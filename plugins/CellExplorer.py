import json
import logging
import math
import os
import shutil
import tempfile
import threading
import gzip
from collections import Counter, OrderedDict, defaultdict
from dataclasses import dataclass
from concurrent.futures import Future, ThreadPoolExecutor, as_completed
from textwrap import dedent
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
from flask import make_response, render_template_string

try:
    import anndata as ad
except ImportError:  # pragma: no cover - runtime guard
    ad = None

from scipy.sparse import load_npz, issparse
from skimage import io, measure
from skimage.segmentation import expand_labels, find_boundaries

from matplotlib import cm, colors
from PySide6.QtCore import QMetaObject, QObject, Qt, Q_ARG, Slot
from PySide6.QtWidgets import QApplication, QFileDialog

LOGGER = logging.getLogger(__name__)

_GLOBAL_STATE: Dict[str, object] = {}
_FILE_DIALOG_HELPER: Optional["FileDialogHelper"] = None

OBS_CATEGORY_LIMIT = 256

_FILETREE_TEMPLATE = dedent(
    """
    <!doctype html>
    <html lang="en">
      <head>
        <meta charset="UTF-8" />
        <meta name="viewport" content="width=device-width, initial-scale=1.0" />
        <title>Select file</title>
        <script src="/static/vendor/jquery-3.5.1/jquery-3.5.1.min.js"></script>
        <link rel="stylesheet" href="/static/vendor/jstree-3.3.16/style.min.css" />
        <script src="/static/vendor/jstree-3.3.16/jstree.min.js"></script>
        <style>
          body { font-family: system-ui, sans-serif; margin: 0; padding: 0.5rem; }
          #b2c-filetree { min-height: 360px; }
        </style>
      </head>
      <body>
        <div id="b2c-filetree"></div>
        <script>
          const pickerField = {{ field|tojson }};
          const slideRoot = {{ base_path|tojson }};
          const startRel = {{ start_rel|tojson }};

          function toAbsolutePath(relPath) {
            const cleanBase = slideRoot.replace(/\\\\/g, "/").replace(/\/+$/, "");
            const cleanRel = relPath.replace(/^\\/+/, "");
            if (!cleanRel) {
              return cleanBase;
            }
            return cleanBase + "/" + cleanRel;
          }

          $(function () {
            const initialRoot = (startRel && startRel.length) ? startRel : ".";
            $("#b2c-filetree")
              .jstree({
                core: {
                  data: {
                    url: function (node) {
                      const target =
                        node.id === "#"
                          ? initialRoot
                          : $("#b2c-filetree").jstree().get_path(node, "/");
                      return "/plugins/CellExplorer/filetree_data?root=" + encodeURIComponent(target);
                    },
                    dataType: "json",
                  },
                  themes: { dots: true, icons: true },
                },
              })
              .on("open_node.jstree", function (e, data) {
                if (data.node.children.length === 0) {
                  $("#b2c-filetree").jstree("load_node", data.node);
                }
              })
              .on("select_node.jstree", function (e, data) {
                if (data.node.data && data.node.data.isdirectory) {
                  data.instance.open_node(data.node);
                  return;
                }
                const relPath = data.instance.get_path(data.node, "/");
                const absolute = toAbsolutePath(relPath);
                window.parent.postMessage({type: "b2c_file_selected", field: pickerField, path: absolute}, "*");
              });
          });
        </script>
      </body>
    </html>
    """
)


def _as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    if isinstance(value, (int, float)):
        return value != 0
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def _safe_int(value: object, default: int) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _safe_float(value: object, default: float) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _normalize_path(path: str) -> str:
    if not path:
        return ""
    return os.path.abspath(os.path.expanduser(os.path.expandvars(path)))


def _rgba_to_css(rgba: Tuple[float, float, float, float], *, alpha_override: Optional[float] = None) -> str:
    r, g, b, a = rgba
    if alpha_override is not None:
        a = alpha_override
    return f"rgba({int(round(r * 255))},{int(round(g * 255))},{int(round(b * 255))},{round(a, 4)})"


def _get_cmap(cmap: object) -> "colors.Colormap":
    if isinstance(cmap, str):
        return cm.get_cmap(cmap)
    return cmap  # type: ignore[return-value]


def _colormap_sample(cmap: object, value: float, alpha: Optional[float] = None) -> str:
    cmap_obj = _get_cmap(cmap)
    rgba = cmap_obj(np.clip(value, 0.0, 1.0))
    return _rgba_to_css(rgba, alpha_override=alpha)


def _sample_gradient(cmap: object, steps: int = 8) -> List[Tuple[float, str]]:
    cmap_obj = _get_cmap(cmap)
    gradient = []
    for i in range(steps):
        pos = i / max(steps - 1, 1)
        gradient.append((pos, _rgba_to_css(cmap_obj(pos))))
    return gradient


def _unique_sorted(iterable: Iterable) -> List:
    return sorted(set(iterable))


def _color_to_css(value: object, *, alpha_override: Optional[float] = None) -> Optional[str]:
    if value is None:
        return None
    try:
        rgba = colors.to_rgba(value)
    except ValueError:
        return None
    return _rgba_to_css(rgba, alpha_override=alpha_override)


@dataclass(frozen=True)
class Tile:
    id: int
    r0: int
    r1: int
    c0: int
    c1: int

    def to_dict(self) -> Dict[str, int]:
        return {"id": self.id, "r0": self.r0, "r1": self.r1, "c0": self.c0, "c1": self.c1}


class B2CContext:
    """
    Lightweight loader for H&E image, sparse nuclei labels, and AnnData container.
    Provides centroid pixel coordinates and convenience accessors.
    """

    def __init__(
        self,
        he_image_path: str,
        labels_npz_path: str,
        adata_path: str,
        *,
        obsm_key: Optional[str] = None,
    ):
        if ad is None:
            raise ImportError("anndata is required (pip install anndata scanpy).")

        he_image_path = os.path.abspath(os.path.expanduser(os.path.expandvars(he_image_path)))
        labels_npz_path = os.path.abspath(os.path.expanduser(os.path.expandvars(labels_npz_path)))
        adata_path = os.path.abspath(os.path.expanduser(os.path.expandvars(adata_path)))

        if not os.path.exists(he_image_path):
            raise FileNotFoundError(f"H&E image not found: {he_image_path}")
        if os.path.isdir(he_image_path):
            raise IsADirectoryError(f"H&E image path is a directory, expected a file: {he_image_path}")
        if not os.path.exists(labels_npz_path):
            raise FileNotFoundError(f"Label NPZ not found: {labels_npz_path}")
        if os.path.isdir(labels_npz_path):
            raise IsADirectoryError(f"Label NPZ path is a directory, expected a file: {labels_npz_path}")
        if not os.path.exists(adata_path):
            raise FileNotFoundError(f"AnnData (.h5ad) not found: {adata_path}")
        if os.path.isdir(adata_path):
            raise IsADirectoryError(f"AnnData path is a directory, expected a file: {adata_path}")

        LOGGER.info("Loading H&E image: %s", he_image_path)
        self.he = io.imread(he_image_path)

        LOGGER.info("Loading sparse labels: %s", labels_npz_path)
        self.lab_sp = load_npz(labels_npz_path)

        if self.lab_sp.shape != self.he.shape[:2]:
            raise ValueError(
                f"Image/label size mismatch: H&E {self.he.shape[:2]} vs labels {self.lab_sp.shape}"
            )

        LOGGER.info("Loading AnnData: %s", adata_path)
        self.adata = ad.read_h5ad(adata_path)

        self._gene_cache: Dict[str, np.ndarray] = {}
        self._obs_cache: Dict[str, np.ndarray] = {}
        self._obs_meta_cache: Dict[str, Dict[str, object]] = {}

        self._available_obsm = [
            k
            for k, arr in self.adata.obsm.items()
            if hasattr(arr, "shape") and arr.shape[1] >= 2
        ]
        if not self._available_obsm:
            raise ValueError("No obsm entries with >=2 columns were found in AnnData.")

        self.obsm_key = self._resolve_obsm_key(obsm_key)
        xy = np.asarray(self.adata.obsm[self.obsm_key])
        self.cols = np.rint(xy[:, 0]).astype(int)
        self.rows = np.rint(xy[:, 1]).astype(int)

        H, W = self.shape
        np.clip(self.cols, 0, W - 1, out=self.cols)
        np.clip(self.rows, 0, H - 1, out=self.rows)

    def _resolve_obsm_key(self, preferred: Optional[str]) -> str:
        if preferred and preferred in self._available_obsm:
            return preferred
        for candidate in ("spatial_cropped_150_buffer", "spatial"):
            if candidate in self._available_obsm:
                return candidate
        return self._available_obsm[0]

    @property
    def shape(self) -> Tuple[int, int]:
        return self.lab_sp.shape

    @property
    def available_obsm(self) -> List[str]:
        return list(self._available_obsm)

    @property
    def gene_names(self) -> List[str]:
        return list(map(str, self.adata.var_names))

    @property
    def obs_columns(self) -> List[str]:
       return list(map(str, self.adata.obs.columns))

    def crop_dense(self, r0: int, r1: int, c0: int, c1: int) -> Tuple[np.ndarray, np.ndarray]:
        he = self.he[r0:r1, c0:c1]
        lab = self.lab_sp[r0:r1, c0:c1].toarray().astype(np.int32, copy=False)
        return he, lab

    def gene_vector(self, gene: str) -> np.ndarray:
        if gene not in self._gene_cache:
            if gene not in self.adata.var_names:
                raise KeyError(f"Gene '{gene}' not found in AnnData.")
            values = self.adata[:, gene].X
            if issparse(values):
                values = values.toarray()
            self._gene_cache[gene] = np.asarray(values).ravel()
        return self._gene_cache[gene]

    def obs_vector(self, column: str) -> np.ndarray:
        if column not in self._obs_cache:
            if column not in self.adata.obs.columns:
                raise KeyError(f"Column '{column}' not found in AnnData.obs.")
            self._obs_cache[column] = self.adata.obs[column].to_numpy()
        return self._obs_cache[column]

    def obs_metadata(self, column: str) -> Dict[str, object]:
        if column not in self.adata.obs.columns:
            raise KeyError(f"Column '{column}' not found in AnnData.obs.")
        cached = self._obs_meta_cache.get(column)
        if cached is not None:
            return cached
        metadata = self._build_obs_metadata(column)
        self._obs_meta_cache[column] = metadata
        return metadata

    def _build_obs_metadata(self, column: str) -> Dict[str, object]:
        series = self.adata.obs[column]
        dtype = getattr(series, "dtype", None)
        categories: List[str] = []
        limit_hit = False

        if dtype is not None and hasattr(dtype, "categories"):
            cat_values = list(getattr(series.cat, "categories", getattr(dtype, "categories", [])))
            categories = [str(cat) for cat in cat_values]
        elif dtype is not None and getattr(dtype, "kind", "") in {"O", "U", "S"}:
            raw_unique = series.dropna().astype(str).unique()
            limit_hit = raw_unique.size > OBS_CATEGORY_LIMIT
            if not limit_hit:
                categories = _unique_sorted(map(str, raw_unique.tolist()))
        else:
            limit_hit = True

        color_key = f"{column}_colors"
        raw_colors = self.adata.uns.get(color_key)
        color_map: Dict[str, str] = {}
        if isinstance(raw_colors, (list, tuple, np.ndarray)) and categories:
            for cat, raw_color in zip(categories, raw_colors):
                if raw_color is None:
                    continue
                raw_color_str = str(raw_color)
                if raw_color_str:
                    color_map[str(cat)] = raw_color_str

        return {
            "categories": categories,
            "color_map": color_map,
            "category_limit_hit": limit_hit,
        }

    def obs_color_map(self, column: str) -> Dict[str, str]:
        meta = self.obs_metadata(column)
        colors_dict = meta.get("color_map") or {}
        return dict(colors_dict)  # shallow copy


def make_tiles(
    ctx: B2CContext, tile_h: int = 1500, tile_w: int = 1500, stride_h: Optional[int] = None, stride_w: Optional[int] = None
) -> List[Tile]:
    H, W = ctx.shape
    stride_h = tile_h if stride_h is None else stride_h
    stride_w = tile_w if stride_w is None else stride_w

    tiles: List[Tile] = []
    tid = 0
    r = 0
    while r < H:
        r1 = min(r + tile_h, H)
        c = 0
        while c < W:
            c1 = min(c + tile_w, W)
            tiles.append(Tile(tid, r, r1, c, c1))
            tid += 1
            if c1 == W:
                break
            c += stride_w
        if r1 == H:
            break
        r += stride_h
    return tiles


def _compute_expansion_distance(
    lab_raw: np.ndarray,
    *,
    mode: str,
    max_bin_distance: float,
    mpp: float,
    bin_um: float,
    volume_ratio: float,
) -> int:
    dist_px_fixed = int(math.ceil(max(1.0, max_bin_distance * (bin_um / mpp))))
    if mode != "volume_ratio":
        return dist_px_fixed

    lab_ids, counts = np.unique(lab_raw[lab_raw > 0], return_counts=True)
    if lab_ids.size == 0:
        return dist_px_fixed

    r_eff = np.sqrt(counts / np.pi)
    delta_r = (np.sqrt(volume_ratio) - 1.0) * r_eff
    dist_px = int(np.clip(np.median(delta_r), 1, 5 * dist_px_fixed))
    return max(1, dist_px)


def _expand_labels_tile(
    ctx: B2CContext,
    tile: Tile,
    *,
    mode: str,
    max_bin_distance: float,
    mpp: float,
    bin_um: float,
    volume_ratio: float,
    pad_factor: int = 2,
) -> Tuple[np.ndarray, np.ndarray, int]:
    r0, r1, c0, c1 = tile.r0, tile.r1, tile.c0, tile.c1
    he_crop, lab_raw = ctx.crop_dense(r0, r1, c0, c1)

    dist_px = _compute_expansion_distance(
        lab_raw,
        mode=mode,
        max_bin_distance=max_bin_distance,
        mpp=mpp,
        bin_um=bin_um,
        volume_ratio=volume_ratio,
    )

    H, W = ctx.shape
    pad = pad_factor * dist_px

    rp0 = max(r0 - pad, 0)
    rp1 = min(r1 + pad, H)
    cp0 = max(c0 - pad, 0)
    cp1 = min(c1 + pad, W)

    lab_pad = ctx.lab_sp[rp0:rp1, cp0:cp1].toarray().astype(np.int32, copy=False)
    lab_exp_pad = expand_labels(lab_pad, distance=dist_px)

    r0_rel, r1_rel = r0 - rp0, r1 - rp0
    c0_rel, c1_rel = c0 - cp0, c1 - cp0

    lab_exp = lab_exp_pad[r0_rel:r1_rel, c0_rel:c1_rel]
    return he_crop, lab_exp.astype(np.int32, copy=False), int(dist_px)


def _polygons_from_labels(lab: np.ndarray, tile: Tile) -> Dict[int, List[List[List[float]]]]:
    """
    FAST polygon extraction using cv2.findContours in C++.
    This is 10-100x faster than the pure Python approach.
    """
    import cv2
    
    r0, c0 = tile.r0, tile.c0
    output: Dict[int, List[List[List[float]]]] = {}
    
    # Get unique labels (excluding background)
    unique_labels = np.unique(lab)
    unique_labels = unique_labels[unique_labels > 0]
    
    if len(unique_labels) == 0:
        return output
    
    # Process each label - cv2.findContours is MUCH faster than skimage
    for lbl in unique_labels:
        mask = (lab == lbl).astype(np.uint8) * 255
        
        # cv2.findContours is implemented in C++ and is very fast
        contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        
        poly_list = []
        for contour in contours:
            if len(contour) < 3:
                continue
            # contour is shape (N, 1, 2) - reshape to (N, 2)
            contour = contour.reshape(-1, 2)
            polygon = [[float(x + c0), float(y + r0)] for x, y in contour]
            if polygon and polygon[0] != polygon[-1]:
                polygon.append(polygon[0])
            if polygon:
                poly_list.append(polygon)
        
        if poly_list:
            output[int(lbl)] = poly_list
    
    return output


def _outline_paths(lab: np.ndarray, tile: Tile) -> List[List[List[float]]]:
    r0, c0 = tile.r0, tile.c0
    edges = find_boundaries(lab, mode="inner")
    contours = measure.find_contours(edges.astype(np.uint8), 0.5)
    paths: List[List[List[float]]] = []
    for contour in contours:
        if contour.shape[0] < 2:
            continue
        path: List[List[float]] = []
        for y, x in contour:
            path.append([float(x + c0), float(y + r0)])
        paths.append(path)
    return paths

def _outline_paths_per_label(lab: np.ndarray, tile: Tile, labels: Optional[List[int]] = None) -> Dict[int, List[List[List[float]]]]:
    """
    Generate outline paths for each individual label using cv2 (fast C++ implementation).
    """
    import cv2
    
    r0, c0 = tile.r0, tile.c0
    result: Dict[int, List[List[List[float]]]] = {}
    
    unique_labels = labels if labels is not None else np.unique(lab).tolist()
    unique_labels = [lbl for lbl in unique_labels if lbl != 0]
    
    if not unique_labels:
        return result
    
    for lbl in unique_labels:
        mask = (lab == lbl).astype(np.uint8) * 255
        contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        
        paths = []
        for contour in contours:
            if len(contour) < 2:
                continue
            contour = contour.reshape(-1, 2)
            path = [[float(x + c0), float(y + r0)] for x, y in contour]
            paths.append(path)
        
        if paths:
            result[int(lbl)] = paths
    
    return result



def _centroids_for_tile(ctx: B2CContext, tile: Tile) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    r0, r1, c0, c1 = tile.r0, tile.r1, tile.c0, tile.c1
    in_tile = (ctx.rows >= r0) & (ctx.rows < r1) & (ctx.cols >= c0) & (ctx.cols < c1)
    indices = np.where(in_tile)[0]
    if indices.size == 0:
        empty_int = np.empty((0,), dtype=np.int32)
        empty_local = np.empty((0, 2), dtype=np.int32)
        return empty_int, empty_int, empty_int, empty_local

    rr_local = (ctx.rows[indices] - r0).astype(np.int32)
    cc_local = (ctx.cols[indices] - c0).astype(np.int32)
    local = np.stack([rr_local, cc_local], axis=1)
    return indices, ctx.rows[indices], ctx.cols[indices], local


class TileCache:
    def __init__(self, max_items: int = 6):
        self._max = max_items
        self._data: "OrderedDict[Tuple, Dict]" = OrderedDict()
        self._lock = threading.RLock()

    def clear(self) -> None:
        with self._lock:
            self._data.clear()

    def get(self, key: Tuple) -> Optional[Dict]:
        with self._lock:
            entry = self._data.get(key)
            if entry is not None:
                self._data.move_to_end(key)
            return entry

    def set(self, key: Tuple, value: Dict) -> None:
        with self._lock:
            self._data[key] = value
            self._data.move_to_end(key)
            while len(self._data) > self._max:
                self._data.popitem(last=False)


class GeometryCache:
    """
    Cache for tile geometry (polygons, outlines, centroids).
    Invalidated ONLY when geometry-affecting parameters change.
    
    Geometry parameters: b2c_mode, max_bin_distance, mpp, bin_um, volume_ratio, pad_factor
    """
    def __init__(self, max_items: int = 100):
        self._max = max_items
        self._data: "OrderedDict[Tuple, Dict]" = OrderedDict()
        self._lock = threading.RLock()
        self._hits = 0
        self._misses = 0

    def _make_key(
        self,
        dataset_id: str,
        tile_id: int,
        b2c_mode: str,
        max_bin_distance: float,
        mpp: float,
        bin_um: float,
        volume_ratio: float,
        pad_factor: int,
    ) -> Tuple:
        """
        Create cache key from ONLY geometry-affecting parameters.
        Round floats to avoid cache misses from floating-point precision.
        """
        return (
            dataset_id,
            tile_id,
            b2c_mode,
            round(max_bin_distance, 3),
            round(mpp, 4),
            round(bin_um, 4),
            round(volume_ratio, 4),
            pad_factor,
        )

    def get(self, key: Tuple) -> Optional[Dict]:
        """Thread-safe cache get with LRU update"""
        with self._lock:
            entry = self._data.get(key)
            if entry is not None:
                self._data.move_to_end(key)
                self._hits += 1
                return entry
            else:
                self._misses += 1
                return None

    def set(self, key: Tuple, value: Dict) -> None:
        """Thread-safe cache set with LRU eviction"""
        with self._lock:
            self._data[key] = value
            self._data.move_to_end(key)

            # Evict oldest if over limit
            while len(self._data) > self._max:
                evicted_key, _ = self._data.popitem(last=False)
                LOGGER.debug(f"GeometryCache: Evicted {evicted_key}")

    def stats(self) -> Dict[str, object]:
        """Get cache statistics"""
        with self._lock:
            total = self._hits + self._misses
            hit_rate = (self._hits / total * 100) if total > 0 else 0
            return {
                "size": len(self._data),
                "max": self._max,
                "hits": self._hits,
                "misses": self._misses,
                "hit_rate": round(hit_rate, 2),
            }

    def clear(self) -> None:
        """Clear all cached geometry"""
        with self._lock:
            self._data.clear()
            self._hits = 0
            self._misses = 0


class ColorCache:
    """
    Cache for colored overlays.
    Invalidated when color parameters OR geometry parameters change.
    Higher capacity than GeometryCache since colors change more frequently.
    
    Color parameters: overlay_type, gene, obs_col, category, color_mode, gradient_color, expr_quantile
    """
    def __init__(self, max_items: int = 200):
        self._max = max_items
        self._data: "OrderedDict[Tuple, Dict]" = OrderedDict()
        self._lock = threading.RLock()
        self._hits = 0
        self._misses = 0

    def _make_key(
        self,
        dataset_id: str,
        tile_id: int,
        overlay_type: str,
        geometry_key: Tuple,
        gene: Optional[str] = None,
        obs_col: Optional[str] = None,
        category: Optional[str] = None,
        color_mode: str = "gradient",
        gradient_color: Optional[str] = None,
        expr_quantile: Optional[float] = None,
    ) -> Tuple:
        """
        Create cache key from color parameters + geometry key.
        
        The geometry_key ensures that color cache is invalidated when
        geometry changes (since we need to recolor the new geometry).
        """
        if overlay_type == "gene":
            return (
                dataset_id,
                tile_id,
                "gene",
                geometry_key,  # Include geometry key for invalidation
                gene,
                color_mode,
                gradient_color,
                round(expr_quantile, 4) if expr_quantile else None,
            )
        else:  # observation
            return (
                dataset_id,
                tile_id,
                "obs",
                geometry_key,
                obs_col,
                category or "__all__",
            )

    def get(self, key: Tuple) -> Optional[Dict]:
        """Thread-safe cache get with LRU update"""
        with self._lock:
            entry = self._data.get(key)
            if entry is not None:
                self._data.move_to_end(key)
                self._hits += 1
                return entry
            else:
                self._misses += 1
                return None

    def set(self, key: Tuple, value: Dict) -> None:
        """Thread-safe cache set with LRU eviction"""
        with self._lock:
            self._data[key] = value
            self._data.move_to_end(key)

            while len(self._data) > self._max:
                self._data.popitem(last=False)

    def invalidate_geometry(self, geometry_key: Tuple) -> int:
        """
        Invalidate all color cache entries for a specific geometry.
        Called when geometry parameters change.
        
        Returns: Number of entries invalidated
        """
        with self._lock:
            keys_to_remove = [
                k for k in self._data.keys() if len(k) > 3 and k[3] == geometry_key
            ]
            for k in keys_to_remove:
                del self._data[k]
            return len(keys_to_remove)

    def stats(self) -> Dict[str, object]:
        """Get cache statistics"""
        with self._lock:
            total = self._hits + self._misses
            hit_rate = (self._hits / total * 100) if total > 0 else 0
            return {
                "size": len(self._data),
                "max": self._max,
                "hits": self._hits,
                "misses": self._misses,
                "hit_rate": round(hit_rate, 2),
            }

    def clear(self) -> None:
        """Clear all cached colors"""
        with self._lock:
            self._data.clear()
            self._hits = 0
            self._misses = 0


class Plugin:
    def __init__(self, app):
        self.app = app
        self.out_dir = os.path.expanduser("~/.tissuumaps/plugins/CellExplorer")
        os.makedirs(self.out_dir, exist_ok=True)
        self.log_path = os.path.join(self.out_dir, "CellExplorer.log")
        if not LOGGER.handlers:
            handler = logging.FileHandler(self.log_path)
            handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
            LOGGER.addHandler(handler)
        LOGGER.setLevel(logging.INFO)

        # CRITICAL FIX: Use app-level storage instead of module-level global
        # This survives module reloading which TissUUmaps does on each request
        if not hasattr(app, '_cell_explorer_state'):
            LOGGER.info("Initializing NEW app-level state (first time or after app restart)")
            app._cell_explorer_state = {
                "context": None,
                "tiles": [],
                "dataset_id": None,
                "tile_cache": TileCache(max_items=32),
                "geometry_cache": GeometryCache(max_items=100),
                "color_cache": ColorCache(max_items=200),
                "tile_cache_size": 32,
                "dataset_config": None,
                "loading": False,
                "cache_warm_job": None,
                "staged_dir": None,
                "staged_dir_managed": False,
                "worker_pool": None,
                "worker_count": 1,
                "tile_jobs": {},
                "tile_jobs_lock": threading.RLock(),
                "single_tile_mode": False,
            }
        else:
            LOGGER.info(f"Reusing EXISTING app-level state (context exists: {app._cell_explorer_state.get('context') is not None})")
        
        self._state = app._cell_explorer_state
        self.current_params: Dict[str, str] = {}
        self.state_path = os.path.join(self.out_dir, "dataset_state.json")
        if self._state.get("dataset_config") is None:
            cached = self._load_cached_config()
            if cached:
                self._state["dataset_config"] = cached

        self.presets_path = os.path.join(self.out_dir, "presets.json")
        self.presets = self._load_presets()

    # Helpers
    def _json_response(self, payload: Dict) -> "flask.Response":
        return make_response(json.dumps(payload, default=_json_default), 200, {"Content-Type": "application/json"})

    def _json_response_compressed(self, payload: Dict) -> "flask.Response":
        """
        Return gzip-compressed JSON response.
        
        Reduces payload size by 70-90% for typical JSON responses.
        Browsers automatically decompress gzip responses.
        """
        json_str = json.dumps(payload, default=_json_default)
        json_bytes = json_str.encode('utf-8')
        
        # Compress (level 6 is good balance of speed vs compression)
        compressed = gzip.compress(json_bytes, compresslevel=6)
        
        response = make_response(compressed, 200)
        response.headers['Content-Type'] = 'application/json'
        response.headers['Content-Encoding'] = 'gzip'
        response.headers['Content-Length'] = len(compressed)
        
        LOGGER.debug(
            "Response: %d bytes → %d bytes (%.1f%% of original)",
            len(json_bytes), len(compressed), len(compressed) / len(json_bytes) *  100
        )
        
        return response

    def _error_response(self, status: int, message: str, exc: Optional[Exception] = None):
        if exc is not None:
            LOGGER.exception(message)
        else:
            LOGGER.error(message)
        payload = {"status": "error", "message": message}
        if exc is not None:
            payload["detail"] = repr(exc)
        return make_response(json.dumps(payload), status, {"Content-Type": "application/json"})

    def _error_from_exception(self, exc: Exception, context: str):
        if isinstance(exc, FileNotFoundError):
            status = 404
            message = str(exc)
        elif isinstance(exc, (ValueError, KeyError, RuntimeError)):
            status = 400
            message = str(exc)
        elif isinstance(exc, ImportError):
            status = 500
            message = f"{context}: missing dependency ({exc})"
        else:
            status = 500
            message = f"{context}: {exc}"
        return self._error_response(status, message, exc)

    def _ensure_context(self) -> B2CContext:
        if self.context is None:
            config = self._state.get("dataset_config") or self._load_cached_config()
            if config:
                LOGGER.warning("⚠️ Context is None, rehydrating dataset (this should only happen once per session!)")
                ctx = B2CContext(
                    he_image_path=config["he_path"],
                    labels_npz_path=config["labels_path"],
                    adata_path=config["adata_path"],
                    obsm_key=config.get("obsm_key"),
                )
                tiles = make_tiles(
                    ctx,
                    tile_h=config.get("tile_h", 1500) or 1500,
                    tile_w=config.get("tile_w", 1500) or 1500,
                    stride_h=config.get("stride_h"),
                    stride_w=config.get("stride_w"),
                )
                self.context = ctx
                self.tiles = tiles
                self.dataset_id = os.path.basename(config["adata_path"])
                
                # IMPORTANT: Only clear old TileCache, NOT the new caches!
                self.tile_cache.clear()
                # DO NOT clear geometry_cache or color_cache here!
                
                self._state["dataset_config"] = config
                LOGGER.info(f"Context rehydrated. Dataset ID: {self.dataset_id}")
            else:
                raise RuntimeError("Load a dataset first.")
        return self.context

    @property
    def context(self) -> Optional[B2CContext]:
        ctx = self._state.get("context")
        if ctx is None:
            LOGGER.debug(f"🔴 context getter: returning None (app._cell_explorer_state id: {id(self._state)})")
        else:
            LOGGER.debug(f"context getter: returning context (ctx id: {id(ctx)}, state id: {id(self._state)})")
        return ctx  # type: ignore[return-value]

    @context.setter
    def context(self, value: Optional[B2CContext]) -> None:
        if value is None:
            LOGGER.warning(f"🔴 context setter: Setting context to None! (state id: {id(self._state)})")
        else:
            LOGGER.info(f"context setter: Storing context (ctx id: {id(value)}, state id: {id(self._state)})")
        self._state["context"] = value

    @property
    def tiles(self) -> List[Tile]:
        return self._state.get("tiles", [])  # type: ignore[return-value]

    @tiles.setter
    def tiles(self, value: List[Tile]) -> None:
        self._state["tiles"] = value

    @property
    def dataset_id(self) -> Optional[str]:
        return self._state.get("dataset_id")  # type: ignore[return-value]

    @dataset_id.setter
    def dataset_id(self, value: Optional[str]) -> None:
        self._state["dataset_id"] = value

    @property
    def tile_cache(self) -> TileCache:
        cache = self._state.get("tile_cache")
        if cache is None:
            cache = TileCache(max_items=self._state.get("tile_cache_size", 32))
            self._state["tile_cache"] = cache
        return cache  # type: ignore[return-value]

    def _load_presets(self) -> Dict[str, Dict]:
        if not os.path.exists(self.presets_path):
            return {}
        try:
            with open(self.presets_path, "r") as fh:
                presets = json.load(fh)
            return presets if isinstance(presets, dict) else {}
        except Exception as exc:  # pragma: no cover - runtime guard
            LOGGER.warning("Failed to read presets: %s", exc)
            return {}

    def _save_presets(self) -> None:
        with open(self.presets_path, "w") as fh:
            json.dump(self.presets, fh, indent=2)

    def _config_paths_valid(self, config: Dict) -> bool:
        required = ["he_path", "labels_path", "adata_path"]
        for key in required:
            path = config.get(key)
            if not path or not os.path.exists(os.path.abspath(os.path.expanduser(str(path)))):
                LOGGER.warning("Cached path for %s is missing: %s", key, path)
                return False
        return True

    def _normalize_tile_cache_size(self, value: object) -> int:
        default_size = self._state.get("tile_cache_size") or 32
        size = _safe_int(value, default_size)
        return max(1, min(512, size))

    def _configure_tile_cache(self, size: int) -> None:
        size = self._normalize_tile_cache_size(size)
        current_size = self._state.get("tile_cache_size", 4)
        cache = self._state.get("tile_cache")
        if cache is None:
            cache = TileCache(max_items=size)
        if current_size == size and cache is self._state.get("tile_cache"):
            self._state["tile_cache_size"] = size
            return
        new_cache = TileCache(max_items=size)
        if hasattr(cache, "_data"):
            for key, value in cache._data.items():
                new_cache.set(key, value)
        self._state["tile_cache"] = new_cache
        self._state["tile_cache_size"] = size

    def _resolve_stage_root(self, override: Optional[str]) -> Optional[str]:
        candidates = [
            override,
            os.environ.get("TISSUUMAPS_STAGE_ROOT"),
            os.environ.get("BIN2CELL_STAGE_ROOT"),
            os.environ.get("SLURM_TMPDIR"),
            os.environ.get("TMPDIR"),
        ]
        for candidate in candidates:
            if not candidate:
                continue
            normalized = _normalize_path(candidate)
            if not normalized:
                continue
            try:
                os.makedirs(normalized, exist_ok=True)
                return normalized
            except OSError as exc:
                LOGGER.warning("Failed to prepare staging root %s: %s", normalized, exc)
        return None

    def _cleanup_staged_dir(self) -> None:
        staged_dir = self._state.get("staged_dir")
        managed = self._state.get("staged_dir_managed", False)
        if staged_dir and managed and os.path.isdir(staged_dir):
            try:
                shutil.rmtree(staged_dir)
                LOGGER.info("Removed previous staged dataset at %s", staged_dir)
            except Exception as exc:
                LOGGER.warning("Failed to remove staged dataset directory %s: %s", staged_dir, exc)
        self._state["staged_dir"] = None
        self._state["staged_dir_managed"] = False

    def _stage_dataset_files(
        self,
        *,
        he_path: str,
        labels_path: str,
        adata_path: str,
        stage_root: str,
    ) -> Tuple[Dict[str, str], Dict[str, object]]:
        os.makedirs(stage_root, exist_ok=True)
        dest_root = tempfile.mkdtemp(prefix="b2c_stage_", dir=stage_root)
        staged_paths: Dict[str, str] = {}
        inputs = {
            "he_path": he_path,
            "labels_path": labels_path,
            "adata_path": adata_path,
        }
        for key, src in inputs.items():
            normalized_src = _normalize_path(src)
            if not normalized_src:
                continue
            # Validate that the path exists and is a file, not a directory
            if not os.path.exists(normalized_src):
                raise FileNotFoundError(f"File not found for {key}: {normalized_src}")
            if os.path.isdir(normalized_src):
                raise IsADirectoryError(f"Path for {key} is a directory, expected a file: {normalized_src}")
            try:
                common = os.path.commonpath([normalized_src, stage_root])
            except ValueError:
                common = ""
            if common == stage_root:
                staged_paths[key] = normalized_src
                continue
            dest_path = os.path.join(dest_root, os.path.basename(normalized_src))
            LOGGER.info("Staging %s -> %s", normalized_src, dest_path)
            shutil.copy2(normalized_src, dest_path)
            staged_paths[key] = dest_path
        self._state["staged_dir"] = dest_root
        self._state["staged_dir_managed"] = True
        return staged_paths, {
            "staged_dir": dest_root,
            "stage_root": stage_root,
            "files": staged_paths,
        }

    def _stop_cache_warmer(self) -> None:
        job = self._state.get("cache_warm_job")
        if not job:
            return
        stop_event = job.get("stop_event")
        if stop_event:
            stop_event.set()
        thread = job.get("thread")
        if thread and thread.is_alive():
            thread.join(timeout=0.1)
        self._state["cache_warm_job"] = None

    def _start_cache_warmer(
        self,
        *,
        enabled: bool,
        tile_limit: int,
        warm_params: Dict[str, object],
    ) -> None:
        self._stop_cache_warmer()
        if not enabled:
            return
        if not self.tiles:
            return
        if self.context is None:
            return
        stop_event = threading.Event()
        total = tile_limit if tile_limit > 0 else len(self.tiles)
        job: Dict[str, object] = {
            "stop_event": stop_event,
            "progress": 0,
            "total": total,
            "params": warm_params,
        }

        def worker() -> None:
            warmed = 0
            future_map: Dict[Future, int] = {}
            scheduled = 0
            for tile in self.tiles:
                if stop_event.is_set():
                    break
                try:
                    future = self._schedule_tile_job(
                        tile.id,
                        b2c_mode=warm_params.get("b2c_mode", "fixed"),
                        max_bin_distance=_safe_float(warm_params.get("max_bin_distance"), 2.0),
                        mpp=_safe_float(warm_params.get("mpp"), 0.3),
                        bin_um=_safe_float(warm_params.get("bin_um"), 2.0),
                        volume_ratio=_safe_float(warm_params.get("volume_ratio"), 4.0),
                        pad_factor=max(1, _safe_int(warm_params.get("pad_factor"), 2)),
                    )
                    future_map[future] = tile.id
                    scheduled += 1
                except Exception as exc:  # pragma: no cover - background log
                    LOGGER.warning("Cache warm scheduling failed for tile %s: %s", tile.id, exc)
                    continue
                if tile_limit and scheduled >= tile_limit:
                    break
            job["total"] = scheduled
            if not future_map:
                job["done"] = True
                return
            for future in as_completed(list(future_map.keys())):
                if stop_event.is_set():
                    break
                tile_id = future_map.get(future)
                try:
                    future.result()
                except Exception as exc:  # pragma: no cover - background log
                    LOGGER.warning("Cache warm failed for tile %s: %s", tile_id, exc)
                warmed += 1
                job["progress"] = warmed
                if tile_limit and warmed >= tile_limit:
                    break
            job["done"] = True

        thread = threading.Thread(target=worker, name="Bin2CellCacheWarmer", daemon=True)
        job["thread"] = thread
        self._state["cache_warm_job"] = job
        thread.start()

    def _cache_warm_status(self) -> Dict[str, object]:
        job = self._state.get("cache_warm_job")
        if not job:
            return {"active": False, "progress": 0, "total": 0}
        thread = job.get("thread")
        active = bool(thread and thread.is_alive())
        status = {
            "active": active,
            "progress": job.get("progress", 0),
            "total": job.get("total") or len(self.tiles),
        }
        if job.get("done"):
            status["active"] = False
            status["done"] = True
        return status

    def _shutdown_worker_pool(self) -> None:
        pool = self._state.get("worker_pool")
        if pool:
            pool.shutdown(wait=False)
        self._state["worker_pool"] = None
        self._state["worker_count"] = 0
        jobs = self._state.get("tile_jobs")
        if jobs:
            for future in jobs.values():
                if isinstance(future, Future):
                    future.cancel()
        self._state["tile_jobs"] = {}

    def _configure_worker_pool(self, requested: Optional[int]) -> ThreadPoolExecutor:
        desired = max(1, min(_safe_int(requested, os.cpu_count() or 1), (os.cpu_count() or 1)))
        pool = self._state.get("worker_pool")
        current = self._state.get("worker_count", 0)
        if pool and current == desired:
            return pool
        if pool:
            pool.shutdown(wait=False)
        pool = ThreadPoolExecutor(
            max_workers=desired,
            thread_name_prefix="Bin2CellTile",
        )
        self._state["worker_pool"] = pool
        self._state["worker_count"] = desired
        if "tile_jobs" not in self._state:
            self._state["tile_jobs"] = {}
        if "tile_jobs_lock" not in self._state:
            self._state["tile_jobs_lock"] = threading.RLock()
        return pool

    def _tile_cache_key(
        self,
        tile_id: int,
        *,
        b2c_mode: str,
        max_bin_distance: float,
        mpp: float,
        bin_um: float,
        volume_ratio: float,
        pad_factor: int,
    ) -> Tuple:
        return (
            self.dataset_id,
            tile_id,
            b2c_mode,
            round(max_bin_distance, 3),
            round(mpp, 4),
            round(bin_um, 4),
            round(volume_ratio, 4),
            pad_factor,
        )

    def _schedule_tile_job(
        self,
        tile_id: int,
        *,
        b2c_mode: str,
        max_bin_distance: float,
        mpp: float,
        bin_um: float,
        volume_ratio: float,
        pad_factor: int,
    ) -> Future:
        cache_key = self._tile_cache_key(
            tile_id,
            b2c_mode=b2c_mode,
            max_bin_distance=max_bin_distance,
            mpp=mpp,
            bin_um=bin_um,
            volume_ratio=volume_ratio,
            pad_factor=pad_factor,
        )
        cached = self.tile_cache.get(cache_key)
        if cached is not None:
            future: Future = Future()
            future.set_result(cached)
            return future

        single_tile_mode = bool(self._state.get("single_tile_mode"))
        if single_tile_mode:
            self.tile_cache.clear()
            self._cancel_tile_jobs(keep_key=cache_key)

        pool = self._state.get("worker_pool")
        if pool is None:
            pool = self._configure_worker_pool(self._state.get("worker_count") or 1)

        jobs = self._state.setdefault("tile_jobs", {})
        lock = self._state.setdefault("tile_jobs_lock", threading.RLock())
        with lock:
            existing = jobs.get(cache_key)
            if existing:
                return existing

            future = pool.submit(
                self._tile_entry,
                tile_id,
                b2c_mode=b2c_mode,
                max_bin_distance=max_bin_distance,
                mpp=mpp,
                bin_um=bin_um,
                volume_ratio=volume_ratio,
                pad_factor=pad_factor,
                cache_key=cache_key,
            )

            def _cleanup(fut: Future) -> None:
                with lock:
                    jobs.pop(cache_key, None)

            future.add_done_callback(_cleanup)
            jobs[cache_key] = future
            return future

    def _cancel_tile_jobs(self, keep_key: Optional[Tuple] = None) -> None:
        jobs = self._state.get("tile_jobs")
        if not jobs:
            return
        lock = self._state.setdefault("tile_jobs_lock", threading.RLock())
        with lock:
            for key in list(jobs.keys()):
                if keep_key is not None and key == keep_key:
                    continue
                future = jobs.pop(key, None)
                if isinstance(future, Future):
                    future.cancel()

    #  Endpoints

    def pick_file(self, params: Dict):
        try:
            field = str(params.get("field") or "")
            if field not in {"he_path", "labels_path", "adata_path"}:
                raise ValueError("Unsupported field for file picker.")

            current = params.get("current_path") or ""
            start_dir = params.get("start_dir") or ""
            if current:
                current = os.path.abspath(os.path.expanduser(os.path.expandvars(str(current))))
                if os.path.isdir(current):
                    start_dir = current
                else:
                    start_dir = os.path.dirname(current)
            if not start_dir:
                start_dir = os.path.expanduser("~")

            slide_root = os.path.abspath(self.app.config.get("SLIDE_DIR", "."))
            start_rel = ""
            try:
                start_abs = os.path.abspath(start_dir)
                if os.path.commonpath([start_abs, slide_root]) == slide_root:
                    rel = os.path.relpath(start_abs, slide_root)
                    start_rel = "." if rel == "." else rel
            except Exception:
                start_rel = ""

            use_web_picker = not self.app.config.get("isStandalone")
            if use_web_picker:
                return self._json_response(
                    {
                        "status": "web",
                        "field": field,
                        "start_rel": start_rel,
                    }
                )

            filters = {
                "he_path": "H&E images (*.tif *.tiff *.TIF *.TIFF);;All files (*)",
                "labels_path": "Label matrices (*.npz *.NPZ);;All files (*)",
                "adata_path": "AnnData (*.h5ad *.H5AD);;All files (*)",
            }
            captions = {
                "he_path": "Select H&E image",
                "labels_path": "Select label matrix (.npz)",
                "adata_path": "Select AnnData (.h5ad)",
            }

            app = QApplication.instance()
            if app is None:
                raise RuntimeError("QApplication instance not found; cannot open file dialog.")

            global _FILE_DIALOG_HELPER
            if _FILE_DIALOG_HELPER is None:
                helper = FileDialogHelper()
                helper.moveToThread(app.thread())
                _FILE_DIALOG_HELPER = helper

            filename, _ = _FILE_DIALOG_HELPER.get_open_file_name(
                captions[field],
                start_dir,
                filters[field],
            )

            if not filename:
                return self._json_response({"status": "cancelled"})

            path = os.path.abspath(filename)
            if not os.path.isfile(path):
                raise FileNotFoundError(f"Selected path is not a file: {path}")
            return self._json_response({"status": "ok", "path": path, "field": field})
        except Exception as exc:
            return self._error_from_exception(exc, "Failed to pick file")

    def _resolve_root(self, rel_path: str) -> str:
        base = os.path.abspath(self.app.config.get("SLIDE_DIR", "."))
        target = os.path.abspath(os.path.join(base, rel_path.lstrip("/")))
        if os.path.commonpath([target, base]) != base:
            raise ValueError("Path outside slide directory")
        if not os.path.exists(target):
            raise FileNotFoundError(target)
        return target

    #clean this up
    def _normalize_path_for_response(self, path: str) -> str:
        """
        Convert a file path to the format expected by TissUUmaps based on mode.
        - In web server mode: convert absolute paths within SLIDE_DIR to relative paths
        - In standalone mode: return absolute paths as-is
        """
        if not path:
            return path
        is_standalone = self.app.config.get("isStandalone", False)
        if is_standalone:
            # In standalone mode, return absolute paths as is
            return os.path.abspath(os.path.expanduser(os.path.expandvars(path)))
        # In web server mode, convert to relative if within SLIDE_DIR
        slide_dir = os.path.abspath(self.app.config.get("SLIDE_DIR", "."))
        try:
            # Try to normalize the path
            abs_path = os.path.abspath(os.path.expanduser(os.path.expandvars(path)))
            if os.path.commonpath([abs_path, slide_dir]) == slide_dir:
                rel_path = os.path.relpath(abs_path, slide_dir)
                if rel_path == ".":
                    # If file is at root of SLIDE_DIR, return just the filename
                    return os.path.basename(abs_path)
                return rel_path
        except (ValueError, OSError):

            pass
        # If path is outside SLIDE_DIR or normalization failed, return as-is
        # (TissUUmaps will handle the error appropriately)
        return path

    def filetree(self, params: Dict):
        field = str(params.get("field") or "")
        if field not in {"he_path", "labels_path", "adata_path"}:
            return self._error_response(400, "Unsupported field")
        start_rel = params.get("start_rel") or ""
        base_path = os.path.abspath(self.app.config.get("SLIDE_DIR", "."))
        html = render_template_string(
            _FILETREE_TEMPLATE,
            field=field,
            base_path=base_path,
            start_rel=start_rel,
        )
        # Return HTML as JSON string for API consumption
        return self._json_response({"html": html})

    def filetree_data(self, params: Dict):
        rel_root = params.get("root") or "."
        try:
            root = self._resolve_root(rel_root)
        except Exception as exc:
            return self._error_from_exception(exc, "Invalid root")
        entries = []
        try:
            with os.scandir(root) as it:
                for entry in sorted(it, key=lambda e: e.name.lower()):
                    if entry.name.startswith("."):
                        continue
                    rel_path = os.path.relpath(entry.path, os.path.abspath(self.app.config.get("SLIDE_DIR", ".")))
                    if entry.is_dir():
                        entries.append(
                            {
                                "text": entry.name,
                                "icon": "jstree-folder",
                                "state": {"opened": False},
                                "data": {"isdirectory": True},
                                "children": True,
                            }
                        )
                    else:
                        entries.append(
                            {
                                "text": entry.name,
                                "icon": "jstree-file",
                                "data": {"isdirectory": False, "relpath": rel_path},
                            }
                        )
        except Exception as exc:
            return self._error_from_exception(exc, "Failed to read directory")
        return self._json_response(entries)

    def _load_cached_config(self) -> Optional[Dict[str, object]]:
        if not hasattr(self, "state_path"):
            return None
        if not os.path.exists(self.state_path):
            return None
        try:
            with open(self.state_path, "r") as fh:
                data = json.load(fh)
            if isinstance(data, dict) and self._config_paths_valid(data):
                return data
        except Exception as exc:
            LOGGER.warning("Failed to read cached dataset config: %s", exc)
        return None

    def load_dataset(self, params: Dict):
        try:
            self._state["loading"] = True
            self._stop_cache_warmer()
            self._shutdown_worker_pool()
            he_path = params.get("he_path", "")
            labels_path = params.get("labels_path", "")
            adata_path = params.get("adata_path", "")
            LOGGER.info(
                "load_dataset request paths: he_path=%r, labels_path=%r, adata_path=%r",
                he_path,
                labels_path,
                adata_path,
            )
            obsm_key = params.get("obsm_key")
            tile_h = int(params.get("tile_h", 1500) or 1500)
            tile_w = int(params.get("tile_w", 1500) or 1500)
            stride_h_param = params.get("stride_h")
            stride_w_param = params.get("stride_w")
            stride_h = int(stride_h_param) if stride_h_param not in (None, "", "None") else None
            stride_w = int(stride_w_param) if stride_w_param not in (None, "", "None") else None
            single_tile_mode = _as_bool(params.get("single_tile_mode"))
            tile_cache_size = self._normalize_tile_cache_size(params.get("tile_cache_size"))
            if single_tile_mode:
                tile_cache_size = 1
            self._configure_tile_cache(tile_cache_size)
            tile_workers = max(1, min(_safe_int(params.get("tile_workers"), os.cpu_count() or 1), os.cpu_count() or 1))
            self._configure_worker_pool(tile_workers)
            stage_to_local = _as_bool(params.get("stage_to_local"))
            stage_root_override = params.get("stage_root")
            stage_root = self._resolve_stage_root(stage_root_override) if stage_to_local else None
            if stage_to_local and not stage_root:
                raise ValueError(
                    "Local staging requested but no scratch directory is available. "
                    "Set 'stage_root' or export TISSUUMAPS_STAGE_ROOT / SLURM_TMPDIR."
                )
            warm_cache_requested = _as_bool(params.get("warm_cache", True))
            warm_cache = warm_cache_requested and not single_tile_mode
            warm_cache_tiles = 0 if single_tile_mode else max(0, _safe_int(params.get("warm_cache_tiles"), 0))
            warm_params = {
                "b2c_mode": params.get("warm_b2c_mode") or params.get("b2c_mode") or "fixed",
                "max_bin_distance": _safe_float(params.get("warm_max_bin_distance"), 2.0),
                "mpp": _safe_float(params.get("warm_mpp"), 0.3),
                "bin_um": _safe_float(params.get("warm_bin_um"), 2.0),
                "volume_ratio": _safe_float(params.get("warm_volume_ratio"), 4.0),
                "pad_factor": max(1, _safe_int(params.get("warm_pad_factor"), 2)),
            }


            source_he_path = _normalize_path(he_path)
            source_labels_path = _normalize_path(labels_path)
            source_adata_path = _normalize_path(adata_path)

            # Validate paths early with clear error messages
            if not source_he_path or not source_labels_path or not source_adata_path:
                missing = []
                if not source_he_path:
                    missing.append("H&E image (.tif/.tiff)")
                if not source_labels_path:
                    missing.append("labels NPZ file")
                if not source_adata_path:
                    missing.append("AnnData (.h5ad) file")

                raise ValueError(f"Missing required files: {', '.join(missing)}. Please use the 'Browse' buttons or enter the paths manually in the text fields.")
            
            if not os.path.exists(source_he_path):
                raise FileNotFoundError(f"H&E image file not found: {source_he_path}")
            if os.path.isdir(source_he_path):
                raise IsADirectoryError(f"H&E image path is a directory, not a file: {source_he_path}. Please select the image file (e.g., he.tiff) inside this directory.")
            
            if not os.path.exists(source_labels_path):
                raise FileNotFoundError(f"Labels NPZ file not found: {source_labels_path}")
            if os.path.isdir(source_labels_path):
                raise IsADirectoryError(f"Labels NPZ path is a directory, not a file: {source_labels_path}. Please select the .npz file inside this directory.")
            
            if not os.path.exists(source_adata_path):
                raise FileNotFoundError(f"AnnData file not found: {source_adata_path}")
            if os.path.isdir(source_adata_path):
                raise IsADirectoryError(f"AnnData path is a directory, not a file: {source_adata_path}. Please select the .h5ad file inside this directory.")

            requested_config = {
                "source_he_path": source_he_path,
                "source_labels_path": source_labels_path,
                "source_adata_path": source_adata_path,
                "obsm_key": obsm_key,
                "tile_h": tile_h,
                "tile_w": tile_w,
                "stride_h": stride_h,
                "stride_w": stride_w,
                "stage_to_local": stage_to_local,
                "stage_root": stage_root if stage_to_local else None,
                "tile_workers": tile_workers,
                "single_tile_mode": single_tile_mode,
            }
            self._state["single_tile_mode"] = single_tile_mode

            existing_config = self._state.get("dataset_config")
            if existing_config and not self._config_paths_valid(existing_config):
                LOGGER.warning("Clearing invalid cached dataset config.")
                existing_config = None
                self._state["dataset_config"] = None

            if existing_config and all(existing_config.get(k) == requested_config.get(k) for k in requested_config):
                LOGGER.info("Dataset already loaded; reusing cached context.")
                ctx = self._ensure_context()
                tiles = self.tiles
                updated_config = dict(existing_config)
                updated_config.update(
                    {
                        "tile_cache_size": tile_cache_size,
                        "warm_cache": warm_cache,
                        "warm_cache_tiles": warm_cache_tiles,
                        "stage_to_local": stage_to_local,
                        "stage_root": stage_root,
                        "tile_workers": tile_workers,
                        "single_tile_mode": single_tile_mode,
                    }
                )
                self._state["dataset_config"] = updated_config
                self._state["single_tile_mode"] = single_tile_mode
                payload = {
                    "status": "ok",
                    "tile_count": len(tiles),
                    "tiles": [tile.to_dict() for tile in tiles],
                    "obsm_key": ctx.obsm_key,
                    "available_obsm": ctx.available_obsm,
                    "obs_columns": ctx.obs_columns,
                    "gene_count": len(ctx.gene_names),
                    "genes_preview": ctx.gene_names[:512],
                    "dataset_id": self.dataset_id,
                    "shape": ctx.shape,
                    **updated_config,
                }
                payload.update(
                    {
                        "tile_cache_size": tile_cache_size,
                        "stage_to_local": stage_to_local,
                        "stage_root": stage_root,
                        "warm_cache": warm_cache,
                        "warm_cache_tiles": warm_cache_tiles,
                        "single_tile_mode": single_tile_mode,
                    }
                )
                # Convert paths to appropriate format for response (relative for web server, absolute for standalone)
                if "he_path" in payload:
                    payload["he_path"] = self._normalize_path_for_response(payload["he_path"])
                self._start_cache_warmer(enabled=warm_cache, tile_limit=warm_cache_tiles, warm_params=warm_params)
                payload["cache_warm_status"] = self._cache_warm_status()
                return self._json_response(payload)

            self._cleanup_staged_dir()
            effective_paths = {
                "he_path": source_he_path,
                "labels_path": source_labels_path,
                "adata_path": source_adata_path,
            }
            staging_info: Dict[str, object] = {
                "enabled": stage_to_local,
                "stage_root": stage_root,
            }
            if stage_to_local and stage_root:
                staged_paths, staging_details = self._stage_dataset_files(
                    he_path=source_he_path,
                    labels_path=source_labels_path,
                    adata_path=source_adata_path,
                    stage_root=stage_root,
                )
                effective_paths.update(staged_paths)
                staging_info.update(staging_details)

            ctx = B2CContext(
                he_image_path=effective_paths["he_path"],
                labels_npz_path=effective_paths["labels_path"],
                adata_path=effective_paths["adata_path"],
                obsm_key=obsm_key,
            )
            tiles = make_tiles(ctx, tile_h=tile_h, tile_w=tile_w, stride_h=stride_h, stride_w=stride_w)

            self.context = ctx
            self.tiles = tiles
            self.tile_cache.clear()
            self.dataset_id = os.path.basename(effective_paths["adata_path"])

            dataset_config = {
                "he_path": effective_paths["he_path"],
                "labels_path": effective_paths["labels_path"],
                "adata_path": effective_paths["adata_path"],
                "obsm_key": ctx.obsm_key,
                "tile_h": tile_h,
                "tile_w": tile_w,
                "stride_h": stride_h,
                "stride_w": stride_w,
                "source_he_path": source_he_path,
                "source_labels_path": source_labels_path,
                "source_adata_path": source_adata_path,
                "stage_to_local": stage_to_local,
                "stage_root": stage_root,
                "staging": staging_info,
                "tile_cache_size": tile_cache_size,
                "warm_cache": warm_cache,
                "warm_cache_tiles": warm_cache_tiles,
                "tile_workers": tile_workers,
                "single_tile_mode": single_tile_mode,
            }

            payload = {
                "status": "ok",
                "tile_count": len(tiles),
                "tiles": [tile.to_dict() for tile in tiles],
                "obsm_key": ctx.obsm_key,
                "available_obsm": ctx.available_obsm,
                "obs_columns": ctx.obs_columns,
                "gene_count": len(ctx.gene_names),
                "genes_preview": ctx.gene_names[:512],
                "dataset_id": self.dataset_id,
                "shape": ctx.shape,
                **dataset_config,
            }
            # Convert paths to appropriate format for response (relative for web server, absolute for standalone)
            if "he_path" in payload:
                payload["he_path"] = self._normalize_path_for_response(payload["he_path"])
            self._state["dataset_config"] = dataset_config
            try:
                with open(self.state_path, "w") as fh:
                    json.dump(dataset_config, fh)
            except Exception as exc:
                LOGGER.warning("Failed to persist dataset config: %s", exc)
            LOGGER.info("Dataset loaded successfully: %d tiles", len(tiles))
            self._start_cache_warmer(enabled=warm_cache, tile_limit=warm_cache_tiles, warm_params=warm_params)
            payload["cache_warm_status"] = self._cache_warm_status()
            return self._json_response(payload)
        except Exception as exc:
            return self._error_from_exception(exc, "Failed to load dataset")
        finally:
            self._state["loading"] = False

    def describe_dataset(self, params: Dict):
        try:
            ctx = self._ensure_context()
            config = self._state.get("dataset_config") or {}
            payload = {
                "dataset_id": self.dataset_id,
                "shape": ctx.shape,
                "tile_count": len(self.tiles),
                "obsm_key": ctx.obsm_key,
                "available_obsm": ctx.available_obsm,
                "obs_columns": ctx.obs_columns,
                "gene_count": len(ctx.gene_names),
                "tile_cache_size": self._state.get("tile_cache_size"),
                "stage_to_local": config.get("stage_to_local"),
                "stage_root": config.get("stage_root"),
                "warm_cache": config.get("warm_cache"),
                "warm_cache_tiles": config.get("warm_cache_tiles"),
                "tile_workers": config.get("tile_workers") or self._state.get("worker_count"),
                "single_tile_mode": config.get("single_tile_mode"),
                "cache_warm_status": self._cache_warm_status(),
            }
            return self._json_response(payload)
        except Exception as exc:
            return self._error_from_exception(exc, "Failed to describe dataset")

    def list_tiles(self, params: Dict):
        try:
            self._ensure_context()
            return self._json_response({"tiles": [tile.to_dict() for tile in self.tiles]})
        except Exception as exc:
            return self._error_from_exception(exc, "Failed to list tiles")

    def list_obs_columns(self, params: Dict):
        try:
            ctx = self._ensure_context()
            cols = ctx.obs_columns
            prefer = []
            fallback = []
            for col in cols:
                series = ctx.obs_vector(col)
                if series.dtype.kind in {"O", "U", "S"} or str(series.dtype).startswith("category"):
                    prefer.append(col)
                else:
                    fallback.append(col)
            ordered = prefer + [c for c in fallback if c not in prefer]
            return self._json_response({"columns": ordered})
        except Exception as exc:
            return self._error_from_exception(exc, "Failed to list observation columns")

    def describe_obs_column(self, params: Dict):
        try:
            ctx = self._ensure_context()
            obs_col = params.get("obs_col")
            if not obs_col:
                raise ValueError("Provide 'obs_col'.")
            meta = ctx.obs_metadata(obs_col)
            payload = {
                "obs_col": obs_col,
                "categories": meta.get("categories", []),
                "color_map": meta.get("color_map", {}),
                "category_limit_hit": meta.get("category_limit_hit", False),
            }
            return self._json_response(payload)
        except Exception as exc:
            return self._error_from_exception(exc, "Failed to describe observation column")

    def list_genes(self, params: Dict):
        try:
            ctx = self._ensure_context()
            limit = int(params.get("limit", 0) or 0)
            genes = ctx.gene_names
            if limit > 0:
                genes = genes[:limit]
            return self._json_response({"genes": genes})
        except Exception as exc:
            return self._error_from_exception(exc, "Failed to list genes")

    #  Overlay core

    def _tile_entry(
        self,
        tile_id: int,
        *,
        b2c_mode: str,
        max_bin_distance: float,
        mpp: float,
        bin_um: float,
        volume_ratio: float,
        pad_factor: int = 2,
        cache_key: Optional[Tuple] = None,
    ) -> Dict:
        """
        Get or compute tile geometry (polygons, outlines, centroids).
        
        This function is PURE - it depends ONLY on geometry parameters.
        Color parameters have no effect here.
        
        Uses GeometryCache for caching instead of the old TileCache.
        """
        ctx = self._ensure_context()
        if tile_id < 0 or tile_id >= len(self.tiles):
            raise ValueError(f"Tile id {tile_id} out of range.")
        tile = self.tiles[tile_id]
        
        # Get geometry cache
        geometry_cache = self._state.get("geometry_cache")
        if geometry_cache is None:
            geometry_cache = GeometryCache(max_items=100)
            self._state["geometry_cache"] = geometry_cache
        
        # Create geometry cache key (ONLY geometry parameters!)
        geom_key = geometry_cache._make_key(
            self.dataset_id,
            tile_id,
            b2c_mode,
            max_bin_distance,
            mpp,
            bin_um,
            volume_ratio,
            pad_factor,
        )
        
        # Check geometry cache
        cached = geometry_cache.get(geom_key)
        if cached is not None:
            LOGGER.debug(f"GeometryCache HIT for tile {tile_id}")
            return cached
        
        LOGGER.info(f"GeometryCache MISS for tile {tile_id}, computing geometry...")
        
        import time
        t0 = time.time()
        
        # Cache miss - do expensive geometry computation
        _, lab_raw = ctx.crop_dense(tile.r0, tile.r1, tile.c0, tile.c1)
        t1 = time.time()
        LOGGER.info(f"Tile {tile_id}: crop_dense took {(t1-t0)*1000:.1f}ms")
        
        # Only do expensive label expansion if actually needed
        if b2c_mode == "none":
            # Skip expansion entirely - huge time savings!
            lab_exp = lab_raw  # Use raw labels directly
            dist_px = 0
            t2 = time.time()
            LOGGER.info(f"Tile {tile_id}: SKIPPED _expand_labels_tile (b2c_mode=none)")
        else:
            he_crop, lab_exp, dist_px = _expand_labels_tile(
                ctx,
                tile,
                mode=b2c_mode,
                max_bin_distance=max_bin_distance,
                mpp=mpp,
                bin_um=bin_um,
                volume_ratio=volume_ratio,
                pad_factor=pad_factor,
            )
            t2 = time.time()
            LOGGER.info(f"Tile {tile_id}: _expand_labels_tile took {(t2-t1)*1000:.1f}ms")

        centroid_idx, rows_abs, cols_abs, local_rc = _centroids_for_tile(ctx, tile)
        labels_at_centroids = lab_exp[local_rc[:, 0], local_rc[:, 1]] if local_rc.size else np.empty((0,), dtype=np.int32)
        t3 = time.time()
        LOGGER.info(f"Tile {tile_id}: _centroids_for_tile took {(t3-t2)*1000:.1f}ms, found {len(centroid_idx)} centroids")

        import cv2
        from scipy import ndimage
        
        def generate_polygons_fast(lab_array, tile_r0, tile_c0):
            """Generate polygons for all labels using bounding box + optimizations.
            
            Optimizations:
            1. Bounding box extraction (only process small regions)
            2. cv2.approxPolyDP for polygon simplification
            3. Integer coordinates (smaller JSON)
            """
            slices = ndimage.find_objects(lab_array.astype(np.int32))
            
            polygons = {}
            
            for lbl_idx, bbox in enumerate(slices):
                if bbox is None:
                    continue
                lbl = lbl_idx + 1
                
                row_slice, col_slice = bbox
                sub_array = lab_array[row_slice, col_slice]
                
                # Optimized mask creation
                mask = (sub_array == lbl).astype(np.uint8) * 255
                
                # findContours with CHAIN_APPROX_SIMPLE already reduces points
                contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
                
                poly_list = []
                offset_x = col_slice.start + tile_c0
                offset_y = row_slice.start + tile_r0
                
                for contour in contours:
                    if len(contour) < 3:
                        continue
                    
                    # Simplify polygon with approxPolyDP (epsilon=2.0 for faster browser rendering)
                    epsilon = 2.0
                    simplified = cv2.approxPolyDP(contour, epsilon, True)
                    
                    if len(simplified) < 3:
                        continue
                    
                    # Reshape and add offsets
                    pts = simplified.reshape(-1, 2)
                    pts[:, 0] += offset_x
                    pts[:, 1] += offset_y
                    
                    # Use integers (smaller JSON, faster)
                    polygon = [[int(x), int(y)] for x, y in pts]
                    
                    # Close polygon if needed
                    if polygon and polygon[0] != polygon[-1]:
                        polygon.append(polygon[0])
                    
                    if polygon:
                        poly_list.append(polygon)
                
                if poly_list:
                    polygons[lbl] = poly_list
            
            return polygons
        
        r0, c0 = tile.r0, tile.c0
        
        # Only compute expanded polygons if expansion is enabled
        need_expanded = (b2c_mode != "none")
        
        if need_expanded:
            polygons_exp = generate_polygons_fast(lab_exp, r0, c0)
            t4 = time.time()
            LOGGER.info(f"Tile {tile_id}: _polygons_from_labels(exp) VECTORIZED took {(t4-t3)*1000:.1f}ms, {len(polygons_exp)} labels")
        else:
            polygons_exp = {}
            t4 = time.time()
            LOGGER.info(f"Tile {tile_id}: SKIPPED expanded polygons (b2c_mode=none)")
        
        polygons_raw = generate_polygons_fast(lab_raw, r0, c0)
        t5 = time.time()
        LOGGER.info(f"Tile {tile_id}: _polygons_from_labels(raw) VECTORIZED took {(t5-t4)*1000:.1f}ms, {len(polygons_raw)} labels")

        # DERIVE outline paths from already-computed polygons (instead of recomputing)
        # This reuses the polygon data we just generated - instant!
        outline_exp = []
        if need_expanded:
            for lbl, poly_list in polygons_exp.items():
                outline_exp.extend(poly_list)
        t6 = time.time()
        LOGGER.info(f"Tile {tile_id}: outline_exp DERIVED from polygons ({len(outline_exp)} paths)")
        
        outline_raw = []
        for lbl, poly_list in polygons_raw.items():
            outline_raw.extend(poly_list)
        t7 = time.time()
        LOGGER.info(f"Tile {tile_id}: outline_raw DERIVED from polygons ({len(outline_raw)} paths)")
        
        # Get unique labels for per-label outlines (skip pre-computation - use on-demand)
        unique_labels = np.unique(labels_at_centroids)
        unique_labels = unique_labels[unique_labels > 0].tolist()
        
        # Pre-compute per-label outlines too for selected labels mode
        per_label_nuclei_outlines = {}
        per_label_expanded_outlines = {}
        
        if unique_labels:
            # Use the already-generated polygons as outlines (same geometry)
            for lbl in unique_labels:
                if lbl in polygons_raw:
                    per_label_nuclei_outlines[lbl] = polygons_raw[lbl]
                if lbl in polygons_exp:
                    per_label_expanded_outlines[lbl] = polygons_exp[lbl]
        
        t9 = time.time()
        LOGGER.info(f"Tile {tile_id}: Per-label outlines from polygons ({len(unique_labels)} labels)")
        LOGGER.info(f"Tile {tile_id}: TOTAL geometry computation took {(t9-t0)*1000:.1f}ms")

        entry = {
            "tile": tile,
            "lab_exp": lab_exp,
            "lab_raw": lab_raw,
            "dist_px": dist_px,
            "polygons_exp": polygons_exp,
            "polygons_raw": polygons_raw,
            "outline_exp": outline_exp,
            "outline_raw": outline_raw,
            "centroid_indices": centroid_idx,
            "centroid_rows": rows_abs,
            "centroid_cols": cols_abs,
            "labels_at_centroids": labels_at_centroids,
            "geometry_key": geom_key,  # NEW: Store for ColorCache linking
            # PRE-COMPUTED per-label outlines (for instant selected-only mode)
            "per_label_nuclei_outlines": per_label_nuclei_outlines,
            "per_label_expanded_outlines": per_label_expanded_outlines,
        }
        
        # Cache in geometry cache
        geometry_cache.set(geom_key, entry)
        
        # Also cache in old TileCache for backward compatibility with cache warmer
        old_cache_key = cache_key or self._tile_cache_key(
            tile_id, b2c_mode=b2c_mode, max_bin_distance=max_bin_distance,
            mpp=mpp, bin_um=bin_um, volume_ratio=volume_ratio, pad_factor=pad_factor
        )
        self.tile_cache.set(old_cache_key, entry)
        
        LOGGER.info(f"Geometry cached for tile {tile_id}")
        return entry

    def get_overlay(self, params: Dict):
        """
        Main API endpoint for getting overlay data.
        
        Flow:
        1. Parse parameters (geometry, color, visual)
        2. Get geometry from GeometryCache (via _tile_entry)
        3. Check ColorCache for colored overlay
        4. If cache miss, compute colors (fast: <50ms)
        5. Return compressed response
        """
        try:
            if self._state.get("loading"):
                return self._error_response(409, "Dataset is still loading; try again in a moment.")
            
            # Parse overlay type and tile
            overlay_type = params.get("overlay_type", "gene")
            tile_id = int(params.get("tile_id", 0))
            
            # Geometry parameters (affect polygon shapes)
            b2c_mode = params.get("b2c_mode", "fixed")
            max_bin_distance = float(params.get("max_bin_distance", 2.0) or 2.0)
            mpp = float(params.get("mpp", 0.3) or 0.3)
            bin_um = float(params.get("bin_um", 2.0) or 2.0)
            volume_ratio = float(params.get("volume_ratio", 4.0) or 4.0)
            pad_factor = int(params.get("pad_factor", 2) or 2)
            
            # Color parameters (affect which cells shown and their colors)
            gene = params.get("gene") or params.get("genes")
            obs_col = params.get("obs_col")
            category = params.get("category")
            color_mode = params.get("color_mode", "gradient")
            gradient_color = params.get("gradient_color")
            expr_quantile_param = params.get("expr_quantile")
            expr_quantile = float(expr_quantile_param) if expr_quantile_param not in (None, "") else None
            
            # Visual parameters (client-side only, not used in backend caching)
            overlay_alpha = float(params.get("overlay_alpha", 0.5) or 0.5)
            render_mode = params.get("render_mode", "fill")
            stroke_width = float(params.get("stroke_width", 1.0) or 1.0)
            show_outlines = _as_bool(params.get("show_outlines", True))
            all_expanded_outline = _as_bool(params.get("all_expanded_outline", False))
            all_nuclei_outline = _as_bool(params.get("all_nuclei_outline", False))
            nuclei_outline_color = params.get("nuclei_outline_color", "#000000")
            nuclei_outline_alpha = float(params.get("nuclei_outline_alpha", 0.6) or 0.6)
            
            # Legacy parameter
            include_geometry = _as_bool(params.get("include_geometry", True))

            # STEP 1: Get geometry (uses GeometryCache via _tile_entry)
            entry_future = self._schedule_tile_job(
                tile_id,
                b2c_mode=b2c_mode,
                max_bin_distance=max_bin_distance,
                mpp=mpp,
                bin_um=bin_um,
                volume_ratio=volume_ratio,
                pad_factor=pad_factor,
            )
            geometry_entry = entry_future.result()
            geometry_key = geometry_entry.get("geometry_key")
            
            # STEP 2: Check color cache
            import time
            t0 = time.time()
            
            color_cache = self._state.get("color_cache")
            if color_cache is None:
                color_cache = ColorCache(max_items=200)
                self._state["color_cache"] = color_cache
            
            color_key = color_cache._make_key(
                self.dataset_id,
                tile_id,
                overlay_type,
                geometry_key,
                gene=gene,
                obs_col=obs_col,
                category=category,
                color_mode=color_mode,
                gradient_color=gradient_color,
                expr_quantile=expr_quantile,
            )
            
            cached_colors = color_cache.get(color_key)
            if cached_colors is not None:
                LOGGER.info(f"ColorCache HIT for tile {tile_id}, overlay_type={overlay_type}")
                colored_overlay = cached_colors
            else:
                LOGGER.info(f"ColorCache MISS for tile {tile_id}, computing colors...")
                t1 = time.time()
                
                # Compute colors (fast: <50ms)
                if overlay_type == "gene":
                    colored_overlay = self._prepare_gene_overlay(geometry_entry, params, include_geometry=include_geometry)
                elif overlay_type == "observation":
                    colored_overlay = self._prepare_obs_overlay(geometry_entry, params, include_geometry=include_geometry)
                else:
                    raise ValueError(f"Unknown overlay_type '{overlay_type}'")
                
                t2 = time.time()
                LOGGER.info(f"Color computation took {(t2-t1)*1000:.1f}ms")
                
                # Cache the colored overlay
                color_cache.set(color_key, colored_overlay)
                LOGGER.info(f"ColorCache stored for tile {tile_id}")

            # STEP 3: Build final response with all parameters
            payload = {
                "status": "ok",
                "overlay_type": overlay_type,
                "tile": geometry_entry["tile"].to_dict(),
                "dist_px": geometry_entry["dist_px"],
                
                # Geometry parameters (for client-side geometry hash)
                "b2c_mode": b2c_mode,
                "max_bin_distance": max_bin_distance,
                "mpp": mpp,
                "bin_um": bin_um,
                "volume_ratio": volume_ratio,
                "pad_factor": pad_factor,
                
                # Visual parameters (for client-side rendering)
                "overlay_alpha": overlay_alpha,
                "render_mode": render_mode,
                "stroke_width": stroke_width,
                "show_outlines": show_outlines,
                "all_expanded_outline": all_expanded_outline,
                "all_nuclei_outline": all_nuclei_outline,
                "nuclei_outline_color": nuclei_outline_color,
                "nuclei_outline_alpha": nuclei_outline_alpha,
                
                # Color parameters (for client-side color hash)
                "color_mode": color_mode,
                "gradient_color": gradient_color,
                "expr_quantile": expr_quantile,
                
                # Legacy
                "geometry_included": include_geometry,
                
                # Colored overlay data
                **colored_overlay
            }

            # STEP 4: Return compressed response
            return self._json_response_compressed(payload)
            
        except Exception as exc:
            return self._error_from_exception(exc, "Failed to build overlay")

    def _prepare_gene_overlay(self, entry: Dict, params: Dict, *, include_geometry: bool) -> Dict:
        ctx = self._ensure_context()
        genes_value = params.get("genes") or params.get("gene")
        if not genes_value:
            raise ValueError("Provide 'gene' or 'genes' for gene overlay.")
        if isinstance(genes_value, str):
            genes = [g.strip() for g in genes_value.split(",") if g.strip()]
        else:
            genes = list(genes_value)

        if not genes:
            raise ValueError("No gene names provided.")

        color_mode = params.get("color_mode", "gradient")
        gene_color = params.get("gene_color", "#ff6b6b") or "#ff6b6b"
        gradient_color = params.get("gradient_color")
        render_mode = params.get("render_mode", "fill")
        expr_quantile = params.get("expr_quantile")
        expr_quantile = float(expr_quantile) if expr_quantile not in (None, "") else None
        top_n = int(params.get("top_n", 0) or 0)
        overlay_alpha = float(params.get("overlay_alpha", 0.5) or 0.5)
        show_centroids = str(params.get("show_centroids", "false")).lower() in ("1", "true", "yes", "on")
        highlight_color = params.get("highlight_color", "#39ff14")
        highlight_width = float(params.get("highlight_width", 2.0) or 2.0)
        cmap_name = params.get("cmap_name", "viridis")
        vmin_param = params.get("vmin")
        vmax_param = params.get("vmax")
        vmin = float(vmin_param) if vmin_param not in (None, "") else None
        vmax = float(vmax_param) if vmax_param not in (None, "") else None
        all_expanded_outline = str(params.get("all_expanded_outline", "false")).lower() in ("1", "true", "yes", "on")
        all_nuclei_outline = str(params.get("all_nuclei_outline", "false")).lower() in ("1", "true", "yes", "on")

        try:
            solid_rgba = colors.to_rgba(gene_color)
        except ValueError:
            gene_color = "#ff6b6b"
            solid_rgba = colors.to_rgba(gene_color)

        custom_cmap = None
        if color_mode == "gradient" and gradient_color:
            try:
                target_rgba = colors.to_rgba(gradient_color)
            except ValueError:
                target_rgba = colors.to_rgba("#4285f4")
            custom_cmap = colors.LinearSegmentedColormap.from_list(
                "bin2cell_custom_gradient",
                [(0.0, (1.0, 1.0, 1.0, 0.0)), (1.0, target_rgba)],
            )
        else:
            target_rgba = None

        lab_exp: np.ndarray = entry["lab_exp"]
        polygons_exp: Dict[int, List[List[List[float]]]] = entry["polygons_exp"]
        polygons_raw: Dict[int, List[List[List[float]]]] = entry["polygons_raw"]
        centroid_indices: np.ndarray = entry["centroid_indices"]
        labels_at_centroids: np.ndarray = entry["labels_at_centroids"]

        overlays = []
        centroid_payload = []

        for gene in genes:
            values = ctx.gene_vector(gene)
            tile_values = values[centroid_indices] if centroid_indices.size else np.empty((0,), dtype=float)

            selected_mask = np.ones_like(tile_values, dtype=bool)
            threshold = None
            # Validate expr_quantile is in valid range [0, 1]
            if expr_quantile is not None and tile_values.size:
                try:
                    q = float(expr_quantile)
                    if 0 <= q <= 1:
                        threshold = float(np.quantile(tile_values, q))
                        selected_mask &= tile_values > threshold
                    else:
                        LOGGER.warning(f"expr_quantile {q} out of range [0,1], ignoring")
                except (ValueError, TypeError) as e:
                    LOGGER.warning(f"Invalid expr_quantile value: {expr_quantile}, ignoring")

            selected_indices = centroid_indices[selected_mask]
            selected_values = tile_values[selected_mask]
            selected_labels = labels_at_centroids[selected_mask]

            label_expr: Dict[int, float] = {}
            for lbl, val in zip(selected_labels, selected_values):
                if lbl <= 0:
                    continue
                if (lbl not in label_expr) or (val > label_expr[lbl]):
                    label_expr[int(lbl)] = float(val)

            if top_n and len(label_expr) > top_n:
                top_labels = sorted(label_expr.keys(), key=lambda l: label_expr[l], reverse=True)[:top_n]
                label_expr = {lbl: label_expr[lbl] for lbl in top_labels}

            if not label_expr:
                overlays.append(
                    {
                        "gene": gene,
                        "features": [],
                        "legend": {
                            "type": "empty",
                            "gene": gene,
                            "message": "No labels passed filters in this tile.",
                        },
                    }
                )
                continue

            values_array = np.array(list(label_expr.values()), dtype=float)
            vmin_local = float(values_array.min()) if vmin is None else vmin
            vmax_local = float(values_array.max()) if vmax is None else vmax
            if vmax_local == vmin_local:
                vmax_local = vmin_local + 1e-6
            norm = colors.Normalize(vmin=vmin_local, vmax=vmax_local)

            # Use pre-computed polygons from geometry cache
            # (on-demand generation removed - polygons now always pre-computed)
            feature_polygons = polygons_exp if all_expanded_outline else polygons_raw
            
            # Generate per-label outlines for selected cells (if needed for selected-only mode)
            per_label_nuclei_outlines: Dict[int, List] = {}
            per_label_expanded_outlines: Dict[int, List] = {}
            if include_geometry:
                # Get labels that will be shown
                shown_labels = set(label_expr.keys())
                
                # Use PRE-CACHED per-label outlines from geometry (instant lookup!)
                cached_nuclei = entry.get("per_label_nuclei_outlines", {})
                cached_expanded = entry.get("per_label_expanded_outlines", {})
                
                # Filter to only shown labels (instant, just dict lookups)
                per_label_nuclei_outlines = {lbl: cached_nuclei.get(lbl, []) for lbl in shown_labels if lbl in cached_nuclei}
                per_label_expanded_outlines = {lbl: cached_expanded.get(lbl, []) for lbl in shown_labels if lbl in cached_expanded}
                
                LOGGER.debug(f"Using cached outlines: {len(per_label_nuclei_outlines)} nuclei, {len(per_label_expanded_outlines)} expanded")
            
            features = []
            for lbl, expr_value in label_expr.items():
                polygons = feature_polygons.get(lbl) or polygons_raw.get(lbl) or polygons_exp.get(lbl)
                if color_mode == "solid":
                    # Use solid gene color for all expressing cells
                    fill_color = _rgba_to_css(solid_rgba, alpha_override=overlay_alpha)
                    stroke_color = gene_color
                else:
                    color_fraction = float(norm(expr_value))
                    cmap_source = custom_cmap or cmap_name
                    fill_color = _colormap_sample(cmap_source, color_fraction, alpha=overlay_alpha)
                    stroke_color = _colormap_sample(cmap_source, color_fraction, alpha=1.0)

                feature = {
                    "label": int(lbl),
                    "polygons": polygons if include_geometry else [],
                    "fill": fill_color if render_mode == "fill" else None,
                    "stroke": stroke_color,
                    "stroke_width": highlight_width,
                    "value": float(expr_value),
                }
                
                # Add per-feature outlines for selected-only mode
                if include_geometry:
                    feature["nuclei_outline_paths"] = per_label_nuclei_outlines.get(lbl, [])
                    feature["expanded_outline_paths"] = per_label_expanded_outlines.get(lbl, [])
                
                if include_geometry or polygons:
                    features.append(feature)
                else:
                    # Preserve style/value for client-side geometry reuse
                    features.append({k: v for k, v in feature.items() if k not in ("polygons", "nuclei_outline_paths", "expanded_outline_paths")})

            legend = {
                "type": "binary" if color_mode == "binary" else ("solid" if color_mode == "solid" else "continuous"),
                "gene": gene,
                "color": gene_color if color_mode == "solid" else highlight_color,
                "overlay_alpha": overlay_alpha,
            }
            if color_mode == "gradient":
                cmap_source = custom_cmap or cmap_name
                legend.update(
                    {
                        "min": vmin_local,
                        "max": vmax_local,
                        "cmap": cmap_name,
                        "gradient": _sample_gradient(cmap_source),
                        "gradient_color": gradient_color,
                    }
                )

            overlays.append(
                {
                    "gene": gene,
                    "features": features,
                    "legend": legend,
                    "render_mode": render_mode,
                    "color_mode": color_mode,
                    "gene_color": gene_color,
                    "gradient_color": gradient_color,
                }
            )

            if show_centroids and centroid_indices.size:
                gene_centroids = []
                for lbl, idx_global, row, col, value in zip(
                    labels_at_centroids[selected_mask],
                    selected_indices,
                    entry["centroid_rows"][selected_mask],
                    entry["centroid_cols"][selected_mask],
                    selected_values,
                ):
                    if lbl <= 0:
                        continue
                    gene_centroids.append(
                        {
                            "index": int(idx_global),
                            "label": int(lbl),
                            "x": float(col),
                            "y": float(row),
                            "value": float(value),
                        }
                    )
                centroid_payload.append({"gene": gene, "points": gene_centroids})

        geometry_block = {}
        if include_geometry:
            geometry_block = {
                "polygons_exp": polygons_exp,
                "polygons_raw": polygons_raw,
                "outline_exp": entry["outline_exp"],
                "outline_raw": entry["outline_raw"],
                # Per-label outlines for selected-only mode
                "per_label_nuclei": entry.get("per_label_nuclei_outlines", {}),
                "per_label_expanded": entry.get("per_label_expanded_outlines", {}),
            }

        payload: Dict = {
            "overlay_type": "gene",
            "overlays": overlays,
            "all_expanded_outline": all_expanded_outline,
            "all_nuclei_outline": all_nuclei_outline,
            "expanded_outline": entry["outline_exp"] if (all_expanded_outline and include_geometry) else [],
            "nuclei_outline": entry["outline_raw"] if (all_nuclei_outline and include_geometry) else [],
        }
        if geometry_block:
            payload["geometry"] = geometry_block
        if show_centroids:
            payload["centroids"] = centroid_payload
        return payload

    def _prepare_obs_overlay(self, entry: Dict, params: Dict, *, include_geometry: bool) -> Dict:
        ctx = self._ensure_context()
        obs_col = params.get("obs_col")
        if not obs_col:
            raise ValueError("Provide 'obs_col' for observation overlay.")

        lab_exp: np.ndarray = entry["lab_exp"]
        polygons_exp: Dict[int, List[List[List[float]]]] = entry["polygons_exp"]
        polygons_raw: Dict[int, List[List[List[float]]]] = entry["polygons_raw"]
        centroid_indices: np.ndarray = entry["centroid_indices"]
        labels_at_centroids: np.ndarray = entry["labels_at_centroids"]

        category_filter = params.get("category")
        render_mode = params.get("render_mode", "fill")
        overlay_alpha = float(params.get("overlay_alpha", 0.5) or 0.5)
        cmap_name = params.get("cmap_name", "tab20")
        highlight_color = params.get("highlight_color", "#39ff14")
        highlight_width = float(params.get("highlight_width", 2.0) or 2.0)
        show_centroids = str(params.get("show_centroids", "false")).lower() in ("1", "true", "yes", "on")
        all_expanded_outline = str(params.get("all_expanded_outline", "true")).lower() in ("1", "true", "yes", "on")
        all_nuclei_outline = str(params.get("all_nuclei_outline", "false")).lower() in ("1", "true", "yes", "on")
        legend_outside = str(params.get("legend_outside", "true")).lower() in ("1", "true", "yes", "on")

        obs_vals = ctx.obs_vector(obs_col)
        tile_vals = obs_vals[centroid_indices] if centroid_indices.size else np.empty((0,), dtype=object)

        label_to_cat: Dict[int, str] = {}
        counts: Dict[int, Counter] = defaultdict(Counter)
        for lbl, val in zip(labels_at_centroids, tile_vals):
            if lbl <= 0:
                continue
            key = "" if val is None or (isinstance(val, float) and math.isnan(val)) else str(val)
            counts[int(lbl)][key] += 1
        for lbl, cnt in counts.items():
            if not cnt:
                continue
            label_to_cat[lbl] = cnt.most_common(1)[0][0]

        if not label_to_cat:
            return {
                "overlay_type": "observation",
                "obs_col": obs_col,
                "features": [],
                "legend": {"type": "empty", "message": "No categories found in tile."},
            }

        meta = ctx.obs_metadata(obs_col)
        ordered_categories = list(meta.get("categories") or [])
        predefined_colors = dict(meta.get("color_map") or {})

        if category_filter:
            allowed_categories = {category_filter}
        else:
            allowed_categories = set(label_to_cat.values())

        if ordered_categories:
            categories_sorted = [cat for cat in ordered_categories if cat in allowed_categories]
            remainder = [cat for cat in allowed_categories if cat not in categories_sorted]
            categories_sorted.extend(sorted(remainder))
        else:
            categories_sorted = sorted(allowed_categories)
        if not categories_sorted and allowed_categories:
            categories_sorted = sorted(allowed_categories)

        cmap_obj = cm.get_cmap(cmap_name, max(1, len(categories_sorted) or 1))
        denom = max(1, len(categories_sorted) - 1)
        category_styles: Dict[str, Dict[str, Optional[str]]] = {}
        for idx, cat in enumerate(categories_sorted):
            base_color = predefined_colors.get(cat)
            if base_color:
                fill_css = _color_to_css(base_color, alpha_override=overlay_alpha) if render_mode == "fill" else None
                stroke_css = _color_to_css(base_color, alpha_override=1.0)
                legend_color = base_color
            else:
                fraction = idx / denom if denom else 0.0
                fill_css = _colormap_sample(cmap_obj, fraction, alpha=overlay_alpha) if render_mode == "fill" else None
                stroke_css = _colormap_sample(cmap_obj, fraction, alpha=1.0)
                legend_color = stroke_css
            category_styles[cat] = {
                "fill": fill_css,
                "stroke": stroke_css,
                "legend": legend_color,
            }

        # Use pre-computed polygons from geometry cache
        # (on-demand generation removed - polygons now always pre-computed)
        feature_polygons = polygons_exp if all_expanded_outline else polygons_raw
        
        # Generate per-label outlines for selected cells (if needed for selected-only mode)
        per_label_nuclei_outlines: Dict[int, List] = {}
        per_label_expanded_outlines: Dict[int, List] = {}
        if include_geometry:
            # Get labels that will be shown
            shown_labels = set(label_to_cat.keys())
            
            # Use PRE-CACHED per-label outlines from geometry (instant lookup!)
            cached_nuclei = entry.get("per_label_nuclei_outlines", {})
            cached_expanded = entry.get("per_label_expanded_outlines", {})
            
            # Filter to only shown labels (instant, just dict lookups)
            per_label_nuclei_outlines = {lbl: cached_nuclei.get(lbl, []) for lbl in shown_labels if lbl in cached_nuclei}
            per_label_expanded_outlines = {lbl: cached_expanded.get(lbl, []) for lbl in shown_labels if lbl in cached_expanded}
            
            LOGGER.debug(f"[OBS] Using cached outlines: {len(per_label_nuclei_outlines)} nuclei, {len(per_label_expanded_outlines)} expanded")
        
        features = []
        for lbl, cat in label_to_cat.items():
            if category_filter and cat != category_filter:
                continue
            polygons = feature_polygons.get(lbl) or polygons_raw.get(lbl) or polygons_exp.get(lbl)
            style = category_styles.get(cat)
            if not style:
                continue
            fill_color = style["fill"] if render_mode == "fill" else None
            stroke_color = style["stroke"] or highlight_color

            feature_entry = {
                "label": int(lbl),
                "category": cat,
                "polygons": polygons if include_geometry else [],
                "fill": fill_color,
                "stroke": stroke_color,
                "stroke_width": highlight_width,
            }
            
            # Add per-feature outlines for selected-only mode
            if include_geometry:
                feature_entry["nuclei_outline_paths"] = per_label_nuclei_outlines.get(lbl, [])
                feature_entry["expanded_outline_paths"] = per_label_expanded_outlines.get(lbl, [])
            
            if include_geometry or polygons:
                features.append(feature_entry)
            else:
                features.append({k: v for k, v in feature_entry.items() if k != "polygons"})

        legend = {
            "type": "categorical",
            "obs_col": obs_col,
            "items": [
                {"label": cat, "color": (category_styles.get(cat) or {}).get("legend", highlight_color)}
                for cat in categories_sorted
            ],
            "legend_outside": legend_outside,
        }

        geometry_block = {}
        if include_geometry:
            geometry_block = {
                "polygons_exp": polygons_exp,
                "polygons_raw": polygons_raw,
                "outline_exp": entry["outline_exp"],
                "outline_raw": entry["outline_raw"],
                # Per-label outlines for selected-only mode
                "per_label_nuclei": entry.get("per_label_nuclei_outlines", {}),
                "per_label_expanded": entry.get("per_label_expanded_outlines", {}),
            }

        payload: Dict = {
            "overlay_type": "observation",
            "obs_col": obs_col,
            "category_filter": category_filter,
            "features": features,
            "legend": legend,
            "render_mode": render_mode,
            "all_expanded_outline": all_expanded_outline,
            "all_nuclei_outline": all_nuclei_outline,
            "expanded_outline": entry["outline_exp"] if (all_expanded_outline and include_geometry) else [],
            "nuclei_outline": entry["outline_raw"] if (all_nuclei_outline and include_geometry) else [],
        }

        if geometry_block:
            payload["geometry"] = geometry_block

        if show_centroids and centroid_indices.size:
            points = []
            for lbl, idx_global, row, col, obs_val in zip(
                labels_at_centroids,
                centroid_indices,
                entry["centroid_rows"],
                entry["centroid_cols"],
                tile_vals,
            ):
                if lbl <= 0:
                    continue
                label_str = label_to_cat.get(int(lbl))
                points.append(
                    {
                        "index": int(idx_global),
                        "label": int(lbl),
                        "category": label_str,
                        "x": float(col),
                        "y": float(row),
                    }
                )
            payload["centroids"] = [{"points": points}]

        return payload

    def export_overlay(self, params: Dict):
        try:
            response = self.get_overlay(params)
            if response.status_code and response.status_code >= 400:
                return response
            payload = json.loads(response.get_data(as_text=True))

            tile_id = payload["tile"]["id"]
            overlay_type = payload["overlay_type"]
            name = params.get("name")
            if not name:
                suffix = overlay_type
                if overlay_type == "gene" and payload.get("overlays"):
                    suffix = "_".join(ov["gene"] for ov in payload["overlays"])
                if overlay_type == "observation":
                    suffix = payload.get("obs_col", "observation")
                name = f"tile{tile_id}_{suffix}"

            safe_name = "".join(ch if ch.isalnum() or ch in ("_", "-", ".") else "_" for ch in name)
            out_path = os.path.join(self.out_dir, f"{safe_name}.geojson")

            features = []
            if overlay_type == "gene":
                for gene_overlay in payload.get("overlays", []):
                    gene_name = gene_overlay.get("gene")
                    for feature in gene_overlay.get("features", []):
                        geom = _feature_polygons_to_geojson(feature["polygons"])
                        props = {
                            "tile_id": tile_id,
                            "overlay_type": "gene",
                            "gene": gene_name,
                            "value": feature.get("value"),
                        }
                        features.append({"type": "Feature", "geometry": geom, "properties": props})
            else:
                obs_col = payload.get("obs_col")
                for feature in payload.get("features", []):
                    geom = _feature_polygons_to_geojson(feature["polygons"])
                    props = {
                        "tile_id": tile_id,
                        "overlay_type": "observation",
                        "obs_col": obs_col,
                        "category": feature.get("category"),
                    }
                    features.append({"type": "Feature", "geometry": geom, "properties": props})

            geojson = {"type": "FeatureCollection", "features": features}
            with open(out_path, "w") as fh:
                json.dump(geojson, fh)

            return self._json_response({"status": "ok", "path": out_path})
        except Exception as exc:
            return self._error_from_exception(exc, "Failed to export overlay")

    def save_preset(self, params: Dict):
        try:
            name = params.get("name")
            if not name:
                raise ValueError("Preset requires 'name'.")
            config = params.get("config")
            if not isinstance(config, dict):
                raise ValueError("Preset requires 'config' JSON object.")
            self.presets[name] = config
            self._save_presets()
            return self._json_response({"status": "ok", "name": name})
        except Exception as exc:
            return self._error_from_exception(exc, "Failed to save preset")

    def list_presets(self, params: Dict):
        try:
            return self._json_response({"presets": self.presets})
        except Exception as exc:
            return self._error_from_exception(exc, "Failed to list presets")

    def delete_preset(self, params: Dict):
        try:
            name = params.get("name")
            if not name or name not in self.presets:
                raise KeyError("Preset not found.")
            del self.presets[name]
            self._save_presets()
            return self._json_response({"status": "ok", "deleted": name})
        except Exception as exc:
            return self._error_from_exception(exc, "Failed to delete preset")


def _feature_polygons_to_geojson(polygons: List[List[List[float]]]) -> Dict:
    if not polygons:
        return {"type": "GeometryCollection", "geometries": []}
    if len(polygons) == 1:
        return {"type": "Polygon", "coordinates": [polygons[0]]}
    return {"type": "MultiPolygon", "coordinates": [[poly] for poly in polygons]}


def _json_default(obj):
    if isinstance(obj, Tile):
        return obj.to_dict()
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")
class FileDialogHelper(QObject):
    def __init__(self) -> None:
        super().__init__()
        self._result: Tuple[str, str] = ("", "")

    @Slot(str, str, str)
    def _open_dialog(self, caption: str, start_dir: str, filters: str) -> None:
        self._result = QFileDialog.getOpenFileName(None, caption, start_dir, filters)

    def get_open_file_name(self, caption: str, start_dir: str, filters: str) -> Tuple[str, str]:
        self._result = ("", "")
        QMetaObject.invokeMethod(
            self,
            "_open_dialog",
            Qt.BlockingQueuedConnection,
            Q_ARG(str, caption),
            Q_ARG(str, start_dir),
            Q_ARG(str, filters),
        )
        return self._result
