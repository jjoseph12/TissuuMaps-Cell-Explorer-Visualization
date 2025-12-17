var Bin2CellExplorer;
Bin2CellExplorer = {
  name: "Cell Explorer",
  parameters: {
    _sec_data: { label: "Dataset", type: "section", collapsed: false },
    he_path: { label: "H&E image (.tif/.tiff)", type: "text", default: "" },
    he_path_browse: { label: "Browse H&E image…", type: "button" },
    labels_path: { label: "Label matrix (.npz)", type: "text", default: "" },
    labels_path_browse: { label: "Browse label matrix…", type: "button" },
    adata_path: { label: "AnnData (.h5ad)", type: "text", default: "" },
    adata_path_browse: { label: "Browse AnnData…", type: "button" },
    obsm_key: { label: "obsm coord key", type: "select", default: "spatial_cropped_150_buffer" },
    tile_h: { label: "Tile height (px)", type: "number", default: 1500 },
    tile_w: { label: "Tile width (px)", type: "number", default: 1500 },
    stride_h: { label: "Stride height (px, optional)", type: "number", default: "" },
    stride_w: { label: "Stride width (px, optional)", type: "number", default: "" },
    single_tile_mode: { label: "Single-tile focus (disable caching)", type: "checkbox", default: false },
    stage_to_local: { label: "Stage dataset to node-local scratch", type: "checkbox", default: false },
    stage_root: { label: "Scratch override (optional)", type: "text", default: "" },
    tile_cache_size: { label: "Tile cache size", type: "number", default: 32 },
    warm_cache: { label: "Warm cache after load", type: "checkbox", default: true },
    warm_cache_tiles: { label: "Tiles to prewarm (0 = all)", type: "number", default: 32 },
    tile_workers: { label: "Tile worker threads", type: "number", default: 4 },
    load_dataset_btn: { label: "Load dataset", type: "button" },

    _sec_overlay: { label: "Overlay", type: "section", collapsed: false },
    tile_id: { label: "Tile ID", type: "select", default: "" },
    overlay_type: { label: "Overlay type", type: "select", default: "gene", options: ["gene", "observation"] },
    genes: { label: "Gene(s) (comma separated)", type: "text", default: "COL1A1" },
    obs_col: { label: "Observation column", type: "select", default: "" },
    category: { label: "Category filter (optional)", type: "select", default: "" },
    render_mode: { label: "Render mode", type: "select", default: "fill", options: ["fill", "outline"] },
    color_mode: { label: "Gene color mode", type: "select", default: "gradient", options: ["gradient", "solid"] },
    gradient_color: { label: "Gradient color (optional)", type: "text", default: "#4285f4", attributes: { type: "color" } },
    gene_color: { label: "Gene color (solid mode)", type: "text", default: "#ff6b6b", attributes: { type: "color" } },
    expr_quantile: { label: "Expr. quantile (0-1)", type: "number", default: "" },
    top_n: { label: "Top N labels", type: "number", default: "" },
    b2c_mode: { label: "Expand mode", type: "select", default: "none", options: ["none", "fixed", "volume_ratio"] },
    max_bin_distance: { label: "Max bin distance", type: "number", default: 2.0 },
    mpp: { label: "Microns per pixel", type: "number", default: 0.3 },
    bin_um: { label: "Bin size (µm)", type: "number", default: 2.0 },
    volume_ratio: { label: "Volume ratio", type: "number", default: 4.0 },
    render_mode: { label: "Render mode", type: "select", default: "fill", options: ["fill", "outline"] },
    overlay_alpha: { label: "Overlay opacity (0-1)", type: "number", default: 0.7 },
    highlight_color: { label: "Outline color", type: "text", default: "#39ff14", attributes: { type: "color" } },
    highlight_width: { label: "Outline width", type: "number", default: 2.0 },
    all_expanded_outline: { label: "Show expanded outlines (all cells)", type: "checkbox", default: false },
    expanded_outlines_selected: { label: "Show expanded outlines (selected only)", type: "checkbox", default: false },
    all_nuclei_outline: { label: "Show nuclei outlines (all cells)", type: "checkbox", default: false },
    nuclei_outlines_selected: { label: "Show nuclei outlines (selected only)", type: "checkbox", default: true },
    nuclei_outline_color: { label: "Nuclei outline color", type: "text", default: "#000000", attributes: { type: "color" } },
    nuclei_outline_alpha: { label: "Nuclei outline opacity (0-1)", type: "number", default: 0.6 },
    apply_overlay_btn: { label: "Update overlay", type: "button" },




    _sec_export: { label: "Export & Presets", type: "section", collapsed: true },
    export_name: { label: "Export name (optional)", type: "text", default: "" },
    export_geojson_btn: { label: "Export overlay to GeoJSON", type: "button" },
    save_preset_name: { label: "Preset name", type: "text", default: "" },
    save_preset_btn: { label: "Save preset", type: "button" },
    preset_select: { label: "Presets", type: "select", default: "" },
    load_preset_btn: { label: "Apply preset", type: "button" },
    delete_preset_btn: { label: "Delete preset", type: "button" }
  }
};

Bin2CellExplorer.state = {
  datasetLoaded: false,
  tiles: [],
  obsColumns: [],
  obsMetadata: {},
  obsMetadataRequests: {},
  obsCategorySelections: {},
  pendingCategoryValue: "",
  selectedObsCol: null,
  genesPreview: [],
  overlays: null,
  layers: {},
  presets: {},
  slideShape: null,
  tileOverviewCanvas: null,
  tileOverviewScale: null,
  selectedTileId: null,
  filePickerModalId: null,
  filePickerField: null,
  lastOverlayPayload: null,
  renderPayload: null,
  geometryCache: {},
  hiddenGenes: {},
  hiddenCategories: {},
  canvasOverlay: null,
  canvasCtx: null,
  viewerHooksInstalled: false,
  overlayRefreshToken: null,
  cacheWarmStatus: null,
  singleTileMode: false,
  // NEW: Track last request parameters for smart updates
  lastRequest: {
    geometry: null,  // { tile_id, b2c_mode, mpp, bin_um, volume_ratio, ... }
    color: null,     // { overlay_type, gene, obs_col, category, color_mode, ... }
    visual: null     // { overlay_alpha, render_mode, all_expanded_outline, ... }
  }
};

Bin2CellExplorer._prefixes = ["Bin2CellExplorer_", "CellExplorer_"];
Bin2CellExplorer._findElement = function (name) {
  for (let i = 0; i < Bin2CellExplorer._prefixes.length; i += 1) {
    const el = document.getElementById(Bin2CellExplorer._prefixes[i] + name);
    if (el) return el;
  }
  return null;
};
Bin2CellExplorer._domId = function (name) {
  for (let i = 0; i < Bin2CellExplorer._prefixes.length; i += 1) {
    const id = Bin2CellExplorer._prefixes[i] + name;
    if (document.getElementById(id)) return id;
  }
  // Prefer the filename-based prefix if nothing is mounted yet.
  return "CellExplorer_" + name;
};

Bin2CellExplorer.get = function (name) {
  let element = Bin2CellExplorer._findElement(name);
  if (!element) {
    return "";
  }
  if (element.matches && element.matches("input, select, textarea")) {
    // element is already the actual form control
  } else {
    const nested = element.querySelector("input, select, textarea");
    if (nested) {
      element = nested;
    }
  }
  if (element.tagName === "SELECT") {
    return element.value || "";
  }
  if (element.type === "checkbox") {
    return element.checked;
  }
  if (element.type === "number" || element.type === "range") {
    const raw = element.value;
    return raw === "" ? "" : Number(raw);
  }
  const value = element.value !== undefined ? element.value : element.textContent;
  return value == null ? "" : value;
};

Bin2CellExplorer.set = function (name, value) {
  let element = Bin2CellExplorer._findElement(name);
  if (!element) return;
  if (!(element.matches && element.matches("input, select, textarea"))) {
    const nested = element.querySelector("input, select, textarea");
    if (nested) {
      element = nested;
    }
  }
  if (!element) return;
  if (element.tagName === "SELECT") {
    element.value = value == null ? "" : String(value);
    element.dispatchEvent(new Event("change", { bubbles: true }));
    return;
  }
  if (element.type === "checkbox") {
    element.checked = !!value;
    element.dispatchEvent(new Event("change", { bubbles: true }));
    return;
  }
  element.value = value == null ? "" : value;
  element.dispatchEvent(new Event("input", { bubbles: true }));
};

Bin2CellExplorer.init = function (container) {
  Bin2CellExplorer.state.container = container;
  container.classList.add("bin2cell-explorer-panel");

  const legend = document.createElement("div");
  legend.id = "Bin2CellExplorer_legend";
  legend.className = "bin2cell-legend mt-2";
  container.appendChild(legend);

  const status = document.createElement("div");
  status.id = "Bin2CellExplorer_status";
  status.className = "bin2cell-status mt-2 small text-muted";
  container.appendChild(status);

  const overviewBlock = document.createElement("div");
  overviewBlock.id = "Bin2CellExplorer_tile_overview_block";
  overviewBlock.className = "tile-overview mt-3";
  const overviewTitle = document.createElement("div");
  overviewTitle.className = "small fw-bold mb-1";
  overviewTitle.textContent = "Tile overview";
  const overviewCanvas = document.createElement("canvas");
  overviewCanvas.id = "Bin2CellExplorer_tile_overview";
  overviewCanvas.width = 320;
  overviewCanvas.height = 320;
  overviewCanvas.style.border = "1px solid rgba(0,0,0,0.1)";
  overviewCanvas.style.borderRadius = "12px";
  overviewCanvas.style.background = "linear-gradient(135deg, #fdfbfb, #ebedee)";
  overviewCanvas.style.boxShadow = "0 6px 18px rgba(15,23,42,0.15)";
  overviewCanvas.style.cursor = "pointer";
  const overviewHint = document.createElement("div");
  overviewHint.className = "small text-muted mt-1";
  overviewHint.id = "Bin2CellExplorer_tile_overview_hint";
  overviewHint.textContent = "Load dataset to display tiles";
  overviewBlock.appendChild(overviewTitle);
  overviewBlock.appendChild(overviewCanvas);
  overviewBlock.appendChild(overviewHint);
  container.appendChild(overviewBlock);

  Bin2CellExplorer.state.tileOverviewCanvas = overviewCanvas;
  overviewCanvas.addEventListener("click", Bin2CellExplorer.onTileOverviewClick);

  interfaceUtils.alert("Bin2Cell Explorer loaded");

  // Listen for file selection messages from web file browser iframe
  window.addEventListener("message", function (event) {
    if (event.data && event.data.type === "b2c_file_selected") {
      Bin2CellExplorer.onWebFilePicked(event.data.field, event.data.path);
    }
  });

  Bin2CellExplorer.toggleOverlayInputs("gene");
  Bin2CellExplorer.updateDatasetMode();

  // Set initial visibility for expansion and render params
  Bin2CellExplorer.toggleExpansionParams("none");  // Default: hide all expansion params
  Bin2CellExplorer.toggleRenderParams("fill");      // Default mode
};

/**
 * Set up dynamic show/hide of parameters based on dropdown values
 */
Bin2CellExplorer.setupConditionalVisibility = function () {
  // Wait a bit for TissUUmaps to create the UI elements
  setTimeout(function () {
    // Helper to get the row element for a parameter
    const getParamRow = function (paramName) {
      // TissUUmaps typically creates inputs with IDs like "Bin2CellExplorer_paramName"
      const prefixes = Bin2CellExplorer._prefixes;
      for (let i = 0; i < prefixes.length; i++) {
        const elem = document.getElementById(prefixes[i] + paramName);
        if (elem && elem.closest) {
          return elem.closest('.form-group, .row, tr');
        }
      }
      return null;
    };

    // Update visibility based on b2c_mode
    const updateExpansionVisibility = function () {
      const mode = Bin2CellExplorer.get("b2c_mode") || "fixed";

      const maxBinRow = getParamRow("max_bin_distance");
      const mppRow = getParamRow("mpp");
      const binUmRow = getParamRow("bin_um");
      const volumeRow = getParamRow("volume_ratio");

      // Show/hide based on mode
      if (mode === "none") {
        // Hide all expansion params
        if (maxBinRow) maxBinRow.style.display = "none";
        if (mppRow) mppRow.style.display = "none";
        if (binUmRow) binUmRow.style.display = "none";
        if (volumeRow) volumeRow.style.display = "none";
      } else if (mode === "fixed") {
        // Show fixed-mode params, hide volume params
        if (maxBinRow) maxBinRow.style.display = "";
        if (mppRow) mppRow.style.display = "";
        if (binUmRow) binUmRow.style.display = "";
        if (volumeRow) volumeRow.style.display = "none";
      } else if (mode === "volume_ratio") {
        // Hide max_bin_distance, show others
        if (maxBinRow) maxBinRow.style.display = "none";
        if (mppRow) mppRow.style.display = "";
        if (binUmRow) binUmRow.style.display = "";
        if (volumeRow) volumeRow.style.display = "";
      }
    };

    // Update visibility based on render_mode
    const updateRenderVisibility = function () {
      const mode = Bin2CellExplorer.get("render_mode") || "fill";

      const colorRow = getParamRow("highlight_color");
      const widthRow = getParamRow("highlight_width");

      // Only show highlight options when render_mode is "outline"
      if (mode === "outline") {
        if (colorRow) colorRow.style.display = "";
        if (widthRow) widthRow.style.display = "";
      } else {
        if (colorRow) colorRow.style.display = "none";
        if (widthRow) widthRow.style.display = "none";
      }
    };

    // Initial update
    updateExpansionVisibility();
    updateRenderVisibility();

    // Add change listeners
    const b2cSelect = document.getElementById("Bin2CellExplorer_b2c_mode");
    if (b2cSelect) {
      b2cSelect.addEventListener("change", updateExpansionVisibility);
    }

    const renderSelect = document.getElementById("Bin2CellExplorer_render_mode");
    if (renderSelect) {
      renderSelect.addEventListener("change", updateRenderVisibility);
    }

    console.log("[CellExplorer] ✨ Conditional visibility initialized");
  }, 500);  // Wait 500ms for UI to be created
};

Bin2CellExplorer.inputTrigger = function (inputName) {
  switch (inputName) {
    case "load_dataset_btn":
      Bin2CellExplorer.loadDataset();
      break;
    case "he_path_browse":
      Bin2CellExplorer.browseForFile("he_path");
      break;
    case "labels_path_browse":
      Bin2CellExplorer.browseForFile("labels_path");
      break;
    case "adata_path_browse":
      Bin2CellExplorer.browseForFile("adata_path");
      break;
    case "single_tile_mode":
      Bin2CellExplorer.updateDatasetMode();
      break;
    case "overlay_type":
      Bin2CellExplorer.toggleOverlayInputs(Bin2CellExplorer.get("overlay_type"));
      break;
    case "color_mode":
      Bin2CellExplorer.updateGeneColorControls();
      break;
    case "b2c_mode":
      Bin2CellExplorer.toggleExpansionParams(Bin2CellExplorer.get("b2c_mode"));
      break;
    case "render_mode":
      Bin2CellExplorer.toggleRenderParams(Bin2CellExplorer.get("render_mode"));
      break;
    case "apply_overlay_btn":
      Bin2CellExplorer.requestOverlay();
      break;
    case "export_geojson_btn":
      Bin2CellExplorer.exportOverlay();
      break;
    case "save_preset_btn":
      Bin2CellExplorer.savePreset();
      break;
    case "load_preset_btn":
      Bin2CellExplorer.applySelectedPreset();
      break;
    case "delete_preset_btn":
      Bin2CellExplorer.deleteSelectedPreset();
      break;
    case "preset_select":
      Bin2CellExplorer.populatePresetPreview();
      break;
    default:
      break;
  }
};

Bin2CellExplorer.toggleOverlayInputs = function (mode) {
  const geneIds = ["genes", "color_mode", "gradient_color", "gene_color", "expr_quantile", "top_n"];
  const obsIds = ["obs_col", "category"];
  geneIds.forEach((id) => Bin2CellExplorer.toggleParam(id, mode === "gene"));
  obsIds.forEach((id) => Bin2CellExplorer.toggleParam(id, mode === "observation"));
  if (mode === "gene") {
    Bin2CellExplorer.updateGeneColorControls();
  } else if (mode === "observation" && Bin2CellExplorer.state.datasetLoaded) {
    Bin2CellExplorer.populateCategorySelect(Bin2CellExplorer.get("obs_col") || "");
  }
};

Bin2CellExplorer.toggleExpansionParams = function (mode) {
  // Show/hide expansion parameters based on b2c_mode
  const fixedOnlyIds = ["max_bin_distance"];
  const volumeOnlyIds = ["volume_ratio"];
  const commonIds = ["mpp", "bin_um"];

  if (mode === "none") {
    // Hide all expansion parameters
    fixedOnlyIds.forEach((id) => Bin2CellExplorer.toggleParam(id, false));
    volumeOnlyIds.forEach((id) => Bin2CellExplorer.toggleParam(id, false));
    commonIds.forEach((id) => Bin2CellExplorer.toggleParam(id, false));
  } else {
    // Show/hide based on mode
    fixedOnlyIds.forEach((id) => Bin2CellExplorer.toggleParam(id, mode === "fixed"));
    volumeOnlyIds.forEach((id) => Bin2CellExplorer.toggleParam(id, mode === "volume_ratio"));
    commonIds.forEach((id) => Bin2CellExplorer.toggleParam(id, true));
  }
};

Bin2CellExplorer.toggleRenderParams = function (mode) {
  // Show/hide highlight parameters based on render_mode
  const outlineOnlyIds = ["highlight_color", "highlight_width"];
  outlineOnlyIds.forEach((id) => Bin2CellExplorer.toggleParam(id, mode === "outline"));
};

Bin2CellExplorer.toggleParam = function (name, visible) {
  const element = Bin2CellExplorer._findElement(name);
  if (!element) return;
  let wrapper = element.closest(".form-group, .row, .input-group");
  if (!wrapper) {
    wrapper = element.parentElement;
  }
  if (wrapper) wrapper.style.display = visible ? "" : "none";
};

Bin2CellExplorer.ensureObject = function (payload) {
  if (!payload) return {};
  if (typeof payload === "string") {
    try {
      return JSON.parse(payload);
    } catch (err) {
      console.warn("Bin2CellExplorer: failed to parse payload", err);
      return {};
    }
  }
  return payload;
};

Bin2CellExplorer.setStatus = function (msg) {
  const status = document.getElementById("Bin2CellExplorer_status") || document.getElementById("CellExplorer_status");
  if (status) status.textContent = msg || "";
};

Bin2CellExplorer.updateGeneColorControls = function () {
  const mode = Bin2CellExplorer.get("color_mode");
  Bin2CellExplorer.toggleParam("gradient_color", mode === "gradient");
  Bin2CellExplorer.toggleParam("gene_color", mode === "solid");
};

Bin2CellExplorer.browseForFile = function (field) {
  const current = Bin2CellExplorer.get(field) || "";
  const payload = {
    field: field,
    current_path: current,
  };
  Bin2CellExplorer.api(
    "pick_file",
    payload,
    function (resp) {
      const data = Bin2CellExplorer.ensureObject(resp);
      if (data && data.status === "web") {
        Bin2CellExplorer.openWebFileBrowser(field, data);
        return;
      }
      if (!data || data.status === "cancelled") {
        Bin2CellExplorer.setStatus("File selection cancelled.");
        return;
      }
      if (data.status !== "ok" || !data.path) {
        Bin2CellExplorer.setStatus("File selection failed.");
        return;
      }
      // Set value robustly
      Bin2CellExplorer.set(field, data.path);
      const el = Bin2CellExplorer._findElement(field);
      if (el) {
        el.value = String(data.path);
        el.setAttribute("value", String(data.path));
      }
      console.log("DEBUG: Browse set " + field + " to:", data.path, "(post-set read:", Bin2CellExplorer.get(field), ")");
      Bin2CellExplorer.setStatus("Selected " + data.path);
    },
    Bin2CellExplorer.handleError
  );
};

Bin2CellExplorer.openWebFileBrowser = function (field, data) {
  const modalUID = "Bin2CellExplorer_filepicker";
  Bin2CellExplorer.state.filePickerField = field;
  Bin2CellExplorer.state.filePickerModalId = modalUID;

  // Use API call instead of direct iframe URL
  console.log("DEBUG: Opening web file browser for field:", field);
  Bin2CellExplorer.api(
    "filetree",
    {
      field: field,
      start_rel: data.start_rel || "",
    },
    function (response) {
      console.log("DEBUG: Got filetree response:", response);
      const data = Bin2CellExplorer.ensureObject(response);
      const iframe = document.createElement("iframe");
      iframe.srcdoc = data.html; // Use srcdoc to set HTML content directly
      iframe.width = "100%";
      iframe.height = "360px";
      iframe.style.border = "0";

      const content = document.createElement("div");
      content.appendChild(iframe);

      const closeBtn = HTMLElementUtils.createButton({
        extraAttributes: { class: "btn btn-primary mx-2" },
      });
      closeBtn.innerText = "Cancel";
      closeBtn.addEventListener("click", function () {
        Bin2CellExplorer.closeWebFileBrowser();
      });
      const buttons = document.createElement("div");
      buttons.appendChild(closeBtn);

      interfaceUtils.generateModal("Select file", content, buttons, modalUID);
      setTimeout(function () {
        $(`#${modalUID}_modal`).on("hidden.bs.modal", function () {
          Bin2CellExplorer.closeWebFileBrowser(true);
        });
      }, 0);
    },
    Bin2CellExplorer.handleError
  );
};

Bin2CellExplorer.closeWebFileBrowser = function (skipHide) {
  const modalUID = Bin2CellExplorer.state.filePickerModalId;
  if (!skipHide && modalUID) {
    $(`#${modalUID}_modal`).modal("hide");
  }
  Bin2CellExplorer.state.filePickerModalId = null;
  Bin2CellExplorer.state.filePickerField = null;
};

Bin2CellExplorer.onWebFilePicked = function (field, absolutePath) {
  if (!absolutePath || field !== Bin2CellExplorer.state.filePickerField) {
    return;
  }
  // Set value robustly
  Bin2CellExplorer.set(field, absolutePath);
  const el = Bin2CellExplorer._findElement(field);
  if (el) {
    el.value = String(absolutePath);
    el.setAttribute("value", String(absolutePath));
  }
  console.log("DEBUG: Web browse set " + field + " to:", absolutePath, "(post-set read:", Bin2CellExplorer.get(field), ")");
  Bin2CellExplorer.closeWebFileBrowser();
};

Bin2CellExplorer.loadDataset = function () {
  const he_path = Bin2CellExplorer.get("he_path");
  const labels_path = Bin2CellExplorer.get("labels_path");
  const adata_path = Bin2CellExplorer.get("adata_path");

  console.log("DEBUG: UI values - he_path:", he_path, "labels_path:", labels_path, "adata_path:", adata_path);

  // Front-end validation: require all three paths before calling backend
  if (!he_path || !labels_path || !adata_path) {
    const missing = [];
    if (!he_path) missing.push("H&E image");
    if (!labels_path) missing.push("labels matrix");
    if (!adata_path) missing.push("AnnData");
    const msg = "Please provide: " + missing.join(", ") + ". Use Browse or paste full paths.";
    console.error("Bin2CellExplorer:", msg);
    Bin2CellExplorer.setStatus(msg);
    interfaceUtils.alert(msg);
    return;
  }

  Bin2CellExplorer.setStatus("Loading dataset…");
  const payload = {
    he_path: he_path,
    labels_path: labels_path,
    adata_path: adata_path,
    obsm_key: Bin2CellExplorer.get("obsm_key"),
    tile_h: Bin2CellExplorer.get("tile_h"),
    tile_w: Bin2CellExplorer.get("tile_w"),
    stride_h: Bin2CellExplorer.get("stride_h"),
    stride_w: Bin2CellExplorer.get("stride_w"),
    single_tile_mode: Bin2CellExplorer.isChecked("single_tile_mode"),
    stage_to_local: Bin2CellExplorer.isChecked("stage_to_local"),
    stage_root: Bin2CellExplorer.get("stage_root"),
    tile_cache_size: Bin2CellExplorer.get("tile_cache_size"),
    warm_cache: Bin2CellExplorer.isChecked("warm_cache"),
    warm_cache_tiles: Bin2CellExplorer.get("warm_cache_tiles"),
    warm_b2c_mode: Bin2CellExplorer.get("b2c_mode"),
    warm_max_bin_distance: Bin2CellExplorer.get("max_bin_distance"),
    warm_mpp: Bin2CellExplorer.get("mpp"),
    warm_bin_um: Bin2CellExplorer.get("bin_um"),
    warm_volume_ratio: Bin2CellExplorer.get("volume_ratio"),
    tile_workers: Bin2CellExplorer.get("tile_workers")
  };
  Bin2CellExplorer.api(
    "load_dataset",
    payload,
    function (resp) {
      Bin2CellExplorer.onDatasetLoaded(Bin2CellExplorer.ensureObject(resp));
    },
    Bin2CellExplorer.handleError
  );
};

Bin2CellExplorer.onDatasetLoaded = function (data) {
  Bin2CellExplorer.state.datasetLoaded = true;
  Bin2CellExplorer.state.tiles = data.tiles || [];
  Bin2CellExplorer.state.obsColumns = data.obs_columns || [];
  Bin2CellExplorer.state.obsMetadata = {};
  Bin2CellExplorer.state.obsMetadataRequests = {};
  Bin2CellExplorer.state.obsCategorySelections = {};
  Bin2CellExplorer.state.pendingCategoryValue = "";
  Bin2CellExplorer.state.selectedObsCol = null;
  Bin2CellExplorer.state.genesPreview = data.genes_preview || [];
  Bin2CellExplorer.state.datasetId = data.dataset_id;
  Bin2CellExplorer.state.hePath = data.he_path;
  Bin2CellExplorer.state.slideShape = data.shape;
  Bin2CellExplorer.state.geometryCache = {};
  Bin2CellExplorer.state.hiddenGenes = {};
  Bin2CellExplorer.state.hiddenCategories = {};
  Bin2CellExplorer.state.renderPayload = null;

  // Update the input field with the normalized path (relative for web server, absolute for standalone)
  // This prevents TissUUmaps from trying to add the layer with the wrong path format
  if (data.he_path) {
    interfaceUtils.setValueForElement(Bin2CellExplorer._domId("he_path"), "value", data.he_path);
  }

  if (typeof data.tile_cache_size !== "undefined") {
    interfaceUtils.setValueForElement(Bin2CellExplorer._domId("tile_cache_size"), "value", data.tile_cache_size);
  }
  if (typeof data.stage_to_local !== "undefined") {
    Bin2CellExplorer.setCheckbox("stage_to_local", !!data.stage_to_local);
  }
  if (typeof data.stage_root !== "undefined" && data.stage_root !== null) {
    interfaceUtils.setValueForElement(Bin2CellExplorer._domId("stage_root"), "value", data.stage_root || "");
  }
  if (typeof data.warm_cache !== "undefined") {
    Bin2CellExplorer.setCheckbox("warm_cache", !!data.warm_cache);
  }
  if (typeof data.warm_cache_tiles !== "undefined") {
    interfaceUtils.setValueForElement(Bin2CellExplorer._domId("warm_cache_tiles"), "value", data.warm_cache_tiles);
  }
  if (typeof data.tile_workers !== "undefined") {
    interfaceUtils.setValueForElement(Bin2CellExplorer._domId("tile_workers"), "value", data.tile_workers);
  }
  if (typeof data.single_tile_mode !== "undefined") {
    Bin2CellExplorer.setCheckbox("single_tile_mode", !!data.single_tile_mode);
  }
  Bin2CellExplorer.updateDatasetMode();

  Bin2CellExplorer.populateSelect("tile_id", Bin2CellExplorer.state.tiles.map(function (tile) { return tile.id; }));
  Bin2CellExplorer.populateSelect("obs_col", Bin2CellExplorer.state.obsColumns);
  Bin2CellExplorer.populateSelect("category", [""]);
  Bin2CellExplorer.populateSelect("obsm_key", data.available_obsm || [], data.obsm_key);

  if (Bin2CellExplorer.state.tiles.length) {
    Bin2CellExplorer.set("tile_id", String(Bin2CellExplorer.state.tiles[0].id));
    Bin2CellExplorer.state.selectedTileId = Bin2CellExplorer.state.tiles[0].id;
  }
  if (Bin2CellExplorer.state.obsColumns.length) {
    const firstObs = Bin2CellExplorer.state.obsColumns[0];
    Bin2CellExplorer.set("obs_col", firstObs);
    Bin2CellExplorer.onObsColumnChange(firstObs);
  } else {
    Bin2CellExplorer.populateCategorySelect("");
  }

  Bin2CellExplorer.requestPresets();

  Bin2CellExplorer.renderTileOverview();
  Bin2CellExplorer.attachTileSelectListener();
  Bin2CellExplorer.attachObsSelectListener();
  Bin2CellExplorer.attachCategorySelectListener();

  Bin2CellExplorer.state.cacheWarmStatus = data.cache_warm_status || null;
  if (Bin2CellExplorer.state.cacheWarmStatus && Bin2CellExplorer.state.cacheWarmStatus.active) {
    Bin2CellExplorer.setStatus(
      "Dataset loaded (" +
      Bin2CellExplorer.state.tiles.length +
      " tiles). Cache warming " +
      Bin2CellExplorer.state.cacheWarmStatus.progress +
      "/" +
      (Bin2CellExplorer.state.cacheWarmStatus.total || "?") +
      "…"
    );
  } else {
    Bin2CellExplorer.setStatus("Dataset loaded (" + Bin2CellExplorer.state.tiles.length + " tiles)");
  }
};

Bin2CellExplorer.populateSelect = function (name, options, selected) {
  const domId = Bin2CellExplorer._domId(name);
  interfaceUtils.cleanSelect(domId);
  options.forEach(function (opt) {
    interfaceUtils.addSingleElementToSelect(domId, String(opt));
  });
  if (selected !== undefined && selected !== null && selected !== "") {
    interfaceUtils.setValueForElement(domId, "value", String(selected));
  }
  if (name === "category") {
    const selectEl = document.getElementById(domId);
    if (selectEl && selectEl.options.length) {
      for (let i = 0; i < selectEl.options.length; i += 1) {
        if (selectEl.options[i].value === "") {
          selectEl.options[i].text = selectEl.options[i].text || "All categories";
          if (!selectEl.options[i].text.trim()) {
            selectEl.options[i].text = "All categories";
          }
          break;
        }
      }
    }
  }
};

Bin2CellExplorer.attachTileSelectListener = function () {
  const select = Bin2CellExplorer._findElement("tile_id");
  if (!select || select.__bin2cell_tile_listener) return;
  select.addEventListener("change", function () {
    const value = parseInt(select.value, 10);
    if (!isNaN(value)) {
      Bin2CellExplorer.state.selectedTileId = value;
      Bin2CellExplorer.renderTileOverview();
      Bin2CellExplorer.panToTile(value, true);
    }
  });
  select.__bin2cell_tile_listener = true;
};

Bin2CellExplorer.attachObsSelectListener = function () {
  const select = Bin2CellExplorer._findElement("obs_col");
  if (!select || select.__bin2cell_obs_listener) return;
  select.addEventListener("change", function () {
    Bin2CellExplorer.onObsColumnChange(select.value || "");
  });
  select.__bin2cell_obs_listener = true;
};

Bin2CellExplorer.attachCategorySelectListener = function () {
  const select = Bin2CellExplorer._findElement("category");
  if (!select || select.__bin2cell_category_listener) return;
  select.addEventListener("change", function () {
    const col = Bin2CellExplorer.get("obs_col");
    if (col) {
      Bin2CellExplorer.state.obsCategorySelections[col] = select.value || "";
    }
  });
  select.__bin2cell_category_listener = true;
};

Bin2CellExplorer.onObsColumnChange = function (column, options) {
  if (!Bin2CellExplorer.state.datasetLoaded) return;
  const value = column || "";
  Bin2CellExplorer.state.selectedObsCol = value || null;
  const opts = options || {};
  if (!opts.keepCategory) {
    const cached = value ? Bin2CellExplorer.state.obsCategorySelections[value] : "";
    Bin2CellExplorer.state.pendingCategoryValue = cached || "";
    interfaceUtils.setValueForElement(Bin2CellExplorer._domId("category"), "value", "");
  }
  Bin2CellExplorer.ensureObsMetadata(value);
};

Bin2CellExplorer.ensureObsMetadata = function (column) {
  const col = column || "";
  Bin2CellExplorer.populateCategorySelect(col);
  if (!col) return;
  if (Bin2CellExplorer.state.obsMetadata[col]) {
    Bin2CellExplorer.populateCategorySelect(col);
    return;
  }
  if (Bin2CellExplorer.state.obsMetadataRequests[col]) {
    return;
  }
  Bin2CellExplorer.state.obsMetadataRequests[col] = true;
  Bin2CellExplorer.api(
    "describe_obs_column",
    { obs_col: col },
    function (resp) {
      delete Bin2CellExplorer.state.obsMetadataRequests[col];
      const data = Bin2CellExplorer.ensureObject(resp) || {};
      Bin2CellExplorer.state.obsMetadata[col] = data;
      Bin2CellExplorer.populateCategorySelect(col);
    },
    function (jqXHR, textStatus, errorThrown) {
      delete Bin2CellExplorer.state.obsMetadataRequests[col];
      Bin2CellExplorer.handleError(jqXHR, textStatus, errorThrown);
    }
  );
};

Bin2CellExplorer.populateCategorySelect = function (column) {
  const metadata = column ? Bin2CellExplorer.state.obsMetadata[column] : null;
  const categories = (metadata && Array.isArray(metadata.categories)) ? metadata.categories.slice() : [];
  const options = [""].concat(categories);
  let desired = Bin2CellExplorer.state.pendingCategoryValue;
  if (!desired) {
    desired = Bin2CellExplorer.get("category") || "";
  }
  if (options.indexOf(desired) === -1) {
    desired = "";
  }
  Bin2CellExplorer.populateSelect("category", options, desired);
  if (desired && Bin2CellExplorer.state.pendingCategoryValue === desired) {
    Bin2CellExplorer.state.pendingCategoryValue = "";
  }
  if (metadata && Array.isArray(metadata.categories) && column) {
    if (desired) {
      Bin2CellExplorer.state.obsCategorySelections[column] = desired;
    } else {
      delete Bin2CellExplorer.state.obsCategorySelections[column];
    }
  }
};

Bin2CellExplorer.requestOverlay = function () {
  if (!Bin2CellExplorer.state.datasetLoaded) {
    interfaceUtils.alert("Load a dataset first.");
    return;
  }

  const overlayType = Bin2CellExplorer.get("overlay_type") || "gene";
  const payload = Bin2CellExplorer.collectOverlayParams();
  payload.overlay_type = overlayType;
  const tileId = Number(payload.tile_id);

  // SMART UPDATE: Detect what changed
  const changeType = Bin2CellExplorer.detectChangeType(payload);

  if (changeType === 'NONE') {
    console.log('[CellExplorer] ⏭️  No changes detected, skipping update');
    return;
  }

  if (changeType === 'VISUAL') {
    console.log('[CellExplorer] ⚡ Visual-only change, updating instantly!');
    const startTime = performance.now();
    Bin2CellExplorer.updateVisualOnly(payload);
    const elapsed = performance.now() - startTime;
    Bin2CellExplorer.setStatus(`Overlay updated (${elapsed.toFixed(1)}ms)`);
    return;
  }

  // For COLOR or GEOMETRY changes, make backend request
  console.log(`[CellExplorer] 🔄 ${changeType} change, requesting from backend...`);
  payload.include_geometry = Bin2CellExplorer.shouldRequestGeometry(tileId, payload);

  Bin2CellExplorer.setStatus("Computing overlay…");
  Bin2CellExplorer.api(
    "get_overlay",
    payload,
    function (resp) {
      const data = Bin2CellExplorer.ensureObject(resp);
      Bin2CellExplorer.renderOverlay(data);
      Bin2CellExplorer.updateLastRequest(payload);  // Track successful request
      Bin2CellExplorer.setStatus("Overlay ready (" + data.overlay_type + ")");
    },
    Bin2CellExplorer.handleError
  );
};

/**
 * Update visual parameters only (instant, no backend request).
 * This is called when ONLY visual parameters changed.
 */
Bin2CellExplorer.updateVisualOnly = function (params) {
  const payload = Bin2CellExplorer.state.lastOverlayPayload;

  if (!payload) {
    console.warn('[CellExplorer] ⚠️  No cached overlay data, making full backend request');
    Bin2CellExplorer.state.lastRequest = { geometry: null, color: null, visual: null };
    return Bin2CellExplorer.requestOverlay();
  }

  // Update visual parameters in cached payload
  payload.overlay_alpha = parseFloat(params.overlay_alpha) || 0.7;
  payload.render_mode = params.render_mode || "fill";
  payload.highlight_color = params.highlight_color || "#39ff14";
  payload.highlight_width = parseFloat(params.highlight_width) || 2.0;
  payload.all_expanded_outline = params.all_expanded_outline;
  payload.expanded_outlines_selected = params.expanded_outlines_selected;
  payload.all_nuclei_outline = params.all_nuclei_outline;
  payload.nuclei_outlines_selected = params.nuclei_outlines_selected;
  payload.nuclei_outline_color = params.nuclei_outline_color || "#000000";
  payload.nuclei_outline_alpha = parseFloat(params.nuclei_outline_alpha) || 0.6;

  // Re-render canvas with updated parameters (fast!)
  Bin2CellExplorer.drawOverlayToCanvas(payload);
  Bin2CellExplorer.state.renderPayload = payload;

  // Update tracked visual params
  Bin2CellExplorer.updateLastRequest(params);

  console.log('[CellExplorer] ✨ Visual update complete');
};

Bin2CellExplorer.renderTileOverview = function () {
  const canvas = Bin2CellExplorer.state.tileOverviewCanvas;
  const hint = document.getElementById("Bin2CellExplorer_tile_overview_hint");
  if (!canvas) return;
  const ctx = canvas.getContext("2d");
  const tiles = Bin2CellExplorer.state.tiles || [];
  const shape = Bin2CellExplorer.state.slideShape;
  if (!tiles.length || !shape) {
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    if (hint) hint.textContent = "Load dataset to display tiles";
    Bin2CellExplorer.state.tileOverviewScale = null;
    return;
  }

  const slideHeight = shape[0];
  const slideWidth = shape[1];
  const margin = 16;
  const maxWidth = 320;
  const maxHeight = 320;
  const scale = Math.min((maxWidth - margin * 2) / slideWidth, (maxHeight - margin * 2) / slideHeight);
  const canvasWidth = Math.ceil(slideWidth * scale + margin * 2);
  const canvasHeight = Math.ceil(slideHeight * scale + margin * 2);
  canvas.width = canvasWidth;
  canvas.height = canvasHeight;
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  const gradient = ctx.createLinearGradient(0, 0, canvas.width, canvas.height);
  gradient.addColorStop(0, "#fdfbfb");
  gradient.addColorStop(1, "#ebedee");
  ctx.fillStyle = gradient;
  ctx.fillRect(0, 0, canvas.width, canvas.height);

  const selectedId = Bin2CellExplorer.state.selectedTileId;
  ctx.lineWidth = 1;
  ctx.strokeStyle = "rgba(60,60,60,0.5)";
  ctx.font = "10px sans-serif";
  ctx.textAlign = "center";
  ctx.textBaseline = "middle";

  tiles.forEach(function (tile) {
    const x = margin + tile.c0 * scale;
    const y = margin + tile.r0 * scale;
    const w = (tile.c1 - tile.c0) * scale;
    const h = (tile.r1 - tile.r0) * scale;
    if (tile.id === selectedId) {
      ctx.fillStyle = "rgba(255, 193, 7, 0.45)";
      ctx.fillRect(x, y, w, h);
    }
    ctx.strokeStyle = "rgba(60,60,60,0.6)";
    ctx.strokeRect(x, y, w, h);
    if (w > 10 && h > 10) {
      ctx.fillStyle = "#202124";
      ctx.fillText(String(tile.id), x + w / 2, y + h / 2);
    }
  });

  if (hint) {
    hint.textContent = "Click a tile to select and center";
  }
  Bin2CellExplorer.state.tileOverviewScale = { scale: scale, margin: margin };
};

Bin2CellExplorer.onTileOverviewClick = function (evt) {
  const canvas = Bin2CellExplorer.state.tileOverviewCanvas;
  const scaleInfo = Bin2CellExplorer.state.tileOverviewScale;
  if (!canvas || !scaleInfo || !Bin2CellExplorer.state.tiles.length) return;
  const rect = canvas.getBoundingClientRect();
  const x = (evt.clientX - rect.left);
  const y = (evt.clientY - rect.top);
  const slideX = (x - scaleInfo.margin) / scaleInfo.scale;
  const slideY = (y - scaleInfo.margin) / scaleInfo.scale;
  if (slideX < 0 || slideY < 0) return;
  const tile = Bin2CellExplorer.state.tiles.find(function (t) {
    return slideX >= t.c0 && slideX < t.c1 && slideY >= t.r0 && slideY < t.r1;
  });
  if (!tile) return;
  Bin2CellExplorer.state.selectedTileId = tile.id;
  Bin2CellExplorer.set("tile_id", String(tile.id));
  Bin2CellExplorer.renderTileOverview();
  Bin2CellExplorer.panToTile(tile.id, true);
};

Bin2CellExplorer.panToTile = function (tileId, animate) {
  const tile = Bin2CellExplorer.state.tiles.find(function (t) { return t.id === tileId; });
  const shape = Bin2CellExplorer.state.slideShape;
  if (!tile || !shape) return;
  const viewer = tmapp[tmapp["object_prefix"] + "_viewer"];
  if (!viewer || !viewer.world.getItemCount()) return;
  const image = viewer.world.getItemAt(0);
  const width = shape[1];
  const height = shape[0];
  const xNorm = tile.c0 / width;
  const yNorm = tile.r0 / height;
  const wNorm = (tile.c1 - tile.c0) / width;
  const hNorm = (tile.r1 - tile.r0) / height;
  const rect = new OpenSeadragon.Rect(xNorm, yNorm, wNorm, hNorm);
  viewer.viewport.fitBounds(rect, animate !== false);
};

Bin2CellExplorer.collectOverlayParams = function () {
  const params = {
    tile_id: Bin2CellExplorer.get("tile_id"),
    overlay_type: Bin2CellExplorer.get("overlay_type"),
    genes: Bin2CellExplorer.get("genes"),
    obs_col: Bin2CellExplorer.get("obs_col"),
    category: Bin2CellExplorer.get("category"),
    render_mode: Bin2CellExplorer.get("render_mode"),
    color_mode: Bin2CellExplorer.get("color_mode"),
    gradient_color: Bin2CellExplorer.get("gradient_color"),
    gene_color: Bin2CellExplorer.get("gene_color"),
    expr_quantile: Bin2CellExplorer.get("expr_quantile"),
    top_n: Bin2CellExplorer.get("top_n"),
    b2c_mode: Bin2CellExplorer.get("b2c_mode"),
    max_bin_distance: Bin2CellExplorer.get("max_bin_distance"),
    mpp: Bin2CellExplorer.get("mpp"),
    bin_um: Bin2CellExplorer.get("bin_um"),
    volume_ratio: Bin2CellExplorer.get("volume_ratio"),
    overlay_alpha: Bin2CellExplorer.get("overlay_alpha"),
    highlight_color: Bin2CellExplorer.get("highlight_color"),
    highlight_width: Bin2CellExplorer.get("highlight_width"),
    all_expanded_outline: Bin2CellExplorer.isChecked("all_expanded_outline"),
    expanded_outlines_selected: Bin2CellExplorer.isChecked("expanded_outlines_selected"),
    all_nuclei_outline: Bin2CellExplorer.isChecked("all_nuclei_outline"),
    nuclei_outlines_selected: Bin2CellExplorer.isChecked("nuclei_outlines_selected"),
    nuclei_outline_color: Bin2CellExplorer.get("nuclei_outline_color"),
    nuclei_outline_alpha: Bin2CellExplorer.get("nuclei_outline_alpha")
  };
  return params;
};

Bin2CellExplorer.geometrySignature = function (payload) {
  if (!payload) return "";
  const tileId = payload.tile_id != null ? Number(payload.tile_id) : -1;
  return [
    tileId,
    payload.b2c_mode,
    payload.max_bin_distance,
    payload.mpp,
    payload.bin_um,
    payload.volume_ratio,
    payload.pad_factor || 2
  ].join("|");
};

/**
 * Detect what type of parameters changed.
 * Returns: 'VISUAL' | 'COLOR' | 'GEOMETRY' | 'NONE'
 */
Bin2CellExplorer.detectChangeType = function (newParams) {
  const last = Bin2CellExplorer.state.lastRequest;

  // First request ever
  if (!last.geometry && !last.color && !last.visual) {
    return 'GEOMETRY';
  }

  // Check geometry parameters
  const geometryParams = ['tile_id', 'b2c_mode', 'max_bin_distance', 'mpp', 'bin_um', 'volume_ratio'];
  for (let i = 0; i < geometryParams.length; i++) {
    const param = geometryParams[i];
    if (String(newParams[param]) !== String(last.geometry[param])) {
      console.log(`[CellExplorer] Geometry changed: ${param} (${last.geometry[param]} → ${newParams[param]})`);
      return 'GEOMETRY';
    }
  }

  // Check color parameters
  const colorParams = ['overlay_type', 'genes', 'obs_col', 'category', 'color_mode', 'gradient_color', 'expr_quantile'];
  for (let i = 0; i < colorParams.length; i++) {
    const param = colorParams[i];
    if (String(newParams[param] || '') !== String(last.color[param] || '')) {
      console.log(`[CellExplorer] Color changed: ${param} (${last.color[param]} → ${newParams[param]})`);
      return 'COLOR';
    }
  }

  // Check visual parameters
  const visualParams = ['overlay_alpha', 'render_mode', 'highlight_color', 'highlight_width',
    'all_expanded_outline', 'expanded_outlines_selected', 'all_nuclei_outline', 'nuclei_outlines_selected',
    'nuclei_outline_color', 'nuclei_outline_alpha'];
  for (let i = 0; i < visualParams.length; i++) {
    const param = visualParams[i];
    if (String(newParams[param] || '') !== String(last.visual[param] || '')) {
      console.log(`[CellExplorer] Visual changed: ${param} (${last.visual[param]} → ${newParams[param]})`);
      return 'VISUAL';
    }
  }

  return 'NONE';  // No changes
};

/**
 * Update cached parameters after successful request.
 */
Bin2CellExplorer.updateLastRequest = function (params) {
  Bin2CellExplorer.state.lastRequest = {
    geometry: {
      tile_id: params.tile_id,
      b2c_mode: params.b2c_mode,
      max_bin_distance: params.max_bin_distance,
      mpp: params.mpp,
      bin_um: params.bin_um,
      volume_ratio: params.volume_ratio,
      pad_factor: params.pad_factor || 2
    },
    color: {
      overlay_type: params.overlay_type,
      genes: params.genes,
      obs_col: params.obs_col,
      category: params.category,
      color_mode: params.color_mode,
      gradient_color: params.gradient_color,
      expr_quantile: params.expr_quantile
    },
    visual: {
      overlay_alpha: params.overlay_alpha,
      render_mode: params.render_mode,
      highlight_color: params.highlight_color,
      highlight_width: params.highlight_width,
      all_expanded_outline: params.all_expanded_outline,
      expanded_outlines_selected: params.expanded_outlines_selected,
      all_nuclei_outline: params.all_nuclei_outline,
      nuclei_outlines_selected: params.nuclei_outlines_selected,
      nuclei_outline_color: params.nuclei_outline_color,
      nuclei_outline_alpha: params.nuclei_outline_alpha
    }
  };
};

Bin2CellExplorer.shouldRequestGeometry = function (tileId, params) {
  if (tileId == null || isNaN(tileId)) return true;
  const sig = Bin2CellExplorer.geometrySignature(params);
  const cached = Bin2CellExplorer.state.geometryCache[tileId];
  if (!cached) return true;
  return cached.signature !== sig;
};

Bin2CellExplorer.exportOverlay = function () {
  if (!Bin2CellExplorer.state.datasetLoaded) {
    interfaceUtils.alert("Load a dataset first.");
    return;
  }
  const params = Bin2CellExplorer.collectOverlayParams();
  params.name = Bin2CellExplorer.get("export_name");
  Bin2CellExplorer.api(
    "export_overlay",
    params,
    function (resp) {
      const data = Bin2CellExplorer.ensureObject(resp);
      Bin2CellExplorer.setStatus("GeoJSON written to " + data.path);
      interfaceUtils.alert("Overlay exported:\n" + data.path);
    },
    Bin2CellExplorer.handleError
  );
};

Bin2CellExplorer.collectPresetConfig = function () {
  return Bin2CellExplorer.collectOverlayParams();
};

Bin2CellExplorer.savePreset = function () {
  const name = (Bin2CellExplorer.get("save_preset_name") || "").trim();
  if (!name) {
    interfaceUtils.alert("Provide a preset name.");
    return;
  }
  const payload = {
    name: name,
    config: Bin2CellExplorer.collectPresetConfig()
  };
  Bin2CellExplorer.api(
    "save_preset",
    payload,
    function (resp) {
      const data = Bin2CellExplorer.ensureObject(resp);
      Bin2CellExplorer.setStatus("Preset saved.");
      Bin2CellExplorer.requestPresets(name);
    },
    Bin2CellExplorer.handleError
  );
};

Bin2CellExplorer.requestPresets = function (selectName) {
  Bin2CellExplorer.api(
    "list_presets",
    {},
    function (resp) {
      const data = Bin2CellExplorer.ensureObject(resp);
      Bin2CellExplorer.state.presets = data.presets || {};
      const names = Object.keys(Bin2CellExplorer.state.presets);
      Bin2CellExplorer.populateSelect("preset_select", [""].concat(names), selectName || "");
    },
    Bin2CellExplorer.handleError
  );
};

Bin2CellExplorer.applySelectedPreset = function () {
  const select = Bin2CellExplorer.get("preset_select");
  if (!select || !(select in Bin2CellExplorer.state.presets)) {
    interfaceUtils.alert("Select a preset first.");
    return;
  }
  const config = Bin2CellExplorer.state.presets[select];
  if (!config) return;
  let pendingCategory = null;
  let pendingObsCol = null;
  Object.keys(config).forEach(function (key) {
    if (key === "category") {
      pendingCategory = config[key] || "";
      return;
    }
    const domId = "Bin2CellExplorer_" + key;
    if (document.getElementById(domId)) {
      interfaceUtils.setValueForElement(domId, "value", config[key]);
      if (key === "obs_col") {
        pendingObsCol = config[key] || "";
      }
    }
  });
  Bin2CellExplorer.toggleOverlayInputs(config.overlay_type || Bin2CellExplorer.get("overlay_type"));
  if (pendingCategory !== null) {
    Bin2CellExplorer.state.pendingCategoryValue = pendingCategory;
  }
  if (pendingObsCol !== null) {
    Bin2CellExplorer.onObsColumnChange(pendingObsCol, { keepCategory: pendingCategory !== null });
  } else if (pendingCategory !== null) {
    Bin2CellExplorer.populateCategorySelect(Bin2CellExplorer.get("obs_col") || "");
  }
  Bin2CellExplorer.setStatus("Preset applied: " + select);
};

Bin2CellExplorer.populatePresetPreview = function () {
  const select = Bin2CellExplorer.get("preset_select");
  if (!select) return;
  const preset = Bin2CellExplorer.state.presets[select];
  if (!preset) return;
  Bin2CellExplorer.setStatus("Preset '" + select + "' ready to apply.");
};

Bin2CellExplorer.deleteSelectedPreset = function () {
  const select = Bin2CellExplorer.get("preset_select");
  if (!select || !(select in Bin2CellExplorer.state.presets)) {
    interfaceUtils.alert("Select a preset to delete.");
    return;
  }
  Bin2CellExplorer.api(
    "delete_preset",
    { name: select },
    function () {
      Bin2CellExplorer.setStatus("Preset deleted: " + select);
      Bin2CellExplorer.requestPresets();
    },
    Bin2CellExplorer.handleError
  );
};

Bin2CellExplorer.handleError = function (jqXHR, textStatus, errorThrown) {
  let message = "";
  if (jqXHR) {
    if (jqXHR.responseJSON) {
      message = jqXHR.responseJSON.message || JSON.stringify(jqXHR.responseJSON);
    } else if (jqXHR.responseText) {
      try {
        const parsed = JSON.parse(jqXHR.responseText);
        message = parsed.message || jqXHR.responseText;
      } catch (parseErr) {
        message = jqXHR.responseText;
      }
    }
  }
  if (!message && errorThrown) {
    message = errorThrown;
  }
  if (!message && textStatus) {
    message = textStatus;
  }
  if (!message) {
    message = "Unknown error";
  }
  Bin2CellExplorer.setStatus("Error: " + message);
  interfaceUtils.alert("Bin2Cell Explorer error:\n" + message);
};

Bin2CellExplorer.isChecked = function (name) {
  const el = Bin2CellExplorer._findElement(name);
  return !!(el && el.checked);
};

Bin2CellExplorer.setCheckbox = function (name, value) {
  const el = Bin2CellExplorer._findElement(name);
  if (el) {
    el.checked = !!value;
  }
};

Bin2CellExplorer.updateDatasetMode = function () {
  const single = Bin2CellExplorer.isChecked("single_tile_mode");
  Bin2CellExplorer.toggleDatasetControl("tile_cache_size", !single);
  Bin2CellExplorer.toggleDatasetControl("warm_cache", !single);
  Bin2CellExplorer.toggleDatasetControl("warm_cache_tiles", !single);
};

Bin2CellExplorer.toggleDatasetControl = function (name, enabled) {
  const el = Bin2CellExplorer._findElement(name);
  if (!el) return;
  el.disabled = !enabled;
  const wrapper = el.closest(".form-group, .row, .input-group");
  if (wrapper) {
    wrapper.style.opacity = enabled ? "" : "0.5";
  }
};

Bin2CellExplorer.ensureCanvasOverlay = function () {
  const viewer = tmapp[tmapp["object_prefix"] + "_viewer"];
  if (!viewer) return null;
  const container = viewer.container || viewer.canvas || viewer.element;
  if (!container) return null;
  if (!Bin2CellExplorer.state.canvasOverlay) {
    const canvas = document.createElement("canvas");
    canvas.id = "Bin2CellExplorer_canvas_overlay";
    canvas.style.position = "absolute";
    canvas.style.top = "0";
    canvas.style.left = "0";
    canvas.style.width = "100%";
    canvas.style.height = "100%";
    canvas.style.pointerEvents = "none";
    canvas.style.zIndex = 10;
    container.appendChild(canvas);
    Bin2CellExplorer.state.canvasOverlay = canvas;
    Bin2CellExplorer.state.canvasCtx = canvas.getContext("2d");
    // Hide legacy SVG layer if present to avoid double rendering
    if (Bin2CellExplorer.state.d3layer && Bin2CellExplorer.state.d3layer.remove) {
      Bin2CellExplorer.state.d3layer.remove();
      Bin2CellExplorer.state.d3layer = null;
    }
  }
  return Bin2CellExplorer.state.canvasOverlay;
};

Bin2CellExplorer.projectPoint = function (image, viewer, pt) {
  const vp = image.imageToViewportCoordinates(pt[0], pt[1], true);
  const px = viewer.viewport.pixelFromPoint(vp, true);
  return px;
};

Bin2CellExplorer.drawPolygonList = function (ctx, polygons, style, projectFn) {
  if (!polygons || !polygons.length) return;
  polygons.forEach(function (poly) {
    if (!poly || poly.length < 3) return;
    const path = new Path2D();
    poly.forEach(function (pt, idx) {
      const screen = projectFn(pt);
      if (idx === 0) {
        path.moveTo(screen.x, screen.y);
      } else {
        path.lineTo(screen.x, screen.y);
      }
    });
    path.closePath();
    if (style.fill) {
      ctx.fillStyle = style.fill;
      ctx.fill(path);
    }
    if (style.stroke) {
      ctx.strokeStyle = style.stroke;
      ctx.lineWidth = style.strokeWidth || 1.0;
      ctx.stroke(path);
    } else if (style.fill) {
      // Add subtle border to filled cells to reduce "floating" effect
      ctx.strokeStyle = "rgba(0,0,0,0.15)";
      ctx.lineWidth = 0.5;
      ctx.stroke(path);
    }
  });
};

Bin2CellExplorer.drawLineList = function (ctx, lines, style, projectFn) {
  if (!lines || !lines.length) return;
  ctx.strokeStyle = style.stroke || "rgba(150,150,150,0.6)";
  ctx.lineWidth = style.strokeWidth || 1.0;
  lines.forEach(function (line) {
    if (!line || line.length < 2) return;
    const path = new Path2D();
    line.forEach(function (pt, idx) {
      const screen = projectFn(pt);
      if (idx === 0) {
        path.moveTo(screen.x, screen.y);
      } else {
        path.lineTo(screen.x, screen.y);
      }
    });
    ctx.stroke(path);
  });
};

Bin2CellExplorer.inflateGeometry = function (data) {
  if (!data || !data.tile) return data;
  const tileId = data.tile.id;
  const sig = Bin2CellExplorer.geometrySignature(data);
  if (data.geometry) {
    Bin2CellExplorer.state.geometryCache[tileId] = {
      signature: sig,
      geometry: data.geometry
    };
  }
  const cached = Bin2CellExplorer.state.geometryCache[tileId];
  if (!cached || !cached.geometry) {
    return data;
  }

  const geom = cached.geometry;
  const result = JSON.parse(JSON.stringify(data));
  const useExpanded = !!result.all_expanded_outline;
  const polySource = useExpanded ? geom.polygons_exp : geom.polygons_raw;

  if (!result.expanded_outline || !result.expanded_outline.length) {
    result.expanded_outline = geom.outline_exp || [];
  }
  if (!result.nuclei_outline || !result.nuclei_outline.length) {
    result.nuclei_outline = geom.outline_raw || [];
  }

  // Per-label outlines cache (for selected-only mode)
  const perLabelNuclei = geom.per_label_nuclei || {};
  const perLabelExpanded = geom.per_label_expanded || {};

  if (result.overlay_type === "gene") {
    (result.overlays || []).forEach(function (overlay) {
      (overlay.features || []).forEach(function (feature) {
        // Restore polygons
        if (!feature.polygons || !feature.polygons.length) {
          const polys = (polySource && polySource[feature.label]) || [];
          feature.polygons = polys;
        }
        // Restore per-feature nuclei outlines
        if (!feature.nuclei_outline_paths || !feature.nuclei_outline_paths.length) {
          feature.nuclei_outline_paths = perLabelNuclei[feature.label] || [];
        }
        // Restore per-feature expanded outlines
        if (!feature.expanded_outline_paths || !feature.expanded_outline_paths.length) {
          feature.expanded_outline_paths = perLabelExpanded[feature.label] || [];
        }
      });
    });
  } else if (result.overlay_type === "observation") {
    (result.features || []).forEach(function (feature) {
      // Restore polygons
      if (!feature.polygons || !feature.polygons.length) {
        const polys = (polySource && polySource[feature.label]) || [];
        feature.polygons = polys;
      }
      // Restore per-feature nuclei outlines
      if (!feature.nuclei_outline_paths || !feature.nuclei_outline_paths.length) {
        feature.nuclei_outline_paths = perLabelNuclei[feature.label] || [];
      }
      // Restore per-feature expanded outlines
      if (!feature.expanded_outline_paths || !feature.expanded_outline_paths.length) {
        feature.expanded_outline_paths = perLabelExpanded[feature.label] || [];
      }
    });
  }
  return result;
};

Bin2CellExplorer.drawOverlayToCanvas = function (payload) {
  const canvas = Bin2CellExplorer.ensureCanvasOverlay();
  const ctx = Bin2CellExplorer.state.canvasCtx;
  const viewer = tmapp[tmapp["object_prefix"] + "_viewer"];
  if (!canvas || !ctx || !viewer) return;
  if (!viewer.world || !viewer.world.getItemCount()) return;
  const image = viewer.world.getItemAt(0);
  if (!image) return;

  const width = viewer.container ? viewer.container.clientWidth : canvas.clientWidth;
  const height = viewer.container ? viewer.container.clientHeight : canvas.clientHeight;
  const dpr = window.devicePixelRatio || 1;
  canvas.width = width * dpr;
  canvas.height = height * dpr;
  canvas.style.width = width + "px";
  canvas.style.height = height + "px";
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, width, height);

  const projectFn = function (pt) {
    return Bin2CellExplorer.projectPoint(image, viewer, pt);
  };

  // Draw outlines based on payload flags
  const showExpandedAll = payload.all_expanded_outline || false;
  const showExpandedSelected = payload.expanded_outlines_selected || false;
  const showNucleiAll = payload.all_nuclei_outline || false;
  const showNucleiSelected = payload.nuclei_outlines_selected || false;

  // Collect label IDs of selected/displayed cells for filtering
  let selectedLabels = null;
  if (showExpandedSelected || showNucleiSelected) {
    selectedLabels = new Set();
    if (payload.overlay_type === "gene") {
      (payload.overlays || []).forEach(function (overlay) {
        if (Bin2CellExplorer.state.hiddenGenes[overlay.gene]) return;
        (overlay.features || []).forEach(function (feature) {
          if (feature.label != null) {
            selectedLabels.add(feature.label);
          }
        });
      });
    } else if (payload.overlay_type === "observation") {
      (payload.features || []).forEach(function (feature) {
        if (Bin2CellExplorer.state.hiddenCategories[feature.category]) return;
        if (feature.label != null) {
          selectedLabels.add(feature.label);
        }
      });
    }
  }

  // Draw expanded outlines (all cells)
  if (showExpandedAll && payload.expanded_outline && payload.expanded_outline.length) {
    Bin2CellExplorer.drawLineList(ctx, payload.expanded_outline, {
      stroke: "rgba(100,100,100,0.9)",  // Darker gray, more opaque - matches selected-only
      strokeWidth: 1.0
    }, projectFn);
  }

  // Draw nuclei outlines (all cells)
  if (showNucleiAll && payload.nuclei_outline && payload.nuclei_outline.length) {
    const nucleiColor = payload.nuclei_outline_color || "#000000";
    const nucleiAlpha = payload.nuclei_outline_alpha != null ? payload.nuclei_outline_alpha : 0.6;
    const nucleiStroke = Bin2CellExplorer.hexToRgba(nucleiColor, nucleiAlpha);
    Bin2CellExplorer.drawLineList(ctx, payload.nuclei_outline, {
      stroke: nucleiStroke,
      strokeWidth: 0.8
    }, projectFn);
  }

  // Get overlay_alpha from payload
  const overlayAlpha = payload.overlay_alpha != null ? payload.overlay_alpha : 0.7;
  const renderMode = payload.render_mode || "fill";
  const highlightColor = payload.highlight_color || "#39ff14";
  const highlightWidth = payload.highlight_width || 2.0;

  if (payload.overlay_type === "gene") {
    (payload.overlays || []).forEach(function (overlay) {
      if (Bin2CellExplorer.state.hiddenGenes[overlay.gene]) return;
      (overlay.features || []).forEach(function (feature) {
        if (!feature.polygons || !feature.polygons.length) return;

        let fillColor = null;
        let strokeColor = null;
        let strokeWidth = 1.0;

        if (renderMode === "fill") {
          // Fill mode: show colored fills with subtle borders
          if (feature.fill) {
            fillColor = Bin2CellExplorer.applyAlpha(feature.fill, overlayAlpha);
          }
        } else if (renderMode === "outline") {
          // Outline mode: no fill, use highlight color for stroke
          fillColor = null;
          strokeColor = Bin2CellExplorer.applyAlpha(highlightColor, overlayAlpha);
          strokeWidth = highlightWidth;
        }

        const style = {
          fill: fillColor,
          stroke: strokeColor,
          strokeWidth: strokeWidth
        };
        Bin2CellExplorer.drawPolygonList(ctx, feature.polygons, style, projectFn);

        // Draw per-feature outlines if selected-only mode is enabled
        if (showExpandedSelected && feature.expanded_outline_paths) {
          Bin2CellExplorer.drawLineList(ctx, feature.expanded_outline_paths, {
            stroke: "rgba(100,100,100,0.9)",  // Darker gray, more opaque
            strokeWidth: 1.0
          }, projectFn);
        }

        if (showNucleiSelected && feature.nuclei_outline_paths) {
          const nucleiColor = payload.nuclei_outline_color || "#000000";
          const nucleiAlpha = payload.nuclei_outline_alpha != null ? payload.nuclei_outline_alpha : 0.6;
          const nucleiStroke = Bin2CellExplorer.hexToRgba(nucleiColor, nucleiAlpha);
          Bin2CellExplorer.drawLineList(ctx, feature.nuclei_outline_paths, {
            stroke: nucleiStroke,
            strokeWidth: 0.8
          }, projectFn);
        }
      });
    });
  } else if (payload.overlay_type === "observation") {
    (payload.features || []).forEach(function (feature) {
      if (Bin2CellExplorer.state.hiddenCategories[feature.category]) return;
      if (!feature.polygons || !feature.polygons.length) return;

      let fillColor = null;
      let strokeColor = null;
      let strokeWidth = 1.0;

      if (renderMode === "fill") {
        // Fill mode: show colored fills
        if (feature.fill) {
          fillColor = Bin2CellExplorer.applyAlpha(feature.fill, overlayAlpha);
        }
      } else if (renderMode === "outline") {
        // Outline mode: no fill, use highlight color for stroke
        fillColor = null;
        strokeColor = Bin2CellExplorer.applyAlpha(highlightColor, overlayAlpha);
        strokeWidth = highlightWidth;
      }

      const style = {
        fill: fillColor,
        stroke: strokeColor,
        strokeWidth: strokeWidth
      };
      Bin2CellExplorer.drawPolygonList(ctx, feature.polygons, style, projectFn);

      // Draw per-feature outlines if selected-only mode is enabled
      if (showExpandedSelected && feature.expanded_outline_paths) {
        Bin2CellExplorer.drawLineList(ctx, feature.expanded_outline_paths, {
          stroke: "rgba(100,100,100,0.9)",  // Darker gray, more opaque
          strokeWidth: 1.0
        }, projectFn);
      }

      if (showNucleiSelected && feature.nuclei_outline_paths) {
        const nucleiColor = payload.nuclei_outline_color || "#000000";
        const nucleiAlpha = payload.nuclei_outline_alpha != null ? payload.nuclei_outline_alpha : 0.6;
        const nucleiStroke = Bin2CellExplorer.hexToRgba(nucleiColor, nucleiAlpha);
        Bin2CellExplorer.drawLineList(ctx, feature.nuclei_outline_paths, {
          stroke: nucleiStroke,
          strokeWidth: 0.8
        }, projectFn);
      }
    });
  }
};

Bin2CellExplorer.renderOverlay = function (data) {
  const resolved = Bin2CellExplorer.inflateGeometry(data);
  Bin2CellExplorer.state.overlays = resolved;
  Bin2CellExplorer.state.renderPayload = resolved;
  Bin2CellExplorer.state.lastOverlayPayload = resolved;
  Bin2CellExplorer.drawOverlayToCanvas(resolved);
  Bin2CellExplorer.renderLegend(resolved);
  Bin2CellExplorer.renderTileOverview();
  Bin2CellExplorer.ensureViewerHooks();
};


Bin2CellExplorer.imageToViewport = function (x, y, tiledImage) {
  const viewer = tmapp[tmapp["object_prefix"] + "_viewer"];
  if (!viewer) return { x: x, y: y };
  const image = tiledImage || viewer.world.getItemAt(0);
  if (!image) return { x: x, y: y };
  let px = x;
  if (image.getFlip && image.getFlip()) {
    px = image.getContentSize().x - px;
  }
  const point = image.imageToViewportCoordinates(px, y, true);
  return point;
};

Bin2CellExplorer.polygonsToPath = function (polygons) {
  if (!polygons || !polygons.length) return "";
  const viewer = tmapp[tmapp["object_prefix"] + "_viewer"];
  if (!viewer || !viewer.world.getItemCount()) return "";
  const image = viewer.world.getItemAt(0);
  const segments = [];
  polygons.forEach(function (poly) {
    if (!poly || poly.length < 3) return;
    let path = "";
    poly.forEach(function (pt, idx) {
      const vp = Bin2CellExplorer.imageToViewport(pt[0], pt[1], image);
      path += (idx === 0 ? "M" : "L") + vp.x + " " + vp.y;
    });
    path += "Z";
    segments.push(path);
  });
  return segments.join(" ");
};

Bin2CellExplorer.lineToPath = function (points) {
  if (!points || !points.length) return "";
  const viewer = tmapp[tmapp["object_prefix"] + "_viewer"];
  if (!viewer || !viewer.world.getItemCount()) return "";
  const image = viewer.world.getItemAt(0);
  let path = "";
  points.forEach(function (pt, idx) {
    const vp = Bin2CellExplorer.imageToViewport(pt[0], pt[1], image);
    path += (idx === 0 ? "M" : "L") + vp.x + " " + vp.y;
  });
  return path;
};

Bin2CellExplorer.renderLegend = function (data) {
  const legend = document.getElementById("Bin2CellExplorer_legend");
  if (!legend) return;
  legend.innerHTML = "";

  if (data.overlay_type === "gene") {
    Bin2CellExplorer.buildGeneLegend(legend, data.overlays || []);
  } else if (data.overlay_type === "observation") {
    Bin2CellExplorer.buildObservationLegend(legend, data.legend || {});
  }
};

Bin2CellExplorer.buildGeneLegend = function (container, overlays) {
  if (!overlays.length) {
    container.textContent = "No genes in overlay.";
    return;
  }
  const title = document.createElement("div");
  title.textContent = "Genes";
  title.className = "fw-bold mb-1";
  container.appendChild(title);

  overlays.forEach(function (overlay) {
    const gene = overlay.gene;
    const legend = overlay.legend || {};

    const row = document.createElement("div");
    row.className = "bin2cell-legend-row";

    const checkbox = document.createElement("input");
    checkbox.type = "checkbox";
    checkbox.checked = true;
    checkbox.dataset.gene = gene;
    checkbox.addEventListener("change", function (evt) {
      Bin2CellExplorer.toggleGeneLayer(gene, evt.target.checked);
    });

    const label = document.createElement("label");
    label.textContent = " " + gene;
    label.className = "ms-1";

    row.appendChild(checkbox);
    row.appendChild(label);

    if (legend.type === "continuous" && legend.gradient) {
      const gradient = document.createElement("div");
      gradient.className = "bin2cell-gradient";
      gradient.style.height = "12px";
      gradient.style.flex = "1";
      gradient.style.marginLeft = "8px";
      gradient.style.borderRadius = "4px";
      gradient.style.background = "linear-gradient(to right," + legend.gradient.map(function (stop) {
        return stop[1];
      }).join(",") + ")";
      gradient.title = "min: " + legend.min + " max: " + legend.max;
      row.appendChild(gradient);
    } else if (legend.type === "solid" && legend.color) {
      const swatch = document.createElement("span");
      swatch.className = "bin2cell-swatch";
      swatch.style.display = "inline-block";
      swatch.style.width = "16px";
      swatch.style.height = "16px";
      swatch.style.marginLeft = "8px";
      swatch.style.borderRadius = "3px";
      swatch.style.border = "1px solid rgba(0,0,0,0.2)";
      swatch.style.background = legend.color;
      row.appendChild(swatch);
    }
    container.appendChild(row);
  });
};

Bin2CellExplorer.buildObservationLegend = function (container, legendData) {
  const items = legendData.items || [];
  if (!items.length) {
    container.textContent = "No observation categories.";
    return;
  }
  const title = document.createElement("div");
  title.textContent = legendData.obs_col || "Observation";
  title.className = "fw-bold mb-1";
  container.appendChild(title);

  items.forEach(function (item) {
    const row = document.createElement("div");
    row.className = "bin2cell-legend-row";

    const checkbox = document.createElement("input");
    checkbox.type = "checkbox";
    checkbox.checked = true;
    checkbox.dataset.category = item.label;
    checkbox.addEventListener("change", function (evt) {
      Bin2CellExplorer.toggleCategoryLayer(item.label, evt.target.checked);
    });

    const swatch = document.createElement("span");
    swatch.className = "bin2cell-swatch";
    swatch.style.display = "inline-block";
    swatch.style.width = "14px";
    swatch.style.height = "14px";
    swatch.style.marginLeft = "8px";
    swatch.style.borderRadius = "3px";
    swatch.style.background = item.color;

    const label = document.createElement("label");
    label.textContent = " " + item.label;
    label.className = "ms-1";

    row.appendChild(checkbox);
    row.appendChild(swatch);
    row.appendChild(label);
    container.appendChild(row);
  });
};

Bin2CellExplorer.toggleGeneLayer = function (gene, visible) {
  if (!gene) return;
  if (visible) {
    delete Bin2CellExplorer.state.hiddenGenes[gene];
  } else {
    Bin2CellExplorer.state.hiddenGenes[gene] = true;
  }
  Bin2CellExplorer.scheduleOverlayRefresh();
};

Bin2CellExplorer.toggleCategoryLayer = function (category, visible) {
  if (!category) return;
  if (visible) {
    delete Bin2CellExplorer.state.hiddenCategories[category];
  } else {
    Bin2CellExplorer.state.hiddenCategories[category] = true;
  }
  Bin2CellExplorer.scheduleOverlayRefresh();
};

Bin2CellExplorer.ensureViewerHooks = function () {
  if (Bin2CellExplorer.state.viewerHooksInstalled) return;
  const viewer = tmapp[tmapp["object_prefix"] + "_viewer"];
  if (!viewer) return;
  const schedule = function () {
    Bin2CellExplorer.scheduleOverlayRefresh();
  };
  viewer.addHandler("animation", schedule);
  viewer.addHandler("resize", schedule);
  viewer.addHandler("flip", schedule);
  viewer.addHandler("rotate", schedule);
  Bin2CellExplorer.state.viewerHooksInstalled = true;
};

Bin2CellExplorer.scheduleOverlayRefresh = function () {
  if (!Bin2CellExplorer.state.renderPayload) return;
  if (Bin2CellExplorer.state.overlayRefreshToken) {
    cancelAnimationFrame(Bin2CellExplorer.state.overlayRefreshToken);
  }
  Bin2CellExplorer.state.overlayRefreshToken = requestAnimationFrame(function () {
    Bin2CellExplorer.state.overlayRefreshToken = null;
    Bin2CellExplorer.rebuildOverlayGeometry();
  });
};

Bin2CellExplorer.rebuildOverlayGeometry = function () {
  const payload = Bin2CellExplorer.state.renderPayload || Bin2CellExplorer.state.lastOverlayPayload;
  if (!payload) return;
  Bin2CellExplorer.drawOverlayToCanvas(payload);
};

// Helper to apply alpha to any color (hex or rgba)
Bin2CellExplorer.applyAlpha = function (color, alpha) {
  if (!color) return null;

  // If it's already rgba, replace the alpha
  if (color.startsWith('rgba(')) {
    return color.replace(/,\s*[\d.]+\)$/, `, ${alpha})`);
  }

  // If it's rgb, convert to rgba
  if (color.startsWith('rgb(')) {
    return color.replace(')', `, ${alpha})`).replace('rgb', 'rgba');
  }

  // If it's hex, convert to rgba
  if (color.startsWith('#')) {
    return Bin2CellExplorer.hexToRgba(color, alpha);
  }

  return color;
};

// Helper to convert #RRGGBB to rgba(r,g,b,a)
Bin2CellExplorer.hexToRgba = function (hex, alpha) {
  try {
    const h = String(hex || "").trim();
    let r = 160, g = 160, b = 160;
    if (/^#?[0-9a-fA-F]{6}$/.test(h)) {
      const s = h.charAt(0) === '#' ? h.substring(1) : h;
      r = parseInt(s.substring(0, 2), 16);
      g = parseInt(s.substring(2, 4), 16);
      b = parseInt(s.substring(4, 6), 16);
    }
    const a = Math.max(0, Math.min(1, Number(alpha)));
    return `rgba(${r},${g},${b},${a})`;
  } catch (e) {
    return `rgba(160,160,160,${Math.max(0, Math.min(1, Number(alpha) || 0.5))})`;
  }
};
// Ensure filename-based global exists for TissUUmaps plugin loader, cause we changed name
var CellExplorer = Bin2CellExplorer;
