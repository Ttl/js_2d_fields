import { Complex } from './complex.js';
import { usableSweepPoints } from './sparameters.js';
import { exportSnP, downloadFile } from './snp_export.js';
import { darkAxis, darkBackground, draw, fieldsHaveWantedView, drawResultsPlot, drawSParamPlot, drawParameterSweepPlot, setGlobals, getScaleRange, setScaleRange, getActualDataRange,
    freeze, unfreeze, isFrozen, conductorFillShapes, dielectricFillShapes, computeGeometryView, displayTop,
    centroidHoverTrace, triMeanE, updateTriImage } from './plot.js';
import { initCustomGeometryEditor, activateCustomGeometry, validateCustomGeometry, getCustomGeometryText, setCustomGeometryText,
         getCustomOverrides, customSweepParams, flushCustomGeometryEdits } from './custom_geometry_editor.js';
import { solverToGeometryText, setParamInText } from './custom_geometry_text.js';
import { initLayoutPanels, syncLogPanel, setLogStatus, logSolveStarted } from './layout_panels.js';
import { buildSolverFromParams as _buildSolverFromParams, platingOptions } from './solver_factory.js';
import { DEFAULT_GRID_N } from './mesher.js';

// Lazy Plotly access - allows app to function while Plotly is loading
const getPlotly = () => window.Plotly;

// The main thread's solver is a GEOMETRY/PLOT view model: it is built here so the
// geometry preview and the field plots have something to read, but it never solves.
// Every solve runs in the worker (see solve_worker.js) and its field arrays are grafted
// back on. Construction is cheap (no meshing or linear algebra).
const buildSolverFromParams = (p) => _buildSolverFromParams(p, log);

let solver = null;
let isSimulating = false;
let frequencySweepResults = null;  // Array of {freq, result} objects
let currentTab = 'geometry';
let geometryChanged = false;  // Track if geometry has changed since last solve
let lastSolvedGeometry = null;  // Hash of geometry params from last solve
let lastSolvedFrequency = null;  // Frequency params from last solve
let isSweeping = false;
let parameterSweepResults = null;
let lastSweepGeometry = null;  // Geometry hash at sweep time (excluding swept param)
let lastSweepParam = null;     // Which parameter was swept
let lastSweepDisplayUnit = null; // Display unit used during last sweep

// Modes tab state
let isSolvingModes = false;
// Field plot state: the view-model solver whose fields the worker's retained solve
// belongs to, the simulate job of that solve, and its highest frequency (the default
// plot frequency).
let plotFieldsSolver = null, plotFieldsJob = null;
let plotFieldsMaxFreq = null;
let plotFieldsSolveKey = null;   // solveInputKey() of that solve
let plotFieldsBusy = false, plotFieldsPending = false;
let modesResult = null;          // last solveModes() summary { modes, nconv, ... }
let modesSelectedIdx = -1;       // index into modesResult.modes of the plotted mode
let modesFieldCache = new Map(); // mode index, field grid fetched from the worker
let modesSolver = null;          // INDEPENDENT solver for the Modes tab — built from the
                                 // sidebar geometry but kept separate from the main `solver`,
                                 // so solving modes never disturbs the main solve's results.
let lastModesGeometry = null;    // geometry hash at the last modes solve (staleness notice)
let lastModesFrequency = null;   // modes frequency at the last modes solve (staleness notice)

// Full parameter config table. Each entry drives input writing + axis labeling.
// fixedUnit: cosmetic axis label for plain-number inputs (sigma). Absent = derive from geometry input.
const SWEEP_PARAM_CONFIG = {
    custom_sigma:       { label: 'Conductivity',          inputId: 'inp_custom_sigma', fixedUnit: 'S/m', group: 'custom' },
    // Always available
    w:                  { label: 'Trace Width',          inputId: 'inp_w',             group: 'always' },
    h:                  { label: 'Substrate Height',     inputId: 'inp_h',             group: 'always' },
    t:                  { label: 'Trace Thickness',      inputId: 'inp_t',             group: 'always' },
    er:                 { label: 'Permittivity',         inputId: 'inp_er',            group: 'always' },
    tand:               { label: 'Loss Tangent',         inputId: 'inp_tand',          group: 'always' },
    sigma:              { label: 'Conductivity',         inputId: 'inp_sigma',         fixedUnit: 'S/m', group: 'always' },
    // 'shared' = applies to EVERY line type, including those (coax) whose own geometry
    // block replaces the 'always' set above.
    rq:                 { label: 'Surface Roughness',    inputId: 'inp_rq',            group: 'shared' },
    // Coaxial only. No "(coax)" suffix: unlike the stripline entries below, these are
    // never listed alongside the 'always' group (updateSweepParamList turns that group
    // off for coax), so the same plain names as microstrip are unambiguous here.
    coax_d:             { label: 'Inner Conductor Diameter', inputId: 'inp_coax_d',    group: 'coax' },
    coax_D:             { label: 'Dielectric Diameter',      inputId: 'inp_coax_D',    group: 'coax' },
    coax_er:            { label: 'Permittivity',             inputId: 'inp_coax_er',   group: 'coax' },
    coax_tand:          { label: 'Loss Tangent',             inputId: 'inp_coax_tand', group: 'coax' },
    coax_sigma:         { label: 'Conductivity',             inputId: 'inp_coax_sigma', fixedUnit: 'S/m', group: 'coax' },
    // Rectangular waveguide only. Same reasoning as the coax block above.
    wg_a:               { label: 'Broad Wall (a)',           inputId: 'inp_wg_a',      group: 'waveguide' },
    wg_b:               { label: 'Narrow Wall (b)',          inputId: 'inp_wg_b',      group: 'waveguide' },
    wg_er:              { label: 'Permittivity',             inputId: 'inp_wg_er',     group: 'waveguide' },
    wg_tand:            { label: 'Loss Tangent',             inputId: 'inp_wg_tand',   group: 'waveguide' },
    wg_sigma:           { label: 'Conductivity',             inputId: 'inp_wg_sigma', fixedUnit: 'S/m', group: 'waveguide' },
    // Differential types only
    trace_spacing:      { label: 'Trace Spacing',        inputId: 'inp_trace_spacing', group: 'diff' },
    // GCPW types only
    gap:                { label: 'GCPW Gap',              inputId: 'inp_gap',           group: 'gcpw' },
    via_gap:            { label: 'Via Gap',               inputId: 'inp_via_gap',       group: 'gcpw' },
    gnd_width:          { label: 'Ground Width',          inputId: 'inp_gnd_width',     group: 'gcpw' },
    // Stripline types only
    stripline_top_h:    { label: 'Top Dielectric Height (stripline)', inputId: 'inp_air_top',    group: 'stripline' },
    er_top:             { label: 'Top Permittivity (stripline)',       inputId: 'inp_er_top',     group: 'stripline' },
    tand_top:           { label: 'Top Loss Tangent (stripline)',       inputId: 'inp_tand_top',   group: 'stripline' },
    // Solder mask (if enabled)
    sm_t_sub:           { label: 'Solder Mask Thickness (substrate side)', inputId: 'inp_sm_t_sub',  group: 'sm' },
    sm_t_trace:         { label: 'Solder Mask Thickness (trace top)',       inputId: 'inp_sm_t_trace', group: 'sm' },
    sm_t_side:          { label: 'Solder Mask Thickness (trace side)',      inputId: 'inp_sm_t_side', group: 'sm' },
    sm_er:              { label: 'Solder Mask Permittivity',            inputId: 'inp_sm_er',     group: 'sm' },
    sm_tand:            { label: 'Solder Mask Loss Tangent',            inputId: 'inp_sm_tand',   group: 'sm' },
    // Top dielectric (if enabled)
    top_diel_h:         { label: 'Top Dielectric Height',   inputId: 'inp_top_diel_h',  group: 'top_diel' },
    top_diel_er:        { label: 'Top Dielectric Permittivity', inputId: 'inp_top_diel_er', group: 'top_diel' },
    top_diel_tand:      { label: 'Top Dielectric Loss Tangent', inputId: 'inp_top_diel_tand', group: 'top_diel' },
    // Ground cutout (if enabled)
    gnd_cut_w:          { label: 'Ground Cutout Width',  inputId: 'inp_gnd_cut_w',    group: 'gnd_cut' },
    gnd_cut_h:          { label: 'Ground Cutout Height', inputId: 'inp_gnd_cut_h',    group: 'gnd_cut' },
    // Enclosure (if enabled)
    enclosure_width:    { label: 'Enclosure Width',  inputId: 'inp_enclosure_width',  group: 'enclosure' },
    enclosure_height:   { label: 'Enclosure Height', inputId: 'inp_enclosure_height', group: 'enclosure' },
    // Plating (if enabled)
    plating_t:          { label: 'Plating Thickness',    inputId: 'inp_plating_t',   group: 'plating' },
    plating_sigma:      { label: 'Plating Conductivity', inputId: 'inp_plating_sigma', fixedUnit: 'S/m', group: 'plating' },
    plating_rq:         { label: 'Plating Roughness',    inputId: 'inp_plating_rq',  group: 'plating' },
    plating_rq_iface:   { label: 'Plating Interface Roughness', inputId: 'inp_plating_rq_iface', group: 'plating' },
};

// --- Unit Parsing Helper ---

/**
 * Get value from input field with unit parsing
 * Returns value in SI base units (meters for length, Hz for frequency)
 * @param {string} id - Input element ID
 * @returns {number} - Parsed value in SI units
 */
function getInputValue(id) {
    const element = document.getElementById(id);
    if (!element) return NaN;

    const defaultUnit = window.getDefaultUnit ? window.getDefaultUnit(id) : '';

    // Use value if present, otherwise fallback to placeholder
    let raw = element.value;
    if (!raw || raw.trim() === '') {
        raw = element.placeholder || '';
    }
    // Placeholder keywords pass through: "auto" (enclosure size), "full" (ground width)
    if (raw === "auto" || raw === "full") {
        return raw;
    }

    return window.parseValueWithUnit
        ? window.parseValueWithUnit(raw, defaultUnit)
        : parseFloat(raw);
}

function getInputValueUnitless(id) {
    const el = document.getElementById(id);
    if (!el) return NaN;

    let raw = el.value;
    if (!raw || raw.trim() === '') {
        raw = el.placeholder || '';
    }

    return parseFloat(raw);
}

// --- URL Parameter Serialization ---

/**
 * Default settings (in display units, matching what getUISettings returns)
 * These are used to filter out default values from URL parameters.
 * Doesn't need to match HTML defaults.
 * DO NOT CHANGE OR ALL EXISTING LINKS WILL BREAK.
 */
// Sweep config of a key: a fixed entry, or a geometry parameter of the custom type.
function sweepParamConfig(key) {
    return SWEEP_PARAM_CONFIG[key] || customSweepParams().find(c => c.key === key) || null;
}

const DEFAULT_SETTINGS = {
    tl_type: 'microstrip',
    mesh_backend: 'rectilinear',
    w: 0.35,           // mm
    h: 0.21,           // mm
    t: 35,             // μm
    er: 4.4,
    tand: 0.02,
    sigma: 5.8e7,
    freq_start: 0.1,   // GHz
    freq_stop: 10,     // GHz
    freq_points: 10,
    trace_spacing: 0.2, // mm
    gap: 0.1,          // mm
    via_gap: 0.1,      // mm
    gnd_width: NaN,    // full domain width
    stripline_top_h: 0.4, // mm
    er_top: 4.5,
    tand_top: 0.02,
    use_sm: 0,
    sm_t_sub: 20,      // μm
    sm_t_trace: 20,    // μm
    sm_t_side: 20,     // μm
    sm_er: 3.5,
    sm_tand: 0.02,
    use_top_diel: 0,
    top_diel_h: 0.2,   // mm
    top_diel_er: 4.5,
    top_diel_tand: 0.02,
    use_gnd_cut: 0,
    gnd_cut_w: 0.5,    // mm
    gnd_cut_h: 0.5,    // mm
    use_enclosure: 0,
    use_side_gnd: 0,
    use_top_gnd: 0,
    enclosure_width: NaN,  // auto
    enclosure_height: NaN, // auto
    max_iters: 10,
    tolerance: 0.01,
    min_converged_passes: 2,
    estimate_error: 1,
    max_nodes: 20,
    rq: 0,             // μm
    use_plating: 0,
    plating_sigma: 1e7,
    plating_t: 4,    // μm
    plating_rq: 0,   // μm
    plating_rq_iface: NaN,  // μm, same as plating_rq
    plating_top: 1,
    plating_sides: 1,
    plating_bottom: 0,
    plating_thick_corners: 1,
    sparam_length: 10, // mm
    sparam_z_ref: 50,
    use_causal_materials: 0,
    interp_sweep: 1,
    interp_tolerance: 0.5,
    // Modes tab (eigenmode viewer)
    modes_freq: 10,    // GHz
    modes_nev: 6,
    modes_mesh_density: 8,    // bulk cells per wavelength (TriBackend wavelengthDensity)
    modes_shrink_domain: true, // mesh only the field region of an auto-sized open domain
    // Broadside coupled stripline (display units: mm, μm)
    bs_w: 0.2,           // mm
    bs_t: 35,            // μm
    bs_x_offset: 0,      // mm
    bs_sigma: 5.8e7,
    bs_h_bottom: 0.2,    // mm
    bs_er_bottom: 4.4,
    bs_tand_bottom: 0.02,
    bs_h_middle: 0.2,    // mm
    bs_er_middle: 4.4,
    bs_tand_middle: 0.02,
    bs_h_top: 0.2,       // mm
    bs_er_top: 4.4,
    bs_tand_top: 0.02,

    // Coaxial (display units: mm). Defaults are RG-402-like semi-rigid, ~48 ohm.
    coax_d: 0.92,        // mm, inner conductor diameter
    coax_D: 2.95,        // mm, dielectric diameter = shield inner diameter
    coax_er: 2.1,
    coax_tand: 0.0002,
    coax_sigma: 5.8e7,
    // Which conductor the plating lands on, coax's counterpart to plating_top/sides/
    // bottom. Named coax_* so TYPE_ONLY_KEYS keeps them out of every other type's link.
    coax_plating_inner: 1,
    coax_plating_outer: 0,

    // Rectangular waveguide (display units: mm). Defaults are WR-90 (X band):
    // fc = 6.557 GHz, second cutoff 13.114 GHz, published band 8.2-12.4 GHz.
    wg_a: 22.86,         // mm, broad inner wall
    wg_b: 10.16,         // mm, narrow inner wall
    wg_er: 1.0,
    wg_tand: 0,
    wg_sigma: 5.8e7,

    // Custom geometry: the geometry text and the conductor conductivity.
    custom_geom: '',
    custom_sigma: 5.8e7,
};

// Certified error (fraction) as a percentage with three significant digits,
// scientific notation below 0.001%.
function fmtErrPct(err) {
    const pct = 100 * err;
    if (pct === 0) return '0';
    return pct < 1e-3 ? pct.toExponential(2) : pct.toPrecision(3);
}

// The sidebar inputs behind the settings (getUISettings, share links, restoreSettings)
// and the solve parameters (getParams), in settings key order.
// kind: unit (a number in the input's unit: display units in the settings, SI in the
//   params), num, int, float (parseFloat of the value, no placeholder), chk (1/0 in the
//   settings, boolean in the params), bool (boolean in the settings), pct (percent in the
//   input, a fraction in the settings and the params).
// use: s = settings, p = params, h = geometry hash, l = edits redraw the geometry.
const SETTINGS_FIELDS = [
    ['custom_sigma', 'inp_custom_sigma', 'num', 'sl'],
    ['w', 'inp_w', 'unit', 'sphl'],
    ['h', 'inp_h', 'unit', 'sphl'],
    ['t', 'inp_t', 'unit', 'sphl'],
    ['er', 'inp_er', 'num', 'sphl'],
    ['tand', 'inp_tand', 'num', 'sphl'],
    ['sigma', 'inp_sigma', 'num', 'sphl'],
    ['freq_start', 'freq-start', 'unit', 'spl'],
    ['freq_stop', 'freq-stop', 'unit', 's'],
    ['freq_points', 'freq-points', 'int', 's'],
    ['trace_spacing', 'inp_trace_spacing', 'unit', 'sphl'],
    ['gap', 'inp_gap', 'unit', 'sphl'],
    ['via_gap', 'inp_via_gap', 'unit', 'sphl'],
    ['gnd_width', 'inp_gnd_width', 'unit', 'sphl'],
    ['stripline_top_h', 'inp_air_top', 'unit', 'sphl'],
    ['er_top', 'inp_er_top', 'num', 'sphl'],
    ['tand_top', 'inp_tand_top', 'num', 'sphl'],
    ['use_sm', 'chk_solder_mask', 'chk', 'sphl'],
    ['sm_t_sub', 'inp_sm_t_sub', 'unit', 'sphl'],
    ['sm_t_trace', 'inp_sm_t_trace', 'unit', 'sphl'],
    ['sm_t_side', 'inp_sm_t_side', 'unit', 'sphl'],
    ['sm_er', 'inp_sm_er', 'num', 'sphl'],
    ['sm_tand', 'inp_sm_tand', 'num', 'sphl'],
    ['use_top_diel', 'chk_top_diel', 'chk', 'sphl'],
    ['top_diel_h', 'inp_top_diel_h', 'unit', 'sphl'],
    ['top_diel_er', 'inp_top_diel_er', 'num', 'sphl'],
    ['top_diel_tand', 'inp_top_diel_tand', 'num', 'sphl'],
    ['use_gnd_cut', 'chk_gnd_cut', 'chk', 'sphl'],
    ['gnd_cut_w', 'inp_gnd_cut_w', 'unit', 'sphl'],
    ['gnd_cut_h', 'inp_gnd_cut_h', 'unit', 'sphl'],
    ['use_enclosure', 'chk_enclosure', 'chk', 'sphl'],
    ['use_side_gnd', 'chk_side_gnd', 'chk', 'sphl'],
    ['use_top_gnd', 'chk_top_gnd', 'chk', 'sphl'],
    ['enclosure_width', 'inp_enclosure_width', 'unit', 'sphl'],
    ['enclosure_height', 'inp_enclosure_height', 'unit', 'sphl'],
    ['max_iters', 'inp_max_iters', 'int', 'sp'],
    ['tolerance', 'inp_tolerance', 'pct', 'sp'],
    ['min_converged_passes', 'inp_min_converged_passes', 'num', 'sp'],
    ['estimate_error', 'chk_estimate_error', 'chk', 'sp'],
    ['max_nodes', 'inp_max_nodes', 'int', 'sp'],
    ['rq', 'inp_rq', 'unit', 'sphl'],
    ['use_plating', 'chk_plating', 'chk', 'sphl'],
    ['plating_sigma', 'inp_plating_sigma', 'num', 'sphl'],
    ['plating_t', 'inp_plating_t', 'unit', 'sphl'],
    ['plating_rq', 'inp_plating_rq', 'unit', 'sphl'],
    ['plating_rq_iface', 'inp_plating_rq_iface', 'unit', 'sphl'],
    ['plating_top', 'chk_plating_top', 'chk', 'sphl'],
    ['plating_sides', 'chk_plating_sides', 'chk', 'sphl'],
    ['plating_bottom', 'chk_plating_bottom', 'chk', 'sphl'],
    ['plating_thick_corners', 'chk_plating_thick_corners', 'chk', 'sphl'],
    ['sparam_length', 'sparam-length', 'unit', 's'],
    ['sparam_z_ref', 'sparam-z-ref', 'num', 's'],
    ['use_causal_materials', 'chk_causal_materials', 'chk', 'sph'],
    ['interp_sweep', 'chk_interp_sweep', 'chk', 's'],
    ['interp_tolerance', 'interp_tolerance', 'float', 's'],
    ['modes_freq', 'modes-freq', 'unit', 's'],
    ['modes_nev', 'modes-nev', 'int', 's'],
    ['modes_mesh_density', 'modes-mesh-density', 'int', 's'],
    ['modes_shrink_domain', 'modes-shrink-domain', 'bool', 's'],
    // Broadside coupled stripline
    ['bs_w', 'inp_bs_w', 'unit', 'sphl'],
    ['bs_t', 'inp_bs_t', 'unit', 'sphl'],
    ['bs_x_offset', 'inp_bs_x_offset', 'unit', 'sphl'],
    ['bs_sigma', 'inp_bs_sigma', 'num', 'sphl'],
    ['bs_h_bottom', 'inp_bs_h_bottom', 'unit', 'sphl'],
    ['bs_er_bottom', 'inp_bs_er_bottom', 'num', 'sphl'],
    ['bs_tand_bottom', 'inp_bs_tand_bottom', 'num', 'sphl'],
    ['bs_h_middle', 'inp_bs_h_middle', 'unit', 'sphl'],
    ['bs_er_middle', 'inp_bs_er_middle', 'num', 'sphl'],
    ['bs_tand_middle', 'inp_bs_tand_middle', 'num', 'sphl'],
    ['bs_h_top', 'inp_bs_h_top', 'unit', 'sphl'],
    ['bs_er_top', 'inp_bs_er_top', 'num', 'sphl'],
    ['bs_tand_top', 'inp_bs_tand_top', 'num', 'sphl'],
    // Coaxial (diameters in; CoaxSolver derives the radii)
    ['coax_d', 'inp_coax_d', 'unit', 'sphl'],
    ['coax_D', 'inp_coax_D', 'unit', 'sphl'],
    ['coax_er', 'inp_coax_er', 'num', 'sphl'],
    ['coax_tand', 'inp_coax_tand', 'num', 'sphl'],
    ['coax_sigma', 'inp_coax_sigma', 'num', 'sphl'],
    ['coax_plating_inner', 'chk_plating_inner', 'chk', 'sphl'],
    ['coax_plating_outer', 'chk_plating_outer', 'chk', 'sphl'],
    // Rectangular waveguide (inner wall dimensions)
    ['wg_a', 'inp_wg_a', 'unit', 'sphl'],
    ['wg_b', 'inp_wg_b', 'unit', 'sphl'],
    ['wg_er', 'inp_wg_er', 'num', 'sphl'],
    ['wg_tand', 'inp_wg_tand', 'num', 'sphl'],
    ['wg_sigma', 'inp_wg_sigma', 'num', 'sphl'],
];
const fieldsFor = (use) => SETTINGS_FIELDS.filter(f => f[3].includes(use));

// Value of an input as stored in the settings (display units).
function readSetting(id, kind) {
    const el = document.getElementById(id);
    switch (kind) {
        case 'unit': {
            if (!el) return NaN;
            const defaultUnit = window.getDefaultUnit ? window.getDefaultUnit(id) : '';
            const siValue = window.parseValueWithUnit ?
                window.parseValueWithUnit(el.value, defaultUnit) :
                parseFloat(el.value);
            // Convert back to display units for serialization
            const unitMap = { 'mm': 1e3, 'μm': 1e6, 'GHz': 1e-9, 'm': 1 };
            return siValue * (unitMap[defaultUnit] || 1);
        }
        case 'num': return getInputValueUnitless(id);
        case 'int': return parseInt(el.value);
        case 'float': return parseFloat(el.value);
        case 'chk': return el.checked ? 1 : 0;
        case 'bool': return el.checked;
        // The input is in percent. The solvers (and saved settings / share links)
        // use the fraction.
        case 'pct': return getInputValueUnitless(id) / 100;
    }
}

// Value of an input as passed to the solvers (SI units).
function readParam(id, kind) {
    switch (kind) {
        case 'unit': return getInputValue(id);
        case 'chk': return document.getElementById(id).checked;
        default: return readSetting(id, kind);
    }
}

/**
 * Get current UI settings as a serializable object (in display units)
 */
function getUISettings() {
    const settings = {
        tl_type: document.getElementById('tl_type').value,
        mesh_backend: (document.getElementById('mesh_backend')?.value) ?? 'rectilinear',
        custom_geom: getCustomGeometryText(),
    };
    for (const [key, id, kind] of fieldsFor('s')) settings[key] = readSetting(id, kind);
    return settings;
}

/**
 * Serialize settings to URL-safe base64 string
 * Only includes non-default parameters to keep URLs short
 */
// Geometry fields belonging to the other transmission-line types (microstrip /
// diff / gcpw / stripline and their solder-mask, top-dielectric and ground-cutout
// options). For broadside coupled stripline these are excluded from the link; it
// has its own bs_* fields. The top-ground enclosure controls are excluded too —
// broadside's top/bottom grounds are intrinsic, so only the side-wall enclosure
// applies (see buildSolverFromParams). Everything not listed here (frequency,
// s-params, side enclosure, plating, surface roughness, solver, modes, material
// and interpolation options) is shared and round-trips when non-default.
const BROADSIDE_EXCLUDED_KEYS = new Set([
    'w', 'h', 't', 'er', 'tand', 'sigma',
    'trace_spacing', 'gap', 'via_gap', 'gnd_width',
    'stripline_top_h', 'er_top', 'tand_top',
    'use_sm', 'sm_t_sub', 'sm_t_trace', 'sm_t_side', 'sm_er', 'sm_tand',
    'use_top_diel', 'top_diel_h', 'top_diel_er', 'top_diel_tand',
    'use_gnd_cut', 'gnd_cut_w', 'gnd_cut_h',
    'use_top_gnd', 'enclosure_height',
]);

// Coax: a fixed cross-section with no board stackup, so no stackup option applies. Plating,
// surface roughness, solver and frequency settings are shared and DO round-trip, so they
// are deliberately absent here (as is mesh_backend — a coax link must carry it).
const COAX_EXCLUDED_KEYS = new Set([
    ...BROADSIDE_EXCLUDED_KEYS,
    'use_enclosure', 'enclosure_width', 'use_side_gnd',
]);

// Rectangular waveguide: the walls are the domain boundary, exactly as for coax, so the
// same board-stackup and enclosure options are meaningless. The interpolating sweep is
// excluded too, it is forced off for this type (see toggleParameterVisibility).
const WAVEGUIDE_EXCLUDED_KEYS = new Set([
    ...COAX_EXCLUDED_KEYS,
    // Both are forced/ignored for this type: the sweep is analytic after one eigensolve
    // so interpolation is off, and the S-parameters are referenced to the modal impedance
    // rather than to sparam_z_ref.
    'interp_sweep', 'interp_tolerance', 'sparam_z_ref',
]);

// Each type's own geometry keys, so a link for one type never carries another's.
const TYPE_ONLY_KEYS = {
    broadside_stripline: Object.keys(DEFAULT_SETTINGS).filter(k => k.startsWith('bs_')),
    coax: Object.keys(DEFAULT_SETTINGS).filter(k => k.startsWith('coax_')),
    rect_waveguide: Object.keys(DEFAULT_SETTINGS).filter(k => k.startsWith('wg_')),
    custom: Object.keys(DEFAULT_SETTINGS).filter(k => k.startsWith('custom_')),
};
const EXCLUDED_BY_TYPE = {
    broadside_stripline: BROADSIDE_EXCLUDED_KEYS,
    coax: COAX_EXCLUDED_KEYS,
    rect_waveguide: WAVEGUIDE_EXCLUDED_KEYS,
    // The text carries the whole stackup, the boundaries and the plated faces. Model
    // Thick Plating stays: it applies to every plated conductor.
    custom: new Set([...COAX_EXCLUDED_KEYS, 'use_plating', 'plating_sigma', 'plating_t', 'plating_rq', 'plating_rq_iface',
        'plating_top', 'plating_sides', 'plating_bottom']),
};

// The settings that differ from the defaults, without the keys of other line types.
function nonDefaultSettings(settings) {
    // Exclusions only ever SHORTEN a newly generated URL; settingsFromURL merges over
    // defaults and ignores unknown keys, so existing links are unaffected.
    const excluded = EXCLUDED_BY_TYPE[settings.tl_type];
    const otherTypeKeys = new Set();
    for (const [type, keys] of Object.entries(TYPE_ONLY_KEYS)) {
        if (type === settings.tl_type) continue;
        for (const k of keys) otherTypeKeys.add(k);
    }

    // Filter out default values (and the other types' geometry fields)
    const out = {};
    for (const key in settings) {
        if (excluded && excluded.has(key)) continue;
        if (otherTypeKeys.has(key)) continue;

        const value = settings[key];
        const defaultValue = DEFAULT_SETTINGS[key];

        // Include if value differs from default
        // Handle NaN specially (NaN !== NaN is true, so we need special comparison)
        const bothNaN = (typeof value === 'number' && isNaN(value)) &&
                        (typeof defaultValue === 'number' && isNaN(defaultValue));

        if (bothNaN) {
            // Both NaN, skip (it's the default)
            continue;
        } else if (value !== defaultValue) {
            out[key] = value;
        }
    }

    return out;
}

// Settings as a URL-safe base64 string.
function settingsToURL(settings) {
    return btoa(encodeURIComponent(JSON.stringify(nonDefaultSettings(settings))));
}

/**
 * Deserialize settings from URL-safe base64 string
 */
function settingsFromURL(encoded) {
    try {
        const json = decodeURIComponent(atob(encoded));
        return JSON.parse(json);
    } catch (e) {
        log(`Failed to parse URL parameters: ${(e && e.message) || e}`);
        return null;
    }
}

// Show the option section of each advanced checkbox to match its state. Restoring a
// setting sets .checked without firing change, so the sections need this afterwards.
function syncCheckboxSections() {
    [['chk_solder_mask', 'solder-mask-params'], ['chk_top_diel', 'top-diel-params'],
     ['chk_gnd_cut', 'gnd-cut-params'], ['chk_enclosure', 'enclosure-params'],
     ['chk_plating', 'plating-params']].forEach(([id, sectionId]) => {
        const checkbox = document.getElementById(id);
        const section = document.getElementById(sectionId);
        if (checkbox && section) {
            section.style.display = checkbox.checked ? 'block' : 'none';
        }
    });
    // An enclosure is a physical boundary: the Modes tab cannot shrink it.
    document.getElementById('modes-shrink-domain').disabled = document.getElementById('chk_enclosure').checked;
}

/**
 * Restore UI settings from a settings object
 * Merges with defaults to explicitly set all values, preventing browser-remembered inputs
 */
function restoreSettings(settings) {
    if (!settings) return false;

    try {
        // Merge with defaults - URL settings override defaults
        // We explicitly set ALL values to prevent browser-remembered inputs
        const fullSettings = { ...DEFAULT_SETTINGS, ...settings };

        // Helper to restore value with unit
        const setValueWithUnit = (id, value) => {
            const element = document.getElementById(id);
            if (!element || value === undefined || value === null) return;
            // NaN is the "auto"/"full" default of a placeholder input: clear it so a
            // browser-remembered value does not survive the restore
            if (isNaN(value)) { element.value = ''; return; }
            // Format number to remove floating point artifacts
            const formattedValue = parseFloat(value.toPrecision(12));
            const unit = window.getDefaultUnit ? window.getDefaultUnit(id) : '';
            if (unit && element.classList.contains('unit-input')) {
                element.value = `${formattedValue} ${unit}`;
            } else {
                element.value = formattedValue;
            }
        };

        // Set input values - now always from fullSettings to override browser memory
        const tlTypeSelect = document.getElementById('tl_type');
        tlTypeSelect.value = fullSettings.tl_type;
        // Trigger change event to update UI visibility for the selected transmission line type
        tlTypeSelect.dispatchEvent(new Event('change', { bubbles: true }));

        const meshBackendSelect = document.getElementById('mesh_backend');
        if (meshBackendSelect && fullSettings.mesh_backend) {
            // Map legacy saved values onto current options ('triangular' and
            // 'fullwave_occ' predate the current dropdown).
            const legacy = { triangular: 'fullwave_mqs', fullwave_occ: 'fullwave_mqs' };
            const v = legacy[fullSettings.mesh_backend] ?? fullSettings.mesh_backend;
            meshBackendSelect.value = v;
            // An unknown value leaves the select EMPTY (selectedIndex -1) and
            // getParams() would then silently fall back to rectilinear — pin to
            // the first (default) option instead.
            if (meshBackendSelect.value !== v) meshBackendSelect.selectedIndex = 0;
            // Fire change so dependent UI (backend-vs-line-type enforcement) updates.
            meshBackendSelect.dispatchEvent(new Event('change', { bubbles: true }));
        }

        for (const [key, id, kind] of fieldsFor('s')) {
            const value = fullSettings[key];
            if (kind === 'unit') { setValueWithUnit(id, value); continue; }
            if (key === 'min_converged_passes' && value === undefined) continue;
            const el = document.getElementById(id);
            if (kind === 'chk') el.checked = !!value;
            else if (kind === 'bool') el.checked = value !== false;
            // Stored as a fraction (settings/link format), displayed in percent.
            // toPrecision strips float noise (0.003 * 100 = 0.30000000000000004).
            else if (kind === 'pct') el.value = parseFloat((100 * value).toPrecision(10));
            else el.value = value;
        }
        setCustomGeometryText(fullSettings.custom_geom || '');

        // tl_type is restored above the backend dropdown and the sweep checkboxes, so both
        // locks have to run once every input is in place, otherwise a stale or
        // hand-edited link can leave a full-wave-only type selected with the quasi-static
        // backend (which cannot mesh it), or leave the interpolating sweep ticked on a
        // medium that does not support it.
        if (window.enforceBackendForType) window.enforceBackendForType();
        if (window.enforceSweepOptionsForType) window.enforceSweepOptionsForType();
        syncCheckboxSections();
        syncLogPanel();

        return true;
    } catch (e) {
        console.error('Failed to restore settings:', e);
        return false;
    }
}

// Custom geometry links carry the geometry text, so they go into the URL fragment (never
// sent to the server, no length limit on its side) as deflate + base64url with a "z."
// prefix. Every other type keeps the original ?params= form, so existing links and the
// links of the fixed types do not change.
async function streamBytes(stream) {
    return new Uint8Array(await new Response(stream).arrayBuffer());
}

async function encodeLongSettings(settings) {
    const json = JSON.stringify(settings);
    if (typeof CompressionStream === 'undefined') return btoa(encodeURIComponent(json));
    const bytes = await streamBytes(new Blob([json]).stream().pipeThrough(new CompressionStream('deflate-raw')));
    let bin = '';
    for (const v of bytes) bin += String.fromCharCode(v);
    return 'z.' + btoa(bin).replace(/\+/g, '-').replace(/\//g, '_').replace(/=+$/, '');
}

// Largest decompressed link accepted. A geometry text is a few kB, while a crafted link
// of a few kB can inflate to gigabytes.
const MAX_LINK_BYTES = 1 << 20;

// Reads a stream into bytes, giving up past `limit` bytes.
async function streamBytesCapped(stream, limit) {
    const reader = stream.getReader(), chunks = [];
    let n = 0;
    for (;;) {
        const { done, value } = await reader.read();
        if (done) break;
        n += value.length;
        if (n > limit) { reader.cancel().catch(() => {}); throw new Error(`the link decompresses to more than ${limit} bytes`); }
        chunks.push(value);
    }
    const out = new Uint8Array(n);
    let o = 0;
    for (const c of chunks) { out.set(c, o); o += c.length; }
    return out;
}

async function decodeLongSettings(encoded) {
    if (!encoded.startsWith('z.')) return settingsFromURL(encoded);
    try {
        const b64 = encoded.slice(2).replace(/-/g, '+').replace(/_/g, '/');
        const bin = atob(b64 + '='.repeat((4 - b64.length % 4) % 4));
        const bytes = Uint8Array.from(bin, c => c.charCodeAt(0));
        const out = await streamBytesCapped(new Blob([bytes]).stream().pipeThrough(new DecompressionStream('deflate-raw')),
            MAX_LINK_BYTES);
        return JSON.parse(new TextDecoder().decode(out));
    } catch (e) {
        log(`Failed to parse URL parameters: ${(e && e.message) || e}`);
        return null;
    }
}

/**
 * Copy current settings as URL to clipboard
 */
async function copySettingsLink() {
    const settings = getUISettings();
    const base = `${window.location.origin}${window.location.pathname}`;
    let url;
    if (settings.tl_type === 'custom') {
        // Same default filtering as settingsToURL, then the compressed fragment form.
        url = `${base}#params=${await encodeLongSettings(nonDefaultSettings(settings))}`;
        if (url.length > 8000) {
            log(`The link is ${url.length} characters long. Some applications truncate long links, ` +
                `the geometry text file (Save file) is the safer way to share this geometry.`);
        }
    } else {
        url = `${base}?params=${settingsToURL(settings)}`;
    }

    navigator.clipboard.writeText(url).then(() => {
        const btn = document.getElementById('copy-link-btn');
        const originalText = btn.textContent;
        btn.textContent = 'Copied!';
        setTimeout(() => { btn.textContent = originalText; }, 2000);
    }).catch(err => {
        console.error('Failed to copy link:', err);
        // Fallback: show prompt with URL
        prompt('Copy this URL:', url);
    });
}

/**
 * Check URL for params and restore if present
 */
function loadSettingsFromURL() {
    const urlParams = new URLSearchParams(window.location.search);
    const paramsStr = urlParams.get('params');
    if (paramsStr) {
        const settings = settingsFromURL(paramsStr);
        if (settings && restoreSettings(settings)) {
            log('Settings restored from URL');
            return true;
        }
    }
    return false;
}

// The fragment form decodes asynchronously, so it is applied after init has drawn the
// default geometry.
async function loadSettingsFromFragment() {
    const m = /[#&]params=([^&]+)/.exec(window.location.hash);
    if (!m) return;
    const settings = await decodeLongSettings(m[1]);
    if (!settings || !restoreSettings(settings)) return;
    log('Settings restored from URL');
    // The link has been applied: without it in the address bar a reload keeps the
    // edits made since (the session copy of the geometry) instead of the link's text.
    history.replaceState(null, '', window.location.pathname + window.location.search);
    updateGeometry();
    draw(true);
    updateSweepParamList();
    autoFillSweepRange();
}

// Solver log.
//
// The naive version (`c.textContent += msg + '\n'; c.scrollTop = c.scrollHeight;`) had
// problems that made the log the slowest thing on the page during a solve:
//   * `textContent +=` reads the whole accumulated text, concatenates, and replaces the
//     single text node. O(n) per line, O(n²) over a run, and it also destroyed any text
//     the user had selected.
//   * reading `scrollHeight` forces a synchronous layout on every line.
//   * the unconditional scroll-to-bottom fought the user.
// Now: lines are queued and flushed once per animation frame (one layout per frame, not
// per line), appended as new text nodes, auto-scroll only when the view is already
// pinned to the bottom, and the buffer is capped so scrollback stays cheap.
const LOG_MAX_LINES = 5000;
let _logQueue = [], _logFlushScheduled = false, _logLineCount = 0;
// Line for the collapsed log bar in place of the last logged line (a result summary).
let _logStatus = null;

function _flushLog() {
    _logFlushScheduled = false;
    if (!_logQueue.length) return;
    const c = document.getElementById('console_out');
    if (!c) { _logQueue = []; return; }
    // Sample the scroll position before mutating: within 4px of the bottom counts as
    // "following the tail" (fractional device pixels never land exactly on 0).
    const atBottom = c.scrollHeight - c.scrollTop - c.clientHeight < 4;
    const chunk = _logQueue.join('');
    _logQueue = [];
    c.appendChild(document.createTextNode(chunk));
    const lines = chunk.split('\n').filter(l => l.trim());
    if (_logStatus !== null) setLogStatus(_logStatus);
    else if (lines.length) setLogStatus(lines[lines.length - 1].trim());
    _logStatus = null;
    _logLineCount += chunk.split('\n').length - 1;
    if (_logLineCount > LOG_MAX_LINES) {
        // Drop the oldest nodes wholesale rather than re-splitting text: keeping the
        // node count and total text bounded is what keeps scrollback responsive.
        while (c.childNodes.length > 1 && _logLineCount > LOG_MAX_LINES) {
            const first = c.firstChild;
            _logLineCount -= (first.textContent.split('\n').length - 1);
            c.removeChild(first);
        }
    }
    if (atBottom) c.scrollTop = c.scrollHeight;
}

// Frequency unit for a sweep, chosen from its highest frequency.
const FREQ_UNITS = [{ name: 'THz', scale: 1e12 }, { name: 'GHz', scale: 1e9 },
                    { name: 'MHz', scale: 1e6 }, { name: 'kHz', scale: 1e3 }, { name: 'Hz', scale: 1 }];
function freqUnit(fMax) {
    return FREQ_UNITS.find(u => fMax >= u.scale) ?? FREQ_UNITS[FREQ_UNITS.length - 1];
}
function formatFreq(f, unit) {
    return String(parseFloat((f / unit.scale).toPrecision(4)));
}

// Im(Zc) is printed when it is at least this fraction of Re(Zc). Dielectric loss alone
// gives Im/Re -> tand/2 at high frequency (1 % on FR4), conductor loss below a few
// hundred MHz gives a negative Im of several percent and more.
const ZC_IM_THRESHOLD = 0.02;
function formatZc(Zc, threshold = ZC_IM_THRESHOLD) {
    const re = Zc.re.toFixed(2);
    if (!(Math.abs(Zc.im) >= threshold * Math.abs(Zc.re))) return re;
    return `${re} ${Zc.im < 0 ? '−' : '+'} j${Math.abs(Zc.im).toFixed(2)}`;
}

// `status` replaces the last line of `msg` in the collapsed log bar.
function log(msg, status = null) {
    _logQueue.push(msg + '\n');
    _logStatus = status;
    if (_logFlushScheduled) return;
    _logFlushScheduled = true;
    // rAF coalesces a burst into one layout. It does not fire in a background tab, so
    // fall back to a timer there rather than letting the queue grow unboundedly.
    if (typeof requestAnimationFrame === 'function' && document.visibilityState !== 'hidden')
        requestAnimationFrame(_flushLog);
    else setTimeout(_flushLog, 100);
}

// Going to the background does not cancel an already-scheduled rAF, it just
// stops it from ever firing so a flush scheduled while visible would strand the
// queue until the tab came back. Drain it on the transition. Lines queued from
// then on pick the timer path above by themselves.
document.addEventListener('visibilitychange', () => {
    if (document.visibilityState === 'hidden' && _logFlushScheduled) _flushLog();
});

// Solve worker client
//
// Every solve runs in solve_worker.js so a multi-second WASM call cannot freeze the tab.
// This is the RPC layer: one long-lived module worker, jobs keyed by id, streamed
// messages routed to per-job handlers.
//
// The worker is created lazily on the first solve rather than at load, so the page (and
// the geometry preview, which needs no worker) is interactive immediately and a browser
// that cannot construct a module worker fails at a point where we can report it.
let _worker = null, _workerJobId = 0;
const _workerJobs = new Map();   // id -> { resolve, reject, handlers }

// Fail every outstanding job. Used on the worker-level failures that can never produce a
// 'done'/'error' reply of their own. Without it the UI would sit in "solving" forever.
function failAllJobs(message) {
    const err = new Error(message);
    for (const [id, job] of _workerJobs) { _workerJobs.delete(id); job.reject(err); }
}

function getWorker() {
    if (_worker) return _worker;
    // The literal "solve_worker.js" is rewritten by build.sh to add a content hash. Keep
    // it a plain string so that sed can find it.
    const w = new Worker(new URL("solve_worker.js", import.meta.url), { type: 'module' });
    _worker = w;
    _worker.onmessage = (e) => {
        const m = e.data;
        // Log lines are printed regardless of whether the job is still being awaited (a
        // stopped job can still emit its trailing messages).
        if (m.type === 'log') { log(m.msg); return; }
        const job = _workerJobs.get(m.id);
        if (!job) return;
        if (m.type === 'done') { _workerJobs.delete(m.id); job.resolve(m); return; }
        if (m.type === 'error') { _workerJobs.delete(m.id); job.reject(new Error(m.message)); return; }
        const h = job.handlers[m.type];
        if (h) h(m);
    };
    _worker.onerror = (e) => {
        // A worker-level failure (module load error, uncaught throw) never resolves the
        // outstanding job, so fail them all rather than hanging the UI in "solving".
        failAllJobs(e.message || 'Solver worker failed');
        // Treat the instance as unusable and let the next solve build a fresh one. It has
        // to be terminated explicitly: dropping the reference does not stop a worker, and
        // this one is holding a multi-hundred-MB WASM heap.
        w.terminate();
        if (_worker === w) _worker = null;
        // The solve kept for plots died with the worker.
        plotFieldsSolver = null; plotFieldsJob = null; plotFieldsSolveKey = null;
    };
    _worker.onmessageerror = () => {
        // A reply that could not be deserialized (structured clone refused something in a
        // result). The job it belonged to is unidentifiable — its id was in the message —
        // so nothing can ever settle it. The worker itself is still healthy, so keep it.
        failAllJobs('Solver worker sent a result that could not be read.');
    };
    return _worker;
}

function workerJob(type, payload, handlers = {}) {
    const w = getWorker();
    const id = ++_workerJobId;
    return new Promise((resolve, reject) => {
        _workerJobs.set(id, { resolve, reject, handlers });
        w.postMessage({ id, type, ...payload });
    });
}

// Cancellation is cooperative: the worker polls this between WASM calls, exactly where
// the old main-thread code polled `stopRequested`. Stop targets the solve or sweep job by
// id, so it also holds when that job is still queued behind a plot job.
let _stoppableJob = null;
function workerStop() {
    if (_worker && _stoppableJob !== null) _worker.postMessage({ type: 'stop', job: _stoppableJob });
}

// Liveness indicator.
//
// A long solve gives no sign of life between phase messages. A ticking elapsed
// time distinguishes the two at a glance.
let _hbTimer = null, _hbStart = 0, _hbPhase = '', _hbEl = null;

function _hbFormat(ms) {
    const s = Math.floor(ms / 1000);
    return s < 60 ? `${s}s` : `${Math.floor(s / 60)}m ${String(s % 60).padStart(2, '0')}s`;
}

function _hbRender() {
    if (!_hbEl) return;
    const t = _hbFormat(Date.now() - _hbStart);
    _hbEl.textContent = _hbPhase ? `${_hbPhase} · ${t}` : t;
}

function heartbeatStart(el) {
    heartbeatStop();
    _hbEl = el; _hbStart = Date.now(); _hbPhase = '';
    _hbRender();
    _hbTimer = setInterval(_hbRender, 1000);
}

// Phase text from the worker; the elapsed clock keeps running underneath it.
function heartbeatPhase(text) {
    _hbPhase = text || '';
    _hbRender();
}

function heartbeatStop(finalText = null) {
    if (_hbTimer) { clearInterval(_hbTimer); _hbTimer = null; }
    if (_hbEl && finalText !== null) _hbEl.textContent = finalText;
    _hbEl = null;
}

// Complex survives structured clone only as a bare {re, im}, the prototype is dropped.
// sparameters.js does genuine complex arithmetic with mode.Zc, so rebuild the instances
// on arrival here.
function reviveResult(result) {
    if (result && Array.isArray(result.modes)) {
        for (const m of result.modes) {
            if (m.Zc && !(m.Zc instanceof Complex)) m.Zc = new Complex(m.Zc.re, m.Zc.im);
        }
    }
    return result;
}
const reviveSweep = (rows) => { for (const r of rows) reviveResult(r.result); return rows; };

// Graft a worker solve's field arrays onto the main thread's view-model solver, so
// plot.js keeps reading solver.x/y/V/Ex/Ey/triMesh unchanged.
function applyFields(target, fields) {
    if (!target || !fields) return;
    target.x = fields.x;
    target.y = fields.y;
    target.V = fields.V;
    target.Ex = fields.Ex;
    target.Ey = fields.Ey;
    target.ExIm = fields.ExIm || null;
    target.EyIm = fields.EyIm || null;
    target.surfaceK = fields.surfaceK || null;
    target.currentJ = fields.currentJ || null;
    target.currentMesh = fields.currentMesh || null;
    target.fieldMesh = fields.fieldMesh || null;
    target.surfaceKSource = fields.surfaceKSource || null;
    target.idealGrounds = fields.idealGrounds || null;
    target.fieldFreq = fields.fieldFreq ?? null;
    target.fieldKind = fields.fieldKind || 'static';
    if (fields.triMesh) target.triMesh = fields.triMesh;
    target.solution_valid = true;
    target.mesh_generated = true;
}

// Field plot frequency in Hz from Plot Options, null when empty (the highest solved).
function plotFrequency() {
    const el = document.getElementById('plot-freq');
    if (!el || !el.value.trim()) return null;
    const f = getInputValue('plot-freq');
    return isFinite(f) && f >= 0 ? f : null;
}

// The frequency the plot fields of the kept solve are drawn at: the Plot Options one,
// limited like the worker limits it (full-wave plots stay within the solve range).
// { f, limited }, limited when the Plot Options frequency was above the limit.
function keptPlotFrequency() {
    const want = plotFrequency() ?? plotFieldsMaxFreq;
    const f = Math.min(want, solver.plotFreqLimit(plotFieldsMaxFreq));
    return { f, limited: f < want };
}

function logPlotLimit(f) {
    log(`Plot frequency is above the solved range: full-wave fields are plotted at most at `
        + `${(f / 1e9).toPrecision(4)} GHz. Raise the stop frequency to plot higher.`);
}

// Re-solve the plot fields of the last solve at the Plot Options frequency. Needs the
// worker's retained solve of the current geometry; requests arriving while one runs
// collapse into one more run after it. afterSolve: called when a solve ends, which
// already logged the plot frequency limit.
async function updatePlotFields(afterSolve = false) {
    if (!solver || solver !== plotFieldsSolver || !solver.solution_valid) return;
    if (isSimulating || isSweeping || isSolvingModes) return;
    if (plotFieldsBusy) { plotFieldsPending = true; return; }
    plotFieldsBusy = true;
    try {
        do {
            plotFieldsPending = false;
            const { f, limited } = keptPlotFrequency();
            if (f === solver.fieldFreq) {
                if (limited && afterSolve !== true) logPlotLimit(f);
                continue;
            }
            const target = solver;
            const { fields } = await workerJob('plotFields', { freq: f });
            if (fields && target === solver && target === plotFieldsSolver) {
                applyFields(target, fields);
                draw();
            }
        } while (plotFieldsPending);
    } catch (e) {
        log('Plot field update failed: ' + (e.message || e));
    } finally {
        plotFieldsBusy = false;
    }
}

// Everything a simulate job depends on except the plot frequency. A Solve with the same
// key as the kept solve only needs new plot fields.
function solveInputKey() {
    try {
        return JSON.stringify({
            params: getParams(),
            frequencies: getFrequencies(),
            interp: !!document.getElementById('chk_interp_sweep')?.checked,
            interpTol: document.getElementById('interp_tolerance')?.value,
        });
    } catch (e) {
        return null;
    }
}

// Solve button: when only the plot frequency changed since the kept solve, update the
// plot fields instead of solving again. Returns true when it handled the click.
function replotInsteadOfSolve() {
    if (!solver || solver !== plotFieldsSolver || !solver.solution_valid) return false;
    // A running solve: runSimulation reports it.
    if (isSimulating || isSweeping || isSolvingModes) return false;
    const key = solveInputKey();
    if (!key || key !== plotFieldsSolveKey) return false;
    const { f, limited } = keptPlotFrequency();
    if (f === solver.fieldFreq) {
        if (limited) logPlotLimit(f);
        else log('Nothing changed since the last solve.');
    } else {
        log(`Only the plot frequency changed: updating the field plots at ${(f / 1e9).toPrecision(4)} GHz.`);
        updatePlotFields();
    }
    return true;
}

function getFrequencies() {
    const start = getInputValue('freq-start');
    const stop = getInputValue('freq-stop');
    let points = parseInt(document.getElementById('freq-points').value);

    // Validate points - default to 1 if invalid
    if (isNaN(points) || points < 1) {
        points = 1;
        document.getElementById('freq-points').value = '1';
    }

    const freqs = [];
    if (points === 1) {
        // Single frequency point - use start frequency
        freqs.push(start);
    } else {
        // Multiple points - linear spacing
        for (let i = 0; i < points; i++) {
            freqs.push(start + (stop - start) * i / (points - 1));
        }
    }
    return freqs;
}

// Interpolating-sweep tolerance, as a fraction. Validated here (on the thread that owns
// the input) so a bad value fails before the worker job starts.
function interpTolerance() {
    const tolPercent = parseFloat(document.getElementById('interp_tolerance')?.value);
    if (isNaN(tolPercent) || tolPercent <= 0) {
        throw new Error("Interpolation tolerance must be a positive number.");
    }
    return tolPercent / 100;
}

// Keys of getParams() that define the geometry, hashed for change tracking.
const GEOMETRY_HASH_KEYS = ['tl_type', 'custom_geom', 'custom_overrides', ...fieldsFor('h').map(f => f[0])];

function getGeometryHash() {
    const p = getParams();
    return JSON.stringify(Object.fromEntries(GEOMETRY_HASH_KEYS.map(k => [k, p[k]])));
}

/**
 * Get a hash of frequency parameters for change tracking
 */
function getFrequencyHash() {
    return JSON.stringify({
        freq_start: getInputValue('freq-start'),
        freq_stop: getInputValue('freq-stop'),
        freq_points: parseInt(document.getElementById('freq-points').value)
    });
}

/**
 * Update notices on Results and S-parameters tabs
 */
function updateResultNotices() {
    // Shows the notice with this text, or hides it when the text is null.
    const setNotice = (id, text) => {
        const notice = document.getElementById(id);
        if (!notice) return;
        if (text) document.getElementById(`${id}-text`).textContent = text;
        notice.style.display = text ? 'block' : 'none';
    };
    const hasResults = frequencySweepResults && frequencySweepResults.length > 0;
    let resultsText = null, sparamText = null;
    // exportTitle stays undefined when the button keeps its current title.
    let exportable = false, exportTitle;

    if (!hasResults) {
        resultsText = 'No results available. Run solver to view results.';
        sparamText = 'No results available. Run solver to view S-parameters.';
    } else {
        const geometryChanged = lastSolvedGeometry && getGeometryHash() !== lastSolvedGeometry;
        const frequencyChanged = lastSolvedFrequency && getFrequencyHash() !== lastSolvedFrequency;
        if (!isSimulating && (geometryChanged || frequencyChanged)) {
            // Keep the old results visible under the notice.
            resultsText = sparamText =
                `${geometryChanged ? 'Geometry' : 'Frequency'} changed. Solve to update results.`;
            exportTitle = 'Cannot export - geometry or frequency changed';
        } else {
            // A self-referenced medium drops its below-cutoff points (they have an
            // imaginary modal impedance), so the S-parameter tab can legitimately end up
            // with nothing to draw while the Results tab is full. Say why, rather than
            // leaving an empty plot and an export button that only fails when clicked.
            exportable = usableSweepPoints(frequencySweepResults).length > 0;
            if (!exportable) {
                sparamText = 'Every sweep point is below the cutoff — ' +
                    'the mode is evanescent there, so there are no propagating ' +
                    'S-parameters. Raise the sweep frequency above the cutoff.';
            }
            exportTitle = exportable ? '' : 'Cannot export - no propagating sweep points';
        }
    }

    setNotice('results-notice', resultsText);
    setNotice('sparam-notice', sparamText);
    const exportBtn = document.getElementById('export-snp');
    if (exportBtn) {
        exportBtn.disabled = !exportable;
        if (exportTitle !== undefined) exportBtn.title = exportTitle;
    }
    // The differential-mode checkboxes apply to two-mode results only.
    const differential = hasResults && frequencySweepResults[0].result.modes.length === 2;
    for (const id of ['results-diff', 'sparam-diff']) {
        const el = document.getElementById(id);
        if (el) el.disabled = !differential;
    }
}

function switchTab(tabName) {
    document.querySelectorAll('.tab-button').forEach(btn =>
        btn.classList.toggle('active', btn.dataset.tab === tabName));
    document.querySelectorAll('.tab-content').forEach(div =>
        div.classList.toggle('active', div.id === `tab-${tabName}`));
    currentTab = tabName;

    if (tabName === 'results') {
        updateResultNotices();
        if (frequencySweepResults) {
            drawResultsPlot();
        }
    } else if (tabName === 'sparams') {
        updateResultNotices();
        if (frequencySweepResults) {
            drawSParamPlot();
        }
    } else if (tabName === 'sweep') {
        updateSweepParamList();
        updateSweepDiffCheckbox();
        updateSweepNotice();
        redrawSweepPlot();
    } else if (tabName === 'geometry') {
        // Refresh the geometry plot when switching back
        draw();
    } else if (tabName === 'modes') {
        const Plotly = getPlotly();
        const c = document.getElementById('modes-plot');
        updateModesNotice();
        if (Plotly && c) {
            if (modesResult) {
                // Already solved: keep the field plot; just resize to the now-visible tab.
                if (c.data) Plotly.Plots.resize(c);
            } else {
                // Nothing solved yet: (re)build an independent solver from the current sidebar
                // geometry and show a live geometry preview of the structure.
                modesSolver = buildSolverFromParams(getParams());
                plotModesGeometry();
            }
        }
    }
}

// ===================== Modes Tab =====================
// Full-wave eigenmode viewer: solve the cross-section's modes at a chosen frequency,
// list them, and plot the selected mode's transverse E-field. Always uses the
// triangular full-wave backend (independent of the sidebar Solver dropdown).

async function runModesSolve() {
    if (isSolvingModes) return;
    if (isSimulating || isSweeping) { setModesStatus('Cannot solve modes while a simulation is running.'); return; }
    flushCustomGeometryEdits();

    const freq = getInputValue('modes-freq');
    let nev = parseInt(document.getElementById('modes-nev').value);
    if (!isFinite(freq) || freq <= 0) { setModesStatus('Enter a valid frequency.'); return; }
    if (!isFinite(nev) || nev < 1) { nev = 1; document.getElementById('modes-nev').value = '1'; }
    nev = Math.min(nev, 30);
    let meshDensity = parseInt(document.getElementById('modes-mesh-density').value);
    if (!isFinite(meshDensity)) meshDensity = 8;
    meshDensity = Math.min(Math.max(meshDensity, 3), 40);
    document.getElementById('modes-mesh-density').value = String(meshDensity);

    // Build an INDEPENDENT solver from the current sidebar geometry. The Modes tab keeps its
    // own solver so solving modes never disturbs the main solve (Geometry field plot, Results,
    // S-parameters all stay intact).
    modesSolver = buildSolverFromParams(getParams());
    if (!modesSolver) { setModesStatus('Geometry is invalid — check parameters.'); return; }
    // Stamp what is being solved now, not what the sidebar reads when the worker
    // returns: the UI stays interactive during the eigensolve, so an edit made while it
    // runs must still trip the staleness notice.
    const solvedGeometry = getGeometryHash();

    const btn = document.getElementById('btn-solve-modes');
    isSolvingModes = true;
    btn.disabled = true;
    modesResult = null; modesSelectedIdx = -1;
    setModesStatus('Meshing…');
    heartbeatStart(document.getElementById('modes-status'));
    log(`Solving ${nev} modes at ${formatValueWithUnit(freq, 'GHz')}…`);
    // The modes solve treats open walls as radiating ABCs — same clearance concern.
    for (const w of modesSolver.openBoundaryWarnings()) log(`⚠ Warning: ${w}`);

    // Drive the adaptive mesh refinement with the sidebar Solver Settings, same as the
    // main solve (Max Nodes is entered in thousands).
    const p = getParams();
    const refineOpts = {
        maxRefineIters: p.max_iters,
        refineTol: p.tolerance,
        maxNodes: p.max_nodes * 1000,
        minConvergedPasses: p.min_converged_passes,
        // Never verify here: the certificate covers only the quasi-TEM static
        // solution, not the higher-order modes this tab displays, so it would add
        // bisected-mesh solves (the tab is slow already) without certifying anything
        // the user sees.
        certify: false,
        wavelengthDensity: meshDensity,
        // Only an auto-sized open domain can be shrunk. An enclosure is a physical
        // boundary (the checkbox is disabled in the UI while an enclosure is on).
        shrinkDomain: document.getElementById('modes-shrink-domain').checked
            && !document.getElementById('chk_enclosure').checked,
    };

    try {
        // The eigensolve is the single worst blocking call in the app, so it
        // runs in the worker. modesFieldCache is dropped here because the
        // worker's mode grids belong to the solve that just started.
        modesFieldCache = new Map();
        const { result } = await workerJob('modes', { params: p, freq, nev, refineOpts }, {
            progress: (m) => { if (m.text) heartbeatPhase(m.text); },
        });
        heartbeatStop();
        modesResult = result;
        // Record what the displayed modes were solved at, so the staleness notice can flag
        // a later geometry/frequency change.
        lastModesGeometry = solvedGeometry;
        lastModesFrequency = freq;
        updateModesNotice();
        showModesWarnings(result.warnings || []);
        if (result.error) {
            setModesStatus('Solve failed: ' + result.error);
            log('Modes solve error: ' + result.error);
        } else if (!result.modes || result.modes.length === 0) {
            setModesStatus('No modes converged. Try increasing the number of modes.');
        } else {
            const nProp = result.modes.filter(m => m.status === 'propagating').length;
            setModesStatus(`${result.modes.length} converged, ${nProp} propagating (N=${result.N}, ${result.nTris} triangles).`);
            log(`Found ${result.modes.length} modes (${nProp} propagating).`);
        }
        renderModesTable();
        // Auto-select the quasi-TEM mode (highest overlap with the static drive).
        if (modesResult && modesResult.modes && modesResult.modes.length) {
            let best = 0, bestOvl = -1;
            modesResult.modes.forEach((m, i) => { if (m.overlap > bestOvl) { bestOvl = m.overlap; best = i; } });
            selectMode(best, true);   // fresh solve → focus the view on the structure
        }
    } catch (e) {
        setModesStatus('Solve failed: ' + (e.message || e));
        log('Modes solve exception: ' + (e.message || e));
        console.error(e);
    } finally {
        heartbeatStop();
        isSolvingModes = false;
        btn.disabled = false;
    }
}

function setModesStatus(msg) {
    const el = document.getElementById('modes-status');
    if (el) el.textContent = msg;
}

// Solve-level warnings from solveModes (e.g. floating grounds), shown above the mode
// list and logged. Cleared when there are none.
function showModesWarnings(warnings) {
    const box = document.getElementById('modes-warning');
    const text = document.getElementById('modes-warning-text');
    if (!box || !text) return;
    if (!warnings.length) { box.style.display = 'none'; text.textContent = ''; return; }
    text.textContent = warnings.map(w => w.message).join(' ');
    box.style.display = 'block';
    for (const w of warnings) log(`⚠ Warning: ${w.message}`);
}

// Show a staleness warning on the Modes tab (mirroring Results / S-parameters) when the
// geometry or the modes frequency has changed since the displayed modes were solved, so the
// user knows the field plot is out of date. Keeps the old plot visible, like the other tabs.
function updateModesNotice() {
    const notice = document.getElementById('modes-notice');
    const noticeText = document.getElementById('modes-notice-text');
    if (!notice || !noticeText) return;
    if (!modesResult || !lastModesGeometry) { notice.style.display = 'none'; return; }
    const geometryChanged = getGeometryHash() !== lastModesGeometry;
    const freqChanged = lastModesFrequency != null && getInputValue('modes-freq') !== lastModesFrequency;
    if (isSolvingModes || !(geometryChanged || freqChanged)) {
        notice.style.display = 'none';
        return;
    }
    noticeText.textContent = geometryChanged
        ? 'Geometry changed. Solve Modes to update.'
        : 'Frequency changed. Solve Modes to update.';
    notice.style.display = 'block';
}

function renderModesTable() {
    const container = document.getElementById('modes-list');
    if (!container) return;
    if (!modesResult || !modesResult.modes || !modesResult.modes.length) {
        container.innerHTML = '<p style="color:var(--text-muted); font-size:0.85em; padding:8px;">No modes to display.</p>';
        return;
    }
    const STATUS_LABEL = { propagating: 'PROP', near_cutoff: 'evan?', evanescent: 'evan', spurious: 'spurious', nullspace: 'null' };
    let html = '<table><thead><tr><th>#</th><th>status</th><th>&epsilon;<sub>eff</sub></th><th>overlap</th></tr></thead><tbody>';
    modesResult.modes.forEach((m, i) => {
        const eeff = m.eps_eff != null ? m.eps_eff.toFixed(4) : '–';
        const star = m.overlap > 0.5 ? ' *' : '';
        // "#" = position in the energy-sorted list (sequential), matching the plot title;
        // the raw eigensolve index (m.idx) is not meaningful to the user.
        html += `<tr class="${m.status}" data-idx="${i}">` +
            `<td>${i}</td><td class="status">${STATUS_LABEL[m.status] || m.status}${star}</td>` +
            `<td>${eeff}</td><td>${m.overlap.toFixed(3)}</td></tr>`;
    });
    html += '</tbody></table>';
    container.innerHTML = html;
    container.querySelectorAll('tr[data-idx]').forEach(tr =>
        tr.addEventListener('click', () => selectMode(parseInt(tr.dataset.idx))));
    highlightSelectedMode();
}

function highlightSelectedMode() {
    document.querySelectorAll('#modes-list tr[data-idx]').forEach(tr =>
        tr.classList.toggle('selected', parseInt(tr.dataset.idx) === modesSelectedIdx));
}

// resetView=true recomputes the focused view (used right after a fresh solve); when
// switching between modes we keep the user's current zoom/pan so the plot doesn't jump.
// The per-mode field grid lives in the worker (it is resampled from the FEM mesh on
// demand). Fetching is async and cached per index: with up to 30 modes, shipping every
// grid with the solve would cost far more than the few the user clicks through.
async function selectMode(idx, resetView = false) {
    if (!modesResult || !modesResult.modes || !modesResult.modes[idx]) return;
    modesSelectedIdx = idx;
    highlightSelectedMode();
    let grid = modesFieldCache.get(idx);
    if (grid === undefined) {
        // A modes solve started while this fetch waits replaces the cache, and the
        // field that comes back belongs to the solve before it.
        const cache = modesFieldCache;
        try {
            ({ grid } = await workerJob('modeField', { idx }));
        } catch (e) {
            log('Mode field fetch failed: ' + (e.message || e));
            grid = null;
        }
        // Cache successes only. Remembering a failure would leave the row permanently
        // blank with no way to retry short of re-solving. Re-asking on the next click
        // costs one worker round trip and lets a transient failure heal itself.
        if (cache !== modesFieldCache) return;
        if (grid) cache.set(idx, grid);
        // A click on another row while this one was in flight wins. Do not paint over it.
        if (modesSelectedIdx !== idx) return;
    }
    plotModesField(grid, modesResult.modes[idx], idx, resetView);
}

function plotModesField(grid, mode, idx, resetView = false) {
    const Plotly = getPlotly();
    const container = document.getElementById('modes-plot');
    if (!Plotly || !container || !modesSolver) return;
    if (!grid) {
        Plotly.purge(container);
        container._modesMesh = null;
        container._modesRelayoutBound = false;
        container.innerHTML = '<p style="color:var(--text-muted); padding:12px;">No field available for this mode.</p>';
        return;
    }

    const maxY = displayTop(modesSolver);
    const M = grid.mesh;

    const eeff = mode.eps_eff != null ? `, ε_eff=${mode.eps_eff.toFixed(3)}` : '';
    const STATUS_LABEL = { propagating: 'propagating', near_cutoff: 'near cutoff (evan?)', evanescent: 'evanescent', spurious: 'spurious', nullspace: 'null-space' };
    const title = `Mode ${idx} transverse |E| (${STATUS_LABEL[mode.status] || mode.status}${eeff})`;
    // A mode's eigenvector has an arbitrary scale, so absolute |E| is meaningless — show
    // the field normalized to its own maximum (0–1). This also keeps the colorbar tick
    // labels a constant width across modes, so switching modes never resizes the plot
    // area (and thus never shifts the equal-aspect x-range). The maximum is that of the
    // field resampled on a grid (grid.zmax): the mesh vertices sit on the conductor
    // corners, where the field is singular, and their maximum would dim the rest of the plot.
    const zmax = grid.zmax;
    const inv = zmax > 0 ? 1 / zmax : 1;
    const colorbar = { title: { text: '|E_t| / max', font: { color: '#aaa' } }, tickfont: { color: '#aaa' } };

    // |E| on the triangles: the color is an image (updateModesImage), the hover and
    // colorbar ride on invisible markers at the triangle centroids.
    container._modesMesh = { blocks: [{ tris: M.tris, Jv: M.E }], zmax: zmax || 1 };
    const traces = [centroidHoverTrace([M], (b, t) => triMeanE(b, t) * inv,
        'x: %{x:.3f} mm<br>y: %{y:.3f} mm<br>|E_t|/max: %{marker.color:.3f}<extra></extra>')];

    // Switching modes (a plot of the same kind exists, no view reset): update only the
    // field data and title in place — leaving the axes untouched preserves the current
    // zoom/pan exactly. (Re-issuing the layout would re-run the equal-aspect scaleanchor
    // solver and drift the x-range a little each time.)
    if (!resetView && container.data && container.data.length && container.data[0].type === traces[0].type) {
        const t = traces[0];
        Plotly.restyle(container, { x: [t.x], y: [t.y], 'marker.color': [t.marker.color] }, [0]);
        Plotly.relayout(container, { 'title.text': title });
        updateModesImage(container);
        return;
    }

    // Fresh plot (first render or after a new solve): full render at the focused view.
    const shapes = buildGeometryShapes(maxY);
    const view = computeModesView(maxY);
    const layout = modesPlotLayout(title, view, shapes);
    layout.coloraxis = { cmin: 0, cmax: 1, colorscale: 'Viridis', colorbar };
    // The image sits below the traces, under the grid lines.
    for (const ax of [layout.xaxis, layout.yaxis]) { ax.showgrid = false; ax.zeroline = false; }

    Plotly.react(container, traces, layout,
        { responsive: true, displayModeBar: true, scrollZoom: true, modeBarButtonsToRemove: ["select2d", "lasso2d"] });
    updateModesImage(container);
    if (!container._modesRelayoutBound) {
        container._modesRelayoutBound = true;
        // A zoom, pan or resize redraws the image for the new view.
        container.on('plotly_relayout', () => {
            if (container._modesMesh && !container._triImageUpdate) updateModesImage(container);
        });
    }
}

// Redraws the |E| image of the mode's triangles for the current axis ranges and plot size.
function updateModesImage(container) {
    updateTriImage(container, () => {
        const mm = container._modesMesh;
        return mm ? { blocks: mm.blocks, zmin: 0, zmax: mm.zmax, db: false, bleed: 2 } : { blocks: null };
    });
}

// Draw a geometry-only preview into the mode plot (before any solve), so the tab shows
// the structure like the reference viewer does.
function plotModesGeometry() {
    const Plotly = getPlotly();
    const container = document.getElementById('modes-plot');
    if (!Plotly || !container || !modesSolver || !modesSolver.conductors) return;
    container._modesMesh = null;
    const maxY = displayTop(modesSolver);
    const view = computeModesView(maxY);
    // Invisible scatter spanning the view so the axes (and shapes) scale correctly.
    const traces = [{ type: 'scatter', x: [view.xRange[0], view.xRange[1]], y: [view.yRange[0], view.yRange[1]],
        mode: 'markers', marker: { size: 0, opacity: 0 }, hoverinfo: 'skip', showlegend: false }];
    const layout = modesPlotLayout('Geometry (click Solve Modes to compute fields)', view, buildGeometryShapes(maxY));
    Plotly.react(container, traces, layout, { responsive: true, displayModeBar: true, scrollZoom: true, modeBarButtonsToRemove: ["select2d", "lasso2d"] });
}

// Focused view (mm) around the signal conductors — the geometry tab's zoom
// (computeGeometryView) with a whole-domain fallback when there are no signals.
function computeModesView(maxY) {
    return computeGeometryView(modesSolver, maxY)
        ?? { xRange: [-modesSolver.domain_width * 500, modesSolver.domain_width * 500], yRange: [0, maxY * 1000] };
}

// Plotly rectangle shapes for the geometry, reusing the geometry tab's shape builders so
// the two tabs render identically (incl. plated-edge markers). Everything sits ABOVE the
// heatmap: dielectrics faint (the field shows through), conductors opaque on top.
function buildGeometryShapes(maxY) {
    return [
        ...dielectricFillShapes(modesSolver, maxY,
            { alpha: 0.12, airAlpha: 0, layer: 'above', lineColor: 'rgba(160,160,160,0.25)' }),
        ...conductorFillShapes(modesSolver, maxY),
    ];
}

// Shared Plotly layout for the Modes tab (field render and geometry-only preview).
function modesPlotLayout(title, view, shapes) {
    return {
        title: { text: title, font: { color: '#fff' } },
        xaxis: darkAxis('Width (mm)', { scaleanchor: 'y', scaleratio: 1, range: view.xRange }),
        yaxis: darkAxis('Height (mm)', { range: view.yRange }),
        margin: { l: 70, r: 90, t: 50, b: 60 },
        showlegend: false, hovermode: 'closest', dragmode: 'pan',
        ...darkBackground(),
        shapes,
    };
}

// Translate the "Solver" dropdown value into the backend + triangular loss options.
// Accepts the legacy 'triangular' value (maps to the accurate MQS full-wave mode).
function getParams() {
    const isCustom = document.getElementById('tl_type').value === 'custom';
    const p = {
        tl_type: document.getElementById('tl_type').value,
        mesh_backend: (document.getElementById('mesh_backend')?.value) ?? 'rectilinear',
        // Custom geometry: the text, plus sidebar parameter values that differ from it
        // (a parameter sweep sets the input without touching the text).
        custom_geom: isCustom ? getCustomGeometryText() : '',
        custom_overrides: isCustom ? getCustomOverrides() : {},
    };
    for (const [key, id, kind] of fieldsFor('p')) {
        if (key === 'sigma') {
            p.sigma = getInputValueUnitless(isCustom ? 'inp_custom_sigma' : 'inp_sigma');
        } else if (key === 'freq_start') {
            p.freq = getInputValue(id);
            p.nx = DEFAULT_GRID_N;
            p.ny = DEFAULT_GRID_N;
        } else if (key === 'estimate_error') {
            // 1/0 like in getUISettings, the same key must not be a boolean in one
            // params object and a number in the other.
            p[key] = readSetting(id, kind);
        } else {
            p[key] = readParam(id, kind);
        }
    }
    return p;
}

// Helper function to add common optional geometry parameters
function updateGeometry() {

    const pbar = document.getElementById('progress_bar');
    pbar.style.width = "0%";

    const p = getParams();
    // The custom geometry editor lists the errors of its own text, so they stay out of
    // the log. An error the editor did not catch (it validates without the solver
    // options) is still logged.
    const isCustom = p.tl_type === 'custom';
    let buildError = null;
    const built = isCustom ? _buildSolverFromParams(p, (msg) => { buildError = msg; })
                           : buildSolverFromParams(p);
    if (buildError && !validateCustomGeometry().errors.length) log(buildError);
    // While the custom geometry text is being edited it is invalid most of the time.
    // Keep the last valid geometry on screen, the editor lists the errors.
    const keep = !built && isCustom && solver && solver.geometry_params;
    if (!keep) solver = built;
    if (plotFieldsSolver && solver !== plotFieldsSolver) releasePlotFields(plotFieldsJob);
}

// Let the worker drop the solve kept for plots at other frequencies.
function releasePlotFields(job) {
    if (job === plotFieldsJob) { plotFieldsSolver = null; plotFieldsJob = null; }
    workerJob('plotRelease', { job }).catch(() => {});
}

// Writes the current geometry (every rectangle the native solver builds, solder mask,
// vias, enclosure walls and all) as custom geometry text and switches to the custom type.
function convertToCustomGeometry() {
    const p = getParams();
    const native = buildSolverFromParams(p, log);
    if (!native) return;
    let text;
    try {
        text = solverToGeometryText(native, { units: 'mm', pinWalls: true });
    } catch (e) { log('ERROR: ' + e.message); return; }
    const names = { microstrip: 'Microstrip', diff_microstrip: 'Differential microstrip', stripline: 'Stripline',
        diff_stripline: 'Differential stripline', gcpw: 'GCPW', diff_gcpw: 'Differential GCPW',
        broadside_stripline: 'Broadside coupled stripline', coax: 'Coaxial line' };
    document.getElementById('inp_custom_sigma').value = native.sigma_cond;
    setCustomGeometryText(`# Converted from: ${names[p.tl_type] || p.tl_type}\n` + text);
    const sel = document.getElementById('tl_type');
    sel.value = 'custom';
    sel.dispatchEvent(new Event('change', { bubbles: true }));
    log('Converted the geometry to custom geometry text.');
}

// Build a FieldSolver from the given sidebar params WITHOUT touching any global state.
// updateGeometry() uses it to (re)build the main `solver`; the Modes tab uses it to build
// its OWN independent solver (modesSolver), so solving modes never disturbs the main
// solve's results. Returns the solver, or null if the parameters are invalid.
async function runSimulation() {
    // One solve at a time. The Modes tab already guards itself this way
    // (runModesSolve). Both can be in flight at once: the worker would
    // serialize them behind its job queue while the button already read "Stop"
    // for a job that had not started, the two solves would fight over the
    // single heartbeat element, and Stop, one worker-wide flag, could land on
    // the wrong job.
    if (isSimulating || isSweeping || isSolvingModes) {
        log('Cannot start a solve while another solve is running.');
        return;
    }

    // Check if solver is valid before attempting to run simulation
    if (!solver) {
        log("ERROR: Cannot run simulation - solver initialization failed due to invalid parameters.");
        return;
    }

    const p = getParams();
    let frequencies = getFrequencies();
    // A waveguide is a DC block: the f=0 limit of its equivalent circuit is correct
    // (Z0 -> 0) but the literal evaluation divides by zero, and computeSParamsSingleEnded
    // has a freq===0 branch that would model it as a lossless through. Drop the point.
    if (solver.allow_dc === false && frequencies.includes(0)) {
        frequencies = frequencies.filter(f => f > 0);
        // The medium names its own reason (allow_dc is a generic flag, not a waveguide one).
        log(`Note: DC (0 Hz) is skipped — ${solver.dc_block_reason || 'this line type does not propagate at DC'}.`);
        if (!frequencies.length) {
            log('ERROR: no non-zero frequencies to solve.');
            return;
        }
    }
    // Calculate hash at the start of the solve. UI can be edited during the
    // solve.
    const solvedGeometry = getGeometryHash();
    const solvedFrequency = getFrequencyHash();
    const solvedKey = solveInputKey();
    // The view model these fields belong to. A mid-solve edit runs
    // updateGeometry(), which replaces `solver` with one describing a different
    // cross-section. Identity is the exact test: updateGeometry() always
    // assigns a fresh object.
    const solvedSolver = solver;
    clearStoredScales();

    const btn = document.getElementById('btn_solve');
    const pbar = document.getElementById('progress_bar');
    const ptext = document.getElementById('progress_text');

    // Monotonic progress bar: the interpolating sweep's completion is only
    // estimated, so different estimators can disagree frame-to-frame. Clamp the
    // displayed width so it only ever increases within a run.
    let displayedProgress = 0;
    const setProgress = (frac) => {
        frac = Math.max(0, Math.min(1, frac));
        if (frac > displayedProgress) displayedProgress = frac;
        pbar.style.width = (displayedProgress * 100) + '%';
    };

    // Change button to "Stop" mode
    btn.textContent = 'Stop';
    btn.classList.add('stop-mode');
    logSolveStarted();
    isSimulating = true;
    // A simulate job replaces the worker's kept solve.
    plotFieldsSolver = null; plotFieldsJob = null;
    updateResultNotices();
    displayedProgress = 0;
    pbar.style.width = '0%';
    if (ptext) ptext.textContent = '';
    heartbeatStart(ptext);
    log("Starting simulation...");
    // The open-wall clearance is measured from the solved field (result warnings).
    for (const w of solver.openBoundaryWarnings({ clearance: false })) log(`⚠ Warning: ${w}`);
    if (solver.mode_type === 'waveguide') {
        // State the single-mode limitation and the usable band up front, every solve.
        log(`Rectangular waveguide: fundamental mode only ` +
            `(cutoff ${(solver.fc / 1e9).toFixed(3)} GHz, single-mode up to ` +
            `${(solver.fc2 / 1e9).toFixed(3)} GHz).`);
        for (const w of solver.waveguideWarnings(Math.max(...frequencies))) log(`⚠ Warning: ${w}`);
    }

    // How the run ended, so the finally block can report it honestly instead of
    // announcing "Done" over a failure or a cancellation.
    let outcome = 'error';

    try {
        // The whole solve pipeline (mesh refinement, causal recompute, interpolating or
        // discrete sweep) runs in solve_worker.js. Parameter validation lives there too,
        // so an invalid combination fails the job and lands in the catch below exactly as
        // it used to. Here we only translate streamed messages into DOM updates.
        const seenModeWarnings = new Set();
        const logModeWarnings = (warnings) => {
            for (const mw of warnings || []) {
                const key = `${mw.type || 'ambiguous'}|${mw.reason || ''}|${mw.mode}`;
                if (seenModeWarnings.has(key)) continue;
                seenModeWarnings.add(key);
                log(`\u26a0 Mode warning: ${mw.message}`);
            }
        };

        const simJob = workerJob('simulate', {
            params: p,
            frequencies,
            opts: {
                useInterpolation: !!document.getElementById('chk_interp_sweep')?.checked,
                interpTolerance: interpTolerance(),
                plotFreq: plotFrequency(),
            },
        }, {
            progress: (m) => {
                setProgress(m.frac);
                if (m.text) heartbeatPhase(m.text);
            },
            // Mesh refinement finished: surface the certificate and mesh-quality notes at
            // the same point in the log as before the sweep output starts.
            meshDone: (m) => {
                if (m.meta.certification && m.meta.certification.pass) {
                    log(`Estimated remaining error ${fmtErrPct(m.meta.certification.err)}% ` +
                        `(tolerance ${(100 * p.tolerance).toFixed(2)}%)`);
                }
                if (m.meta.meshQuality && m.meta.meshQuality.maxQ > 100) {
                    log(`\u26a0 Mesh quality warning: worst triangle Q=${m.meta.meshQuality.maxQ.toFixed(0)} ` +
                        `results may be inaccurate. Try adjusting geometry or enclosure size.`);
                }
                logModeWarnings(m.warnings);
            },
            warnings: (m) => logModeWarnings(m.warnings),
            // Live plots while the sweep runs. The worker sends fields only when they
            // are already the final ones (quasi-static). They carry no surface current
            // yet, so a current view keeps the plot it has until the final fields.
            partial: (m) => {
                if (m.fields && solver === solvedSolver && fieldsHaveWantedView(m.fields)) {
                    applyFields(solver, m.fields);
                    draw();
                }
                if (m.sweepResults) {
                    frequencySweepResults = reviveSweep(m.sweepResults);
                    drawResultsPlot();
                    drawSParamPlot();
                }
            },
        });
        const simJobId = _workerJobId;
        _stoppableJob = simJobId;
        const out = await simJob;

        if (solver === solvedSolver) applyFields(solver, out.fields);
        // The worker keeps this solve for plots at other frequencies while its geometry
        // is on screen.
        if (!out.stopped) {
            if (solver === solvedSolver) {
                plotFieldsSolver = solvedSolver;
                plotFieldsJob = simJobId;
                plotFieldsMaxFreq = Math.max(...frequencies);
                plotFieldsSolveKey = solvedKey;
            } else {
                releasePlotFields(simJobId);
            }
        }
        logModeWarnings(out.sweepWarnings);

        if (out.stopped) {
            outcome = 'stopped';
            log("Simulation stopped by user");
            return;
        }

        frequencySweepResults = reviveSweep(out.sweepResults);
        const results = reviveResult(out.meshResult);

        draw();
        drawResultsPlot();
        drawSParamPlot();

        // Sort results by frequency
        frequencySweepResults.sort((a, b) => a.freq - b.freq);

        // Display summary. Zc = sqrt((R + jwL) / (G + jwC)) and eps_eff = (beta/k0)^2 at the
        // first and last sweep point, Z0 = 1/(c sqrt(C C0)) once from the mesh solve. The
        // first point skips DC (Zc infinite without a conducting dielectric) and, for a
        // waveguide, points below cutoff. A DC-only solve reports its DC row.
        const rows = frequencySweepResults;
        const propagates = r => r.freq > 0 && !Number.isNaN(r.result.modes[0].Z0);
        const dcOnly = rows.length > 0 && rows.every(r => r.freq === 0) && !Number.isNaN(rows[0].result.modes[0].Z0);
        const firstRow = rows.find(propagates) ?? (dcOnly ? rows[0] : undefined);
        const fmtZc = (r, Zc, threshold) => Number.isFinite(Zc.re) ? formatZc(Zc, threshold) : '∞';
        const lastRow = rows[rows.length - 1];
        const ends = !firstRow ? [] : (rows.length === 1 || firstRow === lastRow) ? [firstRow] : [firstRow, lastRow];
        const unit = freqUnit(Math.max(...frequencies));
        const fAt = f => `${formatFreq(f, unit)} ${unit.name}`;
        const freqLine = `Frequency: ${ends.map(r => formatFreq(r.freq, unit)).join(' - ')} ${unit.name}\n`;
        const fRange = ends.length ? ends.map(r => formatFreq(r.freq, unit)).join('-') + ` ${unit.name}` : '';
        const span = (fmt) => ends.map(fmt).join(' - ');
        const loss = r => r.result.modes[0].alpha_total.toFixed(3);
        const lossShort = `${span(loss)} dB/m @ ${fRange}`;
        const lossStr = `Loss: ${span(loss)} dB/m`;
        let cutoffNote = '';
        if (firstRow && rows[0].freq > 0 && firstRow !== rows[0])
            cutoffNote = `\n(from ${fAt(firstRow.freq)}: the sweep starts below cutoff)`;

        // Check if differential results
        if (results.modes.length === 2) {
            const odd = results.modes.find(m => m.mode === 'odd');
            const even = results.modes.find(m => m.mode === 'even');
            const modeOf = (r, name) => r.result.modes.find(m => m.mode === name);
            const zdiff = r => fmtZc(r, modeOf(r, 'odd').Zc.mul(2));
            const zcm = r => fmtZc(r, modeOf(r, 'even').Zc.mul(0.5));
            const eps = (r, name) => modeOf(r, name).eps_eff.toFixed(3);
            const zcFull = (name, scale) => `${span(r => fmtZc(r, modeOf(r, name).Zc.mul(scale), 0))} Ohm`;
            // For an asymmetric pair the two traces are not interchangeable: the physical self
            // terms differ (C11 ≠ C22) and S22 ≠ S11. Surface that here — otherwise the summary
            // looks identical to a symmetric line. odd/even are then only the approximate eigenmodes.
            const mC = results.RLGC_matrix?.C, mL = results.RLGC_matrix?.L;
            const asymStr = (results.physMatrix && mC && mL)
                ? `\n\nAsymmetric pair (S22 ≠ S11; odd/even are the approximate eigenmodes):\n` +
                  `  Self-C:  C11 = ${(mC[0][0] * 1e12).toFixed(2)} pF/m,  C22 = ${(mC[1][1] * 1e12).toFixed(2)} pF/m\n` +
                  `  Self-L:  L11 = ${(mL[0][0] * 1e9).toFixed(2)} nH/m,  L22 = ${(mL[1][1] * 1e9).toFixed(2)} nH/m`
                : '';
            // Traces of different metal or finish: per-line R and L.
            const mR = results.RLGC_matrix?.R;
            const lineStr = (odd.RLGC.dR !== undefined && mR && mL)
                ? `\n\nUnequal traces (line 1 = positive trace, mode conversion included in the S-parameters):\n` +
                  `  R11 = ${mR[0][0].toFixed(2)} Ohm/m,  R22 = ${mR[1][1].toFixed(2)} Ohm/m\n` +
                  `  L11 = ${(mL[0][0] * 1e9).toFixed(2)} nH/m,  L22 = ${(mL[1][1] * 1e9).toFixed(2)} nH/m`
                : '';
            log(`\nDIFFERENTIAL RESULTS:\n` +
                     `======================\n` +
                     freqLine +
                     `Differential Impedance Z_diff: ${zcFull('odd', 2)}  (2 x Zc_odd)\n` +
                     `Common-Mode Impedance Z_common: ${zcFull('even', 0.5)}  (Zc_even / 2)\n` +
                     `\nModal Impedances:\n` +
                     `  Odd-Mode  Z0_odd:  ${odd.Z0.toFixed(2)} Ohm\n` +
                     `  Even-Mode Z0_even: ${even.Z0.toFixed(2)} Ohm\n` +
                     `  Odd-Mode  Zc_odd:  ${zcFull('odd', 1)}\n` +
                     `  Even-Mode Zc_even: ${zcFull('even', 1)}\n` +
                     `  eps_eff odd:  ${span(r => eps(r, 'odd'))}\n` +
                     `  eps_eff even: ${span(r => eps(r, 'even'))}` +
                     `${asymStr}${lineStr}\n` +
                     `\n${lossStr}`,
                `Zdiff ${span(zdiff)} Ω,   Zcm ${span(zcm)} Ω,  ` +
                `εeff odd ${span(r => eps(r, 'odd'))} / even ${span(r => eps(r, 'even'))},  loss ${lossShort}`);
        } else {
            const zc = r => fmtZc(r, r.result.modes[0].Zc);
            const eps = r => r.result.modes[0].eps_eff.toFixed(3);
            log(`\nRESULTS:\n` +
                     `----------------------\n` +
                     (!firstRow
                        ? `Below cutoff across the whole sweep — attenuation only.\n`
                        : freqLine +
                          `Z0: ${results.modes[0].Z0.toFixed(2)} Ohm\n` +
                          `Zc: ${span(r => fmtZc(r, r.result.modes[0].Zc, 0))} Ohm\n` +
                          `eps_eff: ${span(eps)}${cutoffNote}\n`) +
                     `${lossStr}`,
                !firstRow ? 'Below cutoff'
                    : `Zc ${span(zc)} Ω,   εeff ${span(eps)},  loss ${lossShort}`);
        }

        // Update plots
        drawResultsPlot();
        drawSParamPlot();

        // Save geometry and frequency hash for change tracking. Both were sampled when
        // the solve started, so a mid-solve edit correctly reads as a change.
        lastSolvedGeometry = solvedGeometry;
        lastSolvedFrequency = solvedFrequency;
        updateResultNotices();
        outcome = 'done';

    } catch (e) {
        console.error(e);
        log("Error: " + e.message);
    } finally {
        // Restore button to "Solve" mode. The custom geometry editor keeps Solve off while
        // its text has errors, which it may have gained during the solve.
        btn.textContent = 'Solve';
        btn.classList.remove('stop-mode');
        btn.disabled = btn.dataset.customInvalid === '1';
        // A full bar and a "Done" only for a run that actually produced results; a
        // cancelled or failed solve rewinds the bar instead of claiming completion.
        pbar.style.width = outcome === 'done' ? '100%' : '0%';
        // Keep the elapsed total on screen rather than blanking it.
        const elapsed = _hbFormat(Date.now() - _hbStart);
        heartbeatStop(outcome === 'done' ? `Done in ${elapsed}`
                    : outcome === 'stopped' ? `Stopped after ${elapsed}` : '');
        isSimulating = false;
        // A Plot Frequency edited during the run.
        updatePlotFields(true);
    }
}

function getSweepDiffMode() {
    const cb = document.getElementById('sweep-diff');
    return cb ? cb.checked : false;
}

function updateSweepDiffCheckbox() {
    const cb = document.getElementById('sweep-diff');
    if (!cb) return;
    const hasResults = parameterSweepResults && parameterSweepResults.length > 0;
    const isDiff = hasResults && parameterSweepResults[0].result.modes.length === 2;
    cb.disabled = !isDiff;
}

function redrawSweepPlot() {
    if (!parameterSweepResults || parameterSweepResults.length === 0 || !lastSweepParam) return;
    const ySel = document.getElementById('sweep-y-selector').value;
    const cfg = sweepParamConfig(lastSweepParam);
    const xLabel = cfg ? cfg.label + (lastSweepDisplayUnit ? ` (${lastSweepDisplayUnit})` : '') : lastSweepParam;
    drawParameterSweepPlot(parameterSweepResults, xLabel, ySel, getSweepDiffMode());
}

function getGeometryHashExcluding(paramKey) {
    const hash = JSON.parse(getGeometryHash());
    delete hash[paramKey];
    return JSON.stringify(hash);
}

function updateSweepNotice() {
    const notice = document.getElementById('sweep-notice');
    const noticeText = document.getElementById('sweep-notice-text');
    if (!notice) return;

    if (!parameterSweepResults || parameterSweepResults.length === 0) {
        notice.style.display = 'none';
        return;
    }

    if (!isSweeping && lastSweepGeometry) {
        const currentHash = getGeometryHashExcluding(lastSweepParam);
        if (currentHash !== lastSweepGeometry) {
            noticeText.textContent = 'Geometry changed. Run sweep to update results.';
            notice.style.display = 'block';
            return;
        }
    }
    notice.style.display = 'none';
}

function updateSweepParamList() {
    const tlType = document.getElementById('tl_type').value;
    const isDiff = tlType.startsWith('diff_');
    const isGcpw = tlType.includes('gcpw');
    const isStripline = tlType.includes('stripline');
    const isCoax      = tlType === 'coax';
    const isWaveguide = tlType === 'rect_waveguide';
    const isCustom    = tlType === 'custom';
    // All three replace the microstrip stackup entirely with their own geometry block.
    const isSelfBounded = isCoax || isWaveguide || isCustom;
    const useSm       = document.getElementById('chk_solder_mask').checked;
    const useTopDiel  = document.getElementById('chk_top_diel').checked;
    const useGndCut   = document.getElementById('chk_gnd_cut').checked;
    const useEnclosure= document.getElementById('chk_enclosure').checked;
    const usePlating  = document.getElementById('chk_plating').checked;

    const groupEnabled = {
        shared: true,                    // surface roughness, applies to every type
        // Coax and waveguide each have their own geometry block and none of the microstrip
        // stackup inputs, so the 'always' set and the board-stackup groups are off for
        // them. Plating still applies (all-around: the centre conductor / the walls).
        always: !isSelfBounded,
        diff: isDiff,
        gcpw: isGcpw,
        stripline: isStripline,
        coax: isCoax,
        waveguide: isWaveguide,
        custom: isCustom,
        sm: useSm && !isSelfBounded,
        top_diel: useTopDiel && !isSelfBounded,
        gnd_cut: useGndCut && !isSelfBounded,
        enclosure: useEnclosure && !isSelfBounded,
        plating: usePlating && !isCustom,
    };

    const sel = document.getElementById('sweep-x-selector');
    const previousValue = sel.value;
    sel.innerHTML = '';
    for (const [key, cfg] of Object.entries(SWEEP_PARAM_CONFIG)) {
        if (!groupEnabled[cfg.group]) continue;
        const opt = document.createElement('option');
        opt.value = key;
        opt.textContent = cfg.label;
        sel.appendChild(opt);
    }
    // Geometry parameters of the custom type, first in the list.
    if (isCustom) {
        const first = sel.firstChild;
        for (const cfg of customSweepParams()) {
            const opt = document.createElement('option');
            opt.value = cfg.key;
            opt.textContent = cfg.label;
            sel.insertBefore(opt, first);
        }
        if (sel.options.length) sel.selectedIndex = 0;
    }
    // Restore previous selection if still available
    if ([...sel.options].some(o => o.value === previousValue)) sel.value = previousValue;
}

// The sweep tab lists the parameters of the custom geometry text. Activating the type
// fills the text without a geometry change, so it syncs here too.
let customParamKeys = null;
function syncCustomSweepParams() {
    const keys = customSweepParams().map(c => c.key).join();
    if (keys === customParamKeys) return;
    customParamKeys = keys;
    updateSweepParamList();
    autoFillSweepRange();
}

function getZeroDefaultMax(displayUnit) {
    // Return a sensible max in display units for zero-valued params
    const maxInMeters = 2e-6; // 2 μm as reference
    if (!displayUnit) return 1; // unitless
    return +(convertToDisplayUnit(maxInMeters, displayUnit)).toPrecision(4);
}

function autoFillSweepRange() {
    const xSel = document.getElementById('sweep-x-selector').value;
    const cfg = sweepParamConfig(xSel);
    if (!cfg) return;
    const inputEl = document.getElementById(cfg.inputId);
    if (!inputEl) return;
    const displayUnit = getSweepDisplayUnit(cfg);
    const isUnitless = !displayUnit || cfg.fixedUnit;
    const currentVal = isUnitless ? parseFloat(inputEl.value) : getInputValue(cfg.inputId);
    if (isNaN(currentVal) || currentVal < 0) return;
    const minInput = document.getElementById('sweep-x-min');
    const maxInput = document.getElementById('sweep-x-max');
    let minNum, maxNum;
    if (currentVal === 0) {
        // For zero-valued params (e.g. roughness), use a sensible default range
        minNum = 0;
        maxNum = getZeroDefaultMax(isUnitless ? '' : displayUnit);
    } else {
        const displayVal = isUnitless ? currentVal : convertToDisplayUnit(currentVal, displayUnit);
        minNum = +(displayVal * 0.5).toPrecision(4);
        maxNum = +(displayVal * 2.0).toPrecision(4);
    }
    if (isUnitless) {
        minInput.value = minNum;
        maxInput.value = maxNum;
    } else {
        minInput.value = `${minNum} ${displayUnit}`;
        maxInput.value = `${maxNum} ${displayUnit}`;
    }
}

/**
 * Extract the unit suffix the user typed into a geometry input field.
 * Falls back to getDefaultUnit() when the field has no explicit unit.
 */
function extractUnitFromInput(inputId) {
    const el = document.getElementById(inputId);
    if (!el) return '';
    const raw = (el.value || '').trim();
    const match = raw.match(/[+-]?(?:\d+\.?\d*|\.\d+)(?:[e][+-]?\d+)?\s*([a-zμµ]+)$/i);
    if (match && match[1]) return match[1];
    return window.getDefaultUnit ? window.getDefaultUnit(inputId) : '';
}

/**
 * Determine the display unit for a sweep parameter.
 * fixedUnit params (sigma) use their fixed label; others derive from geometry input.
 * Returns '' for unitless params (er, tand).
 */
function getSweepDisplayUnit(cfg) {
    if (cfg.fixedUnit) return cfg.fixedUnit;
    const unit = extractUnitFromInput(cfg.inputId);
    return unit || '';
}

function convertToDisplayUnit(valueSI, unit) {
    const factors = {
        // 'μm' is U+03BC, the second 'µm' is U+00B5.
        'mm': 1e3, 'μm': 1e6, 'µm': 1e6, 'um': 1e6, 'nm': 1e9,
        'cm': 1e2, 'm': 1,
        'mil': 1 / 25.4e-6, 'mils': 1 / 25.4e-6,
        'in': 1 / 25.4e-3, 'inch': 1 / 25.4e-3, 'inches': 1 / 25.4e-3,
        'ft': 1 / 0.3048, 'foot': 1 / 0.3048, 'feet': 1 / 0.3048,
        'GHz': 1e-9, 'MHz': 1e-6,
        'S/m': 1,
    };
    return valueSI * (factors[unit] || 1);
}

async function runParameterSweep() {
    const xSel = document.getElementById('sweep-x-selector').value;
    const ySel = document.getElementById('sweep-y-selector').value;
    const cfg = sweepParamConfig(xSel);
    const displayUnit = getSweepDisplayUnit(cfg);
    const isUnitless = !displayUnit || cfg.fixedUnit;
    const parseVal = (str) => {
        if (isUnitless) return parseFloat(str);
        return window.parseValueWithUnit ? window.parseValueWithUnit(str, displayUnit) : parseFloat(str);
    };
    const minSI = parseVal(document.getElementById('sweep-x-min').value);
    const maxSI = parseVal(document.getElementById('sweep-x-max').value);
    // For unitless/fixedUnit params, min/maxSI are already display values.
    // Unit conversion leaves float noise trim to 12 significant digits so the
    // output is clean.
    const trimFloat = (v) => parseFloat(v.toPrecision(12));
    const minDisplay = trimFloat(isUnitless ? minSI : convertToDisplayUnit(minSI, displayUnit));
    const maxDisplay = trimFloat(isUnitless ? maxSI : convertToDisplayUnit(maxSI, displayUnit));
    const nPoints = parseInt(document.getElementById('sweep-points').value, 10);
    const freqHz = getInputValue('sweep-freq');

    if (isNaN(minDisplay) || isNaN(maxDisplay) || minDisplay >= maxDisplay) { log('ERROR: Invalid sweep range.'); return; }
    if (isNaN(nPoints) || nPoints < 2)                       { log('ERROR: Points must be >= 2.'); return; }
    if (isNaN(freqHz) || freqHz < 0)                         { log('ERROR: Invalid frequency.'); return; }

    const runBtn = document.getElementById('btn-run-sweep');
    const stopBtn = document.getElementById('btn-stop-sweep');
    const solveBtn = document.getElementById('btn_solve');
    const progressText = document.getElementById('sweep-progress-text');
    runBtn.style.display = 'none';
    stopBtn.style.display = '';
    solveBtn.disabled = true;
    isSweeping = true;
    parameterSweepResults = [];

    const inputEl = document.getElementById(cfg.inputId);
    const originalValue = inputEl.value;
    const p = getParams();

    // Save geometry hash excluding the swept parameter
    lastSweepParam = xSel;
    lastSweepDisplayUnit = displayUnit;
    lastSweepGeometry = getGeometryHashExcluding(xSel);
    updateSweepNotice();

    heartbeatStart(progressText);
    log(`Parameter sweep: ${cfg.label} ${minDisplay}–${maxDisplay}${displayUnit ? ' ' + displayUnit : ''} (${nPoints} pts) @ ${(freqHz/1e9).toFixed(3)} GHz`);

    try {
        // Every sweep point's params are resolved HERE, on the thread that owns the
        // sidebar inputs: the worker has no DOM, so it cannot read the swept input
        // itself. It receives a ready-made list and never touches the form.
        const points = [];
        for (let i = 0; i < nPoints; i++) {
            // Interpolation reintroduces float noise (0.105 + 0.315/9 = 0.14000000000000001).
            const displayVal = trimFloat(minDisplay + (maxDisplay - minDisplay) * i / (nPoints - 1));
            inputEl.value = isUnitless ? displayVal : `${displayVal} ${displayUnit}`;
            points.push({ displayVal, params: getParams() });
        }
        inputEl.value = originalValue;
        updateGeometry();

        const sweepJob = workerJob('paramSweep', {
            points, freqHz, opts: { estimateError: !!p.estimate_error },
        }, {
            // Surface the first point's verification outcome once. Later points run the
            // legacy per-pass gate with the same refinement settings.
            firstPointCert: (m) => {
                if (m.certification && m.certification.pass) {
                    log(`Estimated remaining error (first-point) ${fmtErrPct(m.certification.err)}%. ` +
                        `Later sweep points are not re-verified.`);
                }
                // Every accuracy note, not just the first: one solve can carry a failed
                // certificate and a loss-accuracy note (skin-transition,
                // broadside-proximity) at once.
                for (const aw of m.warnings || []) {
                    if (aw.type !== 'accuracy' && aw.type !== 'open-boundary') continue;
                    log(`\u26a0 First point: ${aw.message} Later sweep points are not re-verified.`);
                }
            },
            partial: (m) => {
                parameterSweepResults = m.sweepResults.map(r => ({ ...r, result: reviveResult(r.result) }));
                heartbeatPhase(`${m.index + 1}/${m.total}`);
                if (m.index === 0) updateSweepDiffCheckbox();
                redrawSweepPlot();
            },
        });
        _stoppableJob = _workerJobId;
        const out = await sweepJob;

        parameterSweepResults = out.results.map(r => ({ ...r, result: reviveResult(r.result) }));
        redrawSweepPlot();
        log(`Sweep complete: ${parameterSweepResults.length} points.`);
    } catch(e) {
        console.error(e);
        log('Sweep error: ' + e.message);
    } finally {
        inputEl.value = originalValue;
        updateGeometry();
        runBtn.style.display = '';
        stopBtn.style.display = 'none';
        // The custom geometry editor keeps Solve and Run Sweep off while its text has errors.
        solveBtn.disabled = runBtn.disabled = solveBtn.dataset.customInvalid === '1';
        isSweeping = false;
        heartbeatStop('');
    }
}

function resizeCanvas() {
    const Plotly = getPlotly();
    if (!Plotly) return;
    for (const id of ['sim_canvas', 'modes-plot']) {
        const container = document.getElementById(id);
        if (container && container.data) Plotly.Plots.resize(container);
    }
}

function bindEvents() {
    document.getElementById('btn_solve').onclick = () => {
        const btn = document.getElementById('btn_solve');
        if (btn.textContent === 'Stop') {
            // Cancellation lives entirely in the worker, which polls the flag between
            // WASM calls — the same granularity the main-thread loop had.
            workerStop();
            log("Stop requested...");
        } else {
            flushCustomGeometryEdits();
            // Before updateGeometry(), which replaces the solver the kept solve belongs to.
            if (replotInsteadOfSolve()) return;
            // Start the simulation
            updateGeometry(); // Ensure geometry is updated with latest parameters
            runSimulation();
        }
    };

    // Tab switching
    document.querySelectorAll('.tab-button').forEach(btn => {
        btn.addEventListener('click', () => {
            switchTab(btn.dataset.tab);
        });
    });

    // Modes tab events
    const btnSolveModes = document.getElementById('btn-solve-modes');
    if (btnSolveModes) btnSolveModes.addEventListener('click', runModesSolve);

    // Parameter sweep events
    document.getElementById('btn-run-sweep').addEventListener('click', () => {
        if (isSweeping || isSimulating || isSolvingModes) {
            log('Cannot sweep while another solve is running.');
            return;
        }
        flushCustomGeometryEdits();
        runParameterSweep();
    });
    document.getElementById('btn-stop-sweep').addEventListener('click', () => {
        workerStop();
        log('Sweep stop requested...');
    });

    document.getElementById('sweep-diff').addEventListener('change', redrawSweepPlot);
    document.getElementById('sweep-y-selector').addEventListener('change', redrawSweepPlot);

    document.getElementById('sweep-x-selector').addEventListener('change', autoFillSweepRange);

    document.getElementById('tl_type').addEventListener('change', updateSweepParamList);
    ['chk_solder_mask','chk_top_diel','chk_gnd_cut','chk_enclosure','chk_plating']
        .forEach(id => document.getElementById(id).addEventListener('change', updateSweepParamList));

    updateSweepParamList();
    autoFillSweepRange();

    // Results and S-parameter plot controls redraw their plot once there are results.
    for (const [id, ev, redraw] of [
        ['results-plot-selector', 'change', drawResultsPlot],
        ['sparam-length', 'input', drawSParamPlot],
        ['sparam-z-ref', 'input', drawSParamPlot],
        ['sparam-plot-mode', 'change', drawSParamPlot],
        ['sparam-diff', 'change', drawSParamPlot],
        ['results-log-x', 'change', drawResultsPlot],
        ['results-diff', 'change', drawResultsPlot],
        ['sparam-log-x', 'change', drawSParamPlot],
    ]) {
        document.getElementById(id)?.addEventListener(ev, () => {
            if (frequencySweepResults) redraw();
        });
    }

    // Freeze buttons (linked — both tabs share the same frozen state)
    const freezeResultsBtn = document.getElementById('freeze-results-btn');
    const freezeSParamsBtn = document.getElementById('freeze-sparams-btn');
    const freezeBtns = [freezeResultsBtn, freezeSParamsBtn].filter(Boolean);

    function toggleFreeze() {
        if (isFrozen()) {
            unfreeze();
            for (const btn of freezeBtns) {
                btn.textContent = 'Freeze';
                btn.classList.remove('freeze-active');
            }
        } else {
            if (!frequencySweepResults || frequencySweepResults.length === 0) return;
            freeze();
            for (const btn of freezeBtns) {
                btn.textContent = 'Unfreeze';
                btn.classList.add('freeze-active');
            }
        }
        if (frequencySweepResults) {
            drawResultsPlot();
            drawSParamPlot();
        }
    }

    for (const btn of freezeBtns) {
        btn.addEventListener('click', toggleFreeze);
    }

    // Touchstone header lines of a custom geometry: the defaults for conductors without
    // their own metal, then the geometry text as solved, with the sidebar parameter
    // values that differ from the text written into it.
    const customGeometryLines = p => {
        let text = p.custom_geom;
        for (const [name, v] of Object.entries(p.custom_overrides || {})) {
            text = setParamInText(text, name, String(+v.toPrecision(12)));
        }
        return [
            `!   Default conductivity: ${p.sigma.toExponential(2)} S/m`,
            `!   Default surface roughness RMS: ${(p.rq * 1e6).toFixed(2)} um`,
            ...(p.plating_thick_corners ? ['!   Model thick plating: on'] : []),
            '!   Geometry:',
            ...text.replace(/\s+$/, '').split('\n').map(l => (l.trim() ? `!     ${l}` : '!')),
        ];
    };

    // Export SnP button
    const exportSnpBtn = document.getElementById('export-snp');
    if (exportSnpBtn) {
        exportSnpBtn.addEventListener('click', () => {
            if (isSimulating) {
                log('Cannot export while simulation is running.');
                return;
            }
            if (!frequencySweepResults || frequencySweepResults.length === 0) {
                log('No results to export. Run simulation first.');
                return;
            }

            // Check if geometry or frequency has changed
            const currentGeometry = getGeometryHash();
            const currentFrequency = getFrequencyHash();
            if ((lastSolvedGeometry && currentGeometry !== lastSolvedGeometry) ||
                (lastSolvedFrequency && currentFrequency !== lastSolvedFrequency)) {
                log('Cannot export: Geometry or frequency has changed. Run simulation again.');
                return;
            }
            const length = getInputValue('sparam-length');
            const Z_ref = getInputValueUnitless('sparam-z-ref');
            const isDifferential = solver && solver.is_differential;
            const p = getParams();
            const params = {
                tlType: p.tl_type,
                traceWidth: p.w,
                traceThickness: p.t,
                substrateHeight: p.h,
                epsilonR: p.er,
                tanDelta: p.tand,
                sigma: p.sigma,
                traceSpacing: p.tl_type.startsWith('diff_') ? p.trace_spacing : null,
                surfaceRoughness: p.tl_type === 'custom' ? null : p.rq,
                // Coax, waveguide and custom geometry are not trace-on-substrate
                // stackups, so they describe themselves.
                geometryLines: p.tl_type === 'custom' ? customGeometryLines(p) : p.tl_type === 'coax' ? [
                    `!   Inner conductor diameter: ${(p.coax_d * 1e6).toFixed(1)} um`,
                    `!   Dielectric diameter (shield ID): ${(p.coax_D * 1e6).toFixed(1)} um`,
                    `!   Dielectric permittivity: ${p.coax_er}`,
                    `!   Loss tangent: ${p.coax_tand}`,
                    `!   Conductivity: ${p.coax_sigma.toExponential(2)} S/m`,
                ] : p.tl_type === 'rect_waveguide' ? [
                    `!   Broad wall a: ${(p.wg_a * 1e6).toFixed(1)} um`,
                    `!   Narrow wall b: ${(p.wg_b * 1e6).toFixed(1)} um`,
                    `!   Fill permittivity: ${p.wg_er}`,
                    `!   Loss tangent: ${p.wg_tand}`,
                    `!   Wall conductivity: ${p.wg_sigma.toExponential(2)} S/m`,
                    `!   Fundamental mode only (TE10 when a >= b)`,
                ] : null,
                // Mirror what the solver actually built, not the raw checkboxes: a coax
                // selects conductors and a waveguide plates its whole wall, so the
                // top/sides/bottom boxes (hidden for both) would misdescribe the run.
                // Custom geometry carries its plating per conductor in the text.
                plating: p.tl_type === 'custom' ? null : platingOptions(p,
                    p.tl_type === 'coax' ? { inner: p.coax_plating_inner, outer: p.coax_plating_outer }
                    : p.tl_type === 'rect_waveguide' ? { all: true }
                    : null),
                freqStart: frequencySweepResults[0].freq,
                freqStop: frequencySweepResults[frequencySweepResults.length - 1].freq,
                numPoints: frequencySweepResults.length
            };
            // generateS2P rejects a sweep it cannot represent (e.g. every point below a
            // waveguide's cutoff) with a message written for the user, surface it in the
            // log like every other failure here, rather than letting it escape to the
            // console where the button just appears to do nothing.
            try {
                log(`Exported ${exportSnP(frequencySweepResults, length, Z_ref, isDifferential, params)}`);
            } catch (e) {
                log(`Export failed: ${e && e.message ? e.message : e}`);
            }
        });
    }

    // Export CSV button
    const exportCsvBtn = document.getElementById('export-csv-btn');
    if (exportCsvBtn) {
        exportCsvBtn.addEventListener('click', () => {
            if (isSimulating) {
                log('Cannot export while simulation is running.');
                return;
            }
            if (!frequencySweepResults || frequencySweepResults.length === 0) {
                log('No results to export. Run simulation first.');
                return;
            }
            const isDifferential = frequencySweepResults[0].result.modes.length === 2;
            // A pair gets the odd mode's columns, then the even mode's.
            const suffixes = isDifferential ? ['_odd', '_even'] : [''];
            const header = ['Freq_Hz'];
            for (const m of suffixes) {
                header.push(`Re_Z0${m}_Ohm`, `Im_Z0${m}_Ohm`, `eps_eff${m}`,
                    `conductor_loss${m}_dBpm`, `dielectric_loss${m}_dBpm`, `total_loss${m}_dBpm`,
                    `R${m}_Ohmpm`, `L${m}_Hpm`, `G${m}_Spm`, `C${m}_Fpm`);
            }
            const rows = [header];
            for (const { freq, result } of frequencySweepResults) {
                const row = [freq];
                for (const m of result.modes.slice(0, suffixes.length)) {
                    row.push(m.Zc.re, m.Zc.im, m.eps_eff,
                        m.alpha_c, m.alpha_d, m.alpha_total,
                        m.RLGC.R, m.RLGC.L, m.RLGC.G, m.RLGC.C);
                }
                rows.push(row);
            }
            downloadFile(rows.map(r => r.join(',')).join('\n'), 'results.csv', 'text/csv');
            log('Exported results.csv');
        });
    }

    // Frequency points validation. Default to 1 when empty
    const freqPointsEl = document.getElementById('freq-points');
    if (freqPointsEl) {
        freqPointsEl.addEventListener('blur', () => {
            const val = parseInt(freqPointsEl.value);
            if (isNaN(val) || val < 1 || freqPointsEl.value.trim() === '') {
                freqPointsEl.value = '1';
            }
        });
    }

    // Solver and plot parameter validation
    const validationRules = {
        'freq-start': { default: 0.1, label: 'Start frequency' },
        'freq-stop': { default: 10, label: 'Stop frequency' },
        'inp_max_iters': { min: 1, default: 10, integer: true, label: 'Max iterations' },
        'inp_max_nodes': { min: 1, default: 20, integer: true, label: 'Max nodes' },
        'inp_tolerance': { min: 0, max: 100, default: 1, label: 'Tolerance' },
        'sparam-length': { default: 0.01, label: 'Line length' },
        'sparam-z-ref': { min: 1, default: 50, label: 'Reference impedance' }
    };

    Object.entries(validationRules).forEach(([id, rule]) => {
        const el = document.getElementById(id);
        if (el) {
            el.addEventListener('blur', () => {
                let val = rule.integer ? parseInt(el.value) : parseFloat(el.value);
                if (isNaN(val) || el.value.trim() === '') {
                    el.value = rule.default;
                }
                else if (val < rule.min) {
                    el.value = rule.min;
                }
                else if (val > rule.max) {
                    el.value = rule.max;
                }
            });
        }
    });

    // Rebuild the preview and flag the solved results, sweep and modes as stale.
    const geometryEdited = (resetZoom = false, afterDraw = null) => {
        updateGeometry();
        draw(resetZoom);
        afterDraw?.();
        updateResultNotices();
        updateSweepNotice();
        updateModesNotice();
    };

    // Real-time geometry updates for the parameter inputs and checkboxes
    for (const [, id, kind] of fieldsFor('l')) {
        const ev = kind === 'chk' ? 'change' : 'input';
        document.getElementById(id)?.addEventListener(ev, () => geometryEdited());
    }
    // Hashed settings that leave the preview as it is (causal materials) only flag the
    // results as stale.
    for (const [, id, kind, use] of fieldsFor('h')) {
        if (use.includes('l')) continue;
        const ev = kind === 'chk' ? 'change' : 'input';
        document.getElementById(id)?.addEventListener(ev, () => {
            updateResultNotices();
            updateSweepNotice();
            updateModesNotice();
        });
    }

    // Custom geometry editor: every edit rebuilds the preview, and the parameter list
    // of the sweep tab follows the parameters defined in the text.
    initCustomGeometryEditor({
        log,
        onGeometryChange: (resetZoom = false) => {
            if (document.getElementById('tl_type').value !== 'custom') return;
            geometryEdited(resetZoom === true, syncCustomSweepParams);
        },
        onHighlightChange: () => {
            if (document.getElementById('tl_type').value === 'custom') draw();
        },
    });
    document.getElementById('btn-convert-custom').addEventListener('click', convertToCustomGeometry);

    // Transmission line type selector - reset zoom when type changes
    document.getElementById('tl_type').addEventListener('change', () => {
        if (document.getElementById('tl_type').value === 'custom') {
            activateCustomGeometry();
            syncCustomSweepParams();
        }
        geometryEdited(true);  // Reset zoom/pan for new geometry
    });

    // Frequency inputs - update notices when changed
    ['freq-start', 'freq-stop', 'freq-points'].forEach(id => {
        const el = document.getElementById(id);
        if (el) {
            el.addEventListener('change', () => {
                updateResultNotices();
            });
        }
    });

    // Modes frequency - flag the modes plot as stale when changed after a solve
    const modesFreqEl = document.getElementById('modes-freq');
    if (modesFreqEl) {
        modesFreqEl.addEventListener('input', updateModesNotice);
        modesFreqEl.addEventListener('change', updateModesNotice);
    }

    // Plot options redraw the solved fields.
    const redrawFields = () => { if (solver && solver.solution_valid) draw(); };
    for (const id of ['plot-mode', 'plot-streamlines', 'plot-contours']) {
        document.getElementById(id)?.addEventListener('change', redrawFields);
    }
    const plotFreqEl = document.getElementById('plot-freq');
    if (plotFreqEl) plotFreqEl.addEventListener('change', () => updatePlotFields());
    const plotEfieldDbEl = document.getElementById('plot-efield-db');
    if (plotEfieldDbEl) {
        // dB and linear keep separate scales, so the dialog reloads the new one.
        plotEfieldDbEl.addEventListener('change', () => {
            redrawFields();
            if (scaleDialogOpen) openScaleDialog();
        });
    }

    // Copy link button
    const copyLinkBtn = document.getElementById('copy-link-btn');
    if (copyLinkBtn) {
        copyLinkBtn.addEventListener('click', copySettingsLink);
    }

    // Scale dialog event listeners
    setupScaleDialog();
}

// --- Scale Dialog Management ---

// Store separate scales for each view type
const scaleRanges = {
    potential: { min: null, max: null },
    efield: { min: null, max: null },
    efield_db: { min: null, max: null },
    current: { min: null, max: null },
    current_db: { min: null, max: null },
    density: { min: null, max: null },
    density_db: { min: null, max: null },
    geometry: { min: null, max: null }
};

// A new solve autoscales every view: the field range of the previous geometry means
// nothing for the next one.
function clearStoredScales() {
    for (const r of Object.values(scaleRanges)) { r.min = null; r.max = null; }
    closeScaleDialog();
}

let scaleDialogOpen = false;

function getViewType(view) {
    if (view.startsWith('potential')) return 'potential';
    if (view.startsWith('efield')) return view.endsWith('_db') ? 'efield_db' : 'efield';
    if (view.startsWith('current')) return view.endsWith('_db') ? 'current_db' : 'current';
    if (view.startsWith('density')) return view.endsWith('_db') ? 'density_db' : 'density';
    return 'geometry';
}

function setupScaleDialog() {
    const zMinInput = document.getElementById("zMinInput");
    const zMaxInput = document.getElementById("zMaxInput");
    const zMinSlider = document.getElementById("zMinSlider");
    const zMaxSlider = document.getElementById("zMaxSlider");

    if (zMinInput) {
        zMinInput.addEventListener("input", () => {
            const min = Number(zMinInput.value);
            const max = Number(zMaxInput.value);
            if (zMinSlider) zMinSlider.value = min;
            updateScaleFromDialog();
        });
    }

    if (zMaxInput) {
        zMaxInput.addEventListener("input", () => {
            const min = Number(zMinInput.value);
            const max = Number(zMaxInput.value);
            if (zMaxSlider) zMaxSlider.value = max;
            updateScaleFromDialog();
        });
    }

    if (zMinSlider) {
        zMinSlider.addEventListener("input", (e) => {
            zMinInput.value = Number(e.target.value).toFixed(2);
            updateScaleFromDialog();
        });
    }

    if (zMaxSlider) {
        zMaxSlider.addEventListener("input", (e) => {
            zMaxInput.value = Number(e.target.value).toFixed(2);
            updateScaleFromDialog();
        });
    }
}

function updateScaleFromDialog() {
    const min = Number(document.getElementById("zMinInput").value);
    const max = Number(document.getElementById("zMaxInput").value);

    // Save to current view's scale
    const scaleInfo = getScaleRange();
    const viewType = getViewType(scaleInfo.view);
    scaleRanges[viewType].min = min;
    scaleRanges[viewType].max = max;

    // Apply to plot
    setScaleRange(min, max);
}

function toggleScaleDialog() {
    const dlg = document.getElementById("scaleDialog");
    if (!dlg) return;

    if (scaleDialogOpen) {
        dlg.style.display = "none";
        scaleDialogOpen = false;
    } else {
        openScaleDialog();
    }
}

function openScaleDialog() {
    const dlg = document.getElementById("scaleDialog");
    if (!dlg) return;

    const scaleInfo = getScaleRange();
    const viewType = getViewType(scaleInfo.view);
    const dbRow = document.getElementById("efieldDbRow");
    if (dbRow) dbRow.style.display = /^(efield|current|density)/.test(viewType) ? 'block' : 'none';

    // Get actual data range (before any user scaling)
    const actualDataRange = getActualDataRange();
    let actualMin = actualDataRange.min !== null ? actualDataRange.min : 0;
    let actualMax = actualDataRange.max !== null ? actualDataRange.max : 1;

    // Get stored scale or use current computed scale
    let minVal = scaleRanges[viewType].min;
    let maxVal = scaleRanges[viewType].max;

    // If no stored scale, use current computed values
    if (minVal === null || maxVal === null) {
        minVal = actualMin;
        maxVal = actualMax;
        scaleRanges[viewType].min = minVal;
        scaleRanges[viewType].max = maxVal;
    }

    document.getElementById("zMinInput").value = Number(minVal).toFixed(2);
    document.getElementById("zMaxInput").value = Number(maxVal).toFixed(2);

    const minSlider = document.getElementById("zMinSlider");
    const maxSlider = document.getElementById("zMaxSlider");

    // Determine slider bounds based on view type and actual data
    let sliderMinBound, sliderMaxBound;

    if (viewType === 'potential') {
        // Potential has theoretical bounds: [-1,1] for differential, [0,1] for single-ended
        // Check if differential odd mode by looking at whether actualMin is negative
        const isPotentialOddMode = actualMin < -0.1;
        sliderMinBound = isPotentialOddMode ? -1.0 : 0.0;
        sliderMaxBound = 1.0;
    } else if (viewType.endsWith('_db')) {
        sliderMinBound = actualMin - 40;
        sliderMaxBound = actualMax + 20;
    } else {
        // For E-field and geometry, 1.5x the autoscale range, reaching the true peak
        sliderMinBound = actualMin < -0.1 ? actualMin * 1.5 : 0.0;
        sliderMaxBound = Math.max(actualMax * 1.5, actualDataRange.peak || 0);
    }

    if (minSlider) {
        minSlider.min = sliderMinBound;
        minSlider.max = maxVal;
        minSlider.step = (minSlider.max - minSlider.min) / 200;
        minSlider.value = minVal;
    }

    if (maxSlider) {
        maxSlider.min = minVal;
        maxSlider.max = sliderMaxBound;
        maxSlider.step = (maxSlider.max - maxSlider.min) / 200;
        maxSlider.value = maxVal;
    }

    dlg.style.display = "block";
    scaleDialogOpen = true;
}

function closeScaleDialog() {
    const dlg = document.getElementById("scaleDialog");
    if (dlg) {
        dlg.style.display = "none";
        scaleDialogOpen = false;
    }
}

// Reset color scale to actual data range (called when autoscale is triggered)
function resetColorScale() {
    // Clear stored scale for current view
    const scaleInfo = getScaleRange();
    const viewType = getViewType(scaleInfo.view);
    scaleRanges[viewType].min = null;
    scaleRanges[viewType].max = null;

    // Redraw to apply actual data range
    draw();
}

// Make functions globally accessible for HTML onclick handlers
window.toggleScaleDialog = toggleScaleDialog;
window.closeScaleDialog = closeScaleDialog;
window.resetColorScale = resetColorScale;

// Handle view changes to restore appropriate scale
window.onViewChanged = function(view) {
    // Close scale dialog when switching views
    // The user can reopen it to adjust the scale for the new view
    if (scaleDialogOpen) {
        closeScaleDialog();
    }
};

// Get stored scale override for current view (called by plot.js)
window.getStoredScale = function(view) {
    const viewType = getViewType(view);
    const stored = scaleRanges[viewType];
    if (stored.min !== null && stored.max !== null) {
        return { min: stored.min, max: stored.max };
    }
    return null;
};

function init() {
    // Set up globals for plot.js
    setGlobals({
        getSolver: () => solver,
        getFrequencySweepResults: () => frequencySweepResults,
        getInputValue: getInputValue
    });

    bindEvents();

    // Check for URL parameters and restore settings if present
    const hasURLParams = loadSettingsFromURL();

    // Update checkbox section visibility after settings restore
    if (typeof toggleParameterVisibility === 'function') {
        toggleParameterVisibility();
    }
    // Update checkbox sections
    syncCheckboxSections();

    // Interpolating sweep toggle
    const interpChk = document.getElementById('chk_interp_sweep');
    const interpTolGroup = document.getElementById('interp-tolerance-group');
    if (interpChk && interpTolGroup) {
        const updateInterpVisibility = () => {
            interpTolGroup.style.display = interpChk.checked ? '' : 'none';
        };
        interpChk.addEventListener('change', updateInterpVisibility);
        updateInterpVisibility();
    }

    // A reload keeps the form state: the type may come back as custom, with the
    // geometry text the browser restored or, without one, the first template.
    if (document.getElementById('tl_type').value === 'custom') {
        activateCustomGeometry();
        syncCustomSweepParams();
    }
    updateGeometry();
    draw();
    initLayoutPanels();
    resizeCanvas();
    window.addEventListener('resize', resizeCanvas);
    // Panes shown or hidden without a window resize (custom geometry editor, log)
    // change the plot size too; Plotly.react with an unchanged config does not re-measure.
    if (typeof ResizeObserver !== 'undefined') {
        let queued = false;
        const observer = new ResizeObserver(() => {
            if (queued) return;
            queued = true;
            requestAnimationFrame(() => { queued = false; resizeCanvas(); });
        });
        for (const id of ['sim_canvas', 'modes-plot']) {
            const container = document.getElementById(id);
            if (container) observer.observe(container);
        }
    }
    log("Ready. Click 'Solve' to start simulation.");
    loadSettingsFromFragment();
    window.addEventListener('hashchange', loadSettingsFromFragment);
}

// Start when DOM is ready
window.addEventListener('DOMContentLoaded', init);

// Redraw plots when Plotly finishes loading (in case solver ran before Plotly loaded)
window.addEventListener('plotly-loaded', () => {
    if (solver) {
        draw();
    }
    if (frequencySweepResults && frequencySweepResults.length > 0) {
        drawResultsPlot();
        drawSParamPlot();
    }
});
