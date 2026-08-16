import { Complex } from './complex.js';
import { computeSParamsSingleEnded, computeSParamsDifferential, sParamTodB, usableSweepPoints } from './sparameters.js';
import { exportSnP } from './snp_export.js';
import { draw, drawResultsPlot, drawSParamPlot, drawParameterSweepPlot, setGlobals, setCurrentView, getScaleRange, setScaleRange, getActualDataRange,
    freeze, unfreeze, isFrozen, conductorFillShapes, dielectricFillShapes, computeGeometryView } from './plot.js';
import { buildSolverFromParams as _buildSolverFromParams, platingOptions } from './solver_factory.js';

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
    if (raw === "auto") {
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
    modes_mesh_density: 12,   // bulk cells per wavelength (TriBackend wavelengthDensity)
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
};

/**
 * Get current UI settings as a serializable object (in display units)
 */
function getUISettings() {
    // Helper to get display value (strip unit and return raw number)
    const getDisplayValue = (id) => {
        const element = document.getElementById(id);
        if (!element) return NaN;
        const defaultUnit = window.getDefaultUnit ? window.getDefaultUnit(id) : '';
        const siValue = window.parseValueWithUnit ?
            window.parseValueWithUnit(element.value, defaultUnit) :
            parseFloat(element.value);

        // Convert back to display units for serialization
        const unitMap = {
            'mm': 1e3, 'μm': 1e6, 'GHz': 1e-9, 'm': 1
        };
        const scale = unitMap[defaultUnit] || 1;
        return siValue * scale;
    };

    return {
        tl_type: document.getElementById('tl_type').value,
        mesh_backend: (document.getElementById('mesh_backend')?.value) ?? 'rectilinear',
        w: getDisplayValue('inp_w'),
        h: getDisplayValue('inp_h'),
        t: getDisplayValue('inp_t'),
        er: getInputValueUnitless('inp_er'),
        tand: getInputValueUnitless('inp_tand'),
        sigma: getInputValueUnitless('inp_sigma'),
        freq_start: getDisplayValue('freq-start'),
        freq_stop: getDisplayValue('freq-stop'),
        freq_points: parseInt(document.getElementById('freq-points').value),
        trace_spacing: getDisplayValue('inp_trace_spacing'),
        gap: getDisplayValue('inp_gap'),
        via_gap: getDisplayValue('inp_via_gap'),
        stripline_top_h: getDisplayValue('inp_air_top'),
        er_top: getInputValueUnitless('inp_er_top'),
        tand_top: getInputValueUnitless('inp_tand_top'),
        use_sm: document.getElementById('chk_solder_mask').checked ? 1 : 0,
        sm_t_sub: getDisplayValue('inp_sm_t_sub'),
        sm_t_trace: getDisplayValue('inp_sm_t_trace'),
        sm_t_side: getDisplayValue('inp_sm_t_side'),
        sm_er: getInputValueUnitless('inp_sm_er'),
        sm_tand: getInputValueUnitless('inp_sm_tand'),
        use_top_diel: document.getElementById('chk_top_diel').checked ? 1 : 0,
        top_diel_h: getDisplayValue('inp_top_diel_h'),
        top_diel_er: getInputValueUnitless('inp_top_diel_er'),
        top_diel_tand: getInputValueUnitless('inp_top_diel_tand'),
        use_gnd_cut: document.getElementById('chk_gnd_cut').checked ? 1 : 0,
        gnd_cut_w: getDisplayValue('inp_gnd_cut_w'),
        gnd_cut_h: getDisplayValue('inp_gnd_cut_h'),
        use_enclosure: document.getElementById('chk_enclosure').checked ? 1 : 0,
        use_side_gnd: document.getElementById('chk_side_gnd').checked ? 1 : 0,
        use_top_gnd: document.getElementById('chk_top_gnd').checked ? 1 : 0,
        enclosure_width: getDisplayValue('inp_enclosure_width'),
        enclosure_height: getDisplayValue('inp_enclosure_height'),
        max_iters: parseInt(document.getElementById('inp_max_iters').value),
        // The input is in percent. The solvers (and saved settings / share links)
        // use the fraction.
        tolerance: getInputValueUnitless('inp_tolerance') / 100,
        min_converged_passes: getInputValueUnitless('inp_min_converged_passes'),
        estimate_error: document.getElementById('chk_estimate_error').checked ? 1 : 0,
        max_nodes: parseInt(document.getElementById('inp_max_nodes').value),
        rq: getDisplayValue('inp_rq'),
        use_plating: document.getElementById('chk_plating').checked ? 1 : 0,
        plating_sigma: getInputValueUnitless('inp_plating_sigma'),
        plating_t: getDisplayValue('inp_plating_t'),
        plating_rq: getDisplayValue('inp_plating_rq'),
        plating_top: document.getElementById('chk_plating_top').checked ? 1 : 0,
        plating_sides: document.getElementById('chk_plating_sides').checked ? 1 : 0,
        plating_bottom: document.getElementById('chk_plating_bottom').checked ? 1 : 0,
        plating_thick_corners: document.getElementById('chk_plating_thick_corners').checked ? 1 : 0,
        sparam_length: getDisplayValue('sparam-length'),
        sparam_z_ref: getInputValueUnitless('sparam-z-ref'),
        use_causal_materials: document.getElementById('chk_causal_materials').checked ? 1 : 0,
        interp_sweep: document.getElementById('chk_interp_sweep').checked ? 1 : 0,
        interp_tolerance: parseFloat(document.getElementById('interp_tolerance').value),
        modes_freq: getDisplayValue('modes-freq'),
        modes_nev: parseInt(document.getElementById('modes-nev').value),
        modes_mesh_density: parseInt(document.getElementById('modes-mesh-density').value),
        bs_w: getDisplayValue('inp_bs_w'),
        bs_t: getDisplayValue('inp_bs_t'),
        bs_x_offset: getDisplayValue('inp_bs_x_offset'),
        bs_sigma: getInputValueUnitless('inp_bs_sigma'),
        bs_h_bottom: getDisplayValue('inp_bs_h_bottom'),
        bs_er_bottom: getInputValueUnitless('inp_bs_er_bottom'),
        bs_tand_bottom: getInputValueUnitless('inp_bs_tand_bottom'),
        bs_h_middle: getDisplayValue('inp_bs_h_middle'),
        bs_er_middle: getInputValueUnitless('inp_bs_er_middle'),
        bs_tand_middle: getInputValueUnitless('inp_bs_tand_middle'),
        bs_h_top: getDisplayValue('inp_bs_h_top'),
        bs_er_top: getInputValueUnitless('inp_bs_er_top'),
        bs_tand_top: getInputValueUnitless('inp_bs_tand_top'),

        // Coaxial
        coax_d: getDisplayValue('inp_coax_d'),
        coax_D: getDisplayValue('inp_coax_D'),
        coax_er: getInputValueUnitless('inp_coax_er'),
        coax_tand: getInputValueUnitless('inp_coax_tand'),
        coax_sigma: getInputValueUnitless('inp_coax_sigma'),
        coax_plating_inner: document.getElementById('chk_plating_inner').checked ? 1 : 0,
        coax_plating_outer: document.getElementById('chk_plating_outer').checked ? 1 : 0,

        // Rectangular waveguide
        wg_a: getDisplayValue('inp_wg_a'),
        wg_b: getDisplayValue('inp_wg_b'),
        wg_er: getInputValueUnitless('inp_wg_er'),
        wg_tand: getInputValueUnitless('inp_wg_tand'),
        wg_sigma: getInputValueUnitless('inp_wg_sigma'),
    };
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
    'trace_spacing', 'gap', 'via_gap',
    'stripline_top_h', 'er_top', 'tand_top',
    'use_sm', 'sm_t_sub', 'sm_t_trace', 'sm_t_side', 'sm_er', 'sm_tand',
    'use_top_diel', 'top_diel_h', 'top_diel_er', 'top_diel_tand',
    'use_gnd_cut', 'gnd_cut_w', 'gnd_cut_h',
    'use_top_gnd', 'enclosure_height',
]);

// Coax: the shield IS the domain boundary, so no board-stackup option applies. Plating,
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
};
const EXCLUDED_BY_TYPE = {
    broadside_stripline: BROADSIDE_EXCLUDED_KEYS,
    coax: COAX_EXCLUDED_KEYS,
    rect_waveguide: WAVEGUIDE_EXCLUDED_KEYS,
};

function settingsToURL(settings) {
    // Exclusions only ever SHORTEN a newly generated URL; settingsFromURL merges over
    // defaults and ignores unknown keys, so existing links are unaffected.
    const excluded = EXCLUDED_BY_TYPE[settings.tl_type];
    const otherTypeKeys = new Set();
    for (const [type, keys] of Object.entries(TYPE_ONLY_KEYS)) {
        if (type === settings.tl_type) continue;
        for (const k of keys) otherTypeKeys.add(k);
    }

    // Filter out default values (and the other types' geometry fields)
    const nonDefaultSettings = {};
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
            nonDefaultSettings[key] = value;
        }
    }

    const json = JSON.stringify(nonDefaultSettings);
    // Use base64 encoding for URL-safe serialization
    return btoa(encodeURIComponent(json));
}

/**
 * Deserialize settings from URL-safe base64 string
 */
function settingsFromURL(encoded) {
    try {
        const json = decodeURIComponent(atob(encoded));
        return JSON.parse(json);
    } catch (e) {
        log('Failed to parse URL parameters:', e);
        return null;
    }
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
            if (!element || value === undefined || value === null || isNaN(value)) return;
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

        setValueWithUnit('inp_w', fullSettings.w);
        setValueWithUnit('inp_h', fullSettings.h);
        setValueWithUnit('inp_t', fullSettings.t);
        document.getElementById('inp_er').value = fullSettings.er;
        document.getElementById('inp_tand').value = fullSettings.tand;
        document.getElementById('inp_sigma').value = fullSettings.sigma;
        setValueWithUnit('freq-start', fullSettings.freq_start);
        setValueWithUnit('freq-stop', fullSettings.freq_stop);
        document.getElementById('freq-points').value = fullSettings.freq_points;
        setValueWithUnit('inp_trace_spacing', fullSettings.trace_spacing);
        setValueWithUnit('inp_gap', fullSettings.gap);
        setValueWithUnit('inp_via_gap', fullSettings.via_gap);
        setValueWithUnit('inp_air_top', fullSettings.stripline_top_h);
        document.getElementById('inp_er_top').value = fullSettings.er_top;
        document.getElementById('inp_tand_top').value = fullSettings.tand_top;

        // Checkboxes
        document.getElementById('chk_solder_mask').checked = !!fullSettings.use_sm;
        setValueWithUnit('inp_sm_t_sub', fullSettings.sm_t_sub);
        setValueWithUnit('inp_sm_t_trace', fullSettings.sm_t_trace);
        setValueWithUnit('inp_sm_t_side', fullSettings.sm_t_side);
        document.getElementById('inp_sm_er').value = fullSettings.sm_er;
        document.getElementById('inp_sm_tand').value = fullSettings.sm_tand;

        document.getElementById('chk_top_diel').checked = !!fullSettings.use_top_diel;
        setValueWithUnit('inp_top_diel_h', fullSettings.top_diel_h);
        document.getElementById('inp_top_diel_er').value = fullSettings.top_diel_er;
        document.getElementById('inp_top_diel_tand').value = fullSettings.top_diel_tand;

        document.getElementById('chk_gnd_cut').checked = !!fullSettings.use_gnd_cut;
        setValueWithUnit('inp_gnd_cut_w', fullSettings.gnd_cut_w);
        setValueWithUnit('inp_gnd_cut_h', fullSettings.gnd_cut_h);

        document.getElementById('chk_enclosure').checked = !!fullSettings.use_enclosure;
        document.getElementById('chk_side_gnd').checked = !!fullSettings.use_side_gnd;
        document.getElementById('chk_top_gnd').checked = !!fullSettings.use_top_gnd;
        setValueWithUnit('inp_enclosure_width', fullSettings.enclosure_width);
        setValueWithUnit('inp_enclosure_height', fullSettings.enclosure_height);

        document.getElementById('inp_max_iters').value = fullSettings.max_iters;
        // Stored as a fraction (settings/link format), displayed in percent.
        // toPrecision strips float noise (0.003 * 100 = 0.30000000000000004).
        document.getElementById('inp_tolerance').value =
            parseFloat((100 * fullSettings.tolerance).toPrecision(10));
        if (fullSettings.min_converged_passes !== undefined)
            document.getElementById('inp_min_converged_passes').value = fullSettings.min_converged_passes;
        document.getElementById('chk_estimate_error').checked = !!fullSettings.estimate_error;
        document.getElementById('inp_max_nodes').value = fullSettings.max_nodes;
        setValueWithUnit('inp_rq', fullSettings.rq);

        document.getElementById('chk_plating').checked = !!fullSettings.use_plating;
        document.getElementById('inp_plating_sigma').value = fullSettings.plating_sigma;
        setValueWithUnit('inp_plating_t', fullSettings.plating_t);
        setValueWithUnit('inp_plating_rq', fullSettings.plating_rq);
        document.getElementById('chk_plating_top').checked = !!fullSettings.plating_top;
        document.getElementById('chk_plating_sides').checked = !!fullSettings.plating_sides;
        document.getElementById('chk_plating_bottom').checked = !!fullSettings.plating_bottom;
        document.getElementById('chk_plating_thick_corners').checked = !!fullSettings.plating_thick_corners;

        document.getElementById('chk_causal_materials').checked = !!fullSettings.use_causal_materials;

        document.getElementById('chk_interp_sweep').checked = !!fullSettings.interp_sweep;
        document.getElementById('interp_tolerance').value = fullSettings.interp_tolerance;

        setValueWithUnit('modes-freq', fullSettings.modes_freq);
        document.getElementById('modes-nev').value = fullSettings.modes_nev;
        document.getElementById('modes-mesh-density').value = fullSettings.modes_mesh_density;

        setValueWithUnit('sparam-length', fullSettings.sparam_length);
        document.getElementById('sparam-z-ref').value = fullSettings.sparam_z_ref;

        // Broadside coupled stripline
        setValueWithUnit('inp_bs_w', fullSettings.bs_w);
        setValueWithUnit('inp_bs_t', fullSettings.bs_t);
        setValueWithUnit('inp_bs_x_offset', fullSettings.bs_x_offset);
        document.getElementById('inp_bs_sigma').value = fullSettings.bs_sigma;
        setValueWithUnit('inp_bs_h_bottom', fullSettings.bs_h_bottom);
        document.getElementById('inp_bs_er_bottom').value = fullSettings.bs_er_bottom;
        document.getElementById('inp_bs_tand_bottom').value = fullSettings.bs_tand_bottom;
        setValueWithUnit('inp_bs_h_middle', fullSettings.bs_h_middle);
        document.getElementById('inp_bs_er_middle').value = fullSettings.bs_er_middle;
        document.getElementById('inp_bs_tand_middle').value = fullSettings.bs_tand_middle;
        setValueWithUnit('inp_bs_h_top', fullSettings.bs_h_top);
        document.getElementById('inp_bs_er_top').value = fullSettings.bs_er_top;
        document.getElementById('inp_bs_tand_top').value = fullSettings.bs_tand_top;

        setValueWithUnit('inp_coax_d', fullSettings.coax_d);
        setValueWithUnit('inp_coax_D', fullSettings.coax_D);
        document.getElementById('inp_coax_er').value = fullSettings.coax_er;
        document.getElementById('inp_coax_tand').value = fullSettings.coax_tand;
        document.getElementById('inp_coax_sigma').value = fullSettings.coax_sigma;
        document.getElementById('chk_plating_inner').checked = !!fullSettings.coax_plating_inner;
        document.getElementById('chk_plating_outer').checked = !!fullSettings.coax_plating_outer;

        setValueWithUnit('inp_wg_a', fullSettings.wg_a);
        setValueWithUnit('inp_wg_b', fullSettings.wg_b);
        document.getElementById('inp_wg_er').value = fullSettings.wg_er;
        document.getElementById('inp_wg_tand').value = fullSettings.wg_tand;
        document.getElementById('inp_wg_sigma').value = fullSettings.wg_sigma;

        // tl_type is restored above the backend dropdown and the sweep checkboxes, so both
        // locks have to run once every input is in place, otherwise a stale or
        // hand-edited link can leave a full-wave-only type selected with the quasi-static
        // backend (which cannot mesh it), or leave the interpolating sweep ticked on a
        // medium that does not support it.
        if (window.enforceBackendForType) window.enforceBackendForType();
        if (window.enforceSweepOptionsForType) window.enforceSweepOptionsForType();

        return true;
    } catch (e) {
        console.error('Failed to restore settings:', e);
        return false;
    }
}

/**
 * Copy current settings as URL to clipboard
 */
function copySettingsLink() {
    const settings = getUISettings();
    const encoded = settingsToURL(settings);
    const url = `${window.location.origin}${window.location.pathname}?params=${encoded}`;

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

function log(msg) {
    _logQueue.push(msg + '\n');
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
// the old main-thread code polled `stopRequested`.
function workerStop() {
    if (_worker) _worker.postMessage({ type: 'stop' });
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
    if (fields.triMesh) target.triMesh = fields.triMesh;
    target.solution_valid = true;
    target.mesh_generated = true;
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

/**
 * Get a hash of geometry parameters for change tracking
 */
// Interpolating-sweep tolerance, as a fraction. Validated here (on the thread that owns
// the input) so a bad value fails before the worker job starts.
function interpTolerance() {
    const tolPercent = parseFloat(document.getElementById('interp_tolerance')?.value);
    if (isNaN(tolPercent) || tolPercent <= 0) {
        throw new Error("Interpolation tolerance must be a positive number.");
    }
    return tolPercent / 100;
}

function getGeometryHash() {
    const p = getParams();
    return JSON.stringify({
        tl_type: p.tl_type,
        w: p.w,
        h: p.h,
        t: p.t,
        er: p.er,
        tand: p.tand,
        sigma: p.sigma,
        trace_spacing: p.trace_spacing,
        gap: p.gap,
        via_gap: p.via_gap,
        stripline_top_h: p.stripline_top_h,
        er_top: p.er_top,
        tand_top: p.tand_top,
        use_sm: p.use_sm,
        sm_t_sub: p.sm_t_sub,
        sm_t_trace: p.sm_t_trace,
        sm_t_side: p.sm_t_side,
        sm_er: p.sm_er,
        sm_tand: p.sm_tand,
        use_top_diel: p.use_top_diel,
        top_diel_h: p.top_diel_h,
        top_diel_er: p.top_diel_er,
        top_diel_tand: p.top_diel_tand,
        use_gnd_cut: p.use_gnd_cut,
        gnd_cut_w: p.gnd_cut_w,
        gnd_cut_h: p.gnd_cut_h,
        use_enclosure: p.use_enclosure,
        use_side_gnd: p.use_side_gnd,
        use_top_gnd: p.use_top_gnd,
        enclosure_width: p.enclosure_width,
        enclosure_height: p.enclosure_height,
        rq: p.rq,
        use_plating: p.use_plating,
        plating_sigma: p.plating_sigma,
        plating_t: p.plating_t,
        plating_rq: p.plating_rq,
        plating_top: p.plating_top,
        plating_sides: p.plating_sides,
        plating_bottom: p.plating_bottom,
        bs_w: p.bs_w,
        bs_t: p.bs_t,
        bs_x_offset: p.bs_x_offset,
        bs_sigma: p.bs_sigma,
        bs_h_bottom: p.bs_h_bottom,
        bs_er_bottom: p.bs_er_bottom,
        bs_tand_bottom: p.bs_tand_bottom,
        bs_h_middle: p.bs_h_middle,
        bs_er_middle: p.bs_er_middle,
        bs_tand_middle: p.bs_tand_middle,
        bs_h_top: p.bs_h_top,
        bs_er_top: p.bs_er_top,
        bs_tand_top: p.bs_tand_top,
        coax_d: p.coax_d,
        coax_D: p.coax_D,
        coax_er: p.coax_er,
        coax_tand: p.coax_tand,
        coax_sigma: p.coax_sigma,
        coax_plating_inner: p.coax_plating_inner,
        coax_plating_outer: p.coax_plating_outer,
        wg_a: p.wg_a,
        wg_b: p.wg_b,
        wg_er: p.wg_er,
        wg_tand: p.wg_tand,
        wg_sigma: p.wg_sigma
    });
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
    const resultsNotice = document.getElementById('results-notice');
    const resultsNoticeText = document.getElementById('results-notice-text');
    const sparamNotice = document.getElementById('sparam-notice');
    const sparamNoticeText = document.getElementById('sparam-notice-text');
    const exportBtn = document.getElementById('export-snp');
    const resultsDiffCheckbox = document.getElementById('results-diff');
    const sparamDiffCheckbox = document.getElementById('sparam-diff');

    if (!frequencySweepResults || frequencySweepResults.length === 0) {
        // No results exist
        if (resultsNotice) {
            resultsNoticeText.textContent = 'No results available. Run solver to view results.';
            resultsNotice.style.display = 'block';
        }
        if (sparamNotice) {
            sparamNoticeText.textContent = 'No results available. Run solver to view S-parameters.';
            sparamNotice.style.display = 'block';
        }
        if (exportBtn) {
            exportBtn.disabled = true;
        }
        // Disable differential-mode checkboxes when no results
        if (resultsDiffCheckbox) {
            resultsDiffCheckbox.disabled = true;
        }
        if (sparamDiffCheckbox) {
            sparamDiffCheckbox.disabled = true;
        }
    } else {
        const currentGeometry = getGeometryHash();
        const currentFrequency = getFrequencyHash();
        const geometryChanged = lastSolvedGeometry && currentGeometry !== lastSolvedGeometry;
        const frequencyChanged = lastSolvedFrequency && currentFrequency !== lastSolvedFrequency;

        // Enable/disable differential-mode checkboxes based on whether results are differential
        const resultsAreDifferential = frequencySweepResults[0].result.modes.length === 2;
        if (resultsDiffCheckbox) {
            resultsDiffCheckbox.disabled = !resultsAreDifferential;
        }
        if (sparamDiffCheckbox) {
            sparamDiffCheckbox.disabled = !resultsAreDifferential;
        }

        if (!isSimulating && geometryChanged) {
            // Geometry changed - show notice but keep old results visible
            if (resultsNotice) {
                resultsNoticeText.textContent = 'Geometry changed. Solve to update results.';
                resultsNotice.style.display = 'block';
            }
            if (sparamNotice) {
                sparamNoticeText.textContent = 'Geometry changed. Solve to update results.';
                sparamNotice.style.display = 'block';
            }
            if (exportBtn) {
                exportBtn.disabled = true;
                exportBtn.title = 'Cannot export - geometry or frequency changed';
            }
        } else if (!isSimulating && frequencyChanged) {
            // Only frequency changed
            if (resultsNotice) {
                resultsNoticeText.textContent = 'Frequency changed. Solve to update results.';
                resultsNotice.style.display = 'block';
            }
            if (sparamNotice) {
                sparamNoticeText.textContent = 'Frequency changed. Solve to update results.';
                sparamNotice.style.display = 'block';
            }
            if (exportBtn) {
                exportBtn.disabled = true;
                exportBtn.title = 'Cannot export - geometry or frequency changed';
            }
        } else {
            // No changes - hide notices, enable export
            if (resultsNotice) {
                resultsNotice.style.display = 'none';
            }
            // A self-referenced medium drops its below-cutoff points (they have an
            // imaginary modal impedance), so the S-parameter tab can legitimately end up
            // with nothing to draw while the Results tab is full. Say why, rather than
            // leaving an empty plot and an export button that only fails when clicked.
            const exportable = usableSweepPoints(frequencySweepResults).length > 0;
            if (sparamNotice) {
                if (!exportable) {
                    sparamNoticeText.textContent = 'Every sweep point is below the cutoff — ' +
                        'the mode is evanescent there, so there are no propagating ' +
                        'S-parameters. Raise the sweep frequency above the cutoff.';
                    sparamNotice.style.display = 'block';
                } else {
                    sparamNotice.style.display = 'none';
                }
            }
            if (exportBtn) {
                exportBtn.disabled = !exportable;
                exportBtn.title = exportable ? '' : 'Cannot export - no propagating sweep points';
            }
        }
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

    const freq = getInputValue('modes-freq');
    let nev = parseInt(document.getElementById('modes-nev').value);
    if (!isFinite(freq) || freq <= 0) { setModesStatus('Enter a valid frequency.'); return; }
    if (!isFinite(nev) || nev < 1) { nev = 1; document.getElementById('modes-nev').value = '1'; }
    nev = Math.min(nev, 30);
    let meshDensity = parseInt(document.getElementById('modes-mesh-density').value);
    if (!isFinite(meshDensity)) meshDensity = 12;
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
        try {
            ({ grid } = await workerJob('modeField', { idx }));
        } catch (e) {
            log('Mode field fetch failed: ' + (e.message || e));
            grid = null;
        }
        // Cache successes only. Remembering a failure would leave the row permanently
        // blank with no way to retry short of re-solving. Re-asking on the next click
        // costs one worker round trip and lets a transient failure heal itself.
        if (grid) modesFieldCache.set(idx, grid);
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
        container.innerHTML = '<p style="color:var(--text-muted); padding:12px;">No field available for this mode.</p>';
        return;
    }

    const maxY = Math.max(
        modesSolver.dielectrics.reduce((m, d) => Math.max(m, d.y_max), 0),
        modesSolver.conductors.reduce((m, c) => Math.max(m, c.y_max), 0));

    const yArr = Array.from(grid.y);
    const maxYIdx = yArr.findIndex(y => y > maxY);
    const nyDisp = maxYIdx > 0 ? maxYIdx : yArr.length;
    const xMM = Array.from(grid.x, v => v * 1000);
    const yMM = yArr.slice(0, nyDisp).map(v => v * 1000);
    // A mode's eigenvector has an arbitrary scale, so absolute |E| is meaningless — show
    // the field normalized to its own maximum (0–1). This also keeps the colorbar tick
    // labels a constant width across modes, so switching modes never resizes the plot
    // area (and thus never shifts the equal-aspect x-range).
    let zmax = 0;
    for (const row of grid.E) for (const v of row) if (v > zmax) zmax = v;
    const inv = zmax > 0 ? 1 / zmax : 1;
    const z = grid.E.slice(0, nyDisp).map(row => Array.from(row, v => v * inv));

    const eeff = mode.eps_eff != null ? `, ε_eff=${mode.eps_eff.toFixed(3)}` : '';
    const STATUS_LABEL = { propagating: 'propagating', near_cutoff: 'near cutoff (evan?)', evanescent: 'evanescent', spurious: 'spurious', nullspace: 'null-space' };
    const title = `Mode ${idx} — transverse |E| (${STATUS_LABEL[mode.status] || mode.status}${eeff})`;

    // Switching modes (a plot already exists, no view reset): update only the field data
    // and title in place — leaving the axes untouched preserves the current zoom/pan
    // exactly. (Re-issuing the layout would re-run the equal-aspect scaleanchor solver and
    // drift the x-range a little each time.) The colorbar auto-rescales to the new field.
    if (!resetView && container.data && container.data.length) {
        Plotly.restyle(container, { z: [z], x: [xMM], y: [yMM] }, [0]);
        Plotly.relayout(container, { 'title.text': title });
        return;
    }

    // Fresh plot (first render or after a new solve): full render at the focused view.
    const shapes = buildGeometryShapes(maxY);
    const view = computeModesView(maxY);

    const traces = [{
        type: 'heatmap', x: xMM, y: yMM, z,
        colorscale: 'Viridis', zsmooth: 'best',
        zmin: 0, zmax: 1,
        colorbar: { title: { text: '|E_t| / max', font: { color: '#aaa' } }, tickfont: { color: '#aaa' } },
        hoverinfo: 'skip',
    }];

    Plotly.react(container, traces, modesPlotLayout(title, view, shapes),
        { responsive: true, displayModeBar: true, scrollZoom: true });
}

// Draw a geometry-only preview into the mode plot (before any solve), so the tab shows
// the structure like the reference viewer does.
function plotModesGeometry() {
    const Plotly = getPlotly();
    const container = document.getElementById('modes-plot');
    if (!Plotly || !container || !modesSolver || !modesSolver.conductors) return;
    const maxY = Math.max(
        modesSolver.dielectrics.reduce((m, d) => Math.max(m, d.y_max), 0),
        modesSolver.conductors.reduce((m, c) => Math.max(m, c.y_max), 0));
    const view = computeModesView(maxY);
    // Invisible scatter spanning the view so the axes (and shapes) scale correctly.
    const traces = [{ type: 'scatter', x: [view.xRange[0], view.xRange[1]], y: [view.yRange[0], view.yRange[1]],
        mode: 'markers', marker: { size: 0, opacity: 0 }, hoverinfo: 'skip', showlegend: false }];
    const layout = modesPlotLayout('Geometry — click Solve Modes to compute fields', view, buildGeometryShapes(maxY));
    Plotly.react(container, traces, layout, { responsive: true, displayModeBar: true, scrollZoom: true });
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
        xaxis: { title: { text: 'Width (mm)', font: { color: '#aaa' } }, scaleanchor: 'y', scaleratio: 1, range: view.xRange, color: '#aaa', gridcolor: '#444', zerolinecolor: '#555' },
        yaxis: { title: { text: 'Height (mm)', font: { color: '#aaa' } }, range: view.yRange, color: '#aaa', gridcolor: '#444', zerolinecolor: '#555' },
        margin: { l: 70, r: 90, t: 50, b: 60 },
        showlegend: false, hovermode: 'closest', dragmode: 'pan',
        paper_bgcolor: '#2a2a2a', plot_bgcolor: '#1a1a1a', font: { color: '#fff' },
        shapes,
    };
}

// Translate the "Solver" dropdown value into the backend + triangular loss options.
// Accepts the legacy 'triangular' value (maps to the accurate MQS full-wave mode).
function getParams() {
    return {
        tl_type: document.getElementById('tl_type').value,
        mesh_backend: (document.getElementById('mesh_backend')?.value) ?? 'rectilinear',
        w: getInputValue('inp_w'),
        h: getInputValue('inp_h'),
        t: getInputValue('inp_t'),
        er: getInputValueUnitless('inp_er'),
        tand: getInputValueUnitless('inp_tand'),
        sigma: getInputValueUnitless('inp_sigma'),
        freq: getInputValue('freq-start'),
        nx: 30,  // Fixed initial grid size
        ny: 30,  // Fixed initial grid size
        // Differential parameters
        trace_spacing: getInputValue('inp_trace_spacing'),
        // GCPW specific parameters
        gap: getInputValue('inp_gap'),
        via_gap: getInputValue('inp_via_gap'),
        // Stripline parameters
        stripline_top_h: getInputValue('inp_air_top'),
        er_top: getInputValueUnitless('inp_er_top'),
        tand_top: getInputValueUnitless('inp_tand_top'),
        // Solder mask parameters
        use_sm: document.getElementById('chk_solder_mask').checked,
        sm_t_sub: getInputValue('inp_sm_t_sub'),
        sm_t_trace: getInputValue('inp_sm_t_trace'),
        sm_t_side: getInputValue('inp_sm_t_side'),
        sm_er: getInputValueUnitless('inp_sm_er'),
        sm_tand: getInputValueUnitless('inp_sm_tand'),
        // Top dielectric parameters
        use_top_diel: document.getElementById('chk_top_diel').checked,
        top_diel_h: getInputValue('inp_top_diel_h'),
        top_diel_er: getInputValueUnitless('inp_top_diel_er'),
        top_diel_tand: getInputValueUnitless('inp_top_diel_tand'),
        // Ground cutout parameters
        use_gnd_cut: document.getElementById('chk_gnd_cut').checked,
        gnd_cut_w: getInputValue('inp_gnd_cut_w'),
        gnd_cut_h: getInputValue('inp_gnd_cut_h'),
        // Enclosure parameters
        use_enclosure: document.getElementById('chk_enclosure').checked,
        use_side_gnd: document.getElementById('chk_side_gnd').checked,
        use_top_gnd: document.getElementById('chk_top_gnd').checked,
        enclosure_width: getInputValue('inp_enclosure_width'),
        enclosure_height: getInputValue('inp_enclosure_height'),
        max_iters: parseInt(document.getElementById('inp_max_iters').value),
        // Percent in the UI, fraction in the settings object, keeps saved settings
        // and share links (which diff against DEFAULT_SETTINGS.tolerance) unchanged.
        tolerance: getInputValueUnitless('inp_tolerance') / 100,
        min_converged_passes: getInputValueUnitless('inp_min_converged_passes'),
        // 1/0 like every other checkbox here and in getUISettings, the same key must not
        // be a boolean in one params object and a number in the other.
        estimate_error: document.getElementById('chk_estimate_error').checked ? 1 : 0,
        max_nodes: parseInt(document.getElementById('inp_max_nodes').value),
        // Surface roughness parameter
        rq: getInputValue('inp_rq'),
        // Surface plating parameters
        use_plating: document.getElementById('chk_plating').checked,
        plating_sigma: getInputValueUnitless('inp_plating_sigma'),
        plating_t: getInputValue('inp_plating_t'),
        plating_rq: getInputValue('inp_plating_rq'),
        plating_top: document.getElementById('chk_plating_top').checked,
        plating_sides: document.getElementById('chk_plating_sides').checked,
        plating_bottom: document.getElementById('chk_plating_bottom').checked,
        plating_thick_corners: document.getElementById('chk_plating_thick_corners').checked,
        // Causal material parameters
        use_causal_materials: document.getElementById('chk_causal_materials').checked,
        // Broadside coupled stripline parameters
        bs_w: getInputValue('inp_bs_w'),
        bs_t: getInputValue('inp_bs_t'),
        bs_x_offset: getInputValue('inp_bs_x_offset'),
        bs_sigma: getInputValueUnitless('inp_bs_sigma'),
        bs_h_bottom: getInputValue('inp_bs_h_bottom'),
        bs_er_bottom: getInputValueUnitless('inp_bs_er_bottom'),
        bs_tand_bottom: getInputValueUnitless('inp_bs_tand_bottom'),
        bs_h_middle: getInputValue('inp_bs_h_middle'),
        bs_er_middle: getInputValueUnitless('inp_bs_er_middle'),
        bs_tand_middle: getInputValueUnitless('inp_bs_tand_middle'),
        bs_h_top: getInputValue('inp_bs_h_top'),
        bs_er_top: getInputValueUnitless('inp_bs_er_top'),
        bs_tand_top: getInputValueUnitless('inp_bs_tand_top'),

        // Coaxial (diameters in; CoaxSolver derives the radii)
        coax_d: getInputValue('inp_coax_d'),
        coax_D: getInputValue('inp_coax_D'),
        coax_er: getInputValueUnitless('inp_coax_er'),
        coax_tand: getInputValueUnitless('inp_coax_tand'),
        coax_sigma: getInputValueUnitless('inp_coax_sigma'),
        coax_plating_inner: document.getElementById('chk_plating_inner').checked,
        coax_plating_outer: document.getElementById('chk_plating_outer').checked,

        // Rectangular waveguide (inner wall dimensions)
        wg_a: getInputValue('inp_wg_a'),
        wg_b: getInputValue('inp_wg_b'),
        wg_er: getInputValueUnitless('inp_wg_er'),
        wg_tand: getInputValueUnitless('inp_wg_tand'),
        wg_sigma: getInputValueUnitless('inp_wg_sigma'),
    };
}

// Helper function to add common optional geometry parameters
function updateGeometry() {
    setCurrentView("geometry");

    const pbar = document.getElementById('progress_bar');
    pbar.style.width = "0%";

    solver = buildSolverFromParams(getParams());
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
    // The view model these fields belong to. A mid-solve edit runs
    // updateGeometry(), which replaces `solver` with one describing a different
    // cross-section. Identity is the exact test: updateGeometry() always
    // assigns a fresh object.
    const solvedSolver = solver;

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
    isSimulating = true;
    updateResultNotices();
    displayedProgress = 0;
    pbar.style.width = '0%';
    if (ptext) ptext.textContent = '';
    heartbeatStart(ptext);
    log("Starting simulation...");
    for (const w of solver.openBoundaryWarnings()) log(`⚠ Warning: ${w}`);
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

        const out = await workerJob('simulate', {
            params: p,
            frequencies,
            opts: {
                useInterpolation: !!document.getElementById('chk_interp_sweep')?.checked,
                interpTolerance: interpTolerance(),
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
                    log(`Estimated remaining error ${(100 * m.meta.certification.err).toExponential(2)}% ` +
                        `(tolerance ${(100 * p.tolerance).toFixed(2)}%)`);
                }
                if (m.meta.meshQuality && m.meta.meshQuality.maxQ > 100) {
                    log(`\u26a0 Mesh quality warning: worst triangle Q=${m.meta.meshQuality.maxQ.toFixed(0)} ` +
                        `results may be inaccurate. Try adjusting geometry or enclosure size.`);
                }
                logModeWarnings(m.warnings);
            },
            warnings: (m) => logModeWarnings(m.warnings),
            // Live plots while the sweep runs. Redrawing is cheap next to a solve, and
            // now that the solve is off-thread these actually paint.
            partial: (m) => {
                if (m.fields && solver === solvedSolver) { applyFields(solver, m.fields); draw(); }
                if (m.sweepResults) {
                    frequencySweepResults = reviveSweep(m.sweepResults);
                    drawResultsPlot();
                    drawSParamPlot();
                }
            },
        });

        if (solver === solvedSolver) applyFields(solver, out.fields);
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

        // Display summary
        const f0 = frequencies[0] / 1e9;
        const mode0 = frequencySweepResults[0].result.modes[0];
        const loss0 = mode0.alpha_total;
        const isSingleFreq = frequencies.length === 1;

        // Check if differential results
        if (results.modes.length === 2) {
            const odd = results.modes.find(m => m.mode === 'odd');
            const even = results.modes.find(m => m.mode === 'even');
            let lossStr;
            if (isSingleFreq) {
                lossStr = `Loss: ${loss0.toFixed(3)} dB/m @ ${f0.toFixed(2)} GHz`;
            } else {
                const fn = frequencies[frequencies.length - 1] / 1e9;
                const lossN = frequencySweepResults[frequencySweepResults.length - 1].result.modes[0].alpha_total;
                lossStr = `Loss: ${loss0.toFixed(3)} dB/m @ ${f0.toFixed(2)} GHz - ${lossN.toFixed(3)} dB/m @ ${fn.toFixed(2)} GHz`;
            }
            // For an asymmetric pair the two traces are not interchangeable: the physical self
            // terms differ (C11 ≠ C22) and S22 ≠ S11. Surface that here — otherwise the summary
            // looks identical to a symmetric line. odd/even are then only the approximate eigenmodes.
            const mC = results.RLGC_matrix?.C, mL = results.RLGC_matrix?.L;
            const asymStr = (results.physMatrix && mC && mL)
                ? `\n\nAsymmetric pair (S22 ≠ S11; odd/even are the approximate eigenmodes):\n` +
                  `  Self-C:  C11 = ${(mC[0][0] * 1e12).toFixed(2)} pF/m,  C22 = ${(mC[1][1] * 1e12).toFixed(2)} pF/m\n` +
                  `  Self-L:  L11 = ${(mL[0][0] * 1e9).toFixed(2)} nH/m,  L22 = ${(mL[1][1] * 1e9).toFixed(2)} nH/m`
                : '';
            log(`\nDIFFERENTIAL RESULTS:\n` +
                     `======================\n` +
                     `Differential Impedance Z_diff: ${results.Z_diff.toFixed(2)} Ohm  (2 x Z_odd)\n` +
                     `Common-Mode Impedance Z_common: ${results.Z_common.toFixed(2)} Ohm  (Z_even / 2)\n` +
                     `\nModal Impedances:\n` +
                     `  Odd-Mode  Z_odd:  ${odd.Z0.toFixed(2)} Ohm  (eps_eff = ${odd.eps_eff.toFixed(3)})\n` +
                     `  Even-Mode Z_even: ${even.Z0.toFixed(2)} Ohm  (eps_eff = ${even.eps_eff.toFixed(3)})` +
                     `${asymStr}\n` +
                     `\n${lossStr}`);
        } else {
            let lossStr;
            if (isSingleFreq) {
                lossStr = `Loss: ${loss0.toFixed(3)} dB/m @ ${f0.toFixed(2)} GHz`;
            } else {
                const fn = frequencies[frequencies.length - 1] / 1e9;
                const lossN = frequencySweepResults[frequencySweepResults.length - 1].result.modes[0].alpha_total;
                lossStr = `Loss: ${loss0.toFixed(3)} dB/m @ ${f0.toFixed(2)} GHz - ${lossN.toFixed(3)} dB/m @ ${fn.toFixed(2)} GHz`;
            }
            // Z0/eps_eff are quoted at the first point, which for a waveguide can be below
            // cutoff (where they are NaN by design). Quote the first propagating point
            // instead, and say which one it is, rather than printing "NaN Ohm".
            let sumMode = mode0, sumF = f0, cutoffNote = '';
            if (Number.isNaN(mode0.Z0)) {
                const firstProp = frequencySweepResults.find(r => !Number.isNaN(r.result.modes[0].Z0));
                if (firstProp) {
                    sumMode = firstProp.result.modes[0];
                    sumF = firstProp.freq / 1e9;
                    cutoffNote = ` (at ${sumF.toFixed(2)} GHz — the sweep starts below cutoff)`;
                }
            }
            log(`\nRESULTS:\n` +
                     `----------------------\n` +
                     (Number.isNaN(sumMode.Z0)
                        ? `Below cutoff across the whole sweep — attenuation only.\n`
                        : `Z0: ${sumMode.Z0.toFixed(2)} Ohm${cutoffNote}\n` +
                          `eps_eff: ${sumMode.eps_eff.toFixed(3)}\n`) +
                     `${lossStr}`);
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
        // Restore button to "Solve" mode
        btn.textContent = 'Solve';
        btn.classList.remove('stop-mode');
        // A full bar and a "Done" only for a run that actually produced results; a
        // cancelled or failed solve rewinds the bar instead of claiming completion.
        pbar.style.width = outcome === 'done' ? '100%' : '0%';
        // Keep the elapsed total on screen rather than blanking it.
        const elapsed = _hbFormat(Date.now() - _hbStart);
        heartbeatStop(outcome === 'done' ? `Done in ${elapsed}`
                    : outcome === 'stopped' ? `Stopped after ${elapsed}` : '');
        isSimulating = false;
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
    const cfg = SWEEP_PARAM_CONFIG[lastSweepParam];
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
    // Both replace the microstrip stackup entirely with their own geometry block.
    const isSelfBounded = isCoax || isWaveguide;
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
        sm: useSm && !isSelfBounded,
        top_diel: useTopDiel && !isSelfBounded,
        gnd_cut: useGndCut && !isSelfBounded,
        enclosure: useEnclosure && !isSelfBounded,
        plating: usePlating,
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
    // Restore previous selection if still available
    if ([...sel.options].some(o => o.value === previousValue)) sel.value = previousValue;
}

function getZeroDefaultMax(displayUnit) {
    // Return a sensible max in display units for zero-valued params
    const maxInMeters = 2e-6; // 2 μm as reference
    if (!displayUnit) return 1; // unitless
    return +(convertToDisplayUnit(maxInMeters, displayUnit)).toPrecision(4);
}

function autoFillSweepRange() {
    const xSel = document.getElementById('sweep-x-selector').value;
    const cfg = SWEEP_PARAM_CONFIG[xSel];
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
        'mm': 1e3, 'μm': 1e6, 'um': 1e6, 'nm': 1e9,
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
    const cfg = SWEEP_PARAM_CONFIG[xSel];
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

        const out = await workerJob('paramSweep', {
            points, freqHz, opts: { estimateError: !!p.estimate_error },
        }, {
            // Surface the first point's verification outcome once. Later points run the
            // legacy per-pass gate with the same refinement settings.
            firstPointCert: (m) => {
                if (m.certification && m.certification.pass) {
                    log(`Estimated remaining error (first-point) ${(100 * m.certification.err).toExponential(2)}%. ` +
                        `Later sweep points are not re-verified.`);
                }
                // Every accuracy note, not just the first: one solve can carry a failed
                // certificate and a loss-accuracy note (skin-transition,
                // broadside-proximity) at once.
                for (const aw of m.warnings || []) {
                    if (aw.type !== 'accuracy') continue;
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
        solveBtn.disabled = false;
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

    // Results plot selector change
    const resultsSelector = document.getElementById('results-plot-selector');
    if (resultsSelector) {
        resultsSelector.addEventListener('change', () => {
            if (frequencySweepResults) {
                drawResultsPlot();
            }
        });
    }

    // S-parameter controls
    const sparamLength = document.getElementById('sparam-length');
    const sparamZref = document.getElementById('sparam-z-ref');
    const sparamMode = document.getElementById('sparam-plot-mode');
    const sparamDiff = document.getElementById('sparam-diff');
    if (sparamLength) {
        sparamLength.addEventListener('input', () => {
            if (frequencySweepResults) {
                drawSParamPlot();
            }
        });
    }
    if (sparamZref) {
        sparamZref.addEventListener('input', () => {
            if (frequencySweepResults) {
                drawSParamPlot();
            }
        });
    }
    if (sparamMode) {
        sparamMode.addEventListener('change', () => {
            if (frequencySweepResults) {
                drawSParamPlot();
            }
        });
    }
    if (sparamDiff) {
        sparamDiff.addEventListener('change', () => {
            if (frequencySweepResults) {
                drawSParamPlot();
            }
        });
    }

    // Log checkbox for results plot
    const resultsLogX = document.getElementById('results-log-x');
    if (resultsLogX) {
        resultsLogX.addEventListener('change', () => {
            if (frequencySweepResults) {
                drawResultsPlot();
            }
        });
    }

    // Differential-mode checkbox for results plot
    const resultsDiff = document.getElementById('results-diff');
    if (resultsDiff) {
        resultsDiff.addEventListener('change', () => {
            if (frequencySweepResults) {
                drawResultsPlot();
            }
        });
    }

    // Log checkbox for S-parameter plot
    const sparamLogX = document.getElementById('sparam-log-x');
    if (sparamLogX) {
        sparamLogX.addEventListener('change', () => {
            if (frequencySweepResults) {
                drawSParamPlot();
            }
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
                traceSpacing: p.trace_spacing,
                surfaceRoughness: p.rq,
                // Coax and waveguide are not trace-on-substrate stackups, so they
                // describe themselves.
                geometryLines: p.tl_type === 'coax' ? [
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
                plating: platingOptions(p,
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
            const rows = [];
            if (isDifferential) {
                rows.push([
                    'Freq_Hz',
                    'Re_Z0_odd_Ohm', 'Im_Z0_odd_Ohm', 'eps_eff_odd',
                    'conductor_loss_odd_dBpm', 'dielectric_loss_odd_dBpm', 'total_loss_odd_dBpm',
                    'R_odd_Ohmpm', 'L_odd_Hpm', 'G_odd_Spm', 'C_odd_Fpm',
                    'Re_Z0_even_Ohm', 'Im_Z0_even_Ohm', 'eps_eff_even',
                    'conductor_loss_even_dBpm', 'dielectric_loss_even_dBpm', 'total_loss_even_dBpm',
                    'R_even_Ohmpm', 'L_even_Hpm', 'G_even_Spm', 'C_even_Fpm'
                ]);
                for (const { freq, result } of frequencySweepResults) {
                    const m0 = result.modes[0];
                    const m1 = result.modes[1];
                    rows.push([
                        freq,
                        m0.Zc.re, m0.Zc.im, m0.eps_eff,
                        m0.alpha_c, m0.alpha_d, m0.alpha_total,
                        m0.RLGC.R, m0.RLGC.L, m0.RLGC.G, m0.RLGC.C,
                        m1.Zc.re, m1.Zc.im, m1.eps_eff,
                        m1.alpha_c, m1.alpha_d, m1.alpha_total,
                        m1.RLGC.R, m1.RLGC.L, m1.RLGC.G, m1.RLGC.C
                    ]);
                }
            } else {
                rows.push([
                    'Freq_Hz',
                    'Re_Z0_Ohm', 'Im_Z0_Ohm', 'eps_eff',
                    'conductor_loss_dBpm', 'dielectric_loss_dBpm', 'total_loss_dBpm',
                    'R_Ohmpm', 'L_Hpm', 'G_Spm', 'C_Fpm'
                ]);
                for (const { freq, result } of frequencySweepResults) {
                    const m = result.modes[0];
                    rows.push([
                        freq,
                        m.Zc.re, m.Zc.im, m.eps_eff,
                        m.alpha_c, m.alpha_d, m.alpha_total,
                        m.RLGC.R, m.RLGC.L, m.RLGC.G, m.RLGC.C
                    ]);
                }
            }
            const csv = rows.map(r => r.join(',')).join('\n');
            const blob = new Blob([csv], { type: 'text/csv' });
            const url = URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = 'results.csv';
            a.click();
            URL.revokeObjectURL(url);
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

    // Real-time geometry updates for all parameter inputs
    const geometryInputs = [
        'inp_w', 'inp_h', 'inp_t', 'inp_er', 'inp_tand', 'inp_sigma',
        'inp_trace_spacing',
        'inp_gap', 'inp_via_gap',
        'inp_air_top', 'inp_er_top', 'inp_tand_top',
        'inp_sm_t_sub', 'inp_sm_t_trace', 'inp_sm_t_side', 'inp_sm_er', 'inp_sm_tand',
        'inp_top_diel_h', 'inp_top_diel_er', 'inp_top_diel_tand',
        'inp_gnd_cut_w', 'inp_gnd_cut_h',
        'inp_enclosure_width', 'inp_enclosure_height',
        'inp_rq',
        'inp_plating_sigma', 'inp_plating_t', 'inp_plating_rq',
        'inp_bs_w', 'inp_bs_t', 'inp_bs_x_offset', 'inp_bs_sigma',
        'inp_bs_h_bottom', 'inp_bs_er_bottom', 'inp_bs_tand_bottom',
        'inp_bs_h_middle', 'inp_bs_er_middle', 'inp_bs_tand_middle',
        'inp_bs_h_top', 'inp_bs_er_top', 'inp_bs_tand_top',
        'inp_coax_d', 'inp_coax_D', 'inp_coax_er', 'inp_coax_tand', 'inp_coax_sigma',
        'inp_wg_a', 'inp_wg_b', 'inp_wg_er', 'inp_wg_tand', 'inp_wg_sigma',
        'freq-start'
    ];

    geometryInputs.forEach(id => {
        const el = document.getElementById(id);
        if (el) {
            el.addEventListener('input', () => {
                updateGeometry();
                draw();
                updateResultNotices();
                updateSweepNotice();
                updateModesNotice();
            });
        }
    });

    // Real-time updates for checkboxes
    const geometryCheckboxes = [
        'chk_solder_mask', 'chk_top_diel', 'chk_gnd_cut', 'chk_enclosure', 'chk_side_gnd', 'chk_top_gnd',
        'chk_plating', 'chk_plating_top', 'chk_plating_sides', 'chk_plating_bottom', 'chk_plating_thick_corners',
        'chk_plating_inner', 'chk_plating_outer'
    ];

    geometryCheckboxes.forEach(id => {
        const el = document.getElementById(id);
        if (el) {
            el.addEventListener('change', () => {
                updateGeometry();
                draw();
                updateResultNotices();
                updateSweepNotice();
                updateModesNotice();
            });
        }
    });

    // Transmission line type selector - reset zoom when type changes
    document.getElementById('tl_type').addEventListener('change', () => {
        updateGeometry();
        draw(true);  // Reset zoom/pan for new geometry
        updateResultNotices();
        updateSweepNotice();
        updateModesNotice();
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

    // Plot options - mode selector
    const plotModeEl = document.getElementById('plot-mode');
    if (plotModeEl) {
        plotModeEl.addEventListener('change', () => {
            if (solver && solver.solution_valid) {
                draw();
            }
        });
    }

    // Plot options - streamlines and contours
    const plotStreamlinesEl = document.getElementById('plot-streamlines');
    const plotContoursEl = document.getElementById('plot-contours');
    if (plotStreamlinesEl) {
        plotStreamlinesEl.addEventListener('change', () => {
            if (solver && solver.solution_valid) {
                draw();
            }
        });
    }
    if (plotContoursEl) {
        plotContoursEl.addEventListener('change', () => {
            if (solver && solver.solution_valid) {
                draw();
            }
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
    geometry: { min: null, max: null }
};

let scaleDialogOpen = false;

function getViewType(view) {
    if (view === 'potential') return 'potential';
    if (view.startsWith('efield')) return 'efield';
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
    } else {
        // For E-field and geometry, use 1.5x actual data range for margin
        sliderMinBound = actualMin < -0.1 ? actualMin * 1.5 : 0.0;
        sliderMaxBound = actualMax * 1.5;
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
    ['chk_solder_mask', 'chk_top_diel', 'chk_gnd_cut', 'chk_enclosure', 'chk_plating'].forEach(id => {
        const checkbox = document.getElementById(id);
        if (checkbox) {
            const sectionId = id.replace('chk_', '') + '-params';
            const section = document.getElementById(sectionId);
            if (section) {
                section.style.display = checkbox.checked ? 'block' : 'none';
            }
        }
    });

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

    updateGeometry();
    draw();
    resizeCanvas();
    window.addEventListener('resize', resizeCanvas);
    log("Ready. Click 'Solve' to start simulation.");
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
