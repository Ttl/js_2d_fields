import { MicrostripSolver } from './microstrip.js';
import { BroadsideStriplineSolver } from './broadside_stripline.js';
import { computeSParamsSingleEnded, computeSParamsDifferential, sParamTodB } from './sparameters.js';
import { exportSnP } from './snp_export.js';
import { draw, drawResultsPlot, drawSParamPlot, drawParameterSweepPlot, setGlobals, setCurrentView, getScaleRange, setScaleRange, getActualDataRange,
    freeze, unfreeze, isFrozen } from './plot.js';
import { InterpolatingSweep } from './interpolating_sweep.js';

// Lazy Plotly access - allows app to function while Plotly is loading
const getPlotly = () => window.Plotly;

let solver = null;
let stopRequested = false;
let isSimulating = false;
let frequencySweepResults = null;  // Array of {freq, result} objects
let currentTab = 'geometry';
let geometryChanged = false;  // Track if geometry has changed since last solve
let lastSolvedGeometry = null;  // Hash of geometry params from last solve
let lastSolvedFrequency = null;  // Frequency params from last solve
let sweepStopRequested = false;
let isSweeping = false;
let parameterSweepResults = null;
let lastSweepGeometry = null;  // Geometry hash at sweep time (excluding swept param)
let lastSweepParam = null;     // Which parameter was swept
let lastSweepDisplayUnit = null; // Display unit used during last sweep

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
    rq:                 { label: 'Surface Roughness',    inputId: 'inp_rq',            group: 'always' },
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
        tolerance: getInputValueUnitless('inp_tolerance'),
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
    };
}

/**
 * Serialize settings to URL-safe base64 string
 * Only includes non-default parameters to keep URLs short
 */
function settingsToURL(settings) {
    // Broadside coupled stripline uses an allowlist. Only its own fields plus
    // shared frequency/sparam/enclosure/plating fields are written into the link.
    if (settings.tl_type === 'broadside_stripline') {
        const allow = new Set([
            'tl_type',
            'bs_w', 'bs_t', 'bs_x_offset', 'bs_sigma',
            'bs_h_bottom', 'bs_er_bottom', 'bs_tand_bottom',
            'bs_h_middle', 'bs_er_middle', 'bs_tand_middle',
            'bs_h_top', 'bs_er_top', 'bs_tand_top',
            'freq_start', 'freq_stop', 'freq_points',
            'sparam_length', 'sparam_z_ref',
            'use_enclosure', 'use_side_gnd', 'enclosure_width',
            'use_plating', 'plating_sigma', 'plating_t', 'plating_rq',
            'plating_top', 'plating_sides', 'plating_bottom', 'plating_thick_corners',
            'use_causal_materials', 'interp_sweep', 'interp_tolerance',
            'max_iters', 'tolerance', 'max_nodes',
        ]);
        const out = {};
        for (const key in settings) {
            if (!allow.has(key)) continue;
            const value = settings[key];
            const defaultValue = DEFAULT_SETTINGS[key];
            const bothNaN = (typeof value === 'number' && isNaN(value)) &&
                            (typeof defaultValue === 'number' && isNaN(defaultValue));
            if (bothNaN) continue;
            if (value !== defaultValue) out[key] = value;
        }
        // Always include tl_type so the link round-trips.
        out.tl_type = 'broadside_stripline';
        const json = JSON.stringify(out);
        return btoa(encodeURIComponent(json));
    }

    // Filter out default values
    const nonDefaultSettings = {};
    for (const key in settings) {
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
        document.getElementById('inp_tolerance').value = fullSettings.tolerance;
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

function log(msg) {
    const c = document.getElementById('console_out');
    c.textContent += msg + "\n";
    c.scrollTop = c.scrollHeight;
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
        bs_tand_top: p.bs_tand_top
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
            if (sparamNotice) {
                sparamNotice.style.display = 'none';
            }
            if (exportBtn) {
                exportBtn.disabled = false;
                exportBtn.title = '';
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
    }
}

function getParams() {
    return {
        tl_type: document.getElementById('tl_type').value,
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
        tolerance: getInputValueUnitless('inp_tolerance'),
        min_converged_passes: getInputValueUnitless('inp_min_converged_passes'),
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
    };
}

// Helper function to add common optional geometry parameters
function addCommonOptions(options, p) {
    // Solder mask
    if (p.use_sm) {
        options.use_sm = true;
        options.sm_t_sub = p.sm_t_sub;
        options.sm_t_trace = p.sm_t_trace;
        options.sm_t_side = p.sm_t_side;
        options.sm_er = p.sm_er;
        options.sm_tand = p.sm_tand;
    }

    // Top dielectric
    if (p.use_top_diel) {
        options.top_diel_h = p.top_diel_h;
        options.top_diel_er = p.top_diel_er;
        options.top_diel_tand = p.top_diel_tand;
    }

    // Ground cutout
    if (p.use_gnd_cut) {
        options.gnd_cut_width = p.gnd_cut_w;
        options.gnd_cut_sub_h = p.gnd_cut_h;
    }

    // Enclosure
    if (p.use_enclosure) {
        options.enclosure_width = p.enclosure_width;
        if (options.enclosure_height === undefined) {
            options.enclosure_height = p.enclosure_height;
        }

        const left_bc = p.use_side_gnd ? "gnd" : "open";
        const right_bc = p.use_side_gnd ? "gnd" : "open";
        const top_bc = p.use_top_gnd ? "gnd" : (options.boundaries ? options.boundaries[2] : "open");
        const bottom_bc = options.boundaries ? options.boundaries[3] : "gnd";
        options.boundaries = [left_bc, right_bc, top_bc, bottom_bc];
    }

    // Surface plating
    if (p.use_plating) {
        options.plating = {
            sigma: p.plating_sigma,
            thickness: p.plating_t,
            rq: p.plating_rq,
            top: p.plating_top,
            sides: p.plating_sides,
            bottom: p.plating_bottom,
            thick_corners: p.plating_thick_corners
        };
    }
}

function updateGeometry() {
    const p = getParams();
    setCurrentView("geometry");

    const pbar = document.getElementById('progress_bar');
    pbar.style.width = "0%";

    try {
        if (p.tl_type === 'gcpw') {
            const options = {
                substrate_height: p.h,
                trace_width: p.w,
                trace_thickness: p.t,
                gnd_thickness: 35e-6,
                epsilon_r: p.er,
                tan_delta: p.tand,
                sigma_cond: p.sigma,
                freq: p.freq,
                nx: p.nx,
                ny: p.ny,
                boundaries: ["open", "open", "open", "gnd"],
                // Coplanar-specific
                use_coplanar_gnd: true,
                gap: p.gap,
                via_gap: p.via_gap,
                use_vias: true,
                // Surface roughness
                rq: p.rq,
            };
            addCommonOptions(options, p);
            solver = new MicrostripSolver(options);
        } else if (p.tl_type === 'diff_gcpw') {
            const options = {
                substrate_height: p.h,
                trace_width: p.w,
                trace_thickness: p.t,
                trace_spacing: p.trace_spacing,  // Enables differential mode
                gnd_thickness: 35e-6,
                epsilon_r: p.er,
                tan_delta: p.tand,
                sigma_cond: p.sigma,
                freq: p.freq,
                nx: p.nx,
                ny: p.ny,
                boundaries: ["open", "open", "open", "gnd"],
                // Coplanar-specific
                use_coplanar_gnd: true,
                gap: p.gap,
                via_gap: p.via_gap,
                use_vias: true,
                // Surface roughness
                rq: p.rq,
            };
            addCommonOptions(options, p);
            solver = new MicrostripSolver(options);
        } else if (p.tl_type === 'diff_microstrip') {
            // Differential Microstrip
            const options = {
                trace_width: p.w,
                substrate_height: p.h,
                trace_thickness: p.t,
                trace_spacing: p.trace_spacing,  // Enable differential mode
                epsilon_r: p.er,
                tan_delta: p.tand,
                sigma_cond: p.sigma,
                freq: p.freq,
                nx: p.nx,
                ny: p.ny,
                boundaries: ["open", "open", "open", "gnd"],
                // Surface roughness
                rq: p.rq
            };
            addCommonOptions(options, p);
            solver = new MicrostripSolver(options);
        } else if (p.tl_type === 'stripline') {
            const options = {
                trace_width: p.w,
                substrate_height: p.h,
                trace_thickness: p.t,
                epsilon_r: p.er,
                epsilon_r_top: p.er_top,
                tan_delta_top: p.tand_top,
                enclosure_height: p.stripline_top_h,
                tan_delta: p.tand,
                sigma_cond: p.sigma,
                freq: p.freq,
                nx: p.nx,
                ny: p.ny,
                boundaries: ["open", "open", "gnd", "gnd"],
                // Surface roughness
                rq: p.rq
            };
            addCommonOptions(options, p);
            solver = new MicrostripSolver(options);
        } else if (p.tl_type === 'broadside_stripline') {
            const options = {
                trace_width: p.bs_w,
                trace_thickness: p.bs_t,
                x_offset: p.bs_x_offset,
                sigma_cond: p.bs_sigma,
                h_bottom: p.bs_h_bottom,
                er_bottom: p.bs_er_bottom,
                tand_bottom: p.bs_tand_bottom,
                h_middle: p.bs_h_middle,
                er_middle: p.bs_er_middle,
                tand_middle: p.bs_tand_middle,
                h_top: p.bs_h_top,
                er_top: p.bs_er_top,
                tand_top: p.bs_tand_top,
                freq: p.freq,
                nx: p.nx,
                ny: p.ny,
                rq: p.rq,
                boundaries: ["open", "open", "gnd", "gnd"],
            };
            // Enclosure: only side ground walls apply (top/bottom are intrinsic).
            if (p.use_enclosure) {
                options.enclosure_width = p.enclosure_width;
                if (p.use_side_gnd) {
                    options.boundaries = ["gnd", "gnd", "gnd", "gnd"];
                }
            }
            // Plating
            if (p.use_plating) {
                options.plating = {
                    sigma: p.plating_sigma,
                    thickness: p.plating_t,
                    rq: p.plating_rq,
                    top: p.plating_top,
                    sides: p.plating_sides,
                    bottom: p.plating_bottom,
                    thick_corners: p.plating_thick_corners
                };
            }
            solver = new BroadsideStriplineSolver(options);
        } else if (p.tl_type === 'diff_stripline') {
            const options = {
                trace_width: p.w,
                substrate_height: p.h,
                trace_thickness: p.t,
                trace_spacing: p.trace_spacing,  // Enable differential mode
                epsilon_r: p.er,
                epsilon_r_top: p.er_top,
                enclosure_height: p.stripline_top_h,
                tan_delta: p.tand,
                tan_delta_top: p.tand_top,
                sigma_cond: p.sigma,
                freq: p.freq,
                nx: p.nx,
                ny: p.ny,
                boundaries: ["open", "open", "gnd", "gnd"],
                // Surface roughness
                rq: p.rq
            };
            addCommonOptions(options, p);
            solver = new MicrostripSolver(options);
        } else {
            // Microstrip (with optional solder mask, top dielectric, ground cutout)
            const options = {
                trace_width: p.w,
                substrate_height: p.h,
                trace_thickness: p.t,
                epsilon_r: p.er,
                tan_delta: p.tand,
                sigma_cond: p.sigma,
                freq: p.freq,
                nx: p.nx,
                ny: p.ny,
                boundaries: ["open", "open", "open", "gnd"],
                // Surface roughness
                rq: p.rq
            };
            addCommonOptions(options, p);
            solver = new MicrostripSolver(options);
        }

        // Store causal materials option on solver
        if (solver) {
            solver.use_causal_materials = p.use_causal_materials;
        }
    } catch (error) {
        // Log validation errors to the console
        log('ERROR: ' + error.message);
        // Set solver to null to prevent simulation from running with invalid parameters
        solver = null;
    }
}

async function runSimulation() {
    // Check if solver is valid before attempting to run simulation
    if (!solver) {
        log("ERROR: Cannot run simulation - solver initialization failed due to invalid parameters.");
        return;
    }

    const p = getParams();
    const frequencies = getFrequencies();
    const btn = document.getElementById('btn_solve');
    const pbar = document.getElementById('progress_bar');
    const ptext = document.getElementById('progress_text');

    // Change button to "Stop" mode
    btn.textContent = 'Stop';
    btn.classList.add('stop-mode');
    stopRequested = false;
    isSimulating = true;
    updateResultNotices();
    pbar.style.width = '0%';
    if (ptext) ptext.style.display = 'block';
    log("Starting simulation...");

    try {

        if (p.sigma < 1e4) {
            throw new Error("Signal line conductivity is too low to be considered a conductor.");
        }

        // Validate surface roughness
        if (p.rq < 0) {
            throw new Error("Surface roughness cannot be negative.");
        }

        // Validate plating parameters
        if (p.use_plating) {
            if (p.plating_sigma < 1e4) {
                throw new Error("Plating conductivity is too low to be considered a conductor.");
            }
            if (p.plating_t < 0) {
                throw new Error("Plating thickness must be non-negative.");
            }
            if (p.plating_rq < 0) {
                throw new Error("Plating roughness cannot be negative.");
            }
        }

        // Clear previous results
        frequencySweepResults = [];

        // Use the highest frequency for mesh generation (skin depth calculation)
        const maxFreq = Math.max(...frequencies);
        solver.freq = maxFreq;

        // Ensure mesh is generated before solving
        log("Calculating mesh...");
        solver.ensure_mesh();
        log("Mesh generated: " + solver.x.length + "x" + solver.y.length);

        // Run adaptive refinement at highest frequency first
        // Note: Causal materials will be applied during frequency sweep in computeAtFrequency()
        log(`Running adaptive analysis (max ${p.max_iters} iterations, max ${p.max_nodes} nodes, tolerance ${p.tolerance})...`);

        let results = await solver.solve_adaptive({
            max_iters: p.max_iters,
            energy_tol: p.tolerance,
            param_tol: 0.05,
            max_nodes: p.max_nodes*1000,
            min_converged_passes: p.min_converged_passes,
            onProgress: (info) => {
                const progress = info.iteration / p.max_iters * 0.5;  // First half is for mesh refinement
                pbar.style.width = (progress * 100) + "%";
                if (ptext) ptext.textContent = `Mesh refinement ${info.iteration}/${p.max_iters}: ` +
                                   `Energy err=${info.energy_error.toExponential(2)}, ` +
                                   `Grid=${info.nodes_x}x${info.nodes_y}`;
                log(`Pass ${info.iteration}: Energy error=${info.energy_error.toExponential(3)}, Param error=${info.param_error.toExponential(3)}, Grid=${info.nodes_x}x${info.nodes_y}`);
            },
            shouldStop: () => stopRequested
        });

        if (stopRequested) {
            log("Simulation stopped by user");
            pbar.style.width = "0%";
            return;
        }

        // Use the initial results as cache for frequency-dependent calculations
        const cachedResults = results;

        // If causal materials are enabled, recalculate max frequency with the model applied
        // (solve_adaptive used non-causal materials during mesh refinement)
        if (solver.use_causal_materials) {
            const maxFreqResult = await solver.computeAtFrequency(maxFreq, cachedResults);
            frequencySweepResults.push({ freq: maxFreq, result: maxFreqResult });
        } else {
            // Store the result from solve_adaptive (non-causal)
            frequencySweepResults.push({ freq: maxFreq, result: results });
        }

        // Redraw to show E-field overlay on geometry
        draw();
        drawResultsPlot();
        drawSParamPlot();

        // Check if interpolating sweep is enabled
        const fMax = Math.max(...frequencies);
        const nonZeroFreqs = frequencies.filter(f => f > 0);
        const fMinNonZero = nonZeroFreqs.length > 0 ? Math.min(...nonZeroFreqs) : 0;
        const useInterpolation = document.getElementById('chk_interp_sweep')?.checked
            && nonZeroFreqs.length > 1
            && fMax > fMinNonZero;

        if (useInterpolation) {
            // Interpolating sweep: adaptively sample RLGC, then interpolate
            const tolPercent = parseFloat(document.getElementById('interp_tolerance')?.value);
            if (isNaN(tolPercent) || tolPercent <= 0) {
                throw new Error("Interpolation tolerance must be a positive number.");
            }
            const tolerance = tolPercent / 100;

            // Handle DC point separately (can't use log-frequency axis)
            const hasDC = frequencies.includes(0);
            if (hasDC) {
                const dcResult = await solver.computeAtFrequency(0, cachedResults);
                frequencySweepResults.push({ freq: 0, result: dcResult });
            }

            log(`Interpolating sweep (tolerance ${tolPercent}%)...`);

            const sweep = new InterpolatingSweep(solver, cachedResults, { tolerance });
            const nSamples = await sweep.run(fMinNonZero, fMax, {
                onProgress: (info) => {
                    const progress = 0.5 + 0.5 * Math.min(info.iteration / 4, 0.9);
                    pbar.style.width = (progress * 100) + "%";
                    if (ptext) ptext.textContent = `Interpolating sweep: ${info.totalSamples} samples, ` +
                        `error=${(info.maxError * 100).toFixed(3)}%`;

                    // Update plots in real-time from current interpolation
                    if (info.iteration > 0) {
                        frequencySweepResults = hasDC
                            ? [{ freq: 0, result: frequencySweepResults.find(r => r.freq === 0).result },
                               ...sweep.buildResults(nonZeroFreqs)]
                            : sweep.buildResults(nonZeroFreqs);
                        frequencySweepResults.sort((a, b) => a.freq - b.freq);
                        drawResultsPlot();
                        drawSParamPlot();
                    }
                },
                shouldStop: () => stopRequested
            });

            if (!stopRequested) {
                // Build final results from converged interpolation
                const interpResults = sweep.buildResults(nonZeroFreqs);
                if (hasDC) {
                    const dcEntry = frequencySweepResults.find(r => r.freq === 0);
                    frequencySweepResults = [dcEntry, ...interpResults];
                } else {
                    frequencySweepResults = interpResults;
                }
                log(`Interpolating sweep: ${nSamples + (hasDC ? 1 : 0)} exact solves for ${frequencies.length} output points`);
            }
        } else {
            // Discrete sweep: compute at every frequency point
            log(`Calculating frequency sweep (${frequencies.length} points)...`);

            for (let i = 0; i < frequencies.length; i++) {
                const freq = frequencies[i];

                // Skip if this is the max frequency (already calculated above)
                if (freq === maxFreq) {
                    continue;
                }

                // Yield to event loop to allow UI updates
                await new Promise(resolve => setTimeout(resolve, 0));

                if (stopRequested) {
                    log("Simulation stopped by user");
                    break;
                }

                // Use optimized frequency sweep - only recalculates frequency-dependent losses
                // (or full solve if causal materials are enabled)
                const result = await solver.computeAtFrequency(freq, cachedResults);

                frequencySweepResults.push({ freq, result });

                // Update progress (second half is for frequency sweep)
                const progress = 0.5 + (i + 1) / frequencies.length * 0.5;
                pbar.style.width = (progress * 100) + "%";
                if (ptext) ptext.textContent = `Frequency sweep: ${i + 1}/${frequencies.length} (${(freq / 1e9).toFixed(2)} GHz)`;

                // Yield to event loop periodically and update plots in real time
                if (i % 10 === 0) {
                    await new Promise(resolve => setTimeout(resolve, 0));
                    frequencySweepResults.sort((a, b) => a.freq - b.freq);
                    drawResultsPlot();
                    drawSParamPlot();
                }
            }
        }

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
            log(`\nDIFFERENTIAL RESULTS:\n` +
                     `======================\n` +
                     `Differential Impedance Z_diff: ${results.Z_diff.toFixed(2)} Ohm  (2 x Z_odd)\n` +
                     `Common-Mode Impedance Z_common: ${results.Z_common.toFixed(2)} Ohm  (Z_even / 2)\n` +
                     `\nModal Impedances:\n` +
                     `  Odd-Mode  Z_odd:  ${odd.Z0.toFixed(2)} Ohm  (eps_eff = ${odd.eps_eff.toFixed(3)})\n` +
                     `  Even-Mode Z_even: ${even.Z0.toFixed(2)} Ohm  (eps_eff = ${even.eps_eff.toFixed(3)})\n` +
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
            log(`\nRESULTS:\n` +
                     `----------------------\n` +
                     `Z0: ${mode0.Z0.toFixed(2)} Ohm\n` +
                     `eps_eff: ${mode0.eps_eff.toFixed(3)}\n` +
                     `${lossStr}`);
        }

        // Update plots
        drawResultsPlot();
        drawSParamPlot();

        // Save geometry and frequency hash for change tracking
        lastSolvedGeometry = getGeometryHash();
        lastSolvedFrequency = getFrequencyHash();
        updateResultNotices();

    } catch (e) {
        console.error(e);
        log("Error: " + e.message);
    } finally {
        // Restore button to "Solve" mode
        btn.textContent = 'Solve';
        btn.classList.remove('stop-mode');
        pbar.style.width = '100%';
        if (ptext) ptext.style.display = 'none';
        stopRequested = false;
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
    const useSm       = document.getElementById('chk_solder_mask').checked;
    const useTopDiel  = document.getElementById('chk_top_diel').checked;
    const useGndCut   = document.getElementById('chk_gnd_cut').checked;
    const useEnclosure= document.getElementById('chk_enclosure').checked;
    const usePlating  = document.getElementById('chk_plating').checked;

    const groupEnabled = {
        always: true,
        diff: isDiff,
        gcpw: isGcpw,
        stripline: isStripline,
        sm: useSm,
        top_diel: useTopDiel,
        gnd_cut: useGndCut,
        enclosure: useEnclosure,
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
    // For unitless/fixedUnit params, min/maxSI are already display values
    const minDisplay = isUnitless ? minSI : convertToDisplayUnit(minSI, displayUnit);
    const maxDisplay = isUnitless ? maxSI : convertToDisplayUnit(maxSI, displayUnit);
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
    sweepStopRequested = false;
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

    log(`Parameter sweep: ${cfg.label} ${minDisplay}–${maxDisplay}${displayUnit ? ' ' + displayUnit : ''} (${nPoints} pts) @ ${(freqHz/1e9).toFixed(3)} GHz`);

    try {
        for (let i = 0; i < nPoints; i++) {
            if (sweepStopRequested) { log('Sweep stopped.'); break; }

            const displayVal = minDisplay + (maxDisplay - minDisplay) * i / (nPoints - 1);
            // Temporarily set value for updateGeometry to read, then restore
            inputEl.value = isUnitless ? displayVal : `${displayVal} ${displayUnit}`;
            updateGeometry();
            inputEl.value = originalValue;

            if (!solver) { log(`Point ${i+1}: solver init failed, skipping.`); continue; }

            solver.ensure_mesh();
            const cachedResults = await solver.solve_adaptive({
                max_iters: p.max_iters,
                energy_tol: p.tolerance,
                param_tol: 0.05,
                max_nodes: p.max_nodes * 1000,
                min_converged_passes: p.min_converged_passes,
                onProgress: () => {},
                shouldStop: () => sweepStopRequested
            });

            if (sweepStopRequested) { log('Sweep stopped.'); break; }

            const result = await solver.computeAtFrequency(freqHz, cachedResults);
            parameterSweepResults.push({ paramValue: displayVal, result });
            progressText.textContent = `${i + 1}/${nPoints}`;
            if (i === 0) updateSweepDiffCheckbox();

            redrawSweepPlot();
            await new Promise(r => setTimeout(r, 0)); // yield to UI
        }
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
        sweepStopRequested = false;
        progressText.textContent = '';
    }
}

function resizeCanvas() {
    const container = document.getElementById('sim_canvas');
    const Plotly = getPlotly();
    if (container && Plotly) {
        Plotly.Plots.resize(container);
    }
}

function bindEvents() {
    document.getElementById('btn_solve').onclick = () => {
        const btn = document.getElementById('btn_solve');
        if (btn.textContent === 'Stop') {
            // Stop the simulation
            stopRequested = true;
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

    // Parameter sweep events
    document.getElementById('btn-run-sweep').addEventListener('click', () => {
        if (isSweeping || isSimulating) { log('Cannot sweep while simulation is running.'); return; }
        runParameterSweep();
    });
    document.getElementById('btn-stop-sweep').addEventListener('click', () => {
        sweepStopRequested = true;
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
                plating: p.use_plating ? {
                    sigma: p.plating_sigma,
                    thickness: p.plating_t,
                    rq: p.plating_rq,
                    top: p.plating_top,
                    sides: p.plating_sides,
                    bottom: p.plating_bottom,
                    thick_corners: p.plating_thick_corners
                } : null,
                freqStart: frequencySweepResults[0].freq,
                freqStop: frequencySweepResults[frequencySweepResults.length - 1].freq,
                numPoints: frequencySweepResults.length
            };
            const filename = exportSnP(frequencySweepResults, length, Z_ref, isDifferential, params);
            log(`Exported ${filename}`);
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
        'inp_tolerance': { min: 0, max: 1, default: 0.05, label: 'Tolerance' },
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
            });
        }
    });

    // Real-time updates for checkboxes
    const geometryCheckboxes = [
        'chk_solder_mask', 'chk_top_diel', 'chk_gnd_cut', 'chk_enclosure', 'chk_side_gnd', 'chk_top_gnd',
        'chk_plating', 'chk_plating_top', 'chk_plating_sides', 'chk_plating_bottom', 'chk_plating_thick_corners'
    ];

    geometryCheckboxes.forEach(id => {
        const el = document.getElementById(id);
        if (el) {
            el.addEventListener('change', () => {
                updateGeometry();
                draw();
                updateResultNotices();
                updateSweepNotice();
            });
        }
    });

    // Transmission line type selector - reset zoom when type changes
    document.getElementById('tl_type').addEventListener('change', () => {
        updateGeometry();
        draw(true);  // Reset zoom/pan for new geometry
        updateResultNotices();
        updateSweepNotice();
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
