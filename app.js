import { MicrostripSolver } from './microstrip.js';
import { CONSTANTS } from './field_solver.js';
import { makeStreamlineTraceFromConductors } from './streamlines.js';
import { computeSParamsSingleEnded, computeSParamsDifferential, sParamTodB } from './sparameters.js';
import { exportSnP } from './snp_export.js';
import { draw, drawResultsPlot, drawSParamPlot, setGlobals, setCurrentView } from './plot.js';
const Plotly = window.Plotly;

let solver = null;
let stopRequested = false;
let frequencySweepResults = null;  // Array of {freq, result} objects
let currentTab = 'geometry';
let geometryChanged = false;  // Track if geometry has changed since last solve
let lastSolvedGeometry = null;  // Hash of geometry params from last solve
let lastSolvedFrequency = null;  // Frequency params from last solve

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
    return window.parseValueWithUnit ?
        window.parseValueWithUnit(element.value, defaultUnit) :
        parseFloat(element.value);
}

// --- URL Parameter Serialization ---

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
        er: parseFloat(document.getElementById('inp_er').value),
        tand: parseFloat(document.getElementById('inp_tand').value),
        sigma: parseFloat(document.getElementById('inp_sigma').value),
        freq_start: getDisplayValue('freq-start'),
        freq_stop: getDisplayValue('freq-stop'),
        freq_points: parseInt(document.getElementById('freq-points').value),
        trace_spacing: getDisplayValue('inp_trace_spacing'),
        gap: getDisplayValue('inp_gap'),
        top_gnd_w: getDisplayValue('inp_top_gnd_w'),
        via_gap: getDisplayValue('inp_via_gap'),
        stripline_top_h: getDisplayValue('inp_air_top'),
        er_top: parseFloat(document.getElementById('inp_er_top').value),
        tand_top: parseFloat(document.getElementById('inp_tand_top').value),
        use_sm: document.getElementById('chk_solder_mask').checked ? 1 : 0,
        sm_t_sub: getDisplayValue('inp_sm_t_sub'),
        sm_t_trace: getDisplayValue('inp_sm_t_trace'),
        sm_t_side: getDisplayValue('inp_sm_t_side'),
        sm_er: parseFloat(document.getElementById('inp_sm_er').value),
        sm_tand: parseFloat(document.getElementById('inp_sm_tand').value),
        use_top_diel: document.getElementById('chk_top_diel').checked ? 1 : 0,
        top_diel_h: getDisplayValue('inp_top_diel_h'),
        top_diel_er: parseFloat(document.getElementById('inp_top_diel_er').value),
        top_diel_tand: parseFloat(document.getElementById('inp_top_diel_tand').value),
        use_gnd_cut: document.getElementById('chk_gnd_cut').checked ? 1 : 0,
        gnd_cut_w: getDisplayValue('inp_gnd_cut_w'),
        gnd_cut_h: getDisplayValue('inp_gnd_cut_h'),
        use_enclosure: document.getElementById('chk_enclosure').checked ? 1 : 0,
        use_side_gnd: document.getElementById('chk_side_gnd').checked ? 1 : 0,
        use_top_gnd: document.getElementById('chk_top_gnd').checked ? 1 : 0,
        enclosure_width: getDisplayValue('inp_enclosure_width'),
        enclosure_height: getDisplayValue('inp_enclosure_height'),
        max_iters: parseInt(document.getElementById('inp_max_iters').value),
        tolerance: parseFloat(document.getElementById('inp_tolerance').value),
        max_nodes: parseInt(document.getElementById('inp_max_nodes').value),
        rq: getDisplayValue('inp_rq'),
        sparam_length: getDisplayValue('sparam-length'),
        sparam_z_ref: parseFloat(document.getElementById('sparam-z-ref').value),
    };
}

/**
 * Serialize settings to URL-safe base64 string
 */
function settingsToURL(settings) {
    const json = JSON.stringify(settings);
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
 */
function restoreSettings(settings) {
    if (!settings) return false;

    try {
        // Helper to restore value with unit
        const setValueWithUnit = (id, value) => {
            const element = document.getElementById(id);
            if (!element || value === undefined || value === null || isNaN(value)) return;
            const unit = window.getDefaultUnit ? window.getDefaultUnit(id) : '';
            if (unit && element.classList.contains('unit-input')) {
                element.value = `${value} ${unit}`;
            } else {
                element.value = value;
            }
        };

        // Set input values
        if (settings.tl_type) document.getElementById('tl_type').value = settings.tl_type;
        setValueWithUnit('inp_w', settings.w);
        setValueWithUnit('inp_h', settings.h);
        setValueWithUnit('inp_t', settings.t);
        if (settings.er !== undefined) document.getElementById('inp_er').value = settings.er;
        if (settings.tand !== undefined) document.getElementById('inp_tand').value = settings.tand;
        if (settings.sigma !== undefined) document.getElementById('inp_sigma').value = settings.sigma;
        setValueWithUnit('freq-start', settings.freq_start);
        setValueWithUnit('freq-stop', settings.freq_stop);
        if (settings.freq_points !== undefined) document.getElementById('freq-points').value = settings.freq_points;
        setValueWithUnit('inp_trace_spacing', settings.trace_spacing);
        setValueWithUnit('inp_gap', settings.gap);
        setValueWithUnit('inp_top_gnd_w', settings.top_gnd_w);
        setValueWithUnit('inp_via_gap', settings.via_gap);
        setValueWithUnit('inp_air_top', settings.stripline_top_h);
        if (settings.er_top !== undefined) document.getElementById('inp_er_top').value = settings.er_top;
        if (settings.tand_top !== undefined) document.getElementById('inp_tand_top').value = settings.tand_top;

        // Checkboxes
        if (settings.use_sm !== undefined) document.getElementById('chk_solder_mask').checked = !!settings.use_sm;
        setValueWithUnit('inp_sm_t_sub', settings.sm_t_sub);
        setValueWithUnit('inp_sm_t_trace', settings.sm_t_trace);
        setValueWithUnit('inp_sm_t_side', settings.sm_t_side);
        if (settings.sm_er !== undefined) document.getElementById('inp_sm_er').value = settings.sm_er;
        if (settings.sm_tand !== undefined) document.getElementById('inp_sm_tand').value = settings.sm_tand;

        if (settings.use_top_diel !== undefined) document.getElementById('chk_top_diel').checked = !!settings.use_top_diel;
        setValueWithUnit('inp_top_diel_h', settings.top_diel_h);
        if (settings.top_diel_er !== undefined) document.getElementById('inp_top_diel_er').value = settings.top_diel_er;
        if (settings.top_diel_tand !== undefined) document.getElementById('inp_top_diel_tand').value = settings.top_diel_tand;

        if (settings.use_gnd_cut !== undefined) document.getElementById('chk_gnd_cut').checked = !!settings.use_gnd_cut;
        setValueWithUnit('inp_gnd_cut_w', settings.gnd_cut_w);
        setValueWithUnit('inp_gnd_cut_h', settings.gnd_cut_h);

        if (settings.use_enclosure !== undefined) document.getElementById('chk_enclosure').checked = !!settings.use_enclosure;
        if (settings.use_side_gnd !== undefined) document.getElementById('chk_side_gnd').checked = !!settings.use_side_gnd;
        if (settings.use_top_gnd !== undefined) document.getElementById('chk_top_gnd').checked = !!settings.use_top_gnd;
        setValueWithUnit('inp_enclosure_width', settings.enclosure_width);
        setValueWithUnit('inp_enclosure_height', settings.enclosure_height);

        if (settings.max_iters !== undefined) document.getElementById('inp_max_iters').value = settings.max_iters;
        if (settings.tolerance !== undefined) document.getElementById('inp_tolerance').value = settings.tolerance;
        if (settings.max_nodes !== undefined) document.getElementById('inp_max_nodes').value = settings.max_nodes;
        setValueWithUnit('inp_rq', settings.rq);

        setValueWithUnit('sparam-length', settings.sparam_length);
        if (settings.sparam_z_ref !== undefined) document.getElementById('sparam-z-ref').value = settings.sparam_z_ref;

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

function shouldStop() {
    return stopRequested;
}

function log(msg) {
    const c = document.getElementById('console_out');
    c.textContent += msg + "\n";
    c.scrollTop = c.scrollHeight;
}

function nanToNull(input) {
    return Number.isNaN(input) ? null : input;
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
        top_gnd_w: p.top_gnd_w,
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
        rq: p.rq
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
    } else {
        const currentGeometry = getGeometryHash();
        const currentFrequency = getFrequencyHash();
        const geometryChanged = lastSolvedGeometry && currentGeometry !== lastSolvedGeometry;
        const frequencyChanged = lastSolvedFrequency && currentFrequency !== lastSolvedFrequency;

        if (geometryChanged) {
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
        } else if (frequencyChanged) {
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
        er: parseFloat(document.getElementById('inp_er').value),
        tand: parseFloat(document.getElementById('inp_tand').value),
        sigma: parseFloat(document.getElementById('inp_sigma').value),
        freq: getInputValue('freq-start'),
        nx: 30,  // Fixed initial grid size
        ny: 30,  // Fixed initial grid size
        // Differential parameters
        trace_spacing: getInputValue('inp_trace_spacing'),
        // GCPW specific parameters
        gap: getInputValue('inp_gap'),
        top_gnd_w: getInputValue('inp_top_gnd_w'),
        via_gap: getInputValue('inp_via_gap'),
        // Stripline parameters
        stripline_top_h: getInputValue('inp_air_top'),
        er_top: parseFloat(document.getElementById('inp_er_top').value),
        tand_top: parseFloat(document.getElementById('inp_tand_top').value),
        // Solder mask parameters
        use_sm: document.getElementById('chk_solder_mask').checked,
        sm_t_sub: getInputValue('inp_sm_t_sub'),
        sm_t_trace: getInputValue('inp_sm_t_trace'),
        sm_t_side: getInputValue('inp_sm_t_side'),
        sm_er: parseFloat(document.getElementById('inp_sm_er').value),
        sm_tand: parseFloat(document.getElementById('inp_sm_tand').value),
        // Top dielectric parameters
        use_top_diel: document.getElementById('chk_top_diel').checked,
        top_diel_h: getInputValue('inp_top_diel_h'),
        top_diel_er: parseFloat(document.getElementById('inp_top_diel_er').value),
        top_diel_tand: parseFloat(document.getElementById('inp_top_diel_tand').value),
        // Ground cutout parameters
        use_gnd_cut: document.getElementById('chk_gnd_cut').checked,
        gnd_cut_w: getInputValue('inp_gnd_cut_w'),
        gnd_cut_h: getInputValue('inp_gnd_cut_h'),
        // Enclosure parameters
        use_enclosure: document.getElementById('chk_enclosure').checked,
        use_side_gnd: document.getElementById('chk_side_gnd').checked,
        use_top_gnd: document.getElementById('chk_top_gnd').checked,
        enclosure_width: nanToNull(getInputValue('inp_enclosure_width')),
        enclosure_height: nanToNull(getInputValue('inp_enclosure_height')),
        max_iters: parseInt(document.getElementById('inp_max_iters').value),
        tolerance: parseFloat(document.getElementById('inp_tolerance').value),
        max_nodes: parseInt(document.getElementById('inp_max_nodes').value),
        // Surface roughness parameter
        rq: getInputValue('inp_rq'),
    };
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
                top_gnd_width: p.top_gnd_w,
                via_gap: p.via_gap,
                use_vias: true,
                // Surface roughness
                rq: p.rq,
            };
            // Add solder mask options
            if (p.use_sm) {
                options.use_sm = true;
                options.sm_t_sub = p.sm_t_sub;
                options.sm_t_trace = p.sm_t_trace;
                options.sm_t_side = p.sm_t_side;
                options.sm_er = p.sm_er;
                options.sm_tand = p.sm_tand;
            }
            // Add top dielectric options
            if (p.use_top_diel) {
                options.top_diel_h = p.top_diel_h;
                options.top_diel_er = p.top_diel_er;
                options.top_diel_tand = p.top_diel_tand;
            }
            // Add ground cutout options
            if (p.use_gnd_cut) {
                options.gnd_cut_width = p.gnd_cut_w;
                options.gnd_cut_sub_h = p.gnd_cut_h;
            }
            // Add enclosure options
            if (p.use_enclosure) {
                options.enclosure_width = p.enclosure_width;
                options.enclosure_height = p.enclosure_height;
                // Set boundaries based on side and top ground options
                const left_bc = p.use_side_gnd ? "gnd" : "open";
                const right_bc = p.use_side_gnd ? "gnd" : "open";
                const top_bc = p.use_top_gnd ? "gnd" : "open";
                options.boundaries = [left_bc, right_bc, top_bc, "gnd"];
            }
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
                top_gnd_width: p.top_gnd_w,
                via_gap: p.via_gap,
                use_vias: true,
                // Surface roughness
                rq: p.rq,
            };
            // Add solder mask options
            if (p.use_sm) {
                options.use_sm = true;
                options.sm_t_sub = p.sm_t_sub;
                options.sm_t_trace = p.sm_t_trace;
                options.sm_t_side = p.sm_t_side;
                options.sm_er = p.sm_er;
                options.sm_tand = p.sm_tand;
            }
            // Add top dielectric options
            if (p.use_top_diel) {
                options.top_diel_h = p.top_diel_h;
                options.top_diel_er = p.top_diel_er;
                options.top_diel_tand = p.top_diel_tand;
            }
            // Add ground cutout options
            if (p.use_gnd_cut) {
                options.gnd_cut_width = p.gnd_cut_w;
                options.gnd_cut_sub_h = p.gnd_cut_h;
            }
            // Add enclosure options
            if (p.use_enclosure) {
                options.enclosure_width = p.enclosure_width;
                options.enclosure_height = p.enclosure_height;
                // Set boundaries based on side and top ground options
                const left_bc = p.use_side_gnd ? "gnd" : "open";
                const right_bc = p.use_side_gnd ? "gnd" : "open";
                const top_bc = p.use_top_gnd ? "gnd" : "open";
                options.boundaries = [left_bc, right_bc, top_bc, "gnd"];
            }
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
            // Add solder mask options
            if (p.use_sm) {
                options.use_sm = true;
                options.sm_t_sub = p.sm_t_sub;
                options.sm_t_trace = p.sm_t_trace;
                options.sm_t_side = p.sm_t_side;
                options.sm_er = p.sm_er;
                options.sm_tand = p.sm_tand;
            }
            // Add top dielectric options
            if (p.use_top_diel) {
                options.top_diel_h = p.top_diel_h;
                options.top_diel_er = p.top_diel_er;
                options.top_diel_tand = p.top_diel_tand;
            }
            // Add ground cutout options
            if (p.use_gnd_cut) {
                options.gnd_cut_width = p.gnd_cut_w;
                options.gnd_cut_sub_h = p.gnd_cut_h;
            }
            // Add enclosure options
            if (p.use_enclosure) {
                options.enclosure_width = p.enclosure_width;
                options.enclosure_height = p.enclosure_height;
                // Set boundaries based on side and top ground options
                const left_bc = p.use_side_gnd ? "gnd" : "open";
                const right_bc = p.use_side_gnd ? "gnd" : "open";
                const top_bc = p.use_top_gnd ? "gnd" : "open";
                options.boundaries = [left_bc, right_bc, top_bc, "gnd"];
            }
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
            // Add solder mask options
            if (p.use_sm) {
                options.use_sm = true;
                options.sm_t_sub = p.sm_t_sub;
                options.sm_t_trace = p.sm_t_trace;
                options.sm_t_side = p.sm_t_side;
                options.sm_er = p.sm_er;
                options.sm_tand = p.sm_tand;
            }
            // Add top dielectric options
            if (p.use_top_diel) {
                options.top_diel_h = p.top_diel_h;
                options.top_diel_er = p.top_diel_er;
                options.top_diel_tand = p.top_diel_tand;
            }
            // Add ground cutout options
            if (p.use_gnd_cut) {
                options.gnd_cut_width = p.gnd_cut_w;
                options.gnd_cut_sub_h = p.gnd_cut_h;
            }
            // Add enclosure options (stripline already has top ground)
            if (p.use_enclosure) {
                options.enclosure_width = p.enclosure_width;
                // Set side boundaries based on side ground option
                const left_bc = p.use_side_gnd ? "gnd" : "open";
                const right_bc = p.use_side_gnd ? "gnd" : "open";
                options.boundaries = [left_bc, right_bc, "gnd", "gnd"];
            }
            solver = new MicrostripSolver(options);
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
            // Enclosure options (stripline already has top ground)
            if (p.use_enclosure) {
                options.enclosure_width = p.enclosure_width;
                // Set side boundaries based on side ground option
                const left_bc = p.use_side_gnd ? "gnd" : "open";
                const right_bc = p.use_side_gnd ? "gnd" : "open";
                options.boundaries = [left_bc, right_bc, "gnd", "gnd"];
            }
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

            // Solder mask
            if (p.use_sm) {
                options.use_sm = true;
                options.sm_t_sub = p.sm_t_sub;
                options.sm_t_trace = p.sm_t_trace;
                options.sm_t_side = p.sm_t_side;
                options.sm_er = p.sm_er;
                options.sm_tand = p.sm_tand;
            }

            // Top dielectric (embedded microstrip)
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

            // Enclosure options
            if (p.use_enclosure) {
                options.enclosure_width = p.enclosure_width;
                options.enclosure_height = p.enclosure_height;
                // Set boundaries based on side and top ground options
                const left_bc = p.use_side_gnd ? "gnd" : "open";
                const right_bc = p.use_side_gnd ? "gnd" : "open";
                const top_bc = p.use_top_gnd ? "gnd" : "open";
                options.boundaries = [left_bc, right_bc, top_bc, "gnd"];
            }

            solver = new MicrostripSolver(options);
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
    pbar.style.width = '0%';
    if (ptext) ptext.style.display = 'block';
    log("Starting simulation...");

    try {
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
        log(`Running adaptive analysis (max ${p.max_iters} iterations, max ${p.max_nodes} nodes, tolerance ${p.tolerance})...`);

        let results = await solver.solve_adaptive({
            max_iters: p.max_iters,
            energy_tol: p.tolerance,
            param_tol: 0.05,
            max_nodes: p.max_nodes,
            onProgress: (info) => {
                const progress = info.iteration / p.max_iters * 0.5;  // First half is for mesh refinement
                pbar.style.width = (progress * 100) + "%";
                if (ptext) ptext.textContent = `Mesh refinement ${info.iteration}/${p.max_iters}: ` +
                                   `Energy err=${info.energy_error.toExponential(2)}, ` +
                                   `Grid=${info.nodes_x}x${info.nodes_y}`;
                log(`Pass ${info.iteration}: Energy error=${info.energy_error.toExponential(3)}, Param error=${info.param_error.toExponential(3)}, Grid=${info.nodes_x}x${info.nodes_y}`);
            },
            shouldStop: shouldStop
        });

        if (stopRequested) {
            log("Simulation stopped by user");
            pbar.style.width = "0%";
            return;
        }

        // Store the first result (highest frequency)
        frequencySweepResults.push({ freq: maxFreq, result: results });

        // Redraw to show E-field overlay on geometry
        draw();

        // Now run frequency sweep with cached fields (fast path)
        // The potential distribution and electric fields don't change with frequency,
        // only the losses depend on frequency, so we can reuse the cached results.
        log(`Calculating frequency sweep (${frequencies.length} points)...`);

        // Use the initial results as cache for frequency-dependent calculations
        const cachedResults = results;

        for (let i = 0; i < frequencies.length; i++) {
            const freq = frequencies[i];

            // Skip if this is the max frequency (already calculated)
            if (freq === maxFreq) {
                continue;
            }

            if (stopRequested) {
                log("Simulation stopped by user");
                break;
            }

            // Use optimized frequency sweep - only recalculates frequency-dependent losses
            const result = solver.computeAtFrequency(freq, cachedResults);

            frequencySweepResults.push({ freq, result });

            // Update progress (second half is for frequency sweep)
            const progress = 0.5 + (i + 1) / frequencies.length * 0.5;
            pbar.style.width = (progress * 100) + "%";
            if (ptext) ptext.textContent = `Frequency sweep: ${i + 1}/${frequencies.length} (${(freq / 1e9).toFixed(2)} GHz)`;

            // Yield to event loop periodically to prevent UI freeze
            if (i % 10 === 0) {
                await new Promise(resolve => setTimeout(resolve, 0));
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
    }
}

function resizeCanvas() {
    const container = document.getElementById('sim_canvas');
    if (container) {
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

    // Log checkbox for results plot
    const resultsLogX = document.getElementById('results-log-x');
    if (resultsLogX) {
        resultsLogX.addEventListener('change', () => {
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

    // Export SnP button
    const exportSnpBtn = document.getElementById('export-snp');
    if (exportSnpBtn) {
        exportSnpBtn.addEventListener('click', () => {
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
            const Z_ref = parseFloat(document.getElementById('sparam-z-ref').value);
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
                freqStart: frequencySweepResults[0].freq,
                freqStop: frequencySweepResults[frequencySweepResults.length - 1].freq,
                numPoints: frequencySweepResults.length
            };
            const filename = exportSnP(frequencySweepResults, length, Z_ref, isDifferential, params);
            log(`Exported ${filename}`);
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
        'inp_max_nodes': { min: 1000, default: 20000, integer: true, label: 'Max nodes' },
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
        'inp_gap', 'inp_top_gnd_w', 'inp_via_gap',
        'inp_air_top', 'inp_er_top', 'inp_tand_top',
        'inp_sm_t_sub', 'inp_sm_t_trace', 'inp_sm_t_side', 'inp_sm_er', 'inp_sm_tand',
        'inp_top_diel_h', 'inp_top_diel_er', 'inp_top_diel_tand',
        'inp_gnd_cut_w', 'inp_gnd_cut_h',
        'inp_enclosure_width', 'inp_enclosure_height',
        'freq-start'
    ];

    geometryInputs.forEach(id => {
        const el = document.getElementById(id);
        if (el) {
            el.addEventListener('input', () => {
                updateGeometry();
                draw();
                updateResultNotices();
            });
        }
    });

    // Real-time updates for checkboxes
    const geometryCheckboxes = [
        'chk_solder_mask', 'chk_top_diel', 'chk_gnd_cut', 'chk_enclosure', 'chk_side_gnd', 'chk_top_gnd'
    ];

    geometryCheckboxes.forEach(id => {
        const el = document.getElementById(id);
        if (el) {
            el.addEventListener('change', () => {
                updateGeometry();
                draw();
                updateResultNotices();
            });
        }
    });

    // Transmission line type selector - reset zoom when type changes
    document.getElementById('tl_type').addEventListener('change', () => {
        updateGeometry();
        draw(true);  // Reset zoom/pan for new geometry
        updateResultNotices();
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
}


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
    ['chk_solder_mask', 'chk_top_diel', 'chk_gnd_cut', 'chk_enclosure'].forEach(id => {
        const checkbox = document.getElementById(id);
        if (checkbox) {
            const sectionId = id.replace('chk_', '') + '-params';
            const section = document.getElementById(sectionId);
            if (section) {
                section.style.display = checkbox.checked ? 'block' : 'none';
            }
        }
    });

    updateGeometry();
    draw();
    resizeCanvas();
    window.addEventListener('resize', resizeCanvas);
    log("Ready. Click 'Solve' to start simulation.");
}

// Start when DOM is ready
window.addEventListener('DOMContentLoaded', init);
