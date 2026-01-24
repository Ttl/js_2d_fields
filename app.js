import { MicrostripSolver } from './microstrip.js';
import { CONSTANTS } from './field_solver.js';
import { makeStreamlineTraceFromConductors } from './streamlines.js';
import { computeSParamsSingleEnded, computeSParamsDifferential, sParamTodB } from './sparameters.js';
import { exportSnP } from './snp_export.js';
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
        console.error('Failed to parse URL parameters:', e);
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
            console.log('Settings restored from URL');
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
    currentView = "geometry";

    // Show/hide mode selector based on differential transmission line
    const isDiff = p.tl_type.startsWith('diff_');
    const modeGroup = document.getElementById('plot-mode-group');
    if (modeGroup) {
        modeGroup.style.display = isDiff ? 'block' : 'none';
    }

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
                enclosure_height: p.stripline_top_h + p.t,
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
                enclosure_height: p.stripline_top_h + p.t,
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
    ptext.style.display = 'block';
    log("Starting simulation...");

    try {
        // Clear previous results
        frequencySweepResults = [];

        // Use the highest frequency for mesh generation (skin depth calculation)
        const maxFreq = Math.max(...frequencies);
        solver.freq = maxFreq;
        solver.omega = 2 * Math.PI * maxFreq;

        // Ensure mesh is generated before solving
        log("Calculating mesh...");
        solver.ensure_mesh();
        log("Mesh generated: " + solver.x.length + "x" + solver.y.length);

        // Run adaptive refinement at highest frequency first
        log(`Running adaptive analysis (max ${p.max_iters} iterations, max ${p.max_nodes} nodes, tolerance ${p.tolerance})...`);

        let results = await solver.solve_adaptive({
            max_iters: p.max_iters,
            tolerance: p.tolerance,
            param_tol: 0.001,
            max_nodes: p.max_nodes,
            onProgress: (info) => {
                const progress = info.iteration / p.max_iters * 0.5;  // First half is for mesh refinement
                pbar.style.width = (progress * 100) + "%";
                ptext.textContent = `Mesh refinement ${info.iteration}/${p.max_iters}: ` +
                                   `Energy err=${info.energy_error.toExponential(2)}, ` +
                                   `Grid=${info.nodes_x}x${info.nodes_y}`;
                log(`Pass ${info.iteration}: Energy error=${info.energy_error.toExponential(3)}, Param error=${info.param_error.toExponential(3)}, Grid=${info.nodes_x}x${info.nodes_y}`);
            },
            shouldStop: shouldStop
        });

        if (stopRequested) {
            log("Simulation stopped by user");
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
            ptext.textContent = `Frequency sweep: ${i + 1}/${frequencies.length} (${(freq / 1e9).toFixed(2)} GHz)`;

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
        // Restore button to "Solve Physics" mode
        btn.textContent = 'Solve Physics';
        btn.classList.remove('stop-mode');
        pbar.style.width = '100%';
        ptext.style.display = 'none';
        stopRequested = false;
    }
}

let showMesh = false;
let currentView = "geometry";

// Helper function to check if solver is in differential mode
function isDifferentialMode() {
    if (!solver || !solver.Ex || !solver.Ey) return false;
    // In differential mode, Ex and Ey are arrays of 2 arrays (odd and even modes)
    // Check if Ex[0] and Ex[1] are both arrays
    return Array.isArray(solver.Ex) &&
           solver.Ex.length === 2 &&
           Array.isArray(solver.Ex[0]) &&
           Array.isArray(solver.Ex[1]);
}

// Get the selected mode index from sidebar (0=odd, 1=even)
function getSelectedModeIndex() {
    const modeSelect = document.getElementById('plot-mode');
    return modeSelect && modeSelect.value === 'even' ? 1 : 0;
}

// Helper function to get Ex/Ey fields (handles differential mode)
function getFields() {
    if (!solver || !solver.Ex || !solver.Ey) {
        return { Ex: null, Ey: null };
    }

    if (isDifferentialMode()) {
        const modeIndex = getSelectedModeIndex();
        return { Ex: solver.Ex[modeIndex], Ey: solver.Ey[modeIndex] };
    } else {
        // Single-ended mode
        return { Ex: solver.Ex[0], Ey: solver.Ey[0] };
    }
}

// Helper function to get voltage potential (handles differential mode)
function getPotential() {
    if (!solver || !solver.V) {
        return null;
    }

    if (isDifferentialMode()) {
        const modeIndex = getSelectedModeIndex();
        return solver.V[modeIndex];
    } else {
        // Single-ended mode
        return solver.V[0];
    }
}

// Get plot options from sidebar
function getPlotOptions() {
    const streamlinesEl = document.getElementById('plot-streamlines');
    const contoursEl = document.getElementById('plot-contours');

    const streamlinesVal = streamlinesEl ? streamlinesEl.value.trim() : '';
    const contoursVal = contoursEl ? contoursEl.value.trim() : '';

    return {
        streamlines: streamlinesVal === '' ? 0 : parseInt(streamlinesVal) || 0,
        contours: contoursVal === '' ? 0 : parseInt(contoursVal) || 0
    };
}


// Interpolation functions for higher-resolution plots

// Find the index i such that arr[i] <= val < arr[i+1]
function find_idx(arr, val) {
    // A linear scan is used here. For very large grids, binary search could be an optimization.
    for (let i = 0; i < arr.length - 1; i++) {
        if (arr[i] <= val && val < arr[i + 1]) {
            return i;
        }
    }
    // Handle the edge case where val is exactly the last element
    if (val === arr[arr.length - 1]) {
        return arr.length - 2;
    }
    return -1;
}

// Bilinear interpolation for a point (x, y) within a grid cell
function bilinearInterpolate(x, y, x1, y1, x2, y2, q11, q12, q21, q22) {
    const denom = (x2 - x1) * (y2 - y1);
    if (denom === 0) {
        // Avoid division by zero if the grid cell has no area
        return q11;
    }
    const w11 = (x2 - x) * (y2 - y);
    const w12 = (x2 - x) * (y - y1);
    const w21 = (x - x1) * (y2 - y);
    const w22 = (x - x1) * (y - y1);
    return (w11 * q11 + w12 * q12 + w21 * q21 + w22 * q22) / denom;
}

// Interpolates data from an old grid (z_old) to a new, finer grid (x_new, y_new)
function interpolateGrid(x_old, y_old, z_old, x_new, y_new) {
    const z_new = Array(y_new.length).fill(0).map(() => Array(x_new.length).fill(0));

    for (let j = 0; j < y_new.length; j++) {
        const y_val = y_new[j];
        const y_idx1 = find_idx(y_old, y_val);
        if (y_idx1 === -1) continue;
        const y_idx2 = y_idx1 + 1;

        for (let i = 0; i < x_new.length; i++) {
            const x_val = x_new[i];
            const x_idx1 = find_idx(x_old, x_val);
            if (x_idx1 === -1) continue;
            const x_idx2 = x_idx1 + 1;

            // Values at the four corners of the cell in the old grid
            const q11 = z_old[y_idx1][x_idx1];
            const q12 = z_old[y_idx2][x_idx1];
            const q21 = z_old[y_idx1][x_idx2];
            const q22 = z_old[y_idx2][x_idx2];

            z_new[j][i] = bilinearInterpolate(
                x_val, y_val,
                x_old[x_idx1], y_old[y_idx1], x_old[x_idx2], y_old[y_idx2],
                q11, q12, q21, q22
            );
        }
    }
    return z_new;
}


function draw(resetZoom = false) {
    if (!solver) return;

    const container = document.getElementById('sim_canvas');
    const plotOptions = getPlotOptions();

    // Preserve current view state if plot exists (unless resetZoom is requested)
    let currentXRange = null;
    let currentYRange = null;
    if (!resetZoom && container && container.layout && container.layout.xaxis) {
        currentXRange = container.layout.xaxis.range;
        currentYRange = container.layout.yaxis.range;
    }

    let zData = [];
    let title = "";
    let colorscale = "Viridis";
    let zTitle = "";
    let shapes = [];
    let xMM, yMM, nx, ny, nyDisplay;

    // View selection
    if (currentView === "geometry") {
        title = "Transmission Line Geometry";

        // Determine display bounds using actual domain extent
        const maxY = solver.dielectrics.reduce((max, d) => Math.max(max, d.y_max), 0);

        // Draw dielectrics as rectangles (color by epsilon_r)
        for (const diel of solver.dielectrics) {
            if (diel.y_min > maxY) continue;

            const yMax = Math.min(diel.y_max, maxY);
            const er = diel.epsilon_r;

            // Color mapping: air (1.0) = white, higher er = green shades
            let fillcolor;
            if (er <= 1.01) {
                fillcolor = 'rgba(255, 255, 255, 0.8)';
            } else {
                // Green shades for dielectrics
                const intensity = Math.min(255, 100 + (er - 1) * 30);
                fillcolor = `rgba(100, ${intensity}, 100, 0.8)`;
            }

            shapes.push({
                type: 'rect',
                x0: diel.x_min * 1000,
                y0: diel.y_min * 1000,
                x1: diel.x_max * 1000,
                y1: yMax * 1000,
                fillcolor: fillcolor,
                line: { color: 'rgba(128, 128, 128, 0.3)', width: 0.5 },
                layer: 'below'
            });
        }

        // Draw conductors as rectangles (signal = orange, ground = dark gray)
        for (const cond of solver.conductors) {
            if (cond.y_min > maxY) continue;

            const yMax = Math.min(cond.y_max, maxY);
            const fillcolor = cond.is_signal ?
                'rgba(217, 119, 6, 1.0)' :  // Orange for signal
                'rgba(51, 51, 51, 1.0)';     // Dark gray for ground

            shapes.push({
                type: 'rect',
                x0: cond.x_min * 1000,
                y0: cond.y_min * 1000,
                x1: cond.x_max * 1000,
                y1: yMax * 1000,
                fillcolor: fillcolor,
                line: { color: 'rgba(0, 0, 0, 0.5)', width: 1 },
                layer: 'above'
            });
        }

        // If solution available, overlay E-field contours
        if (solver.solution_valid && solver.mesh_generated) {
            nx = solver.x.length;
            ny = solver.y.length;

            // Limit display Y
            const yArr = Array.from(solver.y);
            const maxYIdx = yArr.findIndex(y => y > maxY);
            nyDisplay = maxYIdx > 0 ? maxYIdx : ny;

            xMM = Array.from(solver.x, v => v * 1000);
            yMM = yArr.slice(0, nyDisplay).map(v => v * 1000);

            // Compute E-field magnitude
            const { Ex, Ey } = getFields();
            if (Ex && Ey && Ex.length >= nyDisplay) {
                for (let i = 0; i < nyDisplay; i++) {
                    const row = [];
                    if (Ex[i] && Ey[i]) {
                        for (let j = 0; j < nx; j++) {
                            row.push(Math.hypot(Ex[i][j], Ey[i][j]));
                        }
                    }
                    zData.push(row);
                }
            }

            // In geometry view with differential solver, always show odd mode
            const modeLabel = isDifferentialMode() ? " (Odd Mode)" : "";
            title = `Transmission Line Geometry with E-field${modeLabel}`;
        } else {
            // No solution - just axis scaling
            xMM = [0, solver.w * 2000];
            yMM = [0, maxY * 1000];
        }
    }

    else if ((currentView === "potential" || currentView === "potential_odd" || currentView === "potential_even") && solver.solution_valid) {
        // Ensure mesh exists for field visualization
        if (!solver.mesh_generated) {
            solver.ensure_mesh();
        }

        nx = solver.x.length;
        ny = solver.y.length;

        // Limit display Y to domain extent
        const yArr = Array.from(solver.y);
        const maxY = yArr[ny - 1];
        const maxYIdx = yArr.findIndex(y => y > maxY);
        nyDisplay = maxYIdx > 0 ? maxYIdx : ny;

        xMM = Array.from(solver.x, v => v * 1000);
        yMM = yArr.slice(0, nyDisplay).map(v => v * 1000);

        let modeLabel = "";
        if (currentView === "potential_odd") {
            modeLabel = " (Odd Mode)";
        } else if (currentView === "potential_even") {
            modeLabel = " (Even Mode)";
        }
        title = `Electric Potential${modeLabel} (V)`;
        zTitle = "Volts";

        const V = getPotential();
        if (V && V.length >= nyDisplay) {
            for (let i = 0; i < nyDisplay; i++) {
                zData.push(V[i].slice(0, nx));
            }
        }
    }

    else if ((currentView === "efield" || currentView === "efield_odd" || currentView === "efield_even") && solver.solution_valid) {
        // Ensure mesh exists for field visualization
        if (!solver.mesh_generated) {
            solver.ensure_mesh();
        }

        nx = solver.x.length;
        ny = solver.y.length;

        // Limit display Y to actual domain extent
        const yArr = Array.from(solver.y);
        const maxY = yArr[ny - 1];
        const maxYIdx = yArr.findIndex(y => y > maxY);
        nyDisplay = maxYIdx > 0 ? maxYIdx : ny;

        xMM = Array.from(solver.x, v => v * 1000);
        yMM = yArr.slice(0, nyDisplay).map(v => v * 1000);

        let modeLabel = "";
        if (currentView === "efield_odd") {
            modeLabel = " (Odd Mode)";
        } else if (currentView === "efield_even") {
            modeLabel = " (Even Mode)";
        }
        title = `|E| Field Magnitude${modeLabel} (V/m)`;
        zTitle = "V/m";

        const { Ex, Ey } = getFields();
        if (Ex && Ey && Ex.length >= nyDisplay) {
            for (let i = 0; i < nyDisplay; i++) {
                const row = [];
                if (Ex[i] && Ey[i]) {
                    for (let j = 0; j < nx; j++) {
                        row.push(Math.hypot(Ex[i][j], Ey[i][j]));
                    }
                }
                zData.push(row);
            }
        }
    }

    else {
        title = "No Data Available";
        // Create minimal dummy data
        xMM = [0, (solver.w || 1) * 2000];
        yMM = [0, (solver.h || 1) * 1000];
    }

    // Interpolate for smoother plots
    const INTERP_ENABLED = true; // Control flag for interpolation
    const INTERP_FACTOR = 4; // Interpolation multiplier

    // Save original mesh coordinates for mesh overlay before interpolation
    let xMM_mesh = xMM;
    let yMM_mesh = yMM;
    let nx_mesh = nx;
    let nyDisplay_mesh = nyDisplay;

    if (INTERP_ENABLED && zData.length > 0 && (currentView.includes("potential") || currentView.includes("efield"))) {

        const x_old = Array.from(solver.x); // Original grid X (meters)
        const y_old = Array.from(solver.y).slice(0, nyDisplay); // Original grid Y (meters)

        if (x_old.length > 1 && y_old.length > 1 && zData.length > 1 && zData[0].length > 1) {
            const nx_interp = (nx - 1) * INTERP_FACTOR + 1;
            const ny_interp = (nyDisplay - 1) * INTERP_FACTOR + 1;

            // Create a new, uniformly spaced, finer grid for interpolation
            const x_new = new Float64Array(nx_interp);
            const y_new = new Float64Array(ny_interp);

            const x_min = x_old[0];
            const x_max = x_old[x_old.length - 1];
            const y_min = y_old[0];
            const y_max = y_old[y_old.length - 1];

            for (let i = 0; i < nx_interp; i++) {
                x_new[i] = x_min + (x_max - x_min) * i / (nx_interp - 1);
            }
            for (let j = 0; j < ny_interp; j++) {
                y_new[j] = y_min + (y_max - y_min) * j / (ny_interp - 1);
            }

            // Perform interpolation from the original zData to the new grid
            const z_interp = interpolateGrid(x_old, y_old, zData, x_new, y_new);

            // Update plot variables with the new, higher-resolution data
            xMM = Array.from(x_new, v => v * 1000);
            yMM = Array.from(y_new, v => v * 1000);
            zData = z_interp;
        }
    }


    // Main field trace
    let traces = [];

    if (currentView === "geometry" && zData.length > 0) {
        const { Ex, Ey } = getFields();

        //const fieldTraces = [];

        //const stepX = Math.max(1, Math.floor(nx / 30));        // density control
        //const stepY = Math.max(1, Math.floor(nyDisplay / 30));
        //const scale = 0.8; // arrow length in mm


        //if (Ex && Ey && Ex.length >= nyDisplay) {
        //    const xLines = [];
        //    const yLines = [];

        //    for (let i = 0; i < nyDisplay; i += stepY) {
        //        for (let j = 0; j < nx; j += stepX) {
        //            const ex = Ex[i]?.[j];
        //            const ey = Ey[i]?.[j];
        //            if (!ex || !ey) continue;

        //            const mag = Math.hypot(ex, ey);
        //            if (mag === 0) continue;

        //            // Base point (mm)
        //            const x0 = xMM[j];
        //            const y0 = yMM[i];

        //            // Direction (normalized)
        //            const dx = (ex / mag) * scale;
        //            const dy = (ey / mag) * scale;

        //            // Line segment
        //            xLines.push(x0, x0 + dx, null);
        //            yLines.push(y0, y0 + dy, null);
        //        }
        //    }

        //    traces.push({
        //        type: "scatter",
        //        mode: "lines",
        //        x: xLines,
        //        y: yLines,
        //        line: {
        //            width: 1.2,
        //            color: "black"
        //        },
        //        hoverinfo: "skip",
        //        name: "E-field lines"
        //    });
        //}
        //
        //traces.push({
        //    type: "contour",
        //    x: xMM,
        //    y: yMM,
        //    z: zData,
        //    colorscale: "Hot",
        //    opacity: 0.6,
        //    contours: {
        //        showlines: true,
        //        coloring: "heatmap",
        //        ncontours: 15
        //    },
        //    colorbar: {
        //        title: "|E| (V/m)",
        //        len: 0.6
        //    },
        //    hovertemplate:
        //        "x: %{x:.2f} mm<br>" +
        //        "y: %{y:.2f} mm<br>" +
        //        "|E|: %{z:.3e} V/m<extra></extra>"
        //});

        // Add streamlines if requested via plot options
        if (plotOptions.streamlines > 0) {
            traces.push(
                makeStreamlineTraceFromConductors(
                    Ex,
                    Ey,
                    solver.x,
                    solver.y,
                    solver.conductors,
                    plotOptions.streamlines
                )
            );
        }

    } else if (currentView === "geometry") {
        // Geometry only. Invisible scatter for axis scaling
        traces.push({
            type: "scatter",
            x: xMM,
            y: yMM,
            mode: "markers",
            marker: { size: 0, opacity: 0 },
            showlegend: false,
            hoverinfo: "skip"
        });
    } else if (zData.length > 0) {
        // Field views only. Use heatmap with optional contour lines
        const contourSettings = {
            coloring: 'heatmap',
            showlines: plotOptions.contours > 0,
        };
        if (plotOptions.contours > 0) {
            contourSettings.size = 0;  // Auto-calculate based on data
            contourSettings.ncontours = plotOptions.contours;
        }
        traces.push({
            type: "contour",
            x: xMM,
            y: yMM,
            z: zData,
            colorscale: colorscale,
            contours: contourSettings,
            line: {
              smoothing: 1.3,
              width: 0.5
            },
            colorbar: {
                title: zTitle,
                len: 0.8
            },
            hovertemplate:
                "x: %{x:.2f} mm<br>" +
                "y: %{y:.2f} mm<br>" +
                "value: %{z:.3e}<extra></extra>"
        });
    }

    // Mesh overlay
    if (showMesh && solver.solution_valid) {
        const stepX = 1;
        const stepY = 1;

        // Use original mesh coordinates (before interpolation)
        for (let j = 0; j < nx_mesh; j += stepX) {
            traces.push({
                type: "scatter",
                x: [xMM_mesh[j], xMM_mesh[j]],
                y: [yMM_mesh[0], yMM_mesh[nyDisplay_mesh - 1]],
                mode: "lines",
                line: { width: 0.2, color: "black" },
                showlegend: false,
                hoverinfo: "skip"
            });
        }

        for (let i = 0; i < nyDisplay_mesh; i += stepY) {
            traces.push({
                type: "scatter",
                x: [xMM_mesh[0], xMM_mesh[nx_mesh - 1]],
                y: [yMM_mesh[i], yMM_mesh[i]],
                mode: "lines",
                line: { width: 0.2, color: "black" },
                showlegend: false,
                hoverinfo: "skip"
            });
        }
    }

    // UI menues
    const layout = {
        title: title,
        xaxis: {
            title: "Width (mm)",
            scaleanchor: "y",
            scaleratio: 1,
            range: currentXRange  // Preserve zoom/pan
        },
        yaxis: {
            title: "Height (mm)",
            range: currentYRange  // Preserve zoom/pan
        },
        margin: { l: 70, r: 90, t: 50, b: 60 },
        hovermode: "closest",
        dragmode: "pan",
        plot_bgcolor: "#f8f9fa",
        shapes: shapes,  // Add vector shapes for geometry

        updatemenus: [
            {
                x: 0.01,
                y: 1.15,
                showactive: true,
                active: (() => {
                    if (currentView === "geometry") return 0;
                    if (currentView === "potential") return 1;
                    if (currentView === "efield") return 2;
                    return 0;
                })(),
                buttons: (() => {
                    const buttons = [
                        {
                            label: "Geometry",
                            method: "skip",
                            args: []
                        }
                    ];

                    if (solver.solution_valid) {
                        buttons.push({
                            label: "Potential",
                            method: "skip",
                            args: []
                        });
                        buttons.push({
                            label: "|E| Field",
                            method: "skip",
                            args: []
                        });
                    }

                    return buttons;
                })()
            }
        ]
    };

    const config = {
        responsive: true,
        displayModeBar: true,
        scrollZoom: true,
        modeBarButtonsToAdd: [
            {
                name: "Toggle Mesh",
                icon: Plotly.Icons.grid,
                click: () => {
                    showMesh = !showMesh;
                    draw();
                }
            },
            {
                name: "Auto Z Scale",
                icon: Plotly.Icons.autoscale,
                click: () => {
                    Plotly.relayout(container, { "zaxis.autorange": true });
                }
            }
        ]
    };

    Plotly.react(container, traces, layout, config);

    if (!container._viewListenerBound) {
        container.on('plotly_buttonclicked', (event) => {
            // View handling: Geometry(0), Potential(1), E-field(2)
            // Mode selector in sidebar controls odd/even for differential lines
            if (event.menu.active === 0) {
                currentView = "geometry";
            } else if (event.menu.active === 1) {
                currentView = "potential";
            } else if (event.menu.active === 2) {
                currentView = "efield";
            }
            draw();
        });
        container._viewListenerBound = true;
    }

}

function resizeCanvas() {
    const container = document.getElementById('sim_canvas');
    if (container) {
        Plotly.Plots.resize(container);
    }
}

function getYAxisLabel(selector) {
    const labels = {
        're_z0': 'Re(Z0) (Ohm)',
        'im_z0': 'Im(Z0) (Ohm)',
        'loss': 'Loss (dB/m)',
        'R': 'R (Ohm/m)',
        'L': 'L (H/m)',
        'C': 'C (F/m)',
        'G': 'G (S/m)'
    };
    return labels[selector] || selector;
}

function drawResultsPlot() {
    if (!frequencySweepResults || frequencySweepResults.length === 0) return;

    const selector = document.getElementById('results-plot-selector').value;
    const isDifferential = solver && solver.is_differential;

    const freqs = frequencySweepResults.map(r => r.freq / 1e9);
    const traces = [];

    // Use lines+markers mode so single frequency points are visible
    const plotMode = freqs.length === 1 ? 'markers' : 'lines+markers';

    if (selector === 're_z0') {
        if (isDifferential) {
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[0].Zc.re),
                name: 'Odd mode',
                type: 'scatter',
                mode: plotMode
            });
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[1].Zc.re),
                name: 'Even mode',
                type: 'scatter',
                mode: plotMode
            });
        } else {
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[0].Zc.re),
                name: 'Re(Z0)',
                type: 'scatter',
                mode: plotMode
            });
        }
    } else if (selector === 'im_z0') {
        if (isDifferential) {
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[0].Zc.im),
                name: 'Odd mode',
                type: 'scatter',
                mode: plotMode
            });
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[1].Zc.im),
                name: 'Even mode',
                type: 'scatter',
                mode: plotMode
            });
        } else {
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[0].Zc.im),
                name: 'Im(Z0)',
                type: 'scatter',
                mode: plotMode
            });
        }
    } else if (selector === 'loss') {
        if (isDifferential) {
            // Odd mode losses (solid lines)
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[0].alpha_c),
                name: 'Conductor (odd)',
                type: 'scatter',
                mode: plotMode
            });
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[0].alpha_d),
                name: 'Dielectric (odd)',
                type: 'scatter',
                mode: plotMode
            });
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[0].alpha_total),
                name: 'Total (odd)',
                type: 'scatter',
                mode: plotMode,
                line: { width: 2 }
            });
            // Even mode losses (dashed lines)
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[1].alpha_c),
                name: 'Conductor (even)',
                type: 'scatter',
                mode: plotMode,
                line: { dash: 'dash' }
            });
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[1].alpha_d),
                name: 'Dielectric (even)',
                type: 'scatter',
                mode: plotMode,
                line: { dash: 'dash' }
            });
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[1].alpha_total),
                name: 'Total (even)',
                type: 'scatter',
                mode: plotMode,
                line: { width: 2, dash: 'dash' }
            });
        } else {
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[0].alpha_c),
                name: 'Conductor',
                type: 'scatter',
                mode: plotMode
            });
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[0].alpha_d),
                name: 'Dielectric',
                type: 'scatter',
                mode: plotMode
            });
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[0].alpha_total),
                name: 'Total',
                type: 'scatter',
                mode: plotMode,
                line: { width: 2 }
            });
        }
    } else {
        // RLGC parameters
        const paramKey = selector; // R, L, G, or C
        if (isDifferential) {
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[0].RLGC[paramKey]),
                name: `${paramKey} (odd)`,
                type: 'scatter',
                mode: plotMode
            });
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[1].RLGC[paramKey]),
                name: `${paramKey} (even)`,
                type: 'scatter',
                mode: plotMode
            });
        } else {
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[0].RLGC[paramKey]),
                name: paramKey,
                type: 'scatter',
                mode: plotMode
            });
        }
    }

    const useLogX = document.getElementById('results-log-x').checked;
    const layout = {
        xaxis: {
            title: 'Frequency (GHz)',
            type: useLogX ? 'log' : 'linear'
        },
        yaxis: {
            title: getYAxisLabel(selector)
        },
        margin: { l: 80, r: 40, t: 40, b: 60 },
        showlegend: true,
        legend: { x: 0.02, y: 0.98 }
    };

    Plotly.newPlot('results-plot', traces, layout, { responsive: true });
}

function drawSParamPlot() {
    if (!frequencySweepResults || frequencySweepResults.length === 0) return;

    const length = getInputValue('sparam-length');
    const Z_ref = parseFloat(document.getElementById('sparam-z-ref').value);
    const isDifferential = solver && solver.is_differential;
    const plotMode = document.getElementById('sparam-plot-mode').value; // 'magnitude' or 'phase'

    const freqs = frequencySweepResults.map(r => r.freq / 1e9);
    const traces = [];

    // Use lines+markers mode so single frequency points are visible
    const lineMode = freqs.length === 1 ? 'markers' : 'lines+markers';

    // Helper to convert complex S-parameter to phase in degrees
    const sParamToPhase = (complexVal) => {
        return complexVal.arg() * 180 / Math.PI;
    };

    if (!isDifferential) {
        // 2-port S-parameters
        const S11_data = [];
        const S21_data = [];

        for (const { freq, result } of frequencySweepResults) {
            const sp = computeSParamsSingleEnded(freq, result.modes[0].RLGC, length, Z_ref);
            if (plotMode === 'magnitude') {
                S11_data.push(sParamTodB(sp.S11));
                S21_data.push(sParamTodB(sp.S21));
            } else {
                S11_data.push(sParamToPhase(sp.S11));
                S21_data.push(sParamToPhase(sp.S21));
            }
        }

        const label = plotMode === 'magnitude' ? '(dB)' : '(deg)';
        traces.push({
            x: freqs,
            y: S11_data,
            name: `S11 ${label}`,
            type: 'scatter',
            mode: lineMode
        });
        traces.push({
            x: freqs,
            y: S21_data,
            name: `S21 ${label}`,
            type: 'scatter',
            mode: lineMode
        });
    } else {
        // 4-port S-parameters (mixed-mode)
        const SDD11_data = [];
        const SDD21_data = [];
        const SCC11_data = [];
        const SCC21_data = [];

        for (const { freq, result } of frequencySweepResults) {
            const oddMode = result.modes.find(m => m.mode === 'odd');
            const evenMode = result.modes.find(m => m.mode === 'even');

            const sp = computeSParamsDifferential(
                freq,
                oddMode.RLGC,
                evenMode.RLGC,
                length,
                Z_ref
            );

            if (plotMode === 'magnitude') {
                SDD11_data.push(sParamTodB(sp.SDD11));
                SDD21_data.push(sParamTodB(sp.SDD21));
                SCC11_data.push(sParamTodB(sp.SCC11));
                SCC21_data.push(sParamTodB(sp.SCC21));
            } else {
                SDD11_data.push(sParamToPhase(sp.SDD11));
                SDD21_data.push(sParamToPhase(sp.SDD21));
                SCC11_data.push(sParamToPhase(sp.SCC11));
                SCC21_data.push(sParamToPhase(sp.SCC21));
            }
        }

        const label = plotMode === 'magnitude' ? '(dB)' : '(deg)';
        traces.push({
            x: freqs,
            y: SDD11_data,
            name: `SDD11 ${label}`,
            type: 'scatter',
            mode: lineMode
        });
        traces.push({
            x: freqs,
            y: SDD21_data,
            name: `SDD21 ${label}`,
            type: 'scatter',
            mode: lineMode
        });
        traces.push({
            x: freqs,
            y: SCC11_data,
            name: `SCC11 ${label}`,
            type: 'scatter',
            mode: lineMode,
            line: { dash: 'dash' }
        });
        traces.push({
            x: freqs,
            y: SCC21_data,
            name: `SCC21 ${label}`,
            type: 'scatter',
            mode: lineMode,
            line: { dash: 'dash' }
        });
    }

    const useLogX = document.getElementById('sparam-log-x').checked;
    const yTitle = plotMode === 'magnitude' ? 'Magnitude (dB)' : 'Phase (degrees)';
    const layout = {
        xaxis: {
            title: 'Frequency (GHz)',
            type: useLogX ? 'log' : 'linear'
        },
        yaxis: {
            title: yTitle
        },
        margin: { l: 80, r: 40, t: 40, b: 60 },
        showlegend: true,
        legend: { x: 0.02, y: 0.02 }
    };

    Plotly.newPlot('sparam-plot', traces, layout, { responsive: true });
}

function viridis(t) {
    // Simple heatmap approximation
    t = Math.max(0, Math.min(1, t));
    // R, G, B interpolation
    const r = Math.floor(255 * Math.sin(t * 2));
    const g = Math.floor(255 * Math.sin(t * 3));
    const b = Math.floor(255 * Math.cos(t * 1.5));
    // Better pseudocolor: (Blue -> Cyan -> Green -> Yellow -> Red)
    // Manual standard mapping for clarity:
    if(t < 0.25) return `rgb(0, ${Math.floor(t*4*255)}, 255)`;
    if(t < 0.5) return `rgb(0, 255, ${Math.floor((0.5-t)*4*255)})`;
    if(t < 0.75) return `rgb(${Math.floor((t-0.5)*4*255)}, 255, 0)`;
    return `rgb(255, ${Math.floor((1-t)*4*255)}, 0)`;
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
        sparamLength.addEventListener('change', () => {
            if (frequencySweepResults) {
                drawSParamPlot();
            }
        });
    }
    if (sparamZref) {
        sparamZref.addEventListener('change', () => {
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
        'freq-start': { min: 0.001, default: 0.1, label: 'Start frequency' },
        'freq-stop': { min: 0.001, default: 10, label: 'Stop frequency' },
        'inp_max_iters': { min: 1, default: 10, integer: true, label: 'Max iterations' },
        'inp_max_nodes': { min: 100, default: 20000, integer: true, label: 'Max nodes' },
        'inp_tolerance': { min: 0.0001, default: 0.05, label: 'Tolerance' },
        'sparam-length': { min: 0.0001, default: 0.01, label: 'Line length' },
        'sparam-z-ref': { min: 1, default: 50, label: 'Reference impedance' }
    };

    Object.entries(validationRules).forEach(([id, rule]) => {
        const el = document.getElementById(id);
        if (el) {
            el.addEventListener('blur', () => {
                let val = rule.integer ? parseInt(el.value) : parseFloat(el.value);
                if (isNaN(val) || val < rule.min || el.value.trim() === '') {
                    el.value = rule.default;
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
}

// Start when DOM is ready
window.addEventListener('DOMContentLoaded', init);
