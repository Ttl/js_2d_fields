// Eye tab: stimulus, coupling and termination settings, the channel sweep request and
// the plots of the time-domain link simulation (eye_worker.js).
import { darkAxis, darkBackground } from '../plot.js';
import { channelTableFromRows, checkCoverage, fmtHz } from './channel_grid.js';
import { requiredBandwidth } from './stimulus.js';
import { patternBits } from './pattern.js';
import { lowFrequencyNotes, BOARD_WIDTH, GAP_ERR_MAX, SKIN_ERR } from './sensitivity.js';

// Accuracy notes the tab quantifies on the waveform instead of logging them.
const QUANTIFIED = new Set(['narrow-gap', 'skin-transition', 'ground-spreading']);
// Waveform change (fraction of the swing) above which such a note is a warning.
const LOW_FREQ_WARN = 0.01;
// Accuracy notes of sweep rows (discrete sweeps carry them per row).
const rowWarnings = (rows) => rows.flatMap(r => (r.result && r.result.warnings) || []);

const getPlotly = () => window.Plotly;
const $ = (id) => document.getElementById(id);

const COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#17becf'];

// Settings of the tab: key, input id, kind. 'eye' inputs are stored in the settings and
// links as SI numbers (qty says how to read and show them), 'text' keeps a select value.
const EYE_FIELDS = [
    ['eye_length', 'eye-length', 'eye', 'len'],
    ['eye_stim', 'eye-stim', 'text'],
    ['eye_delay', 'eye-delay', 'eye', 'time'],
    ['eye_width', 'eye-width', 'eye', 'time'],
    ['eye_clock_freq', 'eye-clock-freq', 'eye', 'freq'],
    ['eye_duty', 'eye-duty', 'eye', 'num'],
    ['eye_clock_mode', 'eye-clock-mode', 'text'],
    ['eye_cycles', 'eye-cycles', 'eye', 'num'],
    ['eye_vlow', 'eye-vlow', 'eye', 'volt'],
    ['eye_vhigh', 'eye-vhigh', 'eye', 'volt'],
    ['eye_tr', 'eye-tr', 'eye', 'time'],
    ['eye_tf', 'eye-tf', 'eye', 'time'],
    ['eye_edge_def', 'eye-edge-def', 'text'],
    ['eye_edge_shape', 'eye-edge-shape', 'text'],
    ['eye_rs', 'eye-rs', 'eye', 'res'],
    ['eye_skew', 'eye-skew', 'eye', 'time'],
    ['eye_tx_ac', 'eye-tx-ac', 'chk'],
    ['eye_ctx', 'eye-ctx', 'eye', 'cap'],
    ['eye_rx_ac', 'eye-rx-ac', 'chk'],
    ['eye_crx', 'eye-crx', 'eye', 'cap'],
    ['eye_term', 'eye-term', 'text'],
    ['eye_rt', 'eye-rt', 'eye', 'res'],
    ['eye_vtt', 'eye-vtt', 'eye', 'volt'],
    ['eye_ct', 'eye-ct', 'eye', 'cap'],
    ['eye_rdiff', 'eye-rdiff', 'eye', 'res'],
    ['eye_bitrate', 'eye-bitrate', 'eye', 'rate'],
    ['eye_pattern', 'eye-pattern', 'text'],
    ['eye_bits', 'eye-bits', 'text'],
    ['eye_nbits', 'eye-nbits', 'eye', 'num'],
    ['eye_spe', 'eye-spe', 'eye', 'num'],
    ['eye_span', 'eye-span', 'eye', 'time'],
];
const QTY = Object.fromEntries(EYE_FIELDS.filter(f => f[2] === 'eye').map(f => [f[1], f[3]]));

// Defaults, the same as the HTML values (links only carry what differs). NaN: empty
// (automatic span).
const EYE_DEFAULTS = {
    eye_length: 0.1, eye_stim: 'pulse', eye_delay: 0, eye_width: 1e-9,
    eye_clock_freq: 1e9, eye_duty: 50, eye_clock_mode: 'steady', eye_cycles: 8,
    eye_vlow: 0, eye_vhigh: 1, eye_tr: 35e-12, eye_tf: 35e-12, eye_edge_def: '10-90',
    eye_edge_shape: 'rcos', eye_rs: 50, eye_skew: 0, eye_tx_ac: 0, eye_ctx: 1e-7,
    eye_rx_ac: 0, eye_crx: 1e-7, eye_term: 'r', eye_rt: 50, eye_vtt: 1.2, eye_ct: 1e-7,
    eye_rdiff: 100, eye_spe: 16, eye_span: NaN,
    eye_bitrate: 10e9, eye_pattern: 'prbs7', eye_bits: '0011010111', eye_nbits: 32767,
};

// Whether a setting applies to the selected options; links leave out the others.
function eyeKeyActive(settings, key) {
    const v = (k) => settings[k] ?? EYE_DEFAULTS[k];
    switch (key) {
        case 'eye_width': return v('eye_stim') === 'pulse';
        case 'eye_clock_freq': case 'eye_duty': case 'eye_clock_mode': return v('eye_stim') === 'clock';
        case 'eye_cycles': return v('eye_stim') === 'clock' && v('eye_clock_mode') === 'burst';
        case 'eye_bitrate': case 'eye_pattern': return v('eye_stim') === 'pattern';
        case 'eye_bits': return v('eye_stim') === 'pattern' && v('eye_pattern') === 'custom';
        case 'eye_nbits': return v('eye_stim') === 'pattern' && ['prbs23', 'prbs31'].includes(v('eye_pattern'));
        case 'eye_delay': case 'eye_span': return v('eye_stim') !== 'pattern';
        case 'eye_ctx': return !!v('eye_tx_ac');
        case 'eye_crx': return !!v('eye_rx_ac');
        case 'eye_rt': return ['vtt', 'r', 'rc'].includes(v('eye_term'));
        case 'eye_vtt': return v('eye_term') === 'vtt';
        case 'eye_ct': return v('eye_term') === 'rc';
        case 'eye_rdiff': return v('eye_term') === 'diff';
        default: return true;
    }
}

// Settings value of an 'eye' input: SI number to 12 significant digits (so a typed
// "35 ps" equals the default 35e-12), NaN when empty or invalid.
function readEyeSetting(id) {
    try {
        const v = parseQuantity(id, QTY[id], id, { optional: true });
        return v === null ? NaN : parseFloat(v.toPrecision(12));
    } catch {
        return NaN;
    }
}

// Display units of each quantity: [threshold, scale, unit], the first whose threshold
// the magnitude reaches.
const SHOW = {
    len: [[0, 1e-3, 'mm']],
    time: [[1e-3, 1e-3, 'ms'], [1e-6, 1e-6, 'µs'], [1e-9, 1e-9, 'ns'], [0, 1e-12, 'ps']],
    freq: [[1e9, 1e9, 'GHz'], [1e6, 1e6, 'MHz'], [1e3, 1e3, 'kHz'], [0, 1, 'Hz']],
    cap: [[1e-6, 1e-6, 'µF'], [1e-9, 1e-9, 'nF'], [0, 1e-12, 'pF']],
    volt: [[0, 1, 'V']],
    rate: [[1e9, 1e9, 'Gb/s'], [1e6, 1e6, 'Mb/s'], [1e3, 1e3, 'kb/s'], [0, 1, 'b/s']],
    res: [[0, 1, '']],
    num: [[0, 1, '']],
};

// Put a settings value into an 'eye' input: a number with a fitting unit, or a string
// as it is (settings saved before the numeric form).
function writeEyeSetting(id, value) {
    const el = $(id);
    if (!el) return;
    if (typeof value === 'string') { el.value = value; return; }
    if (value === null || value === undefined || !isFinite(value)) { el.value = ''; return; }
    const a = Math.abs(value);
    const [, scale, unit] = SHOW[QTY[id]].find(([th]) => a >= th) || SHOW[QTY[id]].at(-1);
    const shown = parseFloat((value / scale).toPrecision(12));
    el.value = unit ? `${shown} ${unit}` : `${shown}`;
}

// ---------------------------------------------------------------- input parsing

const UNIT = {
    time: { s: 1, ms: 1e-3, us: 1e-6, 'μs': 1e-6, ns: 1e-9, ps: 1e-12, fs: 1e-15 },
    freq: { hz: 1, khz: 1e3, mhz: 1e6, ghz: 1e9, thz: 1e12 },
    cap: { f: 1, mf: 1e-3, uf: 1e-6, 'μf': 1e-6, nf: 1e-9, pf: 1e-12, ff: 1e-15 },
    res: { '': 1, ohm: 1, ohms: 1, 'ω': 1, k: 1e3, kohm: 1e3, 'kω': 1e3 },
    volt: { v: 1, mv: 1e-3 },
    rate: { 'b/s': 1, bps: 1, 'kb/s': 1e3, kbps: 1e3, 'mb/s': 1e6, mbps: 1e6, 'gb/s': 1e9, gbps: 1e9, 'tb/s': 1e12, tbps: 1e12 },
    pct: { '%': 1 },
    num: {},
};
const DEFAULT_UNIT = { time: 'ps', freq: 'ghz', cap: 'nf', res: 'ohm', volt: 'v', pct: '%', num: '', rate: 'gb/s' };

// Value of a text input in SI units. Throws with the field label on a bad value.
function parseQuantity(id, kind, label, { optional = false } = {}) {
    const el = $(id);
    const raw = (el ? el.value : '').trim().replace(/µ/g, 'μ');
    if (!raw) {
        if (optional) return null;
        throw new Error(`${label}: enter a value.`);
    }
    if (kind === 'len') {
        const v = window.parseValueWithUnit ? window.parseValueWithUnit(raw, 'mm') : parseFloat(raw) * 1e-3;
        if (!isFinite(v)) throw new Error(`${label}: "${raw}" is not a length.`);
        return v;
    }
    const m = raw.match(/^([+-]?(?:\d+\.?\d*|\.\d+)(?:e[+-]?\d+)?)\s*([a-zμω%/]*)$/i);
    if (!m) throw new Error(`${label}: "${raw}" is not a number.`);
    const unit = (m[2] || DEFAULT_UNIT[kind]).toLowerCase();
    const table = UNIT[kind];
    if (kind === 'num') {
        if (m[2]) throw new Error(`${label}: no unit expected.`);
        return parseFloat(m[1]);
    }
    if (!(unit in table)) throw new Error(`${label}: unknown unit "${m[2]}".`);
    return parseFloat(m[1]) * table[unit];
}

// ---------------------------------------------------------------- state

let app = null;
let eyeWorker = null, eyeJobId = 0;
const eyeJobs = new Map();
let busy = false, stopRequested = false;
let lastResult = null, frozen = null;
let channelCache = null;   // { key, fMax, rows, source }

function isBusy() { return busy; }

function getEyeWorker() {
    if (eyeWorker) return eyeWorker;
    // The literal "eye_worker.js" is rewritten by build.sh to add a content hash.
    const w = new Worker(new URL("eye_worker.js", import.meta.url), { type: 'module' });
    w.onmessage = (e) => {
        const m = e.data, job = eyeJobs.get(m.id);
        if (!job) return;
        eyeJobs.delete(m.id);
        if (m.type === 'done') job.resolve(m.result); else job.reject(new Error(m.message));
    };
    w.onerror = (e) => {
        for (const [id, job] of eyeJobs) { eyeJobs.delete(id); job.reject(new Error(e.message || 'Eye worker failed')); }
        w.terminate();
        if (eyeWorker === w) eyeWorker = null;
    };
    eyeWorker = w;
    return w;
}

function runEyeJob(table, cfg) {
    const w = getEyeWorker(), id = ++eyeJobId;
    return new Promise((resolve, reject) => {
        eyeJobs.set(id, { resolve, reject });
        w.postMessage({ id, type: 'run', table, cfg });
    });
}

// ---------------------------------------------------------------- settings

function readConfig() {
    const diff = app.isDifferential();
    const type = $('eye-stim').value;
    const stimulus = {
        type,
        vLow: parseQuantity('eye-vlow', 'volt', 'Low level'),
        vHigh: parseQuantity('eye-vhigh', 'volt', 'High level'),
        tr: parseQuantity('eye-tr', 'time', 'Rise time'),
        tf: parseQuantity('eye-tf', 'time', 'Fall time'),
        edgeDef: $('eye-edge-def').value,
        shape: $('eye-edge-shape').value,
        delay: parseQuantity('eye-delay', 'time', 'Start delay'),
        skew: diff ? parseQuantity('eye-skew', 'time', 'P/N skew') : 0,
    };
    if (type === 'pulse') stimulus.width = parseQuantity('eye-width', 'time', 'Pulse width');
    if (type === 'pattern') {
        stimulus.bitRate = parseQuantity('eye-bitrate', 'rate', 'Bit rate');
        const nbits = Math.round(parseQuantity('eye-nbits', 'num', 'Bits simulated'));
        if (!(nbits >= 2)) throw new Error('Bits simulated: at least 2.');
        stimulus.bits = patternBits({ pattern: $('eye-pattern').value, bits: $('eye-bits').value, nbits });
    }
    if (type === 'clock') {
        stimulus.clockFreq = parseQuantity('eye-clock-freq', 'freq', 'Clock frequency');
        stimulus.duty = parseQuantity('eye-duty', 'num', 'Duty cycle') / 100;
        stimulus.steady = $('eye-clock-mode').value === 'steady';
        stimulus.cycles = Math.round(parseQuantity('eye-cycles', 'num', 'Cycles'));
        if (!stimulus.steady && !(stimulus.cycles >= 1 && stimulus.cycles <= 2000)) {
            throw new Error('Cycles: between 1 and 2000.');
        }
    }
    const term = { type: $('eye-term').value };
    if (['vtt', 'r', 'rc'].includes(term.type)) term.rt = parseQuantity('eye-rt', 'res', 'Termination resistance');
    if (term.type === 'vtt') term.vtt = parseQuantity('eye-vtt', 'volt', 'V_TT');
    if (term.type === 'rc') term.ct = parseQuantity('eye-ct', 'cap', 'Termination capacitance');
    if (term.type === 'diff') {
        if (!diff) throw new Error('A differential termination needs a differential line.');
        term.rdiff = parseQuantity('eye-rdiff', 'res', 'Differential termination');
    }
    const network = {
        rs: parseQuantity('eye-rs', 'res', 'Driver resistance'),
        tx: { ac: $('eye-tx-ac').checked },
        rx: { ac: $('eye-rx-ac').checked },
        term,
    };
    if (network.tx.ac) network.tx.c = parseQuantity('eye-ctx', 'cap', 'TX coupling capacitor');
    if (network.rx.ac) network.rx.c = parseQuantity('eye-crx', 'cap', 'RX coupling capacitor');
    const spe = parseQuantity('eye-spe', 'num', 'Samples per edge');
    if (!(spe >= 4 && spe <= 256)) throw new Error('Samples per edge: between 4 and 256.');
    const span = parseQuantity('eye-span', 'time', 'Time span', { optional: true });
    const length = parseQuantity('eye-length', 'len', 'Line length');
    if (!(length > 0)) throw new Error('Line length must be positive.');
    return {
        length, network, stimulus, samplesPerEdge: spe, span: span || 0,
        causal: !!($('chk_causal_materials') && $('chk_causal_materials').checked),
    };
}

// Show the inputs that apply to the selected options.
function updateVisibility() {
    const diff = app.isDifferential();
    const type = $('eye-stim').value, term = $('eye-term').value;
    const show = (id, on) => { const el = $(id); if (el) el.style.display = on ? '' : 'none'; };
    show('eye-row-width', type === 'pulse');
    for (const id of ['eye-row-clock-freq', 'eye-row-duty', 'eye-row-clock-mode']) show(id, type === 'clock');
    show('eye-row-cycles', type === 'clock' && $('eye-clock-mode').value === 'burst');
    const pat = $('eye-pattern').value;
    for (const id of ['eye-row-bitrate', 'eye-row-pattern']) show(id, type === 'pattern');
    show('eye-row-bits', type === 'pattern' && pat === 'custom');
    show('eye-row-nbits', type === 'pattern' && (pat === 'prbs23' || pat === 'prbs31'));
    for (const id of ['eye-row-delay', 'eye-row-span']) show(id, type !== 'pattern');
    show('eye-row-skew', diff);
    show('eye-row-rt', ['vtt', 'r', 'rc'].includes(term));
    show('eye-row-vtt', term === 'vtt');
    show('eye-row-ct', term === 'rc');
    show('eye-row-rdiff', term === 'diff');
    $('eye-ctx').disabled = !$('eye-tx-ac').checked;
    $('eye-crx').disabled = !$('eye-rx-ac').checked;
    const diffOpt = $('eye-term').querySelector('option[value="diff"]');
    if (diffOpt) diffOpt.disabled = !diff;
    show('eye-show-legs-label', diff);
    show('eye-show-cm-label', diff);
    $('eye-term-label').textContent = diff ? 'per leg' : '';
    // Level readout and clock data rate.
    let note = '';
    try {
        const lo = parseQuantity('eye-vlow', 'volt', 'Low level'), hi = parseQuantity('eye-vhigh', 'volt', 'High level');
        note = diff ? `${fmtVolt(2 * Math.abs(hi - lo))} pp differential, ${fmtVolt((hi + lo) / 2)} common mode` : '';
    } catch { note = ''; }
    $('eye-level-note').textContent = note;
    let rate = '';
    try { if (type === 'clock') rate = `= ${+(2 * parseQuantity('eye-clock-freq', 'freq', 'f') / 1e9).toPrecision(4)} Gb/s data rate`; } catch { rate = ''; }
    $('eye-rate-note').textContent = rate;
    const hasEye = type === 'clock' || type === 'pattern';
    $('eye-view').querySelector('option[value="eye"]').disabled = !hasEye;
    if (!hasEye && $('eye-view').value === 'eye') $('eye-view').value = 'wave';
    const wg = app.isSelfReferenced();
    $('btn-eye-run').disabled = wg || busy;
    setNotice(wg ? ['Time-domain simulation needs a TEM or quasi-TEM line: not available for a waveguide.'] : null, 'eye-notice');
}

function fmtVolt(v) { return Math.abs(v) >= 1 ? `${+v.toPrecision(4)} V` : `${+(v * 1e3).toPrecision(4)} mV`; }

function setNotice(list, id = 'eye-warnings') {
    const el = $(id);
    if (!el) return;
    if (!list || !list.length) { el.style.display = 'none'; return; }
    el.innerHTML = list.map(w => `<div><strong>⚠</strong> ${escapeHtml(w)}</div>`).join('');
    el.style.display = 'block';
}
const escapeHtml = (s) => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

function setStatus(text) { $('eye-status').textContent = text || ''; }

// ---------------------------------------------------------------- channel

// Channel key: every solver parameter except the plot/solve frequency.
function channelKey(params) {
    const { freq, ...rest } = params;
    void freq;
    return JSON.stringify(rest);
}

// Sweep rows from DC to at least fNeed: the cache, the last solve, or a channel sweep.
async function getChannel(fNeed) {
    const params = app.getParams();
    const key = channelKey(params);
    if (channelCache && channelCache.key === key && channelCache.fMax >= fNeed) {
        return { rows: channelCache.rows, source: channelCache.source, warnings: channelCache.warnings };
    }
    const sw = app.getSweepState();
    let reason = 'there is no solve of the current geometry';
    if (sw.rows && sw.current) {
        const cov = checkCoverage(sw.rows, fNeed, sw.interp);
        if (cov.ok) {
            if (sw.interp && sw.keptJob != null) {
                try {
                    const out = await app.solveJob('eyeChannel', { mode: 'reuse', job: sw.keptJob, fNeed }, {});
                    const source = `the last solve (interpolating sweep, DC to ${fmtHz(out.fMax)})`;
                    channelCache = { key, fMax: out.fMax, rows: out.rows, source, warnings: out.warnings || [] };
                    return channelCache;
                } catch (e) {
                    reason = e.message;
                }
            } else if (!sw.interp) {
                const fMax = Math.max(...sw.rows.map(r => r.freq));
                const source = `the last solve (DC to ${fmtHz(fMax)})`;
                if (cov.warning) app.log(`⚠ Eye: ${cov.warning}.`);
                channelCache = { key, fMax, rows: sw.rows, source, warnings: rowWarnings(sw.rows) };
                return channelCache;
            } else {
                reason = 'the last solve is no longer kept';
            }
        } else {
            reason = cov.reason;
        }
    } else if (sw.rows) {
        reason = 'the geometry or frequencies changed since the last solve';
    }
    app.log(`Eye: running a channel sweep from DC to ${fmtHz(fNeed)}: ${reason}.`);
    setStatus(`Channel sweep to ${fmtHz(fNeed)}...`);
    const out = await app.solveJob('eyeChannel', {
        mode: 'solve', params, fNeed, interpTolerance: app.interpTolerance(),
    }, {
        progress: (m) => setStatus(m.text || 'Channel sweep...'),
    });
    if (out.stopped) return null;
    const logged = (out.warnings || []).filter(w => !QUANTIFIED.has(w && w.reason));
    for (const w of new Set(logged.map(w => w.message || String(w)))) app.log(`⚠ Eye channel: ${w}`);
    const source = `the channel sweep (DC to ${fmtHz(out.fMax)})`;
    channelCache = { key, fMax: out.fMax, rows: out.rows, source, warnings: out.warnings || [] };
    return channelCache;
}

// ---------------------------------------------------------------- low-frequency notes

const NOTE_NAME = { 'narrow-gap': 'narrow-gap loss', 'skin-transition': 'DC-skin transition', 'ground-spreading': 'ground spreading' };

// Warning text of a low-frequency note that changes the waveform noticeably.
function lowFreqText(lf) {
    const pct = `${+(100 * lf.impact).toPrecision(2)} %`;
    if (lf.reason === 'narrow-gap') {
        return `Below ${fmtHz(lf.band)} the quasi-static model overstates the skin-effect loss of closely facing conductors ` +
            `(up to ${100 * GAP_ERR_MAX} %): its estimated error changes this receiver waveform by ${pct} of the swing. ` +
            `The full-wave solver models the current in the gap.`;
    }
    if (lf.reason === 'skin-transition') {
        return `Below ${fmtHz(lf.band)} the quasi-static conductor loss and internal inductance are approximate (DC-skin ` +
            `transition): a ${100 * SKIN_ERR} % error there changes this receiver waveform by ${pct} of the swing. ` +
            `The full-wave solver resolves the transition-region current.`;
    }
    return `Below ${fmtHz(lf.band)} the ground return spreads sideways, modelled as an ideal return: the resulting ` +
        `inductance increase (up to a ${BOARD_WIDTH * 1e3} mm wide board) changes this receiver waveform by ${pct} of the swing. ` +
        `Use custom geometry with finite grounds when this matters.`;
}

// One-line summary of the checked notes for the info line.
function lowFreqSummary(list) {
    if (!list || !list.length) return '';
    return 'low-frequency model notes change the receiver by ' +
        list.map(lf => `${lf.impact < 1e-5 ? '< 0.001' : +(100 * lf.impact).toPrecision(2)} % (${NOTE_NAME[lf.reason]} below ${fmtHz(lf.band)})`).join(', ');
}

// ---------------------------------------------------------------- run

async function runEye() {
    if (busy) return;
    if (app.isBusy()) { setStatus('Wait for the running solve to finish.'); return; }
    let cfg, fNeed;
    try {
        cfg = readConfig();
        fNeed = requiredBandwidth(cfg.stimulus);
    } catch (e) {
        setStatus('');
        setNotice([e.message]);
        return;
    }
    busy = true; stopRequested = false;
    $('btn-eye-run').style.display = 'none';
    $('btn-eye-stop').style.display = '';
    setNotice(null);
    setStatus('Preparing...');
    const t0 = Date.now();
    try {
        const ch = await getChannel(fNeed);
        if (!ch || stopRequested) { setStatus('Stopped.'); return; }
        const table = channelTableFromRows(ch.rows);
        if (table.differential !== app.isDifferential()) throw new Error('The channel does not match the line type, solve again.');
        setStatus('Simulating...');
        const res = await runEyeJob(table, { ...cfg, notes: lowFrequencyNotes(ch.warnings) });
        if (stopRequested) { setStatus('Stopped.'); return; }
        res.source = ch.source;
        res.cfg = cfg;
        for (const lf of res.lowFreq || []) if (lf.impact > LOW_FREQ_WARN) res.warnings.push(lowFreqText(lf));
        lastResult = res;
        $('eye-export-btn').disabled = false;
        setNotice(res.warnings);
        for (const w of res.warnings) app.log(`⚠ Eye: ${w}`);
        draw();
        setStatus(`Done in ${((Date.now() - t0) / 1000).toFixed(1)} s`);
    } catch (e) {
        setStatus('');
        setNotice([e.message]);
        app.log('Eye simulation failed: ' + e.message);
    } finally {
        busy = false;
        $('btn-eye-run').style.display = '';
        $('btn-eye-stop').style.display = 'none';
        updateVisibility();
    }
}

function stopEye() {
    if (!busy) return;
    stopRequested = true;
    app.stopSolveJob();
    // The simulation itself is short. A long one is cut by replacing the worker.
    if (eyeWorker && eyeJobs.size) {
        eyeWorker.terminate(); eyeWorker = null;
        for (const [id, job] of eyeJobs) { eyeJobs.delete(id); job.reject(new Error('Stopped.')); }
    }
}

// ---------------------------------------------------------------- plots

function timeUnit(tMax) {
    if (tMax < 2e-9) return { s: 1e12, name: 'ps' };
    if (tMax < 2e-6) return { s: 1e9, name: 'ns' };
    if (tMax < 2e-3) return { s: 1e6, name: 'μs' };
    return { s: 1e3, name: 'ms' };
}

// Traces to draw, [key, label] by the toggles.
function selectedTraces(res) {
    const on = (id) => $(id) && $(id).checked;
    const out = [];
    if (!res.differential) {
        if (on('eye-show-rx')) out.push(['rx', 'RX']);
        if (on('eye-show-tx')) out.push(['tx', 'TX pad']);
        if (on('eye-show-src')) out.push(['src', 'Source EMF']);
        return out;
    }
    if (on('eye-show-rx')) out.push(['rxDiff', 'RX diff']);
    if (on('eye-show-tx')) out.push(['txDiff', 'TX pad diff']);
    if (on('eye-show-src')) out.push(['srcDiff', 'Source EMF diff']);
    if (on('eye-show-legs')) {
        if (on('eye-show-rx') || !on('eye-show-tx')) out.push(['rxP', 'RX P'], ['rxN', 'RX N']);
        if (on('eye-show-tx')) out.push(['txP', 'TX pad P'], ['txN', 'TX pad N']);
    }
    if (on('eye-show-cm')) out.push(['rxCM', 'RX common mode']);
    return out;
}

const LABEL_ORDER = ['rx', 'tx', 'src', 'rxDiff', 'txDiff', 'srcDiff', 'rxP', 'rxN', 'txP', 'txN', 'rxCM'];
const colorOf = (key) => COLORS[LABEL_ORDER.indexOf(key) % COLORS.length];

function layoutBase(xTitle, yTitle) {
    return {
        xaxis: darkAxis(xTitle), yaxis: darkAxis(yTitle),
        margin: { l: 70, r: 30, t: 30, b: 55 },
        showlegend: true, legend: { x: 1, xanchor: 'right', y: 1, font: { color: '#fff' }, bgcolor: 'rgba(0,0,0,0.3)' },
        ...darkBackground(),
    };
}

function drawWave(res) {
    const Plotly = getPlotly();
    const sel = selectedTraces(res);
    const tMax = Math.max(Math.abs(res.plot.t[0]), Math.abs(res.plot.t[res.plot.t.length - 1]));
    const u = timeUnit(tMax);
    const traces = [];
    const add = (r, faint) => {
        const x = Array.from(r.plot.t, v => v * u.s);
        for (const [key, label] of sel) {
            const y = r.plot.traces[key];
            if (!y) continue;
            traces.push({
                x, y: Array.from(y), name: faint ? `${label} (frozen)` : label, type: 'scatter', mode: 'lines',
                line: { color: colorOf(key), width: faint ? 1.5 : 2, dash: faint ? 'dot' : 'solid' },
                opacity: faint ? 0.5 : 1, legendgroup: key + (faint ? 'f' : ''),
            });
        }
    };
    if (frozen && frozen.kind === res.kind) add(frozen, true);
    add(res, false);
    const layout = layoutBase(`Time (${u.name})`, 'Voltage (V)');
    layout.shapes = (res.markers || []).map(m => ({
        type: 'line', x0: m.t * u.s, x1: m.t * u.s, yref: 'paper', y0: 0, y1: 1,
        line: { color: '#777', width: 1, dash: 'dot' },
    }));
    layout.annotations = (res.markers || []).map((m, i) => ({
        x: m.t * u.s, y: 1 - 0.05 * i, yref: 'paper', text: m.label, showarrow: false, xanchor: 'left', yanchor: 'top',
        font: { color: '#999', size: 10 },
    }));
    Plotly.newPlot('eye-plot', traces, layout, { responsive: true, modeBarButtonsToRemove: ['select2d', 'lasso2d'] });
}

// Density eye of a pattern: log counts over 2 UI, the eye height and width marked.
function drawEyeHeatmap(res) {
    const Plotly = getPlotly();
    const h = res.eyeHist, inf = res.eyeInfo;
    const z = [];
    for (let iy = 0; iy < h.vBins; iy++) {
        const row = new Array(h.tBins);
        for (let ix = 0; ix < h.tBins; ix++) {
            const c = h.counts[iy * h.tBins + ix];
            row[ix] = c > 0 ? Math.log10(1 + c) : null;
        }
        z.push(row);
    }
    const x = Array.from({ length: h.tBins }, (_, i) => 2 * (i + 0.5) / h.tBins);
    const y = Array.from({ length: h.vBins }, (_, i) => h.y0 + (h.y1 - h.y0) * (i + 0.5) / h.vBins);
    const traces = [{
        type: 'heatmap', x, y, z, showscale: false, hoverongaps: false, zsmooth: false,
        colorscale: [[0, '#1d2f6f'], [0.35, '#2c7fb8'], [0.7, '#f2c14e'], [1, '#ffffff']],
        hovertemplate: '%{x:.3f} UI<br>%{y:.4g} V<extra></extra>', name: 'RX',
    }];
    if (frozen && frozen.eyeHist && frozen.eyeHist.tBins === h.tBins) {
        const f = frozen.eyeHist, fz = [];
        for (let iy = 0; iy < f.vBins; iy++) {
            const row = new Array(f.tBins);
            for (let ix = 0; ix < f.tBins; ix++) row[ix] = f.counts[iy * f.tBins + ix] > 0 ? 1 : 0;
            fz.push(row);
        }
        traces.push({ type: 'contour', x, y: Array.from({ length: f.vBins }, (_, i) => f.y0 + (f.y1 - f.y0) * (i + 0.5) / f.vBins),
            z: fz, showscale: false, contours: { coloring: 'lines', start: 0.5, end: 0.5, size: 1 },
            line: { color: '#ddd', width: 1, dash: 'dot' }, name: 'Frozen', hoverinfo: 'skip' });
    }
    const layout = layoutBase(`Time (UI = ${fmtT(inf.ui)})`, 'Voltage (V)');
    layout.xaxis.range = [0, 2];
    layout.showlegend = false;
    layout.shapes = [];
    layout.annotations = [];
    if (isFinite(inf.heightTop)) {
        const open = inf.heightTop > inf.heightBottom;
        const mark = { color: open ? '#4ade80' : '#f87171', width: 2 };
        layout.shapes.push({ type: 'line', x0: 1, x1: 1, y0: inf.heightBottom, y1: inf.heightTop, line: mark });
        if (inf.widthRight > inf.widthLeft) {
            layout.shapes.push({ type: 'line', x0: inf.widthLeft, x1: inf.widthRight, y0: inf.thr, y1: inf.thr, line: mark });
        }
        layout.annotations.push({ x: 1, y: inf.heightTop, text: open ? 'eye height' : 'closed', showarrow: false,
            yanchor: 'bottom', font: { color: mark.color, size: 10 } });
    }
    Plotly.newPlot('eye-plot', traces, layout, { responsive: true, modeBarButtonsToRemove: ['select2d', 'lasso2d'] });
}

function drawEyeDiagram(res) {
    if (res.eyeHist) { drawEyeHeatmap(res); return; }
    const Plotly = getPlotly();
    const traces = [];
    const add = (r, faint) => {
        if (!r.eye) return;
        traces.push({
            x: Array.from(r.eye.x), y: Array.from(r.eye.y), type: 'scatter', mode: 'lines',
            name: faint ? 'Frozen' : (r.differential ? 'RX diff' : 'RX'),
            line: { color: faint ? '#aaa' : COLORS[0], width: 1.5, dash: faint ? 'dot' : 'solid' },
            opacity: faint ? 0.5 : 0.85, connectgaps: false,
        });
    };
    if (frozen) add(frozen, true);
    add(res, false);
    const ui = res.eye ? res.eye.ui : NaN;
    const layout = layoutBase(`Time (UI = ${isFinite(ui) ? fmtT(ui) : '-'})`, 'Voltage (V)');
    layout.xaxis.range = [0, 2];
    Plotly.newPlot('eye-plot', traces, layout, { responsive: true, modeBarButtonsToRemove: ['select2d', 'lasso2d'] });
}

function drawSpectrum(res) {
    const Plotly = getPlotly();
    const traces = [];
    const db = (a) => Array.from(a, v => 20 * Math.log10(Math.max(v, 1e-15)));
    const add = (r, faint) => {
        const x = Array.from(r.hplot.f, v => v / 1e9);
        traces.push({ x, y: db(r.hplot.mag), name: faint ? 'Transfer (frozen)' : (r.differential ? 'Transfer, diff RX / diff EMF' : 'Transfer, RX / EMF'),
            type: 'scatter', mode: 'lines', line: { color: COLORS[0], width: faint ? 1.5 : 2, dash: faint ? 'dot' : 'solid' }, opacity: faint ? 0.5 : 1 });
        if (!faint) traces.push({ x, y: db(r.hplot.stim), name: 'Rising edge spectrum', type: 'scatter', mode: 'lines',
            line: { color: COLORS[1], width: 1.5 } });
    };
    if (frozen) add(frozen, true);
    add(res, false);
    const layout = layoutBase('Frequency (GHz)', 'Magnitude (dB)');
    const vline = (f, label, color) => ({ f: f / 1e9, label, color });
    const lines = [vline(res.fNeed, 'needed bandwidth', '#888'), vline(res.fTop, 'sweep stop (extrapolated above)', '#d97706')];
    layout.shapes = lines.map(l => ({ type: 'line', x0: l.f, x1: l.f, yref: 'paper', y0: 0, y1: 1, line: { color: l.color, width: 1, dash: 'dash' } }));
    layout.annotations = lines.map((l, i) => ({ x: l.f, y: 1 - 0.06 * i, yref: 'paper', text: l.label, showarrow: false, xanchor: 'left', yanchor: 'top', font: { color: l.color, size: 10 } }));
    layout.yaxis.range = [-80, 10];
    // Up to a few times the needed bandwidth, beyond it the edge carries little.
    const fEnd = res.hplot.f.length ? res.hplot.f[res.hplot.f.length - 1] : 0;
    layout.xaxis.range = [0, Math.min(fEnd, 4 * Math.max(res.fNeed, res.fTop)) / 1e9];
    Plotly.newPlot('eye-plot', traces, layout, { responsive: true, modeBarButtonsToRemove: ['select2d', 'lasso2d'] });
}

function fmtT(s) {
    const a = Math.abs(s);
    if (a >= 1e-6) return `${+(s * 1e6).toPrecision(4)} μs`;
    if (a >= 1e-9) return `${+(s * 1e9).toPrecision(4)} ns`;
    return `${+(s * 1e12).toPrecision(4)} ps`;
}

function drawMetrics(res) {
    const el = $('eye-metrics');
    if (!res) { el.innerHTML = ''; return; }
    const rows = res.metrics.list.map(m =>
        `<div${m.note ? ` title="${escapeHtml(m.note)}"` : ''}><span>${escapeHtml(m.label)}</span><span>${escapeHtml(m.value)}</span></div>`).join('');
    const i = res.info;
    const tau = res.tau.map(fmtT).join(' / ');
    const info = `Line delay ${tau}${res.differential ? ' (odd / even)' : ''} · ` +
        (res.kind === 'steady' ? 'clock steady state'
            : res.kind === 'pattern' ? `pattern steady state, ${i.samplesPerBit} samples per bit, step responses over ${fmtT(i.window)}`
            : `window ${fmtT(i.window)}`) +
        ` · ${i.N} samples, step ${fmtT(i.dt)}` + (i.slowModes ? ` · ${i.slowModes} slow coupling modes` : '') +
        ` · line parameters from ${escapeHtml(res.source)}` +
        (res.lowFreq && res.lowFreq.length ? ` · ${escapeHtml(lowFreqSummary(res.lowFreq))}` : '');
    el.innerHTML = `<div class="eye-metric-grid">${rows}</div><div class="eye-info">${info}</div>`;
}

function draw() {
    const Plotly = getPlotly();
    if (!Plotly || !lastResult) return;
    const view = $('eye-view').value;
    if (view === 'eye' && (lastResult.eye || lastResult.eyeHist)) drawEyeDiagram(lastResult);
    else if (view === 'spectrum') drawSpectrum(lastResult);
    else drawWave(lastResult);
    drawMetrics(lastResult);
}

function exportCsv() {
    if (!lastResult) return;
    const { csv, differential } = lastResult;
    const keys = Object.keys(csv.traces);
    const lines = [`# Eye tab transient, ${differential ? 'differential' : 'single-ended'}, voltages in V`,
        ['time_s', ...keys].join(',')];
    for (let i = 0; i < csv.t.length; i++) {
        lines.push([csv.t[i].toExponential(9), ...keys.map(k => csv.traces[k][i].toPrecision(9))].join(','));
    }
    const blob = new Blob([lines.join('\n') + '\n'], { type: 'text/csv' });
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = 'transient.csv';
    a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}

function toggleFreeze() {
    const btn = $('eye-freeze-btn');
    if (frozen) { frozen = null; btn.textContent = 'Freeze'; btn.classList.remove('freeze-active'); }
    else if (lastResult) { frozen = lastResult; btn.textContent = 'Unfreeze'; btn.classList.add('freeze-active'); }
    draw();
}

// ---------------------------------------------------------------- setup

/**
 * ctx: { getParams, getSweepState, solveJob(type, payload, handlers), stopSolveJob,
 *        log, isBusy, isDifferential, isSelfReferenced, interpTolerance }
 */
function initEyeTab(ctx) {
    app = ctx;
    $('btn-eye-run').addEventListener('click', runEye);
    $('btn-eye-stop').addEventListener('click', stopEye);
    $('eye-export-btn').addEventListener('click', exportCsv);
    $('eye-freeze-btn').addEventListener('click', toggleFreeze);
    $('eye-view').addEventListener('change', draw);
    for (const id of ['eye-show-rx', 'eye-show-tx', 'eye-show-src', 'eye-show-legs', 'eye-show-cm']) {
        $(id).addEventListener('change', draw);
    }
    const panel = $('eye-settings');
    panel.addEventListener('change', updateVisibility);
    panel.addEventListener('input', updateVisibility);
    panel.addEventListener('keydown', (e) => { if (e.key === 'Enter' && e.target.tagName === 'INPUT') runEye(); });
}

function onEyeTabShown() {
    updateVisibility();
    const Plotly = getPlotly(), c = $('eye-plot');
    if (lastResult) draw();
    else if (Plotly && c && !c.data) {
        Plotly.newPlot('eye-plot', [], { ...layoutBase('Time', 'Voltage (V)'),
            annotations: [{ text: 'Set the stimulus and the link, then Simulate.', showarrow: false, xref: 'paper', yref: 'paper', x: 0.5, y: 0.5, font: { color: '#888' } }] },
            { responsive: true });
    }
}

export { initEyeTab, onEyeTabShown, isBusy as eyeBusy, updateVisibility as updateEyeVisibility, EYE_FIELDS, EYE_DEFAULTS,
    eyeKeyActive, readEyeSetting, writeEyeSetting, parseQuantity };
