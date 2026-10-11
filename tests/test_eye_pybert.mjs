// Time-domain link simulation (src/eye) against PyBERT.
//
// PyBERT's native channel (Howard Johnson's line model with its differential source and
// load resistances and parasitic capacitances) is the reference. tests/pybert_reference.py
// hands over the line's R, L, G, C from PyBERT's own gamma and Zc, so both simulate the
// same line, and PyBERT's channel impulse response driven by the same Gaussian-edged
// pulse, with the Tx FFE taps and a PCIe 3.0 reference-form CTLE applied on its side.
// The differential channel is a single line here: R_s and R_L as given, each parasitic
// capacitance split in two in series by PyBERT, so C / 2 per end.
//
//   1. Channel pulse response, transient path.
//   2. Pulse responses of the pattern engine (the one equalization uses): no FFE, PCIe
//      P7, P8, a CTLE of -6 dB, P7 with a CTLE of -9 dB.
// Both are compared after the best time shift (PyBERT trims its impulse response). They
// agree to 0.2-0.4 % of the peak, unchanged by finer sampling on either side or a denser
// line-parameter table: the remainder is PyBERT's resampling of its impulse response.
// PyBERT has no PCIe settings of its own and its CTLE-from-file path is marked broken,
// so the CTLE is applied on its side by tests/pybert_reference.py: this checks the Eye
// engine's channel, terminations, FFE and CTLE application, not the PCIe values.
//
// Skipped when PyBERT is not installed.
import { spawnSync } from 'child_process';
import { check, done } from './helpers.js';
import { loadFFT } from '../src/eye/fft.js';
import { simulateLink } from '../src/eye/transient.js';
import { patternBits } from '../src/eye/pattern.js';
import { PCIE_PRESETS, firFromPresetDb } from '../src/eye/eq.js';

if (spawnSync('python3', ['-c', 'import pybert'], { stdio: 'ignore' }).status !== 0) {
    console.log('PyBERT not found: skipped');
    process.exit(0);
}

const fft = await loadFFT();
const BASE = {
    bit_rate_gbps: 10, nspui: 64, f_max_ghz: 160, f_step_mhz: 5, l_ch: 0.5,
    R0: 1.452, w0: 10e6, Rdc: 0.1876, Z0: 100, v0: 0.67, Theta0: 0.02,
    rs: 100, cout_pf: 0.5, rin: 100, cin_pf: 0.5, f_tab_max: 2e11, edge_sigma: 5e-12, t0: 1e-10,
};

function reference(extra) {
    const r = spawnSync('python3', ['tests/pybert_reference.py'], { input: JSON.stringify({ ...BASE, ...extra }), maxBuffer: 1 << 28 });
    if (r.status !== 0) throw new Error('pybert_reference.py failed: ' + r.stderr.toString());
    return JSON.parse(r.stdout.toString());
}

const ref0 = reference({ taps: [0, 1, 0] });
const tb = ref0.table;
const table = {
    differential: false, freqs: Float64Array.from(tb.freqs), keys: ['R', 'L', 'G', 'C'],
    data: { R: Float64Array.from(tb.R), L: Float64Array.from(tb.L), G: Float64Array.from(tb.G), C: Float64Array.from(tb.C) },
    lext: { L: tb.L_ext },
};
const ui = ref0.ui, ts = ref0.ts;
// Gaussian edges of PyBERT's sigma (10-90 % = 2.563 sigma).
const tr = 2 * 1.2815515655446004 * BASE.edge_sigma;
const network = { rs: BASE.rs, padTx: BASE.cout_pf * 1e-12 / 2, padRx: BASE.cin_pf * 1e-12 / 2,
    tx: { ac: false }, rx: { ac: false }, term: { type: 'r', rt: BASE.rin } };

const peakOf = (y) => { let i = 0; for (let k = 1; k < y.length; k++) if (y[k] > y[i]) i = k; return i; };
// Largest difference of (xo, yo) from the reference (xr, yr) over [x0, x1] after the
// best time shift: absolute delays do not compare (PyBERT trims its impulse response),
// and aligning on the largest sample alone is off by up to half a sample.
function bestShiftDiff(xo, yo, xr, yr, x0, x1, step) {
    const diffAt = (sh) => {
        let e = 0;
        for (let i = 0; i < xo.length; i++) {
            const x = xo[i] + sh;
            if (x < x0 || x > x1) continue;
            const v = interp(xr, yr, x);
            if (isFinite(v)) e = Math.max(e, Math.abs(yo[i] - v));
        }
        return e;
    };
    let best = Infinity;
    for (let k = -50; k <= 50; k++) best = Math.min(best, diffAt(k * step / 50));
    return best;
}
const interp = (xs, ys, x) => {
    let lo = 0, hi = xs.length - 1;
    if (x <= xs[0] || x >= xs[hi]) return NaN;
    while (hi - lo > 1) { const m = (lo + hi) >> 1; if (xs[m] <= x) lo = m; else hi = m; }
    return ys[lo] + (ys[hi] - ys[lo]) * (x - xs[lo]) / (xs[hi] - xs[lo]);
};

// 1. Channel pulse, transient path: one UI between 50 % points. PyBERT trims the front of
// its channel impulse response, so both are aligned on their peaks.
{
    const stimulus = { type: 'pulse', width: ui, vLow: 0, vHigh: 1, tr, tf: tr, edgeDef: '10-90', shape: 'gauss', delay: 0 };
    const cfg = { length: BASE.l_ch, samplesPerEdge: 16, network, stimulus, span: 3e-9 };
    const r = simulateLink(table, cfg, fft);
    const p = ref0.pulse, pr = peakOf(p), peak = p[pr], rx = Array.from(r.traces.rx), po = peakOf(rx);
    const xo = Array.from(r.t, t => t - r.t[po]), xr = p.map((_, k) => (k - pr) * ts);
    const e = bestShiftDiff(xo, rx, xr, p, -3 * ui, 20 * ui, ts);
    check('channel pulse response = PyBERT (transient path)', e < 0.01 * peak, `max diff ${(100 * e / peak).toFixed(3)} % of the ${(peak * 1e3).toFixed(1)} mV peak`);
}

// 2. Pattern engine pulse responses (2 q, q the half response used by the equalizers).
const cases = [
    ['no FFE', [0, 1, 0], null],
    ['PCIe P7', firFromPresetDb(...PCIE_PRESETS.P7), null],
    ['PCIe P8', firFromPresetDb(...PCIE_PRESETS.P8), null],
    ['CTLE -6 dB, 2 / 8 GHz', [0, 1, 0], { adcDb: -6, fp1: 2e9, fp2: 8e9 }],
    ['PCIe P7 + CTLE -9 dB', firFromPresetDb(...PCIE_PRESETS.P7), { adcDb: -9, fp1: 2e9, fp2: 8e9 }],
];
for (const [name, taps, ctle] of cases) {
    const ref = reference({ taps, ctle });
    const yRef = ctle ? ref.pulse_ctle : ref.pulse;
    const cfg = { length: BASE.l_ch, samplesPerEdge: 16, network,
        stimulus: { type: 'pattern', bitRate: 1 / ui, bits: patternBits({ pattern: 'prbs7' }), fir: taps,
            vLow: 0, vHigh: 1, tr, tf: tr, edgeDef: '10-90', shape: 'gauss' },
        txEq: { name: null }, rxEq: { ctle, dfeTaps: 1 } };
    const r = simulateLink(table, cfg, fft);
    // Both on a time axis from their own peak, in UI.
    const xo = Array.from(r.pulse.x), yo = Array.from(r.pulse.y, v => 2 * v);
    const po = peakOf(yo), pr = peakOf(yRef);
    const xr = yRef.map((_, k) => (k - pr) * ts / ui);
    const peak = yRef[pr];
    const e = bestShiftDiff(xo.map(x => x - xo[po]), yo, xr, yRef, -3, 30, ts / ui);
    check(`pulse response = PyBERT: ${name}`, e < 0.01 * peak, `max diff ${(100 * e / peak).toFixed(3)} % of the ${(peak * 1e3).toFixed(1)} mV peak`);
}

done();
