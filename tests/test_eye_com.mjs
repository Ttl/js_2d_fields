// Statistical eye (src/eye/stat_eye.js through patternSim) against independent references.
//
// tests/com_reference.py builds its own pulse response of the same line (Howard Johnson's
// model by PyBERT's calc_gamma, S21 between matched 100 ohm ends, the COM forms of the
// Gaussian Tx edge, Tx FFE and CTLE from PyChOpMarg, the pulse by 93A-24) at the Eye
// engine's sample point and DFE taps rule, and evaluates it with:
//   - the COM statistical method (IEEE 802.3 Annex 93A, PyChOpMarg's delta_pmf): eye
//     height at the sample point and eye width at the threshold, DFE taps clipped to their
//     limits, at BER 1e-12, without jitter;
//   - a Monte Carlo of the PCIe jitter model, every Tx FIR level change moving with its
//     boundary's Rj and dual-Dirac jitter through the channel step response, at BER
//     levels it resolves (1e-3, 1e-4).
//
//   1. 16 GT/s, 0.6 m (~16 dB at 8 GHz), P7, CTLE -9 dB, 2 DFE taps.
//   2. 0.75 m, P7, CTLE -12 dB, 3 DFE taps, the first two clipped by their limits.
//   3. 0.4 m, P7, CTLE -9 dB, no DFE (eye width over the whole UI).
//   4. Case 1 with Table 8-14 jitter at 16 GT/s (Rj 0.71 ps, DjDD 3.75 ps): the Monte
//      Carlo. COM's own jitter term (per symbol, on the slope of the equalized pulse) is a
//      different model, printed for comparison only: it does not move the FIR edges, and
//      here takes about a third of the PCIe model's eye height loss.
// Without jitter the eye height agrees to ~0.3 % (the default 2049-bin grid; ~0.15 % at
// 16385) and the width to ~0.002 UI. COM's Rx FFE, packages, noise and crosstalk are not
// used: this checks the Eye engine's statistical eye, not COM or the PCIe values.
//
// Skipped when PyBERT or PyChOpMarg is not installed.
import { spawnSync } from 'child_process';
import { check, done } from './helpers.js';
import { loadFFT } from '../src/eye/fft.js';
import { simulateLink } from '../src/eye/transient.js';
import { patternBits } from '../src/eye/pattern.js';
import { PCIE_PRESETS, firFromPresetDb } from '../src/eye/eq.js';

if (spawnSync('python3', ['-c', 'import pybert, pychopmarg'], { stdio: 'ignore' }).status !== 0) {
    console.log('PyBERT / PyChOpMarg not found: skipped');
    process.exit(0);
}

const python = (input) => {
    const r = spawnSync('python3', ['tests/com_reference.py'], { input: JSON.stringify(input), maxBuffer: 1 << 28 });
    if (r.status !== 0) throw new Error('com_reference.py failed: ' + r.stderr.toString());
    return JSON.parse(r.stdout.toString());
};

const fft = await loadFFT();
const BIT_RATE = 16e9, TR = 20e-12;
const LINE = { R0: 1.452, w0: 10e6, Rdc: 0.1876, Z0: 100, v0: 0.67, Theta0: 0.02 };
const P7 = firFromPresetDb(...PCIE_PRESETS.P7);

function tableOf(line) {
    const tb = python({ table_only: true, line, f_max: 3e11 });
    return { differential: false, freqs: Float64Array.from(tb.freqs), keys: ['R', 'L', 'G', 'C'],
        data: { R: Float64Array.from(tb.R), L: Float64Array.from(tb.L), G: Float64Array.from(tb.G), C: Float64Array.from(tb.C) },
        lext: { L: tb.L_ext } };
}

// The single line stands for the differential pair: 100 ohm source and load, no pads.
function ours(table, c, ber) {
    const cfg = { length: c.line.l_ch, samplesPerEdge: 16,
        network: { rs: 100, padTx: 0, padRx: 0, tx: { ac: false }, rx: { ac: false }, term: { type: 'r', rt: 100 } },
        stimulus: { type: 'pattern', bitRate: BIT_RATE, bits: patternBits({ pattern: 'prbs7' }), fir: c.fir, vLow: 0, vHigh: 1,
            tr: TR, tf: TR, edgeDef: '20-80', shape: 'gauss', jitter: c.jitter || null },
        txEq: { name: null }, rxEq: { ctle: c.ctle, dfeTaps: c.dfe, dfeLimits: c.limits },
        statEye: { ber, offsets: [0] } };
    return simulateLink(table, cfg, fft);
}

function reference(c, r, mc) {
    return python({ bit_rate: BIT_RATE, ber: 1e-12, dfe_taps: c.dfe, dfe_limits: c.limits,
        x: Array.from(r.pulse.x), q: Array.from(r.pulse.y), M: 64, fstep: 10e6, line: c.line, z_ref: 100,
        T_r: TR, A_v: 0.25, fir: c.fir, ctle: c.ctle, jitter: c.jitter || null, mc: mc || null });
}

const mV = (v) => `${(v * 1e3).toFixed(2)} mV`;
const CTLE = (adcDb) => ({ adcDb, fp1: 2e9, fp2: 8e9 });
const CASES = [
    { name: '0.6 m, P7, CTLE -9 dB, 2 DFE taps', line: { ...LINE, l_ch: 0.6 }, fir: P7, ctle: CTLE(-9), dfe: 2, limits: [0.012, 0] },
    { name: '0.75 m, P7, CTLE -12 dB, 3 DFE taps (2 clipped)', line: { ...LINE, l_ch: 0.75 }, fir: P7, ctle: CTLE(-12), dfe: 3, limits: [0.004, 0.002, 0] },
    { name: '0.4 m, P7, CTLE -9 dB, no DFE', line: { ...LINE, l_ch: 0.4 }, fir: P7, ctle: CTLE(-9), dfe: 0, limits: [] },
];

// 1-3. Without jitter, against the COM method on the independent pulse.
for (const c of CASES) {
    const r = ours(tableOf(c.line), c, 1e-12);
    const ref = reference(c, r).ref;
    const se = r.statEye;
    check(`${c.name}: pulse response`, ref.pulse_diff < 1e-3 * ref.peak,
        `max diff ${(100 * ref.pulse_diff / ref.peak).toFixed(3)} % of the ${mV(ref.peak)} peak`);
    const dDfe = Math.max(0, ...r.eqChoice.dfe.map((h, k) => Math.abs(h - ref.dfe[k])));
    check(`${c.name}: DFE taps`, r.eqChoice.dfe.length === c.dfe && dDfe < 2e-3 * ref.peak,
        `[${r.eqChoice.dfe.map(mV).join(', ')}] vs [${ref.dfe.map(mV).join(', ')}]`);
    check(`${c.name}: eye height = COM method`, se.height > 0 && Math.abs(se.height - ref.height) < 0.01 * ref.height,
        `${mV(se.height)} vs ${mV(ref.height)}, ${(100 * (se.height / ref.height - 1)).toFixed(2)} %`);
    check(`${c.name}: eye width = COM method`, Math.abs(se.width - ref.width) < 0.005,
        `${se.width.toFixed(4)} vs ${ref.width.toFixed(4)} UI`);
}

// 4. Jitter: the Monte Carlo of the PCIe model. The eye top at BER b is the lower
// 2 b quantile of the one level (BER = P(V1 < v) / 2 + P(V0 > v) / 2).
{
    const c = { ...CASES[0], name: '0.6 m Gen4 jitter', jitter: { rj: 0.71e-12, djdd: 3.75e-12 } };
    const table = tableOf(c.line);
    const r = ours(table, c, 1e-12);
    const res = reference(c, r, { n: 1000000, probs: [1e-3, 1e-4], ew_ber: 1e-3, n_ew: 100000, post: 60 });
    const tol = [0.15e-3, 0.3e-3];
    res.mc.probs.forEach((p, i) => {
        const top = 0.5 * ours(table, c, p / 2).statEye.height;
        check(`${c.name}: eye top at BER ${p / 2} = Monte Carlo`, Math.abs(top - res.mc.top[i]) < tol[i],
            `${mV(top)} vs ${mV(res.mc.top[i])} (${res.mc.n} one bits)`);
    });
    const w = ours(table, c, 1e-3).statEye.width;
    check(`${c.name}: eye width at BER 1e-3 = Monte Carlo`, Math.abs(w - res.mc.width) < 0.01,
        `${w.toFixed(4)} vs ${res.mc.width.toFixed(4)} UI`);
    console.log(`  BER 1e-12: EH ${mV(r.statEye.height)}, EW ${r.statEye.width.toFixed(3)} UI; COM jitter model on the same pulse `
        + `${mV(res.com_jitter.height)}, ${res.com_jitter.width.toFixed(3)} UI; without jitter ${mV(res.ours.height)}, ${res.ours.width.toFixed(3)} UI`);
}

done();
