// Time-domain link simulation (src/eye) against analytic circuits.
//
//   1. Lossless matched line: the receiver is the stimulus / 2 delayed by l sqrt(LC).
//   2. DC levels of every termination, V_TT included, against a resistive divider.
//   3. Mismatched lossless line: lattice-diagram reflection levels.
//   4. Skin-effect line: step response vs 0.5 erfc(beta / (2 sqrt(t - tau))).
//   5. RC termination on a short line: the lumped RC step response.
//   6. AC coupling: slow-mode split vs the brute-force long window, droop vs first order.
//   7. Clock steady state vs the last period of a long burst from rest.
//   8. Symmetric pair: differential response = twice the odd-mode single-ended one.
//      Line + network transfer vs the MTL / single-ended S-parameters.
//   9. Coverage rule of reused sweeps.
//  10. FFT round trip, non-power-of-two sizes, sampling invariance.
//  11. Non-causal dielectric model: precursor warning.
//  12. Open end, strongly mismatched source, skin-effect line: the low-frequency
//      ringing settles inside the window (no wrap-around before the arrival).
//  13. Low-frequency model notes (narrow gap, skin transition, ground spreading): the bound on the
//      receiver change holds against perturbed steady-state runs and separates GHz
//      links (negligible) from slow, low-impedance ones.
//  14. Bit patterns: PRBS generators maximal length; the periodic pattern engine equals
//      the clock steady state for clock-like patterns (DC and AC coupled, a pair with
//      skew); eye metrics and histogram on a lossless line; AC-coupled baseline wander.
//  15. Drivers, pads, terminations and equalization: LVDS current-mode levels, Thevenin
//      termination levels, RX pad capacitance as an RC; PCIe preset taps, FIR source
//      levels, a CTLE of 0 dB with far poles leaves the signal, the CTLE's DC gain; DFE
//      and the automatic choices beat no equalization.
import { check, done } from './helpers.js';
import { loadFFT, smoothSize } from '../src/eye/fft.js';
import { simulateLink } from '../src/eye/transient.js';
import { buildStimulus, edgeValue, erfc } from '../src/eye/stimulus.js';
import { makeRLGCEvaluator, checkCoverage } from '../src/eye/channel_grid.js';
import { buildNetwork, makeTransfer } from '../src/eye/network.js';
import { lowFrequencyImpact, lowFrequencyNotes, perturbation } from '../src/eye/sensitivity.js';
import { prbsBits, prbsPeriod, patternBits, PRBS } from '../src/eye/pattern.js';
import { PCIE_PRESETS, REDUCED_SWING_PRESETS, boostDb, firFromPresetDb, presetDbFromFir, ctle2Transfer } from '../src/eye/eq.js';
import { computeSParamsSingleEnded, computeSParamsDifferentialMTL } from '../src/sparameters.js';

const fft = await loadFFT();

const FREQS = [0, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12];
// Frequency-dependent parameters need a dense table: splining a decade-spaced one breaks
// the R/L Kramers-Kronig pairing and shows as a small precursor.
const DENSE = [0, ...Array.from({ length: 121 }, (_, i) => 1e6 * Math.pow(10, i / 20))];
function seTable(fn, freqs = FREQS) {
    const data = { R: [], L: [], G: [], C: [] };
    for (const f of freqs) { const v = fn(f); for (const k in data) data[k].push(v[k]); }
    return { differential: false, freqs: Float64Array.from(freqs), keys: ['R', 'L', 'G', 'C'],
        data: Object.fromEntries(Object.entries(data).map(([k, v]) => [k, Float64Array.from(v)])),
        lext: { L: fn(freqs[freqs.length - 1]).L } };
}
function pairTable(fn) {
    const keys = ['R11', 'R12', 'R22', 'L11', 'L12', 'L22', 'G11', 'G12', 'G22', 'C11', 'C12', 'C22'];
    const data = Object.fromEntries(keys.map(k => [k, []]));
    for (const f of FREQS) { const v = fn(f); for (const k of keys) data[k].push(v[k]); }
    const top = fn(FREQS[FREQS.length - 1]);
    return { differential: true, freqs: Float64Array.from(FREQS), keys,
        data: Object.fromEntries(Object.entries(data).map(([k, v]) => [k, Float64Array.from(v)])),
        lext: { L11: top.L11, L12: top.L12, L22: top.L22 } };
}
// Physical 2x2 entries of a symmetric pair from its odd and even mode values.
const fromModes = (o, e) => {
    const r = {};
    for (const k of ['R', 'L', 'G', 'C']) {
        r[k + '11'] = r[k + '22'] = (o[k] + e[k]) / 2;
        r[k + '12'] = (e[k] - o[k]) / 2;
    }
    return r;
};

const L0 = 250e-9, C0 = 100e-12;           // 50 Ohm, 2e8 m/s
const lossless = seTable(() => ({ R: 0, L: L0, G: 0, C: C0 }));
const base = (over = {}) => ({
    length: 0.1,
    network: { rs: 50, tx: { ac: false }, rx: { ac: false }, term: { type: 'r', rt: 50 }, ...(over.network || {}) },
    stimulus: { type: 'step', vLow: 0, vHigh: 1, tr: 30e-12, tf: 30e-12, edgeDef: '10-90', shape: 'gauss', delay: 0, ...(over.stimulus || {}) },
    samplesPerEdge: over.samplesPerEdge || 16, span: over.span, noSlowSplit: over.noSlowSplit, causal: over.causal,
});
const maxErr = (r, name, fn, from = -Infinity, to = Infinity) => {
    let e = 0;
    for (let i = 0; i < r.t.length; i++) {
        if (r.t[i] < from || r.t[i] > to) continue;
        e = Math.max(e, Math.abs(r.traces[name][i] - fn(r.t[i])));
    }
    return e;
};
const valueAt = (r, name, t) => {
    const ts = r.t; let i = 1;
    while (i < ts.length - 1 && ts[i] < t) i++;
    const a = (t - ts[i - 1]) / (ts[i] - ts[i - 1]);
    return r.traces[name][i - 1] + a * (r.traces[name][i] - r.traces[name][i - 1]);
};

// 1. Matched lossless line.
{
    const cfg = base();
    const r = simulateLink(lossless, cfg, fft);
    const st = buildStimulus({ ...cfg.stimulus, differential: false });
    const tau = 0.1 * Math.sqrt(L0 * C0);
    const eRx = maxErr(r, 'rx', t => 0.5 * edgeValue(st.up, t - st.t0 - tau));
    const eTx = maxErr(r, 'tx', t => 0.5 * edgeValue(st.up, t - st.t0));
    check('matched line: RX = stimulus/2 delayed by l sqrt(LC)', eRx < 1e-5, `max err ${eRx.toExponential(2)} V`);
    check('matched line: TX pad = stimulus/2', eTx < 1e-5, `max err ${eTx.toExponential(2)} V`);
    const d = r.metrics.list.find(m => m.label.startsWith('Delay'));
    check('matched line: delay metric', d && Math.abs(parseFloat(d.value) - 500) < 0.5, d && d.value);
    check('matched line: no warnings', r.warnings.length === 0, r.warnings.join(' | '));
    for (const shape of ['ramp', 'rcos']) {
        const c2 = base({ stimulus: { shape } });
        const r2 = simulateLink(lossless, c2, fft);
        const s2 = buildStimulus({ ...c2.stimulus, differential: false });
        const e = maxErr(r2, 'rx', t => 0.5 * edgeValue(s2.up, t - s2.t0 - tau));
        // Band-limit ringing at the ramp corners, the raised cosine has continuous slope.
        const tol = shape === 'ramp' ? 5e-3 : 1e-4;
        check(`matched line, ${shape} edge`, e < tol, `max err ${e.toExponential(2)} V`);
    }
}

// 2. DC levels: lossy line R_dc = 5 Ohm/m (0.5 Ohm over 0.1 m).
{
    const lossy = seTable(() => ({ R: 5, L: L0, G: 0, C: C0 }));
    const Rl = 0.5, Rs = 50;
    const cases = [
        { term: { type: 'r', rt: 50 }, lo: 0, hi: 50 / (Rs + Rl + 50) },
        { term: { type: 'vtt', rt: 50, vtt: 1.2 }, lo: 1.2 * (Rs + Rl) / (Rs + Rl + 50), hi: 1.2 * (Rs + Rl) / (Rs + Rl + 50) + 50 / (Rs + Rl + 50) },
        { term: { type: 'rc', rt: 50, ct: 1e-9 }, lo: 0, hi: 1 },
        { term: { type: 'open' }, lo: 0, hi: 1 },
    ];
    for (const c of cases) {
        // The RC termination charges over 50 ns, the open line settles slowly through the
        // small mismatch of its lossy Zc at low frequencies.
        const span = { rc: 2e-6, open: 50e-9 }[c.term.type];
        const r = simulateLink(lossy, base({ network: { term: c.term }, span }), fft);
        const v = r.traces.rx;
        const e0 = Math.abs(v[0] - c.lo), e1 = Math.abs(v[v.length - 1] - c.hi);
        check(`DC levels, ${c.term.type} termination`, e0 < 1e-9 && e1 < 2e-4, `initial ${v[0].toFixed(6)} (${c.lo.toFixed(6)}), final ${v[v.length - 1].toFixed(6)} (${c.hi.toFixed(6)})`);
    }
}

// 3. Mismatched lossless line: Rs = 25, Z0 = 50, RL = 100.
{
    const r = simulateLink(lossless, base({ network: { rs: 25, term: { type: 'r', rt: 100 } } }), fft);
    const st = buildStimulus({ ...base().stimulus, differential: false });
    const tau = 0.5e-9, t0 = st.t0;
    const gs = -1 / 3, gl = 1 / 3, v1 = 2 / 3;
    // Plateaus between arrivals: RX steps at tau, 3 tau, ...; TX at 0, 2 tau, ...
    const rx = [v1 * (1 + gl), v1 * (1 + gl) * (1 + gl * gs)];
    const tx = [v1, v1 * (1 + gl * (1 + gs))];
    const eRx2 = Math.abs(valueAt(r, 'rx', t0 + 4 * tau) - rx[1]);
    const eT = Math.max(Math.abs(valueAt(r, 'tx', t0 + tau) - tx[0]), Math.abs(valueAt(r, 'tx', t0 + 3 * tau) - tx[1]));
    check('mismatched line: RX first and second plateau', Math.abs(valueAt(r, 'rx', t0 + 2 * tau) - rx[0]) < 1e-5 && eRx2 < 1e-5,
        `${valueAt(r, 'rx', t0 + 2 * tau).toFixed(6)} / ${valueAt(r, 'rx', t0 + 4 * tau).toFixed(6)}`);
    check('mismatched line: TX plateaus', eT < 1e-5, `${valueAt(r, 'tx', t0 + tau).toFixed(6)} / ${valueAt(r, 'tx', t0 + 3 * tau).toFixed(6)}`);
    const fin = r.traces.rx[r.traces.rx.length - 1];
    check('mismatched line: settles to the divider', Math.abs(fin - 0.8) < 1e-4, fin.toFixed(6));
}

// 4. Skin-effect line, R = k sqrt(f) with equal internal reactance, matched to Z0.
{
    const k = 2e-3;   // Ohm/m/sqrt(Hz): 63 Ohm/m at 1 GHz
    const len = 1;
    const skin = seTable(f => {
        const ff = Math.max(f, 1e6);
        return { R: k * Math.sqrt(ff), L: L0 + k * Math.sqrt(ff) / (2 * Math.PI * ff), G: 0, C: C0 };
    });
    const cfg = { ...base({ stimulus: { tr: 5e-12 }, samplesPerEdge: 8 }), length: len };
    const r = simulateLink(skin, cfg, fft);
    const st = buildStimulus({ ...cfg.stimulus, differential: false });
    const tau = len * Math.sqrt(L0 * C0), Z0 = 50;
    const beta = len * k / (2 * Z0 * Math.sqrt(Math.PI));
    let e = 0;
    for (const dtx of [20e-12, 50e-12, 100e-12, 200e-12, 500e-12]) {
        const t = st.t0 + tau + dtx;
        e = Math.max(e, Math.abs(valueAt(r, 'rx', t) - 0.5 * erfc(beta / (2 * Math.sqrt(dtx)))));
    }
    // The erfc form is the high-frequency limit (Zc = Z0, gamma linear in Zs): about
    // 2 % of the swing from the Zc mismatch reflections.
    check('skin-effect line: step response vs erfc', e < 2e-2, `max err ${e.toExponential(2)} V of 0.5 V`);
}

// 5. RC termination, 1 mm line: v_d = 1 - Rs/(Rs+Rt) exp(-t/((Rs+Rt)C)).
{
    const cfg = { ...base({ network: { term: { type: 'rc', rt: 50, ct: 10e-12 } }, stimulus: { tr: 10e-12 } }), length: 1e-3 };
    const r = simulateLink(lossless, cfg, fft);
    const st = buildStimulus({ ...cfg.stimulus, differential: false });
    const tc = 100 * 10e-12;
    const e = maxErr(r, 'rx', t => t < st.t0 ? 0 : 1 - 0.5 * Math.exp(-(t - st.t0) / tc), st.t0 + 100e-12, st.t0 + 5e-9);
    check('RC termination: lumped step response', e < 5e-3, `max err ${e.toExponential(2)} V`);
}

// 6. AC coupling: 1 nF TX and RX caps, 50 Ohm, 1 ns pulse.
{
    const net = { tx: { ac: true, c: 1e-9 }, rx: { ac: true, c: 1e-9 } };
    const stim = { type: 'pulse', width: 1e-9 };
    const span = 40e-9;
    const split = simulateLink(lossless, base({ network: net, stimulus: stim, span }), fft);
    const brute = simulateLink(lossless, base({ network: net, stimulus: stim, span, noSlowSplit: true }), fft);
    // The brute-force run keeps only the 1 GOhm leak of the floating line slow.
    check('AC coupling: split uses slow modes', split.info.slowModes > brute.info.slowModes,
        `${split.info.slowModes} slow modes, windows ${(split.info.window * 1e9).toFixed(1)} / ${(brute.info.window * 1e9).toFixed(1)} ns`);
    let e = 0;
    for (const t of [0.6e-9, 0.8e-9, 1.2e-9, 1.5e-9, 2e-9, 5e-9, 20e-9, 39e-9]) {
        e = Math.max(e, Math.abs(valueAt(split, 'rx', t) - valueAt(brute, 'rx', t)));
    }
    check('AC coupling: split matches the brute-force window', e < 2e-4, `max diff ${e.toExponential(2)} V`);
    // Droop over the pulse: the series loop time constant is 100 Ohm * 0.5 nF.
    const tauAC = 100 * 0.5e-9;
    const st = buildStimulus({ ...base({ stimulus: stim }).stimulus, differential: false });
    const tRx = st.t0 + 0.5e-9;
    const vTop = valueAt(split, 'rx', tRx + 0.8e-9), vStart = valueAt(split, 'rx', tRx + 0.2e-9);
    const expect = 0.5 * (Math.exp(-0.8e-9 / tauAC) - Math.exp(-0.2e-9 / tauAC));
    check('AC coupling: droop over the pulse', Math.abs((vTop - vStart) - expect) < 1e-4,
        `${((vTop - vStart) * 1e3).toFixed(3)} mV vs ${(expect * 1e3).toFixed(3)} mV`);
}

// 7. Clock steady state vs burst, mismatched lossy line.
{
    const lossy = seTable(() => ({ R: 20, L: L0, G: 0, C: C0 }));
    const net = { rs: 30, term: { type: 'r', rt: 80 } };
    const clk = { type: 'clock', clockFreq: 1e9, duty: 0.5 };
    const steady = simulateLink(lossy, base({ network: net, stimulus: { ...clk, steady: true } }), fft);
    const burst = simulateLink(lossy, base({ network: net, stimulus: { ...clk, steady: false, cycles: 30 } }), fft);
    const st = buildStimulus({ ...base({ stimulus: { ...clk, cycles: 30 } }).stimulus, differential: false });
    // The steady state time axis is the source period from 0; the burst one has the
    // same edges shifted by its start (t0 - halfMax of the steady stimulus is identical).
    const sst = buildStimulus({ ...base({ stimulus: { ...clk, steady: true } }).stimulus, differential: false });
    // Mid-plateau times of cycles 25 and 26 (receiver edges every 0.5 ns after the
    // 0.5 ns line delay), away from the edges where the two time grids interpolate
    // differently.
    let e = 0;
    for (let i = 0; i < 4; i++) {
        const tt = st.t0 + 0.5e-9 + 0.25e-9 + 25e-9 + i * 0.5e-9;
        const ts = tt - st.t0 + sst.t0;
        const P = 1e-9, tsm = ((ts % P) + P) % P;
        e = Math.max(e, Math.abs(valueAt(burst, 'rx', tt) - valueAt(steady, 'rx', tsm)));
    }
    check('clock: steady state matches a long burst', e < 1e-6, `max diff ${e.toExponential(2)} V`);
    check('clock: eye folded', steady.eye && steady.eye.x.length > 10 && burst.eye && burst.eye.x.length > 10);
    // Source EMF of the steady state: a 0/1 V square wave with 50 % duty.
    const src = steady.traces.src;
    let mn = Infinity, mx = -Infinity, mean = 0;
    for (const v of src) { mn = Math.min(mn, v); mx = Math.max(mx, v); mean += v / src.length; }
    check('clock: steady source waveform', Math.abs(mn) < 1e-6 && Math.abs(mx - 1) < 1e-6 && Math.abs(mean - 0.5) < 2e-3,
        `min ${mn.toFixed(6)}, max ${mx.toFixed(6)}, mean ${mean.toFixed(4)}`);
    const duty = steady.metrics.list.find(m => m.label === 'Duty cycle');
    check('clock: RX duty cycle 50 %', duty && Math.abs(parseFloat(duty.value) - 50) < 0.1, duty && duty.value);
}

// 8. Symmetric pair vs its odd mode; transfer vs S-parameters.
{
    const odd = { R: 10, L: 300e-9, G: 0, C: 120e-12 }, even = { R: 15, L: 400e-9, G: 0, C: 90e-12 };
    const pair = pairTable(() => fromModes(odd, even));
    const Zo = Math.sqrt(odd.L / odd.C);
    const cfgD = base({ network: { rs: Zo, term: { type: 'diff', rdiff: 2 * Zo } } });
    const rd = simulateLink(pair, cfgD, fft);
    const se = simulateLink(seTable(() => odd), base({ network: { rs: Zo, term: { type: 'r', rt: Zo } } }), fft);
    let e = 0;
    const d0 = rd.traces.rxDiff[0], s0 = se.traces.rx[0];
    for (let i = 0; i < 200; i++) {
        const t = i * 10e-12;
        e = Math.max(e, Math.abs((valueAt(rd, 'rxDiff', t) - d0) - 2 * (valueAt(se, 'rx', t) - s0)));
    }
    check('symmetric pair: differential = 2 x odd mode', e < 1e-5, `max diff ${e.toExponential(2)} V`);
    let cm = 0;
    for (const v of rd.traces.rxCM) cm = Math.max(cm, Math.abs(v - rd.traces.rxCM[0]));
    check('symmetric pair: no common mode from a balanced drive', cm < 1e-6, `${cm.toExponential(2)} V`);

    // Asymmetric pair with 50 Ohm everywhere: H(rx_k, leg_l) = S(far k, near l) / 2.
    const asym = { R11: 8, R12: 2, R22: 12, L11: 320e-9, L12: 60e-9, L22: 280e-9, G11: 1e-3, G12: -2e-4, G22: 1.5e-3, C11: 110e-12, C12: -20e-12, C22: 125e-12 };
    const net = buildNetwork({ differential: true, rs: 50, tx: { ac: false }, rx: { ac: false }, term: { type: 'r', rt: 50 } });
    const tr = makeTransfer(net, 0.1, [net.nodes.d[0], net.nodes.d[1], net.nodes.a[0], net.nodes.a[1]]);
    let eS = 0;
    const m2 = (k) => [[asym[k + '11'], asym[k + '12']], [asym[k + '12'], asym[k + '22']]];
    for (const f of [1e6, 1e8, 1e9, 5e9, 2e10]) {
        const H = tr(f, asym);
        const sp = computeSParamsDifferentialMTL(f, m2('R'), m2('L'), m2('G'), m2('C'), 0.1, 50);
        // Rows: rxP, rxN, txP, txN; columns: leg P, leg N. S port order 1,2 near, 3,4 far.
        const want = [[2, 0], [2, 1], [3, 0], [3, 1], [0, 0], [0, 1], [1, 0], [1, 1]];
        want.forEach(([p, l], idx) => {
            const S = sp.S[p][l];
            // TX pad: V = (1 + S11) / 2 on the driven port, S/2 on the other.
            const exR = p < 2 ? ((p === l ? 1 : 0) + S.re) / 2 : S.re / 2, exI = S.im / 2;
            const row = p >= 2 ? p - 2 : p + 2;
            const i = row * 2 + l;
            eS = Math.max(eS, Math.hypot(H.re[i] - exR, H.im[i] - exI));
            void idx;
        });
    }
    check('asymmetric pair: network transfer = MTL S-parameters', eS < 1e-9, `max diff ${eS.toExponential(2)}`);
    const net1 = buildNetwork({ differential: false, rs: 50, tx: { ac: false }, rx: { ac: false }, term: { type: 'r', rt: 50 } });
    const tr1 = makeTransfer(net1, 0.3, [net1.nodes.d[0]]);
    let e1 = 0;
    for (const f of [0, 1e6, 1e9, 1e10, 1e11]) {
        const v = { R: 30, L: 260e-9, G: 2e-3, C: 105e-12 };
        const H = tr1(f, v), S = computeSParamsSingleEnded(f, v, 0.3, 50).S21;
        e1 = Math.max(e1, Math.hypot(H.re[0] - S.re / 2, H.im[0] - S.im / 2));
    }
    check('single-ended: network transfer = S21 / 2', e1 < 1e-9, `max diff ${e1.toExponential(2)}`);
}

// 9. Coverage rule.
{
    const rows = (fs) => fs.map(f => ({ freq: f }));
    const lin = (a, b, n) => Array.from({ length: n }, (_, i) => a + (b - a) * i / (n - 1));
    check('coverage: no DC point', !checkCoverage(rows([1e8, 1e9, 1e10]), 5e9, true).ok);
    check('coverage: stop too low', !checkCoverage(rows([0, 1e8, 1e9]), 5e9, true).ok);
    check('coverage: interpolating sweep from DC', checkCoverage(rows([0, 1e8, 2e10]), 1e10, true).ok);
    check('coverage: coarse linear sweep rejected', !checkCoverage(rows(lin(0, 2e10, 10)), 1e10, false).ok);
    check('coverage: dense linear sweep accepted', checkCoverage(rows(lin(0, 2e10, 201)), 1e10, false).ok);
}

// 10. FFT and sampling.
{
    for (const n of [30, 64, 360, 1000]) {
        const x = Float64Array.from({ length: n }, (_, i) => Math.sin(0.3 * i) + 0.01 * i * i / n);
        const X = fft.rfft(x);
        let e = 0;
        for (const k of [0, 1, 7, n / 2]) {
            let re = 0, im = 0;
            for (let i = 0; i < n; i++) { re += x[i] * Math.cos(2 * Math.PI * k * i / n); im -= x[i] * Math.sin(2 * Math.PI * k * i / n); }
            e = Math.max(e, Math.hypot(X[2 * k] - re, X[2 * k + 1] - im));
        }
        const y = fft.irfftBatch(n, 1, X);
        let er = 0;
        for (let i = 0; i < n; i++) er = Math.max(er, Math.abs(y[i] - x[i]));
        check(`FFT n=${n}: vs DFT and round trip`, e < 1e-9 && er < 1e-12, `${e.toExponential(1)} / ${er.toExponential(1)}`);
    }
    check('smoothSize', smoothSize(1001) === 1024 && smoothSize(241) === 250 && smoothSize(7) === 8);
    const lossy = seTable(f => ({ R: 5 + 1e-4 * Math.sqrt(f), L: L0, G: 1e-12 * f, C: C0 }));
    const a = simulateLink(lossy, base({ network: { rs: 40 }, stimulus: { type: 'pulse', width: 0.3e-9 } }), fft);
    const b = simulateLink(lossy, base({ network: { rs: 40 }, stimulus: { type: 'pulse', width: 0.3e-9 }, samplesPerEdge: 32 }), fft);
    let e = 0;
    for (let i = 0; i < 300; i++) { const t = i * 5e-12; e = Math.max(e, Math.abs(valueAt(a, 'rx', t) - valueAt(b, 'rx', t))); }
    check('doubling the sampling rate changes the pulse < 1e-3 of its swing', e < 5e-4, `${e.toExponential(2)} V`);
}

// 11. Constant loss tangent (G ~ f, C constant) is non-causal: a precursor warning.
{
    const ncausal = seTable(f => ({ R: 0, L: L0, G: 2 * Math.PI * f * C0 * 0.05, C: C0 }));
    const r = simulateLink(ncausal, { ...base({ causal: false }), length: 0.5 }, fft);
    check('non-causal dielectric: precursor warning', r.warnings.some(w => w.includes('before the signal can arrive')), r.warnings.join(' | '));
    const r2 = simulateLink(lossless, base(), fft);
    check('lossless line: no precursor warning', !r2.warnings.some(w => w.includes('before the signal')));
}

// 12. 0.5 m skin-effect line, 5 ps edge, 5 Ohm source, open end: heavily attenuated at
// the edge frequencies but ringing for tens of round trips at low frequencies.
{
    // Causal skin effect: equal resistance and internal reactance, no dielectric loss.
    const skin = seTable(f => {
        const ff = Math.max(f, 1e6);
        return { R: 3 + 2e-4 * Math.sqrt(ff), L: L0 + 2e-4 * Math.sqrt(ff) / (2 * Math.PI * ff), G: 0, C: C0 };
    }, DENSE);
    for (const rs of [20, 5]) {
        const cfg = { ...base({ network: { rs, term: { type: 'open' } }, stimulus: { tr: 5e-12, tf: 5e-12, vHigh: 3.3 } }), length: 0.5 };
        const r = simulateLink(skin, cfg, fft);
        const st = buildStimulus({ ...cfg.stimulus, differential: false });
        const tArrive = st.t0 + 0.95 * 0.5 * Math.sqrt(L0 * C0);
        let pre = 0;
        for (let i = 0; i < r.t.length && r.t[i] < tArrive; i++) pre = Math.max(pre, Math.abs(r.traces.rx[i]));
        check(`open end, ${rs} Ohm source: nothing before the arrival`, pre < 3.3e-4 && r.warnings.length === 0,
            `max ${pre.toExponential(2)} V of 3.3 V, window ${(r.info.window * 1e9).toFixed(0)} ns${r.warnings.length ? ', ' + r.warnings.join(' | ') : ''}`);
    }
}

// 13. Low-frequency notes.
{
    const sig = 5.8e7, mu0 = 4e-7 * Math.PI, w = 0.35e-3, t = 35e-6;
    const skin = seTable(f => {
        const Rdc = 1 / (sig * w * t), Rac = f > 0 ? Math.sqrt(Math.PI * f * mu0 / sig) / w : 0;
        const R = Math.hypot(Rdc, Rac);
        return { R, L: 300e-9 + (f > 0 ? Rac / (2 * Math.PI * f) : 0), G: 0, C: 120e-12 };
    }, [0, ...Array.from({ length: 141 }, (_, i) => 1e3 * Math.pow(10, i / 20))]);
    skin.lext = { L: 300e-9 };
    const notes = lowFrequencyNotes([
        { reason: 'narrow-gap', gap: 300e-6, sigma: sig }, { reason: 'narrow-gap', gap: 200e-6, sigma: sig },
        { reason: 'skin-transition', fBand: 14e6 }, { reason: 'skin-transition', fBand: 12e6 },
        { reason: 'ground-spreading', fOn: 594e3, gap: 210e-6 }, { reason: 'open-boundary' }]);
    check('notes merged: narrowest gap, highest skin band, one spreading note',
        notes.length === 3 && notes[0].gap === 200e-6 && notes[1].fBand === 14e6);
    const mk = (len, r, stimulus) => ({ length: len, samplesPerEdge: 16,
        network: { rs: r, tx: { ac: false }, rx: { ac: false }, term: { type: 'r', rt: r } },
        stimulus: { vLow: 0, vHigh: 1, edgeDef: '10-90', shape: 'rcos', delay: 0, ...stimulus } });
    const ghz = mk(0.1, 50, { type: 'pulse', width: 1e-9, tr: 35e-12, tf: 35e-12 });
    const lfG = lowFrequencyImpact(skin, ghz, simulateLink(skin, ghz, fft), notes);
    check('GHz pulse on a terminated line: notes negligible', lfG.every(x => x.impact < 1e-4),
        lfG.map(x => `${x.reason} ${(100 * x.impact).toExponential(1)} %`).join(', '));
    const dc = makeRLGCEvaluator(skin)(0);
    for (const fc of [2e5, 1e7]) {
        const cfg = mk(1, 5, { type: 'clock', clockFreq: fc, duty: 0.5, steady: true, tr: fc > 1e6 ? 1e-9 : 10e-9, tf: fc > 1e6 ? 1e-9 : 10e-9 });
        const base = simulateLink(skin, cfg, fft);
        const lf = lowFrequencyImpact(skin, cfg, base, notes);
        notes.forEach((note, k) => {
            const p = simulateLink(skin, { ...cfg, perturb: perturbation(note, dc, skin.lext) }, fft);
            let e = 0, lo = Infinity, hi = -Infinity;
            for (let i = 0; i < base.t.length; i++) {
                e = Math.max(e, Math.abs(base.traces.rx[i] - p.traces.rx[i]));
                lo = Math.min(lo, base.traces.rx[i]); hi = Math.max(hi, base.traces.rx[i]);
            }
            const actual = e / (hi - lo), bound = lf[k].impact;
            check(`${fc / 1e3} kHz clock, 1 m, 5 Ohm, ${note.reason}: bound >= actual, within 3x`,
                bound >= actual * (1 - 1e-6) && (actual === 0 ? bound < 1e-9 : bound <= 3 * actual),
                `actual ${(100 * actual).toPrecision(3)} %, bound ${(100 * bound).toPrecision(3)} %`);
        });
        if (fc === 2e5) check('slow clock into 5 Ohm on 1 m: ground spreading matters (> 1 %)', lf[2].impact > 0.01, `${(100 * lf[2].impact).toPrecision(3)} %`);
    }
}

// 14. Bit patterns.
{
    for (const name of ['prbs7', 'prbs9', 'prbs11', 'prbs15']) {
        const P = prbsPeriod(name), n = PRBS[name][0], b = prbsBits(name, 2 * P);
        let ones = 0, periodic = true, shorter = false, run = 0, r1 = 0;
        for (let i = 0; i < P; i++) { ones += b[i]; if (b[i] !== b[i + P]) periodic = false; }
        for (const q of [1, 3, 7, P / 7, P / 31].filter(q => Number.isInteger(q) && q < P)) {
            let same = true;
            for (let i = 0; i < P; i++) if (b[i] !== b[(i + q) % P]) { same = false; break; }
            if (same) shorter = true;
        }
        for (let i = 0; i < 2 * P; i++) { run = i && b[i] === b[i - 1] ? run + 1 : 1; if (b[i]) r1 = Math.max(r1, run); }
        check(`${name}: period 2^${n}-1, ${2 ** (n - 1)} ones, longest run ${n}`, periodic && !shorter && ones === 2 ** (n - 1) && r1 === n,
            `ones ${ones}, longest run ${r1}`);
    }
    const lossy = seTable(() => ({ R: 20, L: L0, G: 0, C: C0 }));
    const pcfg = (stim, net = {}, len = 0.1) => ({ length: len, samplesPerEdge: 16,
        network: { rs: 30, tx: { ac: false }, rx: { ac: false }, term: { type: 'r', rt: 80 }, ...net },
        stimulus: { vLow: 0, vHigh: 1, tr: 30e-12, tf: 30e-12, edgeDef: '10-90', shape: 'rcos', delay: 0, ...stim } });
    // Pattern at 2 Gb/s vs the clock steady state, compared at mid-bit times of the
    // receiver (the line delay is 0.5 ns, one UI).
    const atPeriodic = (r, name, t, P) => {
        const n = r.t.length, dtc = r.t[1] - r.t[0], x = (((t % P) + P) % P) / dtc, i = Math.floor(x), f = x - i;
        const v = r.traces[name], m = Math.round(P / dtc);
        return v[i % m] + (v[(i + 1) % m] - v[i % m]) * f;
    };
    for (const [bits, fClk, duty] of [['01', 1e9, 0.5], ['0011', 0.5e9, 0.5], ['0111', 0.5e9, 0.75]]) {
        for (const ac of [false, true]) {
            const net = ac ? { tx: { ac: true, c: 1e-9 }, rx: { ac: true, c: 1e-9 } } : {};
            const rp = simulateLink(lossy, pcfg({ type: 'pattern', bitRate: 2e9, bits: patternBits({ pattern: 'custom', bits }) }, net), fft);
            const rc = simulateLink(lossy, pcfg({ type: 'clock', clockFreq: fClk, duty, steady: true }, net), fft);
            const sc = buildStimulus({ ...pcfg({ type: 'clock', clockFreq: fClk, duty, steady: true }).stimulus, differential: false });
            const ui = 0.5e-9, tRise = (bits.indexOf('01') + 1) * ui, P = 1 / fClk;
            let e = 0;
            // The pattern result holds one period of the pattern.
            for (let k = 0; k < bits.length; k++) {
                const t = (k + 0.5) * ui;
                e = Math.max(e, Math.abs(atPeriodic(rp, 'rx', t, bits.length * ui) - atPeriodic(rc, 'rx', t - tRise + sc.t0, P)));
            }
            check(`pattern ${bits} = clock steady state (${ac ? 'AC' : 'DC'} coupled), mid-bit`, e < 2e-4, `max diff ${e.toExponential(2)} V`);
        }
    }
    // Pair with negative skew: pattern '01' vs the clock.
    {
        const odd = { R: 10, L: 300e-9, G: 0, C: 120e-12 }, even = { R: 15, L: 400e-9, G: 0, C: 90e-12 };
        const pair = pairTable(() => fromModes(odd, even));
        const net = { term: { type: 'diff', rdiff: 90 } };
        const rp = simulateLink(pair, pcfg({ type: 'pattern', bitRate: 2e9, skew: -40e-12, bits: patternBits({ pattern: 'custom', bits: '01' }) }, net), fft);
        const rc = simulateLink(pair, pcfg({ type: 'clock', clockFreq: 1e9, duty: 0.5, steady: true, skew: -40e-12 }, net), fft);
        const sc = buildStimulus({ ...pcfg({ type: 'clock', clockFreq: 1e9, duty: 0.5, steady: true, skew: -40e-12 }).stimulus, differential: true });
        const sp = buildStimulus({ ...pcfg({ type: 'pattern', bitRate: 2e9, skew: -40e-12, bits: Uint8Array.of(0, 1) }).stimulus, differential: true });
        let e = 0;
        for (const name of ['rxDiff', 'rxCM']) for (let k = 0; k < 2; k++) {
            const t = (k + 0.5) * 0.5e-9;
            // The pattern's P leg rises at 0.5 ns (bit 1), the clock's at sc.t0.
            e = Math.max(e, Math.abs(atPeriodic(rp, name, t, 1e-9) - atPeriodic(rc, name, t - 0.5e-9 + sc.t0, 1e-9)));
        }
        void sp;
        check('pair pattern with -40 ps skew = clock steady state (diff and common mode)', e < 2e-4, `max diff ${e.toExponential(2)} V`);
    }
    // Lossless matched line: a fully open eye, no jitter, the line delay.
    {
        const r = simulateLink(lossless, pcfg({ type: 'pattern', bitRate: 10e9, bits: patternBits({ pattern: 'prbs7' }) }, { rs: 50, term: { type: 'r', rt: 50 } }), fft);
        const m = Object.fromEntries(r.metrics.list.map(x => [x.label, x.value]));
        check('lossless PRBS7: eye height = amplitude, width = UI, no jitter, 500 ps delay',
            m['Eye height (centre)'].startsWith('500 mV') && m['Eye width (at threshold)'].startsWith('100 ps')
            && m['Crossing jitter pk-pk / rms'] === '0 ps / 0 ps' && m['Delay (50 %)'] === '500 ps',
            JSON.stringify([m['Eye height (centre)'], m['Eye width (at threshold)'], m['Crossing jitter pk-pk / rms'], m['Delay (50 %)']]));
        check('eye histogram filled', r.eyeHist && r.eyeHist.counts.some(c => c > 0));
        // Every time column carries one unit per bit: no brighter columns from the
        // sampling beating against the bins.
        const h = r.eyeHist, col = new Float64Array(h.tBins);
        for (let iy = 0; iy < h.vBins; iy++) for (let ix = 0; ix < h.tBins; ix++) col[ix] += h.counts[iy * h.tBins + ix];
        const cmin = Math.min(...col), cmax = Math.max(...col);
        check('eye histogram: equal column totals (= bits)', Math.abs(cmin - 127) < 1e-3 && Math.abs(cmax - 127) < 1e-3, `${cmin.toFixed(4)} .. ${cmax.toFixed(4)}`);
    }
    // AC coupling: an unbalanced pattern settles around its mean (baseline wander), the
    // mean of the receiver is zero; DC coupled the mean is the divided mean level.
    {
        const bits = patternBits({ pattern: 'custom', bits: '0000000001' });
        const ac = simulateLink(lossy, pcfg({ type: 'pattern', bitRate: 2e9, bits }, { tx: { ac: true, c: 1e-9 }, rx: { ac: true, c: 1e-9 } }), fft);
        const dc = simulateLink(lossy, pcfg({ type: 'pattern', bitRate: 2e9, bits }), fft);
        const mean = (v) => v.reduce((a, x) => a + x, 0) / v.length;
        const mAc = mean(ac.traces.rx), mDc = mean(dc.traces.rx), want = 0.1 * 80 / (30 + 2 + 80);
        check('AC-coupled unbalanced pattern: zero-mean receiver; DC coupled: divided mean', Math.abs(mAc) < 1e-6 && Math.abs(mDc - want) < 1e-4,
            `AC mean ${mAc.toExponential(2)} V, DC mean ${mDc.toFixed(5)} V (want ${want.toFixed(5)})`);
    }
}

// 15. Drivers, pads, terminations, equalization.
{
    const odd = { R: 5, L: 300e-9, G: 0, C: 120e-12 }, even = { R: 8, L: 400e-9, G: 0, C: 90e-12 };
    const pair = pairTable(() => fromModes(odd, even));
    // LVDS: +-3.5 mA steered through 100 Ohm, 10 kOhm output resistance to 1.2 V.
    {
        const cfg = { length: 0.1, samplesPerEdge: 16,
            network: { driver: 'current', rout: 1e4, vcm: 1.2, tx: { ac: false }, rx: { ac: false }, term: { type: 'diff', rdiff: 100 } },
            stimulus: { type: 'step', vLow: -3.5e-3, vHigh: 3.5e-3, tr: 200e-12, tf: 200e-12, edgeDef: '20-80', shape: 'rcos', delay: 0 }, span: 20e-9 };
        const r = simulateLink(pair, cfg, fft);
        const d = r.traces.rxDiff, cm = r.traces.rxCM;
        // Odd half circuit: 3.5 mA into 10 kOhm || (0.5 Ohm line + 50 Ohm), the far end
        // across the 50 Ohm half of the termination.
        const rl = 0.5, leg = 3.5e-3 * (1e4 * (rl + 50)) / (1e4 + rl + 50) * 50 / (rl + 50), want = 2 * leg;
        check('LVDS current mode: +-348 mV differential, 1.2 V common mode',
            Math.abs(d[d.length - 1] - want) < 1e-4 && Math.abs(d[0] + want) < 1e-4 && Math.abs(cm[cm.length - 1] - 1.2) < 1e-6,
            `diff ${(d[0] * 1e3).toFixed(2)} -> ${(d[d.length - 1] * 1e3).toFixed(2)} mV (want ${(want * 1e3).toFixed(2)}), CM ${cm[cm.length - 1].toFixed(5)} V`);
    }
    // Thevenin termination: 130 / 82 Ohm from 3.3 V, driver 7 Ohm, EMF 1.6 -> 2.4 V.
    {
        const lossy = seTable(() => ({ R: 2, L: L0, G: 0, C: C0 }));
        const cfg = { length: 0.1, samplesPerEdge: 16, span: 10e-9,
            network: { rs: 7, tx: { ac: false }, rx: { ac: false }, term: { type: 'thev', rup: 130, rdn: 82, vsup: 3.3 } },
            stimulus: { type: 'step', vLow: 1.6, vHigh: 2.4, tr: 200e-12, tf: 200e-12, edgeDef: '20-80', shape: 'rcos', delay: 0 } };
        const r = simulateLink(lossy, cfg, fft), v = r.traces.rx;
        const rth = 130 * 82 / 212, vth = 3.3 * 82 / 212, rser = 7 + 0.2;
        const lvl = (e) => (e * rth + vth * rser) / (rth + rser);
        check('Thevenin termination: DC levels', Math.abs(v[0] - lvl(1.6)) < 1e-6 && Math.abs(v[v.length - 1] - lvl(2.4)) < 1e-4,
            `${v[0].toFixed(5)} / ${v[v.length - 1].toFixed(5)} V (want ${lvl(1.6).toFixed(5)} / ${lvl(2.4).toFixed(5)})`);
    }
    // RX pad capacitance on a 1 mm line, open end: v = 1 - exp(-t / (Rs C)).
    {
        const cfg = { length: 1e-3, samplesPerEdge: 16, span: 3e-9,
            network: { rs: 50, padRx: 10e-12, tx: { ac: false }, rx: { ac: false }, term: { type: 'open' } },
            stimulus: { type: 'step', vLow: 0, vHigh: 1, tr: 10e-12, tf: 10e-12, edgeDef: '10-90', shape: 'gauss', delay: 0 } };
        const r = simulateLink(lossless, cfg, fft);
        const st = buildStimulus({ ...cfg.stimulus, differential: false });
        // The 1 mm line adds its own 0.1 pF and 5 ps of delay; its 0.25 nH make the RC a
        // second-order circuit, up to L / (R^2 C) = 1 % off the RC curve.
        const tau = 50 * (10e-12 + C0 * 1e-3), td = 1e-3 * Math.sqrt(L0 * C0);
        const e = maxErr(r, 'rx', t => t < st.t0 + td ? 0 : 1 - Math.exp(-(t - st.t0 - td) / tau), st.t0 + 50e-12, st.t0 + 2.5e-9);
        check('RX pad capacitance: RC charging through the source resistance', e < 5e-3, `max err ${e.toExponential(2)} V`);
    }
    // PCIe preset taps (published coefficients) and their dB round trip.
    {
        const want = { P0: [0, 0.75, -0.25], P1: [0, 0.833, -0.167], P4: [0, 1, 0], P5: [-0.1, 0.9, 0], P7: [-0.1, 0.7, -0.2], P8: [-0.125, 0.75, -0.125], P9: [-0.167, 0.833, 0] };
        let e = 0, rt = 0;
        for (const [p, w] of Object.entries(want)) {
            const c = firFromPresetDb(...PCIE_PRESETS[p]);
            e = Math.max(e, ...c.map((x, i) => Math.abs(x - w[i])));
            const db = presetDbFromFir(c);
            rt = Math.max(rt, Math.abs(db.preshoot - PCIE_PRESETS[p][0]), Math.abs(db.deemphasis - PCIE_PRESETS[p][1]));
        }
        check('PCIe presets: taps within 0.003 of the published values, dB round trip', e < 3e-3 && rt < 1e-9, `${e.toExponential(1)}, ${rt.toExponential(1)} dB`);
        // P10: no pre-cursor (Va / Vd = 1) and 9.5 dB nominal maximum boost.
        const p10 = firFromPresetDb(...PCIE_PRESETS.P10);
        check('PCIe P10: c-1 = 0, boost 9.5 dB', p10[0] === 0 && Math.abs(boostDb(p10) - 9.5) < 1e-9, `boost ${boostDb(p10).toFixed(3)} dB`);
    }
    // FIR source levels: P7 on PRBS7, the source takes the levels 1 - 2a, 1 - 2b,
    // 1 - 2a - 2b and 1 (times the half swing).
    {
        const fir = firFromPresetDb(...PCIE_PRESETS.P7);
        const cfg = { length: 0.1, samplesPerEdge: 16, network: { rs: 50, tx: { ac: false }, rx: { ac: false }, term: { type: 'r', rt: 50 } },
            stimulus: { type: 'pattern', bitRate: 5e9, bits: patternBits({ pattern: 'prbs7' }), fir, vLow: 0, vHigh: 1, tr: 30e-12, tf: 30e-12, edgeDef: '10-90', shape: 'rcos' } };
        const r = simulateLink(lossless, cfg, fft);
        const ui = 200e-12, a = -fir[0], b = -fir[2], want = [1, 1 - 2 * a, 1 - 2 * b, 1 - 2 * a - 2 * b];
        const lv = new Set();
        for (let k = 0; k < 127; k++) lv.add((2 * valueAt(r, 'src', (k + 0.5) * ui) - 1).toFixed(4));
        const got = [...lv].map(Number).map(Math.abs);
        const ok = got.every(g => want.some(w => Math.abs(g - w) < 1e-4)) && want.every(w => got.some(g => Math.abs(g - w) < 1e-4));
        check('TX FIR P7: source levels', ok, [...lv].sort().join(' '));
    }
    // CTLE: 0 dB with poles far above the signal leaves it; the DC gain scales the mean.
    {
        const lossy = seTable(() => ({ R: 20, L: L0, G: 0, C: C0 }));
        const base = { length: 0.1, samplesPerEdge: 16, network: { rs: 50, tx: { ac: false }, rx: { ac: false }, term: { type: 'r', rt: 50 } },
            stimulus: { type: 'pattern', bitRate: 2e9, bits: patternBits({ pattern: 'prbs7' }), vLow: 0, vHigh: 1, tr: 50e-12, tf: 50e-12, edgeDef: '10-90', shape: 'rcos' } };
        const r0 = simulateLink(lossy, { ...base, rxEq: { ctle: { adcDb: 0, fp1: 1e13, fp2: 2e13 }, dfeTaps: 0 } }, fft);
        let e = 0;
        for (let i = 0; i < r0.t.length; i++) e = Math.max(e, Math.abs(r0.traces.rxEq[i] - r0.traces.rx[i]));
        check('CTLE at 0 dB with far poles: equalized = RX', e < 2e-3, `max diff ${e.toExponential(2)} V`);
        const r6 = simulateLink(lossy, { ...base, rxEq: { ctle: { adcDb: -6, fp1: 2e9, fp2: 8e9 }, dfeTaps: 0 } }, fft);
        const mean = (v) => v.reduce((s, x) => s + x, 0) / v.length;
        // Over the shown 127-bit period the means relate by the DC gain.
        const ratio = mean(r6.traces.rxEq) / mean(r6.traces.rx);
        check('CTLE: DC gain on the mean level', Math.abs(ratio - Math.pow(10, -6 / 20)) < 1e-6, ratio.toFixed(6));
        const r10 = simulateLink(lossy, { ...base, rxEq: { ctle: { order: 2, adcDb: -10, fz1: 450e6, fp2: 9.5e9, fp3: 28e9 }, dfeTaps: 0 } }, fft);
        const ratio2 = mean(r10.traces.rxEq) / mean(r10.traces.rx);
        check('2nd-order CTLE: DC gain on the mean level', Math.abs(ratio2 - Math.pow(10, -10 / 20)) < 1e-6, ratio2.toFixed(6));
    }
    // 2nd-order CTLE (PCIe 32 / 64 GT/s) against the formula in magnitude and phase form:
    // |H| = k prod |s + z| / prod |s + p|, arg H = sum atan(w / z) - sum atan(w / p).
    {
        const tw = (f) => 2 * Math.PI * f, A = Math.pow(10, -9 / 20);
        const wz1 = tw(450e6), wp1 = tw(1.65 * 450e6), wp2 = tw(9.5e9), wp3 = tw(28e9);
        const H = ctle2Transfer(-9, 450e6, 9.5e9, 28e9);
        let e = 0;
        for (const f of [0, 1e8, 7e8, 3e9, 1.6e10, 6e10]) {
            const w = tw(f), m = (a) => Math.hypot(a, w), ph = (a) => Math.atan2(w, a);
            const mag = (wp1 * wp3 * wp3 / wz1) * m(wz1) * m(A * wp2) / (m(wp1) * m(wp2) * m(wp3) * m(wp3));
            const arg = ph(wz1) + ph(A * wp2) - ph(wp1) - ph(wp2) - 2 * ph(wp3);
            const [re, im] = H(f);
            e = Math.max(e, Math.hypot(re - mag * Math.cos(arg), im - mag * Math.sin(arg)) / mag);
        }
        const [dc] = H(0), peak = Math.max(...Array.from({ length: 400 }, (_, i) => Math.hypot(...H(1e8 * Math.pow(1.02, i)))));
        check('2nd-order CTLE = PCIe 32 GT/s formula, DC gain A', e < 1e-12 && Math.abs(dc - A) < 1e-12,
            `rel diff ${e.toExponential(1)}, DC ${dc.toFixed(6)}, peaking ${(20 * Math.log10(peak / dc)).toFixed(2)} dB over DC`);
    }
    // DFE and automatic choices on a lossy pair at 10 Gb/s.
    {
        const skinPair = pairTable(f => {
            const ff = Math.max(f, 1e6), r = 5 + 4e-4 * Math.sqrt(ff), li = 4e-4 * Math.sqrt(ff) / (2 * Math.PI * ff);
            return { R11: r, R12: 0.1 * r, R22: r, L11: 330e-9 + li, L12: 50e-9, L22: 330e-9 + li, G11: 0, G12: 0, G22: 0, C11: 110e-12, C12: -10e-12, C22: 110e-12 };
        });
        const base = { length: 0.5, samplesPerEdge: 16, network: { rs: 50, tx: { ac: true, c: 220e-9 }, rx: { ac: false }, term: { type: 'r', rt: 50 } },
            stimulus: { type: 'pattern', bitRate: 10e9, bits: patternBits({ pattern: 'prbs9' }), vLow: 0, vHigh: 1, tr: 25e-12, tf: 25e-12, edgeDef: '20-80', shape: 'rcos' } };
        const height = (r) => { const m = r.metrics.list.find(x => x.label === 'Eye height (centre)'); return m.value === 'closed' ? 0 : parseFloat(m.value) * (m.value.includes(' mV') ? 1e-3 : 1); };
        const none = simulateLink(skinPair, base, fft);
        const dfe = simulateLink(skinPair, { ...base, rxEq: { dfeTaps: 4 } }, fft);
        const auto = simulateLink(skinPair, { ...base, txEq: { auto: true }, rxEq: { ctle: { auto: true, adcMin: -12, adcMax: 0, fp1: 2.5e9, fp2: 10e9 }, dfeTaps: 2 } }, fft);
        check('DFE (4 taps) opens the eye beyond no equalization', height(dfe) > height(none), `${(height(none) * 1e3).toFixed(1)} -> ${(height(dfe) * 1e3).toFixed(1)} mV`);
        check('the eye is at least the peak-distortion estimate', height(dfe) >= dfe.eqChoice.pdaHeight - 1e-6 && height(auto) >= auto.eqChoice.pdaHeight - 1e-6,
            `DFE ${(height(dfe) * 1e3).toFixed(1)} >= ${(dfe.eqChoice.pdaHeight * 1e3).toFixed(1)}, auto ${(height(auto) * 1e3).toFixed(1)} >= ${(auto.eqChoice.pdaHeight * 1e3).toFixed(1)} mV`);
        check('automatic preset + CTLE + DFE beats no equalization', height(auto) > height(none),
            `${auto.eqChoice.firName}, CTLE ${auto.eqChoice.adc} dB: ${(height(auto) * 1e3).toFixed(1)} mV`);
        // DFE tap limits: a tap stops at its limit, and the rest of its post-cursor counts
        // in the eye estimate. The 4-tap DFE above, with limits below its taps.
        const t = dfe.eqChoice.dfe, lim = t.map(x => 0.5 * Math.abs(x));
        const lim2 = simulateLink(skinPair, { ...base, rxEq: { dfeTaps: 4, dfeLimits: lim } }, fft);
        const clamped = lim2.eqChoice.dfe.every((x, k) => Math.abs(Math.abs(x) - lim[k]) < 1e-12);
        check('DFE tap limits: taps clamped, eye smaller than unlimited, at least the estimate',
            clamped && height(lim2) < height(dfe) && height(lim2) >= lim2.eqChoice.pdaHeight - 1e-6,
            `taps ${lim2.eqChoice.dfe.map(x => (x * 1e3).toFixed(1)).join('/')} mV (limits ${lim.map(x => (x * 1e3).toFixed(1)).join('/')}), `
            + `eye ${(height(lim2) * 1e3).toFixed(1)} vs ${(height(dfe) * 1e3).toFixed(1)} mV, estimate ${(lim2.eqChoice.pdaHeight * 1e3).toFixed(1)} mV`);
        // A limit relative to the cursor (PCIe 32 GT/s: |d1| <= 0.8 h0): same taps as the
        // equivalent voltage limit.
        const h0 = lim2.pulse.cursors.find(c => c.k === 0).value;
        const rel = simulateLink(skinPair, { ...base, rxEq: { dfeTaps: 4, dfeLimits: lim.map(x => ({ rel: x / h0 })) } }, fft);
        const dRel = Math.max(...rel.eqChoice.dfe.map((x, k) => Math.abs(x - lim2.eqChoice.dfe[k])));
        check('DFE tap limit as a ratio to the cursor = the same limit in volts', dRel < 1e-12 && Math.abs(height(rel) - height(lim2)) < 1e-9,
            `tap diff ${dRel.toExponential(1)} V`);
        // Reduced swing: only P4, P1, P9, P5, P6, P3. A lossier line at 16 Gb/s, where full
        // swing picks a preset reduced swing lacks.
        const lossy = { ...base, length: 3.0, stimulus: { ...base.stimulus, bitRate: 16e9, tr: 15e-12, tf: 15e-12 } };
        const full = simulateLink(skinPair, { ...lossy, txEq: { auto: true }, rxEq: { dfeTaps: 0, ctle: null } }, fft);
        const rs = simulateLink(skinPair, { ...lossy, txEq: { auto: true, presets: REDUCED_SWING_PRESETS }, rxEq: { dfeTaps: 0, ctle: null } }, fft);
        check('reduced-swing preset search stays within P4, P1, P9, P5, P6, P3 and never beats full swing',
            REDUCED_SWING_PRESETS.includes(rs.eqChoice.firName) && !REDUCED_SWING_PRESETS.includes(full.eqChoice.firName)
            && rs.eqChoice.pdaHeight < full.eqChoice.pdaHeight,
            `full swing ${full.eqChoice.firName} ${(full.eqChoice.pdaHeight * 1e3).toFixed(1)} mV, reduced ${rs.eqChoice.firName} ${(rs.eqChoice.pdaHeight * 1e3).toFixed(1)} mV`);
    }
}

void makeRLGCEvaluator;
done();
