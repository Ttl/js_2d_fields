// Statistical eye (src/eye/stat_eye.js).
//
//   1. The chain against brute-force enumeration of every bit sequence and dual-Dirac
//      jitter combination: 3-tap FIR, unequal rise and fall responses, a DFE tap, both
//      conditional distributions at two sampling phases.
//   2. Lossless matched line with Rj, DjDD and PW DDj: the eye width at the threshold
//      against the analytic BER(t) = P(delta > t) / 2 + P(delta < t - UI) / 2, and the
//      eye height (the full amplitude), the mean UI centre in the middle of the width.
//   3. Short-memory lossy line without jitter (few cursors, so the worst sequence is more
//      likely than 1e-12): eye height at the chosen sample point = the peak-distortion worst
//      case = the PRBS15 time-domain eye there, with and without DFE.
//   5. Equalization search: the TX preset / CTLE gain / sample point chosen is the argmax of
//      the full-grid eye height times width (the default search and every pair on the full
//      grid); without the statistical eye the peak-distortion argmax, a different choice here.
//   4. Baseline wander of a slow AC-coupling mode: the passband level times the rms of a
//      first-order low-pass of random +-1 bits, sqrt((1 - a) / (1 + a)), a = exp(-UI / lambda).
import { check, done } from './helpers.js';
import { loadFFT } from '../src/eye/fft.js';
import { simulateLink } from '../src/eye/transient.js';
import { erfc } from '../src/eye/stimulus.js';
import { patternBits } from '../src/eye/pattern.js';
import { makeStatEye } from '../src/eye/stat_eye.js';
import { PACKAGES } from '../src/eye/network.js';

const fft = await loadFFT();
const FREQS = [0, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12];
function seTable(fn) {
    const data = { R: [], L: [], G: [], C: [] };
    for (const f of FREQS) { const v = fn(f); for (const k in data) data[k].push(v[k]); }
    return { differential: false, freqs: Float64Array.from(FREQS), keys: ['R', 'L', 'G', 'C'],
        data: Object.fromEntries(Object.entries(data).map(([k, v]) => [k, Float64Array.from(v)])),
        lext: { L: fn(FREQS[FREQS.length - 1]).L } };
}
const L0 = 250e-9, C0 = 100e-12;
const metric = (r, label) => { const m = r.metrics.list.find(x => x.label === label); return m ? m.value : undefined; };
const volts = (s) => parseFloat(s) * (/ mV/.test(s) ? 1e-3 : / µV| uV/.test(s) ? 1e-6 : 1);

// 1. Brute force.
{
    const spb = 16, N = 160, dt = 1 / spb;
    const edge = (t, s) => 0.5 * erfc(-t / (s * Math.SQRT2));
    // An edge at index 32, a slow tail and a negative reflection 2.4 UI later.
    const mk = (s) => Float64Array.from({ length: N }, (_, n) => {
        const t = (n - 32) * dt;
        return 0.4 * edge(t, s) + 0.08 * edge(t - 0.6, 0.4) - 0.06 * edge(t - 2.4, s);
    });
    const gUp = mk(0.12), gDn = mk(0.2);
    const fir = [-0.1, 0.7, -0.2], dfe = [0.03];
    const p = { gUp, gDn, fir, spb, dt, n0: 41, dfe, jitter: { djdd: 0.1, pwDj: 0.08 }, vmid: 0.01, sigmaSlow: 0 };
    const nodes = [-0.07, -0.03, 0.03, 0.07];
    const se = makeStatEye(p);
    const at = (a, x) => { if (x <= 0) return 0; if (x >= N - 1) return a[N - 1]; const i = Math.floor(x), f = x - i; return a[i] + (a[i + 1] - a[i]) * f; };
    const lvl = (d, b) => fir[0] * d(b + 1) + fir[1] * d(b) + fir[2] * d(b - 1);
    const nV = 1201, bin = 1.8 / (nV - 1), grid = { v0: p.vmid - (nV - 1) / 2 * bin, bin, nV };
    let worst = 0;
    const rows = [];
    for (const phase of [0, 6]) {
        const s = p.n0 + phase;
        // Bits -8..3, jitter on boundaries -6..2 (the others settled or not arrived);
        // older bits held at the settled level of the oldest one.
        const mLo = -8, mHi = 3, nb = mHi - mLo + 1, bLo = -6, bHi = 2;
        const exact = { 1: [], [-1]: [] };
        for (let code = 0; code < (1 << nb); code++) {
            const d = (m) => (m < mLo || m > mHi) ? 0 : ((code >> (m - mLo)) & 1) ? 1 : -1;
            let base = p.vmid + gUp[N - 1] * lvl(d, mLo + 1);
            for (let k = 1; k <= dfe.length; k++) base -= dfe[k - 1] * d(-k);
            const xs = [];
            for (let b = mLo + 2; b <= mHi - 1; b++) {
                const x = lvl(d, b) - lvl(d, b - 1);
                if (x === 0) continue;
                const g = x > 0 ? gUp : gDn;
                if (b >= bLo && b <= bHi) xs.push([b, x, g]); else base += x * at(g, s - b * spb);
            }
            let dist = new Map([[base, 1]]);
            for (const [b, x, g] of xs) {
                const nd = new Map();
                for (const [v, q] of dist) for (const t of nodes) {
                    const k = Math.round((v + x * at(g, s - b * spb - t / dt)) * 1e6) / 1e6;
                    nd.set(k, (nd.get(k) || 0) + q / 4);
                }
                dist = nd;
            }
            for (const [v, q] of dist) exact[d(0)].push([v, q]);
        }
        const eng = se.evaluate(phase, grid);
        for (const [cond, pm] of [[1, eng.p1], [-1, eng.p0]]) {
            const ex = exact[cond];
            let tot = 0, mE = 0, mS = 0, minE = Infinity, maxE = -Infinity, minS = Infinity, maxS = -Infinity, band = 0;
            for (const [, q] of ex) tot += q;
            for (const [v, q] of ex) { mE += v * q / tot; minE = Math.min(minE, v); maxE = Math.max(maxE, v); }
            for (let i = 0; i < nV; i++) { const v = grid.v0 + i * bin; mS += v * pm[i]; if (pm[i] > 0) { minS = Math.min(minS, v); maxS = Math.max(maxS, v); } }
            // The exact CDF lies between the engine's CDF shifted by +-8 bins at every level.
            for (let lv = -0.6; lv <= 0.9; lv += 0.01) {
                let cE = 0, lo = 0, hi = 0;
                for (const [v, q] of ex) if (v < lv) cE += q / tot;
                for (let i = 0; i < nV; i++) { const v = grid.v0 + i * bin; if (v < lv - 8 * bin) lo += pm[i]; if (v < lv + 8 * bin) hi += pm[i]; }
                if (cE < lo - 1e-12 || cE > hi + 1e-12) band++;
            }
            worst = Math.max(worst, Math.abs(mE - mS) / bin, Math.abs(minE - minS) / bin, Math.abs(maxE - maxS) / bin, band ? 99 : 0);
            rows.push(`${phase}/${cond}: mean ${(mE * 1e3).toFixed(2)}/${(mS * 1e3).toFixed(2)} mV, support ${(minE * 1e3).toFixed(1)}..${(maxE * 1e3).toFixed(1)} vs ${(minS * 1e3).toFixed(1)}..${(maxS * 1e3).toFixed(1)}`);
        }
    }
    check('chain = brute-force enumeration (FIR, rise != fall, DFE, dual Dirac): means and supports within 2 bins, CDF within 8 bins',
        worst < 2 && !se.symmetric, `worst ${worst.toFixed(2)} bins; ${rows.join('; ')}`);
}

// 2. Lossless matched line, analytic eye width.
{
    const lossless = seTable(() => ({ R: 0, L: L0, G: 0, C: C0 }));
    const ui = 100e-12, j = { rj: 1e-12, djdd: 10e-12, pwRj: 0.6e-12, pwDj: 4e-12 };
    const cfg = { length: 0.1, samplesPerEdge: 16, network: { rs: 50, tx: { ac: false }, rx: { ac: false }, term: { type: 'r', rt: 50 } },
        stimulus: { type: 'pattern', bitRate: 1 / ui, bits: patternBits({ pattern: 'prbs7' }), vLow: 0, vHigh: 1, tr: 10e-12, tf: 10e-12,
            edgeDef: '10-90', shape: 'gauss', jitter: j }, statEye: { ber: 1e-12 } };
    const r = simulateLink(lossless, cfg, fft);
    const sig = Math.sqrt(j.rj ** 2 + j.pwRj ** 2 / 2), shifts = [];
    for (const a of [-j.djdd / 2, j.djdd / 2]) for (const b of [-j.pwDj / 4, j.pwDj / 4]) shifts.push(a + b);
    const Q = (z) => 0.5 * erfc(z / Math.SQRT2);
    const pLate = (t) => shifts.reduce((s, a) => s + Q((t - a) / sig), 0) / shifts.length;
    // BER(t) = P(delta > t) / 2 + P(delta < t - UI) / 2 = target, t from the crossing.
    let lo = 0, hi = ui / 2;
    for (let i = 0; i < 100; i++) { const m = 0.5 * (lo + hi); if (0.5 * pLate(m) + 0.5 * pLate(ui - m) > 1e-12) lo = m; else hi = m; }
    const want = (ui - 2 * lo) / ui;
    const s = r.statEye;
    check('lossless line: statistical eye width = analytic BER 1e-12 width (Rj, PW Rj, DjDD, PW DDj)', Math.abs(s.width - want) < 1.5e-3,
        `${s.width.toFixed(4)} vs ${want.toFixed(4)} UI`);
    // Symmetric jitter: the middle of the eye width is the mean UI centre, the sample point
    // sits `offset` from it, the height is the amplitude.
    const mid = 0.5 * (s.widthLeft + s.widthRight);
    check('lossless line: eye height = amplitude (0.5 V) within 2 bins, the eye width centred on the mean UI centre',
        Math.abs(s.height - 0.5) < 2 * s.bin && Math.abs(mid + s.offset) < 0.01 && Math.abs(s.offset) <= 0.1 + 1e-9,
        `${s.height.toFixed(4)} V, bin ${(s.bin * 1e3).toFixed(2)} mV; sample point ${s.offset} UI from the mean centre, width middle ${mid.toFixed(4)} UI from the sample point`);
    check('lossless line: symmetric eye (one chain run)', s.symmetric === true);
}

// 3. Short memory, no jitter: the worst case is above 1e-12.
{
    // A lossy line with a resistive source mismatch: a few UI of ISI, one small reflection.
    const line = seTable(() => ({ R: 60, L: L0, G: 0, C: C0 }));
    const base = (rxEq) => ({ length: 0.1, samplesPerEdge: 16, network: { rs: 40, tx: { ac: false }, rx: { ac: false }, term: { type: 'r', rt: 50 } },
        stimulus: { type: 'pattern', bitRate: 4e9, bits: patternBits({ pattern: 'prbs15' }), vLow: 0, vHigh: 1, tr: 60e-12, tf: 60e-12,
            edgeDef: '10-90', shape: 'gauss' }, rxEq, statEye: { ber: 1e-12, nV: 4097 } });
    for (const taps of [0, 2]) {
        // A 0 dB CTLE with far poles makes the equalization path run without changing the signal.
        const r = simulateLink(line, base({ dfeTaps: taps, ctle: taps ? null : { adcDb: 0, fp1: 1e13, fp2: 2e13 } }), fft);
        const s = r.statEye, pda = r.eqChoice.pdaHeight, td = volts(metric(r, 'Eye height (centre)'));
        const bin = s.bin;
        check(`short memory, ${taps} DFE taps: statistical eye height = peak-distortion worst case = PRBS15 eye`,
            Math.abs(s.height - pda) < Math.max(1e-4, 2 * bin) && Math.abs(s.height - td) < Math.max(1e-4, 2 * bin),
            `stat ${(s.height * 1e3).toFixed(3)} mV at ${s.offset} UI, peak distortion ${(pda * 1e3).toFixed(3)} mV, time domain ${(td * 1e3).toFixed(3)} mV, bin ${(bin * 1e3).toFixed(3)} mV`);
    }
}

// 4. Baseline wander.
{
    const line = seTable(() => ({ R: 5, L: L0, G: 0, C: C0 }));
    const ui = 1e-9, c = 100e-9;
    const cfg = { length: 0.1, samplesPerEdge: 16, network: { rs: 50, tx: { ac: true, c }, rx: { ac: false }, term: { type: 'r', rt: 50 } },
        stimulus: { type: 'pattern', bitRate: 1 / ui, bits: patternBits({ pattern: 'prbs7' }), vLow: 0, vHigh: 1, tr: 100e-12, tf: 100e-12,
            edgeDef: '10-90', shape: 'gauss' }, statEye: { ber: 1e-12 } };
    const r = simulateLink(line, cfg, fft);
    const got = volts(metric(r, 'Baseline wander (AC coupling)') || 'NaN');
    // One pole, lambda = C (Rs + R l + Rt); the slow mode's amplitude is the passband level.
    const lam = c * (50 + 0.5 + 50), a = Math.exp(-ui / lam), amp = 0.5 * 50 / 100.5;
    const want = amp * Math.sqrt((1 - a) / (1 + a));
    check('AC coupling: baseline wander rms = passband level x low-pass of random bits', Math.abs(got / want - 1) < 0.01,
        `${(got * 1e3).toFixed(4)} mV vs ${(want * 1e3).toFixed(4)} mV (lambda ${(lam * 1e6).toFixed(2)} us)`);
}

// 5. Equalization search.
{
    const keys = ['R11', 'R12', 'R22', 'L11', 'L12', 'L22', 'G11', 'G12', 'G22', 'C11', 'C12', 'C22'];
    const data = Object.fromEntries(keys.map(k => [k, []]));
    const fn = (f) => {
        const ff = Math.max(f, 1e6), r = 5 + 1.5e-3 * Math.sqrt(ff), li = 1.5e-3 * Math.sqrt(ff) / (2 * Math.PI * ff);
        return { R11: r, R12: 0.1 * r, R22: r, L11: 330e-9 + li, L12: 50e-9, L22: 330e-9 + li, G11: 0, G12: 0, G22: 0, C11: 110e-12, C12: -10e-12, C22: 110e-12 };
    };
    for (const f of FREQS) { const v = fn(f); for (const k of keys) data[k].push(v[k]); }
    const top = fn(FREQS[FREQS.length - 1]);
    const pair = { differential: true, freqs: Float64Array.from(FREQS), keys,
        data: Object.fromEntries(Object.entries(data).map(([k, v]) => [k, Float64Array.from(v)])), lext: { L11: top.L11, L12: top.L12, L22: top.L22 } };
    const cfg = (statEye) => ({ length: 1.0, samplesPerEdge: 16,
        network: { rs: 50, tx: { ac: false }, rx: { ac: false }, term: { type: 'r', rt: 50 }, padTx: 0.5e-12, padRx: 0.5e-12,
            pkgTx: PACKAGES.pcie16Root.model, pkgRx: PACKAGES.pcie16NonRoot.model },
        stimulus: { type: 'pattern', bitRate: 16e9, bits: patternBits({ pattern: 'prbs7' }), vLow: 0, vHigh: 0.804, tr: 21.875e-12, tf: 21.875e-12,
            edgeDef: '10-90', shape: 'gauss', jitter: { rj: 0.71e-12, djdd: 3.75e-12, pwRj: 0.54e-12, pwDj: 5e-12 } },
        txEq: { auto: true, presets: ['P4', 'P5'] }, rxEq: { ctle: { auto: true, adcMin: -10, adcMax: -8, fp1: 2e9, fp2: 16e9 }, dfeTaps: 2, dfeLimits: [0.03, 0.02] },
        statEye });
    const name = (c) => `${c.firName}/${c.adc}`;
    const offsets = [-0.05, 0, 0.05];
    const all = simulateLink(pair, cfg({ ber: 1e-12, offsets, sweep: 6, refine: 18 }), fft);
    const def = simulateLink(pair, cfg({ ber: 1e-12, offsets }), fft);
    const pda = simulateLink(pair, cfg(null), fft);
    const prod = (e) => e.height * e.width;
    const best = all.eqSearch.full.reduce((a, b) => prod(b) > prod(a) ? b : a);
    const pairName = (e) => `${name(e)}@${e.offset ?? e.sampleOffset}`;
    check('equalization search: EH x EW argmax over every candidate and sample point, the default search agrees, the reported eye is the winner\'s',
        all.eqSearch.full.length === 18 && pairName(all.eqChoice) === pairName(best) && pairName(def.eqChoice) === pairName(best)
        && Math.abs(def.statEye.height - best.height) < 1e-12 && Math.abs(def.statEye.width - best.width) < 1e-12,
        `${pairName(best)}: ${(best.height * 1e3).toFixed(1)} mV x ${best.width.toFixed(3)} UI; full grid ${def.eqSearch.full.length} of 18; `
        + all.eqSearch.full.slice().sort((a, b) => prod(b) - prod(a)).slice(0, 6).map(e => `${pairName(e)} ${(prod(e) * 1e3).toFixed(2)}`).join(', '));
    check('without the statistical eye: the peak-distortion choice (different here)',
        name(pda.eqChoice) !== name(best) && !pda.eqSearch && pda.eqChoice.sampleOffset === null && !pda.statEye,
        `${name(pda.eqChoice)} (${(pda.eqChoice.pdaHeight * 1e3).toFixed(1)} mV) vs EH x EW ${name(best)}`);
}

done();
