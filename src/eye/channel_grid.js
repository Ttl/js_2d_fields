// Per-unit-length line parameters for the time-domain simulation.
//
// A channel table holds the RLGC of the solved line at a set of sweep frequencies,
// starting with DC. For a pair the entries are the physical 2x2 matrices (the ones the
// MTL S-parameters use), stored as their three independent entries. The evaluator
// interpolates the table at any frequency:
//   0 < f < f1 (first nonzero row): blended from DC in sqrt(f) (G linearly in f),
//   f1..f_top: cubic splines in log f (R and G in log value when positive),
//   f > f_top: R ~ sqrt(f), internal L ~ 1/sqrt(f), G ~ f, C constant.
import { CubicSpline } from '../interpolating_sweep.js';
import { buildPhysicalRLGC } from '../sparameters.js';

const PAIR_KEYS = ['R11', 'R12', 'R22', 'L11', 'L12', 'L22', 'G11', 'G12', 'G22', 'C11', 'C12', 'C22'];
const SE_KEYS = ['R', 'L', 'G', 'C'];

const sym3 = (m) => [m[0][0], 0.5 * (m[0][1] + m[1][0]), m[1][1]];
const symOE = (o, e) => [(o + e) / 2, (e - o) / 2, (o + e) / 2];

/**
 * Channel table from sweep rows [{ freq, result }]. Rows must include f = 0.
 * Returns { differential, freqs, keys, data: { key: Float64Array }, lext }, lext the
 * external inductance entries (the high-frequency limit of L) at the top row.
 */
function channelTableFromRows(rows) {
    if (!rows || !rows.length) throw new Error('No sweep results for the channel.');
    const sorted = [...rows].filter(r => r && r.result && r.result.modes && r.result.modes.length)
        .sort((a, b) => a.freq - b.freq);
    // One row per frequency.
    const uniq = sorted.filter((r, i) => i === 0 || r.freq !== sorted[i - 1].freq);
    if (!uniq.length || uniq[0].freq !== 0) throw new Error('The channel needs a DC (0 Hz) sweep point.');
    const m0 = uniq[0].result.modes[0];
    if (m0.self_referenced) throw new Error('Time-domain simulation needs a TEM or quasi-TEM line.');
    const differential = uniq[0].result.modes.length === 2;
    const keys = differential ? PAIR_KEYS : SE_KEYS;
    const data = Object.fromEntries(keys.map(k => [k, new Float64Array(uniq.length)]));
    let lext = null;
    uniq.forEach((row, i) => {
        const res = row.result;
        if (!differential) {
            const m = res.modes[0], { R, L, G, C } = m.RLGC;
            data.R[i] = R; data.L[i] = L; data.G[i] = G; data.C[i] = C;
            if (i === uniq.length - 1) lext = { L: Number.isFinite(m.L_external) ? m.L_external : L };
            return;
        }
        const odd = res.modes.find(m => m.mode === 'odd') || res.modes[0];
        const even = res.modes.find(m => m.mode === 'even') || res.modes[1];
        const ph = buildPhysicalRLGC(odd.RLGC, even.RLGC, res.physMatrix || null);
        const put = (name, m) => { const v = sym3(m); data[name + '11'][i] = v[0]; data[name + '12'][i] = v[1]; data[name + '22'][i] = v[2]; };
        put('R', ph.R); put('L', ph.L); put('G', ph.G); put('C', ph.C);
        if (i === uniq.length - 1) {
            let le;
            if (res.physMatrix && res.physMatrix.L) le = sym3(res.physMatrix.L);
            else {
                const lo = Number.isFinite(odd.L_external) ? odd.L_external : odd.RLGC.L;
                const lev = Number.isFinite(even.L_external) ? even.L_external : even.RLGC.L;
                le = symOE(lo, lev);
            }
            lext = { L11: le[0], L12: le[1], L22: le[2] };
        }
    });
    return { differential, freqs: Float64Array.from(uniq.map(r => r.freq)), keys, data, lext };
}

const kind = (key) => key[0];   // R, L, G or C

/**
 * Evaluator of a channel table: f (Hz) -> { key: value }.
 */
function makeRLGCEvaluator(table) {
    const { freqs, keys, data } = table;
    const n = freqs.length;
    const fTop = freqs[n - 1];
    const f1 = n > 1 ? freqs[1] : 0;
    const dc = Object.fromEntries(keys.map(k => [k, data[k][0]]));
    const ts = Array.from(freqs.subarray(1), f => Math.log10(f));
    const splines = {};
    for (const k of keys) {
        const ys = Array.from(data[k].subarray(1));
        if (ys.length === 0) { splines[k] = null; continue; }
        if (ys.length === 1) { const v = ys[0]; splines[k] = { evaluate: () => v }; continue; }
        const useLog = (kind(k) === 'R' || kind(k) === 'G') && ys.every(v => v > 0);
        if (useLog) {
            const s = new CubicSpline(ts, ys.map(Math.log));
            splines[k] = { evaluate: (t) => Math.exp(s.evaluate(t)) };
        } else {
            splines[k] = new CubicSpline(ts, ys);
        }
    }
    const top = Object.fromEntries(keys.map(k => [k, data[k][n - 1]]));
    const first = Object.fromEntries(keys.map(k => [k, n > 1 ? data[k][1] : data[k][0]]));
    const lext = table.lext || {};

    return function evaluate(f) {
        const out = {};
        if (f <= 0 || n === 1) {
            for (const k of keys) out[k] = dc[k];
            return out;
        }
        if (f < f1) {
            const s = Math.sqrt(f / f1), r = f / f1;
            for (const k of keys) out[k] = dc[k] + (first[k] - dc[k]) * (kind(k) === 'G' ? r : s);
            return out;
        }
        if (f <= fTop) {
            const t = Math.log10(f);
            for (const k of keys) out[k] = splines[k].evaluate(t);
            return out;
        }
        const s = Math.sqrt(f / fTop);
        for (const k of keys) {
            const v = top[k];
            switch (kind(k)) {
                case 'R': out[k] = v * s; break;
                case 'G': out[k] = v * (f / fTop); break;
                case 'L': {
                    const le = Number.isFinite(lext[k]) ? lext[k] : v;
                    out[k] = le + (v - le) / s;
                    break;
                }
                default: out[k] = v;
            }
        }
        return out;
    };
}

// Lowest sweep frequency above DC a channel needs. Below the first row the evaluator
// blends from DC in sqrt(f), which misses the fall of the internal inductance through the
// DC-skin transition: an L without the matching R, a receiver response before arrival.
function channelLowFreq(fNeed) {
    return Math.min(10e6, fNeed / 100);
}

/**
 * Whether sweep rows cover a time-domain simulation up to fNeed.
 * interp: the rows come from an interpolating sweep (evaluated on a dense grid).
 * Returns { ok, reason, warning }.
 */
function checkCoverage(rows, fNeed, interp) {
    if (!rows || !rows.length) return { ok: false, reason: 'no sweep results' };
    const fs = [...new Set(rows.map(r => r.freq))].sort((a, b) => a - b);
    if (fs[0] !== 0) return { ok: false, reason: 'the sweep has no DC (0 Hz) point' };
    const top = fs[fs.length - 1];
    if (top < fNeed) {
        return { ok: false, reason: `the sweep stops at ${fmtHz(top)}, below the ${fmtHz(fNeed)} the stimulus needs` };
    }
    if (interp) return { ok: true };
    const nz = fs.filter(f => f > 0);
    if (nz.length < 8) return { ok: false, reason: 'the sweep has too few points' };
    if (nz[0] > 2 * channelLowFreq(fNeed)) {
        return { ok: false, reason: `the first sweep point above DC (${fmtHz(nz[0])}) is too high to resolve the internal inductance at low frequencies` };
    }
    // Large gaps in log f relative to the decade the line parameters change over.
    let worst = 1;
    for (let i = 1; i < nz.length; i++) worst = Math.max(worst, nz[i] / nz[i - 1]);
    return { ok: true, warning: worst > 4 ? 'the sweep points are sparse, the line parameters are splined over large gaps' : null };
}

function fmtHz(f) {
    if (f >= 1e9) return `${+(f / 1e9).toPrecision(4)} GHz`;
    if (f >= 1e6) return `${+(f / 1e6).toPrecision(4)} MHz`;
    if (f >= 1e3) return `${+(f / 1e3).toPrecision(4)} kHz`;
    return `${+f.toPrecision(4)} Hz`;
}

export { channelTableFromRows, makeRLGCEvaluator, checkCoverage, channelLowFreq, fmtHz, PAIR_KEYS, SE_KEYS };
