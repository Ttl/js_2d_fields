// Effect of the low-frequency line-model notes on a time-domain result.
//
// The solvers flag two approximations that only hold above some frequency:
//   narrow-gap: the quasi-static surface-impedance loss overstates the skin-effect
//     loss of conductors facing each other across a gap g by about 25 delta/g
//     (field_solver.js _narrow_gap_note), noted where delta/g > 0.05.
//   skin-transition: below the frequency where the skin depth is half the conductor
//     thickness, the quasi-static R and internal L blend the DC and skin solutions
//     with a calibrated correction (field_solver.js, up to 7 %); taken as a 10 % error
//     of the skin-effect part of R and L.
//   ground-spreading: below f_on the return current spreads sideways in the ground
//     planes, which are modelled as ideal returns; the loop inductance rises by about
//     mu0/(2 pi) ln(l/g) as the return widens to l = g f_on / f (until the board
//     width), the ground resistance is approximate.
// Each note becomes a perturbation of the line parameters by its estimated error; the
// bound on the resulting change of the receiver waveform is the measure of whether the
// note matters for this stimulus and circuit.

import { mainName } from './transient.js';
import { makeRLGCEvaluator } from './channel_grid.js';
import { buildStimulus, legDerivSpectrum } from './stimulus.js';
import { buildNetwork, makeTransfer } from './network.js';

const MU0 = 4e-7 * Math.PI;
// Widest return assumed for the ground spreading.
const BOARD_WIDTH = 0.1;
// Largest relative loss error of the narrow-gap estimate (delta/g = 1).
const GAP_ERR_MAX = 0.25;
// Error of the skin-effect part of R and L in the DC-skin transition.
const SKIN_ERR = 0.1;

// The notes of the channel, merged: the narrowest gap and the highest onset.
function lowFrequencyNotes(warnings) {
    let gap = null, spread = null, skin = null;
    for (const w of warnings || []) {
        if (!w || typeof w !== 'object') continue;
        if (w.reason === 'narrow-gap' && w.gap > 0 && w.sigma > 0) {
            if (!gap || w.gap < gap.gap) gap = { reason: 'narrow-gap', gap: w.gap, sigma: w.sigma };
        } else if (w.reason === 'skin-transition' && w.fBand > 0) {
            if (!skin || w.fBand > skin.fBand) skin = { reason: 'skin-transition', fBand: w.fBand };
        } else if (w.reason === 'ground-spreading' && w.fOn > 0 && w.gap > 0) {
            if (!spread || w.fOn > spread.fOn) spread = { reason: 'ground-spreading', fOn: w.fOn, gap: w.gap };
        }
    }
    return [gap, skin, spread].filter(Boolean);
}

// Frequency below which a note applies.
function noteBand(note) {
    if (note.reason === 'ground-spreading') return note.fOn;
    if (note.reason === 'skin-transition') return note.fBand;
    // delta = 0.05 g
    const d = 0.05 * note.gap;
    return 1 / (Math.PI * MU0 * note.sigma * d * d);
}

/**
 * Perturbation of the evaluated line parameters v at f for a note.
 * dc: the parameters at DC, lext: the external inductance entries (table.lext).
 */
function perturbation(note, dc, lext) {
    if (note.reason === 'narrow-gap') {
        return (v, f) => {
            if (!(f > 0)) return v;
            const r = Math.sqrt(1 / (Math.PI * f * MU0 * note.sigma)) / note.gap;
            if (r <= 0.05) return v;
            const k = 1 - GAP_ERR_MAX * Math.min(r, 1);
            const out = { ...v };
            // The skin-effect part of R and L: R above its DC value, L above the
            // external inductance.
            for (const key of Object.keys(v)) {
                if (key[0] === 'R') out[key] = dc[key] + (v[key] - dc[key]) * k;
                if (key[0] === 'L' && Number.isFinite(lext[key])) out[key] = lext[key] + (v[key] - lext[key]) * k;
            }
            return out;
        };
    }
    if (note.reason === 'skin-transition') {
        return (v, f) => {
            if (!(f > 0) || f >= note.fBand) return v;
            const out = { ...v };
            for (const key of Object.keys(v)) {
                if (key[0] === 'R') out[key] = dc[key] + (v[key] - dc[key]) * (1 + SKIN_ERR);
                if (key[0] === 'L' && Number.isFinite(lext[key])) out[key] = lext[key] + (v[key] - lext[key]) * (1 + SKIN_ERR);
            }
            return out;
        };
    }
    const lnMax = Math.log(Math.max(1, BOARD_WIDTH / note.gap));
    return (v, f) => {
        const ln = f > 0 ? Math.min(lnMax, Math.log(Math.max(1, note.fOn / f))) : lnMax;
        if (!(ln > 0)) return v;
        // A common return: the same inductance on every entry of a pair.
        const dL = MU0 / (2 * Math.PI) * ln;
        const out = { ...v };
        for (const key of Object.keys(v)) if (key[0] === 'L') out[key] = v[key] + dL;
        return out;
    };
}

// Upper bound of the change of the main receiver voltage, over all time, when the line
// parameters carry each note's estimated error: (1/pi) int_0^inf |dO(w)| dw with
// dO = sum over legs of dH_l(w) U_l(w), U the driver spectrum. A steady clock sums
// 2 |dO_n| over its harmonics instead. Unlike a second time-domain run this sees
// changes below the window's frequency step: their slow tails are bounded too.
// Returns [{ reason, band, impact }], impact relative to the receiver swing of base.
function lowFrequencyImpact(table, cfg, base, notes) {
    const diff = table.differential;
    const stim = buildStimulus({ ...cfg.stimulus, differential: diff });
    const net = buildNetwork({ ...cfg.network, differential: diff });
    const ev0 = makeRLGCEvaluator(table), dc = ev0(0), lext = table.lext || {};
    const { d } = net.nodes, nl = net.nl;
    const transfer = makeTransfer(net, cfg.length, diff ? [d[0], d[1]] : [d[0]]);
    const scratch = () => ({ re: new Float64Array(2 * nl), im: new Float64Array(2 * nl) });
    const s0 = scratch(), s1 = scratch();
    const main = base.traces[mainName(diff)];
    let lo = Infinity, hi = -Infinity;
    for (const v of main) { lo = Math.min(lo, v); hi = Math.max(hi, v); }
    const swing = hi - lo;
    // |dO| at f for driver spectra U (one [re, im] per leg).
    const dOut = (f, pert, U) => {
        const a = transfer(f, ev0(f), s0), b = transfer(f, pert(ev0(f), f), s1);
        let re = 0, im = 0;
        for (let l = 0; l < nl; l++) {
            // Main output: rx (single-ended) or rxP - rxN.
            const dr = (b.re[l] - a.re[l]) - (diff ? b.re[nl + l] - a.re[nl + l] : 0);
            const di = (b.im[l] - a.im[l]) - (diff ? b.im[nl + l] - a.im[nl + l] : 0);
            re += dr * U[l][0] - di * U[l][1];
            im += dr * U[l][1] + di * U[l][0];
        }
        return Math.hypot(re, im);
    };
    return notes.map(note => {
        const band = noteBand(note), pert = perturbation(note, dc, lext);
        let bound = 0;
        if (stim.steady) {
            const P = stim.period, nMax = Math.min(100000, Math.floor(band * P));
            for (let n = 1; n <= nMax; n++) {
                const f = n / P, w = 2 * Math.PI * f;
                // Fourier coefficients of the legs: U'_period(f) / (j w P).
                const U = stim.legs.map(leg => { const [r, i] = legDerivSpectrum(stim, leg, f, true); return [i / (w * P), -r / (w * P)]; });
                bound += 2 * dOut(f, pert, U);
            }
        } else {
            // Log grid up to the band edge; the integrand vanishes towards DC (dH -> 0).
            const perDecade = 200, f0 = band * 1e-7;
            const n = Math.ceil(7 * perDecade), step = Math.log(band / f0) / n;
            let prev = null;
            for (let k = 0; k <= n; k++) {
                const f = f0 * Math.exp(k * step), w = 2 * Math.PI * f;
                const U = stim.legs.map(leg => { const [r, i] = legDerivSpectrum(stim, leg, f); return [i / w, -r / w]; });
                // dw = w d(ln f)
                const g = dOut(f, pert, U) * w;
                if (prev !== null) bound += 0.5 * (prev + g) * step;
                prev = g;
            }
            bound /= Math.PI;
        }
        return { reason: note.reason, band, impact: swing > 0 ? bound / swing : 0 };
    });
}

export { lowFrequencyImpact, lowFrequencyNotes, noteBand, perturbation, BOARD_WIDTH, GAP_ERR_MAX, SKIN_ERR };
