// Statistical eye of an NRZ link: the distribution of the slicer signal at each sampling
// phase over all bit sequences (independent, equiprobable bits) and the transmitter jitter
// of every bit boundary, and the BER contours, eye height and eye width it gives.
//
// With the transition responses G+ and G- (the slicer signal per unit FIR transition up
// or down, half the swing included, each with its edge shape), FIR levels
// L_b = sum_j c_j d_(b-j) and transitions x_b = L_b - L_(b-1), the signal at sample time s
// (bit 0 from its boundary at 0) is
//   V(s) = Vmid + sum_b x_b G_sgn(x_b)(s - b UI - delta_b) - sum_k h_k d_-k
// with delta_b the jitter of boundary b and h the DFE taps (ideal decisions). It is split
// into a part linear in the bits, sum_m d_m w_m with w_m the cursors of the jitter-free
// pulse response of G_lin = (G+ + G-) / 2 less the DFE taps, and per-boundary residuals
//   x_b (G_sgn(s - b UI - delta_b) - G_lin(s - b UI)),
// which only differ from zero near an arriving edge (jitter, unequal rise and fall). The
// boundaries with a residual form a window that runs as a Markov chain over the bits each
// transition depends on; the bits outside it are independent +-w_m terms, convolved
// first. Residuals below half a voltage bin (the slow tail of a lossy line moved by the
// jitter) enter as independent Gaussian noise of their mean and variance instead. Slow
// AC-coupling modes add a Gaussian baseline wander.
import { erfc } from './stimulus.js';

const Phi = (z) => 0.5 * erfc(-z / Math.SQRT2);

// Gaussian jitter nodes per sigma, and their extent (one-sided tail beyond 8 sigma 6e-16).
const NODES_PER_SIGMA = 8, SIGMA_SPAN = 8;
// Mass below which a probability bin at the end of a support is dropped.
const PRUNE = 1e-40;

/**
 * Jitter displacement distribution of one bit boundary as nodes and probabilities: a
 * Gaussian (Rj with the PW Rj) around the four dual-Dirac positions (DjDD, PW DDj), the
 * split jitter.js uses. Returns { t: Float64Array, p: Float64Array, ext } (ext: largest
 * |t|).
 */
function jitterNodes(j) {
    const sig = j ? Math.sqrt((j.rj || 0) ** 2 + 0.5 * (j.pwRj || 0) ** 2) : 0;
    const a = j ? 0.5 * (j.djdd || 0) : 0, b = j ? 0.25 * (j.pwDj || 0) : 0;
    let gt = [0], gp = [1];
    if (sig > 0) {
        const n = NODES_PER_SIGMA * SIGMA_SPAN, h = 1 / NODES_PER_SIGMA;
        gt = []; gp = [];
        for (let i = -n; i <= n; i++) {
            // Bin [i - 1/2, i + 1/2] h of the unit Gaussian, the end nodes take the tails.
            const lo = i === -n ? 0 : Phi((i - 0.5) * h), hi = i === n ? 1 : Phi((i + 0.5) * h);
            gt.push(i * h * sig); gp.push(hi - lo);
        }
    }
    const dd = [];
    for (const sa of a > 0 ? [-a, a] : [0]) for (const sb of b > 0 ? [-b, b] : [0]) dd.push(sa + sb);
    const t = new Float64Array(dd.length * gt.length), p = new Float64Array(t.length);
    let k = 0, ext = 0;
    for (const d of dd) for (let i = 0; i < gt.length; i++) {
        t[k] = d + gt[i]; p[k] = gp[i] / dd.length; ext = Math.max(ext, Math.abs(t[k])); k++;
    }
    return { t, p, ext };
}

// Response at a fractional sample index: 0 before the window, settled after it.
function sampler(a) {
    const N = a.length;
    return (x) => {
        if (x <= 0) return x === 0 ? a[0] : 0;
        if (x >= N - 1) return a[N - 1];
        const i = Math.floor(x), f = x - i;
        return a[i] + (a[i + 1] - a[i]) * f;
    };
}

/**
 * Statistical eye engine.
 * p: { gUp, gDn (Float64Array transition responses per unit FIR step, half the swing
 *      included, sampled at dt, the boundary at index 0), fir [c-1, c0, c+1], spb (samples
 *      per bit), dt, n0 (sample index of bit 0's cursor), dfe (taps), jitter ({ rj, djdd,
 *      pwRj, pwDj } or null), vmid (threshold level, the mean signal), sigmaSlow (baseline
 *      wander rms) }
 * Returns { layout(s, eps, weak), evaluate(phase, grid), symmetric }: s a sample time,
 * phase in samples from n0.
 */
function makeStatEye(p) {
    const { fir, spb, dt, n0 } = p;
    const N = p.gUp.length;
    const gUp = sampler(p.gUp), gDn = sampler(p.gDn);
    const lin = new Float64Array(N);
    for (let n = 0; n < N; n++) lin[n] = 0.5 * (p.gUp[n] + p.gDn[n]);
    const gLin = sampler(lin);
    let gMax = 0, dirDiff = 0;
    for (let n = 0; n < N; n++) { gMax = Math.max(gMax, Math.abs(lin[n])); dirDiff = Math.max(dirDiff, Math.abs(p.gUp[n] - p.gDn[n])); }
    // Odd symmetry in the bits: the zero distribution mirrors the one distribution.
    const symmetric = dirDiff <= 1e-12 * gMax;
    const jn = jitterNodes(p.jitter);
    const jExt = jn.ext / dt;
    const jt = Float64Array.from(jn.t, t => t / dt);
    const pre = fir[0] !== 0 ? 1 : 0, post = fir[2] !== 0 ? 1 : 0;
    const xMax = 2 * (Math.abs(fir[0]) + Math.abs(fir[1]) + Math.abs(fir[2]));
    const dfe = p.dfe || [];

    // Linear weights w_m at sample time s: the pulse response cursors and the DFE.
    const weights = (s) => {
        const mLo = Math.floor((s - N) / spb) - 3, mHi = Math.ceil(s / spb) + 3;
        const w = new Map();
        for (let m = mLo; m <= mHi; m++) {
            let v = 0;
            for (let j = -1; j <= 1; j++) {
                const c = fir[j + 1];
                if (c !== 0) v += c * (gLin(s - (m + j) * spb) - gLin(s - (m + j + 1) * spb));
            }
            if (m < 0 && -m <= dfe.length) v -= dfe[-m - 1];
            if (v !== 0) w.set(m, v);
        }
        return w;
    };

    // Largest residual of boundary b, and the largest change the jitter makes to it.
    const residual = (s, b) => {
        const tau = s - b * spb;
        if (tau + jExt <= 0 || tau - jExt >= N - 1) return { max: 0, jitter: 0 };
        const gl = gLin(tau);
        let mx = 0, spread = 0;
        const lo = Math.max(0, Math.ceil(tau - jExt)), hi = Math.min(N - 1, Math.floor(tau + jExt));
        const at = (x) => {
            const u = gUp(x), d = gDn(x);
            mx = Math.max(mx, Math.abs(u - gl), Math.abs(d - gl));
            spread = Math.max(spread, Math.abs(u - gUp(tau)), Math.abs(d - gDn(tau)));
        };
        at(tau - jExt); at(tau + jExt); at(tau);
        for (let n = lo; n <= hi; n++) at(n);
        return { max: xMax * mx, jitter: xMax * spread };
    };

    // Transition values over the bits it depends on (equiprobable): x and its probability.
    const xDist = new Map();
    for (let code = 0; code < 16; code++) {
        const d = (k) => ((code >> (k + 1)) & 1) ? 1 : -1;   // d(-1) = d_(b+1) .. d(2) = d_(b-2)
        const x = fir[0] * (d(-1) - d(0)) + fir[1] * (d(0) - d(1)) + fir[2] * (d(1) - d(2));
        xDist.set(x, (xDist.get(x) || 0) + 1 / 16);
    }
    // Mean and variance of a boundary's residual over the bits and the jitter.
    const residualMoments = (s, b) => {
        const tau = s - b * spb, gl = gLin(tau);
        let m = 0, m2 = 0;
        for (const [x, px] of xDist) {
            if (x === 0) continue;
            const g = x > 0 ? gUp : gDn;
            for (let i = 0; i < jt.length; i++) {
                const v = x * (g(tau - jt[i]) - gl), q = px * jn.p[i];
                m += q * v; m2 += q * v * v;
            }
        }
        return { mean: m, var: Math.max(0, m2 - m * m) };
    };

    // Bound of |V - vmid| over all bits and jitter at s, and the residual window. Residuals
    // above `eps` join the chain; those up to `weak` only as noise moments.
    const layout = (s, eps, weak = 0) => {
        const w = weights(s);
        let bound = 0;
        for (const v of w.values()) bound += Math.abs(v);
        const res = new Map();
        const bLo = Math.floor((s - N - jExt) / spb) - 1, bHi = Math.ceil((s + jExt) / spb) + 1;
        let bA = 0, bB = 0, noiseMean = 0, noiseVar = 0;
        for (let b = bLo; b <= bHi; b++) {
            const r = residual(s, b);
            if (r.max <= eps) continue;
            bound += r.max;
            if (r.max <= weak) {
                const mo = residualMoments(s, b);
                noiseMean += mo.mean; noiseVar += mo.var;
                continue;
            }
            res.set(b, r);
            bA = Math.min(bA, b); bB = Math.max(bB, b);
        }
        return { w, res, bound, bA, bB, noiseMean, noiseVar };
    };

    /**
     * Conditional distributions of V at sample n0 + phase on the grid
     * { v0, bin, nV } (bin centres v0 + i bin). Returns { p1, p0 } (Float64Array, mass 1).
     */
    const evaluate = (phase, grid) => {
        const s = n0 + phase, { bin, nV } = grid;
        const c0 = Math.round((p.vmid - grid.v0) / bin);
        const lay = layout(s, 0.02 * bin, 0.5 * bin);
        const { w, res, bA, bB } = lay;
        // Bits of the chain lo .. hi; a transition depends on K + 1 of them.
        const lo = bA - 1 - post, hi = bB + pre;
        const K = 1 + post + pre;
        // Bits outside the chain: independent +-w terms.
        let base = new Float64Array(nV);
        base[c0] = 1;
        let bLo = c0, bHi = c0;
        const small = [];
        for (const [m, v] of w) {
            if (m >= lo && m <= hi) continue;
            if (Math.abs(v) < 4 * bin) { small.push(v); continue; }
            const k = Math.round(v / bin);
            const next = new Float64Array(nV);
            for (let i = bLo; i <= bHi; i++) {
                const q = 0.5 * base[i];
                if (q === 0) continue;
                next[clampIdx(i - k, nV)] += q;
                next[clampIdx(i + k, nV)] += q;
            }
            base = next;
            bLo = Math.max(0, bLo - Math.abs(k)); bHi = Math.min(nV - 1, bHi + Math.abs(k));
        }
        if (small.length) {
            const f = smallTerms(small, bin);
            const conv = new Float64Array(nV);
            for (let i = bLo; i <= bHi; i++) {
                const q = base[i];
                if (q === 0) continue;
                for (let k = 0; k < f.p.length; k++) conv[clampIdx(i + f.off + k, nV)] += q * f.p[k];
            }
            base = conv;
            bLo = Math.max(0, bLo + f.off); bHi = Math.min(nV - 1, bHi + f.off + f.p.length - 1);
        }

        const run = (cond) => {
            const S = 1 << K;
            let cur = Array.from({ length: S }, () => null), nxt = Array.from({ length: S }, () => null);
            const lim = new Int32Array(2 * S), nlim = new Int32Array(2 * S);
            const bitOf = (state, m, first) => (state >> (m - first)) & 1;
            // Initial states: bits lo .. lo + K - 1.
            for (let st = 0; st < S; st++) {
                let shift = 0, ok = true;
                for (let m = lo; m < lo + K; m++) {
                    const d = bitOf(st, m, lo) ? 1 : -1;
                    if (m === 0 && d !== cond) ok = false;
                    shift += (w.get(m) || 0) * d;
                }
                if (!ok) continue;
                const k = Math.round(shift / bin), a = new Float64Array(nV);
                for (let i = bLo; i <= bHi; i++) if (base[i] !== 0) a[clampIdx(i + k, nV)] += base[i];
                cur[st] = a;
                lim[2 * st] = clampIdx(bLo + k, nV); lim[2 * st + 1] = clampIdx(bHi + k, nV);
            }
            for (let b = bA; b <= bB; b++) {
                const m = b + pre, wm = w.get(m) || 0;
                const r = res.get(b);
                const kern = new Map();
                for (let ns = 0; ns < S; ns++) { if (nxt[ns]) nxt[ns].fill(0); nlim[2 * ns] = nV; nlim[2 * ns + 1] = -1; }
                for (let st = 0; st < S; st++) {
                    const src = cur[st];
                    if (!src) continue;
                    const sLo = lim[2 * st], sHi = lim[2 * st + 1];
                    if (sHi < sLo) continue;
                    for (let v = 0; v < 2; v++) {
                        const d = v ? 1 : -1;
                        if (m === 0 && d !== cond) continue;
                        const W = st | (v << K);
                        // Bits of the window: position t holds d_(lo(b) + t), lo(b) = b - 1 - post.
                        const bit = (q) => ((W >> (q - (b - 1 - post))) & 1) ? 1 : -1;
                        let x = 0;
                        if (fir[0] !== 0) x += fir[0] * (bit(b + 1) - bit(b));
                        x += fir[1] * (bit(b) - bit(b - 1));
                        if (fir[2] !== 0) x += fir[2] * (bit(b - 1) - bit(b - 2));
                        const key = `${x}|${d}`;
                        let kr = kern.get(key);
                        if (!kr) { kr = kernel(s, b, x, wm * d, r); kern.set(key, kr); }
                        const ns = W >> 1;
                        if (!nxt[ns]) nxt[ns] = new Float64Array(nV);
                        const dst = nxt[ns];
                        for (let e = 0; e < kr.off.length; e++) {
                            const o = kr.off[e], q = kr.p[e];
                            const i0 = Math.max(sLo, -o), i1 = Math.min(sHi, nV - 1 - o);
                            for (let i = i0; i <= i1; i++) dst[i + o] += q * src[i];
                            // Beyond the grid: onto its end bins.
                            for (let i = sLo; i < i0; i++) dst[0] += q * src[i];
                            for (let i = Math.max(i1 + 1, sLo); i <= sHi; i++) dst[nV - 1] += q * src[i];
                            nlim[2 * ns] = Math.min(nlim[2 * ns], clampIdx(sLo + o, nV));
                            nlim[2 * ns + 1] = Math.max(nlim[2 * ns + 1], clampIdx(sHi + o, nV));
                        }
                    }
                }
                // Trim negligible ends.
                for (let ns = 0; ns < S; ns++) {
                    const a = nxt[ns];
                    if (!a) continue;
                    let l = nlim[2 * ns], h = nlim[2 * ns + 1];
                    while (l <= h && a[l] < PRUNE) a[l++] = 0;
                    while (h >= l && a[h] < PRUNE) a[h--] = 0;
                    nlim[2 * ns] = l; nlim[2 * ns + 1] = h;
                }
                [cur, nxt] = [nxt, cur];
                lim.set(nlim);
            }
            const out = new Float64Array(nV);
            for (let st = 0; st < S; st++) {
                const a = cur[st];
                if (!a) continue;
                for (let i = lim[2 * st]; i <= lim[2 * st + 1]; i++) out[i] += a[i];
            }
            let tot = 0;
            for (let i = 0; i < nV; i++) tot += out[i];
            for (let i = 0; i < nV; i++) out[i] /= tot;
            // Weak residuals and the baseline wander: one Gaussian.
            const sig = Math.sqrt(lay.noiseVar + (p.sigmaSlow || 0) ** 2);
            return sig > 0.05 * bin || Math.abs(lay.noiseMean) > 0.05 * bin ? gaussConv(out, sig / bin, lay.noiseMean / bin) : out;
        };

        // Residual of boundary b for transition x, with the new bit's linear term wd.
        const kernel = (s_, b, x, wd, r) => {
            const tau = s_ - b * spb;
            const vals = new Map();
            const add = (v, q) => { const k = Math.round(v / bin); vals.set(k, (vals.get(k) || 0) + q); };
            if (!r || x === 0) add(wd, 1);
            else {
                const g = x > 0 ? gUp : gDn, gl = gLin(tau);
                if (!r.jitter || r.jitter < 0.01 * bin) add(wd + x * (g(tau) - gl), 1);
                else for (let i = 0; i < jt.length; i++) add(wd + x * (g(tau - jt[i]) - gl), jn.p[i]);
            }
            const off = Int32Array.from(vals.keys()), pr = Float64Array.from(vals.values());
            return { off, p: pr };
        };

        const p1 = run(1);
        let p0;
        if (symmetric && 2 * c0 === nV - 1) {
            p0 = new Float64Array(nV);
            for (let i = 0; i < nV; i++) p0[i] = p1[nV - 1 - i];
        } else p0 = run(-1);
        return { p1, p0 };
    };

    return { layout, evaluate, symmetric };
}

function clampIdx(i, n) { return i < 0 ? 0 : i >= n ? n - 1 : i; }

// Distribution of a sum of small independent +-w terms on a grid finer than bin, deposited
// onto bin (linear split, the mean kept). Returns { off, p }: p[k] at bin offset off + k.
function smallTerms(ws, bin) {
    let S = 0;
    for (const w of ws) S += Math.abs(w);
    let fine = bin / 16;
    const nMax = 1 << 15;
    if (2 * S / fine > nMax) fine = 2 * S / nMax;
    // Steps of each term on the fine grid; the support is their sum.
    const ks = ws.map(w => Math.round(Math.abs(w) / fine));
    let K = 0, varRest = 0;
    ws.forEach((w, i) => { if (ks[i] === 0) varRest += w * w; else K += ks[i]; });
    const n = 2 * K + 3, c = K + 1;
    let a = new Float64Array(n), b = new Float64Array(n);
    a[c] = 1;
    let lo = c, hi = c;
    for (const k of ks) {
        if (k === 0) continue;
        b.fill(0, lo - k, hi + k + 1);
        for (let i = lo; i <= hi; i++) {
            const q = 0.5 * a[i];
            if (q === 0) continue;
            b[i - k] += q; b[i + k] += q;
        }
        [a, b] = [b, a];
        lo -= k; hi += k;
        while (lo < hi && a[lo] < PRUNE) a[lo++] = 0;
        while (hi > lo && a[hi] < PRUNE) a[hi--] = 0;
    }
    if (varRest > 0) a = gaussConv(a, Math.sqrt(varRest) / fine);
    // Onto the coarse grid: fine index i at (i - c) fine / bin coarse bins.
    const kLo = Math.floor((0 - c) * fine / bin) - 1, kHi = Math.ceil((n - 1 - c) * fine / bin) + 1;
    const out = new Float64Array(kHi - kLo + 1);
    for (let i = 0; i < n; i++) {
        if (a[i] === 0) continue;
        const x = (i - c) * fine / bin - kLo, k = Math.floor(x), f = x - k;
        out[k] += a[i] * (1 - f);
        out[k + 1] += a[i] * f;
    }
    return { off: kLo, p: out };
}

// Convolution with a unit-mass Gaussian of sigma bins and mean mu bins (cut at 8 sigma).
// A sigma below a bin keeps the mean by splitting it over the neighbouring bins.
function gaussConv(a, sig, mu = 0) {
    const n = a.length, k0 = Math.round(mu), h = Math.ceil(8 * sig) + 1;
    const g = new Float64Array(2 * h + 1);
    if (sig > 0.05) {
        for (let k = -h; k <= h; k++) g[k + h] = Phi((k + k0 - mu + 0.5) / sig) - Phi((k + k0 - mu - 0.5) / sig);
    } else {
        const f = mu - k0, j = f >= 0 ? 1 : -1;
        g[h] = 1 - Math.abs(f); g[h + j] = Math.abs(f);
    }
    let gs = 0;
    for (const x of g) gs += x;
    const out = new Float64Array(n);
    for (let i = 0; i < n; i++) {
        const q = a[i];
        if (q === 0) continue;
        for (let k = -h; k <= h; k++) if (g[k + h] !== 0) out[clampIdx(i + k0 + k, n)] += q * g[k + h] / gs;
    }
    return out;
}

/**
 * BER of a decision threshold at the bin boundaries: ber[i] at v0 + (i - 1/2) bin, i = 0..nV,
 * for one and zero bits equally likely. Tails summed from their own ends.
 */
function berCurve(p1, p0) {
    const n = p1.length, ber = new Float64Array(n + 1);
    const c1 = new Float64Array(n + 1), c0 = new Float64Array(n + 1);
    for (let i = 0; i < n; i++) c1[i + 1] = c1[i] + p1[i];          // P(V1 < boundary i+1)
    for (let i = n - 1; i >= 0; i--) c0[i] = c0[i + 1] + p0[i];      // P(V0 > boundary i)
    for (let i = 0; i <= n; i++) ber[i] = 0.5 * (c1[i] + c0[i]);
    return ber;
}

// BER at the level v (linear in log10 BER between bin boundaries).
function berAt(ber, grid, v) {
    const x = (v - grid.v0) / grid.bin + 0.5, i = Math.floor(x), f = x - i;
    if (i < 0) return ber[0];
    if (i >= ber.length - 1) return ber[ber.length - 1];
    const l = (q) => Math.log10(Math.max(q, 1e-300));
    return Math.pow(10, l(ber[i]) + (l(ber[i + 1]) - l(ber[i])) * f);
}

/**
 * Vertical opening at a target BER around the level of least BER: { top, bottom, height }
 * (height 0 and top = bottom = NaN when closed).
 */
function verticalOpening(ber, grid, target) {
    let im = 0;
    for (let i = 1; i < ber.length; i++) if (ber[i] < ber[im]) im = i;
    if (!(ber[im] < target)) return { top: NaN, bottom: NaN, height: 0 };
    const l = (q) => Math.log10(Math.max(q, 1e-300)), lt = Math.log10(target);
    const vAt = (i) => grid.v0 + (i - 0.5) * grid.bin;
    let i = im;
    while (i + 1 < ber.length && ber[i + 1] < target) i++;
    const top = i + 1 < ber.length ? vAt(i) + grid.bin * (lt - l(ber[i])) / (l(ber[i + 1]) - l(ber[i])) : vAt(i);
    i = im;
    while (i - 1 >= 0 && ber[i - 1] < target) i--;
    const bottom = i - 1 >= 0 ? vAt(i) - grid.bin * (lt - l(ber[i])) / (l(ber[i - 1]) - l(ber[i])) : vAt(i);
    return { top, bottom, height: top - bottom };
}

// Rows of the plotted maps.
const PLOT_ROWS = 256;

/**
 * Statistical eye over one UI of sampling phases and its BER-target opening.
 * p: makeStatEye parameters; opts: { ber (target), phases (per UI, even), nV (odd) }.
 * Returns { ber, height (peak within +-0.1 UI of the eye's centre), heightPhase,
 * heightSample (at n0), width, widthLeft, widthRight (phases in UI from n0), berSample
 * (at n0 and the threshold), thr, plot: { x (UI, the sample point at 1), y, density
 * (log10 probability per row), logBer } }.
 */
function statEye(p, opts) {
    const se = makeStatEye(p);
    const { spb } = p, nPh = opts.phases, target = opts.ber;
    // A phase decides bit 0 only inside its own eye, which need not be centred on n0 (a
    // flat pulse top puts n0 at its start): a coarse scan over +-1 UI finds the centre.
    // The scan also sizes the voltage grid.
    const coarse = Array.from({ length: 17 }, (_, k) => (k - 8) / 8);
    let bound = 0;
    for (const u of coarse) bound = Math.max(bound, se.layout(p.n0 + u * spb, 0).bound);
    bound = 1.05 * bound + 8 * (p.sigmaSlow || 0);
    if (!(bound > 0)) bound = 1e-3;
    const nV = opts.nV, bin = 2 * bound / (nV - 1);
    const grid = { v0: p.vmid - (nV - 1) / 2 * bin, bin, nV };
    const thr = p.vmid;
    const at = (u) => {
        const { p1, p0 } = se.evaluate(u * spb, grid);
        const ber = berCurve(p1, p0);
        return { p1, p0, ber, atThr: berAt(ber, grid, thr) };
    };
    const span = (list, open) => {
        // Contiguous run of open entries around the least BER: [first, last] indices.
        let k = 0;
        for (let i = 1; i < list.length; i++) if (list[i].atThr < list[k].atThr) k = i;
        if (!open(list[k])) return null;
        let l = k, r = k;
        while (l > 0 && open(list[l - 1])) l--;
        while (r < list.length - 1 && open(list[r + 1])) r++;
        return { l, r, k };
    };
    // With a DFE the slicer holds bit 0's feedback only within half a UI of the sample point
    // n0 (as the time-domain slicer does), so the eye lies in that window.
    const dfeWindow = (p.dfe || []).length > 0;
    let centre = 0;
    if (!dfeWindow) {
        const cs = coarse.map(u => ({ u, atThr: at(u).atThr }));
        const cr = span(cs, c => c.atThr < target);
        // The centre of the open phases at the threshold (edges interpolated in log BER),
        // else the least BER there.
        const lb = (q) => Math.log10(Math.max(q, 1e-300)), lt = Math.log10(target);
        const cross = (i, o) => {
            if (o < 0 || o >= cs.length) return cs[i].u;
            const a = lb(cs[i].atThr), b = lb(cs[o].atThr);
            return cs[i].u + (cs[o].u - cs[i].u) * Math.min(1, Math.max(0, (lt - a) / (b - a)));
        };
        centre = cr ? 0.5 * (cross(cr.l, cr.l - 1) + cross(cr.r, cr.r + 1)) : cs[span(cs, () => true).k].u;
        centre = Math.round(centre * nPh) / nPh;
    }

    const cols = [];
    for (let k = 0; k < nPh; k++) {
        const u = centre + (k - nPh / 2) / nPh, c = at(u);
        const dens = new Float64Array(nV);
        for (let i = 0; i < nV; i++) dens[i] = 0.5 * (c.p1[i] + c.p0[i]);
        cols.push({ u, ber: c.ber, dens, open: verticalOpening(c.ber, grid, target), atThr: c.atThr });
    }
    // Width at the threshold: the open phases around the least BER, edges by bisection.
    let width = 0, wl = NaN, wr = NaN;
    const ws = span(cols, c => c.atThr < target);
    if (ws) {
        const edge = (inside, outside) => {
            let a = inside, b = outside;
            for (let it = 0; it < 7; it++) {
                const m = 0.5 * (a + b);
                if (at(m).atThr < target) a = m; else b = m;
            }
            return 0.5 * (a + b);
        };
        // Past the end of a DFE window the next decision's feedback applies: the window ends it.
        const step = 1 / nPh;
        wl = dfeWindow && ws.l === 0 ? cols[0].u : edge(cols[ws.l].u, cols[ws.l].u - step);
        wr = dfeWindow && ws.r === nPh - 1 ? centre + 0.5 : edge(cols[ws.r].u, cols[ws.r].u + step);
        width = Math.min(1, wr - wl);
    }
    // Peak height within +-0.1 UI of the centre of the open phases.
    const mid = ws ? 0.5 * (wl + wr) : centre;
    let height = 0, hPhase = NaN;
    for (const c of cols) {
        if (Math.abs(c.u - mid) > 0.1 + 1e-9) continue;
        if (c.open.height > height) { height = c.open.height; hPhase = c.u; }
    }
    // At the sample point n0.
    const s0 = at(0), o0 = verticalOpening(s0.ber, grid, target);

    // Maps over 2 UI, the eye centre at 1 UI, rows over the populated voltages. The
    // neighbouring eyes repeat it (stationary bits; with a DFE the seams are its switching
    // times).
    let iLo = nV, iHi = -1;
    for (const c of cols) for (let i = 0; i < nV; i++) if (c.dens[i] > 1e-18) { iLo = Math.min(iLo, i); iHi = Math.max(iHi, i); }
    const pad = Math.max(2, Math.round(0.04 * (iHi - iLo)));
    iLo = Math.max(0, iLo - pad); iHi = Math.min(nV - 1, iHi + pad);
    const gsz = Math.max(1, Math.ceil((iHi - iLo + 1) / PLOT_ROWS)), rows = Math.ceil((iHi - iLo + 1) / gsz);
    const y = Array.from({ length: rows }, (_, r) => grid.v0 + (iLo + (r + 0.5) * gsz - 0.5) * bin);
    const x = Array.from({ length: 2 * nPh + 1 }, (_, j) => j / nPh);
    const density = [], logBer = [];
    for (let r = 0; r < rows; r++) {
        const dr = new Array(x.length), br = new Array(x.length);
        for (let j = 0; j < x.length; j++) {
            const c = cols[((j - nPh + nPh / 2) % nPh + nPh) % nPh];
            let d = 0;
            for (let i = iLo + r * gsz; i < Math.min(iLo + (r + 1) * gsz, nV); i++) d += c.dens[i];
            dr[j] = d > 1e-30 ? Math.log10(d) : null;
            br[j] = Math.log10(Math.max(berAt(c.ber, grid, y[r]), 1e-30));
        }
        density.push(dr); logBer.push(br);
    }
    return {
        ber: target, thr, bin, centre, height, heightPhase: hPhase, heightSample: o0.height,
        sampleTop: o0.top, sampleBottom: o0.bottom,
        width, widthLeft: wl, widthRight: wr, berSample: s0.atThr, symmetric: se.symmetric,
        open: cols.map(c => ({ u: c.u, top: c.open.top, bottom: c.open.bottom })),
        plot: { x, y, density, logBer },
    };
}

export { makeStatEye, statEye, jitterNodes, berCurve, berAt, verticalOpening, smallTerms, gaussConv };
