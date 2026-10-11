// Periodic steady state of a repeating bit pattern, and its eye.
//
// The fast part of the link (exact transfer minus the slow AC-coupling modes) responds
// to every edge of a leg with the same step response, rising or falling. Those step
// responses come from one FFT window each; their sample differences ds (zero once
// settled) are convolved, circularly over the pattern period, with the trains of edge
// amplitudes on the bit boundaries. The cumulative sum is the periodic response up to
// a constant. The slow modes are integrated exactly over one period with the initial
// state that makes them periodic. The mean of everything is the DC operating point at
// the legs' mean levels.
import { legMean, edgeSpectrum, edgeValue } from './stimulus.js';
import { makeTransfer, dcSolve } from './network.js';
import { smoothSize } from './fft.js';
import { fmtHz } from './channel_grid.js';
import { lineEstimates, reflectionTrips, taper, deriveTraces, firstOrder, slowSplit, packOutput,
    N_MAX, fmtTime, fmtV, fmtPct, mainName, srcName } from './transient.js';

// Samples of one pattern period.
const M_MAX = 1 << 20;
// Bits of the waveform shown (the eye uses the whole period).
const SHOW_BITS = 200;
// Eye histogram bins over 2 UI and over the voltage range.
const EYE_T_BINS = 320, EYE_V_BINS = 240;

/**
 * Step responses of every (output, leg, edge direction) on a window of N samples.
 * offsets[l]: the 50 % time of leg l's edge in the window. Returns { s (Float64Array per
 * response, index (o * nl + l) * 2 + e, e = 0 rising, 1 falling), tail, extrapFraction,
 * hplot }.
 */
function edgeResponses(ctx, transfer, slowRes, N, dt, offsets) {
    const { stim, net, ev, out, fft, fTop, diff } = ctx;
    const nl = net.nl, no = out.nodes.length, nr = no * nl * 2;
    const K = N / 2, df = 1 / (N * dt), fNyq = K * df;
    const spec = new Float64Array(nr * (K + 1) * 2);
    const O0 = new Float64Array(nr);
    const scratch = { re: new Float64Array(no * nl), im: new Float64Array(no * nl) };
    const plotStride = Math.max(1, Math.floor(K / 2000));
    const hf = [], hm = [], hp = [], sm = [];
    let eTot = 0, eExt = 0;
    for (let k = 0; k <= K; k++) {
        const f = k * df, w = 2 * Math.PI * f;
        const H = transfer(f, ev(f), scratch);
        if (k % plotStride === 0 && k > 0) {
            let re, im;
            if (!diff) { re = H.re[0]; im = H.im[0]; } else {
                re = 0.5 * (H.re[0] - H.re[1] - H.re[2] + H.re[3]);
                im = 0.5 * (H.im[0] - H.im[1] - H.im[2] + H.im[3]);
            }
            hf.push(f); hm.push(Math.hypot(re, im)); hp.push(Math.atan2(im, re));
            sm.push(Math.abs(edgeSpectrum(stim.up, f)));
        }
        for (const s of slowRes) {
            const d = 1 + w * w * s.lambda * s.lambda;
            const cr = 1 / d, ci = -w * s.lambda / d;
            for (let i = 0; i < no * nl; i++) { H.re[i] -= s.r[i] * cr; H.im[i] -= s.r[i] * ci; }
        }
        if (k === 0) {
            for (let o = 0; o < no; o++) for (let l = 0; l < nl; l++) for (let e = 0; e < 2; e++) {
                O0[(o * nl + l) * 2 + e] = H.re[o * nl + l];
            }
            continue;
        }
        if (k === K) continue;
        const W = taper(f, fNyq);
        if (W === 0) continue;
        const E = [edgeSpectrum(stim.up, f), edgeSpectrum(stim.down, f)];
        eTot += E[0] * E[0] * W * W;
        if (f > fTop) eExt += E[0] * E[0] * W * W;
        for (let l = 0; l < nl; l++) {
            const ph = -w * offsets[l], c = Math.cos(ph), sn = Math.sin(ph);
            for (let e = 0; e < 2; e++) {
                // Edge derivative spectrum E e^(-j w t_l) W, then O'(f) / (j w dt).
                const ur = E[e] * W * c, ui = E[e] * W * sn, sc = 1 / (w * dt);
                for (let o = 0; o < no; o++) {
                    const hr = H.re[o * nl + l], hi = H.im[o * nl + l];
                    const sr = hr * ur - hi * ui, si = hr * ui + hi * ur;
                    const base = 2 * (((o * nl + l) * 2 + e) * (K + 1) + k);
                    spec[base] = si * sc; spec[base + 1] = -sr * sc;
                }
            }
        }
    }
    const p = fft.irfftBatch(N, nr, spec);
    const s = [];
    let eTail = 0, eAll = 0;
    for (let r = 0; r < nr; r++) {
        const v = new Float64Array(N), p0 = p[r * N];
        for (let n = 0; n < N; n++) v[n] = O0[r] * n / N + p[r * N + n] - p0;
        for (let n = 1; n < N; n++) {
            const d = (v[n] - v[n - 1]) ** 2;
            eAll += d;
            if (n > 0.9 * N) eTail += d;
        }
        s.push(v);
    }
    return {
        s, tail: eAll > 0 ? eTail / eAll : 0, extrapFraction: eTot > 0 ? eExt / eTot : 0,
        hplot: { f: Float64Array.from(hf), mag: Float64Array.from(hm), phase: Float64Array.from(hp), stim: Float64Array.from(sm) },
    };
}

// One period of a driver leg on the grid n * dt (edges at k * UI + shift).
function sampleLegPeriodic(stim, leg, M, dt) {
    const P = M * dt, u = new Float64Array(M);
    // Level between the edges: start from v0 at t = 0 and step at every 50 % point.
    // v0 is the level before the first boundary; an edge moved before t = 0 by a negative
    // skew has happened by t = 0 (and comes again, a period later, inside the period).
    let level = leg.v0, j = 0;
    for (const e of leg.edges) if (e.t < 0) level += e.dv;
    const edges = leg.edges.map(e => ({ ...e, t: ((e.t % P) + P) % P })).sort((a, b) => a.t - b.t);
    for (let n = 0; n < M; n++) {
        while (j < edges.length && edges[j].t <= n * dt) level += edges[j++].dv;
        u[n] = level;
    }
    // The transitions: within half an edge of each 50 % point, periodic.
    for (const e of edges) {
        const h = e.edge.half, n0 = Math.ceil((e.t - h) / dt), n1 = Math.floor((e.t + h) / dt);
        for (let n = n0; n <= n1; n++) {
            const m = ((n % M) + M) % M, tt = n * dt;
            // u[m] holds the level on its side of this edge; replace by the shaped value.
            const after = tt >= e.t;
            u[m] += e.dv * (edgeValue(e.edge, tt - e.t) - (after ? 1 : 0));
        }
    }
    return u;
}

// Periodic steady state of a slow mode driven by one period of s (exact for s linear
// between samples).
function firstOrderPeriodic(lambda, s, dt) {
    const M = s.length;
    const ext = new Float64Array(M + 1);
    ext.set(s); ext[M] = s[0];
    const z = firstOrder(lambda, ext, dt);
    // z(t) + z0 exp(-t / lambda) is periodic for z0 = z(P) / (1 - exp(-P / lambda)).
    const z0 = z[M] / -Math.expm1(-M * dt / lambda);
    const out = new Float64Array(M);
    for (let n = 0; n < M; n++) out[n] = z[n] + z0 * Math.exp(-n * dt / lambda);
    return out;
}

function patternSim(ctx, dt0) {
    const { stim, net, est, warnings, ev, len, out, diff, fft, fTop } = ctx;
    const nl = net.nl, no = out.nodes.length;
    const ui = stim.ui, nb = stim.bits.length;
    let spb = Math.max(8, Math.ceil(ui / dt0 - 1e-9));
    if (nb * spb > M_MAX) {
        spb = Math.floor(M_MAX / nb);
        if (spb < 8) throw new Error(`The pattern is too long: at most ${M_MAX / 8} bits.`);
        if (spb < 16) warnings.push(`${nb} bits fit only ${spb} samples per bit in ${M_MAX} samples: the edges are coarsely sampled.`);
    }
    const dt = ui / spb, M = nb * spb;

    // Step-response window: a quiet lead, the line delay and the reflections.
    const skew = stim.legs.length > 1 ? stim.legs[1].edges[0].t - stim.legs[0].edges[0].t : 0;
    const nLead = Math.ceil((64 * dt + stim.halfMax + Math.abs(skew)) / dt);
    const offsets = stim.legs.map((leg, l) => nLead * dt + (l === 1 ? skew : 0));
    const tauMax = Math.max(...est.modes.map(m => m.tau)), tauMin = Math.min(...est.modes.map(m => m.tau));
    const trips = reflectionTrips(lineEstimates(ev, 1 / (8 * tauMax), len, diff), net);
    const edgeLong = Math.max(stim.up.t1090, stim.down.t1090);
    const tBase = nLead * dt + Math.abs(skew) + tauMax * (1 + 2 * trips) + 20 * edgeLong + 2 * stim.halfMax;
    const split = slowSplit(ctx, tBase, dt);
    let T = split.T;

    // Main step (P rising, N falling), relative response before the first arrival.
    const iRes = (o, l, e) => (o * nl + l) * 2 + e;
    const mainStep = (s) => {
        if (!diff) return (n) => s[iRes(0, 0, 0)][n];
        return (n) => (s[iRes(0, 0, 0)][n] - s[iRes(1, 0, 0)][n]) - (s[iRes(0, 1, 1)][n] - s[iRes(1, 1, 1)][n]);
    };
    const precursor = (v, N) => {
        const iPre = Math.min(N, Math.floor((nLead * dt - stim.halfMax + 0.9 * tauMin) / dt));
        let pre = 0, swing = 0;
        for (let n = 0; n < N; n++) swing = Math.max(swing, Math.abs(v(n) - v(0)));
        for (let n = 0; n < iPre; n++) pre = Math.max(pre, Math.abs(v(n) - v(0)));
        return swing > 0 ? pre / swing : 0;
    };

    const transfer = makeTransfer(net, len, out.nodes);
    let resp = null, N = 0, nonCausal = false, prevPre = Infinity;
    for (let attempt = 0; attempt < 5; attempt++) {
        N = Math.min(N_MAX, smoothSize(T / dt));
        resp = edgeResponses(ctx, transfer, split.slowRes, N, dt, offsets);
        if (resp.tail <= 1e-6 || 2 * T > split.tCap || N >= N_MAX) break;
        const pre = precursor(mainStep(resp.s), N);
        nonCausal = pre > 2e-3 && pre > 0.7 * prevPre;
        if (nonCausal) break;
        prevPre = pre;
        T *= 2;
    }
    if (resp.tail > 1e-6 && !nonCausal) {
        warnings.push('The step responses have not settled within their window (reflections still ringing); the eye is not accurate.');
    }
    const preRel = precursor(mainStep(resp.s), N);
    if (preRel > 2e-3) {
        warnings.push(`The receiver responds ${(100 * preRel).toFixed(2)} % of the swing before the signal can arrive`
            + (ctx.cfg.causal ? '.' : ': enable causal materials for a causal dielectric model.'));
    }
    if (resp.extrapFraction > 0.01) {
        warnings.push(`${(100 * resp.extrapFraction).toFixed(1)} % of the edge spectrum is above the sweep stop `
            + `${fmtHz(fTop)}, where the line parameters are extrapolated.`);
    }

    // Circular convolution over the period: edge trains (impulses on the bit boundaries)
    // with the folded step differences.
    const L = smoothSize(2 * M), K = L / 2;
    const acc = Array.from({ length: no }, () => new Float64Array(2 * (K + 1)));
    const buf = new Float64Array(L);
    for (let l = 0; l < nl; l++) {
        const leg = stim.legs[l];
        // Bit boundary of each edge (the leg's skew is in its responses).
        const boundary = (e) => Math.round((e.t - (l === 1 ? skew : 0) - stim.t0) / ui);
        for (let e = 0; e < 2; e++) {
            buf.fill(0);
            let any = false;
            for (const ed of leg.edges) {
                if ((ed.dv > 0) !== (e === 0)) continue;
                buf[(((boundary(ed) % nb) + nb) % nb) * spb] += ed.dv;
                any = true;
            }
            if (!any) continue;
            const A = fft.rfft(buf);
            for (let o = 0; o < no; o++) {
                const s = resp.s[iRes(o, l, e)];
                buf.fill(0);
                for (let n = 0; n < N; n++) buf[n % M] += n === 0 ? s[0] : s[n] - s[n - 1];
                const D = fft.rfft(buf), a = acc[o];
                for (let k = 0; k <= K; k++) {
                    const ar = A[2 * k], ai = A[2 * k + 1], dr = D[2 * k], di = D[2 * k + 1];
                    a[2 * k] += ar * dr - ai * di;
                    a[2 * k + 1] += ar * di + ai * dr;
                }
            }
        }
    }
    const conv = fft.irfftBatch(L, no, (() => {
        const all = new Float64Array(no * 2 * (K + 1));
        acc.forEach((a, o) => all.set(a, o * 2 * (K + 1)));
        return all;
    })());

    // Driver legs over the period: bit boundaries at k * UI from t = 0.
    const src = stim.legs.map((leg, l) => sampleLegPeriodic(stim, {
        v0: leg.v0, edges: leg.edges.map(e => ({ ...e, t: e.t - stim.t0 })),
    }, M, dt));
    const xdc = dcSolve(net, ctx.Gdc, stim.legs.map(leg => legMean(stim, leg)));
    const o = {};
    for (let k = 0; k < no; k++) {
        // Wrap the linear convolution to the period, integrate, align the edges to
        // t = k * UI, remove the mean.
        const c = new Float64Array(M);
        for (let n = 0; n < L; n++) c[n % M] += conv[k * L + n];
        // G[m] = sum of c up to m is periodic (the edges of a period sum to zero); the
        // responses put edge k at index k * spb + nLead.
        const G = new Float64Array(M);
        let run = 0, mean = 0;
        for (let m = 0; m < M; m++) { run += c[m]; G[m] = run; mean += run / M; }
        const v = new Float64Array(M), base = xdc[out.nodes[k]];
        for (let n = 0; n < M; n++) v[n] = base + G[(n + nLead) % M] - mean;
        o[out.names[k]] = v;
    }
    // Slow modes: periodic steady state, without their mean (in the DC point).
    for (const sm of split.slowRes) {
        for (let l = 0; l < nl; l++) {
            const z = firstOrderPeriodic(sm.lambda, src[l], dt);
            let zm = 0;
            for (let n = 0; n < M; n++) zm += z[n] / M;
            for (let k = 0; k < no; k++) {
                const r = sm.r[k * nl + l];
                if (r === 0) continue;
                const v = o[out.names[k]];
                for (let n = 0; n < M; n++) v[n] += r * (z[n] - zm);
            }
        }
    }

    const full = deriveTraces(diff, o, src);
    const main = full[mainName(diff)];
    const eye = eyeAnalysis(main, dt, ui, spb, nb);
    const delayMetric = eyeDelay(eye, full[srcName(diff)], dt, ui, tauMin);

    // The shown waveform: the first bits of the period.
    const nShow = Math.min(M, SHOW_BITS * spb);
    const t = new Float64Array(nShow);
    for (let n = 0; n < nShow; n++) t[n] = n * dt;
    const traces = Object.fromEntries(Object.entries(full).map(([k, v]) => [k, v.slice(0, nShow)]));
    const metrics = { list: [
        { label: 'Bit rate / UI', value: `${+(stim.bitRate / 1e9).toPrecision(4)} Gb/s / ${fmtTime(ui)}` },
        { label: 'Pattern length', value: `${nb} bits` },
        ...eye.metrics, delayMetric,
    ].filter(Boolean), crossTime: NaN };
    const res = {
        kind: 'pattern', t, traces, metrics, hplot: resp.hplot,
        info: { dt, N, window: N * dt, slowModes: split.slow.length, trips, samplesPerBit: spb, periodSamples: M },
        eyeHist: eye.hist, eyeInfo: eye.info,
    };
    packOutput(res);
    return res;
}

// Eye of one period of the main receiver trace: crossings at the threshold half way
// between the mean one and zero levels at the eye centre, the centre half a UI from the
// mean crossing phase. Returns { metrics, hist, info }.
function eyeAnalysis(v, dt, ui, spb, nb) {
    const M = v.length;
    const at = (t) => {
        const x = ((t / dt) % M + M) % M, n = Math.floor(x), f = x - n;
        return v[n] + (v[(n + 1) % M] - v[n]) * f;
    };
    const crossings = (thr) => {
        const ph = [];
        for (let n = 0; n < M; n++) {
            const a = v[n] - thr, b = v[(n + 1) % M] - thr;
            if ((a < 0) !== (b < 0)) ph.push(((n + a / (a - b)) * dt) % ui);
        }
        let C = 0, S = 0;
        for (const p of ph) { C += Math.cos(2 * Math.PI * p / ui); S += Math.sin(2 * Math.PI * p / ui); }
        const pc = ((Math.atan2(S, C) / (2 * Math.PI)) * ui + ui) % ui;
        const dev = ph.map(p => { let d = p - pc; d -= ui * Math.round(d / ui); return d; });
        return { pc, dev };
    };
    const levels = (pc, thr) => {
        const tc = pc + ui / 2, ones = [], zeros = [];
        for (let k = 0; k < nb; k++) (at(k * ui + tc) > thr ? ones : zeros).push(k);
        return { tc, ones, zeros };
    };
    let mean = 0;
    for (const x of v) mean += x / M;
    let thr = mean, cr = crossings(thr), lv = levels(cr.pc, thr);
    for (let it = 0; it < 2 && lv.ones.length && lv.zeros.length; it++) {
        const m1 = lv.ones.reduce((s, k) => s + at(k * ui + lv.tc), 0) / lv.ones.length;
        const m0 = lv.zeros.reduce((s, k) => s + at(k * ui + lv.tc), 0) / lv.zeros.length;
        thr = 0.5 * (m1 + m0);
        cr = crossings(thr);
        lv = levels(cr.pc, thr);
    }
    const metrics = [];
    const info = { thr, centre: 1, ui };
    if (!lv.ones.length || !lv.zeros.length || !cr.dev.length) {
        metrics.push({ label: 'Eye', value: 'no transitions at the receiver' });
        return { metrics, hist: histogram(v, dt, ui, 0), info };
    }
    // Levels and inner opening at a phase, for the bits decided at the centre.
    const opening = (t) => {
        let min1 = Infinity, max0 = -Infinity, s1 = 0, s0 = 0;
        for (const k of lv.ones) { const x = at(k * ui + t); min1 = Math.min(min1, x); s1 += x; }
        for (const k of lv.zeros) { const x = at(k * ui + t); max0 = Math.max(max0, x); s0 += x; }
        return { min1, max0, m1: s1 / lv.ones.length, m0: s0 / lv.zeros.length };
    };
    const c = opening(lv.tc);
    // The centre wins ties (a flat top).
    let best = c.min1 - c.max0, bestOff = 0;
    for (let j = -32; j <= 32; j++) {
        const off = j * ui / 64, op = opening(lv.tc + off), h = op.min1 - op.max0;
        if (h > best + 1e-9 * Math.abs(c.m1 - c.m0)) { best = h; bestOff = off; }
    }
    let dmin = Infinity, dmax = -Infinity, ss = 0;
    for (const d of cr.dev) { dmin = Math.min(dmin, d); dmax = Math.max(dmax, d); ss += d * d; }
    const jpp = dmax - dmin, jrms = Math.sqrt(ss / cr.dev.length);
    const height = c.min1 - c.max0, amp = c.m1 - c.m0;
    metrics.push({ label: 'Eye height (centre)', value: height > 0 ? `${fmtV(height)} (${fmtPct(height / amp)})` : 'closed', note: 'percentage of the eye amplitude' });
    metrics.push({ label: 'Eye height (best phase)', value: best > 0 ? `${fmtV(best)} at ${bestOff >= 0 ? '+' : ''}${fmtTime(bestOff)}` : 'closed',
        note: 'largest inner opening over the sampling phase, relative to the centre' });
    metrics.push({ label: 'Eye width (at threshold)', value: jpp < ui ? `${fmtTime(ui - jpp)} (${fmtPct(1 - jpp / ui)} UI)` : 'closed' });
    metrics.push({ label: 'Crossing jitter pk-pk / rms', value: `${fmtTime(jpp)} / ${fmtTime(jrms)}`, note: 'data-dependent (ISI) jitter at the threshold' });
    metrics.push({ label: 'One / zero level (centre)', value: `${fmtV(c.m1)} / ${fmtV(c.m0)}` });
    metrics.push({ label: 'Eye amplitude', value: fmtV(amp) });
    // Window of the histogram: the mean crossing at 0.5 UI, the centre at 1 UI.
    const start = cr.pc - ui / 2;
    Object.assign(info, {
        crossPhase: cr.pc, centreTime: lv.tc, start, heightTop: c.min1, heightBottom: c.max0,
        widthLeft: 0.5 + dmax / ui, widthRight: 1.5 + dmin / ui,
    });
    return { metrics, hist: histogram(v, dt, ui, start), info };
}

// 2 UI density histogram of a periodic trace folded every UI, from `start`. Every bit
// gives every time column the same weight, spread over the voltage bins the trace
// crosses within the column (its values at the column edges), so the column totals are
// uniform whatever the samples per bit, and steep edges stay continuous.
function histogram(v, dt, ui, start) {
    const M = v.length, nb = Math.max(1, Math.round(M * dt / ui));
    let lo = Infinity, hi = -Infinity;
    for (const x of v) { lo = Math.min(lo, x); hi = Math.max(hi, x); }
    const pad = 0.06 * (hi - lo || 1);
    const y0 = lo - pad, y1 = hi + pad;
    const counts = new Float32Array(EYE_T_BINS * EYE_V_BINS);
    const half = EYE_T_BINS / 2;
    const at = (t) => {
        const x = ((t / dt) % M + M) % M, n = Math.floor(x), f = x - n;
        return v[n] + (v[(n + 1) % M] - v[n]) * f;
    };
    const bin = (y) => Math.min(EYE_V_BINS - 1, Math.max(0, Math.floor((y - y0) / (y1 - y0) * EYE_V_BINS)));
    const edge = new Float64Array(half + 1);
    for (let k = 0; k < nb; k++) {
        const t0 = start + k * ui;
        for (let j = 0; j <= half; j++) edge[j] = at(t0 + j * ui / half);
        for (let ix = 0; ix < half; ix++) {
            let a = bin(edge[ix]), b = bin(edge[ix + 1]);
            if (a > b) { const tmp = a; a = b; b = tmp; }
            const w = 1 / (b - a + 1);
            for (let iy = a; iy <= b; iy++) {
                counts[iy * EYE_T_BINS + ix] += w;
                counts[iy * EYE_T_BINS + ix + half] += w;
            }
        }
    }
    return { counts, tBins: EYE_T_BINS, vBins: EYE_V_BINS, y0, y1 };
}

// Source-to-receiver delay of the data: the receiver's mean crossing phase against the
// source's, unwrapped with the line delay estimate.
function eyeDelay(eye, src, dt, ui, tau) {
    if (!isFinite(eye.info.crossPhase)) return null;
    let lo = Infinity, hi = -Infinity;
    for (const x of src) { lo = Math.min(lo, x); hi = Math.max(hi, x); }
    const mid = 0.5 * (lo + hi), M = src.length;
    let C = 0, S = 0;
    for (let n = 0; n < M; n++) {
        const a = src[n] - mid, b = src[(n + 1) % M] - mid;
        if ((a < 0) !== (b < 0)) {
            const p = ((n + a / (a - b)) * dt) % ui;
            C += Math.cos(2 * Math.PI * p / ui); S += Math.sin(2 * Math.PI * p / ui);
        }
    }
    const ps = ((Math.atan2(S, C) / (2 * Math.PI)) * ui + ui) % ui;
    let d = ((eye.info.crossPhase - ps) % ui + ui) % ui;
    d += Math.max(0, Math.round((tau - d) / ui)) * ui;
    return { label: 'Delay (50 %)', value: fmtTime(d), note: 'source to receiver crossing, unwrapped with the line delay' };
}

export { patternSim, eyeAnalysis, sampleLegPeriodic, firstOrderPeriodic };
