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
//
// Receiver equalization acts on the main receiver signal (RX, or RX P - RX N): the CTLE
// output is a further output whose step responses are the main ones through the CTLE
// transfer; an ideal DFE then subtracts, over each bit, its post-cursors times the
// previous transmitted bits. TX preset and CTLE gain can be chosen automatically by the
// peak-distortion eye of the pulse response (eq.js).
import { legMean, edgeSpectrum, edgeValue, buildStimulus } from './stimulus.js';
import { ctleOf, pdaEye, dfeTap, PCIE_PRESETS, firFromPresetDb, presetDbFromFir, boostDb } from './eq.js';
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
 * hplot, main }: main holds the spectra of the main receiver signal per (leg, edge) for
 * receiver equalization (mainSteps).
 */
function edgeResponses(ctx, transfer, slowRes, N, dt, offsets) {
    const { stim, net, ev, out, fft, fTop, diff } = ctx;
    const nl = net.nl, no = out.nodes.length, nr = no * nl * 2;
    const K = N / 2, df = 1 / (N * dt), fNyq = K * df;
    const spec = new Float64Array(nr * (K + 1) * 2);
    const O0 = new Float64Array(nr);
    // Main signal (RX, or RX P - RX N) per (leg, edge): the same spectra combined.
    const mspec = new Float64Array(nl * 2 * (K + 1) * 2), mO0 = new Float64Array(nl * 2);
    const mainRow = (arrRe, arrIm, l) => diff
        ? [arrRe[l] - arrRe[nl + l], arrIm[l] - arrIm[nl + l]] : [arrRe[l], arrIm[l]];
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
            for (let l = 0; l < nl; l++) for (let e = 0; e < 2; e++) mO0[l * 2 + e] = mainRow(H.re, H.im, l)[0];
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
                const [hr, hi] = mainRow(H.re, H.im, l);
                const sr = hr * ur - hi * ui, si = hr * ui + hi * ur;
                const base = 2 * ((l * 2 + e) * (K + 1) + k);
                mspec[base] = si * sc; mspec[base + 1] = -sr * sc;
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
        main: { spec: mspec, O0: mO0, N, K, df },
        s, tail: eAll > 0 ? eTail / eAll : 0, extrapFraction: eTot > 0 ? eExt / eTot : 0,
        hplot: { f: Float64Array.from(hf), mag: Float64Array.from(hm), phase: Float64Array.from(hp), stim: Float64Array.from(sm) },
    };
}

// Step responses of the main signal through a filter C (f -> [re, im], C(0) real) per
// (leg, edge), index l * 2 + e. C null: the main signal itself.
function mainSteps(ctx, main, C) {
    const { N, K, df } = main, nr = main.O0.length;
    const spec = new Float64Array(main.spec.length), O0 = new Float64Array(nr);
    for (let r = 0; r < nr; r++) {
        O0[r] = main.O0[r] * (C ? C(0)[0] : 1);
        for (let k = 1; k < K; k++) {
            const i = 2 * (r * (K + 1) + k);
            const [cr, ci] = C ? C(k * df) : [1, 0];
            const a = main.spec[i], b = main.spec[i + 1];
            spec[i] = a * cr - b * ci; spec[i + 1] = a * ci + b * cr;
        }
    }
    const p = ctx.fft.irfftBatch(N, nr, spec);
    const out = [];
    for (let r = 0; r < nr; r++) {
        const v = new Float64Array(N), p0 = p[r * N];
        for (let n = 0; n < N; n++) v[n] = O0[r] * n / N + p[r * N + n] - p0;
        out.push(v);
    }
    return out;
}

// Pulse response q of the main signal (signal = base + sum d_m q(t - m UI), d = +-1)
// from its step responses, with transmitter FIR taps [c-1, c0, c+1].
function pulseResponse(stim, steps, spb, fir) {
    const N = steps[0].length, nl = steps.length / 2, dV = stim.dV;
    const p = new Float64Array(N);
    const at = (a, n) => n >= 0 ? a[n] : 0;
    for (let n = 0; n < N; n++) {
        // One isolated bit: P rises then falls a UI later; N the opposite.
        let v = dV * (steps[0][n] - at(steps[1], n - spb));
        if (nl > 1) v += -dV * (steps[2 * 1 + 1][n] - at(steps[2 * 1 + 0], n - spb));
        p[n] = 0.5 * v;
    }
    const q = new Float64Array(N);
    for (let n = 0; n < N; n++) {
        // Level of bit n = sum_j c_j d_(n-j): c-1 takes the next bit (one UI earlier).
        q[n] = fir[0] * (n + spb < N ? p[n + spb] : 0) + fir[1] * p[n] + fir[2] * (n - spb >= 0 ? p[n - spb] : 0);
    }
    return q;
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
    const tauLine = Math.max(...est.modes.map(m => m.tau));
    const trips = reflectionTrips(lineEstimates(ev, 1 / (8 * tauLine), len, diff), net);
    // Delays of the path: the line modes plus the packages.
    const tauMax = tauLine + net.pkgDelay, tauMin = Math.min(...est.modes.map(m => m.tau)) + net.pkgDelay;
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

    // Settling of a response: the share of its derivative energy in the last 10 %.
    const tailOf = (v, N) => {
        let eAll = 0, eTail = 0;
        for (let n = 1; n < N; n++) { const d = (v(n) - v(n - 1)) ** 2; eAll += d; if (n > 0.9 * N) eTail += d; }
        return eAll > 0 ? eTail / eAll : 0;
    };
    // Amplitude settling, for the common mode: a sharp skew pulse dominates its energy, so
    // a slow ring after it passes the energy measure. Largest deviation from the final
    // value in the last 10 %, relative to the peak.
    const tailAmp = (v, N) => {
        let peak = 0, dev = 0;
        const fin = v(N - 1);
        for (let n = 0; n < N; n++) { peak = Math.max(peak, Math.abs(v(n))); if (n > 0.9 * N) dev = Math.max(dev, Math.abs(v(n) - fin)); }
        return peak > 0 ? dev / peak : 0;
    };
    const CM_SETTLED = 1e-4;
    // Common mode of a P rise with an N fall: zero for a balanced pair; with an open
    // common mode at both ends (a current driver into a differential termination) it
    // rings long, which only matters when the drive excites it.
    const cmStep = (s) => (n) => 0.5 * ((s[iRes(0, 0, 0)][n] + s[iRes(1, 0, 0)][n]) - (s[iRes(0, 1, 1)][n] + s[iRes(1, 1, 1)][n]));

    const cmExcited = (s, N) => {
        const cm = cmStep(s), ms = mainStep(s);
        let cmMax = 0, mMax = 0;
        for (let n = 0; n < N; n++) { cmMax = Math.max(cmMax, Math.abs(cm(n))); mMax = Math.max(mMax, Math.abs(ms(n))); }
        return cmMax > 0.01 * mMax;
    };

    const transfer = makeTransfer(net, len, out.nodes);
    let resp = null, N = 0, nonCausal = false, prevPre = Infinity, tail = 0;
    for (let attempt = 0; attempt < 5; attempt++) {
        N = Math.min(N_MAX, smoothSize(T / dt));
        resp = edgeResponses(ctx, transfer, split.slowRes, N, dt, offsets);
        tail = tailOf(mainStep(resp.s), N);
        const cmOpen = diff && cmExcited(resp.s, N) && tailAmp(cmStep(resp.s), N) > CM_SETTLED;
        if ((tail <= 1e-6 && !cmOpen) || 2 * T > split.tCap || N >= N_MAX) break;
        const pre = precursor(mainStep(resp.s), N);
        nonCausal = pre > 2e-3 && pre > 0.7 * prevPre;
        if (nonCausal) break;
        prevPre = pre;
        T *= 2;
    }
    if (tail > 1e-6 && !nonCausal) {
        warnings.push('The step responses have not settled within their window (reflections still ringing); the eye is not accurate.');
    } else if (diff && cmExcited(resp.s, N) && tailAmp(cmStep(resp.s), N) > CM_SETTLED && !nonCausal) {
        warnings.push('The common mode the drive excites (skew, unequal edges or an asymmetric pair) has not settled within the window: the common-mode traces are not accurate, the differential eye is.');
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

    // Receiver and transmitter equalization: choose the TX preset and the CTLE gain
    // (automatic settings) by the peak-distortion eye of the pulse response.
    const eq = ctx.cfg.rxEq || {}, txEq = ctx.cfg.txEq || {};
    const dfeTaps = Math.max(0, Math.round(eq.dfeTaps || 0));
    // Tap magnitude limits of the DFE taps in order, V or { rel } of the cursor; missing ones
    // are unlimited.
    const dfeLimits = eq.dfeLimits || null;
    const ctleCands = !eq.ctle ? [null]
        : eq.ctle.auto ? Array.from({ length: Math.round(eq.ctle.adcMax - eq.ctle.adcMin) + 1 }, (_, i) => eq.ctle.adcMin + i)
        : [eq.ctle.adcDb];
    // txEq.presets: the presets to try (reduced swing implements only some).
    const firCands = txEq.auto
        ? Object.entries(PCIE_PRESETS).filter(([name]) => !txEq.presets || txEq.presets.includes(name))
            .map(([name, [ps, de]]) => ({ name, fir: firFromPresetDb(ps, de) }))
        : [{ name: txEq.name || null, fir: stim.fir || [0, 1, 0] }];
    let choice = null;
    const eqActive = !!eq.ctle || dfeTaps > 0 || txEq.auto;
    if (eqActive) {
        for (const adc of ctleCands) {
            const C = adc === null ? null : ctleOf(eq.ctle, adc);
            const steps = mainSteps(ctx, resp.main, C);
            for (const fc of firCands) {
                const q = pulseResponse(stim, steps, spb, fc.fir);
                const pda = pdaEye(q, spb, dfeTaps, dfeLimits);
                if (!choice || pda.height > choice.pda.height) choice = { adc, C, steps, fir: fc.fir, firName: fc.name, q, pda };
            }
        }
    }
    let st = stim;
    if (choice && txEq.auto) {
        st = buildStimulus({ ...ctx.cfg.stimulus, differential: diff, fir: choice.fir });
    }

    // Circular convolution over the period: edge trains (impulses on the bit boundaries)
    // with the folded step differences. With a CTLE the equalized main signal is one more
    // output ('eq').
    const useEqOut = !!(choice && choice.C);
    const nOut = no + (useEqOut ? 1 : 0);
    const respOf = (o, l, e) => o < no ? resp.s[iRes(o, l, e)] : choice.steps[l * 2 + e];
    const L = smoothSize(2 * M), K = L / 2;
    const acc = Array.from({ length: nOut }, () => new Float64Array(2 * (K + 1)));
    const buf = new Float64Array(L);
    for (let l = 0; l < nl; l++) {
        const leg = st.legs[l];
        // Bit boundary of each edge (the leg's skew is in its responses).
        const boundary = (e) => Math.round((e.t - (l === 1 ? skew : 0) - st.t0) / ui);
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
            for (let o = 0; o < nOut; o++) {
                const s = respOf(o, l, e);
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
    const conv = fft.irfftBatch(L, nOut, (() => {
        const all = new Float64Array(nOut * 2 * (K + 1));
        acc.forEach((a, o) => all.set(a, o * 2 * (K + 1)));
        return all;
    })());

    // Driver legs over the period: bit boundaries at k * UI from t = 0.
    const src = st.legs.map((leg, l) => sampleLegPeriodic(st, {
        v0: leg.v0, edges: leg.edges.map(e => ({ ...e, t: e.t - st.t0 })),
    }, M, dt));
    const xdc = dcSolve(net, ctx.Gdc, st.legs.map(leg => legMean(st, leg)));
    const o = {};
    // The equalized output: the CTLE's DC gain on the main signal's DC level.
    const adcLin = useEqOut ? choice.C(0)[0] : 1;
    const mainDc = diff ? xdc[out.nodes[0]] - xdc[out.nodes[1]] : xdc[out.nodes[0]];
    const outName = (k) => k < no ? out.names[k] : 'eq';
    const outBase = (k) => k < no ? xdc[out.nodes[k]] : adcLin * mainDc;
    for (let k = 0; k < nOut; k++) {
        // Wrap the linear convolution to the period, integrate, align the edges to
        // t = k * UI, remove the mean.
        const c = new Float64Array(M);
        for (let n = 0; n < L; n++) c[n % M] += conv[k * L + n];
        // G[m] = sum of c up to m is periodic (the edges of a period sum to zero); the
        // responses put edge k at index k * spb + nLead.
        const G = new Float64Array(M);
        let run = 0, mean = 0;
        for (let m = 0; m < M; m++) { run += c[m]; G[m] = run; mean += run / M; }
        const v = new Float64Array(M), base = outBase(k);
        for (let n = 0; n < M; n++) v[n] = base + G[(n + nLead) % M] - mean;
        o[outName(k)] = v;
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
            // The slow modes lie far below the CTLE corners: its DC gain applies.
            if (useEqOut) {
                const r = adcLin * (diff ? sm.r[l] - sm.r[nl + l] : sm.r[l]);
                const v = o.eq;
                for (let n = 0; n < M; n++) v[n] += r * (z[n] - zm);
            }
        }
    }

    const eqTrace = o.eq;
    delete o.eq;
    const full = deriveTraces(diff, o, src);
    const main = full[mainName(diff)];
    // Slicer signal: the CTLE output (or the main signal), less the ideal DFE feedback.
    let slicer = eqTrace || main, centre = null;
    if (choice) {
        // The main cursor of bit m sits at m UI + tau_s; the DFE corrects the bit's UI.
        const tauS = (choice.pda.n0 - nLead) * dt;
        centre = tauS;
        if (dfeTaps > 0) {
            slicer = Float64Array.from(slicer);
            const d = (m) => 2 * stim.bits[((m % nb) + nb) % nb] - 1;
            const h = Array.from({ length: dfeTaps }, (_, k) => {
                const i = choice.pda.n0 + (k + 1) * spb;
                return dfeTap(i < choice.q.length ? choice.q[i] : 0, k + 1, dfeLimits, choice.q[choice.pda.n0]);
            });
            choice.dfeLimited = h.map((x, k) => {
                const i = choice.pda.n0 + (k + 1) * spb;
                return i < choice.q.length && Math.abs(choice.q[i] - x) > 1e-15;
            });
            for (let m = 0; m < nb; m++) {
                let fb = 0;
                for (let k = 0; k < dfeTaps; k++) fb += h[k] * d(m - 1 - k);
                const n0 = Math.round((m * ui + tauS - ui / 2) / dt);
                for (let n = n0; n < n0 + spb; n++) slicer[((n % M) + M) % M] -= fb;
            }
            choice.dfe = h;
        }
        if (eqTrace || dfeTaps > 0) full.rxEq = slicer;
    }
    const eye = choice
        ? eyeAnalysis(slicer, dt, ui, spb, nb, { centre, bits: stim.bits, crossings: dfeTaps === 0 })
        : eyeAnalysis(main, dt, ui, spb, nb);
    const delayMetric = eyeDelay(eye, full[srcName(diff)], dt, ui, tauMin);
    const eqMetrics = [];
    if (choice) {
        const pins = eyeAnalysis(main, dt, ui, spb, nb);
        const ph = pins.metrics.find(m => m.label === 'Eye height (centre)');
        if (ph) eqMetrics.push({ label: 'Eye height at the pins (no RX EQ)', value: ph.value });
        if (txEq.auto || txEq.name || (stim.fir && stim.fir[1] !== 1)) {
            const db = presetDbFromFir(choice.fir);
            // Taps 0 / 1 / 0 (PCIe P4) send plain NRZ.
            const plain = Math.abs(choice.fir[0]) < 1e-9 && Math.abs(choice.fir[2]) < 1e-9;
            eqMetrics.push({ label: 'TX FIR' + (txEq.auto ? ' (chosen)' : ''), value: `${choice.firName ? choice.firName + ': ' : ''}`
                + (plain ? 'no TX equalization' : `${boostDb(choice.fir).toFixed(1)} dB boost`),
                note: (plain ? '' : `pre-shoot ${db.preshoot.toFixed(1)} dB, de-emphasis ${db.deemphasis.toFixed(1)} dB, `)
                    + `taps ${choice.fir.map(c => c.toFixed(3)).join(' / ')}`
                    + (choice.firName === 'P10' ? '; P10 at its nominal maximum boost, the real one follows the LF level the transmitter advertises' : '') });
        }
        if (choice.adc !== null) eqMetrics.push({ label: 'CTLE DC gain' + (eq.ctle.auto ? ' (chosen)' : ''), value: `${choice.adc} dB`,
            note: eq.ctle.order === 2 ? '2nd-order CTLE' : '1st-order CTLE' });
        if (choice.dfe) eqMetrics.push({ label: 'DFE taps', value: choice.dfe.map((x, k) => fmtV(x) + (choice.dfeLimited[k] ? ' (at its limit)' : '')).join(', '),
            note: 'ideal decisions; each tap the post-cursor at the sampling point, within its limit' });
        eqMetrics.push({ label: 'Pulse-response eye estimate', value: fmtV(choice.pda.height), note: 'peak distortion: main cursor minus the residual ISI, worst case' });
    }

    // The shown waveform: the first bits of the period.
    const nShow = Math.min(M, SHOW_BITS * spb);
    const t = new Float64Array(nShow);
    for (let n = 0; n < nShow; n++) t[n] = n * dt;
    const traces = Object.fromEntries(Object.entries(full).map(([k, v]) => [k, v.slice(0, nShow)]));
    const metrics = { list: [
        { label: 'Bit rate / UI', value: `${+(stim.bitRate / 1e9).toPrecision(4)} Gb/s / ${fmtTime(ui)}` },
        { label: 'Pattern length', value: `${nb} bits` },
        ...eye.metrics, delayMetric, ...eqMetrics,
    ].filter(Boolean), crossTime: NaN };
    const res = {
        kind: 'pattern', t, traces, metrics, hplot: resp.hplot,
        info: { dt, N, window: N * dt, slowModes: split.slow.length, trips, samplesPerBit: spb, periodSamples: M },
        eyeHist: eye.hist, eyeInfo: eye.info,
        eyeAt: choice && (eqTrace || dfeTaps > 0) ? 'slicer' : 'pins',
    };
    if (choice) {
        res.pulse = pulseOutput(choice, dt, spb, dfeTaps, nLead);
        res.eqChoice = { firName: choice.firName, fir: choice.fir, adc: choice.adc, pdaHeight: choice.pda.height, dfe: choice.dfe || [] };
    }
    packOutput(res);
    return res;
}

// Eye of one period of the main receiver trace: crossings at the threshold half way
// between the mean one and zero levels at the eye centre, the centre half a UI from the
// mean crossing phase. Returns { metrics, hist, info }.
// opts.centre (time of bit 0's sampling point) and opts.bits (transmitted bits) fix the
// centre and the decisions instead (an equalized slicer signal); the width is then the
// span of phases with an open eye. opts.crossings: also report the crossing jitter.
function eyeAnalysis(v, dt, ui, spb, nb, opts = null) {
    if (opts) return eyeAtCentre(v, dt, ui, nb, opts);
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

function eyeAtCentre(v, dt, ui, nb, { centre, bits, crossings }) {
    const M = v.length;
    const at = (t) => {
        const x = ((t / dt) % M + M) % M, n = Math.floor(x), f = x - n;
        return v[n] + (v[(n + 1) % M] - v[n]) * f;
    };
    const opening = (off) => {
        let min1 = Infinity, max0 = -Infinity, s1 = 0, s0 = 0, n1 = 0, n0 = 0;
        for (let m = 0; m < nb; m++) {
            const x = at(m * ui + centre + off);
            if (bits[m]) { min1 = Math.min(min1, x); s1 += x; n1++; } else { max0 = Math.max(max0, x); s0 += x; n0++; }
        }
        return { min1, max0, m1: s1 / n1, m0: s0 / n0, h: min1 - max0 };
    };
    const c = opening(0), thr = 0.5 * (c.m1 + c.m0), amp = c.m1 - c.m0;
    const metrics = [];
    let best = c.h, bestOff = 0, left = 0, right = 0;
    const steps = 64;
    for (let j = -steps / 2; j <= steps / 2; j++) {
        const off = j * ui / steps, h = opening(off).h;
        if (h > best + 1e-9 * Math.abs(amp)) { best = h; bestOff = off; }
    }
    // Open span around the centre (inner opening > 0), within +-UI/2.
    if (c.h > 0) {
        for (let j = 1; j <= steps / 2 && opening(-j * ui / steps).h > 0; j++) left = -j * ui / steps;
        for (let j = 1; j <= steps / 2 && opening(j * ui / steps).h > 0; j++) right = j * ui / steps;
    }
    metrics.push({ label: 'Eye height (centre)', value: c.h > 0 ? `${fmtV(c.h)} (${fmtPct(c.h / amp)})` : 'closed', note: 'percentage of the eye amplitude' });
    metrics.push({ label: 'Eye height (best phase)', value: best > 0 ? `${fmtV(best)} at ${bestOff >= 0 ? '+' : ''}${fmtTime(bestOff)}` : 'closed',
        note: 'largest inner opening over the sampling phase, relative to the centre' });
    metrics.push({ label: 'Eye width (open phases)', value: c.h > 0 ? `${fmtTime(right - left)} (${fmtPct((right - left) / ui)} UI)` : 'closed',
        note: 'span of sampling phases with an open eye, in steps of UI/64' });
    if (crossings) {
        const ph = [];
        for (let n = 0; n < M; n++) {
            const a = v[n] - thr, b = v[(n + 1) % M] - thr;
            if ((a < 0) !== (b < 0)) { let d = ((n + a / (a - b)) * dt) - centre - ui / 2; d -= ui * Math.round(d / ui); ph.push(d); }
        }
        if (ph.length) {
            const jpp = Math.max(...ph) - Math.min(...ph);
            const mean = ph.reduce((s, x) => s + x, 0) / ph.length;
            const jrms = Math.sqrt(ph.reduce((s, x) => s + (x - mean) ** 2, 0) / ph.length);
            metrics.push({ label: 'Crossing jitter pk-pk / rms', value: `${fmtTime(jpp)} / ${fmtTime(jrms)}`, note: 'data-dependent (ISI) jitter at the threshold' });
        }
    }
    metrics.push({ label: 'One / zero level (centre)', value: `${fmtV(c.m1)} / ${fmtV(c.m0)}` });
    metrics.push({ label: 'Eye amplitude', value: fmtV(amp) });
    const tc = ((centre % ui) + ui) % ui, start = tc - ui;
    const info = { thr, centre: 1, ui, crossPhase: ((tc - ui / 2) % ui + ui) % ui, centreTime: tc, start,
        heightTop: c.min1, heightBottom: c.max0, widthLeft: 1 + left / ui, widthRight: 1 + right / ui };
    return { metrics, hist: histogram(v, dt, ui, start), info };
}

// Pulse response for the plot: q around the main cursor, in UI from it, and its cursors.
function pulseOutput(choice, dt, spb, dfeTaps, nLead) {
    const q = choice.q, n0 = choice.pda.n0;
    const i0 = Math.max(0, n0 - 5 * spb), i1 = Math.min(q.length, n0 + 40 * spb);
    const stride = Math.max(1, Math.floor((i1 - i0) / 4000));
    const x = [], y = [];
    for (let i = i0; i < i1; i += stride) { x.push((i - n0) / spb); y.push(q[i]); }
    return { x: Float64Array.from(x), y: Float64Array.from(y), cursors: choice.pda.cursors, dfeTaps, dfe: choice.dfe || [],
        delay: (n0 - nLead) * dt };
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
