// Time-domain response of the link to a step, a pulse or a clock.
//
// Transient (step, pulse, clock burst): the response is
//   o(t) = DC operating point + slow part + fast part.
// The slow part is the AC coupling: modes of the lumped low-frequency circuit (line as
// its DC network) whose time constants exceed what an FFT window can hold, integrated
// exactly in time. The fast part is the rest of the exact transfer,
// H(f) - H_slow(f), times the stimulus spectrum, transformed with one inverse FFT.
// Its derivative is time limited (it settles), so the window holds it without wrap
// around, and it is integrated spectrally:
//   o(t_n) = O'(0) n/N + p_n - p_0,  p = IFFT[O'(f) / (j w)] (no DC bin).
//
// Clock steady state: the periodic response is exact on the harmonic grid, one inverse
// FFT over one period, any time constant included. With jitter, the clock runs through
// the pattern engine over many cycles (jitteredClock).
import { buildStimulus, legDerivSpectrum, legNetChange, legMean, sampleLeg, legValue, edgeSpectrum, requiredBandwidth } from './stimulus.js';
import { makeRLGCEvaluator, fmtHz } from './channel_grid.js';
import { buildNetwork, makeTransfer, dcMatrix, dcSolve, lumpedModes } from './network.js';
import { smoothSize } from './fft.js';
import { patternSim } from './pattern_sim.js';
import { jitterActive, jitterPerEdge, boundaryJitter, Q_BER12 } from './jitter.js';

const N_MAX = 1 << 21;
const PLOT_BUCKETS = 6000;
const CSV_MAX_ROWS = 100000;
// Cycles of a jittered steady clock.
const CLOCK_JITTER_CYCLES = 4096;

// Delay, impedance and attenuation of the line modes at f (estimates for sizing).
function lineEstimates(ev, f, len, diff) {
    const v = ev(f);
    const w = 2 * Math.PI * f;
    const mode = (R, L, G, C) => {
        const Z = Math.sqrt(L / C);
        return { Z, tau: len * Math.sqrt(L * C), alpha: R / (2 * Z) + G * Z / 2, R, L, G, C, w };
    };
    if (!diff) return { modes: [mode(v.R, v.L, v.G, v.C)] };
    const avg = (k) => 0.5 * (v[k + '11'] + v[k + '22']);
    const odd = mode(avg('R') - v.R12, avg('L') - v.L12, avg('G') - v.G12, avg('C') - v.C12);
    const even = mode(avg('R') + v.R12, avg('L') + v.L12, avg('G') + v.G12, avg('C') + v.C12);
    return { modes: [odd, even] };
}

// Termination impedance a mode sees at high frequency (coupling caps shorted).
function termImpedance(term, modeIdx) {
    switch (term.type) {
        case 'vtt': case 'r': case 'rc': return term.rt;
        case 'thev': return term.rup * term.rdn / (term.rup + term.rdn);
        case 'diff': return modeIdx === 0 ? term.rdiff / 2 : Infinity;
        default: return Infinity;
    }
}

// Round trips until the reflections of the first mode decay below 1e-5. est holds the
// line at a low frequency: the low-frequency content of an edge settles last, with
// the least attenuation per trip.
function reflectionTrips(est, net) {
    const m = est.modes[0];
    const gam = (a, b) => isFinite(a) ? Math.abs((a - b) / (a + b)) : 1;
    const zl = termImpedance(net.term, 0);
    // Each end reflects at its worst junction: the line, a line-form package, the source
    // or the termination (behind an R + L package's resistance).
    const end = (pkg, z) => !pkg ? gam(z, m.Z)
        : pkg.form === 'line' ? Math.max(gam(z, m.Z), gam(pkg.z0, m.Z), gam(z, pkg.z0)) : gam(z + pkg.r, m.Z);
    const pk = (near) => (net.packages || []).find(p => p.near === near);
    const tx = pk(net.nodes.a), rx = net.nodes.q ? pk(net.nodes.q) : null;
    const gs = end(tx && tx.model, net.rs), gl = end(rx && rx.model, zl);
    const rho = Math.abs(gs * gl) * Math.exp(-2 * m.alpha * m.tau / Math.sqrt(m.L * m.C));
    if (rho < 1e-12) return 0;
    if (rho >= 0.999) return 200;
    return Math.min(200, Math.ceil(Math.log(1e-5) / Math.log(rho)));
}

// Raised-cosine taper from half the Nyquist frequency to zero at it.
function taper(f, fNyq) {
    const f0 = 0.5 * fNyq;
    if (f <= f0) return 1;
    if (f >= fNyq) return 0;
    return 0.5 * (1 + Math.cos(Math.PI * (f - f0) / (fNyq - f0)));
}

// Output node voltages and names.
function outputs(net) {
    const { a, d } = net.nodes;
    return net.differential
        ? { nodes: [d[0], d[1], a[0], a[1]], names: ['rxP', 'rxN', 'txP', 'txN'] }
        : { nodes: [d[0], a[0]], names: ['rx', 'tx'] };
}

// Derived traces: differential and common mode of a pair, and the source EMF.
function deriveTraces(diff, o, src) {
    if (!diff) return { rx: o.rx, tx: o.tx, src: src[0] };
    const n = o.rxP.length;
    const mk = (f) => { const r = new Float64Array(n); for (let i = 0; i < n; i++) r[i] = f(i); return r; };
    return {
        rxDiff: mk(i => o.rxP[i] - o.rxN[i]),
        rxCM: mk(i => 0.5 * (o.rxP[i] + o.rxN[i])),
        rxP: o.rxP, rxN: o.rxN,
        txDiff: mk(i => o.txP[i] - o.txN[i]),
        txP: o.txP, txN: o.txN,
        srcDiff: mk(i => src[0][i] - src[1][i]),
        srcP: src[0], srcN: src[1],
    };
}

const mainName = (diff) => diff ? 'rxDiff' : 'rx';
const srcName = (diff) => diff ? 'srcDiff' : 'src';

/**
 * Simulate the link.
 * table: channel table (channel_grid.js). fft: loadFFT() instance.
 * cfg: { length, network, stimulus, samplesPerEdge, span, causal }
 */
function simulateLink(table, cfg, fft) {
    const len = cfg.length;
    if (!(len > 0)) throw new Error('Line length must be positive.');
    const diff = table.differential;
    const spe = Math.min(256, Math.max(4, cfg.samplesPerEdge || 16));
    const dt0 = Math.min(cfg.stimulus.tr, cfg.stimulus.tf) / spe;
    // A quiet lead before the first edge: the spectral integration takes its constant
    // from the first sample. The time axis is shifted back by it.
    const steady = (cfg.stimulus.type === 'clock' && cfg.stimulus.steady) || cfg.stimulus.type === 'pattern';
    const lead = steady ? 0 : 64 * dt0;
    const stim = buildStimulus({ ...cfg.stimulus, differential: diff, delay: (cfg.stimulus.delay || 0) + lead });
    // cfg.perturb(v, f): a change of the line parameters (tests of sensitivity.js).
    const ev0 = makeRLGCEvaluator(table);
    const ev = cfg.perturb ? (f) => cfg.perturb(ev0(f), f) : ev0;
    const net = buildNetwork({ ...cfg.network, differential: diff });
    const fTop = table.freqs[table.freqs.length - 1];
    const warnings = [];
    const Gdc = dcMatrix(net, ev(0), len);
    const modes = lumpedModes(net, Gdc);
    const out = outputs(net);
    const est = lineEstimates(ev, 0.35 / stim.t1090Min, len, diff);
    const ctx = { table, cfg, fft, len, diff, stim, ev, net, fTop, warnings, Gdc, modes, out, est, lead };
    if (net.rxFloating) {
        warnings.push(diff
            ? 'The receiver has no DC path: its common-mode level is set by the receiver bias network, not modelled here (a 1 GOhm leak holds it at 0 V).'
            : 'The receiver has no DC path: its DC level is set by the receiver bias network, not modelled here (a 1 GOhm leak holds it at 0 V).');
    }
    const res = stim.type === 'pattern' ? patternSim(ctx, dt0)
        : stim.steady ? (jitterActive(cfg.stimulus.jitter) ? jitteredClock(ctx, dt0) : steadyState(ctx, dt0))
        : transient(ctx, dt0);
    res.differential = diff;
    // A current-mode driver's transfer is a transimpedance (V/A).
    res.currentDrive = net.current;
    res.stimType = stim.type;
    res.period = stim.period;
    delete res.eyeSignal;
    res.warnings = [...new Set(warnings)];
    res.tau = est.modes.map(m => m.tau);
    res.fTop = fTop;
    res.fNeed = requiredBandwidth(cfg.stimulus);
    return res;
}

// Exact first-order response of a slow mode to a sampled input s on a uniform grid
// (input linear between samples), from rest.
function firstOrder(lambda, s, dt) {
    const n = s.length, z = new Float64Array(n);
    const a = Math.exp(-dt / lambda), k = (lambda / dt) * (1 - a);
    for (let i = 1; i < n; i++) {
        const u0 = s[i - 1], u1 = s[i];
        z[i] = a * z[i - 1] + u1 - a * u0 - (u1 - u0) * k;
    }
    return z;
}

function transient(ctx, dt0) {
    const { stim, net, est, modes, warnings, ev, len, out, fTop, diff } = ctx;
    const nl = net.nl, no = out.nodes.length;
    const tauLine = Math.max(...est.modes.map(m => m.tau));
    const trips = reflectionTrips(lineEstimates(ev, 1 / (8 * tauLine), len, diff), net);
    // Delays of the path: the line modes plus the packages.
    const tauMax = tauLine + net.pkgDelay, tauMin = Math.min(...est.modes.map(m => m.tau)) + net.pkgDelay;
    const edgeLong = Math.max(stim.up.t1090, stim.down.t1090);
    const tBase = stim.tEnd + tauMax * (1 + 2 * trips) + 20 * edgeLong + 2 * stim.halfMax;

    const split = slowSplit(ctx, tBase, dt0);
    const { tCap, slow, slowRes } = split;
    let T = split.T;

    // Receiver response before the first possible arrival, relative to its swing.
    let tEdge0 = Infinity;
    for (const leg of stim.legs) for (const e of leg.edges) tEdge0 = Math.min(tEdge0, e.t - e.edge.half);
    const tArrive = tEdge0 + 0.9 * tauMin;
    const precursor = (v, n, dtv) => {
        let pre = 0, swing = 0;
        const iPre = Math.min(n, Math.floor(tArrive / dtv));
        for (let i = 0; i < n; i++) swing = Math.max(swing, Math.abs(v(i) - v(0)));
        for (let i = 0; i < iPre; i++) pre = Math.max(pre, Math.abs(v(i) - v(0)));
        return swing > 0 ? pre / swing : 0;
    };

    const transfer = makeTransfer(net, len, out.nodes);
    let fast = null, N = 0, dt = dt0, nonCausal = false, prevPre = Infinity;
    for (let attempt = 0; attempt < 5; attempt++) {
        N = smoothSize(T / dt0);
        dt = dt0;
        if (N > N_MAX) { N = N_MAX; dt = T / N; }
        fast = fastResponse(ctx, transfer, slowRes, N, dt);
        if (fast.tail <= 1e-6 || 2 * T > tCap) break;
        // Unsettled ringing wraps around to the start of the window and shrinks as the
        // window grows. A non-causal line model responds before t = 0, which wraps into
        // any window: when a longer window leaves the early response as it was, stop.
        const fo = fast.o;
        const pre = precursor(diff ? (i) => fo[0][i] - fo[1][i] : (i) => fo[0][i], N, dt);
        nonCausal = pre > 2e-3 && pre > 0.7 * prevPre;
        if (nonCausal) break;
        prevPre = pre;
        T *= 2;
    }
    if (dt > dt0) warnings.push(`The window needs more than ${N_MAX} samples: the time step is coarsened to ${fmtTime(dt)}.`);
    if (fast.tail > 1e-6 && !nonCausal) {
        warnings.push('The response has not settled within the simulation window (reflections still ringing); the end of the window is not accurate.');
    }
    T = N * dt;

    // Display grid: the FFT window, extended or cropped to the requested span.
    // Automatic span: the edges and the settling estimate (the window can be longer, it
    // holds coupling time constants and the extra settling found by the tail check).
    const span = cfg_span(ctx, Math.min(T, tBase) - ctx.lead) + ctx.lead;
    const nFine = Math.min(N, Math.max(2, Math.floor(span / dt) + 1));
    const nExt = span > T * (1 + 1e-9) ? 4000 : 0;
    const t = new Float64Array(nFine + nExt);
    for (let i = 0; i < nFine; i++) t[i] = i * dt;
    for (let m = 0; m < nExt; m++) t[nFine + m] = T + (m + 1) * (span - T) / nExt;

    // Driver EMF of each leg on the display grid.
    const src = stim.legs.map(leg => {
        const s = new Float64Array(nFine + nExt);
        s.set(sampleLeg(stim, leg, dt, nFine));
        const fin = leg.v0 + legNetChange(stim, leg);
        for (let m = 0; m < nExt; m++) s[nFine + m] = fin;
        return s;
    });
    const du = src.map((s, l) => {
        const v0 = stim.legs[l].v0, r = new Float64Array(s.length);
        for (let i = 0; i < s.length; i++) r[i] = s[i] - v0;
        return r;
    });

    const xdc = dcSolve(net, ctx.Gdc, stim.legs.map(l => l.v0));
    const o = {};
    for (let k = 0; k < no; k++) {
        const v = new Float64Array(nFine + nExt);
        const f = fast.o[k];
        for (let i = 0; i < nFine; i++) v[i] = xdc[out.nodes[k]] + f[i];
        for (let m = 0; m < nExt; m++) v[nFine + m] = xdc[out.nodes[k]] + fast.final[k];
        o[out.names[k]] = v;
    }
    // Slow modes: exact on the fine grid, then the exponential relaxation to the final
    // input over the extension.
    for (const sm of slowRes) {
        // The residue r carries P * beta, so each leg drives the mode with its own du.
        for (let l = 0; l < nl; l++) {
            const z = firstOrder(sm.lambda, du[l].subarray(0, nFine), dt);
            const zT = z[nFine - 1], uf = du[l][du[l].length - 1];
            for (let k = 0; k < no; k++) {
                const r = sm.r[k * nl + l];
                if (r === 0) continue;
                const v = o[out.names[k]];
                for (let i = 0; i < nFine; i++) v[i] += r * z[i];
                for (let m = 0; m < nExt; m++) {
                    const tt = t[nFine + m] - t[nFine - 1];
                    v[nFine + m] += r * (uf + (zT - uf) * Math.exp(-tt / sm.lambda));
                }
            }
        }
    }

    const traces = deriveTraces(diff, o, src);
    const main = traces[mainName(diff)];

    const preRel = precursor((i) => main[i], main.length, dt);
    if (preRel > 2e-3) {
        warnings.push(`The receiver responds ${(100 * preRel).toFixed(2)} % of the swing before the signal can arrive`
            + (ctx.cfg.causal ? '.' : ': enable causal materials for a causal dielectric model.'));
    }
    if (fast.extrapFraction > 0.01) {
        warnings.push(`${(100 * fast.extrapFraction).toFixed(1)} % of the stimulus spectrum is above the sweep stop `
            + `${fmtHz(fTop)}, where the line parameters are extrapolated.`);
    }

    // Step size of the first data edge at the main output (fast part settled).
    const first = stim.legs.map(l => l.edges[0].dv);
    const hf0 = fast.hf0;
    const stepOf = (name) => {
        const idx = out.names.indexOf(name);
        let s = 0;
        for (let l = 0; l < nl; l++) s += hf0[idx * nl + l] * first[l];
        return s;
    };
    const step = diff ? stepOf('rxP') - stepOf('rxN') : stepOf('rx');

    const metrics = transientMetrics(t, traces, diff, stim, step, slow.length > 0, nFine);
    for (let i = 0; i < t.length; i++) t[i] -= ctx.lead;
    if (isFinite(metrics.crossTime)) metrics.crossTime -= ctx.lead;
    const res = {
        kind: 'transient', t, traces, metrics, hplot: fast.hplot,
        info: { dt, N, window: T, slowModes: slow.length, trips },
        markers: [{ t: stim.t0 - ctx.lead, label: 'source 50 %' }, { t: stim.t0 - ctx.lead + tauMin, label: net.pkgDelay > 0 ? '+ line and package delay' : '+ line delay' }],
    };
    if (stim.type === 'clock') res.eye = foldEye(t, main, stim.period, metrics.crossTime, stim.count, nFine);
    packOutput(res);
    return res;
}

// Lumped modes too slow for the window, handled in the time domain, and the window that
// holds the others. noSlowSplit (tests): every mode below 1 ms in the window, the
// brute-force reference (the 1 GOhm leaks of floating nodes stay slow).
// Returns { T, tCap, slow, slowRes }: slowRes the residues at the outputs,
// r[o][l] / (1 + j w lambda).
function slowSplit(ctx, tBase, dt0) {
    const { modes, net, out } = ctx;
    const nl = net.nl, no = out.nodes.length;
    const tCap = ctx.cfg.noSlowSplit ? 5e-3 : Math.min(64 * tBase, N_MAX * dt0);
    let T = tBase;
    const slow = [];
    const order = modes.lambda.map((l, i) => i).filter(i => modes.lambda[i] > 0)
        .sort((i, j) => modes.lambda[i] - modes.lambda[j]);
    for (const i of order) {
        const lam = modes.lambda[i];
        if (5 * lam <= T) continue;
        if (5 * lam <= tCap) { T = 5 * lam; continue; }
        slow.push(i);
    }
    const slowRes = slow.map(i => {
        const r = new Float64Array(no * nl);
        for (let o = 0; o < no; o++) for (let l = 0; l < nl; l++) {
            r[o * nl + l] = modes.P[out.nodes[o] * net.n + i] * modes.beta[i * nl + l];
        }
        return { lambda: modes.lambda[i], r };
    });
    return { T, tCap, slow, slowRes };
}

function cfg_span(ctx, T) {
    const s = ctx.cfg.span;
    return s > 0 ? s : T;
}

// The fast part on N samples: spectra, inverse FFT, spectral integration.
function fastResponse(ctx, transfer, slowRes, N, dt) {
    const { stim, net, ev, out, fft, fTop, diff } = ctx;
    const nl = net.nl, no = out.nodes.length;
    const K = N / 2, df = 1 / (N * dt), fNyq = K * df;
    const spec = new Float64Array(no * (K + 1) * 2);
    const hf0 = new Float64Array(no * nl);
    const scratch = { re: new Float64Array(no * nl), im: new Float64Array(no * nl) };
    const plotStride = Math.max(1, Math.floor(K / 2000));
    const hf = [], hm = [], hp = [], sm = [];
    let eTot = 0, eExt = 0;
    const net0 = stim.legs.map(leg => legNetChange(stim, leg));
    const O0 = new Float64Array(no);
    for (let k = 0; k <= K; k++) {
        const f = k * df, w = 2 * Math.PI * f;
        const H = transfer(f, ev(f), scratch);
        // Main transfer for the plot, before the slow part is taken out.
        if (k % plotStride === 0 && k > 0) {
            let re, im;
            if (!diff) { re = H.re[0]; im = H.im[0]; } else {
                // Differential out per differential in: (H_PP - H_PN - H_NP + H_NN) / 2.
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
            for (let i = 0; i < no * nl; i++) hf0[i] = H.re[i];
            for (let o = 0; o < no; o++) {
                let s = 0;
                for (let l = 0; l < nl; l++) s += H.re[o * nl + l] * net0[l];
                O0[o] = s;
            }
            continue;
        }
        if (k === K) continue;
        const W = taper(f, fNyq);
        if (W === 0) continue;
        const U = stim.legs.map(leg => legDerivSpectrum(stim, leg, f));
        const e = (U[0][0] ** 2 + U[0][1] ** 2) * W * W;
        eTot += e;
        if (f > fTop) eExt += e;
        for (let o = 0; o < no; o++) {
            let sr = 0, si = 0;
            for (let l = 0; l < nl; l++) {
                const hr = H.re[o * nl + l], hi = H.im[o * nl + l];
                sr += hr * U[l][0] - hi * U[l][1];
                si += hr * U[l][1] + hi * U[l][0];
            }
            // O'(f) W / (j w dt)
            const sc = W / (w * dt);
            const base = 2 * (o * (K + 1) + k);
            spec[base] = si * sc;
            spec[base + 1] = -sr * sc;
        }
    }
    const p = fft.irfftBatch(N, no, spec);
    const o = [], final = new Float64Array(no);
    let eTail = 0, eAll = 0;
    for (let k = 0; k < no; k++) {
        const v = new Float64Array(N);
        const p0 = p[k * N];
        for (let n = 0; n < N; n++) v[n] = O0[k] * n / N + p[k * N + n] - p0;
        for (let n = 1; n < N; n++) {
            const d = (v[n] - v[n - 1]) ** 2;
            eAll += d;
            if (n > 0.9 * N) eTail += d;
        }
        o.push(v);
        final[k] = O0[k];
    }
    return {
        o, final, hf0, tail: eAll > 0 ? eTail / eAll : 0,
        extrapFraction: eTot > 0 ? eExt / eTot : 0,
        hplot: { f: Float64Array.from(hf), mag: Float64Array.from(hm), phase: Float64Array.from(hp), stim: Float64Array.from(sm) },
    };
}

function steadyState(ctx, dt0) {
    const { stim, net, ev, len, out, fft, diff, warnings, fTop } = ctx;
    const nl = net.nl, no = out.nodes.length;
    const P = stim.period;
    let M = smoothSize(Math.max(64, P / dt0));
    if (M > N_MAX) {
        M = N_MAX;
        warnings.push(`One clock period needs more than ${N_MAX} samples: the time step is coarsened.`);
    }
    const dt = P / M, K = M / 2, fNyq = K / P;
    const transfer = makeTransfer(net, len, out.nodes);
    const scratch = { re: new Float64Array(no * nl), im: new Float64Array(no * nl) };
    const spec = new Float64Array(no * (K + 1) * 2);
    const plotStride = Math.max(1, Math.floor(K / 2000));
    const hf = [], hm = [], hp = [], sm = [];
    let eTot = 0, eExt = 0;
    for (let k = 1; k < K; k++) {
        const f = k / P, w = 2 * Math.PI * f;
        const H = transfer(f, ev(f), scratch);
        const W = taper(f, fNyq);
        // Fourier coefficients of the legs: (1/P) U'_period(f) / (j w).
        const U = stim.legs.map(leg => {
            const [r, i] = legDerivSpectrum(stim, leg, f, true);
            const s = W / (w * P);
            return [i * s, -r * s];
        });
        const e = U[0][0] ** 2 + U[0][1] ** 2;
        eTot += e; if (f > fTop) eExt += e;
        let hr0, hi0;
        if (!diff) { hr0 = H.re[0]; hi0 = H.im[0]; } else {
            hr0 = 0.5 * (H.re[0] - H.re[1] - H.re[2] + H.re[3]);
            hi0 = 0.5 * (H.im[0] - H.im[1] - H.im[2] + H.im[3]);
        }
        if (k % plotStride === 0) {
            hf.push(f); hm.push(Math.hypot(hr0, hi0)); hp.push(Math.atan2(hi0, hr0));
            sm.push(Math.abs(edgeSpectrum(stim.up, f)));
        }
        for (let o = 0; o < no; o++) {
            let sr = 0, si = 0;
            for (let l = 0; l < nl; l++) {
                const hr = H.re[o * nl + l], hi = H.im[o * nl + l];
                sr += hr * U[l][0] - hi * U[l][1];
                si += hr * U[l][1] + hi * U[l][0];
            }
            const base = 2 * (o * (K + 1) + k);
            spec[base] = sr * M; spec[base + 1] = si * M;
        }
    }
    if (eTot > 0 && eExt / eTot > 0.01) {
        warnings.push(`${(100 * eExt / eTot).toFixed(1)} % of the clock spectrum is above the sweep stop `
            + `${fmtHz(fTop)}, where the line parameters are extrapolated.`);
    }
    const xbar = dcSolve(net, ctx.Gdc, stim.legs.map(leg => legMean(stim, leg)));
    const p = fft.irfftBatch(M, no, spec);
    // Two periods for display.
    const n2 = 2 * M;
    const t = new Float64Array(n2);
    for (let i = 0; i < n2; i++) t[i] = i * dt;
    const o = {};
    for (let k = 0; k < no; k++) {
        const v = new Float64Array(n2), base = xbar[out.nodes[k]];
        for (let i = 0; i < n2; i++) v[i] = base + p[k * M + (i % M)];
        o[out.names[k]] = v;
    }
    const src = stim.legs.map(leg => {
        const s = new Float64Array(n2);
        for (let i = 0; i < n2; i++) s[i] = legValue(stim, leg, t[i]);
        return s;
    });
    const traces = deriveTraces(diff, o, src);
    const main = traces[mainName(diff)];
    const metrics = clockMetrics(t.subarray(0, M), main.subarray(0, M), traces[srcName(diff)].subarray(0, M), P, ctx);
    const res = {
        kind: 'steady', t, traces, metrics,
        hplot: { f: Float64Array.from(hf), mag: Float64Array.from(hm), phase: Float64Array.from(hp), stim: Float64Array.from(sm) },
        info: { dt, N: M, window: P, slowModes: 0 },
    };
    res.eye = foldSteady(main.subarray(0, M), P, metrics.crossTime);
    packOutput(res);
    return res;
}

// ---------------------------------------------------------------- metrics

// First time in [i0, i1) where v crosses level in direction dir (+1 up, -1 down).
function crossing(t, v, level, i0, i1, dir) {
    for (let i = Math.max(1, i0); i < i1; i++) {
        const a = dir * (v[i - 1] - level), b = dir * (v[i] - level);
        if (a < 0 && b >= 0) return t[i - 1] + (t[i] - t[i - 1]) * (a / (a - b));
    }
    return NaN;
}

// Last time before index i where v crosses level in direction dir.
function crossingBefore(t, v, level, i, dir, iMin = 0) {
    for (let j = Math.min(i, v.length - 1); j > Math.max(0, iMin); j--) {
        const a = dir * (v[j - 1] - level), b = dir * (v[j] - level);
        if (a < 0 && b >= 0) return t[j - 1] + (t[j] - t[j - 1]) * (a / (a - b));
    }
    return NaN;
}

const indexAt = (t, x) => {
    let lo = 0, hi = t.length - 1;
    if (!(x > t[0])) return 0;
    if (x >= t[hi]) return hi;
    while (hi - lo > 1) { const m = (lo + hi) >> 1; if (t[m] <= x) lo = m; else hi = m; }
    return hi;
};

// Rise or fall times and the 50 % time of the edge from lo to hi around index range.
function edgeTimes(t, v, lo, hi, i0, i1) {
    const dir = hi > lo ? 1 : -1, lv = (p) => lo + p * (hi - lo);
    const t50 = crossing(t, v, lv(0.5), i0, i1, dir);
    if (!isFinite(t50)) return { t50 };
    const i50 = indexAt(t, t50);
    const t10 = crossingBefore(t, v, lv(0.1), i50, dir, i0), t20 = crossingBefore(t, v, lv(0.2), i50, dir, i0);
    const t90 = crossing(t, v, lv(0.9), i50 - 1, i1, dir), t80 = crossing(t, v, lv(0.8), i50 - 1, i1, dir);
    return { t50, i50, t1090: t90 - t10, t2080: t80 - t20 };
}

function transientMetrics(t, traces, diff, stim, step, acCoupled, nFine) {
    const v = traces[mainName(diff)], s = traces[srcName(diff)];
    const n = v.length;
    const m = [];
    const lo = v[0], hi = lo + step;
    const out = { list: m, crossTime: NaN };
    const sLo = s[0], sHi = s[0] + (diff ? 2 : 1) * stim.legs[0].edges[0].dv;
    const src = edgeTimes(t, s, sLo, sHi, 0, n);
    if (!(Math.abs(step) > 1e-12 * Math.max(1, Math.abs(lo)))) {
        m.push({ label: 'Receiver step', value: 'none (the coupling blocks it)' });
        return out;
    }
    const rise = edgeTimes(t, v, lo, hi, 0, n);
    const dirWord = step > 0 ? 'Rise' : 'Fall';
    m.push({ label: 'Initial level', value: fmtV(lo) });
    m.push({ label: 'Step level', value: fmtV(hi), note: 'settled fast response to the first edge' });
    if (!isFinite(rise.t50)) {
        m.push({ label: 'Delay (50 %)', value: 'the receiver never crosses 50 %' });
        return out;
    }
    out.crossTime = rise.t50;
    m.push({ label: 'Delay (50 %)', value: fmtTime(rise.t50 - src.t50) });
    m.push({ label: `${dirWord} time 10-90 %`, value: fmtTime(rise.t1090) });
    m.push({ label: `${dirWord} time 20-80 %`, value: fmtTime(rise.t2080) });
    const dir = Math.sign(step);
    if (stim.type === 'step') {
        let peak = -Infinity, iPeak = rise.i50;
        for (let i = rise.i50; i < n; i++) { const x = dir * (v[i] - hi); if (x > peak) { peak = x; iPeak = i; } }
        let ring = Infinity;
        for (let i = iPeak; i < n; i++) ring = Math.min(ring, dir * (v[i] - hi));
        m.push({ label: 'Overshoot', value: fmtPct(Math.max(0, peak) / Math.abs(step)) });
        m.push({ label: 'Ringback', value: fmtPct(Math.max(0, -ring) / Math.abs(step)) });
        let iLast = -1;
        for (let i = rise.i50; i < nFine; i++) if (Math.abs(v[i] - hi) > 0.02 * Math.abs(step)) iLast = i;
        m.push({ label: 'Settling time (2 %)', value: iLast >= nFine - 1 ? 'not settled' : fmtTime(t[Math.max(iLast + 1, rise.i50)] - rise.t50) });
        if (acCoupled) m.push({ label: 'Level at end, relative to step', value: fmtPct((v[n - 1] - lo) / step) });
        return out;
    }
    // Pulse and clock: the first falling edge after the rise.
    const fall = edgeTimes(t, v, hi, lo, rise.i50 + 1, n);
    if (stim.type === 'pulse') {
        let peak = -Infinity;
        const iEnd = isFinite(fall.t50) ? fall.i50 : n;
        for (let i = rise.i50; i < iEnd; i++) peak = Math.max(peak, dir * (v[i] - lo));
        m.push({ label: 'Pulse peak', value: `${fmtV(lo + dir * peak)} (${fmtPct(peak / Math.abs(step))} of step)` });
        if (isFinite(fall.t50)) {
            m.push({ label: 'Width (50 %)', value: fmtTime(fall.t50 - rise.t50) });
            m.push({ label: `${step > 0 ? 'Fall' : 'Rise'} time 10-90 %`, value: fmtTime(fall.t1090) });
            m.push({ label: `${step > 0 ? 'Fall' : 'Rise'} time 20-80 %`, value: fmtTime(fall.t2080) });
            let und = Infinity;
            for (let i = fall.i50; i < n; i++) und = Math.min(und, dir * (v[i] - lo));
            m.push({ label: 'Undershoot after pulse', value: fmtPct(Math.max(0, -und) / Math.abs(step)) });
        } else {
            m.push({ label: 'Width (50 %)', value: 'does not return below 50 %' });
        }
        m.push({ label: 'Overshoot', value: fmtPct(Math.max(0, peak - Math.abs(step)) / Math.abs(step)) });
        return out;
    }
    // Clock burst: swing of the last period.
    const P = stim.period;
    const tLast = rise.t50 + (stim.count - 1) * P;
    const i0 = indexAt(t, tLast - 0.25 * P), i1 = indexAt(t, tLast + 0.75 * P);
    let vmax = -Infinity, vmin = Infinity;
    for (let i = i0; i < i1; i++) { vmax = Math.max(vmax, v[i]); vmin = Math.min(vmin, v[i]); }
    if (i1 > i0) m.push({ label: 'Swing in the last period', value: fmtV(vmax - vmin) });
    if (isFinite(fall.t50)) m.push({ label: 'First high time (50 %)', value: fmtTime(fall.t50 - rise.t50) });
    return out;
}

// Steady clock with transmitter jitter: the clock as a periodic pattern of alternating
// bits (UI = P / 2) over CLOCK_JITTER_CYCLES cycles in the pattern engine, each edge with
// its own jitter draw and the falling ones moved by the duty cycle. The waveform and eye
// are of that run; the clock metrics are of the clock without jitter.
function jitteredClock(ctx, dt0) {
    const nominal = steadyState(ctx, dt0);
    const c = ctx.cfg.stimulus, P = 1 / c.clockFreq, n = 2 * CLOCK_JITTER_CYCLES;
    const pe = jitterPerEdge(c.jitter), duty = c.duty ?? 0.5;
    if (pe.dj / 2 + Q_BER12 * pe.rj >= 0.5 * Math.min(duty, 1 - duty) * P) {
        throw new Error('The transmitter jitter at BER 1e-12 exceeds half of the high or low time.');
    }
    const bits = new Uint8Array(n);
    for (let k = 0; k < n; k += 2) bits[k] = 1;
    const rand = boundaryJitter(n, c.jitter), times = Float64Array.from(rand);
    for (let k = 1; k < n; k += 2) times[k] += (duty - 0.5) * P;
    const pcfg = { ...c, type: 'pattern', bitRate: 2 * c.clockFreq, bits, fir: null, jitter: { times } };
    const stim = buildStimulus({ ...pcfg, differential: ctx.diff });
    const run = patternSim({ ...ctx, stim, cfg: { ...ctx.cfg, stimulus: pcfg, rxEq: null, txEq: null, statEye: null } }, dt0);
    let peak = 0;
    for (const x of rand) peak = Math.max(peak, Math.abs(x));
    // Rising crossings of the receiver at the eye threshold: time interval error against
    // the ideal clock and period jitter.
    const v = run.eyeSignal, M = v.length, dt = run.info.dt, thr = run.eyeInfo.thr;
    const rise = [];
    for (let i = 0; i < M; i++) {
        const a = v[i] - thr, b = v[(i + 1) % M] - thr;
        if (a < 0 && b >= 0) rise.push((i + a / (a - b)) * dt);
    }
    const jitterRows = [];
    if (rise.length > 2) {
        const tie = rise.map((t, k) => t - k * P), mean = tie.reduce((x, y) => x + y, 0) / tie.length;
        const per = rise.map((t, k) => (k + 1 < rise.length ? rise[k + 1] : rise[0] + M * dt) - t - P);
        const ppRms = (a, m = 0) => [Math.max(...a) - Math.min(...a), Math.sqrt(a.reduce((x, y) => x + (y - m) ** 2, 0) / a.length)];
        const [tpp, trms] = ppRms(tie, mean), [ppp, prms] = ppRms(per);
        jitterRows.push({ label: 'Rising-edge TIE pk-pk / rms', value: `${fmtTime(tpp)} / ${fmtTime(trms)}`,
            note: 'time interval error of the receiver rising crossings against an ideal clock, this run' });
        jitterRows.push({ label: 'Period jitter pk-pk / rms', value: `${fmtTime(ppp)} / ${fmtTime(prms)}`,
            note: 'rising crossing to the next, against the period, this run' });
    }
    delete run.eyeSignal;
    const pick = (label) => run.metrics.list.find(m => m.label === label);
    const list = nominal.metrics.list.filter(m => m.label !== 'Eye height (centre)').map(m => ({ ...m,
        note: (m.note ? m.note + '; ' : '') + 'without jitter' }));
    list.push(...['Eye height (centre)', 'Eye width (at threshold)', 'Crossing jitter pk-pk / rms'].map(pick).filter(Boolean), ...jitterRows);
    list.push({ label: 'TX jitter per edge', value: `${fmtTime(pe.rj)} rms + ${fmtTime(pe.dj)} pp`,
        note: `random (Rj with the PW Rj) plus dual Dirac (DjDD with the PW DDj); ${fmtTime(pe.tj)} total at BER 1e-12` });
    list.push({ label: 'TX jitter of this run', value: `${fmtTime(peak)} largest displacement`,
        note: `one random draw per edge over ${CLOCK_JITTER_CYCLES} cycles; the eye shows this draw, not the BER 1e-12 tails` });
    return { ...run, kind: 'steady', metrics: { list, crossTime: nominal.metrics.crossTime }, hplot: nominal.hplot,
        info: { ...run.info, jitterCycles: CLOCK_JITTER_CYCLES } };
}

// Metrics of one steady-state clock period.
function clockMetrics(t, v, s, P, ctx) {
    const n = v.length;
    let vmax = -Infinity, vmin = Infinity;
    for (let i = 0; i < n; i++) { vmax = Math.max(vmax, v[i]); vmin = Math.min(vmin, v[i]); }
    const m = [];
    const out = { list: m, crossTime: NaN };
    m.push({ label: 'Swing', value: fmtV(vmax - vmin) });
    m.push({ label: 'High / low level', value: `${fmtV(vmax)} / ${fmtV(vmin)}` });
    if (!(vmax - vmin > 0)) return out;
    const mid = 0.5 * (vmax + vmin);
    // Two periods so crossings near the period boundary are found.
    const t2 = new Float64Array(2 * n), v2 = new Float64Array(2 * n);
    for (let i = 0; i < 2 * n; i++) { t2[i] = (i < n ? t[i] : t[i - n] + P); v2[i] = v[i % n]; }
    const tr = crossing(t2, v2, mid, 1, 2 * n, 1);
    const tf = crossing(t2, v2, mid, indexAt(t2, tr) + 1, 2 * n, -1);
    if (!isFinite(tr) || !isFinite(tf)) return out;
    const high = tf - tr;
    out.crossTime = tr;
    m.push({ label: 'Duty cycle', value: fmtPct(high / P) });
    m.push({ label: 'Duty-cycle distortion', value: fmtTime(Math.abs(high - P / 2)) });
    const er = edgeTimes(t2, v2, vmin, vmax, Math.max(1, indexAt(t2, tr - 0.5 * P)), 2 * n);
    const ef = edgeTimes(t2, v2, vmax, vmin, indexAt(t2, tr) + 1, 2 * n);
    if (isFinite(er.t1090)) m.push({ label: 'Rise time 10-90 / 20-80 %', value: `${fmtTime(er.t1090)} / ${fmtTime(er.t2080)}` });
    if (isFinite(ef.t1090)) m.push({ label: 'Fall time 10-90 / 20-80 %', value: `${fmtTime(ef.t1090)} / ${fmtTime(ef.t2080)}` });
    // Source 50 % rising time and the delay, unwrapped with the line delay estimate.
    let smin = Infinity, smax = -Infinity;
    for (let i = 0; i < n; i++) { smin = Math.min(smin, s[i]); smax = Math.max(smax, s[i]); }
    const s2 = new Float64Array(2 * n);
    for (let i = 0; i < 2 * n; i++) s2[i] = s[i % n];
    const ts = crossing(t2, s2, 0.5 * (smin + smax), 1, 2 * n, 1);
    if (isFinite(ts)) {
        const tau = Math.min(...ctx.est.modes.map(mm => mm.tau)) + ctx.net.pkgDelay;
        let dly = ((tr - ts) % P + P) % P;
        dly += Math.max(0, Math.round((tau - dly) / P)) * P;
        m.push({ label: 'Delay (50 %)', value: fmtTime(dly), note: 'source to receiver, unwrapped with the line delay' });
    }
    // Peaks against the driver rails (the source levels), and the crossings of the mid
    // level: more than two per period are non-monotonic edges (double clocking).
    if (!ctx.net.current) {
        m.push({ label: 'Maximum / minimum', value: `${fmtV(vmax)} / ${fmtV(vmin)}` });
        m.push({ label: 'Beyond the driver rails (above / below)', value: `${fmtV(Math.max(0, vmax - smax))} / ${fmtV(Math.max(0, smin - vmin))}`,
            note: 'overshoot above the high source level and undershoot below the low one' });
    }
    let nCross = 0;
    for (let i = 0; i < n; i++) if ((v[i] - mid < 0) !== (v[(i + 1) % n] - mid < 0)) nCross++;
    m.push({ label: 'Mid-level crossings per period', value: nCross <= 2 ? `${nCross} (monotonic)` : `${nCross} (non-monotonic edges: double clocking risk)` });
    // Eye of the clock (UI = P/2): the waveform against itself shifted by one UI.
    const ui = P / 2;
    const at = (x) => {
        const xm = ((x % P) + P) % P, i = Math.min(n - 1, Math.floor(xm / (P / n)));
        const j = (i + 1) % n, fr = xm / (P / n) - i;
        return v[i] + (v[j] - v[i]) * fr;
    };
    // Crossing of the rising edge with the shifted (falling) one.
    let tc = tr, best = Infinity;
    for (let k = -200; k <= 200; k++) {
        const x = tr + k * ui / 1000, d = Math.abs(at(x) - at(x + ui));
        if (d < best) { best = d; tc = x; }
    }
    const vc = at(tc);
    m.push({ label: 'Eye crossing level', value: `${fmtV(vc)} (${fmtPct((vc - vmin) / (vmax - vmin))})` });
    const ctr = tc + ui / 2;
    m.push({ label: 'Eye height (centre)', value: fmtV(Math.abs(at(ctr) - at(ctr + ui))) });
    return out;
}

// Clock eye: segments of 2 UI (UI = P/2) starting every UI, aligned so the first
// rising 50 % crossing sits at 0.5 UI. NaN separates the segments.
function foldEye(t, v, P, tCross, count, nFine) {
    if (!isFinite(tCross)) return null;
    const ui = P / 2;
    const xs = [], ys = [];
    const tStart = tCross - 0.5 * ui;
    // A burst folds only its own cycles, not the idle line before and after it.
    let tEnd = t[Math.min(nFine, t.length) - 1];
    if (isFinite(count)) tEnd = Math.min(tEnd, tCross + count * P + 0.5 * ui);
    const dtLocal = t[1] - t[0];
    const stride = Math.max(1, Math.floor(ui / dtLocal / 200));
    // At most the last 400 segments of a long burst (the settled part).
    const kFirst = Math.max(0, Math.floor((tEnd - tStart) / ui) - 401);
    for (let k = kFirst; ; k++) {
        const a = tStart + k * ui, b = a + 2 * ui;
        if (b > tEnd + 1e-15 && k > kFirst) break;
        let i = indexAt(t, a);
        for (; i < nFine && t[i] <= b; i += stride) { xs.push((t[i] - a) / ui); ys.push(v[i]); }
        xs.push(NaN); ys.push(NaN);
        if (b > tEnd) break;
    }
    return { x: Float64Array.from(xs), y: Float64Array.from(ys), ui };
}

// Eye of a steady-state clock: the period from 0.5 UI before the rising crossing and the
// same shifted by one UI, each over 2 UI (UI = P/2).
function foldSteady(v, P, tCross) {
    if (!isFinite(tCross)) return null;
    const n = v.length, ui = P / 2, per = n / 2;
    const xs = [], ys = [];
    const at = (x) => {
        const xm = ((x % P) + P) % P * n / P, i = Math.floor(xm) % n, fr = xm - Math.floor(xm);
        return v[i] + (v[(i + 1) % n] - v[i]) * fr;
    };
    const a = tCross - 0.5 * ui;
    for (const shift of [0, ui]) {
        for (let j = 0; j <= n; j++) { const x = 2 * j / n; xs.push(x); ys.push(at(a + shift + x * ui)); }
        xs.push(NaN); ys.push(NaN);
    }
    return { x: Float64Array.from(xs), y: Float64Array.from(ys), ui, per };
}

// ---------------------------------------------------------------- output packing

// Plot data (min/max decimated on the main trace) and CSV data (uniform stride).
function packOutput(res) {
    const { t, traces } = res;
    const n = t.length;
    const names = Object.keys(traces);
    const main = traces[mainName(res.differential ?? names.includes('rxDiff'))] || traces[names[0]];
    let idx;
    if (n <= 2 * PLOT_BUCKETS) {
        idx = null;
    } else {
        const list = [];
        const per = n / PLOT_BUCKETS;
        for (let b = 0; b < PLOT_BUCKETS; b++) {
            const i0 = Math.floor(b * per), i1 = Math.min(n, Math.floor((b + 1) * per));
            let imin = i0, imax = i0;
            for (let i = i0; i < i1; i++) {
                if (main[i] < main[imin]) imin = i;
                if (main[i] > main[imax]) imax = i;
            }
            if (imin === imax) list.push(imin); else list.push(Math.min(imin, imax), Math.max(imin, imax));
        }
        if (list[list.length - 1] !== n - 1) list.push(n - 1);
        idx = Int32Array.from(list);
    }
    const pick = (a) => idx ? Float64Array.from(idx, i => a[i]) : a;
    res.plot = { t: pick(t), traces: Object.fromEntries(names.map(k => [k, pick(traces[k])])) };
    const stride = Math.max(1, Math.ceil(n / CSV_MAX_ROWS));
    const strided = (a) => { const r = new Float64Array(Math.ceil(n / stride)); for (let i = 0, j = 0; i < n; i += stride, j++) r[j] = a[i]; return r; };
    res.csv = { t: strided(t), traces: Object.fromEntries(names.map(k => [k, strided(traces[k])])) };
}

// ---------------------------------------------------------------- formatting

function fmtTime(s) {
    if (!isFinite(s)) return '-';
    const a = Math.abs(s);
    if (a < 1e-16) return '0 ps';
    if (a >= 1e-3) return `${+(s * 1e3).toPrecision(4)} ms`;
    if (a >= 1e-6) return `${+(s * 1e6).toPrecision(4)} us`;
    if (a >= 1e-9) return `${+(s * 1e9).toPrecision(4)} ns`;
    return `${+(s * 1e12).toPrecision(4)} ps`;
}
function fmtV(v) {
    if (!isFinite(v)) return '-';
    if (Math.abs(v) < 1e-7) return '0 mV';
    return Math.abs(v) >= 1 ? `${+v.toPrecision(4)} V` : `${+(v * 1e3).toPrecision(4)} mV`;
}
const fmtPct = (x) => isFinite(x) ? `${+(100 * x).toPrecision(3)} %` : '-';

export { simulateLink, requiredBandwidth, lineEstimates, reflectionTrips, taper, outputs, deriveTraces, firstOrder,
    slowSplit, packOutput, crossing, N_MAX, fmtTime, fmtV, fmtPct, mainName, srcName };
