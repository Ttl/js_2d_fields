// Driver stimulus: edge shapes, their spectra, and the edge lists of a step, a single
// pulse or a clock.
//
// An edge goes from 0 to 1 with its 50 % point at t = 0. The shapes are odd-symmetric
// about that point, so the spectrum of the edge derivative (unit area) is real.

const RAMP = 'ramp', RCOS = 'rcos', GAUSS = 'gauss';

// Gaussian quantiles: the 10-90 % width is 2 * 1.2816 sigma, the 20-80 % one 2 * 0.8416 sigma.
const Z90 = 1.2815515655446004, Z80 = 0.8416212335729143;

// Complementary error function, fractional error below 1.2e-7 (Numerical Recipes erfcc).
function erfc(x) {
    const z = Math.abs(x), t = 1 / (1 + 0.5 * z);
    const r = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 +
        t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 +
        t * (-0.82215223 + t * 0.17087277)))))))));
    return x >= 0 ? r : 2 - r;
}

/**
 * Edge of a given shape whose 10-90 % (def '10-90') or 20-80 % (def '20-80') transition
 * time is `time`.
 * Returns { shape, w, half, t1090 }: w is the full transition time (ramp, raised cosine)
 * or sigma (Gaussian), half the half-width of the support (the Gaussian is cut at 5 sigma,
 * 3e-7 of the step), t1090 the 10-90 % time.
 */
function makeEdge(shape, time, def = '10-90') {
    if (!(time > 0)) throw new Error('Rise and fall times must be positive.');
    const is1090 = def !== '20-80';
    if (shape === RAMP) {
        const w = time / (is1090 ? 0.8 : 0.6);
        return { shape, w, half: w / 2, t1090: 0.8 * w };
    }
    if (shape === RCOS) {
        // 0.5 + 0.5 sin(pi t / w): the 10 % and 90 % points sit at sin = -+0.8.
        const f1090 = 2 * Math.asin(0.8) / Math.PI, f2080 = 2 * Math.asin(0.6) / Math.PI;
        const w = time / (is1090 ? f1090 : f2080);
        return { shape, w, half: w / 2, t1090: f1090 * w };
    }
    if (shape === GAUSS) {
        const s = time / (2 * (is1090 ? Z90 : Z80));
        return { shape, w: s, half: 5 * s, t1090: 2 * Z90 * s };
    }
    throw new Error(`Unknown edge shape "${shape}".`);
}

// Edge value (0..1) at time t from its 50 % point.
function edgeValue(e, t) {
    if (e.shape === GAUSS) return 0.5 * erfc(-t / (e.w * Math.SQRT2));
    if (t <= -e.half) return 0;
    if (t >= e.half) return 1;
    if (e.shape === RAMP) return 0.5 + t / e.w;
    return 0.5 + 0.5 * Math.sin(Math.PI * t / e.w);
}

// Fourier transform of the edge derivative at f (real, 1 at f = 0).
function edgeSpectrum(e, f) {
    if (e.shape === RAMP) {
        const x = Math.PI * f * e.w;
        return x === 0 ? 1 : Math.sin(x) / x;
    }
    if (e.shape === RCOS) {
        const y = 2 * f * e.w, d = 1 - y * y;
        return Math.abs(d) < 1e-9 ? Math.PI / 4 : Math.cos(Math.PI * f * e.w) / d;
    }
    const s = e.w;
    return Math.exp(-2 * Math.PI * Math.PI * s * s * f * f);
}

/**
 * Stimulus of the driver legs (one leg single-ended, P and N for a pair).
 *
 * cfg: { type: 'step' | 'pulse' | 'clock' | 'pattern', differential, vLow, vHigh, tr, tf,
 *        edgeDef, shape, delay, width (pulse), clockFreq, duty, steady (clock steady
 *        state), cycles (clock burst), bitRate, bits (pattern, Uint8Array), fir (pattern
 *        transmitter taps [c-1, c0, c+1], |c| summing to 1), skew (N leg delay) }
 *
 * Voltages are the open-circuit source levels. Each leg is v0 plus a list of signed
 * steps { t, dv, edge } (t the 50 % time). A clock lists the edges of its first period
 * and repeats them `count` times every `period` (steady: indefinitely). A bit pattern
 * is periodic like a steady clock: its transitions sit on the bit boundaries k / bitRate.
 */
function buildStimulus(cfg) {
    const dV = cfg.vHigh - cfg.vLow;
    if (!(Math.abs(dV) > 0)) throw new Error('The high and low levels must differ.');
    const up = makeEdge(cfg.shape, cfg.tr, cfg.edgeDef);
    const down = makeEdge(cfg.shape, cfg.tf, cfg.edgeDef);
    // A level change dv rises with the rise time when it goes up.
    const edgeFor = (dv) => dv > 0 ? up : down;
    const skew = cfg.differential ? (cfg.skew || 0) : 0;
    const halfMax = Math.max(up.half, down.half);
    const delay = Math.max(0, cfg.delay || 0);
    const t0 = delay + halfMax + Math.max(0, -skew);

    // Data transition times of one period and the level change of the P leg at each.
    let times, period = null, count = 1, steady = false;
    if (cfg.type === 'pattern') return patternStimulus(cfg, { dV, up, down, edgeFor, skew, halfMax });
    if (cfg.type === 'step') {
        times = [t0];
    } else if (cfg.type === 'pulse') {
        if (!(cfg.width > 0)) throw new Error('Pulse width must be positive.');
        times = [t0, t0 + cfg.width];
    } else if (cfg.type === 'clock') {
        if (!(cfg.clockFreq > 0)) throw new Error('Clock frequency must be positive.');
        const duty = cfg.duty ?? 0.5;
        if (!(duty > 0 && duty < 1)) throw new Error('Duty cycle must be between 0 and 100 %.');
        period = 1 / cfg.clockFreq;
        times = [t0, t0 + duty * period];
        steady = !!cfg.steady;
        count = steady ? Infinity : Math.max(1, Math.round(cfg.cycles || 1));
        // Edges overlap when the transitions are longer than the high or low time.
        const minHalfTime = Math.min(duty, 1 - duty) * period;
        if (2 * halfMax > minHalfTime) {
            throw new Error('The rise and fall times are too long for this clock: the edges overlap.');
        }
    } else {
        throw new Error(`Unknown stimulus type "${cfg.type}".`);
    }
    if (cfg.type === 'pulse' && 2 * halfMax > cfg.width) {
        throw new Error('The pulse is shorter than its edges.');
    }

    const legEdges = (v0, sign, shift) => ({
        v0,
        edges: times.map((t, k) => {
            const dv = sign * (k % 2 === 0 ? dV : -dV);
            return { t: t + shift, dv, edge: edgeFor(dv) };
        }),
    });
    const legs = [legEdges(cfg.vLow, 1, 0)];
    if (cfg.differential) legs.push(legEdges(cfg.vHigh, -1, skew));

    // Last edge end of a finite stimulus.
    let tEnd = 0;
    for (const leg of legs) for (const e of leg.edges) tEnd = Math.max(tEnd, e.t + e.edge.half);
    if (count > 1 && isFinite(count)) tEnd += (count - 1) * period;

    return {
        type: cfg.type, differential: !!cfg.differential, legs, period, count, steady,
        t0, tEnd, tStart: delay, dV, up, down,
        edgeMin: Math.min(cfg.tr, cfg.tf),
        t1090Min: Math.min(up.t1090, down.t1090),
        halfMax,
    };
}

// Periodic NRZ bit pattern: bit k occupies [k, k + 1) / bitRate, the transitions at the
// bit boundaries. The period starts at the level of the last bit (v0). A transmitter FIR
// sets the level of bit n to mid + (swing / 2) sum_j c_j d_(n-j), d = +-1, j = -1, 0, 1:
// the low and high levels are its full swing.
function patternStimulus(cfg, { dV, up, down, edgeFor, skew, halfMax }) {
    if (!(cfg.bitRate > 0)) throw new Error('Bit rate must be positive.');
    const bits = cfg.bits;
    if (!bits || bits.length < 2) throw new Error('The pattern needs at least two bits.');
    const ui = 1 / cfg.bitRate, n = bits.length;
    if (2 * halfMax > ui) throw new Error('The rise and fall times are too long for this bit rate: the edges overlap.');
    if (Math.abs(skew) > 0.5 * ui) throw new Error('The P/N skew must stay below half a bit.');
    const fir = cfg.fir || [0, 1, 0];
    const legEdges = (vLo, vHi, shift) => {
        const mid = 0.5 * (vLo + vHi), half = 0.5 * (vHi - vLo);
        const d = (k) => 2 * bits[((k % n) + n) % n] - 1;
        const levels = new Float64Array(n);
        for (let k = 0; k < n; k++) levels[k] = mid + half * (fir[0] * d(k + 1) + fir[1] * d(k) + fir[2] * d(k - 1));
        const edges = [];
        for (let k = 0; k < n; k++) {
            const dv = levels[k] - levels[(k + n - 1) % n];
            if (Math.abs(dv) <= 1e-12 * Math.abs(half)) continue;
            edges.push({ t: k * ui + shift, dv, edge: edgeFor(dv) });
        }
        return { v0: levels[n - 1], edges };
    };
    // Every edge inside [0, period): a negative skew moves the P leg later instead.
    const t0 = Math.max(0, -skew);
    const legs = [legEdges(cfg.vLow, cfg.vHigh, t0)];
    if (cfg.differential) legs.push(legEdges(cfg.vHigh, cfg.vLow, t0 + skew));
    return {
        type: 'pattern', differential: !!cfg.differential, legs, period: n * ui, count: Infinity, steady: true,
        t0, tEnd: n * ui, tStart: 0, dV, up, down, bits, ui, bitRate: cfg.bitRate, fir,
        edgeMin: Math.min(cfg.tr, cfg.tf), t1090Min: Math.min(up.t1090, down.t1090), halfMax,
    };
}

// Level of a leg at time t (finite stimulus; a steady clock is evaluated periodically).
function legValue(stim, leg, t) {
    let v = leg.v0;
    const P = stim.period;
    if (stim.steady) {
        // Periodic: t and the edge times folded into [0, P), each edge with its
        // occurrences one period either side. The level at the period start is v0.
        const tm = ((t % P) + P) % P;
        for (const e of leg.edges) {
            const te = ((e.t % P) + P) % P, d = tm - te;
            v += e.dv * (edgeValue(e.edge, d) + edgeValue(e.edge, d + P) + edgeValue(e.edge, d - P) - 1);
        }
        return v;
    }
    const n = isFinite(stim.count) ? stim.count : 1;
    for (let k = 0; k < n; k++) {
        const shift = k * (P || 0);
        for (const e of leg.edges) {
            const tau = t - e.t - shift;
            if (tau <= -e.edge.half && e.edge.shape !== GAUSS) continue;
            v += e.dv * edgeValue(e.edge, tau);
        }
    }
    return v;
}

// Leg levels on the grid t_i = i * dt, i < n (finite stimulus).
function sampleLeg(stim, leg, dt, n) {
    const out = new Float64Array(n).fill(leg.v0);
    const tail = new Float64Array(n + 1);
    const reps = isFinite(stim.count) ? stim.count : 1;
    for (let k = 0; k < reps; k++) {
        const shift = k * (stim.period || 0);
        for (const e of leg.edges) {
            const tc = e.t + shift, h = e.edge.half;
            const i0 = Math.max(0, Math.ceil((tc - h) / dt));
            const i1 = Math.min(n - 1, Math.floor((tc + h) / dt));
            for (let i = i0; i <= i1; i++) out[i] += e.dv * edgeValue(e.edge, i * dt - tc);
            if (i1 + 1 < n) tail[Math.max(0, i1 + 1)] += e.dv;
        }
    }
    let acc = 0;
    for (let i = 0; i < n; i++) { acc += tail[i]; out[i] += acc; }
    return out;
}

// Level of a leg at the times of an arbitrary increasing grid (finite stimulus).
function sampleLegAt(stim, leg, times) {
    const out = new Float64Array(times.length);
    for (let i = 0; i < times.length; i++) out[i] = legValue(stim, leg, times[i]);
    return out;
}

/**
 * Fourier transform of the leg's derivative at f, as [re, im]: the sum over its edges
 * of dv * E(f) * exp(-j w t), times the period repetition sum of a clock burst.
 * onePeriod: the edges of one period only (steady-state coefficients).
 */
function legDerivSpectrum(stim, leg, f, onePeriod = false) {
    const w = 2 * Math.PI * f;
    let re = 0, im = 0;
    for (const e of leg.edges) {
        const a = e.dv * edgeSpectrum(e.edge, f), ph = -w * e.t;
        re += a * Math.cos(ph); im += a * Math.sin(ph);
    }
    const n = stim.count;
    if (!onePeriod && n > 1 && isFinite(n)) {
        // sum_{k<n} z^k with z = exp(-j w P): (1 - z^n) / (1 - z).
        const th = w * stim.period;
        const dr = 1 - Math.cos(th), di = Math.sin(th);   // 1 - z
        let gr, gi;
        if (dr * dr + di * di < 1e-24) { gr = n; gi = 0; } else {
            const nr = 1 - Math.cos(n * th), ni = Math.sin(n * th);
            const d = dr * dr + di * di;
            gr = (nr * dr + ni * di) / d; gi = (ni * dr - nr * di) / d;
        }
        const r = re * gr - im * gi;
        im = re * gi + im * gr; re = r;
    }
    return [re, im];
}

// Total level change of a leg over the finite stimulus.
function legNetChange(stim, leg) {
    let s = 0;
    for (const e of leg.edges) s += e.dv;
    return s * (isFinite(stim.count) ? stim.count : 1);
}

// Mean level of a steady clock or pattern leg over one period.
function legMean(stim, leg) {
    const P = stim.period;
    // Level changes at the edge times inside the period [0, P): the leg sits at v0
    // before the first and returns to it after the last. The edges are odd-symmetric
    // about their 50 % points, so those points alone fix the mean.
    let m = leg.v0;
    const ts = leg.edges.map(e => ((e.t % P) + P) % P);
    for (let i = 0; i < leg.edges.length; i++) m += leg.edges[i].dv * (P - ts[i]) / P;
    return m;
}

// Highest frequency the sweep has to cover for this stimulus: about three times the
// knee frequency 0.35/t_r of the fastest edge, and five clock harmonics.
function requiredBandwidth(cfg) {
    // Only the edges and the rate matter, a two-bit pattern stands in for the bits.
    const c = cfg.type === 'pattern' ? { ...cfg, bits: Uint8Array.of(0, 1) } : cfg;
    const stim = buildStimulus({ ...c, differential: false, skew: 0 });
    let f = 1 / stim.t1090Min;
    if (stim.type === 'pattern') f = Math.max(f, 2.5 * stim.bitRate);
    else if (stim.period) f = Math.max(f, 5 / stim.period);
    return f;
}

export { requiredBandwidth, makeEdge, edgeValue, edgeSpectrum, buildStimulus, legValue, sampleLeg, sampleLegAt,
    legDerivSpectrum, legNetChange, legMean, erfc, RAMP, RCOS, GAUSS };
