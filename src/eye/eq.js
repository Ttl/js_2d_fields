// Equalization: transmitter FIR (pre-shoot / de-emphasis), receiver CTLE and DFE, and
// the peak-distortion eye estimate used to choose their settings.

// PCIe 3.0+ transmitter presets: [pre-shoot dB, de-emphasis dB]. P10 is not fixed: its
// de-emphasis follows the LF level the transmitter advertises in training; -9.5 dB is the
// nominal maximum boost at full swing (+-1.5 dB). Reduced swing implements only
// REDUCED_SWING_PRESETS.

const PCIE_PRESETS = {
    P0: [0, -6.0], P1: [0, -3.5], P2: [0, -4.4], P3: [0, -2.5], P4: [0, 0],
    P5: [1.9, 0], P6: [2.5, 0], P7: [3.5, -6.0], P8: [3.5, -3.5], P9: [3.5, 0], P10: [0, -9.5],
};

const REDUCED_SWING_PRESETS = ['P4', 'P1', 'P9', 'P5', 'P6', 'P3'];

// Boost 20 log(Vd / Vb): the isolated-bit level against a long run's, the quantity the
// PCIe limits are stated in.
function boostDb([cm, c0, cp]) {
    const a = -cm, b = -cp;
    return 20 * Math.log10((a + c0 + b) / (c0 - a - b));
}

/**
 * Taps [c-1, c0, c+1] of a 3-tap FIR from pre-shoot and de-emphasis (dB), normalized to
 * |c-1| + c0 + |c+1| = 1 (the full swing). With a = -c-1, b = -c+1 the levels are
 *   isolated bit 1, first bit after a transition Va = 1 - 2a, bit before a transition
 *   Vc = 1 - 2b, a long run Vb = 1 - 2a - 2b;
 * de-emphasis = 20 log(Vb / Va), pre-shoot = 20 log(Vc / Vb) (PCIe definitions).
 */
function firFromPresetDb(preshootDb, deemphDb) {
    const r1 = Math.pow(10, deemphDb / 20), r2 = Math.pow(10, preshootDb / 20);
    // Va + Vc = 1 + Vb, Vb = r1 Va, Vc = r2 Vb.
    const va = 1 / (1 + r1 * (r2 - 1));
    const vb = r1 * va, vc = r2 * vb;
    const a = (1 - va) / 2, b = (1 - vc) / 2;
    if (!(a >= -1e-12 && b >= -1e-12 && a + b < 0.5)) throw new Error('These pre-shoot and de-emphasis values give no valid FIR taps.');
    return [-Math.max(0, a), 1 - Math.max(0, a) - Math.max(0, b), -Math.max(0, b)];
}

// Pre-shoot and de-emphasis (dB) of taps [c-1, c0, c+1].
function presetDbFromFir([cm, c0, cp]) {
    const a = -cm, b = -cp, n = a + c0 + b;
    const va = (c0 + b - a) / n, vb = (c0 - a - b) / n, vc = (c0 + a - b) / n;
    return { preshoot: 20 * Math.log10(vc / vb), deemphasis: 20 * Math.log10(vb / va) };
}

/**
 * 1st-order CTLE, the PCIe 8 and 16 GT/s behavioral form:
 * H(s) = wp2 (s + A wp1) / ((s + wp1)(s + wp2)), a DC gain A (adcDb), a zero at A fp1,
 * poles at fp1 and fp2. Returns f -> [re, im].
 */
function ctleTransfer(adcDb, fp1, fp2) {
    const A = Math.pow(10, adcDb / 20), w1 = 2 * Math.PI * fp1, w2 = 2 * Math.PI * fp2;
    return (f) => {
        const w = 2 * Math.PI * f;
        // numerator w2 (A w1 + j w), denominator (w1 + j w)(w2 + j w)
        const nr = w2 * A * w1, ni = w2 * w;
        const dr = w1 * w2 - w * w, di = w * (w1 + w2);
        const den = dr * dr + di * di;
        return [(nr * dr + ni * di) / den, (ni * dr - nr * di) / den];
    };
}

/**
 * 2nd-order CTLE, the PCIe 32 GT/s behavioral form:
 * H(s) = (wp1 wp3 wp4 / wz1) (s + wz1)(s + A wp2) / ((s + wp1)(s + wp2)(s + wp3)(s + wp4)),
 * with fp1 = 1.65 fz1 and fp4 = fp3: a DC gain A, a low-frequency zero / pole pair at fz1 /
 * fp1 and the peaking zero at A fp2. Returns f -> [re, im].
 */
function ctle2Transfer(adcDb, fz1, fp2, fp3) {
    const A = Math.pow(10, adcDb / 20), w = (f) => 2 * Math.PI * f;
    const wz1 = w(fz1), wp1 = w(1.65 * fz1), wp2 = w(fp2), wp3 = w(fp3);
    const k = wp1 * wp3 * wp3 / wz1;
    const zeros = [wz1, A * wp2], poles = [wp1, wp2, wp3, wp3];
    return (f) => {
        const x = w(f);
        // product of (a + j x) factors
        let re = k, im = 0;
        const mul = (a, inv) => {
            const d = inv ? a * a + x * x : 1, cr = a / d, ci = (inv ? -x : x) / d;
            const r = re * cr - im * ci; im = re * ci + im * cr; re = r;
        };
        for (const z of zeros) mul(z, false);
        for (const p of poles) mul(p, true);
        return [re, im];
    };
}

// Transfer function of a CTLE setting ({ order, fp1, fp2, fz1, fp3 }) at DC gain adcDb.
function ctleOf(ctle, adcDb) {
    return ctle.order === 2 ? ctle2Transfer(adcDb, ctle.fz1, ctle.fp2, ctle.fp3) : ctleTransfer(adcDb, ctle.fp1, ctle.fp2);
}

// DFE tap k (1-based) for post-cursor value h: h within the tap's limit, a voltage or
// { rel } for |d_k| <= rel * h0, h0 the cursor (absent or non-positive: unlimited).
function dfeTap(h, k, limits, h0) {
    const l = limits && limits[k - 1];
    const lim = l && typeof l === 'object' ? l.rel * Math.abs(h0) : l;
    return lim > 0 ? Math.max(-lim, Math.min(lim, h)) : h;
}

/**
 * Peak-distortion eye estimate from a sampled pulse response q (the response to +1
 * minus that to -1, halved: the signal is base + sum_m d_m q(t - m UI), d = +-1).
 * The inner eye height at sampling index n0 is 2 (q(n0) - sum |residual cursors|): every
 * cursor, the first `dfeTaps` post-cursors less their DFE taps (clamped to `dfeLimits`).
 * The best index is searched within half a UI of the pulse peak. Returns { height, n0,
 * cursors }.
 */
function pdaEye(q, spb, dfeTaps = 0, dfeLimits = null) {
    const span = Math.ceil(q.length / spb);
    const N = q.length;
    let peak = 0;
    for (let n = 1; n < N; n++) if (q[n] > q[peak]) peak = n;
    let best = -Infinity, bestN = peak;
    const half = Math.floor(spb / 2);
    for (let n0 = Math.max(0, peak - half); n0 <= Math.min(N - 1, peak + half); n0++) {
        let isi = 0;
        for (let k = -span; k <= span; k++) {
            if (k === 0) continue;
            const i = n0 + k * spb;
            if (i < 0 || i >= N) continue;
            isi += k >= 1 && k <= dfeTaps ? Math.abs(q[i] - dfeTap(q[i], k, dfeLimits, q[n0])) : Math.abs(q[i]);
        }
        const h = 2 * (q[n0] - isi);
        if (h > best) { best = h; bestN = n0; }
    }
    const cursors = [];
    for (let k = -3; k <= Math.max(8, dfeTaps); k++) {
        const i = bestN + k * spb;
        cursors.push({ k, value: i >= 0 && i < N ? q[i] : 0 });
    }
    return { height: best, n0: bestN, cursors };
}

// Pulse response at a fractional sample index (linear, zero outside).
function pulseAt(q, x) {
    if (x < 0 || x > q.length - 1) return 0;
    const i = Math.floor(x), f = x - i;
    return i + 1 < q.length ? q[i] + (q[i + 1] - q[i]) * f : q[i];
}

/**
 * Peak-distortion eye height at a fractional sample index s with DFE taps h (the post-
 * cursors 1..h.length less their taps, every other cursor as ISI), and the cursors there.
 */
function pdaAt(q, s, spb, h) {
    const span = Math.ceil(q.length / spb);
    const h0 = pulseAt(q, s);
    let isi = 0;
    for (let k = -span; k <= span; k++) {
        if (k === 0) continue;
        const v = pulseAt(q, s + k * spb);
        isi += Math.abs(k >= 1 && k <= h.length ? v - h[k - 1] : v);
    }
    const cursors = [];
    for (let k = -3; k <= Math.max(8, h.length); k++) cursors.push({ k, value: pulseAt(q, s + k * spb) });
    return { height: 2 * (h0 - isi), cursors };
}

export { pulseAt, pdaAt, PCIE_PRESETS, REDUCED_SWING_PRESETS, boostDb, firFromPresetDb, presetDbFromFir, ctleTransfer, ctle2Transfer, ctleOf, pdaEye, dfeTap };
