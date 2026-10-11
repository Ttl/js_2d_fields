// Equalization: transmitter FIR (pre-shoot / de-emphasis), receiver CTLE and DFE, and
// the peak-distortion eye estimate used to choose their settings.

// PCIe 3.0+ transmitter presets: [pre-shoot dB, de-emphasis dB]. P10's de-emphasis is
// set by the transmitter's full-swing and low-frequency limits; -9.5 dB is typical.
const PCIE_PRESETS = {
    P0: [0, -6.0], P1: [0, -3.5], P2: [0, -4.4], P3: [0, -2.5], P4: [0, 0],
    P5: [1.9, 0], P6: [2.5, 0], P7: [3.5, -6.0], P8: [3.5, -3.5], P9: [3.5, 0], P10: [0, -9.5],
};

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
 * CTLE of the PCIe 3.0 reference form: H(s) = wp2 (s + A wp1) / ((s + wp1)(s + wp2)), a
 * DC gain A (adcDb), a zero at A fp1, poles at fp1 and fp2. Returns f -> [re, im].
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
 * Peak-distortion eye estimate from a sampled pulse response q (the response to +1
 * minus that to -1, halved: the signal is base + sum_m d_m q(t - m UI), d = +-1).
 * The inner eye height at sampling index n0 is 2 (q(n0) - sum |q(n0 + k spb)|) over the
 * cursors not cancelled by a DFE of `dfeTaps` taps (k = 1..dfeTaps). The best index is
 * searched within half a UI of the pulse peak. Returns { height, n0, cursors }.
 */
function pdaEye(q, spb, dfeTaps = 0, span = 64) {
    const N = q.length;
    let peak = 0;
    for (let n = 1; n < N; n++) if (q[n] > q[peak]) peak = n;
    let best = -Infinity, bestN = peak;
    const half = Math.floor(spb / 2);
    for (let n0 = Math.max(0, peak - half); n0 <= Math.min(N - 1, peak + half); n0++) {
        let isi = 0;
        for (let k = -span; k <= span; k++) {
            if (k === 0 || (k >= 1 && k <= dfeTaps)) continue;
            const i = n0 + k * spb;
            if (i >= 0 && i < N) isi += Math.abs(q[i]);
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

export { PCIE_PRESETS, firFromPresetDb, presetDbFromFir, ctleTransfer, pdaEye };
