"""COM-method reference for tests/test_eye_com.mjs.

The statistical method of IEEE 802.3 Annex 93A (COM), as implemented by PyChOpMarg: the
interference distribution is the convolution of one +-h two-point distribution per
cursor (delta_pmf, 93A-39/40), the DFE taps are the post-cursors at the sample point
clipped to their limits (93A-26/27), and the jitter of COM is a Gaussian (sigma_Rj) and a
dual Dirac (A_DD) on the pulse slope per symbol (calc_hJ, 93A-28/32). No noise, no
crosstalk, no packages, no Rx FFE: the parts that compare with the Eye engine's
statistical eye.

Reads a JSON case on stdin, writes JSON to stdout:
  - "ours": the eye height and width at the target BER of the Eye engine's pulse response
    (as given, sample point and DFE taps from its own rules), and the DFE taps;
  - "pulse", "ts": the independent pulse response: Howard Johnson's line model (PyBERT's
    calc_gamma) as S21 between matched 100 ohm ends, a Gaussian Tx edge of the COM form
    (93A-46, T_r the 20-80 % time), the Tx FFE and CTLE of PyChOpMarg, the pulse by
    (93A-24);
  - "ref": eye height and width of that pulse at the same sample point (aligned to ours by
    the best time shift);
  - "jitter": the COM jitter terms on our pulse (when the case has jitter), for comparison
    only: COM moves each symbol's pulse, the PCIe model each transition's edge.
"""

import json
import sys

import numpy as np
from scipy.signal import fftconvolve
from pybert.utility.channel import calc_gamma
from pychopmarg.utility.filter import calc_Hctle, calc_Hffe
from pychopmarg.utility.probability import delta_pmf


def interp(x, y):
    "Linear interpolation of samples y at the points x, 0 outside."
    return lambda t: np.interp(t, x, y, left=0.0, right=0.0)


def open_eye(q_at, k_lo, k_hi, dfe, ber, nbins, dfe_window, jitter=None):
    """Eye height at phase 0 and eye width at the 0 threshold, for the pulse q_at(t) (t in
    UI from the sample point), cursors k_lo..k_hi, fixed DFE taps (subtracted from the
    first post-cursors at every phase, as the slicer holds bit 0's feedback in its window).
    BER(v) = P(V1 < v) / 2 + P(V0 > v) / 2 with V = h0 d0 + interference."""

    def isi_pmf(phase):
        h = np.array([q_at(phase + k) for k in range(k_lo, k_hi + 1)])
        h0 = h[-k_lo]
        h = np.delete(h, -k_lo)
        post = h[-k_lo:]
        post[:len(dfe)] -= dfe
        ymax = 1.02 * (np.abs(h).sum() + abs(h0)) + 1e-9
        y = np.linspace(-ymax, ymax, 2 * nbins + 1)
        _, p = delta_pmf(h, L=2, y=y)
        if jitter:
            p = np.maximum(fftconvolve(p, jitter(phase, y), mode="same"), 0.0)
            p /= p.sum()
        return h0, y, p

    # Eye height: P(h0 + I < top) = 2 BER (the zero distribution's tail is far away).
    h0, y, p = isi_pmf(0.0)
    cdf = np.cumsum(p)
    i = np.searchsorted(cdf, 2 * ber)
    top = h0 + y[i] - 0.5 * (y[1] - y[0])
    height = max(0.0, 2 * top)

    def ber_at(phase):
        h0, y, p = isi_pmf(phase)
        return p[y < -h0].sum() + 0.5 * p[y == -h0].sum()

    def edge(direction):
        reach = 0.5 if dfe_window else 1.0
        if dfe_window and ber_at(direction * (0.5 - 1e-6)) < ber:
            return direction * 0.5
        a, b = 0.0, reach
        for _ in range(14):
            m = 0.5 * (a + b)
            if ber_at(direction * m) < ber:
                a = m
            else:
                b = m
        return direction * 0.5 * (a + b)

    width = edge(1) - edge(-1) if height > 0 else 0.0
    return height, width


def dfe_taps(q_at, n, limits):
    h0 = q_at(0.0)
    taps = []
    for k in range(1, n + 1):
        h = q_at(float(k))
        lim = limits[k - 1] if limits and k - 1 < len(limits) else None
        taps.append(float(np.clip(h, -lim, lim)) if lim else h)
    del h0
    return np.array(taps)


def line_table(ln, f_max, per_dec=40):
    "R, L, G, C per unit length of the Johnson line from its gamma and Zc (1 Hz for DC)."
    f = np.concatenate([[1.0], np.logspace(3, np.log10(f_max), int(per_dec * (np.log10(f_max) - 3)) + 1)])
    w = 2 * np.pi * f
    gamma, zc = calc_gamma(ln["R0"], ln["w0"], ln["Rdc"], ln["Z0"], ln["v0"] * 3.0e8, ln["Theta0"],
                           np.concatenate([[0.0], w]))
    z, y = gamma[1:] * zc[1:], gamma[1:] / zc[1:]
    return {"freqs": [0.0] + list(f[1:]), "R": list(z.real), "L": list(z.imag / w),
            "G": [0.0] + list(y[1:].real), "C": list(y.imag / w), "L_ext": ln["Z0"] / (ln["v0"] * 3.0e8)}


def main() -> None:
    case = json.load(sys.stdin)
    if case.get("table_only"):
        json.dump(line_table(case["line"], case["f_max"]), sys.stdout)
        return
    ui = 1 / case["bit_rate"]
    ber, nbins = case["ber"], case.get("nbins", 20000)
    ndfe, limits = case["dfe_taps"], case.get("dfe_limits")

    # Ours: the Eye engine's pulse response, x in UI from its sample point.
    xo, qo = np.array(case["x"]), np.array(case["q"])
    q_ours = interp(xo, qo)
    k_lo, k_hi = int(np.ceil(xo[0])), int(np.floor(xo[-1]))
    dfe_o = dfe_taps(q_ours, ndfe, limits)
    eh_o, ew_o = open_eye(q_ours, k_lo, k_hi, dfe_o, ber, nbins, ndfe > 0)
    result = {"ours": {"height": eh_o, "width": ew_o, "dfe": list(dfe_o)}}

    # Independent pulse response, by (93A-24).
    M, fstep = case["M"], case["fstep"]
    ts = ui / M
    t = np.arange(0, 1 / fstep, ts)
    f = np.arange(len(t) // 2 + 1) * fstep
    fl = f.copy()
    fl[0] = 1.0
    ln = case["line"]
    # calc_gamma() overwrites its first point, so it gets a dummy one.
    gamma, zc = calc_gamma(ln["R0"], ln["w0"], ln["Rdc"], ln["Z0"], ln["v0"] * 3.0e8, ln["Theta0"],
                           np.concatenate([[0.0], 2 * np.pi * fl]))
    gamma, zc = gamma[1:], zc[1:]
    gl, zr = gamma * ln["l_ch"], case["z_ref"]
    s21 = 2 / (2 * np.cosh(gl) + (zc / zr + zr / zc) * np.sinh(gl))
    ht = np.exp(-2 * (np.pi * f * case["T_r"] / 1.6832) ** 2)              # (93A-46)
    htx = calc_Hffe(f, ui, np.array(case["fir"]), 1, hasCurs=True)         # (93A-21)
    H = case["A_v"] * ht * s21 * htx
    if case.get("ctle"):
        c = case["ctle"]
        # fz = fp1: the PCIe 3.0 form wp2 (s + A wp1) / ((s + wp1)(s + wp2)), A = gDC.
        H = H * calc_Hctle(f, c["fp1"], c["fp1"], c["fp2"], 1e9, c["adcDb"], 0.0)   # (93A-22)
    xsinc = M * np.sinc(ui * f)
    p = np.fft.irfft(xsinc * H, len(t))

    # Align on ours: the shift (UI) minimising the largest difference over -3..20 UI.
    tr = np.arange(len(p)) * ts / ui
    sel = (xo > -3) & (xo < 20)
    peak_r, peak_o = tr[np.argmax(p)], xo[np.argmax(qo)]
    shifts = peak_r - peak_o + np.linspace(-0.05, 0.05, 1001)
    errs = [np.max(np.abs(np.interp(xo[sel] + s, tr, p) - qo[sel])) for s in shifts]
    s_best = shifts[int(np.argmin(errs))]
    q_ref = lambda x: np.interp(np.asarray(x) + s_best, tr, p, left=0.0, right=0.0)
    span = min(case.get("ref_span", 200), int(np.floor((tr[-1] - s_best) * 0.9)))
    dfe_r = dfe_taps(q_ref, ndfe, limits)
    eh_r, ew_r = open_eye(q_ref, k_lo, span, dfe_r, ber, nbins, ndfe > 0)
    tail = float(np.sqrt(sum(q_ref(float(k)) ** 2 for k in range(k_hi + 1, span + 1))))
    result["ref"] = {"height": eh_r, "width": ew_r, "dfe": list(dfe_r), "pulse_diff": float(min(errs)),
                     "peak": float(np.max(p)), "tail_rms": tail}

    # COM jitter terms on our pulse at the sample point: sigma_Rj and A_DD (UI) on the
    # slope of every cursor, 93A-28/32/40/42.
    if case.get("jitter"):
        j = case["jitter"]
        rj, add = j["rj"] / ui, 0.5 * j["djdd"] / ui
        d = 1 / 64

        def jit(phase, y):
            hj = np.array([(q_ours(phase + k + d) - q_ours(phase + k - d)) / (2 * d) for k in range(k_lo, k_hi + 1)])
            var = rj ** 2 * (hj ** 2).sum()
            g = np.exp(-y ** 2 / (2 * var)) if var > 0 else (y == 0).astype(float)
            _, pj = delta_pmf(add * hj, L=2, y=y)
            out = np.maximum(fftconvolve(g / g.sum(), pj, mode="same"), 0.0)
            return out / out.sum()

        eh_j, ew_j = open_eye(q_ours, k_lo, k_hi, dfe_o, ber, nbins, ndfe > 0, jitter=jit)
        result["com_jitter"] = {"height": eh_j, "width": ew_j}

    # Monte Carlo of the PCIe jitter model on the independent pulse: random bits, the Tx FIR
    # levels L_b = c-1 d_(b+1) + c0 d_b + c+1 d_(b-1), every level change moving with the
    # jitter of its boundary (Gaussian Rj and dual Dirac +-DjDD/2, independent per boundary)
    # through the channel step response S without the FIR:
    #   V(s) = sum_b (L_b - L_(b-1)) S(s - b - delta_b) - sum_k h_k d_-k  (s, b, delta in UI).
    if case.get("mc"):
        mc, j = case["mc"], case.get("jitter") or {"rj": 0.0, "djdd": 0.0}
        rj, add = j["rj"] / ui, 0.5 * j["djdd"] / ui
        p_n = np.fft.irfft(xsinc * H / htx, len(t))
        qn = lambda x: np.interp(np.asarray(x) + s_best - 1, tr, p_n, left=0.0, right=0.0)
        b_lo, b_hi = -mc.get("post", 80), -k_lo + 2
        tg = np.arange(-b_hi - 2, -b_lo + 2, 1 / 512)
        Sg = sum(qn(tg - m) for m in range(0, -b_lo + 40))
        s_inf = Sg[-1]
        c = case["fir"]
        nb = b_hi - b_lo + 1

        def sample(phase, n, seed):
            rng = np.random.default_rng(seed)
            d = 2.0 * rng.integers(0, 2, size=(n, nb + 3)) - 1.0      # bits b_lo - 2 .. b_hi + 1
            lev = c[0] * d[:, 2:] + c[1] * d[:, 1:-1] + c[2] * d[:, :-2]   # levels b_lo - 1 .. b_hi
            delta = np.zeros((n, nb))
            if rj:
                delta += rng.normal(0.0, rj, (n, nb))
            if add:
                delta += add * (2.0 * rng.integers(0, 2, size=(n, nb)) - 1.0)
            # S on its uniform grid, linear between nodes.
            x = (phase - np.arange(b_lo, b_hi + 1)[None, :] - delta - tg[0]) * 512
            i = np.clip(x.astype(np.int64), 0, len(Sg) - 2)
            fr = x - i
            s_at = Sg[i] + (Sg[i + 1] - Sg[i]) * fr
            v = lev[:, 0] * s_inf + (np.diff(lev, axis=1) * s_at).sum(axis=1)
            d0 = d[:, -b_lo + 2]
            for k, h in enumerate(dfe_r):
                v -= h * d[:, -b_lo + 1 - k]
            return v, d0

        def chunks(phase, n, seed):
            vs, ds = [], []
            for i in range(0, n, 20000):
                v, d0 = sample(phase, min(20000, n - i), seed + i)
                vs.append(v)
                ds.append(d0)
            return np.concatenate(vs), np.concatenate(ds)

        v, d0 = chunks(0.0, mc["n"], 1)
        ones = v[d0 > 0]
        out = {"probs": mc["probs"], "top": [float(np.quantile(ones, q)) for q in mc["probs"]], "n": int(len(ones))}
        if mc.get("ew_ber"):
            target, n_ew = mc["ew_ber"], mc.get("n_ew", 200000)

            def ber_at(ph):
                v, d0 = chunks(ph, n_ew, 7)
                return float(np.mean(np.sign(v) != d0))

            def edge(direction):
                a, b = 0.0, 0.5 if ndfe > 0 else 1.0
                if ndfe > 0 and ber_at(direction * (0.5 - 1e-6)) < target:
                    return direction * 0.5
                for _ in range(9):
                    m = 0.5 * (a + b)
                    if ber_at(direction * m) < target:
                        a = m
                    else:
                        b = m
                return direction * 0.5 * (a + b)

            out["width"] = edge(1) - edge(-1)
        result["mc"] = out
    json.dump(result, sys.stdout)


if __name__ == "__main__":
    main()
