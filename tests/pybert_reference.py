"""PyBERT reference for tests/test_eye_pybert.mjs.

Reads a JSON case on stdin, runs PyBERT's native channel (Howard Johnson's line model,
differential source/load resistances and parasitic capacitances) and writes JSON to
stdout with:
  - the line's R, L, G, C per unit length on a log frequency grid (from PyBERT's own
    gamma and Zc), for the Eye engine to simulate the same line;
  - PyBERT's channel impulse response driven by the given pulse (Gaussian edges of the
    Eye engine), with the Tx FFE taps and, optionally, a CTLE of the PCIe 3.0 reference
    form applied: the reference pulse responses.
"""

import json
import sys

import numpy as np
from scipy.special import erf

from pybert.pybert import PyBERT
from pybert.utility.channel import calc_gamma


def main() -> None:
    case = json.load(sys.stdin)
    bert = PyBERT(run_simulation=False, gui=False)
    bert.bit_rate = case["bit_rate_gbps"]
    bert.nspui = case["nspui"]
    bert.f_max = case["f_max_ghz"]
    bert.f_step = case["f_step_mhz"]
    bert.l_ch, bert.R0, bert.w0, bert.Rdc = case["l_ch"], case["R0"], case["w0"], case["Rdc"]
    bert.Z0, bert.v0, bert.Theta0 = case["Z0"], case["v0"], case["Theta0"]
    bert.rs, bert.cout, bert.rin, bert.cin = case["rs"], case["cout_pf"], case["rin"], case["cin_pf"]

    chnl_h = bert.calc_chnl_h()          # V/sample, referenced to the driver EMF
    t = bert.t
    ts = t[1]
    ui = bert.ui

    # Line parameters from PyBERT's gamma and Zc: Z = gamma Zc, Y = gamma / Zc.
    per_dec = case.get("tab_per_decade", 40)
    f_tab = np.concatenate([[1.0], np.logspace(3, np.log10(case["f_tab_max"]), int(per_dec * (np.log10(case["f_tab_max"]) - 3)) + 1)])
    w_tab = 2 * np.pi * f_tab
    # calc_gamma() overwrites its first point (DC), so it gets a dummy one.
    gamma, zc = calc_gamma(case["R0"], case["w0"], case["Rdc"], case["Z0"], case["v0"] * 3.0e8, case["Theta0"],
                           np.concatenate([[0.0], w_tab]))
    gamma, zc = gamma[1:], zc[1:]
    z = gamma * zc
    y = gamma / zc
    table = {
        "freqs": [0.0] + list(f_tab[1:]),
        # The 1 Hz point stands in for DC (C(f) has a slowly diverging power law at DC).
        "R": list(np.concatenate([[z[0].real], z[1:].real])),
        "L": list(np.concatenate([[z[0].imag / w_tab[0]], z[1:].imag / w_tab[1:]])),
        "G": list(np.concatenate([[0.0], y[1:].real])),
        "C": list(np.concatenate([[y[0].imag / w_tab[0]], y[1:].imag / w_tab[1:]])),
        "L_ext": case["Z0"] / (case["v0"] * 3.0e8),
    }

    # Driver pulse: one UI between the 50 % points, Gaussian edges of the given sigma,
    # starting at t0.
    sigma, t0 = case["edge_sigma"], case["t0"]
    edge = lambda tt: 0.5 * (1 + erf(tt / (sigma * np.sqrt(2))))
    n_resp = len(chnl_h) + int(round((t0 + 4 * ui) / ts))
    tt = np.arange(n_resp) * ts
    pulse = edge(tt - t0) - edge(tt - t0 - ui)

    # Tx FFE [c-1, c0, c+1]: c-1 weights the next bit (one UI earlier).
    c = case["taps"]
    shift = int(round(ui / ts))
    drive = c[1] * pulse
    drive[:-shift] += c[0] * pulse[shift:]
    drive[shift:] += c[2] * pulse[:-shift]
    out = np.convolve(drive, chnl_h)[:n_resp]

    result = {"ts": ts, "ui": ui, "table": table, "pulse": list(out)}
    if case.get("ctle"):
        adc, fp1, fp2 = case["ctle"]["adcDb"], case["ctle"]["fp1"], case["ctle"]["fp2"]
        a, w1, w2 = 10 ** (adc / 20), 2 * np.pi * fp1, 2 * np.pi * fp2
        nfft = 1 << int(np.ceil(np.log2(4 * n_resp)))
        fw = np.fft.rfftfreq(nfft, ts) * 2 * np.pi
        hc = w2 * (1j * fw + a * w1) / ((1j * fw + w1) * (1j * fw + w2))
        result["pulse_ctle"] = list(np.fft.irfft(np.fft.rfft(out, nfft) * hc, nfft)[:n_resp])
    json.dump(result, sys.stdout)


if __name__ == "__main__":
    main()
