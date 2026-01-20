#!/usr/bin/env python
import sys
import numpy as np
from microstrip_ref_v2 import MicrostripSolver2D

def force_symmetric(arr):
    a = np.asarray(arr, dtype=float).copy()
    n = a.size

    # First half indices
    i = np.arange(n // 2)
    j = n - 1 - i

    avg = 0.5 * (a[i] + a[j])
    a[i] = avg
    a[j] = avg

    return a

def solve_microstrip(plots=True):
    print("Solving Microstrip (v2 with Mesher)...")
    solver = MicrostripSolver2D(
        substrate_height=1.6e-3,
        trace_width=3e-3,
        trace_thickness=35e-6,
        gnd_thickness=35e-6,
        epsilon_r=4.5,
        tan_delta=0.02,
        sigma_cond=5.8e7,
        freq=1e9,
        nx=20, ny=20,
        use_sm=False,
        boundaries=["open", "open", "open", "gnd"]
    )

    # Print X mesh
    print("X mesh", solver.x)
    print("X mesh diff", np.diff(solver.x))
    print("X mesh diff forced symmetric", force_symmetric(np.diff(solver.x)))

    Z0, eps_eff, C, C0 = solver.calculate_parameters()
    Ex, Ey = solver.compute_fields()
    #Z0, eps_eff, C, C0, Ex, Ey = solver.solve_adaptive(param_tol=0.001)

    alpha_cond, J = solver.calculate_conductor_loss(Ex, Ey, Z0)
    alpha_diel = solver.calculate_dielectric_loss(Ex, Ey, Z0)
    alpha_total = alpha_cond + alpha_diel

    z, rlgc, eps_eff_mode = solver.rlgc(alpha_cond, alpha_diel, C, Z0)
    print(f"Z (complex) = {z:.2f} ohm, eps_eff {eps_eff_mode:.3f}, RLGC {rlgc}")

    print(f"\n{'='*50}")
    print(f"MICROSTRIP ANALYSIS RESULTS")
    print(f"{'='*50}")
    print(f"Characteristic Impedance Z0:  {Z0:.2f} Ω")
    print(f"Effective Permittivity εᵣₑff: {eps_eff:.3f}")
    print(f"Losses (dB/m) @ 1GHz:         Diel={alpha_diel:.4f}, Cond={alpha_cond:.4g}")
    print(f"Total Attenuation:            {alpha_total:.4f} dB/m")
    print(f"{'='*50}\n")

    if plots:
        solver.plot_geometry()
        solver.plot(Ex, Ey)


if __name__ == "__main__":
    plots = False
    if len(sys.argv) > 1:
        if "plot" in sys.argv[1]:
            plots = True
    solve_microstrip(plots)
