#!/usr/bin/env python
import sys
import numpy as np
from microstrip_ref_v2 import MicrostripSolver2D
from gcpw_ref_v2 import GroundedCPWSolver2D

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
    if 0:
        solver = MicrostripSolver2D(
            substrate_height=1.6e-3,
            trace_width=3e-3,
            trace_thickness=35e-6,
            gnd_thickness=35e-6,
            epsilon_r=4.5,
            tan_delta=0.02,
            sigma_cond=5.8e7,
            freq=1e9,
            nx=2, ny=2,
            gnd_cut_width=3e-3,
            gnd_cut_sub_h=1e-3,
            use_sm=False,
            boundaries=["open", "open", "open", "gnd"]
        )
    else:
        solver = GroundedCPWSolver2D(
            substrate_height=1.6e-3,
            trace_width=0.3e-3,
            trace_thickness=35e-6,
            gap=0.15e-3,
            top_gnd_width=5e-3,
            via_gap=0.5e-3,
            gnd_thickness=35e-6,
            epsilon_r=4.5,
            tan_delta=0.02,
            sigma_cond=5.8e7,
            freq=1e9,
            nx=2, ny=2,
            use_sm=False,
            boundaries=["open", "open", "open", "gnd"]
        )


    # Print X mesh
    if 1:
        m = float(solver.x[-1] / 2)
        trace_points = [float(x) for x in solver.x if m - 3e-3/2 < x < m + 3e-3/2]
        print(f"Trace points (mm from center): {[round(1e3*(x - m), 3) for x in trace_points]}")
        left_points = [float(x) for x in solver.x if x < m - 3e-3/2]
        print(f"Trace points left of trace (mm from center): {[round(1e3*(x - m), 3) for x in left_points]}")
        print(f"Number of trace points: {len(trace_points)}")
        print("X mesh", solver.x)

    # Trace center
    t = solver.t
    m = solver.y_trace_start + solver.t/2

    trace_y_points = [float(y) for y in solver.y if m - t/2 < y < m + t/2]
    print(f"Trace y-points (mm from center): {[round(1e3*(y - m), 3) for y in trace_y_points]}")

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
