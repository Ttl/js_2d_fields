#!/usr/bin/env python
import sys
from microstrip_ref import *

def solve_microstrip_mask(plots=True):
    print("Solving Microstrip...")
    solver = MicrostripSolver2D(
        substrate_height=0.1e-3,
        trace_width=0.1e-3,
        trace_thickness=35e-6,
        gnd_thickness=35e-6,
        epsilon_r=4.5,
        tan_delta=0.02,
        sigma_cond=5.8e7,
        freq=1e9,
        nx=100, ny=100,
        skin_cells=50,
        use_sm=True,
        boundaries=["open", "open", "open", "gnd"]
    )

    Z0, eps_eff, C, C0 = solver.calculate_parameters()
    Ex, Ey = solver.compute_fields()

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

def solve_microstrip_cut(plots=True):
    print("Solving Microstrip with gnd cut...")
    solver = MicrostripSolver2D(
        substrate_height=1.6e-3,
        trace_width=3e-3,
        trace_thickness=35e-6,
        gnd_thickness=35e-6,
        epsilon_r=4.5,
        tan_delta=0.02,
        sigma_cond=5.8e7,
        freq=1e9,
        nx=100, ny=100,
        skin_cells=50,
        use_sm=False,
        gnd_cut_width=3e-3,
        gnd_cut_sub_h=1e-3,
        boundaries=["open", "open", "open", "gnd"]
    )

    Z0, eps_eff, C, C0 = solver.calculate_parameters()
    Ex, Ey = solver.compute_fields()

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

def solve_microstrip_embed(plots=True):
    print("Solving Embedded Microstrip...")
    solver = MicrostripSolver2D(
        substrate_height=1.6e-3,
        trace_width=3e-3,
        trace_thickness=35e-6,
        gnd_thickness=35e-6,
        epsilon_r=4.5,
        tan_delta=0.02,
        sigma_cond=5.8e7,
        freq=1e9,
        nx=200, ny=200,
        skin_cells=50,
        use_sm=False,
        #gnd_cut_width=3e-3,
        #gnd_cut_sub_h=1e-3,
        top_diel_h=0.2e-3,
        top_diel_er=4.5,
        boundaries=["open", "open", "open", "gnd"]
    )

    Z0, eps_eff, C, C0 = solver.calculate_parameters()
    Ex, Ey = solver.compute_fields()

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

def solve_stripline(plots=True):
    print("Solving Stripline...")
    # Stripline: 0.2mm bottom dielectric, 0.2mm top dielectric
    solver = MicrostripSolver2D(
        substrate_height=0.2e-3,
        trace_width=0.15e-3,
        trace_thickness=35e-6,
        gnd_thickness=16e-6,
        epsilon_r=4.1,
        epsilon_r_top=4.1,
        air_top=0.2e-3,
        tan_delta=0.02,
        sigma_cond=5.8e7,
        freq=1e9,
        nx=100, ny=100,
        skin_cells=50,
        boundaries=["open", "open", "gnd", "gnd"]
    )

    Z0, eps_eff, C, C0 = solver.calculate_parameters()
    Ex, Ey = solver.compute_fields()

    alpha_cond, J = solver.calculate_conductor_loss(Ex, Ey, Z0)
    alpha_diel = solver.calculate_dielectric_loss(Ex, Ey, Z0)
    alpha_total = alpha_cond + alpha_diel

    z, rlgc, eps_eff_mode = solver.rlgc(alpha_cond, alpha_diel, C, Z0)
    print(f"Z (complex) = {z:.2f} ohm, eps_eff {eps_eff_mode:.3f}, RLGC {rlgc}")

    print(f"\n{'='*50}")
    print(f"STRIPLINE ANALYSIS RESULTS")
    print(f"{'='*50}")
    print(f"Characteristic Impedance Z0:  {Z0:.2f} Ω")
    print(f"Effective Permittivity εᵣₑff: {eps_eff:.3f}")
    print(f"Losses (dB/m) @ 1GHz:         Diel={alpha_diel:.4f}, Cond={alpha_cond:.4f}")
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
    solve_microstrip_mask(plots)
    solve_microstrip_cut(plots)
    solve_microstrip_embed(plots)
    solve_stripline(plots)
