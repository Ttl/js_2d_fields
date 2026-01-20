#!/usr/bin/env python
"""
Test suite for microstrip_ref_v2.py

This file tests the new Mesher-based implementation against the original
microstrip_ref.py implementation.
"""

import sys
import numpy as np
from microstrip_ref_v2 import MicrostripSolver2D as MicrostripSolver2D_v2
from microstrip_ref import MicrostripSolver2D as MicrostripSolver2D_v1


def compare_results(name, v1_value, v2_value, tolerance_percent=5.0):
    """Compare two values and print result."""
    diff_percent = abs(v1_value - v2_value) / abs(v1_value) * 100
    status = "✓ PASS" if diff_percent <= tolerance_percent else "✗ FAIL"
    print(f"  {name:25s}: v1={v1_value:8.4f}, v2={v2_value:8.4f}, "
          f"diff={diff_percent:5.2f}% {status}")
    return diff_percent <= tolerance_percent


def test_simple_microstrip(plots=False):
    """Test basic microstrip case."""
    print("\n" + "="*70)
    print("TEST 1: Simple Microstrip (1.6mm substrate, 3mm trace)")
    print("="*70)

    params = dict(
        substrate_height=1.6e-3,
        trace_width=3e-3,
        trace_thickness=35e-6,
        gnd_thickness=35e-6,
        epsilon_r=4.5,
        tan_delta=0.02,
        sigma_cond=5.8e7,
        freq=1e9,
        nx=100, ny=100,
        use_sm=False,
        boundaries=["open", "open", "open", "gnd"]
    )

    print("\nSolving with v1 (original)...")
    solver_v1 = MicrostripSolver2D_v1(**params)
    Z0_v1, eps_eff_v1, C_v1, C0_v1 = solver_v1.calculate_parameters()
    Ex_v1, Ey_v1 = solver_v1.compute_fields()
    alpha_cond_v1, _ = solver_v1.calculate_conductor_loss(Ex_v1, Ey_v1, Z0_v1)
    alpha_diel_v1 = solver_v1.calculate_dielectric_loss(Ex_v1, Ey_v1, Z0_v1)

    print("\nSolving with v2 (Mesher)...")
    params["nx"] = 50
    params["ny"] = 50
    solver_v2 = MicrostripSolver2D_v2(**params)
    #Z0_v2, eps_eff_v2, C_v2, C0_v2 = solver_v2.calculate_parameters()
    Z0_v2, eps_eff_v2, C_v2, C0_v2, Ex, Ey = solver_v2.solve_adaptive()
    Ex_v2, Ey_v2 = solver_v2.compute_fields()
    alpha_cond_v2, _ = solver_v2.calculate_conductor_loss(Ex_v2, Ey_v2, Z0_v2)
    alpha_diel_v2 = solver_v2.calculate_dielectric_loss(Ex_v2, Ey_v2, Z0_v2)

    print("\nComparison (Z0/εeff: 5%, losses: 50%):")
    all_pass = True
    all_pass &= compare_results("Z0 [Ω]", Z0_v1, Z0_v2, 5.0)
    all_pass &= compare_results("ε_eff", eps_eff_v1, eps_eff_v2, 5.0)
    all_pass &= compare_results("α_cond [dB/m]", alpha_cond_v1, alpha_cond_v2, 50.0)
    all_pass &= compare_results("α_diel [dB/m]", alpha_diel_v1, alpha_diel_v2, 50.0)

    if plots:
        print("\nPlotting v1 geometry...")
        solver_v1.plot_geometry()
        print("Plotting v2 geometry...")
        solver_v2.plot_geometry()

    return all_pass


def test_microstrip_with_solder_mask(plots=False):
    """Test microstrip with solder mask."""
    print("\n" + "="*70)
    print("TEST 2: Microstrip with Solder Mask")
    print("="*70)

    params = dict(
        substrate_height=0.1e-3,
        trace_width=0.1e-3,
        trace_thickness=35e-6,
        gnd_thickness=35e-6,
        epsilon_r=4.5,
        tan_delta=0.02,
        sigma_cond=5.8e7,
        freq=1e9,
        nx=100, ny=100,
        use_sm=True,
        boundaries=["open", "open", "open", "gnd"]
    )

    print("\nSolving with v1 (original)...")
    solver_v1 = MicrostripSolver2D_v1(**params)
    Z0_v1, eps_eff_v1, C_v1, C0_v1 = solver_v1.calculate_parameters()
    Ex_v1, Ey_v1 = solver_v1.compute_fields()
    alpha_cond_v1, _ = solver_v1.calculate_conductor_loss(Ex_v1, Ey_v1, Z0_v1)
    alpha_diel_v1 = solver_v1.calculate_dielectric_loss(Ex_v1, Ey_v1, Z0_v1)

    print("\nSolving with v2 (Mesher)...")
    params["nx"] = 50
    params["ny"] = 50
    solver_v2 = MicrostripSolver2D_v2(**params)
    Z0_v2, eps_eff_v2, C_v2, C0_v2, Ex, Ey = solver_v2.solve_adaptive()
    Ex_v2, Ey_v2 = solver_v2.compute_fields()
    alpha_cond_v2, _ = solver_v2.calculate_conductor_loss(Ex_v2, Ey_v2, Z0_v2)
    alpha_diel_v2 = solver_v2.calculate_dielectric_loss(Ex_v2, Ey_v2, Z0_v2)

    print("\nComparison (Z0/εeff: 5%, losses: 50%):")
    all_pass = True
    all_pass &= compare_results("Z0 [Ω]", Z0_v1, Z0_v2, 5.0)
    all_pass &= compare_results("ε_eff", eps_eff_v1, eps_eff_v2, 5.0)
    all_pass &= compare_results("α_cond [dB/m]", alpha_cond_v1, alpha_cond_v2, 50.0)
    all_pass &= compare_results("α_diel [dB/m]", alpha_diel_v1, alpha_diel_v2, 50.0)

    if plots:
        print("\nPlotting v1 geometry...")
        solver_v1.plot_geometry()
        print("Plotting v2 geometry...")
        solver_v2.plot_geometry()

    return all_pass


def test_microstrip_with_gnd_cut(plots=False):
    """Test microstrip with ground cutout."""
    print("\n" + "="*70)
    print("TEST 3: Microstrip with Ground Cutout")
    print("="*70)

    params = dict(
        substrate_height=1.6e-3,
        trace_width=3e-3,
        trace_thickness=35e-6,
        gnd_thickness=35e-6,
        epsilon_r=4.5,
        tan_delta=0.02,
        sigma_cond=5.8e7,
        freq=1e9,
        nx=100, ny=100,
        use_sm=False,
        gnd_cut_width=3e-3,
        gnd_cut_sub_h=1e-3,
        boundaries=["open", "open", "open", "gnd"]
    )

    print("\nSolving with v1 (original)...")
    solver_v1 = MicrostripSolver2D_v1(**params)
    Z0_v1, eps_eff_v1, C_v1, C0_v1 = solver_v1.calculate_parameters()
    Ex_v1, Ey_v1 = solver_v1.compute_fields()
    alpha_cond_v1, _ = solver_v1.calculate_conductor_loss(Ex_v1, Ey_v1, Z0_v1)
    alpha_diel_v1 = solver_v1.calculate_dielectric_loss(Ex_v1, Ey_v1, Z0_v1)

    print("\nSolving with v2 (Mesher)...")
    params["nx"] = 50
    params["ny"] = 50
    solver_v2 = MicrostripSolver2D_v2(**params)
    Z0_v2, eps_eff_v2, C_v2, C0_v2, Ex, Ey = solver_v2.solve_adaptive()
    Ex_v2, Ey_v2 = solver_v2.compute_fields()
    alpha_cond_v2, _ = solver_v2.calculate_conductor_loss(Ex_v2, Ey_v2, Z0_v2)
    alpha_diel_v2 = solver_v2.calculate_dielectric_loss(Ex_v2, Ey_v2, Z0_v2)

    print("\nComparison (Z0/εeff: 5%, losses: 50%):")
    all_pass = True
    all_pass &= compare_results("Z0 [Ω]", Z0_v1, Z0_v2, 5.0)
    all_pass &= compare_results("ε_eff", eps_eff_v1, eps_eff_v2, 5.0)
    all_pass &= compare_results("α_cond [dB/m]", alpha_cond_v1, alpha_cond_v2, 50.0)
    all_pass &= compare_results("α_diel [dB/m]", alpha_diel_v1, alpha_diel_v2, 50.0)

    if plots:
        print("\nPlotting v1 geometry...")
        solver_v1.plot_geometry()
        print("Plotting v2 geometry...")
        solver_v2.plot_geometry()

    return all_pass


def test_embedded_microstrip(plots=False):
    """Test embedded microstrip (trace in dielectric layer)."""
    print("\n" + "="*70)
    print("TEST 4: Embedded Microstrip")
    print("="*70)

    params = dict(
        substrate_height=1.6e-3,
        trace_width=3e-3,
        trace_thickness=35e-6,
        gnd_thickness=35e-6,
        epsilon_r=4.5,
        tan_delta=0.02,
        sigma_cond=5.8e7,
        freq=1e9,
        nx=200, ny=200,
        use_sm=False,
        top_diel_h=0.2e-3,
        top_diel_er=4.5,
        boundaries=["open", "open", "open", "gnd"]
    )

    print("\nSolving with v1 (original)...")
    solver_v1 = MicrostripSolver2D_v1(**params)
    Z0_v1, eps_eff_v1, C_v1, C0_v1 = solver_v1.calculate_parameters()
    Ex_v1, Ey_v1 = solver_v1.compute_fields()
    alpha_cond_v1, _ = solver_v1.calculate_conductor_loss(Ex_v1, Ey_v1, Z0_v1)
    alpha_diel_v1 = solver_v1.calculate_dielectric_loss(Ex_v1, Ey_v1, Z0_v1)

    print("\nSolving with v2 (Mesher)...")
    params["nx"] = 50
    params["ny"] = 50
    solver_v2 = MicrostripSolver2D_v2(**params)
    #Z0_v2, eps_eff_v2, C_v2, C0_v2 = solver_v2.calculate_parameters()
    Z0_v2, eps_eff_v2, C_v2, C0_v2, Ex, Ey = solver_v2.solve_adaptive()
    Ex_v2, Ey_v2 = solver_v2.compute_fields()
    alpha_cond_v2, _ = solver_v2.calculate_conductor_loss(Ex_v2, Ey_v2, Z0_v2)
    alpha_diel_v2 = solver_v2.calculate_dielectric_loss(Ex_v2, Ey_v2, Z0_v2)

    print("\nComparison (Z0/εeff: 5%, losses: 50%):")
    all_pass = True
    all_pass &= compare_results("Z0 [Ω]", Z0_v1, Z0_v2, 5.0)
    all_pass &= compare_results("ε_eff", eps_eff_v1, eps_eff_v2, 5.0)
    all_pass &= compare_results("α_cond [dB/m]", alpha_cond_v1, alpha_cond_v2, 50.0)
    all_pass &= compare_results("α_diel [dB/m]", alpha_diel_v1, alpha_diel_v2, 50.0)

    if plots:
        print("\nPlotting v1 geometry...")
        solver_v1.plot_geometry()
        print("Plotting v2 geometry...")
        solver_v2.plot_geometry()

    return all_pass


def test_stripline(plots=False):
    """Test stripline configuration (top and bottom ground)."""
    print("\n" + "="*70)
    print("TEST 5: Stripline")
    print("="*70)

    params = dict(
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
        boundaries=["open", "open", "gnd", "gnd"]
    )

    print("\nSolving with v1 (original)...")
    solver_v1 = MicrostripSolver2D_v1(**params)
    Z0_v1, eps_eff_v1, C_v1, C0_v1 = solver_v1.calculate_parameters()
    Ex_v1, Ey_v1 = solver_v1.compute_fields()
    alpha_cond_v1, _ = solver_v1.calculate_conductor_loss(Ex_v1, Ey_v1, Z0_v1)
    alpha_diel_v1 = solver_v1.calculate_dielectric_loss(Ex_v1, Ey_v1, Z0_v1)

    print("\nSolving with v2 (Mesher)...")
    solver_v2 = MicrostripSolver2D_v2(**params)
    Z0_v2, eps_eff_v2, C_v2, C0_v2 = solver_v2.calculate_parameters()
    Ex_v2, Ey_v2 = solver_v2.compute_fields()
    alpha_cond_v2, _ = solver_v2.calculate_conductor_loss(Ex_v2, Ey_v2, Z0_v2)
    alpha_diel_v2 = solver_v2.calculate_dielectric_loss(Ex_v2, Ey_v2, Z0_v2)

    print("\nComparison (Z0/εeff: 5%, losses: 50%):")
    all_pass = True
    all_pass &= compare_results("Z0 [Ω]", Z0_v1, Z0_v2, 5.0)
    all_pass &= compare_results("ε_eff", eps_eff_v1, eps_eff_v2, 5.0)
    all_pass &= compare_results("α_cond [dB/m]", alpha_cond_v1, alpha_cond_v2, 50.0)
    all_pass &= compare_results("α_diel [dB/m]", alpha_diel_v1, alpha_diel_v2, 50.0)

    if plots:
        print("\nPlotting v1 geometry...")
        solver_v1.plot_geometry()
        print("Plotting v2 geometry...")
        solver_v2.plot_geometry()

    return all_pass


def run_all_tests(plots=False):
    """Run all tests and summarize results."""
    print("\n" + "#"*70)
    print("# MICROSTRIP_REF_V2 TEST SUITE")
    print("#"*70)

    tests = [
        ("Simple Microstrip", test_simple_microstrip),
        ("Solder Mask", test_microstrip_with_solder_mask),
        ("Ground Cutout", test_microstrip_with_gnd_cut),
        ("Embedded Microstrip", test_embedded_microstrip),
        ("Stripline", test_stripline),
    ]

    results = []
    for name, test_func in tests:
        try:
            passed = test_func(plots=plots)
            results.append((name, passed))
        except Exception as e:
            print(f"\n✗ TEST FAILED with exception: {e}")
            import traceback
            traceback.print_exc()
            results.append((name, False))

    # Summary
    print("\n" + "="*70)
    print("TEST SUMMARY")
    print("="*70)
    for name, passed in results:
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"  {name:30s}: {status}")

    total = len(results)
    passed = sum(1 for _, p in results if p)
    print(f"\nTotal: {passed}/{total} tests passed")
    print("="*70 + "\n")

    return passed == total


if __name__ == "__main__":
    plots = False
    if len(sys.argv) > 1:
        if "plot" in sys.argv[1]:
            plots = True

    success = run_all_tests(plots=plots)
    sys.exit(0 if success else 1)
