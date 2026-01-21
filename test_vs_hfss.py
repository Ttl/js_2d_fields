#!/usr/bin/env python
from microstrip_ref_v2 import MicrostripSolver2D
from gcpw_ref_v2 import GroundedCPWSolver2D

# Alias for backwards compatibility
DifferentialMicrostripSolver2D = MicrostripSolver2D

def test_microstrip_solution(solver_results, reference, test_name="Microstrip"):
    """
    Test microstrip solver results against reference values.

    Args:
        solver_results: dict with keys 'Z0', 'eps_eff', 'alpha_diel', 'alpha_cond', 'C', 'R', 'L', 'G'
        reference: dict with reference values
        test_name: name of the test for printing
    """
    # Global error thresholds (relative error in %)
    MAX_Z0_ERROR = 5.0
    MAX_DIEL_LOSS_ERROR = 10.0
    MAX_COND_LOSS_ERROR = 50.0
    MAX_C_ERROR = 10.0
    MAX_R_ERROR = 15.0
    MAX_L_ERROR = 10.0
    MAX_G_ERROR = 15.0
    MAX_EPS_EFF_ERROR = 5.0

    # Error thresholds mapping
    error_thresholds = {
        'Z0': MAX_Z0_ERROR,
        'diel_loss': MAX_DIEL_LOSS_ERROR,
        'cond_loss': MAX_COND_LOSS_ERROR,
        'C': MAX_C_ERROR,
        'R': MAX_R_ERROR,
        'L': MAX_L_ERROR,
        'G': MAX_G_ERROR,
        'eps_eff': MAX_EPS_EFF_ERROR
    }

    print(f"\n{'='*80}")
    print(f"{test_name.upper()} VALIDATION TEST")
    print(f"{'='*80}")
    print(f"{'Parameter':<15} {'Solved':<15} {'Reference':<15} {'Error (%)':<12} {'Status':<10}")
    print(f"{'-'*80}")

    all_passed = True
    test_results = {}

    for param, ref_value in reference.items():
        if param not in solver_results:
            continue

        solved_value = solver_results[param]

        # Calculate relative error
        if ref_value != 0:
            rel_error = abs((solved_value - ref_value) / ref_value) * 100
        else:
            rel_error = abs(solved_value) * 100

        # Check against threshold
        threshold = error_thresholds.get(param, 10.0)  # default 10%
        passed = rel_error <= threshold
        all_passed = all_passed and passed

        status = "✓ PASS" if passed else "✗ FAIL"

        # Format values based on magnitude
        if param == 'C':
            solved_str = f"{solved_value*1e12:.2f} pF"
            ref_str = f"{ref_value*1e12:.2f} pF"
        elif param == 'L':
            solved_str = f"{solved_value*1e9:.2f} nH"
            ref_str = f"{ref_value*1e9:.2f} nH"
        elif param == 'G':
            solved_str = f"{solved_value*1e3:.2f} mS"
            ref_str = f"{ref_value*1e3:.2f} mS"
        elif param == 'R':
            solved_str = f"{solved_value:.2f} Ω"
            ref_str = f"{ref_value:.2f} Ω"
        elif 'loss' in param:
            solved_str = f"{solved_value:.4f} dB/m"
            ref_str = f"{ref_value:.4f} dB/m"
        else:
            solved_str = f"{solved_value:.3f}"
            ref_str = f"{ref_value:.3f}"

        print(f"{param:<15} {solved_str:<15} {ref_str:<15} {rel_error:<12.2f} {status:<10}")

        test_results[param] = {
            'solved': solved_value,
            'reference': ref_value,
            'error': rel_error,
            'passed': passed
        }

    print(f"{'-'*80}")
    print(f"Overall Result: {'✓ ALL TESTS PASSED' if all_passed else '✗ SOME TESTS FAILED'}")
    print(f"{'='*80}\n")

    assert all_passed, f"{test_name} validation failed - see errors above"

    return test_results


def test_differential_solution(solver_results, reference, test_name="Differential Microstrip"):
    """
    Test differential microstrip solver results against reference values.

    Args:
        solver_results: dict with differential mode parameters
        reference: dict with reference values
        test_name: name of the test for printing
    """
    # Global error thresholds (relative error in %)
    MAX_Z_DIFF_ERROR = 6.0
    MAX_Z_COMMON_ERROR = 6.0
    MAX_Z_ODD_ERROR = 6.0
    MAX_Z_EVEN_ERROR = 6.0
    MAX_EPS_EFF_ERROR = 5.0
    MAX_LOSS_ERROR = 50.0
    MAX_C_ERROR = 10.0

    # Error thresholds mapping
    error_thresholds = {
        'Z_diff': MAX_Z_DIFF_ERROR,
        'Z_common': MAX_Z_COMMON_ERROR,
        'Z_odd': MAX_Z_ODD_ERROR,
        'Z_even': MAX_Z_EVEN_ERROR,
        'eps_eff_odd': MAX_EPS_EFF_ERROR,
        'eps_eff_even': MAX_EPS_EFF_ERROR,
        'C_odd': MAX_C_ERROR,
        'C_even': MAX_C_ERROR,
        'alpha_c_odd': MAX_LOSS_ERROR,
        'alpha_c_even': MAX_LOSS_ERROR,
        'alpha_d_odd': MAX_LOSS_ERROR,
        'alpha_d_even': MAX_LOSS_ERROR,
        'alpha_total_odd': MAX_LOSS_ERROR,
        'alpha_total_even': MAX_LOSS_ERROR,
    }

    print(f"\n{'='*80}")
    print(f"{test_name.upper()} VALIDATION TEST")
    print(f"{'='*80}")
    print(f"{'Parameter':<20} {'Solved':<15} {'Reference':<15} {'Error (%)':<12} {'Status':<10}")
    print(f"{'-'*80}")

    all_passed = True
    test_results = {}

    for param, ref_value in reference.items():
        if param not in solver_results:
            continue

        solved_value = solver_results[param]

        # Calculate relative error
        if ref_value != 0:
            rel_error = abs((solved_value - ref_value) / ref_value) * 100
        else:
            rel_error = abs(solved_value) * 100

        # Check against threshold
        threshold = error_thresholds.get(param, 10.0)
        passed = rel_error <= threshold
        all_passed = all_passed and passed

        status = "✓ PASS" if passed else "✗ FAIL"

        # Format values based on parameter type
        if param.startswith('C_'):
            solved_str = f"{solved_value*1e12:.2f} pF"
            ref_str = f"{ref_value*1e12:.2f} pF"
        elif param.startswith('alpha_'):
            solved_str = f"{solved_value:.4f} dB/m"
            ref_str = f"{ref_value:.4f} dB/m"
        elif param.startswith('Z_'):
            solved_str = f"{solved_value:.2f} Ω"
            ref_str = f"{ref_value:.2f} Ω"
        else:
            solved_str = f"{solved_value:.3f}"
            ref_str = f"{ref_value:.3f}"

        print(f"{param:<20} {solved_str:<15} {ref_str:<15} {rel_error:<12.2f} {status:<10}")

        test_results[param] = {
            'solved': solved_value,
            'reference': ref_value,
            'error': rel_error,
            'passed': passed
        }

    print(f"{'-'*80}")
    print(f"Overall Result: {'✓ ALL TESTS PASSED' if all_passed else '✗ SOME TESTS FAILED'}")
    print(f"{'='*80}\n")

    assert all_passed, f"{test_name} validation failed - see errors above"

    return test_results


# Updated solve_microstrip function with testing
def solve_microstrip():
    solver = MicrostripSolver2D(
        substrate_height=1.6e-3,
        trace_width=3e-3,
        trace_thickness=35e-6,
        gnd_thickness=35e-6,
        epsilon_r=4.5,
        tan_delta=0.02,
        sigma_cond=5.8e7,
        freq=1e9,
        nx=10, ny=10,
        use_sm=False,
        boundaries=["open", "open", "open", "gnd"]
    )
    #Z0, eps_eff, C, C0 = solver.calculate_parameters()
    #Ex, Ey = solver.compute_fields()
    Z0, eps_eff, C, C0 = solver.solve_adaptive()
    alpha_cond, J = solver.calculate_conductor_loss(solver.Ex, solver.Ey, Z0)
    alpha_diel = solver.calculate_dielectric_loss(solver.Ex, solver.Ey, Z0)
    alpha_total = alpha_cond + alpha_diel
    z, rlgc, eps_eff_mode = solver.rlgc(alpha_cond, alpha_diel, C, Z0)

    # Prepare results dictionary
    solver_results = {
        'Z0': Z0,
        'eps_eff': eps_eff,
        'diel_loss': alpha_diel,
        'cond_loss': alpha_cond,
        'C': rlgc['C'],
        'R': rlgc['R'],
        'L': rlgc['L'],
        'G': rlgc['G']
    }

    # Reference values
    reference = {
        "Z0": 49.8,
        "diel_loss": 2.99,
        "cond_loss": 0.285,
        "C": 123e-12,
        "R": 3.26,
        "G": 13.86e-3,
        "L": 307e-9
    }

    # Test against reference
    test_microstrip_solution(solver_results, reference, "Microstrip 50Ω")

    return solver_results

# Updated solve_differential_microstrip function with testing
def solve_differential_microstrip():
    solver = DifferentialMicrostripSolver2D(
        substrate_height=1.6e-3,
        trace_width=3e-3,
        trace_thickness=35e-6,
        gnd_thickness=16e-6,
        epsilon_r=4.5,
        freq=1e9,
        nx=10,
        ny=10,
        trace_spacing=1e-3  # This enables differential mode
    )

    # solve_adaptive now automatically detects differential mode
    results = solver.solve_adaptive()
    rlgc = solver.modal_to_physical_rlgc(results)

    reference = {
        'Z_odd': 40.23,
        'Z_even': 57.65,
        'eps_eff_even': 3.65,
        'eps_eff_odd': 2.98,
        'alpha_c_odd': 0.363,
        'alpha_c_even': 0.269,
        'alpha_d_odd': 2.67,
        'alpha_d_even': 3.21
    }

    # Test against reference
    test_differential_solution(results, reference, "Differential Microstrip")

    return results

def solve_gcpw():
    solver = GroundedCPWSolver2D(
        substrate_height=1.6e-3,
        trace_width=0.3e-3,
        trace_thickness=35e-6,
        gap=0.15e-3,
        top_gnd_width=5e-3,
        via_gap=0.5e-3,
        use_sm=False,
        sm_er=3.5,
        nx=10, ny=10,
    )

    #Z0, eps_eff, C, C0 = solver.calculate_parameters()
    #Ex, Ey = solver.compute_fields()
    Z0, eps_eff, C, C0 = solver.solve_adaptive(param_tol=0.001, max_iters=20, max_nodes=10000)

    alpha_cond, J = solver.calculate_conductor_loss(solver.Ex, solver.Ey, Z0)
    alpha_diel = solver.calculate_dielectric_loss(solver.Ex, solver.Ey, Z0)
    alpha_total = alpha_cond + alpha_diel
    z, rlgc, eps_eff_mode = solver.rlgc(alpha_cond, alpha_diel, C, Z0)

    solver_results = {
        'Z0': Z0,
        'eps_eff': eps_eff,
        'diel_loss': alpha_diel,
        'cond_loss': alpha_cond,
        'C': rlgc['C'],
        'R': rlgc['R'],
        'L': rlgc['L'],
        'G': rlgc['G']
    }

    # Reference values
    reference = {
        "Z0": 64.27,
        "loss": 4.288,
        "eps_eff": 2.5563
    }

    # Test against reference
    test_microstrip_solution(solver_results, reference, "GCPW")

def solve_gcpw_mask():
    solver = GroundedCPWSolver2D(
        substrate_height=1.6e-3,
        trace_width=0.3e-3,
        trace_thickness=35e-6,
        gap=0.15e-3,
        top_gnd_width=5e-3,
        via_gap=0.5e-3,
        use_sm=True,
        sm_er=3.5,
        nx=10, ny=10,
    )

    #Z0, eps_eff, C, C0 = solver.calculate_parameters()
    #Ex, Ey = solver.compute_fields()
    Z0, eps_eff, C, C0 = solver.solve_adaptive(param_tol=0.001, max_iters=20, max_nodes=10000)

    alpha_cond, J = solver.calculate_conductor_loss(solver.Ex, solver.Ey, Z0)
    alpha_diel = solver.calculate_dielectric_loss(solver.Ex, solver.Ey, Z0)
    alpha_total = alpha_cond + alpha_diel
    z, rlgc, eps_eff_mode = solver.rlgc(alpha_cond, alpha_diel, C, Z0)

    solver_results = {
        'Z0': Z0,
        'eps_eff': eps_eff,
        'diel_loss': alpha_diel,
        'cond_loss': alpha_cond,
        'C': rlgc['C'],
        'R': rlgc['R'],
        'L': rlgc['L'],
        'G': rlgc['G']
    }

    # Reference values
    reference = {
        "Z0": 60.52,
        "loss": 4.70,
        "eps_eff": 2.883,
        'C': 93.58e-12,
        #'R': 30.16,
        'G': 9.7e-3,
        'L': 342.7e-9
    }

    # Test against reference
    test_microstrip_solution(solver_results, reference, "GCPW mask")

def solve_microstrip_embed():
    solver = MicrostripSolver2D(
        substrate_height=1.6e-3,
        trace_width=3e-3,
        trace_thickness=35e-6,
        gnd_thickness=35e-6,
        epsilon_r=4.5,
        tan_delta=0.02,
        sigma_cond=5.8e7,
        freq=1e9,
        nx=10, ny=10,
        use_sm=False,
        top_diel_h=0.2e-3,
        top_diel_er=4.5,
        boundaries=["open", "open", "open", "gnd"]
    )

    #Z0, eps_eff, C, C0 = solver.calculate_parameters()
    #Ex, Ey = solver.compute_fields()
    Z0, eps_eff, C, C0 = solver.solve_adaptive()

    alpha_cond, J = solver.calculate_conductor_loss(solver.Ex, solver.Ey, Z0)
    alpha_diel = solver.calculate_dielectric_loss(solver.Ex, solver.Ey, Z0)
    alpha_total = alpha_cond + alpha_diel

    z, rlgc, eps_eff_mode = solver.rlgc(alpha_cond, alpha_diel, C, Z0)

    solver_results = {
        'Z0': Z0,
        'eps_eff': eps_eff,
        'diel_loss': alpha_diel,
        'cond_loss': alpha_cond,
        'loss': alpha_total,
        'C': rlgc['C'],
        'R': rlgc['R'],
        'L': rlgc['L'],
        'G': rlgc['G']
    }

    # Reference values
    reference = {
        "Z0": 48.15,
        "eps_eff": 3.621,
        "loss": 3.48,
        "C": 131.8e-12
    }

    # Test against reference
    test_microstrip_solution(solver_results, reference, "Embedded microstrip")

def solve_microstrip_cut():
    solver = MicrostripSolver2D(
        substrate_height=1.6e-3,
        trace_width=3e-3,
        trace_thickness=35e-6,
        gnd_thickness=35e-6,
        epsilon_r=4.5,
        tan_delta=0.02,
        sigma_cond=5.8e7,
        freq=1e9,
        nx=10, ny=10,
        use_sm=False,
        gnd_cut_width=3e-3,
        gnd_cut_sub_h=1e-3,
        #top_diel_h=0.2e-3,
        #top_diel_er=4.5,
        boundaries=["open", "open", "open", "gnd"]
    )

    #Z0, eps_eff, C, C0 = solver.calculate_parameters()
    #Ex, Ey = solver.compute_fields()
    Z0, eps_eff, C, C0 = solver.solve_adaptive()

    alpha_cond, J = solver.calculate_conductor_loss(solver.Ex, solver.Ey, Z0)
    alpha_diel = solver.calculate_dielectric_loss(solver.Ex, solver.Ey, Z0)
    alpha_total = alpha_cond + alpha_diel

    z, rlgc, eps_eff_mode = solver.rlgc(alpha_cond, alpha_diel, C, Z0)

    solver_results = {
        'Z0': Z0,
        'eps_eff': eps_eff,
        'diel_loss': alpha_diel,
        'cond_loss': alpha_cond,
        'loss': alpha_total,
        'C': rlgc['C'],
        'R': rlgc['R'],
        'L': rlgc['L'],
        'G': rlgc['G']
    }

    # Reference values
    reference = {
        "Z0": 55.84,
        "eps_eff": 3.28,
        "loss": 3.19,
        "C": 108.25e-12
    }

    # Test against reference
    test_microstrip_solution(solver_results, reference, "Microstrip with ground cut")

def solve_differential_stripline():
    solver = DifferentialMicrostripSolver2D(
            substrate_height=0.2e-3,
            trace_width=0.15e-3,
            trace_thickness=35e-6,
            gnd_thickness=16e-6,
            epsilon_r=4.1,
            epsilon_r_top=4.1,
            air_top=0.2e-3,
            freq=1e9,
            nx=10,
            ny=10,
            trace_spacing=0.1e-3,  # This enables differential mode
            boundaries=["open", "open", "gnd", "gnd"]
        )

    # solve_adaptive now automatically detects differential mode
    results = solver.solve_adaptive()
    rlgc = solver.modal_to_physical_rlgc(results)

    reference = {
        'Z_odd': 37.6,
        'Z_even': 61.36,
        'eps_eff_even': 4.162,
        'eps_eff_odd': 4.195,
        'alpha_total_odd': 7.93,
        'alpha_total_even': 6.47,
    }

    # Test against reference
    test_differential_solution(results, reference, "Differential Stripline")

    return results

if __name__ == "__main__":
    solve_microstrip()
    solve_microstrip_embed()
    solve_microstrip_cut()
    solve_differential_stripline()
    solve_differential_microstrip()
    solve_gcpw()
    solve_gcpw_mask()
