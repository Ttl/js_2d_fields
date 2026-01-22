import { MicrostripSolver } from './microstrip.js';

/**
 * Test microstrip solver results against reference values.
 *
 * @param {Object} solver_results - Object with keys 'Z0', 'eps_eff', 'alpha_diel', 'alpha_cond', 'C', 'R', 'L', 'G'
 * @param {Object} reference - Object with reference values
 * @param {string} test_name - Name of the test for printing
 * @returns {Object} Test results with pass/fail status
 */
function test_microstrip_solution(solver_results, reference, test_name = "Microstrip") {
    // Global error thresholds (relative error in %)
    const MAX_Z0_ERROR = 5.0;
    const MAX_DIEL_LOSS_ERROR = 10.0;
    const MAX_COND_LOSS_ERROR = 50.0;
    const MAX_C_ERROR = 10.0;
    const MAX_R_ERROR = 15.0;
    const MAX_L_ERROR = 10.0;
    const MAX_G_ERROR = 15.0;
    const MAX_EPS_EFF_ERROR = 5.0;

    // Error thresholds mapping
    const error_thresholds = {
        'Z0': MAX_Z0_ERROR,
        'diel_loss': MAX_DIEL_LOSS_ERROR,
        'cond_loss': MAX_COND_LOSS_ERROR,
        'C': MAX_C_ERROR,
        'R': MAX_R_ERROR,
        'L': MAX_L_ERROR,
        'G': MAX_G_ERROR,
        'eps_eff': MAX_EPS_EFF_ERROR
    };

    console.log(`\n${'='.repeat(80)}`);
    console.log(`${test_name.toUpperCase()} VALIDATION TEST`);
    console.log(`${'='.repeat(80)}`);
    console.log(`${'Parameter'.padEnd(15)} ${'Solved'.padEnd(15)} ${'Reference'.padEnd(15)} ${'Error (%)'.padEnd(12)} ${'Status'.padEnd(10)}`);
    console.log(`${'-'.repeat(80)}`);

    let all_passed = true;
    const test_results = {};

    for (const [param, ref_value] of Object.entries(reference)) {
        if (!(param in solver_results)) {
            continue;
        }

        const solved_value = solver_results[param];

        // Calculate relative error
        let rel_error;
        if (ref_value !== 0) {
            rel_error = Math.abs((solved_value - ref_value) / ref_value) * 100;
        } else {
            rel_error = Math.abs(solved_value) * 100;
        }

        // Check against threshold
        const threshold = error_thresholds[param] || 10.0; // default 10%
        const passed = rel_error <= threshold;
        all_passed = all_passed && passed;

        const status = passed ? "✓ PASS" : "✗ FAIL";

        // Format values based on magnitude
        let solved_str, ref_str;
        if (param === 'C') {
            solved_str = `${(solved_value * 1e12).toFixed(2)} pF`;
            ref_str = `${(ref_value * 1e12).toFixed(2)} pF`;
        } else if (param === 'L') {
            solved_str = `${(solved_value * 1e9).toFixed(2)} nH`;
            ref_str = `${(ref_value * 1e9).toFixed(2)} nH`;
        } else if (param === 'G') {
            solved_str = `${(solved_value * 1e3).toFixed(2)} mS`;
            ref_str = `${(ref_value * 1e3).toFixed(2)} mS`;
        } else if (param === 'R') {
            solved_str = `${solved_value.toFixed(2)} Ω`;
            ref_str = `${ref_value.toFixed(2)} Ω`;
        } else if (param.includes('loss')) {
            solved_str = `${solved_value.toFixed(4)} dB/m`;
            ref_str = `${ref_value.toFixed(4)} dB/m`;
        } else {
            solved_str = `${solved_value.toFixed(3)}`;
            ref_str = `${ref_value.toFixed(3)}`;
        }

        console.log(`${param.padEnd(15)} ${solved_str.padEnd(15)} ${ref_str.padEnd(15)} ${rel_error.toFixed(2).padEnd(12)} ${status.padEnd(10)}`);

        test_results[param] = {
            'solved': solved_value,
            'reference': ref_value,
            'error': rel_error,
            'passed': passed
        };
    }

    console.log(`${'-'.repeat(80)}`);
    console.log(`Overall Result: ${all_passed ? '✓ ALL TESTS PASSED' : '✗ SOME TESTS FAILED'}`);
    console.log(`${'='.repeat(80)}\n`);

    if (!all_passed) {
        throw new Error(`${test_name} validation failed - see errors above`);
    }

    return test_results;
}

function test_differential_solution(solver_results, reference, test_name = "Differential Microstrip") {
    // Global error thresholds (relative error in %)
    const MAX_Z_DIFF_ERROR = 6.0;
    const MAX_Z_COMMON_ERROR = 6.0;
    const MAX_Z_ODD_ERROR = 6.0;
    const MAX_Z_EVEN_ERROR = 6.0;
    const MAX_EPS_EFF_ERROR = 5.0;
    const MAX_LOSS_ERROR = 50.0;
    const MAX_C_ERROR = 10.0;

    // Error thresholds mapping
    const error_thresholds = {
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
    };

    console.log(`\n${'='.repeat(80)}`);
    console.log(`${test_name.toUpperCase()} VALIDATION TEST`);
    console.log(`${'='.repeat(80)}`);
    console.log(`${'Parameter'.padEnd(20)} ${'Solved'.padEnd(15)} ${'Reference'.padEnd(15)} ${'Error (%)'.padEnd(12)} ${'Status'.padEnd(10)}`);
    console.log(`${'-'.repeat(80)}`);

    let all_passed = true;
    const test_results = {};

    for (const [param, ref_value] of Object.entries(reference)) {
        if (!(param in solver_results)) {
            continue;
        }

        const solved_value = solver_results[param];

        // Calculate relative error
        let rel_error;
        if (ref_value !== 0) {
            rel_error = Math.abs((solved_value - ref_value) / ref_value) * 100;
        } else {
            rel_error = Math.abs(solved_value) * 100;
        }

        // Check against threshold
        const threshold = error_thresholds[param] || 10.0;
        const passed = rel_error <= threshold;
        all_passed = all_passed && passed;

        const status = passed ? "✓ PASS" : "✗ FAIL";

        // Format values based on parameter type
        let solved_str, ref_str;
        if (param.startsWith('C_')) {
            solved_str = `${(solved_value * 1e12).toFixed(2)} pF`;
            ref_str = `${(ref_value * 1e12).toFixed(2)} pF`;
        } else if (param.startsWith('alpha_')) {
            solved_str = `${solved_value.toFixed(4)} dB/m`;
            ref_str = `${ref_value.toFixed(4)} dB/m`;
        } else if (param.startsWith('Z_')) {
            solved_str = `${solved_value.toFixed(2)} Ω`;
            ref_str = `${ref_value.toFixed(2)} Ω`;
        } else {
            solved_str = `${solved_value.toFixed(3)}`;
            ref_str = `${ref_value.toFixed(3)}`;
        }

        console.log(`${param.padEnd(20)} ${solved_str.padEnd(15)} ${ref_str.padEnd(15)} ${rel_error.toFixed(2).padEnd(12)} ${status.padEnd(10)}`);

        test_results[param] = {
            'solved': solved_value,
            'reference': ref_value,
            'error': rel_error,
            'passed': passed
        };
    }

    console.log(`${'-'.repeat(80)}`);
    console.log(`Overall Result: ${all_passed ? '✓ ALL TESTS PASSED' : '✗ SOME TESTS FAILED'}`);
    console.log(`${'='.repeat(80)}\n`);

    if (!all_passed) {
        throw new Error(`${test_name} validation failed - see errors above`);
    }

    return test_results;
}

async function solve_microstrip() {
    const solver = new MicrostripSolver({
        substrate_height: 1.6e-3,
        trace_width: 3e-3,
        trace_thickness: 35e-6,
        gnd_thickness: 35e-6,
        epsilon_r: 4.5,
        tan_delta: 0.02,
        sigma_cond: 5.8e7,
        freq: 1e9,
        nx: 10,
        ny: 10,
        use_sm: false,
        boundaries: ["open", "open", "open", "gnd"]
    });

    const results = await solver.solve_adaptive();
    const mode = results.modes[0];

    // Prepare results dictionary matching Python format
    const solver_results = {
        'Z0': mode.Z0,
        'eps_eff': mode.eps_eff,
        'diel_loss': mode.alpha_d,
        'cond_loss': mode.alpha_c,
        'C': mode.RLGC.C,
        'R': mode.RLGC.R,
        'L': mode.RLGC.L,
        'G': mode.RLGC.G
    };

    // Reference values from HFSS
    const reference = {
        "Z0": 49.8,
        "diel_loss": 2.99,
        "cond_loss": 0.285,
        "C": 123e-12,
        "R": 3.26,
        "G": 13.86e-3,
        "L": 307e-9
    };

    // Test against reference
    test_microstrip_solution(solver_results, reference, "Microstrip 50Ω");

    return solver_results;
}

async function solve_microstrip_embed() {
    const solver = new MicrostripSolver({
        substrate_height: 1.6e-3,
        trace_width: 3e-3,
        trace_thickness: 35e-6,
        gnd_thickness: 35e-6,
        epsilon_r: 4.5,
        tan_delta: 0.02,
        sigma_cond: 5.8e7,
        freq: 1e9,
        nx: 10,
        ny: 10,
        use_sm: false,
        top_diel_h: 0.2e-3,
        top_diel_er: 4.5,
        boundaries: ["open", "open", "open", "gnd"]
    });

    const results = await solver.solve_adaptive();
    const mode = results.modes[0];

    const solver_results = {
        'Z0': mode.Z0,
        'eps_eff': mode.eps_eff,
        'diel_loss': mode.alpha_d,
        'cond_loss': mode.alpha_c,
        'loss': mode.alpha_total,
        'C': mode.RLGC.C,
        'R': mode.RLGC.R,
        'L': mode.RLGC.L,
        'G': mode.RLGC.G
    };

    // Reference values from HFSS
    const reference = {
        "Z0": 48.15,
        "eps_eff": 3.621,
        "loss": 3.48,
        "C": 131.8e-12
    };

    // Test against reference
    test_microstrip_solution(solver_results, reference, "Embedded microstrip");

    return solver_results;
}

async function solve_microstrip_cut() {
    const solver = new MicrostripSolver({
        substrate_height: 1.6e-3,
        trace_width: 3e-3,
        trace_thickness: 35e-6,
        gnd_thickness: 35e-6,
        epsilon_r: 4.5,
        tan_delta: 0.02,
        sigma_cond: 5.8e7,
        freq: 1e9,
        nx: 10,
        ny: 10,
        use_sm: false,
        gnd_cut_width: 3e-3,
        gnd_cut_sub_h: 1e-3,
        boundaries: ["open", "open", "open", "gnd"]
    });

    const results = await solver.solve_adaptive();
    const mode = results.modes[0];

    const solver_results = {
        'Z0': mode.Z0,
        'eps_eff': mode.eps_eff,
        'diel_loss': mode.alpha_d,
        'cond_loss': mode.alpha_c,
        'loss': mode.alpha_total,
        'C': mode.RLGC.C,
        'R': mode.RLGC.R,
        'L': mode.RLGC.L,
        'G': mode.RLGC.G
    };

    // Reference values from HFSS
    const reference = {
        "Z0": 55.84,
        "eps_eff": 3.28,
        "loss": 3.19,
        "C": 108.25e-12
    };

    // Test against reference
    test_microstrip_solution(solver_results, reference, "Microstrip with ground cut");

    return solver_results;
}

async function solve_differential_microstrip() {
    const solver = new MicrostripSolver({
        substrate_height: 1.6e-3,
        trace_width: 3e-3,
        trace_thickness: 35e-6,
        gnd_thickness: 16e-6,
        epsilon_r: 4.5,
        freq: 1e9,
        nx: 10,
        ny: 10,
        trace_spacing: 1e-3  // This enables differential mode
    });

    const results = await solver.solve_adaptive({energy_tol: 0.01});
    const odd = results.modes.find(m => m.mode === 'odd');
    const even = results.modes.find(m => m.mode === 'even');

    // Map to flat structure for test function
    const solver_results = {
        'Z_odd': odd.Z0,
        'Z_even': even.Z0,
        'Z_diff': results.Z_diff,
        'Z_common': results.Z_common,
        'eps_eff_odd': odd.eps_eff,
        'eps_eff_even': even.eps_eff,
        'alpha_c_odd': odd.alpha_c,
        'alpha_c_even': even.alpha_c,
        'alpha_d_odd': odd.alpha_d,
        'alpha_d_even': even.alpha_d,
        'alpha_total_odd': odd.alpha_total,
        'alpha_total_even': even.alpha_total
    };

    const reference = {
        'Z_odd': 40.23,
        'Z_even': 57.65,
        'eps_eff_even': 3.65,
        'eps_eff_odd': 2.98,
        'alpha_c_odd': 0.363,
        'alpha_c_even': 0.269,
        'alpha_d_odd': 2.67,
        'alpha_d_even': 3.21
    };

    // Test against reference
    test_differential_solution(solver_results, reference, "Differential Microstrip");

    return results;
}

async function solve_differential_stripline() {
    const solver = new MicrostripSolver({
        substrate_height: 0.2e-3,
        trace_width: 0.15e-3,
        trace_thickness: 35e-6,
        gnd_thickness: 16e-6,
        epsilon_r: 4.1,
        epsilon_r_top: 4.1,
        enclosure_height: 0.2e-3 + 35e-6,
        freq: 1e9,
        nx: 10,
        ny: 10,
        trace_spacing: 0.1e-3,  // This enables differential mode
        boundaries: ["open", "open", "gnd", "gnd"]
    });

    const results = await solver.solve_adaptive();
    const odd = results.modes.find(m => m.mode === 'odd');
    const even = results.modes.find(m => m.mode === 'even');

    // Map to flat structure for test function
    const solver_results = {
        'Z_odd': odd.Z0,
        'Z_even': even.Z0,
        'Z_diff': results.Z_diff,
        'Z_common': results.Z_common,
        'eps_eff_odd': odd.eps_eff,
        'eps_eff_even': even.eps_eff,
        'alpha_c_odd': odd.alpha_c,
        'alpha_c_even': even.alpha_c,
        'alpha_d_odd': odd.alpha_d,
        'alpha_d_even': even.alpha_d,
        'alpha_total_odd': odd.alpha_total,
        'alpha_total_even': even.alpha_total
    };

    const reference = {
        'Z_odd': 37.6,
        'Z_even': 61.36,
        'eps_eff_even': 4.162,
        'eps_eff_odd': 4.195,
        'alpha_total_odd': 7.93,
        'alpha_total_even': 6.47,
    };

    // Test against reference
    test_differential_solution(solver_results, reference, "Differential Stripline");

    return results;
}

// Run tests
async function runTests() {
    await solve_microstrip();
    await solve_microstrip_embed();
    await solve_microstrip_cut();
    await solve_differential_stripline();
    await solve_differential_microstrip();
}

runTests();
