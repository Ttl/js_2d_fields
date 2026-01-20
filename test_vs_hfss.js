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
        skin_cells: 50,
        use_sm: false,
        boundaries: ["open", "open", "open", "gnd"]
    });

    const results = await solver.solve_adaptive();

    // Prepare results dictionary matching Python format
    const solver_results = {
        'Z0': results.Z0,
        'eps_eff': results.eps_eff,
        'diel_loss': results.alpha_diel_db_m,
        'cond_loss': results.alpha_cond_db_m,
        'C': results.RLGC.C,
        'R': results.RLGC.R,
        'L': results.RLGC.L,
        'G': results.RLGC.G
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
        skin_cells: 50,
        use_sm: false,
        top_diel_h: 0.2e-3,
        top_diel_er: 4.5,
        boundaries: ["open", "open", "open", "gnd"]
    });

    const results = await solver.solve_adaptive();

    const solver_results = {
        'Z0': results.Z0,
        'eps_eff': results.eps_eff,
        'diel_loss': results.alpha_diel_db_m,
        'cond_loss': results.alpha_cond_db_m,
        'loss': results.total_alpha_db_m,
        'C': results.RLGC.C,
        'R': results.RLGC.R,
        'L': results.RLGC.L,
        'G': results.RLGC.G
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
        skin_cells: 50,
        use_sm: false,
        gnd_cut_width: 3e-3,
        gnd_cut_sub_h: 1e-3,
        boundaries: ["open", "open", "open", "gnd"]
    });

    const results = await solver.solve_adaptive();

    const solver_results = {
        'Z0': results.Z0,
        'eps_eff': results.eps_eff,
        'diel_loss': results.alpha_diel_db_m,
        'cond_loss': results.alpha_cond_db_m,
        'loss': results.total_alpha_db_m,
        'C': results.RLGC.C,
        'R': results.RLGC.R,
        'L': results.RLGC.L,
        'G': results.RLGC.G
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

// Run tests
async function runTests() {
    await solve_microstrip();
    await solve_microstrip_embed();
    await solve_microstrip_cut();
}

runTests();
