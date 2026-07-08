import { MicrostripSolver } from '../src/microstrip.js';

/**
 * Test: Verify that poor side plating (lower sigma than bulk) correctly
 * increases loss at bottom corners within plating thickness distance.
 *
 * This validates the intermediate case where side plating is worse than bulk.
 */

async function test_poor_side_plating_at_corners() {
    console.log('Testing: Poor side plating at bottom corners');
    console.log('Expected: Loss increases when side sigma < bulk sigma\n');

    const freq = 1e9;
    const copper_sigma = 5.96e7;  // Good bulk conductor
    const poor_sigma = 1e7;       // Poor plating (e.g., solder, tin)
    const plating_thickness = 10e-6;  // 10 μm

    console.log(`Bulk conductivity (copper): ${(copper_sigma / 1e7).toFixed(2)} × 10^7 S/m`);
    console.log(`Side plating conductivity: ${(poor_sigma / 1e7).toFixed(2)} × 10^7 S/m`);
    console.log(`Plating thickness: ${(plating_thickness * 1e6).toFixed(1)} μm\n`);

    // Case 1: No plating (baseline)
    const solver1 = new MicrostripSolver({
        substrate_height: 0.8e-3,
        trace_width: 100e-6,
        trace_thickness: 35e-6,
        gnd_thickness: 35e-6,
        epsilon_r: 4.5,
        tan_delta: 0.02,
        sigma_cond: copper_sigma,
        freq: freq,
        nx: 50,
        ny: 25,
        plating: null,  // No plating
        boundaries: ["open", "open", "open", "gnd"],
        rq: 0  // No roughness
    });

    // Case 2: Poor side plating (should increase loss at corners)
    const poor_plating = {
        sigma: poor_sigma,  // Much worse than copper
        thickness: plating_thickness,
        rq: 0,  // No roughness
        top: false,
        sides: true,  // Poor plating on sides
        bottom: false,
        thick_corners: true  // Model thick plating covering bottom near edges
    };

    const solver2 = new MicrostripSolver({
        substrate_height: 0.8e-3,
        trace_width: 100e-6,
        trace_thickness: 35e-6,
        gnd_thickness: 35e-6,
        epsilon_r: 4.5,
        tan_delta: 0.02,
        sigma_cond: copper_sigma,
        freq: freq,
        nx: 50,
        ny: 25,
        plating: poor_plating,
        boundaries: ["open", "open", "open", "gnd"],
        rq: 0
    });

    // Case 3: Same poor conductivity as solid bulk (reference)
    const solver3 = new MicrostripSolver({
        substrate_height: 0.8e-3,
        trace_width: 100e-6,
        trace_thickness: 35e-6,
        gnd_thickness: 35e-6,
        epsilon_r: 4.5,
        tan_delta: 0.02,
        sigma_cond: poor_sigma,  // Poor bulk, no plating
        freq: freq,
        nx: 50,
        ny: 25,
        plating: null,
        boundaries: ["open", "open", "open", "gnd"],
        rq: 0
    });

    const result1 = await solver1.solve_adaptive({ max_refinements: 0 });
    const mode1 = result1.modes[0];

    const result2 = await solver2.solve_adaptive({ max_refinements: 0 });
    const mode2 = result2.modes[0];

    const result3 = await solver3.solve_adaptive({ max_refinements: 0 });
    const mode3 = result3.modes[0];

    console.log('Case 1: No plating (good copper bulk)');
    console.log(`  Conductor loss = ${mode1.alpha_c.toFixed(6)} dB/m`);

    console.log('\nCase 2: Poor side plating on good copper bulk');
    console.log(`  Conductor loss = ${mode2.alpha_c.toFixed(6)} dB/m`);

    console.log('\nCase 3: Poor bulk conductor (no plating)');
    console.log(`  Conductor loss = ${mode3.alpha_c.toFixed(6)} dB/m`);

    // Analysis
    console.log('\n' + '='.repeat(60));
    console.log('ANALYSIS:');

    const loss_increase = mode2.alpha_c - mode1.alpha_c;
    const loss_increase_pct = (loss_increase / mode1.alpha_c) * 100;

    console.log(`\nLoss increase from poor side plating: ${loss_increase.toFixed(6)} dB/m`);
    console.log(`Percentage increase: ${loss_increase_pct.toFixed(2)}%`);

    if (loss_increase > 0) {
        console.log('\n✓ CORRECT: Poor side plating increases loss at corners');
        console.log('  Within plating thickness from corners, side plating impedance is used.');
    } else {
        console.log('\n✗ ERROR: Poor side plating should increase loss!');
        console.log('  The corner model may not be applying side plating impedance correctly.');
        return false;
    }

    // Check corner coverage
    const trace_width = 100e-6;
    const coverage_fraction = (2 * plating_thickness) / trace_width;
    console.log(`\nGeometric coverage of bottom surface by side plating:`);
    console.log(`  Coverage width: 2 × ${(plating_thickness * 1e6).toFixed(1)} μm = ${(2 * plating_thickness * 1e6).toFixed(1)} μm`);
    console.log(`  Trace width: ${(trace_width * 1e6).toFixed(1)} μm`);
    console.log(`  Coverage fraction: ${(coverage_fraction * 100).toFixed(1)}%`);

    // Expected behavior: loss should be between baseline and full poor conductor
    const expected_range_min = mode1.alpha_c;
    const expected_range_max = mode3.alpha_c;

    console.log(`\nExpected loss range:`);
    console.log(`  Minimum (no plating): ${expected_range_min.toFixed(6)} dB/m`);
    console.log(`  Maximum (full poor conductor): ${expected_range_max.toFixed(6)} dB/m`);
    console.log(`  Actual (poor side plating): ${mode2.alpha_c.toFixed(6)} dB/m`);

    if (mode2.alpha_c >= expected_range_min && mode2.alpha_c <= expected_range_max) {
        console.log('\n✓ CORRECT: Loss is within expected range');
        console.log('  Poor side plating partially degrades conductor, as expected.');
    } else {
        console.log('\n✗ ERROR: Loss is outside expected range!');
        return false;
    }

    // Verify corner distance rule
    console.log('\n' + '='.repeat(60));
    console.log('PHYSICAL MODEL VERIFICATION:');
    console.log(`\nAt bottom surface position x from left edge:`);
    console.log(`  If x < ${(plating_thickness * 1e6).toFixed(1)} μm: covered by left side plating → use poor sigma`);
    console.log(`  If x > ${((trace_width - plating_thickness) * 1e6).toFixed(1)} μm: covered by right side plating → use poor sigma`);
    console.log(`  Otherwise: exposed copper bulk → use good sigma`);

    return true;
}

// Test with varying plating thickness
async function test_varying_plating_thickness() {
    console.log('\n\n' + '='.repeat(60));
    console.log('Testing: Varying plating thickness');
    console.log('Expected: Loss increases as plating thickness increases\n');

    const freq = 1e9;
    const copper_sigma = 5.96e7;
    const poor_sigma = 1e7;
    const trace_width = 100e-6;

    const thicknesses = [
        { t: 5e-6, label: '5 μm (5% coverage)' },
        { t: 10e-6, label: '10 μm (10% coverage)' },
        { t: 25e-6, label: '25 μm (25% coverage)' },
        { t: 50e-6, label: '50 μm (50% coverage)' }
    ];

    const results = [];

    for (const { t, label } of thicknesses) {
        const plating = {
            sigma: poor_sigma,
            thickness: t,
            rq: 0,
            top: false,
            sides: true,
            bottom: false,
            thick_corners: true  // Model thick plating covering bottom near edges
        };

        const solver = new MicrostripSolver({
            substrate_height: 0.8e-3,
            trace_width: trace_width,
            trace_thickness: 35e-6,
            gnd_thickness: 35e-6,
            epsilon_r: 4.5,
            tan_delta: 0.02,
            sigma_cond: copper_sigma,
            freq: freq,
            nx: 50,
            ny: 25,
            plating: plating,
            boundaries: ["open", "open", "open", "gnd"],
            rq: 0
        });

        const result = await solver.solve_adaptive({ max_refinements: 0 });
        const mode = result.modes[0];

        results.push({ t, label, alpha: mode.alpha_c });
        console.log(`  ${label}: α = ${mode.alpha_c.toFixed(6)} dB/m`);
    }

    // Verify monotonic increase
    console.log('\n' + '='.repeat(60));
    let monotonic = true;
    for (let i = 1; i < results.length; i++) {
        if (results[i].alpha < results[i-1].alpha) {
            monotonic = false;
            console.log(`✗ ERROR: Loss decreased from ${results[i-1].label} to ${results[i].label}`);
        }
    }

    if (monotonic) {
        console.log('✓ CORRECT: Loss increases monotonically with plating thickness');
        console.log('  As expected: thicker poor plating covers more of bottom surface.');
        return true;
    } else {
        console.log('✗ ERROR: Loss should increase monotonically with plating thickness');
        return false;
    }
}

// Run tests
async function run_all_tests() {
    let test1_pass = false;
    let test2_pass = false;

    try {
        console.log('='.repeat(60));
        console.log('POOR SIDE PLATING TEST SUITE');
        console.log('='.repeat(60) + '\n');

        test1_pass = await test_poor_side_plating_at_corners();
        test2_pass = await test_varying_plating_thickness();

        console.log('\n' + '='.repeat(60));
        if (test1_pass && test2_pass) {
            console.log('✓ ALL TESTS PASSED');
            console.log('='.repeat(60));
            process.exit(0);
        } else {
            console.log('✗ SOME TESTS FAILED');
            console.log('='.repeat(60));
            process.exit(1);
        }
    } catch (err) {
        console.error('\n✗ Test error:', err);
        process.exit(1);
    }
}

run_all_tests();
