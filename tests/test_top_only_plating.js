import { MicrostripSolver } from '../src/microstrip.js';

/**
 * Test: Verify that top-only plating correctly extends down the sides
 *
 * When only top plating is enabled:
 * - Top surface: layered plating impedance
 * - Sides: from top to (top - plating.thickness): single-layer plating with plating.rq
 * - Bottom: no plating (bulk impedance)
 */

async function test_top_only_plating() {
    console.log('='.repeat(60));
    console.log('TOP-ONLY PLATING TEST');
    console.log('='.repeat(60));
    console.log('\nPhysical model:');
    console.log('  When only TOP plating enabled:');
    console.log('  - Top surface: layered plating');
    console.log('  - Sides (near top): single-layer plating with plating.rq');
    console.log('  - Bottom: no plating (bulk)\n');

    const freq = 1e9;
    const copper_sigma = 5.96e7;
    const silver_sigma = 6.3e7;  // Silver - better than copper
    const plating_thickness = 10e-6;  // 10 μm

    console.log('Configuration:');
    console.log(`  Frequency: ${freq / 1e9} GHz`);
    console.log(`  Bulk: Copper (${(copper_sigma / 1e7).toFixed(2)} × 10^7 S/m)`);
    console.log(`  Plating: Silver (${(silver_sigma / 1e7).toFixed(2)} × 10^7 S/m - BETTER)`);
    console.log(`  Plating thickness: ${(plating_thickness * 1e6).toFixed(0)} μm`);
    console.log(`  Trace height: 35 μm`);

    // Case 1: Only top plating
    const top_only_plating = {
        sigma: silver_sigma,
        thickness: plating_thickness,
        rq: 0,  // Smooth plating to isolate conductivity effect
        top: true,      // TOP enabled
        sides: false,   // Sides disabled
        bottom: false   // Bottom disabled
    };

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
        plating: top_only_plating,
        boundaries: ["open", "open", "open", "gnd"],
        rq: 0  // Smooth bulk
    });

    // Case 2: No plating (baseline)
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
        plating: null,
        boundaries: ["open", "open", "open", "gnd"],
        rq: 0
    });

    // Case 3: Full plating (all surfaces)
    const full_plating = {
        sigma: silver_sigma,
        thickness: plating_thickness,
        rq: 0,  // Smooth
        top: true,
        sides: true,
        bottom: true
    };

    const solver3 = new MicrostripSolver({
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
        plating: full_plating,
        boundaries: ["open", "open", "open", "gnd"],
        rq: 0
    });

    const result1 = await solver1.solve_adaptive({ max_refinements: 0 });
    const mode1 = result1.modes[0];

    const result2 = await solver2.solve_adaptive({ max_refinements: 0 });
    const mode2 = result2.modes[0];

    const result3 = await solver3.solve_adaptive({ max_refinements: 0 });
    const mode3 = result3.modes[0];

    console.log('\n' + '='.repeat(60));
    console.log('RESULTS:\n');

    console.log('Case 1: Top-only plating');
    console.log(`  Z0 = ${mode1.Z0.toFixed(4)} Ω`);
    console.log(`  Conductor loss = ${mode1.alpha_c.toFixed(6)} dB/m`);

    console.log('\nCase 2: No plating (copper baseline)');
    console.log(`  Z0 = ${mode2.Z0.toFixed(4)} Ω`);
    console.log(`  Conductor loss = ${mode2.alpha_c.toFixed(6)} dB/m`);

    console.log('\nCase 3: Full plating (all surfaces)');
    console.log(`  Z0 = ${mode3.Z0.toFixed(4)} Ω`);
    console.log(`  Conductor loss = ${mode3.alpha_c.toFixed(6)} dB/m`);

    console.log('\n' + '='.repeat(60));
    console.log('ANALYSIS:\n');

    const loss_reduction_top_only = mode2.alpha_c - mode1.alpha_c;
    const loss_reduction_full = mode2.alpha_c - mode3.alpha_c;

    console.log(`Loss reduction from top-only plating: ${loss_reduction_top_only.toFixed(6)} dB/m`);
    console.log(`Loss reduction from full plating: ${loss_reduction_full.toFixed(6)} dB/m`);

    const coverage_ratio = (plating_thickness / 35e-6) * 100;
    console.log(`\nTop plating coverage on sides: ${coverage_ratio.toFixed(1)}%`);
    console.log(`  (${(plating_thickness * 1e6).toFixed(0)} μm out of 35 μm trace height)`);

    if (loss_reduction_top_only > 0) {
        console.log('\n✓ CORRECT: Top plating extends down sides and reduces loss');
    } else {
        console.log('\n✗ ERROR: Top plating should extend down sides!');
        return false;
    }

    if (mode1.alpha_c > mode3.alpha_c) {
        console.log('✓ CORRECT: Top-only has higher loss than full plating');
        console.log('  (Bottom and lower sides are unplated copper)');
    } else {
        console.log('✗ ERROR: Top-only should have higher loss than full plating');
        return false;
    }

    if (mode1.alpha_c < mode2.alpha_c) {
        console.log('✓ CORRECT: Top-only has lower loss than no plating');
        console.log('  (Top and upper sides are plated gold)');
    } else {
        console.log('✗ ERROR: Top-only should have lower loss than no plating');
        return false;
    }

    console.log('\n' + '='.repeat(60));
    console.log('EXPECTED BEHAVIOR:\n');

    console.log('Top surface (y = y_max):');
    console.log('  → Layered plating impedance');

    console.log('\nSide surfaces:');
    console.log(`  → From y_max down to (y_max - ${(plating_thickness * 1e6).toFixed(0)} μm):`);
    console.log('    • Single-layer plating with plating.rq');
    console.log(`  → Below (y_max - ${(plating_thickness * 1e6).toFixed(0)} μm):`);
    console.log('    • Bulk copper (no plating)');

    console.log('\nBottom surface (y = y_min):');
    console.log('  → Bulk copper (no plating)');

    console.log('\n' + '='.repeat(60));
    return true;
}

// Test with varying plating thickness
async function test_top_plating_thickness_variation() {
    console.log('\n' + '='.repeat(60));
    console.log('TOP PLATING THICKNESS VARIATION');
    console.log('='.repeat(60));
    console.log('\nExpected: Loss should decrease as top plating thickness increases\n');

    const freq = 1e9;
    const copper_sigma = 5.96e7;
    const silver_sigma = 6.3e7;  // Silver - better than copper
    const trace_height = 35e-6;

    const thicknesses = [
        { t: 5e-6, label: '5 μm (14% coverage)', coverage: 5 / 35 },
        { t: 10e-6, label: '10 μm (29% coverage)', coverage: 10 / 35 },
        { t: 20e-6, label: '20 μm (57% coverage)', coverage: 20 / 35 },
        { t: 35e-6, label: '35 μm (100% coverage)', coverage: 35 / 35 }
    ];

    const results = [];

    for (const { t, label } of thicknesses) {
        const plating = {
            sigma: silver_sigma,
            thickness: t,
            rq: 0,  // Smooth
            top: true,
            sides: false,  // Only top
            bottom: false
        };

        const solver = new MicrostripSolver({
            substrate_height: 0.8e-3,
            trace_width: 100e-6,
            trace_thickness: trace_height,
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

    // Verify monotonic decrease
    console.log('\n' + '='.repeat(60));
    let monotonic = true;
    for (let i = 1; i < results.length; i++) {
        if (results[i].alpha > results[i - 1].alpha) {
            monotonic = false;
            console.log(`✗ ERROR: Loss increased from ${results[i-1].label} to ${results[i].label}`);
        }
    }

    if (monotonic) {
        console.log('✓ CORRECT: Loss decreases monotonically with plating thickness');
        console.log('  As expected: more side coverage → lower loss');
        return true;
    } else {
        console.log('✗ ERROR: Loss should decrease monotonically');
        return false;
    }
}

// Run tests
async function run_all_tests() {
    let test1_pass = false;
    let test2_pass = false;

    try {
        test1_pass = await test_top_only_plating();
        test2_pass = await test_top_plating_thickness_variation();

        console.log('\n' + '='.repeat(60));
        if (test1_pass && test2_pass) {
            console.log('✓ ALL TESTS PASSED');
            console.log('Top-only plating correctly extends down sides');
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
