import { MicrostripSolver } from '../src/microstrip.js';

/**
 * Test: Verify that very thick plating (>> skin depth) with no roughness
 * equals a solid conductor made of the plating material
 */

async function test_thick_plating_equivalence() {
    console.log('Testing thick plating equivalence...\n');

    const freq = 1e9;  // 1 GHz
    const gold_sigma = 4.1e7;  // Gold conductivity
    const copper_sigma = 5.96e7;  // Copper conductivity

    // Calculate skin depth in gold at 1 GHz
    const mu0 = 4 * Math.PI * 1e-7;
    const skin_depth_gold = Math.sqrt(2 / (2 * Math.PI * freq * mu0 * gold_sigma));
    console.log(`Gold skin depth at 1 GHz: ${(skin_depth_gold * 1e6).toFixed(3)} μm`);

    // Case 1: Very thick gold plating (10× skin depth) on copper, no roughness
    // Top and sides plated, bottom not plated
    const thick_plating = {
        sigma: gold_sigma,
        thickness: 10 * skin_depth_gold,  // Much thicker than skin depth
        rq: 0,  // No roughness
        top: true,
        sides: true,
        bottom: false
    };

    const solver1 = new MicrostripSolver({
        substrate_height: 0.8e-3,
        trace_width: 100e-6,
        trace_thickness: 35e-6,
        gnd_thickness: 35e-6,
        epsilon_r: 4.5,
        tan_delta: 0.02,
        sigma_cond: copper_sigma,  // Copper bulk
        freq: freq,
        nx: 50,
        ny: 25,
        plating: thick_plating,
        boundaries: ["open", "open", "open", "gnd"],
        rq: 0  // No bulk roughness
    });

    // Case 2: Solid gold conductor, no plating, no roughness
    const solver2 = new MicrostripSolver({
        substrate_height: 0.8e-3,
        trace_width: 100e-6,
        trace_thickness: 35e-6,
        gnd_thickness: 35e-6,
        epsilon_r: 4.5,
        tan_delta: 0.02,
        sigma_cond: gold_sigma,  // Gold bulk (no plating)
        freq: freq,
        nx: 50,
        ny: 25,
        plating: null,  // No plating
        boundaries: ["open", "open", "open", "gnd"],
        rq: 0  // No roughness
    });

    const result1 = await solver1.solve_adaptive({ max_refinements: 0 });
    const mode1 = result1.modes[0];

    const result2 = await solver2.solve_adaptive({ max_refinements: 0 });
    const mode2 = result2.modes[0];

    console.log('\nCase 1: Thick gold plating on copper (top + sides)');
    console.log(`  Z0 = ${mode1.Z0.toFixed(4)} Ω`);
    console.log(`  Conductor loss = ${mode1.alpha_c.toFixed(6)} dB/m`);

    console.log('\nCase 2: Solid gold conductor (no plating)');
    console.log(`  Z0 = ${mode2.Z0.toFixed(4)} Ω`);
    console.log(`  Conductor loss = ${mode2.alpha_c.toFixed(6)} dB/m`);

    // Check Z0 equivalence
    const z0_diff = Math.abs(mode1.Z0 - mode2.Z0);
    const z0_rel_diff = z0_diff / mode2.Z0;
    console.log(`\nZ0 difference: ${z0_diff.toFixed(4)} Ω (${(z0_rel_diff * 100).toFixed(3)}%)`);

    // Check conductor loss equivalence
    const alpha_diff = Math.abs(mode1.alpha_c - mode2.alpha_c);
    const alpha_rel_diff = alpha_diff / mode2.alpha_c;
    console.log(`Conductor loss difference: ${alpha_diff.toFixed(6)} dB/m (${(alpha_rel_diff * 100).toFixed(3)}%)`);

    // Tolerance check
    // 2%: the equivalence compares alpha_c, which blends in R_dc — and a
    // copper-CORE conductor genuinely has lower R_dc than solid gold (~1.1%
    // of R_total here). The gap was hidden while the loss integrand
    // over-predicted R_ac; with the calibrated vacuum-field integrand the
    // real R_dc difference shows.
    const tolerance = 0.02;

    if (z0_rel_diff > tolerance) {
        console.log(`\n⚠ WARNING: Z0 differs by ${(z0_rel_diff * 100).toFixed(3)}% (> ${tolerance * 100}% tolerance)`);
        console.log('   Thick plating should equal solid conductor of plating material');
    } else {
        console.log(`\n✓ Z0 match within ${(tolerance * 100).toFixed(1)}% tolerance`);
    }

    if (alpha_rel_diff > tolerance) {
        console.log(`⚠ WARNING: Conductor loss differs by ${(alpha_rel_diff * 100).toFixed(3)}% (> ${tolerance * 100}% tolerance)`);
        console.log('   Thick plating should equal solid conductor of plating material');
    } else {
        console.log(`✓ Conductor loss match within ${(tolerance * 100).toFixed(1)}% tolerance`);
    }

    // Additional test: check that bottom surface (unplated) doesn't affect result
    // when plating is very thick on top and sides
    console.log('\n' + '='.repeat(60));
    console.log('Checking bottom surface contribution...');

    // Case 3: Same as Case 1 but with bottom also plated
    const thick_plating_all = {
        sigma: gold_sigma,
        thickness: 10 * skin_depth_gold,
        rq: 0,
        top: true,
        sides: true,
        bottom: true  // Now bottom is also plated
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
        plating: thick_plating_all,
        boundaries: ["open", "open", "open", "gnd"],
        rq: 0
    });

    const result3 = await solver3.solve_adaptive({ max_refinements: 0 });
    const mode3 = result3.modes[0];

    console.log('\nCase 3: Thick gold plating on all surfaces (top + sides + bottom)');
    console.log(`  Z0 = ${mode3.Z0.toFixed(4)} Ω`);
    console.log(`  Conductor loss = ${mode3.alpha_c.toFixed(6)} dB/m`);

    const alpha_diff_bottom = Math.abs(mode1.alpha_c - mode3.alpha_c);
    const alpha_rel_diff_bottom = alpha_diff_bottom / mode1.alpha_c;
    console.log(`\nLoss difference (bottom plated vs not): ${alpha_diff_bottom.toFixed(6)} dB/m (${(alpha_rel_diff_bottom * 100).toFixed(3)}%)`);

    if (alpha_rel_diff_bottom > 0.05) {
        console.log(`⚠ WARNING: Bottom surface affects loss significantly even with thick plating on top/sides`);
        console.log('   With very thick plating, bottom surface should have minimal impact');
    } else {
        console.log(`✓ Bottom surface has minimal impact (< 5%) with thick top/side plating`);
    }

    // Most important check: ALL surfaces plated should equal solid conductor
    console.log('\n' + '='.repeat(60));
    console.log('Checking equivalence: All surfaces plated vs solid conductor...');

    const alpha_diff_solid = Math.abs(mode3.alpha_c - mode2.alpha_c);
    const alpha_rel_diff_solid = alpha_diff_solid / mode2.alpha_c;
    console.log(`\nCase 3 (all plated) vs Case 2 (solid gold):`);
    console.log(`  Loss difference: ${alpha_diff_solid.toFixed(6)} dB/m (${(alpha_rel_diff_solid * 100).toFixed(3)}%)`);

    if (alpha_rel_diff_solid > tolerance) {
        console.log(`\n⚠ WARNING: All-surface plating differs from solid conductor by ${(alpha_rel_diff_solid * 100).toFixed(3)}%`);
        console.log('   This should match within tolerance when plating >> skin depth');
        return false;
    } else {
        console.log(`\n✓ All-surface thick plating equals solid conductor within ${(tolerance * 100).toFixed(1)}% tolerance`);
        return true;
    }
}

// Run test
test_thick_plating_equivalence()
    .then(passed => {
        console.log('\n' + '='.repeat(60));
        if (passed) {
            console.log('✓ Test PASSED: Thick plating equivalence verified');
        } else {
            console.log('✗ Test FAILED: Thick plating does not equal solid conductor');
            process.exit(1);
        }
        console.log('='.repeat(60));
    })
    .catch(err => {
        console.error('\nTest failed:', err);
        process.exit(1);
    });
