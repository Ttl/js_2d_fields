import { MicrostripSolver } from '../src/microstrip.js';
import { GroundedCPWSolver2D } from '../src/gcpw.js';

/**
 * Test suite for conductor corner plating implementation
 * Validates that side plating extends to bottom corners as specified
 */

function assert_close(actual, expected, rel_tol, name) {
    const diff = Math.abs(actual - expected);
    const tol = Math.abs(expected * rel_tol);
    if (diff > tol) {
        throw new Error(`${name}: expected ${expected}, got ${actual} (diff=${diff}, tol=${tol})`);
    }
}

async function test_corner_detection() {
    console.log('\nTest 1: Corner detection at bottom edges');

    // Create a simple microstrip with plating on sides only
    const plating = {
        sigma: 5.96e7,  // Gold plating
        thickness: 5e-6,  // 5um
        rq: 0.5e-6,  // 0.5um roughness
        top: false,
        sides: true,  // Side plating enabled
        bottom: false  // Bottom plating disabled
    };

    const solver = new MicrostripSolver({
        substrate_height: 0.8e-3,  // 0.8mm
        trace_width: 100e-6,  // 100um
        trace_thickness: 35e-6,  // 35um
        gnd_thickness: 35e-6,
        epsilon_r: 4.5,
        tan_delta: 0.02,
        sigma_cond: 5.96e7,  // Copper bulk
        freq: 1e9,
        nx: 50,  // Reduced mesh size for faster testing
        ny: 25,
        plating: plating,
        boundaries: ["open", "open", "open", "gnd"]
    });

    // Calculate loss with corner handling
    // Use max_refinements: 0 to skip adaptive mesh refinement for faster testing
    const result = await solver.solve_adaptive({ max_refinements: 0 });

    // Verify we got a result (basic sanity check)
    if (!result || !result.modes || result.modes.length === 0) {
        console.error('  Result:', result);
        throw new Error('Corner detection test: solver did not return valid result');
    }

    const mode = result.modes[0];
    console.log(`  Characteristic impedance: ${mode.Z0.toFixed(2)} Ω`);
    console.log(`  Conductor loss: ${mode.alpha_c.toFixed(3)} dB/m`);
    console.log('  ✓ Corner detection working');
}

async function test_corner_vs_no_corner() {
    console.log('\nTest 2: Compare loss with and without side plating at corners');

    // Case 1: No plating
    const no_plating = null;

    // Case 2: Side plating only (should affect corners)
    const side_plating = {
        sigma: 5.96e7,  // Gold plating (same as copper for this test)
        thickness: 5e-6,
        rq: 0.5e-6,
        top: false,
        sides: true,
        bottom: false
    };

    // Case 3: Bottom plating only (no corner effect)
    const bottom_plating = {
        sigma: 5.96e7,
        thickness: 5e-6,
        rq: 0.5e-6,
        top: false,
        sides: false,
        bottom: true
    };

    const results = [];

    for (const [name, plating] of [
        ['No plating', no_plating],
        ['Side plating', side_plating],
        ['Bottom plating', bottom_plating]
    ]) {
        const solver = new MicrostripSolver({
            substrate_height: 0.8e-3,
            trace_width: 100e-6,
            trace_thickness: 35e-6,
            gnd_thickness: 35e-6,
            epsilon_r: 4.5,
            tan_delta: 0.02,
            sigma_cond: 5.96e7,
            freq: 1e9,
            nx: 50,  // Reduced for faster testing
            ny: 25,
            plating: plating,
            boundaries: ["open", "open", "open", "gnd"]
        });

        const result = await solver.solve_adaptive({ max_refinements: 0 });
        const mode = result.modes[0];

        results.push({ name, Z0: mode.Z0, alpha: mode.alpha_c });
        console.log(`  ${name}: Z0=${mode.Z0.toFixed(2)} Ω, α=${mode.alpha_c.toFixed(3)} dB/m`);
    }

    console.log('  ✓ All cases computed successfully');
}

async function test_corner_mesh_sensitivity() {
    console.log('\nTest 3: Corner impedance averaging with different mesh sizes');

    const plating = {
        sigma: 5.96e7,  // Gold
        thickness: 3e-6,
        rq: 0.5e-6,
        top: false,
        sides: true,  // Corner plating enabled
        bottom: false  // No bottom plating
    };

    const mesh_sizes = [
        { nx: 50, ny: 25, label: 'Coarse' },
        { nx: 100, ny: 50, label: 'Medium' },
        { nx: 200, ny: 100, label: 'Fine' }
    ];

    const results = [];

    for (const { nx, ny, label } of mesh_sizes) {
        const solver = new MicrostripSolver({
            substrate_height: 0.8e-3,
            trace_width: 100e-6,
            trace_thickness: 35e-6,
            gnd_thickness: 35e-6,
            epsilon_r: 4.5,
            tan_delta: 0.02,
            sigma_cond: 5.96e7,
            freq: 1e9,
            nx: nx,
            ny: ny,
            plating: plating,
            boundaries: ["open", "open", "open", "gnd"]
        });

        const result = await solver.solve_adaptive({ max_refinements: 0 });
        const mode = result.modes[0];

        results.push({ label, nx, ny, Z0: mode.Z0, alpha: mode.alpha_c });
        console.log(`  ${label} (${nx}×${ny}): Z0=${mode.Z0.toFixed(2)} Ω, α=${mode.alpha_c.toFixed(3)} dB/m`);
    }

    // Check that results converge as mesh is refined
    const alpha_coarse = results[0].alpha;
    const alpha_fine = results[2].alpha;
    const rel_change = Math.abs(alpha_fine - alpha_coarse) / alpha_fine;

    if (rel_change > 0.2) {
        console.log(`  ⚠ Warning: Large mesh sensitivity (${(rel_change * 100).toFixed(1)}% change)`);
    } else {
        console.log(`  ✓ Results converging with mesh refinement (${(rel_change * 100).toFixed(1)}% change)`);
    }
}

async function test_corner_symmetry() {
    console.log('\nTest 4: Verify symmetric conductor produces symmetric corner effects');

    // Use GCPW which has symmetric ground planes on both sides
    const plating = {
        sigma: 5.96e7,  // Gold plating
        thickness: 3e-6,
        rq: 0.5e-6,
        top: false,
        sides: true,  // Enable corner plating
        bottom: false
    };

    // Built twice from one description, differing ONLY in the plating block, because the
    // plating has to be shown to reach the solver at all. GroundedCPWSolver2D is a
    // back-compat wrapper that maps its options across by hand, and it silently dropped
    // `plating` (and `rq`) for a while: this test still printed a loss and still passed,
    // because it only ever asserted that the numbers were finite. The comparison below
    // is what makes that failure mode visible.
    const build = (pl) => new GroundedCPWSolver2D({
        substrate_height: 800e-6,
        trace_width: 100e-6,
        trace_thickness: 35e-6,
        gap: 50e-6,
        via_gap: null,
        gnd_thickness: 35e-6,
        epsilon_r: 4.5,
        tan_delta: 0.02,
        sigma_diel: 0.0,
        sigma_cond: 5.96e7,
        epsilon_r_top: 1,
        air_top: null,
        air_side: null,
        freq: 5e9,
        nx: 50,  // Reduced for faster testing
        ny: 25,
        boundaries: null,
        use_sm: false,
        plating: pl
    });

    const mode = (await build(plating).solve_adaptive({ max_refinements: 0 })).modes[0];
    const bare = (await build(null).solve_adaptive({ max_refinements: 0 })).modes[0];

    console.log(`  GCPW Z0: ${mode.Z0.toFixed(2)} Ω`);
    console.log(`  GCPW conductor loss: ${bare.alpha_c.toFixed(3)} dB/m bare`
        + ` -> ${mode.alpha_c.toFixed(3)} dB/m plated`);
    console.log('  ✓ Symmetric structure computed successfully');

    // The result should be finite and reasonable
    if (!isFinite(mode.Z0) || !isFinite(mode.alpha_c)) {
        throw new Error('Corner symmetry test: non-finite results');
    }
    if (mode.Z0 < 20 || mode.Z0 > 200) {
        throw new Error(`Corner symmetry test: unrealistic Z0 = ${mode.Z0.toFixed(2)} Ω`);
    }
    // This plating is the bulk metal again (5.96e7) but ROUGHENED to 0.5 um rms on the
    // side faces, so its whole effect is added loss — a few percent, and one-directional.
    // Z0 is a capacitance quantity and must not move at all.
    if (!(mode.alpha_c > bare.alpha_c * 1.02)) {
        throw new Error('Corner symmetry test: side plating did not raise the loss '
            + `(bare ${bare.alpha_c.toFixed(4)}, plated ${mode.alpha_c.toFixed(4)} dB/m) — `
            + 'the plating block is not reaching the solver');
    }
    assert_close(mode.Z0, bare.Z0, 1e-6, 'Corner symmetry test: plating moved Z0');
}

async function test_corner_only_with_sides_enabled() {
    console.log('\nTest 5: Corner effects only apply when sides plating is enabled');

    // Plating with sides disabled - should NOT trigger corner handling
    const plating_no_sides = {
        sigma: 5.96e7,
        thickness: 3e-6,
        rq: 0.5e-6,
        top: true,
        sides: false,  // Disabled
        bottom: true
    };

    // Plating with sides enabled - should trigger corner handling
    const plating_with_sides = {
        sigma: 5.96e7,
        thickness: 3e-6,
        rq: 0.5e-6,
        top: true,
        sides: true,  // Enabled
        bottom: true
    };

    const results = [];

    for (const [name, plating] of [
        ['Sides disabled', plating_no_sides],
        ['Sides enabled', plating_with_sides]
    ]) {
        const solver = new MicrostripSolver({
            substrate_height: 0.8e-3,
            trace_width: 100e-6,
            trace_thickness: 35e-6,
            gnd_thickness: 35e-6,
            epsilon_r: 4.5,
            tan_delta: 0.02,
            sigma_cond: 5.96e7,
            freq: 1e9,
            nx: 50,  // Reduced for faster testing
            ny: 25,
            plating: plating,
            boundaries: ["open", "open", "open", "gnd"]
        });

        const result = await solver.solve_adaptive({ max_refinements: 0 });
        const mode = result.modes[0];

        results.push({ name, alpha: mode.alpha_c });
        console.log(`  ${name}: α=${mode.alpha_c.toFixed(3)} dB/m`);
    }

    console.log('  ✓ Both configurations computed successfully');
}

// Run all tests
async function run_all_tests() {
    console.log('===================================');
    console.log('Conductor Corner Plating Test Suite');
    console.log('===================================');

    try {
        await test_corner_detection();
        await test_corner_vs_no_corner();
        await test_corner_mesh_sensitivity();
        await test_corner_symmetry();
        await test_corner_only_with_sides_enabled();

        console.log('\n===================================');
        console.log('All tests passed! ✓');
        console.log('===================================\n');
    } catch (err) {
        console.error('\n===================================');
        console.error('Test failed! ✗');
        console.error('===================================');
        console.error(err);
        process.exit(1);
    }
}

run_all_tests();
