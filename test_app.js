/**
 * Command-line test for web app solver integration.
 * Tests the same code paths as app.js without requiring a browser.
 * Run with: node test_app.js
 */

import { MicrostripSolver } from './microstrip.js';
import { GroundedCPWSolver2D } from './gcpw.js';

// Simulate the getParams() function from app.js with typical default values
function getDefaultParams() {
    return {
        tl_type: 'microstrip',
        w: 3.0 * 1e-3,      // 3mm trace width
        h: 1.6 * 1e-3,      // 1.6mm substrate height
        t: 35 * 1e-6,       // 35um trace thickness
        er: 4.5,            // FR4 permittivity
        tand: 0.02,         // Loss tangent
        sigma: 5.8e7,       // Copper conductivity
        freq: 1.0 * 1e9,    // 1 GHz
        nx: 100,
        ny: 100,
        // GCPW specific parameters
        gap: 0.2 * 1e-3,
        top_gnd_w: 1.0 * 1e-3,
        via_gap: 0.5 * 1e-3,
        via_d: 0.3 * 1e-3,
    };
}

// Test microstrip solver creation (mirrors app.js updateGeometry)
async function testMicrostripCreation() {
    console.log("Test 1: Microstrip solver creation...");
    const p = getDefaultParams();
    p.tl_type = 'microstrip';

    const solver = new MicrostripSolver({
        trace_width: p.w,
        substrate_height: p.h,
        trace_thickness: p.t,
        epsilon_r: p.er,
        tan_delta: p.tand,
        sigma_cond: p.sigma,
        freq: p.freq,
        nx: p.nx,
        ny: p.ny
    });

    console.log(`  Grid: ${solver.x.length}x${solver.y.length}`);

    if (solver.x.length < 10 || solver.y.length < 10) {
        throw new Error("Grid too small");
    }
    if (!solver.V || !solver.epsilon_r || !solver.conductor_mask) {
        throw new Error("Solver arrays not initialized");
    }

    console.log("  PASS: Microstrip solver created successfully");
    return solver;
}

// Test full analysis (mirrors app.js runSimulation)
async function testMicrostripAnalysis(solver) {
    console.log("Test 2: Microstrip full analysis...");

    const results = await solver.perform_analysis();

    console.log(`  Z0: ${results.Z0.toFixed(2)} Ω`);
    console.log(`  eps_eff: ${results.eps_eff.toFixed(3)}`);
    console.log(`  Diel loss: ${results.alpha_diel_db_m.toFixed(4)} dB/m`);
    console.log(`  Cond loss: ${results.alpha_cond_db_m.toFixed(4)} dB/m`);

    // Validate results are in reasonable range
    if (results.Z0 < 10 || results.Z0 > 200) {
        throw new Error(`Z0 out of range: ${results.Z0}`);
    }
    if (results.eps_eff < 1 || results.eps_eff > 10) {
        throw new Error(`eps_eff out of range: ${results.eps_eff}`);
    }
    if (!results.RLGC || typeof results.RLGC.R !== 'number') {
        throw new Error("RLGC not computed");
    }
    if (!solver.solution_valid) {
        throw new Error("solution_valid not set");
    }
    if (!solver.Ex || !solver.Ey) {
        throw new Error("E-field arrays not computed");
    }

    console.log("  PASS: Analysis completed successfully");
    return results;
}

// Test GCPW solver creation
async function testGCPWCreation() {
    console.log("Test 3: GCPW solver creation...");
    const p = getDefaultParams();
    p.tl_type = 'gcpw';

    const solver = new GroundedCPWSolver2D(
        p.h, p.w, p.t, p.gap, p.top_gnd_w, p.via_gap, p.via_d,
        35e-6, p.er, p.tand, 0.0, p.sigma, 1,
        null, null, p.freq, p.nx, p.ny
    );

    console.log(`  Grid: ${solver.x.length}x${solver.y.length}`);

    if (solver.x.length < 10 || solver.y.length < 10) {
        throw new Error("Grid too small");
    }

    console.log("  PASS: GCPW solver created successfully");
    return solver;
}

// Test Stripline solver creation and analysis
async function testStriplineCreation() {
    console.log("Test 5: Stripline solver creation...");

    const solver = new MicrostripSolver({
        substrate_height: 0.2e-3,
        trace_width: 0.15e-3,
        trace_thickness: 35e-6,
        gnd_thickness: 16e-6,
        epsilon_r: 4.1,
        epsilon_r_top: 4.1,
        air_top: 0.2e-3,
        tan_delta: 0.02,
        sigma_cond: 5.8e7,
        freq: 1e9,
        nx: 100,
        ny: 100,
        boundaries: ["open", "open", "gnd", "gnd"]
    });

    console.log(`  Grid: ${solver.x.length}x${solver.y.length}`);

    const results = await solver.perform_analysis();
    console.log(`  Z0: ${results.Z0.toFixed(2)} Ω, eps_eff: ${results.eps_eff.toFixed(3)}`);

    // Stripline should have eps_eff close to substrate permittivity
    if (results.eps_eff < 3.5 || results.eps_eff > 4.5) {
        throw new Error(`Stripline eps_eff out of range: ${results.eps_eff}`);
    }

    console.log("  PASS: Stripline solver created successfully");
    return solver;
}

// Test Microstrip with Solder Mask
async function testSolderMaskCreation() {
    console.log("Test 6: Microstrip with Solder Mask...");

    const solver = new MicrostripSolver({
        substrate_height: 0.1e-3,
        trace_width: 0.1e-3,
        trace_thickness: 35e-6,
        gnd_thickness: 35e-6,
        epsilon_r: 4.5,
        tan_delta: 0.02,
        sigma_cond: 5.8e7,
        freq: 1e9,
        nx: 100,
        ny: 100,
        use_sm: true,
        sm_t_sub: 20e-6,
        sm_t_trace: 20e-6,
        sm_t_side: 20e-6,
        sm_er: 3.5,
        sm_tand: 0.02,
        boundaries: ["open", "open", "open", "gnd"]
    });

    console.log(`  Grid: ${solver.x.length}x${solver.y.length}`);

    const results = await solver.perform_analysis();
    console.log(`  Z0: ${results.Z0.toFixed(2)} Ω, eps_eff: ${results.eps_eff.toFixed(3)}`);

    if (results.Z0 < 10 || results.Z0 > 200) {
        throw new Error(`Solder mask Z0 out of range: ${results.Z0}`);
    }

    console.log("  PASS: Solder Mask microstrip created successfully");
    return solver;
}

// Test Embedded Microstrip (Top Dielectric)
async function testTopDielectricCreation() {
    console.log("Test 7: Embedded Microstrip (Top Dielectric)...");

    const solver = new MicrostripSolver({
        substrate_height: 1.6e-3,
        trace_width: 3e-3,
        trace_thickness: 35e-6,
        gnd_thickness: 35e-6,
        epsilon_r: 4.5,
        tan_delta: 0.02,
        sigma_cond: 5.8e7,
        freq: 1e9,
        nx: 100,
        ny: 100,
        top_diel_h: 0.2e-3,
        top_diel_er: 4.5,
        boundaries: ["open", "open", "open", "gnd"]
    });

    console.log(`  Grid: ${solver.x.length}x${solver.y.length}`);

    const results = await solver.perform_analysis();
    console.log(`  Z0: ${results.Z0.toFixed(2)} Ω, eps_eff: ${results.eps_eff.toFixed(3)}`);

    // Embedded microstrip should have higher eps_eff than regular microstrip
    if (results.eps_eff < 3.0 || results.eps_eff > 5.0) {
        throw new Error(`Top dielectric eps_eff out of range: ${results.eps_eff}`);
    }

    console.log("  PASS: Top Dielectric microstrip created successfully");
    return solver;
}

// Test Microstrip with Ground Cutout
async function testGroundCutoutCreation() {
    console.log("Test 8: Microstrip with Ground Cutout...");

    const solver = new MicrostripSolver({
        substrate_height: 1.6e-3,
        trace_width: 3e-3,
        trace_thickness: 35e-6,
        gnd_thickness: 35e-6,
        epsilon_r: 4.5,
        tan_delta: 0.02,
        sigma_cond: 5.8e7,
        freq: 1e9,
        nx: 100,
        ny: 100,
        gnd_cut_width: 3e-3,
        gnd_cut_sub_h: 1e-3,
        boundaries: ["open", "open", "open", "gnd"]
    });

    console.log(`  Grid: ${solver.x.length}x${solver.y.length}`);

    const results = await solver.perform_analysis();
    console.log(`  Z0: ${results.Z0.toFixed(2)} Ω, eps_eff: ${results.eps_eff.toFixed(3)}`);

    // Ground cutout should increase Z0 compared to regular microstrip
    if (results.Z0 < 40 || results.Z0 > 100) {
        throw new Error(`Ground cutout Z0 out of range: ${results.Z0}`);
    }

    console.log("  PASS: Ground Cutout microstrip created successfully");
    return solver;
}

// Test draw() data access patterns (mirrors app.js draw() exactly)
function testDrawDataAccess(solver) {
    console.log("Test 4: Draw data access patterns...");

    const nx = solver.x.length;
    const ny = solver.y.length;

    // Exactly mirror draw() function from app.js
    const yArr = Array.from(solver.y);
    const maxY = Math.min(yArr[ny - 1], solver.h + 5 * solver.w);
    const maxYIdx = yArr.findIndex(y => y > maxY);
    const nyDisplay = maxYIdx > 0 ? maxYIdx : ny;

    if (nyDisplay <= 0) {
        throw new Error(`Invalid nyDisplay: ${nyDisplay}`);
    }

    // Convert Float64Arrays to regular arrays (as in app.js)
    const xMM = Array.from(solver.x, v => v * 1000);
    const yMM = yArr.slice(0, nyDisplay).map(v => v * 1000);

    if (xMM.length !== nx) {
        throw new Error(`xMM length mismatch: ${xMM.length} vs ${nx}`);
    }
    if (yMM.length !== nyDisplay) {
        throw new Error(`yMM length mismatch: ${yMM.length} vs ${nyDisplay}`);
    }

    // Test geometry view data access
    let zData = [];
    for (let i = 0; i < nyDisplay; i++) {
        const row = [];
        for (let j = 0; j < nx; j++) {
            if (solver.conductor_mask[i][j]) {
                row.push(solver.ground_mask[i][j] ? -2 : -1);
            } else {
                row.push(solver.epsilon_r[i][j]);
            }
        }
        zData.push(row);
    }

    if (zData.length !== nyDisplay) {
        throw new Error(`zData rows mismatch: ${zData.length} vs ${nyDisplay}`);
    }
    if (zData[0].length !== nx) {
        throw new Error(`zData cols mismatch: ${zData[0].length} vs ${nx}`);
    }

    // Test potential view (if solution valid)
    if (solver.solution_valid) {
        zData = [];
        for (let i = 0; i < nyDisplay; i++) {
            zData.push(Array.from(solver.V[i].slice(0, nx)));
        }

        // Test E-field view
        zData = [];
        for (let i = 0; i < nyDisplay; i++) {
            const row = [];
            for (let j = 0; j < nx; j++) {
                const E = Math.hypot(solver.Ex[i][j], solver.Ey[i][j]);
                if (isNaN(E)) {
                    throw new Error(`E-field contains NaN at [${i}][${j}]`);
                }
                row.push(E);
            }
            zData.push(row);
        }
    }

    console.log(`  nyDisplay: ${nyDisplay}, xMM: ${xMM.length}, yMM: ${yMM.length}`);
    console.log("  PASS: Draw data access successful");
}

// Main test runner
async function runTests() {
    console.log("=== Web App Integration Tests ===\n");

    let passed = 0;
    let failed = 0;

    const tests = [
        async () => {
            const solver = await testMicrostripCreation();
            await testMicrostripAnalysis(solver);
            testDrawDataAccess(solver);
            return 3; // 3 tests passed
        },
        async () => { await testGCPWCreation(); return 1; },
        async () => { await testStriplineCreation(); return 1; },
        async () => { await testSolderMaskCreation(); return 1; },
        async () => { await testTopDielectricCreation(); return 1; },
        async () => { await testGroundCutoutCreation(); return 1; },
    ];

    for (const test of tests) {
        try {
            passed += await test();
        } catch (e) {
            console.error(`  FAIL: ${e.message}`);
            failed++;
        }
    }

    console.log(`\n=== Results: ${passed} passed, ${failed} failed ===`);

    if (failed > 0) {
        process.exit(1);
    }
}

runTests();
