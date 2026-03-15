/**
 * Test interpolating sweep against discrete sweep.
 * Verifies that interpolated RLGC values match exact values within tolerance.
 */

import { MicrostripSolver } from '../src/microstrip.js';
import { CubicSpline, InterpolatingSweep } from '../src/interpolating_sweep.js';

// --- CubicSpline unit tests ---

console.log('Testing CubicSpline');
console.log('====================\n');

// Test 1: Linear function should be interpolated exactly
{
    const x = [0, 1, 2, 3, 4];
    const y = x.map(v => 2 * v + 1);
    const spline = new CubicSpline(x, y);

    let maxErr = 0;
    for (let xv = 0; xv <= 4; xv += 0.1) {
        const exact = 2 * xv + 1;
        const interp = spline.evaluate(xv);
        maxErr = Math.max(maxErr, Math.abs(exact - interp));
    }
    console.log(`  Linear function max error: ${maxErr.toExponential(2)} ${maxErr < 1e-10 ? 'PASS' : 'FAIL'}`);
}

// Test 2: Quadratic function
{
    const x = [0, 1, 2, 3, 4, 5];
    const y = x.map(v => v * v);
    const spline = new CubicSpline(x, y);

    let maxErr = 0;
    for (let xv = 0; xv <= 5; xv += 0.1) {
        const exact = xv * xv;
        const interp = spline.evaluate(xv);
        maxErr = Math.max(maxErr, Math.abs(exact - interp));
    }
    console.log(`  Quadratic function max error: ${maxErr.toExponential(2)} ${maxErr < 0.15 ? 'PASS' : 'FAIL'}`);
}

// Test 3: sin(x) with 10 points - should be quite accurate
{
    const n = 10;
    const x = [];
    const y = [];
    for (let i = 0; i < n; i++) {
        const xv = i * Math.PI / (n - 1);
        x.push(xv);
        y.push(Math.sin(xv));
    }
    const spline = new CubicSpline(x, y);

    let maxErr = 0;
    for (let xv = 0; xv <= Math.PI; xv += 0.01) {
        const exact = Math.sin(xv);
        const interp = spline.evaluate(xv);
        maxErr = Math.max(maxErr, Math.abs(exact - interp));
    }
    console.log(`  sin(x) with ${n} points max error: ${maxErr.toExponential(2)} ${maxErr < 1e-4 ? 'PASS' : 'FAIL'}`);
}

// Test 4: Two points (linear fallback)
{
    const spline = new CubicSpline([0, 1], [0, 1]);
    const val = spline.evaluate(0.5);
    console.log(`  Two-point linear: ${Math.abs(val - 0.5) < 1e-10 ? 'PASS' : 'FAIL'}`);
}

// --- InterpolatingSweep test ---

console.log('\nTesting InterpolatingSweep');
console.log('==========================\n');

const options = {
    trace_width: 0.2e-3,
    substrate_height: 0.2e-3,
    trace_thickness: 35e-6,
    epsilon_r: 4.4,
    tan_delta: 0.02,
    sigma_cond: 5.8e7,
    freq: 10e9,
    nx: 20,
    ny: 20,
};

const solver = new MicrostripSolver(options);

(async () => {
    // First solve to get cached results
    const results = await solver.solve_adaptive({
        max_iters: 3,
        energy_tol: 0.05,
        param_tol: 0.05,
        max_nodes: 5000
    });

    const cachedResults = results;
    const frequencies = [];
    const fStart = 0.1e9;
    const fStop = 10e9;
    const nPoints = 100;

    for (let i = 0; i < nPoints; i++) {
        frequencies.push(fStart + (fStop - fStart) * i / (nPoints - 1));
    }

    // Compute discrete sweep (reference)
    console.log('  Computing discrete sweep...');
    const discreteResults = [];
    for (const freq of frequencies) {
        const result = await solver.computeAtFrequency(freq, cachedResults);
        discreteResults.push({ freq, result });
    }

    // Compute interpolating sweep
    console.log('  Computing interpolating sweep...');
    const sweep = new InterpolatingSweep(solver, cachedResults, {
        tolerance: 0.001,
        initialPoints: 8
    });

    const nSamples = await sweep.run(fStart, fStop);
    const interpResults = sweep.buildResults(frequencies);

    console.log(`  Interpolating sweep used ${nSamples} exact solves for ${nPoints} output points`);

    // Compare RLGC values
    let maxRelErr = { R: 0, L: 0, G: 0, C: 0 };
    let maxRelErrOverall = 0;

    for (let i = 0; i < nPoints; i++) {
        const disc = discreteResults[i].result.modes[0].RLGC;
        const interp = interpResults[i].result.modes[0].RLGC;

        for (const key of ['R', 'L', 'G', 'C']) {
            const floor = key === 'G' ? 1e-6 : 1e-12;
            const denom = Math.max(Math.abs(disc[key]), floor);
            const err = Math.abs(disc[key] - interp[key]) / denom;
            if (err > maxRelErr[key]) maxRelErr[key] = err;
            if (err > maxRelErrOverall) maxRelErrOverall = err;
        }
    }

    console.log('\n  Max relative errors:');
    for (const key of ['R', 'L', 'G', 'C']) {
        const pct = (maxRelErr[key] * 100).toFixed(4);
        console.log(`    ${key}: ${pct}% ${maxRelErr[key] < 0.01 ? 'PASS' : 'FAIL'}`);
    }
    console.log(`\n  Overall max relative error: ${(maxRelErrOverall * 100).toFixed(4)}% ${maxRelErrOverall < 0.01 ? 'PASS' : 'FAIL'}`);

    // Compare derived quantities
    let maxZcErr = 0;
    let maxEpsErr = 0;
    for (let i = 0; i < nPoints; i++) {
        const disc = discreteResults[i].result.modes[0];
        const interp = interpResults[i].result.modes[0];

        const zcErr = Math.abs(disc.Zc.re - interp.Zc.re) / Math.max(Math.abs(disc.Zc.re), 1e-6);
        const epsErr = Math.abs(disc.eps_eff - interp.eps_eff) / Math.max(Math.abs(disc.eps_eff), 1e-6);
        if (zcErr > maxZcErr) maxZcErr = zcErr;
        if (epsErr > maxEpsErr) maxEpsErr = epsErr;
    }

    console.log(`\n  Zc(re) max relative error: ${(maxZcErr * 100).toFixed(4)}% ${maxZcErr < 0.01 ? 'PASS' : 'FAIL'}`);
    console.log(`  eps_eff max relative error: ${(maxEpsErr * 100).toFixed(4)}% ${maxEpsErr < 0.01 ? 'PASS' : 'FAIL'}`);

    console.log('\nDone.');
})();
