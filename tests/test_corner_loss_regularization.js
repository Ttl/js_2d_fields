// QS conductor-loss integrand regression.
//
// Production (rect-based solvers) uses the VACUUM-field integrand: H_t from the
// vacuum (C0) solve — the quasi-TEM surface-current pattern, ε-independent —
// scaled by Z0_vac² and the global VACUUM_LOSS_CAL. This replaced the
// dielectric-field integrand whose substrate-interface corner singularity made
// the sum practically mesh-divergent and +21–26% high vs the tri-backend MQS
// volume loss. This test checks the three properties that formulation bought:
//   1. mesh stability — R moves < 2.5% between a coarse and a dense FDM mesh;
//   2. accuracy — R within 8% of the tri-backend MQS reference (measured ~+4%);
//   3. ε-invariance — R at εr 2.2 and 9.8 agrees within 1% (the loss integrand
//      must not see the dielectric; a leak of the dielectric fields or the √εr
//      factor into the loss path breaks this immediately).

import { MicrostripSolver } from '../src/microstrip.js';

const GEOM = {
    substrate_height: 0.21e-3, trace_width: 0.35e-3, trace_thickness: 35e-6,
    gnd_thickness: 35e-6, epsilon_r: 4.4, tan_delta: 0.02,
    sigma_cond: 5.8e7, rq: 0, freq: 1e9,
};
const SOLVE_OPTS = { max_iters: 12, energy_tol: 0.005, param_tol: 0.02,
    max_nodes: 40000, min_converged_passes: 2, certify: false };

const silence = async fn => {
    const log = console.log;
    console.log = () => {};
    try { return await fn(); } finally { console.log = log; }
};

async function solveFdm(extra, max_nodes) {
    const s = new MicrostripSolver({ ...GEOM, mesh_backend: 'rectilinear', ...extra });
    const base = await silence(() => s.solve_adaptive({ ...SOLVE_OPTS, max_nodes }));
    const r = await silence(() => s.computeAtFrequency(1e9, base));
    return r.modes[0].RLGC.R;
}

async function main() {
    console.log('='.repeat(60));
    console.log('QS CONDUCTOR-LOSS INTEGRAND TEST (default microstrip @ 1 GHz)');
    console.log('='.repeat(60));
    let pass = true;
    const check = (ok, msg) => {
        console.log(`  ${ok ? '✓' : '✗'} ${msg}`);
        if (!ok) pass = false;
    };

    // 1. Mesh stability (the dielectric-field integrand climbed without bound here)
    const coarse = await solveFdm({ nx: 150, ny: 150 }, 70000);
    const dense = await solveFdm({ nx: 300, ny: 300 }, 250000);
    const drift = Math.abs(dense / coarse - 1);
    console.log(`\nFDM R: coarse=${coarse.toFixed(3)} dense=${dense.toFixed(3)} Ω/m (drift ${(drift * 100).toFixed(2)}%)`);
    check(drift < 0.025, `R mesh-stable (drift ${(drift * 100).toFixed(2)}% < 2.5%)`);

    // 2. Accuracy vs tri-backend MQS (matches HFSS within ~0.5%)
    const tri = new MicrostripSolver({ ...GEOM, mesh_backend: 'triangular' });
    const triBase = await silence(() => tri.solve_adaptive(SOLVE_OPTS));
    const triR = (await silence(() => tri.computeAtFrequency(1e9, triBase))).modes[0].RLGC.R;
    const err = dense / triR - 1;
    console.log(`tri MQS reference R=${triR.toFixed(3)} Ω/m; dense FDM error ${(err * 100).toFixed(2)}%`);
    check(Math.abs(err) < 0.08, `FDM R within 8% of tri MQS (${(err * 100).toFixed(1)}%)`);

    // 3. ε-invariance of R (quasi-TEM: loss follows current, not charge)
    const R22 = await solveFdm({ epsilon_r: 2.2, tan_delta: 0.001, nx: 150, ny: 150 }, 70000);
    const R98 = await solveFdm({ epsilon_r: 9.8, tan_delta: 0.001, nx: 150, ny: 150 }, 70000);
    const spread = Math.abs(R98 / R22 - 1);
    console.log(`R(εr=2.2)=${R22.toFixed(3)}  R(εr=9.8)=${R98.toFixed(3)}  spread ${(spread * 100).toFixed(3)}%`);
    check(spread < 0.01, `R is ε-invariant (spread ${(spread * 100).toFixed(2)}% < 1%)`);

    console.log('\n' + '='.repeat(60));
    console.log(pass ? '✓ ALL TESTS PASSED' : '✗ SOME TESTS FAILED');
    console.log('='.repeat(60));
    process.exit(pass ? 0 : 1);
}

main().catch(err => { console.error('✗ Test error:', err); process.exit(1); });
