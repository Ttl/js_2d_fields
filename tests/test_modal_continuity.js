// test_modal_continuity.js — regression tests for the general two-conductor (asymmetric)
// modal decomposition used by the differential solvers.
//
// The modal path replaces the odd/even drive results with the genuine eigenmodes of an
// asymmetric pair. Two properties must hold, and both guard against bugs that were once live:
//
//   1. CONTINUITY — as the asymmetry → 0, the modal modes must approach the odd/even modes
//      smoothly. A bug in the eigenvector normalisation (normalising to the max component)
//      made the modal Z0 JUMP ~11% for a 2% perturbation, because the rescaling is
//      discontinuous when the dominant component switches at symmetry. Fixed by |v|=√2.
//
//   2. DEGENERACY FALLBACK — for a homogeneous-εr (velocity-degenerate) pair
//      the eigenvectors are arbitrary, so the modal path must fall back to the
//      odd/even decomposition. A missing guard once produced >100% Z0 errors.
//
// Also checks that the two backends agree on a strongly-asymmetric case.
//
// Backend: `MESH_BACKEND=triangular node tests/test_modal_continuity.js`; default rectilinear.

import { BroadsideStriplineSolver } from '../src/broadside_stripline.js';

const MESH_BACKEND = process.env.MESH_BACKEND || 'rectilinear';
const SOLVE_OPTS = { max_iters: 8, max_nodes: 20000, min_converged_passes: 2, param_tol: 0.05, energy_tol: 0.01 };

// Inhomogeneous broadside (εr_bottom=5, εr_middle=3) so the two modes are non-degenerate and
// the modal path is meaningful; εr_top sets the dielectric asymmetry (εr_top=5 ⇒ symmetric).
function makeSolver(erTop, backend = MESH_BACKEND, extra = {}) {
    return new BroadsideStriplineSolver({
        trace_width: 0.2e-3, trace_thickness: 35e-6, x_offset: 0, sigma_cond: 5.8e7,
        h_bottom: 0.2e-3, er_bottom: 5, tand_bottom: 0.005,
        h_middle: 0.3e-3, er_middle: 3, tand_middle: 0.005,
        h_top: 0.2e-3, er_top: erTop, tand_top: 0.005,
        enclosure_width: 3e-3, freq: 5e9, nx: 30, ny: 30,
        boundaries: ['gnd', 'gnd', 'gnd', 'gnd'], mesh_backend: backend, ...extra,
    });
}

// Solve and return the two modes' eps_eff / Z0. `force` ('on'|'off'|undefined) toggles the
// modal path via the test hook in the guard so we can compare modal vs odd/even at the same point.
//
// Results are MEMOIZED per (geometry, backend, force): the sweep re-solves the εr_top
// 5.0 / 5.1 / 7.0 points that other tests also need, and each duplicated solve is a full
// adaptive mesh build (the triangular backend cannot use half-domain symmetry for a
// broadside pair, so these are its most expensive solves). Only RESULTS are cached — every
// actual solve still uses a fresh solver, because the tri backend evaluates the modal force
// flag once at build time (in _prepareStatic), so a mesh may never be reused across force
// values.
const _resultCache = new Map();
async function solveModes(solver, force) {
    const key = JSON.stringify([solver.mesh_backend, force ?? 'auto',
        solver.er_bottom, solver.er_middle, solver.er_top,
        solver.h_bottom, solver.h_middle, solver.h_top, solver.x_offset]);
    if (_resultCache.has(key)) return _resultCache.get(key);
    const result = await _solveModesUncached(solver, force);
    _resultCache.set(key, result);
    return result;
}
async function _solveModesUncached(solver, force) {
    globalThis.__MODAL_FORCE__ = force;
    solver.use_causal_materials = false;
    const log = console.log, warn = console.warn;
    console.log = () => {}; console.warn = () => {};
    try {
        const r = await solver.solve_adaptive(SOLVE_OPTS);
        return {
            eps: [r.modes[0].eps_eff, r.modes[1].eps_eff],
            Z0: [r.modes[0].Z0, r.modes[1].Z0],
        };
    } finally {
        console.log = log; console.warn = warn;
        globalThis.__MODAL_FORCE__ = undefined;
    }
}

const relDiff = (a, b) => Math.abs(a - b) / Math.max(Math.abs(a), Math.abs(b), 1e-30);
function assert(cond, msg) { if (!cond) throw new Error(msg); }

// --- Tests -----------------------------------------------------------------

// At exact symmetry the modal modes must equal the odd/even modes (they ARE odd/even). A small
// tolerance covers the per-conductor-vs-direct-drive discretisation difference only.
async function test_symmetric_matches_oddeven() {
    const modal = await solveModes(makeSolver(5.0), 'on');
    const oe = await solveModes(makeSolver(5.0), 'off');
    for (let k = 0; k < 2; k++) {
        assert(relDiff(modal.Z0[k], oe.Z0[k]) < 0.015,
            `symmetric mode ${k}: modal Z0=${modal.Z0[k].toFixed(2)} vs odd/even ${oe.Z0[k].toFixed(2)} (>1.5%)`);
        assert(relDiff(modal.eps[k], oe.eps[k]) < 0.015,
            `symmetric mode ${k}: modal eps=${modal.eps[k].toFixed(3)} vs odd/even ${oe.eps[k].toFixed(3)} (>1.5%)`);
    }
}

// CONTINUITY: a small (≈2%) dielectric perturbation must move the modal Z0 only slightly off the
// odd/even Z0 — the normalisation-discontinuity bug produced ~11% here, so require < 2%.
async function test_continuity_small_asymmetry() {
    const modal = await solveModes(makeSolver(5.1), 'on');   // εr_top 5→5.1, ~2% perturbation
    const oe = await solveModes(makeSolver(5.1), 'off');
    for (let k = 0; k < 2; k++) {
        const d = relDiff(modal.Z0[k], oe.Z0[k]);
        assert(d < 0.02,
            `slight asymmetry mode ${k}: modal Z0=${modal.Z0[k].toFixed(2)} deviates ${(d * 100).toFixed(1)}% ` +
            `from odd/even ${oe.Z0[k].toFixed(2)} (should be <2% — a larger jump means the eigenvector ` +
            `normalisation is discontinuous)`);
    }
}

// SMOOTHNESS: sweeping the asymmetry, the modal Z0 must vary smoothly — each step's change is
// bounded and the deviation from odd/even grows monotonically (no jumps).
async function test_smooth_divergence() {
    const sweep = [5.0, 5.1, 5.3, 5.6, 6.0, 7.0];
    const modal = [], oe = [];
    for (const et of sweep) { modal.push(await solveModes(makeSolver(et), 'on')); oe.push(await solveModes(makeSolver(et), 'off')); }
    for (let k = 0; k < 2; k++) {
        let prevDev = -1;
        for (let i = 0; i < sweep.length; i++) {
            const dev = relDiff(modal[i].Z0[k], oe[i].Z0[k]);
            if (i > 0) {
                const step = relDiff(modal[i].Z0[k], modal[i - 1].Z0[k]);
                assert(step < 0.10,
                    `mode ${k} at εr_top=${sweep[i]}: modal Z0 step ${(step * 100).toFixed(1)}% from previous ` +
                    `point is too large (geometry changed only slightly) — non-smooth modal Z0`);
                assert(dev >= prevDev - 0.005,
                    `mode ${k} at εr_top=${sweep[i]}: deviation from odd/even (${(dev * 100).toFixed(1)}%) ` +
                    `decreased as asymmetry grew — non-monotonic`);
            }
            prevDev = dev;
        }
    }
}

// DEGENERACY FALLBACK: a geometrically asymmetric but homogeneous-εr pair is velocity-degenerate,
// so the modal path must fall back to odd/even — forcing modal 'on' must give the same answer as
// odd/even (the guard returns the odd/even result regardless of the force flag when degenerate is
// NOT the path; here we verify the *default* path equals the odd/even path for the homogeneous case).
async function test_degenerate_falls_back() {
    // Homogeneous εr=4.4, asymmetric heights + offset.
    const geom = {
        er_bottom: 4.4, er_middle: 4.4, er_top: 4.4,
        h_bottom: 0.2e-3, h_middle: 0.235e-3, h_top: 0.165e-3, x_offset: 0.3e-3,
    };
    const def = await solveModes(makeSolver(4.4, MESH_BACKEND, geom), undefined);  // guard active → fallback
    const oe = await solveModes(makeSolver(4.4, MESH_BACKEND, geom), 'off');        // forced odd/even
    for (let k = 0; k < 2; k++) {
        assert(relDiff(def.Z0[k], oe.Z0[k]) < 1e-6,
            `homogeneous (degenerate) mode ${k}: default Z0=${def.Z0[k].toFixed(3)} did not fall back to ` +
            `odd/even ${oe.Z0[k].toFixed(3)} — the degeneracy guard is not engaging`);
    }
}

// BACKEND CONSISTENCY: a strongly-asymmetric inhomogeneous pair must give the same modes on both
// backends (independent implementations of the same modal decomposition).
async function test_backend_consistency() {
    const q = await solveModes(makeSolver(7.0, 'rectilinear'), 'on');
    const f = await solveModes(makeSolver(7.0, 'triangular'), 'on');
    for (let k = 0; k < 2; k++) {
        assert(relDiff(q.eps[k], f.eps[k]) < 0.03,
            `backend eps mode ${k}: QS=${q.eps[k].toFixed(3)} vs FW=${f.eps[k].toFixed(3)} (>3%)`);
        assert(relDiff(q.Z0[k], f.Z0[k]) < 0.03,
            `backend Z0 mode ${k}: QS=${q.Z0[k].toFixed(2)} vs FW=${f.Z0[k].toFixed(2)} (>3%)`);
    }
}

async function runTests() {
    const cases = [
        test_symmetric_matches_oddeven,
        test_continuity_small_asymmetry,
        test_smooth_divergence,
        test_degenerate_falls_back,
        test_backend_consistency,
    ];
    console.log(`\n### Modal continuity tests (mesh_backend = "${MESH_BACKEND}") ###`);
    let passed = 0;
    const failed = [];
    for (const c of cases) {
        try { await c(); passed++; console.log(`  ✓ ${c.name}`); }
        catch (e) { failed.push(`${c.name}: ${e.message}`); console.log(`  ✗ ${c.name}`); }
    }
    console.log(`\nSUMMARY [${MESH_BACKEND}]: ${passed}/${cases.length} passed, ${failed.length} failed`);
    failed.forEach(f => console.log(`  ✗ ${f}`));
    if (failed.length) process.exitCode = 1;
}

runTests();
