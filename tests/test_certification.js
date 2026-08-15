// Verification certificate (adaptive refinement, BOTH backends).
//
// The refinement loop's pass-to-pass convergence gate measures the RATE of approach,
// not the DISTANCE to the converged answer — on the default UI microstrip it stopped
// with a reported "error" 5–20× below the actual one (C was 1.2% high at a 1% or even
// 0.5% tolerance setting). The certificate (tri_backend._certifyStatic for the
// triangular backend, FieldSolver2D._certifyStatic for the rectilinear one) fixes the
// semantics: a tripped gate is verified against a uniform-bisection solve, and the
// solve only stops when the Richardson-extrapolated remaining error passes the
// tolerance (or ends with an explicit accuracy warning when it can't). Both backends
// implement the same contract, so the UI Tolerance means the same thing on both.
//
// These tests pin the CONTRACT, reference-free, per backend:
//   1. Tolerance is honest: two solves at tolerances t1 > t2 must agree on C/eps_eff
//      within ~(t1 + t2) — if the gate were optimistic again, the coarse solve would
//      sit several t1 away from the fine one (the original bug: 1.2% apart at 0.5%).
//   2. The certificate is exposed and passing on a normal solve.
//   3. A solve that CANNOT meet the tolerance (tiny iteration/node budget) fails
//      certification and surfaces a type:'accuracy' warning instead of silence.
//   4. certify = false restores the legacy gate (no certification object).
//
// Run: node tests/test_certification.js
import { MicrostripSolver } from '../src/microstrip.js';

const GEOM = {
    trace_width: 0.35e-3, substrate_height: 0.21e-3, trace_thickness: 35e-6,
    epsilon_r: 4.4, tan_delta: 0.02, sigma_cond: 5.8e7, freq: 1e9, rq: 0,
    boundaries: ['open', 'open', 'open', 'gnd'], mesh_backend: 'triangular',
    nx: 30, ny: 30,
};

async function solveTol(tol, { maxIters = 12, maxNodes = 40000, triOpts = {} } = {}) {
    const s = new MicrostripSolver(GEOM);
    s.tri_opts = { lossMethod: 'auto', ...triOpts };
    const log = console.log;
    console.log = () => {};
    let r;
    try {
        r = await s.solve_adaptive({
            max_iters: maxIters, energy_tol: tol, max_nodes: maxNodes, min_converged_passes: 2,
        });
    } finally { console.log = log; }
    return { s, r, m: r.modes[0] };
}

let failures = 0;
function check(name, cond, detail = '') {
    console.log(`${cond ? '✓ PASS' : '✗ FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
    if (!cond) failures++;
}

// 1+2: honest tolerance + certificate exposed
const coarse = await solveTol(0.01);
const fine = await solveTol(0.001);
const dC = Math.abs(coarse.m.RLGC.C - fine.m.RLGC.C) / fine.m.RLGC.C;
const dEps = Math.abs(coarse.m.eps_eff - fine.m.eps_eff) / fine.m.eps_eff;
check('1% and 0.1% solves agree on C within 1.1%', dC < 0.011, `dC=${(100 * dC).toFixed(3)}%`);
check('1% and 0.1% solves agree on eps_eff within 1.1%', dEps < 0.011, `dEps=${(100 * dEps).toFixed(3)}%`);
check('certificate exposed on solver', !!coarse.s.certification, JSON.stringify(coarse.s.certification));
check('certificate passed at 1% tolerance', coarse.s.certification && coarse.s.certification.pass,
    coarse.s.certification ? `err=${(100 * coarse.s.certification.err).toFixed(3)}%` : 'missing');
check('certified error below the tolerance', coarse.s.certification && coarse.s.certification.err < 0.01);

// 3: unreachable tolerance → failed certificate + accuracy warning, not silence
const capped = await solveTol(0.0005, { maxIters: 3, maxNodes: 6000 });
const warn = (capped.r.warnings || capped.s.modeWarnings || []).find(w => w.type === 'accuracy');
check('capped solve fails certification', capped.s.certification && !capped.s.certification.pass,
    JSON.stringify(capped.s.certification));
check('capped solve surfaces an accuracy warning', !!warn, warn ? warn.message.slice(0, 60) + '…' : 'no warning');

// 4: opt-out restores the legacy gate
const legacy = await solveTol(0.01, { triOpts: { certify: false } });
check('certify:false leaves no certification', legacy.s.certification == null,
    JSON.stringify(legacy.s.certification ?? null));

// ---- Rectilinear (quasi-static) backend: the same contract ---------------------
// Coarse mesh hints + low frequency keep the base grid small enough that the
// certificate can afford its second bisection level (a genuine measured r).
const QS_GEOM = {
    ...GEOM, mesh_backend: 'rectilinear', freq: 1e8, nx: 40, ny: 40,
};

async function solveTolQS(tol, { maxIters = 12, maxNodes = 60000, extra = {} } = {}) {
    const s = new MicrostripSolver(QS_GEOM);
    const log = console.log;
    console.log = () => {};
    let r;
    try {
        r = await s.solve_adaptive({
            max_iters: maxIters, energy_tol: tol, param_tol: 0.05, max_nodes: maxNodes,
            min_converged_passes: 1, ...extra,
        });
    } finally { console.log = log; }
    return { s, r, m: r.modes[0] };
}

console.log('\nrectilinear backend:');

// 1+2: honest tolerance + certificate exposed
const qsCoarse = await solveTolQS(0.01);
const qsFine = await solveTolQS(0.005);
const qsDC = Math.abs(qsCoarse.m.C - qsFine.m.C) / qsFine.m.C;
const qsDEps = Math.abs(qsCoarse.m.eps_eff - qsFine.m.eps_eff) / qsFine.m.eps_eff;
check('1% and 0.5% solves agree on C within 1.5%', qsDC < 0.015, `dC=${(100 * qsDC).toFixed(3)}%`);
check('1% and 0.5% solves agree on eps_eff within 1.5%', qsDEps < 0.015, `dEps=${(100 * qsDEps).toFixed(3)}%`);
check('certificate exposed on solver', !!qsCoarse.s.certification, JSON.stringify(qsCoarse.s.certification));
check('certificate passed at 1% tolerance', qsCoarse.s.certification && qsCoarse.s.certification.pass,
    qsCoarse.s.certification ? `err=${(100 * qsCoarse.s.certification.err).toFixed(3)}%` : 'missing');
check('certified error below the tolerance', qsCoarse.s.certification && qsCoarse.s.certification.err < 0.01);

// 3: unreachable tolerance → failed certificate + accuracy warning, not silence
const qsCapped = await solveTolQS(0.001, { maxIters: 2, maxNodes: 3000 });
const qsWarn = (qsCapped.r.warnings || []).find(w => w.type === 'accuracy');
check('capped solve fails certification', qsCapped.s.certification && !qsCapped.s.certification.pass,
    JSON.stringify(qsCapped.s.certification));
check('capped solve surfaces an accuracy warning', !!qsWarn, qsWarn ? qsWarn.message.slice(0, 60) + '…' : 'no warning');

// 4: opt-out restores the legacy gate
const qsLegacy = await solveTolQS(0.01, { extra: { certify: false } });
check('certify:false leaves no certification', qsLegacy.s.certification == null,
    JSON.stringify(qsLegacy.s.certification ?? null));

console.log(failures ? `\n${failures} FAILURE(S)` : '\nAll certification tests passed.');
process.exit(failures ? 1 : 0);
