// Coaxial line against CLOSED FORMS — the analytic-truth anchor for the shaped
// (non-rectangular) geometry path, in the same spirit as test_parallel_plate.js.
//
// A coax is the strongest validation target the suite has: every reported quantity has
// an exact closed form, and the geometry has NO corners, so it isolates formula and
// discretization error from corner-singularity error (which dominates every rectangular
// trace). It is also the geometry the SIBC saturation constant in conductor_loss.js was
// originally calibrated against.
//
//   Z0    = eta0/(2*pi*sqrt(er)) * ln(b/a)
//   C     = 2*pi*eps0*er / ln(b/a)          L_ext = mu0/(2*pi) * ln(b/a)
//   eps_eff (static) = er EXACTLY on any mesh — the homogeneous-fill invariant, since
//           the eps-weighted and vacuum solves differ only by a scalar on the stiffness
//           matrix and so have the identical potential.
//   alpha_d = pi*f*sqrt(er)/c0 * tand       (mesh-independent: C*Z0 = sqrt(er)/c0, so
//           the geometry cancels out of alpha_d = omega*tand*C*Z0/2 entirely)
//   alpha_c = Rs*(1/a + 1/b) / (2*eta*ln(b/a)),  Rs = 1/(sigma*delta), eta = eta0/sqrt(er)
//
// Frequencies are chosen to exercise BOTH conductor-loss estimators: 50 MHz is below
// F_STATIC_MAX so it takes the static-field perturbation path, while 1/10 GHz are deep
// in the skin regime (delta/2a ~ 2e-3) and take the SIBC-blended path.
//
// alpha_c is the one quantity that converges slowly, because it is a pure surface
// integral: it lands ~1% low on the default mesh and is the only thing that mesh size
// is chosen for (see CoaxSolver._hFine). Everything else — Z0, C, L, eps_eff, alpha_d —
// is at or below 0.05% here and stays under 0.1% on meshes several times coarser, so the
// 0.5% tolerances below are loose by design and the 2% on alpha_c is the real bound.
//
// Geometry construction, shape predicates and parameter validation live in
// tests/test_geometry.js (pure JS, milliseconds); this file owns the physics.
//
// Run: node tests/test_coax.js
import { CoaxSolver } from '../src/coax.js';

const EPS0 = 8.854187817e-12, MU0 = 4e-7 * Math.PI;
const C0 = 1 / Math.sqrt(EPS0 * MU0), ETA0 = Math.sqrt(MU0 / EPS0);
const NP_TO_DB = 8.685889638;

// RG-402-like semi-rigid: 0.92 mm centre conductor, 2.95 mm PTFE dielectric.
const D_IN = 0.92e-3, D_OUT = 2.95e-3, ER = 2.1, TAND = 2e-4, SIG = 5.8e7;
const a = D_IN / 2, b = D_OUT / 2;
const LNBA = Math.log(b / a), ETA = ETA0 / Math.sqrt(ER);

const Z0_x = ETA0 / (2 * Math.PI * Math.sqrt(ER)) * LNBA;
const C_x = 2 * Math.PI * EPS0 * ER / LNBA;
const L_x = MU0 * LNBA / (2 * Math.PI);
const alphaD_x = (f) => Math.PI * f * Math.sqrt(ER) / C0 * TAND * NP_TO_DB;
const alphaC_x = (f) => {
    const delta = Math.sqrt(2 / (2 * Math.PI * f * MU0 * SIG));
    return (1 / (SIG * delta)) * (1 / a + 1 / b) / (2 * ETA * LNBA) * NP_TO_DB;
};
// Series resistance implied by the analytic alpha_c, and the characteristic impedance
// including it. In the skin regime Zs = Rs(1+j), so L_int = R/omega — at 50 MHz that is
// a real +1.1% on L and lifts Re(Zc) 0.57% above the lossless Z0. Comparing Zc against
// Z0 alone would be wrong physics, not a tolerance problem.
const R_x = (f) => 2 * (alphaC_x(f) / NP_TO_DB) * Z0_x;
const Zc_x = (f) => Math.sqrt((L_x + R_x(f) / (2 * Math.PI * f)) / C_x);

let failures = 0;
function check(name, cond, detail = '') {
    console.log(`${cond ? '✓ PASS' : '✗ FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
    if (!cond) failures++;
}
const pct = (a_, b_) => ((a_ - b_) / b_ * 100).toFixed(3) + '%';
function near(name, got, want, tolPct) {
    check(name, Math.abs(got / want - 1) * 100 <= tolPct,
        `${got.toExponential(5)} vs ${want.toExponential(5)}, ${pct(got, want)}, tol ${tolPct}%`);
}

function makeSolver(extra = {}) {
    const s = new CoaxSolver({
        inner_diameter: D_IN, dielectric_diameter: D_OUT, epsilon_r: ER,
        tan_delta: TAND, sigma_cond: SIG, freq: 1e9, ...extra,
    });
    // Djordjevic-Sarkar would shift er away from the nominal value and invalidate every
    // closed form above.
    s.use_causal_materials = false;
    return s;
}

async function solveCoax(extra = {}, opts = {}) {
    const s = makeSolver(extra);
    const log = console.log; console.log = () => {};
    try {
        const first = await s.solve_adaptive({ max_nodes: 20000, ...opts });
        return { s, first };
    } finally { console.log = log; }
}
async function at(s, f) {
    const log = console.log; console.log = () => {};
    try { return (await s.computeAtFrequency(f)).modes[0]; } finally { console.log = log; }
}

// ---------------------------------------------------------------- closed forms
console.log('\n=== closed forms ===');
const { s, first } = await solveCoax();
console.log(`  mesh: ${s.triMesh.nTris} triangles`);
check('exactly one mode (single-ended)', first.modes.length === 1);
check('no solver warnings', !first.warnings || first.warnings.length === 0,
    JSON.stringify(first.warnings || []));

for (const f of [50e6, 1e9, 10e9]) {
    const m = await at(s, f);
    console.log(`\n-- ${(f / 1e9).toFixed(3)} GHz --`);
    near('Z0', m.Z0, Z0_x, 0.5);
    near('Zc.re', m.Zc.re, Zc_x(f), 0.5);
    near('R', m.RLGC.R, R_x(f), 2.0);
    near('C', m.RLGC.C, C_x, 0.5);
    near('L_external', m.L_external, L_x, 0.5);
    near('eps_eff_mode', m.eps_eff_mode, ER, 0.5);
    // Geometry-independent: alpha_d = omega*tand*C*Z0/2 and C*Z0 = sqrt(er)/c0.
    near('alpha_d', m.alpha_d, alphaD_x(f), 0.5);
    // Surface integral — the slowest-converging quantity (see header).
    near('alpha_c', m.alpha_c, alphaC_x(f), 2.0);
}

// The homogeneous-fill invariant, and the sharpest assertion in the file. Below
// F_STATIC_MAX the reported eps_eff_mode IS the static W_eps/W_air, and with a uniform
// epsMap the eps-weighted and vacuum stiffness matrices differ by a scalar — so that
// ratio is er to machine precision on ANY mesh, at any refinement level.
//
// It doubles as the material-tagging check: if even a few triangles were tagged air
// (the failure mode if the mesher and tagMaterials ever materialized different
// polygons, leaving vacuum slivers at the vertex "horns"), the ratio would drift off er.
{
    const m = await at(s, 50e6);
    near('eps_eff_mode is EXACTLY er (homogeneous fill / uniform tagging)', m.eps_eff_mode, ER, 1e-6);
}

// ---------------------------------------------------------------- symmetry
// A coax is mirror-symmetric about x=0, so the backend meshes only the x>=0 half and
// applies a PMC (natural) wall on the plane — worth ~2.6x on solve time. This is exact
// rather than approximate because n % 4 == 0 and phase 0 put polygon vertices precisely
// on the y axis, so the half really is half the n-gon.
//
// The full-domain reference below is deliberately run on a COARSER mesh: it is a
// different mesh either way, so agreement between them is a mesh-independent statement,
// and this keeps the fast tier affordable.
console.log('\n=== half-domain symmetry ===');
{
    check('the solve used the half domain', s._triBackend.symmetry === true);
    let minX = Infinity;
    for (let i = 0; i < s.triMesh.nodes.length; i += 2) minX = Math.min(minX, s.triMesh.nodes[i]);
    check('half-domain mesh is x >= 0', minX > -b * 1e-6, `min x = ${minX.toExponential(2)}`);

    const fullS = makeSolver();
    // Clearing the flag is exactly what isXSymmetric keys on, so this exercises the gate.
    for (const o of [...fullS.conductors, ...fullS.dielectrics]) o.shape.xSymmetric = false;
    fullS.tri_mesh_hints.hFine *= 2;
    const log = console.log; console.log = () => {};
    await fullS.solve_adaptive({ max_nodes: 20000 });
    console.log = log;
    check('the reference solve used the full domain', fullS._triBackend.symmetry === false);

    const h = await at(s, 1e9), fu = await at(fullS, 1e9);
    near('half vs full: Z0', h.Z0, fu.Z0, 1.0);
    near('half vs full: C', h.RLGC.C, fu.RLGC.C, 1.0);
    near('half vs full: L', h.RLGC.L, fu.RLGC.L, 1.0);
    near('half vs full: alpha_c', h.alpha_c, fu.alpha_c, 2.0);
    near('half vs full: alpha_d', h.alpha_d, fu.alpha_d, 1.0);
}

console.log(`\n${failures === 0 ? 'ALL TESTS PASSED' : failures + ' TEST(S) FAILED'}`);
process.exit(failures === 0 ? 0 : 1);
