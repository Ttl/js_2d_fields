// Regression tests for full-wave (tri) backend correctness fixes:
//   1. R_dc convention parity with the FDM backend (signal + ground return in
//      series; coplanar GCPW grounds are NOT parallel signal metal).
//   2. R(f=0) = R_dc (used to return 0).
//   3. Differential R from the static perturbation path is per-trace (used to be
//      ~2× high — pair loss divided by single-trace power).
//   4. Plating with every face unchecked is bare metal (effectiveSurface used to
//      swap the whole conductor to the plating σ on the mere presence of plating).
//   5. An eigensolve failure at a sweep point surfaces a user-visible warning and
//      falls back to the quasi-static ε instead of failing silently.
//   6. The Robin ABC is applied on the bottom wall too (abc.bottom used to be
//      silently ignored → PMC).
//   7. An air edge bridging TWO different conductors across a one-element gap is
//      NOT marked PEC (isCondEdge used to key on endpoints only → gap shorted).
import { MicrostripSolver } from '../src/microstrip.js';
import { initTriBackend, TriBackend } from '../src/tri_solver/tri_backend.js';
import { buildTriFreedomMap, assembleTriFEM } from '../src/tri_solver/tri_fem.js';

const ctx = await initTriBackend();
let failures = 0;
function check(name, ok, detail = '') {
    console.log(`${ok ? '✓ PASS' : '✗ FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
    if (!ok) failures++;
}
const relDiff = (a, b) => Math.abs(a - b) / Math.max(Math.abs(b), 1e-30);

async function tri(opts, triOpts = {}) {
    const s = new MicrostripSolver(opts);
    const b = new TriBackend(ctx, s, { maxNodes: 12000, ...triOpts });
    await b.buildMesh();
    return { s, b };
}
async function fdm(opts) {
    const s = new MicrostripSolver(opts);
    s.ensure_mesh();
    const cached = await s.solve_adaptive({ max_iters: 6, energy_tol: 0.005, param_tol: 0.02, max_nodes: 40000 });
    return { s, cached };
}
function expectedRdc(s) {
    let sig = 0, gnd = 0;
    for (const c of s.conductors) {
        const a = Math.abs(c.width * c.height);
        if (c.is_signal) sig += a; else gnd += a;
    }
    return (sig > 0 ? 1 / (s.sigma_cond * sig) : 0) + (gnd > 0 ? 1 / (s.sigma_cond * gnd) : 0);
}

// ---- 1+2: R_dc — GCPW-like geometry (coplanar ground strips) ----
{
    const gcpwOpts = {
        trace_width: 0.3e-3, substrate_height: 0.2e-3, trace_thickness: 35e-6,
        epsilon_r: 3.66, tan_delta: 0.003, sigma_cond: 5.8e7, rq: 0,
        use_coplanar_gnd: true, gap: 0.2e-3, via_gap: 1.0e-3, use_vias: true,
        gnd_thickness: 35e-6, freq: 1e5,
    };
    const { s, b } = await tri(gcpwOpts);
    const hasCoplanarGnd = s.conductors.some(c => !c.is_signal && c.width < s.domain_width * 0.99);
    check('GCPW geometry has coplanar ground rects', hasCoplanarGnd);

    // f = 0: exact geometric DC resistance, FDM convention (series ground return).
    const r0 = b.solveAt(0).modes[0];
    const rdcExp = expectedRdc(s);
    check('R(f=0) = R_dc (signal + ground series)', relDiff(r0.RLGC.R, rdcExp) < 1e-9,
        `tri ${r0.RLGC.R.toFixed(4)} vs expected ${rdcExp.toFixed(4)} Ω/m`);
    check('R(f=0) > 0', r0.RLGC.R > 0);

    // Low frequency (δ ≫ t): R ≈ R_dc — must match the FDM backend's value.
    const rLow = b.solveAt(1e5).modes[0];
    const { s: sf, cached } = await fdm(gcpwOpts);
    const rf = (await sf.computeAtFrequency(1e5, cached)).modes[0];
    check('GCPW low-f R matches FDM (was ~10× low: grounds counted as signal metal)',
        relDiff(rLow.RLGC.R, rf.RLGC.R) < 0.25,
        `tri ${rLow.RLGC.R.toFixed(4)} vs fdm ${rf.RLGC.R.toFixed(4)} Ω/m`);
}

// ---- 3: differential R from the static perturbation path ----
{
    const diffOpts = {
        trace_width: 0.3e-3, trace_spacing: 0.2e-3, substrate_height: 0.2e-3,
        trace_thickness: 35e-6, epsilon_r: 3.66, tan_delta: 0.003,
        sigma_cond: 5.8e7, rq: 0, freq: 1e9,
    };
    const { b } = await tri(diffOpts, { lossMethod: 'static' });
    const rs = b.solveAt(1e9).modes;
    const { s: sf, cached } = await fdm(diffOpts);
    const rf = (await sf.computeAtFrequency(1e9, cached)).modes;
    for (const mode of ['odd', 'even']) {
        const t = rs.find(m => m.mode === mode), f = rf.find(m => m.mode === mode);
        check(`differential ${mode} static-path R is per-trace (was 2× high)`,
            relDiff(t.RLGC.R, f.RLGC.R) < 0.30,
            `tri ${t.RLGC.R.toFixed(2)} vs fdm ${f.RLGC.R.toFixed(2)} Ω/m`);
    }
}

// ---- 4: plating with all faces unchecked == bare metal ----
{
    const base = {
        trace_width: 0.3e-3, substrate_height: 0.2e-3, trace_thickness: 35e-6,
        epsilon_r: 3.66, tan_delta: 0.003, sigma_cond: 5.8e7, rq: 0, freq: 10e9,
    };
    const { b: bBare } = await tri(base);
    const { b: bNoFace } = await tri({ ...base,
        plating: { sigma: 8.7e6, thickness: 4e-6, rq: 0.3e-6, top: false, sides: false, bottom: false } });
    const { b: bPlated } = await tri({ ...base,
        plating: { sigma: 8.7e6, thickness: 4e-6, rq: 0.3e-6, top: true, sides: true, bottom: true } });
    const aBare = bBare.solveAt(10e9).modes[0].alpha_c;
    const aNoFace = bNoFace.solveAt(10e9).modes[0].alpha_c;
    const aPlated = bPlated.solveAt(10e9).modes[0].alpha_c;
    check('plating with no faces selected == bare metal', relDiff(aNoFace, aBare) < 1e-6,
        `no-face ${aNoFace.toFixed(4)} vs bare ${aBare.toFixed(4)} dB/m`);
    check('plating with faces selected increases loss (lower σ)', aPlated > aBare * 1.2,
        `plated ${aPlated.toFixed(4)} vs bare ${aBare.toFixed(4)} dB/m`);
}

// ---- 5: eigensolve failure → warning + static-ε fallback (not silence) ----
{
    const base = {
        trace_width: 0.3e-3, substrate_height: 0.2e-3, trace_thickness: 35e-6,
        epsilon_r: 3.66, tan_delta: 0.003, sigma_cond: 5.8e7, rq: 0, freq: 10e9,
    };
    const s = new MicrostripSolver(base);
    const brokenCtx = { ...ctx, helpers: { ...ctx.helpers,
        solveGeneralized: () => { throw new Error('injected eigensolver failure'); } } };
    const b = new TriBackend(brokenCtx, s, { maxNodes: 12000 });
    await b.buildMesh();   // refinement falls back to the static metric internally
    const res = b.solveAt(10e9);
    const warn = (res.warnings || []).find(w => w.type === 'eigensolve');
    check('eigensolve failure produces a user-visible warning', !!warn,
        warn ? warn.message.slice(0, 60) + '…' : 'no warning');
    const eps_static = b._static[b.modeNames[0]].eps_eff_static;
    check('eigensolve failure falls back to quasi-static ε (finite result)',
        isFinite(res.modes[0].eps_eff_mode) && relDiff(res.modes[0].eps_eff_mode, eps_static) < 1e-9,
        `eps_mode ${res.modes[0].eps_eff_mode.toFixed(4)} vs static ${eps_static.toFixed(4)}`);
}

// ---- 6: bottom-wall Robin ABC is assembled ----
{
    const base = {
        trace_width: 0.3e-3, substrate_height: 0.2e-3, trace_thickness: 35e-6,
        epsilon_r: 3.66, tan_delta: 0.003, sigma_cond: 5.8e7, rq: 0, freq: 10e9,
    };
    const { b } = await tri(base, { symmetry: false });
    const k2 = (2 * Math.PI * 10e9 / 299792458) ** 2;
    const sumAbsIm = (csr) => { let sum = 0; for (const v of csr.valIm) sum += Math.abs(v); return sum; };
    const fmPmc = buildTriFreedomMap(b.mesh, b.condRect, { bottom: 'pmc' });
    const femPmc = assembleTriFEM(b.mesh, fmPmc, k2, b.mesh.epsMap, { bottom: 'pmc' }, b.condRect, null);
    const femAbc = assembleTriFEM(b.mesh, fmPmc, k2, b.mesh.epsMap, { bottom: true }, b.condRect, null);
    check('abc.bottom = true adds the Robin (imaginary) boundary term',
        sumAbsIm(femAbc.csrA) > sumAbsIm(femPmc.csrA) + 1e-12,
        `|Im| abc ${sumAbsIm(femAbc.csrA).toExponential(2)} vs pmc ${sumAbsIm(femPmc.csrA).toExponential(2)}`);
}

// ---- 7: air edge bridging two conductors must stay FREE (not PEC) ----
{
    // Synthetic 2-triangle mesh: n0 sits on conductor A, n1 on conductor B, and
    // edge (n0,n1) crosses the air gap between them (midpoint on neither rect).
    const mesh = {
        nodes: new Float64Array([0, 0, 1, 0, 0.5, 1, 0.5, -1]),
        tris: new Int32Array([0, 1, 2, 0, 3, 1]),
        edges: new Int32Array([0, 1, 1, 2, 2, 0, 0, 3, 3, 1]),
        nNodes: 4, nTris: 2, nEdges: 5,
    };
    const condRect = {
        rects: [
            { xmin: -0.1, xmax: 0.1, ymin: -0.1, ymax: 0.1 },   // conductor A (contains n0)
            { xmin: 0.9, xmax: 1.1, ymin: -0.1, ymax: 0.1 },    // conductor B (contains n1)
        ],
        xmin_domain: -2, xmax_domain: 2, ymin_domain: -2, ymax_domain: 2,
    };
    const fm = buildTriFreedomMap(mesh, condRect, { bottom: true });
    check('bridging air edge is not a conductor edge', fm.isCondEdge[0] === 0);
    check('bridging air edge keeps its transverse DOFs', fm.edgeF[0] >= 0 && fm.edgeF[1] >= 0);
    check('bridging air edge keeps its longitudinal midpoint DOF', fm.edgeNodeF[0] >= 0);
    check('conductor groups are distinct on the two endpoints',
        fm.condNodeGroup[0] === 1 && fm.condNodeGroup[1] === 2);
}

console.log(failures === 0 ? '\nALL CORRECTNESS TESTS PASSED' : `\n${failures} TEST(S) FAILED`);
process.exit(failures === 0 ? 0 : 1);
