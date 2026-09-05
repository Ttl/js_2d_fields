// Tree-cotree gauge of the full-wave pencil (tri_fem.js buildTriGauge /
// assembleTriFEMDecomposed(..., gauge)) — the low-frequency stability fix.
//
// The plain Lee–Sun–Cendes pencil pairs curl-curl entries ~1/h² with O(k²) mass terms.
// The quasi-TEM eigenvector is nearly a pure gradient (S·∇ψ = 0), so once (k₀h)² reaches
// ~1e-14 on the smallest elements the mode's own Ritz residual sits at the 1e-4 gate and
// the solve returns nothing: on the 500 mm differential-microstrip mesh this happened at
// 100 MHz with 12k triangles (h_min 7e-8 m), and on the coarse ~1k-triangle mesh below
// it happens at 1 MHz — same mechanism, cheap to test. In the gauge the gradient part is
// an explicit nodal variable and the pencil is well-conditioned at any frequency.
//
// Checks, on enclosed differential (PEC and PMC symmetry planes), open microstrip and
// GCPW (several PEC components) meshes, for the closed real pencil, the complex lossy-ε
// pencil (the mode viewer's) and, structurally, the radiating-ABC pencil (Robin
// template; its pick converges nothing on these coarse meshes with or without the gauge):
//   1. the gauge is a congruence: yᵀX'y == (Ty)ᵀX(Ty) for every k²/Robin template,
//      and the curl-curl template vanishes exactly on the gradient subspace;
//   2. the gauged patterns are exactly symmetric (the eigensolver's LDLᵀ path);
//   3. gauged and plain eigensolves agree at 1 GHz, the eigenvector maps back onto the
//      static drive with the same overlap;
//   4. the gauged solve converges the quasi-TEM at 1 MHz and 100 kHz, where the plain
//      one does not, with the 1 GHz eigenvalue (dispersion ~3e-5);
//   5. end to end: TriBackend._eigenPick returns the quasi-TEM at 1 MHz.
//
// Run: node src/tri_solver/tests/gauge_test.mjs
import { MicrostripSolver } from '../../microstrip.js';
import { buildTriFreedomMap, assembleTriFEMDecomposed, femFromDecomposition, staticToEdgeDofs,
         buildTriGauge, gaugeExpand, gaugeSeed } from '../tri_fem.js';

const c0 = 299792458;
let failures = 0;
function check(name, cond, detail = '') {
    console.log(`${cond ? '✓ PASS' : '✗ FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
    if (!cond) failures++;
}
const silence = async fn => { const l = console.log; console.log = () => {}; try { return await fn(); } finally { console.log = l; } };
const matvec = (M, val, x) => {
    const y = new Float64Array(x.length);
    for (let i = 0; i < x.length; i++) { let a = 0; for (let p = M.rowPtr[i]; p < M.rowPtr[i+1]; p++) a += val[p] * x[M.colIdx[p]]; y[i] = a; }
    return y;
};
const dot = (a, b) => { let s = 0; for (let i = 0; i < a.length; i++) s += a[i] * b[i]; return s; };
const quad = (M, val, x) => dot(x, matvec(M, val, x));
const patternSymmetric = M => {
    const N = M.rowPtr.length - 1, set = new Set();
    for (let i = 0; i < N; i++) for (let p = M.rowPtr[i]; p < M.rowPtr[i+1]; p++) set.add(i * N + M.colIdx[p]);
    for (let i = 0; i < N; i++) for (let p = M.rowPtr[i]; p < M.rowPtr[i+1]; p++) if (!set.has(M.colIdx[p] * N + i)) return false;
    return true;
};
let rs = 20260905;
const rnd = () => { rs = (rs * 1103515245 + 12345) & 0x7fffffff; return rs / 0x7fffffff - 0.5; };

const CASES = {
    'diff-ms (enclosed)': {
        substrate_height: 0.2104e-3, trace_width: 0.35e-3, trace_thickness: 50e-6, gnd_thickness: 35e-6,
        epsilon_r: 4.4, tan_delta: 0.002, sigma_cond: 5.8e7, freq: 1e9, trace_spacing: 0.5e-3,
        enclosure_width: 3e-3, enclosure_height: 1.05e-3, boundaries: ['gnd', 'gnd', 'gnd', 'gnd'],
    },
    'microstrip (open)': {
        substrate_height: 0.254e-3, trace_width: 0.3e-3, trace_thickness: 35e-6, gnd_thickness: 35e-6,
        epsilon_r: 3.66, tan_delta: 0.003, sigma_cond: 5.8e7, freq: 5e9, boundaries: ['open', 'open', 'open', 'gnd'],
    },
    'gcpw (open)': {
        substrate_height: 0.2e-3, trace_width: 0.2e-3, trace_thickness: 35e-6, gnd_thickness: 35e-6,
        epsilon_r: 3.66, tan_delta: 0.003, sigma_cond: 5.8e7, freq: 5e9, use_coplanar_gnd: true, gap: 0.15e-3,
        boundaries: ['open', 'open', 'open', 'gnd'],
    },
};

// One quasi-TEM eigensolve of a pencil; returns { eps, x (original layout), ovl } or null.
function solveQuasiTEM(helpers, N, fem, f, epsStatic, seed, staticSeed, expand, nFT) {
    const k2 = (2 * Math.PI * f / c0) ** 2;
    const res = helpers.solveGeneralized(N, fem.csrA, fem.csrB, [-k2 * epsStatic, 0], 8, 17, seed);
    if (!res || res.nconv <= 0) return null;
    let best = null;
    for (let i = 0; i < res.nconv; i++) {
        if (res.evalsRe[i] >= 0) continue;
        const eps = -res.evalsRe[i] / k2;
        if (eps < 0.5 * epsStatic) continue;
        const x = expand(res.evecsRe.subarray(i * N, (i + 1) * N));
        let d = 0, nS = 0, nV = 0;
        for (let k = 0; k < nFT; k++) { d += staticSeed[k] * x[k]; nS += staticSeed[k] ** 2; nV += x[k] ** 2; }
        const ovl = Math.abs(d) / Math.sqrt(nS * nV);
        if (!best || ovl > best.ovl) best = { eps, x, ovl };
    }
    return best;
}

for (const [name, geo] of Object.entries(CASES)) {
    console.log(`\n=== ${name} ===`);
    const s = new MicrostripSolver({ ...geo, nx: 10, ny: 10, mesh_backend: 'triangular' });
    s.use_causal_materials = false;
    const tri = await silence(() => s._ensureTriBackend(null, { maxNodes: 8000, refineTol: 1e-3, maxRefineIters: 3, certify: false }));
    const mesh = tri.mesh, helpers = tri.ctx.helpers;
    console.log(`  mesh: ${mesh.nTris} triangles`);
    const lossyEps = mesh.epsMap.map(e => ({ re: e.re, im: -0.02 * e.re }));
    for (const mode of tri.modeNames) {
        const st = tri._static[mode];
        const variants = [{ label: 'closed', abc: st.abc, epsMap: mesh.epsMap, solve: true },
                          { label: 'closed, lossy ε', abc: st.abc, epsMap: lossyEps, solve: true }];
        if (Object.values(st.abc).some(v => v === 'pmc'))
            variants.push({ label: 'radiating ABC', solve: false, epsMap: mesh.epsMap,
                            abc: Object.fromEntries(Object.entries(st.abc).map(([k, v]) => [k, v === 'pmc' ? true : v])) });
        let eps1GHz = null;
        for (const { label, abc, epsMap, solve } of variants) {
            const tag = `${name} ${mode} [${label}]`;
            const fm = buildTriFreedomMap(mesh, tri.condRect, abc);
            const N = fm.nFreeTransverse + fm.nFreeLongitudinal, nFT = fm.nFreeTransverse;
            const gauge = buildTriGauge(mesh, fm);
            const decP = assembleTriFEMDecomposed(mesh, fm, epsMap, abc, tri.condRect);
            const decG = assembleTriFEMDecomposed(mesh, fm, epsMap, abc, tri.condRect, gauge);
            // 1. congruence
            const y = new Float64Array(N); for (let i = 0; i < N; i++) y[i] = rnd();
            const x = gaugeExpand(gauge, y);
            let worst = 0;
            for (const [M, Mg, tpl] of [[decP.A, decG.A, 'v1Re'], [decP.A, decG.A, 'v1Im'], [decP.A, decG.A, 'vrIm'],
                                        [decP.B, decG.B, 'v0Re'], [decP.B, decG.B, 'v1Re'], [decP.B, decG.B, 'v1Im']]) {
                const a = quad(M, M[tpl], x), b = quad(Mg, Mg[tpl], y);
                if (a === 0 && b === 0) continue;
                worst = Math.max(worst, Math.abs(a - b) / Math.abs(a));
            }
            check(`${tag}: congruence of the k²/Robin templates`, worst < 1e-10, `worst rel ${worst.toExponential(2)}`);
            const yg = new Float64Array(N); for (let i = gauge.nC; i < nFT; i++) yg[i] = rnd();
            const xg = gaugeExpand(gauge, yg);
            const sGrad = quad(decG.A, decG.A.v0Re, yg), sPlain = quad(decP.A, decP.A.v0Re, xg), mPlain = -quad(decP.A, decP.A.v1Re, xg);
            check(`${tag}: curl-curl vanishes exactly on the gradient subspace`, sGrad === 0,
                  `gauged ${sGrad}; plain roundoff ${(Math.abs(sPlain) / mPlain).toExponential(1)} of the mass term`);
            const cc = quad(decP.A, decP.A.v0Re, x), cg = quad(decG.A, decG.A.v0Re, y);
            check(`${tag}: curl-curl energy preserved`, Math.abs(cc - cg) < 1e-10 * Math.abs(cc), `rel ${(Math.abs(cc - cg) / Math.abs(cc)).toExponential(2)}`);
            // 2. symmetric patterns
            check(`${tag}: gauged patterns symmetric`, patternSymmetric(decG.A) && patternSymmetric(decG.B));
            if (!solve) continue;
            // 3./4. eigensolves
            const staticSeed = staticToEdgeDofs(st.phiEps, mesh, fm);
            const seedG = gaugeSeed(gauge, st.phiEps, mesh, fm);
            const solveP = f => solveQuasiTEM(helpers, N, femFromDecomposition(decP, (2 * Math.PI * f / c0) ** 2), f, st.eps_eff_static, staticSeed, staticSeed, v => v, nFT);
            const solveG = f => solveQuasiTEM(helpers, N, femFromDecomposition(decG, (2 * Math.PI * f / c0) ** 2), f, st.eps_eff_static, seedG, staticSeed, v => gaugeExpand(gauge, v), nFT);
            const p1 = solveP(1e9), g1 = solveG(1e9);
            check(`${tag}: plain and gauged agree at 1 GHz`, p1 && g1 && Math.abs(p1.eps - g1.eps) < 1e-6 * p1.eps,
                  p1 && g1 ? `eps ${p1.eps.toFixed(7)} vs ${g1.eps.toFixed(7)}, overlap ${p1.ovl.toFixed(3)} vs ${g1.ovl.toFixed(3)}` : 'a solve returned nothing');
            if (label === 'closed' && g1) eps1GHz = g1.eps;
            check(`${tag}: mapped eigenvector overlaps the static drive`, g1 && p1 && g1.ovl > 0.8 && Math.abs(g1.ovl - p1.ovl) < 1e-3);
            for (const f of [1e6, 1e5]) {
                const g = solveG(f);
                check(`${tag}: gauged converges at ${f / 1e6} MHz`, g && g1 && Math.abs(g.eps / g1.eps - 1) < 1e-3 && g.ovl > 0.8,
                      g ? `eps ${g.eps.toFixed(7)} (1 GHz ${g1.eps.toFixed(7)}), overlap ${g.ovl.toFixed(3)}` : 'no converged quasi-TEM');
            }
            if (label === 'closed' && mode === tri.modeNames[0]) {
                const p = solveP(1e6);
                console.log(`  (plain pencil at 1 MHz: ${p ? `converged eps ${p.eps.toFixed(6)}` : 'no converged quasi-TEM — the failure the gauge removes'})`);
            }
        }
        // 5. end to end through the backend's pick at 1 MHz (its closed pick, the one the
        // sweep and the static anchor use)
        const { fw, fwErr } = tri._eigenPick(st, 1e6, st.phiEps, st.eps_eff_static);
        check(`${name} ${mode}: TriBackend._eigenPick converges at 1 MHz`, !!fw && eps1GHz !== null && Math.abs(fw.eps / eps1GHz - 1) < 1e-3 && fw.bestOvl > 0.8,
              fw ? `eps ${fw.eps.toFixed(6)} (1 GHz ${eps1GHz?.toFixed(6)}, static ${st.eps_eff_static.toFixed(4)}), overlap ${fw.bestOvl.toFixed(3)}` : (fwErr ? String(fwErr.message || fwErr) : 'null'));
    }
}

console.log(failures ? `\n${failures} check(s) FAILED` : '\nAll checks passed');
process.exit(failures ? 1 : 0);
