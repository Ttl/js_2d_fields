// Phase 3 prototype: full-wave eigenmode solve on the generalized mesh.
// static (Z0, eps_eff_static) + eigensolve (eps_eff_mode dispersion) + analyzeTriMode.
import createModule from '../eigen_solver.js';
import { createWasmHelpers } from '../fem_core.js';
import { MicrostripSolver } from '../../microstrip.js';
import { initGmsh } from '../gmsh_mesh.js';
import { buildOccMeshFromGeometry } from '../occ_to_mesh.js';
import { buildTriFreedomMap, solveTriStatic, computeTriEnergy, staticToEdgeDofs, assembleTriFEM } from '../tri_fem.js';
import { analyzeTriMode } from '../tri_ms_solver.js';

const M = await createModule();
const { solveGeneralized } = createWasmHelpers(M);
const G = await initGmsh();
const c0 = 299792458, eps0 = 8.854187817e-12;

function meshOf(solver) {
    const dom = { x_min: -solver.domain_width / 2, x_max: solver.domain_width / 2,
                  y_min: -solver.t_gnd, y_max: solver.domain_height };
    return buildOccMeshFromGeometry(G, {
        conductors: solver.conductors, dielectrics: solver.dielectrics,
        domain: dom, boundaries: solver.boundaries,
        hFine: Math.min(Math.abs(solver.t), solver.w / 4) / 3,
        hCoarse: (dom.y_max - dom.y_min) / 18,
    });
}

function fullwave(solver, f, condPot = [1.0]) {
    const mesh = meshOf(solver);
    const fm = buildTriFreedomMap(mesh, mesh.condRect, {});
    // static
    const phiEps = solveTriStatic(mesh, fm, mesh.epsMap, condPot);
    const phiAir = solveTriStatic(mesh, fm, null, condPot);
    const W_eps = computeTriEnergy(phiEps, mesh, mesh.epsMap);
    const W_air = computeTriEnergy(phiAir, mesh, null);
    const C = eps0 * 2 * W_eps, L = 1 / (c0 * c0 * eps0 * 2 * W_air);
    const Z_static = Math.sqrt(L / C), eps_eff_static = W_eps / W_air;

    // eigensolve
    const k2 = (2 * Math.PI * f / c0) ** 2;
    const fem = assembleTriFEM(mesh, fm, k2, mesh.epsMap, {}, mesh.condRect, null);
    const N = fem.N;
    const seed = staticToEdgeDofs(phiEps, mesh, fm);
    const nev = 6, ncv = Math.min(2 * nev + 1, N - 1);
    const shift = [-k2 * eps_eff_static, 0];
    const res = solveGeneralized(N, fem.csrA, fem.csrB, shift, nev, ncv, seed);
    if (!res.nconv) return { Z_static, eps_eff_static, err: 'no eigen convergence' };
    // pick quasi-TEM mode = max overlap with static seed (transverse part)
    let bestIdx = 0, bestOvl = -1;
    for (let i = 0; i < res.nconv; i++) {
        const vR = res.evecsRe.slice(i * N, (i + 1) * N);
        let dot = 0, nS = 0, nV = 0;
        for (let k = 0; k < fm.nFreeTransverse; k++) { dot += seed[k] * vR[k]; nS += seed[k] ** 2; nV += vR[k] ** 2; }
        const ovl = nS > 0 && nV > 0 ? Math.abs(dot) / Math.sqrt(nS * nV) : 0;
        if (ovl > bestOvl) { bestOvl = ovl; bestIdx = i; }
    }
    const vRe = res.evecsRe.slice(bestIdx * N, (bestIdx + 1) * N);
    const vIm = res.evecsIm.slice(bestIdx * N, (bestIdx + 1) * N);
    const g2Re = res.evalsRe[bestIdx], g2Im = res.evalsIm[bestIdx];
    const eps_eff_mode = -g2Re / k2;
    const Z_freq = Z_static * Math.sqrt(eps_eff_static / eps_eff_mode);
    const an = analyzeTriMode(vRe, vIm, g2Re, g2Im, mesh, fm, k2, f, mesh.epsMap, mesh.condRect);
    return { Z_static, eps_eff_static, eps_eff_mode, Z_freq, ovl: bestOvl, Zpv: an.Zpv, ZpiLineAvg: an.ZpiLineAvg, nconv: res.nconv };
}

const ms = new MicrostripSolver({ trace_width: 1.5e-3, substrate_height: 0.5e-3, trace_thickness: 35e-6, epsilon_r: 4.4, tan_delta: 0.02 });
for (const f of [1e9, 10e9, 30e9]) {
    const r = fullwave(ms, f);
    console.log(`f=${(f/1e9).toFixed(0).padStart(2)}GHz  Z_static=${r.Z_static?.toFixed(2)}  eps_static=${r.eps_eff_static?.toFixed(3)}  eps_mode=${r.eps_eff_mode?.toFixed(3)}  Z(f)=${r.Z_freq?.toFixed(2)}  Zpv=${r.Zpv?.toFixed(2)}  ovl=${r.ovl?.toFixed(3)}  nconv=${r.nconv}${r.err?'  ERR:'+r.err:''}`);
}
console.log('\nFULLWAVE DONE');
