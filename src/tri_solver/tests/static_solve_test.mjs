// Phase 2/3 validation: static FEM (C/L energy) on the generalized mesh.
// Z0 = sqrt(L/C), eps_eff = W_eps/W_air. Compares tri result to the FDM solver.
import { MicrostripSolver } from '../../microstrip.js';
import { initGmsh } from '../gmsh_mesh.js';
import { buildOccMeshFromGeometry } from '../occ_to_mesh.js';
import { buildTriFreedomMap, solveTriStatic, computeTriEnergy } from '../tri_fem.js';

const G = await initGmsh();
const c0 = 299792458, eps0 = 8.854187817e-12;

function triStatic(solver, condPotentials, hCoarseDiv = 16) {
    const dom = { x_min: -solver.domain_width / 2, x_max: solver.domain_width / 2,
                  y_min: -solver.t_gnd, y_max: solver.domain_height };
    const mesh = buildOccMeshFromGeometry(G, {
        conductors: solver.conductors, dielectrics: solver.dielectrics,
        domain: dom, boundaries: solver.boundaries,
        hFine: Math.min(solver.t, solver.w / 4) / 2,
        hCoarse: (dom.y_max - dom.y_min) / hCoarseDiv,
    });
    const fm = buildTriFreedomMap(mesh, mesh.condRect, {});
    const phiEps = solveTriStatic(mesh, fm, mesh.epsMap, condPotentials);
    const phiAir = solveTriStatic(mesh, fm, null, condPotentials);
    const W_eps = computeTriEnergy(phiEps, mesh, mesh.epsMap);
    const W_air = computeTriEnergy(phiAir, mesh, null);
    const C = eps0 * 2 * W_eps;
    const L = 1 / (c0 * c0 * eps0 * 2 * W_air);
    const Z0 = Math.sqrt(L / C);
    const eps_eff = W_eps / W_air;
    return { Z0, eps_eff, nTris: mesh.nTris };
}

function report(name, r) {
    console.log(`  ${name.padEnd(22)} Z0=${r.Z0.toFixed(2)} Ω  eps_eff=${r.eps_eff.toFixed(4)}  (nTris=${r.nTris})`);
}

console.log('=== Tri static FEM (generalized mesher) ===');
// Microstrip: w=1.5 h=0.5 er=4.4 → ~ Z0 50 Ω, eps_eff ~3.3 (hammerstad)
report('microstrip', triStatic(new MicrostripSolver({
    trace_width: 1.5e-3, substrate_height: 0.5e-3, trace_thickness: 35e-6,
    epsilon_r: 4.4, tan_delta: 0.02,
}), [1.0]));

// Stripline: signal between two grounds
report('stripline', triStatic(new MicrostripSolver({
    trace_width: 0.3e-3, substrate_height: 0.5e-3, trace_thickness: 35e-6,
    epsilon_r: 4.4, tan_delta: 0.02,
    epsilon_r_top: 4.4, enclosure_height: 0.5e-3, boundaries: ['open', 'open', 'gnd', 'gnd'],
}), [1.0]));

// GCPW: signal=1, two grounds=0
report('gcpw', triStatic(new MicrostripSolver({
    trace_width: 0.5e-3, substrate_height: 0.5e-3, trace_thickness: 35e-6,
    epsilon_r: 4.4, tan_delta: 0.02, use_coplanar_gnd: true, gap: 0.2e-3,
}), [1.0, 0.0, 0.0]));

// Differential — odd ([1,-1]) and even ([1,1])
const diff = new MicrostripSolver({
    trace_width: 0.35e-3, substrate_height: 0.21e-3, trace_thickness: 35e-6,
    trace_spacing: 0.1e-3, epsilon_r: 4.4, tan_delta: 0.02,
});
// rects order from mesher: by construction order (sig-1 then sig1)
const odd = triStatic(diff, [-1.0, 1.0]);
const even = triStatic(diff, [1.0, 1.0]);
report('diff odd', odd);
report('diff even', even);
console.log(`  Z_diff=${(2*odd.Z0).toFixed(2)} Ω  Z_common=${(even.Z0/2).toFixed(2)} Ω`);

console.log('\nSOLVE TEST DONE');
