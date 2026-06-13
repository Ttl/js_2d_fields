// Phase 2 acceptance: cross-check tri static FEM vs the production FDM solver.
import { MicrostripSolver } from '../../microstrip.js';
import { initGmsh } from '../gmsh_mesh.js';
import { buildGmshMeshFromGeometry } from '../geom_to_mesh.js';
import { buildTriFreedomMap, solveTriStatic, computeTriEnergy } from '../tri_fem.js';

const G = await initGmsh();
const c0 = 299792458, eps0 = 8.854187817e-12;

function triStatic(solver, condPotentials) {
    const dom = { x_min: -solver.domain_width / 2, x_max: solver.domain_width / 2,
                  y_min: -solver.t_gnd, y_max: solver.domain_height };
    const mesh = buildGmshMeshFromGeometry(G, {
        conductors: solver.conductors, dielectrics: solver.dielectrics,
        domain: dom, boundaries: solver.boundaries,
        hFine: Math.min(solver.t, solver.w / 4) / 3,
        hCoarse: (dom.y_max - dom.y_min) / 20,
    });
    const fm = buildTriFreedomMap(mesh, mesh.condRect, {});
    const phiEps = solveTriStatic(mesh, fm, mesh.epsMap, condPotentials);
    const phiAir = solveTriStatic(mesh, fm, null, condPotentials);
    const W_eps = computeTriEnergy(phiEps, mesh, mesh.epsMap);
    const W_air = computeTriEnergy(phiAir, mesh, null);
    const C = eps0 * 2 * W_eps, L = 1 / (c0 * c0 * eps0 * 2 * W_air);
    return { Z0: Math.sqrt(L / C), eps_eff: W_eps / W_air };
}

async function fdm(opts) {
    const s = new MicrostripSolver(opts);
    s.freq = 1e9;
    s.ensure_mesh();
    const res = await s.solve_adaptive({ max_iters: 8, energy_tol: 0.005, param_tol: 0.02, max_nodes: 40000 });
    const m = res.modes[0];
    return { Z0: m.Z0, eps_eff: m.eps_eff };
}

function cmp(name, tri, ref) {
    const dZ = 100 * (tri.Z0 - ref.Z0) / ref.Z0;
    const dE = 100 * (tri.eps_eff - ref.eps_eff) / ref.eps_eff;
    console.log(`${name.padEnd(14)} tri: Z0=${tri.Z0.toFixed(2)} ee=${tri.eps_eff.toFixed(3)} | fdm: Z0=${ref.Z0.toFixed(2)} ee=${ref.eps_eff.toFixed(3)} | ΔZ0=${dZ.toFixed(1)}% Δee=${dE.toFixed(1)}%`);
}

const msOpts = { trace_width: 1.5e-3, substrate_height: 0.5e-3, trace_thickness: 35e-6, epsilon_r: 4.4, tan_delta: 0.02 };
cmp('microstrip', triStatic(new MicrostripSolver(msOpts), [1.0]), await fdm(msOpts));

const slOpts = { trace_width: 0.3e-3, substrate_height: 0.5e-3, trace_thickness: 35e-6, epsilon_r: 4.4, tan_delta: 0.02, epsilon_r_top: 4.4, enclosure_height: 0.5e-3, boundaries: ['open', 'open', 'gnd', 'gnd'] };
cmp('stripline', triStatic(new MicrostripSolver(slOpts), [1.0]), await fdm(slOpts));

const gcpwOpts = { trace_width: 0.5e-3, substrate_height: 0.5e-3, trace_thickness: 35e-6, epsilon_r: 4.4, tan_delta: 0.02, use_coplanar_gnd: true, gap: 0.2e-3 };
cmp('gcpw', triStatic(new MicrostripSolver(gcpwOpts), [1.0, 0.0, 0.0]), await fdm(gcpwOpts));

console.log('\nXCHECK DONE');
