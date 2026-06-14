// Phase 2 coverage: all geometry-affecting features mesh + static-solve sanely.
import { MicrostripSolver } from '../../microstrip.js';
import { BroadsideStriplineSolver } from '../../broadside_stripline.js';
import { initGmsh } from '../gmsh_mesh.js';
import { buildOccMeshFromGeometry } from '../occ_to_mesh.js';
import { buildTriFreedomMap, solveTriStatic, computeTriEnergy } from '../tri_fem.js';

const G = await initGmsh();
const c0 = 299792458, eps0 = 8.854187817e-12;

function meshOf(solver) {
    const dom = { x_min: -solver.domain_width / 2, x_max: solver.domain_width / 2,
                  y_min: -solver.t_gnd, y_max: solver.domain_height };
    return buildOccMeshFromGeometry(G, {
        conductors: solver.conductors, dielectrics: solver.dielectrics,
        domain: dom, boundaries: solver.boundaries,
        hFine: Math.min(Math.abs(solver.t || 35e-6), (solver.w || 0.3e-3) / 4) / 2.5,
        hCoarse: (dom.y_max - dom.y_min) / 16,
    });
}

// drive: signal rects → polarity value, ground rects → 0
function drivePotentials(mesh, mode = 'single') {
    return mesh.condRect.rectRoles.map(r => {
        if (!r.is_signal) return 0;
        if (mode === 'even') return 1;
        if (mode === 'odd') return r.polarity;   // ±1
        return 1;                                 // single
    });
}

function staticZ(mesh, mode) {
    const fm = buildTriFreedomMap(mesh, mesh.condRect, {});
    const pot = drivePotentials(mesh, mode);
    const phiEps = solveTriStatic(mesh, fm, mesh.epsMap, pot);
    const phiAir = solveTriStatic(mesh, fm, null, pot);
    const W_eps = computeTriEnergy(phiEps, mesh, mesh.epsMap);
    const W_air = computeTriEnergy(phiAir, mesh, null);
    const C = eps0 * 2 * W_eps, L = 1 / (c0 * c0 * eps0 * 2 * W_air);
    return { Z0: Math.sqrt(L / C), eps_eff: W_eps / W_air };
}

function checkMesh(name, mesh) {
    let bad = 0;
    for (let t = 0; t < mesh.nTris; t++) {
        const v0 = mesh.tris[3*t], v1 = mesh.tris[3*t+1], v2 = mesh.tris[3*t+2];
        const ax=mesh.nodes[2*v0],ay=mesh.nodes[2*v0+1],bx=mesh.nodes[2*v1],by=mesh.nodes[2*v1+1],cx=mesh.nodes[2*v2],cy=mesh.nodes[2*v2+1];
        if (0.5*Math.abs((bx-ax)*(cy-ay)-(cx-ax)*(by-ay)) < 1e-20) bad++;
    }
    let emin=9e9,emax=-9e9; for(const e of mesh.epsMap){emin=Math.min(emin,e.re);emax=Math.max(emax,e.re);}
    const roles = mesh.condRect.rectRoles.map(r => r.is_signal?`s${r.polarity}`:'g').join(',');
    console.log(`${name.padEnd(20)} nTris=${String(mesh.nTris).padStart(6)} cond=${mesh.condRect.rects.length} [${roles}] ε[${emin.toFixed(2)},${emax.toFixed(2)}] degen=${bad} wallPEC=${JSON.stringify(mesh.condRect.wallPEC)}`);
    return bad === 0;
}

// --- top dielectric coating ---
let m = meshOf(new MicrostripSolver({ trace_width: 1.5e-3, substrate_height: 0.5e-3, trace_thickness: 35e-6, epsilon_r: 4.4, tan_delta: 0.02, top_diel_h: 0.05e-3, top_diel_er: 3.5 }));
checkMesh('top_diel', m); console.log('   ', JSON.stringify(staticZ(m, 'single')));

// --- solder mask ---
m = meshOf(new MicrostripSolver({ trace_width: 1.5e-3, substrate_height: 0.5e-3, trace_thickness: 35e-6, epsilon_r: 4.4, tan_delta: 0.02, use_sm: true, sm_t_sub: 20e-6, sm_t_trace: 20e-6, sm_t_side: 20e-6, sm_er: 3.3, sm_tand: 0.02 }));
checkMesh('solder_mask', m); console.log('   ', JSON.stringify(staticZ(m, 'single')));

// --- gnd_cut (partial bottom grounds — not absorbed) ---
m = meshOf(new MicrostripSolver({ trace_width: 1.5e-3, substrate_height: 0.5e-3, trace_thickness: 35e-6, epsilon_r: 4.4, tan_delta: 0.02, gnd_cut_width: 3e-3, gnd_cut_sub_h: 0.3e-3 }));
checkMesh('gnd_cut', m); console.log('   ', JSON.stringify(staticZ(m, 'single')));

// --- enclosure (side + top grounds) ---
m = meshOf(new MicrostripSolver({ trace_width: 1.5e-3, substrate_height: 0.5e-3, trace_thickness: 35e-6, epsilon_r: 4.4, tan_delta: 0.02, enclosure_width: 6e-3, enclosure_height: 3e-3, boundaries: ['gnd','gnd','gnd','gnd'] }));
checkMesh('enclosure', m); console.log('   ', JSON.stringify(staticZ(m, 'single')));

// --- GCPW + vias ---
m = meshOf(new MicrostripSolver({ trace_width: 0.5e-3, substrate_height: 0.5e-3, trace_thickness: 35e-6, epsilon_r: 4.4, tan_delta: 0.02, use_coplanar_gnd: true, use_vias: true, gap: 0.2e-3, via_gap: 0.8e-3 }));
checkMesh('gcpw+vias', m); console.log('   ', JSON.stringify(staticZ(m, 'single')));

// --- broadside coupled stripline (asymmetric, embedded traces at different heights) ---
m = meshOf(new BroadsideStriplineSolver({ trace_width: 0.3e-3, trace_thickness: 35e-6, h_bottom: 0.2e-3, h_middle: 0.2e-3, h_top: 0.2e-3, er_bottom: 4.4, er_middle: 4.4, er_top: 4.4, tand_bottom: 0.02, tand_middle: 0.02, tand_top: 0.02, x_offset: 0.1e-3 }));
checkMesh('broadside', m);
console.log('    odd ', JSON.stringify(staticZ(m, 'odd')), ' even', JSON.stringify(staticZ(m, 'even')));

console.log('\nCOVERAGE DONE');
