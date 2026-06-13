// Phase 2 test: drive the generalized mesher from real app geometry.
import { MicrostripSolver } from '../../microstrip.js';
import { initGmsh } from '../gmsh_mesh.js';
import { buildGmshMeshFromGeometry } from '../geom_to_mesh.js';

const G = await initGmsh();

function domainOf(s) {
    return { x_min: -s.domain_width / 2, x_max: s.domain_width / 2, y_min: -s.t_gnd, y_max: s.domain_height };
}

function testCase(name, solver) {
    const dom = domainOf(solver);
    const h = solver.h;
    const mesh = buildGmshMeshFromGeometry(G, {
        conductors: solver.conductors,
        dielectrics: solver.dielectrics,
        domain: dom,
        boundaries: solver.boundaries,
        hFine: Math.min(solver.t, solver.w / 4) / 2,
        hCoarse: (dom.y_max - dom.y_min) / 12,
    });
    const cr = mesh.condRect;
    console.log(`\n=== ${name} ===`);
    console.log(`  domain raw: x[${dom.x_min.toExponential(2)},${dom.x_max.toExponential(2)}] y[${dom.y_min.toExponential(2)},${dom.y_max.toExponential(2)}]`);
    console.log(`  meshed:     x[${cr.xmin_domain.toExponential(2)},${cr.xmax_domain.toExponential(2)}] y[${cr.ymin_domain.toExponential(2)},${cr.ymax_domain.toExponential(2)}] wallPEC=${JSON.stringify(cr.wallPEC)}`);
    console.log(`  nNodes=${mesh.nNodes} nTris=${mesh.nTris} nEdges=${mesh.nEdges}`);
    console.log(`  conductors in: ${cr.rects.length} (roles: ${cr.rectRoles.map(r => r.is_signal ? `sig${r.polarity}` : 'gnd').join(',')})`);
    // sanity: epsMap range
    let emin = 1e9, emax = -1e9;
    for (const e of mesh.epsMap) { emin = Math.min(emin, e.re); emax = Math.max(emax, e.re); }
    console.log(`  epsMap range: [${emin.toFixed(3)}, ${emax.toFixed(3)}]  (substrate er=${solver.er})`);
    // sanity: every triangle non-degenerate
    let bad = 0, minArea = Infinity;
    for (let t = 0; t < mesh.nTris; t++) {
        const v0 = mesh.tris[3*t], v1 = mesh.tris[3*t+1], v2 = mesh.tris[3*t+2];
        const ax = mesh.nodes[2*v0], ay = mesh.nodes[2*v0+1];
        const bx = mesh.nodes[2*v1], by = mesh.nodes[2*v1+1];
        const cx = mesh.nodes[2*v2], cy = mesh.nodes[2*v2+1];
        const area = 0.5 * Math.abs((bx-ax)*(cy-ay) - (cx-ax)*(by-ay));
        if (area < 1e-20) bad++;
        minArea = Math.min(minArea, area);
    }
    console.log(`  degenerate tris: ${bad}, minArea=${minArea.toExponential(2)}`);
    return mesh;
}

// Plain microstrip (FR4)
testCase('microstrip', new MicrostripSolver({
    trace_width: 1.5e-3, substrate_height: 0.5e-3, trace_thickness: 35e-6,
    epsilon_r: 4.4, tan_delta: 0.02, sigma_cond: 5.8e7,
}));

// Differential microstrip
testCase('diff_microstrip', new MicrostripSolver({
    trace_width: 0.35e-3, substrate_height: 0.21e-3, trace_thickness: 35e-6,
    trace_spacing: 0.1e-3, epsilon_r: 4.4, tan_delta: 0.02,
}));

// Stripline (top + bottom ground)
testCase('stripline', new MicrostripSolver({
    trace_width: 0.3e-3, substrate_height: 0.5e-3, trace_thickness: 35e-6,
    epsilon_r: 4.4, tan_delta: 0.02,
    epsilon_r_top: 4.4, enclosure_height: 0.5e-3, boundaries: ['open', 'open', 'gnd', 'gnd'],
}));

// GCPW (coplanar grounds)
testCase('gcpw', new MicrostripSolver({
    trace_width: 0.5e-3, substrate_height: 0.5e-3, trace_thickness: 35e-6,
    epsilon_r: 4.4, tan_delta: 0.02, use_coplanar_gnd: true, gap: 0.2e-3,
}));

console.log('\nMESH TEST DONE');
