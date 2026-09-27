// The full-wave mesh follows every axis-aligned dielectric face, thin layers
// included: no triangle crosses one (the refiner's smoothing and edge swaps used to
// move the faces of layers thinner than two element sizes, and the centroid tagging
// then gave triangles across an interface one material). Checked on the main and the
// Modes-tab mesh of a slotline on a thin finite substrate, and on solder-masked lines,
// whose mask slivers must not trip the mesh-quality warning (Q > 100).
import { CustomGeometrySolver } from '../src/custom_geometry.js';
import { MicrostripSolver } from '../src/microstrip.js';

let failures = 0;
function check(name, ok, detail = '') {
    console.log(`${ok ? '✓ PASS' : '✗ FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
    if (!ok) failures++;
}
async function quiet(fn) {
    const log = console.log, warn = console.warn;
    console.log = () => {}; console.warn = () => {};
    try { return await fn(); } finally { console.log = log; console.warn = warn; }
}
const SOLVE = { max_iters: 12, energy_tol: 0.01, param_tol: 0.05, max_nodes: 20000, min_converged_passes: 2 };

// Triangles with vertices strictly on both sides of a dielectric rect's face, over the face's extent.
function crossing(mesh, dielectrics) {
    const { nodes, tris, nTris } = mesh, tol = 1e-12;
    const cut = (v, c) => Math.min(...v) < c - tol && Math.max(...v) > c + tol;
    let n = 0;
    for (let t = 0; t < nTris; t++) {
        const xs = [0, 1, 2].map(k => nodes[2 * tris[3 * t + k]]), ys = [0, 1, 2].map(k => nodes[2 * tris[3 * t + k] + 1]);
        if (dielectrics.some(d => !d.shape && (
            (Math.max(...ys) > d.y_min + tol && Math.min(...ys) < d.y_max - tol && (cut(xs, d.x_min) || cut(xs, d.x_max)))
            || (Math.max(...xs) > d.x_min + tol && Math.min(...xs) < d.x_max - tol && (cut(ys, d.y_min) || cut(ys, d.y_max)))))) n++;
    }
    return n;
}

const SLOT = 'units mm\nw = 1; s = 0.3; t = 0.05; h = 0.2104; wsub = 4\nbounds open open open open\n' +
    'diel x=-wsub/2 y=-h w=wsub h=h er=4.4 tand=0.02\nsig+ x=-s/2-w y=0 w=w h=t\ngnd x=s/2 y=0 w=w h=t\n';
{
    const s = new CustomGeometrySolver({ text: SLOT, sigma_cond: 5.8e7, freq: 1e9, nx: 10, ny: 10, mesh_backend: 'triangular' });
    s.tri_opts = { lossMethod: 'auto' };
    await quiet(() => s.solve_adaptive(SOLVE));
    const tb = s._triBackend;
    check('slotline on a thin substrate: main mesh follows the dielectric faces', crossing(tb.mesh, s.dielectrics) === 0,
        `${crossing(tb.mesh, s.dielectrics)} of ${tb.mesh.nTris} triangles cross one, max Q ${tb.meshQuality.maxQ.toFixed(1)}`);
    const sm = new CustomGeometrySolver({ text: SLOT, sigma_cond: 5.8e7, freq: 1e9, nx: 10, ny: 10 });
    await quiet(() => sm.solveModes(60e9, 12, null, { shrinkDomain: true, maxNodes: 20000, wavelengthDensity: 8 }));
    const mm = sm._modesBackend.mesh;
    check('slotline on a thin substrate: Modes-tab mesh follows the dielectric faces', crossing(mm, sm.dielectrics) === 0,
        `${crossing(mm, sm.dielectrics)} of ${mm.nTris} triangles cross one`);
}
const B = { trace_thickness: 35e-6, gnd_thickness: 35e-6, epsilon_r: 4.4, tan_delta: 0.02, sigma_cond: 5.8e7, freq: 1e9,
    nx: 10, ny: 10, substrate_height: 0.2e-3, trace_width: 0.35e-3, boundaries: ['open', 'open', 'open', 'gnd'],
    use_sm: true, mesh_backend: 'triangular' };
for (const [name, o] of [['solder-masked microstrip', B],
    ['solder-masked GCPW', { ...B, trace_width: 0.3e-3, use_coplanar_gnd: true, gap: 0.2e-3, via_gap: 0.5e-3, use_vias: true }]]) {
    const s = new MicrostripSolver(o);
    s.tri_opts = { lossMethod: 'auto' };
    await quiet(() => s.solve_adaptive(SOLVE));
    const tb = s._triBackend, n = crossing(tb.mesh, s.dielectrics);
    check(`${name}: mesh follows the mask faces, no quality warning`, n === 0 && tb.meshQuality.maxQ < 100,
        `${n} of ${tb.mesh.nTris} triangles cross a face, max Q ${tb.meshQuality.maxQ.toFixed(1)}`);
}

if (failures) { console.log(`\n${failures} CHECK(S) FAILED`); process.exit(1); }
console.log('\nALL CHECKS PASSED');
