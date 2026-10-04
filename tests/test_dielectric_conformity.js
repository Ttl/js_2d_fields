// The full-wave mesh follows every axis-aligned dielectric face, thin layers
// included: no triangle crosses one (the refiner's smoothing and edge swaps used to
// move the faces of layers thinner than two element sizes, and the centroid tagging
// then gave triangles across an interface one material). Checked on the main and the
// Modes-tab mesh of a slotline on a thin finite substrate, on solder-masked lines,
// whose mask slivers must not trip the mesh-quality warning (Q > 100), and on layers of
// one er that differ only in conductivity or loss tangent.
import { CustomGeometrySolver } from '../src/custom_geometry.js';
import { MicrostripSolver } from '../src/microstrip.js';
import { check, quiet, APP, done } from './helpers.js';
import { buildTriRegions } from '../src/tri_solver/resample.js';

const SOLVE = { ...APP, max_iters: 12 };

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
    await quiet(() => s.solve_adaptive(SOLVE));
    const tb = s._triBackend, n = crossing(tb.mesh, s.dielectrics);
    check(`${name}: mesh follows the mask faces, no quality warning`, n === 0 && tb.meshQuality.maxQ < 100,
        `${n} of ${tb.mesh.nTris} triangles cross a face, max Q ${tb.meshQuality.maxQ.toFixed(1)}`);
}

// Layers of one er are still separate materials: the face between two that differ in
// sigma or tand only is meshed, and every triangle carries its own layer's loss.
for (const [name, lower, upper] of [['conductivity', 'sigma=10', ''], ['loss tangent', 'tand=0.05', 'tand=0.001']]) {
    const text = 'units um\nbounds open open open gnd\n' +
        `diel x=-inf w=inf y=0 h=150 er=11.9 ${lower}\ndiel x=-inf w=inf y=150 h=50 er=11.9 ${upper}\n` +
        'sig+ x=-10 w=20 y=200 h=2\n';
    const s = new CustomGeometrySolver({ text, sigma_cond: 5.8e7, freq: 1e9, nx: 10, ny: 10, mesh_backend: 'triangular' });
    // Nominal materials: the causal model rescales the loss map from er tand.
    s.use_causal_materials = false;
    await quiet(() => s.solve_adaptive(SOLVE));
    const mesh = s._triBackend.mesh, { nodes, tris, nTris, lossMap } = mesh;
    const want = yc => (yc < 150e-6 ? s.dielectrics[0] : s.dielectrics[1]);
    let wrong = 0, onFace = 0;
    for (let t = 0; t < nTris; t++) {
        const ys = [0, 1, 2].map(k => nodes[2 * tris[3 * t + k] + 1]);
        if (ys.filter(y => Math.abs(y - 150e-6) < 1e-12).length === 2) onFace++;
        const yc = (ys[0] + ys[1] + ys[2]) / 3;
        if (!(yc > 0 && yc < 200e-6) || mesh.epsMap[t].re < 11) continue;
        const d = want(yc), l = lossMap[t];
        if ((l.sigma || 0) !== (d.sigma || 0) || Math.abs(l.re - d.epsilon_r * (d.tan_delta || 0)) > 1e-12) wrong++;
    }
    check(`layers of one er differing in ${name}: the face between them is meshed`,
        crossing(mesh, s.dielectrics) === 0 && onFace > 4, `${crossing(mesh, s.dielectrics)} crossing, ${onFace} triangles on the face`);
    check(`layers of one er differing in ${name}: every triangle has its own layer's loss`, wrong === 0, `${wrong} wrong`);
    // The plotted field is recovered per material region, which must not span the face.
    const { regionOf } = buildTriRegions(mesh);
    // Region of the triangle nearest the face on the side `below` or above it.
    const regionNext = below => {
        let best = -1, dist = Infinity;
        for (let t = 0; t < nTris; t++) {
            const yc = [0, 1, 2].reduce((a, k) => a + nodes[2 * tris[3 * t + k] + 1], 0) / 3;
            const d = below ? 150e-6 - yc : yc - 150e-6;
            if (d > 0 && d < dist && yc > 0 && yc < 200e-6) { dist = d; best = regionOf[t]; }
        }
        return best;
    };
    const lo = regionNext(true), hi = regionNext(false);
    check(`layers of one er differing in ${name}: the plot regions split at the face`,
        lo >= 0 && hi >= 0 && lo !== hi, `regions ${lo} / ${hi}`);
}

done();

console.log('\nALL CHECKS PASSED');
