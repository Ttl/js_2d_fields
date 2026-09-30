// Anisotropic MQS skin band (refineSkinBand aniso, the default): elements stretched
// along the conductor surface, sized across it by the skin depth.
//   - the band converges without degenerate triangles where the normal is undefined (a
//     microstrip at 100 MHz: the band reaches the medial axis of the trace, where a
//     per-edge metric used to bisect short edges down to needles and fail the solve)
//   - R matches the isotropic band, uncapped, with a fraction of its triangles
//   - the coax as n-gons at 40 GHz: the isotropic band runs into its triangle budget
//     there (0.33 um skin depth on 12 mm of perimeter), the anisotropic one resolves it
//     and R matches the closed form
//   - a corner rounded far below the skin depth is still a corner (the 10 nm radius
//     against the sharp trace is pinned in test_custom_geometry_shapes.js)
//
// Run: node tests/test_mqs_aniso_band.js
import { CustomGeometrySolver } from '../src/custom_geometry.js';
import { check, quiet, rel, done } from './helpers.js';

const MU0 = 4e-7 * Math.PI, SIG = 5.8e7;
const MS = `units mm
w = 0.35; t = 0.035; h = 0.21
bounds open open open gnd
domain auto
diel  x=-inf  y=0  w=inf  h=h  er=4.4  tand=0.02
sig+  x=-w/2  y=h  w=w  h=t`;

async function solveAt(text, f, triOpts = {}) {
    const s = new CustomGeometrySolver({ text, freq: f, mesh_backend: 'triangular', sigma_cond: SIG });
    s.use_causal_materials = false;
    s.tri_opts = { mqsInterpTol: 0, dispTol: 0, ...triOpts };
    const first = await quiet(() => s.solve_adaptive({ max_nodes: 20000 }));
    const r = await quiet(() => s.computeAtFrequency(f, first));
    return { R: r.modes[0].RLGC.R, mesh: s._triBackend._skinCache.mesh, warn: (r.warnings || []).map(w => w.type) };
}

// Smallest area / longest-edge^2 over the triangles: 0 for a degenerate one.
function minShape(m) {
    let q = Infinity;
    for (let t = 0; t < m.nTris; t++) {
        const p = [0, 1, 2].map(k => [m.nodes[2 * m.tris[3 * t + k]], m.nodes[2 * m.tris[3 * t + k] + 1]]);
        const A = Math.abs((p[1][0] - p[0][0]) * (p[2][1] - p[0][1]) - (p[2][0] - p[0][0]) * (p[1][1] - p[0][1])) / 2;
        const L = Math.max(...[0, 1, 2].map(k => Math.hypot(p[(k + 1) % 3][0] - p[k][0], p[(k + 1) % 3][1] - p[k][1])));
        q = Math.min(q, A / (L * L));
    }
    return q;
}
const failed = w => w.filter(t => t === 'mqs-band-capped' || t === 'mqs-solve-failed');

// --- band reaching the medial axis of the trace ---
{
    const a = await solveAt(MS, 1e8);
    check('microstrip 100 MHz: band converges, MQS solve runs', failed(a.warn).length === 0, JSON.stringify(a.warn));
    check('microstrip 100 MHz: no degenerate triangles', minShape(a.mesh) > 1e-6, minShape(a.mesh).toExponential(2));
}

// --- against the isotropic band ---
{
    const iso = await solveAt(MS, 1e10, { mqsBandAniso: false, mqsMaxTris: 400000 });
    const an = await solveAt(MS, 1e10);
    check('microstrip 10 GHz: R = isotropic band', rel(an.R, iso.R) < 1e-3, `${an.R.toFixed(3)} / ${iso.R.toFixed(3)} ohm/m`);
    check('microstrip 10 GHz: under half the triangles', an.mesh.nTris < 0.5 * iso.mesh.nTris,
        `${an.mesh.nTris} / ${iso.mesh.nTris}`);
}

// --- coax as n-gons in the deep skin regime ---
{
    const a = 0.46e-3, b = 1.475e-3, f = 4e10;
    const text = `units mm
r1 = 0.460179; r2 = 1.4753; r3 = 1.62283
bounds open open open open
domain -1.1*r3 1.1*r3 -1.1*r3 1.1*r3
diel  ngon  x=0  y=0  r=r2  n=128  er=2.1  tand=0
sig+  ngon  x=0  y=0  r=r1  n=92
gnd   ngon  x=0  y=0  r=r3  r_in=r2  n=128`;
    const iso = await solveAt(text, f, { mqsBandAniso: false });
    check('coax 40 GHz: the isotropic band is capped', iso.warn.includes('mqs-band-capped'), JSON.stringify(iso.warn));
    const an = await solveAt(text, f);
    check('coax 40 GHz: anisotropic band not capped', failed(an.warn).length === 0, JSON.stringify(an.warn));
    const Rs = Math.sqrt(Math.PI * f * MU0 / SIG);
    const Rx = Rs / (2 * Math.PI) * (1 / a + 1 / b);
    check('coax 40 GHz: R = closed form', rel(an.R, Rx) < 3e-3, `${an.R.toFixed(3)} vs ${Rx.toFixed(3)} ohm/m`);
    check('coax 40 GHz: no degenerate triangles', minShape(an.mesh) > 1e-6, minShape(an.mesh).toExponential(2));
}

done();
