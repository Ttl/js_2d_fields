// Full-wave mesh budget on a wide geometry with long grounds (fuzzer seed-2 custom case
// 91: a pair on five full-width layers, a coplanar ground running to the domain edge 31 um
// from one trace, a finite bottom ground, 13 mm domain):
//   1. ground surfaces far from the signals are painted coarser with the distance: every
//      ground curve used to take the signal surface size along its whole length, and 43%
//      of the triangles sat along grounds more than 1 mm from any signal.
//   2. the initial mesh leaves budget for the adaptive passes: it was accepted at 1.3x the
//      budget, no pass ran, the 31 um gap had 49 triangles and the certificate failed.
//      The pair is asymmetric, so the certificate covers its modal C and C0.
//   3. the triangle-count estimate of the graded ground field stays within 35% of the
//      mesh gmsh builds (the pre-mesh sizing relies on it; the layer interfaces it does
//      not count make it read low on this many-layer stack).
//   4. an asymmetric pair (unequal traces, a ground on one side) reports its line modes,
//      C = 1/2 v^T C v with the modal vector v of the mesh: the certificate measures those
//      quantities, vector included. It measured the +-1 drives, and certified a mesh whose
//      odd C was 1.3% off at 0.6% estimated error.
//
// Run: node tests/test_fullwave_mesh_budget.js
import { CustomGeometrySolver } from '../src/custom_geometry.js';
import { estimateOccTriCount } from '../src/tri_solver/occ_to_mesh.js';
import { check, quiet, rel, APP, done } from './helpers.js';

const text = [
    'units um', 'bounds gnd gnd open open', 'domain -6676.1 6676.1 auto auto',
    'diel x=-inf w=inf y=0 h=304.04 er=6.7129 tand=0.0029469',
    'diel x=-inf w=inf y=304.04 h=360.2 er=9.162 tand=0.01252',
    'diel x=-inf w=inf y=664.25 h=599.7 er=6.7373 tand=0.0034598',
    'diel x=-inf w=inf y=1263.9 h=218.13 er=6.5341 tand=0.024448',
    'diel x=-inf w=inf y=1482.1 h=90.575 er=5.7959 tand=0.0062508',
    'diel x=-2768 w=3334 y=1077.4 h=76.463 er=6.8897 tand=0.011592',
    'sig- x=-812.18 w=696.81 y=1482.1 h=43.513 sigma=3746700 rq=0.14317',
    'sig+ x=115.37 w=696.81 y=1482.1 h=43.513',
    'gnd x=-1883.9 w=3767.9 y=-6.8247 h=6.8247 sigma=641000 rq=0.51154',
    'gnd x=843.13 w=inf y=1482.1 h=43.513 rq=0.38213',
].join('\n');
// Even-mode C (pF/m) on a 160k-node mesh at 0.05% tolerance.
const C_EVEN_REF = 206.0;

const s = new CustomGeometrySolver({ text, sigma_cond: 11e6, freq: 4.69e6, mesh_backend: 'triangular', nx: 30, ny: 30 });
const r = await quiet(() => s.solve_adaptive(APP));
const tb = s._triBackend, mesh = tb.mesh;

// Triangles near a ground surface and far from both traces.
const gnds = [[-1883.9, 1883.9, -6.8247, 0], [843.13, 7000, 1482.1, 1525.6]];
let far = 0, gap = 0;
for (let t = 0; t < mesh.nTris; t++) {
    let x = 0, y = 0;
    for (let k = 0; k < 3; k++) { const v = mesh.tris[3 * t + k]; x += mesh.nodes[2 * v] / 3e-6; y += mesh.nodes[2 * v + 1] / 3e-6; }
    const dSig = Math.hypot(Math.max(0, Math.abs(x) - 812.18), Math.max(0, Math.abs(y - 1503.9) - 21.8));
    const dGnd = Math.min(...gnds.map(g => Math.hypot(Math.max(0, g[0] - x, x - g[1]), Math.max(0, g[2] - y, y - g[3]))));
    if (dGnd < 150 && dSig > 1000) far++;
    if (Math.hypot(x - 827.6, y - 1503.9) < 100) gap++;
}
check('1. ground surfaces far from the signals hold few triangles', far < 0.15 * mesh.nTris,
    `${far} of ${mesh.nTris}`);
check('2. the initial mesh leaves budget for adaptive passes', tb.initialTris < 0.5 * tb.maxTris && mesh.nTris > tb.initialTris,
    `initial ${tb.initialTris}, final ${mesh.nTris}, budget ${tb.maxTris}`);
check('2. the gap between trace and coplanar ground is resolved', gap > 150, `${gap} triangles within 100 um`);
check('2. the static certificate passes', !!(tb.certification && tb.certification.pass),
    JSON.stringify(tb.certification && { pass: tb.certification.pass, err: tb.certification.err }));
const even = r.modes.find(m => m.mode === 'even');
check('2. even-mode C matches the fine-mesh reference', rel(even.RLGC.C * 1e12, C_EVEN_REF) < 0.01,
    `${(even.RLGC.C * 1e12).toFixed(2)} vs ${C_EVEN_REF} pF/m`);

// 3. The estimate for the sizing the initial mesh was built with, against its count.
const st = tb.initialSizeStats;
const est = estimateOccTriCount(st, st);
check('3. the graded ground estimate is within 35% of the built mesh', Math.abs(est / tb.initialTris - 1) < 0.35,
    `estimate ${est.toFixed(0)}, built ${tb.initialTris}`);

// 4. Certificate quantities of an asymmetric pair, at the 100 MHz anchor (no dispersion).
{
    const text = [
        'units um', 'bounds open open open gnd',
        'diel x=-inf w=inf y=0 h=364.21 er=4.4645', 'diel x=-inf w=inf y=364.21 h=276.12 er=7.4634',
        'diel x=-inf w=inf y=640.32 h=41.471 er=6.122', 'diel x=-369.68 w=939.02 y=340.87 h=48.29 er=4.8587',
        'sig+ x=-163.15 w=152.47 y=364.21 h=35.985', 'sig- x=10.682 w=300.99 y=364.21 h=35.985',
        'gnd x=-615.85 w=-inf y=364.21 h=35.985', 'gnd x=-406.21 w=960.95 y=681.79 h=6.5416',
    ].join('\n');
    const s4 = new CustomGeometrySolver({ text, freq: 1e8, mesh_backend: 'triangular', nx: 30, ny: 30 });
    const r4 = await quiet(() => s4.solve_adaptive({ ...APP, max_nodes: 10000 }));
    const t4 = s4._triBackend;
    const vecs = t4._certModal && t4._certModal.vecs;
    check('4. the pair is decomposed into line modes', !!vecs, JSON.stringify(vecs));
    if (vecs) {
        const q = t4._staticEnergies(t4.mesh, vecs), q1 = t4._staticEnergies(t4.mesh);
        const eps0 = 8.854187817e-12;
        r4.modes.forEach((m, k) => {
            check(`4. ${m.mode}: the certified energy is the reported C`, rel(m.RLGC.C, eps0 * q[2 * k]) < 1e-6,
                `${m.RLGC.C} vs ${eps0 * q[2 * k]}`);
        });
        check('4. the +-1 drives are not the reported odd mode', rel(r4.modes[0].RLGC.C, eps0 * q1[0]) > 0.01,
            `${(100 * rel(r4.modes[0].RLGC.C, eps0 * q1[0])).toFixed(1)}%`);
    }
}
done();
