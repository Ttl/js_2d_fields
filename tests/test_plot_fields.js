// Plot fields at a chosen frequency (plotFieldsAt) and the surface current plot data.
//
//   - coax |E| from the resampler stays at the analytic peak 1/(a ln(b/a)): grid points
//     inside the unmeshed inner conductor must not enter the difference stencil
//   - |K| per ampere: the driven trace carries 1 A on both backends, half and full
//     domain, single and differential (per trace), coax K = 1/(2 pi r) on both surfaces
//   - full-wave backend: the E-field above the static limit is the static field plus the
//     eigenmode's change from F_STATIC_MAX (continuous there, as smooth as the static
//     field, sharp at dielectric interfaces), with fieldFreq / fieldKind reported
//   - MQS |J|: 1 A through the trace at low frequency, skin-confined at high frequency
//   - MQS |K| (tangential H): 1 A on the trace, exactly 1 A returned by the ground wall,
//     no collapse at the trace corners when the current fills the metal; the ideal-return
//     ground model is taken below the blend midpoint
//   - MQS |J| with meshed plating: one block per conductor, sharp at the plating/bulk interface
//   - MQS |J| of shaped conductors (elliptical coax): each conductor's own triangles, the
//     shield's block never picks up the centre conductor inside its bounding box
//   - the solve at the plot frequency keeps its fields: plotting there solves nothing and
//     matches a fresh solve, on both backends
//
// Run: node tests/test_plot_fields.js
import { CoaxSolver } from '../src/coax.js';
import { MicrostripSolver } from '../src/microstrip.js';
import { GroundedCPWSolver2D } from '../src/gcpw.js';
import { CustomGeometrySolver } from '../src/custom_geometry.js';
import { check, done, quiet } from './helpers.js';

const MS = { trace_width: 0.35e-3, substrate_height: 0.21e-3, trace_thickness: 35e-6, epsilon_r: 4.4,
    tan_delta: 0.02, sigma_cond: 5.8e7, freq: 1e9, boundaries: ['open', 'open', 'open', 'gnd'] };

// Current carried by the segments whose midpoint passes `sel`.
function current(K, sel) {
    let I = 0;
    for (let i = 0; i < K.K.length; i++) {
        const mx = (K.x0[i] + K.x1[i]) / 2, my = (K.y0[i] + K.y1[i]) / 2;
        if (sel(mx, my)) I += K.K[i] * Math.hypot(K.x1[i] - K.x0[i], K.y1[i] - K.y0[i]);
    }
    return I;
}

async function solved(s, f) {
    return quiet(async () => {
        const cached = await s.solve_adaptive({ max_nodes: 20000 });
        await s.plotFieldsAt(f, cached);
        return s.getPlotFields();
    });
}

// ---------------------------------------------------------------- coax
{
    const a = 0.46e-3, b = 1.475e-3;
    const s = new CoaxSolver({ inner_diameter: 2 * a, dielectric_diameter: 2 * b, epsilon_r: 2.1, freq: 1e9 });
    s.use_causal_materials = false;
    const pf = await solved(s, 1e9);
    let peak = 0;
    for (let j = 0; j < pf.y.length; j++) for (let i = 0; i < pf.x.length; i++)
        peak = Math.max(peak, Math.hypot(pf.Ex[0][j][i], pf.Ey[0][j][i]));
    const Emax = 1 / (a * Math.log(b / a));
    check('coax |E| peak at the analytic inner-surface field', Math.abs(peak / Emax - 1) < 0.05,
        `${peak.toFixed(0)} vs ${Emax.toFixed(0)} V/m`);
    const K = pf.surfaceK[0];
    let kin = 0, kout = Infinity;
    for (let i = 0; i < K.K.length; i++) {
        const r = Math.hypot((K.x0[i] + K.x1[i]) / 2, (K.y0[i] + K.y1[i]) / 2);
        // The shield's outer surface faces field-free air.
        if (r < (a + b) / 2) kin = Math.max(kin, K.K[i]); else if (r < 1.01 * b) kout = Math.min(kout, K.K[i]);
    }
    check('coax inner K = 1/(2 pi a)', Math.abs(kin * 2 * Math.PI * a - 1) < 0.03, `${kin.toFixed(1)} A/m`);
    check('coax outer K = 1/(2 pi b)', Math.abs(kout * 2 * Math.PI * b - 1) < 0.03, `${kout.toFixed(1)} A/m`);
    check('coax |K| from the MQS solve, skin mesh for the |J| overlay', pf.surfaceKSource[0] === 'mqs'
        && !!(pf.currentMesh && pf.currentMesh[0] && pf.currentMesh[0].nTris > 0), pf.surfaceKSource[0]);
}

// ---------------------------------------------------------------- microstrip
for (const backend of ['rectilinear', 'triangular']) {
    const s = new MicrostripSolver(MS);
    s.mesh_backend = backend;
    s.use_causal_materials = false;
    const pf = await solved(s, 20e9);
    const I = current(pf.surfaceK[0], (x, y) => y > 0.1e-3);
    check(`${backend} microstrip: trace carries 1 A`, Math.abs(I - 1) < 0.01, I.toFixed(4));
    const Ig = current(pf.surfaceK[0], (x, y) => y < 0.1e-3);
    check(`${backend} microstrip: ground return near 1 A`, Math.abs(Ig - 1) < 0.06, Ig.toFixed(4));
    check(`${backend} microstrip: field frequency reported`, pf.fieldFreq === 20e9, String(pf.fieldFreq));
    check(`${backend} microstrip: field kind`, pf.fieldKind === (backend === 'triangular' ? 'fullwave' : 'static'), pf.fieldKind);
}

// ---------------------------------------------------------------- full-wave E quality
// The full-wave E plot is the static field plus the change of the mode field from
// F_STATIC_MAX: continuous with the static field there, and as smooth as it at high
// frequency (the mode field sampled per element kinks the contours at element edges).
{
    // Mean relative deviation of |E| from the linear interpolation of its neighbours
    // along x and y, over samples above 1 % of the peak: kinks show as local bumps.
    const roughness = (pf) => {
        const { x, y } = pf, Ex = pf.Ex[0], Ey = pf.Ey[0];
        const E = (j, i) => Math.hypot(Ex[j][i], Ey[j][i]);
        let peak = 0, sum = 0, n = 0;
        for (let j = 0; j < y.length; j++) for (let i = 0; i < x.length; i++) peak = Math.max(peak, E(j, i));
        for (let j = 2; j < y.length - 2; j++) for (let i = 2; i < x.length - 2; i++) {
            const c = E(j, i), l = E(j, i - 1), r = E(j, i + 1), d = E(j - 1, i), u = E(j + 1, i);
            if (![c, l, r, d, u].every(v => v > 0.01 * peak)) continue;
            const ax = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]), ay = (y[j] - y[j - 1]) / (y[j + 1] - y[j - 1]);
            sum += Math.abs(c - (l * (1 - ax) + r * ax)) / c + Math.abs(c - (d * (1 - ay) + u * ay)) / c;
            n += 2;
        }
        return sum / n;
    };
    const relDiff = (a, b) => {
        let num = 0, den = 0;
        for (let j = 0; j < a.y.length; j++) for (let i = 0; i < a.x.length; i++) {
            num += (a.Ex[0][j][i] - b.Ex[0][j][i]) ** 2 + (a.Ey[0][j][i] - b.Ey[0][j][i]) ** 2;
            den += b.Ex[0][j][i] ** 2 + b.Ey[0][j][i] ** 2;
        }
        return Math.sqrt(num / den);
    };
    const s = new MicrostripSolver({ ...MS, mesh_backend: 'triangular' });
    s.use_causal_materials = false;
    const cached = await quiet(() => s.solve_adaptive({ max_nodes: 20000 }));
    await quiet(() => s.plotFieldsAt(50e6, cached));
    const st = s.getPlotFields();
    await quiet(() => s.plotFieldsAt(150e6, cached));
    const lo = s.getPlotFields();
    const dLo = relDiff(lo, st);
    check('full-wave E just above 100 MHz continues the static field', lo.fieldKind === 'fullwave' && dLo < 0.005,
        `${lo.fieldKind}, ${(dLo * 100).toFixed(3)} %`);
    // At exactly 100 MHz the mode is its own reference: as smooth as the static field.
    await quiet(() => s.plotFieldsAt(100e6, cached));
    const at = s.getPlotFields();
    const rAt = roughness(at) / roughness(st);
    check('full-wave E at exactly 100 MHz as smooth as the static field', at.fieldKind === 'fullwave' && rAt < 1.05,
        `${at.fieldKind}, x${rAt.toFixed(2)}`);
    await quiet(() => s.plotFieldsAt(20e9, cached));
    const hi = s.getPlotFields();
    const r = roughness(hi) / roughness(st);
    check('full-wave E at 20 GHz as smooth as the static field', hi.fieldKind === 'fullwave' && r < 1.15, `x${r.toFixed(2)}`);
    // Normal D is continuous across the substrate top at any frequency: E_y just below
    // and just above it (clear of the trace) differ by the permittivity ratio, in one
    // grid row, the differencing baseline must not blend the two dielectrics.
    for (const [label, pf] of [['static', st], ['full-wave 20 GHz', hi]]) {
        const i = pf.x.findIndex(v => v > 0.45e-3);
        const h = MS.substrate_height;
        const jb = pf.y.findIndex(v => v >= h) - 1, ja = pf.y.findIndex(v => v > h);
        const ratio = pf.Ey[0][ja][i] / pf.Ey[0][jb][i];
        check(`${label} E_y jumps by er across the substrate top`, Math.abs(ratio / MS.epsilon_r - 1) < 0.05,
            `${ratio.toFixed(3)} at y ${(pf.y[jb] * 1e3).toFixed(4)} / ${(pf.y[ja] * 1e3).toFixed(4)} mm`);
    }
    // |E| on the mesh for the image and traced contours: at the substrate top, clear of
    // the trace, a vertex read from the air triangles and from the substrate triangles
    // jumps as the grid rows on either side of the interface do (normal D continuous,
    // tangential E continuous), not smeared across it.
    {
        const M = hi.fieldMesh && hi.fieldMesh[0];
        check('full-wave E on the mesh present', !!(M && M.E.length && M.E.every(isFinite)), M ? `${M.E.length / 3} triangles` : 'none');
        // A point shared by triangles of one material reads one value in all of them, or a
        // contour line ends on the shared edge (the trace corners, where the field is
        // singular, used to differ by tens of percent).
        {
            const seen = new Map();
            let worst = 0, at = '';
            for (let t = 0; M && t < M.E.length / 3; t++) {
                const cy = (M.tris[6 * t + 1] + M.tris[6 * t + 3] + M.tris[6 * t + 5]) / 3;
                for (let a = 0; a < 3; a++) {
                    const x = M.tris[6 * t + 2 * a], y = M.tris[6 * t + 2 * a + 1];
                    const key = `${Math.round(x * 1e10)},${Math.round(y * 1e10)},${cy > MS.substrate_height}`;
                    const v = M.E[3 * t + a], u = seen.get(key);
                    if (u === undefined) { seen.set(key, v); continue; }
                    const d = Math.abs(v - u) / Math.max(u, v);
                    if (d > worst) { worst = d; at = `(${(x * 1e3).toFixed(4)}, ${(y * 1e3).toFixed(4)}) mm`; }
                }
            }
            check('full-wave E on the mesh: one value per shared point within a material', worst < 1e-9,
                `worst ${(worst * 100).toFixed(2)} % at ${at}`);
        }
        const h = MS.substrate_height, above = [], below = [];
        for (let t = 0; M && t < M.E.length / 3; t++) {
            const cy = (M.tris[6 * t + 1] + M.tris[6 * t + 3] + M.tris[6 * t + 5]) / 3;
            for (let a = 0; a < 3; a++) {
                const x = M.tris[6 * t + 2 * a], y = M.tris[6 * t + 2 * a + 1];
                if (Math.abs(y - h) > 1e-9 || x < 0.4e-3 || x > 0.6e-3) continue;
                (cy > h ? above : below).push([x, M.E[3 * t + a]]);
            }
        }
        const jA = hi.y.findIndex(v => v >= h), jB = jA - 1;
        const Eg = (j, i) => Math.hypot(hi.Ex[0][j][i], hi.Ey[0][j][i]);
        const errs = [];
        for (const [x, e] of above) {
            const b = below.find(([xb]) => Math.abs(xb - x) < 1e-9);
            if (!b) continue;
            const i = hi.x.findIndex(v => v >= x);
            errs.push(Math.abs((e / b[1]) / (Eg(jA, i) / Eg(jB, i)) - 1));
        }
        errs.sort((a, b) => a - b);
        const med = errs[errs.length >> 1];
        check('full-wave E on the mesh jumps at the substrate top like the grid', errs.length > 0 && med < 0.05,
            `${errs.length} vertices, median deviation ${med !== undefined ? (med * 100).toFixed(1) : '-'} %`);
    }
    // At 100 GHz the mode differs strongly from the static field: the grid row on the
    // substrate top must read the same (air) side for both parts, and the outermost
    // columns, where the static field has no sample, the mode itself.
    await quiet(() => s.plotFieldsAt(100e9, cached));
    const vh = s.getPlotFields();
    const E = (j, i) => Math.hypot(vh.Ex[0][j][i], vh.Ey[0][j][i]);
    const iC = vh.x.findIndex(v => v > 0.45e-3);
    const jI = vh.y.findIndex(v => v >= MS.substrate_height);
    const rowRatio = E(jI, iC) / E(jI + 1, iC);
    check('full-wave E at 100 GHz: the interface row continues the air side', Math.abs(rowRatio - 1) < 0.05, rowRatio.toFixed(3));
    const nx = vh.x.length, jM = vh.y.findIndex(v => v > 0.4e-3);
    const edge = E(jM, nx - 1), inner = E(jM, nx - 2);
    check('full-wave E at 100 GHz: no spike on the open side wall', edge <= 2 * inner + 1e-9 * E(jI, iC),
        `${edge.toFixed(2)} vs ${inner.toFixed(2)} V/m`);
    // The mode does change: the field pulls into the substrate at 20 GHz.
    const dHi = relDiff(hi, st);
    check('full-wave E at 20 GHz differs from the static field', dHi > 0.02, `${(dHi * 100).toFixed(2)} %`);
}

// ---------------------------------------------------------------- MQS current density
{
    const s = new MicrostripSolver({ ...MS, mesh_backend: 'triangular' });
    s.use_causal_materials = false;
    const cached = await quiet(() => s.solve_adaptive({ max_nodes: 20000 }));
    // Trapezoid integral of |J| over the sampled blocks (both halves of the trace).
    const integral = (blocks) => {
        let I = 0;
        for (const b of blocks) {
            for (let j = 0; j < b.y.length - 1; j++) for (let i = 0; i < b.x.length - 1; i++) {
                const c = [b.J[j][i], b.J[j][i + 1], b.J[j + 1][i], b.J[j + 1][i + 1]];
                if (c.some(v => v == null)) continue;
                I += (c[0] + c[1] + c[2] + c[3]) / 4 * Math.abs(b.x[i + 1] - b.x[i]) * (b.y[j + 1] - b.y[j]);
            }
        }
        return I;
    };
    await quiet(() => s.plotFieldsAt(1e6, cached));
    const J = s.getPlotFields().currentJ;
    check('MQS |J| present on the full-wave backend', !!(J && J[0] && J[0].length === 2), String(J && J[0] && J[0].length));
    // At 1 MHz the skin depth is twice the thickness: the current is nearly in phase,
    // so the integral of |J| is the 1 A trace current.
    const I1 = integral(J[0]);
    check('MQS |J| at 1 MHz integrates to 1 A', Math.abs(I1 - 1) < 0.03, I1.toFixed(4));
    await quiet(() => s.plotFieldsAt(10e9, cached));
    const blocks = s.getPlotFields().currentJ[0];
    // At 10 GHz the current sits in the skin: the mid-thickness row of the trace centre
    // carries far less than its faces.
    const b = blocks[0], iMid = Math.floor(b.x.length / 2);
    const col = b.J.map(row => row[iMid]).filter(v => v != null);
    const mid = col[Math.floor(col.length / 2)], face = Math.max(col[0], col[col.length - 1]);
    check('MQS |J| at 10 GHz confined to the skin', mid < 1e-3 * face, `${mid.toExponential(2)} vs ${face.toExponential(2)}`);

    // |K| from the same solve: the tangential H on every metal surface, whose contour
    // integral is the enclosed current (Ampere) at any frequency.
    for (const f of [1e5, 10e9]) {
        await quiet(() => s.plotFieldsAt(f, cached));
        const pf = s.getPlotFields();
        const K = pf.surfaceK[0];
        const It = current(K, (x, y) => y > 0.1e-3), Ig = current(K, (x, y) => y < 0.1e-3);
        check(`MQS |K| at ${f / 1e9} GHz: from the MQS solve`, pf.surfaceKSource[0] === 'mqs', pf.surfaceKSource[0]);
        check(`MQS |K| at ${f / 1e9} GHz: trace 1 A`, Math.abs(It - 1) < 0.01, It.toFixed(4));
        check(`MQS |K| at ${f / 1e9} GHz: ground returns 1 A`, Math.abs(Ig - 1) < 0.01, Ig.toFixed(4));
        check(`MQS |K| at ${f / 1e9} GHz: skin mesh for the overlay`, !!(pf.currentMesh && pf.currentMesh[0] && pf.currentMesh[0].nTris > 0));
        if (f === 1e5) {
            // Near DC the current fills the trace: the surface field falls towards the
            // corners but stays of the order of the centre value.
            const top = MS.substrate_height + MS.trace_thickness, edge = MS.trace_width / 2;
            const near = (x0) => {
                let best = null, d = Infinity;
                for (let i = 0; i < K.K.length; i++) {
                    const mx = (K.x0[i] + K.x1[i]) / 2, my = (K.y0[i] + K.y1[i]) / 2;
                    if (Math.abs(my - top) > 1e-7) continue;
                    if (Math.abs(mx - x0) < d) { d = Math.abs(mx - x0); best = K.K[i]; }
                }
                return best;
            };
            const ratio = near(edge - 1e-6) / near(0);
            check('MQS |K| near DC: no collapse at the trace corner', ratio > 0.3, ratio.toFixed(3));
        }
    }
}

// ---------------------------------------------------------------- shaped conductor |J|
{
    const text = `units mm
d = 0.9; a = 2.2; b = 1.4; t_sh = 0.15
bounds open open open open
domain -1.1*(a+t_sh) 1.1*(a+t_sh) -1.1*(b+t_sh) 1.1*(b+t_sh)
diel  ellipse  x=0  y=0  rx=a  ry=b  n=128  er=2.1  tand=0.0002
sig+  ngon  x=0  y=0  r=d/2  n=64
gnd   ellipse  x=0  y=0  rx=a+t_sh  ry=b+t_sh  rx_in=a  ry_in=b  n=128
`;
    const s = new CustomGeometrySolver({ text, freq: 1e9, mesh_backend: 'triangular' });
    s.use_causal_materials = false;
    const cached = await quiet(() => s.solve_adaptive({ max_nodes: 20000 }));
    await quiet(() => s.plotFieldsAt(1e9, cached));
    const blocks = s.getPlotFields().currentJ[0];
    const radius = b => {
        let lo = Infinity, hi = 0;
        for (let i = 0; i < b.tris.length; i += 2) { const r = Math.hypot(b.tris[i], b.tris[i + 1]); lo = Math.min(lo, r); hi = Math.max(hi, r); }
        return [lo, hi];
    };
    const r = blocks.map(radius);
    const centre = r.filter(([, hi]) => hi < 0.46e-3).length, shield = r.filter(([lo]) => lo > 1.3e-3).length;
    check('shaped |J|: triangle blocks, one per conductor and mirror half', blocks.length === 4 && blocks.every(b => b.tris && b.Jv),
        `${blocks.length} blocks`);
    check('shaped |J|: shield block holds only shield metal', centre === 2 && shield === 2,
        r.map(([lo, hi]) => `${(lo * 1e3).toFixed(3)}..${(hi * 1e3).toFixed(3)}`).join(' '));
}

// ---------------------------------------------------------------- MQS |J| with meshed plating
// The plating layer and its core are one block: no unsampled hole where the core is,
// and the doubled interface lines carry the jump of sigma.
{
    const pl = { sigma: 8.7e6, thickness: 3e-6, rq: 0, top: true, sides: true, bottom: false, thick_corners: true };
    const s = new MicrostripSolver({ ...MS, plating: pl, mesh_backend: 'triangular' });
    s.use_causal_materials = false;
    const cached = await quiet(() => s.solve_adaptive({ max_nodes: 20000 }));
    await quiet(() => s.plotFieldsAt(1e9, cached));
    const blocks = s.getPlotFields().currentJ[0];
    const nulls = blocks.reduce((a, b) => a + b.J.flat().filter(v => v == null).length, 0);
    check('plated |J|: one block per trace half', blocks.length === 2, `${blocks.length} blocks`);
    check('plated |J|: the block covers the core', nulls === 0, `${nulls} unsampled`);
    const b = blocks[0], j = Math.floor(b.y.length / 2);
    const i = b.x.findIndex((v, k) => v === b.x[k + 1]);
    const ratio = i >= 0 ? b.J[j][i] / b.J[j][i + 1] : NaN;
    check('plated |J|: sigma jump at the side interface', Math.abs(ratio - 5.8e7 / 8.7e6) < 0.05 * 5.8e7 / 8.7e6,
        ratio.toFixed(3));
}

// ---------------------------------------------------------------- plot frequency cache
// The solve at the plot frequency keeps its fields: plotting there after the solve runs
// no MQS or eigen solve (no static re-solve on the rectilinear backend), and gives the
// fields a fresh solve at that frequency gives.
{
    const f = 5e9;
    const s = new MicrostripSolver({ ...MS, freq: f, mesh_backend: 'triangular' });
    s.plot_freq_target = f;
    const cached = await quiet(() => s.solve_adaptive({ max_nodes: 20000 }));
    await quiet(() => s.computeAtFrequency(1e9, cached));   // a sweep point after it
    const tri = await s._ensureTriBackend();
    let calls = 0;
    for (const k of ['_mqsSolve', '_eigenPick', '_modeAtFreq']) {
        const o = tri[k].bind(tri);
        tri[k] = (...a) => { calls++; return o(...a); };
    }
    await quiet(() => s.plotFieldsAt(f, cached));
    check('tri plot at the solved plot frequency: no solve', calls === 0, `${calls} calls`);
    const a = s.getPlotFields();
    await quiet(() => s.plotFieldsAt(2e9, cached));
    await quiet(() => s.plotFieldsAt(f, cached));
    const b = s.getPlotFields();
    const maxRel = (p, q) => {
        let m = 0, r = 0;
        for (let i = 0; i < p.length; i++) { m = Math.max(m, Math.abs(p[i] - q[i])); r = Math.max(r, Math.abs(q[i])); }
        return m / r;
    };
    const dK = maxRel(a.surfaceK[0].K, b.surfaceK[0].K);
    const dE = maxRel(a.Ex[0].flatMap(r => Array.from(r)), b.Ex[0].flatMap(r => Array.from(r)));
    check('tri cached plot fields match a fresh solve', a.surfaceK[0].K.length === b.surfaceK[0].K.length && dK < 1e-9 && dE < 1e-6,
        `K ${dK.toExponential(1)}, Ex ${dE.toExponential(1)}`);
    check('tri cached plot: full-wave field at the plot frequency', a.fieldKind === 'fullwave' && a.fieldFreq === f, `${a.fieldKind} ${a.fieldFreq}`);
}
{
    const f = 5e9;
    const s = new MicrostripSolver({ ...MS, freq: f, mesh_backend: 'rectilinear' });
    s.plot_freq_target = f;
    const cached = await quiet(() => s.solve_adaptive({ max_nodes: 20000 }));
    await quiet(() => s.computeAtFrequency(f, cached));
    await quiet(() => s.computeAtFrequency(1e9, cached));
    let calls = 0;
    const o = s.computeAtFrequency.bind(s);
    s.computeAtFrequency = (...a) => { calls++; return o(...a); };
    await quiet(() => s.plotFieldsAt(f, cached));
    check('rectilinear causal plot at the solved plot frequency: no solve', calls === 0, `${calls} calls`);
}

// ---------------------------------------------------------------- ideal-ground blend
// Coplanar grounds reaching the domain edge are ideal returns towards DC: the plot takes
// the field of the ground model with the larger blend weight, and says so.
{
    const s = new GroundedCPWSolver2D({ trace_width: 0.35e-3, gap: 0.1e-3, via_gap: 0.2e-3, gnd_thickness: 35e-6,
        substrate_height: 0.21e-3, trace_thickness: 35e-6, epsilon_r: 4.4, tan_delta: 0.02, sigma_cond: 5.8e7, freq: 1e9 });
    s.mesh_backend = 'triangular';
    s.use_causal_materials = false;
    const cached = await quiet(() => s.solve_adaptive({ max_nodes: 20000 }));
    for (const [f, ideal] of [[1e5, true], [1e8, false]]) {
        await quiet(() => s.plotFieldsAt(f, cached));
        const pf = s.getPlotFields();
        check(`GCPW at ${f / 1e6} MHz: ${ideal ? 'ideal' : 'finite'} ground field`, pf.idealGrounds && pf.idealGrounds[0] === ideal,
            JSON.stringify(pf.idealGrounds));
        const K = pf.surfaceK[0];
        const It = current(K, (x, y) => Math.abs(x) < 0.176e-3 && y > 0.2e-3);
        const Io = current(K, (x, y) => !(Math.abs(x) < 0.176e-3 && y > 0.2e-3));
        check(`GCPW at ${f / 1e6} MHz: return balances the trace`, Math.abs(Io / It - 1) < 0.03, `${It.toFixed(3)} / ${Io.toFixed(3)}`);
    }
}

// ---------------------------------------------------------------- differential pair
for (const backend of ['rectilinear', 'triangular']) {
    for (const full of [false, true]) {
        const s = new MicrostripSolver({ ...MS, trace_width: 0.2e-3, trace_spacing: 0.15e-3,
            mesh_backend: backend, ...(full ? { symmetry: false } : {}) });
        s.use_causal_materials = false;
        check(`${backend}${full ? ' full' : ''} diff: domain`, !!s.sym_half === (backend === 'rectilinear' && !full));
        const pf = await solved(s, 5e9);
        ['odd', 'even'].forEach((mode, m) => {
            const K = pf.surfaceK[m];
            const Ip = current(K, (x, y) => y > 0.1e-3 && x > 0), In = current(K, (x, y) => y > 0.1e-3 && x < 0);
            check(`${backend}${full ? ' full' : ''} diff ${mode}: 1 A per trace`,
                Math.abs(Ip - 1) < 0.01 && Math.abs(In - 1) < 0.01, `${Ip.toFixed(3)} / ${In.toFixed(3)}`);
        });
    }
}

done();
