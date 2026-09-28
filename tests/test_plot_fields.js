// Plot fields at a chosen frequency (plotFieldsAt) and the surface current plot data.
//
//   - coax |E| from the resampler stays at the analytic peak 1/(a ln(b/a)): grid points
//     inside the unmeshed inner conductor must not enter the difference stencil
//   - |K| per ampere: the driven trace carries 1 A on both backends, half and full
//     domain, single and differential (per trace), coax K = 1/(2 pi r) on both surfaces
//   - full-wave backend: the E-field above the static limit is the eigenmode field,
//     with fieldFreq / fieldKind reported
//   - MQS |J|: 1 A through the trace at low frequency, skin-confined at high frequency
//   - MQS |K|: 1 A on the trace, exactly 1 A returned by the ground wall; the ideal-return
//     ground model is taken below the blend midpoint
//
// Run: node tests/test_plot_fields.js
import { CoaxSolver } from '../src/coax.js';
import { MicrostripSolver } from '../src/microstrip.js';
import { GroundedCPWSolver2D } from '../src/gcpw.js';
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
        if (r < (a + b) / 2) kin = Math.max(kin, K.K[i]); else kout = Math.min(kout, K.K[i]);
    }
    check('coax inner K = 1/(2 pi a)', Math.abs(kin * 2 * Math.PI * a - 1) < 0.03, `${kin.toFixed(1)} A/m`);
    check('coax outer K = 1/(2 pi b)', Math.abs(kout * 2 * Math.PI * b - 1) < 0.03, `${kout.toFixed(1)} A/m`);
}

// ---------------------------------------------------------------- microstrip
for (const backend of ['rectilinear', 'triangular']) {
    const s = new MicrostripSolver(MS);
    s.mesh_backend = backend;
    s.use_causal_materials = false;
    const pf = await solved(s, 20e9);
    const I = current(pf.surfaceK[0], (x, y) => y > 0.1e-3);
    // Full-wave |K| comes from the MQS current: its magnitude summed along the surface
    // exceeds the 1 A net current by the phase spread across the trace.
    check(`${backend} microstrip: trace carries 1 A`, Math.abs(I - 1) < (backend === 'triangular' ? 0.03 : 0.01), I.toFixed(4));
    const Ig = current(pf.surfaceK[0], (x, y) => y < 0.1e-3);
    check(`${backend} microstrip: ground return near 1 A`, Math.abs(Ig - 1) < 0.06, Ig.toFixed(4));
    check(`${backend} microstrip: field frequency reported`, pf.fieldFreq === 20e9, String(pf.fieldFreq));
    check(`${backend} microstrip: field kind`, pf.fieldKind === (backend === 'triangular' ? 'fullwave' : 'static'), pf.fieldKind);
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

    // |K| from the same solve: integral of J to the midline on the trace, the surface H
    // on the ground wall, which returns exactly the 1 A (Ampere) at any frequency.
    for (const f of [1e5, 10e9]) {
        await quiet(() => s.plotFieldsAt(f, cached));
        const pf = s.getPlotFields();
        const K = pf.surfaceK[0];
        const It = current(K, (x, y) => y > 0.1e-3), Ig = current(K, (x, y) => y < 0.1e-3);
        check(`MQS |K| at ${f / 1e9} GHz: from the MQS solve`, pf.surfaceKSource[0] === 'mqs', pf.surfaceKSource[0]);
        check(`MQS |K| at ${f / 1e9} GHz: trace 1 A`, Math.abs(It - 1) < 0.03, It.toFixed(4));
        check(`MQS |K| at ${f / 1e9} GHz: ground returns 1 A`, Math.abs(Ig - 1) < 0.01, Ig.toFixed(4));
    }
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
                Math.abs(Ip - 1) < 0.03 && Math.abs(In - 1) < 0.03, `${Ip.toFixed(3)} / ${In.toFixed(3)}`);
        });
    }
}

done();
