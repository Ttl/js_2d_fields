// Kramers-Kronig pairing of the full-wave MQS R(f) and L(f) with meshed plating.
//
// A causal series impedance has L(f) - X_KK(f)/w constant, X_KK the Hilbert transform
// of R(f). The case is a diff microstrip with thick-corner plating (meshed layer) of
// another conductivity than the bulk and a rough bulk: the trace bottom the plating
// leaves bare is the only rough face. Its roughness once scaled the whole core volume
// loss, most of which sits under the plating and moves as delta passes the layer, and
// the Eye tab warned of a precursor (L - X_KK/w drifted 3 % over 0.1-20 GHz). The bare
// bottom strips of the layer under the side plating took no roughness either when the
// bulk roughness was the solver-wide one (0.46 % with plating of the bulk sigma). The
// rough faces of a meshed layer are weighted by the perfect-conductor |K|^2: the
// finite-sigma surface integral rises ~15 % over the band, and with a rough layer's
// reactance several times its resistance that drift showed in L(f) alone (1.2 %
// with plating on every face, rq 3 um under 2 um).
//
// Run: node tests/test_plating_kk.js
import { buildSolverFromParams } from '../src/solver_factory.js';
import { InterpolatingSweep } from '../src/interpolating_sweep.js';
import { channelTableFromRows, makeRLGCEvaluator } from '../src/eye/channel_grid.js';
import { check, quiet, done } from './helpers.js';

const BASE = {
    tl_type: 'diff_microstrip', mesh_backend: 'fullwave_mqs', custom_geom: '', custom_overrides: {},
    w: 0.35e-3, h: 0.21e-3, t: 35e-6, er: 4.4, tand: 0.002, sigma: 5.8e7, freq: 1e9, nx: 30, ny: 30,
    trace_spacing: 0.2e-3, gap: 0.1e-3, via_gap: 0.1e-3, gnd_width: NaN,
    max_iters: 10, tolerance: 0.01, min_converged_passes: 2, estimate_error: 1, max_nodes: 20,
    rq: 1e-6, use_plating: true, plating_sigma: 1e7, plating_t: 4e-6, plating_rq: 0, plating_rq_iface: NaN,
    plating_top: true, plating_sides: true, plating_bottom: false, plating_thick_corners: true,
    use_causal_materials: true,
};
const F_TOP = 40e9;

async function channel(p) {
    const s = buildSolverFromParams(p, () => {});
    s.freq = F_TOP;
    s.ensure_mesh();
    const cached = await s.solve_adaptive({ max_iters: 10, energy_tol: 0.01, param_tol: 0.05,
        max_nodes: 20000, min_converged_passes: 2 });
    const dc = await s.computeAtFrequency(0, cached);
    const sweep = new InterpolatingSweep(s, cached, { tolerance: 0.005 });
    await sweep.run(10e6, F_TOP, {});
    const n = Math.ceil(Math.log10(F_TOP / 10e6) * 40) + 1;
    const grid = Array.from({ length: n }, (_, i) => 10e6 * Math.pow(F_TOP / 10e6, i / (n - 1)));
    return makeRLGCEvaluator(channelTableFromRows([{ freq: 0, result: dc }, ...sweep.buildResults(grid)]));
}

// L - X_KK/w of the odd mode at each f, X_KK(f) = (2f/pi) int (R(f') - R(f)) f' / (f'^2 - f^2) dln f'.
function kkOffset(ev, freqs) {
    const odd = (v, k) => v[k + '11'] - v[k + '12'];
    const u0 = 0, u1 = Math.log(1e19), nU = 40000, du = (u1 - u0) / nU;
    const fs = new Float64Array(nU + 1), Rs = new Float64Array(nU + 1);
    for (let i = 0; i <= nU; i++) { fs[i] = Math.exp(u0 + i * du); Rs[i] = odd(ev(fs[i]), 'R'); }
    return freqs.map(f => {
        const v = ev(f), R0 = odd(v, 'R');
        let s = 0;
        for (let i = 0; i <= nU; i++) {
            const fp = fs[i];
            const g = Math.abs(fp - f) < 1e-9 * f ? 0 : (Rs[i] - R0) * fp / (fp * fp - f * f);
            s += (i === 0 || i === nU ? 0.5 : 1) * g;
        }
        const X = (2 * f / Math.PI) * s * du;
        return odd(v, 'L') - X / (2 * Math.PI * f);
    });
}

const FREQS = [1e8, 3e8, 1e9, 3e9, 1e10, 2e10];
const spread = (d) => (Math.max(...d) - Math.min(...d)) / d[d.length - 1];
for (const [label, o, tol] of [
    ['meshed plating, rough bare bottom', {}, 0.002],
    ['meshed plating of the bulk sigma, rough bare bottom', { plating_sigma: 5.8e7 }, 0.002],
    ['meshed plating on every face, rough plating and bulk',
        { rq: 3e-6, plating_rq: 2e-6, plating_bottom: true }, 0.002],
    ['unplated rough reference', { use_plating: false }, 0.003],
]) {
    const ev = await quiet(() => channel({ ...BASE, ...o }));
    const d = kkOffset(ev, FREQS);
    const sp = spread(d);
    check(`${label}: L - X_KK/w flat over 0.1-20 GHz`, sp < tol,
        `spread ${(100 * sp).toFixed(3)} %, limit ${100 * tol} %; ${d.map(x => (x * 1e9).toFixed(2)).join(' ')} nH/m`);
}
done();
