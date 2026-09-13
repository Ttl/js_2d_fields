// Quasi-static (rectilinear FDM) edge cases outside the fuzzer's draw ranges.
// Each case pins one item from docs/untested_edge_cases.md (Q1..Q6):
//
//   Q1 uncovered air nodes must be lossless: a ground cutout with gnd_cut_sub_h = 0
//      leaves a void no dielectric rect covers, and its loss tangent must be 0, not
//      the array's fill value. Cross-checked against the full-wave alpha_d.
//   Q2 internal inductance must stay bounded below the skin regime: L_int(f) is
//      finite, non-increasing from 1 kHz to 10 GHz, and within a factor of two of
//      the full-wave MQS value at 1 kHz (the FDM has no lateral-crowding model, so
//      it is not expected to match closely; it must not diverge as 1/sqrt(f)).
//   Q3 the skin-transition warning fires for embedded (negative-thickness) traces
//      exactly as it does for surface traces of the same |t|.
//   Q4 solder-mask bracket lines depend on the mask, not on the trace: R must be
//      continuous when the trace thickness crosses the mask thickness.
//   Q5 degenerate geometries are rejected at construction: gap = 0 on a GCPW and
//      enclosure_height == trace_thickness on a stripline both paint ground over
//      trace nodes if allowed through.
//   Q6 the mesher's conductor-dimension helpers use |height| so an embedded trace
//      gets the same bracket treatment as a surface trace.
//
// Run: node tests/test_qs_edge_cases.js
import { MicrostripSolver } from '../src/microstrip.js';

let pass = 0, fail = 0;
function check(name, cond, detail = '') {
    console.log(`  ${cond ? 'PASS' : 'FAIL'}: ${name}${detail ? ` (${detail})` : ''}`);
    cond ? pass++ : fail++;
}
const quiet = async (fn) => {
    const log = console.log, warn = console.warn;
    console.log = () => {}; console.warn = () => {};
    try { return await fn(); } finally { console.log = log; console.warn = warn; }
};
const APP = { max_iters: 10, energy_tol: 0.01, param_tol: 0.05, max_nodes: 20000, min_converged_passes: 2 };
const rel = (a, b) => Math.abs(a - b) / Math.max(Math.abs(a), Math.abs(b), 1e-30);

const BASE = {
    trace_width: 0.2e-3, substrate_height: 0.3e-3, trace_thickness: 35e-6, gnd_thickness: 35e-6,
    epsilon_r: 4.4, tan_delta: 0.02, sigma_cond: 5.8e7, freq: 1e9, rq: 0, nx: 30, ny: 30,
    boundaries: ['open', 'open', 'open', 'gnd'],
};
function make(opts, backend = 'rectilinear') {
    const s = new MicrostripSolver({ ...BASE, ...opts, mesh_backend: backend });
    if (backend === 'triangular') s.tri_opts = { lossMethod: 'auto' };
    return s;
}
async function solve(opts, backend) {
    const s = make(opts, backend);
    const r = await quiet(() => s.solve_adaptive({ ...APP }));
    return { s, r, m: r.modes[0] };
}

// Q1: ground cutout with zero substrate extension.
{
    console.log('Q1 ground cutout with gnd_cut_sub_h = 0');
    const cut = { gnd_cut_width: 0.3e-3, gnd_cut_sub_h: 0 };
    const { s, m } = await solve(cut, 'rectilinear');
    let maxTand = 0, offenders = 0;
    const tandMax = Math.max(s.tan_delta ?? 0, ...s.dielectrics.map(d => d.tan_delta));
    for (let i = 0; i < s.y.length; i++)
        for (let j = 0; j < s.x.length; j++) {
            if (s.conductor_mask[i][j]) continue;
            maxTand = Math.max(maxTand, s.tand[i][j]);
            if (s.tand[i][j] > tandMax + 1e-12) offenders++;
        }
    check('no non-conductor node exceeds the largest dielectric tan_delta', offenders === 0,
        `${offenders} nodes, max tand ${maxTand}`);
    const fw = await solve(cut, 'triangular');
    check('alpha_d agrees with full-wave within 5%', rel(m.alpha_d, fw.m.alpha_d) < 0.05,
        `fdm ${m.alpha_d.toFixed(3)} vs tri ${fw.m.alpha_d.toFixed(3)} dB/m`);
    check('Z0 agrees with full-wave within 2%', rel(m.Z0, fw.m.Z0) < 0.02,
        `fdm ${m.Z0.toFixed(2)} vs tri ${fw.m.Z0.toFixed(2)}`);
}

// Q2: internal inductance below the skin regime.
{
    console.log('Q2 internal inductance from 1 kHz to 10 GHz');
    const geo = { substrate_height: 0.1e-3, epsilon_r: 4, tan_delta: 0.01, freq: 10e9 };
    const { s, r } = await solve(geo, 'rectilinear');
    const freqs = [1e3, 1e4, 1e5, 3e5, 1e6, 3e6, 1e7, 3e7, 1e8, 1e9, 1e10];
    const lint = [];
    for (const f of freqs) {
        const rf = await quiet(() => s.computeAtFrequency(f, r));
        lint.push(rf.modes[0].L_internal);
    }
    check('L_int finite at every frequency', lint.every(v => Number.isFinite(v) && v >= 0),
        lint.map(v => (v * 1e9).toFixed(1)).join(' '));
    let monotone = true;
    for (let k = 1; k < lint.length; k++) if (lint[k] > lint[k - 1] * 1.001) monotone = false;
    check('L_int non-increasing with frequency', monotone,
        `nH/m: ${lint.map(v => (v * 1e9).toFixed(1)).join(' ')}`);
    const fw = await solve(geo, 'triangular');
    const fw1k = await quiet(() => fw.s.computeAtFrequency(1e3, fw.r));
    const ratio = lint[0] / fw1k.modes[0].L_internal;
    check('L_int(1 kHz) within a factor of two of full-wave MQS', ratio > 0.5 && ratio < 2,
        `fdm ${(lint[0] * 1e9).toFixed(1)} vs tri ${(fw1k.modes[0].L_internal * 1e9).toFixed(1)} nH/m`);
}

// Q3: embedded trace skin-transition warning.
{
    console.log('Q3 skin-transition warning for negative trace thickness');
    // 10 MHz: skin depth 21 um in copper, t = 20 um, well inside the transition.
    const hasSkinWarn = r => (r.warnings || []).some(w => w.reason === 'skin-transition');
    const pos = await solve({ trace_thickness: 20e-6, freq: 1e7 }, 'rectilinear');
    const neg = await solve({ trace_thickness: -20e-6, freq: 1e7 }, 'rectilinear');
    check('surface trace at delta ~ t warns', hasSkinWarn(pos.r));
    check('embedded trace at delta ~ |t| warns', hasSkinWarn(neg.r));
    const posHi = await solve({ trace_thickness: 20e-6, freq: 20e9 }, 'rectilinear');
    const negHi = await solve({ trace_thickness: -20e-6, freq: 20e9 }, 'rectilinear');
    check('neither warns at 20 GHz', !hasSkinWarn(posHi.r) && !hasSkinWarn(negHi.r));
}

// Q4: R continuity when the trace thickness crosses the solder-mask thickness.
{
    console.log('Q4 solder-mask bracket lines across t = sm thickness');
    const sm = { use_sm: true, sm_t_sub: 20e-6, sm_t_trace: 20e-6, sm_t_side: 20e-6, sm_er: 3.5, sm_tand: 0.02 };
    const a = await solve({ ...sm, trace_thickness: 20e-6 * (1 + 1e-4) }, 'rectilinear');
    const b = await solve({ ...sm, trace_thickness: 20e-6 * (1 - 1e-4) }, 'rectilinear');
    // R at the app budget carries ~2% mesh noise across unrelated thickness changes
    // (R is not in the convergence gate); the 1% gate catches the bracket-line and
    // sliver-cell jumps (2-7% before the fix) without tripping on that noise.
    check('R continuous across the threshold (< 1%)', rel(a.m.RLGC.R, b.m.RLGC.R) < 0.01,
        `${a.m.RLGC.R.toFixed(3)} vs ${b.m.RLGC.R.toFixed(3)} ohm/m`);
    check('C continuous across the threshold (< 0.05%)', rel(a.m.RLGC.C, b.m.RLGC.C) < 5e-4,
        `${(a.m.RLGC.C * 1e12).toFixed(3)} vs ${(b.m.RLGC.C * 1e12).toFixed(3)} pF/m`);
    const thin = make({ ...sm, trace_thickness: 10e-6 });
    thin.ensure_mesh();
    const ySm = thin.y_sub_end + sm.sm_t_sub;
    const near = thin.y.filter(y => Math.abs(y - ySm) < sm.sm_t_sub / 5 && Math.abs(y - ySm) > 1e-12).length;
    check('mask thicker than the trace still gets bracket lines', near >= 2, `${near} lines within sm_t/5 of the mask top`);
}

// Q5: degenerate geometries rejected at construction.
{
    console.log('Q5 validation of touching conductors');
    const throws = (opts) => { try { make(opts); return false; } catch (e) { return /gap|enclosure_height/.test(e.message); } };
    const gcpw = { use_coplanar_gnd: true, use_vias: true, via_gap: 0.2e-3 };
    check('GCPW gap = 0 rejected', throws({ ...gcpw, gap: 0 }));
    check('GCPW gap > 0 accepted', !throws({ ...gcpw, gap: 0.1e-3 }));
    const sl = { epsilon_r_top: 4.4, tan_delta_top: 0.02, boundaries: ['open', 'open', 'gnd', 'gnd'] };
    check('stripline enclosure_height == t rejected', throws({ ...sl, enclosure_height: 35e-6 }));
    check('stripline enclosure_height < t rejected', throws({ ...sl, enclosure_height: 30e-6 }));
    check('stripline enclosure_height > t accepted', !throws({ ...sl, enclosure_height: 0.1e-3 }));
}

// Q6: mesher helpers on an embedded trace.
{
    console.log('Q6 mesher conductor-dimension helpers with negative height');
    const s = make({ trace_thickness: -20e-6 });
    s.ensure_mesh();
    check('_min_conductor_dimension positive', s.mesher._min_conductor_dimension() > 0,
        `${s.mesher._min_conductor_dimension()}`);
    // The trace faces at y_trace_start (top) and y_trace_start - 20 um (bottom) must each
    // carry an inside and an outside bracket line at |t|/20.
    const off = 20e-6 / 20;
    const top = s.y_trace_start, bot = s.y_trace_start - 20e-6;
    const has = v => s.y.some(y => Math.abs(y - v) < 1e-12);
    check('bracket lines present around both faces',
        has(top + off) && has(top - off) && has(bot + off) && has(bot - off));
}

console.log(`\n${pass} passed, ${fail} failed`);
process.exit(fail ? 1 : 0);
