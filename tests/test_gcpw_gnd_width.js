// Finite coplanar ground width on GCPW and differential GCPW.
//
// coplanar_gnd_width is the width of each coplanar ground measured outward from
// the gap edge; null/'full' keeps the pre-feature layout (grounds and via fences
// run to the domain wall). A finite ground puts the via fence via_gap inside
// both of its edges, so the fence is coplanar_gnd_width - 2 * via_gap wide and
// must be positive (a vialess ground is rejected). Pins the geometry lists, validation,
// the factory/link default, the QS physics (wide ground -> full-width result,
// narrower ground -> higher Z0), the half-domain identity and QS-vs-tri agreement.
//
// Run: node tests/test_gcpw_gnd_width.js
import { MicrostripSolver } from '../src/microstrip.js';
import { buildSolverFromParams } from '../src/solver_factory.js';

let failures = 0;
function check(name, cond, detail = '') {
    console.log(`${cond ? '✓ PASS' : '✗ FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
    if (!cond) failures++;
}
const near = (x, y, tol = 1e-12) => Math.abs(x - y) < tol;
const rel = (a, b) => Math.abs(a - b) / Math.max(Math.abs(a), Math.abs(b), 1e-30);
const quiet = async (fn) => {
    const log = console.log, warn = console.warn;
    console.log = () => {}; console.warn = () => {};
    try { return await fn(); } finally { console.log = log; console.warn = warn; }
};
const APP = { max_iters: 10, energy_tol: 0.01, param_tol: 0.05, max_nodes: 20000, min_converged_passes: 2 };

const W = 0.35e-3, H = 0.21e-3, T = 35e-6, GAP = 0.1e-3, VIA_GAP = 0.1e-3;
const BASE = {
    trace_width: W, substrate_height: H, trace_thickness: T, gnd_thickness: 35e-6,
    epsilon_r: 4.4, tan_delta: 0.02, sigma_cond: 5.8e7, freq: 1e9, nx: 30, ny: 30,
    use_coplanar_gnd: true, gap: GAP, via_gap: VIA_GAP, use_vias: true,
    boundaries: ['open', 'open', 'open', 'gnd'],
};
const build = (o = {}) => new MicrostripSolver({ ...BASE, ...o });
const grounds = (s) => s.conductors.filter(c => !c.is_signal);
// Coplanar grounds sit on the trace layer, via fences span bottom ground to trace top
const coplanar = (s) => grounds(s).filter(c => near(c.y_min, H) && near(c.y_max, H + T));
const fences = (s) => grounds(s).filter(c => near(c.y_min, 0) && near(c.y_max, H + T));
const right = (list) => list.find(c => c.x_min > 0), left = (list) => list.find(c => c.x_max < 0);

// ---------- default: full-width grounds (pre-feature layout) ----------
for (const [label, o] of [['omitted', {}], ["'full'", { coplanar_gnd_width: 'full' }], ['null', { coplanar_gnd_width: null }]]) {
    const s = build(o);
    const xw = s.domain_width / 2;
    check(`${label}: coplanar grounds run to the domain wall`,
        near(right(coplanar(s)).x_max, xw) && near(left(coplanar(s)).x_min, -xw));
    check(`${label}: via fences run from gap + via_gap to the wall`,
        near(right(fences(s)).x_min, W / 2 + GAP + VIA_GAP) && near(right(fences(s)).x_max, xw)
        && near(left(fences(s)).x_max, -W / 2 - GAP - VIA_GAP) && near(left(fences(s)).x_min, -xw));
}

// ---------- finite ground, single-ended ----------
{
    const GW = 1e-3;
    const s = build({ coplanar_gnd_width: GW });
    const gr = right(coplanar(s)), gl = left(coplanar(s));
    check('finite: coplanar ground spans gap edge .. gap edge + width',
        near(gr.x_min, W / 2 + GAP) && near(gr.x_max, W / 2 + GAP + GW)
        && near(gl.x_max, -W / 2 - GAP) && near(gl.x_min, -W / 2 - GAP - GW));
    const fr = right(fences(s)), fl = left(fences(s));
    check('finite: via fence via_gap inside both ground edges',
        near(fr.x_min, W / 2 + GAP + VIA_GAP) && near(fr.x_max, W / 2 + GAP + GW - VIA_GAP)
        && near(fl.x_max, -W / 2 - GAP - VIA_GAP) && near(fl.x_min, -W / 2 - GAP - GW + VIA_GAP));
    check('finite: via thickness = width - 2 * via_gap', near(fr.x_max - fr.x_min, GW - 2 * VIA_GAP));
    check('finite: bare substrate between the ground edge and the wall',
        s.domain_width / 2 - gr.x_max > 15 * H - 1e-12);
    check('finite: domain narrower than 1.5x the active width or the microstrip clearance',
        s.domain_width <= Math.max(1.5, 1) * (W + 2 * (GAP + GW)) + 2 * Math.max(8 * W, 15 * H) + 1e-12);
    check('finite: x=0 half-domain symmetry still detected', s.sym_half === true);
    check('finite: conductor count unchanged (bottom gnd, 2 fences, trace, 2 grounds)', s.conductors.length === 6);
}

// ---------- finite ground, differential ----------
{
    const GW = 0.8e-3, SP = 0.2e-3;
    const s = build({ coplanar_gnd_width: GW, trace_spacing: SP });
    const outer = W + SP / 2;   // outer trace edge
    const gr = right(coplanar(s)), fr = right(fences(s));
    check('diff finite: ground from outer gap edge, width GW',
        near(gr.x_min, outer + GAP) && near(gr.x_max, outer + GAP + GW));
    check('diff finite: fence via_gap inside both ground edges',
        near(fr.x_min, outer + GAP + VIA_GAP) && near(fr.x_max, outer + GAP + GW - VIA_GAP));
    check('diff finite: mirror-symmetric', near(left(coplanar(s)).x_min, -gr.x_max) && near(left(fences(s)).x_min, -fr.x_max));
}

// ---------- via thickness zero / negative ----------
{
    let err = null;
    try { build({ coplanar_gnd_width: 2 * VIA_GAP }); } catch (e) { err = e.message; }
    check('width == 2 * via_gap (vialess): rejected', !!err && /via thickness/.test(err), err ? err.split('\n').pop().trim() : 'no error');
    const s = build({ coplanar_gnd_width: 2 * VIA_GAP + 1e-6 });
    check('width just above 2 * via_gap: 1 um fence on both sides', fences(s).length === 2
        && near(right(fences(s)).x_max - right(fences(s)).x_min, 1e-6, 1e-15));
    err = null;
    try { build({ coplanar_gnd_width: 2 * VIA_GAP - 1e-6 }); } catch (e) { err = e.message; }
    check('width < 2 * via_gap: rejected as negative via thickness', !!err && /via thickness/.test(err));
    err = null;
    try { build({ coplanar_gnd_width: 0 }); } catch (e) { err = e.message; }
    check('width 0: rejected', !!err && /coplanar_gnd_width/.test(err));
    err = null;
    try { build({ coplanar_gnd_width: -1e-3 }); } catch (e) { err = e.message; }
    check('negative width: rejected', !!err && /coplanar_gnd_width/.test(err));
}

// ---------- enclosure fit ----------
{
    const GW = 1e-3, active = W + 2 * (GAP + GW);
    let err = null;
    try { build({ coplanar_gnd_width: GW, enclosure_width: active - 1e-6, boundaries: ['gnd', 'gnd', 'open', 'gnd'] }); }
    catch (e) { err = e.message; }
    check('enclosure narrower than trace + gaps + grounds: rejected', !!err && /exceeds enclosure/.test(err));
    const s = build({ coplanar_gnd_width: GW, enclosure_width: active + 0.5e-3, boundaries: ['gnd', 'gnd', 'open', 'gnd'] });
    const gr = right(coplanar(s));
    check('enclosure wider: ground ends 0.25 mm inside the side wall',
        near(s.domain_width / 2 - s.t_gnd - gr.x_max, 0.25e-3));
    const flush = build({ coplanar_gnd_width: GW, enclosure_width: active, boundaries: ['gnd', 'gnd', 'open', 'gnd'] });
    check('enclosure exactly the active width: ground touches the side wall',
        near(flush.domain_width / 2 - flush.t_gnd, right(coplanar(flush)).x_max));
    // Full-width grounds in an enclosure keep the old fit rule (trace + gaps + via_gap)
    // and the old layout: ground and fence run to the domain edge, through the side slab
    const full = build({ enclosure_width: W + 2 * (GAP + VIA_GAP) + 0.1e-3, boundaries: ['gnd', 'gnd', 'open', 'gnd'] });
    check('full-width grounds in an enclosure run to the domain edge',
        near(right(coplanar(full)).x_max, full.domain_width / 2) && near(right(fences(full)).x_max, full.domain_width / 2));
}

// ---------- solder mask ----------
{
    const GW = 1e-3, SM = { use_sm: true, sm_t_sub: 20e-6, sm_t_trace: 20e-6, sm_t_side: 20e-6, sm_er: 3.5, sm_tand: 0.02 };
    const s = build({ coplanar_gnd_width: GW, ...SM });
    const masks = s.dielectrics.filter(d => d.epsilon_r === 3.5);
    const gr = right(coplanar(s));
    check('sm: mask on the ground top spans exactly the ground',
        !!masks.find(d => near(d.y_min, H + T) && near(d.x_min, gr.x_min) && near(d.x_max, gr.x_max)));
    check('sm: side band on the ground outer edge',
        !!masks.find(d => near(d.x_min, gr.x_max) && near(d.x_max, gr.x_max + SM.sm_t_side) && near(d.y_min, H)));
    check('sm: substrate mask from the side band to the wall',
        !!masks.find(d => near(d.x_min, gr.x_max + SM.sm_t_side) && near(d.x_max, s.domain_width / 2) && near(d.y_min, H) && near(d.y_max, H + SM.sm_t_sub)));
    check('sm: no mask rect overlaps a conductor interior',
        !masks.some(d => s.conductors.some(c => d.x_min < c.x_max - 1e-12 && c.x_min < d.x_max - 1e-12
                                              && d.y_min < c.y_max - 1e-12 && c.y_min < d.y_max - 1e-12)));
    const full = build(SM);
    check('sm full-width: mask on the ground top runs to the wall',
        !!full.dielectrics.find(d => d.epsilon_r === 3.5 && near(d.y_min, H + T) && near(d.x_max, full.domain_width / 2)));
    check('sm: differential finite ground builds', !!build({ coplanar_gnd_width: GW, trace_spacing: 0.2e-3, ...SM }));
}

// ---------- factory / link default ----------
{
    const P = { tl_type: 'gcpw', w: W, h: H, t: T, er: 4.4, tand: 0.02, sigma: 5.8e7, gap: GAP, via_gap: VIA_GAP,
        freq: 1e9, nx: 30, ny: 30, rq: 0, mesh_backend: 'rectilinear' };
    const errs = [];
    const s0 = buildSolverFromParams(P, e => errs.push(e));
    check('factory: gnd_width absent (old links) -> full-width grounds', !!s0 && s0.coplanar_gnd_width === null, errs.join(' '));
    const s1 = buildSolverFromParams({ ...P, gnd_width: 'full' }, e => errs.push(e));
    check("factory: gnd_width 'full' (empty UI input) -> full-width grounds", !!s1 && s1.coplanar_gnd_width === null);
    const s2 = buildSolverFromParams({ ...P, gnd_width: 1e-3 }, e => errs.push(e));
    check('factory: numeric gnd_width -> finite grounds', !!s2 && s2.coplanar_gnd_width === 1e-3
        && near(right(coplanar(s2)).x_max, W / 2 + GAP + 1e-3));
    const s3 = buildSolverFromParams({ ...P, tl_type: 'diff_gcpw', trace_spacing: 0.2e-3, gnd_width: 1e-3 }, e => errs.push(e));
    check('factory: diff_gcpw takes gnd_width too', !!s3 && s3.coplanar_gnd_width === 1e-3);
    const bad = [];
    check('factory: negative via thickness surfaces as an error',
        buildSolverFromParams({ ...P, gnd_width: 0.1e-3 }, e => bad.push(e)) === null && /via thickness/.test(bad.join(' ')));
}

// ---------- QS physics ----------
{
    const solve = async (o, extra = {}) => {
        const s = build(o);
        const r = await quiet(() => s.solve_adaptive({ ...APP, ...extra }));
        return { s, m: r.modes[0] };
    };
    const full = await solve({});
    const wide = await solve({ coplanar_gnd_width: 5e-3 });
    check('QS: 5 mm grounds reproduce the full-width Z0 (< 1%)', rel(wide.m.Z0, full.m.Z0) < 0.01,
        `${wide.m.Z0.toFixed(3)} vs ${full.m.Z0.toFixed(3)} ohm`);
    check('QS: 5 mm grounds reproduce the full-width eps_eff (< 1%)', rel(wide.m.eps_eff, full.m.eps_eff) < 0.01,
        `${wide.m.eps_eff.toFixed(4)} vs ${full.m.eps_eff.toFixed(4)}`);
    // The ground is at 0 V whatever its width, so only the fringe field past its
    // outer edge changes: a strip narrower than the substrate height reads a
    // higher Z0 than a wide ground (~0.5% here), wider strips are within mesh noise.
    const zWide = (await solve({ via_gap: 0.02e-3, coplanar_gnd_width: 3e-3 })).m.Z0;
    const zStrip = (await solve({ via_gap: 0.02e-3, coplanar_gnd_width: 0.06e-3 })).m.Z0;
    check('QS: 60 um ground strips raise Z0 above 3 mm grounds (> 0.2%)', (zStrip - zWide) / zWide > 0.002,
        `${zStrip.toFixed(3)} vs ${zWide.toFixed(3)} ohm`);
    // The strips still lower Z0 well below the bare microstrip
    const ms = new MicrostripSolver({ ...BASE, use_coplanar_gnd: false, use_vias: false });
    const rMs = await quiet(() => ms.solve_adaptive({ ...APP }));
    check('QS: narrow ground strips stay below the microstrip Z0', zStrip < rMs.modes[0].Z0,
        `${zStrip.toFixed(2)} vs microstrip ${rMs.modes[0].Z0.toFixed(2)} ohm`);
    // Half-domain identity
    const halfS = build({ coplanar_gnd_width: 0.6e-3 }), fullS = build({ coplanar_gnd_width: 0.6e-3, symmetry: false });
    const rh = await quiet(() => halfS.solve_adaptive({ ...APP })), rf = await quiet(() => fullS.solve_adaptive({ ...APP }));
    check('QS: half-domain solve == full-domain solve', halfS.sym_half && !fullS.sym_half
        && rel(rh.modes[0].RLGC.C, rf.modes[0].RLGC.C) < 1e-6 && rel(rh.modes[0].RLGC.R, rf.modes[0].RLGC.R) < 1e-6,
        `C ${rh.modes[0].RLGC.C.toExponential(6)} vs ${rf.modes[0].RLGC.C.toExponential(6)}`);
}

// ---------- full-wave backend ----------
{
    const o = { coplanar_gnd_width: 0.6e-3 };
    const qs = build(o);
    const rq = await quiet(() => qs.solve_adaptive({ ...APP }));
    const tri = build({ ...o, mesh_backend: 'triangular' });
    tri.use_causal_materials = false;
    tri.tri_opts = { lossMethod: 'auto' };
    let rt = null, err = '';
    try { rt = await quiet(() => tri.solve_adaptive({ ...APP })); } catch (e) { err = e.message; }
    check('tri: finite-ground GCPW meshes and solves', !!rt, err);
    if (rt) {
        check('tri: Z0 agrees with QS (< 3%)', rel(rt.modes[0].Z0, rq.modes[0].Z0) < 0.03,
            `tri ${rt.modes[0].Z0.toFixed(3)} vs qs ${rq.modes[0].Z0.toFixed(3)} ohm`);
        check('tri: eps_eff agrees with QS (< 3%)', rel(rt.modes[0].eps_eff, rq.modes[0].eps_eff) < 0.03,
            `tri ${rt.modes[0].eps_eff.toFixed(4)} vs qs ${rq.modes[0].eps_eff.toFixed(4)}`);
        const warns = [...(rt.warnings || []), ...(tri.modeWarnings || [])].map(w => w.type || w);
        check('tri: grounds are bonded through the fence (no floating-grounds warning)', !warns.includes('floating-grounds'), warns.join(','));
    }
}

console.log(failures === 0 ? '\nALL GCPW GROUND WIDTH TESTS PASSED' : `\n${failures} TEST(S) FAILED`);
process.exit(failures === 0 ? 0 : 1);
