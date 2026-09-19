// A native solver converted to custom geometry text must solve identically.
//
// native -> solverToGeometryText -> CustomGeometrySolver. With units m the text carries
// the doubles exactly, so the rectangle lists, the domain and the FDM results have to be
// bit-identical. That pins every scalar the backends read besides the lists: a field the
// custom solver derives differently shows up as a changed number here.
//
// The pinWalls form (-inf / inf edges) rebuilds wall-touching rectangles from the domain,
// which may differ by an ulp, so it is compared to 1e-9.
//
// TRI=1 adds the triangular backend on a subset (slow).
import { MicrostripSolver } from '../src/microstrip.js';
import { BroadsideStriplineSolver } from '../src/broadside_stripline.js';
import { CustomGeometrySolver } from '../src/custom_geometry.js';
import { solverToGeometryText } from '../src/custom_geometry_text.js';

let failures = 0;
function check(name, ok, detail = '') {
    console.log(`${ok ? '✓ PASS' : '✗ FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
    if (!ok) failures++;
}

const base = {
    trace_width: 0.3e-3, substrate_height: 0.254e-3, trace_thickness: 35e-6,
    epsilon_r: 3.66, tan_delta: 0.003, sigma_cond: 5.8e7, rq: 0.4e-6, gnd_thickness: 35e-6,
    freq: 5e9, nx: 30, ny: 30, boundaries: ['open', 'open', 'open', 'gnd'],
};
const plating = { sigma: 1.45e7, thickness: 4e-6, rq: 0.2e-6, top: true, sides: true, bottom: false, thick_corners: true };
const bs = {
    trace_width: 0.2e-3, trace_thickness: 18e-6, x_offset: 0, sigma_cond: 5.8e7,
    h_bottom: 0.2e-3, er_bottom: 4.1, tand_bottom: 0.015,
    h_middle: 0.15e-3, er_middle: 3.6, tand_middle: 0.01,
    h_top: 0.2e-3, er_top: 4.1, tand_top: 0.015,
    freq: 5e9, nx: 30, ny: 30, rq: 0, boundaries: ['open', 'open', 'gnd', 'gnd'],
};

const CASES = [
    ['microstrip', MicrostripSolver, base],
    ['diff microstrip', MicrostripSolver, { ...base, trace_spacing: 0.2e-3 }],
    ['embedded microstrip', MicrostripSolver, { ...base, trace_thickness: -35e-6 }],
    ['microstrip + solder mask + plating', MicrostripSolver, { ...base, use_sm: true, plating }],
    ['microstrip + top dielectric + gnd cutout', MicrostripSolver,
        { ...base, top_diel_h: 0.1e-3, top_diel_er: 2.9, top_diel_tand: 0.01, gnd_cut_width: 0.6e-3, gnd_cut_sub_h: 0.2e-3 }],
    ['stripline', MicrostripSolver,
        { ...base, boundaries: ['open', 'open', 'gnd', 'gnd'], enclosure_height: 0.4e-3, epsilon_r_top: 3.66, tan_delta_top: 0.003 }],
    ['GCPW with vias', MicrostripSolver,
        { ...base, use_coplanar_gnd: true, gap: 0.2e-3, via_gap: 0.3e-3, use_vias: true }],
    ['diff GCPW, finite grounds, solder mask', MicrostripSolver,
        { ...base, trace_spacing: 0.2e-3, use_coplanar_gnd: true, gap: 0.2e-3, via_gap: 0.2e-3, use_vias: true,
          coplanar_gnd_width: 1e-3, use_sm: true }],
    ['enclosed microstrip', MicrostripSolver,
        { ...base, boundaries: ['gnd', 'gnd', 'gnd', 'gnd'], enclosure_height: 1e-3, enclosure_width: 3e-3 }],
    ['lid only', MicrostripSolver, { ...base, boundaries: ['open', 'open', 'gnd', 'gnd'], enclosure_height: 1e-3 }],
    ['broadside stripline', BroadsideStriplineSolver, bs],
    ['broadside stripline, offset', BroadsideStriplineSolver, { ...bs, x_offset: 0.1e-3 }],
];

const SOLVE = { max_iters: 4, energy_tol: 0.01, param_tol: 0.05, max_nodes: 6000, min_converged_passes: 2 };

async function quiet(fn) {
    const log = console.log;
    console.log = () => {};
    try { return await fn(); } finally { console.log = log; }
}

function rectKey(o) {
    return [o.x, o.y, o.width, o.height, o.epsilon_r, o.tan_delta, o.is_signal, o.polarity,
        o.plating ? JSON.stringify(Object.entries(o.plating).sort()) : null].join('|');
}

function fingerprint(s, r) {
    return r.modes.map(m => [m.Z0?.re ?? m.Z0, m.eps_eff, m.RLGC.R, m.RLGC.L, m.RLGC.G, m.RLGC.C]).flat();
}

function common(opts, extra = {}) {
    return { sigma_cond: opts.sigma_cond, freq: opts.freq, nx: opts.nx, ny: opts.ny, rq: opts.rq, ...extra };
}

for (const [name, Cls, opts] of CASES) {
    const native = new Cls(opts);
    const text = solverToGeometryText(native);
    let custom;
    try {
        custom = new CustomGeometrySolver({ text, ...common(opts) });
    } catch (e) {
        check(`${name}: converts`, false, e.message);
        continue;
    }
    const sameLists = native.conductors.length === custom.conductors.length
        && native.dielectrics.length === custom.dielectrics.length
        && native.conductors.every((c, i) => rectKey(c) === rectKey(custom.conductors[i]))
        && native.dielectrics.every((d, i) => rectKey(d) === rectKey(custom.dielectrics[i]));
    check(`${name}: identical rectangle lists`, sameLists);
    const sameDomain = native.domain_width === custom.domain_width && native.domain_height === custom.domain_height
        && native.domain_y_min === custom.domain_y_min && !!native.sym_half === !!custom.sym_half
        && native.is_differential === custom.is_differential
        && JSON.stringify(native.boundaries) === JSON.stringify(custom.boundaries);
    check(`${name}: identical domain, symmetry and mode setup`, sameDomain,
        `sym_half ${native.sym_half}/${custom.sym_half}`);

    const rn = await quiet(() => native.solve_adaptive(SOLVE));
    const rc = await quiet(() => custom.solve_adaptive(SOLVE));
    const fn = fingerprint(native, rn), fc = fingerprint(custom, rc);
    const same = fn.length === fc.length && fn.every((v, i) => v === fc[i]);
    check(`${name}: bit-identical FDM results`, same, same ? `Z0 ${fn[0].toFixed(3)}` : `${fn}\n    vs ${fc}`);
    const wn = (rn.warnings || []).map(w => w.reason ?? w.type).sort().join(',');
    const wc = (rc.warnings || []).map(w => w.reason ?? w.type).sort().join(',');
    check(`${name}: same warnings`, wn === wc, `[${wn}] vs [${wc}]`);

    // Wall-pinned form in um: rounding and rebuilt wall edges, compared loosely.
    const pinned = new CustomGeometrySolver({
        text: solverToGeometryText(native, { units: 'um', pinWalls: true }), ...common(opts) });
    const rp = await quiet(() => pinned.solve_adaptive(SOLVE));
    const fp = fingerprint(pinned, rp);
    const rel = Math.max(...fn.map((v, i) => Math.abs(v - fp[i]) / (Math.abs(v) + 1e-300)));
    check(`${name}: pinned-wall um text agrees`, rel < 1e-6, `max rel ${rel.toExponential(2)}`);
}

if (process.env.TRI) {
    for (const [name, Cls, opts] of [CASES[0], CASES[1], CASES[6], CASES[10]]) {
        const native = new Cls({ ...opts, mesh_backend: 'triangular' });
        const custom = new CustomGeometrySolver({
            text: solverToGeometryText(native), ...common(opts), mesh_backend: 'triangular' });
        const rn = await quiet(() => native.solve_adaptive({ ...SOLVE, max_nodes: 20000 }));
        const rc = await quiet(() => custom.solve_adaptive({ ...SOLVE, max_nodes: 20000 }));
        const fn = fingerprint(native, rn), fc = fingerprint(custom, rc);
        const rel = Math.max(...fn.map((v, i) => Math.abs(v - fc[i]) / (Math.abs(v) + 1e-300)));
        check(`${name}: triangular backend agrees`, rel < 1e-6, `max rel ${rel.toExponential(2)}`);
    }
}

console.log(failures === 0 ? '\nALL CUSTOM GEOMETRY EQUIVALENCE TESTS PASSED' : `\n${failures} TEST(S) FAILED`);
process.exit(failures === 0 ? 0 : 1);
