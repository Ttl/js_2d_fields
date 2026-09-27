// DC line parameters of both backends against the exact values of dcLineParameters
// (src/dc_inductance.js), with the grounds of unlimited width (walls, and grounds
// reaching an open domain edge) as ideal returns:
//   1. R at f = 0 on both backends, which drops the wall grounds, over native lines
//      with a wall ground, an enclosure, stacked and standing walls, plated traces,
//      blocks of two metals and a finite ground in open space
//   2. both backends classify the same grounds as walls and as unlimited; a full-width
//      ground drawn as the lowest object under a gnd boundary is the wall itself (no
//      extra slab)
//   3. the free-space perfect-conductor inductance L_pec against the external
//      inductance of both backends
//   4. the full-wave internal inductance at f = 0 is the exact value
//   5. at 100 Hz both backends reach the DC values: R within 0.2% (the wall resistance
//      vanishes as its return current spreads), the full-wave internal L within 4% of
//      the exact one or 1% of the line L (on a finite ground in open space an
//      independent check: no wall, the ground meshed), the quasi-static total L within 2%.
import { CustomGeometrySolver } from '../src/custom_geometry.js';
import { MicrostripSolver } from '../src/microstrip.js';
import { dcLineParameters } from '../src/dc_inductance.js';

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
const rel = (a, b) => Math.abs(a - b) / Math.abs(b);

const MS = { trace_thickness: 35e-6, gnd_thickness: 35e-6, epsilon_r: 4.5, tan_delta: 0, sigma_cond: 5.8e7, freq: 1e6, nx: 10, ny: 10 };
const MS_B = ['open', 'open', 'open', 'gnd'];
const COPLANAR = { use_coplanar_gnd: true, gap: 0.2e-3, via_gap: 0.5e-3, use_vias: true };
const custom = text => B => new CustomGeometrySolver({ text, sigma_cond: 5.8e7, freq: 1e6, nx: 10, ny: 10, ...B });
const cases = {
    microstrip: B => new MicrostripSolver({ ...MS, substrate_height: 0.21e-3, trace_width: 0.35e-3, boundaries: MS_B, ...B }),
    'plated microstrip': B => new MicrostripSolver({ ...MS, substrate_height: 0.21e-3, trace_width: 0.35e-3, boundaries: MS_B,
        plating: { sigma: 1e7, thickness: 4e-6, rq: 0, top: true, sides: true, bottom: false, thick_corners: true }, ...B }),
    GCPW: B => new MicrostripSolver({ ...MS, substrate_height: 1.6e-3, trace_width: 1e-3, ...COPLANAR, boundaries: MS_B, ...B }),
    'differential GCPW': B => new MicrostripSolver({ ...MS, substrate_height: 0.2e-3, trace_width: 0.2e-3, trace_spacing: 0.2e-3,
        ...COPLANAR, boundaries: MS_B, ...B }),
    enclosure: B => new MicrostripSolver({ ...MS, sigma_cond: 1e7, substrate_height: 1.6e-3, trace_width: 3e-3,
        enclosure_width: 10e-3, enclosure_height: 3e-3, boundaries: ['gnd', 'gnd', 'gnd', 'gnd'], ...B }),
    stripline: B => new MicrostripSolver({ ...MS, trace_width: 0.2e-3, substrate_height: 0.2e-3, enclosure_height: 0.235e-3,
        boundaries: ['open', 'open', 'gnd', 'gnd'], epsilon_r_top: 4.5, ...B }),
    'stacked ground slabs': custom('units mm\nbounds open open open gnd\ndiel x=-inf w=inf y=-0.07 h=1.67 er=4.5\n' +
        'gnd x=-inf w=inf y=-0.035 h=0.035\ngnd x=-inf w=inf y=-0.07 h=0.035\nsig+ x=-1.5 w=3 y=1.6 h=0.035'),
    'ground cutout': B => new MicrostripSolver({ ...MS, substrate_height: 0.2e-3, trace_width: 0.3e-3, boundaries: MS_B,
        gnd_cut_width: 1e-3, gnd_cut_sub_h: 0.3e-3, ...B }),
    'ground posts under a lid': custom('units mm\nbounds open open gnd gnd\ndiel x=-2 w=4 y=0 h=1 er=4.5\n' +
        'gnd x=-2.035 w=0.035 y=0 h=inf\ngnd x=2 w=0.035 y=0 h=inf\nsig+ x=-0.5 w=1 y=1 h=0.035'),
    'nickel block on copper': custom('units mm\nbounds open open open gnd\ndiel x=-inf w=inf y=0 h=0.21 er=4.4\n' +
        'sig+ x=-0.175 w=0.35 y=0.21 h=0.031 sigma=5.8e7\nsig+ x=-0.175 w=0.35 y=0.241 h=0.004 sigma=1e7'),
    // A finite plate over the auto ground wall far below: the open side and top walls
    // must not take return current (natural BC), or the backends part ways. The DC
    // return runs in the far wall, a loop as large as the domain, so its inductance at
    // 100 Hz is the domain's, not the free-space value: the backends are checked
    // against each other instead.
    'finite plate over a far wall': custom('units mm\nbounds open open open gnd\ndiel x=-5 w=10 y=0 h=0.21 er=4.4\n' +
        'gnd x=-5 w=10 y=-0.035 h=0.035\nsig+ x=-0.175 w=0.35 y=0.21 h=0.035'),
    'finite ground in open space': custom('units mm\nbounds open open open open\ndiel x=-3 w=6 y=0 h=0.21 er=4.4\n' +
        'gnd x=-3 w=6 y=-0.035 h=0.035\nsig+ x=-0.175 w=0.35 y=0.21 h=0.035'),
};
// A gnd boundary snaps onto a full-width ground at the edge of the stack.
{
    const s = custom('units mm\nbounds open open open gnd\ndiel x=-inf w=inf y=0 h=1.6 er=4.5\n' +
        'gnd x=-inf w=inf y=-0.035 h=0.035\ngnd x=-inf w=inf y=-0.07 h=0.035\nsig+ x=-1.5 w=3 y=1.6 h=0.035')({});
    const gnds = s.conductors.map((c, i) => i).filter(i => !s.conductors[i].is_signal);
    check('ground planes drawn below the substrate are the wall, with no slab added below them',
        gnds.length === 2 && gnds.every(i => s._wall_grounds().has(i)) && Math.abs(s.domain_y_min + 0.07e-3) < 1e-12,
        `${gnds.length} grounds, walls [${[...s._wall_grounds()]}], domain bottom ${(s.domain_y_min * 1e3).toFixed(3)} mm`);
}

const SOLVE = { max_iters: 12, energy_tol: 0.01, param_tol: 0.05, max_nodes: 20000, min_converged_passes: 2 };

async function run(mk, tri, freqs) {
    const s = mk(tri ? { mesh_backend: 'triangular' } : {});
    if (tri) s.tri_opts = { lossMethod: 'auto' };
    return quiet(async () => {
        const r0 = await s.solve_adaptive(SOLVE);
        const at = [];
        for (const f of freqs) at.push((await s.computeAtFrequency(f, r0)).modes);
        return { s, at };
    });
}

const DOMAIN_LOOP = new Set(['finite plate over a far wall']);
for (const [name, mk] of Object.entries(cases)) {
    const [q, t] = await Promise.all([run(mk, false, [0, 100]), run(mk, true, [0, 100])]);
    const wq = [...q.s._wall_grounds()].sort(), wt = [...t.s._wall_grounds()].sort();
    const uq = [...q.s._unlimited_grounds()].sort(), ut = [...t.s._unlimited_grounds()].sort();
    check(`${name}: both backends take the same grounds as walls and as unlimited`,
        wq.join() === wt.join() && uq.join() === ut.join(), `walls [${wq}] / [${wt}], unlimited [${uq}] / [${ut}]`);
    const s = t.s, box = { x_min: -s.domain_width / 2, x_max: s.domain_width / 2, y_min: s.domain_y_min, y_max: s.domain_height };
    t.at[0].forEach((m, k) => {
        const ex = dcLineParameters(s.conductors, m.mode, { sigmaDefault: s.sigma_cond,
            unlimited: s._unlimited_grounds(), walls: s._wall_grounds(), box });
        const mq = q.at[0][k], tag = `${name} ${m.mode}`;
        check(`${tag}: R at DC, quasi-static`, rel(mq.RLGC.R, ex.R) < 1e-9, `${mq.RLGC.R.toFixed(5)} vs ${ex.R.toFixed(5)} ohm/m`);
        check(`${tag}: R at DC, full-wave`, rel(m.RLGC.R, ex.R) < 1e-9, `${m.RLGC.R.toFixed(5)} vs ${ex.R.toFixed(5)} ohm/m`);
        check(`${tag}: free-space L_pec within 1% of the external inductance of both backends`,
            rel(mq.L_external, ex.Lpec) < 0.01 && rel(m.L_external, ex.Lpec) < 0.01,
            `${(ex.Lpec * 1e9).toFixed(1)} vs ${(mq.L_external * 1e9).toFixed(1)} / ${(m.L_external * 1e9).toFixed(1)} nH/m`);
        check(`${tag}: full-wave internal L at DC is the exact value`, rel(m.L_internal, ex.Lint) < 1e-9,
            `${(m.L_internal * 1e9).toFixed(2)} vs ${(ex.Lint * 1e9).toFixed(2)} nH/m`);
        const t100 = t.at[1][k], q100 = q.at[1][k];
        if (DOMAIN_LOOP.has(name)) {
            check(`${tag}: L at 100 Hz, the backends within 3%`, rel(q100.RLGC.L, t100.RLGC.L) < 0.03,
                `${(q100.RLGC.L * 1e9).toFixed(1)} / ${(t100.RLGC.L * 1e9).toFixed(1)} nH/m`);
            check(`${tag}: R at 100 Hz within 0.2% of DC on both backends`,
                rel(q100.RLGC.R, ex.R) < 2e-3 && rel(t100.RLGC.R, ex.R) < 2e-3,
                `${q100.RLGC.R.toFixed(5)} / ${t100.RLGC.R.toFixed(5)} vs ${ex.R.toFixed(5)} ohm/m`);
            return;
        }
        check(`${tag}: R at 100 Hz within 0.2% of DC on both backends`,
            rel(q100.RLGC.R, ex.R) < 2e-3 && rel(t100.RLGC.R, ex.R) < 2e-3,
            `${q100.RLGC.R.toFixed(5)} / ${t100.RLGC.R.toFixed(5)} vs ${ex.R.toFixed(5)} ohm/m`);
        check(`${tag}: full-wave internal L at 100 Hz within 4% of the exact DC value (or 1% of L)`,
            Math.abs(t100.L_internal - ex.Lint) < Math.max(0.04 * ex.Lint, 0.01 * t100.RLGC.L),
            `${(t100.L_internal * 1e9).toFixed(2)} vs ${(ex.Lint * 1e9).toFixed(2)} nH/m`);
        const Lq = q100.L_external + ex.Lint;
        check(`${tag}: quasi-static L at 100 Hz within 2% of the exact DC value`, rel(q100.RLGC.L, Lq) < 0.02,
            `${(q100.RLGC.L * 1e9).toFixed(1)} vs ${(Lq * 1e9).toFixed(1)} nH/m`);
    });
}

if (failures) { console.log(`\n${failures} CHECK(S) FAILED`); process.exit(1); }
console.log('\nALL CHECKS PASSED');
