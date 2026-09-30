// Trapezoids and n-gons on the full-wave solver, each against a reference:
//   1. a trapezoid with a vanishing angle solves like the rectangle
//   2. a square drawn as an n-gon (a polygon path through the mesher and the MQS loss)
//      solves like the same square drawn as a rectangle
//   3. the coax converted to n-gons against the closed forms, and the native coax
//   4. a differential pair of etched traces: half domain against full domain
//   5. per-face plating on a trapezoid: plating=all equals top,sides,bottom, a bottom
//      face differs from the top face
//   6. a confocal elliptic coax against its closed form
//   7. rounded trace corners: a vanishing radius is the sharp trace, a real one lowers R
//   8. twinax (stadium shell around two wires): half domain against full domain
import { CustomGeometrySolver } from '../src/custom_geometry.js';
import { buildSolverFromParams } from '../src/solver_factory.js';
import { solverToGeometryText } from '../src/custom_geometry_text.js';
import { check, quiet, rel, APP, done } from './helpers.js';

const pct = v => `${(100 * v).toFixed(2)}%`;
const SOLVE = { ...APP, max_iters: 8 };
const modes = r => r.modes.map(m => ({ Z0: m.Z0?.re ?? m.Z0, eps: m.eps_eff, R: m.RLGC.R, L: m.RLGC.L, lossVia: m.lossVia }));
async function solve(text, extra = {}) {
    const s = new CustomGeometrySolver({ text, freq: 1e9, mesh_backend: 'triangular', ...extra });
    const r = await quiet(() => s.solve_adaptive(SOLVE));
    return { m: modes(r), warnings: [...(r.warnings || []), ...(s.modeWarnings || [])] };
}
function agree(name, a, b, zTol, rTol) {
    a.forEach((m, i) => {
        const n = a.length > 1 ? `${name} mode ${i}` : name;
        check(`${n}: Z0 and eps_eff`, rel(m.Z0, b[i].Z0) < zTol && rel(m.eps, b[i].eps) < zTol,
            `Z0 ${m.Z0.toFixed(3)} / ${b[i].Z0.toFixed(3)}, eps ${m.eps.toFixed(4)} / ${b[i].eps.toFixed(4)}`);
        check(`${n}: R`, rel(m.R, b[i].R) < rTol, `${m.R.toFixed(3)} / ${b[i].R.toFixed(3)} ohm/m, ${pct(rel(m.R, b[i].R))}`);
    });
}

const MS = trace => `units mm\nbounds open open open gnd\ndiel x=-inf w=inf y=0 h=0.2 er=4.3 tand=0.02\n${trace}\n`;

// --- 1. Trapezoid with a vanishing angle ---
const rect = await solve(MS('sig+ x=-0.15 y=0.2 w=0.3 h=0.035'));
{
    const trap = await solve(MS('sig+ trap x=-0.15 y=0.2 w=0.3 h=0.035 angle=0.0001'));
    check('trapezoid solves with the MQS loss', trap.m[0].lossVia === 'mqs', trap.m[0].lossVia);
    agree('trapezoid at 1e-4 degrees = rectangle', trap.m, rect.m, 1e-3, 0.005);
    const etched = await solve(MS('sig+ trap x=-0.15 y=0.2 w=0.3 h=0.035 angle=30'));
    check('etched trace: narrower top raises Z0 and R', etched.m[0].Z0 > rect.m[0].Z0 && etched.m[0].R > rect.m[0].R,
        `Z0 ${etched.m[0].Z0.toFixed(2)} vs ${rect.m[0].Z0.toFixed(2)}, R ${etched.m[0].R.toFixed(2)} vs ${rect.m[0].R.toFixed(2)}`);
}

// --- 2. Square as an n-gon ---
{
    // n = 4 turned by 45 degrees is an axis-aligned square of side r*sqrt(2).
    const a = 0.2, r = a / Math.SQRT2;
    const sq = await solve(MS(`sig+ ngon x=0 y=${0.3 + a / 2} r=${r} n=4 rot=45`));
    const box = await solve(MS(`sig+ x=${-a / 2} y=0.3 w=${a} h=${a}`));
    agree('square n-gon = square rectangle', sq.m, box.m, 3e-3, 0.02);
}

// --- 3. Coax as n-gons ---
{
    const p = { tl_type: 'coax', coax_d: 0.92e-3, coax_D: 2.95e-3, coax_er: 2.1, coax_tand: 2e-4, coax_sigma: 5.8e7,
        rq: 0, freq: 1e9, mesh_backend: 'fullwave_mqs', use_plating: false };
    const text = solverToGeometryText(buildSolverFromParams(p, () => {}), { units: 'mm' });
    const c = await solve(text, { sigma_cond: 5.8e7 });
    const eta = 376.730313668 / Math.sqrt(2.1);
    const a = 0.46e-3, b = 1.475e-3;
    const Z0 = eta / (2 * Math.PI) * Math.log(b / a);
    const Rs = Math.sqrt(Math.PI * 1e9 * 4e-7 * Math.PI / 5.8e7);
    // The centre conductor's curvature adds about delta/(2a) to its share.
    const R = Rs / (2 * Math.PI) * (1 / a + 1 / b);
    check('coax n-gons: Z0 = closed form', rel(c.m[0].Z0, Z0) < 1e-3, `${c.m[0].Z0.toFixed(3)} vs ${Z0.toFixed(3)}`);
    check('coax n-gons: R = closed form (MQS)', c.m[0].lossVia === 'mqs' && rel(c.m[0].R, R) < 0.01,
        `${c.m[0].R.toFixed(4)} vs ${R.toFixed(4)}, ${pct(rel(c.m[0].R, R))}`);
    check('coax n-gons: shield band not capped', !c.warnings.some(w => w.type === 'mqs-band-capped'));
}

// --- 4. Etched pair, half domain against full domain ---
{
    const pair = MS('sig+ trap x=0.1 y=0.2 w=0.3 h=0.035 angle=30 mirror=1');
    const half = await solve(pair), full = await solve(pair, { symmetry: false });
    agree('etched pair, half vs full domain', half.m, full.m, 5e-3, 0.01);
}

// --- 5. Per-face plating on a trapezoid ---
{
    const plated = faces => MS(`plating sigma=1e7 t=3um\nsig+ trap x=-0.15 y=0.2 w=0.3 h=0.035 angle=30 plating=${faces}`);
    const all = await solve(plated('all')), each = await solve(plated('top,sides,bottom'));
    check('plating=all equals every face', rel(all.m[0].R, each.m[0].R) < 1e-9, `${all.m[0].R} / ${each.m[0].R}`);
    const top = await solve(plated('top')), bottom = await solve(plated('bottom'));
    const bare = await solve(MS('sig+ trap x=-0.15 y=0.2 w=0.3 h=0.035 angle=30'));
    // A microstrip carries most of its current on the face towards the ground.
    check('plating: every face plated > bottom > top > bare', all.m[0].R > bottom.m[0].R && bottom.m[0].R > top.m[0].R
        && top.m[0].R > bare.m[0].R,
        `${all.m[0].R.toFixed(2)} > ${bottom.m[0].R.toFixed(2)} > ${top.m[0].R.toFixed(2)} > ${bare.m[0].R.toFixed(2)}`);
}

// --- 6. Confocal elliptic coax ---
{
    // Inner and outer ellipses with the same foci: Z0 = eta/(2 pi) ln((a2+b2)/(a1+b1)).
    const a1 = 0.6, b1 = 0.3, a2 = 1.5, b2 = Math.sqrt(a2 * a2 - (a1 * a1 - b1 * b1)), t = 0.15;
    const text = `units mm\nbounds open open open open\ndomain ${-1.1 * (a2 + t)} ${1.1 * (a2 + t)} ${-1.1 * (b2 + t)} ${1.1 * (b2 + t)}\n` +
        `diel ellipse x=0 y=0 rx=${a2} ry=${b2} n=192 er=2.1\nsig+ ellipse x=0 y=0 rx=${a1} ry=${b1} n=96\n` +
        `gnd ellipse x=0 y=0 rx=${a2 + t} ry=${b2 + t} rx_in=${a2} ry_in=${b2} n=192\n`;
    const e = await solve(text);
    const Z0 = 376.730313668 / Math.sqrt(2.1) / (2 * Math.PI) * Math.log((a2 + b2) / (a1 + b1));
    check('confocal elliptic coax: Z0 = closed form', rel(e.m[0].Z0, Z0) < 3e-3, `${e.m[0].Z0.toFixed(3)} vs ${Z0.toFixed(3)}`);
}

// --- 7. Rounded trace corners ---
{
    const tiny = await solve(MS('sig+ x=-0.15 y=0.2 w=0.3 h=0.035 radius=0.01um'));
    agree('corner radius 10 nm = sharp rectangle', tiny.m, rect.m, 1e-3, 0.01);
    const round = await solve(MS('sig+ x=-0.15 y=0.2 w=0.3 h=0.035 radius=10um radius_bottom=0'));
    check('rounded top corners lower R', round.m[0].R < rect.m[0].R && round.m[0].R > 0.95 * rect.m[0].R,
        `${round.m[0].R.toFixed(3)} vs ${rect.m[0].R.toFixed(3)} ohm/m`);
}

// --- 8. Twinax, half against full domain ---
{
    const tw = `units mm\nd = 0.4; di = 1.2; s = 1.25; g = 0.02; t_sh = 0.03\nhw = s/2+di/2+g; hh = di/2+g\n` +
        'bounds open open open open\ndomain -1.2*(hw+t_sh) 1.2*(hw+t_sh) -1.5*(hh+t_sh) 1.5*(hh+t_sh)\n' +
        'diel ngon x=s/2 y=0 r=di/2 n=64 er=2.1 tand=0.0003 mirror=1\nsig+ ngon x=s/2 y=0 r=d/2 n=32 mirror=1\n' +
        'gnd x=-hw-t_sh w=2*(hw+t_sh) y=-hh-t_sh h=2*(hh+t_sh) radius=hh+t_sh wall=t_sh\n';
    const half = await solve(tw), full = await solve(tw, { symmetry: false });
    agree('twinax, half vs full domain', half.m, full.m, 5e-3, 0.02);
    check('twinax: about 100 ohm differential', half.m[0].Z0 > 45 && half.m[0].Z0 < 55, `Z_odd ${half.m[0].Z0.toFixed(2)}`);
}

done();
