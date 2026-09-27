// Exact DC line parameters of rectangular conductors (src/dc_inductance.js):
//   1. the closed-form rectangle log integral against brute-force quadrature, and the
//      self geometric mean distance of a square (0.44705 a)
//   2. the internal inductance of a square conductor and of a thin strip with a far
//      return on a wall ground, against mu0/2pi ln(r_pec / GMD): 0.59017 a over
//      0.44705 a for the square, (w/4) over 0.22313 w for the strip; and the slab
//      term mu0 d/3 int K^2 of a thick wall under a microstrip
//   3. the boundary-element grading is converged at the default alpha
//   4. the DC resistance: the traces in series with the grounds, which drop out when
//      a wall ground is there
import { rectLogIntegral, dcLineParameters } from '../src/dc_inductance.js';
import { check, relErr as rel, done } from './helpers.js';

const MU = 2e-7;

// 1. Closed form.
function brute(a, b, n) {
    let s = 0;
    const ha = [(a.x1 - a.x0) / n, (a.y1 - a.y0) / n], hb = [(b.x1 - b.x0) / n, (b.y1 - b.y0) / n];
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
        const x = a.x0 + (i + 0.5) * ha[0], y = a.y0 + (j + 0.5) * ha[1];
        for (let k = 0; k < n; k++) for (let l = 0; l < n; l++) {
            const u = x - (b.x0 + (k + 0.5) * hb[0]), v = y - (b.y0 + (l + 0.5) * hb[1]);
            s += 0.5 * Math.log(u * u + v * v);
        }
    }
    return s * ha[0] * ha[1] * hb[0] * hb[1];
}
const A = { x0: 0, x1: 1, y0: 0, y1: 0.3 }, B = { x0: 2, x1: 2.5, y0: -1, y1: 1 }, C = { x0: 0.2, x1: 1.7, y0: 0.5, y1: 0.6 };
for (const [name, a, b] of [['separate', A, B], ['overlapping in x', A, C], ['diagonal', B, C]]) {
    const exact = rectLogIntegral(a, b), q = brute(a, b, 40);
    check(`rect log integral, ${name}: closed form matches quadrature`, Math.abs(exact - q) < 2e-4 * Math.max(1, Math.abs(q)),
        `${exact.toFixed(6)} vs ${q.toFixed(6)}`);
}
const unit = { x0: 0, x1: 1, y0: 0, y1: 1 };
check('square self integral is ln of its GMD 0.44705', Math.abs(rectLogIntegral(unit, unit) - Math.log(0.447049)) < 1e-5,
    `${rectLogIntegral(unit, unit).toFixed(6)}`);

// 2. Internal inductance with a far return on a wall ground, thin enough that its own
// slab term is negligible.
const wallAndTrace = (trace, d, wallW) => [
    { x_min: -wallW / 2, x_max: wallW / 2, y_min: -d - 1e-6, y_max: -d, is_signal: false },
    { ...trace, is_signal: true, polarity: 1 }];
const walls = new Set([0]);
{
    const a = 1e-3, ref = MU * Math.log(0.59017 / 0.44705);
    const r = dcLineParameters(wallAndTrace({ x_min: -a / 2, x_max: a / 2, y_min: 0, y_max: a }, 0.1, 2), 'single', { unlimited: walls });
    check('square conductor: DC internal inductance', rel(r.Lint, ref) < 0.005, `${(r.Lint * 1e9).toFixed(2)} vs ${(ref * 1e9).toFixed(2)} nH/m`);
}
{
    const ref = MU * Math.log(0.25 / 0.22313);
    const r = dcLineParameters(wallAndTrace({ x_min: -0.5e-3, x_max: 0.5e-3, y_min: 0, y_max: 1e-6 }, 0.1, 2), 'single', { unlimited: walls });
    check('thin strip (t/w = 1e-3): DC internal inductance near the zero-thickness limit', rel(r.Lint, ref) < 0.03,
        `${(r.Lint * 1e9).toFixed(2)} vs ${(ref * 1e9).toFixed(2)} nH/m`);
}

// The wall slab term grows with the wall thickness as mu0 d/3 int K^2: four times the
// thickness adds four times the term of a thin wall.
{
    const trace = { x_min: -0.175e-3, x_max: 0.175e-3, y_min: 0.21e-3, y_max: 0.245e-3, is_signal: true, polarity: 1 };
    const Lint = d => dcLineParameters([{ x_min: -3e-3, x_max: 3e-3, y_min: -d, y_max: 0, is_signal: false }, trace],
        'single', { unlimited: walls }).Lint;
    const thin = Lint(1e-7), t1 = Lint(35e-6), t4 = Lint(140e-6);
    check('wall slab term proportional to the wall thickness', rel(t4 - thin, 4 * (t1 - thin)) < 0.01,
        `${((t1 - thin) * 1e9).toFixed(2)} / ${((t4 - thin) * 1e9).toFixed(2)} nH/m at 35 / 140 um`);
}

// 3. Grading: a microstrip, with its plane as a wall ground and as a finite ground.
const ms = [{ x_min: -3e-3, x_max: 3e-3, y_min: -35e-6, y_max: 0, is_signal: false },
            { x_min: -0.175e-3, x_max: 0.175e-3, y_min: 0.21e-3, y_max: 0.245e-3, is_signal: true, polarity: 1 }];
for (const [name, w] of [['wall ground', walls], ['finite ground', new Set()]]) {
    const coarse = dcLineParameters(ms, 'single', { unlimited: w }), fine = dcLineParameters(ms, 'single', { unlimited: w, alpha: 0.1 });
    check(`microstrip, ${name}: L and L_pec converged at the default grading`,
        rel(coarse.L, fine.L) < 2e-4 && rel(coarse.Lpec, fine.Lpec) < 2e-4,
        `L ${(coarse.L * 1e9).toFixed(3)} / ${(fine.L * 1e9).toFixed(3)}, L_pec ${(coarse.Lpec * 1e9).toFixed(3)} / ${(fine.Lpec * 1e9).toFixed(3)} nH/m`);
}

// 4. DC resistance.
{
    const sigma = 5.8e7, Rt = 1 / (sigma * 0.35e-3 * 35e-6), Rg = 1 / (sigma * 6e-3 * 35e-6);
    const finite = dcLineParameters(ms, 'single'), wall = dcLineParameters(ms, 'single', { unlimited: walls });
    check('microstrip on a finite ground: trace and ground in series', rel(finite.R, Rt + Rg) < 1e-12, `${finite.R.toFixed(6)} ohm/m`);
    check('microstrip on a wall ground: the trace alone', rel(wall.R, Rt) < 1e-12, `${wall.R.toFixed(6)} ohm/m`);
    const pair = [ms[0], { ...ms[1], x_min: 0.1e-3, x_max: 0.3e-3 }, { ...ms[1], x_min: -0.3e-3, x_max: -0.1e-3, polarity: -1 }];
    const Rp = 1 / (sigma * 0.2e-3 * 35e-6);
    const odd = dcLineParameters(pair, 'odd'), even = dcLineParameters(pair, 'even');
    check('differential pair, odd mode: the ground carries no DC current', rel(odd.R, Rp) < 1e-12, `${odd.R.toFixed(6)} ohm/m`);
    check('differential pair, even mode: the ground returns 2I', rel(even.R, Rp + 2 * Rg) < 1e-12, `${even.R.toFixed(6)} ohm/m`);
}

done();

console.log('\nALL CHECKS PASSED');
