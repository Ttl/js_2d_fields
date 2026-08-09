// Coax surface roughness and plating — behavioural invariants, not closed forms.
//
// Split out of tests/test_coax.js (which owns the analytic closed forms and stays in the
// fast tier) because each case is a full mesh build + solve. Lives in the slow tier
// alongside the other plating suites.
//
// The point of these checks is that a circle has ONE continuous surface: the per-face
// top/sides/bottom plating model has nothing to select between, so CoaxSolver normalizes
// any enabled face to all-around and the surface-impedance classifiers resolve it through
// a dedicated 'all' face. Roughness needs no coax-specific handling at all — bare metal
// already flows through the single uniform Zbare group.
//
// Run: node tests/test_coax_plating.js
import { CoaxSolver } from '../src/coax.js';

const D_IN = 0.92e-3, D_OUT = 2.95e-3, ER = 2.1, TAND = 2e-4, SIG = 5.8e7;
const a = D_IN / 2, b = D_OUT / 2;

let failures = 0;
function check(name, cond, detail = '') {
    console.log(`${cond ? '\u2713 PASS' : '\u2717 FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
    if (!cond) failures++;
}
const pct = (a_, b_) => ((a_ - b_) / b_ * 100).toFixed(3) + '%';
function near(name, got, want, tolPct) {
    check(name, Math.abs(got / want - 1) * 100 <= tolPct,
        `${got.toExponential(5)} vs ${want.toExponential(5)}, ${pct(got, want)}, tol ${tolPct}%`);
}

function makeSolver(extra = {}) {
    const s = new CoaxSolver({
        inner_diameter: D_IN, dielectric_diameter: D_OUT, epsilon_r: ER,
        tan_delta: TAND, sigma_cond: SIG, freq: 1e9, ...extra,
    });
    s.use_causal_materials = false;
    return s;
}
// NOTE: this file runs five full mesh-build + solve cycles (~50 s), which is why it is
// in the slow tier rather than beside the closed-form checks. Lowering max_nodes to cut
// that is counterproductive — a budget below the natural mesh size makes the
// triangle-budget retry loop in buildMesh re-run gmsh several times, which measured
// SLOWER (68 s) than just letting the default mesh through.
async function solveCoax(extra = {}, opts = {}) {
    const s = makeSolver(extra);
    const log = console.log; console.log = () => {};
    try {
        const first = await s.solve_adaptive({ max_nodes: 20000, ...opts });
        return { s, first };
    } finally { console.log = log; }
}
async function at(s, f) {
    const log = console.log; console.log = () => {};
    try { return (await s.computeAtFrequency(f)).modes[0]; } finally { console.log = log; }
}

const { s } = await solveCoax();

// ---------------------------------------------------------------- roughness
console.log('\n=== surface roughness ===');
{
    // Both conductors are bare, so roughness reaches the whole surface through the
    // single uniform Zbare group — no per-face bookkeeping is involved on a circle.
    // The model is the Grujic gradient model (see surface_roughness.js), which is not
    // capped at the Hammerstad 2x, so the assertion is monotonicity in Rq rather than a
    // fixed bound.
    const smooth = await at(s, 10e9);
    const r05 = await at((await solveCoax({ rq: 0.5e-6 })).s, 10e9);
    const r10 = await at((await solveCoax({ rq: 1e-6 })).s, 10e9);
    check('alpha_c increases monotonically with Rq',
        smooth.alpha_c < r05.alpha_c && r05.alpha_c < r10.alpha_c,
        `${smooth.alpha_c.toFixed(4)} < ${r05.alpha_c.toFixed(4)} < ${r10.alpha_c.toFixed(4)} dB/m`);
    near('roughness barely moves Z0', r10.Z0, smooth.Z0, 0.1);
    near('roughness does not move alpha_d', r10.alpha_d, smooth.alpha_d, 0.1);
}

// ---------------------------------------------------------------- plating
console.log('\n=== plating ===');
{
    // A circle has one continuous surface, so CoaxSolver normalizes any enabled face to
    // all-around; the shield stays bare.
    const plated = makeSolver({ plating: { sigma: 6.3e7, thickness: 4e-6, rq: 0, top: true, sides: false, bottom: false } });
    check('plating normalized to all-around',
        plated.plating.all === true && plated.plating.top && plated.plating.sides && plated.plating.bottom);
    check('plating has no thick-corner wrap (no corners)', plated.plating.thick_corners === false);
    check('plating applies to the inner conductor only',
        plated.conductors[0].plating !== null && plated.conductors[1].plating === null);

    // Silver (6.3e7) over copper (5.8e7) must lower loss; tin (8.7e6) must raise it.
    // Both are bounded by the loss of a solid conductor of the plating metal, since the
    // plating is several skin depths thick at 10 GHz (delta_Ag = 0.63 um << 4 um).
    const ag = await solveCoax({ plating: { sigma: 6.3e7, thickness: 4e-6, rq: 0, top: true, sides: true, bottom: true } });
    const sn = await solveCoax({ plating: { sigma: 8.7e6, thickness: 4e-6, rq: 0, top: true, sides: true, bottom: true } });
    const bare = await at(s, 10e9);
    const mAg = await at(ag.s, 10e9), mSn = await at(sn.s, 10e9);
    check('silver plating lowers alpha_c', mAg.alpha_c < bare.alpha_c,
        `${bare.alpha_c.toFixed(4)} -> ${mAg.alpha_c.toFixed(4)} dB/m`);
    check('tin plating raises alpha_c', mSn.alpha_c > bare.alpha_c,
        `${bare.alpha_c.toFixed(4)} -> ${mSn.alpha_c.toFixed(4)} dB/m`);
    // Only the inner conductor is plated, and it carries 1/a of the 1/a+1/b loss, so the
    // shift is bounded by that share of the pure-metal ratio.
    const share = (1 / a) / (1 / a + 1 / b);
    const agLimit = 1 - share * (1 - Math.sqrt(SIG / 6.3e7));
    check('silver shift is bounded by the inner conductor share',
        mAg.alpha_c / bare.alpha_c > agLimit - 1e-3,
        `ratio ${(mAg.alpha_c / bare.alpha_c).toFixed(4)}, limit ${agLimit.toFixed(4)}`);
    near('plating barely moves Z0', mAg.Z0, bare.Z0, 0.1);
}


console.log(`\n${failures === 0 ? 'ALL TESTS PASSED' : failures + ' TEST(S) FAILED'}`);
process.exit(failures === 0 ? 0 : 1);
