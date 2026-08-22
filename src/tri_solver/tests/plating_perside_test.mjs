// Per-face plating: verify the full-wave backend applies plating per selected
// face (top/sides/bottom) and tracks the FDM backend across configurations.
import { MicrostripSolver } from '../../microstrip.js';
import { initTriBackend, TriBackend } from '../tri_backend.js';

const ctx = await initTriBackend();
const F = 20e9;

// thin tin-like plating (lower sigma than copper) on a copper trace
const base = {
    trace_width: 0.3e-3, substrate_height: 0.2e-3, trace_thickness: 35e-6,
    epsilon_r: 3.66, tan_delta: 0.003, sigma_cond: 5.8e7, rq: 0,
};
const plating = (faces) => ({
    sigma: 8.7e6, thickness: 4e-6, rq: 0.3e-6,
    top: !!faces.top, sides: !!faces.sides, bottom: !!faces.bottom,
    thick_corners: false,
});

// The FDM side needs a LOT of nodes here. Its conductor loss is a surface integral
// of the vacuum field, which converges from below and much more slowly than the
// field energy the adaptive loop actually gates on: at 20 GHz on this geometry it
// runs ~12-16% low at 40k nodes (what this test used to use), ~6% low at 137k, and
// ~3% low at 280k. The original wide 10-12% bands were reading that mesh error, not
// a model difference, so this now solves to ~280k and the bands reflect real
// cross-backend spread.
async function fdm(opts) {
    const s = new MicrostripSolver({ ...opts, freq: F });
    s.ensure_mesh();
    const cached = await s.solve_adaptive({ max_iters: 15, energy_tol: 0.001, param_tol: 0.02,
        max_nodes: 250000, min_converged_passes: 2, certify: false });
    return (await s.computeAtFrequency(F, cached)).modes[0];
}
async function tri(opts) {
    const s = new MicrostripSolver({ ...opts, freq: F });
    const b = new TriBackend(ctx, s);
    await b.buildMesh();
    return b.solveAt(F).modes[0];
}
const pct = (a, b) => (100 * (a - b) / b).toFixed(1) + '%';

const cases = [
    ['no plating',          {}],
    ['all faces',           { plating: plating({ top: 1, sides: 1, bottom: 1 }) }],
    ['top only',            { plating: plating({ top: 1, sides: 0, bottom: 0 }) }],
    ['top+sides (no bot)',  { plating: plating({ top: 1, sides: 1, bottom: 0 }) }],
    ['bottom only',         { plating: plating({ top: 0, sides: 0, bottom: 1 }) }],
];

console.log(`microstrip @ ${F/1e9} GHz — base Cu σ=5.8e7, plating σ=8.7e6 (tin), t=4µm, Rq=0.3µm\n`);
console.log('config                 tri α_c    fdm α_c    Δ        tri R    fdm R');
console.log('-'.repeat(78));
const results = {};
for (const [name, opts] of cases) {
    const o = { ...base, ...opts };
    const t = await tri(o);
    const f = await fdm(o);
    results[name] = { t, f };
    console.log(
        `${name.padEnd(22)} ${t.alpha_c.toFixed(3).padStart(7)}    ${f.alpha_c.toFixed(3).padStart(7)}` +
        `    ${pct(t.alpha_c, f.alpha_c).padStart(6)}   ${t.RLGC.R.toFixed(1).padStart(6)}   ${f.RLGC.R.toFixed(1).padStart(6)}`
    );
}

// Per-face loss ADDED by plating each face (vs bare copper) — the physical
// signature of per-face plating. Bottom (facing ground) carries most current and
// adds the most loss; top the least. (Default solver = auto → MQS: each face's
// smooth current is weighted by its own plating impedance. The per-face split
// differs from FDM by method spread — MQS is a volume eddy solve, FDM a surface-
// impedance grid — so agreement is a band, not tight; ~4% on uniform plating.)
const aT = n => results[n].t.alpha_c, aF = n => results[n].f.alpha_c;
const tBare = aT('no plating'), fBare = aF('no plating');
console.log('\nloss added per face (α_c − bare):  tri      fdm');
console.log(`  top    ${(aT('top only')-tBare).toFixed(2).padStart(7)}  ${(aF('top only')-fBare).toFixed(2).padStart(7)}`);
console.log(`  bottom ${(aT('bottom only')-tBare).toFixed(2).padStart(7)}  ${(aF('bottom only')-fBare).toFixed(2).padStart(7)}`);
console.log(`  all    ${(aT('all faces')-tBare).toFixed(2).padStart(7)}  ${(aF('all faces')-fBare).toFixed(2).padStart(7)}`);

// Assertions: robust properties + a method-spread band vs FDM (FDM is itself a
// model, not ground truth; MQS's volume trace/ground split is arguably more
// accurate, so we assert a 20% band rather than tight agreement on plated cases).
const checks = [
    // Bare baseline: measured ~2.9% at this budget (the residual is the FDM's
    // remaining surface under-resolution, not a model difference); band leaves ~2x.
    ['no-plating tri ≈ FDM (<7%)',        Math.abs(aT('no plating')/aF('no plating') - 1) < 0.07],
    // AGGREGATE plating effect (all faces) — should track FDM closely once the
    // trace body uses bulk σ and every face is correctly classified+scaled
    ['all-faces tri ≈ FDM (<7%)',         Math.abs(aT('all faces')/aF('all faces') - 1) < 0.07],
    // per-face selectivity actually does something and is correctly ordered
    ['plating raises loss (all > none)',  aT('all faces') > aT('no plating')],
    ['top-only between none and all',     tBare < aT('top only') && aT('top only') < aT('all faces')],
    ['bottom adds more than top',         (aT('bottom only')-tBare) > (aT('top only')-tBare)],
    ['top+sides ≥ top only',              aT('top+sides (no bot)') >= aT('top only') - 1e-9],
    ['top+sides < all faces',             aT('top+sides (no bot)') < aT('all faces')],
    // individual per-face configs: the bottom/top/side split near the singular
    // bottom corners is method-dependent (full-wave vs FDM grid), so a wider band
    ['per-face configs within 10% of FDM',
        ['top only','top+sides (no bot)','bottom only']
            .every(n => Math.abs(aT(n)/aF(n) - 1) < 0.10)],
];
console.log('\nchecks:');
let ok = true;
for (const [label, pass] of checks) { console.log(`  ${pass ? '✓' : '✗'} ${label}`); ok = ok && pass; }
console.log(`\n${ok ? '✓ ALL CHECKS PASSED' : '✗ SOME CHECKS FAILED'}`);
process.exit(ok ? 0 : 1);
