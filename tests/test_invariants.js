// Identities the solver must satisfy EXACTLY, with no reference model anywhere.
//
// Every other physics test in the suite compares against something: a closed form for one
// specific geometry (coax, waveguide, parallel plate), HFSS, a measurement, or the other
// backend. This file compares the solver against ITSELF, by choosing configurations where
// the answer is known a priori — so a failure localizes to a formula rather than to "the
// reference disagrees by 3%, is that the mesh?".
//
// Four families, in the order they run (the pure-arithmetic ones first, so a broken
// S-parameter path fails in milliseconds instead of after a minute of solving):
//
//   D. S-PARAMETERS — a lossless line is unitary, a line referenced to its own Zc has
//      S11 = 0, a zero-length line is the identity, two half-length sections cascade to
//      one full-length section, and a SYMMETRIC pair pushed through the 4-port MTL path
//      must reproduce the modal odd/even formulas with exactly zero mode conversion.
//      All of these hold to machine precision — they are algebra, not discretization.
//
//   A. HOMOGENEOUS FILL — with one permittivity everywhere the mode is exactly TEM at
//      every frequency, so eps_eff_mode == er on ANY mesh with ZERO dispersion, and
//      alpha_d = pi*f*sqrt(er)*tand/c0 with the geometry cancelled out entirely. This is
//      the only place in the suite where eps_eff and alpha_d have an exact expected value
//      on a rectangular geometry, which makes it the suite's absolute-accuracy meter for
//      the full-wave eigensolve. Measured today: the quasi-static path is exact to 5e-6%,
//      the eigensolve carries a +0.63% bias on the default mesh that converges away as
//      O(h) (0.34% / 0.17% at 986 / 2932 triangles).
//
//   B. DECOUPLING LIMIT — a pair spaced 20 substrate heights apart IS two independent
//      single-ended lines, so both modes must land on the single-ended solve. Validates
//      the whole differential path (mode decomposition, per-trace power normalization,
//      symmetry handling) against a path with no coupling in it at all.
//
//   E. MTL MATRICES — for a homogeneous fill [L][C] = mu0*eps0*er*I exactly, including
//      the OFF-diagonal, which is a cancellation between (Lo-Le) and (Co-Ce) that no
//      single-mode check can see. Plus the structural facts: [C] and [L] symmetric,
//      mutual capacitance negative, diagonally dominant, [L] positive definite.
//
// Two things must be switched off or accounted for or the identities fail for legitimate
// reasons, and both cost real debugging time to rediscover:
//
//   • Djordjevic-Sarkar moves er away from its nominal value, so the stripline sections
//     (er = 4.3, tand = 0.02) are wrong by the DS shift unless use_causal_materials is
//     false: measured -2.05% on eps and +0.72% on alpha_d at 5 GHz, on both backends,
//     which would fail the gates below by 2x. The er = 1 microstrip is NOT affected —
//     both backends skip the model when |er - 1| < 1e-6 or tand == 0 (djordjevic_sarkar.js
//     applyDjordjevicSarkar, tri_backend.js _applyCausal), and that case trips both
//     guards, so causal on/off is bit-identical there. The flag is set uniformly anyway,
//     so the er = 1 geometry stays safe if its materials are ever edited.
//   • The reported eps_eff is c^2*L*C with the INTERNAL inductance included, so it sits
//     above er by L_int/L_ext (~1.8% at 1 GHz here) even for a perfect solver. The
//     quantity that must equal er is eps_eff_mode (= c^2*L_external*C). The relation
//     between the two is itself asserted below, which ties this file to the surface
//     reactance work in tests/test_surface_reactance.js.
//
// Runtime ~60 s, dominated by seven mesh builds — hence the slow tier. Section D alone is
// instant and could be split into the fast tier if that is ever wanted.
//
// Run: node tests/test_invariants.js
import { MicrostripSolver } from '../src/microstrip.js';
import { computeSParamsSingleEnded, computeSParamsDifferential,
         computeSParamsDifferentialMTL, buildPhysicalRLGC } from '../src/sparameters.js';

const EPS0 = 8.854187817e-12, MU0 = 4e-7 * Math.PI;
const C0 = 1 / Math.sqrt(EPS0 * MU0);
const NP_TO_DB = 8.685889638;
const SIG = 5.8e7;

let failures = 0;
function check(name, cond, detail = '') {
    console.log(`${cond ? '✓ PASS' : '✗ FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
    if (!cond) failures++;
}
const pct = (got, want) => ((got - want) / want * 100).toFixed(4) + '%';
function near(name, got, want, tolPct) {
    check(name, Math.abs(got / want - 1) * 100 <= tolPct,
        `${got.toExponential(6)} vs ${want.toExponential(6)}, ${pct(got, want)}, tol ${tolPct}%`);
}
const quiet = async (fn) => {
    const log = console.log, warn = console.warn;
    console.log = () => {}; console.warn = () => {};
    try { return await fn(); } finally { console.log = log; console.warn = warn; }
};

// One solve, returning the whole result so differential callers can reach both modes and
// the physical matrix. `hints` overrides the triangular mesher's element sizes, which is
// how the convergence sweep in section A refines without touching the geometry.
async function solve(opts, freq, { backend = 'triangular', hints = null } = {}) {
    const s = new MicrostripSolver({ ...opts, mesh_backend: backend });
    s.use_causal_materials = false;         // DS would move er off its nominal value
    if (hints) s.tri_mesh_hints = hints;
    const cached = await quiet(() => s.solve_adaptive(backend === 'triangular'
        ? { max_nodes: 200000 }
        : { max_iters: 6, energy_tol: 0.005, param_tol: 0.02, max_nodes: 40000 }));
    const result = await quiet(() => s.computeAtFrequency(freq, cached));
    return { s, result, mode: result.modes[0],
             nTris: s._triBackend?.mesh?.nTris ?? s.triMesh?.nTris ?? null };
}
// Several frequencies off ONE mesh — the homogeneous invariants are frequency statements,
// so they must be read on a fixed discretization or the mesh moves under them.
async function sweep(opts, freqs, { backend = 'triangular', hints = null } = {}) {
    const s = new MicrostripSolver({ ...opts, mesh_backend: backend });
    s.use_causal_materials = false;
    if (hints) s.tri_mesh_hints = hints;
    const cached = await quiet(() => s.solve_adaptive(backend === 'triangular'
        ? { max_nodes: 200000 }
        : { max_iters: 6, energy_tol: 0.005, param_tol: 0.02, max_nodes: 40000 }));
    const out = new Map();
    for (const f of freqs) {
        out.set(f, (await quiet(() => s.computeAtFrequency(f, cached))).modes[0]);
    }
    return { s, at: out, nTris: s._triBackend?.mesh?.nTris ?? s.triMesh?.nTris ?? null };
}

// ============================================================================ D
// S-parameters. Pure arithmetic on an RLGC set — no solver, no mesh, so every tolerance
// here is machine epsilon and any failure is a formula error.
console.log('\n=== D. S-parameter identities ===');
{
    const f = 10e9;
    const L = 3.5e-7, C = 1.2e-10;
    const lossless = { R: 0, L, G: 0, C };
    const Zc = Math.sqrt(L / C);
    const mag = (z) => Math.hypot(z.re, z.im);
    const p = (z) => z.re * z.re + z.im * z.im;

    for (const len of [1e-3, 25e-3, 137e-3]) {
        // A lossless line can only redistribute power between reflection and
        // transmission. True at ANY reference impedance — 50 Ω is deliberately
        // mismatched here so S11 is large and the identity is not trivially satisfied.
        const s50 = computeSParamsSingleEnded(f, lossless, len, 50);
        near(`lossless ${(len * 1e3).toFixed(0)}mm: |S11|^2+|S21|^2 == 1 at Zref=50`,
            p(s50.S11) + p(s50.S21), 1, 1e-9);
        // Referenced to its own characteristic impedance a uniform line cannot reflect,
        // at any length — this catches a sign or sinh/cosh slip that unitarity alone
        // would absorb.
        const sZc = computeSParamsSingleEnded(f, lossless, len, Zc);
        check(`lossless ${(len * 1e3).toFixed(0)}mm: S11 == 0 at Zref=Zc`,
            mag(sZc.S11) < 1e-12, `|S11| = ${mag(sZc.S11).toExponential(2)}`);
        // Reciprocity: a passive isotropic line has S12 = S21 identically.
        check(`lossless ${(len * 1e3).toFixed(0)}mm: S12 == S21 (reciprocity)`,
            mag({ re: s50.S12.re - s50.S21.re, im: s50.S12.im - s50.S21.im }) < 1e-12);
    }

    // Zero length is the identity two-port, whatever the reference impedance.
    const z = computeSParamsSingleEnded(f, lossless, 0, 50);
    check('zero-length line is the identity two-port',
        mag(z.S11) < 1e-12 && Math.abs(mag(z.S21) - 1) < 1e-12,
        `|S11| = ${mag(z.S11).toExponential(2)}, |S21| = ${mag(z.S21).toFixed(12)}`);

    // A lossy line must LOSE power — the same expression that is exactly 1 above has to
    // come out below 1, so the unitarity check cannot be passing for a degenerate reason.
    const lossy = { R: 80, L, G: 0.05, C };
    const sl = computeSParamsSingleEnded(f, lossy, 25e-3, 50);
    const sum = p(sl.S11) + p(sl.S21);
    check('lossy line is strictly sub-unitary', sum < 1 && sum > 0.5,
        `|S11|^2+|S21|^2 = ${sum.toFixed(6)}`);

    // Cascade: the line is uniform, so two half-length sections chained through their
    // T-matrices must reproduce the full-length section exactly. This is the sharpest
    // available check that gamma and Zc enter the two-port consistently — a wrong Zc
    // reflects at the artificial junction and the cascade stops matching.
    const cx = (re, im) => ({ re, im });
    const cmul = (a, b) => cx(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
    const cdiv = (a, b) => { const d = b.re * b.re + b.im * b.im;
        return cx((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d); };
    const csub = (a, b) => cx(a.re - b.re, a.im - b.im);
    const cneg = (a) => cx(-a.re, -a.im);
    const sToT = (S) => {
        const det = csub(cmul(S.S11, S.S22), cmul(S.S12, S.S21));
        return { T11: cdiv(cneg(det), S.S21), T12: cdiv(S.S11, S.S21),
                 T21: cdiv(cneg(S.S22), S.S21), T22: cdiv(cx(1, 0), S.S21) };
    };
    const tMul = (A, B) => {
        const e = (a, b, c, d) => cx(cmul(a, b).re + cmul(c, d).re, cmul(a, b).im + cmul(c, d).im);
        return { T11: e(A.T11, B.T11, A.T12, B.T21), T12: e(A.T11, B.T12, A.T12, B.T22),
                 T21: e(A.T21, B.T11, A.T22, B.T21), T22: e(A.T21, B.T12, A.T22, B.T22) };
    };
    const tToS = (T) => ({ S11: cdiv(T.T12, T.T22), S21: cdiv(cx(1, 0), T.T22),
        S12: csub(T.T11, cdiv(cmul(T.T12, T.T21), T.T22)), S22: cneg(cdiv(T.T21, T.T22)) });
    const full = computeSParamsSingleEnded(f, lossy, 40e-3, 50);
    const half = computeSParamsSingleEnded(f, lossy, 20e-3, 50);
    const casc = tToS(tMul(sToT(half), sToT(half)));
    check('two half-length sections cascade to the full-length section',
        mag(csub(casc.S11, full.S11)) < 1e-12 && mag(csub(casc.S21, full.S21)) < 1e-12,
        `dS11 ${mag(csub(casc.S11, full.S11)).toExponential(2)},`
        + ` dS21 ${mag(csub(casc.S21, full.S21)).toExponential(2)}`);

    // The 4-port MTL path against the modal path, on a SYMMETRIC pair where both are
    // valid. computeSParamsDifferential builds the mixed-mode terms from the odd/even
    // two-ports and hardcodes the conversion terms to zero; computeSParamsDifferentialMTL
    // builds a genuine 4x4 from the physical [R][L][G][C] and TRANSFORMS to mixed mode,
    // so its SDC/SCD are computed, not assumed. On a symmetric pair the two must agree
    // and the conversion must vanish — this is the only check in the suite that the
    // asymmetric-capable path degenerates correctly onto the symmetric one.
    const odd = { R: 90, L: 3.2e-7, G: 0.06, C: 1.5e-10 };
    const even = { R: 70, L: 3.8e-7, G: 0.04, C: 1.1e-10 };
    const P = buildPhysicalRLGC(odd, even, null);
    const mtl = computeSParamsDifferentialMTL(f, P.R, P.L, P.G, P.C, 40e-3, 50);
    const modal = computeSParamsDifferential(f, odd, even, 40e-3, 50);
    const dd = (a, b) => Math.hypot(a.re - b.re, a.im - b.im);
    for (const k of ['SDD11', 'SDD21', 'SCC11', 'SCC21']) {
        check(`symmetric pair: MTL 4-port reproduces the modal ${k}`,
            dd(mtl[k], modal[k]) < 1e-12, `delta ${dd(mtl[k], modal[k]).toExponential(2)}`);
    }
    for (const k of ['SDC11', 'SDC21', 'SCD11', 'SCD21']) {
        check(`symmetric pair: ${k} == 0 (no mode conversion by symmetry)`,
            mtl[k].abs() < 1e-12, `|${k}| = ${mtl[k].abs().toExponential(2)}`);
    }
    // ...and the pair is genuinely coupled, so the zeros above are not zero-input zeros.
    check('the test pair is actually coupled (|SDD21| differs from |SCC21|)',
        Math.abs(mtl.SDD21.abs() - mtl.SCC21.abs()) > 1e-3,
        `|SDD21| ${mtl.SDD21.abs().toFixed(6)} vs |SCC21| ${mtl.SCC21.abs().toFixed(6)}`);
}

// ============================================================================ A
// Homogeneous fill. er = 1 everywhere makes the whole cross-section vacuum, so the exact
// answer for eps_eff_mode is 1.000000000 at every frequency on every mesh — no closed
// form for the geometry is needed, and none of it depends on the trace shape.
console.log('\n=== A. homogeneous fill: microstrip with er = 1 (vacuum cross-section) ===');
const MS_VAC = {
    trace_width: 0.3e-3, substrate_height: 0.2e-3, trace_thickness: 35e-6,
    epsilon_r: 1, tan_delta: 0, sigma_cond: SIG, rq: 0, freq: 1e9,
    boundaries: ['open', 'open', 'open', 'gnd'],
};
// F_STATIC_MAX sits between these two groups: below it the reported eps_eff_mode IS the
// quasi-static W_eps/W_air ratio, above it the full-wave eigenvalue. They are separate
// code paths that must land on the same exact answer, so both are checked.
const F_STATIC = [1e6, 5e7], F_WAVE = [1e9, 40e9];
{
    const { at, nTris } = await sweep(MS_VAC, [...F_STATIC, ...F_WAVE]);
    console.log(`  default mesh: ${nTris} triangles`);

    // The quasi-static path: a uniform epsMap makes the eps-weighted and vacuum stiffness
    // matrices differ by a scalar, so the ratio is exactly er on any mesh. Same statement
    // tests/test_coax.js makes for a shaped conductor; this is the rectangular case, and
    // it doubles as a material-tagging check (a stray air sliver would drift the ratio).
    for (const f of F_STATIC) {
        near(`quasi-static @ ${(f / 1e6).toFixed(0)} MHz: eps_eff_mode is EXACTLY 1`,
            at.get(f).eps_eff_mode, 1, 1e-4);
    }

    // The full-wave path: a homogeneous guide has NO dispersion, beta = k0*sqrt(er) at
    // every frequency, so the two eigensolves must agree with each other far more tightly
    // than either agrees with 1. Splitting the statement this way separates "the
    // eigensolve is biased" (next check) from "the eigensolve drifts with frequency".
    const wave = F_WAVE.map((f) => at.get(f).eps_eff_mode);
    near('full-wave: eps_eff_mode at 40 GHz == at 1 GHz (homogeneous fill cannot disperse)',
        wave[1], wave[0], 0.02);

    // The absolute bias of the eigensolve on the default mesh. Measured 0.632% high, and
    // it is discretization, not a formula error — the convergence sweep below shows it
    // going away as O(h). The gate exists to catch it GROWING (a coarser default mesh, a
    // worse element, a broken size field), which nothing else in the suite would see.
    for (const f of F_WAVE) {
        near(`full-wave @ ${(f / 1e9).toFixed(0)} GHz: eps_eff_mode == 1 on the default mesh`,
            at.get(f).eps_eff_mode, 1, 1.0);
    }

    // Algebraic ties between the reported quantities, at machine precision. These are
    // what make the checks above meaningful for the OTHER outputs: if eps_eff_mode is
    // right then so are Z0 and eps_eff, because they are the same three numbers.
    const m = at.get(1e9);
    near('eps_eff_mode == c^2*L_external*C', m.eps_eff_mode, C0 * C0 * m.L_external * m.RLGC.C, 1e-4);
    near('Z0 == sqrt(eps_eff_mode)/(c*C)', m.Z0, Math.sqrt(m.eps_eff_mode) / (C0 * m.RLGC.C), 1e-4);
    // The reported eps_eff carries the internal inductance and eps_eff_mode does not, so
    // their ratio is exactly L/L_external. This is where a regression in the surface
    // reactance (tests/test_surface_reactance.js) would surface as an eps_eff error.
    near('eps_eff / eps_eff_mode == L / L_external (the internal-inductance share)',
        m.eps_eff / m.eps_eff_mode, (m.L_external + m.L_internal) / m.L_external, 1e-4);
    check('the internal-inductance share is large enough for that ratio to mean something',
        m.eps_eff / m.eps_eff_mode - 1 > 1e-3,
        `eps_eff ${m.eps_eff.toFixed(6)} vs eps_eff_mode ${m.eps_eff_mode.toFixed(6)}`);
}

// Mesh convergence of the eigensolve bias. The exact answer is known, so this measures
// the ORDER of the method rather than just its value: halving the element size must halve
// the error. Measured 0.337% at 986 triangles and 0.172% at 2932 — first order, ratio
// 1.95. A method that stopped converging (or converged to the wrong constant) would keep
// the ratio near 1 while every reference-based test in the suite still passed.
console.log('\n--- A. eigensolve convergence against the exact homogeneous answer ---');
{
    const errs = [];
    for (const [hFine, hCoarse] of [[52e-6, 260e-6], [26e-6, 130e-6]]) {
        const { at, nTris } = await sweep(MS_VAC, [1e9], { hints: { hFine, hCoarse } });
        const err = Math.abs(at.get(1e9).eps_eff_mode - 1);
        errs.push(err);
        console.log(`  hFine=${(hFine * 1e6).toFixed(0)}um  ${nTris} triangles  `
            + `eps_eff_mode ${at.get(1e9).eps_eff_mode.toFixed(7)}  error ${(err * 100).toFixed(4)}%`);
    }
    check('refining halves the error (first-order convergence to the exact answer)',
        errs[0] / errs[1] > 1.6 && errs[0] / errs[1] < 2.6,
        `error ratio ${(errs[0] / errs[1]).toFixed(3)} for a 2x element-size reduction`);
    check('the refined mesh is within 0.25% of the exact answer',
        errs[1] < 2.5e-3, `${(errs[1] * 100).toFixed(4)}%`);
}

// The same idea with loss, on a homogeneous STRIPLINE (er = er_top, fully enclosed): the
// dielectric attenuation of a TEM line is alpha_d = pi*f*sqrt(er)*tand/c0, in which the
// geometry has cancelled out completely — alpha_d = omega*tand*C*Z0/2 and C*Z0 =
// sqrt(er)/c0. So this is an absolute check on the dielectric loss integral with no mesh
// term in the expected value at all. Measured: triangular 0.30%, rectilinear 1.85%.
console.log('\n=== A. homogeneous fill: stripline, geometry-free alpha_d ===');
const ER_SL = 4.3, TAND_SL = 0.02;
const SL_HOMO = {
    trace_width: 0.2e-3, substrate_height: 0.2e-3, trace_thickness: 35e-6,
    epsilon_r: ER_SL, epsilon_r_top: ER_SL, tan_delta: TAND_SL, tan_delta_top: TAND_SL,
    enclosure_height: 0.2e-3 + 35e-6, sigma_cond: SIG, rq: 0, freq: 5e9,
    boundaries: ['open', 'open', 'gnd', 'gnd'],
};
const alphaD_x = (f) => Math.PI * f * Math.sqrt(ER_SL) / C0 * TAND_SL * NP_TO_DB;
{
    for (const [backend, tolAlpha, tolZ0] of [['triangular', 0.6, 0.8], ['rectilinear', 2.5, 1e-4]]) {
        const { mode: m } = await solve(SL_HOMO, 5e9, { backend });
        near(`${backend}: alpha_d == pi*f*sqrt(er)*tand/c0 (geometry cancels)`,
            m.alpha_d, alphaD_x(5e9), tolAlpha);
        // Z0 = sqrt(er)/(c*C) for a homogeneous line. On the rectilinear backend
        // L_external is defined as 1/(c^2*C_vacuum), which makes this exact by
        // construction; on the triangular backend L_external comes from its own solve, so
        // there it is a real check that the two agree.
        near(`${backend}: Z0 == sqrt(er)/(c*C)`, m.Z0, Math.sqrt(ER_SL) / (C0 * m.RLGC.C), tolZ0);
        if (m.eps_eff_mode !== undefined) {
            near(`${backend}: eps_eff_mode == er`, m.eps_eff_mode, ER_SL, 1.0);
        }
    }
}

// ============================================================================ B
// A pair spaced far enough apart is two single-ended lines. Nothing about the differential
// machinery — modal decomposition, per-trace power normalization, the symmetry gate — has
// a single-ended counterpart to be checked against anywhere else in the suite; this makes
// one. Measured at 20 substrate heights: Z 0.27% / 0.12%, eps_eff 0.59% / 0.11%,
// alpha_c 0.29% / 0.10%.
console.log('\n=== B. decoupling limit: a widely spaced pair == two single-ended lines ===');
{
    const SE = {
        trace_width: 0.2e-3, substrate_height: 0.2e-3, trace_thickness: 35e-6,
        epsilon_r: 3.66, tan_delta: 0.003, sigma_cond: SIG, rq: 0, freq: 5e9,
        boundaries: ['open', 'open', 'open', 'gnd'],
    };
    const SPACING = 20 * SE.substrate_height;
    const { mode: se } = await solve(SE, 5e9);
    const { result } = await solve({ ...SE, trace_spacing: SPACING }, 5e9);
    const odd = result.modes.find((m) => m.mode === 'odd');
    const even = result.modes.find((m) => m.mode === 'even');
    check('the spaced solve produced both modes', !!odd && !!even);

    for (const [tag, m] of [['odd', odd], ['even', even]]) {
        near(`${tag} mode Z0 -> single-ended Z0`, m.Z0, se.Z0, 1.0);
        near(`${tag} mode eps_eff -> single-ended eps_eff`, m.eps_eff, se.eps_eff, 1.2);
        near(`${tag} mode alpha_c -> single-ended alpha_c`, m.alpha_c, se.alpha_c, 1.0);
        near(`${tag} mode alpha_d -> single-ended alpha_d`, m.alpha_d, se.alpha_d, 1.0);
    }
    // The residual coupling must still be there, with the right sign: the odd mode packs
    // more capacitance between the traces (lower Z), the even mode less (higher Z). If
    // this ordering were absent the checks above could be passing because the pair solve
    // silently collapsed to a single trace.
    check('residual coupling has the right sign: Z_odd < Z_single-ended < Z_even',
        odd.Z0 < se.Z0 && se.Z0 < even.Z0,
        `${odd.Z0.toFixed(3)} < ${se.Z0.toFixed(3)} < ${even.Z0.toFixed(3)} ohm`);
    const spread = (even.Z0 - odd.Z0) / se.Z0;
    check('the modal spread is small but nonzero at 20h', spread > 1e-4 && spread < 0.01,
        `(Z_even - Z_odd)/Z_se = ${(spread * 100).toFixed(3)}%`);
}

// ============================================================================ E
// The multiconductor identity. For a homogeneous fill the per-unit-length matrices satisfy
// [L][C] = mu0*eps0*er*I EXACTLY — diagonal AND off-diagonal. The off-diagonal is the
// interesting half: it is a cancellation between the inductive and capacitive coupling
// (Lo-Le against Co-Ce) that no single-mode check can reach, and it is the sharpest test
// available of the modal-to-physical transform that feeds the exported 4-port.
console.log('\n=== E. MTL matrices: [L][C] == mu0*eps0*er*I for a homogeneous fill ===');
{
    const { result } = await solve({ ...SL_HOMO, trace_spacing: 0.2e-3 }, 5e9);
    const odd = result.modes.find((m) => m.mode === 'odd');
    const even = result.modes.find((m) => m.mode === 'even');
    check('the homogeneous pair produced both modes', !!odd && !!even);

    // Per mode first, so a matrix failure below can be attributed. Both modes of a
    // homogeneous line travel at the same speed: c^2*L_ext*C = er for each.
    for (const [tag, m] of [['odd', odd], ['even', even]]) {
        near(`${tag} mode: c^2*L_external*C == er`,
            C0 * C0 * m.L_external * m.RLGC.C, ER_SL, 1.2);
        // alpha_d has no geometry in it, so BOTH modes must land on the same closed form
        // as the single-ended stripline above — the same number from three solves.
        near(`${tag} mode: alpha_d == pi*f*sqrt(er)*tand/c0`, m.alpha_d, alphaD_x(5e9), 1.0);
    }

    // The physical 2x2 matrices exactly as the S-parameter export builds them. L uses
    // L_external: the internal inductance is metal, not field in the dielectric, and it
    // is not part of the TEM identity.
    const P = buildPhysicalRLGC(
        { ...odd.RLGC, L: odd.L_external }, { ...even.RLGC, L: even.L_external },
        result.physMatrix ?? null);
    const { L, C } = P;
    const mul2 = (A, B) => [[A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1]],
                            [A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]]];
    const LC = mul2(L, C);
    const want = ER_SL / (C0 * C0);
    console.log(`  [C] = [[${C[0][0].toExponential(4)}, ${C[0][1].toExponential(4)}], ...]`
        + `  [L] = [[${L[0][0].toExponential(4)}, ${L[0][1].toExponential(4)}], ...]`);

    near('[L][C] diagonal == er/c^2', LC[0][0], want, 1.2);
    near('[L][C] is diagonally uniform (both traces identical)', LC[1][1], LC[0][0], 1e-6);
    // The cancellation. Measured 8.6e-4 of the diagonal, an order of magnitude below the
    // ~0.6% error carried by either mode on its own — which is the point: the off-diagonal
    // is right even where the individual modes are still converging.
    check('[L][C] off-diagonal vanishes (inductive and capacitive coupling cancel)',
        Math.abs(LC[0][1]) / LC[0][0] < 5e-3,
        `|off|/diag = ${(Math.abs(LC[0][1]) / LC[0][0]).toExponential(2)}`);
    check('[L][C] off-diagonal is symmetric', Math.abs(LC[0][1] - LC[1][0]) <= 1e-12 * Math.abs(LC[0][0]));

    // Structural facts about the Maxwell matrices, independent of the fill.
    check('[C] is symmetric', Math.abs(C[0][1] - C[1][0]) <= 1e-12 * Math.abs(C[0][0]));
    check('[L] is symmetric', Math.abs(L[0][1] - L[1][0]) <= 1e-12 * Math.abs(L[0][0]));
    check('mutual capacitance is negative (Maxwell convention)', C[0][1] < 0,
        `C12 = ${C[0][1].toExponential(4)} F/m`);
    check('[C] is diagonally dominant', Math.abs(C[0][0]) > Math.abs(C[0][1]),
        `|C11| ${Math.abs(C[0][0]).toExponential(4)} vs |C12| ${Math.abs(C[0][1]).toExponential(4)}`);
    check('mutual inductance is positive', L[0][1] > 0, `L12 = ${L[0][1].toExponential(4)} H/m`);
    check('[L] is positive definite', L[0][0] > 0 && L[0][0] * L[1][1] - L[0][1] * L[1][0] > 0,
        `det = ${(L[0][0] * L[1][1] - L[0][1] * L[1][0]).toExponential(4)}`);
    check('the pair is genuinely coupled (off-diagonals are not negligible)',
        Math.abs(C[0][1]) / Math.abs(C[0][0]) > 0.01 && Math.abs(L[0][1]) / L[0][0] > 0.01,
        `|C12/C11| = ${(Math.abs(C[0][1]) / Math.abs(C[0][0])).toFixed(4)},`
        + ` L12/L11 = ${(L[0][1] / L[0][0]).toFixed(4)}`);

    // Modal round-trip: the physical matrices must decompose back to the modal values the
    // solver reported. Pins the sign convention (odd = X11 - X12, even = X11 + X12) that
    // the exported 4-port depends on.
    near('round-trip: C11 - C12 == C_odd', C[0][0] - C[0][1], odd.RLGC.C, 1e-6);
    near('round-trip: C11 + C12 == C_even', C[0][0] + C[0][1], even.RLGC.C, 1e-6);
    near('round-trip: L11 - L12 == L_odd', L[0][0] - L[0][1], odd.L_external, 1e-6);
    near('round-trip: L11 + L12 == L_even', L[0][0] + L[0][1], even.L_external, 1e-6);
}

console.log(`\n${failures === 0 ? 'ALL TESTS PASSED' : failures + ' CHECK(S) FAILED'}`);
process.exit(failures === 0 ? 0 : 1);
