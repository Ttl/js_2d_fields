// Parallel-plate exact test — the loss path against CLOSED FORMS, not references.
//
// Geometry: wide trace (W/H = 50) centered between two ground planes, gap H above
// and below, homogeneous fill εr, on the TRIANGULAR backend whose open side walls
// are natural-Neumann (PMC) for the static solve. The field between the plates is
// uniform, so everything has a closed form:
//   C  = 2·ε0·εr·W/H            L_ext = μ0·H/(2W)         Z0 = sqrt(L/C)
//   C/C0 = εr (EXACT on any mesh — homogeneous-fill invariant)
//   G  = ω·C·tanδ
//   R(f) + jωL_int(f) from the exact slab impedance Z(d) = (1+j)/(σδ)·coth((1+j)d/δ):
//     trace splits at its mid-plane (H=0 there by symmetry) → Z(T/2);
//     each ground is a slab with a field-free back → Z(TG);
//     the two trace-face/ground pairs are in parallel → Z_int = [Z(T/2)+Z(TG)]/(2W).
//   The coth form is exact at EVERY frequency, so this pins the MQS volume
//   eddy-current loss solve through the whole skin-effect transition — the only
//   loss test in the suite with analytic truth (the rest compare to other solvers).
//
// Fringing at the trace edges (side clearance 3H) adds a real ~2.5% to C over the
// plate formula (upper-bounded by 6H/W = 12%, converged under mesh refinement), and
// Z0 shifts correspondingly — those comparisons get a ±4% band. The sharp
// assertions are the ones fringing cancels out of: C/C0, G/(ωC), and R(f)/coth.
// KNOWN LIMIT (asserted only loosely): near DC (t ≪ δ) the MQS R comes out ~16%
// below R_DC and L_int ~25% high — the near-DC internal-impedance estimation is a
// separate open issue; this test pins today's behavior from drifting further.
//
// Run: node tests/test_parallel_plate.js
import { MicrostripSolver } from '../src/microstrip.js';

const EPS0 = 8.8541878128e-12, MU0 = 4e-7 * Math.PI;
const H = 0.1e-3, W = 5e-3, T = 35e-6, TG = 35e-6;
const ER = 4.0, TAND = 0.02, SIG = 5.8e7;

let failures = 0;
function check(name, cond, detail = '') {
    console.log(`${cond ? '✓ PASS' : '✗ FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
    if (!cond) failures++;
}
const relErr = (a, b) => Math.abs(a - b) / Math.abs(b);
const pct = (a, b) => ((a - b) / b * 100).toFixed(2) + '%';

// Exact slab internal impedance (1+j)/(σδ)·coth((1+j)·d/δ) → [Re, Im]
function zSlab(d, delta) {
    const x = d / delta;
    const den = Math.cosh(2 * x) - Math.cos(2 * x);
    const cothRe = Math.sinh(2 * x) / den, cothIm = -Math.sin(2 * x) / den;
    const k = 1 / (SIG * delta);
    return [k * (cothRe - cothIm), k * (cothRe + cothIm)];
}
function exactRL(f) {
    const delta = 1 / Math.sqrt(Math.PI * f * MU0 * SIG);
    const [aR, aI] = zSlab(T / 2, delta), [bR, bI] = zSlab(TG, delta);
    return { R: (aR + bR) / (2 * W), Lint: (aI + bI) / (2 * W) / (2 * Math.PI * f) };
}

async function solveAt(f) {
    const s = new MicrostripSolver({
        trace_width: W, substrate_height: H, trace_thickness: T,
        gnd_thickness: TG, epsilon_r: ER, epsilon_r_top: ER,
        tan_delta: TAND, tan_delta_top: TAND, sigma_cond: SIG,
        enclosure_width: W + 6 * H, enclosure_height: H + T,
        freq: f, boundaries: ['open', 'open', 'gnd', 'gnd'],
    });
    s.mesh_backend = 'triangular';
    const log = console.log; console.log = () => {};
    try { return (await s.solve_adaptive()).modes[0]; }
    finally { console.log = log; }
}

const Cx = 2 * EPS0 * ER * W / H;
const Lx = MU0 * H / (2 * W);
const Z0x = Math.sqrt(Lx / Cx);

// ---- static / TEM quantities (frequency-independent; take the 1 GHz solve) ----
const m1 = await solveAt(1e9);
check('C/C0 == εr (homogeneous-fill static invariant, exact on any mesh)',
    relErr(m1.RLGC.C / m1.C0, ER) < 0.005, `C/C0=${(m1.RLGC.C / m1.C0).toFixed(5)}`);
check('C matches the plate formula (+fringing band)', relErr(m1.RLGC.C, Cx) < 0.04,
    `${(m1.RLGC.C * 1e12).toFixed(1)} vs ${(Cx * 1e12).toFixed(1)} pF/m, ${pct(m1.RLGC.C, Cx)}`);
check('Z0 matches sqrt(L/C) of the plates (fringing band)', relErr(m1.Z0, Z0x) < 0.04,
    `${m1.Z0.toFixed(4)} vs ${Z0x.toFixed(4)} Ω, ${pct(m1.Z0, Z0x)}`);
check('G/(ωC) == tanδ', relErr(m1.RLGC.G / (2 * Math.PI * 1e9 * m1.RLGC.C), TAND) < 0.01,
    `${(m1.RLGC.G / (2 * Math.PI * 1e9 * m1.RLGC.C)).toFixed(5)} vs ${TAND}`);

// ---- R(f), L(f) through the skin-effect transition vs the exact coth model ----
const CASES = [
    { f: 1e8, tolR: 0.03, tolL: 0.04, note: 'transition, t/δ≈2.7' },
    { f: 1e9, tolR: 0.03, tolL: 0.04, note: 'skin onset, t/δ≈8.4', mode: m1 },
    { f: 1e10, tolR: 0.05, tolL: 0.04, note: 'skin regime, t/δ≈27' },
];
const solved = [];
for (const c of CASES) {
    const m = c.mode ?? await solveAt(c.f);
    const ex = exactRL(c.f);
    solved.push({ f: c.f, R: m.RLGC.R, L: m.RLGC.L });
    check(`R at ${c.f / 1e9} GHz matches coth model (${c.note})`, relErr(m.RLGC.R, ex.R) < c.tolR,
        `${m.RLGC.R.toFixed(4)} vs ${ex.R.toFixed(4)} Ω/m, ${pct(m.RLGC.R, ex.R)}`);
    check(`L at ${c.f / 1e9} GHz matches plate + internal`, relErr(m.RLGC.L, Lx + ex.Lint) < c.tolL,
        `${(m.RLGC.L * 1e9).toFixed(3)} vs ${((Lx + ex.Lint) * 1e9).toFixed(3)} nH/m, ${pct(m.RLGC.L, Lx + ex.Lint)}`);
}
check('R increases with frequency, L decreases (skin effect direction)',
    solved[0].R < solved[1].R && solved[1].R < solved[2].R &&
    solved[0].L > solved[1].L && solved[1].L >= solved[2].L);

// ---- near-DC: known-imperfect regime, pinned loosely so it can't silently worsen ----
{
    const f = 1e6;
    const m = await solveAt(f);
    const ex = exactRL(f);
    console.log(`(info)  near-DC (1 MHz): R ${m.RLGC.R.toFixed(4)} vs exact ${ex.R.toFixed(4)} Ω/m (${pct(m.RLGC.R, ex.R)}), ` +
        `L ${(m.RLGC.L * 1e9).toFixed(2)} vs ${((Lx + ex.Lint) * 1e9).toFixed(2)} nH/m (${pct(m.RLGC.L, Lx + ex.Lint)}) — known MQS DC-limit bias`);
    check('near-DC R within the known ~±25% envelope', relErr(m.RLGC.R, ex.R) < 0.25);
    check('near-DC L within the known ~±35% envelope', relErr(m.RLGC.L, Lx + ex.Lint) < 0.35);
}

console.log(failures === 0 ? '\nPARALLEL PLATE OK' : `\nPARALLEL PLATE: ${failures} FAILURE(S)`);
process.exit(failures === 0 ? 0 : 1);
