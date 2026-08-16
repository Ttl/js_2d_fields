// QS (rectilinear) wide-gap GCPW → microstrip continuity.
//
// A GCPW whose ground gap is many trace widths IS a microstrip — the tri backend
// already pins this on its MQS path (test_gcpw_mqs.js section 2b, ±5%). The QS
// backend must agree with itself the same way, INCLUDING in the DC→skin
// transition: the transition-notch calibration in calculate_conductor_loss used
// to be gated on !use_coplanar_gnd, which put a ~7% line-type step between a
// microstrip and the same line built as a wide-gap GCPW at δ/t ≈ 0.3–0.4.
//
// The residual gap under a uniform formula is ~7%: the QS surface integrand
// still collects loss on the (far, current-free per MQS) coplanar pour and via
// faces — a known wide-gap GCPW bias, tracked in docs/backend_agreement_handover.md.
// The ±10% band catches the gate regression (13.6% pre-fix) without pinning that
// separate bias tighter than it deserves.
import { MicrostripSolver } from '../src/microstrip.js';

let failures = 0;
function check(name, ok, detail = '') {
    console.log(`${ok ? '✓ PASS' : '✗ FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
    if (!ok) failures++;
}

const MU0 = 4e-7 * Math.PI;
const SIGMA = 5.8e7, T = 35e-6, W = 0.3e-3;
// f chosen so δ/t = 0.3 — the peak of the old gate's discontinuity.
const DELTA = 0.3 * T;
const FREQ = 2 / (2 * Math.PI * MU0 * SIGMA * DELTA * DELTA);

const base = {
    trace_width: W, substrate_height: 0.254e-3, trace_thickness: T,
    epsilon_r: 3.66, tan_delta: 0.003, sigma_cond: SIGMA, rq: 0,
    gnd_thickness: T, freq: FREQ, nx: 30, ny: 30,
    boundaries: ['open', 'open', 'open', 'gnd'],
};

async function solveR(opts) {
    const s = new MicrostripSolver(opts);
    const log = console.log;
    console.log = () => {};
    try {
        const r = await s.solve_adaptive({
            max_iters: 10, energy_tol: 0.01, param_tol: 0.05, max_nodes: 20000, min_converged_passes: 2,
        });
        return r.modes[0].RLGC.R;
    } finally {
        console.log = log;
    }
}

const rMs = await solveR(base);
const rW = await solveR({ ...base, use_coplanar_gnd: true, gap: 10 * W, via_gap: 0.3e-3, use_vias: true });
const rel = Math.abs(rW - rMs) / rMs;
check('wide-gap GCPW R → microstrip R at δ/t = 0.3 (±10%)', rel < 0.10,
    `gcpw ${rW.toFixed(2)} vs microstrip ${rMs.toFixed(2)} Ω/m, ${(100 * rel).toFixed(1)}%`);

console.log(failures === 0 ? '\nALL QS GCPW CONTINUITY TESTS PASSED' : `\n${failures} TEST(S) FAILED`);
process.exit(failures === 0 ? 0 : 1);
