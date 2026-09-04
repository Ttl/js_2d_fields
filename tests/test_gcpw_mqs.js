// GCPW / differential-GCPW MQS conductor loss (tri backend).
//
// Explicit ground rects (coplanar grounds, via slabs) used to force the
// perturbation path; MQS now treats them as passive C = 0 return conductors
// with a distance-graded skin band (fine at the slot edge, coarsening with
// distance from the signal, so the band cost is independent of how far the
// ground pour extends). These tests pin:
//   1. The 'auto' loss method takes the MQS path for GCPW (no applicability
//      warning) and produces R in a sane band around the perturbation estimate.
//   2. Physical limits: R(f=0) = geometric R_dc; wide-gap GCPW converges to
//      the microstrip R (the coplanar grounds stop mattering); R(f) increases.
//   3. The passive grounds actually carry return current: R_gnd > 0 with a
//      minority share of R_total.
//   4. The two-pass band respects its budget: trace budget + ground budget.
//   5. Differential GCPW solves both modes on the MQS path with sane R.
import { MicrostripSolver } from '../src/microstrip.js';
import { initTriBackend, TriBackend } from '../src/tri_solver/tri_backend.js';
import { mqsConductorLoss } from '../src/tri_solver/mqs_loss.js';

const ctx = await initTriBackend();
let failures = 0;
function check(name, ok, detail = '') {
    console.log(`${ok ? '✓ PASS' : '✗ FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
    if (!ok) failures++;
}
const relDiff = (a, b) => Math.abs(a - b) / Math.max(Math.abs(b), 1e-30);

const SIGMA = 5.8e7;
const base = {
    trace_width: 0.3e-3, substrate_height: 0.254e-3, trace_thickness: 35e-6,
    epsilon_r: 3.66, tan_delta: 0.003, sigma_cond: SIGMA, rq: 0,
    gnd_thickness: 35e-6, freq: 1e9,
};
const gcpw = { ...base, use_coplanar_gnd: true, gap: 0.15e-3, via_gap: 0.3e-3, use_vias: true };

async function tri(opts, triOpts = {}) {
    const s = new MicrostripSolver(opts);
    const b = new TriBackend(ctx, s, { maxNodes: 18000, ...triOpts });
    await b.buildMesh();
    return b;
}

// ---- 1+2+3+4: single-ended GCPW ----
{
    const b = await tri(gcpw);
    const bP = await tri(gcpw, { lossMethod: 'perturbation' });

    const r0 = b.solveAt(0).modes[0];
    let sig = 0, gnd = 0;
    for (const c of b.solver.conductors) {
        const a = Math.abs(c.width * c.height);
        if (c.is_signal) sig += a; else gnd += a;
    }
    const rdc = 1 / (SIGMA * sig) + 1 / (SIGMA * gnd);
    check('GCPW R(f=0) = geometric R_dc', relDiff(r0.RLGC.R, rdc) < 1e-9,
        `${r0.RLGC.R.toFixed(4)} vs ${rdc.toFixed(4)} Ω/m`);

    const rLow = b.solveAt(1e5).modes[0].RLGC.R;
    const r1 = b.solveAt(1e9).modes[0].RLGC.R;
    const r10 = b.solveAt(10e9).modes[0].RLGC.R;
    const warns = (b._modeWarnings || []).map(w => w.type);
    check('GCPW auto path is MQS (no applicability/fallback warning)',
        !warns.includes('mqs-shape') && !warns.includes('mqs-no-interior'), warns.join(',') || 'no warnings');
    check('GCPW R(f) increases with f', rLow < r1 && r1 < r10,
        `${rLow.toFixed(2)} < ${r1.toFixed(2)} < ${r10.toFixed(2)} Ω/m`);

    // MQS resolves the slot-corner current the perturbation SIBC over-counts, so
    // it should land BELOW the perturbation estimate at skin-regime frequencies,
    // but not far below (both model the same line).
    const p10 = bP.solveAt(10e9).modes[0].RLGC.R;
    check('GCPW MQS R within (-25%, +5%) of perturbation at 10 GHz',
        r10 > 0.75 * p10 && r10 < 1.05 * p10, `mqs ${r10.toFixed(2)} vs pert ${p10.toFixed(2)} Ω/m`);

    // Ground rects carry return current: direct MQS call on the cached skin mesh.
    const mqs = mqsConductorLoss(b._skinCache.mesh, b.condRect, 10e9, SIGMA,
        ctx.helpers.solveComplexSymmetric, 0, { cache: {} });
    const share = mqs.R_gnd / mqs.R_total;
    check('GCPW ground loss share in (5%, 60%)', share > 0.05 && share < 0.6,
        `R_trace ${mqs.R_trace.toFixed(2)}, R_gnd ${mqs.R_gnd.toFixed(2)} Ω/m → ${(100 * share).toFixed(1)}%`);

    // Two-pass band budget: trace cap + ground budget (default mqsMaxTris/2),
    // with slack for the conformity closure of the final refinement pass.
    const budget = 40000 + 20000;
    check('GCPW skin mesh respects trace+ground budget', b._skinCache.mesh.nTris <= 1.3 * budget,
        `${b._skinCache.mesh.nTris} tris vs budget ${budget}`);
}

// ---- 2b: wide-gap limit → microstrip ----
{
    const bM = await tri(base);
    const rM = bM.solveAt(10e9).modes[0].RLGC.R;
    const bW = await tri({ ...gcpw, gap: 0.6e-3, domain_width: 8e-3 });
    const rW = bW.solveAt(10e9).modes[0].RLGC.R;
    check('wide-gap GCPW R → microstrip R at 10 GHz (±5%)', relDiff(rW, rM) < 0.05,
        `gcpw ${rW.toFixed(2)} vs microstrip ${rM.toFixed(2)} Ω/m`);
}

// ---- 5: differential GCPW, both modes ----
{
    const b = await tri({ ...gcpw, trace_spacing: 0.2e-3 });
    const bP = await tri({ ...gcpw, trace_spacing: 0.2e-3 }, { lossMethod: 'perturbation' });
    const rA = b.solveAt(10e9), rP = bP.solveAt(10e9);
    for (let mi = 0; mi < rA.modes.length; mi++) {
        const m = rA.modes[mi], p = rP.modes[mi];
        const name = m.mode ?? `mode${mi}`;
        check(`diff GCPW ${name} R sane vs perturbation (±30%)`,
            m.RLGC.R > 0 && relDiff(m.RLGC.R, p.RLGC.R) < 0.30,
            `mqs ${m.RLGC.R.toFixed(2)} vs pert ${p.RLGC.R.toFixed(2)} Ω/m`);
        check(`diff GCPW ${name} Z0 matches perturbation backend (±1%)`,
            relDiff(m.Z0, p.Z0) < 0.01, `${m.Z0.toFixed(2)} vs ${p.Z0.toFixed(2)} Ω`);
    }
}

console.log(failures === 0 ? '\nALL GCPW MQS TESTS PASSED' : `\n${failures} TEST(S) FAILED`);
process.exit(failures === 0 ? 0 : 1);
