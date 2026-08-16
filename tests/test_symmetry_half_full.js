// Half-domain ≡ full-domain identity across feature-stacked families (tri backend).
//
// The symmetry-plane void bug produced ~23%-wrong C with certificates claiming
// ~1%: Richardson certifies convergence of the meshing process it is given, not
// correctness of the limit. The cheap structural guard is identity testing —
// solve the same geometry on the half domain (default for mirror-symmetric
// geometries) and the full domain (tri_opts.symmetry: false) and demand the same
// answer. test_symmetry_plane_mask.js pins the original void-bug geometry; this
// file extends the identity to the other families and features that exercise the
// symmetry machinery differently:
//   - single-ended: the plane cuts the TRACE itself (MQS symmetry-plane
//     edge/corner exclusion, half-domain integral scaling)
//   - GCPW: passive pours + via slabs as C = 0 return conductors
//   - stripline cover, closed enclosure walls, ground cutout centered on the plane
//
// Loss methods are matched between the two sides (auto takes MQS only WITH
// symmetry, so a naive R comparison would compare MQS against perturbation):
//   - single-ended: half {auto → MQS} vs full {lossMethod 'mqs'} — MQS identity
//   - differential: half {lossMethod 'perturbation'} vs full {auto →
//     perturbation} — perturbation identity (full-domain diff MQS is refused,
//     see below)
// Also pins the mqs-asym-pair refusal: forcing 'mqs' on a full-domain
// differential pair must warn and fall back, not return the silently ~4×-low
// parallel-drive answer it used to.
import { MicrostripSolver } from '../src/microstrip.js';

let failures = 0;
function check(name, ok, detail = '') {
    console.log(`${ok ? '✓ PASS' : '✗ FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
    if (!ok) failures++;
}
const relDiff = (a, b) => Math.abs(a - b) / Math.max(Math.abs(a), Math.abs(b), 1e-30);
const quiet = async (fn) => {
    const log = console.log, warn = console.warn;
    console.log = () => {}; console.warn = () => {};
    try { return await fn(); } finally { console.log = log; console.warn = warn; }
};

const FAMILIES = [
    {
        name: 'microstrip + sm + plating', diff: false, cTol: 0.01, rTol: 0.03,
        geom: {
            trace_width: 0.35e-3, substrate_height: 0.21e-3, trace_thickness: 35e-6, gnd_thickness: 35e-6,
            epsilon_r: 4.4, tan_delta: 0.02, sigma_cond: 5.8e7, freq: 1e9, rq: 0,
            boundaries: ['open', 'open', 'open', 'gnd'],
            use_sm: true, sm_t_sub: 20e-6, sm_t_trace: 20e-6, sm_t_side: 20e-6, sm_er: 3.5, sm_tand: 0.02,
            plating: { sigma: 1e7, thickness: 4e-6, rq: 0, top: true, sides: true, bottom: false, thick_corners: true },
        },
    },
    {
        name: 'gcpw + vias + sm', diff: false, cTol: 0.01, rTol: 0.03,
        geom: {
            trace_width: 0.3e-3, substrate_height: 0.254e-3, trace_thickness: 35e-6, gnd_thickness: 35e-6,
            epsilon_r: 3.66, tan_delta: 0.003, sigma_cond: 5.8e7, freq: 1e9, rq: 0,
            boundaries: ['open', 'open', 'open', 'gnd'],
            use_coplanar_gnd: true, gap: 0.15e-3, via_gap: 0.3e-3, use_vias: true,
            use_sm: true, sm_t_sub: 15e-6, sm_t_trace: 15e-6, sm_t_side: 15e-6, sm_er: 3.4, sm_tand: 0.02,
        },
    },
    {
        name: 'diff stripline', diff: true, cTol: 0.01, rTol: 0.08,
        geom: {
            trace_width: 0.15e-3, trace_spacing: 0.15e-3, substrate_height: 0.2e-3, trace_thickness: 17e-6,
            gnd_thickness: 17e-6, epsilon_r: 4.1, epsilon_r_top: 4.1, tan_delta: 0.02, tan_delta_top: 0.02,
            enclosure_height: 0.4e-3, sigma_cond: 5.8e7, freq: 2e9, rq: 0,
            boundaries: ['open', 'open', 'gnd', 'gnd'],
        },
    },
    {
        name: 'diff microstrip + enclosure', diff: true, cTol: 0.01, rTol: 0.08,
        geom: {
            trace_width: 0.2e-3, trace_spacing: 0.2e-3, substrate_height: 0.2e-3, trace_thickness: 35e-6,
            gnd_thickness: 20e-6, epsilon_r: 4.4, tan_delta: 0.02, sigma_cond: 5.8e7, freq: 1e9, rq: 0,
            enclosure_width: 2.4e-3, enclosure_height: 0.8e-3,
            boundaries: ['gnd', 'gnd', 'gnd', 'gnd'],
        },
    },
    {
        name: 'diff microstrip + gnd cutout', diff: true, cTol: 0.03, rTol: 0.08,
        geom: {
            trace_width: 0.2e-3, trace_spacing: 0.15e-3, substrate_height: 0.25e-3, trace_thickness: 35e-6,
            gnd_thickness: 20e-6, epsilon_r: 4.4, tan_delta: 0.02, sigma_cond: 5.8e7, freq: 1e9, rq: 0,
            boundaries: ['open', 'open', 'open', 'gnd'],
            gnd_cut_width: 0.25e-3, gnd_cut_sub_h: 0.1e-3,
        },
    },
];

async function solve(geom, triOpts) {
    const s = new MicrostripSolver({ ...geom, mesh_backend: 'triangular', nx: 30, ny: 30 });
    s.tri_opts = triOpts;
    const r = await quiet(() => s.solve_adaptive({
        max_iters: 12, energy_tol: 0.01, param_tol: 0.05, max_nodes: 20000, min_converged_passes: 2,
    }));
    return { modes: r.modes, warns: (r.warnings || []).map(w => w.type) };
}

for (const fam of FAMILIES) {
    if (!fam.diff) {
        // 2 solves: C/Z0 and MQS-R identities from the same pair
        const half = await solve(fam.geom, {});                                 // auto → MQS on half domain
        const full = await solve(fam.geom, { symmetry: false, lossMethod: 'mqs' });
        const h = half.modes[0], f = full.modes[0];
        check(`${fam.name}: half C == full C (±${100 * fam.cTol}%)`,
            relDiff(h.RLGC.C, f.RLGC.C) < fam.cTol,
            `${(h.RLGC.C * 1e12).toFixed(2)} vs ${(f.RLGC.C * 1e12).toFixed(2)} pF/m`);
        check(`${fam.name}: half Z0 == full Z0 (±${100 * fam.cTol}%)`,
            relDiff(h.Z0, f.Z0) < fam.cTol, `${h.Z0.toFixed(2)} vs ${f.Z0.toFixed(2)} Ω`);
        check(`${fam.name}: half MQS R == full MQS R (±${100 * fam.rTol}%)`,
            relDiff(h.RLGC.R, f.RLGC.R) < fam.rTol,
            `${h.RLGC.R.toFixed(2)} vs ${f.RLGC.R.toFixed(2)} Ω/m`);
    } else {
        // 3 solves: half/full auto for C/Z0; half-perturbation vs full-auto for R
        const half = await solve(fam.geom, {});
        const full = await solve(fam.geom, { symmetry: false });                // auto → perturbation on full
        const halfP = await solve(fam.geom, { lossMethod: 'perturbation' });
        for (let i = 0; i < half.modes.length; i++) {
            const mode = half.modes[i].mode;
            check(`${fam.name} ${mode}: half C == full C (±${100 * fam.cTol}%)`,
                relDiff(half.modes[i].RLGC.C, full.modes[i].RLGC.C) < fam.cTol,
                `${(half.modes[i].RLGC.C * 1e12).toFixed(2)} vs ${(full.modes[i].RLGC.C * 1e12).toFixed(2)} pF/m`);
            check(`${fam.name} ${mode}: half pert R == full pert R (±${100 * fam.rTol}%)`,
                relDiff(halfP.modes[i].RLGC.R, full.modes[i].RLGC.R) < fam.rTol,
                `${halfP.modes[i].RLGC.R.toFixed(2)} vs ${full.modes[i].RLGC.R.toFixed(2)} Ω/m`);
        }
    }
}

// Forced 'mqs' on a full-domain differential pair: must warn and fall back to
// perturbation (it used to silently return both modes ~4× low).
{
    const fam = FAMILIES[2];   // diff stripline
    const full = await solve(fam.geom, { symmetry: false });                    // perturbation reference
    const fullM = await solve(fam.geom, { symmetry: false, lossMethod: 'mqs' });
    check('full-domain diff + forced mqs warns (mqs-asym-pair)',
        fullM.warns.includes('mqs-asym-pair'), fullM.warns.join(',') || 'no warnings');
    check('full-domain diff + forced mqs falls back to perturbation R',
        relDiff(fullM.modes[0].RLGC.R, full.modes[0].RLGC.R) < 0.01,
        `${fullM.modes[0].RLGC.R.toFixed(2)} vs ${full.modes[0].RLGC.R.toFixed(2)} Ω/m`);
}

console.log(failures === 0 ? '\nALL HALF≡FULL IDENTITY TESTS PASSED' : `\n${failures} TEST(S) FAILED`);
process.exit(failures === 0 ? 0 : 1);
