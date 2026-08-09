// The internal inductance is the surface REACTANCE, not the surface resistance.
//
// A conductor's surface impedance has two halves and the transmission line uses both:
//
//     R_ac    = ∮ Re(Zs)|K|²dl / (2P)        L_int = ∮ Im(Zs)|K|²dl / (2P·ω)
//
// Smooth metal in the skin regime has Zs = Rs(1 + j), so the two integrals coincide and
// L_int = R_ac/ω. That closed form is a TRAP: it is exact for bare metal and silently
// wrong for every rough or plated surface, where Zs sits off 45°. At 1 µm rms copper
// Im(Zs)/Re(Zs) reaches 4 at 1 GHz and 13 at 10 GHz, so rebuilding L_int from R_ac
// recovers a fifth of the true increment and the line comes out fast — ~1% in β at 2 µm
// rms, measured against scikit-rf's Coaxial in test_dev/test_coax_vs_skrf.py section 3c.
// tri_backend.js accumulated exactly that closed form until X_ac was added beside R_ac.
//
// These checks need no reference model. With ONE surface impedance over the whole
// boundary, the geometric integral is a constant that CANCELS in the ratio between two
// solves of the same geometry (plating and roughness change no geometry, so both land
// on the identical mesh):
//
//     L_int(rough) / L_int(smooth)  ==  Im(Zs_rough) / Im(Zs_smooth)      exactly
//
// A solver deriving L_int from R_ac lands on Re(Zs_rough)/Re(Zs_smooth) instead, which
// is 3-5x smaller here — the two candidate answers are far apart, so the identity can be
// gated tightly and still discriminate. Each case asserts the wrong answer too, so the
// separation is visible in the output rather than assumed.
//
// tri_backend.js builds L_int at THREE independent sites, and each one had to learn
// this separately — so each gets a case here:
//
//   perturbation  (X_ac beside R_ac)         geometries with explicit ground rects:
//                                            coax, GCPW, broadside stripline
//   MQS           (X_total from mqs_loss.js) symmetric, wall-absorbed ground: microstrip
//   waveguide TE  (X_ac from the wall Zs)    rectangular waveguide
//
// The rectilinear (FDM) backend has always done this correctly — field_solver.js
// accumulates Im(Zs) into its own L_int — so it needs no case of its own; what it does
// provide is a second opinion, and microstrip below is chosen to run the MQS path so a
// regression there shows up as a backend disagreement too.
//
// Run: node tests/test_surface_reactance.js
import { CoaxSolver } from '../src/coax.js';
import { MicrostripSolver } from '../src/microstrip.js';
import { RectWaveguideSolver } from '../src/rect_waveguide.js';
import { calculate_Zrough, calculate_Zrough_layered } from '../src/surface_roughness.js';

const SIG = 5.8e7;
const D_IN = 0.92e-3, D_OUT = 2.95e-3;
const a = D_IN / 2, b = D_OUT / 2;
const RQ = 1e-6;
// Tin: 6.7x worse than the copper bulk and thinner than its own skin depth here, so the
// layered impedance is a genuine two-metal mix rather than either metal on its own.
const TIN = { sigma: 8.7e6, thickness: 0.5e-6, rq: 0 };
const FREQS = [1e9, 1e10];
// WR-90 cuts off at 6.56 GHz; below that the guide reports NaN for everything but alpha.
const WG_FREQS = [1e10, 1.2e10];

let failures = 0;
function check(name, cond, detail = '') {
    console.log(`${cond ? '✓ PASS' : '✗ FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
    if (!cond) failures++;
}
const pct = (got, want) => ((got - want) / want * 100).toFixed(3) + '%';
function near(name, got, want, tolPct) {
    check(name, Math.abs(got / want - 1) * 100 <= tolPct,
        `${got.toExponential(5)} vs ${want.toExponential(5)}, ${pct(got, want)}, tol ${tolPct}%`);
}

const quiet = async (fn) => {
    const log = console.log, warn = console.warn;
    console.log = () => {}; console.warn = () => {};
    try { return await fn(); } finally { console.log = log; console.warn = warn; }
};

// L_internal at each frequency of interest, keyed by frequency.
async function internalL(solver, freqs = FREQS) {
    const cached = await quiet(() => solver.solve_adaptive({ max_nodes: 20000 }));
    const out = new Map();
    for (const f of freqs) {
        const m = await quiet(() => solver.computeAtFrequency(f, cached).then((r) => r.modes[0]));
        out.set(f, m.L_internal);
    }
    return out;
}

const coax = (extra = {}) => {
    const s = new CoaxSolver({
        inner_diameter: D_IN, dielectric_diameter: D_OUT, epsilon_r: 2.1,
        tan_delta: 2e-4, sigma_cond: SIG, freq: 1e9, rq: 0, ...extra,
    });
    s.use_causal_materials = false;
    return s;
};
// GCPW as app_solver.js builds it: MicrostripSolver + use_coplanar_gnd. Its coplanar
// grounds are explicit ground rects, which is what routes it to the perturbation path.
const gcpw = (extra = {}) => new MicrostripSolver({
    substrate_height: 200e-6, trace_width: 300e-6, trace_thickness: 35e-6,
    gnd_thickness: 35e-6, epsilon_r: 3.5, tan_delta: 0.005, sigma_cond: SIG,
    freq: 1e9, boundaries: ['open', 'open', 'open', 'gnd'], rq: 0,
    use_coplanar_gnd: true, gap: 150e-6, via_gap: null, use_vias: true,
    mesh_backend: 'triangular', ...extra,
});
// The same solver WITHOUT coplanar grounds: symmetric, ground absorbed into the wall,
// so there is no ground rect and MQS applies. That is the only difference from gcpw().
const microstrip = (extra = {}) => new MicrostripSolver({
    substrate_height: 200e-6, trace_width: 400e-6, trace_thickness: 35e-6,
    gnd_thickness: 35e-6, epsilon_r: 4.3, tan_delta: 0.02, sigma_cond: SIG,
    freq: 1e9, nx: 200, ny: 100, boundaries: ['open', 'open', 'open', 'gnd'], rq: 0,
    mesh_backend: 'triangular', ...extra,
});
const waveguide = (extra = {}) => new RectWaveguideSolver({
    width: 22.86e-3, height: 10.16e-3, epsilon_r: 1, tan_delta: 0,
    sigma_cond: SIG, freq: 1e10, rq: 0, ...extra,
});

// One ratio check: L_int moved by exactly the ratio of the surface REACTANCES, and
// demonstrably NOT by the ratio of the surface resistances.
function checkRatio(label, bareL, caseL, zsOf, freqs = FREQS) {
    for (const f of freqs) {
        const ghz = (f / 1e9).toFixed(0);
        const zBare = calculate_Zrough(f, SIG, 0), zCase = zsOf(f);
        const got = caseL.get(f) / bareL.get(f);
        const wantIm = zCase.im / zBare.im, wantRe = zCase.re / zBare.re;
        // The mesh is shared, so the geometric integral cancels outright — this is an
        // identity, not an approximation, and 0.5% is loose for it by three decades.
        near(`${label} @ ${ghz} GHz: L_int ratio == Im(Zs) ratio`, got, wantIm, 0.5);
        // ...and the resistance ratio is somewhere else entirely. Without this the gate
        // above could be satisfied by a coincidence at one frequency.
        check(`${label} @ ${ghz} GHz: Re(Zs) ratio is a clearly different answer`,
            Math.abs(got / wantRe - 1) > 0.25,
            `L_int ${got.toFixed(4)}, Im(Zs) ${wantIm.toFixed(4)}, Re(Zs) ${wantRe.toFixed(4)}`
            + ` — R_ac/omega would be off by ${pct(wantRe, wantIm)}`);
    }
}

// ------------------------------------------------------------------- coax
console.log('\n=== coax: one continuous surface per conductor ===');
const coaxBare = await internalL(coax());

// Absolute anchor first: for SMOOTH metal both conventions agree, so this pins that
// L_int is the surface reactance in absolute terms and not merely proportional to it.
// The residual is the mesh discretization the loss integral already carries (~1%, the
// same figure R lands on in tests/test_coax.js), not a model difference.
for (const f of FREQS) {
    const w = 2 * Math.PI * f;
    const closed = calculate_Zrough(f, SIG, 0).im * (1 / a + 1 / b) / (2 * Math.PI);
    near(`smooth @ ${(f / 1e9).toFixed(0)} GHz: omega*L_int == Im(Zs)*(1/a + 1/b)/(2pi)`,
        coaxBare.get(f) * w, closed, 2.0);
}

// Roughness: both conductors bare metal, so ONE impedance covers the whole boundary.
checkRatio('coax rough 1um', coaxBare, await internalL(coax({ rq: RQ })),
    (f) => calculate_Zrough(f, SIG, RQ));

// Plating, on BOTH conductors so the boundary again carries a single impedance — this
// time the layered plating-over-bulk one, a different code path to the same identity.
checkRatio('coax tin-plated', coaxBare,
    await internalL(coax({ plating: { ...TIN, inner: true, outer: true } })),
    (f) => calculate_Zrough_layered(f, SIG, TIN.rq, TIN.sigma, TIN.thickness));

// ------------------------------------------------------------------- gcpw
// The ratio is a property of Zs alone, so a geometry with rectangular faces and coplanar
// grounds must produce the SAME number as the coax above. That it does is the evidence
// that the reactance is carried through the per-face surface groups, not special-cased.
console.log('\n=== gcpw: rectangular faces, coplanar grounds (perturbation path) ===');
checkRatio('gcpw rough 1um', await internalL(gcpw()), await internalL(gcpw({ rq: RQ })),
    (f) => calculate_Zrough(f, SIG, RQ));

// ------------------------------------------------------------------- mqs
// A centred microstrip over a wall-absorbed ground is symmetric with no ground rect, so
// it takes the MQS volume eddy solve instead — a different assembly of L_int (X_total
// out of mqs_loss.js) reaching the same identity. The R(f)/X(f) anchor caches inside
// that path interpolate BOTH, so a mismatch there would surface here as well.
console.log('\n=== microstrip: MQS path ===');
checkRatio('microstrip rough 1um', await internalL(microstrip()),
    await internalL(microstrip({ rq: RQ })), (f) => calculate_Zrough(f, SIG, RQ));

// ------------------------------------------------------------------- waveguide
// The TE mode has no conductors at all — the loss is on the enclosing wall, and its
// equivalent circuit builds L_int from the same wall surface impedance. Swept above the
// WR-90 cutoff (6.56 GHz), since below it every reported quantity but alpha is NaN.
console.log('\n=== rectangular waveguide: TE wall impedance ===');
checkRatio('WR-90 rough 1um', await internalL(waveguide(), WG_FREQS),
    await internalL(waveguide({ rq: RQ }), WG_FREQS),
    (f) => calculate_Zrough(f, SIG, RQ), WG_FREQS);

console.log(failures === 0 ? '\nALL TESTS PASSED' : `\n${failures} CHECK(S) FAILED`);
process.exit(failures === 0 ? 0 : 1);
