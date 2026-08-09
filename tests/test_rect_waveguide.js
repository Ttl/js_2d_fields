// Rectangular waveguide against CLOSED FORMS — the analytic-truth anchor for the
// SOURCE-FREE (non-TEM) solve path, in the same spirit as tests/test_coax.js.
//
// A waveguide is the sharpest validation target in the suite after coax, and it exercises
// a path nothing else does: there is no driven conductor, so no static solve, no quasi-TEM
// overlap pick and no capacitance-derived RLGC. Everything is anchored on the eigensolve's
// own cutoff wavenumber, and every reported quantity has an exact closed form:
//
//   kc      = pi/max(a,b)                                   (geometric — no er)
//   fc      = c0*kc/(2*pi*sqrt(er))
//   beta    = sqrt(k0^2*er - kc^2)
//   eps_eff = er - (kc/k0)^2
//   Z_TE    = omega*mu0/beta = k0*eta0/beta                 (independent of er)
//   Zpv     = 2*min(a,b)/max(a,b) * Z_TE                    (power-voltage definition)
//   a_d     = k0^2*er*tand/(2*beta)
//   a_c     = Rs*(1 + 2*(b/a)*(fc/f)^2) / (b*eta*sqrt(1-(fc/f)^2)),  eta = eta0/sqrt(er)
//   a_ev    = sqrt(kc^2 - k0^2*er)                          (below cutoff)
//
// (a_c is algebraically identical to Pozar 3.96, Rs(2*b*pi^2 + a^3*k^2)/(a^3*b*beta*k*eta).)
//
// Unlike every other medium here, alpha_c is NOT the slow-converging quantity: the walls
// are straight and the mode is a smooth sinusoid, so there is no corner singularity and no
// polygonized curve. Measured, alpha_c is converged to five figures on the coarsest usable
// mesh while kc is still moving — see the table in RectWaveguideSolver's constructor. The
// tolerances below are set with margin over the measured errors (all under 0.05%), so a
// failure means a real regression, not mesh noise.
//
// Geometry construction and parameter validation live in tests/test_geometry.js (pure JS,
// milliseconds); this file owns the physics.
//
// Run: node tests/test_rect_waveguide.js
import { RectWaveguideSolver } from '../src/rect_waveguide.js';
import { calculate_Zrough, calculate_Zrough_layered } from '../src/surface_roughness.js';
import { computeSParamsSingleEnded, sParamTodB } from '../src/sparameters.js';

const EPS0 = 8.854187817e-12, MU0 = 4e-7 * Math.PI;
const C0 = 1 / Math.sqrt(EPS0 * MU0), ETA0 = Math.sqrt(MU0 / EPS0);
const NP_TO_DB = 8.685889638;

// WR-90 (X band): fc = 6.557 GHz, second cutoff 13.114 GHz, published band 8.2-12.4 GHz.
const A = 22.86e-3, B = 10.16e-3, SIG = 5.8e7;
const KAPPA = 2 * B / A;                     // Zpv / Z_TE

let failures = 0;
function check(name, cond, detail = '') {
    console.log(`${cond ? '✓ PASS' : '✗ FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
    if (!cond) failures++;
}
const pct = (a_, b_) => ((a_ - b_) / b_ * 100).toFixed(4) + '%';
function near(name, got, want, tolPct) {
    check(name, Math.abs(got / want - 1) * 100 <= tolPct,
        `${got.toExponential(5)} vs ${want.toExponential(5)}, ${pct(got, want)}, tol ${tolPct}%`);
}

// --- closed forms, parameterised by fill ---
const kc = Math.PI / A;
const forms = (er, tand) => {
    const fc = C0 * kc / (2 * Math.PI * Math.sqrt(er));
    const eta = ETA0 / Math.sqrt(er);
    return {
        fc, kc,
        beta: (f) => Math.sqrt((2 * Math.PI * f / C0) ** 2 * er - kc * kc),
        epsEff: (f) => er - (kc / (2 * Math.PI * f / C0)) ** 2,
        ZTE: (f) => 2 * Math.PI * f * MU0 / Math.sqrt((2 * Math.PI * f / C0) ** 2 * er - kc * kc),
        alphaD: (f) => (2 * Math.PI * f / C0) ** 2 * er * tand
            / (2 * Math.sqrt((2 * Math.PI * f / C0) ** 2 * er - kc * kc)),
        alphaC: (f) => {
            const Rs = 1 / (SIG * Math.sqrt(2 / (2 * Math.PI * f * MU0 * SIG)));
            const r = (fc / f) ** 2;
            return Rs * (1 + 2 * (B / A) * r) / (B * eta * Math.sqrt(1 - r));
        },
        alphaEv: (f) => Math.sqrt(kc * kc - (2 * Math.PI * f / C0) ** 2 * er),
    };
};

const quiet = async (fn) => {
    const log = console.log; console.log = () => {};
    try { return await fn(); } finally { console.log = log; }
};

async function build(opts = {}, hScale = 1) {
    const s = new RectWaveguideSolver({ width: A, height: B, sigma_cond: SIG, freq: 10e9, ...opts });
    // Djordjevic-Sarkar would shift er away from the nominal value and invalidate every
    // closed form above.
    s.use_causal_materials = false;
    if (hScale !== 1) {
        s.tri_mesh_hints = { hFine: s.tri_mesh_hints.hFine * hScale,
                             hCoarse: s.tri_mesh_hints.hCoarse * hScale };
    }
    await quiet(() => s.solve_adaptive({ max_nodes: 20000 }));
    return s;
}
const at = async (s, f) => (await quiet(() => s.computeAtFrequency(f))).modes[0];
const warnsOf = async (s, f) => (await quiet(() => s.computeAtFrequency(f))).warnings || [];

// ---------------------------------------------------------------- air fill, in band
console.log('\n=== WR-90, air fill ===');
const air = forms(1.0, 0);
const s = await build();
console.log(`  mesh: ${s.triMesh.nTris} triangles, half-domain=${s._triBackend.symmetry}`);

// The eigensolve's own cutoff must agree with pi/a, or the analytic fallback would have
// engaged and every number below would be analytic-vs-analytic rather than a real check.
{
    const raw = s._triBackend._wg.kcRaw;
    check('kc from the field solve is not gated', s._triBackend._wg.gated === false);
    near('kc (field solve vs pi/a)', raw, kc, 0.05);
}

for (const f of [8.2e9, 10e9, 12.4e9]) {
    const m = await at(s, f);
    console.log(`\n-- ${(f / 1e9).toFixed(2)} GHz --`);
    check(`${f / 1e9} GHz: propagating`, m.belowCutoff === false);
    near('fc', m.fc, air.fc, 0.05);
    near('beta', m.beta, air.beta(f), 0.2);
    near('eps_eff', m.eps_eff, air.epsEff(f), 0.3);
    near('Z_TE', m.Z_TE, air.ZTE(f), 0.5);
    near('Zc.re (= Zpv)', m.Zc.re, KAPPA * air.ZTE(f), 0.5);
    near('alpha_c', m.alpha_c, air.alphaC(f) * NP_TO_DB, 1.0);
    // Equivalent-circuit elements, normalized to Zpv. L carries the wall surface reactance
    // on top of kappa*mu0, which is a ~3e-4 relative increment for copper.
    check('L_external is exactly kappa*mu0',
        Math.abs(m.L_external / (KAPPA * MU0) - 1) < 1e-12, m.L_external.toExponential(6));
    near('RLGC.L', m.RLGC.L, KAPPA * MU0, 0.5);
    near('RLGC.C', m.RLGC.C, (EPS0 - kc * kc / ((2 * Math.PI * f) ** 2 * MU0)) / KAPPA, 0.5);
    check('RLGC.G is zero for a lossless fill', m.RLGC.G === 0);
    near('RLGC.R = 2*alpha_c*Zpv', m.RLGC.R,
        2 * (m.alpha_c / NP_TO_DB) * m.Zc.re, 0.5);

    // THE S-parameter contract: computeSParamsSingleEnded rebuilds gamma and Zc from RLGC
    // alone, so if this identity fails the exported Touchstone silently disagrees with the
    // Results tab. Checked here rather than assumed.
    const { R, L, G, C } = m.RLGC, w = 2 * Math.PI * f;
    const zr = R, zi = w * L, yr = G, yi = w * C;
    const gr = zr * yr - zi * yi, gi = zr * yi + zi * yr;
    const mag = Math.hypot(gr, gi);
    const gRe = Math.sqrt((mag + gr) / 2), gIm = Math.sqrt((mag - gr) / 2);
    near('gamma from RLGC: alpha', gRe, (m.alpha_c + m.alpha_d) / NP_TO_DB, 1.0);
    near('gamma from RLGC: beta', gIm, m.beta, 0.2);
    const zcMag = Math.hypot(zr, zi) / Math.hypot(yr, yi);
    near('Zc from RLGC', Math.sqrt(zcMag), m.Zc.abs(), 1e-6);
}

// ---------------------------------------------------------------- dielectric fill
// Exercises alpha_d and the 1/sqrt(er) scaling of the cutoff. kc itself must NOT move —
// that invariant is what lets the backend solve the eigenproblem once and propagate
// beta(f) analytically even under causal materials.
console.log('\n=== WR-90, PTFE fill (er=2.1, tand=2e-4) ===');
{
    const ER = 2.1, TAND = 2e-4, f = 10e9;
    const ptfe = forms(ER, TAND);
    const sd = await build({ epsilon_r: ER, tan_delta: TAND });
    const m = await at(sd, f);
    near('kc is unchanged by the fill', m.kc, kc, 0.05);
    near('fc scales as 1/sqrt(er)', m.fc, ptfe.fc, 0.05);
    near('beta', m.beta, ptfe.beta(f), 0.2);
    near('eps_eff', m.eps_eff, ptfe.epsEff(f), 0.3);
    near('alpha_d', m.alpha_d, ptfe.alphaD(f) * NP_TO_DB, 0.5);
    near('alpha_c', m.alpha_c, ptfe.alphaC(f) * NP_TO_DB, 1.0);
    near('RLGC.G = omega*eps0*er*tand/kappa', m.RLGC.G,
        2 * Math.PI * f * EPS0 * ER * TAND / KAPPA, 0.5);
    near('alpha_total', m.alpha_total, (ptfe.alphaC(f) + ptfe.alphaD(f)) * NP_TO_DB, 1.0);
}

// ---------------------------------------------------------------- below cutoff
// The reported quantities are NaN by design (it keeps the Results plots' autoscale clean
// and breaks the trace at cutoff); only the attenuation is meaningful. The S-parameter
// paths key on belowCutoff to drop these points entirely.
console.log('\n=== below cutoff ===');
{
    const f = 5e9;
    const m = await at(s, f);
    check('belowCutoff flag set', m.belowCutoff === true);
    near('alpha_total = 8.686*sqrt(kc^2-k0^2)', m.alpha_total, air.alphaEv(f) * NP_TO_DB, 0.5);
    check('eps_eff / Zc / RLGC are NaN below cutoff',
        Number.isNaN(m.eps_eff) && Number.isNaN(m.Zc.re) && Number.isNaN(m.Z0) &&
        Number.isNaN(m.RLGC.R) && Number.isNaN(m.RLGC.L) &&
        Number.isNaN(m.RLGC.G) && Number.isNaN(m.RLGC.C));
    check('below-cutoff warning raised',
        (await warnsOf(s, f)).some(w => w.type === 'wg-cutoff'));
}

// ---------------------------------------------------------------- warnings
console.log('\n=== single-mode band ===');
{
    check('no warnings in band', (await warnsOf(s, 10e9)).length === 0,
        JSON.stringify((await warnsOf(s, 10e9)).map(w => w.type)));
    // Above the TE20/TE01 cutoff the guide is over-moded and this solver still reports only
    // the fundamental — the limitation the help text states, asserted here.
    check('over-moded warning above the second cutoff',
        (await warnsOf(s, 15e9)).some(w => w.type === 'wg-overmoded'));
}

// ---------------------------------------------------------------- S-parameters
// Referenced to its own characteristic impedance a uniform line is matched BY DEFINITION,
// so S11 must vanish. That only comes out of the ABCD->S conversion when the reference is
// the true COMPLEX Zc: Zc of any lossy line has a small imaginary part, and referencing to
// Re(Zc) alone leaves |S11| ~ |Im(Zc)|/(2*Re(Zc)) — around -80 dB here, which reads as a
// physical reflection when there is none. This pins the general behaviour (complex Zr),
// not a waveguide special case.
console.log('\n=== self-referenced S-parameters ===');
{
    const LEN = 0.01;   // 10 mm
    for (const f of [8.2e9, 10e9, 12.4e9]) {
        const g = (f / 1e9).toFixed(1) + ' GHz';
        const m = await at(s, f);
        const sp = computeSParamsSingleEnded(f, m.RLGC, LEN, m.Zc);
        check(`${g}: S11 vanishes to the -300 dB floor when referenced to the complex Zc`,
            sParamTodB(sp.S11) === -300, `${sParamTodB(sp.S11).toFixed(1)} dB, |S11|=${sp.S11.abs().toExponential(2)}`);
        check(`${g}: S22 likewise`, sParamTodB(sp.S22) === -300,
            sParamTodB(sp.S22).toFixed(1) + ' dB');
        // Referencing to Re(Zc) instead must be visibly WORSE — this is the regression the
        // fix addresses, so assert the difference rather than only the fixed value.
        const real = computeSParamsSingleEnded(f, m.RLGC, LEN, m.Zc.re);
        check(`${g}: a real-only reference would leave a spurious S11`,
            sParamTodB(real.S11) > -120, sParamTodB(real.S11).toFixed(1) + ' dB');
        // The matched-section response: attenuation and phase over the length, nothing else.
        near(`${g}: |S21| = exp(-alpha*L)`, sp.S21.abs(),
            Math.exp(-(m.alpha_total / NP_TO_DB) * LEN), 0.001);
        // beta*L < pi at these frequencies and this length, so no phase wrap to undo.
        near(`${g}: arg(S21) = -beta*L`, sp.S21.arg(), -m.beta * LEN, 0.1);
        check(`${g}: reciprocal (S12 === S21)`, sp.S12 === sp.S21);
    }
}

// ---------------------------------------------------------------- roughness and plating
// The wall is ONE continuous surface, so both enter as a single scalar Zs.re/Rs scaling of
// alpha_c. This is the assertion that catches a silent fall-through to bare metal: the
// per-(rect, face) surface-group machinery finds no rects here, so the waveguide path must
// apply the surface impedance itself.
console.log('\n=== surface finish ===');
{
    const f = 10e9, RQ = 1e-6;
    const Rs = 1 / (SIG * Math.sqrt(2 / (2 * Math.PI * f * MU0 * SIG)));
    const base = await at(s, f);

    const rough = await at(await build({ rq: RQ }), f);
    near('roughness scales alpha_c by Zs.re/Rs', rough.alpha_c / base.alpha_c,
        calculate_Zrough(f, SIG, RQ).re / Rs, 0.5);
    check('roughness raises the internal inductance (slows the mode)',
        rough.L_internal > base.L_internal && rough.eps_eff > base.eps_eff,
        `L_int ${base.L_internal.toExponential(3)} -> ${rough.L_internal.toExponential(3)}`);

    const PL = { sigma: 4.5e7, thickness: 2e-6, rq: 0.5e-6, top: true, sides: true, bottom: true };
    const plated = await at(await build({ rq: RQ, plating: PL }), f);
    near('plating scales alpha_c by the layered Zs.re/Rs', plated.alpha_c / base.alpha_c,
        calculate_Zrough_layered(f, SIG, PL.rq, PL.sigma, PL.thickness).re / Rs, 0.5);
}

// ---------------------------------------------------------------- mesh independence
// A different mesh either way, so agreement between them is a mesh-independent statement.
console.log('\n=== mesh independence ===');
{
    const f = 10e9;
    const coarse = await build({}, 2);
    const mf = await at(s, f), mc = await at(coarse, f);
    console.log(`  ${s.triMesh.nTris} vs ${coarse.triMesh.nTris} triangles`);
    check('the coarse mesh really is coarser', coarse.triMesh.nTris < s.triMesh.nTris / 2);
    near('fine vs coarse: eps_eff', mc.eps_eff, mf.eps_eff, 0.1);
    near('fine vs coarse: alpha_c', mc.alpha_c, mf.alpha_c, 2.0);
    near('fine vs coarse: Zc', mc.Zc.re, mf.Zc.re, 0.5);
}

// ---------------------------------------------------------------- b > a  (TE01)
// max(a, b) has to pick the HEIGHT here, in both the analytic cutoff and the seed field.
// Rotating the guide 90 degrees must reproduce the same physics exactly.
console.log('\n=== b > a (fundamental is TE01) ===');
{
    const f = 10e9;
    const tall = await build({ width: B, height: A });
    const m = await at(tall, f);
    near('kc = pi/b', m.kc, kc, 0.05);
    near('fc matches the a > b guide', m.fc, air.fc, 0.05);
    near('eps_eff matches the a > b guide', m.eps_eff, air.epsEff(f), 0.3);
    near('Zpv matches (kappa = 2*min/max either way)', m.Zc.re, KAPPA * air.ZTE(f), 0.5);
}

console.log(`\n${failures === 0 ? 'ALL TESTS PASSED' : failures + ' TEST(S) FAILED'}`);
process.exit(failures === 0 ? 0 : 1);
