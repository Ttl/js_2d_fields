// sparameters.js - S-Parameter Calculations for Transmission Lines
import { Complex } from './complex.js';
import { Matrix2x2, mat2Mul, mat2Inv, mat2T } from './matrix.js';

/**
 * Compute 2-port S-parameters for a single-ended transmission line
 *
 * @param {number} freq - Frequency in Hz
 * @param {object} rlgc - RLGC parameters {R, L, G, C} in SI units
 * @param {number} length - Line length in meters
 * @param {number|Complex} Z_ref - Reference impedance in Ohms (typically the real 50).
 *     A Complex is accepted, for a line referenced to its own slightly-complex Zc.
 * @returns {object} - {S11, S21, S12, S22} as Complex numbers
 */
function computeSParamsSingleEnded(freq, rlgc, length, Z_ref) {
    const omega = 2 * Math.PI * freq;
    const { R, L, G, C } = rlgc;

    // The reference impedance may be complex. It usually is not, a 50 Ohm port is real,
    // but a line referenced to its own Zc needs the true complex value. Zc of any lossy
    // line has a small imaginary part, and the B/Zr and C*Zr terms below only cancel (so a
    // matched line gives S11 = 0) when Zr is that exact value. Referencing to Re(Zc)
    // instead leaves |S11| ~ |Im(Zc)|/(2*Re(Zc)), which looks like a real reflection.
    //
    // WAVE DEFINITION. This is the pseudo-wave / traveling-wave conversion, and
    // those two match here for a line, so Gamma is (Z - Zr)/(Z + Zr) either
    // way. Both are the same diagonal similarity transform of
    // M = (Z - ZR)(Z + ZR)^-1, differing only in which diagonal matrix
    // conjugates it:
    //     pseudo    (Marks & Williams 1992)   S = U    M U^-1,    U    = diag(sqrt(Re z0)/|z0|)
    //     traveling                           S = Y^.5 M Y^-.5,   Y^.5 = diag(1/sqrt(z0))
    // When every port shares one reference, both conjugating matrices are a scalar times
    // the identity, the similarity collapses, and S = M exactly, for any complex Zr and
    // any network, not just a uniform line. With per-port references they separate, but
    // only in the transmission terms: S[m][n] picks up U_m/U_n vs sqrt(z0_n/z0_m), and
    // those ratios are both 1 on the diagonal, so reflection terms agree regardless.
    //
    // Power waves (Kurokawa 1965) are the genuinely different convention: b conjugates the
    // reference, Gamma = (Z - Zr*)/(Z + Zr). That is what you want for conjugate-matching
    // problems, and it is wrong here. Referenced to its own Z0 a perfectly uniform line
    // would report |Gamma| = |Im Z0|/|Z0| rather than zero.
    //
    // Cross-checked against scikit-rf z2s(..., s_def=...) on a lossy waveguide
    // 2-port with Zr = Z0 complex: 'pseudo' and 'traveling' each agree with
    // this function to float precision, while 'power' differs by |Im Z0/Re Z0|.
    // With a real Zr all three agree exactly — skrf's own s2s() short-circuits
    // on real z0 for precisely that reason.
    //
    // Duck-typed, not `instanceof Complex`: build.sh ships the lazily-imported tri_solver
    // tree with its own copy of complex.js alongside the one esbuild inlines into the app
    // bundle, so a Zc built inside the backend is not an instance of the class this module
    // sees. `instanceof` would silently take the number branch there and produce NaN
    // S-parameters in the built app while dev (single module instance) looked perfect.
    // Rebuilding as a local Complex also guarantees the arithmetic below is this copy's.
    const Zr = (typeof Z_ref === 'number')
        ? new Complex(Z_ref, 0)
        : new Complex(Z_ref.re, Z_ref.im);

    // Handle DC case (frequency = 0)
    // At DC, the transmission line behaves as a simple series resistance R*length
    if (freq === 0 || omega === 0) {
        const R_total = R * length;

        // ABCD matrix for series resistance: A = 1, B = R_total, C = 0, D = 1
        // den = A + B/Zr + C*Zr + D = 2 + R_total/Zr
        // The real case stays in scalar arithmetic: complex division computes
        // (a*c)/(c^2+d^2) rather than a/c, which differs by an ULP and would perturb every
        // existing DC result for no reason.
        const BoverZ = Zr.im === 0
            ? new Complex(R_total / Zr.re, 0)
            : new Complex(R_total, 0).div(Zr);
        const den = BoverZ.add(new Complex(2, 0));

        // S11 = (A + B/Zr - C*Zr - D) / den = (R_total/Zr) / den
        const S11 = BoverZ.div(den);

        // S21 = 2 / den
        const S21 = new Complex(2, 0).div(den);

        const S12 = S21;  // Reciprocal

        // S22 = (-A + B/Zr - C*Zr + D) / den = (R_total/Zr) / den
        const S22 = BoverZ.div(den);

        return { S11, S21, S12, S22 };
    }

    // Series impedance per unit length: Z = R + jwL
    const Z_per_length = new Complex(R, omega * L);

    // Shunt admittance per unit length: Y = G + jwC
    const Y_per_length = new Complex(G, omega * C);

    // Propagation constant: gamma = sqrt(Z * Y)
    const gamma = Z_per_length.mul(Y_per_length).sqrt();

    // Characteristic impedance: Z0 = sqrt(Z / Y)
    const Z0 = Z_per_length.div(Y_per_length).sqrt();

    // gamma * length
    const gl = gamma.mul(length);

    // ABCD matrix elements
    // A = cosh(gamma * l)
    // B = Z0 * sinh(gamma * l)
    // C = sinh(gamma * l) / Z0
    // D = cosh(gamma * l)
    const A = gl.cosh();
    const B = Z0.mul(gl.sinh());
    const C_abcd = gl.sinh().div(Z0);
    const D = gl.cosh();

    // Convert ABCD to S-parameters (Zr built above, possibly complex)
    // Common denominator: A + B/Zr + C*Zr + D
    const den = A.add(B.div(Zr)).add(C_abcd.mul(Zr)).add(D);

    // S11 = (A + B/Zr - C*Zr - D) / den
    const S11 = A.add(B.div(Zr)).sub(C_abcd.mul(Zr)).sub(D).div(den);

    // S21 = 2 / den (for reciprocal network)
    const S21 = new Complex(2, 0).div(den);

    // S12 = 2 * (AD - BC) / den = 2 / den (for reciprocal lossless or lossy uniform line, AD-BC=1)
    // For a uniform transmission line, det(ABCD) = 1, so S12 = S21
    const S12 = S21;

    // S22 = (-A + B/Zr - C*Zr + D) / den
    const S22 = A.neg().add(B.div(Zr)).sub(C_abcd.mul(Zr)).add(D).div(den);

    return { S11, S21, S12, S22 };
}

/**
 * Compute modal characteristic impedance from RLGC parameters
 * @param {number} freq - Frequency in Hz
 * @param {object} rlgc - RLGC parameters {R, L, G, C}
 * @returns {number} - Modal Z0 magnitude in Ohms
 */
function computeZ0(freq, rlgc) {
    const omega = 2 * Math.PI * freq;
    const { R, L, G, C } = rlgc;

    // Handle DC case (frequency = 0)
    // At DC, Z0 = sqrt(R/G) with G=0 gives infinity
    if (freq === 0 || omega === 0) {
        return 1e12;  // Return very large impedance for DC
    }

    // Series impedance per unit length: Z = R + jwL
    const Z_per_length = new Complex(R, omega * L);

    // Shunt admittance per unit length: Y = G + jwC
    const Y_per_length = new Complex(G, omega * C);

    // Characteristic impedance: Z0 = sqrt(Z / Y)
    const Z0 = Z_per_length.div(Y_per_length).sqrt();

    return Z0.abs();
}

/**
 * Compute 4-port S-parameters for a differential transmission line
 * Ports 1,3 are at one end (positive and negative), Ports 2,4 at the other end
 *
 * This function uses modal decomposition with odd and even mode parameters.
 *
 * @param {number} freq - Frequency in Hz
 * @param {object} rlgc_odd - RLGC parameters for odd mode {R, L, G, C}
 * @param {object} rlgc_even - RLGC parameters for even mode {R, L, G, C}
 * @param {number} length - Line length in meters
 * @param {number} Z_ref - Reference impedance in Ohms (typically 50)
 * @returns {object} - {S: 4x4 array of Complex, SDD11, SDD21, SCC11, SCC21, SDC11, SDC21, SCD11, SCD21}
 */
function computeSParamsDifferential(freq, rlgc_odd, rlgc_even, length, Z_ref) {
    // Compute single-ended S-params for each mode using system reference impedance
    const S_odd = computeSParamsSingleEnded(freq, rlgc_odd, length, Z_ref);
    const S_even = computeSParamsSingleEnded(freq, rlgc_even, length, Z_ref);

    // For ideal symmetric differential pairs with no coupling between modes,
    // the 4-port S-matrix can be constructed from odd and even mode responses.
    //
    // Port assignment:
    //   Port 1 = near end, trace + (in+)
    //   Port 2 = near end, trace - (in-)
    //   Port 3 = far end, trace + (out+)
    //   Port 4 = far end, trace - (out-)
    //
    // Modal decomposition for symmetric coupled lines:
    //   When exciting port 1 (a1=1, others=0):
    //   - Odd mode excitation: a_odd = (a1 - a2)/√2 = 1/√2
    //   - Even mode excitation: a_even = (a1 + a2)/√2 = 1/√2
    //
    // S-parameter formulas (for symmetric coupled lines):
    //   S11 = S22 = (S_odd_11 + S_even_11) / 2  (reflection at near end)
    //   S33 = S44 = (S_odd_22 + S_even_22) / 2  (reflection at far end)
    //   S21 = S12 = (S_even_11 - S_odd_11) / 2  (near-end coupling, NEXT)
    //   S43 = S34 = (S_even_22 - S_odd_22) / 2  (far-end reflection coupling)
    //   S31 = S13 = S42 = S24 = (S_odd_21 + S_even_21) / 2  (through transmission)
    //   S41 = S14 = S32 = S23 = (S_even_21 - S_odd_21) / 2  (far-end coupling, FEXT)

    const half = new Complex(0.5, 0);

    // Calculate 4-port S-parameters
    const Snn = S_odd.S11.add(S_even.S11).mul(half);   // Near-end reflection (S11, S22)
    const Sff = S_odd.S22.add(S_even.S22).mul(half);   // Far-end reflection (S33, S44)
    const Snext = S_even.S11.sub(S_odd.S11).mul(half); // Near-end coupling (S21, S12)
    const Sfref = S_even.S22.sub(S_odd.S22).mul(half); // Far-end reflection coupling (S43, S34)
    const Sthru = S_odd.S21.add(S_even.S21).mul(half); // Through transmission (S31, S13, S42, S24)
    const Sfext = S_even.S21.sub(S_odd.S21).mul(half); // Far-end coupling (S41, S14, S32, S23)

    // Build 4x4 matrix (symmetric coupled line)
    const S = [
        [Snn,   Snext, Sthru, Sfext],  // Row 1: S11, S12, S13, S14
        [Snext, Snn,   Sfext, Sthru],  // Row 2: S21, S22, S23, S24
        [Sthru, Sfext, Sff,   Sfref],  // Row 3: S31, S32, S33, S34
        [Sfext, Sthru, Sfref, Sff  ]   // Row 4: S41, S42, S43, S44
    ];

    // Mixed-mode S-parameters (for plotting)
    // SDD11 = S_odd_11 (differential reflection)
    // SDD21 = S_odd_21 (differential transmission)
    // SCC11 = S_even_11 (common-mode reflection)
    // SCC21 = S_even_21 (common-mode transmission)
    // SDC, SCD = 0 for symmetric line

    return {
        S,
        SDD11: S_odd.S11,
        SDD21: S_odd.S21,
        SCC11: S_even.S11,
        SCC21: S_even.S21,
        SDC11: new Complex(0, 0),
        SDC21: new Complex(0, 0),
        SCD11: new Complex(0, 0),
        SCD21: new Complex(0, 0)
    };
}

// Modal→physical for a per-unit-length quantity, given the unit-norm voltage-eigenvector matrix
// Tv (columns = modes). Series quantities (R, L) transform as Tv·diag·Tvᵀ; shunt quantities
// (G, C) as Tv⁻ᵀ·diag·Tv⁻¹ (verified to reproduce the Maxwell [C] exactly). Both reduce to the
// symmetric reconstruction X11=X22=(Xo+Xe)/2, X12=(Xe−Xo)/2 when Tv = [1,∓1]/√2.
function _modalToPhys2x2(Xo, Xe, Tv, shunt) {
    const D = [[Xo, 0], [0, Xe]];
    const A = shunt ? mat2T(mat2Inv(Tv)) : Tv;        // left factor
    const B = shunt ? mat2Inv(Tv) : mat2T(Tv);        // right factor
    return mat2Mul(mat2Mul(A, D), B);
}

/**
 * The physical 2×2 per-unit-length [R][L][G][C] matrices for a coupled pair, EXACTLY as fed to
 * the MTL 4-port S-parameter computation. With a physMatrix ({C,L,Tv?}) the reactive [C]/[L] are
 * the genuine (asymmetric) Maxwell matrices, and [R]/[G] use the modal-eigenvector transform when
 * Tv is present (else the symmetric reconstruction, which is what a velocity-degenerate pair
 * gets). Without a physMatrix every quantity uses the symmetric odd/even reconstruction
 * (X11=X22=(Xo+Xe)/2, X12=(Xe−Xo)/2). Shared by computeSParamsDiffAuto and the RLGC_matrix result
 * field so the displayed/exported matrix matches the one that actually drives the S-parameters.
 * @param {object} rlgc_odd  per-mode scalars {R,L,G,C} for the odd mode
 * @param {object} rlgc_even per-mode scalars {R,L,G,C} for the even mode
 * @param {object|null} physMatrix {C, L, Tv?:[[vo],[ve]]} physical p.u.l. matrices + eigenvectors
 * @returns {{R:number[][], L:number[][], G:number[][], C:number[][]}} physical 2×2 matrices
 */
function buildPhysicalRLGC(rlgc_odd, rlgc_even, physMatrix) {
    const sym = (o, e) => [[(o + e) / 2, (e - o) / 2], [(e - o) / 2, (o + e) / 2]];
    if (!physMatrix || !physMatrix.C || !physMatrix.L) {
        return { R: sym(rlgc_odd.R, rlgc_even.R), L: sym(rlgc_odd.L, rlgc_even.L),
                 G: sym(rlgc_odd.G, rlgc_even.G), C: sym(rlgc_odd.C, rlgc_even.C) };
    }
    let R, G;
    if (physMatrix.Tv) {
        const [vo, ve] = physMatrix.Tv;               // eigenvectors (any norm); use unit-norm columns
        const no = Math.hypot(vo[0], vo[1]) || 1, ne = Math.hypot(ve[0], ve[1]) || 1;
        const Tv = [[vo[0] / no, ve[0] / ne], [vo[1] / no, ve[1] / ne]];
        R = _modalToPhys2x2(rlgc_odd.R, rlgc_even.R, Tv, false);    // series
        G = _modalToPhys2x2(rlgc_odd.G, rlgc_even.G, Tv, true);     // shunt
    } else {
        R = sym(rlgc_odd.R, rlgc_even.R);
        G = sym(rlgc_odd.G, rlgc_even.G);
    }
    return { R, L: physMatrix.L, G, C: physMatrix.C };
}

/**
 * Differential S-parameters, automatically using the full MTL 4-port when the physical 2×2
 * [C]/[L] matrices are available (asymmetric-correct) and falling back to the odd/even
 * combination otherwise.
 * @param {object} physMatrix {C, L, Tv?:[[vo],[ve]]} physical p.u.l. matrices + eigenvectors, or null
 */
function computeSParamsDiffAuto(freq, rlgc_odd, rlgc_even, physMatrix, length, Z_ref) {
    if (!physMatrix || !physMatrix.C || !physMatrix.L) {
        return computeSParamsDifferential(freq, rlgc_odd, rlgc_even, length, Z_ref);
    }
    const { R, L, G, C } = buildPhysicalRLGC(rlgc_odd, rlgc_even, physMatrix);
    return computeSParamsDifferentialMTL(freq, R, L, G, C, length, Z_ref);
}

// Build a complex 2×2 from real part [[..]] and imaginary part [[..]].
function mComplex(re, im) {
    return new Matrix2x2(
        new Complex(re[0][0], im[0][0]), new Complex(re[0][1], im[0][1]),
        new Complex(re[1][0], im[1][0]), new Complex(re[1][1], im[1][1]));
}

/**
 * Full multiconductor (MTL) 4-port S-parameters for a coupled pair from the physical 2×2
 * per-unit-length matrices [R],[L],[G],[C]. Unlike computeSParamsDifferential (which combines
 * odd/even single-ended responses and is exact only for a SYMMETRIC pair), this handles an
 * asymmetric pair: it yields S11≠S22 and non-zero mode-conversion SDC/SCD.
 *
 * Method: [Z]=[R]+jω[L], [Y]=[G]+jω[C]; [γ]=√([Z][Y]); chain blocks A=cosh(γℓ),
 * B=γ⁻¹sinh(γℓ)·[Z], C=[Y]·γ⁻¹sinh(γℓ), D=Aᵀ (reciprocal); block ABCD→4-port S; mixed-mode
 * transform → SDD/SCC/SDC/SCD. Port order: {1,2}=near ends of trace1/2, {3,4}=far ends.
 *
 * @param {number} freq Hz; @param {number[][]} R2,L2,G2,C2 physical 2×2 p.u.l. matrices (SI)
 * @param {number} length m; @param {number} Z_ref Ω
 */
function computeSParamsDifferentialMTL(freq, R2, L2, G2, C2, length, Z_ref) {
    const omega = 2 * Math.PI * freq;
    const Z = mComplex(R2, L2.map(row => row.map(v => v * omega)));   // [R] + jω[L]
    const Y = mComplex(G2, C2.map(row => row.map(v => v * omega)));   // [G] + jω[C]
    const gamma = Z.mul(Y).sqrt();                                     // [γ] = √([Z][Y])
    const gl = gamma.mul(length);
    const sinch = gamma.inv().mul(gl.sinh());                          // γ⁻¹·sinh(γℓ)
    const A = gl.cosh();
    const B = sinch.mul(Z);   // γ⁻¹sinh(γℓ)·[Z] — γ and [Z] do NOT commute for an asymmetric pair
    const Cb = Y.mul(sinch);
    const D = A.transpose();                                           // reciprocal: D = Aᵀ
    // Block ABCD → S, Z_ref·I on every port. Operand order matters for the non-commuting
    // blocks: from the wave equations, S11 = (A+BY−CZ−D)·den⁻¹ (den⁻¹ on the RIGHT) while
    // S22 = den⁻¹·(−A+BY−CZ+D) (den⁻¹ on the LEFT), with den = A+BY+CZ+D.
    const Y0 = 1 / Z_ref;
    const BY = B.mul(Y0), CZ = Cb.mul(Z_ref);
    const den = A.add(BY).add(CZ).add(D);
    const denI = den.inv();
    const S11 = A.add(BY).sub(CZ).sub(D).mul(denI);                   // near-near
    const S22 = denI.mul(A.mul(-1).add(BY).sub(CZ).add(D));          // far-far
    const S21 = denI.mul(2);                                           // far-near: 2·den⁻¹
    const S12 = S21;   // = S21ᵀ by reciprocity; den is symmetric (A+Aᵀ, B, C all symmetric)
    // Assemble single-ended 4×4 (ports 1,2 near; 3,4 far), then mixed-mode transform.
    const blk = (m) => [[m.a, m.b], [m.c, m.d]];
    const b11 = blk(S11), b12 = blk(S12), b21 = blk(S21), b22 = blk(S22);
    const S = [
        [b11[0][0], b11[0][1], b12[0][0], b12[0][1]],
        [b11[1][0], b11[1][1], b12[1][0], b12[1][1]],
        [b21[0][0], b21[0][1], b22[0][0], b22[0][1]],
        [b21[1][0], b21[1][1], b22[1][0], b22[1][1]],
    ];
    // M: [V1,V2,V3,V4] → [d1,d2,c1,c2], rows = (1/√2)·{[1,-1,0,0],[0,0,1,-1],[1,1,0,0],[0,0,1,1]}.
    const r = 1 / Math.SQRT2;
    const M = [[r, -r, 0, 0], [0, 0, r, -r], [r, r, 0, 0], [0, 0, r, r]];
    const mm = Array.from({ length: 4 }, () => Array.from({ length: 4 }, () => new Complex(0, 0)));
    for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) {
        let s = new Complex(0, 0);
        for (let a = 0; a < 4; a++) for (let b = 0; b < 4; b++)
            s = s.add(S[a][b].mul(M[i][a] * M[j][b]));                // (M·S·Mᵀ)_ij
        mm[i][j] = s;
    }
    // mm order [d1,d2,c1,c2].
    return {
        S,
        SDD11: mm[0][0], SDD21: mm[1][0],
        SCC11: mm[2][2], SCC21: mm[3][2],
        SDC11: mm[0][2], SDC21: mm[1][2],
        SCD11: mm[2][0], SCD21: mm[3][0],
    };
}

/**
 * Convert Complex S-parameter to dB magnitude
 * @param {Complex} s - S-parameter
 * @returns {number} - Magnitude in dB
 */
function sParamTodB(s) {
    const mag = s.abs();
    if (mag < 1e-15) return -300;  // Prevent log(0)
    return 20 * Math.log10(mag);
}

/**
 * Convert Complex S-parameter to phase in degrees
 * @param {Complex} s - S-parameter
 * @returns {number} - Phase in degrees
 */
function sParamToPhase(s) {
    return s.arg() * 180 / Math.PI;
}

// Self-referenced (non-TEM) media
//
// A rectangular waveguide has no external TEM reference to normalize against.
// Referenced to Z0 instead.
//
// Marked on the mode object by the backend (`self_referenced`), so nothing here has to
// know which media are non-TEM.
function isSelfReferenced(sweepResults) {
    const m = sweepResults && sweepResults[0] && sweepResults[0].result
        && sweepResults[0].result.modes && sweepResults[0].result.modes[0];
    return !!(m && m.self_referenced);
}

// Single-ended S-parameters for one sweep point, picking the reference from the mode
// itself: a non-TEM medium has no external reference, so it normalizes to the mode's own
// characteristic impedance. Every other line uses the caller's fixed Z_ref.
//
// Zc is passed complex. Referencing to Re(Zc) instead leaves a spurious |S11| of roughly
// |Im(Zc)|/(2*Re(Zc)), about -80 dB on a copper waveguide, because the B/Zr and C*Zr
// terms in the T to S conversion only cancel when Zr is the true (complex) Zc.
function sparamsForPoint(freq, result, length, Z_ref) {
    const mode = result.modes[0];
    return computeSParamsSingleEnded(freq, mode.RLGC, length,
        mode.self_referenced ? mode.Zc : Z_ref);
}

// Sweep points that can actually produce S-parameters. Below cutoff a waveguide
// mode is evanescent and its modal impedance is imaginary. Those points are
// dropped. Non-self-referenced media are returned unchanged.
function usableSweepPoints(sweepResults) {
    if (!isSelfReferenced(sweepResults)) return sweepResults;
    return sweepResults.filter(({ result }) => {
        const m = result.modes[0];
        return !m.belowCutoff && m.Zc && isFinite(m.Zc.re) && m.Zc.re > 0;
    });
}

export {
    isSelfReferenced,
    sparamsForPoint,
    usableSweepPoints,
    computeSParamsSingleEnded,
    computeSParamsDifferential,
    computeSParamsDifferentialMTL,
    computeSParamsDiffAuto,
    buildPhysicalRLGC,
    computeZ0,
    sParamTodB,
    sParamToPhase
};
