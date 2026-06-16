import { Complex } from './complex.js';

/**
 * 2x2 Complex Matrix class
 * Matrix is stored as [[a, b], [c, d]] where each element is a Complex number
 */
class Matrix2x2 {
    constructor(a, b, c, d) {
        // Ensure all elements are Complex
        this.a = a instanceof Complex ? a : new Complex(a);
        this.b = b instanceof Complex ? b : new Complex(b);
        this.c = c instanceof Complex ? c : new Complex(c);
        this.d = d instanceof Complex ? d : new Complex(d);
    }

    static identity() {
        return new Matrix2x2(
            new Complex(1, 0), new Complex(0, 0),
            new Complex(0, 0), new Complex(1, 0)
        );
    }

    static zero() {
        return new Matrix2x2(
            new Complex(0, 0), new Complex(0, 0),
            new Complex(0, 0), new Complex(0, 0)
        );
    }

    static scalar(s) {
        // Create a scalar matrix (s * I)
        const sc = s instanceof Complex ? s : new Complex(s);
        return new Matrix2x2(sc, new Complex(0, 0), new Complex(0, 0), sc);
    }

    add(other) {
        return new Matrix2x2(
            this.a.add(other.a), this.b.add(other.b),
            this.c.add(other.c), this.d.add(other.d)
        );
    }

    sub(other) {
        return new Matrix2x2(
            this.a.sub(other.a), this.b.sub(other.b),
            this.c.sub(other.c), this.d.sub(other.d)
        );
    }

    mul(other) {
        if (other instanceof Matrix2x2) {
            // Matrix multiplication
            return new Matrix2x2(
                this.a.mul(other.a).add(this.b.mul(other.c)),
                this.a.mul(other.b).add(this.b.mul(other.d)),
                this.c.mul(other.a).add(this.d.mul(other.c)),
                this.c.mul(other.b).add(this.d.mul(other.d))
            );
        } else {
            // Scalar multiplication
            const s = other instanceof Complex ? other : new Complex(other);
            return new Matrix2x2(
                this.a.mul(s), this.b.mul(s),
                this.c.mul(s), this.d.mul(s)
            );
        }
    }

    det() {
        // Determinant: ad - bc
        return this.a.mul(this.d).sub(this.b.mul(this.c));
    }

    inv() {
        // Inverse: (1/det) * [[d, -b], [-c, a]]
        const det = this.det();
        const invDet = new Complex(1, 0).div(det);
        return new Matrix2x2(
            this.d.mul(invDet), this.b.neg().mul(invDet),
            this.c.neg().mul(invDet), this.a.mul(invDet)
        );
    }

    trace() {
        // Trace: a + d
        return this.a.add(this.d);
    }

    transpose() {
        // Transpose: swap b and c
        return new Matrix2x2(this.a, this.c, this.b, this.d);
    }

    /**
     * Matrix square root using eigendecomposition
     * For a 2x2 matrix A, sqrt(A) such that sqrt(A) * sqrt(A) = A
     */
    sqrt() {
        // For 2x2 matrix, use the formula:
        // sqrt(A) = (1/s) * (A + sqrt(det(A)) * I)
        // where s = sqrt(trace(A) + 2*sqrt(det(A)))

        const det = this.det();
        const sqrtDet = det.sqrt();
        const tr = this.trace();

        // s^2 = trace + 2*sqrt(det)
        const s2 = tr.add(sqrtDet.mul(2));
        const s = s2.sqrt();

        // sqrt(A) = (A + sqrt(det)*I) / s
        const sqrtDetI = Matrix2x2.scalar(sqrtDet);
        const numerator = this.add(sqrtDetI);

        // Divide by s
        const invS = new Complex(1, 0).div(s);
        return numerator.mul(invS);
    }

    /**
     * Matrix exponential using Pad approximation for 2x2
     * For ABCD/S-parameter calculations, we typically use sinh/cosh instead
     */
    exp() {
        // For 2x2 matrices, exp(A) = e^(tr/2) * [cosh(q)*I + sinh(q)/q * (A - tr/2*I)]
        // where q = sqrt((tr^2/4 - det))
        const tr = this.trace();
        const det = this.det();
        const halfTr = tr.mul(0.5);

        // q^2 = tr^2/4 - det
        const q2 = halfTr.mul(halfTr).sub(det);
        const q = q2.sqrt();

        // exp(tr/2)
        const expHalfTr = halfTr.exp();

        // Handle case where q is very small
        const qAbs = q.abs();
        let coshQ, sinhQoverQ;

        if (qAbs < 1e-10) {
            // Taylor expansion for small q
            coshQ = new Complex(1, 0);
            sinhQoverQ = new Complex(1, 0);
        } else {
            coshQ = q.cosh();
            sinhQoverQ = q.sinh().div(q);
        }

        // A - (tr/2)*I
        const Ashift = this.sub(Matrix2x2.scalar(halfTr));

        // Result = exp(tr/2) * [cosh(q)*I + (sinh(q)/q)*(A - tr/2*I)]
        const term1 = Matrix2x2.identity().mul(coshQ);
        const term2 = Ashift.mul(sinhQoverQ);

        return term1.add(term2).mul(expHalfTr);
    }

    /**
     * Matrix sinh for transmission line calculations
     * sinh(A) = (exp(A) - exp(-A)) / 2
     */
    sinh() {
        const expA = this.exp();
        const negA = this.mul(-1);
        const expNegA = negA.exp();
        return expA.sub(expNegA).mul(0.5);
    }

    /**
     * Matrix cosh for transmission line calculations
     * cosh(A) = (exp(A) + exp(-A)) / 2
     */
    cosh() {
        const expA = this.exp();
        const negA = this.mul(-1);
        const expNegA = negA.exp();
        return expA.add(expNegA).mul(0.5);
    }

    toString() {
        return `[[${this.a}, ${this.b}], [${this.c}, ${this.d}]]`;
    }
}

// --- Real 2×2 helpers (matrices as [[a,b],[c,d]] of plain numbers) ----------
const mat2Mul = (A, B) => [
    [A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1]],
    [A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]],
];
const mat2Inv = (M) => { const d = M[0][0] * M[1][1] - M[0][1] * M[1][0]; return [[M[1][1] / d, -M[0][1] / d], [-M[1][0] / d, M[0][0] / d]]; };
const mat2T = (M) => [[M[0][0], M[1][0]], [M[0][1], M[1][1]]];

// Eigenpairs of a real 2×2 M (here M = [C0]⁻¹[C], real eigenvalues since C, C0 are SPD).
// Eigenvectors are normalised to |v|=√2 — equals the ±1 odd/even drive for a symmetric pair
// (so symmetric results are reproduced) AND is continuous through symmetry (normalising to the
// max component would jump when the dominant component switches, making the modal Z0
// discontinuous in the asymmetry). Returns { vals:[λ1,λ2], vecs:[v1,v2] }.
function eig2x2(M) {
    const a = M[0][0], b = M[0][1], c = M[1][0], d = M[1][1];
    const tr = a + d, disc = Math.sqrt(Math.max(0, tr * tr - 4 * (a * d - b * c)));
    const vecFor = (l) => {
        const v = Math.abs(b) > 1e-300 ? [b, l - a] : (Math.abs(c) > 1e-300 ? [l - d, c] : [1, 0]);
        const n = Math.hypot(v[0], v[1]);
        return n > 0 ? [v[0] * Math.SQRT2 / n, v[1] * Math.SQRT2 / n] : v;
    };
    return { vals: [(tr + disc) / 2, (tr - disc) / 2], vecs: [vecFor((tr + disc) / 2), vecFor((tr - disc) / 2)] };
}

export { Matrix2x2, mat2Mul, mat2Inv, mat2T, eig2x2 };
