#include <emscripten.h>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include <Eigen/Dense>
#include <Spectra/SymEigsSolver.h>
#include <complex>
#include <vector>
#include <cmath>

using Complex = std::complex<double>;
using SpMat = Eigen::SparseMatrix<Complex>;
using RSpMat = Eigen::SparseMatrix<double>;
using CVector = Eigen::VectorXcd;
using RVector = Eigen::VectorXd;
using CMatrix = Eigen::MatrixXcd;
using RMatrix = Eigen::MatrixXd;

// Build complex sparse matrix from CSR arrays
static SpMat buildSparseMatrix(int N, int nnz,
    const int* rowPtr, const int* colIdx,
    const double* vals_re, const double* vals_im) {
    std::vector<Eigen::Triplet<Complex>> triplets;
    triplets.reserve(nnz);
    for (int i = 0; i < N; i++) {
        for (int j = rowPtr[i]; j < rowPtr[i + 1]; j++) {
            triplets.emplace_back(i, colIdx[j], Complex(vals_re[j], vals_im[j]));
        }
    }
    SpMat mat(N, N);
    mat.setFromTriplets(triplets.begin(), triplets.end());
    mat.makeCompressed();
    return mat;
}

// Build real sparse matrix from CSR arrays (ignores imaginary part)
static RSpMat buildRealSparseMatrix(int N, int nnz,
    const int* rowPtr, const int* colIdx, const double* vals_re) {
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(nnz);
    for (int i = 0; i < N; i++) {
        for (int j = rowPtr[i]; j < rowPtr[i + 1]; j++) {
            triplets.emplace_back(i, colIdx[j], vals_re[j]);
        }
    }
    RSpMat mat(N, N);
    mat.setFromTriplets(triplets.begin(), triplets.end());
    mat.makeCompressed();
    return mat;
}

// Combined shift-invert-B operator: y = (A-σB)⁻¹ B x
// This is a standard (non-generalized) operator for SymEigsSolver.
// Eigenvalues of this operator are ν = 1/(λ-σ), so λ = σ + 1/ν.
class ShiftInvertBOp {
public:
    using Scalar = double;
private:
    const RSpMat& m_B;
    int m_n;
    Eigen::SimplicialLDLT<RSpMat> m_solver;
    bool m_factored = false;
public:
    ShiftInvertBOp(const RSpMat& A, const RSpMat& B, double sigma)
        : m_B(B), m_n(A.rows())
    {
        RSpMat C = A - sigma * B;
        m_solver.compute(C);
        m_factored = (m_solver.info() == Eigen::Success);
    }
    bool factored() const { return m_factored; }
    Eigen::Index rows() const { return m_n; }
    Eigen::Index cols() const { return m_n; }
    // y = (A - σB)⁻¹ B x
    void perform_op(const double* x_in, double* y_out) const {
        Eigen::Map<const RVector> x(x_in, m_n);
        Eigen::Map<RVector> y(y_out, m_n);
        RVector Bx = m_B * x;
        y.noalias() = m_solver.solve(Bx);
    }
};

// Real-valued shift-invert eigensolver using Spectra (implicit restart Lanczos).
static int solve_real(
    int N,
    int nnz_a,
    int* a_rowPtr, int* a_colIdx, double* a_vals_re,
    int nnz_b,
    int* b_rowPtr, int* b_colIdx, double* b_vals_re,
    double sigma_re,
    int num_eigenvalues, int krylov_size,
    double* out_evals_re, double* out_evals_im,
    double* out_evecs_re, double* out_evecs_im,
    double* init_re
) {
    int n = N;
    int nev = num_eigenvalues;
    int ncv = krylov_size;

    if (ncv <= nev) ncv = 2 * nev + 1;
    if (ncv > n) ncv = n;

    RSpMat A = buildRealSparseMatrix(n, nnz_a, a_rowPtr, a_colIdx, a_vals_re);
    RSpMat B = buildRealSparseMatrix(n, nnz_b, b_rowPtr, b_colIdx, b_vals_re);

    ShiftInvertBOp op(A, B, sigma_re);
    if (!op.factored()) return -1;

    Spectra::SymEigsSolver<ShiftInvertBOp> eigs(op, nev, ncv);

    if (init_re != nullptr) {
        double norm2 = 0;
        for (int i = 0; i < n; i++) norm2 += init_re[i] * init_re[i];
        if (norm2 > 1e-60) {
            eigs.init(init_re);
        } else {
            eigs.init();
        }
    } else {
        eigs.init();
    }

    int nconv = eigs.compute(Spectra::SortRule::LargestMagn, 1000, 1e-10);

    if (eigs.info() != Spectra::CompInfo::Successful &&
        eigs.info() != Spectra::CompInfo::NotConverging) {
        return -1;
    }

    // Extract eigenvalues and transform back: λ = σ + 1/ν
    RVector nu = eigs.eigenvalues();
    RMatrix evecs = eigs.eigenvectors();

    int nout = std::min(nconv, nev);
    for (int k = 0; k < nout; k++) {
        double lambda = sigma_re + 1.0 / nu(k);
        out_evals_re[k] = lambda;
        out_evals_im[k] = 0.0;
        for (int j = 0; j < n; j++) {
            out_evecs_re[k * n + j] = evecs(j, k);
            out_evecs_im[k * n + j] = 0.0;
        }
    }

    return nout;
}

// Complex-valued shift-invert Arnoldi for generalized Ax = λBx.
// Uses complex Arnoldi with (A-σB)⁻¹B operator and NaN-safe Hessenberg decomposition.
static int solve_complex(
    int N,
    int nnz_a,
    int* a_rowPtr, int* a_colIdx, double* a_vals_re, double* a_vals_im,
    int nnz_b,
    int* b_rowPtr, int* b_colIdx, double* b_vals_re, double* b_vals_im,
    double sigma_re, double sigma_im,
    int num_eigenvalues, int krylov_size,
    double* out_evals_re, double* out_evals_im,
    double* out_evecs_re, double* out_evecs_im,
    double* init_re, double* init_im
) {
    int n = N;
    int nev = num_eigenvalues;
    int m = krylov_size;

    if (m <= nev) m = 2 * nev + 1;
    if (m > n) m = n;

    Complex sigma(sigma_re, sigma_im);

    SpMat A = buildSparseMatrix(n, nnz_a, a_rowPtr, a_colIdx, a_vals_re, a_vals_im);
    SpMat B = buildSparseMatrix(n, nnz_b, b_rowPtr, b_colIdx, b_vals_re, b_vals_im);

    // Factor C = A - σB
    SpMat C = A - sigma * B;
    Eigen::SparseLU<SpMat> solver;
    solver.compute(C);
    if (solver.info() != Eigen::Success) {
        return -1; // Factorization failed
    }

    // Arnoldi iteration with modified Gram-Schmidt (double orthogonalization)
    CMatrix V(n, m + 1);
    CMatrix H = CMatrix::Zero(m + 1, m);

    // Initial vector
    CVector v0(n);
    if (init_re != nullptr) {
        for (int i = 0; i < n; i++) {
            double re = init_re[i];
            double im = (init_im != nullptr) ? init_im[i] : 0.0;
            v0(i) = Complex(re, im);
        }
        if (v0.norm() < 1e-30) {
            for (int i = 0; i < n; i++)
                v0(i) = Complex(1.0 / std::sqrt((double)n), 0.0);
        }
    } else {
        for (int i = 0; i < n; i++) {
            v0(i) = Complex(1.0 / std::sqrt((double)n), 0.0);
        }
    }
    V.col(0) = v0 / v0.norm();

    int actual_m = m;

    for (int j = 0; j < m; j++) {
        CVector Bv = B * V.col(j);
        CVector w = solver.solve(Bv);

        // Check for NaN/Inf from solver
        if (!w.allFinite()) {
            actual_m = j;
            break;
        }

        // Modified Gram-Schmidt with double orthogonalization
        for (int i = 0; i <= j; i++) {
            Complex h = V.col(i).dot(w);
            H(i, j) = h;
            w -= h * V.col(i);
        }
        for (int i = 0; i <= j; i++) {
            Complex h = V.col(i).dot(w);
            H(i, j) += h;
            w -= h * V.col(i);
        }

        double wnorm = w.norm();
        H(j + 1, j) = Complex(wnorm, 0.0);

        if (wnorm < 1e-14 || !std::isfinite(wnorm)) {
            actual_m = j + 1;
            break;
        }

        if (j + 1 < m) {
            V.col(j + 1) = w / wnorm;
        }
    }

    if (actual_m < 1) return 0;

    CMatrix Hm = H.topLeftCorner(actual_m, actual_m);

    // Check Hessenberg matrix for NaN before eigendecomposition
    if (!Hm.allFinite()) return -2;

    Eigen::ComplexEigenSolver<CMatrix> eig_solver;
    eig_solver.compute(Hm);

    if (eig_solver.info() != Eigen::Success) {
        return -2;
    }

    CVector theta = eig_solver.eigenvalues();
    CMatrix Y = eig_solver.eigenvectors();

    std::vector<int> indices(actual_m);
    for (int i = 0; i < actual_m; i++) indices[i] = i;

    std::sort(indices.begin(), indices.end(), [&](int a, int b) {
        return std::abs(theta(a)) > std::abs(theta(b));
    });

    int nout = std::min(nev, actual_m);

    int nconv = 0;
    for (int k = 0; k < nout; k++) {
        int idx = indices[k];
        Complex th = theta(idx);

        if (std::abs(th) < 1e-15) continue;

        Complex lambda = sigma + 1.0 / th;

        CVector x = V.leftCols(actual_m) * Y.col(idx);

        CVector Ax = A * x;
        CVector Bx = B * x;
        CVector residual = Ax - lambda * Bx;
        double rel_res = residual.norm() / (x.norm() * (Ax.norm() + std::abs(lambda) * Bx.norm()));

        if (rel_res < 1e-4 || k < nev) {
            out_evals_re[nconv] = lambda.real();
            out_evals_im[nconv] = lambda.imag();
            for (int j = 0; j < n; j++) {
                out_evecs_re[nconv * n + j] = x(j).real();
                out_evecs_im[nconv * n + j] = x(j).imag();
            }
            nconv++;
        }
    }

    return nconv;
}

// Dispatch: use real path when sigma_im == 0 and matrices have no imaginary parts
static bool hasImaginaryParts(int nnz, const double* vals_im) {
    for (int i = 0; i < nnz; i++) {
        if (vals_im[i] != 0.0) return true;
    }
    return false;
}

static int solve_core(
    int N,
    int nnz_a,
    int* a_rowPtr, int* a_colIdx, double* a_vals_re, double* a_vals_im,
    int nnz_b,
    int* b_rowPtr, int* b_colIdx, double* b_vals_re, double* b_vals_im,
    double sigma_re, double sigma_im,
    int num_eigenvalues, int krylov_size,
    double* out_evals_re, double* out_evals_im,
    double* out_evecs_re, double* out_evecs_im,
    double* init_re, double* init_im
) {
    // Use fast real path when shift is real and matrices are real
    if (sigma_im == 0.0 &&
        !hasImaginaryParts(nnz_a, a_vals_im) &&
        !hasImaginaryParts(nnz_b, b_vals_im)) {
        return solve_real(N, nnz_a, a_rowPtr, a_colIdx, a_vals_re,
                          nnz_b, b_rowPtr, b_colIdx, b_vals_re,
                          sigma_re, num_eigenvalues, krylov_size,
                          out_evals_re, out_evals_im, out_evecs_re, out_evecs_im,
                          init_re);
    }

    return solve_complex(N, nnz_a, a_rowPtr, a_colIdx, a_vals_re, a_vals_im,
                         nnz_b, b_rowPtr, b_colIdx, b_vals_re, b_vals_im,
                         sigma_re, sigma_im, num_eigenvalues, krylov_size,
                         out_evals_re, out_evals_im, out_evecs_re, out_evecs_im,
                         init_re, init_im);
}

// Public API
extern "C" {

// Rayleigh quotient iteration for complex generalized eigenvalue problem.
// Refines a known eigenvalue (sigma) of Ax = λBx where A,B are complex sparse.
// Uses inverse iteration: w = (A-σB)⁻¹Bx, then Rayleigh quotient update.
// Returns 1 on success (1 eigenvalue), 0 or negative on failure.
// Much simpler than full Arnoldi — no Krylov subspace, no spurious modes.
EMSCRIPTEN_KEEPALIVE
int solve_rqi(
    int N,
    int nnz_a,
    int* a_rowPtr, int* a_colIdx, double* a_vals_re, double* a_vals_im,
    int nnz_b,
    int* b_rowPtr, int* b_colIdx, double* b_vals_re, double* b_vals_im,
    double sigma_re, double sigma_im,
    int max_iter,
    double* init_re, double* init_im,
    double* out_eval_re, double* out_eval_im,
    double* out_evec_re, double* out_evec_im
) {
    int n = N;
    SpMat A = buildSparseMatrix(n, nnz_a, a_rowPtr, a_colIdx, a_vals_re, a_vals_im);
    SpMat B = buildSparseMatrix(n, nnz_b, b_rowPtr, b_colIdx, b_vals_re, b_vals_im);

    Complex mu(sigma_re, sigma_im);

    // Initial vector
    CVector x(n);
    for (int i = 0; i < n; i++) {
        double re = init_re[i];
        double im = (init_im != nullptr) ? init_im[i] : 0.0;
        x(i) = Complex(re, im);
    }
    double xnorm = x.norm();
    if (xnorm < 1e-30) return -3; // Zero initial vector
    x /= xnorm;

    for (int iter = 0; iter < max_iter; iter++) {
        // Factor (A - μB)
        SpMat C = A - mu * B;
        Eigen::SparseLU<SpMat> solver;
        solver.compute(C);
        if (solver.info() != Eigen::Success) return -1;

        // Inverse iteration: w = (A - μB)⁻¹ B x
        CVector Bx = B * x;
        CVector w = solver.solve(Bx);
        if (!w.allFinite()) return -2;

        // Normalize
        double wnorm = w.norm();
        if (wnorm < 1e-30) return -2;
        x = w / wnorm;

        // Rayleigh quotient: μ = x^H A x / (x^H B x)
        CVector Ax = A * x;
        Bx = B * x;
        Complex num = x.conjugate().dot(Ax);
        Complex den = x.conjugate().dot(Bx);
        if (std::abs(den) < 1e-30) return -2;
        mu = num / den;
    }

    // Output
    *out_eval_re = mu.real();
    *out_eval_im = mu.imag();
    for (int i = 0; i < n; i++) {
        out_evec_re[i] = x(i).real();
        out_evec_im[i] = x(i).imag();
    }
    return 1;
}

// Quadratic eigenvalue problem RQI: refine eigenvalue of (A0 + λA1 + λ²A2)x = 0.
// Uses inverse iteration with Q(σ) = A0+σA1+σ²A2, then quadratic Rayleigh quotient.
// Returns 1 on success, negative on failure.
EMSCRIPTEN_KEEPALIVE
int solve_qep_rqi(
    int N,
    int nnz0, int* r0, int* c0, double* v0re, double* v0im,
    int nnz1, int* r1, int* c1, double* v1re, double* v1im,
    int nnz2, int* r2, int* c2, double* v2re, double* v2im,
    double sigma_re, double sigma_im,
    int max_iter,
    double* init_re, double* init_im,
    double* out_eval_re, double* out_eval_im,
    double* out_evec_re, double* out_evec_im
) {
    int n = N;
    SpMat A0 = buildSparseMatrix(n, nnz0, r0, c0, v0re, v0im);
    SpMat A1 = buildSparseMatrix(n, nnz1, r1, c1, v1re, v1im);
    SpMat A2 = buildSparseMatrix(n, nnz2, r2, c2, v2re, v2im);

    Complex mu(sigma_re, sigma_im);

    CVector x(n);
    for (int i = 0; i < n; i++) {
        double re = init_re[i];
        double im = (init_im != nullptr) ? init_im[i] : 0.0;
        x(i) = Complex(re, im);
    }
    double xnorm = x.norm();
    if (xnorm < 1e-30) return -3;
    x /= xnorm;

    for (int iter = 0; iter < max_iter; iter++) {
        // Q(mu) = A0 + mu*A1 + mu^2*A2
        SpMat Q = A0 + mu * A1 + (mu * mu) * A2;
        Eigen::SparseLU<SpMat> solver;
        solver.compute(Q);
        if (solver.info() != Eigen::Success) return -1;

        // Inverse iteration: solve Q(mu)*w = A2*x
        CVector rhs = A2 * x;
        CVector w = solver.solve(rhs);
        if (!w.allFinite()) return -2;

        double wnorm = w.norm();
        if (wnorm < 1e-30) return -2;
        x = w / wnorm;

        // Quadratic Rayleigh quotient: find mu such that x^H * Q(mu) * x = 0
        // c0 + c1*mu + c2*mu^2 = 0
        Complex qc0 = x.conjugate().dot(A0 * x);
        Complex qc1 = x.conjugate().dot(A1 * x);
        Complex qc2 = x.conjugate().dot(A2 * x);

        if (std::abs(qc2) < 1e-30) {
            // Linear: mu = -qc0/qc1
            if (std::abs(qc1) < 1e-30) return -2;
            mu = -qc0 / qc1;
        } else {
            Complex disc = qc1 * qc1 - 4.0 * qc0 * qc2;
            Complex sq = std::sqrt(disc);
            Complex mu1 = (-qc1 + sq) / (2.0 * qc2);
            Complex mu2 = (-qc1 - sq) / (2.0 * qc2);
            mu = (std::abs(mu1 - mu) < std::abs(mu2 - mu)) ? mu1 : mu2;
        }
    }

    *out_eval_re = mu.real();
    *out_eval_im = mu.imag();
    for (int i = 0; i < n; i++) {
        out_evec_re[i] = x(i).real();
        out_evec_im[i] = x(i).imag();
    }
    return 1;
}

EMSCRIPTEN_KEEPALIVE
int solve_generalized_eigen(
    int N,
    int nnz_a,
    int* a_rowPtr, int* a_colIdx, double* a_vals_re, double* a_vals_im,
    int nnz_b,
    int* b_rowPtr, int* b_colIdx, double* b_vals_re, double* b_vals_im,
    double sigma_re, double sigma_im,
    int num_eigenvalues, int krylov_size,
    double* out_evals_re, double* out_evals_im,
    double* out_evecs_re, double* out_evecs_im
) {
    return solve_core(N, nnz_a, a_rowPtr, a_colIdx, a_vals_re, a_vals_im,
                      nnz_b, b_rowPtr, b_colIdx, b_vals_re, b_vals_im,
                      sigma_re, sigma_im, num_eigenvalues, krylov_size,
                      out_evals_re, out_evals_im, out_evecs_re, out_evecs_im,
                      nullptr, nullptr);
}

EMSCRIPTEN_KEEPALIVE
int solve_generalized_eigen_with_init(
    int N,
    int nnz_a,
    int* a_rowPtr, int* a_colIdx, double* a_vals_re, double* a_vals_im,
    int nnz_b,
    int* b_rowPtr, int* b_colIdx, double* b_vals_re, double* b_vals_im,
    double sigma_re, double sigma_im,
    int num_eigenvalues, int krylov_size,
    double* out_evals_re, double* out_evals_im,
    double* out_evecs_re, double* out_evecs_im,
    double* init_re, double* init_im
) {
    return solve_core(N, nnz_a, a_rowPtr, a_colIdx, a_vals_re, a_vals_im,
                      nnz_b, b_rowPtr, b_colIdx, b_vals_re, b_vals_im,
                      sigma_re, sigma_im, num_eigenvalues, krylov_size,
                      out_evals_re, out_evals_im, out_evecs_re, out_evecs_im,
                      init_re, init_im);
}

// Direct sparse solver: factorize once, solve for nRhs right-hand sides.
// Uses SimplicialLDLT for SPD matrices, falls back to SparseLU.
// rhs_arr and x_arr are nRhs×N column-major (each RHS is contiguous N doubles).
// Returns 0 on success, negative on error.
EMSCRIPTEN_KEEPALIVE
int solve_sparse_multi(
    int N, int nnz,
    int* rowPtr, int* colIdx, double* values,
    int nRhs, double* rhs_arr, double* x_arr
) {
    RSpMat A = buildRealSparseMatrix(N, nnz, rowPtr, colIdx, values);

    // SimplicialLDLT reads only the lower triangle: for a NONSYMMETRIC input it
    // silently solves the wrong (symmetrized) system while reporting Success.
    // Only take the fast path when the matrix is actually symmetric.
    bool symmetric = true;
    {
        double scale = 0;
        for (int k = 0; k < A.nonZeros(); k++) scale = std::max(scale, std::abs(A.valuePtr()[k]));
        RSpMat D = RSpMat(A - RSpMat(A.transpose()));
        for (int k = 0; k < D.nonZeros(); k++) {
            if (std::abs(D.valuePtr()[k]) > 1e-12 * scale) { symmetric = false; break; }
        }
    }

    // Try LDLT first (fast for SPD/symmetric indefinite)
    Eigen::SimplicialLDLT<RSpMat> ldlt;
    if (symmetric) ldlt.compute(A);
    if (symmetric && ldlt.info() == Eigen::Success) {
        for (int r = 0; r < nRhs; r++) {
            Eigen::Map<RVector> b(rhs_arr + r * N, N);
            Eigen::Map<RVector> x(x_arr + r * N, N);
            x = ldlt.solve(b);
        }
        return 0;
    }

    // Fallback to SparseLU
    Eigen::SparseLU<RSpMat> lu;
    lu.compute(A);
    if (lu.info() != Eigen::Success) return -1;
    for (int r = 0; r < nRhs; r++) {
        Eigen::Map<RVector> b(rhs_arr + r * N, N);
        Eigen::Map<RVector> x(x_arr + r * N, N);
        x = lu.solve(b);
    }
    return 0;
}

} // extern "C"
