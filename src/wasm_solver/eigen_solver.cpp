#include <emscripten.h>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include <Eigen/Dense>
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

    // Walk ALL Ritz pairs in closest-to-shift order and keep the converged ones,
    // up to nev — a rejected (unconverged) pair's slot can be filled by a
    // converged pair further from the shift.
    int nconv = 0;
    for (int k = 0; k < actual_m && nconv < nev; k++) {
        int idx = indices[k];
        Complex th = theta(idx);

        if (std::abs(th) < 1e-15) continue;

        Complex lambda = sigma + 1.0 / th;

        CVector x = V.leftCols(actual_m) * Y.col(idx);

        CVector Ax = A * x;
        CVector Bx = B * x;
        CVector residual = Ax - lambda * Bx;
        double rel_res = residual.norm() / (x.norm() * (Ax.norm() + std::abs(lambda) * Bx.norm()));

        // Accept only Ritz pairs that actually converged within the Krylov space.
        // (This gate used to be `rel_res < 1e-4 || k < nev`, and since the loop is
        // bounded by nout = min(nev, actual_m) the `k < nev` arm was always true —
        // unconverged junk pairs were returned as converged modes.)
        if (rel_res < 1e-4) {
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

// Real-arithmetic shift-invert Arnoldi for generalized Ax = λBx with REAL A, B, σ.
// This is the SAME algorithm as solve_complex — on real data with a real start
// vector the complex Arnoldi stays exactly in the real subspace, so this twin
// produces the same Krylov space and Ritz pairs while paying real-arithmetic
// costs (the sparse LU alone is ~4× cheaper). Unlike the earlier removed "fast
// path" (Spectra SymEigsSolver = Lanczos, which wrongly assumed the operator is
// symmetric in the standard inner product — it is only B-self-adjoint, and B is
// indefinite), plain Arnoldi makes NO symmetry assumption: complex-pair Ritz
// values from the real Hessenberg are handled, and every returned pair passes
// the explicit residual gate against the true pencil.
static int solve_real(
    int N,
    int nnz_a, int* a_rowPtr, int* a_colIdx, double* a_vals_re,
    int nnz_b, int* b_rowPtr, int* b_colIdx, double* b_vals_re,
    double sigma_re,
    int num_eigenvalues, int krylov_size,
    double* out_evals_re, double* out_evals_im,
    double* out_evecs_re, double* out_evecs_im,
    double* init_re
) {
    int n = N;
    int nev = num_eigenvalues;
    int m = krylov_size;
    if (m <= nev) m = 2 * nev + 1;
    if (m > n) m = n;

    RSpMat A = buildRealSparseMatrix(n, nnz_a, a_rowPtr, a_colIdx, a_vals_re);
    RSpMat B = buildRealSparseMatrix(n, nnz_b, b_rowPtr, b_colIdx, b_vals_re);

    // Factor C = A - σB (real general LU — pivoted, safe for the indefinite pencil)
    RSpMat C = A - sigma_re * B;
    Eigen::SparseLU<RSpMat> solver;
    solver.compute(C);
    if (solver.info() != Eigen::Success) {
        return -1;
    }

    RMatrix V(n, m + 1);
    RMatrix H = RMatrix::Zero(m + 1, m);

    RVector v0(n);
    if (init_re != nullptr) {
        for (int i = 0; i < n; i++) v0(i) = init_re[i];
        if (v0.norm() < 1e-30) v0.setConstant(1.0 / std::sqrt((double)n));
    } else {
        v0.setConstant(1.0 / std::sqrt((double)n));
    }
    V.col(0) = v0 / v0.norm();

    int actual_m = m;
    for (int j = 0; j < m; j++) {
        RVector Bv = B * V.col(j);
        RVector w = solver.solve(Bv);
        if (!w.allFinite()) { actual_m = j; break; }
        // Modified Gram-Schmidt with double orthogonalization
        for (int i = 0; i <= j; i++) {
            double h = V.col(i).dot(w);
            H(i, j) = h;
            w -= h * V.col(i);
        }
        for (int i = 0; i <= j; i++) {
            double h = V.col(i).dot(w);
            H(i, j) += h;
            w -= h * V.col(i);
        }
        double wnorm = w.norm();
        H(j + 1, j) = wnorm;
        if (wnorm < 1e-14 || !std::isfinite(wnorm)) { actual_m = j + 1; break; }
        if (j + 1 < m) V.col(j + 1) = w / wnorm;
    }

    if (actual_m < 1) return 0;

    RMatrix Hm = H.topLeftCorner(actual_m, actual_m);
    if (!Hm.allFinite()) return -2;

    // Real Schur eigendecomposition — Ritz values may come in complex-conjugate
    // pairs (the real indefinite pencil admits them); they are handled like any
    // other candidate and must pass the residual gate below.
    Eigen::EigenSolver<RMatrix> eig_solver;
    eig_solver.compute(Hm);
    if (eig_solver.info() != Eigen::Success) return -2;

    CVector theta = eig_solver.eigenvalues();
    CMatrix Y = eig_solver.eigenvectors();

    std::vector<int> indices(actual_m);
    for (int i = 0; i < actual_m; i++) indices[i] = i;
    std::sort(indices.begin(), indices.end(), [&](int a, int b) {
        return std::abs(theta(a)) > std::abs(theta(b));
    });

    int nconv = 0;
    for (int k = 0; k < actual_m && nconv < nev; k++) {
        int idx = indices[k];
        Complex th = theta(idx);
        if (std::abs(th) < 1e-15) continue;
        Complex lambda = sigma_re + 1.0 / th;

        // Ritz vector: complex combination of the REAL basis
        RVector xr = V.leftCols(actual_m) * Y.col(idx).real();
        RVector xi = V.leftCols(actual_m) * Y.col(idx).imag();

        // Residual ‖Ax − λBx‖ with real A, B and complex λ, x (componentwise)
        RVector Axr = A * xr, Axi = A * xi;
        RVector Bxr = B * xr, Bxi = B * xi;
        double lr = lambda.real(), li = lambda.imag();
        double res2 = 0, xn2 = 0, Axn2 = 0, Bxn2 = 0;
        for (int i2 = 0; i2 < n; i2++) {
            const double rr = Axr(i2) - (lr * Bxr(i2) - li * Bxi(i2));
            const double ri = Axi(i2) - (lr * Bxi(i2) + li * Bxr(i2));
            res2 += rr * rr + ri * ri;
            xn2 += xr(i2) * xr(i2) + xi(i2) * xi(i2);
            Axn2 += Axr(i2) * Axr(i2) + Axi(i2) * Axi(i2);
            Bxn2 += Bxr(i2) * Bxr(i2) + Bxi(i2) * Bxi(i2);
        }
        double rel_res = std::sqrt(res2) /
            (std::sqrt(xn2) * (std::sqrt(Axn2) + std::abs(lambda) * std::sqrt(Bxn2)));

        if (rel_res < 1e-4) {
            out_evals_re[nconv] = lambda.real();
            out_evals_im[nconv] = lambda.imag();
            for (int j = 0; j < n; j++) {
                out_evecs_re[nconv * n + j] = xr(j);
                out_evecs_im[nconv * n + j] = xi(j);
            }
            nconv++;
        }
    }

    return nconv;
}

static bool allZero(int nnz, const double* v) {
    if (v == nullptr) return true;
    for (int i = 0; i < nnz; i++) if (v[i] != 0.0) return false;
    return true;
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
    // Fully real problem (closed PEC/PMC boundaries, lossless ε, real shift and
    // seed): run the real-arithmetic Arnoldi twin — the same iteration the
    // complex path would perform (a real start keeps the complex Arnoldi exactly
    // real), at real-arithmetic cost. Anything imaginary anywhere (radiating-ABC
    // Robin term, lossy ε, complex shift/seed) goes through the complex path.
    // If the real path fails outright (factorization/decomposition failure or
    // zero converged pairs) the complex path is retried as a safety net before
    // reporting failure.
    const bool realProblem = sigma_im == 0.0 &&
        allZero(nnz_a, a_vals_im) && allZero(nnz_b, b_vals_im) &&
        (init_re == nullptr || init_im == nullptr || allZero(N, init_im));
    if (realProblem) {
        int nconv = solve_real(N, nnz_a, a_rowPtr, a_colIdx, a_vals_re,
                               nnz_b, b_rowPtr, b_colIdx, b_vals_re,
                               sigma_re, num_eigenvalues, krylov_size,
                               out_evals_re, out_evals_im, out_evecs_re, out_evecs_im,
                               init_re);
        if (nconv > 0) return nconv;
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
  try {
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
  } catch (...) {
    return -9;
  }
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
  try {
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
  } catch (...) {
    return -9;
  }
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
    // Catch everything (Spectra std::invalid_argument, Eigen std::bad_alloc, …):
    // an uncaught exception would abort() and permanently poison the module.
    try {
        return solve_core(N, nnz_a, a_rowPtr, a_colIdx, a_vals_re, a_vals_im,
                          nnz_b, b_rowPtr, b_colIdx, b_vals_re, b_vals_im,
                          sigma_re, sigma_im, num_eigenvalues, krylov_size,
                          out_evals_re, out_evals_im, out_evecs_re, out_evecs_im,
                          nullptr, nullptr);
    } catch (...) {
        return -9;
    }
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
    try {
        return solve_core(N, nnz_a, a_rowPtr, a_colIdx, a_vals_re, a_vals_im,
                          nnz_b, b_rowPtr, b_colIdx, b_vals_re, b_vals_im,
                          sigma_re, sigma_im, num_eigenvalues, krylov_size,
                          out_evals_re, out_evals_im, out_evecs_re, out_evecs_im,
                          init_re, init_im);
    } catch (...) {
        return -9;
    }
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
  try {
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
  } catch (...) {
    return -9;
  }
}

} // extern "C"
