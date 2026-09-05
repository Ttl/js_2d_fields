#include <emscripten.h>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include <Eigen/Dense>
#include <complex>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>

using Complex = std::complex<double>;
using SpMat = Eigen::SparseMatrix<Complex>;
using RSpMat = Eigen::SparseMatrix<double>;
using CVector = Eigen::VectorXcd;
using RVector = Eigen::VectorXd;
using CMatrix = Eigen::MatrixXcd;
using RMatrix = Eigen::MatrixXd;

// Complex-symmetric scalar for Eigen's simplicial LDL^T
//
// The MQS eddy-current system K = S + j*β*M (S the real P2 stiffness, M the
// conductor mass, β = ωμ_0σ) is complex symmetric, K^T = K, not Hermitian.
// Eigen's SimplicialLDLT computes LDL^H, the wrong factorization for such a
// matrix. Its only complex-aware operations, though, are numext::conj() and
// numext::real(), and for a scalar type whose NumTraits says IsComplex = 0 both
// are the identity: the same code then computes LDL^T with a complex diagonal
// D, which is exactly the complex-symmetric factorization. `csym` is
// std::complex<double> wearing that disguise.
//
// No pivoting is needed: S is SPD (every MQS domain wall is Dirichlet) and βM is
// positive semidefinite, so K has a positive-definite Hermitian part, and so does
// every Schur complement, i.e. every pivot has Re(d) > 0 (the Higham 1998 class;
// measured min Re(d) ≈ 1 with a |d| spread of ~1e-3 on the reference systems).
// solve_complex_symmetric still verifies the residual and falls back to the
// pivoted complex SparseLU if anything is off.
struct csym {
    Complex v;
    csym() : v(0.0) {}
    csym(double r) : v(r) {}
    csym(int r) : v((double)r) {}
    csym(double r, double i) : v(r, i) {}
    csym(const Complex& c) : v(c) {}
    csym operator+(const csym& o) const { return csym(v + o.v); }
    csym operator-(const csym& o) const { return csym(v - o.v); }
    csym operator*(const csym& o) const { return csym(v * o.v); }
    csym operator/(const csym& o) const { return csym(v / o.v); }
    csym operator-() const { return csym(-v); }
    csym& operator+=(const csym& o) { v += o.v; return *this; }
    csym& operator-=(const csym& o) { v -= o.v; return *this; }
    csym& operator*=(const csym& o) { v *= o.v; return *this; }
    csym& operator/=(const csym& o) { v /= o.v; return *this; }
    bool operator==(const csym& o) const { return v == o.v; }
    bool operator!=(const csym& o) const { return v != o.v; }
    // Orderings exist only because Eigen's "real" scalar concept requires them;
    // LDLᵀ never orders values. Compare magnitudes so they are at least sane.
    bool operator<(const csym& o) const { return std::abs(v) < std::abs(o.v); }
    bool operator>(const csym& o) const { return std::abs(v) > std::abs(o.v); }
    bool operator<=(const csym& o) const { return std::abs(v) <= std::abs(o.v); }
    bool operator>=(const csym& o) const { return std::abs(v) >= std::abs(o.v); }
};
static inline csym operator*(double a, const csym& b) { return csym(a * b.v); }
static inline csym abs(const csym& a) { return csym(std::abs(a.v)); }
static inline csym sqrt(const csym& a) { return csym(std::sqrt(a.v)); }
namespace Eigen {
template<> struct NumTraits<csym> {
    typedef csym Real; typedef csym NonInteger; typedef csym Nested; typedef csym Literal;
    enum { IsComplex = 0, IsInteger = 0, IsSigned = 1, RequireInitialization = 0,
           ReadCost = 2, AddCost = 2, MulCost = 6 };
    static inline csym epsilon() { return csym(std::numeric_limits<double>::epsilon()); }
    static inline csym dummy_precision() { return csym(1e-12); }
    static inline csym highest() { return csym(std::numeric_limits<double>::max()); }
    static inline csym lowest() { return csym(-std::numeric_limits<double>::max()); }
    static inline int digits10() { return std::numeric_limits<double>::digits10; }
    static inline int digits() { return std::numeric_limits<double>::digits; }
    static inline int max_digits10() { return std::numeric_limits<double>::max_digits10; }
    static inline int min_exponent() { return std::numeric_limits<double>::min_exponent; }
    static inline int max_exponent() { return std::numeric_limits<double>::max_exponent; }
    static inline csym infinity() { return csym(std::numeric_limits<double>::infinity()); }
    static inline csym quiet_NaN() { return csym(std::numeric_limits<double>::quiet_NaN()); }
};
}
using SSpMat = Eigen::SparseMatrix<csym>;
using SVector = Eigen::Matrix<csym, Eigen::Dynamic, 1>;
using CSymLDLT = Eigen::SimplicialLDLT<SSpMat>;

// ---- Sparse pattern cache ----------------------------------------------------
// A frequency sweep re-solves systems with the SAME sparsity pattern and
// different VALUES at every point (the MQS block system: only β = ωμσ changes;
// the eigensolver's C = A − k²·(…): only k² changes). All the pattern-dependent
// work — the transpose-partner map for the value-symmetry check, the AMD
// ordering and the symbolic elimination analysis of SimplicialLDLT — is reused
// across calls; only the numeric factorization and the triangular solves run
// per call. Entries are keyed on the EXACT pattern (N + pointers + indices
// compared in full, no hashing) and kept in small per-call-site LRU stores (the
// odd and even differential modes alternate every sweep point, so a store needs
// at least two slots to avoid thrashing).
struct SparsePatternCache {
    int N = 0;
    std::vector<int> rowPtr, colIdx;      // the pattern this entry serves
    std::vector<int> transposeIdx;        // index of entry (j,i) for each entry k=(i,j)
    bool patternSymmetric = false;
    Eigen::SimplicialLDLT<RSpMat>* ldlt = nullptr;   // analyzePattern() already done
    CSymLDLT* cldlt = nullptr;                        // complex-symmetric twin (solve_complex_symmetric)
    // Largest |σ| at which the unpivoted LDLT failed the accuracy probe for this
    // pattern (the FEM pencil's curl null-space clusters eigenvalues near 0, so
    // LDLT reliably fails at SMALL |σ| and works at large |σ|). Future calls with
    // |σ| in the failed region skip straight to LU instead of paying a doomed
    // factorize+probe first.
    double ldltFailSigmaMax = 0.0;
    unsigned long lastUse = 0;
    ~SparsePatternCache() { delete ldlt; delete cldlt; }
};

static SparsePatternCache* findPatternCache(std::vector<SparsePatternCache*>& cache,
                                            int N, int nnz, const int* rowPtr, const int* colIdx) {
    static unsigned long useCounter = 0;
    const size_t MAX_ENTRIES = 4;
    useCounter++;
    for (SparsePatternCache* e : cache) {
        if (e->N == N && (int)e->colIdx.size() == nnz
            && std::equal(e->rowPtr.begin(), e->rowPtr.end(), rowPtr)
            && std::equal(e->colIdx.begin(), e->colIdx.end(), colIdx)) {
            e->lastUse = useCounter;
            return e;
        }
    }
    // Miss: build a new entry (transpose map + pattern-symmetry classification).
    // The arrays may be CSR or CSC — the logic is storage-order agnostic.
    SparsePatternCache* e = new SparsePatternCache();
    e->N = N;
    e->rowPtr.assign(rowPtr, rowPtr + N + 1);
    e->colIdx.assign(colIdx, colIdx + nnz);
    e->lastUse = useCounter;
    // CSC of A == CSR of Aᵀ, built by a counting pass. The pattern is symmetric
    // iff CSC pointers/indices equal the CSR ones (requires ascending indices per
    // row, which every producer in this codebase emits; a violation only demotes
    // the matrix to the un-cached LU path — never wrong results).
    std::vector<int> cscPtr(N + 1, 0), cscRow(nnz), cscK(nnz);
    for (int k = 0; k < nnz; k++) cscPtr[colIdx[k] + 1]++;
    for (int j = 0; j < N; j++) cscPtr[j + 1] += cscPtr[j];
    {
        std::vector<int> fill(cscPtr.begin(), cscPtr.end() - 1);
        for (int i = 0; i < N; i++)
            for (int k = rowPtr[i]; k < rowPtr[i + 1]; k++) {
                int p = fill[colIdx[k]]++;
                cscRow[p] = i; cscK[p] = k;
            }
    }
    e->patternSymmetric =
        std::equal(cscPtr.begin(), cscPtr.end(), e->rowPtr.begin())
        && std::equal(cscRow.begin(), cscRow.end(), e->colIdx.begin());
    if (e->patternSymmetric) e->transposeIdx = std::move(cscK);
    // LRU-evict beyond the cap (the factor objects hold real memory).
    if (cache.size() >= MAX_ENTRIES) {
        size_t victim = 0;
        for (size_t i = 1; i < cache.size(); i++)
            if (cache[i]->lastUse < cache[victim]->lastUse) victim = i;
        delete cache[victim];
        cache.erase(cache.begin() + victim);
    }
    cache.push_back(e);
    return e;
}

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
//
// If krylov_max > krylov_size, the SAME Krylov sequence is extended (subspace
// doubled, re-extracting Ritz pairs each time) until `nev` pairs pass the residual
// gate or the cap is reached. A single fixed-size pass converges only the Ritz
// pairs closest to the shift, and the strict gate below then silently DROPPED
// genuine far-from-shift modes (enclosure cavity modes with the shift at the
// quasi-TEM ε: only the quasi-TEM survived). Extension reuses the sparse LU (the
// dominant cost), so growing 20→80 vectors costs far less than a second solve.
// Callers that only need near-shift modes pass krylov_max == krylov_size and get
// the old single-pass behavior and cost.
static int solve_complex(
    int N,
    int nnz_a,
    int* a_rowPtr, int* a_colIdx, double* a_vals_re, double* a_vals_im,
    int nnz_b,
    int* b_rowPtr, int* b_colIdx, double* b_vals_re, double* b_vals_im,
    double sigma_re, double sigma_im,
    int num_eigenvalues, int krylov_size, int krylov_max,
    double* out_evals_re, double* out_evals_im,
    double* out_evecs_re, double* out_evecs_im,
    double* init_re, double* init_im
) {
    int n = N;
    int nev = num_eigenvalues;
    int m = krylov_size;

    if (m <= nev) m = 2 * nev + 1;
    if (m > n) m = n;

    // Subspace growth ceiling: caller's krylov_max, floored at m, capped by the
    // problem size and by keeping the basis V under ~128 MB for very large N.
    const int memCols = (int)(128.0e6 / (16.0 * n));
    const int mMax = std::min(n, std::max(m, std::min(krylov_max, memCols)));

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

    // Arnoldi iteration with modified Gram-Schmidt (double orthogonalization).
    // V grows column-blocks on demand; H is small, allocate the ceiling up front.
    CMatrix V(n, m + 1);
    CMatrix H = CMatrix::Zero(mMax + 1, mMax);

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

    // Extract the Ritz pairs of the leading actual_m×actual_m Hessenberg block,
    // walk them closest-to-shift first, and emit into the out arrays those whose
    // TRUE pencil residual passes the gate — up to nev. A rejected (unconverged)
    // pair's slot can be filled by a converged pair further from the shift. Only
    // the nearest 3·nev+8 candidates are residual-tested (each test is a V·y and
    // two spmv): anything further from the shift than that is junk.
    // Returns the converged count, or -2 on a failed dense decomposition.
    auto extract = [&](int actual_m) -> int {
        CMatrix Hm = H.topLeftCorner(actual_m, actual_m);
        if (!Hm.allFinite()) return -2;
        Eigen::ComplexEigenSolver<CMatrix> eig_solver;
        eig_solver.compute(Hm);
        if (eig_solver.info() != Eigen::Success) return -2;

        CVector theta = eig_solver.eigenvalues();
        CMatrix Y = eig_solver.eigenvectors();

        std::vector<int> indices(actual_m);
        for (int i = 0; i < actual_m; i++) indices[i] = i;
        std::sort(indices.begin(), indices.end(), [&](int a, int b) {
            return std::abs(theta(a)) > std::abs(theta(b));
        });

        const int kmax = std::min(actual_m, 3 * nev + 8);
        int nconv = 0;
        for (int k = 0; k < kmax && nconv < nev; k++) {
            int idx = indices[k];
            Complex th = theta(idx);

            if (std::abs(th) < 1e-15) continue;

            Complex lambda = sigma + 1.0 / th;

            CVector x = V.leftCols(actual_m) * Y.col(idx);

            CVector Ax = A * x;
            CVector Bx = B * x;
            CVector residual = Ax - lambda * Bx;
            double rel_res = residual.norm() / (x.norm() * (Ax.norm() + std::abs(lambda) * Bx.norm()));

            // Accept only Ritz pairs that actually converged within the Krylov
            // space. (This gate used to be `rel_res < 1e-4 || k < nev`, which
            // returned unconverged junk pairs as converged modes.)
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
    };

    int mCur = m;
    int j = 0;
    while (true) {
        int actual_m = mCur;
        bool breakdown = false;
        for (; j < mCur; j++) {
            CVector Bv = B * V.col(j);
            CVector w = solver.solve(Bv);

            // Check for NaN/Inf from solver
            if (!w.allFinite()) {
                actual_m = j;
                breakdown = true;
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

            if (!std::isfinite(wnorm)) {
                actual_m = j + 1;
                breakdown = true;
                break;
            }

            // Happy breakdown: V spans an invariant subspace (e.g. the start
            // vector was an exact eigenvector — the mode viewer's static seed is
            // nearly one). The wanted pairs OUTSIDE that subspace are not in it,
            // so deflate and continue: zero the subdiagonal (Hm becomes block
            // triangular, preserving the invariant block's exact Ritz pairs) and
            // restart the basis with a fresh vector orthogonal to everything so
            // the iteration can reach the rest of the spectrum.
            if (wnorm < 1e-14) {
                H(j + 1, j) = Complex(0.0, 0.0);
                CVector r(n);
                uint64_t st = 0x9E3779B97F4A7C15ull + (uint64_t)j;
                for (int i = 0; i < n; i++) {
                    st ^= st << 13; st ^= st >> 7; st ^= st << 17;
                    r(i) = Complex((double)(st & 0xFFFFFFull) / (double)0xFFFFFFull - 0.5, 0.0);
                }
                for (int pass = 0; pass < 2; pass++)
                    for (int i = 0; i <= j; i++) r -= V.col(i).dot(r) * V.col(i);
                double rn = r.norm();
                if (rn < 1e-12) {   // subspace genuinely exhausted
                    actual_m = j + 1;
                    breakdown = true;
                    break;
                }
                V.col(j + 1) = r / rn;
                continue;
            }

            V.col(j + 1) = w / wnorm;
        }

        if (actual_m < 1) return 0;

        int nconv = extract(actual_m);
        if (nconv < 0) return nconv;
        if (nconv >= nev || breakdown || mCur >= mMax) return nconv;

        mCur = std::min(2 * mCur, mMax);
        V.conservativeResize(Eigen::NoChange, mCur + 1);
    }
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
    int num_eigenvalues, int krylov_size, int krylov_max,
    double* out_evals_re, double* out_evals_im,
    double* out_evecs_re, double* out_evecs_im,
    double* init_re
) {
    int n = N;
    int nev = num_eigenvalues;
    int m = krylov_size;
    if (m <= nev) m = 2 * nev + 1;
    if (m > n) m = n;

    // Subspace growth ceiling (see solve_complex): real basis is 8 bytes/entry.
    const int memCols = (int)(128.0e6 / (8.0 * n));
    const int mMax = std::min(n, std::max(m, std::min(krylov_max, memCols)));

    RSpMat A = buildRealSparseMatrix(n, nnz_a, a_rowPtr, a_colIdx, a_vals_re);
    RSpMat B = buildRealSparseMatrix(n, nnz_b, b_rowPtr, b_colIdx, b_vals_re);

    // Factor C = A - σB. The FEM pencil is SYMMETRIC (indefinite), so the primary
    // factorization is SimplicialLDLT with the symbolic analysis cached across
    // calls (a frequency sweep re-factors the same pattern at every point — only
    // k² changes the values): considerably cheaper than general LU in both
    // factorize and the m Arnoldi back-solves. LDLT is unpivoted and can lose
    // accuracy on an indefinite matrix, so the factorization is verified with a
    // random-RHS residual probe; any doubt falls back to the pivoted SparseLU
    // (the previous always-on path). The final residual gate below checks every
    // returned Ritz pair against the TRUE pencil either way.
    RSpMat C = A - sigma_re * B;
    C.makeCompressed();
    bool ldltOK = false;
    // Per-solve iterative refinement is only worth its residual matvec + second
    // solve on every Arnoldi step when the raw factor is marginal. The probe
    // below measures the raw solve first: at ~1e-8 (typical for this pencil at
    // the quasi-TEM shift, measured 1e-9 to 2e-8) the operator is already four
    // orders inside the 1e-4 Ritz gate that verifies every returned pair against
    // the true pencil, so refinement is switched off for the whole call.
    bool ldltRefine = true;
    SparsePatternCache* pc = nullptr;
    {
        static std::vector<SparsePatternCache*> store;
        pc = findPatternCache(store, n, (int)C.nonZeros(), C.outerIndexPtr(), C.innerIndexPtr());
        if (pc->patternSymmetric && std::abs(sigma_re) > 1.5 * pc->ldltFailSigmaMax) {
            const double* cv = C.valuePtr();
            const int cnnz = (int)C.nonZeros();
            double scale = 0;
            for (int k = 0; k < cnnz; k++) scale = std::max(scale, std::abs(cv[k]));
            const double tol = 1e-12 * scale;
            bool valueSym = true;
            for (int k = 0; k < cnnz; k++) {
                if (std::abs(cv[k] - cv[pc->transposeIdx[k]]) > tol) { valueSym = false; break; }
            }
            if (valueSym) {
                if (!pc->ldlt) {
                    pc->ldlt = new Eigen::SimplicialLDLT<RSpMat>();
                    pc->ldlt->analyzePattern(C);
                }
                pc->ldlt->factorize(C);
                if (pc->ldlt->info() == Eigen::Success) {
                    // Probe with the same solve+iterative-refinement procedure the
                    // Arnoldi will use: a mildly inaccurate unpivoted factor (moderate
                    // element growth at an unlucky shift) is fully rescued by one or
                    // two refinement steps; only genuine breakdown falls back to LU.
                    RVector bp(n);
                    uint64_t st = 0xDA3E39CB94B95BDBull;
                    for (int i = 0; i < n; i++) {
                        st ^= st << 13; st ^= st >> 7; st ^= st << 17;
                        bp(i) = (double)(st & 0xFFFFFFull) / (double)0xFFFFFFull - 0.5;
                    }
                    RVector xp = pc->ldlt->solve(bp);
                    if (xp.allFinite()) {
                        double bn = bp.norm();
                        RVector r = bp - C * xp;
                        double rn = r.norm() / bn;
                        if (std::isfinite(rn) && rn < 1e-7) {
                            ldltOK = true;
                            ldltRefine = false;   // raw factor is good: no per-solve refinement
                        } else {
                            for (int it = 0; it < 3; it++) {
                                if (!std::isfinite(rn) || rn < 1e-11) break;
                                xp += pc->ldlt->solve(r);
                                r = bp - C * xp;
                                rn = r.norm() / bn;
                            }
                            // Same 1e-8 operator accuracy the plain-LU probe
                            // required, refinement just gets marginal factors
                            // there.
                            if (std::isfinite(rn) && rn < 1e-8) ldltOK = true;
                        }
                    }
                }
                if (!ldltOK)
                    pc->ldltFailSigmaMax = std::max(pc->ldltFailSigmaMax, std::abs(sigma_re));
            }
        }
    }
    Eigen::SparseLU<RSpMat> solver;
    if (!ldltOK) {
        solver.compute(C);
        if (solver.info() != Eigen::Success) {
            return -1;
        }
    }
    // Operator solve: LDLT + iterative refinement to LU-grade accuracy, or LU.
    auto solveC = [&](const RVector& v) -> RVector {
        if (!ldltOK) return RVector(solver.solve(v));
        RVector x = pc->ldlt->solve(v);
        if (!ldltRefine) return x;
        const double vn = v.norm();
        double prev = std::numeric_limits<double>::infinity();
        for (int it = 0; it < 2; it++) {
            RVector r = v - C * x;
            const double rn = r.norm();
            if (!(rn > 1e-10 * vn) || !(rn < 0.5 * prev)) break;  // done / stagnated / NaN
            prev = rn;
            x += pc->ldlt->solve(r);
        }
        return x;
    };

    RMatrix V(n, m + 1);
    RMatrix H = RMatrix::Zero(mMax + 1, mMax);

    RVector v0(n);
    if (init_re != nullptr) {
        for (int i = 0; i < n; i++) v0(i) = init_re[i];
        if (v0.norm() < 1e-30) v0.setConstant(1.0 / std::sqrt((double)n));
    } else {
        v0.setConstant(1.0 / std::sqrt((double)n));
    }
    V.col(0) = v0 / v0.norm();

    // Residual-gated Ritz extraction — the real twin of solve_complex's extract.
    auto extract = [&](int actual_m) -> int {
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

        const int kmax = std::min(actual_m, 3 * nev + 8);
        int nconv = 0;
        for (int k = 0; k < kmax && nconv < nev; k++) {
            int idx = indices[k];
            Complex th = theta(idx);
            if (std::abs(th) < 1e-15) continue;
            Complex lambda = sigma_re + 1.0 / th;

            // Ritz vector: complex combination of the real basis. A real Ritz
            // value (the common case on this real symmetric pencil, the
            // conjugate-pair case is the exception) comes with an exactly real
            // eigenvector from the real Schur form, so its imaginary half is
            // identically zero and its two sparse products are skipped. Half the
            // extraction cost for bit-identical result.
            const bool realPair = th.imag() == 0.0 && Y.col(idx).imag().isZero(0.0);
            RVector xr = V.leftCols(actual_m) * Y.col(idx).real();
            RVector xi = realPair ? RVector(RVector::Zero(n)) : RVector(V.leftCols(actual_m) * Y.col(idx).imag());

            // Residual ‖Ax − λBx‖ with real A, B and complex λ, x (componentwise)
            RVector Axr = A * xr, Bxr = B * xr;
            RVector Axi = realPair ? RVector(RVector::Zero(n)) : RVector(A * xi);
            RVector Bxi = realPair ? RVector(RVector::Zero(n)) : RVector(B * xi);
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
    };

    int mCur = m;
    int j = 0;
    while (true) {
        int actual_m = mCur;
        bool breakdown = false;
        for (; j < mCur; j++) {
            RVector Bv = B * V.col(j);
            RVector w = solveC(Bv);
            if (!w.allFinite()) { actual_m = j; breakdown = true; break; }
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
            if (!std::isfinite(wnorm)) { actual_m = j + 1; breakdown = true; break; }
            // Happy breakdown → deflate and continue with a fresh orthogonal
            // vector (see solve_complex for the rationale).
            if (wnorm < 1e-14) {
                H(j + 1, j) = 0.0;
                RVector r(n);
                uint64_t st = 0x9E3779B97F4A7C15ull + (uint64_t)j;
                for (int i = 0; i < n; i++) {
                    st ^= st << 13; st ^= st >> 7; st ^= st << 17;
                    r(i) = (double)(st & 0xFFFFFFull) / (double)0xFFFFFFull - 0.5;
                }
                for (int pass = 0; pass < 2; pass++)
                    for (int i = 0; i <= j; i++) r -= V.col(i).dot(r) * V.col(i);
                double rn = r.norm();
                if (rn < 1e-12) { actual_m = j + 1; breakdown = true; break; }
                V.col(j + 1) = r / rn;
                continue;
            }
            V.col(j + 1) = w / wnorm;
        }

        if (actual_m < 1) return 0;

        int nconv = extract(actual_m);
        if (nconv < 0) return nconv;
        if (nconv >= nev || breakdown || mCur >= mMax) return nconv;

        mCur = std::min(2 * mCur, mMax);
        V.conservativeResize(Eigen::NoChange, mCur + 1);
    }
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
    int num_eigenvalues, int krylov_size, int krylov_max,
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
                               sigma_re, num_eigenvalues, krylov_size, krylov_max,
                               out_evals_re, out_evals_im, out_evecs_re, out_evecs_im,
                               init_re);
        if (nconv > 0) return nconv;
    }
    return solve_complex(N, nnz_a, a_rowPtr, a_colIdx, a_vals_re, a_vals_im,
                         nnz_b, b_rowPtr, b_colIdx, b_vals_re, b_vals_im,
                         sigma_re, sigma_im, num_eigenvalues, krylov_size, krylov_max,
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

// krylov_max ≤ krylov_size → single fixed-size pass (old behavior). Larger →
// the Krylov subspace doubles (reusing the factorization) until num_eigenvalues
// pairs pass the residual gate or krylov_max columns are reached.
EMSCRIPTEN_KEEPALIVE
int solve_generalized_eigen(
    int N,
    int nnz_a,
    int* a_rowPtr, int* a_colIdx, double* a_vals_re, double* a_vals_im,
    int nnz_b,
    int* b_rowPtr, int* b_colIdx, double* b_vals_re, double* b_vals_im,
    double sigma_re, double sigma_im,
    int num_eigenvalues, int krylov_size, int krylov_max,
    double* out_evals_re, double* out_evals_im,
    double* out_evecs_re, double* out_evecs_im
) {
    // Catch everything (Spectra std::invalid_argument, Eigen std::bad_alloc, …):
    // an uncaught exception would abort() and permanently poison the module.
    try {
        return solve_core(N, nnz_a, a_rowPtr, a_colIdx, a_vals_re, a_vals_im,
                          nnz_b, b_rowPtr, b_colIdx, b_vals_re, b_vals_im,
                          sigma_re, sigma_im, num_eigenvalues, krylov_size, krylov_max,
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
    int num_eigenvalues, int krylov_size, int krylov_max,
    double* out_evals_re, double* out_evals_im,
    double* out_evecs_re, double* out_evecs_im,
    double* init_re, double* init_im
) {
    try {
        return solve_core(N, nnz_a, a_rowPtr, a_colIdx, a_vals_re, a_vals_im,
                          nnz_b, b_rowPtr, b_colIdx, b_vals_re, b_vals_im,
                          sigma_re, sigma_im, num_eigenvalues, krylov_size, krylov_max,
                          out_evals_re, out_evals_im, out_evecs_re, out_evecs_im,
                          init_re, init_im);
    } catch (...) {
        return -9;
    }
}

// Direct sparse solver: factorize once, solve for nRhs right-hand sides.
// Uses SimplicialLDLT for symmetric matrices, falls back to SparseLU.
// rhs_arr and x_arr are nRhs×N column-major (each RHS is contiguous N doubles).
// Returns 0 on success, negative on error.
EMSCRIPTEN_KEEPALIVE
int solve_sparse_multi(
    int N, int nnz,
    int* rowPtr, int* colIdx, double* values,
    int nRhs, double* rhs_arr, double* x_arr
) {
  try {
    static std::vector<SparsePatternCache*> store;
    SparsePatternCache* pc = findPatternCache(store, N, nnz, rowPtr, colIdx);

    // SimplicialLDLT reads only the lower triangle: for a NONSYMMETRIC input it
    // silently solves the wrong (symmetrized) system while reporting Success.
    // Only take the fast path when the matrix is actually symmetric — the
    // structural transpose map is cached with the pattern, so the per-call check
    // is a single O(nnz) value comparison (no sparse transpose/subtraction).
    bool symmetric = pc->patternSymmetric;
    if (symmetric) {
        double scale = 0;
        for (int k = 0; k < nnz; k++) scale = std::max(scale, std::abs(values[k]));
        const double tol = 1e-12 * scale;
        for (int k = 0; k < nnz; k++) {
            if (std::abs(values[k] - values[pc->transposeIdx[k]]) > tol) { symmetric = false; break; }
        }
    }

    RSpMat A = buildRealSparseMatrix(N, nnz, rowPtr, colIdx, values);

    // LDLT first (fast for SPD/symmetric indefinite), reusing the cached
    // symbolic analysis: analyzePattern() runs once per pattern, factorize()
    // (numeric only) per call.
    if (symmetric) {
        if (!pc->ldlt) {
            pc->ldlt = new Eigen::SimplicialLDLT<RSpMat>();
            pc->ldlt->analyzePattern(A);
        }
        pc->ldlt->factorize(A);
        if (pc->ldlt->info() == Eigen::Success) {
            for (int r = 0; r < nRhs; r++) {
                Eigen::Map<RVector> b(rhs_arr + r * N, N);
                Eigen::Map<RVector> x(x_arr + r * N, N);
                x = pc->ldlt->solve(b);
            }
            return 0;
        }
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


// Complex-symmetric direct solver: factor K = valRe + j·valIm (K^T = K) once and
// solve nRhs right-hand sides. Every RHS and every solution is 2N doubles: N
// real parts followed by N imaginary parts, the [Ar; Ai] layout the MQS
// post-processing indexes (mqs_loss.js). Unpivoted LDL^T on the csym scalar (see
// above), with the AMD ordering + symbolic analysis cached per pattern: a sweep
// re-solves the same skin mesh at every frequency with only β changed. The
// solution's residual is checked and a failure falls back to the pivoted complex
// SparseLU, so the caller always gets a correct answer or an error code.
// Returns 0 (LDL^T), 1 (solved by the SparseLU fallback), negative on failure.
EMSCRIPTEN_KEEPALIVE
int solve_complex_symmetric(
    int N, int nnz,
    int* rowPtr, int* colIdx, double* valRe, double* valIm,
    int nRhs, double* rhs_arr, double* x_arr
) {
  try {
    for (int i = 0; i < N; i++)
        for (int p = rowPtr[i]; p < rowPtr[i + 1]; p++) {
            if (colIdx[p] < 0 || colIdx[p] >= N) return -10;               // invalid column
            if (p > rowPtr[i] && colIdx[p] < colIdx[p - 1]) return -11;   // unsorted row
        }
    static std::vector<SparsePatternCache*> store;
    SparsePatternCache* pc = findPatternCache(store, N, nnz, rowPtr, colIdx);

    // LDL^T reads only the lower triangle: require value symmetry of BOTH parts
    // (cheap with the cached transpose map), else go straight to LU.
    bool symmetric = pc->patternSymmetric;
    if (symmetric) {
        double scale = 0;
        for (int k = 0; k < nnz; k++) scale = std::max(scale, std::abs(valRe[k]) + std::abs(valIm[k]));
        const double tol = 1e-12 * scale;
        for (int k = 0; k < nnz; k++) {
            const int t = pc->transposeIdx[k];
            if (std::abs(valRe[k] - valRe[t]) > tol || std::abs(valIm[k] - valIm[t]) > tol) { symmetric = false; break; }
        }
    }
    const size_t stride = (size_t)2 * N;
    if (symmetric) {
        SSpMat A(N, N);
        {
            std::vector<Eigen::Triplet<csym>> triplets;
            triplets.reserve(nnz);
            for (int i = 0; i < N; i++)
                for (int p = rowPtr[i]; p < rowPtr[i + 1]; p++)
                    triplets.emplace_back(i, colIdx[p], csym(valRe[p], valIm[p]));
            A.setFromTriplets(triplets.begin(), triplets.end());
        }
        A.makeCompressed();
        if (!pc->cldlt) {
            pc->cldlt = new CSymLDLT();
            pc->cldlt->analyzePattern(A);
        }
        CSymLDLT* ldlt = pc->cldlt;
        ldlt->factorize(A);
        bool ok = ldlt->info() == Eigen::Success;
        if (ok) {
            SVector b(N), x(N), res(N);
            for (int r = 0; r < nRhs && ok; r++) {
                const double* br = rhs_arr + r * stride;
                for (int i = 0; i < N; i++) b[i] = csym(br[i], br[N + i]);
                x = ldlt->solve(b);
                // Residual guard: the unpivoted factorization has no other way to
                // report a bad pivot (NaN or a lost digit), and it costs one
                // sparse product. 1e-8 relative is far above what a healthy
                // factorization produces (~1e-15) and far below any physical use.
                res = A * x - b;
                double rn = 0, bn = 0;
                for (int i = 0; i < N; i++) { rn += std::norm(res[i].v); bn += std::norm(b[i].v); }
                if (bn > 0 && !(rn <= 1e-16 * bn)) { ok = false; break; }
                double* xr = x_arr + r * stride;
                for (int i = 0; i < N; i++) { xr[i] = x[i].v.real(); xr[N + i] = x[i].v.imag(); }
            }
            if (ok) return 0;
        }
        // Drop the factor (it holds memory the LU below needs) and fall through.
        delete pc->cldlt;
        pc->cldlt = nullptr;
    }

    // Fallback: pivoted complex LU, robust, several times slower.
    SpMat A = buildSparseMatrix(N, nnz, rowPtr, colIdx, valRe, valIm);
    Eigen::SparseLU<SpMat> lu;
    lu.compute(A);
    if (lu.info() != Eigen::Success) return -1;
    CVector b(N), x(N);
    for (int r = 0; r < nRhs; r++) {
        const double* br = rhs_arr + r * stride;
        for (int i = 0; i < N; i++) b[i] = Complex(br[i], br[N + i]);
        x = lu.solve(b);
        if (lu.info() != Eigen::Success) return -2;
        double* xr = x_arr + r * stride;
        for (int i = 0; i < N; i++) { xr[i] = x[i].real(); xr[N + i] = x[i].imag(); }
    }
    return 1;
  } catch (...) {
    return -9;
  }
}

} // extern "C"
