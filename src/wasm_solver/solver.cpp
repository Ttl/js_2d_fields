#include <emscripten.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>

using namespace Eigen;

typedef SparseMatrix<double> RSpMat;

namespace {

// Sparse pattern cache
//
// The quasi-static solve hits the same sparsity pattern several times in a row
// with different values. At one grid level the vacuum (eps = 1) and the
// dielectric operator share grid, conductor mask and stencil, and a frequency
// sweep re-solves the converged grid at every point. The pattern-dependent
// work, the transpose-partner map used by the value-symmetry check, plus
// SimplicialLDLT's AMD ordering and symbolic elimination analysis, is reused
// across those calls. Only the numeric factorization and the triangular solves
// run per call.
//
// Entries are keyed on the exact pattern (N and the index arrays compared in
// full) and evicted LRU. Two slots is the working set of one adaptive pass. A
// half-domain differential solve alternates a PEC and a PMC pattern, each used
// twice back to back. A cached entry holds a numeric factor, so the cap stays
// small.
struct SparsePatternCache {
    int N = 0;
    std::vector<int> rowPtr, colIdx;   // the pattern this entry serves
    std::vector<int> transposeIdx;     // index of entry (j,i) for each entry k=(i,j)
    bool patternSymmetric = false;
    SimplicialLDLT<RSpMat>* ldlt = nullptr;   // analyzePattern() already done
    unsigned long lastUse = 0;
    ~SparsePatternCache() { delete ldlt; }
};

SparsePatternCache* findPatternCache(std::vector<SparsePatternCache*>& cache,
                                     int N, int nnz, const int* rowPtr, const int* colIdx) {
    static unsigned long useCounter = 0;
    const size_t MAX_ENTRIES = 2;
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
    SparsePatternCache* e = new SparsePatternCache();
    e->N = N;
    e->rowPtr.assign(rowPtr, rowPtr + N + 1);
    e->colIdx.assign(colIdx, colIdx + nnz);
    e->lastUse = useCounter;
    // CSC of A == CSR of A^T, built by a counting pass. The pattern is symmetric
    // iff the CSC pointers/indices equal the CSR ones (requires ascending
    // indices per row, which the callers emit. A violation only demotes the
    // matrix to the LU path, never wrong results).
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
    // LRU-evict beyond the cap. The eviction happens before the new factor is
    // built, so the peak is MAX_ENTRIES numeric factors.
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

}  // namespace

extern "C" {

// Factor A once and solve nRhs right-hand sides. b and x_out are nRhs
// contiguous vectors of length N (column r at offset r*N). The certificate's
// odd/even mode solves share one operator, so amortizing the factorization
// (which dominates: the back-substitutions are cheap by comparison) halves the
// differential-pair cost.
//
// force_lu = 0 picks the factorization automatically: Cholesky (SimplicialLDLT,
// ~1.4x faster than LU on the FDM Laplace operator and roughly half the factor
// memory) when the matrix is symmetric, LU otherwise, and LU as a fallback if
// the Cholesky factorization reports failure. force_lu = 1 forces LU.
//
// The symmetry test is a real test, not an assumption: SimplicialLDLT reads
// only the lower triangle, so it would silently solve the wrong (symmetrized)
// system for an unsymmetric input while reporting Success.
EMSCRIPTEN_KEEPALIVE
int solve_sparse_multi(
    int N,
    int nnz,
    int* rowPtr,
    int* colIdx,
    double* values,
    int nRhs,
    double* b,
    double* x_out,
    int force_lu
) {
    try {
        // Validate the CSR arrays BEFORE building the matrix — setFromTriplets
        // performs the out-of-bounds writes a bad colIdx would cause, so the
        // check must precede it.
        for (int i = 0; i < N; i++) {
            for (int p = rowPtr[i]; p < rowPtr[i + 1]; p++) {
                if (colIdx[p] >= N || colIdx[p] < 0) {
                    return 10; // Invalid column index
                }
                if (p > rowPtr[i] && colIdx[p] < colIdx[p-1]) {
                    return 11; // Unsorted columns
                }
            }
        }

        static std::vector<SparsePatternCache*> store;
        SparsePatternCache* pc = nullptr;
        bool symmetric = false;
        if (!force_lu) {
            pc = findPatternCache(store, N, nnz, rowPtr, colIdx);
            symmetric = pc->patternSymmetric;
            if (symmetric) {
                double scale = 0;
                for (int k = 0; k < nnz; k++) scale = std::max(scale, std::abs(values[k]));
                const double tol = 1e-12 * scale;
                for (int k = 0; k < nnz; k++) {
                    if (std::abs(values[k] - values[pc->transposeIdx[k]]) > tol) { symmetric = false; break; }
                }
            }
        }

        // Build Eigen sparse matrix from CSR format
        RSpMat A(N, N);
        {
            std::vector<Triplet<double>> triplets;
            triplets.reserve(nnz);
            for (int i = 0; i < N; i++) {
                for (int p = rowPtr[i]; p < rowPtr[i + 1]; p++) {
                    triplets.push_back(Triplet<double>(i, colIdx[p], values[p]));
                }
            }
            A.setFromTriplets(triplets.begin(), triplets.end());
        }
        A.makeCompressed();

        if (symmetric) {
            if (!pc->ldlt) {
                pc->ldlt = new SimplicialLDLT<RSpMat>();
                pc->ldlt->analyzePattern(A);
            }
            pc->ldlt->factorize(A);
            if (pc->ldlt->info() == Success) {
                for (int r = 0; r < nRhs; r++) {
                    Map<VectorXd> b_vec(b + (size_t)r * N, N);
                    Map<VectorXd> x_vec(x_out + (size_t)r * N, N);
                    x_vec = pc->ldlt->solve(b_vec);
                    if (pc->ldlt->info() != Success) return 4; // Solving failed
                }
                return 0;
            }
            // Cholesky failed on this matrix, drop the factor (it holds memory
            // the LU below is about to need) and fall through.
            delete pc->ldlt;
            pc->ldlt = nullptr;
        }

        SparseLU<RSpMat> solver;
        solver.compute(A);
        if (solver.info() != Success) {
            return 1; // Decomposition failed
        }
        for (int r = 0; r < nRhs; r++) {
            Map<VectorXd> b_vec(b + (size_t)r * N, N);
            Map<VectorXd> x_vec(x_out + (size_t)r * N, N);
            x_vec = solver.solve(b_vec);
            if (solver.info() != Success) {
                return 2; // Solving failed
            }
        }

        return 0;
    } catch (...) {
        return 99; // Unknown error
    }
}

EMSCRIPTEN_KEEPALIVE
int solve_sparse(
    int N,
    int nnz,
    int* rowPtr,
    int* colIdx,
    double* values,
    double* b,
    double* x_out,
    int force_lu
) {
    return solve_sparse_multi(N, nnz, rowPtr, colIdx, values, 1, b, x_out, force_lu);
}

}
