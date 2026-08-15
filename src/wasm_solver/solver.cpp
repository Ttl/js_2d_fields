#include <emscripten.h>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>

using namespace Eigen;

extern "C" {

// Factor A once and solve nRhs right-hand sides. b and x_out are nRhs
// contiguous vectors of length N (column r at offset r*N). The certificate's
// odd/even mode solves share one operator, so amortizing the factorization
// (which dominates: the back-substitutions are cheap by comparison) halves the
// differential-pair cost.
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
    int use_lu
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

        // Build Eigen sparse matrix from CSR format
        SparseMatrix<double> A(N, N);

        std::vector<Triplet<double>> triplets;
        triplets.reserve(nnz);

        for (int i = 0; i < N; i++) {
            for (int p = rowPtr[i]; p < rowPtr[i + 1]; p++) {
                triplets.push_back(Triplet<double>(i, colIdx[p], values[p]));
            }
        }

        A.setFromTriplets(triplets.begin(), triplets.end());
        A.makeCompressed();

        if (use_lu) {
            SparseLU<SparseMatrix<double>> solver;
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
        } else {
            SimplicialLDLT<SparseMatrix<double>> solver;
            solver.compute(A);
            if (solver.info() != Success) {
                return 3; // Decomposition failed
            }
            for (int r = 0; r < nRhs; r++) {
                Map<VectorXd> b_vec(b + (size_t)r * N, N);
                Map<VectorXd> x_vec(x_out + (size_t)r * N, N);
                x_vec = solver.solve(b_vec);
                if (solver.info() != Success) {
                    return 4; // Solving failed
                }
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
    int use_lu
) {
    return solve_sparse_multi(N, nnz, rowPtr, colIdx, values, 1, b, x_out, use_lu);
}

}
