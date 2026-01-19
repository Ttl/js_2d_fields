#include <emscripten.h>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>

using namespace Eigen;

extern "C" {
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
    try {
        // Build Eigen sparse matrix from CSR format
        SparseMatrix<double> A(N, N);
        A.reserve(nnz);
        
        std::vector<Triplet<double>> triplets;
        triplets.reserve(nnz);
        
        for (int i = 0; i < N; i++) {
            for (int p = rowPtr[i]; p < rowPtr[i + 1]; p++) {
                triplets.push_back(Triplet<double>(i, colIdx[p], values[p]));
            }
        }
        
        A.setFromTriplets(triplets.begin(), triplets.end());
        A.makeCompressed();
        
        // Map input/output arrays
        Map<VectorXd> b_vec(b, N);
        Map<VectorXd> x_vec(x_out, N);

        for (int i = 0; i < N; i++) {
            for (int p = rowPtr[i]; p < rowPtr[i + 1] - 1; p++) {
                if (colIdx[p] >= N || colIdx[p] < 0) {
                    return 10; // Invalid column index
                }
                if (colIdx[p] >= colIdx[p + 1]) {
                    return 11; // Unsorted columns
                }
            }
            if (rowPtr[i + 1] - 1 >= 0 && colIdx[rowPtr[i + 1] - 1] >= N) {
                return 10; // Invalid column index
            }
        }
                
        if (use_lu) {
            SparseLU<SparseMatrix<double>> solver;
            solver.compute(A);
            if (solver.info() != Success) {
                return 1; // Decomposition failed
            }
            x_vec = solver.solve(b_vec);
            if (solver.info() != Success) {
                return 2; // Solving failed
            }
        } else {
            SimplicialLDLT<SparseMatrix<double>> solver;
            solver.compute(A);
            if (solver.info() != Success) {
                return 3; // Decomposition failed
            }
            x_vec = solver.solve(b_vec);
            if (solver.info() != Success) {
                return 4; // Solving failed
            }
        }
        
        return 0;
    } catch (...) {
        return 99; // Unknown error
    }
}
}
