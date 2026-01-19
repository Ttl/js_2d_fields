const CONSTANTS = {
    EPS0: 8.854187817e-12,
    MU0: 4 * Math.PI * 1e-7,
    C: 299792458,
    PI: Math.PI
};

// --- Math Utils ---

function linspace(start, end, n) {
    if (n <= 1) return [start];
    const arr = new Float64Array(n);
    const step = (end - start) / (n - 1);
    for (let i = 0; i < n; i++) arr[i] = start + i * step;
    return arr;
}

function diff(arr) {
    const res = new Float64Array(arr.length - 1);
    for (let i = 0; i < arr.length - 1; i++) res[i] = arr[i+1] - arr[i];
    return res;
}

function smooth_transition(start, end, n_points, curve_end='end', beta=4.0) {
    if (n_points <= 1) return new Float64Array([start]);
    if (n_points === 2) return new Float64Array([start, end]);

    const xi = linspace(0, 1, n_points);
    const res = new Float64Array(n_points);

    for(let i=0; i<n_points; i++) {
        let eta = 0;
        if (curve_end === 'end') {
            eta = Math.tanh(beta * xi[i]) / Math.tanh(beta);
        } else if (curve_end === 'both') {
            eta = (Math.tanh(beta * (xi[i] - 0.5)) / Math.tanh(beta * 0.5) + 1) / 2;
        } else {
            eta = 1 - Math.tanh(beta * (1 - xi[i])) / Math.tanh(beta);
        }
        res[i] = start + eta * (end - start);
    }
    return res;
}

function buildCSR(colLists, valLists, N) {
    let nnz = 0;
    for (let i = 0; i < N; i++) nnz += colLists[i].length;

    const rowPtr = new Int32Array(N + 1);
    const colIdx = new Int32Array(nnz);
    const values = new Float64Array(nnz);

    let p = 0;
    for (let i = 0; i < N; i++) {
        rowPtr[i] = p;
        const cols = colLists[i];
        const vals = valLists[i];
        
        // Create array of (col, val) pairs and sort by column index
        const pairs = [];
        for (let k = 0; k < cols.length; k++) {
            pairs.push({ col: cols[k], val: vals[k] });
        }
        pairs.sort((a, b) => a.col - b.col);
        
        // Write sorted data
        for (let k = 0; k < pairs.length; k++) {
            colIdx[p] = pairs[k].col;
            values[p] = pairs[k].val;
            p++;
        }
    }
    rowPtr[N] = p;

    return { rowPtr, colIdx, values };
}

async function solveWithWASM(csr, B, useLU = false) {
    await Module.ready;
    
    const N = B.length;
    const nnz = csr.values.length;
    
    console.log(`Solving: N=${N}, nnz=${nnz}, useLU=${useLU}`);
    
    // Allocate memory
    const pRow = Module._malloc(4 * (N + 1));
    const pCol = Module._malloc(4 * nnz);
    const pVal = Module._malloc(8 * nnz);
    const pB   = Module._malloc(8 * N);
    const pX   = Module._malloc(8 * N);
    
    try {
        // Copy data - create views AFTER malloc
        const rowView = new Int32Array(Module.HEAP32.buffer, pRow, N + 1);
        const colView = new Int32Array(Module.HEAP32.buffer, pCol, nnz);
        const valView = new Float64Array(Module.HEAPF64.buffer, pVal, nnz);
        const bView = new Float64Array(Module.HEAPF64.buffer, pB, N);
        
        rowView.set(csr.rowPtr);
        colView.set(csr.colIdx);
        valView.set(csr.values);
        bView.set(B);

        if (!Module._solve_sparse) {
            throw new Error("WASM function solve_sparse not found. Module not loaded properly.");
        }
        
        // Call solver
        const status = Module._solve_sparse(
            N, nnz,
            pRow, pCol, pVal,
            pB, pX,
            useLU ? 1 : 0
        );
        
        console.log(`Solver status: ${status}`);
        
        if (status !== 0) {
            const errors = {
                1: "LU decomposition failed",
                2: "LU solving failed", 
                3: "Cholesky decomposition failed (matrix may not be positive definite)",
                4: "Cholesky solving failed",
                99: "Unknown C++ exception"
            };
            throw new Error(errors[status] || `WASM solver failed with code: ${status}`);
        }
        
        // Copy result
        const xView = new Float64Array(Module.HEAPF64.buffer, pX, N);
        const x = new Float64Array(N);
        x.set(xView);
        
        return x;
    } finally {
        // Always free memory
        //Module._free(pRow);
        //Module._free(pCol);
        //Module._free(pVal);
        //Module._free(pB);
        //Module._free(pX);
    }
}

// --- Solver Class ---

class FieldSolver2D {
    constructor() {
        this.x = null;
        this.y = null;
        this.V = null;
        this.epsilon_r = null;
        this.conductor_mask = null; // 1 if conductor, 0 if dielectric
        this.bc_mask = null; // Boundary conditions
        this.solution_valid = false;

        // Computed fields
        this.Ex = null;
        this.Ey = null;
    }

    async solve_laplace_iterative(vacuum = false, onProgress = null) {
        const ny = this.y.length, nx = this.x.length;
        const dx = diff(this.x), dy = diff(this.y);
        const N = nx * ny;
        const idx = (i, j) => i * nx + j;

        const get_er = (i, j) => vacuum ? 1.0 : this.epsilon_r[i][j];
        const is_cond = (i, j) => this.conductor_mask[i][j];

        // --- Reduced system mapping ---
        const is_unknown = new Int8Array(N);
        let N_unknown = 0;

        for (let i = 0; i < ny; i++)
            for (let j = 0; j < nx; j++) {
                const n = idx(i, j);
                if (!is_cond(i, j)) {
                    is_unknown[n] = 1;
                    N_unknown++;
                }
            }

        const full_to_red = new Int32Array(N).fill(-1);
        const red_to_full = new Int32Array(N_unknown);

        let k = 0;
        for (let n = 0; n < N; n++) {
            if (is_unknown[n]) {
                full_to_red[n] = k;
                red_to_full[k] = n;
                k++;
            }
        }

        // --- Build sparse system (lists first) ---
        const B = new Float64Array(N_unknown);
        const diag = new Float64Array(N_unknown);

        const colLists = Array(N_unknown);  // Column indices for each row
        const valLists = Array(N_unknown);  // Values for each row

        for (let i = 0; i < N_unknown; i++) {
            colLists[i] = [];
            valLists[i] = [];
        }

        const addA = (r, c, v) => {
            colLists[r].push(c);
            valLists[r].push(v);
            if (r === c) diag[r] += v;
        };

        // --- Matrix assembly (Python parity) ---
        for (let i = 0; i < ny; i++) {
            for (let j = 0; j < nx; j++) {
                if (is_cond(i, j)) continue;

                const fn = idx(i, j);
                const n = full_to_red[fn];

                const boundary =
                    i === 0 || i === ny - 1 || j === 0 || j === nx - 1;

                let dxr, dxl, dyu, dyd;
                if (boundary) {
                    dxr = j < nx - 1 ? dx[j] : dx[j - 1];
                    dxl = j > 0 ? dx[j - 1] : dx[j];
                    dyu = i < ny - 1 ? dy[i] : dy[i - 1];
                    dyd = i > 0 ? dy[i - 1] : dy[i];
                } else {
                    dxr = dx[j];
                    dxl = dx[j - 1];
                    dyu = dy[i];
                    dyd = dy[i - 1];
                }

                let err, erl, eru, erd;
                if (vacuum) {
                    err = erl = eru = erd = 1.0;
                } else {
                    const erc = get_er(i, j);
                    if (boundary) {
                        err = 0.5 * (erc + get_er(i, Math.min(j + 1, nx - 1)));
                        erl = 0.5 * (erc + get_er(i, Math.max(j - 1, 0)));
                        eru = 0.5 * (erc + get_er(Math.min(i + 1, ny - 1), j));
                        erd = 0.5 * (erc + get_er(Math.max(i - 1, 0), j));
                    } else {
                        err = is_cond(i, j + 1) ? erc : 0.5 * (erc + get_er(i, j + 1));
                        erl = is_cond(i, j - 1) ? erc : 0.5 * (erc + get_er(i, j - 1));
                        eru = is_cond(i + 1, j) ? erc : 0.5 * (erc + get_er(i + 1, j));
                        erd = is_cond(i - 1, j) ? erc : 0.5 * (erc + get_er(i - 1, j));
                    }
                }

                const area_i = 0.5 * (dyd + dyu);
                const area_j = 0.5 * (dxl + dxr);

                let cr, cl, cu, cd;
                if (boundary) {
                    cr = j < nx - 1 ? -err * area_i / dxr : 0;
                    cl = j > 0 ? -erl * area_i / dxl : 0;
                    cu = i < ny - 1 ? -eru * area_j / dyu : 0;
                    cd = i > 0 ? -erd * area_j / dyd : 0;
                } else {
                    cr = -err * area_i / dxr;
                    cl = -erl * area_i / dxl;
                    cu = -eru * area_j / dyu;
                    cd = -erd * area_j / dyd;
                }

                const cc = -(cr + cl + cu + cd);
                addA(n, n, cc);

                const handle = (ii, jj, c) => {
                    const fn2 = idx(ii, jj);
                    if (!is_cond(ii, jj)) {
                        addA(n, full_to_red[fn2], c);
                    } else {
                        B[n] -= c * this.V[ii][jj];
                    }
                };

                if (j < nx - 1) handle(i, j + 1, cr);
                if (j > 0) handle(i, j - 1, cl);
                if (i < ny - 1) handle(i + 1, j, cu);
                if (i > 0) handle(i - 1, j, cd);
            }
        }

        // --- Convert to CSR ---
        const { rowPtr, colIdx, values } = buildCSR(colLists, valLists, N_unknown);

        const csr = { rowPtr, colIdx, values };
        const x = await solveWithWASM(csr, B, true); // false = Cholesky

            /*
        // --- Preconditioned CG ---
        const x = new Float64Array(N_unknown);
        const r = new Float64Array(N_unknown);
        const z = new Float64Array(N_unknown);
        const p = new Float64Array(N_unknown);
        const Ap = new Float64Array(N_unknown);
        const invDiag = new Float64Array(N_unknown);

        for (let i = 0; i < N_unknown; i++) invDiag[i] = 1.0 / diag[i];

        let bnorm = 0;
        for (let i = 0; i < N_unknown; i++) {
            r[i] = B[i];
            bnorm += B[i] * B[i];
        }
        bnorm = Math.sqrt(bnorm);

        const rtol = 1e-7;
        const atol = 1e-12;
        const maxIter = 10000;

        for (let i = 0; i < N_unknown; i++) {
            z[i] = r[i] * invDiag[i];
            p[i] = z[i];
        }

        let rz_old = 0;
        for (let i = 0; i < N_unknown; i++) rz_old += r[i] * z[i];

        for (let iter = 0; iter < maxIter; iter++) {
            // Ap = A * p
            for (let i = 0; i < N_unknown; i++) {
                let sum = 0;
                for (let k = rowPtr[i]; k < rowPtr[i + 1]; k++) {
                    sum += values[k] * p[colIdx[k]];
                }
                Ap[i] = sum;
            }

            let alpha_den = 0;
            for (let i = 0; i < N_unknown; i++) {
                alpha_den += p[i] * Ap[i];
            }

            const alpha = rz_old / alpha_den;

            let err = 0;
            for (let i = 0; i < N_unknown; i++) {
                x[i] += alpha * p[i];
                r[i] -= alpha * Ap[i];
                err += r[i] * r[i];
            }

            const rnorm = Math.sqrt(err);
            if (rnorm <= rtol * bnorm + atol) break;

            for (let i = 0; i < N_unknown; i++) {
                z[i] = r[i] * invDiag[i];
            }

            let rz_new = 0;
            for (let i = 0; i < N_unknown; i++) rz_new += r[i] * z[i];

            const beta = rz_new / rz_old;
            rz_old = rz_new;

            for (let i = 0; i < N_unknown; i++) {
                p[i] = z[i] + beta * p[i];
            }

            if (iter % 100 === 0 && onProgress) {
                onProgress(iter, maxIter, rnorm / bnorm);
                await new Promise(r => setTimeout(r, 0));
            }
        }
       */

        // --- Reconstruct solution ---
        for (let k = 0; k < N_unknown; k++) {
            const n = red_to_full[k];
            const i = (n / nx) | 0;
            const j = n % nx;
            this.V[i][j] = x[k];
        }
    }

    compute_fields() {
        const ny = this.y.length;
        const nx = this.x.length;
        const dx = diff(this.x);
        const dy = diff(this.y);

        this.Ex = Array(ny).fill().map(() => new Float64Array(nx));
        this.Ey = Array(ny).fill().map(() => new Float64Array(nx));

        for(let i=1; i<ny-1; i++) {
            for(let j=1; j<nx-1; j++) {
                if (this.conductor_mask[i][j]) continue;

                const dxl = dx[j-1];
                const dxr = dx[j];
                const dyd = dy[i-1];
                const dyu = dy[i];

                this.Ex[i][j] = -(
                    (dxl / (dxr * (dxl + dxr))) * this.V[i][j+1] +
                    ((dxr - dxl) / (dxl * dxr)) * this.V[i][j] -
                    (dxr / (dxl * (dxl + dxr))) * this.V[i][j-1]
                );

                this.Ey[i][j] = -(
                    (dyd / (dyu * (dyd + dyu))) * this.V[i+1][j] +
                    ((dyu - dyd) / (dyd * dyu)) * this.V[i][j] -
                    (dyu / (dyd * (dyd + dyu))) * this.V[i-1][j]
                );
            }
        }
        this.solution_valid = true;
    }

    calculate_capacitance(vacuum=false) {
        let Q = 0.0;
        const ny = this.y.length;
        const nx = this.x.length;
        const dx = diff(this.x);
        const dy = diff(this.y);

        const get_dx = j => j < dx.length ? dx[j] : dx[dx.length-1];
        const get_dy = i => i < dy.length ? dy[i] : dy[dy.length-1];

        // Iterate over signal trace interface
        for (let i = 1; i < ny - 1; i++) {
            for (let j = 1; j < nx - 1; j++) {
                if (!this.signal_mask[i][j]) continue;

                // Check 4 neighbors
                const check_neighbor = (ni, nj, is_vertical_flux) => {
                    if (this.signal_mask[ni][nj]) return; // Internal to conductor

                    // E-field Normal
                    let En;
                    let dist;
                    let area;

                    if (is_vertical_flux) {
                         // Neighbor is Top/Bottom
                         dist = Math.abs(this.y[i] - this.y[ni]);
                         En = (this.V[i][j] - this.V[ni][nj]) / dist;
                         // Average dx for area
                         area = (get_dx(j-1) + get_dx(j)) / 2;
                    } else {
                        // Neighbor is Left/Right
                        dist = Math.abs(this.x[j] - this.x[nj]);
                        En = (this.V[i][j] - this.V[ni][nj]) / dist;
                        // Average dy for area
                        area = (get_dy(i-1) + get_dy(i)) / 2;
                    }

                    const er = vacuum ? 1 : this.epsilon_r[ni][nj];
                    Q += CONSTANTS.EPS0 * er * En * area;
                };

                check_neighbor(i, j+1, false); // Right
                check_neighbor(i, j-1, false); // Left
                check_neighbor(i+1, j, true);  // Top
                check_neighbor(i-1, j, true);  // Bottom
            }
        }
        return Math.abs(Q);
    }
}

