import createWASMModule from './wasm_solver/solver.js';

export const CONSTANTS = {
    EPS0: 8.854187817e-12,
    MU0: 4 * Math.PI * 1e-7,
    C: 299792458,
    PI: Math.PI
};

// --- Complex Number Class ---
class Complex {
    constructor(re, im = 0) {
        this.re = re;
        this.im = im;
    }

    add(other) {
        return new Complex(this.re + other.re, this.im + other.im);
    }

    sub(other) {
        return new Complex(this.re - other.re, this.im - other.im);
    }

    mul(other) {
        if (typeof other === 'number') {
            return new Complex(this.re * other, this.im * other);
        }
        return new Complex(
            this.re * other.re - this.im * other.im,
            this.re * other.im + this.im * other.re
        );
    }

    div(other) {
        if (typeof other === 'number') {
            return new Complex(this.re / other, this.im / other);
        }
        const denominator = other.re * other.re + other.im * other.im;
        return new Complex(
            (this.re * other.re + this.im * other.im) / denominator,
            (this.im * other.re - this.re * other.im) / denominator
        );
    }

    sqrt() {
        const r = Math.sqrt(this.re * this.re + this.im * this.im);
        const theta = Math.atan2(this.im, this.re);
        return new Complex(
            Math.sqrt(r) * Math.cos(theta / 2),
            Math.sqrt(r) * Math.sin(theta / 2)
        );
    }

    toString() {
        if (this.im === 0) return `${this.re.toFixed(2)}`;
        if (this.re === 0) return `${this.im.toFixed(2)}j`;
        return `${this.re.toFixed(2)}${this.im > 0 ? '+' : ''}${this.im.toFixed(2)}j`;
    }
}

// --- Math Utils ---

export function diff(arr) {
    const res = new Float64Array(arr.length - 1);
    for (let i = 0; i < arr.length - 1; i++) res[i] = arr[i+1] - arr[i];
    return res;
}

function linspace(start, end, n) {
    if (n <= 1) return [start];
    const arr = new Float64Array(n);
    const step = (end - start) / (n - 1);
    for (let i = 0; i < n; i++) arr[i] = start + i * step;
    return arr;
}

export function _smooth_transition(x0, x1, n, curve_end, beta=4.0) {
    if (n <= 1) {
        if (n === 0) return new Float64Array(0); // If 0 points, return empty array
        return new Float64Array([x0, x1]); // Python returns [start, end] for n_points <= 1
    }
    
    const pts = new Float64Array(n);
    const xi = linspace(0, 1, n);

    for(let i=0; i<n; i++) {
        let eta;
        if (curve_end === 'end') {
            eta = Math.tanh(beta * xi[i]) / Math.tanh(beta);
        } else if (curve_end === 'both') {
            eta = (Math.tanh(beta * (xi[i] - 0.5)) / Math.tanh(beta * 0.5) + 1) / 2;
        } else { // 'start'
            eta = 1 - Math.tanh(beta * (1 - xi[i])) / Math.tanh(beta);
        }
        pts[i] = x0 + eta * (x1 - x0);
    }
    return pts;
}

export function _enforce_interfaces(arr, interfaces) {
    let list = Array.from(arr);

    // Always add all interface values, then use Set to handle uniqueness and sort.
    // This directly mirrors the Python np.append and np.sort behavior (with unique implied).
    for (let val of interfaces) {
        list.push(val);
    }
    
    return Float64Array.from(new Set(list)).sort((a, b) => a - b);
}

export function _concat_arrays(arrays) {
    let totalLen = 0;
    for(let a of arrays) totalLen += a.length;
    let res = new Float64Array(totalLen);
    let offset = 0;
    for(let a of arrays) {
        res.set(a, offset);
        offset += a.length;
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

// Store the initialized WASM module (singleton pattern)
let WASMModuleInstance = null;

async function solveWithWASM(csr, B, useLU = false) {
    if (!WASMModuleInstance) {
        // Initialize the module if it hasn't been already
        WASMModuleInstance = await createWASMModule();
    }
    
    const N = B.length;
    const nnz = csr.values.length;
    
    // Allocate memory
    const pRow = WASMModuleInstance._malloc(4 * (N + 1));
    const pCol = WASMModuleInstance._malloc(4 * nnz);
    const pVal = WASMModuleInstance._malloc(8 * nnz);
    const pB   = WASMModuleInstance._malloc(8 * N);
    const pX   = WASMModuleInstance._malloc(8 * N);

    try {
        // Re-acquire HEAP views to ensure they are current in case memory grew
        const currentHEAP32 = WASMModuleInstance.HEAP32;
        const currentHEAPF64 = WASMModuleInstance.HEAPF64;

        // Copy data - create views AFTER malloc
        const rowView = new Int32Array(currentHEAP32.buffer, pRow, N + 1);
        const colView = new Int32Array(currentHEAP32.buffer, pCol, nnz);
        const valView = new Float64Array(currentHEAPF64.buffer, pVal, nnz);
        const bView = new Float64Array(currentHEAPF64.buffer, pB, N);
        
        rowView.set(csr.rowPtr);
        colView.set(csr.colIdx);
        valView.set(csr.values);
        bView.set(B);

        if (!WASMModuleInstance._solve_sparse) {
            throw new Error("WASM function solve_sparse not found. Module not loaded properly.");
        }
        
        // Call solver
        const status = WASMModuleInstance._solve_sparse(
            N, nnz,
            pRow, pCol, pVal,
            pB, pX,
            useLU ? 1 : 0
        );
        
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
        const xView = new Float64Array(WASMModuleInstance.HEAPF64.buffer, pX, N);
        const x = new Float64Array(N);
        x.set(xView);
        
        return x;
    } finally {
        // Always free memory
        WASMModuleInstance._free(pRow);
        WASMModuleInstance._free(pCol);
        WASMModuleInstance._free(pVal);
        WASMModuleInstance._free(pB);
        WASMModuleInstance._free(pX);
    }
}

// --- Solver Class ---

export class FieldSolver2D {
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
        const x = await solveWithWASM(csr, B, true);

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
                    // Only add flux if the neighbor is NOT part of the signal conductor
                    if (this.signal_mask[ni][nj]) return;

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

                // Changed these to use the bounds checking from the python side
                // No need to check for signal_mask[i,j+1] etc. here, as check_neighbor does it.
                // The python code calculates the flux across faces *between* a signal and a dielectric region.

                // Right neighbor
                if (!this.signal_mask[i][j + 1]) {
                    check_neighbor(i, j + 1, false);
                }
                // Left neighbor
                if (!this.signal_mask[i][j - 1]) {
                    check_neighbor(i, j - 1, false);
                }
                // Top neighbor
                if (!this.signal_mask[i + 1][j]) {
                    check_neighbor(i + 1, j, true);
                }
                // Bottom neighbor
                if (!this.signal_mask[i - 1][j]) {
                    check_neighbor(i - 1, j, true);
                }
            }
        }
        return Math.abs(Q);
    }

    calculate_conductor_loss(Ex, Ey, Z0) {
        if (!this.solution_valid) {
            throw new Error("Fields (Ex, Ey) are not valid. Run compute_fields() first.");
        }
        const ny = this.y.length;
        const nx = this.x.length;
        const dx_array = diff(this.x);
        const dy_array = diff(this.y);

        const Rs = Math.sqrt(this.omega * CONSTANTS.MU0 / (2 * this.sigma_cond));
        const delta = Math.sqrt(2 / (this.omega * CONSTANTS.MU0 * this.sigma_cond));

        const get_dx = j => (j >= 0 && j < dx_array.length) ? dx_array[j] : dx_array[dx_array.length - 1];
        const get_dy = i => (i >= 0 && i < dy_array.length) ? dy_array[i] : dy_array[dy_array.length - 1];

        let Pc = 0.0;

        for (let i = 1; i < ny - 1; i++) {
            for (let j = 1; j < nx - 1; j++) {
                if (!(this.signal_mask[i][j] || this.ground_mask[i][j])) {
                    continue;
                }

                const neighbors = [
                    { ni: i, nj: j + 1, direction: 'r', dl_func: get_dy, idx: i },
                    { ni: i, nj: j - 1, direction: 'l', dl_func: get_dy, idx: i },
                    { ni: i + 1, nj: j, direction: 'u', dl_func: get_dx, idx: j },
                    { ni: i - 1, nj: j, direction: 'd', dl_func: get_dx, idx: j },
                ];

                let cell_K_sq = 0.0;
                let cell_dl = 0.0;

                for (const { ni, nj, direction, dl_func, idx: dl_idx } of neighbors) {
                    if (ni < 0 || ni >= ny || nj < 0 || nj >= nx) continue; // Out of bounds
                    if (this.signal_mask[ni][nj] || this.ground_mask[ni][nj]) {
                        continue; // Neighbor is also conductor
                    }

                    const eps_diel = this.epsilon_r[ni][nj];
                    const Ex_diel = Ex[ni][nj];
                    const Ey_diel = Ey[ni][nj];

                    let E_norm;
                    if (direction === 'r') {
                        E_norm = Ex_diel;
                    } else if (direction === 'l') {
                        E_norm = -Ex_diel;
                    } else if (direction === 'u') {
                        E_norm = Ey_diel;
                    } else if (direction === 'd') {
                        E_norm = -Ey_diel;
                    } else {
                        E_norm = 0; // Should not happen
                    }

                    const Z0_freespace = Math.sqrt(CONSTANTS.MU0 / CONSTANTS.EPS0);
                    const H_tan = Math.abs(E_norm) * Math.sqrt(eps_diel) / Z0_freespace;
                    const K = H_tan;
                    const dl = dl_func(dl_idx);

                    cell_K_sq += K * K * dl;
                    cell_dl += dl;
                }

                if (cell_dl > 0) {
                    const dP = 0.5 * Rs * cell_K_sq;
                    Pc += dP;
                }
            }
        }

        const Pc_1W = Pc * 2 * Z0; // Normalized to 1W power flow
        const alpha_db_per_m = 8.686 * Pc_1W / 2.0;

        return alpha_db_per_m;
    }

    calculate_dielectric_loss(Ex, Ey, Z0) {
        if (!this.solution_valid) {
            throw new Error("Fields (Ex, Ey) are not valid. Run compute_fields() first.");
        }
        const ny = this.y.length;
        const nx = this.x.length;
        const dx_array = diff(this.x);
        const dy_array = diff(this.y);

        const get_dx = j => (j >= 0 && j < dx_array.length) ? dx_array[j] : dx_array[dx_array.length - 1];
        const get_dy = i => (i >= 0 && i < dy_array.length) ? dy_array[i] : dy_array[dy_array.length - 1];

        let Pd = 0.0;

        for (let i = 0; i < ny - 1; i++) {
            for (let j = 0; j < nx - 1; j++) {
                if (this.conductor_mask[i][j]) continue;
                if (this.epsilon_r[i][j] <= 1.01) continue; // Air or vacuum, no loss

                const E2 = Ex[i][j] * Ex[i][j] + Ey[i][j] * Ey[i][j];
                const dA = get_dx(j) * get_dy(i);
                
                // Assuming tan_delta is a property of the class (MicrostripSolver2D or FieldSolver2D)
                // and omega is also available (2 * pi * freq)
                if (this.tan_delta && this.omega) {
                    Pd += 0.5 * this.omega * CONSTANTS.EPS0 * this.epsilon_r[i][j] * this.tan_delta * E2 * dA;
                }
            }
        }

        const P_flow = 1.0 / (2 * Z0);
        return 8.686 * (Pd / (2 * P_flow));
    }

    rlgc(alpha_cond, alpha_diel, C_mode, Z0_mode) {
        // Convert dB/m to Np/m
        const alpha_c_np = alpha_cond / 8.686;
        const alpha_d_np = alpha_diel / 8.686;

        const R = 2 * Z0_mode * alpha_c_np;
        const G = 2 * alpha_d_np / Z0_mode;

        const L = (Z0_mode * Z0_mode) * C_mode;

        const eps_eff_mode = (L * C_mode) / (CONSTANTS.MU0 * CONSTANTS.EPS0);

        // Complex characteristic impedance Zc = sqrt((R + jwL) / (G + jwC))
        const num = new Complex(R, this.omega * L);
        const den = new Complex(G, this.omega * C_mode);
        const Zc = (num.div(den)).sqrt();

        const rlgc = {
            R: R,
            L: L,
            G: G,
            C: C_mode
        };

        return { Zc, rlgc, eps_eff_mode };
    }

    async perform_analysis(onProgress = null) {
        const totalSteps = 2; // Two main solve_laplace_iterative calls
        let currentStep = 0;

        const updateProgress = (stepFraction, overallStep) => {
            if (onProgress) {
                const totalProgress = ((overallStep - 1) / totalSteps) + (stepFraction / totalSteps);
                onProgress(totalProgress);
            }
        };

        // 1. Calculate C0 (vacuum capacitance)
        currentStep = 1;
        await this.solve_laplace_iterative(true, (i, max) => updateProgress(i / max, currentStep));
        const C0 = this.calculate_capacitance(true);

        // 2. Calculate C (with dielectric capacitance)
        currentStep = 2;
        await this.solve_laplace_iterative(false, (i, max) => updateProgress(i / max, currentStep));
        const C_with_diel = this.calculate_capacitance(false);

        // 3. Calculate Z0 and effective permittivity
        const eps_eff = C_with_diel / C0;
        const Z0 = 1 / (CONSTANTS.C * Math.sqrt(C_with_diel * C0));

        // 4. Compute fields Ex, Ey
        this.compute_fields(); // This will set this.Ex and this.Ey

        // 5. Calculate losses
        // Conductor loss depends on Ex, Ey, Z0, omega, sigma_cond
        const alpha_cond_db_m = this.calculate_conductor_loss(this.Ex, this.Ey, Z0);

        // Dielectric loss depends on Ex, Ey, Z0, omega, epsilon_r, tan_delta
        const alpha_diel_db_m = this.calculate_dielectric_loss(this.Ex, this.Ey, Z0);

        // 6. Calculate RLGC and complex Z0
        const { Zc, rlgc, eps_eff_mode } = this.rlgc(alpha_cond_db_m, alpha_diel_db_m, C_with_diel, Z0);

        const total_alpha_db_m = alpha_cond_db_m + alpha_diel_db_m;

        return {
            Z0: Z0, // Characteristic Impedance (real part approximation)
            Zc: Zc, // Complex Characteristic Impedance
            eps_eff: eps_eff,
            RLGC: rlgc,
            alpha_cond_db_m: alpha_cond_db_m,
            alpha_diel_db_m: alpha_diel_db_m,
            total_alpha_db_m: total_alpha_db_m
        };
    }
}
