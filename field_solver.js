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
        this.V = null;  // Stored as array: [V] for single-ended, [V_odd, V_even] for differential
        this.epsilon_r = null;
        this.conductor_mask = null; // 1 if conductor, 0 if dielectric
        this.bc_mask = null; // Boundary conditions
        this.solution_valid = false;

        // Computed fields - stored as array: [fields] for single-ended, [odd, even] for differential
        this.Ex = null;
        this.Ey = null;
    }

    /**
     * Create a voltage array based on conductor masks and solve mode.
     * @param {string} mode - 'single', 'odd', or 'even'
     * @returns {Array<Float64Array>} - 2D voltage array
     */
    _create_voltage_array(mode = 'single') {
        const ny = this.y.length;
        const nx = this.x.length;
        const V = Array(ny).fill().map(() => new Float64Array(nx));

        for (let i = 0; i < ny; i++) {
            for (let j = 0; j < nx; j++) {
                if (this.ground_mask[i][j]) {
                    V[i][j] = 0.0;
                } else if (mode === 'odd' && this.is_differential) {
                    // Odd mode: positive trace = +1V, negative trace = -1V
                    if (this.signal_p_mask[i][j]) V[i][j] = 1.0;
                    else if (this.signal_n_mask[i][j]) V[i][j] = -1.0;
                } else if (mode === 'even' && this.is_differential) {
                    // Even mode: both traces = +1V
                    if (this.signal_mask[i][j]) V[i][j] = 1.0;
                } else {
                    // Single-ended: signal = +1V
                    if (this.signal_mask[i][j]) V[i][j] = 1.0;
                }
            }
        }
        return V;
    }

    /**
     * Solve Laplace equation for the given voltage array.
     * @param {Array<Float64Array>} V - 2D voltage array with conductor boundary conditions set
     * @param {boolean} vacuum - If true, solve with vacuum permittivity
     * @param {function} onProgress - Optional progress callback
     * @returns {Array<Float64Array>} - The solved voltage array (same reference as input)
     */
    async solve_laplace_iterative(V, vacuum = false, onProgress = null) {
        // Ensure mesh is generated
        if (this.ensure_mesh) {
            this.ensure_mesh();
        }

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
                        B[n] -= c * V[ii][jj];
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

        // --- Reconstruct solution ---
        for (let k = 0; k < N_unknown; k++) {
            const n = red_to_full[k];
            const i = (n / nx) | 0;
            const j = n % nx;
            V[i][j] = x[k];
        }

        return V;
    }

    /**
     * Compute E-field from voltage distribution.
     * @param {Array<Float64Array>} V - 2D voltage array
     * @returns {{Ex: Array<Float64Array>, Ey: Array<Float64Array>}} - E-field components
     */
    compute_fields(V) {
        const ny = this.y.length;
        const nx = this.x.length;
        const dx = diff(this.x);
        const dy = diff(this.y);

        const Ex = Array(ny).fill().map(() => new Float64Array(nx));
        const Ey = Array(ny).fill().map(() => new Float64Array(nx));

        for(let i=1; i<ny-1; i++) {
            for(let j=1; j<nx-1; j++) {
                if (this.conductor_mask[i][j]) continue;

                const dxl = dx[j-1];
                const dxr = dx[j];
                const dyd = dy[i-1];
                const dyu = dy[i];

                Ex[i][j] = -(
                    (dxl / (dxr * (dxl + dxr))) * V[i][j+1] +
                    ((dxr - dxl) / (dxl * dxr)) * V[i][j] -
                    (dxr / (dxl * (dxl + dxr))) * V[i][j-1]
                );

                Ey[i][j] = -(
                    (dyd / (dyu * (dyd + dyu))) * V[i+1][j] +
                    ((dyu - dyd) / (dyd * dyu)) * V[i][j] -
                    (dyu / (dyd * (dyd + dyu))) * V[i-1][j]
                );
            }
        }
        this.solution_valid = true;
        return { Ex, Ey };
    }

    /**
     * Calculate capacitance from voltage distribution.
     * @param {Array<Float64Array>} V - 2D voltage array
     * @param {boolean} vacuum - If true, use vacuum permittivity
     * @returns {number} - Capacitance in F/m
     */
    calculate_capacitance(V, vacuum=false) {
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
                         En = (V[i][j] - V[ni][nj]) / dist;
                         // Average dx for area
                         area = (get_dx(j-1) + get_dx(j)) / 2;
                    } else {
                        // Neighbor is Left/Right
                        dist = Math.abs(this.x[j] - this.x[nj]);
                        En = (V[i][j] - V[ni][nj]) / dist;
                        // Average dy for area
                        area = (get_dy(i-1) + get_dy(i)) / 2;
                    }

                    const er = vacuum ? 1 : this.epsilon_r[ni][nj];
                    Q += CONSTANTS.EPS0 * er * En * area;
                };

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
        let V = this._create_voltage_array('single');
        V = await this.solve_laplace_iterative(V, true, (i, max) => updateProgress(i / max, currentStep));
        const C0 = this.calculate_capacitance(V, true);

        // 2. Calculate C (with dielectric capacitance)
        currentStep = 2;
        V = this._create_voltage_array('single');
        V = await this.solve_laplace_iterative(V, false, (i, max) => updateProgress(i / max, currentStep));
        const C_with_diel = this.calculate_capacitance(V, false);

        // 3. Calculate Z0 and effective permittivity
        const eps_eff = C_with_diel / C0;
        const Z0 = 1 / (CONSTANTS.C * Math.sqrt(C_with_diel * C0));

        // 4. Compute fields Ex, Ey
        const { Ex, Ey } = this.compute_fields(V);

        // 5. Calculate losses
        // Conductor loss depends on Ex, Ey, Z0, omega, sigma_cond
        const alpha_cond_db_m = this.calculate_conductor_loss(Ex, Ey, Z0);

        // Dielectric loss depends on Ex, Ey, Z0, omega, epsilon_r, tan_delta
        const alpha_diel_db_m = this.calculate_dielectric_loss(Ex, Ey, Z0);

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
            total_alpha_db_m: total_alpha_db_m,
            V: V,  // Return V for storage
            Ex: Ex,
            Ey: Ey
        };
    }

    // ============================================================================
    // ADAPTIVE MESHING METHODS
    // ============================================================================

    _compute_refine_metrics(V, Ex, Ey) {
        /**
         * For each grid interval, compute a metric indicating how much refinement
         * would help, based on voltage gradients and field energy in adjacent cells.
         */
        const ny = V.length;
        const nx = V[0].length;

        // Metric for splitting interval [x[j], x[j+1]]
        const x_metrics = new Float64Array(this.x.length - 1);
        // Metric for splitting interval [y[i], y[i+1]]
        const y_metrics = new Float64Array(this.y.length - 1);

        for (let i = 0; i < ny - 1; i++) {
            for (let j = 0; j < nx - 1; j++) {
                // Skip cells fully inside conductors
                if (this.conductor_mask[i][j] &&
                    this.conductor_mask[Math.min(i + 1, ny - 1)][j] &&
                    this.conductor_mask[i][Math.min(j + 1, nx - 1)]) {
                    continue;
                }

                const eps = this.epsilon_r[i][j];

                // Voltage differences across this cell
                const dV_x = j < nx - 1 ? Math.abs(V[i][j + 1] - V[i][j]) : 0;
                const dV_y = i < ny - 1 ? Math.abs(V[i + 1][j] - V[i][j]) : 0;

                // Field magnitude for weighting
                const E2 = Ex[i][j] ** 2 + Ey[i][j] ** 2;
                const E_mag = E2 > 0 ? Math.sqrt(E2) : 1e-12;

                // Boundary detection
                const is_boundary = (!this.conductor_mask[i][j] && (
                    (i > 0 && this.conductor_mask[i - 1][j]) ||
                    (i < ny - 1 && this.conductor_mask[i + 1][j]) ||
                    (j > 0 && this.conductor_mask[i][j - 1]) ||
                    (j < nx - 1 && this.conductor_mask[i][j + 1])));
                const boundary_mult = is_boundary ? 2.0 : 1.0;

                // Weight by field strength, permittivity, and boundary importance
                const weight = E_mag * eps * boundary_mult;

                // Accumulate to the interval metrics
                if (j < x_metrics.length) {
                    x_metrics[j] += dV_x * weight;
                }
                if (i < y_metrics.length) {
                    y_metrics[i] += dV_y * weight;
                }
            }
        }

        return { x_metrics, y_metrics };
    }

    _check_symmetry(coords, center, tol = 1e-10) {
        /**
         * Check if coordinate array is symmetric about center.
         */
        const n = coords.length;
        for (let k = 0; k < Math.floor(n / 2); k++) {
            const left = coords[k];
            const right = coords[n - 1 - k];
            if (Math.abs((left - center) + (right - center)) > tol) {
                return false;
            }
        }
        return true;
    }

    _symmetrize_metrics(metrics) {
        /**
         * Average metrics for symmetric pairs.
         */
        const n = metrics.length;
        const result = new Float64Array(n);
        for (let k = 0; k < n; k++) {
            result[k] = metrics[k];
        }
        for (let k = 0; k < Math.floor(n / 2); k++) {
            const avg = 0.5 * (metrics[k] + metrics[n - 1 - k]);
            result[k] = avg;
            result[n - 1 - k] = avg;
        }
        return result;
    }

    _select_lines_to_refine(x_metrics, y_metrics, frac = 0.15) {
        /**
         * Select which grid intervals to split, respecting left-right symmetry.
         */
        const x_center = (this.x[0] + this.x[this.x.length - 1]) / 2;
        const x_symmetric = this._check_symmetry(this.x, x_center);

        let x_metrics_proc = x_metrics;
        if (x_symmetric) {
            x_metrics_proc = this._symmetrize_metrics(x_metrics);
        }

        // Decide how many x vs y lines based on relative total metric
        let total_x = 0;
        let total_y = 0;
        for (let i = 0; i < x_metrics_proc.length; i++) total_x += x_metrics_proc[i];
        for (let i = 0; i < y_metrics.length; i++) total_y += y_metrics[i];
        const total = total_x + total_y;

        if (total < 1e-15) {
            return { selected_x: new Set(), selected_y: new Set() };
        }

        let n_total = Math.floor(frac * (x_metrics_proc.length + y_metrics.length));
        n_total = Math.max(1, n_total);

        // Allocate proportionally to where the error is
        const n_x = Math.floor(n_total * total_x / total);
        const n_y = n_total - n_x;

        // Select top intervals
        const x_ranked = Array.from(x_metrics_proc.keys()).sort((a, b) => x_metrics_proc[b] - x_metrics_proc[a]);
        const y_ranked = Array.from(y_metrics.keys()).sort((a, b) => y_metrics[b] - y_metrics[a]);

        const selected_x = new Set();
        const selected_y = new Set();

        for (let idx = 0; idx < Math.min(n_x, x_ranked.length); idx++) {
            const j = x_ranked[idx];
            if (x_metrics_proc[j] > 0) {
                selected_x.add(j);
                if (x_symmetric) {
                    const partner = x_metrics_proc.length - 1 - j;
                    if (partner >= 0 && partner < x_metrics_proc.length) {
                        selected_x.add(partner);
                    }
                }
            }
        }

        for (let idx = 0; idx < Math.min(n_y, y_ranked.length); idx++) {
            const i = y_ranked[idx];
            if (y_metrics[i] > 0) {
                selected_y.add(i);
            }
        }

        return { selected_x, selected_y };
    }

    _refine_selected_lines(selected_x, selected_y) {
        /**
         * Add new grid lines at midpoints of selected intervals.
         */
        const x_center = (this.x[0] + this.x[this.x.length - 1]) / 2;
        const x_symmetric = this._check_symmetry(this.x, x_center);

        const new_x = new Set();
        const new_y = new Set();

        for (const j of selected_x) {
            const midpoint = 0.5 * (this.x[j] + this.x[j + 1]);

            if (x_symmetric) {
                if (midpoint <= x_center) {
                    new_x.add(midpoint);
                    const symmetric_point = 2 * x_center - midpoint;
                    if (symmetric_point > this.x[0] && symmetric_point < this.x[this.x.length - 1]) {
                        new_x.add(symmetric_point);
                    }
                }
            } else {
                new_x.add(midpoint);
            }
        }

        for (const i of selected_y) {
            const midpoint = 0.5 * (this.y[i] + this.y[i + 1]);
            new_y.add(midpoint);
        }

        // Merge and sort
        const all_x = new Set([...this.x, ...new_x]);
        const all_y = new Set([...this.y, ...new_y]);

        this.x = Float64Array.from([...all_x].sort((a, b) => a - b));
        this.y = Float64Array.from([...all_y].sort((a, b) => a - b));
    }

    refine_mesh(V, Ex, Ey, frac = 0.15) {
        /**
         * Main refinement routine.
         */
        const { x_metrics, y_metrics } = this._compute_refine_metrics(V, Ex, Ey);
        const { selected_x, selected_y } = this._select_lines_to_refine(x_metrics, y_metrics, frac);
        this._refine_selected_lines(selected_x, selected_y);

        // Invalidate solution since mesh has changed
        this.solution_valid = false;
        this.Ex = null;
        this.Ey = null;
    }

    _compute_energy_error(Ex, Ey, prev_energy) {
        /**
         * Compute relative change in stored electromagnetic energy.
         * Fields must be interpolated to cell centers to match epsilon grid.
         */
        const ny = this.y.length;
        const nx = this.x.length;
        const dx_array = diff(this.x);
        const dy_array = diff(this.y);

        let energy = 0.0;
        for (let i = 0; i < ny - 1; i++) {
            for (let j = 0; j < nx - 1; j++) {
                if (this.conductor_mask[i][j]) {
                    continue;
                }

                const E2 = Ex[i][j] ** 2 + Ey[i][j] ** 2;
                const dA = dx_array[j] * dy_array[i];
                energy += 0.5 * CONSTANTS.EPS0 * this.epsilon_r[i][j] * E2 * dA;
            }
        }

        if (prev_energy === null) {
            return { energy, rel_error: 1.0 };
        }

        const rel_error = Math.abs(energy - prev_energy) / Math.max(Math.abs(prev_energy), 1e-12);
        return { energy, rel_error };
    }

    _compute_parameter_error(Z0, C, prev_Z0, prev_C) {
        /**
         * Track convergence of the quantities you actually care about.
         */
        if (prev_Z0 === null) {
            return 1.0;
        }

        const z_err = Math.abs(Z0 - prev_Z0) / Math.max(Math.abs(prev_Z0), 1e-12);
        const c_err = Math.abs(C - prev_C) / Math.max(Math.abs(prev_C), 1e-12);
        return Math.max(z_err, c_err);
    }

    async _solve_single_mode(mode, vacuum_first = true) {
        /**
         * Solve a single mode (odd or even) for differential pair.
         *
         * Parameters:
         * -----------
         * mode : string - 'odd' or 'even'
         * vacuum_first : boolean - Whether to solve vacuum case first for C0 calculation
         *
         * Returns:
         * --------
         * {Z0, eps_eff, C, C0, Ex, Ey, V}
         */
        let C0;
        let V;

        if (vacuum_first) {
            // Calculate C0 (vacuum capacitance)
            V = this._create_voltage_array(mode);
            V = await this.solve_laplace_iterative(V, true);

            // Use only positive trace for charge calculation
            const orig_signal_mask = this.signal_mask.map(row => row.slice());
            this.signal_mask = this.signal_p_mask.map(row => row.slice());
            C0 = this.calculate_capacitance(V, true);
            this.signal_mask = orig_signal_mask;
        }

        // Solve with dielectric
        V = this._create_voltage_array(mode);
        V = await this.solve_laplace_iterative(V, false);

        // Use only positive trace for charge calculation
        const orig_signal_mask = this.signal_mask.map(row => row.slice());
        this.signal_mask = this.signal_p_mask.map(row => row.slice());
        const C = this.calculate_capacitance(V, false);
        this.signal_mask = orig_signal_mask;

        // Calculate fields
        const { Ex, Ey } = this.compute_fields(V);

        // Calculate impedance
        let eps_eff, Z0;
        if (C0 !== undefined) {
            eps_eff = C / C0;
            Z0 = 1 / (CONSTANTS.C * Math.sqrt(C * C0));
        }

        return { Z0, eps_eff, C, C0, Ex, Ey, V };
    }

    _calculate_differential_conductor_loss(Ex, Ey, Z0, mode) {
        /**
         * Calculate conductor loss for differential pair.
         */
        const Rs = Math.sqrt(this.omega * CONSTANTS.MU0 / (2 * this.sigma_cond));
        const delta = Math.sqrt(2 / (this.omega * CONSTANTS.MU0 * this.sigma_cond));

        const nx = this.x.length;
        const ny = this.y.length;
        const dx_array = diff(this.x);
        const dy_array = diff(this.y);

        const get_dx = j => (j >= 0 && j < dx_array.length) ? dx_array[j] : dx_array[dx_array.length - 1];
        const get_dy = i => (i >= 0 && i < dy_array.length) ? dy_array[i] : dy_array[dy_array.length - 1];

        let Pc = 0.0;

        for (let i = 1; i < ny - 1; i++) {
            for (let j = 1; j < nx - 1; j++) {
                if (!this.signal_p_mask[i][j] && !this.signal_n_mask[i][j] && !this.ground_mask[i][j]) {
                    continue;
                }

                const neighbors = [
                    [i, j + 1, 'r', i, get_dy],
                    [i, j - 1, 'l', i, get_dy],
                    [i + 1, j, 'u', j, get_dx],
                    [i - 1, j, 'd', j, get_dx]
                ];

                let cell_K_sq = 0.0;
                let cell_dl = 0.0;

                for (const [ni, nj, direction, dl_idx, dl_func] of neighbors) {
                    if (this.signal_p_mask[ni][nj] || this.signal_n_mask[ni][nj] || this.ground_mask[ni][nj]) {
                        continue;
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

        // Power normalization for 1W
        const Pc_1W = 0.5 * Pc * 2 * Z0;
        const alpha_db_per_m = 8.686 * Pc_1W / 2.0;

        return alpha_db_per_m;
    }

    _calculate_differential_dielectric_loss(Ex, Ey, Z0) {
        /**
         * Calculate dielectric loss for differential pair.
         */
        const nx = this.x.length;
        const ny = this.y.length;
        const dx_array = diff(this.x);
        const dy_array = diff(this.y);

        let Pd = 0.0;

        for (let i = 0; i < ny - 1; i++) {
            for (let j = 0; j < nx - 1; j++) {
                if (this.signal_p_mask[i][j] || this.signal_n_mask[i][j] || this.ground_mask[i][j]) {
                    continue;
                }
                if (this.epsilon_r[i][j] <= 1.01) {
                    continue;
                }

                const E2 = Ex[i][j] * Ex[i][j] + Ey[i][j] * Ey[i][j];
                const dA = dx_array[j] * dy_array[i];
                Pd += 0.5 * this.omega * CONSTANTS.EPS0 * this.epsilon_r[i][j] * this.tan_delta * E2 * dA;
            }
        }

        // Power flow for 1W in the transmission line
        const P_flow = 1.0 / (2 * Z0);

        // Attenuation constant
        return 8.686 * (0.5 * Pd / (2 * P_flow));
    }

    async solve_adaptive(options = {}) {
        /**
         * Adaptive mesh solve with robust convergence criteria.
         * Automatically handles both single-ended and differential modes.
         */
        // Ensure mesh is generated
        if (this.ensure_mesh) {
            this.ensure_mesh();
        }

        const {
            max_iters = 10,
            refine_frac,
            energy_tol = 0.05,
            param_tol = 0.1,
            max_nodes = 20000,
            min_converged_passes = 2,
            onProgress = null,
            shouldStop = null
        } = options;

        // Set default refine_frac based on mode
        const refineFrac = refine_frac !== undefined ? refine_frac : (this.is_differential ? 0.15 : 0.2);

        let converged_count = 0;
        let results = null;

        if (this.is_differential) {
            // ============================================================================
            // DIFFERENTIAL MODE
            // ============================================================================
            let prev_energy_odd = null;
            let prev_energy_even = null;
            let prev_Z_odd = null;
            let prev_Z_even = null;

            let Z_odd, eps_eff_odd, C_odd, C0_odd, Ex_odd, Ey_odd, V_odd;
            let Z_even, eps_eff_even, C_even, C0_even, Ex_even, Ey_even, V_even;

            for (let it = 0; it < max_iters; it++) {
                // Solve even mode (+1V and +1V)
                const even_results = await this._solve_single_mode('even', true);
                Z_even = even_results.Z0;
                eps_eff_even = even_results.eps_eff;
                C_even = even_results.C;
                C0_even = even_results.C0;
                Ex_even = even_results.Ex;
                Ey_even = even_results.Ey;
                V_even = even_results.V;

                // Solve odd mode (+1V and -1V)
                const odd_results = await this._solve_single_mode('odd', true);
                Z_odd = odd_results.Z0;
                eps_eff_odd = odd_results.eps_eff;
                C_odd = odd_results.C;
                C0_odd = odd_results.C0;
                Ex_odd = odd_results.Ex;
                Ey_odd = odd_results.Ey;
                V_odd = odd_results.V;

                // Energy-based error
                const { energy: energy_odd, rel_error: energy_err_odd } =
                    this._compute_energy_error(Ex_odd, Ey_odd, prev_energy_odd);
                const { energy: energy_even, rel_error: energy_err_even } =
                    this._compute_energy_error(Ex_even, Ey_even, prev_energy_even);
                const energy_err = Math.max(energy_err_odd, energy_err_even);

                // Parameter-based error
                let param_err;
                if (prev_Z_odd !== null) {
                    const z_odd_err = Math.abs(Z_odd - prev_Z_odd) / Math.max(Math.abs(prev_Z_odd), 1e-12);
                    const z_even_err = Math.abs(Z_even - prev_Z_even) / Math.max(Math.abs(prev_Z_even), 1e-12);
                    param_err = Math.max(z_odd_err, z_even_err);
                } else {
                    param_err = 1.0;
                }

                console.log(`Pass ${it + 1}: Energy err=${energy_err.toExponential(3)}, Param err=${param_err.toExponential(3)}, Grid=${this.x.length}x${this.y.length}`);

                // Call progress callback and yield to event loop
                if (onProgress) {
                    onProgress({
                        iteration: it + 1,
                        max_iterations: max_iters,
                        energy_error: energy_err,
                        param_error: param_err,
                        nodes_x: this.x.length,
                        nodes_y: this.y.length
                    });
                }

                // Yield to event loop to allow UI updates
                await new Promise(resolve => setTimeout(resolve, 0));

                // Check convergence
                if (prev_Z_odd !== null) {
                    if (energy_err < energy_tol && param_err < param_tol) {
                        converged_count++;
                        if (converged_count >= min_converged_passes) {
                            console.log(`Converged after ${it + 1} passes`);
                            break;
                        }
                    } else {
                        converged_count = 0;
                    }
                }

                prev_energy_odd = energy_odd;
                prev_energy_even = energy_even;
                prev_Z_odd = Z_odd;
                prev_Z_even = Z_even;

                // Node budget check
                if (this.x.length * this.y.length > max_nodes) {
                    console.log("Node budget reached");
                    break;
                }

                // Check if stop was requested
                if (shouldStop && shouldStop()) {
                    console.log("Adaptive solve stopped by user");
                    break;
                }

                // Refine mesh (use odd mode fields for refinement)
                if (it !== max_iters - 1) {
                    this.refine_mesh(V_odd, Ex_odd, Ey_odd, refineFrac);
                    this._setup_geometry();
                }
            }

            // Calculate losses
            const alpha_c_odd = this._calculate_differential_conductor_loss(Ex_odd, Ey_odd, Z_odd, 'odd');
            const alpha_d_odd = this._calculate_differential_dielectric_loss(Ex_odd, Ey_odd, Z_odd);

            const alpha_c_even = this._calculate_differential_conductor_loss(Ex_even, Ey_even, Z_even, 'even');
            const alpha_d_even = this._calculate_differential_dielectric_loss(Ex_even, Ey_even, Z_even);

            // Store fields and potentials as arrays for differential mode
            this.V = [V_odd, V_even];
            this.Ex = [Ex_odd, Ex_even];
            this.Ey = [Ey_odd, Ey_even];

            // Differential and common mode impedances
            const Z_diff = 2 * Z_odd;
            const Z_common = Z_even / 2;

            results = {
                Z_odd: Z_odd,
                Z_even: Z_even,
                Z_diff: Z_diff,
                Z_common: Z_common,
                eps_eff_odd: eps_eff_odd,
                eps_eff_even: eps_eff_even,
                C_odd: C_odd,
                C_even: C_even,
                alpha_c_odd: alpha_c_odd,
                alpha_c_even: alpha_c_even,
                alpha_d_odd: alpha_d_odd,
                alpha_d_even: alpha_d_even,
                alpha_total_odd: alpha_c_odd + alpha_d_odd,
                alpha_total_even: alpha_c_even + alpha_d_even
            };

        } else {
            // ============================================================================
            // SINGLE-ENDED MODE
            // ============================================================================
            let prev_energy = null;
            let prev_Z0 = null;
            let prev_C = null;
            let V, Ex, Ey;

            for (let it = 0; it < max_iters; it++) {
                // Calculate parameters
                results = await this.perform_analysis();
                const Z0 = results.Z0;
                const eps_eff = results.eps_eff;
                const C = results.RLGC.C;

                // Get fields and potential from perform_analysis
                V = results.V;
                Ex = results.Ex;
                Ey = results.Ey;

                // Energy-based error (primary criterion)
                const { energy, rel_error: energy_err } = this._compute_energy_error(Ex, Ey, prev_energy);

                // Parameter-based error (what you actually care about)
                const param_err = this._compute_parameter_error(Z0, C, prev_Z0, prev_C);

                const nodes = this.x.length * this.y.length;

                console.log(`Pass ${it + 1}: Energy err=${energy_err.toExponential(3)}, Param err=${param_err.toExponential(3)}, Grid=${this.x.length}x${this.y.length}`);

                // Call progress callback
                if (onProgress) {
                    onProgress({
                        iteration: it + 1,
                        total_iterations: max_iters,
                        energy_error: energy_err,
                        param_error: param_err,
                        nodes_x: this.x.length,
                        nodes_y: this.y.length,
                        Z0: Z0,
                        eps_eff: eps_eff
                    });
                }

                // Yield control periodically
                await new Promise(r => setTimeout(r, 0));

                // Check convergence (both criteria must be met)
                if (prev_energy !== null) {
                    if (energy_err < energy_tol && param_err < param_tol) {
                        converged_count++;
                        if (converged_count >= min_converged_passes) {
                            console.log(`Converged after ${it + 1} passes`);
                            break;
                        }
                    } else {
                        converged_count = 0; // Reset if we diverge
                    }
                }

                prev_energy = energy;
                prev_Z0 = Z0;
                prev_C = C;

                // Node budget check
                if (nodes > max_nodes) {
                    console.log("Node budget reached");
                    break;
                }

                // Check if stop was requested
                if (shouldStop && shouldStop()) {
                    console.log("Adaptive solve stopped by user");
                    break;
                }

                // Refine mesh (skip on last iteration)
                if (it !== max_iters - 1) {
                    this.refine_mesh(V, Ex, Ey, refineFrac);
                    this._setup_geometry();
                }
            }

            // Wrap fields and potential in array for consistency with differential mode
            this.V = [V];
            this.Ex = [Ex];
            this.Ey = [Ey];
        }

        return results;
    }
}
