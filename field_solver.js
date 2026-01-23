import createWASMModule from './wasm_solver/solver.js';
import { Complex } from "./complex.js";
import { calculate_Zrough } from './surface_roughness.js';

export const CONSTANTS = {
    EPS0: 8.854187817e-12,
    MU0: 4 * Math.PI * 1e-7,
    C: 299792458,
    PI: Math.PI
};

// --- Math Utils ---

export function diff(arr) {
    const res = new Float64Array(arr.length - 1);
    for (let i = 0; i < arr.length - 1; i++) res[i] = arr[i+1] - arr[i];
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

function isArrayLike2D(arr, ny, nx) {
    if (!arr || typeof arr !== "object") return false;
    if (arr.length !== ny) return false;

    for (let i = 0; i < ny; i++) {
        const row = arr[i];
        if (!row || typeof row !== "object") return false;
        if (typeof row.length !== "number") return false;
        if (row.length !== nx) return false;
    }
    return true;
}

function validate_laplace_inputs(V, x, y, epsilon_r, conductor_mask, vacuum = false) {
    const errors = [];

    const ny = y.length;
    const nx = x.length;

    // --- Shape checks ---
    if (!isArrayLike2D(V, ny, nx)) {
        errors.push("V must be a (ny, nx) 2D array matching mesh dimensions")
    }

    if (!isArrayLike2D(conductor_mask, ny, nx)) {
        errors.push("conductor_mask must be a (ny, nx) 2D array matching mesh dimensions")
    }

    if (!vacuum && !isArrayLike2D(epsilon_r, ny, nx)) {
        errors.push("epsilon_r must be a (ny, nx) 2D array matching mesh dimensions when vacuum=false")
    }

    // --- dx/dy checks ---
    const dx = diff(x);
    const dy = diff(y);

    const check_spacing = (d, name) => {
        const min = Math.min(...d);
        const max = Math.max(...d);

        if (Number.isNaN(min)) {
            errors.push(`NaN in ${name})`);
        }
        if (!(min > 1e-15)) {
            errors.push(`${name}: min spacing <= 1e-15 (min=${min})`);
        }
        if (!(max / min < 1e12)) {
            errors.push(`${name}: spacing ratio too large (max/min = ${max / min})`);
        }
    };

    if (dx.length > 0) check_spacing(dx, "dx");
    if (dy.length > 0) check_spacing(dy, "dy");

    // --- conductor presence ---
    let has_conductor = false;
    for (let i = 0; i < ny && !has_conductor; i++) {
        for (let j = 0; j < nx; j++) {
            if (conductor_mask[i][j]) {
                has_conductor = true;
                break;
            }
        }
    }
    if (!has_conductor) {
        errors.push("No conductor cells found in conductor_mask");
    }

    // --- V validity ---
    for (let i = 0; i < ny; i++) {
        for (let j = 0; j < nx; j++) {
            const v = V[i][j];
            if (!Number.isFinite(v)) {
                errors.push(`V contains non-finite value at (${i}, ${j}): ${v}`);
                break;
            }
        }
    }

    // --- epsilon_r validity ---
    if (!vacuum) {
        for (let i = 0; i < ny; i++) {
            for (let j = 0; j < nx; j++) {
                const er = epsilon_r[i][j];
                if (!Number.isFinite(er)) {
                    errors.push(`epsilon_r contains non-finite value at (${i}, ${j}): ${er}`);
                    return errors;
                }
                if (!(er > 0)) {
                    errors.push(`epsilon_r must be > 0 at (${i}, ${j}), got ${er}`);
                    return errors;
                }
            }
        }
    }

    return errors;
}

// --- Solver Class ---
export class FieldSolver2D {
    constructor() {
        this.x = null;
        this.y = null;
        this.V = null;  // Stored as array: [V] for single-ended, [V_odd, V_even] for differential
        this.epsilon_r = null;
        this.tand = null;
        this.conductor_mask = null; // 1 if conductor, 0 if dielectric
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
    async solve_laplace(V, vacuum = false, onProgress = null) {
        // Ensure mesh is generated
        if (this.ensure_mesh) {
            this.ensure_mesh();
        }

        const errors = validate_laplace_inputs(
                V,
                this.x,
                this.y,
                this.epsilon_r,
                this.conductor_mask,
                vacuum
            );

        if (errors.length > 0) {
                throw new Error(
                    "Laplace solver input validation failed:\n" +
                    errors.map(e => " - " + e).join("\n")
                );
            }

        const ny = this.y.length, nx = this.x.length;
        const dx = diff(this.x), dy = diff(this.y);
        const N = nx * ny;
        const idx = (i, j) => i * nx + j;

        const get_er = (i, j) => vacuum ? 1.0 : this.epsilon_r[i][j];
        const is_cond = (i, j) => this.conductor_mask[i][j];

        // Remove mesh nodes internal to conductors
        // E-field inside conductors is known to be 0.
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

        // Build sparse system
        const B = new Float64Array(N_unknown);
        const diag = new Float64Array(N_unknown);

        const colLists = Array(N_unknown);
        const valLists = Array(N_unknown);

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

        const { rowPtr, colIdx, values } = buildCSR(colLists, valLists, N_unknown);
        const csr = { rowPtr, colIdx, values };
        const x = await solveWithWASM(csr, B, true);

        // Reconstruct solution for full mesh
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

    /**
     * Calculate conductor cross-sectional area from conductor dimensions.
     * Uses the Conductor class dimensions directly (width * height) rather than
     * summing mesh elements for accurate DC resistance calculation.
     *
     * For differential mode, includes both signal traces in signal_area.
     * Ground area includes all ground conductors (bottom, top, sides, vias).
     *
     * @returns {{signal_area: number, ground_area: number}} - Cross-sectional areas in m^2
     */
    _calculate_conductor_area() {
        if (!this.conductors) {
            throw new Error("Conductors array not available");
        }

        let signal_area = 0;
        let ground_area = 0;

        for (const cond of this.conductors) {
            const area = Math.abs(cond.width * cond.height);
            if (cond.is_signal) {
                signal_area += area;
            } else {
                ground_area += area;
            }
        }

        return { signal_area, ground_area };
    }

    /**
     * Calculate conductor losses including both DC and AC (skin effect) contributions.
     *
     * The total resistance is calculated as R_total = sqrt(R_dc^2 + R_ac^2) where:
     * - R_dc: DC resistance from conductor cross-sectional area (uses Conductor dimensions directly)
     * - R_ac: AC resistance from skin effect and surface roughness
     *
     * For transmission lines:
     * - Current flows through signal conductor and returns through ground
     *
     * For differential traces:
     * - Both signal traces contribute to signal_area
     * - Power normalization factor of 0.5 is applied to AC components
     *
     * @param {Array<Array<number>>} Ex - Electric field x-component
     * @param {Array<Array<number>>} Ey - Electric field y-component
     * @param {number} Z0 - Characteristic impedance (real part)
     * @returns {{R_ac: number, R_dc: number, R_total: number, L_internal: number}}
     */
    calculate_conductor_loss(Ex, Ey, Z0) {
        if (!this.solution_valid) throw new Error("Fields invalid");

        // Calculate DC resistance from conductor dimensions
        const { signal_area, ground_area } = this._calculate_conductor_area();

        // DC resistance per unit length for transmission line
        // Current flows through signal conductor and returns through ground (series connection)
        // R_dc = R_signal + R_ground = 1/(σ×A_signal) + 1/(σ×A_ground)
        const R_signal = 1.0 / (this.sigma_cond * signal_area);
        const R_ground = 1.0 / (this.sigma_cond * ground_area);
        const R_dc = R_signal + R_ground;

        // Handle DC case (frequency = 0)
        if (this.freq === 0 || this.omega === 0) {
            return {
                R_ac: 0,
                R_dc: R_dc,
                R_total: R_dc,
                L_internal: 0
            };
        }

        // Use roughness from constructor
        const rq = this.rq || 0;

        // Pre-calculate Surface Impedance (Z_rough) using gradient model
        // Returns Complex object {re, im}
        const Z_surf = calculate_Zrough(this.freq, this.sigma_cond, rq);

        const ny = this.y.length;
        const nx = this.x.length;
        const dx_array = diff(this.x);
        const dy_array = diff(this.y);

        const get_dx = j => (j >= 0 && j < dx_array.length) ? dx_array[j] : dx_array[dx_array.length - 1];
        const get_dy = i => (i >= 0 && i < dy_array.length) ? dy_array[i] : dy_array[dy_array.length - 1];

        let sum_H2_dl_R = 0.0; // Sum for Resistance
        let sum_H2_dl_L = 0.0; // Sum for Inductance

        // Helper to check mask type
        const isSignal = (i, j) => this.signal_mask[i][j];
        const isGround = (i, j) => this.ground_mask[i][j];
        const isConductor = (i, j) => isSignal(i,j) || isGround(i,j);

        for (let i = 1; i < ny - 1; i++) {
            for (let j = 1; j < nx - 1; j++) {

                // We only care about the boundary inside the dielectric
                // but adjacent to a conductor.
                // Current loop iterates all cells. Look for Dielectric cells
                // that have a conductor neighbor.

                if (isConductor(i, j)) continue; // Skip if we are inside metal

                // Check neighbors to see if we are on a boundary
                const neighbors = [
                    { ni: i, nj: j + 1, direction: 'r', dl_func: get_dy, idx: i },
                    { ni: i, nj: j - 1, direction: 'l', dl_func: get_dy, idx: i },
                    { ni: i + 1, nj: j, direction: 'u', dl_func: get_dx, idx: j },
                    { ni: i - 1, nj: j, direction: 'd', dl_func: get_dx, idx: j },
                ];

                for (const { ni, nj, direction, dl_func, idx: dl_idx } of neighbors) {
                    if (ni < 0 || ni >= ny || nj < 0 || nj >= nx) continue;

                    // If neighbor is a conductor, we are on the surface
                    if (isConductor(ni, nj)) {

                        const eps_diel = this.epsilon_r[i][j]; // Use permittivity of current (dielectric) cell
                        const Ex_val = Ex[i][j];
                        const Ey_val = Ey[i][j];

                        // Calculate tangential Magnetic Field H_tan
                        // H_tan = (1/Z_wave) * E_norm
                        // Z_wave = Z0 / sqrt(eps_r)

                        let E_norm = 0.0;
                        if (direction === 'r' || direction === 'l') E_norm = Math.abs(Ex_val);
                        else E_norm = Math.abs(Ey_val);

                        const Z0_freespace = 376.73; // sqrt(mu0/eps0)
                        const H_tan = E_norm * Math.sqrt(eps_diel) / Z0_freespace;

                        const dl = dl_func(dl_idx);
                        const H2_dl = H_tan * H_tan * dl;

                        // Accumulate for all conductors (signal and ground use same roughness)
                        sum_H2_dl_R += Z_surf.re * H2_dl;
                        sum_H2_dl_L += Z_surf.im * H2_dl;
                    }
                }
            }
        }

        // Power normalization factor: differential has 0.5 factor
        // This is because we integrate over both traces but report normalized loss
        const power_factor = this.is_differential ? 0.5 : 1.0;

        const Z0_sq = Z0 * Z0;

        // AC Resistance per unit length from skin effect (Ohm/m)
        const R_ac = power_factor * sum_H2_dl_R * Z0_sq;

        // This doesn't hold if conductor thickness is smaller than skin depth
        // Need to solve magnetic field for accurate L_internal at low frequency
        // but is not a problem at even moderately high frequency >1 MHz.
        // In practice very minimal error since DC can be solved correctly.
        const L_internal = power_factor * sum_H2_dl_L * Z0_sq / this.omega;

        const R_total = Math.sqrt(R_dc * R_dc + R_ac * R_ac);

        return { R_ac, R_dc, R_total, L_internal };
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

        // Helper function for conductor detection based on mode
        const isConductor = this.is_differential
            ? (i, j) => this.signal_p_mask[i][j] || this.signal_n_mask[i][j] || this.ground_mask[i][j]
            : (i, j) => this.conductor_mask[i][j];

        let Pd = 0.0;

        for (let i = 0; i < ny - 1; i++) {
            for (let j = 0; j < nx - 1; j++) {
                if (isConductor(i, j)) continue;

                const E2 = Ex[i][j] * Ex[i][j] + Ey[i][j] * Ey[i][j];
                const dA = get_dx(j) * get_dy(i);

                Pd += 0.5 * this.omega * CONSTANTS.EPS0 * this.epsilon_r[i][j] * this.tand[i][j] * E2 * dA;
            }
        }

        // Power normalization: differential has 0.5 factor
        const power_factor = this.is_differential ? 0.5 : 1.0;
        const P_flow = 1.0 / (2 * Z0);
        return 8.686 * (power_factor * Pd / (2 * P_flow));
    }

    rlgc(R_total, L_internal, alpha_diel, C_mode, Z0_mode) {

        // Dielectric loss conductance
        const alpha_d_np = alpha_diel / 8.686;
        // alpha_d = G * Z0 / 2  => G = 2 * alpha_d / Z0
        const G = 2 * alpha_d_np / Z0_mode;

        // External Inductance (Geometric)
        const L_ext = (Z0_mode * Z0_mode) * C_mode;

        // Total Inductance
        const L_total = L_ext + L_internal;

        // Re-calculate complex Zc and Epsilon_eff with the new L and R
        const omega = this.omega;

        // Zc = sqrt( (R + jwL) / (G + jwC) )
        const Z_num = new Complex(R_total, omega * L_total);
        const Z_den = new Complex(G, omega * C_mode);
        const Zc = Z_num.div(Z_den).sqrt();

        // Effective Permittivity
        // gamma = sqrt( (R+jwL)(G+jwC) ) = alpha + j*beta
        // beta = Im(gamma)
        // eps_eff = (beta / k0)^2  where k0 = omega/c0
        const gamma = Z_num.mul(Z_den).sqrt();
        const beta = gamma.im;
        const k0 = omega / 299792458.0;
        const eps_eff_new = Math.pow(beta / k0, 2);

        return {
            Zc: Zc,
            rlgc: {
                R: R_total,
                L: L_total,
                G: G,
                C: C_mode
            },
            eps_eff_mode: eps_eff_new,
            L_internal: L_internal,
            L_external: L_ext
        };
    }

    rlgc_from_alpha(alpha_cond, alpha_diel, C_mode, Z0_mode) {
        // Convert dB/m to Np/m
        const alpha_c_np = alpha_cond / 8.686;
        const alpha_d_np = alpha_diel / 8.686;

        const R = 2 * Z0_mode * alpha_c_np;
        const G = 2 * alpha_d_np / Z0_mode;

        const L = (Z0_mode * Z0_mode) * C_mode;

        const eps_eff_mode = (L * C_mode) / (CONSTANTS.MU0 * CONSTANTS.EPS0);

        // Handle DC case (omega = 0)
        if (this.omega === 0) {
            if (G === 0) {
                return {
                    Zc: new Complex(Infinity, 0),
                    rlgc: { R, L, G: 0, C: C_mode },
                    eps_eff_mode
                };
            }
            return {
                Zc: new Complex(Math.sqrt(R / G), 0),
                rlgc: { R, L, G, C: C_mode },
                eps_eff_mode
            };
        }

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
        const totalSteps = 2; // Two main solve_laplace calls
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
        V = await this.solve_laplace(V, true, (i, max) => updateProgress(i / max, currentStep));
        const C0 = this.calculate_capacitance(V, true);

        // 2. Calculate C (with dielectric capacitance)
        currentStep = 2;
        V = this._create_voltage_array('single');
        V = await this.solve_laplace(V, false, (i, max) => updateProgress(i / max, currentStep));
        const C_with_diel = this.calculate_capacitance(V, false);

        // 3. Calculate Z0 and effective permittivity
        const eps_eff = C_with_diel / C0;
        const Z0 = 1 / (CONSTANTS.C * Math.sqrt(C_with_diel * C0));

        // 4. Compute fields Ex, Ey
        const { Ex, Ey } = this.compute_fields(V);

        // 5. Calculate losses
        // Conductor loss with surface roughness and DC resistance
        const { R_ac, R_dc, R_total, L_internal } = this.calculate_conductor_loss(Ex, Ey, Z0);

        // Dielectric loss depends on Ex, Ey, Z0, omega, epsilon_r, tan_delta
        const alpha_diel_db_m = this.calculate_dielectric_loss(Ex, Ey, Z0);

        // 6. Calculate RLGC and complex Z0 using surface roughness aware approach
        const { Zc, rlgc, eps_eff_mode } = this.rlgc(R_total, L_internal, alpha_diel_db_m, C_with_diel, Z0);

        // Calculate conductor loss alpha from R_total for reporting
        const alpha_cond_db_m = 8.686 * R_total / (2 * Z0);
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

    // Adaptive Meshing
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

            // Ensure symmetry by only considering the left side.
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
         * Main mesh refinement routine.
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

        if (prev_energy === null || prev_energy === undefined) {
            return { energy, rel_error: 1.0 };
        }

        const rel_error = Math.abs(energy - prev_energy) / Math.max(Math.abs(prev_energy), 1e-12);
        return { energy, rel_error };
    }

    _compute_parameter_error(Z0, C, prev_Z0, prev_C) {
        /**
         * Track convergence of Z0 and C parameters.
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
         * Solve a single mode and return full results.
         *
         * Parameters:
         * -----------
         * mode : string - 'single', 'odd', or 'even'
         * vacuum_first : boolean - Whether to solve vacuum case first for C0 calculation
         *
         * Returns:
         * --------
         * {mode, Z0, eps_eff, C, C0, RLGC, Zc, alpha_c, alpha_d, alpha_total, V, Ex, Ey}
         */
        let C0;
        let V;

        if (vacuum_first) {
            // Calculate C0 (vacuum capacitance)
            V = this._create_voltage_array(mode);
            V = await this.solve_laplace(V, true);

            if (this.is_differential) {
                // Use only positive trace for charge calculation
                const orig_signal_mask = this.signal_mask.map(row => row.slice());
                this.signal_mask = this.signal_p_mask.map(row => row.slice());
                C0 = this.calculate_capacitance(V, true);
                this.signal_mask = orig_signal_mask;
            } else {
                C0 = this.calculate_capacitance(V, true);
            }
        }

        // Solve with dielectric
        V = this._create_voltage_array(mode);
        V = await this.solve_laplace(V, false);

        let C;
        if (this.is_differential) {
            // Use only positive trace for charge calculation
            const orig_signal_mask = this.signal_mask.map(row => row.slice());
            this.signal_mask = this.signal_p_mask.map(row => row.slice());
            C = this.calculate_capacitance(V, false);
            this.signal_mask = orig_signal_mask;
        } else {
            C = this.calculate_capacitance(V, false);
        }

        // Calculate fields
        const { Ex, Ey } = this.compute_fields(V);

        // Calculate impedance
        let eps_eff, Z0;
        if (C0 !== undefined) {
            eps_eff = C / C0;
            Z0 = 1 / (CONSTANTS.C * Math.sqrt(C * C0));
        }

        // Calculate conductor losses with surface roughness and DC resistance
        const { R_ac, R_dc, R_total, L_internal } = this.calculate_conductor_loss(Ex, Ey, Z0);

        // Calculate dielectric loss (returns alpha in dB/m)
        const alpha_d = this.calculate_dielectric_loss(Ex, Ey, Z0);

        // Calculate RLGC using new surface roughness aware approach
        const { Zc, rlgc, eps_eff_mode, L_external } = this.rlgc(R_total, L_internal, alpha_d, C, Z0);

        // Calculate conductor loss alpha from R_total for reporting
        const alpha_c = 8.686 * R_total / (2 * Zc.re);
        const alpha_total = alpha_c + alpha_d;

        return {
            mode,
            Z0,
            eps_eff: eps_eff_mode,
            C, C0,
            RLGC: rlgc, Zc,
            alpha_c, alpha_d, alpha_total,
            L_internal, L_external,
            V, Ex, Ey
        };
    }

    async solve_adaptive(options = {}) {
        /**
         * Adaptive mesh solve with robust convergence criteria.
         * Automatically handles both single-ended and differential modes.
         *
         * Options:
         * --------
         * skip_mesh: boolean - If true, skip mesh refinement (use existing mesh)
         *
         * Returns:
         * --------
         * {
         *   modes: [{mode, Z0, eps_eff, C, C0, RLGC, Zc, alpha_c, alpha_d, alpha_total, V, Ex, Ey}, ...],
         *   Z_diff: (only for differential) 2 * Z_odd,
         *   Z_common: (only for differential) Z_even / 2
         * }
         */
        // Ensure mesh is generated
        if (this.ensure_mesh) {
            this.ensure_mesh();
        }

        const {
            max_iters = 10,
            refine_frac,
            energy_tol = 0.02,
            param_tol = 0.1,
            max_nodes = 20000,
            min_converged_passes = 2,
            onProgress = null,
            shouldStop = null,
            skip_mesh = false
        } = options;

        // If skip_mesh is true, just solve once with existing mesh
        if (skip_mesh) {
            const modeNames = this.is_differential ? ['odd', 'even'] : ['single'];
            const modeResults = [];
            for (const modeName of modeNames) {
                const result = await this._solve_single_mode(modeName, true);
                modeResults.push(result);
            }
            // Store fields as arrays
            this.V = modeResults.map(r => r.V);
            this.Ex = modeResults.map(r => r.Ex);
            this.Ey = modeResults.map(r => r.Ey);
            return this._build_results(modeResults);
        }

        // Set default refine_frac based on mode
        const refineFrac = refine_frac !== undefined ? refine_frac : (this.is_differential ? 0.15 : 0.2);

        // Define modes to solve
        const modeNames = this.is_differential ? ['odd', 'even'] : ['single'];

        // Tracking variables for convergence
        const prevEnergy = {};
        const prevZ0 = {};
        let converged_count = 0;
        let modeResults = null;

        for (let it = 0; it < max_iters; it++) {
            // Solve all modes
            modeResults = [];
            for (const modeName of modeNames) {
                const result = await this._solve_single_mode(modeName, true);
                modeResults.push(result);
            }

            // Compute max errors across all modes
            let max_energy_err = 0;
            let max_param_err = 0;

            for (let i = 0; i < modeNames.length; i++) {
                const modeName = modeNames[i];
                const r = modeResults[i];

                const { energy, rel_error: energy_err } = this._compute_energy_error(r.Ex, r.Ey, prevEnergy[modeName]);

                const param_err = prevZ0[modeName] !== undefined
                    ? Math.abs(r.Z0 - prevZ0[modeName]) / Math.max(Math.abs(prevZ0[modeName]), 1e-12)
                    : 1.0;

                max_energy_err = Math.max(max_energy_err, energy_err);
                max_param_err = Math.max(max_param_err, param_err);

                prevEnergy[modeName] = energy;
                prevZ0[modeName] = r.Z0;
            }

            console.log(`Pass ${it + 1}: Energy err=${max_energy_err.toExponential(3)}, Param err=${max_param_err.toExponential(3)}, Grid=${this.x.length}x${this.y.length}`);

            // Call progress callback
            if (onProgress) {
                onProgress({
                    iteration: it + 1,
                    max_iterations: max_iters,
                    energy_error: max_energy_err,
                    param_error: max_param_err,
                    nodes_x: this.x.length,
                    nodes_y: this.y.length
                });
            }

            // Yield to event loop to allow UI updates
            await new Promise(resolve => setTimeout(resolve, 0));

            // Check convergence
            const hasPrevious = Object.keys(prevZ0).length === modeNames.length &&
                                Object.values(prevZ0).every(v => v !== undefined);
            if (hasPrevious && it > 0) {
                if (max_energy_err < energy_tol && max_param_err < param_tol) {
                    converged_count++;
                    if (converged_count >= min_converged_passes) {
                        console.log(`Converged after ${it + 1} passes`);
                        break;
                    }
                } else {
                    converged_count = 0;
                }
            }

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

            // Refine mesh using first mode's fields (odd for differential, single for single-ended)
            if (it !== max_iters - 1) {
                const refineMode = modeResults[0];
                this.refine_mesh(refineMode.V, refineMode.Ex, refineMode.Ey, refineFrac);
                this._setup_geometry();
            }
        }

        // Store fields as arrays
        this.V = modeResults.map(r => r.V);
        this.Ex = modeResults.map(r => r.Ex);
        this.Ey = modeResults.map(r => r.Ey);

        // Build unified result structure
        return this._build_results(modeResults);
    }

    _build_results(modeResults) {
        /**
         * Build the unified result structure from mode results.
         */
        const result = { modes: modeResults };

        if (this.is_differential) {
            const odd = modeResults.find(m => m.mode === 'odd');
            const even = modeResults.find(m => m.mode === 'even');
            result.Z_diff = 2 * odd.Z0;
            result.Z_common = even.Z0 / 2;
        }

        return result;
    }
}
