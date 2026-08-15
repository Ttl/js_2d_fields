import createWASMModule from './wasm_solver/solver.js';
import { Complex } from "./complex.js";
import { calculate_Zrough, calculate_Zrough_layered } from './surface_roughness.js';
import { applyDjordjevicSarkar } from './djordjevic_sarkar.js';
import { classifyModalDecomposition } from './geometry_symmetry.js';
import { buildPhysicalRLGC } from './sparameters.js';

export const CONSTANTS = {
    EPS0: 8.854187817e-12,
    MU0: 4 * Math.PI * 1e-7,
    C: 299792458,
    PI: Math.PI
};

// Global calibration of the vacuum-field conductor-loss integrand, applied to
// the AC sums (both Re and Im parts, the internal-inductance/roughness
// reactance uses the same |H_t|^2 integral). Absorbs the uniform,
// frequency-flat underestimate of the discrete surface integral. Minimax-fit vs
// the tri-backend MQS. Recalibrate against the tri-backend MQS loss if the
// integrand changes.
export const VACUUM_LOSS_CAL = 1.093;

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

// Solve one matrix against several right-hand sides with a single
// factorization. The factorization dominates the solve cost, so k systems that
// share an operator (e.g. the odd and even modes of a differential pair, which
// differ only in their Dirichlet values) cost barely more than one.
async function solveWithWASMMulti(csr, Bs, useLU = false) {
    if (!WASMModuleInstance) {
        // Initialize the module if it hasn't been already
        WASMModuleInstance = await createWASMModule();
    }

    const nRhs = Bs.length;
    const N = Bs[0].length;
    const nnz = csr.values.length;

    const bytesNeeded = 10 * (12 * nnz + 20 * N) + 16 * N * (nRhs - 1);

    if (bytesNeeded > 1e9) {
      throw new Error(`Problem too large. Tried to allocate ${bytesNeeded/1e9} GB.`);
    }

    // Allocate memory
    const pRow = WASMModuleInstance._malloc(4 * (N + 1));
    const pCol = WASMModuleInstance._malloc(4 * nnz);
    const pVal = WASMModuleInstance._malloc(8 * nnz);
    const pB   = WASMModuleInstance._malloc(8 * N * nRhs);
    const pX   = WASMModuleInstance._malloc(8 * N * nRhs);

    try {
        // Re-acquire HEAP views to ensure they are current in case memory grew
        const currentHEAP32 = WASMModuleInstance.HEAP32;
        const currentHEAPF64 = WASMModuleInstance.HEAPF64;

        // Copy data - create views AFTER malloc
        const rowView = new Int32Array(currentHEAP32.buffer, pRow, N + 1);
        const colView = new Int32Array(currentHEAP32.buffer, pCol, nnz);
        const valView = new Float64Array(currentHEAPF64.buffer, pVal, nnz);
        const bView = new Float64Array(currentHEAPF64.buffer, pB, N * nRhs);

        rowView.set(csr.rowPtr);
        colView.set(csr.colIdx);
        valView.set(csr.values);
        for (let r = 0; r < nRhs; r++) bView.set(Bs[r], r * N);

        if (!WASMModuleInstance._solve_sparse_multi) {
            throw new Error("WASM function solve_sparse_multi not found. Module not loaded properly.");
        }

        // Call solver
        const status = WASMModuleInstance._solve_sparse_multi(
            N, nnz,
            pRow, pCol, pVal,
            nRhs, pB, pX,
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

        // Copy results
        const xView = new Float64Array(WASMModuleInstance.HEAPF64.buffer, pX, N * nRhs);
        const xs = [];
        for (let r = 0; r < nRhs; r++) {
            const x = new Float64Array(N);
            x.set(xView.subarray(r * N, (r + 1) * N));
            xs.push(x);
        }

        return xs;
    } finally {
        // Always free memory
        WASMModuleInstance._free(pRow);
        WASMModuleInstance._free(pCol);
        WASMModuleInstance._free(pVal);
        WASMModuleInstance._free(pB);
        WASMModuleInstance._free(pX);
    }
}

async function solveWithWASM(csr, B, useLU = false) {
    return (await solveWithWASMMulti(csr, [B], useLU))[0];
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

    if (!isArrayLike2D(V, ny, nx)) {
        errors.push("V must be a (ny, nx) 2D array matching mesh dimensions")
    }

    if (!isArrayLike2D(conductor_mask, ny, nx)) {
        errors.push("conductor_mask must be a (ny, nx) 2D array matching mesh dimensions")
    }

    if (!vacuum && !isArrayLike2D(epsilon_r, ny, nx)) {
        errors.push("epsilon_r must be a (ny, nx) 2D array matching mesh dimensions when vacuum=false")
    }

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

    // Check for at least one conductor
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

    // V
    for (let i = 0; i < ny; i++) {
        for (let j = 0; j < nx; j++) {
            const v = V[i][j];
            if (!Number.isFinite(v)) {
                errors.push(`V contains non-finite value at (${i}, ${j}): ${v}`);
                break;
            }
        }
    }

    // epsilon_r
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

export class FieldSolver2D {
    constructor() {
        this.x = null;
        this.y = null;
        this.V = null;  // Stored as array: [V] for single-ended, [V_odd, V_even] for differential
        this.epsilon_r = null;
        this.tand = null;
        this.conductor_mask = null; // 1 if conductor, 0 if dielectric
        this.solution_valid = false;
        this.use_causal_materials = true;

        // Computed fields - stored as array: [fields] for single-ended, [odd, even] for differential
        this.Ex = null;
        this.Ey = null;
    }

    /**
     * Pre-mesh sanity check: can this geometry be meshed at fine enough detail to be
     * resolved while keeping the mesh under the node budget?
     *
     * Two distinct requirements set the mesh size, and the binding one varies:
     *   1. FEATURE resolution — the smallest geometric feature (thin trace, coating,
     *      conductor edge) needs a few cells across it, setting the FINE cell size hFine.
     *   2. WAVELENGTH resolution — at high frequency the field-concentration ("active")
     *      region near the conductors must be sampled at a fraction of the in-medium
     *      wavelength. This is a GLOBAL, ungradeable floor over that region: a graded mesh
     *      can only be so coarse there. A genuinely electrically-large cross-section (e.g.
     *      a 1 m structure at 100 GHz → λ≈3 mm → ~hundreds of cells per side) needs an
     *      impossibly dense mesh. We apply this to the active region, NOT the auto-padded
     *      air domain, and use only ~3 cells/λ (just above Nyquist), so the gate fires only
     *      for the surely-unsolvable — never a normal mm-scale line at any frequency.
     *
     * The two backends pay for this differently, so the estimate is per-backend:
     *   - Rectilinear FDM is a TENSOR grid: nodes = nx·ny, and a fine feature forces fine
     *     graded LINES spanning the whole opposite axis. Budget = max_nodes mesh lines².
     *   - Triangular FEM grades in 2D: a fine feature costs only a LOCAL patch (so it can't
     *     be unsolvable on its own), and the full-wave eigensolve is ~4× heavier per entity,
     *     so its triangle budget is max_nodes/4 (memory parity — see tri_backend.js).
     *
     * Throws an Error (with the dominant cause and remedies) when the geometry cannot fit
     * the budget; returns silently when it can or when the domain isn't known yet.
     *
     * @param {number} maxNodes - The node budget (UI "Max Nodes"), default 20000.
     * @param {number} [freqOverride] - Frequency to evaluate the wavelength term at
     *     (solveModes solves at its own frequency, not this.freq).
     */
    // Warnings for 'open' boundaries that sit too close to the conductors. The open
    // boundary conditions (the FDM open stencil, the triangular backend's first-order
    // radiating ABC) approximate an unbounded exterior, which only holds where the
    // fringing field has mostly decayed — and its decay scale is the substrate
    // thickness. Require OPEN_CLEARANCE (3) substrate heights between each open wall
    // and the nearest conductor (this covers all line types: for a microstrip it is
    // ≥3·h from the trace to a side wall, and ≥3·h of air above the substrate
    // interface). Conductors touching a wall (coplanar ground pours, full-width
    // ground planes) are part of the boundary structure and are ignored.
    // Returns an array of human-readable warning strings (empty when fine).
    openBoundaryWarnings() {
        const out = [];
        const b = this.boundaries;
        if (!b || !this.conductors || !this.conductors.length) return out;
        const xMin = -this.domain_width / 2, xMax = this.domain_width / 2;
        const yMin = -(this.t_gnd ?? 0), yMax = this.domain_height;
        if (!(xMax > xMin) || !(yMax > yMin)) return out;
        // Decay scale: substrate stack thickness (non-air dielectrics); if there is
        // none (all-air line), fall back to the conductor stack height.
        let lo = Infinity, hi = -Infinity;
        for (const d of (this.dielectrics || [])) {
            if ((d.epsilon_r || 1) <= 1.001) continue;
            lo = Math.min(lo, d.y_min); hi = Math.max(hi, d.y_max);
        }
        if (!(hi > lo)) {
            for (const c of this.conductors) { lo = Math.min(lo, c.y_min); hi = Math.max(hi, c.y_max); }
        }
        const h = hi - lo;
        if (!(h > 0)) return out;
        const OPEN_CLEARANCE = 3;
        const tol = Math.max(xMax - xMin, yMax - yMin) * 1e-9;
        const walls = [
            ['left', b[0], c => c.x_min - xMin],
            ['right', b[1], c => xMax - c.x_max],
            ['top', b[2], c => yMax - c.y_max],
            ['bottom', b[3], c => c.y_min - yMin],
        ];
        const mm = (v) => `${(v * 1000).toPrecision(3)} mm`;
        const close = new Map();   // wall name → distance of nearest (non-touching) conductor
        for (const [name, bc, dist] of walls) {
            if (bc !== 'open') continue;
            let dMin = Infinity;
            for (const c of this.conductors) {
                const d = dist(c);
                if (d > tol && d < dMin) dMin = d;   // ≤ tol: touches this wall — intentional
            }
            if (dMin < OPEN_CLEARANCE * h) close.set(name, dMin);
        }
        if (!close.size) return out;
        // One combined warning for everything that tripped: left+right merge into
        // "sides" (a symmetric line trips both together), and the remaining wall
        // names are listed in one sentence ("sides and top", "left and top", …)
        // with the smallest clearance reported.
        if (close.has('left') && close.has('right')) {
            close.set('sides', Math.min(close.get('left'), close.get('right')));
            close.delete('left'); close.delete('right');
        }
        const parts = ['sides', 'left', 'right', 'top', 'bottom'].filter(n => close.has(n));
        const list = parts.length > 1
            ? parts.slice(0, -1).join(', ') + ' and ' + parts[parts.length - 1]
            : parts[0];
        const dMin = Math.min(...close.values());
        out.push(`The open boundary is too close to the conductors on the ${list}: only ${mm(dMin)} ` +
            `to the nearest conductor, less than ${OPEN_CLEARANCE}× the ${mm(h)} substrate height. ` +
            `The open-boundary approximation may be inaccurate this close to the fields; enlarge the ` +
            `enclosure or use grounded walls there.`);
        return out;
    }

    // modesOpts (truthy = guarding a Modes-tab solve, which ALWAYS runs the triangular
    // backend regardless of the Solver dropdown): {wavelengthDensity} — the Modes mesh
    // density, because that solve wavelength-caps the bulk of the WHOLE domain (see
    // TriBackend._wavelengthCap) rather than just the active patch.
    _check_meshability(maxNodes = 20000, freqOverride = undefined, modesOpts = null) {
        const freq = freqOverride ?? this.freq;
        const W = this.domain_width;
        const yBottom = -(this.t_gnd ?? 0);
        const H = (this.domain_height ?? NaN) - yBottom;
        // Subclasses that haven't set up a domain yet (or degenerate inputs) — nothing to check.
        if (!(W > 0) || !(H > 0) || !this.conductors || this.conductors.length === 0) return;

        // --- 1. Smallest geometric feature → the fine cell size hFine ---
        // Mirror what the real meshers key off: conductor widths/heights and dielectric
        // layer thicknesses (FDM corner_size = min_conductor_dimension/10; tri hFine =
        // min(thickness, width/4)). A few cells per feature.
        const featDims = [];
        for (const c of this.conductors) {
            if (Math.abs(c.width) > 0) featDims.push(Math.abs(c.width));
            if (Math.abs(c.height) > 0) featDims.push(Math.abs(c.height));
        }
        for (const d of (this.dielectrics || [])) {
            if (Math.abs(d.height) > 0) featDims.push(Math.abs(d.height));
        }
        const dMin = featDims.length ? Math.min(...featDims) : Math.min(W, H);
        const CELLS_PER_FEATURE = 3;
        const hFine = Math.max(dMin / CELLS_PER_FEATURE, 1e-12);

        // --- 2. The field-concentration ("active") region that must be wavelength-resolved ---
        // A bound quasi-TEM mode's field lives within ~a few substrate heights of the
        // conductors; the auto-sized domain pads far beyond that with near-field-free air
        // that only needs a coarse mesh. So the WAVELENGTH-resolution requirement applies
        // to the active region, NOT the full (often heavily over-padded) air domain — that
        // is what separates a normal mm-scale line at 100 GHz (active region ≪ λ, fine) from
        // a genuinely electrically-large 1 m cross-section at 100 GHz (active region ≫ λ).
        // Build the structure box from the conductors and the real substrate dielectrics
        // (ε_r > 1). The air fill is often modelled as a full-height ε_r=1 dielectric that
        // spans the whole padded domain — including it would wrongly make the structure look
        // domain-sized. The substrate thickness + conductor extent is the true field scale.
        const stackYs = [...this.conductors,
            ...(this.dielectrics || []).filter(d => (d.epsilon_r || 1) > 1.001)];
        let yLo = Infinity, yHi = -Infinity;
        for (const r of stackYs) { yLo = Math.min(yLo, r.y_min); yHi = Math.max(yHi, r.y_max); }
        if (!(yHi > yLo)) { yLo = 0; yHi = dMin; }       // degenerate fallback
        const G = Math.max(yHi - yLo, dMin);             // transverse structure scale (substrate stack)
        // Horizontal span of the conductor cluster, ignoring full-domain-width grounds
        // (which are absorbed into PEC walls and carry no localized field structure).
        let cxLo = Infinity, cxHi = -Infinity;
        for (const c of this.conductors) {
            if (Math.abs(c.width) >= W * 0.99) continue;
            cxLo = Math.min(cxLo, c.x_min); cxHi = Math.max(cxHi, c.x_max);
        }
        const clusterW = (cxHi > cxLo) ? (cxHi - cxLo) : dMin;
        // Field decays ~exponentially over a few G; pad the cluster by 3·G each side
        // vertically and horizontally, capped at the actual domain.
        const Lx = Math.min(W, clusterW + 6 * G);
        const Ly = Math.min(H, (yHi - yLo) + 6 * G);

        // --- 3. Cell sizes: coarse geometric background + wavelength-resolved active patch ---
        const epsMax = Math.max(1, ...(this.dielectrics || []).map(d => d.epsilon_r || 1));
        const geomCoarse = Math.max(Math.min(W, H) / 5, hFine);  // mesher bulk target
        let hWave = geomCoarse;                                   // active-region cell size
        let lambdaLimited = false;
        if (freq > 0) {
            const lambdaMin = CONSTANTS.C / (freq * Math.sqrt(epsMax));
            // ~3 cells per wavelength — just above the Nyquist floor. These quasi-TEM
            // backends solve bound modes on much coarser meshes than a full-wave λ/10
            // rule would demand, so this only flags structures so electrically large that
            // even a crude (alias-level) mesh blows the budget — i.e. surely unsolvable.
            const N_LAMBDA = 3;
            const hLambda = lambdaMin / N_LAMBDA;
            if (hLambda < hWave) { hWave = Math.max(hLambda, hFine); lambdaLimited = true; }
        }

        // --- 4. Estimate the mesh size for the active backend and compare to budget ---
        const fmtL = (m) => m >= 1e-3 ? `${(m * 1e3).toFixed(3)} mm` : `${(m * 1e6).toFixed(1)} µm`;
        const lambdaMm = freq > 0 ? fmtL(CONSTANTS.C / (freq * Math.sqrt(epsMax))) : '';
        const cause = lambdaLimited
            ? `the ${fmtL(Lx)}×${fmtL(Ly)} field region is electrically large at ${(freq / 1e9).toFixed(2)} GHz ` +
              `(λ≈${lambdaMm} in ε_r=${epsMax.toFixed(1)}), needing ~${fmtL(hWave)} cells to resolve the wavelength`
            : `resolving the ${fmtL(dMin)} feature across the ${fmtL(W)}×${fmtL(H)} domain`;
        const remedy = lambdaLimited
            ? 'reduce the structure/enclosure size, lower the maximum frequency, or raise Max Nodes'
            : 'enlarge the smallest feature, shrink the domain, or raise Max Nodes';

        if (this.mesh_backend === 'triangular' || modesOpts) {
            // Triangles grade in 2D, so a fine feature only costs a small local patch, never
            // enough on its own to be "surely unsolvable" (the suite meshes few µm features
            // in cm-scale domains fine). The only triangular blow-up that is genuinely
            // unsolvable is an electrically-large field region: a coarse background plus a
            // wavelength-resolved active patch that exceeds the triangle budget.
            const FW_NODES_PER_TRI = 4;
            const triBudget = Math.max(800, maxNodes / FW_NODES_PER_TRI);
            // Background floor is optimistic (cells ~half the thin domain dimension):
            // the mesher grades field-free air far coarser than the geomCoarse used for
            // the FDM line estimate, so min(W,H)/5 falsely rejected high-aspect (wide,
            // thin) domains the real mesher handles within budget (e.g. a 38 mm×0.27 mm
            // stripline meshes in ~4.1k tris vs the ~6.8k that formula claimed). Reject
            // only surely-unsolvable cases; the backend's coarsen-and-rebuild loop and
            // the eigenSolveBytes guard catch anything the estimate misses, cleanly.
            let triCoarse = Math.max(Math.min(W, H) / 2, hFine);
            let tris, triCause = cause, triRemedy = remedy;
            if (modesOpts && freq > 0) {
                // Modes solve: the bulk of the WHOLE domain is wavelength-capped at the
                // Modes mesh density (default 12 cells/λ) so cavity/higher-order modes
                // stay resolved — coarsening below that mis-classifies them, so an
                // electrically huge domain must be rejected here rather than degraded
                // through the coarsen-and-rebuild loop (whose FIRST gmsh build would
                // also be enormous).
                const nLambda = modesOpts.wavelengthDensity > 0 ? modesOpts.wavelengthDensity : 12;
                const hBulk = CONSTANTS.C / (freq * Math.sqrt(epsMax)) / nLambda;
                if (hBulk < triCoarse) {
                    triCoarse = Math.max(hBulk, hFine);
                    triCause = `the ${fmtL(W)}×${fmtL(H)} domain is electrically large at ` +
                        `${(freq / 1e9).toFixed(2)} GHz (λ≈${lambdaMm} in ε_r=${epsMax.toFixed(1)}); ` +
                        `the Modes solve resolves the whole domain at ${nLambda} cells/λ ` +
                        `(~${fmtL(triCoarse)} cells) to keep cavity/higher-order modes trustworthy`;
                    triRemedy = 'lower the Modes frequency or Mesh density, shrink the enclosure, or raise Max Nodes';
                }
                tris = 2 * (W / triCoarse) * (H / triCoarse);
            } else {
                tris = 2 * (W / triCoarse) * (H / triCoarse);     // coarse background
                if (lambdaLimited) tris += 2 * (Lx / hWave) * (Ly / hWave);   // active wavelength patch
            }
            if (tris > triBudget) {
                throw new Error(
                    `Geometry cannot be meshed for the full-wave (triangular) solver within the node budget: ` +
                    `${triCause}, needing ~${Math.round(tris).toLocaleString()} triangles vs a budget of ` +
                    `${Math.round(triBudget).toLocaleString()} (Max Nodes ${maxNodes.toLocaleString()}). To proceed, ${triRemedy}.`);
            }
        } else {
            // Rectilinear tensor grid: nodes = nx·ny. The coarse geometric cell sets the
            // baseline line count per axis; the wavelength-resolved active region adds fine
            // lines over its span; a fine feature adds graded transition lines (logarithmic).
            const GROWTH = 1.3;                                    // graded neighbour ratio
            const gradeLines = Math.max(0, Math.log(geomCoarse / hFine) / Math.log(GROWTH));
            const axisLines = (L, La) => L / geomCoarse + (lambdaLimited ? La / hWave : 0) + 2 * gradeLines;
            const nx = axisLines(W, Lx), ny = axisLines(H, Ly);
            const nodes = nx * ny;
            if (nodes > maxNodes) {
                throw new Error(
                    `Geometry cannot be meshed for the rectilinear (FDM) solver within the node budget: ` +
                    `${cause}, needing ~${Math.round(nx)}×${Math.round(ny)} ≈ ${Math.round(nodes).toLocaleString()} mesh nodes vs a budget of ` +
                    `${maxNodes.toLocaleString()} (Max Nodes). To proceed, ${remedy}.`);
            }
        }
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
     * Voltage array driving the positive trace at vp and the negative trace at vn
     * (grounds at 0). Used for per-conductor / modal-eigenvector excitation.
     */
    _create_voltage_array_drive(vp, vn) {
        const ny = this.y.length, nx = this.x.length;
        const V = Array(ny).fill().map(() => new Float64Array(nx));
        for (let i = 0; i < ny; i++) {
            for (let j = 0; j < nx; j++) {
                if (this.ground_mask[i][j]) V[i][j] = 0.0;
                else if (this.signal_p_mask[i][j]) V[i][j] = vp;
                else if (this.signal_n_mask[i][j]) V[i][j] = vn;
            }
        }
        return V;
    }

    // Charge on the positive / negative trace for a given potential field (uses the
    // validated Gauss-flux integrator with the per-trace mask). Returns |Q|.
    _trace_charge(V, mask, vacuum) {
        const orig = this.signal_mask;
        this.signal_mask = mask;
        const Q = this.calculate_capacitance(V, vacuum);
        this.signal_mask = orig;
        return Q;
    }

    /**
     * Solve Laplace equation for the given voltage array.
     * @param {Array<Float64Array>} V - 2D voltage array with conductor boundary conditions set
     * @param {boolean} vacuum - If true, solve with vacuum permittivity
     * @param {function} onProgress - Optional progress callback
     * @returns {Array<Float64Array>} - The solved voltage array (same reference as input)
     */
    async solve_laplace(V, vacuum = false, onProgress = null) {
        return (await this.solve_laplace_multi([V], vacuum))[0];
    }

    // Solve the same operator (grid + epsilon + conductor mask) for several sets
    // of Dirichlet drive voltages at once. The matrix does not depend on the
    // drive, only the right-hand side does, so all systems share a single
    // assembly and a single factorization (see solveWithWASMMulti). Each Vs[m]
    // is filled in place and the array of solutions is returned.
    async solve_laplace_multi(Vs, vacuum = false) {
        // Ensure mesh is generated
        if (this.ensure_mesh) {
            this.ensure_mesh();
        }

        for (const V of Vs) {
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
        }

        const ny = this.y.length, nx = this.x.length;
        const dx = diff(this.x), dy = diff(this.y);
        const N = nx * ny;
        const idx = (i, j) => i * nx + j;

        const get_er = (i, j) => vacuum ? 1.0 : this.epsilon_r[i][j];
        const is_cond = (i, j) => this.conductor_mask[i][j];

        // Remove mesh nodes internal to conductors
        // E-field inside conductors is 0.
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
        const Bs = Vs.map(() => new Float64Array(N_unknown));
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
                        for (let m = 0; m < Vs.length; m++)
                            Bs[m][n] -= c * Vs[m][ii][jj];
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
        const xs = await solveWithWASMMulti(csr, Bs, true);

        // Reconstruct solutions for full mesh
        for (let m = 0; m < Vs.length; m++) {
            const V = Vs[m], x = xs[m];
            for (let k = 0; k < N_unknown; k++) {
                const n = red_to_full[k];
                const i = (n / nx) | 0;
                const j = n % nx;
                V[i][j] = x[k];
            }
        }

        return Vs;
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
     * - R_dc: DC resistance from conductor cross-sectional area
     * - R_ac: AC resistance from skin effect and surface roughness
     *
     * Two integrand variants:
     *
     * vacuum_fields = true (production for rect-based solvers): Ex/Ey are the
     * vacuum (C0) solve fields and Z0 is Z0_vac = 1/(c*C0). In the quasi-TEM
     * skin-effect limit the H pattern is the harmonic conjugate of the vacuum
     * potential (the same identity as L_ext = 1/(c^2*C0)), so H_t = E_n(vac)/η0,
     * the surface current distribution, which is permittivity independent. The
     * AC sums are scaled by VACUUM_LOSS_CAL. Algebraically identical to the
     * legacy variant for homogeneous fill apart from the calibration.
     *
     * vacuum_fields = false (legacy): Ex/Ey are the dielectric solve fields,
     * H_t = E_n*√εr(local)/η0, Z0 is the line impedance. For mixed dielectric
     * this uses the charge distribution as a current proxy, it overestimates
     * corner-dominated microstrip loss and its substrate-interface corner
     * singularity makes the sum mesh-divergent.
     *
     * @param {Array<Array<number>>} Ex - Electric field x-component
     * @param {Array<Array<number>>} Ey - Electric field y-component
     * @param {number} Z0 - Line impedance (legacy) or vacuum impedance 1/(c*C0)
     * @param {boolean} vacuum_fields - Ex/Ey are vacuum-solve fields
     * @returns {{R_ac: number, R_dc: number, R_total: number, L_internal: number}}
     */
    calculate_conductor_loss(Ex, Ey, Z0, vacuum_fields = false) {
        if (!this.solution_valid) throw new Error("Fields invalid");

        const { signal_area, ground_area } = this._calculate_conductor_area();

        // DC resistance per unit length for transmission line
        // Current flows through signal conductor and returns through ground (series connection)
        const R_signal = 1.0 / (this.sigma_cond * signal_area);
        const R_ground = 1.0 / (this.sigma_cond * ground_area);
        const R_dc = R_signal + R_ground;

        // Handle DC case (frequency = 0)
        if (this.freq === 0) {
            return {
                R_ac: 0,
                R_dc: R_dc,
                R_total: R_dc,
                L_internal: 0
            };
        }

        // Use roughness from constructor
        const rq = this.rq || 0;

        // Default surface impedance (no plating)
        const Z_surf_default = calculate_Zrough(this.freq, this.sigma_cond, rq);

        // Cache for per-conductor surface impedances to avoid redundant layered solves
        // Key: "ci_surface" e.g. "3_top", value: Complex Z_surf
        const Z_cache = new Map();

        // Helper: check if a point is at a conductor corner
        // Returns corner type: 'bottom-left', 'bottom-right', or null
        const getCornerType = (i, j, ci, direction) => {
            if (!this.conductors || ci < 0 || !this.conductor_id) return null;
            const cond = this.conductors[ci];
            if (!cond) return null;

            // Only check for bottom corners (direction 'u' = bottom face)
            if (direction !== 'u') return null;

            // Check if there's also a horizontal neighbor with the same conductor
            const has_left = (j > 0) && this.conductor_id[i] && this.conductor_id[i][j-1] === ci;
            const has_right = (j < this.x.length - 1) && this.conductor_id[i] && this.conductor_id[i][j+1] === ci;

            if (has_left) return 'bottom-left';
            if (has_right) return 'bottom-right';
            return null;
        };

        // Helper: get surface impedance for a conductor boundary segment
        // Now with corner detection and geometric coverage from thick side plating
        const getZsurf = (ci, direction, i, j, dl) => {
            if (!this.conductors || ci < 0) return Z_surf_default;
            const cond = this.conductors[ci];
            if (!cond || !cond.plating) return Z_surf_default;

            // Map direction to surface face:
            // 'u' = dielectric is below conductor neighbor = conductor's BOTTOM face
            // 'd' = dielectric is above conductor neighbor = conductor's TOP face
            // 'l','r' = SIDE faces
            let surface;
            if (direction === 'd') surface = 'top';
            else if (direction === 'u') surface = 'bottom';
            else surface = 'sides';

            // Geometric coverage from TOP plating extending down the sides
            // If only top plating (not sides), plating extends down from top edge
            // Uses fractional coverage for smooth parameter sweeps
            if ((direction === 'l' || direction === 'r') && cond.plating.top && !cond.plating.sides) {
                const t = cond.plating.thickness;
                // Cell spans [y[i], y[i] + dl] in y-direction
                const y_start = this.y[i];
                const y_end = y_start + dl;
                // Top plating covers [y_max - t, y_max]
                const overlap = Math.max(0, Math.min(cond.y_max, y_end) - Math.max(cond.y_max - t, y_start));
                const fraction = dl > 0 ? Math.min(overlap / dl, 1.0) : 0;

                if (fraction > 0) {
                    const key_top_side = `${ci}_top_side_plating`;
                    let Z_plating;
                    if (Z_cache.has(key_top_side)) {
                        Z_plating = Z_cache.get(key_top_side);
                    } else {
                        Z_plating = calculate_Zrough(
                            this.freq, cond.plating.sigma, cond.plating.rq
                        );
                        Z_cache.set(key_top_side, Z_plating);
                    }

                    if (fraction >= 1.0) return Z_plating;

                    // Weighted average with bulk side impedance for uncovered part
                    return new Complex(
                        fraction * Z_plating.re + (1 - fraction) * Z_surf_default.re,
                        fraction * Z_plating.im + (1 - fraction) * Z_surf_default.im
                    );
                }
            }

            // Geometric coverage of bottom surface by thick side plating
            // Uses fractional coverage for smooth parameter sweeps
            if (direction === 'u' && cond.plating.sides && !cond.plating.bottom && cond.plating.thick_corners) {
                const t = cond.plating.thickness;
                // Cell spans [x[j], x[j] + dl] in x-direction
                const x_start = this.x[j];
                const x_end = x_start + dl;
                // Left side plating covers [x_min, x_min + t]
                const left_overlap = Math.max(0, Math.min(cond.x_min + t, x_end) - Math.max(cond.x_min, x_start));
                // Right side plating covers [x_max - t, x_max]
                const right_overlap = Math.max(0, Math.min(cond.x_max, x_end) - Math.max(cond.x_max - t, x_start));
                const fraction = dl > 0 ? Math.min((left_overlap + right_overlap) / dl, 1.0) : 0;

                if (fraction > 0) {
                    const key_corner = `${ci}_corner_plating`;
                    let Z_plating;
                    if (Z_cache.has(key_corner)) {
                        Z_plating = Z_cache.get(key_corner);
                    } else {
                        // Side plating material with bulk surface roughness
                        Z_plating = calculate_Zrough(
                            this.freq, cond.plating.sigma, rq
                        );
                        Z_cache.set(key_corner, Z_plating);
                    }

                    if (fraction >= 1.0) return Z_plating;

                    // Weighted average: covered part uses plating, rest uses bulk
                    return new Complex(
                        fraction * Z_plating.re + (1 - fraction) * Z_surf_default.re,
                        fraction * Z_plating.im + (1 - fraction) * Z_surf_default.im
                    );
                }
            }

            // Check for bottom corner
            const cornerType = getCornerType(i, j, ci, direction);

            // At bottom corners with side plating enabled, use single-layer plating impedance
            // This models plating extending from sides to bottom at corners
            if (cornerType && cond.plating.sides && cond.plating.thick_corners) {
                // Determine corner size (characteristic dimension)
                const corner_size = Math.min(cond.width, Math.abs(cond.height)) / 10;

                // Get corner plating impedance (single-layer, no bulk)
                const key_corner = `${ci}_corner_plating`;
                let Z_corner;
                if (Z_cache.has(key_corner)) {
                    Z_corner = Z_cache.get(key_corner);
                } else {
                    // At corners: single-layer with plating sigma and bulk rq
                    // - sigma: plating material (extends from sides)
                    // - rq: bulk surface roughness (bottom surface preparation)
                    Z_corner = calculate_Zrough(
                        this.freq, cond.plating.sigma, rq  // Use bulk rq, not plating.rq
                    );
                    Z_cache.set(key_corner, Z_corner);
                }

                // If mesh cell is small (pure corner region), use pure corner plating impedance
                if (dl < corner_size) {
                    return Z_corner;
                }

                // If mesh cell is large and includes both corner and bulk,
                // average based on corner_size fraction
                const corner_fraction = corner_size / dl;

                // Get bottom surface impedance
                let Z_bottom;
                if (cond.plating.bottom) {
                    const key_bottom = `${ci}_bottom`;
                    if (Z_cache.has(key_bottom)) {
                        Z_bottom = Z_cache.get(key_bottom);
                    } else {
                        Z_bottom = calculate_Zrough_layered(
                            this.freq, this.sigma_cond,
                            cond.plating.rq, cond.plating.sigma, cond.plating.thickness
                        );
                        Z_cache.set(key_bottom, Z_bottom);
                    }
                } else {
                    Z_bottom = Z_surf_default;
                }

                // Weighted average: corner region uses corner plating impedance, bulk uses bottom impedance
                const Z_avg_re = corner_fraction * Z_corner.re + (1 - corner_fraction) * Z_bottom.re;
                const Z_avg_im = corner_fraction * Z_corner.im + (1 - corner_fraction) * Z_bottom.im;
                return new Complex(Z_avg_re, Z_avg_im);
            }

            // Standard surface impedance (no corner effects)
            if (!cond.plating[surface]) return Z_surf_default;

            const key = `${ci}_${surface}`;
            if (Z_cache.has(key)) return Z_cache.get(key);

            const Z = calculate_Zrough_layered(
                this.freq, this.sigma_cond,
                cond.plating.rq, cond.plating.sigma, cond.plating.thickness
            );
            Z_cache.set(key, Z);
            return Z;
        };

        const ny = this.y.length;
        const nx = this.x.length;
        const dx_array = diff(this.x);
        const dy_array = diff(this.y);

        const get_dx = j => (j >= 0 && j < dx_array.length) ? dx_array[j] : dx_array[dx_array.length - 1];
        const get_dy = i => (i >= 0 && i < dy_array.length) ? dy_array[i] : dy_array[dy_array.length - 1];

        let sum_H2_dl_R = 0.0; // Sum for Resistance
        let sum_H2_dl_L = 0.0; // Sum for Inductance

        const isSignal = (i, j) => this.signal_mask[i][j];
        const isGround = (i, j) => this.ground_mask[i][j];
        const isConductor = (i, j) => isSignal(i,j) || isGround(i,j);

        for (let i = 1; i < ny - 1; i++) {
            for (let j = 1; j < nx - 1; j++) {
                if (isConductor(i, j)) continue;

                const neighbors = [
                    { ni: i, nj: j + 1, direction: 'r', dl_func: get_dy, idx: i },
                    { ni: i, nj: j - 1, direction: 'l', dl_func: get_dy, idx: i },
                    { ni: i + 1, nj: j, direction: 'u', dl_func: get_dx, idx: j },
                    { ni: i - 1, nj: j, direction: 'd', dl_func: get_dx, idx: j },
                ];

                for (const { ni, nj, direction, dl_func, idx: dl_idx } of neighbors) {
                    if (ni < 0 || ni >= ny || nj < 0 || nj >= nx) continue;

                    if (isConductor(ni, nj)) {
                        const Ex_val = Ex[i][j];
                        const Ey_val = Ey[i][j];

                        let E_norm = 0.0;
                        if (direction === 'r' || direction === 'l') E_norm = Math.abs(Ex_val);
                        else E_norm = Math.abs(Ey_val);

                        const Z0_freespace = 376.73;
                        // Vacuum fields: H pattern is the vacuum dual, no eps factor.
                        // Legacy: local plane-wave relation on the dielectric field.
                        const eps_fac = vacuum_fields ? 1.0 : Math.sqrt(this.epsilon_r[i][j]);
                        const H_tan = E_norm * eps_fac / Z0_freespace;

                        const dl = dl_func(dl_idx);

                        // Look up per-surface impedance (with plating if applicable)
                        const ci = this.conductor_id ? this.conductor_id[ni][nj] : -1;
                        const Z_surf = getZsurf(ci, direction, i, j, dl);

                        const H2_dl = H_tan * H_tan * dl;
                        sum_H2_dl_R += Z_surf.re * H2_dl;
                        sum_H2_dl_L += Z_surf.im * H2_dl;
                    }
                }
            }
        }

        // Power normalization factor: differential has 0.5 factor
        // This is because we integrate over both traces but report normalized loss
        const power_factor = this.is_differential ? 0.5 : 1.0;

        // Vacuum variant: |H|^2 per unit current is |H_vac per 1V|^2*Z0_vac^2, since the
        // vacuum drive at 1V carries I_vac = 1/Z0_vac (legacy: same algebra with the
        // line Z0). The calibration applies to both sums (same integral and same bias).
        const cal = vacuum_fields ? VACUUM_LOSS_CAL : 1.0;
        const Z0_sq = Z0 * Z0 * cal;

        // AC Resistance per unit length from skin effect (Ohm/m)
        const R_ac = power_factor * sum_H2_dl_R * Z0_sq;

        // This doesn't hold if conductor thickness is smaller than skin depth
        // Need to solve magnetic field for accurate L_internal at low frequency
        // but is not a problem at even moderately high frequency >1 MHz.
        // In practice very minimal error since DC can be solved correctly.
        const L_internal = power_factor * sum_H2_dl_L * Z0_sq / (2 * Math.PI * this.freq);

        const R_total = Math.sqrt(R_dc * R_dc + R_ac * R_ac);

        return { R_ac, R_dc, R_total, L_internal };
    }

    // Conductor loss for a solved mode, choosing the integrand variant:
    // rect-based solvers (conductor_id present) use the vacuum-field integrand
    // when the mode's vacuum fields are available.
    // The only production caller lacking vacuum fields on a rect solver is
    // _solve_single_mode(vacuum_first=false), whose loss output is discarded
    // and recomputed by the caller with the cached vacuum fields.
    _mode_conductor_loss(Ex, Ey, Z0, C0, Ex0, Ey0) {
        if (this.conductor_id && Ex0 && Ey0 && C0 > 0) {
            const Z0_vac = 1 / (CONSTANTS.C * C0);
            return this.calculate_conductor_loss(Ex0, Ey0, Z0_vac, true);
        }
        return this.calculate_conductor_loss(Ex, Ey, Z0);
    }

    calculate_dielectric_loss(Ex, Ey, Z0) {
        if (!this.solution_valid) {
            throw new Error("Fields (Ex, Ey) are not valid. Run compute_fields() first.");
        }

        // No dielectric loss at DC
        // If material conductivity is implemented this is not true
        if (this.freq === 0) {
            return 0;
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

                Pd += 0.5 * (2 * Math.PI * this.freq) * CONSTANTS.EPS0 * this.epsilon_r[i][j] * this.tand[i][j] * E2 * dA;
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

        // Handle DC case (frequency = 0)
        if (this.freq === 0) {
            // At DC, Zc = sqrt(R/G) = sqrt(R/0) = infinity
            // For S-parameter calculations, use a very large impedance
            const Zc = new Complex(1e12, 0);  // Effectively infinite impedance

            // eps_eff at DC is calculated from C/C0
            // From Z0 = 1/(c*sqrt(C*C0)), we get C0 = 1/(c^2*Z0^2*C)
            // Therefore eps_eff = C/C0 = c^2 * Z0^2 * C^2
            const c2 = CONSTANTS.C * CONSTANTS.C;
            const eps_eff_dc = c2 * Z0_mode * Z0_mode * C_mode * C_mode;

            return {
                Zc: Zc,
                rlgc: {
                    R: R_total,
                    L: L_total,
                    G: G,
                    C: C_mode
                },
                eps_eff_mode: eps_eff_dc,
                L_internal: L_internal,
                L_external: L_ext
            };
        }

        // Re-calculate complex Zc and Epsilon_eff with the new L and R
        const omega = 2 * Math.PI * this.freq;

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

        // 6. Calculate RLGC and complex Z0
        const { Zc, rlgc, eps_eff_mode } = this.rlgc(R_total, L_internal, alpha_diel_db_m, C_with_diel, Z0);

        // Calculate conductor loss alpha from R_total for reporting
        const alpha_cond_db_m = 8.686 * R_total / (2 * Z0);
        const total_alpha_db_m = alpha_cond_db_m + alpha_diel_db_m;

        return {
            Z0: Z0, // Characteristic Impedance (static approximation)
            Zc: Zc, // Complex Characteristic Impedance (includes loss)
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

    refine_mesh_multi(modes, frac = 0.15) {
        /**
         * Mesh refinement using combined metrics from multiple modes.
         * Each mode's metrics are normalized to equal total weight before summing
         * so that a mode with weaker absolute fields (e.g. even mode) still gets
         * equal refinement budget relative to the dominant odd mode.
         */
        const nx_intervals = this.x.length - 1;
        const ny_intervals = this.y.length - 1;
        const x_combined = new Float64Array(nx_intervals);
        const y_combined = new Float64Array(ny_intervals);

        for (const { V, Ex, Ey } of modes) {
            const { x_metrics, y_metrics } = this._compute_refine_metrics(V, Ex, Ey);
            const total = x_metrics.reduce((s, v) => s + v, 0) +
                          y_metrics.reduce((s, v) => s + v, 0);
            const scale = total > 0 ? 1 / total : 1;
            for (let j = 0; j < nx_intervals; j++) x_combined[j] += x_metrics[j] * scale;
            for (let i = 0; i < ny_intervals; i++) y_combined[i] += y_metrics[i] * scale;
        }

        const { selected_x, selected_y } = this._select_lines_to_refine(x_combined, y_combined, frac);
        this._refine_selected_lines(selected_x, selected_y);

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

    // General two-conductor modal analysis for a differential pair, run on the converged
    // mesh as a post-pass. Drives each trace independently to assemble the full 2×2
    // capacitance matrices [C] (dielectric) and [C0] (vacuum), diagonalises [C0]⁻¹[C] to
    // get the two genuine line modes (eps_eff = eigenvalues, modal voltage ratios =
    // eigenvectors), and builds each mode's fields/loss from the eigenvector combination of
    // the per-trace fields. Replaces the odd/even drive results, which only coincide with
    // the modes for a SYMMETRIC pair (where the eigenvectors come out as [1,∓1], so this
    // reproduces the existing odd/even basis exactly — no change to symmetric results).
    async _solve_modal_differential() {
        // Per-trace static fields (dielectric + vacuum). The Laplace solve is linear in the
        // drive, so the field for any drive [vp,vn] is vp·(A field) + vn·(B field).
        const solveDrive = async (vp, vn, vac) => {
            const V = await this.solve_laplace(this._create_voltage_array_drive(vp, vn), vac);
            return { V, ...this.compute_fields(V) };
        };
        const A = await solveDrive(1, 0, false), Av = await solveDrive(1, 0, true);
        const B = await solveDrive(0, 1, false), Bv = await solveDrive(0, 1, true);
        // Full 2×2 Maxwell capacitance matrices from per-trace charges (Gauss flux). Diagonal
        // is the self charge (driven trace), off-diagonal the induced charge on the other
        // trace (negative). Symmetrised for reciprocity. Building the matrix this way (rather
        // than the energy integral) keeps the validated absolute scale, so C_k = ½·vᵀ·Cm·v
        // reproduces the existing C_odd/C_even exactly for a symmetric pair.
        const sp = this.signal_p_mask, sn = this.signal_n_mask;
        const maxwell = (DA, DB, vac) => {
            const m12 = -0.5 * (this._trace_charge(DA.V, sn, vac) + this._trace_charge(DB.V, sp, vac));
            return [[this._trace_charge(DA.V, sp, vac), m12],
                    [m12, this._trace_charge(DB.V, sn, vac)]];
        };
        const Cm = maxwell(A, B, false), Cm0 = maxwell(Av, Bv, true);
        const quad = (M, v) => 0.5 * (v[0] * v[0] * M[0][0] + 2 * v[0] * v[1] * M[0][1] + v[1] * v[1] * M[1][1]);
        // Shared symmetric/degenerate/modal decision (thresholds, ordering — see
        // classifyModalDecomposition; the triangular backend uses the identical guard).
        // Symmetric: null physMatrix, keep the odd/even drive results (returning null
        // signals the caller). Degenerate: physMatrix without Tv still drives the
        // asymmetric 4-port S-matrix while the odd/even per-mode results are kept.
        const { physMatrix, modalVecs } = classifyModalDecomposition(Cm, Cm0, this.conductors, this.dielectrics);
        this._modalPhys = physMatrix;
        if (!modalVecs) return null;
        const comb = (a, P, b, Q) => P.map((row, i) => row.map((val, j) => a * val + b * Q[i][j]));
        const results = [];
        ['odd', 'even'].forEach((label, li) => {
            const v = modalVecs[li];
            const V = comb(v[0], A.V, v[1], B.V);
            const Ex = comb(v[0], A.Ex, v[1], B.Ex), Ey = comb(v[0], A.Ey, v[1], B.Ey);
            // Same eigenvector combination of the vacuum drive fields: the mode's
            // vacuum field, feeding the conductor-loss integrand (_mode_conductor_loss).
            const Ex0 = comb(v[0], Av.Ex, v[1], Bv.Ex), Ey0 = comb(v[0], Av.Ey, v[1], Bv.Ey);
            const Ck = quad(Cm, v);          // ½·vᵀ·Cm·v  (mode capacitance, matches C_odd/C_even)
            const C0k = quad(Cm0, v);
            const eps_eff = Ck / C0k;
            const Z0 = 1 / (CONSTANTS.C * Math.sqrt(Ck * C0k));
            const { R_total, L_internal } = this._mode_conductor_loss(Ex, Ey, Z0, C0k, Ex0, Ey0);
            const alpha_d = this.calculate_dielectric_loss(Ex, Ey, Z0);
            const { Zc, rlgc, eps_eff_mode, L_external } = this.rlgc(R_total, L_internal, alpha_d, Ck, Z0);
            const alpha_c = 8.686 * R_total / (2 * Zc.re);
            results.push({
                mode: label, Z0, eps_eff: eps_eff_mode, C: Ck, C0: C0k, RLGC: rlgc, Zc,
                alpha_c, alpha_d, alpha_total: alpha_c + alpha_d, L_internal, L_external,
                V, Ex, Ey, Ex0, Ey0, modalVec: v,
            });
        });
        return results;
    }

    // Signal-trace capacitance for a solved potential. For a differential pair the
    // charge is averaged over the two traces: an asymmetric pair (e.g. broadside
    // stripline with unequal top/bottom dielectric heights) carries a different charge
    // on each trace, while for a symmetric pair the average changes nothing.
    _signal_capacitance(V, vacuum) {
        if (!this.is_differential) {
            return this.calculate_capacitance(V, vacuum);
        }
        const orig_signal_mask = this.signal_mask;
        this.signal_mask = this.signal_p_mask;
        const Cp = this.calculate_capacitance(V, vacuum);
        this.signal_mask = this.signal_n_mask;
        const Cn = this.calculate_capacitance(V, vacuum);
        this.signal_mask = orig_signal_mask;
        return 0.5 * (Cp + Cn);
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
        let Ex0, Ey0;

        if (vacuum_first) {
            // Calculate C0 (vacuum capacitance)
            V = this._create_voltage_array(mode);
            V = await this.solve_laplace(V, true);
            C0 = this._signal_capacitance(V, true);
            // Vacuum fields: the conductor-loss integrand for rect-based solvers
            // (the quasi-TEM H pattern, see _mode_conductor_loss). Frequency- and
            // ε-independent, so cached mode results reuse them across the sweep.
            if (this.conductor_id) ({ Ex: Ex0, Ey: Ey0 } = this.compute_fields(V));
        }

        // Solve with dielectric
        V = this._create_voltage_array(mode);
        V = await this.solve_laplace(V, false);
        const C = this._signal_capacitance(V, false);

        // Calculate fields
        const { Ex, Ey } = this.compute_fields(V);

        // Calculate impedance
        let eps_eff, Z0;
        if (C0 !== undefined) {
            eps_eff = C / C0;
            Z0 = 1 / (CONSTANTS.C * Math.sqrt(C * C0));
        }

        // Calculate conductor losses with surface roughness and DC resistance
        const { R_ac, R_dc, R_total, L_internal } = this._mode_conductor_loss(Ex, Ey, Z0, C0, Ex0, Ey0);

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
            V, Ex, Ey,
            Ex0, Ey0
        };
    }

    // Lazily create the triangular FEM backend (and build its mesh + cached
    // static solve). The import is dynamic so the gmsh/eigen WASM is only loaded
    // when the user selects the triangular backend. Persists across the sweep.
    async _ensureTriBackend(onProgress = null, extraOpts = null, shouldStop = null) {
        // Merge the solver-mode opts (lossMethod) with the UI adaptive controls.
        // A call WITH extraOpts (solve_adaptive) pins the effective opts: a cached
        // backend built with different opts (e.g. the user changed Max Nodes) is
        // discarded and rebuilt. A call WITHOUT extraOpts (mid-sweep
        // computeAtFrequency) reuses whatever the initial solve built.
        const wantOpts = extraOpts ? { ...(this.tri_opts || {}), ...extraOpts } : null;
        const wantKey = wantOpts ? JSON.stringify(wantOpts) : null;
        if (this._triBackend &&
            (wantKey === null || wantKey === this._triBackendOptsKey)) {
            return this._triBackend;
        }
        const { initTriBackend, TriBackend } = await import('./tri_solver/tri_backend.js');
        const ctx = await initTriBackend();
        const opts = wantOpts ?? { ...(this.tri_opts || {}) };
        const tri = new TriBackend(ctx, this, opts);
        // Cache only AFTER the mesh builds: a throw here must not leave a broken
        // (mesh-less) backend behind for the next call to reuse.
        await tri.buildMesh(onProgress, shouldStop);   // emits real adaptive-refinement passes (live)
        this._triBackend = tri;
        this._triBackendOptsKey = wantKey ?? JSON.stringify(opts);
        return this._triBackend;
    }

    // ---- Mode viewer ----------------------------------------------------------
    // Solve the full-wave eigenproblem for the lowest `nev` modes at `freq` and return
    // a classified list (see TriBackend.solveModes). Always uses the triangular full-wave
    // backend on the FULL domain (symmetry off) so symmetric AND antisymmetric higher-order
    // modes both appear. A dedicated backend instance is cached on `_modesBackend` so it
    // does not disturb the main (possibly half-domain) solve in `_triBackend`.
    // refineOpts wires the sidebar adaptive-mesh controls (maxRefineIters ← Max
    // Iterations, refineTol ← Tolerance, maxNodes ← Max Nodes) into buildMesh's
    // refinement loop, exactly like solve_adaptive does for the main solve, plus
    // wavelengthDensity ← the Modes tab's Mesh density (cells/λ for the bulk
    // wavelength cap — see TriBackend._wavelengthCap).
    async solveModes(freq, nev = 4, onProgress = null, refineOpts = {}) {
        // Same pre-mesh guard as solve_adaptive, evaluated at the MODES frequency and
        // ALWAYS on the triangular estimate — this solve runs the triangular backend
        // regardless of the sidebar Solver selection, so branching on this.mesh_backend
        // here would apply the FDM tensor-grid rules (falsely rejecting fine features
        // the tri mesher grades locally, and skipping the whole-domain wavelength check
        // that stops an electrically huge modes mesh from hanging gmsh).
        this._check_meshability(refineOpts.maxNodes ?? 20000, freq,
            { wavelengthDensity: refineOpts.wavelengthDensity });
        const { initTriBackend, TriBackend } = await import('./tri_solver/tri_backend.js');
        const ctx = await initTriBackend();
        // modesFreq lets buildMesh size the bulk to the wavelength at this frequency, so
        // high-frequency cavity/higher-order modes are resolved (and not mis-flagged as
        // spurious by the mesh-convergence test).
        // Mesh.Algorithm 1 (MeshAdapt) gives the cleanest eigenmode spectrum for the mode
        // viewer (fewer spurious low-ε_eff artifacts than the default frontal-Delaunay).
        const { shouldStop = null, ...meshOpts } = refineOpts;
        const opts = { ...(this.tri_opts || {}), symmetry: false, modesFreq: freq,
            gmshOptions: { 'Mesh.Algorithm': 1 }, ...meshOpts };
        const tri = new TriBackend(ctx, this, opts);
        await tri.buildMesh(onProgress, shouldStop);   // adaptive refinement passes, emitted via onProgress
        this._modesBackend = tri;
        return tri.solveModes(freq, nev);
    }

    // Resample the sortedIdx-th mode's field from the last solveModes() for plotting.
    getModeField(sortedIdx) {
        return this._modesBackend ? this._modesBackend.getModeField(sortedIdx) : null;
    }

    // Verification certificate (rectilinear backend)
    // Mirror of TriBackend._certifyStatic, so the UI Tolerance means the same thing
    // on both backends: the verified remaining error of the reported static
    // quantities, not the pass-to-pass change the refinement gate measures.
    //
    // Insert a grid line at the midpoint of every x and y interval (uniform bisection
    // of the tensor grid, exact projected size (2nx-1)(2ny-1) ~= 4x nodes per level)
    // and re-solve the statics:
    //   d1 = max rel change of (C, C0) per mode, grid -> bisect
    //   d2 = same, bisect -> bisect^2,   r = d2/d1
    //   certified error = d1/(1−r)
    // Every reported static quantity (C, C0, eps_eff, Z0) is a combination of C and
    // C0, so the static capacitances certify the reported numbers. Losses are not
    // covered (the triangular certificate does not cover them either). Gates match
    // the triangular backend:
    //   * r > rMax  -> pre-asymptotic, cannot certify.
    //   * d1 > tol  -> already failed. Level 2 skipped (cost control).
    //   * d1 < noise floor -> pass outright.
    //   * level-2 grid over its (smaller) node cap -> fall back to r = 0.5. Measured
    //     r on microstrip geometries (single-ended and differential) is 0.28-0.41,
    //     so the fallback overestimates the remaining error (conservative) and the
    //     x1.5 safety additionally covers a true r up to 2/3. Level 2 gets its own
    //     cap well below the level-1 cap because its value (a genuine measured r and
    //     the pre-asymptotic gate) matters most on coarse grids, while on mid-size
    //     grids a 16x-node LU dominates the whole solve time for little sharpening.
    //     A base grid whose bisection is already over the level-1 cap returns null
    //     (the caller keeps the legacy gate). The caps also guard the WASM LU's
    //     ~1 GB allocation limit.
    // The pass decision applies `safety` (x1.5) on top of the estimate; `err` stays
    // the un-inflated best estimate (it is what warnings report).
    //
    // Coarsening (solving on decimated grids, ~0.3x base cost instead of ~6x) was
    // measured as an alternative and rejected: the coarse-level convergence ratio
    // does not transfer to the base level, and the resulting estimate came out 1.2x
    // to 3.5x optimistic vs the bisection reference.
    //
    // Cost reductions that keep the certificate's meaning intact:
    //   * The base level reuses the (C, C0) the refinement pass just computed
    //     (q0 option) instead of re-solving the base grid.
    //   * Odd/even modes share each operator, so every grid level factors the
    //     signal and vacuum matrices once and back-substitutes per mode
    //     (solve_laplace_multi).
    //   * r is measured once per solve and reused by later certificates (knownR
    //     option), deleting the 16x level-2 solve from every call but the first.
    //   * After a failed certificate the loop predicts how many refinement passes
    //     the measured error trend needs before a pass is plausible and skips
    //     certifying until then (certSkipPasses in solve_adaptive).

    // Run fn on a temporarily swapped grid. Every grid-derived array (masks,
    // epsilon_r, conductor ids) is rebuilt from this.x/this.y by _setup_geometry, so
    // swapping the line arrays and rebuilding is a complete state switch, and the
    // finally-rebuild restores the caller's exact state. (this.dx/dy are left
    // untouched, matching _refine_selected_lines, nothing reads them after the
    // initial mesh generation.)
    async _withGrid(x, y, fn) {
        const saved_x = this.x, saved_y = this.y;
        this.x = x;
        this.y = y;
        this._setup_geometry();
        try {
            return await fn();
        } finally {
            this.x = saved_x;
            this.y = saved_y;
            this._setup_geometry();
        }
    }

    // Static (C, C0) per mode — the quantities _solve_single_mode reports, without
    // the loss post-processing the certificate does not cover. All modes share the
    // signal-dielectric operator and the vacuum operator, so each matrix is
    // factored once and back-substituted per mode (solve_laplace_multi): a
    // differential pair does 2 factorizations per grid level instead of 4.
    async _staticCapacitances() {
        const modeNames = this.is_differential ? ['odd', 'even'] : ['single'];
        const signal = await this.solve_laplace_multi(
            modeNames.map(m => this._create_voltage_array(m)), false);
        const vacuum = await this.solve_laplace_multi(
            modeNames.map(m => this._create_voltage_array(m)), true);
        const out = [];
        for (let i = 0; i < modeNames.length; i++) {
            out.push(this._signal_capacitance(signal[i], false));
            out.push(this._signal_capacitance(vacuum[i], true));
        }
        return out;
    }

    // Options beyond the caps/gates:
    //   * q0: base-grid (C, C0) per mode in _staticCapacitances order, when the
    //     caller already has them (the refinement pass that tripped the gate just
    //     computed exactly these). Skips re-solving the base level.
    //   * knownR: convergence ratio measured by an earlier certificate in the same
    //     solve. Skips the 16x level-2 solve entirely: r decreases (or holds) under
    //     refinement, so an earlier genuine measurement used on a finer grid errs
    //     conservative, and the x1.5 safety still covers a moderately larger true r.
    //     Only genuinely measured, non-pre-asymptotic r values are reused.
    async _certifyStatic(tol, { maxNodes = 600000, l2MaxNodes = 150000, rMax = 0.7, safety = 1.5,
                                q0 = null, knownR = null } = {}) {
        const maxDiff = (a, b) => {
            let m = 0;
            for (let i = 0; i < a.length; i++)
                m = Math.max(m, Math.abs(a[i] - b[i]) / Math.max(Math.abs(b[i]), 1e-300));
            return m;
        };
        const bisect = (arr) => {
            const out = new Float64Array(2 * arr.length - 1);
            for (let i = 0; i < arr.length - 1; i++) {
                out[2 * i] = arr[i];
                out[2 * i + 1] = 0.5 * (arr[i] + arr[i + 1]);
            }
            out[out.length - 1] = arr[arr.length - 1];
            return out;
        };
        const bisectedNodes = (x, y) => (2 * x.length - 1) * (2 * y.length - 1);
        if (bisectedNodes(this.x, this.y) > maxNodes) return null;
        if (!q0) q0 = await this._staticCapacitances();
        const x1 = bisect(this.x), y1 = bisect(this.y);
        const q1 = await this._withGrid(x1, y1, () => this._staticCapacitances());
        const d1 = maxDiff(q0, q1);
        const base = { nodes: this.x.length * this.y.length, d1, safety };
        if (d1 < 5e-5) return { ...base, pass: true, err: d1, r: 0, levels: 1 };
        if (d1 * safety >= tol) return { ...base, pass: false, err: d1, r: null, levels: 1 };
        let r = 0.5, levels = 1, rSource = 'fallback';
        if (knownR != null) {
            r = knownR;
            rSource = 'reused';
        } else if (bisectedNodes(x1, y1) <= l2MaxNodes) {
            const x2 = bisect(x1), y2 = bisect(y1);
            const q2 = await this._withGrid(x2, y2, () => this._staticCapacitances());
            r = maxDiff(q1, q2) / d1;
            levels = 2;
            rSource = 'measured';
            if (r >= rMax)
                return { ...base, pass: false, err: d1, r, levels, rSource, preAsymptotic: true };
        }
        const err = d1 / (1 - Math.min(r, rMax));
        return { ...base, pass: err * safety < tol, err, r, levels, rSource };
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
         *   Z_common: (only for differential) Z_even / 2,
         *   RLGC_matrix: (only for differential) {
         *     R: [[R11, R12], [R21, R22]],  // Resistance matrix (Ohm/m)
         *     L: [[L11, L12], [L21, L22]],  // Inductance matrix (H/m)
         *     G: [[G11, G12], [G21, G22]],  // Conductance matrix (S/m)
         *     C: [[C11, C12], [C21, C22]]   // Capacitance matrix (F/m)
         *   }
         * }
         *
         * Note: For differential pairs, RLGC_matrix represents the physical 2x2 per-unit-length
         * parameter matrices relating the voltages and currents on the two traces. The diagonal
         * elements (11, 22) are self-parameters, and off-diagonal elements (12, 21) are mutual
         * coupling parameters. For L and C, coupling terms are negative.
         */
        // Reject geometries that can't be meshed finely enough to resolve features and
        // the wavelength while staying under the node budget, before either backend
        // starts building a (possibly enormous) mesh.
        this._check_meshability(options.max_nodes ?? 20000);

        // Triangular FEM backend: delegate the whole solve (mesh + static +
        // full-wave eigenmode + loss) to TriBackend, lazily loaded so the gmsh
        // WASM is only fetched when this backend is selected.
        if (this.mesh_backend === 'triangular') {
            // Wire the UI adaptive controls (max iterations, tolerance, max nodes) into
            // the triangular backend's refinement loop; buildMesh reports each pass via
            // onProgress. Undefined values fall back to the backend defaults.
            const triOpts = {
                maxRefineIters: options.max_iters,
                refineTol: options.energy_tol,
                maxNodes: options.max_nodes,
                minConvergedPasses: options.min_converged_passes,
            };
            // Only forward certify when the caller decided it (the UI "Estimate
            // solution error" checkbox): an unconditional `certify: undefined` would
            // shadow a certify set through tri_opts.
            if (options.certify !== undefined) triOpts.certify = options.certify;
            const tri = await this._ensureTriBackend(options.onProgress, triOpts, options.shouldStop);
            return tri.solveAt(this.freq);
        }

        // Ensure mesh is generated
        if (this.ensure_mesh) {
            this.ensure_mesh();
        }

        const {
            max_iters = 10,
            refine_frac,
            energy_tol = 0.01,
            param_tol = 0.1,
            max_nodes = 20000,
            min_converged_passes = 1,
            certify = true,
            certify_max_nodes = 600000,
            certify_l2_max_nodes = 150000,
            onProgress = null,
            shouldStop = null,
            skip_mesh = false
        } = options;

        // If skip_mesh is true, just solve once with existing mesh
        if (skip_mesh) {
            let modeResults;
            if (this.is_differential) {
                modeResults = [await this._solve_single_mode('odd', true), await this._solve_single_mode('even', true)];
                const modal = await this._solve_modal_differential();   // null for symmetric/degenerate
                if (modal) modeResults = modal;
            } else {
                modeResults = [await this._solve_single_mode('single', true)];
            }
            // Store fields as arrays
            this.V = modeResults.map(r => r.V);
            this.Ex = modeResults.map(r => r.Ex);
            this.Ey = modeResults.map(r => r.Ey);
            return this._build_results(modeResults);
        }

        // Set default refine_frac based on mode
        const refineFrac = refine_frac !== undefined ? refine_frac : (this.is_differential ? 0.15 : 0.2);

        // Verification certificate state (see _certifyStatic). Mirrors the triangular
        // backend: a tripped pass-to-pass gate is only a CANDIDATE stop until the
        // certified remaining error passes the tolerance.
        const certOpts = { maxNodes: certify_max_nodes, l2MaxNodes: certify_l2_max_nodes };
        this.certification = null;
        this._certWarn = null;

        // The refinement pass that trips the gate has already computed exactly the
        // quantities the certificate compares against (C, C0 per mode via
        // _signal_capacitance), so the certificate's base level is free.
        const certQ0 = () => modeResults ? modeResults.flatMap(r => [r.C, r.C0]) : null;

        // Predictive re-certification: after a failed certificate the error must
        // still shrink by err*safety/tol before a pass is possible, so certifying
        // again at the very next gate trip mostly wastes a 4x-node solve. Estimate
        // the per-pass error shrink (measured from consecutive failed certificates
        // when available, else assume 30% per adaptive pass) and hold off
        // certifying for the predicted number of passes, capped so a bad estimate
        // cannot postpone convergence detection for long. The post-loop
        // certification is unconditional, so a hold can never lose the final
        // verdict, only defer it.
        const certSkipPasses = (cert, prevFail, itNow) => {
            const need = (cert.err * (cert.safety ?? 1.5)) / energy_tol;
            if (!(need > 1)) return 1;
            let rho = 0.7;
            if (prevFail && cert.err < prevFail.err && itNow > prevFail.it) {
                rho = Math.pow(cert.err / prevFail.err, 1 / (itNow - prevFail.it));
                rho = Math.min(Math.max(rho, 0.3), 0.95);
            }
            return Math.min(Math.max(Math.ceil(Math.log(need) / Math.log(1 / rho)), 1), 4);
        };
        let certR = null;          // measured convergence ratio, reused by later certificates
        let certHoldUntil = -1;    // certification is skipped while it < certHoldUntil
        let lastFailedCert = null;

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
                    if (converged_count >= min_converged_passes && it >= certHoldUntil) {
                        // The pass-to-pass gate above measures the RATE of approach,
                        // not the absolute error. Certify the actual remaining error
                        // (see _certifyStatic) and keep refining if the certificate
                        // fails. A null certificate (grid over the certification cap)
                        // keeps the legacy behavior.
                        let cert = null;
                        if (certify) {
                            try {
                                cert = await this._certifyStatic(energy_tol,
                                    { ...certOpts, q0: certQ0(), knownR: certR });
                            } catch { cert = null; }
                            if (cert) {
                                this.certification = cert;
                                if (cert.rSource === 'measured' && !cert.preAsymptotic)
                                    certR = cert.r;
                            }
                        }
                        if (!certify || !cert || cert.pass) {
                            console.log(`Converged after ${it + 1} passes`);
                            break;
                        }
                        // gate was optimistic, keep refining and hold off the next
                        // certificate until the predicted refinement has accumulated
                        certHoldUntil = it + certSkipPasses(cert, lastFailedCert, it);
                        lastFailedCert = { err: cert.err, it };
                        converged_count = 0;
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

            // Refine mesh using combined fields from all modes so that regions
            // important for any mode (e.g. even mode near ground planes) get refined.
            if (it !== max_iters - 1) {
                if (modeResults.length > 1) {
                    this.refine_mesh_multi(modeResults, refineFrac);
                } else {
                    const refineMode = modeResults[0];
                    this.refine_mesh(refineMode.V, refineMode.Ex, refineMode.Ey, refineFrac);
                }
                this._setup_geometry();
            }
        }

        // Accuracy report: if refinement ended without a passing certificate
        // (iteration cap, node budget), measure the final grid once so the user gets
        // an estimated error. A certificate that already covered this grid (pass or
        // fail) is not recomputed. Mirrors TriBackend.buildMesh.
        if (certify) {
            let cert = this.certification;
            const nNodes = this.x.length * this.y.length;
            if (!cert || (!cert.pass && cert.nodes !== nNodes)) {
                try {
                    cert = await this._certifyStatic(energy_tol,
                        { ...certOpts, q0: certQ0(), knownR: certR });
                    this.certification = cert;
                }
                catch { /* keep whatever we had */ }
            }
            if (cert && !cert.pass) {
                // Failure regimes: pre-asymptotic (estimate is only a lower
                // bound), estimate over the tolerance, and estimate under the
                // tolerance but without the safety margin the certificate
                // requires (saying "stopped before reaching 1%: error is 0.7%"
                // would read as a contradiction).
                const pct = (x) => (100 * x).toFixed(2);
                const estStr = cert.preAsymptotic
                    ? `at least ${pct(cert.err)}% (grid still pre-asymptotic, the estimate is a lower bound)`
                    : cert.err < energy_tol
                        ? `about ${pct(cert.err)}%, within the tolerance but without the ` +
                          `×${cert.safety ?? 1.5} margin needed to certify it`
                        : `about ${pct(cert.err)}%`;
                this._certWarn = { type: 'accuracy', mode: 'all', message:
                    `Quasi-static mesh refinement could not verify the requested tolerance ` +
                    `(${pct(energy_tol)}%): estimated remaining error is ${estStr}. ` +
                    `Increase Max Nodes / Max Iterations, or relax Tolerance.` };
            }
        }

        // For an ASYMMETRIC differential pair, replace the odd/even drive results with the
        // genuine two-conductor modal decomposition on the converged mesh. _solve_modal_differential
        // returns null for a symmetric or degenerate pair, in which case the odd/even results
        // (already in modeResults) are exact and are kept.
        if (this.is_differential) {
            const modal = await this._solve_modal_differential();
            if (modal) modeResults = modal;
        }

        // Store fields as arrays
        this.V = modeResults.map(r => r.V);
        this.Ex = modeResults.map(r => r.Ex);
        this.Ey = modeResults.map(r => r.Ey);

        // Build unified result structure
        return this._build_results(modeResults);
    }

    _modal_to_physical_rlgc(odd, even) {
        /**
         * Physical 2×2 RLGC matrices for the differential pair (relating the per-unit-length
         * voltages/currents on the two traces: V1 = Z11·I1 + Z12·I2, V2 = Z21·I1 + Z22·I2).
         *
         * @param {object} odd - Odd mode results with RLGC
         * @param {object} even - Even mode results with RLGC
         * @returns {object} - Physical 2x2 matrices { R, L, G, C }
         */
        // Delegate to the shared builder: the genuine asymmetric [C]/[L] (with Tv-reconstructed
        // [R]/[G]) when an asymmetric physMatrix was computed for this pair, else the symmetric
        // odd/even reconstruction below. Sharing the builder with the S-parameter path guarantees
        // RLGC_matrix matches the matrices that actually drive the S-parameters.
        //
        //   Symmetric reconstruction: X11 = X22 = (X_odd + X_even)/2,  X12 = X21 = (X_even − X_odd)/2
        //   - Odd mode = opposite-polarity (differential) drive; even mode = same-polarity (common).
        //   - L12 = (L_even − L_odd)/2 > 0 (positive mutual L); C12 = (C_even − C_odd)/2 < 0.
        // For an asymmetric pair the diagonal self terms differ (C11 ≠ C22) and come from the
        // physical Maxwell matrix rather than the odd/even average.
        return buildPhysicalRLGC(odd.RLGC, even.RLGC, this._modalPhys);
    }

    _build_results(modeResults) {
        /**
         * Build the unified result structure from mode results.
         */
        const result = { modes: modeResults };

        // Failed verification certificate: surface as an accuracy warning on every
        // result built from this grid, the same per-solve channel the triangular
        // backend uses (result.warnings / solver.modeWarnings), so the UI logs it for
        // the adaptive solve and for every sweep point alike.
        if (this._certWarn) {
            result.warnings = [this._certWarn];
            this.modeWarnings = result.warnings;
        }

        if (this.is_differential) {
            const odd = modeResults.find(m => m.mode === 'odd');
            const even = modeResults.find(m => m.mode === 'even');
            result.Z_diff = 2 * odd.Z0;
            result.Z_common = even.Z0 / 2;

            // Add physical 2x2 RLGC matrix
            result.RLGC_matrix = this._modal_to_physical_rlgc(odd, even);
            // True physical 2×2 [C]/[L] for the asymmetric MTL 4-port S-parameter path.
            if (this._modalPhys) result.physMatrix = this._modalPhys;
        }

        return result;
    }

    /**
     * Perform a frequency sweep with automatic mesh generation at optimal frequency.
     * This is the recommended single-entry-point API for frequency sweeps.
     *
     * @param {object} options - Sweep configuration
     * @param {number[]} options.frequencies - Array of frequencies in Hz
     * @param {number} [options.energy_tol=0.02] - Energy convergence tolerance for adaptive mesh
     * @param {number} [options.max_nodes=20000] - Maximum mesh nodes
     * @param {function} [options.onProgress] - Progress callback
     * @param {function} [options.shouldStop] - Stop check callback
     * @returns {Promise<object>} - Results organized for plotting:
     *   {
     *     frequencies: [...],
     *     modes: [{
     *       mode: 'single'|'odd'|'even',
     *       Z0: [...], Zc_re: [...], Zc_im: [...],
     *       eps_eff: [...],
     *       alpha_c: [...], alpha_d: [...], alpha_total: [...],
     *       RLGC: { R: [...], L: [...], G: [...], C: [...] },  // modal parameters
     *       C: number, C0: number  // static values
     *     }, ...],
     *     Z_diff: [...],    // differential only
     *     Z_common: [...],  // differential only
     *     RLGC_matrix: {    // differential only - physical 2x2 matrices
     *       R: { R11: [...], R12: [...], R21: [...], R22: [...] },
     *       L: { L11: [...], L12: [...], L21: [...], L22: [...] },
     *       G: { G11: [...], G12: [...], G21: [...], G22: [...] },
     *       C: { C11: [...], C12: [...], C21: [...], C22: [...] }
     *     },
     *     mesh: { nx, ny }
     *   }
     */
    async solve_sweep(options = {}) {
        const {
            frequencies,
            energy_tol = 0.02,
            max_nodes = 20000,
            onProgress = null,
            shouldStop = null
        } = options;

        // Validate frequencies
        if (!frequencies || !Array.isArray(frequencies) || frequencies.length === 0) {
            throw new Error('frequencies must be a non-empty array');
        }

        // Sort frequencies and find max for optimal meshing
        const sortedFreqs = [...frequencies].sort((a, b) => a - b);
        const maxFreq = sortedFreqs[sortedFreqs.length - 1];

        // Set frequency to max for finest skin depth mesh
        this.freq = maxFreq;
        // Sweep-max hint for the triangular backend: lets it size the MQS skin
        // band for the WHOLE sweep up front, so an ascending sweep reuses one
        // skin mesh (and its cached assembly) instead of rebuilding per point.
        this._sweepFmax = maxFreq;

        // Fail fast (before building the FDM mesh below) if the geometry can't be meshed
        // finely enough to resolve features and the wavelength within the node budget.
        this._check_meshability(max_nodes);

        // Force mesh regeneration
        this.mesh_generated = false;

        // Generate mesh and run adaptive refinement
        if (this.ensure_mesh) {
            this.ensure_mesh();
        }

        const initResult = await this.solve_adaptive({
            energy_tol,
            max_nodes,
            onProgress,
            shouldStop
        });

        // Initialize result arrays
        const modeNames = this.is_differential ? ['odd', 'even'] : ['single'];
        const resultModes = modeNames.map(modeName => {
            const initMode = initResult.modes.find(m => m.mode === modeName);
            return {
                mode: modeName,
                Z0: [],
                Zc_re: [],
                Zc_im: [],
                eps_eff: [],
                alpha_c: [],
                alpha_d: [],
                alpha_total: [],
                RLGC: { R: [], L: [], G: [], C: [] },
                C: initMode.C,
                C0: initMode.C0
            };
        });

        const result = {
            frequencies: [],
            modes: resultModes,
            mesh: { nx: this.x.length, ny: this.y.length }
        };

        if (this.is_differential) {
            result.Z_diff = [];
            result.Z_common = [];
            // Initialize RLGC_matrix arrays for 2x2 physical matrices
            result.RLGC_matrix = {
                R: { R11: [], R12: [], R21: [], R22: [] },
                L: { L11: [], L12: [], L21: [], L22: [] },
                G: { G11: [], G12: [], G21: [], G22: [] },
                C: { C11: [], C12: [], C21: [], C22: [] }
            };
        }

        // Compute at each frequency
        for (const freq of sortedFreqs) {
            const freqResult = await this.computeAtFrequency(freq, initResult);

            result.frequencies.push(freq);

            // Extract mode results
            for (let i = 0; i < modeNames.length; i++) {
                const modeName = modeNames[i];
                const modeResult = freqResult.modes.find(m => m.mode === modeName);
                const outMode = resultModes[i];

                outMode.Z0.push(modeResult.Z0);
                outMode.Zc_re.push(modeResult.Zc.re);
                outMode.Zc_im.push(modeResult.Zc.im);
                outMode.eps_eff.push(modeResult.eps_eff);
                outMode.alpha_c.push(modeResult.alpha_c);
                outMode.alpha_d.push(modeResult.alpha_d);
                outMode.alpha_total.push(modeResult.alpha_total);
                outMode.RLGC.R.push(modeResult.RLGC.R);
                outMode.RLGC.L.push(modeResult.RLGC.L);
                outMode.RLGC.G.push(modeResult.RLGC.G);
                outMode.RLGC.C.push(modeResult.RLGC.C);
            }

            // Differential-specific results
            if (this.is_differential) {
                result.Z_diff.push(freqResult.Z_diff);
                result.Z_common.push(freqResult.Z_common);

                // Add physical 2x2 RLGC matrix values
                const rlgc_mat = freqResult.RLGC_matrix;
                result.RLGC_matrix.R.R11.push(rlgc_mat.R[0][0]);
                result.RLGC_matrix.R.R12.push(rlgc_mat.R[0][1]);
                result.RLGC_matrix.R.R21.push(rlgc_mat.R[1][0]);
                result.RLGC_matrix.R.R22.push(rlgc_mat.R[1][1]);

                result.RLGC_matrix.L.L11.push(rlgc_mat.L[0][0]);
                result.RLGC_matrix.L.L12.push(rlgc_mat.L[0][1]);
                result.RLGC_matrix.L.L21.push(rlgc_mat.L[1][0]);
                result.RLGC_matrix.L.L22.push(rlgc_mat.L[1][1]);

                result.RLGC_matrix.G.G11.push(rlgc_mat.G[0][0]);
                result.RLGC_matrix.G.G12.push(rlgc_mat.G[0][1]);
                result.RLGC_matrix.G.G21.push(rlgc_mat.G[1][0]);
                result.RLGC_matrix.G.G22.push(rlgc_mat.G[1][1]);

                result.RLGC_matrix.C.C11.push(rlgc_mat.C[0][0]);
                result.RLGC_matrix.C.C12.push(rlgc_mat.C[0][1]);
                result.RLGC_matrix.C.C21.push(rlgc_mat.C[1][0]);
                result.RLGC_matrix.C.C22.push(rlgc_mat.C[1][1]);
            }
        }

        // Drop the sweep-max hint: a later single-frequency solve on this same
        // solver shouldn't keep building the whole-sweep skin band.
        this._sweepFmax = null;
        return result;
    }

    /**
     * Compute frequency-dependent results using cached fields.
     * This is a fast path for frequency sweeps where only frequency changes,
     * not the geometry or dielectric distribution.
     *
     * @param {number} freq - Frequency in Hz
     * @param {object} cachedResults - Results from a previous solve containing V, Ex, Ey, C, C0, Z0
     * @returns {object} - New results with updated frequency-dependent parameters
     */
    async computeAtFrequency(freq, cachedResults) {
        // Update frequency
        this.freq = freq;

        // Triangular FEM backend: re-run the per-frequency solve (eigenmode +
        // loss) on the cached mesh/static solution. skipFieldResample: sweep
        // points don't need the plot-field resample (the displayed fields come
        // from the main solve; per-point resamples are overwritten anyway).
        if (this.mesh_backend === 'triangular') {
            const tri = await this._ensureTriBackend();
            return tri.solveAt(freq, { skipFieldResample: true });
        }

        // If causal materials are enabled, we must re-solve the Laplace equation
        // because epsilon_r changes with frequency, which changes the field distribution
        if (this.use_causal_materials) {
            // Apply the causal model to update epsilon_r and tand
            applyDjordjevicSarkar(this);

            // Re-solve at this frequency with updated material parameters
            const modeResults = [];

            if (this.is_differential) {
                // ASYMMETRIC pair: the odd/even drive fields are mode MIXTURES, not the
                // line's modes — use the modal decomposition under the causal materials,
                // exactly as solve_adaptive does on the converged mesh (it also refreshes
                // _modalPhys, so the 4-port S-matrix tracks the causal ε). Returns null
                // for a velocity-degenerate pair → keep the drive results. Gated on
                // _modalPhys from the initial solve: a symmetric pair (null) skips the
                // four extra Laplace solves — its odd/even drives are exact.
                const modal = this._modalPhys ? await this._solve_modal_differential() : null;
                if (modal) {
                    modeResults.push(...modal);
                } else {
                // Solve both odd and even modes
                const oddMode = await this._solve_single_mode('odd', false);
                const evenMode = await this._solve_single_mode('even', false);

                // Use cached C0 values from initial solve (vacuum doesn't
                // change).  Same for the vacuum fields: permittivity- and
                // frequency-independent, so the causal re-solve (which only
                // shifts the dielectric fields) reuses them for the
                // conductor-loss integrand.
                const cachedOdd = cachedResults.modes.find(m => m.mode === 'odd');
                const cachedEven = cachedResults.modes.find(m => m.mode === 'even');
                oddMode.C0 = cachedOdd.C0;
                evenMode.C0 = cachedEven.C0;
                oddMode.Ex0 = cachedOdd.Ex0; oddMode.Ey0 = cachedOdd.Ey0;
                evenMode.Ex0 = cachedEven.Ex0; evenMode.Ey0 = cachedEven.Ey0;

                // Recalculate eps_eff and Z0 with new C and cached C0
                oddMode.eps_eff = oddMode.C / oddMode.C0;
                oddMode.Z0 = 1 / (CONSTANTS.C * Math.sqrt(oddMode.C * oddMode.C0));
                evenMode.eps_eff = evenMode.C / evenMode.C0;
                evenMode.Z0 = 1 / (CONSTANTS.C * Math.sqrt(evenMode.C * evenMode.C0));

                // Recalculate RLGC parameters with corrected Z0
                const recalc = (mode) => {
                    const { R_ac, R_dc, R_total, L_internal } = this._mode_conductor_loss(
                        mode.Ex, mode.Ey, mode.Z0, mode.C0, mode.Ex0, mode.Ey0);
                    const alpha_d = this.calculate_dielectric_loss(mode.Ex, mode.Ey, mode.Z0);
                    const { Zc, rlgc, eps_eff_mode, L_external } = this.rlgc(R_total, L_internal, alpha_d, mode.C, mode.Z0);
                    mode.RLGC = rlgc;
                    mode.Zc = Zc;
                    mode.eps_eff = eps_eff_mode;
                    mode.alpha_c = 8.686 * R_total / (2 * Zc.re);
                    mode.alpha_d = alpha_d;
                    mode.alpha_total = mode.alpha_c + alpha_d;
                    mode.L_internal = L_internal;
                    mode.L_external = L_external;
                };

                recalc(oddMode);
                recalc(evenMode);

                modeResults.push(oddMode, evenMode);
                }
            } else {
                // Solve single mode
                const result = await this._solve_single_mode('single', false);

                // Use cached C0 from initial solve and the cached vacuum fields
                // (ε- and frequency-independent) for the conductor-loss integrand.
                result.C0 = cachedResults.modes[0].C0;
                result.Ex0 = cachedResults.modes[0].Ex0;
                result.Ey0 = cachedResults.modes[0].Ey0;

                // Recalculate eps_eff and Z0 with new C and cached C0
                result.eps_eff = result.C / result.C0;
                result.Z0 = 1 / (CONSTANTS.C * Math.sqrt(result.C * result.C0));

                // Recalculate RLGC parameters with corrected Z0
                const { R_ac, R_dc, R_total, L_internal } = this._mode_conductor_loss(
                    result.Ex, result.Ey, result.Z0, result.C0, result.Ex0, result.Ey0);
                const alpha_d = this.calculate_dielectric_loss(result.Ex, result.Ey, result.Z0);
                const { Zc, rlgc, eps_eff_mode, L_external } = this.rlgc(R_total, L_internal, alpha_d, result.C, result.Z0);

                result.RLGC = rlgc;
                result.Zc = Zc;
                result.eps_eff = eps_eff_mode;
                result.alpha_c = 8.686 * R_total / (2 * Zc.re);
                result.alpha_d = alpha_d;
                result.alpha_total = result.alpha_c + alpha_d;
                result.L_internal = L_internal;
                result.L_external = L_external;

                modeResults.push(result);
            }

            return this._build_results(modeResults);
        }

        // Fast path: Non-causal materials - use cached fields
        const modeResults = [];

        for (const cached of cachedResults.modes) {
            const { mode, V, Ex, Ey, Ex0, Ey0, C, C0, Z0 } = cached;

            // Recalculate conductor losses with new frequency (affects skin depth)
            const { R_ac, R_dc, R_total, L_internal } = this._mode_conductor_loss(Ex, Ey, Z0, C0, Ex0, Ey0);

            // Recalculate dielectric loss (affects omega)
            const alpha_d = this.calculate_dielectric_loss(Ex, Ey, Z0);

            // Recalculate RLGC with new frequency
            const { Zc, rlgc, eps_eff_mode, L_external } = this.rlgc(R_total, L_internal, alpha_d, C, Z0);

            // Calculate conductor loss alpha from R_total
            const alpha_c = 8.686 * R_total / (2 * Zc.re);
            const alpha_total = alpha_c + alpha_d;

            modeResults.push({
                mode,
                Z0,
                eps_eff: eps_eff_mode,
                C, C0,
                RLGC: rlgc, Zc,
                alpha_c, alpha_d, alpha_total,
                L_internal, L_external,
                V, Ex, Ey,
                Ex0, Ey0
            });
        }

        return this._build_results(modeResults);
    }
}
