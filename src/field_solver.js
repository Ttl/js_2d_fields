import createWASMModule from './wasm_solver/solver.js';
import { Complex } from "./complex.js";
import { calculate_Zrough, calculate_Zrough_layered, wallSpreadFactor } from './surface_roughness.js';
import { applyDjordjevicSarkar } from './djordjevic_sarkar.js';
import { classifyModalDecomposition, conductorFinishKey } from './geometry_symmetry.js';
import { buildPhysicalRLGC } from './sparameters.js';
import { visibleAreas, platedThrough, platingArea, insideRingHole, bodyDistance } from './shapes.js';
import { unlimitedGrounds } from './wall_grounds.js';

export const CONSTANTS = {
    EPS0: 8.854187817e-12,
    MU0: 4 * Math.PI * 1e-7,
    C: 299792458,
    PI: Math.PI
};

// Node cap for the certificate's convergence-ratio (level-2) measurement, on
// the peak grid it solves (see _certifyStatic). It runs at most once per
// adaptive solve (knownR caches the result for every later certificate) and
// only when the conservative fallback ratio is what stands between the current
// grid and a pass, so it can afford a grid well above the level-1 working size.
const CERT_L2_MAX_NODES = 400000;

// --- Math Utils ---

export function diff(arr) {
    const res = new Float64Array(arr.length - 1);
    for (let i = 0; i < arr.length - 1; i++) res[i] = arr[i+1] - arr[i];
    return res;
}

// Store the initialized WASM module (singleton pattern)
let WASMModuleInstance = null;

// Solve one matrix against several right-hand sides with a single
// factorization. The factorization dominates the solve cost, so k systems that
// share an operator (e.g. the odd and even modes of a differential pair, which
// differ only in their Dirichlet values) cost barely more than one.
//
// forceLU forces the LU path. The default lets the WASM side pick, which means
// Cholesky whenever the matrix tests symmetric (the FDM Laplace operator always
// is, half-domain included. See solve_laplace_multi row scaling).
async function solveWithWASMMulti(csr, Bs, forceLU = false) {
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
            forceLU ? 1 : 0
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

async function solveWithWASM(csr, B, forceLU = false) {
    return (await solveWithWASMMulti(csr, [B], forceLU))[0];
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

    // Bottom of the solve domain. Rectangular line types put their bottom ground
    // below y=0, so it defaults to -t_gnd. Solvers without a bottom ground set it.
    get domain_y_min() { return this._domain_y_min ?? -(this.t_gnd ?? 0); }
    set domain_y_min(v) { this._domain_y_min = v; }

    ensure_mesh() {
        // Triangular backend builds & owns its own mesh in TriBackend.
        if (this.mesh_backend === 'triangular') return;
        if (this.mesh_generated) {
            return;
        }

        // Generate mesh
        [this.x, this.y] = this.mesher.generate_mesh();

        // Calculate spacing arrays
        this.dx = new Float64Array(this.x.length - 1);
        for (let i = 0; i < this.x.length - 1; i++) {
            this.dx[i] = this.x[i + 1] - this.x[i];
        }

        this.dy = new Float64Array(this.y.length - 1);
        for (let i = 0; i < this.y.length - 1; i++) {
            this.dy[i] = this.y[i + 1] - this.y[i];
        }

        // Setup geometry
        this._setup_geometry();
        this.mesh_generated = true;
    }

    _setup_geometry() {
        const tol = 1e-11;
        const nx = this.x.length;
        const ny = this.y.length;

        // Initialize mask and material arrays
        // Nodes no dielectric rect covers are air: lossless, epsilon_r = 1.
        this.epsilon_r = Array(ny).fill().map(() => new Float64Array(nx).fill(1));
        this.tand = Array(ny).fill().map(() => new Float64Array(nx).fill(0));
        this.signal_mask = Array(ny).fill().map(() => new Uint8Array(nx));
        this.ground_mask = Array(ny).fill().map(() => new Uint8Array(nx));

        // For differential mode, track positive and negative traces separately
        if (this.is_differential) {
            this.signal_p_mask = Array(ny).fill().map(() => new Uint8Array(nx));
            this.signal_n_mask = Array(ny).fill().map(() => new Uint8Array(nx));
        }

        // Apply dielectrics (last overwrites)
        for (const diel of this.dielectrics) {
            for (let i = 0; i < ny; i++) {
                const yc = this.y[i];
                if (yc >= diel.y_min - tol && yc <= diel.y_max + tol) {
                    for (let j = 0; j < nx; j++) {
                        const xc = this.x[j];
                        if (xc >= diel.x_min - tol && xc <= diel.x_max + tol) {
                            this.epsilon_r[i][j] = diel.epsilon_r;
                            this.tand[i][j] = diel.tan_delta;
                        }
                    }
                }
            }
        }

        // Apply conductors (also build conductor_id map for per-conductor loss calc)
        this.conductor_id = Array(ny).fill().map(() => new Int16Array(nx).fill(-1));
        for (let ci = 0; ci < this.conductors.length; ci++) {
            const cond = this.conductors[ci];
            for (let i = 0; i < ny; i++) {
                const yc = this.y[i];
                if (yc >= cond.y_min - tol && yc <= cond.y_max + tol) {
                    for (let j = 0; j < nx; j++) {
                        const xc = this.x[j];
                        if (xc >= cond.x_min - tol && xc <= cond.x_max + tol) {
                            this.conductor_id[i][j] = ci;
                            if (cond.is_signal) {
                                this.signal_mask[i][j] = 1;
                                if (this.is_differential) {
                                    if (cond.polarity > 0) {
                                        this.signal_p_mask[i][j] = 1;
                                    } else {
                                        this.signal_n_mask[i][j] = 1;
                                    }
                                }
                            } else {
                                this.ground_mask[i][j] = 1;
                            }
                        }
                    }
                }
            }
        }

        // Finalize conductor mask
        this.conductor_mask = Array(ny).fill().map((_, i) => {
            const row = new Uint8Array(nx);
            for (let j = 0; j < nx; j++) {
                row[j] = this.signal_mask[i][j] | this.ground_mask[i][j];
            }
            return row;
        });

        // Cell-centred permittivity for the flux coefficients (see
        // FieldSolver2D._paint_cell_materials).
        this._paint_cell_materials();
    }

    /**
     * Cell-centred material arrays: epsilon_cell[i][j] / tand_cell[i][j] describe
     * the rectangle [x[j], x[j+1]] x [y[i], y[i+1]], sampled at its centre.
     *
     * The nodal epsilon_r samples material at grid points, which
     * leaves a dielectric interface ambiguous. A node sitting exactly on the
     * interface (which is where the mesher puts its lines) belongs to whichever
     * region was painted last, and averaging two nodal values for a face
     * permittivity then smears the interface by half a cell. A cell never
     * straddles an interface that lies on a grid line, so the cell array is
     * unambiguous, and the face coefficients built from it (see
     * solve_laplace_multi) are exact.
     *
     * Called by each solver's _setup_geometry after the nodal arrays are painted.
     */
    _paint_cell_materials() {
        const nx = this.x.length, ny = this.y.length;
        const nxc = Math.max(nx - 1, 0), nyc = Math.max(ny - 1, 0);
        const eps = Array(nyc).fill().map(() => new Float64Array(nxc).fill(1));
        const tand = Array(nyc).fill().map(() => new Float64Array(nxc));
        const xc = new Float64Array(nxc), yc = new Float64Array(nyc);
        for (let j = 0; j < nxc; j++) xc[j] = 0.5 * (this.x[j] + this.x[j + 1]);
        for (let i = 0; i < nyc; i++) yc[i] = 0.5 * (this.y[i] + this.y[i + 1]);
        // Same "later dielectric wins" order as the nodal painting, so a stack
        // (substrate, then solder mask over it) resolves identically.
        for (const d of (this.dielectrics || [])) {
            for (let i = 0; i < nyc; i++) {
                if (yc[i] < d.y_min || yc[i] > d.y_max) continue;
                const er = eps[i], td = tand[i];
                for (let j = 0; j < nxc; j++) {
                    if (xc[j] < d.x_min || xc[j] > d.x_max) continue;
                    er[j] = d.epsilon_r;
                    td[j] = d.tan_delta;
                }
            }
        }
        this.epsilon_cell = eps;
        this.tand_cell = tand;
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
    // clearance: false drops this geometric rule, for callers that go on to solve and
    // get the measured check (openBoundaryFieldWarning) instead.
    openBoundaryWarnings({ clearance = true } = {}) {
        const out = [];
        const b = this.boundaries;
        if (!clearance || !b || !this.conductors || !this.conductors.length) return out;
        const xMin = -this.domain_width / 2, xMax = this.domain_width / 2;
        const yMin = this.domain_y_min, yMax = this.domain_height;
        if (!(xMax > xMin) || !(yMax > yMin)) return out;
        // Decay scale: substrate stack thickness (non-air dielectrics); if there is
        // none (all-air line), fall back to the conductor stack height.
        // A dielectric that runs into an open top or bottom wall is exterior fill, not
        // a substrate, and sets no decay scale.
        let lo = Infinity, hi = -Infinity;
        const yTol = (yMax - yMin) * 1e-9;
        for (const d of (this.dielectrics || [])) {
            if ((d.epsilon_r || 1) <= 1.001) continue;
            if ((b[2] === 'open' && d.y_max >= yMax - yTol) || (b[3] === 'open' && d.y_min <= yMin + yTol)) continue;
            lo = Math.min(lo, d.y_min); hi = Math.max(hi, d.y_max);
        }
        if (!(hi > lo)) {
            for (const c of this.conductors) { lo = Math.min(lo, c.y_min); hi = Math.max(hi, c.y_max); }
        }
        const h = hi - lo;
        if (!(h > 0)) return out;
        // A ground ring around every signal conductor (a coax shield) keeps the field
        // off all four walls.
        const signals = this.conductors.filter(c => c.is_signal);
        if (this.conductors.some(g => !g.is_signal && g.shape && g.shape.type === 'ring'
            && signals.every(c => insideRingHole(g.shape, c)))) return out;
        const OPEN_CLEARANCE = 3;
        const tol = Math.max(xMax - xMin, yMax - yMin) * 1e-9;
        // Each wall: [name, bc, distance-to-wall, "g lies strictly between c and
        // this wall", "g covers c's extent along the wall direction"]. The last
        // two implement shielding: a grounded conductor that touches the wall and
        // sits between a conductor and that wall (GCPW), screens the fields and
        // the open boundary behind it is fine.
        const walls = [
            ['left', b[0], c => c.x_min - xMin,
                (g, c) => g.x_max <= c.x_min + tol, (g, c) => g.y_min <= c.y_min + tol && g.y_max >= c.y_max - tol],
            ['right', b[1], c => xMax - c.x_max,
                (g, c) => g.x_min >= c.x_max - tol, (g, c) => g.y_min <= c.y_min + tol && g.y_max >= c.y_max - tol],
            ['top', b[2], c => yMax - c.y_max,
                (g, c) => g.y_min >= c.y_max - tol, (g, c) => g.x_min <= c.x_min + tol && g.x_max >= c.x_max - tol],
            ['bottom', b[3], c => c.y_min - yMin,
                (g, c) => g.y_max <= c.y_min + tol, (g, c) => g.x_min <= c.x_min + tol && g.x_max >= c.x_max - tol],
        ];
        const mm = (v) => `${(v * 1000).toPrecision(3)} mm`;
        const close = new Map();   // wall name → distance of nearest (non-touching, unshielded) conductor
        for (const [name, bc, dist, between, covers] of walls) {
            if (bc !== 'open') continue;
            const shields = this.conductors.filter(g => !g.is_signal && dist(g) <= tol);
            let dMin = Infinity;
            for (const c of this.conductors) {
                const d = dist(c);
                if (d <= tol || d >= dMin) continue;   // ≤ tol: touches this wall — intentional
                if (shields.some(g => between(g, c) && covers(g, c))) continue;
                dMin = d;
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

    // Post-solve check of the open walls from the solved potential. An open wall is a
    // natural (zero normal flux) boundary, exact only where the field has decayed. The
    // field left on it is tangential, and the energy it would carry beyond the wall
    // is about eps * int(E_t^2) ds * d, with d the distance from the wall to the
    // signal conductors (the lateral scale of the wall field). Relative to the stored
    // energy this tracks the truncation error: Z0 reads about twice that fraction
    // high on microstrip, coupled microstrip (both modes) and GCPW, over 0.1..7 %.
    // The wall potential is not usable for this: an air channel between grounded
    // pours holds a constant potential that carries no energy.
    // Vmodes is one potential grid V[iy][ix] per mode on this.x, this.y. Returns a
    // warning object when the estimated Z0 error exceeds OPEN_WALL_TOL, else null.
    openBoundaryFieldWarning(Vmodes) {
        const b = this.boundaries, x = this.x, y = this.y;
        if (!b || !x || !y || !Vmodes || !b.includes('open')) return null;
        const OPEN_WALL_TOL = 0.02;
        const nx = x.length, ny = y.length;
        let sx0 = Infinity, sx1 = -Infinity, sy0 = Infinity, sy1 = -Infinity;
        for (const c of (this.conductors || [])) {
            if (!c.is_signal) continue;
            sx0 = Math.min(sx0, c.x_min); sx1 = Math.max(sx1, c.x_max);
            sy0 = Math.min(sy0, c.y_min); sy1 = Math.max(sy1, c.y_max);
        }
        if (!(sx1 > sx0)) return null;
        const diels = this.dielectrics || [];
        const epsAt = (px, py) => {
            let e = 1;
            for (const d of diels) {
                if (px >= d.x_min && px <= d.x_max && py >= d.y_min && py <= d.y_max) e = d.epsilon_r || 1;
            }
            return e;
        };
        if (!this._openWallCache) this._openWallCache = new WeakMap();
        // On a half-domain grid the first column is the symmetry plane and the right
        // wall stands for both sides.
        const half = x[0] > -this.domain_width / 4;
        const nudge = 1e-9 * Math.max(x[nx - 1] - x[0], y[ny - 1] - y[0]);
        let worst = null;
        for (const V of Vmodes) {
            if (!V || V.length !== ny || !V[0] || V[0].length !== nx) continue;
            let walls = this._openWallCache.get(V);
            if (!walls) {
                let W = 0;
                for (let i = 0; i + 1 < ny; i++) for (let j = 0; j + 1 < nx; j++) {
                    const dx = x[j + 1] - x[j], dy = y[i + 1] - y[i];
                    const ex = 0.5 * ((V[i][j + 1] - V[i][j]) + (V[i + 1][j + 1] - V[i + 1][j])) / dx;
                    const ey = 0.5 * ((V[i + 1][j] - V[i][j]) + (V[i + 1][j + 1] - V[i][j + 1])) / dy;
                    const e2 = ex * ex + ey * ey;
                    if (isFinite(e2)) W += epsAt(0.5 * (x[j] + x[j + 1]), 0.5 * (y[i] + y[i + 1])) * e2 * dx * dy;
                }
                walls = {};
                const column = (j, px, d, name) => {
                    if (!(d > 0) || !(W > 0)) return;
                    let t = 0;
                    for (let i = 0; i + 1 < ny; i++) {
                        const dy = y[i + 1] - y[i], e = (V[i + 1][j] - V[i][j]) / dy;
                        if (isFinite(e)) t += epsAt(px, 0.5 * (y[i] + y[i + 1])) * e * e * dy;
                    }
                    walls[name] = 2 * t * d / W;
                };
                const row = (i, py, d, name) => {
                    if (!(d > 0) || !(W > 0)) return;
                    let t = 0;
                    for (let j = 0; j + 1 < nx; j++) {
                        const dx = x[j + 1] - x[j], e = (V[i][j + 1] - V[i][j]) / dx;
                        if (isFinite(e)) t += epsAt(0.5 * (x[j] + x[j + 1]), py) * e * e * dx;
                    }
                    walls[name] = 2 * t * d / W;
                };
                if (b[0] === 'open' && !half) column(0, x[0] + nudge, sx0 - x[0], 'left');
                if (b[1] === 'open') column(nx - 1, x[nx - 1] - nudge, x[nx - 1] - sx1,
                    half && b[0] === 'open' ? 'sides' : 'right');
                if (b[2] === 'open') row(ny - 1, y[ny - 1] - nudge, y[ny - 1] - sy1, 'top');
                if (b[3] === 'open') row(0, y[0] + nudge, sy0 - y[0], 'bottom');
                this._openWallCache.set(V, walls);
            }
            const total = Object.values(walls).reduce((a, v) => a + v, 0);
            if (!worst || total > worst.total) worst = { total, walls };
        }
        if (!worst || !(worst.total > OPEN_WALL_TOL)) return null;
        // Name the walls that matter; left + right read as "sides".
        const w = { ...worst.walls };
        if (w.left !== undefined && w.right !== undefined) {
            w.sides = w.left + w.right;
            delete w.left; delete w.right;
        }
        const parts = ['sides', 'left', 'right', 'top', 'bottom'].filter(n => w[n] >= 0.2 * worst.total);
        const list = parts.length > 1
            ? parts.slice(0, -1).join(', ') + ' and ' + parts[parts.length - 1]
            : parts[0];
        return { type: 'open-boundary', mode: 'all', estimate: worst.total, walls: parts, message:
            `The field has not decayed at the open boundary on the ${list}: truncating it there is ` +
            `estimated to raise Z0 by about ${(100 * worst.total).toPrecision(2)}%. Enlarge the enclosure ` +
            `or use grounded walls there.` };
    }

    // modesOpts (truthy = guarding a Modes-tab solve, which ALWAYS runs the triangular
    // backend regardless of the Solver dropdown): {wavelengthDensity} — the Modes mesh
    // density, because that solve wavelength-caps the bulk of the WHOLE domain (see
    // TriBackend._wavelengthCap) rather than just the active patch.
    _check_meshability(maxNodes = 20000, freqOverride = undefined, modesOpts = null) {
        const freq = freqOverride ?? this.freq;
        // modesOpts.domainBox: the Modes tab's shrunken open domain (see _modes_domain_box).
        const box = modesOpts && modesOpts.domainBox;
        const W = box ? box.x_max - box.x_min : this.domain_width;
        const yBottom = box ? box.y_min : this.domain_y_min;
        const H = (box ? box.y_max : (this.domain_height ?? NaN)) - yBottom;
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
                const nLambda = modesOpts.wavelengthDensity > 0 ? modesOpts.wavelengthDensity : 8;
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
            // Half-domain symmetry solve meshes only x in [0, W/2].
            const symK = this.sym_half ? 0.5 : 1;
            const nx = axisLines(symK * W, symK * Lx), ny = axisLines(H, Ly);
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
     * @returns {Array<Float64Array>} - The solved voltage array (same reference as input)
     */
    async solve_laplace(V, vacuum = false, planeBC = null) {
        return (await this.solve_laplace_multi([V], vacuum, planeBC))[0];
    }

    // Symmetry-plane boundary condition at x=0 for a half-domain solve: the odd
    // mode sees an electric wall (V=0), everything else the magnetic wall the
    // FDM's natural boundary treatment already provides.
    _plane_bc(mode) {
        return this.sym_half ? (mode === 'odd' ? 'pec' : 'pmc') : null;
    }

    // Solve the same operator (grid + epsilon + conductor mask) for several sets
    // of Dirichlet drive voltages at once. The matrix does not depend on the
    // drive, only the right-hand side does, so all systems share a single
    // assembly and a single factorization (see solveWithWASMMulti). Each Vs[m]
    // is filled in place and the array of solutions is returned.
    //
    // planeBC (half-domain solves only): 'pec' x=0 column set to V=0 (odd
    // mode). Anything else leaves the natural zero-flux treatment (magnetic
    // wall). The matrix depends on it, so drives batched into one call must
    // share the same planeBC.
    async solve_laplace_multi(Vs, vacuum = false, planeBC = null) {
        // Ensure mesh is generated
        if (this.ensure_mesh) {
            this.ensure_mesh();
        }

        for (const V of Vs) {
            const errors = validate_laplace_inputs(
                V, this.x, this.y, this.epsilon_r, this.conductor_mask, vacuum);
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

        const is_cond = (i, j) => this.conductor_mask[i][j];
        // Cell-centred permittivity, the operator's material of record (see
        // FieldSolver2D._paint_cell_materials). The vacuum solve reads 1 everywhere.
        const ec = vacuum ? null : this.epsilon_cell;
        const epsC = ec ? (ci, cj) => ec[ci][cj] : () => 1.0;
        // PEC symmetry plane: the j=0 column is set to V=0 like a conductor.
        const pec = planeBC === 'pec';
        const pinned = (i, j) => is_cond(i, j) || (pec && j === 0);

        // Remove mesh nodes internal to conductors
        // E-field inside conductors is 0.
        const is_unknown = new Int8Array(N);
        let N_unknown = 0;

        for (let i = 0; i < ny; i++)
            for (let j = 0; j < nx; j++) {
                const n = idx(i, j);
                if (!pinned(i, j)) {
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

        // Build sparse system. The 5-point stencil gives at most 5 entries per
        // row, and the natural (row-major) node numbering already emits them in
        // ascending column order (down, left, self, right, up), so the CSR is
        // written straight into typed arrays.
        const Bs = Vs.map(() => new Float64Array(N_unknown));
        const rowPtr = new Int32Array(N_unknown + 1);
        const colIdx = new Int32Array(5 * N_unknown);
        const values = new Float64Array(5 * N_unknown);
        let nnz = 0;
        const addA = (c, v) => {
            colIdx[nnz] = c;
            values[nnz] = v;
            nnz++;
        };

        for (let i = 0; i < ny; i++) {
            for (let j = 0; j < nx; j++) {
                if (pinned(i, j)) continue;

                const fn = idx(i, j);
                const n = full_to_red[fn];

                const [cd, cl, cr, cu] = this._stencil(i, j, dx, dy, epsC, pec);

                const cc = -(cr + cl + cu + cd);

                // Emit the row in ascending column order: (i-1,j), (i,j-1),
                // (i,j), (i,j+1), (i+1,j). Row-major numbering makes that the
                // sorted order, which is what the WASM entry point requires.
                const handle = (ii, jj, c) => {
                    if (!pinned(ii, jj)) {
                        addA(full_to_red[idx(ii, jj)], c);
                    } else {
                        for (let m = 0; m < Vs.length; m++)
                            Bs[m][n] -= c * Vs[m][ii][jj];
                    }
                };

                if (i > 0) handle(i - 1, j, cd);
                if (j > 0) handle(i, j - 1, cl);
                addA(n, cc);
                if (j < nx - 1) handle(i, j + 1, cr);
                if (i < ny - 1) handle(i + 1, j, cu);
                rowPtr[n + 1] = nnz;
            }
        }

        const csr = { rowPtr, colIdx: colIdx.subarray(0, nnz), values: values.subarray(0, nnz) };
        const xs = await solveWithWASMMulti(csr, Bs);

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

    // Flux coefficients [down, left, right, up] of node (i, j) in the Laplace
    // operator of solve_laplace_multi, then the cells and half widths of its control
    // volume. epsC(ci, cj) is the cell permittivity, pec pins the half-domain symmetry
    // plane. The coefficients are negative and symmetric between neighbouring rows
    // (see the notes below).
    _stencil(i, j, dx, dy, epsC, pec) {
        const ny = this.y.length, nx = this.x.length;
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

        // The control volume around node (i,j) spans half a cell in each
        // direction. At a domain edge the missing half is the mirror of
        // the present one (V[-1] = V[1]),
        // which is what the dxl/dyd fallbacks above encode. Each of its
        // four faces is therefore split across two cells, and the flux
        // through it is the cell permittivities weighted by how much of
        // the face each cell covers:
        //
        //     cr = -( eps[cid][j] * hd + eps[ciu][j] * hu ) / dxr
        //
        // and so on. With interfaces on grid lines (which is where the
        // mesher puts them) no cell straddles a material boundary, so
        // this is exact, unlike averaging the two nodal permittivities,
        // which smears every interface by half a cell. It is also
        // symmetric by construction: node (i,j+1)'s left face reads the
        // same cells with the same weights, which is what lets the
        // system be factored by Cholesky rather than LU.
        //
        // Conductors need no special case: a cell fully inside a
        // conductor is only ever read by a node inside that conductor,
        // and those nodes are pinned.
        const cid = i > 0 ? i - 1 : 0;                 // cell row below the node
        const ciu = i < ny - 1 ? i : ny - 2;           // cell row above
        const cjl = j > 0 ? j - 1 : 0;                 // cell column left of the node
        const cjr = j < nx - 1 ? j : nx - 2;           // cell column right
        const hd = 0.5 * dyd, hu = 0.5 * dyu;          // face halves in y
        const wl = 0.5 * dxl, wr = 0.5 * dxr;          // face halves in x

        let cr = 0, cl = 0, cu = 0, cd = 0;
        if (j < nx - 1) cr = -(epsC(cid, j) * hd + epsC(ciu, j) * hu) / dxr;
        if (j > 0)      cl = -(epsC(cid, cjl) * hd + epsC(ciu, cjl) * hu) / dxl;
        if (i < ny - 1) cu = -(epsC(i, cjl) * wl + epsC(i, cjr) * wr) / dyu;
        if (i > 0)      cd = -(epsC(cid, cjl) * wl + epsC(cid, cjr) * wr) / dyd;
        // Magnetic-wall symmetry plane. The mirrored ghost (V[-1] = V[1],
        // dxl = dxr) folds the left flux into the right one, so the
        // full-domain row restricted to x >= 0 is 2*cr*(V1-V0) plus
        // vertical terms over the full width dx[0] (cjl == cjr == 0 there
        // already gives that width).
        if (this.sym_half && j === 0 && !pec) cr *= 2;

        // Half-domain rows at j >= 1 each stand for a mirror pair of
        // full-domain nodes while the plane row (j = 0) stands for one,
        // so the restriction above (cr *= 2) leaves the operator
        // unsymmetric: A[0,1] = 2*cr but A[1,0] = cr. Scaling every
        // j >= 1 row by 2 makes it exactly P^T A P (P = the half -> full
        // prolongation that duplicates x > 0 nodes), which is symmetric
        // positive definite. Same solution, scaling a row and its
        // right-hand side entry together is an identity, but Cholesky
        // can factor it. Cholesky is ~1.4x faster than LU here.
        // (In the PEC case every unknown row has j >= 1, so this is a
        // uniform x2 and the matrix was symmetric either way.)
        if (this.sym_half && j > 0) {
            cr *= 2; cl *= 2; cu *= 2; cd *= 2;
        }
        // Also the node's control volume: the cells around it and the half widths.
        return [cd, cl, cr, cu, cid, ciu, cjl, cjr, hd, hu, wl, wr];
    }

    // The cached _dc_signal_inductance result for a mode's symmetry-plane condition,
    // or null when it has not been computed on the current grid.
    _dc_signal_entry(mode) {
        const c = this._dcSig, key = this._plane_bc(mode) || 'none';
        return c && c.x === this.x && c.y === this.y ? c.m[key] ?? null : null;
    }
    _dc_signal_matrix(mode) { return this._dc_signal_entry(mode)?.M ?? null; }

    async _ensure_dc_signal_inductance(planeBC = null) {
        const key = planeBC || 'none';
        if (!this._dcSig || this._dcSig.x !== this.x || this._dcSig.y !== this.y) this._dcSig = { x: this.x, y: this.y, m: {} };
        const cache = this._dcSig.m;
        if (!(key in cache)) cache[key] = await this._dc_signal_inductance(planeBC).catch(() => null);
        const entry = cache[key];
        // The ground sheet setup only where the return current can spread.
        if (entry && entry.sheet === undefined && this._ground_may_spread()) {
            entry.sheet = await this._ground_sheet_setup(entry.dc, entry.st, entry.nets, entry.pec).catch(() => null);
        }
        return entry;
    }

    // Whether some ground's spreading length delta^2/d can reach the blend range of the
    // thin-sheet model at this frequency (see calculate_conductor_loss), with the return
    // width W_K bounded below by a quarter of the narrowest trace.
    _ground_may_spread() {
        if (!(this.freq > 0)) return true;
        const sig = (this.conductors || []).filter(c => c.is_signal);
        const wLo = 0.25 * Math.min(...sig.map(c => Math.abs(c.width)));
        const omega = 2 * Math.PI * this.freq;
        const ideal = this._unlimited_grounds();
        return (this.conductors || []).some((c, ci) => {
            if (c.is_signal || ideal.has(ci)) return false;
            const d = Math.min(Math.abs(c.width), Math.abs(c.height));
            return 2 / (omega * CONSTANTS.MU0 * this._bulk_sigma(c)) / d / wLo > 0.02;
        });
    }

    // DC internal inductance of the signal conductors with perfectly conducting grounds,
    // per unit current: a 2x2 matrix over the positive (or single-ended) and negative
    // traces seen on this grid, zero where a trace is absent, or null. The traces carry
    // uniform current, split by conductivity, and the result is the loop inductance
    // (integral of J A) less the external inductance of perfectly conducting traces on
    // the same operator, so their discretization errors largely cancel. A half-domain
    // trace stands for itself and its mirror image: the currents are the full-domain
    // conjugates of the drives. Frequency independent, cached per grid and planeBC.
    //
    // Returns { M, ... } with what the thin-sheet ground model (_ground_sheet_setup)
    // needs of the same system.
    async _dc_signal_inductance(planeBC) {
        if (!this.conductors || !this.ground_mask) return null;
        const nx = this.x.length, ny = this.y.length, ncx = nx - 1;
        const dx = diff(this.x), dy = diff(this.y);
        const pec = planeBC === 'pec';
        const mult = this.sym_half ? 2 : 1;
        const netOf = c => (this.is_differential && c.polarity < 0 ? 1 : 0);
        // Signal net and conductivity of each cell, the later conductor winning.
        const cellNet = new Int8Array((ny - 1) * ncx).fill(-1), cellSig = new Float64Array((ny - 1) * ncx);
        for (const c of this.conductors) {
            if (!c.is_signal) continue;
            const sig = this._bulk_sigma(c);
            for (let i = 0; i < ny - 1; i++) {
                const yc = 0.5 * (this.y[i] + this.y[i + 1]);
                if (yc < c.y_min || yc > c.y_max) continue;
                for (let j = 0; j < ncx; j++) {
                    const xc = 0.5 * (this.x[j] + this.x[j + 1]);
                    if (xc < c.x_min || xc > c.x_max) continue;
                    cellNet[i * ncx + j] = netOf(c); cellSig[i * ncx + j] = sig;
                }
            }
        }
        const G = [0, 0];
        for (let i = 0; i < ny - 1; i++) for (let j = 0; j < ncx; j++) {
            const k = i * ncx + j;
            if (cellNet[k] >= 0) G[cellNet[k]] += cellSig[k] * dx[j] * dy[i] * mult;
        }
        const nets = [0, 1].filter(k => G[k] > 0);
        if (!nets.length) return null;
        const isSig = (i, j) => (this.is_differential ? this.signal_p_mask[i][j] || this.signal_n_mask[i][j] : this.signal_mask[i][j]);
        const sigNet = (i, j) => (this.is_differential && this.signal_n_mask[i][j] ? 1 : 0);
        let grounded = pec;
        for (let i = 0; i < ny && !grounded; i++) for (let j = 0; j < nx; j++) if (this.ground_mask[i][j]) { grounded = true; break; }
        if (!grounded) return null;

        const one = () => 1.0;
        const st = [];
        for (let i = 0; i < ny; i++) for (let j = 0; j < nx; j++) st.push(this._stencil(i, j, dx, dy, one, pec));
        // Solves the vacuum operator with pinned nodes at fixed(m, i, j) and the source
        // density J(m, cell) (per unit mu0). Returns the full node arrays and the
        // right-hand sides.
        const solve = async (pinned, fixed, J, nDrives) => {
            const idx = new Int32Array(nx * ny).fill(-1);
            let N = 0;
            for (let n = 0; n < nx * ny; n++) if (!pinned((n / nx) | 0, n % nx)) idx[n] = N++;
            const Bs = Array.from({ length: nDrives }, () => new Float64Array(N));
            const rowPtr = new Int32Array(N + 1), colIdx = new Int32Array(5 * N), values = new Float64Array(5 * N);
            let nnz = 0;
            for (let i = 0; i < ny; i++) for (let j = 0; j < nx; j++) {
                const r = idx[i * nx + j];
                if (r < 0) continue;
                const [cd, cl, cr, cu, cid, ciu, cjl, cjr, hd, hu, wl, wr] = st[i * nx + j];
                const scale = this.sym_half && j > 0 ? 2 : 1;
                for (let m = 0; m < nDrives; m++) {
                    Bs[m][r] += scale * (J(m, cid * ncx + cjl) * wl * hd + J(m, cid * ncx + cjr) * wr * hd
                        + J(m, ciu * ncx + cjl) * wl * hu + J(m, ciu * ncx + cjr) * wr * hu);
                }
                const nb = (ii, jj, c) => {
                    const q = idx[ii * nx + jj];
                    if (q >= 0) { colIdx[nnz] = q; values[nnz++] = c; }
                    else for (let m = 0; m < nDrives; m++) Bs[m][r] -= c * fixed(m, ii, jj);
                };
                if (i > 0) nb(i - 1, j, cd);
                if (j > 0) nb(i, j - 1, cl);
                colIdx[nnz] = r; values[nnz++] = -(cd + cl + cr + cu);
                if (j < nx - 1) nb(i, j + 1, cr);
                if (i < ny - 1) nb(i + 1, j, cu);
                rowPtr[r + 1] = nnz;
            }
            const csr = { rowPtr, colIdx: colIdx.subarray(0, nnz), values: values.subarray(0, nnz) };
            const xs = await solveWithWASMMulti(csr, Bs);
            const As = xs.map((x, m) => {
                const A = new Float64Array(nx * ny);
                for (let n = 0; n < nx * ny; n++) A[n] = idx[n] >= 0 ? x[idx[n]] : fixed(m, (n / nx) | 0, n % nx);
                return A;
            });
            return { As, Bs, idx, csr, N };
        };

        // Uniform current in each net, unit full-domain current per drive.
        const J = (m, k) => (cellNet[k] === nets[m] ? cellSig[k] / G[nets[m]] : 0);
        const dcPinned = (i, j) => this.ground_mask[i][j] || (pec && j === 0);
        const dc = await solve(dcPinned, () => 0, J, nets.length);
        const MU0 = CONSTANTS.MU0;
        const Ldc = nets.map(() => nets.map(() => 0));
        for (let a = 0; a < nets.length; a++) for (let b = 0; b < nets.length; b++) {
            let w = 0;
            for (let n = 0; n < nx * ny; n++) { const r = dc.idx[n]; if (r >= 0) w += dc.Bs[a][r] * dc.As[b][n]; }
            Ldc[a][b] = MU0 * w;
        }

        // Perfectly conducting traces at unit potential, one drive per net: the energy
        // form K of the same operator gives the external inductance mu0 K^-1.
        const pecPinned = (i, j) => dcPinned(i, j) || isSig(i, j);
        const ext = await solve(pecPinned, (m, i, j) => (isSig(i, j) && sigNet(i, j) === nets[m] ? 1 : 0), () => 0, nets.length);
        const K = nets.map(() => nets.map(() => 0));
        for (let i = 0; i < ny; i++) for (let j = 0; j < nx; j++) {
            const n = i * nx + j;
            const edges = [];
            if (j < nx - 1) edges.push([n + 1, -st[n + 1][1]]);   // right edge, the neighbour's left coefficient
            if (i < ny - 1) edges.push([n + nx, -st[n][3]]);      // up edge
            for (const [q, g] of edges) {
                for (let a = 0; a < nets.length; a++) for (let b = 0; b < nets.length; b++) {
                    K[a][b] += g * (ext.As[a][n] - ext.As[a][q]) * (ext.As[b][n] - ext.As[b][q]);
                }
            }
        }
        let Lext;
        if (nets.length === 1) Lext = [[MU0 / K[0][0]]];
        else {
            const det = K[0][0] * K[1][1] - K[0][1] * K[1][0];
            Lext = [[MU0 * K[1][1] / det, -MU0 * K[0][1] / det], [-MU0 * K[1][0] / det, MU0 * K[0][0] / det]];
        }
        const M = [[0, 0], [0, 0]];
        nets.forEach((a, ia) => nets.forEach((b, ib) => { M[a][b] = Ldc[ia][ib] - Lext[ia][ib]; }));
        // The thin-sheet ground setup reuses the DC solve, see _ensure_dc_signal_inductance.
        return { M, dc, st, nets, pec, sheet: undefined };
    }

    // Ground slabs absorbed into the domain walls, and the grounds of unlimited width:
    // the walls and every ground reaching an open domain edge (unlimitedGrounds, the
    // rule of both backends). The unlimited grounds are ideal returns: at DC they carry
    // the return current with no resistance, their surface model spreads it without
    // limit towards DC, and the thin-sheet model below is left to the other grounds,
    // conductors in open space.
    _ground_classes() {
        const key = this.conductors;
        const dom = { x_min: -this.domain_width / 2, x_max: this.domain_width / 2,
            y_min: this.domain_y_min ?? -this.t_gnd, y_max: this.domain_height };
        const c = this._gndClasses;
        if (c && c.key === key && c.x_max === dom.x_max && c.y_min === dom.y_min && c.y_max === dom.y_max) return c;
        const cls = (this.domain_shape || !key) ? { walls: new Set(), unlimited: new Set() }
            : unlimitedGrounds(dom, key, this.boundaries, this.domain_width * 1e-9);
        this._gndClasses = { key, ...dom, ...cls };
        return this._gndClasses;
    }
    _wall_grounds() { return this._ground_classes().walls; }
    _unlimited_grounds() { return this._ground_classes().unlimited; }

    // Thin-sheet model of the grounds for the return current spreading sideways at low
    // frequency. Each ground conductor is tied across its thin dimension into groups
    // (a column of a horizontal ground, a row of a vertical one) that carry the sheet
    // current Y_s (E0 - j omega A) with Y_s the admittance of a one-sided slab. A thick
    // ground (a via slab) is split across its thickness into layers no thicker than
    // half its distance to the nearest trace: towards DC it is transparent to the
    // field, which one potential across its whole thickness would stiffen. The rest
    // of the grid is eliminated on the factorization of the DC trace solve (grounds
    // pinned there too), leaving the dense Schur complement S on the groups and the
    // right-hand side r of each trace drive. Per frequency (S + j omega mu0 D) u = r
    // then gives the ground's resistance and its inductance above perfectly conducting
    // grounds, see _ground_sheet_impedance.
    async _ground_sheet_setup(dc, st, nets, pec) {
        const nx = this.x.length, ny = this.y.length;
        // Unlimited grounds stay perfectly conducting here (their nodes stay pinned at 0).
        const ideal = this._unlimited_grounds();
        const signals = this.conductors.filter(c => c.is_signal);
        const layers = this.conductors.map((c, ci) => {
            if (c.is_signal || ideal.has(ci)) return null;
            const d = Math.min(Math.abs(c.width), Math.abs(c.height));
            const gap = Math.min(...signals.map(sc => bodyDistance(sc, c)));
            const n = gap > 0 ? Math.max(1, Math.ceil(d / (0.5 * gap))) : 1;
            return { n, d: d / n };
        });
        const groupOf = new Int32Array(nx * ny).fill(-1);
        const groups = [], byKey = new Map();
        for (let i = 0; i < ny; i++) for (let j = 0; j < nx; j++) {
            if (!this.ground_mask[i][j] || (pec && j === 0)) continue;
            const ci = this.conductor_id[i][j];
            const c = this.conductors[ci];
            if (!c || c.is_signal || ideal.has(ci)) continue;
            const horiz = Math.abs(c.height) <= Math.abs(c.width);
            const { n, d } = layers[ci];
            const across = horiz ? (this.y[i] - c.y_min) : (this.x[j] - c.x_min);
            const layer = Math.min(n - 1, Math.max(0, Math.floor(across / d)));
            const key = `${ci}:${horiz ? 'c' : 'r'}${horiz ? j : i}:${layer}`;
            let g = byKey.get(key);
            if (g === undefined) {
                g = groups.length;
                byKey.set(key, g);
                groups.push({ ci, d, sigma: this._bulk_sigma(c), area: 0 });
            }
            groupOf[i * nx + j] = g;
            // Conductor area of the node's control volume, with its row weight.
            const [, , , , , , , , hd, hu, wl, wr] = st[i * nx + j];
            const ox = Math.min(this.x[j] + wr, c.x_max) - Math.max(this.x[j] - wl, c.x_min);
            const oy = Math.min(this.y[i] + hu, c.y_max) - Math.max(this.y[i] - hd, c.y_min);
            if (ox > 0 && oy > 0) groups[g].area += (this.sym_half && j > 0 ? 2 : 1) * ox * oy;
        }
        const m = groups.length;
        if (m === 0 || m > 1500) return null;

        // Couplings: C (eliminated node -> group, entries of the eliminated node's row)
        // and the group block P^T K_gg P.
        const S = Array.from({ length: m }, () => new Float64Array(m));
        const cEntries = [];                 // [reduced row of the eliminated node, group, coefficient]
        const nbs = (n, i, j) => {
            const [cd, cl, cr, cu] = st[n];
            const out = [];
            if (i > 0) out.push([n - nx, cd]);
            if (j > 0) out.push([n - 1, cl]);
            if (j < nx - 1) out.push([n + 1, cr]);
            if (i < ny - 1) out.push([n + nx, cu]);
            return out;
        };
        for (let i = 0; i < ny; i++) for (let j = 0; j < nx; j++) {
            const n = i * nx + j, g = groupOf[n];
            if (g >= 0) {
                const [cd, cl, cr, cu] = st[n];
                S[g][g] -= cd + cl + cr + cu;
                for (const [q, c] of nbs(n, i, j)) if (groupOf[q] >= 0) S[g][groupOf[q]] += c;
            } else if (dc.idx[n] >= 0) {
                for (const [q, c] of nbs(n, i, j)) if (groupOf[q] >= 0) cEntries.push([dc.idx[n], groupOf[q], c]);
            }
        }
        // r_b = -C^T z_b with z_b the DC trace solution of drive b.
        const r = nets.map((_, b) => {
            const v = new Float64Array(m);
            const z = dc.As[b];
            for (let i = 0; i < ny; i++) for (let j = 0; j < nx; j++) {
                const n = i * nx + j;
                if (dc.idx[n] < 0 || groupOf[n] >= 0) continue;
                for (const [q, c] of nbs(n, i, j)) if (groupOf[q] >= 0) v[groupOf[q]] -= c * z[n];
            }
            return v;
        });
        // S -= C^T K_aa^-1 C, in chunks of right-hand sides (the solver keeps the factor
        // of the matrix across the chunks), keeping only the rows C touches.
        const colEntries = Array.from({ length: m }, () => []);
        for (const e of cEntries) colEntries[e[1]].push(e);
        const chunk = Math.max(1, Math.min(m, Math.floor(4e6 / dc.N)));
        for (let k0 = 0; k0 < m; k0 += chunk) {
            const ks = [];
            for (let k = k0; k < Math.min(m, k0 + chunk); k++) ks.push(k);
            const Bs = ks.map(k => {
                const B = new Float64Array(dc.N);
                for (const [row, , c] of colEntries[k]) B[row] += c;
                return B;
            });
            const xs = await solveWithWASMMulti(dc.csr, Bs);
            ks.forEach((l, t) => {
                const X = xs[t];
                for (const [row, k, c] of cEntries) S[k][l] -= c * X[row];
            });
        }
        return { S, r, groups, nets };
    }

    // _ground_sheet_impedance of a mode's setup, computed once per frequency.
    _ground_sheet_cached(mode) {
        const sheet = this._dc_signal_entry(mode)?.sheet;
        if (!sheet) return null;
        if (sheet.freq !== this.freq) { sheet.freq = this.freq; sheet.Z = this._ground_sheet_impedance(sheet); }
        return sheet.Z;
    }

    // Ground contribution of the thin-sheet model at the current frequency, as 2x2
    // matrices over the traces per unit full-domain current: { R, L } with L the
    // inductance added to perfectly conducting grounds. null without a setup.
    _ground_sheet_impedance(sheet) {
        const omega = 2 * Math.PI * this.freq, MU0 = CONSTANTS.MU0;
        if (!sheet || !(omega > 0)) return null;
        const { S, r, groups, nets } = sheet;
        const m = groups.length;
        // Sheet admittance of each group: 1/Zs of a slab driven from one side, sigma d
        // for a thin one, times the group's sheet width.
        const Y = groups.map(g => {
            const delta = Math.sqrt(2 / (omega * MU0 * g.sigma));
            const q = new Complex(g.d / delta, g.d / delta);
            // coth(q) = (e^2q + 1) / (e^2q - 1), capped where it is 1 to double precision.
            const coth = q.re > 20 ? new Complex(1, 0) : (() => {
                const e = new Complex(Math.exp(2 * q.re) * Math.cos(2 * q.im), Math.exp(2 * q.re) * Math.sin(2 * q.im));
                return new Complex(e.re + 1, e.im).div(new Complex(e.re - 1, e.im));
            })();
            const Zs = new Complex(1 / (g.sigma * delta), 1 / (g.sigma * delta)).mul(coth);
            return new Complex(1, 0).div(Zs).mul(new Complex(g.area / g.d, 0));
        });
        // Complex symmetric LDL^T of S + j omega mu0 D, no pivoting: its real part is
        // positive semidefinite and the imaginary part a positive diagonal.
        const Are = S.map(row => Float64Array.from(row)), Aim = S.map(() => new Float64Array(m));
        for (let k = 0; k < m; k++) { Are[k][k] -= omega * MU0 * Y[k].im; Aim[k][k] += omega * MU0 * Y[k].re; }
        // Lower triangle only: L[i][k] overwrites A[i][k], D the diagonal.
        const tr = new Float64Array(m), ti = new Float64Array(m);
        for (let k = 0; k < m; k++) {
            const pr = Are[k][k], pi = Aim[k][k], pd = pr * pr + pi * pi;
            for (let i = k + 1; i < m; i++) { tr[i] = Are[i][k]; ti[i] = Aim[i][k]; }
            for (let i = k + 1; i < m; i++) {
                if (tr[i] === 0 && ti[i] === 0) continue;
                const fr = (tr[i] * pr + ti[i] * pi) / pd, fi = (ti[i] * pr - tr[i] * pi) / pd;
                const Ri = Are[i], Ii = Aim[i];
                for (let j = k + 1; j <= i; j++) {
                    Ri[j] -= fr * tr[j] - fi * ti[j];
                    Ii[j] -= fr * ti[j] + fi * tr[j];
                }
                Ri[k] = fr; Ii[k] = fi;
            }
        }
        const solve = b => {
            const ur = Float64Array.from(b), ui = new Float64Array(m);
            for (let i = 0; i < m; i++) for (let k = 0; k < i; k++) {
                const lr = Are[i][k], li = Aim[i][k];
                ur[i] -= lr * ur[k] - li * ui[k]; ui[i] -= lr * ui[k] + li * ur[k];
            }
            for (let i = 0; i < m; i++) {
                const dr = Are[i][i], di = Aim[i][i], dd = dr * dr + di * di;
                const xr = (ur[i] * dr + ui[i] * di) / dd, xi = (ui[i] * dr - ur[i] * di) / dd;
                ur[i] = xr; ui[i] = xi;
            }
            for (let i = m - 1; i >= 0; i--) for (let k = i + 1; k < m; k++) {
                const lr = Are[k][i], li = Aim[k][i];
                ur[i] -= lr * ur[k] - li * ui[k]; ui[i] -= lr * ui[k] + li * ur[k];
            }
            return [ur, ui];
        };
        const U = r.map(solve);
        // Per pair of drives: the magnetic energy change outside the grounds,
        // mu0 Re[-conj(u_a) . (-r_b) - j omega mu0 sum D conj(u_a) u_b] (the DC energy
        // of perfectly conducting grounds cancels), and the complex power into the
        // sheets, omega^2 mu0^2 sum conj(Y) conj(u_a) u_b.
        const R = [[0, 0], [0, 0]], L = [[0, 0], [0, 0]];
        nets.forEach((na, a) => nets.forEach((nb, b) => {
            const [ar, ai] = U[a], [br, bi] = U[b];
            let eRe = 0, pRe = 0, pIm = 0;
            for (let k = 0; k < m; k++) {
                // conj(u_a) u_b
                const cr = ar[k] * br[k] + ai[k] * bi[k], cim = ar[k] * bi[k] - ai[k] * br[k];
                eRe += ar[k] * r[b][k];                               // Re(conj(u_a) r_b)
                // Re(-j omega mu0 D conj(u_a) u_b) = omega mu0 Im(Y conj(u_a) u_b), D = Y_k
                eRe += omega * MU0 * (Y[k].re * cim + Y[k].im * cr);
                // conj(Y) conj(u_a) u_b
                pRe += Y[k].re * cr + Y[k].im * cim;
                pIm += Y[k].re * cim - Y[k].im * cr;
            }
            const w2 = omega * omega * MU0 * MU0;
            R[na][nb] = w2 * pRe;
            L[na][nb] = MU0 * eRe + w2 * pIm / omega;
        }));
        // Symmetric parts.
        for (const X of [R, L]) { const o = 0.5 * (X[0][1] + X[1][0]); X[0][1] = X[1][0] = o; }
        return { R, L };
    }

    /**
     * Compute E-field from voltage distribution.
     * @param {Array<Float64Array>} V - 2D voltage array
     * @param {string|null} planeBC - symmetry-plane BC at x=0 ('pec'|'pmc'|null).
     * @returns {{Ex: Array<Float64Array>, Ey: Array<Float64Array>}} - E-field components
     */
    compute_fields(V, planeBC = null) {
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
        // Symmetry plane at x=0: the loss/dielectric integrands read the j=0
        // column, so fill it with the parity-exact values instead of zeros.
        // PMC (even/single): V mirrors evenly -> Ex=0, Ey from the same 3-point
        // y-formula. PEC (odd): V mirrors oddly (V[-1]=-V[1], V[0]=0) -> Ey=0,
        // and the central x-difference collapses to -V[1]/dx[0].
        if (this.sym_half && planeBC) {
            for (let i = 1; i < ny - 1; i++) {
                if (this.conductor_mask[i][0]) continue;
                if (planeBC === 'pec') {
                    Ex[i][0] = -V[i][1] / dx[0];
                } else {
                    const dyd = dy[i-1];
                    const dyu = dy[i];
                    Ey[i][0] = -(
                        (dyd / (dyu * (dyd + dyu))) * V[i+1][0] +
                        ((dyu - dyd) / (dyd * dyu)) * V[i][0] -
                        (dyu / (dyd * (dyd + dyu))) * V[i-1][0]
                    );
                }
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

        // Iterate over signal trace interface. On a half-domain (sym_half) solve
        // the plane column j=0 is included: a mirrored signal has
        // real top/bottom flux there (area = the half cell dx[0]/2), and no flux
        // crosses the plane itself.
        const j0 = this.sym_half ? 0 : 1;
        for (let i = 1; i < ny - 1; i++) {
            for (let j = j0; j < nx - 1; j++) {
                if (!this.signal_mask[i][j]) continue;

                // Half-widths of the node's control volume, and the cells its
                // faces cut through. The same decomposition the operator uses
                // (see solve_laplace_multi). Reading the face permittivity off
                // the cells instead of the neighbour node makes this contour the
                // discrete Gauss law of the system that was actually solved, so
                // Q is exact for the computed potential rather than an
                // independent quadrature that disagrees at every interface.
                const hd = get_dy(i - 1) / 2, hu = get_dy(i) / 2;
                const wl = j > 0 ? get_dx(j - 1) / 2 : get_dx(0) / 2;
                const wr = get_dx(j) / 2;
                const cjl = j > 0 ? j - 1 : 0;

                // Check 4 neighbors
                const check_neighbor = (ni, nj, is_vertical_flux) => {
                    // Only add flux if the neighbor is NOT part of the signal conductor
                    if (this.signal_mask[ni][nj]) return;

                    // E-field Normal
                    let En;
                    let dist;
                    let area;
                    let er = 1;

                    if (is_vertical_flux) {
                         // Neighbor is Top/Bottom
                         dist = Math.abs(this.y[i] - this.y[ni]);
                         En = (V[i][j] - V[ni][nj]) / dist;
                         // Average dx for area (half cell at the symmetry plane)
                         area = j > 0 ? (get_dx(j-1) + get_dx(j)) / 2 : get_dx(0) / 2;
                         if (!vacuum) {
                             const row = this.epsilon_cell[ni > i ? i : i - 1];
                             er = (row[cjl] * wl + row[j] * wr) / (wl + wr);
                         }
                    } else {
                        // Neighbor is Left/Right
                        dist = Math.abs(this.x[j] - this.x[nj]);
                        En = (V[i][j] - V[ni][nj]) / dist;
                        // Average dy for area
                        area = (get_dy(i-1) + get_dy(i)) / 2;
                        if (!vacuum) {
                            const col = nj > j ? j : j - 1;
                            er = (this.epsilon_cell[i - 1][col] * hd +
                                  this.epsilon_cell[i][col] * hu) / (hd + hu);
                        }
                    }

                    Q += CONSTANTS.EPS0 * er * En * area;
                };

                // Right neighbor
                if (!this.signal_mask[i][j + 1]) {
                    check_neighbor(i, j + 1, false);
                }
                // Left neighbor
                if (j > 0 && !this.signal_mask[i][j - 1]) {
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
        // A signal cut in half by the symmetry plane captures exactly half the
        // charge. A signal entirely inside x > 0 (one trace of a differential
        // pair) keeps its full contour.
        const scale = (this.sym_half && this._sym_signal_straddles) ? 2 : 1;
        return scale * Math.abs(Q);
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
        // Overlapping conductors of one kind count the shared area once.
        const visible = visibleAreas(this.conductors);

        for (const [i, cond] of this.conductors.entries()) {
            const area = visible && visible.has(i) ? visible.get(i) : Math.abs(cond.width * cond.height);
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
     * Signal and ground are conductors in series: R_total = R_signal + R_ground.
     * - R_signal = sqrt(R_dc^2 + R_ac^2) of the signal conductors (DC resistance of
     *   the cross-section, skin-effect surface integral with roughness)
     * - R_ground = the surface integral on the ground conductors with the slab
     *   resistance of their finite thickness and the lateral spreading of the return
     *   current, floored at the geometric DC resistance
     * R_ac and R_dc are returned as the sums over both.
     *
     * Two integrand variants:
     *
     * vacuum_fields = true (production for rect-based solvers): Ex/Ey are the
     * vacuum (C0) solve fields and Z0 is Z0_vac = 1/(c*C0). In the quasi-TEM
     * skin-effect limit the H pattern is the harmonic conjugate of the vacuum
     * potential (the same identity as L_ext = 1/(c^2*C0)), so H_t = E_n(vac)/η0,
     * the surface current distribution, which is permittivity independent.
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
    // line (0 = positive trace, 1 = negative trace, 2 = both): Ex/Ey are the fields of a
    // pair with unit current in that trace only, or in both (see _line_asymmetry). The
    // result is then the quadratic form I^T X I of R and of the internal L for those
    // currents: no differential power factor, the DC resistance and skin depth of the
    // driven trace(s).
    calculate_conductor_loss(Ex, Ey, Z0, vacuum_fields = false, mode = null, line = null) {
        if (!this.solution_valid) throw new Error("Fields invalid");

        const { signal_area, ground_area } = this._calculate_conductor_area();

        // DC resistance per unit length, in the same per-line convention as the AC
        // integral (power_factor below). For a differential pair signal_area sums
        // both traces, and the two modes have different DC return paths:
        //   odd  - current returns through the partner trace, the ground carries no
        //          net current at DC: R_dc = R_one_trace = 2/(σ*signal_area).
        //   even - both traces carry I, the ground returns 2I:
        //          P = 2I^2*R_even = 2I^2*R_one_trace + (2I)^2*R_gnd
        //          => R_dc = R_one_trace + 2*R_gnd.
        // The mode-blind form 1/(σ*signal_area) + 1/(σ*ground_area) is ~2x low
        // per line.
        // Signal metal is the plating metal when the plating is at least as thick
        // as the trace (the whole cross-section is plating).
        const sigma_sig = this._signal_sigma(line === 2 ? null : line);
        const ownSigma = this._own_sigma();
        // Per-conductor conductivities: conductances add within the positive traces,
        // the negative traces and the grounds.
        // Per-conductor conductances also when a plating layer conducts beside the bulk.
        const thinPlated = (this.conductors || []).some(c => platingArea(c) > 0 && !this._solid_plating(c));
        const dcG = (ownSigma || thinPlated) ? this._dc_conductances() : null;
        // A ground of unlimited width is an ideal return: with one, the return current
        // has no DC resistance and the other grounds carry none of it.
        const R_gnd0 = this._unlimited_grounds().size ? 0
            : dcG ? (dcG.gnd > 0 ? 1.0 / dcG.gnd : 0)
            : ground_area > 0 ? 1.0 / (this.sigma_cond * ground_area) : 0;
        // Signal and ground are separate conductors in series, so each keeps its own
        // DC term and the two resistances add.
        let R_dc_sig, R_dc_gnd;
        if (line !== null) {
            const g = this._dc_conductances();
            // Both traces driven: the two trace resistances add and the ground returns 2 I.
            R_dc_sig = line === 2 ? 1.0 / g.pos + 1.0 / g.neg : 1.0 / (line === 0 ? g.pos : g.neg);
            R_dc_gnd = line === 2 ? 4.0 * R_gnd0 : R_gnd0;
        } else if (this.is_differential && (mode === 'odd' || mode === 'even')) {
            // Both modes put the same current magnitude through each trace, so the
            // per-line value is the mean of the two trace resistances.
            R_dc_sig = dcG ? 0.5 * (1.0 / dcG.pos + 1.0 / dcG.neg) : 2.0 / (sigma_sig * signal_area);
            R_dc_gnd = mode === 'odd' ? 0 : 2.0 * R_gnd0;
        } else {
            R_dc_sig = dcG ? 1.0 / (dcG.pos + dcG.neg) : 1.0 / (sigma_sig * signal_area);
            R_dc_gnd = R_gnd0;
        }
        const R_dc = R_dc_sig + R_dc_gnd;

        if (this.freq === 0) {
            // DC point: R is the geometric DC resistance and the internal inductance
            // is the low-frequency plateau of the surface integral below (slab
            // reactance mu0 d/3 per face once delta >> d). Returning before the
            // skin-transition block, so its warning is cleared explicitly.
            const L_internal = this._dc_internal_inductance(Ex, Ey, Z0, vacuum_fields, mode, line);
            this._skinTransitionWarn = null;
            this._platingTransitionWarn = null;
            return { R_ac: 0, R_dc, R_total: R_dc, L_internal };
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
        // xStart overrides the segment's x origin (half-domain solves evaluate the
        // mirror image of a node's left segment [x[j-1], x[j]]).
        // Bare-metal surface impedance of conductor ci: its own conductivity and
        // roughness when it has them (custom geometry), the solver-wide ones otherwise.
        const sigmaOf = ci => {
            const c = this.conductors && ci >= 0 ? this.conductors[ci] : null;
            return (c && c.sigma > 0) ? c.sigma : this.sigma_cond;
        };
        const bareRq = ci => {
            const c = this.conductors && ci >= 0 ? this.conductors[ci] : null;
            return (c && c.rq !== undefined && c.rq !== null) ? c.rq : rq;
        };
        // A conductor drawn as touching blocks of different metals: a block thinner than
        // a few skin depths does not hide the block behind it, so its outer face takes
        // the layered impedance (the block as a plating on the metal behind). Returns
        // the conductor touching the face opposite to `direction`, or null.
        const backingOf = (ci, direction) => {
            const c = this.conductors[ci];
            const tol = 1e-9 * this.domain_width;
            const covers = (a0, a1, b0, b1) => Math.min(a1, b1) - Math.max(a0, b0) > 0.5 * (a1 - a0);
            for (let k = 0; k < this.conductors.length; k++) {
                const o = this.conductors[k];
                if (k === ci || o.is_signal !== c.is_signal || (o.polarity || 0) !== (c.polarity || 0)) continue;
                const behind = direction === 'd' ? Math.abs(o.y_max - c.y_min) < tol
                    : direction === 'u' ? Math.abs(o.y_min - c.y_max) < tol
                    : direction === 'r' ? Math.abs(o.x_min - c.x_max) < tol
                    : Math.abs(o.x_max - c.x_min) < tol;
                if (!behind) continue;
                const horiz = direction === 'd' || direction === 'u';
                if (horiz ? covers(c.x_min, c.x_max, o.x_min, o.x_max) : covers(c.y_min, c.y_max, o.y_min, o.y_max)) return k;
            }
            return -1;
        };
        const bareZ = (ci, direction = null, r = bareRq(ci)) => {
            const sg = sigmaOf(ci);
            if (ownSigma && direction && !this.conductors[ci].plating) {
                const key = `${ci}_bare_${direction}_${r}`;
                if (!Z_cache.has(key)) {
                    const k = backingOf(ci, direction);
                    const c = this.conductors[ci];
                    const d = (direction === 'd' || direction === 'u') ? Math.abs(c.height) : Math.abs(c.width);
                    let z = (k >= 0 && sigmaOf(k) !== sg)
                        ? calculate_Zrough_layered(this.freq, sigmaOf(k), r, sg, d) : null;
                    if (!z) {
                        // End face of a thin layer (a face shorter than two skin depths,
                        // with the other metal along one of its ends): the current there
                        // spreads into that metal, blend towards its impedance.
                        const ends = (direction === 'd' || direction === 'u') ? ['r', 'l'] : ['d', 'u'];
                        const len = (direction === 'd' || direction === 'u') ? Math.abs(c.width) : Math.abs(c.height);
                        const dlt = Math.sqrt(2 / (2 * Math.PI * this.freq * 4e-7 * Math.PI * sg));
                        const w = Math.min(1, len / (2 * dlt));
                        const ke = ends.map(e => backingOf(ci, e)).find(v => v >= 0 && sigmaOf(v) !== sg);
                        if (w < 1 && ke !== undefined) {
                            const zb = calculate_Zrough(this.freq, sg, r), zo = calculate_Zrough(this.freq, sigmaOf(ke), r);
                            z = new Complex(w * zb.re + (1 - w) * zo.re, w * zb.im + (1 - w) * zo.im);
                        }
                    }
                    Z_cache.set(key, z);
                }
                const z = Z_cache.get(key);
                if (z) return z;
            }
            if (r === rq && sg === this.sigma_cond) return Z_surf_default;
            const key = `${ci}_bare_${r}`;
            if (!Z_cache.has(key)) Z_cache.set(key, calculate_Zrough(this.freq, sg, r));
            return Z_cache.get(key);
        };

        // Plating over the bulk metal. With `dc` set it also stands in for the faces
        // modelled as plating metal alone (top plating down the sides, thick corners):
        // their semi-infinite plating reactance never returns to the bulk value, so
        // its excess over the bulk would grow as 1/sqrt(f) in the low-frequency blend.
        const layeredZ = (ci, platingRq) => {
            const cond = this.conductors[ci];
            const key = `${ci}_layered_${platingRq}`;
            if (!Z_cache.has(key)) Z_cache.set(key, calculate_Zrough_layered(
                this.freq, sigmaOf(ci), platingRq, cond.plating.sigma, cond.plating.thickness));
            return Z_cache.get(key);
        };
        const getZsurf = (ci, direction, i, j, dl, xStart = null, dc = false) => {
            if (!this.conductors || ci < 0) return Z_surf_default;
            const Z_bare = bareZ(ci, direction);
            const cond = this.conductors[ci];
            if (!cond || !cond.plating) return Z_bare;

            // Plating at least as thick as the conductor: the whole cross-section is
            // plating metal and every face sees the solid plating impedance.
            if (this._solid_plating(cond)) {
                const key = `${ci}_solid`;
                if (!Z_cache.has(key)) {
                    Z_cache.set(key, calculate_Zrough(this.freq, cond.plating.sigma, cond.plating.rq ?? 0));
                }
                return Z_cache.get(key);
            }

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
                    if (dc) {
                        Z_plating = layeredZ(ci, cond.plating.rq);
                    } else if (Z_cache.has(key_top_side)) {
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
                        fraction * Z_plating.re + (1 - fraction) * Z_bare.re,
                        fraction * Z_plating.im + (1 - fraction) * Z_bare.im
                    );
                }
            }

            // Geometric coverage of bottom surface by thick side plating
            // Uses fractional coverage for smooth parameter sweeps
            if (direction === 'u' && cond.plating.sides && !cond.plating.bottom && cond.plating.thick_corners) {
                const t = cond.plating.thickness;
                // Cell spans [x[j], x[j] + dl] in x-direction
                const x_start = xStart ?? this.x[j];
                const x_end = x_start + dl;
                // Left side plating covers [x_min, x_min + t]
                const left_overlap = Math.max(0, Math.min(cond.x_min + t, x_end) - Math.max(cond.x_min, x_start));
                // Right side plating covers [x_max - t, x_max]
                const right_overlap = Math.max(0, Math.min(cond.x_max, x_end) - Math.max(cond.x_max - t, x_start));
                const fraction = dl > 0 ? Math.min((left_overlap + right_overlap) / dl, 1.0) : 0;

                if (fraction > 0) {
                    const key_corner = `${ci}_corner_plating`;
                    let Z_plating;
                    if (dc) {
                        Z_plating = layeredZ(ci, bareRq(ci));
                    } else if (Z_cache.has(key_corner)) {
                        Z_plating = Z_cache.get(key_corner);
                    } else {
                        // Side plating material with bulk surface roughness
                        Z_plating = calculate_Zrough(
                            this.freq, cond.plating.sigma, bareRq(ci)
                        );
                        Z_cache.set(key_corner, Z_plating);
                    }

                    if (fraction >= 1.0) return Z_plating;

                    // Weighted average: covered part uses plating, rest uses bulk
                    return new Complex(
                        fraction * Z_plating.re + (1 - fraction) * Z_bare.re,
                        fraction * Z_plating.im + (1 - fraction) * Z_bare.im
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
                if (dc) {
                    Z_corner = layeredZ(ci, bareRq(ci));
                } else if (Z_cache.has(key_corner)) {
                    Z_corner = Z_cache.get(key_corner);
                } else {
                    // At corners: single-layer with plating sigma and bulk rq
                    // - sigma: plating material (extends from sides)
                    // - rq: bulk surface roughness (bottom surface preparation)
                    Z_corner = calculate_Zrough(
                        this.freq, cond.plating.sigma, bareRq(ci)  // Use bulk rq, not plating.rq
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
                const Z_bottom = cond.plating.bottom ? layeredZ(ci, cond.plating.rq) : Z_bare;

                // Weighted average: corner region uses corner plating impedance, bulk uses bottom impedance
                const Z_avg_re = corner_fraction * Z_corner.re + (1 - corner_fraction) * Z_bottom.re;
                const Z_avg_im = corner_fraction * Z_corner.im + (1 - corner_fraction) * Z_bottom.im;
                return new Complex(Z_avg_re, Z_avg_im);
            }

            // Standard surface impedance (no corner effects)
            if (!cond.plating[surface]) return Z_bare;
            return layeredZ(ci, cond.plating.rq);
        };

        const ny = this.y.length;
        const nx = this.x.length;
        const dx_array = diff(this.x);
        const dy_array = diff(this.y);

        const get_dx = j => (j >= 0 && j < dx_array.length) ? dx_array[j] : dx_array[dx_array.length - 1];
        const get_dy = i => (i >= 0 && i < dy_array.length) ? dy_array[i] : dy_array[dy_array.length - 1];

        let sum_H2_dl_R = 0.0; // Sum for Resistance
        let sum_H2_dl_L = 0.0; // Sum for Inductance

        // Surface impedances above assume a semi-infinite metal, whose reactance
        // Im(Zs) = Rs grows as sqrt(f) and makes Im(Zs)/omega diverge as f -> 0.
        // A conductor of finite thickness d has Zs = Rs (1+j) coth((1+j) d/delta),
        // whose reactance tends to omega mu0 d/3 (a finite internal inductance)
        // once delta > d. Only the reactance takes the slab factor: the resistance
        // path has its own DC limit (R_total below) and a transition calibration
        // fitted against the semi-infinite R_ac. A signal trace carries current on
        // both faces: with surface fields a (bottom) and b (top) its DC internal
        // inductance is mu t (a^2 - ab + b^2) / 3 per unit width, the per-face slab
        // of thickness t (1 - ab / (a^2 + b^2)). That is t/2 when both faces carry
        // the same field (stripline) and t when one does (microstrip). Ground planes
        // and pours are treated as one-sided slabs of their full thickness.
        const deltaOf = (sigma) => Math.sqrt(2 / (2 * Math.PI * this.freq * 4e-7 * Math.PI * sigma));
        // delta is the skin depth of the signal metal (transition calibration and
        // warning below), deltaCond the one of each conductor.
        const delta = deltaOf((ownSigma || line !== null) ? sigma_sig : this.sigma_cond);
        const deltaCond = c => deltaOf(this._solid_plating(c) ? c.plating.sigma : (c.sigma > 0 ? c.sigma : this.sigma_cond));
        const slabReactanceFactor = (d, dlt = delta) => {
            const x = d / dlt;
            if (!(x > 0)) return 1;
            if (x > 20) return 1;
            const den = Math.cosh(2 * x) - Math.cos(2 * x);
            return (Math.sinh(2 * x) - Math.sin(2 * x)) / den;
        };
        // A block stacked on another metal of the same conductor is as thick as the stack.
        const stackH = (c, ci) => {
            let h = Math.abs(c.height);
            if (ownSigma) for (const dir of ['d', 'u']) {
                const k = backingOf(ci, dir);
                if (k >= 0) h += Math.abs(this.conductors[k].height);
            }
            return h;
        };
        // Reactance integrand Im(Zs)|H|^2 dl per conductor (the slab factor is applied
        // after the loop, once the face split is known) and |H|^2 dl on the bottom and
        // top faces of each conductor.
        // For the DC blend of signal conductors: the current (the signed tangential
        // H, dl) of each signal net, the reactance integrand of the smooth metal the
        // DC solve describes (sumLref) and the one with the low-frequency face
        // impedances (sumLdc, see layeredZ).
        const nCond = (this.conductors || []).length;
        const sumL = new Float64Array(nCond), faceBot = new Float64Array(nCond), faceTop = new Float64Array(nCond);
        const sumLref = new Float64Array(nCond), sumLdc = new Float64Array(nCond), netI = [0, 0];
        let sumLDefault = 0;
        // Smooth reference reactance of a signal face: the bulk metal under any
        // plating, the smooth face impedance for blocks of different metals (which
        // the DC solve resolves).
        const refIm = (ci, direction) => {
            const c = this.conductors[ci];
            return c.plating ? 1 / (this._bulk_sigma(c) * deltaCond(c)) : bareZ(ci, direction, 0).im;
        };
        // direction points from the dielectric node to the conductor: 'u' is a bottom face.
        // span is the segment getZsurf sees for coverage fractions, dl its quadrature weight.
        const addFace = (ci, direction, i, j, span, xStart, dl, H_tan, H_out) => {
            const Zs = getZsurf(ci, direction, i, j, span, xStart);
            const H2dl = H_tan * H_tan * dl;
            if (isGroundCond(ci)) addGnd(ci, Zs.re, H_tan, dl);
            else sum_H2_dl_R += Zs.re * H2dl;
            if (!(ci >= 0 && ci < nCond)) { sumLDefault += Zs.im * H2dl; return; }
            sumL[ci] += Zs.im * H2dl;
            const c = this.conductors[ci];
            if (c.is_signal && vacuum_fields) {
                netI[this.is_differential && c.polarity < 0 ? 1 : 0] += H_out * dl;
                sumLref[ci] += refIm(ci, direction) * H2dl;
                sumLdc[ci] += getZsurf(ci, direction, i, j, span, xStart, true).im * H2dl;
            }
            if (direction === 'u') faceBot[ci] += H2dl;
            else if (direction === 'd') faceTop[ci] += H2dl;
        };
        // Effective slab thickness of a signal conductor from the split of the surface
        // field between its bottom and top faces, taken over the stack it belongs to.
        const signalSlab = (c, ci) => {
            const members = [ci];
            if (ownSigma) for (const dir of ['d', 'u']) { const k = backingOf(ci, dir); if (k >= 0) members.push(k); }
            let A = 0, B = 0;
            for (const k of members) { A += faceBot[k]; B += faceTop[k]; }
            const share = A + B > 0 ? Math.sqrt(A * B) / (A + B) : 0.5;
            return stackH(c, ci) * (1 - share);
        };
        // Ground conductors take the resistance twin of that factor,
        // Re[(1+j) coth((1+j) d/delta)]: 1 for a thick ground, delta/d (the sheet
        // resistance 1/(sigma d)) once delta passes its thickness. The vacuum-field
        // |H|^2 it weights is the return current confined under the traces, which holds
        // far below the frequency where delta reaches the ground thickness, while the
        // geometric DC resistance spreads the return over the whole ground width.
        // Signal conductors keep the semi-infinite Rs: their DC limit and transition
        // are handled by R_total below.
        const slabResistanceFactor = (d, dlt = delta) => {
            const x = d / dlt;
            if (!(x > 0) || x > 20) return 1;
            return (Math.sinh(2 * x) + Math.sin(2 * x)) / (Math.cosh(2 * x) - Math.cos(2 * x));
        };
        const kR = (this.conductors || []).map((c, ci) => (c.is_signal || !vacuum_fields) ? 1 : slabResistanceFactor(
            Math.min(Math.abs(c.width), stackH(c, ci)), deltaCond(c)));
        const isGroundCond = ci => ci >= 0 && ci < kR.length && !this.conductors[ci].is_signal;
        // Per ground conductor: Re(Zs)-weighted |H|^2 and the plain |H|, |H|^2 moments
        // that give the width of its return current.
        const gndR = kR.map(() => 0), gndS1 = kR.map(() => 0), gndS2 = kR.map(() => 0);
        const addGnd = (ci, zre, H, dl) => {
            gndR[ci] += zre * kR[ci] * H * H * dl; gndS1[ci] += H * dl; gndS2[ci] += H * H * dl;
        };

        const isSignal = (i, j) => this.signal_mask[i][j];
        const isGround = (i, j) => this.ground_mask[i][j];
        const isConductor = (i, j) => isSignal(i,j) || isGround(i,j);

        // Half-domain solve: include the plane column j=0 so the first
        // top/bottom face segment next to the plane is not dropped. The nj
        // bounds check below keeps j=0 from probing a left neighbor, and the
        // plane itself is not a conductor face (PEC pinning is not in
        // conductor_mask), so no spurious cut faces can enter the integral.
        const jlo = this.sym_half ? 0 : 1;
        for (let i = 1; i < ny - 1; i++) {
            for (let j = jlo; j < nx - 1; j++) {
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
                        // Signed: the outward normal field, whose sign is the current's.
                        const E_out = direction === 'r' ? -Ex_val : direction === 'l' ? Ex_val
                            : direction === 'u' ? -Ey_val : Ey_val;
                        const H_out = E_out * eps_fac / Z0_freespace;

                        // Look up per-surface impedance (with plating if applicable)
                        const ci = this.conductor_id ? this.conductor_id[ni][nj] : -1;

                        if (this.sym_half && (direction === 'u' || direction === 'd')) {
                            // Half-domain solve, horizontal faces: the full-domain
                            // quadrature weights each node by the segment to its
                            // right. Mirroring flips that orientation, so on a
                            // symmetric grid a node collects its own right segment
                            // and the mirror image of its left segment. Each at
                            // half weight (the global x2 compensates it),
                            // with the true segment span passed to getZsurf so
                            // plating-coverage fractions see the correct geometry.
                            // The plane node (j=0) has only its right segment,
                            // the mirror of [-x1, 0] belongs to node 1's left term.
                            const segs = j > 0
                                ? [[j, get_dx(j)], [j - 1, get_dx(j - 1)]]
                                : [[0, get_dx(0)]];
                            for (const [js, dseg] of segs)
                                addFace(ci, direction, i, j, dseg, this.x[js], dseg / 2, H_tan, H_out);
                        } else if (this.centred_loss_quadrature && !this.sym_half) {
                            // Full-domain solve with different surface finishes: the one-sided
                            // rule below gives mirrored conductors unequal shares of the
                            // integral (they only add up right as a pair), so each node
                            // takes half of the segment on either side instead. The total
                            // of a mirror-symmetric geometry is the same either way, and on a
                            // half domain mirrored conductors share their finish by construction.
                            const horiz = direction === 'u' || direction === 'd';
                            const k = horiz ? j : i, get = horiz ? get_dx : get_dy;
                            const segs = k > 0 ? [[k, get(k)], [k - 1, get(k - 1)]] : [[0, get(0)]];
                            for (const [ks, dseg] of segs)
                                addFace(ci, direction, i, j, dseg, horiz ? this.x[ks] : null, dseg / 2, H_tan, H_out);
                        } else {
                            const dl = dl_func(dl_idx);
                            addFace(ci, direction, i, j, dl, null, dl, H_tan, H_out);
                        }
                    }
                }
            }
        }

        // Signal conductors with the DC internal inductance of uniform trace current
        // (_dc_signal_inductance) blend it with the semi-infinite surface value of the
        // smooth metal, 1/L^2 = 1/L_dc^2 + 1/L_ac^2, the counterpart of
        // R = sqrt(R_dc^2 + R_ac^2): the current spreads over the cross-section once the
        // surface reactance passes the DC value, before delta reaches the thickness. What
        // roughness and plating add to the surface reactance is a surface layer and
        // stays on top. Without the DC value each signal face takes the slab factor of
        // its effective thickness (signalSlab).
        const omega = 2 * Math.PI * this.freq;
        const dcM = vacuum_fields && omega > 0 ? this._dc_signal_matrix(line === null ? mode : null) : null;
        let Ldc = 0;
        if (dcM) {
            // Full-domain currents: a half domain holds half of each.
            const m = this.sym_half ? 2 : 1, I0 = m * netI[0], I1 = m * netI[1];
            Ldc = (I0 * I0 * dcM[0][0] + 2 * I0 * I1 * dcM[0][1] + I1 * I1 * dcM[1][1]) / m;
        }
        const idealGnd = this._unlimited_grounds();
        // Lateral spreading of the return current in a ground thinner than delta
        // (wallSpreadFactor): the current of effective width W_K = (int|K|)^2 / int|K|^2
        // diffuses sideways over delta^2/d. A ground conductor cannot spread wider than
        // itself; a ground of unlimited width has no such limit, so its resistance
        // vanishes towards DC, where it is an ideal return. Its reactance keeps the slab
        // value of the confined current (the DC convention of dcLineParameters): the
        // inductance the spreading adds outside the ground is not modelled, and removing
        // the one inside would make L rise with frequency.
        const spreadG = kR.map(() => 1), spreadU = kR.map(() => 0);
        for (let c = 0; c < gndR.length; c++) {
            if (!(gndS2[c] > 0) || !vacuum_fields) continue;
            const cond = this.conductors[c];
            const d = Math.min(Math.abs(cond.width), Math.abs(cond.height));
            const wMax = Math.max(Math.abs(cond.width), Math.abs(cond.height));
            // Half-domain moments cover x >= 0: a conductor straddling the plane has
            // twice the width seen, one beside it is complete (its mirror image is a
            // separate conductor).
            const x0 = Math.min(cond.x, cond.x + cond.width), x1 = Math.max(cond.x, cond.x + cond.width);
            const straddles = this.sym_half && x0 < 0 && x1 > 0;
            const Wk = (straddles ? 2 : 1) * gndS1[c] * gndS1[c] / gndS2[c];
            const dlt = deltaCond(cond);
            spreadU[c] = dlt * dlt / d / Wk;
            const g = wallSpreadFactor(2 * Math.PI * spreadU[c]);
            spreadG[c] = idealGnd.has(c) ? g : Math.max(g, Math.min(1, Wk / wMax));
        }
        let sumSig = 0, sumExcess = 0, sumLgnd = 0, sumLsheet = 0;
        (this.conductors || []).forEach((c, ci) => {
            if (!(sumL[ci] !== 0)) return;
            if (c.is_signal && Ldc > 0) {
                sumSig += sumLref[ci]; sumExcess += sumLdc[ci] - sumLref[ci];
                return;
            }
            if (!c.is_signal) {
                const v = sumL[ci] * slabReactanceFactor(Math.min(Math.abs(c.width), stackH(c, ci)), deltaCond(c));
                if (idealGnd.has(ci)) sumLgnd += v; else sumLsheet += v;
                return;
            }
            sum_H2_dl_L += sumL[ci] * slabReactanceFactor(signalSlab(c, ci), deltaCond(c));
        });
        if (sumSig > 0) {
            const Lac = sumSig / omega;
            sum_H2_dl_L += omega / Math.sqrt(1 / (Ldc * Ldc) + 1 / (Lac * Lac)) + sumExcess;
        }
        sum_H2_dl_L += sumLDefault * slabReactanceFactor(Math.abs(this.t) / 2);

        // Half-domain solve: the surface integral covered only x >= 0. The
        // mirror half contributes the same power. After this the sums mean
        // "full domain" again, so the existing power_factor convention applies
        // unchanged.
        if (this.sym_half) {
            sum_H2_dl_R *= 2;
            for (let c = 0; c < gndR.length; c++) gndR[c] *= 2;
            sum_H2_dl_L *= 2;
            sumLgnd *= 2;
            sumLsheet *= 2;
        }

        // Power normalization factor: differential has 0.5 factor
        // This is because we integrate over both traces but report normalized loss
        const power_factor = (this.is_differential && line === null) ? 0.5 : 1.0;

        // Vacuum variant: |H|^2 per unit current is |H_vac per 1V|^2*Z0_vac^2, since the
        // vacuum drive at 1V carries I_vac = 1/Z0_vac (legacy: same algebra with the
        // line Z0).
        const Z0_sq = Z0 * Z0;

        // AC Resistance per unit length from skin effect (Ohm/m)
        const R_ac_sig = power_factor * sum_H2_dl_R * Z0_sq;
        // Ground resistance with the spreading factors above.
        let sum_H2_dl_Rgnd = 0.0, sumRsheet = 0.0, spread = 0;
        for (let c = 0; c < gndR.length; c++) {
            if (!(gndS2[c] > 0)) continue;
            if (idealGnd.has(c)) { sum_H2_dl_Rgnd += gndR[c] * spreadG[c]; continue; }
            sumRsheet += gndR[c] * spreadG[c];
            spread = Math.max(spread, spreadU[c]);
        }
        // Once the spreading length delta^2/d of a ground passes a fraction of its return
        // width W_K, the thin-sheet solve (_ground_sheet_impedance) takes over the
        // resistance and inductance of the grounds that are not walls (_wall_grounds) from
        // the surface model: fully above 0.2, blended on log(spread) down to 0.02, where
        // the two agree to a few tenths of a percent.
        const wSheet = Math.min(1, Math.max(0, Math.log10(spread / 0.02)));
        const sheetZ = vacuum_fields && omega > 0 && wSheet > 0
            ? this._ground_sheet_cached(line === null ? mode : null) : null;
        if (sheetZ) {
            // Full-domain currents: a half domain holds half of each.
            const m = this.sym_half ? 2 : 1, I0 = m * netI[0], I1 = m * netI[1];
            const quad = X => I0 * I0 * X[0][0] + 2 * I0 * I1 * X[0][1] + I1 * I1 * X[1][1];
            sum_H2_dl_Rgnd += (1 - wSheet) * sumRsheet + wSheet * quad(sheetZ.R);
            sum_H2_dl_L += sumLgnd + (1 - wSheet) * sumLsheet + wSheet * omega * quad(sheetZ.L);
        } else {
            sum_H2_dl_Rgnd += sumRsheet;
            sum_H2_dl_L += sumLgnd + sumLsheet;
        }
        const R_ac_gnd = power_factor * sum_H2_dl_Rgnd * Z0_sq;
        const R_ac = R_ac_sig + R_ac_gnd;

        // Bounded below the skin regime by the DC solve of the traces, and by the slab
        // factor or the thin-sheet solve of the grounds above.
        const L_internal = power_factor * sum_H2_dl_L * Z0_sq / (2 * Math.PI * this.freq);

		// DC-skin transition correction (vacuum-field path only): against
		// tri-MQS the sqrt(R_dc^2+R_ac^2) is consistently high, a log-normal
		// notch in δ/t, = −7% at δ/t = 0.4, gone below δ/t = 0.12 and decaying
		// by δ/t = 1.3 (where R_dc takes over). Calibrated on ms/sl (w/h
		// 0.125-1.9, εr 4.4-9.8, t 17-70 µm, σ 1e6-5.8e7, f 0.25-4 GHz). A δ/t
		// sweep (0.12-1.3) over microstrip, stripline, diff microstrip and
		// narrow-/wide-gap GCPW showed the same ~7% bump at δ/t 0.3-0.4 in
		// every family, so the notch applies to all line types.
        // Alternatives measured and rejected: no
        // notch (microstrip +10% mid-transition), and slab/p-norm blends
        // R_dc*Re[q*coth q] or (R_dc^4+R_ac^5)^0.25.
        let transitionCal = 1.0;
        // Cleared unconditionally: the warning describes this call's frequency, and
        // the branch below is skipped on the legacy integrand and on solvers without
        // a rectangular conductor thickness. Leaving it set would carry a stale note
        // into the next sweep point.
        this._skinTransitionWarn = null;
        this._platingTransitionWarn = this._plating_transition_note(this.freq);
        // |t|: an embedded trace (negative thickness) has the same skin transition.
        const tAbs = Math.abs(this.t);
        if (vacuum_fields && tAbs > 0 && this.freq > 0) {
            const lx = Math.log(delta / tAbs / 0.4);
            transitionCal = 1 - 0.07 * Math.exp(-(lx * lx) / (2 * 0.45 * 0.45));
            const tMin = this.t_gnd > 0 ? Math.min(tAbs, this.t_gnd) : tAbs;
            // reason distinguishes loss-accuracy notes from certificate notes for
            // machine consumers (the fuzzer relaxes its R gate on loss reasons).
            this._skinTransitionWarn = (delta > 0.5 * tMin)
                ? { type: 'accuracy', reason: 'skin-transition', mode: 'all', message:
                    `Conductor loss and internal-inductance accuracy is reduced in the DC-skin ` +
                    `transition (skin depth ${(delta * 1e6).toFixed(1)} µm vs conductor thickness ` +
                    `${(tMin * 1e6).toFixed(1)} µm). The full-wave solver resolves the ` +
                    `transition-region current accurately.` }
                : null;
        }
        // The signal blends its DC and skin terms; the ground term is already valid at
        // every delta and only floors at its geometric DC resistance.
        const R_total = Math.max(R_dc_sig, transitionCal * Math.sqrt(R_dc_sig * R_dc_sig + R_ac_sig * R_ac_sig))
            + Math.max(R_dc_gnd, R_ac_gnd);

        return { R_ac, R_dc, R_total, L_internal };
    }

    // Conductor loss for a solved mode, choosing the integrand variant:
    // rect-based solvers (conductor_id present) use the vacuum-field integrand
    // when the mode's vacuum fields are available.
    // The only production caller lacking vacuum fields on a rect solver is
    // _solve_single_mode(vacuum_first=false), whose loss output is discarded
    // and recomputed by the caller with the cached vacuum fields.
    // Plating through the whole cross-section (at least as thick as the conductor, or
    // filling its width) makes it plating metal; the layered plating-over-bulk
    // impedance has no bulk to stand on.
    _solid_plating(cond) {
        return !!cond && platedThrough(cond);
    }

    // Accuracy note for plated conductors in the skin transition: the layered
    // plating-over-bulk surface impedance assumes a bulk thick against its skin
    // depth. Returns null when no plated rectangular conductor has less than two
    // bulk skin depths under its plating (solid plating is exact by convention).
    //   meshedThick - thick plating was solved as meshed metal and is left out
    //   fullWave    - the full-wave solver, which can mesh the plating
    _plating_transition_note(f, { meshedThick = false, fullWave = false } = {}) {
        if (!(f > 0) || !this.conductors) return null;
        let worst = null;
        for (const c of this.conductors) {
            const pl = c.plating;
            if (c.shape || !pl || !(pl.sigma > 0) || !(pl.top || pl.sides || pl.bottom)) continue;
            if (meshedThick && pl.thick_corners) continue;
            if (this._solid_plating(c)) continue;
            const tp = pl.thickness ?? 0, t = Math.abs(c.height);
            const bulk = t - tp;
            const delta = Math.sqrt(2 / (2 * Math.PI * f * 4e-7 * Math.PI * (c.sigma > 0 ? c.sigma : this.sigma_cond)));
            if (bulk < 2 * delta && (!worst || bulk < worst.bulk)) worst = { bulk, t, tp, delta, thick: !!pl.thick_corners };
        }
        if (!worst) return null;
        const delta = worst.delta;
        return { type: 'accuracy', reason: 'plating-transition', mode: 'all', message:
            `Plating accuracy is reduced: the ${(worst.tp * 1e6).toFixed(2)} µm plating on a ` +
            `${(worst.t * 1e6).toFixed(2)} µm conductor leaves ${(worst.bulk * 1e6).toFixed(2)} µm of bulk metal ` +
            `under it, less than two skin depths (${(delta * 1e6).toFixed(2)} µm) at this frequency. ` +
            `The layered surface impedance assumes a thick bulk, so conductor loss and internal ` +
            `inductance can be off by tens of percent here.` +
            (fullWave && !worst.thick ? ' Model Thick Plating (Advanced Options) solves the plating as a ' +
                'layer of metal on the full-wave solver.' : '') };
    }

    // DC conductivity of the signal metal: the plating's when every signal
    // conductor is solid plating, else the bulk.
    _signal_sigma(line = null) {
        const sig = (this.conductors || []).filter(c => c.is_signal
            && (line === null || (c.polarity < 0) === (line === 1)));
        if (sig.length && sig.every(c => this._solid_plating(c))) return sig[0].plating.sigma;
        if (!this._own_sigma()) return this.sigma_cond;
        // Signal conductors of different metals: the area-weighted mean.
        let a = 0, ga = 0;
        for (const c of sig) {
            const area = Math.abs(c.width * c.height);
            a += area; ga += this._bulk_sigma(c) * area;
        }
        return a > 0 ? ga / a : this.sigma_cond;
    }

    // True when a conductor carries a conductivity of its own (custom geometry).
    _own_sigma() {
        return (this.conductors || []).some(c => c.sigma > 0);
    }

    // Conductivity of a conductor's cross-section.
    _bulk_sigma(c) {
        if (this._solid_plating(c)) return c.plating.sigma;
        return c.sigma > 0 ? c.sigma : this.sigma_cond;
    }

    // Per-unit-length DC conductance of the positive traces, the negative traces and
    // the grounds.
    _dc_conductances() {
        const g = { pos: 0, neg: 0, gnd: 0 };
        const visible = visibleAreas(this.conductors || []);
        for (const [i, c] of (this.conductors || []).entries()) {
            // A plating layer inside the outline conducts at its own sigma.
            const a = visible && visible.has(i) ? visible.get(i) : Math.abs(c.width * c.height);
            const ap = this._solid_plating(c) ? 0 : Math.min(platingArea(c), a);
            const v = this._bulk_sigma(c) * (a - ap) + (ap > 0 ? c.plating.sigma * ap : 0);
            if (!c.is_signal) g.gnd += v; else if (c.polarity < 0) g.neg += v; else g.pos += v;
        }
        return g;
    }

    // Internal inductance at DC: the surface integral evaluated where the skin depth
    // is 100x the thickest conductor, so every face sits on its mu0 d/3 plateau.
    _dc_internal_inductance(...args) {
        let tMax = 0;
        for (const c of this.conductors || []) tMax = Math.max(tMax, Math.abs(c.height));
        if (!(tMax > 0)) return 0;
        const delta = 100 * tMax;
        let sigma = this.sigma_cond;
        for (const c of this.conductors || []) if (c.sigma > sigma) sigma = c.sigma;
        const fDc = 2 / (2 * Math.PI * 4e-7 * Math.PI * sigma * delta * delta);
        const f0 = this.freq;
        this.freq = fDc;
        try { return this.calculate_conductor_loss(...args).L_internal; }
        finally { this.freq = f0; }
    }

    // Per-line loss data of a pair (line 1 = the positive trace), null when not needed or
    // the vacuum fields are missing. The surface integral of the vacuum field of given
    // trace currents I is I^T R I, and I^T L_int I for the reactance.
    //   - mirror-symmetric geometry whose traces differ in metal or finish: { dR, dL } =
    //     R11 - R22 and L11 - L22. The field of unit current in one trace is half the sum
    //     (difference) of the even and odd vacuum fields, each scaled to unit current per
    //     trace by its vacuum impedance. The modes supply R11 + R22 and R12.
    //   - asymmetric geometry (physMatrix): { Rm, Lim } = [X11, X12, X22]. The per-trace
    //     vacuum solves are combined with the voltages V = Cm0^-1 I / c of the currents
    //     [1, 0], [0, 1] and [1, 1], the last giving X12.
    _line_asymmetry(odd, even) {
        if (!this.is_differential || this.sym_half || !this.conductor_id) return null;
        const warn = [this._skinTransitionWarn, this._platingTransitionWarn];
        const mix = (a, A, b, B) => A.map((row, i) => {
            const out = new Float64Array(row.length), rb = B[i];
            for (let j = 0; j < row.length; j++) out[j] = a * row[j] + b * rb[j];
            return out;
        });
        try {
            const tv = this._modalPhys ? this._traceVac : null;
            if (tv) {
                const M = tv.Cm0, det = M[0][0] * M[1][1] - M[0][1] * M[1][0], k = 1 / (CONSTANTS.C * det);
                const X = [[1, 0], [0, 1], [1, 1]].map((I, line) => {
                    const va = k * (M[1][1] * I[0] - M[0][1] * I[1]), vb = k * (M[0][0] * I[1] - M[1][0] * I[0]);
                    return this.calculate_conductor_loss(mix(va, tv.Av.Ex, vb, tv.Bv.Ex), mix(va, tv.Av.Ey, vb, tv.Bv.Ey),
                        1, true, null, line);
                });
                const m3 = key => [X[0][key], (X[2][key] - X[0][key] - X[1][key]) / 2, X[1][key]];
                return { Rm: m3('R_total'), Lim: m3('L_internal') };
            }
            if (this._modalPhys || !this._pair_finish_differs()) return null;
            if (!odd || !even || !odd.Ex0 || !even.Ex0 || !(odd.C0 > 0) || !(even.C0 > 0)) return null;
            const zo = 1 / (CONSTANTS.C * odd.C0), ze = 1 / (CONSTANTS.C * even.C0);
            const lines = [1, -1].map((sb, line) => this.calculate_conductor_loss(
                mix(0.5 * ze, even.Ex0, 0.5 * sb * zo, odd.Ex0), mix(0.5 * ze, even.Ey0, 0.5 * sb * zo, odd.Ey0), 1, true, null, line));
            return { dR: lines[0].R_total - lines[1].R_total, dL: lines[0].L_internal - lines[1].L_internal };
        } finally {
            [this._skinTransitionWarn, this._platingTransitionWarn] = warn;
        }
    }

    _pair_finish_differs() {
        if (this._finish_differs === undefined) {
            const keys = neg => [...new Set((this.conductors || []).filter(c => c.is_signal && (c.polarity < 0) === neg)
                .map(conductorFinishKey))].sort().join(';');
            this._finish_differs = keys(true) !== keys(false);
        }
        return this._finish_differs;
    }

    _mode_conductor_loss(Ex, Ey, Z0, C0, Ex0, Ey0, mode = null) {
        if (this.conductor_id && Ex0 && Ey0 && C0 > 0) {
            const Z0_vac = 1 / (CONSTANTS.C * C0);
            return this.calculate_conductor_loss(Ex0, Ey0, Z0_vac, true, mode);
        }
        return this.calculate_conductor_loss(Ex, Ey, Z0, false, mode);
    }

    calculate_dielectric_loss(V, Z0) {
        if (!this.solution_valid) {
            throw new Error("Potential V is not valid. Run the solve first.");
        }

        // No dielectric loss at DC
        // If material conductivity is implemented this is not true
        if (this.freq === 0) {
            return 0;
        }
        const Pd = this._dielectric_power(V, 2 * Math.PI * this.freq);

        // Power normalization factor
        const power_factor = this.is_differential ? 0.5 : 1.0;
        const P_flow = 1.0 / (2 * Z0);
        return 8.686 * (power_factor * Pd / (2 * P_flow));
    }

    // Power dissipated in the dielectrics by the potential V at angular frequency omega
    // (full domain). For a pair, 2 * P of the potential of drive v is v^T G v.
    _dielectric_power(V, omega) {
        const ny = this.y.length;
        const nx = this.x.length;
        const dx_array = diff(this.x);
        const dy_array = diff(this.y);

        // Cell-wise integral over the operator's cell-centred materials (see
        // _paint_cell_materials): the field in a cell is the mean of its edge
        // differences of V, so a cell next to a dielectric interface weights its own
        // side's field with its own eps*tand. Sampling nodal material at the cell
        // corner instead puts one region's eps*tand on an interface node against a
        // field that belongs mostly to the other side (a lossy cover over an air layer
        // read G 30% high). Cells inside a conductor have a constant V and drop out.
        let Pd = 0.0;
        const wPd = 0.5 * omega * CONSTANTS.EPS0;
        for (let i = 0; i < ny - 1; i++) {
            const V0 = V[i], V1 = V[i + 1];
            const ec = this.epsilon_cell[i], tc = this.tand_cell[i];
            const dy = dy_array[i];
            for (let j = 0; j < nx - 1; j++) {
                const w = ec[j] * tc[j];
                if (w === 0) continue;
                const dx = dx_array[j];
                const Ex = -0.5 * ((V0[j + 1] - V0[j]) + (V1[j + 1] - V1[j])) / dx;
                const Ey = -0.5 * ((V1[j] - V0[j]) + (V1[j + 1] - V0[j + 1])) / dy;
                Pd += wPd * w * (Ex * Ex + Ey * Ey) * dx * dy;
            }
        }
        // Half-domain solve: the cells cover x >= 0 of a mirror-symmetric field.
        if (this.sym_half) Pd *= 2;
        return Pd;
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

    // Adaptive Meshing
    // Discrete flux-jump refinement indicator (FDM analog of the tri backend's
    // Kelly marker). At each interior node the jump in ε*dV/dh between the two
    // adjacent intervals measures local discretization error. It vanishes where
    // the discrete solution resolves the field and concentrates where curvature
    // is under-resolved. The legacy dV*|E|*ε metric is a field intensity
    // indicator, it keeps refining strong-field intervals. Jumps at
    // conductor-adjacent nodes are physical surface charge and are skipped, the
    // same rule as the tri Kelly indicator's Dirichlet-edge skip.
    //
    // planeBC (half-domain solves): the symmetry-plane BC of this field, see the
    // j = 0 term below.
    _compute_refine_metrics_jump(V, vacuum = false, planeBC = null) {
        const ny = this.y.length, nx = this.x.length;
        const x_metrics = new Float64Array(nx - 1);
        const y_metrics = new Float64Array(ny - 1);
        const dx = diff(this.x), dy = diff(this.y);
        const cm = this.conductor_mask;
        // Face permittivities exactly as the operator forms them (see
        // solve_laplace_multi). The cell-weighted average over the half-faces
        // above and below (x faces) or left and right (y faces) of the node. The
        // nodal average this used to take differs from the operator's face value
        // at any interface running along the jump direction, which shows up as a
        // spurious jump where the discrete flux is in fact continuous, the
        // indicator then spends refinement on interfaces that are already exact.
        const ec = vacuum ? null : this.epsilon_cell;
        const epsX = (i, cj) => {           // face at column cj, node row i
            if (!ec) return 1.0;
            const hd = i > 0 ? dy[i - 1] : dy[0], hu = i < ny - 1 ? dy[i] : dy[ny - 2];
            const cid = i > 0 ? i - 1 : 0, ciu = i < ny - 1 ? i : ny - 2;
            return (ec[cid][cj] * hd + ec[ciu][cj] * hu) / (hd + hu);
        };
        const epsY = (ci, j) => {           // face at row ci, node column j
            if (!ec) return 1.0;
            const wl = j > 0 ? dx[j - 1] : dx[0], wr = j < nx - 1 ? dx[j] : dx[nx - 2];
            const cjl = j > 0 ? j - 1 : 0, cjr = j < nx - 1 ? j : nx - 2;
            return (ec[ci][cjl] * wl + ec[ci][cjr] * wr) / (wl + wr);
        };

        for (let i = 0; i < ny; i++) {
            // transverse control width so the accumulation approximates an area integral
            const wy = ((i + 1 < ny ? this.y[i + 1] : this.y[i]) - (i > 0 ? this.y[i - 1] : this.y[i])) / 2;
            for (let j = 1; j < nx - 1; j++) {
                if (cm[i][j] || cm[i][j - 1] || cm[i][j + 1]) continue;
                const eL = epsX(i, j - 1);
                const eR = epsX(i, j);
                const J = Math.abs(eR * (V[i][j + 1] - V[i][j]) / dx[j] -
                                   eL * (V[i][j] - V[i][j - 1]) / dx[j - 1]);
                x_metrics[j - 1] += J * dx[j - 1] * wy;
                x_metrics[j] += J * dx[j] * wy;
            }
            // Half-domain symmetry plane (j = 0) with a magnetic wall: the mirrored
            // ghost column (V[-1] = V[1], same spacing and epsilon) gives the
            // full-domain centre node's jump, credited to the one interval it
            // borders here. The electric wall (odd mode) pins V[0] = 0 with an odd
            // mirror, so its jump is 0.
            if (planeBC === 'pmc' && nx > 1 && !cm[i][0] && !cm[i][1]) {
                const J = 2 * epsX(i, 0) * Math.abs(V[i][1] - V[i][0]) / dx[0];
                x_metrics[0] += J * dx[0] * wy;
            }
        }
        for (let j = 0; j < nx; j++) {
            const wx = ((j + 1 < nx ? this.x[j + 1] : this.x[j]) - (j > 0 ? this.x[j - 1] : this.x[j])) / 2;
            for (let i = 1; i < ny - 1; i++) {
                if (cm[i][j] || cm[i - 1][j] || cm[i + 1][j]) continue;
                const eD = epsY(i - 1, j);
                const eU = epsY(i, j);
                const J = Math.abs(eU * (V[i + 1][j] - V[i][j]) / dy[i] -
                                   eD * (V[i][j] - V[i - 1][j]) / dy[i - 1]);
                y_metrics[i - 1] += J * dy[i - 1] * wx;
                y_metrics[i] += J * dy[i] * wx;
            }
        }
        return { x_metrics, y_metrics };
    }

    _compute_refine_metrics(V, Ex, Ey, vacuum = false, planeBC = null) {
        /**
         * For each grid interval, compute a metric indicating how much
         * refinement would help. Default ('blend'): the flux-jump error
         * indicator plus the legacy dV*|E|*ε field-intensity metric. The jump
         * metric targets the static discretization error (C, C0).  But the
         * intensity metric's conductor-surface emphasis is important for the
         * conductor-loss integral. The blend keeps both. Opt-outs:
         * solver.refine_metric = 'intensity' (legacy) or 'jump' (static-only).
         *
         * vacuum: the fields are from the vacuum (C0) solve.
         * planeBC: the field's symmetry-plane BC on a half-domain solve.
         */
        const metric = this.refine_metric ?? 'blend';
        if (metric === 'jump') {
            return this._compute_refine_metrics_jump(V, vacuum, planeBC);
        }
        if (metric === 'blend') {
            // Surface component restricted to conductor-adjacent cells.
            const alpha = this.refine_surface_weight ?? 0.35;
            const j = this._compute_refine_metrics_jump(V, vacuum, planeBC);
            const f = this._compute_refine_metrics_intensity(V, Ex, Ey, vacuum, true);
            const norm = (m) => {
                const t = m.x_metrics.reduce((s, v) => s + v, 0) + m.y_metrics.reduce((s, v) => s + v, 0);
                return t > 0 ? 1 / t : 0;
            };
            const nj = norm(j), nf = norm(f) * alpha;
            const x_metrics = j.x_metrics.map((v, k) => v * nj + f.x_metrics[k] * nf);
            const y_metrics = j.y_metrics.map((v, k) => v * nj + f.y_metrics[k] * nf);
            return { x_metrics, y_metrics };
        }
        return this._compute_refine_metrics_intensity(V, Ex, Ey, vacuum);
    }

    // Legacy field-intensity metric: dV * |E| * ε with a x2 conductor-boundary
    // boost. Not an error indicator (it keeps refining strong-field intervals
    // after they converge) but its surface emphasis resolves the conductor-loss
    // integrand.
    _compute_refine_metrics_intensity(V, Ex, Ey, vacuum = false, boundaryOnly = false) {
        const ny = V.length;
        const nx = V[0].length;
        const cm = this.conductor_mask;

        // Metric for splitting interval [x[j], x[j+1]]
        const x_metrics = new Float64Array(this.x.length - 1);
        // Metric for splitting interval [y[i], y[i+1]]
        const y_metrics = new Float64Array(this.y.length - 1);

        // Boundary detection: a dielectric node next to a conductor
        const isBoundary = (i, j) => !cm[i][j] && (
            (i > 0 && cm[i - 1][j]) || (i < ny - 1 && cm[i + 1][j]) ||
            (j > 0 && cm[i][j - 1]) || (j < nx - 1 && cm[i][j + 1]));

        // Sample of the cell with lower corners (i, js) [the sampled one] and (i, jo):
        // the field-intensity weight and the vertical voltage step at the sample
        // node, or null when the cell is skipped.
        const sample = (i, js, jo) => {
            // Skip cells fully inside conductors
            if (cm[i][js] && cm[i + 1][js] && cm[i][jo]) return null;
            if (boundaryOnly && !isBoundary(i, js)) return null;
            const eps = vacuum ? 1.0 : this.epsilon_r[i][js];
            const E2 = Ex[i][js] ** 2 + Ey[i][js] ** 2;
            const E_mag = E2 > 0 ? Math.sqrt(E2) : 1e-12;
            // Weight by field strength, permittivity, and boundary importance
            const weight = E_mag * eps * (isBoundary(i, js) ? 2.0 : 1.0);
            return { weight, dV_y: Math.abs(V[i + 1][js] - V[i][js]) };
        };

        for (let i = 0; i < ny - 1; i++) {
            for (let j = 0; j < nx - 1; j++) {
                // Voltage difference across this cell, the cell is sampled at its
                // lower-left corner.
                const dV_x = Math.abs(V[i][j + 1] - V[i][j]);
                const sR = sample(i, j, j + 1);
                if (!this.sym_half) {
                    if (!sR) continue;
                    x_metrics[j] += dV_x * sR.weight;
                    y_metrics[i] += sR.dV_y * sR.weight;
                    continue;
                }
                // Half-domain solve: the mirror image of this cell samples the mirror
                // of the right corner. Averaging both corners (each with its own skip
                // test) makes the x metric equal the full domain's pair-averaged value
                // and 2x the y metric the full-domain sum, so the refinement trajectory
                // is the full-domain one (see _compute_refine_metrics_jump's plane term).
                const sL = sample(i, j + 1, j);
                if (sR) { x_metrics[j] += 0.5 * dV_x * sR.weight; y_metrics[i] += 0.5 * sR.dV_y * sR.weight; }
                if (sL) { x_metrics[j] += 0.5 * dV_x * sL.weight; y_metrics[i] += 0.5 * sL.dV_y * sL.weight; }
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
        // A half-domain grid is not mirror-symmetric about its own center. A
        // coincidental match would mis-symmetrize the metrics about W/4.
        const x_symmetric = !this.sym_half && this._check_symmetry(this.x, x_center);

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

        // Size the pass as the full-domain grid would: a half-domain grid stands for
        // twice its x intervals, and each x interval selected there is a mirror pair
        // (the symmetric full-domain path selects both partners), so x gets half the
        // count. Keeps the half-domain refinement on the full-domain trajectory.
        const nxEq = this.sym_half ? 2 * x_metrics_proc.length : x_metrics_proc.length;
        let n_total = Math.floor(frac * (nxEq + y_metrics.length));
        n_total = Math.max(1, n_total);

        // Allocate proportionally to where the error is
        let n_x = Math.floor(n_total * total_x / total);
        const n_y = n_total - n_x;
        if (this.sym_half) n_x = Math.ceil(n_x / 2);

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
        const x_symmetric = !this.sym_half && this._check_symmetry(this.x, x_center);

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

    // Mesh refinement. solve_adaptive always goes through the multi-set form (one
    // set per mode plus one per mode's vacuum solve), so there is no single-set
    // variant, pass [{V, Ex, Ey}] for one field.
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

        for (const { V, Ex, Ey, vacuum, planeBC } of modes) {
            const { x_metrics, y_metrics } = this._compute_refine_metrics(V, Ex, Ey, vacuum, planeBC);
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

    _compute_energy_error(Ex, Ey, prev_energy, vacuum = false) {
        /**
         * Compute relative change in stored electromagnetic energy.
         * vacuum: score the vacuum (C0) field, which carries no permittivity.
         */
        const ny = this.y.length;
        const nx = this.x.length;
        const dx_array = diff(this.x);
        const dy_array = diff(this.y);

        // Each cell is sampled at its lower-left corner. On a half-domain solve the
        // mirrored cells sample the mirror of the right corner, so both corners are
        // averaged and the half sum doubled = the full-domain energy exactly (the
        // same rule as calculate_dielectric_loss), so the convergence decisions
        // follow the full-domain trajectory.
        const term = (i, j) => this.conductor_mask[i][j] ? 0
            : (vacuum ? 1 : this.epsilon_r[i][j]) * (Ex[i][j] ** 2 + Ey[i][j] ** 2);
        let energy = 0.0;
        for (let i = 0; i < ny - 1; i++) {
            for (let j = 0; j < nx - 1; j++) {
                const dA = dx_array[j] * dy_array[i];
                const t = this.sym_half ? 0.5 * (term(i, j) + term(i, j + 1)) : term(i, j);
                energy += 0.5 * CONSTANTS.EPS0 * t * dA;
            }
        }
        if (this.sym_half) energy *= 2;

        if (prev_energy === null || prev_energy === undefined) {
            return { energy, rel_error: 1.0 };
        }

        const rel_error = Math.abs(energy - prev_energy) / Math.max(Math.abs(prev_energy), 1e-12);
        return { energy, rel_error };
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
        // Half-domain pair: symmetric by construction and only one trace is meshed, so
        // the odd/even drives are the modes and the per-trace solves are impossible.
        if (this.sym_half) { this._modalPhys = null; this._traceVac = null; return null; }
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
        // Per-trace vacuum fields for the per-line loss matrices (_line_asymmetry), and
        // the dielectric loss matrix over omega, [G11, G12, G22] / omega: v^T G v is twice
        // the power of the potential of drive v, so G12 comes from the drive [1, 1].
        const pA = this._dielectric_power(A.V, 1), pB = this._dielectric_power(B.V, 1);
        const pAB = this._dielectric_power(A.V.map((row, i) => row.map((v, j) => v + B.V[i][j])), 1);
        this._traceVac = { Av, Bv, Cm0, Gw: [2 * pA, pAB - pA - pB, 2 * pB] };
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
        await this._ensure_dc_signal_inductance(null);
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
            const { R_total, L_internal } = this._mode_conductor_loss(Ex, Ey, Z0, C0k, Ex0, Ey0, label);
            const alpha_d = this.calculate_dielectric_loss(V, Z0);
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
        // Half-domain pair: only one trace is painted, so signal_mask already is that
        // trace's contour and carries the per-trace charge directly (the 0.5*(Cp+Cn)
        // average below would see an empty partner mask and silently halve C).
        if (!this.is_differential || this.sym_half) {
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

    async _solve_single_mode(mode, vacuum_first = true, withLoss = true) {
        /**
         * Solve a single mode and return full results.
         *
         * Parameters:
         * -----------
         * mode : string - 'single', 'odd', or 'even'
         * vacuum_first : boolean - Whether to solve vacuum case first for C0 calculation
         * withLoss : boolean - false returns the fields, C, C0 and Z0 only
         *                      (_mode_loss_results completes them)
         *
         * Returns:
         * --------
         * {mode, Z0, eps_eff, C, C0, RLGC, Zc, alpha_c, alpha_d, alpha_total, V, Ex, Ey}
         */
        let C0;
        let V;
        let V0, Ex0, Ey0;
        const planeBC = this._plane_bc(mode);

        if (vacuum_first) {
            // Calculate C0 (vacuum capacitance)
            V0 = this._create_voltage_array(mode);
            V0 = await this.solve_laplace(V0, true, planeBC);
            C0 = this._signal_capacitance(V0, true);
            // Vacuum fields: the conductor-loss integrand for rect-based solvers
            // (the quasi-TEM H pattern, see _mode_conductor_loss). Frequency- and
            // ε-independent, so cached mode results reuse them across the sweep.
            // V0 is kept as well: adaptive refinement feeds the vacuum solve
            // into its metrics so C0 (which the certificate certifies)
            // converges alongside C, see solve_adaptive.
            if (this.conductor_id) ({ Ex: Ex0, Ey: Ey0 } = this.compute_fields(V0, planeBC));
        }

        // Solve with dielectric
        V = this._create_voltage_array(mode);
        V = await this.solve_laplace(V, false, planeBC);
        const C = this._signal_capacitance(V, false);

        // Calculate fields
        const { Ex, Ey } = this.compute_fields(V, planeBC);

        // Calculate impedance
        let Z0;
        if (C0 !== undefined) Z0 = 1 / (CONSTANTS.C * Math.sqrt(C * C0));

        const r = { mode, Z0, C, C0, V, Ex, Ey, V0, Ex0, Ey0 };
        return withLoss ? this._mode_loss_results(r) : r;
    }

    // Losses, RLGC and the reported parameters of a mode solved by _solve_single_mode.
    async _mode_loss_results(r) {
        const { mode, Z0, C, C0, V, Ex, Ey, V0, Ex0, Ey0 } = r;
        // Calculate conductor losses with surface roughness and DC resistance
        await this._ensure_dc_signal_inductance(this._plane_bc(mode));
        const { R_total, L_internal } = this._mode_conductor_loss(Ex, Ey, Z0, C0, Ex0, Ey0, mode);

        // Calculate dielectric loss (returns alpha in dB/m)
        const alpha_d = this.calculate_dielectric_loss(V, Z0);

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
            V0, Ex0, Ey0
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
    // Field region of an auto-sized open domain for the Modes solve: the signal
    // cluster padded by four substrate-stack heights on the open sides, the ground
    // side kept. Returns null when the box would not be smaller than the domain, or
    // when the domain is user-sized (an enclosure is a physical boundary).
    _modes_domain_box() {
        if (this.enclosure_width != null || this.enclosure_height != null) return null;
        if (!this.conductors || !this.domain_width || !(this.domain_height > 0)) return null;
        const W = this.domain_width, H = this.domain_height, yBottom = this.domain_y_min;
        const stack = [...this.conductors,
            ...(this.dielectrics || []).filter(d => (d.epsilon_r || 1) > 1.001)];
        let yLo = Infinity, yHi = -Infinity, xLo = Infinity, xHi = -Infinity;
        for (const r of stack) { yLo = Math.min(yLo, r.y_min); yHi = Math.max(yHi, r.y_max); }
        for (const c of this.conductors) {
            if (Math.abs(c.width) >= W * 0.99) continue;
            xLo = Math.min(xLo, c.x_min); xHi = Math.max(xHi, c.x_max);
        }
        if (!(yHi > yLo) || !(xHi > xLo)) return null;
        const G = yHi - yLo, pad = 4 * G;
        const b = (this.boundaries || ['open', 'open', 'open', 'gnd']);
        const box = {
            x_min: b[0] === 'open' ? Math.max(-W / 2, xLo - pad) : -W / 2,
            x_max: b[1] === 'open' ? Math.min(W / 2, xHi + pad) : W / 2,
            y_min: yBottom,
            y_max: b[2] === 'open' ? Math.min(H, yHi + pad) : H,
        };
        const shrunk = (box.x_max - box.x_min) * (box.y_max - box.y_min) < 0.95 * W * (H - yBottom);
        return shrunk ? box : null;
    }

    async solveModes(freq, nev = 4, onProgress = null, refineOpts = {}) {
        // Optional shrink of an auto-sized open domain to the field region: the modes
        // solve wavelength-resolves the whole meshed box, and the padded far field of an
        // open line costs triangles for modes of the artificial box only.
        const domainBox = refineOpts.shrinkDomain ? this._modes_domain_box() : null;
        // Same pre-mesh guard as solve_adaptive, evaluated at the MODES frequency and
        // ALWAYS on the triangular estimate — this solve runs the triangular backend
        // regardless of the sidebar Solver selection, so branching on this.mesh_backend
        // here would apply the FDM tensor-grid rules (falsely rejecting fine features
        // the tri mesher grades locally, and skipping the whole-domain wavelength check
        // that stops an electrically huge modes mesh from hanging gmsh).
        this._check_meshability(refineOpts.maxNodes ?? 20000, freq,
            { wavelengthDensity: refineOpts.wavelengthDensity, domainBox });
        const { initTriBackend, TriBackend } = await import('./tri_solver/tri_backend.js');
        const ctx = await initTriBackend();
        // modesFreq lets buildMesh size the bulk to the wavelength at this frequency, so
        // high-frequency cavity/higher-order modes are resolved (and not mis-flagged as
        // spurious by the mesh-convergence test).
        // Mesh.Algorithm 1 (MeshAdapt) gives the cleanest eigenmode spectrum for the mode
        // viewer (fewer spurious low-ε_eff artifacts than the default frontal-Delaunay).
        const { shouldStop = null, shrinkDomain = false, ...meshOpts } = refineOpts;
        const opts = { ...(this.tri_opts || {}), symmetry: false, modesFreq: freq,
            gmshOptions: { 'Mesh.Algorithm': 1 }, ...meshOpts, domainBox };
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
        this._repaint_geometry();
        try {
            return await fn();
        } finally {
            this.x = saved_x;
            this.y = saved_y;
            this._repaint_geometry();
        }
    }

    // _setup_geometry repaints epsilon_r/tand from the nominal material values, so
    // every rebuild (initial mesh, each refinement pass, the certificate's grid
    // swaps) drops the causal dispersion and has to re-apply it. Always go through
    // this instead of calling _setup_geometry directly, or the solve silently runs
    // at the f_ref permittivity while a sweep point at the same frequency does not.
    // applyDjordjevicSarkar recomputes from cached nominals, so it is idempotent.
    _repaint_geometry() {
        this._setup_geometry();
        this._apply_causal_materials();
    }

    _apply_causal_materials() {
        if (this.use_causal_materials && this.epsilon_r && this.tand) {
            applyDjordjevicSarkar(this);
        }
    }

    // Static (C, C0) per mode — the quantities _solve_single_mode reports, without
    // the loss post-processing the certificate does not cover. All modes share the
    // signal-dielectric operator and the vacuum operator, so each matrix is
    // factored once and back-substituted per mode (solve_laplace_multi): a
    // differential pair does 2 factorizations per grid level instead of 4.
    async _staticCapacitances() {
        const modeNames = this.is_differential ? ['odd', 'even'] : ['single'];
        // One factorization for every mode, except on a half-domain solve where
        // the matrix carries the mode's symmetry-plane BC (odd: PEC, even/single:
        // PMC) and each mode solves alone.
        const solveSet = async (vacuum) => {
            const Vs = modeNames.map(m => this._create_voltage_array(m));
            if (!this.sym_half) return this.solve_laplace_multi(Vs, vacuum);
            const out = [];
            for (let i = 0; i < modeNames.length; i++)
                out.push(await this.solve_laplace(Vs[i], vacuum, this._plane_bc(modeNames[i])));
            return out;
        };
        const signal = await solveSet(false);
        const vacuum = await solveSet(true);
        const out = [];
        for (let i = 0; i < modeNames.length; i++) {
            out.push(this._signal_capacitance(signal[i], false));
            out.push(this._signal_capacitance(vacuum[i], true));
        }
        return out;
    }

    // Richardson certificate for the static quantities (C, C0 per mode).
    //
    // The comparison level halves one axis at a time rather than both at once.
    // For a tensor-product discretisation the error separates as
    // A_x*h_x^p + A_y*h_y^p, so the two single-axis differences add up to the
    // both-axes difference. Verified in tests to within 0.7% on a microstrip, a
    // stripline and a differential GCPW + solder mask, while the peak grid is
    // 2N instead of 4N. That is ~0.7x the solve time at level 1 and a
    // quarter of it at level 2 (two 4N solves instead of one 16N), and it halves
    // the peak factorization, which is what the node caps below really guard.
    //
    // Options beyond the caps/gates:
    //   * q0: base-grid (C, C0) per mode in _staticCapacitances order, when the
    //     caller already has them (the refinement pass that tripped the gate just
    //     computed exactly these). Skips re-solving the base level.
    //   * knownR: convergence ratio measured by an earlier certificate in the same
    //     solve. Skips the level-2 solves entirely: r decreases (or holds) under
    //     refinement, so an earlier genuine measurement used on a finer grid errs
    //     conservative, and the x1.5 safety still covers a moderately larger true r.
    //     Only genuinely measured, non-pre-asymptotic r values are reused.
    //   * l2MaxNodes: peak-grid cap for the level-2 (convergence-ratio) solves.
    //
    // All caps are on the peak grid a level actually solves, so they bound peak
    // memory directly.
    async _certifyStatic(tol, { maxNodes = 600000, l2MaxNodes = CERT_L2_MAX_NODES,
                                rMax = 0.7, safety = 1.5,
                                q0 = null, knownR = null } = {}) {
        // Relative distance from a base level to the pair of once-refined levels,
        // per quantity, worst quantity wins. The x and y components come out with
        // the same sign in practice, so adding their magnitudes is tight, and
        // conservative when they ever disagree.
        const splitDiff = (qBase, qx, qy) => {
            let m = 0;
            for (let i = 0; i < qBase.length; i++) {
                const dx = Math.abs(qBase[i] - qx[i]) / Math.max(Math.abs(qx[i]), 1e-300);
                const dy = Math.abs(qBase[i] - qy[i]) / Math.max(Math.abs(qy[i]), 1e-300);
                m = Math.max(m, dx + dy);
            }
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
        // Largest grid the split pair solves: one axis doubled, the other as is.
        const peakNodes = (x, y) => Math.max((2 * x.length - 1) * y.length,
                                             x.length * (2 * y.length - 1));
        if (peakNodes(this.x, this.y) > maxNodes) return null;
        if (!q0) q0 = await this._staticCapacitances();
        const x1 = bisect(this.x), y1 = bisect(this.y);
        const qx = await this._withGrid(x1, this.y, () => this._staticCapacitances());
        const qy = await this._withGrid(this.x, y1, () => this._staticCapacitances());
        const d1 = splitDiff(q0, qx, qy);
        const base = { nodes: this.x.length * this.y.length, d1, safety };
        if (d1 < 5e-5) return { ...base, pass: true, err: d1, r: 0, rSource: null, levels: 1 };

        // Level 2 refines each axis once more, again one at a time: qxx is the
        // x-halved grid with x halved again, qyy likewise in y, so each component
        // is a clean three-level sequence in its own direction.
        const l2Peak = Math.max((4 * this.x.length - 3) * this.y.length,
                                this.x.length * (4 * this.y.length - 3));
        const measureR = async () => {
            const qxx = await this._withGrid(bisect(x1), this.y, () => this._staticCapacitances());
            const qyy = await this._withGrid(this.x, bisect(y1), () => this._staticCapacitances());
            let d2 = 0;
            for (let i = 0; i < q0.length; i++) {
                const dx = Math.abs(qx[i] - qxx[i]) / Math.max(Math.abs(qxx[i]), 1e-300);
                const dy = Math.abs(qy[i] - qyy[i]) / Math.max(Math.abs(qyy[i]), 1e-300);
                d2 = Math.max(d2, dx + dy);
            }
            return d2 / d1;
        };

        if (d1 * safety >= tol) {
            // Already failing on d1 alone. Still measure r when the level-2 grids
            // are affordable: err = d1/(1-r) is more accurate than the d1 lower
            // bound, and the measured r is reused (knownR) by every later
            // certificate in this solve, otherwise the final certificate on a
            // large grid is stuck with the conservative r=0.5 fallback even when
            // the true ratio is ~0.3. It runs at most once per adaptive solve
            // (knownR caches the result).
            if (knownR == null && l2Peak <= l2MaxNodes) {
                const r = await measureR();
                if (r >= rMax)
                    return { ...base, pass: false, err: d1, r, levels: 2, rSource: 'measured', preAsymptotic: true };
                return { ...base, pass: false, err: d1 / (1 - r), r, levels: 2, rSource: 'measured' };
            }
            const rEst = knownR != null ? knownR : null;
            return {
                ...base, pass: false,
                err: rEst != null ? d1 / (1 - Math.min(rEst, rMax)) : d1,
                r: rEst, levels: 1, rSource: rEst != null ? 'reused' : null,
            };
        }
        let r = knownR != null ? knownR : 0.5;
        let rSource = knownR != null ? 'reused' : 'fallback';
        let levels = 1;
        let err = d1 / (1 - Math.min(r, rMax));
        // Measuring r can only shrink the estimate, so when the conservative
        // r = 0.5 fallback already clears the tolerance the level-2 solves would
        // only improve the already passing estimate. Pay for it only when the
        // alternative is continuing with extra refinement pass.
        if (err * safety >= tol && knownR == null && l2Peak <= l2MaxNodes) {
            r = await measureR();
            levels = 2;
            rSource = 'measured';
            if (r >= rMax)
                return { ...base, pass: false, err: d1, r, levels, rSource, preAsymptotic: true };
            err = d1 / (1 - Math.min(r, rMax));
        }
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
        this._apply_causal_materials();

        if (this.sym_half) {
            console.log('x=0 mirror symmetry detected: meshing right half only');
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
            certify_l2_max_nodes = CERT_L2_MAX_NODES,
            onProgress = null,
            shouldStop = null,
            skip_mesh = false
        } = options;

        // If skip_mesh is true, just solve once with existing mesh
        if (skip_mesh) {
            let modeResults;
            if (this.is_differential) {
                modeResults = [await this._solve_single_mode('odd', true), await this._solve_single_mode('even', true)];
                const modal = await this._solve_modal_differential();   // null for symmetric/degenerate/half-domain
                if (modal) modeResults = modal;
            } else {
                modeResults = [await this._solve_single_mode('single', true)];
            }
            // Store fields as arrays
            this.V = modeResults.map(r => r.V);
            this.Ex = modeResults.map(r => r.Ex);
            this.Ey = modeResults.map(r => r.Ey);
            this._plotModes = modeResults.map(r => r.mode);
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
        const prevEnergy0 = {};
        const prevZ0 = {};
        let converged_count = 0;
        let modeResults = null;

        for (let it = 0; it < max_iters; it++) {
            // Solve all modes. Refinement and its convergence gate need the fields,
            // C, C0 and Z0 only: the losses, with the DC trace solve and the ground
            // sheet setup of each grid, follow once on the final grid.
            modeResults = [];
            for (const modeName of modeNames) {
                const result = await this._solve_single_mode(modeName, true, false);
                modeResults.push(result);
            }

            // Compute max errors across all modes
            let max_energy_err = 0;
            let max_param_err = 0;

            for (let i = 0; i < modeNames.length; i++) {
                const modeName = modeNames[i];
                const r = modeResults[i];

                const { energy, rel_error: energy_err } = this._compute_energy_error(r.Ex, r.Ey, prevEnergy[modeName]);
                // The vacuum field is scored too. It is not a duplicate of the
                // dielectric one, it carries C0 (hence Z0) and it is the entire
                // conductor-loss integrand (see _mode_conductor_loss), which is a
                // surface quantity that settles later than the dielectric energy.
                let energy0_err = 0;
                if (r.Ex0 && r.Ey0) {
                    const e0 = this._compute_energy_error(r.Ex0, r.Ey0, prevEnergy0[modeName], true);
                    energy0_err = e0.rel_error;
                    prevEnergy0[modeName] = e0.energy;
                }

                const param_err = prevZ0[modeName] !== undefined
                    ? Math.abs(r.Z0 - prevZ0[modeName]) / Math.max(Math.abs(prevZ0[modeName]), 1e-12)
                    : 1.0;

                max_energy_err = Math.max(max_energy_err, energy_err, energy0_err);
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
            // important for any mode (e.g. even mode near ground planes) get
            // refined. Each mode contributes both its dielectric solve and its
            // vacuum (C0) solve: the certificate certifies C and C0.
            if (it !== max_iters - 1) {
                const refineSets = modeResults.flatMap(m => {
                    const planeBC = this._plane_bc(m.mode);
                    const sets = [{ V: m.V, Ex: m.Ex, Ey: m.Ey, planeBC }];
                    if (m.V0 && m.Ex0 && m.Ey0) {
                        sets.push({ V: m.V0, Ex: m.Ex0, Ey: m.Ey0, vacuum: true, planeBC });
                    }
                    return sets;
                });
                this.refine_mesh_multi(refineSets, refineFrac);
                this._repaint_geometry();
            }
        }

        for (let i = 0; i < modeResults.length; i++) modeResults[i] = await this._mode_loss_results(modeResults[i]);

        // V0 is only needed by the refinement metrics above. Ex0/Ey0 are what
        // the conductor-loss calculation reuses across the frequency sweep.
        // Drop the vacuum potential so a converged mode result to save some
        // memory. (~5 MB/mode on a 600k-node grid).
        for (const m of modeResults) m.V0 = null;

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
                this._certWarn = { type: 'accuracy', reason: 'certificate', mode: 'all', message:
                    `Quasi-static mesh refinement could not verify the requested tolerance ` +
                    `(${pct(energy_tol)}%): estimated remaining error is ${estStr}. ` +
                    `Increase Max Nodes / Max Iterations, or relax Tolerance.` };
            }
        }

        // For an ASYMMETRIC differential pair, replace the odd/even drive results with the
        // genuine two-conductor modal decomposition on the converged mesh. _solve_modal_differential
        // returns null for a symmetric, degenerate or half-domain pair, in which case the
        // odd/even results (already in modeResults) are exact and are kept.
        if (this.is_differential) {
            const modal = await this._solve_modal_differential();
            if (modal) modeResults = modal;
        }

        // Store fields as arrays
        this.V = modeResults.map(r => r.V);
        this.Ex = modeResults.map(r => r.Ex);
        this.Ey = modeResults.map(r => r.Ey);
        this._plotModes = modeResults.map(r => r.mode);

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

        // Failed verification certificate and/or DC-skin transition note: surface as
        // accuracy warnings on every result built from this grid, the same per-solve
        // channel the triangular backend uses (result.warnings / solver.modeWarnings),
        // so the UI logs them for the adaptive solve and for every sweep point alike.
        const warns = [];
        if (this._certWarn) warns.push(this._certWarn);
        if (this._skinTransitionWarn) warns.push(this._skinTransitionWarn);
        if (this._platingTransitionWarn) warns.push(this._platingTransitionWarn);
        // Geometry-level accuracy note a subclass may set at construction (e.g. the
        // broadside strong-coupling warning). Like the two above, this only reaches
        // rectilinear results: the triangular backend returns from solveAt before
        // _build_results and models these regimes accurately (MQS).
        if (this._proximityWarn) warns.push(this._proximityWarn);
        if (this._causalWarn) warns.push(this._causalWarn);
        const openWarn = this.openBoundaryFieldWarning(modeResults.map(m => m.V));
        if (openWarn) warns.push(openWarn);
        if (warns.length) {
            result.warnings = warns;
            this.modeWarnings = result.warnings;
        }

        if (this.is_differential) {
            const odd = modeResults.find(m => m.mode === 'odd');
            const even = modeResults.find(m => m.mode === 'even');
            result.Z_diff = 2 * odd.Z0;
            result.Z_common = even.Z0 / 2;

            // Traces of different metal or finish: R11 - R22 and L11 - L22, carried on
            // the mode RLGC (see buildPhysicalRLGC).
            const asym = this._line_asymmetry(odd, even);
            if (asym) for (const m of [odd, even]) Object.assign(m.RLGC, asym);
            // Asymmetric geometry: the exact dielectric loss matrix of the per-trace solves.
            if (this._modalPhys && this._traceVac && this._traceVac.Gw) {
                const w = 2 * Math.PI * this.freq;
                for (const m of [odd, even]) m.RLGC.Gm = this._traceVac.Gw.map(v => v * w);
            }

            // Add physical 2x2 RLGC matrix
            result.RLGC_matrix = this._modal_to_physical_rlgc(odd, even);
            // True physical 2×2 [C]/[L] for the asymmetric MTL 4-port S-parameter path.
            if (this._modalPhys) result.physMatrix = this._modalPhys;
        }

        return result;
    }

    // Plotting payload: the grid and per-mode fields the app's field view reads.
    // On a half-domain (sym_half) solve the solver-internal arrays cover only
    // x >= 0. This mirrors fresh copies onto the full domain with the correct
    // parity per mode so plots (and the mesh-line overlay, which draws the same
    // x array) span the full cross-section. The internal half-domain arrays are
    // never touched, the sweep paths keep recomputing losses from them.
    //
    //   V(-x) = sV * V(x) with sV = -1 for 'odd' (electric wall) else +1
    //   Ex(-x) = -sV * Ex(x). Ey(-x) = sV * Ey(x).
    getPlotFields() {
        const base = { x: this.x, y: this.y, V: this.V, Ex: this.Ex, Ey: this.Ey,
                       triMesh: this.triMesh || null };
        // Mirror only a genuine half grid (x[0] === 0). Fields grafted by the
        // triangular backend already span the full domain (x[0] < 0).
        if (!this.sym_half || !this.V || !this.x || this.x.length < 2 || this.x[0] < 0)
            return base;
        const n = this.x.length;
        const nf = 2 * n - 1;
        const xs = new Float64Array(nf);
        for (let j = 0; j < n - 1; j++) xs[j] = -this.x[n - 1 - j];
        for (let j = 0; j < n; j++) xs[n - 1 + j] = this.x[j];
        const mirror = (A, sign) => A.map(row => {
            const out = new Float64Array(nf);
            for (let j = 0; j < n - 1; j++) out[j] = sign * row[n - 1 - j];
            for (let j = 0; j < n; j++) out[n - 1 + j] = row[j];
            return out;
        });
        const V = [], Ex = [], Ey = [];
        for (let m = 0; m < this.V.length; m++) {
            const sV = (this._plotModes && this._plotModes[m] === 'odd') ? -1 : 1;
            V.push(mirror(this.V[m], sV));
            Ex.push(mirror(this.Ex[m], -sV));
            Ey.push(mirror(this.Ey[m], sV));
        }
        return { x: xs, y: this.y, V, Ex, Ey, triMesh: this.triMesh || null };
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

        // DC internal inductance of the traces, cached per grid.
        for (const mode of this.is_differential ? ['odd', 'even'] : ['single']) {
            await this._ensure_dc_signal_inductance(this._plane_bc(mode));
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
                        mode.Ex, mode.Ey, mode.Z0, mode.C0, mode.Ex0, mode.Ey0, mode.mode);
                    const alpha_d = this.calculate_dielectric_loss(mode.V, mode.Z0);
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
                    result.Ex, result.Ey, result.Z0, result.C0, result.Ex0, result.Ey0, result.mode);
                const alpha_d = this.calculate_dielectric_loss(result.V, result.Z0);
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
            const { R_ac, R_dc, R_total, L_internal } = this._mode_conductor_loss(Ex, Ey, Z0, C0, Ex0, Ey0, mode);

            // Recalculate dielectric loss (affects omega)
            const alpha_d = this.calculate_dielectric_loss(V, Z0);

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
