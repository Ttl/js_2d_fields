import { FieldSolver2D, CONSTANTS, diff } from './field_solver.js';

class MicrostripSolver extends FieldSolver2D {
    constructor(options) {
        super();

        // Required parameters
        this.h = options.substrate_height;
        this.w = options.trace_width;
        this.t = options.trace_thickness;
        this.er = options.epsilon_r;

        // Optional parameters with defaults
        this.t_gnd = options.gnd_thickness ?? 35e-6;
        this.tan_delta = options.tan_delta ?? 0.02;
        this.sigma_cond = options.sigma_cond ?? 5.8e7;
        this.freq = options.freq ?? 1e9;
        this.er_top = options.epsilon_r_top ?? 1;

        // Grid parameters
        this.nx_req = options.nx ?? 300;
        this.ny_req = options.ny ?? 300;
        this.skin_cells = options.skin_cells ?? 100;

        // Ground cut parameters
        this.gnd_cut_width = options.gnd_cut_width ?? 0.0;
        this.gnd_cut_sub_h = options.gnd_cut_sub_h ?? 0.0;

        // Top dielectric (embedded microstrip)
        this.top_diel_h = options.top_diel_h ?? 0.0;
        this.top_diel_er = options.top_diel_er ?? 1.0;
        this.top_diel_tand = options.top_diel_tand ?? 0.0;

        // Solder mask parameters
        this.use_sm = options.use_sm ?? false;
        this.sm_t_sub = options.sm_t_sub ?? 20e-6;
        this.sm_t_trace = options.sm_t_trace ?? 20e-6;
        this.sm_t_side = options.sm_t_side ?? 20e-6;
        this.sm_er = options.sm_er ?? 3.5;
        this.sm_tand = options.sm_tand ?? 0.02;

        // Domain and boundary parameters
        const air_side = options.air_side ?? null;
        const air_top = options.air_top ?? null;
        this.boundaries = options.boundaries ?? ["open", "open", "open", "gnd"];

        this.omega = 2 * Math.PI * this.freq;

        // Domain width calculation
        if (air_side === null) {
            this.domain_width = 2 * Math.max(this.w * 8, this.h * 15);
        } else {
            this.domain_width = this.w + 2 * air_side;
        }

        // --- Physical Y-Coordinates ---
        // Bottom extension for cut ground
        this.y_ext_start = this.t_gnd;
        this.y_ext_end = this.t_gnd + this.gnd_cut_sub_h;

        // New bottom ground plane location
        this.y_gnd_bot_start = this.y_ext_end;
        this.y_gnd_bot_end = this.y_gnd_bot_start + this.t_gnd;
        if (this.gnd_cut_width === 0) {
            this.y_gnd_bot_end = this.y_gnd_bot_start;
        }

        this.y_sub_start = this.y_gnd_bot_end;
        this.y_sub_end = this.y_sub_start + this.h;

        // --- TOP DIELECTRIC ---
        this.y_top_diel_start = this.y_sub_end;
        this.y_top_diel_end = this.y_top_diel_start + this.top_diel_h;

        // Trace is embedded in top dielectric
        this.y_trace_start = this.y_top_diel_start;
        this.y_trace_end = this.y_trace_start + this.t;

        // Solder mask extents
        this.y_sm_sub_end = this.y_top_diel_end + this.sm_t_sub;
        this.y_sm_trace_end = this.y_trace_end + this.sm_t_trace;

        this.y_top_start = this.y_top_diel_end;
        if (this.use_sm) {
            this.y_top_start = Math.max(this.y_sm_sub_end, this.y_sm_trace_end);
        }

        if (air_top === null) {
            this.top_dielectric_h = this.h * 15;
            this.has_top_gnd = false;
        } else {
            this.top_dielectric_h = air_top + this.t;
            this.has_top_gnd = (this.boundaries[2] === "gnd");
        }

        this.y_top_end = this.y_top_start + this.top_dielectric_h;

        if (this.has_top_gnd) {
            this.y_gnd_top_start = this.y_top_end;
            this.y_gnd_top_end = this.y_gnd_top_start + this.t_gnd;
            this.domain_height = this.y_gnd_top_end;
        } else {
            this.y_gnd_top_start = null;
            this.y_gnd_top_end = null;
            this.domain_height = this.y_top_end;
        }

        this.generate_grid();
    }



    generate_grid() {
        this.delta_s = Math.sqrt(2 / (this.omega * CONSTANTS.MU0 * this.sigma_cond)); // Use 'omega' and 'sigma_cond'

        // Generate Grids
        this.x = this._grid_x(this.nx_req);
        this.y = this._grid_y(this.ny_req);

        // Compute differentials (dx, dy) for solver use
        this.dx = new Float64Array(this.x.length - 1);
        for(let i=0; i<this.x.length-1; i++) this.dx[i] = this.x[i+1] - this.x[i];

        this.dy = new Float64Array(this.y.length - 1);
        for(let i=0; i<this.y.length-1; i++) this.dy[i] = this.y[i+1] - this.y[i];

        this.init_matrices();
    }


    // --- X Grid Generation Port ---
    _grid_x(n) {
        const cx = this.domain_width / 2;
        const xl = cx - this.w / 2;
        const xr = cx + this.w / 2;

        const ds = this.delta_s;
        const corner = Math.min(3 * ds, this.w / 4);
        const ncorner = Math.max(3, Math.floor(this.skin_cells / 4));

        // Collect all critical x-interfaces
        let x_if = [0, xl, xr, this.domain_width];

        if (this.gnd_cut_width > 0) {
            const cut_l = cx - this.gnd_cut_width / 2;
            const cut_r = cx + this.gnd_cut_width / 2;
            x_if.push(cut_l, cut_r);
        }

        if (this.use_sm) {
            const xsl = xl - this.sm_t_side;
            const xsr = xr + this.sm_t_side;
            x_if.push(xsl, xsr);
        }

        // Sort and unique
        x_if = Array.from(new Set(x_if)).sort((a, b) => a - b);
        const n_regions = x_if.length - 1;

        // Adaptive point allocation based on field importance
        let region_weights = [];
        for (let k = 0; k < n_regions; k++) {
            const x0 = x_if[k];
            const x1 = x_if[k + 1];
            const width = x1 - x0;
            let weight;

            // High weight for trace region and near-trace
            if (x0 >= xl - 1e-15 && x1 <= xr + 1e-15) {
                weight = 6.0;  // Trace itself - highest field
            } else {
                weight = 1.0;
            }

            // Ground cut gets extra points (field discontinuity)
            if (this.gnd_cut_width > 0) {
                const cut_l = cx - this.gnd_cut_width / 2;
                const cut_r = cx + this.gnd_cut_width / 2;
                if (Math.abs(x0 - cut_l) < 1e-15 || Math.abs(x1 - cut_r) < 1e-15) {
                    weight *= 1.5;
                }
            }

            region_weights.push(weight * width);
        }

        // Normalize and allocate points
        const total_weight = region_weights.reduce((a, b) => a + b, 0);
        let region_points = [];
        let allocated = 0;

        for (let k = 0; k < n_regions; k++) {
            let pts;
            if (k === n_regions - 1) {
                pts = Math.max(3, n - allocated);
            } else {
                pts = Math.max(3, Math.floor(n * region_weights[k] / total_weight));
            }
            region_points.push(pts);
            allocated += pts;
        }

        // Generate grid segments
        let x_parts = [];

        for (let k = 0; k < n_regions; k++) {
            const x0 = x_if[k];
            const x1 = x_if[k + 1];
            const npts = region_points[k];
            let seg;

            if (x0 >= xl - 1e-15 && x1 <= xr + 1e-15) {
                // TRACE REGION - fine mesh with corner grading
                const nc = Math.min(ncorner, Math.floor(npts / 3));
                const n_mid = Math.max(3, npts - 2 * nc);

                const corner_left = Math.min(corner, (x1 - x0) / 3);
                const corner_right = Math.min(corner, (x1 - x0) / 3);

                let cl, cr;

                if (nc > 0 && corner_left > 1e-15) {
                    cl = this._smooth_transition(x0, x0 + corner_left, nc, 'start', 4.0);
                } else {
                    cl = new Float64Array([x0]);
                }

                const mid0 = x0 + corner_left;
                const mid1 = Math.max(x1 - corner_right, mid0);
                const mid = this._smooth_transition(mid0, mid1, n_mid, 'both', 4.0);

                if (nc > 0 && corner_right > 1e-15) {
                    cr = this._smooth_transition(mid1, x1, nc, 'end', 4.0);
                } else {
                    cr = new Float64Array([x1]);
                }

                seg = this._concat_arrays([cl, mid.subarray(1), cr.subarray(1)]);

            } else {
                // DIELECTRIC/AIR REGION - grade toward interfaces
                let end_curve = "both";
                let beta_val = 3.0;

                // Grade AWAY from domain edges (fewer points at edges)
                if (Math.abs(x0) < 1e-15) {
                    end_curve = "end";  // Dense at right (toward trace)
                    beta_val = 3.0;
                } else if (Math.abs(x1 - this.domain_width) < 1e-15) {
                    end_curve = "start";  // Dense at left (toward trace)
                    beta_val = 3.0;
                }

                // Grade toward trace edges (more points near trace)
                if (x1 <= xl + 1e-15 || x0 >= xr - 1e-15) {
                    if (x1 <= xl + 1e-15) {
                        end_curve = "end";  // Dense toward trace (right side)
                    } else {
                        end_curve = "start";  // Dense toward trace (left side)
                    }
                    beta_val = 4.0;
                }

                seg = this._smooth_transition(x0, x1, npts, end_curve, beta_val);
            }

            if (k > 0) {
                seg = seg.subarray(1);
            }
            x_parts.push(seg);
        }

        let x = this._concat_arrays(x_parts);

        // Enforce exact interface locations
        for (const xi of x_if) {
            let min_dist = Infinity;
            for (let i = 0; i < x.length; i++) {
                min_dist = Math.min(min_dist, Math.abs(x[i] - xi));
            }
            if (min_dist > 1e-12) {
                const newX = new Float64Array(x.length + 1);
                newX.set(x);
                newX[x.length] = xi;
                x = newX.sort();
            }
        }

        return x;
    }

    // --- Y Grid Generation Port ---
    _grid_y(n) {
        const ds = this.delta_s;
        const corner = Math.min(3 * ds, this.t / 4);
        const ncorner = Math.max(3, Math.floor(this.skin_cells / 4));

        // Collect all critical y-interfaces
        let y_if = [
            0.0,
            this.y_gnd_bot_start,
            this.y_gnd_bot_end,
            this.y_sub_start,
            this.y_sub_end,
        ];

        if (this.top_diel_h > 0) {
            y_if.push(this.y_top_diel_start, this.y_top_diel_end);
        }

        y_if.push(
            this.y_trace_start,
            this.y_trace_end,
            this.y_top_start,
            this.y_top_end
        );

        if (this.use_sm) {
            y_if.push(this.y_sm_sub_end, this.y_sm_trace_end);
        }

        if (this.has_top_gnd) {
            y_if.push(this.y_gnd_top_start, this.y_gnd_top_end);
        }

        // Sort and unique
        y_if = Array.from(new Set(y_if)).sort((a, b) => a - b);
        const n_regions = y_if.length - 1;

        // Adaptive point allocation based on field importance
        let region_weights = [];
        for (let k = 0; k < n_regions; k++) {
            const y0 = y_if[k];
            const y1 = y_if[k + 1];
            const height = y1 - y0;
            let weight;

            // TRACE - highest field concentration
            if (y0 >= this.y_trace_start - 1e-15 && y1 <= this.y_trace_end + 1e-15) {
                weight = 10;
            }
            // GROUND PLANES - high current density
            else if (y1 <= this.y_ext_start + 1e-15 ||
                     (this.has_top_gnd &&
                      y0 >= this.y_gnd_top_start - 1e-15 &&
                      y1 <= this.y_gnd_top_end + 1e-15)) {
                weight = 0;
            }
            // SUBSTRATE - high field region
            else if (y0 >= this.y_sub_start - 1e-15 && y1 <= this.y_sub_end + 1e-15) {
                weight = 0.75;
            }
            // TOP DIELECTRIC - field transitions
            else if (this.top_diel_h > 0 &&
                     y0 >= this.y_top_diel_start - 1e-15 &&
                     y1 <= this.y_top_diel_end + 1e-15) {
                weight = 1;
            }
            // SOLDER MASK - moderate field
            else if (this.use_sm && (
                (y0 >= this.y_top_diel_end - 1e-15 && y1 <= this.y_sm_sub_end + 1e-15) ||
                (y0 >= this.y_trace_end - 1e-15 && y1 <= this.y_sm_trace_end + 1e-15))) {
                weight = 1;
            }
            // AIR/FAR-FIELD - lower field
            else {
                if (this.has_top_gnd) {
                    weight = 0.75;
                } else {
                    weight = 0.1;
                }
            }

            region_weights.push(weight * height);
        }

        // Normalize and allocate points
        const total_weight = region_weights.reduce((a, b) => a + b, 0);
        let region_points = [];
        let allocated = 0;

        for (let k = 0; k < n_regions; k++) {
            let pts;
            if (k === n_regions - 1) {
                pts = Math.max(3, n - allocated);
            } else {
                pts = Math.max(3, Math.floor(n * region_weights[k] / total_weight));
            }
            region_points.push(pts);
            allocated += pts;
        }

        // Generate grid segments
        let y_parts = [];

        for (let k = 0; k < n_regions; k++) {
            const y0 = y_if[k];
            const y1 = y_if[k + 1];
            const npts = region_points[k];
            let seg;

            // Determine grading strategy
            let end_curve = "both";
            let beta_val = 3.0;

            // Domain boundaries - grade AWAY from edges
            if (Math.abs(y0) < 1e-15) {
                end_curve = "end";  // Dense at top (away from bottom edge)
                beta_val = 4;
            }
            if (this.has_top_gnd && Math.abs(y1 - this.y_gnd_top_end) < 1e-15) {
                end_curve = "start";  // Dense at bottom (away from top edge)
                beta_val = 4;
            } else if (!this.has_top_gnd && Math.abs(y1 - this.y_top_end) < 1e-15) {
                end_curve = "start";  // Dense at bottom (away from top edge)
                beta_val = 4;
            }

            // TRACE REGION - special corner grading
            if (y0 >= this.y_trace_start - 1e-15 && y1 <= this.y_trace_end + 1e-15) {
                const nc = Math.min(ncorner, Math.floor(npts / 3));
                const n_mid = Math.max(4, npts - 2 * nc);

                const corner_bot = Math.min(corner, (y1 - y0) / 3);
                const corner_top = Math.min(corner, (y1 - y0) / 3);

                let cb, ct;

                if (nc > 0 && corner_bot > 1e-15) {
                    cb = this._smooth_transition(y0, y0 + corner_bot, nc, 'start', 4.0);
                } else {
                    cb = new Float64Array([y0]);
                }

                const mid0 = y0 + corner_bot;
                const mid1 = Math.max(y1 - corner_top, mid0);
                const mid = this._smooth_transition(mid0, mid1, n_mid, 'both', 2.0);

                if (nc > 0 && corner_top > 1e-15) {
                    ct = this._smooth_transition(mid1, y1, nc, 'end', 4.0);
                } else {
                    ct = new Float64Array([y1]);
                }

                seg = this._concat_arrays([cb, mid.subarray(1), ct.subarray(1)]);
            }
            // GROUND PLANES - skin depth grading
            else if ((y0 >= this.y_gnd_bot_start - 1e-15 && y1 <= this.y_gnd_bot_end + 1e-15) ||
                     (this.has_top_gnd &&
                      y0 >= this.y_gnd_top_start - 1e-15 &&
                      y1 <= this.y_gnd_top_end + 1e-15)) {
                seg = this._smooth_transition(y0, y1, npts, end_curve, 4.0);
            }
            // SUBSTRATE & TOP DIELECTRIC - grade toward trace interface
            else if ((y0 >= this.y_sub_start - 1e-15 && y1 <= this.y_sub_end + 1e-15) ||
                     (this.top_diel_h > 0 &&
                      y0 >= this.y_top_diel_start - 1e-15 &&
                      y1 <= this.y_top_diel_end + 1e-15)) {
                // Grade toward trace (more points near trace)
                if (y0 >= this.y_top_diel_start - 1e-15 && y1 <= this.y_top_diel_end + 1e-15) {
                    // Top dielectric - grade toward trace at bottom
                    seg = this._smooth_transition(y0, y1, npts, 'start', 3.5);
                } else {
                    // Substrate - grade toward trace at top
                    seg = this._smooth_transition(y0, y1, npts, 'end', 3.5);
                }
            }
            // AIR/FAR-FIELD - less aggressive grading
            else {
                seg = this._smooth_transition(y0, y1, npts, end_curve, beta_val);
            }

            if (k > 0) {
                seg = seg.subarray(1);
            }
            y_parts.push(seg);
        }

        let y = this._concat_arrays(y_parts);

        // Enforce exact interface locations
        for (const yi of y_if) {
            let min_dist = Infinity;
            for (let i = 0; i < y.length; i++) {
                min_dist = Math.min(min_dist, Math.abs(y[i] - yi));
            }
            if (min_dist > 1e-12) {
                const newY = new Float64Array(y.length + 1);
                newY.set(y);
                newY[y.length] = yi;
                y = newY.sort();
            }
        }

        return y;
    }

    _smooth_transition(x0, x1, n, curve_end, beta) {
        // Matches Python field_solver_ref.py _smooth_transition exactly
        if (n <= 1) {
            return new Float64Array([x0, x1]);
        }

        const pts = new Float64Array(n);
        const tanhBeta = Math.tanh(beta);
        const tanhBetaHalf = Math.tanh(beta * 0.5);

        for (let i = 0; i < n; i++) {
            const xi = i / (n - 1);
            let eta;

            if (curve_end === 'end') {
                // Dense at end: tanh(beta * xi) / tanh(beta)
                eta = Math.tanh(beta * xi) / tanhBeta;
            } else if (curve_end === 'both') {
                // Dense at both ends
                eta = (Math.tanh(beta * (xi - 0.5)) / tanhBetaHalf + 1) / 2;
            } else {
                // 'start' - Dense at start: 1 - tanh(beta * (1 - xi)) / tanh(beta)
                eta = 1 - Math.tanh(beta * (1 - xi)) / tanhBeta;
            }

            pts[i] = x0 + (x1 - x0) * eta;
        }
        return pts;
    }

    _enforce_interfaces(arr, interfaces) {
        // If an interface point is missing (distance > tol), inject it
        const tol = 1e-12;
        let list = Array.from(arr);
        let changed = false;

        for (let val of interfaces) {
            let min_dist = Number.MAX_VALUE;
            for(let p of list) min_dist = Math.min(min_dist, Math.abs(p - val));
            
            if (min_dist > tol) {
                list.push(val);
                changed = true;
            }
        }
        if(changed) {
            return Float64Array.from(new Set(list)).sort();
        }
        return arr;
    }

    _concat_arrays(arrays) {
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

    init_matrices() {
        const nx = this.x.length;
        const ny = this.y.length;

        // Initialize 2D Arrays (in Row-Major: [ny][nx])
        this.epsilon_r = Array(ny).fill().map(() => new Float64Array(nx).fill(1.0));
        this.conductor_mask = Array(ny).fill().map(() => new Uint8Array(nx).fill(0));
        this.signal_mask = Array(ny).fill().map(() => new Uint8Array(nx).fill(0));
        this.ground_mask = Array(ny).fill().map(() => new Uint8Array(nx).fill(0));
        this.V = Array(ny).fill().map(() => new Float64Array(nx).fill(0.0));

        const tol = 1e-11;
        const cx = this.domain_width / 2;
        const xl = cx - this.w / 2;
        const xr = cx + this.w / 2;
        const xsl = xl - this.sm_t_side;
        const xsr = xr + this.sm_t_side;

        const cut_l = cx - this.gnd_cut_width / 2;
        const cut_r = cx + this.gnd_cut_width / 2;

        for (let i = 0; i < ny; i++) {
            const yc = this.y[i];

            // --- PERMITTIVITY ---
            if (yc - tol <= this.y_sub_end + tol) {
                for (let j = 0; j < nx; j++) {
                    this.epsilon_r[i][j] = this.er;
                }
            } else if (yc - tol <= this.y_top_diel_end + tol) {
                for (let j = 0; j < nx; j++) {
                    this.epsilon_r[i][j] = this.top_diel_er;
                }
            } else {
                for (let j = 0; j < nx; j++) {
                    this.epsilon_r[i][j] = this.er_top;
                }

                if (this.use_sm) {
                    for (let j = 0; j < nx; j++) {
                        const xc = this.x[j];
                        if (yc <= this.y_sm_sub_end - tol) {
                            this.epsilon_r[i][j] = this.sm_er;
                        }
                        if ((xsl - tol <= xc && xc <= xsr + tol) && (yc - tol <= this.y_sm_trace_end + tol)) {
                            this.epsilon_r[i][j] = this.sm_er;
                        }
                    }
                }
            }

            // --- BOTTOM GROUND ---
            if (yc >= this.y_gnd_bot_start - tol && yc <= this.y_gnd_bot_end + tol) {
                for (let j = 0; j < nx; j++) {
                    const xc = this.x[j];
                    if (this.gnd_cut_width > 0 && cut_l <= xc && xc <= cut_r) {
                        continue;
                    }
                    this.ground_mask[i][j] = 1;
                }
            }

            if (yc <= this.y_ext_start + tol) {
                for (let j = 0; j < nx; j++) {
                    this.ground_mask[i][j] = 1;
                }
            }

            // --- TOP GROUND ---
            if (this.has_top_gnd && yc >= this.y_gnd_top_start - tol) {
                for (let j = 0; j < nx; j++) {
                    this.ground_mask[i][j] = 1;
                }
            }

            // --- TRACE ---
            if (this.y_trace_start - tol <= yc && yc <= this.y_trace_end + tol) {
                for (let j = 0; j < nx; j++) {
                    const xc = this.x[j];
                    if (xc >= xl && xc <= xr) {
                        this.signal_mask[i][j] = 1;
                    }
                }
            }
        }

        // Finalize potentials
        for (let i = 0; i < ny; i++) {
            for (let j = 0; j < nx; j++) {
                if (this.signal_mask[i][j]) {
                    this.V[i][j] = 1.0;
                    this.conductor_mask[i][j] = 1;
                } else if (this.ground_mask[i][j]) {
                    this.V[i][j] = 0.0;
                    this.conductor_mask[i][j] = 1;
                }
            }
        }
    }

}

export { MicrostripSolver };
