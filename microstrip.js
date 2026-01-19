import { FieldSolver2D, CONSTANTS, diff } from './field_solver.js';

class MicrostripSolver extends FieldSolver2D {
    constructor(w, h, t, er, tan_delta, sigma_cond, freq, nx, ny) { // Renamed params for clarity
        super();
        this.w = w;
        this.h = h;
        this.t = t;
        this.er = er;
        this.tan_delta = tan_delta; // Set for FieldSolver2D methods
        this.sigma_cond = sigma_cond; // Set for FieldSolver2D methods
        this.freq = freq;
        this.omega = 2 * Math.PI * freq; // Renamed to 'omega' for FieldSolver2D methods
        
        // Requirements
        this.nx_req = nx;
        this.ny_req = ny;

        // Derived geometry
        this.domain_width = 2 * Math.max(this.w * 8, this.h * 15);
        this.top_air = this.h * 15; // Matched Python 'h * 15' default

        // Coordinates (Y)
        this.y_gnd_bot_start = 0;
        this.y_gnd_bot_end = 0; // Assuming thin ground for this model unless t_gnd specified
        this.y_sub_start = 0;
        this.y_sub_end = this.h;
        this.y_trace_start = this.h;
        this.y_trace_end = this.h + this.t;
        this.y_top_start = this.y_trace_end;
        this.y_top_end = this.y_top_start + this.top_air;

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
        const ncorner = Math.max(3, Math.floor(n / 40));

        // Critical X Interfaces
        let x_if = [0, xl, xr, this.domain_width];
        // Sort and Unique
        x_if = Float64Array.from(new Set(x_if)).sort();

        const n_regions = x_if.length - 1;

        // 1. Calculate Weights
        let region_weights = [];
        for (let k = 0; k < n_regions; k++) {
            const x0 = x_if[k];
            const x1 = x_if[k + 1];
            const width = x1 - x0;
            let weight = 1.0;

            // Trace Region check (using epsilon tolerance)
            if (x0 >= xl - 1e-15 && x1 <= xr + 1e-15) {
                weight = 6.0; // High weight for trace
            }
            region_weights.push(weight * width);
        }

        // 2. Allocate Points
        const total_weight = region_weights.reduce((a, b) => a + b, 0);
        let region_points = [];
        let allocated = 0;

        for (let k = 0; k < n_regions; k++) {
            let pts;
            if (k === n_regions - 1) {
                pts = n - allocated;
            } else {
                pts = Math.max(3, Math.floor(n * region_weights[k] / total_weight));
            }
            region_points.push(pts);
            allocated += pts;
        }

        // 3. Generate Segments
        let x_parts = [];
        for (let k = 0; k < n_regions; k++) {
            const x0 = x_if[k];
            const x1 = x_if[k + 1];
            const npts = region_points[k];
            let seg;

            if (x0 >= xl - 1e-15 && x1 <= xr + 1e-15) {
                // --- TRACE REGION (Complex Grading) ---
                const nc = Math.min(ncorner, Math.floor(npts / 3));
                const n_mid = Math.max(3, npts - 2 * nc);

                const corner_left = Math.min(corner, (x1 - x0) / 3);
                const corner_right = Math.min(corner, (x1 - x0) / 3);

                let cl, cr;

                // Left Corner
                if (nc > 0 && corner_left > 1e-15) {
                    cl = this._smooth_transition(x0, x0 + corner_left, nc, 'start', 4.0);
                } else {
                    cl = new Float64Array([x0]);
                }

                // Middle
                const mid0 = x0 + corner_left;
                const mid1 = Math.max(x1 - corner_right, mid0);
                const mid = this._smooth_transition(mid0, mid1, n_mid, 'both', 4.0);

                // Right Corner
                if (nc > 0 && corner_right > 1e-15) {
                    cr = this._smooth_transition(mid1, x1, nc, 'end', 4.0);
                } else {
                    cr = new Float64Array([x1]);
                }

                // Concatenate Trace parts: cl + mid[1:] + cr[1:]
                seg = this._concat_arrays([cl, mid.subarray(1), cr.subarray(1)]);

            } else {
                // --- AIR / DIELECTRIC REGION ---
                let end_curve = "both";
                let beta_val = 3.0;

                // Edges of Domain
                if (Math.abs(x0) < 1e-15) {
                    end_curve = "end"; // Dense right
                } else if (Math.abs(x1 - this.domain_width) < 1e-15) {
                    end_curve = "start"; // Dense left
                }

                // Edges of Trace
                if (x1 <= xl + 1e-15 || x0 >= xr - 1e-15) {
                    if (x1 <= xl + 1e-15) end_curve = "end";
                    else end_curve = "start";
                    beta_val = 4.0;
                }

                seg = this._smooth_transition(x0, x1, npts, end_curve, beta_val);
            }

            // Stitching: Skip first point for segments after the first one
            if (k > 0) {
                seg = seg.subarray(1);
            }
            x_parts.push(seg);
        }

        let x = this._concat_arrays(x_parts);

        // 4. Enforce Interfaces (Inject points if missed)
        x = this._enforce_interfaces(x, x_if);

        return x;
    }

    // --- Y Grid Generation Port ---
    _grid_y(n) {
        const ds = this.delta_s;
        const corner = Math.min(3 * ds, this.t / 4);
        const ncorner = Math.max(3, Math.floor(n / 40));

        // Critical Y Interfaces
        let y_if = [
            this.y_gnd_bot_start,
            this.y_gnd_bot_end, // likely 0 if t_gnd=0, unique will handle it
            this.y_sub_start,
            this.y_sub_end,
            this.y_trace_start,
            this.y_trace_end,
            this.y_top_start,
            this.y_top_end
        ];
        y_if = Float64Array.from(new Set(y_if)).sort();
        const n_regions = y_if.length - 1;

        // 1. Weights
        let region_weights = [];
        for (let k = 0; k < n_regions; k++) {
            const y0 = y_if[k];
            const y1 = y_if[k + 1];
            const height = y1 - y0;
            let weight = 0.1; // Default Air

            if (y0 >= this.y_trace_start - 1e-15 && y1 <= this.y_trace_end + 1e-15) {
                weight = 10.0; // Trace
            } else if (y1 <= this.y_gnd_bot_end + 1e-15) {
                weight = 0.0; // Ground (metal thickness)
            } else if (y0 >= this.y_sub_start - 1e-15 && y1 <= this.y_sub_end + 1e-15) {
                weight = 0.75; // Substrate
            }
            
            region_weights.push(weight * height);
        }

        // 2. Allocate
        const total_weight = region_weights.reduce((a, b) => a + b, 0);
        let region_points = [];
        let allocated = 0;
        
        for (let k = 0; k < n_regions; k++) {
            let pts;
            if (k === n_regions - 1) pts = n - allocated;
            else pts = Math.max(3, Math.floor(n * region_weights[k] / total_weight));
            region_points.push(pts);
            allocated += pts;
        }

        // 3. Segments
        let y_parts = [];
        for (let k = 0; k < n_regions; k++) {
            const y0 = y_if[k];
            const y1 = y_if[k + 1];
            const npts = region_points[k];
            let seg;
            let end_curve = "both";
            let beta_val = 3.0;

            if (Math.abs(y0) < 1e-15) {
                end_curve = "end";
                beta_val = 4.0;
            }

            if (y0 >= this.y_trace_start - 1e-15 && y1 <= this.y_trace_end + 1e-15) {
                // --- TRACE (Vertical) ---
                const nc = Math.min(ncorner, Math.floor(npts / 3));
                const n_mid = Math.max(4, npts - 2 * nc);
                
                const corner_bot = Math.min(corner, (y1 - y0)/3);
                const corner_top = Math.min(corner, (y1 - y0)/3);
                
                let cb, ct;

                // Bottom Corner
                if (nc > 0 && corner_bot > 1e-15) {
                    cb = this._smooth_transition(y0, y0 + corner_bot, nc, 'start', 4.0);
                } else {
                    cb = new Float64Array([y0]);
                }

                // Mid
                const mid0 = y0 + corner_bot;
                const mid1 = Math.max(y1 - corner_top, mid0);
                const mid = this._smooth_transition(mid0, mid1, n_mid, 'both', 2.0);

                // Top Corner
                if (nc > 0 && corner_top > 1e-15) {
                    ct = this._smooth_transition(mid1, y1, nc, 'end', 4.0);
                } else {
                    ct = new Float64Array([y1]);
                }

                seg = this._concat_arrays([cb, mid.subarray(1), ct.subarray(1)]);

            } else if (y0 >= this.y_sub_start - 1e-15 && y1 <= this.y_sub_end + 1e-15) {
                // Substrate
                seg = this._smooth_transition(y0, y1, npts, 'end', 3.5);
            } else {
                // Air / Others
                seg = this._smooth_transition(y0, y1, npts, end_curve, beta_val);
            }

            if (k > 0) seg = seg.subarray(1);
            y_parts.push(seg);
        }

        let y = this._concat_arrays(y_parts);
        y = this._enforce_interfaces(y, y_if);

        return y;
    }

    _smooth_transition(x0, x1, n, curve_end, beta) {
        // Standard geometric/hyperbolic grading logic
        // Matches typical Python solvers' behavior given parameters
        const pts = new Float64Array(n);
        for(let i=0; i<n; i++) {
            let u = i / (n - 1);
            let val;
            
            // Note: Different solvers use different map functions (sinh, tanh, power).
            // A robust "beta" based map often implies a Tanh or geometric progression.
            // Using Tanh map to support 'both' grading effectively.
            if(n <= 1) { pts[i] = x0; continue; }

            if (curve_end === 'start') {
                // Dense at start: u^beta
                val = Math.pow(u, beta);
            } else if (curve_end === 'end') {
                // Dense at end: 1 - (1-u)^beta
                val = 1.0 - Math.pow(1.0 - u, beta);
            } else {
                // Both: Tanh sigmoid map 
                // Map u [0,1] to [-1, 1], apply tanh, remap back
                const u_map = 2 * u - 1;
                const top = Math.tanh(beta * u_map);
                const bot = Math.tanh(beta);
                val = 0.5 * (1 + top / bot);
            }
            pts[i] = x0 + (x1 - x0) * val;
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
        const xl = cx - this.w/2;
        const xr = cx + this.w/2;

        // Iterate exactly as Python _setup_geometry
        for (let i = 0; i < ny; i++) {
            const yc = this.y[i];

            // 1. Permittivity (Row operation)
            // Python: if yc <= y_sub_end: er, else: er_top
            const is_sub = (yc - tol <= this.y_sub_end + tol);
            const val = is_sub ? this.er : (this.er_top || 1.0); // Handle er_top if exists
            
            for(let j=0; j<nx; j++) {
                this.epsilon_r[i][j] = val;
            }

            // 2. Bottom Ground
            if (yc >= this.y_gnd_bot_start - tol && yc <= this.y_gnd_bot_end + tol) {
                for(let j=0; j<nx; j++) this.ground_mask[i][j] = 1;
            }

            // 3. Trace
            if (yc >= this.y_trace_start - tol && yc <= this.y_trace_end + tol) {
                for(let j=0; j<nx; j++) {
                    const xc = this.x[j];
                    if (xc >= xl - tol && xc <= xr + tol) {
                        this.signal_mask[i][j] = 1;
                    }
                }
            }
        }

        // Finalize Potentials
        for(let i=0; i<ny; i++) {
            for(let j=0; j<nx; j++) {
                if(this.signal_mask[i][j]) {
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
