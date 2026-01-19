import { FieldSolver2D, CONSTANTS, diff } from './field_solver.js';

// Helper functions (copied from microstrip.js)
function _smooth_transition(x0, x1, n, curve_end, beta=4.0) {
    if (n === 0) {
        return new Float64Array(0);
    }
    if (n === 1) {
        return new Float64Array([x0]);
    }
    
    const pts = new Float64Array(n);

    for(let i=0; i<n; i++) {
        let u = i / (n - 1); // This will no longer be 0/0=NaN due to n > 1 check
        let val;
        
        if (curve_end === 'start') {
            val = Math.pow(u, beta);
        } else if (curve_end === 'end') {
            val = 1.0 - Math.pow(1.0 - u, beta);
        } else { // 'both'
            const u_map = 2 * u - 1;
            const top = Math.tanh(beta * u_map);
            const bot = Math.tanh(beta);
            val = 0.5 * (1 + top / bot);
        }
        pts[i] = x0 + (x1 - x0) * val;
    }
    return pts;
}

function _enforce_interfaces(arr, interfaces) {
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

function _concat_arrays(arrays) {
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


export class GroundedCPWSolver2D extends FieldSolver2D {
    constructor(
        substrate_height, trace_width, trace_thickness,
        gap, top_gnd_width, via_gap, via_diameter,
        gnd_thickness = 35e-6, epsilon_r = 4.5, tan_delta = 0.02,
        sigma_diel = 0.0, sigma_cond = 5.8e7, epsilon_r_top = 1,
        air_top = null, air_side = null, freq = 1e9,
        nx = 300, ny = 300, skin_cells = 50, boundaries = null,
        use_sm = false, sm_t_sub = 20e-6, sm_t_trace = 20e-6,
        sm_t_side = 20e-6, sm_er = 3.5, sm_tand = 0.02
    ) {
        super();

        this.h = substrate_height;
        this.w = trace_width;
        this.t = trace_thickness;
        this.gap = gap;
        this.top_gnd_w = top_gnd_width;
        this.via_gap = via_gap;
        this.via_d = via_diameter;
        this.t_gnd = gnd_thickness;
        this.er = epsilon_r;
        this.er_top = epsilon_r_top;
        this.tan_delta = tan_delta;
        this.sigma_diel = sigma_diel;
        this.sigma_cond = sigma_cond;

        // Solder Mask
        this.use_sm = use_sm;
        this.sm_t_sub = sm_t_sub;
        this.sm_t_trace = sm_t_trace;
        this.sm_t_side = sm_t_side;
        this.sm_er = sm_er;
        this.sm_tand = sm_tand;

        this.freq = freq;
        this.omega = 2 * CONSTANTS.PI * freq;
        this.nx_req = nx;
        this.ny_req = ny;
        this.skin_cells = skin_cells;

        this.active_width = this.w + 2 * Math.max(this.gap + this.top_gnd_w, this.via_gap + this.via_d);

        if (air_side === null) {
            this.domain_width = this.active_width * 1.5;
        } else {
            this.domain_width = this.active_width + 2 * air_side;
        }

        this.boundaries = boundaries || ["open", "open", "open", "gnd"];

        // --- Y-Coordinates ---
        this.y_gnd_bot_start = 0.0;
        this.y_gnd_bot_end   = this.t_gnd;
        this.y_sub_start     = this.y_gnd_bot_end;
        this.y_sub_end       = this.y_sub_start + this.h;
        this.y_trace_start   = this.y_sub_end;
        this.y_trace_end     = this.y_trace_start + this.t;
        this.y_sm_sub_end    = this.y_sub_end + this.sm_t_sub;
        this.y_sm_trace_end  = this.y_trace_end + this.sm_t_trace;
        this.y_top_start     = this.y_trace_end;
        if (use_sm) {
             this.y_top_start = Math.max(this.y_sm_sub_end, this.y_sm_trace_end);
        }
        
        this.top_dielectric_h = air_top !== null ? air_top : this.h * 10;
        this.y_top_end       = this.y_top_start + this.top_dielectric_h;
        this.has_top_gnd     = (this.boundaries[2] === "gnd");
        this.domain_height   = (this.y_top_end + this.t_gnd);
        if (!this.has_top_gnd) {
            this.domain_height = this.y_top_end;
        }

        if (this.has_top_gnd) {
            this.y_gnd_top_start = this.y_top_end;
            this.y_gnd_top_end   = this.y_gnd_top_start + this.t_gnd;
        } else {
            this.y_gnd_top_start = null;
            this.y_gnd_top_end = null;
        }

        // Skin depth
        this.delta_s = Math.sqrt(2 / (this.omega * CONSTANTS.MU0 * this.sigma_cond));

        this.generate_grid();
    }

    async calculate_parameters() {
        await this.solve_laplace_iterative(true);
        const C0 = this.calculate_capacitance(true);

        await this.solve_laplace_iterative(false);
        const C = this.calculate_capacitance(false);

        const eps_eff = C / C0;
        const Z0 = 1 / (CONSTANTS.C * Math.sqrt(C * C0));
        return { Z0, eps_eff, C, C0 };
    }

    generate_grid() {
        this._calculate_geometry_x();
        this.x = this._grid_x(this.nx_req);
        this.y = this._grid_y(this.ny_req);

        this.dx = diff(this.x);
        this.dy = diff(this.y);

        this.init_matrices();
    }

    // This needs to be implemented
    _calculate_geometry_x() {
        const cx = this.domain_width / 2;
        // Signal
        this.x_tr_l = cx - this.w / 2;
        this.x_tr_r = cx + this.w / 2;

        // Top Grounds
        this.x_gap_l = this.x_tr_l - this.gap;
        this.x_gap_r = this.x_tr_r + this.gap;

        // Vias (Single via on each side)
        // Inner edge of via is (trace_edge - via_gap)
        this.via_x_left_inner = this.x_tr_l - this.via_gap;
        this.via_x_left_center = this.via_x_left_inner - (this.via_d / 2);

        this.via_x_right_inner = this.x_tr_r + this.via_gap;
        this.via_x_right_center = this.via_x_right_inner + (this.via_d / 2);
    }

    // This needs to be implemented
    _grid_x(n) {
        const cx = this.domain_width / 2;

        const xl_trace = this.x_tr_l;
        const xr_trace = this.x_tr_r;
        const xl_gap = this.x_gap_l;
        const xr_gap = this.x_gap_r;

        const ds = this.delta_s;
        const corner_w = Math.min(3 * ds, this.w / 4);
        const corner_gap = Math.min(3 * ds, this.gap / 4);

        const sm = this.use_sm ? this.sm_t_side : 0.0;

        const n_corner = Math.max(5, Math.floor(n / 40));
        const n_outer_gnd = Math.max(Math.floor(n * 0.12), 8);
        const n_gap_bulk = Math.max(Math.floor(n * 0.14), 10);
        const n_trace_center = Math.max(Math.floor(n * 0.15), 8);

        // Solder mask sidewall resolution
        const n_sm = this.use_sm ? Math.max(5, Math.floor(n / 40)) : 0;

        // --- LEFT SIDE ---

        // 1. Left outer ground
        const a = _smooth_transition(0, xl_gap - corner_gap, n_outer_gnd, 'end');

        // 2. Left ground corner
        const a_corner = _smooth_transition(xl_gap - corner_gap, xl_gap, n_corner, 'end');

        let b_sm1 = null;
        let gap_start = xl_gap;
        // 2.5 Left gap SM sidewall (ground edge)
        if (this.use_sm && sm > 0) {
            b_sm1 = _smooth_transition(xl_gap, xl_gap + sm, n_sm, 'both', 4.0);
            gap_start = xl_gap + sm;
        }

        // 3. Left gap bulk
        const b = _smooth_transition(gap_start, xl_trace - sm, n_gap_bulk, 'both', 6.0);

        let b_sm2 = null;
        // 3.5 Left gap SM sidewall (trace edge)
        if (this.use_sm && sm > 0) {
            b_sm2 = _smooth_transition(xl_trace - sm, xl_trace, n_sm, 'both', 4.0);
        }

        // 4. Trace left corner
        const c_cl = _smooth_transition(xl_trace, xl_trace + corner_w, n_corner, 'start');

        // 5. Trace center
        const c = _smooth_transition(xl_trace + corner_w, xr_trace - corner_w, n_trace_center, 'both');

        // 6. Trace right corner
        const c_cr = _smooth_transition(xr_trace - corner_w, xr_trace, n_corner, 'end');

        // --- RIGHT SIDE ---

        let d_sm1 = null;
        let gap_r_start = xr_trace;
        // 7. Right gap SM sidewall (trace edge)
        if (this.use_sm && sm > 0) {
            d_sm1 = _smooth_transition(xr_trace, xr_trace + sm, n_sm, 'both', 4.0);
            gap_r_start = xr_trace + sm;
        }

        // 8. Right gap bulk
        const d = _smooth_transition(gap_r_start, xr_gap - sm, n_gap_bulk, 'both', 6.0);

        let d_sm2 = null;
        // 8.5 Right gap SM sidewall (ground edge)
        if (this.use_sm && sm > 0) {
            d_sm2 = _smooth_transition(xr_gap - sm, xr_gap, n_sm, 'both', 4.0);
        }

        // 9. Right ground corner
        const e_corner = _smooth_transition(xr_gap, xr_gap + corner_gap, n_corner, 'start');

        // Remaining outer ground
        let used_length = a.length + (a_corner.length > 1 ? a_corner.length - 1 : 0);
        if (b_sm1) used_length += (b_sm1.length > 1 ? b_sm1.length - 1 : 0);
        used_length += (b.length > 1 ? b.length - 1 : 0);
        if (b_sm2) used_length += (b_sm2.length > 1 ? b_sm2.length - 1 : 0);
        used_length += (c_cl.length > 1 ? c_cl.length - 1 : 0);
        used_length += (c.length > 1 ? c.length - 1 : 0);
        used_length += (c_cr.length > 1 ? c_cr.length - 1 : 0);
        if (d_sm1) used_length += (d_sm1.length > 1 ? d_sm1.length - 1 : 0);
        used_length += (d.length > 1 ? d.length - 1 : 0);
        if (d_sm2) used_length += (d_sm2.length > 1 ? d_sm2.length - 1 : 0);
        used_length += (e_corner.length > 1 ? e_corner.length - 1 : 0);
        
        const remaining_n = Math.max(n - used_length, 6);

        const e = _smooth_transition(xr_gap + corner_gap, this.domain_width, remaining_n + 1, 'start');

        // --- Combine all segments ---
        const parts = [a]; // Start with 'a' fully

        // Add subarrays, checking length before slicing
        if (a_corner.length > 1) parts.push(a_corner.subarray(1));
        if (b_sm1 && b_sm1.length > 1) parts.push(b_sm1.subarray(1));
        if (b.length > 1) parts.push(b.subarray(1));
        if (b_sm2 && b_sm2.length > 1) parts.push(b_sm2.subarray(1));

        if (c_cl.length > 1) parts.push(c_cl.subarray(1));
        if (c.length > 1) parts.push(c.subarray(1));
        if (c_cr.length > 1) parts.push(c_cr.subarray(1));

        if (d_sm1 && d_sm1.length > 1) parts.push(d_sm1.subarray(1));
        if (d.length > 1) parts.push(d.subarray(1));
        if (d_sm2 && d_sm2.length > 1) parts.push(d_sm2.subarray(1));

        if (e_corner.length > 1) parts.push(e_corner.subarray(1));
        if (e.length > 1) parts.push(e.subarray(1));


        let x = _concat_arrays(parts);
        x = _enforce_interfaces(x, [0, xl_trace, xr_trace, xl_gap, xr_gap, this.domain_width]); // Ensure critical points are there
        return x;
    }

    _grid_x(n) {
        const cx = this.domain_width / 2;

        const xl_trace = this.x_tr_l;
        const xr_trace = this.x_tr_r;
        const xl_gap = this.x_gap_l;
        const xr_gap = this.x_gap_r;

        const ds = this.delta_s;
        const corner_w = Math.min(3 * ds, this.w / 4);
        const corner_gap = Math.min(3 * ds, this.gap / 4);

        const sm = this.use_sm ? this.sm_t_side : 0.0;

        const n_corner = Math.max(5, Math.floor(n / 40));
        const n_outer_gnd = Math.max(Math.floor(n * 0.12), 8);
        const n_gap_bulk = Math.max(Math.floor(n * 0.14), 10);
        const n_trace_center = Math.max(Math.floor(n * 0.15), 8);

        // Solder mask sidewall resolution
        const n_sm = this.use_sm ? Math.max(5, Math.floor(n / 40)) : 0;

        // --- LEFT SIDE ---

        // 1. Left outer ground
        const a = _smooth_transition(0, xl_gap - corner_gap, n_outer_gnd, 'end');

        // 2. Left ground corner
        const a_corner = _smooth_transition(xl_gap - corner_gap, xl_gap, n_corner, 'end');

        let b_sm1 = null;
        let gap_start = xl_gap;
        // 2.5 Left gap SM sidewall (ground edge)
        if (this.use_sm && sm > 0) {
            b_sm1 = _smooth_transition(xl_gap, xl_gap + sm, n_sm, 'both', 4.0);
            gap_start = xl_gap + sm;
        }

        // 3. Left gap bulk
        const b = _smooth_transition(gap_start, xl_trace - sm, n_gap_bulk, 'both', 6.0);

        let b_sm2 = null;
        // 3.5 Left gap SM sidewall (trace edge)
        if (this.use_sm && sm > 0) {
            b_sm2 = _smooth_transition(xl_trace - sm, xl_trace, n_sm, 'both', 4.0);
        }

        // 4. Trace left corner
        const c_cl = _smooth_transition(xl_trace, xl_trace + corner_w, n_corner, 'start');

        // 5. Trace center
        const c = _smooth_transition(xl_trace + corner_w, xr_trace - corner_w, n_trace_center, 'both');

        // 6. Trace right corner
        const c_cr = _smooth_transition(xr_trace - corner_w, xr_trace, n_corner, 'end');

        // --- RIGHT SIDE ---

        let d_sm1 = null;
        let gap_r_start = xr_trace;
        // 7. Right gap SM sidewall (trace edge)
        if (this.use_sm && sm > 0) {
            d_sm1 = _smooth_transition(xr_trace, xr_trace + sm, n_sm, 'both', 4.0);
            gap_r_start = xr_trace + sm;
        }

        // 8. Right gap bulk
        const d = _smooth_transition(gap_r_start, xr_gap - sm, n_gap_bulk, 'both', 6.0);

        let d_sm2 = null;
        // 8.5 Right gap SM sidewall (ground edge)
        if (this.use_sm && sm > 0) {
            d_sm2 = _smooth_transition(xr_gap - sm, xr_gap, n_sm, 'both', 4.0);
        }

        // 9. Right ground corner
        const e_corner = _smooth_transition(xr_gap, xr_gap + corner_gap, n_corner, 'start');

        // Remaining outer ground
        let used_length = a.length + (a_corner.length > 1 ? a_corner.length - 1 : 0);
        if (b_sm1) used_length += (b_sm1.length > 1 ? b_sm1.length - 1 : 0);
        used_length += (b.length > 1 ? b.length - 1 : 0);
        if (b_sm2) used_length += (b_sm2.length > 1 ? b_sm2.length - 1 : 0);
        used_length += (c_cl.length > 1 ? c_cl.length - 1 : 0);
        used_length += (c.length > 1 ? c.length - 1 : 0);
        used_length += (c_cr.length > 1 ? c_cr.length - 1 : 0);
        if (d_sm1) used_length += (d_sm1.length > 1 ? d_sm1.length - 1 : 0);
        used_length += (d.length > 1 ? d.length - 1 : 0);
        if (d_sm2) used_length += (d_sm2.length > 1 ? d_sm2.length - 1 : 0);
        used_length += (e_corner.length > 1 ? e_corner.length - 1 : 0);
        
        const remaining_n = Math.max(n - used_length, 6);

        const e = _smooth_transition(xr_gap + corner_gap, this.domain_width, remaining_n + 1, 'start');

        // --- Combine all segments ---
        const parts = [a]; // Start with 'a' fully

        // Add subarrays, checking length before slicing
        if (a_corner.length > 1) parts.push(a_corner.subarray(1));
        if (b_sm1 && b_sm1.length > 1) parts.push(b_sm1.subarray(1));
        if (b.length > 1) parts.push(b.subarray(1));
        if (b_sm2 && b_sm2.length > 1) parts.push(b_sm2.subarray(1));

        if (c_cl.length > 1) parts.push(c_cl.subarray(1));
        if (c.length > 1) parts.push(c.subarray(1));
        if (c_cr.length > 1) parts.push(c_cr.subarray(1));

        if (d_sm1 && d_sm1.length > 1) parts.push(d_sm1.subarray(1));
        if (d.length > 1) parts.push(d.subarray(1));
        if (d_sm2 && d_sm2.length > 1) parts.push(d_sm2.subarray(1));

        if (e_corner.length > 1) parts.push(e_corner.subarray(1));
        if (e.length > 1) parts.push(e.subarray(1));


        let x = _concat_arrays(parts);
        x = _enforce_interfaces(x, [0, xl_trace, xr_trace, xl_gap, xr_gap, this.domain_width]); // Ensure critical points are there
        return x;
    }

    _grid_y(n) {
        const n_gnd = Math.max(3, Math.floor(n * 0.05));

        const ds = this.delta_s;
        const corner = Math.min(3 * ds, this.t / 4);
        const ncorner = Math.floor(this.skin_cells / 4);
        const n_trace = this.skin_cells + 4 - 2 * ncorner;
        const n_sm_val = 20; // This was hardcoded in python example for else branch.

        let grid_parts = [];
        let n_remain, n_sub_bot, n_sub_top;

        if (!this.use_sm) {
            n_remain = n - (this.has_top_gnd ? n_gnd * 2 : n_gnd) - n_trace;
            n_sub_bot = Math.floor(n_remain * 0.4);
            n_sub_top = n_remain - n_sub_bot;

            const l1 = _smooth_transition(this.y_gnd_bot_start, this.y_gnd_bot_end, n_gnd, 'both');
            const l2 = _smooth_transition(this.y_sub_start, this.y_sub_end, n_sub_bot, 'both');
            const l3_cb = _smooth_transition(this.y_trace_start, this.y_trace_start + corner, ncorner, 'start');
            const l3 = _smooth_transition(this.y_trace_start + corner, this.y_trace_end - corner, n_trace, 'both', 1);
            const l3_ct = _smooth_transition(this.y_trace_end - corner, this.y_trace_end, ncorner, 'end');
            const l4 = _smooth_transition(this.y_top_start, this.y_top_end, n_sub_top, 'start');
            
            grid_parts = [l1];
            if (l2.length > 1) grid_parts.push(l2.subarray(1));
            if (l3_cb.length > 1) grid_parts.push(l3_cb.subarray(1));
            if (l3.length > 1) grid_parts.push(l3.subarray(1));
            if (l3_ct.length > 1) grid_parts.push(l3_ct.subarray(1));
            if (l4.length > 1) grid_parts.push(l4.subarray(1));
        } else {
            n_remain = n - (this.has_top_gnd ? n_gnd * 2 : n_gnd) - n_trace - n_sm_val;
            n_sub_bot = Math.floor(n_remain * 0.4);
            n_sub_top = n_remain - n_sub_bot;

            const l1 = _smooth_transition(this.y_gnd_bot_start, this.y_gnd_bot_end, n_gnd, 'both');
            const l2 = _smooth_transition(this.y_sub_start, this.y_sub_end, n_sub_bot, 'both');
            const l3_cb = _smooth_transition(this.y_trace_start, this.y_trace_start + corner, ncorner, 'start');
            const l3 = _smooth_transition(this.y_trace_start + corner, this.y_trace_end - corner, n_trace, 'both', 1);
            const l3_ct = _smooth_transition(this.y_trace_end - corner, this.y_trace_end, ncorner, 'end');
            const l_sm = _smooth_transition(this.y_trace_end, this.y_top_start, n_sm_val, 'both');
            const l4 = _smooth_transition(this.y_top_start, this.y_top_end, n_sub_top, 'start');
            
            grid_parts = [l1];
            if (l2.length > 1) grid_parts.push(l2.subarray(1));
            if (l3_cb.length > 1) grid_parts.push(l3_cb.subarray(1));
            if (l3.length > 1) grid_parts.push(l3.subarray(1));
            if (l3_ct.length > 1) grid_parts.push(l3_ct.subarray(1));
            if (l_sm.length > 1) grid_parts.push(l_sm.subarray(1));
            if (l4.length > 1) grid_parts.push(l4.subarray(1));
        }

        if (this.has_top_gnd) {
            const l5 = _smooth_transition(this.y_gnd_top_start, this.y_gnd_top_end, n_gnd, 'start');
            if (l5.length > 1) grid_parts.push(l5.subarray(1));
        }

        // Filter out empty arrays before concatenating
        const filtered_parts = grid_parts.filter(arr => arr.length > 0);
        let y = _concat_arrays(filtered_parts);
        
        y = _enforce_interfaces(y, [
            this.y_gnd_bot_start, this.y_gnd_bot_end,
            this.y_sub_start, this.y_sub_end,
            this.y_trace_start, this.y_trace_end,
            this.y_top_start, this.y_top_end
        ]);
        if (this.has_top_gnd) {
            y = _enforce_interfaces(y, [this.y_gnd_top_start, this.y_gnd_top_end]);
        }

        return y;
    }
    init_matrices() {
        const nx = this.x.length;
        const ny = this.y.length;
        const tol = 1e-11; // Tolerance for floating point comparisons

        this.epsilon_r = Array(ny).fill().map(() => new Float64Array(nx).fill(1.0));
        this.conductor_mask = Array(ny).fill().map(() => new Uint8Array(nx).fill(0));
        this.signal_mask = Array(ny).fill().map(() => new Uint8Array(nx).fill(0));
        this.ground_mask = Array(ny).fill().map(() => new Uint8Array(nx).fill(0));
        this.V = Array(ny).fill().map(() => new Float64Array(nx).fill(0.0));

        const cx = this.domain_width / 2;
        const xl_trace = this.x_tr_l;
        const xr_trace = this.x_tr_r;
        const xl_gap = this.x_gap_l;
        const xr_gap = this.x_gap_r;
        
        for (let i = 0; i < ny; i++) {
            const yc = this.y[i];
            for (let j = 0; j < nx; j++) {
                const xc = this.x[j];

                // --- 1. Permittivity (Dielectrics) ---
                if (yc <= this.y_sub_end + tol) {
                    this.epsilon_r[i][j] = this.er;
                } else {
                    this.epsilon_r[i][j] = this.er_top;
                }

                // Apply Solder Mask if enabled
                if (this.use_sm) {
                    // --- Zone A: Over substrate in gaps (between grounds and signal) ---
                    const is_in_gap_x = (xc > xl_gap - tol && xc < xl_trace + tol) ||
                                        (xc > xr_trace - tol && xc < xr_gap + tol);
                    if (is_in_gap_x && (yc > this.y_sub_end - tol && yc <= this.y_sub_end + this.sm_t_sub + tol)) {
                        this.epsilon_r[i][j] = this.sm_er;
                    }

                    // --- Zone B: On top of metals (signal and top grounds) ---
                    // On top of signal trace
                    if ((xc >= xl_trace - tol && xc <= xr_trace + tol) &&
                        (yc > this.y_trace_end - tol && yc <= this.y_trace_end + this.sm_t_trace + tol)) {
                        this.epsilon_r[i][j] = this.sm_er;
                    }
                    // On top of left ground
                    if ((xc <= xl_gap + tol) &&
                        (yc > this.y_trace_end - tol && yc <= this.y_trace_end + this.sm_t_trace + tol)) {
                        this.epsilon_r[i][j] = this.sm_er;
                    }
                    // On top of right ground
                    if ((xc >= xr_gap - tol) &&
                        (yc > this.y_trace_end - tol && yc <= this.y_trace_end + this.sm_t_trace + tol)) {
                        this.epsilon_r[i][j] = this.sm_er;
                    }

                    // --- Zone C: Conformal coating at metal edges (sidewalls) ---
                    // Left edge of signal trace
                    if ((xc >= xl_trace - this.sm_t_side - tol && xc < xl_trace + tol) &&
                        (yc > this.y_sub_end - tol && yc <= this.y_trace_end + this.sm_t_trace + tol)) {
                        this.epsilon_r[i][j] = this.sm_er;
                    }
                    // Right edge of signal trace
                    if ((xc > xr_trace - tol && xc <= xr_trace + this.sm_t_side + tol) &&
                        (yc > this.y_sub_end - tol && yc <= this.y_trace_end + this.sm_t_trace + tol)) {
                        this.epsilon_r[i][j] = this.sm_er;
                    }
                    // Right edge of left ground (inner edge)
                    if ((xc > xl_gap - tol && xc <= xl_gap + this.sm_t_side + tol) &&
                        (yc > this.y_sub_end - tol && yc <= this.y_trace_end + this.sm_t_trace + tol)) {
                        this.epsilon_r[i][j] = this.sm_er;
                    }
                    // Left edge of right ground (inner edge)
                    if ((xc >= xr_gap - this.sm_t_side - tol && xc < xr_gap + tol) &&
                        (yc > this.y_sub_end - tol && yc <= this.y_trace_end + this.sm_t_trace + tol)) {
                        this.epsilon_r[i][j] = this.sm_er;
                    }
                }

                // --- 2. Conductor Masks ---
                // Bottom Ground Plane
                if (yc <= this.y_gnd_bot_end) {
                    this.ground_mask[i][j] = 1;
                }

                // Top Grounds (Extend from gap edges to domain edges)
                if (yc >= this.y_trace_start && yc <= this.y_trace_end) {
                    if (xc <= xl_gap || xc >= xr_gap) {
                        this.ground_mask[i][j] = 1;
                    }
                }

                // Vias (One per side, from via_gap edge to domain edges)
                if (yc >= this.y_sub_start && yc <= this.y_sub_end) {
                    if (xc <= this.via_x_left_inner || xc >= this.via_x_right_inner) {
                        this.ground_mask[i][j] = 1;
                    }
                }

                // Signal Trace
                if (yc >= this.y_trace_start && yc <= this.y_trace_end) {
                    if (xc >= xl_trace && xc <= xr_trace) {
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
