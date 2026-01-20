import { FieldSolver2D, CONSTANTS, diff } from './field_solver.js';
import { Dielectric, Conductor, Mesher } from './mesher.js';

export class GroundedCPWSolver2D extends FieldSolver2D {
    constructor(options = {}) {
        super();

        // Extract parameters with defaults
        const {
            substrate_height,
            trace_width,
            trace_thickness,
            gap,
            top_gnd_width,
            via_gap,
            gnd_thickness = 35e-6,
            epsilon_r = 4.5,
            tan_delta = 0.02,
            sigma_diel = 0.0,
            sigma_cond = 5.8e7,
            epsilon_r_top = 1,
            air_top = null,
            air_side = null,
            freq = 1e9,
            nx = 300,
            ny = 300,
            boundaries = null,
            use_sm = false,
            sm_t_sub = 20e-6,
            sm_t_trace = 20e-6,
            sm_t_side = 20e-6,
            sm_er = 3.5,
            sm_tand = 0.02
        } = options;

        // Store parameters
        this.h = substrate_height;
        this.w = trace_width;
        this.t = trace_thickness;
        this.gap = gap;
        this.top_gnd_w = top_gnd_width;
        this.via_gap = via_gap;
        this.t_gnd = gnd_thickness;
        this.er = epsilon_r;
        this.er_top = epsilon_r_top;
        this.tan_delta = tan_delta;
        this.sigma_diel = sigma_diel;
        this.sigma_cond = sigma_cond;

        // Solder Mask Parameters
        this.use_sm = use_sm;
        this.sm_t_sub = sm_t_sub;
        this.sm_t_trace = sm_t_trace;
        this.sm_t_side = sm_t_side;
        this.sm_er = sm_er;
        this.sm_tand = sm_tand;

        this.freq = freq;
        this.omega = 2 * CONSTANTS.PI * freq;
        this.nx = nx;
        this.ny = ny;

        // Store air parameters
        this.air_top = air_top;
        this.air_side = air_side;

        this.active_width = this.w + 2 * Math.max(this.gap + this.top_gnd_w, this.via_gap);

        // Domain sizing
        if (this.air_side === null) {
            this.domain_width = this.active_width * 1.5;
        } else {
            this.domain_width = this.active_width + 2 * this.air_side;
        }

        this.boundaries = boundaries || ["open", "open", "open", "gnd"];

        // Calculate physical coordinates
        this._calculate_coordinates();

        // Skin depth
        this.delta_s = Math.sqrt(2 / (this.omega * CONSTANTS.MU0 * this.sigma_cond));

        // Calculate geometry X coordinates
        this._calculate_geometry_x();

        // Build geometry lists
        const [dielectrics, conductors] = this._build_geometry_lists();
        this.dielectrics = dielectrics;
        this.conductors = conductors;

        // Create mesher and generate mesh
        const mesher = new Mesher(
            this.domain_width, this.domain_height,
            this.nx, this.ny, this.delta_s,
            this.conductors, this.dielectrics,
            true  // symmetric
        );

        [this.x, this.y] = mesher.generate_mesh();

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
    }

    _calculate_coordinates() {
        // Bottom ground
        this.y_gnd_bot_start = 0.0;
        this.y_gnd_bot_end = this.t_gnd;

        // Substrate
        this.y_sub_start = this.y_gnd_bot_end;
        this.y_sub_end = this.y_sub_start + this.h;

        // Trace
        this.y_trace_start = this.y_sub_end;
        this.y_trace_end = this.y_trace_start + this.t;

        // Solder mask extents
        this.y_sm_sub_end = this.y_sub_end + this.sm_t_sub;
        this.y_sm_trace_end = this.y_trace_end + this.sm_t_trace;

        this.y_top_start = this.y_trace_end;
        if (this.use_sm) {
            this.y_top_start = Math.max(this.y_sm_sub_end, this.y_sm_trace_end);
        }

        // Top air/dielectric region
        if (this.air_top === null) {
            this.top_dielectric_h = this.h * 10;
            this.has_top_gnd = false;
        } else {
            this.top_dielectric_h = this.air_top + this.t;
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
    }

    _calculate_geometry_x() {
        const cx = this.domain_width / 2;

        // Signal trace
        this.x_tr_l = cx - this.w / 2;
        this.x_tr_r = cx + this.w / 2;

        // Top ground gaps
        this.x_gap_l = this.x_tr_l - this.gap;
        this.x_gap_r = this.x_tr_r + this.gap;

        // Via positions
        this.via_x_left_inner = this.x_tr_l - this.via_gap;
        this.via_x_right_inner = this.x_tr_r + this.via_gap;
    }

    _build_geometry_lists() {
        const dielectrics = [];
        const conductors = [];

        const cx = this.domain_width / 2;
        const xl = this.x_tr_l;
        const xr = this.x_tr_r;
        const xl_gap = this.x_gap_l;
        const xr_gap = this.x_gap_r;

        // --- DIELECTRICS ---

        // Substrate (full width)
        dielectrics.push(new Dielectric(
            0, this.y_sub_start,
            this.domain_width, this.h,
            this.er, this.tan_delta
        ));

        // Top air/dielectric region
        dielectrics.push(new Dielectric(
            0, this.y_top_start,
            this.domain_width, this.top_dielectric_h,
            this.er_top, 0.0
        ));

        // Solder mask regions (overwrites previous)
        if (this.use_sm) {
            // Solder mask on substrate in gaps (between grounds and signal)
            // Left gap
            dielectrics.push(new Dielectric(
                xl_gap, this.y_sub_end,
                xl - xl_gap, this.sm_t_sub,
                this.sm_er, this.sm_tand
            ));
            // Right gap
            dielectrics.push(new Dielectric(
                xr, this.y_sub_end,
                xr_gap - xr, this.sm_t_sub,
                this.sm_er, this.sm_tand
            ));

            // Solder mask on left side of trace
            const xsl = xl - this.sm_t_side;
            if (xsl >= 0) {
                dielectrics.push(new Dielectric(
                    xsl, this.y_trace_start,
                    this.sm_t_side, this.t + this.sm_t_trace,
                    this.sm_er, this.sm_tand
                ));
            }

            // Solder mask on right side of trace
            const xsr = xr;
            if (xsr + this.sm_t_side <= this.domain_width) {
                dielectrics.push(new Dielectric(
                    xsr, this.y_trace_start,
                    this.sm_t_side, this.t + this.sm_t_trace,
                    this.sm_er, this.sm_tand
                ));
            }

            // Solder mask on top of signal trace
            dielectrics.push(new Dielectric(
                xl, this.y_trace_end,
                this.w, this.sm_t_trace,
                this.sm_er, this.sm_tand
            ));

            // Solder mask on top of left ground
            dielectrics.push(new Dielectric(
                0, this.y_trace_end,
                xl_gap, this.sm_t_trace,
                this.sm_er, this.sm_tand
            ));

            // Solder mask on top of right ground
            dielectrics.push(new Dielectric(
                xr_gap, this.y_trace_end,
                this.domain_width - xr_gap, this.sm_t_trace,
                this.sm_er, this.sm_tand
            ));
        }

        // --- CONDUCTORS ---

        // Bottom ground plane (full width)
        if (this.t_gnd > 0) {
            conductors.push(new Conductor(
                0, this.y_gnd_bot_start,
                this.domain_width, this.t_gnd,
                false
            ));
        }

        // Left via (from inner edge to left boundary)
        if (this.via_x_left_inner > 0) {
            conductors.push(new Conductor(
                0, this.y_sub_start,
                this.via_x_left_inner, this.h,
                false
            ));
        }

        // Right via (from inner edge to right boundary)
        if (this.via_x_right_inner < this.domain_width) {
            conductors.push(new Conductor(
                this.via_x_right_inner, this.y_sub_start,
                this.domain_width - this.via_x_right_inner, this.h,
                false
            ));
        }

        // Left top ground (from left edge to gap edge)
        conductors.push(new Conductor(
            0, this.y_trace_start,
            xl_gap, this.t,
            false
        ));

        // Right top ground (from gap edge to right edge)
        conductors.push(new Conductor(
            xr_gap, this.y_trace_start,
            this.domain_width - xr_gap, this.t,
            false
        ));

        // Signal trace
        conductors.push(new Conductor(
            xl, this.y_trace_start,
            this.w, this.t,
            true
        ));

        // Top ground plane (if present)
        if (this.has_top_gnd) {
            conductors.push(new Conductor(
                0, this.y_gnd_top_start,
                this.domain_width, this.t_gnd,
                false
            ));
        }

        return [dielectrics, conductors];
    }

    _setup_geometry() {
        const tol = 1e-11;
        const nx = this.x.length;
        const ny = this.y.length;

        const dielectrics = this.dielectrics;
        const conductors = this.conductors;

        // Initialize field arrays
        this.V = Array(ny).fill().map(() => new Float64Array(nx).fill(0.0));
        this.epsilon_r = Array(ny).fill().map(() => new Float64Array(nx).fill(1.0));
        this.signal_mask = Array(ny).fill().map(() => new Uint8Array(nx).fill(0));
        this.ground_mask = Array(ny).fill().map(() => new Uint8Array(nx).fill(0));

        // Apply dielectrics (last overwrites)
        for (const diel of dielectrics) {
            for (let i = 0; i < ny; i++) {
                const yc = this.y[i];
                if (yc >= diel.y_min - tol && yc <= diel.y_max + tol) {
                    for (let j = 0; j < nx; j++) {
                        const xc = this.x[j];
                        if (xc >= diel.x_min - tol && xc <= diel.x_max + tol) {
                            this.epsilon_r[i][j] = diel.epsilon_r;
                        }
                    }
                }
            }
        }

        // Apply conductors
        for (const cond of conductors) {
            for (let i = 0; i < ny; i++) {
                const yc = this.y[i];
                if (yc >= cond.y_min - tol && yc <= cond.y_max + tol) {
                    for (let j = 0; j < nx; j++) {
                        const xc = this.x[j];
                        if (xc >= cond.x_min - tol && xc <= cond.x_max + tol) {
                            if (cond.is_signal) {
                                this.signal_mask[i][j] = 1;
                            } else {
                                this.ground_mask[i][j] = 1;
                            }
                            this.V[i][j] = cond.voltage;
                        }
                    }
                }
            }
        }

        // Finalize conductor mask
        this.conductor_mask = Array(ny).fill().map(() => new Uint8Array(nx).fill(0));
        for (let i = 0; i < ny; i++) {
            for (let j = 0; j < nx; j++) {
                this.conductor_mask[i][j] = this.signal_mask[i][j] | this.ground_mask[i][j];
            }
        }
    }
}
