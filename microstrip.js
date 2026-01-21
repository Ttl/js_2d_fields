import { FieldSolver2D, CONSTANTS, diff } from './field_solver.js';
import { Dielectric, Conductor, Mesher } from './mesher.js';

// ============================================================================
// MICROSTRIP SOLVER V2
// ============================================================================

class MicrostripSolver extends FieldSolver2D {
    constructor(options) {
        super();

        // Store parameters
        this.h = options.substrate_height;
        this.w = options.trace_width;
        this.t = options.trace_thickness;
        this.t_gnd = options.gnd_thickness ?? 35e-6;
        this.er = options.epsilon_r;
        this.er_top = options.epsilon_r_top ?? 1;
        this.tan_delta = options.tan_delta ?? 0.02;
        this.sigma_diel = options.sigma_diel ?? 0.0;
        this.sigma_cond = options.sigma_cond ?? 5.8e7;

        // Differential mode parameters
        this.trace_spacing = options.trace_spacing ?? null;
        this.is_differential = (this.trace_spacing !== null && this.trace_spacing > 0);

        this.gnd_cut_width = options.gnd_cut_width ?? 0.0;
        this.gnd_cut_sub_h = options.gnd_cut_sub_h ?? 0.0;
        this.top_diel_h = options.top_diel_h ?? 0.0;
        this.top_diel_er = options.top_diel_er ?? 1.0;
        this.top_diel_tand = options.top_diel_tand ?? 0.0;

        // Solder Mask Parameters
        this.use_sm = options.use_sm ?? false;
        this.sm_t_sub = options.sm_t_sub ?? 20e-6;
        this.sm_t_trace = options.sm_t_trace ?? 20e-6;
        this.sm_t_side = options.sm_t_side ?? 20e-6;
        this.sm_er = options.sm_er ?? 3.5;
        this.sm_tand = options.sm_tand ?? 0.02;

        this.freq = options.freq ?? 1e9;
        this.omega = 2 * Math.PI * this.freq;
        this.nx = options.nx ?? 300;
        this.ny = options.ny ?? 300;

        // Store air parameters
        const air_side = options.air_side ?? null;
        const air_top = options.air_top ?? null;

        // Domain sizing
        if (this.is_differential) {
            // For differential, span includes both traces and spacing
            const trace_span = 2 * this.w + this.trace_spacing;
            if (air_side === null) {
                this.domain_width = 2 * Math.max(trace_span * 4, this.h * 15);
            } else {
                this.domain_width = trace_span + 2 * air_side;
            }
        } else {
            // Single-ended
            if (air_side === null) {
                this.domain_width = 2 * Math.max(this.w * 8, this.h * 15);
            } else {
                this.domain_width = this.w + 2 * air_side;
            }
        }

        this.boundaries = options.boundaries ?? ["open", "open", "open", "gnd"];

        // Calculate physical coordinates
        this._calculate_coordinates(air_top);

        // Skin depth
        this.delta_s = Math.sqrt(2 / (this.omega * CONSTANTS.MU0 * this.sigma_cond));

        // Build geometry lists
        const [dielectrics, conductors] = this._build_geometry_lists();
        this.dielectrics = dielectrics;
        this.conductors = conductors;

        // Create mesher but don't generate mesh yet
        this.mesher = new Mesher(
            this.domain_width, this.domain_height,
            this.nx, this.ny, this.delta_s,
            this.conductors, this.dielectrics,
            true  // symmetric
        );

        // Mesh will be generated when needed
        this.x = null;
        this.y = null;
        this.dx = null;
        this.dy = null;
        this.mesh_generated = false;
    }

    _calculate_coordinates(air_top) {
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

        // Top dielectric
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

        // Top air/dielectric region
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
    }

    _build_geometry_lists() {
        const dielectrics = [];
        const conductors = [];

        const cx = this.domain_width / 2;
        const xl = cx - this.w / 2;
        const xr = cx + this.w / 2;

        // Substrate (covers both cutout extension and main substrate)
        if (this.gnd_cut_sub_h > 0) {
            dielectrics.push(new Dielectric(
                0, this.y_ext_start,
                this.domain_width, this.y_trace_start - this.y_ext_start,
                this.er, this.tan_delta
            ));
        } else {
            dielectrics.push(new Dielectric(
                0, this.y_sub_start,
                this.domain_width, this.h,
                this.er, this.tan_delta
            ));
        }

        // Top dielectric (if present)
        if (this.top_diel_h > 0) {
            dielectrics.push(new Dielectric(
                0, this.y_top_diel_start,
                this.domain_width, this.top_diel_h,
                this.top_diel_er, this.top_diel_tand
            ));
        }

        // Top air/dielectric region
        dielectrics.push(new Dielectric(
            0, this.y_top_start,
            this.domain_width, this.top_dielectric_h,
            this.er_top, 0.0
        ));

        // Solder mask regions (overwrites previous)
        if (this.use_sm) {
            // Solder mask on substrate (full width)
            dielectrics.push(new Dielectric(
                0, this.y_top_diel_end,
                this.domain_width, this.sm_t_sub,
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
            const xsr = xr + this.sm_t_side;
            if (xsr <= this.domain_width) {
                dielectrics.push(new Dielectric(
                    xr, this.y_trace_start,
                    this.sm_t_side, this.t + this.sm_t_trace,
                    this.sm_er, this.sm_tand
                ));
            }

            // Solder mask on top of trace
            dielectrics.push(new Dielectric(
                xl, this.y_trace_end,
                this.w, this.sm_t_trace,
                this.sm_er, this.sm_tand
            ));
        }

        // --- CONDUCTORS ---

        // Bottom ground (beneath everything)
        if (this.t_gnd > 0) {
            conductors.push(new Conductor(
                0, 0,
                this.domain_width, this.t_gnd,
                false
            ));
        }

        // Bottom ground plane (above cutout extension)
        if (this.gnd_cut_width === 0) {
            // No cutout - full ground plane
            if (this.y_gnd_bot_end > this.y_gnd_bot_start) {
                conductors.push(new Conductor(
                    0, this.y_gnd_bot_start,
                    this.domain_width, this.t_gnd,
                    false
                ));
            }
        } else {
            // With cutout - ground on sides only
            const cut_l = cx - this.gnd_cut_width / 2;
            const cut_r = cx + this.gnd_cut_width / 2;

            // Left ground
            if (cut_l > 0) {
                conductors.push(new Conductor(
                    0, this.y_gnd_bot_start,
                    cut_l, this.t_gnd,
                    false
                ));
            }

            // Right ground
            if (cut_r < this.domain_width) {
                conductors.push(new Conductor(
                    cut_r, this.y_gnd_bot_start,
                    this.domain_width - cut_r, this.t_gnd,
                    false
                ));
            }
        }

        // Signal trace(s)
        if (this.is_differential) {
            // Left trace (negative in odd mode)
            const xl_left = cx - this.w - this.trace_spacing / 2;
            conductors.push(new Conductor(
                xl_left, this.y_trace_start,
                this.w, this.t,
                true
            ));
            // Right trace (positive in odd mode)
            const xl_right = cx + this.trace_spacing / 2;
            conductors.push(new Conductor(
                xl_right, this.y_trace_start,
                this.w, this.t,
                true
            ));
        } else {
            // Single trace
            conductors.push(new Conductor(
                xl, this.y_trace_start,
                this.w, this.t,
                true
            ));
        }

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

    ensure_mesh() {
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

        // Initialize field arrays
        this.V = Array(ny).fill().map(() => new Float64Array(nx));
        this.epsilon_r = Array(ny).fill().map(() => new Float64Array(nx).fill(1.0));
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
                        }
                    }
                }
            }
        }

        // Apply conductors
        let signal_cond_idx = 0;
        for (const cond of this.conductors) {
            for (let i = 0; i < ny; i++) {
                const yc = this.y[i];
                if (yc >= cond.y_min - tol && yc <= cond.y_max + tol) {
                    for (let j = 0; j < nx; j++) {
                        const xc = this.x[j];
                        if (xc >= cond.x_min - tol && xc <= cond.x_max + tol) {
                            if (cond.is_signal) {
                                this.signal_mask[i][j] = 1;
                                // For differential: first signal is negative, second is positive
                                if (this.is_differential) {
                                    if (signal_cond_idx === 0) {
                                        this.signal_n_mask[i][j] = 1;
                                    } else {
                                        this.signal_p_mask[i][j] = 1;
                                    }
                                }
                            } else {
                                this.ground_mask[i][j] = 1;
                            }
                            this.V[i][j] = cond.voltage;
                        }
                    }
                }
            }
            if (cond.is_signal) {
                signal_cond_idx++;
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
    }

}

export { MicrostripSolver };
