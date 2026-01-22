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

        // Coplanar ground options (from GCPW)
        this.use_coplanar_gnd = options.use_coplanar_gnd ?? false;
        this.gap = options.gap ?? 0;                    // Gap from signal to top ground
        this.top_gnd_width = options.top_gnd_width ?? 0; // Width of top ground planes
        this.via_gap = options.via_gap ?? 0;            // Gap from ground edge to via
        this.use_vias = options.use_vias ?? false;      // Enable via generation

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
            if (this.use_coplanar_gnd) {
                // Coplanar: active width includes gaps, top grounds, vias
                const active_width = trace_span + 2 * (this.gap + Math.max(this.top_gnd_width, this.via_gap));
                if (air_side === null) {
                    this.domain_width = Math.max(active_width * 1.5, this.h * 10);
                } else {
                    this.domain_width = active_width + 2 * air_side;
                }
            } else {
                if (air_side === null) {
                    this.domain_width = 2 * Math.max(trace_span * 4, this.h * 15);
                } else {
                    this.domain_width = trace_span + 2 * air_side;
                }
            }
        } else {
            // Single-ended
            if (this.use_coplanar_gnd) {
                // Coplanar: active width includes gaps, top grounds, vias
                const active_width = this.w + 2 * (this.gap + Math.max(this.top_gnd_width, this.via_gap));
                if (air_side === null) {
                    this.domain_width = Math.max(active_width * 1.5, this.h * 10);
                } else {
                    this.domain_width = active_width + 2 * air_side;
                }
            } else {
                if (air_side === null) {
                    this.domain_width = 2 * Math.max(this.w * 8, this.h * 15);
                } else {
                    this.domain_width = this.w + 2 * air_side;
                }
            }
        }

        this.boundaries = options.boundaries ?? ["open", "open", "open", "gnd"];

        // Calculate physical coordinates
        this._calculate_coordinates(air_top);

        // Calculate coplanar geometry if enabled
        if (this.use_coplanar_gnd) {
            this._calculate_coplanar_geometry_x();
        }

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

    _calculate_coplanar_geometry_x() {
        // Calculate x-coordinates for gaps, top grounds, vias
        // Handle both single-ended and differential layouts
        const cx = this.domain_width / 2;

        if (this.is_differential) {
            // Differential: two traces with spacing between
            // Layout: Via|TopGnd|Gap|Sig(-)|Space|Sig(+)|Gap|TopGnd|Via
            const half_spacing = this.trace_spacing / 2;

            // Left trace (negative polarity)
            this.x_tr_left_l = cx - this.w - half_spacing;
            this.x_tr_left_r = cx - half_spacing;

            // Right trace (positive polarity)
            this.x_tr_right_l = cx + half_spacing;
            this.x_tr_right_r = cx + this.w + half_spacing;

            // Outer gaps (from outer edges of traces)
            this.x_gap_outer_l = this.x_tr_left_l - this.gap;
            this.x_gap_outer_r = this.x_tr_right_r + this.gap;

            // Via positions (via_gap is distance from ground edge to via edge)
            this.via_x_left_inner = this.x_gap_outer_l - this.via_gap;
            this.via_x_right_inner = this.x_gap_outer_r + this.via_gap;
        } else {
            // Single-ended: one trace centered
            this.x_tr_l = cx - this.w / 2;
            this.x_tr_r = cx + this.w / 2;

            // Gaps from signal to top ground
            this.x_gap_l = this.x_tr_l - this.gap;
            this.x_gap_r = this.x_tr_r + this.gap;

            // Via positions (via_gap is distance from ground edge to via edge)
            this.via_x_left_inner = this.x_gap_l - this.via_gap;
            this.via_x_right_inner = this.x_gap_r + this.via_gap;
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
            if (this.use_coplanar_gnd) {
                // Coplanar solder mask: in gaps between signals and grounds
                this._add_coplanar_solder_mask(dielectrics);
            } else {
                // Standard microstrip solder mask
                // Solder mask on substrate (full width)
                dielectrics.push(new Dielectric(
                    0, this.y_sub_end,
                    this.domain_width, this.sm_t_sub,
                    this.sm_er, this.sm_tand
                ));

                if (this.is_differential) {
                    // Differential: solder mask for both traces
                    const half_spacing = this.trace_spacing / 2;
                    const xl_left = cx - this.w - half_spacing;
                    const xr_left = cx - half_spacing;
                    const xl_right = cx + half_spacing;
                    const xr_right = cx + this.w + half_spacing;

                    // Left trace solder mask
                    dielectrics.push(new Dielectric(
                        xl_left - this.sm_t_side, this.y_trace_start,
                        this.sm_t_side, this.t + this.sm_t_trace,
                        this.sm_er, this.sm_tand
                    ));
                    dielectrics.push(new Dielectric(
                        xr_left, this.y_trace_start,
                        this.sm_t_side, this.t + this.sm_t_trace,
                        this.sm_er, this.sm_tand
                    ));
                    dielectrics.push(new Dielectric(
                        xl_left, this.y_trace_end,
                        this.w, this.sm_t_trace,
                        this.sm_er, this.sm_tand
                    ));

                    // Right trace solder mask
                    dielectrics.push(new Dielectric(
                        xl_right - this.sm_t_side, this.y_trace_start,
                        this.sm_t_side, this.t + this.sm_t_trace,
                        this.sm_er, this.sm_tand
                    ));
                    dielectrics.push(new Dielectric(
                        xr_right, this.y_trace_start,
                        this.sm_t_side, this.t + this.sm_t_trace,
                        this.sm_er, this.sm_tand
                    ));
                    dielectrics.push(new Dielectric(
                        xl_right, this.y_trace_end,
                        this.w, this.sm_t_trace,
                        this.sm_er, this.sm_tand
                    ));
                } else {
                    // Single-ended: solder mask for one trace
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
            }
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

        // Coplanar vias through substrate (if enabled)
        if (this.use_coplanar_gnd && this.use_vias) {
            // Vias should extend from bottom ground to top coplanar grounds
            // Start from top of bottom ground (y_ext_start = t_gnd)
            // End at top of coplanar grounds (y_trace_end)
            const via_y_start = this.y_ext_start;
            const via_height = this.y_trace_end - this.y_ext_start;

            // Left via (from inner edge to left boundary)
            if (this.via_x_left_inner > 0) {
                conductors.push(new Conductor(
                    0, via_y_start,
                    this.via_x_left_inner, via_height,
                    false
                ));
            }

            // Right via (from inner edge to right boundary)
            if (this.via_x_right_inner < this.domain_width) {
                conductors.push(new Conductor(
                    this.via_x_right_inner, via_y_start,
                    this.domain_width - this.via_x_right_inner, via_height,
                    false
                ));
            }
        }

        // Signal trace(s)
        if (this.is_differential) {
            // Left trace (negative in odd mode, polarity = -1)
            const xl_left = cx - this.w - this.trace_spacing / 2;
            conductors.push(new Conductor(
                xl_left, this.y_trace_start,
                this.w, this.t,
                true, -1
            ));
            // Right trace (positive in odd mode, polarity = +1)
            const xl_right = cx + this.trace_spacing / 2;
            conductors.push(new Conductor(
                xl_right, this.y_trace_start,
                this.w, this.t,
                true, 1
            ));
        } else {
            // Single trace (polarity = +1)
            conductors.push(new Conductor(
                xl, this.y_trace_start,
                this.w, this.t,
                true, 1
            ));
        }

        // Coplanar top grounds (on same layer as signal traces)
        if (this.use_coplanar_gnd) {
            if (this.is_differential) {
                // Differential: grounds on outer edges only
                // Left top ground (from left edge to outer gap edge)
                conductors.push(new Conductor(
                    0, this.y_trace_start,
                    this.x_gap_outer_l, this.t,
                    false
                ));

                // Right top ground (from outer gap edge to right edge)
                conductors.push(new Conductor(
                    this.x_gap_outer_r, this.y_trace_start,
                    this.domain_width - this.x_gap_outer_r, this.t,
                    false
                ));
            } else {
                // Single-ended: grounds on both sides of the trace
                // Left top ground (from left edge to gap edge)
                conductors.push(new Conductor(
                    0, this.y_trace_start,
                    this.x_gap_l, this.t,
                    false
                ));

                // Right top ground (from gap edge to right edge)
                conductors.push(new Conductor(
                    this.x_gap_r, this.y_trace_start,
                    this.domain_width - this.x_gap_r, this.t,
                    false
                ));
            }
        }

        // Top ground plane (if present - for stripline)
        if (this.has_top_gnd) {
            conductors.push(new Conductor(
                0, this.y_gnd_top_start,
                this.domain_width, this.t_gnd,
                false
            ));
        }

        return [dielectrics, conductors];
    }

    _add_coplanar_solder_mask(dielectrics) {
        // Coplanar solder mask: covers gaps and tops of conductors
        // Adapted from gcpw.js lines 194-283

        if (this.is_differential) {
            // Differential coplanar solder mask
            const xl = this.x_tr_left_l;
            const xr_left = this.x_tr_left_r;
            const xl_right = this.x_tr_right_l;
            const xr = this.x_tr_right_r;
            const xl_gap = this.x_gap_outer_l;
            const xr_gap = this.x_gap_outer_r;

            // Solder mask on substrate in outer gaps
            const xl_sub_start = xl_gap + this.sm_t_side;
            const xl_sub_end = xl - this.sm_t_side;
            if (xl_sub_end > xl_sub_start) {
                dielectrics.push(new Dielectric(
                    xl_sub_start, this.y_sub_end,
                    xl_sub_end - xl_sub_start, this.sm_t_sub,
                    this.sm_er, this.sm_tand
                ));
            }

            const xr_sub_start = xr + this.sm_t_side;
            const xr_sub_end = xr_gap - this.sm_t_side;
            if (xr_sub_end > xr_sub_start) {
                dielectrics.push(new Dielectric(
                    xr_sub_start, this.y_sub_end,
                    xr_sub_end - xr_sub_start, this.sm_t_sub,
                    this.sm_er, this.sm_tand
                ));
            }

            // Solder mask in center gap (between traces)
            const center_start = xr_left + this.sm_t_side;
            const center_end = xl_right - this.sm_t_side;
            if (center_end > center_start) {
                dielectrics.push(new Dielectric(
                    center_start, this.y_sub_end,
                    center_end - center_start, this.sm_t_sub,
                    this.sm_er, this.sm_tand
                ));
            }

            // Solder mask on sides of left trace
            dielectrics.push(new Dielectric(
                xl - this.sm_t_side, this.y_trace_start,
                this.sm_t_side, this.t + this.sm_t_trace,
                this.sm_er, this.sm_tand
            ));
            dielectrics.push(new Dielectric(
                xr_left, this.y_trace_start,
                this.sm_t_side, this.t + this.sm_t_trace,
                this.sm_er, this.sm_tand
            ));

            // Solder mask on sides of right trace
            dielectrics.push(new Dielectric(
                xl_right - this.sm_t_side, this.y_trace_start,
                this.sm_t_side, this.t + this.sm_t_trace,
                this.sm_er, this.sm_tand
            ));
            dielectrics.push(new Dielectric(
                xr, this.y_trace_start,
                this.sm_t_side, this.t + this.sm_t_trace,
                this.sm_er, this.sm_tand
            ));

            // Solder mask on outer gap sides (ground side)
            dielectrics.push(new Dielectric(
                xl_gap, this.y_trace_start,
                this.sm_t_side, this.t + this.sm_t_trace,
                this.sm_er, this.sm_tand
            ));
            dielectrics.push(new Dielectric(
                xr_gap - this.sm_t_side, this.y_trace_start,
                this.sm_t_side, this.t + this.sm_t_trace,
                this.sm_er, this.sm_tand
            ));

            // Solder mask on top of traces
            dielectrics.push(new Dielectric(
                xl, this.y_trace_end,
                this.w, this.sm_t_trace,
                this.sm_er, this.sm_tand
            ));
            dielectrics.push(new Dielectric(
                xl_right, this.y_trace_end,
                this.w, this.sm_t_trace,
                this.sm_er, this.sm_tand
            ));

            // Solder mask on top of grounds
            dielectrics.push(new Dielectric(
                0, this.y_trace_end,
                xl_gap, this.sm_t_trace,
                this.sm_er, this.sm_tand
            ));
            dielectrics.push(new Dielectric(
                xr_gap, this.y_trace_end,
                this.domain_width - xr_gap, this.sm_t_trace,
                this.sm_er, this.sm_tand
            ));
        } else {
            // Single-ended coplanar solder mask (from original gcpw.js)
            const xl = this.x_tr_l;
            const xr = this.x_tr_r;
            const xl_gap = this.x_gap_l;
            const xr_gap = this.x_gap_r;

            // Trace side positions
            const xsl = xl - this.sm_t_side;
            const xsr = xr + this.sm_t_side;

            // Ground side mask positions
            const xl_gnd_side_end = Math.min(xl_gap + this.sm_t_side, xl);
            const xr_gnd_side_start = Math.max(xr_gap - this.sm_t_side, xr);

            // Solder mask on substrate in gaps (between grounds and signal)
            // Left gap
            const xl_sub_start = xl_gnd_side_end;
            const xl_sub_end = Math.min(xl, Math.max(xl_sub_start, xsl));
            if (xl_sub_end > xl_sub_start) {
                dielectrics.push(new Dielectric(
                    xl_sub_start, this.y_sub_end,
                    xl_sub_end - xl_sub_start, this.sm_t_sub,
                    this.sm_er, this.sm_tand
                ));
            }

            // Right gap
            const xr_sub_end = xr_gnd_side_start;
            const xr_sub_start_calc = Math.max(xr, Math.min(xr_sub_end, xsr));
            if (xr_sub_end > xr_sub_start_calc) {
                dielectrics.push(new Dielectric(
                    xr_sub_start_calc, this.y_sub_end,
                    xr_sub_end - xr_sub_start_calc, this.sm_t_sub,
                    this.sm_er, this.sm_tand
                ));
            }

            // Solder mask on left side of trace
            if (xsl >= 0) {
                dielectrics.push(new Dielectric(
                    xsl, this.y_trace_start,
                    this.sm_t_side, this.t + this.sm_t_trace,
                    this.sm_er, this.sm_tand
                ));
            }

            // Solder mask on right side of trace
            if (xsr <= this.domain_width) {
                dielectrics.push(new Dielectric(
                    xr, this.y_trace_start,
                    this.sm_t_side, this.t + this.sm_t_trace,
                    this.sm_er, this.sm_tand
                ));
            }

            // Solder mask on ground side of left gap
            if (xl_gnd_side_end > xl_gap) {
                dielectrics.push(new Dielectric(
                    xl_gap, this.y_trace_start,
                    xl_gnd_side_end - xl_gap, this.t + this.sm_t_trace,
                    this.sm_er, this.sm_tand
                ));
            }

            // Solder mask on ground side of right gap
            if (xr_gap > xr_gnd_side_start) {
                dielectrics.push(new Dielectric(
                    xr_gnd_side_start, this.y_trace_start,
                    xr_gap - xr_gnd_side_start, this.t + this.sm_t_trace,
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

        // Initialize mask and material arrays (V is created by solver based on mode)
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

        // Apply conductors - use polarity to determine signal_p vs signal_n
        for (const cond of this.conductors) {
            for (let i = 0; i < ny; i++) {
                const yc = this.y[i];
                if (yc >= cond.y_min - tol && yc <= cond.y_max + tol) {
                    for (let j = 0; j < nx; j++) {
                        const xc = this.x[j];
                        if (xc >= cond.x_min - tol && xc <= cond.x_max + tol) {
                            if (cond.is_signal) {
                                this.signal_mask[i][j] = 1;
                                // Use polarity to determine positive/negative trace
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
    }

}

export { MicrostripSolver };
