import { FieldSolver2D, CONSTANTS, diff } from './field_solver.js';
import { Dielectric, Conductor, Mesher } from './mesher.js';

class MicrostripSolver extends FieldSolver2D {
    constructor(options) {
        super();

        // Validate parameters before proceeding
        this._validate_parameters(options);

        // Store parameters
        this.h = options.substrate_height;
        this.w = options.trace_width;
        this.t = options.trace_thickness;
        this.t_gnd = options.gnd_thickness ?? 35e-6;
        this.er = options.epsilon_r;
        this.er_top = options.epsilon_r_top ?? 1;
        this.tand_top = options.tan_delta_top ?? 0;
        this.tan_delta = options.tan_delta ?? 0.02;
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
        this.gap = options.gap ?? 0; // Gap from signal to top ground
        this.top_gnd_width = options.top_gnd_width ?? 0; // Width of top ground planes
        this.via_gap = options.via_gap ?? 0; // Gap from ground edge to via
        this.use_vias = options.use_vias ?? false; // Enable via generation

        // Enclosure options
        this.enclosure_width = options.enclosure_width ?? null;
        this.enclosure_width = this.enclosure_width === "auto" ? null : this.enclosure_width;
        this.enclosure_height = options.enclosure_height ?? null;
        this.enclosure_height = this.enclosure_height === "auto" ? null : this.enclosure_height;

        // Solder Mask Parameters
        this.use_sm = options.use_sm ?? false;
        this.sm_t_sub = options.sm_t_sub ?? 20e-6;
        this.sm_t_trace = options.sm_t_trace ?? 20e-6;
        this.sm_t_side = options.sm_t_side ?? 20e-6;
        this.sm_er = options.sm_er ?? 3.5;
        this.sm_tand = options.sm_tand ?? 0.02;

        this.freq = options.freq ?? 1e9;
        this.nx = options.nx ?? 300;
        this.ny = options.ny ?? 300;

        // Surface roughness parameter (RMS roughness in meters)
        this.rq = options.rq ?? 0;

        // Domain sizing
        if (this.enclosure_width !== null) {
            // Use explicit enclosure width
            // If side ground is enabled, this is the distance between inner walls
            // Otherwise, it's the total domain width
            if (this.has_side_gnd) {
                // enclosure_width is the air box size (inner wall to inner wall)
                // Add ground thickness on both sides
                this.domain_width = this.enclosure_width + 2 * this.t_gnd;
            } else {
                this.domain_width = this.enclosure_width;
            }
        } else if (this.is_differential) {
            // For differential, span includes both traces and spacing
            const trace_span = 2 * this.w + this.trace_spacing;
            if (this.use_coplanar_gnd) {
                // Coplanar: active width includes gaps, top grounds, vias
                const active_width = trace_span + 2 * (this.gap + Math.max(this.top_gnd_width, this.via_gap));
                this.domain_width = Math.max(active_width * 1.5, this.h * 10);
            } else {
                this.domain_width = 2 * Math.max(trace_span * 4, this.h * 15);
            }
        } else {
            // Single-ended
            if (this.use_coplanar_gnd) {
                // Coplanar: active width includes gaps, top grounds, vias
                const active_width = this.w + 2 * (this.gap + Math.max(this.top_gnd_width, this.via_gap));
                this.domain_width = Math.max(active_width * 1.5, this.h * 10);
            } else {
                this.domain_width = 2 * Math.max(this.w * 8, this.h * 15);
            }
        }

        this.boundaries = options.boundaries ?? ["open", "open", "open", "gnd"];

        // Determine if side grounds are present based on boundaries
        this.has_side_gnd = (this.boundaries[0] === "gnd" || this.boundaries[1] === "gnd");

        // Calculate physical coordinates
        this._calculate_coordinates();

        // Calculate coplanar geometry if enabled
        if (this.use_coplanar_gnd) {
            this._calculate_coplanar_geometry_x();
        }

        // Skin depth
        this.delta_s = Math.sqrt(2 / (2 * Math.PI * Math.max(1e6, this.freq) * CONSTANTS.MU0 * this.sigma_cond));

        // Build geometry lists
        const [dielectrics, conductors] = this._build_geometry_lists();
        this.dielectrics = dielectrics;
        this.conductors = conductors;

        // Create mesher but don't generate mesh yet
        // Geometry is centered at x=0, so domain spans from -domain_width/2 to +domain_width/2
        // Y-coordinate system: bottom ground extends from -t_gnd to 0 (top surface at y=0)
        this.mesher = new Mesher(
            this.domain_width, this.domain_height,
            this.nx, this.ny, this.delta_s,
            this.conductors, this.dielectrics,
            true,  // symmetric
            -this.domain_width / 2,  // x_min
            this.domain_width / 2,   // x_max
            -this.t_gnd,             // y_min (bottom of bottom ground)
            this.domain_height       // y_max (top of domain)
        );

        // Mesh will be generated when needed
        this.x = null;
        this.y = null;
        this.dx = null;
        this.dy = null;
        this.mesh_generated = false;
    }

    _validate_parameters(options) {
        const errors = [];

        const isValidNumber = (val) => typeof val === 'number' && !isNaN(val) && isFinite(val);

        const checkPositive = (val, name) => {
            if (val === undefined || val === null) return;
            if (!isValidNumber(val)) {
                errors.push(`${name} must be a valid number (got ${val})`);
            } else if (val <= 0) {
                errors.push(`${name} must be positive (> 0), got ${val}`);
            }
        };

        const checkNonNegative = (val, name) => {
            if (val === undefined || val === null) return;
            if (!isValidNumber(val)) {
                errors.push(`${name} must be a valid number (got ${val})`);
            } else if (val < 0) {
                errors.push(`${name} must be non-negative (>= 0), got ${val}`);
            }
        };

        // Required positive parameters
        checkPositive(options.substrate_height, 'substrate_height');
        checkPositive(options.trace_width, 'trace_width');
        checkPositive(options.epsilon_r, 'epsilon_r');

        // Allow negative thickness for conductors inside the substrate
        isValidNumber(options.trace_thickness, 'trace_thickness');

        checkNonNegative(options.freq, 'frequency');

        // Optional positive parameters
        if (options.gnd_thickness !== undefined) checkPositive(options.gnd_thickness, 'gnd_thickness');
        if (options.epsilon_r_top !== undefined) checkPositive(options.epsilon_r_top, 'epsilon_r_top');
        if (options.sigma_cond !== undefined) checkPositive(options.sigma_cond, 'sigma_cond');
        if (options.trace_spacing !== undefined && options.trace_spacing !== null && options.trace_spacing > 0) {
            checkPositive(options.trace_spacing, 'trace_spacing');
        }

        // Non-negative parameters
        checkNonNegative(options.tan_delta, 'tan_delta');
        checkNonNegative(options.gnd_cut_width, 'gnd_cut_width');
        checkNonNegative(options.gnd_cut_sub_h, 'gnd_cut_sub_h');
        checkNonNegative(options.top_diel_h, 'top_diel_h');
        checkNonNegative(options.top_diel_tand, 'top_diel_tand');
        checkNonNegative(options.gap, 'gap');
        checkNonNegative(options.top_gnd_width, 'top_gnd_width');
        checkNonNegative(options.via_gap, 'via_gap');

        if (options.trace_thickness <= -options.substrate_height) {
            errors.push("trace_thickness must be > -substrate_height");
        }

        if (options.enclosure_height != null && (options.trace_thickness > options.enclosure_height)) {
            errors.push("trace_thickness > enclosure_height");
        }

        // Top dielectric epsilon_r should be positive if top dielectric is used
        if (options.top_diel_h > 0 && options.top_diel_er !== undefined) {
            checkPositive(options.top_diel_er, 'top_diel_er');
        }

        // Solder mask parameters - sm_t_sub and sm_t_trace can be negative
        // But sm_er must be positive if solder mask is used
        if (options.use_sm) {
            if (options.sm_t_sub !== undefined && !isValidNumber(options.sm_t_sub)) {
                errors.push(`sm_t_sub must be a valid number (got ${options.sm_t_sub})`);
            }
            if (options.sm_t_trace !== undefined && !isValidNumber(options.sm_t_trace)) {
                errors.push(`sm_t_trace must be a valid number (got ${options.sm_t_trace})`);
            }
            checkNonNegative(options.sm_t_side, 'sm_t_side');
            checkPositive(options.sm_er, 'sm_er');
            checkNonNegative(options.sm_tand, 'sm_tand');
        }

        // Enclosure parameters
        if (options.enclosure_width !== undefined && options.enclosure_width !== "auto") {
            checkPositive(options.enclosure_width, 'enclosure_width');
        }
        if (options.enclosure_height !== undefined && options.enclosure_height !== "auto") {
            checkPositive(options.enclosure_height, 'enclosure_height');
        }

        // Grid size parameters
        if (options.nx !== undefined) {
            if (!Number.isInteger(options.nx) || options.nx <= 0) {
                errors.push(`nx must be a positive integer (got ${options.nx})`);
            }
        }
        if (options.ny !== undefined) {
            if (!Number.isInteger(options.ny) || options.ny <= 0) {
                errors.push(`ny must be a positive integer (got ${options.ny})`);
            }
        }

        // Check that active area fits in enclosure if enclosure is enabled
        const has_side_gnd_temp = (options.boundaries && (options.boundaries[0] === "gnd" || options.boundaries[1] === "gnd"));
        if (has_side_gnd_temp && options.enclosure_width !== undefined && options.enclosure_width !== "auto") {
            // Calculate active area width based on configuration
            let active_width = 0;
            const w = options.trace_width || 0;
            const trace_spacing = options.trace_spacing || 0;
            const is_differential = (trace_spacing !== null && trace_spacing > 0);

            if (is_differential) {
                const trace_span = 2 * w + trace_spacing;
                if (options.use_coplanar_gnd) {
                    const gap = options.gap || 0;
                    const top_gnd_width = options.top_gnd_width || 0;
                    const via_gap = options.via_gap || 0;
                    active_width = trace_span + 2 * (gap + Math.max(top_gnd_width, via_gap));
                } else {
                    active_width = trace_span;
                }
            } else {
                if (options.use_coplanar_gnd) {
                    const gap = options.gap || 0;
                    const top_gnd_width = options.top_gnd_width || 0;
                    const via_gap = options.via_gap || 0;
                    active_width = w + 2 * (gap + Math.max(top_gnd_width, via_gap));
                } else {
                    active_width = w;
                }
            }

            // enclosure_width is the inner width (between ground walls), so just check active area
            if (active_width > options.enclosure_width) {
                errors.push(`Active area width (${(active_width * 1000).toFixed(3)} mm) exceeds enclosure inner width (${(options.enclosure_width * 1000).toFixed(3)} mm)`);
            }
        }

        // If there are any errors, throw them
        if (errors.length > 0) {
            throw new Error('Parameter validation failed:\n' + errors.map(e => '  - ' + e).join('\n'));
        }
    }

    _calculate_coordinates() {
        // Bottom extension for cut ground
        this.y_ext_start = 0;
        this.y_ext_end = this.gnd_cut_sub_h;

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

        // Trace position - starts at top of top dielectric (or substrate if no top diel)
        // Negative thickness means trace extends down into substrate
        this.y_trace_start = this.y_top_diel_start;
        this.y_trace_end = this.y_trace_start + this.t;

        // Solder mask extents
        this.y_sm_sub_end = this.y_top_diel_end + this.sm_t_sub;
        // For negative trace thickness, solder mask top should be at substrate surface
        const y_trace_top = Math.max(this.y_trace_start, this.y_trace_end);
        this.y_sm_trace_end = y_trace_top + this.sm_t_trace;

        this.y_top_start = this.y_top_diel_end;

        // Top air/dielectric region
        if (this.enclosure_height !== null) {
            // Use enclosure height (distance from highest dielectric to top of domain)
            // Highest dielectric is y_top_start
            this.top_dielectric_h = this.enclosure_height;
            this.has_top_gnd = (this.boundaries[2] === "gnd");
        } else {
            this.top_dielectric_h = (this.h + this.t) * 15;
            this.has_top_gnd = false;
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
        // Geometry is centered at x=0

        if (this.is_differential) {
            // Differential: two traces with spacing between
            // Layout: Via|TopGnd|Gap|Sig(-)|Space|Sig(+)|Gap|TopGnd|Via
            const half_spacing = this.trace_spacing / 2;

            // Left trace (negative polarity)
            this.x_tr_left_l = -this.w - half_spacing;
            this.x_tr_left_r = -half_spacing;

            // Right trace (positive polarity)
            this.x_tr_right_l = half_spacing;
            this.x_tr_right_r = this.w + half_spacing;

            // Outer gaps (from outer edges of traces)
            this.x_gap_outer_l = this.x_tr_left_l - this.gap;
            this.x_gap_outer_r = this.x_tr_right_r + this.gap;

            // Via positions (via_gap is distance from ground edge to via edge)
            this.via_x_left_inner = this.x_gap_outer_l - this.via_gap;
            this.via_x_right_inner = this.x_gap_outer_r + this.via_gap;
        } else {
            // Single-ended: one trace centered at x=0
            this.x_tr_l = -this.w / 2;
            this.x_tr_r = this.w / 2;

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

        // Geometry is centered at x=0
        const xl = -this.w / 2;
        const xr = this.w / 2;
        const x_min = -this.domain_width / 2;
        const x_max = this.domain_width / 2;

        // Substrate (covers both cutout extension and main substrate)
        if (this.gnd_cut_sub_h > 0) {
            dielectrics.push(new Dielectric(
                x_min, this.y_ext_start,
                this.domain_width, this.y_trace_start - this.y_ext_start,
                this.er, this.tan_delta
            ));
        } else {
            dielectrics.push(new Dielectric(
                x_min, this.y_sub_start,
                this.domain_width, this.h,
                this.er, this.tan_delta
            ));
        }

        // Top dielectric (if present)
        if (this.top_diel_h > 0) {
            dielectrics.push(new Dielectric(
                x_min, this.y_top_diel_start,
                this.domain_width, this.top_diel_h,
                this.top_diel_er, this.top_diel_tand
            ));
        }

        // Top air/dielectric region
        dielectrics.push(new Dielectric(
            x_min, this.y_top_start,
            this.domain_width, this.top_dielectric_h,
            this.er_top, this.tand_top
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
                    x_min, this.y_sub_end,
                    this.domain_width, this.sm_t_sub,
                    this.sm_er, this.sm_tand
                ));

                if (this.is_differential) {
                    // Differential: solder mask for both traces
                    const half_spacing = this.trace_spacing / 2;
                    const xl_left = -this.w - half_spacing;
                    const xr_left = -half_spacing;
                    const xl_right = half_spacing;
                    const xr_right = this.w + half_spacing;

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
                    if (xsl >= x_min) {
                        dielectrics.push(new Dielectric(
                            xsl, this.y_trace_start,
                            this.sm_t_side, this.t + this.sm_t_trace,
                            this.sm_er, this.sm_tand
                        ));
                    }

                    // Solder mask on right side of trace
                    const xsr = xr + this.sm_t_side;
                    if (xsr <= x_max) {
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

        // Conductors

        // Bottom ground (beneath everything)
        if (this.t_gnd > 0) {
            conductors.push(new Conductor(
                x_min, -this.t_gnd,
                this.domain_width, this.t_gnd,
                false
            ));
        }

        // Bottom ground plane (above cutout extension)
        if (this.gnd_cut_width === 0) {
            // No cutout - full ground plane
            if (this.y_gnd_bot_end > this.y_gnd_bot_start) {
                conductors.push(new Conductor(
                    x_min, this.y_gnd_bot_start,
                    this.domain_width, this.t_gnd,
                    false
                ));
            }
        } else {
            // With cutout - ground on sides only
            const cut_l = -this.gnd_cut_width / 2;
            const cut_r = this.gnd_cut_width / 2;

            // Left ground
            if (cut_l > x_min) {
                conductors.push(new Conductor(
                    x_min, this.y_gnd_bot_start,
                    cut_l - x_min, this.t_gnd,
                    false
                ));
            }

            // Right ground
            if (cut_r < x_max) {
                conductors.push(new Conductor(
                    cut_r, this.y_gnd_bot_start,
                    x_max - cut_r, this.t_gnd,
                    false
                ));
            }
        }

        // Coplanar vias through substrate (if enabled)
        if (this.use_coplanar_gnd && this.use_vias) {
            // Vias should extend from bottom ground to top coplanar grounds
            // Start from top of bottom ground (y_ext_start = t_gnd)
            // End at top of coplanar grounds (max of trace start/end for negative thickness)
            const via_y_start = this.y_ext_start;
            const via_y_end = Math.max(this.y_trace_start, this.y_trace_end);
            const via_height = via_y_end - this.y_ext_start;

            // Left via (from inner edge to left boundary)
            if (this.via_x_left_inner > x_min) {
                conductors.push(new Conductor(
                    x_min, via_y_start,
                    this.via_x_left_inner - x_min, via_height,
                    false
                ));
            }

            // Right via (from inner edge to right boundary)
            if (this.via_x_right_inner < x_max) {
                conductors.push(new Conductor(
                    this.via_x_right_inner, via_y_start,
                    x_max - this.via_x_right_inner, via_height,
                    false
                ));
            }
        }

        // Signal trace(s)
        if (this.is_differential) {
            // Left trace (negative in odd mode, polarity = -1)
            const xl_left = -this.w - this.trace_spacing / 2;
            conductors.push(new Conductor(
                xl_left, this.y_trace_start,
                this.w, this.t,
                true, -1
            ));
            // Right trace (positive in odd mode, polarity = +1)
            const xl_right = this.trace_spacing / 2;
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
                    x_min, this.y_trace_start,
                    this.x_gap_outer_l - x_min, this.t,
                    false
                ));

                // Right top ground (from outer gap edge to right edge)
                conductors.push(new Conductor(
                    this.x_gap_outer_r, this.y_trace_start,
                    x_max - this.x_gap_outer_r, this.t,
                    false
                ));
            } else {
                // Single-ended: grounds on both sides of the trace
                // Left top ground (from left edge to gap edge)
                conductors.push(new Conductor(
                    x_min, this.y_trace_start,
                    this.x_gap_l - x_min, this.t,
                    false
                ));

                // Right top ground (from gap edge to right edge)
                conductors.push(new Conductor(
                    this.x_gap_r, this.y_trace_start,
                    x_max - this.x_gap_r, this.t,
                    false
                ));
            }
        }

        // Top ground plane (if present, for stripline)
        if (this.has_top_gnd) {
            conductors.push(new Conductor(
                x_min, this.y_gnd_top_start,
                this.domain_width, this.t_gnd,
                false
            ));
        }

        // Side ground planes (based on boundaries array)
        const side_gnd_thickness = this.t_gnd;
        const side_gnd_height = this.has_top_gnd ?
            (this.y_gnd_top_start + this.t_gnd) :
            (this.y_top_start + this.top_dielectric_h);

        // Left side ground (if boundaries[0] === "gnd")
        if (this.boundaries[0] === "gnd") {
            conductors.push(new Conductor(
                x_min, -this.t_gnd,
                side_gnd_thickness, side_gnd_height + this.t_gnd,
                false
            ));
        }

        // Right side ground (if boundaries[1] === "gnd")
        if (this.boundaries[1] === "gnd") {
            conductors.push(new Conductor(
                x_max - side_gnd_thickness, -this.t_gnd,
                side_gnd_thickness, side_gnd_height + this.t_gnd,
                false
            ));
        }

        return [dielectrics, conductors];
    }

    _add_coplanar_solder_mask(dielectrics) {
        // Coplanar solder mask: covers gaps and tops of conductors

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
            const x_min = -this.domain_width / 2;
            const x_max = this.domain_width / 2;
            dielectrics.push(new Dielectric(
                x_min, this.y_trace_end,
                xl_gap - x_min, this.sm_t_trace,
                this.sm_er, this.sm_tand
            ));
            dielectrics.push(new Dielectric(
                xr_gap, this.y_trace_end,
                x_max - xr_gap, this.sm_t_trace,
                this.sm_er, this.sm_tand
            ));
        } else {
            // Single-ended coplanar solder mask
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
            if (xsl > -this.domain_width/2) {
                dielectrics.push(new Dielectric(
                    xsl, this.y_trace_start,
                    this.sm_t_side, this.t + this.sm_t_trace,
                    this.sm_er, this.sm_tand
                ));
            }

            // Solder mask on right side of trace
            if (xsr <= this.domain_width/2) {
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
            const x_min = -this.domain_width / 2;
            const x_max = this.domain_width / 2;
            dielectrics.push(new Dielectric(
                x_min, this.y_trace_end,
                xl_gap - x_min, this.sm_t_trace,
                this.sm_er, this.sm_tand
            ));

            // Solder mask on top of right ground
            dielectrics.push(new Dielectric(
                xr_gap, this.y_trace_end,
                x_max - xr_gap, this.sm_t_trace,
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

        // Initialize mask and material arrays
        this.epsilon_r = Array(ny).fill().map(() => new Float64Array(nx).fill(1));
        this.tand = Array(ny).fill().map(() => new Float64Array(nx).fill(1));
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

        // Apply conductors
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
