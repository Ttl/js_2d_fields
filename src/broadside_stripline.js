import { FieldSolver2D } from './field_solver.js';
import { Dielectric, Conductor, Mesher } from './mesher.js';

/**
 * Broadside-coupled stripline.
 *
 * Two signal traces stacked vertically inside three dielectric layers
 * (bottom, middle, top), enclosed top and bottom by ground planes.
 * The upper trace can be horizontally offset relative to the lower one.
 *
 * Always differential. Lower trace polarity = -1, upper = +1.
 *
 * Substrate stack (y, growing upward):
 *   bottom ground (-t_gnd .. 0)
 *   bottom dielectric (h_bottom)   top-of-bottom-ground → bottom of lower trace
 *   lower trace (|t|)              grows up when t > 0, down when t < 0
 *   middle dielectric              gap between the two traces
 *   upper trace (|t|)              grows down when t > 0, up when t < 0
 *   top dielectric (h_top)         top of upper trace → bottom of top ground
 *   top ground (t_gnd)
 *
 * h_bottom: top of bottom ground → bottom reference of lower trace (bottom reference = lowest point when t>0)
 * h_middle: bottom of lower trace → top of upper trace (includes both trace thicknesses + gap)
 * h_top:    top of upper trace → bottom of top ground
 * Ground-to-ground spacing = h_bottom + h_middle + h_top (trace thickness excluded)
 *
 * Negative t reverses both conductor directions (lower grows down, upper grows up).
 */
class BroadsideStriplineSolver extends FieldSolver2D {
    constructor(options) {
        super();

        this._validate_parameters(options);

        this.w = options.trace_width;
        this.t = options.trace_thickness;
        this.t_gnd = options.gnd_thickness ?? 35e-6;

        this.h_bottom = options.h_bottom;
        this.h_middle = options.h_middle;
        this.h_top = options.h_top;

        this.er_bottom = options.er_bottom;
        this.er_middle = options.er_middle;
        this.er_top = options.er_top;

        this.tand_bottom = options.tand_bottom ?? 0;
        this.tand_middle = options.tand_middle ?? 0;
        this.tand_top = options.tand_top ?? 0;

        this.sigma_cond = options.sigma_cond ?? 5.8e7;
        this.x_offset = options.x_offset ?? 0;

        this.freq = options.freq ?? 1e9;
        this.nx = options.nx ?? 300;
        this.ny = options.ny ?? 300;

        this.rq = options.rq ?? 0;
        this.plating = options.plating ?? null;

        // Always differential — base FieldSolver2D mode decomposition relies on this.
        this.is_differential = true;

        // Boundaries: top/bottom always grounded by the ground planes.
        // Sides default open; enclosure makes them ground.
        this.boundaries = options.boundaries ?? ["open", "open", "gnd", "gnd"];
        this.has_side_gnd = (this.boundaries[0] === "gnd" || this.boundaries[1] === "gnd");

        // Domain width
        const total_substrate_h = this.h_bottom + this.h_middle + this.h_top;
        const trace_extent = this.w + Math.abs(this.x_offset);
        if (options.enclosure_width != null && options.enclosure_width !== "auto") {
            this.enclosure_width = options.enclosure_width;
            if (this.has_side_gnd) {
                this.domain_width = this.enclosure_width + 2 * this.t_gnd;
            } else {
                this.domain_width = this.enclosure_width;
            }
        } else {
            this.enclosure_width = null;
            this.domain_width = 2 * Math.max(trace_extent * 8, total_substrate_h * 4);
        }

        this._calculate_coordinates();
        const [dielectrics, conductors] = this._build_geometry_lists();
        this.dielectrics = dielectrics;
        this.conductors = conductors;

        // Symmetric mesh only when x_offset is zero — otherwise geometry isn't mirror-symmetric.
        const symmetric = (this.x_offset === 0);

        this.mesher = new Mesher(
            this.domain_width, this.domain_height,
            this.nx, this.ny,
            this.conductors, this.dielectrics,
            symmetric,
            -this.domain_width / 2,
            this.domain_width / 2,
            -this.t_gnd,
            this.domain_height
        );

        this.x = null;
        this.y = null;
        this.dx = null;
        this.dy = null;
        this.mesh_generated = false;
    }

    _validate_parameters(options) {
        const errors = [];
        const isNum = (v) => typeof v === 'number' && !isNaN(v) && isFinite(v);
        const positive = (v, name) => {
            if (!isNum(v)) errors.push(`${name} must be a valid number (got ${v})`);
            else if (v <= 0) errors.push(`${name} must be positive, got ${v}`);
        };
        const nonneg = (v, name) => {
            if (v == null) return;
            if (!isNum(v)) errors.push(`${name} must be a valid number (got ${v})`);
            else if (v < 0) errors.push(`${name} must be non-negative, got ${v}`);
        };
        positive(options.trace_width, 'trace_width');
        const nonzero = (v, name) => {
            if (!isNum(v)) errors.push(`${name} must be a valid number (got ${v})`);
            else if (v === 0) errors.push(`${name} must be non-zero`);
        };
        nonzero(options.trace_thickness, 'trace_thickness');
        positive(options.h_bottom, 'h_bottom');
        positive(options.h_middle, 'h_middle');
        positive(options.h_top, 'h_top');
        if (isNum(options.trace_thickness) && isNum(options.h_middle) && isNum(options.h_bottom) && isNum(options.h_top)) {
            const t = options.trace_thickness;
            if (t > 0) {
                if (options.h_middle <= 2 * t)
                    errors.push(`h_middle (${options.h_middle}) must be greater than 2 × trace_thickness (${2 * t}) — conductors would collide`);
            } else {
                const abs_t = Math.abs(t);
                if (options.h_bottom <= abs_t)
                    errors.push(`h_bottom (${options.h_bottom}) must be greater than |trace_thickness| (${abs_t}) — lower conductor would penetrate bottom ground`);
                if (options.h_top <= abs_t)
                    errors.push(`h_top (${options.h_top}) must be greater than |trace_thickness| (${abs_t}) — upper conductor would penetrate top ground`);
            }
        }
        positive(options.er_bottom, 'er_bottom');
        positive(options.er_middle, 'er_middle');
        positive(options.er_top, 'er_top');
        nonneg(options.tand_bottom, 'tand_bottom');
        nonneg(options.tand_middle, 'tand_middle');
        nonneg(options.tand_top, 'tand_top');
        if (options.x_offset != null && !isNum(options.x_offset)) {
            errors.push(`x_offset must be a valid number (got ${options.x_offset})`);
        }
        if (errors.length > 0) {
            throw new Error('Parameter validation failed:\n' + errors.map(e => '  - ' + e).join('\n'));
        }
    }

    _calculate_coordinates() {
        const t = this.t;

        // Lower trace: reference at h_bottom (bottom of lower trace for t>0).
        // Positive t → grows upward; negative t → grows downward.
        this.y_lower_trace_start = this.h_bottom + Math.min(t, 0);
        this.y_lower_trace_end   = this.h_bottom + Math.max(t, 0);

        // Upper trace: reference at h_bottom + h_middle (top of upper trace for t>0).
        // Positive t → grows downward; negative t → grows upward.
        const y_upper_ref = this.h_bottom + this.h_middle;
        this.y_upper_trace_start = y_upper_ref - Math.max(t, 0);
        this.y_upper_trace_end   = y_upper_ref - Math.min(t, 0);

        // Dielectric layers aligned to substrate definitions (independent of trace direction).
        this.y_bot_diel_start = 0;
        this.y_bot_diel_end   = this.h_bottom;
        this.y_mid_diel_start = this.h_bottom;
        this.y_mid_diel_end   = this.h_bottom + this.h_middle;
        this.y_top_diel_start = this.h_bottom + this.h_middle;
        this.y_gnd_top_start  = this.h_bottom + this.h_middle + this.h_top;
        this.y_top_diel_end   = this.y_gnd_top_start;

        this.y_gnd_top_end  = this.y_gnd_top_start + this.t_gnd;
        this.domain_height  = this.y_gnd_top_end;
    }

    _build_geometry_lists() {
        const dielectrics = [];
        const conductors = [];

        const x_min = -this.domain_width / 2;
        const x_max = this.domain_width / 2;

        // Three dielectric layers spanning full width.
        // Cells inside the trace conductors are overwritten by the conductor mask;
        // permittivity values inside are unused since E-field is zero there.
        dielectrics.push(new Dielectric(
            x_min, this.y_bot_diel_start,
            this.domain_width, this.h_bottom,
            this.er_bottom, this.tand_bottom
        ));
        // Middle dielectric spans the full h_middle extent (h_bottom → h_bottom + h_middle).
        dielectrics.push(new Dielectric(
            x_min, this.y_mid_diel_start,
            this.domain_width, this.h_middle,
            this.er_middle, this.tand_middle
        ));
        dielectrics.push(new Dielectric(
            x_min, this.y_top_diel_start,
            this.domain_width, this.h_top,
            this.er_top, this.tand_top
        ));

        // Bottom ground plane
        conductors.push(new Conductor(
            x_min, -this.t_gnd,
            this.domain_width, this.t_gnd,
            false
        ));

        // Top ground plane
        conductors.push(new Conductor(
            x_min, this.y_gnd_top_start,
            this.domain_width, this.t_gnd,
            false
        ));

        const abs_t = Math.abs(this.t);

        // Lower trace (negative polarity)
        const xl_lower = -this.w / 2;
        conductors.push(new Conductor(
            xl_lower, this.y_lower_trace_start,
            this.w, abs_t,
            true, -1, this.plating
        ));

        // Upper trace (positive polarity), shifted by x_offset
        const xl_upper = -this.w / 2 + this.x_offset;
        conductors.push(new Conductor(
            xl_upper, this.y_upper_trace_start,
            this.w, abs_t,
            true, 1, this.plating
        ));

        // Side ground walls if enclosure enabled
        const side_full_height = this.y_gnd_top_end + this.t_gnd;
        if (this.boundaries[0] === "gnd") {
            conductors.push(new Conductor(
                x_min, -this.t_gnd,
                this.t_gnd, side_full_height,
                false
            ));
        }
        if (this.boundaries[1] === "gnd") {
            conductors.push(new Conductor(
                x_max - this.t_gnd, -this.t_gnd,
                this.t_gnd, side_full_height,
                false
            ));
        }

        return [dielectrics, conductors];
    }

    ensure_mesh() {
        if (this.mesh_generated) return;

        [this.x, this.y] = this.mesher.generate_mesh();

        this.dx = new Float64Array(this.x.length - 1);
        for (let i = 0; i < this.x.length - 1; i++) this.dx[i] = this.x[i + 1] - this.x[i];

        this.dy = new Float64Array(this.y.length - 1);
        for (let i = 0; i < this.y.length - 1; i++) this.dy[i] = this.y[i + 1] - this.y[i];

        this._setup_geometry();
        this.mesh_generated = true;
    }

    _setup_geometry() {
        const tol = 1e-11;
        const nx = this.x.length;
        const ny = this.y.length;

        this.epsilon_r = Array(ny).fill().map(() => new Float64Array(nx).fill(1));
        this.tand = Array(ny).fill().map(() => new Float64Array(nx).fill(0));
        this.signal_mask = Array(ny).fill().map(() => new Uint8Array(nx));
        this.ground_mask = Array(ny).fill().map(() => new Uint8Array(nx));
        this.signal_p_mask = Array(ny).fill().map(() => new Uint8Array(nx));
        this.signal_n_mask = Array(ny).fill().map(() => new Uint8Array(nx));

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
                                if (cond.polarity > 0) this.signal_p_mask[i][j] = 1;
                                else this.signal_n_mask[i][j] = 1;
                            } else {
                                this.ground_mask[i][j] = 1;
                            }
                        }
                    }
                }
            }
        }

        this.conductor_mask = Array(ny).fill().map((_, i) => {
            const row = new Uint8Array(nx);
            for (let j = 0; j < nx; j++) {
                row[j] = this.signal_mask[i][j] | this.ground_mask[i][j];
            }
            return row;
        });
    }
}

export { BroadsideStriplineSolver };
