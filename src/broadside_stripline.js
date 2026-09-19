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
        // Numerical backend: 'rectilinear' (FDM, default) or 'triangular' (FEM).
        this.mesh_backend = options.mesh_backend ?? 'rectilinear';

        // Domain width
        const total_substrate_h = this.h_bottom + this.h_middle + this.h_top;
        if (options.enclosure_width != null && options.enclosure_width !== "auto") {
            this.enclosure_width = options.enclosure_width;
            if (this.has_side_gnd) {
                this.domain_width = this.enclosure_width + 2 * this.t_gnd;
            } else {
                this.domain_width = this.enclosure_width;
            }
        } else {
            this.enclosure_width = null;
            // Same trace-to-wall clearance rule as the other line types: a margin set by
            // the trace width and the substrate stack.
            // x_offset only translated the upper trace, so it widens the domain
            // by the translation instead of scaling the margin with it.
            const margin = Math.max(this.w * 8, total_substrate_h * 4);
            this.domain_width = 2 * (margin + Math.abs(this.x_offset));
        }

        this._calculate_coordinates();
        const [dielectrics, conductors] = this._build_geometry_lists();
        this.dielectrics = dielectrics;
        this.conductors = conductors;

        // Symmetric mesh only when x_offset is zero — otherwise geometry isn't mirror-symmetric.
        // No half-domain (sym_half) meshing. A broadside pair traces are on top
        // of each other so the plane can't separate the modes.
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

        // Strong broadside coupling: reduced conductor-loss accuracy on this
        // (rectilinear) backend: the pertrubation surface intergral
        // over-weights the facing surfaces relative to the true current
        // distribution once the traces couple strongly. Measured vs the
        // tri-backend multi-drive MQS reference (2026-08-16, 20 fuzzer specs):
        // rows with w/facing-gap > 1.75 read R +10 to +45% high (growing with
        // coupling). The full trace width is deliberately used even with
        // x_offset, the worst measured case had zero vertical overlap.
        // Threshold lowered 1.75 -> 1.5 after a fuzzer sweep. The bias
        // already reaches +21% at w/gap 1.68 and +26% at 1.70
        // (seeds 1/4), while 1.5-1.6 rows sit ~5-11% (conservative warns).
        const facing_gap = this.h_middle - 2 * Math.abs(this.t);
        this._proximityWarn = (facing_gap > 0 && this.w / facing_gap >= 1.5)
            ? { type: 'accuracy', reason: 'broadside-proximity', mode: 'all', message:
                `Strongly coupled broadside pair (trace width ${(this.w * 1e6).toFixed(0)} µm vs ` +
                `${(facing_gap * 1e6).toFixed(0)} µm facing gap): conductor loss accuracy is reduced. ` +
                `R typically reads up to 50% high in this regime. The full-wave solver models the ` +
                `broadside proximity effect accurately.` }
            : null;
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
        // Both traces must sit strictly inside the enclosure: the lower one is centred,
        // the upper one is shifted by x_offset, and a trace reaching the side wall is
        // shorted to it (or cut by an open truncation).
        if (options.enclosure_width != null && options.enclosure_width !== "auto") {
            positive(options.enclosure_width, 'enclosure_width');
            if (isNum(options.enclosure_width) && isNum(options.trace_width)) {
                const active_width = options.trace_width + 2 * Math.abs(isNum(options.x_offset) ? options.x_offset : 0);
                if (active_width >= options.enclosure_width)
                    errors.push(`Active area width (${(active_width * 1000).toFixed(3)} mm, trace_width + 2 × |x_offset|) must be smaller than the enclosure inner width (${(options.enclosure_width * 1000).toFixed(3)} mm)`);
            }
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
}

export { BroadsideStriplineSolver };
