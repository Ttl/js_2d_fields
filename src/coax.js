import { FieldSolver2D } from './field_solver.js';
import { Dielectric, Conductor } from './mesher.js';
import { shapeBBox } from './shapes.js';

/**
 * Coaxial transmission line.
 *
 * A round centre conductor of diameter d inside a round shield of inner diameter D,
 * with a homogeneous dielectric filling the annulus. Single-ended, always TEM.
 *
 * FULL-WAVE BACKEND ONLY. The rectilinear quasi-static backend cannot represent this:
 * its tensor-product grid mesher and point-in-box material rasterization would
 * staircase both circles. The constructor throws rather than silently producing a
 * square "coax".
 *
 * Geometry model — the meshed domain IS the dielectric disk:
 *
 *   • the centre conductor is a solid disk of radius a = d/2 (meshed, PEC);
 *   • the dielectric fills the disk of radius b = D/2 (this is the whole domain);
 *   • the SHIELD is the complement of that disk — everything at radius >= b. It has
 *     zero meshed area, so no shield metal is meshed and no dead air cavities exist
 *     outside it, but it still owns every node and edge on the outer boundary, which
 *     is what makes that boundary PEC and gives it conductor-loss edges.
 *
 * Modelling the shield this way (rather than as a ring inside a bounding box) is why
 * the mesh is minimal: only the annulus and the centre disk carry triangles. The cost
 * is that the shield has no thickness, so it contributes no DC resistance — accurate
 * at any frequency where the skin depth is below the real shield wall, which is the
 * entire useful range for a coaxial cable.
 *
 * Plating selects a conductor (options.plating.inner / .outer), not a face, see the
 * plating block in the constructor.
 *
 * Exact closed forms this geometry must reproduce (used by tests/test_coax.js):
 *   Z0     = eta0/(2*pi*sqrt(er)) * ln(b/a)
 *   C      = 2*pi*eps0*er / ln(b/a)
 *   L      = mu0/(2*pi) * ln(b/a)
 *   eps_eff= er exactly (homogeneous fill)
 *   a_d    = pi*f*sqrt(er)/c0 * tand
 *   a_c    = Rs*(1/a + 1/b) / (2*eta*ln(b/a)),  eta = eta0/sqrt(er)
 *   with each conductor's Rs set by its own surface: plating the inner conductor scales
 *   the 1/a term, plating the shield scales the 1/b term.
 */
class CoaxSolver extends FieldSolver2D {
    constructor(options) {
        super();

        this._validate_parameters(options);

        // Datasheets specify coax by DIAMETER, so that is what the UI collects; the
        // solver works in radii throughout.
        this.d = options.inner_diameter;
        this.D = options.dielectric_diameter;
        this.a = this.d / 2;
        this.b = this.D / 2;

        this.epsilon_r = options.epsilon_r;
        this.tan_delta = options.tan_delta ?? 0;
        this.sigma_cond = options.sigma_cond ?? 5.8e7;
        this.rq = options.rq ?? 0;
        this.freq = options.freq ?? 1e9;

        // Full-wave only. The UI locks the solver dropdown for coax; this is the
        // programmatic backstop for a hand-edited link or a direct API caller.
        this.mesh_backend = options.mesh_backend ?? 'triangular';
        if (this.mesh_backend !== 'triangular') {
            throw new Error(
                'Coaxial lines require the full-wave solver — the quasi-static backend ' +
                'meshes on a rectilinear grid and cannot represent circular conductors.');
        }

        this.is_differential = false;
        // Fully enclosed by the shield: every wall is PEC. modeConfig() then emits an
        // empty abc, so the eigensolve stays all-real (the fast Arnoldi path) and never
        // escalates to the radiating ABC — a closed structure has no open wall.
        this.boundaries = ['gnd', 'gnd', 'gnd', 'gnd'];

        // --- Circle discretization ---------------------------------------------
        // Both circles are materialized as regular polygons (see shapes.js for why
        // that is the geometry rather than an approximation of it). N is decided ONCE,
        // here, and never changes afterwards: the same polygon is consumed by the OCC
        // mesher, tagMaterials, buildTriFreedomMap, the loss integrals and the causal
        // re-tag, and if any of them materialized a different N they would disagree
        // about where the metal is — e.g. tagMaterials would label the slivers near
        // each polygon vertex as air.
        //
        // N does NOT limit mesh resolution: gmsh subdivides each polygon side as finely
        // as the background size field asks. It is purely geometric fidelity, plus a
        // floor on the boundary node count, so it is tied to the conductor-surface
        // element size only to avoid forcing one element per side.
        // The CEILINGS are set by geometric fidelity, not by hSurf, and they matter more
        // than they look: every polygon side forces at least one mesh edge, so n_inner +
        // n_outer is a hard FLOOR on the element count. Ceilings of 192/256 put that
        // floor near 800 triangles, which a small "Max Nodes" setting then cannot get
        // under — the triangle-budget loop coarsens repeatedly, never succeeds, and the
        // solve ends up slower than with a generous budget.
        //
        // An area-matched n-gon's perimeter error is pi^2/(6 n^2): 4e-4 at n=64, 1e-4 at
        // n=128. Since alpha_c itself is only good to ~1%, anything past ~128 is far
        // below the noise floor and buys nothing but forced edges.
        const hFine = this._hFine();
        const hSurf = 0.35 * hFine;                 // matches occSurfScale in occ_to_mesh
        this.n_inner = clamp4(2 * Math.PI * this.a / hSurf, 32, 128);
        this.n_outer = clamp4(2 * Math.PI * this.b / hSurf, 32, 128);

        // phase 0 with n % 4 === 0 puts vertices exactly on both axes, which is what
        // makes the x >= 0 half an EXACT half of the polygon — required by the
        // half-domain symmetry solve.
        const outer = { type: 'circle', cx: 0, cy: 0, r: this.b, n: this.n_outer, phase: 0, xSymmetric: true };
        const inner = { type: 'circle', cx: 0, cy: 0, r: this.a, n: this.n_inner, phase: 0, xSymmetric: true };
        // The shield reuses the outer polygon so the domain outline and the shield
        // boundary are the same vertices by construction, not merely by agreement.
        const shield = { type: 'outside_circle', cx: 0, cy: 0, r: this.b, n: this.n_outer, phase: 0, xSymmetric: true };

        this.domain_shape = outer;

        // --- Plating ------------------------------------------------------------
        // A circle has ONE continuous surface, so the per-face top/sides/bottom model
        // has nothing to select between. What a coax selects between instead is WHICH
        // CONDUCTOR carries the layer, inner or outer conductor. Each selected
        // conductor is plated all the way around: `all: true` is what the
        // surface-impedance classifiers key on, and thick_corners is dropped
        // since there are no corners to wrap.
        //
        // Both conductors share the same plating object (it is read-only downstream,
        // and makePlatingZs caches its impedance by sigma/rq/thickness), so
        // `this.plating` stays "the layer in effect" and the per-conductor
        // assignment below carries the selection.
        const pl = options.plating;
        const named = pl && (pl.inner !== undefined || pl.outer !== undefined);
        const selected = pl && (named ? (pl.inner || pl.outer)
                                      : (pl.all || pl.top || pl.sides || pl.bottom));
        const platingOn = !!(pl && pl.sigma > 0 && selected);
        const platesInner = platingOn && (named ? !!pl.inner : true);
        const platesOuter = platingOn && (named ? !!pl.outer : false);
        this.plating = platingOn
            ? { ...pl, top: true, sides: true, bottom: true, all: true, thick_corners: false,
                inner: platesInner, outer: platesOuter }
            : null;

        // The shield is modelled as a zero-thickness shell of semi-infinite metal (see the
        // header), which is exactly the substrate the layered plating-over-bulk impedance
        // assumes, so plating it needs no extra geometry, only the surface impedance on
        // its boundary edges changes.
        const innerPlating = platesInner ? this.plating : null;
        const outerPlating = platesOuter ? this.plating : null;

        // --- Geometry lists ------------------------------------------------------
        // Bounding boxes are still supplied because bbox-only consumers (plot framing,
        // meshability estimate, open-boundary check) read them; every containment test
        // goes through the shape.
        const ib = shapeBBox(inner);
        // The shield is drawn as a ring out to Rd; that extent is cosmetic and only
        // sets its bbox. The meshed domain stops at b.
        const Rd = this.b * 1.10;

        this.dielectrics = [
            new Dielectric(-this.b, -this.b, 2 * this.b, 2 * this.b,
                this.epsilon_r, this.tan_delta, outer),
        ];
        this.conductors = [
            new Conductor(ib.xmin, ib.ymin, ib.xmax - ib.xmin, ib.ymax - ib.ymin,
                true, 1, innerPlating, inner),
            new Conductor(-Rd, -Rd, 2 * Rd, 2 * Rd, false, 0, outerPlating, shield),
        ];

        // Domain box must strictly ENCLOSE the meshed disk (tri_backend derives its
        // plotting/resample box from it). t_gnd is the distance from y=0 to the box
        // floor in the backend's convention.
        this.domain_width = 2 * Rd;
        this.domain_height = Rd;
        this.t_gnd = Rd;

        // Element sizing for the full-wave mesher. Passed explicitly so coax sizing
        // does not depend on w/t, which have no meaning here.
        this.tri_mesh_hints = { hFine, hCoarse: (this.b - this.a) / 2 };
        // w/t are read by generic code (mesh sizing fallbacks, plot axis fallback).
        this.w = 2 * this.a;
        this.t = 2 * this.a;

        this.x = null;
        this.y = null;
        this.dx = null;
        this.dy = null;
        this.mesh_generated = false;
    }

    // Base element size. tri_backend's own formula is min(t, w/4)*1.5, keyed on a thin
    // rectangular trace — a solid round conductor has no thin dimension, and using it
    // leaves the surfaces badly under-resolved (measured: alpha_c 3% low, because the
    // conductor-loss line integral is driven entirely by surface resolution while Z0
    // and C converge long before it).
    //
    // Two scales actually matter here:
    //   • the annulus width (b - a), which the bulk field varies over;
    //   • the inner circumference, which sets how finely the surface current — and so
    //     the ∮|Ht|² integral behind alpha_c — is sampled.
    //
    // The divisors were picked by measuring, across b/a = 1.5, 3.2 and 10, how far this
    // can be coarsened before anything degrades. alpha_c is the ONLY quantity that
    // constrains the mesh: it is a pure surface integral, so it converges far later than
    // everything else. Z0, C and alpha_d stay within 0.1% on every mesh tried, including
    // ones several times coarser than this, and never enter the choice.
    //
    // Measured through solve_adaptive() at these divisors — alpha_c error / triangles:
    //   b/a = 1.5   -0.55%   1815
    //   b/a = 3.2   -0.88%   2641
    //   b/a = 10    -1.25%   4400
    // Halving them again (a 2x finer mesh) buys only ~0.2 pp and costs ~1.5x the
    // triangles and ~1.6x the wall time; doubling them costs ~0.5 pp and starts pushing
    // alpha_c past 1.5%. This sits at that knee.
    //
    // Two traps when re-tuning:
    //   • alpha_c is NOT monotonic in h. The adaptive refinement decides where surface
    //     resolution actually lands, so neighbouring settings scatter by a few tenths of
    //     a percentage point. Don't fine-tune against a single geometry.
    //   • Measure through solve_adaptive(), not by driving TriBackend directly. The
    //     refinement runs at the solver's reference frequency, so a hand-built backend
    //     at a different frequency converges to a different mesh and flatters the result.
    _hFine() {
        return Math.min((this.b - this.a) / 4, 2 * Math.PI * this.a / 32);
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

        positive(options.inner_diameter, 'inner_diameter');
        positive(options.dielectric_diameter, 'dielectric_diameter');
        if (isNum(options.inner_diameter) && isNum(options.dielectric_diameter)) {
            const ratio = options.dielectric_diameter / options.inner_diameter;
            if (ratio <= 1) {
                errors.push(`dielectric_diameter (${options.dielectric_diameter}) must be greater ` +
                    `than inner_diameter (${options.inner_diameter})`);
            } else if (ratio < 1.02) {
                // Below this the annulus is thinner than the elements the mesher would
                // need, and the fragment collapses instead of failing cleanly.
                errors.push(`dielectric_diameter / inner_diameter = ${ratio.toFixed(4)} is too close ` +
                    `to 1 to mesh — the dielectric annulus would be degenerate (need >= 1.02)`);
            }
        }
        if (!isNum(options.epsilon_r)) errors.push(`epsilon_r must be a valid number (got ${options.epsilon_r})`);
        else if (options.epsilon_r < 1) errors.push(`epsilon_r must be at least 1, got ${options.epsilon_r}`);
        nonneg(options.tan_delta, 'tan_delta');
        nonneg(options.rq, 'rq');
        if (options.sigma_cond != null) positive(options.sigma_cond, 'sigma_cond');

        if (errors.length > 0) {
            throw new Error('Parameter validation failed:\n' + errors.map(e => '  - ' + e).join('\n'));
        }
    }

    // The triangular backend builds and owns its own mesh in TriBackend, and it is the
    // only backend coax supports, so there is nothing to do here (and no Mesher).
    ensure_mesh() {
        return;
    }
}

// Round to a multiple of 4 within [lo, hi]. Multiples of 4 put polygon vertices exactly
// on both axes, which is what makes the x >= 0 half an exact half.
function clamp4(v, lo, hi) {
    const n = Math.round(v / 4) * 4;
    return Math.max(lo, Math.min(hi, n));
}

export { CoaxSolver };
