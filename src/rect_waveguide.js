import { FieldSolver2D, CONSTANTS } from './field_solver.js';
import { Dielectric } from './mesher.js';

const C0 = CONSTANTS.C;

/**
 * Hollow rectangular waveguide.
 *
 * A rectangular metal pipe of inner dimensions a (broad wall) x b (narrow wall), filled
 * with a homogeneous dielectric. Single-ended, and the FIRST medium in this solver that
 * is not TEM or quasi-TEM.
 *
 * FULL-WAVE BACKEND ONLY, and for a different reason than coax. Coax cannot use the
 * rectilinear backend because a tensor-product grid staircases a circle; a waveguide
 * meshes fine on any grid, but the quasi-static backend solves Laplace for a TEM mode
 * and a hollow guide HAS no TEM mode — there is no second conductor, so no potential
 * difference and no static drive. The constructor throws rather than returning zeros.
 *
 * Geometry model — the meshed domain IS the guide interior:
 *
 *   • the domain is the rectangle x in [-a/2, a/2], y in [0, b];
 *   • it is filled by one homogeneous dielectric;
 *   • there are NO conductors. The walls live entirely in `boundaries` (all 'gnd'),
 *     which _clipDomain turns into wallPEC on all four sides, isEdgePEC turns into
 *     Dirichlet edges, and buildLossEdges turns into conductor-loss surfaces.
 *
 * That last point is what makes this medium source-free, and it is why the backend needs
 * a waveguide path at all: everything downstream of _prepareStatic assumes a driven
 * conductor. See TriBackend._prepareWaveguide / waveguideEigen.
 *
 * ONLY THE FUNDAMENTAL MODE IS SOLVED — TE10 when a >= b, TE01 otherwise. Above the
 * second cutoff the guide is over-moded and the higher modes are simply not in the
 * reported results (the Modes tab can show them). A warning is raised there.
 *
 * Exact closed forms this geometry must reproduce (used by tests/test_rect_waveguide.js),
 * with L = max(a, b), kc = pi/L, k0 = 2*pi*f/c0, eta = eta0/sqrt(er):
 *   fc      = c0*kc/(2*pi*sqrt(er))
 *   beta    = sqrt(k0^2*er - kc^2)
 *   eps_eff = er - (kc/k0)^2
 *   Z_TE    = k0*eta0/beta
 *   Zpv     = 2*(b/a)*Z_TE
 *   a_d     = k0^2*er*tand/(2*beta)
 *   a_c     = Rs*(1 + 2*(b/a)*(fc/f)^2) / (b*eta*sqrt(1 - (fc/f)^2))     [TE10]
 */
class RectWaveguideSolver extends FieldSolver2D {
    constructor(options) {
        super();

        this._validate_parameters(options);

        // Inner dimensions. `a` is the broad wall by convention, but nothing requires
        // a >= b: max(a, b) picks the fundamental either way (TE10 or TE01).
        this.a = options.width;
        this.b = options.height;

        this.epsilon_r = options.epsilon_r ?? 1.0;
        this.tan_delta = options.tan_delta ?? 0;
        this.sigma_cond = options.sigma_cond ?? 5.8e7;
        this.rq = options.rq ?? 0;
        this.freq = options.freq ?? 10e9;

        // Full-wave only. The UI locks the solver dropdown and buildSolverFromParams
        // overrides it; this is the backstop for a hand-edited link or an API caller.
        this.mesh_backend = options.mesh_backend ?? 'triangular';
        if (this.mesh_backend !== 'triangular') {
            throw new Error(
                'Rectangular waveguides require the full-wave solver — the quasi-static ' +
                'backend solves Laplace for a TEM mode, and a hollow waveguide has no TEM ' +
                'mode (no second conductor, so no static drive).');
        }

        this.is_differential = false;
        // Fully enclosed: every wall is PEC. modeConfig() then emits an empty abc, so the
        // eigensolve stays all-real (the fast Arnoldi path) and never escalates to the
        // radiating ABC — a closed structure has no open wall.
        this.boundaries = ['gnd', 'gnd', 'gnd', 'gnd'];

        // --- Analytic reference -------------------------------------------------
        // kc is purely GEOMETRIC: it does not depend on er, so causal materials can move
        // eps_r across the sweep without ever invalidating it. That is what lets the
        // backend solve the eigenproblem ONCE and propagate beta(f) analytically.
        const L = Math.max(this.a, this.b);
        this.kc_analytic = Math.PI / L;
        this.fc = C0 * this.kc_analytic / (2 * Math.PI * Math.sqrt(this.epsilon_r));
        // Second cutoff: TE20 (2*pi/a) or TE01 (pi/b), whichever is lower. Above it the
        // guide is over-moded and this solver still reports only the fundamental.
        this.kc2_analytic = Math.min(2 * Math.PI / this.a, Math.PI / this.b);
        this.fc2 = C0 * this.kc2_analytic / (2 * Math.PI * Math.sqrt(this.epsilon_r));

        // --- Backend dispatch flags ---------------------------------------------
        this.mode_type = 'waveguide';    // TriBackend switches on this
        // Half-domain symmetry would be VALID here (the fundamental's symmetry plane is a
        // magnetic wall), but resampleModeField has no parity/mirroring support the way
        // resampleStatic does, so the plotted mode field would be blank for x < 0. An
        // empty box meshes cheaply, so the full domain costs little.
        this.tri_symmetry = false;
        this.has_potential = false;      // no static potential to plot
        this.allow_dc = false;           // a guide is a DC block; the f=0 limit is degenerate
        // The shunt C crosses zero at cutoff, which blows up the interpolator's relative
        // error test — and after the single eigensolve every sweep point is analytic
        // anyway, so there is nothing to save.
        this.allow_interp_sweep = false;

        // --- Geometry ------------------------------------------------------------
        // Domain (tri_backend derives it as x in [-w/2, w/2], y in [-t_gnd, h]).
        this.domain_width = this.a;
        this.domain_height = this.b;
        this.t_gnd = 0;

        this.dielectrics = [
            new Dielectric(-this.a / 2, 0, this.a, this.b, this.epsilon_r, this.tan_delta),
        ];
        // No conductors: the walls are boundary conditions, not meshed bodies.
        this.conductors = [];

        // COSMETIC ONLY — the walls have no physical thickness in this model (they are a
        // lossy surface impedance on the domain boundary, contributing surface loss but no
        // DC resistance), so nothing but the plot reads this. Without it the geometry view
        // would show a bare rectangle of dielectric with no indication of the metal that
        // makes it a waveguide. Same role as CoaxSolver's Rd = 1.10*b shield extent.
        const tWall = 0.06 * Math.min(this.a, this.b);
        this.enclosure_walls = {
            x_min: -this.a / 2, x_max: this.a / 2, y_min: 0, y_max: this.b,
            thickness: tWall,
        };

        // --- Plating --------------------------------------------------------------
        // The wall is ONE continuous surface, so the per-face top/sides/bottom model has
        // nothing to select between: any enabled face means "plate the whole perimeter".
        // `all: true` is what the surface-impedance classifiers key on. thick_corners is
        // dropped — the interior corners are 90-degree wedges with no field singularity.
        const pl = options.plating;
        const platingOn = pl && pl.sigma > 0 && (pl.all || pl.top || pl.sides || pl.bottom);
        this.plating = platingOn
            ? { ...pl, top: true, sides: true, bottom: true, all: true, thick_corners: false }
            : null;

        // Element sizing. With no conductor curves the mesher installs no background size
        // field (buildOccMeshFromGeometry only builds one from conductor boundary curves),
        // so hFine and hCoarse are set EQUAL and the box meshes uniformly. That is the right
        // shape here: there is no feature to grade towards, only a half-sine to resolve.
        //
        // Measured on WR-90 in air through solve_adaptive() — h = min(a,b)/div:
        //   div  tris  build   per-point   kc err     eps_eff err   alpha_c
        //     8   528  1040ms      37ms   -0.0186%       0.0439%    0.10841
        //    12  1169  1178ms      51ms   -0.0083%       0.0283%    0.10841
        //    16  2042  1716ms      91ms   -0.0047%       0.0228%    0.10839
        //    24  4633  3979ms     228ms   -0.0021%       0.0189%    0.10839
        //
        // Unlike every other medium here, alpha_c is NOT the mesh-limited quantity: the
        // walls are straight and the mode is a smooth sinusoid, so there is no corner
        // singularity and no polygonized curve — alpha_c is converged to five figures on
        // the coarsest mesh tried. kc is the one that moves, and it is already at 1e-4 by
        // div=12. Past that the cost roughly doubles per step for a few 1e-5 of kc.
        const h = Math.min(this.a, this.b) / 12;
        this.tri_mesh_hints = { hFine: h, hCoarse: h };

        // w/t are read by generic code (mesh sizing fallbacks, plot axis fallback).
        this.w = this.a;
        this.t = this.b;

        this.x = null;
        this.y = null;
        this.dx = null;
        this.dy = null;
        this.mesh_generated = false;
    }

    // Non-fatal conditions worth telling the user about, evaluated at `freq`. Returned as
    // strings so app_solver can log them alongside openBoundaryWarnings().
    // The over-moded case also stands in for the electrically-large mesh guard, which
    // _check_meshability skips entirely when there are no conductors.
    waveguideWarnings(freq = this.freq) {
        const out = [];
        const g = (f) => (f / 1e9).toFixed(3) + ' GHz';
        if (freq <= this.fc) {
            out.push(`Frequency ${g(freq)} is at or below the ${g(this.fc)} cutoff — ` +
                `the mode is evanescent and carries no power.`);
        } else if (freq >= this.fc2) {
            out.push(`Frequency ${g(freq)} is at or above the ${g(this.fc2)} second cutoff — ` +
                `the guide is over-moded. Only the fundamental mode is reported; use the ` +
                `Modes tab to inspect the higher-order modes.`);
        }
        return out;
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

        positive(options.width, 'width');
        positive(options.height, 'height');
        if (isNum(options.width) && isNum(options.height) &&
            options.width > 0 && options.height > 0) {
            const ratio = Math.max(options.width, options.height) /
                Math.min(options.width, options.height);
            // Past this the uniform mesh needs an impractical number of elements across
            // the long dimension to keep the short one resolved.
            if (ratio > 100) {
                errors.push(`width / height aspect ratio ${ratio.toFixed(1)} is too extreme ` +
                    `to mesh (need <= 100)`);
            }
        }
        if (options.epsilon_r != null) {
            if (!isNum(options.epsilon_r)) errors.push(`epsilon_r must be a valid number (got ${options.epsilon_r})`);
            else if (options.epsilon_r < 1) errors.push(`epsilon_r must be at least 1, got ${options.epsilon_r}`);
        }
        nonneg(options.tan_delta, 'tan_delta');
        nonneg(options.rq, 'rq');
        if (options.sigma_cond != null) positive(options.sigma_cond, 'sigma_cond');

        if (errors.length > 0) {
            throw new Error('Parameter validation failed:\n' + errors.map(e => '  - ' + e).join('\n'));
        }
    }

    // The triangular backend builds and owns its own mesh in TriBackend, and it is the
    // only backend a waveguide supports, so there is nothing to do here (and no Mesher).
    ensure_mesh() {
        return;
    }
}

export { RectWaveguideSolver };
