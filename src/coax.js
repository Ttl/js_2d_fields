import { CustomGeometrySolver } from './custom_geometry.js';
import { paramChecks, throwIfErrors } from './field_solver.js';
import { coaxToGeometryText } from './custom_geometry_text.js';

/**
 * Coaxial transmission line.
 *
 * A round centre conductor of diameter d inside a round shield of inner diameter D,
 * with a homogeneous dielectric filling the annulus. Single-ended, always TEM.
 *
 * FULL-WAVE BACKEND ONLY. The rectilinear quasi-static backend would staircase both
 * circles, so the constructor throws on it.
 *
 * Geometry model: the custom-geometry coax (coaxToGeometryText). Each circle is a
 * regular n-gon of the circle's area: the centre conductor, the dielectric disk and
 * the shield as a ring of shield_thickness around it, in an open domain. The conductors
 * are polygons, so the MQS conductor loss applies and the |J| plot has a skin mesh.
 *
 * Plating selects a conductor (options.plating.inner / .outer), not a face.
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
class CoaxSolver extends CustomGeometrySolver {
    constructor(options) {
        validateCoax(options);
        const mesh_backend = options.mesh_backend ?? 'triangular';
        if (mesh_backend !== 'triangular') {
            throw new Error(
                'Coaxial lines require the full-wave solver. The quasi-static backend ' +
                'meshes on a rectilinear grid and cannot represent circular conductors.');
        }
        const g = coaxModel(options);
        super({
            text: coaxToGeometryText(g, 'm'),
            sigma_cond: options.sigma_cond ?? 5.8e7,
            rq: options.rq ?? 0,
            freq: options.freq ?? 1e9,
            mesh_backend,
            symmetry: options.symmetry,
        });
        Object.assign(this, g);
        this.is_coax = true;
        this.d = 2 * g.a;
        this.D = 2 * g.b;
    }
}

// Vertex count of every n-gon. The n-gons have the circles' areas, which leaves a
// perimeter error of pi^2/(6n^2) = 1e-4 at n = 128, below the loss accuracy.
const N_GON = 128;

// Radii, vertex counts, material and plating of the coax model.
function coaxModel(options) {
    // Datasheets specify coax by diameter, the model works in radii.
    const a = options.inner_diameter / 2, b = options.dielectric_diameter / 2;
    const pl = options.plating;
    const named = pl && (pl.inner !== undefined || pl.outer !== undefined);
    const selected = pl && (named ? (pl.inner || pl.outer)
                                  : (pl.all || pl.top || pl.sides || pl.bottom));
    const platingOn = !!(pl && pl.sigma > 0 && selected);
    return {
        a, b,
        n_inner: N_GON,
        n_outer: N_GON,
        shield_thickness: options.shield_thickness ?? 0.10 * b,
        epsilon_r: options.epsilon_r,
        tan_delta: options.tan_delta ?? 0,
        // Without a conductor selection the layer goes on the centre conductor.
        plating: platingOn
            ? { sigma: pl.sigma, thickness: pl.thickness, rq: pl.rq ?? 0, rq_interface: pl.rq_interface,
                inner: named ? !!pl.inner : true, outer: named ? !!pl.outer : false }
            : null,
    };
}

function validateCoax(options) {
    const errors = [];
    const { isNum, positive, nonneg } = paramChecks(errors);

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

    throwIfErrors(errors);
}

export { CoaxSolver };
