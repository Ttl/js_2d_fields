// Solver construction from a plain params object.
//
// Deliberately DOM-free: this module is imported both by app_solver.js (main thread, to
// build the cheap geometry-preview solver) and by solve_worker.js (worker thread, where
// there is no document at all). Any getElementById here breaks the worker. `p`
// is the plain object produced by app_solver's getParams().
import { MicrostripSolver } from './microstrip.js';
import { BroadsideStriplineSolver } from './broadside_stripline.js';
import { CoaxSolver } from './coax.js';
import { RectWaveguideSolver } from './rect_waveguide.js';

// Line types with no quasi-static implementation. Coax: a rectilinear grid staircases a
// circle. Rectangular waveguide: a hollow guide has no TEM mode at all, so there is
// nothing for the Laplace solve to find. Both constructors throw on any other backend;
// this set is the app-level guard, and field_solver.html keeps its own copy for the
// Solver-dropdown lock (separate script scope).
export const FULLWAVE_ONLY_TYPES = new Set(['coax', 'rect_waveguide']);

export function solverModeConfig(mode) {
    switch (mode) {
        case 'fullwave_pert':
            return { mesh_backend: 'triangular', tri_opts: { lossMethod: 'perturbation' } };
        case 'fullwave_mqs':
        case 'fullwave_occ':   // legacy saved value (the triangular backend now always uses OCC)
        case 'triangular':     // legacy saved value
            return { mesh_backend: 'triangular', tri_opts: { lossMethod: 'auto' } };
        default:             // 'rectilinear'
            return { mesh_backend: 'rectilinear', tri_opts: null };
    }
}

function addCommonOptions(options, p) {
    // Solder mask
    if (p.use_sm) {
        options.use_sm = true;
        options.sm_t_sub = p.sm_t_sub;
        options.sm_t_trace = p.sm_t_trace;
        options.sm_t_side = p.sm_t_side;
        options.sm_er = p.sm_er;
        options.sm_tand = p.sm_tand;
    }

    // Top dielectric
    if (p.use_top_diel) {
        options.top_diel_h = p.top_diel_h;
        options.top_diel_er = p.top_diel_er;
        options.top_diel_tand = p.top_diel_tand;
    }

    // Ground cutout
    if (p.use_gnd_cut) {
        options.gnd_cut_width = p.gnd_cut_w;
        options.gnd_cut_sub_h = p.gnd_cut_h;
    }

    // Enclosure
    if (p.use_enclosure) {
        options.enclosure_width = p.enclosure_width;
        if (options.enclosure_height === undefined) {
            options.enclosure_height = p.enclosure_height;
        }

        const left_bc = p.use_side_gnd ? "gnd" : "open";
        const right_bc = p.use_side_gnd ? "gnd" : "open";
        const top_bc = p.use_top_gnd ? "gnd" : (options.boundaries ? options.boundaries[2] : "open");
        const bottom_bc = options.boundaries ? options.boundaries[3] : "gnd";
        options.boundaries = [left_bc, right_bc, top_bc, bottom_bc];
    }

    // Surface plating
    const plating = platingOptions(p);
    if (plating) options.plating = plating;
}

// The plating block on its own, so line types that take plating but none of the other
// board-stackup options (coax) can reuse it instead of duplicating it. `extra` carries a
// type's own surface selection where top/sides/bottom has no meaning, coax passes
// inner/outer (which conductor), and CoaxSolver reads those in place of the faces.
export function platingOptions(p, extra = null) {
    if (!p.use_plating) return null;
    return {
        sigma: p.plating_sigma,
        thickness: p.plating_t,
        rq: p.plating_rq,
        top: p.plating_top,
        sides: p.plating_sides,
        bottom: p.plating_bottom,
        thick_corners: p.plating_thick_corners,
        ...extra,
    };
}

// Build a FieldSolver from the given sidebar params without touching any global
// state. updateGeometry() uses it to (re)build the main `solver`. The Modes
// tab builds its own independent solver with it, the worker builds its
// solve-side solver with it.
//
// Returns the solver, or null if the parameters are invalid. `onError` receives
// the validation message, the caller decides where it goes (the app logs it,
// the worker forwards it), because this module cannot reach the DOM.
export function buildSolverFromParams(p, onError = null) {
    let solver = null;
    try {
        if (p.tl_type === 'gcpw') {
            const options = {
                substrate_height: p.h,
                trace_width: p.w,
                trace_thickness: p.t,
                gnd_thickness: 35e-6,
                epsilon_r: p.er,
                tan_delta: p.tand,
                sigma_cond: p.sigma,
                freq: p.freq,
                nx: p.nx,
                ny: p.ny,
                boundaries: ["open", "open", "open", "gnd"],
                // Coplanar-specific
                use_coplanar_gnd: true,
                gap: p.gap,
                via_gap: p.via_gap,
                use_vias: true,
                // Surface roughness
                rq: p.rq,
            };
            addCommonOptions(options, p);
            solver = new MicrostripSolver(options);
        } else if (p.tl_type === 'diff_gcpw') {
            const options = {
                substrate_height: p.h,
                trace_width: p.w,
                trace_thickness: p.t,
                trace_spacing: p.trace_spacing,  // Enables differential mode
                gnd_thickness: 35e-6,
                epsilon_r: p.er,
                tan_delta: p.tand,
                sigma_cond: p.sigma,
                freq: p.freq,
                nx: p.nx,
                ny: p.ny,
                boundaries: ["open", "open", "open", "gnd"],
                // Coplanar-specific
                use_coplanar_gnd: true,
                gap: p.gap,
                via_gap: p.via_gap,
                use_vias: true,
                // Surface roughness
                rq: p.rq,
            };
            addCommonOptions(options, p);
            solver = new MicrostripSolver(options);
        } else if (p.tl_type === 'diff_microstrip') {
            // Differential Microstrip
            const options = {
                trace_width: p.w,
                substrate_height: p.h,
                trace_thickness: p.t,
                trace_spacing: p.trace_spacing,  // Enable differential mode
                epsilon_r: p.er,
                tan_delta: p.tand,
                sigma_cond: p.sigma,
                freq: p.freq,
                nx: p.nx,
                ny: p.ny,
                boundaries: ["open", "open", "open", "gnd"],
                // Surface roughness
                rq: p.rq
            };
            addCommonOptions(options, p);
            solver = new MicrostripSolver(options);
        } else if (p.tl_type === 'stripline') {
            const options = {
                trace_width: p.w,
                substrate_height: p.h,
                trace_thickness: p.t,
                epsilon_r: p.er,
                epsilon_r_top: p.er_top,
                tan_delta_top: p.tand_top,
                enclosure_height: p.stripline_top_h,
                tan_delta: p.tand,
                sigma_cond: p.sigma,
                freq: p.freq,
                nx: p.nx,
                ny: p.ny,
                boundaries: ["open", "open", "gnd", "gnd"],
                // Surface roughness
                rq: p.rq
            };
            addCommonOptions(options, p);
            solver = new MicrostripSolver(options);
        } else if (p.tl_type === 'coax') {
            // Full-wave only. CoaxSolver throws on any other backend. addCommonOptions is
            // deliberately not used. Solder mask, top dielectric, ground cutout and the
            // enclosure have no meaning inside a coax.
            // Plating and surface roughness do apply and are passed through. The plating
            // is selected per conductor here rather than per face (each conductor is one
            // continuous circular surface).
            solver = new CoaxSolver({
                inner_diameter: p.coax_d,
                dielectric_diameter: p.coax_D,
                epsilon_r: p.coax_er,
                tan_delta: p.coax_tand,
                sigma_cond: p.coax_sigma,
                rq: p.rq,
                plating: platingOptions(p, { inner: p.coax_plating_inner,
                                             outer: p.coax_plating_outer }),
                freq: p.freq,
                mesh_backend: 'triangular',
            });
        } else if (p.tl_type === 'rect_waveguide') {
            // Full-wave only. No TEM mode for the quasi-static Laplace solve to
            // find. addCommonOptions is deliberately not used, solder mask, top
            // dielectric, ground cutout and the enclosure have no meaning
            // inside a waveguide. Plating and surface roughness do apply and
            // are passed through.
            solver = new RectWaveguideSolver({
                width: p.wg_a,
                height: p.wg_b,
                epsilon_r: p.wg_er,
                tan_delta: p.wg_tand,
                sigma_cond: p.wg_sigma,
                rq: p.rq,
                plating: platingOptions(p),
                freq: p.freq,
                mesh_backend: 'triangular',
            });
        } else if (p.tl_type === 'broadside_stripline') {
            const options = {
                trace_width: p.bs_w,
                trace_thickness: p.bs_t,
                x_offset: p.bs_x_offset,
                sigma_cond: p.bs_sigma,
                h_bottom: p.bs_h_bottom,
                er_bottom: p.bs_er_bottom,
                tand_bottom: p.bs_tand_bottom,
                h_middle: p.bs_h_middle,
                er_middle: p.bs_er_middle,
                tand_middle: p.bs_tand_middle,
                h_top: p.bs_h_top,
                er_top: p.bs_er_top,
                tand_top: p.bs_tand_top,
                freq: p.freq,
                nx: p.nx,
                ny: p.ny,
                rq: p.rq,
                boundaries: ["open", "open", "gnd", "gnd"],
            };
            // Enclosure: only side ground walls apply (top/bottom are intrinsic).
            if (p.use_enclosure) {
                options.enclosure_width = p.enclosure_width;
                if (p.use_side_gnd) {
                    options.boundaries = ["gnd", "gnd", "gnd", "gnd"];
                }
            }
            // Plating
            if (p.use_plating) {
                options.plating = {
                    sigma: p.plating_sigma,
                    thickness: p.plating_t,
                    rq: p.plating_rq,
                    top: p.plating_top,
                    sides: p.plating_sides,
                    bottom: p.plating_bottom,
                    thick_corners: p.plating_thick_corners
                };
            }
            solver = new BroadsideStriplineSolver(options);
        } else if (p.tl_type === 'diff_stripline') {
            const options = {
                trace_width: p.w,
                substrate_height: p.h,
                trace_thickness: p.t,
                trace_spacing: p.trace_spacing,  // Enable differential mode
                epsilon_r: p.er,
                epsilon_r_top: p.er_top,
                enclosure_height: p.stripline_top_h,
                tan_delta: p.tand,
                tan_delta_top: p.tand_top,
                sigma_cond: p.sigma,
                freq: p.freq,
                nx: p.nx,
                ny: p.ny,
                boundaries: ["open", "open", "gnd", "gnd"],
                // Surface roughness
                rq: p.rq
            };
            addCommonOptions(options, p);
            solver = new MicrostripSolver(options);
        } else {
            // Microstrip (with optional solder mask, top dielectric, ground cutout)
            const options = {
                trace_width: p.w,
                substrate_height: p.h,
                trace_thickness: p.t,
                epsilon_r: p.er,
                tan_delta: p.tand,
                sigma_cond: p.sigma,
                freq: p.freq,
                nx: p.nx,
                ny: p.ny,
                boundaries: ["open", "open", "open", "gnd"],
                // Surface roughness
                rq: p.rq
            };
            addCommonOptions(options, p);
            solver = new MicrostripSolver(options);
        }

        // Store causal materials option on solver
        if (solver) {
            solver.use_causal_materials = p.use_causal_materials;
            // Some media have no quasi-static implementation.
            // Force the full-wave backend. The UI disables the option and the
            // constructors throw, but a stale link can still arrive with 'rectilinear'.
            const backend = FULLWAVE_ONLY_TYPES.has(p.tl_type) ? 'fullwave_mqs' : p.mesh_backend;
            // Solver mode = numerical backend + triangular loss method:
            //   'rectilinear'   = quasi-static FDM (fastest)
            //   'fullwave_pert' = triangular full-wave, perturbation loss (~2× faster)
            //   'fullwave_mqs'  = triangular full-wave, MQS volume loss (most accurate)
            const cfg = solverModeConfig(backend);
            solver.mesh_backend = cfg.mesh_backend;
            solver.tri_opts = cfg.tri_opts;
        }
    } catch (error) {
        // Log validation errors to the console
        if (onError) onError('ERROR: ' + error.message);
        // Set solver to null to prevent simulation from running with invalid parameters
        solver = null;
    }
    return solver;
}

