import { MicrostripSolver } from './microstrip.js';

/**
 * GroundedCPWSolver2D - Thin wrapper around MicrostripSolver for backward compatibility.
 * Maps the old GCPW constructor parameters to the unified MicrostripSolver interface.
 * Kept for test compatibility; app.js uses MicrostripSolver directly.
 */
export class GroundedCPWSolver2D extends MicrostripSolver {
    constructor(options = {}) {
        // Map old GCPW params to new unified MicrostripSolver params
        super({
            substrate_height: options.substrate_height,
            trace_width: options.trace_width,
            trace_thickness: options.trace_thickness,
            trace_spacing: options.trace_spacing ?? null,
            gnd_thickness: options.gnd_thickness ?? 35e-6,
            epsilon_r: options.epsilon_r ?? 4.5,
            epsilon_r_top: options.epsilon_r_top ?? 1,
            tan_delta: options.tan_delta ?? 0.02,
            sigma_diel: options.sigma_diel ?? 0.0,
            sigma_cond: options.sigma_cond ?? 5.8e7,
            freq: options.freq ?? 1e9,
            nx: options.nx ?? 300,
            ny: options.ny ?? 300,
            air_top: options.air_top ?? null,
            air_side: options.air_side ?? null,
            boundaries: options.boundaries ?? ["open", "open", "open", "gnd"],
            // Coplanar-specific (enable by default for GCPW)
            use_coplanar_gnd: true,
            gap: options.gap,
            top_gnd_width: options.top_gnd_width,
            via_gap: options.via_gap,
            use_vias: true,
            // Solder mask
            use_sm: options.use_sm ?? false,
            sm_t_sub: options.sm_t_sub ?? 20e-6,
            sm_t_trace: options.sm_t_trace ?? 20e-6,
            sm_t_side: options.sm_t_side ?? 20e-6,
            sm_er: options.sm_er ?? 3.5,
            sm_tand: options.sm_tand ?? 0.02,
        });
    }
}
