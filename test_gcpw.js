import { GroundedCPWSolver2D } from './gcpw.js';

async function test_gcpw_js() {
    console.log("Solving GCPW (JS version)...");
    const solver = new GroundedCPWSolver2D(
        1.6e-3,  // substrate_height
        0.3e-3,  // trace_width
        35e-6,   // trace_thickness
        0.15e-3, // gap
        5e-3,    // top_gnd_width
        0.5e-3,  // via_gap
        0.3e-3,  // via_diameter
        35e-6,   // gnd_thickness
        4.5,     // epsilon_r
        0.02,    // tan_delta
        0.0,     // sigma_diel
        5.8e7,   // sigma_cond
        1,       // epsilon_r_top
        null,    // air_top
        null,    // air_side
        1e9,     // freq
        200,     // nx
        200,     // ny
        50,      // skin_cells
        null,    // boundaries
        false,   // use_sm
        20e-6,   // sm_t_sub
        20e-6,   // sm_t_trace
        20e-6,   // sm_t_side
        3.5,     // sm_er
        0.02     // sm_tand
    );

    const results = await solver.perform_analysis();

    console.log(`Z (complex) = ${results.Zc.toString()} ohm, eps_eff ${results.eps_eff.toFixed(3)}, RLGC {'R': ${results.RLGC.R.toExponential(3)}, 'L': ${results.RLGC.L.toExponential(3)}, 'G': ${results.RLGC.G.toExponential(3)}, 'C': ${results.RLGC.C.toExponential(3)}}`);

    console.log(`
==================================================`);
    console.log(`GCPW ANALYSIS RESULTS`);
    console.log(`==================================================`);
    console.log(`Characteristic Impedance Z0:  ${results.Z0.toFixed(2)} Ω`);
    console.log(`Effective Permittivity εᵣₑff: ${results.eps_eff.toFixed(3)}`);
    console.log(`Losses (dB/m) @ 1GHz:         Diel=${results.alpha_diel_db_m.toFixed(4)}, Cond=${results.alpha_cond_db_m.toPrecision(4)}`);
    console.log(`Total Attenuation:            ${results.total_alpha_db_m.toFixed(4)} dB/m`);
    console.log(`==================================================\n`);
}

test_gcpw_js();
