import { MicrostripSolver } from './microstrip.js';

async function solve_microstrip_js_test() {
    console.log("Solving Microstrip (JS version) for test...");

    const solver = new MicrostripSolver(
        3e-3,        // trace_width
        1.6e-3,      // substrate_height
        35e-6,       // trace_thickness
        4.5,         // epsilon_r
        0.02,        // tan_delta
        5.8e7,       // sigma_cond
        1e9,         // freq
        100,         // nx
        100,         // ny
        50           // skin_cells
    );

    const results = await solver.perform_analysis();

    console.log(`Z (complex) = ${results.Zc.toString()} ohm, eps_eff ${results.eps_eff.toFixed(3)}, RLGC {'R': ${results.RLGC.R.toExponential(3)}, 'L': ${results.RLGC.L.toExponential(3)}, 'G': ${results.RLGC.G.toExponential(3)}, 'C': ${results.RLGC.C.toExponential(3)}}`);

    console.log(`\n==================================================`);
    console.log(`MICROSTRIP ANALYSIS RESULTS`);
    console.log(`==================================================`);
    console.log(`Characteristic Impedance Z0:  ${results.Z0.toFixed(2)} Ω`);
    console.log(`Effective Permittivity εᵣₑff: ${results.eps_eff.toFixed(3)}`);
    console.log(`Losses (dB/m) @ 1GHz:         Diel=${results.alpha_diel_db_m.toFixed(4)}, Cond=${results.alpha_cond_db_m.toPrecision(4)}`);
    console.log(`Total Attenuation:            ${results.total_alpha_db_m.toFixed(4)} dB/m`);
    console.log(`==================================================\n`);
}

solve_microstrip_js_test();