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

    const { Z0, eps_eff, C, C0 } = await solver.calculate_parameters();
    
    // Compute fields for loss calculations
    solver.compute_fields();

    // Calculate losses
    const alpha_cond = solver.calculate_conductor_loss(solver.Ex, solver.Ey, Z0);
    const alpha_diel = solver.calculate_dielectric_loss(solver.Ex, solver.Ey, Z0);
    const alpha_total = alpha_cond + alpha_diel;

    // Calculate RLGC and complex Z0
    const { Zc, rlgc, eps_eff_mode } = solver.rlgc(alpha_cond, alpha_diel, C, Z0);

    console.log(`Z (complex) = ${Zc.toString()} ohm, eps_eff ${eps_eff_mode.toFixed(3)}, RLGC {'R': ${rlgc.R.toExponential(3)}, 'L': ${rlgc.L.toExponential(3)}, 'G': ${rlgc.G.toExponential(3)}, 'C': ${rlgc.C.toExponential(3)}}`);

    console.log(`\n==================================================`);
    console.log(`MICROSTRIP ANALYSIS RESULTS`);
    console.log(`==================================================`);
    console.log(`Characteristic Impedance Z0:  ${Z0.toFixed(2)} Ω`);
    console.log(`Effective Permittivity εᵣₑff: ${eps_eff.toFixed(3)}`);
    console.log(`Losses (dB/m) @ 1GHz:         Diel=${alpha_diel.toFixed(4)}, Cond=${alpha_cond.toPrecision(4)}`);
    console.log(`Total Attenuation:            ${alpha_total.toFixed(4)} dB/m`);
    console.log(`==================================================\n`);
}

solve_microstrip_js_test();