
const { MicrostripSolver } = require('./microstrip.js');
const TL = { MicrostripSolver };

async function solve_microstrip_js_test() {
    console.log("Solving Microstrip (JS version) for test...");

    const solver = new TL.MicrostripSolver(
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

    console.log(`Characteristic Impedance Z0: ${Z0.toFixed(2)} Ω`);
    console.log(`Effective Permittivity εᵣₑff: ${eps_eff.toFixed(3)}`);
}

solve_microstrip_js_test();