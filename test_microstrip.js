import { MicrostripSolver } from './microstrip.js';

async function solve_microstrip_js_test() {
    console.log("Solving Microstrip...");

    const solver = new MicrostripSolver({
        substrate_height: 1.6e-3,
        trace_width: 3e-3,
        trace_thickness: 35e-6,
        gnd_thickness: 35e-6,
        epsilon_r: 4.5,
        tan_delta: 0.02,
        sigma_cond: 5.8e7,
        freq: 1e9,
        nx: 100,
        ny: 100,
        skin_cells: 50,
        use_sm: false,
        boundaries: ["open", "open", "open", "gnd"]
    });

    const results = await solver.solve_adaptive();

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

async function solve_stripline_js_test() {
    console.log("Solving Stripline...");

    const solver = new MicrostripSolver({
        substrate_height: 0.2e-3,
        trace_width: 0.15e-3,
        trace_thickness: 35e-6,
        gnd_thickness: 16e-6,
        epsilon_r: 4.1,
        epsilon_r_top: 4.1,
        air_top: 0.2e-3,
        tan_delta: 0.02,
        sigma_cond: 5.8e7,
        freq: 1e9,
        nx: 100,
        ny: 100,
        skin_cells: 50,
        boundaries: ["open", "open", "gnd", "gnd"]
    });

    const results = await solver.solve_adaptive();

    console.log(`Z (complex) = ${results.Zc.toString()} ohm, eps_eff ${results.eps_eff.toFixed(3)}, RLGC {'R': ${results.RLGC.R.toExponential(3)}, 'L': ${results.RLGC.L.toExponential(3)}, 'G': ${results.RLGC.G.toExponential(3)}, 'C': ${results.RLGC.C.toExponential(3)}}`);

    console.log(`\n==================================================`);
    console.log(`STRIPLINE ANALYSIS RESULTS`);
    console.log(`==================================================`);
    console.log(`Characteristic Impedance Z0:  ${results.Z0.toFixed(2)} Ω`);
    console.log(`Effective Permittivity εᵣₑff: ${results.eps_eff.toFixed(3)}`);
    console.log(`Losses (dB/m) @ 1GHz:         Diel=${results.alpha_diel_db_m.toFixed(4)}, Cond=${results.alpha_cond_db_m.toFixed(4)}`);
    console.log(`Total Attenuation:            ${results.total_alpha_db_m.toFixed(4)} dB/m`);
    console.log(`==================================================\n`);
}

async function solve_microstrip_embed_js_test() {
    console.log("Solving Embedded Microstrip...");

    const solver = new MicrostripSolver({
        substrate_height: 1.6e-3,
        trace_width: 3e-3,
        trace_thickness: 35e-6,
        gnd_thickness: 35e-6,
        epsilon_r: 4.5,
        tan_delta: 0.02,
        sigma_cond: 5.8e7,
        freq: 1e9,
        nx: 200,
        ny: 200,
        skin_cells: 50,
        use_sm: false,
        top_diel_h: 0.2e-3,
        top_diel_er: 4.5,
        boundaries: ["open", "open", "open", "gnd"]
    });

    const results = await solver.solve_adaptive();

    console.log(`Z (complex) = ${results.Zc.toString()} ohm, eps_eff ${results.eps_eff.toFixed(3)}, RLGC {'R': ${results.RLGC.R.toExponential(3)}, 'L': ${results.RLGC.L.toExponential(3)}, 'G': ${results.RLGC.G.toExponential(3)}, 'C': ${results.RLGC.C.toExponential(3)}}`);

    console.log(`\n==================================================`);
    console.log(`EMBEDDED MICROSTRIP ANALYSIS RESULTS`);
    console.log(`==================================================`);
    console.log(`Characteristic Impedance Z0:  ${results.Z0.toFixed(2)} Ω`);
    console.log(`Effective Permittivity εᵣₑff: ${results.eps_eff.toFixed(3)}`);
    console.log(`Losses (dB/m) @ 1GHz:         Diel=${results.alpha_diel_db_m.toFixed(4)}, Cond=${results.alpha_cond_db_m.toPrecision(4)}`);
    console.log(`Total Attenuation:            ${results.total_alpha_db_m.toFixed(4)} dB/m`);
    console.log(`==================================================\n`);
}

async function solve_microstrip_cut_js_test() {
    console.log("Solving Microstrip with gnd cut...");

    const solver = new MicrostripSolver({
        substrate_height: 1.6e-3,
        trace_width: 3e-3,
        trace_thickness: 35e-6,
        gnd_thickness: 35e-6,
        epsilon_r: 4.5,
        tan_delta: 0.02,
        sigma_cond: 5.8e7,
        freq: 1e9,
        nx: 300,
        ny: 300,
        skin_cells: 50,
        use_sm: false,
        gnd_cut_width: 3e-3,
        gnd_cut_sub_h: 1e-3,
        boundaries: ["open", "open", "open", "gnd"]
    });

    const results = await solver.solve_adaptive();

    console.log(`Z (complex) = ${results.Zc.toString()} ohm, eps_eff ${results.eps_eff.toFixed(3)}, RLGC {'R': ${results.RLGC.R.toExponential(3)}, 'L': ${results.RLGC.L.toExponential(3)}, 'G': ${results.RLGC.G.toExponential(3)}, 'C': ${results.RLGC.C.toExponential(3)}}`);

    console.log(`\n==================================================`);
    console.log(`MICROSTRIP WITH GND CUT ANALYSIS RESULTS`);
    console.log(`==================================================`);
    console.log(`Characteristic Impedance Z0:  ${results.Z0.toFixed(2)} Ω`);
    console.log(`Effective Permittivity εᵣₑff: ${results.eps_eff.toFixed(3)}`);
    console.log(`Losses (dB/m) @ 1GHz:         Diel=${results.alpha_diel_db_m.toFixed(4)}, Cond=${results.alpha_cond_db_m.toPrecision(4)}`);
    console.log(`Total Attenuation:            ${results.total_alpha_db_m.toFixed(4)} dB/m`);
    console.log(`==================================================\n`);
}

// Run tests
async function runTests() {
    await solve_microstrip_js_test();
    await solve_stripline_js_test();
    await solve_microstrip_embed_js_test();
    await solve_microstrip_cut_js_test();
}

runTests();
