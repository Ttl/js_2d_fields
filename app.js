import { MicrostripSolver } from './microstrip.js';
import { CONSTANTS } from './field_solver.js';
import { makeStreamlineTraceFromConductors } from './streamlines.js';
const Plotly = window.Plotly;

let solver = null;
let stopRequested = false;

function shouldStop() {
    return stopRequested;
}

function log(msg) {
    const c = document.getElementById('console_out');
    c.textContent += msg + "\n";
    c.scrollTop = c.scrollHeight;
}

function nanToNull(input) {
    return Number.isNaN(input) ? null : input;
}

function getParams() {
    return {
        tl_type: document.getElementById('tl_type').value,
        w: parseFloat(document.getElementById('inp_w').value) * 1e-3,
        h: parseFloat(document.getElementById('inp_h').value) * 1e-3,
        t: parseFloat(document.getElementById('inp_t').value) * 1e-6,
        er: parseFloat(document.getElementById('inp_er').value),
        tand: parseFloat(document.getElementById('inp_tand').value),
        sigma: parseFloat(document.getElementById('inp_sigma').value),
        freq: parseFloat(document.getElementById('inp_freq').value) * 1e9,
        nx: 30,  // Fixed initial grid size
        ny: 30,  // Fixed initial grid size
        // Differential parameters
        trace_spacing: parseFloat(document.getElementById('inp_trace_spacing').value) * 1e-3,
        // GCPW specific parameters
        gap: parseFloat(document.getElementById('inp_gap').value) * 1e-3,
        top_gnd_w: parseFloat(document.getElementById('inp_top_gnd_w').value) * 1e-3,
        via_gap: parseFloat(document.getElementById('inp_via_gap').value) * 1e-3,
        // Stripline parameters
        air_top: parseFloat(document.getElementById('inp_air_top').value) * 1e-3,
        er_top: parseFloat(document.getElementById('inp_er_top').value),
        // Solder mask parameters
        use_sm: document.getElementById('chk_solder_mask').checked,
        sm_t_sub: parseFloat(document.getElementById('inp_sm_t_sub').value) * 1e-6,
        sm_t_trace: parseFloat(document.getElementById('inp_sm_t_trace').value) * 1e-6,
        sm_t_side: parseFloat(document.getElementById('inp_sm_t_side').value) * 1e-6,
        sm_er: parseFloat(document.getElementById('inp_sm_er').value),
        sm_tand: parseFloat(document.getElementById('inp_sm_tand').value),
        // Top dielectric parameters
        use_top_diel: document.getElementById('chk_top_diel').checked,
        top_diel_h: parseFloat(document.getElementById('inp_top_diel_h').value) * 1e-3,
        top_diel_er: parseFloat(document.getElementById('inp_top_diel_er').value),
        top_diel_tand: parseFloat(document.getElementById('inp_top_diel_tand').value),
        // Ground cutout parameters
        use_gnd_cut: document.getElementById('chk_gnd_cut').checked,
        gnd_cut_w: parseFloat(document.getElementById('inp_gnd_cut_w').value) * 1e-3,
        gnd_cut_h: parseFloat(document.getElementById('inp_gnd_cut_h').value) * 1e-3,
        // Enclosure parameters
        use_enclosure: document.getElementById('chk_enclosure').checked,
        use_side_gnd: document.getElementById('chk_side_gnd').checked,
        use_top_gnd: document.getElementById('chk_top_gnd').checked,
        enclosure_width: nanToNull(parseFloat(document.getElementById('inp_enclosure_width').value) * 1e-3),
        enclosure_height: nanToNull(parseFloat(document.getElementById('inp_enclosure_height').value) * 1e-3),
        max_iters: parseInt(document.getElementById('inp_max_iters').value),
        tolerance: parseFloat(document.getElementById('inp_tolerance').value),
        max_nodes: parseInt(document.getElementById('inp_max_nodes').value),
    };
}

function updateGeometry() {
    const p = getParams();
    currentView = "geometry";

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
                top_gnd_width: p.top_gnd_w,
                via_gap: p.via_gap,
                use_vias: true,
            };
            // Add solder mask options
            if (p.use_sm) {
                options.use_sm = true;
                options.sm_t_sub = p.sm_t_sub;
                options.sm_t_trace = p.sm_t_trace;
                options.sm_t_side = p.sm_t_side;
                options.sm_er = p.sm_er;
                options.sm_tand = p.sm_tand;
            }
            // Add top dielectric options
            if (p.use_top_diel) {
                options.top_diel_h = p.top_diel_h;
                options.top_diel_er = p.top_diel_er;
                options.top_diel_tand = p.top_diel_tand;
            }
            // Add ground cutout options
            if (p.use_gnd_cut) {
                options.gnd_cut_width = p.gnd_cut_w;
                options.gnd_cut_sub_h = p.gnd_cut_h;
            }
            // Add enclosure options
            if (p.use_enclosure) {
                options.use_side_gnd = p.use_side_gnd;
                options.enclosure_width = p.enclosure_width;
                options.enclosure_height = p.enclosure_height;
                if (p.use_top_gnd) {
                    options.boundaries = ["open", "open", "gnd", "gnd"];
                    options.air_top = p.enclosure_height;
                }
            }
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
                top_gnd_width: p.top_gnd_w,
                via_gap: p.via_gap,
                use_vias: true,
            };
            // Add solder mask options
            if (p.use_sm) {
                options.use_sm = true;
                options.sm_t_sub = p.sm_t_sub;
                options.sm_t_trace = p.sm_t_trace;
                options.sm_t_side = p.sm_t_side;
                options.sm_er = p.sm_er;
                options.sm_tand = p.sm_tand;
            }
            // Add top dielectric options
            if (p.use_top_diel) {
                options.top_diel_h = p.top_diel_h;
                options.top_diel_er = p.top_diel_er;
                options.top_diel_tand = p.top_diel_tand;
            }
            // Add ground cutout options
            if (p.use_gnd_cut) {
                options.gnd_cut_width = p.gnd_cut_w;
                options.gnd_cut_sub_h = p.gnd_cut_h;
            }
            // Add enclosure options
            if (p.use_enclosure) {
                options.use_side_gnd = p.use_side_gnd;
                options.enclosure_width = p.enclosure_width;
                options.enclosure_height = p.enclosure_height;
                if (p.use_top_gnd) {
                    options.boundaries = ["open", "open", "gnd", "gnd"];
                    options.air_top = p.enclosure_height;
                }
            }
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
                boundaries: ["open", "open", "open", "gnd"]
            };
            // Add solder mask options
            if (p.use_sm) {
                options.use_sm = true;
                options.sm_t_sub = p.sm_t_sub;
                options.sm_t_trace = p.sm_t_trace;
                options.sm_t_side = p.sm_t_side;
                options.sm_er = p.sm_er;
                options.sm_tand = p.sm_tand;
            }
            // Add top dielectric options
            if (p.use_top_diel) {
                options.top_diel_h = p.top_diel_h;
                options.top_diel_er = p.top_diel_er;
                options.top_diel_tand = p.top_diel_tand;
            }
            // Add ground cutout options
            if (p.use_gnd_cut) {
                options.gnd_cut_width = p.gnd_cut_w;
                options.gnd_cut_sub_h = p.gnd_cut_h;
            }
            // Add enclosure options
            if (p.use_enclosure) {
                options.use_side_gnd = p.use_side_gnd;
                options.enclosure_width = p.enclosure_width;
                options.enclosure_height = p.enclosure_height;
                if (p.use_top_gnd) {
                    options.boundaries = ["open", "open", "gnd", "gnd"];
                    options.air_top = p.enclosure_height;
                }
            }
            solver = new MicrostripSolver(options);
        } else if (p.tl_type === 'stripline') {
            const options = {
                trace_width: p.w,
                substrate_height: p.h,
                trace_thickness: p.t,
                epsilon_r: p.er,
                epsilon_r_top: p.er_top,
                air_top: p.air_top,
                tan_delta: p.tand,
                sigma_cond: p.sigma,
                freq: p.freq,
                nx: p.nx,
                ny: p.ny,
                boundaries: ["open", "open", "gnd", "gnd"]
            };
            // Add solder mask options
            if (p.use_sm) {
                options.use_sm = true;
                options.sm_t_sub = p.sm_t_sub;
                options.sm_t_trace = p.sm_t_trace;
                options.sm_t_side = p.sm_t_side;
                options.sm_er = p.sm_er;
                options.sm_tand = p.sm_tand;
            }
            // Add top dielectric options
            if (p.use_top_diel) {
                options.top_diel_h = p.top_diel_h;
                options.top_diel_er = p.top_diel_er;
                options.top_diel_tand = p.top_diel_tand;
            }
            // Add ground cutout options
            if (p.use_gnd_cut) {
                options.gnd_cut_width = p.gnd_cut_w;
                options.gnd_cut_sub_h = p.gnd_cut_h;
            }
            // Add enclosure options (stripline already has top ground)
            if (p.use_enclosure) {
                options.use_side_gnd = p.use_side_gnd;
                options.enclosure_width = p.enclosure_width;
            }
            solver = new MicrostripSolver(options);
        } else if (p.tl_type === 'diff_stripline') {
            const options = {
                trace_width: p.w,
                substrate_height: p.h,
                trace_thickness: p.t,
                trace_spacing: p.trace_spacing,  // Enable differential mode
                epsilon_r: p.er,
                epsilon_r_top: p.er_top,
                air_top: p.air_top,
                tan_delta: p.tand,
                sigma_cond: p.sigma,
                freq: p.freq,
                nx: p.nx,
                ny: p.ny,
                boundaries: ["open", "open", "gnd", "gnd"]
            };
            // Enclosure options (stripline already has top ground)
            if (p.use_enclosure) {
                options.use_side_gnd = p.use_side_gnd;
                options.enclosure_width = p.enclosure_width;
            }
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
                boundaries: ["open", "open", "open", "gnd"]
            };

            // Solder mask
            if (p.use_sm) {
                options.use_sm = true;
                options.sm_t_sub = p.sm_t_sub;
                options.sm_t_trace = p.sm_t_trace;
                options.sm_t_side = p.sm_t_side;
                options.sm_er = p.sm_er;
                options.sm_tand = p.sm_tand;
            }

            // Top dielectric (embedded microstrip)
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

            // Enclosure options
            if (p.use_enclosure) {
                options.use_side_gnd = p.use_side_gnd;
                options.enclosure_width = p.enclosure_width;
                options.enclosure_height = p.enclosure_height;
                if (p.use_top_gnd) {
                    options.boundaries = ["open", "open", "gnd", "gnd"];
                    options.air_top = p.enclosure_height;
                }
            }

            solver = new MicrostripSolver(options);
        }
    } catch (error) {
        // Log validation errors to the console
        log('ERROR: ' + error.message);
        // Set solver to null to prevent simulation from running with invalid parameters
        solver = null;
    }
}

async function runSimulation() {
    // Check if solver is valid before attempting to run simulation
    if (!solver) {
        log("ERROR: Cannot run simulation - solver initialization failed due to invalid parameters.");
        return;
    }

    const p = getParams();
    const btn = document.getElementById('btn_solve');
    const pbar = document.getElementById('progress_bar');
    const pcont = document.getElementById('progress_container');
    const ptext = document.getElementById('progress_text');

    // Change button to "Stop" mode
    btn.textContent = 'Stop';
    btn.classList.add('stop-mode');
    stopRequested = false;
    pcont.style.display = 'block';
    log("Starting simulation...");

    try {
        // Ensure mesh is generated before solving
        if (!solver.mesh_generated) {
            log("Generating mesh...");
            solver.ensure_mesh();
            log("Mesh generated: " + solver.x.length + "x" + solver.y.length);
        }

        let results;

        log(`Running adaptive analysis (max ${p.max_iters} iterations, max ${p.max_nodes} nodes, tolerance ${p.tolerance})...`);
        ptext.style.display = 'block';

        results = await solver.solve_adaptive({
            max_iters: p.max_iters,
            tolerance: p.tolerance,
            param_tol: 0.001,  // Fixed parameter tolerance
            max_nodes: p.max_nodes,
            onProgress: (info) => {
                const progress = info.iteration / p.max_iters;
                pbar.style.width = (progress * 100) + "%";
                ptext.textContent = `Pass ${info.iteration}/${p.max_iters}: ` +
                                   `Energy err=${info.energy_error.toExponential(2)}, ` +
                                   `Param err=${info.param_error.toExponential(2)}, ` +
                                   `Grid=${info.nodes_x}x${info.nodes_y}`;
                log(`Pass ${info.iteration}: Energy error=${info.energy_error.toExponential(3)}, Param error=${info.param_error.toExponential(3)}, Grid=${info.nodes_x}x${info.nodes_y}`);
            },
            shouldStop: shouldStop
        });

        if (stopRequested) {
            log("Simulation stopped by user");
        }

        // Redraw to show E-field overlay on geometry
        draw();

        // Check if differential results (modes array has 2 elements)
        if (results.modes.length === 2) {
            // Differential microstrip results
            const odd = results.modes.find(m => m.mode === 'odd');
            const even = results.modes.find(m => m.mode === 'even');
            log(`\nDIFFERENTIAL RESULTS:\n` +
                     `======================\n` +
                     `Differential Impedance Z_diff: ${results.Z_diff.toFixed(2)} Ω  (2 × Z_odd)\n` +
                     `Common-Mode Impedance Z_common: ${results.Z_common.toFixed(2)} Ω  (Z_even / 2)\n` +
                     `\nModal Impedances:\n` +
                     `  Odd-Mode  Z_odd:  ${odd.Z0.toFixed(2)} Ω  (εᵣₑff = ${odd.eps_eff.toFixed(3)})\n` +
                     `  Even-Mode Z_even: ${even.Z0.toFixed(2)} Ω  (εᵣₑff = ${even.eps_eff.toFixed(3)})\n` +
                     `\nLosses @ ${(p.freq / 1e9).toFixed(2)} GHz:\n` +
                     `  Odd-Mode:  Diel=${odd.alpha_d.toFixed(4)} dB/m, Cond=${odd.alpha_c.toFixed(4)} dB/m, Total=${odd.alpha_total.toFixed(4)} dB/m\n` +
                     `  Even-Mode: Diel=${even.alpha_d.toFixed(4)} dB/m, Cond=${even.alpha_c.toFixed(4)} dB/m, Total=${even.alpha_total.toFixed(4)} dB/m`);
        } else {
            // Single-ended results
            const mode = results.modes[0];
            log(`\nRESULTS:\n` +
                     `----------------------\n` +
                     `Characteristic Impedance Z0:  ${mode.Z0.toFixed(2)} Ω\n` +
                     `Z0 (complex):  ${mode.Zc.toString()} Ω\n` +
                     `Effective Permittivity: ${mode.eps_eff.toFixed(3)}\n` +
                     `R:             ${mode.RLGC.R.toExponential(3)} Ω/m\n` +
                     `L:             ${mode.RLGC.L.toExponential(3)} H/m\n` +
                     `G:             ${mode.RLGC.G.toExponential(3)} S/m\n` +
                     `C:             ${mode.RLGC.C.toExponential(3)} F/m\n` +
                     `Dielectric Loss: ${mode.alpha_d.toFixed(4)} dB/m\n` +
                     `Conductor Loss:  ${mode.alpha_c.toFixed(4)} dB/m\n` +
                     `Total Loss:      ${mode.alpha_total.toFixed(4)} dB/m`);
        }

    } catch (e) {
        console.error(e);
        log("Error: " + e.message);
    } finally {
        // Restore button to "Solve Physics" mode
        btn.textContent = 'Solve Physics';
        btn.classList.remove('stop-mode');
        pcont.style.display = 'none';
        ptext.style.display = 'none';
        stopRequested = false;
    }
}

let showMesh = false;
let currentView = "geometry";

// Helper function to check if solver is in differential mode
function isDifferentialMode() {
    if (!solver || !solver.Ex || !solver.Ey) return false;
    // In differential mode, Ex and Ey are arrays of 2 arrays (odd and even modes)
    // Check if Ex[0] and Ex[1] are both arrays
    return Array.isArray(solver.Ex) &&
           solver.Ex.length === 2 &&
           Array.isArray(solver.Ex[0]) &&
           Array.isArray(solver.Ex[1]);
}

// Helper function to get Ex/Ey fields (handles differential mode)
function getFields() {
    if (!solver || !solver.Ex || !solver.Ey) {
        return { Ex: null, Ey: null };
    }

    if (isDifferentialMode()) {
        // Differential mode: determine which mode based on currentView
        const modeIndex = currentView === "efield_even" ? 1 : 0;  // Default to odd (0)
        return { Ex: solver.Ex[modeIndex], Ey: solver.Ey[modeIndex] };
    } else {
        // Single-ended mode
        return { Ex: solver.Ex[0], Ey: solver.Ey[0] };
    }
}

// Helper function to get voltage potential (handles differential mode)
function getPotential() {
    if (!solver || !solver.V) {
        return null;
    }

    if (isDifferentialMode()) {
        // Differential mode: determine which mode based on currentView
        const modeIndex = currentView === "potential_even" ? 1 : 0;  // Default to odd (0)
        return solver.V[modeIndex];
    } else {
        // Single-ended mode
        return solver.V[0];
    }
}


// --- Interpolation functions for higher-resolution plots ---

// Find the index i such that arr[i] <= val < arr[i+1]
function find_idx(arr, val) {
    // A linear scan is used here. For very large grids, binary search could be an optimization.
    for (let i = 0; i < arr.length - 1; i++) {
        if (arr[i] <= val && val < arr[i + 1]) {
            return i;
        }
    }
    // Handle the edge case where val is exactly the last element
    if (val === arr[arr.length - 1]) {
        return arr.length - 2;
    }
    return -1;
}

// Bilinear interpolation for a point (x, y) within a grid cell
function bilinearInterpolate(x, y, x1, y1, x2, y2, q11, q12, q21, q22) {
    const denom = (x2 - x1) * (y2 - y1);
    if (denom === 0) {
        // Avoid division by zero if the grid cell has no area
        return q11;
    }
    const w11 = (x2 - x) * (y2 - y);
    const w12 = (x2 - x) * (y - y1);
    const w21 = (x - x1) * (y2 - y);
    const w22 = (x - x1) * (y - y1);
    return (w11 * q11 + w12 * q12 + w21 * q21 + w22 * q22) / denom;
}

// Interpolates data from an old grid (z_old) to a new, finer grid (x_new, y_new)
function interpolateGrid(x_old, y_old, z_old, x_new, y_new) {
    const z_new = Array(y_new.length).fill(0).map(() => Array(x_new.length).fill(0));

    for (let j = 0; j < y_new.length; j++) {
        const y_val = y_new[j];
        const y_idx1 = find_idx(y_old, y_val);
        if (y_idx1 === -1) continue;
        const y_idx2 = y_idx1 + 1;

        for (let i = 0; i < x_new.length; i++) {
            const x_val = x_new[i];
            const x_idx1 = find_idx(x_old, x_val);
            if (x_idx1 === -1) continue;
            const x_idx2 = x_idx1 + 1;

            // Values at the four corners of the cell in the old grid
            const q11 = z_old[y_idx1][x_idx1];
            const q12 = z_old[y_idx2][x_idx1];
            const q21 = z_old[y_idx1][x_idx2];
            const q22 = z_old[y_idx2][x_idx2];

            z_new[j][i] = bilinearInterpolate(
                x_val, y_val,
                x_old[x_idx1], y_old[y_idx1], x_old[x_idx2], y_old[y_idx2],
                q11, q12, q21, q22
            );
        }
    }
    return z_new;
}


function draw(resetZoom = false) {
    if (!solver) return;

    const container = document.getElementById('sim_canvas');

    // Preserve current view state if plot exists (unless resetZoom is requested)
    let currentXRange = null;
    let currentYRange = null;
    if (!resetZoom && container && container.layout && container.layout.xaxis) {
        currentXRange = container.layout.xaxis.range;
        currentYRange = container.layout.yaxis.range;
    }

    let zData = [];
    let title = "";
    let colorscale = "Viridis";
    let zTitle = "";
    let shapes = [];
    let xMM, yMM, nx, ny, nyDisplay;

    // --------------------
    // DATA SELECTION
    // --------------------
    if (currentView === "geometry") {
        title = "Transmission Line Geometry";

        // Determine display bounds
        const maxY = Math.min(solver.h + 5 * solver.w, solver.dielectrics.reduce((max, d) => Math.max(max, d.y_max), 0));

        // Draw dielectrics as rectangles (color by epsilon_r)
        for (const diel of solver.dielectrics) {
            if (diel.y_min > maxY) continue;

            const yMax = Math.min(diel.y_max, maxY);
            const er = diel.epsilon_r;

            // Color mapping: air (1.0) = white, higher er = green shades
            let fillcolor;
            if (er <= 1.01) {
                fillcolor = 'rgba(255, 255, 255, 0.8)';
            } else {
                // Green shades for dielectrics
                const intensity = Math.min(255, 100 + (er - 1) * 30);
                fillcolor = `rgba(100, ${intensity}, 100, 0.8)`;
            }

            shapes.push({
                type: 'rect',
                x0: diel.x_min * 1000,
                y0: diel.y_min * 1000,
                x1: diel.x_max * 1000,
                y1: yMax * 1000,
                fillcolor: fillcolor,
                line: { color: 'rgba(128, 128, 128, 0.3)', width: 0.5 },
                layer: 'below'
            });
        }

        // Draw conductors as rectangles (signal = orange, ground = dark gray)
        for (const cond of solver.conductors) {
            if (cond.y_min > maxY) continue;

            const yMax = Math.min(cond.y_max, maxY);
            const fillcolor = cond.is_signal ?
                'rgba(217, 119, 6, 1.0)' :  // Orange for signal
                'rgba(51, 51, 51, 1.0)';     // Dark gray for ground

            shapes.push({
                type: 'rect',
                x0: cond.x_min * 1000,
                y0: cond.y_min * 1000,
                x1: cond.x_max * 1000,
                y1: yMax * 1000,
                fillcolor: fillcolor,
                line: { color: 'rgba(0, 0, 0, 0.5)', width: 1 },
                layer: 'above'
            });
        }

        // If solution available, overlay E-field contours
        if (solver.solution_valid && solver.mesh_generated) {
            nx = solver.x.length;
            ny = solver.y.length;

            // Limit display Y
            const yArr = Array.from(solver.y);
            const maxYIdx = yArr.findIndex(y => y > maxY);
            nyDisplay = maxYIdx > 0 ? maxYIdx : ny;

            xMM = Array.from(solver.x, v => v * 1000);
            yMM = yArr.slice(0, nyDisplay).map(v => v * 1000);

            // Compute E-field magnitude
            const { Ex, Ey } = getFields();
            if (Ex && Ey && Ex.length >= nyDisplay) {
                for (let i = 0; i < nyDisplay; i++) {
                    const row = [];
                    if (Ex[i] && Ey[i]) {
                        for (let j = 0; j < nx; j++) {
                            row.push(Math.hypot(Ex[i][j], Ey[i][j]));
                        }
                    }
                    zData.push(row);
                }
            }

            // In geometry view with differential solver, always show odd mode
            const modeLabel = isDifferentialMode() ? " (Odd Mode)" : "";
            title = `Transmission Line Geometry with E-field${modeLabel}`;
        } else {
            // No solution - just axis scaling
            xMM = [0, solver.w * 2000];
            yMM = [0, maxY * 1000];
        }
    }

    else if ((currentView === "potential" || currentView === "potential_odd" || currentView === "potential_even") && solver.solution_valid) {
        // Ensure mesh exists for field visualization
        if (!solver.mesh_generated) {
            solver.ensure_mesh();
        }

        nx = solver.x.length;
        ny = solver.y.length;

        // Limit display Y
        const yArr = Array.from(solver.y);
        const maxY = Math.min(yArr[ny - 1], solver.h + 5 * solver.w);
        const maxYIdx = yArr.findIndex(y => y > maxY);
        nyDisplay = maxYIdx > 0 ? maxYIdx : ny;

        xMM = Array.from(solver.x, v => v * 1000);
        yMM = yArr.slice(0, nyDisplay).map(v => v * 1000);

        let modeLabel = "";
        if (currentView === "potential_odd") {
            modeLabel = " (Odd Mode)";
        } else if (currentView === "potential_even") {
            modeLabel = " (Even Mode)";
        }
        title = `Electric Potential${modeLabel} (V)`;
        zTitle = "Volts";

        const V = getPotential();
        if (V && V.length >= nyDisplay) {
            for (let i = 0; i < nyDisplay; i++) {
                zData.push(V[i].slice(0, nx));
            }
        }
    }

    else if ((currentView === "efield" || currentView === "efield_odd" || currentView === "efield_even") && solver.solution_valid) {
        // Ensure mesh exists for field visualization
        if (!solver.mesh_generated) {
            solver.ensure_mesh();
        }

        nx = solver.x.length;
        ny = solver.y.length;

        // Limit display Y
        const yArr = Array.from(solver.y);
        const maxY = Math.min(yArr[ny - 1], solver.h + 5 * solver.w);
        const maxYIdx = yArr.findIndex(y => y > maxY);
        nyDisplay = maxYIdx > 0 ? maxYIdx : ny;

        xMM = Array.from(solver.x, v => v * 1000);
        yMM = yArr.slice(0, nyDisplay).map(v => v * 1000);

        let modeLabel = "";
        if (currentView === "efield_odd") {
            modeLabel = " (Odd Mode)";
        } else if (currentView === "efield_even") {
            modeLabel = " (Even Mode)";
        }
        title = `|E| Field Magnitude${modeLabel} (V/m)`;
        zTitle = "V/m";

        const { Ex, Ey } = getFields();
        if (Ex && Ey && Ex.length >= nyDisplay) {
            for (let i = 0; i < nyDisplay; i++) {
                const row = [];
                if (Ex[i] && Ey[i]) {
                    for (let j = 0; j < nx; j++) {
                        row.push(Math.hypot(Ex[i][j], Ey[i][j]));
                    }
                }
                zData.push(row);
            }
        }
    }

    else {
        title = "No Data Available";
        // Create minimal dummy data
        xMM = [0, (solver.w || 1) * 2000];
        yMM = [0, (solver.h || 1) * 1000];
    }

    // --------------------
    // INTERPOLATION FOR SMOOTHER PLOTS
    // --------------------
    const INTERP_ENABLED = true; // Control flag for interpolation
    const INTERP_FACTOR = 5;     // Interpolation multiplier (e.g., 5x resolution)

    if (INTERP_ENABLED && zData.length > 0 && (currentView.includes("potential") || currentView.includes("efield"))) {
        
        const x_old = Array.from(solver.x); // Original grid X (meters)
        const y_old = Array.from(solver.y).slice(0, nyDisplay); // Original grid Y (meters)

        if (x_old.length > 1 && y_old.length > 1 && zData.length > 1 && zData[0].length > 1) {
            const nx_interp = (nx - 1) * INTERP_FACTOR + 1;
            const ny_interp = (nyDisplay - 1) * INTERP_FACTOR + 1;

            // Create a new, uniformly spaced, finer grid for interpolation
            const x_new = new Float64Array(nx_interp);
            const y_new = new Float64Array(ny_interp);

            const x_min = x_old[0];
            const x_max = x_old[x_old.length - 1];
            const y_min = y_old[0];
            const y_max = y_old[y_old.length - 1];

            for (let i = 0; i < nx_interp; i++) {
                x_new[i] = x_min + (x_max - x_min) * i / (nx_interp - 1);
            }
            for (let j = 0; j < ny_interp; j++) {
                y_new[j] = y_min + (y_max - y_min) * j / (ny_interp - 1);
            }

            // Perform interpolation from the original zData to the new grid
            const z_interp = interpolateGrid(x_old, y_old, zData, x_new, y_new);

            // Update plot variables with the new, higher-resolution data
            xMM = Array.from(x_new, v => v * 1000);
            yMM = Array.from(y_new, v => v * 1000);
            zData = z_interp;
        }
    }


    // --------------------
    // MAIN FIELD TRACE
    // --------------------
    let traces = [];

    if (currentView === "geometry" && zData.length > 0) {
        const { Ex, Ey } = getFields();

        //const fieldTraces = [];

        //const stepX = Math.max(1, Math.floor(nx / 30));        // density control
        //const stepY = Math.max(1, Math.floor(nyDisplay / 30));
        //const scale = 0.8; // arrow length in mm


        //if (Ex && Ey && Ex.length >= nyDisplay) {
        //    const xLines = [];
        //    const yLines = [];

        //    for (let i = 0; i < nyDisplay; i += stepY) {
        //        for (let j = 0; j < nx; j += stepX) {
        //            const ex = Ex[i]?.[j];
        //            const ey = Ey[i]?.[j];
        //            if (!ex || !ey) continue;

        //            const mag = Math.hypot(ex, ey);
        //            if (mag === 0) continue;

        //            // Base point (mm)
        //            const x0 = xMM[j];
        //            const y0 = yMM[i];

        //            // Direction (normalized)
        //            const dx = (ex / mag) * scale;
        //            const dy = (ey / mag) * scale;

        //            // Line segment
        //            xLines.push(x0, x0 + dx, null);
        //            yLines.push(y0, y0 + dy, null);
        //        }
        //    }

        //    traces.push({
        //        type: "scatter",
        //        mode: "lines",
        //        x: xLines,
        //        y: yLines,
        //        line: {
        //            width: 1.2,
        //            color: "black"
        //        },
        //        hoverinfo: "skip",
        //        name: "E-field lines"
        //    });
        //}
        //
        //traces.push({
        //    type: "contour",
        //    x: xMM,
        //    y: yMM,
        //    z: zData,
        //    colorscale: "Hot",
        //    opacity: 0.6,
        //    contours: {
        //        showlines: true,
        //        coloring: "heatmap",
        //        ncontours: 15
        //    },
        //    colorbar: {
        //        title: "|E| (V/m)",
        //        len: 0.6
        //    },
        //    hovertemplate:
        //        "x: %{x:.2f} mm<br>" +
        //        "y: %{y:.2f} mm<br>" +
        //        "|E|: %{z:.3e} V/m<extra></extra>"
        //});

        traces.push(
            makeStreamlineTraceFromConductors(
                Ex,
                Ey,
                solver.x,
                solver.y,
                solver.conductors,
                50
            )
        );

    } else if (currentView === "geometry") {
        // Geometry only - invisible scatter for axis scaling
        traces.push({
            type: "scatter",
            x: xMM,
            y: yMM,
            mode: "markers",
            marker: { size: 0, opacity: 0 },
            showlegend: false,
            hoverinfo: "skip"
        });
    } else if (zData.length > 0) {
        // Field views only - use heatmap
        traces.push({
            type: "contour",
            x: xMM,
            y: yMM,
            z: zData,
            colorscale: colorscale,
            contours: {
                coloring: 'heatmap',
                showlines: true,
            },
            line: {
              smoothing: 1.3,
              width: 0.5
            },
            colorbar: {
                title: zTitle,
                len: 0.8
            },
            hovertemplate:
                "x: %{x:.2f} mm<br>" +
                "y: %{y:.2f} mm<br>" +
                "value: %{z:.3e}<extra></extra>"
        });
    }

    // MESH OVERLAY
    if (showMesh && solver.solution_valid) {
        const stepX = 1;
        const stepY = 1;

        for (let j = 0; j < nx; j += stepX) {
            traces.push({
                type: "scatter",
                x: [xMM[j], xMM[j]],
                y: [yMM[0], yMM[nyDisplay - 1]],
                mode: "lines",
                line: { width: 0.2, color: "white" },
                showlegend: false,
                hoverinfo: "skip"
            });
        }

        for (let i = 0; i < nyDisplay; i += stepY) {
            traces.push({
                type: "scatter",
                x: [xMM[0], xMM[nx - 1]],
                y: [yMM[i], yMM[i]],
                mode: "lines",
                line: { width: 0.2, color: "white" },
                showlegend: false,
                hoverinfo: "skip"
            });
        }
    }

    // --------------------
    // UI MENUS
    // --------------------
    const layout = {
        title: title,
        xaxis: {
            title: "Width (mm)",
            scaleanchor: "y",
            scaleratio: 1,
            range: currentXRange  // Preserve zoom/pan
        },
        yaxis: {
            title: "Height (mm)",
            range: currentYRange  // Preserve zoom/pan
        },
        margin: { l: 70, r: 90, t: 50, b: 60 },
        hovermode: "closest",
        dragmode: "pan",
        plot_bgcolor: "#f8f9fa",
        shapes: shapes,  // Add vector shapes for geometry

        updatemenus: [
            {
                x: 0.01,
                y: 1.15,
                showactive: true,
                active: (() => {
                    if (currentView === "geometry") return 0;
                    if (isDifferentialMode()) {
                        // Differential: Geometry, Potential(odd), Potential(even), E-field(odd), E-field(even)
                        if (currentView === "potential_odd" || currentView === "potential") return 1;
                        if (currentView === "potential_even") return 2;
                        if (currentView === "efield_odd" || currentView === "efield") return 3;
                        if (currentView === "efield_even") return 4;
                    } else {
                        // Single-ended: Geometry, Potential, E-field
                        if (currentView === "potential") return 1;
                        if (currentView === "efield") return 2;
                    }
                    return 0;
                })(),
                buttons: (() => {
                    const buttons = [
                        {
                            label: "Geometry",
                            method: "skip",
                            args: []
                        }
                    ];

                    if (solver.solution_valid) {
                        if (isDifferentialMode()) {
                            // Add separate buttons for odd and even modes
                            buttons.push({
                                label: "Potential (odd)",
                                method: "skip",
                                args: []
                            });
                            buttons.push({
                                label: "Potential (even)",
                                method: "skip",
                                args: []
                            });
                            buttons.push({
                                label: "E-field (odd)",
                                method: "skip",
                                args: []
                            });
                            buttons.push({
                                label: "E-field (even)",
                                method: "skip",
                                args: []
                            });
                        } else {
                            // Single-ended: just one Potential and one E-field button
                            buttons.push({
                                label: "Potential",
                                method: "skip",
                                args: []
                            });
                            buttons.push({
                                label: "|E| Field",
                                method: "skip",
                                args: []
                            });
                        }
                    }

                    return buttons;
                })()
            }
        ]
    };

    const config = {
        responsive: true,
        displayModeBar: true,
        scrollZoom: true,
        modeBarButtonsToAdd: [
            {
                name: "Toggle Mesh",
                icon: Plotly.Icons.grid,
                click: () => {
                    showMesh = !showMesh;
                    draw();
                }
            },
            {
                name: "Auto Z Scale",
                icon: Plotly.Icons.autoscale,
                click: () => {
                    Plotly.relayout(container, { "zaxis.autorange": true });
                }
            }
        ]
    };

    Plotly.react(container, traces, layout, config);

    if (!container._viewListenerBound) {
        container.on('plotly_buttonclicked', (event) => {
            if (event.menu.active === 0) {
                currentView = "geometry";
            } else if (isDifferentialMode()) {
                // Differential mode: Geometry(0), Potential(odd)(1), Potential(even)(2), E-field(odd)(3), E-field(even)(4)
                if (event.menu.active === 1) {
                    currentView = "potential_odd";
                } else if (event.menu.active === 2) {
                    currentView = "potential_even";
                } else if (event.menu.active === 3) {
                    currentView = "efield_odd";
                } else if (event.menu.active === 4) {
                    currentView = "efield_even";
                }
            } else {
                // Single-ended mode: Geometry(0), Potential(1), E-field(2)
                if (event.menu.active === 1) {
                    currentView = "potential";
                } else if (event.menu.active === 2) {
                    currentView = "efield";
                }
            }
            draw();
        });
        container._viewListenerBound = true;
    }

}

function resizeCanvas() {
    const container = document.getElementById('sim_canvas');
    if (container) {
        Plotly.Plots.resize(container);
    }
}

function viridis(t) {
    // Simple heatmap approximation
    t = Math.max(0, Math.min(1, t));
    // R, G, B interpolation
    const r = Math.floor(255 * Math.sin(t * 2));
    const g = Math.floor(255 * Math.sin(t * 3));
    const b = Math.floor(255 * Math.cos(t * 1.5));
    // Better pseudocolor: (Blue -> Cyan -> Green -> Yellow -> Red)
    // Manual standard mapping for clarity:
    if(t < 0.25) return `rgb(0, ${Math.floor(t*4*255)}, 255)`;
    if(t < 0.5) return `rgb(0, 255, ${Math.floor((0.5-t)*4*255)})`;
    if(t < 0.75) return `rgb(${Math.floor((t-0.5)*4*255)}, 255, 0)`;
    return `rgb(255, ${Math.floor((1-t)*4*255)}, 0)`;
}

function bindEvents() {
    document.getElementById('btn_solve').onclick = () => {
        const btn = document.getElementById('btn_solve');
        if (btn.textContent === 'Stop') {
            // Stop the simulation
            stopRequested = true;
            log("Stop requested...");
        } else {
            // Start the simulation
            updateGeometry(); // Ensure geometry is updated with latest parameters
            runSimulation();
        }
    };

    // Real-time geometry updates for all parameter inputs
    const geometryInputs = [
        'inp_w', 'inp_h', 'inp_t', 'inp_er', 'inp_tand', 'inp_sigma', 'inp_freq',
        'inp_trace_spacing',
        'inp_gap', 'inp_top_gnd_w', 'inp_via_gap',
        'inp_air_top', 'inp_er_top',
        'inp_sm_t_sub', 'inp_sm_t_trace', 'inp_sm_t_side', 'inp_sm_er', 'inp_sm_tand',
        'inp_top_diel_h', 'inp_top_diel_er', 'inp_top_diel_tand',
        'inp_gnd_cut_w', 'inp_gnd_cut_h',
        'inp_enclosure_width', 'inp_enclosure_height'
    ];

    geometryInputs.forEach(id => {
        const el = document.getElementById(id);
        if (el) {
            el.addEventListener('input', () => {
                updateGeometry();
                draw();
            });
        }
    });

    // Real-time updates for checkboxes
    const geometryCheckboxes = [
        'chk_solder_mask', 'chk_top_diel', 'chk_gnd_cut', 'chk_enclosure', 'chk_side_gnd', 'chk_top_gnd'
    ];

    geometryCheckboxes.forEach(id => {
        const el = document.getElementById(id);
        if (el) {
            el.addEventListener('change', () => {
                updateGeometry();
                draw();
            });
        }
    });

    // Transmission line type selector - reset zoom when type changes
    document.getElementById('tl_type').addEventListener('change', () => {
        updateGeometry();
        draw(true);  // Reset zoom/pan for new geometry
    });
}


function init() {
    bindEvents();
    updateGeometry();
    draw();
    resizeCanvas();
    window.addEventListener('resize', resizeCanvas);
}

// Start when DOM is ready
window.addEventListener('DOMContentLoaded', init);
