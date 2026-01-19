import { MicrostripSolver } from './microstrip.js';
import { CONSTANTS } from './field_solver.js';

let solver = null;
const canvas = document.getElementById('sim_canvas');
const ctx = document.getElementById('sim_canvas').getContext('2d');

function log(msg) {
    const c = document.getElementById('console_out');
    c.textContent += msg + "\n";
    c.scrollTop = c.scrollHeight;
}

function getParams() {
    return {
        w: parseFloat(document.getElementById('inp_w').value) * 1e-3,
        h: parseFloat(document.getElementById('inp_h').value) * 1e-3,
        t: parseFloat(document.getElementById('inp_t').value) * 1e-6,
        er: parseFloat(document.getElementById('inp_er').value),
        tand: parseFloat(document.getElementById('inp_tand').value),
        sigma: parseFloat(document.getElementById('inp_sigma').value),
        freq: parseFloat(document.getElementById('inp_freq').value) * 1e9,
        nx: parseInt(document.getElementById('inp_nx').value),
        ny: parseInt(document.getElementById('inp_ny').value),
    };
}

function updateGeometry() {
    const p = getParams();
    solver = new MicrostripSolver(p.w, p.h, p.t, p.er, p.tand, p.sigma, p.freq, p.nx, p.ny);
    log("Geometry updated. Grid: " + solver.x.length + "x" + solver.y.length);
}

async function runSimulation() {
    const btn = document.getElementById('btn_solve');
    const pbar = document.getElementById('progress_bar');
    const pcont = document.getElementById('progress_container');

    btn.disabled = true;
    pcont.style.display = 'block';
    log("Starting simulation...");

    try {
        // 1. Solve with Dielectric
        log("Step 1/2: Solving field with dielectric...");
        await solver.solve_laplace_iterative(false, (i, max, diff) => {
            const pct = (i / max) * 50;
            pbar.style.width = pct + "%";
        });
        const C = solver.calculate_capacitance(false);

        // Compute fields for loss calculation and visualization
        solver.compute_fields();
        draw(); // Update view with fields

        // 2. Solve with Vacuum
        log("Step 2/2: Solving field with vacuum...");

        // Backup Voltage from Step 1 to restore for visualization later
        const V_diel = solver.V.map(row => new Float64Array(row));
        const Ex_diel = solver.Ex;
        const Ey_diel = solver.Ey;

        // Reset V for vacuum solve
        solver.V = solver.V.map((row, i) => 
            row.map((val, j) => solver.conductor_mask[i][j] ? val : 0)
        );

        await solver.solve_laplace_iterative(true, (i, max, diff) => {
            const pct = 50 + (i / max) * 50;
            pbar.style.width = pct + "%";
        });
        const C0 = solver.calculate_capacitance(true);

        // 3. Post Process
        const eps_eff = C / C0;
        const Z0_real = 1 / (CONSTANTS.C * Math.sqrt(C * C0)); // Renamed to avoid conflict with complex Z0

        // Restore dielectric fields for loss calc
        solver.V = V_diel;
        solver.Ex = Ex_diel;
        solver.Ey = Ey_diel;

        // Calculate losses using the new FieldSolver2D methods
        const alpha_cond = solver.calculate_conductor_loss(solver.Ex, solver.Ey, Z0_real);
        const alpha_diel = solver.calculate_dielectric_loss(solver.Ex, solver.Ey, Z0_real);
        const alpha_total = alpha_cond + alpha_diel;

        // Calculate RLGC and complex Z0
        const { Zc, rlgc, eps_eff_mode } = solver.rlgc(alpha_cond, alpha_diel, C, Z0_real);

        log(`\nRESULTS:\n` +
                 `----------------------\n` +
                 `Capacitance:   ${(C*1e12).toFixed(2)} pF/m\n` +
                 `Z0 (real):     ${Z0_real.toFixed(2)} Ω\n` +
                 `Z0 (complex):  ${Zc.toString()} Ω\n` +
                 `Eps_eff:       ${eps_eff.toFixed(3)}\n` +
                 `Eps_eff (mode):${eps_eff_mode.toFixed(3)}\n` +
                 `R:             ${rlgc.R.toExponential(3)} Ω/m\n` +
                 `L:             ${rlgc.L.toExponential(3)} H/m\n` +
                 `G:             ${rlgc.G.toExponential(3)} S/m\n` +
                 `C:             ${rlgc.C.toExponential(3)} F/m\n` +
                 `Dielectric Loss: ${alpha_diel.toFixed(4)} dB/m\n` +
                 `Conductor Loss:  ${alpha_cond.toFixed(4)} dB/m\n` +
                 `Total Loss:      ${alpha_total.toFixed(4)} dB/m`);

    } catch (e) {
        console.error(e);
        log("Error: " + e.message);
    } finally {
        btn.disabled = false;
        pcont.style.display = 'none';
    }
}

function resizeCanvas() {
    const container = canvas.parentElement;
    canvas.width = container.clientWidth;
    canvas.height = container.clientHeight;
    if(solver) draw();
}

function draw() {
    if (!solver) return;
    const w = canvas.width;
    const h = canvas.height;
    ctx.clearRect(0, 0, w, h);

    const viewMode = document.getElementById('sel_view_mode').value;
    const showMesh = document.getElementById('chk_mesh').checked;

    // Scaling
    // Find bounding box to fit in canvas with margin
    const domainW = solver.domain_width;
    // Limit display Y to relevant area (up to 5x trace height)
    const maxY = Math.min(solver.y.slice(-1)[0], solver.h + 5 * solver.w);
    const domainH = maxY;

    const scaleX = (w - 40) / domainW;
    const scaleY = (h - 40) / domainH;
    const scale = Math.min(scaleX, scaleY);

    const offsetX = (w - domainW * scale) / 2;
    const offsetY = h - (h - domainH * scale) / 2; // Bottom origin

    const toScreenX = (x) => offsetX + x * scale;
    const toScreenY = (y) => offsetY - y * scale;

    // Helper to draw Rects
    const drawCell = (i, j, color) => {
        const x0 = toScreenX(solver.x[j]);
        const y0 = toScreenY(solver.y[i]);
        // dx/dy are widths
        const x1 = toScreenX(solver.x[j+1]);
        const y1 = toScreenY(solver.y[i+1]);

        ctx.fillStyle = color;
        ctx.fillRect(x0, y1, x1-x0, y0-y1); // Canvas Y inverted
    };

    // Draw Field / Geometry
    const ny = solver.y.length;
    const nx = solver.x.length;

    // Determine min/max for colormapping
    let maxVal = 0;
    if (viewMode === 'efield' && solver.solution_valid) {
         for(let i=0; i<ny; i++) 
            for(let j=0; j<nx; j++) 
                 maxVal = Math.max(maxVal, Math.sqrt(solver.Ex[i][j]**2 + solver.Ey[i][j]**2));
    }

    for (let i = 0; i < ny - 1; i++) {
        if (solver.y[i] > maxY) break;

        for (let j = 0; j < nx - 1; j++) {
            let color;

            if (viewMode === 'geometry') {
                if (solver.conductor_mask[i][j]) {
                    color = solver.ground_mask[i][j] ? '#333' : '#d97706'; // Dark gnd, Orange trace
                } else {
                    // Substrate vs Air
                    const er = solver.epsilon_r[i][j];
                    color = er > 1.1 ? `rgba(100, 200, 100, ${0.2 + er/10})` : '#fff';
                }
            } else if (viewMode === 'potential' && solver.solution_valid) {
                if (solver.conductor_mask[i][j]) {
                     color = solver.ground_mask[i][j] ? '#000' : '#d97706';
                } else {
                    const v = solver.V[i][j];
                    color = viridis(v);
                }
            } else if (viewMode === 'efield' && solver.solution_valid) {
                const emag = Math.sqrt(solver.Ex[i][j]**2 + solver.Ey[i][j]**2);
                color = viridis(Math.pow(emag / (maxVal || 1), 0.5)); // Gamma correct
            } else {
                color = '#eee';
            }

            drawCell(i, j, color);
        }
    }

    // Draw Mesh Lines
    if (showMesh) {
        ctx.strokeStyle = 'rgba(0,0,0,0.1)';
        ctx.lineWidth = 1;
        ctx.beginPath();

        // V lines
        for(let j=0; j<nx; j++) {
            const x = toScreenX(solver.x[j]);
            ctx.moveTo(x, toScreenY(0));
            ctx.lineTo(x, toScreenY(maxY));
        }
        // H lines
        for(let i=0; i<ny; i++) {
            if (solver.y[i] > maxY) break;
            const y = toScreenY(solver.y[i]);
            ctx.moveTo(toScreenX(0), y);
            ctx.lineTo(toScreenX(domainW), y);
        }
        ctx.stroke();
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
    document.getElementById('btn_update_geo').onclick = () => {
        updateGeometry();
        draw();
    };
    document.getElementById('btn_solve').onclick = () => runSimulation();

    ['chk_mesh', 'sel_view_mode'].forEach(id => {
        document.getElementById(id).onchange = () => draw();
    });
}


// Initial setup and event binding
function init() {
    bindEvents();
    updateGeometry();
    resizeCanvas();
    window.addEventListener('resize', () => resizeCanvas());
}

// Start when DOM is ready
window.addEventListener('DOMContentLoaded', init);
