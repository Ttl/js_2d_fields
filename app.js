import { MicrostripSolver } from './microstrip.js';
import { CONSTANTS } from './field_solver.js';
const Plotly = window.Plotly;

let solver = null;

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
        log("Running full analysis...");
        const results = await solver.perform_analysis((progress) => {
            pbar.style.width = (progress * 100) + "%";
        });

        currentView = "efield";
        draw();

        log(`\nRESULTS:\n` +
                 `----------------------\n` +
                 `Characteristic Impedance Z0:  ${results.Z0.toFixed(2)} Ω\n` +
                 `Z0 (complex):  ${results.Zc.toString()} Ω\n` +
                 `Effective Permittivity εᵣₑff: ${results.eps_eff.toFixed(3)}\n` +
                 `R:             ${results.RLGC.R.toExponential(3)} Ω/m\n` +
                 `L:             ${results.RLGC.L.toExponential(3)} H/m\n` +
                 `G:             ${results.RLGC.G.toExponential(3)} S/m\n` +
                 `C:             ${results.RLGC.C.toExponential(3)} F/m\n` +
                 `Dielectric Loss: ${results.alpha_diel_db_m.toFixed(4)} dB/m\n` +
                 `Conductor Loss:  ${results.alpha_cond_db_m.toFixed(4)} dB/m\n` +
                 `Total Loss:      ${results.total_alpha_db_m.toFixed(4)} dB/m`);

    } catch (e) {
        console.error(e);
        log("Error: " + e.message);
    } finally {
        btn.disabled = false;
        pcont.style.display = 'none';
    }
}

let contourCount = 20;
let showMesh = false;
let currentView = "geometry";

function draw() {
    if (!solver) return;

    const container = document.getElementById('sim_canvas');

    const nx = solver.x.length;
    const ny = solver.y.length;

    // Limit display Y
    const maxY = Math.min(solver.y[ny - 1], solver.h + 5 * solver.w);
    const maxYIdx = solver.y.findIndex(y => y > maxY);
    const nyDisplay = maxYIdx > 0 ? maxYIdx : ny;

    // Coordinates in mm
    const xMM = solver.x.map(v => v * 1000);
    const yMM = solver.y.slice(0, nyDisplay).map(v => v * 1000);

    let zData = [];
    let title = "";
    let colorscale = "Viridis";
    let zTitle = "";

    // --------------------
    // DATA SELECTION
    // --------------------
    if (currentView === "geometry") {
        title = "Microstrip Geometry";
        zTitle = "εᵣ / Conductor";

        colorscale = [
            [0.0, "#ffffff"],
            [0.3, "#64c864"],
            [0.7, "#333333"],
            [1.0, "#d97706"]
        ];

        for (let i = 0; i < nyDisplay; i++) {
            const row = [];
            for (let j = 0; j < nx; j++) {
                if (solver.conductor_mask[i][j]) {
                    row.push(solver.ground_mask[i][j] ? -2 : -1);
                } else {
                    row.push(solver.epsilon_r[i][j]);
                }
            }
            zData.push(row);
        }
    }

    else if (currentView === "potential" && solver.solution_valid) {
        title = "Electric Potential (V)";
        zTitle = "Volts";

        for (let i = 0; i < nyDisplay; i++) {
            zData.push(solver.V[i].slice(0, nx));
        }
    }

    else if (currentView === "efield" && solver.solution_valid) {
        title = "|E| Field Magnitude (V/m)";
        zTitle = "V/m";

        for (let i = 0; i < nyDisplay; i++) {
            const row = [];
            for (let j = 0; j < nx; j++) {
                row.push(Math.hypot(solver.Ex[i][j], solver.Ey[i][j]));
            }
            zData.push(row);
        }
    }

    else {
        title = "No Data Available";
        zData = Array(nyDisplay).fill(0).map(() => Array(nx).fill(0));
    }

    // --------------------
    // MAIN FIELD TRACE
    // --------------------
    const traces = [{
        type: currentView === "geometry" ? "contour" : "heatmap",
        x: xMM,
        y: yMM,
        z: zData,
        colorscale: colorscale,
        contours: {
            showlines: true,
            coloring: "heatmap",
            ncontours: contourCount
        },
        colorbar: {
            title: zTitle,
            len: 0.8
        },
        hovertemplate:
            "x: %{x:.2f} mm<br>" +
            "y: %{y:.2f} mm<br>" +
            "value: %{z:.3e}<extra></extra>"
    }];

    // --------------------
    // MESH OVERLAY
    // --------------------
    if (showMesh) {
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
            scaleratio: 1
        },
        yaxis: {
            title: "Height (mm)"
        },
        margin: { l: 70, r: 90, t: 50, b: 60 },
        hovermode: "closest",
        dragmode: "pan",
        plot_bgcolor: "#f8f9fa",

        updatemenus: [
            {
                x: 0.01,
                y: 1.15,
                showactive: true,
                active: currentView === "geometry" ? 0 : currentView === "potential" ? 1 : 2,
                buttons: [
                    {
                        label: "Geometry",
                        method: "skip",
                        args: []
                    },
                    {
                        label: "Potential",
                        method: "skip",
                        args: [],
                        visible: solver.solution_valid
                    },
                    {
                        label: "|E| Field",
                        method: "skip",
                        args: [],
                        visible: solver.solution_valid
                    }
                ]
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
                name: "More Contours",
                icon: Plotly.Icons.zoomIn,
                click: () => {
                    contourCount += 5;
                    draw();
                }
            },
            {
                name: "Fewer Contours",
                icon: Plotly.Icons.zoomOut,
                click: () => {
                    contourCount = Math.max(5, contourCount - 5);
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
            if (event.menu.active === 0) currentView = "geometry";
            else if (event.menu.active === 1) currentView = "potential";
            else if (event.menu.active === 2) currentView = "efield";
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
    document.getElementById('btn_update_geo').onclick = () => {
        updateGeometry();
        draw();
    };
    document.getElementById('btn_solve').onclick = () => runSimulation();

    [
      'inp_contours',
      'inp_vmin',
      'inp_vmax',
      'chk_auto_range',
      'chk_mesh',
      'sel_view_mode'
    ].forEach(id => {
        const el = document.getElementById(id);
        if (!el) return;
        el.onchange = () => draw();
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
