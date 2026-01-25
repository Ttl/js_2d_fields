import { makeStreamlineTraceFromConductors } from './streamlines.js';
import { computeSParamsSingleEnded, computeSParamsDifferential, sParamTodB } from './sparameters.js';

const Plotly = window.Plotly;

let showMesh = false;
let currentView = "geometry";
let zMin = null;
let zMax = null;

// Globals imported from app.js
let getSolver = () => null;
let getFrequencySweepResults = () => null;
let getInputValue = () => NaN;

// Function to set globals from app.js
function setGlobals(globals) {
    getSolver = globals.getSolver || (() => null);
    getFrequencySweepResults = globals.getFrequencySweepResults || (() => null);
    getInputValue = globals.getInputValue || (() => NaN);
}

// Helper to access globals
const get = {
    solver: () => getSolver(),
    frequencySweepResults: () => getFrequencySweepResults(),
    inputValue: (id) => getInputValue(id)
};

function contourScaledB(min, max, n) {
    let eMin = Math.max(Math.max(1, max*1e-2), min);
    eMin = Math.log10(Math.max(eMin, 0.1));
    let eMax = Math.log10(Math.max(eMin + 0.1, Math.max(max, 0.1)));
    eMax = Math.max(eMin + 0.1, eMax);
    const logStep = n == 0 ? 1 : Math.abs((eMax - eMin)) / n;
    return [eMin, eMax, logStep];
}

// Export functions to get/set scale range for current view
function getScaleRange() {
    return { min: zMin, max: zMax, view: currentView };
}

function setScaleRange(min, max) {
    zMin = min;
    zMax = max;

    const container = document.getElementById('sim_canvas');
    if (!container || !container.data) return;

    // For geometry view with contours, update contour properties
    if (currentView === "geometry") {
        // Find the contour trace (if it exists)
        const contourTraceIdx = container.data.findIndex(trace =>
            trace.type === 'contour' && trace.name === 'E-field contours'
        );

        if (contourTraceIdx !== -1) {
            // Get number of contours from plotOptions
            const plotOptions = getPlotOptions();
            const n = plotOptions.contours;
            const limits = contourScaledB(min, max, n);

            if (n > 0) {
                Plotly.restyle(container, {
                    'contours.start': limits[0],
                    'contours.end': limits[1],
                    'contours.size': limits[2]
                }, [contourTraceIdx]);
            }
        }
    } else {
        // For E-field and potential views, update heatmap/contour zmin/zmax
        Plotly.restyle(container, {
            zmin: min,
            zmax: max
        });
    }
}

function draw(resetZoom = false) {
    const solver = get.solver();
    if (!solver) return;

    const container = document.getElementById('sim_canvas');
    const plotOptions = getPlotOptions();

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

    // View selection
    if (currentView === "geometry") {
        title = "Transmission Line Geometry";

        // Determine display bounds using actual domain extent
        const maxY = Math.max(
            solver.dielectrics.reduce((max, d) => Math.max(max, d.y_max), 0),
            solver.conductors.reduce((max, c) => Math.max(max, c.y_max), 0)
        );

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
            if (zData.length > 0) {
                const flatZ = zData.flat();
                zMin = Math.min(...flatZ);
                zMax = Math.max(...flatZ);
            }
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

        // Limit display Y to domain extent
        const yArr = Array.from(solver.y);
        const maxY = yArr[ny - 1];
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
                zData.push(Array.from(V[i].slice(0, nx)));
            }
        }
        const flatZ = zData.flat();
        zMin = Math.min(...flatZ);
        zMax = Math.max(...flatZ);
    }

    else if ((currentView === "efield" || currentView === "efield_odd" || currentView === "efield_even") && solver.solution_valid) {
        // Ensure mesh exists for field visualization
        if (!solver.mesh_generated) {
            solver.ensure_mesh();
        }

        nx = solver.x.length;
        ny = solver.y.length;

        // Limit display Y to actual domain extent
        const yArr = Array.from(solver.y);
        const maxY = yArr[ny - 1];
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
        const flatZ = zData.flat();
        zMin = Math.min(...flatZ);
        zMax = Math.max(...flatZ);
    }

    else {
        title = "No Data Available";
        // Create minimal dummy data
        xMM = [0, (solver.w || 1) * 2000];
        yMM = [0, (solver.h || 1) * 1000];
    }

    // Save original mesh coordinates for mesh overlay before interpolation
    let xMM_mesh = xMM;
    let yMM_mesh = yMM;
    let nx_mesh = nx;
    let nyDisplay_mesh = nyDisplay;

    // Main field trace
    let traces = [];

    if (currentView === "geometry" && zData.length > 0) {
        const { Ex, Ey } = getFields();

        let eMax = Math.max(...zData.flat());
        let eMin = Math.min(...zData.flat());

        // Check if there's a user-defined scale override
        if (window.getStoredScale) {
            const override = window.getStoredScale(currentView);
            if (override) {
                eMin = override.min;
                eMax = override.max;
            }
        }

        const n = plotOptions.contours;
        const limits = contourScaledB(eMin, eMax, n);

        // Add E-field contours if requested
        if (plotOptions.contours > 0) {
            traces.push({
                type: "contour",
                x: xMM,
                y: yMM,
                contours: {
                    showlines: true,
                    coloring: "none",
                    start: limits[0],
                    end: limits[1],
                    size: limits[2]
                },
                z: zData.map(row => row.map(v => Math.log10(Math.max(v, 1e-3)))),
                line: {
                    smoothing: 1.3,
                    width: 1,
                    color: "rgba(0, 0, 0, 0.4)"
                },
                showscale: false,  // No colorbar
                name: "E-field contours",
                hoverinfo: "skip"
            });
        }

        // Add streamlines if requested via plot options
        if (plotOptions.streamlines > 0) {
            const modeIndex = getSelectedModeIndex();
            const mode = modeIndex === 1 ? 'even' : 'odd';

            traces.push(
                makeStreamlineTraceFromConductors(
                    Ex,
                    Ey,
                    solver.x,
                    solver.y,
                    solver.conductors,
                    plotOptions.streamlines,
                    mode
                )
            );
        }

    } else if (currentView === "geometry") {
        // Geometry only. Invisible scatter for axis scaling
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
        // Field views only. Use heatmap with optional contour lines

        // Check if there's a user-defined scale override
        if (window.getStoredScale) {
            const override = window.getStoredScale(currentView);
            if (override) {
                zMin = override.min;
                zMax = override.max;
            }
        }

        const n = plotOptions.contours;

        const contourSettings = {
            coloring: 'heatmap',
            showlines: n > 0
        };

        if (n > 0) {
            const step = (zMax - zMin) / n;
            contourSettings.start = zMin;
            contourSettings.end = zMax;
            contourSettings.size = step;
        }

        traces.push({
            type: n > 0 ? "contour" : "heatmap",
            zsmooth: "best",
            x: xMM,
            y: yMM,
            z: zData,
            zmin: zMin,
            zmax: zMax,
            colorscale: colorscale,
            contours: contourSettings,
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

    // Mesh overlay
    if (showMesh && solver.solution_valid) {
        const stepX = 1;
        const stepY = 1;

        // Use original mesh coordinates (before interpolation)
        for (let j = 0; j < nx_mesh; j += stepX) {
            traces.push({
                type: "scatter",
                x: [xMM_mesh[j], xMM_mesh[j]],
                y: [yMM_mesh[0], yMM_mesh[nyDisplay_mesh - 1]],
                mode: "lines",
                line: { width: 0.2, color: "black" },
                showlegend: false,
                hoverinfo: "skip"
            });
        }

        for (let i = 0; i < nyDisplay_mesh; i += stepY) {
            traces.push({
                type: "scatter",
                x: [xMM_mesh[0], xMM_mesh[nx_mesh - 1]],
                y: [yMM_mesh[i], yMM_mesh[i]],
                mode: "lines",
                line: { width: 0.2, color: "black" },
                showlegend: false,
                hoverinfo: "skip"
            });
        }
    }

    // UI menues
    const layout = {
        title: { text: title, font: { color: '#fff' } },
        xaxis: {
            title: { text: "Width (mm)", font: { color: '#aaa' } },
            scaleanchor: "y",
            scaleratio: 1,
            range: currentXRange,  // Preserve zoom/pan
            color: '#aaa',
            gridcolor: '#444',
            zerolinecolor: '#555'
        },
        yaxis: {
            title: { text: "Height (mm)", font: { color: '#aaa' } },
            range: currentYRange,  // Preserve zoom/pan
            color: '#aaa',
            gridcolor: '#444',
            zerolinecolor: '#555'
        },
        margin: { l: 70, r: 90, t: 50, b: 60 },
        showlegend: false,
        hovermode: "closest",
        dragmode: "pan",
        paper_bgcolor: '#2a2a2a',
        plot_bgcolor: '#1a1a1a',
        font: { color: '#fff' },
        shapes: shapes,  // Add vector shapes for geometry

        updatemenus: (() => {
            const menus = [];

            // View selector (Geometry/Potential/E-field)
            menus.push({
                x: 0.01,
                y: 1.15,
                showactive: true,
                active: (() => {
                    if (currentView === "geometry") return 0;
                    if (currentView === "potential") return 1;
                    if (currentView === "efield") return 2;
                    return 0;
                })(),
                bgcolor: '#2a2a2a',
                bordercolor: '#444',
                font: { color: '#aaa' },
                buttons: (() => {
                    const buttons = [
                        {
                            label: "Geometry",
                            method: "skip",
                            args: []
                        }
                    ];

                    if (solver.solution_valid) {
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

                    return buttons;
                })()
            });

            // Mode selector (Odd/Even) - only for differential lines
            if (isDifferentialMode()) {
                const modeIndex = getSelectedModeIndex();
                menus.push({
                    x: 0.25,
                    y: 1.15,
                    showactive: true,
                    active: modeIndex,
                    bgcolor: '#2a2a2a',
                    bordercolor: '#444',
                    font: { color: '#aaa' },
                    buttons: [
                        {
                            label: "Odd Mode",
                            method: "skip",
                            args: []
                        },
                        {
                            label: "Even Mode",
                            method: "skip",
                            args: []
                        }
                    ]
                });
            }

            return menus;
        })()
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
                name: "Scale Range",
                icon: Plotly.Icons.autoscale,
                click: () => window.toggleScaleDialog && window.toggleScaleDialog()
            }
        ]
    };

    Plotly.react(container, traces, layout, config);

    if (!container._viewListenerBound) {
        container.on('plotly_buttonclicked', (event) => {
            // Determine which menu was clicked based on x position
            // First menu (x=0.01): View selector (Geometry/Potential/E-field)
            // Second menu (x=0.25): Mode selector (Odd/Even) - only for differential

            if (event.menu.x < 0.2) {
                // View selector clicked
                if (event.menu.active === 0) {
                    currentView = "geometry";
                } else if (event.menu.active === 1) {
                    currentView = "potential";
                } else if (event.menu.active === 2) {
                    currentView = "efield";
                }
            } else {
                // Mode selector clicked (differential lines only)
                const plotModeEl = document.getElementById('plot-mode');
                if (plotModeEl) {
                    plotModeEl.value = event.menu.active === 0 ? 'odd' : 'even';
                }
            }
            draw();
        });
        container._viewListenerBound = true;
    }

}

function getYAxisLabel(selector) {
    const labels = {
        're_z0': 'Re(Z0) (Ohm)',
        'im_z0': 'Im(Z0) (Ohm)',
        'eps_eff': 'Effective permittivity',
        'loss': 'Loss (dB/m)',
        'R': 'R (Ohm/m)',
        'L': 'L (H/m)',
        'C': 'C (F/m)',
        'G': 'G (S/m)'
    };
    return labels[selector] || selector;
}

function drawResultsPlot() {
    const frequencySweepResults = get.frequencySweepResults();
    if (!frequencySweepResults || frequencySweepResults.length === 0) return;

    const selector = document.getElementById('results-plot-selector').value;
    const solver = get.solver();
    const isDifferential = solver && solver.is_differential;

    const freqs = frequencySweepResults.map(r => r.freq / 1e9);
    const traces = [];

    // Use lines+markers mode so single frequency points are visible
    const plotMode = freqs.length === 1 ? 'markers' : 'lines+markers';

    if (selector === 're_z0') {
        if (isDifferential) {
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[0].Zc.re),
                name: 'Odd mode',
                type: 'scatter',
                mode: plotMode
            });
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[1].Zc.re),
                name: 'Even mode',
                type: 'scatter',
                mode: plotMode
            });
        } else {
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[0].Zc.re),
                name: 'Re(Z0)',
                type: 'scatter',
                mode: plotMode
            });
        }
    } else if (selector === 'im_z0') {
        if (isDifferential) {
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[0].Zc.im),
                name: 'Odd mode',
                type: 'scatter',
                mode: plotMode
            });
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[1].Zc.im),
                name: 'Even mode',
                type: 'scatter',
                mode: plotMode
            });
        } else {
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[0].Zc.im),
                name: 'Im(Z0)',
                type: 'scatter',
                mode: plotMode
            });
        }
    } else if (selector === 'eps_eff') {
        if (isDifferential) {
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[0].eps_eff),
                name: 'Odd mode',
                type: 'scatter',
                mode: plotMode
            });
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[1].eps_eff),
                name: 'Even mode',
                type: 'scatter',
                mode: plotMode
            });
        } else {
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[0].eps_eff),
                name: 'eps_eff',
                type: 'scatter',
                mode: plotMode
            });
        }
    } else if (selector === 'loss') {
        if (isDifferential) {
            // Odd mode losses (solid lines)
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[0].alpha_c),
                name: 'Conductor (odd)',
                type: 'scatter',
                mode: plotMode
            });
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[0].alpha_d),
                name: 'Dielectric (odd)',
                type: 'scatter',
                mode: plotMode
            });
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[0].alpha_total),
                name: 'Total (odd)',
                type: 'scatter',
                mode: plotMode,
                line: { width: 2 }
            });
            // Even mode losses (dashed lines)
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[1].alpha_c),
                name: 'Conductor (even)',
                type: 'scatter',
                mode: plotMode,
                line: { dash: 'dash' }
            });
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[1].alpha_d),
                name: 'Dielectric (even)',
                type: 'scatter',
                mode: plotMode,
                line: { dash: 'dash' }
            });
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[1].alpha_total),
                name: 'Total (even)',
                type: 'scatter',
                mode: plotMode,
                line: { width: 2, dash: 'dash' }
            });
        } else {
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[0].alpha_c),
                name: 'Conductor',
                type: 'scatter',
                mode: plotMode
            });
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[0].alpha_d),
                name: 'Dielectric',
                type: 'scatter',
                mode: plotMode
            });
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[0].alpha_total),
                name: 'Total',
                type: 'scatter',
                mode: plotMode,
                line: { width: 2 }
            });
        }
    } else {
        // RLGC parameters
        const paramKey = selector; // R, L, G, or C
        if (isDifferential) {
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[0].RLGC[paramKey]),
                name: `${paramKey} (odd)`,
                type: 'scatter',
                mode: plotMode
            });
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[1].RLGC[paramKey]),
                name: `${paramKey} (even)`,
                type: 'scatter',
                mode: plotMode
            });
        } else {
            traces.push({
                x: freqs,
                y: frequencySweepResults.map(r => r.result.modes[0].RLGC[paramKey]),
                name: paramKey,
                type: 'scatter',
                mode: plotMode
            });
        }
    }

    const useLogX = document.getElementById('results-log-x').checked;
    const layout = {
        xaxis: {
            title: { text: 'Frequency (GHz)', font: { color: '#aaa' } },
            type: useLogX ? 'log' : 'linear',
            color: '#aaa',
            gridcolor: '#444',
            zerolinecolor: '#555'
        },
        yaxis: {
            title: { text: getYAxisLabel(selector), font: { color: '#aaa' } },
            color: '#aaa',
            gridcolor: '#444',
            zerolinecolor: '#555'
        },
        margin: { l: 80, r: 40, t: 40, b: 60 },
        showlegend: true,
        legend: { x: 0.02, y: 0.98, font: { color: '#fff' } },
        paper_bgcolor: '#2a2a2a',
        plot_bgcolor: '#1a1a1a',
        font: { color: '#fff' }
    };

    Plotly.newPlot('results-plot', traces, layout, { responsive: true });
}

function drawSParamPlot() {
    const frequencySweepResults = get.frequencySweepResults();
    if (!frequencySweepResults || frequencySweepResults.length === 0) return;

    const length = get.inputValue('sparam-length');
    const Z_ref = parseFloat(document.getElementById('sparam-z-ref').value);

    // Check for invalid inputs
    if (isNaN(length) || length <= 0 || isNaN(Z_ref) || Z_ref <= 0) {
        return;
    }

    const solver = get.solver();
    const isDifferential = solver && solver.is_differential;
    const plotMode = document.getElementById('sparam-plot-mode').value; // 'magnitude' or 'phase'

    const freqs = frequencySweepResults.map(r => r.freq / 1e9);
    const traces = [];

    // Use lines+markers mode so single frequency points are visible
    const lineMode = freqs.length === 1 ? 'markers' : 'lines+markers';

    // Helper to convert complex S-parameter to phase in degrees
    const sParamToPhase = (complexVal) => {
        return complexVal.arg() * 180 / Math.PI;
    };

    if (!isDifferential) {
        // 2-port S-parameters
        const S11_data = [];
        const S21_data = [];

        for (const { freq, result } of frequencySweepResults) {
            const sp = computeSParamsSingleEnded(freq, result.modes[0].RLGC, length, Z_ref);
            if (plotMode === 'magnitude') {
                S11_data.push(sParamTodB(sp.S11));
                S21_data.push(sParamTodB(sp.S21));
            } else {
                S11_data.push(sParamToPhase(sp.S11));
                S21_data.push(sParamToPhase(sp.S21));
            }
        }

        const label = plotMode === 'magnitude' ? '(dB)' : '(deg)';
        traces.push({
            x: freqs,
            y: S11_data,
            name: `S11 ${label}`,
            type: 'scatter',
            mode: lineMode
        });
        traces.push({
            x: freqs,
            y: S21_data,
            name: `S21 ${label}`,
            type: 'scatter',
            mode: lineMode
        });
    } else {
        // 4-port S-parameters (mixed-mode)
        const SDD11_data = [];
        const SDD21_data = [];
        const SCC11_data = [];
        const SCC21_data = [];

        for (const { freq, result } of frequencySweepResults) {
            const oddMode = result.modes.find(m => m.mode === 'odd');
            const evenMode = result.modes.find(m => m.mode === 'even');

            const sp = computeSParamsDifferential(
                freq,
                oddMode.RLGC,
                evenMode.RLGC,
                length,
                Z_ref
            );

            if (plotMode === 'magnitude') {
                SDD11_data.push(sParamTodB(sp.SDD11));
                SDD21_data.push(sParamTodB(sp.SDD21));
                SCC11_data.push(sParamTodB(sp.SCC11));
                SCC21_data.push(sParamTodB(sp.SCC21));
            } else {
                SDD11_data.push(sParamToPhase(sp.SDD11));
                SDD21_data.push(sParamToPhase(sp.SDD21));
                SCC11_data.push(sParamToPhase(sp.SCC11));
                SCC21_data.push(sParamToPhase(sp.SCC21));
            }
        }

        const label = plotMode === 'magnitude' ? '(dB)' : '(deg)';
        traces.push({
            x: freqs,
            y: SDD11_data,
            name: `SDD11 ${label}`,
            type: 'scatter',
            mode: lineMode
        });
        traces.push({
            x: freqs,
            y: SDD21_data,
            name: `SDD21 ${label}`,
            type: 'scatter',
            mode: lineMode
        });
        traces.push({
            x: freqs,
            y: SCC11_data,
            name: `SCC11 ${label}`,
            type: 'scatter',
            mode: lineMode,
            line: { dash: 'dash' }
        });
        traces.push({
            x: freqs,
            y: SCC21_data,
            name: `SCC21 ${label}`,
            type: 'scatter',
            mode: lineMode,
            line: { dash: 'dash' }
        });
    }

    const useLogX = document.getElementById('sparam-log-x').checked;
    const yTitle = plotMode === 'magnitude' ? 'Magnitude (dB)' : 'Phase (degrees)';
    const layout = {
        xaxis: {
            title: { text: 'Frequency (GHz)', font: { color: '#aaa' } },
            type: useLogX ? 'log' : 'linear',
            color: '#aaa',
            gridcolor: '#444',
            zerolinecolor: '#555'
        },
        yaxis: {
            title: { text: yTitle, font: { color: '#aaa' } },
            color: '#aaa',
            gridcolor: '#444',
            zerolinecolor: '#555'
        },
        margin: { l: 80, r: 40, t: 40, b: 60 },
        showlegend: true,
        legend: { x: 0.02, y: 0.02, font: { color: '#fff' } },
        paper_bgcolor: '#2a2a2a',
        plot_bgcolor: '#1a1a1a',
        font: { color: '#fff' }
    };

    Plotly.newPlot('sparam-plot', traces, layout, { responsive: true });
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

// Helper function to check if solver is in differential mode
function isDifferentialMode() {
    const solver = get.solver();
    if (!solver || !solver.Ex || !solver.Ey) return false;
    // In differential mode, Ex and Ey are arrays of 2 arrays (odd and even modes)
    // Check if Ex[0] and Ex[1] are both arrays
    return Array.isArray(solver.Ex) &&
           solver.Ex.length === 2 &&
           Array.isArray(solver.Ex[0]) &&
           Array.isArray(solver.Ex[1]);
}

// Get the selected mode index from sidebar (0=odd, 1=even)
function getSelectedModeIndex() {
    const modeSelect = document.getElementById('plot-mode');
    return modeSelect && modeSelect.value === 'even' ? 1 : 0;
}

// Helper function to get Ex/Ey fields (handles differential mode)
function getFields() {
    const solver = get.solver();
    if (!solver || !solver.Ex || !solver.Ey) {
        return { Ex: null, Ey: null };
    }

    if (isDifferentialMode()) {
        const modeIndex = getSelectedModeIndex();
        return { Ex: solver.Ex[modeIndex], Ey: solver.Ey[modeIndex] };
    } else {
        // Single-ended mode
        return { Ex: solver.Ex[0], Ey: solver.Ey[0] };
    }
}

// Helper function to get voltage potential (handles differential mode)
function getPotential() {
    const solver = get.solver();
    if (!solver || !solver.V) {
        return null;
    }

    if (isDifferentialMode()) {
        const modeIndex = getSelectedModeIndex();
        return solver.V[modeIndex];
    } else {
        // Single-ended mode
        return solver.V[0];
    }
}

// Get plot options from sidebar
function getPlotOptions() {
    const streamlinesEl = document.getElementById('plot-streamlines');
    const contoursEl = document.getElementById('plot-contours');

    const streamlinesVal = streamlinesEl ? streamlinesEl.value.trim() : '';
    const contoursVal = contoursEl ? contoursEl.value.trim() : '';

    return {
        streamlines: streamlinesVal === '' ? 0 : parseInt(streamlinesVal) || 0,
        contours: contoursVal === '' ? 0 : parseInt(contoursVal) || 0
    };
}

// Function to set current view
function setCurrentView(view) {
    currentView = view;
    // Notify app.js that view changed so it can restore the appropriate scale
    if (window.onViewChanged) {
        window.onViewChanged(view);
    }
}

export { draw, drawResultsPlot, drawSParamPlot, setGlobals, setCurrentView, getScaleRange, setScaleRange };
