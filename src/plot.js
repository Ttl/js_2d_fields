import { makeStreamlineTraceFromConductors } from './streamlines.js';
import { computeSParamsSingleEnded, computeSParamsDiffAuto, sParamTodB,
         isSelfReferenced, sparamsForPoint, usableSweepPoints } from './sparameters.js';
import { isComplement, svgRingPath, svgShapePath, shapePoly, isPolyShape, shapeContains } from './shapes.js';

// A polygon or ring shape as a Plotly path shape with `style`.
const polyPathShape = (shape, style) => ({ type: 'path', path: svgShapePath(shape), fillrule: 'evenodd', ...style });

// Outline of a conductor or dielectric in `style`: its polygon, or its rectangle cut at maxY.
const bodyShape = (o, maxY, style) => (isPolyShape(o.shape) ? polyPathShape(o.shape, style) : {
    type: 'rect', x0: o.x_min * 1000, y0: o.y_min * 1000, x1: o.x_max * 1000, y1: Math.min(o.y_max, maxY) * 1000, ...style });

// Plated faces of a polygon shape as gold lines (every edge when it names none).
function platedEdgeLines(shape, plating, line) {
    const poly = shapePoly(shape);
    const n = poly.length >> 1;
    const out = [];
    const loops = shape.type === 'ring' ? [shape.poly, shape.hole] : [poly];
    loops.forEach((p, li) => {
        const m = p.length >> 1;
        for (let i = 0; i < m; i++) {
            const face = (li === 0 && shape.faces && m === n) ? shape.faces[i] : 'all';
            if (face !== 'all' && !plating[face]) continue;
            const j = (i + 1) % m;
            out.push({ type: 'line', x0: p[2 * i] * 1000, y0: p[2 * i + 1] * 1000, x1: p[2 * j] * 1000, y1: p[2 * j + 1] * 1000,
                line, layer: 'above' });
        }
    });
    return out;
}

// Lazy Plotly access - allows app to function while Plotly is loading
const getPlotly = () => window.Plotly;

let showMesh = false;
let currentView = "geometry";
// The view the user picked. currentView is it while its data exists, else the geometry.
let wantedView = "geometry";
let zMin = null;
let zMax = null;
// Store actual data range (before any user scaling)
let actualDataMin = null;
let actualDataMax = null;
// Largest sample, above the autoscale max when a corner singularity is clipped
let actualDataPeak = null;
// Lower |E| bound for the log-spaced contours (smallest field outside the conductors)
let contourFloor = 0;

// Geometry view zoom constants
const SIGNAL_CONDUCTOR_VIEW_FRACTION = 1/3;  // Signal conductors take up this fraction of X-axis view

// Frozen trace state
let frozenResultsData = null;   // Deep copy of frequencySweepResults
let frozenSParamData = null;    // { results: deepCopy, length, zRef }

// Plotly default color cycle (colorway)
const PLOTLY_COLORS = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'
];

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
    const logStep = n === 0 ? 1 : Math.abs((eMax - eMin)) / n;
    return [eMin, eMax, logStep];
}

// Contour levels in log10(|E|) for a dB color range, one per n-th of the range.
function contourLimitsDb(minDb, maxDb, n) {
    const lo = minDb / 20, hi = Math.max(maxDb / 20, lo + 0.005);
    return [lo, hi, n === 0 ? 1 : (hi - lo) / n];
}

// Area-weighted quantiles of the nonzero |E| samples. The plot grid is graded, dense at
// conductor corners, so plain sample quantiles over-weight the singular corner field.
function fieldQuantiles(zData, xMM, yMM, qs) {
    const pts = [];
    const nyz = zData.length, nxz = xMM.length;
    for (let i = 0; i < nyz; i++) {
        const row = zData[i];
        if (!row || !row.length) continue;
        const dy = (yMM[Math.min(i + 1, nyz - 1)] - yMM[Math.max(i - 1, 0)]) / 2;
        for (let j = 0; j < nxz; j++) {
            const v = row[j];
            if (!(v > 0)) continue;
            const dx = (xMM[Math.min(j + 1, nxz - 1)] - xMM[Math.max(j - 1, 0)]) / 2;
            pts.push([v, dx * dy]);
        }
    }
    if (!pts.length) return qs.map(() => 0);
    pts.sort((a, b) => a[0] - b[0]);
    let total = 0;
    for (const p of pts) total += p[1];
    const out = [];
    let k = 0, acc = 0;
    for (const q of qs) {
        while (k < pts.length - 1 && acc + pts[k][1] < q * total) acc += pts[k++][1];
        out.push(pts[k][0]);
    }
    return out;
}

// Autoscale of an |E| grid: the linear max is the area-weighted 99.99th percentile so a
// singular corner does not leave the rest of the plot dark; peak is the true maximum.
function efieldAutoscale(zData, xMM, yMM) {
    const [floor, max] = fieldQuantiles(zData, xMM, yMM, [0.001, 0.9999]);
    let peak = 0;
    for (const row of zData) for (const v of row) if (v > peak) peak = v;
    return { floor, max: max || peak, peak };
}

const toDb = (v) => 20 * Math.log10(v);

// dB color range: from the peak down to the weakest field, at most 60 dB.
function efieldDbRange(auto) {
    const max = Math.ceil(toDb(auto.peak || 1));
    const min = Math.floor(Math.max(max - 60, auto.floor > 0 ? toDb(auto.floor) : max - 60));
    return { min, max: Math.max(max, min + 1) };
}

// Closed rectangular loop as an SVG path, in mm. Two of these in one path with the evenodd
// fill rule give a rectangular ring, the frame drawn around an enclosed medium's domain.
function rectLoopPath(x0, y0, x1, y1) {
    const m = (v) => v * 1000;
    return `M ${m(x0)},${m(y0)} L ${m(x1)},${m(y0)} L ${m(x1)},${m(y1)} L ${m(x0)},${m(y1)} Z`;
}

// Opaque conductor fills (+ yellow plating edges), drawn above the field. Shared by the geometry
// view and the field views so conductors look identical and the field views don't show heatmap/
// contour bleed inside the PEC (the resampled field is identically zero in the conductor interior,
// this just masks the zsmooth/contour interpolation that spills the steep boundary field inward).
function conductorFillShapes(solver, maxY) {
    const out = [];
    const FILL = 'rgba(217, 119, 6, 1.0)';
    const EDGE = { color: 'rgba(0, 0, 0, 0.5)', width: 1 };
    const GOLD = { color: 'rgba(255, 215, 0, 1.0)', width: 3 };

    // Enclosing metal walls of a source-free medium (rectangular waveguide). They are not
    // conductors in `solver.conductors`, physically they live in the boundary conditions,
    // with no meshed thickness, so without this the geometry view would render a bare
    // rectangle of dielectric with nothing to show it is a waveguide. Drawn as a ring
    // (outer loop + inner loop, evenodd) exactly like the coax shield.
    const w = solver.enclosure_walls;
    if (w) {
        const t = w.thickness;
        out.push({
            type: 'path',
            path: rectLoopPath(w.x_min - t, w.y_min - t, w.x_max + t, w.y_max + t) + ' ' +
                  rectLoopPath(w.x_min, w.y_min, w.x_max, w.y_max),
            fillrule: 'evenodd', fillcolor: FILL, line: EDGE, layer: 'above',
        });
        // Plating covers the whole inner surface (one continuous wall, no separate faces),
        // so the indicator is an outline of that surface rather than per-face lines.
        if (solver.plating) {
            out.push({
                type: 'path', path: rectLoopPath(w.x_min, w.y_min, w.x_max, w.y_max),
                fillcolor: 'rgba(0,0,0,0)', line: GOLD, layer: 'above',
            });
        }
    }
    for (const cond of (solver.conductors || [])) {
        const sh = cond.shape;
        if (isPolyShape(sh)) {
            out.push(polyPathShape(sh, { fillcolor: FILL, line: EDGE, layer: 'above' }));
            if (cond.plating) out.push(...platedEdgeLines(sh, cond.plating, GOLD));
            continue;
        }
        if (sh) {
            // A round conductor is drawn with Plotly's ellipse shape. The enclosing
            // shield is an annulus, which needs an SVG path because Plotly's shape path
            // grammar has no arc command (M/L/H/V/Q/C/T/S/Z only), svgRingPath emits
            // two polygonal loops filled with the evenodd rule.
            const cx = sh.cx * 1000, cy = sh.cy * 1000, r = sh.r * 1000;
            if (isComplement(sh)) {
                out.push({
                    type: 'path', path: svgRingPath(cx, cy, r, cond.x_max * 1000),
                    fillrule: 'evenodd', fillcolor: FILL, line: EDGE, layer: 'above',
                });
            } else {
                out.push({
                    type: 'circle', xref: 'x', yref: 'y',
                    x0: cx - r, y0: cy - r, x1: cx + r, y1: cy + r,
                    fillcolor: FILL, line: EDGE, layer: 'above',
                });
            }
            // Plating covers the whole circumference (a circle has no separate faces),
            // so the indicator is a gold outline rather than per-face lines.
            if (cond.plating) {
                out.push({
                    type: 'circle', xref: 'x', yref: 'y',
                    x0: cx - r, y0: cy - r, x1: cx + r, y1: cy + r,
                    fillcolor: 'rgba(0,0,0,0)', line: GOLD, layer: 'above',
                });
            }
            continue;
        }
        if (cond.y_min > maxY) continue;
        const yMax = Math.min(cond.y_max, maxY);
        out.push({
            type: 'rect',
            x0: cond.x_min * 1000, y0: cond.y_min * 1000,
            x1: cond.x_max * 1000, y1: yMax * 1000,
            fillcolor: 'rgba(217, 119, 6, 1.0)',
            line: { color: 'rgba(0, 0, 0, 0.5)', width: 1 },
            layer: 'above'
        });
        if (cond.plating) {   // yellow lines on plated edges
            const x0 = cond.x_min * 1000, x1 = cond.x_max * 1000, y0 = cond.y_min * 1000, y1 = yMax * 1000;
            const plateLine = { color: 'rgba(255, 215, 0, 1.0)', width: 3 };
            if (cond.plating.top) out.push({ type: 'line', x0, y0: y1, x1, y1: y1, line: plateLine, layer: 'above' });
            if (cond.plating.bottom) out.push({ type: 'line', x0, y0: y0, x1, y1: y0, line: plateLine, layer: 'above' });
            if (cond.plating.sides) {
                out.push({ type: 'line', x0: x0, y0: y0, x1: x0, y1: y1, line: plateLine, layer: 'above' });
                out.push({ type: 'line', x0: x1, y0: y0, x1: x1, y1: y1, line: plateLine, layer: 'above' });
            }
        }
    }
    return out;
}

// Dielectric rect shapes, colored by ε_r (air ≈1 → white/transparent, higher ε_r → green
// shades). Shared by the geometry view (opaque, below the contours) and the Modes tab
// (faint, above the field heatmap) so the two tabs use the same color mapping.
function dielectricFillShapes(solver, maxY, { alpha = 0.8, airAlpha = alpha, layer = 'below',
    lineColor = 'rgba(128, 128, 128, 0.3)' } = {}) {
    const out = [];
    for (const diel of (solver.dielectrics || [])) {
        if (!diel.shape && diel.y_min > maxY) continue;
        const yMax = Math.min(diel.y_max, maxY);
        const er = diel.epsilon_r;
        let fillcolor;
        if (er <= 1.01) {
            fillcolor = `rgba(255, 255, 255, ${airAlpha})`;
        } else {
            const intensity = Math.min(255, 100 + (er - 1) * 30);
            fillcolor = `rgba(100, ${intensity}, 100, ${alpha})`;
        }
        const sh = diel.shape;
        if (isPolyShape(sh)) {
            out.push(polyPathShape(sh, { fillcolor, line: { color: lineColor, width: 0.5 }, layer }));
            continue;
        }
        if (sh && !isComplement(sh)) {
            const cx = sh.cx * 1000, cy = sh.cy * 1000, r = sh.r * 1000;
            out.push({
                type: 'circle', xref: 'x', yref: 'y',
                x0: cx - r, y0: cy - r, x1: cx + r, y1: cy + r,
                fillcolor, line: { color: lineColor, width: 0.5 }, layer
            });
            continue;
        }
        out.push({
            type: 'rect',
            x0: diel.x_min * 1000, y0: diel.y_min * 1000,
            x1: diel.x_max * 1000, y1: yMax * 1000,
            fillcolor, line: { color: lineColor, width: 0.5 }, layer
        });
    }
    return out;
}

// Top of what the geometry and field plots show: the highest rectangle, which for the
// fixed line types is their air box. A custom geometry has no air rectangle of its own,
// its solved region is the air, so the plots run to the top of that.
function displayTop(solver) {
    if (solver.user_domain) return solver.domain_height;
    return Math.max(
        solver.dielectrics.reduce((max, d) => Math.max(max, d.y_max), 0),
        solver.conductors.reduce((max, c) => Math.max(max, c.y_max), 0));
}

// Outline of the rectangles written on the text line the custom geometry editor's cursor
// is on (window.customHighlightLine, 0 for none).
function sourceLineHighlightShapes(solver, maxY) {
    const line = window.customHighlightLine;
    if (!line) return [];
    const style = { fillcolor: 'rgba(56, 189, 248, 0.25)', line: { color: 'rgba(56, 189, 248, 1)', width: 2 }, layer: 'above' };
    return [...(solver.dielectrics || []), ...(solver.conductors || [])]
        .filter(o => o.src_line === line && o.y_min <= maxY)
        .map(o => bodyShape(o, maxY, style));
}

// Rectangles generated by mirror=1 in a custom geometry: a dashed outline on each (no
// fill, which would read as another material in the permittivity shading), and the
// mirror plane x=0 as a dash-dot line.
function mirrorImageShapes(solver, maxY) {
    const images = [...(solver.dielectrics || []), ...(solver.conductors || [])]
        .filter(o => o.src_image && o.y_min <= maxY);
    if (!images.length) return [];
    const u = solver.user_domain;
    const x0 = -(solver.x_shift || 0) * 1000;
    const style = { fillcolor: 'rgba(0, 0, 0, 0)', line: { color: 'rgba(60, 60, 60, 0.7)', width: 1, dash: 'dash' }, layer: 'above' };
    return [
        ...images.map(o => bodyShape(o, maxY, style)),
        { type: 'line', x0, x1: x0, y0: (u ? u.y_min : solver.domain_y_min) * 1000, y1: maxY * 1000,
          line: { color: 'rgba(56, 189, 248, 0.7)', width: 1, dash: 'dashdot' }, layer: 'above' },
    ];
}

// Focused view (mm) around the signal conductors: the signal cluster fills `fraction` of
// the x-axis; with a top ground the full stack height is shown, otherwise the conductors
// sit in the bottom `fraction` of the view. Shared by the geometry tab's initial zoom and
// the Modes tab so both frame the structure identically. Returns null when there are no
// signal conductors to frame (caller picks its own fallback).
function computeGeometryView(solver, maxY, fraction = SIGNAL_CONDUCTOR_VIEW_FRACTION) {
    const signal = solver.conductors.filter(c => c.is_signal);
    const grounds = solver.conductors.filter(c => !c.is_signal);
    if (!signal.length) {
        // A source-free enclosed medium (rectangular waveguide) has no signal cluster to
        // centre on, the structure is the domain. Frame the whole cross-section walls
        // included with a margin. Without this the caller falls back to Plotly autoscale,
        // which ignores shapes and so leaves the guide off-centre in a default range.
        const w = solver.enclosure_walls;
        if (!w) return null;
        const outer = w.thickness;
        const pad = 0.10 * Math.max(w.x_max - w.x_min, w.y_max - w.y_min);
        return {
            xRange: [(w.x_min - outer - pad) * 1000, (w.x_max + outer + pad) * 1000],
            yRange: [(w.y_min - outer - pad) * 1000, (w.y_max + outer + pad) * 1000],
        };
    }
    // A signal edge on the domain wall (a slotline half plane) says nothing about where
    // the fields are, so the view is framed by the other signal edges and the ground
    // edges facing them.
    const W2 = (solver.domain_width || 0) / 2, tolW = W2 * 1e-9;
    const inner = v => !(W2 > 0) || Math.abs(Math.abs(v) - W2) > tolW;
    const sx0 = Math.min(...signal.map(c => c.x_min)), sx1 = Math.max(...signal.map(c => c.x_max));
    const sy0 = Math.min(...signal.map(c => c.y_min)), sy1 = Math.max(...signal.map(c => c.y_max));
    // A ground around the signals (coax or twinax shield) holds the fields: frame it.
    const shield = grounds.find(g => g.x_min <= sx0 && g.x_max >= sx1 && g.y_min <= sy0 && g.y_max >= sy1);
    if (shield) {
        const pad = 0.10 * Math.max(shield.x_max - shield.x_min, shield.y_max - shield.y_min);
        // The margin stops at the edge of a custom geometry's solved region.
        const u = solver.user_domain || { x_min: -Infinity, x_max: Infinity, y_min: -Infinity, y_max: Infinity };
        return {
            xRange: [Math.max(shield.x_min - pad, u.x_min) * 1000, Math.min(shield.x_max + pad, u.x_max) * 1000],
            yRange: [Math.max(shield.y_min - pad, u.y_min) * 1000, Math.min(shield.y_max + pad, u.y_max) * 1000],
        };
    }
    let xs = signal.flatMap(c => [c.x_min, c.x_max]).filter(inner);
    if (xs.length < 2) xs = xs.concat(grounds.flatMap(c => [c.x_min, c.x_max]).filter(inner));
    if (xs.length < 2) xs = signal.flatMap(c => [c.x_min, c.x_max]);
    // A return strip no wider than the signal cluster (coplanar strips) belongs to the
    // structure, wide grounds (CPW, planes) are framed by the signals alone.
    const clusterW = Math.max(...xs) - Math.min(...xs);
    for (const g of grounds) {
        if (inner(g.x_min) && inner(g.x_max) && g.x_max - g.x_min <= clusterW) xs.push(g.x_min, g.x_max);
    }
    let xl = Math.min(...xs), xr = Math.max(...xs);
    if (!(xr > xl)) { xl -= W2 / 20; xr += W2 / 20; }
    const center = (xl + xr) / 2;
    const viewWidth = (xr - xl) / fraction;
    const xRange = [(center - viewWidth / 2) * 1000, (center + viewWidth / 2) * 1000];

    const bottomY = grounds.length ? Math.min(...grounds.map(g => g.y_min)) : 0;
    const hasTopGround = grounds.some(c => c.y_max >= maxY * 0.9);
    let yRange;
    if (hasTopGround) {
        yRange = [bottomY * 1000, maxY * 1000];
    } else {
        const topOfConductors = Math.max(...solver.conductors.map(c => c.y_max));
        const viewHeight = (topOfConductors - bottomY) / fraction;
        // Air below the lowest ground (a custom geometry with an open bottom boundary):
        // show part of it, the structure is not sitting on the edge of the solved region.
        const airBelow = solver.user_domain && solver.user_domain.y_min < bottomY - 1e-9 * viewHeight;
        const yLow = airBelow ? Math.max(solver.user_domain.y_min, bottomY - viewHeight / 3) : bottomY;
        yRange = [yLow * 1000, (bottomY + viewHeight) * 1000];
    }
    return { xRange, yRange };
}

// The log-spaced |E| contour-LINE trace (lines only, no fill), shared by the geometry overlay and
// the |E| field view so their contours are identical. z = log10(|E|) with log-spaced levels keeps
// the lines evenly spaced instead of crowding at the singular trace corners. Named "E-field
// contours" so setScaleRange can rescale it live in either view.
function efieldContourTrace(xMM, yMM, zData, limits) {
    return {
        type: "contour",
        x: xMM, y: yMM,
        z: zData.map(row => row.map(v => Math.log10(Math.max(v, 1e-3)))),
        contours: { showlines: true, coloring: "none", start: limits[0], end: limits[1], size: limits[2] },
        line: { smoothing: 1.3, width: 1, color: "rgba(0, 0, 0, 0.4)" },
        showscale: false,
        name: "E-field contours",
        hoverinfo: "skip"
    };
}

// Export functions to get/set scale range for current view
function getScaleRange() {
    return { min: zMin, max: zMax, view: scaleView() };
}

// The E-field and surface current views in dB are separate scales: their range is in dB.
function isDbView() {
    return (currentView.startsWith("efield") || currentView === "current" || currentView === "density")
        && getPlotOptions().efieldDb;
}
function scaleView() {
    return isDbView() ? currentView + "_db" : currentView;
}

// Get actual data range (before any user scaling)
function getActualDataRange() {
    return { min: actualDataMin, max: actualDataMax, peak: actualDataPeak };
}

function setScaleRange(min, max) {
    zMin = min;
    zMax = max;
    // The surface current colors are binned by value, and the triangle-mesh contours
    // are traced per level, so a new range redraws them.
    if (currentView === "current" || currentView === "density"
        || ((currentView.startsWith("efield") || currentView === "geometry") && getFieldMesh())) {
        draw(); return;
    }

    const container = document.getElementById('sim_canvas');
    const Plotly = getPlotly();
    if (!container || !container.data || !Plotly) return;

    const n = getPlotOptions().contours;

    // The shared log-spaced "E-field contours" line trace appears in BOTH the geometry overlay
    // and the |E| field view — rescale its log levels identically wherever it is.
    const cIdx = container.data.findIndex(t => t.type === 'contour' && t.name === 'E-field contours');
    if (cIdx !== -1 && n > 0) {
        const limits = isDbView() ? contourLimitsDb(min, max, n) : contourScaledB(min, max, n);
        Plotly.restyle(container, {
            'contours.start': limits[0], 'contours.end': limits[1], 'contours.size': limits[2]
        }, [cIdx]);
    }

    // Field views also carry the color in a heatmap (|E|) or a linear contour (potential) — update
    // its zmin/zmax (and linear levels for the potential contour). The geometry view has no such trace.
    if (currentView !== "geometry") {
        const hIdx = container.data.findIndex(t =>
            t.type === 'heatmap' || (t.type === 'contour' && t.name !== 'E-field contours'));
        if (hIdx !== -1) {
            const restyle = { zmin: min, zmax: max };
            if (container.data[hIdx].type === 'contour' && n > 0) {
                restyle['contours.start'] = min;
                restyle['contours.end'] = max;
                restyle['contours.size'] = (max - min) / n;
            }
            Plotly.restyle(container, restyle, [hIdx]);
        }
    }
}

// The solver's plot grid in mm, for the mesh overlay of the views with no field grid
// of their own (|K|, |J|).
function solverGridMM(solver) {
    if (!solver.x || !solver.y) return { xMM: [0, 1], yMM: [0, 1], nx: 0, nyDisplay: 0 };
    return { xMM: Array.from(solver.x, v => v * 1000), yMM: Array.from(solver.y, v => v * 1000),
             nx: solver.x.length, nyDisplay: solver.y.length };
}

// Help topic (field_solver.html helpContent) of the view on screen.
function plotHelpTopic() {
    return currentView.startsWith("potential") ? "plot_potential"
        : currentView.startsWith("efield") ? "plot_efield"
        : currentView === "current" ? "plot_current"
        : currentView === "density" ? "plot_density" : "plot_geometry";
}

// Modebar icon for the mesh toggle, a 3x3 grid (Plotly has no grid icon).
const GRID_ICON = {
    width: 1000, height: 1000,
    path: [0, 460, 920].map(p =>
        `M${p} 0h80v1000h-80z M0 ${p}h1000v80h-1000z`).join(' '),
};

// Title suffix naming the frequency and kind of the plotted fields.
function fieldFreqLabel(solver) {
    const f = solver.fieldFreq;
    if (typeof f !== "number" || !(f >= 0)) return "";
    const fs = f >= 1e9 ? `${+(f / 1e9).toPrecision(4)} GHz` : f >= 1e6 ? `${+(f / 1e6).toPrecision(4)} MHz`
        : f >= 1e3 ? `${+(f / 1e3).toPrecision(4)} kHz` : `${+f.toPrecision(4)} Hz`;
    // Below a waveguide cutoff the mode keeps its transverse pattern but decays along z.
    const evanescent = solver.fc > 0 && f <= solver.fc ? ", below cutoff (evanescent)" : "";
    return ` at ${fs}` + (solver.fieldKind === 'fullwave' && currentView.startsWith("efield") ? ", full-wave mode" : "") + evanescent;
}

// Surface current segments of the displayed mode, or null.
function getSurfaceK() {
    const solver = get.solver();
    const K = solver && solver.surfaceK;
    if (!K) return null;
    return K[isDifferentialMode() ? getSelectedModeIndex() : 0] || null;
}

// |E| on the triangles of the displayed mode (meshFieldBlock), or null: the triangular
// backend's plot fields, drawn as an image and traced contours instead of the grid.
function getFieldMesh() {
    const solver = get.solver();
    const M = solver && solver.fieldMesh;
    if (!M) return null;
    return M[isDifferentialMode() ? getSelectedModeIndex() : 0] || null;
}

// Surface current of the displayed mode's MQS solve on the metal with no |J| block: the
// walls and the grounds absorbed into them (a surface impedance), and the grounds held
// as ideal returns (see _mqsSolve), or null. A segment belongs to them when its
// midpoint lies on no conductor with a block.
function bareMetalK(solver, blocks) {
    const mi = isDifferentialMode() ? getSelectedModeIndex() : 0;
    const K = getSurfaceK();
    if (!K || !solver.surfaceKSource || solver.surfaceKSource[mi] !== 'mqs') return null;
    const withJ = (solver.conductors || []).filter(hasDensityBlock(blocks));
    const tol = 1e-9 * (solver.domain_width || 1);
    const out = { x0: [], x1: [], y0: [], y1: [], K: [] };
    for (let i = 0; i < K.K.length; i++) {
        const mx = (K.x0[i] + K.x1[i]) / 2, my = (K.y0[i] + K.y1[i]) / 2;
        if (withJ.some(c => shapeContains(c, mx, my, tol))) continue;
        out.x0.push(K.x0[i]); out.x1.push(K.x1[i]); out.y0.push(K.y0[i]); out.y1.push(K.y1[i]); out.K.push(K.K[i]);
    }
    return out.K.length ? out : null;
}

// Current density blocks of the displayed mode, or null.
function getCurrentJ() {
    const solver = get.solver();
    const J = solver && solver.currentJ;
    if (!J) return null;
    return J[isDifferentialMode() ? getSelectedModeIndex() : 0] || null;
}

// The triangle mesh of the displayed view: the skin mesh of the MQS solve for the
// current views that come from it, else the solve mesh (null on the rectilinear backend).
function viewTriMesh(solver) {
    const mi = isDifferentialMode() ? getSelectedModeIndex() : 0;
    const fromMqs = currentView === "density"
        || (currentView === "current" && solver.surfaceKSource && solver.surfaceKSource[mi] === 'mqs');
    return (fromMqs && solver.currentMesh && solver.currentMesh[mi]) || solver.triMesh || null;
}

// Conductor shapes of the |J| view. A conductor under a current density block gets only
// an outline above the heatmap: its fill (and plating edges) below would show through
// at the block edges as a false hot rim. Metal without a block (walls, wall grounds)
// keeps its fill.
function densityConductorShapes(solver, blocks, maxY) {
    const covered = hasDensityBlock(blocks);
    const conductors = solver.conductors || [];
    const rest = Object.create(solver);
    rest.conductors = conductors.filter(c => !covered(c));
    const out = conductorFillShapes(rest, maxY).map(sh => ({ ...sh, layer: 'below' }));
    const OUTLINE = { fillcolor: 'rgba(0,0,0,0)', line: { color: 'rgba(255, 255, 255, 0.5)', width: 1 }, layer: 'above' };
    for (const c of conductors.filter(covered)) {
        const sh = c.shape;
        if (sh && !isPolyShape(sh)) {
            const cx = sh.cx * 1000, cy = sh.cy * 1000, r = sh.r * 1000;
            out.push({ type: 'circle', xref: 'x', yref: 'y', x0: cx - r, y0: cy - r, x1: cx + r, y1: cy + r, ...OUTLINE });
        } else {
            out.push(bodyShape(c, maxY, OUTLINE));
        }
    }
    return out;
}

// Test of whether a conductor lies under one of the |J| blocks. A complement shell (coax
// shield) has no cross-section and never gets one.
function hasDensityBlock(blocks) {
    const bb = blocks.map(blockBox);
    return c => !(c.shape && isComplement(c.shape)) && bb.some(b => Math.min(b.x1, c.x_max) > Math.max(b.x0, c.x_min)
        && Math.min(b.y1, c.y_max) > Math.max(b.y0, c.y_min));
}

// Bounding box { x0, x1, y0, y1 } of a |J| block, grid or triangles.
function blockBox(b) {
    if (b.tris) {
        const o = { x0: Infinity, x1: -Infinity, y0: Infinity, y1: -Infinity };
        for (let i = 0; i < b.tris.length; i += 2) {
            o.x0 = Math.min(o.x0, b.tris[i]); o.x1 = Math.max(o.x1, b.tris[i]);
            o.y0 = Math.min(o.y0, b.tris[i + 1]); o.y1 = Math.max(o.y1, b.tris[i + 1]);
        }
        return o;
    }
    return { x0: Math.min(b.x[0], b.x[b.x.length - 1]), x1: Math.max(b.x[0], b.x[b.x.length - 1]),
             y0: b.y[0], y1: b.y[b.y.length - 1] };
}

// Hover and color axis of the triangle blocks (shaped conductors): invisible markers at
// the centroids of at most HOVER_MAX triangles, a marker per triangle (~1e5 on a skin
// mesh) makes every pan and zoom slow. The |J| itself is the image of rasterizeDensity.
const HOVER_MAX = 5000;
function densityHoverTrace(blocks, db) {
    const mx = [], my = [], mv = [];
    const n = blocks.reduce((a, b) => a + b.J.length, 0);
    const stride = Math.max(1, Math.ceil(n / HOVER_MAX));
    for (const b of blocks) {
        const T = b.tris;
        for (let t = 0; t < b.J.length; t += stride) {
            const k = 6 * t;
            mx.push((T[k] + T[k + 2] + T[k + 4]) * 1000 / 3); my.push((T[k + 1] + T[k + 3] + T[k + 5]) * 1000 / 3);
            mv.push(db ? (b.J[t] > 0 ? 20 * Math.log10(b.J[t]) : null) : b.J[t]);
        }
    }
    return {
        type: "scattergl", mode: "markers", x: mx, y: my,
        marker: { size: 4, opacity: 0, color: mv, coloraxis: "coloraxis" },
        hovertemplate: `x: %{x:.4f} mm<br>y: %{y:.4f} mm<br>|J|: %{marker.color:${db ? ".1f} dB(A/m²)" : ".4g} A/m²"}<extra></extra>`,
        showlegend: false,
    };
}

// Hover and color axis of |E| on the mesh: invisible markers at the centroids of at most
// HOVER_MAX triangles.
function fieldMeshHoverTrace(M, db) {
    const n = M.E.length / 3, stride = Math.max(1, Math.ceil(n / HOVER_MAX));
    const mx = [], my = [], mv = [];
    for (let t = 0; t < n; t += stride) {
        const k = 6 * t, v = (M.E[3 * t] + M.E[3 * t + 1] + M.E[3 * t + 2]) / 3;
        mx.push((M.tris[k] + M.tris[k + 2] + M.tris[k + 4]) * 1000 / 3);
        my.push((M.tris[k + 1] + M.tris[k + 3] + M.tris[k + 5]) * 1000 / 3);
        mv.push(db ? (v > 0 ? toDb(v) : null) : v);
    }
    return {
        type: "scattergl", mode: "markers", x: mx, y: my,
        marker: { size: 4, opacity: 0, color: mv, coloraxis: "coloraxis" },
        hovertemplate: `x: %{x:.3f} mm<br>y: %{y:.3f} mm<br>|E|: %{marker.color:${db ? ".1f} dB(V/m)" : ".4g} V/m"}<extra></extra>`,
        showlegend: false,
    };
}

// Contour lines of |E| on the mesh at the log10 levels limits = [start, end, step]
// (efieldContourTrace's), traced through each triangle of the linear field.
function fieldMeshContourTrace(M, limits) {
    const [lo, hi, step] = limits;
    const levels = [];
    for (let L = lo; L <= hi + 1e-9 && levels.length < 500; L += step) levels.push(L);
    const x = [], y = [];
    const T = M.tris, n = T.length / 6;
    const lv = new Float64Array(3);
    for (let t = 0; t < n; t++) {
        for (let a = 0; a < 3; a++) lv[a] = Math.log10(Math.max(M.E[3 * t + a], 1e-3));
        const vmin = Math.min(lv[0], lv[1], lv[2]), vmax = Math.max(lv[0], lv[1], lv[2]);
        if (vmax <= lo || vmin >= hi + 1e-9) continue;
        for (const L of levels) {
            if (L <= vmin || L >= vmax) continue;
            let cnt = 0;
            for (let a = 0; a < 3; a++) {
                const b = (a + 1) % 3, da = lv[a] - L, db = lv[b] - L;
                if ((da < 0) !== (db < 0)) {
                    const u = da / (da - db);
                    x.push((T[6 * t + 2 * a] + u * (T[6 * t + 2 * b] - T[6 * t + 2 * a])) * 1000);
                    y.push((T[6 * t + 2 * a + 1] + u * (T[6 * t + 2 * b + 1] - T[6 * t + 2 * a + 1])) * 1000);
                    cnt++;
                }
            }
            if (cnt === 2) { x.push(null); y.push(null); }
            else if (cnt) { x.length -= cnt; y.length -= cnt; }
        }
    }
    return {
        type: "scattergl", mode: "lines", x, y,
        line: { width: 1, color: "rgba(0, 0, 0, 0.4)" },
        name: "E-field contours", showlegend: false, hoverinfo: "skip",
    };
}

// Triangle blocks { tris, Jv } (vertex values, |J| or |E|) rasterized for the visible
// axis ranges (mm) at w x h pixels: each triangle Gouraud-shaded from its vertex values,
// so a skin layer or a field jump on a slanted or curved face stays sharp and smooth at
// any zoom. Returns a layout image above the shapes of layer 'below' (the dielectric
// fills) and under the traces (mesh overlay, contour lines), or null.
function rasterizeDensity(blocks, xr, yr, w, h, zmin, zmax, db) {
    if (!(w > 0 && h > 0)) return null;
    const canvas = document.createElement('canvas');
    canvas.width = w; canvas.height = h;
    const ctx = canvas.getContext('2d');
    const img = ctx.createImageData(w, h);
    const px = img.data;
    const LUT = Array.from({ length: 256 }, (_, i) => viridisAt(i / 255).match(/\d+/g).map(Number));
    const span = zmax - zmin || 1;
    const sx = w / (xr[1] - xr[0]), sy = h / (yr[1] - yr[0]);
    const val = v => db ? (v > 0 ? 20 * Math.log10(v) : zmin) : v;
    let any = false;
    for (const b of blocks) {
        const T = b.tris, Jv = b.Jv, nT = T.length / 6;
        for (let t = 0; t < nT; t++) {
            const k = 6 * t;
            const ax = (T[k] * 1000 - xr[0]) * sx, ay = (yr[1] - T[k + 1] * 1000) * sy;
            const bx = (T[k + 2] * 1000 - xr[0]) * sx, by = (yr[1] - T[k + 3] * 1000) * sy;
            const cx = (T[k + 4] * 1000 - xr[0]) * sx, cy = (yr[1] - T[k + 5] * 1000) * sy;
            const i0 = Math.max(0, Math.ceil(Math.min(ax, bx, cx) - 0.5)), i1 = Math.min(w - 1, Math.floor(Math.max(ax, bx, cx) - 0.5));
            const j0 = Math.max(0, Math.ceil(Math.min(ay, by, cy) - 0.5)), j1 = Math.min(h - 1, Math.floor(Math.max(ay, by, cy) - 0.5));
            if (i0 > i1 || j0 > j1) continue;
            const det = (bx - ax) * (cy - ay) - (cx - ax) * (by - ay);
            if (!det) continue;
            const va = val(Jv[3 * t]), vb = val(Jv[3 * t + 1]), vc = val(Jv[3 * t + 2]);
            for (let j = j0; j <= j1; j++) {
                const qy = j + 0.5;
                for (let i = i0; i <= i1; i++) {
                    const qx = i + 0.5;
                    const l1 = ((qx - ax) * (cy - ay) - (cx - ax) * (qy - ay)) / det;
                    const l2 = ((bx - ax) * (qy - ay) - (qx - ax) * (by - ay)) / det;
                    const l0 = 1 - l1 - l2;
                    if (l0 < -1e-9 || l1 < -1e-9 || l2 < -1e-9) continue;
                    const u = (l0 * va + l1 * vb + l2 * vc - zmin) / span;
                    const c = LUT[Math.max(0, Math.min(255, Math.round(u * 255)))];
                    const o = 4 * (j * w + i);
                    px[o] = c[0]; px[o + 1] = c[1]; px[o + 2] = c[2]; px[o + 3] = 255;
                    any = true;
                }
            }
        }
    }
    if (!any) return null;
    ctx.putImageData(img, 0, 0);
    return { source: canvas.toDataURL(), xref: 'x', yref: 'y', x: xr[0], y: yr[1],
             sizex: xr[1] - xr[0], sizey: yr[1] - yr[0], sizing: 'stretch', layer: 'below' };
}

// The triangle blocks the current view draws as an image: the shaped conductors' |J|,
// or |E| on the mesh.
function imageBlocks() {
    if (currentView === "density") return (getCurrentJ() || []).filter(b => b.tris);
    const M = currentView.startsWith("efield") ? getFieldMesh() : null;
    return M ? [{ tris: M.tris, Jv: M.E }] : [];
}

// Redraws the image of the triangle blocks for the current axis ranges and plot size.
let densityImageFrame = 0;
function updateDensityImage(container) {
    cancelAnimationFrame(densityImageFrame);
    densityImageFrame = requestAnimationFrame(() => {
        const blocks = imageBlocks();
        const fl = container._fullLayout;
        if (!blocks.length || !fl || !fl.xaxis || !fl._size) return;
        const ratio = Math.min(window.devicePixelRatio || 1, 2);
        const xr = fl.xaxis.range.slice().sort((a, b) => a - b), yr = fl.yaxis.range.slice().sort((a, b) => a - b);
        const im = rasterizeDensity(blocks, xr, yr, Math.round(fl._size.w * ratio), Math.round(fl._size.h * ratio),
            zMin, zMax, getPlotOptions().efieldDb);
        container._densityImageUpdate = true;
        getPlotly().relayout(container, { images: im ? [im] : [] }).finally(() => { container._densityImageUpdate = false; });
    });
}

// A |J| block extended by a hair past each face with copies of its edge samples, so the
// heatmap reaches the conductor outline without an antialiasing seam.
function padBlock(b) {
    const w = Math.abs(b.x[b.x.length - 1] - b.x[0]), h = b.y[b.y.length - 1] - b.y[0];
    const pad = 2e-3 * Math.min(w, h);
    const dir = b.x[b.x.length - 1] >= b.x[0] ? 1 : -1;
    const x = [b.x[0] - dir * pad, ...b.x, b.x[b.x.length - 1] + dir * pad];
    const y = [b.y[0] - pad, ...b.y, b.y[b.y.length - 1] + pad];
    const rows = b.J.map(row => [row[0], ...row, row[row.length - 1]]);
    return { x, y, J: [rows[0], ...rows, rows[rows.length - 1]] };
}

// Autoscale of |J|: the linear max is the level 95 % of the current flows below (the
// corner current is singular; weighting by current, not area, keeps wide ground pours
// with little current from pulling it down), dB spans the peak down 40 dB.
function densityAutoscale(blocks, db) {
    const pts = [];
    for (const b of blocks) {
        if (b.tris) {
            const T = b.tris;
            for (let t = 0; t < b.J.length; t++) {
                const k = 6 * t;
                const area = Math.abs((T[k + 2] - T[k]) * (T[k + 5] - T[k + 1]) - (T[k + 4] - T[k]) * (T[k + 3] - T[k + 1])) / 2;
                if (b.J[t] > 0) pts.push([b.J[t], area]);
            }
            continue;
        }
        const nx = b.x.length, ny = b.y.length;
        for (let j = 0; j < ny; j++) {
            const dy = (b.y[Math.min(j + 1, ny - 1)] - b.y[Math.max(j - 1, 0)]) / 2;
            for (let i = 0; i < nx; i++) {
                const v = b.J[j][i];
                if (v > 0) pts.push([v, dy * (b.x[Math.min(i + 1, nx - 1)] - b.x[Math.max(i - 1, 0)]) / 2]);
            }
        }
    }
    if (!pts.length) return { min: 0, max: 1, peak: 1 };
    pts.sort((a, b) => a[0] - b[0]);
    const peak = pts[pts.length - 1][0];
    if (db) { const max = Math.ceil(20 * Math.log10(peak)); return { min: max - 40, max, peak: max }; }
    let total = 0;
    for (const p of pts) total += p[0] * p[1];
    let acc = 0, q = peak;
    for (const p of pts) { acc += p[0] * p[1]; if (acc >= 0.95 * total) { q = p[0]; break; } }
    return { min: 0, max: q, peak };
}

// Viridis, sampled for the binned segment colors (the colorbar uses the same stops).
const VIRIDIS = ['#440154', '#482878', '#3e4989', '#31688e', '#26828e', '#1f9e89', '#35b779', '#6ece58', '#b5de2b', '#fde725'];
const VIRIDIS_SCALE = VIRIDIS.map((c, i) => [i / (VIRIDIS.length - 1), c]);
function viridisAt(u) {
    const t = Math.min(1, Math.max(0, u)) * (VIRIDIS.length - 1);
    const i = Math.min(VIRIDIS.length - 2, Math.floor(t)), f = t - i;
    const a = VIRIDIS[i], b = VIRIDIS[i + 1];
    const ch = (c, k) => parseInt(c.slice(1 + 2 * k, 3 + 2 * k), 16);
    return `rgb(${[0, 1, 2].map(k => Math.round(ch(a, k) + f * (ch(b, k) - ch(a, k)))).join(',')})`;
}

// Grid samples inside a conductor or outside the mesh (|E| exactly 0) of the solver's
// plot grid, cached on the grid arrays.
function metalMask(solver, z) {
    const X = solver.x, Y = solver.y;
    const c = solver._metalMask;
    if (c && c.x === X && c.y === Y && c.ny === z.length) return c.mask;
    const conductors = solver.conductors || [];
    const mask = z.map((row, j) => Uint8Array.from(row, (v, i) =>
        v === 0 || conductors.some(o => shapeContains(o, X[i], Y[j])) ? 1 : 0));
    solver._metalMask = { x: X, y: Y, ny: z.length, mask };
    return mask;
}

// |E| of the display grid carried two samples into the metal (and past the mesh edge):
// the samples next to the dielectric take the mean of their dielectric neighbours. A
// contour cell across a curved or slanted surface would otherwise pack every contour
// level into a staircase along the grid; the conductor fill hides the metal.
function extendIntoMetal(z, solver) {
    const ny = z.length, nx = ny ? z[0].length : 0;
    if (!nx || !solver.x || solver.x.length !== nx) return z;
    const mask = metalMask(solver, z).map(row => row.slice());
    for (let pass = 0; pass < 2; pass++) {
        const fill = [];
        for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
            if (!mask[j][i]) continue;
            let sum = 0, n = 0;
            for (let dj = -1; dj <= 1; dj++) for (let di = -1; di <= 1; di++) {
                const m = mask[j + dj];
                if (m && m[i + di] === 0) { sum += z[j + dj][i + di]; n++; }
            }
            if (n) fill.push(j, i, sum / n);
        }
        for (let k = 0; k < fill.length; k += 3) { z[fill[k]][fill[k + 1]] = fill[k + 2]; mask[fill[k]][fill[k + 1]] = 0; }
    }
    return z;
}

// Autoscale of the surface current: the linear max is the length-weighted 99th
// percentile (the corner current is singular), dB spans the peak down 40 dB.
function currentAutoscale(K, db) {
    const n = K.K.length;
    const idx = Array.from({ length: n }, (_, i) => i).sort((a, b) => K.K[a] - K.K[b]);
    const len = i => Math.hypot(K.x1[i] - K.x0[i], K.y1[i] - K.y0[i]);
    let total = 0;
    for (let i = 0; i < n; i++) total += len(i);
    let acc = 0, q99 = K.K[idx[n - 1]];
    for (const i of idx) { acc += len(i); if (acc >= 0.99 * total) { q99 = K.K[i]; break; } }
    const peak = K.K[idx[n - 1]];
    if (!db) return { min: 0, max: q99, peak };
    const max = Math.ceil(20 * Math.log10(peak));
    return { min: max - 40, max, peak: max };
}

// Surface current traces: the segments binned by color into one line trace per bin,
// plus invisible midpoint markers carrying the hover text and the colorbar.
function surfaceCurrentTraces(K, zmin, zmax, db, colorbar = { len: 0.8 }) {
    const NB = 48;
    const bins = Array.from({ length: NB }, () => ({ x: [], y: [] }));
    const mx = [], my = [], mv = [];
    const span = zmax - zmin || 1;
    for (let i = 0; i < K.K.length; i++) {
        const v = db ? 20 * Math.log10(Math.max(K.K[i], 1e-30)) : K.K[i];
        const b = bins[Math.min(NB - 1, Math.max(0, Math.floor((v - zmin) / span * NB)))];
        b.x.push(K.x0[i] * 1000, K.x1[i] * 1000, null);
        b.y.push(K.y0[i] * 1000, K.y1[i] * 1000, null);
        mx.push((K.x0[i] + K.x1[i]) * 500); my.push((K.y0[i] + K.y1[i]) * 500); mv.push(v);
    }
    const traces = bins.map((b, k) => ({
        type: "scattergl", mode: "lines", x: b.x, y: b.y,
        line: { width: 5, color: viridisAt((k + 0.5) / NB) },
        hoverinfo: "skip", showlegend: false,
    })).filter(t => t.x.length);
    traces.push({
        type: "scattergl", mode: "markers", x: mx, y: my,
        marker: { size: 6, opacity: 0, color: mv, cmin: zmin, cmax: zmax, colorscale: VIRIDIS_SCALE,
                  showscale: true, colorbar: { title: { text: db ? "dB(A/m)" : "A/m" }, ...colorbar } },
        hovertemplate: `x: %{x:.3f} mm<br>y: %{y:.3f} mm<br>|K|: %{marker.color:${db ? ".1f} dB(A/m)" : ".4g} A/m"}<extra></extra>`,
        showlegend: false,
    });
    return traces;
}

// Whether `view` has data to show: a new geometry or a solve still in progress shows
// the geometry until its fields arrive.
function viewAvailable(view, solver) {
    if (view === "geometry") return true;
    if (!solver.solution_valid) return false;
    if (view.startsWith("potential")) return solver.has_potential !== false && !!solver.V;
    if (view.startsWith("efield")) return !!solver.Ex;
    if (view === "current") return !!getSurfaceK();
    if (view === "density") return !!getCurrentJ();
    return false;
}

function draw(resetZoom = false) {
    const solver = get.solver();
    const Plotly = getPlotly();
    if (!solver || !Plotly) return;
    currentView = viewAvailable(wantedView, solver) ? wantedView : "geometry";

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
        const maxY = displayTop(solver);

        // Calculate intelligent zoom ranges for initial view (only if no current view exists)
        if (!currentXRange || resetZoom) {
            const view = computeGeometryView(solver, maxY);
            if (view) { currentXRange = view.xRange; currentYRange = view.yRange; }
        }

        // The solved region of a custom geometry, drawn as air under everything else.
        if (solver.user_domain) {
            const u = solver.user_domain;
            shapes.push({ type: 'rect', x0: u.x_min * 1000, y0: u.y_min * 1000, x1: u.x_max * 1000, y1: u.y_max * 1000,
                fillcolor: 'rgba(255, 255, 255, 0.8)', line: { color: 'rgba(128, 128, 128, 0.6)', width: 1, dash: 'dot' }, layer: 'below' });
        }
        // Dielectrics (opaque, below the field contours) + conductors above.
        shapes.push(...dielectricFillShapes(solver, maxY));
        shapes.push(...conductorFillShapes(solver, maxY));
        shapes.push(...mirrorImageShapes(solver, maxY));
        shapes.push(...sourceLineHighlightShapes(solver, maxY));

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
                const auto = efieldAutoscale(zData, xMM, yMM);
                extendIntoMetal(zData, solver);
                zMin = 0;
                zMax = auto.max;
                contourFloor = auto.floor;
                actualDataMin = zMin;
                actualDataMax = zMax;
                actualDataPeak = auto.peak;
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
        title = `Electric Potential${modeLabel} (V)${fieldFreqLabel(solver)}`;
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
        // Store actual data range for potential view
        actualDataMin = zMin;
        actualDataMax = zMax;
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
        const db = plotOptions.efieldDb;
        title = `|E| Field Magnitude${modeLabel} (${db ? "dB V/m" : "V/m"})${fieldFreqLabel(solver)}`;
        zTitle = db ? "dB(V/m)" : "V/m";

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
        const auto = efieldAutoscale(zData, xMM, yMM);
        extendIntoMetal(zData, solver);
        contourFloor = auto.floor;
        if (db) {
            ({ min: zMin, max: zMax } = efieldDbRange(auto));
            actualDataPeak = zMax;
        } else {
            zMin = 0;
            zMax = auto.max;
            actualDataPeak = auto.peak;
        }
        actualDataMin = zMin;
        actualDataMax = zMax;
        // Mask the conductor interior (field is 0 inside the PEC) so the heatmap/contour bleed
        // across the boundary is hidden, like the geometry view.
        shapes.push(...conductorFillShapes(solver, yArr[nyDisplay - 1]));
    }

    else if (currentView === "density" && solver.solution_valid && getCurrentJ()) {
        const maxY = displayTop(solver);
        if (!currentXRange || resetZoom) {
            const view = computeGeometryView(solver, maxY);
            if (view) { currentXRange = view.xRange; currentYRange = view.yRange; }
        }
        const db = plotOptions.efieldDb;
        let modeLabel = "";
        if (isDifferentialMode()) modeLabel = getSelectedModeIndex() === 1 ? " (Even Mode)" : " (Odd Mode)";
        const mi = isDifferentialMode() ? getSelectedModeIndex() : 0;
        const ideal = (solver.idealGrounds && solver.idealGrounds[mi] ? ", ideal edge grounds" : "")
            + (bareMetalK(solver, getCurrentJ()) ? ", surface |K| (A/m) on metal without |J|" : "");
        title = `Current Density |J| per 1 A${modeLabel} (${db ? "dB A/m²" : "A/m²"})${fieldFreqLabel(solver)}${ideal}`;
        shapes.push(...dielectricFillShapes(solver, maxY).map(s => ({ ...s, layer: 'below' })));
        shapes.push(...densityConductorShapes(solver, getCurrentJ(), maxY));
        const auto = densityAutoscale(getCurrentJ(), db);
        zMin = auto.min; zMax = auto.max;
        actualDataMin = zMin; actualDataMax = zMax; actualDataPeak = auto.peak;
        ({ xMM, yMM, nx, nyDisplay } = solverGridMM(solver));
    }

    else if (currentView === "current" && solver.solution_valid && getSurfaceK()) {
        const maxY = displayTop(solver);
        if (!currentXRange || resetZoom) {
            const view = computeGeometryView(solver, maxY);
            if (view) { currentXRange = view.xRange; currentYRange = view.yRange; }
        }
        const db = plotOptions.efieldDb;
        let modeLabel = "";
        if (isDifferentialMode()) modeLabel = getSelectedModeIndex() === 1 ? " (Even Mode)" : " (Odd Mode)";
        const mi = isDifferentialMode() ? getSelectedModeIndex() : 0;
        const fromMqs = solver.surfaceKSource && solver.surfaceKSource[mi] === 'mqs';
        const ideal = fromMqs && solver.idealGrounds && solver.idealGrounds[mi] ? ", ideal edge grounds" : "";
        // Both plot the tangential H at the surface. In the perfect-conductor limit it is the
        // current sheet, from the MQS solve only in the skin regime.
        const wg = solver.surfaceKSource && solver.surfaceKSource[mi] === 'waveguide';
        title = wg
            ? `Surface Current |K| = |H| per 1 W (${db ? "dB A/m" : "A/m"})${fieldFreqLabel(solver)}`
            : `Surface Current |K| = |H<sub>t</sub>| per 1 A${modeLabel} (${db ? "dB A/m" : "A/m"})` +
                (fromMqs ? `${fieldFreqLabel(solver)}, MQS${ideal}` : ", perfect-conductor limit");
        const below = s => ({ ...s, layer: 'below' });
        shapes.push(...dielectricFillShapes(solver, maxY).map(below));
        shapes.push(...conductorFillShapes(solver, maxY).map(below));
        const auto = currentAutoscale(getSurfaceK(), db);
        zMin = auto.min; zMax = auto.max;
        actualDataMin = zMin; actualDataMax = zMax; actualDataPeak = auto.peak;
        ({ xMM, yMM, nx, nyDisplay } = solverGridMM(solver));
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

    let colorAxis = null;
    if (currentView === "density" && getCurrentJ()) {
        if (window.getStoredScale) {
            const override = window.getStoredScale(scaleView());
            if (override) { zMin = override.min; zMax = override.max; }
        }
        const db = plotOptions.efieldDb;
        // Walls and ideal edge grounds carry no current inside: their surface |K| is
        // drawn on their faces, on its own colorbar below the |J| one.
        const gK = bareMetalK(solver, getCurrentJ());
        colorAxis = { cmin: zMin, cmax: zMax, colorscale: colorscale,
                      colorbar: gK ? { title: { text: db ? "|J| dB(A/m²)" : "|J| A/m²" }, len: 0.45, y: 1, yanchor: "top" }
                                   : { title: { text: db ? "dB(A/m²)" : "A/m²" }, len: 0.8 } };
        if (gK) {
            const k = currentAutoscale(gK, db);
            traces.push(...surfaceCurrentTraces(gK, k.min, k.max, db,
                { title: { text: db ? "|K| dB(A/m)" : "|K| A/m" }, len: 0.45, y: 0, yanchor: "bottom" }));
        }
        const blocks = getCurrentJ();
        const triBlocks = blocks.filter(b => b.tris);
        if (triBlocks.length) traces.push(densityHoverTrace(triBlocks, db));
        for (const b of blocks.filter(b => !b.tris).map(padBlock)) {
            traces.push({
                type: "heatmap", coloraxis: "coloraxis", zsmooth: "best",
                x: Array.from(b.x, v => v * 1000), y: Array.from(b.y, v => v * 1000),
                z: db ? b.J.map(row => row.map(v => v > 0 ? 20 * Math.log10(v) : null)) : b.J,
                hovertemplate: `x: %{x:.4f} mm<br>y: %{y:.4f} mm<br>|J|: %{z:${db ? ".1f} dB(A/m²)" : ".4g} A/m²"}<extra></extra>`,
            });
        }
    } else if (currentView === "current" && getSurfaceK()) {
        if (window.getStoredScale) {
            const override = window.getStoredScale(scaleView());
            if (override) { zMin = override.min; zMax = override.max; }
        }
        traces.push(...surfaceCurrentTraces(getSurfaceK(), zMin, zMax, plotOptions.efieldDb));
    } else if (currentView === "geometry" && zData.length > 0) {
        const { Ex, Ey } = getFields();

        let eMax = zMax;
        let eMin = contourFloor;

        // Check if there's a user-defined scale override
        if (window.getStoredScale) {
            const override = window.getStoredScale(scaleView());
            if (override) {
                eMin = override.min;
                eMax = override.max;
            }
        }

        const n = plotOptions.contours;

        // Add E-field contours if requested (shared with the |E| field view)
        if (n > 0) {
            const M = getFieldMesh();
            traces.push(M ? fieldMeshContourTrace(M, contourScaledB(eMin, eMax, n))
                : efieldContourTrace(xMM, yMM, zData, contourScaledB(eMin, eMax, n)));
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
        // Field views. Heatmap with optional contour lines.

        // Check if there's a user-defined scale override
        let autoscaled = true;
        if (window.getStoredScale) {
            const override = window.getStoredScale(scaleView());
            if (override) {
                zMin = override.min;
                zMax = override.max;
                autoscaled = false;
            }
        }

        const n = plotOptions.contours;
        const hoverTpl = "x: %{x:.2f} mm<br>y: %{y:.2f} mm<br>value: %{z:.3e}<extra></extra>";

        const fieldMesh = currentView.startsWith("efield") ? getFieldMesh() : null;
        if (fieldMesh) {
            // |E| on the triangles: the color is an image (updateDensityImage), the hover
            // and colorbar ride on invisible markers, the contours are traced per triangle.
            const db = plotOptions.efieldDb;
            colorAxis = { cmin: zMin, cmax: zMax, colorscale: colorscale, colorbar: { title: { text: zTitle }, len: 0.8 } };
            traces.push(fieldMeshHoverTrace(fieldMesh, db));
            if (n > 0) {
                const limits = db ? contourLimitsDb(zMin, zMax, n)
                    : contourScaledB(autoscaled ? contourFloor : Math.max(zMin, 0), zMax, n);
                traces.push(fieldMeshContourTrace(fieldMesh, limits));
            }
        } else if (currentView.startsWith("efield")) {
            const db = plotOptions.efieldDb;
            // |E| heatmap (linear or dB) for the color + colorbar...
            traces.push({
                type: "heatmap",
                zsmooth: "best",
                x: xMM, y: yMM,
                z: db ? zData.map(row => row.map(v => v > 0 ? toDb(v) : null)) : zData,
                zmin: zMin, zmax: zMax,
                colorscale: colorscale,
                colorbar: { title: { text: zTitle }, len: 0.8 },
                hovertemplate: db ? "x: %{x:.2f} mm<br>y: %{y:.2f} mm<br>value: %{z:.1f} dB(V/m)<extra></extra>" : hoverTpl
            });
            // ...overlaid with log-spaced contour lines: the geometry view's in linear, one per
            // n-th of the color range in dB.
            if (n > 0) {
                const limits = db ? contourLimitsDb(zMin, zMax, n)
                    : contourScaledB(autoscaled ? contourFloor : Math.max(zMin, 0), zMax, n);
                traces.push(efieldContourTrace(xMM, yMM, zData, limits));
            }
        } else {
            // Potential (and any other field view): linear heatmap + linear contour lines.
            const contourSettings = { coloring: 'heatmap', showlines: n > 0 };
            if (n > 0) {
                contourSettings.start = zMin;
                contourSettings.end = zMax;
                contourSettings.size = (zMax - zMin) / n;
            }
            traces.push({
                type: n > 0 ? "contour" : "heatmap",
                zsmooth: "best",
                x: xMM, y: yMM, z: zData,
                zmin: zMin, zmax: zMax,
                colorscale: colorscale,
                contours: contourSettings,
                line: { smoothing: 1.3, width: 0.5 },
                colorbar: { title: { text: zTitle }, len: 0.8 },
                hovertemplate: hoverTpl
            });
        }
    }

    // Mesh overlay
    const triMesh = viewTriMesh(solver);
    if (showMesh && solver.solution_valid && triMesh) {
        // Triangular backend: draw triangle edges (deduped) as one batched trace.
        const { nodes, tris, nTris } = triMesh;
        const seen = new Set();
        const ex = [], ey = [];
        const nNodesTri = nodes.length / 2;
        const addEdge = (a, b) => {
            const n0 = a < b ? a : b, n1 = a < b ? b : a;
            const key = n0 * (nNodesTri + 1) + n1;
            if (seen.has(key)) return;
            seen.add(key);
            ex.push(nodes[2 * n0] * 1000, nodes[2 * n1] * 1000, null);
            ey.push(nodes[2 * n0 + 1] * 1000, nodes[2 * n1 + 1] * 1000, null);
        };
        for (let t = 0; t < nTris; t++) {
            const v0 = tris[3 * t], v1 = tris[3 * t + 1], v2 = tris[3 * t + 2];
            addEdge(v0, v1); addEdge(v1, v2); addEdge(v2, v0);
        }
        traces.push({
            type: "scattergl", x: ex, y: ey, mode: "lines",
            line: { width: 0.3, color: "rgba(0,0,0,0.5)" },
            showlegend: false, hoverinfo: "skip"
        });
    } else if (showMesh && solver.solution_valid) {
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

    // An |E| image below the traces sits under the grid lines, which the heatmap covered.
    const gridOff = currentView.startsWith("efield") && !!getFieldMesh();

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
            zerolinecolor: '#555',
            showgrid: !gridOff, zeroline: !gridOff
        },
        yaxis: {
            title: { text: "Height (mm)", font: { color: '#aaa' } },
            range: currentYRange,  // Preserve zoom/pan
            color: '#aaa',
            gridcolor: '#444',
            zerolinecolor: '#555',
            showgrid: !gridOff, zeroline: !gridOff
        },
        margin: { l: 70, r: 90, t: 50, b: 60 },
        showlegend: false,
        hovermode: "closest",
        dragmode: "pan",
        paper_bgcolor: '#2a2a2a',
        plot_bgcolor: '#1a1a1a',
        font: { color: '#fff' },
        shapes: shapes,  // Add vector shapes for geometry
        ...(colorAxis ? { coloraxis: colorAxis } : {}),

        updatemenus: (() => {
            const menus = [];

            // View selector (Geometry/Potential/E-field)
            const viewButtons = [{ label: "Geometry", method: "skip", args: [] }];
            if (solver.solution_valid) {
                // A source-free medium (rectangular waveguide) has no static potential to
                // show, its field is the mode field, so the Potential button is omitted
                // rather than left to render a blank heatmap.
                if (solver.has_potential !== false) {
                    viewButtons.push({ label: "Potential", method: "skip", args: [] });
                }
                viewButtons.push({ label: "|E| Field", method: "skip", args: [] });
                if (solver.surfaceK && solver.surfaceK.some(Boolean)) {
                    viewButtons.push({ label: "|K| Current", method: "skip", args: [] });
                }
                if (solver.currentJ && solver.currentJ.some(Boolean)) {
                    viewButtons.push({ label: "|J| Density", method: "skip", args: [] });
                }
            }
            // Both the highlighted button and the click handler key off the LABEL, never a
            // fixed index, with Potential absent, "|E| Field" is at index 1, not 2.
            // Prefix match so the differential "_odd"/"_even" view variants land on their
            // own button rather than falling through to the first one.
            const activeLabel = currentView.startsWith("geometry") ? "Geometry"
                : currentView.startsWith("potential") ? "Potential"
                : currentView === "current" ? "|K| Current"
                : currentView === "density" ? "|J| Density" : "|E| Field";
            menus.push({
                x: 0.01,
                y: 1.15,
                showactive: true,
                active: Math.max(0, viewButtons.findIndex(b => b.label === activeLabel)),
                bgcolor: '#2a2a2a',
                bordercolor: '#444',
                font: { color: '#aaa' },
                buttons: viewButtons
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
        modeBarButtonsToRemove: ["select2d", "lasso2d"],
        modeBarButtonsToAdd: [
            {
                name: "Toggle Mesh",
                icon: GRID_ICON,
                click: () => {
                    showMesh = !showMesh;
                    draw();
                }
            },
            {
                name: "Scale Range",
                icon: Plotly.Icons.autoscale,
                click: () => window.toggleScaleDialog && window.toggleScaleDialog()
            },
            {
                name: "About this plot",
                icon: Plotly.Icons.question,
                click: () => window.showHelpModal && window.showHelpModal(plotHelpTopic())
            }
        ]
    };

    Plotly.react(container, traces, layout, config);
    if (currentView === "density" || currentView.startsWith("efield")) updateDensityImage(container);

    if (!container._viewListenerBound) {
        container.on('plotly_buttonclicked', (event) => {
            // Determine which menu was clicked based on x position
            // First menu (x=0.01): View selector (Geometry/Potential/E-field)
            // Second menu (x=0.25): Mode selector (Odd/Even) - only for differential

            if (event.menu.x < 0.2) {
                // View selector clicked. Key off the LABEL, not the index: the Potential
                // button is absent for a source-free medium (see the button list above),
                // so index 1 is not always "potential".
                const btn = event.menu.buttons[event.menu.active];
                const label = btn && btn.label;
                setCurrentView(label === "Geometry" ? "geometry"
                    : label === "Potential" ? "potential"
                    : label === "|K| Current" ? "current"
                    : label === "|J| Density" ? "density" : "efield");
            } else {
                // Mode selector clicked (differential lines only)
                const plotModeEl = document.getElementById('plot-mode');
                if (plotModeEl) {
                    plotModeEl.value = event.menu.active === 0 ? 'odd' : 'even';
                }
                // Trigger view change notification for mode switch
                if (window.onViewChanged) {
                    window.onViewChanged(currentView);
                }
            }
            draw();
        });
        container._viewListenerBound = true;
    }

    // Listen for autoscale events to reset color scale
    if (!container._autoscaleListenerBound) {
        let ignoreNextAutoscale = false;

        // Track double-clicks to distinguish from autoscale button
        container.on('plotly_doubleclick', () => {
            ignoreNextAutoscale = true;
            // Clear the flag after a short delay in case the autoscale event doesn't fire
            setTimeout(() => {
                ignoreNextAutoscale = false;
            }, 200);
        });

        // Handle autoscale button click
        container.on('plotly_relayout', (eventData) => {
            // A zoom or pan redraws the |J| image of the shaped conductors for the new view.
            if ((currentView === "density" || currentView.startsWith("efield")) && !container._densityImageUpdate)
                updateDensityImage(container);
            // Check if this is an autoscale event (both axes autoscaling)
            if (eventData && eventData['xaxis.autorange'] === true && eventData['yaxis.autorange'] === true) {
                // Only reset color scale if this is from the autoscale button, not double-click
                if (!ignoreNextAutoscale && window.resetColorScale) {
                    window.resetColorScale();
                }
                ignoreNextAutoscale = false;
            }
        });

        container._autoscaleListenerBound = true;
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

/**
 * Extract a single value from a mode result for a given selector and scaling.
 * Used by both frequency sweep results and parameter sweep plots.
 */
function extractModeValue(mode, selector, scale) {
    switch (selector) {
        case 're_z0':   return scale * mode.Zc.re;
        case 'im_z0':   return scale * mode.Zc.im;
        case 'eps_eff': return mode.eps_eff;
        case 'loss':    return mode.alpha_total;
        default:        return scale * mode.RLGC[selector];
    }
}

function buildResultsTraces(sweepResults, selector, useDiffMode) {
    const resultsAreDifferential = sweepResults[0].result.modes.length === 2;
    const freqs = sweepResults.map(r => r.freq / 1e9);
    const plotMode = freqs.length === 1 ? 'markers' : 'lines+markers';
    const traces = [];

    // Mode labels
    const mode0 = useDiffMode ? 'Differential' : 'Odd';
    const mode1 = useDiffMode ? 'Common' : 'Even';

    if (selector === 'loss') {
        if (resultsAreDifferential) {
            const suffix0 = useDiffMode ? 'diff' : 'odd';
            const suffix1 = useDiffMode ? 'common' : 'even';
            // Mode 0 losses (solid lines)
            traces.push({
                x: freqs,
                y: sweepResults.map(r => r.result.modes[0].alpha_c),
                name: `Conductor (${suffix0})`, type: 'scatter', mode: plotMode
            });
            traces.push({
                x: freqs,
                y: sweepResults.map(r => r.result.modes[0].alpha_d),
                name: `Dielectric (${suffix0})`, type: 'scatter', mode: plotMode
            });
            traces.push({
                x: freqs,
                y: sweepResults.map(r => r.result.modes[0].alpha_total),
                name: `Total (${suffix0})`, type: 'scatter', mode: plotMode,
                line: { width: 2 }
            });
            // Mode 1 losses (dashed lines)
            traces.push({
                x: freqs,
                y: sweepResults.map(r => r.result.modes[1].alpha_c),
                name: `Conductor (${suffix1})`, type: 'scatter', mode: plotMode,
                line: { dash: 'dash' }
            });
            traces.push({
                x: freqs,
                y: sweepResults.map(r => r.result.modes[1].alpha_d),
                name: `Dielectric (${suffix1})`, type: 'scatter', mode: plotMode,
                line: { dash: 'dash' }
            });
            traces.push({
                x: freqs,
                y: sweepResults.map(r => r.result.modes[1].alpha_total),
                name: `Total (${suffix1})`, type: 'scatter', mode: plotMode,
                line: { width: 2, dash: 'dash' }
            });
        } else {
            traces.push({
                x: freqs,
                y: sweepResults.map(r => r.result.modes[0].alpha_c),
                name: 'Conductor', type: 'scatter', mode: plotMode
            });
            traces.push({
                x: freqs,
                y: sweepResults.map(r => r.result.modes[0].alpha_d),
                name: 'Dielectric', type: 'scatter', mode: plotMode
            });
            traces.push({
                x: freqs,
                y: sweepResults.map(r => r.result.modes[0].alpha_total),
                name: 'Total', type: 'scatter', mode: plotMode,
                line: { width: 2 }
            });
        }
    } else {
        // Z0, eps_eff, RLGC parameters
        const scale0 = useDiffMode ? 2 : 1;
        const scale1 = useDiffMode ? 0.5 : 1;
        if (resultsAreDifferential) {
            traces.push({
                x: freqs,
                y: sweepResults.map(r => extractModeValue(r.result.modes[0], selector, scale0)),
                name: `${mode0} mode`, type: 'scatter', mode: plotMode
            });
            traces.push({
                x: freqs,
                y: sweepResults.map(r => extractModeValue(r.result.modes[1], selector, scale1)),
                name: `${mode1} mode`, type: 'scatter', mode: plotMode
            });
        } else {
            traces.push({
                x: freqs,
                y: sweepResults.map(r => extractModeValue(r.result.modes[0], selector, 1)),
                name: getYAxisLabel(selector), type: 'scatter', mode: plotMode
            });
        }
    }

    return traces;
}

function drawResultsPlot() {
    const frequencySweepResults = get.frequencySweepResults();
    const Plotly = getPlotly();
    if (!frequencySweepResults || frequencySweepResults.length === 0 || !Plotly) return;

    const selector = document.getElementById('results-plot-selector').value;
    const resultsAreDifferential = frequencySweepResults[0].result.modes.length === 2;
    const useDiffMode = document.getElementById('results-diff').checked && resultsAreDifferential;

    const activeTraces = buildResultsTraces(frequencySweepResults, selector, useDiffMode);

    // Assign explicit colors and legend groups so frozen traces don't shift the color cycle
    for (let i = 0; i < activeTraces.length; i++) {
        const color = PLOTLY_COLORS[i % PLOTLY_COLORS.length];
        activeTraces[i].line = { ...activeTraces[i].line, color };
        activeTraces[i].marker = { color };
        activeTraces[i].legendgroup = `group${i}`;
    }

    const allTraces = [];

    if (frozenResultsData) {
        const frozenDiff = frozenResultsData[0].result.modes.length === 2;
        const frozenUseDiff = document.getElementById('results-diff').checked && frozenDiff;
        const frozen = buildResultsTraces(frozenResultsData, selector, frozenUseDiff);
        for (let i = 0; i < frozen.length; i++) {
            const color = PLOTLY_COLORS[i % PLOTLY_COLORS.length];
            frozen[i].line = { color };
            frozen[i].opacity = 0.35;
            frozen[i].showlegend = false;
            frozen[i].hoverinfo = 'skip';
            frozen[i].mode = 'lines';
            frozen[i].legendgroup = `group${i}`;
        }
        allTraces.push(...frozen);
    }

    allTraces.push(...activeTraces);

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

    Plotly.newPlot('results-plot', allTraces, layout, { responsive: true });
}

function buildSParamTraces(sweepResults, length, Z_ref, plotMode, useMixedMode) {
    // A self-referenced medium (waveguide) drops its below-cutoff points and normalizes
    // each remaining one to its own modal impedance; every other medium passes through.
    sweepResults = usableSweepPoints(sweepResults);
    if (!sweepResults.length) return [];
    const resultsAreDifferential = sweepResults[0].result.modes.length === 2;
    const freqs = sweepResults.map(r => r.freq / 1e9);
    const lineMode = freqs.length === 1 ? 'markers' : 'lines+markers';
    const traces = [];

    const sParamToPhase = (complexVal) => complexVal.arg() * 180 / Math.PI;

    if (!resultsAreDifferential) {
        const S11_data = [];
        const S21_data = [];

        for (const { freq, result } of sweepResults) {
            const sp = sparamsForPoint(freq, result, length, Z_ref);
            if (plotMode === 'magnitude') {
                S11_data.push(sParamTodB(sp.S11));
                S21_data.push(sParamTodB(sp.S21));
            } else {
                S11_data.push(sParamToPhase(sp.S11));
                S21_data.push(sParamToPhase(sp.S21));
            }
        }

        const label = plotMode === 'magnitude' ? '(dB)' : '(deg)';
        traces.push({ x: freqs, y: S11_data, name: `S11 ${label}`, type: 'scatter', mode: lineMode });
        traces.push({ x: freqs, y: S21_data, name: `S21 ${label}`, type: 'scatter', mode: lineMode });
    } else if (useMixedMode) {
        // An asymmetric pair converts between differential and common mode: SDC/SCD are non-zero
        // (zero for a symmetric line). Plot those terms for an asymmetric pair so the mixed-mode
        // signature of the asymmetry is visible, not just the pure SDD/SCC responses.
        const isAsymmetric = sweepResults.some(({ result }) => pairIsAsymmetric(result));
        const SDD11_data = [], SDD21_data = [], SCC11_data = [], SCC21_data = [];
        const SDC11_data = [], SCD11_data = [], SCD21_data = [];
        const conv = plotMode === 'magnitude' ? sParamTodB : sParamToPhase;

        for (const { freq, result } of sweepResults) {
            const oddMode = result.modes.find(m => m.mode === 'odd');
            const evenMode = result.modes.find(m => m.mode === 'even');
            const sp = computeSParamsDiffAuto(freq, oddMode.RLGC, evenMode.RLGC, result.physMatrix, length, Z_ref);

            SDD11_data.push(conv(sp.SDD11));
            SDD21_data.push(conv(sp.SDD21));
            SCC11_data.push(conv(sp.SCC11));
            SCC21_data.push(conv(sp.SCC21));
            if (isAsymmetric) {
                SDC11_data.push(conv(sp.SDC11));   // common→differential conversion (reflection)
                SCD11_data.push(conv(sp.SCD11));   // differential→common conversion (reflection)
                SCD21_data.push(conv(sp.SCD21));   // differential→common conversion (transmission)
            }
        }

        const label = plotMode === 'magnitude' ? '(dB)' : '(deg)';
        traces.push({ x: freqs, y: SDD11_data, name: `SDD11 ${label}`, type: 'scatter', mode: lineMode });
        traces.push({ x: freqs, y: SDD21_data, name: `SDD21 ${label}`, type: 'scatter', mode: lineMode });
        traces.push({ x: freqs, y: SCC11_data, name: `SCC11 ${label}`, type: 'scatter', mode: lineMode, line: { dash: 'dash' } });
        traces.push({ x: freqs, y: SCC21_data, name: `SCC21 ${label}`, type: 'scatter', mode: lineMode, line: { dash: 'dash' } });
        if (isAsymmetric) {
            traces.push({ x: freqs, y: SDC11_data, name: `SDC11 ${label}`, type: 'scatter', mode: lineMode, line: { dash: 'dot' } });
            traces.push({ x: freqs, y: SCD11_data, name: `SCD11 ${label}`, type: 'scatter', mode: lineMode, line: { dash: 'dot' } });
            traces.push({ x: freqs, y: SCD21_data, name: `SCD21 ${label}`, type: 'scatter', mode: lineMode, line: { dash: 'dot' } });
        }
    } else {
        // An asymmetric coupled pair has a non-degenerate second column:
        // S22≠S11, S32≠S41, S42≠S31. Plot those extra terms so the asymmetry is visible.
        const isAsymmetric = sweepResults.some(({ result }) => pairIsAsymmetric(result));

        const S11_data = [], S21_data = [], S31_data = [], S41_data = [];
        const S22_data = [], S32_data = [], S42_data = [];
        const conv = plotMode === 'magnitude' ? sParamTodB : sParamToPhase;

        for (const { freq, result } of sweepResults) {
            const oddMode = result.modes.find(m => m.mode === 'odd');
            const evenMode = result.modes.find(m => m.mode === 'even');
            const sp = computeSParamsDiffAuto(freq, oddMode.RLGC, evenMode.RLGC, result.physMatrix, length, Z_ref);

            S11_data.push(conv(sp.S[0][0]));
            S21_data.push(conv(sp.S[1][0]));
            S31_data.push(conv(sp.S[2][0]));
            S41_data.push(conv(sp.S[3][0]));
            if (isAsymmetric) {
                S22_data.push(conv(sp.S[1][1]));
                S32_data.push(conv(sp.S[2][1]));
                S42_data.push(conv(sp.S[3][1]));
            }
        }

        const label = plotMode === 'magnitude' ? '(dB)' : '(deg)';
        traces.push({ x: freqs, y: S11_data, name: `S11 ${label}`, type: 'scatter', mode: lineMode });
        traces.push({ x: freqs, y: S21_data, name: `S21 ${label}`, type: 'scatter', mode: lineMode });
        traces.push({ x: freqs, y: S31_data, name: `S31 ${label}`, type: 'scatter', mode: lineMode });
        traces.push({ x: freqs, y: S41_data, name: `S41 ${label}`, type: 'scatter', mode: lineMode, line: { dash: 'dash' } });
        if (isAsymmetric) {
            traces.push({ x: freqs, y: S22_data, name: `S22 ${label}`, type: 'scatter', mode: lineMode, line: { dash: 'dot' } });
            traces.push({ x: freqs, y: S32_data, name: `S32 ${label}`, type: 'scatter', mode: lineMode, line: { dash: 'dot' } });
            traces.push({ x: freqs, y: S42_data, name: `S42 ${label}`, type: 'scatter', mode: lineMode, line: { dash: 'dot' } });
        }
    }

    return traces;
}

// A pair whose two lines differ: by geometry (physMatrix), or by the metal or finish of
// its traces on a mirror-symmetric geometry (dR / dL on the mode RLGC).
function pairIsAsymmetric(result) {
    if (result.physMatrix) return true;
    const odd = result.modes && result.modes.find(m => m.mode === 'odd');
    return !!(odd && odd.RLGC && (odd.RLGC.dR || odd.RLGC.dL));
}

function drawSParamPlot() {
    const frequencySweepResults = get.frequencySweepResults();
    const Plotly = getPlotly();
    if (!frequencySweepResults || frequencySweepResults.length === 0 || !Plotly) return;

    const length = get.inputValue('sparam-length');
    const Z_ref = parseFloat(document.getElementById('sparam-z-ref').value);
    const useMixedMode = document.getElementById('sparam-diff').checked;

    // A self-referenced medium ignores the reference-impedance box entirely (the UI hides
    // it and shows "Z0" instead), so it must not gate on parsing that field.
    const selfRef = isSelfReferenced(frequencySweepResults);
    if (isNaN(length) || length <= 0 || (!selfRef && (isNaN(Z_ref) || Z_ref <= 0))) {
        return;
    }

    const plotMode = document.getElementById('sparam-plot-mode').value;

    const activeTraces = buildSParamTraces(frequencySweepResults, length, Z_ref, plotMode, useMixedMode);

    // Assign explicit colors and legend groups so frozen traces don't shift the color cycle
    for (let i = 0; i < activeTraces.length; i++) {
        const color = PLOTLY_COLORS[i % PLOTLY_COLORS.length];
        activeTraces[i].line = { ...activeTraces[i].line, color };
        activeTraces[i].marker = { color };
        activeTraces[i].legendgroup = `group${i}`;
    }

    const allTraces = [];

    if (frozenSParamData) {
        const frozen = buildSParamTraces(
            frozenSParamData.results, frozenSParamData.length,
            frozenSParamData.zRef, plotMode, useMixedMode
        );
        for (let i = 0; i < frozen.length; i++) {
            const color = PLOTLY_COLORS[i % PLOTLY_COLORS.length];
            frozen[i].line = { color };
            frozen[i].opacity = 0.35;
            frozen[i].showlegend = false;
            frozen[i].hoverinfo = 'skip';
            frozen[i].mode = 'lines';
            frozen[i].legendgroup = `group${i}`;
        }
        allTraces.push(...frozen);
    }

    allTraces.push(...activeTraces);

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

    Plotly.newPlot('sparam-plot', allTraces, layout, { responsive: true });
}

function drawParameterSweepPlot(sweepData, xLabel, ySelector, useDiffMode) {
    const Plotly = getPlotly();
    if (!sweepData || sweepData.length === 0 || !Plotly) return;
    const xVals = sweepData.map(d => d.paramValue);
    const isDiff = sweepData[0].result.modes.length === 2;

    const name0 = !isDiff ? getYAxisLabel(ySelector) : (useDiffMode ? 'Differential' : 'Odd');
    const name1 = useDiffMode ? 'Common' : 'Even';
    const scale0 = isDiff && useDiffMode ? 2 : 1;
    const scale1 = isDiff && useDiffMode ? 0.5 : 1;

    const yVals0 = sweepData.map(d => extractModeValue(d.result.modes[0], ySelector, scale0));
    const yVals1 = isDiff ? sweepData.map(d => extractModeValue(d.result.modes[1], ySelector, scale1)) : null;

    const traces = [];

    if (xVals.length >= 2) {
        // Dense interpolated traces for hover-anywhere-on-line capability
        const interpPts = 500;
        const addInterpTrace = (xArr, yArr, name, color) => {
            const xInterp = [];
            const yInterp = [];
            for (let i = 0; i < xArr.length - 1; i++) {
                const nSeg = Math.max(2, Math.round(interpPts / (xArr.length - 1)));
                for (let j = 0; j < nSeg; j++) {
                    const t = j / nSeg;
                    xInterp.push(xArr[i] + t * (xArr[i + 1] - xArr[i]));
                    yInterp.push(yArr[i] + t * (yArr[i + 1] - yArr[i]));
                }
            }
            // Add final point
            xInterp.push(xArr[xArr.length - 1]);
            yInterp.push(yArr[yArr.length - 1]);

            // Interpolated line trace (hoverable, no visible markers)
            traces.push({
                x: xInterp, y: yInterp,
                name, type: 'scatter', mode: 'lines',
                line: { color, width: 2 },
                hoverinfo: 'x+y+name',
                showlegend: false
            });
            // Markers at actual computed points
            traces.push({
                x: xArr, y: yArr,
                name, type: 'scatter', mode: 'markers',
                marker: { color, size: 7 },
                hoverinfo: 'x+y+name',
                legendgroup: name,
                showlegend: true
            });
        };

        addInterpTrace(xVals, yVals0, name0, PLOTLY_COLORS[0]);
        if (isDiff) addInterpTrace(xVals, yVals1, name1, PLOTLY_COLORS[1]);
    } else {
        // Single point - just markers
        traces.push({ x: xVals, y: yVals0,
            name: name0, type: 'scatter', mode: 'markers',
            marker: { color: PLOTLY_COLORS[0] } });
        if (isDiff) {
            traces.push({ x: xVals, y: yVals1,
                name: name1, type: 'scatter', mode: 'markers',
                marker: { color: PLOTLY_COLORS[1] } });
        }
    }

    const layout = {
        xaxis: { title: { text: xLabel, font: { color: '#aaa' } }, color: '#aaa', gridcolor: '#444', zerolinecolor: '#555' },
        yaxis: { title: { text: getYAxisLabel(ySelector), font: { color: '#aaa' } }, color: '#aaa', gridcolor: '#444', zerolinecolor: '#555' },
        margin: { l: 80, r: 40, t: 40, b: 60 },
        hovermode: 'closest',
        showlegend: isDiff,
        legend: { x: 0.02, y: 0.98, font: { color: '#fff' } },
        paper_bgcolor: '#2a2a2a', plot_bgcolor: '#1a1a1a', font: { color: '#fff' }
    };
    Plotly.newPlot('sweep-plot', traces, layout, { responsive: true });
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
        contours: contoursVal === '' ? 0 : parseInt(contoursVal) || 0,
        efieldDb: !!document.getElementById('plot-efield-db')?.checked
    };
}

// Function to set current view
function setCurrentView(view) {
    currentView = view;
    wantedView = view;
    // Notify app.js that view changed so it can restore the appropriate scale
    if (window.onViewChanged) {
        window.onViewChanged(view);
    }
}

// Unified freeze/unfreeze for both Results and S-Parameters tabs
function freeze() {
    const data = get.frequencySweepResults();
    if (data && data.length > 0) {
        frozenResultsData = JSON.parse(JSON.stringify(data));
        frozenSParamData = {
            results: JSON.parse(JSON.stringify(data)),
            length: get.inputValue('sparam-length'),
            zRef: parseFloat(document.getElementById('sparam-z-ref').value)
        };
    }
}
function unfreeze() {
    frozenResultsData = null;
    frozenSParamData = null;
}
function isFrozen() { return frozenResultsData !== null; }

export { draw, drawResultsPlot, drawSParamPlot, drawParameterSweepPlot, setGlobals, setCurrentView, getScaleRange, setScaleRange, getActualDataRange,
    freeze, unfreeze, isFrozen, conductorFillShapes, dielectricFillShapes, computeGeometryView, displayTop,
    rasterizeDensity };
