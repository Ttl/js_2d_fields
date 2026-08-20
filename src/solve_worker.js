// Solve worker.
//
// Every solve needs multiple synchronous WASM calls, and a single one of them
// is not interruptible. One call can take tens of seconds freezing the UI.
// Firefox's "this tab is slowing down the browser" dialog (which fires around
// 5 s) can be triggered during a solve. The only fix is to run the solve on
// a different thread, which is this file.
//
// The main thread keeps a solver instance too, but purely as a geometry/plot view model
// (construction is cheap: no meshing, no linear algebra). After a solve the worker ships
// back the field arrays and app_solver grafts them onto that instance, so plot.js keeps
// reading `solver.x/y/V/Ex/Ey/triMesh` exactly as before.
//
// Protocol
// main to worker:
//   { id, type: 'simulate'  , params, frequencies, opts }
//   { id, type: 'modes'     , params, freq, nev, refineOpts }
//   { id, type: 'modeField' , idx }          (after a 'modes' job, same worker)
//   { id, type: 'paramSweep', points, freqHz, opts }
//   {     type: 'stop' }                     (no id: cancels whatever is running)
// worker to main:
//   { id, type: 'log'     , msg }
//   { id, type: 'progress', frac, text }
//   { id, type: 'partial' , ... }            (live plot data mid-sweep)
//   { id, type: 'done'    , ... }
//   { id, type: 'error'   , message }
//
// Cancellation: `stop` sets a flag that the existing shouldStop callbacks poll. The
// worker only sees the message between WASM calls.

import { buildSolverFromParams } from './solver_factory.js';
import { InterpolatingSweep } from './interpolating_sweep.js';

let stopRequested = false;
let currentId = null;
let modesSolver = null;   // retained between 'modes' and its follow-up 'modeField' calls

const post = (msg, transfer) => self.postMessage(msg, transfer || []);
const log = (msg) => post({ id: currentId, type: 'log', msg });
const progress = (frac, text) => post({ id: currentId, type: 'progress', frac, text });
const shouldStop = () => stopRequested;

// Serialization

// Per-mode field arrays are the biggest thing a solve produces (V, Ex, Ey, and
// the vacuum-drive Ex0/Ey0/V0 the loss integrand caches). They are needed
// inside the worker, computeAtFrequency reuses them across the sweep, but the
// main thread only ever plots one solve. Sending them with all ~30 sweep points
// would move tens of MB per solve for nothing, so every result crossing the
// boundary is stripped, and the fields travel once, separately, via
// fieldPayload().
const HEAVY = ['V', 'Ex', 'Ey', 'V0', 'Ex0', 'Ey0'];

function stripResult(result) {
    if (!result) return result;
    const out = { ...result };
    if (Array.isArray(result.modes)) {
        out.modes = result.modes.map(m => {
            const c = { ...m };
            for (const k of HEAVY) delete c[k];
            // Zc is a Complex instance. Structured clone drops the prototype, so it
            // arrives as a bare {re, im}. Send it as such deliberately and let the main
            // thread revive it.
            if (c.Zc) c.Zc = { re: c.Zc.re, im: c.Zc.im };
            return c;
        });
    }
    return out;
}

const stripSweep = (rows) => rows.map(({ freq, result }) => ({ freq, result: stripResult(result) }));

// Exactly the set of properties both backends graft onto the solver for plotting
// (FieldSolver2D.solve_adaptive and TriBackend.solveAt agree on this list).
// getPlotFields mirrors half-domain (sym_half) FDM solves onto the full domain,
// for everything else it returns the same properties unchanged.
function fieldPayload(solver) {
    if (!solver || !solver.solution_valid) return null;
    if (solver.getPlotFields) return solver.getPlotFields();
    return {
        x: solver.x, y: solver.y,
        V: solver.V, Ex: solver.Ex, Ey: solver.Ey,
        triMesh: solver.triMesh || null,
    };
}

// The accuracy notes the main thread cannot derive from its own (never-solved) view
// model. Reported once, from the 'meshDone' message, at the point in the log where the
// main-thread version used to print them.
const solveMeta = (solver) => ({
    certification: solver.certification || null,
    meshQuality: solver.meshQuality || null,
});

function makeSolver(params) {
    const s = buildSolverFromParams(params, log);
    if (!s) throw new Error('Solver initialization failed due to invalid parameters.');
    return s;
}

// Simulate

// Progress-bar work model, moved here verbatim from runSimulation so the worker can emit
// a ready-made fraction: mesh refinement is a handful of passes on a coarse mesh and is
// much cheaper than the sweep that follows, so it owns only the small leading slice.
const EST_MESH_PASSES = 4;
const EST_SWEEP_POINTS = 30;
const MESH_PASS_COST = 0.5;
const MESH_FRACTION = (EST_MESH_PASSES * MESH_PASS_COST) /
    (EST_MESH_PASSES * MESH_PASS_COST + EST_SWEEP_POINTS);

async function jobSimulate({ params, frequencies, opts }) {
    const solver = makeSolver(params);
    const p = params;

    if (p.sigma < 1e4) throw new Error('Signal line conductivity is too low to be considered a conductor.');
    if (p.rq < 0) throw new Error('Surface roughness cannot be negative.');
    if (p.use_plating) {
        if (p.plating_sigma < 1e4) throw new Error('Plating conductivity is too low to be considered a conductor.');
        if (p.plating_t < 0) throw new Error('Plating thickness must be non-negative.');
        if (p.plating_rq < 0) throw new Error('Plating roughness cannot be negative.');
    }

    let sweepResults = [];
    const maxFreq = Math.max(...frequencies);
    solver.freq = maxFreq;

    log('Calculating mesh...');
    solver.ensure_mesh();
    log(solver.x ? ('Mesh generated: ' + solver.x.length + 'x' + solver.y.length)
                 : 'Triangular mesh will be generated by the solver...');
    log(`Running adaptive analysis (max ${p.max_iters} iterations, max ${p.max_nodes}k nodes, ` +
        `tolerance ${(100 * p.tolerance).toFixed(2)}%)...`);

    const results = await solver.solve_adaptive({
        max_iters: p.max_iters,
        energy_tol: p.tolerance,
        param_tol: 0.05,
        max_nodes: p.max_nodes * 1000,
        min_converged_passes: p.min_converged_passes,
        certify: !!p.estimate_error,
        onProgress: (info) => {
            const frac = Math.min(info.iteration / EST_MESH_PASSES, 1) * MESH_FRACTION;
            const meshStr = info.n_tris != null ? `Tris=${info.n_tris}` : `Grid=${info.nodes_x}x${info.nodes_y}`;
            progress(frac, `Pass ${info.iteration}/${p.max_iters} · error ${info.energy_error.toExponential(2)}`);
            log(`Pass ${info.iteration}: Energy error=${info.energy_error.toExponential(3)}, ` +
                `Param error=${info.param_error.toExponential(3)}, ${meshStr}`);
        },
        shouldStop,
    });

    post({ id: currentId, type: 'meshDone', meta: solveMeta(solver),
           warnings: (results && results.warnings) || solver.modeWarnings || null });

    if (stopRequested) return { stopped: true, sweepResults: [], fields: null };

    const cachedResults = results;
    if (solver.use_causal_materials) {
        sweepResults.push({ freq: maxFreq, result: await solver.computeAtFrequency(maxFreq, cachedResults) });
    } else {
        sweepResults.push({ freq: maxFreq, result: results });
    }

    // First plottable state: the converged mesh solve. Ship the fields now so the
    // geometry tab paints its E-field overlay while the sweep is still running.
    post({ id: currentId, type: 'partial', fields: fieldPayload(solver),
           sweepResults: stripSweep(sweepResults) });

    const fMax = maxFreq;
    const nonZeroFreqs = frequencies.filter(f => f > 0);
    const fMinNonZero = nonZeroFreqs.length > 0 ? Math.min(...nonZeroFreqs) : 0;
    const useInterpolation = opts.useInterpolation
        && solver.allow_interp_sweep !== false
        && nonZeroFreqs.length > 1
        && fMax > fMinNonZero;

    let sweepWarnings = null;

    if (useInterpolation) {
        const tolerance = opts.interpTolerance;
        const hasDC = frequencies.includes(0);
        if (hasDC) {
            sweepResults.push({ freq: 0, result: await solver.computeAtFrequency(0, cachedResults) });
        }
        log(`Interpolating sweep (tolerance ${(tolerance * 100)}%)...`);

        const sweep = new InterpolatingSweep(solver, cachedResults, { tolerance });

        const SWEEP_START = MESH_FRACTION, EST_TOTAL_POINTS = EST_SWEEP_POINTS, POINT_CAP = 0.9;
        const setSweep = (frac, text) =>
            progress(SWEEP_START + (1 - SWEEP_START) * Math.max(0, Math.min(1, frac)), text);
        let refineErr0 = null, lastErrFrac = 0, lastFinalErr = null;

        const nSamples = await sweep.run(fMinNonZero, fMax, {
            onProgress: (info) => {
                const pointFloor = Math.min(POINT_CAP, info.totalSamples / EST_TOTAL_POINTS);
                if (info.phase === 'initial') {
                    setSweep(pointFloor, `Solving initial points ${info.pointsComputed}/${info.initialPoints}`);
                } else {
                    const tol = info.tolerance, err = info.maxError;
                    if (info.final && isFinite(err)) {
                        lastFinalErr = err;
                        if (refineErr0 === null && err > tol) refineErr0 = err;
                        if (err <= tol) {
                            lastErrFrac = 1;
                        } else if (refineErr0 !== null && refineErr0 > tol) {
                            const num = Math.log10(refineErr0) - Math.log10(err);
                            const den = Math.log10(refineErr0) - Math.log10(tol);
                            if (den > 0) lastErrFrac = Math.max(0, Math.min(1, num / den));
                        }
                    }
                    const maxIt = info.maxIterations || 8;
                    const within = info.midpointsTotal > 0
                        ? Math.min(1, info.midpointsDone / info.midpointsTotal) : 1;
                    const bandLo = Math.min(1, (info.iteration - 1) / maxIt);
                    const bandHi = Math.min(1, info.iteration / maxIt);
                    const iterFloor = bandLo + (bandHi - bandLo) * within;
                    const convFrac = Math.max(lastErrFrac, iterFloor * 0.9) * 0.98;
                    const errPart = lastFinalErr !== null ? ` · error ${(lastFinalErr * 100).toFixed(3)}%` : '';
                    setSweep(Math.max(pointFloor, convFrac), `${info.totalSamples} points solved${errPart}`);
                }
                if (info.iteration > 0) {
                    const dc = hasDC ? sweepResults.find(r => r.freq === 0) : null;
                    const live = dc ? [dc, ...sweep.buildResults(nonZeroFreqs)]
                                    : sweep.buildResults(nonZeroFreqs);
                    live.sort((a, b) => a.freq - b.freq);
                    post({ id: currentId, type: 'partial', sweepResults: stripSweep(live) });
                }
            },
            shouldStop,
        });

        if (!stopRequested) {
            const interpResults = sweep.buildResults(nonZeroFreqs);
            const dc = hasDC ? sweepResults.find(r => r.freq === 0) : null;
            sweepResults = dc ? [dc, ...interpResults] : interpResults;
            log(`Interpolating sweep: ${nSamples + (hasDC ? 1 : 0)} exact solves for ${frequencies.length} output points`);
        }
        sweepWarnings = sweep.warnings || null;
    } else {
        log(`Calculating frequency sweep (${frequencies.length} points)...`);
        for (let i = 0; i < frequencies.length; i++) {
            const freq = frequencies[i];
            if (freq === maxFreq) continue;
            // Hand the worker's event loop a turn so the pending `stop` message
            // is actually delivered. `stop` arrives as a message and nothing in
            // this loop yields to the task queue on its own: the FDM backend's
            // only yield lives inside solve_adaptive, and TriBackend.solveAt is
            // synchronous.  Without this the flag is not read until the whole
            // sweep has finished, and Stop silently does nothing. (The
            // interpolating branch is safe sweep.run yields once per exact
            // solve.)
            await new Promise(r => setTimeout(r, 0));
            // No log line here: the main thread prints "Simulation stopped by
            // user" for every stopped job, from the `stopped` flag in the
            // reply.
            if (stopRequested) break;

            const result = await solver.computeAtFrequency(freq, cachedResults);
            if (result && result.warnings)
                post({ id: currentId, type: 'warnings', warnings: result.warnings });
            sweepResults.push({ freq, result });

            progress(MESH_FRACTION + (i + 1) / frequencies.length * (1 - MESH_FRACTION),
                     `Frequency sweep: ${i + 1}/${frequencies.length} (${(freq / 1e9).toFixed(2)} GHz)`);

            if (i % 10 === 0) {
                const live = sweepResults.slice().sort((a, b) => a.freq - b.freq);
                post({ id: currentId, type: 'partial', sweepResults: stripSweep(live) });
            }
        }
    }

    sweepResults.sort((a, b) => a.freq - b.freq);
    return {
        stopped: stopRequested,
        sweepResults: stripSweep(sweepResults),
        // The converged mesh solve itself: the summary block reports Z_diff / Z_common /
        // RLGC_matrix from it, and an interpolated sweep row is not a substitute (it is
        // rebuilt from the interpolant, not from a full solve).
        meshResult: stripResult(results),
        fields: fieldPayload(solver),
        sweepWarnings,
    };
}

// ---------------------------------------------------------------- modes

async function jobModes({ params, freq, nev, refineOpts }) {
    modesSolver = makeSolver(params);
    const result = await modesSolver.solveModes(freq, nev, (pp) => {
        const line = `Refinement pass ${pp.iteration}/${pp.max_iterations}: ` +
            `param err=${(isFinite(pp.param_error) ? pp.param_error : 1).toExponential(2)}, ` +
            `${pp.n_tris || 0} triangles…`;
        progress(0, line);
        log(`Modes pass ${pp.iteration}: param error=` +
            `${(isFinite(pp.param_error) ? pp.param_error : 1).toExponential(3)}, Tris=${pp.n_tris || 0}`);
    }, refineOpts);
    // The per-mode field grids are fetched lazily (jobModeField): with up to 30 modes,
    // eagerly resampling and shipping every one of them would cost far more than the
    // handful the user actually clicks on.
    return { result };
}

function jobModeField({ idx }) {
    if (!modesSolver) return { grid: null };
    return { grid: modesSolver.getModeField(idx) || null };
}

// Parameter Sweep

// `points` is precomputed on the main thread: one fully-resolved params object per sweep
// value, because deriving them means reading the sidebar inputs (DOM), which the worker
// has no access to.
async function jobParamSweep({ points, freqHz, opts }) {
    const out = [];
    for (let i = 0; i < points.length; i++) {
        if (stopRequested) { log('Sweep stopped.'); break; }
        const { params, displayVal } = points[i];
        let solver;
        try {
            solver = makeSolver(params);
        } catch {
            log(`Point ${i + 1}: solver init failed, skipping.`);
            continue;
        }
        solver.ensure_mesh();
        const cached = await solver.solve_adaptive({
            max_iters: params.max_iters,
            energy_tol: params.tolerance,
            param_tol: 0.05,
            max_nodes: params.max_nodes * 1000,
            min_converged_passes: params.min_converged_passes,
            // Verify the first point only. Adjacent sweep points share the geometry
            // family and mesh behavior, so one certificate is representative, while
            // verifying every point would multiply the whole sweep by the overhead.
            certify: !!opts.estimateError && i === 0,
            onProgress: () => {},
            shouldStop,
        });
        if (stopRequested) { log('Sweep stopped.'); break; }

        if (i === 0 && opts.estimateError) {
            post({ id: currentId, type: 'firstPointCert',
                   certification: solver.certification || null,
                   warnings: (cached && cached.warnings) || null });
        }

        const result = await solver.computeAtFrequency(freqHz, cached);
        out.push({ paramValue: displayVal, result: stripResult(result) });
        post({ id: currentId, type: 'partial', sweepResults: out, index: i, total: points.length });
    }
    return { stopped: stopRequested, results: out };
}

// Dispatch

const JOBS = {
    simulate: jobSimulate,
    modes: jobModes,
    modeField: jobModeField,
    paramSweep: jobParamSweep,
};

async function runJob(msg) {
    const job = JOBS[msg.type];
    if (!job) { post({ id: msg.id, type: 'error', message: `unknown job "${msg.type}"` }); return; }

    currentId = msg.id;
    stopRequested = false;
    try {
        const payload = await job(msg);
        post({ id: msg.id, type: 'done', ...payload });
    } catch (err) {
        post({ id: msg.id, type: 'error', message: (err && err.message) || String(err) });
    } finally {
        currentId = null;
    }
}

// Jobs run one at a time. The handler is async, so without this chain a message
// arriving while a job is awaiting would start a second handler concurrently, and since
// each job begins by clearing `stopRequested` and claiming `currentId`, an overlap would
// silently cancel a pending Stop and misroute the first job's log/progress messages. Now
// that the UI stays interactive during a solve, overlapping requests are reachable (e.g.
// clicking a mode row while the main solve runs), so the ordering has to be explicit.
let jobQueue = Promise.resolve();

self.onmessage = (e) => {
    const msg = e.data;
    // Stop is control, not work: it must take effect on the job that is running now, so
    // it bypasses the queue entirely.
    if (msg.type === 'stop') { stopRequested = true; return; }
    jobQueue = jobQueue.then(() => runJob(msg));
};
