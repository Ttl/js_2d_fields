import { Complex } from './complex.js';

/**
 * Natural cubic spline interpolation.
 * Given data points (x_i, y_i), builds a C2-continuous piecewise cubic.
 */
class CubicSpline {
    /**
     * @param {number[]} x - Sorted x values (strictly increasing)
     * @param {number[]} y - Corresponding y values
     */
    constructor(x, y) {
        const n = x.length;
        if (n < 2) throw new Error('CubicSpline requires at least 2 points');

        this.x = x;
        this.y = y;
        this.n = n;

        if (n === 2) {
            // Linear interpolation
            this.a = [y[0]];
            this.b = [(y[1] - y[0]) / (x[1] - x[0])];
            this.c = [0];
            this.d = [0];
            return;
        }

        // Compute coefficients using the tridiagonal algorithm
        const h = new Float64Array(n - 1);
        for (let i = 0; i < n - 1; i++) {
            h[i] = x[i + 1] - x[i];
        }

        // Set up tridiagonal system for second derivatives (natural spline: S''(x0) = S''(xn) = 0)
        const alpha = new Float64Array(n);
        for (let i = 1; i < n - 1; i++) {
            alpha[i] = (3 / h[i]) * (y[i + 1] - y[i]) - (3 / h[i - 1]) * (y[i] - y[i - 1]);
        }

        // Solve tridiagonal system
        const l = new Float64Array(n);
        const mu = new Float64Array(n);
        const z = new Float64Array(n);

        l[0] = 1;
        for (let i = 1; i < n - 1; i++) {
            l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
        }
        l[n - 1] = 1;

        // Back substitution
        this.c = new Float64Array(n);
        this.b = new Float64Array(n - 1);
        this.d = new Float64Array(n - 1);
        this.a = new Float64Array(n - 1);

        for (let j = n - 2; j >= 0; j--) {
            this.c[j] = z[j] - mu[j] * this.c[j + 1];
            this.b[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (this.c[j + 1] + 2 * this.c[j]) / 3;
            this.d[j] = (this.c[j + 1] - this.c[j]) / (3 * h[j]);
            this.a[j] = y[j];
        }
    }

    /**
     * Evaluate the spline at a given x value.
     * Clamps to the boundary values outside the data range.
     * @param {number} xv - x value to evaluate at
     * @returns {number}
     */
    evaluate(xv) {
        const { x, n } = this;

        // Clamp to range
        if (xv <= x[0]) return this.y[0];
        if (xv >= x[n - 1]) return this.y[n - 1];

        // Binary search for interval
        let lo = 0, hi = n - 2;
        while (lo < hi) {
            const mid = (lo + hi) >> 1;
            if (x[mid + 1] < xv) lo = mid + 1;
            else hi = mid;
        }

        const dx = xv - x[lo];
        return this.a[lo] + this.b[lo] * dx + this.c[lo] * dx * dx + this.d[lo] * dx * dx * dx;
    }
}

/**
 * Interpolating frequency sweep.
 *
 * Adaptively samples exact RLGC at a small number of frequencies,
 * builds cubic spline interpolation over log-frequency, then evaluates
 * interpolated RLGC at all output frequencies.
 */
class InterpolatingSweep {
    /**
     * @param {object} solver - MicrostripSolver instance
     * @param {object} cachedResults - Results from initial solve (for computeAtFrequency)
     * @param {object} options
     * @param {number} [options.tolerance=0.001] - Max relative error tolerance
     * @param {number} [options.maxPoints=200] - Safety cap on number of sample points
     * @param {number} [options.maxIterations=8] - Max refinement iterations
     * @param {number} [options.initialPoints=12] - Initial number of sample points
     */
    constructor(solver, cachedResults, options = {}) {
        this.solver = solver;
        this.cachedResults = cachedResults;
        this.tolerance = options.tolerance ?? 0.001;
        this.maxPoints = options.maxPoints ?? 200;
        this.maxIterations = options.maxIterations ?? 8;
        this.initialPoints = options.initialPoints ?? 12;

        // Map from log10(freq) -> { modes: [{ mode, R, L, G, C }] }
        this.samplePoints = new Map();
        this.splines = null; // Built after adaptive refinement
    }

    /**
     * Compute exact RLGC at a frequency and store it.
     * @param {number} freq - Frequency in Hz
     * @returns {object} - The result from computeAtFrequency
     */
    async _computeExact(freq) {
        const result = await this.solver.computeAtFrequency(freq, this.cachedResults);
        // Collect per-point solver warnings (e.g. the full-wave backend's eigensolve
        // failure / mode-ambiguity warnings) so the app can surface them after the
        // sweep — they would otherwise be silently discarded with the result object.
        if (result && result.warnings) {
            if (!this.warnings) this.warnings = [];
            this.warnings.push(...result.warnings);
        }
        const t = Math.log10(freq);
        const modeData = result.modes.map(m => ({
            mode: m.mode,
            R: m.RLGC.R,
            L: m.RLGC.L,
            G: m.RLGC.G,
            C: m.RLGC.C
        }));
        this.samplePoints.set(t, modeData);
        return result;
    }

    /**
     * Build cubic splines from current sample points.
     * Creates one spline per (mode, RLGC component) combination.
     */
    _buildSplines() {
        // Sort sample points by log-frequency
        const entries = [...this.samplePoints.entries()].sort((a, b) => a[0] - b[0]);
        const ts = entries.map(e => e[0]);
        const numModes = entries[0][1].length;

        this.splines = [];
        for (let mi = 0; mi < numModes; mi++) {
            const modeSplines = {};
            for (const key of ['R', 'L', 'G', 'C']) {
                const values = entries.map(e => e[1][mi][key]);
                modeSplines[key] = new CubicSpline(ts, values);
            }
            modeSplines.mode = entries[0][1][mi].mode;
            this.splines.push(modeSplines);
        }
    }

    /**
     * Evaluate interpolated RLGC at a log-frequency.
     * @param {number} t - log10(freq)
     * @returns {Array} - [{mode, R, L, G, C}]
     */
    _evaluateSplines(t) {
        return this.splines.map(ms => ({
            mode: ms.mode,
            R: ms.R.evaluate(t),
            L: ms.L.evaluate(t),
            G: ms.G.evaluate(t),
            C: ms.C.evaluate(t)
        }));
    }

    /**
     * Compute max relative error between exact and interpolated RLGC.
     * @param {Array} exact - [{mode, R, L, G, C}]
     * @param {Array} interpolated - [{mode, R, L, G, C}]
     * @returns {number} - Maximum relative error across all components and modes
     */
    _computeError(exact, interpolated) {
        let maxErr = 0;
        for (let mi = 0; mi < exact.length; mi++) {
            for (const key of ['R', 'L', 'G', 'C']) {
                const ex = exact[mi][key];
                const interp = interpolated[mi][key];
                // Floor prevents division by near-zero
                const floor = key === 'G' ? 1e-6 : 1e-12;
                const denom = Math.max(Math.abs(ex), floor);
                const err = Math.abs(ex - interp) / denom;
                if (err > maxErr) maxErr = err;
            }
        }
        return maxErr;
    }

    /**
     * Run the adaptive interpolating sweep.
     * @param {number} fMin - Minimum frequency in Hz (must be > 0)
     * @param {number} fMax - Maximum frequency in Hz
     * @param {object} [callbacks]
     * @param {function} [callbacks.onProgress] - Called with {iteration, totalSamples, maxError}
     * @param {function} [callbacks.shouldStop] - Returns true to abort
     * @returns {number} - Number of exact solves performed
     */
    async run(fMin, fMax, callbacks = {}) {
        const { onProgress, shouldStop } = callbacks;
        // Sweep-max hint for the triangular backend (see solve_sweep): sizes the
        // MQS skin band for the whole sweep so every sample reuses one skin mesh.
        this.solver._sweepFmax = fMax;
        const tMin = Math.log10(fMin);
        const tMax = Math.log10(fMax);

        // Step 1: Initial geometrically spaced points
        const nInit = Math.min(this.initialPoints, this.maxPoints);
        const initialTs = [];
        for (let i = 0; i < nInit; i++) {
            initialTs.push(tMin + (tMax - tMin) * i / (nInit - 1));
        }

        // Compute exact RLGC at initial points
        for (const t of initialTs) {
            if (shouldStop && shouldStop()) return this.samplePoints.size;
            const freq = Math.pow(10, t);
            await this._computeExact(freq);
            await new Promise(resolve => setTimeout(resolve, 0)); // Yield to UI
        }

        if (onProgress) {
            onProgress({ iteration: 0, totalSamples: this.samplePoints.size, maxError: Infinity });
        }

        // Step 2: Adaptive refinement with selective midpoint computation
        // Track which intervals need testing. Initially all intervals are untested.
        // After each iteration, only subdivisions of failed intervals are retested.
        let intervalsToTest = null; // null means "test all"

        for (let iter = 0; iter < this.maxIterations; iter++) {
            if (shouldStop && shouldStop()) break;

            // Build splines from current points
            this._buildSplines();

            // Get sorted sample log-frequencies
            const sortedTs = [...this.samplePoints.keys()].sort((a, b) => a - b);

            // Determine which midpoints to compute
            const midpointsToCompute = [];
            if (intervalsToTest === null) {
                // First iteration: test all intervals
                for (let i = 0; i < sortedTs.length - 1; i++) {
                    midpointsToCompute.push((sortedTs[i] + sortedTs[i + 1]) / 2);
                }
            } else {
                // Subsequent iterations: only test sub-intervals of previously failed intervals
                for (const [tLo, tHi] of intervalsToTest) {
                    // Find the sample points within [tLo, tHi] and compute midpoints
                    // of the sub-intervals created by the previous insertion
                    const subPoints = sortedTs.filter(t => t >= tLo - 1e-15 && t <= tHi + 1e-15);
                    for (let i = 0; i < subPoints.length - 1; i++) {
                        const tMid = (subPoints[i] + subPoints[i + 1]) / 2;
                        if (!this.samplePoints.has(tMid)) {
                            midpointsToCompute.push(tMid);
                        }
                    }
                }
            }

            if (midpointsToCompute.length === 0) break;

            // Compute exact RLGC at needed midpoints and check against current spline
            let maxError = 0;
            const failedIntervals = [];

            for (const tMid of midpointsToCompute) {
                if (shouldStop && shouldStop()) return this.samplePoints.size;
                // Enforce the safety cap DURING the loop, not just between iterations:
                // an iteration starting below the cap can otherwise queue hundreds of
                // midpoints and blow past maxPoints (e.g. 293 solves for a 200 cap).
                if (this.samplePoints.size >= this.maxPoints) break;

                // Evaluate interpolated value BEFORE computing exact (using current spline)
                const interpolated = this._evaluateSplines(tMid);

                // Compute exact
                const freq = Math.pow(10, tMid);
                await this._computeExact(freq);
                await new Promise(resolve => setTimeout(resolve, 0)); // Yield to UI

                const exact = this.samplePoints.get(tMid);
                const err = this._computeError(exact, interpolated);
                if (err > maxError) maxError = err;
                if (err > this.tolerance) {
                    // Find the interval boundaries around this midpoint
                    const sortedNow = [...this.samplePoints.keys()].sort((a, b) => a - b);
                    // Binary search for tMid in the sorted array
                    let idx = -1;
                    for (let j = 0; j < sortedNow.length; j++) {
                        if (Math.abs(sortedNow[j] - tMid) < 1e-15) { idx = j; break; }
                    }
                    if (idx > 0 && idx < sortedNow.length - 1) {
                        failedIntervals.push([sortedNow[idx - 1], sortedNow[idx + 1]]);
                    }
                }
            }

            if (onProgress) {
                onProgress({ iteration: iter + 1, totalSamples: this.samplePoints.size, maxError });
            }

            // Check convergence
            if (failedIntervals.length === 0) break;

            // Safety cap
            if (this.samplePoints.size >= this.maxPoints) break;

            intervalsToTest = failedIntervals;
        }

        // Build final splines
        this._buildSplines();
        return this.samplePoints.size;
    }

    /**
     * Evaluate interpolated RLGC at an array of frequencies and build
     * results in the same format as the discrete sweep.
     *
     * @param {number[]} frequencies - Output frequencies in Hz
     * @returns {Array<{freq: number, result: object}>} - Same format as frequencySweepResults
     */
    buildResults(frequencies) {
        if (!this.splines) {
            throw new Error('Must call run() before buildResults()');
        }

        const results = [];
        const solver = this.solver;

        for (const freq of frequencies) {
            const t = Math.log10(freq);
            const interpolatedModes = this._evaluateSplines(t);

            // Build a full result object matching computeAtFrequency output
            const modeResults = interpolatedModes.map(im => {
                const { R, L, G, C } = im;
                const omega = 2 * Math.PI * freq;

                // Compute derived quantities from RLGC (same formulas as field_solver.js rlgc())
                // Note: freq=0 is excluded by the useInterpolation guard (fMin > 0)
                const Z_num = new Complex(R, omega * L);
                const Z_den = new Complex(G, omega * C);
                const Zc = Z_num.div(Z_den).sqrt();

                // gamma = sqrt((R+jwL)(G+jwC)), eps_eff = (beta/k0)^2
                const gamma = Z_num.mul(Z_den).sqrt();
                const beta = gamma.im;
                const k0 = omega / 299792458.0;
                const eps_eff = (beta / k0) * (beta / k0);

                // alpha_c = R / (2 * Re(Zc))  in Np/m, convert to dB/m
                const alpha_c = 8.686 * R / (2 * Zc.re);

                // alpha_d from G: alpha_d = G * Re(Zc) / 2  in Np/m -> dB/m
                const alpha_d = 8.686 * G * Zc.re / 2;

                const alpha_total = alpha_c + alpha_d;

                // L_external/L_internal split is not available from interpolated RLGC
                // (would require the cached field data). Not used for S-params or plots.
                const L_external = L;
                const L_internal = 0;

                return {
                    mode: im.mode,
                    Z0: Zc.re, // Use real part of complex Zc as Z0
                    eps_eff,
                    C, C0: C, // C0 not meaningful for interpolated; use C as placeholder
                    RLGC: { R, L, G, C },
                    Zc,
                    alpha_c, alpha_d, alpha_total,
                    L_internal,
                    L_external
                };
            });

            const result = { modes: modeResults };

            if (solver.is_differential) {
                const odd = modeResults.find(m => m.mode === 'odd');
                const even = modeResults.find(m => m.mode === 'even');
                if (odd && even) {
                    result.Z_diff = 2 * odd.Z0;
                    result.Z_common = even.Z0 / 2;
                }
            }

            results.push({ freq, result });
        }

        return results;
    }
}

export { CubicSpline, InterpolatingSweep };
