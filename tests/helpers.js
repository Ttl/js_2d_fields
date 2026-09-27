// Shared test helpers: pass/fail bookkeeping, console silencing, error metrics
// and the app's default adaptive-solve settings.

let passes = 0, failures = 0;

export function check(name, ok, detail = '') {
    if (ok) passes++; else failures++;
    console.log(`${ok ? '✓ PASS' : '✗ FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
    return ok;
}

export const failureCount = () => failures;

// Print the summary and exit nonzero if any check failed.
export function done() {
    console.log(failures ? `\n✗ ${failures} of ${passes + failures} checks failed`
        : `\n✓ all ${passes} checks passed`);
    process.exit(failures ? 1 : 0);
}

// Run fn with console.log/warn silenced (solver progress output).
export async function quiet(fn) {
    const log = console.log, warn = console.warn;
    console.log = () => {}; console.warn = () => {};
    try { return await fn(); } finally { console.log = log; console.warn = warn; }
}

// Symmetric relative difference.
export const rel = (a, b) => Math.abs(a - b) / Math.max(Math.abs(a), Math.abs(b), 1e-30);
// Relative error against a reference b.
export const relErr = (a, b) => Math.abs(a - b) / Math.abs(b);

// Adaptive solve settings used by the app.
export const APP = { max_iters: 10, energy_tol: 0.01, param_tol: 0.05, max_nodes: 20000, min_converged_passes: 2 };
