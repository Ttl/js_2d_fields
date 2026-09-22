// Complex-symmetric direct solver unit test — no FEM, no mesh: drives
// solveComplexSymmetric (solve_complex_symmetric in eigen_solver.cpp, the
// unpivoted LDLᵀ on Eigen's SimplicialLDLT with the "fake-real" csym scalar) on
// systems shaped like the MQS eddy-current problem, K = S + jβM with S SPD and M
// positive semidefinite on a "metal" subset, and verifies every solution by an
// independent complex residual in JS. Also pins the safety net: an unsymmetric
// input and a zero-pivot symmetric input both come back correct through the
// SparseLU fallback (rc = 1), and a singular / malformed system is an error,
// never a silent garbage vector.
//
// Run: node src/tri_solver/tests/complex_symmetric_test.mjs
process.env.TRI_STATS = process.env.TRI_STATS || '1';   // read at fem_core import time
const { default: createModule } = await import('../../wasm_solver/eigen_solver.js');
const { createWasmHelpers } = await import('../fem_core.js');

let failures = 0;
function check(name, cond, detail = '') {
    console.log(`${cond ? '✓ PASS' : '✗ FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
    if (!cond) failures++;
}

// CSR from a Map of "i,j" -> {re, im}; both triangles stored, rows sorted.
function toCSR(N, entries) {
    const rows = Array.from({ length: N }, () => []);
    for (const [key, v] of entries) { const [i, j] = key.split(',').map(Number); rows[i].push([j, v]); }
    const rowPtr = new Int32Array(N + 1), colIdx = [], valRe = [], valIm = [];
    for (let i = 0; i < N; i++) {
        rows[i].sort((a, b) => a[0] - b[0]);
        for (const [j, v] of rows[i]) { colIdx.push(j); valRe.push(v.re); valIm.push(v.im); }
        rowPtr[i + 1] = colIdx.length;
    }
    return { rowPtr, colIdx: Int32Array.from(colIdx), valRe: Float64Array.from(valRe), valIm: Float64Array.from(valIm) };
}
function addEntry(entries, i, j, re, im) {
    const k = `${i},${j}`, e = entries.get(k) || { re: 0, im: 0 };
    e.re += re; e.im += im; entries.set(k, e);
}
// ||K x - b|| / ||b|| with x, b in the [re(N); im(N)] layout.
function relResidual(N, csr, x, b) {
    let rn = 0, bn = 0;
    for (let i = 0; i < N; i++) {
        let sr = 0, si = 0;
        for (let p = csr.rowPtr[i]; p < csr.rowPtr[i + 1]; p++) {
            const j = csr.colIdx[p], kr = csr.valRe[p], ki = csr.valIm[p];
            sr += kr * x[j] - ki * x[N + j];
            si += kr * x[N + j] + ki * x[j];
        }
        rn += (sr - b[i]) ** 2 + (si - b[N + i]) ** 2;
        bn += b[i] ** 2 + b[N + i] ** 2;
    }
    return Math.sqrt(rn / bn);
}
// Deterministic pseudo-random numbers (LCG) so failures reproduce.
let seed = 12345;
const rnd = () => (seed = (seed * 1103515245 + 12345) % 2147483648) / 2147483648 - 0.5;

// MQS-shaped system on a 2D grid graph: S = 5-point Laplacian (SPD with the
// Dirichlet "walls" folded into the diagonal), M = lumped mass on a metal block,
// β a skin-effect-sized ratio so the imaginary part dominates inside the metal.
function buildMQSLike(nx, ny, beta) {
    const N = nx * ny, entries = new Map();
    const id = (i, j) => i * nx + j;
    for (let i = 0; i < ny; i++) for (let j = 0; j < nx; j++) {
        const n = id(i, j);
        addEntry(entries, n, n, 4, 0);          // Dirichlet neighbours off-grid stay in the diagonal
        if (i > 0) addEntry(entries, n, id(i - 1, j), -1, 0);
        if (i < ny - 1) addEntry(entries, n, id(i + 1, j), -1, 0);
        if (j > 0) addEntry(entries, n, id(i, j - 1), -1, 0);
        if (j < nx - 1) addEntry(entries, n, id(i, j + 1), -1, 0);
        const metal = i >= ny / 3 && i < 2 * ny / 3 && j >= nx / 3 && j < 2 * nx / 3;
        if (metal) addEntry(entries, n, n, 0, beta * (0.5 + rnd() * 0.2));
    }
    return { N, csr: toCSR(N, entries) };
}

const M = await createModule();
const helpers = createWasmHelpers(M);
const stats = globalThis.__TRI_STATS__;
const lastLin = () => stats.lin[stats.lin.length - 1];

// ---- 1: MQS-like system, two RHS, twice (second call hits the pattern cache) ----
{
    const { N, csr } = buildMQSLike(60, 50, 80);
    const rhs = [0, 1].map(() => { const b = new Float64Array(2 * N); for (let i = 0; i < 2 * N; i++) b[i] = rnd(); return b; });
    for (const pass of [1, 2]) {
        const xs = helpers.solveComplexSymmetric(N, csr, rhs);
        const res = xs.map((x, r) => relResidual(N, csr, x, rhs[r]));
        check(`MQS-like N=${N}, pass ${pass}: both RHS solved to residual < 1e-12`, res.every(r => r < 1e-12),
            res.map(r => r.toExponential(1)).join(', '));
        check(`MQS-like pass ${pass}: solved by the LDLᵀ path (no LU fallback)`, lastLin().luFallback === false);
    }
    // Same pattern, different β (a sweep point): cached symbolic analysis, new values.
    const csr2 = { ...csr, valIm: csr.valIm.map(v => 3.3 * v) };
    const [x2] = helpers.solveComplexSymmetric(N, csr2, [rhs[0]]);
    check('MQS-like, re-weighted β on the cached pattern: residual < 1e-12',
        relResidual(N, csr2, x2, rhs[0]) < 1e-12, relResidual(N, csr2, x2, rhs[0]).toExponential(1));
}

// ---- 2: unsymmetric input is refused by LDLᵀ and solved by the LU fallback ----
{
    const { N, csr } = buildMQSLike(20, 20, 30);
    const p = csr.rowPtr[5];          // perturb one off-diagonal of row 5 only
    for (let q = p; q < csr.rowPtr[6]; q++) if (csr.colIdx[q] !== 5) { csr.valRe[q] += 0.3; break; }
    const b = new Float64Array(2 * N); for (let i = 0; i < 2 * N; i++) b[i] = rnd();
    const [x] = helpers.solveComplexSymmetric(N, csr, [b]);
    check('unsymmetric input: correct through the LU fallback', relResidual(N, csr, x, b) < 1e-12 && lastLin().luFallback === true,
        `res ${relResidual(N, csr, x, b).toExponential(1)}, luFallback=${lastLin().luFallback}`);
}

// ---- 3: symmetric but zero leading pivot: LDLᵀ breaks down, LU fallback is exact ----
{
    const N = 2, entries = new Map();
    addEntry(entries, 0, 1, 1, 0); addEntry(entries, 1, 0, 1, 0);
    const csr = toCSR(N, entries);
    const b = Float64Array.from([1, 2, 0, 0]);
    const [x] = helpers.solveComplexSymmetric(N, csr, [b]);
    check('zero-pivot symmetric system: x = [2, 1] via the LU fallback',
        Math.abs(x[0] - 2) < 1e-14 && Math.abs(x[1] - 1) < 1e-14 && lastLin().luFallback === true,
        `x=[${x[0]}, ${x[1]}], luFallback=${lastLin().luFallback}`);
}

// ---- 4: errors are errors, never garbage ----
{
    const N = 2, entries = new Map();
    for (const [i, j] of [[0, 0], [0, 1], [1, 0], [1, 1]]) addEntry(entries, i, j, 1, 0);   // singular
    const csr = toCSR(N, entries);
    let threw = false;
    try { helpers.solveComplexSymmetric(N, csr, [Float64Array.from([1, 0, 0, 0])]); } catch { threw = true; }
    check('singular system throws', threw);
    const bad = { ...csr, colIdx: Int32Array.from([0, 7, 0, 1]) };
    threw = false;
    try { helpers.solveComplexSymmetric(N, bad, [Float64Array.from([1, 0, 0, 0])]); } catch { threw = true; }
    check('out-of-range column index throws', threw);
    threw = false;
    try { helpers.solveComplexSymmetric(N, csr, [Float64Array.from([1, 0])]); } catch { threw = true; }
    check('RHS of the wrong length (not 2N) throws', threw);
}

console.log(failures === 0 ? '\nAll complex-symmetric solver checks passed.' : `\n${failures} check(s) FAILED.`);
process.exit(failures === 0 ? 0 : 1);
