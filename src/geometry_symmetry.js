// Geometric symmetry test + modal-decomposition classifier for a differential pair.
//
// Decides whether swapping the two signal conductors is an exact symmetry of the whole
// cross-section (conductor geometry + dielectric stack). That is precisely the condition
// under which the odd/even modal decomposition is exact and the 4-port S-matrix has
// S22≡S11 — so when it holds we keep the odd/even path and expose no physical [C]/[L],
// and when it fails we expose the asymmetric physical matrices.
//
// This replaces the |C11−C22| capacitance-imbalance heuristic that was previously used to
// make the same decision. The symmetric/asymmetric distinction is a property of the INPUT
// geometry, not of the discretised solve: a truly symmetric line has C11=C22 exactly, but
// the FDM/FEM mesh breaks that equality by a resolution-dependent noise floor, so the
// numeric test needed a hand-tuned threshold that was mesh-dependent and blind to small
// real asymmetries. Reading symmetry from the geometry is exact and mesh-independent.
//
// Returns:
//   true  — symmetric: swapping the two signals leaves geometry + dielectric invariant
//   false — asymmetric: no such swap symmetry exists
//   null  — cannot decide from geometry (not exactly two signals, or missing data);
//           the caller should fall back to the numeric capacitance test.

import { eig2x2, mat2Mul, mat2Inv } from './matrix.js';

// All conductors/dielectrics are axis-aligned rectangles exposing x_min/x_max/y_min/y_max.
// A swap symmetry must be a mirror across an axis-aligned plane (vertical for an edge-coupled
// pair, horizontal for a broadside pair) — the perpendicular bisector of the two signals.

function rectsBBoxScale(...lists) {
    let lo = Infinity, hi = -Infinity, blo = Infinity, bhi = -Infinity;
    for (const list of lists) {
        for (const r of list) {
            if (r.x_min < lo) lo = r.x_min;
            if (r.x_max > hi) hi = r.x_max;
            if (r.y_min < blo) blo = r.y_min;
            if (r.y_max > bhi) bhi = r.y_max;
        }
    }
    return Math.max(hi - lo, bhi - blo);
}

// Every rect must have a mirror partner across the line {axis}={coord}, matching in material
// (keyOf) and in its span on the preserved (non-mirrored) axis. Generalises the x=0 mirror
// test used for half-domain meshing to an arbitrary axis-aligned mirror plane.
function mirrorInvariant(rects, keyOf, axis, coord, tol) {
    const q = v => Math.round(v / tol);   // tolerance-quantised bucket key
    const buckets = new Map();
    for (const r of rects) {
        const pLo = axis === 'x' ? r.y_min : r.x_min;   // preserved axis (unchanged by the mirror)
        const pHi = axis === 'x' ? r.y_max : r.x_max;
        const k = `${keyOf(r)}|${q(pLo)}|${q(pHi)}`;
        if (!buckets.has(k)) buckets.set(k, []);
        const mLo = axis === 'x' ? r.x_min : r.y_min;   // mirrored axis
        const mHi = axis === 'x' ? r.x_max : r.y_max;
        buckets.get(k).push([mLo, mHi]);
    }
    for (const spans of buckets.values()) {
        for (const [lo, hi] of spans) {
            const rLo = 2 * coord - hi, rHi = 2 * coord - lo;   // reflected span
            if (!spans.some(([p, s]) => Math.abs(p - rLo) < tol && Math.abs(s - rHi) < tol)) return false;
        }
    }
    return true;
}

export function conductorSwapSymmetric(conductors, dielectrics) {
    if (!conductors || !dielectrics) return null;
    const signals = conductors.filter(c => c.is_signal);
    if (signals.length !== 2) return null;            // not a two-conductor pair → can't define a swap
    const tol = rectsBBoxScale(conductors, dielectrics) * 1e-6;
    if (!(tol > 0)) return null;

    const [a, b] = signals;
    const wa = a.x_max - a.x_min, wb = b.x_max - b.x_min;
    const ha = a.y_max - a.y_min, hb = b.y_max - b.y_min;
    if (Math.abs(wa - wb) > tol || Math.abs(ha - hb) > tol) return false;   // unequal traces → asymmetric

    const acx = (a.x_min + a.x_max) / 2, acy = (a.y_min + a.y_max) / 2;
    const bcx = (b.x_min + b.x_max) / 2, bcy = (b.y_min + b.y_max) / 2;
    const dx = Math.abs(acx - bcx), dy = Math.abs(acy - bcy);

    // The swap mirror is the perpendicular bisector of the two signal centres. It must be
    // axis-aligned: the pair is separated along x (edge-coupled → vertical mirror) or along y
    // (broadside → horizontal mirror). Separation on both axes admits no axis-aligned mirror
    // that swaps them (a point reflection would, but that does not preserve a layered dielectric
    // stack) → asymmetric.
    let axis, coord;
    if (dx > tol && dy <= tol) { axis = 'x'; coord = (acx + bcx) / 2; }
    else if (dy > tol && dx <= tol) { axis = 'y'; coord = (acy + bcy) / 2; }
    else return false;

    const condKey = c => (c.is_signal ? 's' : 'g');   // signals are interchangeable; grounds map to grounds
    const dielKey = d => `${(d.epsilon_r ?? 0).toFixed(6)}:${(d.tan_delta ?? 0).toFixed(6)}`;
    return mirrorInvariant(conductors, condKey, axis, coord, tol)
        && mirrorInvariant(dielectrics, dielKey, axis, coord, tol);
}

// Is the geometry mirror-symmetric about x=0? (Mesh only the right half)
export function isXSymmetric(conductors, dielectrics, domainW) {
    const tol = domainW * 1e-6;
    // Shaped geometry can't be compared by rect spans, so a shape declares its own mirror
    // symmetry instead (CoaxSolver sets xSymmetric on polygons built with n % 4 === 0 and
    // phase 0, which puts vertices exactly on the y axis so the x >= 0 half is an exact
    // half).
    const shaped = (o) => !!o.shape;
    if (![...conductors, ...dielectrics].filter(shaped).every(o => o.shape.xSymmetric === true))
        return false;
    conductors = conductors.filter(o => !shaped(o));
    dielectrics = dielectrics.filter(o => !shaped(o));
    const mirrorOf = (list, extra) => {
        // every rect must have a mirror partner (xmin <=> -xmax) with same y + material
        const key = (r) => `${extra(r)}|${r.y_min.toFixed(12)}|${r.y_max.toFixed(12)}`;
        const buckets = new Map();
        for (const o of list) {
            const k = key(o);
            if (!buckets.has(k)) buckets.set(k, []);
            buckets.get(k).push([o.x_min, o.x_max]);
        }
        for (const spans of buckets.values()) {
            for (const [xmin, xmax] of spans) {
                const mlo = -xmax, mhi = -xmin;
                const found = spans.some(([a, b]) => Math.abs(a - mlo) < tol && Math.abs(b - mhi) < tol);
                if (!found) return false;
            }
        }
        return true;
    };
    const condOk = mirrorOf(conductors, c => (c.is_signal ? 's' + Math.abs(c.polarity || 0) : 'g'));
    const dielOk = mirrorOf(dielectrics, d => `${d.epsilon_r.toFixed(6)}:${(d.tan_delta || 0).toFixed(6)}`);
    return condOk && dielOk;
}

// Half-domain symmetry: mirror-symmetric about x=0 and, for a differential
// pair, no signal conductor cutting the plane. A broadside pair's traces are
// stacked at x=0, so both modes are x-symmetric there and the plane cannot
// separate even from odd those pairs must solve on the full domain.
export function halfDomainSymmetry(conductors, dielectrics, domainW, isDifferential) {
    const tol = domainW * 1e-6;
    const diffStraddles = !!isDifferential && conductors.some(
        c => c.is_signal && c.x_min < -tol && c.x_max > tol);
    return { diffStraddles,
             ok: !diffStraddles && isXSymmetric(conductors, dielectrics, domainW) };
}

// Modal-decomposition classifier for a differential pair, shared by BOTH backends
// (FieldSolver2D._solve_modal_differential and TriBackend._prepareStaticModal) so the
// symmetric/asymmetric decision, thresholds, and odd/even mode ordering can never drift
// between them — the same geometry must get the same 4-port S-matrix shape from either
// solver.
//
// Inputs: the physical 2×2 Maxwell capacitance matrices [C] (dielectric) and [C0]
// (vacuum), BOTH in F/m, plus the geometry lists for the exact swap-symmetry test.
// Diagonalises [C0]⁻¹[C] (eps_eff = eigenvalues, modal voltage ratios = eigenvectors)
// and applies the guard chain:
//
//   1. SYMMETRIC — swap symmetry from the geometry (exact, mesh-independent); when the
//      geometric test cannot decide (null: not exactly two signals / missing data) fall
//      back to the |C11−C22| capacitance imbalance against a mesh-noise floor (≈3e-6
//      fine … ≈1.3e-3 coarse; a 5% broadside height imbalance ≈7e-3, 10% ≈1.5e-2;
//      5e-3 separates the two). The odd/even decomposition is exact: no physMatrix.
//   2. DEGENERATE — asymmetric but velocity-degenerate modes (homogeneous dielectric,
//      relative eigenvalue separation < 2%): the eigenvectors of [C0]⁻¹[C] are
//      arbitrary, so no per-mode rebuild is possible (modalVecs null, physMatrix
//      without Tv). The asymmetric [C]/[L] still drive the MTL 4-port S-matrix; the
//      caller keeps its odd/even per-mode results (loss [R]/[G] use the symmetric
//      reconstruction).
//   3. MODAL — non-degenerate asymmetric pair: physMatrix {C, L, Tv} plus
//      modalVecs = [v_odd, v_even], ordered by differential character
//      ((v0·v1)/|v|² — opposite-sign components → odd), for the caller to rebuild
//      per-mode fields from the eigenvector combination of its per-trace solves.
//
// globalThis.__MODAL_FORCE__ = 'on' | 'off' is a test override: 'off' forces SYMMETRIC,
// 'on' skips the symmetry AND degeneracy gates.
//
// Returns { physMatrix, modalVecs }:
//   physMatrix — null (symmetric) or {C, L: (1/c²)[C0]⁻¹, Tv?: [v_odd, v_even]}
//   modalVecs  — null (symmetric or degenerate) or [v_odd, v_even]
export function classifyModalDecomposition(C, C0, conductors, dielectrics) {
    const { vals, vecs } = eig2x2(mat2Mul(mat2Inv(C0), C));
    const sep = Math.abs(vals[0] - vals[1]) / (Math.abs(vals[0]) + Math.abs(vals[1]) + 1e-30);
    const force = globalThis.__MODAL_FORCE__;
    const symGeo = conductorSwapSymmetric(conductors, dielectrics);
    const asym = Math.abs(C[0][0] - C[1][1]) / (Math.abs(C[0][0]) + Math.abs(C[1][1]) + 1e-30);
    const isSymmetric = force === 'off'
        || (force !== 'on' && (symGeo === null ? asym < 5e-3 : symGeo));
    if (isSymmetric) return { physMatrix: null, modalVecs: null };
    const C_LIGHT = 299792458;
    const kL = 1 / (C_LIGHT * C_LIGHT);
    const C0inv = mat2Inv(C0);
    const physMatrix = { C, L: [[C0inv[0][0] * kL, C0inv[0][1] * kL], [C0inv[1][0] * kL, C0inv[1][1] * kL]] };
    if (force !== 'on' && sep < 0.02) return { physMatrix, modalVecs: null };
    const diffScore = v => { const n = v[0] * v[0] + v[1] * v[1]; return n > 0 ? v[0] * v[1] / n : 0; };
    const order = diffScore(vecs[0]) <= diffScore(vecs[1]) ? [0, 1] : [1, 0];
    const modalVecs = [vecs[order[0]], vecs[order[1]]];
    physMatrix.Tv = modalVecs;
    return { physMatrix, modalVecs };
}
