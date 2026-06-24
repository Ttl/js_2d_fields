// Geometric symmetry test for a differential conductor pair.
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
