// Transmitter edge jitter of a bit pattern, as the PCIe channel tolerancing parameters
// (6.0 base spec 8.5.1.4, Table 8-14) define it at the die pad: uncorrelated Rj and
// DjDD, and pulse-width Rj and DDj. All terms act per edge and are independent between
// edges:
//   Rj      Gaussian, sigma rj
//   DjDD    dual Dirac, +-djdd / 2
//   PW Rj   Gaussian, sigma pwRj / sqrt(2): the width of a UI (the difference of two
//           edges) then has sigma pwRj
//   PW DDj  dual Dirac, +-pwDj / 4: the width of a UI then spans pwDj peak to peak
// This is the split the spec's derivation uses (Fig 8-78: the PW Rj enters the edge Rj
// times 0.707, the PW DDj the edge DjDD halved). Every edge at one bit boundary (both
// legs of a pair, the FIR steps) moves by the same displacement.

// Q of BER 1e-12 (one side): TJ = DJ + 2 Q RJ = DJ + 14.07 RJ.
const Q_BER12 = 7.0344;

function jitterActive(j) {
    return !!j && (!!j.times || j.rj > 0 || j.djdd > 0 || j.pwRj > 0 || j.pwDj > 0);
}

// Per-edge totals: Gaussian sigma, dual-Dirac peak-to-peak, and the BER 1e-12 total
// (dual-Dirac estimate).
function jitterPerEdge(j) {
    const rj = Math.sqrt((j.rj || 0) ** 2 + 0.5 * (j.pwRj || 0) ** 2);
    const dj = (j.djdd || 0) + 0.5 * (j.pwDj || 0);
    return { rj, dj, tj: dj + 2 * Q_BER12 * rj };
}

// Seeded generator (mulberry32): the same pattern gives the same jitter every run.
function rng(seed) {
    let a = seed >>> 0;
    return () => {
        a = (a + 0x6D2B79F5) >>> 0;
        let t = a;
        t = Math.imul(t ^ (t >>> 15), t | 1);
        t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}

// Displacement of each of n bit boundaries (s, positive later). j.times overrides.
function boundaryJitter(n, j, seed = 1) {
    if (j.times) {
        if (j.times.length !== n) throw new Error('Jitter: one displacement per bit.');
        return Float64Array.from(j.times);
    }
    const u = rng(seed), out = new Float64Array(n);
    const sRj = Math.sqrt((j.rj || 0) ** 2 + 0.5 * (j.pwRj || 0) ** 2);
    const a = 0.5 * (j.djdd || 0), b = 0.25 * (j.pwDj || 0);
    const sign = () => (u() < 0.5 ? -1 : 1);
    for (let k = 0; k < n; k++) {
        // Box-Muller; 1 - u() keeps the logarithm finite.
        const g = Math.sqrt(-2 * Math.log(1 - u())) * Math.cos(2 * Math.PI * u());
        out[k] = sRj * g + a * sign() + b * sign();
    }
    return out;
}

export { Q_BER12, jitterActive, jitterPerEdge, boundaryJitter };
