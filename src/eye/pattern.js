// Bit patterns for the time-domain simulation: PRBS sequences and custom bit strings.

// PRBS generator polynomials x^n + x^k + 1 (ITU-T O.150 and common SerDes use).
const PRBS = {
    prbs7: [7, 6], prbs9: [9, 5], prbs11: [11, 9], prbs15: [15, 14], prbs23: [23, 18], prbs31: [31, 28],
};

/**
 * Bits of a PRBS: nbits of the sequence of x^n + x^k + 1 from the all-ones state (a
 * Fibonacci LFSR, output = the feedback bit). The full period is 2^n - 1 bits.
 */
function prbsBits(name, nbits) {
    const p = PRBS[name];
    if (!p) throw new Error(`Unknown pattern "${name}".`);
    const [n, k] = p;
    const out = new Uint8Array(nbits);
    // state bit i = stage i + 1 of the register; 32-bit integer is enough up to n = 31.
    let s = (n === 31 ? 0x7fffffff : (1 << n) - 1) >>> 0;
    for (let i = 0; i < nbits; i++) {
        const b = ((s >>> (n - 1)) ^ (s >>> (k - 1))) & 1;
        s = ((s << 1) | b) & (n === 31 ? 0x7fffffff : (1 << n) - 1);
        out[i] = b;
    }
    return out;
}

// Full period of a PRBS.
const prbsPeriod = (name) => PRBS[name] ? 2 ** PRBS[name][0] - 1 : 0;

// A string of 0 and 1 (spaces and underscores ignored) as bits.
function parseBitString(str) {
    const s = String(str || '').replace(/[\s_]/g, '');
    if (!s) throw new Error('Custom pattern: enter the bits, e.g. 0011010111.');
    if (!/^[01]+$/.test(s)) throw new Error('Custom pattern: only 0 and 1 are allowed.');
    return Uint8Array.from(s, c => c === '1' ? 1 : 0);
}

// PRBS too long for a full period: the number of bits simulated is set by the user.
const PRBS_TRUNCATED = new Set(['prbs23', 'prbs31']);

/**
 * Bits of the pattern setting: { pattern: 'prbs7'..'prbs31' | 'custom', bits (custom
 * string), nbits (bits simulated of PRBS23 / PRBS31, the shorter ones run a full period) }.
 */
function patternBits({ pattern, bits, nbits }) {
    let b;
    if (pattern === 'custom') b = parseBitString(bits);
    else if (!PRBS_TRUNCATED.has(pattern) || !nbits) b = prbsBits(pattern, prbsPeriod(pattern));
    else b = prbsBits(pattern, Math.min(prbsPeriod(pattern), Math.max(2, Math.round(nbits))));
    let ones = 0;
    for (const x of b) ones += x;
    if (ones === 0 || ones === b.length) throw new Error('The pattern needs both 0 and 1 bits.');
    return b;
}

export { PRBS, PRBS_TRUNCATED, prbsBits, prbsPeriod, parseBitString, patternBits };
