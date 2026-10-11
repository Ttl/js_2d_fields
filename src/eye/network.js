// Circuit model of the link: driver, AC coupling capacitors, the line and the receiver
// termination, as a small nodal admittance (MNA) system.
//
// Per leg (one single-ended, P and N for a pair):
//
//   V_drv --R_s--(a)--[C_tx]--(b)== line ==(c)--[C_rx]--(d)-- termination
//
// The driver is a Norton equivalent at (a). A DC-coupled side merges the nodes on both
// sides of its capacitor. The line enters as its exact 2- or 4-port admittance matrix.
import { cholesky, forwardSub, backwardSubT, jacobiEigen, solveReal, solveComplexInPlace } from './dense.js';

// Conductance to ground of a node without any DC path (1 GOhm).
const LEAK = 1e-9;
// Floor of the total series resistance of a conductor, so a lossless line stays a
// finite conductance at DC.
const R_FLOOR = 1e-6;

// ---------------------------------------------------------------- complex helpers

function csqrt(a, b) {
    const r = Math.hypot(a, b);
    if (r === 0) return [0, 0];
    const re = Math.sqrt(0.5 * (r + Math.abs(a)));
    const im = 0.5 * b / re;
    return a >= 0 ? [re, im] : [Math.abs(im), b >= 0 ? re : -re];
}
const cmul = (a, b, c, d) => [a * c - b * d, a * d + b * c];
function cdiv(a, b, c, d) {
    const den = c * c + d * d;
    return [(a * c + b * d) / den, (b * c - a * d) / den];
}
function cexp(a, b) {
    const m = Math.exp(a);
    return [m * Math.cos(b), m * Math.sin(b)];
}

// x coth x and x csch x as functions of q = x^2 (entire in q, so no branch choice), and
// their derivatives with respect to q. Re(x) >= 0 keeps the exponentials bounded.
function hcoth(qr, qi) {
    if (qr * qr + qi * qi < 1e-8) {
        // 1 + q/3 - q^2/45 + 2 q^3/945
        const [q2r, q2i] = cmul(qr, qi, qr, qi), [q3r, q3i] = cmul(q2r, q2i, qr, qi);
        return [1 + qr / 3 - q2r / 45 + 2 * q3r / 945, qi / 3 - q2i / 45 + 2 * q3i / 945];
    }
    const [xr, xi] = csqrt(qr, qi);
    const [er, ei] = cexp(-2 * xr, -2 * xi);
    const [cr, ci] = cdiv(1 + er, ei, 1 - er, -ei);
    return cmul(xr, xi, cr, ci);
}

function hcsch(qr, qi) {
    if (qr * qr + qi * qi < 1e-8) {
        // 1 - q/6 + 7 q^2/360 - 31 q^3/15120
        const [q2r, q2i] = cmul(qr, qi, qr, qi), [q3r, q3i] = cmul(q2r, q2i, qr, qi);
        return [1 - qr / 6 + 7 * q2r / 360 - 31 * q3r / 15120, -qi / 6 + 7 * q2i / 360 - 31 * q3i / 15120];
    }
    const [xr, xi] = csqrt(qr, qi);
    const [e1r, e1i] = cexp(-xr, -xi);
    const [e2r, e2i] = cmul(e1r, e1i, e1r, e1i);
    const [sr, si] = cdiv(2 * e1r, 2 * e1i, 1 - e2r, -e2i);   // csch x
    return cmul(xr, xi, sr, si);
}

function hcothDeriv(qr, qi) {
    if (qr * qr + qi * qi < 1e-8) {
        const [q2r, q2i] = cmul(qr, qi, qr, qi);
        return [1 / 3 - 2 * qr / 45 + 2 * q2r / 315, -2 * qi / 45 + 2 * q2i / 315];
    }
    // (coth x - x csch^2 x) / (2x)
    const [xr, xi] = csqrt(qr, qi);
    const [er, ei] = cexp(-2 * xr, -2 * xi);
    const [cr, ci] = cdiv(1 + er, ei, 1 - er, -ei);
    const [d2r, d2i] = cmul(1 - er, -ei, 1 - er, -ei);
    const [c2r, c2i] = cdiv(4 * er, 4 * ei, d2r, d2i);       // csch^2 x
    const [xc2r, xc2i] = cmul(xr, xi, c2r, c2i);
    return cdiv(cr - xc2r, ci - xc2i, 2 * xr, 2 * xi);
}

function hcschDeriv(qr, qi) {
    if (qr * qr + qi * qi < 1e-8) {
        const [q2r, q2i] = cmul(qr, qi, qr, qi);
        return [-1 / 6 + 7 * qr / 180 - 31 * q2r / 5040, 7 * qi / 180 - 31 * q2i / 5040];
    }
    // csch x (1 - x coth x) / (2x)
    const [xr, xi] = csqrt(qr, qi);
    const [e1r, e1i] = cexp(-xr, -xi);
    const [e2r, e2i] = cmul(e1r, e1i, e1r, e1i);
    const [sr, si] = cdiv(2 * e1r, 2 * e1i, 1 - e2r, -e2i);
    const [hr, hi] = hcoth(qr, qi);
    const [nr, ni] = cmul(sr, si, 1 - hr, -hi);
    return cdiv(nr, ni, 2 * xr, 2 * xi);
}

// f(Q) for a 2x2 complex matrix Q = [[a, b], [c, d]] (arrays [re, im]), by the Newton
// form f(l2) I + f[l1, l2] (Q - l2 I) over its eigenvalues.
function matFun2(Q, f, fp) {
    const [a, b, c, d] = Q;
    const hr = 0.5 * (a[0] + d[0]), hi = 0.5 * (a[1] + d[1]);
    const det = [a[0] * d[0] - a[1] * d[1] - (b[0] * c[0] - b[1] * c[1]),
                 a[0] * d[1] + a[1] * d[0] - (b[0] * c[1] + b[1] * c[0])];
    const [h2r, h2i] = cmul(hr, hi, hr, hi);
    const [sr, si] = csqrt(h2r - det[0], h2i - det[1]);
    const l1 = [hr + sr, hi + si], l2 = [hr - sr, hi - si];
    const f2 = f(l2[0], l2[1]);
    let dd;
    const gap = Math.hypot(sr, si) * 2, scale = Math.max(Math.hypot(l1[0], l1[1]), Math.hypot(l2[0], l2[1]));
    if (gap > 1e-5 * scale && gap > 1e-300) {
        const f1 = f(l1[0], l1[1]);
        dd = cdiv(f1[0] - f2[0], f1[1] - f2[1], 2 * sr, 2 * si);
    } else {
        dd = fp(hr, hi);
    }
    const m = (x, y) => cmul(dd[0], dd[1], x, y);
    const ma = m(a[0] - l2[0], a[1] - l2[1]), mb = m(b[0], b[1]);
    const mc = m(c[0], c[1]), md = m(d[0] - l2[0], d[1] - l2[1]);
    return [[f2[0] + ma[0], f2[1] + ma[1]], mb, mc, [f2[0] + md[0], f2[1] + md[1]]];
}

const m2mul = (A, B) => [
    add(cm(A[0], B[0]), cm(A[1], B[2])), add(cm(A[0], B[1]), cm(A[1], B[3])),
    add(cm(A[2], B[0]), cm(A[3], B[2])), add(cm(A[2], B[1]), cm(A[3], B[3])),
];
const cm = (x, y) => cmul(x[0], x[1], y[0], y[1]);
const add = (x, y) => [x[0] + y[0], x[1] + y[1]];
function m2inv(A) {
    const det = add(cm(A[0], A[3]), cm(A[1], A[2]).map(v => -v));
    const inv = (x) => cdiv(x[0], x[1], det[0], det[1]);
    return [inv(A[3]), inv(A[1].map(v => -v)), inv(A[2].map(v => -v)), inv(A[0])];
}

// ---------------------------------------------------------------- line admittance

/**
 * Admittance blocks of a uniform line of length len at frequency f.
 * v: line parameters at f ({R, L, G, C}, or the 2x2 entries R11 R12 R22 ...).
 * Returns { y11, y12 }: single-ended complex numbers [re, im], or for a pair 2x2
 * matrices [Y_00, Y_01, Y_10, Y_11] of complex numbers. The matrix is
 * [[y11, y12], [y12, y11]] over (near, far) ports.
 *
 * With Q = (G + jwC)(R + jwL) len^2, y11 = sqrt(Q) coth(sqrt(Q)) (R + jwL)^-1 / len and
 * y12 = -sqrt(Q) csch(sqrt(Q)) (R + jwL)^-1 / len, exact for any loss and at DC.
 */
function lineAdmittance(f, v, len, differential) {
    const w = 2 * Math.PI * f;
    const rmin = R_FLOOR / len;
    if (!differential) {
        const zr = Math.max(v.R, rmin), zi = w * v.L;
        const [qr0, qi0] = cmul(v.G, w * v.C, zr, zi);
        const qr = qr0 * len * len, qi = qi0 * len * len;
        const [zlr, zli] = [zr * len, zi * len];
        const h1 = hcoth(qr, qi), h2 = hcsch(qr, qi);
        return { y11: cdiv(h1[0], h1[1], zlr, zli), y12: cdiv(-h2[0], -h2[1], zlr, zli) };
    }
    const Z = [[Math.max(v.R11, rmin), w * v.L11], [v.R12, w * v.L12], [v.R12, w * v.L12], [Math.max(v.R22, rmin), w * v.L22]];
    const Y = [[v.G11, w * v.C11], [v.G12, w * v.C12], [v.G12, w * v.C12], [v.G22, w * v.C22]];
    const l2 = len * len;
    const Q = m2mul(Y, Z).map(x => [x[0] * l2, x[1] * l2]);
    const Zi = m2inv(Z).map(x => [x[0] / len, x[1] / len]);
    const H1 = matFun2(Q, hcoth, hcothDeriv), H2 = matFun2(Q, hcsch, hcschDeriv);
    const y11 = m2mul(H1, Zi), y12 = m2mul(H2, Zi).map(x => [-x[0], -x[1]]);
    return { y11, y12 };
}

// ---------------------------------------------------------------- netlist

/**
 * Circuit of the link.
 * cfg: { differential, rs, tx: { ac, c }, rx: { ac, c },
 *        term: { type: 'vtt' | 'rc' | 'r' | 'diff' | 'open', rt, vtt, ct, rdiff } }
 * Returns the node map, the real conductance and capacitance matrices without the
 * line, the input matrix (driver EMF of each leg to node currents) and the constant
 * V_TT source vector.
 */
function buildNetwork(cfg) {
    const diff = !!cfg.differential, nl = diff ? 2 : 1;
    const pos = (x, what) => { if (!(x > 0) || !isFinite(x)) throw new Error(`${what} must be a positive number.`); return x; };
    const rs = pos(cfg.rs, 'Driver resistance');
    const term = cfg.term || { type: 'open' };
    let n = 0;
    const a = [], b = [], c = [], d = [], mid = [];
    const gst = [], cst = [];   // [i, j, value], j = -1 is ground
    const vttI = [];            // [node, current] of the V_TT source
    for (let l = 0; l < nl; l++) {
        a[l] = n++;
        b[l] = cfg.tx && cfg.tx.ac ? n++ : a[l];
        c[l] = n++;
        d[l] = cfg.rx && cfg.rx.ac ? n++ : c[l];
        gst.push([a[l], -1, 1 / rs]);
        if (cfg.tx && cfg.tx.ac) cst.push([a[l], b[l], pos(cfg.tx.c, 'TX coupling capacitance')]);
        if (cfg.rx && cfg.rx.ac) cst.push([c[l], d[l], pos(cfg.rx.c, 'RX coupling capacitance')]);
        switch (term.type) {
            case 'vtt': {
                const g = 1 / pos(term.rt, 'Termination resistance');
                gst.push([d[l], -1, g]);
                vttI.push([d[l], (term.vtt || 0) * g]);
                break;
            }
            case 'r':
                gst.push([d[l], -1, 1 / pos(term.rt, 'Termination resistance')]);
                break;
            case 'rc':
                mid[l] = n++;
                gst.push([d[l], mid[l], 1 / pos(term.rt, 'Termination resistance')]);
                cst.push([mid[l], -1, pos(term.ct, 'Termination capacitance')]);
                break;
            case 'diff':
                if (!diff) throw new Error('A differential termination needs a differential line.');
                break;
            case 'open':
                break;
            default:
                throw new Error(`Unknown termination "${term.type}".`);
        }
    }
    if (term.type === 'diff') gst.push([d[0], d[1], 1 / pos(term.rdiff, 'Differential termination resistance')]);

    // Nodes with a resistive path to ground or a source: the line conducts end to end.
    const parent = Array.from({ length: n + 1 }, (_, i) => i);   // n is ground
    const find = (i) => { while (parent[i] !== i) i = parent[i] = parent[parent[i]]; return i; };
    const unite = (i, j) => { parent[find(i)] = find(j); };
    for (const [i, j] of gst) unite(i, j < 0 ? n : j);
    for (let l = 0; l < nl; l++) unite(b[l], c[l]);
    const floating = [];
    for (let i = 0; i < n; i++) {
        if (find(i) !== find(n)) { floating.push(i); gst.push([i, -1, LEAK]); }
    }

    const G = new Float64Array(n * n), C = new Float64Array(n * n);
    const stamp = (M, [i, j, v]) => {
        M[i * n + i] += v;
        if (j >= 0) { M[j * n + j] += v; M[i * n + j] -= v; M[j * n + i] -= v; }
    };
    for (const s of gst) stamp(G, s);
    for (const s of cst) stamp(C, s);
    const B = new Float64Array(n * nl);
    for (let l = 0; l < nl; l++) B[a[l] * nl + l] = 1 / rs;
    const vtt = new Float64Array(n);
    for (const [i, cur] of vttI) vtt[i] += cur;

    return {
        n, nl, differential: diff, G, C, B, vtt,
        nodes: { a, b, c, d, mid },
        rxFloating: d.some(i => floating.includes(i)),
        floating,
        term, rs,
    };
}

// Stamp the line blocks (complex) into Ar, Ai.
function stampLine(net, Ar, Ai, Y) {
    const { n, nl } = net, { b, c } = net.nodes;
    const put = (i, j, z) => { Ar[i * n + j] += z[0]; Ai[i * n + j] += z[1]; };
    if (nl === 1) {
        put(b[0], b[0], Y.y11); put(c[0], c[0], Y.y11);
        put(b[0], c[0], Y.y12); put(c[0], b[0], Y.y12);
        return;
    }
    for (let p = 0; p < 2; p++) for (let q = 0; q < 2; q++) {
        const k = 2 * p + q;
        put(b[p], b[q], Y.y11[k]); put(c[p], c[q], Y.y11[k]);
        put(b[p], c[q], Y.y12[k]); put(c[p], b[q], Y.y12[k]);
    }
}

/**
 * Transfer from the driver EMFs to the node voltages `outNodes` at frequency f, with the
 * line parameters v. Returns { re, im }: Float64Arrays of outNodes.length * nl entries
 * (row = output node, column = driver leg).
 */
function makeTransfer(net, len, outNodes) {
    const { n, nl } = net;
    const Ar = new Float64Array(n * n), Ai = new Float64Array(n * n);
    const Br = new Float64Array(n * nl), Bi = new Float64Array(n * nl);
    const no = outNodes.length;
    return function transfer(f, v, out) {
        const w = 2 * Math.PI * f;
        for (let i = 0; i < n * n; i++) { Ar[i] = net.G[i]; Ai[i] = w * net.C[i]; }
        stampLine(net, Ar, Ai, lineAdmittance(f, v, len, net.differential));
        Br.set(net.B); Bi.fill(0);
        solveComplexInPlace(Ar, Ai, n, Br, Bi, nl);
        const re = out ? out.re : new Float64Array(no * nl), im = out ? out.im : new Float64Array(no * nl);
        for (let o = 0; o < no; o++) for (let l = 0; l < nl; l++) {
            re[o * nl + l] = Br[outNodes[o] * nl + l];
            im[o * nl + l] = Bi[outNodes[o] * nl + l];
        }
        return { re, im };
    };
}

// Real nodal matrix at DC: conductances plus the DC line, symmetrized.
function dcMatrix(net, vdc, len) {
    const { n } = net;
    const Ar = Float64Array.from(net.G), Ai = new Float64Array(n * n);
    stampLine(net, Ar, Ai, lineAdmittance(0, vdc, len, net.differential));
    const S = new Float64Array(n * n);
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) S[i * n + j] = 0.5 * (Ar[i * n + j] + Ar[j * n + i]);
    return S;
}

// Node voltages at DC with the driver legs at the EMFs u (array of nl).
function dcSolve(net, Gdc, u) {
    const { n, nl } = net;
    const rhs = Float64Array.from(net.vtt);
    for (let i = 0; i < n; i++) for (let l = 0; l < nl; l++) rhs[i] += net.B[i * nl + l] * u[l];
    return solveReal(Gdc, n, rhs);
}

/**
 * Modes of the lumped low-frequency circuit (the line as its DC network, the
 * capacitors kept): C x' + G x = B u. With G = L L^T and L^-1 C L^-T = W diag(lambda) W^T,
 * mode i obeys lambda_i z_i' + z_i = beta_i u and x = P z.
 * Returns { lambda, P (n x n, column i = mode i), beta (n x nl) }.
 */
function lumpedModes(net, Gdc) {
    const { n, nl } = net;
    const L = cholesky(Gdc, n);
    if (!L) throw new Error('The circuit has no DC solution.');
    // Y = L^-1 C (column by column), then M = L^-1 Y^T = L^-1 C L^-T.
    const Y = new Float64Array(n * n);
    for (let j = 0; j < n; j++) {
        const col = new Float64Array(n);
        for (let i = 0; i < n; i++) col[i] = net.C[i * n + j];
        const y = forwardSub(L, n, col);
        for (let i = 0; i < n; i++) Y[i * n + j] = y[i];
    }
    const M = new Float64Array(n * n);
    for (let j = 0; j < n; j++) {
        const col = new Float64Array(n);
        for (let i = 0; i < n; i++) col[i] = Y[j * n + i];
        const y = forwardSub(L, n, col);
        for (let i = 0; i < n; i++) M[i * n + j] = y[i];
    }
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
        const s = 0.5 * (M[i * n + j] + M[j * n + i]); M[i * n + j] = M[j * n + i] = s;
    }
    const { values, vectors } = jacobiEigen(M, n);
    const lambda = Array.from(values, v => Math.max(0, v));
    const P = new Float64Array(n * n);
    for (let j = 0; j < n; j++) {
        const w = new Float64Array(n);
        for (let i = 0; i < n; i++) w[i] = vectors[i * n + j];
        const x = backwardSubT(L, n, w);
        for (let i = 0; i < n; i++) P[i * n + j] = x[i];
    }
    const beta = new Float64Array(n * nl);
    for (let l = 0; l < nl; l++) {
        const col = new Float64Array(n);
        for (let i = 0; i < n; i++) col[i] = net.B[i * nl + l];
        const y = forwardSub(L, n, col);
        for (let j = 0; j < n; j++) {
            let s = 0;
            for (let k = 0; k < n; k++) s += vectors[k * n + j] * y[k];
            beta[j * nl + l] = s;
        }
    }
    return { lambda, P, beta };
}

export { buildNetwork, lineAdmittance, makeTransfer, dcMatrix, dcSolve, lumpedModes, stampLine,
    hcoth, hcsch, hcothDeriv, hcschDeriv, matFun2, LEAK };
