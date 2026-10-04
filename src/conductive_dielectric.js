// Dielectrics with a conductivity (sigma, S/m), such as a doped silicon substrate.
//
// Both backends treat them in the quasi-static limit with the complex permittivity
//   eps* = eps_r (1 - j tand) - j sigma / (omega eps0)
// The static solve with eps* gives the complex line capacitance C* = C' - j C'', so the
// shunt admittance is Y = j omega C* = omega C'' + j omega C': C is frequency dependent
// (a substrate that conducts screens the field like a floating conductor at low
// frequency) and G = omega C''. The magnetic field is left alone: the material carries
// no eddy currents, which holds while the skin depth in it is large next to its size
// (conductiveDielectricWarning).

const EPS0 = 8.854187817e-12;
const MU0 = 4 * Math.PI * 1e-7;

export const hasConductiveDielectric = dielectrics => (dielectrics || []).some(d => d.sigma > 0);

// Lowest relaxation frequency sigma / (2 pi eps0 eps_r) of the conductive dielectrics.
function minRelaxationFrequency(dielectrics) {
    let fr = Infinity;
    for (const d of dielectrics || []) {
        if (d.sigma > 0) fr = Math.min(fr, d.sigma / (2 * Math.PI * EPS0 * d.epsilon_r));
    }
    return fr;
}

// Angular frequency the complex solve runs at. At f = 0 eps* is singular, the DC
// limit is taken a million times below the relaxation frequency instead.
export function conductiveOmega(dielectrics, f) {
    const fe = f > 0 ? f : 1e-6 * minRelaxationFrequency(dielectrics);
    return 2 * Math.PI * fe;
}

// Complex relative permittivity { re, im } of a material at angular frequency omega.
// dc drops the loss tangent: the DC limit solves at a small omega > 0, where omega times
// the polarization loss would still add to G.
export function complexEps(er, tand, sigma, omega, dc = false) {
    return { re: er, im: -((dc ? 0 : er * (tand || 0)) + (sigma > 0 ? sigma / (omega * EPS0) : 0)) };
}

// Accuracy warning for a conductive dielectric whose skin depth at f is not large next
// to its size in the solved domain (box: { x_min, x_max, y_min, y_max }). Eddy currents
// in it would screen the magnetic field and lower L, which the dielectric model does
// not include. Null when every conductive dielectric is fine.
export function conductiveDielectricWarning(dielectrics, f, box) {
    if (!(f > 0)) return null;
    let worst = null;
    for (const d of dielectrics || []) {
        if (!(d.sigma > 0)) continue;
        const w = Math.min(d.x_max, box.x_max) - Math.max(d.x_min, box.x_min);
        const h = Math.min(d.y_max, box.y_max) - Math.max(d.y_min, box.y_min);
        if (!(w > 0 && h > 0)) continue;
        const size = Math.min(w, h);
        const delta = Math.sqrt(2 / (2 * Math.PI * f * MU0 * d.sigma));
        const ratio = size / delta;
        if (ratio > SKIN_RATIO_MAX && (!worst || ratio > worst.ratio)) worst = { d, size, delta, ratio };
    }
    if (!worst) return null;
    const um = v => `${+(v * 1e6).toPrecision(3)} um`;
    return { type: 'accuracy', reason: 'conductive-dielectric', mode: 'all', freq: f, message:
        `The dielectric with sigma = ${worst.d.sigma} S/m is too conductive to be solved as a lossy ` +
        `dielectric at ${+(f / 1e9).toPrecision(4)} GHz: its skin depth ${um(worst.delta)} is not large ` +
        `next to its size ${um(worst.size)}. The eddy currents it carries are not modelled, so L, R and ` +
        `the losses may be inaccurate. Model it as a conductor if it conducts like one.` };
}

// Size over skin depth above which the eddy currents are no longer negligible: the
// field they induce is of order (size / delta)^2, about 4 % at 1/3.
const SKIN_RATIO_MAX = 1 / 3;
