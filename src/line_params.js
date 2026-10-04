// Line quantities from the per-unit-length R, L, G, C that both backends and the
// interpolating sweep report.
import { Complex } from './complex.js';

const C0 = 299792458;

// DC row of the line: the f -> 0 limit of Zc and of the phase eps_eff = (beta/k0)^2
// reported at f > 0, so the DC point joins the sweep.
//
// With a shunt conductance (a conducting dielectric) gamma = sqrt((R+jwL)(G+jwC)) tends to
// sqrt(RG) (1 + jw (L/R + C/G) / 2), so
//   Zc -> sqrt(R/G)
//   (beta/k0)^2 -> (c^2/4) (L/Zc + C Zc)^2
// The last is at least c^2 L C, equal on a distortionless line (R/L = G/C).
// Without G both run off to infinity as f -> 0: Zc is reported infinite and eps_eff takes
// the lossless value c^2 L C, which the phase value approaches while wL >> R.
// Returns { Zc, eps_eff }, Zc in ohm (Infinity without G).
export function dcLineLimit(R, L, G, C) {
    if (G > 0 && R > 0) {
        const Zc = Math.sqrt(R / G);
        const v = 0.5 * C0 * (L / Zc + C * Zc);
        return { Zc, eps_eff: v * v };
    }
    return { Zc: Infinity, eps_eff: C0 * C0 * L * C };
}

// Attenuation in Np/m, { alpha_c, alpha_d }: together the exact Re(gamma). alpha_d is the
// attenuation of the same line with ideal conductors (no R, the external L only), alpha_c
// what the conductors add to it. alpha_d is then independent of the metal and its finish,
// G Z0 / 2 on a low-loss line, and exact on an R-G line (a conducting substrate below its
// relaxation frequency), where G Z0 / 2 would overstate the loss several times. alpha_c
// is not negative: Re(gamma) grows with R and L. At DC the attenuation sqrt(RG) needs
// the conductor resistance and is all alpha_c.
export function lineAttenuation(R, L, G, C, omega, L_ext = L) {
    if (!(omega > 0)) return { alpha_c: G > 0 && R > 0 ? Math.sqrt(R * G) : 0, alpha_d: 0 };
    const Y = new Complex(G, omega * C);
    const alpha = new Complex(R, omega * L).mul(Y).sqrt().re;
    const alphaD = new Complex(0, omega * L_ext).mul(Y).sqrt().re;
    if (!(alpha > 0)) return { alpha_c: 0, alpha_d: 0 };
    const alpha_d = Math.min(Math.max(alphaD, 0), alpha);
    return { alpha_c: alpha - alpha_d, alpha_d };
}
