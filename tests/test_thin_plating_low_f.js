// Thin plating on a thin, poorly conducting trace below the skin transition: a
// 1.2 um 7.5e6 S/m plating on a 19 um 5e5 S/m trace, delta far above both.
//
//   1 the layered model on a finite slab (H = 0 at its back) matches the analytic
//     slab Rs (1+j) coth((1+j) d/delta) for one metal, and the DC inductance of a
//     sheet over a slab with the plating's own field energy counted.
//   2 the internal inductance the plating adds or removes is bounded by the slab:
//     with the plating, L at 100 kHz stays above the 30 GHz (skin limit) value on
//     both backends, and the two agree. The quasi-static excess was taken against a
//     semi-infinite bulk and read 150 nH/m below the external inductance.
//   3 plating roughness adds the inductance of a surface layer: 0.45 um on the plating
//     adds less than 1 um on the bare trace. The full-wave meshed plating layer took
//     its surface current from its volume loss at the whole trace's slab thickness,
//     8x too large.
//
// Run: node tests/test_thin_plating_low_f.js
import { calculate_Zrough_layered, slabCoth } from '../src/surface_roughness.js';
import { MicrostripSolver } from '../src/microstrip.js';
import { check, quiet, APP, done } from './helpers.js';

const MU0 = 4e-7 * Math.PI;

// --- 1 ---
{
    console.log('1 layered model on a finite slab');
    let worst = 0;
    for (const f of [1e3, 1e6, 1e9, 1e11]) {
        const s = 5e7, d = 10e-6, w = 2 * Math.PI * f, dlt = Math.sqrt(2 / (w * MU0 * s)), Rs = 1 / (s * dlt);
        const z = calculate_Zrough_layered(f, s, 0, s, 1e-6, 0, d), c = slabCoth(d / dlt);
        worst = Math.max(worst, Math.abs(z.re / (Rs * c.re) - 1), Math.abs(z.im / (Rs * c.im) - 1));
    }
    check('one metal: slab Zs within 0.1% of the analytic coth', worst < 1e-3, worst.toExponential(2));
    const sb = 5e5, sp = 7.5e6, tp = 1.17e-6, d = 9.6e-6, f = 1e3;
    const L = calculate_Zrough_layered(f, sb, 0, sp, tp, 0, d).im / (2 * Math.PI * f);
    // Linear H across the plating (K to H1) and the bulk (H1 to 0), H1/K the bulk share.
    const fr = sb * (d - tp) / (sp * tp + sb * (d - tp));
    const Lref = MU0 * ((d - tp) / 3 * fr * fr + tp * (1 + fr + fr * fr) / 3);
    check('sheet over slab: DC inductance within 1%', Math.abs(L / Lref - 1) < 0.01,
        `${(L * 1e12).toFixed(3)} vs ${(Lref * 1e12).toFixed(3)} pH`);
    const half = calculate_Zrough_layered(1e9, 5.8e7, 1e-6, 8.7e6, 3e-6);
    const inf = calculate_Zrough_layered(1e9, 5.8e7, 1e-6, 8.7e6, 3e-6, 1e-6, Infinity);
    check('slab = Infinity is the semi-infinite model, bit-identical', half.re === inf.re && half.im === inf.im);
}

// --- 2, 3 ---
const MS = { trace_width: 56e-6, substrate_height: 260e-6, trace_thickness: 19.3e-6, gnd_thickness: 28e-6,
    epsilon_r: 2.94, tan_delta: 0.01, sigma_cond: 5.07e5, freq: 1e9, rq: 0, nx: 30, ny: 30,
    boundaries: ['open', 'open', 'open', 'gnd'] };
const PL = { sigma: 7.5e6, thickness: 1.17e-6, rq: 0, top: true, sides: true, bottom: false, thick_corners: true };
const Lof = async (backend, extra) => {
    const s = new MicrostripSolver({ ...MS, ...extra, mesh_backend: backend });
    s.use_causal_materials = false;
    return quiet(async () => {
        const r = await s.solve_adaptive({ ...APP });
        const out = {};
        for (const f of [1e5, 30e9]) out[f] = (await s.computeAtFrequency(f, r)).modes[0].RLGC.L;
        return out;
    });
};
const res = {};
for (const backend of ['rectilinear', 'triangular']) {
    res[backend] = {
        bare: await Lof(backend, {}),
        bareRough: await Lof(backend, { rq: 1e-6 }),
        plated: await Lof(backend, { plating: PL }),
        platedRough: await Lof(backend, { plating: { ...PL, rq: 0.45e-6 } }),
    };
}
const nH = v => (v * 1e9).toFixed(1);
{
    console.log('2 plating inductance bounded by the slab');
    for (const [backend, r] of Object.entries(res)) {
        check(`${backend}: plated L at 100 kHz above its 30 GHz value`, r.plated[1e5] > r.plated[30e9],
            `${nH(r.plated[1e5])} vs ${nH(r.plated[30e9])} nH/m`);
        check(`${backend}: plating changes L at 100 kHz by less than 15%`,
            Math.abs(r.plated[1e5] / r.bare[1e5] - 1) < 0.15, `${nH(r.plated[1e5])} vs bare ${nH(r.bare[1e5])} nH/m`);
    }
    const q = res.rectilinear.plated[1e5], w = res.triangular.plated[1e5];
    check('quasi-static and full-wave plated L at 100 kHz within 5%', Math.abs(q / w - 1) < 0.05,
        `${nH(q)} vs ${nH(w)} nH/m`);
}
{
    console.log('3 plating roughness inductance');
    for (const [backend, r] of Object.entries(res)) {
        const dPl = r.platedRough[1e5] - r.plated[1e5], dBare = r.bareRough[1e5] - r.bare[1e5];
        check(`${backend}: 0.45 um plating rq adds less than 0.6x what 1 um bare rq adds`,
            dPl > 0 && dPl < 0.6 * dBare, `${nH(dPl)} vs ${nH(dBare)} nH/m`);
    }
}
done();
