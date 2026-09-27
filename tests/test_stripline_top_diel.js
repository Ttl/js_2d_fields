// Stripline with a top-dielectric layer (air beside the trace), both backends.
//
// The app carves the layer out of the bottom of the stripline cover, so the
// ground-to-ground spacing is unchanged (tests/test_geometry.js pins that layout).
// This test checks the physics of the resulting geometry:
//   1. a layer with the cover's own material is the plain stripline: C and R agree
//      with the layerless solve (only the mesh differs, an extra interface line)
//   2. an air layer taller than the trace lowers eps_eff below the homogeneous value
//      and raises Z0, and both backends agree on C and Z0
//   3. a thicker air layer lowers eps_eff further (monotone in the layer height)
//   4. full-wave even mode of a diff stripline with an air cover (no layer, top
//      permittivity 1) agrees with quasi-static. The closed-wall eigen pick used to
//      return a strip/parallel-plate hybrid here (C 13% high, eigen-anchor warning);
//      floating grounds now route the eigensolves to radiating walls, see
//      TriBackend._eigenAbc, and the Modes tab reports them (checked here too)
//
// Run: node tests/test_stripline_top_diel.js
import { buildSolverFromParams } from '../src/solver_factory.js';
import { check, quiet, rel, APP, done } from './helpers.js';


const ER = 4.4, H = 0.3e-3, H_TOP = 0.3e-3, T = 35e-6;
const P = { tl_type: 'stripline', w: 0.2e-3, h: H, t: T, er: ER, tand: 0.02, er_top: ER, tand_top: 0.02,
    stripline_top_h: H_TOP, sigma: 5.8e7, freq: 1e9, nx: 30, ny: 30, rq: 0, use_causal_materials: false,
    use_top_diel: 0, top_diel_h: 0.1e-3, top_diel_er: 1, top_diel_tand: 0 };

async function solveAll(over, backend = 'rectilinear') {
    const errs = [];
    const s = buildSolverFromParams({ ...P, ...over, mesh_backend: backend }, e => errs.push(e));
    if (!s) throw new Error(errs.join('\n'));
    return quiet(() => s.solve_adaptive({ ...APP }));
}
const solve = async (over, backend) => (await solveAll(over, backend)).modes[0];

{
    console.log('1. layer of the cover material == plain stripline (QS)');
    const plain = await solve({});
    const same = await solve({ use_top_diel: 1, top_diel_h: 0.1e-3, top_diel_er: ER, top_diel_tand: 0.02 });
    check('C agrees within 0.5%', rel(plain.RLGC.C, same.RLGC.C) < 5e-3,
        `${(plain.RLGC.C * 1e12).toFixed(3)} vs ${(same.RLGC.C * 1e12).toFixed(3)} pF/m`);
    check('R agrees within 2%', rel(plain.RLGC.R, same.RLGC.R) < 0.02,
        `${plain.RLGC.R.toFixed(3)} vs ${same.RLGC.R.toFixed(3)} ohm/m`);
    check('G agrees within 0.5%', rel(plain.RLGC.G, same.RLGC.G) < 5e-3);
}

{
    console.log('2. air layer taller than the trace, both backends');
    const plain = await solve({});
    const air = { use_top_diel: 1, top_diel_h: 0.1e-3, top_diel_er: 1, top_diel_tand: 0 };
    const qs = await solve(air, 'rectilinear');
    const fw = await solve(air, 'triangular');
    check('QS eps_eff below the homogeneous stripline', qs.eps_eff < plain.eps_eff - 0.1,
        `${qs.eps_eff.toFixed(3)} vs ${plain.eps_eff.toFixed(3)}`);
    check('QS eps_eff above 1', qs.eps_eff > 1);
    check('QS Z0 above the homogeneous stripline', qs.Z0 > plain.Z0,
        `${qs.Z0.toFixed(2)} vs ${plain.Z0.toFixed(2)} ohm`);
    check('backends agree on C within 2%', rel(qs.RLGC.C, fw.RLGC.C) < 0.02,
        `qs ${(qs.RLGC.C * 1e12).toFixed(3)} vs tri ${(fw.RLGC.C * 1e12).toFixed(3)} pF/m`);
    check('backends agree on Z0 within 2%', rel(qs.Z0, fw.Z0) < 0.02,
        `qs ${qs.Z0.toFixed(2)} vs tri ${fw.Z0.toFixed(2)} ohm`);
    check('backends agree on R within 10%', rel(qs.RLGC.R, fw.RLGC.R) < 0.10,
        `qs ${qs.RLGC.R.toFixed(3)} vs tri ${fw.RLGC.R.toFixed(3)} ohm/m`);
    // Air fills the region beside the trace, so the shunt loss drops with it.
    check('G below the homogeneous stripline', qs.RLGC.G < plain.RLGC.G);
}

{
    console.log('3. eps_eff monotone in the air-layer height (QS)');
    const at = async (hL) => (await solve({ use_top_diel: 1, top_diel_h: hL, top_diel_er: 1, top_diel_tand: 0 })).eps_eff;
    const e1 = await at(0.02e-3), e2 = await at(0.1e-3), e3 = await at(0.2e-3);
    check('thin (trace pokes through) > tall > taller', e1 > e2 && e2 > e3,
        `${e1.toFixed(3)} > ${e2.toFixed(3)} > ${e3.toFixed(3)}`);
}

{
    console.log('4. full-wave even mode of a diff stripline with an air cover');
    const air = { tl_type: 'diff_stripline', trace_spacing: 0.2e-3, er_top: 1, tand_top: 0 };
    const qs = await solveAll(air, 'rectilinear');
    const fw = await solveAll(air, 'triangular');
    for (let k = 0; k < 2; k++) {
        check(`mode ${k} C agrees within 2%`, rel(qs.modes[k].RLGC.C, fw.modes[k].RLGC.C) < 0.02,
            `qs ${(qs.modes[k].RLGC.C * 1e12).toFixed(2)} vs tri ${(fw.modes[k].RLGC.C * 1e12).toFixed(2)} pF/m`);
    }
    check('no eigen-anchor warning', !(fw.warnings || []).some(w => w.type === 'eigen-anchor'),
        (fw.warnings || []).map(w => w.type).join(',') || 'none');
    // Modes tab: the same cross-section reports its floating grounds.
    const errs = [];
    const s = buildSolverFromParams({ ...P, ...air, mesh_backend: 'fullwave_mqs' }, e => errs.push(e));
    const modes = await quiet(() => s.solveModes(1e9, 4, null, { maxNodes: 8000 }));
    check('Modes tab warns about floating grounds', (modes.warnings || []).some(w => w.type === 'floating-grounds'));
    const boxed = buildSolverFromParams({ ...P, ...air, mesh_backend: 'fullwave_mqs', use_enclosure: 1,
        enclosure_width: 3e-3, use_side_gnd: 1, use_top_gnd: 0 }, e => errs.push(e));
    const modesBoxed = await quiet(() => boxed.solveModes(1e9, 4, null, { maxNodes: 8000 }));
    check('no floating-grounds warning with enclosure side walls', !(modesBoxed.warnings || []).some(w => w.type === 'floating-grounds'));
}

done();
