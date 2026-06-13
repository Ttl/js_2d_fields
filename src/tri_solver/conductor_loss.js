// Conductor loss via perturbation with Galerkin-projected H (OpenParEM2D approach).
// Projects Ht and Hz into FEM spaces using ALL DOFs (no PEC elimination for H).
// PEC eliminates tangential E, not tangential H.
// Loss integral: α_c = Rs · ∮(|Ht|² + |Hz|²) dl / (4P)

import { csqrt, tripletsToCSR } from './fem_core.js';
import { ne1, ne2, ne3, nf1, nf2, nf3, nf4, nf5, nf6,
         lv, le, le3, lf3, lvGrad, leGrad, le3Grad, lf3Grad, triCoefficients,
         ne1Curl, ne2Curl, ne3Curl, nf1Curl, nf2Curl, nf3Curl, nf4Curl, nf5Curl, nf6Curl,
         QW, QL1, QL2, QL3, NQ, QW12, QL1_12, QL2_12, QL3_12, NQ12,
         getLzOffsets, computeTriP2Matrices, computeTriP3Matrices } from './tri_fem.js';

// --- Helper: evaluate Nedelec basis functions at a point for P2 or P3 ---
// Returns { Wx, Wy, curlW, nNed } with nNed = 8 (P2) or 15 (P3)
function evalNedelecBasis(coeff, order, edgeVerts, xq, yq) {
    const nNed = order >= 3 ? 15 : 8;
    const Wx = new Float64Array(nNed), Wy = new Float64Array(nNed), curlW = new Float64Array(nNed);
    for (let k = 0; k < 3; k++) {
        const [p, q] = edgeVerts[k];
        [Wx[k], Wy[k]] = ne1(coeff, p, q, xq, yq); curlW[k] = ne1Curl(coeff, p, q);
        [Wx[k+4], Wy[k+4]] = ne2(coeff, p, q, xq, yq); curlW[k+4] = ne2Curl(coeff, p, q, xq, yq);
        if (order >= 3) {
            [Wx[k+8], Wy[k+8]] = ne3(coeff, p, q, xq, yq); curlW[k+8] = ne3Curl(coeff, p, q, xq, yq);
        }
    }
    [Wx[3], Wy[3]] = nf1(coeff, xq, yq); curlW[3] = nf1Curl(coeff, xq, yq);
    [Wx[7], Wy[7]] = nf2(coeff, xq, yq); curlW[7] = nf2Curl(coeff, xq, yq);
    if (order >= 3) {
        [Wx[11], Wy[11]] = nf3(coeff, xq, yq); curlW[11] = nf3Curl(coeff, xq, yq);
        [Wx[12], Wy[12]] = nf4(coeff, xq, yq); curlW[12] = nf4Curl(coeff, xq, yq);
        [Wx[13], Wy[13]] = nf5(coeff, xq, yq); curlW[13] = nf5Curl(coeff, xq, yq);
        [Wx[14], Wy[14]] = nf6(coeff, xq, yq); curlW[14] = nf6Curl(coeff, xq, yq);
    }
    return { Wx, Wy, curlW, nNed };
}

// --- Helper: evaluate Lagrange basis functions at a point for P2 or P3 ---
// Returns { Nz, Gx, Gy, nLag } with nLag = 6 (P2) or 10 (P3)
function evalLagrangeBasis(coeff, order, edgeVerts, xq, yq) {
    const nLag = order >= 3 ? 10 : 6;
    const Nz = new Float64Array(nLag), Gx = new Float64Array(nLag), Gy = new Float64Array(nLag);
    for (let k = 0; k < 3; k++) {
        Nz[k] = lv(coeff, k, xq, yq);
        [Gx[k], Gy[k]] = lvGrad(coeff, k, xq, yq);
    }
    for (let k = 0; k < 3; k++) {
        const [p, q] = edgeVerts[k];
        Nz[k+3] = le(coeff, p, q, xq, yq);
        [Gx[k+3], Gy[k+3]] = leGrad(coeff, p, q, xq, yq);
        if (order >= 3) {
            Nz[k+6] = le3(coeff, p, q, xq, yq);
            [Gx[k+6], Gy[k+6]] = le3Grad(coeff, p, q, xq, yq);
        }
    }
    if (order >= 3) {
        Nz[9] = lf3(coeff, xq, yq);
        [Gx[9], Gy[9]] = lf3Grad(coeff, xq, yq);
    }
    return { Nz, Gx, Gy, nLag };
}

// --- Helper: compute full H DOF counts for variable-order mesh ---
// Returns { NH, NHz, htEdgeOff, htFaceOff, hzVertOff, hzEdgeMidOff, hzEdge3Off, hzFaceOff }
function computeHDofCounts(mesh, fm) {
    const { nEdges, nTris, nNodes } = mesh;
    const elemOrder = fm.elemOrder;

    // Ht DOFs: per edge (2 for P2, 3 for P3) + per face (2 for P2, 6 for P3)
    // Use max order touching each edge for full H space (H is NOT subject to minimum rule)
    const edgeMaxOrder = new Uint8Array(nEdges).fill(2);
    if (elemOrder) {
        for (let t = 0; t < nTris; t++) {
            const ord = elemOrder[t];
            for (let k = 0; k < 3; k++) {
                const e = mesh.triEdges[3*t+k];
                if (ord > edgeMaxOrder[e]) edgeMaxOrder[e] = ord;
            }
        }
    }

    // Edge Ht DOFs: variable per edge
    const htEdgeDofStart = new Int32Array(nEdges);
    let htEdgeTotal = 0;
    for (let e = 0; e < nEdges; e++) {
        htEdgeDofStart[e] = htEdgeTotal;
        htEdgeTotal += edgeMaxOrder[e] >= 3 ? 3 : 2;
    }

    // Face Ht DOFs: variable per triangle
    const htFaceDofStart = new Int32Array(nTris);
    let htFaceTotal = 0;
    for (let t = 0; t < nTris; t++) {
        htFaceDofStart[t] = htFaceTotal;
        const ord = elemOrder ? elemOrder[t] : 2;
        htFaceTotal += ord >= 3 ? 6 : 2;
    }
    const NH = htEdgeTotal + htFaceTotal;

    // Hz DOFs: vertex + edge midpoint + (P3: edge extra + face bubble)
    const hzEdgeMidStart = nNodes; // after vertices
    let hzEdge3Start = hzEdgeMidStart + nEdges; // after edge midpoints
    let hzEdge3Total = 0;
    const hzEdge3DofStart = new Int32Array(nEdges).fill(-1);
    if (elemOrder) {
        for (let e = 0; e < nEdges; e++) {
            if (edgeMaxOrder[e] >= 3) {
                hzEdge3DofStart[e] = hzEdge3Start + hzEdge3Total;
                hzEdge3Total++;
            }
        }
    }
    let hzFaceStart = hzEdge3Start + hzEdge3Total;
    let hzFaceTotal = 0;
    const hzFaceDofStart = new Int32Array(nTris).fill(-1);
    if (elemOrder) {
        for (let t = 0; t < nTris; t++) {
            if (elemOrder[t] >= 3) {
                hzFaceDofStart[t] = hzFaceStart + hzFaceTotal;
                hzFaceTotal++;
            }
        }
    }
    const NHz = nNodes + nEdges + hzEdge3Total + hzFaceTotal;

    return { NH, NHz, htEdgeTotal, htEdgeDofStart, htFaceDofStart, edgeMaxOrder,
             hzEdge3DofStart, hzFaceDofStart };
}

const MU0 = 4 * Math.PI * 1e-7;

// Corner saturation radius r₀ = SAT_KAPPA·δ used by the perturbation-loss corner
// correction (see satCornerIntegral). The SIBC regularization length is
// ℓ = |Zs|/(ωμ₀) = δ/√2; κ = 1 calibrated against the QEP SIBC eigenvalue on the
// microstrip and rectangular-coax test geometries.
const SAT_KAPPA = 1.0;

// --- Galerkin projection of Ht and Hz (also used for Poynting power) ---
// Ht: jωμ·Mt·Ht = ẑ×∇Ez + γ·(ẑ×Et)  →  Mt·Ht_s = RHS_t  (Ht_s = jωμ·Ht)
// Hz: jωμ·Mz·Hz = curl_z(Et)          →  Mz·Hz_s = RHS_z  (Hz_s = jωμ·Hz)
// ALL DOFs — H is not subject to PEC BCs.
// Uses faceF from freedom map to identify conductor-interior triangles (geometry-independent).
export function projectH(mesh, fm, vecRe, vecIm, gamma, freq, enrichment, wasmSolver) {
    const { nodes, tris, triEdges, triSigns, nTris, nEdges, nNodes } = mesh;
    const { edgeF, faceF, nodeF, edgeNodeF, nFreeTransverse, nFreeVertexDof,
            elemOrder, edgeF3, faceF3, edgeNodeF3, faceNodeF,
            nFreeEdgeNodeDof, nFreeEdgeNode3Dof, nFreeFaceNodeDof } = fm;
    const omega = 2 * Math.PI * freq;
    const omu = omega * MU0;
    const { lzOff, lzEdgeMidOff, lzEdge3Off, lzFaceNodeOff } = getLzOffsets(fm);
    const edgeVerts = [[0, 1], [1, 2], [2, 0]];

    // Full H DOF counts (variable for hp)
    const hDofs = computeHDofCounts(mesh, fm);
    const { NH, NHz, htEdgeDofStart, htFaceDofStart, edgeMaxOrder,
            hzEdge3DofStart, hzFaceDofStart } = hDofs;

    // Ht system
    const MtR = [], MtC = [], MtV = [];
    const rhsRe = new Float64Array(NH), rhsIm = new Float64Array(NH);

    // Hz system
    const MzR = [], MzC = [], MzV = [];
    const rhsHzRe = new Float64Array(NHz), rhsHzIm = new Float64Array(NHz);

    for (let t = 0; t < nTris; t++) {
        if (faceF[2*t] < 0) continue;

        const v0 = tris[3*t], v1 = tris[3*t+1], v2 = tris[3*t+2];
        const { coeff, Area } = triCoefficients(nodes, v0, v1, v2);
        const txs = [nodes[2*v0], nodes[2*v1], nodes[2*v2]];
        const tys = [nodes[2*v0+1], nodes[2*v1+1], nodes[2*v2+1]];
        const verts = [v0, v1, v2];
        const order = elemOrder ? elemOrder[t] : 2;
        const nNed = order >= 3 ? 15 : 8;
        const nLag = order >= 3 ? 10 : 6;

        // Ht DOFs: full numbering (variable per edge/face)
        const hD = new Int32Array(nNed), hS = new Float64Array(nNed).fill(1);
        for (let k = 0; k < 3; k++) {
            const eIdx = triEdges[3*t+k], base = htEdgeDofStart[eIdx];
            hD[k] = base; hS[k] = triSigns[3*t+k]; // ne1
            hD[k+4] = base + 1; // ne2
            if (order >= 3) { hD[k+8] = base + 2; hS[k+8] = triSigns[3*t+k]; } // ne3
        }
        const faceBlockStart = hDofs.htEdgeTotal;
        hD[3] = faceBlockStart + htFaceDofStart[t]; // nf1
        hD[7] = faceBlockStart + htFaceDofStart[t] + 1; // nf2
        if (order >= 3) {
            for (let k = 0; k < 4; k++)
                hD[11+k] = faceBlockStart + htFaceDofStart[t] + 2 + k; // nf3..6
        }

        // Hz DOFs: full numbering
        const hzD = new Int32Array(nLag);
        for (let k = 0; k < 3; k++) hzD[k] = verts[k];
        for (let k = 0; k < 3; k++) hzD[k+3] = nNodes + triEdges[3*t+k];
        if (order >= 3) {
            for (let k = 0; k < 3; k++) {
                const d3 = hzEdge3DofStart[triEdges[3*t+k]];
                if (d3 >= 0) hzD[k+6] = d3;
            }
            const fnd = hzFaceDofStart[t];
            if (fnd >= 0) hzD[9] = fnd;
        }

        // E DOFs: free (PEC eliminated) — for building RHS
        const eD = new Int32Array(nNed), eS = new Float64Array(nNed).fill(1);
        for (let k = 0; k < 3; k++) {
            const eIdx = triEdges[3*t+k];
            eD[k] = edgeF[2*eIdx]; eS[k] = triSigns[3*t+k];
            eD[k+4] = edgeF[2*eIdx+1];
            if (order >= 3) { eD[k+8] = edgeF3 ? edgeF3[eIdx] : -1; eS[k+8] = triSigns[3*t+k]; }
        }
        eD[3] = faceF[2*t]; eD[7] = faceF[2*t+1];
        if (order >= 3 && faceF3) {
            for (let k = 0; k < 4; k++) eD[11+k] = faceF3[4*t+k];
        }

        // Gather Et from eigenvector
        const etR = new Float64Array(nNed), etI = new Float64Array(nNed);
        for (let k = 0; k < nNed; k++) { const g=eD[k]; if(g>=0){etR[k]=eS[k]*vecRe[g]; etI[k]=eS[k]*vecIm[g];} }

        // Gather Ez from eigenvector
        const ezR = new Float64Array(nLag), ezI = new Float64Array(nLag);
        for (let k = 0; k < 3; k++) { const nf=nodeF[verts[k]]; if(nf>=0){ezR[k]=vecRe[lzOff+nf]; ezI[k]=vecIm[lzOff+nf];} }
        for (let k = 0; k < 3; k++) { const enf=edgeNodeF[triEdges[3*t+k]]; if(enf>=0){ezR[k+3]=vecRe[lzEdgeMidOff+enf]; ezI[k+3]=vecIm[lzEdgeMidOff+enf];} }
        if (order >= 3 && edgeNodeF3) {
            const signs = triSigns;
            for (let k = 0; k < 3; k++) {
                const enf3 = edgeNodeF3[triEdges[3*t+k]];
                if (enf3 >= 0) { ezR[k+6] = signs[3*t+k]*vecRe[lzEdge3Off+enf3]; ezI[k+6] = signs[3*t+k]*vecIm[lzEdge3Off+enf3]; }
            }
            if (faceNodeF) {
                const fnf = faceNodeF[t];
                if (fnf >= 0) { ezR[9] = vecRe[lzFaceNodeOff+fnf]; ezI[9] = vecIm[lzFaceNodeOff+fnf]; }
            }
        }

        const MtEl = new Float64Array(nNed * nNed);
        const rR = new Float64Array(nNed), rI = new Float64Array(nNed);
        const MzEl = new Float64Array(nLag * nLag);
        const rzR = new Float64Array(nLag), rzI = new Float64Array(nLag);

        const useQ12 = order >= 3;
        const nqp = useQ12 ? NQ12 : NQ;
        const qw = useQ12 ? QW12 : QW;
        const ql1 = useQ12 ? QL1_12 : QL1;
        const ql2 = useQ12 ? QL2_12 : QL2;
        const ql3 = useQ12 ? QL3_12 : QL3;

        for (let q = 0; q < nqp; q++) {
            const w = qw[q];
            const xq = txs[0]*ql1[q]+txs[1]*ql2[q]+txs[2]*ql3[q];
            const yq = tys[0]*ql1[q]+tys[1]*ql2[q]+tys[2]*ql3[q];

            const { Wx, Wy, curlW } = evalNedelecBasis(coeff, order, edgeVerts, xq, yq);
            const { Nz, Gx: LGx, Gy: LGy } = evalLagrangeBasis(coeff, order, edgeVerts, xq, yq);

            // ∇Ez and Et interpolation from eigenvector
            let dxR=0,dxI=0,dyR=0,dyI=0;
            for (let k=0;k<nLag;k++){dxR+=LGx[k]*ezR[k];dxI+=LGx[k]*ezI[k];dyR+=LGy[k]*ezR[k];dyI+=LGy[k]*ezI[k];}

            let exR=0,exI=0,eyR=0,eyI=0;
            for (let k=0;k<nNed;k++){exR+=Wx[k]*etR[k];exI+=Wx[k]*etI[k];eyR+=Wy[k]*etR[k];eyI+=Wy[k]*etI[k];}

            // Enrichment contributions
            if (enrichment) {
                const { corners, evalFn, vecRe: evR, vecIm: evI, nStd, nCorners } = enrichment;
                for (let ci = 0; ci < nCorners; ci++) {
                    const { dphidx, dphidy } = evalFn(xq, yq, corners[ci]);
                    if (dphidx === 0 && dphidy === 0) continue;
                    const ctR = evR[nStd + ci], ctI = evI[nStd + ci];
                    exR += dphidx*ctR; exI += dphidx*ctI;
                    eyR += dphidy*ctR; eyI += dphidy*ctI;
                    const czR = evR[nStd + nCorners + ci], czI = evI[nStd + nCorners + ci];
                    dxR += dphidx*czR; dxI += dphidx*czI;
                    dyR += dphidy*czR; dyI += dphidy*czI;
                }
            }

            // Ht RHS (γ-scaled convention)
            const fxR=-dyR+gamma.re*(-eyR)-gamma.im*(-eyI);
            const fxI=-dyI+gamma.re*(-eyI)+gamma.im*(-eyR);
            const fyR=dxR+gamma.re*exR-gamma.im*exI;
            const fyI=dxI+gamma.re*exI+gamma.im*exR;

            // curl_z(Et)
            let curlR=0, curlI=0;
            for (let k=0;k<nNed;k++){ curlR+=curlW[k]*etR[k]; curlI+=curlW[k]*etI[k]; }

            // Ht mass matrix and RHS
            for (let i=0;i<nNed;i++){
                rR[i]+=w*(Wx[i]*fxR+Wy[i]*fyR);
                rI[i]+=w*(Wx[i]*fxI+Wy[i]*fyI);
                for (let j=0;j<nNed;j++) MtEl[i*nNed+j]+=w*(Wx[i]*Wx[j]+Wy[i]*Wy[j]);
            }

            // Hz mass matrix and RHS
            for (let i=0;i<nLag;i++){
                rzR[i]+=w*Nz[i]*curlR;
                rzI[i]+=w*Nz[i]*curlI;
                for (let j=0;j<nLag;j++) MzEl[i*nLag+j]+=w*Nz[i]*Nz[j];
            }
        }
        for (let k=0;k<nNed*nNed;k++) MtEl[k]*=Area;
        for (let k=0;k<nNed;k++){rR[k]*=Area; rI[k]*=Area;}
        for (let k=0;k<nLag*nLag;k++) MzEl[k]*=Area;
        for (let k=0;k<nLag;k++){rzR[k]*=Area; rzI[k]*=Area;}

        // Assemble Ht system
        for (let li=0;li<nNed;li++){
            const gi=hD[li], si=hS[li];
            rhsRe[gi]+=si*rR[li]; rhsIm[gi]+=si*rI[li];
            for (let lj=0;lj<nNed;lj++){
                const gj=hD[lj], v=si*hS[lj]*MtEl[li*nNed+lj];
                if(v!==0){MtR.push(gi);MtC.push(gj);MtV.push(v);}
            }
        }

        // Assemble Hz system
        const hzS = new Float64Array(nLag).fill(1);
        if (order >= 3) { for (let k=0;k<3;k++) hzS[k+6] = triSigns[3*t+k]; }
        for (let li=0;li<nLag;li++){
            const gi=hzD[li];
            rhsHzRe[gi]+=hzS[li]*rzR[li]; rhsHzIm[gi]+=hzS[li]*rzI[li];
            for (let lj=0;lj<nLag;lj++){
                const gj=hzD[lj], v=hzS[li]*hzS[lj]*MzEl[li*nLag+lj];
                if(v!==0){MzR.push(gi);MzC.push(gj);MzV.push(v);}
            }
        }
    }

    // H DOFs are allocated for ALL edges/faces/nodes, but conductor-interior
    // triangles are skipped in the assembly. On meshes with conductor interiors
    // (meshConductorInterior), interior-only DOFs would leave empty rows →
    // structurally singular mass matrices. Pin them with unit diagonals (rhs is
    // zero there, so H = 0 inside the conductor — physically correct).
    {
        const touchedT = new Uint8Array(NH);
        for (const r of MtR) touchedT[r] = 1;
        for (let i = 0; i < NH; i++) if (!touchedT[i]) { MtR.push(i); MtC.push(i); MtV.push(1); }
        const touchedZ = new Uint8Array(NHz);
        for (const r of MzR) touchedZ[r] = 1;
        for (let i = 0; i < NHz; i++) if (!touchedZ[i]) { MzR.push(i); MzC.push(i); MzV.push(1); }
    }

    const csrHt = tripletsToCSR(MtR, MtC, MtV, NH);
    const csrHz = tripletsToCSR(MzR, MzC, MzV, NHz);

    const [hR, hI] = wasmSolver(NH, csrHt, [rhsRe, rhsIm]);
    const [hzR, hzI] = wasmSolver(NHz, csrHz, [rhsHzRe, rhsHzIm]);
    return { htRe: hR, htIm: hI, hzRe: hzR, hzIm: hzI, omu, NH, NHz, hDofs };
}

// --- Gather H DOFs for a triangle using full numbering ---
function gatherHt(triIdx, mesh, htRe, htIm, hDofs) {
    const { triEdges, triSigns, nEdges } = mesh;
    const order = hDofs ? (mesh.triEdges ? (() => {
        // Determine element order from face DOF count
        const faceDofs = hDofs.htFaceDofStart[triIdx + 1] !== undefined
            ? hDofs.htFaceDofStart[triIdx + 1] - hDofs.htFaceDofStart[triIdx]
            : 2;
        return faceDofs > 2 ? 3 : 2;
    })() : 2) : 2;
    // Simpler: check edgeMaxOrder of first edge
    const elemOrd = hDofs && hDofs.edgeMaxOrder
        ? Math.max(...[0,1,2].map(k => hDofs.edgeMaxOrder[triEdges[3*triIdx+k]])) : 2;
    const nNed = elemOrd >= 3 ? 15 : 8;
    const hR = new Float64Array(nNed), hI = new Float64Array(nNed);
    if (hDofs) {
        for (let k = 0; k < 3; k++) {
            const eIdx = triEdges[3*triIdx+k], s = triSigns[3*triIdx+k];
            const base = hDofs.htEdgeDofStart[eIdx];
            hR[k] = s*htRe[base]; hI[k] = s*htIm[base];
            hR[k+4] = htRe[base+1]; hI[k+4] = htIm[base+1];
            if (elemOrd >= 3 && hDofs.edgeMaxOrder[eIdx] >= 3) {
                hR[k+8] = s*htRe[base+2]; hI[k+8] = s*htIm[base+2];
            }
        }
        const fb = hDofs.htEdgeTotal + hDofs.htFaceDofStart[triIdx];
        hR[3]=htRe[fb]; hI[3]=htIm[fb];
        hR[7]=htRe[fb+1]; hI[7]=htIm[fb+1];
        if (elemOrd >= 3) {
            for (let k = 0; k < 4; k++) { hR[11+k]=htRe[fb+2+k]; hI[11+k]=htIm[fb+2+k]; }
        }
    } else {
        for (let k = 0; k < 3; k++) {
            const eIdx = triEdges[3*triIdx+k], s = triSigns[3*triIdx+k];
            hR[k] = s*htRe[2*eIdx]; hI[k] = s*htIm[2*eIdx];
            hR[k+4] = htRe[2*eIdx+1]; hI[k+4] = htIm[2*eIdx+1];
        }
        const fo = 2*nEdges;
        hR[3]=htRe[fo+2*triIdx]; hI[3]=htIm[fo+2*triIdx];
        hR[7]=htRe[fo+2*triIdx+1]; hI[7]=htIm[fo+2*triIdx+1];
    }
    return { hR, hI, nNed };
}

// --- Poynting power from projected Ht (consistent with loss integral) ---
export function computePoyntingFromProjectedH(mesh, fm, vecRe, vecIm, htRe, htIm, omu, hDofs) {
    const { nodes, tris, triEdges, triSigns, nTris, nEdges } = mesh;
    const { edgeF, faceF, nodeF, edgeNodeF, nFreeTransverse, nFreeVertexDof,
            elemOrder, edgeF3, faceF3, edgeNodeF3, faceNodeF,
            nFreeEdgeNodeDof, nFreeEdgeNode3Dof } = fm;
    const lzOff = nFreeTransverse;
    const edgeVerts = [[0,1],[1,2],[2,0]];
    let Pz = 0;

    for (let t = 0; t < nTris; t++) {
        if (faceF[2*t] < 0) continue;

        const v0=tris[3*t],v1=tris[3*t+1],v2=tris[3*t+2];
        const order = elemOrder ? elemOrder[t] : 2;
        const nNed = order >= 3 ? 15 : 8;
        const {coeff,Area}=triCoefficients(nodes,v0,v1,v2);
        const txs=[nodes[2*v0],nodes[2*v1],nodes[2*v2]];
        const tys=[nodes[2*v0+1],nodes[2*v1+1],nodes[2*v2+1]];

        // Gather Et from eigenvector (P2/P3 aware)
        const eD=new Int32Array(nNed), eS=new Float64Array(nNed).fill(1);
        for(let k=0;k<3;k++){
            const eIdx=triEdges[3*t+k]; eD[k]=edgeF[2*eIdx]; eS[k]=triSigns[3*t+k];
            eD[k+4]=edgeF[2*eIdx+1];
            if (order >= 3) { eD[k+8] = edgeF3 ? edgeF3[eIdx] : -1; eS[k+8] = triSigns[3*t+k]; }
        }
        eD[3]=faceF[2*t]; eD[7]=faceF[2*t+1];
        if (order >= 3 && faceF3) { for(let k=0;k<4;k++) eD[11+k]=faceF3[4*t+k]; }
        const etR=new Float64Array(nNed),etI=new Float64Array(nNed);
        for(let k=0;k<nNed;k++){const g=eD[k]; if(g>=0){etR[k]=eS[k]*vecRe[g]; etI[k]=eS[k]*vecIm[g];}}

        // Gather Ht from projection — use projected DOF order for basis evaluation
        const {hR,hI,nNed:hNed}=gatherHt(t, mesh, htRe, htIm, hDofs);
        const hOrder = hNed >= 15 ? 3 : 2;

        const useQ12 = order >= 3 || hOrder >= 3;
        const nqp = useQ12 ? NQ12 : NQ;
        const qw = useQ12 ? QW12 : QW;
        const ql1 = useQ12 ? QL1_12 : QL1, ql2 = useQ12 ? QL2_12 : QL2, ql3 = useQ12 ? QL3_12 : QL3;

        for(let q=0;q<nqp;q++){
            const w=qw[q];
            const xq=txs[0]*ql1[q]+txs[1]*ql2[q]+txs[2]*ql3[q];
            const yq=tys[0]*ql1[q]+tys[1]*ql2[q]+tys[2]*ql3[q];
            const {Wx,Wy}=evalNedelecBasis(coeff, hOrder, edgeVerts, xq, yq);

            let exR=0,exI=0,eyR=0,eyI=0;
            for(let k=0;k<nNed;k++){exR+=Wx[k]*etR[k];exI+=Wx[k]*etI[k];eyR+=Wy[k]*etR[k];eyI+=Wy[k]*etI[k];}

            let hxR=0,hxI=0,hyR=0,hyI=0;
            for(let k=0;k<nNed;k++){hxR+=Wx[k]*hR[k];hxI+=Wx[k]*hI[k];hyR+=Wy[k]*hR[k];hyI+=Wy[k]*hI[k];}
            const htxR=hxI/omu, htxI=-hxR/omu;
            const htyR=hyI/omu, htyI=-hyR/omu;

            const sz = 0.5*( exR*htyR+exI*htyI - eyR*htxR-eyI*htxI );
            Pz += w*Area*sz;
        }
    }
    return Pz;
}

// --- Direct evaluation of ∮(|Ht|² + |Hz|²) dl from E eigenvector ---
// Computes H at boundary quadrature points directly from E (including enrichment),
// bypassing the Galerkin projection that loses the corner singularity.
// Physical formula: jωμHt = ẑ×∇ez + ẑ×et (no γ, since et = γ·Et and γ cancels).
// Hz = curl_z(et)/(jωμ·γ) — enrichment doesn't contribute (curl of gradient = 0).
function evalH2dlDirect(mesh, fm, vecRe, vecIm, gamma, omu, isLossEdge, enrichment) {
    const { nodes, tris, edges, triEdges, triSigns, nEdges, nTris, nNodes } = mesh;
    const { edgeF, faceF, nodeF, edgeNodeF, nFreeTransverse, nFreeVertexDof } = fm;
    const lzOff = nFreeTransverse;
    const edgeVerts = [[0,1],[1,2],[2,0]];
    const gMag2 = gamma.re*gamma.re + gamma.im*gamma.im;

    const edgeToTri = new Array(nEdges);
    for (let e=0;e<nEdges;e++) edgeToTri[e]=[];
    for (let t=0;t<nTris;t++) for (let le=0;le<3;le++) edgeToTri[triEdges[3*t+le]].push(t);

    const GL3p=[0.11270,0.50000,0.88730], GL3w=[0.27778,0.44444,0.27778];
    const omu2 = omu*omu;
    let h2dl = 0;

    for (let e=0;e<nEdges;e++){
        if (!isLossEdge[e]) continue;

        let adj=-1;
        for(const t of edgeToTri[e]){
            if (faceF[2*t] >= 0) { adj=t; break; }
        }
        if(adj<0) continue;

        const n0=edges[2*e],n1=edges[2*e+1];
        const v0=tris[3*adj],v1=tris[3*adj+1],v2=tris[3*adj+2];
        const {coeff}=triCoefficients(nodes,v0,v1,v2);
        const x0=nodes[2*n0],y0=nodes[2*n0+1],tx=nodes[2*n1]-x0,ty=nodes[2*n1+1]-y0;
        const L=Math.sqrt(tx*tx+ty*ty);

        // Gather Et DOFs from eigenvector (PEC eliminated)
        const eD=new Int32Array(8), eS=new Float64Array(8).fill(1);
        for(let k=0;k<3;k++){
            const eIdx=triEdges[3*adj+k];
            eD[k]=edgeF[2*eIdx]; eS[k]=triSigns[3*adj+k];
            eD[k+4]=edgeF[2*eIdx+1]; eS[k+4]=1;
        }
        eD[3]=faceF[2*adj]; eD[7]=faceF[2*adj+1];
        const etR=new Float64Array(8),etI=new Float64Array(8);
        for(let k=0;k<8;k++){const g=eD[k]; if(g>=0){etR[k]=eS[k]*vecRe[g]; etI[k]=eS[k]*vecIm[g];}}

        // Gather Ez DOFs
        const verts=[v0,v1,v2];
        const ezR=new Float64Array(6), ezI=new Float64Array(6);
        for(let k=0;k<3;k++){const nf=nodeF[verts[k]]; if(nf>=0){ezR[k]=vecRe[lzOff+nf]; ezI[k]=vecIm[lzOff+nf];}}
        for(let k=0;k<3;k++){const enf=edgeNodeF[triEdges[3*adj+k]]; if(enf>=0){ezR[k+3]=vecRe[lzOff+nFreeVertexDof+enf]; ezI[k+3]=vecIm[lzOff+nFreeVertexDof+enf];}}

        for(let q=0;q<3;q++){
            const px=x0+GL3p[q]*tx, py=y0+GL3p[q]*ty;

            // Et at quadrature point
            const Wx=new Float64Array(8),Wy=new Float64Array(8);
            for(let k=0;k<3;k++){const [p,qq]=edgeVerts[k];[Wx[k],Wy[k]]=ne1(coeff,p,qq,px,py);[Wx[k+4],Wy[k+4]]=ne2(coeff,p,qq,px,py);}
            [Wx[3],Wy[3]]=nf1(coeff,px,py);[Wx[7],Wy[7]]=nf2(coeff,px,py);

            let exR=0,exI=0,eyR=0,eyI=0;
            for(let k=0;k<8;k++){exR+=Wx[k]*etR[k];exI+=Wx[k]*etI[k];eyR+=Wy[k]*etR[k];eyI+=Wy[k]*etI[k];}

            // ∇Ez at quadrature point
            let dxR=0,dxI=0,dyR=0,dyI=0;
            for(let k=0;k<3;k++){const [gx,gy]=lvGrad(coeff,k,px,py); dxR+=gx*ezR[k];dxI+=gx*ezI[k];dyR+=gy*ezR[k];dyI+=gy*ezI[k];}
            for(let k=0;k<3;k++){const [p,qq]=edgeVerts[k]; const [gx,gy]=leGrad(coeff,p,qq,px,py); dxR+=gx*ezR[k+3];dxI+=gx*ezI[k+3];dyR+=gy*ezR[k+3];dyI+=gy*ezI[k+3];}

            // Add enrichment contributions to Et and ∇Ez
            if (enrichment) {
                const { corners, evalFn, vecRe: evR, vecIm: evI, nStd, nCorners } = enrichment;
                for (let ci = 0; ci < nCorners; ci++) {
                    const { dphidx, dphidy } = evalFn(px, py, corners[ci]);
                    if (dphidx === 0 && dphidy === 0) continue;
                    const ctR = evR[nStd + ci], ctI = evI[nStd + ci];
                    exR += dphidx*ctR; exI += dphidx*ctI;
                    eyR += dphidy*ctR; eyI += dphidy*ctI;
                    const czR = evR[nStd + nCorners + ci], czI = evI[nStd + nCorners + ci];
                    dxR += dphidx*czR; dxI += dphidx*czI;
                    dyR += dphidy*czR; dyI += dphidy*czI;
                }
            }

            // Ht: jωμHt = ẑ×∇ez + ẑ×et = (-∂ez/∂y - ey, ∂ez/∂x + ex)
            const nHxR = -(dyR + eyR), nHxI = -(dyI + eyI);
            const nHyR = dxR + exR, nHyI = dxI + exI;
            const ht2 = (nHxR*nHxR + nHxI*nHxI + nHyR*nHyR + nHyI*nHyI) / omu2;

            // Hz: curl_z(et)/(jωμ·γ) — no enrichment contribution (curl∇ = 0)
            let curlR=0, curlI=0;
            for(let k=0;k<3;k++){
                const [p,qq]=edgeVerts[k];
                curlR+=ne1Curl(coeff,p,qq)*etR[k]; curlI+=ne1Curl(coeff,p,qq)*etI[k];
                const c2=ne2Curl(coeff,p,qq,px,py);
                curlR+=c2*etR[k+4]; curlI+=c2*etI[k+4];
            }
            curlR+=nf1Curl(coeff,px,py)*etR[3]+nf2Curl(coeff,px,py)*etR[7];
            curlI+=nf1Curl(coeff,px,py)*etI[3]+nf2Curl(coeff,px,py)*etI[7];
            // |Hz|² = |curl_z(et)|² / (|γ|² · (ωμ)²)
            const hz2 = (curlR*curlR + curlI*curlI) / (gMag2 * omu2);

            h2dl += GL3w[q]*L*(ht2 + hz2);
        }
    }
    return h2dl;
}

// --- Evaluate ∮(|Ht|² + |Hz|²) dl on loss boundary edges (from projected H) ---
function evalH2dl(mesh, fm, htRe, htIm, hzRe, hzIm, omu, isLossEdge, hDofs) {
    const { nodes, tris, edges, triEdges, nEdges, nTris, nNodes } = mesh;
    const { faceF, elemOrder } = fm;
    const edgeVerts = [[0,1],[1,2],[2,0]];

    const edgeToTri = new Array(nEdges);
    for (let e=0;e<nEdges;e++) edgeToTri[e]=[];
    for (let t=0;t<nTris;t++) for (let le=0;le<3;le++) edgeToTri[triEdges[3*t+le]].push(t);

    const GL3p=[0.11270,0.50000,0.88730], GL3w=[0.27778,0.44444,0.27778];
    const omu2 = omu*omu;
    let h2dl = 0;

    for (let e=0;e<nEdges;e++){
        if (!isLossEdge[e]) continue;

        let adj=-1;
        for(const t of edgeToTri[e]){
            if (faceF[2*t] >= 0) { adj=t; break; }
        }
        if(adj<0) continue;

        const n0=edges[2*e],n1=edges[2*e+1];
        // Use projected DOF counts (from gatherHt) for basis evaluation — the projected
        // H lives in the full space, so P2 elements adjacent to P3 need P3 evaluation.
        const {hR,hI,nNed:hNed} = gatherHt(adj, mesh, htRe, htIm, hDofs);
        const hOrder = hNed >= 15 ? 3 : 2;
        const nNed = hNed;
        const nLag = hOrder >= 3 ? 10 : 6;
        const v0=tris[3*adj],v1=tris[3*adj+1],v2=tris[3*adj+2];
        const {coeff}=triCoefficients(nodes,v0,v1,v2);
        const x0=nodes[2*n0],y0=nodes[2*n0+1],tx=nodes[2*n1]-x0,ty=nodes[2*n1+1]-y0;
        const L=Math.sqrt(tx*tx+ty*ty);

        // Gather Hz DOFs for the adjacent triangle
        const adjVerts = [tris[3*adj], tris[3*adj+1], tris[3*adj+2]];
        const hzR_loc = new Float64Array(nLag), hzI_loc = new Float64Array(nLag);
        for (let k=0;k<3;k++) { hzR_loc[k]=hzRe[adjVerts[k]]; hzI_loc[k]=hzIm[adjVerts[k]]; }
        for (let k=0;k<3;k++) { hzR_loc[k+3]=hzRe[nNodes+triEdges[3*adj+k]]; hzI_loc[k+3]=hzIm[nNodes+triEdges[3*adj+k]]; }
        if (hOrder >= 3 && hDofs) {
            const signs = mesh.triSigns;
            for (let k=0;k<3;k++) {
                const d3 = hDofs.hzEdge3DofStart[triEdges[3*adj+k]];
                if (d3 >= 0) { hzR_loc[k+6]=signs[3*adj+k]*hzRe[d3]; hzI_loc[k+6]=signs[3*adj+k]*hzIm[d3]; }
            }
            const fnd = hDofs.hzFaceDofStart[adj];
            if (fnd >= 0) { hzR_loc[9]=hzRe[fnd]; hzI_loc[9]=hzIm[fnd]; }
        }

        for(let q=0;q<3;q++){
            const px=x0+GL3p[q]*tx, py=y0+GL3p[q]*ty;

            const {Wx,Wy}=evalNedelecBasis(coeff, hOrder, edgeVerts, px, py);
            let hxR=0,hxI=0,hyR=0,hyI=0;
            for(let k=0;k<nNed;k++){hxR+=Wx[k]*hR[k];hxI+=Wx[k]*hI[k];hyR+=Wy[k]*hR[k];hyI+=Wy[k]*hI[k];}
            const ht2=(hxR*hxR+hxI*hxI+hyR*hyR+hyI*hyI)/omu2;

            const {Nz}=evalLagrangeBasis(coeff, hOrder, edgeVerts, px, py);
            let hzValR=0, hzValI=0;
            for (let k=0;k<nLag;k++) { hzValR+=Nz[k]*hzR_loc[k]; hzValI+=Nz[k]*hzI_loc[k]; }
            const hz2=(hzValR*hzValR+hzValI*hzValI)/omu2;

            h2dl += GL3w[q]*L*(ht2 + hz2);
        }
    }
    return h2dl;
}

// --- Analytical corner correction for ∮|Ht|²dl ---
// Near PEC corners, H ~ C·r^{ν-1} where ν depends on wedge angle/dielectric.
// The polynomial FEM overshoots |Ht|² on boundary edges touching corners.
// Fix: replace the FEM integral on corner-adjacent edges with the analytical
// integral ∫C²·r^{2(ν-1)}dr, fitting C from a neighboring edge where FEM is accurate.
//
// For corners at a dielectric interface (substrate/air), ν is computed from:
//   cos²(νπ/2) = εr / (2(1 + εr))
// For corners in free space: ν = 2/3.

// Compute singularity exponent for PEC corner at dielectric interface
function cornerNu(epsR) {
    if (!epsR || epsR <= 1.01) return 2/3;
    const target = epsR / (2 * (1 + epsR));
    return (2 / Math.PI) * Math.acos(Math.sqrt(target));
}

// Analytic integral of the saturated corner model over the corner edge [0, L]:
//   |H(r)|² = C²·r^{2(ν-1)}  for r > r₀   (PEC singular field)
//   |H(r)|² = C²·r₀^{2(ν-1)} for r ≤ r₀   (physical saturation at skin-depth scale)
// where u = 2ν-1. The PEC perturbation integral (r₀ = 0) overestimates the physical
// loss: the real field stops growing at r ~ δ while the PEC field continues as
// r^{ν-1}. With u as small as 0.117 (substrate corners at εr=4.4) the sub-δ region
// carries a large share of the PEC integral, which is why the uncorrected loss
// climbs far past the true (SIBC-validated) value on dense meshes.
function satCornerIntegral(C2, u, L, r0) {
    if (r0 <= 0) return C2 * Math.pow(L, u) / u;                       // PEC limit
    if (L <= r0) return C2 * L * Math.pow(r0, u - 1);                  // fully saturated
    return C2 * (Math.pow(r0, u) + (Math.pow(L, u) - Math.pow(r0, u)) / u);
}

export function evalH2dlCorrected(mesh, fm, htRe, htIm, hzRe, hzIm, omu, isLossEdge, hDofs,
                                  condRect, hSub, epsR, satR0 = 0) {
    const { nodes, edges, tris, triEdges, nEdges, nTris, nNodes } = mesh;
    const { faceF, isCondEdge } = fm;
    const edgeVertsLocal = [[0,1],[1,2],[2,0]];
    const TOL = 1e-10;
    const omu2 = omu*omu;

    // Build edge-to-triangle map
    const edgeToTri = new Array(nEdges);
    for (let e=0;e<nEdges;e++) edgeToTri[e]=[];
    for (let t=0;t<nTris;t++) for (let le=0;le<3;le++) edgeToTri[triEdges[3*t+le]].push(t);

    // --- Step 1: Compute per-edge h2dl and |Ht|² at 3 quadrature points ---
    const GL3p=[0.11270,0.50000,0.88730], GL3w=[0.27778,0.44444,0.27778];
    const edgeH2dl = new Float64Array(nEdges);
    // Store |Ht|² at each of the 3 Gauss points per edge (for fitting)
    const edgeHt2 = new Array(nEdges); // edgeHt2[e] = [ht2_q0, ht2_q1, ht2_q2] or null

    for (let e=0;e<nEdges;e++){
        if (!isLossEdge[e]) continue;
        let adj=-1;
        for(const t of edgeToTri[e]) if (faceF[2*t] >= 0) { adj=t; break; }
        if(adj<0) continue;

        const n0=edges[2*e],n1=edges[2*e+1];
        const {hR,hI,nNed:hNed} = gatherHt(adj, mesh, htRe, htIm, hDofs);
        const hOrder = hNed >= 15 ? 3 : 2;
        const nLag = hOrder >= 3 ? 10 : 6;
        const v0=tris[3*adj],v1=tris[3*adj+1],v2=tris[3*adj+2];
        const {coeff}=triCoefficients(nodes,v0,v1,v2);
        const x0=nodes[2*n0],y0=nodes[2*n0+1],tx=nodes[2*n1]-x0,ty=nodes[2*n1+1]-y0;
        const L=Math.sqrt(tx*tx+ty*ty);

        const adjVerts=[v0,v1,v2];
        const hzR_loc=new Float64Array(nLag), hzI_loc=new Float64Array(nLag);
        for(let k=0;k<3;k++){hzR_loc[k]=hzRe[adjVerts[k]]; hzI_loc[k]=hzIm[adjVerts[k]];}
        for(let k=0;k<3;k++){hzR_loc[k+3]=hzRe[nNodes+triEdges[3*adj+k]]; hzI_loc[k+3]=hzIm[nNodes+triEdges[3*adj+k]];}
        if(hOrder>=3&&hDofs){
            const signs=mesh.triSigns;
            for(let k=0;k<3;k++){const d3=hDofs.hzEdge3DofStart[triEdges[3*adj+k]]; if(d3>=0){hzR_loc[k+6]=signs[3*adj+k]*hzRe[d3]; hzI_loc[k+6]=signs[3*adj+k]*hzIm[d3];}}
            const fnd=hDofs.hzFaceDofStart[adj]; if(fnd>=0){hzR_loc[9]=hzRe[fnd]; hzI_loc[9]=hzIm[fnd];}
        }

        const ht2q = [0,0,0];
        for(let q=0;q<3;q++){
            const px=x0+GL3p[q]*tx, py=y0+GL3p[q]*ty;
            const {Wx,Wy}=evalNedelecBasis(coeff, hOrder, edgeVertsLocal, px, py);
            let hxR=0,hxI=0,hyR=0,hyI=0;
            for(let k=0;k<hNed;k++){hxR+=Wx[k]*hR[k];hxI+=Wx[k]*hI[k];hyR+=Wy[k]*hR[k];hyI+=Wy[k]*hI[k];}
            ht2q[q]=(hxR*hxR+hxI*hxI+hyR*hyR+hyI*hyI)/omu2;
            const {Nz}=evalLagrangeBasis(coeff, hOrder, edgeVertsLocal, px, py);
            let hzR2=0,hzI2=0;
            for(let k=0;k<nLag;k++){hzR2+=Nz[k]*hzR_loc[k]; hzI2+=Nz[k]*hzI_loc[k];}
            const hz2=(hzR2*hzR2+hzI2*hzI2)/omu2;
            edgeH2dl[e] += GL3w[q]*L*(ht2q[q]+hz2);
        }
        edgeHt2[e] = ht2q;
    }

    let h2dl = 0;
    for (let e = 0; e < nEdges; e++) h2dl += edgeH2dl[e];

    // --- Step 2: Identify PEC corners with correct ν ---
    // No condRect (e.g. circular geometry): no rectangular corners to correct.
    const condRects = condRect ? (condRect.rects || [condRect]) : [];
    const symX = (condRect && condRect.symmetry > 1) ? condRect.xmin_domain : null;
    const corners = [];
    for (const cr of condRects) {
        for (const p of [{x:cr.xmin,y:cr.ymin},{x:cr.xmax,y:cr.ymin},{x:cr.xmax,y:cr.ymax},{x:cr.xmin,y:cr.ymax}]) {
            if (symX !== null && Math.abs(p.x - symX) < TOL) continue;
            let nIdx = -1;
            for (let n=0;n<nNodes;n++) if(Math.abs(nodes[2*n]-p.x)<TOL && Math.abs(nodes[2*n+1]-p.y)<TOL){nIdx=n;break;}
            if (nIdx < 0) continue;
            // Corners on the dielectric interface (trace bottom at y=hSub) have the
            // dielectric-modified exponent ν = cornerNu(epsR) ≈ 0.559 for εr=4.4;
            // free-space corners (trace top) have ν = 2/3.
            const onInterface = (hSub !== undefined && epsR > 1.01 &&
                                 Math.abs(p.y - hSub) < TOL);
            const nu = onInterface ? cornerNu(epsR) : 2/3;
            corners.push({ x: p.x, y: p.y, nodeIdx: nIdx, nu });
        }
    }

    // --- Step 3: For each corner, correct edges on each surface ---
    let totalCorr = 0;
    for (const corner of corners) {
        // Find all signal conductor loss edges touching this corner node
        const touchingEdges = [];
        for (let e = 0; e < nEdges; e++) {
            if (!isLossEdge[e] || !isCondEdge[e]) continue;
            const n0 = edges[2*e], n1 = edges[2*e+1];
            if (n0 !== corner.nodeIdx && n1 !== corner.nodeIdx) continue;
            const other = n0 === corner.nodeIdx ? n1 : n0;
            const ox = nodes[2*other], oy = nodes[2*other+1];
            const dist = Math.sqrt((ox-corner.x)**2 + (oy-corner.y)**2);
            // Which endpoint of the edge is at the corner? n0 or n1?
            const cornerAtN0 = (n0 === corner.nodeIdx);
            touchingEdges.push({ e, other, dist, cornerAtN0 });
        }

        // Group by surface (same y = horizontal, same x = vertical)
        const surfaceGroups = new Map();
        for (const te of touchingEdges) {
            const cy = nodes[2*corner.nodeIdx+1], oy = nodes[2*te.other+1];
            const cx = nodes[2*corner.nodeIdx], ox = nodes[2*te.other];
            let key;
            if (Math.abs(cy-oy) < TOL) key = `y=${cy.toExponential(6)}`;
            else if (Math.abs(cx-ox) < TOL) key = `x=${cx.toExponential(6)}`;
            else key = `d=${te.e}`;
            if (!surfaceGroups.has(key)) surfaceGroups.set(key, []);
            surfaceGroups.get(key).push(te);
        }

        for (const [surfKey, surfEdges] of surfaceGroups) {
            surfEdges.sort((a,b) => a.dist - b.dist);
            const ce = surfEdges[0]; // Corner edge (touches the corner)
            const L = ce.dist;       // = distance from corner to far endpoint

            // Energy-based fitting: use the edge integral of the NEXT edge on this
            // surface to predict what the corner edge integral should be.
            //
            // If H ~ C·r^{ν-1}, the integral on edge [r₁,r₂] is:
            //   h2dl_edge = C²·(r₂^{2ν-1} - r₁^{2ν-1})/(2ν-1)
            //
            // The corner edge integral [0, L]:
            //   h2dl_corner = C²·L^{2ν-1}/(2ν-1)
            //
            // Ratio (1/(2ν-1) cancels!):
            //   h2dl_corner = h2dl_fit · L^{2ν-1} / (r₂^{2ν-1} - r₁^{2ν-1})

            // Two-parameter model: |H(r)|² = C²·r^{2(ν-1)} + C_reg²
            // (singular + regular parts). Use TWO fitting edges to solve a 2×2
            // system for C² and C_reg², then predict the corner edge integral.
            //
            // Edge [r₁,r₂] integral: I = C²·(r₂^u-r₁^u)/u + C_reg²·(r₂-r₁)
            // where u = 2ν-1.
            //
            // Corner edge [0,L]: I_corner = C²·L^u/u + C_reg²·L

            let fitNode = ce.other;
            let walkDist = ce.dist;
            const fitEdges = []; // {e, r1, r2} — non-corner edges on this surface

            for (let step = 0; step < 4; step++) {
                let foundNext = false;
                for (let e = 0; e < nEdges; e++) {
                    if (!isLossEdge[e] || !isCondEdge[e] || e === ce.e) continue;
                    const en0 = edges[2*e], en1 = edges[2*e+1];
                    if (en0 !== fitNode && en1 !== fitNode) continue;
                    const enOther = en0 === fitNode ? en1 : en0;
                    const eox = nodes[2*enOther], eoy = nodes[2*enOther+1];
                    if (surfKey.startsWith('y=') && Math.abs(eoy - nodes[2*fitNode+1]) > TOL) continue;
                    if (surfKey.startsWith('x=') && Math.abs(eox - nodes[2*fitNode]) > TOL) continue;
                    const eL = Math.sqrt((eox-nodes[2*fitNode])**2 + (eoy-nodes[2*fitNode+1])**2);
                    if (edgeH2dl[e] > 0) fitEdges.push({ e, r1: walkDist, r2: walkDist + eL });
                    walkDist += eL;
                    fitNode = enOther;
                    foundNext = true;
                    break;
                }
                if (!foundNext) break;
                if (fitEdges.length >= 2) break;
            }

            if (fitEdges.length < 1) continue;

            const nu = corner.nu;
            const u = 2*nu - 1;

            let h2dl_analytical;

            if (fitEdges.length >= 2) {
                // 2×2 system: [a₁ b₁; a₂ b₂]·[C²; C_reg²] = [I₁; I₂]
                const f0 = fitEdges[0], f1 = fitEdges[1];
                const a1 = (Math.pow(f0.r2, u) - Math.pow(f0.r1, u)) / u;
                const b1 = f0.r2 - f0.r1;
                const I1 = edgeH2dl[f0.e];
                const a2 = (Math.pow(f1.r2, u) - Math.pow(f1.r1, u)) / u;
                const b2 = f1.r2 - f1.r1;
                const I2 = edgeH2dl[f1.e];
                const det = a1*b2 - a2*b1;
                if (Math.abs(det) < 1e-30) continue;
                const C2 = (I1*b2 - I2*b1) / det;
                const Creg2 = (a1*I2 - a2*I1) / det;

                // Corner edge: singular part saturated below r₀ (skin-depth scale),
                // regular part integrates to C_reg²·L.
                // Ensure C2 >= 0 (physical singularity must be non-negative)
                const C2_eff = Math.max(C2, 0);
                const Creg2_eff = Math.max(Creg2, 0);
                h2dl_analytical = satCornerIntegral(C2_eff, u, L, satR0) + Creg2_eff * L;
            } else {
                // Fallback: single fitting edge, pure singularity model
                const f0 = fitEdges[0];
                const fitFactor = Math.pow(f0.r2, u) - Math.pow(f0.r1, u);
                if (Math.abs(fitFactor) < 1e-30) continue;
                const C2_fb = edgeH2dl[f0.e] * u / fitFactor;
                h2dl_analytical = satCornerIntegral(C2_fb, u, L, satR0);
            }

            const corr = h2dl_analytical - edgeH2dl[ce.e];
            totalCorr += corr;
        }
    }

    const h2dl_corrected = h2dl + totalCorr;
    if (Math.abs(totalCorr) > 1e-30 && h2dl > 0) {
        globalThis.__TRI_DEBUG__ && console.log(`    Corner correction: ${(totalCorr/h2dl*100).toFixed(1)}%, ` +
            `corrected h2dl=${h2dl_corrected.toExponential(3)} (${corners.length} corners)`);
    }
    return h2dl_corrected;
}

// --- Build isLossEdge array for microstrip geometry ---
// Marks conductor boundary edges + ground (y=0) edges.
// condRect: if provided and condRect.symmetry > 1, excludes edges on the
// symmetry plane (x = xmin_domain) — these are not physical conductor surfaces.
export function buildMicrostripLossEdges(mesh, fm, condRect) {
    const { nodes, edges, nEdges, nTris, triEdges } = mesh;
    const symX = (condRect && condRect.symmetry > 1) ? condRect.xmin_domain : null;
    const isLossEdge = new Uint8Array(nEdges);
    // On meshes with conductor interiors (meshConductorInterior), isCondEdge also
    // marks edges strictly inside the metal. A loss edge is a conductor SURFACE
    // edge: it must touch at least one exterior (non-conductor) triangle.
    const touchesExterior = new Uint8Array(nEdges);
    for (let t = 0; t < nTris; t++) {
        if (fm.faceF && fm.faceF[2*t] < 0) continue;
        for (let k = 0; k < 3; k++) touchesExterior[triEdges[3*t+k]] = 1;
    }
    for (let e = 0; e < nEdges; e++) {
        const n0=edges[2*e], n1=edges[2*e+1];
        // Conductor boundary edges (from freedom map)
        if (fm.isCondEdge && fm.isCondEdge[e] && touchesExterior[e]) {
            // Exclude symmetry plane edges — not physical conductor surfaces
            if (symX !== null && Math.abs(nodes[2*n0] - symX) < 1e-12 &&
                Math.abs(nodes[2*n1] - symX) < 1e-12) continue;
            isLossEdge[e] = 1; continue;
        }
        // Ground plane (y=0)
        if (Math.abs(nodes[2*n0+1]) < 1e-12 && Math.abs(nodes[2*n1+1]) < 1e-12) {
            if (!(fm.isCondEdge && fm.isCondEdge[e])) isLossEdge[e] = 1;
        }
    }
    return isLossEdge;
}

// --- Main entry ---
// isLossEdge: optional Uint8Array marking which edges contribute to conductor loss.
//   If not provided, uses buildMicrostripLossEdges (backward compatible).
// condArea: conductor cross-section area for DC resistance (default: computed from condRects).
export function solveConductorLoss(condRects, freq, sigma, extMesh, fm, vecRe, vecIm,
                                    gamma2Re, gamma2Im, P_poynting, Z0, epsMap, isLossEdge, enrichment, projectedH, wasmSolver) {
    const omega = 2*Math.PI*freq;
    const delta = Math.sqrt(2/(omega*MU0*sigma));
    const Rs = 1/(sigma*delta);
    const P = Math.abs(P_poynting);

    let gamma = csqrt(gamma2Re, gamma2Im);
    if(gamma.im<0) gamma={re:-gamma.re,im:-gamma.im};
    const omu = 2 * Math.PI * freq * MU0;

    // Build loss edge mask if not provided
    if (!isLossEdge) isLossEdge = buildMicrostripLossEdges(extMesh, fm, condRects && condRects[0]);

    // Galerkin-projected H with γ-scaled convention: both h2dl and P_proj carry the same
    // γ² factor from the e_t term, which cancels exactly in the ratio αc = h2dl/(4P).
    // The enrichment improves the RHS accuracy (enriched ∇Ez and Et), and the projection
    // distributes the singularity information over the boundary through the mass matrix.
    const projResult = projectedH || projectH(extMesh, fm, vecRe, vecIm, gamma, freq, enrichment, wasmSolver);
    const {htRe,htIm,hzRe,hzIm,hDofs} = projResult;
    // Use corner-corrected h2dl: analytically corrects corner-adjacent edges using a
    // two-parameter model (singular + regular) fitted from neighboring edges.
    // Derive hSub and epsR from condRects/epsMap for dielectric-aware corner correction.
    const condRect0 = condRects && condRects[0];
    const hSub_loss = condRect0 ? condRect0.ymin : undefined;
    let epsR_loss = 1;
    if (epsMap) for (const e of epsMap) if (e.re > epsR_loss) epsR_loss = e.re;
    // Saturate the corner singular model below r₀ = SAT_KAPPA·δ: the physical field
    // stops growing at skin-depth scale, while the PEC model grows as r^{ν-1} to r=0.
    // Without saturation the corrected integral converges to the PEC perturbation
    // limit, which overestimates the true (SIBC) loss. κ calibrated against the
    // QEP SIBC eigenvalue (validated to −0.03% on analytic circular coax).
    const h2dl = evalH2dlCorrected(extMesh, fm, htRe, htIm, hzRe, hzIm, omu, isLossEdge,
        hDofs, condRect0, hSub_loss, epsR_loss, SAT_KAPPA * delta);
    const Pproj = computePoyntingFromProjectedH(extMesh, fm, vecRe, vecIm, htRe, htIm, omu, hDofs);
    const Peff = Math.abs(Pproj) > 1e-30 ? Math.abs(Pproj) : P;

    const alpha_c_ac = Peff>1e-30 ? Rs*h2dl/(4*Peff) : 0;

    // For half-domain, signal_area from condRects is half the physical conductor area
    const sym = (condRects && condRects[0] && condRects[0].symmetry) || 1;
    let signal_area=0;
    if (condRects) for(const cr of condRects) signal_area+=Math.abs((cr.xmax-cr.xmin)*(cr.ymax-cr.ymin));
    signal_area *= sym;
    const R_dc = signal_area>0 ? 1/(sigma*signal_area) : 0;
    const R_ac = 2*alpha_c_ac*Z0;
    const R_combined = Math.sqrt(R_dc*R_dc+R_ac*R_ac);
    const alpha_c = Z0>0 ? R_combined/(2*Z0) : 0;
    const L_int = Peff>1e-30 ? Rs*h2dl/(4*Peff*omega) : 0;

    globalThis.__TRI_DEBUG__ && console.log(`  Conductor loss: Rs=${Rs.toFixed(4)} Ohm/sq, delta=${(delta*1e6).toFixed(2)} um, P=${P.toExponential(3)}`);
    globalThis.__TRI_DEBUG__ && console.log(`    h2dl=${h2dl.toExponential(3)}, alpha_c_ac=${(alpha_c_ac*8.686).toFixed(4)} dB/m`);
    globalThis.__TRI_DEBUG__ && console.log(`    R_ac=${R_ac.toFixed(4)} Ohm/m, R_dc=${R_dc.toFixed(4)} Ohm/m, R=${R_combined.toFixed(4)} Ohm/m`);

    return {R_ac, R_dc, R_combined, L_int, alpha_c, alpha_c_dBm:alpha_c*8.686, Rs, delta, h2dl};
}

// --- Zienkiewicz-Zhu error estimator on projected Ht ---
// Computes per-element error from the discontinuity of curl_z(Ht) between elements.
// Algorithm: compute raw curl_z(Ht) per element at centroid, smooth to vertices via
// area-weighted averaging, then error = ||curl_raw - curl_smooth||² × Area per element.
// Returns Float64Array(nTris) of per-element error estimates.
export function computeHtZZMetric(mesh, fm, vecRe, vecIm, gamma2Re, gamma2Im, freq, condRects, enrichment, projectedH, wasmSolver) {
    const { nodes, tris, triEdges, triSigns, nTris, nEdges, nNodes } = mesh;
    const { elemOrder } = fm;
    const edgeVerts = [[0,1],[1,2],[2,0]];
    const TOL = 1e-12;

    let gamma = csqrt(gamma2Re, gamma2Im);
    if (gamma.im < 0) gamma = { re: -gamma.re, im: -gamma.im };

    const { htRe, htIm, omu, hDofs } = projectedH || projectH(mesh, fm, vecRe, vecIm, gamma, freq, enrichment || null, wasmSolver);

    const curlRe = new Float64Array(nTris);
    const curlIm = new Float64Array(nTris);
    const triArea = new Float64Array(nTris);

    for (let t = 0; t < nTris; t++) {
        const v0 = tris[3*t], v1 = tris[3*t+1], v2 = tris[3*t+2];
        const xc = (nodes[2*v0]+nodes[2*v1]+nodes[2*v2])/3;
        const yc = (nodes[2*v0+1]+nodes[2*v1+1]+nodes[2*v2+1])/3;

        let inCond = false;
        if (condRects) for (const cr of condRects)
            if (xc>=cr.xmin-TOL && xc<=cr.xmax+TOL && yc>=cr.ymin-TOL && yc<=cr.ymax+TOL) { inCond=true; break; }
        if (inCond) continue;

        const { coeff, Area } = triCoefficients(nodes, v0, v1, v2);
        triArea[t] = Area;

        // Use the projected H DOF count (from gatherHt) for curl evaluation,
        // not elemOrder — the projected H lives in the full space and P2 elements
        // adjacent to P3 elements have ne3 DOFs from the projection.
        const {hR,hI,nNed} = gatherHt(t, mesh, htRe, htIm, hDofs);
        const evalOrder = nNed >= 15 ? 3 : 2;
        const {curlW} = evalNedelecBasis(coeff, evalOrder, edgeVerts, xc, yc);

        let cR = 0, cI = 0;
        for (let k = 0; k < nNed; k++) { cR += curlW[k]*hR[k]; cI += curlW[k]*hI[k]; }

        curlRe[t] = cI / omu;
        curlIm[t] = -cR / omu;
    }

    // Step 3: smooth curl to vertices (area-weighted average of adjacent element curls)
    const vertCurlRe = new Float64Array(nNodes);
    const vertCurlIm = new Float64Array(nNodes);
    const vertWeight = new Float64Array(nNodes);

    for (let t = 0; t < nTris; t++) {
        if (triArea[t] <= 0) continue;
        const A = triArea[t];
        for (let k = 0; k < 3; k++) {
            const v = tris[3*t+k];
            vertCurlRe[v] += curlRe[t] * A;
            vertCurlIm[v] += curlIm[t] * A;
            vertWeight[v] += A;
        }
    }
    for (let v = 0; v < nNodes; v++) {
        if (vertWeight[v] > 0) {
            vertCurlRe[v] /= vertWeight[v];
            vertCurlIm[v] /= vertWeight[v];
        }
    }

    // Step 4: per-element error = ∫|curl_raw - curl_smooth|² dA
    // Approximate: evaluate smoothed curl at centroid (average of 3 vertex values),
    // then error = |curl_raw - curl_smooth|² × Area
    const metric = new Float64Array(nTris);
    for (let t = 0; t < nTris; t++) {
        if (triArea[t] <= 0) continue;
        const v0 = tris[3*t], v1 = tris[3*t+1], v2 = tris[3*t+2];
        const smoothRe = (vertCurlRe[v0] + vertCurlRe[v1] + vertCurlRe[v2]) / 3;
        const smoothIm = (vertCurlIm[v0] + vertCurlIm[v1] + vertCurlIm[v2]) / 3;
        const dRe = curlRe[t] - smoothRe;
        const dIm = curlIm[t] - smoothIm;
        metric[t] = (dRe*dRe + dIm*dIm) * triArea[t];
    }

    return metric;
}

// --- Static conductor loss from ∮|∇φ|²dl ---
// Uses the static potential φ (P2 Lagrange) to compute conductor loss directly.
// For TEM: Ht = ẑ×∇φ/Z₀, so ∮|Ht|²dl = ∮|∇φ|²dl/Z₀².
// α_c = Rs · ∮|∇φ|²dl / (2·Z₀)  (with φ normalized so V=1).
// For quasi-TEM at frequency f: scale by √(ε_eff(f)/ε_eff(static)) via Z₀(f).
export function staticConductorLoss(condRects, freq, sigma, mesh, fm, phi, Z0, epsEff, epsEffMode, enrichment = null, isLossEdgeArg = null) {
    const omega = 2*Math.PI*freq;
    const delta = Math.sqrt(2/(omega*MU0*sigma));
    const Rs = 1/(sigma*delta);

    const { nodes, tris, edges, triEdges, nEdges, nTris, nNodes } = mesh;
    const { faceF, nodeF, edgeNodeF, nFreeVertexDof } = fm;
    const edgeVerts = [[0,1],[1,2],[2,0]];

    // Build loss edge mask (caller may supply a generalized one)
    const isLossEdge = isLossEdgeArg || buildMicrostripLossEdges(mesh, fm, condRects && condRects[0]);

    const edgeToTri = new Array(nEdges);
    for (let e=0;e<nEdges;e++) edgeToTri[e]=[];
    for (let t=0;t<nTris;t++) for (let le=0;le<3;le++) edgeToTri[triEdges[3*t+le]].push(t);

    // Standard 3-point Gauss-Legendre on [0,1]
    const GL3p=[0.11270,0.50000,0.88730], GL3w=[0.27778,0.44444,0.27778];

    // Graded quadrature for edges touching PEC corners: use substitution s = t^β
    // to cluster points toward t=0 (the corner end). With β = 1/(2ν), this cancels
    // the r^{2(ν-1)} singularity in |∇φ|². Use 5-point Gauss on the s variable.
    // GL5 on [0,1]:
    const GL5p=[0.04691008,0.23076534,0.50000000,0.76923466,0.95308992];
    const GL5w=[0.11846345,0.23931434,0.28444444,0.23931434,0.11846345];

    // Identify corner node indices for graded quadrature
    const cornerNodes = new Map(); // nodeIdx → corner info
    if (enrichment && enrichment.corners) {
        for (const c of enrichment.corners) {
            if (c.nodeIdx >= 0) cornerNodes.set(c.nodeIdx, c);
        }
    } else if (condRects) {
        // Even without enrichment, detect corners for graded quadrature
        const TOL = 1e-10;
        for (const cr of condRects) {
            const symX = (cr.symmetry > 1) ? cr.xmin_domain : null;
            const pts = [
                { x: cr.xmin, y: cr.ymin }, { x: cr.xmax, y: cr.ymin },
                { x: cr.xmax, y: cr.ymax }, { x: cr.xmin, y: cr.ymax },
            ];
            for (const p of pts) {
                if (symX !== null && Math.abs(p.x - symX) < TOL) continue;
                for (let n = 0; n < nNodes; n++) {
                    if (Math.abs(nodes[2*n] - p.x) < TOL && Math.abs(nodes[2*n+1] - p.y) < TOL) {
                        cornerNodes.set(n, { nu: 2/3 }); // default ν for 90° corner
                        break;
                    }
                }
            }
        }
    }

    let gradPhi2dl = 0;

    for (let e=0;e<nEdges;e++){
        if (!isLossEdge[e]) continue;

        // Find adjacent exterior triangle
        let adj=-1;
        for(const t of edgeToTri[e]){
            if (faceF[2*t] >= 0) { adj=t; break; }
        }
        if(adj<0) continue;

        const v0=tris[3*adj],v1=tris[3*adj+1],v2=tris[3*adj+2];
        const {coeff}=triCoefficients(nodes,v0,v1,v2);
        const n0=edges[2*e],n1=edges[2*e+1];
        const x0=nodes[2*n0],y0=nodes[2*n0+1],tx=nodes[2*n1]-x0,ty=nodes[2*n1+1]-y0;
        const L=Math.sqrt(tx*tx+ty*ty);

        // Gather φ DOFs for the adjacent triangle
        const verts=[v0,v1,v2];
        const elemOrd = fm.elemOrder ? fm.elemOrder[adj] : 2;
        const nLocalZ = elemOrd >= 3 ? 10 : 6;
        const phiLoc = new Float64Array(nLocalZ);
        for (let k=0;k<3;k++) phiLoc[k] = phi.phiVertex[verts[k]];
        for (let k=0;k<3;k++) phiLoc[k+3] = phi.phiEdge[triEdges[3*adj+k]];
        if (elemOrd >= 3 && phi.phiEdge3) {
            const signs = mesh.triSigns;
            for (let k=0;k<3;k++) phiLoc[k+6] = signs[3*adj+k] * phi.phiEdge3[triEdges[3*adj+k]];
            phiLoc[9] = phi.phiFaceNode ? phi.phiFaceNode[adj] : 0;
        }

        // Check if this edge touches a PEC corner vertex
        const c0 = cornerNodes.get(n0), c1 = cornerNodes.get(n1);
        let qp, qw, nq;
        if (c0 || c1) {
            // Use graded quadrature: s = t^β where β = 1/(2ν) ≈ 0.75 for ν=2/3.
            // If corner is at n0: t goes from 0 (corner) to 1, s = t^β, dt = s^(1/β-1)/(β) ds.
            // If corner is at n1: reverse the direction.
            const nu = (c0 || c1).nu || 2/3;
            const beta = 1 / (2 * nu); // grading exponent
            nq = 5;
            qp = new Float64Array(nq);
            qw = new Float64Array(nq);
            for (let q = 0; q < nq; q++) {
                const s = GL5p[q]; // uniform in s ∈ [0,1]
                const t = Math.pow(s, beta); // graded: clusters toward t=0
                const dtds = beta * Math.pow(s, beta - 1); // dt/ds
                if (c0) {
                    qp[q] = t; // corner at n0, t=0 is at n0
                } else {
                    qp[q] = 1 - t; // corner at n1, t=0 is at n1
                }
                qw[q] = GL5w[q] * dtds;
            }
        } else {
            nq = 3;
            qp = GL3p;
            qw = GL3w;
        }

        for(let q=0;q<nq;q++){
            const px=x0+qp[q]*tx, py=y0+qp[q]*ty;

            // ∇φ at quadrature point (P2/P3 Lagrange gradient + enrichment)
            let gx=0, gy=0;
            for (let k=0;k<3;k++){
                const [dgx,dgy]=lvGrad(coeff,k,px,py);
                gx+=dgx*phiLoc[k]; gy+=dgy*phiLoc[k];
            }
            for (let k=0;k<3;k++){
                const [p,qq]=edgeVerts[k];
                const [dgx,dgy]=leGrad(coeff,p,qq,px,py);
                gx+=dgx*phiLoc[k+3]; gy+=dgy*phiLoc[k+3];
            }
            if (elemOrd >= 3) {
                for (let k=0;k<3;k++){
                    const [p,qq]=edgeVerts[k];
                    const [dgx,dgy]=le3Grad(coeff,p,qq,px,py);
                    gx+=dgx*phiLoc[k+6]; gy+=dgy*phiLoc[k+6];
                }
                { const [dgx,dgy]=lf3Grad(coeff,px,py); gx+=dgx*phiLoc[9]; gy+=dgy*phiLoc[9]; }
            }

            // Add enrichment gradient: ∇φ_enrich = Σ c_i · ∇ψ_i
            if (enrichment && enrichment.coeffs) {
                const { corners, evalFn, coeffs } = enrichment;
                for (let ci = 0; ci < corners.length; ci++) {
                    if (Math.abs(coeffs[ci]) < 1e-30) continue;
                    const { dphidx, dphidy } = evalFn(px, py, corners[ci]);
                    gx += coeffs[ci] * dphidx;
                    gy += coeffs[ci] * dphidy;
                }
            }

            gradPhi2dl += qw[q]*L*(gx*gx + gy*gy);
        }
    }

    // For half-domain symmetry, the line integral covers half the perimeter.
    // condRects[0].symmetry = 2 for half domain. Scale gradPhi2dl and signal_area.
    const sym = (condRects && condRects[0] && condRects[0].symmetry) || 1;
    gradPhi2dl *= sym;

    // Quasi-TEM: |Ht| = √ε_eff/η₀ · |∇φ_eps| where η₀ = √(μ₀/ε₀) ≈ 377Ω.
    // Uses the ε-weighted static potential φ_eps (not φ_air).
    // ∮|Ht|²dl = ε_eff/η₀² · ∮|∇φ_eps|²dl. Power: P = V²/(2Z₀).
    // α_c = Rs · ε_eff · Z₀ · ∮|∇φ_eps|²dl / (2η₀²)
    // Frequency correction: use ε_eff(mode) and Z₀(f).
    const EPS0 = 8.854187817e-12;
    const ETA0 = Math.sqrt(MU0 / EPS0); // ~376.73 Ω
    const ee = (epsEffMode && epsEff) ? epsEffMode : (epsEff || 1);
    const Z0f = (epsEffMode && epsEff) ? Z0 * Math.sqrt(epsEff / epsEffMode) : Z0;
    const alpha_c_static = Rs * ee * Z0f * gradPhi2dl / (2 * ETA0 * ETA0);

    let signal_area=0;
    if (condRects) for(const cr of condRects) signal_area+=Math.abs((cr.xmax-cr.xmin)*(cr.ymax-cr.ymin));
    signal_area *= sym;
    const R_dc = signal_area>0 ? 1/(sigma*signal_area) : 0;
    const R_ac = 2*alpha_c_static*Z0f;
    const R_combined = Math.sqrt(R_dc*R_dc+R_ac*R_ac);
    const alpha_c = Z0f>0 ? R_combined/(2*Z0f) : 0;

    globalThis.__TRI_DEBUG__ && console.log(`  Static loss: gradPhi2dl=${gradPhi2dl.toExponential(3)}, alpha_c_static=${(alpha_c_static*8.686).toFixed(4)} dB/m`);

    return { alpha_c, alpha_c_dBm: alpha_c*8.686, Rs, delta, gradPhi2dl, R_ac, R_dc, R_combined };
}

// --- QEP SIBC conductor loss via quadratic eigenvalue RQI ---
// Formulates the waveguide problem WITHOUT Lee-Jin substitution:
//   Q(γ)·[Et; Ez] = (A₀ + γ·A₁ + γ²·A₂)·[Et; Ez] = 0
// where A₀_zz = Sz - k²εMz ≠ 0 (no spectral singularity for SIBC).
// Adds SIBC terms to A₀, then refines PEC eigenvalue via QEP RQI.
export function qepSibcConductorLoss(mesh, fm_pec, vecRe, vecIm, gamma2Re, gamma2Im,
                                     freq, sigma_cond, epsMap, k2, condRect, abc, solveQepRQI,
                                     projectedH, wasmSolver, rqiIters = 3, isLossEdgeIn = null) {
    const { nodes, edges, tris, triEdges, triSigns, nEdges, nTris, nNodes } = mesh;
    const { edgeF, faceF, nodeF, edgeNodeF, nFreeTransverse, nFreeVertexDof,
            isCondEdge, isCondNode,
            edgeF3, faceF3, edgeNodeF3, faceNodeF, elemOrder, edgeOrder,
            nFreeEdgeNodeDof, nFreeEdgeNode3Dof, nFreeFaceNodeDof } = fm_pec;
    const hasP3 = !!elemOrder;
    const omega = 2 * Math.PI * freq;
    const delta_skin = Math.sqrt(2 / (omega * MU0 * sigma_cond));
    const Rs = 1 / (sigma_cond * delta_skin);
    const Zs_re = Rs, Zs_im = Rs;
    const k0 = Math.sqrt(k2);
    const TOL = 1e-10;

    // Build loss edge mask and SIBC freedom map (reuse from previous SIBC code)
    const isLossEdge = isLossEdgeIn || buildMicrostripLossEdges(mesh, fm_pec, condRect);
    const isCondBoundaryNode = new Uint8Array(nNodes);
    for (let e = 0; e < nEdges; e++) {
        if (!isLossEdge[e]) continue;
        isCondBoundaryNode[edges[2*e]] = 1;
        isCondBoundaryNode[edges[2*e+1]] = 1;
    }

    const xmin_d = condRect.xmin_domain, xmax_d = condRect.xmax_domain;
    const ymax_d = condRect.ymax_domain;
    const condRects = condRect.rects || [condRect];

    // SIBC edge PEC check: loss edges (signal + ground) → FREE with impedance BC.
    // Conductor-interior → PEC. Domain boundary → PEC.
    function isEdgePEC_s(e) {
        if (isLossEdge[e]) return false; // Signal boundary + ground ��� FREE
        if (isCondEdge[e]) return true;  // Conductor interior
        const n0=edges[2*e],n1=edges[2*e+1];
        const y0=nodes[2*n0+1],y1=nodes[2*n1+1],x0=nodes[2*n0],x1=nodes[2*n1];
        if (Math.abs(y0)<TOL && Math.abs(y1)<TOL) return false; // Ground → FREE
        if (!abc.left && Math.abs(x0-xmin_d)<TOL && Math.abs(x1-xmin_d)<TOL) return true;
        if (!abc.right && Math.abs(x0-xmax_d)<TOL && Math.abs(x1-xmax_d)<TOL) return true;
        if (!abc.top && Math.abs(y0-ymax_d)<TOL && Math.abs(y1-ymax_d)<TOL) return true;
        return false;
    }
    function isNodePEC_s(n) {
        if (isCondBoundaryNode[n]) return false; // Signal boundary node → FREE
        if (Math.abs(nodes[2*n+1])<TOL) return false; // Ground node → FREE
        if (isCondNode[n]) return true; // Conductor interior
        const x=nodes[2*n],y=nodes[2*n+1];
        if (!abc.left && Math.abs(x-xmin_d)<TOL) return true;
        if (!abc.right && Math.abs(x-xmax_d)<TOL) return true;
        if (!abc.top && Math.abs(y-ymax_d)<TOL) return true;
        return false;
    }

    // Compute SIBC edge order: use edgeOrder from PEC map if available,
    // otherwise max of adjacent element orders for boundary edges.
    let sEdgeOrder = null;
    if (hasP3) {
        sEdgeOrder = new Uint8Array(nEdges).fill(2);
        if (edgeOrder) {
            for (let e = 0; e < nEdges; e++) sEdgeOrder[e] = edgeOrder[e];
        }
        // For loss edges (boundary, 1 adjacent element), use that element's order
        for (let t = 0; t < nTris; t++) {
            const order = elemOrder[t];
            for (let k = 0; k < 3; k++) {
                const e = triEdges[3*t+k];
                if (isLossEdge[e] && order > sEdgeOrder[e]) sEdgeOrder[e] = order;
            }
        }
    }

    // Build SIBC DOF maps (P2 + P3)
    const sEdgeF = new Int32Array(2*nEdges).fill(-1);
    const sEdgeF3 = hasP3 ? new Int32Array(nEdges).fill(-1) : null;
    let sNEdge = 0;
    for (let e=0;e<nEdges;e++) {
        if (isEdgePEC_s(e)) continue;
        sEdgeF[2*e]=sNEdge++; sEdgeF[2*e+1]=sNEdge++;
        if (hasP3 && sEdgeOrder[e] >= 3) sEdgeF3[e] = sNEdge++;
    }
    const sFaceF = new Int32Array(2*nTris).fill(-1);
    const sFaceF3 = hasP3 ? new Int32Array(4*nTris).fill(-1) : null;
    let sNFace = 0;
    for (let t=0;t<nTris;t++) {
        if (faceF[2*t]<0) continue;
        sFaceF[2*t]=sNEdge+sNFace++; sFaceF[2*t+1]=sNEdge+sNFace++;
        if (hasP3 && elemOrder[t] >= 3) {
            for (let k=0;k<4;k++) sFaceF3[4*t+k]=sNEdge+sNFace++;
        }
    }
    const sNTrans = sNEdge + sNFace;
    const sNodeF = new Int32Array(nNodes).fill(-1);
    let sNVert = 0;
    for (let n=0;n<nNodes;n++) { if(isNodePEC_s(n)) continue; sNodeF[n]=sNVert++; }
    const sEdgeNodeF = new Int32Array(nEdges).fill(-1);
    let sNEdgeNode = 0;
    for (let e=0;e<nEdges;e++) {
        if (isEdgePEC_s(e)) continue;
        const ym=(nodes[2*edges[2*e]+1]+nodes[2*edges[2*e+1]+1])/2;
        const xm=(nodes[2*edges[2*e]]+nodes[2*edges[2*e+1]])/2;
        if (!abc.left && Math.abs(xm-xmin_d)<TOL) continue;
        if (!abc.right && Math.abs(xm-xmax_d)<TOL) continue;
        if (!abc.top && Math.abs(ym-ymax_d)<TOL) continue;
        sEdgeNodeF[e]=sNEdgeNode++;
    }
    const sEdgeNodeF3 = hasP3 ? new Int32Array(nEdges).fill(-1) : null;
    let sNEdgeNode3 = 0;
    if (hasP3) {
        for (let e=0;e<nEdges;e++) {
            if (sEdgeOrder[e] < 3) continue;
            if (sEdgeNodeF[e] < 0) continue;
            sEdgeNodeF3[e] = sNEdgeNode3++;
        }
    }
    const sFaceNodeF = hasP3 ? new Int32Array(nTris).fill(-1) : null;
    let sNFaceNode = 0;
    if (hasP3) {
        for (let t=0;t<nTris;t++) {
            if (elemOrder[t] < 3) continue;
            if (faceF[2*t] < 0) continue;
            sFaceNodeF[t] = sNFaceNode++;
        }
    }
    const sNLong = sNVert + sNEdgeNode + sNEdgeNode3 + sNFaceNode;
    const sN = sNTrans + sNLong;
    const sLzOff = sNTrans;
    const sLzMidOff = sLzOff + sNVert;
    const sLzEdge3Off = sLzMidOff + sNEdgeNode;
    const sLzFaceNodeOff = sLzEdge3Off + sNEdgeNode3;
    const N_pec = fm_pec.nFreeTransverse + fm_pec.nFreeLongitudinal;

    const nP3elem = hasP3 ? elemOrder.reduce((s,v) => s+(v>=3?1:0), 0) : 0;
    globalThis.__TRI_DEBUG__ && console.log(`  QEP-SIBC: N_pec=${N_pec}, N_sibc=${sN} (+${sN-N_pec} conductor DOFs), P3=${nP3elem}`);

    // --- Assemble A₀, A₁, A₂ matrices ---
    // A₀ = [[St-k²εMt, 0], [0, Sz-k²εMz]]  (element ARe for tt, element B_zz for zz)
    // A₁ = [[0, -G], [Gᵀ, 0]]  (gradient coupling with signs)
    // A₂ = [[-Mt, 0], [0, 0]]  (negative transverse mass)
    const a0R=[],a0C=[],a0VRe=[],a0VIm=[];
    const a1R=[],a1C=[],a1VRe=[],a1VIm=[];
    const a2R=[],a2C=[],a2VRe=[],a2VIm=[];

    for (let t=0;t<nTris;t++) {
        if (sFaceF[2*t] < 0 && faceF[2*t] < 0) continue; // Conductor interior

        const v0=tris[3*t],v1=tris[3*t+1],v2=tris[3*t+2];
        const eps = epsMap[t];
        const order = hasP3 ? elemOrder[t] : 2;

        let ARe, AIm, BRe, BIm, Dtt, nLocal, nTrans;
        if (order >= 3) {
            const m = computeTriP3Matrices(nodes, v0, v1, v2, eps.re, eps.im, k0);
            ARe=m.ARe; AIm=m.AIm; BRe=m.BRe; BIm=m.BIm; Dtt=m.Dtt;
            nLocal=25; nTrans=15;
        } else {
            const m = computeTriP2Matrices(nodes, v0, v1, v2, eps.re, eps.im, k0);
            ARe=m.ARe; AIm=m.AIm; BRe=m.BRe; BIm=m.BIm; Dtt=m.Dtt;
            nLocal=14; nTrans=8;
        }

        // Local-to-global DOF mapping (SIBC freedom map)
        const gd = new Int32Array(nLocal).fill(-1);
        const sg = new Float64Array(nLocal).fill(1);
        for (let le=0;le<3;le++) {
            const eIdx=triEdges[3*t+le], s=triSigns[3*t+le];
            gd[le]=sEdgeF[2*eIdx]; sg[le]=s;       // ne1 (antisymmetric)
            gd[le+4]=sEdgeF[2*eIdx+1];              // ne2 (symmetric)
            if (order >= 3) {
                // ne3 (antisymmetric like ne1)
                gd[le+8] = sEdgeF3 ? sEdgeF3[eIdx] : -1;
                sg[le+8] = s;
                // le midpoint at P3 index 18
                const enf=sEdgeNodeF[eIdx];
                gd[le+18]=enf>=0?sLzMidOff+enf:-1;
                // le3 (antisymmetric)
                const enf3 = sEdgeNodeF3 ? sEdgeNodeF3[eIdx] : -1;
                gd[le+21]=enf3>=0?sLzEdge3Off+enf3:-1;
                sg[le+21]=s;
            } else {
                // le midpoint at P2 index 11
                const enf=sEdgeNodeF[eIdx];
                gd[le+11]=enf>=0?sLzMidOff+enf:-1;
            }
        }
        gd[3]=sFaceF[2*t]; gd[7]=sFaceF[2*t+1];    // nf1,nf2
        if (order >= 3 && sFaceF3) {
            for (let k=0;k<4;k++) gd[11+k]=sFaceF3[4*t+k]; // nf3..nf6
        }
        const verts=[v0,v1,v2];
        const lvOff = order >= 3 ? 15 : 8;
        for (let k=0;k<3;k++) {
            const nf=sNodeF[verts[k]];
            gd[lvOff+k]=nf>=0?sLzOff+nf:-1;
        }
        if (order >= 3 && sFaceNodeF) {
            const fnf=sFaceNodeF[t];
            gd[24]=fnf>=0?sLzFaceNodeOff+fnf:-1;
        }

        // Assemble A₀, A₁, A₂ from element matrices
        for (let li=0;li<nLocal;li++) {
            const gi=gd[li]; if (gi<0) continue;
            const si=sg[li];
            for (let lj=0;lj<nLocal;lj++) {
                const gj=gd[lj]; if (gj<0) continue;
                const sij=si*sg[lj];
                const idx=li*nLocal+lj;

                // A₀: transverse from ARe/AIm (St-k²εMt), longitudinal from BRe/BIm (Sz-k²εMz)
                if (li<nTrans && lj<nTrans) {
                    const re=sij*ARe[idx], im=sij*AIm[idx];
                    if (re!==0||im!==0) { a0R.push(gi);a0C.push(gj);a0VRe.push(re);a0VIm.push(im); }
                } else if (li>=nTrans && lj>=nTrans) {
                    const re=sij*BRe[idx], im=sij*BIm[idx];
                    if (re!==0||im!==0) { a0R.push(gi);a0C.push(gj);a0VRe.push(re);a0VIm.push(im); }
                }

                // A₁: coupling [[0, -Gᵀ],[G, 0]]
                if (li<nTrans && lj>=nTrans) {
                    const re=-sij*BRe[idx], im=-sij*BIm[idx];
                    if (re!==0||im!==0) { a1R.push(gi);a1C.push(gj);a1VRe.push(re);a1VIm.push(im); }
                } else if (li>=nTrans && lj<nTrans) {
                    const re=sij*BRe[idx], im=sij*BIm[idx];
                    if (re!==0||im!==0) { a1R.push(gi);a1C.push(gj);a1VRe.push(re);a1VIm.push(im); }
                }

                // A₂: [[-Mt,0],[0,0]]
                if (li<nTrans && lj<nTrans) {
                    const re=-sij*Dtt[li*nTrans+lj];
                    if (re!==0) { a2R.push(gi);a2C.push(gj);a2VRe.push(re);a2VIm.push(0); }
                }
            }
        }
    }

    const nBaseA0 = a0R.length; // Track where base entries end

    // --- Add SIBC boundary terms to A₀ ---
    const omu = omega * MU0;
    const Zs2 = Zs_re*Zs_re + Zs_im*Zs_im;
    const C_re = omu * Zs_im / Zs2;
    const C_im = omu * Zs_re / Zs2;

    for (let e=0;e<nEdges;e++) {
        if (!isLossEdge[e]) continue;
        const c0=sEdgeF[2*e], c1=sEdgeF[2*e+1];
        if (c0<0) continue;
        const n0=edges[2*e],n1=edges[2*e+1];
        const L=Math.sqrt((nodes[2*n1]-nodes[2*n0])**2+(nodes[2*n1+1]-nodes[2*n0+1])**2);
        // Transverse tangential trace integrals on edge (p,q):
        // ne1_tang = 1/L, ne2_tang = (λ_p-λ_q)/L, ne3_tang = (λ_p-λ_q)²/L
        // ∫ne1²dl = 1/L, ∫ne2²dl = 1/(3L), ∫ne3²dl = 1/(5L)
        // Cross: ∫ne1·ne2 = 0 (odd), ∫ne2·ne3 = 0 (odd), ∫ne1·ne3 = 1/(3L)
        a0R.push(c0);a0C.push(c0);a0VRe.push(C_re/L);a0VIm.push(C_im/L);
        a0R.push(c1);a0C.push(c1);a0VRe.push(C_re/(3*L));a0VIm.push(C_im/(3*L));
        // P3 ne3 transverse impedance
        const c2 = sEdgeF3 ? sEdgeF3[e] : -1;
        if (c2 >= 0) {
            a0R.push(c2);a0C.push(c2);a0VRe.push(C_re/(5*L));a0VIm.push(C_im/(5*L));
            a0R.push(c0);a0C.push(c2);a0VRe.push(C_re/(3*L));a0VIm.push(C_im/(3*L));
            a0R.push(c2);a0C.push(c0);a0VRe.push(C_re/(3*L));a0VIm.push(C_im/(3*L));
        }
        // Longitudinal: P2 1D mass on conductor edges
        const cn0=sNodeF[n0]>=0?sLzOff+sNodeF[n0]:-1;
        const cn1=sNodeF[n1]>=0?sLzOff+sNodeF[n1]:-1;
        const ce=sEdgeNodeF[e]>=0?sLzMidOff+sEdgeNodeF[e]:-1;
        const ml=L/30;
        for (const [ri,ci,v] of [[cn0,cn0,4],[cn1,cn1,4],[ce,ce,16],[cn0,cn1,-1],[cn1,cn0,-1],[cn0,ce,2],[ce,cn0,2],[cn1,ce,2],[ce,cn1,2]]) {
            if (ri>=0&&ci>=0) { a0R.push(ri);a0C.push(ci);a0VRe.push(C_re*ml*v);a0VIm.push(C_im*ml*v); }
        }
        // P3 le3 longitudinal impedance: le3(p,q) = λ_p·λ_q·(λ_p-λ_q) = t(1-t)(1-2t) on edge
        // ∫le3²dl = L/210, ∫lv0·le3 = L/210, ∫lv1·le3 = -L/210, ∫le·le3 = 0 (odd)
        const ce3 = sEdgeNodeF3 ? sEdgeNodeF3[e] : -1;
        if (ce3 >= 0) {
            const cle3 = sLzEdge3Off + ce3;
            const ml3 = L/210;
            a0R.push(cle3);a0C.push(cle3);a0VRe.push(C_re*ml3);a0VIm.push(C_im*ml3);
            if (cn0>=0) { a0R.push(cn0);a0C.push(cle3);a0VRe.push(C_re*ml3);a0VIm.push(C_im*ml3);
                          a0R.push(cle3);a0C.push(cn0);a0VRe.push(C_re*ml3);a0VIm.push(C_im*ml3); }
            if (cn1>=0) { a0R.push(cn1);a0C.push(cle3);a0VRe.push(-C_re*ml3);a0VIm.push(-C_im*ml3);
                          a0R.push(cle3);a0C.push(cn1);a0VRe.push(-C_re*ml3);a0VIm.push(-C_im*ml3); }
        }
    }

    // --- No A₁ SIBC boundary correction needed ---
    // The longitudinal equation requires IBP of BOTH the Laplacian and the divergence:
    //   Laplacian IBP: -∮N(∂Ez/∂n)ds → gives +C∮NEz ds + γ∮N(n̂·Et)ds
    //   Divergence IBP: -γ∮N(n̂·Et)ds (from -γ∫N(∇·Et)dA → γ∫∇N·Et - γ∮N(n̂·Et)ds)
    // The γ-dependent terms CANCEL EXACTLY, leaving only +C∮NEz ds in A₀_zz.
    // Equivalently: the natural flux of the combined weak form is n̂·(∇Ez+γEt) = jωμ₀Hτ,
    // which the SIBC converts wholly to (jωμ₀/Zs)·Ez. Verified against the analytic
    // circular-coax SIBC eigenvalue (−0.03% error) and √ε_r scaling on uniform fill.

    const csr0 = tripletsToCSR(a0R,a0C,a0VRe,sN,a0VIm);
    const csr1 = tripletsToCSR(a1R,a1C,a1VRe,sN,a1VIm);
    const csr2 = tripletsToCSR(a2R,a2C,a2VRe,sN,a2VIm);

    // --- Map PEC eigenvector to physical fields [Et; Ez] ---
    // Lee-Jin: x = [e_t; ez] where e_t = γ·Et, ez = Ez
    // Physical: x_phys = [Et; Ez] = [e_t/γ; ez]
    let gamma0 = csqrt(gamma2Re, gamma2Im);
    if (gamma0.im < 0) gamma0 = { re: -gamma0.re, im: -gamma0.im };
    const gMag2 = gamma0.re*gamma0.re + gamma0.im*gamma0.im;
    // 1/γ = conj(γ)/|γ|²
    const invG_re = gamma0.re / gMag2, invG_im = -gamma0.im / gMag2;

    const initRe = new Float64Array(sN);
    const initIm = new Float64Array(sN);

    // Map transverse DOFs: Et = e_t / γ (complex division)
    const lzP = fm_pec.nFreeTransverse;
    for (let e=0;e<nEdges;e++) {
        // ne1
        const p0=edgeF[2*e], s0=sEdgeF[2*e];
        if (p0>=0 && s0>=0) {
            const etRe=vecRe[p0], etIm=vecIm[p0];
            initRe[s0] = etRe*invG_re - etIm*invG_im;
            initIm[s0] = etRe*invG_im + etIm*invG_re;
        }
        // ne2
        const p1=edgeF[2*e+1], s1=sEdgeF[2*e+1];
        if (p1>=0 && s1>=0) {
            const etRe=vecRe[p1], etIm=vecIm[p1];
            initRe[s1] = etRe*invG_re - etIm*invG_im;
            initIm[s1] = etRe*invG_im + etIm*invG_re;
        }
        // ne3 (P3)
        if (hasP3 && edgeF3 && sEdgeF3) {
            const p3=edgeF3[e], s3=sEdgeF3[e];
            if (p3>=0 && s3>=0) {
                const etRe=vecRe[p3], etIm=vecIm[p3];
                initRe[s3] = etRe*invG_re - etIm*invG_im;
                initIm[s3] = etRe*invG_im + etIm*invG_re;
            }
        }
    }
    // nf1, nf2
    for (let t=0;t<nTris;t++) {
        for (let k=0;k<2;k++) {
            const p=faceF[2*t+k], s=sFaceF[2*t+k];
            if (p>=0 && s>=0) {
                const etRe=vecRe[p], etIm=vecIm[p];
                initRe[s] = etRe*invG_re - etIm*invG_im;
                initIm[s] = etRe*invG_im + etIm*invG_re;
            }
        }
        // nf3-6 (P3)
        if (hasP3 && faceF3 && sFaceF3) {
            for (let k=0;k<4;k++) {
                const p=faceF3[4*t+k], s=sFaceF3[4*t+k];
                if (p>=0 && s>=0) {
                    const etRe=vecRe[p], etIm=vecIm[p];
                    initRe[s] = etRe*invG_re - etIm*invG_im;
                    initIm[s] = etRe*invG_im + etIm*invG_re;
                }
            }
        }
    }
    // Map longitudinal DOFs: Ez = ez (no scaling)
    // lv (vertex)
    for (let n=0;n<nNodes;n++) {
        const pn=nodeF[n], sn=sNodeF[n];
        if (pn>=0 && sn>=0) { initRe[sLzOff+sn]=vecRe[lzP+pn]; initIm[sLzOff+sn]=vecIm[lzP+pn]; }
    }
    // le (edge midpoint)
    const emP=lzP+nFreeVertexDof;
    for (let e=0;e<nEdges;e++) {
        const pe=edgeNodeF[e], se=sEdgeNodeF[e];
        if (pe>=0 && se>=0) { initRe[sLzMidOff+se]=vecRe[emP+pe]; initIm[sLzMidOff+se]=vecIm[emP+pe]; }
    }
    // le3 (P3 edge extra)
    if (hasP3 && edgeNodeF3 && sEdgeNodeF3) {
        const e3P = emP + nFreeEdgeNodeDof;
        for (let e=0;e<nEdges;e++) {
            const pe=edgeNodeF3[e], se=sEdgeNodeF3[e];
            if (pe>=0 && se>=0) { initRe[sLzEdge3Off+se]=vecRe[e3P+pe]; initIm[sLzEdge3Off+se]=vecIm[e3P+pe]; }
        }
    }
    // lf3 (P3 face bubble)
    if (hasP3 && faceNodeF && sFaceNodeF) {
        const fnP = emP + nFreeEdgeNodeDof + (nFreeEdgeNode3Dof || 0);
        for (let t=0;t<nTris;t++) {
            const pf=faceNodeF[t], sf=sFaceNodeF[t];
            if (pf>=0 && sf>=0) { initRe[sLzFaceNodeOff+sf]=vecRe[fnP+pf]; initIm[sLzFaceNodeOff+sf]=vecIm[fnP+pf]; }
        }
    }

    // --- Seed conductor boundary DOFs from projected H ---
    // The mapped PEC eigenvector has E=0 on the conductor boundary. Seeding those DOFs
    // with the SIBC relation E = Z_s·H gives RQI a better starting vector (the converged
    // eigenvalue is unaffected — verified identical with 3 vs 6 RQI iterations).
    // E_tang = Z_s·H_z (z-component of SIBC) and E_z = Z_s·H_t_tang (t-component of SIBC)
    // For the QEP physical field: Et = e_t/γ, so we set DOF = Z_s·H/(jωμ₀·γ)
    // (projectH returns jωμ₀·H, so divide by jωμ₀ to get physical H, then multiply by Z_s/γ)
    if (projectedH || wasmSolver) {
        if (!projectedH) {
            projectedH = projectH(mesh, fm_pec, vecRe, vecIm, gamma0, freq, null, wasmSolver);
        }
        const { htRe: pHtRe, htIm: pHtIm, hzRe: pHzRe, hzIm: pHzIm, hDofs: pHDofs } = projectedH;

        // Complex factor: Z_s / (jωμ₀ · γ) for transverse, Z_s / (jωμ₀) for longitudinal
        // Z_s/(jωμ₀) = Rs(1+j)/(jωμ₀) = Rs(1+j)(−j)/(ωμ₀) = Rs(−j+1)/(ωμ₀) = Rs(1−j)/(ωμ₀)
        const omu = omega * MU0;
        const ZomuRe = Zs_re / omu;  // Re(Zs/(jωμ₀)) ... actually need full complex division
        // Zs/(jωμ₀) = (Zs_re + j·Zs_im) / (j·ωμ₀) = (Zs_re + j·Zs_im)·(-j)/(ωμ₀)
        //            = (Zs_im - j·Zs_re) / (ωμ₀)
        const ZjmuRe = Zs_im / omu;
        const ZjmuIm = -Zs_re / omu;
        // Zs/(jωμ₀·γ) = ZjmuRe/γ complex division
        const ZgRe = ZjmuRe * invG_re - ZjmuIm * invG_im;
        const ZgIm = ZjmuRe * invG_im + ZjmuIm * invG_re;

        // Seed transverse conductor boundary DOFs: Et = Z_s·Hz/(jωμ₀·γ)
        // (SIBC z-component: Et_τ = Zs·Hz, but here we use Hz projected DOFs)
        // For each loss edge with SIBC DOF but NO PEC DOF:
        for (let e = 0; e < nEdges; e++) {
            if (!isLossEdge[e]) continue;
            // ne1 DOF on conductor boundary
            const s0 = sEdgeF[2*e], p0 = edgeF[2*e];
            if (s0 >= 0 && p0 < 0) {
                // This DOF is freed by SIBC but was PEC → seed from Hz
                // Hz DOF for this edge's vertices: average
                const n0 = edges[2*e], n1 = edges[2*e+1];
                const hzR = (pHzRe[n0] + pHzRe[n1]) / 2;
                const hzI = (pHzIm[n0] + pHzIm[n1]) / 2;
                initRe[s0] = ZgRe * hzR - ZgIm * hzI;
                initIm[s0] = ZgRe * hzI + ZgIm * hzR;
            }
            // ne2 DOF
            const s1 = sEdgeF[2*e+1], p1 = edgeF[2*e+1];
            if (s1 >= 0 && p1 < 0) {
                const n0 = edges[2*e], n1 = edges[2*e+1];
                const hzR = (pHzRe[n0] + pHzRe[n1]) / 2;
                const hzI = (pHzIm[n0] + pHzIm[n1]) / 2;
                // ne2 is the (λp-λq) moment — scale by 0 at midpoint
                initRe[s1] = 0; initIm[s1] = 0;
            }
        }
        // Seed longitudinal conductor boundary DOFs: Ez = Z_s·Ht_τ/(jωμ₀)
        // For conductor boundary nodes, use the projected Ht evaluated at the node
        // Approximate: use the average of adjacent conductor edge Ht DOFs
        for (let n = 0; n < nNodes; n++) {
            const sn = sNodeF[n], pn = nodeF[n];
            if (sn >= 0 && pn < 0) {
                // Freed conductor boundary node — seed from Ht tangential
                // Find adjacent loss edges and average their Ht ne1 DOFs
                let htSumR = 0, htSumI = 0, cnt = 0;
                for (let e = 0; e < nEdges; e++) {
                    if (!isLossEdge[e]) continue;
                    if (edges[2*e] !== n && edges[2*e+1] !== n) continue;
                    const base = pHDofs.htEdgeDofStart[e];
                    htSumR += pHtRe[base]; htSumI += pHtIm[base]; cnt++;
                }
                if (cnt > 0) {
                    htSumR /= cnt; htSumI /= cnt;
                    // Ez = Z_s·Ht_τ, no γ division for longitudinal (ez = Ez directly)
                    initRe[sLzOff + sn] = ZjmuRe * htSumR - ZjmuIm * htSumI;
                    initIm[sLzOff + sn] = ZjmuRe * htSumI + ZjmuIm * htSumR;
                }
            }
        }
        // Edge midpoint longitudinal DOFs
        for (let e = 0; e < nEdges; e++) {
            const se = sEdgeNodeF[e], pe = edgeNodeF[e];
            if (se >= 0 && pe < 0 && isLossEdge[e]) {
                const base = pHDofs.htEdgeDofStart[e];
                const htR = pHtRe[base], htI = pHtIm[base];
                initRe[sLzMidOff + se] = ZjmuRe * htR - ZjmuIm * htI;
                initIm[sLzMidOff + se] = ZjmuRe * htI + ZjmuIm * htR;
            }
        }
    }

    // --- QEP RQI ---
    // Shift = γ₀ (NOT γ₀²)
    const sigma = [gamma0.re, gamma0.im];

    // Build A₀ with very large impedance (≈PEC) for baseline.
    // PEC limit: Z_s → 0, C = jωμ/Z_s → ∞. Use C_pec = 1000 × C_real.
    const pecScale = 1000;
    const pecR=[...a0R.slice(0,nBaseA0)], pecC=[...a0C.slice(0,nBaseA0)];
    const pecVRe=[...a0VRe.slice(0,nBaseA0)], pecVIm=[...a0VIm.slice(0,nBaseA0)];
    for (let i = nBaseA0; i < a0R.length; i++) {
        pecR.push(a0R[i]); pecC.push(a0C[i]);
        pecVRe.push(a0VRe[i]*pecScale); pecVIm.push(a0VIm[i]*pecScale);
    }
    const csr0_pec = tripletsToCSR(pecR, pecC, pecVRe, sN, pecVIm);

    let pecResult;
    try {
        pecResult = solveQepRQI(sN, csr0_pec, csr1, csr2, sigma, initRe, initIm, rqiIters);
        globalThis.__TRI_DEBUG__ && console.log(`    QEP PEC-limit: γ=(${pecResult.evalRe.toExponential(4)}, ${pecResult.evalIm.toExponential(4)}) [Lee-Jin PEC: (${gamma0.re.toExponential(4)}, ${gamma0.im.toExponential(4)})]`);
    } catch(e) {
        globalThis.__TRI_DEBUG__ && console.log(`    QEP PEC-limit failed: ${e.message}`);
    }

    // Verify: compute QEP residual |Q(γ₀)·x₀| for the PEC eigenvalue
    // This tests whether the QEP matrices are correctly derived from Lee-Jin
    function csrMatVec(csr, xRe, xIm, N) {
        const yRe = new Float64Array(N), yIm = new Float64Array(N);
        for (let i = 0; i < N; i++) {
            let sr=0,si=0;
            for (let j = csr.rowPtr[i]; j < csr.rowPtr[i+1]; j++) {
                const c=csr.colIdx[j], vre=csr.valRe[j], vim=csr.valIm?csr.valIm[j]:0;
                sr += vre*xRe[c] - vim*xIm[c];
                si += vre*xIm[c] + vim*xRe[c];
            }
            yRe[i]=sr; yIm[i]=si;
        }
        return {re:yRe, im:yIm};
    }
    // Q(γ₀) = A₀ + γ₀·A₁ + γ₀²·A₂
    const A0x = csrMatVec(csr0_pec, initRe, initIm, sN);
    const A1x = csrMatVec(csr1, initRe, initIm, sN);
    const A2x = csrMatVec(csr2, initRe, initIm, sN);
    const g0r=gamma0.re, g0i=gamma0.im;
    const g02r=g0r*g0r-g0i*g0i, g02i=2*g0r*g0i;
    let resInt=0,xInt=0,resCond=0;
    for (let i=0;i<sN;i++) {
        const qr = A0x.re[i] + (g0r*A1x.re[i]-g0i*A1x.im[i]) + (g02r*A2x.re[i]-g02i*A2x.im[i]);
        const qi = A0x.im[i] + (g0r*A1x.im[i]+g0i*A1x.re[i]) + (g02r*A2x.im[i]+g02i*A2x.re[i]);
        const r2 = qr*qr + qi*qi;
        const x2 = initRe[i]*initRe[i] + initIm[i]*initIm[i];
        if (x2 > 1e-30) { resInt += r2; xInt += x2; }
        else resCond += r2;
    }
    globalThis.__TRI_DEBUG__ && console.log(`    QEP residual: interior=${(Math.sqrt(resInt/(xInt+1e-30))).toExponential(3)}, conductor=${Math.sqrt(resCond).toExponential(3)}`);

    let result;
    try {
        // QEP with SIBC
        result = solveQepRQI(sN, csr0, csr1, csr2, sigma, initRe, initIm, rqiIters);
    } catch (e) {
        globalThis.__TRI_DEBUG__ && console.log(`  QEP-SIBC RQI failed: ${e.message}`);
        return { alpha_c: NaN, alpha_c_dBm: NaN, Rs, delta: delta_skin };
    }

    // Conductor loss = α_SIBC - α_PEC_QEP (both from same QEP DOF set)
    const alpha_total = result.evalRe;
    const alpha_pec_qep = pecResult ? pecResult.evalRe : gamma0.re;
    const alpha_c = alpha_total - alpha_pec_qep;

    globalThis.__TRI_DEBUG__ && console.log(`  QEP-SIBC: γ_sibc=(${result.evalRe.toExponential(4)}, ${result.evalIm.toExponential(4)})`);
    globalThis.__TRI_DEBUG__ && console.log(`    α_total=${(alpha_total*8.686).toFixed(3)} dB/m, α_diel=${(alpha_pec_qep*8.686).toFixed(3)} dB/m, α_cond=${(alpha_c*8.686).toFixed(3)} dB/m`);

    return {
        alpha_c,
        alpha_c_dBm: alpha_c * 8.686,
        gamma_re: result.evalRe, gamma_im: result.evalIm,
        Rs, delta: delta_skin, N_sibc: sN, N_pec,
    };
}

// --- H-field QEP SIBC conductor loss ---
// Uses H as primary variable: curl((1/ε)curl(H)) - k₀²μ₀H = 0
// Key advantage: H_tang ≠ 0 on PEC → SIBC perturbation is first-order in Z_s
// → correct ε_r scaling (unlike E-field QEP which is second-order).
// SIBC coefficient: Z_s (not jωμ₀/Z_s as in E-field).
export function qepSibcConductorLossH(mesh, fm_pec, vecRe, vecIm, gamma2Re, gamma2Im,
                                       freq, sigma_cond, epsMap, k2, condRect, abc,
                                       solveQepRQI, projectedH, wasmSolver) {
    const { nodes, edges, tris, triEdges, triSigns, nEdges, nTris, nNodes } = mesh;
    const { isCondEdge, isCondNode } = fm_pec;
    const omega = 2 * Math.PI * freq;
    const delta_skin = Math.sqrt(2 / (omega * MU0 * sigma_cond));
    const Rs = 1 / (sigma_cond * delta_skin);
    const Zs_re = Rs, Zs_im = Rs; // Z_s = (1+j)Rs
    const k0 = Math.sqrt(k2);
    const TOL = 1e-10;

    // Get projected H (initial vector source)
    let gamma0 = csqrt(gamma2Re, gamma2Im);
    if (gamma0.im < 0) gamma0 = { re: -gamma0.re, im: -gamma0.im };

    if (!projectedH) {
        projectedH = projectH(mesh, fm_pec, vecRe, vecIm, gamma0, freq, null, wasmSolver);
    }
    const { htRe, htIm, hzRe, hzIm, omu, hDofs } = projectedH;
    const { NH, NHz, htEdgeDofStart, htFaceDofStart, edgeMaxOrder } = hDofs;

    const isLossEdge = buildMicrostripLossEdges(mesh, fm_pec, condRect);
    const xmin_d = condRect.xmin_domain, xmax_d = condRect.xmax_domain;
    const ymax_d = condRect.ymax_domain;

    // --- Build H-field SIBC freedom map ---
    // Transverse Ht: ALL edges free (H_tang ≠ 0 on PEC). Use projectH's full numbering.
    // Longitudinal Hz: Eliminate on domain PEC walls (ground, side walls). Keep on conductor.
    // Transverse DOFs: same numbering as projectH (htEdgeDofStart, htFaceDofStart)
    const hNTrans = NH; // all transverse DOFs free

    // Longitudinal Hz: eliminate domain wall nodes
    function isHzEliminated(n) {
        // Ground plane (y=0): Hz = 0 (PEC wall)
        if (Math.abs(nodes[2*n+1]) < TOL) return true;
        // Domain walls (if not ABC)
        const x = nodes[2*n], y = nodes[2*n+1];
        if (!abc.left && abc.left !== 'pmc' && Math.abs(x - xmin_d) < TOL) return true;
        if (!abc.right && Math.abs(x - xmax_d) < TOL) return true;
        if (!abc.top && Math.abs(y - ymax_d) < TOL) return true;
        // PMC symmetry at x=0: Hz has Neumann BC (natural), keep free
        if (abc.left === 'pmc' && Math.abs(x - xmin_d) < TOL) return false;
        return false;
    }
    function isHzEdgeEliminated(e) {
        const n0 = edges[2*e], n1 = edges[2*e+1];
        const ym = (nodes[2*n0+1] + nodes[2*n1+1]) / 2;
        const xm = (nodes[2*n0] + nodes[2*n1]) / 2;
        if (Math.abs(ym) < TOL) return true; // ground
        if (!abc.left && abc.left !== 'pmc' && Math.abs(xm - xmin_d) < TOL) return true;
        if (!abc.right && Math.abs(xm - xmax_d) < TOL) return true;
        if (!abc.top && Math.abs(ym - ymax_d) < TOL) return true;
        return false;
    }

    // Hz vertex DOFs
    const hNodeF = new Int32Array(nNodes).fill(-1);
    let hNVert = 0;
    for (let n = 0; n < nNodes; n++) {
        if (isHzEliminated(n)) continue;
        hNodeF[n] = hNVert++;
    }
    // Hz edge midpoint DOFs
    const hEdgeNodeF = new Int32Array(nEdges).fill(-1);
    let hNEdgeNode = 0;
    for (let e = 0; e < nEdges; e++) {
        if (isHzEdgeEliminated(e)) continue;
        hEdgeNodeF[e] = hNEdgeNode++;
    }
    const hNLong = hNVert + hNEdgeNode;
    const hN = hNTrans + hNLong;
    const hLzOff = hNTrans;
    const hLzMidOff = hLzOff + hNVert;

    globalThis.__TRI_DEBUG__ && console.log(`  QEP-SIBC-H: N_h=${hN} (Ht=${hNTrans}, Hz=${hNLong})`);

    // --- Assemble H-field QEP matrices (Form 1: same as E-field!) ---
    // For non-magnetic media: curl(curl(H)) = k₀²ε_r·H (same structure as E-field).
    // A₀ = [[Att - k₀²ε·Dtt, 0], [0, Dzz1 - k₀²ε·Dzz2]]  (= ARe/AIm from E-field)
    // A₁ = [[0, -G^T], [G, 0]]  (= BRe off-diagonal, same as E-field)
    // A₂ = [[-Dtt, 0], [0, 0]]  (same as E-field)
    // Only freedom map, SIBC coefficient, and initial vector differ.
    const a0R=[],a0C=[],a0VRe=[],a0VIm=[];
    const a1R=[],a1C=[],a1VRe=[],a1VIm=[];
    const a2R=[],a2C=[],a2VRe=[],a2VIm=[];

    const hasP3 = !!fm_pec.elemOrder;

    for (let t = 0; t < nTris; t++) {
        // Skip conductor-interior triangles
        const v0=tris[3*t], v1=tris[3*t+1], v2=tris[3*t+2];
        const yc = (nodes[2*v0+1]+nodes[2*v1+1]+nodes[2*v2+1])/3;
        const xc = (nodes[2*v0]+nodes[2*v1]+nodes[2*v2])/3;
        let inCond = false;
        const condRects = condRect.rects || [condRect];
        for (const cr of condRects) {
            if (xc >= cr.xmin-TOL && xc <= cr.xmax+TOL &&
                yc >= cr.ymin-TOL && yc <= cr.ymax+TOL) { inCond = true; break; }
        }
        if (inCond) continue;

        const eps = epsMap[t];
        const order = hasP3 ? fm_pec.elemOrder[t] : 2;

        // Use standard E-field element matrices (Form 1 — identical for H-field)
        let ARe, AIm, BRe, BIm, Dtt, nTrans, nLong;
        if (order >= 3) {
            const m = computeTriP3Matrices(nodes, v0, v1, v2, eps.re, eps.im, k0);
            ARe=m.ARe; AIm=m.AIm; BRe=m.BRe; BIm=m.BIm; Dtt=m.Dtt;
            nTrans=15; nLong=10;
        } else {
            const m = computeTriP2Matrices(nodes, v0, v1, v2, eps.re, eps.im, k0);
            ARe=m.ARe; AIm=m.AIm; BRe=m.BRe; BIm=m.BIm; Dtt=m.Dtt;
            nTrans=8; nLong=6;
        }
        const nLocal = nTrans + nLong;

        // Build local-to-global DOF mapping for H-field
        // Transverse: use projectH numbering (htEdgeDofStart, htFaceDofStart)
        const gd_t = new Int32Array(nTrans).fill(-1); // transverse local→global
        const sg_t = new Float64Array(nTrans).fill(1);
        for (let le = 0; le < 3; le++) {
            const eIdx = triEdges[3*t+le], s = triSigns[3*t+le];
            const base = htEdgeDofStart[eIdx];
            gd_t[le] = base;       sg_t[le] = s;  // ne1
            gd_t[le+4] = base + 1;                 // ne2
            if (order >= 3 && edgeMaxOrder[eIdx] >= 3) {
                gd_t[le+8] = base + 2; sg_t[le+8] = s; // ne3
            }
        }
        const fBase = NH - (htFaceDofStart[nTris-1] + (fm_pec.elemOrder ? (fm_pec.elemOrder[nTris-1]>=3?6:2) : 2));
        // Actually, face DOF global = htEdgeTotal + htFaceDofStart[t] + local_face_idx
        const htEdgeTotal = hDofs.htEdgeTotal;
        const fOff = htEdgeTotal + htFaceDofStart[t];
        gd_t[3] = fOff;     // nf1
        gd_t[7] = fOff + 1; // nf2
        if (order >= 3) {
            for (let k = 0; k < 4; k++) gd_t[11+k] = fOff + 2 + k; // nf3-6
        }

        // Longitudinal: Hz vertex + edge midpoint DOFs
        const gd_z = new Int32Array(nLong).fill(-1);
        const sg_z = new Float64Array(nLong).fill(1);
        const verts = [v0, v1, v2];
        for (let k = 0; k < 3; k++) {
            const nf = hNodeF[verts[k]];
            gd_z[k] = nf >= 0 ? hLzOff + nf : -1;
        }
        for (let le = 0; le < 3; le++) {
            const eIdx = triEdges[3*t+le];
            const enf = hEdgeNodeF[eIdx];
            gd_z[le+3] = enf >= 0 ? hLzMidOff + enf : -1;
            if (order >= 3) {
                sg_z[le+6] = triSigns[3*t+le]; // le3 antisymmetric
                // le3 and lf3 DOFs: skip for now (P2 only in initial implementation)
            }
        }

        // Assemble A₀, A₁, A₂ — SAME as E-field QEP (Form 1)
        for (let li = 0; li < nLocal; li++) {
            const gi = li < nTrans ? gd_t[li] : gd_z[li - nTrans];
            if (gi < 0) continue;
            const si = li < nTrans ? sg_t[li] : sg_z[li - nTrans];

            for (let lj = 0; lj < nLocal; lj++) {
                const gj = lj < nTrans ? gd_t[lj] : gd_z[lj - nTrans];
                if (gj < 0) continue;
                const sij = si * (lj < nTrans ? sg_t[lj] : sg_z[lj - nTrans]);
                const idx = li * nLocal + lj;

                // A₀: tt from ARe (St-k²εMt), zz from BRe (Sz-k²εMz)
                if (li < nTrans && lj < nTrans) {
                    const re=sij*ARe[idx], im=sij*AIm[idx];
                    if (re!==0||im!==0) { a0R.push(gi);a0C.push(gj);a0VRe.push(re);a0VIm.push(im); }
                } else if (li >= nTrans && lj >= nTrans) {
                    const re=sij*BRe[idx], im=sij*BIm[idx];
                    if (re!==0||im!==0) { a0R.push(gi);a0C.push(gj);a0VRe.push(re);a0VIm.push(im); }
                }

                // A₁: coupling [[0, -G^T], [G, 0]]
                if (li < nTrans && lj >= nTrans) {
                    const re=-sij*BRe[idx], im=-sij*BIm[idx];
                    if (re!==0||im!==0) { a1R.push(gi);a1C.push(gj);a1VRe.push(re);a1VIm.push(im); }
                } else if (li >= nTrans && lj < nTrans) {
                    const re=sij*BRe[idx], im=sij*BIm[idx];
                    if (re!==0||im!==0) { a1R.push(gi);a1C.push(gj);a1VRe.push(re);a1VIm.push(im); }
                }

                // A₂: [[-Dtt, 0], [0, 0]]
                if (li < nTrans && lj < nTrans) {
                    const re=-sij*Dtt[li*nTrans+lj];
                    if (re!==0) { a2R.push(gi);a2C.push(gj);a2VRe.push(re);a2VIm.push(0); }
                }
            }
        }
    }

    const nBaseA0 = a0R.length;

    // --- Add H-field SIBC boundary terms to A₀ ---
    // From curl-curl IBP: -(1/ε)∮curl_z(H)·V_τ ds = -jωZ_s∮H_τ·V_τ ds
    // Coefficient: C_H = -jωZ_s = ω(Z_s_im - j·Z_s_re) = ωR_s(1-j) for Z_s=R_s(1+j)
    // H_tang ≠ 0 on PEC → first-order perturbation → correct ε_r scaling
    // SIBC coefficient for H-field (Form 1):
    // From curl-curl IBP: -∮curl_z(H)·V_τ ds, with curl_z(H) = jωε₀ε_adj·Ez
    // and SIBC Ez = Z_s·H_τ: C_H = -jωε₀ε_adj·Z_s per conductor edge.
    // Since ε_adj varies per edge, compute per-edge below.
    const EPS0 = 8.854187817e-12;

    // Build edge-to-triangle adjacency for finding ε_adj per loss edge
    const edgeToTri = new Array(nEdges);
    for (let e = 0; e < nEdges; e++) edgeToTri[e] = [];
    for (let t = 0; t < nTris; t++)
        for (let le = 0; le < 3; le++) edgeToTri[triEdges[3*t+le]].push(t);

    for (let e = 0; e < nEdges; e++) {
        if (!isLossEdge[e]) continue;
        const base = htEdgeDofStart[e];
        const n0 = edges[2*e], n1 = edges[2*e+1];
        const L = Math.sqrt((nodes[2*n1]-nodes[2*n0])**2 + (nodes[2*n1+1]-nodes[2*n0+1])**2);

        // Find adjacent non-conductor element's ε for this edge
        let epsAdj = 1.0;
        for (const t of edgeToTri[e]) {
            const v0=tris[3*t], v1=tris[3*t+1], v2=tris[3*t+2];
            const yc=(nodes[2*v0+1]+nodes[2*v1+1]+nodes[2*v2+1])/3;
            const xc=(nodes[2*v0]+nodes[2*v1]+nodes[2*v2])/3;
            let inC = false;
            for (const cr of (condRect.rects || [condRect])) {
                if (xc>=cr.xmin-TOL&&xc<=cr.xmax+TOL&&yc>=cr.ymin-TOL&&yc<=cr.ymax+TOL) { inC=true; break; }
            }
            if (!inC) { epsAdj = epsMap[t].re; break; }
        }

        // C_H = -jωε₀ε_adj·Z_s = ωε₀ε_adj(Z_s_im - j·Z_s_re) = ωε₀ε_adj·Rs·(1-j)
        const CH_re = omega * EPS0 * epsAdj * Zs_im;
        const CH_im = -omega * EPS0 * epsAdj * Zs_re;

        // Transverse: C_H · ∫W_tang·V_tang dl
        const c0 = base, c1 = base + 1;
        a0R.push(c0);a0C.push(c0);a0VRe.push(CH_re/L);a0VIm.push(CH_im/L);
        a0R.push(c1);a0C.push(c1);a0VRe.push(CH_re/(3*L));a0VIm.push(CH_im/(3*L));

        // P3 ne3
        if (edgeMaxOrder[e] >= 3) {
            const c2 = base + 2;
            a0R.push(c2);a0C.push(c2);a0VRe.push(CH_re/(5*L));a0VIm.push(CH_im/(5*L));
            a0R.push(c0);a0C.push(c2);a0VRe.push(CH_re/(3*L));a0VIm.push(CH_im/(3*L));
            a0R.push(c2);a0C.push(c0);a0VRe.push(CH_re/(3*L));a0VIm.push(CH_im/(3*L));
        }

        // Longitudinal: C_H · ∫N·N dl (double-IBP cancellation gives same C_H)
        const cn0 = hNodeF[n0] >= 0 ? hLzOff + hNodeF[n0] : -1;
        const cn1 = hNodeF[n1] >= 0 ? hLzOff + hNodeF[n1] : -1;
        const ce = hEdgeNodeF[e] >= 0 ? hLzMidOff + hEdgeNodeF[e] : -1;
        const ml = L / 30;
        for (const [ri,ci,v] of [[cn0,cn0,4],[cn1,cn1,4],[ce,ce,16],[cn0,cn1,-1],[cn1,cn0,-1],[cn0,ce,2],[ce,cn0,2],[cn1,ce,2],[ce,cn1,2]]) {
            if (ri >= 0 && ci >= 0) { a0R.push(ri);a0C.push(ci);a0VRe.push(CH_re*ml*v);a0VIm.push(CH_im*ml*v); }
        }
    }

    const csr0 = tripletsToCSR(a0R, a0C, a0VRe, hN, a0VIm);
    const csr1 = tripletsToCSR(a1R, a1C, a1VRe, hN, a1VIm);
    const csr2 = tripletsToCSR(a2R, a2C, a2VRe, hN, a2VIm);

    // --- Map projected H to initial vector ---
    // projectH returns jωμ₀·H. Physical H = htRe/(jωμ₀). QEP uses Ht (not γHt).
    // For the QEP initial vector, we need Ht (physical H, not Lee-Jin scaled).
    // Since the QEP eigenvalue is γ (not γ²), and the substitution is h_t = γHt:
    //   Physical Ht = htRe_projected / (jωμ₀)
    // But the QEP initial vector should be [Ht; Hz] in physical units.
    // The factor 1/(jωμ₀) is a complex constant that scales all DOFs uniformly
    // → doesn't affect the eigenvector direction. We can use htRe/htIm directly
    //   (the RQI normalizes the vector).

    const initRe = new Float64Array(hN);
    const initIm = new Float64Array(hN);

    // Transverse Ht: direct copy from projectH (same DOF numbering)
    for (let i = 0; i < NH; i++) {
        initRe[i] = htRe[i];
        initIm[i] = htIm[i];
    }
    // Longitudinal Hz: map from projectH full numbering to freedom map
    for (let n = 0; n < nNodes; n++) {
        const hn = hNodeF[n];
        if (hn >= 0) { initRe[hLzOff + hn] = hzRe[n]; initIm[hLzOff + hn] = hzIm[n]; }
    }
    for (let e = 0; e < nEdges; e++) {
        const he = hEdgeNodeF[e];
        if (he >= 0) { initRe[hLzMidOff + he] = hzRe[nNodes + e]; initIm[hLzMidOff + he] = hzIm[nNodes + e]; }
    }

    // --- QEP RQI ---
    const sigma = [gamma0.re, gamma0.im];

    // PEC-limit baseline: for H-field, PEC is the NATURAL BC (Z_s → 0).
    // Use A₀ WITHOUT boundary terms (just the volume integrals).
    const csr0_pec = tripletsToCSR(a0R.slice(0,nBaseA0), a0C.slice(0,nBaseA0),
        a0VRe.slice(0,nBaseA0), hN, a0VIm.slice(0,nBaseA0));

    let pecResult;
    try {
        pecResult = solveQepRQI(hN, csr0_pec, csr1, csr2, sigma, initRe, initIm, 3);
        globalThis.__TRI_DEBUG__ && console.log(`    QEP-H PEC-limit: γ=(${pecResult.evalRe.toExponential(4)}, ${pecResult.evalIm.toExponential(4)}) [Lee-Jin PEC: (${gamma0.re.toExponential(4)}, ${gamma0.im.toExponential(4)})]`);
    } catch(e) {
        globalThis.__TRI_DEBUG__ && console.log(`    QEP-H PEC-limit failed: ${e.message}`);
    }

    let result;
    try {
        result = solveQepRQI(hN, csr0, csr1, csr2, sigma, initRe, initIm, 3);
    } catch (e) {
        globalThis.__TRI_DEBUG__ && console.log(`  QEP-SIBC-H RQI failed: ${e.message}`);
        return { alpha_c: NaN, alpha_c_dBm: NaN, Rs, delta: delta_skin };
    }

    const alpha_total = result.evalRe;
    const alpha_pec_qep = pecResult ? pecResult.evalRe : gamma0.re;
    const alpha_c = alpha_total - alpha_pec_qep;

    globalThis.__TRI_DEBUG__ && console.log(`  QEP-SIBC-H: γ_sibc=(${result.evalRe.toExponential(4)}, ${result.evalIm.toExponential(4)})`);
    globalThis.__TRI_DEBUG__ && console.log(`    α_total=${(alpha_total*8.686).toFixed(3)} dB/m, α_diel=${(alpha_pec_qep*8.686).toFixed(3)} dB/m, α_cond=${(alpha_c*8.686).toFixed(3)} dB/m`);

    return {
        alpha_c,
        alpha_c_dBm: alpha_c * 8.686,
        gamma_re: result.evalRe, gamma_im: result.evalIm,
        Rs, delta: delta_skin, N_sibc: hN,
    };
}

// --- First-order SIBC conductor loss from projected H ---
// Computes the conductor loss directly from the H-field boundary integral,
// using the SIBC boundary mass matrix and the QEP derivative normalization.
// No iterative solve needed — just inner products with projected H DOFs.
//
// δγ = x_H^H · δA₀ · x_H / (x_H^H · Q'(γ₀) · x_H)
// where Q'(γ) = A₁ + 2γA₂, δA₀ = SIBC boundary mass (coefficient Z_s).
//
// The numerator gives Z_s·∮|H_tang|²dl (first-order, correct ε scaling).
// The denominator gives the normalization from the waveguide power flow.
export function sibcConductorLossFirstOrder(mesh, fm_pec, vecRe, vecIm, gamma2Re, gamma2Im,
                                             freq, sigma_cond, epsMap, k2, condRect, abc,
                                             projectedH, wasmSolver) {
    const { nodes, edges, tris, triEdges, triSigns, nEdges, nTris, nNodes } = mesh;
    const { isCondEdge, isCondNode } = fm_pec;
    const omega = 2 * Math.PI * freq;
    const delta_skin = Math.sqrt(2 / (omega * MU0 * sigma_cond));
    const Rs = 1 / (sigma_cond * delta_skin);
    const Zs_re = Rs, Zs_im = Rs;

    let gamma0 = csqrt(gamma2Re, gamma2Im);
    if (gamma0.im < 0) gamma0 = { re: -gamma0.re, im: -gamma0.im };

    if (!projectedH) {
        projectedH = projectH(mesh, fm_pec, vecRe, vecIm, gamma0, freq, null, wasmSolver);
    }
    const { htRe, htIm, hzRe, hzIm, hDofs } = projectedH;
    const { htEdgeDofStart, htFaceDofStart, edgeMaxOrder } = hDofs;
    const NH = hDofs.NH;

    const isLossEdge = buildMicrostripLossEdges(mesh, fm_pec, condRect);

    // --- Numerator: x^H · δA₀ · x = Z_s · ∮|H_tang|² dl + Z_s · ∮|Hz|² dl ---
    // The projected H stores jωμ₀·H. The boundary integral is:
    //   ∮ |jωμ₀·H_tang|² / (ωμ₀)² dl = ∮|H_tang|² dl  (the (ωμ₀)² cancels in ratio)
    // So we compute: Z_s · Σ_{loss edges} [∫|Ht_tang|²dl + ∫|Hz|²dl] using projected H DOFs.
    //
    // Ht tangential trace integrals (same formulas as SIBC):
    //   ne1: ∫ne1_tang² dl = 1/L,  ne2: 1/(3L),  ne1×ne3: 1/(3L), ne3: 1/(5L)

    let numRe = 0, numIm = 0;

    for (let e = 0; e < nEdges; e++) {
        if (!isLossEdge[e]) continue;
        const n0 = edges[2*e], n1 = edges[2*e+1];
        const L = Math.sqrt((nodes[2*n1]-nodes[2*n0])**2 + (nodes[2*n1+1]-nodes[2*n0+1])**2);
        const base = htEdgeDofStart[e];

        // Ht tangential: ne1 DOF
        const h0r = htRe[base], h0i = htIm[base];
        const h0sq = h0r*h0r + h0i*h0i; // |ne1 DOF|²
        numRe += h0sq / L; // ∫ne1²dl = 1/L

        // ne2 DOF
        const h1r = htRe[base+1], h1i = htIm[base+1];
        const h1sq = h1r*h1r + h1i*h1i;
        numRe += h1sq / (3*L); // ∫ne2²dl = 1/(3L)

        // ne3 DOF (P3)
        if (edgeMaxOrder[e] >= 3) {
            const h2r = htRe[base+2], h2i = htIm[base+2];
            const h2sq = h2r*h2r + h2i*h2i;
            numRe += h2sq / (5*L);
            // ne1×ne3 cross: complex conjugate product
            numRe += 2*(h0r*h2r + h0i*h2i) / (3*L);
        }

        // Hz on this edge: P2 1D mass ∫|Hz|²dl
        // Hz DOFs: vertex + edge midpoint
        const hz0r = hzRe[n0], hz0i = hzIm[n0];
        const hz1r = hzRe[n1], hz1i = hzIm[n1];
        const hzMr = hzRe[nNodes + e], hzMi = hzIm[nNodes + e];
        const ml = L / 30;
        // Diagonal: 4·|v0|² + 4·|v1|² + 16·|mid|²
        numRe += ml * (4*(hz0r*hz0r+hz0i*hz0i) + 4*(hz1r*hz1r+hz1i*hz1i) + 16*(hzMr*hzMr+hzMi*hzMi));
        // Off-diagonal: -2·Re(v0·v1*) + 4·Re(v0·mid*) + 4·Re(v1·mid*)
        numRe += ml * (-2*(hz0r*hz1r+hz0i*hz1i) + 4*(hz0r*hzMr+hz0i*hzMi) + 4*(hz1r*hzMr+hz1i*hzMi));
    }

    // Multiply by Z_s (complex): num_total = (Zs_re + j·Zs_im) · numRe
    // (numRe is real because it's |x|² terms)
    const numTotRe = Zs_re * numRe;
    const numTotIm = Zs_im * numRe;

    // --- Denominator: use Poynting power P (already available from projected H) ---
    // The perturbation formula: α_c = Rs · ∮|H|² dl / (4P)
    // P is computed from the projected H and the E-field eigenvector.
    const omu = omega * MU0;
    const P = Math.abs(computePoyntingFromProjectedH(mesh, fm_pec, vecRe, vecIm,
        htRe, htIm, omu, hDofs));

    // α_c = Rs · (∮|H_phys|² dl) / (4P)
    // The projected H stores jωμ₀·H. So |H_phys|² = |h_projected|²/(ωμ₀)².
    // numRe = ∮|h_projected_tang|² dl (from DOF-based computation above)
    const h2dl_phys = numRe / (omu * omu);
    const alpha_c = Rs * h2dl_phys / (4 * P);
    globalThis.__TRI_DEBUG__ && console.log(`  SIBC-1st: ∮|H|²_dof=${h2dl_phys.toExponential(4)}, P=${P.toExponential(4)}, α_c=${(alpha_c*8.686).toFixed(4)} dB/m`);

    return {
        alpha_c,
        alpha_c_dBm: alpha_c * 8.686,
        Rs, delta: delta_skin,
    };
}
