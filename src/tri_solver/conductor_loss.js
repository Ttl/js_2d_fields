// Conductor loss via perturbation with Galerkin-projected H (OpenParEM2D approach).
// Projects Ht and Hz into FEM spaces using ALL DOFs (no PEC elimination for H).
// PEC eliminates tangential E, not tangential H.
// Loss integral: α_c = Rs · ∮(|Ht|² + |Hz|²) dl / (4P)

import { csqrt, tripletsToCSR, GL3p, GL3w } from './fem_core.js';
import { ne1, ne2, nf1, nf2,
         lv, le, lvGrad, leGrad, triCoefficients,
         ne1Curl, ne2Curl, nf1Curl, nf2Curl,
         QW, QL1, QL2, QL3, NQ,
         getLzOffsets } from './tri_fem.js';

// --- Helper: edge → list-of-adjacent-triangles adjacency map ---
function buildEdgeToTri(nEdges, nTris, triEdges) {
    const edgeToTri = new Array(nEdges);
    for (let e = 0; e < nEdges; e++) edgeToTri[e] = [];
    for (let t = 0; t < nTris; t++) for (let le = 0; le < 3; le++) edgeToTri[triEdges[3*t+le]].push(t);
    return edgeToTri;
}

// --- Helper: evaluate P2 Nedelec basis functions at a point ---
// Returns { Wx, Wy, curlW, nNed } with nNed = 8.
function evalNedelecBasis(coeff, edgeVerts, xq, yq) {
    const nNed = 8;
    const Wx = new Float64Array(nNed), Wy = new Float64Array(nNed), curlW = new Float64Array(nNed);
    for (let k = 0; k < 3; k++) {
        const [p, q] = edgeVerts[k];
        [Wx[k], Wy[k]] = ne1(coeff, p, q, xq, yq); curlW[k] = ne1Curl(coeff, p, q);
        [Wx[k+4], Wy[k+4]] = ne2(coeff, p, q, xq, yq); curlW[k+4] = ne2Curl(coeff, p, q, xq, yq);
    }
    [Wx[3], Wy[3]] = nf1(coeff, xq, yq); curlW[3] = nf1Curl(coeff, xq, yq);
    [Wx[7], Wy[7]] = nf2(coeff, xq, yq); curlW[7] = nf2Curl(coeff, xq, yq);
    return { Wx, Wy, curlW, nNed };
}

// --- Helper: evaluate P2 Lagrange basis functions at a point ---
// Returns { Nz, Gx, Gy, nLag } with nLag = 6.
function evalLagrangeBasis(coeff, edgeVerts, xq, yq) {
    const nLag = 6;
    const Nz = new Float64Array(nLag), Gx = new Float64Array(nLag), Gy = new Float64Array(nLag);
    for (let k = 0; k < 3; k++) {
        Nz[k] = lv(coeff, k, xq, yq);
        [Gx[k], Gy[k]] = lvGrad(coeff, k, xq, yq);
    }
    for (let k = 0; k < 3; k++) {
        const [p, q] = edgeVerts[k];
        Nz[k+3] = le(coeff, p, q, xq, yq);
        [Gx[k+3], Gy[k+3]] = leGrad(coeff, p, q, xq, yq);
    }
    return { Nz, Gx, Gy, nLag };
}

// --- Helper: compute full H DOF counts (P2) ---
// Returns { NH, NHz, htEdgeTotal, htEdgeDofStart, htFaceDofStart }
function computeHDofCounts(mesh, fm) {
    const { nEdges, nTris, nNodes } = mesh;

    // Ht DOFs: 2 per edge (ne1, ne2) + 2 per face (nf1, nf2).
    const htEdgeDofStart = new Int32Array(nEdges);
    let htEdgeTotal = 0;
    for (let e = 0; e < nEdges; e++) {
        htEdgeDofStart[e] = htEdgeTotal;
        htEdgeTotal += 2;
    }

    const htFaceDofStart = new Int32Array(nTris);
    let htFaceTotal = 0;
    for (let t = 0; t < nTris; t++) {
        htFaceDofStart[t] = htFaceTotal;
        htFaceTotal += 2;
    }
    const NH = htEdgeTotal + htFaceTotal;

    // Hz DOFs: vertex + edge midpoint.
    const NHz = nNodes + nEdges;

    return { NH, NHz, htEdgeTotal, htEdgeDofStart, htFaceDofStart };
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
export function projectH(mesh, fm, vecRe, vecIm, gamma, freq, wasmSolver) {
    const { nodes, tris, triEdges, triSigns, nTris, nEdges, nNodes } = mesh;
    const { edgeF, faceF, nodeF, edgeNodeF } = fm;
    const omega = 2 * Math.PI * freq;
    const omu = omega * MU0;
    const { lzOff, lzEdgeMidOff } = getLzOffsets(fm);
    const edgeVerts = [[0, 1], [1, 2], [2, 0]];

    // Full H DOF counts
    const hDofs = computeHDofCounts(mesh, fm);
    const { NH, NHz, htEdgeDofStart, htFaceDofStart } = hDofs;

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
        const nNed = 8;
        const nLag = 6;

        // Ht DOFs: full numbering
        const hD = new Int32Array(nNed), hS = new Float64Array(nNed).fill(1);
        for (let k = 0; k < 3; k++) {
            const eIdx = triEdges[3*t+k], base = htEdgeDofStart[eIdx];
            hD[k] = base; hS[k] = triSigns[3*t+k]; // ne1
            hD[k+4] = base + 1; // ne2
        }
        const faceBlockStart = hDofs.htEdgeTotal;
        hD[3] = faceBlockStart + htFaceDofStart[t]; // nf1
        hD[7] = faceBlockStart + htFaceDofStart[t] + 1; // nf2

        // Hz DOFs: full numbering
        const hzD = new Int32Array(nLag);
        for (let k = 0; k < 3; k++) hzD[k] = verts[k];
        for (let k = 0; k < 3; k++) hzD[k+3] = nNodes + triEdges[3*t+k];

        // E DOFs: free (PEC eliminated) — for building RHS
        const eD = new Int32Array(nNed), eS = new Float64Array(nNed).fill(1);
        for (let k = 0; k < 3; k++) {
            const eIdx = triEdges[3*t+k];
            eD[k] = edgeF[2*eIdx]; eS[k] = triSigns[3*t+k];
            eD[k+4] = edgeF[2*eIdx+1];
        }
        eD[3] = faceF[2*t]; eD[7] = faceF[2*t+1];

        // Gather Et from eigenvector
        const etR = new Float64Array(nNed), etI = new Float64Array(nNed);
        for (let k = 0; k < nNed; k++) { const g=eD[k]; if(g>=0){etR[k]=eS[k]*vecRe[g]; etI[k]=eS[k]*vecIm[g];} }

        // Gather Ez from eigenvector
        const ezR = new Float64Array(nLag), ezI = new Float64Array(nLag);
        for (let k = 0; k < 3; k++) { const nf=nodeF[verts[k]]; if(nf>=0){ezR[k]=vecRe[lzOff+nf]; ezI[k]=vecIm[lzOff+nf];} }
        for (let k = 0; k < 3; k++) { const enf=edgeNodeF[triEdges[3*t+k]]; if(enf>=0){ezR[k+3]=vecRe[lzEdgeMidOff+enf]; ezI[k+3]=vecIm[lzEdgeMidOff+enf];} }

        const MtEl = new Float64Array(nNed * nNed);
        const rR = new Float64Array(nNed), rI = new Float64Array(nNed);
        const MzEl = new Float64Array(nLag * nLag);
        const rzR = new Float64Array(nLag), rzI = new Float64Array(nLag);

        const nqp = NQ;
        const qw = QW;
        const ql1 = QL1;
        const ql2 = QL2;
        const ql3 = QL3;

        for (let q = 0; q < nqp; q++) {
            const w = qw[q];
            const xq = txs[0]*ql1[q]+txs[1]*ql2[q]+txs[2]*ql3[q];
            const yq = tys[0]*ql1[q]+tys[1]*ql2[q]+tys[2]*ql3[q];

            const { Wx, Wy, curlW } = evalNedelecBasis(coeff, edgeVerts, xq, yq);
            const { Nz, Gx: LGx, Gy: LGy } = evalLagrangeBasis(coeff, edgeVerts, xq, yq);

            // ∇Ez and Et interpolation from eigenvector
            let dxR=0,dxI=0,dyR=0,dyI=0;
            for (let k=0;k<nLag;k++){dxR+=LGx[k]*ezR[k];dxI+=LGx[k]*ezI[k];dyR+=LGy[k]*ezR[k];dyI+=LGy[k]*ezI[k];}

            let exR=0,exI=0,eyR=0,eyI=0;
            for (let k=0;k<nNed;k++){exR+=Wx[k]*etR[k];exI+=Wx[k]*etI[k];eyR+=Wy[k]*etR[k];eyI+=Wy[k]*etI[k];}

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

// --- Gather H DOFs (P2) for a triangle using full numbering ---
function gatherHt(triIdx, mesh, htRe, htIm, hDofs) {
    const { triEdges, triSigns, nEdges } = mesh;
    const nNed = 8;
    const hR = new Float64Array(nNed), hI = new Float64Array(nNed);
    if (hDofs) {
        for (let k = 0; k < 3; k++) {
            const eIdx = triEdges[3*triIdx+k], s = triSigns[3*triIdx+k];
            const base = hDofs.htEdgeDofStart[eIdx];
            hR[k] = s*htRe[base]; hI[k] = s*htIm[base];
            hR[k+4] = htRe[base+1]; hI[k+4] = htIm[base+1];
        }
        const fb = hDofs.htEdgeTotal + hDofs.htFaceDofStart[triIdx];
        hR[3]=htRe[fb]; hI[3]=htIm[fb];
        hR[7]=htRe[fb+1]; hI[7]=htIm[fb+1];
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
    const { edgeF, faceF, nFreeTransverse } = fm;
    const lzOff = nFreeTransverse;
    const edgeVerts = [[0,1],[1,2],[2,0]];
    let Pz = 0;

    for (let t = 0; t < nTris; t++) {
        if (faceF[2*t] < 0) continue;

        const v0=tris[3*t],v1=tris[3*t+1],v2=tris[3*t+2];
        const nNed = 8;
        const {coeff,Area}=triCoefficients(nodes,v0,v1,v2);
        const txs=[nodes[2*v0],nodes[2*v1],nodes[2*v2]];
        const tys=[nodes[2*v0+1],nodes[2*v1+1],nodes[2*v2+1]];

        // Gather Et from eigenvector
        const eD=new Int32Array(nNed), eS=new Float64Array(nNed).fill(1);
        for(let k=0;k<3;k++){
            const eIdx=triEdges[3*t+k]; eD[k]=edgeF[2*eIdx]; eS[k]=triSigns[3*t+k];
            eD[k+4]=edgeF[2*eIdx+1];
        }
        eD[3]=faceF[2*t]; eD[7]=faceF[2*t+1];
        const etR=new Float64Array(nNed),etI=new Float64Array(nNed);
        for(let k=0;k<nNed;k++){const g=eD[k]; if(g>=0){etR[k]=eS[k]*vecRe[g]; etI[k]=eS[k]*vecIm[g];}}

        // Gather Ht from projection
        const {hR,hI}=gatherHt(t, mesh, htRe, htIm, hDofs);

        const nqp = NQ;
        const qw = QW;
        const ql1 = QL1, ql2 = QL2, ql3 = QL3;

        for(let q=0;q<nqp;q++){
            const w=qw[q];
            const xq=txs[0]*ql1[q]+txs[1]*ql2[q]+txs[2]*ql3[q];
            const yq=tys[0]*ql1[q]+tys[1]*ql2[q]+tys[2]*ql3[q];
            const {Wx,Wy}=evalNedelecBasis(coeff, edgeVerts, xq, yq);

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
    const edgeToTri = buildEdgeToTri(nEdges, nTris, triEdges);

    // --- Step 1: Compute per-edge h2dl and |Ht|² at 3 quadrature points ---
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
        const nLag = 6;
        const v0=tris[3*adj],v1=tris[3*adj+1],v2=tris[3*adj+2];
        const {coeff}=triCoefficients(nodes,v0,v1,v2);
        const x0=nodes[2*n0],y0=nodes[2*n0+1],tx=nodes[2*n1]-x0,ty=nodes[2*n1+1]-y0;
        const L=Math.sqrt(tx*tx+ty*ty);

        const adjVerts=[v0,v1,v2];
        const hzR_loc=new Float64Array(nLag), hzI_loc=new Float64Array(nLag);
        for(let k=0;k<3;k++){hzR_loc[k]=hzRe[adjVerts[k]]; hzI_loc[k]=hzIm[adjVerts[k]];}
        for(let k=0;k<3;k++){hzR_loc[k+3]=hzRe[nNodes+triEdges[3*adj+k]]; hzI_loc[k+3]=hzIm[nNodes+triEdges[3*adj+k]];}

        const ht2q = [0,0,0];
        for(let q=0;q<3;q++){
            const px=x0+GL3p[q]*tx, py=y0+GL3p[q]*ty;
            const {Wx,Wy}=evalNedelecBasis(coeff, edgeVertsLocal, px, py);
            let hxR=0,hxI=0,hyR=0,hyI=0;
            for(let k=0;k<hNed;k++){hxR+=Wx[k]*hR[k];hxI+=Wx[k]*hI[k];hyR+=Wy[k]*hR[k];hyI+=Wy[k]*hI[k];}
            ht2q[q]=(hxR*hxR+hxI*hxI+hyR*hyR+hyI*hyI)/omu2;
            const {Nz}=evalLagrangeBasis(coeff, edgeVertsLocal, px, py);
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
                                    gamma2Re, gamma2Im, P_poynting, Z0, epsMap, isLossEdge, projectedH, wasmSolver) {
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
    // The projection distributes the singularity information over the boundary through
    // the mass matrix.
    const projResult = projectedH || projectH(extMesh, fm, vecRe, vecIm, gamma, freq, wasmSolver);
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
export function computeHtZZMetric(mesh, fm, vecRe, vecIm, gamma2Re, gamma2Im, freq, condRects, projectedH, wasmSolver) {
    const { nodes, tris, triEdges, triSigns, nTris, nEdges, nNodes } = mesh;
    const edgeVerts = [[0,1],[1,2],[2,0]];
    const TOL = 1e-12;

    let gamma = csqrt(gamma2Re, gamma2Im);
    if (gamma.im < 0) gamma = { re: -gamma.re, im: -gamma.im };

    const { htRe, htIm, omu, hDofs } = projectedH || projectH(mesh, fm, vecRe, vecIm, gamma, freq, wasmSolver);

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

        const {hR,hI,nNed} = gatherHt(t, mesh, htRe, htIm, hDofs);
        const {curlW} = evalNedelecBasis(coeff, edgeVerts, xc, yc);

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
export function staticConductorLoss(condRects, freq, sigma, mesh, fm, phi, Z0, epsEff, epsEffMode, isLossEdgeArg = null) {
    const omega = 2*Math.PI*freq;
    const delta = Math.sqrt(2/(omega*MU0*sigma));
    const Rs = 1/(sigma*delta);

    const { nodes, tris, edges, triEdges, nEdges, nTris, nNodes } = mesh;
    const { faceF, nodeF, edgeNodeF, nFreeVertexDof } = fm;
    const edgeVerts = [[0,1],[1,2],[2,0]];

    // Build loss edge mask (caller may supply a generalized one)
    const isLossEdge = isLossEdgeArg || buildMicrostripLossEdges(mesh, fm, condRects && condRects[0]);

    const edgeToTri = buildEdgeToTri(nEdges, nTris, triEdges);

    // Standard 3-point Gauss-Legendre on [0,1]

    // Graded quadrature for edges touching PEC corners: use substitution s = t^β
    // to cluster points toward t=0 (the corner end). With β = 1/(2ν), this cancels
    // the r^{2(ν-1)} singularity in |∇φ|². Use 5-point Gauss on the s variable.
    // GL5 on [0,1]:
    const GL5p=[0.04691008,0.23076534,0.50000000,0.76923466,0.95308992];
    const GL5w=[0.11846345,0.23931434,0.28444444,0.23931434,0.11846345];

    // Identify corner node indices for graded quadrature
    const cornerNodes = new Map(); // nodeIdx → corner info
    if (condRects) {
        // Detect conductor corners for graded quadrature
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
        const nLocalZ = 6;
        const phiLoc = new Float64Array(nLocalZ);
        for (let k=0;k<3;k++) phiLoc[k] = phi.phiVertex[verts[k]];
        for (let k=0;k<3;k++) phiLoc[k+3] = phi.phiEdge[triEdges[3*adj+k]];

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

            // ∇φ at quadrature point (P2 Lagrange gradient)
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

