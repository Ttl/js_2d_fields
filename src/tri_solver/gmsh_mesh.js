// gmsh WASM mesh generator for 2D transmission line cross-sections.
// Builds geometry from domain bounds, dielectric interface, and PEC conductor rectangles.
// Supports symmetry (half-domain with PMC at x=0).
//
// Usage:
//   const G = await initGmsh();
//   const mesh = generateGmshMesh(G, { domainW, domainH, conductors, interfaceY, hFine, hCoarse });

// ---- Internal helpers ----

function _cstr(G, str) {
    const len = G.lengthBytesUTF8(str) + 1;
    const ptr = G.stackAlloc(len);
    G.stringToUTF8(str, ptr, len);
    return ptr;
}

// ---- Public API ----

// Load and initialize the gmsh WASM module. Call once; pass the returned G
// to generateGmshMesh for every mesh.
export async function initGmsh() {
    const gmshUrl = new URL('../wasm_solver/gmsh.js', import.meta.url).href;
    const G = await (await import(gmshUrl)).default();
    const stack = G.stackSave();
    const ierr  = G.stackAlloc(4);
    G._gmshInitialize(0, 0, 1, 0, ierr);
    G._gmshOptionSetNumber(_cstr(G, 'General.Verbosity'), 0, ierr);
    G.stackRestore(stack);
    return G;
}

