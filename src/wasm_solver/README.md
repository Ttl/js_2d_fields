# WASM solver libraries

C++/WebAssembly numerics used by the field solver. All three runtime binaries
are checked in; rebuilding is only needed when the sources or pinned dependency
commits change.

| Output                  | Source            | Used by                  |
| ----------------------- | ----------------- | ------------------------ |
| `solver.{js,wasm}`      | `solver.cpp`      | quasi-static FDM backend |
| `eigen_solver.{js,wasm}`| `eigen_solver.cpp`| full-wave (triangular) backend |
| `gmsh.{js,wasm}`        | `gmsh/` submodule + `gmsh.patch` | triangular-mesh generation |

The triangular backend's JS modules live in `../tri_solver/` and load
`eigen_solver.js` and `gmsh.js` from this directory.

## Dependencies (git submodules)

- `eigen/`   — Eigen (linear algebra), used by both solvers
- `spectra/` — Spectra (eigensolvers), used by `eigen_solver.cpp`
- `gmsh/`    — Gmsh meshing kernel (pinned to `gmsh_4_15_0-179-g77cc784cc`)

```sh
git submodule update --init src/wasm_solver/eigen src/wasm_solver/spectra src/wasm_solver/gmsh
```

## Building

Requires [emscripten](https://emscripten.org/) (`emcc` / `emcmake` on `PATH`).

```sh
make solver        # solver.{js,wasm}
make eigen_solver  # eigen_solver.{js,wasm}
make all           # both of the above
```

### gmsh (OCC-enabled)

Rebuilding gmsh is heavier: it needs **OpenCASCADE (OCCT) cross-compiled to
WASM**, which gmsh links statically for the OCC-fragment mesher.

1. Cross-compile OCCT 7.5.1 to WASM and install it (default `/tmp/occt-wasm-install`).
   gmsh's configure only looks for the STEP/IGES base toolkits, and the rest of
   DataExchange (XDE, VRML, STL) would pull in Visualization and FreeType, so trim
   the module list first, then configure with the same exception/optimization
   flags as gmsh:

   ```sh
   git clone --depth 1 --branch V7_5_1 https://github.com/Open-Cascade-SAS/OCCT.git occt
   sed -i 's/^DataExchange .*/DataExchange TKXSBase TKSTEPBase TKSTEPAttr TKSTEP209 TKSTEP TKIGES/' occt/adm/MODULES
   emcmake cmake -S occt -B occt/build -G Ninja -DCMAKE_BUILD_TYPE=Release -DBUILD_LIBRARY_TYPE=Static \
     -DBUILD_MODULE_FoundationClasses=ON -DBUILD_MODULE_ModelingData=ON \
     -DBUILD_MODULE_ModelingAlgorithms=ON -DBUILD_MODULE_DataExchange=ON \
     -DBUILD_MODULE_Visualization=OFF -DBUILD_MODULE_ApplicationFramework=OFF -DBUILD_MODULE_Draw=OFF \
     -DUSE_FREETYPE=OFF -DUSE_TK=OFF -DUSE_TBB=OFF -DUSE_FREEIMAGE=OFF -DUSE_RAPIDJSON=OFF \
     -DUSE_OPENGL=OFF -DUSE_GLES2=OFF -DBUILD_DOC_Overview=OFF \
     -DINSTALL_DIR=/tmp/occt-wasm-install -DINSTALL_DIR_LAYOUT=Unix \
     -DCMAKE_CXX_FLAGS=-fwasm-exceptions -DCMAKE_C_FLAGS=-fwasm-exceptions \
     -DCMAKE_CXX_FLAGS_RELEASE="-O2 -DNDEBUG" -DCMAKE_C_FLAGS_RELEASE="-O2 -DNDEBUG"
   ninja -C occt/build -j16 install
   ```

   About 20 minutes on 16 cores.

2. Build gmsh against it:

   ```sh
   OCCT_WASM=/tmp/occt-wasm-install make gmsh
   ```

   This applies `gmsh.patch` (disables the unused STEP/IGES readers to trim the
   binary), CMake-configures gmsh with OCC, compiles `libgmsh.a`, then links and
   dead-strips via `link-occ-wasm.sh`, producing `gmsh.{js,wasm}` (~11.1 MB).

   Everything is compiled and linked with `-fwasm-exceptions` at `-O2`. The
   default JS exception handling routes calls through JS trampolines, which was
   about a third of a mesh build's time; the wasm-exception build meshes ~1.6x
   faster and produces bit-identical meshes.

`gmsh.patch` is a small source patch against the pinned gmsh commit; it only
removes unused code (size optimization) and is not required for correctness.

`link-occ-wasm.sh` also post-processes `gmsh.js` when it detects a JS-exception
build (`invoke_*` trampolines present): it replaces Emscripten's
`getWasmTableEntry` with a JS-array mirror of the function table, since the
`WebAssembly.Table.get` lookups in the trampolines dominate otherwise. A
wasm-exception build has no trampolines and is left untouched.
