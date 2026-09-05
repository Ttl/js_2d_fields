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

1. Cross-compile OCCT to WASM and install it (e.g. to `/tmp/occt-wasm-install`).
   Use OCCT 7.5.1 with `emcmake`, static libraries, and the modeling toolkits
   (FoundationClasses + ModelingData + ModelingAlgorithms + DataExchange);
   `-DBUILD_LIBRARY_TYPE=Static -DBUILD_MODULE_Visualization=OFF
   -DBUILD_MODULE_ApplicationFramework=OFF`. Build with `-Os` for a smaller
   binary.

2. Build gmsh against it:

   ```sh
   OCCT_WASM=/tmp/occt-wasm-install make gmsh
   ```

   This applies `gmsh.patch` (disables the unused STEP/IGES readers to trim the
   binary), CMake-configures gmsh with OCC, compiles `libgmsh.a`, then links and
   dead-strips via `link-occ-wasm.sh`, producing `gmsh.{js,wasm}` (~13.7 MB).

`gmsh.patch` is a small source patch against the pinned gmsh commit; it only
removes unused code (size optimization) and is not required for correctness.

`link-occ-wasm.sh` also post-processes the generated `gmsh.js`: it replaces
Emscripten's `getWasmTableEntry` with a JS-array mirror of the function table.
The build uses JS exception handling (`-fexceptions`), so most calls inside OCC
go through `invoke_*` trampolines that look their target up with
`WebAssembly.Table.get`, which was ~35% of a mesh build's wall time. The patched
`gmsh.js` is what is checked in; rebuilding without the script loses it (a
rebuild with `-fwasm-exceptions` and `-O2` instead of `-Oz` would remove the
trampolines altogether, at the cost of a larger binary, and is untested).
