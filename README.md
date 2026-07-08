# 2D Transmission Line Field Solver

![Header](https://github.com/Ttl/js_2d_fields/blob/master/docs/header.png?raw=true)

A browser-based quasi-static and full-wave 2D field solver for transmission line analysis. Computes characteristic impedance, effective permittivity, RLGC parameters, losses, and S-parameters.

Try it online: https://hforsten.com/field_solver.html

## Features

- Transmission Line Types: Microstrip, stripline, GCPW (single-ended and differential)
- Electric Field Solving: Calculates characteristic impedance
- Full RLGC Extraction: Resistance, inductance, capacitance, conductance per unit length
- Loss Modeling: Conductor losses (skin effect, surface roughness) and dielectric losses
- S-Parameter Export: Touchstone .s2p and .s4p file generation
- Visualization: 2D potential plots, E-field streamlines, frequency-dependent plots
- Adaptive Meshing: Automatic mesh refinement for accurate field solutions

## Quick Start

1. Build (or obtain) the WASM solver binaries (see below). They are not checked into git.
2. Host src folder on a web server (for example `python -m http.server 8000`). Open `src/field_solver.html` in a browser.

### Validity

- Designed for microstrip, stripline, and coplanar waveguide structures commonly used in PCB RF and high-speed digital designs.
- Provides good results when the transmission line supports TEM or quasi-TEM propagation, which covers most practical PCB geometries below the onset of higher-order modes.
- Results have been checked against EM solver and actual measurement data with different geometries and transmission line types, showing close agreement with small error in typical use cases (`tests` folder).
- Accurate from RF through microwave and high-speed digital frequencies where return currents are confined and skin effect is significant.
- Sufficient for Most Practical Designs. Suitable for impedance control, loss estimation, and S-parameter generation in the vast majority of PCB transmission line applications.
- Dispersion and higher-order modes are support with full-wave solver.

## Common Tasks

### Run Tests

Tests are organized into tiers and run through `tests/run.mjs`:

```bash
npm run test:fast   # ~1 min, run before every commit
npm run test:slow   # full solver validation, ~5 min
npm test            # fast + slow
npm run test:e2e    # browser tests (see below)
npm run test:fuzz   # QS-vs-fullwave fuzzer (~10 min, deterministic seed)
```

Fast tier: Tests small parts of the solver.
Slow tier: Test against reference structures.
Browser (e2e) tier, needs `playwright-core` (`npm install`).

### Build WASM Solvers

The solver numerics run as WebAssembly (WASM) modules built with Emscripten. All
of them live in `src/wasm_solver/`:

| Output                   | Source             | Used by                  |
| ------------------------ | ------------------ | ------------------------ |
| `solver.{js,wasm}`       | `solver.cpp`       | quasi-static FDM backend |
| `eigen_solver.{js,wasm}` | `eigen_solver.cpp` | full-wave backend        |
| `gmsh.{js,wasm}`         | `gmsh/` submodule  | full-wave mesh generation |

#### Prerequisites

**Emscripten Compiler (emcc)**

Install the Emscripten SDK:

```bash
# Clone the emsdk repository
git clone https://github.com/emscripten-core/emsdk.git
cd emsdk

# Install and activate the latest SDK
./emsdk install latest
./emsdk activate latest

# Add to PATH (add this to your .bashrc or .zshrc for permanent use)
source ./emsdk_env.sh
```

Verify installation:
```bash
emcc --version
```

**Submodules (Eigen, Spectra, gmsh)**

Third-party sources are git submodules under `src/wasm_solver/`. Initialize them:

```bash
git submodule update --init --recursive
```

#### Build Steps

```bash
cd src/wasm_solver
make all          # builds solver.{js,wasm} and eigen_solver.{js,wasm}
```

Building gmsh is heavier: it additionally needs OpenCASCADE cross-compiled to
WASM. See `src/wasm_solver/README.md` for the OCCT prerequisite, then:

```bash
OCCT_WASM=/path/to/occt-wasm-install make gmsh   # builds gmsh.{js,wasm}
```
