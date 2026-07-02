# Build the WebAssembly solver libraries for js_2d_fields.
#
#   make solver        -> solver.{js,wasm}        quasi-static FDM sparse solver
#   make eigen_solver  -> eigen_solver.{js,wasm}  full-wave FEM eigen / QEP solver
#   make all           -> solver + eigen_solver
#   make gmsh          -> gmsh.{js,wasm}          OCC-enabled 2D mesher (see README)
#
# Requirements: emscripten (emcc / emcmake on PATH). `make gmsh` additionally
# needs OpenCASCADE cross-compiled to WASM -- see README.md (OCCT_WASM below).
#
# Third-party sources are git submodules: eigen/, spectra/, gmsh/. Run
# `git submodule update --init` first.

EMCC ?= emcc

all: solver eigen_solver

# --- quasi-static solver -----------------------------------------------------
solver: solver.js
solver.js: solver.cpp | eigen
	$(EMCC) solver.cpp -I eigen -O3 -fwasm-exceptions \
	  -s WASM=1 -s EXPORT_ES6=1 -s ENVIRONMENT=web,node \
	  -s EXPORTED_FUNCTIONS='["_solve_sparse","_malloc","_free"]' \
	  -s EXPORTED_RUNTIME_METHODS='["HEAP32","HEAPF64"]' \
	  -s ALLOW_MEMORY_GROWTH=1 -msimd128 -s STACK_SIZE=33554432 \
	  -o solver.js

# --- full-wave eigen solver --------------------------------------------------
# -fwasm-exceptions: C++ throws (Spectra invalid_argument, bad_alloc) must be
# catchable inside the entry points' try/catch — without it any throw abort()s
# and permanently poisons the module for the rest of the session.
eigen_solver: eigen_solver.js
eigen_solver.js: eigen_solver.cpp | eigen spectra
	$(EMCC) eigen_solver.cpp -I eigen -I spectra/include -O3 -fwasm-exceptions \
	  -s WASM=1 -s EXPORT_ES6=1 -s ENVIRONMENT=web,node \
	  -s EXPORTED_FUNCTIONS='["_solve_generalized_eigen","_solve_generalized_eigen_with_init","_solve_sparse_multi","_solve_rqi","_solve_qep_rqi","_malloc","_free"]' \
	  -s EXPORTED_RUNTIME_METHODS='["HEAP32","HEAPF64"]' \
	  -s ALLOW_MEMORY_GROWTH=1 -msimd128 -s STACK_SIZE=33554432 \
	  -o eigen_solver.js

# --- gmsh mesher (OCC-enabled) ----------------------------------------------
# gmsh.{js,wasm} is vendored (checked in) because rebuilding it needs OpenCASCADE
# cross-compiled to WASM. To rebuild from the gmsh submodule:
#
#   OCCT_WASM=/path/to/occt-wasm-install make gmsh
#
# This applies the size-trim patch, CMake-configures gmsh with OCC, compiles
# libgmsh.a, then links + dead-strips via link-occ-wasm.sh. See README.md for
# building the OCCT WASM prerequisite.

OCCT_WASM  ?= /tmp/occt-wasm-install
GMSH_BUILD := gmsh/build-wasm-occ

gmsh-patch:
	cd gmsh && (git apply --check ../gmsh.patch 2>/dev/null && git apply ../gmsh.patch \
	  && echo "applied gmsh.patch" || echo "gmsh.patch already applied / not applicable; skipping")

gmsh-configure: gmsh-patch
	mkdir -p $(GMSH_BUILD)
	cd $(GMSH_BUILD) && emcmake cmake .. \
	  -DCMAKE_BUILD_TYPE=Release \
	  -DENABLE_BUILD_LIB=ON -DENABLE_BUILD_SHARED=OFF -DENABLE_BUILD_DYNAMIC=OFF \
	  -DENABLE_OCC=ON -DOCC_STATIC=ON \
	  -DCMAKE_FIND_ROOT_PATH=$(OCCT_WASM) \
	  -DCMAKE_FIND_ROOT_PATH_MODE_INCLUDE=BOTH \
	  -DCMAKE_FIND_ROOT_PATH_MODE_LIBRARY=BOTH \
	  -DENABLE_FLTK=OFF -DENABLE_OPENGL=OFF -DENABLE_OCC_CAF=OFF -DENABLE_CAIRO=OFF \
	  -DENABLE_MPEG_ENCODE=OFF -DENABLE_OPENMP=OFF -DENABLE_NETGEN=OFF -DENABLE_HXT=OFF \
	  -DENABLE_DOMHEX=OFF -DENABLE_NII2MESH=OFF -DENABLE_BAMG=OFF -DENABLE_BLOSSOM=OFF \
	  -DENABLE_QUADTRI=OFF -DENABLE_METIS=OFF -DENABLE_PARSER=OFF -DENABLE_PLUGINS=OFF \
	  -DENABLE_ONELAB=OFF -DENABLE_ONELAB_METAMODEL=OFF -DENABLE_GETDP=OFF -DENABLE_KBIPACK=OFF \
	  -DENABLE_GMP=OFF -DENABLE_DINTEGRATION=OFF -DENABLE_TINYXML2=OFF -DENABLE_PRO=OFF \
	  -DENABLE_TESTS=OFF -DENABLE_GMM=OFF -DENABLE_OPTHOM=OFF -DENABLE_UNTANGLE=OFF \
	  -DENABLE_ALGLIB=OFF -DENABLE_WINSLOWUNTANGLER=OFF -DENABLE_VOROPP=OFF

gmsh-lib: gmsh-configure
	cd $(GMSH_BUILD) && emmake make -j$$(nproc) lib

gmsh: gmsh-lib
	OCCT_WASM=$(OCCT_WASM) ./link-occ-wasm.sh
	cp $(GMSH_BUILD)/gmsh.js $(GMSH_BUILD)/gmsh.wasm .

clean:
	rm -f solver.js solver.wasm eigen_solver.js eigen_solver.wasm

.PHONY: all solver eigen_solver gmsh gmsh-patch gmsh-configure gmsh-lib clean
