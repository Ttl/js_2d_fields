#!/bin/bash
# Link gmsh.{js,wasm} (OCC-enabled, size-trimmed) from the compiled libgmsh.a.
#
# Drops the DataExchange (STEP/IGES) OCC libs + the file-I/O exports
# (gmshOpen/gmshWrite) that are the only paths reaching the STEP/IGES readers, so
# the linker dead-strips that code (gmsh.patch additionally #if-0's the direct
# references in GModelIO_OCC.cpp). The OCC-fragment mesher uses none of it. This
# trims the binary from ~24.5 MB to ~13.7 MB.
#
# Requirements:
#   - emcc on PATH (or set EMSDK to an emsdk dir and its env will be sourced)
#   - OpenCASCADE cross-compiled to WASM (OCCT_WASM, default /tmp/occt-wasm-install)
# Run from this directory (src/wasm_solver/), after `make gmsh-lib`.
set -e
[ -n "$EMSDK" ] && [ -f "$EMSDK/emsdk_env.sh" ] && source "$EMSDK/emsdk_env.sh"
OCCT_WASM="${OCCT_WASM:-/tmp/occt-wasm-install}"
cd "$(dirname "$0")/gmsh/build-wasm-occ"

OCC="$OCCT_WASM/lib"
OCCLIBS="-L$OCC -lTKOffset -lTKFeat -lTKFillet -lTKBool -lTKMesh -lTKHLR -lTKBO -lTKPrim \
 -lTKShHealing -lTKTopAlgo -lTKGeomAlgo -lTKBRep -lTKGeomBase -lTKG3d -lTKG2d -lTKMath -lTKernel"
EXPORTS="['_gmshInitialize','_gmshFinalize','_gmshClear','_gmshModelAdd',\
'_gmshModelGeoAddPoint','_gmshModelGeoAddLine','_gmshModelGeoAddCurveLoop','_gmshModelGeoAddPlaneSurface',\
'_gmshModelGeoSynchronize','_gmshModelGeoMeshSetSize',\
'_gmshModelOccAddPoint','_gmshModelOccAddLine','_gmshModelOccAddCurveLoop','_gmshModelOccAddPlaneSurface',\
'_gmshModelOccAddRectangle','_gmshModelOccFragment','_gmshModelOccSynchronize','_gmshModelOccRemove',\
'_gmshModelGetEntities','_gmshModelAddPhysicalGroup','_gmshModelGetEntitiesForPhysicalGroup',\
'_gmshModelGetBoundary','_gmshModelGetEntitiesInBoundingBox','_gmshModelMeshSetSize','_gmshModelGetBoundingBox',\
'_gmshModelMeshFieldAdd','_gmshModelMeshFieldSetNumber','_gmshModelMeshFieldSetNumbers',\
'_gmshModelMeshFieldSetString','_gmshModelMeshFieldSetAsBackgroundMesh',\
'_gmshModelMeshGenerate','_gmshModelMeshGetNodes','_gmshModelMeshGetElements','_gmshModelMeshGetElementsByType',\
'_gmshOptionSetNumber','_gmshFree','_gmshMalloc']"
emcc libgmsh.a $OCCLIBS -o gmsh.js \
  -Oz -s MODULARIZE=1 -s EXPORT_NAME="GmshModule" -s EXPORT_ES6=1 \
  -s ALLOW_MEMORY_GROWTH=1 -fexceptions -s ASSERTIONS=0 -s INITIAL_MEMORY=67108864 \
  -s EXPORTED_FUNCTIONS="$EXPORTS" \
  -s EXPORTED_RUNTIME_METHODS="['ccall','cwrap','getValue','setValue','UTF8ToString','stringToUTF8','lengthBytesUTF8','stackAlloc','stackSave','stackRestore']"
echo "=== link exit: $? ==="
# Emscripten's getWasmTableEntry calls WebAssembly.Table.get on every indirect
# call routed through an invoke_* trampoline (this build uses JS exceptions,
# -fexceptions, so that is most calls inside OCC). Table.get is slow in V8
# (~1 µs): it was 10 s of a 28 s mesh build. Mirror the table in a JS array, as
# older Emscripten versions did (wasmTableMirror). Static build, no addFunction,
# so the entries never change.
python3 - gmsh.js <<'PY'
import sys
p = sys.argv[1]; s = open(p).read()
old = 'var getWasmTableEntry=funcPtr=>wasmTable.get(funcPtr);'
new = ('var wasmTableMirror=[];var getWasmTableEntry=funcPtr=>{var f=wasmTableMirror[funcPtr];'
       'if(!f){f=wasmTableMirror[funcPtr]=wasmTable.get(funcPtr)}return f};')
assert s.count(old) == 1, 'getWasmTableEntry pattern not found: update the mirror patch'
open(p, 'w').write(s.replace(old, new)); print('patched getWasmTableEntry (table mirror)')
PY
ls -la gmsh.wasm | awk '{printf ">>> gmsh+OCC WASM: %.2f MB\n",$5/1048576}'
