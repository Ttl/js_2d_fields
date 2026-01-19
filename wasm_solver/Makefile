all:
	emcc solver.cpp -I eigen -O3   -s WASM=1   -s EXPORT_ES6=1   -s ENVIRONMENT=web,node   -s EXPORTED_FUNCTIONS='["_solve_sparse", "_malloc", "_free"]'   -s EXPORTED_RUNTIME_METHODS='["HEAP32","HEAPF64"]'   -s ALLOW_MEMORY_GROWTH=1   -msimd128   -s STACK_SIZE=33554432   -o solver.js
