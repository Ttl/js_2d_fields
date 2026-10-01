// The eigen_solver WASM module and its helpers, one instance shared by the full-wave
// backend and the quasi-static complex solve (conductive dielectrics).
import createModule from '../wasm_solver/eigen_solver.js';
import { createWasmHelpers } from './fem_core.js';

let _promise = null;
export function eigenWasm() {
    if (!_promise) _promise = createModule().then(M => ({ M, helpers: createWasmHelpers(M) }));
    return _promise;
}
