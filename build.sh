#!/usr/bin/env bash
# Optional build script to minify files
set -e

rm -rf dist
mkdir dist

#make -C src/wasm_solver/

# JS
# tri_solver/ is kept external (not bundled): it contains emscripten modules that
# load their .wasm via import.meta.url and a computed dynamic import() for gmsh,
# neither of which survive bundling into a single outfile. It is lazy-imported
# from app_solver.js only when the triangular backend is selected, and shipped as
# a verbatim tree (copied below).
npx esbuild src/app_solver.js \
  --bundle \
  --minify \
  --format=esm \
  --sourcemap \
  --platform=node \
  --external:./tri_solver/* \
  --outfile=dist/app_solver.js

# Solve worker: a second entry point, deliberately not bundled into app_solver.js.
# app_solver constructs it with `new Worker(new URL("solve_worker.js", import.meta.url))`,
# which esbuild leaves alone, the URL is resolved by the browser at runtime, relative to
# the app bundle, so the two files must sit side by side in dist/. It carries its own copy
# of the solver core (both threads need it). The tri_solver tree stays external and shared,
# exactly as for the app bundle.
npx esbuild src/solve_worker.js \
  --bundle \
  --minify \
  --format=esm \
  --sourcemap \
  --platform=node \
  --external:./tri_solver/* \
  --outfile=dist/solve_worker.js

# CSS
npx postcss src/solver-style.css \
  --use cssnano \
  -o dist/solver-style.css

# HTML - swap to CDN first, then minify (order matters!)
# Copy and swap CDN URL first
cp src/field_solver.html dist/field_solver.html
sed -i "s|const PRIMARY_PLOTLY_SRC = 'plotly-3.3.0.min.js';|const PRIMARY_PLOTLY_SRC = 'https://cdn.plot.ly/plotly-3.3.0.min.js';|g" dist/field_solver.html

# Then minify (which will inline the CDN URL)
npx html-minifier-terser \
  --collapse-whitespace \
  --remove-comments \
  --minify-js true \
  --minify-css true \
  dist/field_solver.html \
  -o dist/field_solver.html

# Copy plotly as fallback in case CDN is unreachable
cp src/plotly-3.3.0.min.js dist/

cp src/wasm_solver/solver.wasm dist/solver.wasm
cp src/.htaccess dist/.htaccess

# Full-wave (triangular) backend runtime: shipped as a verbatim ES-module tree —
# NOT bundled, because the emscripten modules locate their .wasm via
# import.meta.url and gmsh is loaded through a computed dynamic import, neither
# of which survives bundling. The WHOLE runtime (tri_solver modules, the shared
# `../` sibling modules, and the eigen/gmsh WASM) lives under ONE content-hashed
# directory, dist/tri-<hash>/, so the immutable cache policy can never serve a
# mixed-version module graph after a deployment: when any file changes, every
# URL in the subtree changes. The only reference from outside is the bundle's
# lazy import('./tri_solver/tri_backend.js'), rewritten below. The layout
# mirrors src/ shifted one level down, so all relative imports resolve inside
# the versioned tree:
#   tri-<hash>/complex.js …                          ← `../foo.js` from tri_solver/*
#   tri-<hash>/tri_solver/*.js
#   tri-<hash>/wasm_solver/{eigen_solver,gmsh}.{js,wasm}
TRI_STAGE=dist/.tri-stage
mkdir -p "$TRI_STAGE/wasm_solver"
cp -r src/tri_solver "$TRI_STAGE/tri_solver"
# dev-only tests (node/playwright) are not shipped
rm -rf "$TRI_STAGE/tri_solver/tests"
rm -f "$TRI_STAGE/tri_solver/_smoke_test.mjs"
# Shared src modules imported as `../*.js` from tri_solver/*.js — keep this list
# in sync with those imports. (They are also inlined into the app bundle; the
# external tree needs standalone copies.)
cp src/complex.js src/surface_roughness.js src/djordjevic_sarkar.js src/matrix.js \
   src/sparameters.js src/geometry_symmetry.js src/shapes.js "$TRI_STAGE/"
# Each emscripten module locates its .wasm next to its own .js.
cp src/wasm_solver/eigen_solver.js src/wasm_solver/eigen_solver.wasm \
   src/wasm_solver/gmsh.js src/wasm_solver/gmsh.wasm "$TRI_STAGE/wasm_solver/"
TRI_HASH=$(cd "$TRI_STAGE" && find . -type f | LC_ALL=C sort | xargs cat | sha256sum | cut -c1-8)
mv "$TRI_STAGE" "dist/tri-${TRI_HASH}"
# Point both bundles' lazy import at the versioned tree. The worker sits next to the app
# bundle, so the same relative path resolves for both.
sed -i "s|\./tri_solver/tri_backend\.js|./tri-${TRI_HASH}/tri_solver/tri_backend.js|g" \
  dist/app_solver.js dist/solve_worker.js

# Cache-busting: append content hashes as query strings so browsers fetch
# new versions when files change. The HTML must be served with no-cache headers
# (e.g. Cache-Control: no-cache) so users always get the latest asset URLs.
BUILD_DATE=$(date -u '+%Y-%m-%d')
GIT_HASH=$(git rev-parse --short HEAD)
sed -i "s|BUILD_VERSION|${BUILD_DATE} (${GIT_HASH})|g" dist/field_solver.html

CSS_HASH=$(sha256sum dist/solver-style.css | cut -c1-8)
WASM_HASH=$(sha256sum dist/solver.wasm | cut -c1-8)

# Inject WASM hash into bundled JS (solver.wasm is fetched at runtime).
# This must happen BEFORE JS_HASH is computed: a solver.wasm-only change must
# alter the bundle (and therefore its URL), or returning users would keep the
# old immutable-cached bundle pointing at the old ?v= and never fetch the new
# WASM.
sed -i "s|\"solver\.wasm\"|\"solver.wasm?v=${WASM_HASH}\"|g" dist/app_solver.js dist/solve_worker.js

# The worker is fetched by URL from inside the app bundle, so its hash must be injected
# BEFORE the app bundle is hashed — same ordering argument as solver.wasm above.
WORKER_HASH=$(sha256sum dist/solve_worker.js | cut -c1-8)
sed -i "s|\"solve_worker\.js\"|\"solve_worker.js?v=${WORKER_HASH}\"|g" dist/app_solver.js

JS_HASH=$(sha256sum dist/app_solver.js | cut -c1-8)

# Inject JS and CSS hashes into HTML
sed -i \
  -e "s|\"app_solver\.js\"|\"app_solver.js?v=${JS_HASH}\"|g" \
  -e "s|\"solver-style\.css\"|\"solver-style.css?v=${CSS_HASH}\"|g" \
  dist/field_solver.html
