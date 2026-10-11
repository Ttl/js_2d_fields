// Development server: serves a directory with caching disabled, so every reload gets
// the current modules, workers and WASM (python -m http.server lets the browser cache
// them heuristically).
//
//   node serve.mjs [port] [dir]      default 8000, the current directory
import http from 'http';
import fs from 'fs';
import path from 'path';

const port = parseInt(process.argv[2] || '8000', 10);
const root = path.resolve(process.argv[3] || '.');

const TYPES = {
    '.html': 'text/html; charset=utf-8', '.js': 'text/javascript; charset=utf-8',
    '.mjs': 'text/javascript; charset=utf-8', '.css': 'text/css; charset=utf-8',
    '.json': 'application/json', '.map': 'application/json', '.wasm': 'application/wasm',
    '.svg': 'image/svg+xml', '.png': 'image/png', '.ico': 'image/x-icon', '.txt': 'text/plain; charset=utf-8',
};

http.createServer((req, res) => {
    const url = new URL(req.url, 'http://localhost');
    let file = path.join(root, decodeURIComponent(url.pathname));
    if (file !== root && !file.startsWith(root + path.sep)) { res.writeHead(403).end(); return; }
    fs.stat(file, (err, st) => {
        if (!err && st.isDirectory()) file = path.join(file, 'index.html');
        fs.readFile(file, (err2, data) => {
            if (err2) { res.writeHead(404, { 'Content-Type': 'text/plain' }).end('Not found'); return; }
            res.writeHead(200, {
                'Content-Type': TYPES[path.extname(file).toLowerCase()] || 'application/octet-stream',
                'Cache-Control': 'no-store',
            });
            res.end(data);
        });
    });
}).listen(port, () => {
    // Point at the app wherever it sits under the served directory.
    const page = ['field_solver.html', 'src/field_solver.html'].find(p => fs.existsSync(path.join(root, p))) || '';
    console.log(`Serving ${root} without caching: http://localhost:${port}/${page}`);
});
