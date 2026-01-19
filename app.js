const App = {
    solver: null,
    canvas: document.getElementById('sim_canvas'),
    ctx: document.getElementById('sim_canvas').getContext('2d'),

    init() {
        this.bindEvents();
        this.updateGeometry();
        this.resizeCanvas();
        window.addEventListener('resize', () => this.resizeCanvas());
    },

    bindEvents() {
        document.getElementById('btn_update_geo').onclick = () => {
            this.updateGeometry();
            this.draw();
        };
        document.getElementById('btn_solve').onclick = () => this.runSimulation();

        ['chk_mesh', 'sel_view_mode'].forEach(id => {
            document.getElementById(id).onchange = () => this.draw();
        });
    },

    getParams() {
        return {
            w: parseFloat(document.getElementById('inp_w').value) * 1e-3,
            h: parseFloat(document.getElementById('inp_h').value) * 1e-3,
            t: parseFloat(document.getElementById('inp_t').value) * 1e-6,
            er: parseFloat(document.getElementById('inp_er').value),
            tand: parseFloat(document.getElementById('inp_tand').value),
            sigma: parseFloat(document.getElementById('inp_sigma').value),
            freq: parseFloat(document.getElementById('inp_freq').value) * 1e9,
            nx: parseInt(document.getElementById('inp_nx').value),
            ny: parseInt(document.getElementById('inp_ny').value),
        };
    },

    updateGeometry() {
        const p = this.getParams();
        this.solver = new MicrostripSolver(p.w, p.h, p.t, p.er, p.tand, p.sigma, p.freq, p.nx, p.ny);
        this.log("Geometry updated. Grid: " + this.solver.x.length + "x" + this.solver.y.length);
    },

    async runSimulation() {
        const btn = document.getElementById('btn_solve');
        const pbar = document.getElementById('progress_bar');
        const pcont = document.getElementById('progress_container');

        btn.disabled = true;
        pcont.style.display = 'block';
        this.log("Starting simulation...");

        try {
            // 1. Solve with Dielectric
            this.log("Step 1/2: Solving field with dielectric...");
            await this.solver.solve_laplace_iterative(false, (i, max, diff) => {
                const pct = (i / max) * 50;
                pbar.style.width = pct + "%";
            });
            const C = this.solver.calculate_capacitance(false);

            // Compute fields for loss calculation and visualization
            this.solver.compute_fields();
            this.draw(); // Update view with fields

            // 2. Solve with Vacuum
            this.log("Step 2/2: Solving field with vacuum...");

            // Backup Voltage from Step 1 to restore for visualization later
            const V_diel = this.solver.V.map(row => new Float64Array(row));
            const Ex_diel = this.solver.Ex;
            const Ey_diel = this.solver.Ey;

            // Reset V for vacuum solve
            this.solver.V = this.solver.V.map((row, i) => 
                row.map((val, j) => this.solver.conductor_mask[i][j] ? val : 0)
            );

            await this.solver.solve_laplace_iterative(true, (i, max, diff) => {
                const pct = 50 + (i / max) * 50;
                pbar.style.width = pct + "%";
            });
            const C0 = this.solver.calculate_capacitance(true);

            // 3. Post Process
            const eps_eff = C / C0;
            const Z0 = 1 / (CONSTANTS.C * Math.sqrt(C * C0));

            // Restore dielectric fields for loss calc
            this.solver.V = V_diel;
            this.solver.Ex = Ex_diel;
            this.solver.Ey = Ey_diel;

            const losses = this.solver.calculate_losses(Z0);

            this.log(`\nRESULTS:\n` +
                     `----------------------\n` +
                     `Capacitance:   ${(C*1e12).toFixed(2)} pF/m\n` +
                     `Z0:            ${Z0.toFixed(2)} Ω\n` +
                     `Eps_eff:       ${eps_eff.toFixed(3)}\n` +
                     `Dielectric Loss: ${losses.alpha_diel.toFixed(4)} dB/m\n` +
                     `Conductor Loss:  ${losses.alpha_cond.toFixed(4)} dB/m\n` +
                     `Total Loss:      ${(losses.alpha_diel + losses.alpha_cond).toFixed(4)} dB/m`);

        } catch (e) {
            console.error(e);
            this.log("Error: " + e.message);
        } finally {
            btn.disabled = false;
            pcont.style.display = 'none';
        }
    },

    resizeCanvas() {
        const container = this.canvas.parentElement;
        this.canvas.width = container.clientWidth;
        this.canvas.height = container.clientHeight;
        if(this.solver) this.draw();
    },

    draw() {
        if (!this.solver) return;
        const w = this.canvas.width;
        const h = this.canvas.height;
        this.ctx.clearRect(0, 0, w, h);

        const viewMode = document.getElementById('sel_view_mode').value;
        const showMesh = document.getElementById('chk_mesh').checked;

        // Scaling
        // Find bounding box to fit in canvas with margin
        const domainW = this.solver.domain_width;
        // Limit display Y to relevant area (up to 5x trace height)
        const maxY = Math.min(this.solver.y.slice(-1)[0], this.solver.h + 5 * this.solver.w);
        const domainH = maxY;

        const scaleX = (w - 40) / domainW;
        const scaleY = (h - 40) / domainH;
        const scale = Math.min(scaleX, scaleY);

        const offsetX = (w - domainW * scale) / 2;
        const offsetY = h - (h - domainH * scale) / 2; // Bottom origin

        const toScreenX = (x) => offsetX + x * scale;
        const toScreenY = (y) => offsetY - y * scale;

        // Helper to draw Rects
        const drawCell = (i, j, color) => {
            const x0 = toScreenX(this.solver.x[j]);
            const y0 = toScreenY(this.solver.y[i]);
            // dx/dy are widths
            const x1 = toScreenX(this.solver.x[j+1]);
            const y1 = toScreenY(this.solver.y[i+1]);

            this.ctx.fillStyle = color;
            this.ctx.fillRect(x0, y1, x1-x0, y0-y1); // Canvas Y inverted
        };

        // Draw Field / Geometry
        const ny = this.solver.y.length;
        const nx = this.solver.x.length;

        // Determine min/max for colormapping
        let maxVal = 0;
        if (viewMode === 'efield' && this.solver.solution_valid) {
             for(let i=0; i<ny; i++) 
                for(let j=0; j<nx; j++) 
                     maxVal = Math.max(maxVal, Math.sqrt(this.solver.Ex[i][j]**2 + this.solver.Ey[i][j]**2));
        }

        for (let i = 0; i < ny - 1; i++) {
            if (this.solver.y[i] > maxY) break;

            for (let j = 0; j < nx - 1; j++) {
                let color;

                if (viewMode === 'geometry') {
                    if (this.solver.conductor_mask[i][j]) {
                        color = this.solver.ground_mask[i][j] ? '#333' : '#d97706'; // Dark gnd, Orange trace
                    } else {
                        // Substrate vs Air
                        const er = this.solver.epsilon_r[i][j];
                        color = er > 1.1 ? `rgba(100, 200, 100, ${0.2 + er/10})` : '#fff';
                    }
                } else if (viewMode === 'potential' && this.solver.solution_valid) {
                    if (this.solver.conductor_mask[i][j]) {
                         color = this.solver.ground_mask[i][j] ? '#000' : '#d97706';
                    } else {
                        const v = this.solver.V[i][j];
                        color = this.viridis(v);
                    }
                } else if (viewMode === 'efield' && this.solver.solution_valid) {
                    const emag = Math.sqrt(this.solver.Ex[i][j]**2 + this.solver.Ey[i][j]**2);
                    color = this.viridis(Math.pow(emag / (maxVal || 1), 0.5)); // Gamma correct
                } else {
                    color = '#eee';
                }

                drawCell(i, j, color);
            }
        }

        // Draw Mesh Lines
        if (showMesh) {
            this.ctx.strokeStyle = 'rgba(0,0,0,0.1)';
            this.ctx.lineWidth = 1;
            this.ctx.beginPath();

            // V lines
            for(let j=0; j<nx; j++) {
                const x = toScreenX(this.solver.x[j]);
                this.ctx.moveTo(x, toScreenY(0));
                this.ctx.lineTo(x, toScreenY(maxY));
            }
            // H lines
            for(let i=0; i<ny; i++) {
                if (this.solver.y[i] > maxY) break;
                const y = toScreenY(this.solver.y[i]);
                this.ctx.moveTo(toScreenX(0), y);
                this.ctx.lineTo(toScreenX(domainW), y);
            }
            this.ctx.stroke();
        }
    },

    viridis(t) {
        // Simple heatmap approximation
        t = Math.max(0, Math.min(1, t));
        // R, G, B interpolation
        const r = Math.floor(255 * Math.sin(t * 2));
        const g = Math.floor(255 * Math.sin(t * 3));
        const b = Math.floor(255 * Math.cos(t * 1.5));
        // Better pseudocolor: (Blue -> Cyan -> Green -> Yellow -> Red)
        // Manual standard mapping for clarity:
        if(t < 0.25) return `rgb(0, ${Math.floor(t*4*255)}, 255)`;
        if(t < 0.5) return `rgb(0, 255, ${Math.floor((0.5-t)*4*255)})`;
        if(t < 0.75) return `rgb(${Math.floor((t-0.5)*4*255)}, 255, 0)`;
        return `rgb(255, ${Math.floor((1-t)*4*255)}, 0)`;
    },

    log(msg) {
        const c = document.getElementById('console_out');
        c.textContent += msg + "\n";
        c.scrollTop = c.scrollHeight;
    }
};


// Start
window.onload = () => App.init();

