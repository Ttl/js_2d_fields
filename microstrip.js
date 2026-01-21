import { FieldSolver2D, CONSTANTS, diff } from './field_solver.js';
import { Dielectric, Conductor, Mesher } from './mesher.js';

// ============================================================================
// MICROSTRIP SOLVER V2
// ============================================================================

class MicrostripSolver extends FieldSolver2D {
    constructor(options) {
        super();

        // Store parameters
        this.h = options.substrate_height;
        this.w = options.trace_width;
        this.t = options.trace_thickness;
        this.t_gnd = options.gnd_thickness ?? 35e-6;
        this.er = options.epsilon_r;
        this.er_top = options.epsilon_r_top ?? 1;
        this.tan_delta = options.tan_delta ?? 0.02;
        this.sigma_diel = options.sigma_diel ?? 0.0;
        this.sigma_cond = options.sigma_cond ?? 5.8e7;

        // Differential mode parameters
        this.trace_spacing = options.trace_spacing ?? null;
        this.is_differential = (this.trace_spacing !== null && this.trace_spacing > 0);

        this.gnd_cut_width = options.gnd_cut_width ?? 0.0;
        this.gnd_cut_sub_h = options.gnd_cut_sub_h ?? 0.0;
        this.top_diel_h = options.top_diel_h ?? 0.0;
        this.top_diel_er = options.top_diel_er ?? 1.0;
        this.top_diel_tand = options.top_diel_tand ?? 0.0;

        // Solder Mask Parameters
        this.use_sm = options.use_sm ?? false;
        this.sm_t_sub = options.sm_t_sub ?? 20e-6;
        this.sm_t_trace = options.sm_t_trace ?? 20e-6;
        this.sm_t_side = options.sm_t_side ?? 20e-6;
        this.sm_er = options.sm_er ?? 3.5;
        this.sm_tand = options.sm_tand ?? 0.02;

        this.freq = options.freq ?? 1e9;
        this.omega = 2 * Math.PI * this.freq;
        this.nx = options.nx ?? 300;
        this.ny = options.ny ?? 300;

        // Store air parameters
        const air_side = options.air_side ?? null;
        const air_top = options.air_top ?? null;

        // Domain sizing
        if (this.is_differential) {
            // For differential, span includes both traces and spacing
            const trace_span = 2 * this.w + this.trace_spacing;
            if (air_side === null) {
                this.domain_width = 2 * Math.max(trace_span * 4, this.h * 15);
            } else {
                this.domain_width = trace_span + 2 * air_side;
            }
        } else {
            // Single-ended
            if (air_side === null) {
                this.domain_width = 2 * Math.max(this.w * 8, this.h * 15);
            } else {
                this.domain_width = this.w + 2 * air_side;
            }
        }

        this.boundaries = options.boundaries ?? ["open", "open", "open", "gnd"];

        // Calculate physical coordinates
        this._calculate_coordinates(air_top);

        // Skin depth
        this.delta_s = Math.sqrt(2 / (this.omega * CONSTANTS.MU0 * this.sigma_cond));

        // Build geometry lists
        const [dielectrics, conductors] = this._build_geometry_lists();
        this.dielectrics = dielectrics;
        this.conductors = conductors;

        // Create mesher but don't generate mesh yet
        this.mesher = new Mesher(
            this.domain_width, this.domain_height,
            this.nx, this.ny, this.delta_s,
            this.conductors, this.dielectrics,
            true  // symmetric
        );

        // Mesh will be generated when needed
        this.x = null;
        this.y = null;
        this.dx = null;
        this.dy = null;
        this.mesh_generated = false;
    }

    _calculate_coordinates(air_top) {
        // Bottom extension for cut ground
        this.y_ext_start = this.t_gnd;
        this.y_ext_end = this.t_gnd + this.gnd_cut_sub_h;

        // New bottom ground plane location
        this.y_gnd_bot_start = this.y_ext_end;
        this.y_gnd_bot_end = this.y_gnd_bot_start + this.t_gnd;
        if (this.gnd_cut_width === 0) {
            this.y_gnd_bot_end = this.y_gnd_bot_start;
        }

        this.y_sub_start = this.y_gnd_bot_end;
        this.y_sub_end = this.y_sub_start + this.h;

        // Top dielectric
        this.y_top_diel_start = this.y_sub_end;
        this.y_top_diel_end = this.y_top_diel_start + this.top_diel_h;

        // Trace is embedded in top dielectric
        this.y_trace_start = this.y_top_diel_start;
        this.y_trace_end = this.y_trace_start + this.t;

        // Solder mask extents
        this.y_sm_sub_end = this.y_top_diel_end + this.sm_t_sub;
        this.y_sm_trace_end = this.y_trace_end + this.sm_t_trace;

        this.y_top_start = this.y_top_diel_end;
        if (this.use_sm) {
            this.y_top_start = Math.max(this.y_sm_sub_end, this.y_sm_trace_end);
        }

        // Top air/dielectric region
        if (air_top === null) {
            this.top_dielectric_h = this.h * 15;
            this.has_top_gnd = false;
        } else {
            this.top_dielectric_h = air_top + this.t;
            this.has_top_gnd = (this.boundaries[2] === "gnd");
        }

        this.y_top_end = this.y_top_start + this.top_dielectric_h;

        if (this.has_top_gnd) {
            this.y_gnd_top_start = this.y_top_end;
            this.y_gnd_top_end = this.y_gnd_top_start + this.t_gnd;
            this.domain_height = this.y_gnd_top_end;
        } else {
            this.y_gnd_top_start = null;
            this.y_gnd_top_end = null;
            this.domain_height = this.y_top_end;
        }
    }

    _build_geometry_lists() {
        const dielectrics = [];
        const conductors = [];

        const cx = this.domain_width / 2;
        const xl = cx - this.w / 2;
        const xr = cx + this.w / 2;

        // Substrate (covers both cutout extension and main substrate)
        if (this.gnd_cut_sub_h > 0) {
            dielectrics.push(new Dielectric(
                0, this.y_ext_start,
                this.domain_width, this.y_trace_start - this.y_ext_start,
                this.er, this.tan_delta
            ));
        } else {
            dielectrics.push(new Dielectric(
                0, this.y_sub_start,
                this.domain_width, this.h,
                this.er, this.tan_delta
            ));
        }

        // Top dielectric (if present)
        if (this.top_diel_h > 0) {
            dielectrics.push(new Dielectric(
                0, this.y_top_diel_start,
                this.domain_width, this.top_diel_h,
                this.top_diel_er, this.top_diel_tand
            ));
        }

        // Top air/dielectric region
        dielectrics.push(new Dielectric(
            0, this.y_top_start,
            this.domain_width, this.top_dielectric_h,
            this.er_top, 0.0
        ));

        // Solder mask regions (overwrites previous)
        if (this.use_sm) {
            // Solder mask on substrate (full width)
            dielectrics.push(new Dielectric(
                0, this.y_top_diel_end,
                this.domain_width, this.sm_t_sub,
                this.sm_er, this.sm_tand
            ));

            // Solder mask on left side of trace
            const xsl = xl - this.sm_t_side;
            if (xsl >= 0) {
                dielectrics.push(new Dielectric(
                    xsl, this.y_trace_start,
                    this.sm_t_side, this.t + this.sm_t_trace,
                    this.sm_er, this.sm_tand
                ));
            }

            // Solder mask on right side of trace
            const xsr = xr + this.sm_t_side;
            if (xsr <= this.domain_width) {
                dielectrics.push(new Dielectric(
                    xr, this.y_trace_start,
                    this.sm_t_side, this.t + this.sm_t_trace,
                    this.sm_er, this.sm_tand
                ));
            }

            // Solder mask on top of trace
            dielectrics.push(new Dielectric(
                xl, this.y_trace_end,
                this.w, this.sm_t_trace,
                this.sm_er, this.sm_tand
            ));
        }

        // --- CONDUCTORS ---

        // Bottom ground (beneath everything)
        if (this.t_gnd > 0) {
            conductors.push(new Conductor(
                0, 0,
                this.domain_width, this.t_gnd,
                false
            ));
        }

        // Bottom ground plane (above cutout extension)
        if (this.gnd_cut_width === 0) {
            // No cutout - full ground plane
            if (this.y_gnd_bot_end > this.y_gnd_bot_start) {
                conductors.push(new Conductor(
                    0, this.y_gnd_bot_start,
                    this.domain_width, this.t_gnd,
                    false
                ));
            }
        } else {
            // With cutout - ground on sides only
            const cut_l = cx - this.gnd_cut_width / 2;
            const cut_r = cx + this.gnd_cut_width / 2;

            // Left ground
            if (cut_l > 0) {
                conductors.push(new Conductor(
                    0, this.y_gnd_bot_start,
                    cut_l, this.t_gnd,
                    false
                ));
            }

            // Right ground
            if (cut_r < this.domain_width) {
                conductors.push(new Conductor(
                    cut_r, this.y_gnd_bot_start,
                    this.domain_width - cut_r, this.t_gnd,
                    false
                ));
            }
        }

        // Signal trace(s)
        if (this.is_differential) {
            // Left trace (negative in odd mode)
            const xl_left = cx - this.w - this.trace_spacing / 2;
            conductors.push(new Conductor(
                xl_left, this.y_trace_start,
                this.w, this.t,
                true
            ));
            // Right trace (positive in odd mode)
            const xl_right = cx + this.trace_spacing / 2;
            conductors.push(new Conductor(
                xl_right, this.y_trace_start,
                this.w, this.t,
                true
            ));
        } else {
            // Single trace
            conductors.push(new Conductor(
                xl, this.y_trace_start,
                this.w, this.t,
                true
            ));
        }

        // Top ground plane (if present)
        if (this.has_top_gnd) {
            conductors.push(new Conductor(
                0, this.y_gnd_top_start,
                this.domain_width, this.t_gnd,
                false
            ));
        }

        return [dielectrics, conductors];
    }

    ensure_mesh() {
        if (this.mesh_generated) {
            return;
        }

        // Generate mesh
        [this.x, this.y] = this.mesher.generate_mesh();

        // Calculate spacing arrays
        this.dx = new Float64Array(this.x.length - 1);
        for (let i = 0; i < this.x.length - 1; i++) {
            this.dx[i] = this.x[i + 1] - this.x[i];
        }

        this.dy = new Float64Array(this.y.length - 1);
        for (let i = 0; i < this.y.length - 1; i++) {
            this.dy[i] = this.y[i + 1] - this.y[i];
        }

        // Setup geometry
        this._setup_geometry();
        this.mesh_generated = true;
    }

    _setup_geometry() {
        const tol = 1e-11;
        const nx = this.x.length;
        const ny = this.y.length;

        // Initialize field arrays
        this.V = Array(ny).fill().map(() => new Float64Array(nx));
        this.epsilon_r = Array(ny).fill().map(() => new Float64Array(nx).fill(1.0));
        this.signal_mask = Array(ny).fill().map(() => new Uint8Array(nx));
        this.ground_mask = Array(ny).fill().map(() => new Uint8Array(nx));

        // For differential mode, track positive and negative traces separately
        if (this.is_differential) {
            this.signal_p_mask = Array(ny).fill().map(() => new Uint8Array(nx));
            this.signal_n_mask = Array(ny).fill().map(() => new Uint8Array(nx));
        }

        // Apply dielectrics (last overwrites)
        for (const diel of this.dielectrics) {
            for (let i = 0; i < ny; i++) {
                const yc = this.y[i];
                if (yc >= diel.y_min - tol && yc <= diel.y_max + tol) {
                    for (let j = 0; j < nx; j++) {
                        const xc = this.x[j];
                        if (xc >= diel.x_min - tol && xc <= diel.x_max + tol) {
                            this.epsilon_r[i][j] = diel.epsilon_r;
                        }
                    }
                }
            }
        }

        // Apply conductors
        let signal_cond_idx = 0;
        for (const cond of this.conductors) {
            for (let i = 0; i < ny; i++) {
                const yc = this.y[i];
                if (yc >= cond.y_min - tol && yc <= cond.y_max + tol) {
                    for (let j = 0; j < nx; j++) {
                        const xc = this.x[j];
                        if (xc >= cond.x_min - tol && xc <= cond.x_max + tol) {
                            if (cond.is_signal) {
                                this.signal_mask[i][j] = 1;
                                // For differential: first signal is negative, second is positive
                                if (this.is_differential) {
                                    if (signal_cond_idx === 0) {
                                        this.signal_n_mask[i][j] = 1;
                                    } else {
                                        this.signal_p_mask[i][j] = 1;
                                    }
                                }
                            } else {
                                this.ground_mask[i][j] = 1;
                            }
                            this.V[i][j] = cond.voltage;
                        }
                    }
                }
            }
            if (cond.is_signal) {
                signal_cond_idx++;
            }
        }

        // Finalize conductor mask
        this.conductor_mask = Array(ny).fill().map((_, i) => {
            const row = new Uint8Array(nx);
            for (let j = 0; j < nx; j++) {
                row[j] = this.signal_mask[i][j] | this.ground_mask[i][j];
            }
            return row;
        });
    }

    async solve_adaptive(options = {}) {
        /**
         * Override solve_adaptive to handle both single-ended and differential modes.
         *
         * For single-ended: Returns standard results object
         * For differential: Returns differential_results object with odd/even mode parameters
         */
        if (!this.is_differential) {
            // Single-ended mode - use parent implementation
            return super.solve_adaptive(options);
        } else {
            // Differential mode
            return this._solve_adaptive_differential(options);
        }
    }

    async _solve_adaptive_differential(options = {}) {
        /**
         * Adaptive mesh solve for differential mode with odd and even modes.
         *
         * Returns: object with differential parameters
         */
        // Ensure mesh is generated
        if (this.ensure_mesh) {
            this.ensure_mesh();
        }

        const {
            max_iters = 10,
            refine_frac = 0.15,
            energy_tol = 0.05,
            param_tol = 0.1,
            max_nodes = 20000,
            min_converged_passes = 2,
            onProgress = null
        } = options;

        let prev_energy_odd = null;
        let prev_energy_even = null;
        let prev_Z_odd = null;
        let prev_Z_even = null;
        let converged_count = 0;

        let Z_odd, eps_eff_odd, C_odd, C0_odd, Ex_odd, Ey_odd;
        let Z_even, eps_eff_even, C_even, C0_even, Ex_even, Ey_even;

        for (let it = 0; it < max_iters; it++) {
            // Solve odd mode (+1V and -1V)
            const odd_results = await this._solve_single_mode('odd', true);
            Z_odd = odd_results.Z0;
            eps_eff_odd = odd_results.eps_eff;
            C_odd = odd_results.C;
            C0_odd = odd_results.C0;
            Ex_odd = odd_results.Ex;
            Ey_odd = odd_results.Ey;

            // Solve even mode (+1V and +1V)
            const even_results = await this._solve_single_mode('even', true);
            Z_even = even_results.Z0;
            eps_eff_even = even_results.eps_eff;
            C_even = even_results.C;
            C0_even = even_results.C0;
            Ex_even = even_results.Ex;
            Ey_even = even_results.Ey;

            // Energy-based error
            const { energy: energy_odd, rel_error: energy_err_odd } =
                this._compute_energy_error(Ex_odd, Ey_odd, prev_energy_odd);
            const { energy: energy_even, rel_error: energy_err_even } =
                this._compute_energy_error(Ex_even, Ey_even, prev_energy_even);
            const energy_err = Math.max(energy_err_odd, energy_err_even);

            // Parameter-based error
            let param_err;
            if (prev_Z_odd !== null) {
                const z_odd_err = Math.abs(Z_odd - prev_Z_odd) / Math.max(Math.abs(prev_Z_odd), 1e-12);
                const z_even_err = Math.abs(Z_even - prev_Z_even) / Math.max(Math.abs(prev_Z_even), 1e-12);
                param_err = Math.max(z_odd_err, z_even_err);
            } else {
                param_err = 1.0;
            }

            console.log(`Pass ${it + 1}: Energy err=${energy_err.toExponential(3)}, Param err=${param_err.toExponential(3)}, Grid=${this.x.length}x${this.y.length}`);

            // Call progress callback
            if (onProgress) {
                onProgress({
                    iteration: it + 1,
                    max_iterations: max_iters,
                    energy_error: energy_err,
                    param_error: param_err,
                    nodes_x: this.x.length,
                    nodes_y: this.y.length
                });
            }

            // Check convergence
            if (prev_Z_odd !== null) {
                if (energy_err < energy_tol && param_err < param_tol) {
                    converged_count++;
                    if (converged_count >= min_converged_passes) {
                        console.log(`Converged after ${it + 1} passes`);
                        break;
                    }
                } else {
                    converged_count = 0;
                }
            }

            prev_energy_odd = energy_odd;
            prev_energy_even = energy_even;
            prev_Z_odd = Z_odd;
            prev_Z_even = Z_even;

            // Node budget check
            if (this.x.length * this.y.length > max_nodes) {
                console.log("Node budget reached");
                break;
            }

            // Refine mesh (use odd mode fields for refinement)
            if (it !== max_iters - 1) {
                this.refine_mesh(Ex_odd, Ey_odd, refine_frac);
                this._setup_geometry_after_refinement();
            }
        }

        // Calculate losses
        const alpha_c_odd = this._calculate_differential_conductor_loss(Ex_odd, Ey_odd, Z_odd, 'odd');
        const alpha_d_odd = this._calculate_differential_dielectric_loss(Ex_odd, Ey_odd, Z_odd);

        const alpha_c_even = this._calculate_differential_conductor_loss(Ex_even, Ey_even, Z_even, 'even');
        const alpha_d_even = this._calculate_differential_dielectric_loss(Ex_even, Ey_even, Z_even);

        // Store fields as class members (as arrays for differential)
        this.Ex = [Ex_odd, Ex_even];
        this.Ey = [Ey_odd, Ey_even];

        // Differential and common mode impedances
        const Z_diff = 2 * Z_odd;
        const Z_common = Z_even / 2;

        return {
            Z_odd: Z_odd,
            Z_even: Z_even,
            Z_diff: Z_diff,
            Z_common: Z_common,
            eps_eff_odd: eps_eff_odd,
            eps_eff_even: eps_eff_even,
            C_odd: C_odd,
            C_even: C_even,
            alpha_c_odd: alpha_c_odd,
            alpha_c_even: alpha_c_even,
            alpha_d_odd: alpha_d_odd,
            alpha_d_even: alpha_d_even,
            alpha_total_odd: alpha_c_odd + alpha_d_odd,
            alpha_total_even: alpha_c_even + alpha_d_even
        };
    }

    async _solve_single_mode(mode, vacuum_first = true) {
        /**
         * Solve a single mode (odd or even) for differential pair.
         *
         * Parameters:
         * -----------
         * mode : string - 'odd' or 'even'
         * vacuum_first : boolean - Whether to solve vacuum case first for C0 calculation
         *
         * Returns:
         * --------
         * {Z0, eps_eff, C, C0, Ex, Ey}
         */
        const nx = this.x.length;
        const ny = this.y.length;

        // Set voltages based on mode
        this.V = Array(ny).fill().map(() => new Float64Array(nx));

        if (mode === 'odd') {
            // Odd mode: left=-1V, right=+1V
            for (let i = 0; i < ny; i++) {
                for (let j = 0; j < nx; j++) {
                    if (this.signal_n_mask[i][j]) {
                        this.V[i][j] = -1.0;
                    } else if (this.signal_p_mask[i][j]) {
                        this.V[i][j] = 1.0;
                    } else if (this.ground_mask[i][j]) {
                        this.V[i][j] = 0.0;
                    }
                }
            }
        } else {
            // Even mode: both=+1V
            for (let i = 0; i < ny; i++) {
                for (let j = 0; j < nx; j++) {
                    if (this.signal_n_mask[i][j] || this.signal_p_mask[i][j]) {
                        this.V[i][j] = 1.0;
                    } else if (this.ground_mask[i][j]) {
                        this.V[i][j] = 0.0;
                    }
                }
            }
        }

        let C0;
        if (vacuum_first) {
            // Calculate C0 (vacuum capacitance)
            await this.solve_laplace_iterative(true);

            // Use only positive trace for charge calculation
            const orig_signal_mask = this.signal_mask.map(row => row.slice());
            this.signal_mask = this.signal_p_mask.map(row => row.slice());
            C0 = this.calculate_capacitance(true);
            this.signal_mask = orig_signal_mask;
        }

        // Solve with dielectric
        await this.solve_laplace_iterative(false);

        // Use only positive trace for charge calculation
        const orig_signal_mask = this.signal_mask.map(row => row.slice());
        this.signal_mask = this.signal_p_mask.map(row => row.slice());
        const C = this.calculate_capacitance(false);
        this.signal_mask = orig_signal_mask;

        // Calculate fields
        this.compute_fields();
        const Ex = this.Ex.map(row => row.slice());
        const Ey = this.Ey.map(row => row.slice());

        // Calculate impedance
        let eps_eff, Z0;
        if (C0 !== undefined) {
            eps_eff = C / C0;
            Z0 = 1 / (CONSTANTS.C * Math.sqrt(C * C0));
        }

        return { Z0, eps_eff, C, C0, Ex, Ey };
    }

    _calculate_differential_conductor_loss(Ex, Ey, Z0, mode) {
        /**
         * Calculate conductor loss for differential pair.
         */
        const Rs = Math.sqrt(this.omega * CONSTANTS.MU0 / (2 * this.sigma_cond));
        const delta = Math.sqrt(2 / (this.omega * CONSTANTS.MU0 * this.sigma_cond));

        const nx = this.x.length;
        const ny = this.y.length;
        const dx_array = diff(this.x);
        const dy_array = diff(this.y);

        const get_dx = j => (j >= 0 && j < dx_array.length) ? dx_array[j] : dx_array[dx_array.length - 1];
        const get_dy = i => (i >= 0 && i < dy_array.length) ? dy_array[i] : dy_array[dy_array.length - 1];

        let Pc = 0.0;

        for (let i = 1; i < ny - 1; i++) {
            for (let j = 1; j < nx - 1; j++) {
                if (!this.signal_p_mask[i][j] && !this.signal_n_mask[i][j] && !this.ground_mask[i][j]) {
                    continue;
                }

                const neighbors = [
                    [i, j + 1, 'r', i, get_dy],
                    [i, j - 1, 'l', i, get_dy],
                    [i + 1, j, 'u', j, get_dx],
                    [i - 1, j, 'd', j, get_dx]
                ];

                let cell_K_sq = 0.0;
                let cell_dl = 0.0;

                for (const [ni, nj, direction, dl_idx, dl_func] of neighbors) {
                    if (this.signal_p_mask[ni][nj] || this.signal_n_mask[ni][nj] || this.ground_mask[ni][nj]) {
                        continue;
                    }

                    const eps_diel = this.epsilon_r[ni][nj];
                    const Ex_diel = Ex[ni][nj];
                    const Ey_diel = Ey[ni][nj];

                    let E_norm;
                    if (direction === 'r') {
                        E_norm = Ex_diel;
                    } else if (direction === 'l') {
                        E_norm = -Ex_diel;
                    } else if (direction === 'u') {
                        E_norm = Ey_diel;
                    } else if (direction === 'd') {
                        E_norm = -Ey_diel;
                    }

                    const Z0_freespace = Math.sqrt(CONSTANTS.MU0 / CONSTANTS.EPS0);
                    const H_tan = Math.abs(E_norm) * Math.sqrt(eps_diel) / Z0_freespace;
                    const K = H_tan;
                    const dl = dl_func(dl_idx);

                    cell_K_sq += K * K * dl;
                    cell_dl += dl;
                }

                if (cell_dl > 0) {
                    const dP = 0.5 * Rs * cell_K_sq;
                    Pc += dP;
                }
            }
        }

        // Power normalization for 1W
        const Pc_1W = 0.5 * Pc * 2 * Z0;
        const alpha_db_per_m = 8.686 * Pc_1W / 2.0;

        return alpha_db_per_m;
    }

    _calculate_differential_dielectric_loss(Ex, Ey, Z0) {
        /**
         * Calculate dielectric loss for differential pair.
         */
        const nx = this.x.length;
        const ny = this.y.length;
        const dx_array = diff(this.x);
        const dy_array = diff(this.y);

        let Pd = 0.0;

        for (let i = 0; i < ny - 1; i++) {
            for (let j = 0; j < nx - 1; j++) {
                if (this.signal_p_mask[i][j] || this.signal_n_mask[i][j] || this.ground_mask[i][j]) {
                    continue;
                }
                if (this.epsilon_r[i][j] <= 1.01) {
                    continue;
                }

                const E2 = Ex[i][j] * Ex[i][j] + Ey[i][j] * Ey[i][j];
                const dA = dx_array[j] * dy_array[i];
                Pd += 0.5 * this.omega * CONSTANTS.EPS0 * this.epsilon_r[i][j] * this.tan_delta * E2 * dA;
            }
        }

        // Power flow for 1W in the transmission line
        const P_flow = 1.0 / (2 * Z0);

        // Attenuation constant
        return 8.686 * (0.5 * Pd / (2 * P_flow));
    }
}

export { MicrostripSolver };
