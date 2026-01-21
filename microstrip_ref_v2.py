#!/usr/bin/env python
import sys
import numpy as np
from field_solver_ref import FieldSolver2D, epsilon_0, mu_0, pi, c
from mesher import Dielectric, Conductor, Mesher

# ============================================================================
# MICROSTRIP SOLVER V2
# ============================================================================

class MicrostripSolver2D(FieldSolver2D):
    """
    Microstrip transmission line solver using the new Mesher class.

    This version uses lists of dielectrics and conductors for more flexible
    geometry specification.
    """

    def __init__(self, substrate_height, trace_width, trace_thickness,
                 gnd_thickness=35e-6, epsilon_r=4.5, tan_delta=0.02,
                 sigma_diel=0.0, sigma_cond=5.8e7, epsilon_r_top=1,
                 air_top=None, air_side=None, freq=1e9,
                 nx=300, ny=300, boundaries=None,
                 use_sm=False, sm_t_sub=20e-6, sm_t_trace=20e-6,
                 sm_t_side=20e-6, sm_er=3.5, sm_tand=0.02,
                 gnd_cut_width=0.0,
                 gnd_cut_sub_h=0.0,
                 top_diel_h=0.0,
                 top_diel_er=1.0,
                 top_diel_tand=0.0,
                 trace_spacing=None):

        # Store parameters
        self.h = substrate_height
        self.w = trace_width
        self.t = trace_thickness
        self.t_gnd = gnd_thickness
        self.er = epsilon_r
        self.er_top = epsilon_r_top
        self.tan_delta = tan_delta
        self.sigma_diel = sigma_diel
        self.sigma_cond = sigma_cond

        # Differential mode parameters
        self.trace_spacing = trace_spacing
        self.is_differential = (trace_spacing is not None and trace_spacing > 0)

        self.gnd_cut_width = gnd_cut_width
        self.gnd_cut_sub_h = gnd_cut_sub_h
        self.top_diel_h = top_diel_h
        self.top_diel_er = top_diel_er
        self.top_diel_tand = top_diel_tand

        # Solder Mask Parameters
        self.use_sm = use_sm
        self.sm_t_sub = sm_t_sub
        self.sm_t_trace = sm_t_trace
        self.sm_t_side = sm_t_side
        self.sm_er = sm_er
        self.sm_tand = sm_tand

        self.f = freq
        self.omega = 2 * pi * freq
        self.nx, self.ny = nx, ny

        # Store air parameters
        self.air_top = air_top
        self.air_side = air_side

        # Domain sizing
        if self.is_differential:
            # For differential, span includes both traces and spacing
            trace_span = 2 * self.w + self.trace_spacing
            if air_side is None:
                self.domain_width = 2 * max(trace_span * 4, self.h * 15)
            else:
                self.domain_width = trace_span + 2 * air_side
        else:
            # Single-ended
            if air_side is None:
                self.domain_width = 2 * max(self.w * 8, self.h * 15)
            else:
                self.domain_width = self.w + 2 * air_side

        self.boundaries = boundaries if boundaries else ["open", "open", "open", "gnd"]

        # Calculate physical coordinates
        self._calculate_coordinates()

        # Skin depth
        self.delta_s = np.sqrt(2 / (self.omega * mu_0 * self.sigma_cond))

        # Build geometry lists
        self.dielectrics, self.conductors = self._build_geometry_lists()

        # Create mesher and generate mesh
        mesher = Mesher(
            self.domain_width, self.domain_height,
            nx, ny, self.delta_s,
            self.conductors, self.dielectrics,
            symmetric=True
        )

        self.x, self.y = mesher.generate_mesh()
        self.dx_array = np.diff(self.x)
        self.dy_array = np.diff(self.y)

        # Setup geometry
        self._setup_geometry()

    def _calculate_coordinates(self):
        """Calculate all physical y-coordinates for geometry layers."""
        # Bottom extension for cut ground
        self.y_ext_start = self.t_gnd
        self.y_ext_end = self.t_gnd + self.gnd_cut_sub_h

        # New bottom ground plane location
        self.y_gnd_bot_start = self.y_ext_end
        self.y_gnd_bot_end = self.y_gnd_bot_start + self.t_gnd
        if self.gnd_cut_width == 0:
            self.y_gnd_bot_end = self.y_gnd_bot_start

        self.y_sub_start = self.y_gnd_bot_end
        self.y_sub_end = self.y_sub_start + self.h

        # Top dielectric
        self.y_top_diel_start = self.y_sub_end
        self.y_top_diel_end = self.y_top_diel_start + self.top_diel_h

        # Trace is embedded in top dielectric
        self.y_trace_start = self.y_top_diel_start
        self.y_trace_end = self.y_trace_start + self.t

        # Solder mask extents
        self.y_sm_sub_end = self.y_top_diel_end + self.sm_t_sub
        self.y_sm_trace_end = self.y_trace_end + self.sm_t_trace

        self.y_top_start = self.y_top_diel_end
        if self.use_sm:
            self.y_top_start = max(self.y_sm_sub_end, self.y_sm_trace_end)

        # Top air/dielectric region
        if self.air_top is None:
            self.top_dielectric_h = self.h * 15
            self.has_top_gnd = False
        else:
            self.top_dielectric_h = self.air_top + self.t
            self.has_top_gnd = (self.boundaries[2] == "gnd")

        self.y_top_end = self.y_top_start + self.top_dielectric_h

        if self.has_top_gnd:
            self.y_gnd_top_start = self.y_top_end
            self.y_gnd_top_end = self.y_gnd_top_start + self.t_gnd
            self.domain_height = self.y_gnd_top_end
        else:
            self.y_gnd_top_start = self.y_gnd_top_end = None
            self.domain_height = self.y_top_end

    def _build_geometry_lists(self):
        """Build lists of dielectrics and conductors from legacy parameters."""
        dielectrics = []
        conductors = []

        cx = self.domain_width / 2
        xl, xr = cx - self.w / 2, cx + self.w / 2

        # Substrate (covers both cutout extension and main substrate)
        if self.gnd_cut_sub_h > 0:
            # Substrate covers from extension start to substrate end
            dielectrics.append(Dielectric(
                0, self.y_ext_start,
                self.domain_width, self.y_trace_start,
                self.er, self.tan_delta
            ))
        else:
            # No cutout, just main substrate
            dielectrics.append(Dielectric(
                0, self.y_sub_start,
                self.domain_width, self.h,
                self.er, self.tan_delta
            ))

        # Top dielectric (if present)
        if self.top_diel_h > 0:
            dielectrics.append(Dielectric(
                0, self.y_top_diel_start,
                self.domain_width, self.top_diel_h,
                self.top_diel_er, self.top_diel_tand
            ))

        # Top air/dielectric region
        dielectrics.append(Dielectric(
            0, self.y_top_start,
            self.domain_width, self.top_dielectric_h,
            self.er_top, 0.0
        ))

        # Solder mask regions (overwrites previous)
        if self.use_sm:
            # Solder mask on substrate (full width, but will be overwritten by trace)
            dielectrics.append(Dielectric(
                0, self.y_top_diel_end,
                self.domain_width, self.sm_t_sub,
                self.sm_er, self.sm_tand
            ))

            # Solder mask on left side of trace
            xsl = xl - self.sm_t_side
            if xsl >= 0:
                dielectrics.append(Dielectric(
                    xsl, self.y_trace_start,
                    self.sm_t_side, self.t + self.sm_t_trace,
                    self.sm_er, self.sm_tand
                ))

            # Solder mask on right side of trace
            xsr = xr + self.sm_t_side
            if xsr <= self.domain_width:
                dielectrics.append(Dielectric(
                    xr, self.y_trace_start,
                    self.sm_t_side, self.t + self.sm_t_trace,
                    self.sm_er, self.sm_tand
                ))

            # Solder mask on top of trace
            dielectrics.append(Dielectric(
                xl, self.y_trace_end,
                self.w, self.sm_t_trace,
                self.sm_er, self.sm_tand
            ))

        # --- CONDUCTORS ---

        # Bottom ground (beneath everything)
        if self.t_gnd > 0:
            conductors.append(Conductor(
                0, 0,
                self.domain_width, self.t_gnd,
                is_signal=False
            ))

        # Bottom ground plane (above cutout extension)
        if self.gnd_cut_width == 0:
            # No cutout - full ground plane
            if self.y_gnd_bot_end > self.y_gnd_bot_start:
                conductors.append(Conductor(
                    0, self.y_gnd_bot_start,
                    self.domain_width, self.t_gnd,
                    is_signal=False
                ))
        else:
            # With cutout - ground on sides only
            cut_l = cx - self.gnd_cut_width / 2
            cut_r = cx + self.gnd_cut_width / 2

            # Left ground
            if cut_l > 0:
                conductors.append(Conductor(
                    0, self.y_gnd_bot_start,
                    cut_l, self.t_gnd,
                    is_signal=False
                ))

            # Right ground
            if cut_r < self.domain_width:
                conductors.append(Conductor(
                    cut_r, self.y_gnd_bot_start,
                    self.domain_width - cut_r, self.t_gnd,
                    is_signal=False
                ))

        # Signal trace(s)
        if self.is_differential:
            # Left trace (negative in odd mode)
            xl_left = cx - self.w - self.trace_spacing / 2
            conductors.append(Conductor(
                xl_left, self.y_trace_start,
                self.w, self.t,
                is_signal=True
            ))
            # Right trace (positive in odd mode)
            xl_right = cx + self.trace_spacing / 2
            conductors.append(Conductor(
                xl_right, self.y_trace_start,
                self.w, self.t,
                is_signal=True
            ))
        else:
            # Single trace
            conductors.append(Conductor(
                xl, self.y_trace_start,
                self.w, self.t,
                is_signal=True
            ))

        # Top ground plane (if present)
        if self.has_top_gnd:
            conductors.append(Conductor(
                0, self.y_gnd_top_start,
                self.domain_width, self.t_gnd,
                is_signal=False
            ))

        return dielectrics, conductors

    def _setup_geometry(self ):
        """Setup geometry from dielectric and conductor lists."""
        tol = 1e-11

        dielectrics = self.dielectrics
        conductors = self.conductors

        self.X, self.Y = np.meshgrid(self.x, self.y)

        # Initialize field arrays
        self.V = np.zeros(self.X.shape)
        self.epsilon_r = np.ones(self.X.shape)
        self.signal_mask = np.zeros(self.X.shape, dtype=bool)
        self.ground_mask = np.zeros(self.X.shape, dtype=bool)

        # For differential mode, track positive and negative traces separately
        if self.is_differential:
            self.signal_p_mask = np.zeros(self.X.shape, dtype=bool)
            self.signal_n_mask = np.zeros(self.X.shape, dtype=bool)

        # Initialize with air (epsilon_r = 1)
        self.epsilon_r = np.ones(self.X.shape)

        # Apply dielectrics (last overwrites)
        for diel in dielectrics:
            for i, yc in enumerate(self.y):
                if diel.y_min - tol <= yc <= diel.y_max + tol:
                    for j, xc in enumerate(self.x):
                        if diel.x_min - tol <= xc <= diel.x_max + tol:
                            self.epsilon_r[i, j] = diel.epsilon_r

        # Apply conductors
        signal_cond_idx = 0
        for cond in conductors:
            for i, yc in enumerate(self.y):
                if cond.y_min - tol <= yc <= cond.y_max + tol:
                    for j, xc in enumerate(self.x):
                        if cond.x_min - tol <= xc <= cond.x_max + tol:
                            if cond.is_signal:
                                self.signal_mask[i, j] = True
                                # For differential: first signal is negative, second is positive
                                if self.is_differential:
                                    if signal_cond_idx == 0:
                                        self.signal_n_mask[i, j] = True
                                    else:
                                        self.signal_p_mask[i, j] = True
                            else:
                                self.ground_mask[i, j] = True
                            self.V[i, j] = cond.voltage
            if cond.is_signal:
                signal_cond_idx += 1

        # Finalize conductor mask
        self.conductor_mask = self.signal_mask | self.ground_mask

    def solve_adaptive(self, **kwargs):
        """
        Override solve_adaptive to handle both single-ended and differential modes.

        For single-ended: Returns Z0, eps_eff, C, C0
        For differential: Returns differential_results dict with odd/even mode parameters
        """
        if not self.is_differential:
            # Single-ended mode
            return super().solve_adaptive(**kwargs)
        else:
            # Differential mode
            return self._solve_adaptive_differential(**kwargs)

    def _solve_adaptive_differential(self,
                                      max_iters=10,
                                      refine_frac=0.15,
                                      energy_tol=0.05,
                                      param_tol=0.1,
                                      max_nodes=20000,
                                      min_converged_passes=2):
        """
        Adaptive mesh solve for differential mode with odd and even modes.

        Returns: dict with differential parameters
        """
        prev_energy_odd = None
        prev_energy_even = None
        prev_Z_odd = None
        prev_Z_even = None
        converged_count = 0

        for it in range(max_iters):
            # Solve odd mode (+1V and -1V)
            Z_odd, eps_eff_odd, C_odd, C0_odd, Ex_odd, Ey_odd = self._solve_single_mode(
                'odd', vacuum_first=True
            )

            # Solve even mode (+1V and +1V)
            Z_even, eps_eff_even, C_even, C0_even, Ex_even, Ey_even = self._solve_single_mode(
                'even', vacuum_first=True
            )

            # Energy-based error
            energy_odd, energy_err_odd = self._compute_energy_error(Ex_odd, Ey_odd, prev_energy_odd)
            energy_even, energy_err_even = self._compute_energy_error(Ex_even, Ey_even, prev_energy_even)
            energy_err = max(energy_err_odd, energy_err_even)

            # Parameter-based error
            if prev_Z_odd is not None:
                z_odd_err = abs(Z_odd - prev_Z_odd) / max(abs(prev_Z_odd), 1e-12)
                z_even_err = abs(Z_even - prev_Z_even) / max(abs(prev_Z_even), 1e-12)
                param_err = max(z_odd_err, z_even_err)
            else:
                param_err = 1.0

            print(f"Pass {it+1}: Energy err={energy_err:.3g}, Param err={param_err:.3g}, "
                  f"Nodes={self.nx}x{self.ny}")

            # Check convergence
            if prev_Z_odd is not None:
                if energy_err < energy_tol and param_err < param_tol:
                    converged_count += 1
                    if converged_count >= min_converged_passes:
                        print(f"Converged after {it+1} passes")
                        break
                else:
                    converged_count = 0

            prev_energy_odd = energy_odd
            prev_energy_even = energy_even
            prev_Z_odd = Z_odd
            prev_Z_even = Z_even

            # Node budget check
            if self.nx * self.ny > max_nodes:
                print("Node budget reached")
                break

            # Refine mesh (use odd mode fields for refinement)
            if it != max_iters - 1:
                self.refine_mesh(Ex_odd, Ey_odd, frac=refine_frac)
                self._setup_geometry()

        # Calculate losses
        alpha_c_odd = self._calculate_differential_conductor_loss(Ex_odd, Ey_odd, Z_odd, 'odd')
        alpha_d_odd = self._calculate_differential_dielectric_loss(Ex_odd, Ey_odd, Z_odd)

        alpha_c_even = self._calculate_differential_conductor_loss(Ex_even, Ey_even, Z_even, 'even')
        alpha_d_even = self._calculate_differential_dielectric_loss(Ex_even, Ey_even, Z_even)

        # Store fields as class members (as lists for differential)
        self.Ex = [Ex_odd, Ex_even]
        self.Ey = [Ey_odd, Ey_even]

        # Differential and common mode impedances
        Z_diff = 2 * Z_odd
        Z_common = Z_even / 2

        return {
            'Z_odd': Z_odd,
            'Z_even': Z_even,
            'Z_diff': Z_diff,
            'Z_common': Z_common,
            'eps_eff_odd': eps_eff_odd,
            'eps_eff_even': eps_eff_even,
            'C_odd': C_odd,
            'C_even': C_even,
            'alpha_c_odd': alpha_c_odd,
            'alpha_c_even': alpha_c_even,
            'alpha_d_odd': alpha_d_odd,
            'alpha_d_even': alpha_d_even,
            'alpha_total_odd': alpha_c_odd + alpha_d_odd,
            'alpha_total_even': alpha_c_even + alpha_d_even
        }

    def _solve_single_mode(self, mode, vacuum_first=True):
        """
        Solve a single mode (odd or even) for differential pair.

        Parameters:
        -----------
        mode : str
            'odd' or 'even'
        vacuum_first : bool
            Whether to solve vacuum case first for C0 calculation

        Returns:
        --------
        Z0, eps_eff, C, C0, Ex, Ey
        """
        # Set voltages based on mode
        self.V = np.zeros(self.X.shape)
        if mode == 'odd':
            self.V[self.signal_n_mask] = -1.0  # Left trace: -1V
            self.V[self.signal_p_mask] = 1.0   # Right trace: +1V
        else:  # even
            self.V[self.signal_n_mask] = 1.0   # Left trace: +1V
            self.V[self.signal_p_mask] = 1.0   # Right trace: +1V
        self.V[self.ground_mask] = 0.0

        # Calculate C0 (vacuum capacitance) if requested
        if vacuum_first:
            self.solve_laplace(vacuum=True)
            # Use only positive trace for charge calculation
            orig_signal_mask = self.signal_mask.copy()
            self.signal_mask = self.signal_p_mask.copy()
            C0 = self.calculate_capacitance(vacuum=True)
            self.signal_mask = orig_signal_mask
        else:
            C0 = None

        # Solve with dielectric
        self.solve_laplace()
        # Use only positive trace for charge calculation
        orig_signal_mask = self.signal_mask.copy()
        self.signal_mask = self.signal_p_mask.copy()
        C = self.calculate_capacitance()
        self.signal_mask = orig_signal_mask

        # Calculate fields
        Ex, Ey = self.compute_fields()

        # Calculate impedance
        if C0 is not None:
            eps_eff = C / C0
            Z0 = 1 / (c * np.sqrt(C * C0))
        else:
            eps_eff = None
            Z0 = None

        return Z0, eps_eff, C, C0, Ex, Ey

    def _calculate_differential_conductor_loss(self, Ex, Ey, Z0, mode):
        """Calculate conductor loss for differential pair."""
        Rs = np.sqrt(self.omega * mu_0 / (2 * self.sigma_cond))
        delta = np.sqrt(2 / (self.omega * mu_0 * self.sigma_cond))

        def get_dx(j):
            return self.dx_array[j] if 0 <= j < len(self.dx_array) else self.dx_array[-1]

        def get_dy(i):
            return self.dy_array[i] if 0 <= i < len(self.dy_array) else self.dy_array[-1]

        ny, nx = self.ny, self.nx
        Pc = 0.0

        for i in range(1, ny - 1):
            for j in range(1, nx - 1):
                if not (self.signal_p_mask[i, j] or self.signal_n_mask[i, j] or self.ground_mask[i, j]):
                    continue

                neighbors = [
                    (i, j+1, 'r', get_dy(i)),
                    (i, j-1, 'l', get_dy(i)),
                    (i+1, j, 'u', get_dx(j)),
                    (i-1, j, 'd', get_dx(j)),
                ]

                cell_K_sq = 0.0
                cell_dl = 0.0

                for ni, nj, direction, dl in neighbors:
                    if (self.signal_p_mask[ni, nj] or self.signal_n_mask[ni, nj] or
                        self.ground_mask[ni, nj]):
                        continue

                    eps_diel = self.epsilon_r[ni, nj]
                    Ex_diel = Ex[ni, nj]
                    Ey_diel = Ey[ni, nj]

                    if direction == 'r':
                        E_norm = Ex_diel
                    elif direction == 'l':
                        E_norm = -Ex_diel
                    elif direction == 'u':
                        E_norm = Ey_diel
                    elif direction == 'd':
                        E_norm = -Ey_diel

                    Z0_freespace = np.sqrt(mu_0 / epsilon_0)
                    H_tan = abs(E_norm) * np.sqrt(eps_diel) / Z0_freespace
                    K = H_tan

                    cell_K_sq += K**2 * dl
                    cell_dl += dl

                if cell_dl > 0:
                    dP = 0.5 * Rs * cell_K_sq
                    Pc += dP

        # Power normalization for 1W
        Pc_1W = 0.5 * Pc * 2 * Z0
        alpha_db_per_m = 8.686 * Pc_1W / 2.0

        return alpha_db_per_m

    def _calculate_differential_dielectric_loss(self, Ex, Ey, Z0):
        """Calculate dielectric loss for differential pair."""
        Pd = 0.0
        for i in range(self.ny - 1):
            for j in range(self.nx - 1):
                if (self.signal_p_mask[i, j] or self.signal_n_mask[i, j] or
                    self.ground_mask[i, j]):
                    continue
                if self.epsilon_r[i, j] <= 1.01:
                    continue

                E2 = Ex[i, j]**2 + Ey[i, j]**2
                dA = self.dx_array[j] * self.dy_array[i]
                Pd += 0.5 * self.omega * epsilon_0 * self.epsilon_r[i, j] * self.tan_delta * E2 * dA

        # Power flow for 1W in the transmission line
        P_flow = 1.0 / (2 * Z0)

        # Attenuation constant
        return 8.686 * (0.5 * Pd / (2 * P_flow))


# ============================================================================
# TEST FUNCTION
# ============================================================================

def solve_microstrip(spacing=0, plots=True):
    print(f"\nSolving {'Differential' if spacing > 0 else 'Single-Ended'} Microstrip (v2 with Mesher)...")
    solver = MicrostripSolver2D(
        substrate_height=1.6e-3,
        trace_width=3e-3,
        trace_thickness=35e-6,
        gnd_thickness=35e-6,
        epsilon_r=4.5,
        tan_delta=0.02,
        sigma_cond=5.8e7,
        freq=1e9,
        nx=10, ny=10,
        use_sm=False,
        trace_spacing=spacing if spacing > 0 else None,
        boundaries=["open", "open", "open", "gnd"]
    )

    result = solver.solve_adaptive(param_tol=0.01)

    # Handle both single-ended and differential results
    if isinstance(result, dict):
        # Differential mode - extract results
        print(f"\n{'='*60}")
        print(f"DIFFERENTIAL MICROSTRIP ANALYSIS RESULTS")
        print(f"{'='*60}")
        print(f"Differential Impedance Z_diff: {result['Z_diff']:.2f} Ω  (2 × Z_odd)")
        print(f"Common-Mode Impedance Z_common: {result['Z_common']:.2f} Ω  (Z_even / 2)")
        print(f"\nModal Impedances:")
        print(f"  Odd-Mode  Z_odd:  {result['Z_odd']:.2f} Ω  (εᵣₑff = {result['eps_eff_odd']:.3f})")
        print(f"  Even-Mode Z_even: {result['Z_even']:.2f} Ω  (εᵣₑff = {result['eps_eff_even']:.3f})")
        print(f"\nLosses @ 1GHz:")
        print(f"  Odd-Mode:  Diel={result['alpha_d_odd']:.4f} dB/m, Cond={result['alpha_c_odd']:.4f} dB/m, Total={result['alpha_total_odd']:.4f} dB/m")
        print(f"  Even-Mode: Diel={result['alpha_d_even']:.4f} dB/m, Cond={result['alpha_c_even']:.4f} dB/m, Total={result['alpha_total_even']:.4f} dB/m")
        print(f"{'='*60}\n")

        # For plotting, use odd-mode fields
        if plots:
            solver.plot_geometry()
            solver.plot(solver.Ex[0], solver.Ey[0])  # Odd mode
            solver.plot(solver.Ex[1], solver.Ey[1])  # Odd mode
    else:
        # Single-ended mode
        Z0, eps_eff, C, C0 = result

        alpha_cond, J = solver.calculate_conductor_loss(solver.Ex, solver.Ey, Z0)
        alpha_diel = solver.calculate_dielectric_loss(solver.Ex, solver.Ey, Z0)
        alpha_total = alpha_cond + alpha_diel

        z, rlgc, eps_eff_mode = solver.rlgc(alpha_cond, alpha_diel, C, Z0)
        print(f"Z (complex) = {z:.2f} ohm, eps_eff {eps_eff_mode:.3f}, RLGC {rlgc}")

        print(f"\n{'='*60}")
        print(f"SINGLE-ENDED MICROSTRIP ANALYSIS RESULTS")
        print(f"{'='*60}")
        print(f"Characteristic Impedance Z0:  {Z0:.2f} Ω")
        print(f"Effective Permittivity εᵣₑff: {eps_eff:.3f}")
        print(f"Losses (dB/m) @ 1GHz:         Diel={alpha_diel:.4f}, Cond={alpha_cond:.4g}")
        print(f"Total Attenuation:            {alpha_total:.4f} dB/m")
        print(f"{'='*60}\n")

        if plots:
            solver.plot_geometry()
            solver.plot(solver.Ex, solver.Ey)


if __name__ == "__main__":
    plots = False
    if len(sys.argv) > 1:
        if "plot" in sys.argv[1]:
            plots = True
    solve_microstrip(0, plots)
    solve_microstrip(0.5e-3, plots)
