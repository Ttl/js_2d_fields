#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import epsilon_0, mu_0, c, pi
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve, bicgstab, spilu, LinearOperator, cg, norm
from numpy import pi
import matplotlib.patches as mpatches

mu_0 = 4 * pi * 1e-7

class FieldSolver2D:
    def __init__(self):
        raise ValueError("Call geometry specific initializer")

    def _smooth_transition(self, start, end, n_points, curve_end='end', beta=4.0):
        if n_points <= 1: return np.array([start, end])
        xi = np.linspace(0, 1, n_points)
        if curve_end == 'end':
            eta = np.tanh(beta * xi) / np.tanh(beta)
        elif curve_end == 'both':
            eta = (np.tanh(beta * (xi - 0.5)) / np.tanh(beta * 0.5) + 1) / 2
        else:
            eta = 1 - np.tanh(beta * (1 - xi)) / np.tanh(beta)
        return start + eta * (end - start)

    def solve_laplace(self, vacuum=False):
        self.ny = len(self.y)
        self.nx = len(self.x)
        self.dx_array = np.diff(self.x)
        self.dy_array = np.diff(self.y)

        N = self.nx * self.ny

        # Create mapping: unknown nodes only
        is_unknown = ~self.conductor_mask.flatten()
        unknown_indices = np.where(is_unknown)[0]
        N_unknown = len(unknown_indices)

        # Map from full index to reduced index
        full_to_reduced = -np.ones(N, dtype=int)
        full_to_reduced[unknown_indices] = np.arange(N_unknown)

        A = lil_matrix((N_unknown, N_unknown))
        b = np.zeros(N_unknown)

        def idx(i, j):
            return i * self.nx + j

        for i in range(self.ny):
            for j in range(self.nx):
                full_n = idx(i, j)

                # Skip conductor nodes entirely
                if self.conductor_mask[i, j]:
                    continue

                n = full_to_reduced[full_n]  # Reduced system index

                is_boundary = (i == 0 or i == self.ny - 1 or j == 0 or j == self.nx - 1)

                if is_boundary:
                    dxr = self.dx_array[j] if j < self.nx - 1 else self.dx_array[j - 1]
                    dxl = self.dx_array[j - 1] if j > 0 else self.dx_array[j]
                    dyu = self.dy_array[i] if i < self.ny - 1 else self.dy_array[i - 1]
                    dyd = self.dy_array[i - 1] if i > 0 else self.dy_array[i]
                else:
                    dxr = self.dx_array[j]
                    dxl = self.dx_array[j - 1]
                    dyu = self.dy_array[i]
                    dyd = self.dy_array[i - 1]

                if vacuum:
                    err = erl = eru = erd = 1
                else:
                    erc = self.epsilon_r[i, j]
                    if is_boundary:
                        err = 0.5 * (erc + self.epsilon_r[i, min(j + 1, self.nx - 1)])
                        erl = 0.5 * (erc + self.epsilon_r[i, max(j - 1, 0)])
                        eru = 0.5 * (erc + self.epsilon_r[min(i + 1, self.ny - 1), j])
                        erd = 0.5 * (erc + self.epsilon_r[max(i - 1, 0), j])
                    else:
                        err = erc if self.conductor_mask[i, j + 1] else 0.5 * (erc + self.epsilon_r[i, j + 1])
                        erl = erc if self.conductor_mask[i, j - 1] else 0.5 * (erc + self.epsilon_r[i, j - 1])
                        eru = erc if self.conductor_mask[i + 1, j] else 0.5 * (erc + self.epsilon_r[i + 1, j])
                        erd = erc if self.conductor_mask[i - 1, j] else 0.5 * (erc + self.epsilon_r[i - 1, j])

                area_i = 0.5 * (dyd + dyu)
                area_j = 0.5 * (dxl + dxr)

                if is_boundary:
                    cr = -err * area_i / dxr if j < self.nx - 1 else 0
                    cl = -erl * area_i / dxl if j > 0 else 0
                    cu = -eru * area_j / dyu if i < self.ny - 1 else 0
                    cd = -erd * area_j / dyd if i > 0 else 0
                else:
                    cr = -err * area_i / dxr
                    cl = -erl * area_i / dxl
                    cu = -eru * area_j / dyu
                    cd = -erd * area_j / dyd

                A[n, n] = -(cr + cl + cu + cd)

                # Right
                if j < self.nx - 1:
                    if not self.conductor_mask[i, j + 1]:
                        n_right = full_to_reduced[idx(i, j + 1)]
                        A[n, n_right] = cr
                    else:
                        b[n] -= cr * self.V[i, j + 1]

                # Left
                if j > 0:
                    if not self.conductor_mask[i, j - 1]:
                        n_left = full_to_reduced[idx(i, j - 1)]
                        A[n, n_left] = cl
                    else:
                        b[n] -= cl * self.V[i, j - 1]

                # Up
                if i < self.ny - 1:
                    if not self.conductor_mask[i + 1, j]:
                        n_up = full_to_reduced[idx(i + 1, j)]
                        A[n, n_up] = cu
                    else:
                        b[n] -= cu * self.V[i + 1, j]

                # Down
                if i > 0:
                    if not self.conductor_mask[i - 1, j]:
                        n_down = full_to_reduced[idx(i - 1, j)]
                        A[n, n_down] = cd
                    else:
                        b[n] -= cd * self.V[i - 1, j]

        A_csr = A.tocsr()
        # A is symmetric positive definite

        # Solve reduced system
        if 0:
            diag = A_csr.diagonal()
            M = LinearOperator(A_csr.shape, lambda x: x / diag)
            solution_reduced, info = cg(A_csr, b, rtol=1e-6, M=M)
        elif 0:
            ilu = spilu(A_csr, drop_tol=1e-4, fill_factor=10)

            # Create a LinearOperator that applies the inverse of the ILU
            M = LinearOperator(A_csr.shape, ilu.solve)

            # 3. Solve with Preconditioner (M)
            solution_reduced, info = bicgstab(A_csr, b, rtol=1e-7, maxiter=1000, M=M)

            if info > 0:
                print(f"Warning: Convergence not reached after {info} iterations.")
            elif info < 0:
                print(f"Error: Illegal input or breakdown (info={info}).")
            else:
                print("Success: bicgstab converged.")
        else:
            solution_reduced = spsolve(A_csr, b)

        # Reconstruct full solution
        solution_full = np.zeros(N)
        solution_full[unknown_indices] = solution_reduced
        solution_full[~is_unknown] = self.V.flatten()[~is_unknown]

        self.V = solution_full.reshape((self.ny, self.nx))

    def compute_fields(self):
        Ex = np.zeros_like(self.V)
        Ey = np.zeros_like(self.V)

        for i in range(1, self.ny - 1):
            for j in range(1, self.nx - 1):
                if self.conductor_mask[i, j]:
                    continue

                dxl = self.dx_array[j-1]
                dxr = self.dx_array[j]
                dyd = self.dy_array[i-1]
                dyu = self.dy_array[i]

                Ex[i, j] = -(
                     (dxl / (dxr * (dxl + dxr))) * self.V[i, j + 1] +
                     ((dxr - dxl) / (dxl * dxr)) * self.V[i, j] -
                     (dxr / (dxl * (dxl + dxr))) * self.V[i, j - 1]
                )

                Ey[i, j] = -(
                     (dyd / (dyu * (dyd + dyu))) * self.V[i + 1, j] +
                     ((dyu - dyd) / (dyd * dyu)) * self.V[i, j] -
                     (dyu / (dyd * (dyd + dyu))) * self.V[i - 1, j]
                )
        return Ex, Ey

    def calculate_capacitance(self, vacuum=False):
        Q = 0.0
        def get_dx(j): return self.dx_array[j] if j < len(self.dx_array) else self.dx_array[-1]
        def get_dy(i): return self.dy_array[i] if i < len(self.dy_array) else self.dy_array[-1]

        for i in range(1, self.ny - 1):
            for j in range(1, self.nx - 1):
                if not self.signal_mask[i, j]:
                    continue

                # Check neighbors. If dielectric, calculate flux.
                if not self.signal_mask[i, j+1]: # Right
                    En = (self.V[i, j] - self.V[i, j+1]) / get_dx(j)
                    er = 1 if vacuum else self.epsilon_r[i, j+1]
                    Q += epsilon_0 * er * En * ((get_dy(i-1) + get_dy(i))/2)

                if not self.signal_mask[i, j-1]: # Left
                    En = (self.V[i, j] - self.V[i, j-1]) / get_dx(j-1)
                    er = 1 if vacuum else self.epsilon_r[i, j-1]
                    Q += epsilon_0 * er * En * ((get_dy(i-1) + get_dy(i))/2)

                if not self.signal_mask[i+1, j]: # Top
                    En = (self.V[i, j] - self.V[i+1, j]) / get_dy(i)
                    er = 1 if vacuum else self.epsilon_r[i+1, j]
                    Q += epsilon_0 * er * En * ((get_dx(j-1) + get_dx(j))/2)

                if not self.signal_mask[i-1, j]: # Bottom
                    En = (self.V[i, j] - self.V[i-1, j]) / get_dy(i-1)
                    er = 1 if vacuum else self.epsilon_r[i-1, j]
                    Q += epsilon_0 * er * En * ((get_dx(j-1) + get_dx(j))/2)
        return abs(Q)

    def calculate_parameters(self):
        """
        Calculate Z0, eps_eff, C and C0 (vacuum capacitance
        """
        self.solve_laplace(vacuum=True)
        C0 = self.calculate_capacitance(vacuum=True)

        self.solve_laplace()
        C = self.calculate_capacitance()

        eps_eff = C / C0
        Z0 = 1 / (c * np.sqrt(C * C0))
        return Z0, eps_eff, C, C0

    def calculate_conductor_loss(self, Ex, Ey, Z0):
        """
        Conductor loss using perturbation method.

        For quasi-TEM mode:
        - Surface current K = n × H_tangential
        - H from Ampere's law: ∇ × H = jωD = jωε₀εᵣE
        - Power loss: P = (1/2) Rs ∫|K|² dl
        """

        Rs = np.sqrt(self.omega * mu_0 / (2 * self.sigma_cond))
        delta = np.sqrt(2 / (self.omega * mu_0 * self.sigma_cond))

        def get_dx(j):
            return self.dx_array[j] if 0 <= j < len(self.dx_array) else self.dx_array[-1]

        def get_dy(i):
            return self.dy_array[i] if 0 <= i < len(self.dy_array) else self.dy_array[-1]

        ny, nx = self.ny, self.nx

        # Build current density map for visualization
        J_map = np.zeros((ny, nx), dtype=np.complex128)

        # Calculate surface current and power loss
        Pc = 0.0

        for i in range(1, ny - 1):
            for j in range(1, nx - 1):
                if not (self.signal_mask[i, j] or self.ground_mask[i, j]):
                    continue

                # Check each face for dielectric neighbor
                neighbors = [
                    (i, j+1, 'r', get_dy(i)),
                    (i, j-1, 'l', get_dy(i)),
                    (i+1, j, 'u', get_dx(j)),
                    (i-1, j, 'd', get_dx(j)),
                ]

                cell_K_sq = 0.0
                cell_dl = 0.0

                for ni, nj, direction, dl in neighbors:
                    if self.signal_mask[ni, nj] or self.ground_mask[ni, nj]:
                        continue  # Neighbor is also conductor

                    # Get E-field in dielectric cell
                    eps_diel = self.epsilon_r[ni, nj]
                    Ex_diel = Ex[ni, nj]
                    Ey_diel = Ey[ni, nj]

                    # Normal E-field component (pointing from metal into dielectric)
                    if direction == 'r':
                        E_norm = Ex_diel
                    elif direction == 'l':
                        E_norm = -Ex_diel
                    elif direction == 'u':
                        E_norm = Ey_diel
                    elif direction == 'd':
                        E_norm = -Ey_diel

                    # For quasi-TEM transmission line:
                    # Surface current density K relates to tangential H-field
                    # From boundary condition: K = H_tangential
                    # From Maxwell: ∇ × H = jωε₀εᵣE
                    # For TEM wave: H_tan/E_norm = Y₀√εᵣ where Y₀ = 1/Z₀_freespace
                    # So: H_tan = E_norm × √(ε₀εᵣ/μ₀) = E_norm × √εᵣ / Z₀_freespace

                    Z0_freespace = np.sqrt(mu_0 / epsilon_0)  # ~377 Ω
                    H_tan = abs(E_norm) * np.sqrt(eps_diel) / Z0_freespace

                    # Surface current K = H_tan (A/m)
                    K = H_tan

                    # Accumulate for power calculation
                    cell_K_sq += K**2 * dl
                    cell_dl += dl

                if cell_dl > 0:
                    # Power dissipation: P = (1/2) Rs ∫|K|² dl
                    dP = 0.5 * Rs * cell_K_sq
                    Pc += dP

                    # For visualization: approximate volume current density
                    K_avg = np.sqrt(cell_K_sq / cell_dl)
                    J_surface = K_avg / delta
                    J_map[i, j] = J_surface

        # Scale to 1W power flow
        Pc_1W = Pc * 2 * Z0
        J_map_1W = J_map * np.sqrt(2 * Z0)

        alpha_db_per_m = 8.686 * Pc_1W / 2.0

        # Store for visualization
        metal_cells = np.argwhere(self.signal_mask | self.ground_mask)
        J_viz = np.array([J_map_1W[int(i), int(j)] for i, j in metal_cells])

        self._last_J_solution = J_viz
        self._last_J_cells = metal_cells

        #print(f"Power loss (1W):    {Pc_1W:.6e} W")
        #print(f"Max |J| (1W):       {np.max(np.abs(J_viz)):.3e} A/m²")
        #print(f"Conductor loss:     {alpha_db_per_m:.4f} dB/m\n")

        return alpha_db_per_m, J_viz

    def calculate_dielectric_loss(self, Ex, Ey, Z0):
        Pd = 0.0
        for i in range(self.ny - 1):
            for j in range(self.nx - 1):
                if self.conductor_mask[i, j]:
                    continue
                if self.epsilon_r[i, j] <= 1.01:
                    continue

                E2 = Ex[i, j]**2 + Ey[i, j]**2
                dA = self.dx_array[j] * self.dy_array[i]
                Pd += 0.5 * self.omega * epsilon_0 * self.epsilon_r[i, j] * self.tan_delta * E2 * dA

        P_flow = 1.0 / (2 * Z0)
        return 8.686 * (Pd / (2 * P_flow))

    def plot(self, Ex, Ey):
        fig, axs = plt.subplots(1, 2, figsize=(12, 5))

        E = np.sqrt(Ex**2 + Ey**2)
        im1 = axs[0].contourf(self.X * 1e3, self.Y * 1e3, E, 10, cmap='viridis')
        axs[0].set_title("|E| Field")
        axs[0].set_xlabel("x (mm)")
        axs[0].set_ylabel("y (mm)")
        axs[0].set_aspect('equal')
        plt.colorbar(im1, ax=axs[0], label="Log(V/m)")

        im2 = axs[1].contourf(self.X * 1e3, self.Y * 1e3, self.V, 10, cmap='viridis')
        mask = self.conductor_mask
        axs[1].scatter(self.X[mask] * 1e3, self.Y[mask] * 1e3, s=1, c='black', alpha=0.5)
        axs[1].set_title("Potential V")
        axs[1].set_xlabel("x (mm)")
        axs[1].set_aspect('equal')
        plt.colorbar(im2, ax=axs[1], label="Volts")

        plt.tight_layout()
        plt.show()

    def plot_geometry(self, unit='um', show_mesh=True, mesh_stride=1,
                      mesh_alpha=0.5, mesh_color='black', mesh_linewidth=0.3):
        """
        Plots the permittivity and conductors with optional mesh overlay.

        Parameters:
        -----------
        unit : str
            Unit for display ('um', 'mm', 'm')
        show_mesh : bool
            Whether to overlay mesh grid lines
        mesh_stride : int or tuple
            Stride for mesh lines. If int, same for both directions.
            If tuple (stride_x, stride_y). If None, auto-calculated.
        mesh_alpha : float
            Transparency of mesh lines
        mesh_color : str
            Color of mesh lines
        mesh_linewidth : float
            Width of mesh lines
        """
        scales = {'um': 1e6, 'mm': 1e3, 'm': 1.0}
        scale = scales.get(unit, 1e6)

        fig, ax = plt.subplots(figsize=(12, 8))

        # 1. Create a display array for permittivity
        plot_data = np.copy(self.epsilon_r)

        # 2. Plot the dielectric background
        im = ax.pcolormesh(self.x * scale, self.y * scale, plot_data,
                           shading='auto', cmap='cool', alpha=0.8)

        # 3. Overlay the Conductors
        cond_display = np.full(self.epsilon_r.shape, np.nan)
        cond_display[self.ground_mask] = 0  # Index for ground color
        cond_display[self.signal_mask] = 1  # Index for signal color

        # Define a specific colormap for conductors
        cond_cmap = plt.cm.colors.ListedColormap(['#404040', '#FF8C00']) # Dark Gray and Orange
        ax.pcolormesh(self.x * scale, self.y * scale, cond_display,
                      shading='auto', cmap=cond_cmap)

        # 4. Overlay mesh grid lines if requested
        if show_mesh:
            # Auto-calculate stride if not provided
            if mesh_stride is None:
                # Show at most 100 lines in each direction
                stride_x = max(1, len(self.x) // 100)
                stride_y = max(1, len(self.y) // 100)
            elif isinstance(mesh_stride, (int, float)):
                stride_x = stride_y = int(mesh_stride)
            else:
                stride_x, stride_y = mesh_stride

            # Plot vertical grid lines
            for i in range(0, len(self.x), stride_x):
                ax.plot([self.x[i] * scale, self.x[i] * scale],
                        [self.y[0] * scale, self.y[-1] * scale],
                        color=mesh_color, alpha=mesh_alpha, linewidth=mesh_linewidth,
                        zorder=10)

            # Plot horizontal grid lines
            for j in range(0, len(self.y), stride_y):
                ax.plot([self.x[0] * scale, self.x[-1] * scale],
                        [self.y[j] * scale, self.y[j] * scale],
                        color=mesh_color, alpha=mesh_alpha, linewidth=mesh_linewidth,
                        zorder=10)

        # 5. Formatting
        ax.set_aspect('equal')
        title = f"2D Field Solver Geometry (Conductors & Permittivity)\nScale: {unit}"
        if show_mesh:
            title += f"\nMesh: {len(self.x)} × {len(self.y)} points"
        ax.set_title(title)
        ax.set_xlabel(f"X [{unit}]")
        ax.set_ylabel(f"Y [{unit}]")

        # Add Colorbar for Permittivity
        cbar = fig.colorbar(im, ax=ax, label='Relative Permittivity ($\epsilon_r$)')

        # Add custom Legend for conductors
        gnd_patch = mpatches.Patch(color='#404040', label='Ground', linewidth=0)
        sig_patch = mpatches.Patch(color='#FF8C00', label='Signal Trace', linewidth=0)
        ax.legend(handles=[gnd_patch, sig_patch], loc='upper right')

        plt.tight_layout()
        plt.show()

    def rlgc(self, alpha_cond, alpha_diel, C_mode, Z0_mode):
        # dB/m to Np/m
        alpha_c = alpha_cond / 8.686
        alpha_d = alpha_diel / 8.686

        R = 2 * Z0_mode * alpha_c
        G = 2 * alpha_d / Z0_mode

        L = (Z0_mode**2) * C_mode

        # Modal epsilon_eff from LC
        eps_eff_mode = (L * C_mode) / (mu_0 * epsilon_0)

        Zc = np.sqrt((R + 1j * self.omega * L) /
                     (G + 1j * self.omega * C_mode))

        rlgc = {"R": float(R), "L": float(L), "G": float(G), "C": float(C_mode)}
        return Zc, rlgc, eps_eff_mode

    def modal_to_physical_rlgc(self, results):
        """
        Convert modal (odd/even) parameters to physical 2x2 RLGC matrices
        """
        # Modal parameters
        R_odd = 2 * results['Z_odd'] * (results['alpha_c_odd'] / 8.686)
        R_even = 2 * results['Z_even'] * (results['alpha_c_even'] / 8.686)

        L_odd = results['Z_odd']**2 * results['C_odd']
        L_even = results['Z_even']**2 * results['C_even']

        G_odd = 2 * (results['alpha_d_odd'] / 8.686) / results['Z_odd']
        G_even = 2 * (results['alpha_d_even'] / 8.686) / results['Z_even']

        C_odd = results['C_odd']
        C_even = results['C_even']

        # Transform to physical domain
        # For L and C: coupling terms are NEGATIVE
        R11 = R22 = (R_odd + R_even) / 2
        L11 = L22 = (L_odd + L_even) / 2
        G11 = G22 = (G_odd + G_even) / 2
        C11 = C22 = (C_odd + C_even) / 2

        R12 = R21 = (R_even - R_odd) / 2
        L12 = L21 = (L_even - L_odd) / 2
        G12 = G21 = (G_even - G_odd) / 2
        C12 = C21 = (C_even - C_odd) / 2

        return {
            'R': np.array([[R11, R12], [R21, R22]]),
            'L': np.array([[L11*1e9, L12*1e9], [L21*1e9, L22*1e9]]),
            'G': np.array([[G11*1e3, G12*1e3], [G21*1e3, G22*1e3]]),
            'C': np.array([[C11*1e12, C12*1e12], [C21*1e12, C22*1e12]])
        }


class DifferentialMicrostripSolver2D(FieldSolver2D):
    """
    Differential microstrip solver for calculating odd and even mode impedances.

    Odd mode: Traces driven with opposite polarity (+V and -V)
    Even mode: Traces driven with same polarity (+V and +V)

    Z_diff = 2 * Z_odd
    Z_common = Z_even / 2
    """

    def __init__(self, substrate_height, trace_width, trace_spacing, trace_thickness,
                 gnd_thickness=35e-6, epsilon_r=4.5, tan_delta=0.02,
                 sigma_diel=0.0, sigma_cond=5.8e7, epsilon_r_top=1,
                 air_top=None, air_side=None, freq=1e9,
                 nx=400, ny=300, skin_cells=100, boundaries=None,
                 use_sm=False, sm_t_sub=20e-6, sm_t_trace=20e-6,
                 sm_t_side=20e-6, sm_er=3.5, sm_tand=0.02):
        """
        Initialize differential microstrip geometry.

        Parameters:
        -----------
        substrate_height : float
            PCB substrate thickness (m)
        trace_width : float
            Width of each trace (m)
        trace_spacing : float
            Edge-to-edge spacing between traces (m)
        trace_thickness : float
            Copper thickness (m)
        """
        self.h = substrate_height
        self.w = trace_width
        self.s = trace_spacing  # Edge-to-edge spacing
        self.t = trace_thickness
        self.t_gnd = gnd_thickness
        self.er = epsilon_r
        self.er_top = epsilon_r_top
        self.tan_delta = tan_delta
        self.sigma_diel = sigma_diel
        self.sigma_cond = sigma_cond

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
        self.skin_cells = skin_cells

        # Domain width needs to accommodate both traces
        trace_span = 2 * self.w + self.s
        if air_side is None:
            self.domain_width = 2 * max(trace_span * 4, self.h * 15)
        else:
            self.domain_width = trace_span + 2 * air_side

        self.boundaries = boundaries if boundaries else ["open", "open", "open", "gnd"]

        # --- Physical Y-Coordinates ---
        self.y_gnd_bot_start = 0.0
        self.y_gnd_bot_end = self.t_gnd
        self.y_sub_start = self.y_gnd_bot_end
        self.y_sub_end = self.y_sub_start + self.h
        self.y_trace_start = self.y_sub_end
        self.y_trace_end = self.y_trace_start + self.t

        # Solder Mask Boundaries
        self.y_sm_sub_end = self.y_sub_end + self.sm_t_sub
        self.y_sm_trace_end = self.y_trace_end + self.sm_t_trace
        self.y_top_start = self.y_trace_end

        if self.use_sm:
            self.y_top_start = max(self.y_sm_sub_end, self.y_sm_trace_end)

        if air_top is None:
            self.top_dielectric_h = self.h * 15
            self.has_top_gnd = False
        else:
            self.top_dielectric_h = air_top
            self.has_top_gnd = (self.boundaries[2] == "gnd")

        self.y_top_end = self.y_top_start + self.top_dielectric_h

        if self.has_top_gnd:
            self.y_gnd_top_start = self.y_top_end
            self.y_gnd_top_end = self.y_gnd_top_start + self.t_gnd
            self.domain_height = self.y_gnd_top_end
        else:
            self.y_gnd_top_start = self.y_gnd_top_end = None
            self.domain_height = self.y_top_end

        # Skin depth
        self.delta_s = np.sqrt(2 / (self.omega * mu_0 * self.sigma_cond))

        # --- Grid Generation ---
        self.x, self.dx_array = self._grid_x(nx)
        self.y, self.dy_array = self._grid_y(ny)
        self.X, self.Y = np.meshgrid(self.x, self.y)

        # --- Field Arrays ---
        self.V = np.zeros(self.X.shape)
        self.epsilon_r = np.ones(self.X.shape)
        self.signal_mask = np.zeros(self.X.shape, dtype=bool)
        self.signal_p_mask = np.zeros(self.X.shape, dtype=bool)  # Positive trace
        self.signal_n_mask = np.zeros(self.X.shape, dtype=bool)  # Negative trace
        self.ground_mask = np.zeros(self.X.shape, dtype=bool)

        self._setup_geometry()

    def _grid_x(self, n):
        """Generate non-uniform grid with refinement around both traces."""
        cx = self.domain_width / 2

        # Left trace (negative in odd mode)
        xl_left = cx - self.w - self.s/2
        xr_left = cx - self.s/2

        # Right trace (positive in odd mode)
        xl_right = cx + self.s/2
        xr_right = cx + self.w + self.s/2

        ds = self.delta_s
        corner = min(3*ds, self.w/4)
        ncorner = max(3, n // 20)

        if not self.use_sm:
            # 5-segment grid: left-side | left-trace | gap | right-trace | right-side
            n_seg = n // 5

            a = self._smooth_transition(0, xl_left, n_seg, curve_end="end")
            b = self._smooth_transition(xl_left, xr_left, n_seg, curve_end="both", beta=6)
            c = self._smooth_transition(xr_left, xl_right, n_seg, curve_end="both")
            d = self._smooth_transition(xl_right, xr_right, n_seg, curve_end="both", beta=6)
            e = self._smooth_transition(xr_right, self.domain_width, n - 4*n_seg + 4, curve_end="start")

            x = np.concatenate([a, b[1:], c[1:], d[1:], e[1:]])
        else:
            # With solder mask, extend grid
            xsl_left = xl_left - self.sm_t_side
            xsr_left = xr_left + self.sm_t_side
            xsl_right = xl_right - self.sm_t_side
            xsr_right = xr_right + self.sm_t_side

            n_seg = n // 9

            a = self._smooth_transition(0, xsl_left, n_seg, curve_end="end")
            b = self._smooth_transition(xsl_left, xl_left, n_seg//2, curve_end="both")
            c = self._smooth_transition(xl_left, xr_left, n_seg, curve_end="both")
            d = self._smooth_transition(xr_left, xsr_left, n_seg//2, curve_end="both")
            e = self._smooth_transition(xsr_left, xsl_right, n_seg, curve_end="both")
            f = self._smooth_transition(xsl_right, xl_right, n_seg//2, curve_end="both")
            g = self._smooth_transition(xl_right, xr_right, n_seg, curve_end="both")
            h = self._smooth_transition(xr_right, xsr_right, n_seg//2, curve_end="both")
            i = self._smooth_transition(xsr_right, self.domain_width, n - 8*n_seg + 8, curve_end="start")

            x = np.concatenate([a, b[1:], c[1:], d[1:], e[1:], f[1:], g[1:], h[1:], i[1:]])

        return x, np.diff(x)

    def _grid_y(self, n):
        """Generate non-uniform grid with refinement in vertical direction."""
        n_gnd = max(3, int(n * 0.05))

        ds = self.delta_s
        corner = min(3*ds, self.t/4)
        ncorner = self.skin_cells // 4
        n_trace = self.skin_cells + 4 - 2*ncorner

        if not self.use_sm:
            n_remain = n - (n_gnd * 2 if self.has_top_gnd else n_gnd) - n_trace
            n_sub_bot = int(n_remain * 0.4)
            n_sub_top = n_remain - n_sub_bot

            l1 = self._smooth_transition(self.y_gnd_bot_start, self.y_gnd_bot_end, n_gnd, curve_end="both")
            l2 = self._smooth_transition(self.y_sub_start, self.y_sub_end, n_sub_bot, curve_end="both")
            l3_cb = self._smooth_transition(self.y_trace_start, self.y_trace_start+corner, ncorner, curve_end="start")
            l3 = self._smooth_transition(self.y_trace_start+corner, self.y_trace_end-corner, n_trace, curve_end="both", beta=2)
            l3_ct = self._smooth_transition(self.y_trace_end-corner, self.y_trace_end, ncorner, curve_end="end")
            l4 = self._smooth_transition(self.y_top_start, self.y_top_end, n_sub_top, curve_end="start")
            grid_parts = [l1, l2[1:], l3_cb[1:], l3[1:], l3_ct[1:], l4[1:]]
        else:
            n_sm = 10
            n_remain = n - (n_gnd * 2 if self.has_top_gnd else n_gnd) - n_trace - n_sm
            n_sub_bot = int(n_remain * 0.4)
            n_sub_top = n_remain - n_sub_bot

            l1 = self._smooth_transition(self.y_gnd_bot_start, self.y_gnd_bot_end, n_gnd, curve_end="both")
            l2 = self._smooth_transition(self.y_sub_start, self.y_sub_end, n_sub_bot, curve_end="both")
            l3 = self._smooth_transition(self.y_trace_start, self.y_trace_end, n_trace, curve_end="both")
            l_sm = self._smooth_transition(self.y_trace_end, self.y_top_start, n_sm, curve_end="both", beta=8)
            l4 = self._smooth_transition(self.y_top_start, self.y_top_end, n_sub_top, curve_end="start")
            grid_parts = [l1, l2[1:], l3[1:], l_sm[1:], l4[1:]]

        if self.has_top_gnd:
            l5 = self._smooth_transition(self.y_gnd_top_start, self.y_gnd_top_end, n_gnd, curve_end="start")
            grid_parts.append(l5[1:])

        y = np.concatenate(grid_parts)
        return y, np.diff(y)

    def _setup_geometry(self):
        """Set up differential microstrip geometry with two traces."""
        cx = self.domain_width / 2

        # Left trace
        xl_left = cx - self.w - self.s/2
        xr_left = cx - self.s/2

        # Right trace
        xl_right = cx + self.s/2
        xr_right = cx + self.w + self.s/2

        # Solder mask extensions
        xsl_left = xl_left - self.sm_t_side
        xsr_left = xr_left + self.sm_t_side
        xsl_right = xl_right - self.sm_t_side
        xsr_right = xr_right + self.sm_t_side

        for i, yc in enumerate(self.y):
            # --- Permittivity ---
            if yc <= self.y_sub_end:
                self.epsilon_r[i, :] = self.er
            else:
                self.epsilon_r[i, :] = self.er_top

                if self.use_sm:
                    for j, xc in enumerate(self.x):
                        # Mask over substrate
                        if yc <= self.y_sm_sub_end:
                            self.epsilon_r[i, j] = self.sm_er

                        # Mask around left trace
                        if (xsl_left <= xc <= xsr_left) and (yc <= self.y_sm_trace_end):
                            self.epsilon_r[i, j] = self.sm_er

                        # Mask around right trace
                        if (xsl_right <= xc <= xsr_right) and (yc <= self.y_sm_trace_end):
                            self.epsilon_r[i, j] = self.sm_er

            # --- Conductor Masks ---
            if yc <= self.y_gnd_bot_end and self.boundaries[3] == "gnd":
                self.ground_mask[i, :] = True
            if self.has_top_gnd and yc >= self.y_gnd_top_start:
                self.ground_mask[i, :] = True

            # Traces (Signal)
            if self.y_trace_start <= yc <= self.y_trace_end:
                # Left trace (negative)
                idx_left = (self.x >= xl_left) & (self.x <= xr_left)
                self.signal_n_mask[i, idx_left] = True

                # Right trace (positive)
                idx_right = (self.x >= xl_right) & (self.x <= xr_right)
                self.signal_p_mask[i, idx_right] = True

        # Combined signal mask
        self.signal_mask = self.signal_p_mask | self.signal_n_mask
        self.conductor_mask = self.signal_mask | self.ground_mask

    def calculate_differential_parameters(self):
        """
        Calculate both odd-mode and even-mode impedances with losses.

        Returns:
        --------
        dict with keys:
            Z_odd : Odd-mode impedance (Ohms)
            Z_even : Even-mode impedance (Ohms)
            Z_diff : Differential impedance = 2*Z_odd (Ohms)
            Z_common : Common-mode impedance = Z_even/2 (Ohms)
            eps_eff_odd : Odd-mode effective permittivity
            eps_eff_even : Even-mode effective permittivity
            alpha_c_odd : Odd-mode conductor loss (dB/m)
            alpha_c_even : Even-mode conductor loss (dB/m)
            alpha_d_odd : Odd-mode dielectric loss (dB/m)
            alpha_d_even : Even-mode dielectric loss (dB/m)
        """
        self.V = np.zeros(self.X.shape)
        self.V[self.signal_p_mask] = 1.0   # Right trace: +1V
        self.V[self.signal_n_mask] = -1.0  # Left trace: -1V
        self.V[self.ground_mask] = 0.0

        self.solve_laplace()
        Ex_odd, Ey_odd = self.compute_fields()

        # For odd mode, calculate capacitance for one trace (they're symmetric)
        # Use only positive trace for charge calculation
        self.signal_mask = self.signal_p_mask.copy()
        C_odd = self.calculate_capacitance()

        self.solve_laplace(vacuum=True)
        C0_odd = self.calculate_capacitance(vacuum=True)

        eps_eff_odd = C_odd / C0_odd
        Z_odd = 1 / (c * np.sqrt(C_odd * C0_odd))

        # Calculate losses for odd mode
        # For losses, we need to consider both traces
        self.signal_mask = self.signal_p_mask | self.signal_n_mask

        alpha_c_odd = self._calculate_differential_conductor_loss(Ex_odd, Ey_odd, Z_odd, 'odd')
        alpha_d_odd = self._calculate_differential_dielectric_loss(Ex_odd, Ey_odd, Z_odd)

        # --- EVEN MODE: Same polarity (+1V and +1V) ---
        self.V = np.zeros(self.X.shape)
        self.V[self.signal_p_mask] = 1.0   # Right trace: +1V
        self.V[self.signal_n_mask] = 1.0   # Left trace: +1V
        self.V[self.ground_mask] = 0.0

        self.solve_laplace()

        # For even mode, use positive trace for charge calculation
        self.signal_mask = self.signal_p_mask.copy()
        C_even = self.calculate_capacitance()
        Ex_even, Ey_even = self.compute_fields()

        # Calculate C0 (air dielectric)
        self.solve_laplace(vacuum=True)
        C0_even = self.calculate_capacitance(vacuum=True)

        eps_eff_even = C_even / C0_even
        Z_even = 1 / (c * np.sqrt(C_even * C0_even))

        # Calculate losses for even mode
        self.signal_mask = self.signal_p_mask | self.signal_n_mask

        alpha_c_even = self._calculate_differential_conductor_loss(Ex_even, Ey_even, Z_even, 'even')
        alpha_d_even = self._calculate_differential_dielectric_loss(Ex_even, Ey_even, Z_even)

        eps_eff_1 = self.calculate_eps_eff_from_fields(Ex_even, Ey_even)
        eps_eff_2 = self.calculate_eps_eff_from_fields(Ex_odd, Ey_odd)

        # --- DIFFERENTIAL AND COMMON MODE ---
        Z_diff = 2 * Z_odd
        Z_common = Z_even / 2

        # Restore signal mask
        self.signal_mask = self.signal_p_mask | self.signal_n_mask

        return Ex_odd, Ey_odd, Ex_even, Ey_even, {
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

    def _calculate_differential_conductor_loss(self, Ex, Ey, Z0, mode):
        """
        Calculate conductor loss for differential pair.
        For differential pairs, surface current is distributed on both traces.
        """
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

        #print(f"Power loss (1W):    {Pc_1W:.6e} W")
        #print(f"Conductor loss:     {alpha_db_per_m:.4f} dB/m")

        return alpha_db_per_m

    def _calculate_differential_dielectric_loss(self, Ex, Ey, Z0):
        """
        Calculate dielectric loss for differential pair.
        For differential mode, integrate over entire domain but power is normalized per trace.
        """
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
        # For odd/even mode, this is the power in one trace
        P_flow = 1.0 / (2 * Z0)

        # Attenuation constant: α = Pd / (2*P_flow)
        # Convert to dB/m: α_dB = 8.686 * α
        return 8.686 * (0.5 * Pd / (2 * P_flow))

def solve_differential_microstrip():
    solver = DifferentialMicrostripSolver2D(
            substrate_height=1.6e-3,
            trace_width=3e-3,
            trace_spacing=1e-3,
            trace_thickness=35e-6,
            gnd_thickness=16e-6,
            epsilon_r=4.5,
            freq=1e9,
            nx=200,
            ny=100
        )

    Ex_odd, Ey_odd, Ex_even, Ey_even, results = solver.calculate_differential_parameters()

    rlgc = solver.modal_to_physical_rlgc(results)
    print({i:np.array(rlgc[i]) for i in rlgc.keys()})

    solver.plot_geometry(unit="mm")
    solver.plot(Ex_odd, Ey_odd)
    solver.plot(Ex_even, Ey_even)

def solve_differential_stripline():
    solver = DifferentialMicrostripSolver2D(
            substrate_height=0.2e-3,
            trace_width=0.15e-3,
            trace_spacing=0.1e-3,
            trace_thickness=35e-6,
            gnd_thickness=16e-6,
            epsilon_r=4.1,
            epsilon_r_top=4.1,
            air_top=0.2e-3,
            freq=1e9,
            nx=200,
            ny=200,
            boundaries=["open", "open", "gnd", "gnd"]
        )

    Ex_odd, Ey_odd, Ex_even, Ey_even, results = solver.calculate_differential_parameters()

    rlgc = solver.modal_to_physical_rlgc(results)
    print({i:np.array(rlgc[i]) for i in rlgc.keys()})

    solver.plot_geometry(unit="mm")
    solver.plot(Ex_odd, Ey_odd)
    solver.plot(Ex_even, Ey_even)

if __name__ == "__main__":
    #solve_differential_stripline()
    #solve_differential_microstrip()
    #solve_gcpw()
    #solve_microstrip_embed()
    solve_microstrip()
    #solve_stripline()
