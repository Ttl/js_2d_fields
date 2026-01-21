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

    def _compute_refine_metrics(self, Ex, Ey):
        """
        For each grid interval, compute a metric indicating how much refinement
        would help, based on voltage gradients and field energy in adjacent cells.
        """
        ny, nx = self.V.shape

        # Metric for splitting interval [x[j], x[j+1]]
        x_metrics = np.zeros(len(self.x) - 1)
        # Metric for splitting interval [y[i], y[i+1]]
        y_metrics = np.zeros(len(self.y) - 1)

        for i in range(ny - 1):
            for j in range(nx - 1):
                # Skip cells fully inside conductors
                if (self.conductor_mask[i, j] and
                    self.conductor_mask[min(i+1, ny-1), j] and
                    self.conductor_mask[i, min(j+1, nx-1)]):
                    continue

                eps = self.epsilon_r[i, j]

                # Voltage differences across this cell
                dV_x = abs(self.V[i, j+1] - self.V[i, j]) if j < nx - 1 else 0
                dV_y = abs(self.V[i+1, j] - self.V[i, j]) if i < ny - 1 else 0

                # Field magnitude for weighting
                E2 = Ex[i, j]**2 + Ey[i, j]**2
                E_mag = np.sqrt(E2) if E2 > 0 else 1e-12

                # Boundary detection
                is_boundary = (not self.conductor_mask[i, j] and (
                    (i > 0 and self.conductor_mask[i-1, j]) or
                    (i < ny-1 and self.conductor_mask[i+1, j]) or
                    (j > 0 and self.conductor_mask[i, j-1]) or
                    (j < nx-1 and self.conductor_mask[i, j+1])))
                boundary_mult = 2.0 if is_boundary else 1.0

                # Weight by field strength, permittivity, and boundary importance
                weight = E_mag * eps * boundary_mult

                # Accumulate to the interval metrics
                if j < len(x_metrics):
                    x_metrics[j] += dV_x * weight
                if i < len(y_metrics):
                    y_metrics[i] += dV_y * weight

        return x_metrics, y_metrics


    def _select_lines_to_refine(self, x_metrics, y_metrics, frac=0.15):
        """
        Select which grid intervals to split, respecting left-right symmetry.
        """
        x_center = (self.x[0] + self.x[-1]) / 2
        x_symmetric = self._check_symmetry(self.x, x_center)

        if x_symmetric:
            x_metrics = self._symmetrize_metrics(x_metrics)

        # Decide how many x vs y lines based on relative total metric
        total_x = np.sum(x_metrics)
        total_y = np.sum(y_metrics)
        total = total_x + total_y

        if total < 1e-15:
            return set(), set()

        n_total = int(frac * (len(x_metrics) + len(y_metrics)))
        n_total = max(1, n_total)

        # Allocate proportionally to where the error is
        n_x = int(n_total * total_x / total)
        n_y = n_total - n_x

        # Select top intervals
        x_ranked = np.argsort(x_metrics)[::-1]
        y_ranked = np.argsort(y_metrics)[::-1]

        selected_x = set()
        selected_y = set()

        for j in x_ranked[:n_x]:
            if x_metrics[j] > 0:
                selected_x.add(j)
                if x_symmetric:
                    partner = len(x_metrics) - 1 - j
                    if 0 <= partner < len(x_metrics):
                        selected_x.add(partner)

        for i in y_ranked[:n_y]:
            if y_metrics[i] > 0:
                selected_y.add(i)

        return selected_x, selected_y


    def _check_symmetry(self, coords, center, tol=1e-10):
        """Check if coordinate array is symmetric about center."""
        n = len(coords)
        for k in range(n // 2):
            left = coords[k]
            right = coords[n - 1 - k]
            if abs((left - center) + (right - center)) > tol:
                return False
        return True


    def _symmetrize_metrics(self, metrics):
        """Average metrics for symmetric pairs."""
        n = len(metrics)
        result = metrics.copy()
        for k in range(n // 2):
            avg = 0.5 * (metrics[k] + metrics[n - 1 - k])
            result[k] = avg
            result[n - 1 - k] = avg
        return result


    def _refine_selected_lines(self, selected_x, selected_y):
        """Add new grid lines at midpoints of selected intervals."""
        x_center = (self.x[0] + self.x[-1]) / 2
        x_symmetric = self._check_symmetry(self.x, x_center)

        new_x = set()
        new_y = set()

        for j in selected_x:
            midpoint = 0.5 * (self.x[j] + self.x[j + 1])

            if x_symmetric:
                # Only generate midpoints on one side of center
                if midpoint <= x_center:
                    new_x.add(midpoint)
                    symmetric_point = 2 * x_center - midpoint
                    if self.x[0] < symmetric_point < self.x[-1]:
                        new_x.add(symmetric_point)
            else:
                new_x.add(midpoint)

        for i in selected_y:
            midpoint = 0.5 * (self.y[i] + self.y[i + 1])
            new_y.add(midpoint)

        self.x = np.array(sorted(set(self.x) | new_x))
        self.y = np.array(sorted(set(self.y) | new_y))


    def refine_mesh(self, Ex, Ey, frac=0.15):
        """Main refinement routine."""
        x_metrics, y_metrics = self._compute_refine_metrics(Ex, Ey)
        selected_x, selected_y = self._select_lines_to_refine(x_metrics, y_metrics, frac)
        self._refine_selected_lines(selected_x, selected_y)

    def _compute_energy_error(self, Ex, Ey, prev_energy):
        """
        Compute relative change in stored electromagnetic energy.
        Fields must be interpolated to cell centers to match epsilon grid.
        """
        energy = 0.0
        for i in range(self.ny - 1):
            for j in range(self.nx - 1):
                if self.conductor_mask[i, j]:
                    continue

                E2 = Ex[i, j]**2 + Ey[i, j]**2
                dA = self.dx_array[j] * self.dy_array[i]
                energy += 0.5 * epsilon_0 * self.epsilon_r[i, j] * E2 * dA

        if prev_energy is None:
            return energy, 1.0

        rel_error = abs(energy - prev_energy) / max(abs(prev_energy), 1e-12)
        return energy, rel_error

    def _compute_parameter_error(self, Z0, C, prev_Z0, prev_C):
        """Track convergence of the quantities you actually care about."""
        if prev_Z0 is None:
            return 1.0

        z_err = abs(Z0 - prev_Z0) / max(abs(prev_Z0), 1e-12)
        c_err = abs(C - prev_C) / max(abs(prev_C), 1e-12)
        return max(z_err, c_err)

    def solve_adaptive(self,
                       max_iters=10,
                       refine_frac=0.15,
                       energy_tol=0.05,
                       param_tol=0.1,# Error in Z0/C
                       max_nodes=20000,
                       min_converged_passes=2):
        """
        Adaptive mesh solve with robust convergence criteria.

        Stores computed fields in self.Ex and self.Ey.
        Returns: Z0, eps_eff, C, C0
        """
        prev_energy = None
        prev_Z0 = None
        prev_C = None
        converged_count = 0

        for it in range(max_iters):
            Z0, eps_eff, C, C0 = self.calculate_parameters()
            Ex, Ey = self.compute_fields()

            # Energy-based error (primary criterion)
            energy, energy_err = self._compute_energy_error(Ex, Ey, prev_energy)

            # Parameter-based error (what you actually care about)
            param_err = self._compute_parameter_error(Z0, C, prev_Z0, prev_C)

            print(f"Pass {it+1}: Energy err={energy_err:.3g}, Param err={param_err:.3g}, "
                  f"Nodes={self.nx}x{self.ny}")

            # Check convergence (both criteria must be met)
            if prev_energy is not None:
                if energy_err < energy_tol and param_err < param_tol:
                    converged_count += 1
                    if converged_count >= min_converged_passes:
                        print(f"Converged after {it+1} passes")
                        break
                else:
                    converged_count = 0  # Reset if we diverge

            prev_energy = energy
            prev_Z0 = Z0
            prev_C = C

            # Node budget check
            if self.nx * self.ny > max_nodes:
                print("Node budget reached")
                break

            # Refine mesh
            if it != max_iters - 1:
                self.refine_mesh(Ex, Ey, frac=refine_frac)
                self._setup_geometry()

        # Store fields as class members
        self.Ex = Ex
        self.Ey = Ey

        return Z0, eps_eff, C, C0

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

        mesh_alpha = 0.5
        mesh_color = "white"
        mesh_linewidth = 0.2
        # Plot vertical grid lines
        scale = 1e3
        for i in range(0, len(self.x)):
            axs[0].plot([self.x[i] * scale, self.x[i] * scale],
                    [self.y[0] * scale, self.y[-1] * scale],
                    color=mesh_color, alpha=mesh_alpha, linewidth=mesh_linewidth,
                    zorder=10)

        # Plot horizontal grid lines
        for j in range(0, len(self.y)):
            axs[0].plot([self.x[0] * scale, self.x[-1] * scale],
                    [self.y[j] * scale, self.y[j] * scale],
                    color=mesh_color, alpha=mesh_alpha, linewidth=mesh_linewidth,
                    zorder=10)

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

