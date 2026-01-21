#!/usr/bin/env python
import sys
from field_solver_ref import *

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


def solve_differential_microstrip(plots=False):
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

    print(results)
    rlgc = solver.modal_to_physical_rlgc(results)
    print({i:np.array(rlgc[i]) for i in rlgc.keys()})

    if plots:
        solver.plot_geometry(unit="mm")
        solver.plot(Ex_odd, Ey_odd)
        solver.plot(Ex_even, Ey_even)

def solve_differential_stripline(plots=False):
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

    print(results)
    rlgc = solver.modal_to_physical_rlgc(results)
    print({i:np.array(rlgc[i]) for i in rlgc.keys()})

    if plots:
        solver.plot_geometry(unit="mm")
        solver.plot(Ex_odd, Ey_odd)
        solver.plot(Ex_even, Ey_even)

if __name__ == "__main__":
    plots = False
    if len(sys.argv) > 1:
        if "plot" in sys.argv[1]:
            plots = True
    solve_differential_microstrip(plots)
