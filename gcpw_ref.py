#!/usr/bin/env python
from field_solver_ref import *

class GroundedCPWSolver2D(FieldSolver2D):
    def __init__(self, substrate_height, trace_width, trace_thickness,
                 gap, top_gnd_width, via_gap, via_diameter,
                 gnd_thickness=35e-6, epsilon_r=4.5, tan_delta=0.02,
                 sigma_diel=0.0, sigma_cond=5.8e7, epsilon_r_top=1,
                 air_top=None, air_side=None, freq=1e9,
                 nx=300, ny=300, skin_cells=50, boundaries=None,
                 use_sm=False, sm_t_sub=20e-6, sm_t_trace=20e-6,
                 sm_t_side=20e-6, sm_er=3.5, sm_tand=0.02):

        self.h = substrate_height
        self.w = trace_width
        self.t = trace_thickness
        self.gap = gap # Gap to top ground
        self.top_gnd_w = top_gnd_width
        self.via_gap = via_gap # Edge-to-edge distance: signal edge to via edge
        self.via_d = via_diameter
        self.t_gnd = gnd_thickness
        self.er = epsilon_r
        self.er_top = epsilon_r_top
        self.tan_delta = tan_delta
        self.sigma_diel = sigma_diel
        self.sigma_cond = sigma_cond

        # Solder Mask
        self.use_sm = use_sm
        self.sm_t_sub = sm_t_sub
        self.sm_t_trace = sm_t_trace
        self.sm_t_side = sm_t_side
        self.sm_er = sm_er
        self.sm_tand = sm_tand

        self.f = freq
        self.omega = 2 * np.pi * freq
        self.nx, self.ny = nx, ny
        self.skin_cells = skin_cells

        # Total width: focus on the active area, but we will extend grounds to edges
        self.active_width = self.w + 2 * max(self.gap + self.top_gnd_w, self.via_gap + self.via_d)

        if air_side is None:
            self.domain_width = self.active_width * 1.5
        else:
            self.domain_width = self.active_width + 2 * air_side

        self.boundaries = boundaries if boundaries else ["open", "open", "open", "gnd"]

        # --- Y-Coordinates (Unchanged logic) ---
        self.y_gnd_bot_start = 0.0
        self.y_gnd_bot_end   = self.t_gnd
        self.y_sub_start     = self.y_gnd_bot_end
        self.y_sub_end       = self.y_sub_start + self.h
        self.y_trace_start   = self.y_sub_end
        self.y_trace_end     = self.y_trace_start + self.t
        self.y_sm_sub_end    = self.y_sub_end + self.sm_t_sub
        self.y_sm_trace_end  = self.y_trace_end + self.sm_t_trace
        self.y_top_start     = self.y_trace_end if not use_sm else max(self.y_sm_sub_end, self.y_sm_trace_end)
        self.top_dielectric_h = air_top if air_top is not None else self.h * 10
        self.y_top_end       = self.y_top_start + self.top_dielectric_h
        self.has_top_gnd     = (self.boundaries[2] == "gnd")
        self.domain_height   = (self.y_top_end + self.t_gnd) if self.has_top_gnd else self.y_top_end

        if self.has_top_gnd:
            self.y_gnd_top_start = self.y_top_end
            self.y_gnd_top_end   = self.y_gnd_top_start + self.t_gnd
        else:
            self.y_gnd_top_start = self.y_gnd_top_end = None

        # Skin depth
        self.delta_s = np.sqrt(2 / (self.omega * 4e-7 * np.pi * self.sigma_cond))

        # --- Via & Conductor X-Coordinates ---
        self._calculate_geometry_x()

        # --- Grid and Arrays ---
        self.x, self.dx_array = self._grid_x(nx)
        self.y, self.dy_array = self._grid_y(ny)
        self.X, self.Y = np.meshgrid(self.x, self.y)
        self.V = np.zeros(self.X.shape)
        self.epsilon_r = np.ones(self.X.shape)
        self.signal_mask = np.zeros(self.X.shape, dtype=bool)
        self.ground_mask = np.zeros(self.X.shape, dtype=bool)

        self._setup_geometry()

    def _calculate_geometry_x(self):
        cx = self.domain_width / 2
        # Signal
        self.x_tr_l = cx - self.w/2
        self.x_tr_r = cx + self.w/2

        # Top Grounds
        self.x_gap_l = self.x_tr_l - self.gap
        self.x_gap_r = self.x_tr_r + self.gap

        # Vias (Single via on each side)
        # Inner edge of via is (trace_edge - via_gap)
        self.via_x_left_inner = self.x_tr_l - self.via_gap
        self.via_x_left_center = self.via_x_left_inner - (self.via_d / 2)

        self.via_x_right_inner = self.x_tr_r + self.via_gap
        self.via_x_right_center = self.via_x_right_inner + (self.via_d / 2)

    def _grid_x(self, n):
        cx = self.domain_width / 2

        xl_trace, xr_trace = self.x_tr_l, self.x_tr_r
        xl_gap, xr_gap = self.x_gap_l, self.x_gap_r

        ds = self.delta_s
        corner_w = min(3 * ds, self.w / 4)
        corner_gap = min(3 * ds, self.gap / 4)

        sm = self.sm_t_side if self.use_sm else 0.0

        n_corner = max(5, n // 40)
        n_outer_gnd = max(int(n * 0.12), 8)
        n_gap_bulk = max(int(n * 0.14), 10)
        n_trace_center = max(int(n * 0.15), 8)

        # Solder mask sidewall resolution
        n_sm = max(5, n // 40) if self.use_sm else 0

        # --- LEFT SIDE ---

        # 1. Left outer ground
        a = self._smooth_transition(
            0, xl_gap - corner_gap,
            n_outer_gnd, curve_end="end"
        )

        # 2. Left ground corner
        a_corner = self._smooth_transition(
            xl_gap - corner_gap, xl_gap,
            n_corner, curve_end="end"
        )

        # 2.5 Left gap SM sidewall (ground edge)
        if self.use_sm and sm > 0:
            b_sm1 = self._smooth_transition(
                xl_gap, xl_gap + sm,
                n_sm, curve_end="both", beta=4.0
            )
            gap_start = xl_gap + sm
        else:
            b_sm1 = None
            gap_start = xl_gap

        # 3. Left gap bulk
        b = self._smooth_transition(
            gap_start, xl_trace - sm,
            n_gap_bulk, curve_end="both", beta=6.0
        )

        # 3.5 Left gap SM sidewall (trace edge)
        if self.use_sm and sm > 0:
            b_sm2 = self._smooth_transition(
                xl_trace - sm, xl_trace,
                n_sm, curve_end="both", beta=4.0
            )
        else:
            b_sm2 = None

        # 4. Trace left corner
        c_cl = self._smooth_transition(
            xl_trace, xl_trace + corner_w,
            n_corner, curve_end="start"
        )

        # 5. Trace center
        c = self._smooth_transition(
            xl_trace + corner_w,
            xr_trace - corner_w,
            n_trace_center,
            curve_end="both"
        )

        # 6. Trace right corner
        c_cr = self._smooth_transition(
            xr_trace - corner_w, xr_trace,
            n_corner, curve_end="end"
        )

        # --- RIGHT SIDE ---

        # 7. Right gap SM sidewall (trace edge)
        if self.use_sm and sm > 0:
            d_sm1 = self._smooth_transition(
                xr_trace, xr_trace + sm,
                n_sm, curve_end="both", beta=4.0
            )
            gap_r_start = xr_trace + sm
        else:
            d_sm1 = None
            gap_r_start = xr_trace

        # 8. Right gap bulk
        d = self._smooth_transition(
            gap_r_start, xr_gap - sm,
            n_gap_bulk, curve_end="both", beta=6.0
        )

        # 8.5 Right gap SM sidewall (ground edge)
        if self.use_sm and sm > 0:
            d_sm2 = self._smooth_transition(
                xr_gap - sm, xr_gap,
                n_sm, curve_end="both", beta=4.0
            )
        else:
            d_sm2 = None

        # 9. Right ground corner
        e_corner = self._smooth_transition(
            xr_gap, xr_gap + corner_gap,
            n_corner, curve_end="start"
        )

        # Remaining outer ground
        used = (
            len(a) +
            len(a_corner[1:]) +
            (len(b_sm1[1:]) if b_sm1 is not None else 0) +
            len(b[1:]) +
            (len(b_sm2[1:]) if b_sm2 is not None else 0) +
            len(c_cl[1:]) +
            len(c[1:]) +
            len(c_cr[1:]) +
            (len(d_sm1[1:]) if d_sm1 is not None else 0) +
            len(d[1:]) +
            (len(d_sm2[1:]) if d_sm2 is not None else 0) +
            len(e_corner[1:])
        )

        remaining_n = max(n - used, 6)

        e = self._smooth_transition(
            xr_gap + corner_gap,
            self.domain_width,
            remaining_n + 1,
            curve_end="start"
        )

        # --- Combine all segments ---
        parts = [
            a,
            a_corner[1:]
        ]

        if b_sm1 is not None:
            parts.append(b_sm1[1:])
        parts.append(b[1:])
        if b_sm2 is not None:
            parts.append(b_sm2[1:])

        parts += [
            c_cl[1:],
            c[1:],
            c_cr[1:]
        ]

        if d_sm1 is not None:
            parts.append(d_sm1[1:])
        parts.append(d[1:])
        if d_sm2 is not None:
            parts.append(d_sm2[1:])

        parts += [
            e_corner[1:],
            e[1:]
        ]

        x = np.concatenate(parts)

        return x, np.diff(x)

    def _grid_y(self, n):
        # Same as microstrip - reuse the vertical grid logic
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
            l3 = self._smooth_transition(self.y_trace_start+corner, self.y_trace_end-corner, n_trace, curve_end="both", beta=1)
            l3_ct = self._smooth_transition(self.y_trace_end-corner, self.y_trace_end, ncorner, curve_end="end")
            l4 = self._smooth_transition(self.y_top_start, self.y_top_end, n_sub_top, curve_end="start")
            grid_parts = [l1, l2[1:], l3_cb[1:], l3[1:], l3_ct[1:], l4[1:]]
        else:
            n_sm = 20
            n_remain = n - (n_gnd * 2 if self.has_top_gnd else n_gnd) - n_trace - n_sm
            n_sub_bot = int(n_remain * 0.4)
            n_sub_top = n_remain - n_sub_bot

            l1 = self._smooth_transition(self.y_gnd_bot_start, self.y_gnd_bot_end, n_gnd, curve_end="both")
            l2 = self._smooth_transition(self.y_sub_start, self.y_sub_end, n_sub_bot, curve_end="both")
            l3_cb = self._smooth_transition(self.y_trace_start, self.y_trace_start+corner, ncorner, curve_end="start")
            l3 = self._smooth_transition(self.y_trace_start+corner, self.y_trace_end-corner, n_trace, curve_end="both", beta=1)
            l3_ct = self._smooth_transition(self.y_trace_end-corner, self.y_trace_end, ncorner, curve_end="end")
            l_sm = self._smooth_transition(self.y_trace_end, self.y_top_start, n_sm, curve_end="both")
            l4 = self._smooth_transition(self.y_top_start, self.y_top_end, n_sub_top, curve_end="start")
            grid_parts = [l1, l2[1:], l3_cb[1:], l3[1:], l3_ct[1:], l_sm[1:], l4[1:]]

        if self.has_top_gnd:
            l5 = self._smooth_transition(self.y_gnd_top_start, self.y_gnd_top_end, n_gnd, curve_end="start")
            grid_parts.append(l5[1:])

        y = np.concatenate(grid_parts)
        return y, np.diff(y)

    def _setup_geometry(self):
        cx = self.domain_width / 2

        for i, yc in enumerate(self.y):
            for j, xc in enumerate(self.x):
                # --- 1. Permittivity (Dielectrics) ---
                if yc <= self.y_sub_end:
                    # Substrate Core
                    self.epsilon_r[i, j] = self.er
                else:
                    # Default to top dielectric (Air)
                    self.epsilon_r[i, j] = self.er_top

                # Apply Solder Mask if enabled
                if self.use_sm:
                    # Define metal regions for easier checking
                    is_in_signal = (self.x_tr_l <= xc <= self.x_tr_r) and \
                                  (self.y_trace_start <= yc <= self.y_trace_end)

                    is_in_top_gnd = ((xc <= self.x_gap_l) or (xc >= self.x_gap_r)) and \
                                   (self.y_trace_start <= yc <= self.y_trace_end)

                    # Check if point is in metal (signal or top ground)
                    is_in_metal = is_in_signal or is_in_top_gnd

                    # --- Zone A: Over substrate in gaps (between grounds and signal) ---
                    is_in_gap = (self.x_gap_l < xc < self.x_tr_l) or \
                               (self.x_tr_r < xc < self.x_gap_r)
                    if is_in_gap and (self.y_sub_end < yc <= self.y_sub_end + self.sm_t_sub):
                        self.epsilon_r[i, j] = self.sm_er

                    # --- Zone B: On top of metals (signal and top grounds) ---
                    # On top of signal trace
                    if (self.x_tr_l <= xc <= self.x_tr_r) and \
                       (self.y_trace_end < yc <= self.y_trace_end + self.sm_t_trace):
                        self.epsilon_r[i, j] = self.sm_er

                    # On top of left ground
                    if (xc <= self.x_gap_l) and \
                       (self.y_trace_end < yc <= self.y_trace_end + self.sm_t_trace):
                        self.epsilon_r[i, j] = self.sm_er

                    # On top of right ground
                    if (xc >= self.x_gap_r) and \
                       (self.y_trace_end < yc <= self.y_trace_end + self.sm_t_trace):
                        self.epsilon_r[i, j] = self.sm_er

                    # --- Zone C: Conformal coating at metal edges (sidewalls) ---
                    # Left edge of signal trace
                    if (self.x_tr_l - self.sm_t_side <= xc < self.x_tr_l) and \
                       (self.y_sub_end < yc <= self.y_trace_end + self.sm_t_trace):
                        self.epsilon_r[i, j] = self.sm_er

                    # Right edge of signal trace
                    if (self.x_tr_r < xc <= self.x_tr_r + self.sm_t_side) and \
                       (self.y_sub_end < yc <= self.y_trace_end + self.sm_t_trace):
                        self.epsilon_r[i, j] = self.sm_er

                    # Right edge of left ground (inner edge)
                    if (self.x_gap_l < xc <= self.x_gap_l + self.sm_t_side) and \
                       (self.y_sub_end < yc <= self.y_trace_end + self.sm_t_trace):
                        self.epsilon_r[i, j] = self.sm_er

                    # Left edge of right ground (inner edge)
                    if (self.x_gap_r - self.sm_t_side <= xc < self.x_gap_r) and \
                       (self.y_sub_end < yc <= self.y_trace_end + self.sm_t_trace):
                        self.epsilon_r[i, j] = self.sm_er

                # --- 2. Conductor Masks ---
                # Bottom Ground Plane
                if yc <= self.y_gnd_bot_end:
                    self.ground_mask[i, j] = True

                # Top Grounds (Extend from gap edges to domain edges)
                if self.y_trace_start <= yc <= self.y_trace_end:
                    if xc <= self.x_gap_l or xc >= self.x_gap_r:
                        self.ground_mask[i, j] = True

                # Vias (One per side, from via_gap edge to domain edges)
                if self.y_sub_start <= yc <= self.y_sub_end:
                    if xc <= self.via_x_left_inner or xc >= self.via_x_right_inner:
                        self.ground_mask[i, j] = True

                # Signal Trace
                if self.y_trace_start <= yc <= self.y_trace_end:
                    if self.x_tr_l <= xc <= self.x_tr_r:
                        self.signal_mask[i, j] = True

        # Finalize potentials
        self.V[self.signal_mask] = 1.0
        self.V[self.ground_mask] = 0.0
        self.conductor_mask = self.signal_mask | self.ground_mask

def solve_gcpw(plots=True):
    print("Solving GCPW...")
    solver = GroundedCPWSolver2D(
        substrate_height=1.6e-3,
        trace_width=0.3e-3,
        trace_thickness=35e-6,
        gap=0.15e-3,
        top_gnd_width=5e-3,
        via_gap=0.5e-3,
        via_diameter=0.3e-3,
        use_sm=False,
        sm_er=3.5,
        nx=200, ny=200,
    )

    Z0, eps_eff, C, C0 = solver.calculate_parameters()
    Ex, Ey = solver.compute_fields()

    alpha_cond, J = solver.calculate_conductor_loss(Ex, Ey, Z0)
    alpha_diel = solver.calculate_dielectric_loss(Ex, Ey, Z0)
    alpha_total = alpha_cond + alpha_diel

    z, rlgc, eps_eff_mode = solver.rlgc(alpha_cond, alpha_diel, C, Z0)
    print(f"Z (complex) = {z:.2f} ohm, eps_eff {eps_eff_mode:.3f}, RLGC {rlgc}")

    print(f"\n{'='*50}")
    print(f"GCPW ANALYSIS RESULTS")
    print(f"{'='*50}")
    print(f"Characteristic Impedance Z0:  {Z0:.2f} Ω")
    print(f"Effective Permittivity εᵣₑff: {eps_eff:.3f}")
    print(f"Losses (dB/m) @ 1GHz:         Diel={alpha_diel:.4f}, Cond={alpha_cond:.4g}")
    print(f"Total Attenuation:            {alpha_total:.4f} dB/m")
    print(f"{'='*50}\n")

    if plots:
        solver.plot_geometry(mesh_stride=1)
        solver.plot(Ex, Ey)

if __name__ == "__main__":
    solve_gcpw(False)
