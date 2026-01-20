#!/usr/bin/env python
import sys
from field_solver_ref import *

class MicrostripSolver2D(FieldSolver2D):
    def __init__(self, substrate_height, trace_width, trace_thickness,
                 gnd_thickness=35e-6, epsilon_r=4.5, tan_delta=0.02,
                 sigma_diel=0.0, sigma_cond=5.8e7, epsilon_r_top=1,
                 air_top=None, air_side=None, freq=1e9,
                 nx=300, ny=300, skin_cells=100, boundaries=None,
                 use_sm=False, sm_t_sub=20e-6, sm_t_trace=20e-6,
                 sm_t_side=20e-6, sm_er=3.5, sm_tand=0.02,
                 gnd_cut_width=0.0,
                 gnd_cut_sub_h=0.0,
                 top_diel_h=0.0,
                 top_diel_er=1.0,
                 top_diel_tand=0.0):

        self.h = substrate_height
        self.w = trace_width
        self.t = trace_thickness
        self.t_gnd = gnd_thickness
        self.er = epsilon_r
        self.er_top = epsilon_r_top
        self.tan_delta = tan_delta
        self.sigma_diel = sigma_diel
        self.sigma_cond = sigma_cond

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
        self.skin_cells = skin_cells

        if air_side is None:
            self.domain_width = 2 * max(self.w * 8, self.h * 15)
        else:
            self.domain_width = self.w + 2 * air_side

        self.boundaries = boundaries if boundaries else ["open", "open", "open", "gnd"]

        # --- Physical Y-Coordinates ---
        # Bottom extension for cut ground
        self.y_ext_start = self.t_gnd
        self.y_ext_end = self.t_gnd + self.gnd_cut_sub_h

        # New bottom ground plane location
        self.y_gnd_bot_start = self.y_ext_end
        self.y_gnd_bot_end = self.y_gnd_bot_start + self.t_gnd
        if gnd_cut_width == 0:
            self.y_gnd_bot_end = self.y_gnd_bot_start

        self.y_sub_start = self.y_gnd_bot_end
        self.y_sub_end = self.y_sub_start + self.h

        # --- TOP DIELECTRIC ---
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

        if air_top is None:
            self.top_dielectric_h = self.h * 15
            self.has_top_gnd = False
        else:
            self.top_dielectric_h = air_top + self.t
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
        self.ground_mask = np.zeros(self.X.shape, dtype=bool)

        self._setup_geometry()

    # --------------------------------------------------
    # GEOMETRY SETUP
    # --------------------------------------------------
    def _setup_geometry(self):
        cx = self.domain_width / 2
        xl, xr = cx - self.w / 2, cx + self.w / 2
        xsl, xsr = xl - self.sm_t_side, xr + self.sm_t_side

        cut_l = cx - self.gnd_cut_width / 2
        cut_r = cx + self.gnd_cut_width / 2
        tol = 1e-11

        for i, yc in enumerate(self.y):
            # --- PERMITTIVITY ---
            if yc - tol <= self.y_sub_end + tol:
                self.epsilon_r[i, :] = self.er

            elif yc - tol <= self.y_top_diel_end + tol:
                self.epsilon_r[i, :] = self.top_diel_er

            else:
                self.epsilon_r[i, :] = self.er_top

                if self.use_sm:
                    for j, xc in enumerate(self.x):
                        if yc <= self.y_sm_sub_end - tol:
                            self.epsilon_r[i, j] = self.sm_er
                        if (xsl - tol <= xc <= xsr + tol) and (yc - tol <= self.y_sm_trace_end + tol):
                            self.epsilon_r[i, j] = self.sm_er

            # --- BOTTOM GROUND ---
            if yc >= self.y_gnd_bot_start - tol and yc <= self.y_gnd_bot_end + tol:
                for j, xc in enumerate(self.x):
                    if self.gnd_cut_width > 0 and cut_l <= xc <= cut_r:
                        continue
                    self.ground_mask[i, j] = True

            if yc <= self.y_ext_start + tol:
                self.ground_mask[i, :] = True

            # --- TOP GROUND ---
            if self.has_top_gnd and yc >= self.y_gnd_top_start - tol:
                self.ground_mask[i, :] = True

            # --- TRACE ---
            if self.y_trace_start - tol <= yc <= self.y_trace_end + tol:
                idx_x = (self.x >= xl) & (self.x <= xr)
                self.signal_mask[i, idx_x] = True

        # Finalize potentials
        self.V[self.signal_mask] = 1.0
        self.V[self.ground_mask] = 0.0
        self.conductor_mask = self.signal_mask | self.ground_mask

    def _grid_x(self, n):
        """Generate x-grid with adaptive grading focused on field-critical regions."""
        cx = self.domain_width / 2
        xl, xr = cx - self.w / 2, cx + self.w / 2

        ds = self.delta_s
        corner = min(3 * ds, self.w / 4)
        ncorner = max(3, self.skin_cells // 4)

        # Collect all critical x-interfaces
        x_if = [0, xl, xr, self.domain_width]

        if self.gnd_cut_width > 0:
            cut_l = cx - self.gnd_cut_width / 2
            cut_r = cx + self.gnd_cut_width / 2
            x_if += [cut_l, cut_r]

        if self.use_sm:
            xsl = xl - self.sm_t_side
            xsr = xr + self.sm_t_side
            x_if += [xsl, xsr]

        x_if = np.array(sorted(set(x_if)))
        n_regions = len(x_if) - 1

        # Adaptive point allocation based on field importance
        region_weights = []
        for k in range(n_regions):
            x0, x1 = x_if[k], x_if[k + 1]
            width = x1 - x0

            # High weight for trace region and near-trace
            if x0 >= xl - 1e-15 and x1 <= xr + 1e-15:
                weight = 6.0  # Trace itself - highest field
            else:
                weight = 1.0

            # Ground cut gets extra points (field discontinuity)
            if self.gnd_cut_width > 0:
                cut_l = cx - self.gnd_cut_width / 2
                cut_r = cx + self.gnd_cut_width / 2
                if abs(x0 - cut_l) < 1e-15 or abs(x1 - cut_r) < 1e-15:
                    weight *= 1.5

            region_weights.append(weight * width)

        # Normalize and allocate points
        total_weight = sum(region_weights)
        region_points = []
        allocated = 0

        for k in range(n_regions):
            if k == n_regions - 1:
                pts = n - allocated
            else:
                pts = max(3, int(n * region_weights[k] / total_weight))
            region_points.append(pts)
            allocated += pts

        # Generate grid segments
        x_parts = []

        for k in range(n_regions):
            x0, x1 = x_if[k], x_if[k + 1]
            npts = region_points[k]

            # Determine grading strategy
            if x0 >= xl - 1e-15 and x1 <= xr + 1e-15:
                # TRACE REGION - fine mesh with corner grading
                nc = min(ncorner, npts // 3)
                n_mid = max(3, npts - 2 * nc)

                corner_left = min(corner, (x1 - x0) / 3)
                corner_right = min(corner, (x1 - x0) / 3)

                if nc > 0 and corner_left > 1e-15:
                    cl = self._smooth_transition(x0, x0 + corner_left, nc,
                                                curve_end="start", beta=4.0)
                else:
                    cl = np.array([x0])

                mid0 = x0 + corner_left
                mid1 = max(x1 - corner_right, mid0)

                mid = self._smooth_transition(mid0, mid1, n_mid,
                                             curve_end="both", beta=4.0)

                if nc > 0 and corner_right > 1e-15:
                    cr = self._smooth_transition(mid1, x1, nc,
                                                curve_end="end", beta=4.0)
                else:
                    cr = np.array([x1])

                seg = np.concatenate([cl, mid[1:], cr[1:]])

            else:
                # DIELECTRIC/AIR REGION - grade toward interfaces
                end_curve = "both"
                beta_val = 3.0

                # Grade AWAY from domain edges (fewer points at edges)
                if abs(x0) < 1e-15:
                    end_curve = "end"  # Dense at right (toward trace)
                    beta_val = 3.0
                elif abs(x1 - self.domain_width) < 1e-15:
                    end_curve = "start"  # Dense at left (toward trace)
                    beta_val = 3.0

                # Grade toward trace edges (more points near trace)
                if x1 <= xl + 1e-15 or x0 >= xr - 1e-15:
                    if x1 <= xl + 1e-15:
                        end_curve = "end"  # Dense toward trace (right side)
                    else:
                        end_curve = "start"  # Dense toward trace (left side)
                    beta_val = 4.0

                seg = self._smooth_transition(x0, x1, npts,
                                             curve_end=end_curve, beta=beta_val)

            if k > 0:
                seg = seg[1:]
            x_parts.append(seg)

        x = np.concatenate(x_parts)

        # Enforce exact interface locations
        for xi in x_if:
            if np.min(np.abs(x - xi)) > 1e-12:
                x = np.sort(np.append(x, xi))

        return x, np.diff(x)

    def _grid_y(self, n):
        """Generate y-grid with adaptive grading focused on field-critical regions."""

        ds = self.delta_s
        corner = min(3 * ds, self.t / 4)
        ncorner = max(3, self.skin_cells // 4)

        # Collect all critical y-interfaces
        y_if = [
            0.0,
            self.y_gnd_bot_start,
            self.y_gnd_bot_end,
            self.y_sub_start,
            self.y_sub_end,
        ]

        # Add ground cut interfaces
        if self.gnd_cut_width > 0:
            y_if += [self.y_ext_start, self.y_ext_end]

        if self.top_diel_h > 0:
            y_if += [self.y_top_diel_start, self.y_top_diel_end]

        y_if += [
            self.y_trace_start,
            self.y_trace_end,
            self.y_top_start,
            self.y_top_end
        ]

        if self.use_sm:
            y_if += [self.y_sm_sub_end, self.y_sm_trace_end]

        if self.has_top_gnd:
            y_if += [self.y_gnd_top_start, self.y_gnd_top_end]

        y_if = np.array(sorted(set(y_if)))
        n_regions = len(y_if) - 1

        # Adaptive point allocation based on field importance
        region_weights = []
        for k in range(n_regions):
            y0, y1 = y_if[k], y_if[k + 1]
            height = y1 - y0

            # TRACE - highest field concentration
            if (y0 >= self.y_trace_start - 1e-15 and
                y1 <= self.y_trace_end + 1e-15):
                weight = 10

            # GROUND PLANES - high current density
            elif (y1 <= self.y_ext_start + 1e-15 or
                  (self.has_top_gnd and
                   y0 >= self.y_gnd_top_start - 1e-15 and
                   y1 <= self.y_gnd_top_end + 1e-15)):
                weight = 0

            # SUBSTRATE BELOW CUTOUT - field region when ground cut enabled
            elif (self.gnd_cut_width > 0 and
                  y0 >= self.y_ext_start - 1e-15 and
                  y1 <= self.y_ext_end + 1e-15):
                weight = 0.75

            # SUBSTRATE - high field region
            elif (y0 >= self.y_sub_start - 1e-15 and
                  y1 <= self.y_sub_end + 1e-15):
                weight = 0.75

            # TOP DIELECTRIC - field transitions
            elif (self.top_diel_h > 0 and
                  y0 >= self.y_top_diel_start - 1e-15 and
                  y1 <= self.y_top_diel_end + 1e-15):
                weight = 1

            # SOLDER MASK - moderate field
            elif self.use_sm and (
                (y0 >= self.y_top_diel_end - 1e-15 and y1 <= self.y_sm_sub_end + 1e-15) or
                (y0 >= self.y_trace_end - 1e-15 and y1 <= self.y_sm_trace_end + 1e-15)):
                weight = 1

            # AIR/FAR-FIELD - lower field
            else:
                if self.has_top_gnd:
                    weight = 0.75
                else:
                    weight = 0.1

            region_weights.append(weight * height)

        # Normalize and allocate points
        total_weight = sum(region_weights)
        region_points = []
        allocated = 0

        for k in range(n_regions):
            if k == n_regions - 1:
                pts = n - allocated
            else:
                pts = max(3, int(n * region_weights[k] / total_weight))
            region_points.append(pts)
            allocated += pts

        # Generate grid segments
        y_parts = []

        for k in range(n_regions):
            y0, y1 = y_if[k], y_if[k + 1]
            npts = region_points[k]

            # Determine grading strategy
            end_curve = "both"
            beta_val = 3.0

            # Domain boundaries - grade AWAY from edges
            if abs(y0) < 1e-15:
                end_curve = "end"  # Dense at top (away from bottom edge)
                beta_val = 4
            if self.has_top_gnd and abs(y1 - self.y_gnd_top_end) < 1e-15:
                end_curve = "start"  # Dense at bottom (away from top edge)
                beta_val = 4
            elif not self.has_top_gnd and abs(y1 - self.y_top_end) < 1e-15:
                end_curve = "start"  # Dense at bottom (away from top edge)
                beta_val = 4

            # TRACE REGION - special corner grading
            if (y0 >= self.y_trace_start - 1e-15 and
                y1 <= self.y_trace_end + 1e-15):

                nc = min(ncorner, npts // 3)
                n_mid = max(4, npts - 2 * nc)

                corner_bot = min(corner, (y1 - y0) / 3)
                corner_top = min(corner, (y1 - y0) / 3)

                if nc > 0 and corner_bot > 1e-15:
                    cb = self._smooth_transition(y0, y0 + corner_bot, nc,
                                                curve_end="start", beta=4.0)
                else:
                    cb = np.array([y0])

                mid0 = y0 + corner_bot
                mid1 = max(y1 - corner_top, mid0)

                mid = self._smooth_transition(mid0, mid1, n_mid,
                                             curve_end="both", beta=2.0)

                if nc > 0 and corner_top > 1e-15:
                    ct = self._smooth_transition(mid1, y1, nc,
                                                curve_end="end", beta=4.0)
                else:
                    ct = np.array([y1])

                seg = np.concatenate([cb, mid[1:], ct[1:]])

            # GROUND PLANES - skin depth grading
            elif ((y0 >= self.y_gnd_bot_start - 1e-15 and
                   y1 <= self.y_gnd_bot_end + 1e-15) or
                  (self.has_top_gnd and
                   y0 >= self.y_gnd_top_start - 1e-15 and
                   y1 <= self.y_gnd_top_end + 1e-15)):

                seg = self._smooth_transition(y0, y1, npts,
                                             curve_end=end_curve, beta=4.0)

            # SUBSTRATE & TOP DIELECTRIC - grade toward trace interface
            elif ((y0 >= self.y_sub_start - 1e-15 and
                   y1 <= self.y_sub_end + 1e-15) or
                  (self.top_diel_h > 0 and
                   y0 >= self.y_top_diel_start - 1e-15 and
                   y1 <= self.y_top_diel_end + 1e-15)):

                # Grade toward trace (more points near trace)
                if y0 >= self.y_top_diel_start - 1e-15 and y1 <= self.y_top_diel_end + 1e-15:
                    # Top dielectric - grade toward trace at bottom
                    seg = self._smooth_transition(y0, y1, npts,
                                                 curve_end="start", beta=3.5)
                else:
                    # Substrate - grade toward trace at top
                    seg = self._smooth_transition(y0, y1, npts,
                                                 curve_end="end", beta=3.5)

            # AIR/FAR-FIELD - less aggressive grading
            else:
                seg = self._smooth_transition(y0, y1, npts,
                                             curve_end=end_curve, beta=beta_val)

            if k > 0:
                seg = seg[1:]
            y_parts.append(seg)

        y = np.concatenate(y_parts)

        # Enforce exact interface locations
        for yi in y_if:
            if np.min(np.abs(y - yi)) > 1e-12:
                y = np.sort(np.append(y, yi))

        return y, np.diff(y)


def solve_microstrip(plots=True):
    print("Solving Microstrip...")
    solver = MicrostripSolver2D(
        substrate_height=1.6e-3,
        trace_width=3e-3,
        trace_thickness=35e-6,
        gnd_thickness=35e-6,
        epsilon_r=4.5,
        tan_delta=0.02,
        sigma_cond=5.8e7,
        freq=1e9,
        nx=100, ny=100,
        skin_cells=50,
        use_sm=False,
        boundaries=["open", "open", "open", "gnd"]
    )

    Z0, eps_eff, C, C0 = solver.calculate_parameters()
    Ex, Ey = solver.compute_fields()

    alpha_cond, J = solver.calculate_conductor_loss(Ex, Ey, Z0)
    alpha_diel = solver.calculate_dielectric_loss(Ex, Ey, Z0)
    alpha_total = alpha_cond + alpha_diel

    z, rlgc, eps_eff_mode = solver.rlgc(alpha_cond, alpha_diel, C, Z0)
    print(f"Z (complex) = {z:.2f} ohm, eps_eff {eps_eff_mode:.3f}, RLGC {rlgc}")

    print(f"\n{'='*50}")
    print(f"MICROSTRIP ANALYSIS RESULTS")
    print(f"{'='*50}")
    print(f"Characteristic Impedance Z0:  {Z0:.2f} Ω")
    print(f"Effective Permittivity εᵣₑff: {eps_eff:.3f}")
    print(f"Losses (dB/m) @ 1GHz:         Diel={alpha_diel:.4f}, Cond={alpha_cond:.4g}")
    print(f"Total Attenuation:            {alpha_total:.4f} dB/m")
    print(f"{'='*50}\n")

    if plots:
        solver.plot_geometry()
        solver.plot(Ex, Ey)

def solve_microstrip_cut(plots=True):
    print("Solving Microstrip with gnd cut...")
    solver = MicrostripSolver2D(
        substrate_height=1.6e-3,
        trace_width=3e-3,
        trace_thickness=35e-6,
        gnd_thickness=35e-6,
        epsilon_r=4.5,
        tan_delta=0.02,
        sigma_cond=5.8e7,
        freq=1e9,
        nx=300, ny=300,
        skin_cells=50,
        use_sm=False,
        gnd_cut_width=3e-3,
        gnd_cut_sub_h=1e-3,
        boundaries=["open", "open", "open", "gnd"]
    )

    Z0, eps_eff, C, C0 = solver.calculate_parameters()
    Ex, Ey = solver.compute_fields()

    alpha_cond, J = solver.calculate_conductor_loss(Ex, Ey, Z0)
    alpha_diel = solver.calculate_dielectric_loss(Ex, Ey, Z0)
    alpha_total = alpha_cond + alpha_diel

    z, rlgc, eps_eff_mode = solver.rlgc(alpha_cond, alpha_diel, C, Z0)
    print(f"Z (complex) = {z:.2f} ohm, eps_eff {eps_eff_mode:.3f}, RLGC {rlgc}")

    print(f"\n{'='*50}")
    print(f"MICROSTRIP ANALYSIS RESULTS")
    print(f"{'='*50}")
    print(f"Characteristic Impedance Z0:  {Z0:.2f} Ω")
    print(f"Effective Permittivity εᵣₑff: {eps_eff:.3f}")
    print(f"Losses (dB/m) @ 1GHz:         Diel={alpha_diel:.4f}, Cond={alpha_cond:.4g}")
    print(f"Total Attenuation:            {alpha_total:.4f} dB/m")
    print(f"{'='*50}\n")

    if plots:
        solver.plot_geometry()
        solver.plot(Ex, Ey)

def solve_microstrip_embed(plots=True):
    print("Solving Embedded Microstrip...")
    solver = MicrostripSolver2D(
        substrate_height=1.6e-3,
        trace_width=3e-3,
        trace_thickness=35e-6,
        gnd_thickness=35e-6,
        epsilon_r=4.5,
        tan_delta=0.02,
        sigma_cond=5.8e7,
        freq=1e9,
        nx=200, ny=200,
        skin_cells=50,
        use_sm=False,
        #gnd_cut_width=3e-3,
        #gnd_cut_sub_h=1e-3,
        top_diel_h=0.2e-3,
        top_diel_er=4.5,
        boundaries=["open", "open", "open", "gnd"]
    )

    Z0, eps_eff, C, C0 = solver.calculate_parameters()
    Ex, Ey = solver.compute_fields()

    alpha_cond, J = solver.calculate_conductor_loss(Ex, Ey, Z0)
    alpha_diel = solver.calculate_dielectric_loss(Ex, Ey, Z0)
    alpha_total = alpha_cond + alpha_diel

    z, rlgc, eps_eff_mode = solver.rlgc(alpha_cond, alpha_diel, C, Z0)
    print(f"Z (complex) = {z:.2f} ohm, eps_eff {eps_eff_mode:.3f}, RLGC {rlgc}")

    print(f"\n{'='*50}")
    print(f"MICROSTRIP ANALYSIS RESULTS")
    print(f"{'='*50}")
    print(f"Characteristic Impedance Z0:  {Z0:.2f} Ω")
    print(f"Effective Permittivity εᵣₑff: {eps_eff:.3f}")
    print(f"Losses (dB/m) @ 1GHz:         Diel={alpha_diel:.4f}, Cond={alpha_cond:.4g}")
    print(f"Total Attenuation:            {alpha_total:.4f} dB/m")
    print(f"{'='*50}\n")

    if plots:
        solver.plot_geometry()
        solver.plot(Ex, Ey)

def solve_stripline(plots=True):
    print("Solving Stripline...")
    # Stripline: 0.2mm bottom dielectric, 0.2mm top dielectric
    solver = MicrostripSolver2D(
        substrate_height=0.2e-3,
        trace_width=0.15e-3,
        trace_thickness=35e-6,
        gnd_thickness=16e-6,
        epsilon_r=4.1,
        epsilon_r_top=4.1,
        air_top=0.2e-3,
        tan_delta=0.02,
        sigma_cond=5.8e7,
        freq=1e9,
        nx=100, ny=100,
        skin_cells=50,
        boundaries=["open", "open", "gnd", "gnd"]
    )

    Z0, eps_eff, C, C0 = solver.calculate_parameters()
    Ex, Ey = solver.compute_fields()

    alpha_cond, J = solver.calculate_conductor_loss(Ex, Ey, Z0)
    alpha_diel = solver.calculate_dielectric_loss(Ex, Ey, Z0)
    alpha_total = alpha_cond + alpha_diel

    z, rlgc, eps_eff_mode = solver.rlgc(alpha_cond, alpha_diel, C, Z0)
    print(f"Z (complex) = {z:.2f} ohm, eps_eff {eps_eff_mode:.3f}, RLGC {rlgc}")

    print(f"\n{'='*50}")
    print(f"STRIPLINE ANALYSIS RESULTS")
    print(f"{'='*50}")
    print(f"Characteristic Impedance Z0:  {Z0:.2f} Ω")
    print(f"Effective Permittivity εᵣₑff: {eps_eff:.3f}")
    print(f"Losses (dB/m) @ 1GHz:         Diel={alpha_diel:.4f}, Cond={alpha_cond:.4f}")
    print(f"Total Attenuation:            {alpha_total:.4f} dB/m")
    print(f"{'='*50}\n")

    if plots:
        solver.plot_geometry()
        solver.plot(Ex, Ey)


if __name__ == "__main__":
    plots = False
    if len(sys.argv) > 1:
        if "plot" in sys.argv[1]:
            plots = True
    solve_microstrip(plots)
