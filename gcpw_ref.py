#!/usr/bin/env python
import sys
from field_solver_ref import *

class GroundedCPWSolver2D(FieldSolver2D):
    def __init__(self, substrate_height, trace_width, trace_thickness,
                 gap, top_gnd_width, via_gap,
                 gnd_thickness=35e-6, epsilon_r=4.5, tan_delta=0.02,
                 sigma_diel=0.0, sigma_cond=5.8e7, epsilon_r_top=1,
                 air_top=None, air_side=None, freq=1e9,
                 nx=300, ny=300, skin_cells=50, boundaries=None,
                 use_sm=False, sm_t_sub=20e-6, sm_t_trace=20e-6,
                 sm_t_side=20e-6, sm_er=3.5, sm_tand=0.02):

        # Geometry
        # domain edge (via inside) - via edge - top gnd edge - gap - signal
        # Top gnd is from domain edge to gap
        # Via is under the top gnd from domain edge to via edge
        # via_gap is distance from top gnd gap to via edge
        self.h = substrate_height
        self.w = trace_width
        self.t = trace_thickness
        self.gap = gap # Gap to top ground
        self.top_gnd_w = top_gnd_width
        self.via_gap = via_gap # Edge-to-edge distance: signal edge to via edge
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
        self.active_width = self.w + 2 * max(self.gap + self.top_gnd_w, self.via_gap + self.gap)

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
        self.via_x_left_inner = self.x_tr_l - self.gap - self.via_gap

        self.via_x_right_inner = self.x_tr_r + self.gap + self.via_gap

    def _grid_x(self, n):
        """Generate x-grid with adaptive grading focused on field-critical regions."""
        cx = self.domain_width / 2
        
        # Define all critical x-interfaces
        x_if = [
            0,                          # Left domain edge
            self.via_x_left_inner,      # Left via inner edge
            self.x_gap_l,               # Left ground inner edge
            self.x_tr_l,                # Left trace edge
            self.x_tr_r,                # Right trace edge
            self.x_gap_r,               # Right ground inner edge
            self.via_x_right_inner,     # Right via inner edge
            self.domain_width           # Right domain edge
        ]
        
        # Add solder mask interfaces if enabled
        if self.use_sm:
            sm = self.sm_t_side
            x_if += [
                self.x_gap_l + sm,      # Left gap SM (ground side)
                self.x_tr_l - sm,       # Left gap SM (trace side)
                self.x_tr_r + sm,       # Right gap SM (trace side)
                self.x_gap_r - sm       # Right gap SM (ground side)
            ]
        
        x_if = np.array(sorted(set(x_if)))
        n_regions = len(x_if) - 1
        
        # Adaptive point allocation based on field importance
        # Use density per unit width rather than total weight
        region_densities = []
        
        for k in range(n_regions):
            x0, x1 = x_if[k], x_if[k + 1]
            width = x1 - x0
            
            # Assign density (points per unit width) based on field importance
            
            # Highest field: trace region
            if x0 >= self.x_tr_l - 1e-15 and x1 <= self.x_tr_r + 1e-15:
                density = 6.0  # Trace itself - many points
            
            # High field: gaps between trace and ground
            elif ((self.x_gap_l - 1e-15 <= x0 and x1 <= self.x_tr_l + 1e-15) or
                  (self.x_tr_r - 1e-15 <= x0 and x1 <= self.x_gap_r + 1e-15)):
                density = 6.0  # Gap regions - high field gradients
            
            # Medium field: ground plane regions (between via and gap)
            elif ((self.via_x_left_inner - 1e-15 <= x0 and x1 <= self.x_gap_l + 1e-15) or
                  (self.x_gap_r - 1e-15 <= x0 and x1 <= self.via_x_right_inner + 1e-15)):
                density = 2.0  # Top ground regions
            
            # Lower field: domain edge to via (low field, uniform)
            else:
                density = 0.5  # Outer regions - minimal points needed
            
            # Extra density for solder mask sidewall regions (field discontinuities)
            if self.use_sm:
                sm = self.sm_t_side
                # Check if this region is a SM sidewall
                is_sm_sidewall = (
                    abs(x0 - (self.x_gap_l)) < 1e-15 or
                    abs(x1 - (self.x_gap_l + sm)) < 1e-15 or
                    abs(x1 - (self.x_tr_l - sm)) < 1e-15 or
                    abs(x0 - (self.x_tr_l)) < 1e-15 or
                    abs(x0 - (self.x_tr_r)) < 1e-15 or
                    abs(x1 - (self.x_tr_r + sm)) < 1e-15 or
                    abs(x1 - (self.x_gap_r - sm)) < 1e-15 or
                    abs(x0 - (self.x_gap_r)) < 1e-15
                )
                if is_sm_sidewall:
                    density *= 1.5
            
            print(x0, x1, density)
            region_densities.append(density)
        
        # Calculate total "weighted width" and allocate points
        region_widths = [x_if[k+1] - x_if[k] for k in range(n_regions)]
        weighted_widths = [region_densities[k] * region_widths[k] for k in range(n_regions)]
        total_weighted = sum(weighted_widths)
        
        region_points = []
        allocated = 0
        
        for k in range(n_regions):
            if k == n_regions - 1:
                pts = n - allocated
            else:
                pts = max(3, int(n * weighted_widths[k] / total_weighted))
            region_points.append(pts)
            allocated += pts
        
        # Generate grid segments
        x_parts = []
        ds = self.delta_s
        corner_trace = min(3 * ds, self.w / 4)
        corner_gap = min(3 * ds, self.gap / 4)
        ncorner = max(3, self.skin_cells // 4)
        
        for k in range(n_regions):
            x0, x1 = x_if[k], x_if[k + 1]
            npts = region_points[k]
            
            # Determine grading strategy based on region type
            
            # TRACE REGION - fine mesh with corner grading
            if x0 >= self.x_tr_l - 1e-15 and x1 <= self.x_tr_r + 1e-15:
                nc = min(ncorner, npts // 3)
                n_mid = max(3, npts - 2 * nc)
                
                corner_left = min(corner_trace, (x1 - x0) / 3)
                corner_right = min(corner_trace, (x1 - x0) / 3)
                
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
            
            # GAP REGIONS - grade toward both interfaces
            elif ((self.x_gap_l <= x0 and x1 <= self.x_tr_l + 1e-15) or
                  (self.x_tr_r - 1e-15 <= x0 and x1 <= self.x_gap_r)):
                seg = self._smooth_transition(x0, x1, npts,
                                             curve_end="both", beta=4.0)
            
            # GROUND REGIONS - grade toward gap (inner edge)
            elif ((self.via_x_left_inner < x0 and x1 <= self.x_gap_l + 1e-15) or
                  (self.x_gap_r - 1e-15 <= x0 and x1 < self.via_x_right_inner)):
                # Grade toward trace (inner edge has higher field)
                if x1 <= self.x_gap_l + 1e-15:
                    # Left ground: dense at right (toward trace)
                    seg = self._smooth_transition(x0, x1, npts,
                                                 curve_end="end", beta=3.0)
                else:
                    # Right ground: dense at left (toward trace)
                    seg = self._smooth_transition(x0, x1, npts,
                                                 curve_end="start", beta=3.0)
            
            # OUTER REGIONS (domain edge to via) - uniform, no grading (low field)
            else:
                end = "end" if x0 < self.x_tr_l else "start"
                seg = self._smooth_transition(x0, x1, npts,
                                             curve_end=end, beta=1.5)
            
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
    plots = False
    if len(sys.argv) > 1:
        if "plot" in sys.argv[1]:
            plots = True
    solve_gcpw(plots)
