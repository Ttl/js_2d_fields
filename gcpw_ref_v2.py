#!/usr/bin/env python
import sys
import numpy as np
from field_solver_ref import FieldSolver2D, epsilon_0, mu_0, pi
from mesher import Dielectric, Conductor, Mesher

# ============================================================================
# GROUNDED COPLANAR WAVEGUIDE SOLVER V2
# ============================================================================

class GroundedCPWSolver2D(FieldSolver2D):
    """
    Grounded Coplanar Waveguide (GCPW) transmission line solver using the new Mesher class.

    This version uses lists of dielectrics and conductors for more flexible
    geometry specification.
    """

    def __init__(self, substrate_height, trace_width, trace_thickness,
                 gap, top_gnd_width, via_gap,
                 gnd_thickness=35e-6, epsilon_r=4.5, tan_delta=0.02,
                 sigma_diel=0.0, sigma_cond=5.8e7, epsilon_r_top=1,
                 air_top=None, air_side=None, freq=1e9,
                 nx=300, ny=300, boundaries=None,
                 use_sm=False, sm_t_sub=20e-6, sm_t_trace=20e-6,
                 sm_t_side=20e-6, sm_er=3.5, sm_tand=0.02):

        # Store parameters
        self.h = substrate_height
        self.w = trace_width
        self.t = trace_thickness
        self.gap = gap  # Gap to top ground
        self.top_gnd_w = top_gnd_width
        self.via_gap = via_gap  # Edge-to-edge distance: signal edge to via edge
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

        # Store air parameters
        self.air_top = air_top
        self.air_side = air_side

        # Domain sizing
        self.active_width = self.w + 2 * max(self.gap + self.top_gnd_w, self.via_gap + self.gap)

        if air_side is None:
            self.domain_width = self.active_width * 2
        else:
            self.domain_width = self.active_width + 2 * air_side

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
        # Y-coordinates
        self.y_gnd_bot_start = 0.0
        self.y_gnd_bot_end = self.t_gnd
        self.y_sub_start = self.y_gnd_bot_end
        self.y_sub_end = self.y_sub_start + self.h
        self.y_trace_start = self.y_sub_end
        self.y_trace_end = self.y_trace_start + self.t

        # Solder mask extents
        self.y_sm_sub_end = self.y_sub_end + self.sm_t_sub
        self.y_sm_trace_end = self.y_trace_end + self.sm_t_trace

        self.y_top_start = self.y_trace_end if not self.use_sm else max(self.y_sm_sub_end, self.y_sm_trace_end)

        # Top air/dielectric region
        if self.air_top is None:
            self.top_dielectric_h = self.h * 10
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

        # X-coordinates
        cx = self.domain_width / 2
        self.x_tr_l = cx - self.w / 2
        self.x_tr_r = cx + self.w / 2

        # Top Grounds
        self.x_gap_l = self.x_tr_l - self.gap
        self.x_gap_r = self.x_tr_r + self.gap

        # Vias (Single via on each side)
        self.via_x_left_inner = self.x_tr_l - self.gap - self.via_gap
        self.via_x_right_inner = self.x_tr_r + self.gap + self.via_gap

    def _build_geometry_lists(self):
        """Build lists of dielectrics and conductors from parameters."""
        dielectrics = []
        conductors = []

        # --- DIELECTRICS ---

        # Substrate (full width)
        dielectrics.append(Dielectric(
            0, self.y_sub_start,
            self.domain_width, self.h,
            self.er, self.tan_delta
        ))

        # Top air/dielectric region
        dielectrics.append(Dielectric(
            0, self.y_top_start,
            self.domain_width, self.top_dielectric_h,
            self.er_top, 0.0
        ))

        # Solder mask regions (if enabled)
        if self.use_sm:
            # Solder mask over substrate in gaps (Zone A)
            # Left gap
            if self.x_gap_l > 0:
                gap_width_left = self.x_tr_l - self.x_gap_l
                dielectrics.append(Dielectric(
                    self.x_gap_l, self.y_sub_end,
                    gap_width_left, self.sm_t_sub,
                    self.sm_er, self.sm_tand
                ))

            # Right gap
            if self.x_gap_r < self.domain_width:
                gap_width_right = self.x_gap_r - self.x_tr_r
                dielectrics.append(Dielectric(
                    self.x_tr_r, self.y_sub_end,
                    gap_width_right, self.sm_t_sub,
                    self.sm_er, self.sm_tand
                ))

            # Solder mask on top of signal trace (Zone B)
            dielectrics.append(Dielectric(
                self.x_tr_l, self.y_trace_end,
                self.w, self.sm_t_trace,
                self.sm_er, self.sm_tand
            ))

            # Solder mask on top of left ground
            if self.x_gap_l > 0:
                dielectrics.append(Dielectric(
                    0, self.y_trace_end,
                    self.x_gap_l, self.sm_t_trace,
                    self.sm_er, self.sm_tand
                ))

            # Solder mask on top of right ground
            if self.x_gap_r < self.domain_width:
                dielectrics.append(Dielectric(
                    self.x_gap_r, self.y_trace_end,
                    self.domain_width - self.x_gap_r, self.sm_t_trace,
                    self.sm_er, self.sm_tand
                ))

            # Solder mask on sidewalls (Zone C)
            # Left edge of signal trace
            if self.x_tr_l - self.sm_t_side >= 0:
                dielectrics.append(Dielectric(
                    self.x_tr_l - self.sm_t_side, self.y_sub_end,
                    self.sm_t_side, self.t + self.sm_t_trace,
                    self.sm_er, self.sm_tand
                ))

            # Right edge of signal trace
            if self.x_tr_r + self.sm_t_side <= self.domain_width:
                dielectrics.append(Dielectric(
                    self.x_tr_r, self.y_sub_end,
                    self.sm_t_side, self.t + self.sm_t_trace,
                    self.sm_er, self.sm_tand
                ))

            # Right edge of left ground (inner edge)
            if self.x_gap_l + self.sm_t_side <= self.x_tr_l:
                dielectrics.append(Dielectric(
                    self.x_gap_l, self.y_sub_end,
                    self.sm_t_side, self.t + self.sm_t_trace,
                    self.sm_er, self.sm_tand
                ))

            # Left edge of right ground (inner edge)
            if self.x_gap_r - self.sm_t_side >= self.x_tr_r:
                dielectrics.append(Dielectric(
                    self.x_gap_r - self.sm_t_side, self.y_sub_end,
                    self.sm_t_side, self.t + self.sm_t_trace,
                    self.sm_er, self.sm_tand
                ))

        # --- CONDUCTORS ---

        # Bottom ground plane (full width)
        if self.t_gnd > 0:
            conductors.append(Conductor(
                0, 0,
                self.domain_width, self.t_gnd,
                is_signal=False
            ))

        # Vias (connect bottom ground to top grounds)
        # Left via
        if self.via_x_left_inner > 0:
            conductors.append(Conductor(
                0, self.y_sub_start,
                self.via_x_left_inner, self.h,
                is_signal=False
            ))

        # Right via
        if self.via_x_right_inner < self.domain_width:
            conductors.append(Conductor(
                self.via_x_right_inner, self.y_sub_start,
                self.domain_width - self.via_x_right_inner, self.h,
                is_signal=False
            ))

        # Top ground planes (on sides)
        # Left ground
        if self.x_gap_l > 0:
            conductors.append(Conductor(
                0, self.y_trace_start,
                self.x_gap_l, self.t,
                is_signal=False
            ))

        # Right ground
        if self.x_gap_r < self.domain_width:
            conductors.append(Conductor(
                self.x_gap_r, self.y_trace_start,
                self.domain_width - self.x_gap_r, self.t,
                is_signal=False
            ))

        # Signal trace
        conductors.append(Conductor(
            self.x_tr_l, self.y_trace_start,
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

    def _setup_geometry(self):
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
        for cond in conductors:
            for i, yc in enumerate(self.y):
                if cond.y_min - tol <= yc <= cond.y_max + tol:
                    for j, xc in enumerate(self.x):
                        if cond.x_min - tol <= xc <= cond.x_max + tol:
                            if cond.is_signal:
                                self.signal_mask[i, j] = True
                            else:
                                self.ground_mask[i, j] = True
                            self.V[i, j] = cond.voltage

        # Finalize conductor mask
        self.conductor_mask = self.signal_mask | self.ground_mask


# ============================================================================
# TEST FUNCTION
# ============================================================================

def solve_gcpw(plots=True):
    print("Solving GCPW (v2 with Mesher)...")
    #solver = GroundedCPWSolver2D(
    #    substrate_height=1.6e-3,
    #    trace_width=0.3e-3,
    #    trace_thickness=35e-6,
    #    gap=0.15e-3,
    #    top_gnd_width=5e-3,
    #    via_gap=0.5e-3,
    #    gnd_thickness=35e-6,
    #    epsilon_r=4.5,
    #    tan_delta=0.02,
    #    sigma_cond=5.8e7,
    #    freq=1e9,
    #    nx=10, ny=10,
    #    use_sm=False,
    #    boundaries=["open", "open", "open", "gnd"]
    #)
    solver = GroundedCPWSolver2D(
        substrate_height=1.6e-3,
        trace_width=1.6e-3,
        trace_thickness=35e-6,
        gap=1e-3,
        top_gnd_width=0.5e-3,
        via_gap=1.6e-3,
        gnd_thickness=35e-6,
        epsilon_r=4.5,
        tan_delta=0.02,
        sigma_cond=5.8e7,
        freq=1e9,
        nx=30, ny=30,
        use_sm=False,
        boundaries=["open", "open", "open", "gnd"]
    )

    #Z0, eps_eff, C, C0 = solver.calculate_parameters()
    #Ex, Ey = solver.compute_fields()
    Z0, eps_eff, C, C0 = solver.solve_adaptive()
    Ex, Ey = solver.Ex, solver.Ey

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
        solver.plot_geometry()
        solver.plot(Ex, Ey)


if __name__ == "__main__":
    plots = False
    if len(sys.argv) > 1:
        if "plot" in sys.argv[1]:
            plots = True
    solve_gcpw(plots)
