#!/usr/bin/env python
import sys
import numpy as np
from field_solver_ref import FieldSolver2D, epsilon_0, mu_0, pi
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
                 top_diel_tand=0.0):

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
                self.domain_width, self.gnd_cut_sub_h + self.h,
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
                    self.sm_t_side, self.sm_t_trace,
                    self.sm_er, self.sm_tand
                ))

            # Solder mask on right side of trace
            xsr = xr + self.sm_t_side
            if xsr <= self.domain_width:
                dielectrics.append(Dielectric(
                    xr, self.y_trace_start,
                    self.sm_t_side, self.sm_t_trace,
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

        # Signal trace
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

def solve_microstrip(plots=True):
    print("Solving Microstrip (v2 with Mesher)...")
    solver = MicrostripSolver2D(
        substrate_height=1.6e-3,
        trace_width=3e-3,
        trace_thickness=35e-6,
        gnd_thickness=35e-6,
        epsilon_r=4.5,
        tan_delta=0.02,
        sigma_cond=5.8e7,
        freq=1e9,
        nx=1, ny=1,
        use_sm=False,
        boundaries=["open", "open", "open", "gnd"]
    )

    #Z0, eps_eff, C, C0 = solver.calculate_parameters()
    #Ex, Ey = solver.compute_fields()
    Z0, eps_eff, C, C0, Ex, Ey = solver.solve_adaptive(param_tol=0.001)

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


if __name__ == "__main__":
    plots = False
    if len(sys.argv) > 1:
        if "plot" in sys.argv[1]:
            plots = True
    solve_microstrip(plots)
