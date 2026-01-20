import numpy as np

# ============================================================================
# GEOMETRY CLASSES
# ============================================================================

class Dielectric:
    """Represents a rectangular dielectric region."""
    def __init__(self, x, y, width, height, epsilon_r, tan_delta=0.0):
        """
        Parameters:
        -----------
        x, y : float
            Bottom-left corner coordinates
        width, height : float
            Dimensions of the rectangle
        epsilon_r : float
            Relative permittivity
        tan_delta : float
            Loss tangent
        """
        self.x = x
        self.y = y
        self.width = width
        self.height = height
        self.epsilon_r = epsilon_r
        self.tan_delta = tan_delta

    @property
    def x_min(self):
        return self.x

    @property
    def x_max(self):
        return self.x + self.width

    @property
    def y_min(self):
        return self.y

    @property
    def y_max(self):
        return self.y + self.height


class Conductor:
    """Represents a rectangular conductor region."""
    def __init__(self, x, y, width, height, is_signal=False, voltage=0.0):
        """
        Parameters:
        -----------
        x, y : float
            Bottom-left corner coordinates
        width, height : float
            Dimensions of the rectangle
        is_signal : bool
            True for signal conductor, False for ground
        voltage : float
            Voltage value (1.0 for signal, 0.0 for ground typically)
        """
        self.x = x
        self.y = y
        self.width = width
        self.height = height
        self.is_signal = is_signal
        self.voltage = voltage if voltage != 0.0 else (1.0 if is_signal else 0.0)

    @property
    def x_min(self):
        return self.x

    @property
    def x_max(self):
        return self.x + self.width

    @property
    def y_min(self):
        return self.y

    @property
    def y_max(self):
        return self.y + self.height


# ============================================================================
# MESHER CLASS
# ============================================================================

class Mesher:
    """
    Generates adaptive graded mesh based on conductor and dielectric locations.

    The mesh is graded such that:
    - More mesh lines near signal conductors (high field regions)
    - More mesh lines near conductor surfaces (for accurate loss calculation)
    - More mesh lines at dielectric interfaces
    - Fewer mesh lines in far-field regions
    """

    def __init__(self, domain_width, domain_height, nx, ny, skin_depth,
                 conductors, dielectrics):
        """
        Parameters:
        -----------
        domain_width, domain_height : float
            Physical domain dimensions
        nx, ny : int
            Approximate number of mesh points in each direction
        skin_depth : float
            Skin depth for conductor meshing
        conductors : list of Conductor
            List of conductor objects
        dielectrics : list of Dielectric
            List of dielectric objects
        """
        self.domain_width = domain_width
        self.domain_height = domain_height
        self.nx = nx
        self.ny = ny
        self.skin_depth = skin_depth
        self.conductors = conductors
        self.dielectrics = dielectrics

        # Calculate corner mesh parameters
        # With 100 mesh lines, we want about 2 extra lines near conductor edges
        self.ncorner = max(2, nx // 40)  # About 2-3 lines for nx=100
        self.corner_size = min(3 * skin_depth, self._min_conductor_dimension() / 4)

    def _min_conductor_dimension(self):
        """Find minimum conductor dimension."""
        min_dim = float('inf')
        for cond in self.conductors:
            min_dim = min(min_dim, cond.width, cond.height)
        return min_dim if min_dim != float('inf') else 1e-3

    def _smooth_transition(self, start, end, n_points, curve_end='end', beta=4.0):
        """Create smooth transition with specified grading."""
        if n_points <= 1:
            return np.array([start, end])
        xi = np.linspace(0, 1, n_points)
        if curve_end == 'end':
            eta = np.tanh(beta * xi) / np.tanh(beta)
        elif curve_end == 'both':
            eta = (np.tanh(beta * (xi - 0.5)) / np.tanh(beta * 0.5) + 1) / 2
        else:  # 'start'
            eta = 1 - np.tanh(beta * (1 - xi)) / np.tanh(beta)
        return start + eta * (end - start)

    def _collect_interfaces_x(self):
        """Collect all critical x-coordinates (conductor edges, dielectric edges)."""
        x_if = [0.0, self.domain_width]

        for cond in self.conductors:
            x_if.extend([cond.x_min, cond.x_max])

        for diel in self.dielectrics:
            # Only add if not at domain boundaries
            if diel.x_min > 0:
                x_if.append(diel.x_min)
            if diel.x_max < self.domain_width:
                x_if.append(diel.x_max)

        return np.array(sorted(set(x_if)))

    def _collect_interfaces_y(self):
        """Collect all critical y-coordinates (conductor edges, dielectric edges)."""
        y_if = [0.0, self.domain_height]

        for cond in self.conductors:
            y_if.extend([cond.y_min, cond.y_max])

        for diel in self.dielectrics:
            # Only add if not at domain boundaries
            if diel.y_min > 0:
                y_if.append(diel.y_min)
            if diel.y_max < self.domain_height:
                y_if.append(diel.y_max)

        return np.array(sorted(set(y_if)))

    def _is_inside_conductor(self, x0, x1, y0, y1):
        """Check if region overlaps with any conductor."""
        tol = 1e-15
        for cond in self.conductors:
            if (x0 >= cond.x_min - tol and x1 <= cond.x_max + tol and
                y0 >= cond.y_min - tol and y1 <= cond.y_max + tol):
                return True
        return False

    def _region_weight_x(self, x0, x1):
        """Calculate mesh density weight for x-region."""
        tol = 1e-15

        # Check if region is inside a conductor
        for cond in self.conductors:
            if x0 >= cond.x_min - tol and x1 <= cond.x_max + tol:
                # Inside conductor - high weight for signal, moderate for ground
                return 10.0 if cond.is_signal else 5.0

        # Check if region is near a conductor
        min_dist_signal = float('inf')
        min_dist_ground = float('inf')
        for cond in self.conductors:
            dist = min(abs(x0 - cond.x_min), abs(x0 - cond.x_max),
                      abs(x1 - cond.x_min), abs(x1 - cond.x_max))
            if cond.is_signal:
                min_dist_signal = min(min_dist_signal, dist)
            else:
                min_dist_ground = min(min_dist_ground, dist)

        # Weight based on distance from signal conductors (highest priority)
        if min_dist_signal < self.skin_depth * 5:
            return 5.0  # Very close to signal conductor
        elif min_dist_signal < self.skin_depth * 20:
            return 2.5  # Near signal conductor
        elif min_dist_signal < self.skin_depth * 50:
            return 1.0  # Moderate distance
        # Check ground conductor proximity (for cutout case)
        elif min_dist_ground < self.skin_depth * 5:
            return 1.5  # Near ground conductor edges
        else:
            return 0.2  # Far field - minimize mesh density

    def _region_weight_y(self, y0, y1):
        """Calculate mesh density weight for y-region."""
        tol = 1e-15

        # Check if region is inside a conductor
        for cond in self.conductors:
            if y0 >= cond.y_min - tol and y1 <= cond.y_max + tol:
                # Inside conductor - very high weight for signal
                return 20.0 if cond.is_signal else 6.0

        # Check if region is near a conductor
        min_dist_signal = float('inf')
        min_dist_ground = float('inf')
        for cond in self.conductors:
            dist = min(abs(y0 - cond.y_min), abs(y0 - cond.y_max),
                      abs(y1 - cond.y_min), abs(y1 - cond.y_max))
            if cond.is_signal:
                min_dist_signal = min(min_dist_signal, dist)
            else:
                min_dist_ground = min(min_dist_ground, dist)

        # Check for dielectric interfaces
        at_interface = False
        for diel in self.dielectrics:
            if abs(y0 - diel.y_min) < tol or abs(y0 - diel.y_max) < tol or \
               abs(y1 - diel.y_min) < tol or abs(y1 - diel.y_max) < tol:
                at_interface = True
                break

        # Weight based on distance from signal (highest priority)
        if min_dist_signal < self.skin_depth * 5:
            return 6.0  # Very close to signal conductor
        elif min_dist_signal < self.skin_depth * 20:
            return 3.0  # Near signal conductor
        elif at_interface and min_dist_signal < self.skin_depth * 50:
            return 1.5  # At dielectric interface near signal
        elif min_dist_signal < self.skin_depth * 50:
            return 1.0  # Moderate distance from signal
        # Check ground conductor proximity (for cutout case)
        elif min_dist_ground < self.skin_depth * 5:
            return 1.5  # Near ground conductor edges
        else:
            return 0.15  # Far field - minimize mesh density

    def _mesh_conductor_region(self, start, end, npts, direction='x'):
        """Create fine mesh for conductor region with corner grading."""
        if npts < 3:
            return np.linspace(start, end, max(npts, 2))

        length = end - start
        # Reduce corner meshing - only use 1-2 points per corner for nx=100
        nc = max(1, min(max(1, self.ncorner // 2), npts // 4))
        n_mid = max(npts - 2 * nc, npts // 2)

        corner = min(self.corner_size, length / 4)

        # Generate segments
        parts = []

        # Left/bottom corner
        if nc > 0 and corner > 1e-15:
            parts.append(self._smooth_transition(start, start + corner, nc,
                                                curve_end='start', beta=3.0))
        else:
            parts.append(np.array([start]))

        # Middle section - use uniform distribution for better coverage
        mid_start = start + corner
        mid_end = max(end - corner, mid_start)
        if n_mid > 0:
            parts.append(self._smooth_transition(mid_start, mid_end, n_mid,
                                                curve_end='end', beta=3.0)[1:])


        # Right/top corner
        if nc > 0 and corner > 1e-15 and mid_end < end - 1e-15:
            parts.append(self._smooth_transition(mid_end, end, nc,
                                                curve_end='end', beta=3.0)[1:])
        elif mid_end < end - 1e-15:
            parts.append(np.array([end]))

        return np.concatenate(parts)

    def generate_mesh(self):
        """Generate the x and y mesh arrays."""
        x = self._generate_axis_mesh('x')
        y = self._generate_axis_mesh('y')
        return x, y

    def _generate_axis_mesh(self, axis):
        """Generate mesh for specified axis ('x' or 'y')."""
        if axis == 'x':
            interfaces = self._collect_interfaces_x()
            n_points = self.nx
            domain_size = self.domain_width
            weight_func = self._region_weight_x
        else:
            interfaces = self._collect_interfaces_y()
            n_points = self.ny
            domain_size = self.domain_height
            weight_func = self._region_weight_y

        n_regions = len(interfaces) - 1

        # Calculate weights for each region
        region_weights = []
        for k in range(n_regions):
            i0, i1 = interfaces[k], interfaces[k + 1]
            width = i1 - i0
            weight = weight_func(i0, i1)
            region_weights.append(weight * width)

        # Allocate points based on weights
        # First, ensure minimum points in SIGNAL conductors (not ground planes)
        MIN_CONDUCTOR_POINTS = 5
        region_points = []
        reserved_points = 0
        non_conductor_weight = 0

        for k in range(n_regions):
            i0, i1 = interfaces[k], interfaces[k + 1]
            # Check if this region is inside a SIGNAL conductor
            is_signal_conductor = False
            for cond in self.conductors:
                if not cond.is_signal:
                    continue  # Skip ground planes
                tol = 1e-15
                if axis == 'x':
                    if i0 >= cond.x_min - tol and i1 <= cond.x_max + tol:
                        is_signal_conductor = True
                        break
                else:
                    if i0 >= cond.y_min - tol and i1 <= cond.y_max + tol:
                        is_signal_conductor = True
                        break

            if is_signal_conductor:
                region_points.append(MIN_CONDUCTOR_POINTS)
                reserved_points += MIN_CONDUCTOR_POINTS
            else:
                region_points.append(0)  # Will be allocated later
                non_conductor_weight += region_weights[k]

        # Allocate remaining points based on weights (only for non-conductor regions)
        remaining_points = n_points - reserved_points
        allocated = reserved_points

        for k in range(n_regions):
            if region_points[k] > 0:
                # Already has minimum conductor points
                continue

            if k == n_regions - 1 and allocated < n_points:
                pts = n_points - allocated
            else:
                if non_conductor_weight > 0:
                    pts = max(5, int(remaining_points * region_weights[k] / non_conductor_weight))
                else:
                    pts = 5
            region_points[k] = pts
            allocated += pts

        # Generate mesh segments
        mesh_parts = []

        for k in range(n_regions):
            i0, i1 = interfaces[k], interfaces[k + 1]
            npts = region_points[k]

            # Check if this region is inside a SIGNAL conductor
            is_signal_conductor = False
            for cond in self.conductors:
                if not cond.is_signal:
                    continue  # Skip ground planes for special meshing
                tol = 1e-15
                if axis == 'x':
                    if i0 >= cond.x_min - tol and i1 <= cond.x_max + tol:
                        is_signal_conductor = True
                        break
                else:
                    if i0 >= cond.y_min - tol and i1 <= cond.y_max + tol:
                        is_signal_conductor = True
                        break

            if is_signal_conductor:
                # Use fine conductor meshing with corner grading for signal
                seg = self._mesh_conductor_region(i0, i1, npts, axis)
            else:
                # Determine grading strategy for dielectric regions
                # Grade toward conductors
                end_curve = 'both'
                beta_val = 1.0

                # Check if at domain boundary
                if abs(i0) < 1e-15:
                    end_curve = 'end'  # Dense toward right/top
                    beta_val = 2.0
                elif abs(i1 - domain_size) < 1e-15:
                    end_curve = 'start'  # Dense toward left/bottom
                    beta_val = 2.0

                # Check proximity to conductors
                for cond in self.conductors:
                    if axis == 'x':
                        if abs(i1 - cond.x_min) < 1e-12:
                            end_curve = 'end'  # Dense toward conductor
                            beta_val = 1.0
                        elif abs(i0 - cond.x_max) < 1e-12:
                            end_curve = 'start'  # Dense toward conductor
                            beta_val = 1.0
                    else:
                        if abs(i1 - cond.y_min) < 1e-12:
                            end_curve = 'end'  # Dense toward conductor
                            beta_val = 1.0
                        elif abs(i0 - cond.y_max) < 1e-12:
                            end_curve = 'start'  # Dense toward conductor
                            beta_val = 1.0

                seg = self._smooth_transition(i0, i1, npts,
                                            curve_end=end_curve, beta=beta_val)

            # Avoid duplicating interface points
            if k > 0:
                seg = seg[1:]
            mesh_parts.append(seg)

        mesh = np.concatenate(mesh_parts)

        # Ensure exact interface locations
        for interface in interfaces:
            if np.min(np.abs(mesh - interface)) > 1e-12:
                mesh = np.sort(np.append(mesh, interface))

        # Remove duplicate or near-duplicate points
        # This is critical to avoid division by zero in the solver
        mesh_unique = [mesh[0]]
        min_spacing = domain_size * 1e-10  # Minimum allowed spacing

        for i in range(1, len(mesh)):
            if mesh[i] - mesh_unique[-1] > min_spacing:
                mesh_unique.append(mesh[i])

        return np.array(mesh_unique)


