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
                 conductors, dielectrics, symmetric=False):
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
        self.symmetric = symmetric

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
        """Create fine mesh for conductor region with symmetric distribution about center."""
        center = (start + end) / 2
        length = end - start

        if npts < 3:
            # Always include center for small npts
            return np.array([start, center, end])

        # Generate symmetric mesh by creating left half and mirroring
        # Always include boundaries and center
        mesh_points = [start, end, center]

        # Calculate how many additional points we need (excluding start, end, center)
        n_additional = max(0, npts - 3)

        if n_additional > 0:
            # Divide additional points between left and right halves
            n_half = n_additional // 2

            if n_half > 0:
                # Generate points in left half using smooth transition
                # Points go from start to center (excluding both endpoints)
                left_half_length = length / 2

                # Create graded distribution in left half
                xi = np.linspace(0, 1, n_half + 2)[1:-1]  # Exclude 0 and 1
                # Use tanh grading for finer mesh near edges
                beta = 2.0
                eta = np.tanh(beta * (xi - 0.5)) / np.tanh(beta * 0.5) / 2 + 0.5

                # Map to actual coordinates in left half
                left_points = start + eta * left_half_length
                mesh_points.extend(left_points)

                # Mirror to right half
                right_points = 2 * center - left_points
                mesh_points.extend(right_points)

            # If we have an odd number of additional points, add one more near center
            if n_additional % 2 == 1:
                # Add a point slightly to the left of center
                offset = length / (4 * npts)
                mesh_points.append(center - offset)
                mesh_points.append(center + offset)

        return np.sort(np.array(mesh_points))

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
                            beta_val = 2.0
                        elif abs(i0 - cond.x_max) < 1e-12:
                            end_curve = 'start'  # Dense toward conductor
                            beta_val = 2.0
                    else:
                        if abs(i1 - cond.y_min) < 1e-12:
                            end_curve = 'end'  # Dense toward conductor
                            beta_val = 2.0
                        elif abs(i0 - cond.y_max) < 1e-12:
                            end_curve = 'start'  # Dense toward conductor
                            beta_val = 2.0

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

        # Add center points for all conductors and dielectrics
        center_points = []

        # Add conductor centers
        for cond in self.conductors:
            if axis == 'x':
                center = (cond.x_min + cond.x_max) / 2
            else:
                center = (cond.y_min + cond.y_max) / 2

            # Only add if not at domain boundaries
            if 0 < center < domain_size:
                center_points.append(center)

        # Add dielectric centers
        for diel in self.dielectrics:
            if axis == 'x':
                center = (diel.x_min + diel.x_max) / 2
            else:
                center = (diel.y_min + diel.y_max) / 2

            # Only add if not at domain boundaries
            if 0 < center < domain_size:
                center_points.append(center)

        # Add center points to mesh if they don't already exist
        for center in center_points:
            min_dist = np.min(np.abs(mesh - center))
            # Only add if there's no point very close to it
            if min_dist > domain_size / 1000:
                mesh = np.sort(np.append(mesh, center))

        # Add boundary lines adjacent to all conductor edges
        # Add lines on BOTH sides of each conductor edge (inside and outside)
        boundary_offset = min(self.skin_depth * 3, domain_size / 200)
        boundary_lines = []

        for cond in self.conductors:
            if axis == 'x':
                edges = [cond.x_min, cond.x_max]
                cond_min, cond_max = cond.x_min, cond.x_max
            else:
                edges = [cond.y_min, cond.y_max]
                cond_min, cond_max = cond.y_min, cond.y_max

            for edge in edges:
                # Skip domain boundaries
                if abs(edge) < 1e-15 or abs(edge - domain_size) < 1e-15:
                    continue

                # Determine if this is the left/bottom or right/top edge
                is_left_edge = abs(edge - cond_min) < 1e-15

                if is_left_edge:
                    # Left/bottom edge: add line outside (to the left) and inside (to the right)
                    outside_line = edge - boundary_offset
                    inside_line = edge + boundary_offset
                else:
                    # Right/top edge: add line inside (to the left) and outside (to the right)
                    inside_line = edge - boundary_offset
                    outside_line = edge + boundary_offset

                # Add lines if they're within the domain and within/adjacent to the conductor
                if 0 < outside_line < domain_size:
                    boundary_lines.append(outside_line)
                if cond_min < inside_line < cond_max:
                    boundary_lines.append(inside_line)

        # Add boundary lines to mesh if they don't already exist
        for line in boundary_lines:
            # Only add if there's no existing point very close to it
            if np.min(np.abs(mesh - line)) > boundary_offset / 3:
                mesh = np.sort(np.append(mesh, line))

        # Check if geometry is symmetric and enforce symmetry
        if self.symmetric:
            is_symmetric = self._check_symmetry(axis)
            if is_symmetric:
                mesh = self._enforce_symmetry(mesh, domain_size)

        # Remove duplicate or near-duplicate points
        # This is critical to avoid division by zero in the solver
        mesh_unique = [mesh[0]]
        min_spacing = domain_size * 1e-10  # Minimum allowed spacing

        for i in range(1, len(mesh)):
            if mesh[i] - mesh_unique[-1] > min_spacing:
                mesh_unique.append(mesh[i])

        return np.array(mesh_unique)

    def _check_symmetry(self, axis):
        """Check if geometry is symmetric about the center line."""
        tol = 1e-12

        if axis == 'x':
            center = self.domain_width / 2
            # Check if all conductors are symmetric about center
            for cond in self.conductors:
                x_min_mirror = center - (cond.x_max - center)
                x_max_mirror = center - (cond.x_min - center)

                # Check if there's a matching conductor
                found_match = False
                for other in self.conductors:
                    if (abs(other.x_min - x_min_mirror) < tol and
                        abs(other.x_max - x_max_mirror) < tol and
                        other.is_signal == cond.is_signal):
                        found_match = True
                        break

                # Also check if conductor is centered
                if abs(cond.x_min + cond.x_max - 2 * center) < tol:
                    found_match = True

                if not found_match:
                    return False

            return True
        else:
            return False

    def _enforce_symmetry(self, mesh, domain_size):
        """Enforce symmetry in mesh by averaging symmetric pairs."""
        center = domain_size / 2
        tol = 1e-12

        symmetric_mesh = []

        # Process points from left to right
        used = np.zeros(len(mesh), dtype=bool)

        for i, point in enumerate(mesh):
            if used[i]:
                continue

            # Check if point is at center
            if abs(point - center) < tol:
                symmetric_mesh.append(center)
                used[i] = True
                continue

            # Find mirror point
            mirror_pos = 2 * center - point
            mirror_idx = None
            min_dist = float('inf')

            for j in range(len(mesh)):
                if used[j]:
                    continue
                dist = abs(mesh[j] - mirror_pos)
                if dist < min_dist:
                    min_dist = dist
                    mirror_idx = j

            if mirror_idx is not None and min_dist < domain_size / 100:
                # Found mirror - average the positions
                avg_pos = center + (point - center)  # Keep left point as-is
                mirror_avg_pos = 2 * center - avg_pos

                symmetric_mesh.append(avg_pos)
                symmetric_mesh.append(mirror_avg_pos)
                used[i] = True
                used[mirror_idx] = True
            else:
                # No mirror found - create one
                symmetric_mesh.append(point)
                if 0 < mirror_pos < domain_size:
                    symmetric_mesh.append(mirror_pos)
                used[i] = True

        return np.sort(np.unique(symmetric_mesh))


