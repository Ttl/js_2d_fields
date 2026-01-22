// ============================================================================
// GEOMETRY CLASSES
// ============================================================================

class Dielectric {
    /**
     * Represents a rectangular dielectric region.
     * @param {number} x - Bottom-left corner x-coordinate
     * @param {number} y - Bottom-left corner y-coordinate
     * @param {number} width - Width of the rectangle
     * @param {number} height - Height of the rectangle
     * @param {number} epsilon_r - Relative permittivity
     * @param {number} tan_delta - Loss tangent (default: 0.0)
     */
    constructor(x, y, width, height, epsilon_r, tan_delta = 0.0) {
        this.x = x;
        this.y = y;
        this.width = width;
        this.height = height;
        this.epsilon_r = epsilon_r;
        this.tan_delta = tan_delta;
    }

    get x_min() { return this.x; }
    get x_max() { return this.x + this.width; }
    get y_min() { return this.y; }
    get y_max() { return this.y + this.height; }
}

class Conductor {
    /**
     * Represents a rectangular conductor region.
     * @param {number} x - Bottom-left corner x-coordinate
     * @param {number} y - Bottom-left corner y-coordinate
     * @param {number} width - Width of the rectangle
     * @param {number} height - Height of the rectangle
     * @param {boolean} is_signal - True for signal conductor, False for ground
     * @param {number} polarity - Signal polarity: +1 (positive), -1 (negative), 0 (ground)
     */
    constructor(x, y, width, height, is_signal = false, polarity = 0) {
        this.x = x;
        this.y = y;
        this.width = width;
        this.height = height;
        this.is_signal = is_signal;
        this.polarity = is_signal ? (polarity || 1) : 0;
    }

    get x_min() { return this.x; }
    get x_max() { return this.x + this.width; }
    get y_min() { return this.y; }
    get y_max() { return this.y + this.height; }
}

// ============================================================================
// MESHER CLASS
// ============================================================================

class Mesher {
    /**
     * Generates adaptive graded mesh based on conductor and dielectric locations.
     * @param {number} domain_width - Physical domain width
     * @param {number} domain_height - Physical domain height
     * @param {number} nx - Approximate number of mesh points in x-direction
     * @param {number} ny - Approximate number of mesh points in y-direction
     * @param {number} skin_depth - Skin depth for conductor meshing
     * @param {Array<Conductor>} conductors - List of conductor objects
     * @param {Array<Dielectric>} dielectrics - List of dielectric objects
     * @param {boolean} symmetric - Whether to enforce symmetry (default: false)
     * @param {number} x_min - Minimum x-coordinate (default: 0)
     * @param {number} x_max - Maximum x-coordinate (default: domain_width)
     * @param {number} y_min - Minimum y-coordinate (default: 0)
     * @param {number} y_max - Maximum y-coordinate (default: domain_height)
     */
    constructor(domain_width, domain_height, nx, ny, skin_depth, conductors, dielectrics, symmetric = false, x_min = 0, x_max = null, y_min = 0, y_max = null) {
        this.domain_width = domain_width;
        this.domain_height = domain_height;
        this.nx = nx;
        this.ny = ny;
        this.skin_depth = skin_depth;
        this.conductors = conductors;
        this.dielectrics = dielectrics;
        this.symmetric = symmetric;
        this.x_min = x_min;
        this.x_max = x_max !== null ? x_max : domain_width;
        this.y_min = y_min;
        this.y_max = y_max !== null ? y_max : domain_height;

        // Calculate corner mesh parameters
        this.ncorner = Math.max(2, Math.floor(nx / 40));  // About 2-3 lines for nx=100
        this.corner_size = Math.min(3 * skin_depth, this._min_conductor_dimension() / 4);
    }

    _min_conductor_dimension() {
        let min_dim = Infinity;
        for (const cond of this.conductors) {
            min_dim = Math.min(min_dim, cond.width, cond.height);
        }
        return min_dim !== Infinity ? min_dim : 1e-3;
    }

    _smooth_transition(start, end, n_points, curve_end = 'end', beta = 4.0) {
        if (n_points <= 1) {
            return new Float64Array([start, end]);
        }

        const result = new Float64Array(n_points);
        const tanhBeta = Math.tanh(beta);
        const tanhBetaHalf = Math.tanh(beta * 0.5);

        for (let i = 0; i < n_points; i++) {
            const xi = i / (n_points - 1);
            let eta;

            if (curve_end === 'end') {
                eta = Math.tanh(beta * xi) / tanhBeta;
            } else if (curve_end === 'both') {
                eta = (Math.tanh(beta * (xi - 0.5)) / tanhBetaHalf + 1) / 2;
            } else {  // 'start'
                eta = 1 - Math.tanh(beta * (1 - xi)) / tanhBeta;
            }

            result[i] = start + eta * (end - start);
        }
        return result;
    }

    _collect_interfaces_x() {
        const x_if = new Set([this.x_min, this.x_max]);

        for (const cond of this.conductors) {
            x_if.add(cond.x_min);
            x_if.add(cond.x_max);
        }

        for (const diel of this.dielectrics) {
            if (diel.x_min > this.x_min) {
                x_if.add(diel.x_min);
            }
            if (diel.x_max < this.x_max) {
                x_if.add(diel.x_max);
            }
        }

        return Array.from(x_if).sort((a, b) => a - b);
    }

    _collect_interfaces_y() {
        const y_if = new Set([this.y_min, this.y_max]);

        for (const cond of this.conductors) {
            y_if.add(cond.y_min);
            y_if.add(cond.y_max);
        }

        for (const diel of this.dielectrics) {
            if (diel.y_min > this.y_min) {
                y_if.add(diel.y_min);
            }
            if (diel.y_max < this.y_max) {
                y_if.add(diel.y_max);
            }
        }

        return Array.from(y_if).sort((a, b) => a - b);
    }

    _region_weight_x(x0, x1) {
        const tol = 1e-15;

        // Check if region is inside a conductor
        for (const cond of this.conductors) {
            if (x0 >= cond.x_min - tol && x1 <= cond.x_max + tol) {
                return cond.is_signal ? 10.0 : 5.0;
            }
        }

        // Check if region is near a conductor
        let min_dist_signal = Infinity;
        let min_dist_ground = Infinity;

        for (const cond of this.conductors) {
            const dist = Math.min(
                Math.abs(x0 - cond.x_min), Math.abs(x0 - cond.x_max),
                Math.abs(x1 - cond.x_min), Math.abs(x1 - cond.x_max)
            );
            if (cond.is_signal) {
                min_dist_signal = Math.min(min_dist_signal, dist);
            } else {
                min_dist_ground = Math.min(min_dist_ground, dist);
            }
        }

        // Weight based on distance from signal conductors
        if (min_dist_signal < this.skin_depth * 5) {
            return 5.0;  // Very close to signal conductor
        } else if (min_dist_signal < this.skin_depth * 20) {
            return 2.5;  // Near signal conductor
        } else if (min_dist_signal < this.skin_depth * 50) {
            return 1.0;  // Moderate distance
        } else if (min_dist_ground < this.skin_depth * 5) {
            return 1.5;  // Near ground conductor edges
        } else {
            return 0.2;  // Far field
        }
    }

    _region_weight_y(y0, y1) {
        const tol = 1e-15;

        // Check if region is inside a conductor
        for (const cond of this.conductors) {
            if (y0 >= cond.y_min - tol && y1 <= cond.y_max + tol) {
                return cond.is_signal ? 20.0 : 6.0;
            }
        }

        // Check if region is near a conductor
        let min_dist_signal = Infinity;
        let min_dist_ground = Infinity;

        for (const cond of this.conductors) {
            const dist = Math.min(
                Math.abs(y0 - cond.y_min), Math.abs(y0 - cond.y_max),
                Math.abs(y1 - cond.y_min), Math.abs(y1 - cond.y_max)
            );
            if (cond.is_signal) {
                min_dist_signal = Math.min(min_dist_signal, dist);
            } else {
                min_dist_ground = Math.min(min_dist_ground, dist);
            }
        }

        // Check for dielectric interfaces
        let at_interface = false;
        for (const diel of this.dielectrics) {
            if (Math.abs(y0 - diel.y_min) < tol || Math.abs(y0 - diel.y_max) < tol ||
                Math.abs(y1 - diel.y_min) < tol || Math.abs(y1 - diel.y_max) < tol) {
                at_interface = true;
                break;
            }
        }

        // Weight based on distance from signal
        if (min_dist_signal < this.skin_depth * 5) {
            return 6.0;  // Very close to signal conductor
        } else if (min_dist_signal < this.skin_depth * 20) {
            return 3.0;  // Near signal conductor
        } else if (at_interface && min_dist_signal < this.skin_depth * 50) {
            return 1.5;  // At dielectric interface near signal
        } else if (min_dist_signal < this.skin_depth * 50) {
            return 1.0;  // Moderate distance from signal
        } else if (min_dist_ground < this.skin_depth * 5) {
            return 1.5;  // Near ground conductor edges
        } else {
            return 0.15;  // Far field
        }
    }

    _mesh_conductor_region(start, end, npts, direction = 'x') {
        const center = (start + end) / 2;
        const length = end - start;

        if (npts < 3) {
            return new Float64Array([start, center, end]);
        }

        const mesh_points = [start, end, center];
        const n_additional = Math.max(0, npts - 3);

        if (n_additional > 0) {
            const n_half = Math.floor(n_additional / 2);

            if (n_half > 0) {
                const left_half_length = length / 2;
                const beta = 2.0;
                const left_points = [];

                for (let i = 1; i <= n_half; i++) {
                    const xi = i / (n_half + 1);
                    const eta = Math.tanh(beta * (xi - 0.5)) / Math.tanh(beta * 0.5) / 2 + 0.5;
                    const pt = start + eta * left_half_length;
                    left_points.push(pt);
                    mesh_points.push(pt);
                    mesh_points.push(2 * center - pt);  // Mirror to right
                }
            }

            if (n_additional % 2 === 1) {
                const offset = length / (4 * npts);
                mesh_points.push(center - offset);
                mesh_points.push(center + offset);
            }
        }

        return Float64Array.from(mesh_points.sort((a, b) => a - b));
    }

    generate_mesh() {
        const x = this._generate_axis_mesh('x');
        const y = this._generate_axis_mesh('y');

        // Validate symmetry if requested
        if (this.symmetric) {
            this._validate_mesh_symmetry(x);
        }

        return [x, y];
    }

    _generate_axis_mesh(axis) {
        const interfaces = axis === 'x' ? this._collect_interfaces_x() : this._collect_interfaces_y();
        const n_points = axis === 'x' ? this.nx : this.ny;
        const domain_size = axis === 'x' ? this.domain_width : this.domain_height;
        const weight_func = axis === 'x' ?
            (a, b) => this._region_weight_x(a, b) :
            (a, b) => this._region_weight_y(a, b);

        const n_regions = interfaces.length - 1;

        // Calculate weights for each region
        const region_weights = [];
        for (let k = 0; k < n_regions; k++) {
            const i0 = interfaces[k];
            const i1 = interfaces[k + 1];
            const width = i1 - i0;
            const weight = weight_func(i0, i1);
            region_weights.push(weight * width);
        }

        // Allocate points - ensure minimum points in SIGNAL conductors
        const MIN_CONDUCTOR_POINTS = 5;
        const region_points = [];
        let reserved_points = 0;
        let non_conductor_weight = 0;

        for (let k = 0; k < n_regions; k++) {
            const i0 = interfaces[k];
            const i1 = interfaces[k + 1];
            let is_signal_conductor = false;

            for (const cond of this.conductors) {
                if (!cond.is_signal) continue;
                const tol = 1e-15;
                const cond_min = axis === 'x' ? cond.x_min : cond.y_min;
                const cond_max = axis === 'x' ? cond.x_max : cond.y_max;

                if (i0 >= cond_min - tol && i1 <= cond_max + tol) {
                    is_signal_conductor = true;
                    break;
                }
            }

            if (is_signal_conductor) {
                region_points.push(MIN_CONDUCTOR_POINTS);
                reserved_points += MIN_CONDUCTOR_POINTS;
            } else {
                region_points.push(0);
                non_conductor_weight += region_weights[k];
            }
        }

        // Allocate remaining points based on weights
        const remaining_points = n_points - reserved_points;
        let allocated = reserved_points;

        for (let k = 0; k < n_regions; k++) {
            if (region_points[k] > 0) continue;

            let pts;
            if (k === n_regions - 1 && allocated < n_points) {
                pts = n_points - allocated;
            } else {
                if (non_conductor_weight > 0) {
                    pts = Math.max(5, Math.floor(remaining_points * region_weights[k] / non_conductor_weight));
                } else {
                    pts = 5;
                }
            }
            region_points[k] = pts;
            allocated += pts;
        }

        // Generate mesh segments
        const mesh_parts = [];

        for (let k = 0; k < n_regions; k++) {
            const i0 = interfaces[k];
            const i1 = interfaces[k + 1];
            const npts = region_points[k];

            // Check if this region is inside a SIGNAL conductor
            let is_signal_conductor = false;
            for (const cond of this.conductors) {
                if (!cond.is_signal) continue;
                const tol = 1e-15;
                const cond_min = axis === 'x' ? cond.x_min : cond.y_min;
                const cond_max = axis === 'x' ? cond.x_max : cond.y_max;

                if (i0 >= cond_min - tol && i1 <= cond_max + tol) {
                    is_signal_conductor = true;
                    break;
                }
            }

            let seg;
            if (is_signal_conductor) {
                seg = this._mesh_conductor_region(i0, i1, npts, axis);
            } else {
                // Determine grading strategy
                let end_curve = 'both';
                let beta_val = 1.0;

                if (Math.abs(i0) < 1e-15) {
                    end_curve = 'end';
                    beta_val = 2.0;
                } else if (Math.abs(i1 - domain_size) < 1e-15) {
                    end_curve = 'start';
                    beta_val = 2.0;
                }

                // Check proximity to conductors
                for (const cond of this.conductors) {
                    const cond_min = axis === 'x' ? cond.x_min : cond.y_min;
                    const cond_max = axis === 'x' ? cond.x_max : cond.y_max;

                    if (Math.abs(i1 - cond_min) < 1e-12) {
                        end_curve = 'end';
                        beta_val = 2.0;
                    } else if (Math.abs(i0 - cond_max) < 1e-12) {
                        end_curve = 'start';
                        beta_val = 2.0;
                    }
                }

                seg = this._smooth_transition(i0, i1, npts, end_curve, beta_val);
            }

            if (k > 0) {
                seg = seg.slice(1);
            }
            mesh_parts.push(seg);
        }

        let mesh = this._concat_arrays(mesh_parts);

        // Ensure exact interface locations
        for (const interface_val of interfaces) {
            let min_dist = Infinity;
            for (const val of mesh) {
                min_dist = Math.min(min_dist, Math.abs(val - interface_val));
            }
            if (min_dist > 1e-12) {
                const temp = Array.from(mesh);
                temp.push(interface_val);
                mesh = Float64Array.from(temp.sort((a, b) => a - b));
            }
        }

        // Add center points for all conductors and dielectrics
        const center_points = [];
        const axis_min_for_centers = axis === 'x' ? this.x_min : 0;
        const axis_max_for_centers = axis === 'x' ? this.x_max : this.y_max;

        for (const cond of this.conductors) {
            const center = axis === 'x' ?
                (cond.x_min + cond.x_max) / 2 :
                (cond.y_min + cond.y_max) / 2;
            if (axis_min_for_centers < center && center < axis_max_for_centers) {
                center_points.push(center);
            }
        }

        for (const diel of this.dielectrics) {
            const center = axis === 'x' ?
                (diel.x_min + diel.x_max) / 2 :
                (diel.y_min + diel.y_max) / 2;
            if (axis_min_for_centers < center && center < axis_max_for_centers) {
                center_points.push(center);
            }
        }

        for (const center of center_points) {
            let min_dist = Infinity;
            for (const val of mesh) {
                min_dist = Math.min(min_dist, Math.abs(val - center));
            }
            if (min_dist > domain_size / 1000) {
                const temp = Array.from(mesh);
                temp.push(center);
                mesh = Float64Array.from(temp.sort((a, b) => a - b));
            }
        }

        // Add boundary lines adjacent to conductor edges
        const boundary_offset = Math.min(this.skin_depth * 3, domain_size / 200);
        const boundary_lines = [];
        const axis_min = axis === 'x' ? this.x_min : this.y_min;
        const axis_max = axis === 'x' ? this.x_max : this.y_max;

        for (const cond of this.conductors) {
            const cond_min = axis === 'x' ? cond.x_min : cond.y_min;
            const cond_max = axis === 'x' ? cond.x_max : cond.y_max;
            const edges = [cond_min, cond_max];

            for (const edge of edges) {
                if (Math.abs(edge - axis_min) < 1e-15 || Math.abs(edge - axis_max) < 1e-15) {
                    continue;
                }

                const is_left_edge = Math.abs(edge - cond_min) < 1e-15;
                let outside_line, inside_line;

                if (is_left_edge) {
                    outside_line = edge - boundary_offset;
                    inside_line = edge + boundary_offset;
                } else {
                    inside_line = edge - boundary_offset;
                    outside_line = edge + boundary_offset;
                }

                if (axis_min < outside_line && outside_line < axis_max) {
                    boundary_lines.push(outside_line);
                }
                if (cond_min < inside_line && inside_line < cond_max) {
                    boundary_lines.push(inside_line);
                }
            }
        }

        for (const line of boundary_lines) {
            let min_dist = Infinity;
            for (const val of mesh) {
                min_dist = Math.min(min_dist, Math.abs(val - line));
            }
            if (min_dist > boundary_offset / 3) {
                const temp = Array.from(mesh);
                temp.push(line);
                mesh = Float64Array.from(temp.sort((a, b) => a - b));
            }
        }

        // Remove duplicate or near-duplicate points first
        const mesh_unique = [mesh[0]];
        const min_spacing = domain_size * 1e-10;

        for (let i = 1; i < mesh.length; i++) {
            if (mesh[i] - mesh_unique[mesh_unique.length - 1] > min_spacing) {
                mesh_unique.push(mesh[i]);
            }
        }

        mesh = Float64Array.from(mesh_unique);

        // Check symmetry and enforce if needed (do this AFTER all additions)
        if (this.symmetric && axis === 'x') {
            const is_symmetric = this._check_symmetry(axis);
            if (is_symmetric) {
                const axis_min = axis === 'x' ? this.x_min : this.y_min;
                const axis_max = axis === 'x' ? this.x_max : this.y_max;
                mesh = this._enforce_symmetry(mesh, axis_min, axis_max);
            } else {
                console.warn('Geometry is not symmetric, but symmetric meshing was requested');
            }
        }

        return mesh;
    }

    _check_symmetry(axis) {
        if (axis !== 'x') return false;

        const center = (this.x_min + this.x_max) / 2;
        const tol = 1e-12;

        for (const cond of this.conductors) {
            const x_min_mirror = center - (cond.x_max - center);
            const x_max_mirror = center - (cond.x_min - center);

            let found_match = false;
            for (const other of this.conductors) {
                if (Math.abs(other.x_min - x_min_mirror) < tol &&
                    Math.abs(other.x_max - x_max_mirror) < tol &&
                    other.is_signal === cond.is_signal) {
                    found_match = true;
                    break;
                }
            }

            if (Math.abs(cond.x_min + cond.x_max - 2 * center) < tol) {
                found_match = true;
            }

            if (!found_match) {
                return false;
            }
        }

        return true;
    }

    _enforce_symmetry(mesh, axis_min, axis_max) {
        const center = (axis_min + axis_max) / 2;
        const domain_size = axis_max - axis_min;
        const tol = 1e-12;
        const symmetric_mesh = [];
        const used = new Array(mesh.length).fill(false);

        for (let i = 0; i < mesh.length; i++) {
            if (used[i]) continue;

            const point = mesh[i];

            if (Math.abs(point - center) < tol) {
                symmetric_mesh.push(center);
                used[i] = true;
                continue;
            }

            const mirror_pos = 2 * center - point;
            let mirror_idx = null;
            let min_dist = Infinity;

            for (let j = 0; j < mesh.length; j++) {
                if (used[j]) continue;
                const dist = Math.abs(mesh[j] - mirror_pos);
                if (dist < min_dist) {
                    min_dist = dist;
                    mirror_idx = j;
                }
            }

            if (mirror_idx !== null && min_dist < domain_size / 100) {
                const avg_pos = center + (point - center);
                const mirror_avg_pos = 2 * center - avg_pos;
                symmetric_mesh.push(avg_pos);
                symmetric_mesh.push(mirror_avg_pos);
                used[i] = true;
                used[mirror_idx] = true;
            } else {
                symmetric_mesh.push(point);
                if (axis_min < mirror_pos && mirror_pos < axis_max) {
                    symmetric_mesh.push(mirror_pos);
                }
                used[i] = true;
            }
        }

        return Float64Array.from(Array.from(new Set(symmetric_mesh)).sort((a, b) => a - b));
    }

    _validate_mesh_symmetry(mesh) {
        const center = (this.x_min + this.x_max) / 2;
        const tol = 1e-10;
        const errors = [];

        // Build a set of mesh points for fast lookup
        const mesh_set = new Set(Array.from(mesh).map(x => x.toFixed(15)));

        for (const point of mesh) {
            const mirror_pos = 2 * center - point;
            const mirror_key = mirror_pos.toFixed(15);

            // Skip the center point itself
            if (Math.abs(point - center) < tol) {
                continue;
            }

            // Check if mirror point exists
            if (!mesh_set.has(mirror_key)) {
                // Check for near matches
                let found_near = false;
                for (const other of mesh) {
                    if (Math.abs(other - mirror_pos) < tol) {
                        found_near = true;
                        break;
                    }
                }

                if (!found_near) {
                    errors.push(`Point ${point} has no symmetric counterpart (expected ${mirror_pos})`);
                }
            }
        }

        if (errors.length > 0) {
            console.error('Mesh symmetry validation failed:');
            errors.forEach(err => console.error('  ' + err));
            throw new Error(`Mesh is not symmetric: ${errors.length} symmetry violations detected`);
        }
    }

    _concat_arrays(arrays) {
        let total_len = 0;
        for (const arr of arrays) {
            total_len += arr.length;
        }
        const result = new Float64Array(total_len);
        let offset = 0;
        for (const arr of arrays) {
            result.set(arr, offset);
            offset += arr.length;
        }
        return result;
    }
}

export { Dielectric, Conductor, Mesher };
