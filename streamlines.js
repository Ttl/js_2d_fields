function isInsideConductor(x, y, conductors) {
    for (const c of conductors) {
        if (
            x >= c.x_min &&
            x <= c.x_max &&
            y >= c.y_min &&
            y <= c.y_max
        ) {
            return true;
        }
    }
    return false;
}

function findConductorIntersection(x1, y1, x2, y2, conductors) {
    // Check if line segment from (x1,y1) to (x2,y2) crosses into a conductor
    // Returns the closest intersection point if multiple exist
    
    let closestHit = null;
    let closestT = Infinity;
    
    for (const c of conductors) {
        // Check top edge (entering from above)
        if (y1 >= c.y_max && y2 < c.y_max) {
            const t = (c.y_max - y1) / (y2 - y1);
            const xHit = x1 + t * (x2 - x1);
            if (xHit >= c.x_min && xHit <= c.x_max && t < closestT) {
                closestT = t;
                closestHit = { x: xHit, y: c.y_max, conductor: c };
            }
        }
        
        // Check bottom edge (entering from below)
        if (y1 <= c.y_min && y2 > c.y_min) {
            const t = (c.y_min - y1) / (y2 - y1);
            const xHit = x1 + t * (x2 - x1);
            if (xHit >= c.x_min && xHit <= c.x_max && t < closestT) {
                closestT = t;
                closestHit = { x: xHit, y: c.y_min, conductor: c };
            }
        }
        
        // Check left edge (entering from the left)
        if (x1 <= c.x_min && x2 > c.x_min) {
            const t = (c.x_min - x1) / (x2 - x1);
            const yHit = y1 + t * (y2 - y1);
            if (yHit >= c.y_min && yHit <= c.y_max && t < closestT) {
                closestT = t;
                closestHit = { x: c.x_min, y: yHit, conductor: c };
            }
        }
        
        // Check right edge (entering from the right)
        if (x1 >= c.x_max && x2 < c.x_max) {
            const t = (c.x_max - x1) / (x2 - x1);
            const yHit = y1 + t * (y2 - y1);
            if (yHit >= c.y_min && yHit <= c.y_max && t < closestT) {
                closestT = t;
                closestHit = { x: c.x_max, y: yHit, conductor: c };
            }
        }
    }
    
    return closestHit;
}

function distToConductor(x, y, conductors) {
    let minDist = Infinity;
    for (const c of conductors) {
        // Distance to nearest edge
        const dx = Math.max(c.x_min - x, 0, x - c.x_max);
        const dy = Math.max(c.y_min - y, 0, y - c.y_max);
        minDist = Math.min(minDist, Math.hypot(dx, dy));
    }
    return minDist;
}

function generateConductorSeedsWeighted(
    conductors, spacing,
    xArr, yArr, Ex, Ey,
    yOffset = 1e-7
) {
    const seeds = [];
    const eps = 1e-12;

    for (const c of conductors) {
        if (!c.is_signal) continue;
        if (c.polarity < 0) continue;
        const width = c.x_max - c.x_min;
        const nSamples = Math.max(50, Math.floor(width / spacing) * 4);

        // Process both top and bottom edges
        const edges = [
            { y: c.y_max + yOffset, normalSign: 1 },   // Top edge: +y normal
            { y: c.y_min - yOffset, normalSign: -1 }   // Bottom edge: -y normal
        ];

        for (const edge of edges) {
            // Sample points along edge
            const xs = [];
            const w = [];

            for (let i = 0; i < nSamples; i++) {
                const x = c.x_min + (i / (nSamples - 1)) * width;
                xs.push(x);
                const f = sampleField(x, edge.y, xArr, yArr, Ex, Ey, conductors);
                if (!f) {
                    w.push(eps);
                } else {
                    // Normal field component (use abs of Ey for weighting)
                    w.push(Math.abs(f.ey) + eps);
                }
            }

            // Build CDF
            const C = [];
            let sum = 0;
            for (let i = 0; i < w.length; i++) {
                sum += w[i];
                C.push(sum);
            }

            // Number of streamlines to launch (similar to non-weighted version)
            const nSeeds = Math.max(10, Math.floor(width / spacing));

            // Uniform sampling in CDF space
            for (let k = 0; k <= nSeeds; k++) {
                const target = (k / nSeeds) * sum;
                let i = C.findIndex(v => v >= target);
                if (i < 0) i = C.length - 1;
                seeds.push([xs[i], edge.y]);
            }
        }
    }

    return seeds;
}

function makeStreamlineTraceFromConductors(
    Ex, Ey,
    xSolver, ySolver,
    conductors
) {
    const xLines = [];
    const yLines = [];

    const spacing = 0.5 * (xSolver[xSolver.length - 1] - xSolver[0]) / xSolver.length;
    const ds = spacing / 4;
    const maxSteps = 800;

    const seeds = generateConductorSeedsWeighted(
        conductors, spacing,
        xSolver, ySolver,
        Ex, Ey,
    );

    for (const [x0, y0] of seeds) {
        const line = traceStreamline(
            x0, y0,
            xSolver, ySolver,
            Ex, Ey,
            ds, maxSteps,
            conductors
        );

        if (line.length < 2) continue;

        xLines.push(line[0][0] * 1000);
        yLines.push(line[0][1] * 1000);

        for (let k = 1; k < line.length; k++) {
            xLines.push(line[k][0] * 1000);
            yLines.push(line[k][1] * 1000);
        }

        xLines.push(null);
        yLines.push(null);
    }

    return {
        type: "scatter",
        mode: "lines",
        x: xLines,
        y: yLines,
        line: { width: 1, color: "black" },
        hoverinfo: "skip",
        name: "E-field lines"
    };
}

function sampleV(x, y, xArr, yArr, V) {
    // Find surrounding indices
    let i = yArr.findIndex(v => v > y) - 1;
    let j = xArr.findIndex(v => v > x) - 1;

    if (i < 0 || j < 0 || i >= yArr.length - 1 || j >= xArr.length - 1) {
        return null;
    }

    const x1 = xArr[j], x2 = xArr[j + 1];
    const y1 = yArr[i], y2 = yArr[i + 1];

    const tx = (x - x1) / (x2 - x1);
    const ty = (y - y1) / (y2 - y1);

    function lerp(a, b, t) {
        return a * (1 - t) + b * t;
    }

    const v =
        lerp(
            lerp(V[i][j], V[i][j + 1], tx),
            lerp(V[i + 1][j], V[i + 1][j + 1], tx),
            ty
        );
    return v;
}

function sampleField(x, y, xArr, yArr, Ex, Ey, conductors) {
    if (isInsideConductor(x, y, conductors)) return null;

    // Find surrounding indices
    let i = yArr.findIndex(v => v > y) - 1;
    let j = xArr.findIndex(v => v > x) - 1;

    if (i < 0 || j < 0 || i >= yArr.length - 1 || j >= xArr.length - 1) {
        return null;
    }

    const x1 = xArr[j], x2 = xArr[j + 1];
    const y1 = yArr[i], y2 = yArr[i + 1];

    const tx = (x - x1) / (x2 - x1);
    const ty = (y - y1) / (y2 - y1);

    function lerp(a, b, t) {
        return a * (1 - t) + b * t;
    }

    const ex =
        lerp(
            lerp(Ex[i][j], Ex[i][j + 1], tx),
            lerp(Ex[i + 1][j], Ex[i + 1][j + 1], tx),
            ty
        );

    const ey =
        lerp(
            lerp(Ey[i][j], Ey[i][j + 1], tx),
            lerp(Ey[i + 1][j], Ey[i + 1][j + 1], tx),
            ty
        );

    return { ex, ey };
}

function snapToNormal(x, y, conductors, f) {
    for (const c of conductors) {
        const tol = 1e-9;

        if (Math.abs(y - c.y_max) < tol &&
            x >= c.x_min && x <= c.x_max) {
            return { ex: 0, ey: Math.sign(f.ey) };
        }
    }
    return f;
}

function backtrackToConductor(x0, y0, xArr, yArr, Ex, Ey, ds, conductors) {
    // Trace backwards (against field direction) to find conductor surface
    const maxSteps = 100;
    const points = [[x0, y0]];
    let x = x0;
    let y = y0;
    
    for (let n = 0; n < maxSteps; n++) {
        const f = sampleField(x, y, xArr, yArr, Ex, Ey, conductors);
        if (!f) break;
        
        const m = Math.hypot(f.ex, f.ey);
        if (m === 0) break;

        const dist = distToConductor(x, y, conductors);
        const dsLocal = Math.min(ds, dist / 2); // Smaller steps near conductors
        
        // Step BACKWARDS (negative direction)
        const xNew = x - dsLocal * f.ex / m;
        const yNew = y - dsLocal * f.ey / m;
        
        // Check if we hit a conductor
        const hit = findConductorIntersection(x, y, xNew, yNew, conductors);
        if (hit) {
            points.unshift([hit.x, hit.y]);
            return points;
        }
        
        // Check if inside conductor (fallback)
        if (isInsideConductor(xNew, yNew, conductors)) {
            // Use current point as best approximation
            return points;
        }
        
        points.unshift([xNew, yNew]);
        x = xNew;
        y = yNew;
    }
    
    return points;
}

function traceStreamline(x0, y0, xArr, yArr, Ex, Ey, ds, maxSteps, conductors) {
    // First, backtrack to conductor surface
    const backtrack = backtrackToConductor(x0, y0, xArr, yArr, Ex, Ey, ds / 2, conductors);
    console.log(backtrack.length);
    
    // Start with backtracked points
    const line = backtrack;

    let x = x0;
    let y = y0;

    for (let n = 0; n < maxSteps; n++) {
        if (isInsideConductor(x, y, conductors)) break;

        const f1_init = sampleField(x, y, xArr, yArr, Ex, Ey, conductors);
        const f1 = snapToNormal(x, y, conductors, f1_init);
        if (!f1) break;

        const m1 = Math.hypot(f1.ex, f1.ey);
        if (m1 === 0) break;

        const k1x = f1.ex / m1;
        const k1y = f1.ey / m1;

        const dist = distToConductor(x, y, conductors);
        const dsLocal = Math.min(ds, Math.max(dist / 2, ds / 10)); // Smaller steps near conductors

        const f2 = sampleField(x + 0.5 * ds * k1x, y + 0.5 * dsLocal * k1y, xArr, yArr, Ex, Ey, conductors);
        if (!f2) break;
        const m2 = Math.hypot(f2.ex, f2.ey);
        if (m2 === 0) break;

        const k2x = f2.ex / m2;
        const k2y = f2.ey / m2;

        const f3 = sampleField(x + 0.5 * ds * k2x, y + 0.5 * dsLocal * k2y, xArr, yArr, Ex, Ey, conductors);
        if (!f3) break;
        const m3 = Math.hypot(f3.ex, f3.ey);
        if (m3 === 0) break;
        const k3x = f3.ex / m3;
        const k3y = f3.ey / m3;

        const f4 = sampleField(x + ds * k3x, y + dsLocal * k3y, xArr, yArr, Ex, Ey, conductors);
        if (!f4) break;
        const m4 = Math.hypot(f4.ex, f4.ey);
        if (m4 === 0) break;
        const k4x = f4.ex / m4;
        const k4y = f4.ey / m4;

        const xNew = x + ds * (k1x + 2*k2x + 2*k3x + k4x) / 6;
        const yNew = y + ds * (k1y + 2*k2y + 2*k3y + k4y) / 6;

        // Check if we crossed into a conductor
        const hit = findConductorIntersection(x, y, xNew, yNew, conductors);
        if (hit) {
            line.push([hit.x, hit.y]); // Terminate exactly at surface
            break;
        }

        x = xNew;
        y = yNew;
        line.push([x, y]);
    }

    return line;
}

export { makeStreamlineTraceFromConductors };
