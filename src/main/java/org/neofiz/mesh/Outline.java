package org.neofiz.mesh;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * A closed polygon, and the mesher that turns one into quadrilaterals.
 *
 * <p>Until now there were two meshes: a cylinder wall and a solid cylinder, each written out
 * by hand. That is the right amount of mesher for a validation harness, where the geometry is
 * part of the experiment and changing it would change what is being measured. It is not
 * enough for a sandbox, where the geometry is whatever the player drew.
 *
 * <h2>Cells, not a fitted mesh</h2>
 *
 * The outline is rasterised onto a square lattice: every cell whose centre lies inside becomes
 * one element, shared corners are welded, and the boundary comes out stair-stepped. A fitted
 * mesher would follow the outline exactly and is the obvious thing to want, so it is worth
 * being clear that the blocky answer is the deliberate one.
 *
 * <p>The reason is the timestep. This is an explicit solve, so the stable step is set by the
 * <em>shortest edge anywhere on the mesh</em> -- one sliver in ten thousand elements slows
 * every other element down by the same factor. A fitted mesher trims cells against the
 * outline and produces exactly those slivers, at every corner, unpredictably, as a function of
 * where the player happened to draw. A lattice cannot: every edge is the cell size, so
 * {@link QuadMesh#minimumEdgeLength()} is known before the mesh is built and the cost of a
 * shape is proportional to its area and nothing else. Stair-steps are a resolution error that
 * halves when the cell size halves. A sliver is a performance cliff that does not.
 *
 * <h2>The lattice is anchored at the origin</h2>
 *
 * Cell {@code (i, j)} always spans {@code [i*h, (i+1)*h]} by {@code [j*h, (j+1)*h]}, whatever
 * shape is being meshed. Two shapes meshed together therefore land on the same lattice, and
 * where they touch they share nodes rather than merely coinciding. Shared nodes are a welded
 * joint -- the two bodies are one body, which is what a bracket bolted to a plate should be,
 * and is not yet a substitute for contact between bodies meant to be able to separate.
 *
 * <p>Shapes that do not touch become disconnected components of one mesh. The solver has no
 * difficulty with that; they simply do not interact until there is contact.
 *
 * <h2>Coordinates</h2>
 *
 * Named {@code (r, z)} to match {@link QuadMesh}, but under
 * {@link org.neofiz.core.Formulation#PLANE_STRAIN} they are an ordinary Cartesian
 * {@code (x, y)} and nothing is revolved. Under the axisymmetric formulation {@code r} is a
 * radius and must not be negative; that is the caller's business.
 */
public final class Outline {

    /** The outer ring, {@code r0, z0, r1, z1, ...}, closed implicitly. */
    private final double[] ring;
    /** Rings subtracted from the outer one. */
    private final List<double[]> holes;

    private Outline(double[] ring, List<double[]> holes) {
        this.ring = ring;
        this.holes = holes;
    }

    /**
     * A closed polygon from interleaved coordinates, {@code r0, z0, r1, z1, ...}.
     *
     * <p>Winding does not matter: containment is decided by the even-odd rule, which is
     * insensitive to it, and the mesher emits its own counter-clockwise cells regardless.
     */
    public static Outline of(double... rz) {
        if (rz.length < 6 || rz.length % 2 != 0) {
            throw new IllegalArgumentException(
                    "a polygon needs at least three (r, z) pairs; got " + rz.length + " numbers");
        }
        return new Outline(rz.clone(), new ArrayList<>());
    }

    /** An axis-aligned rectangle. */
    public static Outline rectangle(double r0, double z0, double r1, double z1) {
        if (r1 <= r0 || z1 <= z0) {
            throw new IllegalArgumentException("rectangle needs r1 > r0 and z1 > z0");
        }
        return of(r0, z0, r1, z0, r1, z1, r0, z1);
    }

    /** A regular polygon inscribed in a circle. Use enough sides and it is a disc. */
    public static Outline circle(double centreR, double centreZ, double radius, int sides) {
        if (radius <= 0.0) throw new IllegalArgumentException("radius must be positive");
        if (sides < 3) throw new IllegalArgumentException("need at least three sides");
        final double[] rz = new double[2 * sides];
        for (int k = 0; k < sides; k++) {
            final double a = 2.0 * Math.PI * k / sides;
            rz[2 * k] = centreR + radius * Math.cos(a);
            rz[2 * k + 1] = centreZ + radius * Math.sin(a);
        }
        return new Outline(rz, new ArrayList<>());
    }

    /** This outline moved. Rings and holes travel together. */
    public Outline translated(double dr, double dz) {
        return mapped((r, z) -> new double[] {r + dr, z + dz});
    }

    /**
     * This outline turned about a point, angle in radians, positive anticlockwise.
     *
     * <p>Worth having as a primitive rather than as trigonometry at the call site, because the
     * orientation of a body is not cosmetic. A block landing flat and the same block landing on
     * a corner are different events, and only the second one can tip over.
     */
    public Outline rotated(double radians, double aboutR, double aboutZ) {
        final double c = Math.cos(radians), s = Math.sin(radians);
        return mapped((r, z) -> {
            final double dr = r - aboutR, dz = z - aboutZ;
            return new double[] {aboutR + c * dr - s * dz, aboutZ + s * dr + c * dz};
        });
    }

    /** This outline moved so its lowest point sits at {@code z}. */
    public Outline restingOn(double z) {
        double lowest = Double.MAX_VALUE;
        for (int i = 1; i < ring.length; i += 2) lowest = Math.min(lowest, ring[i]);
        return translated(0.0, z - lowest);
    }

    /** The centre of the outer ring's bounding box, {@code {r, z}}. */
    public double[] centre() {
        final double[] b = bounds();
        return new double[] {0.5 * (b[0] + b[2]), 0.5 * (b[1] + b[3])};
    }

    /** A point-to-point map, for the rigid transforms above. */
    private interface Move {
        double[] at(double r, double z);
    }

    private Outline mapped(Move move) {
        final double[] moved = ring.clone();
        for (int i = 0; i < moved.length; i += 2) {
            final double[] p = move.at(moved[i], moved[i + 1]);
            moved[i] = p[0];
            moved[i + 1] = p[1];
        }
        final List<double[]> movedHoles = new ArrayList<>(holes.size());
        for (double[] hole : holes) {
            final double[] h = hole.clone();
            for (int i = 0; i < h.length; i += 2) {
                final double[] p = move.at(h[i], h[i + 1]);
                h[i] = p[0];
                h[i + 1] = p[1];
            }
            movedHoles.add(h);
        }
        return new Outline(moved, movedHoles);
    }

    /** This outline with a region removed from it. The hole's winding is irrelevant. */
    public Outline withHole(Outline hole) {
        final List<double[]> more = new ArrayList<>(holes);
        more.add(hole.ring.clone());
        more.addAll(hole.holes);
        return new Outline(ring, more);
    }

    /**
     * Whether a point is inside, by the even-odd rule.
     *
     * <p>Even-odd is what makes holes free: a point inside a hole has crossed the outer ring
     * once and the hole ring once, so it is outside on an even count, with no need to know
     * which ring is which or which way either of them winds.
     */
    public boolean contains(double r, double z) {
        boolean in = crosses(ring, r, z);
        for (double[] hole : holes) {
            if (crosses(hole, r, z)) in = !in;
        }
        return in;
    }

    private static boolean crosses(double[] p, double r, double z) {
        boolean in = false;
        final int n = p.length / 2;
        for (int i = 0, j = n - 1; i < n; j = i++) {
            final double ri = p[2 * i], zi = p[2 * i + 1];
            final double rj = p[2 * j], zj = p[2 * j + 1];
            if ((zi > z) != (zj > z) && r < (rj - ri) * (z - zi) / (zj - zi) + ri) {
                in = !in;
            }
        }
        return in;
    }

    /** Enclosed area, holes subtracted, square metres. Always positive. */
    public double area() {
        double a = Math.abs(shoelace(ring));
        for (double[] hole : holes) a -= Math.abs(shoelace(hole));
        return a;
    }

    private static double shoelace(double[] p) {
        double a = 0.0;
        final int n = p.length / 2;
        for (int i = 0, j = n - 1; i < n; j = i++) {
            a += p[2 * j] * p[2 * i + 1] - p[2 * i] * p[2 * j + 1];
        }
        return 0.5 * a;
    }

    /** {@code rMin, zMin, rMax, zMax} of the outer ring. */
    public double[] bounds() {
        double rlo = Double.MAX_VALUE, zlo = Double.MAX_VALUE;
        double rhi = -Double.MAX_VALUE, zhi = -Double.MAX_VALUE;
        for (int i = 0; i < ring.length; i += 2) {
            rlo = Math.min(rlo, ring[i]);
            rhi = Math.max(rhi, ring[i]);
            zlo = Math.min(zlo, ring[i + 1]);
            zhi = Math.max(zhi, ring[i + 1]);
        }
        return new double[] {rlo, zlo, rhi, zhi};
    }

    /** This one shape, meshed at the given cell size. */
    public QuadMesh mesh(double cellSize) {
        return mesh(cellSize, this);
    }

    /**
     * Meshes any number of shapes onto one lattice, producing one {@link QuadMesh}.
     *
     * <p>A cell is kept when its centre is inside <em>any</em> of the shapes, so overlapping
     * shapes become their union and are meshed once rather than twice. Which shape an element
     * came from is recoverable from its centroid, which is how a material gets painted onto it.
     *
     * @param cellSize element edge length, metres. Sets both the resolution and the timestep.
     */
    /**
     * Meshes each shape as a body of its own, on the same lattice but with its own nodes.
     *
     * <p>The difference from {@link #mesh(double, Outline...)} is the whole difference between
     * a bolted assembly and a stack of objects. Shared nodes are a weld: two shapes that touch
     * are one body and cannot come apart however hard they are hit. Separate nodes are two
     * bodies that happen to be in contact, and what holds them together is
     * {@code Contact} -- so they press on one another, slide, and separate.
     *
     * <p>Both are wanted, and neither is the default for the other. A clad plate is welded. A
     * wall of bricks is not, and meshing one with {@link #mesh(double, Outline...)} produces a
     * single monolithic wall that behaves nothing like bricks.
     *
     * <p>Because the lattice is still shared, shapes meshed this way can be placed <em>exactly
     * touching</em> -- a block whose top is at z and another whose bottom is at z give two
     * coincident rows of nodes, zero gap, and a contact that is live from the first step. That
     * matters more than it sounds: an explicit solve runs in microseconds and gravity works in
     * milliseconds, so anything that has to fall into place before the interesting part starts
     * is unaffordable. Things have to begin already stacked.
     */
    public static QuadMesh assemble(double cellSize, Outline... shapes) {
        if (shapes.length == 0) throw new IllegalArgumentException("nothing to mesh");
        final List<double[]> nr = new ArrayList<>();
        final List<double[]> nz = new ArrayList<>();
        final List<int[]> conn = new ArrayList<>();
        int nodes = 0;
        for (Outline shape : shapes) {
            final QuadMesh body = mesh(cellSize, shape);
            nr.add(body.r);
            nz.add(body.z);
            final int[] c = body.conn.clone();
            for (int k = 0; k < c.length; k++) c[k] += nodes;
            conn.add(c);
            nodes += body.nodeCount;
        }

        final double[] r = new double[nodes];
        final double[] z = new double[nodes];
        int at = 0;
        for (int b = 0; b < nr.size(); b++) {
            System.arraycopy(nr.get(b), 0, r, at, nr.get(b).length);
            System.arraycopy(nz.get(b), 0, z, at, nz.get(b).length);
            at += nr.get(b).length;
        }
        int total = 0;
        for (int[] c : conn) total += c.length;
        final int[] all = new int[total];
        at = 0;
        for (int[] c : conn) {
            System.arraycopy(c, 0, all, at, c.length);
            at += c.length;
        }
        return new QuadMesh(r, z, all, new int[0]);
    }

    public static QuadMesh mesh(double cellSize, Outline... shapes) {
        if (!(cellSize > 0.0)) {
            throw new IllegalArgumentException("cell size must be positive, was " + cellSize);
        }
        if (shapes.length == 0) throw new IllegalArgumentException("nothing to mesh");

        double rlo = Double.MAX_VALUE, zlo = Double.MAX_VALUE;
        double rhi = -Double.MAX_VALUE, zhi = -Double.MAX_VALUE;
        for (Outline s : shapes) {
            final double[] b = s.bounds();
            rlo = Math.min(rlo, b[0]);
            zlo = Math.min(zlo, b[1]);
            rhi = Math.max(rhi, b[2]);
            zhi = Math.max(zhi, b[3]);
        }

        // Lattice indices spanning the bounding box. Anchored at the origin, so the cell a
        // point falls in is a property of the point and not of what else is being meshed.
        final int i0 = (int) Math.floor(rlo / cellSize);
        final int i1 = (int) Math.ceil(rhi / cellSize);
        final int j0 = (int) Math.floor(zlo / cellSize);
        final int j1 = (int) Math.ceil(zhi / cellSize);
        if ((long) (i1 - i0) * (j1 - j0) > 20_000_000L) {
            throw new IllegalArgumentException(
                    "a cell size of " + cellSize + " m over this bounding box would need more "
                            + "than twenty million cells; the shape is large or the cell small");
        }

        final Map<Long, Integer> nodeAt = new HashMap<>();
        final List<Double> nr = new ArrayList<>();
        final List<Double> nz = new ArrayList<>();
        final List<Integer> conn = new ArrayList<>();

        for (int j = j0; j < j1; j++) {
            for (int i = i0; i < i1; i++) {
                final double cr = (i + 0.5) * cellSize;
                final double cz = (j + 0.5) * cellSize;
                boolean inside = false;
                for (Outline s : shapes) {
                    if (s.contains(cr, cz)) {
                        inside = true;
                        break;
                    }
                }
                if (!inside) continue;

                // Counter-clockwise in (r, z): right along the bottom, up, left, back down.
                conn.add(node(nodeAt, nr, nz, i, j, cellSize));
                conn.add(node(nodeAt, nr, nz, i + 1, j, cellSize));
                conn.add(node(nodeAt, nr, nz, i + 1, j + 1, cellSize));
                conn.add(node(nodeAt, nr, nz, i, j + 1, cellSize));
            }
        }

        if (conn.isEmpty()) {
            throw new IllegalArgumentException(
                    "no cell centre fell inside the outline; the shape is smaller than one "
                            + cellSize + " m cell, or the coordinates are not in metres");
        }

        final double[] r = new double[nr.size()];
        final double[] z = new double[nz.size()];
        for (int k = 0; k < r.length; k++) {
            r[k] = nr.get(k);
            z[k] = nz.get(k);
        }
        final int[] c = new int[conn.size()];
        for (int k = 0; k < c.length; k++) c[k] = conn.get(k);

        return new QuadMesh(r, z, c, new int[0]);
    }

    private static int node(Map<Long, Integer> nodeAt, List<Double> nr, List<Double> nz,
                            int i, int j, double h) {
        final long key = ((long) i << 32) ^ (j & 0xFFFFFFFFL);
        final Integer have = nodeAt.get(key);
        if (have != null) return have;
        final int index = nr.size();
        nr.add(i * h);
        nz.add(j * h);
        nodeAt.put(key, index);
        return index;
    }
}
