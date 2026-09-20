package org.neofiz.mesh;

/**
 * A mapped mesh of 4-node quadrilaterals in the (r, z) half-plane.
 *
 * <p>Deliberately structure-of-arrays over flat primitive arrays rather than an object
 * graph. This is the layout that ports to a CUDA or HLSL kernel without restructuring,
 * and it is what makes the batched-ensemble architecture affordable later: N lanes are
 * N strides into the same buffers.
 *
 * <p>Node ordering within an element is counter-clockwise in (r, z), with r playing the
 * role of x and z of y. A positive Jacobian determinant is therefore expected everywhere
 * and is checked at construction.
 */
public final class QuadMesh {

    public final int nodeCount;
    public final int elementCount;

    /** Reference radial coordinate per node, metres. */
    public final double[] r;
    /** Reference axial coordinate per node, metres. */
    public final double[] z;
    /** Element connectivity, 4 node indices per element, counter-clockwise. */
    public final int[] conn;

    /**
     * Boundary edges carrying pressure, 2 node indices per edge.
     *
     * <p>Ordering convention: traversing the edge from the first node to the second,
     * rotating the tangent by -90 degrees yields the direction a positive pressure
     * pushes the solid. For a bore surface at constant r, traversed in increasing z,
     * that direction is +r, which is outward from the fluid and into the wall.
     */
    public final int[] pressureEdges;

    public QuadMesh(double[] r, double[] z, int[] conn, int[] pressureEdges) {
        if (r.length != z.length) {
            throw new IllegalArgumentException("r and z must have the same length");
        }
        if (conn.length % 4 != 0) {
            throw new IllegalArgumentException("connectivity must be a multiple of 4");
        }
        if (pressureEdges.length % 2 != 0) {
            throw new IllegalArgumentException("pressure edges must be a multiple of 2");
        }
        this.r = r;
        this.z = z;
        this.conn = conn;
        this.pressureEdges = pressureEdges;
        this.nodeCount = r.length;
        this.elementCount = conn.length / 4;
        checkOrientation();
    }

    /**
     * Shortest edge over the whole mesh. This is the characteristic length that sets the
     * CFL timestep. Using the minimum edge rather than area/diagonal is the conservative
     * choice and stays correct on stretched meshes.
     */
    public double minimumEdgeLength() {
        double min = Double.MAX_VALUE;
        for (int e = 0; e < elementCount; e++) {
            int b = e * 4;
            for (int k = 0; k < 4; k++) {
                int i = conn[b + k];
                int j = conn[b + (k + 1) % 4];
                double dr = r[j] - r[i];
                double dz = z[j] - z[i];
                double len = Math.sqrt(dr * dr + dz * dz);
                if (len < min) min = len;
            }
        }
        return min;
    }

    /** Signed area of an element in the (r, z) plane; positive for correct ordering. */
    public double signedArea(int element) {
        int b = element * 4;
        double a = 0.0;
        for (int k = 0; k < 4; k++) {
            int i = conn[b + k];
            int j = conn[b + (k + 1) % 4];
            a += r[i] * z[j] - r[j] * z[i];
        }
        return 0.5 * a;
    }

    /**
     * Radial coordinate of an element's area centroid, in the reference configuration.
     *
     * <p>The true polygon centroid, not the mean of the four corners. On the mapped meshes
     * here those agree, because every element is a rectangle; on a graded or distorted mesh
     * they do not, and the difference goes straight into the ring volume below.
     */
    public double centroidRadius(int element) {
        return centroid(element, r);
    }

    /** Axial coordinate of an element's area centroid, in the reference configuration. */
    public double centroidZ(int element) {
        return centroid(element, z);
    }

    /** Polygon area centroid along one axis, {@code sum (c_i + c_j) cross / (3 sum cross)}. */
    private double centroid(int element, double[] coordinate) {
        int b = element * 4;
        double a = 0.0;
        double moment = 0.0;
        for (int k = 0; k < 4; k++) {
            int i = conn[b + k];
            int j = conn[b + (k + 1) % 4];
            double cross = r[i] * z[j] - r[j] * z[i];
            a += cross;
            moment += (coordinate[i] + coordinate[j]) * cross;
        }
        return moment / (3.0 * a);
    }

    /**
     * Volume of the ring an element sweeps when revolved about the axis, cubic metres.
     *
     * <p>Pappus's theorem: {@code 2 pi rc A}. This is the volume a defect field has to
     * normalise by, and it is the reason an axisymmetric element cannot be treated as a small
     * piece of material with one flaw in it. An element at r = 60 mm with a 1 mm square
     * section is a torus holding 377 cubic millimetres -- not one flaw, but the weakest of
     * however many the circumference contains.
     *
     * <p>Exactly additive, so the sum over a mesh is the analytic volume of the body to
     * round-off, which is the cheapest available check that a mesh is what it claims to be.
     */
    public double ringVolume(int element) {
        return 2.0 * Math.PI * centroidRadius(element) * signedArea(element);
    }

    /** Total revolved volume of the mesh, cubic metres. */
    public double totalRingVolume() {
        double v = 0.0;
        for (int e = 0; e < elementCount; e++) v += ringVolume(e);
        return v;
    }

    /**
     * The free surface: every edge that belongs to exactly one element.
     *
     * @param edges 2 node indices per edge, in the winding of the element that owns it, so
     *              that rotating the tangent {@code (dr, dz)} by -90 degrees to
     *              {@code (dz, -dr)} points <em>out</em> of the solid
     * @param owner the element each edge belongs to, same order
     * @param nodes every node appearing on the surface, ascending and without repeats
     */
    public record Surface(int[] edges, int[] owner, int[] nodes) {

        public int edgeCount() {
            return owner.length;
        }
    }

    /**
     * Finds the free surface from connectivity alone.
     *
     * <p>Topological rather than geometric, which is what makes it work on a mesh nobody
     * designed: an interior edge is shared by two elements and a surface edge is not, and that
     * is true for a hole, a crack, a body of several disconnected pieces, or whatever a
     * polygon mesher produced. A traversal that walked the outline instead would need to know
     * in advance how many outlines there were.
     *
     * <p>Emitted in element order rather than in map order, so two runs on the same mesh
     * produce the same surface in the same sequence. Contact forces are accumulated in this
     * order and floating-point addition is not associative, so the ordering is part of the
     * solver being reproducible rather than a tidiness preference.
     */
    public Surface surface() {
        return surface(null);
    }

    /**
     * The free surface of what is left when some elements are ignored.
     *
     * <p>The reason this takes an argument at all is erosion. Deleting an element exposes the
     * faces it was hiding: a body eroded through the middle has two new free surfaces that
     * were interior a moment earlier, and anything that kept using the original surface would
     * let the two halves pass straight through one another. Recomputing with the dead elements
     * skipped is the whole of the update.
     *
     * @param skip one flag per element, or null to include them all
     */
    public Surface surface(boolean[] skip) {
        final java.util.HashMap<Long, Integer> count = new java.util.HashMap<>();
        for (int e = 0; e < elementCount; e++) {
            if (skip != null && skip[e]) continue;
            for (int k = 0; k < 4; k++) {
                count.merge(edgeKey(conn[e * 4 + k], conn[e * 4 + (k + 1) % 4]), 1, Integer::sum);
            }
        }

        int free = 0;
        for (int v : count.values()) if (v == 1) free++;

        final int[] edges = new int[2 * free];
        final int[] owner = new int[free];
        final boolean[] onSurface = new boolean[nodeCount];
        int at = 0;
        for (int e = 0; e < elementCount; e++) {
            if (skip != null && skip[e]) continue;
            for (int k = 0; k < 4; k++) {
                final int a = conn[e * 4 + k];
                final int b = conn[e * 4 + (k + 1) % 4];
                if (count.get(edgeKey(a, b)) != 1) continue;
                edges[2 * at] = a;
                edges[2 * at + 1] = b;
                owner[at] = e;
                at++;
                onSurface[a] = true;
                onSurface[b] = true;
            }
        }

        int n = 0;
        for (boolean b : onSurface) if (b) n++;
        final int[] nodes = new int[n];
        int m = 0;
        for (int i = 0; i < nodeCount; i++) if (onSurface[i]) nodes[m++] = i;

        return new Surface(edges, owner, nodes);
    }

    private static long edgeKey(int a, int b) {
        return (long) Math.min(a, b) << 32 | Math.max(a, b);
    }

    private void checkOrientation() {
        for (int e = 0; e < elementCount; e++) {
            if (signedArea(e) <= 0.0) {
                throw new IllegalStateException(
                        "element " + e + " has non-positive area; node ordering must be "
                                + "counter-clockwise in (r, z)");
            }
        }
    }

    /**
     * A mapped annular mesh: the wall of a cylinder, seen as a rectangle in (r, z).
     *
     * <p>This is the M0 geometry and the charter's worked example. The inner face carries
     * pressure; the caller decides what to do with the axial ends.
     *
     * @param innerRadius bore radius, metres
     * @param outerRadius outside radius, metres
     * @param height      axial extent of the modelled slice, metres
     * @param nr          elements through the wall
     * @param nz          elements along the axis
     */
    public static QuadMesh cylinderWall(double innerRadius, double outerRadius, double height,
                                        int nr, int nz) {
        if (outerRadius <= innerRadius) {
            throw new IllegalArgumentException("outer radius must exceed inner radius");
        }
        if (nr < 1 || nz < 1) {
            throw new IllegalArgumentException("need at least one element in each direction");
        }

        int nodesR = nr + 1;
        int nodesZ = nz + 1;
        double[] r = new double[nodesR * nodesZ];
        double[] z = new double[nodesR * nodesZ];

        for (int j = 0; j < nodesZ; j++) {
            for (int i = 0; i < nodesR; i++) {
                int n = j * nodesR + i;
                r[n] = innerRadius + (outerRadius - innerRadius) * i / nr;
                z[n] = height * j / nz;
            }
        }

        int[] conn = new int[nr * nz * 4];
        int c = 0;
        for (int j = 0; j < nz; j++) {
            for (int i = 0; i < nr; i++) {
                conn[c++] = j * nodesR + i;
                conn[c++] = j * nodesR + i + 1;
                conn[c++] = (j + 1) * nodesR + i + 1;
                conn[c++] = (j + 1) * nodesR + i;
            }
        }

        // Bore face: the i = 0 column, traversed in increasing z so that the -90 degree
        // rotation of the tangent points in +r, into the wall.
        int[] edges = new int[nz * 2];
        for (int j = 0; j < nz; j++) {
            edges[j * 2] = j * nodesR;
            edges[j * 2 + 1] = (j + 1) * nodesR;
        }

        return new QuadMesh(r, z, conn, edges);
    }

    /**
     * A mapped mesh of a solid cylinder: the full radius, seen as a rectangle in (r, z).
     *
     * <p>This is the Taylor-impact geometry. It differs from {@link #cylinderWall} in the one
     * way that matters -- there is a column of nodes <em>on the axis</em>, at r = 0 -- and it
     * carries no pressure edges, because the r = 0 column is a symmetry line rather than a
     * surface and pressure on it is meaningless.
     *
     * <p><b>Why r = 0 is safe.</b> The axisymmetric formulation divides by radius to form the
     * hoop strain, eps_theta = u_r / r, which looks fatal on the axis. It is not, because that
     * quotient is only ever evaluated at a <em>quadrature</em> point -- the centroid under the
     * reduced rule, the four Gauss points under the full one -- and every quadrature point of
     * an element with non-zero radial extent lies strictly inside it, at r &gt; 0. The axis
     * nodes contribute through shape functions and never appear in a denominator.
     *
     * <p>What the axis <em>does</em> require is that its radial degree of freedom be pinned:
     * material on the centreline cannot move off it without tearing the body open, so
     * u_r(0, z) = 0 is a symmetry condition, not a boundary condition the caller may choose.
     * {@link #axisNodes} exists to make that easy to get right, and forgetting it does not
     * crash -- it quietly lets the axis open into a hole.
     *
     * @param radius outside radius, metres
     * @param height axial length, metres
     * @param nr     elements along the radius
     * @param nz     elements along the axis
     */
    public static QuadMesh solidCylinder(double radius, double height, int nr, int nz) {
        if (radius <= 0.0) throw new IllegalArgumentException("radius must be positive");
        if (height <= 0.0) throw new IllegalArgumentException("height must be positive");
        if (nr < 1 || nz < 1) {
            throw new IllegalArgumentException("need at least one element in each direction");
        }

        int nodesR = nr + 1;
        int nodesZ = nz + 1;
        double[] r = new double[nodesR * nodesZ];
        double[] z = new double[nodesR * nodesZ];

        for (int j = 0; j < nodesZ; j++) {
            for (int i = 0; i < nodesR; i++) {
                int n = j * nodesR + i;
                r[n] = radius * i / nr;
                z[n] = height * j / nz;
            }
        }

        int[] conn = new int[nr * nz * 4];
        int c = 0;
        for (int j = 0; j < nz; j++) {
            for (int i = 0; i < nr; i++) {
                conn[c++] = j * nodesR + i;
                conn[c++] = j * nodesR + i + 1;
                conn[c++] = (j + 1) * nodesR + i + 1;
                conn[c++] = (j + 1) * nodesR + i;
            }
        }

        return new QuadMesh(r, z, conn, new int[0]);
    }

    /**
     * Node indices on the symmetry axis, r = 0, ordered by increasing z. Every one of these
     * needs its radial degree of freedom pinned. See {@link #solidCylinder}.
     */
    public static int[] axisNodes(int nr, int nz) {
        return innerSurfaceNodes(nr, nz);
    }

    /** Node indices on the z = 0 face, ordered by increasing r. */
    public static int[] nearFaceNodes(int nr, int nz) {
        int[] out = new int[nr + 1];
        for (int i = 0; i <= nr; i++) out[i] = i;
        return out;
    }

    /** Node indices on the z = height face, ordered by increasing r. */
    public static int[] farFaceNodes(int nr, int nz) {
        int[] out = new int[nr + 1];
        for (int i = 0; i <= nr; i++) out[i] = nz * (nr + 1) + i;
        return out;
    }

    /** Node indices on the inner (bore) surface, ordered by increasing z. */
    public static int[] innerSurfaceNodes(int nr, int nz) {
        int[] out = new int[nz + 1];
        for (int j = 0; j <= nz; j++) out[j] = j * (nr + 1);
        return out;
    }

    /** Node indices on the outer surface, ordered by increasing z. */
    public static int[] outerSurfaceNodes(int nr, int nz) {
        int[] out = new int[nz + 1];
        for (int j = 0; j <= nz; j++) out[j] = j * (nr + 1) + nr;
        return out;
    }

    /** Node indices on the z = 0 and z = height faces. */
    public static int[] axialEndNodes(int nr, int nz) {
        int[] out = new int[2 * (nr + 1)];
        for (int i = 0; i <= nr; i++) {
            out[i] = i;
            out[nr + 1 + i] = nz * (nr + 1) + i;
        }
        return out;
    }
}
