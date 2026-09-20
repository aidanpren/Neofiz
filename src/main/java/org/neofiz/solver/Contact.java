package org.neofiz.solver;

import org.neofiz.mesh.QuadMesh;

/**
 * Penalty contact between deformable surfaces, with a Coulomb cone.
 *
 * <p>{@link RigidWall} is the other contact model here and it is not superseded by this one.
 * A plane that is rigid and cannot move is best handled kinematically -- non-penetration is
 * exact, no stiffness has to be chosen, and no frequency is added to the system. None of that
 * survives a second deformable body: there is no prescribed position to snap a node to, both
 * sides accelerate, and the constraint has to be solved for rather than imposed. So this is
 * the penalty method, and the parts of it that are usually free parameters are pinned down
 * below rather than left for a caller to guess at.
 *
 * <h2>Node against segment</h2>
 *
 * The free surface is a closed polyline of two-node segments ({@link QuadMesh#surface()}).
 * Every surface node is checked against every nearby segment; a node that has crossed a
 * segment is pushed back along that segment's outward normal, and the segment's two nodes take
 * the reaction, split by the shape functions at the point of contact. Momentum is therefore
 * conserved pair by pair rather than approximately -- what one side gains the other loses, in
 * the same step, along the same line.
 *
 * <p>Deciding what counts as a contact takes three tests, and only the first is obvious.
 * <b>Projection</b>: the node has to project inside the segment rather than past its end,
 * because clamping to the ends invents a contact in the quarter-plane diagonally beyond a
 * convex corner, where the node is outside the body altogether. <b>Facing</b>: the outward
 * normal at the node, averaged from the segments meeting there, has to oppose the segment's.
 * <b>Depth</b>: the node has to be behind the segment by less than half its length.
 *
 * <p>The second and third are not belt and braces. The region behind a segment is a half-plane
 * and not a box, so a node on the bottom face of a block projects squarely inside the slab of
 * its own left face and reads as two millimetres deep behind it. On an undeformed lattice the
 * depth cap hides that, because two faces of a body are a whole cell apart and the cap is half
 * of one -- but a deforming body brings them inside it, and since the master is chosen by
 * greatest depth, the invented pair would beat the real one and drive the node sideways with
 * several times the right force. The facing test is what excludes it, and the depth cap is
 * then left doing the job it is named for.
 *
 * <p>Two further rings of the surface are excluded topologically. A closed surface is always
 * exactly touching itself, and one ring is not enough at a sharp corner, where the segment on
 * the far side of a neighbour is close enough to register while being the same material.
 *
 * <p>A node in reach of several segments takes the deepest one only, ties broken by the lower
 * segment index. One master per node keeps the force from being counted twice in the overlap
 * behind a concave corner, and the tie-break is what makes that choice reproducible.
 *
 * <h2>The stiffness is not a free parameter</h2>
 *
 * The usual penalty stiffness is built from a bulk modulus and an element size, which makes it
 * a material property -- and then a steel node meeting a foam segment has two answers and the
 * pair needs a rule for combining them, while the added frequency can quietly exceed what the
 * timestep can integrate. This uses the mass-based form instead:
 *
 * <pre>  k = SCALE * m* / dt0^2       m* = reduced mass of the node and the segment point</pre>
 *
 * whose frequency is {@code omega = sqrt(SCALE)/dt0}, independent of the materials, the mesh
 * and the step. Central difference is stable while {@code omega*dt < 2}, and a pair is found
 * from <em>both</em> sides -- every surface here is slave and master -- so the effective
 * stiffness is twice this and the bound is {@code SCALE < 2}. Inside it the contact can never
 * be the thing that breaks the CFL limit. The cost is penetration: a node arriving at speed v
 * turns around in a quarter period and reaches a depth of about {@code v*dt/sqrt(SCALE)}.
 *
 * <p>{@code dt0} is the step in the <em>reference</em> configuration, not the current one. The
 * two part company as soon as an element compresses, and a spring whose constant rises while it
 * is loaded is not a conservative element -- it manufactures energy at the rate the stiffness
 * climbs. Stability survives the change because the current step only falls below the reference
 * one as the mesh crushes, and a falling step makes the contact frequency safer rather than
 * riskier.
 *
 * <p>Being stable is not the same as being usable. Measured on the head-on collision the gate
 * runs -- two elastic blocks at 60 m/s, meshed at 1 mm:
 *
 * <pre>
 *   SCALE    energy closes to    deepest penetration
 *   0.01         -3.2 %              25.5 % of a cell
 *   0.02         +0.4 %              16.6 %
 *   0.05         +2.8 %               8.7 %
 *   0.10         +2.2 %               4.6 %
 *   0.20         +6.0 %               3.5 %
 *   0.40         +3.2 %               1.8 %
 *   0.80         +4.1 %               1.6 %
 *   1.50       +553   %               0.4 %
 *   3.00         diverges
 * </pre>
 *
 * At the soft end the surfaces sink far enough into one another that nodes slide out of the
 * projection windows they were resolving; the springs holding them vanish without doing the
 * work of releasing, and energy is <em>lost</em>. Through the middle two decades the closure
 * is a few per cent with no clean trend -- it is not converging to anything, it is the size of
 * the discrete release impulse and it wanders. Past 1 the double-sided stiffness overruns what
 * the step can integrate and the run leaves the rails. The default sits where penetration is
 * already under five per cent of a cell and there is still an order of magnitude of headroom.
 *
 * <p>A <b>resting</b> contact is far better behaved than that table suggests: a block dropped
 * on a block and left to slide across it closes to 0.15 % over six thousand steps. This is
 * worth saying because it was not true until the facing test above existed, and the way it
 * failed is instructive -- a single spurious pair, admitted because two perpendicular surfaces
 * have a dot product of exactly zero, was worth more energy than the entire impact. The
 * symptom looked exactly like friction pumping energy in, and a normal damper was added to
 * suppress it before the cause was found. The damper is still here ({@link #withDamping}) and
 * defaults to nothing, because with the real bug fixed it makes every measurement above
 * slightly worse.
 *
 * <p>One thing an undamped penalty contact does not do is <em>settle</em>. A stack of blocks
 * under gravity finds its equilibrium at the penetration where the spring carries the weight
 * -- for steel at these meshes, about two picometres -- and then oscillates about it for ever,
 * because nothing removes the energy. It is harmless and invisible, being twelve orders of
 * magnitude below the element size, but it does mean {@link #pairCount()} flickers on a scene
 * that looks perfectly still. A little {@link #withDamping} is the cure where that matters.
 *
 * <p>Momentum, unlike energy, is exact to round-off, and that is what the gate leans on.
 *
 * <p>Because the stiffness is built from nodal mass, this is also formulation-agnostic. Under
 * the axisymmetric rule a node's mass is already the mass of the ring it represents and every
 * assembled force is already a ring force, so k comes out as a ring stiffness with no special
 * case anywhere.
 *
 * <h2>Friction, as a return map</h2>
 *
 * Same shape as {@link J2} and as {@link RigidWall}: compute the tangential force that would
 * remove the relative sliding velocity this step, and if it lies outside the cone
 * {@code |f_t| <= mu*f_n}, return it to the cone surface. Sticking is then exact rather than
 * approached through a tangential spring, and there is no stored slip to maintain for pairs
 * that appear and vanish as surfaces move across one another.
 *
 * <h2>Serial, deliberately</h2>
 *
 * Contact pairs are discovered per step and accumulate into shared nodal force arrays and into
 * running energy totals, so a threaded version would make the answer depend on the schedule.
 * The rest of this solver is bit-for-bit identical at any thread count and that property is
 * worth more than the few per cent this costs -- the work here is over the surface, which is
 * O(sqrt(n)) of a mesh, while the element assembly that stays threaded is O(n).
 */
public final class Contact {

    /**
     * Penalty stiffness as a fraction of {@code m/dt^2}. See the class note: the added
     * frequency is {@code sqrt(SCALE)/dt} regardless of anything else, so this is stable for
     * any value below 4 and the choice is purely a penetration-versus-smoothness trade.
     */
    public static final double DEFAULT_SCALE = 0.1;

    /**
     * A node deeper than this many times the surrounding segment length is ignored.
     *
     * <p>Not a tolerance -- a guard against the worst failure penalty contact has. A node that
     * has passed entirely through a thin body is geometrically "behind" the far surface by a
     * large amount, and pushing it back out the way it came would fire it off at a speed that
     * has nothing to do with the collision. Dropping the pair instead lets it leave, which is
     * wrong but bounded, and {@link #escapedNodes()} counts it so the run says so.
     */
    private static final double MAX_DEPTH = 0.5;

    /**
     * Fraction of the relative sliding a sticking pair removes in one step.
     *
     * <p>Not a fudge factor. A contact between two deformable surfaces is found <b>twice</b>,
     * once from each side, because both surfaces are made of slave nodes and master segments.
     * For the normal spring that is harmless -- two springs in parallel are a stiffer spring.
     * For a tangential force computed to remove <em>all</em> the relative sliding it is not: 
     * applied twice it removes twice as much, which reverses the slip rather than stopping it,
     * and applied every step that is an energy source. It took a block sliding on a block from
     * closing its budget to gaining 245 %.
     *
     * <p>Half, so that the double application lands on exactly full arrest and anything less
     * than double decays geometrically towards it. A sticking pair is still stuck within a few
     * steps and every step of it is dissipative.
     */
    private static final double STICK_RELAXATION = 0.5;
    /**
     * How squarely two surfaces have to face one another to be touching.
     *
     * <p>A cosine, so this admits everything up to about 104 degrees of opposition: two flat
     * faces meet at 180 and read -1, a right-angled corner indenting a flat face reads -0.707,
     * and two perpendicular surfaces read 0.
     *
     * <p>Zero would be the natural threshold and it is not safe, because the perpendicular
     * case sits exactly on it. A block landing on another while sliding sideways puts a
     * top-face node of the lower body (normal straight up) against the vertical side face of
     * the upper one (normal straight out) -- a dot product of exactly zero, no contact at all
     * -- and a hair of elastic deformation tilts the node normal enough to make it very
     * slightly negative. The pair is then admitted with a depth equal to the whole sideways
     * overlap rather than to any real penetration, which on a 1 mm mesh was one spurious pair
     * worth 2.4 J against a 2.5 J impact. It read as though friction were pumping energy in.
     *
     * <p>Rejecting the grazing pairing costs nothing, because every surface here is both slave
     * and master: a contact that one side sees at a glancing angle, the other side sees
     * head-on, and it is found there instead.
     */
    private static final double FACING = 0.25;

    /**
     * Normal contact damping, as a fraction of critical for the penalty pair.
     *
     * <p>Zero, and the history is worth keeping. A damper was added here on the theory that a
     * pure penalty spring cannot survive a resting contact, because a block sliding across
     * another was gaining 245 % of its energy. The real cause was a spurious contact pair
     * between two perpendicular surfaces; with that fixed, the same case closes to 0.15 %
     * undamped, and the damper makes every other measurement slightly worse.
     *
     * <p>It stays because it is the coefficient of restitution in disguise and a caller may
     * want a deader bounce than elasticity gives. It is not needed for stability and not
     * needed for resting contact.
     */
    public static final double DEFAULT_DAMPING = 0.0;

    private final QuadMesh mesh;
    private QuadMesh.Surface surface;
    /** Elements deleted from the solve; the surface is rebuilt around them. */
    private final boolean[] dead;
    private boolean stale;

    /** Segment endpoints, 2 per segment; the same array the surface reported. */
    private int[] seg;
    /** Surface nodes, ascending. */
    private int[] nodes;
    /**
     * Segments incident on each surface node, as a compressed row: {@code segStart[i]} to
     * {@code segStart[i+1]} index into {@code segOf}. Used to exclude a node from its own
     * segments and from its neighbours', which is what stops a surface contacting itself
     * everywhere it is merely continuous.
     */
    private int[] segStart;
    private int[] segOf;
    /** Position of a node in {@link #nodes}, or -1 for an interior node. */
    private final int[] slot;

    private double friction;
    private double scale = DEFAULT_SCALE;
    private double damping = DEFAULT_DAMPING;

    /** The current step, and the one the stiffness is built from. See {@link #resolve}. */
    private double step;
    private double reference;

    // Scratch, sized once. Current coordinates of the surface only.
    private final double[] cr;
    private final double[] cz;
    private double[] nr;
    private double[] nz;
    private double[] length;
    /** Outward normal at each surface node, the mean of the segments meeting there. */
    private double[] pnr;
    private double[] pnz;
    /** Master segment each surface node had last step, or -1. Only used to spot tunnelling. */
    private int[] lastMaster;

    // A uniform grid over the surface, rebuilt each step. Buckets hold segment indices.
    private double cellSize;
    private int[] bucketStart = new int[0];
    private int[] bucketOf = new int[0];
    private int buckets;
    private double gridR0;
    private double gridZ0;
    private int gridCols;
    private int gridRows;

    private int pairs;
    private double energy;
    private double dissipation;
    private double maxPenetration;
    private int escaped;
    /**
     * Per node: the master the search chose, and whether it had to drop a candidate for depth.
     * One slot each, written only by the search and read only by the accumulation, which is
     * what lets the search be threaded without the answer depending on the schedule.
     */
    private int[] master;
    private boolean[] dropped;

    public Contact(QuadMesh mesh) {
        this.mesh = mesh;
        this.dead = new boolean[mesh.elementCount];
        this.cr = new double[mesh.nodeCount];
        this.cz = new double[mesh.nodeCount];
        this.slot = new int[mesh.nodeCount];
        adopt();
    }

    /**
     * Rebuilds everything derived from the free surface.
     *
     * <p>Called once at construction and again whenever an element is deleted, because
     * deleting one exposes the faces it was hiding -- a body eroded through the middle has two
     * new surfaces that were interior a moment earlier, and contact that did not know about
     * them would let the halves pass through each other.
     */
    private void adopt() {
        this.surface = mesh.surface(dead);
        this.seg = surface.edges();
        this.nodes = surface.nodes();

        java.util.Arrays.fill(slot, -1);
        for (int i = 0; i < nodes.length; i++) slot[nodes[i]] = i;

        final int segments = surface.edgeCount();
        this.segStart = new int[nodes.length + 1];
        for (int s = 0; s < segments; s++) {
            segStart[slot[seg[2 * s]] + 1]++;
            segStart[slot[seg[2 * s + 1]] + 1]++;
        }
        for (int i = 0; i < nodes.length; i++) segStart[i + 1] += segStart[i];
        this.segOf = new int[segStart[nodes.length]];
        final int[] fill = segStart.clone();
        for (int s = 0; s < segments; s++) {
            segOf[fill[slot[seg[2 * s]]]++] = s;
            segOf[fill[slot[seg[2 * s + 1]]]++] = s;
        }

        this.nr = new double[segments];
        this.nz = new double[segments];
        this.length = new double[segments];
        this.pnr = new double[nodes.length];
        this.pnz = new double[nodes.length];
        this.lastMaster = new int[nodes.length];
        java.util.Arrays.fill(lastMaster, -1);
        this.master = new int[nodes.length];
        this.dropped = new boolean[nodes.length];
        this.stale = false;
    }

    /**
     * Tells this that an element has been deleted, so the surface is rebuilt before the next
     * step. Marked rather than rebuilt immediately, because a step that erodes a hundred
     * elements should rebuild once.
     */
    void forget(int element) {
        dead[element] = true;
        stale = true;
    }

    /** Coulomb coefficient between surfaces. Zero is frictionless. */
    public Contact withFriction(double mu) {
        if (mu < 0.0) throw new IllegalArgumentException("friction must not be negative");
        this.friction = mu;
        return this;
    }

    /** Penalty stiffness as a fraction of {@code m/dt^2}. See {@link #DEFAULT_SCALE}. */
    /**
     * Normal contact damping as a fraction of critical. Zero is a perfectly elastic pair and
     * is not recommended; see {@link #DEFAULT_DAMPING}.
     */
    public Contact withDamping(double zeta) {
        if (zeta < 0.0 || zeta > 1.0) {
            throw new IllegalArgumentException(
                    "contact damping is a fraction of critical, so it lies in [0, 1]; was " + zeta);
        }
        this.damping = zeta;
        return this;
    }

    public Contact withStiffnessScale(double s) {
        if (!(s > 0.0) || s >= 2.0) {
            throw new IllegalArgumentException(
                    "the stiffness scale must lie in (0, 2); a pair is found from both sides, "
                            + "so the effective stiffness is twice this and above 2 the contact "
                            + "frequency exceeds what central difference can integrate, was " + s);
        }
        this.scale = s;
        return this;
    }

    public double friction() {
        return friction;
    }

    /** Node-segment pairs in contact at the last step. */
    public int pairCount() {
        return pairs;
    }

    /** Energy currently stored in the penalty springs, joules. Recoverable, not lost. */
    public double energy() {
        return energy;
    }

    /** Energy taken out by friction over the whole run, joules. */
    public double dissipation() {
        return dissipation;
    }

    /** Deepest penetration at the last step, metres. The honesty check on the stiffness. */
    public double maxPenetration() {
        return maxPenetration;
    }

    /**
     * Nodes dropped over the whole run for having passed too far through a surface.
     *
     * <p>Cumulative, not per step, and that is the point: a node escapes for one step and is
     * gone, so a counter reset each step would read zero everywhere except in the single step
     * it happened, and a test asserting on it would pass by missing it.
     */
    public int escapedNodes() {
        return escaped;
    }

    /** Surface segments this is watching. */
    public int segmentCount() {
        return surface.edgeCount();
    }

    /**
     * Adds contact forces for the current configuration.
     *
     * <p>Called after the internal forces have been assembled and before they are integrated,
     * so the penalty acts over the same step the material does.
     */
    /**
     * Prepares the current configuration, ready for {@link #search} and {@link #accumulate}.
     *
     * <p>Split in three because the middle part is most of the cost and is the only part that
     * can be threaded. A search writes nothing but its own node's answer, so any number of
     * workers can do it in any order; accumulating adds into shared force arrays and running
     * totals, where the order of a floating-point sum is part of the answer. Threading the
     * search and keeping the sum in node order buys the speed without giving up the property
     * that two runs of the same scene agree bit for bit.
     *
     * @return how many surface nodes there are to search over, or 0 when there is nothing
     */
    int begin(double[] ur, double[] uz, double dt, double referenceStep) {
        this.step = dt;
        this.reference = referenceStep;
        if (stale) adopt();
        pairs = 0;
        energy = 0.0;
        maxPenetration = 0.0;
        if (surface.edgeCount() == 0) return 0;

        geometry(ur, uz);
        grid();
        return nodes.length;
    }

    /** The master segment for a range of surface nodes. Writes only their own slots. */
    void search(int from, int to) {
        for (int i = from; i < to; i++) {
            master[i] = nearestPenetrated(i, nodes[i]);
        }
    }

    /** Applies what the search found, in node order. Serial, and deliberately so. */
    void accumulate(double[] vr, double[] vz, double[] mass, double[] fr, double[] fz) {
        for (int i = 0; i < nodes.length; i++) {
            final int n = nodes[i];
            final int best = master[i];
            if (best < 0) {
                // Tunnelling, and the only way to tell it from an ordinary separation: this
                // node was resolving a contact last step and now has a candidate it is too
                // deep behind to resolve. Counting every too-deep candidate instead would
                // count the far side of every body, every step, for ever.
                if (lastMaster[i] >= 0 && dropped[i]) escaped++;
                lastMaster[i] = -1;
                continue;
            }
            lastMaster[i] = best;
            resolve(i, n, best, vr, vz, mass, fr, fz);
        }
    }

    /** Current surface coordinates, segment normals and lengths. */
    private void geometry(double[] ur, double[] uz) {
        for (int n : nodes) {
            cr[n] = mesh.r[n] + ur[n];
            cz[n] = mesh.z[n] + uz[n];
        }
        for (int s = 0; s < surface.edgeCount(); s++) {
            final int a = seg[2 * s], b = seg[2 * s + 1];
            final double dr = cr[b] - cr[a], dz = cz[b] - cz[a];
            final double len = Math.sqrt(dr * dr + dz * dz);
            length[s] = len;
            if (len > 0.0) {
                // The edge is wound as its owning element winds, counter-clockwise, so the
                // tangent turned by -90 degrees points out of the solid. Getting this backwards
                // does not look like a bug -- it looks like an attractive force.
                nr[s] = dz / len;
                nz[s] = -dr / len;
            } else {
                nr[s] = 0.0;
                nz[s] = 0.0;
            }
        }

        // Which way each node faces: the mean of the segments meeting there, so a node on a
        // flat face points straight out and a node on a corner points into the corner.
        for (int i = 0; i < nodes.length; i++) {
            double sr = 0.0, sz = 0.0;
            for (int k = segStart[i]; k < segStart[i + 1]; k++) {
                sr += nr[segOf[k]];
                sz += nz[segOf[k]];
            }
            final double len = Math.sqrt(sr * sr + sz * sz);
            pnr[i] = len > 0.0 ? sr / len : 0.0;
            pnz[i] = len > 0.0 ? sz / len : 0.0;
        }
    }

    /**
     * A uniform grid over the surface, rebuilt each step by counting sort.
     *
     * <p>Rebuilt rather than updated because the surface is a small fraction of the mesh and
     * two passes over it cost less than deciding what has moved.
     *
     * <p>The cell is <b>1.5</b> times the longest segment, which is not a round number chosen
     * for comfort. A node close enough to a segment to be touching it is at most half a
     * segment behind it and within its projection span, so its distance to that segment's
     * endpoints is at most sqrt(0.5^2 + 1^2) = 1.12 segment lengths. The nine cells around a
     * node therefore have to reach at least that far, and 1.5 does with margin. Two also
     * works and searches 78 % more area for nothing: dropping to 1.5 left every number in the
     * slug scene bit-identical and took 5 % off the threaded runtime.
     */
    private void grid() {
        final int segments = surface.edgeCount();
        double rlo = Double.MAX_VALUE, zlo = Double.MAX_VALUE;
        double rhi = -Double.MAX_VALUE, zhi = -Double.MAX_VALUE;
        double longest = 0.0;
        for (int n : nodes) {
            if (cr[n] < rlo) rlo = cr[n];
            if (cr[n] > rhi) rhi = cr[n];
            if (cz[n] < zlo) zlo = cz[n];
            if (cz[n] > zhi) zhi = cz[n];
        }
        for (int s = 0; s < segments; s++) if (length[s] > longest) longest = length[s];

        cellSize = longest > 0.0 ? 1.5 * longest : 1.0;
        gridR0 = rlo;
        gridZ0 = zlo;
        gridCols = Math.max(1, (int) ((rhi - rlo) / cellSize) + 1);
        gridRows = Math.max(1, (int) ((zhi - zlo) / cellSize) + 1);

        // A segment is filed under both of its endpoints' cells, so a query that visits the
        // cells around a node sees every segment either of whose ends could be near it.
        buckets = gridCols * gridRows;
        if (bucketStart.length < buckets + 1) bucketStart = new int[buckets + 1];
        else java.util.Arrays.fill(bucketStart, 0, buckets + 1, 0);
        if (bucketOf.length < 2 * segments) bucketOf = new int[2 * segments];

        for (int s = 0; s < segments; s++) {
            bucketStart[cellOf(cr[seg[2 * s]], cz[seg[2 * s]]) + 1]++;
            bucketStart[cellOf(cr[seg[2 * s + 1]], cz[seg[2 * s + 1]]) + 1]++;
        }
        for (int c = 0; c < buckets; c++) bucketStart[c + 1] += bucketStart[c];
        final int[] fill = java.util.Arrays.copyOf(bucketStart, buckets + 1);
        for (int s = 0; s < segments; s++) {
            bucketOf[fill[cellOf(cr[seg[2 * s]], cz[seg[2 * s]])]++] = s;
            bucketOf[fill[cellOf(cr[seg[2 * s + 1]], cz[seg[2 * s + 1]])]++] = s;
        }
    }

    private int cellOf(double r, double z) {
        final int i = clamp((int) ((r - gridR0) / cellSize), gridCols);
        final int j = clamp((int) ((z - gridZ0) / cellSize), gridRows);
        return j * gridCols + i;
    }

    private static int clamp(int v, int n) {
        return v < 0 ? 0 : v >= n ? n - 1 : v;
    }

    /**
     * The deepest segment this node has crossed, or -1.
     *
     * <p>Visits the nine cells around the node in a fixed order and, within a cell, the
     * segments in the order the counting sort filed them, which is segment order. Both are
     * deterministic, so the tie-break below picks the same segment on every run.
     */
    private int nearestPenetrated(int i, int n) {
        final int ci = clamp((int) ((cr[n] - gridR0) / cellSize), gridCols);
        final int cj = clamp((int) ((cz[n] - gridZ0) / cellSize), gridRows);

        int best = -1;
        double deepest = 0.0;
        boolean fell = false;

        for (int dj = -1; dj <= 1; dj++) {
            final int j = cj + dj;
            if (j < 0 || j >= gridRows) continue;
            for (int di = -1; di <= 1; di++) {
                final int ii = ci + di;
                if (ii < 0 || ii >= gridCols) continue;
                final int cell = j * gridCols + ii;
                for (int k = bucketStart[cell]; k < bucketStart[cell + 1]; k++) {
                    final int s = bucketOf[k];
                    if (excluded(i, n, s)) continue;
                    final double depth = depth(n, s);
                    if (depth <= 0.0) continue;
                    if (depth > MAX_DEPTH * length[s]) {
                        // Only the segment this node was actually resolving counts. Every
                        // other too-deep candidate is the far side of some body, which is
                        // behind every node of it, every step, for ever.
                        if (s == lastMaster[i]) fell = true;
                        continue;
                    }
                    if (depth > deepest) {
                        deepest = depth;
                        best = s;
                    }
                }
            }
        }
        this.dropped[i] = fell;
        return best;
    }

    /**
     * Whether a node is too close to a segment in the surface's own topology to contact it.
     *
     * <p>A node is on its own segments and is one step around the outline from its neighbours'.
     * Neither is a collision; both would fire constantly, because a surface is always exactly
     * touching itself. Excluding two rings rather than one matters at a sharp corner, where
     * the segment on the far side of the neighbour can be close enough to register while still
     * being the same piece of material.
     */
    private boolean excluded(int i, int n, int s) {
        // Two surfaces in contact face each other, so a node can only be touching a segment
        // whose outward normal opposes its own. Without this the test is far too generous:
        // the region behind a segment is a half-plane, not a box, so a node sitting on the
        // bottom face of a block is "two millimetres deep" behind its own left face, at a
        // projection that lands squarely inside a real segment. The depth cap hides that on an
        // undeformed lattice -- the faces are a whole cell apart and the cap is half of one --
        // but a deforming body brings them inside it, and because the master is chosen by
        // greatest depth, the spurious pair would win over the real one and drive the node
        // sideways with several times the right force.
        final int a = seg[2 * s], b = seg[2 * s + 1];
        if (pnr[i] * nr[s] + pnz[i] * nz[s] > -FACING) return true;
        if (a == n || b == n) return true;
        for (int k = segStart[i]; k < segStart[i + 1]; k++) {
            final int own = segOf[k];
            final int p = seg[2 * own], q = seg[2 * own + 1];
            if (p == a || p == b || q == a || q == b) return true;
        }
        return false;
    }

    /** How far a node lies behind a segment, or 0 if it is outside or off the end. */
    private double depth(int n, int s) {
        final int a = seg[2 * s], b = seg[2 * s + 1];
        final double len = length[s];
        if (len <= 0.0) return 0.0;
        final double tr = (cr[b] - cr[a]) / len, tz = (cz[b] - cz[a]) / len;
        final double dr = cr[n] - cr[a], dz = cz[n] - cz[a];
        final double along = dr * tr + dz * tz;
        if (along < 0.0 || along > len) return 0.0;
        final double gap = dr * nr[s] + dz * nz[s];
        return gap < 0.0 ? -gap : 0.0;
    }

    /** Where along a segment a node projects, as a fraction from the first endpoint. */
    private double parameter(int n, int s) {
        final int a = seg[2 * s];
        final double len = length[s];
        final double tr = (cr[seg[2 * s + 1]] - cr[a]) / len;
        final double tz = (cz[seg[2 * s + 1]] - cz[a]) / len;
        return ((cr[n] - cr[a]) * tr + (cz[n] - cz[a]) * tz) / len;
    }

    private void resolve(int i, int n, int s, double[] vr, double[] vz,
                         double[] mass, double[] fr, double[] fz) {
        final double dt = step;
        final int a = seg[2 * s], b = seg[2 * s + 1];
        final double depth = depth(n, s);
        final double t = parameter(n, s);
        final double wa = 1.0 - t, wb = t;

        // The segment point moves with a mass that depends on where along it the contact sits.
        // At an end it is one node's mass; in the middle the two share, and the shape-function
        // weights are what say how. Squaring them is what makes the pair conserve energy: the
        // same weights distribute the force and gather the velocity.
        // Fixed half weights, not the live ones. The mass a distributed force sees on the
        // segment side is 1 / (wa^2/m_a + wb^2/m_b), which runs from one node's mass at an end
        // to twice that in the middle -- and the stiffness below is built from it, so using
        // the live weights would let a node change its own spring constant by a factor of two
        // while the spring was loaded. A spring whose constant moves under load is an energy
        // source. On the aligned lattices this mesher produces, contact points sit near
        // segment ends and a whisker of lateral motion flips them, so that was worth 4 % of a
        // head-on collision. The live weights still distribute the force, which is what keeps
        // momentum exact.
        final double invMaster = 0.25 / mass[a] + 0.25 / mass[b];
        if (!(invMaster > 0.0)) return;
        final double mMaster = 1.0 / invMaster;
        final double mSlave = mass[n];
        final double mStar = mSlave * mMaster / (mSlave + mMaster);

        // Built from the REFERENCE step, not the current one. The two differ as soon as an
        // element compresses, and a spring whose constant rises while it is loaded is not a
        // conservative element: it manufactures energy at the rate the stiffness climbs. That
        // was worth 6.5 % of the collision before it was pinned. Stability survives the change
        // because the current step only ever falls below the reference one as the mesh
        // crushes, and a falling step makes the contact frequency safer rather than riskier.
        final double k = scale * mStar / (reference * reference);

        final double nR = nr[s], nZ = nz[s];
        // The reduced mass again, but with the live weights this time. Relative motion at this
        // particular contact point is what the damper and the friction are about, not a spring
        // constant that has to hold still while it is loaded, so the reason for freezing the
        // weights above does not apply to either and the accurate value is the right one.
        final double mTrue = 1.0 / (wa * wa / mass[a] + wb * wb / mass[b]);
        final double mTangent = mSlave * mTrue / (mSlave + mTrue);

        // Closing speed along the normal. Negative while the surfaces are still driving into
        // one another, positive as they come apart.
        final double closing = (vr[n] - wa * vr[a] - wb * vr[b]) * nR
                + (vz[n] - wa * vz[a] - wb * vz[b]) * nZ;

        // An optional dashpot across the same pair, as a fraction of critical:
        // c = 2*zeta*sqrt(k*m*), which with the mass-based stiffness above is
        // 2*zeta*m*sqrt(SCALE)/dt0 and is again independent of the materials. Off by default;
        // see DEFAULT_DAMPING for why it is here at all.
        final double damper = 2.0 * damping * mStar * Math.sqrt(scale) / reference;
        // Never pull. A damper that outruns the spring on the way out would suck a separating
        // node back, which is the one thing contact must not do.
        final double normal = Math.max(0.0, k * depth - damper * closing);

        fr[n] += normal * nR;
        fz[n] += normal * nZ;
        fr[a] -= wa * normal * nR;
        fz[a] -= wa * normal * nZ;
        fr[b] -= wb * normal * nR;
        fz[b] -= wb * normal * nZ;

        pairs++;
        energy += 0.5 * k * depth * depth;
        if (depth > maxPenetration) maxPenetration = depth;
        // What the dashpot took out, over the mean closing speed of the step for the same
        // reason the friction term below uses a mean: the force changes the speed it is acting
        // against, and the starting value would claim twice the work.
        final double damperForce = normal - k * depth;
        if (damperForce != 0.0) {
            final double closingAfter = closing + damperForce * dt / mTangent;
            dissipation -= damperForce * 0.5 * (closing + closingAfter) * dt;
        }

        if (friction <= 0.0) return;

        // Coulomb, as a return map. The trial is the force that would remove all relative
        // sliding this step; if it is outside the cone, it comes back to the cone's surface
        // and the surfaces slide.
        final double tR = -nZ, tZ = nR;
        final double slide = (vr[n] - wa * vr[a] - wb * vr[b]) * tR
                + (vz[n] - wa * vz[a] - wb * vz[b]) * tZ;
        double tangent = -STICK_RELAXATION * mTangent * slide / dt;
        final double cone = friction * normal;
        if (tangent > cone) tangent = cone;
        else if (tangent < -cone) tangent = -cone;

        fr[n] += tangent * tR;
        fz[n] += tangent * tZ;
        fr[a] -= wa * tangent * tR;
        fz[a] -= wa * tangent * tZ;
        fr[b] -= wb * tangent * tR;
        fz[b] -= wb * tangent * tZ;

        // Work done against the relative motion, over the mean sliding speed of the step
        // rather than its starting value. With the mean, a pair that sticks accounts for
        // exactly the relative kinetic energy the stick removed; with the starting value it
        // would claim twice that, and the energy audit would never close.
        final double after = slide + tangent * dt / mTangent;
        dissipation += -tangent * 0.5 * (slide + after) * dt;
    }
}
