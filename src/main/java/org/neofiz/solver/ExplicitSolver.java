package org.neofiz.solver;

import org.neofiz.core.Formulation;
import org.neofiz.core.Material;
import org.neofiz.core.Materials;
import org.neofiz.mesh.QuadMesh;

import java.util.Arrays;

/**
 * Two-dimensional Lagrangian explicit finite element solver. Four-node quadrilaterals,
 * central-difference time integration, full or reduced Gauss integration, J2 plasticity.
 *
 * <p><b>No mass scaling, ever.</b> It is the standard trick for making explicit dynamics
 * affordable and it corrupts exactly the inertial behaviour this project is about.
 * Fragment velocities and bulge dynamics both depend on real mass. The offline-solve
 * decision is what buys the right to refuse it.
 *
 * <p><b>Why B is recomputed every step.</b> For small-strain physics the strain-displacement
 * matrix could be cached, which would be several times faster. It is deliberately not
 * cached, because the solver this becomes recomputes B from the deformed configuration every
 * step, and a timing number obtained from a cached B would flatter the cost model that M0
 * exists to test. The same applies to the hourglass shape vectors.
 *
 * <h2>Why the solver is incremental</h2>
 *
 * M0 computed stress from total displacement: state in, stress out, no memory. That stops
 * working the moment the material can yield, because a plastic stress depends on the path
 * taken and not merely on where the material ended up. So stress, back stress and
 * accumulated plastic strain are stored per Gauss point and integrated forward along the
 * strain path, which is what an explicit code actually is.
 *
 * <p>The step order matters and is the standard one: integrate velocity, integrate
 * displacement, take the strain increment from the displacement just applied, update the
 * material, then assemble forces for the next step. Assembling before the material update
 * would advance the mesh and the stress out of phase by a step.
 */
public final class ExplicitSolver {

    /** 2x2 Gauss points in natural coordinates, all weights 1. */
    private static final double G = 0.5773502691896258; // 1/sqrt(3)
    private static final double[] GP_XI = {-G, G, G, -G};
    private static final double[] GP_ETA = {-G, -G, G, G};

    /**
     * The 1-point rule's weight. The natural domain is [-1,1] squared, so a single centroid
     * sample carries the whole area of 4; the 2x2 rule splits the same 4 into four weights
     * of 1. Getting this wrong quarters every element volume and is silent.
     */
    private static final double W_REDUCED = 4.0;

    /**
     * The hourglass mode of a bilinear quad: the nodal pattern that alternates sign around
     * the element. It is the only displacement pattern a centroid sample cannot see.
     */
    private static final double[] H = {1.0, -1.0, 1.0, -1.0};

    /** Default Flanagan-Belytschko stiffness scaling. In the usual 0.01 to 0.10 range. */
    private static final double DEFAULT_HOURGLASS_COEFFICIENT = 0.05;

    private final QuadMesh mesh;
    private final Materials materials;
    private final Formulation formulation;
    private final Integration integration;
    private final Kinematics kinematics;
    /** Out-of-plane thickness for planar formulations, metres. Unused when axisymmetric. */
    private final double thickness;

    /**
     * Elastic constants and flow law, one entry per element.
     *
     * <p>Unrolled out of {@link Materials} at construction rather than read through it in the
     * kernel, because the inner loop wants a double out of an array and not a record
     * dereference behind a branch. A uniform mesh fills every entry with the same number, so
     * the arithmetic is unchanged and the cost is three doubles and a reference per element.
     */
    private final double[] lambdaOf;
    private final double[] muOf;
    private final double[] bulkOf;
    private final double[] densityOf;
    private final J2.Flow[] flowOf;

    /** The fastest dilatational wave on the mesh, m/s. What the CFL step is sized against. */
    private final double fastestWave;

    /**
     * The timestep in the reference configuration. Fixed, unlike {@link #dt}, which falls as
     * elements compress -- and a penalty stiffness has to be built from something that does
     * not move while a contact is loaded, or it manufactures energy. See {@link Contact}.
     */
    private double dt0;

    /**
     * Per-point failure strain and softening state, or null for a material that cannot fail.
     * The one piece of solver state that is not uniform across the mesh, which is the point:
     * it is where the defect field enters.
     */
    private Damage damage;

    // Nodal state. Structure-of-arrays.
    private final double[] ur;
    private final double[] uz;
    private final double[] vr;
    private final double[] vz;
    private final double[] fr;
    private final double[] fz;
    private final double[] mass;

    private final boolean[] fixedR;
    private final boolean[] fixedZ;

    /**
     * Nodes currently resting on the rigid wall. Not a boundary condition but a state: the
     * contact set, which changes every step as the mushroom spreads and again as the
     * specimen lets go. Held as nodal state alongside the rest of it rather than inside
     * {@link RigidWall} so that the structure-of-arrays layout stays intact.
     */
    private final boolean[] onWall;

    // Material state, one set per Gauss point. This is the memory plasticity requires.
    private final int pointsPerElement;
    private final double[] sig;   // stride 4: rr, zz, tt, rz
    private final double[] back;  // stride 4, deviatoric
    private final double[] epsP;  // stride 1, accumulated equivalent plastic strain
    private final double[] wPlastic; // stride 1, plastic work per unit volume, J/m^3

    /**
     * Hourglass restoring force per element, accumulated incrementally.
     *
     * <p>Held as a force rather than recomputed from total displacement each step, because
     * under finite strain the shape vector is rebuilt on a mesh that has moved, so there is
     * no fixed total to multiply. At small strain the two are the same quantity.
     */
    private final double[] hgR;
    private final double[] hgZ;
    /** Hourglass stiffness per element from the last assembly, for the energy audit. */
    private final double[] hourglassStiffness;

    /**
     * Every force any source produces, before anything is summed. Stride 8: four radial
     * values then four axial. An element owns slot {@code e} and fills all four corners; a
     * pressure edge owns slot {@code elementCount + k} and fills the first two.
     *
     * <p>This is what makes the assembly safe to thread. Two elements sharing a node would
     * otherwise write the same nodal slot from two threads, and the usual answers -- an atomic
     * add, or a private force vector per thread reduced afterwards -- both change the
     * <em>order</em> the contributions to that node are summed in. Floating-point addition is
     * not associative, so both give an answer that depends on the thread count, and a
     * nonlinear plastic run amplifies a last-bit difference over 200,000 steps until it is
     * visible. Writing per source and gathering per node fixes the summation order, which is
     * what {@link #setThreads} can then be held to bit-for-bit.
     *
     * <p>The pressure load shares the array rather than keeping its own, so that it is
     * assembled in the same parallel region as the elements and summed in the same pass.
     * Kept separate it would be a third of the threaded step, because a thousand edges are
     * too little work to be worth a barrier and too much to leave on one thread.
     */
    private final double[] contribution;

    /** Where the pressure edges start in {@link #contribution}, in slots. */
    private final int edgeSlotBase;

    /**
     * Sources incident on each node, as offsets into {@link #contribution}: the elements in
     * ascending element order, then the pressure edges in ascending edge order. That is the
     * order the serial assembly added them in, term for term. Compressed-row storage: node
     * {@code n} owns {@code [nodeSlotStart[n], nodeSlotStart[n + 1])}.
     */
    private final int[] nodeSlotStart;
    private final int[] nodeSlot;

    /** Shortest edge seen by each worker, squared. Stride {@link #PAD} to avoid false sharing. */
    private double[] minEdgeLocal;

    private static final int PAD = 8;

    /** Worker pool, or null when the solver runs on the calling thread alone. */
    private Parallel pool;

    /** Published to the workers through the pool, under the ticket's release. */
    private double kernelStrainDt;
    private double kernelPressure;
    private double nodalDecay;
    private double nodalGain;
    private double nodalDtAvg;

    private final Parallel.Kernel assembleKernel = this::assemble;
    private final Parallel.Kernel gatherKernel = (from, to, worker) -> gather(from, to);
    private final Parallel.Kernel nodalKernel = (from, to, worker) -> nodalUpdate(from, to);

    private final double cflSafety;
    private double dt;
    private double dtPrevious;
    /** Shortest current element edge, squared, tracked during the element loop. */
    private double minEdgeSquared;
    private double time;
    private long stepCount;

    /** Mass-proportional damping coefficient, 1/s. Used to relax toward a static answer. */
    private double damping;

    private double targetPressure;
    private double rampDuration;

    /** Flanagan-Belytschko stiffness scaling. Zero disables hourglass control entirely. */
    private double hourglassCoefficient = DEFAULT_HOURGLASS_COEFFICIENT;

    /** Rigid frictionless anvil, or null for an unconstrained body. */
    private RigidWall wall;

    /** Penalty contact between deformable surfaces, or null if nothing may touch. */
    private Contact contact;

    /** The threaded half of the contact pass; see Contact for why only half. */
    private final Parallel.Kernel contactKernel = (from, to, worker) -> contact.search(from, to);

    /** Uniform body acceleration, m/s^2. Zero unless a scene asks for it. */
    private double gravityR;
    private double gravityZ;

    /**
     * Which elements have been deleted, or null when erosion is off.
     *
     * <p>Null rather than an all-false array so that the check in the hot loop is a reference
     * test against a field that is almost always null, and every validation gate in this
     * project -- none of which erodes -- runs the same arithmetic it always did.
     */
    private boolean[] gone;
    /** How much heavier the mesh has been made; see setMassScaling. */
    private double massScale = 1.0;
    private double erosionDamage = Double.POSITIVE_INFINITY;
    private int erodedCount;
    private double erodedEnergy;
    private int invertedCount;

    /** Full 2x2 integration, small strain. */
    public ExplicitSolver(QuadMesh mesh, Material material, Formulation formulation,
                          double thickness, double cflSafety) {
        this(mesh, material, formulation, Integration.FULL, thickness, cflSafety);
    }

    /** Small strain. */
    public ExplicitSolver(QuadMesh mesh, Material material, Formulation formulation,
                          Integration integration, double thickness, double cflSafety) {
        this(mesh, material, formulation, integration, Kinematics.SMALL_STRAIN,
                thickness, cflSafety);
    }

    public ExplicitSolver(QuadMesh mesh, Material material, Formulation formulation,
                          Integration integration, Kinematics kinematics,
                          double thickness, double cflSafety) {
        this(mesh, Materials.uniform(material), formulation, integration, kinematics,
                thickness, cflSafety);
    }

    /**
     * The general form: a mesh whose elements need not be made of the same thing.
     *
     * <p>Everything downstream of here is already per element -- stress, plastic strain, back
     * stress and damage all live per integration point -- so a second substance costs the
     * kernel one array index and nothing else. What it buys is the whole point of a sandbox:
     * a steel bar dropped on a concrete slab is one solve, not two coupled ones.
     */
    public ExplicitSolver(QuadMesh mesh, Materials materials, Formulation formulation,
                          Integration integration, Kinematics kinematics,
                          double thickness, double cflSafety) {
        materials.requireCovers(mesh.elementCount);
        for (Material m : materials.distinct()) {
            if (formulation == Formulation.PLANE_STRESS && !m.isElastic()) {
                // Enforcing sigma_zz = 0 through a return map is a different algorithm, not a
                // special case of this one. Refusing is better than silently solving plane
                // strain and labelling it plane stress.
                throw new IllegalArgumentException(
                        "plane stress plasticity needs its own return map; not implemented");
            }
        }
        this.mesh = mesh;
        this.materials = materials;
        this.formulation = formulation;
        this.integration = integration;
        this.kinematics = kinematics;
        this.thickness = thickness;
        this.cflSafety = cflSafety;
        this.lambdaOf = new double[mesh.elementCount];
        this.muOf = new double[mesh.elementCount];
        this.bulkOf = new double[mesh.elementCount];
        this.densityOf = new double[mesh.elementCount];
        this.flowOf = new J2.Flow[mesh.elementCount];
        for (int e = 0; e < mesh.elementCount; e++) {
            final Material m = materials.at(e);
            final double lam = formulation == Formulation.PLANE_STRESS
                    ? 2.0 * m.lambda() * m.mu() / (m.lambda() + 2.0 * m.mu())
                    : m.lambda();
            lambdaOf[e] = lam;
            muOf[e] = m.mu();
            bulkOf[e] = lam + 2.0 * m.mu() / 3.0;
            densityOf[e] = m.density();
            flowOf[e] = J2.Flow.of(m);
        }
        this.fastestWave = materials.fastestWave();

        int n = mesh.nodeCount;
        this.ur = new double[n];
        this.uz = new double[n];
        this.vr = new double[n];
        this.vz = new double[n];
        this.fr = new double[n];
        this.fz = new double[n];
        this.mass = new double[n];
        this.fixedR = new boolean[n];
        this.fixedZ = new boolean[n];
        this.onWall = new boolean[n];

        this.pointsPerElement = integration.gaussPoints();
        int points = mesh.elementCount * pointsPerElement;
        this.sig = new double[J2.COMPONENTS * points];
        this.back = new double[J2.COMPONENTS * points];
        this.epsP = new double[points];
        this.wPlastic = new double[points];
        this.hgR = new double[mesh.elementCount];
        this.hgZ = new double[mesh.elementCount];
        this.hourglassStiffness = new double[mesh.elementCount];

        this.edgeSlotBase = mesh.elementCount;
        this.contribution = new double[8 * (mesh.elementCount + mesh.pressureEdges.length / 2)];
        this.nodeSlotStart = new int[n + 1];
        this.nodeSlot = new int[mesh.conn.length + mesh.pressureEdges.length];
        buildNodeSlots();
        this.minEdgeLocal = new double[PAD];

        lumpMass();

        // CFL: a dilatational wave must not cross the shortest element edge in one step.
        double edge = mesh.minimumEdgeLength();
        this.minEdgeSquared = edge * edge;
        this.dt = cflSafety * edge / fastestWave;
        this.dt0 = this.dt;
        this.dtPrevious = this.dt;
    }

    /**
     * Inverts the connectivity: for each node, the sources incident on it -- element corners
     * in ascending element order, then pressure edges in ascending edge order. That ordering
     * is the whole point. It is the order the serial assembly added them in, so the gather
     * reproduces the serial sum term for term rather than merely to a tolerance.
     */
    private void buildNodeSlots() {
        final int[] conn = mesh.conn;
        final int[] edges = mesh.pressureEdges;

        for (int k = 0; k < conn.length; k++) nodeSlotStart[conn[k] + 1]++;
        for (int k = 0; k < edges.length; k++) nodeSlotStart[edges[k] + 1]++;
        for (int i = 0; i < mesh.nodeCount; i++) nodeSlotStart[i + 1] += nodeSlotStart[i];

        final int[] cursor = nodeSlotStart.clone();
        for (int e = 0; e < mesh.elementCount; e++) {
            for (int c = 0; c < 4; c++) {
                nodeSlot[cursor[conn[e * 4 + c]]++] = e * 8 + c;
            }
        }
        for (int k = 0; k < edges.length; k += 2) {
            final int slot = (edgeSlotBase + k / 2) * 8;
            nodeSlot[cursor[edges[k]]++] = slot;
            nodeSlot[cursor[edges[k + 1]]++] = slot + 1;
        }
    }

    // ------------------------------------------------------------------ setup

    /**
     * Splits the element kernel and the nodal gather across {@code n} threads. One means the
     * calling thread does everything and no pool exists.
     *
     * <p><b>The answer does not change.</b> Not "changes by a rounding error" -- the element
     * kernel writes only its own element's slot, the gather reads each node's contributions in
     * ascending element order whatever thread produced them, and the chunking is fixed rather
     * than stolen. So the stress field is bit-for-bit identical at any thread count, which is
     * asserted rather than assumed. The nodal integration, the contact set and the pressure
     * load stay on the calling thread: together they are a few per cent of a step, and contact
     * in particular carries an impulse accumulator that would need an ordering rule of its own
     * for a gain that is not there.
     *
     * <p>Workers spin between regions rather than sleeping, so a pool left alive occupies its
     * cores. Call {@link #shutdown()} when the run is done.
     */
    public void setThreads(int n) {
        if (n < 1) throw new IllegalArgumentException("need at least one thread");
        if (pool != null && pool.threads() == n) return;
        shutdown();
        this.minEdgeLocal = new double[PAD * n];
        if (n > 1) this.pool = new Parallel(n);
    }

    /** Stops the worker pool. Safe to call more than once, and on a solver that never had one. */
    public void shutdown() {
        if (pool != null) {
            pool.close();
            pool = null;
        }
    }

    /** Threads the element kernel is split across. One unless {@link #setThreads} says otherwise. */
    public int threads() {
        return pool == null ? 1 : pool.threads();
    }

    /** Pin the radial degree of freedom at a node. */
    public void fixRadial(int node) {
        fixedR[node] = true;
    }

    /** Pin the axial degree of freedom at a node. */
    public void fixAxial(int node) {
        fixedZ[node] = true;
    }

    /**
     * Impose plane-strain conditions on an axisymmetric mesh by pinning every axial degree
     * of freedom. With u_z identically zero the axisymmetric formulation reduces exactly to
     * plane strain, which is what makes the closed-form Lame comparison exact rather than
     * approximate.
     */
    public void fixAllAxial() {
        Arrays.fill(fixedZ, true);
    }

    /**
     * Ramp bore pressure from zero to {@code pascals} over {@code rampSeconds} using a
     * smoothstep, then hold. A ramp long compared with the structure's breathing period
     * keeps the response quasi-static; the damping term removes what ringing remains.
     */
    public void setPressureRamp(double pascals, double rampSeconds) {
        this.targetPressure = pascals;
        this.rampDuration = rampSeconds;
    }

    /**
     * Mass-proportional damping, as a fraction of critical for a mode of the given angular
     * frequency. This is dynamic relaxation: a device for reaching a static answer with an
     * explicit code, not a physical model.
     *
     * <p>The frequency is supplied by the caller rather than inferred here, because the
     * right one is a property of the <em>structure</em> and this class only knows the
     * material. Getting it wrong is not a small error in either direction: too low and the
     * ringing never dies, too high and the system is over-damped and creeps toward
     * equilibrium more slowly than an undamped one would ring down. Both look like the same
     * symptom -- a run that has not settled -- which is why {@link #kineticEnergy()} over
     * {@link #strainEnergy()} has to be checked rather than assumed.
     */
    public void setRelaxationDamping(double fractionOfCritical, double angularFrequency) {
        if (angularFrequency <= 0.0) {
            throw new IllegalArgumentException("angular frequency must be positive");
        }
        this.damping = 2.0 * fractionOfCritical * angularFrequency;
    }

    /**
     * Flanagan-Belytschko hourglass stiffness scaling, dimensionless. Zero disables control,
     * which is useful only for demonstrating that the mode really is free without it.
     *
     * <p>Too small and the hourglass mode grows; too large and it stiffens the element in
     * bending, which is the locking that reduced integration exists to avoid. The usual
     * range is 0.01 to 0.10 and the answer should be insensitive across it -- if it is not,
     * the mesh is being carried by hourglass stiffness rather than by physics.
     */
    public void setHourglassCoefficient(double q) {
        if (q < 0.0) throw new IllegalArgumentException("hourglass coefficient must be >= 0");
        this.hourglassCoefficient = q;
    }

    /**
     * Places a rigid frictionless anvil. Enforced by {@link #step()} only: {@link #advance()}
     * drives a prescribed deformation and has no business being overruled by a constraint.
     */
    public void setRigidWall(RigidWall wall) {
        this.wall = wall;
    }

    /**
     * Lets the body touch itself and anything else in the same mesh. See {@link Contact}.
     *
     * <p>Independent of {@link #setRigidWall}: an anvil is a boundary condition and this is a
     * pair of surfaces, and a scene may have both. Null switches it off, and off is the
     * default, because a mesh with no second body in it pays nothing for a surface it will
     * never touch.
     */
    public void setContact(Contact contact) {
        this.contact = contact;
    }

    /** The contact model, or null if none was set. */
    public Contact contact() {
        return contact;
    }

    /**
     * A uniform body acceleration, m/s^2. Earth is {@code setGravity(0, -9.81)}.
     *
     * <p>Off by default, and that is not an oversight: every validation case in this project
     * runs at stresses where a g is nothing. A tube bursting at 20 MPa is carrying six orders
     * of magnitude more than its own weight, and switching gravity on would change the answer
     * in the eighth decimal place while making every gate depend on the direction the model
     * happens to be drawn in. It matters for a sandbox, where things are supposed to fall,
     * and nowhere else here.
     *
     * <p>Applied to the lumped nodal mass, so it is exact for the rigid-body part of the
     * motion and consistent with whatever quadrature the mass came from.
     */
    public void setGravity(double accelerationR, double accelerationZ) {
        this.gravityR = accelerationR;
        this.gravityZ = accelerationZ;
    }

    /**
     * Lets fully softened elements be removed, so that a body can come apart.
     *
     * <h2>What is thrown away and what is not</h2>
     *
     * A deleted element stops contributing internal force and stops constraining the timestep.
     * Its <b>nodes and their mass stay</b>, which is what makes this conservative rather than
     * merely convenient: an element's internal forces sum to zero, so removing them cannot
     * move the total momentum, and keeping the mass means none of it vanishes. A node left
     * with no live element becomes a free particle carrying whatever momentum it had, which is
     * physically what a fragment too small to mesh is.
     *
     * <p>What <em>is</em> thrown away is the elastic energy the element still held. That is
     * reported by {@link #erodedEnergy()} rather than quietly dropped, because it is the one
     * term that stops an energy audit closing across a deletion, and a run that erodes a lot
     * should be able to say how much it lost.
     *
     * <h2>Why this is not a physical model</h2>
     *
     * Nothing in continuum mechanics says a piece of material ceases to exist. Deletion is a
     * numerical device for getting rid of elements that can no longer be integrated -- they
     * carry no stress, their shape functions are meaningless, and they hold the timestep down
     * for nothing. The <em>physics</em> of separation is the softening that got them there,
     * which is regularised by fracture energy and is the same at any mesh size; deletion is
     * only the disposal. That distinction is what keeps the charter's ban on binary failure
     * intact: an element does not switch off, it softens to nothing and is then swept up.
     *
     * <p>It still costs accuracy, and honestly: the deleted volume is gone, so a body that
     * erodes loses section, and a crack opened this way is exactly one element wide however
     * fine the mesh. This is the trade the sandbox licenses and the validation programme does
     * not.
     *
     * @param damageThreshold element damage at which to delete, 0 to 1; 1.0 is "completely
     *                        softened", which is the only value with a physical reading
     */
    public void setErosion(double damageThreshold) {
        if (damage == null) {
            throw new IllegalStateException(
                    "erosion needs a failure model; call setDamage first");
        }
        if (!(damageThreshold > 0.0) || damageThreshold > 1.0) {
            throw new IllegalArgumentException(
                    "the erosion threshold is a damage, so it lies in (0, 1]; was "
                            + damageThreshold);
        }
        this.erosionDamage = damageThreshold;
        if (gone == null) gone = new boolean[mesh.elementCount];
    }

    /**
     * Inflates the mass of the whole mesh so the stable step reaches a target, m/s^2 of
     * honesty traded for wall-clock time.
     *
     * <h2>Why a scene needs it</h2>
     *
     * Gravity works in milliseconds and an explicit solve works in microseconds. A block
     * falling one millimetre takes 14 ms, which at the CFL step of a 1 mm copper mesh is a
     * hundred and sixty thousand steps -- to watch something drop the height of its own
     * skin. Every scene in this project that involves gravity has had to begin already
     * stacked, because settling was unaffordable. That is the constraint this lifts.
     *
     * <h2>Why the simple form is the right one here</h2>
     *
     * Textbook mass scaling is <em>selective</em>: find the elements whose own CFL step is
     * below the target and inflate only those, so the mesh's inertia is disturbed as little as
     * possible. That matters on a mesh with a few small elements among many large ones -- which
     * is what a fitted mesher produces and is exactly what {@link org.neofiz.mesh.Outline}
     * refuses to make. On a lattice every element is the same size and the same material can
     * be given the same factor, so selective scaling and uniform scaling are <b>the same
     * thing</b>. The mesher's blockiness pays for itself twice.
     *
     * <p>Uniformity is also what makes this cheap. Scaling every mass by {@code f} divides
     * every wave speed by {@code sqrt(f)} and multiplies the stable step by {@code sqrt(f)},
     * so the CFL rule needs no per-element wave speed and the hot loop is untouched.
     *
     * <h2>What it costs</h2>
     *
     * Momentum. A body with {@code f} times its mass carries {@code f} times the momentum at
     * the same speed, so anything whose answer depends on inertia -- an impact, a ricochet, a
     * momentum transfer -- is wrong by that factor. Gravity is not: weight scales with mass,
     * so a scaled body still falls at g and still rests at the same load. Use it for settling,
     * stacking and toppling; do not use it for the moment something arrives at speed.
     * {@link #massScale()} reports the factor so a run can say what it bought.
     *
     * @param targetStep the timestep wanted, seconds. A target below the current step is
     *                   refused rather than silently ignored: mass is never removed.
     */
    public void setMassScaling(double targetStep) {
        if (!(targetStep > 0.0)) {
            throw new IllegalArgumentException("the target step must be positive");
        }
        if (targetStep < dt) {
            throw new IllegalArgumentException(
                    "the step is already " + dt + " s, longer than the target " + targetStep
                            + "; mass scaling adds mass and cannot shorten a step");
        }
        final double factor = (targetStep / dt) * (targetStep / dt);
        for (int i = 0; i < mass.length; i++) mass[i] *= factor;
        massScale *= factor;
        dt = targetStep;
        dtPrevious = dt;
        // The reference step is what the contact stiffness is built from, so it has to move
        // with the real one. Leaving it behind would make the penalty frequency sqrt(f) times
        // what the step can integrate, which for a hundredfold scaling is a divergent run.
        dt0 = dt;
    }

    /** How much heavier the mesh has been made, 1 when nothing was scaled. */
    public double massScale() {
        return massScale;
    }

    /** Whether an element has been deleted. False everywhere when erosion is off. */
    public boolean isEroded(int element) {
        return gone != null && gone[element];
    }

    /** How many elements have been deleted over the run. */
    public int erodedElements() {
        return erodedCount;
    }

    /** Elastic energy discarded with deleted elements, joules. The audit's missing term. */
    public double erodedEnergy() {
        return erodedEnergy;
    }

    /**
     * Elements deleted for having turned inside out rather than for having failed.
     *
     * <p>Separate from the damage count because it means something different. An inverted
     * element is a solve that has gone wrong -- the step was too long, or the material was
     * driven somewhere it cannot go -- and deleting it keeps the run alive but does not make
     * the answer right. A handful at the end of a violent scene is ordinary; a lot of them
     * early is a result to distrust.
     */
    public int invertedElements() {
        return invertedCount;
    }

    /**
     * Deletes whatever has failed or inverted since the last step.
     *
     * <p>Runs once per step and only when erosion is on. Two criteria: damage past the
     * threshold, which is failure, and a non-positive current area, which is not a failure but
     * an element that can no longer be integrated at all -- its Jacobian has changed sign, so
     * the volume it reports is negative and every force it assembles has the wrong sign.
     * Keeping one of those is worse than losing it.
     */
    private void erode() {
        for (int e = 0; e < mesh.elementCount; e++) {
            if (gone[e]) continue;
            final boolean failed = damage != null && elementDamage(e) >= erosionDamage;
            final boolean inverted = currentArea(e) <= 0.0;
            if (!failed && !inverted) continue;

            gone[e] = true;
            erodedCount++;
            if (inverted && !failed) invertedCount++;
            for (int g = 0; g < pointsPerElement; g++) {
                erodedEnergy += pointStrainEnergy(e, g);
            }
            if (contact != null) contact.forget(e);
        }
    }

    /** Signed area of an element as it is now, square metres. Negative means inside out. */
    private double currentArea(int element) {
        final int b = element * 4;
        double a = 0.0;
        for (int k = 0; k < 4; k++) {
            final int i = mesh.conn[b + k];
            final int j = mesh.conn[b + (k + 1) % 4];
            a += (mesh.r[i] + ur[i]) * (mesh.z[j] + uz[j])
                    - (mesh.r[j] + ur[j]) * (mesh.z[i] + uz[i]);
        }
        return 0.5 * a;
    }

    /**
     * Gives every element a failure strain and lets it soften past it. See {@link Damage}.
     *
     * <p>The geometry handed over is the <b>reference</b> one, not the current one. By the time
     * softening runs the element has stretched by tens of percent, so the two differ
     * materially, and the reference is the right choice because a fracture energy measured on a
     * real specimen is quoted per unit of <em>original</em> crack area. Using the current area
     * would make G_f mean something that no test reports.
     *
     * <p>The two mean edge vectors are the natural-coordinate directions of the quad,
     * {@code u = ((p1-p0) + (p2-p3))/2} and {@code v = ((p3-p0) + (p2-p1))/2}, which for a
     * parallelogram are its sides exactly and for a general quad are the sides of the
     * parallelogram with the same area and the same diagonals.
     *
     * @param fractureEnergy energy per unit crack area, J/m^2
     * @param failureStrain  plastic strain at onset, one per <em>element</em>. Under full
     *                       integration all four Gauss points inherit the element's value: the
     *                       defect field's correlation length is far larger than an element, so
     *                       sub-element variation is not something it resolves.
     */
    public void setDamage(double fractureEnergy, double[] failureStrain) {
        final double[] gf = new double[mesh.elementCount];
        java.util.Arrays.fill(gf, fractureEnergy);
        setDamage(gf, failureStrain);
    }

    /**
     * The per-element form, for a mesh holding more than one substance.
     *
     * <p>G_f is the property that separates materials most sharply -- a window pane and a
     * pressure-vessel steel are four orders of magnitude apart -- so a mixed mesh with one
     * fracture energy on it is not a mixed mesh in the way that matters.
     *
     * @param fractureEnergy energy per unit crack area, J/m^2, one per element
     * @param failureStrain  plastic strain at onset, one per element
     */
    public void setDamage(double[] fractureEnergy, double[] failureStrain) {
        if (fractureEnergy.length != mesh.elementCount) {
            throw new IllegalArgumentException(
                    "fracture energy has " + fractureEnergy.length + " entries for a mesh of "
                            + mesh.elementCount + " elements");
        }
        if (failureStrain.length != mesh.elementCount) {
            throw new IllegalArgumentException(
                    "failure strain has " + failureStrain.length + " entries for a mesh of "
                            + mesh.elementCount + " elements");
        }
        // The return map divides by 2 mu + (2/3)(H_soft + H_kin) and Damage clamps H_soft at
        // -E, so this is the worst denominator the solve can ever see. It is positive for every
        // nu < 0.5, which Material already enforces -- so this cannot fire, and it is here to
        // say out loud which two constraints are holding each other up. Checked per material,
        // because on a mixed mesh it is the softest one that would fail first and the mesh
        // average would hide it.
        for (Material m : materials.distinct()) {
            final double worst = 2.0 * m.mu() + (2.0 / 3.0)
                    * (m.kinematicHardening() - m.youngsModulus());
            if (!(worst > 0.0)) {
                throw new IllegalArgumentException(
                        "a fully clamped softening modulus would make the return map singular "
                                + "for " + m.name() + ": 2mu + (2/3)(H_kin - E) = " + worst
                                + " <= 0");
            }
        }

        // Every point of an element is handed that element's whole geometry, so all of them
        // arrive at the same band width. See Damage: the band is one element wide because that
        // is the narrowest strain discontinuity a bilinear quad can represent, and the point
        // volumes already sum to the element volume, so the energy comes out the same under
        // either quadrature rule.
        final int points = mesh.elementCount * pointsPerElement;
        final double[] perPoint = new double[points];
        final double[] gfOf = new double[points];
        final double[] eOf = new double[points];
        final double[] geometry = new double[Damage.GEOMETRY_STRIDE * points];
        for (int e = 0; e < mesh.elementCount; e++) {
            final int b = e * 4;
            final int n0 = mesh.conn[b], n1 = mesh.conn[b + 1];
            final int n2 = mesh.conn[b + 2], n3 = mesh.conn[b + 3];
            final double uR = 0.5 * (mesh.r[n1] - mesh.r[n0] + mesh.r[n2] - mesh.r[n3]);
            final double uZ = 0.5 * (mesh.z[n1] - mesh.z[n0] + mesh.z[n2] - mesh.z[n3]);
            final double vR = 0.5 * (mesh.r[n3] - mesh.r[n0] + mesh.r[n2] - mesh.r[n1]);
            final double vZ = 0.5 * (mesh.z[n3] - mesh.z[n0] + mesh.z[n2] - mesh.z[n1]);
            final double area = Math.abs(mesh.signedArea(e));

            for (int g = 0; g < pointsPerElement; g++) {
                final int p = e * pointsPerElement + g;
                perPoint[p] = failureStrain[e];
                gfOf[p] = fractureEnergy[e];
                eOf[p] = materials.at(e).youngsModulus();
                final int q = Damage.GEOMETRY_STRIDE * p;
                geometry[q] = area;
                geometry[q + 1] = uR;
                geometry[q + 2] = uZ;
                geometry[q + 3] = vR;
                geometry[q + 4] = vZ;
            }
        }
        this.damage = new Damage(gfOf, eOf, perPoint, geometry);
    }

    /** The damage state, or null if none was set. */
    public Damage damage() {
        return damage;
    }

    /** Prescribe a nodal velocity. Also the diagnostic entry point, with {@link #advance()}. */
    public void setVelocity(int node, double velocityR, double velocityZ) {
        vr[node] = velocityR;
        vz[node] = velocityZ;
    }

    /** Give every node the same velocity, as a projectile in flight has. */
    public void setUniformVelocity(double velocityR, double velocityZ) {
        Arrays.fill(vr, velocityR);
        Arrays.fill(vz, velocityZ);
    }

    public Integration integration() {
        return integration;
    }

    /**
     * Which configuration this solver writes equilibrium in. Anything computing a global
     * equilibrium check from the stress field has to ask: under small strain the internal
     * forces are assembled on the reference mesh, so a check that measures the geometry where
     * the nodes have actually moved to is comparing two different configurations and will
     * report an imbalance that is not there.
     */
    public Kinematics kinematics() {
        return kinematics;
    }

    /** Bore pressure currently applied, Pa -- the value in the last force assembly. */
    public double borePressure() {
        return pressureAt(time);
    }

    public double timestep() {
        return dt;
    }

    public double time() {
        return time;
    }

    public long steps() {
        return stepCount;
    }

    public double[] radialDisplacement() {
        return ur;
    }

    public double[] axialDisplacement() {
        return uz;
    }

    /**
     * Kinetic energy with the velocity brought forward to the same instant as the
     * displacement, joules.
     *
     * <p>{@link #kineticEnergy()} is half a step ahead of everything it is ever compared
     * against. Central difference carries velocity at {@code t + dt/2} and displacement at
     * {@code t + dt}, so an audit that adds kinetic to stored energy is adding two numbers
     * from two different instants, and the mismatch is
     *
     * <pre>  KE(t + dt/2) - KE(t + dt) ~ -(dt/2) * d(KE)/dt</pre>
     *
     * which is <b>first order</b> in the step, not second. On a body in rigid translation it
     * is nothing, because the kinetic energy is not changing. On a vibrating body it is a
     * swing of about {@code omega*dt/2} of the vibrational energy -- 1.4 % on a bar meshed at
     * the usual CFL safety factor, with no error anywhere in the solve. That is large enough
     * to be mistaken for a defect in whatever was added last, and was: a collision turns
     * translation into vibration, so the total-energy audit steps up as the bodies meet and
     * the contact model looks like the source.
     *
     * <p>The correction is exact and free. The velocity half a step further on is
     * {@code v + dt*f/m} with the force already assembled at the current displacement, so the
     * centred velocity is {@code v + (dt/2)*f/m} and no extra state is needed. Constrained
     * components are left alone: a pinned node has no velocity to centre, and a node the wall
     * is holding carries a force the wall cancels rather than one that will accelerate it.
     */
    public double centredKineticEnergy() {
        double ke = 0.0;
        final double h = 0.5 * dt;
        for (int i = 0; i < mesh.nodeCount; i++) {
            final double invM = 1.0 / mass[i];
            final boolean held = wall != null && onWall[i] && wall.holds(fz[i]);
            final double cvr = fixedR[i] ? vr[i] : vr[i] + h * fr[i] * invM;
            final double cvz = fixedZ[i] || held ? vz[i] : vz[i] + h * fz[i] * invM;
            ke += 0.5 * mass[i] * (cvr * cvr + cvz * cvz);
        }
        return ke;
    }

    /**
     * Mean speed of an element's four corners, m/s.
     *
     * <p>The field that makes a scene of several bodies legible. Stress and plastic strain say
     * what a body is going through; speed says which piece is going where, which is the
     * question a picture of a collapse is actually being asked.
     */
    public double elementSpeed(int element) {
        double sum = 0.0;
        for (int k = 0; k < 4; k++) {
            final int n = mesh.conn[element * 4 + k];
            sum += Math.sqrt(vr[n] * vr[n] + vz[n] * vz[n]);
        }
        return 0.25 * sum;
    }

    /** Nodal velocity, radial component. Live, not a copy. */
    public double[] velocityR() {
        return vr;
    }

    /** Nodal velocity, axial component. Live, not a copy. */
    public double[] velocityZ() {
        return vz;
    }

    /** Net nodal force from the last assembly, radial component. */
    public double[] forceR() {
        return fr;
    }

    /** Net nodal force from the last assembly, axial component. */
    public double[] forceZ() {
        return fz;
    }

    /** Accumulated equivalent plastic strain, one per Gauss point. */
    public double[] plasticStrain() {
        return epsP;
    }

    /** Stored stress, stride 4 per Gauss point: {rr, zz, tt, rz}. */
    public double[] stress() {
        return sig;
    }

    // ------------------------------------------------------------------ run

    public void run(int steps) {
        for (int i = 0; i < steps; i++) step();
    }

    public void step() {
        // Forces at t = 0, from the stored stress and the initial pressure. Passing a zero
        // strain interval makes this the same code path with the material update turned into
        // a no-op, rather than a second assembly routine that can drift out of step with the
        // real one.
        if (stepCount == 0) advanceMaterialAndAssemble(0.0);

        // Central difference with a timestep that may have changed: the velocity is centred
        // between two steps, so it advances over their mean. Under small strain the geometry
        // is frozen, dt never moves, and this reduces exactly to the constant-step form.
        final double dtAvg = 0.5 * (dtPrevious + dt);
        final double half = 0.5 * damping * dtAvg;
        nodalDecay = (1.0 - half) / (1.0 + half);
        nodalGain = dtAvg / (1.0 + half);
        nodalDtAvg = dtAvg;

        if (wall != null) wall.beginStep();

        // Threaded only without a wall. Contact is not thread-safe here and would not be
        // worth making so: the wall carries running totals -- impulse, destroyed kinetic
        // energy, frictional dissipation -- and those are sums, so splitting them across
        // workers would make the audits depend on the schedule. Since the schedule is now
        // dynamic, that is not merely a different answer but a different answer each run,
        // which is exactly the property the rest of this class is built to avoid. The audits
        // are what caught the one real bug in the contact work, so they are not something to
        // trade for speed on a case whose mesh is two orders of magnitude smaller than the
        // one the cost model is about.
        if (pool != null && wall == null) {
            pool.run(mesh.nodeCount, nodalKernel);
        } else {
            nodalUpdate(0, mesh.nodeCount);
        }

        time += dt;
        stepCount++;

        // Before the assembly, so a deleted element never contributes again -- and after the
        // nodal update, so the geometry the inversion test reads is the one the step produced.
        if (gone != null) erode();

        advanceMaterialAndAssemble(dt);

        // Elements that have compressed carry a shorter wave transit, so the stable step
        // shrinks with them. Holding dt fixed while the mesh crushes is how an explicit run
        // goes unstable late, long after the setup that looked fine.
        if (kinematics.isFinite()) {
            dtPrevious = dt;
            // sqrt(massScale) because inflating every mass by f divides every wave speed by
            // sqrt(f). One factor here is the whole of mass scaling in the hot path.
            dt = cflSafety * Math.sqrt(minEdgeSquared * massScale) / fastestWave;
        }
    }

    /**
     * Central-difference velocity and displacement update over a range of nodes, including
     * the contact branches. One method rather than a fast path beside a general one: the
     * range is the only thing threading changes, so there is no second copy to drift.
     */
    private void nodalUpdate(int from, int to) {
        final double decay = nodalDecay;
        final double gain = nodalGain;
        final double dtAvg = nodalDtAvg;

        for (int i = from; i < to; i++) {
            final double invM = 1.0 / mass[i];

            // The contact state is settled once, before either component is touched. The
            // axial branch below rewrites onWall, so asking again afterwards would apply
            // friction using this step's answer for the normal direction and the previous
            // step's for the tangential one.
            final boolean resting = wall != null && onWall[i] && wall.holds(fz[i]);

            if (fixedR[i]) {
                vr[i] = 0.0;
            } else if (resting && wall.hasFriction()) {
                vr[i] = wall.slide(fr[i], fz[i], vr[i], decay, gain, mass[i], dt);
                ur[i] += vr[i] * dt;
            } else {
                vr[i] = decay * vr[i] + gain * fr[i] * invM;
                ur[i] += vr[i] * dt;
            }
            if (fixedZ[i]) {
                vz[i] = 0.0;
            } else if (resting) {
                // Resting on the anvil with the material still pressing into it. The wall
                // supplies a reaction that cancels fz exactly, so this node does not
                // accelerate and does not move -- and because it does not move, the reaction
                // does no work. Letting the force act and then zeroing the velocity instead
                // would be an energy leak that no amount of refinement removes; see RigidWall.
                //
                // vz is already exactly zero here, established by the arrival below and never
                // disturbed while held. It is deliberately not re-zeroed: an assignment would
                // paper over any momentum this branch failed to account for, which is exactly
                // what the balance audit exists to catch.
                //
                // dtAvg, not dt: the reaction cancels the force over the same interval the
                // velocity update would have applied it over, and under adaptive CFL those
                // differ. Using dt leaves a momentum error invisible at small strain -- where
                // dt never moves -- that grows to 1e-3 on a run whose mesh is crushing.
                wall.hold(fz[i], dtAvg);
            } else {
                final double trial = decay * vz[i] + gain * fz[i] * invM;
                final double zc = mesh.z[i] + uz[i];
                if (wall != null && wall.crosses(zc, trial, dt)) {
                    // Arrival, and it has to be one event: land exactly on the plane and stop
                    // there. Landing with the velocity that reaches the plane and stopping on
                    // the following step instead splits it in two, and the momentum destroyed
                    // by the second half belongs to no impulse -- a leak that is silent at
                    // first impact, where a node already touching the plane lands with zero
                    // velocity anyway, and appears only once the contact patch starts
                    // spreading to nodes arriving from above.
                    wall.arrive(trial, mass[i]);
                    uz[i] += wall.plane() - zc;
                    vz[i] = 0.0;
                    onWall[i] = true;
                } else {
                    vz[i] = trial;
                    uz[i] += vz[i] * dt;
                    if (wall != null) onWall[i] = false;
                }
            }
        }
    }

    /**
     * Applies the strain increment implied by the current velocities, updates the material
     * state, and assembles forces -- without touching the velocities. Diagnostic entry point
     * for driving a prescribed deformation; {@link #step()} runs the same code.
     */
    public void advance() {
        for (int i = 0; i < mesh.nodeCount; i++) {
            if (!fixedR[i]) ur[i] += vr[i] * dt;
            if (!fixedZ[i]) uz[i] += vz[i] * dt;
        }
        advanceMaterialAndAssemble(dt);
    }

    /**
     * @param strainDt interval over which the current velocities act. Zero assembles forces
     *                 from the stored stress without advancing the material.
     */
    private void advanceMaterialAndAssemble(double strainDt) {
        kernelStrainDt = strainDt;
        kernelPressure = pressureAt(time);
        Arrays.fill(minEdgeLocal, Double.MAX_VALUE);

        final int sources = edgeSlotBase + mesh.pressureEdges.length / 2;
        if (pool == null) {
            assemble(0, sources, 0);
            gather(0, mesh.nodeCount);
        } else {
            pool.run(sources, assembleKernel);
            pool.run(mesh.nodeCount, gatherKernel);
        }

        // After the gather, because contact adds to the summed nodal force rather than
        // contributing a slot to it: its pairs are discovered per step, so there is no fixed
        // source list for the gather to index. Before the step that integrates them, so the
        // penalty acts over the same interval the material does.
        if (contact != null) {
            // The search is threaded and the accumulation is not; see Contact. It is the most
            // expensive thing in a scene with several bodies in it -- two and a half times
            // the element kernel on the benchmark -- and almost all of that is the search.
            final int slaves = contact.begin(ur, uz, dt, dt0);
            if (slaves > 0) {
                if (pool != null) pool.run(slaves, contactKernel);
                else contact.search(0, slaves);
                contact.accumulate(vr, vz, mass, fr, fz);
            }
        }

        minEdgeSquared = Double.MAX_VALUE;
        for (int w = 0; w < minEdgeLocal.length; w += PAD) {
            minEdgeSquared = Math.min(minEdgeSquared, minEdgeLocal[w]);
        }
    }

    /**
     * One chunk of the assembly. Elements and pressure edges share a single index space so
     * that they share a single parallel region, and a chunk that straddles the boundary does
     * part of each.
     *
     * <p>A worker may be handed several chunks, so the shortest edge accumulates into its slot
     * rather than overwriting it -- and minimum being exactly associative is why the slots can
     * be combined afterwards in any order without moving the timestep by a bit.
     */
    private void assemble(int from, int to, int worker) {
        final int ne = edgeSlotBase;
        if (from < ne) {
            final int end = Math.min(to, ne);
            final double edge = integration == Integration.FULL
                    ? fullIntegration(kernelStrainDt, from, end)
                    : reducedIntegration(kernelStrainDt, from, end);
            final int slot = worker * PAD;
            if (edge < minEdgeLocal[slot]) minEdgeLocal[slot] = edge;
        }
        if (to > ne) applyPressure(kernelPressure, Math.max(from, ne) - ne, to - ne);
    }

    /**
     * Sums each node's contributions. Replaces the scatter that used to end the element loop
     * and the one that used to end the pressure load, and is the reason the assembly can be
     * threaded without an atomic anywhere.
     *
     * <p>Assigns rather than accumulates, so no zeroing pass is needed -- and so a node that
     * belongs to no source reads exactly zero rather than whatever was there last step.
     */
    private void gather(int from, int to) {
        final double[] c = contribution;
        final int[] slots = nodeSlot;
        for (int n = from; n < to; n++) {
            double sr = 0.0, sz = 0.0;
            final int end = nodeSlotStart[n + 1];
            for (int k = nodeSlotStart[n]; k < end; k++) {
                final int slot = slots[k];
                sr += c[slot];
                sz += c[slot + 4];
            }
            // Gravity is added here rather than contributed as a source, because it is the
            // one force that belongs to a node rather than to anything the node is part of.
            // Adding it after the gather also keeps it out of the parallel assembly, where it
            // would have needed a slot of its own for no benefit.
            fr[n] = sr + mass[n] * gravityR;
            fz[n] = sz + mass[n] * gravityZ;
        }
    }

    /**
     * Radial coordinate of a node in the configuration the kernel works in.
     *
     * <p>The displacement has already been advanced when the kernel runs, so a shift of
     * {@code -dt/2} steps back to the mid-step configuration. That is where the velocity
     * gradient has to be evaluated: it is the one configuration in which rigid body motion
     * yields an exactly skew velocity gradient, and therefore exactly zero strain. Evaluating
     * at either end instead makes rotation leak strain at second order per step.
     */
    private double cr(int node, double shift) {
        return kinematics.isFinite() ? mesh.r[node] + ur[node] + shift * vr[node] : mesh.r[node];
    }

    private double cz(int node, double shift) {
        return kinematics.isFinite() ? mesh.z[node] + uz[node] + shift * vz[node] : mesh.z[node];
    }

    /**
     * Shortest current edge of an element against a running minimum, for the adaptive CFL
     * limit. Squared, so the whole mesh costs one square root per step rather than four per
     * element. Carried as a return value rather than written to a field because each worker
     * keeps its own, and minimum is exactly associative in floating point so combining them
     * afterwards is order-independent.
     */
    private static double shortestEdge(double running,
                                       double r0, double r1, double r2, double r3,
                                       double z0, double z1, double z2, double z3) {
        return Math.min(running, Math.min(
                Math.min(sq(r1 - r0, z1 - z0), sq(r2 - r1, z2 - z1)),
                Math.min(sq(r3 - r2, z3 - z2), sq(r0 - r3, z0 - z3))));
    }

    private static double sq(double a, double b) {
        return a * a + b * b;
    }

    private double pressureAt(double t) {
        if (rampDuration <= 0.0) return targetPressure;
        double s = Math.min(1.0, t / rampDuration);
        return targetPressure * s * s * (3.0 - 2.0 * s);
    }

    // ------------------------------------------------------------------ kernels

    /**
     * Full 2x2 integration. Four stress evaluations and four sets of history per element.
     *
     * <p>Accurate and free of spurious modes, but it volumetrically locks as the material
     * approaches incompressibility -- which is where developed plastic flow lives. Kept
     * because it is the reference the cheap rule is validated against, not because it is the
     * rule to run production in.
     */
    private double fullIntegration(double strainDt, int from, int to) {
        final int[] conn = mesh.conn;
        final double[] ef = contribution;
        final boolean axi = formulation.isAxisymmetric();
        final boolean finite = kinematics.isFinite();
        final double mid = -0.5 * strainDt;   // step back to the mid-step configuration
        final Damage dmg = damage;
        double minEdge = Double.MAX_VALUE;

        for (int e = from; e < to; e++) {
            final int b = e * 4;
            // A deleted element contributes nothing, and has to write that nothing rather
            // than skip: the gather sums fixed slots, so a stale force left in one would go
            // on pushing its nodes for the rest of the run. It also stops being counted in
            // the shortest-edge search, which is most of the point -- an element crushed to
            // nothing would otherwise hold the whole scene to its timestep.
            if (gone != null && gone[e]) {
                for (int k = 0; k < 8; k++) ef[b * 2 + k] = 0.0;
                continue;
            }
            final J2.Flow fl = flowOf[e];
            final double lambda = lambdaOf[e], mu = muOf[e];
            final int n0 = conn[b], n1 = conn[b + 1], n2 = conn[b + 2], n3 = conn[b + 3];

            final double r0 = cr(n0, mid), r1 = cr(n1, mid), r2 = cr(n2, mid), r3 = cr(n3, mid);
            final double z0 = cz(n0, mid), z1 = cz(n1, mid), z2 = cz(n2, mid), z3 = cz(n3, mid);
            // Displacement increment over this step, from the mid-step velocity.
            final double a0 = vr[n0] * strainDt, a1 = vr[n1] * strainDt;
            final double a2 = vr[n2] * strainDt, a3 = vr[n3] * strainDt;
            final double c0 = vz[n0] * strainDt, c1 = vz[n1] * strainDt;
            final double c2 = vz[n2] * strainDt, c3 = vz[n3] * strainDt;

            if (finite) minEdge = shortestEdge(minEdge, r0, r1, r2, r3, z0, z1, z2, z3);

            double g0r = 0, g1r = 0, g2r = 0, g3r = 0;
            double g0z = 0, g1z = 0, g2z = 0, g3z = 0;

            for (int g = 0; g < 4; g++) {
                final double xi = GP_XI[g], eta = GP_ETA[g];
                final double xm = 1.0 - xi, xp = 1.0 + xi;
                final double em = 1.0 - eta, ep = 1.0 + eta;

                final double N0 = 0.25 * xm * em, N1 = 0.25 * xp * em;
                final double N2 = 0.25 * xp * ep, N3 = 0.25 * xm * ep;

                final double dN0x = -0.25 * em, dN1x = 0.25 * em, dN2x = 0.25 * ep, dN3x = -0.25 * ep;
                final double dN0e = -0.25 * xm, dN1e = -0.25 * xp, dN2e = 0.25 * xp, dN3e = 0.25 * xm;

                final double J00 = dN0x * r0 + dN1x * r1 + dN2x * r2 + dN3x * r3;
                final double J01 = dN0x * z0 + dN1x * z1 + dN2x * z2 + dN3x * z3;
                final double J10 = dN0e * r0 + dN1e * r1 + dN2e * r2 + dN3e * r3;
                final double J11 = dN0e * z0 + dN1e * z1 + dN2e * z2 + dN3e * z3;
                final double det = J00 * J11 - J01 * J10;
                final double inv = 1.0 / det;

                final double d0r = (J11 * dN0x - J01 * dN0e) * inv;
                final double d1r = (J11 * dN1x - J01 * dN1e) * inv;
                final double d2r = (J11 * dN2x - J01 * dN2e) * inv;
                final double d3r = (J11 * dN3x - J01 * dN3e) * inv;

                final double d0z = (-J10 * dN0x + J00 * dN0e) * inv;
                final double d1z = (-J10 * dN1x + J00 * dN1e) * inv;
                final double d2z = (-J10 * dN2x + J00 * dN2e) * inv;
                final double d3z = (-J10 * dN3x + J00 * dN3e) * inv;

                final double dEpsR = d0r * a0 + d1r * a1 + d2r * a2 + d3r * a3;
                final double dEpsZ = d0z * c0 + d1z * c1 + d2z * c2 + d3z * c3;
                final double dRdZ = d0z * a0 + d1z * a1 + d2z * a2 + d3z * a3;
                final double dZdR = d0r * c0 + d1r * c1 + d2r * c2 + d3r * c3;
                final double dGamRZ = dRdZ + dZdR;

                double rg = 0.0, dEpsT = 0.0;
                if (axi) {
                    rg = N0 * r0 + N1 * r1 + N2 * r2 + N3 * r3;
                    dEpsT = (N0 * a0 + N1 * a1 + N2 * a2 + N3 * a3) / rg;
                }

                final int p = e * 4 + g;

                // Co-rotate before the constitutive update. The stored components describe a
                // state that has physically turned; rotating them first is what makes the
                // increment objective.
                if (finite) {
                    final double spinDt = 0.5 * (dZdR - dRdZ);   // the spin W_zr times dt
                    if (spinDt != 0.0) {
                        final double cos = Corotational.cos(spinDt), sin = Corotational.sin(spinDt);
                        Corotational.rotate(sig, p, cos, sin);
                        Corotational.rotate(back, p, cos, sin);
                    }
                }

                J2.update(sig, back, epsP, wPlastic, p, dEpsR, dEpsZ, dEpsT, dGamRZ,
                        lambda, mu, fl, strainDt, dmg);

                final int s = J2.COMPONENTS * p;
                final double sR = sig[s], sZ = sig[s + 1], sT = sig[s + 2], sRZ = sig[s + 3];
                final double dV = axi ? 2.0 * Math.PI * rg * det : thickness * det;

                if (axi) {
                    final double h = sT / rg * dV;
                    g0r -= (d0r * sR + d0z * sRZ) * dV + N0 * h;
                    g1r -= (d1r * sR + d1z * sRZ) * dV + N1 * h;
                    g2r -= (d2r * sR + d2z * sRZ) * dV + N2 * h;
                    g3r -= (d3r * sR + d3z * sRZ) * dV + N3 * h;
                } else {
                    g0r -= (d0r * sR + d0z * sRZ) * dV;
                    g1r -= (d1r * sR + d1z * sRZ) * dV;
                    g2r -= (d2r * sR + d2z * sRZ) * dV;
                    g3r -= (d3r * sR + d3z * sRZ) * dV;
                }
                g0z -= (d0z * sZ + d0r * sRZ) * dV;
                g1z -= (d1z * sZ + d1r * sRZ) * dV;
                g2z -= (d2z * sZ + d2r * sRZ) * dV;
                g3z -= (d3z * sZ + d3r * sRZ) * dV;
            }

            ef[b * 2] = g0r; ef[b * 2 + 1] = g1r; ef[b * 2 + 2] = g2r; ef[b * 2 + 3] = g3r;
            ef[b * 2 + 4] = g0z; ef[b * 2 + 5] = g1z; ef[b * 2 + 6] = g2z; ef[b * 2 + 7] = g3z;
        }
        return minEdge;
    }

    /**
     * One-point integration with Flanagan-Belytschko hourglass control. The same physics as
     * {@link #fullIntegration} at a quarter of the stress work and a quarter of the
     * history storage.
     *
     * <p><b>Why control is needed.</b> A single sample at the centroid sees only the
     * constant part of the strain field. The hourglass mode -- nodal displacements
     * alternating in sign around the element -- has exactly zero strain there, so it
     * produces exactly zero internal force and nothing resists it. Left alone it grows until
     * the mesh is a checkerboard of folded elements, and because the mode carries no energy
     * the usual energy audits do not flag it.
     *
     * <p><b>Why the shape vector is not just the alternating pattern.</b> On a distorted
     * element the raw pattern h = (1,-1,1,-1) is contaminated by genuine linear deformation,
     * so resisting it would resist real strain. Flanagan and Belytschko project that
     * contamination out:
     *
     * <pre>  gamma_i = h_i - (h . r) b_i^r - (h . z) b_i^z</pre>
     *
     * where b is the centroid gradient operator. The result is orthogonal to every linear
     * displacement field, translations and constant strains alike, so the hourglass force
     * vanishes identically on anything the element is supposed to represent. That is the
     * property that makes this a correction rather than a fudge, and it is asserted in the
     * tests to machine precision rather than to a tolerance.
     *
     * <p><b>Stiffness form, not viscous.</b> The original viscous form resists hourglass
     * <em>velocity</em>. That is fine for shocks and useless here: dynamic relaxation drives
     * velocity to zero, at which point a viscous term stops resisting and any hourglass
     * <em>displacement</em> already accumulated is frozen in. The stiffness form resists
     * displacement, so it holds at equilibrium, and being conservative it has an energy that
     * can be audited against strain energy.
     */
    private double reducedIntegration(double strainDt, int from, int to) {
        final int[] conn = mesh.conn;
        final double[] ef = contribution;
        final boolean axi = formulation.isAxisymmetric();
        final boolean finite = kinematics.isFinite();
        final double mid = -0.5 * strainDt;   // step back to the mid-step configuration
        final double qhg = hourglassCoefficient;
        final Damage dmg = damage;
        double minEdge = Double.MAX_VALUE;

        for (int e = from; e < to; e++) {
            final int b = e * 4;
            // A deleted element contributes nothing, and has to write that nothing rather
            // than skip: the gather sums fixed slots, so a stale force left in one would go
            // on pushing its nodes for the rest of the run.
            if (gone != null && gone[e]) {
                for (int k = 0; k < 8; k++) ef[b * 2 + k] = 0.0;
                continue;
            }
            final J2.Flow fl = flowOf[e];
            final double lambda = lambdaOf[e], mu = muOf[e];
            final int n0 = conn[b], n1 = conn[b + 1], n2 = conn[b + 2], n3 = conn[b + 3];

            final double r0 = cr(n0, mid), r1 = cr(n1, mid), r2 = cr(n2, mid), r3 = cr(n3, mid);
            final double z0 = cz(n0, mid), z1 = cz(n1, mid), z2 = cz(n2, mid), z3 = cz(n3, mid);
            final double a0 = vr[n0] * strainDt, a1 = vr[n1] * strainDt;
            final double a2 = vr[n2] * strainDt, a3 = vr[n3] * strainDt;
            final double c0 = vz[n0] * strainDt, c1 = vz[n1] * strainDt;
            final double c2 = vz[n2] * strainDt, c3 = vz[n3] * strainDt;

            if (finite) minEdge = shortestEdge(minEdge, r0, r1, r2, r3, z0, z1, z2, z3);

            // Centroid: N_i = 1/4, dN/dxi = (-1,1,1,-1)/4, dN/deta = (-1,-1,1,1)/4.
            final double J00 = 0.25 * (-r0 + r1 + r2 - r3);
            final double J01 = 0.25 * (-z0 + z1 + z2 - z3);
            final double J10 = 0.25 * (-r0 - r1 + r2 + r3);
            final double J11 = 0.25 * (-z0 - z1 + z2 + z3);
            final double det = J00 * J11 - J01 * J10;
            final double inv = 1.0 / det;

            final double d0r = 0.25 * (J01 - J11) * inv;
            final double d1r = 0.25 * (J01 + J11) * inv;
            final double d2r = -d0r;
            final double d3r = -d1r;

            final double d0z = 0.25 * (J10 - J00) * inv;
            final double d1z = -0.25 * (J10 + J00) * inv;
            final double d2z = -d0z;
            final double d3z = -d1z;

            final double dEpsR = d0r * a0 + d1r * a1 + d2r * a2 + d3r * a3;
            final double dEpsZ = d0z * c0 + d1z * c1 + d2z * c2 + d3z * c3;
            final double dRdZ = d0z * a0 + d1z * a1 + d2z * a2 + d3z * a3;
            final double dZdR = d0r * c0 + d1r * c1 + d2r * c2 + d3r * c3;
            final double dGamRZ = dRdZ + dZdR;

            double rg = 0.0, dEpsT = 0.0;
            if (axi) {
                rg = 0.25 * (r0 + r1 + r2 + r3);
                dEpsT = 0.25 * (a0 + a1 + a2 + a3) / rg;
            }

            // Co-rotate before the constitutive update, so the increment is objective.
            if (finite) {
                final double spinDt = 0.5 * (dZdR - dRdZ);   // the spin W_zr times dt
                if (spinDt != 0.0) {
                    final double cos = Corotational.cos(spinDt), sin = Corotational.sin(spinDt);
                    Corotational.rotate(sig, e, cos, sin);
                    Corotational.rotate(back, e, cos, sin);
                }
            }

            J2.update(sig, back, epsP, wPlastic, e, dEpsR, dEpsZ, dEpsT, dGamRZ,
                    lambda, mu, fl, strainDt, dmg);

            final int s = J2.COMPONENTS * e;
            final double sR = sig[s], sZ = sig[s + 1], sT = sig[s + 2], sRZ = sig[s + 3];
            final double dV = W_REDUCED * (axi ? 2.0 * Math.PI * rg * det : thickness * det);

            double g0r, g1r, g2r, g3r;
            if (axi) {
                final double h = sT / rg * dV * 0.25;   // N_i = 1/4 at the centroid
                g0r = -((d0r * sR + d0z * sRZ) * dV + h);
                g1r = -((d1r * sR + d1z * sRZ) * dV + h);
                g2r = -((d2r * sR + d2z * sRZ) * dV + h);
                g3r = -((d3r * sR + d3z * sRZ) * dV + h);
            } else {
                g0r = -(d0r * sR + d0z * sRZ) * dV;
                g1r = -(d1r * sR + d1z * sRZ) * dV;
                g2r = -(d2r * sR + d2z * sRZ) * dV;
                g3r = -(d3r * sR + d3z * sRZ) * dV;
            }
            double g0z = -(d0z * sZ + d0r * sRZ) * dV;
            double g1z = -(d1z * sZ + d1r * sRZ) * dV;
            double g2z = -(d2z * sZ + d2r * sRZ) * dV;
            double g3z = -(d3z * sZ + d3r * sRZ) * dV;

            if (qhg > 0.0) {
                final double hr = r0 - r1 + r2 - r3;
                final double hz = z0 - z1 + z2 - z3;

                final double y0 = H[0] - hr * d0r - hz * d0z;
                final double y1 = H[1] - hr * d1r - hz * d1z;
                final double y2 = H[2] - hr * d2r - hz * d2z;
                final double y3 = H[3] - hr * d3r - hz * d3z;

                // Shear modulus scaling: the hourglass mode is shear-like, so mu is the
                // right stiffness to borrow. Using the bulk modulus over-stiffens it.
                final double k = qhg * mu * dV
                        * (d0r * d0r + d1r * d1r + d2r * d2r + d3r * d3r
                        + d0z * d0z + d1z * d1z + d2z * d2z + d3z * d3z);

                // Accumulated as a force rather than recomputed from total displacement.
                // Under finite strain the shape vector is rebuilt on a mesh that has moved,
                // so there is no fixed total to multiply; at small strain the two agree.
                hgR[e] += k * (y0 * a0 + y1 * a1 + y2 * a2 + y3 * a3);
                hgZ[e] += k * (y0 * c0 + y1 * c1 + y2 * c2 + y3 * c3);
                final double qr = hgR[e], qz = hgZ[e];

                g0r -= y0 * qr; g1r -= y1 * qr; g2r -= y2 * qr; g3r -= y3 * qr;
                g0z -= y0 * qz; g1z -= y1 * qz; g2z -= y2 * qz; g3z -= y3 * qz;

                // Stiffness k and the shape vector are both needed to convert the stored
                // force back into an energy; cache k so the audit does not rebuild them.
                hourglassStiffness[e] = k;
            }

            ef[b * 2] = g0r; ef[b * 2 + 1] = g1r; ef[b * 2 + 2] = g2r; ef[b * 2 + 3] = g3r;
            ef[b * 2 + 4] = g0z; ef[b * 2 + 5] = g1z; ef[b * 2 + 6] = g2z; ef[b * 2 + 7] = g3z;
        }
        return minEdge;
    }

    /**
     * Consistent-edge pressure load. For an axisymmetric edge the shape functions are
     * weighted by 2*pi*r, which puts more of the load on the larger-radius node; the
     * closed forms below are the exact integrals of N_i * 2*pi*r along a straight edge.
     */
    private void applyPressure(double p, int from, int to) {
        final int[] edges = mesh.pressureEdges;
        final double[] c = contribution;

        for (int e = from; e < to; e++) {
            final int k = 2 * e;
            final int slot = (edgeSlotBase + e) * 8;
            final int a = edges[k], bnode = edges[k + 1];
            // A follower load under finite strain: pressure acts on the surface where it
            // now is, at the area it now has, along the normal it now has. A bore that has
            // opened presents more area to the same pressure, and that feedback is part of
            // why a pressurised vessel goes unstable rather than merely stretching.
            final double ra = cr(a, 0.0), rb = cr(bnode, 0.0);
            final double tr = rb - ra;
            final double tz = cz(bnode, 0.0) - cz(a, 0.0);
            final double len = Math.sqrt(tr * tr + tz * tz);
            if (len == 0.0) {
                // A degenerate edge carries no load, but its slot is read by the gather every
                // step, so it has to be cleared rather than skipped.
                c[slot] = 0.0;
                c[slot + 1] = 0.0;
                c[slot + 4] = 0.0;
                c[slot + 5] = 0.0;
                continue;
            }

            // Tangent rotated by -90 degrees: the direction positive pressure pushes.
            final double nr = tz / len;
            final double nz = -tr / len;

            final double wa, wb;
            if (formulation.isAxisymmetric()) {
                final double circumference = 2.0 * Math.PI * len;
                wa = circumference * (ra / 3.0 + rb / 6.0);
                wb = circumference * (ra / 6.0 + rb / 3.0);
            } else {
                wa = thickness * len * 0.5;
                wb = wa;
            }

            c[slot] = p * nr * wa;
            c[slot + 1] = p * nr * wb;
            c[slot + 4] = p * nz * wa;
            c[slot + 5] = p * nz * wb;
        }
    }

    /** Row-sum lumped mass: m_i = rho * integral(N_i dV). */
    private void lumpMass() {
        final int[] conn = mesh.conn;
        final double[] mr = mesh.r;
        final double[] mz = mesh.z;
        final boolean axi = formulation.isAxisymmetric();

        for (int e = 0; e < mesh.elementCount; e++) {
            final double rho = densityOf[e];
            final int b = e * 4;
            final int n0 = conn[b], n1 = conn[b + 1], n2 = conn[b + 2], n3 = conn[b + 3];
            final double r0 = mr[n0], r1 = mr[n1], r2 = mr[n2], r3 = mr[n3];
            final double z0 = mz[n0], z1 = mz[n1], z2 = mz[n2], z3 = mz[n3];

            for (int g = 0; g < 4; g++) {
                final double xi = GP_XI[g], eta = GP_ETA[g];
                final double xm = 1.0 - xi, xp = 1.0 + xi;
                final double em = 1.0 - eta, ep = 1.0 + eta;
                final double N0 = 0.25 * xm * em, N1 = 0.25 * xp * em;
                final double N2 = 0.25 * xp * ep, N3 = 0.25 * xm * ep;

                final double dN0x = -0.25 * em, dN1x = 0.25 * em, dN2x = 0.25 * ep, dN3x = -0.25 * ep;
                final double dN0e = -0.25 * xm, dN1e = -0.25 * xp, dN2e = 0.25 * xp, dN3e = 0.25 * xm;

                final double J00 = dN0x * r0 + dN1x * r1 + dN2x * r2 + dN3x * r3;
                final double J01 = dN0x * z0 + dN1x * z1 + dN2x * z2 + dN3x * z3;
                final double J10 = dN0e * r0 + dN1e * r1 + dN2e * r2 + dN3e * r3;
                final double J11 = dN0e * z0 + dN1e * z1 + dN2e * z2 + dN3e * z3;
                final double det = J00 * J11 - J01 * J10;

                final double rg = N0 * r0 + N1 * r1 + N2 * r2 + N3 * r3;
                final double dV = axi ? 2.0 * Math.PI * rg * det : thickness * det;

                mass[n0] += rho * N0 * dV;
                mass[n1] += rho * N1 * dV;
                mass[n2] += rho * N2 * dV;
                mass[n3] += rho * N3 * dV;
            }
        }
    }

    // ------------------------------------------------------------------ audits

    /**
     * Total axial momentum, kg m/s.
     *
     * <p>The audit this exists for: every element's axial internal forces sum to zero,
     * because the shape function derivatives are the gradient of a partition of unity and so
     * sum to zero themselves. Nothing internal to the mesh can therefore move this number,
     * and with no pressure load and no axial restraint its entire history is the contact
     * impulse. See {@link RigidWall}.
     */
    public double axialMomentum() {
        double p = 0.0;
        for (int i = 0; i < mesh.nodeCount; i++) p += mass[i] * vz[i];
        return p;
    }

    /**
     * Total radial momentum, kg m/s.
     *
     * <p>Only a conserved quantity under {@link Formulation#PLANE_STRAIN} and
     * {@link Formulation#PLANE_STRESS}, where {@code r} is an ordinary Cartesian direction. An
     * axisymmetric body has no net radial momentum to conserve -- every ring's motion cancels
     * around the circumference, and the hoop term in the internal force is a real source in
     * this sum rather than a bookkeeping artefact. It is the audit {@link Contact} is checked
     * against, because a penalty pair puts equal and opposite forces on the two surfaces and
     * so cannot move it at all.
     */
    public double radialMomentum() {
        double p = 0.0;
        for (int i = 0; i < mesh.nodeCount; i++) p += mass[i] * vr[i];
        return p;
    }

    /** Total kinetic energy, joules. */
    public double kineticEnergy() {
        double ke = 0.0;
        for (int i = 0; i < mesh.nodeCount; i++) {
            ke += mass[i] * (vr[i] * vr[i] + vz[i] * vz[i]);
        }
        return 0.5 * ke;
    }

    /**
     * Recoverable elastic strain energy, joules, computed from the stored stress rather than
     * from displacement: {@code ||s||^2 / 4mu + p^2 / 2K}.
     *
     * <p>Taking it from stress is what keeps it meaningful once the material yields. Half of
     * the work put into a plastic element is not stored, it is dissipated, and an energy
     * computed from total strain would count the dissipated part as recoverable and quietly
     * inflate the denominator of every audit that uses it.
     *
     * <p>Compared against kinetic energy this is the quasi-static audit: when KE / SE has
     * fallen far below one, the dynamic answer has settled onto the static one.
     */
    public double strainEnergy() {
        double se = 0.0;
        for (int e = 0; e < mesh.elementCount; e++) {
            if (gone != null && gone[e]) continue;
            for (int g = 0; g < pointsPerElement; g++) {
                se += pointStrainEnergy(e, g);
            }
        }
        return se;
    }

    /** Recoverable energy stored at one integration point, joules. */
    private double pointStrainEnergy(int e, int g) {
        final int p = e * pointsPerElement + g;
        final int s = J2.COMPONENTS * p;
        final double press = (sig[s] + sig[s + 1] + sig[s + 2]) / 3.0;
        final double a = sig[s] - press, bb = sig[s + 1] - press, c = sig[s + 2] - press;
        final double d = sig[s + 3];
        final double density = (a * a + bb * bb + c * c + 2.0 * d * d) / (4.0 * muOf[e])
                + press * press / (2.0 * bulkOf[e]);
        return density * gaussVolume(e, g);
    }

    /**
     * Energy dissipated by plastic flow, joules.
     *
     * <p>In an impact this is where nearly all the input energy goes -- on the reference
     * Taylor case, around 98 % of it -- which is what makes it the term a dynamic energy
     * balance cannot be written without. Omitting it leaves
     * {@link #kineticEnergy()} and {@link #strainEnergy()} to be compared against each other
     * as small residuals, which will agree to a few parts in a thousand whether or not the
     * run is sound.
     *
     * <p>The work density is Cauchy stress against the rate of deformation, so it is work per
     * unit <em>current</em> volume, while the integration below weights it by the reference
     * volume. The two agree because this return map is exactly isochoric -- the flow direction
     * is deviatoric by construction, so the only volume change anywhere is elastic, a fraction
     * of a percent. That is an argument, not an accident, and it would stop being true for a
     * model with volumetric plasticity.
     *
     * <p>{@link Damage} does not break it, and that is worth saying because most damage models
     * would. Softening scales the size of the yield surface and leaves the flow direction
     * alone, so a damaged point is every bit as isochoric as an undamaged one. The price is
     * that a fully damaged element is a fluid rather than a crack -- see {@link Damage} -- and
     * the price of the alternative would have been this identity.
     */
    public double plasticDissipation() {
        double w = 0.0;
        for (int e = 0; e < mesh.elementCount; e++) {
            for (int g = 0; g < pointsPerElement; g++) {
                final int p = e * pointsPerElement + g;
                w += wPlastic[p] * gaussVolume(e, g);
            }
        }
        return w;
    }

    /** Plastic work per unit volume, J/m^3, one per Gauss point. What thermal softening reads. */
    public double[] plasticWorkDensity() {
        return wPlastic;
    }

    /**
     * Energy consumed by softening, joules, from {@link Damage#dissipation}'s closed form.
     *
     * <p><b>This is not a separate energy.</b> It is already inside
     * {@link #plasticDissipation()} -- softening work is plastic work done at a falling flow
     * stress, and the solver stores it in the same accumulator. What this adds is the ability
     * to say how much of the total went into breaking the part rather than into deforming it,
     * and to say it analytically, so that comparing the two measures the integration error
     * rather than reporting it twice.
     *
     * <p>The number to check it against is {@code G_f} times the crack area. Once a band is
     * fully damaged, this equals that exactly -- which is the mesh independence the
     * regularisation exists for, stated as an identity.
     */
    public double fractureDissipation() {
        if (damage == null) return 0.0;
        double w = 0.0;
        for (int e = 0; e < mesh.elementCount; e++) {
            for (int g = 0; g < pointsPerElement; g++) {
                final int p = e * pointsPerElement + g;
                w += damage.dissipation(p, epsP[p]) * gaussVolume(e, g);
            }
        }
        return w;
    }

    /** Largest damage anywhere in the mesh, 0 to 1. Zero if nothing can fail. */
    public double maxDamage() {
        return damage == null ? 0.0 : damage.maximum(epsP);
    }

    /** Fraction of Gauss points that have started to soften. */
    public double damagedFraction() {
        return damage == null ? 0.0 : damage.begunFraction();
    }

    /** Fraction of Gauss points with no deviatoric strength left. */
    public double failedFraction() {
        return damage == null ? 0.0 : damage.failedFraction(epsP);
    }

    /**
     * Damage at an element, averaged over its Gauss points, 0 to 1. Zero if nothing can fail.
     *
     * <p>The per-element companion to {@link #maxDamage()}, which exists so that damage can
     * be drawn rather than only summarised. Averaging matches {@link #elementPlasticStrain},
     * and under reduced integration there is only one point to average.
     */
    public double elementDamage(int element) {
        if (damage == null) return 0.0;
        double sum = 0.0;
        for (int g = 0; g < pointsPerElement; g++) {
            final int p = element * pointsPerElement + g;
            sum += damage.at(p, epsP[p]);
        }
        return sum / pointsPerElement;
    }

    /**
     * Element whose damage is highest, or -1 if nothing has softened. What the burst harness
     * reads to say <em>where</em> a tube failed, as distinct from at what pressure.
     */
    public int mostDamagedElement() {
        if (damage == null) return -1;
        int worst = -1;
        double max = 0.0;
        for (int e = 0; e < mesh.elementCount; e++) {
            for (int g = 0; g < pointsPerElement; g++) {
                final int p = e * pointsPerElement + g;
                final double d = damage.at(p, epsP[p]);
                if (d > max) {
                    max = d;
                    worst = e;
                }
            }
        }
        return worst;
    }

    /**
     * Energy stored in the hourglass restoring springs, joules. Zero under full integration.
     *
     * <p>This is the audit that makes reduced integration honest. Hourglass stiffness is not
     * physics -- it is a numerical prop holding up a mode the quadrature cannot see -- so the
     * question is always how much of the answer it is carrying. Against strain energy it
     * should be a rounding error. When it climbs into the percent range the mesh is being
     * held together by the prop, and the result is an artefact of the hourglass coefficient
     * rather than of the material. Reporting it is what stops that from going unnoticed,
     * which is exactly how hourglassing normally ships: silently, because the mode carries
     * no strain energy and every other audit stays clean.
     */
    public double hourglassEnergy() {
        if (integration == Integration.FULL || hourglassCoefficient == 0.0) return 0.0;

        // The stored quantity is the restoring force Q = k q, so the energy of the spring
        // that produced it is Q^2 / 2k.
        double energy = 0.0;
        for (int e = 0; e < mesh.elementCount; e++) {
            final double k = hourglassStiffness[e];
            if (k > 0.0) energy += 0.5 * (hgR[e] * hgR[e] + hgZ[e] * hgZ[e]) / k;
        }
        return energy;
    }

    /**
     * Adiabatic temperature at a Gauss point, K.
     *
     * <p>Derived from the plastic work stored there rather than integrated as a field of its
     * own, which is exact rather than a shortcut: over the tens of microseconds an impact
     * lasts, the thermal diffusion length in steel is well under one element, so no point
     * exchanges heat with any other.
     *
     * @throws IllegalStateException if the material has no thermal model, because returning
     *                               room temperature would be indistinguishable from a
     *                               specimen that never heated up
     */
    public double temperature(int point) {
        final Material material = materials.at(point / pointsPerElement);
        if (material.johnsonCook() == null) {
            throw new IllegalStateException(
                    "this material has no thermal model; temperature is not defined");
        }
        return material.johnsonCook().temperature(wPlastic[point], material.density());
    }

    /** Highest temperature anywhere in the mesh, K. */
    public double maxTemperature() {
        double max = Double.NEGATIVE_INFINITY;
        for (int p = 0; p < wPlastic.length; p++) max = Math.max(max, temperature(p));
        return max;
    }

    /** Largest accumulated equivalent plastic strain anywhere in the mesh. */
    public double maxPlasticStrain() {
        double max = 0.0;
        for (double v : epsP) if (v > max) max = v;
        return max;
    }

    /** Fraction of Gauss points that have yielded at any point in their history. */
    public double yieldedFraction() {
        int n = 0;
        for (double v : epsP) if (v > 0.0) n++;
        return (double) n / epsP.length;
    }

    /**
     * Radius of the elastic-plastic boundary: the largest radius at which anything has
     * yielded. Meaningful for the pressurised-cylinder case, where plasticity spreads
     * outward from the bore as a front.
     */
    public double plasticFrontRadius() {
        double front = 0.0;
        for (int e = 0; e < mesh.elementCount; e++) {
            boolean yielded = false;
            for (int g = 0; g < pointsPerElement; g++) {
                if (epsP[e * pointsPerElement + g] > 0.0) yielded = true;
            }
            if (yielded) front = Math.max(front, outerRadius(e));
        }
        return front;
    }

    /**
     * Stress at an element centroid, as {sigma_r, sigma_z, sigma_theta, tau_rz}.
     *
     * <p>Under reduced integration the centroid is the stored point, so this is exact. Under
     * full integration it is the mean of the four Gauss points, which is the centroid value
     * of the bilinear field they define.
     */
    public double[] centroidStress(int element) {
        double[] out = new double[J2.COMPONENTS];
        for (int g = 0; g < pointsPerElement; g++) {
            final int s = J2.COMPONENTS * (element * pointsPerElement + g);
            for (int c = 0; c < J2.COMPONENTS; c++) out[c] += sig[s + c];
        }
        for (int c = 0; c < J2.COMPONENTS; c++) out[c] /= pointsPerElement;
        return out;
    }

    /** Equivalent plastic strain at an element, averaged over its Gauss points. */
    public double elementPlasticStrain(int element) {
        double sum = 0.0;
        for (int g = 0; g < pointsPerElement; g++) sum += epsP[element * pointsPerElement + g];
        return sum / pointsPerElement;
    }

    /** Radius of an element centroid. */
    public double centroidRadius(int element) {
        final int b = element * 4;
        return 0.25 * (mesh.r[mesh.conn[b]] + mesh.r[mesh.conn[b + 1]]
                + mesh.r[mesh.conn[b + 2]] + mesh.r[mesh.conn[b + 3]]);
    }

    /**
     * Axial position of an element centroid, in the <b>current</b> configuration.
     *
     * <p>Deliberately not the reference one, unlike {@link #centroidRadius}. A scene of bodies
     * that move wants to know where a piece is now, and colouring a film by where a piece
     * started produces a picture that never changes.
     */
    public double centroidZ(int element) {
        final int b = element * 4;
        double sum = 0.0;
        for (int k = 0; k < 4; k++) {
            final int n = mesh.conn[b + k];
            sum += mesh.z[n] + uz[n];
        }
        return 0.25 * sum;
    }

    /** Reference axial position of a node, metres. For a caller measuring displacement. */
    public double meshZ(int node) {
        return mesh.z[node];
    }

    private double outerRadius(int element) {
        final int b = element * 4;
        return Math.max(Math.max(mesh.r[mesh.conn[b]], mesh.r[mesh.conn[b + 1]]),
                Math.max(mesh.r[mesh.conn[b + 2]], mesh.r[mesh.conn[b + 3]]));
    }

    /** Volume associated with one Gauss point, for energy integration. */
    private double gaussVolume(int element, int g) {
        final int b = element * 4;
        final int n0 = mesh.conn[b], n1 = mesh.conn[b + 1];
        final int n2 = mesh.conn[b + 2], n3 = mesh.conn[b + 3];
        final double r0 = mesh.r[n0], r1 = mesh.r[n1], r2 = mesh.r[n2], r3 = mesh.r[n3];
        final double z0 = mesh.z[n0], z1 = mesh.z[n1], z2 = mesh.z[n2], z3 = mesh.z[n3];

        final double xi = pointsPerElement == 1 ? 0.0 : GP_XI[g];
        final double eta = pointsPerElement == 1 ? 0.0 : GP_ETA[g];
        final double weight = pointsPerElement == 1 ? W_REDUCED : 1.0;

        final double xm = 1.0 - xi, xp = 1.0 + xi, em = 1.0 - eta, ep = 1.0 + eta;
        final double N0 = 0.25 * xm * em, N1 = 0.25 * xp * em;
        final double N2 = 0.25 * xp * ep, N3 = 0.25 * xm * ep;
        final double dN0x = -0.25 * em, dN1x = 0.25 * em, dN2x = 0.25 * ep, dN3x = -0.25 * ep;
        final double dN0e = -0.25 * xm, dN1e = -0.25 * xp, dN2e = 0.25 * xp, dN3e = 0.25 * xm;

        final double J00 = dN0x * r0 + dN1x * r1 + dN2x * r2 + dN3x * r3;
        final double J01 = dN0x * z0 + dN1x * z1 + dN2x * z2 + dN3x * z3;
        final double J10 = dN0e * r0 + dN1e * r1 + dN2e * r2 + dN3e * r3;
        final double J11 = dN0e * z0 + dN1e * z1 + dN2e * z2 + dN3e * z3;
        final double det = J00 * J11 - J01 * J10;
        final double rg = N0 * r0 + N1 * r1 + N2 * r2 + N3 * r3;

        return weight * (formulation.isAxisymmetric()
                ? 2.0 * Math.PI * rg * det : thickness * det);
    }
}
