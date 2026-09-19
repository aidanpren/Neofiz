package org.neofiz.solver;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.neofiz.core.Formulation;
import org.neofiz.core.JohnsonCook;
import org.neofiz.core.Material;
import org.neofiz.mesh.QuadMesh;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Threading the element loop, held to the only standard worth holding it to: the answer is
 * <b>bit-for-bit</b> the same at any thread count.
 *
 * <p>A tolerance would not do here, and the reason is specific rather than fastidious. The
 * obvious way to thread an explicit kernel -- accumulate nodal forces into a private vector
 * per thread and sum them afterwards, or use an atomic add -- changes the <em>order</em> the
 * contributions to a node are summed in. Floating-point addition is not associative, so the
 * force differs in its last bit, the difference feeds a nonlinear return map, and over the
 * 200,000 steps of a real run it grows. A test with a tolerance would pass at every thread
 * count while the answers silently drifted apart, and the tolerance would have to be widened
 * as the runs got longer -- which is exactly the failure that cannot be distinguished from a
 * physics bug.
 *
 * <p>So the kernel does not scatter at all. Each element writes its own slot, and a second
 * pass has each node sum the elements incident on it in ascending element order -- the order
 * the serial scatter used. Identical arithmetic, not equivalent arithmetic.
 *
 * <p>Two configurations are driven, because they exercise disjoint halves of the step. The
 * pressurised cylinder runs the pressure load, hourglass accumulation and a plastic return map
 * on a mesh that never moves. The Taylor specimen runs corotation, adaptive CFL, contact and
 * a Johnson-Cook solve on a mesh that is crushing -- and its timestep is a reduction over
 * every element, which is the other thing threading has to combine.
 */
class ThreadingTest {

    /** Enough elements to give every thread count a non-trivial chunk to disagree over. */
    private static final int WALL = 12;
    private static final int ALONG = 20;

    // ------------------------------------------------------------------ the two runs

    /** Pressurised plastic cylinder: pressure edges, hourglass, plasticity, static mesh. */
    private static ExplicitSolver cylinder() {
        QuadMesh mesh = QuadMesh.cylinderWall(25.0e-3, 31.0e-3, 6.0e-3, WALL, ALONG);
        ExplicitSolver s = new ExplicitSolver(mesh, Material.STEEL_4340.yielding(400e6, 2.0e9),
                Formulation.AXISYMMETRIC, Integration.REDUCED, 1.0, 0.5);
        s.fixAllAxial();
        s.setPressureRamp(180e6, 2.0e-5);
        s.setRelaxationDamping(0.5, 1.0e6);
        return s;
    }

    /** Taylor specimen: finite strain, contact with friction, Johnson-Cook, adaptive CFL. */
    private static ExplicitSolver taylor() {
        final int nr = 6, nz = 30;
        QuadMesh mesh = QuadMesh.solidCylinder(3.7975e-3, 37.97e-3, nr, nz);
        Material steel = Material.STEEL_4340.withJohnsonCook(JohnsonCook.STEEL_4340);
        ExplicitSolver s = new ExplicitSolver(mesh, steel, Formulation.AXISYMMETRIC,
                Integration.FULL, Kinematics.FINITE_STRAIN, 1.0, 0.5);
        for (int n : QuadMesh.axisNodes(nr, nz)) s.fixRadial(n);
        s.setRigidWall(RigidWall.atZ(0.0, 0.15));
        s.setUniformVelocity(0.0, -181.0);
        return s;
    }

    /**
     * Everything the step writes, flattened. Displacement and velocity are downstream of the
     * force, and the force is what threading touches, so a difference anywhere upstream lands
     * in one of these.
     */
    private static double[][] state(ExplicitSolver s) {
        return new double[][]{
                s.stress(), s.plasticStrain(), s.plasticWorkDensity(),
                s.radialDisplacement(), s.axialDisplacement(),
                s.forceR(), s.forceZ(), {s.timestep()}
        };
    }

    private static double[][] snapshot(ExplicitSolver s, int steps, int threads) {
        s.setThreads(threads);
        try {
            s.run(steps);
            double[][] live = state(s);
            double[][] copy = new double[live.length][];
            for (int i = 0; i < live.length; i++) copy[i] = live[i].clone();
            return copy;
        } finally {
            s.shutdown();
        }
    }

    private static void assertIdenticalBits(double[][] one, double[][] many, int threads) {
        for (int a = 0; a < one.length; a++) {
            assertEquals(one[a].length, many[a].length, "array " + a + " changed length");
            for (int i = 0; i < one[a].length; i++) {
                long serial = Double.doubleToRawLongBits(one[a][i]);
                long parallel = Double.doubleToRawLongBits(many[a][i]);
                if (serial != parallel) {
                    throw new AssertionError(String.format(
                            "array %d entry %d differs at %d threads: %.17e vs %.17e "
                                    + "(bits %016x vs %016x). The gather is not reproducing the "
                                    + "serial summation order.",
                            a, i, threads, one[a][i], many[a][i], serial, parallel));
                }
            }
        }
    }

    // ------------------------------------------------------------------ the claim

    @ParameterizedTest(name = "{0} threads")
    @ValueSource(ints = {2, 3, 5, 8})
    @DisplayName("a pressurised plastic cylinder is bit-for-bit identical at any thread count")
    void cylinderIsBitIdentical(int threads) {
        final int steps = 400;
        assertIdenticalBits(snapshot(cylinder(), steps, 1),
                snapshot(cylinder(), steps, threads), threads);
    }

    @ParameterizedTest(name = "{0} threads")
    @ValueSource(ints = {2, 3, 5, 8})
    @DisplayName("so is a Taylor impact, where the mesh moves and the timestep is a reduction")
    void taylorIsBitIdentical(int threads) {
        final int steps = 600;
        assertIdenticalBits(snapshot(taylor(), steps, 1),
                snapshot(taylor(), steps, threads), threads);
    }

    /**
     * The identity has to survive the pool being reconfigured mid-run, because that is how a
     * solver gets reused -- and because a stale chunk boundary or a stale scratch array would
     * otherwise show up only in whichever run happened to change its thread count.
     */
    @Test
    @DisplayName("changing the thread count mid-run changes nothing")
    void threadCountMayChangeMidRun() {
        ExplicitSolver serial = cylinder();
        ExplicitSolver shifting = cylinder();
        try {
            serial.run(400);
            for (int t : new int[]{1, 4, 2, 7, 1, 3}) {
                shifting.setThreads(t);
                shifting.run(400 / 6);
            }
            shifting.run(400 - 6 * (400 / 6));
            assertIdenticalBits(state(serial), state(shifting), -1);
        } finally {
            serial.shutdown();
            shifting.shutdown();
        }
    }

    /**
     * More threads than elements. The chunking must hand some workers an empty range rather
     * than an inverted one, and the timestep reduction must still see a real minimum -- an
     * empty chunk reports positive infinity, and one of those leaking into the minimum would
     * be invisible until the run went unstable.
     */
    @Test
    @DisplayName("more threads than elements is an empty chunk, not a wrong answer")
    void moreThreadsThanElements() {
        QuadMesh mesh = QuadMesh.cylinderWall(25.0e-3, 27.0e-3, 1.0e-3, 2, 1);
        ExplicitSolver s = new ExplicitSolver(mesh, Material.STEEL_4340,
                Formulation.AXISYMMETRIC, Integration.REDUCED, Kinematics.FINITE_STRAIN, 1.0, 0.5);
        s.fixAllAxial();
        s.setPressureRamp(50e6, 1.0e-6);
        try {
            s.setThreads(16);
            double dt0 = s.timestep();
            s.run(50);
            assertTrue(s.timestep() > 0.0 && s.timestep() < 2.0 * dt0,
                    "the timestep reduction picked up an empty chunk's infinity");
            assertTrue(Double.isFinite(s.stress()[0]), "the kernel did not run on every element");
        } finally {
            s.shutdown();
        }
    }

    @Test
    @DisplayName("the pool refuses a thread count below one and is safe to stop twice")
    void poolLifecycle() {
        ExplicitSolver s = cylinder();
        assertEquals(1, s.threads());
        assertThrows(IllegalArgumentException.class, () -> s.setThreads(0));
        s.setThreads(4);
        assertEquals(4, s.threads());
        s.run(10);
        s.shutdown();
        s.shutdown();
        assertEquals(1, s.threads(), "a stopped pool must fall back to the calling thread");
        s.run(10);   // and must still run
    }
}
