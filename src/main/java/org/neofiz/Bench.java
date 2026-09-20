package org.neofiz;

import org.neofiz.core.Formulation;
import org.neofiz.core.JohnsonCook;
import org.neofiz.core.Material;
import org.neofiz.mesh.Outline;
import org.neofiz.mesh.QuadMesh;
import org.neofiz.report.Csv;
import org.neofiz.solver.Contact;
import org.neofiz.solver.ExplicitSolver;
import org.neofiz.solver.Integration;
import org.neofiz.solver.Kinematics;
import org.neofiz.solver.RigidWall;

import java.util.Arrays;
import java.util.Locale;

/**
 * Where the time actually goes.
 *
 * <p>The cost model in the README was an estimate, and estimates about explicit solvers are
 * usually wrong in the same direction: the element kernel is the part that looks expensive and
 * is not, because it is a fixed amount of arithmetic per element with no branches and no
 * pointer chasing. Everything else in a scene -- contact search, the wall, the nodal update --
 * is where the surprises live.
 *
 * <p>So this measures instead. The unit is the <b>element-step</b>, which is the only quantity
 * that lets a 600-element scene and a 60 000-element scene be compared: doubling the mesh
 * doubles the work per step <em>and</em> halves the step, so cost goes as the square of
 * resolution and quoting a frame rate for one mesh says nothing about another.
 *
 * <pre>  ./gradlew bench</pre>
 *
 * <p>Writes {@code runs/bench.csv}.
 */
public final class Bench {

    private static final double MM = 1.0e-3;
    private static final double DEPTH = 10.0 * MM;
    private static final Material COPPER =
            Material.COPPER_OFHC.withJohnsonCook(JohnsonCook.COPPER_OFHC);

    /** What a configuration cost, in element-steps per second. */
    private record Result(String name, int elements, int steps, double seconds) {

        double rate() {
            return elements * (double) steps / seconds;
        }
    }

    private Bench() {
    }

    public static void main(String[] args) {
        System.out.println("NEOFIZ  where the time goes");
        System.out.println("=".repeat(76));
        System.out.printf(Locale.ROOT, "%d cores reported%n%n",
                Runtime.getRuntime().availableProcessors());

        final Result[] results = {
                run("bare elements", 40, false, false, false, 1),
                run("+ damage", 40, true, false, false, 1),
                run("+ rigid floor", 40, true, true, false, 1),
                run("+ contact", 40, true, true, true, 1),
                run("bare, 4 threads", 40, false, false, false, 4),
                run("+ contact, 4 threads", 40, true, true, true, 4),
                run("bare elements, 4x mesh", 80, false, false, false, 1),
                run("+ contact, 4x mesh", 80, true, true, true, 1),
        };

        System.out.printf(Locale.ROOT, "%-24s %9s %8s %9s %14s%n",
                "configuration", "elements", "steps", "seconds", "element-steps/s");
        System.out.println("-".repeat(76));
        for (Result r : results) {
            System.out.printf(Locale.ROOT, "%-24s %9d %8d %9.2f %14.3e%n",
                    r.name(), r.elements(), r.steps(), r.seconds(), r.rate());
        }

        final double bare = results[0].rate();
        System.out.println();
        System.out.printf(Locale.ROOT,
                "Contact costs %.1fx the bare element kernel on this scene.%n",
                bare / results[3].rate());
        System.out.printf(Locale.ROOT,
                "Four threads buy %.2fx without contact and %.2fx with it.%n",
                results[4].rate() / bare, results[5].rate() / results[3].rate());
        System.out.printf(Locale.ROOT,
                "A 4x mesh runs at %.2fx the per-element rate, so the kernel scales.%n%n",
                results[6].rate() / bare);

        // What a second of simulated time would cost, at this rate, for a scene of a given
        // size. This is the number that decides whether a sandbox is interactive, and it is
        // the one an estimate gets wrong: the step shrinks with the mesh, so the cost of a
        // second of physics grows as the cube of the linear resolution.
        System.out.printf(Locale.ROOT, "%-14s %12s %14s %16s%n",
                "cell size", "elements", "step", "hours per sim. second");
        System.out.println("-".repeat(76));
        for (double cell : new double[] {8.0, 4.0, 2.0, 1.0, 0.5}) {
            final double side = 0.2;                       // a 200 mm square of material
            final double elements = (side / (cell * MM)) * (side / (cell * MM));
            final double step = 0.5 * cell * MM / COPPER.dilatationalWaveSpeed();
            final double stepsPerSecond = 1.0 / step;
            final double hours = elements * stepsPerSecond / results[3].rate() / 3600.0;
            System.out.printf(Locale.ROOT, "%10.1f mm %12.0f %12.1f ns %16.1f%n",
                    cell, elements, 1e9 * step, hours);
        }

        final Csv csv = Csv.of("configuration", "elements", "steps", "seconds",
                "elementStepsPerSecond");
        for (Result r : results) {
            csv.row(r.name(), r.elements(), r.steps(), r.seconds(), r.rate());
        }
        System.out.printf(Locale.ROOT, "%nwritten to %s%n", csv.write("bench.csv"));
    }

    /**
     * Times one configuration.
     *
     * <p>Every run is preceded by a shorter one on a fresh solver, because the first few
     * thousand steps of a JVM run are spent compiling the kernel and are not what anyone wants
     * to know. The returned time is the second run only.
     */
    private static Result run(String name, int side, boolean damage, boolean floor,
                              boolean contact, int threads) {
        time(build(side, damage, floor, contact, threads), 400);
        final ExplicitSolver s = build(side, damage, floor, contact, threads);
        final int steps = 2000;
        final double seconds = time(s, steps);
        final int elements = side * side;
        s.shutdown();
        return new Result(name, elements, steps, seconds);
    }

    private static double time(ExplicitSolver s, int steps) {
        final long t0 = System.nanoTime();
        s.run(steps);
        final double seconds = (System.nanoTime() - t0) / 1e9;
        s.shutdown();
        return seconds;
    }

    private static ExplicitSolver build(int side, boolean damage, boolean floor,
                                        boolean contact, int threads) {
        // Two bodies rather than one, so that the contact case has something to find and the
        // others carry exactly the same mesh. Comparing configurations on different meshes
        // would measure the meshes.
        final double half = side * MM / 2.0;
        final QuadMesh mesh = Outline.assemble(MM,
                Outline.rectangle(0.0, 0.0, side * MM, half),
                Outline.rectangle(0.0, half, side * MM, side * MM));

        final ExplicitSolver s = new ExplicitSolver(mesh, COPPER, Formulation.PLANE_STRAIN,
                Integration.REDUCED, Kinematics.FINITE_STRAIN, DEPTH, 0.5);
        s.setThreads(threads);
        if (damage) {
            final double[] eps = new double[mesh.elementCount];
            final double[] gf = new double[mesh.elementCount];
            Arrays.fill(eps, 0.8);
            Arrays.fill(gf, 100.0e3);
            s.setDamage(gf, eps);
            s.setErosion(1.0);
        }
        if (floor) s.setRigidWall(RigidWall.atZ(0.0, 0.3));
        if (contact) s.setContact(new Contact(mesh).withFriction(0.3));
        for (int i = 0; i < mesh.nodeCount; i++) {
            s.setVelocity(i, 0.0, mesh.z[i] > half ? -120.0 : 0.0);
        }
        return s;
    }
}
