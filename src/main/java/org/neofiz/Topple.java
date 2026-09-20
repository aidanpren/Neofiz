package org.neofiz;

import org.neofiz.core.Formulation;
import org.neofiz.core.Material;
import org.neofiz.mesh.Outline;
import org.neofiz.mesh.QuadMesh;
import org.neofiz.render.Film;
import org.neofiz.report.RunLog;
import org.neofiz.solver.Contact;
import org.neofiz.solver.ExplicitSolver;
import org.neofiz.solver.Integration;
import org.neofiz.solver.Kinematics;
import org.neofiz.solver.RigidWall;

import java.nio.file.Path;
import java.util.Locale;

/**
 * A leaning stack of blocks, left alone under gravity until it falls over.
 *
 * <p>The first scene here where <b>nothing is thrown</b>. Every other film in this project
 * starts with something travelling at a few hundred metres a second, because that is what an
 * explicit solver can afford: it runs in microseconds and gravity works in milliseconds, so
 * anything that has to fall before the interesting part begins costs a hundred thousand steps
 * of watching nothing happen. {@link ExplicitSolver#setMassScaling} is what makes this
 * affordable, and this scene is what it was for.
 *
 * <h2>Why toppling is the right case for scaled mass</h2>
 *
 * Mass scaling is a deliberate lie -- a body with {@code f} times its mass carries {@code f}
 * times the momentum -- so it matters that the motion here does not notice. A block rotating
 * about its corner under its own weight turns at
 *
 * <pre>  alpha = torque / inertia = m*g*d / (m*k^2) = g*d / k^2</pre>
 *
 * and the mass cancels <em>exactly</em>. Every purely gravity-driven motion has that property:
 * the force and the inertia are the same mass, so scaling it changes nothing about the
 * trajectory. What scaling does change is the acoustic impedance, {@code rho*c}, which grows
 * as {@code sqrt(f)} -- so a block <em>landing</em> hits harder than it should, by that
 * factor. The scaling here is chosen to keep that below the yield stress rather than for the
 * largest step it could get away with.
 *
 * <h2>Why it falls</h2>
 *
 * Each block is set back by two fifths of its width. Every block is individually over the one
 * below it, which is what makes the stack look like it should stand; the centre of mass of
 * everything above the first joint is not, which is what makes it fall. The overhang problem,
 * built as a scene instead of as an inequality.
 *
 * <pre>  ./gradlew topple</pre>
 *
 * <p>Writes {@code runs/topple.js} for the scene viewer.
 */
public final class Topple {

    private static final double MM = 1.0e-3;
    private static final double CELL = 4.0 * MM;
    private static final double DEPTH = 40.0 * MM;

    private static final int BLOCKS = 5;
    private static final double SIDE = 40.0 * MM;
    /** How far each block is set back from the one below, as a fraction of its width. */
    private static final double OVERHANG = 0.4;

    /**
     * How much longer a step to buy.
     *
     * <p>Twenty, which is four hundred times the mass. The limit is not stability -- the CFL
     * rule takes care of that -- but the acoustic impedance, which grows as the square root of
     * the factor and so makes every landing twenty times harder than it should be. At the
     * speeds a 200 mm stack reaches falling over, that is still well below the yield stress of
     * the steel it is made of; at a hundred it would not be.
     */
    private static final double STEP_FACTOR = 20.0;

    private static final double DURATION = 320.0e-3;
    private static final int FRAMES = 140;

    /** Elastic, because a toppling stack should read as rigid bodies, not as putty. */
    private static final Material STEEL = Material.STEEL_4340;

    private Topple() {
    }

    public static void main(String[] args) {
        final Outline[] blocks = new Outline[BLOCKS];
        for (int i = 0; i < BLOCKS; i++) {
            final double left = i * OVERHANG * SIDE;
            blocks[i] = Outline.rectangle(left, i * SIDE, left + SIDE, (i + 1) * SIDE);
        }

        final QuadMesh mesh = Outline.assemble(CELL, blocks);
        final ExplicitSolver solver = new ExplicitSolver(mesh, STEEL,
                Formulation.PLANE_STRAIN, Integration.REDUCED, Kinematics.FINITE_STRAIN,
                DEPTH, 0.5);

        final double plain = solver.timestep();
        solver.setMassScaling(plain * STEP_FACTOR);
        solver.setGravity(0.0, -9.81);
        solver.setRigidWall(RigidWall.atZ(0.0, 0.5));
        solver.setContact(new Contact(mesh).withFriction(0.5).withDamping(0.05));

        // Where the centre of mass of the stack above the first joint sits, against the block
        // it is standing on. This is the whole reason the thing falls, and it is cheaper to
        // print it than to explain it.
        double above = 0.0;
        for (int i = 1; i < BLOCKS; i++) above += i * OVERHANG * SIDE + 0.5 * SIDE;
        above /= (BLOCKS - 1);

        System.out.println("NEOFIZ  a stack that falls over");
        System.out.println("=".repeat(76));
        System.out.printf(Locale.ROOT,
                "Stack     %d blocks of %s, %.0f mm on a side, each set back %.0f %%%n",
                BLOCKS, STEEL.name(), 1e3 * SIDE, 100.0 * OVERHANG);
        System.out.printf(Locale.ROOT,
                "          %d elements at %.0f mm, %.0f mm tall%n",
                mesh.elementCount, 1e3 * CELL, 1e3 * BLOCKS * SIDE);
        System.out.printf(Locale.ROOT,
                "Balance   everything above the first joint has its centre of mass at "
                        + "%.0f mm,%n          over a block that ends at %.0f mm -- so it "
                        + "cannot stand%n", 1e3 * above, 1e3 * SIDE);
        System.out.printf(Locale.ROOT,
                "Mass      scaled %.0fx for a %.1f us step, up from %.3f us%n",
                solver.massScale(), 1e6 * solver.timestep(), 1e6 * plain);
        System.out.printf(Locale.ROOT,
                "          landings therefore hit %.1fx too hard; gravity-driven motion is "
                        + "exact%n%n", Math.sqrt(solver.massScale()));

        final Film film = new Film(mesh).planar()
                .colourBy("height", "mm", Film.AUTO, Film.AUTO,
                        (s, e) -> 1e3 * s.centroidZ(e))
                .colourBy("speed", "m/s", Film.AUTO, Film.AUTO,
                        (s, e) -> s.elementSpeed(e))
                .colourBy("von Mises", "MPa", Film.AUTO, Film.AUTO,
                        (s, e) -> mises(s.centroidStress(e)) / 1e6)
                .readout("time, ms", s -> 1e3 * s.time())
                .readout("kinetic energy, J", ExplicitSolver::centredKineticEnergy)
                .readout("contact pairs", s -> s.contact().pairCount())
                .readout("highest point, mm", Topple::highest)
                .readout("friction lost, J", s -> s.contact().dissipation());

        final long t0 = System.nanoTime();
        final double interval = DURATION / FRAMES;
        film.capture(solver);
        for (int f = 0; f < FRAMES; f++) {
            final double until = (f + 1) * interval;
            while (solver.time() < until) {
                solver.run(Math.max(1, (int) ((until - solver.time()) / solver.timestep())));
            }
            film.capture(solver);
        }
        final double wall = (System.nanoTime() - t0) / 1e9;

        System.out.printf(Locale.ROOT, "Fell      from %.0f mm tall to %.0f mm%n",
                1e3 * BLOCKS * SIDE, highest(solver));
        System.out.printf(Locale.ROOT,
                "Energy    %.2f J still kinetic, %.2f J lost to friction, %.3f J stored%n",
                solver.centredKineticEnergy(), solver.contact().dissipation(),
                solver.strainEnergy());
        System.out.printf(Locale.ROOT,
                "Contact   %d pairs at the end, deepest penetration %.2f %% of a cell%n",
                solver.contact().pairCount(),
                100.0 * solver.contact().maxPenetration() / CELL);
        System.out.printf(Locale.ROOT, "Run       %d steps, %.0f ms simulated, %.1f s of wall "
                + "clock%n%n", solver.steps(), 1e3 * solver.time(), wall);

        RunLog.of("topple")
                .parameter("blocks", BLOCKS)
                .parameter("overhang", OVERHANG)
                .parameter("cell", CELL)
                .parameter("stepFactor", STEP_FACTOR)
                .result("finalHeight", highest(solver))
                .result("kinetic", solver.centredKineticEnergy())
                .result("friction", solver.contact().dissipation())
                .result("steps", solver.steps())
                .result("wallClock", wall)
                .commit();

        final Path file = film.write("topple.js", "A leaning stack, left alone");
        System.out.printf(Locale.ROOT, "%d frames written to %s (%.1f MB)%n",
                film.size(), file, file.toFile().length() / (1024.0 * 1024.0));
        solver.shutdown();
    }

    /** The highest node in the scene, millimetres. How far the stack has come down. */
    private static double highest(ExplicitSolver s) {
        final double[] uz = s.axialDisplacement();
        double top = 0.0;
        for (int i = 0; i < uz.length; i++) top = Math.max(top, s.meshZ(i) + uz[i]);
        return 1e3 * top;
    }

    private static double mises(double[] s) {
        final double a = s[0] - s[1], b = s[1] - s[2], c = s[2] - s[0];
        return Math.sqrt(0.5 * (a * a + b * b + c * c) + 3.0 * s[3] * s[3]);
    }
}
