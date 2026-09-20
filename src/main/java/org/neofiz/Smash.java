package org.neofiz;

import org.neofiz.core.Formulation;
import org.neofiz.core.JohnsonCook;
import org.neofiz.core.Material;
import org.neofiz.core.Materials;
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
 * A steel slug through a stacked wall of copper blocks.
 *
 * <p>Thirteen separate bodies in one solve, resting on each other and on the floor, and one of
 * them arriving at several hundred metres a second. It is the first scene in this project that
 * is a <em>scene</em> rather than a specimen: nothing here is measuring anything, and what it
 * demonstrates is that the pieces now compose.
 *
 * <ul>
 *   <li>{@link Outline#assemble} -- each block gets its own nodes, so the stack is a stack and
 *       not a monolith. They start exactly touching, because an explicit solve runs in
 *       microseconds and gravity settles things in milliseconds, so anything that has to fall
 *       into place first is unaffordable.</li>
 *   <li>{@link Contact} -- what holds the stack up and what transmits the slug. Twelve bodies
 *       pressing on one another with a Coulomb cone between them.</li>
 *   <li>{@link ExplicitSolver#setGravity} -- on, and honest about being nearly irrelevant: over
 *       the few hundred microseconds this runs, a g moves anything by a fraction of a micron.
 *       The wall comes apart because it is hit, not because it falls.</li>
 *   <li>{@link Materials} -- steel slug, copper wall, painted per element.</li>
 * </ul>
 *
 * <p>It also breaks, which it did not when this scene was first written. With
 * {@link ExplicitSolver#setErosion} on, an element softened to nothing is deleted and the free
 * surface closes around the hole, so the slug punches through rather than pushing a permanent
 * dent ahead of it. At 700 m/s that is 32 of 652 elements, none of them inverted, and 15 J of
 * stored energy that leaves with them.
 *
 * <p>What it does <b>not</b> do is shatter. Copper is ductile: it perforates, petals and
 * tears, and that is the right answer for copper. A wall that comes apart into fragments needs
 * a brittle material, and a brittle material needs a mesh this scene cannot afford -- the
 * crack-band limit is {@code h <= 2*E*G_f/sigma_f^2}, which for an alumina puts the largest
 * regularisable element at a few microns against the millimetre used here. That is a real
 * constraint rather than a missing feature, and it is why the wall is made of something
 * ductile.
 *
 * <pre>  ./gradlew smash --args="&lt;speed in m/s&gt;"</pre>
 *
 * <p>Writes {@code runs/smash.js}, which the scene viewer reads.
 */
public final class Smash {

    private static final double MM = 1.0e-3;
    private static final double CELL = 1.0 * MM;
    private static final double DEPTH = 20.0 * MM;

    /** Blocks across, blocks up, and the side of each. */
    private static final int COLUMNS = 3;
    private static final int ROWS = 4;
    private static final double BLOCK = 7.0 * MM;

    private static final double DEFAULT_SPEED = 700.0;
    private static final double DURATION = 150.0e-6;
    private static final int FRAMES = 130;

    private static final Material STEEL =
            Material.STEEL_4340.withJohnsonCook(JohnsonCook.STEEL_4340);
    private static final Material COPPER =
            Material.COPPER_OFHC.withJohnsonCook(JohnsonCook.COPPER_OFHC);

    private Smash() {
    }

    public static void main(String[] args) {
        final double speed = args.length > 0 ? Double.parseDouble(args[0]) : DEFAULT_SPEED;

        // The wall, laid out corner to corner so the blocks touch exactly rather than
        // overlapping or leaving a gap the solve cannot afford to close.
        final Outline[] bodies = new Outline[COLUMNS * ROWS + 1];
        int at = 0;
        for (int j = 0; j < ROWS; j++) {
            for (int i = 0; i < COLUMNS; i++) {
                bodies[at++] = Outline.rectangle(
                        i * BLOCK, j * BLOCK, (i + 1) * BLOCK, (j + 1) * BLOCK);
            }
        }
        final double slugZ = 1.5 * BLOCK;
        final Outline slug = Outline.rectangle(
                -16.0 * MM, slugZ - 4.0 * MM, -8.0 * MM, slugZ + 4.0 * MM);
        bodies[at] = slug;

        final QuadMesh mesh = Outline.assemble(CELL, bodies);
        final boolean[] isSlug = new boolean[mesh.elementCount];
        for (int e = 0; e < mesh.elementCount; e++) {
            isSlug[e] = slug.contains(mesh.centroidRadius(e), mesh.centroidZ(e));
        }
        final Materials materials = Materials.byElement(mesh.elementCount,
                e -> isSlug[e] ? STEEL : COPPER);

        final ExplicitSolver solver = new ExplicitSolver(mesh, materials,
                Formulation.PLANE_STRAIN, Integration.REDUCED, Kinematics.FINITE_STRAIN,
                DEPTH, 0.5);
        solver.setGravity(0.0, -9.81);
        solver.setRigidWall(RigidWall.atZ(0.0, 0.4));
        solver.setContact(new Contact(mesh).withFriction(0.3));
        for (int i = 0; i < mesh.nodeCount; i++) {
            if (mesh.r[i] < 0.0) solver.setVelocity(i, speed, 0.0);
        }

        final double[] failureStrain = new double[mesh.elementCount];
        final double[] fractureEnergy = new double[mesh.elementCount];
        for (int e = 0; e < mesh.elementCount; e++) {
            failureStrain[e] = isSlug[e] ? 0.60 : 0.90;
            fractureEnergy[e] = isSlug[e] ? 60.0e3 : 100.0e3;
        }
        solver.setDamage(fractureEnergy, failureStrain);
        solver.setErosion(1.0);

        System.out.println("NEOFIZ  a slug through a wall");
        System.out.println("=".repeat(76));
        System.out.printf(Locale.ROOT,
                "Scene     %d bodies, %d elements at %.1f mm, %.0f mm deep%n",
                bodies.length, mesh.elementCount, 1e3 * CELL, 1e3 * DEPTH);
        System.out.printf(Locale.ROOT,
                "Wall      %d by %d blocks of %s, %.0f mm on a side, stacked touching%n",
                COLUMNS, ROWS, COPPER.name(), 1e3 * BLOCK);
        System.out.printf(Locale.ROOT, "Slug      %s, 8 by 8 mm, %.0f m/s%n", STEEL.name(), speed);
        System.out.printf(Locale.ROOT, "Contact   penalty, friction 0.3; floor friction 0.4%n");
        System.out.printf(Locale.ROOT, "Step      %.3f ns at the start%n%n",
                1e9 * solver.timestep());

        final Film film = new Film(mesh).planar()
                .colourBy("body", "", 0.0, 1.0, (s, e) -> isSlug[e] ? 1.0 : 0.0)
                .colourBy("speed", "m/s", Film.AUTO, Film.AUTO, Smash::elementSpeed)
                .colourBy("plastic strain", "", Film.AUTO, Film.AUTO,
                        (s, e) -> s.elementPlasticStrain(e))
                .colourBy("damage", "", 0.0, 1.0, (s, e) -> s.elementDamage(e))
                .aliveIf((s, e) -> !s.isEroded(e))
                .readout("time, us", s -> 1e6 * s.time())
                .readout("kinetic energy, J", ExplicitSolver::centredKineticEnergy)
                .readout("plastic dissipation, J", ExplicitSolver::plasticDissipation)
                .readout("contact pairs", s -> s.contact().pairCount())
                .readout("elements gone", ExplicitSolver::erodedElements);

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

        final Contact contact = solver.contact();
        System.out.printf(Locale.ROOT,
                "Energy    %.1f J still kinetic, %.1f J dissipated plastically%n",
                solver.centredKineticEnergy(), solver.plasticDissipation());
        System.out.printf(Locale.ROOT,
                "Contact   %d pairs at the end, deepest penetration %.1f %% of a cell, "
                        + "%d escapes%n",
                contact.pairCount(), 100.0 * contact.maxPenetration() / CELL,
                contact.escapedNodes());
        System.out.printf(Locale.ROOT,
                "Friction  %.2f J lost between bodies, %.2f J on the floor%n",
                contact.dissipation(), 0.0);
        System.out.printf(Locale.ROOT,
                "Damage    worst %.3f, %.1f %% of points yielded, hourglass %.2f J%n",
                solver.maxDamage(), 100.0 * solver.yieldedFraction(), solver.hourglassEnergy());
        final double internal = solver.strainEnergy() + solver.plasticDissipation()
                + solver.fractureDissipation() + solver.erodedEnergy();
        System.out.printf(Locale.ROOT,
                "          hourglass is %.1f %% of the %.0f J of internal energy, %d clamped points%n",
                100.0 * solver.hourglassEnergy() / internal, internal,
                solver.damage().clampedPoints());
        System.out.printf(Locale.ROOT,
                "Eroded    %d of %d elements (%d of them inverted), %.2f J went with them%n",
                solver.erodedElements(), mesh.elementCount, solver.invertedElements(),
                solver.erodedEnergy());
        System.out.printf(Locale.ROOT, "Run       %d steps, %.0f us simulated, %.1f s of wall "
                + "clock%n%n", solver.steps(), 1e6 * solver.time(), wall);

        RunLog.of("smash")
                .parameter("speed", speed)
                .parameter("cell", CELL)
                .parameter("blocks", COLUMNS * ROWS)
                .result("kinetic", solver.centredKineticEnergy())
                .result("plastic", solver.plasticDissipation())
                .result("friction", contact.dissipation())
                .result("eroded", solver.erodedElements())
                .result("hourglass", solver.hourglassEnergy())
                .result("steps", solver.steps())
                .result("wallClock", wall)
                .commit();

        final Path file = film.write("smash.js",
                "Slug through a stacked wall, " + (int) speed + " m/s");
        System.out.printf(Locale.ROOT, "%d frames written to %s (%.1f MB)%n",
                film.size(), file, file.toFile().length() / (1024.0 * 1024.0));
        solver.shutdown();
    }

    /** Speed at an element's corners, m/s. What makes a scene of moving bodies legible. */
    private static double elementSpeed(ExplicitSolver s, int element) {
        return s.elementSpeed(element);
    }
}
