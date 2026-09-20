package org.neofiz;

import org.neofiz.core.Formulation;
import org.neofiz.core.JohnsonCook;
import org.neofiz.core.Material;
import org.neofiz.core.Materials;
import org.neofiz.mesh.Outline;
import org.neofiz.mesh.QuadMesh;
import org.neofiz.render.Film;
import org.neofiz.report.RunLog;
import org.neofiz.solver.ExplicitSolver;
import org.neofiz.solver.Integration;
import org.neofiz.solver.Kinematics;
import org.neofiz.solver.RigidWall;

import java.nio.file.Path;
import java.util.Locale;

/**
 * An L-shaped bracket, steel-footed and copper-legged, thrown at an anvil on its corner.
 *
 * <p>Nothing in this scene is a new physical model. The return map, the corotational update,
 * the hourglass control, the crack-band softening and the kinematic wall are the same code the
 * burst gates run on. What is new is that the scene can be <em>described</em> at all: it is
 * not a solid of revolution, its outline is not one of two hardcoded generators, and it is not
 * made of one substance.
 *
 * <p>Those three refusals were the whole distance between a validation harness and a sandbox,
 * and each is removed by a different piece:
 *
 * <ul>
 *   <li>{@link Formulation#PLANE_STRAIN} -- a cross-section of a long part rather than a body
 *       spun about an axis. The solver has carried the flag since M0 and nothing had ever used
 *       it. An axisymmetric tube cannot fall over sideways; this can.</li>
 *   <li>{@link Outline} -- the L is six corners, turned by {@link Outline#rotated} and set down
 *       on the anvil. The old mesher could express a rectangle and a rectangle touching the
 *       axis.</li>
 *   <li>{@link Materials} -- steel below, copper above, painted by asking the foot's own
 *       outline which elements are inside it. The interface is welded, which is what sharing
 *       nodes means, and is the honest limit here: two bodies that should be able to
 *       <em>separate</em> still need contact.</li>
 * </ul>
 *
 * <h2>The tilt is the experiment, and it was not obvious</h2>
 *
 * The first version of this scene landed the bracket flat, on the theory that a stiff foot
 * under a soft off-centre leg would tip on its own. It does not. Measured off the film, the top
 * of the leg turned <b>1.08 degrees</b> over the whole event and the section's centroid moved
 * sideways by 16 microns: it squats, it spreads by 17 %, and it stays upright. The reason is
 * that the entire foot arrives at the anvil in the same instant, so the reaction is distributed
 * under the whole base and there is no moment arm anywhere. Mass being off-centre is not a
 * torque when every part of the contact is supported.
 *
 * <p>Landing it on a corner supplies what was missing. One cell of the foot takes the whole
 * reaction while the rest of the body is still travelling, so there is a real arm, and the
 * bracket rotates about the contact as it deforms. That is a rigid-body rotation superposed on
 * a large plastic deformation -- the case the corotational update exists for, and a case the
 * axisymmetric harness could not pose. Pass a tilt of zero to see the flat landing again.
 *
 * <pre>  ./gradlew impact --args="&lt;speed in m/s&gt; &lt;tilt in degrees&gt;"</pre>
 *
 * <p>Writes {@code runs/bracket.js} for the scene viewer.
 */
public final class Impact {

    /** Element edge, metres. Sets the resolution and, through the CFL rule, the timestep. */
    private static final double CELL = 1.0e-3;

    /** Out-of-plane depth of the section, metres. Plane strain treats it as nominal. */
    private static final double DEPTH = 10.0e-3;

    /** The bracket, in metres: a 30 mm foot 12 mm deep, with a 12 mm leg 40 mm tall on it. */
    private static final double SPAN = 30.0e-3;
    private static final double FOOT = 12.0e-3;
    private static final double LEG = 12.0e-3;
    private static final double TALL = 40.0e-3;

    private static final double DEFAULT_SPEED = 150.0;
    private static final double DEFAULT_TILT = 18.0;
    private static final double DURATION = 400.0e-6;
    private static final int FRAMES = 140;

    /**
     * Coulomb friction at the anvil.
     *
     * <p>Not a detail once the landing is on a corner. A frictionless anvil lets the contact
     * point slide, so the bracket skates instead of pivoting and the rotation never develops.
     * 0.3 is dry steel on steel.
     */
    private static final double FRICTION = 0.3;

    /** Hard, strong, barely deforms. */
    private static final Material STEEL =
            Material.STEEL_4340.withJohnsonCook(JohnsonCook.STEEL_4340);
    /** Soft, strongly hardening, mushrooms. */
    private static final Material COPPER =
            Material.COPPER_OFHC.withJohnsonCook(JohnsonCook.COPPER_OFHC);

    private Impact() {
    }

    public static void main(String[] args) {
        final double speed = args.length > 0 ? Double.parseDouble(args[0]) : DEFAULT_SPEED;
        final double tilt = Math.toRadians(
                args.length > 1 ? Double.parseDouble(args[1]) : DEFAULT_TILT);

        // Six corners, anticlockwise from the bottom left, and the foot as a shape of its own
        // so that it can say later which elements are steel. Both get the same turn and the
        // same lift, so the foot stays where it belongs on the bracket.
        final Outline flat = Outline.of(
                0.0, 0.0, SPAN, 0.0, SPAN, FOOT, LEG, FOOT, LEG, TALL, 0.0, TALL);
        final double[] pivot = flat.centre();
        final Outline turned = flat.rotated(tilt, pivot[0], pivot[1]);
        final double lift = -turned.bounds()[1];
        final Outline bracket = turned.translated(0.0, lift);
        final Outline steelFoot = Outline.rectangle(0.0, 0.0, SPAN, FOOT)
                .rotated(tilt, pivot[0], pivot[1]).translated(0.0, lift);

        final QuadMesh mesh = bracket.mesh(CELL);
        final boolean[] isSteel = new boolean[mesh.elementCount];
        for (int e = 0; e < mesh.elementCount; e++) {
            isSteel[e] = steelFoot.contains(mesh.centroidRadius(e), mesh.centroidZ(e));
        }
        final Materials materials = Materials.byElement(mesh.elementCount,
                e -> isSteel[e] ? STEEL : COPPER);

        final ExplicitSolver solver = new ExplicitSolver(mesh, materials,
                Formulation.PLANE_STRAIN, Integration.REDUCED, Kinematics.FINITE_STRAIN,
                DEPTH, 0.6);
        solver.setRigidWall(RigidWall.atZ(0.0, FRICTION));
        solver.setUniformVelocity(0.0, -speed);

        final double[] failureStrain = new double[mesh.elementCount];
        final double[] fractureEnergy = new double[mesh.elementCount];
        for (int e = 0; e < mesh.elementCount; e++) {
            failureStrain[e] = isSteel[e] ? 0.55 : 1.10;
            fractureEnergy[e] = isSteel[e] ? 60.0e3 : 120.0e3;
        }
        solver.setDamage(fractureEnergy, failureStrain);

        int steelCount = 0;
        for (boolean b : isSteel) if (b) steelCount++;

        System.out.println("NEOFIZ  plane-strain impact");
        System.out.println("=".repeat(76));
        System.out.printf(Locale.ROOT,
                "Shape     L bracket, %.0f x %.0f mm, %d elements at %.1f mm, %.0f mm deep%n",
                1e3 * SPAN, 1e3 * TALL, mesh.elementCount, 1e3 * CELL, 1e3 * DEPTH);
        System.out.printf(Locale.ROOT,
                "Material  %s in the foot (%d elements), %s in the leg (%d)%n",
                STEEL.name(), steelCount, COPPER.name(), mesh.elementCount - steelCount);
        System.out.printf(Locale.ROOT,
                "Impact    %.0f m/s onto a rigid anvil, tilted %.1f deg, friction %.1f%n",
                speed, Math.toDegrees(tilt), FRICTION);
        System.out.printf(Locale.ROOT, "Step      %.3f ns at the start%n%n",
                1e9 * solver.timestep());

        final Film film = new Film(mesh).planar()
                // Which substance an element is, so the two halves are legible before anything
                // has happened. The only field here that is a property of the scene rather
                // than of its state.
                .colourBy("material", "", 0.0, 1.0, (s, e) -> isSteel[e] ? 0.0 : 1.0)
                .colourBy("plastic strain", "", Film.AUTO, Film.AUTO,
                        (s, e) -> s.elementPlasticStrain(e))
                .colourBy("damage", "", 0.0, 1.0, (s, e) -> s.elementDamage(e))
                .colourBy("von Mises", "MPa", Film.AUTO, Film.AUTO,
                        (s, e) -> mises(s.centroidStress(e)) / 1e6)
                .readout("time, us", s -> 1e6 * s.time())
                .readout("kinetic energy, J", ExplicitSolver::kineticEnergy)
                .readout("plastic dissipation, J", ExplicitSolver::plasticDissipation)
                .readout("largest plastic strain", ExplicitSolver::maxPlasticStrain)
                .readout("worst damage", ExplicitSolver::maxDamage);

        // Run to a duration rather than a step count. The step shrinks as elements compress,
        // so a fixed number of steps would cover a different amount of time on every scene and
        // the film would end wherever the arithmetic happened to land.
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

        final double ke0 = 0.5 * massOf(mesh, materials) * speed * speed;
        System.out.printf(Locale.ROOT,
                "Energy    %.1f J in, %.1f J dissipated plastically, %.1f J still kinetic%n",
                ke0, solver.plasticDissipation(), solver.kineticEnergy());
        System.out.printf(Locale.ROOT,
                "Strain    largest %.3f, worst damage %.3f, %.1f %% of points yielded%n",
                solver.maxPlasticStrain(), solver.maxDamage(),
                100.0 * solver.yieldedFraction());
        // Against internal energy, not elastic strain energy. Once most of the mesh has
        // yielded the elastic store is a rounding error, and dividing by it would report a
        // catastrophe every time. The usual acceptance limit is a few per cent.
        final double internal = solver.strainEnergy() + solver.plasticDissipation()
                + solver.fractureDissipation();
        System.out.printf(Locale.ROOT,
                "Hourglass %.2f J, %.1f %% of the %.1f J of internal energy%n",
                solver.hourglassEnergy(), 100.0 * solver.hourglassEnergy() / internal, internal);
        System.out.printf(Locale.ROOT, "Run       %d steps, %.1f us simulated, %.1f s of wall "
                        + "clock%n%n",
                solver.steps(), 1e6 * solver.time(), wall);

        RunLog.of("impact")
                .parameter("speed", speed)
                .parameter("tilt", Math.toDegrees(tilt))
                .parameter("cell", CELL)
                .result("kinetic", solver.centredKineticEnergy())
                .result("plastic", solver.plasticDissipation())
                .result("maxPlasticStrain", solver.maxPlasticStrain())
                .result("hourglass", solver.hourglassEnergy())
                .result("steps", solver.steps())
                .result("wallClock", wall)
                .commit();

        final Path file = film.write("bracket.js",
                "Bracket on its corner, " + (int) speed + " m/s");
        System.out.printf(Locale.ROOT, "%d frames written to %s (%.1f MB)%n",
                film.size(), file, file.toFile().length() / (1024.0 * 1024.0));
        solver.shutdown();
    }

    /** Von Mises equivalent stress from the four stored components, Pa. */
    private static double mises(double[] s) {
        final double a = s[0] - s[1], b = s[1] - s[2], c = s[2] - s[0];
        return Math.sqrt(0.5 * (a * a + b * b + c * c) + 3.0 * s[3] * s[3]);
    }

    /** Mass of the section, kg, for the energy the bracket arrives with. */
    private static double massOf(QuadMesh mesh, Materials materials) {
        double m = 0.0;
        for (int e = 0; e < mesh.elementCount; e++) {
            m += materials.at(e).density() * mesh.signedArea(e) * DEPTH;
        }
        return m;
    }
}
