package org.neofiz;

import org.neofiz.render.Film;
import org.neofiz.validate.BurstCase;
import org.neofiz.validate.DefectBurst;

import java.nio.file.Path;
import java.util.Locale;

/**
 * A tube pressurised until it gives way, captured frame by frame so it can be watched.
 *
 * <p>Every number this project has produced so far came out of a run nobody has ever seen.
 * The solver has been tracking the position of every node through several hundred thousand
 * timesteps and discarding all of it except a peak pressure and a station along the axis. This
 * keeps the shape.
 *
 * <p>The case is {@link DefectBurst} rather than {@link BurstCase} because it has somewhere to
 * fail. {@code BurstCase} is a two-element-long slice of wall, which is the right harness for
 * measuring a pressure and has nothing to look at; a defect tube is forty millimetres long
 * with a mesh-independent flaw field in it, so the seed decides where the weak material is.
 *
 * <p>What the film then shows is <b>not</b> a neck: the bulge ratio on the reference tube is
 * about 1.009, so the wall is still swelling very nearly evenly when it gives way, and damage
 * runs from 0.77 to 1.00 across the <em>whole</em> wall rather than reaching one at one
 * station. That is the honest picture and it is the same fact from two other directions --
 * that a pressurised shell averages its defects over its shear-lag length rather than failing
 * at its weakest link, and that the burst pressure does not converge under refinement because
 * the load maximum arrives while the damage is still diffuse. A film that showed a tidy crack
 * would be showing something this model does not currently do.
 *
 * <pre>  ./gradlew burstFilm --args="&lt;seed&gt;"</pre>
 *
 * <p>Writes {@code runs/film.js}, which is a viewer's input rather than a report. The numbers
 * printed here are only enough to say which run is on the film.
 */
public final class BurstFilm {

    /** The tube modelled, in millimetres of length. Long enough to have somewhere to fail. */
    private static final double LENGTH = 40.0e-3;

    private BurstFilm() {
    }

    public static void main(String[] args) {
        final long seed = args.length > 0 ? Long.parseLong(args[0]) : 1L;
        final DefectBurst.Setup setup = DefectBurst.Setup.nominal().withLength(LENGTH);

        System.out.println("NEOFIZ  filming a burst");
        System.out.println("=".repeat(76));
        System.out.printf(Locale.ROOT,
                "Tube      OFHC Cu, mean radius %.1f mm, wall %.2f mm, %.0f mm long%n",
                1e3 * DefectBurst.MEAN_RADIUS, 1e3 * setup.thickness(), 1e3 * setup.length());
        System.out.printf(Locale.ROOT,
                "Mesh      %d through the wall, %d along the axis%n",
                setup.throughWall(), setup.alongAxis());
        System.out.printf(Locale.ROOT,
                "Defects   Weibull m = %.0f, scale %.3f, correlation %.1f mm, seed %d%n%n",
                DefectBurst.MODULUS, DefectBurst.SCALE, 1e3 * DefectBurst.CORRELATION, seed);

        // Built before the run so that the mesh the film indexes is the mesh the run uses.
        // The observer hands both back at every sample; the capture only needs the solver.
        final Film[] film = new Film[1];
        final double reference = DefectBurst.referencePressure(setup.slenderness());

        final DefectBurst.Shot shot = DefectBurst.fire(seed, setup, (solver, mesh, capacity) -> {
            if (film[0] == null) {
                film[0] = new Film(mesh)
                        // Damage keeps its absolute scale: 0 is intact and 1 is gone, and
                        // rescaling that to whatever this run reached would make an
                        // undamaged element look damaged. The other two are fitted to the
                        // run, because their useful range is a property of the tube rather
                        // than of the quantity.
                        .colourBy("plastic strain", "", Film.AUTO, Film.AUTO,
                                (s, e) -> s.elementPlasticStrain(e))
                        .colourBy("damage", "", 0.0, 1.0,
                                (s, e) -> s.elementDamage(e))
                        .colourBy("hoop stress", "MPa", Film.AUTO, Film.AUTO,
                                (s, e) -> s.centroidStress(e)[2] / 1e6)
                        .readout("time, ms", s -> 1e3 * s.time())
                        .readout("applied pressure, MPa", s -> s.borePressure() / 1e6)
                        // Recomputed per frame rather than taken from the observer's
                        // argument, so that the film carries the same quantity the M1 gate
                        // measures: the pressure the wall could hold in this state. It is
                        // the curve the viewer draws its scrub track from, and its maximum
                        // is the burst.
                        .readout("wall capacity, MPa",
                                s -> BurstCase.loadCapacity(s, mesh) / 1e6)
                        .readout("largest plastic strain", s -> s.maxPlasticStrain())
                        .readout("worst damage", s -> s.maxDamage());
            }
            film[0].capture(solver);
        });

        if (film[0] == null) throw new IllegalStateException("the run took no samples");

        System.out.printf(Locale.ROOT, "Burst     %.3f MPa, %.1f %% below the defect-free %.3f%n",
                shot.peakCapacity() / 1e6,
                100.0 * (1.0 - shot.peakCapacity() / reference), reference / 1e6);
        System.out.printf(Locale.ROOT,
                "Failed    %.1f mm along the tube, worst damage %.3f%n",
                1e3 * shot.failureZ(), shot.maxDamage());
        System.out.printf(Locale.ROOT,
                "Bulge     %.3f -- the widest station over the mean, so %s%n",
                shot.bulgeRatio(), shot.bulgeRatio() > 1.02
                        ? "it picked a place rather than swelling evenly"
                        : "it was still swelling evenly when the run ended");
        System.out.printf(Locale.ROOT, "Run       %d steps, %.1f s%n%n",
                shot.steps(), shot.wallClock());

        final Path file = film[0].write("film.js", "A copper tube, pressurised until it gives way");
        System.out.printf(Locale.ROOT, "%d frames written to %s%n", film[0].size(), file);
        System.out.printf(Locale.ROOT, "%.1f MB%n",
                file.toFile().length() / (1024.0 * 1024.0));
    }
}
