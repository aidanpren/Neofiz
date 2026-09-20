package org.neofiz;

import org.neofiz.report.Csv;
import org.neofiz.report.Plot;
import org.neofiz.solver.Integration;
import org.neofiz.solver.Kinematics;
import org.neofiz.sweep.PowerLaw;
import org.neofiz.sweep.Sweep;
import org.neofiz.validate.Burst;
import org.neofiz.validate.BurstCase;

import java.nio.file.Path;
import java.util.List;
import java.util.Locale;

/**
 * Two sweeps, and the law that falls out of them.
 *
 * <p>The M3 exit criterion in the build plan is that a player rediscovers Barlow from a sweep
 * without being told. This is the machine half of that: vary one dimension at a time, measure
 * the burst pressure of each tube, and fit a power law to the result. Nothing here knows what
 * Barlow's formula is. {@link org.neofiz.validate.Burst#barlow} exists and is deliberately
 * not consulted -- the exponents below are measured off the finite element solver, and the
 * comparison against the rule of thumb happens at the end, after the numbers are in.
 *
 * <pre>
 *   wall thickness  -&gt;  burst pressure     a straight line        P proportional to t
 *   bore diameter   -&gt;  burst pressure     inverse relationship   P proportional to 1/d
 * </pre>
 *
 * <p>Together those two exponents are {@code P = 2 sigma t / D}, which is the entire content
 * of the formula, arrived at without writing it down. That is the whole claim the experiment
 * layer makes, on the smallest case that can support it.
 *
 * <h2>What the residual is doing in the report</h2>
 *
 * The thickness exponent is not exactly one and should not be. Barlow is a thin-wall
 * statement, and the sweep deliberately runs from {@code D/t = 10}, which is not thin, to
 * {@code D/t = 100}, which is. Where the fit stops describing the data is where the
 * assumption stops holding, and the report says where that is rather than choosing a range
 * that hides it.
 *
 * <h2>The bias that does not matter here</h2>
 *
 * Each point is a single run at {@link BurstCase#OVERSHOOT}, not the three-run extrapolation
 * the M1 gate uses, so every point carries the same 0.2 % creep bias low. A systematic factor
 * that is the same at every point moves the fitted <em>coefficient</em> and leaves the
 * <em>exponent</em> untouched, and the exponent is the answer. Paying three times the compute
 * to move a number the fit discards would be the wrong trade.
 */
public final class BarlowSweep {

    /** Points per sweep. Twelve over a decade is plenty to see a line and to see a bend. */
    private static final int POINTS = 12;

    /** Mean-diameter-over-thickness range: thick enough to break Barlow, thin enough to obey. */
    private static final double SLENDERNESS_LOW = 10.0;
    private static final double SLENDERNESS_HIGH = 100.0;

    /** Mean radii for the diameter sweep, metres. A factor of eight. */
    private static final double RADIUS_LOW = 12.5e-3;
    private static final double RADIUS_HIGH = 100.0e-3;

    /** Wall thickness held fixed while the diameter varies, metres. */
    private static final double FIXED_WALL = 2.0e-3;

    private BarlowSweep() {
    }

    public static void main(String[] args) {
        System.out.println("NEOFIZ  two sweeps");
        System.out.println("=".repeat(76));
        System.out.printf(Locale.ROOT,
                "Material  %s on its quasi-static flow curve%n", BurstCase.COPPER.name());
        System.out.printf(Locale.ROOT,
                "Sweeps    %d points each, %d threads, one thread per run%n%n",
                POINTS, Math.min(Sweep.DEFAULT_THREADS, POINTS));

        final long t0 = System.nanoTime();
        final PowerLaw thickness = thicknessSweep();
        final PowerLaw diameter = diameterSweep();
        final double elapsed = (System.nanoTime() - t0) * 1e-9;

        System.out.println("=".repeat(76));
        System.out.println("What the two sweeps say, together");
        System.out.println();
        System.out.printf(Locale.ROOT, "   %s%n",
                thickness.describe("burst pressure", "wall thickness"));
        System.out.printf(Locale.ROOT, "   %s%n",
                diameter.describe("burst pressure", "mean diameter"));
        System.out.println();
        System.out.println("   which is P = k t / D. The exponents are the shape of the law;");
        System.out.println("   the coefficient is the material in it.");
        System.out.println();

        // Closing the loop. The thickness sweep ran at a fixed mean diameter, so its fitted
        // coefficient is k sigma / D and the strength falls straight out of it. Nothing in
        // the sweep was told the flow curve: this is the solver's own answer for how strong
        // the copper is, recovered from twelve tubes of different thicknesses bursting.
        final double meanDiameter = 2e3 * BurstCase.MEAN_RADIUS;
        final double fromSweep = 0.5 * thickness.coefficient() * meanDiameter;
        final double effective = Burst.effectiveStress(BurstCase.COPPER.johnsonCook()) / 1e6;
        System.out.printf(Locale.ROOT,
                "   solving P = 2 sigma t / D for sigma gives      %8.2f MPa%n", fromSweep);
        System.out.printf(Locale.ROOT,
                "   the flow curve's own effective stress is       %8.2f MPa   (%+.2f %%)%n",
                effective, 100.0 * (fromSweep - effective) / effective);
        System.out.printf(Locale.ROOT,
                "   its engineering UTS, which a data sheet quotes %8.2f MPa%n",
                Burst.ultimateTensileStrength(BurstCase.COPPER.johnsonCook()) / 1e6);
        System.out.println();
        System.out.println("   The gap between the last two is why Barlow is conservative,");
        System.out.println("   and it is the only number here anybody had to be told.");
        System.out.printf(Locale.ROOT, "%n%d runs in %.1f s.%n", 2 * POINTS, elapsed);
    }

    // ---------------------------------------------------------------- wall thickness

    private static PowerLaw thicknessSweep() {
        final double[] slenderness =
                Sweep.logarithmic(SLENDERNESS_HIGH, SLENDERNESS_LOW, POINTS);
        final List<BurstCase.Result> runs =
                Sweep.over(slenderness, BurstCase::atSlenderness);

        final double[] wall = new double[POINTS];
        final double[] pressure = new double[POINTS];
        for (int i = 0; i < POINTS; i++) {
            runs.get(i).requirePeak();
            wall[i] = 1e3 * runs.get(i).thickness();
            pressure[i] = runs.get(i).burstPressureFem() / 1e6;
        }

        final PowerLaw law = PowerLaw.fit(wall, pressure);
        report("Sweep 1: wall thickness, with the mean radius held at 25 mm",
                "wall thickness, mm", wall, pressure, law, runs);

        // Barlow's thin-wall statement is an exact proportionality, so the departure from an
        // exponent of one is the thick-wall correction and nothing else. Naming its size is
        // more useful than asserting the fit is good.
        System.out.printf(Locale.ROOT,
                "   the exponent is %.4f, which is %+.2f %% away from the thin-wall 1%n%n",
                law.exponent(), 100.0 * (law.exponent() - 1.0));
        return law;
    }

    // ---------------------------------------------------------------- bore diameter

    private static PowerLaw diameterSweep() {
        final double[] radii = Sweep.logarithmic(RADIUS_LOW, RADIUS_HIGH, POINTS);
        final List<BurstCase.Result> runs = Sweep.over(radii, BarlowSweep::atRadius);

        final double[] bore = new double[POINTS];
        final double[] pressure = new double[POINTS];
        for (int i = 0; i < POINTS; i++) {
            runs.get(i).requirePeak();
            bore[i] = 1e3 * (2.0 * runs.get(i).meanRadius() - runs.get(i).thickness());
            pressure[i] = runs.get(i).burstPressureFem() / 1e6;
        }

        final PowerLaw law = PowerLaw.fit(bore, pressure);
        report("Sweep 2: bore diameter, with the wall held at 2 mm",
                "bore diameter, mm", bore, pressure, law, runs);
        System.out.printf(Locale.ROOT,
                "   the exponent is %.4f, which is %+.2f %% away from the thin-wall -1%n%n",
                law.exponent(), 100.0 * (law.exponent() + 1.0) / -1.0);

        // And here is what that residual was. The wall does not know what the bore diameter
        // is; equilibrium is written on the mean radius, so with the wall held fixed the
        // burst pressure goes as 1/(d + t) rather than as 1/d -- which is not a power law in
        // d at all, and is why the fit came back bent and short of -1. Replot the identical
        // runs against the mean diameter and the bend is gone.
        //
        // This is the same statement as Burst's note that Barlow written on the outside
        // diameter is the right expression with the wrong radius, and it is most of why the
        // rule of thumb is conservative. Nobody was told it here. It arrived as curvature in
        // a residual plot.
        final double[] mean = new double[POINTS];
        for (int i = 0; i < POINTS; i++) mean[i] = 2e3 * runs.get(i).meanRadius();
        final PowerLaw onMean = PowerLaw.fit(mean, pressure);

        System.out.println("   The same runs, against the mean diameter instead of the bore:");
        System.out.printf(Locale.ROOT, "   %s%n", onMean);
        System.out.printf(Locale.ROOT,
                "   the exponent is %.4f, and the worst residual fell from %.3f %% to "
                        + "%.3f %%.%n", onMean.exponent(), law.maxResidualPercent(),
                onMean.maxResidualPercent());
        System.out.println("   Which diameter Barlow is written on is not a convention.");
        System.out.println();
        return onMean;
    }

    /** The reference tube at a chosen mean radius, with the wall thickness held fixed. */
    private static BurstCase.Result atRadius(double meanRadius) {
        return BurstCase.run(BurstCase.COPPER, meanRadius, 2.0 * meanRadius / FIXED_WALL,
                BurstCase.ELEMENTS_THROUGH_WALL, Integration.REDUCED,
                Kinematics.FINITE_STRAIN, BurstCase.OVERSHOOT);
    }

    // ---------------------------------------------------------------- reporting

    private static void report(String title, String xLabel, double[] x, double[] y,
                               PowerLaw law, List<BurstCase.Result> runs) {
        System.out.println("-".repeat(76));
        System.out.println(title);
        System.out.println();

        System.out.printf(Locale.ROOT, "   %12s %14s %14s %12s%n",
                xLabel, "burst, MPa", "closed form", "difference");
        for (int i = 0; i < x.length; i++) {
            System.out.printf(Locale.ROOT, "   %12.4f %14.4f %14.4f %11.3f %%%n",
                    x[i], y[i], runs.get(i).burstPressureExact() / 1e6,
                    runs.get(i).burstPressureErrorPercent());
        }
        System.out.println();

        // A power law is a straight line in log-log and a curve in anything else, so this is
        // the only pair of axes on which the eye can check the fit rather than take it.
        final double[] logX = new double[x.length];
        final double[] logY = new double[y.length];
        final double[] logFit = new double[x.length];
        for (int i = 0; i < x.length; i++) {
            logX[i] = Math.log10(x[i]);
            logY[i] = Math.log10(y[i]);
            logFit[i] = Math.log10(law.at(x[i]));
        }

        System.out.println(Plot.of(62, 16)
                .yLabel("log10( burst pressure / MPa )")
                .xLabel("log10( " + xLabel.replace(", mm", " / mm") + " )")
                .series("fitted power law", '-', logX, logFit)
                .series("measured", '#', logX, logY)
                .render());
        System.out.println();
        System.out.printf(Locale.ROOT, "   %s%n", law);
        System.out.println();

        // On log-log axes spanning a decade, a 1 % departure from the fit is a fraction of a
        // character cell and the eye reports a perfect straight line whether or not there is
        // one. The residuals are where a power law that is nearly right differs from one that
        // is right, and "nearly right, with structure" is the more interesting outcome: it
        // means the exponent is flattening something out.
        final double[] residual = new double[x.length];
        final double[] zero = new double[x.length];
        for (int i = 0; i < x.length; i++) {
            residual[i] = 100.0 * (y[i] - law.at(x[i])) / law.at(x[i]);
        }
        System.out.println(Plot.of(62, 9)
                .yLabel("departure from the fitted law, %")
                .xLabel("log10( " + xLabel.replace(", mm", " / mm") + " )")
                .series("zero", '.', logX, zero)
                .series("residual", '#', logX, residual)
                .render());

        final Csv csv = Csv.of(xLabel.replace(", ", "_").replace(" ", "_"),
                "burst_Pa", "closed_form_Pa", "error_percent", "hoop_strain_at_burst");
        for (int i = 0; i < x.length; i++) {
            csv.row(x[i], runs.get(i).burstPressureFem(), runs.get(i).burstPressureExact(),
                    runs.get(i).burstPressureErrorPercent(), runs.get(i).burstHoopStrainFem());
        }
        final Path file = csv.write(xLabel.startsWith("wall")
                ? "sweep-wall-thickness.csv" : "sweep-bore-diameter.csv");
        System.out.printf(Locale.ROOT, "   %d points written to %s%n%n", csv.size(), file);
    }
}
