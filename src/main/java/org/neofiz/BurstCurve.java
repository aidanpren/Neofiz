package org.neofiz;

import org.neofiz.core.JohnsonCook;
import org.neofiz.report.Csv;
import org.neofiz.report.Plot;
import org.neofiz.validate.Burst;
import org.neofiz.validate.BurstCase;

import java.nio.file.Path;
import java.util.List;
import java.util.Locale;

/**
 * The pressure-expansion curve of a bursting tube, as an output rather than as one number
 * off the top of it.
 *
 * <p>The burst harness has always traced {@code P(eps_theta)} through the maximum and thrown
 * all of it away but the peak, because the peak is what the M1 gate is about. The curve is
 * the more interesting object. It is the charter's bulge-as-positive-feedback story in one
 * picture: pressure and expansion rise together while strain hardening is winning, flatten
 * as the wall thins and the radius grows, and turn over at the point where geometry starts
 * winning -- with nothing having broken, no failure criterion consulted, and no threshold
 * anywhere in the model.
 *
 * <p>Three curves are drawn on the same axes and the comparison between them is the report:
 *
 * <ul>
 *   <li><b>The wall's capacity</b>, from {@link BurstCase#loadCapacity} -- the pressure this
 *       tube can hold in its current state. Defined on both sides of the maximum, which is
 *       exactly why it is the thing measured.</li>
 *   <li><b>The applied bore pressure</b>, which tracks the capacity up the stable branch and
 *       then sits on its ceiling while the capacity falls away beneath it. The gap that
 *       opens between them past the peak is the unbalanced force accelerating the wall
 *       outward. It is the instability, drawn.</li>
 *   <li><b>The closed form</b> from {@link Burst#pressure}, rigid-plastic and thin-walled.
 *       It has no elastic strain in it, so it sits a little left of the measured curve and a
 *       little above it -- see {@link Burst#elasticState}, which predicts that offset to
 *       0.004 % from {@code E} and {@code nu} alone.</li>
 * </ul>
 *
 * <p>The terminal plot is a reading aid. The file this also writes is the curve at full
 * precision, for plotting properly.
 */
public final class BurstCurve {

    private BurstCurve() {
    }

    public static void main(String[] args) {
        final BurstCase.Result run = BurstCase.reference();
        final JohnsonCook law = BurstCase.COPPER.johnsonCook();

        System.out.println("NEOFIZ  the burst curve");
        System.out.println("=".repeat(76));
        System.out.printf(Locale.ROOT,
                "Material  %s on its quasi-static flow curve%n", BurstCase.COPPER.name());
        System.out.printf(Locale.ROOT,
                "Geometry  mean radius %.1f mm, wall %.2f mm, D/t = %.0f%n",
                1e3 * run.meanRadius(), 1e3 * run.thickness(), run.slenderness());
        System.out.printf(Locale.ROOT,
                "Mesh      %d elements through the wall, %d total, %s integration%n",
                run.elementsThroughWall(), run.elementCount(),
                run.integration().name().toLowerCase(Locale.ROOT));
        System.out.printf(Locale.ROOT,
                "Load      ramp to %.2f MPa, %.0f%% above the closed-form burst%n",
                run.appliedCeiling() / 1e6, 100.0 * (run.overshoot() - 1.0));
        System.out.printf(Locale.ROOT,
                "Run       %d steps, %d samples, %.1f s wall clock%n%n",
                run.steps(), run.trace().size(), run.wallClockSeconds());

        run.requirePeak();

        final List<BurstCase.Sample> trace = run.trace();
        final int peak = run.peakIndex();

        final double[] strain = new double[trace.size()];
        final double[] capacity = new double[trace.size()];
        final double[] applied = new double[trace.size()];
        for (int i = 0; i < trace.size(); i++) {
            strain[i] = 100.0 * trace.get(i).hoopStrain();
            capacity[i] = trace.get(i).capacity() / 1e6;
            applied[i] = trace.get(i).applied() / 1e6;
        }

        // Where the bore surface first flows. The elastic rise ends here, and everything the
        // curve does afterwards is the hardening-against-geometry competition.
        int firstYield = 0;
        while (firstYield < trace.size() - 1
                && trace.get(firstYield).maxPlasticStrain() <= 0.0) {
            firstYield++;
        }

        // The closed form over the same expansion, sampled finely enough that its own peak is
        // not an artefact of where the samples happened to fall.
        final int points = 200;
        final double[] exactStrain = new double[points];
        final double[] exactPressure = new double[points];
        for (int i = 0; i < points; i++) {
            final double hoop = strain[strain.length - 1] * i / (points - 1.0) / 100.0;
            exactStrain[i] = 100.0 * hoop;
            exactPressure[i] =
                    Burst.pressure(law, run.thickness(), run.meanRadius(), hoop) / 1e6;
        }

        // Drawn back to front: the closed form is the reference, the applied pressure is the
        // boundary condition, and the measured capacity is the answer, so the answer is the
        // one that survives where they overlap.
        System.out.println(Plot.standard()
                .title("A tube that has not broken, and cannot hold any more")
                .yLabel("pressure, MPa")
                .xLabel("mid-wall hoop strain, %")
                .series("closed form", '~', exactStrain, exactPressure)
                .series("applied", '.', strain, applied)
                .series("capacity", '#', strain, capacity)
                .marker("first yield", 'y', strain[firstYield], capacity[firstYield])
                .marker("burst", '*', strain[peak], capacity[peak])
                .render());
        System.out.println();

        // The turnover is a 0.5 % feature on a curve that spans 19 MPa, so on axes that show
        // the elastic rise it is a flat top and nothing more. This is the same curve with the
        // pressure axis opened up around the maximum, and it is where the instability is
        // actually legible: the capacity peaks and falls while the applied pressure, which
        // reached its ceiling long before, does not.
        // The run stops once the capacity has fallen PEAK_MARGIN below its maximum, so a
        // window of twice that either side of the peak holds the whole turnover and nothing
        // else. The closed form is left off these axes deliberately: the comparison has
        // already been made above, and what is worth looking at here is one curve going over
        // the top.
        final double window = 2.0 * BurstCase.PEAK_MARGIN;
        int from = peak;
        while (from > 0 && capacity[from] > (1.0 - window) * capacity[peak]) from--;
        System.out.println(Plot.of(64, 14)
                .title("The same maximum, with the pressure axis opened up")
                .yLabel("pressure, MPa")
                .xLabel("mid-wall hoop strain, %")
                .yRange((1.0 - window) * capacity[peak], 1.0002 * capacity[peak])
                .xRange(strain[from], strain[strain.length - 1])
                .series("capacity", '#', strain, capacity)
                .marker("burst", '*', strain[peak], capacity[peak])
                .render());
        System.out.println();

        System.out.println("-".repeat(76));
        System.out.printf(Locale.ROOT, "burst pressure, this run       %10.4f MPa   %+.3f %%%n",
                run.burstPressureFem() / 1e6, run.burstPressureErrorPercent());
        System.out.printf(Locale.ROOT, "burst pressure, closed form    %10.4f MPa%n",
                run.burstPressureExact() / 1e6);
        System.out.printf(Locale.ROOT, "hoop strain at burst           %10.4f %%    "
                        + "(closed form %.4f %%)%n",
                100.0 * run.burstHoopStrainFem(), 100.0 * run.burstHoopStrainExact());
        System.out.printf(Locale.ROOT, "Barlow on the mean diameter    %10.4f MPa   %+.2f %%%n",
                run.barlowPressure() / 1e6, run.barlowErrorPercent());
        System.out.println();
        System.out.printf(Locale.ROOT, "kinetic / strain energy at the peak        %.4f %%   "
                        + "(quasi-static)%n", 100.0 * run.kineticOverStrainAtPeak());
        System.out.printf(Locale.ROOT, "hourglass / plastic work                   %.4f %%%n",
                100.0 * run.hourglassOverPlasticWork());
        System.out.printf(Locale.ROOT, "ceiling above capacity at the peak         %.2f %%   "
                        + "(the creep imbalance)%n", 100.0 * run.creepImbalance());
        System.out.println();

        // A single run's error is not the gate figure and should not be read as one. The
        // ceiling has to sit above burst for a peak to exist at all, and sitting above burst
        // is what makes the wall creep, which biases the reading low in proportion to the
        // overshoot. Three ceilings and a straight line remove it; see
        // BurstCase.extrapolateToZeroOvershoot for why that extrapolation is exact rather
        // than fitted. It costs three more runs, which on this case is a few seconds.
        System.out.println("-".repeat(76));
        System.out.println("With the harness's own overshoot extrapolated away:");
        final BurstCase.OvershootFit fit = BurstCase.extrapolated();
        System.out.printf(Locale.ROOT, "burst pressure, extrapolated   %10.4f MPa   %+.3f %%%n",
                fit.burstPressure() / 1e6, fit.errorPercent());
        System.out.printf(Locale.ROOT, "  elasticity alone predicts                %+.3f %%%n",
                -run.elasticKnockdownPercent());
        System.out.printf(Locale.ROOT, "  leaving unaccounted for                  %+.3f %%%n",
                fit.errorPercent() + run.elasticKnockdownPercent());
        System.out.printf(Locale.ROOT, "  the fit's own linearity witness           %.5f %%%n",
                fit.linearityPercent());
        System.out.println();

        final Csv csv = Csv.of("time_s", "hoop_strain", "capacity_Pa", "applied_Pa",
                "kinetic_over_strain", "max_plastic_strain", "yielded_fraction");
        for (BurstCase.Sample s : trace) {
            csv.row(s.time(), s.hoopStrain(), s.capacity(), s.applied(),
                    s.kineticOverStrain(), s.maxPlasticStrain(), s.yieldedFraction());
        }
        final Path file = csv.write("burst-curve.csv");
        System.out.printf(Locale.ROOT, "%d samples written to %s%n", csv.size(), file);
    }
}
