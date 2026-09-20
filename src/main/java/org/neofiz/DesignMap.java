package org.neofiz;

import org.neofiz.report.Csv;
import org.neofiz.report.Heatmap;
import org.neofiz.report.Runs;
import org.neofiz.sweep.Sweep;
import org.neofiz.validate.DefectBurst;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;

/**
 * The design map: wall thickness against applied pressure, every cell coloured by how often
 * the tube gives way.
 *
 * <p>The build plan calls this the single most valuable image in the product, and the reason
 * is not the map -- it is the <b>band</b>. Run one tube per cell and the boundary between
 * "holds" and "bursts" is a line, which is what a textbook draws and what a safety factor
 * pretends to be protecting against. Run several nominally identical tubes per cell, with the
 * defect field of {@link DefectBurst} supplying the only difference between them, and the
 * boundary comes out as a band of finite width. That band is the entire reason safety factors
 * exist, rendered as an image, derived from the player's own hardware, requiring no writing at
 * all.
 *
 * <h2>The pressure axis is free</h2>
 *
 * A map of {@code n} thicknesses by {@code m} pressures looks like {@code n * m} runs and is
 * not. A tube's burst pressure is a property of the tube: it is the maximum of the wall's own
 * load capacity, measured off an equilibrium identity, and it does not depend on what ceiling
 * the harness used to find it -- which is the point of {@link DefectBurst} measuring capacity
 * rather than applied pressure, and is checked by the overshoot sweep in
 * {@code BurstCase.extrapolateToZeroOvershoot}. So a tube bursts under an applied pressure
 * {@code P} exactly when {@code P} exceeds its capacity, and one population per thickness
 * answers every pressure at once.
 *
 * <p>The map therefore costs {@code thicknesses x samples} runs, and the pressure resolution
 * is limited by nothing but how many rows fit on a terminal. That is worth saying plainly
 * because the obvious implementation is two orders of magnitude more expensive and produces
 * the same picture.
 *
 * <h2>Arguments</h2>
 *
 * <pre>  ./gradlew designMap --args="&lt;thicknesses&gt; &lt;samples per thickness&gt;"</pre>
 *
 * <p>Defaults are deliberately modest. A shot here is a 160-element tube taken past its peak
 * and then sixty more ring periods to let the damage pick a station, so it is a good deal more
 * expensive than the 16-element tube {@code BurstCase} bursts, and the report prints what it
 * actually cost.
 */
public final class DesignMap {

    /** Thickness columns. */
    private static final int THICKNESSES = 8;

    /** Nominally identical tubes per column. Two is a smoke test; eight is a band. */
    private static final int SAMPLES = 6;

    /** Pressure rows. Free, so there are enough of them to resolve the band. */
    private static final int PRESSURE_ROWS = 25;

    /** Mean diameter over wall thickness, thick end and thin end. */
    private static final double SLENDERNESS_THICK = 12.5;
    private static final double SLENDERNESS_THIN = 40.0;

    private DesignMap() {
    }

    /** One tube: which column it belongs to, and which realisation of the defect field. */
    private record Tube(int column, double slenderness, long seed) {
    }

    public static void main(String[] args) {
        final int columns = args.length > 0 ? Integer.parseInt(args[0]) : THICKNESSES;
        final int samples = args.length > 1 ? Integer.parseInt(args[1]) : SAMPLES;
        if (columns < 2 || samples < 1) {
            throw new IllegalArgumentException("need at least two thicknesses and one sample");
        }

        final double[] slenderness =
                Sweep.logarithmic(SLENDERNESS_THIN, SLENDERNESS_THICK, columns);
        final double[] wall = new double[columns];
        for (int c = 0; c < columns; c++) {
            wall[c] = 2e3 * DefectBurst.MEAN_RADIUS / slenderness[c];
        }

        System.out.println("NEOFIZ  the design map");
        System.out.println("=".repeat(76));
        System.out.printf(Locale.ROOT,
                "Tube      %s, mean radius %.1f mm, %.0f mm of it modelled%n",
                "OFHC Cu", 1e3 * DefectBurst.MEAN_RADIUS, 1e3 * DefectBurst.TUBE_LENGTH);
        System.out.printf(Locale.ROOT,
                "Defects   Weibull m = %.0f, scale %.3f, correlation length %.1f mm%n",
                DefectBurst.MODULUS, DefectBurst.SCALE, 1e3 * DefectBurst.CORRELATION);
        System.out.printf(Locale.ROOT,
                "Map       %d wall thicknesses from %.2f to %.2f mm, %d tubes each%n",
                columns, wall[0], wall[columns - 1], samples);
        System.out.printf(Locale.ROOT,
                "          %d pressure rows, which cost nothing -- see the class note%n",
                PRESSURE_ROWS);
        System.out.printf(Locale.ROOT, "Runs      %d, on %d threads%n%n",
                columns * samples, Math.min(Sweep.DEFAULT_THREADS, columns * samples));

        // One flat list of every tube to be fired, so the whole map is a single ensemble and
        // the machine stays busy to the last run. Nesting a sweep per column inside a sweep
        // over columns would idle most of the cores on the last column.
        final List<Tube> tubes = new ArrayList<>(columns * samples);
        for (int c = 0; c < columns; c++) {
            for (int s = 0; s < samples; s++) tubes.add(new Tube(c, slenderness[c], s + 1L));
        }

        final long t0 = System.nanoTime();
        final List<DefectBurst.Shot> shots = Sweep.over(tubes, DesignMap::fire);
        final double elapsed = (System.nanoTime() - t0) * 1e-9;

        final double[][] capacity = new double[columns][samples];
        final int[] filled = new int[columns];
        int refused = 0;
        for (int i = 0; i < tubes.size(); i++) {
            final DefectBurst.Shot shot = shots.get(i);
            if (!shot.traversed()) {
                refused++;
                continue;
            }
            final int c = tubes.get(i).column();
            capacity[c][filled[c]++] = shot.peakCapacity();
        }
        if (refused > 0) {
            // A tube that never turned over has no burst pressure to place on the map, and
            // quietly treating it as "held" would put it on the safe side of the boundary
            // for a reason that is about the harness rather than the tube.
            System.out.printf(Locale.ROOT,
                    "%d of %d tubes did not traverse a peak and are excluded.%n%n",
                    refused, tubes.size());
        }
        for (int c = 0; c < columns; c++) {
            if (filled[c] == 0) {
                throw new IllegalStateException("no tube at " + wall[c] + " mm burst");
            }
            capacity[c] = Arrays.copyOf(capacity[c], filled[c]);
            Arrays.sort(capacity[c]);
        }

        absoluteMap(wall, capacity);
        normalisedMap(wall, slenderness, capacity);
        columnTable(wall, slenderness, capacity);
        write(wall, slenderness, capacity);

        System.out.printf(Locale.ROOT, "%n%d runs in %.1f s, %.1f s of solver per tube.%n",
                tubes.size(), elapsed, shots.stream()
                        .mapToDouble(DefectBurst.Shot::wallClock).average().orElse(0.0));
    }

    private static DefectBurst.Shot fire(Tube tube) {
        return DefectBurst.fire(tube.seed(),
                DefectBurst.Setup.nominal().withSlenderness(tube.slenderness()));
    }

    // ---------------------------------------------------------------- the map

    /** The map in engineering units: what pressure a tube of this wall will take. */
    private static void absoluteMap(double[] wall, double[][] capacity) {
        double low = Double.MAX_VALUE;
        double high = 0.0;
        for (double[] column : capacity) {
            low = Math.min(low, column[0]);
            high = Math.max(high, column[column.length - 1]);
        }
        // A little air either side, so the all-holds and all-bursts regions are visibly
        // regions rather than the edges of the picture.
        final double span = high - low;
        low -= 0.12 * span;
        high += 0.12 * span;

        final double[] pressure = new double[PRESSURE_ROWS];
        final double[][] fraction = new double[PRESSURE_ROWS][wall.length];
        for (int r = 0; r < PRESSURE_ROWS; r++) {
            pressure[r] = (low + (high - low) * r / (PRESSURE_ROWS - 1.0)) / 1e6;
            for (int c = 0; c < wall.length; c++) {
                fraction[r][c] = burstFraction(capacity[c], pressure[r] * 1e6);
            }
        }

        System.out.println(Heatmap.of(wall, pressure, fraction)
                .title("Where a tube of this wall holds this pressure, and where it does not")
                .yLabel("applied pressure, MPa")
                .xLabel("wall thickness, mm")
                .outcomes("all hold", "all burst")
                .cellWidth(6)
                .render());
        System.out.println();

        // Say what the picture cannot. On axes that span the whole design range the boundary
        // is a clean staircase, and it would be easy to read that as the band not existing.
        final double rowHeight = pressure[1] - pressure[0];
        double widest = 0.0;
        for (double[] column : capacity) {
            widest = Math.max(widest, column[column.length - 1] - column[0]);
        }
        System.out.printf(Locale.ROOT,
                "   The boundary looks like a line here, and is not one. The widest column%n"
                        + "   spans %.3f MPa from its weakest tube to its strongest, against a%n"
                        + "   row height of %.3f MPa -- so the band is narrower than a pixel on%n"
                        + "   these axes. The next map is the same data with that resolved.%n",
                widest / 1e6, rowHeight);
        System.out.println();
    }

    /**
     * The same runs with the pressure axis divided through by each column's own defect-free
     * burst pressure.
     *
     * <p>The knockdown and the scatter are both roughly constant <em>fractions</em>, so in
     * absolute pressure the band is a sliver at the bottom of the map and a wedge at the top,
     * and on any axis wide enough to hold the whole design range it is invisible. Divide
     * through and it becomes a horizontal band of a single width, which is the shape of the
     * statement it actually is: derate by this much, and be this unsure about it.
     */
    private static void normalisedMap(double[] wall, double[] slenderness,
                                      double[][] capacity) {
        final double[][] ratio = new double[wall.length][];
        double low = Double.MAX_VALUE;
        double high = 0.0;
        for (int c = 0; c < wall.length; c++) {
            final double reference = DefectBurst.referencePressure(slenderness[c]);
            ratio[c] = new double[capacity[c].length];
            for (int s = 0; s < capacity[c].length; s++) {
                ratio[c][s] = capacity[c][s] / reference;
            }
            low = Math.min(low, ratio[c][0]);
            high = Math.max(high, ratio[c][ratio[c].length - 1]);
        }
        final double span = Math.max(high - low, 1e-6);
        low -= 0.25 * span;
        high += 0.25 * span;

        final double[] level = new double[PRESSURE_ROWS];
        final double[][] fraction = new double[PRESSURE_ROWS][wall.length];
        for (int r = 0; r < PRESSURE_ROWS; r++) {
            level[r] = low + (high - low) * r / (PRESSURE_ROWS - 1.0);
            for (int c = 0; c < wall.length; c++) {
                fraction[r][c] = burstFraction(ratio[c], level[r]);
            }
        }

        System.out.println(Heatmap.of(wall, level, fraction)
                .title("The same tubes, against their own defect-free burst pressure")
                .yLabel("applied pressure / defect-free burst pressure")
                .xLabel("wall thickness, mm")
                .outcomes("all hold", "all burst")
                .cellWidth(6)
                .render());
        System.out.println();
        System.out.println("   There is the band. Nearly horizontal, which says the knockdown");
        System.out.println("   and the scatter are both fixed fractions rather than fixed");
        System.out.println("   pressures -- a safety factor, in the literal sense, derived from");
        System.out.println("   tubes rather than looked up. Nothing was told to produce it: the");
        System.out.println("   only difference between the tubes in a column is which draw of");
        System.out.println("   the same defect field they got.");
        System.out.println();
        System.out.println("   It is narrow, and that is the M2 result showing up again: along");
        System.out.println("   its axis a tube is a bundle of parallel rings rather than a");
        System.out.println("   chain of links, so a pressurised shell averages its defects");
        System.out.println("   instead of failing at its weakest one.");
        System.out.println();

        // And it slopes, which is the more interesting half. A thicker wall on the same mean
        // radius is proportionally more material, more material is more chances at a bad
        // patch, and that is a size effect. But it is a far smaller one than weakest-link
        // theory asks for, and the gap is the point rather than a discrepancy: it is the M2
        // finding measured a second way, on a different experiment, without being looked for.
        final double first = mean(ratio[0]);
        final double last = mean(ratio[wall.length - 1]);
        final double volumeRatio = wall[wall.length - 1] / wall[0];
        final double weakestLink = Math.pow(volumeRatio, -1.0 / DefectBurst.MODULUS);
        final double observed = last / first;

        System.out.printf(Locale.ROOT,
                "   And it slopes: %.2f %% knockdown at %.2f mm of wall against %.2f %% at%n"
                        + "   %.2f mm, over a factor of %.1f in material. A size effect --%n"
                        + "   more material, more chances at a bad patch.%n",
                100.0 * (first - 1.0), wall[0], 100.0 * (last - 1.0),
                wall[wall.length - 1], volumeRatio);
        System.out.println();
        System.out.printf(Locale.ROOT,
                "   Weakest-link theory says V^(-1/m) over that range, which for m = %.0f is%n"
                        + "   a %.2f %% fall. The measured fall is %.2f %%, smaller by a factor%n"
                        + "   of %.0f. The size effect is real and it is nothing like as strong%n"
                        + "   as a chain of links would make it -- which is the M2 result,%n"
                        + "   arrived at from a different experiment: a pressurised shell%n"
                        + "   averages its defects over the shear-lag length rather than%n"
                        + "   failing at its weakest one.%n",
                DefectBurst.MODULUS, 100.0 * (1.0 - weakestLink), 100.0 * (1.0 - observed),
                (1.0 - weakestLink) / (1.0 - observed));
        System.out.println();
    }

    /**
     * Fraction of a sorted population whose capacity is below the applied pressure.
     *
     * <p>The empirical distribution function of the batch, which is what a map cell is: the
     * proportion of nominally identical tubes that would not have survived this pressure.
     * Strictly below, so a pressure equal to a tube's capacity counts as held -- the tube
     * reached that pressure, which is the definition of its capacity.
     */
    static double burstFraction(double[] sorted, double applied) {
        int below = 0;
        while (below < sorted.length && sorted[below] < applied) below++;
        return (double) below / sorted.length;
    }

    // ---------------------------------------------------------------- the numbers

    private static void columnTable(double[] wall, double[] slenderness, double[][] capacity) {
        System.out.println("-".repeat(76));
        System.out.printf(Locale.ROOT, "   %9s %6s %10s %10s %10s %9s %9s%n",
                "wall, mm", "D/t", "lowest", "mean", "defect-free", "knockdown", "scatter");
        for (int c = 0; c < wall.length; c++) {
            final double mean = mean(capacity[c]);
            final double reference = DefectBurst.referencePressure(slenderness[c]);
            System.out.printf(Locale.ROOT,
                    "   %9.3f %6.1f %8.3f M %8.3f M %8.3f M %8.2f %% %8.2f %%%n",
                    wall[c], slenderness[c], capacity[c][0] / 1e6, mean / 1e6,
                    reference / 1e6, 100.0 * (mean / reference - 1.0),
                    100.0 * sd(capacity[c], mean) / mean);
        }
        System.out.println();
        System.out.println("   'lowest' is the weakest tube of the batch and 'defect-free' is");
        System.out.println("   the closed form for the same tube with no defects in it. The");
        System.out.println("   knockdown column is what a defect field costs you; the scatter");
        System.out.println("   column is what it costs you in confidence. Both are close to");
        System.out.println("   constant across a factor of three in wall thickness, which is");
        System.out.println("   what makes the normalised map flat and a derating rule possible.");
    }

    private static double mean(double[] a) {
        double m = 0.0;
        for (double x : a) m += x;
        return m / a.length;
    }

    private static double sd(double[] a, double mean) {
        if (a.length < 2) return 0.0;
        double v = 0.0;
        for (double x : a) v += (x - mean) * (x - mean);
        return Math.sqrt(v / (a.length - 1));
    }

    // ---------------------------------------------------------------- the file

    private static void write(double[] wall, double[] slenderness, double[][] capacity) {
        final Csv csv = Csv.of("wall_mm", "slenderness", "sample", "capacity_Pa",
                "defect_free_Pa");
        for (int c = 0; c < wall.length; c++) {
            for (int s = 0; s < capacity[c].length; s++) {
                csv.row(wall[c], slenderness[c], s, capacity[c][s],
                        DefectBurst.referencePressure(slenderness[c]));
            }
        }
        final Path file = csv.write("design-map.csv");
        System.out.printf(Locale.ROOT, "%n   %d tubes written to %s   (%s)%n",
                csv.size(), file, Runs.dir().toAbsolutePath());
    }
}
