package org.neofiz.report;

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * A line chart drawn in characters, for showing a curve in the same terminal report
 * everything else in this project prints to.
 *
 * <h2>Why this exists at all</h2>
 *
 * The build plan's case for the experiment layer is that people do not extract physical law
 * from single observations, they extract it from contrast -- and that the project's numbers
 * only start teaching once they are curves rather than table rows. This is the cheapest
 * possible version of that: no window, no dependency, no graphics stack, and available from
 * the first run. It is deliberately not the product's plotting; it is what keeps a curve
 * visible while the product's plotting does not exist yet.
 *
 * <p>Plain ASCII glyphs, not box-drawing characters, because this has to survive whatever
 * code page a Windows console happens to be in. The pretty version of any of these plots
 * comes from the CSV that {@link Csv} writes alongside.
 *
 * <h2>Units</h2>
 *
 * The plot is unit-agnostic and does no scaling: pass numbers already in the units the axis
 * label claims. Pascals plotted against a label saying MPa is the caller's bug, and hiding a
 * factor of a million inside a chart library is how that bug gets made.
 */
public final class Plot {

    /** Columns of plotting area, excluding the y-axis gutter. */
    private final int width;

    /** Rows of plotting area, excluding the x-axis labels. */
    private final int height;

    private static final int GUTTER = 11;

    private String title = "";
    private String xLabel = "";
    private String yLabel = "";

    private final List<Series> series = new ArrayList<>();
    private final List<Marker> markers = new ArrayList<>();

    private double xMin = Double.NaN;
    private double xMax = Double.NaN;
    private double yMin = Double.NaN;
    private double yMax = Double.NaN;

    private record Series(String name, char glyph, double[] x, double[] y) {
    }

    private record Marker(String name, char glyph, double x, double y) {
    }

    private Plot(int width, int height) {
        if (width < 20 || height < 5) {
            throw new IllegalArgumentException("a plot smaller than 20x5 is not a plot");
        }
        this.width = width;
        this.height = height;
    }

    /** A plot of the given interior size, in characters. */
    public static Plot of(int width, int height) {
        return new Plot(width, height);
    }

    /** A plot sized for an 80-column terminal. */
    public static Plot standard() {
        return new Plot(64, 20);
    }

    public Plot title(String text) {
        this.title = text;
        return this;
    }

    public Plot xLabel(String text) {
        this.xLabel = text;
        return this;
    }

    public Plot yLabel(String text) {
        this.yLabel = text;
        return this;
    }

    /**
     * A curve. Points are joined, so a series with gaps in it draws lines across them --
     * split it into two series instead.
     */
    public Plot series(String name, char glyph, double[] x, double[] y) {
        if (x.length != y.length) {
            throw new IllegalArgumentException("x and y differ in length: "
                    + x.length + " and " + y.length);
        }
        series.add(new Series(name, glyph, x.clone(), y.clone()));
        return this;
    }

    /** A single labelled point, drawn over the curves. */
    public Plot marker(String name, char glyph, double x, double y) {
        markers.add(new Marker(name, glyph, x, y));
        return this;
    }

    /** Fixes the x range instead of taking it from the data. */
    public Plot xRange(double min, double max) {
        this.xMin = min;
        this.xMax = max;
        return this;
    }

    /** Fixes the y range instead of taking it from the data. */
    public Plot yRange(double min, double max) {
        this.yMin = min;
        this.yMax = max;
        return this;
    }

    /**
     * The chart, as lines of text with no trailing newline.
     *
     * <p>Axis ranges come from the data unless {@link #xRange} or {@link #yRange} fixed them.
     * A degenerate range -- one point, or a series that never moves -- is padded rather than
     * divided by, so a flat line renders as a flat line instead of a crash.
     */
    public String render() {
        double x0 = xMin;
        double x1 = xMax;
        double y0 = yMin;
        double y1 = yMax;

        if (Double.isNaN(x0) || Double.isNaN(x1)) {
            final double[] bounds = dataBounds(true);
            x0 = bounds[0];
            x1 = bounds[1];
        }
        if (Double.isNaN(y0) || Double.isNaN(y1)) {
            final double[] bounds = dataBounds(false);
            y0 = bounds[0];
            y1 = bounds[1];
        }
        x1 = padded(x0, x1);
        y1 = padded(y0, y1);

        final char[][] grid = new char[height][width];
        for (char[] row : grid) java.util.Arrays.fill(row, ' ');

        for (Series s : series) {
            int previousColumn = Integer.MIN_VALUE;
            int previousRow = Integer.MIN_VALUE;
            for (int i = 0; i < s.x().length; i++) {
                final int column = column(s.x()[i], x0, x1);
                final int row = row(s.y()[i], y0, y1);
                if (previousColumn != Integer.MIN_VALUE) {
                    line(grid, previousColumn, previousRow, column, row, s.glyph());
                } else {
                    put(grid, column, row, s.glyph());
                }
                previousColumn = column;
                previousRow = row;
            }
        }
        for (Marker m : markers) {
            put(grid, column(m.x(), x0, x1), row(m.y(), y0, y1), m.glyph());
        }

        final StringBuilder out = new StringBuilder();
        if (!title.isEmpty()) {
            out.append(" ".repeat(GUTTER)).append(title).append('\n');
        }
        if (!yLabel.isEmpty()) {
            out.append(" ".repeat(GUTTER)).append(yLabel).append('\n');
        }

        for (int r = 0; r < height; r++) {
            // Label the top and bottom rows and roughly every fifth in between, so the axis
            // stays readable at any height without the labels colliding.
            final boolean labelled = r == 0 || r == height - 1 || (height - 1 - r) % 5 == 0;
            final double value = y1 - (y1 - y0) * r / (height - 1.0);
            out.append(labelled
                    ? String.format(Locale.ROOT, "%9s |", axisFormat(value, y0, y1))
                    : " ".repeat(9) + " |");
            out.append(new String(grid[r])).append('\n');
        }

        out.append(" ".repeat(9)).append(" +").append("-".repeat(width)).append('\n');
        out.append(xAxisLabels(x0, x1)).append('\n');
        if (!xLabel.isEmpty()) {
            out.append(" ".repeat(GUTTER + Math.max(0, (width - xLabel.length()) / 2)))
                    .append(xLabel).append('\n');
        }

        final String legend = legend();
        if (!legend.isEmpty()) {
            out.append(legend).append('\n');
        }
        return out.substring(0, out.length() - 1);
    }

    private String legend() {
        final StringBuilder out = new StringBuilder();
        for (Series s : series) {
            if (out.length() > 0) out.append("   ");
            out.append(s.glyph()).append(' ').append(s.name());
        }
        for (Marker m : markers) {
            if (out.length() > 0) out.append("   ");
            out.append(m.glyph()).append(' ').append(m.name());
        }
        return out.length() == 0 ? "" : " ".repeat(GUTTER) + out;
    }

    private String xAxisLabels(double x0, double x1) {
        final char[] axis = new char[GUTTER + width];
        java.util.Arrays.fill(axis, ' ');
        final int ticks = Math.max(2, Math.min(6, width / 12));
        for (int t = 0; t < ticks; t++) {
            final double value = x0 + (x1 - x0) * t / (ticks - 1.0);
            final String text = axisFormat(value, x0, x1);
            // The last tick is right-aligned on the axis end and the rest left-aligned on
            // their own column, which is the arrangement that does not run off either edge.
            int at = GUTTER + (int) Math.round((width - 1.0) * t / (ticks - 1.0));
            if (t == ticks - 1) at -= text.length() - 1;
            else if (t > 0) at -= text.length() / 2;
            for (int k = 0; k < text.length() && at + k < axis.length; k++) {
                if (at + k >= 0) axis[at + k] = text.charAt(k);
            }
        }
        return new String(axis);
    }

    /**
     * Enough significant figures to tell neighbouring tick labels apart, which on a curve
     * whose interesting part is a 2 % window is not the same as "a sensible number of
     * decimals".
     */
    private static String axisFormat(double value, double low, double high) {
        final double span = Math.abs(high - low);
        if (span == 0.0) return String.format(Locale.ROOT, "%.3g", value);
        final double magnitude = Math.max(Math.abs(low), Math.abs(high));
        if (magnitude >= 1e5 || (magnitude > 0.0 && magnitude < 1e-3)) {
            return String.format(Locale.ROOT, "%.2e", value);
        }
        final int decimals = Math.max(0, Math.min(6,
                (int) Math.ceil(-Math.log10(span / 8.0)) + 1));
        return String.format(Locale.ROOT, "%." + decimals + "f", value);
    }

    private double[] dataBounds(boolean horizontal) {
        double low = Double.POSITIVE_INFINITY;
        double high = Double.NEGATIVE_INFINITY;
        for (Series s : series) {
            for (double v : horizontal ? s.x() : s.y()) {
                if (!Double.isFinite(v)) continue;
                low = Math.min(low, v);
                high = Math.max(high, v);
            }
        }
        for (Marker m : markers) {
            final double v = horizontal ? m.x() : m.y();
            if (!Double.isFinite(v)) continue;
            low = Math.min(low, v);
            high = Math.max(high, v);
        }
        if (low > high) throw new IllegalStateException("nothing finite to plot");
        return new double[] {low, high};
    }

    /** A range with no width in it is given some, so the mapping stays invertible. */
    private static double padded(double low, double high) {
        if (high > low) return high;
        return low == 0.0 ? 1.0 : low + Math.abs(low) * 1e-9;
    }

    private int column(double x, double x0, double x1) {
        return (int) Math.round((width - 1.0) * (x - x0) / (x1 - x0));
    }

    private int row(double y, double y0, double y1) {
        return (int) Math.round((height - 1.0) * (y1 - y) / (y1 - y0));
    }

    private void put(char[][] grid, int column, int row, char glyph) {
        if (row < 0 || row >= height || column < 0 || column >= width) return;
        grid[row][column] = glyph;
    }

    /** Bresenham, so that a steep segment is a line rather than two dots. */
    private void line(char[][] grid, int c0, int r0, int c1, int r1, char glyph) {
        final int dc = Math.abs(c1 - c0);
        final int dr = -Math.abs(r1 - r0);
        final int stepC = c0 < c1 ? 1 : -1;
        final int stepR = r0 < r1 ? 1 : -1;
        int error = dc + dr;
        int c = c0;
        int r = r0;
        while (true) {
            put(grid, c, r, glyph);
            if (c == c1 && r == r1) return;
            final int doubled = 2 * error;
            if (doubled >= dr) {
                error += dr;
                c += stepC;
            }
            if (doubled <= dc) {
                error += dc;
                r += stepR;
            }
        }
    }
}
