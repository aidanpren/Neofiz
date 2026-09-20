package org.neofiz.report;

import java.util.Locale;

/**
 * A two-parameter map drawn in characters, one cell per pair of parameter values.
 *
 * <p>This is the shape of the build plan's design map: two dimensions swept at once, each
 * cell coloured by outcome. What makes it worth drawing rather than tabulating is that the
 * boundary between outcomes is the interesting object, and a boundary is a thing the eye
 * finds instantly in a picture and never finds in a table of numbers.
 *
 * <p>Values are fractions in {@code [0, 1]} and the caller does the normalising, because what
 * zero and one mean is the caller's business -- here it is the fraction of nominally
 * identical tubes that burst, which is a probability and needs no scaling, but a map of burst
 * pressure or of muzzle velocity would need its own.
 *
 * <h2>The ramp</h2>
 *
 * Eight steps from blank to {@code #}, in increasing visual weight. Blank for zero is
 * deliberate: on a map whose whole point is where the outcome changes, the region where
 * nothing happens should be quiet and the boundary should be the thing with ink in it. The
 * intermediate glyphs exist for exactly one reason -- a boundary that is a <b>band</b> rather
 * than a line is the entire content of the image, and a two-colour map would throw it away.
 */
public final class Heatmap {

    /**
     * Blank through to solid. Chosen for monotone visual weight in a monospaced terminal
     * font rather than for any information-theoretic property: the reader has to be able to
     * see which way is "more" without consulting the legend.
     */
    private static final char[] RAMP = {' ', '.', ':', '-', '=', '+', '*', '#'};

    private final double[] x;
    private final double[] y;
    private final double[][] values;

    private String title = "";
    private String xLabel = "";
    private String yLabel = "";
    private String lowLabel = "0";
    private String highLabel = "1";
    private int cellWidth = 4;

    private Heatmap(double[] x, double[] y, double[][] values) {
        if (values.length != y.length) {
            throw new IllegalArgumentException("values has " + values.length
                    + " rows for " + y.length + " y values");
        }
        for (double[] row : values) {
            if (row.length != x.length) {
                throw new IllegalArgumentException("a row has " + row.length
                        + " cells for " + x.length + " x values");
            }
        }
        this.x = x.clone();
        this.y = y.clone();
        this.values = new double[values.length][];
        for (int i = 0; i < values.length; i++) this.values[i] = values[i].clone();
    }

    /**
     * A map over the given axes. {@code values[row][column]} is indexed by {@code y} then
     * {@code x}, and is drawn with the first {@code y} at the bottom.
     */
    public static Heatmap of(double[] x, double[] y, double[][] values) {
        return new Heatmap(x, y, values);
    }

    public Heatmap title(String text) {
        this.title = text;
        return this;
    }

    public Heatmap xLabel(String text) {
        this.xLabel = text;
        return this;
    }

    public Heatmap yLabel(String text) {
        this.yLabel = text;
        return this;
    }

    /** What the two ends of the ramp mean, for the legend. */
    public Heatmap outcomes(String low, String high) {
        this.lowLabel = low;
        this.highLabel = high;
        return this;
    }

    /** Characters per cell. Wider cells make a map of few columns readable. */
    public Heatmap cellWidth(int width) {
        if (width < 1) throw new IllegalArgumentException("a cell is at least one character");
        this.cellWidth = width;
        return this;
    }

    public String render() {
        final int gutter = 11;
        final StringBuilder out = new StringBuilder();
        if (!title.isEmpty()) out.append(" ".repeat(gutter)).append(title).append('\n');
        if (!yLabel.isEmpty()) out.append(" ".repeat(gutter)).append(yLabel).append('\n');

        for (int r = values.length - 1; r >= 0; r--) {
            // Label the ends and roughly every fourth row between them, which keeps a tall
            // pressure axis readable without the labels running together.
            final boolean labelled =
                    r == 0 || r == values.length - 1 || r % 4 == 0;
            out.append(labelled
                    ? String.format(Locale.ROOT, "%9s |", format(y[r]))
                    : " ".repeat(9) + " |");
            for (int c = 0; c < x.length; c++) {
                out.append(String.valueOf(glyph(values[r][c])).repeat(cellWidth));
            }
            out.append('\n');
        }

        out.append(" ".repeat(9)).append(" +")
                .append("-".repeat(x.length * cellWidth)).append('\n');
        out.append(xAxis(gutter)).append('\n');
        if (!xLabel.isEmpty()) {
            final int span = x.length * cellWidth;
            out.append(" ".repeat(gutter + Math.max(0, (span - xLabel.length()) / 2)))
                    .append(xLabel).append('\n');
        }
        out.append(legend(gutter));
        return out.toString();
    }

    /**
     * One tick per cell where they fit, otherwise every other cell or every fourth, so that
     * the labels never overlap and always sit under the cell they describe.
     */
    private String xAxis(int gutter) {
        final int span = x.length * cellWidth;
        final char[] axis = new char[gutter + span];
        java.util.Arrays.fill(axis, ' ');

        int every = 1;
        while (every < x.length) {
            final int widest = widestLabel(every);
            if (widest + 1 <= cellWidth * every) break;
            every *= 2;
        }

        for (int c = 0; c < x.length; c += every) {
            final String text = format(x[c]);
            final int centre = gutter + c * cellWidth + cellWidth / 2;
            int at = centre - text.length() / 2;
            for (int k = 0; k < text.length(); k++) {
                final int position = at + k;
                if (position >= 0 && position < axis.length) axis[position] = text.charAt(k);
            }
        }
        return new String(axis);
    }

    private int widestLabel(int every) {
        int widest = 0;
        for (int c = 0; c < x.length; c += every) {
            widest = Math.max(widest, format(x[c]).length());
        }
        return widest;
    }

    private String legend(int gutter) {
        final StringBuilder out = new StringBuilder(" ".repeat(gutter));
        out.append(lowLabel).append("  ");
        for (char c : RAMP) out.append(c).append(c);
        out.append("  ").append(highLabel);
        return out.toString();
    }

    /** The ramp step a fraction falls in. Out-of-range values clamp rather than throw. */
    private static char glyph(double fraction) {
        if (Double.isNaN(fraction)) return '?';
        final int step = (int) Math.round(fraction * (RAMP.length - 1));
        return RAMP[Math.max(0, Math.min(RAMP.length - 1, step))];
    }

    /** Three significant figures, which is as much as an axis tick can carry. */
    private static String format(double value) {
        final double magnitude = Math.abs(value);
        if (magnitude >= 1000.0 || (magnitude > 0.0 && magnitude < 0.01)) {
            return String.format(Locale.ROOT, "%.1e", value);
        }
        if (magnitude >= 100.0) return String.format(Locale.ROOT, "%.0f", value);
        if (magnitude >= 10.0) return String.format(Locale.ROOT, "%.1f", value);
        return String.format(Locale.ROOT, "%.2f", value);
    }
}
