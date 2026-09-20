package org.neofiz.report;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The terminal plot.
 *
 * <p>Nothing here is about whether a chart looks nice, which is not testable and not the
 * point. What is testable is that the mapping from numbers to characters is faithful: that a
 * maximum in the data is a maximum on the screen, that a curve is connected rather than a
 * scatter of dots, and that the degenerate inputs a solver will eventually hand it -- a flat
 * series, a single sample, a value outside a fixed range -- render instead of throwing.
 * A plotting routine that falls over on the first unusual run is worse than no plot, because
 * it takes the run down with it.
 */
class PlotTest {

    /** The plotting area of a plot of this size, without the gutter or the axis rows. */
    private static String[] interior(String rendered, int width) {
        return rendered.lines()
                .filter(line -> line.contains(" |"))
                .map(line -> line.substring(line.indexOf(" |") + 2))
                .map(line -> line.length() >= width ? line.substring(0, width) : line)
                .toArray(String[]::new);
    }

    @Test
    @DisplayName("a maximum in the data is a maximum on the screen")
    void peakIsWhereThePeakIs() {
        // A parabola with its top at x = 5, which is the middle of the range, so the highest
        // drawn row has to be nearer the middle column than either edge. This is the whole
        // contract: if it fails, every curve this thing draws is lying about where its
        // interesting point is.
        final int n = 101;
        final double[] x = new double[n];
        final double[] y = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = 10.0 * i / (n - 1.0);
            y[i] = -(x[i] - 5.0) * (x[i] - 5.0);
        }

        final int width = 60;
        final String[] rows = interior(
                Plot.of(width, 15).series("p", '#', x, y).render(), width);

        int topmost = -1;
        for (int r = 0; r < rows.length && topmost < 0; r++) {
            if (rows[r].indexOf('#') >= 0) topmost = r;
        }
        assertEquals(0, topmost, "the maximum should reach the top row");

        final int first = rows[0].indexOf('#');
        final int last = rows[0].lastIndexOf('#');
        final int centre = (first + last) / 2;
        assertTrue(Math.abs(centre - (width - 1) / 2) <= 2,
                "the peak is drawn at column " + centre + " of " + width);
    }

    @Test
    @DisplayName("a curve is connected, not a scatter of dots")
    void steepSegmentsAreJoined() {
        // Two samples that land eleven rows apart. Plotting points alone would leave nine
        // blank rows between them and read as two unrelated marks; the whole visual claim of
        // a line chart is that consecutive samples are consecutive.
        final String[] rows = interior(
                Plot.of(20, 12).series("s", '#', new double[] {0, 1}, new double[] {0, 1})
                        .render(), 20);
        for (String row : rows) {
            assertTrue(row.indexOf('#') >= 0, "a row of the traverse is empty: [" + row + "]");
        }
    }

    @Test
    @DisplayName("a series that never moves renders as a flat line")
    void flatSeriesDoesNotDivideByZero() {
        // The zero-width range is the case that turns a plot into a crash, and a solver that
        // has not yielded yet produces exactly it.
        final String rendered = Plot.of(30, 8)
                .series("flat", '-', new double[] {0, 1, 2}, new double[] {7, 7, 7})
                .render();
        assertTrue(rendered.contains("-"), "nothing was drawn");
        assertTrue(rendered.contains("7"), "the axis should still say what the value is");
    }

    @Test
    @DisplayName("a single point is a plot")
    void singlePoint() {
        // Both ranges are degenerate here, which is the same padding path as the flat series
        // taken in two dimensions at once.
        final String[] rows = interior(
                Plot.of(20, 5).series("sample", 'o', new double[] {3}, new double[] {4})
                        .render(), 20);
        final long drawn = java.util.Arrays.stream(rows)
                .flatMapToInt(String::chars).filter(c -> c == 'o').count();
        assertEquals(1, drawn, "one sample should be exactly one glyph");
    }

    @Test
    @DisplayName("data outside a fixed range is clipped, not thrown")
    void fixedRangeClips() {
        // Fixed ranges exist so that a zoomed panel can show the top of a curve. The rest of
        // the curve does not stop existing, and it must not take the render down with it.
        final String rendered = Plot.of(24, 6)
                .yRange(0.0, 1.0)
                .series("s", '#', new double[] {0, 1, 2}, new double[] {-50, 0.5, 50})
                .render();
        assertTrue(rendered.contains("#"), "the in-range sample should still be drawn");
    }

    @Test
    @DisplayName("markers are drawn over the curves they annotate")
    void markersWin() {
        // The burst point sits on the capacity curve by construction, so if series glyphs
        // won, the marker would never be visible on any real plot.
        final String rendered = Plot.of(21, 5)
                .series("s", '#', new double[] {0, 1, 2}, new double[] {0, 1, 2})
                .marker("peak", '*', 1, 1)
                .render();
        final String[] rows = interior(rendered, 21);
        assertEquals('*', rows[2].charAt(10), "the marker should own the middle cell");
    }

    @Test
    @DisplayName("axis labels carry enough digits to tell the ticks apart")
    void axisResolvesNarrowRanges() {
        // The zoomed burst panel spans about 0.5 % of 18 MPa. Two decimal places would print
        // the same number on every tick and the panel would be unreadable.
        final String rendered = Plot.of(60, 10)
                .yRange(17.9, 18.05)
                .series("s", '#', new double[] {0, 1}, new double[] {17.95, 18.0})
                .render();
        assertTrue(rendered.contains("18.05") || rendered.contains("18.050"),
                "the top of the y axis is not labelled distinguishably:\n" + rendered);
        assertTrue(rendered.contains("17.90") || rendered.contains("17.900"),
                "the bottom of the y axis is not labelled distinguishably:\n" + rendered);
    }

    @Test
    @DisplayName("the legend names every series and marker")
    void legendIsComplete() {
        final String rendered = Plot.standard()
                .series("capacity", '#', new double[] {0, 1}, new double[] {0, 1})
                .series("applied", '.', new double[] {0, 1}, new double[] {1, 1})
                .marker("burst", '*', 1, 1)
                .render();
        assertTrue(rendered.contains("# capacity"), rendered);
        assertTrue(rendered.contains(". applied"), rendered);
        assertTrue(rendered.contains("* burst"), rendered);
    }

    @Test
    @DisplayName("mismatched or empty inputs are refused rather than half-drawn")
    void refusesNonsense() {
        assertThrows(IllegalArgumentException.class,
                () -> Plot.standard().series("s", '#', new double[] {0, 1}, new double[] {0}));
        assertThrows(IllegalArgumentException.class, () -> Plot.of(4, 4));
        assertThrows(IllegalStateException.class, () -> Plot.standard().render());
    }
}
