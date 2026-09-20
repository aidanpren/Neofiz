package org.neofiz.report;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The two-parameter map.
 *
 * <p>What has to be right is the geometry: the first row of data at the bottom, the first
 * column at the left, and the ramp running the way the legend says it does. A map that is
 * upside down still looks like a map, which is exactly why it needs asserting rather than
 * eyeballing -- the design map's whole job is to show which corner is safe.
 */
class HeatmapTest {

    /** The cell rows of a rendered map, with the gutter stripped. */
    private static String[] cells(String rendered) {
        return rendered.lines()
                .filter(line -> line.contains(" |"))
                .map(line -> line.substring(line.indexOf(" |") + 2))
                .toArray(String[]::new);
    }

    @Test
    @DisplayName("the first row of data is drawn at the bottom")
    void originIsBottomLeft() {
        // values[0] is the lowest y. If this inverted, a design map would put the safe region
        // at the top and read as the exact opposite of what it means.
        final String[] rows = cells(Heatmap.of(
                new double[] {1, 2}, new double[] {10, 20},
                new double[][] {{0.0, 0.0}, {1.0, 1.0}}).cellWidth(1).render());

        assertEquals(2, rows.length);
        assertEquals("##", rows[0], "the high-y row should be drawn first, at the top");
        assertEquals("  ", rows[1], "the low-y row should be drawn last, at the bottom");
    }

    @Test
    @DisplayName("the first column of data is drawn at the left")
    void columnsAreInOrder() {
        final String row = cells(Heatmap.of(
                new double[] {1, 2, 3}, new double[] {10},
                new double[][] {{0.0, 0.5, 1.0}}).cellWidth(1).render())[0];

        assertEquals(' ', row.charAt(0), "the zero column should be leftmost: [" + row + "]");
        assertEquals('#', row.charAt(2), "the one column should be rightmost: [" + row + "]");
        assertTrue(row.charAt(1) != ' ' && row.charAt(1) != '#',
                "a half should be neither end of the ramp: [" + row + "]");
    }

    @Test
    @DisplayName("the ramp is monotone from blank to solid")
    void rampIsMonotone() {
        // Blank for zero is the design decision under test: on a map whose point is where the
        // outcome changes, the region where nothing happens should carry no ink.
        final double[] values = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
        final double[][] row = {values};
        final String drawn = cells(Heatmap.of(
                new double[] {1, 2, 3, 4, 5, 6}, new double[] {1}, row)
                .cellWidth(1).render())[0].substring(0, 6);

        assertEquals(' ', drawn.charAt(0), "zero should be blank");
        assertEquals('#', drawn.charAt(5), "one should be solid");
        for (int i = 1; i < 6; i++) {
            assertTrue(drawn.charAt(i) != drawn.charAt(i - 1) || i == 1,
                    "the ramp repeated a glyph at " + i + " in [" + drawn + "]");
        }
    }

    @Test
    @DisplayName("cells are as wide as they were asked to be")
    void cellWidthIsHonoured() {
        final String[] rows = cells(Heatmap.of(
                new double[] {1, 2}, new double[] {1},
                new double[][] {{1.0, 0.0}}).cellWidth(5).render());
        assertTrue(rows[0].startsWith("#####     "),
                "expected five solid then five blank; got [" + rows[0] + "]");
    }

    @Test
    @DisplayName("values outside [0, 1] clamp, and a hole is marked rather than drawn as safe")
    void outOfRangeAndHoles() {
        // A NaN cell is a cell with no data in it. Rendering it as blank would put it in the
        // "everything held" region, which is the most dangerous possible default.
        final String[] rows = cells(Heatmap.of(
                new double[] {1, 2, 3}, new double[] {1},
                new double[][] {{-5.0, Double.NaN, 5.0}}).cellWidth(1).render());
        assertEquals(' ', rows[0].charAt(0));
        assertEquals('?', rows[0].charAt(1));
        assertEquals('#', rows[0].charAt(2));
    }

    @Test
    @DisplayName("the legend says which end is which")
    void legendNamesTheEnds() {
        final String rendered = Heatmap.of(new double[] {1}, new double[] {1},
                new double[][] {{0.5}}).outcomes("all hold", "all burst").render();
        assertTrue(rendered.contains("all hold"), rendered);
        assertTrue(rendered.contains("all burst"), rendered);
    }

    @Test
    @DisplayName("a grid whose shape does not match its axes is refused")
    void refusesMismatchedShape() {
        assertThrows(IllegalArgumentException.class, () -> Heatmap.of(
                new double[] {1, 2}, new double[] {1, 2}, new double[][] {{0, 0}}));
        assertThrows(IllegalArgumentException.class, () -> Heatmap.of(
                new double[] {1, 2}, new double[] {1}, new double[][] {{0, 0, 0}}));
        assertThrows(IllegalArgumentException.class, () -> Heatmap.of(
                new double[] {1}, new double[] {1}, new double[][] {{0}}).cellWidth(0));
    }

    @Test
    @DisplayName("the map does not alias its caller's arrays")
    void defensivelyCopies() {
        final double[][] values = {{0.0}};
        final Heatmap map = Heatmap.of(new double[] {1}, new double[] {1}, values);
        values[0][0] = 1.0;
        assertEquals(' ', cells(map.cellWidth(1).render())[0].charAt(0),
                "the map re-read its caller's array after being built");
    }
}
