package org.neofiz.report;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Locale;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The trace file.
 *
 * <p>Two properties matter and neither is obvious enough to leave untested. A trace has to
 * survive the round trip <em>exactly</em>, because several of this project's answers are
 * differences in the fourth decimal place of a percentage and a file written at six
 * significant figures cannot be used to recompute one. And it has to be written the same way
 * on every machine, because a decimal comma inside a comma-separated file is a corrupted file
 * that still parses.
 */
class CsvTest {

    @Test
    @DisplayName("a double survives the round trip bit for bit")
    void exactRoundTrip(@TempDir Path dir) throws IOException {
        // Values chosen to break anything short of full precision: a repeating binary
        // fraction, a burst pressure to the pascal, a strain increment near the epsilon of
        // the numbers it is added to, and the extremes of the format.
        final double[] awkward = {
                0.1, 1.0 / 3.0, 18067432.198765431, 1e-300, 1e300,
                Math.PI, Math.nextUp(1.0), -0.0, Double.MIN_VALUE, Double.MAX_VALUE,
        };

        final Csv csv = Csv.of("v");
        for (double v : awkward) csv.row(v);
        final Path file = csv.writeTo(dir.resolve("round-trip.csv"));

        final List<String> lines = Files.readAllLines(file);
        assertEquals(awkward.length + 1, lines.size(), "one header and one row per value");
        for (int i = 0; i < awkward.length; i++) {
            final double read = Double.parseDouble(lines.get(i + 1));
            assertEquals(Double.doubleToLongBits(awkward[i]), Double.doubleToLongBits(read),
                    "value " + i + " came back as " + read + " instead of " + awkward[i]);
        }
    }

    @Test
    @DisplayName("a machine with a decimal comma still writes a readable file")
    void localeIndependent(@TempDir Path dir) throws IOException {
        // The failure this prevents is silent: the file is written, it has the right number
        // of commas plus a few, and every reader of it gets different numbers.
        final Locale original = Locale.getDefault();
        try {
            Locale.setDefault(Locale.GERMANY);
            final Path file = Csv.of("a", "b")
                    .row(1.5, 2.25)
                    .writeTo(dir.resolve("locale.csv"));
            final String row = Files.readAllLines(file).get(1);
            assertEquals(1, row.chars().filter(c -> c == ',').count(),
                    "a decimal comma got into the separators: " + row);
            assertEquals(1.5, Double.parseDouble(row.split(",")[0]), 0.0);
            assertEquals(2.25, Double.parseDouble(row.split(",")[1]), 0.0);
        } finally {
            Locale.setDefault(original);
        }
    }

    @Test
    @DisplayName("the header is the header")
    void headerIsWritten(@TempDir Path dir) throws IOException {
        final Path file = Csv.of("time_s", "capacity_Pa")
                .row(0.0, 1.0)
                .writeTo(dir.resolve("header.csv"));
        assertEquals("time_s,capacity_Pa", Files.readAllLines(file).get(0));
    }

    @Test
    @DisplayName("non-finite values are named rather than spelled in Java")
    void nonFinite(@TempDir Path dir) throws IOException {
        // A kinetic-over-strain ratio on an unloaded first sample used to be one of these.
        // Writing "NaN" is not wrong so much as unhelpful to everything downstream.
        final Path file = Csv.of("v")
                .row(Double.NaN).row(Double.POSITIVE_INFINITY).row(Double.NEGATIVE_INFINITY)
                .writeTo(dir.resolve("edges.csv"));
        assertEquals(List.of("v", "nan", "inf", "-inf"), Files.readAllLines(file));
    }

    @Test
    @DisplayName("a row that does not fit the columns is refused")
    void refusesRaggedRows() {
        final Csv csv = Csv.of("a", "b");
        final IllegalArgumentException e =
                assertThrows(IllegalArgumentException.class, () -> csv.row(1.0));
        assertTrue(e.getMessage().contains("1 values for 2 columns"), e.getMessage());
        assertThrows(IllegalArgumentException.class, () -> csv.row(1.0, 2.0, 3.0));
        assertThrows(IllegalArgumentException.class, () -> Csv.of());
    }

    @Test
    @DisplayName("the parent directory is created rather than assumed")
    void createsParent(@TempDir Path dir) {
        final Path file = Csv.of("v").row(1.0).writeTo(dir.resolve("a/b/c/deep.csv"));
        assertTrue(Files.exists(file), "the file was not created at " + file);
    }

    @Test
    @DisplayName("the run directory is not the IDE's compiler output")
    void runsDirectoryIsItsOwn(@TempDir Path dir) {
        // out/ is claimed by this project's IDE configuration as a compiler output root, and
        // a directory something else empties is not a place to keep run traces. This is the
        // guard on a decision that is otherwise only recorded in a comment.
        final String original = System.getProperty(Runs.PROPERTY);
        try {
            System.clearProperty(Runs.PROPERTY);
            final String name = Runs.dir().getFileName().toString();
            assertEquals("runs", name, "the default run directory moved to " + name);

            System.setProperty(Runs.PROPERTY, dir.resolve("elsewhere").toString());
            final Path redirected = Runs.file("trace.csv");
            assertEquals(dir.resolve("elsewhere").resolve("trace.csv"), redirected);
            assertTrue(Files.isDirectory(dir.resolve("elsewhere")),
                    "the redirected directory was not created");
        } finally {
            if (original == null) System.clearProperty(Runs.PROPERTY);
            else System.setProperty(Runs.PROPERTY, original);
        }
    }

    @Test
    @DisplayName("a labelled row puts a name in front of the numbers")
    void labelledRows() throws Exception {
        // The one text case that has come up: a table of measurements, one row per named
        // configuration. Anything more general needs quoting rules for every column and a
        // decision about missing values, neither of which has a caller.
        Path file = Files.createTempFile("neofiz-label", ".csv");
        Csv.of("configuration", "rate")
                .row("bare elements", 9.0)
                .row("with, commas", 3.5)
                .writeTo(file);

        List<String> lines = Files.readAllLines(file);
        assertEquals("configuration,rate", lines.get(0));
        assertTrue(lines.get(1).startsWith("bare elements,9.0"), lines.get(1));
        assertTrue(lines.get(2).startsWith("\"with, commas\",3.5"),
                "a comma in a label has to be quoted: " + lines.get(2));
        Files.delete(file);
    }

    @Test
    @DisplayName("a labelled row still has to fill the table")
    void labelledRowWidthIsChecked() {
        Csv csv = Csv.of("name", "a", "b");
        assertThrows(IllegalArgumentException.class, () -> csv.row("x", 1.0));
        assertThrows(IllegalArgumentException.class, () -> csv.row("x", 1.0, 2.0, 3.0));
    }
}
