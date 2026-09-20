package org.neofiz.report;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * Numeric tables out, for anything that wants to plot them properly.
 *
 * <p>{@link Plot} keeps a curve visible in the terminal; this is the same curve in a form a
 * spreadsheet, a notebook or the eventual instrument panel can read. Together they are the
 * whole of the project's output story until the experiment layer exists, and the split is
 * intentional: the terminal plot is for the person running the solver right now, the file is
 * for everything that comes later.
 *
 * <h2>Formatting</h2>
 *
 * Seventeen significant digits, which is {@code %.17g} -- enough that a {@code double}
 * survives the round trip exactly. This is not fussiness. Several of the numbers in this
 * project are differences in the fourth decimal place of a percentage, and a trace written at
 * six digits cannot be used to recompute one.
 *
 * <p>{@link Locale#ROOT} throughout, so a machine with a comma decimal separator does not
 * write a comma-separated file whose numbers contain commas.
 */
public final class Csv {

    private final List<String> headers = new ArrayList<>();
    private final List<double[]> rows = new ArrayList<>();
    /** The first column of each row when it is text, or null when every column is a number. */
    private final List<String> labels = new ArrayList<>();
    private boolean labelled;

    private Csv(String... headers) {
        if (headers.length == 0) throw new IllegalArgumentException("a table needs columns");
        this.headers.addAll(List.of(headers));
    }

    /** A table with the given column headings. Include units in them. */
    public static Csv of(String... headers) {
        return new Csv(headers);
    }

    /** One row. Must have exactly as many values as there are columns. */
    public Csv row(double... values) {
        if (values.length != headers.size()) {
            throw new IllegalArgumentException("row has " + values.length
                    + " values for " + headers.size() + " columns");
        }
        rows.add(values.clone());
        labels.add(null);
        return this;
    }

    /**
     * One row whose first column is a name rather than a number.
     *
     * <p>A separate method rather than an overload taking objects, because a table where any
     * column might be text needs quoting rules, an escaping convention and a decision about
     * what a missing value looks like -- and the only case that has come up is a label in
     * front of a row of measurements. Commas and quotes in the label are escaped; nothing
     * else about the format changes.
     */
    public Csv row(String label, double... values) {
        if (values.length != headers.size() - 1) {
            throw new IllegalArgumentException("row has " + (values.length + 1)
                    + " values for " + headers.size() + " columns");
        }
        rows.add(values.clone());
        labels.add(label);
        labelled = true;
        return this;
    }

    public int size() {
        return rows.size();
    }

    /** Writes to a named file in the {@link Runs} directory and returns its path. */
    public Path write(String name) {
        return writeTo(Runs.file(name));
    }

    /** Writes to an explicit path, creating the parent directory, and returns the path. */
    public Path writeTo(Path file) {
        try {
            final Path parent = file.getParent();
            if (parent != null) Files.createDirectories(parent);
            try (BufferedWriter out = Files.newBufferedWriter(file, StandardCharsets.UTF_8)) {
                out.write(String.join(",", headers));
                out.newLine();
                for (int k = 0; k < rows.size(); k++) {
                    final double[] row = rows.get(k);
                    final String label = labels.get(k);
                    if (labelled) out.write(quote(label == null ? "" : label));
                    for (int i = 0; i < row.length; i++) {
                        if (i > 0 || labelled) out.write(',');
                        out.write(format(row[i]));
                    }
                    out.newLine();
                }
            }
        } catch (IOException e) {
            throw new UncheckedIOException("cannot write " + file, e);
        }
        return file;
    }

    /** A label, quoted only when it has to be. */
    private static String quote(String text) {
        final boolean plain = text.indexOf(',') < 0
                && text.indexOf('"') < 0
                && text.indexOf('\n') < 0
                && text.indexOf('\r') < 0;
        return plain ? text : '"' + text.replace("\"", "\"\"") + '"';
    }

    /**
     * A double at full precision, with the non-finite values written as the bare words a
     * reader will recognize rather than as Java's {@code NaN}/{@code Infinity} spelling.
     */
    private static String format(double value) {
        if (Double.isNaN(value)) return "nan";
        if (Double.isInfinite(value)) return value > 0.0 ? "inf" : "-inf";
        return String.format(Locale.ROOT, "%.17g", value);
    }
}
