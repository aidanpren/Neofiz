package org.neofiz.report;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.time.Instant;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Optional;

/**
 * What was run, with what, and what came out.
 *
 * <p>Every program in this project prints a few lines and writes a file, and then the terminal
 * scrolls. That is fine while there is one run; it stops being fine at the point where the
 * interesting question becomes "is this number different from last time, and if so which of
 * the eleven things I changed did it". The answer needs the parameters to travel with the
 * results, and it needs the last run to still be there.
 *
 * <h2>Why not a database, and why not JSON</h2>
 *
 * A line per run, appended, in a file a person can read. Appending is the only write, so two
 * runs going at once cannot corrupt each other beyond interleaving whole lines, and there is
 * no schema to migrate when a program starts recording one more number. The cost is that
 * querying is a linear scan of a text file, which for a log that grows by one line per run is
 * not a cost at all.
 *
 * <p>Parameters and results are {@code key=value} pairs joined by semicolons, so the file is
 * still a four-column CSV and a spreadsheet can open it. Keys and values containing either
 * character are refused rather than escaped: the alternative is an escaping convention that
 * has to be got right in two places, to support names nothing here wants.
 *
 * <h2>The parameters are the identity</h2>
 *
 * {@link #lastMatching} asks whether a run has been done before, by comparing the whole
 * parameter set rather than a name or a hash. That is deliberate. A hash would be shorter and
 * would answer the same question, and would also make the file useless to read -- and the
 * point of a log nobody can read is not clear.
 */
public final class RunLog {

    /** The file, in the {@link Runs} directory. */
    public static final String FILE = "log.csv";

    private static final String[] HEADERS = {"when", "program", "parameters", "results"};

    /** One line of the log. */
    public record Entry(Instant when, String program,
                        Map<String, String> parameters, Map<String, String> results) {

        /** A result as a number, or NaN if it is missing or not one. */
        public double number(String key) {
            final String text = results.get(key);
            if (text == null) return Double.NaN;
            try {
                return Double.parseDouble(text);
            } catch (NumberFormatException e) {
                return Double.NaN;
            }
        }
    }

    private final String program;
    private final Map<String, String> parameters = new LinkedHashMap<>();
    private final Map<String, String> results = new LinkedHashMap<>();

    private RunLog(String program) {
        this.program = check("program", program);
    }

    /** Starts recording a run of a named program. */
    public static RunLog of(String program) {
        return new RunLog(program);
    }

    /** An input. Everything that would change the answer belongs here. */
    public RunLog parameter(String key, double value) {
        return parameter(key, format(value));
    }

    public RunLog parameter(String key, String value) {
        parameters.put(check("key", key), check("value", value));
        return this;
    }

    /** An output. Everything worth comparing against the next run belongs here. */
    public RunLog result(String key, double value) {
        results.put(check("key", key), format(value));
        return this;
    }

    public RunLog result(String key, String value) {
        results.put(check("key", key), check("value", value));
        return this;
    }

    /** Appends the run to the log and returns what was written. */
    public Entry commit() {
        final Entry entry = new Entry(Instant.now(), program,
                new LinkedHashMap<>(parameters), new LinkedHashMap<>(results));
        final Path file = Runs.file(FILE);
        try {
            final boolean fresh = !Files.exists(file) || Files.size(file) == 0;
            final StringBuilder line = new StringBuilder();
            if (fresh) line.append(String.join(",", HEADERS)).append(System.lineSeparator());
            line.append(entry.when()).append(',').append(entry.program()).append(',')
                    .append(pack(entry.parameters())).append(',')
                    .append(pack(entry.results())).append(System.lineSeparator());
            Files.writeString(file, line.toString(), StandardCharsets.UTF_8,
                    StandardOpenOption.CREATE, StandardOpenOption.APPEND);
        } catch (IOException e) {
            throw new UncheckedIOException("cannot append to " + file, e);
        }
        return entry;
    }

    /** Every run recorded so far, oldest first. Empty when nothing has been logged. */
    public static List<Entry> all() {
        final Path file = Runs.file(FILE);
        if (!Files.exists(file)) return List.of();
        final List<Entry> out = new ArrayList<>();
        try {
            for (String line : Files.readAllLines(file, StandardCharsets.UTF_8)) {
                if (line.isBlank() || line.startsWith(HEADERS[0] + ",")) continue;
                final String[] parts = line.split(",", 4);
                if (parts.length < 4) continue;
                out.add(new Entry(Instant.parse(parts[0]), parts[1],
                        unpack(parts[2]), unpack(parts[3])));
            }
        } catch (IOException e) {
            throw new UncheckedIOException("cannot read " + file, e);
        }
        return out;
    }

    /** Runs of one program, oldest first. */
    public static List<Entry> historyOf(String program) {
        final List<Entry> out = new ArrayList<>();
        for (Entry e : all()) if (e.program().equals(program)) out.add(e);
        return out;
    }

    /**
     * The most recent run of this program with exactly these parameters, if there is one.
     *
     * <p>Exactly, in both directions: a run recorded with an extra parameter does not match a
     * query without it. Anything looser would answer "yes, this has been run" for a run that
     * was not the same run, which is worse than answering nothing.
     */
    public Optional<Entry> lastMatching() {
        Entry best = null;
        for (Entry e : all()) {
            if (e.program().equals(program) && e.parameters().equals(parameters)) best = e;
        }
        return Optional.ofNullable(best);
    }

    private static String pack(Map<String, String> map) {
        final StringBuilder out = new StringBuilder();
        for (Map.Entry<String, String> e : map.entrySet()) {
            if (out.length() > 0) out.append(';');
            out.append(e.getKey()).append('=').append(e.getValue());
        }
        return out.toString();
    }

    private static Map<String, String> unpack(String text) {
        final Map<String, String> out = new LinkedHashMap<>();
        if (text.isEmpty()) return out;
        for (String pair : text.split(";")) {
            final int at = pair.indexOf('=');
            if (at > 0) out.put(pair.substring(0, at), pair.substring(at + 1));
        }
        return out;
    }

    /**
     * Enough digits to tell two runs apart, which is what this file is for.
     *
     * <p>Not {@link Csv}'s seventeen. A log is read by people and compared by eye, and a
     * column of twenty-character numbers is unreadable; a trace that has to be recomputed
     * from belongs in a Csv beside it.
     */
    private static String format(double value) {
        if (Double.isNaN(value)) return "nan";
        if (Double.isInfinite(value)) return value > 0.0 ? "inf" : "-inf";
        if (value == Math.rint(value) && Math.abs(value) < 1e15) {
            return String.format(Locale.ROOT, "%.0f", value);
        }
        return String.format(Locale.ROOT, "%.6g", value);
    }

    private static String check(String what, String text) {
        if (text == null || text.isEmpty()) {
            throw new IllegalArgumentException(what + " must not be empty");
        }
        if (text.indexOf(',') >= 0 || text.indexOf(';') >= 0 || text.indexOf('=') >= 0
                || text.indexOf('\n') >= 0 || text.indexOf('\r') >= 0) {
            throw new IllegalArgumentException(
                    what + " must not contain a comma, semicolon, equals or newline; the log "
                            + "keeps them as separators rather than escaping them: " + text);
        }
        return text;
    }
}
