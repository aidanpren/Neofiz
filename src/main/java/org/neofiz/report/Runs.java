package org.neofiz.report;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

/**
 * Where run output lands.
 *
 * <p>The build plan is firm that runs accumulate -- every run stays in the log with its
 * parameters, traces and outcome, because the accumulated record is the player's experience.
 * This is the floor of that: a directory, and a rule about which one. The log's structure,
 * its index and its parameter provenance belong to the experiment layer and are not invented
 * here.
 *
 * <p>The directory is {@code runs/} beside the build file, overridable with the
 * {@code neofiz.runs} system property. It is deliberately <b>not</b> {@code out/}, which this
 * project's IDE configuration already claims as a compiler output root -- writing traces into
 * a directory something else empties is a way to lose data quietly.
 */
public final class Runs {

    /** System property naming an alternative output directory. */
    public static final String PROPERTY = "neofiz.runs";

    private static final String DEFAULT = "runs";

    private Runs() {
    }

    /** The run-output directory, created if it does not exist. */
    public static Path dir() {
        final Path dir = Paths.get(System.getProperty(PROPERTY, DEFAULT));
        try {
            Files.createDirectories(dir);
        } catch (IOException e) {
            throw new UncheckedIOException("cannot create the run directory " + dir, e);
        }
        return dir;
    }

    /** A file in the run-output directory. */
    public static Path file(String name) {
        return dir().resolve(name);
    }
}
