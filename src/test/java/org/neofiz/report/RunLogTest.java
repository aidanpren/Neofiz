package org.neofiz.report;

import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The log of what has been run.
 *
 * <p>Every one of these writes to a temporary run directory rather than the real one, because
 * a test that appends to the project's own log would make the log a record of the test suite.
 * That is what the {@code neofiz.runs} property is for and it is worth exercising here: if
 * redirecting the run directory did not work, a test would say so by polluting the one beside
 * it.
 */
class RunLogTest {

    private Path directory;
    private String previous;

    @BeforeEach
    void redirectRunsDirectory() throws IOException {
        directory = Files.createTempDirectory("neofiz-log");
        previous = System.getProperty(Runs.PROPERTY);
        System.setProperty(Runs.PROPERTY, directory.toString());
    }

    @AfterEach
    void restore() throws IOException {
        if (previous == null) System.clearProperty(Runs.PROPERTY);
        else System.setProperty(Runs.PROPERTY, previous);
        try (var walk = Files.walk(directory)) {
            for (Path p : walk.sorted(Comparator.reverseOrder()).toList()) Files.deleteIfExists(p);
        }
    }

    @Test
    @DisplayName("a run comes back with its parameters and its results")
    void roundTrip() {
        RunLog.of("smash").parameter("speed", 700.0).parameter("cell", 0.001)
                .result("eroded", 32).result("plastic", 1618.9).commit();

        List<RunLog.Entry> all = RunLog.all();
        assertEquals(1, all.size());
        RunLog.Entry e = all.get(0);
        assertEquals("smash", e.program());
        assertEquals("700", e.parameters().get("speed"));
        assertEquals("0.00100000", e.parameters().get("cell"));
        assertEquals(32.0, e.number("eroded"));
        assertEquals(1618.90, e.number("plastic"), 0.01);
    }

    @Test
    @DisplayName("the same run twice is recognised, a different one is not")
    void matchingIsOnTheWholeParameterSet() {
        // Exactly, in both directions. A looser rule would answer "yes, this has been run"
        // for a run that was not the same run, which is worse than answering nothing.
        RunLog.of("smash").parameter("speed", 700.0).result("eroded", 32).commit();

        assertTrue(RunLog.of("smash").parameter("speed", 700.0).lastMatching().isPresent());
        assertFalse(RunLog.of("smash").parameter("speed", 420.0).lastMatching().isPresent());
        assertFalse(RunLog.of("topple").parameter("speed", 700.0).lastMatching().isPresent());
        assertFalse(RunLog.of("smash").parameter("speed", 700.0)
                .parameter("cell", 0.001).lastMatching().isPresent(),
                "an extra parameter makes it a different run");
    }

    @Test
    @DisplayName("the newest matching run wins")
    void matchingReturnsTheLatest() {
        RunLog.of("smash").parameter("speed", 700.0).result("eroded", 32).commit();
        RunLog.of("smash").parameter("speed", 700.0).result("eroded", 41).commit();

        Optional<RunLog.Entry> found = RunLog.of("smash").parameter("speed", 700.0).lastMatching();
        assertTrue(found.isPresent());
        assertEquals(41.0, found.get().number("eroded"));
    }

    @Test
    @DisplayName("runs of one program can be picked out of the log")
    void historyIsPerProgram() {
        RunLog.of("smash").parameter("speed", 700.0).commit();
        RunLog.of("topple").parameter("blocks", 5).commit();
        RunLog.of("smash").parameter("speed", 420.0).commit();

        assertEquals(3, RunLog.all().size());
        assertEquals(2, RunLog.historyOf("smash").size());
        assertEquals(1, RunLog.historyOf("topple").size());
        assertEquals(0, RunLog.historyOf("nothing").size());
    }

    @Test
    @DisplayName("the file is a readable four-column CSV with one header")
    void theFileIsReadable() throws IOException {
        // The whole reason for a text format rather than a database. Two runs, one header,
        // and a person can open it.
        RunLog.of("smash").parameter("speed", 700.0).result("eroded", 32).commit();
        RunLog.of("smash").parameter("speed", 420.0).result("eroded", 6).commit();

        List<String> lines = Files.readAllLines(directory.resolve(RunLog.FILE));
        assertEquals(3, lines.size());
        assertEquals("when,program,parameters,results", lines.get(0));
        assertTrue(lines.get(1).endsWith(",smash,speed=700,eroded=32"), lines.get(1));
        assertTrue(lines.get(2).endsWith(",smash,speed=420,eroded=6"), lines.get(2));
    }

    @Test
    @DisplayName("a name that would break the format is refused rather than escaped")
    void separatorsAreRefused() {
        // Escaping needs the convention to be right in two places, to support names nothing
        // here wants. Refusing is one place and it cannot be subtly wrong.
        RunLog log = RunLog.of("smash");
        assertThrows(IllegalArgumentException.class, () -> log.parameter("a,b", "1"));
        assertThrows(IllegalArgumentException.class, () -> log.parameter("a;b", "1"));
        assertThrows(IllegalArgumentException.class, () -> log.parameter("a=b", "1"));
        assertThrows(IllegalArgumentException.class, () -> log.parameter("ok", "a;b"));
        assertThrows(IllegalArgumentException.class, () -> RunLog.of("with,comma"));
        assertThrows(IllegalArgumentException.class, () -> log.parameter("", "1"));
    }

    @Test
    @DisplayName("an empty log is empty rather than absent")
    void emptyLogReadsAsEmpty() {
        assertEquals(List.of(), RunLog.all());
        assertFalse(RunLog.of("smash").lastMatching().isPresent());
    }
}
