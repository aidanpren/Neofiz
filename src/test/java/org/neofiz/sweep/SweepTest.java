package org.neofiz.sweep;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicInteger;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The ensemble runner.
 *
 * <p>One property carries all the weight: results come back in parameter order whatever order
 * the runs finished in. A sweep's output is fitted, plotted, and compared against a previous
 * sweep, and every one of those becomes irreproducible if the order depends on thread
 * scheduling. So the ordering test deliberately makes the completion order the reverse of the
 * parameter order, rather than hoping a race shows up.
 */
class SweepTest {

    @Test
    @DisplayName("results are in parameter order, not completion order")
    void orderedByParameter() {
        // Later parameters finish first, by a margin far larger than scheduling noise. If the
        // gather were racing, this is the arrangement that would expose it every time rather
        // than one run in a hundred.
        final ConcurrentLinkedQueue<Integer> finished = new ConcurrentLinkedQueue<>();
        final List<Integer> parameters = List.of(0, 1, 2, 3, 4, 5, 6, 7);

        final List<Integer> results = Sweep.over(parameters, i -> {
            try {
                Thread.sleep(10L * (parameters.size() - i));
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
            }
            finished.add(i);
            return i * i;
        }, 8);

        assertEquals(List.of(0, 1, 4, 9, 16, 25, 36, 49), results);
        assertNotEquals(parameters, List.copyOf(finished),
                "the runs finished in parameter order, so this proved nothing; "
                        + "the test's own timing is broken");
    }

    @Test
    @DisplayName("every parameter is run exactly once")
    void runsEachPointOnce() {
        final AtomicInteger calls = new AtomicInteger();
        final double[] values = Sweep.linear(0.0, 1.0, 16);
        final List<Double> results = Sweep.over(values, v -> {
            calls.incrementAndGet();
            return v * 2.0;
        });
        assertEquals(16, calls.get());
        assertEquals(16, results.size());
        for (int i = 0; i < 16; i++) assertEquals(2.0 * values[i], results.get(i), 0.0);
    }

    @Test
    @DisplayName("one thread is a plain loop, and gives the same answers")
    void singleThreadedPathAgrees() {
        final List<Integer> parameters = List.of(3, 1, 4, 1, 5, 9, 2, 6);
        assertEquals(Sweep.over(parameters, i -> i * 7, 8),
                Sweep.over(parameters, i -> i * 7, 1));
    }

    @Test
    @DisplayName("a sweep with nothing in it is not an error")
    void emptySweep() {
        assertEquals(List.of(), Sweep.over(List.<Integer>of(), i -> i));
    }

    @Test
    @DisplayName("a failed point arrives as what it threw")
    void exceptionsArriveUnwrapped() {
        // A case that fails is a real result -- a tube that did not burst, a material the
        // harness refused -- and the caller has to be able to catch it as itself. Wrapping it
        // in an ExecutionException turns a legible refusal into plumbing.
        final IllegalArgumentException e = assertThrows(IllegalArgumentException.class,
                () -> Sweep.over(List.of(1, 2, 3), i -> {
                    if (i == 2) throw new IllegalArgumentException("no instability to find");
                    return i;
                }));
        assertEquals("no instability to find", e.getMessage());
    }

    @Test
    @DisplayName("a nonsense thread count is refused")
    void refusesZeroThreads() {
        assertThrows(IllegalArgumentException.class,
                () -> Sweep.over(List.of(1), i -> i, 0));
    }

    // ------------------------------------------------------------------ spacing

    @Test
    @DisplayName("a linear sweep hits both endpoints exactly")
    void linearEndpoints() {
        // Accumulating an increment would leave the last point a digit off the number that
        // was asked for, and two sweeps over the same range would then disagree at their ends.
        final double[] values = Sweep.linear(0.5, 5.0, 10);
        assertEquals(0.5, values[0], 0.0);
        assertEquals(5.0, values[9], 0.0);
        for (int i = 1; i < values.length; i++) {
            assertTrue(values[i] > values[i - 1], "not increasing at " + i);
        }
    }

    @Test
    @DisplayName("a logarithmic sweep is evenly spaced in the logarithm, endpoints exact")
    void logarithmicSpacing() {
        final double[] values = Sweep.logarithmic(10.0, 100.0, 11);
        assertEquals(10.0, values[0], 0.0);
        assertEquals(100.0, values[10], 0.0);

        final double ratio = values[1] / values[0];
        for (int i = 1; i < values.length; i++) {
            assertEquals(ratio, values[i] / values[i - 1], 1e-12,
                    "the ratio between neighbours moved at " + i);
        }
    }

    @Test
    @DisplayName("a descending range sweeps downward")
    void descendingRange() {
        // The thickness sweep is written as slenderness from high to low, because that is
        // thickness from low to high.
        final double[] values = Sweep.logarithmic(100.0, 10.0, 6);
        assertEquals(100.0, values[0], 0.0);
        assertEquals(10.0, values[5], 0.0);
        for (int i = 1; i < values.length; i++) {
            assertTrue(values[i] < values[i - 1], "not decreasing at " + i);
        }
    }

    @Test
    @DisplayName("spacings that cannot mean anything are refused")
    void refusesBadRanges() {
        assertThrows(IllegalArgumentException.class, () -> Sweep.linear(0.0, 1.0, 1));
        assertThrows(IllegalArgumentException.class, () -> Sweep.logarithmic(0.0, 1.0, 5));
        assertThrows(IllegalArgumentException.class, () -> Sweep.logarithmic(1.0, -1.0, 5));
    }

    @Test
    @DisplayName("the default thread count is fixed for the life of the process")
    void threadCountIsStable() {
        // Read once, so that two sweeps in one report are comparable as timings. This is a
        // guard on a decision, not on arithmetic.
        final List<Integer> observed = new ArrayList<>();
        for (int i = 0; i < 3; i++) observed.add(Sweep.DEFAULT_THREADS);
        assertEquals(List.of(observed.get(0), observed.get(0), observed.get(0)), observed);
        assertTrue(Sweep.DEFAULT_THREADS >= 1);
    }
}
