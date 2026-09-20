package org.neofiz.sweep;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.function.DoubleFunction;
import java.util.function.Function;

/**
 * One parameter varied, every value run, the results in order.
 *
 * <p>The build plan's whole argument for the experiment layer is that people do not extract
 * physical law from single observations -- they extract it from contrast, from varying one
 * quantity and watching the response. This is the smallest object that provides contrast:
 * take a case, take a list of values for one of its parameters, and return a result per
 * value. Everything the experiment layer eventually does -- design maps, scatter, scaling
 * laws, inverse design -- is built on the ensemble being the unit of work rather than the
 * single run.
 *
 * <h2>Where the parallelism goes, and why it goes there</h2>
 *
 * {@link org.neofiz.solver.ExplicitSolver} can split its element kernel across threads, and
 * for a sweep it should not. Two reasons, and the second is the one that bites.
 *
 * <ul>
 *   <li>A sweep is embarrassingly parallel across <em>runs</em>, which have no barrier
 *       between them at all, while threading inside a run pays a barrier twice per step. One
 *       run per thread is the better split whenever there are at least as many runs as
 *       cores, which is every sweep worth taking.</li>
 *   <li>The solver's worker threads <b>spin</b> rather than sleeping between parallel
 *       regions, which is the right trade for one offline run on an idle machine and a
 *       catastrophic one under oversubscription: eight solvers each asking for eight
 *       spinning workers puts sixty-four runnable threads on eight cores, and most of the
 *       machine goes into waiting for a barrier held by a thread that is not scheduled.</li>
 * </ul>
 *
 * <p>So the contract is that a sweep runs each case on one thread and the case must not ask
 * for more. Nothing here enforces that -- the solver is constructed inside the caller's
 * function, where this class cannot see it -- but the default is one thread, so a case has
 * to go out of its way to break it.
 *
 * <h2>Determinism</h2>
 *
 * Results come back in the order the parameters were given, whatever order the runs actually
 * finished in. That is not cosmetic: a sweep's output is fitted, plotted and compared against
 * a previous sweep, and a list whose order depends on thread scheduling makes every one of
 * those operations non-reproducible. The runs themselves are already deterministic -- the
 * solver is bit-for-bit identical at any thread count -- and they share no state, so an
 * ordered gather is the whole of what is needed.
 */
public final class Sweep {

    /**
     * Threads a sweep uses by default: one per available core.
     *
     * <p>Read once. {@link Runtime#availableProcessors()} can change between calls, and a
     * sweep whose thread count moved underneath it would be reporting timings that cannot be
     * compared with each other.
     */
    public static final int DEFAULT_THREADS = Runtime.getRuntime().availableProcessors();

    private Sweep() {
    }

    /**
     * One result per parameter, in parameter order.
     *
     * @param parameters the values to run, in the order they should be reported
     * @param run        the case; must be safe to call from several threads at once, which
     *                   for this project means it must not share a solver or a mesh
     */
    public static <P, R> List<R> over(List<P> parameters, Function<P, R> run) {
        return over(parameters, run, DEFAULT_THREADS);
    }

    /** As {@link #over(List, Function)}, with the thread count named. */
    public static <P, R> List<R> over(List<P> parameters, Function<P, R> run, int threads) {
        if (threads < 1) throw new IllegalArgumentException("need at least one thread");
        if (parameters.isEmpty()) return List.of();

        // More threads than runs is not wrong, only wasteful, and on a spin-heavy solver it
        // is worth not doing: the extra workers would sit in the pool doing nothing.
        final int width = Math.min(threads, parameters.size());
        if (width == 1) {
            final List<R> results = new ArrayList<>(parameters.size());
            for (P p : parameters) results.add(run.apply(p));
            return List.copyOf(results);
        }

        try (ExecutorService pool = Executors.newFixedThreadPool(width)) {
            final List<Future<R>> futures = new ArrayList<>(parameters.size());
            for (P p : parameters) {
                final Callable<R> task = () -> run.apply(p);
                futures.add(pool.submit(task));
            }
            final List<R> results = new ArrayList<>(parameters.size());
            for (Future<R> f : futures) results.add(join(f));
            return List.copyOf(results);
        }
    }

    /** One result per value of a numeric parameter, in order. */
    public static <R> List<R> over(double[] values, DoubleFunction<R> run) {
        return over(values, run, DEFAULT_THREADS);
    }

    /** As {@link #over(double[], DoubleFunction)}, with the thread count named. */
    public static <R> List<R> over(double[] values, DoubleFunction<R> run, int threads) {
        return over(Arrays.stream(values).boxed().toList(),
                (Double value) -> run.apply(value), threads);
    }

    /**
     * {@code count} values from {@code low} to {@code high} inclusive, evenly spaced.
     *
     * <p>The endpoints are exact rather than accumulated, so a sweep's last point is the
     * number that was asked for and two sweeps over the same range agree at their ends.
     */
    public static double[] linear(double low, double high, int count) {
        if (count < 2) throw new IllegalArgumentException("a sweep needs at least two points");
        final double[] values = new double[count];
        for (int i = 0; i < count; i++) {
            values[i] = low + (high - low) * i / (count - 1.0);
        }
        return values;
    }

    /**
     * {@code count} values from {@code low} to {@code high} inclusive, evenly spaced in the
     * logarithm.
     *
     * <p>This is the spacing a scaling law wants. A power law is a straight line in log-log,
     * and fitting one from points bunched at the top of the range weights the fit by where
     * the samples happened to fall rather than by the decades they cover.
     */
    public static double[] logarithmic(double low, double high, int count) {
        if (low <= 0.0 || high <= 0.0) {
            throw new IllegalArgumentException("a logarithmic sweep needs positive endpoints");
        }
        final double[] values = linear(Math.log(low), Math.log(high), count);
        for (int i = 0; i < values.length; i++) values[i] = Math.exp(values[i]);
        // Exact at the ends, where exp(log(x)) is otherwise a digit out.
        values[0] = low;
        values[values.length - 1] = high;
        return values;
    }

    /**
     * Unwraps a finished run, rethrowing what it threw.
     *
     * <p>A case that fails is a real result -- a tube that did not burst, a material that was
     * refused -- and it should reach the caller as the exception it was, not wrapped in an
     * {@link ExecutionException} whose message is the word "ExecutionException".
     */
    private static <R> R join(Future<R> future) {
        try {
            return future.get();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new IllegalStateException("the sweep was interrupted", e);
        } catch (ExecutionException e) {
            final Throwable cause = e.getCause();
            if (cause instanceof RuntimeException r) throw r;
            if (cause instanceof Error err) throw err;
            throw new IllegalStateException("a sweep point failed", cause);
        }
    }
}
