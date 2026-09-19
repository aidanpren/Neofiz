package org.neofiz.solver;

import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.LockSupport;

/**
 * A fixed pool of workers that splits an index range across threads and waits for all of
 * them, with the master thread taking the first chunk itself.
 *
 * <p>Built rather than taken from {@code java.util.concurrent} for one reason: the granularity.
 * An explicit run of the plan's worked example is 234,000 steps, and each step has two
 * parallel regions, so a dispatch that costs 20 microseconds -- an ordinary figure for
 * {@code ExecutorService.invokeAll}, which parks and unparks a thread per task -- would add
 * nine seconds to a run whose entire budget is four. The useful work in one region is around
 * 200 microseconds, so the synchronisation has to cost single-digit microseconds to be
 * invisible, and that means the workers do not sleep between regions.
 *
 * <h2>The barrier</h2>
 *
 * A ticket counter the master increments to publish work, and a completion counter the
 * workers increment when they finish. Both are atomics, which supplies the happens-before
 * edge: everything the master wrote before incrementing the ticket is visible to a worker
 * that reads it, and everything a worker wrote before incrementing the completion count is
 * visible to the master that reads it. No other synchronisation is needed and no field
 * carrying the work itself has to be volatile.
 *
 * <p>Workers spin on the ticket. Spinning burns a core while idle, which is the right trade
 * for an offline solver on a machine that is running nothing else -- the whole architecture
 * is built on the decision not to solve in a frame budget -- but it is a bad neighbour if the
 * pool is left alive and unused, so the spin decays: a tight {@link Thread#onSpinWait()} loop
 * first, then {@link Thread#yield()}, then a short park. The fast path is the first one and
 * costs about a microsecond.
 *
 * <h2>Why the schedule is dynamic, and why it is allowed to be</h2>
 *
 * Chunks are handed out from a shared cursor rather than divided up in advance, so a thread
 * that finishes early comes back for more. That is normally a trade -- better balance, at the
 * cost of an answer that depends on which thread got which work -- and normally the trade is
 * refused in a solver that has to be reproducible.
 *
 * <p>Here there is no trade, because {@link ExplicitSolver} does not let a chunk produce a
 * sum. Each element writes its own slot and a separate pass gathers them in element order, so
 * the result is independent of the partition entirely: equal chunks, unequal chunks, one
 * thread or twelve, all produce the same bits. The determinism is a property of the data
 * layout, not of the schedule, which leaves the schedule free to be chosen purely for speed.
 *
 * <p>And it needs to be. Three separate things make an equal split the wrong one: elements
 * that are yielding cost several times an elastic one and the plastic zone moves during a
 * run; the operating system preempts one worker and not the others; and on a machine with
 * both performance and efficiency cores -- every current Apple part, and increasingly
 * everything else -- an equal split runs at the speed of the slowest core.
 *
 * <p>That third one is not a small effect and it is not hypothetical. Measured here, on four
 * performance cores and six efficiency ones, an equal split of the reduced kernel:
 *
 * <pre>
 *    threads     equal split   shared cursor
 *      2            0.87x          2.03x
 *      8            1.79x          3.74x
 *     10            2.20x          3.94x
 * </pre>
 *
 * <p>The first row is the one worth keeping. <b>Two threads with an equal split are slower
 * than one thread</b>, because the second half of the work lands on a core that takes longer
 * to do half of it than the first core takes to do all of it, and the barrier waits. Nothing
 * about that is visible in a profile of the kernel; it is entirely a property of the schedule.
 */
final class Parallel implements AutoCloseable {

    /** What a worker runs over its slice. The worker index is its slot in any scratch array. */
    interface Kernel {
        void run(int from, int to, int worker);
    }

    /** Tight spins before yielding. About a microsecond on current hardware. */
    private static final int SPINS_BEFORE_YIELD = 2000;
    /** Yields before parking. Beyond this the pool is idle rather than between regions. */
    private static final int YIELDS_BEFORE_PARK = 100;
    private static final long PARK_NANOS = 50_000L;

    /**
     * Chunks aimed at per thread. Too few and a straggler holds up the barrier; too many and
     * the cursor is contended for no benefit. Eight leaves enough tail to absorb a core
     * running at half speed while costing one atomic per few hundred elements.
     */
    private static final int CHUNKS_PER_THREAD = 8;
    /** Smallest chunk worth an atomic and a call. */
    private static final int MIN_GRAIN = 16;

    private final int threads;
    private final Thread[] workers;

    private final AtomicInteger ticket = new AtomicInteger();
    private final AtomicInteger completed = new AtomicInteger();
    private final AtomicInteger cursor = new AtomicInteger();

    /**
     * Published under the ticket's release, read under its acquire. Plain fields: the atomic
     * carries them across, and making them volatile would add a fence per read for nothing.
     */
    private Kernel kernel;
    private int count;
    private int grain;
    private volatile boolean shutdown;

    Parallel(int threads) {
        if (threads < 1) throw new IllegalArgumentException("need at least one thread");
        this.threads = threads;
        this.workers = new Thread[threads - 1];
        for (int i = 0; i < workers.length; i++) {
            final int worker = i + 1;
            Thread t = new Thread(() -> loop(worker), "neofiz-worker-" + worker);
            t.setDaemon(true);
            workers[i] = t;
            t.start();
        }
    }

    int threads() {
        return threads;
    }

    /**
     * Runs {@code kernel} over {@code [0, count)} in chunks taken from a shared cursor, and
     * returns once every chunk has finished. The calling thread takes chunks too rather than
     * waiting, so a one-thread pool is the kernel called directly.
     */
    void run(int count, Kernel kernel) {
        if (threads == 1) {
            kernel.run(0, count, 0);
            return;
        }
        this.kernel = kernel;
        this.count = count;
        this.grain = Math.max(MIN_GRAIN, count / (threads * CHUNKS_PER_THREAD));
        cursor.set(0);
        completed.set(0);
        ticket.incrementAndGet();      // release: publishes kernel, count and grain

        drain(0);

        int spins = 0;
        while (completed.get() < threads - 1) spins = wait(spins);
    }

    /** Takes chunks until there are none left. */
    private void drain(int worker) {
        final int n = count;
        final int g = grain;
        final Kernel k = kernel;
        while (true) {
            final int from = cursor.getAndAdd(g);
            if (from >= n) return;
            k.run(from, Math.min(from + g, n), worker);
        }
    }

    private void loop(int worker) {
        int seen = 0;
        int spins = 0;
        while (true) {
            final int now = ticket.get();   // acquire: sees kernel, count and grain
            if (now == seen) {
                if (shutdown) return;
                spins = wait(spins);
                continue;
            }
            seen = now;
            spins = 0;
            try {
                drain(worker);
            } finally {
                completed.incrementAndGet();   // release: publishes everything written above
            }
        }
    }

    private static int wait(int spins) {
        if (spins < SPINS_BEFORE_YIELD) {
            Thread.onSpinWait();
        } else if (spins < SPINS_BEFORE_YIELD + YIELDS_BEFORE_PARK) {
            Thread.yield();
        } else {
            LockSupport.parkNanos(PARK_NANOS);
            return spins;
        }
        return spins + 1;
    }

    @Override
    public void close() {
        shutdown = true;
        for (Thread t : workers) {
            LockSupport.unpark(t);
            try {
                t.join(1000);
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                return;
            }
        }
    }
}
