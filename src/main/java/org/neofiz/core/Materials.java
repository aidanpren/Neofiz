package org.neofiz.core;

import java.util.ArrayList;
import java.util.List;
import java.util.function.IntFunction;

/**
 * Which material each element is made of.
 *
 * <p>Until now a solve was one mesh of one substance, which is the right shape for a
 * validation harness -- a burst test is a copper tube and nothing else -- and the wrong shape
 * for a sandbox, where the interesting question is what a steel bar does to a concrete slab.
 * This is the indirection that removes the restriction: everything that used to read a
 * constant now reads {@code at(element)}.
 *
 * <p>The uniform case is not a degenerate one-entry array. It holds a single reference and
 * answers every element with it, so a million-element mesh of one substance costs one pointer
 * and the old arithmetic is reproduced <em>bit for bit</em> -- the same {@code Material}
 * object reaches the same kernel by a different route. That matters more than the memory: it
 * is what lets every existing gate keep passing unchanged, which is the only evidence that
 * this refactor did not move an answer.
 */
public final class Materials {

    /** The single material, or null when this map is per-element. */
    private final Material one;
    /** One material per element, or null when this map is uniform. */
    private final Material[] each;

    private Materials(Material one, Material[] each) {
        this.one = one;
        this.each = each;
    }

    /** One substance everywhere, on a mesh of any size. */
    public static Materials uniform(Material material) {
        if (material == null) throw new IllegalArgumentException("material must not be null");
        return new Materials(material, null);
    }

    /** One entry per element, in element order. The array is copied. */
    public static Materials of(Material... byElement) {
        if (byElement.length == 0) {
            throw new IllegalArgumentException("need at least one material");
        }
        final Material[] copy = byElement.clone();
        for (int e = 0; e < copy.length; e++) {
            if (copy[e] == null) {
                throw new IllegalArgumentException("element " + e + " has no material");
            }
        }
        return new Materials(null, copy);
    }

    /**
     * One entry per element, chosen by a rule. The usual way to paint a mesh: the rule reads
     * the element's centroid off the mesh and decides which body it belongs to.
     */
    public static Materials byElement(int elementCount, IntFunction<Material> pick) {
        if (elementCount < 1) throw new IllegalArgumentException("need at least one element");
        final Material[] out = new Material[elementCount];
        for (int e = 0; e < elementCount; e++) out[e] = pick.apply(e);
        return of(out);
    }

    /** The material of one element. */
    public Material at(int element) {
        return one != null ? one : each[element];
    }

    /** Whether one substance covers the whole mesh. */
    public boolean isUniform() {
        return one != null;
    }

    /** How many elements this map covers, or -1 when it is uniform and so covers any mesh. */
    public int elementCount() {
        return each == null ? -1 : each.length;
    }

    /**
     * Refuses a map that does not fit the mesh it is about to be used on.
     *
     * <p>A short map would silently throw on the first element past its end, in a parallel
     * worker, several hundred timesteps in. Checking at construction turns that into a
     * sentence.
     */
    public void requireCovers(int elementCount) {
        if (each != null && each.length != elementCount) {
            throw new IllegalArgumentException(
                    "materials cover " + each.length + " elements but the mesh has "
                            + elementCount);
        }
    }

    /** The distinct materials present, in first-appearance order. Compared by identity. */
    public List<Material> distinct() {
        final List<Material> out = new ArrayList<>();
        if (one != null) {
            out.add(one);
            return out;
        }
        for (Material m : each) {
            boolean seen = false;
            for (Material k : out) if (k == m) { seen = true; break; }
            if (!seen) out.add(m);
        }
        return out;
    }

    /**
     * The fastest dilatational wave anywhere on the mesh, m/s.
     *
     * <p>The CFL condition is per element -- a step is stable when no element's wave crosses
     * its own shortest edge -- but the solver carries one global step, so it has to assume the
     * worst pairing: the shortest edge on the mesh and the fastest material on it, whether or
     * not they are the same element. Conservative by construction, and the cost of that
     * conservatism is exactly what mass scaling is for later.
     */
    public double fastestWave() {
        double fastest = 0.0;
        if (one != null) return one.dilatationalWaveSpeed();
        for (Material m : each) {
            final double c = m.dilatationalWaveSpeed();
            if (c > fastest) fastest = c;
        }
        return fastest;
    }
}
