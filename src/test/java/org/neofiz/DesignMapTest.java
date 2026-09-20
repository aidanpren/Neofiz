package org.neofiz;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * The one piece of arithmetic in the design map.
 *
 * <p>Every cell of the map is {@link DesignMap#burstFraction}, so a sign or an off-by-one
 * here would not produce a broken picture -- it would produce a plausible one, shifted by a
 * row or inverted, and the map's entire purpose is to be believed at a glance.
 */
class DesignMapTest {

    /** A sorted batch, in pascals, of the shape a column of the map actually holds. */
    private static final double[] BATCH = {10.0e6, 10.2e6, 10.4e6, 10.6e6, 10.8e6};

    @Test
    @DisplayName("below the weakest tube nothing bursts, above the strongest everything does")
    void endsAreSaturated() {
        assertEquals(0.0, DesignMap.burstFraction(BATCH, 9.0e6), 0.0);
        assertEquals(1.0, DesignMap.burstFraction(BATCH, 12.0e6), 0.0);
    }

    @Test
    @DisplayName("the fraction counts the tubes that did not reach this pressure")
    void countsBelow() {
        assertEquals(0.2, DesignMap.burstFraction(BATCH, 10.1e6), 1e-12);
        assertEquals(0.4, DesignMap.burstFraction(BATCH, 10.3e6), 1e-12);
        assertEquals(0.8, DesignMap.burstFraction(BATCH, 10.7e6), 1e-12);
    }

    @Test
    @DisplayName("a tube counts as holding the pressure that is its own capacity")
    void capacityIsHeld() {
        // Strictly below, because a tube's capacity is the pressure it reached. Counting it
        // as burst would shift the whole boundary down by one sample and would do it
        // invisibly.
        assertEquals(0.0, DesignMap.burstFraction(BATCH, 10.0e6), 0.0);
        assertEquals(0.2, DesignMap.burstFraction(BATCH, 10.2e6), 1e-12);
    }

    @Test
    @DisplayName("one tube gives a map with a hard edge, which is the point of using more")
    void singleSampleIsAStep() {
        final double[] one = {10.0e6};
        assertEquals(0.0, DesignMap.burstFraction(one, 9.9e6), 0.0);
        assertEquals(1.0, DesignMap.burstFraction(one, 10.1e6), 0.0);
    }

    @Test
    @DisplayName("ties are counted once each")
    void ties() {
        // Nominally identical tubes can land on the same capacity to the bit when the defect
        // draw barely moves the answer, and a duplicate must not be double-counted or skipped.
        final double[] tied = {10.0e6, 10.0e6, 10.5e6};
        assertEquals(2.0 / 3.0, DesignMap.burstFraction(tied, 10.2e6), 1e-12);
    }
}
