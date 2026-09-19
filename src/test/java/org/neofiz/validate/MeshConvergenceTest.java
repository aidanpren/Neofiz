package org.neofiz.validate;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Tests for Richardson extrapolation, driven by sequences whose answer is known exactly.
 *
 * <p>This is a tool for deciding how much to believe every other number in the project, so it
 * is tested against constructed sequences rather than against solver output: sample a function
 * with a known convergence order at three mesh sizes and the extrapolation has to recover the
 * limit to machine precision. If it cannot do that on {@code 10 - 2h^2}, nothing it says about
 * a mushroom diameter is worth reading.
 */
class MeshConvergenceTest {

    @Test
    @DisplayName("a second-order sequence extrapolates to its exact limit")
    void secondOrder() {
        // f(h) = 10 - 2h^2 at h = 4, 2, 1.
        MeshConvergence c = MeshConvergence.halving(10 - 2 * 16, 10 - 2 * 4, 10 - 2 * 1);

        assertEquals(2.0, c.observedOrder(), 1e-12, "order");
        assertEquals(10.0, c.extrapolated(), 1e-12, "the limit is exact, not approximate");
        assertTrue(c.isConverging());
        assertTrue(c.isMonotone());
    }

    @Test
    @DisplayName("a first-order sequence extrapolates to its exact limit")
    void firstOrder() {
        // f(h) = 10 - 2h at h = 4, 2, 1. The rate a quantity read next to a geometric
        // singularity shows, which is the case this tool exists to handle.
        MeshConvergence c = MeshConvergence.halving(2.0, 6.0, 8.0);

        assertEquals(1.0, c.observedOrder(), 1e-12);
        assertEquals(10.0, c.extrapolated(), 1e-12);
    }

    @Test
    @DisplayName("the order is recovered at refinement ratios other than two")
    void otherRefinementRatios() {
        // f(h) = 5 + h^2 at h = 9, 3, 1, refining by three each time.
        MeshConvergence c = new MeshConvergence(5 + 81, 5 + 9, 5 + 1, 3.0);

        assertEquals(2.0, c.observedOrder(), 1e-12);
        assertEquals(5.0, c.extrapolated(), 1e-12);
    }

    @Test
    @DisplayName("the uncertainty band contains the error it is a band on")
    void theBandContainsTheError() {
        // The property that makes a GCI worth quoting. Asserted rather than assumed, because
        // a band that is narrower than the error it describes is worse than no band: it turns
        // an unknown into a confident wrong number.
        MeshConvergence c = MeshConvergence.halving(10 - 2 * 16, 10 - 2 * 4, 10 - 2 * 1);

        double actualError = Math.abs(c.extrapolated() - c.fine()) / Math.abs(c.fine());
        assertTrue(c.gci() > actualError,
                "GCI of " + c.gci() + " does not cover the actual relative error of "
                        + actualError);
        assertEquals(c.gci() * Math.abs(c.fine()), c.uncertainty(), 1e-12,
                "the absolute band must agree with the fractional one");
    }

    @Test
    @DisplayName("an oscillatory sequence is refused, not extrapolated")
    void oscillatoryIsRefused() {
        // The formula would happily return a number here. A sequence that steps up then down
        // has no order to observe, and quoting one would invent a convergence that is not
        // happening.
        MeshConvergence c = MeshConvergence.halving(10.0, 12.0, 11.0);

        assertFalse(c.isMonotone());
        assertFalse(c.isConverging());
        IllegalStateException e = assertThrows(IllegalStateException.class, c::extrapolated);
        assertTrue(e.getMessage().contains("oscillatory"), "got: " + e.getMessage());
        assertThrows(IllegalStateException.class, c::observedOrder);
        assertThrows(IllegalStateException.class, c::gci);
    }

    @Test
    @DisplayName("a sequence whose steps are growing is refused")
    void divergingIsRefused() {
        // Monotone but getting worse under refinement. The observed order comes out negative
        // and the extrapolation runs away from the answer rather than towards it.
        MeshConvergence c = MeshConvergence.halving(10.0, 11.0, 13.0);

        assertTrue(c.isMonotone());
        assertFalse(c.isConverging());
        IllegalStateException e = assertThrows(IllegalStateException.class, c::extrapolated);
        assertTrue(e.getMessage().contains("not converging"), "got: " + e.getMessage());
    }

    @Test
    @DisplayName("an already-converged quantity is reported, not extrapolated")
    void anAlreadyConvergedQuantityIsRefused() {
        // The Taylor final length ratio at 8 / 16 / 32. It moves one part in ten thousand per
        // refinement, which is residual dynamics rather than a mesh trend -- but the steps are
        // very slightly shrinking, so a shrinking-steps test alone lets it through and
        // Richardson relocates it to 0.9056, twelve times further than the whole study moved
        // it. The extrapolation is not merely imprecise here, it is confidently wrong.
        MeshConvergence c = MeshConvergence.halving(0.908089322, 0.907987693, 0.907890227);

        assertTrue(c.isMonotone());
        assertTrue(c.isConverging(), "the steps really are shrinking, which is the trap");
        assertFalse(c.isAsymptotic(), "but not enough to extrapolate on");

        IllegalStateException e = assertThrows(IllegalStateException.class, c::extrapolated);
        assertTrue(e.getMessage().contains("asymptotic"), "got: " + e.getMessage());

        assertTrue(c.summary("").contains("converged"),
                "the summary should report it as converged rather than invent a limit: "
                        + c.summary(""));
        assertTrue(c.lastChange() < 2e-4, "last step " + c.lastChange());
    }

    @Test
    @DisplayName("the refinement ratio has to be a refinement")
    void validation() {
        assertThrows(IllegalArgumentException.class,
                () -> new MeshConvergence(1, 2, 3, 1.0));
        assertThrows(IllegalArgumentException.class,
                () -> new MeshConvergence(1, 2, 3, 0.5));
    }

    // ---------------------------------------------------------- the measured Taylor sequences
    //
    // Pinned as literals rather than re-run. Each of these took three meshes up to 10,240
    // elements, and the conclusions drawn from them are load-bearing enough that they should
    // be asserted on every build rather than only when someone has 40 seconds to spare.

    @Test
    @DisplayName("the Taylor mushroom converges at first order, not the element's second")
    void theMushroomConvergesAtFirstOrder() {
        // Frictionless, Johnson-Cook, reduced integration, at 8 / 16 / 32 across the radius.
        // A bilinear quad gives second order on a smooth field, so an observed order near one
        // says this quantity is being read beside a singularity -- the contact edge, where the
        // free surface meets the anvil -- and refinement will not restore the missing order.
        // That is what makes the extrapolated value the right thing to quote.
        MeshConvergence c = MeshConvergence.halving(9.7758, 9.9653, 10.0545);

        assertTrue(c.observedOrder() > 0.9 && c.observedOrder() < 1.3,
                "observed order " + c.observedOrder() + "; if this has risen towards two the "
                        + "singularity argument no longer holds and the mushroom should simply "
                        + "be refined instead of extrapolated");
        assertEquals(10.134, c.extrapolated(), 0.005);
        assertTrue(c.gci() < 0.02, "the band should be tight enough to be useful");
    }

    @Test
    @DisplayName("gripping the face restores the convergence order, which confirms the diagnosis")
    void grippingTheFaceRestoresTheOrder() {
        // The same quantity at mu = 0.20, where the impact face is largely held. If the slow
        // convergence were a property of the mesh or of the mushroom as a feature, this would
        // look the same. It does not: the order climbs back towards two and the band collapses
        // by an order of magnitude, because material is no longer sliding past the corner that
        // was generating the singular field.
        MeshConvergence gripped = MeshConvergence.halving(8.9872, 9.0263, 9.0373);
        MeshConvergence sliding = MeshConvergence.halving(9.7758, 9.9653, 10.0545);

        assertTrue(gripped.observedOrder() > 1.5,
                "a gripped face should converge close to the element's own order; got "
                        + gripped.observedOrder());
        assertTrue(gripped.gci() < sliding.gci() / 5.0,
                "gripping should tighten the band by a lot: " + gripped.gci() + " against "
                        + sliding.gci());
    }

    @Test
    @DisplayName("the extrapolated sweep crosses the measurement at a plausible friction")
    void theExtrapolatedSweepCrossesTheMeasurement() {
        // mu = 0.09, extrapolated, against the measured 9.5 mm. The number that matters is not
        // that it matches -- a coefficient was chosen to make it match -- but that the
        // coefficient which does so is one a real steel interface could have. The un-converged
        // n = 8 mesh asked for 0.05, which was not.
        MeshConvergence c = MeshConvergence.halving(9.2992, 9.4056, 9.4584);

        assertEquals(9.5, c.extrapolated(), c.uncertainty(),
                "the extrapolated value should land on the measurement within its own band");
        assertTrue(c.observedOrder() > 0.9 && c.observedOrder() < 1.3,
                "the same first-order behaviour should hold at this coefficient too");
    }
}
