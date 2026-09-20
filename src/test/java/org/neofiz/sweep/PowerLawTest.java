package org.neofiz.sweep;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The log-log fit.
 *
 * <p>The fit's job is to hand back an exponent a person can read as physics, and the failure
 * mode that matters is not imprecision -- it is confidence. A power law fitted to data that
 * is not a power law returns an exponent anyway, and that exponent is worse than no answer
 * because it looks like one. So half of what is tested here is that
 * {@link PowerLaw#maxResidualPercent} actually notices.
 */
class PowerLawTest {

    private static double[] powers(double[] x, double coefficient, double exponent) {
        final double[] y = new double[x.length];
        for (int i = 0; i < x.length; i++) y[i] = coefficient * Math.pow(x[i], exponent);
        return y;
    }

    @Test
    @DisplayName("an exact power law comes back exactly")
    void recoversExactly() {
        final double[] x = Sweep.logarithmic(0.5, 5.0, 12);
        final PowerLaw law = PowerLaw.fit(x, powers(x, 9.02, 1.0));

        assertEquals(1.0, law.exponent(), 1e-12);
        assertEquals(9.02, law.coefficient(), 1e-10);
        assertTrue(law.maxResidualPercent() < 1e-10,
                "an exact law should have no residual; got " + law.maxResidualPercent());
        assertEquals(1.0, law.rSquared(), 1e-12);
        assertEquals(12, law.points());
    }

    @Test
    @DisplayName("a negative exponent is recovered as readily as a positive one")
    void negativeExponent() {
        final double[] x = Sweep.logarithmic(23.0, 198.0, 12);
        final PowerLaw law = PowerLaw.fit(x, powers(x, 897.7, -1.0));
        assertEquals(-1.0, law.exponent(), 1e-12);
        assertTrue(law.maxResidualPercent() < 1e-10);
    }

    @Test
    @DisplayName("the residual notices a law that is not a power law")
    void residualCatchesStructure() {
        // 1/(d + t) is what a tube with a fixed wall actually does when it is plotted
        // against its bore diameter rather than its mean diameter. It is nearly a power law
        // and it is not one, and the difference is the whole of the lesson in that sweep.
        final double[] d = Sweep.logarithmic(23.0, 198.0, 12);
        final double[] y = new double[d.length];
        for (int i = 0; i < d.length; i++) y[i] = 1794.0 / (d[i] + 2.0);

        final PowerLaw law = PowerLaw.fit(d, y);

        assertTrue(law.exponent() > -1.0,
                "the fit should fall short of -1; got " + law.exponent());
        assertTrue(law.maxResidualPercent() > 0.5,
                "the residual should see the structure; got " + law.maxResidualPercent() + " %");

        // And r^2 should not, which is the reason the residual is the statistic reported
        // first. Over this range it stays above 0.999 on data the exponent is visibly wrong
        // about.
        assertTrue(law.rSquared() > 0.999,
                "r2 was expected to be unhelpfully high; got " + law.rSquared());
    }

    @Test
    @DisplayName("the same data on the right abscissa fits cleanly")
    void rightAbscissaIsClean() {
        // The companion to the test above: nothing about the numbers changed except which
        // length they are plotted against.
        final double[] d = Sweep.logarithmic(23.0, 198.0, 12);
        final double[] y = new double[d.length];
        final double[] mean = new double[d.length];
        for (int i = 0; i < d.length; i++) {
            y[i] = 1794.0 / (d[i] + 2.0);
            mean[i] = d[i] + 2.0;
        }

        final PowerLaw law = PowerLaw.fit(mean, y);
        assertEquals(-1.0, law.exponent(), 1e-12);
        assertTrue(law.maxResidualPercent() < 1e-10);
    }

    @Test
    @DisplayName("the fitted curve is the one the exponent describes")
    void atMatchesTheFit() {
        final double[] x = Sweep.logarithmic(1.0, 10.0, 6);
        final PowerLaw law = PowerLaw.fit(x, powers(x, 3.0, 0.5));
        assertEquals(3.0 * Math.sqrt(4.0), law.at(4.0), 1e-9);
    }

    @Test
    @DisplayName("data a logarithm cannot be taken of is refused, not filtered")
    void refusesNonPositive() {
        // Dropping a run silently would mean the sweep describes a different experiment than
        // the one that was asked for, and would do it without saying so.
        final IllegalArgumentException e = assertThrows(IllegalArgumentException.class,
                () -> PowerLaw.fit(new double[] {1, 2, 3}, new double[] {1, 0, 3}));
        assertTrue(e.getMessage().contains("point 1"), e.getMessage());

        assertThrows(IllegalArgumentException.class,
                () -> PowerLaw.fit(new double[] {1, -2}, new double[] {1, 2}));
    }

    @Test
    @DisplayName("a fit needs points that differ")
    void refusesDegenerateInput() {
        assertThrows(IllegalArgumentException.class,
                () -> PowerLaw.fit(new double[] {1, 2}, new double[] {1}));
        assertThrows(IllegalArgumentException.class,
                () -> PowerLaw.fit(new double[] {1}, new double[] {1}));

        final IllegalArgumentException e = assertThrows(IllegalArgumentException.class,
                () -> PowerLaw.fit(new double[] {3, 3, 3}, new double[] {1, 2, 3}));
        assertTrue(e.getMessage().contains("no sweep here"), e.getMessage());
    }

    @Test
    @DisplayName("the description reads as a statement about physics")
    void describes() {
        final double[] x = Sweep.logarithmic(0.5, 5.0, 8);
        assertEquals("burst pressure proportional to wall thickness^1.00",
                PowerLaw.fit(x, powers(x, 9.0, 1.0))
                        .describe("burst pressure", "wall thickness"));
    }
}
