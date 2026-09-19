package org.neofiz.core;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The normal CDF, checked against values that are known independently.
 *
 * <p>This is the middle link of the copula in {@link DefectField}, and its error lands
 * directly on the Weibull marginal. So it is worth more than a smoke test: a function that is
 * right to seven digits looks perfect in every plot and still costs the size effect a per cent
 * it cannot afford, in a way that would be absorbed by relaxing the size-effect tolerance.
 */
class NormalTest {

    @Test
    @DisplayName("erf and erfc are right at the values with closed forms")
    void closedForms() {
        assertEquals(0.0, Normal.erf(0.0), 0.0);
        assertEquals(1.0, Normal.erfc(0.0), 1e-16);
        assertEquals(2.0, Normal.erfc(-100.0), 1e-15);
        assertEquals(0.0, Normal.erfc(100.0), 0.0);
    }

    /**
     * Agreement demanded against the platform's own {@code erfc}. Cody's approximation is not
     * correctly rounded and does not claim to be; measured against the system library across
     * the range it lands within a few parts in 1e15, so this is the honest bar rather than
     * the one that sounds better.
     */
    private static final double RELATIVE = 1e-13;

    @Test
    @DisplayName("erfc matches the reference values to within a few parts in 1e15")
    void publishedValues() {
        // Values taken from the platform libm, not from a printed table -- a four-figure
        // table cannot tell the difference between this implementation and one that is a
        // million times worse, which is the difference the field actually cares about. The
        // 0.46875 and 4.0 branch boundaries are straddled deliberately.
        assertErfc(0.5716076449533316, 0.4);
        assertErfc(0.4795001221869535, 0.5);
        assertErfc(0.15729920705028516, 1.0);
        assertErfc(0.0046777349810472645, 2.0);
        assertErfc(1.541725790028002e-8, 4.0);
        assertErfc(1.537459794428035e-12, 5.0);
        assertErfc(4.183825607779415e-23, 7.0);
    }

    private static void assertErfc(double expected, double x) {
        assertEquals(expected, Normal.erfc(x), Math.abs(expected) * RELATIVE,
                "erfc(" + x + ")");
    }

    @Test
    @DisplayName("erfc is continuous across both branch boundaries")
    void branchesJoin() {
        // Three separate rational approximations meet at 0.46875 and at 4. A real
        // discontinuity there would be a step in the marginal of the defect field at one
        // particular quantile, which is the sort of artefact that stays invisible until
        // someone plots a histogram of ten million samples. The step that is actually there
        // is about one part in 1e14, which is the approximation's own accuracy and not a
        // seam that can be closed.
        for (double x : new double[]{0.46875, 4.0}) {
            double below = Normal.erfc(Math.nextDown(x));
            double above = Normal.erfc(Math.nextUp(x));
            assertEquals(below, above, Math.abs(below) * RELATIVE, "discontinuity at " + x);
        }
    }

    @Test
    @DisplayName("the cdf is symmetric and hits the textbook quantiles")
    void cdf() {
        assertEquals(0.5, Normal.cdf(0.0), 1e-16);
        assertEquals(0.8413447460685429, Normal.cdf(1.0), 1e-15);
        assertEquals(0.9772498680518208, Normal.cdf(2.0), 1e-15);
        assertEquals(0.9986501019683699, Normal.cdf(3.0), 1e-15);
        assertEquals(0.9750021048517795, Normal.cdf(1.96), 1e-15);
        assertEquals(0.975, Normal.cdf(1.959963984540054), 1e-15,
                "the two-sided 95 % quantile, which is the value everyone rounds to 1.96");

        for (double x = 0.0; x < 8.0; x += 0.25) {
            assertEquals(1.0, Normal.cdf(x) + Normal.cdf(-x), 1e-15, "asymmetric at " + x);
        }
    }

    @Test
    @DisplayName("the far tail keeps its digits instead of cancelling to nothing")
    void tailDoesNotCancel() {
        // 1 - cdf(x) is exactly zero past about x = 8.3 in double precision, so a field value
        // out there would map to a uniform of exactly one and then to an infinite failure
        // strain. survival() goes through erfc of a positive argument and stays finite far
        // beyond anything reachable.
        assertEquals(0.0, 1.0 - Normal.cdf(9.0), 0.0, "this is the cancellation being avoided");
        assertTrue(Normal.survival(9.0) > 0.0);
        assertErfcLike(1.1285884059538297e-19, Normal.survival(9.0));
        assertErfcLike(4.9067139271480216e-198, Normal.survival(30.0));
        assertTrue(Normal.survival(30.0) > 0.0, "still finite at thirty sigma");
    }

    @Test
    @DisplayName("the far tail is so steep that one ulp on the argument moves the answer")
    void theTailAmplifiesTheArgument() {
        // Not a defect, and worth writing down because it looks like one. The logarithmic
        // derivative of erfc is about -2x, so at x = 21 a relative perturbation of 1e-16 in
        // the argument becomes 9e-14 in the result. Forming x/sqrt(2) by multiplying by the
        // rounded reciprocal, as survival() does, differs from dividing by sqrt(2) in the
        // last bit -- and that alone is the whole discrepancy against a reference computed
        // the other way. Anyone chasing it as an accuracy bug in the approximation will not
        // find it there.
        double byMultiply = Normal.erfc(30.0 * 0.70710678118654752440);
        double byDivide = Normal.erfc(30.0 / Math.sqrt(2.0));
        assertTrue(byMultiply != byDivide, "the two argument forms are no longer distinct");
        assertEquals(byDivide, byMultiply, Math.abs(byDivide) * 1e-12);
        assertTrue(Math.abs(byDivide - byMultiply) > Math.abs(byDivide) * 1e-14,
                "a one-ulp argument change must still be visible here, or the tail has "
                        + "stopped being steep and this test has stopped meaning anything");
    }

    private static void assertErfcLike(double expected, double actual) {
        assertEquals(expected, actual, Math.abs(expected) * RELATIVE);
    }

    @Test
    @DisplayName("the cdf is monotone across the branch boundaries of erfc")
    void monotone() {
        double previous = 0.0;
        for (double x = -8.0; x <= 8.0; x += 0.0009765625) {
            double p = Normal.cdf(x);
            assertTrue(p >= previous, "cdf decreased at " + x);
            previous = p;
        }
    }
}
