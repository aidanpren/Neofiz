package org.neofiz.core;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The Johnson-Cook flow stress as a pure function, tested away from the return map.
 *
 * <p>Worth separating, because the three terms have limits and clamps that are easy to get
 * wrong in ways the Taylor case would absorb silently: a rate term that softens instead of
 * hardening, or a thermal term that goes negative past melt, produces a plausible-looking
 * mushroom rather than an obvious failure.
 */
class JohnsonCookTest {

    private static final JohnsonCook STEEL = JohnsonCook.STEEL_4340;
    private static final double ROOM = 293.0;

    @Test
    @DisplayName("at reference conditions the flow stress is exactly A")
    void referenceConditionsGiveA() {
        // Zero plastic strain, at or below the reference rate, at room temperature: every
        // term is on its identity and the quasi-static yield stress is left. This is the
        // value Material keeps its yieldStress in step with, so the elastic check and the
        // return map agree on where the surface starts.
        assertEquals(STEEL.a(), STEEL.flowStress(0.0, 0.0, ROOM), 0.0);
        assertEquals(STEEL.a(), STEEL.flowStress(0.0, STEEL.referenceStrainRate(), ROOM), 0.0);
    }

    @Test
    @DisplayName("strain hardening is A + B eps^n")
    void strainTerm() {
        for (double eps : new double[]{0.01, 0.1, 0.5, 1.0, 2.0}) {
            assertEquals(STEEL.a() + STEEL.b() * Math.pow(eps, STEEL.n()),
                    STEEL.strainTerm(eps), 1e-6);
        }
        assertEquals(STEEL.a() + STEEL.b(), STEEL.strainTerm(1.0), 1e-6,
                "at unit plastic strain the hardening term is exactly B");
    }

    @Test
    @DisplayName("the strain slope is the derivative of the strain term")
    void strainSlopeIsTheDerivative() {
        // The return map's Newton step is only as good as this, and an analytic derivative
        // that disagrees with its own function is the classic cause of a solve that converges
        // to slightly the wrong root -- slowly enough to look like a tolerance problem.
        for (double eps : new double[]{0.05, 0.2, 0.6, 1.5}) {
            double d = eps * 1e-7;
            double numeric = (STEEL.strainTerm(eps + d) - STEEL.strainTerm(eps - d)) / (2 * d);
            assertEquals(numeric, STEEL.strainSlope(eps), Math.abs(numeric) * 1e-6,
                    "analytic and numeric slope disagree at eps = " + eps);
        }
    }

    @Test
    @DisplayName("the pow-free slope is the same number, not merely a close one")
    void strainSlopeFromIsAlgebra() {
        // The return map takes the slope from the strain term it has already computed --
        // B n eps^(n-1) = n (B eps^n) / eps -- which halves the transcendentals in the hot
        // loop. That is an identity, so it is held to round-off rather than to a tolerance,
        // and over a range wide enough to include the cancellation regime where B eps^n is
        // small against A.
        for (JohnsonCook law : new JohnsonCook[]{JohnsonCook.STEEL_4340, JohnsonCook.COPPER_OFHC}) {
            for (double eps = 1.0e-10; eps < 2.0; eps *= 1.7) {
                double exact = law.strainSlope(eps);
                double cheap = law.strainSlopeFrom(law.strainTerm(eps), eps);
                assertEquals(exact, cheap, Math.abs(exact) * 1e-12,
                        "slope forms disagree at eps = " + eps);
            }
        }
    }

    @Test
    @DisplayName("the strain slope is infinite at zero plastic strain, and says so")
    void strainSlopeIsSingularAtZero() {
        // n < 1, so the first increment of plastic strain is infinitely stiff. Reported
        // rather than quietly floored, because the return map has to safeguard against it
        // and cannot do that if the law hides it.
        assertTrue(Double.isInfinite(STEEL.strainSlope(0.0)),
                "a floored slope here would let an unguarded Newton step look reasonable");
    }

    @ParameterizedTest
    @ValueSource(doubles = {0.0, 1e-12, 0.001, 0.5, 1.0})
    @DisplayName("rate hardening is clamped at the reference rate, never softening")
    void rateTermIsClamped(double rate) {
        // Unclamped, 1 + C ln(epsdot) passes through zero strength at about 1e-31 per second
        // and goes negative below it. Unreachable physically, immediately reachable in the
        // first Newton iteration where the trial rate is zero.
        assertEquals(1.0, STEEL.rateTerm(rate), 0.0,
                "the rate term must not soften below the reference rate");
        assertTrue(STEEL.rateIsClamped(rate));
    }

    @Test
    @DisplayName("rate hardening is worth about 16 % at Taylor-impact rates")
    void rateTermAtImpactRates() {
        // The number the whole rate term exists for. 1e4 to 1e5 per second is what a Taylor
        // mushroom imposes.
        assertEquals(1.0 + 0.014 * Math.log(1e4), STEEL.rateTerm(1e4), 1e-12);
        assertEquals(1.129, STEEL.rateTerm(1e4), 5e-4);
        assertEquals(1.161, STEEL.rateTerm(1e5), 5e-4);
        assertTrue(!STEEL.rateIsClamped(1e4));
    }

    @Test
    @DisplayName("thermal softening runs from 1 at room temperature to 0 at melt")
    void thermalTerm() {
        assertEquals(1.0, STEEL.thermalTerm(ROOM), 0.0);
        assertEquals(0.0, STEEL.thermalTerm(STEEL.meltTemperature()), 0.0);

        // Monotone in between, and floored past melt rather than going negative -- a
        // negative flow stress is a material that resists deformation backwards.
        double previous = 1.0;
        for (double t = ROOM; t <= 2200.0; t += 50.0) {
            double theta = STEEL.thermalTerm(t);
            assertTrue(theta >= 0.0, "thermal term went negative at " + t + " K");
            assertTrue(theta <= previous + 1e-12, "thermal term rose with temperature at " + t);
            previous = theta;
        }
        assertEquals(0.0, STEEL.thermalTerm(3000.0), 0.0, "past melt there is no strength");
        assertEquals(1.0, STEEL.thermalTerm(200.0), 0.0, "below room, no thermal hardening");
    }

    @Test
    @DisplayName("temperature comes from plastic work, adiabatically")
    void temperatureFromWork() {
        final double density = 7850.0;
        // 1 GJ/m^3 of plastic work, which is about what this test's mushroom reaches.
        final double work = 1.0e9;
        final double expected = ROOM + 0.9 * work / (density * 477.0);

        assertEquals(expected, STEEL.temperature(work, density), 1e-9);
        assertEquals(ROOM, STEEL.temperature(0.0, density), 0.0,
                "an undeformed point must sit at room temperature");
        assertEquals(STEEL.temperature(work, density),
                ROOM + STEEL.thermalCoefficient(density) * work, 1e-9,
                "the precomputed coefficient must agree with the full expression");
    }

    @Test
    @DisplayName("rate hardening and thermal softening nearly cancel at Taylor-impact conditions")
    void theTwoTermsNearlyCancel() {
        // This is the substantive physical finding of the M1 work and it belongs in a test,
        // because it is the reason adding Johnson-Cook barely moved the mushroom. At the
        // hottest, most strained point of the reference case -- about 517 K and 0.75 plastic
        // strain at a rate near 4e4 per second -- the two corrections are the same size and
        // opposite in sign.
        final double rate = 0.75 / 20e-6;
        final double temperature = 517.0;

        double hardening = STEEL.rateTerm(rate);
        double softening = STEEL.thermalTerm(temperature);

        assertTrue(hardening > 1.14 && hardening < 1.16,
                "rate hardening was " + hardening + ", expected about 1.147");
        assertTrue(softening > 0.85 && softening < 0.87,
                "thermal softening was " + softening + ", expected about 0.859");
        assertEquals(1.0, hardening * softening, 0.03,
                "the product was " + hardening * softening + ". If these stop cancelling, the "
                        + "explanation for why Johnson-Cook left the mushroom almost unchanged "
                        + "no longer holds and the M1 write-up needs revisiting");
    }

    @Test
    @DisplayName("nonsense constants are refused")
    void validation() {
        assertThrows(IllegalArgumentException.class, () -> new JohnsonCook(
                -1.0, 510e6, 0.26, 0.014, 1.03, 1.0, 293, 1793, 477, 0.9), "negative A");
        assertThrows(IllegalArgumentException.class, () -> new JohnsonCook(
                792e6, 510e6, 1.5, 0.014, 1.03, 1.0, 293, 1793, 477, 0.9), "n above 1");
        assertThrows(IllegalArgumentException.class, () -> new JohnsonCook(
                792e6, 510e6, 0.26, 0.014, 1.03, 1.0, 1793, 293, 477, 0.9), "melt below room");
        assertThrows(IllegalArgumentException.class, () -> new JohnsonCook(
                792e6, 510e6, 0.26, 0.014, 1.03, 0.0, 293, 1793, 477, 0.9), "zero ref rate");
        assertThrows(IllegalArgumentException.class, () -> new JohnsonCook(
                792e6, 510e6, 0.26, 0.014, 1.03, 1.0, 293, 1793, 477, 1.5), "beta above 1");
    }

    @Test
    @DisplayName("a material refuses a yield stress that disagrees with its Johnson-Cook A")
    void materialKeepsTheSurfacesInStep() {
        // The elastic trial checks against Material.yieldStress and the return map solves
        // against the law. If those disagree, a point can be flagged as yielding by one
        // surface and solved against another, which shows up as a tiny spurious plastic
        // strain everywhere rather than as an error.
        Material good = Material.STEEL_4340.withJohnsonCook(STEEL);
        assertEquals(STEEL.a(), good.yieldStress(), 0.0);
        assertEquals(0.0, good.isotropicHardening(), 0.0,
                "linear hardening must be cleared, not left to read as if it still applied");
        assertTrue(good.isRateDependent());

        assertThrows(IllegalArgumentException.class, () -> new Material(
                "mismatched", 205e9, 0.29, 7850.0, 500.0e6, 0.0, 0.0, STEEL));
    }
}
