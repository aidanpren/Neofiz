package org.neofiz.validate;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.neofiz.core.JohnsonCook;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Tests for the burst closed forms, against answers known independently of them.
 *
 * <p>Three of these are exact identities rather than tolerances, which is the point: a
 * bisection that lands on {@code n/sqrt(3)} to eight digits for four different exponents is not
 * being graded on a curve, and neither is one whose non-hardening limit reproduces a formula
 * derived somewhere else in this package from entirely different reasoning.
 */
class BurstTest {

    /**
     * A power law {@code K eps^n}, as closely as Johnson-Cook can express one.
     *
     * <p>{@code A} must be positive -- a material with no yield stress at all is not a material
     * -- so it is set to 1 Pa against a hardening term of order 1e8, eight orders down. The
     * identities below therefore hold to about 1e-8 rather than exactly, and asserting them at
     * 1e-7 is asserting that nothing but that deliberate perturbation is in the way.
     */
    private static JohnsonCook powerLaw(double n) {
        return new JohnsonCook(1.0, 292.0e6, n, 0.0, 1.09, 1.0, 293.0, 1356.0, 383.0, 0.0);
    }

    @Test
    @DisplayName("a power law bursts at exactly n / sqrt(3)")
    void powerLawInstability() {
        for (double n : new double[] {0.1, 0.2, 0.31, 0.5, 1.0}) {
            assertEquals(n / Math.sqrt(3.0), Burst.instabilityStrain(powerLaw(n)), 1e-7,
                    "biaxial instability strain at n = " + n);
        }
    }

    @Test
    @DisplayName("and necks in tension at exactly n, which is Considere")
    void powerLawTensileInstability() {
        for (double n : new double[] {0.1, 0.2, 0.31, 0.5}) {
            assertEquals(n, Burst.tensileInstabilityStrain(powerLaw(n)), 1e-7,
                    "uniaxial instability strain at n = " + n);
        }
    }

    @Test
    @DisplayName("a tube goes unstable at 1/sqrt(3) of the strain the same material necks at")
    void biaxialGoesFirst() {
        JohnsonCook law = JohnsonCook.COPPER_OFHC.quasiStatic();
        double ratio = Burst.instabilityStrain(law) / Burst.tensileInstabilityStrain(law);

        // Not exactly 1/sqrt(3) for a law with a yield stress in it -- that identity belongs
        // to the pure power law -- but the ordering is the physics and it is not close.
        assertTrue(ratio < 0.6 && ratio > 0.5,
                "biaxial tension weakens the geometry twice as fast, so the instability "
                        + "arrives far earlier; got a ratio of " + ratio);
    }

    @Test
    @DisplayName("with no hardening, burst reduces to the fully plastic limit pressure")
    void nonHardeningLimit() {
        // B = 0: nothing can outrun the geometric softening, so the tube is unstable from
        // first yield and its burst pressure is the pressure at which it goes fully plastic.
        // ElasticPlastic derives that from equilibrium plus the yield condition integrated
        // through the wall, with no thin-wall assumption anywhere; this one derives it from a
        // thin-wall force balance and an instability condition. Nothing is shared, so their
        // agreement is evidence, and the *rate* at which they agree is more evidence still:
        //
        //   limit / burst = (2 sy/sqrt(3)) ln(b/a) / ((2/sqrt(3)) sy t/r)
        //                 = 2 artanh(t/2r) / (t/r)  =  1 + (t/r)^2 / 12 + ...
        //
        // so the relative gap is not merely small, it is (t/r)^2 / 12 and must fall by exactly
        // four for every halving of the wall. That is a prediction with no free constant in
        // it, which is worth more than any tolerance.
        JohnsonCook flat = new JohnsonCook(200.0e6, 0.0, 0.31, 0.0, 1.09,
                1.0, 293.0, 1356.0, 383.0, 0.0);
        assertEquals(0.0, Burst.instabilityStrain(flat), 0.0,
                "a material that cannot harden is unstable from the first increment");

        double previous = Double.NaN;
        for (double slenderness : new double[] {25.0, 50.0, 100.0, 200.0}) {
            double meanRadius = 25.0e-3;
            double thickness = 2.0 * meanRadius / slenderness;
            double a = meanRadius - 0.5 * thickness;
            double b = meanRadius + 0.5 * thickness;
            double ratio = thickness / meanRadius;

            double mine = Burst.burstPressure(flat, thickness, meanRadius);
            double theirs = ElasticPlastic.limitPressure(a, b, flat.a());
            double gap = (theirs - mine) / theirs;

            assertTrue(gap > 0.0, "the thin-wall form must sit below the exact one, not "
                    + "either side of it; D/t = " + slenderness);
            // To within the next term of the series, which is (t/r)^4 / 80. Bounding the
            // residual by (t/r)^4 rather than by a fixed tolerance keeps this an assertion
            // about the *order* of the agreement at every slenderness, instead of one that
            // gets easier to pass as the wall thins.
            assertEquals(ratio * ratio / 12.0, gap, Math.pow(ratio, 4),
                    "the gap at D/t = " + slenderness + " is not (t/r)^2 / 12");
            if (!Double.isNaN(previous)) {
                assertEquals(4.0, previous / gap, 0.01,
                        "halving the wall must quarter the gap, at D/t = " + slenderness);
            }
            previous = gap;
        }
    }

    @Test
    @DisplayName("Barlow on the UTS is conservative, by a few per cent")
    void barlowIsConservative() {
        JohnsonCook law = JohnsonCook.COPPER_OFHC.quasiStatic();
        double effective = Burst.effectiveStress(law);
        double uts = Burst.ultimateTensileStrength(law);

        assertTrue(uts < effective,
                "the rule of thumb has to land low or it is not a safe rule of thumb");
        double gap = (effective - uts) / effective;
        assertTrue(gap > 0.02 && gap < 0.05,
                "two large errors nearly cancelling should leave a few per cent; got " + gap);

        // And the same statement in the form the formula is actually used in.
        assertEquals(Burst.barlow(effective, 2.0e-3, 50.0e-3),
                Burst.burstPressure(law, 2.0e-3, 25.0e-3), 1e-6,
                "the effective stress is defined as the one that makes Barlow exact");
    }

    @Test
    @DisplayName("no admissible law bursts beyond unit equivalent strain")
    void boundedByConstruction() {
        for (double n : new double[] {0.05, 0.31, 0.99, 1.0}) {
            for (double b : new double[] {1.0e6, 292.0e6, 50.0e9}) {
                JohnsonCook law = new JohnsonCook(90.0e6, b, n, 0.0, 1.09,
                        1.0, 293.0, 1356.0, 383.0, 0.0);
                double strain = Burst.instabilityStrain(law);
                assertTrue(strain >= 0.0 && strain < 1.0,
                        "B = " + b + ", n = " + n + " gave " + strain);
            }
        }
    }

    @Test
    @DisplayName("a rate-dependent law is refused, not approximated")
    void refusesRateDependence() {
        // The instability strain of a rate-hardening tube depends on how fast the pump runs,
        // so there is no single number for a closed form to predict.
        IllegalArgumentException e = assertThrows(IllegalArgumentException.class,
                () -> Burst.instabilityStrain(JohnsonCook.COPPER_OFHC));
        assertTrue(e.getMessage().contains("flow curve"), e.getMessage());

        assertThrows(IllegalArgumentException.class,
                () -> Burst.burstPressure(JohnsonCook.COPPER_OFHC, 2.0e-3, 25.0e-3));
    }

    @Test
    @DisplayName("the quasi-static form is a flow curve and nothing else")
    void quasiStaticIsAFlowCurve() {
        JohnsonCook law = JohnsonCook.COPPER_OFHC.quasiStatic();
        assertTrue(law.isQuasiStatic());

        for (double rate : new double[] {0.0, 1.0, 1.0e3, 1.0e6}) {
            assertEquals(law.a() + law.b() * Math.pow(0.2, law.n()),
                    law.flowStress(0.2, rate, law.roomTemperature()), 1e-9,
                    "rate " + rate + " must not move the flow stress");
        }
        // Isothermal by removing the heat source, not by removing the thermal term: with no
        // Taylor-Quinney fraction nothing in a run can raise a point off room temperature.
        for (double work : new double[] {0.0, 1.0e6, 1.0e9}) {
            assertEquals(law.roomTemperature(), law.temperature(work, 8960.0), 0.0,
                    "plastic work " + work + " must not heat anything");
        }
        assertEquals(JohnsonCook.COPPER_OFHC.strainTerm(0.3), law.strainTerm(0.3), 0.0,
                "the strain term is the one thing that must survive untouched");
    }

    @Test
    @DisplayName("the pressure curve peaks where the instability condition says it does")
    void curvePeaksAtTheInstability() {
        JohnsonCook law = JohnsonCook.COPPER_OFHC.quasiStatic();
        double peak = Burst.burstHoopStrain(law);
        double best = Burst.burstPressure(law, 2.0e-3, 25.0e-3);

        // Brute force, independent of the derivative condition that produced it.
        for (int i = 1; i < 4000; i++) {
            double strain = i * 1.0e-4;
            assertTrue(Burst.pressure(law, 2.0e-3, 25.0e-3, strain) <= best,
                    "found a higher pressure at " + strain + " than at the claimed peak "
                            + peak);
        }
    }
}
