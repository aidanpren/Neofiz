package org.neofiz.validate;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.EnumSource;
import org.neofiz.core.Material;
import org.neofiz.solver.Integration;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The elastic-plastic cylinder gate: J2 plasticity wired into the solver, against a closed
 * form that cannot be satisfied by accident.
 *
 * <p>A thick wall, b/a = 2. Thin-walled cylinders are nearly useless for this: at b/a = 1.2
 * the burst pressure is only 1.19x the elastic limit, so there is almost no contained-plastic
 * regime between "entirely elastic" and "no equilibrium exists" to validate anything in.
 */
class PlasticCylinderTest {

    private static final double A = 0.025;
    private static final double B = 0.050;
    private static final double SY = 800.0e6;
    private static final Material STEEL = Material.STEEL_4340.yielding(SY, 0.0);

    private static final double P_ELASTIC_LIMIT = ElasticPlastic.elasticLimitPressure(A, B, SY);

    @Test
    @DisplayName("the closed form is self-consistent at both ends of its range")
    void closedFormEndpoints() {
        assertEquals(A, ElasticPlastic.frontRadius(P_ELASTIC_LIMIT, A, B, SY), 1e-12,
                "at the elastic limit the front must sit exactly on the bore");
        assertEquals(B, ElasticPlastic.frontRadius(
                        ElasticPlastic.limitPressure(A, B, SY), A, B, SY), 1e-9,
                "at the limit pressure the front must have reached the outer surface");

        // The plastic and elastic branches must agree where they meet.
        double p = 1.4 * P_ELASTIC_LIMIT;
        double c = ElasticPlastic.frontRadius(p, A, B, SY);
        assertEquals(ElasticPlastic.radialStressPlastic(c, A, p, SY),
                ElasticPlastic.radialStressElastic(c, c, B, SY), Math.abs(p) * 1e-9,
                "radial stress must be continuous across the elastic-plastic front");

        assertEquals(-p, ElasticPlastic.radialStress(A, A, B, p, SY), Math.abs(p) * 1e-12,
                "radial stress at the bore must still equal minus the applied pressure");
        assertEquals(0.0, ElasticPlastic.radialStress(B, A, B, p, SY), Math.abs(p) * 1e-9,
                "the outer surface is free");
    }

    @ParameterizedTest
    @EnumSource(Integration.class)
    @DisplayName("GATE M1a: contained plastic flow matches the closed form to 1%")
    void containedPlasticFlowMatchesClosedForm(Integration rule) {
        for (double factor : new double[]{1.10, 1.35, 1.60}) {
            var r = PlasticCylinderCase.run(STEEL, A, B, P_ELASTIC_LIMIT * factor, 20, rule);
            String at = rule + " at p/pe = " + factor + ": ";

            assertTrue(r.kineticOverStrainEnergy() < 1e-9,
                    at + "not settled; KE/SE = " + r.kineticOverStrainEnergy());

            assertTrue(r.maxYieldConditionErrorPercent() < 1.0,
                    at + "yielded elements are off the yield surface by "
                            + r.maxYieldConditionErrorPercent() + "%");
            assertTrue(r.maxRadialStressErrorPercent() < 1.0,
                    at + "radial stress error " + r.maxRadialStressErrorPercent() + "%");
            assertTrue(r.maxHoopStressErrorPercent() < 1.0,
                    at + "hoop stress error " + r.maxHoopStressErrorPercent() + "%");

            assertEquals(r.frontRadiusExact(), r.frontRadiusFem(), 1.01 * r.elementSize(),
                    at + "the plastic front is more than one element from where the closed "
                            + "form puts it");
        }
    }

    @ParameterizedTest
    @EnumSource(Integration.class)
    @DisplayName("below the elastic limit nothing yields at all")
    void nothingYieldsBelowTheElasticLimit(Integration rule) {
        var r = PlasticCylinderCase.run(STEEL, A, B, 0.95 * P_ELASTIC_LIMIT, 20, rule);

        assertEquals(0.0, r.yieldedFraction(), 0.0,
                "the closed form says first yield happens at the elastic limit; anything "
                        + "yielding below it means the yield criterion is miscalibrated, most "
                        + "likely by the 2/sqrt(3) that separates von Mises from Tresca");
        assertEquals(0.0, r.maxPlasticStrain(), 0.0);
        assertTrue(r.maxRadialStressErrorPercent() < 0.2,
                "and the elastic answer must still be right: " + r.maxRadialStressErrorPercent() + "%");
    }

    @Test
    @DisplayName("first yield happens at the bore, and lags the exact limit by exactly O(h)")
    void firstYieldLagsByTheSamplingOffset() {
        // A centroid-sampled element never sees the stress at the bore, only the stress half
        // an element inside it, and the peak is at the bore. So the discrete cylinder yields
        // *late*, by a predictable amount: the elastic hoop-radial difference falls off as
        // 1/r^2, so first yield should land at pe * ((a + h/2)/a)^2.
        //
        // This is worth pinning down rather than tolerating. It is a real O(h) bias in every
        // yield-onset prediction the solver will ever make, it is always in the
        // non-conservative direction, and the only reason it is acceptable is that it
        // converges away. A test that merely asserted "something yields eventually" would
        // hide both the size of the bias and the fact that it has the right scaling.
        double previousGap = Double.MAX_VALUE;

        for (int nr : new int[]{10, 20, 40}) {
            double h = (B - A) / nr;
            double predicted = P_ELASTIC_LIMIT * Math.pow((A + 0.5 * h) / A, 2.0);

            var below = PlasticCylinderCase.run(STEEL, A, B, 0.98 * predicted, nr, Integration.REDUCED);
            var above = PlasticCylinderCase.run(STEEL, A, B, 1.02 * predicted, nr, Integration.REDUCED);

            assertEquals(0.0, below.yieldedFraction(), 0.0,
                    "nr=" + nr + ": yielded below the predicted onset of "
                            + (predicted / 1e6) + " MPa");
            assertTrue(above.yieldedFraction() > 0.0,
                    "nr=" + nr + ": had not yielded above the predicted onset");
            assertTrue(above.frontRadiusFem() <= A + 1.01 * h,
                    "nr=" + nr + ": first yield must appear in the bore element, not elsewhere; "
                            + "front came out at " + (above.frontRadiusFem() * 1000) + " mm");

            double gap = predicted / P_ELASTIC_LIMIT - 1.0;
            assertTrue(gap < previousGap,
                    "the onset bias must shrink under refinement, not merely be small");
            previousGap = gap;
        }
    }

    @Test
    @DisplayName("the plastic zone grows monotonically with pressure")
    void plasticZoneGrowsWithPressure() {
        double previousFront = 0.0;
        double previousBore = 0.0;
        for (double factor : new double[]{1.05, 1.20, 1.40, 1.60, 1.75}) {
            var r = PlasticCylinderCase.run(STEEL, A, B, P_ELASTIC_LIMIT * factor, 20,
                    Integration.REDUCED);
            assertTrue(r.frontRadiusFem() >= previousFront,
                    "the plastic front retreated when the pressure rose, at p/pe = " + factor);
            assertTrue(r.boreDisplacement() > previousBore,
                    "the bore must keep opening as pressure rises, at p/pe = " + factor);
            previousFront = r.frontRadiusFem();
            previousBore = r.boreDisplacement();
        }
    }

    @Test
    @DisplayName("above the limit pressure the case refuses rather than returning a number")
    void refusesWhenNoEquilibriumExists() {
        double beyond = 1.01 * ElasticPlastic.limitPressure(A, B, SY);
        var thrown = assertThrows(IllegalArgumentException.class,
                () -> PlasticCylinderCase.run(STEEL, A, B, beyond, 20, Integration.REDUCED),
                "a fully plastic cylinder has no static answer to relax onto. Returning a "
                        + "stress field anyway is the dangerous outcome: it looks reasonable "
                        + "and the kinetic energy grows too slowly to be obvious");
        assertTrue(thrown.getMessage().contains("limit pressure"));
    }

    @ParameterizedTest
    @EnumSource(Integration.class)
    @DisplayName("refining the mesh reduces the stress error")
    void convergesUnderRefinement(Integration rule) {
        var coarse = PlasticCylinderCase.run(STEEL, A, B, 1.35 * P_ELASTIC_LIMIT, 10, rule);
        var fine = PlasticCylinderCase.run(STEEL, A, B, 1.35 * P_ELASTIC_LIMIT, 40, rule);

        assertTrue(fine.maxRadialStressErrorPercent() < coarse.maxRadialStressErrorPercent(),
                rule + ": radial stress error must fall under refinement; coarse "
                        + coarse.maxRadialStressErrorPercent() + "% vs fine "
                        + fine.maxRadialStressErrorPercent() + "%");
    }

    @Test
    @DisplayName("hardening raises the pressure a given plastic zone needs")
    void hardeningStiffensTheResponse() {
        // Not against a closed form -- there is not a simple one for a hardening cylinder --
        // but the direction is unambiguous, and it checks that the hardening modulus is
        // actually reaching the return map rather than being carried and ignored.
        Material hardening = Material.STEEL_4340.yielding(SY, 8.0e9);
        double p = 1.5 * P_ELASTIC_LIMIT;

        var perfect = PlasticCylinderCase.run(STEEL, A, B, p, 20, Integration.REDUCED);
        var hardened = PlasticCylinderCase.run(hardening, A, B, p, 20, Integration.REDUCED);

        assertTrue(hardened.boreDisplacement() < perfect.boreDisplacement(),
                "a hardening material must open less at the same pressure");
        assertTrue(hardened.maxPlasticStrain() < perfect.maxPlasticStrain(),
                "and must accumulate less plastic strain");
        assertTrue(hardened.frontRadiusFem() <= perfect.frontRadiusFem(),
                "and its plastic zone must not be larger");
    }
}
