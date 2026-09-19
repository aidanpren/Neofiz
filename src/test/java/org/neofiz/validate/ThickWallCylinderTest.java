package org.neofiz.validate;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.EnumSource;
import org.junit.jupiter.params.provider.MethodSource;
import org.neofiz.core.Formulation;
import org.neofiz.core.Material;
import org.neofiz.mesh.QuadMesh;
import org.neofiz.solver.ExplicitSolver;
import org.neofiz.solver.Integration;
import org.neofiz.validate.ThickWallCylinderCase.EndCondition;

import java.util.Arrays;
import java.util.stream.Stream;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Tier 1 validation. These have closed-form answers, so a deviation is a bug rather than a
 * modelling choice, and they run on every build.
 *
 * <p>Validation is a gate, not an aspiration, and not something done in a final phase.
 * Late validation is how simulation projects die: discovering at M6 that the constitutive
 * models cannot hit their targets wastes a year.
 */
class ThickWallCylinderTest {

    private static final double A = 0.025;
    private static final double B = 0.030;
    private static final double P = 50.0e6;

    @Test
    @DisplayName("Lame satisfies its own boundary conditions")
    void lameBoundaryConditions() {
        assertEquals(-P, Lame.radialStress(A, A, B, P), P * 1e-12,
                "radial stress at the bore must equal minus the applied pressure");
        assertEquals(0.0, Lame.radialStress(B, A, B, P), P * 1e-12,
                "radial stress at a free outer surface must vanish");
        assertTrue(Lame.hoopStress(A, A, B, P) > Lame.hoopStress(B, A, B, P),
                "hoop stress must peak at the bore, which is where a vessel fails");
    }

    @Test
    @DisplayName("the two end conditions are distinguishable only by displacement")
    void endConditionsShareStressesButNotDisplacements() {
        double strain = Lame.radialDisplacementPlaneStrain(A, A, B, P, 205e9, 0.29);
        double stress = Lame.radialDisplacementPlaneStress(A, A, B, P, 205e9, 0.29);
        assertTrue(Math.abs(stress - strain) / strain > 0.05,
                "the closed forms must differ enough for the displacement check to have teeth");
    }

    /** Every end condition against every quadrature rule. Four independent gate runs. */
    static Stream<Arguments> cases() {
        return Arrays.stream(EndCondition.values()).flatMap(end ->
                Arrays.stream(Integration.values()).map(rule -> Arguments.of(end, rule)));
    }

    @ParameterizedTest(name = "{0} / {1}")
    @MethodSource("cases")
    @DisplayName("GATE 1: axisymmetric FEM matches Lame to 1% at production resolution")
    void matchesLameAtProductionResolution(EndCondition end, Integration rule) {
        var r = ThickWallCylinderCase.run(Material.STEEL_4340, A, B, P, 10, end, rule);

        // A genuinely settled run on this problem reaches KE/SE around 1e-20, so 1e-9 still
        // leaves eleven orders of headroom. The looser 1e-3 this started at was the reason a
        // mistuned relaxation schedule once passed the audit while carrying 1.4% error.
        assertTrue(r.kineticOverStrainEnergy() < 1e-9,
                "solution has not relaxed to a static state; KE/SE = " + r.kineticOverStrainEnergy()
                        + ", so comparison against a static closed form is meaningless");

        assertTrue(r.maxDisplacementErrorPercent() < 1.0,
                end + " / " + rule + ": max radial displacement error "
                        + r.maxDisplacementErrorPercent() + "% exceeds 1%");
        assertTrue(r.maxHoopStressErrorPercent() < 1.0,
                end + " / " + rule + ": max hoop stress error "
                        + r.maxHoopStressErrorPercent() + "% exceeds 1%");
    }

    @ParameterizedTest
    @EnumSource(EndCondition.class)
    @DisplayName("reduced integration agrees with full integration, not merely with Lame")
    void reducedTracksFullIntegration(EndCondition end) {
        var full = ThickWallCylinderCase.run(Material.STEEL_4340, A, B, P, 10, end, Integration.FULL);
        var reduced = ThickWallCylinderCase.run(Material.STEEL_4340, A, B, P, 10, end, Integration.REDUCED);

        double relative = Math.abs(reduced.boreDisplacementFem() - full.boreDisplacementFem())
                / Math.abs(full.boreDisplacementFem());
        assertTrue(relative < 1e-3,
                end + ": the two quadrature rules disagree by " + (100 * relative)
                        + "% at the bore. Both may still pass a 1% gate against Lame while "
                        + "differing from each other, so this is the tighter statement");
    }

    @ParameterizedTest
    @EnumSource(EndCondition.class)
    @DisplayName("hourglass energy is a rounding error, not load-bearing")
    void hourglassEnergyStaysNegligible(EndCondition end) {
        var r = ThickWallCylinderCase.run(Material.STEEL_4340, A, B, P, 10, end, Integration.REDUCED);

        assertTrue(r.hourglassOverStrainEnergy() < 1e-4,
                end + ": hourglass energy is " + (100 * r.hourglassOverStrainEnergy())
                        + "% of strain energy. The answer is being propped up by a numerical "
                        + "stiffness rather than by the material");
    }

    @Test
    @DisplayName("reduced integration does not volumetrically lock; full integration does")
    void reducedIntegrationCuresVolumetricLocking() {
        // nu = 0.499 is a stand-in for something this project walks into at M1 by a
        // different route: plastic flow is incompressible, so the tangent response of a
        // fully yielded element approaches exactly this limit. An element that locks here
        // will lock in developed plasticity, where there is no closed form to catch it.
        Material nearlyIncompressible = new Material("nu -> 0.5", 205e9, 0.499, 7850.0);

        var full = ThickWallCylinderCase.run(nearlyIncompressible, A, B, P, 3,
                EndCondition.RESTRAINED, Integration.FULL);
        var reduced = ThickWallCylinderCase.run(nearlyIncompressible, A, B, P, 3,
                EndCondition.RESTRAINED, Integration.REDUCED);

        assertTrue(full.kineticOverStrainEnergy() < 1e-9 && reduced.kineticOverStrainEnergy() < 1e-9,
                "both runs must be settled or this compares noise, not locking");

        assertTrue(full.maxDisplacementErrorPercent() > 5.0,
                "full integration is expected to lock on this coarse mesh; it reported only "
                        + full.maxDisplacementErrorPercent() + "%. If that is now small, "
                        + "something changed and the justification for reduced integration "
                        + "needs rechecking rather than this bound being relaxed");

        assertTrue(reduced.maxDisplacementErrorPercent() < 1.0,
                "reduced integration must stay accurate as nu approaches 0.5; got "
                        + reduced.maxDisplacementErrorPercent() + "%");
    }

    @Test
    @DisplayName("the answer does not depend on the hourglass coefficient")
    void answerIsInsensitiveToTheHourglassCoefficient() {
        // If the result moves when this dial moves, the dial is doing physics, and the value
        // chosen for it silently becomes a material property.
        double[] boreDisplacements = new double[3];
        double[] qs = {0.01, 0.05, 0.10};
        for (int i = 0; i < qs.length; i++) {
            boreDisplacements[i] = boreWithHourglassCoefficient(qs[i]);
        }
        for (int i = 1; i < qs.length; i++) {
            double relative = Math.abs(boreDisplacements[i] - boreDisplacements[0])
                    / Math.abs(boreDisplacements[0]);
            assertTrue(relative < 1e-4,
                    "bore displacement moved " + (100 * relative) + "% between hourglass "
                            + "coefficients " + qs[0] + " and " + qs[i]);
        }
    }

    private static double boreWithHourglassCoefficient(double q) {
        int nr = 10, nz = 4;
        double h = (B - A) / nr;
        QuadMesh mesh = QuadMesh.cylinderWall(A, B, h * nz, nr, nz);
        ExplicitSolver s = new ExplicitSolver(mesh, Material.STEEL_4340,
                Formulation.AXISYMMETRIC, Integration.REDUCED, 1.0, 0.5);
        s.setHourglassCoefficient(q);
        s.fixAllAxial();

        double rMid = 0.5 * (A + B);
        double omega = Material.STEEL_4340.barWaveSpeed() / rMid;
        double period = 2.0 * Math.PI / omega;
        s.setPressureRamp(P, 8.0 * period);
        s.setRelaxationDamping(0.5, omega);
        s.run((int) Math.ceil(13.0 * period / s.timestep()));
        return s.radialDisplacement()[0];
    }

    @Test
    @DisplayName("refining the mesh reduces the error")
    void convergesUnderRefinement() {
        var coarse = ThickWallCylinderCase.run(Material.STEEL_4340, A, B, P, 2, EndCondition.RESTRAINED);
        var fine = ThickWallCylinderCase.run(Material.STEEL_4340, A, B, P, 16, EndCondition.RESTRAINED);

        assertTrue(fine.maxHoopStressErrorPercent() < coarse.maxHoopStressErrorPercent(),
                "hoop stress error must fall under refinement: coarse "
                        + coarse.maxHoopStressErrorPercent() + "% vs fine "
                        + fine.maxHoopStressErrorPercent() + "%");
    }

    @Test
    @DisplayName("the CFL timestep scales with element size")
    void timestepFollowsCfl() {
        var coarse = ThickWallCylinderCase.run(Material.STEEL_4340, A, B, P, 5, EndCondition.RESTRAINED);
        var fine = ThickWallCylinderCase.run(Material.STEEL_4340, A, B, P, 10, EndCondition.RESTRAINED);

        assertEquals(2.0, coarse.timestep() / fine.timestep(), 1e-9,
                "halving the element size must halve the stable timestep; no mass scaling is allowed");
    }
}
