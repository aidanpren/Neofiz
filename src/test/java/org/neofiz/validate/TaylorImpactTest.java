package org.neofiz.validate;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.EnumSource;
import org.neofiz.core.Formulation;
import org.neofiz.core.JohnsonCook;
import org.neofiz.core.Material;
import org.neofiz.mesh.QuadMesh;
import org.neofiz.solver.ExplicitSolver;
import org.neofiz.solver.Integration;
import org.neofiz.solver.Kinematics;
import org.neofiz.solver.RigidWall;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The M1 gate: Taylor cylinder impact.
 *
 * <p>There is no closed form here, which changes what a test can assert. The audits are exact
 * and are asserted as such. The answer is compared against experiment, and the shape of that
 * comparison is the substantive content: the two outputs of this case behave completely
 * differently, and the tests are arranged around that rather than around a single tolerance.
 *
 * <p><b>Final length is a constitutive measure</b> -- converged, rule-independent, and almost
 * untouched by the interface condition. <b>Mushroom diameter is an interface measure</b>,
 * bracketed between free sliding and fully stuck by about a millimetre, with the measurement
 * inside that bracket. So the gate is asserted as containment rather than as a distance from
 * the sliding limit, because the sliding limit is a boundary condition this case does not
 * have evidence for.
 */
class TaylorImpactTest {

    /** Cheap enough to run in a test, coarse enough to be honest about it. */
    private static final int COARSE = 4;
    /** The production configuration below: reduced integration, square elements. */
    private static final int MEDIUM = 8;

    // ------------------------------------------------------------------ the mesh

    @Test
    @DisplayName("the solid cylinder mesh puts nodes on the axis and carries no pressure edges")
    void solidCylinderGeometry() {
        final int nr = 4, nz = 6;
        QuadMesh mesh = QuadMesh.solidCylinder(2.0e-3, 6.0e-3, nr, nz);

        assertEquals((nr + 1) * (nz + 1), mesh.nodeCount);
        assertEquals(nr * nz, mesh.elementCount);
        assertEquals(0, mesh.pressureEdges.length,
                "the r = 0 column is a symmetry line, not a surface; pressure on it is "
                        + "meaningless and an edge there would be silently mis-weighted");

        for (int n : QuadMesh.axisNodes(nr, nz)) {
            assertEquals(0.0, mesh.r[n], 0.0, "axis node " + n + " is not at r = 0");
        }
        for (int n : QuadMesh.nearFaceNodes(nr, nz)) {
            assertEquals(0.0, mesh.z[n], 0.0, "near-face node " + n + " is not at z = 0");
        }
        for (int n : QuadMesh.farFaceNodes(nr, nz)) {
            assertEquals(6.0e-3, mesh.z[n], 1e-15, "far-face node " + n + " is not at the end");
        }

        // Elements touching the axis are the ones that could divide by zero. They do not,
        // because the hoop strain is only ever formed at a quadrature point, and every
        // quadrature point of an element with non-zero radial extent lies at r > 0.
        for (int e = 0; e < mesh.elementCount; e++) {
            assertTrue(mesh.signedArea(e) > 0.0, "element " + e + " is inverted or degenerate");
        }
    }

    @Test
    @DisplayName("material on the axis stays on the axis, exactly")
    void theAxisHolds() {
        final int nr = 3;
        final int nz = 8;
        QuadMesh mesh = QuadMesh.solidCylinder(0.5 * TaylorImpactCase.DIAMETER, 8.0e-3, nr, nz);
        ExplicitSolver s = new ExplicitSolver(mesh, TaylorImpactCase.johnsonCook4340(),
                Formulation.AXISYMMETRIC, Integration.REDUCED, Kinematics.FINITE_STRAIN, 1.0, 0.5);
        for (int n : QuadMesh.axisNodes(nr, nz)) s.fixRadial(n);
        s.setRigidWall(RigidWall.atZ(0.0));
        s.setUniformVelocity(0.0, -TaylorImpactCase.IMPACT_SPEED);
        while (s.time() < 20e-6) s.step();

        // u_r = 0 on the centreline is a symmetry condition, not a boundary condition, and
        // omitting it does not crash -- it quietly lets the axis open into a hole down the
        // middle of the specimen while every other number stays plausible.
        for (int n : QuadMesh.axisNodes(nr, nz)) {
            assertEquals(0.0, s.radialDisplacement()[n], 0.0,
                    "axis node " + n + " moved off the centreline");
        }
        assertTrue(s.maxPlasticStrain() > 0.05, "the specimen should have deformed by now");
    }

    @Test
    @DisplayName("an elastic specimen is refused rather than measured")
    void elasticIsRefused() {
        // Elastic material bounces off with no permanent set, so "final length" and "mushroom
        // diameter" are not measurements of anything. The run would succeed and report
        // numbers, which is the failure mode worth refusing.
        IllegalArgumentException e = assertThrows(IllegalArgumentException.class,
                () -> TaylorImpactCase.run(Material.STEEL_4340, TaylorImpactCase.LENGTH,
                        TaylorImpactCase.DIAMETER, TaylorImpactCase.IMPACT_SPEED,
                        COARSE, Integration.REDUCED, TaylorImpactCase.FRICTIONLESS, 10e-6));
        assertTrue(e.getMessage().contains("yield"), "got: " + e.getMessage());
    }

    // ------------------------------------------------------------------ the audits

    @ParameterizedTest
    @EnumSource(Integration.class)
    @DisplayName("axial momentum change equals the anvil impulse, to floating point")
    void momentumBalanceIsExact(Integration integration) {
        var r = TaylorImpactCase.reference(COARSE, integration);
        assertTrue(r.momentumBalanceError() < 1e-11,
                "momentum balance broken by " + r.momentumBalanceError() + " relative");
    }

    @ParameterizedTest
    @EnumSource(Integration.class)
    @DisplayName("the energy balance closes: everything that left the kinetic account arrived somewhere")
    void energyBalanceCloses(Integration integration) {
        var r = TaylorImpactCase.reference(COARSE, integration);

        // Plastic work, residual kinetic and elastic energy, the arrival loss and the
        // hourglass springs. Nothing else is allowed to absorb the impact.
        assertTrue(Math.abs(r.energyBalanceError()) < 0.01,
                "energy balance is off by " + r.energyBalanceError() * 100
                        + " %; the impact energy is going somewhere unaccounted for");

        // Almost all of it is plastic work. Worth asserting, because it is what makes the
        // other fractions interpretable: if this were not dominant, the balance above would
        // be an identity between small residuals and would close whether or not the run was
        // sound.
        assertTrue(r.plasticFraction() > 0.95,
                "only " + r.plasticFraction() * 100 + " % of the impact energy became plastic "
                        + "work; for a Taylor test that should be around 98 %");
    }

    @Test
    @DisplayName("hourglass control carries a rounding error of the impact energy, not a share of it")
    void hourglassIsNegligible() {
        var r = TaylorImpactCase.reference(MEDIUM, Integration.REDUCED);

        // Against the energy that went in -- not against recoverable strain energy, which is
        // the right denominator for the quasi-static cylinder gates and the wrong one here.
        // By the time an impacted specimen has unloaded, strain energy is a small residual
        // rather than the scale of the problem, and hourglass-over-strain-energy reads about
        // 1.0 on this very run while hourglass control is in fact carrying 0.2 %.
        assertTrue(r.hourglassFraction() < 1e-2,
                "hourglass springs hold " + r.hourglassFraction() * 100 + " % of the impact "
                        + "energy; the mushroom is being shaped by the hourglass coefficient "
                        + "rather than by the material");

        assertEquals(0.0, TaylorImpactCase.reference(COARSE, Integration.FULL).hourglassFraction(),
                0.0, "full integration has no hourglass term to hold anything");
    }

    // ------------------------------------------------------------------ Johnson-Cook

    @Test
    @DisplayName("the specimen heats by a couple of hundred kelvin, from its own plastic work")
    void temperatureRiseIsPhysical() {
        var r = TaylorImpactCase.reference(MEDIUM, Integration.REDUCED);
        double rise = r.maxTemperature() - JohnsonCook.STEEL_4340.roomTemperature();

        // Not a free parameter: the rise is fixed by the plastic work the run actually did,
        // divided by rho c_p. Two hundred kelvin is what 4340 steel does when it is deformed
        // to unit plastic strain adiabatically, and it is enough to matter -- at this
        // temperature the thermal term has taken 14 % off the flow stress.
        assertTrue(rise > 100.0 && rise < 400.0,
                "temperature rose by " + rise + " K, which is outside anything this much "
                        + "plastic work can produce in steel");
        assertTrue(JohnsonCook.STEEL_4340.thermalTerm(r.maxTemperature()) < 0.92,
                "at " + r.maxTemperature() + " K the thermal term should have taken a visible "
                        + "bite out of the flow stress");
    }

    @Test
    @DisplayName("Johnson-Cook stiffens the response: the specimen shortens less")
    void johnsonCookStiffensTheResponse() {
        var jc = TaylorImpactCase.reference(MEDIUM, Integration.REDUCED);
        var surrogate = TaylorImpactCase.referenceSurrogate(MEDIUM, Integration.REDUCED);

        // Length is the clean read on the flow stress, and it moves by two and a half points
        // of the original length -- far outside the mesh convergence of either. This is the
        // measurement that shows Johnson-Cook did something.
        assertTrue(jc.lengthRatio() > surrogate.lengthRatio() + 0.015,
                "Johnson-Cook gave L/L0 = " + jc.lengthRatio() + " against the surrogate's "
                        + surrogate.lengthRatio() + "; a stiffer material must shorten less");
        assertTrue(jc.maxPlasticStrain() < surrogate.maxPlasticStrain(),
                "a stiffer material must reach less plastic strain");
    }

    @Test
    @DisplayName("Johnson-Cook barely moves the mushroom, because its two corrections cancel")
    void johnsonCookBarelyMovesTheMushroom() {
        var jc = TaylorImpactCase.reference(MEDIUM, Integration.REDUCED);
        var surrogate = TaylorImpactCase.referenceSurrogate(MEDIUM, Integration.REDUCED);

        // The finding that reframed the M1 gate. Rate hardening was expected to be worth
        // about 16 % on flow stress and to narrow the mushroom accordingly; thermal softening
        // turns out to give nearly all of it back, so the net effect on the mushroom is under
        // a per cent. See JohnsonCookTest#theTwoTermsNearlyCancel for the two factors
        // measured separately.
        //
        // If this assertion ever fails because the gap widened, the conclusion drawn from it
        // -- that the residual over-prediction is an interface effect and not a constitutive
        // one -- has to be re-derived rather than assumed.
        double change = Math.abs(jc.mushroomDiameter() - surrogate.mushroomDiameter())
                / surrogate.mushroomDiameter();
        assertTrue(change < 0.02,
                "Johnson-Cook moved the mushroom by " + change * 100 + " %, from "
                        + surrogate.mushroomDiameter() * 1e3 + " to "
                        + jc.mushroomDiameter() * 1e3 + " mm");
    }

    // ------------------------------------------------------------------ the gate

    @Test
    @DisplayName("GATE M1: the measured mushroom lies inside the interface bracket")
    void theMeasurementLiesInsideTheInterfaceBracket() {
        var sliding = TaylorImpactCase.reference(MEDIUM, Integration.REDUCED);
        var stuck = TaylorImpactCase.referenceWelded(MEDIUM, Integration.REDUCED);

        // This is the gate, and it is stated as containment rather than as a distance from
        // the frictionless answer on purpose. The two interface limits are a millimetre apart
        // -- an order of magnitude more than Johnson-Cook was worth -- so quoting the sliding
        // number against a 5 % tolerance would be reporting the precision of a boundary
        // condition this case has no evidence for. Coulomb friction with a real coefficient
        // is what turns this from a bracket into a number.
        assertTrue(stuck.mushroomDiameter() < TaylorImpactCase.MEASURED_MUSHROOM_DIAMETER,
                "the fully stuck limit came out at " + stuck.mushroomDiameter() * 1e3
                        + " mm, above the measured 9.5 mm. If even a fully gripped face "
                        + "spreads this far, the bracket no longer contains the answer and "
                        + "something other than friction is wrong");
        assertTrue(sliding.mushroomDiameter() > TaylorImpactCase.MEASURED_MUSHROOM_DIAMETER,
                "the frictionless limit came out at " + sliding.mushroomDiameter() * 1e3
                        + " mm, below the measured 9.5 mm");

        // And the bracket is wide, which is the point.
        double width = (sliding.mushroomDiameter() - stuck.mushroomDiameter())
                / TaylorImpactCase.MEASURED_MUSHROOM_DIAMETER;
        assertTrue(width > 0.05,
                "the interface bracket is only " + width * 100 + " % wide; if free sliding and "
                        + "a fully gripped face really agree this closely, friction is not the "
                        + "explanation for the residual and the gate can be asserted directly");
    }

    @Test
    @DisplayName("the mushroom crosses the measurement at a physically plausible friction")
    void frictionSweepCrossesTheMeasurement() {
        // The bracket above says the answer is reachable; this says it is reachable at a
        // coefficient a real steel-on-steel interface could have, which is the difference
        // between a model that spans the measurement and one that explains it.
        var loose = TaylorImpactCase.referenceAt(MEDIUM, Integration.REDUCED, 0.02);
        var tight = TaylorImpactCase.referenceAt(MEDIUM, Integration.REDUCED, 0.10);

        assertTrue(loose.mushroomDiameter() > TaylorImpactCase.MEASURED_MUSHROOM_DIAMETER,
                "mu = 0.02 already narrows the mushroom past the measurement");
        assertTrue(tight.mushroomDiameter() < TaylorImpactCase.MEASURED_MUSHROOM_DIAMETER,
                "mu = 0.10 does not reach the measurement");

        // Monotone in between, which is what makes the crossing a crossing rather than a
        // coincidence of two samples.
        double previous = Double.MAX_VALUE;
        for (double mu : new double[]{0.0, 0.02, 0.05, 0.10, 0.15}) {
            double d = TaylorImpactCase.referenceAt(MEDIUM, Integration.REDUCED, mu)
                    .mushroomDiameter();
            assertTrue(d < previous,
                    "the mushroom did not narrow between the previous coefficient and " + mu);
            previous = d;
        }
    }

    @Test
    @DisplayName("the coefficient that matches is not mesh-converged, and is not a measurement of friction")
    void theMatchingCoefficientAbsorbsMeshError() {
        // The caveat that stops the sweep above from being over-read. The frictionless
        // baseline still drifts upward with refinement, so the coefficient needed to pull it
        // back onto 9.5 mm drifts with it -- about 0.05 at this mesh and 0.07 at the next.
        // Quoting a fitted mu as "the" friction coefficient would be reporting a
        // discretisation error as a material property.
        double atMedium = TaylorImpactCase.referenceAt(MEDIUM, Integration.REDUCED, 0.05)
                .mushroomErrorVsMeasuredPercent();
        double atCoarse = TaylorImpactCase.referenceAt(COARSE, Integration.REDUCED, 0.05)
                .mushroomErrorVsMeasuredPercent();

        assertTrue(Math.abs(atMedium - atCoarse) > 1.0,
                "the same coefficient gave " + atCoarse + " % and " + atMedium + " % on two "
                        + "meshes. If those have converged, the fitted coefficient has become "
                        + "meaningful and this caveat can be dropped");
    }

    @Test
    @DisplayName("friction dissipates energy, and the balance still closes with it")
    void frictionalDissipationIsAccountedFor() {
        var r = TaylorImpactCase.referenceAt(MEDIUM, Integration.REDUCED, 0.15);

        assertTrue(r.frictionDissipationFraction() > 0.005,
                "a gripped interface at this speed should absorb a percent or two of the "
                        + "impact energy; got " + r.frictionDissipationFraction());

        // Sliding friction is physics the interface is supposed to dissipate, so it enters the
        // balance as a term of its own. Folding it into the arrival loss would make a
        // frictional anvil indistinguishable from a leaky constraint.
        assertTrue(Math.abs(r.energyBalanceError()) < 0.01,
                "with friction on, the energy balance is off by " + r.energyBalanceError() * 100
                        + " %; the dissipation term is missing or double-counted");
        assertTrue(r.momentumBalanceError() < 1e-11,
                "friction disturbed the axial momentum balance");
    }

    @Test
    @DisplayName("the mushroom converges more slowly than the length, and both still converge")
    void theTwoOutputsConvergeDifferently() {
        // A live guard on the expensive study pinned in MeshConvergenceTest. Run at 4 / 8 / 16
        // rather than 8 / 16 / 32 so the suite stays quick, which means it is below the
        // asymptotic range -- the observed order here comes out near 0.7 rather than the 1.09
        // the finer triple gives. That is fine for what this asserts: the sequence still has
        // to be converging, and the mushroom still has to be slower than the element's own
        // second order.
        MeshConvergence mushroom =
                TaylorImpactCase.mushroomConvergence(Integration.REDUCED, 0.0, COARSE);
        MeshConvergence length =
                TaylorImpactCase.lengthConvergence(Integration.REDUCED, 0.0, COARSE);

        assertTrue(mushroom.isConverging(),
                "the mushroom sequence stopped converging: " + mushroom.coarse() * 1e3 + ", "
                        + mushroom.medium() * 1e3 + ", " + mushroom.fine() * 1e3);
        assertTrue(mushroom.observedOrder() < 1.6,
                "the mushroom converged at order " + mushroom.observedOrder() + "; if it has "
                        + "reached the element's second order then the contact-edge "
                        + "singularity argument no longer holds and the extrapolation in the "
                        + "write-up should be replaced by straightforward refinement");

        // Length is the contrast, and it is already so converged that Richardson has nothing
        // to work with -- which the tool says out loud rather than extrapolating anyway.
        assertFalse(length.isAsymptotic(),
                "the length ratio is expected to be converged, not extrapolable");
        assertTrue(length.lastChange() < 1e-3,
                "length moved " + length.lastChange() + " on the last refinement");
    }

    @Test
    @DisplayName("the frictionless limit alone drifts past 5 % under refinement")
    void theSlidingLimitDriftsUnderRefinement() {
        double coarse = TaylorImpactCase.reference(COARSE, Integration.REDUCED)
                .mushroomErrorVsMeasuredPercent();
        double medium = TaylorImpactCase.reference(MEDIUM, Integration.REDUCED)
                .mushroomErrorVsMeasuredPercent();

        // A characterisation test, not an aspiration. Taken alone the frictionless mushroom
        // is inside 5 % at these meshes and drifts out of it by the next one, which is
        // exactly why the gate above is stated as containment. Recording the direction means
        // that when friction lands, the change can be attributed: a real friction coefficient
        // narrows the mushroom at fixed mesh, and this drift stays.
        assertTrue(medium > coarse,
                "the frictionless mushroom no longer grows under refinement (" + coarse
                        + " % then " + medium + " %). If that is now converged, this test has "
                        + "outlived its purpose");
    }

    @Test
    @DisplayName("final length converges, agrees between rules, and ignores the interface")
    void finalLengthIsAConstitutiveMeasure() {
        var coarse = TaylorImpactCase.reference(COARSE, Integration.REDUCED);
        var medium = TaylorImpactCase.reference(MEDIUM, Integration.REDUCED);
        var full = TaylorImpactCase.reference(MEDIUM, Integration.FULL);
        var stuck = TaylorImpactCase.referenceWelded(MEDIUM, Integration.REDUCED);

        // Everything the mushroom is not. Three independent things are varied here -- the
        // mesh, the quadrature rule and the interface condition -- and the length ratio moves
        // by less than half a percent against any of them. That is what makes it the output
        // this case can actually pin the flow stress with.
        assertEquals(coarse.lengthRatio(), medium.lengthRatio(), 0.005, "mesh");
        assertEquals(medium.lengthRatio(), full.lengthRatio(), 0.005, "quadrature rule");
        assertEquals(medium.lengthRatio(), stuck.lengthRatio(), 0.005,
                "the interface condition moves the mushroom by a millimetre and the length by "
                        + "almost nothing; if that stops being true the two outputs are no "
                        + "longer measuring different things");

        assertTrue(medium.lengthRatio() > 0.88 && medium.lengthRatio() < 0.93,
                "length ratio " + medium.lengthRatio() + " is outside anything plausible for "
                        + "4340 steel at 181 m/s");
    }
}
