package org.neofiz.validate;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.neofiz.core.Formulation;
import org.neofiz.core.JohnsonCook;
import org.neofiz.core.Material;
import org.neofiz.mesh.QuadMesh;
import org.neofiz.solver.ExplicitSolver;
import org.neofiz.solver.Integration;
import org.neofiz.solver.Kinematics;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The M1 burst gate, and the audits that decide how much of it to believe.
 *
 * <p>Every run here is expensive by the standards of this suite -- a burst traverse is a few
 * hundred thousand steps -- so the shared ones are computed once into static fields and the
 * tests interrogate them from different directions. That is not only a speed measure: several
 * of these assertions are about two numbers from the <em>same</em> run having to agree, and
 * that only means something if it is in fact the same run.
 */
class BurstCaseTest {

    /** Reference tube, reference mesh, one ceiling. */
    private static final BurstCase.Result REFERENCE = BurstCase.reference();

    /** The gate figure, with the ceiling extrapolated out. */
    private static final BurstCase.OvershootFit FIT = BurstCase.extrapolated();

    private static final double GATE = 10.0;

    // ------------------------------------------------------------------ the gate

    @Test
    @DisplayName("burst pressure is inside the 10 % exit criterion, by a factor of twenty")
    void gate() {
        assertTrue(REFERENCE.traversedPeak(), "the reference run has to actually burst");
        assertTrue(Math.abs(REFERENCE.burstPressureErrorPercent()) < GATE,
                "M1 exit criterion; got " + REFERENCE.burstPressureErrorPercent() + " %");
        assertTrue(Math.abs(REFERENCE.burstPressureErrorPercent()) < 0.6,
                "and the gate is not where the real accuracy is; got "
                        + REFERENCE.burstPressureErrorPercent() + " %");
        assertTrue(Math.abs(FIT.errorPercent()) < 0.35,
                "with the ceiling extrapolated out; got " + FIT.errorPercent() + " %");
    }

    @Test
    @DisplayName("the ceiling is extrapolated away rather than chosen")
    void overshootExtrapolates() {
        // The perturbation being removed is 0.2 % of the answer. The middle measurement,
        // which is held out of the fit, sits on the line to four parts in ten million. A fit
        // that clean is not evidence of care, it is evidence that the perturbation really is
        // linear in the ceiling with no constant term, which is what the mechanism predicts.
        assertTrue(FIT.linearityPercent() < 1.0e-3,
                "the withheld point is off the line by " + FIT.linearityPercent() + " %");
        double removed = 100.0 * Math.abs(FIT.nearPressure() - FIT.burstPressure())
                / FIT.burstPressure();
        assertTrue(removed > 100.0 * FIT.linearityPercent(),
                "an extrapolation is only worth doing if it removes more than it invents: "
                        + "removed " + removed + " %, residual " + FIT.linearityPercent() + " %");
        assertTrue(FIT.slope() < 0.0,
                "a higher ceiling drives a faster creep and reads a lower peak");
    }

    // ------------------------------------------------------------------ the measurement

    @Test
    @DisplayName("the equilibrium identity holds to eight digits on a settled tube")
    void equilibriumIdentity() {
        // Below burst the tube reaches a genuine static equilibrium, and there P a and the
        // hoop stress resultant are the same number exactly. This is the audit of
        // loadCapacity: it says the integral is right, that it is being evaluated on the
        // configuration the solver assembles forces in, and that the run settled.
        BurstCase.Result settled = BurstCase.overshootSweep(0.90);

        assertFalse(settled.traversedPeak(), "0.90 of burst must not burst");
        assertTrue(settled.kineticOverStrainAtPeak() < 1.0e-12,
                "not settled: KE/SE = " + settled.kineticOverStrainAtPeak());
        assertTrue(settled.equilibriumErrorPercent() < 1.0e-4,
                "P a = integral(sigma_theta dr) failed by "
                        + settled.equilibriumErrorPercent() + " %");
    }

    @Test
    @DisplayName("an unloaded mesh has no capacity and no strain")
    void unloadedIsZero() {
        QuadMesh mesh = QuadMesh.cylinderWall(24.0e-3, 26.0e-3, 0.5e-3, 8, 2);
        ExplicitSolver solver = new ExplicitSolver(mesh, BurstCase.COPPER,
                Formulation.AXISYMMETRIC, Integration.REDUCED, Kinematics.FINITE_STRAIN,
                1.0, 0.5);

        assertEquals(0.0, BurstCase.loadCapacity(solver, mesh), 0.0);
        assertEquals(0.0, BurstCase.hoopStrain(solver, mesh, 25.0e-3), 0.0);
    }

    @Test
    @DisplayName("the traverse stays quasi-static")
    void quasiStatic() {
        assertTrue(REFERENCE.kineticOverStrainAtPeak() < 1.0e-2,
                "inertia is participating at the peak: KE/SE = "
                        + REFERENCE.kineticOverStrainAtPeak());
        assertTrue(REFERENCE.hourglassOverPlasticWork() < 1.0e-6,
                "hourglass energy against plastic work: "
                        + REFERENCE.hourglassOverPlasticWork());
    }

    // ------------------------------------------------------------------ the residual

    @Test
    @DisplayName("the residual against a rigid-plastic oracle is the elasticity it leaves out")
    void residualIsElasticity() {
        // The strong claim: two terms computed from E and nu before the run -- the elastic
        // hoop strain thinning the wall without hardening it, and elastic dilatation working
        // against that thinning -- predict the entire measured deficit. Nothing here is
        // fitted, and the ceiling has been extrapolated out of the measurement rather than
        // chosen, so there is no knob left to have tuned.
        double predicted = REFERENCE.elasticKnockdownPercent();
        assertEquals(predicted, -FIT.errorPercent(), 0.06,
                "the elastic correction predicts a " + predicted + " % deficit and the "
                        + "measurement shows " + (-FIT.errorPercent()) + " %");

        // Keeping only the obvious term would over-predict by a third and leave a residual
        // with nowhere to live. This is why the dilatation is not optional.
        Burst.ElasticState state = Burst.elasticState(BurstCase.COPPER.johnsonCook(),
                REFERENCE.thickness(), REFERENCE.meanRadius(),
                BurstCase.COPPER.youngsModulus(), BurstCase.COPPER.poissonRatio());
        assertTrue(state.dilatation() > 0.2 * 2.0 * state.hoopStrain(),
                "dilatation should take back a substantial share of the thinning; got "
                        + state.dilatation() + " against " + (2.0 * state.hoopStrain()));
        assertTrue(200.0 * state.hoopStrain() - REFERENCE.elasticKnockdownPercent() > 0.05,
                "the two terms should differ by more than the tolerance above, or this test "
                        + "is not distinguishing them");

        // And the weaker, independent corroboration: the peak strain shift measures the same
        // elastic hoop strain without using an elastic constant at all.
        double observed = 100.0 * REFERENCE.observedKnockdown();
        assertTrue(observed > 0.0, "the peak must sit at a higher strain than rigid-plastic, "
                + "because elastic hoop strain thins the wall without hardening it");
        assertEquals(200.0 * state.hoopStrain(), observed, 0.1,
                "the strain shift should see the hoop term: predicted "
                        + (200.0 * state.hoopStrain()) + " %, measured " + observed + " %");
    }

    @Test
    @DisplayName("the thin-wall oracle earns its name: the residual shrinks with the wall")
    void thinWallLimit() {
        // Both of these carry the same ceiling perturbation, which is why the assertion is
        // about the difference between them and not about either one's absolute size. Taking
        // the ceiling out needs three runs per slenderness and belongs in the report, where
        // the remainder runs from +0.21 % at D/t = 10 to -0.03 % at D/t = 80 -- shrinking
        // sixfold, but through a sign change, so the reference wall at D/t = 25 sits close to
        // where it vanishes and flatters the accounting there.
        BurstCase.Result thick = BurstCase.atSlenderness(10.0);
        BurstCase.Result thin = BurstCase.atSlenderness(80.0);

        double thickResidual = Math.abs(thick.burstPressureErrorPercent()
                + thick.elasticKnockdownPercent());
        double thinResidual = Math.abs(thin.burstPressureErrorPercent()
                + thin.elasticKnockdownPercent());

        assertTrue(thinResidual < thickResidual,
                "D/t = 80 residual " + thinResidual + " % is not below D/t = 10's "
                        + thickResidual + " %");
        assertTrue(thickResidual - thinResidual > 0.15,
                "and the improvement should be worth having, not noise: " + thickResidual
                        + " % to " + thinResidual + " %");
    }

    @Test
    @DisplayName("the answer is not mesh-limited")
    void meshIndependent() {
        // Unlike the Taylor mushroom, which is read next to a singularity and converges at
        // first order, this one is an integral of a smooth field and is done at four elements
        // through the wall. Reporting that is as much a result as a convergence rate would be.
        BurstCase.Result coarse = BurstCase.referenceAt(4, Integration.REDUCED);
        double gap = 100.0 * Math.abs(coarse.burstPressureFem() - REFERENCE.burstPressureFem())
                / REFERENCE.burstPressureFem();

        assertTrue(gap < 0.01, "doubling the mesh moved the burst pressure by " + gap + " %");
    }

    @Test
    @DisplayName("a coarse fully integrated mesh gets the right answer for the wrong reason")
    void lockingHidesBehindTheHeadline() {
        // This is the case the second measurement exists for. Full integration locks under
        // the near-incompressible plastic flow, which stiffens the wall and delays the
        // instability. Both errors land on the reported pressure and very nearly cancel, so
        // the headline number comes out better than the converged one -- while the peak
        // strain, which the locking also moves and nothing cancels, is wrong by half again.
        BurstCase.Result locked = BurstCase.referenceAt(4, Integration.FULL);

        assertTrue(Math.abs(locked.burstPressureErrorPercent()) < 0.1,
                "the accident this test is about did not happen; got "
                        + locked.burstPressureErrorPercent() + " %");
        assertTrue(Math.abs(locked.burstPressureErrorPercent())
                        < Math.abs(REFERENCE.burstPressureErrorPercent()),
                "the locked mesh is supposed to look better than the good one");

        Burst.ElasticState state = Burst.elasticState(BurstCase.COPPER.johnsonCook(),
                locked.thickness(), locked.meanRadius(),
                BurstCase.COPPER.youngsModulus(), BurstCase.COPPER.poissonRatio());
        double observed = 100.0 * locked.observedKnockdown();
        double sound = 100.0 * REFERENCE.observedKnockdown();

        assertTrue(observed > 1.4 * sound,
                "and the strain is supposed to give it away: the locked mesh reads "
                        + observed + " % where the sound one reads " + sound + " %");
        assertTrue(observed - 200.0 * state.hoopStrain() > 0.15,
                "the locked mesh should miss the predicted hoop term by far more than the "
                        + "sound one does; got " + observed + " % against "
                        + (200.0 * state.hoopStrain()) + " %");
    }

    // ------------------------------------------------------------------ falsification

    @Test
    @DisplayName("with the geometry frozen there is no burst at all")
    void smallStrainCannotBurst() {
        // The sharpest statement available that the finite-strain path is doing real work.
        // Burst is not a strength, it is a feedback between load and shape, and a solver that
        // refuses to update the shape cannot produce it: this tube sits in settled equilibrium
        // at half again its burst pressure, having quietly expanded by a quarter of its radius.
        BurstCase.Result frozen = BurstCase.smallStrain(1.5);

        assertFalse(frozen.traversedPeak(), "small strain must not find a load maximum");
        assertTrue(frozen.finalHoopStrain() > 0.20,
                "and should reach an absurd strain while doing it; got "
                        + frozen.finalHoopStrain());
        assertTrue(frozen.kineticOverStrainAtPeak() < 1.0e-4,
                "settled, not still climbing: KE/SE = " + frozen.kineticOverStrainAtPeak());

        IllegalStateException e = assertThrows(IllegalStateException.class,
                frozen::requirePeak);
        assertTrue(e.getMessage().contains("small-strain"), e.getMessage());
        assertTrue(e.getMessage().contains("no geometric self-weakening"), e.getMessage());
    }

    @Test
    @DisplayName("a material with no hardening law is refused")
    void refusesNonHardening() {
        Material perfect = Material.COPPER_OFHC.yielding(90.0e6, 0.0);
        IllegalArgumentException e = assertThrows(IllegalArgumentException.class,
                () -> BurstCase.run(perfect, 25.0e-3, 25.0, 8, Integration.REDUCED,
                        Kinematics.FINITE_STRAIN, 1.05));
        assertTrue(e.getMessage().contains("instability to find"), e.getMessage());
    }

    @Test
    @DisplayName("the closed form and the harness agree about which material they are running")
    void sameMaterial() {
        JohnsonCook law = BurstCase.COPPER.johnsonCook();
        assertTrue(law.isQuasiStatic(), "the harness must run the flow curve the oracle assumes");
        assertEquals(Burst.burstHoopStrain(law), REFERENCE.burstHoopStrainExact(), 0.0);
        assertEquals(Burst.burstPressure(law, REFERENCE.thickness(), REFERENCE.meanRadius()),
                REFERENCE.burstPressureExact(), 0.0);
    }
}
