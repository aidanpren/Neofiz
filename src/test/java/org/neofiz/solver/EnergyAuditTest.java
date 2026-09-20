package org.neofiz.solver;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.neofiz.core.Formulation;
import org.neofiz.core.Material;
import org.neofiz.mesh.Outline;
import org.neofiz.mesh.QuadMesh;

import java.util.function.ToDoubleFunction;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * What a total-energy audit can and cannot say.
 *
 * <p>Every gate in this project that closes an energy budget adds kinetic energy to stored
 * energy, and central difference does not hold those two at the same instant: velocity is
 * carried at half steps and displacement at whole ones. For a body in rigid translation that
 * costs nothing, because the kinetic energy is not changing. For a body <b>vibrating</b> it is
 * a real bias of order {@code omega*dt/2} of the vibrational energy -- <em>first</em> order in
 * the step, not second, because the leading term is {@code (dt/2) * d(KE)/dt}.
 *
 * <p>This file exists because the contact work ran into exactly that and nearly blamed it on
 * the contact model. A collision turns translation into vibration, so the audit reads clean
 * before the impact and biased after it, and the step change looks precisely like energy
 * created by whatever was added last. Measuring the bias on a bar that touches nothing is what
 * separates the two.
 *
 * <p>The correction is exact and costs nothing:
 * {@link ExplicitSolver#centredKineticEnergy()}. Measured here it takes the audit on a
 * vibrating bar from 1.36 % to 0.16 %, and turns first-order convergence into second.
 */
class EnergyAuditTest {

    private static final double CELL = 1.0e-3;
    private static final double DEPTH = 10.0e-3;

    private static final ToDoubleFunction<ExplicitSolver> HALF = ExplicitSolver::kineticEnergy;
    private static final ToDoubleFunction<ExplicitSolver> CENTRED =
            ExplicitSolver::centredKineticEnergy;

    private static ExplicitSolver vibrating(double cfl) {
        QuadMesh mesh = Outline.rectangle(0.0, 0.0, 20.0 * CELL, 8.0 * CELL).mesh(CELL);
        ExplicitSolver s = new ExplicitSolver(mesh, Material.STEEL_4340,
                Formulation.PLANE_STRAIN, Integration.REDUCED, Kinematics.FINITE_STRAIN,
                DEPTH, cfl);
        // A standing wave along the bar: half a wavelength, so both ends are free. Nothing
        // touches anything, so the total energy is a constant of the continuum problem and
        // every departure the audit shows belongs to the audit.
        for (int i = 0; i < mesh.nodeCount; i++) {
            s.setVelocity(i, 60.0 * Math.sin(Math.PI * mesh.r[i] / (20.0 * CELL)), 0.0);
        }
        return s;
    }

    /** Worst relative departure of the energy budget from its starting value over the run. */
    private static double drift(double cfl, int samples, int per,
                                ToDoubleFunction<ExplicitSolver> kinetic) {
        ExplicitSolver s = vibrating(cfl);
        double initial = kinetic.applyAsDouble(s) + s.strainEnergy() + s.hourglassEnergy();
        double worst = 0.0;
        for (int k = 0; k < samples; k++) {
            s.run(per);
            double total = kinetic.applyAsDouble(s) + s.strainEnergy() + s.hourglassEnergy();
            worst = Math.max(worst, Math.abs(total / initial - 1.0));
        }
        return worst;
    }

    @Test
    @DisplayName("a freely vibrating bar already fails a tight audit, with nothing to touch")
    void vibrationBiasesTheHalfStepAudit() {
        // The measurement that exonerated the contact model. A bar with nothing to hit, given
        // a standing wave and left alone, still moves the budget by over a per cent.
        double bias = drift(0.5, 200, 20, HALF);
        assertTrue(bias > 0.005,
                "a vibrating bar should show the half-step bias; saw " + (100.0 * bias) + " %");
    }

    @Test
    @DisplayName("centring the velocity removes most of it")
    void centringRemovesTheBias() {
        double half = drift(0.5, 200, 20, HALF);
        double centred = drift(0.5, 200, 20, CENTRED);
        assertTrue(centred < half / 5.0,
                "centring should improve the audit fivefold at least: " + half + " -> " + centred);
    }

    @Test
    @DisplayName("the half-step bias is first order in the step and the centred one is second")
    void ordersOfConvergence() {
        // The distinction that says which of the two is the artefact and which is the answer.
        // Quartering the step takes about four off a first-order term and about sixteen off a
        // second-order one; measured, this run gives 4.3 and 10.4. A defect in the solve would
        // not care about the step at all.
        double halfCoarse = drift(0.5, 200, 20, HALF);
        double halfFine = drift(0.125, 200, 80, HALF);
        double centredCoarse = drift(0.5, 200, 20, CENTRED);
        double centredFine = drift(0.125, 200, 80, CENTRED);

        double halfRatio = halfCoarse / halfFine;
        double centredRatio = centredCoarse / centredFine;
        assertTrue(halfRatio > 3.0 && halfRatio < 7.0,
                "the half-step bias should fall about fourfold, saw " + halfRatio);
        assertTrue(centredRatio > 8.0,
                "the centred bias should fall far faster, saw " + centredRatio);
        assertTrue(centredRatio > 2.0 * halfRatio,
                "centred should converge at a visibly higher order: "
                        + centredRatio + " vs " + halfRatio);
    }
}
