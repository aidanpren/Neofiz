package org.neofiz.solver;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.EnumSource;
import org.neofiz.core.Formulation;
import org.neofiz.core.Material;
import org.neofiz.mesh.QuadMesh;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Element-level tests for reduced integration and Flanagan-Belytschko hourglass control.
 *
 * <p>These drive single elements through one prescribed deformation increment rather than
 * running a simulation, which is what lets them assert exact values instead of tolerances.
 * Hourglass control has two properties worth having, and they pull against each other: it
 * must be stiff enough to suppress a mode the quadrature cannot see, and it must be
 * completely invisible to every deformation the element is supposed to represent. The second
 * is the one that is easy to get wrong and hard to notice, because a hourglass term that
 * leaks into real strain shows up only as a quiet loss of accuracy.
 */
class HourglassTest {

    private static final Material STEEL = Material.STEEL_4340;

    /** A deliberately skewed quad. On a rectangle the hourglass projection is trivial. */
    private static QuadMesh distortedQuad() {
        return new QuadMesh(
                new double[]{1.00, 2.10, 1.90, 0.85},
                new double[]{0.00, 0.15, 1.20, 0.95},
                new int[]{0, 1, 2, 3},
                new int[0]);
    }

    private static QuadMesh unitSquare() {
        return new QuadMesh(
                new double[]{1.0, 2.0, 2.0, 1.0},
                new double[]{0.0, 0.0, 1.0, 1.0},
                new int[]{0, 1, 2, 3},
                new int[0]);
    }

    private static ExplicitSolver solver(QuadMesh mesh, Formulation f, Integration i, double q) {
        ExplicitSolver s = new ExplicitSolver(mesh, STEEL, f, i, 1.0, 0.5);
        s.setHourglassCoefficient(q);
        return s;
    }

    /**
     * Imposes a displacement field through one increment of the real solve and returns the
     * resulting nodal forces.
     *
     * <p>The solver is incremental now: stress is integrated along a strain path rather than
     * evaluated from total displacement. So a prescribed displacement is expressed as the
     * velocity that produces it in one timestep. Driving the production path is the point --
     * a test that pokes state in directly stops testing the code that ships.
     */
    private static double[] impose(ExplicitSolver s, double[] uR, double[] uZ) {
        for (int i = 0; i < 4; i++) s.setVelocity(i, uR[i] / s.timestep(), uZ[i] / s.timestep());
        s.advance();
        double[] out = new double[8];
        for (int i = 0; i < 4; i++) {
            out[2 * i] = s.forceR()[i];
            out[2 * i + 1] = s.forceZ()[i];
        }
        return out;
    }

    /** The alternating hourglass pattern in the radial direction, amplitude e. */
    private static double[] hourglassField(double e) {
        return new double[]{e, -e, e, -e};
    }

    private static double[] zeros() {
        return new double[4];
    }

    private static double norm(double[] v) {
        double n = 0.0;
        for (double x : v) n += x * x;
        return Math.sqrt(n);
    }

    @ParameterizedTest
    @EnumSource(value = Formulation.class, names = {"AXISYMMETRIC", "PLANE_STRAIN"})
    @DisplayName("hourglass control is invisible to any linear displacement field")
    void linearFieldsProduceNoHourglassForce(Formulation formulation) {
        QuadMesh mesh = distortedQuad();

        // An arbitrary linear field: translation, rotation and constant strain all at once.
        // gamma is constructed orthogonal to every one of these, so the hourglass term must
        // contribute exactly nothing -- on a skewed element, not just a rectangle.
        double[] uR = new double[4];
        double[] uZ = new double[4];
        for (int i = 0; i < 4; i++) {
            uR[i] = 1.0e-3 + 2.0e-3 * mesh.r[i] - 0.7e-3 * mesh.z[i];
            uZ[i] = -0.4e-3 + 1.1e-3 * mesh.r[i] + 1.5e-3 * mesh.z[i];
        }

        double[] without = impose(solver(mesh, formulation, Integration.REDUCED, 0.0), uR, uZ);
        ExplicitSolver on = solver(mesh, formulation, Integration.REDUCED, 0.05);
        double[] with = impose(on, uR, uZ);

        double scale = norm(without);
        assertTrue(scale > 0.0, "the linear field must actually load the element");

        for (int i = 0; i < 8; i++) {
            assertEquals(without[i], with[i], scale * 1e-12,
                    "hourglass control perturbed the response to a linear field at dof " + i
                            + "; the shape vector is not orthogonal to constant strain and the "
                            + "correction is corrupting real deformation");
        }

        assertEquals(0.0, on.hourglassEnergy(), on.strainEnergy() * 1e-12,
                "a linear field must store no hourglass energy");
    }

    @Test
    @DisplayName("without control the hourglass mode is a genuine zero-energy mode")
    void theModeIsFreeWithoutControl() {
        QuadMesh mesh = unitSquare();
        final double e = 1.0e-4;

        double[] cheap = impose(
                solver(mesh, Formulation.PLANE_STRAIN, Integration.REDUCED, 0.0),
                hourglassField(e), zeros());
        double[] expensive = impose(
                new ExplicitSolver(mesh, STEEL, Formulation.PLANE_STRAIN, Integration.FULL, 1.0, 0.5),
                hourglassField(e), zeros());

        // The reference scale: what a comparable amount of real strain costs on this element.
        double[] stretch = {e * mesh.r[0], e * mesh.r[1], e * mesh.r[2], e * mesh.r[3]};
        double scale = norm(impose(
                solver(mesh, Formulation.PLANE_STRAIN, Integration.REDUCED, 0.0), stretch, zeros()));

        assertTrue(norm(cheap) < scale * 1e-12,
                "one-point integration must see zero strain in the hourglass mode; got "
                        + norm(cheap) + " against a reference of " + scale
                        + ". If this is nonzero the mode is not actually free and the rest of "
                        + "this test is measuring something else");

        assertTrue(norm(expensive) > scale * 0.1,
                "full integration must resist the same mode strongly -- that contrast is the "
                        + "whole reason hourglass control exists");
    }

    @Test
    @DisplayName("control resists the mode, linearly and in the right direction")
    void controlOpposesTheMode() {
        QuadMesh mesh = unitSquare();
        final double e = 1.0e-4;

        ExplicitSolver s = solver(mesh, Formulation.PLANE_STRAIN, Integration.REDUCED, 0.05);
        double[] f = impose(s, hourglassField(e), zeros());
        double energy = s.hourglassEnergy();

        // On a rectangle gamma is exactly the alternating pattern, so the restoring force
        // must be antiparallel to it.
        double projection = 0.0;
        for (int i = 0; i < 4; i++) projection += ((i % 2 == 0) ? 1.0 : -1.0) * f[2 * i];
        assertTrue(projection < 0.0,
                "hourglass force must oppose the mode, not drive it; projection " + projection);

        assertTrue(energy > 0.0, "the restoring spring must store energy");

        // Stiffness form: force linear in displacement, energy quadratic.
        ExplicitSolver twice = solver(mesh, Formulation.PLANE_STRAIN, Integration.REDUCED, 0.05);
        double[] f2 = impose(twice, hourglassField(2 * e), zeros());

        assertEquals(2.0, norm(f2) / norm(f), 1e-9,
                "doubling the hourglass amplitude must double the force");
        assertEquals(4.0, twice.hourglassEnergy() / energy, 1e-9,
                "hourglass energy must be quadratic in amplitude; if it is not, the reported "
                        + "energy does not match the force actually applied and the audit lies");
    }

    @Test
    @DisplayName("the hourglass coefficient scales the stiffness it claims to scale")
    void coefficientScalesStiffness() {
        QuadMesh mesh = distortedQuad();
        final double e = 1.0e-4;

        double[] previous = null;
        for (double q : new double[]{0.02, 0.04, 0.08}) {
            double[] with = impose(
                    solver(mesh, Formulation.PLANE_STRAIN, Integration.REDUCED, q),
                    hourglassField(e), zeros());
            double[] without = impose(
                    solver(mesh, Formulation.PLANE_STRAIN, Integration.REDUCED, 0.0),
                    hourglassField(e), zeros());

            double[] delta = new double[8];
            for (int i = 0; i < 8; i++) delta[i] = with[i] - without[i];

            if (previous != null) {
                assertEquals(2.0, norm(delta) / norm(previous), 1e-9,
                        "doubling the coefficient must double the hourglass force");
            }
            previous = delta;
        }
    }
}
