package org.neofiz.solver;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.neofiz.core.Formulation;
import org.neofiz.core.JohnsonCook;
import org.neofiz.core.Material;
import org.neofiz.mesh.Outline;
import org.neofiz.mesh.QuadMesh;

import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Deleting elements that have failed, so that a body can come apart.
 *
 * <p>Erosion is the one thing in this solver that is frankly not physics: continuum mechanics
 * does not let material cease to exist. What makes it defensible is the division of labour --
 * the <em>softening</em> that gets an element to zero strength is regularised by fracture
 * energy and is mesh-independent, and deletion is only the disposal of something that can no
 * longer be integrated. So what these assert is not that the answer is right, which erosion
 * cannot promise, but that the disposal is clean: mass stays, momentum does not move, the
 * energy that went missing is reported rather than dropped, and the surfaces that were hidden
 * behind the deleted element become real surfaces that other bodies can hit.
 */
class ErosionTest {

    private static final double MM = 1.0e-3;
    private static final double DEPTH = 10.0 * MM;
    private static final Material COPPER =
            Material.COPPER_OFHC.withJohnsonCook(JohnsonCook.COPPER_OFHC);

    /** A bar pulled apart from both ends until the middle fails. */
    private static ExplicitSolver stretching(QuadMesh mesh, double failureStrain) {
        ExplicitSolver s = new ExplicitSolver(mesh, COPPER, Formulation.PLANE_STRAIN,
                Integration.REDUCED, Kinematics.FINITE_STRAIN, DEPTH, 0.5);
        double[] eps = new double[mesh.elementCount];
        double[] gf = new double[mesh.elementCount];
        Arrays.fill(eps, failureStrain);
        Arrays.fill(gf, 6.0e3);
        s.setDamage(gf, eps);
        for (int i = 0; i < mesh.nodeCount; i++) {
            s.setVelocity(i, mesh.r[i] < 10.0 * MM ? -80.0 : 80.0, 0.0);
        }
        return s;
    }

    private static QuadMesh bar() {
        return Outline.rectangle(0.0, 0.0, 20.0 * MM, 4.0 * MM).mesh(MM);
    }

    @Test
    @DisplayName("nothing erodes unless erosion is switched on")
    void offByDefault() {
        // Every validation gate in this project runs without it, and has to keep running the
        // arithmetic it always did -- which is why the flag is a null array rather than an
        // all-false one.
        QuadMesh mesh = bar();
        ExplicitSolver s = stretching(mesh, 0.05);
        s.run(4000);

        assertEquals(0, s.erodedElements());
        assertTrue(s.maxDamage() > 0.99, "the bar should have failed somewhere by now");
        for (int e = 0; e < mesh.elementCount; e++) {
            assertTrue(!s.isEroded(e), "element " + e + " vanished with erosion off");
        }
    }

    @Test
    @DisplayName("a bar pulled apart loses elements and keeps its momentum")
    void barSeparates() {
        // The invariant that makes deletion safe. An element's internal forces sum to zero, so
        // taking them away cannot move the total momentum -- and because the nodes and their
        // mass stay behind, none of the mass goes either. If momentum moved, the deletion
        // would be pushing on something.
        QuadMesh mesh = bar();
        ExplicitSolver s = stretching(mesh, 0.05);
        s.setErosion(1.0);

        double before = s.radialMomentum();
        s.run(4000);

        assertTrue(s.erodedElements() > 0, "nothing failed, so nothing was deleted");
        assertEquals(before, s.radialMomentum(), 1e-9 * Math.max(1e-12, Math.abs(before)) + 1e-14,
                "deleting an element moved the total momentum");
    }

    @Test
    @DisplayName("what was thrown away is reported rather than dropped")
    void erodedEnergyIsAccountedFor() {
        // A deleted element takes its stored elastic energy with it. That is the one term
        // stopping an energy audit closing across a deletion, so it is reported: a run can
        // then say how much of its budget went out with the debris.
        QuadMesh mesh = bar();
        ExplicitSolver s = stretching(mesh, 0.05);
        s.setErosion(1.0);
        s.run(4000);

        assertTrue(s.erodedEnergy() > 0.0,
                "elements were deleted but nothing was booked as lost");
        double budget = s.centredKineticEnergy() + s.strainEnergy() + s.hourglassEnergy()
                + s.plasticDissipation() + s.fractureDissipation() + s.erodedEnergy();
        assertTrue(budget > 0.0);
        // The elastic energy in an element about to fail is small beside the plastic work it
        // took to get there; if this were not so, deletion would be losing the run's energy
        // rather than its rounding.
        assertTrue(s.erodedEnergy() < 0.05 * s.plasticDissipation(),
                "erosion threw away " + s.erodedEnergy() + " J against "
                        + s.plasticDissipation() + " J of plastic work");
    }

    @Test
    @DisplayName("the surface follows, so a hole is a hole")
    void surfaceIsRebuiltAroundTheHole() {
        // The reason Contact has to be told. An element deleted from the middle of a body
        // exposes four faces that were interior a step ago; contact that kept the original
        // surface would let anything entering the hole pass straight through its walls.
        QuadMesh mesh = Outline.rectangle(0.0, 0.0, 5.0 * MM, 5.0 * MM).mesh(MM);
        assertEquals(20, mesh.surface().edgeCount());

        boolean[] dead = new boolean[mesh.elementCount];
        int middle = -1;
        for (int e = 0; e < mesh.elementCount; e++) {
            if (Math.abs(mesh.centroidRadius(e) - 2.5 * MM) < 0.6 * MM
                    && Math.abs(mesh.centroidZ(e) - 2.5 * MM) < 0.6 * MM) {
                middle = e;
            }
        }
        assertTrue(middle >= 0, "no element in the middle of the block");
        dead[middle] = true;

        QuadMesh.Surface after = mesh.surface(dead);
        assertEquals(20 + 4, after.edgeCount(),
                "removing an interior element should add its four faces to the surface");
    }

    @Test
    @DisplayName("an element that turns inside out is removed and counted separately")
    void invertedElementsAreCountedApart() {
        // Not a failure -- a solve that has gone wrong. Keeping an inverted element is worse
        // than losing it, because its Jacobian has changed sign and every force it assembles
        // points the wrong way. Counting it apart from the failures is what lets a run say
        // which of the two happened.
        QuadMesh mesh = bar();
        ExplicitSolver s = stretching(mesh, 0.05);
        s.setErosion(1.0);
        s.run(4000);

        assertTrue(s.invertedElements() <= s.erodedElements());
        assertTrue(s.invertedElements() < 0.5 * mesh.elementCount,
                "half the mesh inverted, which is a result to distrust: "
                        + s.invertedElements());
    }

    @Test
    @DisplayName("erosion without a failure model is refused")
    void refusesWithoutDamage() {
        // There is nothing to erode on. Refusing is better than silently deleting only the
        // elements that happen to invert, which would look like it was working.
        QuadMesh mesh = bar();
        ExplicitSolver s = new ExplicitSolver(mesh, COPPER, Formulation.PLANE_STRAIN,
                Integration.REDUCED, Kinematics.FINITE_STRAIN, DEPTH, 0.5);
        assertThrows(IllegalStateException.class, () -> s.setErosion(1.0));

        ExplicitSolver t = stretching(bar(), 0.05);
        assertThrows(IllegalArgumentException.class, () -> t.setErosion(0.0));
        assertThrows(IllegalArgumentException.class, () -> t.setErosion(1.5));
    }
}
