package org.neofiz.solver;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.neofiz.core.Formulation;
import org.neofiz.core.Material;
import org.neofiz.mesh.Outline;
import org.neofiz.mesh.QuadMesh;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * A uniform body force.
 *
 * <p>Trivial to add and easy to get subtly wrong, which is why it is worth three assertions.
 * A body force applied to anything other than the lumped nodal mass gives the wrong
 * acceleration on a graded mesh; applied to the wrong sign it reads as a very convincing
 * simulation of the ceiling; and applied inside the parallel assembly it would need a source
 * slot and a share of the gather. It is none of those.
 */
class GravityTest {

    private static final double MM = 1.0e-3;
    private static final double G = 9.81;

    private static ExplicitSolver falling(QuadMesh mesh) {
        ExplicitSolver s = new ExplicitSolver(mesh, Material.STEEL_4340,
                Formulation.PLANE_STRAIN, Integration.REDUCED, Kinematics.SMALL_STRAIN,
                10.0 * MM, 0.5);
        s.setGravity(0.0, -G);
        return s;
    }

    @Test
    @DisplayName("a free body falls at g, whatever it is made of or meshed at")
    void freeFallIsExact() {
        // The one case with a closed form, and it has to hold to round-off: a body with
        // nothing touching it has no internal forces at all, so every node sees exactly its
        // own weight over its own mass.
        QuadMesh mesh = Outline.rectangle(0.0, 0.0, 6.0 * MM, 4.0 * MM).mesh(MM);
        ExplicitSolver s = falling(mesh);
        s.run(2000);

        double t = s.time();
        double[] uz = s.axialDisplacement();
        // The discrete answer, not the continuum one. Leapfrog from rest gives
        // u_n = a*dt^2*n(n+1)/2, which is 0.5*a*t*(t + dt) -- half a step of lead over the
        // textbook 0.5*a*t^2, and at 2000 steps that is a tenth of a per cent. Asserting the
        // continuum form to round-off would have failed for the right reason and been fixed
        // by loosening the tolerance, which would have thrown away the exactness this has.
        double expected = -0.5 * G * t * (t + s.timestep());
        for (int i = 0; i < mesh.nodeCount; i++) {
            assertEquals(expected, uz[i], 1e-12 * Math.abs(expected),
                    "node " + i + " did not fall at g");
        }
        assertEquals(-G * t, s.axialMomentum() / totalMass(s, mesh), 1e-9 * G * t,
                "momentum should be the whole weight integrated over the run");
    }

    @Test
    @DisplayName("a body resting on the anvil stays there")
    void restingBodyDoesNotSink() {
        // Gravity plus the kinematic wall, which is the combination every sandbox scene uses
        // and the one where a sign error is invisible in the first frame.
        QuadMesh mesh = Outline.rectangle(0.0, 0.0, 20.0 * MM, 5.0 * MM).mesh(MM);
        ExplicitSolver s = falling(mesh);
        s.setRigidWall(RigidWall.atZ(0.0));
        s.run(4000);

        double lowest = Double.MAX_VALUE;
        double[] uz = s.axialDisplacement();
        for (int i = 0; i < mesh.nodeCount; i++) lowest = Math.min(lowest, mesh.z[i] + uz[i]);
        assertTrue(lowest >= -1e-15, "the body sank through the anvil to " + lowest);
        // And it is held rather than bouncing: a slab of steel under one g compresses by
        // nanometres, so nothing should have moved perceptibly.
        assertTrue(Math.abs(s.kineticEnergy()) < 1e-6,
                "a resting slab should be still, kinetic energy " + s.kineticEnergy());
    }

    @Test
    @DisplayName("off by default, because every gate here would rather it were")
    void defaultIsNone() {
        // A tube bursting at 20 MPa carries six orders of magnitude more than its own weight.
        // Switching gravity on by default would move every validation number in this project
        // in its eighth decimal place and make each one depend on which way the mesh is drawn.
        QuadMesh mesh = Outline.rectangle(0.0, 0.0, 6.0 * MM, 4.0 * MM).mesh(MM);
        ExplicitSolver s = new ExplicitSolver(mesh, Material.STEEL_4340,
                Formulation.PLANE_STRAIN, Integration.REDUCED, Kinematics.SMALL_STRAIN,
                10.0 * MM, 0.5);
        s.run(500);
        assertEquals(0.0, s.axialMomentum(), 0.0, "something pulled on an unforced body");
    }

    private static double totalMass(ExplicitSolver s, QuadMesh mesh) {
        return Material.STEEL_4340.density() * meshArea(mesh) * 10.0 * MM;
    }

    private static double meshArea(QuadMesh mesh) {
        double a = 0.0;
        for (int e = 0; e < mesh.elementCount; e++) a += mesh.signedArea(e);
        return a;
    }
}
