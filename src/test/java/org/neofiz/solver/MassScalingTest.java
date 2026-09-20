package org.neofiz.solver;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.neofiz.core.Formulation;
import org.neofiz.core.Material;
import org.neofiz.mesh.Outline;
import org.neofiz.mesh.QuadMesh;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Buying timestep with mass.
 *
 * <p>Mass scaling is a lie told on purpose, so what these assert is the shape of the lie: that
 * the things it is supposed to leave alone are left exactly alone, and that the thing it
 * breaks is broken by the factor advertised and no more. Gravity is the case it is for --
 * weight scales with mass, so a scaled body falls at the same g and rests at the same load.
 * Inertia is the case it is not for, and a test that did not show that would be hiding it.
 */
class MassScalingTest {

    private static final double MM = 1.0e-3;
    private static final double G = 9.81;

    private static ExplicitSolver block(double scale) {
        QuadMesh mesh = Outline.rectangle(0.0, 0.0, 8.0 * MM, 8.0 * MM).mesh(MM);
        ExplicitSolver s = new ExplicitSolver(mesh, Material.STEEL_4340,
                Formulation.PLANE_STRAIN, Integration.REDUCED, Kinematics.SMALL_STRAIN,
                10.0 * MM, 0.5);
        if (scale > 1.0) s.setMassScaling(s.timestep() * Math.sqrt(scale));
        return s;
    }

    @Test
    @DisplayName("a hundredfold mass buys a tenfold step")
    void stepGrowsAsTheSquareRoot() {
        // The whole arithmetic of it: inflating mass by f divides the wave speed by sqrt(f)
        // and so multiplies the stable step by sqrt(f). Asserting the factor rather than the
        // step keeps this independent of the mesh.
        ExplicitSolver plain = block(1.0);
        ExplicitSolver heavy = block(100.0);

        assertEquals(100.0, heavy.massScale(), 1e-9);
        assertEquals(10.0 * plain.timestep(), heavy.timestep(), 1e-12 * plain.timestep(),
                "a hundred times the mass should be ten times the step");
    }

    @Test
    @DisplayName("gravity is untouched, which is the case it exists for")
    void fallingIsUnaffected() {
        // Weight scales with mass, so the acceleration does not. A scaled body falls at g,
        // lands at the same speed, and rests on the floor at the same load -- which is why
        // settling, stacking and toppling are what this is licensed for.
        ExplicitSolver plain = block(1.0);
        ExplicitSolver heavy = block(400.0);
        plain.setGravity(0.0, -G);
        heavy.setGravity(0.0, -G);

        final double until = 200.0e-6;
        while (plain.time() < until) plain.run(1);
        while (heavy.time() < until) heavy.run(1);

        double plainDrop = lowest(plain);
        double heavyDrop = lowest(heavy);
        assertTrue(Math.abs(heavyDrop / plainDrop - 1.0) < 0.01,
                "a scaled body should fall the same distance: " + plainDrop + " vs " + heavyDrop);
        assertTrue(plainDrop < 0.0, "nothing fell");
    }

    @Test
    @DisplayName("momentum is wrong by exactly the factor, which is the price")
    void inertiaIsWrongByTheFactor() {
        // Stated as a test rather than only as a warning in a comment. A body at a given
        // speed carries f times the momentum, so anything that depends on inertia -- an
        // impact, a ricochet, a transfer -- is wrong by f. Knowing it is wrong by exactly f
        // is what makes it usable.
        ExplicitSolver plain = block(1.0);
        ExplicitSolver heavy = block(25.0);
        setAll(plain, -50.0);
        setAll(heavy, -50.0);

        assertEquals(25.0, heavy.axialMomentum() / plain.axialMomentum(), 1e-9,
                "momentum should be off by exactly the mass factor");
        assertEquals(25.0, heavy.kineticEnergy() / plain.kineticEnergy(), 1e-9,
                "and so should the kinetic energy");
    }

    @Test
    @DisplayName("a target shorter than the current step is refused")
    void willNotRemoveMass() {
        // The failure this stops is a caller passing a step they wanted rather than one they
        // can have, and silently getting nothing. Mass scaling only ever adds.
        ExplicitSolver s = block(1.0);
        double now = s.timestep();
        assertThrows(IllegalArgumentException.class, () -> s.setMassScaling(now * 0.5));
        assertThrows(IllegalArgumentException.class, () -> s.setMassScaling(0.0));
        assertEquals(1.0, s.massScale(), 0.0, "a refused call must change nothing");
        assertEquals(now, s.timestep(), 0.0);
    }

    @Test
    @DisplayName("the contact stiffness follows the step, or the run diverges")
    void contactStaysStable() {
        // The trap. Contact stiffness is built from the reference step, so a step that grows
        // tenfold without the stiffness following makes the penalty frequency ten times what
        // the integrator can carry. This ran away before setMassScaling moved dt0 with it.
        QuadMesh mesh = Outline.assemble(MM,
                Outline.rectangle(0.0, 0.0, 8.0 * MM, 4.0 * MM),
                Outline.rectangle(0.0, 4.0 * MM, 8.0 * MM, 8.0 * MM));
        ExplicitSolver s = new ExplicitSolver(mesh, Material.STEEL_4340,
                Formulation.PLANE_STRAIN, Integration.REDUCED, Kinematics.FINITE_STRAIN,
                10.0 * MM, 0.5);
        s.setMassScaling(s.timestep() * 10.0);
        s.setContact(new Contact(mesh));
        s.setGravity(0.0, -G);
        s.setRigidWall(RigidWall.atZ(0.0));
        s.run(4000);

        assertTrue(Double.isFinite(s.kineticEnergy()), "the run diverged");
        assertTrue(s.kineticEnergy() < 1.0e-3,
                "two blocks resting under gravity should be nearly still, not "
                        + s.kineticEnergy() + " J");
        assertTrue(lowest(s) > -MM, "the stack sank through the floor");
    }

    private static void setAll(ExplicitSolver s, double vz) {
        for (int i = 0; i < s.axialDisplacement().length; i++) s.setVelocity(i, 0.0, vz);
    }

    private static double lowest(ExplicitSolver s) {
        double low = 0.0;
        for (double u : s.axialDisplacement()) low = Math.min(low, u);
        return low;
    }
}
