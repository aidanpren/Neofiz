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
 * Penalty contact between deformable surfaces.
 *
 * <p>There is no closed form for two meshed blocks hitting one another, so what this leans on
 * instead are the invariants that hold whatever the collision does. Momentum is the strong one
 * and it is exact: every pair applies equal and opposite forces, so the sum over the mesh
 * cannot move at all, and a normal pointing the wrong way, a reaction on the wrong node or a
 * weight that does not sum to one all break it immediately.
 *
 * <p>Energy is the weak one, and the thresholds here say so. Explicit penalty contact closes
 * the budget to a few per cent, and the amount is set by the stiffness rather than by the mesh
 * -- see {@link Contact} for the measured table. The audit is also read with
 * {@link ExplicitSolver#centredKineticEnergy()} rather than the half-step kinetic energy,
 * because a collision turns translation into vibration and the half-step form is biased for
 * anything that vibrates. {@link EnergyAuditTest} pins that down separately, and it exists
 * because this gate nearly blamed the bias on the contact model.
 */
class ContactTest {

    private static final double MM = 1.0e-3;
    private static final double CELL = 1.0 * MM;
    private static final double DEPTH = 10.0 * MM;
    /** Elastic, so that anything the energy audit sees is contact rather than plasticity. */
    private static final Material STEEL = Material.STEEL_4340;

    /** Two blocks with a gap between them, the left one twice as tall as the right. */
    private static QuadMesh twoBlocks() {
        return Outline.mesh(CELL,
                Outline.rectangle(0.0, 0.0, 8.0 * MM, 12.0 * MM),
                Outline.rectangle(13.0 * MM, 2.0 * MM, 21.0 * MM, 8.0 * MM));
    }

    private static ExplicitSolver solver(QuadMesh mesh) {
        return new ExplicitSolver(mesh, STEEL, Formulation.PLANE_STRAIN,
                Integration.REDUCED, Kinematics.FINITE_STRAIN, DEPTH, 0.5);
    }

    /** Drives everything left of the gap right at {@code v}, everything right of it left. */
    private static void approach(ExplicitSolver s, QuadMesh mesh, double v) {
        for (int i = 0; i < mesh.nodeCount; i++) {
            s.setVelocity(i, mesh.r[i] < 10.0 * MM ? v : -v, 0.0);
        }
    }

    private static double total(ExplicitSolver s) {
        return s.centredKineticEnergy() + s.strainEnergy() + s.hourglassEnergy()
                + s.contact().energy() + s.contact().dissipation();
    }

    @Test
    @DisplayName("two blocks in one mesh pass through each other without contact")
    void withoutContactTheyInterpenetrate() {
        // The state of affairs contact exists to end, asserted rather than assumed. Nothing in
        // the force assembly couples separate element groups, so before this feature the two
        // blocks were simulated correctly and independently, and met by overlapping.
        QuadMesh mesh = twoBlocks();
        ExplicitSolver s = solver(mesh);
        approach(s, mesh, 60.0);
        s.run(3000);

        double[] ur = s.radialDisplacement();
        double leftFront = -Double.MAX_VALUE, rightBack = Double.MAX_VALUE;
        for (int i = 0; i < mesh.nodeCount; i++) {
            if (mesh.r[i] == 8.0 * MM) leftFront = Math.max(leftFront, mesh.r[i] + ur[i]);
            if (mesh.r[i] == 13.0 * MM) rightBack = Math.min(rightBack, mesh.r[i] + ur[i]);
        }
        assertTrue(leftFront > rightBack,
                "with no contact the faces should have crossed: " + leftFront + " vs " + rightBack);
    }

    @Test
    @DisplayName("a collision conserves momentum to round-off")
    void momentumIsConserved() {
        // The strongest statement available about contact, and the one that fails loudly for
        // almost any mistake. Asymmetric blocks, so the total is not zero by construction and
        // a sign error cannot hide inside it.
        QuadMesh mesh = twoBlocks();
        ExplicitSolver s = solver(mesh);
        s.setContact(new Contact(mesh));
        approach(s, mesh, 60.0);

        double before = s.radialMomentum();
        s.run(4000);

        assertEquals(before, s.radialMomentum(), 1e-9 * Math.abs(before),
                "radial momentum moved during a collision");
        assertEquals(0.0, s.axialMomentum(), 1e-12,
                "a head-on collision must not produce axial momentum");
    }

    @Test
    @DisplayName("the blocks bounce, nothing tunnels, and the energy is accounted for")
    void elasticCollisionClosesTheEnergyBudget() {
        QuadMesh mesh = twoBlocks();
        ExplicitSolver s = solver(mesh);
        s.setContact(new Contact(mesh));
        approach(s, mesh, 60.0);

        double initial = total(s);
        boolean touched = false;
        for (int k = 0; k < 400; k++) {
            s.run(20);
            if (s.contact().pairCount() > 0) touched = true;
        }

        assertTrue(touched, "the blocks never came into contact");
        assertEquals(0, s.contact().escapedNodes(),
                "a node it was resolving went too deep to resolve; see Contact");
        // Four per cent, not a threshold chosen to pass: the measured closure at the default
        // stiffness is 2.8 %, and the shape of that number is in Contact's class note.
        assertTrue(Math.abs(total(s) / initial - 1.0) < 0.04,
                "energy closed to " + (100.0 * (total(s) / initial - 1.0)) + " %");

        // And they separated: the left block is no longer moving right.
        double leftV = 0.0, leftM = 0.0;
        for (int i = 0; i < mesh.nodeCount; i++) {
            if (mesh.r[i] < 10.0 * MM) {
                leftV += s.velocityR()[i];
                leftM++;
            }
        }
        assertTrue(leftV / leftM < 0.0,
                "the left block should have been turned around, mean v = " + leftV / leftM);
    }

    @Test
    @DisplayName("penetration stays a small fraction of a cell")
    void penetrationIsBounded() {
        // What the mass-based stiffness is chosen for, and the reason the contact geometry
        // means anything: a surface sunk a quarter of a cell into another is not a surface.
        QuadMesh mesh = twoBlocks();
        ExplicitSolver s = solver(mesh);
        s.setContact(new Contact(mesh));
        approach(s, mesh, 60.0);

        double worst = 0.0;
        for (int k = 0; k < 400; k++) {
            s.run(20);
            worst = Math.max(worst, s.contact().maxPenetration());
        }
        assertTrue(worst > 0.0, "nothing ever touched");
        assertTrue(worst < 0.10 * CELL,
                "penetration reached " + (100.0 * worst / CELL) + " % of a cell");
    }

    @Test
    @DisplayName("a surface at rest does not contact itself")
    void noSelfContactAtRest() {
        // A closed surface is always exactly touching itself, so without the exclusions every
        // segment would register against its neighbours and the body would leap apart on the
        // first step. An L rather than a rectangle, because the re-entrant corner is where one
        // ring of topological exclusion is not enough.
        QuadMesh mesh = Outline.of(0.0, 0.0, 20.0 * MM, 0.0, 20.0 * MM, 6.0 * MM,
                6.0 * MM, 6.0 * MM, 6.0 * MM, 20.0 * MM, 0.0, 20.0 * MM).mesh(CELL);
        ExplicitSolver s = solver(mesh);
        s.setContact(new Contact(mesh));
        s.run(200);

        assertEquals(0, s.contact().pairCount(), "a body at rest touched itself");
        assertEquals(0, s.contact().escapedNodes(),
                "a body at rest reported nodes tunnelling through it");
        assertEquals(0.0, s.kineticEnergy(), 1e-18, "a body at rest started moving");
    }

    @Test
    @DisplayName("a node is not caught by a surface it is not facing")
    void facingIsRequired() {
        // The bug this exists to stop, stated as a property. The region behind a segment is a
        // half-plane, so every node of a block is arbitrarily deep behind the slabs of its own
        // perpendicular faces, at projections landing squarely inside real segments. A block
        // squeezed hard enough to bring those faces within the depth cap would then contact
        // itself -- and, because the master is the deepest candidate, would prefer the
        // imaginary pair to any real one.
        QuadMesh mesh = Outline.rectangle(0.0, 0.0, 3.0 * MM, 3.0 * MM).mesh(CELL);
        ExplicitSolver s = solver(mesh);
        s.setContact(new Contact(mesh));
        for (int i = 0; i < mesh.nodeCount; i++) {
            if (mesh.r[i] == 0.0) s.setVelocity(i, 120.0, 0.0);
            if (mesh.r[i] == 3.0 * MM) s.setVelocity(i, -120.0, 0.0);
        }
        s.run(1500);

        assertEquals(0, s.contact().pairCount(),
                "a block crushed against itself found an imaginary contact");
    }

    @Test
    @DisplayName("friction opposes sliding and is paid for out of the energy budget")
    void frictionOpposesSliding() {
        QuadMesh mesh = Outline.mesh(CELL,
                Outline.rectangle(0.0, 0.0, 20.0 * MM, 5.0 * MM),
                Outline.rectangle(5.0 * MM, 6.0 * MM, 13.0 * MM, 11.0 * MM));

        double free = slideSpeed(mesh, 0.0);
        double rough = slideSpeed(mesh, 0.4);
        assertTrue(rough < free, "friction did not slow the slider: " + rough + " vs " + free);
        assertTrue(rough > 0.0, "friction reversed the slider, which a cone cannot do");
    }

    /** Mean sideways speed of the upper block after it is pressed down onto the lower one. */
    private static double slideSpeed(QuadMesh mesh, double mu) {
        ExplicitSolver s = solver(mesh);
        Contact c = new Contact(mesh).withFriction(mu);
        s.setContact(c);
        for (int i = 0; i < mesh.nodeCount; i++) {
            if (mesh.z[i] > 5.5 * MM) s.setVelocity(i, 30.0, -40.0);
            if (mesh.z[i] == 0.0) s.fixAxial(i);
        }

        double initial = total(s);
        for (int q = 0; q < 30; q++) {
            s.run(200);
            System.out.printf("SLIDE mu=%.1f t=%5d pairs=%3d ke=%.3f se=%.3f cn=%.3f dis=%.3f tot=%.4f pen=%.3e esc=%d%n",
                mu, (q + 1) * 200, c.pairCount(), s.centredKineticEnergy(), s.strainEnergy(),
                c.energy(), c.dissipation(), total(s) / initial, c.maxPenetration(), c.escapedNodes());
        }

        assertTrue(c.dissipation() >= 0.0, "friction produced energy");
        assertTrue(Math.abs(total(s) / initial - 1.0) < 0.10,
                "the sliding run closed to " + (100.0 * (total(s) / initial - 1.0)) + " %");

        double v = 0.0, n = 0.0;
        for (int i = 0; i < mesh.nodeCount; i++) {
            if (mesh.z[i] > 5.5 * MM) {
                v += s.velocityR()[i];
                n++;
            }
        }
        return v / n;
    }

    @Test
    @DisplayName("a stiffness that could break the timestep is refused")
    void refusesUnstableStiffness() {
        // The bound is exact rather than advisory: omega*dt = sqrt(scale), and central
        // difference is stable to 2. What it does not say is that the model is *useful* all
        // the way there; Contact's table shows it is not.
        QuadMesh mesh = twoBlocks();
        assertThrows(IllegalArgumentException.class,
                () -> new Contact(mesh).withStiffnessScale(4.0));
        assertThrows(IllegalArgumentException.class,
                () -> new Contact(mesh).withStiffnessScale(0.0));
        assertThrows(IllegalArgumentException.class,
                () -> new Contact(mesh).withFriction(-0.1));
    }

    @Test
    @DisplayName("the surface is the free surface, and both blocks are on it")
    void surfaceCoversEveryFreeEdge() {
        // 8x12 and 8x6 blocks: two perimeters, 40 and 28 edges. Counting them by hand is the
        // check that a mesh of several disconnected pieces is handled as one surface rather
        // than as whichever piece a traversal happened to start on.
        QuadMesh mesh = twoBlocks();
        assertEquals(2 * (8 + 12) + 2 * (8 + 6), new Contact(mesh).segmentCount());
    }
}
