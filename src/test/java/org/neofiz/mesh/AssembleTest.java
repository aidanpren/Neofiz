package org.neofiz.mesh;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.neofiz.core.Formulation;
import org.neofiz.core.Material;
import org.neofiz.solver.Contact;
import org.neofiz.solver.ExplicitSolver;
import org.neofiz.solver.Integration;
import org.neofiz.solver.Kinematics;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Welded against stacked.
 *
 * <p>The same two rectangles, meshed two ways, are two different physical objects, and the
 * distinction is the whole difference between a clad plate and a wall of bricks. Sharing nodes
 * is a weld that cannot come apart; separate nodes are two bodies whose only relationship is
 * whatever contact gives them. Getting this wrong does not throw -- it quietly simulates a
 * monolith and calls it a stack.
 */
class AssembleTest {

    private static final double MM = 1.0e-3;

    private static final Outline LOWER = Outline.rectangle(0.0, 0.0, 6.0 * MM, 6.0 * MM);
    private static final Outline UPPER = Outline.rectangle(0.0, 6.0 * MM, 6.0 * MM, 12.0 * MM);

    @Test
    @DisplayName("touching shapes weld when meshed together and do not when assembled")
    void nodeCountsSayWhichItIs() {
        QuadMesh welded = Outline.mesh(MM, LOWER, UPPER);
        QuadMesh stacked = Outline.assemble(MM, LOWER, UPPER);

        assertEquals(72, welded.elementCount);
        assertEquals(72, stacked.elementCount, "the same cells either way");

        // Welded: 7 columns by 13 rows, the shared row counted once.
        assertEquals(7 * 13, welded.nodeCount);
        // Assembled: two blocks of 7 by 7, the shared row counted twice.
        assertEquals(2 * 7 * 7, stacked.nodeCount);
        assertEquals(welded.nodeCount + 7, stacked.nodeCount);
    }

    @Test
    @DisplayName("a weld carries tension and a stack does not")
    void aWeldCarriesTensionAndAStackDoesNot() {
        // The property, rather than the node count. Pull the upper block straight up: welded,
        // it drags the lower one along; stacked, it simply leaves.
        double weldedPull = lowerBlockRise(Outline.mesh(MM, LOWER, UPPER));
        double stackedPull = lowerBlockRise(Outline.assemble(MM, LOWER, UPPER));

        assertTrue(weldedPull > 1.0e-6,
                "a welded lower block should be dragged upward, rose " + weldedPull);
        assertTrue(Math.abs(stackedPull) < 1.0e-9,
                "a stacked lower block should be left behind, moved " + stackedPull);
    }

    /** How far the lower block's base rises when the upper block is pulled away. */
    private static double lowerBlockRise(QuadMesh mesh) {
        ExplicitSolver s = new ExplicitSolver(mesh, Material.STEEL_4340,
                Formulation.PLANE_STRAIN, Integration.REDUCED, Kinematics.FINITE_STRAIN,
                10.0 * MM, 0.5);
        s.setContact(new Contact(mesh));
        for (int i = 0; i < mesh.nodeCount; i++) {
            if (mesh.z[i] > 6.0 * MM) s.setVelocity(i, 0.0, 40.0);
        }
        s.run(1500);

        double rise = 0.0;
        double[] uz = s.axialDisplacement();
        for (int i = 0; i < mesh.nodeCount; i++) {
            if (mesh.z[i] == 0.0) rise = Math.max(rise, uz[i]);
        }
        return rise;
    }

    @Test
    @DisplayName("an assembled stack is one mesh with several free surfaces")
    void surfacesAreCountedPerBody() {
        // Twelve separate blocks have twelve perimeters, and contact has to see all of them.
        // A traversal-based surface finder would have found whichever one it started on.
        Outline[] blocks = new Outline[6];
        for (int i = 0; i < blocks.length; i++) {
            blocks[i] = Outline.rectangle(0.0, i * 5.0 * MM, 4.0 * MM, (i + 1) * 5.0 * MM);
        }
        QuadMesh mesh = Outline.assemble(MM, blocks);
        assertEquals(6 * 2 * (4 + 5), mesh.surface().edgeCount());
        assertEquals(6 * 4 * 5, mesh.elementCount);
    }
}
