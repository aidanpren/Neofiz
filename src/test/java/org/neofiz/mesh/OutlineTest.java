package org.neofiz.mesh;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The polygon mesher.
 *
 * <p>A mesher is the one component whose output is checkable by arithmetic rather than by
 * inspection: a shape of known area, meshed at a known cell size, has a known element count,
 * and every node has to land on the lattice. What cannot be checked that way -- that the
 * elements are wound the right way, that shared corners are welded rather than duplicated --
 * is what {@link QuadMesh} refuses to be built without, so those are asserted here as the
 * properties they are rather than left to the constructor to catch by accident.
 */
class OutlineTest {

    private static final double H = 1.0e-3;

    @Test
    @DisplayName("a rectangle on the lattice meshes to exactly the cells it covers")
    void rectangleIsExact() {
        // 10 mm by 6 mm at a 1 mm cell is sixty elements and seventy-seven nodes, and both
        // numbers can be counted by hand. Anything else means cells were dropped, duplicated,
        // or welded wrongly.
        QuadMesh mesh = Outline.rectangle(0.0, 0.0, 10.0 * H, 6.0 * H).mesh(H);
        assertEquals(60, mesh.elementCount);
        assertEquals(11 * 7, mesh.nodeCount, "shared corners must be welded, not duplicated");
        assertEquals(60.0 * H * H, area(mesh), 1e-15);
    }

    @Test
    @DisplayName("every edge is the cell size, which is what the timestep is bought with")
    void everyEdgeIsTheCellSize() {
        // The whole reason for a lattice rather than a fitted mesher. If this ever stops
        // holding, the stable timestep stops being predictable from the cell size alone.
        QuadMesh mesh = Outline.circle(20.0 * H, 20.0 * H, 8.0 * H, 40).mesh(H);
        assertEquals(H, mesh.minimumEdgeLength(), 1e-15);
    }

    @Test
    @DisplayName("a disc comes out within one cell of its analytic area")
    void discAreaConverges() {
        // Rasterising is an approximation and this is how big it is: at a cell size of a
        // twentieth of the radius the area is within a couple of per cent, and it is the
        // convergence that matters rather than the number.
        double radius = 20.0 * H;
        Outline disc = Outline.circle(0.0, 0.0, radius, 180);
        double coarse = Math.abs(area(disc.mesh(2.0 * H)) / (Math.PI * radius * radius) - 1.0);
        double fine = Math.abs(area(disc.mesh(0.5 * H)) / (Math.PI * radius * radius) - 1.0);
        assertTrue(fine < coarse,
                "halving the cell twice should improve the area: " + coarse + " -> " + fine);
        assertTrue(fine < 0.01, "a fine disc should be within one per cent, was " + fine);
    }

    @Test
    @DisplayName("a hole is a hole, and its winding does not matter")
    void holesAreSubtracted() {
        Outline plate = Outline.rectangle(0.0, 0.0, 20.0 * H, 20.0 * H);
        Outline bored = plate.withHole(Outline.rectangle(8.0 * H, 8.0 * H, 12.0 * H, 12.0 * H));

        assertEquals(400.0 * H * H, plate.area(), 1e-18);
        assertEquals(384.0 * H * H, bored.area(), 1e-18);
        assertEquals(400 - 16, bored.mesh(H).elementCount);
        assertFalse(bored.contains(10.0 * H, 10.0 * H), "the middle of the hole is outside");
        assertTrue(bored.contains(2.0 * H, 2.0 * H), "the corner of the plate is inside");
    }

    @Test
    @DisplayName("two shapes meshed together weld where they touch and not where they do not")
    void sharedLatticeWelds() {
        Outline left = Outline.rectangle(0.0, 0.0, 5.0 * H, 5.0 * H);
        Outline touching = Outline.rectangle(5.0 * H, 0.0, 10.0 * H, 5.0 * H);
        Outline apart = Outline.rectangle(9.0 * H, 0.0, 14.0 * H, 5.0 * H);

        // Welded: the shared column of six nodes is counted once, so 11x6 rather than 2x(6x6).
        QuadMesh joined = Outline.mesh(H, left, touching);
        assertEquals(50, joined.elementCount);
        assertEquals(11 * 6, joined.nodeCount);

        // Apart: two components, and the node count is the two blocks with nothing shared.
        QuadMesh separate = Outline.mesh(H, left, apart);
        assertEquals(50, separate.elementCount);
        assertEquals(2 * 6 * 6, separate.nodeCount);
    }

    @Test
    @DisplayName("a turned shape keeps its area and lands on the anvil")
    void rigidTransformsPreserveArea() {
        Outline bracket = Outline.of(0.0, 0.0, 30.0 * H, 0.0, 30.0 * H, 12.0 * H,
                12.0 * H, 12.0 * H, 12.0 * H, 40.0 * H, 0.0, 40.0 * H);
        double[] c = bracket.centre();
        Outline turned = bracket.rotated(Math.toRadians(18.0), c[0], c[1]).restingOn(0.0);

        assertEquals(bracket.area(), turned.area(), 1e-12 * bracket.area());
        assertEquals(0.0, turned.bounds()[1], 1e-18, "the lowest point must sit on z = 0");

        // Turning a shape changes how many cells its centres fall in, but not by much, and
        // never to nothing.
        assertTrue(Math.abs(turned.mesh(H).elementCount - bracket.mesh(H).elementCount) < 30,
                "a rotation should not change the element count materially");
    }

    @Test
    @DisplayName("a shape smaller than a cell is refused rather than meshed to nothing")
    void refusesWhatItCannotMesh() {
        // The failure this catches in practice is coordinates in millimetres rather than
        // metres, which produces a shape a thousand times too small and an empty mesh.
        assertThrows(IllegalArgumentException.class,
                () -> Outline.rectangle(0.0, 0.0, 1.0e-6, 1.0e-6).mesh(H));
        assertThrows(IllegalArgumentException.class,
                () -> Outline.rectangle(0.0, 0.0, H, H).mesh(0.0));
        assertThrows(IllegalArgumentException.class, () -> Outline.of(0.0, 0.0, 1.0, 1.0));
    }

    private static double area(QuadMesh mesh) {
        double a = 0.0;
        for (int e = 0; e < mesh.elementCount; e++) a += mesh.signedArea(e);
        return a;
    }
}
