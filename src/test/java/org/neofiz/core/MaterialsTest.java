package org.neofiz.core;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.neofiz.mesh.QuadMesh;
import org.neofiz.solver.ExplicitSolver;
import org.neofiz.solver.Integration;
import org.neofiz.solver.Kinematics;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * One mesh, more than one substance.
 *
 * <p>The property that matters most here is the one that is easiest to lose: a mesh of one
 * material has to give <b>bit-for-bit</b> the same answer through the new indirection as it
 * did through the old constant. Every validation gate in this project is a single-material
 * run, so if that identity holds then the gates are still testing what they were written to
 * test, and if it does not then every number in the README has quietly moved.
 */
class MaterialsTest {

    private static final Material SOFT =
            new Material("soft", 70.0e9, 0.33, 2700.0, 200.0e6, 1.0e9, 0.0);
    private static final Material HARD =
            new Material("hard", 210.0e9, 0.30, 7850.0, 700.0e6, 3.0e9, 0.0);

    private static QuadMesh wall() {
        return QuadMesh.cylinderWall(24.0e-3, 26.0e-3, 4.0e-3, 4, 4);
    }

    private static ExplicitSolver run(Object materials, int steps) {
        QuadMesh mesh = wall();
        ExplicitSolver s = materials instanceof Materials m
                ? new ExplicitSolver(mesh, m, Formulation.AXISYMMETRIC,
                        Integration.REDUCED, Kinematics.FINITE_STRAIN, 1.0, 0.5)
                : new ExplicitSolver(mesh, (Material) materials, Formulation.AXISYMMETRIC,
                        Integration.REDUCED, Kinematics.FINITE_STRAIN, 1.0, 0.5);
        s.setPressureRamp(300.0e6, 2.0e-6);
        s.run(steps);
        return s;
    }

    @Test
    @DisplayName("a uniform map reproduces the single-material solve exactly")
    void uniformIsBitIdentical() {
        // Not "to a tolerance". The same Material object reaches the same kernel by a
        // different route, so any difference at all would mean the route changed the
        // arithmetic, and every gate in the project is a single-material run.
        ExplicitSolver direct = run(SOFT, 400);
        ExplicitSolver mapped = run(Materials.uniform(SOFT), 400);

        assertEquals(direct.timestep(), mapped.timestep(), 0.0, "the CFL step must not move");
        double[] a = direct.stress(), b = mapped.stress();
        for (int i = 0; i < a.length; i++) {
            assertEquals(a[i], b[i], 0.0, "stress component " + i + " differs");
        }
        assertEquals(direct.strainEnergy(), mapped.strainEnergy(), 0.0);
        assertEquals(direct.plasticDissipation(), mapped.plasticDissipation(), 0.0);
    }

    @Test
    @DisplayName("a mixed mesh behaves like neither of its materials alone")
    void mixingChangesTheAnswer() {
        // The cheapest possible evidence that the paint is reaching the kernel rather than
        // being stored and ignored: half a wall of each has to land between the two.
        QuadMesh mesh = wall();
        Materials half = Materials.byElement(mesh.elementCount,
                e -> mesh.centroidRadius(e) < 25.0e-3 ? HARD : SOFT);

        double soft = run(Materials.uniform(SOFT), 400).strainEnergy();
        double hard = run(Materials.uniform(HARD), 400).strainEnergy();
        double mixed = run(half, 400).strainEnergy();

        assertTrue(Math.min(soft, hard) < mixed && mixed < Math.max(soft, hard),
                "mixed " + mixed + " should lie between " + soft + " and " + hard);
    }

    @Test
    @DisplayName("the timestep is sized by the fastest material on the mesh")
    void cflTakesTheFastestWave() {
        // One global step means the worst pairing has to be assumed: the shortest edge on the
        // mesh and the fastest wave on it, whether or not they are the same element. A map
        // that averaged, or took the first element's material, would run unstably wherever the
        // fast material actually is.
        //
        // Which of these two is faster is worked out rather than assumed, because the obvious
        // guess is wrong: c = sqrt((lambda + 2mu) / rho), and the soft material here is the
        // aluminium-like one, whose density falls further than its stiffness does. The stiffer
        // material is the slower one, and a test written around the stiff-means-fast intuition
        // would have passed for the wrong reason or failed for a right one.
        QuadMesh mesh = wall();
        boolean softIsFaster = SOFT.dilatationalWaveSpeed() > HARD.dilatationalWaveSpeed();
        Material fastest = softIsFaster ? SOFT : HARD;
        Material slower = softIsFaster ? HARD : SOFT;
        Materials mixed = Materials.byElement(mesh.elementCount, e -> e == 0 ? fastest : slower);

        assertEquals(fastest.dilatationalWaveSpeed(), mixed.fastestWave(), 0.0);
        assertEquals(run(Materials.uniform(fastest), 1).timestep(), run(mixed, 1).timestep(),
                0.0, "one fast element must set the step for the whole mesh");
        assertTrue(run(mixed, 1).timestep() < run(Materials.uniform(slower), 1).timestep(),
                "a single fast element has to cost the whole mesh its longer step");
    }

    @Test
    @DisplayName("a map that does not fit its mesh is refused at construction")
    void shortMapsAreRefused() {
        // Otherwise the failure is an array index out of bounds, in a parallel worker, several
        // hundred timesteps into a run that has already cost a minute.
        QuadMesh mesh = wall();
        Materials tooFew = Materials.of(SOFT, HARD);
        assertThrows(IllegalArgumentException.class,
                () -> new ExplicitSolver(mesh, tooFew, Formulation.AXISYMMETRIC,
                        Integration.REDUCED, Kinematics.SMALL_STRAIN, 1.0, 0.5));
        assertThrows(IllegalArgumentException.class, () -> Materials.of(SOFT, null));
        assertThrows(IllegalArgumentException.class, Materials::of);
    }

    @Test
    @DisplayName("distinct materials are reported once each, in first-appearance order")
    void distinctIsOrderedAndDeduplicated() {
        // What the plane-stress and softening checks iterate over. A duplicate would only cost
        // time; a missing entry would let an unchecked material through.
        Materials m = Materials.of(SOFT, SOFT, HARD, SOFT, HARD);
        assertEquals(2, m.distinct().size());
        assertEquals(SOFT, m.distinct().get(0));
        assertEquals(HARD, m.distinct().get(1));
        assertEquals(5, m.elementCount());
        assertFalse(m.isUniform());

        Materials one = Materials.uniform(HARD);
        assertTrue(one.isUniform());
        assertEquals(-1, one.elementCount(), "a uniform map fits a mesh of any size");
        assertEquals(1, one.distinct().size());
    }
}
