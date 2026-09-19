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
 * Finite-strain kinematics: objectivity under rotation, and agreement with the small-strain
 * formulation where the two must agree.
 *
 * <p>Objectivity is the property most worth testing to machine precision, because its failure
 * mode is invisible. A non-objective update produces stress out of pure rotation -- no load,
 * no deformation, just spin -- and it grows with accumulated rotation rather than with
 * anything a modeller would think to check. It shows up in folding and spalling regions,
 * late in a run, in exactly the places where there is no closed form to compare against and
 * the answer looks plausible either way.
 */
class FiniteStrainTest {

    private static final Material STEEL = Material.STEEL_4340;

    private static QuadMesh unitSquare() {
        return new QuadMesh(
                new double[]{0.0, 1.0, 1.0, 0.0},
                new double[]{0.0, 0.0, 1.0, 1.0},
                new int[]{0, 1, 2, 3},
                new int[0]);
    }

    private static ExplicitSolver finiteSolver(Material m) {
        return new ExplicitSolver(unitSquare(), m, Formulation.PLANE_STRAIN,
                Integration.REDUCED, Kinematics.FINITE_STRAIN, 1.0, 0.5);
    }

    /**
     * Drives the mesh to a prescribed set of node positions over one step, by setting the
     * velocity that gets there. This is how a rigid rotation is imposed exactly: the nodes
     * land on the true rotated positions, and the solver is left to work out for itself that
     * no stretching happened.
     */
    private static void moveTo(ExplicitSolver s, QuadMesh mesh, double[] targetR, double[] targetZ) {
        double dt = s.timestep();
        for (int i = 0; i < 4; i++) {
            double nowR = mesh.r[i] + s.radialDisplacement()[i];
            double nowZ = mesh.z[i] + s.axialDisplacement()[i];
            s.setVelocity(i, (targetR[i] - nowR) / dt, (targetZ[i] - nowZ) / dt);
        }
        s.advance();
    }

    @Test
    @DisplayName("rigid rotation transports the stress tensor exactly and invents no strain")
    void rigidRotationIsExactlyObjective() {
        QuadMesh mesh = unitSquare();
        ExplicitSolver s = finiteSolver(STEEL);

        // Put the element under a stress state with a clear orientation: stretch it along z
        // only, so sigma_zz and sigma_rr differ and a rotation has something to mix.
        double[] r = {0.0, 1.0, 1.0, 0.0};
        double[] z = {0.0, 0.0, 1.002, 1.002};
        moveTo(s, mesh, r, z);

        double[] before = s.centroidStress(0).clone();
        double misesBefore = vonMises(before);
        double pressureBefore = pressure(before);
        assertTrue(misesBefore > 1.0e6, "the element must actually be stressed");

        // Now rotate rigidly through a quarter turn, in many small steps. No stretching at
        // any point, so the rate of deformation is identically zero and the only thing that
        // may change is the frame the components are written in.
        final int steps = 360;
        final double total = Math.PI / 2.0;
        for (int k = 1; k <= steps; k++) {
            double a = total * k / steps;
            double c = Math.cos(a), sn = Math.sin(a);
            double[] tr = new double[4], tz = new double[4];
            for (int i = 0; i < 4; i++) {
                tr[i] = c * r[i] - sn * z[i];
                tz[i] = sn * r[i] + c * z[i];
            }
            moveTo(s, mesh, tr, tz);
        }

        double[] after = s.centroidStress(0);

        // Invariants are frame-independent: they must not have moved at all.
        assertEquals(misesBefore, vonMises(after), misesBefore * 1e-9,
                "equivalent stress changed under pure rotation; the update is not objective "
                        + "and is manufacturing stress out of spin");
        assertEquals(pressureBefore, pressure(after), Math.abs(pressureBefore) * 1e-9 + 1.0,
                "pressure changed under pure rotation");

        // And the components must be exactly the quarter-turn transform of the originals:
        // a quarter turn swaps rr and zz and flips the sign of the shear.
        assertEquals(before[1], after[0], misesBefore * 1e-9, "sigma_rr should be the old sigma_zz");
        assertEquals(before[0], after[1], misesBefore * 1e-9, "sigma_zz should be the old sigma_rr");
        assertEquals(before[2], after[2], misesBefore * 1e-9, "the hoop component does not rotate");

        assertEquals(0.0, s.plasticStrain()[0], 0.0,
                "rotation must not yield anything; if it does, the spurious stress got large "
                        + "enough to cross the yield surface and is now permanent");
    }

    @Test
    @DisplayName("a non-objective update would visibly fail the same test")
    void smallStrainIsNotObjectiveAndThisMatters() {
        // The contrast is the point. The small-strain formulation is not wrong -- it is
        // correct to first order in strain and the entire M0 gate set rests on it -- but it
        // has no notion of the mesh turning, so the same rotation corrupts it. Without this
        // comparison the objectivity test above only proves the code does something.
        QuadMesh mesh = unitSquare();
        ExplicitSolver s = new ExplicitSolver(mesh, STEEL, Formulation.PLANE_STRAIN,
                Integration.REDUCED, Kinematics.SMALL_STRAIN, 1.0, 0.5);

        double[] r = {0.0, 1.0, 1.0, 0.0};
        double[] z = {0.0, 0.0, 1.002, 1.002};
        moveTo(s, mesh, r, z);
        double misesBefore = vonMises(s.centroidStress(0));

        for (int k = 1; k <= 360; k++) {
            double a = (Math.PI / 2.0) * k / 360;
            double c = Math.cos(a), sn = Math.sin(a);
            double[] tr = new double[4], tz = new double[4];
            for (int i = 0; i < 4; i++) {
                tr[i] = c * r[i] - sn * z[i];
                tz[i] = sn * r[i] + c * z[i];
            }
            moveTo(s, mesh, tr, tz);
        }

        double misesAfter = vonMises(s.centroidStress(0));
        assertTrue(Math.abs(misesAfter - misesBefore) > 0.1 * misesBefore,
                "small strain was expected to mangle this rotation, and did not. Either the "
                        + "kinematics flag is not reaching the kernel, or the rotation imposed "
                        + "here is too small to be a real test of the finite-strain path");
    }

    @Test
    @DisplayName("the quarter turn is not an arbitrary choice: a full turn hides the defect")
    void afullTurnWouldHideIt() {
        // Measured at a quarter turn the non-objective error is about 500x the stress. Taken
        // all the way round it is 1e-12, indistinguishable from exact -- the linearised
        // strain increments of a rigid rotation cancel over a complete revolution, so the
        // element arrives back where it started with its error undone.
        //
        // The error is periodic in rotation, not monotonic. An objectivity test that happened
        // to use 360 degrees would certify a broken update as perfect. Pinning that here
        // stops the angle in the test above from being "tidied up" later.
        double quarter = rotationError(Kinematics.SMALL_STRAIN, 0.25);
        double full = rotationError(Kinematics.SMALL_STRAIN, 1.00);

        assertTrue(quarter > 1.0, "a quarter turn must expose the defect; got " + quarter);
        assertTrue(full < 1e-6,
                "a full turn must hide it, which is the point; got " + full);

        assertTrue(rotationError(Kinematics.FINITE_STRAIN, 0.25) < 1e-9,
                "finite strain must be exact at every angle, not merely at the convenient ones");
        assertTrue(rotationError(Kinematics.FINITE_STRAIN, 3.30) < 1e-9,
                "including after several turns and a fraction");
    }

    /** Relative change in equivalent stress after a rigid rotation through {@code turns}. */
    private static double rotationError(Kinematics k, double turns) {
        QuadMesh mesh = unitSquare();
        ExplicitSolver s = new ExplicitSolver(mesh, STEEL, Formulation.PLANE_STRAIN,
                Integration.REDUCED, k, 1.0, 0.5);

        double[] r = {0.0, 1.0, 1.0, 0.0};
        double[] z = {0.0, 0.0, 1.002, 1.002};
        moveTo(s, mesh, r, z);
        double before = vonMises(s.centroidStress(0));

        final int steps = 400;
        for (int i = 1; i <= steps; i++) {
            double a = 2.0 * Math.PI * turns * i / steps;
            double c = Math.cos(a), sn = Math.sin(a);
            double[] tr = new double[4], tz = new double[4];
            for (int j = 0; j < 4; j++) {
                tr[j] = c * r[j] - sn * z[j];
                tz[j] = sn * r[j] + c * z[j];
            }
            moveTo(s, mesh, tr, tz);
        }
        return Math.abs(vonMises(s.centroidStress(0)) - before) / before;
    }

    @Test
    @DisplayName("isochoric compression to 50% yields exactly sy/sqrt(3) at exactly zero pressure")
    void largeIsochoricDeformationHitsTheExactYieldState() {
        // Squash in z by s, stretch in r by 1/s, so the Jacobian is exactly one at every
        // instant. Two exact consequences follow, and between them they pin down the whole
        // finite-strain stack:
        //
        //   pressure = 0     the volumetric strain is ln(det F), and det F is exactly 1.
        //                    Small strain gets this wrong, because its linearised trace of
        //                    a large isochoric deformation is not zero.
        //   sigma_zz = -sy / sqrt(3)
        //                    a deviatoric plane-strain state at yield. Only the deviator can
        //                    be nonzero, symmetry forces sigma_tt to sit at the mean, and
        //                    von Mises then fixes the magnitude with no freedom left.
        final double sy = 800.0e6;
        Material steel = Material.STEEL_4340.yielding(sy, 0.0);
        QuadMesh mesh = unitSquare();
        ExplicitSolver s = new ExplicitSolver(mesh, steel, Formulation.PLANE_STRAIN,
                Integration.REDUCED, Kinematics.FINITE_STRAIN, 1.0, 0.5);

        final int n = 4000;
        final double finalScale = 0.5;      // 50% height, 100% wider
        for (int i = 1; i <= n; i++) {
            double c = 1.0 + (finalScale - 1.0) * i / n;
            moveTo(s, mesh, new double[]{0, 1 / c, 1 / c, 0}, new double[]{0, 0, c, c});
        }

        double[] sig = s.centroidStress(0);
        assertTrue(s.plasticStrain()[0] > 0.3,
                "50% compression must drive substantial plastic strain; got " + s.plasticStrain()[0]);

        assertEquals(0.0, pressure(sig), sy * 1e-3,
                "an exactly isochoric path must build no pressure. A nonzero value means the "
                        + "volumetric strain is being taken as the linearised trace rather than "
                        + "the log of the Jacobian, or that plastic flow is leaking volume");
        assertEquals(-sy / Math.sqrt(3.0), sig[1], sy * 5e-3,
                "axial stress must sit at the plane-strain yield value");
        assertEquals(sy, vonMises(sig), sy * 5e-3,
                "and the state must be exactly on the yield surface");
    }

    @ParameterizedTest
    @EnumSource(Integration.class)
    @DisplayName("finite and small strain agree when the strain really is small")
    void formulationsAgreeAtSmallStrain(Integration rule) {
        // 0.02% strain. The two formulations differ at order strain, so they must agree to
        // far better than the gates care about. If they do not, the finite-strain path has a
        // defect that has nothing to do with large deformation.
        double stretch = 1.0 + 2.0e-4;

        double[] small = stressAfterStretch(Kinematics.SMALL_STRAIN, rule, stretch);
        double[] finite = stressAfterStretch(Kinematics.FINITE_STRAIN, rule, stretch);

        double scale = vonMises(small);
        assertTrue(scale > 0.0, "the stretch must load the element");
        for (int c = 0; c < 4; c++) {
            assertEquals(small[c], finite[c], scale * 1e-3,
                    "component " + c + " differs between formulations at 0.02% strain");
        }
    }

    @Test
    @DisplayName("and diverge, in the right direction, when it is not")
    void formulationsDivergeAtLargeStrain() {
        // 20% stretch. The finite-strain answer is the smaller one: it measures strain
        // against the configuration the material is actually in, and reports the true stress
        // on the area that actually exists, where small strain keeps dividing by the original.
        double[] small = stressAfterStretch(Kinematics.SMALL_STRAIN, Integration.REDUCED, 1.20);
        double[] finite = stressAfterStretch(Kinematics.FINITE_STRAIN, Integration.REDUCED, 1.20);

        assertTrue(vonMises(finite) < vonMises(small),
                "at 20% stretch the two must disagree, with finite strain lower: got finite "
                        + vonMises(finite) / 1e6 + " MPa against small " + vonMises(small) / 1e6);
        assertTrue(Math.abs(vonMises(finite) - vonMises(small)) > 0.05 * vonMises(small),
                "and the gap must be substantial at 20%, not a rounding difference");
    }

    private static double[] stressAfterStretch(Kinematics k, Integration rule, double stretch) {
        QuadMesh mesh = unitSquare();
        ExplicitSolver s = new ExplicitSolver(mesh, STEEL, Formulation.PLANE_STRAIN,
                rule, k, 1.0, 0.5);

        // Pull along z in many increments, so the finite-strain path integrates a real path
        // rather than taking one enormous step.
        final int n = 500;
        for (int i = 1; i <= n; i++) {
            double h = 1.0 + (stretch - 1.0) * i / n;
            moveTo(s, mesh, new double[]{0, 1, 1, 0}, new double[]{0, 0, h, h});
        }
        return s.centroidStress(0);
    }

    private static double vonMises(double[] sig) {
        double p = (sig[0] + sig[1] + sig[2]) / 3.0;
        double a = sig[0] - p, b = sig[1] - p, c = sig[2] - p;
        return Math.sqrt(1.5 * (a * a + b * b + c * c + 2.0 * sig[3] * sig[3]));
    }

    private static double pressure(double[] sig) {
        return -(sig[0] + sig[1] + sig[2]) / 3.0;
    }
}
