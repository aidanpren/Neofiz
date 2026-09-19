package org.neofiz.solver;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.neofiz.core.Formulation;
import org.neofiz.core.JohnsonCook;
import org.neofiz.core.Material;
import org.neofiz.mesh.QuadMesh;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Tests for progressive softening and its regularisation.
 *
 * <p>The claim under test is a single sentence -- <em>breaking a part costs G_f per unit of
 * crack area, whatever the mesh</em> -- and it is the kind of claim that fails quietly. A
 * softening model with the slope wrong is not unstable and does not produce implausible
 * stresses; it produces a part that gets weaker every time the mesh is refined, which reads as
 * convergence failing to arrive rather than as a bug. So the tests here are identities where
 * an identity is available and refinement studies where one is not.
 *
 * <p>The strain path used throughout stretches in z and contracts equally in r and theta. It
 * is exactly volume-preserving, so no pressure builds and the deviator stays on one axis,
 * which fixes the plastic flow direction and therefore fixes the band orientation -- the band
 * comes out normal to z and its width is the element's z extent. That makes the answer a
 * number that can be written down rather than one that has to be believed.
 */
class DamageTest {

    private static final double E = 205.0e9;
    private static final double NU = 0.29;
    private static final double SY = 800.0e6;
    private static final double MU = E / (2.0 * (1.0 + NU));
    private static final double LAMBDA = E * NU / ((1.0 + NU) * (1.0 - 2.0 * NU));

    /** 50 kJ/m^2, a ductile structural steel. */
    private static final double GF = 50.0e3;
    private static final double EPS_F = 0.20;

    private static final double DT = 1.0e-8;

    // ---------------------------------------------------------------- single point

    /** Geometry for one rectangular element of the given side lengths. */
    private static double[] box(double dr, double dz) {
        return new double[] {dr * dz, dr, 0.0, 0.0, dz};
    }

    /** One material point driven along a fixed deviatoric axis. */
    private static final class Point {
        final double[] sig = new double[4];
        final double[] back = new double[4];
        final double[] epsP = new double[1];
        final double[] work = new double[1];
        final J2.Flow flow;
        final Damage damage;

        Point(J2.Flow flow, Damage damage) {
            this.flow = flow;
            this.damage = damage;
        }

        static Point linear(double hIso, double dr, double dz) {
            return new Point(J2.Flow.linear(SY, hIso, 0.0),
                    new Damage(GF, E, new double[] {EPS_F}, box(dr, dz)));
        }

        /** One increment of axial stretch at constant volume. */
        void pull(double d) {
            J2.update(sig, back, epsP, work, 0, -0.5 * d, d, -0.5 * d, 0.0,
                    LAMBDA, MU, flow, DT, damage);
        }

        double mises() {
            return J2.vonMises(sig, 0);
        }

        double damageNow() {
            return damage.at(0, epsP[0]);
        }

        /** Plastic work and plastic strain as they stood on the step damage began. */
        double workAtOnset, strainAtOnset;
        /** Steps spent on the softening branch. */
        int softeningSteps;

        /** Drives to full damage and returns the plastic work spent from onset onwards. */
        double softenToFailure(double step) {
            while (!damage.hasBegun(0)) pull(step);
            workAtOnset = work[0];
            strainAtOnset = epsP[0];
            softeningSteps = 0;
            while (damageNow() < 1.0) {
                pull(step);
                softeningSteps++;
            }
            // One more, to prove a fully damaged point stays put rather than going somewhere.
            pull(step);
            return work[0] - workAtOnset;
        }
    }

    @Test
    @DisplayName("the softening branch is a triangle whose area is exactly G_f / h")
    void theTriangleIsExact() {
        Point p = Point.linear(2.0e9, 0.5e-3, 0.5e-3);
        while (!p.damage.hasBegun(0)) p.pull(1.0e-5);

        double h = p.damage.bandWidth(0);
        double sf = p.damage.strengthAtOnset(0);
        double range = Damage.softeningRange(GF, sf, h);

        // The closed form at full damage, from the state alone.
        double full = p.damage.dissipation(0, EPS_F + range);
        assertEquals(GF / h, full, GF / h * 1e-14,
                "a fully softened point must have spent exactly G_f per unit band width");
        assertEquals(0.5 * sf * range, full, full * 1e-14,
                "and that must be the area of the triangle");

        // Past full damage it stops accumulating rather than going negative.
        assertEquals(full, p.damage.dissipation(0, EPS_F + 10.0 * range), full * 1e-14,
                "dissipation must stop at the end of the softening branch");
        assertEquals(0.0, p.damage.dissipation(0, EPS_F - 0.01), 0.0,
                "and must be zero before onset");
    }

    @Test
    @DisplayName("the band is oriented by the flow, so its width is the z extent, not sqrt(A)")
    void theBandIsOrientedByTheFlow() {
        // Same z extent, wildly different aspect ratios. An axial pull opens a band normal to
        // z, whose width is dz whatever dr does.
        double dz = 0.4e-3;
        for (double dr : new double[] {0.05e-3, 0.4e-3, 4.0e-3}) {
            Point p = Point.linear(2.0e9, dr, dz);
            while (!p.damage.hasBegun(0)) p.pull(1.0e-5);
            assertEquals(dz, p.damage.bandWidth(0), dz * 1e-12,
                    "band width must follow the crack normal, not the element area");
        }

        // What the usual shortcut would have said, for the record: sqrt(A) is right only on
        // the square and is out by sqrt(aspect ratio) either side of it.
        assertEquals(Math.sqrt(0.05e-3 * dz) / dz, Math.sqrt(0.05 / 0.4), 1e-12);
        assertTrue(Math.sqrt(4.0e-3 * dz) / dz > 3.1,
                "on a 10:1 element sqrt(A) would over-estimate the band by more than three");
    }

    @Test
    @DisplayName("the softening return lands exactly on the softened surface")
    void softeningIsExactlyOnTheSurface() {
        Point p = Point.linear(2.0e9, 0.5e-3, 0.5e-3);
        while (!p.damage.hasBegun(0)) p.pull(1.0e-5);

        double sf = p.damage.strengthAtOnset(0);
        double slope = p.damage.softeningModulus(0);
        assertTrue(slope < 0.0, "the softening modulus must be negative");

        int checked = 0;
        while (p.damageNow() < 1.0) {
            p.pull(1.0e-5);
            double expected = Math.max(0.0, sf + slope * (p.epsP[0] - EPS_F));
            assertEquals(expected, p.mises(), SY * 1e-13,
                    "damage is negative hardening; the return map is still exact");
            checked++;
        }
        assertTrue(checked > 100, "the softening branch must take real steps, took " + checked);
    }

    @Test
    @DisplayName("a fully damaged point has no deviatoric strength and still carries pressure")
    void fullDamageLeavesPressureOnly() {
        Point p = Point.linear(2.0e9, 0.5e-3, 0.5e-3);
        p.softenToFailure(1.0e-5);

        assertEquals(1.0, p.damageNow(), 0.0, "must be fully damaged");
        assertEquals(0.0, p.mises(), SY * 1e-12,
                "a fully damaged point must carry no von Mises stress at all");

        // Now compress it. The volumetric response is untouched, which is the documented
        // limitation: this is a fluid, not a crack.
        double d = -1.0e-4;
        J2.update(p.sig, p.back, p.epsP, p.work, 0, d, d, d, 0.0, LAMBDA, MU, p.flow, DT,
                p.damage);
        assertEquals(-(3.0 * LAMBDA + 2.0 * MU) * d, J2.pressure(p.sig, 0),
                Math.abs((3.0 * LAMBDA + 2.0 * MU) * d) * 1e-12,
                "hydrostatic stiffness must survive full damage");
        assertEquals(0.0, p.mises(), SY * 1e-12, "and deviatoric strength must not come back");
    }

    @Test
    @DisplayName("the solver undershoots the closed form by exactly one step in n")
    void theSolverUndershootsByOneStepInN() {
        // The plastic work accumulator takes the stress at the end of each step, which on a
        // falling branch is a right Riemann sum of a straight line. That has a known answer:
        // the sum is short of the integral by one step out of the whole branch, so the
        // relative deficit is 1/n and nothing else. Measuring it against the closed form over
        // the *same* strain interval is what isolates it -- comparing against G_f/h instead
        // would fold in the overshoot at onset, where the crossing increment is taken on the
        // undamaged surface, and the two errors are not the same size.
        double previous = Double.NaN;
        for (double step : new double[] {4.0e-5, 2.0e-5, 1.0e-5}) {
            Point p = Point.linear(2.0e9, 0.5e-3, 0.5e-3);
            double spent = p.softenToFailure(step);

            double exact = p.damage.dissipation(0, p.epsP[0])
                    - p.damage.dissipation(0, p.strainAtOnset);
            double deficit = (exact - spent) / exact;

            assertEquals(1.0 / p.softeningSteps, deficit, 0.2 / p.softeningSteps,
                    "the deficit must be one step out of " + p.softeningSteps);
            if (!Double.isNaN(previous)) {
                assertTrue(deficit < 0.6 * previous,
                        "halving the step must roughly halve the deficit: " + previous
                                + " then " + deficit);
            }
            previous = deficit;
        }
        assertTrue(previous < 0.01, "and it must be small: " + previous);
    }

    @Test
    @DisplayName("energy per unit crack area is fixed; energy per unit volume is not")
    void theFractureEnergyIsPerUnitArea() {
        // Four band widths over a 27-fold range. The regularised quantity is the one that
        // holds still; the unregularised one is shown moving beside it, because "the energy
        // did not change" is only evidence if something else did.
        double[] widths = {0.1e-3, 0.3e-3, 0.9e-3, 2.7e-3};
        double[] perVolume = new double[widths.length];
        double[] perArea = new double[widths.length];

        for (int i = 0; i < widths.length; i++) {
            Point p = Point.linear(2.0e9, 0.5e-3, widths[i]);
            double spent = p.softenToFailure(2.0e-6);
            assertEquals(0, p.damage.clampedPoints(), "no element here may hit the clamp");
            perVolume[i] = spent;
            perArea[i] = spent * p.damage.bandWidth(i == 0 ? 0 : 0);
        }

        for (int i = 1; i < widths.length; i++) {
            assertEquals(perArea[0], perArea[i], perArea[0] * 0.02,
                    "energy per unit crack area must not move with the band width");
        }
        assertEquals(GF, perArea[0], GF * 0.02, "and it must be G_f");

        // The same runs, read per unit volume: a 27-fold spread, which is exactly what would
        // have gone into the answer if the slope had been a material constant.
        assertEquals(widths[3] / widths[0], perVolume[0] / perVolume[3],
                widths[3] / widths[0] * 0.03,
                "energy per unit volume must fall as the band widens, in inverse proportion");
    }

    @Test
    @DisplayName("an element too large to regularise is clamped and counted, not tolerated")
    void snapBackIsClampedAndCounted() {
        // The limit depends on the flow stress at onset, which includes the hardening picked
        // up on the way, so it has to be measured rather than assumed from sy0.
        Point probe = Point.linear(2.0e9, 0.5e-3, 0.5e-3);
        while (!probe.damage.hasBegun(0)) probe.pull(1.0e-5);
        double sf = probe.damage.strengthAtOnset(0);
        double limit = Damage.largestRegularisableElement(GF, E, sf);

        Point fine = Point.linear(2.0e9, 0.5e-3, 0.5 * limit);
        while (!fine.damage.hasBegun(0)) fine.pull(1.0e-5);
        assertEquals(0, fine.damage.clampedPoints(), "half the limit must not clamp");
        assertNotEquals(-E, fine.damage.softeningModulus(0));

        Point coarse = Point.linear(2.0e9, 0.5e-3, 2.0 * limit);
        while (!coarse.damage.hasBegun(0)) coarse.pull(1.0e-5);
        assertEquals(1, coarse.damage.clampedPoints(), "twice the limit must clamp");
        assertEquals(-E, coarse.damage.softeningModulus(0), 0.0,
                "the clamp must land on exactly -E, so counting it needs no tolerance");

        // And the clamped point still returns to a sane stress: the denominator
        // 2mu + (2/3)(-E) is small but positive, which is the whole reason for the clamp.
        double before = coarse.mises();
        coarse.pull(1.0e-5);
        assertTrue(Double.isFinite(coarse.mises()) && coarse.mises() < before,
                "a clamped point must still soften monotonically, got " + coarse.mises());
        assertTrue(2.0 * MU - (2.0 / 3.0) * E > 0.0,
                "and the clamped return map must be non-singular");
    }

    @Test
    @DisplayName("Johnson-Cook softens through the closed form, with no iteration left to do")
    void johnsonCookSoftensWithoutIterating() {
        JohnsonCook law = JohnsonCook.STEEL_4340;
        Material steel = Material.STEEL_4340.withJohnsonCook(law);
        Damage damage = new Damage(GF, E, new double[] {EPS_F}, box(0.5e-3, 0.5e-3));
        Point p = new Point(J2.Flow.of(steel), damage);

        while (!damage.hasBegun(0)) p.pull(1.0e-5);
        double sf = damage.strengthAtOnset(0);
        double slope = damage.softeningModulus(0);

        int checked = 0;
        while (p.damageNow() < 1.0) {
            p.pull(1.0e-5);
            double expected = Math.max(0.0, sf + slope * (p.epsP[0] - EPS_F));
            assertEquals(expected, p.mises(), sf * 1e-13,
                    "a frozen flow stress removes the nonlinearity the solve existed for");
            checked++;
        }
        assertTrue(checked > 50, "the softening branch must take real steps");

        // The frozen strength must be the quasi-static flow stress the law had reached, not
        // A and not the rate-hardened value.
        assertTrue(sf > law.a(), "the material must have hardened before it failed");
    }

    @Test
    @DisplayName("a failure strain that is never reached changes nothing, bit for bit")
    void anUnreachableFailureStrainChangesNothing() {
        Point plain = new Point(J2.Flow.linear(SY, 2.0e9, 1.0e9), null);
        Point armed = new Point(J2.Flow.linear(SY, 2.0e9, 1.0e9),
                new Damage(GF, E, new double[] {10.0}, box(0.5e-3, 0.5e-3)));

        for (int i = 0; i < 4000; i++) {
            plain.pull(1.0e-5);
            armed.pull(1.0e-5);
        }
        assertTrue(armed.epsP[0] > 0.03, "the path must have yielded hard");
        for (int c = 0; c < 4; c++) {
            assertEquals(plain.sig[c], armed.sig[c], 0.0, "stress component " + c + " drifted");
            assertEquals(plain.back[c], armed.back[c], 0.0, "back stress " + c + " drifted");
        }
        assertEquals(plain.epsP[0], armed.epsP[0], 0.0);
        assertEquals(plain.work[0], armed.work[0], 0.0);
        assertEquals(0.0, armed.damage.begunFraction(), 0.0);
    }

    @Test
    @DisplayName("the constructor refuses states it cannot regularise")
    void validation() {
        double[] one = {EPS_F};
        assertThrows(IllegalArgumentException.class,
                () -> new Damage(0.0, E, one, box(1e-3, 1e-3)));
        assertThrows(IllegalArgumentException.class,
                () -> new Damage(GF, 0.0, one, box(1e-3, 1e-3)));
        assertThrows(IllegalArgumentException.class,
                () -> new Damage(GF, E, one, new double[] {1.0, 1.0}));
        assertThrows(IllegalArgumentException.class,
                () -> new Damage(GF, E, new double[] {0.0}, box(1e-3, 1e-3)));
        assertThrows(IllegalArgumentException.class,
                () -> new Damage(GF, E, one, box(0.0, 1e-3)));
    }

    // ---------------------------------------------------------------- whole mesh

    private static final double BORE = 20.0e-3;
    private static final double OUTER = 30.0e-3;
    private static final double HEIGHT = 4.0e-3;

    /**
     * Stretches a tube axially at constant prescribed velocity, with one axial row of elements
     * given a low failure strain, and returns the energy the softening took.
     *
     * <p>Prescribed rather than dynamic on purpose. Localisation is what a softening solve does
     * on its own, and letting it happen here would make the crack area an outcome of the run
     * instead of something the test fixes -- which is fine for the burst case and useless for
     * measuring an energy.
     */
    private static double ringCrackEnergy(int nr, int nz, int weakRow, Integration rule) {
        QuadMesh mesh = QuadMesh.cylinderWall(BORE, OUTER, HEIGHT, nr, nz);
        Material steel = Material.STEEL_4340.yielding(SY, 2.0e9);
        ExplicitSolver solver = new ExplicitSolver(mesh, steel, Formulation.AXISYMMETRIC,
                rule, 1.0, 0.5);

        double[] failure = new double[mesh.elementCount];
        java.util.Arrays.fill(failure, 10.0);
        for (int i = 0; i < nr; i++) failure[weakRow * nr + i] = EPS_F;
        solver.setDamage(GF, failure);

        // Uniform axial stretch: v_z proportional to z, so every element sees the same strain
        // increment and only the weak row can soften.
        double rate = 1.0e-3 / HEIGHT;          // per second, scaled below by the timestep
        for (int n = 0; n < mesh.nodeCount; n++) {
            solver.setVelocity(n, 0.0, rate * mesh.z[n] / solver.timestep());
        }

        for (int i = 0; i < 400000 && solver.failedFraction() < nr / (double) mesh.elementCount;
             i++) {
            solver.advance();
        }
        assertEquals(nr / (double) mesh.elementCount, solver.failedFraction(), 1e-12,
                "exactly the weak row must have failed, on a " + nr + "x" + nz + " mesh");
        assertEquals(0, solver.damage().clampedPoints(), "no point may hit the snap-back clamp");
        return solver.fractureDissipation();
    }

    @Test
    @DisplayName("a through-wall ring crack costs G_f times its area, on any mesh")
    void theRingCrackCostsItsArea() {
        // The crack is an annulus between bore and outside radius. Its area does not depend on
        // how many elements are used to span it, so neither may the energy -- and the 2 pi r
        // that makes every element here a ring has to cancel for that to be true.
        double exact = GF * Math.PI * (OUTER * OUTER - BORE * BORE);

        double[] measured = new double[3];
        int[] through = {2, 4, 8};
        for (int i = 0; i < through.length; i++) {
            measured[i] = ringCrackEnergy(through[i], 3, 1, Integration.REDUCED);
            assertEquals(exact, measured[i], exact * 1e-10,
                    "energy must be G_f times the annulus area at nr = " + through[i]);
        }
        assertEquals(measured[0], measured[2], measured[0] * 1e-10,
                "and must not move across a fourfold refinement through the wall");
    }

    @Test
    @DisplayName("refining along the axis narrows the band without changing the crack")
    void refiningAlongTheAxisChangesTheBandNotTheCrack() {
        // This is the refinement that would break an unregularised model. Halving the axial
        // element halves the band width, so the softening slope has to double to keep the
        // energy where it is. Nothing about the crack itself changes: same annulus, same area.
        double exact = GF * Math.PI * (OUTER * OUTER - BORE * BORE);
        for (int nz : new int[] {2, 4, 8, 16}) {
            assertEquals(exact, ringCrackEnergy(4, nz, nz / 2, Integration.REDUCED),
                    exact * 1e-10, "an eightfold axial refinement must not move the energy");
        }
    }

    @Test
    @DisplayName("full and reduced integration dissipate the same, exactly")
    void theQuadratureRuleIsNotAMaterialProperty() {
        // Four points per element instead of one, all sharing the element's band width. If
        // each point sized its own band the element would host two parallel cracks in each
        // direction and this would come out doubled.
        double reduced = ringCrackEnergy(4, 3, 1, Integration.REDUCED);
        double full = ringCrackEnergy(4, 3, 1, Integration.FULL);
        assertEquals(reduced, full, reduced * 1e-10,
                "the quadrature rule must not change the fracture toughness");
    }
}
