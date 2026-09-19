package org.neofiz.solver;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.neofiz.core.JohnsonCook;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Single-point tests for the return map when the flow stress is Johnson-Cook, where the
 * plastic multiplier has to be solved for rather than written down.
 *
 * <p>The exact closed form was one of the properties worth having about the linear law, so
 * the central test here is not that the solve produces a plausible answer -- it is that where
 * a closed form still exists, the solve <b>reproduces it to machine precision</b>. A
 * degenerate Johnson-Cook law with no hardening, no rate sensitivity and no heating is
 * perfect plasticity, for which the answer is known exactly. If the iteration cannot land on
 * that, nothing it produces on the real law is worth anything.
 */
class RateDependentReturnMapTest {

    private static final double E = 205.0e9;
    private static final double NU = 0.29;
    private static final double MU = E / (2.0 * (1.0 + NU));
    private static final double LAMBDA = E * NU / ((1.0 + NU) * (1.0 - 2.0 * NU));
    private static final double DENSITY = 7850.0;
    private static final double SY = 800.0e6;
    private static final double ROOM = 293.0;

    /** A Johnson-Cook law stripped down to whatever the test needs. */
    private static JohnsonCook law(double a, double b, double n, double c, double beta) {
        return new JohnsonCook(a, b, n, c, 1.03, 1.0, ROOM, 1793.0, 477.0, beta);
    }

    /** One material point driven through prescribed strain increments. */
    private static final class Point {
        final double[] sig = new double[4];
        final double[] back = new double[4];
        final double[] epsP = new double[1];
        final double[] work = new double[1];
        final J2.Flow flow;

        Point(J2.Flow flow) {
            this.flow = flow;
        }

        static Point johnsonCook(JohnsonCook jc) {
            return new Point(new J2.Flow(jc.a(), 0.0, 0.0, jc, jc.thermalCoefficient(DENSITY)));
        }

        static Point linear(double sy, double hIso) {
            return new Point(J2.Flow.linear(sy, hIso, 0.0));
        }

        double shear(double dGam, double dt) {
            return J2.update(sig, back, epsP, work, 0, 0, 0, 0, dGam, LAMBDA, MU, flow, dt);
        }

        double tau() {
            return sig[3];
        }

        double pressure() {
            return J2.pressure(sig, 0);
        }
    }

    @Test
    @DisplayName("with no hardening, rate sensitivity or heating, the solve matches the exact return")
    void degeneratesToTheClosedForm() {
        // B = 0 removes strain hardening, C = 0 removes rate hardening, beta = 0 stops the
        // plastic work from raising the temperature, so the thermal term stays pinned at 1.
        // What is left is perfect plasticity with flow stress exactly A -- a case whose
        // multiplier is (||xi|| - sqrt(2/3) sy) / 2mu with no iteration at all.
        Point solved = Point.johnsonCook(law(SY, 0.0, 0.26, 0.0, 0.0));
        Point exact = Point.linear(SY, 0.0);

        final double dt = 1.0e-8;
        for (int step = 0; step < 40; step++) {
            double a = solved.shear(2.0e-4, dt);
            double b = exact.shear(2.0e-4, dt);

            // Machine precision, not a tolerance. The safeguarded Newton has a convergence
            // tolerance, but on this law the residual is linear in the multiplier, so the
            // first Newton step is the answer and the tolerance never comes into it.
            assertEquals(b, a, Math.abs(b) * 1e-12 + 1e-300,
                    "multiplier differs from the closed form at step " + step);
            assertEquals(exact.tau(), solved.tau(), Math.abs(exact.tau()) * 1e-12,
                    "stress differs from the closed form at step " + step);
            assertEquals(exact.epsP[0], solved.epsP[0], exact.epsP[0] * 1e-12 + 1e-300,
                    "plastic strain differs from the closed form at step " + step);
        }
        assertTrue(solved.epsP[0] > 0.001, "the path must actually have gone plastic");
    }

    @Test
    @DisplayName("the solved state lands exactly on the rate- and temperature-dependent surface")
    void consistencyHolds() {
        // The return map's own consistency condition, read back out afterwards. With the full
        // law this is the only check that the iteration converged to the right root rather
        // than to a nearby one: the equivalent stress has to equal the flow stress evaluated
        // at the strain, rate and temperature the step itself produced.
        JohnsonCook jc = JohnsonCook.STEEL_4340;
        Point p = Point.johnsonCook(jc);
        final double dt = 2.0e-8;
        final double dGam = 5.0e-4;

        double previousEps = 0.0;
        for (int step = 0; step < 600; step++) {
            // The temperature the solve will have used: the one standing at the start of the
            // step, before this increment's plastic work is added. Reading it afterwards
            // instead compares the answer against a surface that moved after it was computed,
            // which is how the explicit temperature treatment looks like an error.
            double temperature = jc.temperature(p.work[0], DENSITY);

            double dGamma = p.shear(dGam, dt);
            if (dGamma == 0.0) continue;

            double dEps = p.epsP[0] - previousEps;
            previousEps = p.epsP[0];

            double rate = dEps / dt;
            double expected = jc.flowStress(p.epsP[0], rate, temperature);

            assertEquals(expected, J2.vonMises(p.sig, 0), expected * 1e-9,
                    "at step " + step + " the equivalent stress is off the surface the law "
                            + "defines at this strain, rate and temperature");
        }
        assertTrue(previousEps > 0.01, "the path must reach appreciable plastic strain");
    }

    @Test
    @DisplayName("the explicit temperature lags by one step, and the lag is negligible at real steps")
    void temperatureLagIsBoundedByTheStep() {
        // The temperature is taken at the start of the step rather than solved for alongside
        // the multiplier, so the stress can end a step fractionally outside the surface
        // evaluated at the *end* of it -- by exactly the amount the point heated while
        // flowing. That is a real approximation and it is worth measuring rather than
        // asserting away, because it is the one place this return map is not exact.
        //
        // It scales with the plastic strain taken in a single step. The reference Taylor case
        // reaches about 0.75 plastic strain over roughly five hundred CFL-limited steps, so
        // an increment near 1.5e-3 is what it actually uses -- not the far smaller number it
        // is tempting to assume from the 40 ns timestep alone.
        JohnsonCook jc = JohnsonCook.STEEL_4340;

        double runLike = overshoot(jc, 2.0e-3);
        double brutal = overshoot(jc, 2.0e-2);

        // A part in ten thousand, against a gate of five parts in a hundred. Small enough to
        // ignore and large enough that it should be stated rather than called zero.
        assertTrue(runLike < 1.0e-3,
                "at the increment the reference case actually takes, the end-of-step surface "
                        + "is violated by " + runLike + " relative. Beyond about 1e-3 the "
                        + "explicit temperature stops being a free approximation and the "
                        + "thermal term belongs inside the solve");
        assertTrue(brutal > runLike * 2.0,
                "the lag must grow with the step taken (" + runLike + " then " + brutal
                        + "), or it is not the lag being measured");
    }

    /** Relative amount by which the stress ends a step outside the end-of-step surface. */
    private static double overshoot(JohnsonCook jc, double dGam) {
        Point p = Point.johnsonCook(jc);
        final double dt = 1.0e-8;
        double worst = 0.0;
        double previousEps = 0.0;
        for (int step = 0; step < 400; step++) {
            p.shear(dGam, dt);
            double dEps = p.epsP[0] - previousEps;
            previousEps = p.epsP[0];
            if (dEps == 0.0) continue;
            double surface = jc.flowStress(p.epsP[0], dEps / dt,
                    jc.temperature(p.work[0], DENSITY));
            worst = Math.max(worst, (J2.vonMises(p.sig, 0) - surface) / surface);
        }
        return worst;
    }

    @Test
    @DisplayName("a faster strain path raises the flow stress")
    void rateHardeningIsObservable() {
        // The whole point of the rate term, expressed as the thing it must do. Identical
        // strain increments, different intervals, so the only difference is the rate.
        JohnsonCook jc = law(SY, 0.0, 0.26, 0.014, 0.0);
        double slowStress = drive(jc, 1.0e-5);
        double fastStress = drive(jc, 1.0e-9);

        assertTrue(fastStress > slowStress * 1.02,
                "the fast path reached " + fastStress / 1e6 + " MPa and the slow one "
                        + slowStress / 1e6 + " MPa; four decades of strain rate must be worth "
                        + "more than this");

        // And quantitatively: four decades of rate at C = 0.014 is 1 + 0.014 ln(1e4).
        assertEquals(1.0 + 0.014 * Math.log(1.0e4), fastStress / slowStress, 0.02,
                "the rate hardening is not the logarithm the law specifies");
    }

    /**
     * Shears a point well into the plastic regime and returns the flow stress it settles at.
     *
     * <p>The step count is not arbitrary. An increment of 2e-4 engineering shear is worth
     * about 16 MPa elastically, so roughly the first thirty steps are spent merely reaching
     * the yield surface at sy/sqrt(3); anything that stops soon after is comparing two
     * barely-yielded points and will report that rate and temperature do almost nothing.
     */
    private static double drive(JohnsonCook jc, double dt) {
        Point p = Point.johnsonCook(jc);
        for (int step = 0; step < 2000; step++) p.shear(2.0e-4, dt);
        assertTrue(p.epsP[0] > 0.1,
                "the driver only reached " + p.epsP[0] + " plastic strain; it is measuring the "
                        + "approach to yield rather than flow");
        return J2.vonMises(p.sig, 0);
    }

    @Test
    @DisplayName("plastic heating softens the material as it deforms")
    void thermalSofteningIsObservable() {
        // Same law twice, differing only in whether plastic work becomes heat. This isolates
        // the thermal term from everything else, which matters because on the Taylor case it
        // is very nearly cancelled by the rate term and so is otherwise hard to see at all.
        double adiabatic = drive(law(SY, 0.0, 0.26, 0.0, 0.9), 1.0e-8);
        double isothermal = drive(law(SY, 0.0, 0.26, 0.0, 0.0), 1.0e-8);

        assertTrue(adiabatic < isothermal,
                "heating must soften: adiabatic " + adiabatic / 1e6 + " MPa against isothermal "
                        + isothermal / 1e6 + " MPa");
        assertEquals(SY, isothermal, SY * 1e-9,
                "with no heating and no hardening the flow stress must sit exactly at A");
    }

    @Test
    @DisplayName("plastic flow still changes no volume")
    void flowIsIncompressible() {
        // Unchanged by the new solve, and worth re-asserting on this path: the flow direction
        // is deviatoric by construction, so nothing the multiplier does can touch the
        // pressure. The plastic dissipation audit integrates work density against reference
        // volumes and is only correct because of this.
        Point p = Point.johnsonCook(JohnsonCook.STEEL_4340);
        for (int step = 0; step < 50; step++) {
            p.shear(5.0e-4, 1.0e-8);
            assertEquals(0.0, p.pressure(), 1.0,
                    "pure shear produced a pressure of " + p.pressure() + " Pa at step " + step);
        }
        assertTrue(p.epsP[0] > 0.0, "the path must have yielded for this to mean anything");
    }

    @Test
    @DisplayName("the stress never ends a step outside the yield surface")
    void neverOutsideTheSurface() {
        // The one thing a return map exists to guarantee. Under an iterative solve this is
        // also a statement about the solve: a root found to the wrong side, or an iteration
        // that gave up early, shows up here as a state the surface does not contain.
        JohnsonCook jc = JohnsonCook.STEEL_4340;
        Point p = Point.johnsonCook(jc);

        // Deliberately brutal increments -- far larger than any CFL-limited step -- so the
        // elastic trial overshoots the surface by a long way and the solve has real work.
        for (double dGam : new double[]{1e-3, 5e-3, 2e-2, 1e-1, 5e-2, 1e-3}) {
            double before = p.epsP[0];
            // The surface the solve targeted, which is the one at the temperature standing
            // when the step began. Temperature is explicit here by design; see
            // temperatureLagIsBoundedByTheStep for what that costs and why it is acceptable.
            double temperature = jc.temperature(p.work[0], DENSITY);

            p.shear(dGam, 1.0e-8);
            double rate = (p.epsP[0] - before) / 1.0e-8;
            double surface = jc.flowStress(p.epsP[0], rate, temperature);

            assertTrue(J2.vonMises(p.sig, 0) <= surface * (1.0 + 1e-8),
                    "equivalent stress " + J2.vonMises(p.sig, 0) / 1e6 + " MPa exceeds the "
                            + "surface at " + surface / 1e6 + " MPa after a " + dGam
                            + " shear increment");
        }
    }

    @Test
    @DisplayName("plastic work is the flow stress times the plastic strain increment")
    void plasticWorkIsConsistent() {
        // The work accumulator feeds the temperature, which feeds the flow stress, so an
        // error here is self-reinforcing rather than merely wrong. Checked against the
        // independent expression sigma_y * d(eps_p) rather than against the expansion the
        // kernel actually uses.
        JohnsonCook jc = JohnsonCook.STEEL_4340;
        Point p = Point.johnsonCook(jc);
        final double dt = 1.0e-8;

        double expected = 0.0;
        double previousEps = 0.0;
        for (int step = 0; step < 40; step++) {
            p.shear(5.0e-4, dt);
            double dEps = p.epsP[0] - previousEps;
            previousEps = p.epsP[0];
            // Flow stress after the update, which is where the state ends up sitting.
            expected += J2.vonMises(p.sig, 0) * dEps;
        }

        assertEquals(expected, p.work[0], expected * 1e-3,
                "accumulated plastic work " + p.work[0] + " against " + expected
                        + " from the flow stress and the strain increments");
        assertTrue(p.work[0] > 0.0);
    }
}
