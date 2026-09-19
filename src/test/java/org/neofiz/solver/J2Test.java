package org.neofiz.solver;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Single-point tests for the J2 return map, driven by prescribed strain paths.
 *
 * <p>A constitutive model is the easiest part of a solver to get subtly wrong and the
 * hardest to catch later, because its errors arrive as plausible stresses rather than as
 * blow-ups. So these assert exact identities -- consistency, incompressibility, the
 * hardening law itself -- to machine precision, rather than comparing against a reference
 * curve to three digits. An identity that holds to 1e-15 is evidence; one that holds to 1e-3
 * is a coincidence waiting to be discovered at M5.
 */
class J2Test {

    private static final double E = 205.0e9;
    private static final double NU = 0.29;
    private static final double SY = 800.0e6;          // 4340 in the right ballpark
    private static final double MU = E / (2.0 * (1.0 + NU));
    private static final double LAMBDA = E * NU / ((1.0 + NU) * (1.0 - 2.0 * NU));

    /** One material point: stress, back stress, equivalent plastic strain. */
    private static final class Point {
        final double[] sig = new double[4];
        final double[] back = new double[4];
        final double[] epsP = new double[1];
        final double[] work = new double[1];
        final double hIso, hKin, sy;
        final J2.Flow flow;

        Point(double sy, double hIso, double hKin) {
            this.sy = sy;
            this.hIso = hIso;
            this.hKin = hKin;
            this.flow = J2.Flow.linear(sy, hIso, hKin);
        }

        /** The interval is irrelevant to a rate-independent law, and asserted so below. */
        private static final double DT = 1.0e-8;

        double shear(double dGam) {
            return J2.update(sig, back, epsP, work, 0, 0, 0, 0, dGam, LAMBDA, MU, flow, DT);
        }

        double strain(double dRR, double dZZ, double dTT, double dGam) {
            return J2.update(sig, back, epsP, work, 0, dRR, dZZ, dTT, dGam,
                    LAMBDA, MU, flow, DT);
        }

        double mises() {
            return J2.vonMises(sig, 0);
        }

        double f() {
            return J2.yieldFunction(sig, back, epsP, 0, sy, hIso);
        }
    }

    @Test
    @DisplayName("shear yields at sy/sqrt(3), the von Mises shear strength")
    void shearYieldsAtTheMisesLimit() {
        Point p = new Point(SY, 0.0, 0.0);
        double step = 1.0e-6;

        double last = 0.0;
        for (int i = 0; i < 20000; i++) {
            if (p.shear(step) > 0.0) break;
            last = p.sig[3];
        }
        assertTrue(last > 0.0, "the path must reach yield");

        // Elastic shear stress just before yield must be within one step of sy/sqrt(3).
        double expected = SY / Math.sqrt(3.0);
        assertEquals(expected, p.sig[3], MU * step * 1.01,
                "von Mises predicts shear yield at sy/sqrt(3), not at sy");
    }

    /**
     * Exact equivalent plastic strain after monotonic shear to total engineering strain
     * gamma, for linear isotropic hardening. Splitting gamma into its elastic and plastic
     * parts and eliminating the stress gives
     *
     * <pre>  gamma = (sy + H epsP) / (sqrt(3) mu) + sqrt(3) epsP</pre>
     *
     * which inverts in closed form. Radial return is exact under proportional loading -- the
     * trial deviator and the returned deviator point the same way -- so the integrated curve
     * should match this to roundoff, not to the step size.
     */
    private static double exactShearPlasticStrain(double gamma, double h) {
        double root3mu = Math.sqrt(3.0) * MU;
        return (gamma - SY / root3mu) / (h / root3mu + Math.sqrt(3.0));
    }

    @Test
    @DisplayName("while flowing, equivalent stress equals the hardening law exactly")
    void flowStressFollowsTheHardeningLaw() {
        // This is the consistency condition itself. If it holds to machine precision at
        // every step of a long plastic path, the return map is landing exactly on the
        // updated surface rather than near it.
        double hIso = 2.0e9;
        Point p = new Point(SY, hIso, 0.0);

        final int steps = 6000;
        final double dGam = 2.0e-6;
        int plasticSteps = 0;
        for (int i = 0; i < steps; i++) {
            double dGamma = p.shear(dGam);
            if (dGamma > 0.0) {
                plasticSteps++;
                assertEquals(SY + hIso * p.epsP[0], p.mises(), (SY + hIso) * 1e-14,
                        "flow stress left the hardening law at step " + i);
                assertEquals(0.0, p.f(), SY * 1e-14,
                        "stress is not on the yield surface at step " + i);
            }
        }
        assertTrue(plasticSteps > 1000, "the path must spend real time flowing");

        // The integrated curve, not merely a spot check that something yielded.
        double expected = exactShearPlasticStrain(steps * dGam, hIso);
        assertEquals(expected, p.epsP[0], expected * 1e-10,
                "accumulated plastic strain does not match the closed-form shear curve");
    }

    @Test
    @DisplayName("plastic flow changes no volume, to machine precision")
    void plasticFlowIsExactlyIncompressible() {
        Point p = new Point(SY, 1.0e9, 1.0e9);

        // A path with a large hydrostatic component, so any leak from the deviatoric return
        // into the pressure has something big to contaminate.
        double pressureBefore = 0.0;
        for (int i = 0; i < 4000; i++) {
            double[] trial = p.sig.clone();
            double tr = -3.0e-6;                       // pure compaction each step
            double lt = LAMBDA * tr;
            double expected = -((trial[0] + lt + 2 * MU * -1.0e-6)
                    + (trial[1] + lt + 2 * MU * -1.0e-6)
                    + (trial[2] + lt + 2 * MU * -1.0e-6)) / 3.0;

            double dGamma = p.strain(-1.0e-6, -1.0e-6, -1.0e-6, 4.0e-6);
            double actual = J2.pressure(p.sig, 0);

            assertEquals(expected, actual, Math.max(1.0, Math.abs(expected)) * 1e-13,
                    "the return map moved the pressure at step " + i + " (dgamma " + dGamma
                            + "). Plastic flow is deviatoric; metals do not change volume when "
                            + "they yield, and a leak here manufactures energy over a long run");
            pressureBefore = actual;
        }
        assertTrue(Math.abs(pressureBefore) > 1.0e8, "the path must build real pressure");
        assertTrue(p.epsP[0] > 0.0, "and must actually yield");
    }

    @Test
    @DisplayName("perfect plasticity caps the flow stress at sy and holds it there")
    void perfectPlasticityCaps() {
        Point p = new Point(SY, 0.0, 0.0);
        final int steps = 20000;
        final double dGam = 2.0e-6;
        for (int i = 0; i < steps; i++) p.shear(dGam);

        assertEquals(SY, p.mises(), SY * 1e-12,
                "with no hardening the equivalent stress must sit exactly on sy forever");

        // With the stress pinned, every further increment of strain is plastic.
        double expected = exactShearPlasticStrain(steps * dGam, 0.0);
        assertEquals(expected, p.epsP[0], expected * 1e-10,
                "all strain beyond yield must go into plastic flow");
    }

    @Test
    @DisplayName("unloading from the yield surface is elastic and follows the elastic slope")
    void unloadingIsElastic() {
        Point p = new Point(SY, 2.0e9, 0.0);
        for (int i = 0; i < 4000; i++) p.shear(2.0e-6);
        assertTrue(p.epsP[0] > 0.0, "must be flowing before the reversal means anything");

        double tauAtReversal = p.sig[3];
        double epsPAtReversal = p.epsP[0];

        double reverse = -2.0e-6;
        double dGamma = p.shear(reverse);

        assertEquals(0.0, dGamma, 0.0,
                "the first reversed step must be elastic; a plastic multiplier here means the "
                        + "model is flowing on unloading, which no metal does");
        assertEquals(tauAtReversal + MU * reverse, p.sig[3], Math.abs(tauAtReversal) * 1e-14,
                "unloading must follow the elastic shear modulus exactly");
        assertEquals(epsPAtReversal, p.epsP[0], 0.0,
                "no plastic strain may accumulate during unloading");
    }

    @Test
    @DisplayName("kinematic hardening reproduces Bauschinger; isotropic hardening cannot")
    void bauschingerSeparatesTheTwoHardeningKinds() {
        double h = 4.0e9;
        Point kinematic = new Point(SY, 0.0, h);
        Point isotropic = new Point(SY, h, 0.0);

        // Identical monotonic path. Under forward loading the two are indistinguishable.
        for (int i = 0; i < 6000; i++) {
            kinematic.shear(2.0e-6);
            isotropic.shear(2.0e-6);
        }
        assertEquals(isotropic.mises(), kinematic.mises(), SY * 1e-10,
                "the two hardening kinds must agree under monotonic loading -- that is exactly "
                        + "why a monotonic test cannot tell them apart");

        double forwardTau = kinematic.sig[3];

        // Now reverse, and find where each one yields again.
        double reverseKin = reverseYieldShear(kinematic);
        double reverseIso = reverseYieldShear(isotropic);

        // Pure kinematic: the surface translated but did not grow, so the elastic span is
        // still the original 2 * sy/sqrt(3) wide.
        double span = forwardTau - reverseKin;
        assertEquals(2.0 * SY / Math.sqrt(3.0), span, span * 1e-3,
                "pure kinematic hardening must leave the elastic range exactly 2*sy wide");

        // Pure isotropic: the surface grew, so reverse yield happens later, not earlier.
        assertTrue(reverseIso < reverseKin,
                "isotropic hardening must predict a later reverse yield than kinematic: "
                        + reverseIso + " vs " + reverseKin + ". A model without a kinematic "
                        + "component over-predicts reverse yield on every cycle, which is the "
                        + "error that matters for anything loaded repeatedly");
    }

    private static double reverseYieldShear(Point p) {
        for (int i = 0; i < 40000; i++) {
            double tauBefore = p.sig[3];
            if (p.shear(-2.0e-6) > 0.0) return tauBefore;
        }
        throw new AssertionError("never reached reverse yield");
    }

    @Test
    @DisplayName("an infinite yield stress is exactly linear elasticity")
    void elasticMaterialTakesTheSamePath() {
        Point p = new Point(Double.POSITIVE_INFINITY, 0.0, 0.0);

        double dRR = 3.0e-4, dZZ = -1.0e-4, dTT = 2.0e-4, dGam = 5.0e-4;
        p.strain(dRR, dZZ, dTT, dGam);

        double lt = LAMBDA * (dRR + dZZ + dTT);
        assertEquals(lt + 2 * MU * dRR, p.sig[0], 1e-6);
        assertEquals(lt + 2 * MU * dZZ, p.sig[1], 1e-6);
        assertEquals(lt + 2 * MU * dTT, p.sig[2], 1e-6);
        assertEquals(MU * dGam, p.sig[3], 1e-6);
        assertEquals(0.0, p.epsP[0], 0.0, "an elastic material must accumulate no plastic strain");
    }

    @Test
    @DisplayName("the stress never ends a step outside the yield surface")
    void stressNeverEscapesTheSurface() {
        // Deliberately brutal: huge steps, so the elastic trial overshoots the surface by a
        // long way. An explicit run takes tiny steps, but a code that only stays admissible
        // for small increments is one step size away from silently drifting outside.
        Point p = new Point(SY, 1.0e9, 1.0e9);
        java.util.Random rng = new java.util.Random(20260918L);

        for (int i = 0; i < 5000; i++) {
            p.strain(1e-3 * rng.nextGaussian(), 1e-3 * rng.nextGaussian(),
                    1e-3 * rng.nextGaussian(), 1e-3 * rng.nextGaussian());
            assertTrue(p.f() <= SY * 1e-12,
                    "stress escaped the yield surface at step " + i + "; f = " + p.f());
            assertTrue(Double.isFinite(p.mises()), "stress went non-finite at step " + i);
        }
        assertTrue(p.epsP[0] > 0.1, "the path must have yielded hard");
    }
}
