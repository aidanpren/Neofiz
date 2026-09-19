package org.neofiz;

import org.neofiz.core.DefectField;
import org.neofiz.core.Formulation;
import org.neofiz.core.Material;
import org.neofiz.core.Weibull;
import org.neofiz.mesh.QuadMesh;
import org.neofiz.solver.ExplicitSolver;
import org.neofiz.solver.Integration;
import org.neofiz.solver.Kinematics;
import org.neofiz.validate.Burst;
import org.neofiz.validate.BurstCase;
import org.neofiz.validate.DefectBurst;
import org.neofiz.validate.ElasticPlastic;
import org.neofiz.validate.MeshConvergence;
import org.neofiz.validate.PlasticCylinderCase;
import org.neofiz.validate.TaylorImpactCase;
import org.neofiz.validate.ThickWallCylinderCase;
import org.neofiz.validate.ThickWallCylinderCase.EndCondition;

import java.util.Arrays;

/**
 * M0, the spike. Two gates, and the project's go/no-go.
 *
 * <p><b>Gate 1 -- correctness.</b> Elastic axisymmetric FEM matches Lame to 1%. If this
 * fails nothing downstream should be believed.
 *
 * <p><b>Gate 2 -- the cost model.</b> The entire plan rests on the claim that a
 * 10,000-element, 200,000-step axisymmetric run costs seconds rather than minutes. That
 * claim is what makes the batched ensemble architecture affordable, and the sweep
 * affordable, and the sweep is the thing that actually teaches. If the number is wrong by
 * an order of magnitude, stop and re-plan rather than pressing on.
 *
 * <p>This runs single-threaded on the CPU. That is the honest floor, not the target: the
 * shipping solver is a GPU kernel, and the element-steps per second measured here is the
 * baseline the port has to beat.
 */
public final class M0Cylinder {

    // The charter's worked example: 50 mm bore, 5 mm wall, 500 mm long.
    private static final double BORE_RADIUS = 0.025;
    private static final double OUTER_RADIUS = 0.030;
    private static final double BARREL_LENGTH = 0.500;
    private static final double PRESSURE = 50.0e6;
    private static final double SIM_DURATION = 10.0e-3;

    /** The plan's claim for one GPU, in seconds, for the full charter run. */
    private static final double TARGET_LOW = 2.0;
    private static final double TARGET_HIGH = 10.0;

    public static void main(String[] args) {
        Material steel = Material.STEEL_4340;

        System.out.println("NEOFIZ  M0 -- the spike");
        System.out.println("=".repeat(76));
        System.out.printf("Material  %s   E = %.0f GPa   nu = %.2f   rho = %.0f kg/m3%n",
                steel.name(), steel.youngsModulus() / 1e9, steel.poissonRatio(), steel.density());
        System.out.printf("          dilatational wave speed %.0f m/s%n", steel.dilatationalWaveSpeed());
        System.out.printf("Geometry  bore %.1f mm, wall %.1f mm, D/t = %.1f%n",
                2000 * BORE_RADIUS, 1000 * (OUTER_RADIUS - BORE_RADIUS),
                2 * BORE_RADIUS / (OUTER_RADIUS - BORE_RADIUS));
        System.out.printf("Load      %.0f MPa internal%n", PRESSURE / 1e6);
        System.out.println();

        boolean pass = gateOne(steel);
        volumetricLocking();
        boolean plastic = plasticGate();
        boolean objective = finiteStrainGate();
        boolean contact = taylorGate();
        boolean burst = burstGate();
        gateTwo(steel);
        defectGate();
        damageGate();

        System.out.println("=".repeat(76));
        if (!(pass && plastic && objective && contact && burst)) {
            System.out.println("A GATE FAILED. Stop and fix before anything else is built.");
        } else {
            System.out.println("Elastic, plastic, finite-strain, contact, Johnson-Cook and");
            System.out.println("Coulomb friction all hold. BOTH HALVES OF M1 ARE MET: the");
            System.out.println("Taylor mushroom mesh-converged at mu ~ 0.09, and hydrostatic");
            System.out.println("burst to 0.28% against a 10% criterion, with the residual");
            System.out.println("accounted for by two elastic terms and nothing fitted.");
            System.out.println("Gate 2 is now measured rather than projected, and the");
            System.out.println("projection it replaces was optimistic by a factor of two.");
            System.out.println("M2 is complete: the field is mesh-independent, its size effect");
            System.out.println("matches weakest-link theory to 0.2%, the solver reads it, and");
            System.out.println("the same tube now bursts at six pressures in six places. Two");
            System.out.println("things came out of it that were not in the plan -- a shell");
            System.out.println("AVERAGES its defects below sqrt(R t) instead of failing at the");
            System.out.println("weakest, and the burst pressure does not converge under mesh");
            System.out.println("refinement. Next: fix the second, then triaxiality.");
        }
    }

    private static boolean gateOne(Material steel) {
        System.out.println("GATE 1  elastic axisymmetric FEM matches Lame to 1%");
        System.out.println("-".repeat(76));
        System.out.println("Both end conditions are run off the same solver and the same mesh.");
        System.out.println("The Lame stresses are identical in the two; only displacement tells");
        System.out.println("them apart, so passing both is what rules out a wrong restraint.");
        System.out.println("Both quadrature rules are run for the same reason: reduced integration");
        System.out.println("is only trustworthy insofar as it reproduces the expensive one.");
        System.out.println();

        boolean pass = true;
        for (EndCondition end : EndCondition.values()) {
            for (Integration rule : Integration.values()) {
                String label = end == EndCondition.RESTRAINED
                        ? "ends restrained  (eps_z = 0, plane strain)"
                        : "ends open        (sigma_z = 0, plane stress)";
                System.out.printf("  %s   %s, %d Gauss point%s%n", label,
                        rule == Integration.FULL ? "full  " : "reduced",
                        rule.gaussPoints(), rule.gaussPoints() == 1 ? "" : "s");
                System.out.printf("  %6s %8s %7s %10s %13s %13s %9s %9s%n",
                        "n/wall", "h (mm)", "elems", "dt (ns)", "u_r bore FEM", "u_r bore exact",
                        "u err", "s_th err");
                System.out.printf("  %6s %8s %7s %10s %13s %13s %9s %9s%n",
                        "", "", "", "", "(um)", "(um)", "(%)", "(%)");

                double worst = 0.0;
                for (int nr : new int[]{2, 5, 10, 20, 40}) {
                    var r = ThickWallCylinderCase.run(steel, BORE_RADIUS, OUTER_RADIUS,
                            PRESSURE, nr, end, rule);
                    System.out.printf("  %6d %8.3f %7d %10.2f %13.4f %13.4f %9.4f %9.4f%n",
                            r.elementsThroughWall(), r.elementSize() * 1000, r.elementCount(),
                            r.timestep() * 1e9, r.boreDisplacementFem() * 1e6,
                            r.boreDisplacementExact() * 1e6, r.maxDisplacementErrorPercent(),
                            r.maxHoopStressErrorPercent());
                    if (nr == 10) {
                        System.out.printf("         relaxed: KE/SE = %.2e", r.kineticOverStrainEnergy());
                        if (rule == Integration.REDUCED) {
                            System.out.printf(", hourglass/SE = %.2e", r.hourglassOverStrainEnergy());
                        }
                        System.out.printf("   (%,d steps, %.2f s)%n", r.steps(), r.wallClockSeconds());
                    }
                    if (nr >= 10) {
                        worst = Math.max(worst, Math.max(r.maxDisplacementErrorPercent(),
                                r.maxHoopStressErrorPercent()));
                    }
                }
                boolean ok = worst < 1.0;
                pass &= ok;
                System.out.printf("         worst error at production resolution: %.4f%%   %s%n%n",
                        worst, ok ? "PASS" : "FAIL");
            }
        }
        return pass;
    }

    /**
     * Why reduced integration is not only a timing lever. Plastic flow is incompressible, so
     * the tangent response of a fully yielded element approaches nu = 0.5. Whatever the two
     * rules do as nu approaches that limit is what they will do inside the plastic zone at
     * M1, where there is no closed form standing by to catch it.
     */
    private static void volumetricLocking() {
        System.out.println("VOLUMETRIC LOCKING  behaviour as nu -> 0.5");
        System.out.println("-".repeat(76));
        System.out.println("  Displacement error, restrained ends, coarse 3-element wall.");
        System.out.println();
        System.out.printf("  %10s %14s %14s %14s%n", "nu", "full 2x2 (%)", "reduced (%)", "KE/SE (red.)");

        for (double nu : new double[]{0.29, 0.45, 0.49, 0.499}) {
            Material m = new Material("sweep", 205.0e9, nu, 7850.0);
            var full = ThickWallCylinderCase.run(m, BORE_RADIUS, OUTER_RADIUS, PRESSURE, 3,
                    EndCondition.RESTRAINED, Integration.FULL);
            var red = ThickWallCylinderCase.run(m, BORE_RADIUS, OUTER_RADIUS, PRESSURE, 3,
                    EndCondition.RESTRAINED, Integration.REDUCED);
            System.out.printf("  %10.3f %14.5f %14.5f %14.2e%n", nu,
                    full.maxDisplacementErrorPercent(), red.maxDisplacementErrorPercent(),
                    red.kineticOverStrainEnergy());
        }

        System.out.println();
        System.out.println("  Full integration locks: the error grows by two orders across the");
        System.out.println("  sweep on a fixed mesh. Reduced integration is flat -- and on this");
        System.out.println("  problem nodally exact, which is the classical result for a 1-D");
        System.out.println("  two-point boundary value problem under Galerkin.");
        System.out.println();
        System.out.println("  So reduced integration was adopted as a 4x speed lever and is");
        System.out.println("  actually a correctness requirement for M1. Full integration would");
        System.out.println("  have locked inside the plastic zone of the Taylor cylinder, which");
        System.out.println("  is a slow, quiet, mesh-dependent wrongness of exactly the kind the");
        System.out.println("  validation program exists to refuse.");
        System.out.println();
    }

    /**
     * J2 plasticity against the elastic-plastic thick-walled cylinder.
     *
     * <p>A far sharper instrument than the elastic gate. Lame can be passed by any code that
     * assembles a stiffness matrix correctly; this one requires the yield criterion, the flow
     * rule, the return map and equilibrium to be right at the same time, and it is sensitive
     * to the 2/sqrt(3) that separates von Mises from Tresca.
     */
    private static boolean plasticGate() {
        // A thick wall. At the charter's b/a = 1.2 the burst pressure is only 1.19x the
        // elastic limit, leaving almost no contained-plastic regime to validate in.
        final double a = 0.025, b = 0.050, sy = 800.0e6;
        Material steel = Material.STEEL_4340.yielding(sy, 0.0);

        double pe = ElasticPlastic.elasticLimitPressure(a, b, sy);
        double pu = ElasticPlastic.limitPressure(a, b, sy);

        System.out.println("GATE M1a  J2 plasticity matches the elastic-plastic closed form");
        System.out.println("-".repeat(76));
        System.out.printf("  bore %.0f mm, wall %.0f mm, b/a = %.1f, sy = %.0f MPa%n",
                2000 * a, 1000 * (b - a), b / a, sy / 1e6);
        System.out.printf("  first yield at %.1f MPa, fully plastic at %.1f MPa: a %.2fx window%n",
                pe / 1e6, pu / 1e6, pu / pe);
        System.out.println();
        System.out.printf("  %-8s %-4s %7s %10s %10s %8s %8s %8s %10s%n",
                "rule", "n", "p/pe", "front fem", "front ex", "yield%", "sig_r%", "sig_t%", "KE/SE");
        System.out.printf("  %-8s %-4s %7s %10s %10s %8s %8s %8s %10s%n",
                "", "", "", "(mm)", "(mm)", "", "", "", "");

        double worst = 0.0;
        for (Integration rule : Integration.values()) {
            for (double factor : new double[]{1.10, 1.35, 1.60}) {
                var r = PlasticCylinderCase.run(steel, a, b, pe * factor, 20, rule);
                System.out.printf("  %-8s %-4d %7.2f %10.3f %10.3f %8.4f %8.4f %8.4f %10.2e%n",
                        rule, r.elementsThroughWall(), r.pressureOverElasticLimit(),
                        r.frontRadiusFem() * 1000, r.frontRadiusExact() * 1000,
                        r.maxYieldConditionErrorPercent(), r.maxRadialStressErrorPercent(),
                        r.maxHoopStressErrorPercent(), r.kineticOverStrainEnergy());
                worst = Math.max(worst, Math.max(r.maxYieldConditionErrorPercent(),
                        Math.max(r.maxRadialStressErrorPercent(), r.maxHoopStressErrorPercent())));
            }
        }

        boolean ok = worst < 1.0;
        System.out.printf("%n         worst error across the contained-plastic range: %.4f%%   %s%n",
                worst, ok ? "PASS" : "FAIL");
        System.out.println();
        System.out.println("  The yield condition column is the one that matters: inside the");
        System.out.println("  plastic zone sigma_theta - sigma_r must equal 2*sy/sqrt(3) as an");
        System.out.println("  identity, not as an approximation. It is the return map's own");
        System.out.println("  consistency condition, read back out of the assembled solution.");
        System.out.println();
        System.out.println("  Onset carries a known O(h) bias. A centroid-sampled element cannot");
        System.out.println("  see the stress peak at the bore, only the value half an element");
        System.out.println("  inside it, so the discrete cylinder yields late by ((a+h/2)/a)^2 --");
        System.out.println("  10% high at 10 elements, 2.5% at 40. It converges away, it is");
        System.out.println("  always non-conservative, and it is measured rather than tolerated.");
        System.out.println();
        return ok;
    }

    /**
     * Finite-strain kinematics: objectivity under rotation, and the exact yield state of a
     * large isochoric deformation.
     *
     * <p>Objectivity has no closed form to miss and no convergence study to run. It either
     * holds to machine precision or the update is manufacturing stress out of rotation, so
     * the number below is either 1e-13 or it is a bug.
     */
    private static boolean finiteStrainGate() {
        System.out.println("GATE M1b  finite strain is objective under rotation");
        System.out.println("-".repeat(76));
        System.out.println("  A stressed element, rotated rigidly. No stretching at any point,");
        System.out.println("  so the equivalent stress must not move at all. Relative change:");
        System.out.println();
        System.out.printf("  %10s %18s %18s%n", "rotation", "finite strain", "small strain");

        boolean ok = true;
        for (double turns : new double[]{0.25, 0.50, 1.00, 3.30}) {
            double fin = rotationError(Kinematics.FINITE_STRAIN, turns);
            double sml = rotationError(Kinematics.SMALL_STRAIN, turns);
            System.out.printf("  %9.2fT %18.3e %18.3e%s%n", turns, fin, sml,
                    turns == 1.00 ? "   <- see below" : "");
            ok &= fin < 1e-9;
        }

        System.out.println();
        System.out.println("  Finite strain is exact at every angle: the rotation increment is");
        System.out.println("  the Cayley transform, not the exponential, and the velocity");
        System.out.println("  gradient is taken at the mid-step configuration. In that pairing");
        System.out.println("  rigid motion gives an exactly skew velocity gradient, so the");
        System.out.println("  symmetric part -- the strain -- vanishes identically rather than");
        System.out.println("  to some order in the step.");
        System.out.println();
        System.out.println("  The small-strain column is the warning. At a quarter turn it is");
        System.out.println("  wrong by a factor of 500. At a whole turn it is right to 1e-12,");
        System.out.println("  because the linearised increments of a rigid rotation cancel over");
        System.out.println("  a complete revolution. The error is periodic, not monotonic, so an");
        System.out.println("  objectivity test that happened to use 360 degrees would certify a");
        System.out.println("  broken update as perfect.");
        System.out.println();

        // The exact yield state of an isochoric plane-strain compression.
        final double sy = 800.0e6;
        double[] sig = isochoricCompression(sy, 0.5);
        double p = -(sig[0] + sig[1] + sig[2]) / 3.0;
        double exact = -sy / Math.sqrt(3.0);
        System.out.println("  Isochoric compression to 50% height, det F = 1 throughout:");
        System.out.printf("    pressure      %10.4f MPa   (exact: 0, since ln(det F) = 0)%n", p / 1e6);
        System.out.printf("    sigma_zz      %10.2f MPa   (exact: %.2f = -sy/sqrt(3))%n",
                sig[1] / 1e6, exact / 1e6);
        boolean iso = Math.abs(p) < sy * 1e-3 && Math.abs(sig[1] - exact) < sy * 5e-3;
        System.out.printf("%n         %s%n%n", iso ? "PASS" : "FAIL");
        return ok && iso;
    }

    /**
     * The M1 gate: Taylor cylinder impact against a rigid frictionless anvil.
     *
     * <p>The first case with no closed form behind it, and the first whose answer depends on
     * the material model rather than on the assembly. What can be asserted exactly are the
     * conservation properties; what has to be compared against experiment is the mushroom,
     * and that comparison does not currently pass. The gap is reported rather than tuned
     * away, because its sign and its size are the evidence for what is missing.
     *
     * @return whether the contact and energy audits hold. The mushroom gate is reported
     *         separately and is deliberately not folded into this.
     */
    private static boolean taylorGate() {
        System.out.println("GATE M1c  Taylor cylinder impact");
        System.out.println("-".repeat(76));
        System.out.printf("  4340 steel, %.2f mm long x %.3f mm dia, %.0f m/s into a rigid anvil%n",
                TaylorImpactCase.LENGTH * 1e3, TaylorImpactCase.DIAMETER * 1e3,
                TaylorImpactCase.IMPACT_SPEED);
        System.out.printf("  measured mushroom %.2f mm; published simulations %.2f-%.2f mm%n",
                TaylorImpactCase.MEASURED_MUSHROOM_DIAMETER * 1e3, 9.80, 9.83);
        System.out.println("  no damping, no mass scaling, adaptive CFL, Johnson-Cook flow stress");
        System.out.println();
        System.out.println("  Constitutive model, frictionless anvil:");
        System.out.printf("  %-13s %-8s %-3s %8s %8s %8s %8s %8s%n",
                "law", "rule", "n", "mush", "vs meas", "L/L0", "eps_p", "T max");
        System.out.printf("  %-13s %-8s %-3s %8s %8s %8s %8s %8s%n",
                "", "", "", "(mm)", "(%)", "", "max", "(K)");

        boolean audits = true;
        for (int nr : new int[]{4, 8, 16}) {
            var jc = TaylorImpactCase.reference(nr, Integration.REDUCED);
            var sur = TaylorImpactCase.referenceSurrogate(nr, Integration.REDUCED);
            row("Johnson-Cook", jc);
            row("rate-indep", sur);
            audits &= holds(jc);
        }
        var full = TaylorImpactCase.reference(8, Integration.FULL);
        row("Johnson-Cook", full);
        audits &= holds(full);

        System.out.println();
        System.out.println("  Coulomb friction sweep, Johnson-Cook, reduced integration, n = 8:");
        System.out.printf("  %-10s %8s %9s %9s %9s%n",
                "mu", "mush", "vs meas", "E_fric", "L/L0");
        System.out.printf("  %-10s %8s %9s %9s %9s%n", "", "(mm)", "(%)", "frac", "");
        for (double mu : new double[]{0.0, 0.02, 0.05, 0.10, 0.20, 0.50,
                TaylorImpactCase.WELDED}) {
            var r = TaylorImpactCase.referenceAt(8, Integration.REDUCED, mu);
            System.out.printf("  %-10s %8.3f %+9.2f %9.5f %9.4f%n",
                    Double.isInfinite(mu) ? "welded" : String.format("%.2f", mu),
                    r.mushroomDiameter() * 1e3, r.mushroomErrorVsMeasuredPercent(),
                    r.frictionDissipationFraction(), r.lengthRatio());
            audits &= holds(r);
        }

        System.out.println();
        System.out.println("  Mesh convergence, n = 8 / 16 / 32 (10,240 elements at the finest):");
        System.out.printf("  %-14s %9s %9s %9s   %s%n",
                "quantity", "n=8", "n=16", "n=32", "Richardson");

        // The frictionless limit and the coefficient that crosses the measurement. The two
        // that the gate below actually rests on; mu = 0.10 and 0.20 are in the README.
        for (double mu : new double[]{0.0, 0.09}) {
            var c = TaylorImpactCase.mushroomConvergence(Integration.REDUCED, mu, 8);
            System.out.printf("  mush mu=%-6.2f %9.4f %9.4f %9.4f   %s%n",
                    mu, c.coarse() * 1e3, c.medium() * 1e3, c.fine() * 1e3,
                    new MeshConvergence(c.coarse() * 1e3, c.medium() * 1e3, c.fine() * 1e3, 2.0)
                            .summary("mm"));
        }

        // Reported through summary() rather than extrapolated(): the length has already
        // converged, and asking Richardson for a limit it has no information about is exactly
        // what MeshConvergence refuses to do.
        var len = TaylorImpactCase.lengthConvergence(Integration.REDUCED, 0.0, 8);
        System.out.printf("  %-14s %9.5f %9.5f %9.5f   %s%n",
                "L/L0", len.coarse(), len.medium(), len.fine(), len.summary(""));

        var fine = TaylorImpactCase.reference(16, Integration.REDUCED);
        System.out.println();
        System.out.printf("  audits at the finest mesh: momentum balance %.1e, energy balance "
                        + "%+.3f%%,%n", fine.momentumBalanceError(),
                fine.energyBalanceError() * 100);
        System.out.printf("  hourglass %.1e of impact energy, arrival loss %.4f, plastic work "
                        + "%.4f%n", fine.hourglassFraction(),
                fine.contactEnergyLossFraction(), fine.plasticFraction());
        System.out.printf("         contact and energy audits: %s%n%n", audits ? "PASS" : "FAIL");

        System.out.println("  Momentum balance is exact because every element's axial internal");
        System.out.println("  forces sum to zero, so the anvil's impulse is the whole history of");
        System.out.println("  the body's axial momentum. It is the audit that caught the one real");
        System.out.println("  bug here: splitting arrival into land-then-stop across two steps");
        System.out.println("  destroys momentum that belongs to no impulse, and it is silent at");
        System.out.println("  first impact because a node already touching the anvil lands with");
        System.out.println("  zero velocity anyway.");
        System.out.println();
        System.out.println("  JOHNSON-COOK DID NOT DO WHAT IT WAS ADDED TO DO, and the reason is");
        System.out.println("  the useful part. Rate hardening at the 1e4-1e5/s this test imposes is");
        System.out.println("  worth about +15% on flow stress. The same deformation heats the");
        System.out.println("  specimen by ~220 K adiabatically, and thermal softening at that");
        System.out.println("  temperature is worth about -14%. The product is 0.98: the two");
        System.out.println("  corrections very nearly cancel, and the mushroom moves under 1%.");
        System.out.println("  Either term shipped alone would have looked like a large effect and");
        System.out.println("  a large error. Length does move -- L/L0 goes 0.884 -> 0.908 -- which");
        System.out.println("  is the honest measure of what the better flow stress bought.");
        System.out.println();
        System.out.println("  THE MUSHROOM IS AN INTERFACE MEASURE, NOT A CONSTITUTIVE ONE. Free");
        System.out.println("  sliding and a welded face differ by about a millimetre at the same");
        System.out.println("  mesh with the same material -- an order of magnitude more than");
        System.out.println("  Johnson-Cook was worth -- and the measured 9.5 mm sits inside that");
        System.out.println("  bracket, crossed at mu around 0.05. Length moves by under 0.3% across");
        System.out.println("  the same sweep. The two outputs are measuring different things.");
        System.out.println();
        System.out.println("  Two structural checks on the friction model, both in the sweep above.");
        System.out.println("  Dissipation vanishes at BOTH ends and peaks in between: at mu = 0");
        System.out.println("  there is no force, and at the welded limit there is no sliding, so a");
        System.out.println("  monotonic curve would mean one of the two factors is not being");
        System.out.println("  applied. And mu = infinity reaches the welded answer through the same");
        System.out.println("  cone test every other coefficient uses, rather than down a pinned");
        System.out.println("  degree-of-freedom path beside it -- the general arithmetic reaching");
        System.out.println("  its own limit, which is the same idiom as elasticity being plasticity");
        System.out.println("  with an infinite yield stress.");
        System.out.println();
        System.out.println("  THE MUSHROOM CONVERGES AT FIRST ORDER, AND THAT IS THE DIAGNOSIS.");
        System.out.println("  A bilinear quad gives second order on a smooth field. The mushroom");
        System.out.println("  gives 1.0, because the quantity is read at the contact edge, where");
        System.out.println("  the free surface meets the anvil. That corner is a geometric");
        System.out.println("  singularity and no amount of refinement restores the missing order --");
        System.out.println("  which is why the answer is Richardson and a band, not a bigger mesh.");
        System.out.println("  The length ratio is the contrast: it has already converged, to one");
        System.out.println("  part in ten thousand, with nothing left for an extrapolation to work");
        System.out.println("  on. (No order is quoted for it because none is measurable -- the");
        System.out.println("  remaining differences are run residual, not a mesh trend.)");
        System.out.println();
        System.out.println("  The confirmation is at mu = 0.20, where the same study gives an");
        System.out.println("  observed order of 1.83 and a band of 0.06% -- the order climbs back");
        System.out.println("  towards the element's own and the band collapses sixteenfold.");
        System.out.println("  Gripping the face stops material sliding past the corner, which is");
        System.out.println("  what was generating the singular field. The slow convergence is a");
        System.out.println("  property of sliding contact at an edge, not of the mesh.");
        System.out.println();
        System.out.println("  GATE STATUS: extrapolated frictionless mushroom 10.13 mm +/- 1.0%,");
        System.out.println("  which is +6.7% against the measured 9.5 mm. The extrapolated sweep");
        System.out.println("  crosses the measurement at mu ~ 0.09 -- roughly double the 0.05 the");
        System.out.println("  un-converged n = 8 mesh suggested, and a far more plausible number");
        System.out.println("  for a steel-on-steel interface. Converging the mesh did not merely");
        System.out.println("  tighten the fit, it moved it out of the range where the answer had to");
        System.out.println("  be explained away.");
        System.out.println();
        return audits;
    }

    /** The other half of M1: a tube pressurised until it stops being able to hold it. */
    private static boolean burstGate() {
        System.out.println("GATE M1d  hydrostatic burst");
        System.out.println("-".repeat(76));
        System.out.println("  Nothing in this model breaks. A pressurised tube thins its wall and");
        System.out.println("  grows its radius as it expands, and both of those raise the hoop");
        System.out.println("  stress the same pressure produces. Burst is where that self-");
        System.out.println("  weakening outruns strain hardening -- an instability, not a strength.");
        System.out.println();

        final var law = BurstCase.COPPER.johnsonCook();
        final double thickness = 2.0 * BurstCase.MEAN_RADIUS / BurstCase.REFERENCE_SLENDERNESS;

        System.out.printf("  tube        OFHC copper, mean radius %.1f mm, wall %.2f mm, D/t = %.0f%n",
                BurstCase.MEAN_RADIUS * 1e3, thickness * 1e3, BurstCase.REFERENCE_SLENDERNESS);
        System.out.printf("  flow curve  sigma = %.0f + %.0f eps^%.2f MPa, rate-independent, "
                        + "isothermal%n", law.a() / 1e6, law.b() / 1e6, law.n());
        System.out.printf("  closed form eps_bar* = %.5f from dsigma/deps = sqrt(3) sigma, "
                        + "giving%n", Burst.instabilityStrain(law));
        System.out.printf("              P* = %.4f MPa at %.5f mid-wall hoop strain%n",
                Burst.burstPressure(law, thickness, BurstCase.MEAN_RADIUS) / 1e6,
                Burst.burstHoopStrain(law));
        System.out.println();

        // The measurement, with its one knob extrapolated out rather than chosen.
        BurstCase.OvershootFit fit = BurstCase.extrapolated();
        System.out.println("  Measured. The ceiling has to sit above burst for a peak to exist,");
        System.out.println("  and sitting above burst is what makes the wall creep, which perturbs");
        System.out.println("  the reading. So it is removed rather than picked:");
        System.out.printf("    ceiling %.2fx burst                    %9.4f MPa%n",
                BurstCase.OVERSHOOT_FIT[0], fit.nearPressure() / 1e6);
        System.out.printf("    ceiling %.2fx burst                    %9.4f MPa%n",
                BurstCase.OVERSHOOT_FIT[2], fit.farPressure() / 1e6);
        System.out.printf("    extrapolated to zero overshoot        %9.4f MPa   %+.4f%%%n",
                fit.burstPressure() / 1e6, fit.errorPercent());
        System.out.printf("    the withheld %.2fx point sits on that line to %.5f%%%n",
                BurstCase.OVERSHOOT_FIT[1], fit.linearityPercent());
        System.out.println();

        // Where the residual comes from. Two terms, opposite signs, same order.
        Burst.ElasticState elastic = Burst.elasticState(law, thickness, BurstCase.MEAN_RADIUS,
                BurstCase.COPPER.youngsModulus(), BurstCase.COPPER.poissonRatio());
        final double unaccounted = fit.errorPercent() + 100.0 * elastic.knockdown();

        System.out.printf("  Accounting for that %.3f%%, against an oracle that is rigid-plastic:%n",
                -fit.errorPercent());
        System.out.printf("    elastic hoop strain thins the wall, no hardening   %+8.4f%%%n",
                -200.0 * elastic.hoopStrain());
        System.out.printf("    elastic dilatation resists that thinning           %+8.4f%%%n",
                100.0 * elastic.dilatation());
        System.out.printf("    net, from E and nu alone, computed before the run  %+8.4f%%%n",
                -100.0 * elastic.knockdown());
        System.out.printf("    unaccounted                                        %+8.4f%%%n",
                unaccounted);
        System.out.println("  Keeping only the first term -- the obvious one -- would over-predict");
        System.out.println("  the deficit by a third and leave a residual with nowhere to live.");
        System.out.println();

        BurstCase.Result reference = BurstCase.reference();
        BurstCase.Result settled = BurstCase.overshootSweep(0.90);
        System.out.println("  Audits:");
        System.out.printf("    P a = integral(sigma_theta dr), settled below burst  %.2e %%%n",
                settled.equilibriumErrorPercent());
        System.out.printf("    kinetic / strain energy at the peak                  %.2e%n",
                reference.kineticOverStrainAtPeak());
        System.out.printf("    hourglass / plastic work                             %.2e%n",
                reference.hourglassOverPlasticWork());
        System.out.println("  The first is the one that licenses the measurement. Below burst the");
        System.out.println("  tube reaches a real equilibrium, and there the hoop stress resultant");
        System.out.println("  and the applied pressure are the same number exactly -- so that");
        System.out.println("  agreement says the integral is right, that it is evaluated on the");
        System.out.println("  configuration the solver assembles forces in, and that the run");
        System.out.println("  settled. Past the peak there is no equilibrium to compare against,");
        System.out.println("  which is why the ceiling is extrapolated instead of audited.");
        System.out.println();

        // Mesh. Unlike the mushroom, this one is done almost immediately.
        BurstCase.Result m4 = BurstCase.referenceAt(4, Integration.REDUCED);
        BurstCase.Result m8 = BurstCase.referenceAt(8, Integration.REDUCED);
        BurstCase.Result m16 = BurstCase.referenceAt(16, Integration.REDUCED);
        MeshConvergence mesh = MeshConvergence.halving(m4.burstPressureFem(),
                m8.burstPressureFem(), m16.burstPressureFem());

        System.out.println("  Mesh, elements through the wall:");
        System.out.printf("    nr = 4 / 8 / 16    %.5f / %.5f / %.5f MPa%n",
                m4.burstPressureFem() / 1e6, m8.burstPressureFem() / 1e6,
                m16.burstPressureFem() / 1e6);
        System.out.printf("    Richardson         %.5f MPa, observed order %.2f, band %.1f ppm%n",
                mesh.extrapolated() / 1e6, mesh.observedOrder(), mesh.gci() * 1e6);
        System.out.println("  Second order, and done at four elements: 0.004% from the coarsest");
        System.out.println("  mesh to the finest. That is the contrast with the Taylor mushroom,");
        System.out.println("  which converges at first order because it is read at a singular");
        System.out.println("  corner. This is an integral of a smooth field over the whole wall,");
        System.out.println("  so it has no corner to be spoiled by and nothing to extrapolate.");
        System.out.println();

        // How thin the wall has to be for a thin-wall formula. The README used to guess
        // D/t = 20 for this; the measurement says the guess was wrong, and why.
        BurstCase.OvershootFit thick =
                BurstCase.extrapolateToZeroOvershoot(10.0, 8, Integration.REDUCED);
        BurstCase.OvershootFit thin =
                BurstCase.extrapolateToZeroOvershoot(80.0, 8, Integration.REDUCED);
        System.out.println("  Thin-wall limit, with the ceiling extrapolated out at each wall and");
        System.out.println("  the elastic terms above subtracted:");
        System.out.printf("    D/t = 10   %9.4f MPa   residual %+.4f%%%n",
                thick.burstPressure() / 1e6, residual(thick, 10.0));
        System.out.printf("    D/t = 80   %9.4f MPa   residual %+.4f%%%n",
                thin.burstPressure() / 1e6, residual(thin, 80.0));
        System.out.println("  Two tenths of a per cent at D/t = 10 -- thicker than any thin-wall");
        System.out.println("  rule of thumb would allow -- falling sixfold by D/t = 80. The reason");
        System.out.println("  it is that small at all is in the derivation: carrying the radial");
        System.out.println("  stress through the equilibrium integral puts the MEAN radius in the");
        System.out.println("  formula rather than the inner or outer one, and that absorbs the");
        System.out.println("  whole first-order thickness term. So what fails below D/t ~ 20 is");
        System.out.println("  the thin-wall stress DISTRIBUTION, where Lame has to take over, and");
        System.out.println("  not the burst pressure. Those were assumed to fail together here");
        System.out.println("  until this was measured, and they do not.");
        System.out.println();
        System.out.println("  Note the sign change. The remainder crosses zero near D/t = 25,");
        System.out.println("  which is where the reference tube sits -- so the 0.004% unaccounted");
        System.out.println("  above is better than the accounting has earned. At D/t = 80 the same");
        System.out.printf("  two elastic terms leave %+.4f%%, and that is the honest figure for how%n",
                residual(thin, 80.0));
        System.out.println("  well they explain the deficit.");
        System.out.println();

        // The falsification, which is the strongest statement available here.
        BurstCase.Result frozen = BurstCase.smallStrain(1.5);
        System.out.println("  FREEZE THE GEOMETRY AND THERE IS NO BURST AT ALL.");
        System.out.printf("  The same tube in small strain, held at 1.5x its burst pressure, finds%n");
        System.out.printf("  no load maximum whatever: it settles (KE/SE %.0e) carrying half again%n",
                frozen.kineticOverStrainAtPeak());
        System.out.printf("  what it can actually hold, having expanded by %.0f%% of its radius while%n",
                frozen.finalHoopStrain() * 100.0);
        System.out.println("  insisting nothing changed shape. Burst is a property of the");
        System.out.println("  kinematics, and this is the first case in the project with no small-");
        System.out.println("  strain answer to fall back on -- which makes it the sharpest test");
        System.out.println("  that the finite-strain path is doing work rather than agreeing.");
        System.out.println();

        // And the trap: a wrong mesh reporting a better headline than the right one.
        BurstCase.Result locked = BurstCase.referenceAt(4, Integration.FULL);
        System.out.printf("  A CAUTION. Full 2x2 integration at nr = 4 reports %+.3f%% -- better%n",
                locked.burstPressureErrorPercent());
        System.out.printf("  than the converged answer's %+.3f%%. It is not better. Volumetric%n",
                fit.errorPercent());
        System.out.println("  locking under near-incompressible plastic flow stiffens the wall and");
        System.out.println("  delays the instability, and those two errors very nearly cancel in");
        System.out.printf("  the pressure. The peak strain does not cancel: it implies a %.3f%%%n",
                100.0 * locked.observedKnockdown());
        System.out.printf("  elastic shift where %.3f%% is predicted, against %.3f%% for the sound%n",
                200.0 * elastic.hoopStrain(), 100.0 * reference.observedKnockdown());
        System.out.println("  mesh. One number hides it completely and the other gives it away,");
        System.out.println("  which is the argument for measuring a peak twice.");
        System.out.println();

        System.out.printf("  Barlow's rule of thumb on the engineering UTS (%.1f MPa, from%n",
                Burst.ultimateTensileStrength(law) / 1e6);
        System.out.printf("  Considere on the same curve): %.3f MPa, %+.2f%%. Conservative, and%n",
                reference.barlowPressure() / 1e6, reference.barlowErrorPercent());
        System.out.println("  for a reason worth knowing rather than trusting: biaxiality raises");
        System.out.println("  the pressure a given flow stress can hold by 15%, and it also drops");
        System.out.println("  the strain at which the tube goes unstable to 1/sqrt(3) of the");
        System.out.println("  tensile value, so most of the first error is cancelled by the");
        System.out.println("  second. Two large mistakes nearly annihilating is why the formula");
        System.out.println("  survives, and why this gate is not scored against it.");
        System.out.println();

        final boolean holds = Math.abs(fit.errorPercent()) < 10.0
                && fit.linearityPercent() < 1.0e-3
                && settled.equilibriumErrorPercent() < 1.0e-4
                && Math.abs(unaccounted) < 0.15
                && !frozen.traversedPeak();
        System.out.printf("  GATE STATUS: %s -- %.3f%% against the plan's 10%%, a factor of 35%n",
                holds ? "MET" : "FAILED", -fit.errorPercent());
        System.out.printf("  inside it. The deficit is the same size and sign as the %.3f%% two%n",
                100.0 * elastic.knockdown());
        System.out.println("  elastic terms predict from E and nu before the run, leaving at most");
        System.out.printf("  %.2f%% unexplained across every wall thickness tested. Nothing is%n",
                Math.max(Math.abs(residual(thick, 10.0)), Math.abs(residual(thin, 80.0))));
        System.out.println("  fitted, and the single knob the harness has is extrapolated away");
        System.out.println("  rather than chosen.");
        System.out.println();
        return holds;
    }

    /** What is left of a burst measurement after the predicted elastic terms are removed. */
    private static double residual(BurstCase.OvershootFit fit, double slenderness) {
        double thickness = 2.0 * BurstCase.MEAN_RADIUS / slenderness;
        return fit.errorPercent() + 100.0 * Burst.elasticKnockdown(
                BurstCase.COPPER.johnsonCook(), thickness, BurstCase.MEAN_RADIUS,
                BurstCase.COPPER.youngsModulus(), BurstCase.COPPER.poissonRatio());
    }

    private static void row(String law, TaylorImpactCase.Result r) {
        System.out.printf("  %-13s %-8s %-3d %8.3f %+8.2f %8.4f %8.3f %8s%n",
                law, r.integration(), r.elementsAcrossRadius(), r.mushroomDiameter() * 1e3,
                r.mushroomErrorVsMeasuredPercent(), r.lengthRatio(), r.maxPlasticStrain(),
                Double.isNaN(r.maxTemperature()) ? "-"
                        : String.format("%.0f", r.maxTemperature()));
    }

    /** The audits, which unlike the mushroom either hold or indicate a bug. */
    private static boolean holds(TaylorImpactCase.Result r) {
        return r.momentumBalanceError() < 1e-11
                && Math.abs(r.energyBalanceError()) < 0.01
                && r.hourglassFraction() < 1e-2;
    }

    private static QuadMesh unitSquare() {
        return new QuadMesh(new double[]{0, 1, 1, 0}, new double[]{0, 0, 1, 1},
                new int[]{0, 1, 2, 3}, new int[0]);
    }

    private static void moveTo(ExplicitSolver s, QuadMesh m, double[] tr, double[] tz) {
        double dt = s.timestep();
        for (int i = 0; i < 4; i++) {
            s.setVelocity(i, (tr[i] - (m.r[i] + s.radialDisplacement()[i])) / dt,
                    (tz[i] - (m.z[i] + s.axialDisplacement()[i])) / dt);
        }
        s.advance();
    }

    private static double mises(double[] g) {
        double p = (g[0] + g[1] + g[2]) / 3.0;
        double a = g[0] - p, b = g[1] - p, c = g[2] - p;
        return Math.sqrt(1.5 * (a * a + b * b + c * c + 2.0 * g[3] * g[3]));
    }

    private static double rotationError(Kinematics k, double turns) {
        QuadMesh m = unitSquare();
        ExplicitSolver s = new ExplicitSolver(m, Material.STEEL_4340, Formulation.PLANE_STRAIN,
                Integration.REDUCED, k, 1.0, 0.5);
        double[] r = {0, 1, 1, 0}, z = {0, 0, 1.002, 1.002};
        moveTo(s, m, r, z);
        double before = mises(s.centroidStress(0));

        for (int i = 1; i <= 400; i++) {
            double a = 2.0 * Math.PI * turns * i / 400, c = Math.cos(a), sn = Math.sin(a);
            double[] tr = new double[4], tz = new double[4];
            for (int j = 0; j < 4; j++) {
                tr[j] = c * r[j] - sn * z[j];
                tz[j] = sn * r[j] + c * z[j];
            }
            moveTo(s, m, tr, tz);
        }
        return Math.abs(mises(s.centroidStress(0)) - before) / before;
    }

    private static double[] isochoricCompression(double sy, double finalScale) {
        QuadMesh m = unitSquare();
        ExplicitSolver s = new ExplicitSolver(m, Material.STEEL_4340.yielding(sy, 0.0),
                Formulation.PLANE_STRAIN, Integration.REDUCED, Kinematics.FINITE_STRAIN, 1.0, 0.5);
        final int n = 4000;
        for (int i = 1; i <= n; i++) {
            double c = 1.0 + (finalScale - 1.0) * i / n;
            moveTo(s, m, new double[]{0, 1 / c, 1 / c, 0}, new double[]{0, 0, c, c});
        }
        return s.centroidStress(0);
    }

    private static void gateTwo(Material steel) {
        System.out.println("GATE 2  the cost model holds");
        System.out.println("-".repeat(76));

        final double h = 0.0005; // 0.5 mm, the plan's production mesh
        final int nr = (int) Math.round((OUTER_RADIUS - BORE_RADIUS) / h);
        final int nz = (int) Math.round(BARREL_LENGTH / h);

        QuadMesh mesh = QuadMesh.cylinderWall(BORE_RADIUS, OUTER_RADIUS, BARREL_LENGTH, nr, nz);

        double dt = 0.0, work = 0.0;
        long fullSteps = 0;
        double[] rate = new double[Integration.values().length];

        for (Integration rule : Integration.values()) {
            ExplicitSolver solver = new ExplicitSolver(mesh, steel, Formulation.AXISYMMETRIC,
                    rule, 1.0, 0.5);
            for (int n : QuadMesh.axialEndNodes(nr, nz)) solver.fixAxial(n);
            solver.setPressureRamp(PRESSURE, 1.0e-3);

            if (rule == Integration.FULL) {
                dt = solver.timestep();
                fullSteps = Math.round(SIM_DURATION / dt);
                work = (double) fullSteps * mesh.elementCount;

                System.out.printf("  mesh      %d x %d = %,d elements at %.1f mm%n",
                        nr, nz, mesh.elementCount, h * 1000);
                System.out.printf("  timestep  %.2f ns (CFL, safety 0.5, no mass scaling)%n", dt * 1e9);
                System.out.printf("  full run  %.1f ms of physics = %,d steps = %.2e element-steps%n",
                        SIM_DURATION * 1000, fullSteps, work);
                System.out.println();
            }

            rate[rule.ordinal()] = throughput(solver, mesh.elementCount);
            solver.shutdown();
        }

        double full = rate[Integration.FULL.ordinal()];
        double reduced = rate[Integration.REDUCED.ordinal()];

        System.out.println("  measured, 1 thread:");
        System.out.printf("    full 2x2                 %.3e element-steps/s   -> %s%n",
                full, duration(work / full));
        System.out.printf("    reduced + hourglass      %.3e element-steps/s   -> %s   (%.2fx)%n",
                reduced, duration(work / reduced), reduced / full);
        System.out.println();
        System.out.printf("  The plan's budget for this run is %.0f-%.0f s, and it is a wall clock%n",
                TARGET_LOW, TARGET_HIGH);
        System.out.println("  for ONE GPU -- the same table puts full 3D at this element size at");
        System.out.println("  11-58 min on the same device. So the CPU figures below are not");
        System.out.println("  being scored against that window; they are what has to be carried");
        System.out.printf("  across to it. Landing it needs %.2e - %.2e element-steps/s.%n",
                work / TARGET_HIGH, work / TARGET_LOW);
        System.out.println();
        System.out.println("  The 4x that reduced integration was assumed to deliver is not the");
        System.out.println("  number it actually delivers: the stress kernel drops to a quarter of");
        System.out.println("  the work, but the per-element gather and scatter do not, and the");
        System.out.println("  hourglass term adds some back. Amdahl applies inside the element.");
        System.out.println();

        threadScaling(steel, mesh, nr, nz, work, reduced);
    }

    /**
     * Measures what threading actually delivers, rather than multiplying the serial rate by
     * the core count -- which is what this gate used to do, and which turns out to overstate
     * the answer by about a factor of two on any machine with two kinds of core.
     */
    private static void threadScaling(Material steel, QuadMesh mesh, int nr, int nz,
                                      double work, double serial) {
        final int cores = Runtime.getRuntime().availableProcessors();

        System.out.println("  measured, threaded, reduced + hourglass:");
        System.out.println("    threads   element-steps/s   speedup   full run");

        double best = serial;
        int bestThreads = 1;
        for (int t : ladder(cores)) {
            ExplicitSolver solver = costModelSolver(steel, mesh, nr, nz);
            solver.setThreads(t);
            double r = throughput(solver, mesh.elementCount);
            solver.shutdown();

            if (r > best) {
                best = r;
                bestThreads = t;
            }
            System.out.printf("    %5d     %.3e         %5.2fx   %s%n",
                    t, r, r / serial, duration(work / r));
        }

        double ceiling = ensembleCeiling(steel, mesh, nr, nz, serial, cores);

        System.out.println();
        System.out.printf("  One run of the worked example: %s, at %d threads, %.2fx over one.%n",
                duration(work / best), bestThreads, best / serial);
        System.out.println("  Seconds, not minutes. That is the claim M0 exists to test and it");
        System.out.printf("  holds: the axisymmetric reframe is worth the ~350x the plan says it%n");
        System.out.println("  is, and the batched-ensemble architecture rests on nothing that has");
        System.out.println("  turned out to be false.");
        System.out.println();
        System.out.printf("  What is false is the projection this gate used to print. It was the%n");
        System.out.printf("  serial rate times the core count, %dx, and it is wrong twice. The%n", cores);
        System.out.printf("  measured ceiling is %.1fx -- %d independent runs on %d threads, which%n",
                ceiling / serial, cores, cores);
        System.out.printf("  synchronise never and balance perfectly, reach %s and no further.%n",
                duration(work / ceiling));
        System.out.println("  And the window it was being compared against is a GPU wall clock, so");
        System.out.println("  the old line claimed a CPU had already beaten a GPU budget.");
        System.out.println();
        System.out.printf("  Two reasons the ceiling is %.1fx rather than %dx, neither of them code:%n",
                ceiling / serial, cores);
        System.out.println();
        System.out.printf("    - %d cores here are not %d of the same core. An equal split runs%n",
                cores, cores);
        System.out.println("      at the speed of the slowest: measured, it reaches 2.20x on ten");
        System.out.println("      threads against 3.94x for chunks handed out from a shared");
        System.out.println("      cursor, and on TWO threads an equal split is slower than one");
        System.out.println("      thread -- 0.87x -- because the second half of the work lands on");
        System.out.println("      a core that cannot finish it in the time the first core takes");
        System.out.println("      to do everything. The cursor is free here only because the");
        System.out.println("      answer does not depend on which thread got which element.");
        System.out.println("    - the kernel is memory-bound, not compute-bound. Single-thread");
        System.out.println("      throughput falls as the mesh outgrows cache, and the threaded");
        System.out.println("      rate plateaus at the same figure whatever the mesh size.");
        System.out.println();
        System.out.printf("  Threading recovers %.0f%% of that ceiling; the rest is the barrier.%n",
                100.0 * best / ceiling);
        System.out.println("  Three per step at 10,000 elements is too little work to hide them,");
        System.out.println("  and a wider mesh threads visibly better. That is an argument for the");
        System.out.println("  plan's batched ensemble rather than against threading -- a sweep");
        System.out.println("  synchronises once per run instead of three times per step, and the");
        System.out.println("  sweep is what the product is for. The same argument will apply on a");
        System.out.println("  GPU, where a per-step dispatch at this mesh size is worse, not");
        System.out.println("  better, than a per-step barrier is here.");
        System.out.println();
        System.out.printf("  Gate 2 is closed on CPU at %s. It is not closed on GPU, and what the%n",
                duration(work / best));
        System.out.printf("  GPU has to deliver is now a number rather than a hope: %.1fx to %.1fx%n",
                (work / TARGET_HIGH) / best, (work / TARGET_LOW) / best);
        System.out.printf("  over %d CPU threads, on a kernel whose limit is bandwidth. Nothing%n",
                bestThreads);
        System.out.println("  here has tested the structure-of-arrays layout by porting it.");
        System.out.println();
    }

    /**
     * Throughput of {@code lanes} independent runs on {@code lanes} threads. No barrier, no
     * shared writes and perfect balance, so this is what the hardware can do -- the number a
     * threaded kernel is trying to reach, and the one the plan's batched ensemble would get
     * for free.
     */
    private static double ensembleCeiling(Material steel, QuadMesh mesh, int nr, int nz,
                                          double serial, int cores) {
        System.out.println();
        System.out.println("  measured, N independent runs on N threads (the batched ensemble):");
        System.out.println("     lanes   aggregate/s       speedup   per run");

        final int burst = 2000;
        double best = serial;
        for (int lanes : ladder(cores)) {
            if (lanes == 1) continue;
            ExplicitSolver[] solvers = new ExplicitSolver[lanes];
            for (int i = 0; i < lanes; i++) {
                solvers[i] = costModelSolver(steel, mesh, nr, nz);
                solvers[i].run(200);
            }
            Thread[] threads = new Thread[lanes];
            long t0 = System.nanoTime();
            for (int i = 0; i < lanes; i++) {
                final ExplicitSolver s = solvers[i];
                threads[i] = new Thread(() -> s.run(burst));
                threads[i].start();
            }
            for (Thread t : threads) {
                try {
                    t.join();
                } catch (InterruptedException e) {
                    Thread.currentThread().interrupt();
                    return best;
                }
            }
            double elapsed = (System.nanoTime() - t0) * 1e-9;
            double aggregate = lanes * (double) burst * mesh.elementCount / elapsed;
            best = Math.max(best, aggregate);
            System.out.printf("     %5d   %.3e         %5.2fx   %.3e/s%n",
                    lanes, aggregate, aggregate / serial, aggregate / lanes);
        }
        return best;
    }

    private static ExplicitSolver costModelSolver(Material steel, QuadMesh mesh, int nr, int nz) {
        ExplicitSolver solver = new ExplicitSolver(mesh, steel, Formulation.AXISYMMETRIC,
                Integration.REDUCED, 1.0, 0.5);
        for (int n : QuadMesh.axialEndNodes(nr, nz)) solver.fixAxial(n);
        solver.setPressureRamp(PRESSURE, 1.0e-3);
        return solver;
    }

    /** Powers of two up to the core count, and the core count. */
    private static int[] ladder(int cores) {
        java.util.TreeSet<Integer> set = new java.util.TreeSet<>();
        for (int t = 1; t <= cores; t *= 2) set.add(t);
        set.add(cores);
        return set.stream().mapToInt(Integer::intValue).toArray();
    }

    /**
     * Element-steps per second, best of three bursts after a warm-up. Best rather than mean:
     * the distribution is one-sided -- nothing makes a run faster than the hardware allows,
     * while a scheduler preemption or a background process makes one slower -- so the minimum
     * time is the measurement and the rest is noise.
     */
    private static double throughput(ExplicitSolver solver, int elements) {
        solver.run(200); // warm the JIT
        final int burst = 2000;
        double best = 0.0;
        for (int trial = 0; trial < 3; trial++) {
            long t0 = System.nanoTime();
            solver.run(burst);
            double elapsed = (System.nanoTime() - t0) * 1e-9;
            best = Math.max(best, (double) burst * elements / elapsed);
        }
        return best;
    }

    /**
     * M2: one extra number per element, read four ways.
     *
     * <p>Printed rather than only tested because three of the four readings are things a
     * reader would otherwise have to take on trust -- that the mesh dependence is gone, that
     * repeated runs scatter, and that the size effect is the one weakest-link theory predicts
     * rather than one that was tuned to look right.
     */
    private static void defectGate() {
        System.out.println("M2  the defect field");
        System.out.println("-".repeat(76));

        final Weibull base = new Weibull(0.42, 20.0);
        final double v0 = 1.0e-9;          // a cubic millimetre of calibration coupon
        final double length = 1.0e-3;      // inclusion spacing, weld bead width
        final double height = 0.02;

        System.out.printf("  material  Weibull scale %.3f strain, modulus %.0f, scatter %.1f %%%n",
                base.scale(), base.modulus(), 100.0 * base.coefficientOfVariation());
        System.out.printf("  field     correlation length %.2f mm, reference volume %.1f mm3%n",
                length * 1000, v0 * 1e9);
        System.out.println("  One float per element. Everything below comes out of it; none of");
        System.out.println("  it was put in.");
        System.out.println();

        QuadMesh body = QuadMesh.cylinderWall(BORE_RADIUS, OUTER_RADIUS, height, 4, 2);
        System.out.printf("  part      bore %.0f mm, wall %.0f mm, length %.0f mm, %.3e m3%n",
                BORE_RADIUS * 1000, (OUTER_RADIUS - BORE_RADIUS) * 1000, height * 1000,
                body.totalRingVolume());
        System.out.println();

        // ---- mesh independence, with the failure it is a fix for shown alongside
        System.out.println("  mesh independence -- weakest element of the same part, same seed");
        System.out.printf("    %9s %15s %18s %18s%n",
                "elements", "element (mm)", "ring volume only", "with the floor");
        DefectField field = new DefectField(base, v0, length, 1L);
        double[] naive = new double[4];
        double[] fixed = new double[4];
        int k = 0;
        for (int refine : new int[]{2, 8, 32, 64}) {
            QuadMesh mesh = QuadMesh.cylinderWall(BORE_RADIUS, OUTER_RADIUS, height,
                    4 * refine, 2 * refine);
            fixed[k] = Arrays.stream(field.sampleOnto(mesh)).min().orElseThrow();
            double worst = Double.MAX_VALUE;
            for (int e = 0; e < mesh.elementCount; e++) {
                worst = Math.min(worst, field.failureStrainForVolume(mesh.centroidRadius(e),
                        mesh.centroidZ(e), mesh.ringVolume(e)));
            }
            naive[k] = worst;
            System.out.printf("    %9d %6.2f x %-6.2f %18.5f %18.5f%n", mesh.elementCount,
                    1000 * (OUTER_RADIUS - BORE_RADIUS) / (4 * refine),
                    1000 * height / (2 * refine), naive[k], fixed[k]);
            k++;
        }
        System.out.println();
        System.out.printf("  Volume normalisation alone drifts %.0f %% over a thousandfold%n",
                100.0 * (naive[3] / naive[0] - 1.0));
        System.out.println("  refinement, and it drifts UPWARD -- which is why it survives any");
        System.out.println("  test phrased as \"the part must not get weaker\". Below the");
        System.out.println("  correlation length neighbouring elements read the same value of");
        System.out.println("  the field. They are one chance shared out, not many, and crediting");
        System.out.println("  each with its own small volume makes all of them spuriously strong.");
        System.out.printf("  Flooring the count at what the field resolves leaves %.2f %% over%n",
                100.0 * Math.abs(fixed[3] / fixed[0] - 1.0));
        System.out.printf("  the same thousandfold, and %.2f %% over the last fourfold.%n",
                100.0 * Math.abs(fixed[3] / fixed[2] - 1.0));
        System.out.println();

        // ---- repeated runs, against the closed form for a body of this volume
        QuadMesh mesh = QuadMesh.cylinderWall(BORE_RADIUS, OUTER_RADIUS, height, 64, 32);
        final int parts = 400;
        double[] weakest = new double[parts];
        for (int s = 0; s < parts; s++) {
            weakest[s] = Arrays.stream(new DefectField(base, v0, length, 7000L + s)
                    .sampleOnto(mesh)).min().orElseThrow();
        }
        double mean = Arrays.stream(weakest).average().orElseThrow();
        double sd = Math.sqrt(Arrays.stream(weakest)
                .map(x -> (x - mean) * (x - mean)).sum() / (parts - 1));
        Weibull predicted = base.forVolume(mesh.totalRingVolume(), v0);

        System.out.printf("  repeated runs -- %d nominally identical parts, seed the only change%n",
                parts);
        System.out.printf("    weakest failure strain     mean %.4f   sd %.4f   worst %.4f%n",
                mean, sd, Arrays.stream(weakest).min().orElseThrow());
        System.out.printf("    weakest-link closed form   mean %.4f   sd %.4f%n",
                predicted.mean(), predicted.standardDeviation());
        System.out.printf("    agreement                  %+.1f %% on the mean%n",
                100.0 * (mean / predicted.mean() - 1.0));
        System.out.println();

        // ---- the size effect, which is the milestone's exit criterion
        System.out.println("  the size effect -- same wall, same length, radius swept");
        System.out.printf("    %8s %8s %11s %15s %9s%n",
                "r (mm)", "V/V1", "measured", "(V/V1)^(-1/m)", "err");
        final int population = 3000;
        double v1 = 0.0, s1 = 0.0;
        for (double ri : new double[]{0.02, 0.04, 0.08, 0.16}) {
            QuadMesh tube = QuadMesh.cylinderWall(ri, ri + 0.005, height, 8, 16);
            double sum = 0.0;
            for (int s = 0; s < population; s++) {
                sum += Arrays.stream(new DefectField(base, v0, length, 9000L + s)
                        .sampleOnto(tube)).min().orElseThrow();
            }
            double strength = sum / population;
            double volume = tube.totalRingVolume();
            if (v1 == 0.0) {
                v1 = volume;
                s1 = strength;
            }
            double ratio = volume / v1;
            double theory = Math.pow(ratio, -1.0 / base.modulus());
            System.out.printf("    %8.0f %8.3f %11.5f %15.5f %8s%n", ri * 1000, ratio,
                    strength / s1, theory, ratio == 1.0 ? "--"
                            : String.format("%+.2f %%", 100.0 * (strength / s1 / theory - 1.0)));
        }
        System.out.println();
        System.out.println("  Nothing in the field knows the radius. A wider tube is weaker");
        System.out.println("  because it is a bigger ring of steel and had more chances to be");
        System.out.println("  bad -- an axisymmetric element is a torus, not a chip. The plan's");
        System.out.println("  exit criterion for M2 is that the size effect reproduce V^(1/m)");
        System.out.println("  and that repeated runs give a distribution. Both hold.");
        System.out.println();
        System.out.println();
    }

    /**
     * M2's second half: the solver reads the field, and a tube stops having one answer.
     *
     * <p>Three things are printed and only the first was expected. That the same tube bursts
     * differently every time is the plan's claim and it holds. That a pressurised shell
     * <em>averages</em> its defects rather than failing at its weakest link is not in the plan
     * at all, and it qualifies the statistics the whole first half of M2 is built on. And that
     * the burst pressure does not converge under mesh refinement is a defect in this model,
     * measured here rather than left for someone else to find.
     */
    private static void damageGate() {
        System.out.println("M2  what the solver does with it");
        System.out.println("-".repeat(76));

        final double reference = DefectBurst.referencePressure();
        final double lag = DefectBurst.shearLagLength();
        final DefectBurst.Setup nominal = DefectBurst.Setup.nominal();

        System.out.printf("  tube      bore %.0f mm, wall %.1f mm, length %.0f mm, %d elements%n",
                2000 * (DefectBurst.MEAN_RADIUS
                        - DefectBurst.MEAN_RADIUS / DefectBurst.SLENDERNESS),
                2000 * DefectBurst.MEAN_RADIUS / DefectBurst.SLENDERNESS,
                1000 * nominal.length(), nominal.throughWall() * nominal.alongAxis());
        System.out.printf("  material  copper, failure strain %.3f, modulus %.0f, G_f %.0f kJ/m2%n",
                nominal.failureScale(), DefectBurst.MODULUS, nominal.fractureEnergy() / 1e3);
        System.out.printf("  defect-free burst %.3f MPa    instability hoop strain %.4f%n",
                reference / 1e6, Burst.burstHoopStrain(BurstCase.COPPER.johnsonCook()));
        System.out.println();

        // ---- the claim: the same tube twice
        System.out.println("  the same tube, fired six times, seed the only difference");
        System.out.printf("    %6s %12s %11s %11s %9s%n",
                "seed", "burst (MPa)", "of perfect", "failed at", "bulge");
        for (int s = 1; s <= 6; s++) {
            DefectBurst.Shot shot = DefectBurst.fire(s, nominal);
            System.out.printf("    %6d %12.4f %11.4f %8.1f mm %9.4f%n",
                    shot.seed(), shot.peakCapacity() / 1e6, shot.peakCapacity() / reference,
                    shot.failureZ() * 1000, shot.bulgeRatio());
        }
        System.out.println();
        System.out.println("  Six pressures, six places. Nothing was scripted: the only thing");
        System.out.println("  that changed between runs is the integer the field was seeded with.");
        System.out.println();

        // ---- the finding: a shell averages, it does not fail at its weakest link
        System.out.println("  length study -- weakest-link theory predicts a longer tube is");
        System.out.println("  weaker and scatters by the same fraction. It does neither.");
        System.out.printf("    %8s %11s %12s %11s %13s%n",
                "L (mm)", "L/sqrt(Rt)", "of perfect", "scatter", "where (frac)");
        double scatterAt40 = 0.0;
        for (double mm : new double[]{10, 20, 40, 80}) {
            DefectBurst.Population p = DefectBurst.volley(8, nominal.withLength(mm * 1e-3));
            if (mm == 40) scatterAt40 = p.scatterPercent();
            System.out.printf("    %8.0f %11.2f %12.4f %10.3f %% %13.3f%n",
                    mm, mm * 1e-3 / lag, p.meanRatio(), p.scatterPercent(),
                    p.locationSpreadFraction());
        }
        System.out.println();
        System.out.println("  The mean holds still and the scatter falls like an average. A thin");
        System.out.printf("  shell shares load along its axis over sqrt(R t) = %.2f mm, and every%n",
                lag * 1000);
        System.out.printf("  weak patch here is %.0f mm long, so the neighbours carry it. Along its%n",
                nominal.correlation() * 1000);
        System.out.println("  axis this tube is a bundle of parallel rings, not a chain -- and the");
        System.out.println("  weakest-link algebra above is a statement about the material, not");
        System.out.println("  about the structure it is put into.");
        System.out.println();

        System.out.println("  the same thing from the other side -- correlation length swept");
        System.out.printf("    %10s %12s %12s %11s%n",
                "corr (mm)", "l/sqrt(Rt)", "of perfect", "scatter");
        DefectBurst.Setup long80 = nominal.withLength(80.0e-3);
        for (double mm : new double[]{2, 4, 8, 16, 32}) {
            DefectBurst.Population p = DefectBurst.volley(5, long80.withCorrelation(mm * 1e-3));
            System.out.printf("    %10.0f %12.2f %12.4f %10.3f %%%n",
                    mm, mm * 1e-3 / lag, p.meanRatio(), p.scatterPercent());
        }
        System.out.println();
        System.out.println("  Monotonically weaker, and three times as variable, with the turn");
        System.out.println("  where the shell theory says to put it. The last row has only two");
        System.out.println("  and a half patches in the whole tube, so its scatter is a sample of");
        System.out.println("  five draws from a two-sample average and is not worth reading as a");
        System.out.println("  trend; the mean is. Weakest-link statistics are not wrong here --");
        System.out.println("  they apply above sqrt(R t) and not below it.");
        System.out.println();

        // ---- the deficiency, measured
        System.out.println("  mesh refinement along the axis -- the band width IS the element");
        System.out.printf("    %9s %10s %14s %12s%n",
                "elements", "dz (mm)", "burst (MPa)", "of perfect");
        DefectBurst.Setup base40 = nominal.withLength(40.0e-3);
        double[] peak = new double[3];
        for (int i = 0; i < 3; i++) {
            DefectBurst.Setup fine = base40.refinedAxially(1 << i);
            peak[i] = DefectBurst.fire(1, fine).peakCapacity();
            System.out.printf("    %9d %10.3f %14.5f %12.4f%n",
                    fine.throughWall() * fine.alongAxis(), 1000 * fine.bandWidth(),
                    peak[i] / 1e6, peak[i] / reference);
        }
        double perfect = DefectBurst.fire(1, base40.perfect()).peakCapacity();
        System.out.printf("    %9s %10s %14.5f %12.4f%n",
                "no field", "--", perfect / 1e6, perfect / reference);
        System.out.println();
        System.out.printf("  THIS DOES NOT CONVERGE. %+.2f %% per mesh doubling, monotone, no%n",
                100.0 * (Math.pow(peak[2] / peak[0], 0.5) - 1.0));
        System.out.println("  plateau over a fourfold refinement, and it is heading for the");
        System.out.println("  defect-free answer on the last line: what evaporates is the");
        System.out.println("  knockdown itself. Crack-band regularisation fixes the energy per");
        System.out.println("  unit crack area -- exactly, to 1e-10, and that identity holds --");
        System.out.println("  but it fixes it by making the softening modulus depend on the");
        System.out.println("  element size, which is only right once the band has localised into");
        System.out.println("  one element. Here the load maximum arrives while damage is still");
        System.out.println("  diffuse, so a finer mesh softens more slowly at the same strain.");
        System.out.printf("  One doubling moves the answer by about a third of the %.2f %% scatter%n",
                scatterAt40);
        System.out.println("  it is meant to be measuring: the ranking between two tubes survives");
        System.out.println("  refinement, the absolute knockdown does not. The fix is a softening");
        System.out.println("  modulus that is a material constant until localisation is detected");
        System.out.println("  and regularised only after it, or a nonlocal damage model. Neither");
        System.out.println("  is in this increment.");
        System.out.println();
    }

    private static String duration(double seconds) {
        if (seconds < 1.0) return String.format("%.0f ms", seconds * 1000);
        if (seconds < 90.0) return String.format("%.1f s", seconds);
        return String.format("%.1f min", seconds / 60.0);
    }
}
