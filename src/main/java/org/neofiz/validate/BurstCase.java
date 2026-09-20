package org.neofiz.validate;

import org.neofiz.core.Formulation;
import org.neofiz.core.JohnsonCook;
import org.neofiz.core.Material;
import org.neofiz.mesh.QuadMesh;
import org.neofiz.solver.ExplicitSolver;
import org.neofiz.solver.Integration;
import org.neofiz.solver.Kinematics;

import java.util.ArrayList;
import java.util.List;

/**
 * The second half of the M1 gate: a pressurised tube taken past its burst pressure, compared
 * against {@link Burst}.
 *
 * <p>Same geometry family and same relaxation harness as {@link ThickWallCylinderCase} and
 * {@link PlasticCylinderCase}, and deliberately so -- the elastic and contained-plastic gates
 * keep running on the same mesh, so a regression here can be attributed. What is new is that
 * the answer being measured is a <b>maximum of a curve</b> rather than a point on one, and
 * that changes how the run has to be driven.
 *
 * <h2>You cannot walk up to a load maximum under load control</h2>
 *
 * Burst is where the pressure a tube can hold stops increasing with expansion. Ramp the
 * pressure toward it and there is no difficulty until the last instant, and then no
 * equilibrium at all: past the peak every state the tube can reach carries less than what is
 * being applied, so it accelerates outward and keeps accelerating. Reporting the applied
 * pressure at the moment a run visibly diverges would work, but it is biased high and the bias
 * depends on how fast the ramp was and how sensitive the divergence detector is -- three
 * knobs, one answer.
 *
 * <p>So the pressure is not what gets measured. What gets measured is the <b>wall's
 * capacity</b>, read straight off the stress field through the exact axisymmetric equilibrium
 * identity of {@link #loadCapacity}. That is a function of state, so it is defined on both
 * sides of the peak: as the tube creeps outward past burst, the capacity turns over and falls
 * while the applied pressure sits on its ceiling. The maximum of the capacity is the burst
 * pressure, and it is insensitive to how far above burst the ceiling was set -- which is
 * checkable, and is checked.
 *
 * <p>Two things make the traverse gentle enough to resolve. The ceiling is only a few per cent
 * above the expected burst, so the force imbalance driving the runaway starts near zero; and
 * the relaxation damping that exists to reach static answers also turns the runaway into a
 * creep. <b>The damping is doing the job of a displacement-control device</b>, without any of
 * the machinery, and the audit that says it worked is that the capacity and the applied
 * pressure agree all the way up the stable branch.
 *
 * <h2>Where the equilibrium check can and cannot be applied</h2>
 *
 * The obvious audit for {@link #loadCapacity} is that it reproduce the applied pressure, and
 * the obvious place to apply it is all the way up the stable branch. That is wrong twice over,
 * and both ways of being wrong are instructive.
 *
 * <p><b>During a ramp, a plastically flowing structure is never in equilibrium.</b> It is in a
 * steady creep, and the imbalance driving that creep is whatever the damping needs in order to
 * move the wall as fast as the ramp is loading it -- 20 % of the applied pressure on this case,
 * which is neither an error nor reducible by refinement. <b>And past the peak there is no
 * equilibrium by definition</b>: the applied pressure exceeds what the tube can hold, and the
 * excess is exactly the overshoot, so the same comparison at the peak returns the overshoot
 * and audits nothing.
 *
 * <p>What saves the measurement is that the drag does not have to be modelled out by hand. The
 * damping appears in the radial balance as a body force, so
 * {@code P a = integral(sigma_theta dr) + integral(rho c v r dr)}: the hoop stress resultant is
 * the applied pressure <em>minus</em> the drag, which is to say it is the pressure the wall
 * would hold if it were standing still. That is the quantity wanted, and it is why the peak is
 * measurable on a creeping tube at all. The residual contamination is second order -- the drag
 * also perturbs the through-wall {@code sigma_r} profile, worth {@code O(t/r)} of itself -- and
 * it scales with the overshoot, which is what makes {@link #overshootSweep} the audit that
 * matters rather than a nicety.
 *
 * <p>So {@code equilibriumErrorPercent} is reported from the <em>final</em> sample, where it
 * means something only for a run that did not burst: ramp below burst, let it settle, and the
 * identity has to hold tightly. That is the form the test takes.
 *
 * <h2>What the run can refuse to answer</h2>
 *
 * If the ceiling never exceeds the tube's capacity, there is no peak, no burst, and nothing to
 * report -- so {@link Result#requirePeak} throws rather than returning the last sample and
 * letting it read as an answer. That is not a hypothetical: it is exactly what a
 * <em>small-strain</em> run does, at any ceiling whatever, because with the geometry frozen
 * the self-weakening that causes burst does not exist. See {@link Burst}.
 */
public final class BurstCase {

    // ---------------------------------------------------------------- reference geometry

    /** Mean radius of the reference tube, metres. */
    public static final double MEAN_RADIUS = 25.0e-3;

    /** Mean diameter over wall thickness for the reference tube: thin enough for Burst. */
    public static final double REFERENCE_SLENDERNESS = 25.0;

    /** Elements through the wall at the reference resolution. */
    public static final int ELEMENTS_THROUGH_WALL = 8;

    /**
     * How far above the closed-form burst pressure the applied ceiling is set, for a single
     * run. Two things are asked of it: that it be above the tube's real capacity, so a peak
     * exists to find, and that it be close enough for the traverse to stay quasi-static.
     *
     * <p>It is not, however, a number anyone should have to defend, which is why the gate
     * figure comes from {@link #extrapolateToZeroOvershoot} instead.
     */
    public static final double OVERSHOOT = 1.06;

    /**
     * The three ceilings the extrapolation is fitted through: two to define the line, one to
     * check that it is one.
     *
     * <p>Spread over a factor of four in the perturbation, which is what makes the linearity
     * check worth anything. Tighter than about 1.02 and the post-peak creep no longer clears
     * {@link #PEAK_MARGIN} inside {@link #MAX_PERIODS}; wider than about 1.15 and the traverse
     * stops being quasi-static, with kinetic energy climbing through 4 % of strain energy.
     */
    public static final double[] OVERSHOOT_FIT = {1.03, 1.06, 1.12};

    /** Smoothstep ramp length, in ring-breathing periods. */
    public static final double RAMP_PERIODS = 20.0;

    /** Total run length, in ring-breathing periods. Generous: the post-peak creep is slow. */
    public static final double MAX_PERIODS = 300.0;

    /** Mass-proportional relaxation damping, as a fraction of critical for the breathing mode. */
    public static final double DAMPING = 0.5;

    /** Drop below the running maximum that counts as having traversed the peak. */
    public static final double PEAK_MARGIN = 0.005;

    /** Steps between samples of the capacity curve. */
    public static final int SAMPLE_STEPS = 250;

    private static final double CFL_SAFETY = 0.5;
    private static final int ELEMENTS_ALONG_AXIS = 2;

    /**
     * Annealed OFHC copper on its quasi-static flow curve. Low yield and strong hardening,
     * which is what puts the instability out at 11 % equivalent strain where the geometric
     * term is worth 17 % -- on a high-yield, weakly-hardening steel the same test bursts after
     * 3 % strain and barely exercises the finite-strain path it exists to test.
     */
    public static final Material COPPER =
            Material.COPPER_OFHC.withJohnsonCook(JohnsonCook.COPPER_OFHC.quasiStatic());

    /**
     * One point on the pressure-expansion curve, sampled every {@link #SAMPLE_STEPS} steps.
     *
     * <p>The run has always computed these and kept only the peak. It is the curve, not the
     * peak, that carries the physics: an elastic rise, a knee where the bore surface yields,
     * a plateau where hardening and geometric self-weakening are trading, and a turnover
     * where hardening loses. The peak is one number off the top of that story.
     *
     * <p>{@code capacity} is the wall's own load capacity from {@link #loadCapacity}, which
     * is what "the pressure this tube can hold" means on both sides of the maximum;
     * {@code applied} is the ramp's bore pressure, which tracks it up the stable branch and
     * then sits on its ceiling while the capacity falls away. Plotting both is what makes
     * the instability legible rather than asserted -- the gap between them past the peak is
     * the unbalanced force that is accelerating the wall outward.
     *
     * @param time              solver time, seconds
     * @param hoopStrain        logarithmic mid-wall hoop strain, the closed form's abscissa
     * @param capacity          the wall's load capacity, Pa
     * @param applied           bore pressure the ramp is applying, Pa
     * @param kineticOverStrain kinetic energy over strain energy; the quasi-static audit
     * @param maxPlasticStrain  largest equivalent plastic strain anywhere in the wall
     * @param yieldedFraction   fraction of the wall currently flowing; locates the knee
     */
    public record Sample(double time, double hoopStrain, double capacity, double applied,
                         double kineticOverStrain, double maxPlasticStrain,
                         double yieldedFraction) {
    }

    public record Result(
            Integration integration,
            Kinematics kinematics,
            double meanRadius,
            double thickness,
            double slenderness,
            int elementsThroughWall,
            int elementCount,
            double overshoot,
            boolean traversedPeak,
            double burstPressureFem,
            double burstPressureExact,
            double burstPressureErrorPercent,
            double burstHoopStrainFem,
            double burstHoopStrainExact,
            double appliedCeiling,
            double appliedAtPeak,
            double equilibriumErrorPercent,
            double kineticOverStrainAtPeak,
            double hourglassOverPlasticWork,
            double elasticKnockdownPercent,
            double barlowPressure,
            double barlowErrorPercent,
            double maxPlasticStrain,
            double finalHoopStrain,
            int steps,
            double wallClockSeconds,
            List<Sample> trace) {

        /**
         * The measured burst pressure, or an explanation of why there is not one.
         *
         * <p>A run whose ceiling never exceeded the tube's capacity has no maximum in it. The
         * last sample it took is a pressure the tube <em>held</em>, which is the opposite of
         * the quantity being asked for, and returning it would be worse than returning
         * nothing.
         */
        /**
         * How far the ceiling sat above the wall's capacity at the instant the peak was read,
         * as a fraction. This is the size of the viscous body force the relaxation damping was
         * applying at that moment, and therefore the size of the perturbation contaminating
         * the reading: it is not an error estimate, it is the knob whose irrelevance
         * {@link #overshootSweep} has to demonstrate.
         *
         * <p>It cannot be made small by running longer or refining anything, because past the
         * peak there is by definition no equilibrium -- the applied pressure exceeds what the
         * tube can hold, and that excess <em>is</em> the overshoot. Which is why the honest
         * response is a sweep and not a tolerance.
         */
        public double creepImbalance() {
            return (appliedCeiling - burstPressureFem) / appliedCeiling;
        }

        /**
         * The elastic hoop strain, read off this run's own peak <em>strain</em> and expressed
         * as the pressure knockdown it would cause on its own, as a fraction.
         *
         * <p>Elasticity translates the whole pressure-expansion curve along the strain axis by
         * the elastic hoop strain {@code e}, so the shift in the measured peak strain
         * <b>is</b> {@code e}, obtainable without touching an elastic constant. See
         * {@link Burst#elasticState} for why that translation is exact.
         *
         * <p>What this is <em>not</em> is the pressure knockdown, because the curve is also
         * scaled, by {@code exp(D - 2e)} rather than {@code exp(-2e)}: elastic dilatation
         * works against the thinning and takes back about a quarter of it. So this number
         * deliberately reports only the term the strain can see, and comparing it against
         * {@code 2e} from {@link Burst#elasticState} is a consistency check on two
         * measurements of one quantity -- agreement to about 15 % on the reference tube, which
         * is as well as a mid-wall estimate does on a wall with a strain gradient through it.
         *
         * <p>Loose as that is, it is the check that earns its keep. A coarse fully integrated
         * mesh locks, which stiffens the wall and delays the instability; the two errors very
         * nearly cancel in the reported <em>pressure</em> and do not cancel at all here, so
         * this reads half again too high while the headline looks better than the converged
         * answer. One number agreeing to 15 % catches what the other one hides completely.
         */
        public double observedKnockdown() {
            return 1.0 - Math.exp(-2.0 * (burstHoopStrainFem - burstHoopStrainExact));
        }

        public double requirePeak() {
            if (!traversedPeak) {
                throw new IllegalStateException(String.format(
                        "no burst: the tube held the %.3f MPa ceiling to %.1f%% hoop strain "
                                + "without its capacity ever turning over%s. There is no load "
                                + "maximum in this run to report",
                        appliedCeiling / 1e6, 100.0 * finalHoopStrain,
                        kinematics.isFinite() ? ""
                                : ", which is what small-strain kinematics always does: with "
                                + "the wall thickness and bore radius frozen at their "
                                + "reference values there is no geometric self-weakening and "
                                + "so no instability"));
            }
            return burstPressureFem;
        }

        /**
         * Index into {@link #trace} of the sample the burst pressure was read from, or
         * {@code -1} for a run with no peak in it.
         *
         * <p>Found by scanning rather than recorded during the run, so that it cannot drift
         * out of agreement with the trace it indexes.
         */
        public int peakIndex() {
            if (!traversedPeak) return -1;
            int best = -1;
            double most = Double.NEGATIVE_INFINITY;
            for (int i = 0; i < trace.size(); i++) {
                if (trace.get(i).capacity() > most) {
                    most = trace.get(i).capacity();
                    best = i;
                }
            }
            return best;
        }
    }

    private BurstCase() {
    }

    // ---------------------------------------------------------------- entry points

    /** The M1 burst gate: reference tube, reference mesh, reduced integration. */
    public static Result reference() {
        return referenceAt(ELEMENTS_THROUGH_WALL, Integration.REDUCED);
    }

    /** The reference tube at a chosen resolution and integration rule. */
    public static Result referenceAt(int elementsThroughWall, Integration integration) {
        return run(COPPER, MEAN_RADIUS, REFERENCE_SLENDERNESS, elementsThroughWall,
                integration, Kinematics.FINITE_STRAIN, OVERSHOOT);
    }

    /** The reference mesh on a tube of chosen slenderness, for the thin-wall limit study. */
    public static Result atSlenderness(double slenderness) {
        return run(COPPER, MEAN_RADIUS, slenderness, ELEMENTS_THROUGH_WALL,
                Integration.REDUCED, Kinematics.FINITE_STRAIN, OVERSHOOT);
    }

    /** The reference case at a chosen ceiling, for the insensitivity check. */
    public static Result overshootSweep(double overshoot) {
        return run(COPPER, MEAN_RADIUS, REFERENCE_SLENDERNESS, ELEMENTS_THROUGH_WALL,
                Integration.REDUCED, Kinematics.FINITE_STRAIN, overshoot);
    }

    /**
     * Burst pressure with the harness's own overshoot extrapolated away.
     *
     * @param nearPressure   measured peak at the lowest ceiling
     * @param farPressure    measured peak at the highest ceiling
     * @param burstPressure  the fit's intercept at a ceiling of exactly the burst pressure
     * @param hoopStrain     peak hoop strain, extrapolated the same way
     * @param slope          d(peak) / d(ceiling), Pa per unit of overshoot
     * @param linearityPercent how far the middle measurement sits off the line through the
     *                         outer two, relative to the extrapolated pressure
     * @param exact          the closed-form burst pressure being compared against
     * @param errorPercent   the extrapolated measurement against that closed form
     */
    public record OvershootFit(double nearPressure, double farPressure, double burstPressure,
                               double hoopStrain, double slope, double linearityPercent,
                               double exact, double errorPercent) {
    }

    /**
     * Removes the overshoot from the answer instead of choosing it.
     *
     * <p>The ceiling has to sit above burst for a peak to exist, and sitting above burst is
     * precisely what makes the wall creep, which perturbs the reading. There is no setting
     * that avoids this -- the perturbation goes to zero only in the limit where the traverse
     * takes forever. So it is extrapolated out, and the justification is not the fit quality
     * but the mechanism: the drag is proportional to the overshoot and enters the through-wall
     * {@code sigma_r} profile linearly, so the perturbation is <b>linear in the overshoot with
     * no constant term</b>. It should therefore extrapolate exactly, and
     * {@code linearityPercent} is the check -- measured at 4 parts in 100 000, against a
     * perturbation the fit is removing 20 times larger than that.
     *
     * <p>Fitted through the outer two of {@link #OVERSHOOT_FIT} with the middle one held back
     * as the witness, which is the only arrangement in which the linearity number is evidence
     * rather than a restatement of the fit.
     */
    public static OvershootFit extrapolateToZeroOvershoot(double slenderness,
                                                          int elementsThroughWall,
                                                          Integration integration) {
        final double[] o = OVERSHOOT_FIT;
        final Result near = geometry(slenderness, elementsThroughWall, integration, o[0]);
        final Result mid = geometry(slenderness, elementsThroughWall, integration, o[1]);
        final Result far = geometry(slenderness, elementsThroughWall, integration, o[2]);
        near.requirePeak();
        mid.requirePeak();
        far.requirePeak();

        final double slope = (far.burstPressureFem() - near.burstPressureFem()) / (o[2] - o[0]);
        final double pressure = near.burstPressureFem() + slope * (1.0 - o[0]);
        final double strainSlope =
                (far.burstHoopStrainFem() - near.burstHoopStrainFem()) / (o[2] - o[0]);
        final double strain = near.burstHoopStrainFem() + strainSlope * (1.0 - o[0]);

        final double witness = near.burstPressureFem() + slope * (o[1] - o[0]);
        final double exact = near.burstPressureExact();

        return new OvershootFit(near.burstPressureFem(), far.burstPressureFem(),
                pressure, strain, slope,
                100.0 * Math.abs(mid.burstPressureFem() - witness) / pressure,
                exact, 100.0 * (pressure - exact) / exact);
    }

    /** The M1 gate figure: reference tube and mesh, ceiling extrapolated out. */
    public static OvershootFit extrapolated() {
        return extrapolateToZeroOvershoot(REFERENCE_SLENDERNESS, ELEMENTS_THROUGH_WALL,
                Integration.REDUCED);
    }

    private static Result geometry(double slenderness, int elementsThroughWall,
                                   Integration integration, double overshoot) {
        return run(COPPER, MEAN_RADIUS, slenderness, elementsThroughWall, integration,
                Kinematics.FINITE_STRAIN, overshoot);
    }

    /**
     * The same tube with the geometry frozen. Expected to hold its ceiling indefinitely and
     * to report no peak at all; see {@link Result#requirePeak}.
     */
    public static Result smallStrain(double overshoot) {
        return run(COPPER, MEAN_RADIUS, REFERENCE_SLENDERNESS, ELEMENTS_THROUGH_WALL,
                Integration.REDUCED, Kinematics.SMALL_STRAIN, overshoot);
    }

    // ---------------------------------------------------------------- the run

    public static Result run(Material material, double meanRadius, double slenderness,
                             int elementsThroughWall, Integration integration,
                             Kinematics kinematics, double overshoot) {
        final JohnsonCook law = material.johnsonCook();
        if (law == null) {
            throw new IllegalArgumentException(
                    "burst needs a hardening flow curve; an elastic-perfectly-plastic tube "
                            + "bursts the instant it goes fully plastic and has no instability "
                            + "to find");
        }
        final double thickness = 2.0 * meanRadius / slenderness;
        final double inner = meanRadius - 0.5 * thickness;
        final double outer = meanRadius + 0.5 * thickness;

        final int nr = elementsThroughWall;
        final int nz = ELEMENTS_ALONG_AXIS;
        final double h = thickness / nr;

        QuadMesh mesh = QuadMesh.cylinderWall(inner, outer, h * nz, nr, nz);
        ExplicitSolver solver = new ExplicitSolver(mesh, material, Formulation.AXISYMMETRIC,
                integration, kinematics, 1.0, CFL_SAFETY);

        // eps_z = 0. Not an approximation of the closed form but a restatement of it: the
        // 1 : 1/2 : 0 stress state of a thin pressurised tube has zero deviatoric axial
        // component, so von Mises flow produces no axial strain whether the ends are held or
        // free. Holding them makes the same solution reachable on a two-element slice.
        solver.fixAllAxial();

        final double exact = Burst.burstPressure(law, thickness, meanRadius);
        final double ceiling = overshoot * exact;

        final double omega = material.barWaveSpeed() / meanRadius;
        final double period = 2.0 * Math.PI / omega;
        solver.setPressureRamp(ceiling, RAMP_PERIODS * period);
        solver.setRelaxationDamping(DAMPING, omega);

        final List<Sample> trace = new ArrayList<>();
        boolean traversed = false;
        double peakCapacity = 0.0;
        double peakStrain = 0.0;
        double peakApplied = 0.0;
        double peakKinetic = 0.0;
        double capacity = 0.0;
        double applied = 0.0;
        double strain = 0.0;

        final double deadline = MAX_PERIODS * period;
        final long t0 = System.nanoTime();
        while (solver.time() < deadline) {
            solver.run(SAMPLE_STEPS);

            capacity = loadCapacity(solver, mesh);
            applied = solver.borePressure();
            strain = hoopStrain(solver, mesh, meanRadius);

            // Strain energy is zero on the first sample of a run that has not been loaded
            // yet, and the ratio is reported rather than tested, so it is left at zero there
            // instead of arriving in the trace as an infinity.
            final double strainEnergy = solver.strainEnergy();
            final double kinetic = strainEnergy > 0.0
                    ? solver.kineticEnergy() / strainEnergy : 0.0;
            trace.add(new Sample(solver.time(), strain, capacity, applied, kinetic,
                    solver.maxPlasticStrain(), solver.yieldedFraction()));

            if (capacity > peakCapacity) {
                peakCapacity = capacity;
                peakStrain = strain;
                peakApplied = applied;
                peakKinetic = kinetic;
            } else if (applied >= ceiling && capacity < peakCapacity * (1.0 - PEAK_MARGIN)) {
                // Past the peak, and past the end of the ramp so that a transient dip on the
                // way up cannot be mistaken for the turnover.
                traversed = true;
                break;
            }
        }
        final double wallClock = (System.nanoTime() - t0) * 1e-9;

        final double barlow = Burst.barlow(Burst.ultimateTensileStrength(law),
                thickness, 2.0 * meanRadius);

        // Against plastic dissipation, not strain energy: by burst almost all the work put
        // into the tube has been dissipated, so the residual elastic energy is the wrong and
        // far too flattering denominator. Compare TaylorImpactCase, which had to make the
        // same correction for the same reason.
        final double plastic = solver.plasticDissipation();
        final double hourglass = plastic > 0.0 ? solver.hourglassEnergy() / plastic : 0.0;

        return new Result(integration, kinematics, meanRadius, thickness, slenderness,
                nr, mesh.elementCount, overshoot, traversed,
                peakCapacity, exact, 100.0 * (peakCapacity - exact) / exact,
                peakStrain, Burst.burstHoopStrain(law),
                ceiling, peakApplied,
                100.0 * Math.abs(capacity - applied) / applied, peakKinetic, hourglass,
                100.0 * Burst.elasticKnockdown(law, thickness, meanRadius,
                        material.youngsModulus(), material.poissonRatio()),
                barlow, 100.0 * (barlow - exact) / exact,
                solver.maxPlasticStrain(), strain,
                (int) solver.steps(), wallClock, List.copyOf(trace));
    }

    // ---------------------------------------------------------------- measurement

    /**
     * The pressure the wall is currently able to hold, from the exact axisymmetric equilibrium
     * identity.
     *
     * <p>Radial equilibrium with no body force is
     * {@code d(sigma_r)/dr + (sigma_r - sigma_theta)/r = 0}, and multiplying by {@code r}
     * turns the left side into a perfect derivative:
     *
     * <pre>
     *   d(r sigma_r)/dr = sigma_theta
     *   [r sigma_r] from a to b  =  integral(sigma_theta dr)
     *   P a = integral(sigma_theta dr)                        with sigma_r(b) = 0
     * </pre>
     *
     * <p>No approximation anywhere: not thin-wall, not elastic, not small-strain. It holds in
     * the current configuration with Cauchy stress, which is exactly what the solver stores,
     * so this is a <b>measurement of the same field the solver assembles forces from</b>
     * rather than an independent model of it. Along the stable branch it must reproduce the
     * applied pressure, and the fact that it does is what licenses trusting it past the peak,
     * where there is nothing to compare against.
     *
     * <p>The integral is taken over the whole two-dimensional mesh and divided by the height,
     * which averages over {@code z} instead of picking a row out of the connectivity. The
     * answer is the same -- plane strain makes every row identical -- and it costs nothing,
     * but it keeps this function independent of how {@link QuadMesh#cylinderWall} happens to
     * order its nodes.
     */
    public static double loadCapacity(ExplicitSolver solver, QuadMesh mesh) {
        final boolean deformed = solver.kinematics().isFinite();
        final double[] ur = solver.radialDisplacement();
        final double[] uz = solver.axialDisplacement();

        double hoopArea = 0.0;
        for (int e = 0; e < mesh.elementCount; e++) {
            hoopArea += solver.centroidStress(e)[2] * area(mesh, ur, uz, e, deformed);
        }

        double bore = Double.MAX_VALUE;
        double low = Double.MAX_VALUE;
        double high = -Double.MAX_VALUE;
        for (int i = 0; i < mesh.nodeCount; i++) {
            final double r = deformed ? mesh.r[i] + ur[i] : mesh.r[i];
            final double z = deformed ? mesh.z[i] + uz[i] : mesh.z[i];
            if (r < bore) bore = r;
            if (z < low) low = z;
            if (z > high) high = z;
        }
        return hoopArea / (bore * (high - low));
    }

    /**
     * Logarithmic hoop strain at mid-wall, the strain the closed form is written in.
     *
     * <p>Always from the displacement, in both kinematics. A small-strain run still
     * accumulates the displacement -- it only declines to assemble forces on it -- so this
     * remains the expansion the run produced, and reporting it is how a small-strain run gets
     * to say out loud that it reached 27 % hoop strain while insisting nothing had changed
     * shape.
     */
    public static double hoopStrain(ExplicitSolver solver, QuadMesh mesh, double meanRadius) {
        final double[] ur = solver.radialDisplacement();
        double bore = Double.MAX_VALUE;
        double rim = -Double.MAX_VALUE;
        for (int i = 0; i < mesh.nodeCount; i++) {
            final double r = mesh.r[i] + ur[i];
            if (r < bore) bore = r;
            if (r > rim) rim = r;
        }
        return Math.log(0.5 * (bore + rim) / meanRadius);
    }

    /** Signed area of an element in the (r, z) plane, by the shoelace formula. */
    private static double area(QuadMesh mesh, double[] ur, double[] uz, int element,
                              boolean deformed) {
        final int b = element * 4;
        double sum = 0.0;
        for (int k = 0; k < 4; k++) {
            final int i = mesh.conn[b + k];
            final int j = mesh.conn[b + (k + 1) % 4];
            final double ri = deformed ? mesh.r[i] + ur[i] : mesh.r[i];
            final double rj = deformed ? mesh.r[j] + ur[j] : mesh.r[j];
            final double zi = deformed ? mesh.z[i] + uz[i] : mesh.z[i];
            final double zj = deformed ? mesh.z[j] + uz[j] : mesh.z[j];
            sum += ri * zj - rj * zi;
        }
        return 0.5 * sum;
    }
}
