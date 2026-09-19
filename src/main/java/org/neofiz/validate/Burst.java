package org.neofiz.validate;

import org.neofiz.core.JohnsonCook;

/**
 * Closed-form hydrostatic burst of a thin-walled tube: the pressure at which internal pressure
 * stops being containable. Rigid-plastic, von Mises, thin wall, finite strain.
 *
 * <h2>Burst is an instability, not a strength</h2>
 *
 * Nothing breaks here. There is no failure stress, no critical strain, no criterion of any
 * kind -- and there had better not be, because the charter forbids binary failure. What
 * happens instead is that a pressurised tube is <b>geometrically self-weakening</b>: as it
 * expands the wall thins and the radius grows, and both of those raise the hoop stress that
 * the same pressure produces. Strain hardening pushes back. Burst is the crossing point,
 * where the material stops being able to keep up with its own geometry.
 *
 * <p>That is why this case cannot be run in small strain at all. Freeze the geometry and the
 * feedback disappears: wall thickness and radius are constants, hoop stress is proportional
 * to pressure forever, and the tube quietly carries any pressure the hardening law can
 * eventually match. The instability is a property of the <em>kinematics</em>. It is the first
 * case in this project with no small-strain answer to fall back on, which also makes it the
 * sharpest available test that the finite-strain path is doing real work rather than
 * reproducing the small-strain one to three digits.
 *
 * <h2>Deriving it</h2>
 *
 * A thin closed-end tube under internal pressure carries hoop stress {@code sigma_theta},
 * axial {@code sigma_z = sigma_theta / 2} and radial {@code ~0}. That stress state has zero
 * deviatoric axial component, so von Mises flow puts <b>no axial strain rate</b> anywhere:
 * the tube deforms in plane strain whether or not its ends are held. Two consequences, and
 * both are load-bearing below. Equivalent stress and strain follow from the 1 : 1/2 : 0 state
 *
 * <pre>
 *   sigma_bar = (sqrt(3)/2) sigma_theta        eps_bar = (2/sqrt(3)) eps_theta
 * </pre>
 *
 * <p>and incompressibility with {@code eps_z = 0} gives {@code eps_r = -eps_theta}, so the
 * wall thins at exactly the rate the radius grows:
 *
 * <pre>
 *   t = t0 exp(-eps_theta)      r = r0 exp(+eps_theta)      t/r = (t0/r0) exp(-2 eps_theta)
 * </pre>
 *
 * <p>Equilibrium of a half-ring is {@code P = sigma_theta t / r}, which assembles into the
 * pressure the tube can hold at a given expansion:
 *
 * <pre>
 *   P(eps_bar) = (2/sqrt(3)) (t0/r0) sigma_bar(eps_bar) exp(-sqrt(3) eps_bar)
 * </pre>
 *
 * <p>and burst is that function's maximum:
 *
 * <pre>
 *   dsigma_bar/deps_bar = sqrt(3) sigma_bar
 * </pre>
 *
 * <p><b>Which radius.</b> The exact axisymmetric equilibrium identity is
 * {@code P a = integral(sigma_theta dr)} over the wall, with {@code a} the <em>inner</em>
 * radius -- see {@link BurstCase#loadCapacity}. Keeping the radial stress, which runs from
 * {@code -P} at the bore to zero outside, and using the yield condition
 * {@code sigma_theta - sigma_r = (2/sqrt(3)) sigma_bar} turns that into
 * {@code P (a + t/2) = (2/sqrt(3)) sigma_bar t}. So the radius that belongs in the thin-wall
 * formula is the <b>mean</b> radius, and using it moves the leading error from first order in
 * {@code t/r} to second. Barlow's formula written on the outside diameter is the same
 * expression with the wrong radius, which is most of why it is conservative.
 *
 * <h2>The one memorable number</h2>
 *
 * For a pure power law {@code sigma = K eps^n} the condition collapses to
 *
 * <pre>  eps_bar* = n / sqrt(3)</pre>
 *
 * <p>exactly -- the burst strain is the hardening exponent, scaled by the biaxiality. Compare
 * the uniaxial tensile instability, where Considere's construction gives
 * {@code dsigma/deps = sigma} and therefore {@code eps* = n}: a pressurised tube goes unstable
 * at 58 % of the strain the same material would reach in a tension test, because in biaxial
 * tension the geometry weakens twice as fast. {@link #instabilityStrain} reproduces
 * {@code n/sqrt(3)} to nine digits, which is what certifies the bisection.
 *
 * <p>All arguments SI: metres, pascals.
 */
public final class Burst {

    /** {@code sigma_theta / sigma_bar} for the 1 : 1/2 : 0 biaxial state: {@code 2/sqrt(3)}. */
    public static final double BIAXIAL = 2.0 / Math.sqrt(3.0);

    /** Logarithmic slope the flow curve must sustain to hold a pressurised tube stable. */
    public static final double BIAXIAL_SOFTENING = Math.sqrt(3.0);

    /** The same for uniaxial tension: Considere's {@code dsigma/deps = sigma}. */
    public static final double UNIAXIAL_SOFTENING = 1.0;

    private Burst() {
    }

    /**
     * Equivalent plastic strain at which {@code dsigma/deps = softening * sigma}, by bisection.
     *
     * <p>Unique and bracketed by construction, which is worth spelling out because it is what
     * makes bisection the right tool rather than a lazy one. The flow curve's slope
     * {@code B n eps^(n-1)} is decreasing and the curve itself is increasing, so the residual
     * is strictly monotone and has at most one root. At {@code eps = 1} it is
     * {@code B n - k (A + B)}, negative for every admissible law, since {@code n <= 1},
     * {@code k >= 1} and {@code A > 0}. <b>So no Johnson-Cook tube bursts beyond unit
     * equivalent strain</b>, and the upper bracket never has to be searched for.
     *
     * <p>Returns zero when the residual is negative everywhere, which is not a failure to
     * converge but a material statement: a law whose hardening cannot outrun the geometric
     * softening even at the first increment of plastic strain -- {@code B = 0}, or {@code n = 1}
     * with {@code B < k A} -- is unstable from first yield, and its burst pressure is its
     * elastic limit. {@link #burstPressure} then returns
     * {@code (2/sqrt(3)) A t/r}, which agrees with
     * {@link ElasticPlastic#limitPressure} in the thin-wall limit; they are independent
     * derivations and that agreement is a test.
     */
    public static double instabilityStrain(JohnsonCook law, double softening) {
        requireFlowCurve(law);
        double hi = 1.0;
        double lo = hi;
        while (residual(law, softening, lo) <= 0.0) {
            lo *= 0.5;
            if (lo < 1e-15) return 0.0;
        }
        for (int i = 0; i < 200; i++) {
            final double mid = 0.5 * (lo + hi);
            if (residual(law, softening, mid) > 0.0) lo = mid; else hi = mid;
        }
        return 0.5 * (lo + hi);
    }

    private static double residual(JohnsonCook law, double softening, double strain) {
        return law.strainSlope(strain) - softening * law.strainTerm(strain);
    }

    /** Equivalent plastic strain at burst. For a power law, exactly {@code n / sqrt(3)}. */
    public static double instabilityStrain(JohnsonCook law) {
        return instabilityStrain(law, BIAXIAL_SOFTENING);
    }

    /** Equivalent plastic strain at the uniaxial tensile instability. For a power law, {@code n}. */
    public static double tensileInstabilityStrain(JohnsonCook law) {
        return instabilityStrain(law, UNIAXIAL_SOFTENING);
    }

    /**
     * Engineering ultimate tensile strength: true stress at the tensile instability, converted
     * to the original cross-section. This is the number a data sheet quotes and the number
     * {@link #barlow} is meant to be fed.
     */
    public static double ultimateTensileStrength(JohnsonCook law) {
        final double strain = tensileInstabilityStrain(law);
        return law.strainTerm(strain) * Math.exp(-strain);
    }

    /** Logarithmic mid-wall hoop strain at burst: {@code (sqrt(3)/2)} of the equivalent strain. */
    public static double burstHoopStrain(JohnsonCook law) {
        return instabilityStrain(law) / BIAXIAL;
    }

    /**
     * Pressure a tube of initial thickness {@code t0} and mean radius {@code r0} can hold at a
     * given logarithmic mid-wall hoop strain. Rising while hardening wins, falling after.
     *
     * <p>The strain is taken as entirely plastic, which is what rigid-plastic means. See
     * {@link #elasticKnockdown} for the size of what that leaves out.
     */
    public static double pressure(JohnsonCook law, double thickness, double meanRadius,
                                  double hoopStrain) {
        requireFlowCurve(law);
        if (thickness <= 0.0 || meanRadius <= thickness) {
            throw new IllegalArgumentException("need a positive thickness inside the radius");
        }
        final double equivalent = BIAXIAL * hoopStrain;
        return BIAXIAL * law.strainTerm(equivalent) * (thickness / meanRadius)
                * Math.exp(-2.0 * hoopStrain);
    }

    /** Burst pressure: the maximum of {@link #pressure} over expansion. */
    public static double burstPressure(JohnsonCook law, double thickness, double meanRadius) {
        return pressure(law, thickness, meanRadius, burstHoopStrain(law));
    }

    /**
     * Barlow's formula, the engineering rule of thumb: {@code P = 2 sigma t / D}. Quote it on
     * the mean diameter and feed it the engineering UTS and it lands within a few per cent of
     * the instability solution, for reasons that are worth knowing rather than trusting --
     * see {@link #effectiveStress}.
     */
    public static double barlow(double stress, double thickness, double diameter) {
        return 2.0 * stress * thickness / diameter;
    }

    /**
     * The stress that would make {@link #barlow} on the mean diameter exactly right:
     * {@code (2/sqrt(3)) sigma_bar* exp(-2 eps_theta*)}.
     *
     * <p>Comparing this against {@link #ultimateTensileStrength} is the whole explanation of
     * why a formula with no geometry change and no instability in it works at all. Two large
     * errors very nearly cancel. Biaxiality <em>raises</em> the pressure a given flow stress
     * can hold, by {@code 2/sqrt(3)}, some 15 %; but it also makes the tube go unstable at
     * {@code 1/sqrt(3)} of the uniaxial strain, so it reaches that pressure with far less
     * hardening banked and, in exchange, far less thinning spent. What is left over is a few
     * per cent, and on the low side, so the rule of thumb is conservative. That is luck with
     * a reason, not a derivation, and it is why this class does not validate against Barlow.
     */
    public static double effectiveStress(JohnsonCook law) {
        final double strain = burstHoopStrain(law);
        return BIAXIAL * law.strainTerm(BIAXIAL * strain) * Math.exp(-2.0 * strain);
    }

    /**
     * The elastic part of the state at burst, which is the whole of the difference between
     * this oracle and the model a finite-strain solver actually integrates.
     *
     * <h2>Why elasticity moves a rigid-plastic answer at all</h2>
     *
     * Write the rigid-plastic curve as {@code P_rp(x) = C sigma(kappa x) exp(-2x)} in total
     * hoop strain {@code x}. Two things change when the material is allowed to be elastic, and
     * <b>they pull in opposite directions</b>:
     *
     * <ul>
     *   <li><b>The elastic hoop strain thins the wall without hardening the material.</b> Flow
     *       stress is set by the plastic strain {@code x - e}, while the geometry thins by the
     *       total {@code x}. Worth {@code -2e}.</li>
     *   <li><b>Elastic dilatation resists that thinning.</b> The {@code exp(-2x)} factor comes
     *       from incompressibility: {@code eps_r = -eps_theta}. A real material under a mean
     *       stress of {@code +p} swells by {@code D = (1-2nu) 3p/E}, and with {@code eps_z = 0}
     *       that goes entirely into {@code eps_r + eps_theta = D}, so the wall is thicker than
     *       incompressible kinematics predicts. Worth {@code +D}.</li>
     * </ul>
     *
     * <p>Together, {@code P(x) = exp(D - 2e) P_rp(x - e)}: the rigid-plastic curve
     * <b>translated by {@code e} and scaled by {@code exp(D - 2e)}</b>. Since neither
     * {@code e} nor {@code D} moves appreciably across the flat top of the curve, that is not
     * a bound but two statements: the peak sits at a hoop strain higher by {@code e}, and the
     * peak pressure is lower by {@code 1 - exp(D - 2e)}.
     *
     * <p>For copper the two terms are 0.385 % down and 0.105 % up, netting <b>0.279 %</b>, and
     * the measured deficit on the reference tube is 0.283 %. Keeping only the first term -- the
     * obvious one -- would over-predict the deficit by a third and leave a residual with no
     * home. That is the reason this is worth computing properly rather than waving at: two
     * terms of the same order, one of which is easy to forget.
     *
     * <h2>What it does not cover</h2>
     *
     * About 0.04 % of the measured deficit is left over after this, and it is worth saying so
     * rather than widening a tolerance until it disappears. It is independent of wall thickness
     * and independent of the loading rate -- both checked -- so it is neither the thin-wall
     * assumption nor the harness. The remaining suspect is that the solver's elastic response
     * is hypoelastic, a rate form integrated along the path, whose departure from an exact
     * hyperelastic law is of order the elastic strain times the total strain: 0.002 by 0.1,
     * which is the size observed. That is a hypothesis with the right magnitude and the right
     * invariances, not a measurement.
     *
     * @param hoopStress     {@code sigma_theta} at burst, Pa
     * @param axialStress    {@code sigma_z}, taken at its fully plastic value
     * @param radialStress   {@code sigma_r}, wall-averaged as {@code -P/2}
     * @param hoopStrain     the elastic hoop strain {@code e}
     * @param dilatation     the elastic volume strain {@code D}
     * @param knockdown      {@code 1 - exp(D - 2e)}, the fraction the peak pressure drops by
     */
    public record ElasticState(double hoopStress, double axialStress, double radialStress,
                               double hoopStrain, double dilatation, double knockdown) {
    }

    /**
     * The elastic state at the burst point of a given tube. See {@link ElasticState}.
     *
     * <p>Assembled from the yield condition {@code sigma_theta - sigma_r = (2/sqrt(3)) sigma_bar}
     * with {@code sigma_r} at its wall-average {@code -P/2} and {@code sigma_z} at the fully
     * plastic {@code (sigma_r + sigma_theta)/2}. The last of those is an approximation -- the
     * real plane-strain axial stress drifts up from {@code nu (sigma_r + sigma_theta)} as flow
     * accumulates -- but a harmless one, because {@code sigma_bar} is stationary with respect
     * to {@code sigma_z} exactly there, so the error it causes is second order.
     */
    public static ElasticState elasticState(JohnsonCook law, double thickness,
                                            double meanRadius, double youngsModulus,
                                            double poissonRatio) {
        final double pressure = burstPressure(law, thickness, meanRadius);
        final double radial = -0.5 * pressure;
        final double hoop = BIAXIAL * law.strainTerm(instabilityStrain(law)) + radial;
        final double axial = 0.5 * (radial + hoop);

        final double hoopStrain = (hoop - poissonRatio * (axial + radial)) / youngsModulus;
        final double dilatation =
                (1.0 - 2.0 * poissonRatio) * (radial + hoop + axial) / youngsModulus;

        return new ElasticState(hoop, axial, radial, hoopStrain, dilatation,
                1.0 - Math.exp(dilatation - 2.0 * hoopStrain));
    }

    /** Fraction by which elasticity lowers the burst pressure. See {@link #elasticState}. */
    public static double elasticKnockdown(JohnsonCook law, double thickness, double meanRadius,
                                          double youngsModulus, double poissonRatio) {
        return elasticState(law, thickness, meanRadius, youngsModulus, poissonRatio)
                .knockdown();
    }

    private static void requireFlowCurve(JohnsonCook law) {
        if (!law.isQuasiStatic()) {
            // With a rate term in the law the instability strain depends on how fast the
            // pump runs and there is no single number to predict; with self-heating it
            // depends on the loading history too. Neither belongs in a burst test that takes
            // minutes. See JohnsonCook.quasiStatic.
            throw new IllegalArgumentException(
                    "burst is a property of a flow curve: pass a rate-independent, isothermal "
                            + "law (JohnsonCook.quasiStatic), not one whose strength depends on "
                            + "the pumping rate");
        }
    }
}
