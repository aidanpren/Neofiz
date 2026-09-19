package org.neofiz.core;

/**
 * A two-parameter Weibull distribution, and the weakest-link algebra that comes with it.
 *
 * <pre>  P(a sample of unit reference volume survives x) = exp(-(x / lambda)^m)</pre>
 *
 * <p>The scale {@code lambda} is the characteristic value -- the level at which survival has
 * fallen to {@code 1/e} -- and the modulus {@code m} is the scatter parameter. A large modulus
 * is a tight distribution; a small one is a broad one. This carries the difference in
 * <em>character</em> between materials that a single strength number cannot: ductile structural
 * steels sit around m = 20-50, brittle ceramics and cast irons around m = 5-15.
 *
 * <h2>The one identity the rest of M2 rests on</h2>
 *
 * The minimum of {@code n} independent Weibull samples is <b>again Weibull</b>, with the same
 * modulus and a scale reduced by {@code n^(-1/m)}:
 *
 * <pre>  P(all n survive x) = exp(-(x/lambda)^m)^n = exp(-(x / (lambda n^(-1/m)))^m)</pre>
 *
 * <p>That is exact, not an approximation, and it is the reason the whole construction is cheap.
 * It means three statements that look like three separate features are one calculation:
 *
 * <ul>
 *   <li><b>Bigger parts are weaker.</b> A body of volume V contains V/V0 reference volumes, so
 *       its strength is the minimum of that many samples: {@code (V2/V1)^(1/m)}. Eight times
 *       the volume at m = 10 is about 19 % weaker.</li>
 *   <li><b>Volume normalisation.</b> An element holding volume V must be given the
 *       distribution of the minimum over V, not a fresh draw from the base law -- otherwise
 *       refining the mesh takes more samples of the same distribution, finds a worse extreme,
 *       and the part gets weaker without bound.</li>
 *   <li><b>An axisymmetric element is a ring.</b> A real flaw is a point; an element here is a
 *       full 360-degree torus, so its stored value has to stand for the weakest of the roughly
 *       {@code 2 pi r / l} independent flaws around that circumference. Using the true ring
 *       volume {@code 2 pi r A} in the normalisation <em>is</em> that statement -- there is no
 *       second correction to apply.</li>
 * </ul>
 *
 * <p>Because those are the same calculation, {@link #forVolume} is implemented as
 * {@link #minimumOf} rather than alongside it. A test asserting they agree would otherwise be
 * comparing two roundings of the same algebra, which is a test that can pass while the intent
 * has drifted.
 *
 * <h2>What the distribution is over</h2>
 *
 * Here it is over <b>local failure strain</b>, following the plan: one extra number per element,
 * sampled at mesh time. Nothing in this class knows that. It is equally the distribution of a
 * strength, a toughness or a time to failure, and the size effect reads the same way in each.
 *
 * @param scale   the characteristic value, at the reference volume the caller is using
 * @param modulus the Weibull modulus m; small is scattered, large is tight
 */
public record Weibull(double scale, double modulus) {

    /** Ductile structural steel: tight scatter, weak size effect. */
    public static final double MODULUS_DUCTILE_STEEL = 30.0;
    /** Cast iron and ceramics: broad scatter, strong size effect. */
    public static final double MODULUS_BRITTLE = 8.0;

    public Weibull {
        // Written as a negated positive test so that NaN is refused too. A NaN modulus
        // survives "modulus <= 0" and then silently turns every sampled failure strain into
        // NaN, which reaches the solver as an element that never yields and never fails.
        if (!(scale > 0.0)) {
            throw new IllegalArgumentException("scale must be positive, was " + scale);
        }
        if (!(modulus > 0.0)) {
            throw new IllegalArgumentException("modulus must be positive, was " + modulus);
        }
    }

    /**
     * {@code (x / lambda)^m}, the cumulative hazard: the expected number of failures a body
     * would suffer at level {@code x} if it could fail more than once.
     *
     * <p>Worth naming rather than leaving inside the exponential, because this is the
     * <b>additive</b> quantity and additivity is the structural fact underneath the whole
     * defect field. Cut a body into pieces however you like, give each piece the distribution
     * its own volume earns it, and the hazards sum back to the hazard of the whole:
     *
     * <pre>  sum_e (V_e/V0) (x/lambda)^m  =  (V/V0) (x/lambda)^m</pre>
     *
     * <p>So the survival probability of the body is <b>independent of the partition</b> --
     * which is to say, independent of the mesh. That is the same shape of guarantee as the
     * one the threaded assembly rests on, and for the same reason: nothing is allowed to
     * depend on how the work was divided up.
     */
    public double cumulativeHazard(double x) {
        if (x <= 0.0) return 0.0;
        return Math.pow(x / scale, modulus);
    }

    /** Probability that a sample of the reference volume exceeds {@code x}. */
    public double survival(double x) {
        if (x <= 0.0) return 1.0;
        return Math.exp(-cumulativeHazard(x));
    }

    /** Probability that a sample of the reference volume has failed by {@code x}. */
    public double cdf(double x) {
        if (x <= 0.0) return 0.0;
        return -Math.expm1(-cumulativeHazard(x));
    }

    /**
     * The value at cumulative probability {@code p}. This is the sampling transform: feed it a
     * uniform and get a Weibull.
     *
     * <p>Uses {@code log1p} rather than {@code log(1 - p)} because the interesting end is the
     * one where p is close to 1, and {@code 1 - p} throws away every digit of it. A uniform of
     * 0.999999 is a sample from the strong tail, and getting it wrong would make the strong
     * tail look like a hard ceiling.
     */
    public double quantile(double p) {
        if (!(p >= 0.0 && p <= 1.0)) {
            throw new IllegalArgumentException("probability must lie in [0, 1], was " + p);
        }
        return scale * Math.pow(-Math.log1p(-p), 1.0 / modulus);
    }

    /**
     * The value at survival probability {@code s}, which is {@code quantile(1 - s)} written so
     * that the weak tail keeps its digits. The weak tail is the one that fails first, so it is
     * the one a defect field spends its time in.
     */
    public double quantileFromSurvival(double s) {
        if (!(s >= 0.0 && s <= 1.0)) {
            throw new IllegalArgumentException("probability must lie in [0, 1], was " + s);
        }
        return scale * Math.pow(-Math.log(s), 1.0 / modulus);
    }

    /**
     * The distribution of the weakest of {@code n} independent samples. Same modulus, scale
     * multiplied by {@code n^(-1/m)}.
     *
     * <p>{@code n} is a double rather than an int on purpose: the count of independent flaws
     * around a ring, or in an element of arbitrary volume, is not a whole number and rounding
     * it would put a staircase into a size effect that is smooth.
     */
    public Weibull minimumOf(double n) {
        if (!(n > 0.0)) {
            throw new IllegalArgumentException("sample count must be positive, was " + n);
        }
        return new Weibull(scale * Math.pow(n, -1.0 / modulus), modulus);
    }

    /**
     * The distribution governing a body of volume {@code volume}, given that {@code scale} was
     * measured on specimens of {@code referenceVolume}.
     *
     * <p>This is {@link #minimumOf} with {@code n = volume / referenceVolume}, and is written
     * as a call to it rather than as the same formula a second time.
     */
    public Weibull forVolume(double volume, double referenceVolume) {
        if (!(volume > 0.0)) {
            throw new IllegalArgumentException("volume must be positive, was " + volume);
        }
        if (!(referenceVolume > 0.0)) {
            throw new IllegalArgumentException(
                    "reference volume must be positive, was " + referenceVolume);
        }
        return minimumOf(volume / referenceVolume);
    }

    /**
     * The factor by which a body of {@code volume} is weaker than one of
     * {@code referenceVolume}: {@code (V0/V)^(1/m)}, below one for the larger body.
     *
     * <p>The plan's headline number. It is quoted the other way round there -- strength scales
     * as {@code (V2/V1)^(1/m)} with V1 the larger -- which is the same statement.
     */
    public double sizeFactor(double volume, double referenceVolume) {
        return forVolume(volume, referenceVolume).scale() / scale;
    }

    /** Expected value, {@code lambda * Gamma(1 + 1/m)}. */
    public double mean() {
        return scale * gamma(1.0 + 1.0 / modulus);
    }

    /** Standard deviation, {@code lambda * sqrt(Gamma(1 + 2/m) - Gamma(1 + 1/m)^2)}. */
    public double standardDeviation() {
        double g1 = gamma(1.0 + 1.0 / modulus);
        double g2 = gamma(1.0 + 2.0 / modulus);
        return scale * Math.sqrt(Math.max(0.0, g2 - g1 * g1));
    }

    /**
     * Standard deviation over mean, which depends on the modulus alone. This is the number to
     * quote when explaining what a modulus means to someone who has never met one: m = 8 is
     * about 15 % scatter, m = 30 about 4 %.
     */
    public double coefficientOfVariation() {
        double g1 = gamma(1.0 + 1.0 / modulus);
        double g2 = gamma(1.0 + 2.0 / modulus);
        return Math.sqrt(Math.max(0.0, g2 - g1 * g1)) / g1;
    }

    // ---------------------------------------------------------------- gamma

    /**
     * Lanczos approximation, g = 7, nine coefficients. Good to about 15 significant figures
     * over the range that matters here, which is {@code 1 + 1/m} and {@code 1 + 2/m} for a
     * positive modulus and therefore always greater than one.
     *
     * <p>Present because the JDK has no gamma function and because the alternative -- quoting
     * a mean and a scatter from a Monte Carlo estimate -- would give the class a test that
     * agrees with itself to three digits and nothing more.
     */
    static double gamma(double x) {
        if (x < 0.5) {
            // Reflection. Not reachable from mean() or standardDeviation(), but a gamma that
            // is wrong outside the range it is called from is a trap for the next caller.
            return Math.PI / (Math.sin(Math.PI * x) * gamma(1.0 - x));
        }
        final double[] c = LANCZOS;
        double z = x - 1.0;
        double a = c[0];
        for (int i = 1; i < c.length; i++) a += c[i] / (z + i);
        double t = z + 7.5;
        return Math.sqrt(2.0 * Math.PI) * Math.pow(t, z + 0.5) * Math.exp(-t) * a;
    }

    private static final double[] LANCZOS = {
            0.99999999999980993,
            676.5203681218851,
            -1259.1392167224028,
            771.32342877765313,
            -176.61502916214059,
            12.507343278686905,
            -0.13857109526572012,
            9.9843695780195716e-6,
            1.5056327351493116e-7,
    };
}
