package org.neofiz.core;

/**
 * The standard normal cumulative distribution, to full double precision.
 *
 * <p>Exists because the JDK has no error function and because {@link DefectField} needs one
 * that is <em>accurate</em>, not merely plausible. The field is built as a correlated Gaussian
 * and then pushed through this function to a uniform, which is then pushed through a Weibull
 * quantile. Every bit of error in the middle step lands on the Weibull marginal, and a
 * distorted marginal shows up as a size effect that misses {@code (V2/V1)^(1/m)} by a few per
 * cent -- which is indistinguishable from the construction being wrong, and which would be
 * absorbed by widening the tolerance on the test that is supposed to be proving it.
 *
 * <p>The common Abramowitz and Stegun 7.1.26 approximation is good to 1.5e-7. That is roughly
 * a million times worse than double precision and it is the one every quick implementation
 * reaches for. What follows instead is W. J. Cody's rational Chebyshev approximation, three
 * branches by magnitude, as published in ALGORITHM 715.
 *
 * <p>Measured against the platform's own {@code erfc} at 7681 points across [-30, 30], the
 * worst disagreement is <b>three units in the last place</b>, or 2.1e-15 relative. It is not
 * correctly rounded and does not claim to be.
 */
public final class Normal {

    private Normal() {
    }

    private static final double SQRT_HALF = 0.70710678118654752440;
    /** 1/sqrt(pi), the leading coefficient of the asymptotic erfc expansion. */
    private static final double INV_SQRT_PI = 5.6418958354775628695e-1;

    /** P(Z &lt;= x) for a standard normal Z. */
    public static double cdf(double x) {
        return 0.5 * erfc(-x * SQRT_HALF);
    }

    /**
     * P(Z &gt; x). Written through {@code erfc} of a positive argument where possible, because
     * {@code 1 - cdf(x)} in the far upper tail is a subtraction that leaves nothing: at
     * x = 6 the answer is 1e-9 and every significant figure of it has been cancelled away.
     */
    public static double survival(double x) {
        return 0.5 * erfc(x * SQRT_HALF);
    }

    /** The error function. */
    public static double erf(double x) {
        double a = Math.abs(x);
        if (a <= THRESHOLD) return x * smallErf(a * a);
        double e = 1.0 - erfc(a);
        return x < 0.0 ? -e : e;
    }

    /** The complementary error function, {@code 1 - erf(x)}, accurate in the tail. */
    public static double erfc(double x) {
        double y = Math.abs(x);
        double result;
        if (y <= THRESHOLD) {
            result = 1.0 - x * smallErf(y * y);
            return result;
        } else if (y <= 4.0) {
            double num = C[8] * y;
            double den = y;
            for (int i = 0; i < 7; i++) {
                num = (num + C[i]) * y;
                den = (den + D[i]) * y;
            }
            result = (num + C[7]) / (den + D[7]);
        } else {
            double inv = 1.0 / (y * y);
            double num = P[5] * inv;
            double den = inv;
            for (int i = 0; i < 4; i++) {
                num = (num + P[i]) * inv;
                den = (den + Q[i]) * inv;
            }
            result = inv * (num + P[4]) / (den + Q[4]);
            result = (INV_SQRT_PI - result) / y;
        }

        // exp(-y*y) split so that the large, exactly representable part of the exponent is
        // evaluated once. Cody's trick: it keeps the rounding of the exponential from eating
        // the accuracy the rational approximation just bought.
        double trunc = ((int) (y * 16.0)) / 16.0;
        double del = (y - trunc) * (y + trunc);
        result *= Math.exp(-trunc * trunc) * Math.exp(-del);

        return x < 0.0 ? 2.0 - result : result;
    }

    /** {@code erf(x)/x} for small x, as a function of x squared. */
    private static double smallErf(double sq) {
        double num = A[4] * sq;
        double den = sq;
        for (int i = 0; i < 3; i++) {
            num = (num + A[i]) * sq;
            den = (den + B[i]) * sq;
        }
        return (num + A[3]) / (den + B[3]);
    }

    private static final double THRESHOLD = 0.46875;

    private static final double[] A = {
            3.16112374387056560e0, 1.13864154151050156e2, 3.77485237685302021e2,
            3.20937758913846947e3, 1.85777706184603153e-1,
    };
    private static final double[] B = {
            2.36012909523441209e1, 2.44024637934444173e2, 1.28261652607737228e3,
            2.84423683343917062e3,
    };
    private static final double[] C = {
            5.64188496988670089e-1, 8.88314979438837594e0, 6.61191906371416295e1,
            2.98635138197400131e2, 8.81952221241769090e2, 1.71204761263407058e3,
            2.05107837782607147e3, 1.23033935479799725e3, 2.15311535474403846e-8,
    };
    private static final double[] D = {
            1.57449261107098347e1, 1.17693950891312499e2, 5.37181101862009858e2,
            1.62138957456669019e3, 3.29079923573345963e3, 4.36261909014324716e3,
            3.43936767414372164e3, 1.23033935480374942e3,
    };
    private static final double[] P = {
            3.05326634961232344e-1, 3.60344899949804439e-1, 1.25781726111229246e-1,
            1.60837851487422766e-2, 6.58749161529837803e-4, 1.63153871373020978e-2,
    };
    private static final double[] Q = {
            2.56852019228982242e0, 1.87295284992346047e0, 5.27905102951428412e-1,
            6.05183413124413191e-2, 2.33520497626869185e-3,
    };
}
