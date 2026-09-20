package org.neofiz.sweep;

import java.util.Locale;

/**
 * A power law fitted to a sweep: {@code y = a x^b}, by least squares in log-log.
 *
 * <p>This is the step that turns a column of numbers into something a person has learned.
 * The build plan is explicit that a log-log fit is sufficient and robust for this and that
 * symbolic regression is a later luxury -- the value is not in the sophistication of the fit
 * but in handing back an <b>exponent</b>, because the exponent is the physics. A sweep of
 * wall thickness against burst pressure that comes back with {@code b = 1.00} is Barlow's
 * formula, measured rather than read, by someone who was never told to expect it.
 *
 * <h2>What the residual is for</h2>
 *
 * A power law fitted to data that is not a power law still returns an exponent, and that
 * exponent is worse than no answer because it looks like one. So {@link #maxResidualPercent}
 * is reported beside {@code b} and is meant to be looked at first: it is the largest distance
 * any point sits from the fitted line, as a percentage of the value there. A clean scaling
 * law over a decade lands in the tenths of a per cent. Anything that does not is telling you
 * that the relationship has structure the exponent is flattening out -- which is itself worth
 * seeing, and is exactly what a thick-walled tube does to Barlow.
 *
 * <p>{@code rSquared} is reported too because people expect it, and it is the weaker of the
 * two statistics here: over a wide range in {@code x} a power law fit has an {@code r^2} of
 * 0.999 or better almost regardless of how badly it is wrong in the middle.
 */
public record PowerLaw(double exponent, double coefficient, double maxResidualPercent,
                       double rSquared, int points) {

    /**
     * Fits {@code y = a x^b} to the given points.
     *
     * <p>Both arrays must be strictly positive: a logarithm is taken of each, and a sweep
     * that produced a zero or a negative value did not produce a point on a power law. That
     * is refused rather than filtered, because silently dropping a run is how a sweep comes
     * to describe a different experiment than the one that was asked for.
     */
    public static PowerLaw fit(double[] x, double[] y) {
        if (x.length != y.length) {
            throw new IllegalArgumentException("x and y differ in length: "
                    + x.length + " and " + y.length);
        }
        if (x.length < 2) {
            throw new IllegalArgumentException("a power law needs at least two points");
        }

        final int n = x.length;
        final double[] logX = new double[n];
        final double[] logY = new double[n];
        for (int i = 0; i < n; i++) {
            if (!(x[i] > 0.0) || !(y[i] > 0.0)) {
                throw new IllegalArgumentException(String.format(Locale.ROOT,
                        "point %d is (%g, %g); a power law is fitted in the logarithm and "
                                + "needs both coordinates positive", i, x[i], y[i]));
            }
            logX[i] = Math.log(x[i]);
            logY[i] = Math.log(y[i]);
        }

        double meanX = 0.0;
        double meanY = 0.0;
        for (int i = 0; i < n; i++) {
            meanX += logX[i];
            meanY += logY[i];
        }
        meanX /= n;
        meanY /= n;

        double covariance = 0.0;
        double variance = 0.0;
        for (int i = 0; i < n; i++) {
            covariance += (logX[i] - meanX) * (logY[i] - meanY);
            variance += (logX[i] - meanX) * (logX[i] - meanX);
        }
        if (variance == 0.0) {
            throw new IllegalArgumentException(
                    "every point has the same x; there is no sweep here to fit");
        }

        final double exponent = covariance / variance;
        final double coefficient = Math.exp(meanY - exponent * meanX);

        double residualSum = 0.0;
        double totalSum = 0.0;
        double worst = 0.0;
        for (int i = 0; i < n; i++) {
            final double predicted = coefficient * Math.pow(x[i], exponent);
            residualSum += (logY[i] - Math.log(predicted)) * (logY[i] - Math.log(predicted));
            totalSum += (logY[i] - meanY) * (logY[i] - meanY);
            worst = Math.max(worst, Math.abs(y[i] - predicted) / y[i]);
        }

        return new PowerLaw(exponent, coefficient, 100.0 * worst,
                totalSum == 0.0 ? 1.0 : 1.0 - residualSum / totalSum, n);
    }

    /** The fitted value at {@code x}. */
    public double at(double x) {
        return coefficient * Math.pow(x, exponent);
    }

    /**
     * The law as a person would write it, with the exponent to two decimals because that is
     * the precision at which it is a statement about physics rather than about arithmetic.
     */
    public String describe(String yName, String xName) {
        return String.format(Locale.ROOT, "%s proportional to %s^%.2f", yName, xName, exponent);
    }

    @Override
    public String toString() {
        return String.format(Locale.ROOT,
                "y = %.4g x^%.4f  (%d points, worst residual %.3f %%, r2 = %.6f)",
                coefficient, exponent, points, maxResidualPercent, rSquared);
    }
}
