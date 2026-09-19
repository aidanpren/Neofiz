package org.neofiz.validate;

/**
 * Closed-form elastic-plastic thick-walled cylinder under internal pressure. Plane strain,
 * von Mises, elastic-perfectly-plastic, small strain.
 *
 * <p>Tier 1 validation, and a much sharper instrument than {@link Lame}. The elastic case
 * tests the stiffness assembly and little else -- any code that integrates B transpose sigma
 * correctly will pass it. This one cannot be passed by accident: it requires the yield
 * criterion, the flow rule, the return map and equilibrium to all be right simultaneously,
 * and it is sensitive to the difference between von Mises and Tresca, which differ here by
 * the factor of 2/sqrt(3) that shows up in every formula below.
 *
 * <h2>The solution</h2>
 *
 * Plasticity spreads outward from the bore as a front. Inside it, equilibrium plus the yield
 * condition sigma_theta - sigma_r = 2 sy / sqrt(3) integrate directly; outside it, the
 * material is still elastic and carries a Lame field that happens to be exactly at yield on
 * its inner face. Matching the two at the front radius c gives
 *
 * <pre>  p(c) = (sy / sqrt(3)) * [ 2 ln(c/a) + (1 - c^2/b^2) ]</pre>
 *
 * <p>Note what is <em>not</em> in that expression: neither E nor nu. The stress field of a
 * contained plastic zone is statically determinate, so it depends only on the yield stress
 * and the geometry. That is the same trap {@link Lame} carries -- a stress-only comparison
 * cannot see an error in the elastic constants -- so the front radius is checked too, since
 * where the front sits does depend on the elastic response that put it there.
 *
 * <p>All arguments SI: metres, pascals.
 */
public final class ElasticPlastic {

    /** Plane-strain von Mises effective yield in the hoop-radial difference: 2 sy / sqrt(3). */
    public static final double MISES_PLANE_STRAIN = 2.0 / Math.sqrt(3.0);

    private ElasticPlastic() {
    }

    /** Internal pressure that puts the elastic-plastic front at radius c. */
    public static double pressureForFront(double c, double a, double b, double sy) {
        if (c < a || c > b) throw new IllegalArgumentException("front must lie within the wall");
        return sy / Math.sqrt(3.0) * (2.0 * Math.log(c / a) + (1.0 - c * c / (b * b)));
    }

    /** Pressure at which the bore first yields. Below this the response is entirely elastic. */
    public static double elasticLimitPressure(double a, double b, double sy) {
        return pressureForFront(a, a, b, sy);
    }

    /**
     * Pressure at which the plastic zone reaches the outer surface. Beyond this the cylinder
     * has no remaining elastic constraint and bursts; nothing here models what follows.
     */
    public static double limitPressure(double a, double b, double sy) {
        return pressureForFront(b, a, b, sy);
    }

    /**
     * Front radius for a given pressure, by bisection on the monotone {@link
     * #pressureForFront}. Returns a when the pressure is below the elastic limit.
     */
    public static double frontRadius(double p, double a, double b, double sy) {
        if (p <= elasticLimitPressure(a, b, sy)) return a;
        if (p >= limitPressure(a, b, sy)) return b;
        double lo = a, hi = b;
        for (int i = 0; i < 200; i++) {
            double mid = 0.5 * (lo + hi);
            if (pressureForFront(mid, a, b, sy) < p) lo = mid; else hi = mid;
        }
        return 0.5 * (lo + hi);
    }

    /** Radial stress inside the plastic zone, r &lt;= c. */
    public static double radialStressPlastic(double r, double a, double p, double sy) {
        return -p + MISES_PLANE_STRAIN * sy * Math.log(r / a);
    }

    /** Hoop stress inside the plastic zone. Exactly one yield stress above the radial one. */
    public static double hoopStressPlastic(double r, double a, double p, double sy) {
        return radialStressPlastic(r, a, p, sy) + MISES_PLANE_STRAIN * sy;
    }

    /** Radial stress in the elastic annulus outside the front, r &gt;= c. */
    public static double radialStressElastic(double r, double c, double b, double sy) {
        return sy * c * c / (Math.sqrt(3.0) * b * b) * (1.0 - b * b / (r * r));
    }

    /** Hoop stress in the elastic annulus outside the front. */
    public static double hoopStressElastic(double r, double c, double b, double sy) {
        return sy * c * c / (Math.sqrt(3.0) * b * b) * (1.0 + b * b / (r * r));
    }

    /**
     * Radial stress anywhere in the wall, choosing the branch by comparison with the front.
     *
     * <p>Below the elastic limit this is plain {@link Lame}, not the elastic-annulus formula
     * with the front collapsed onto the bore. Those are not the same thing: the annulus
     * formula assumes its inner face is exactly at yield, so with c = a it returns the field
     * at the elastic limit pressure regardless of what pressure was asked for. Blending the
     * two at the boundary is what makes this continuous across first yield.
     */
    public static double radialStress(double r, double a, double b, double p, double sy) {
        if (p <= elasticLimitPressure(a, b, sy)) return Lame.radialStress(r, a, b, p);
        double c = frontRadius(p, a, b, sy);
        return r <= c ? radialStressPlastic(r, a, p, sy) : radialStressElastic(r, c, b, sy);
    }

    /** Hoop stress anywhere in the wall. */
    public static double hoopStress(double r, double a, double b, double p, double sy) {
        if (p <= elasticLimitPressure(a, b, sy)) return Lame.hoopStress(r, a, b, p);
        double c = frontRadius(p, a, b, sy);
        return r <= c ? hoopStressPlastic(r, a, p, sy) : hoopStressElastic(r, c, b, sy);
    }
}
