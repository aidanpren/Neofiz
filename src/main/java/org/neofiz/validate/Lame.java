package org.neofiz.validate;

/**
 * Lame's closed-form solution for a thick-walled cylinder under internal pressure.
 *
 * <p>Tier 1 of the validation program: exact answers, run on every build. A deviation here
 * is a bug rather than a modelling choice, which is what makes this the right M0 gate.
 *
 * <p>Note for later: Barlow's thin-wall formula, sigma = P*d/(2*t), is the companion check
 * above a diameter-to-thickness ratio of about 20 and degrades below it, where this
 * solution takes over. Barlow's 1836 derivation is actually wrong as a thick-wall hoop
 * stress -- it violates the three-dimensional constitutive law -- and survives because it
 * happens to coincide with the correct thin-wall plane-stress yield pressure.
 *
 * <p>All arguments SI: metres, pascals.
 */
public final class Lame {

    private Lame() {
    }

    /** Radial stress at radius r. Equals -p at the bore and zero at the outside. */
    public static double radialStress(double r, double a, double b, double p) {
        double a2 = a * a, b2 = b * b;
        return p * a2 / (b2 - a2) * (1.0 - b2 / (r * r));
    }

    /** Hoop stress at radius r. Maximum at the bore, which is where a vessel fails. */
    public static double hoopStress(double r, double a, double b, double p) {
        double a2 = a * a, b2 = b * b;
        return p * a2 / (b2 - a2) * (1.0 + b2 / (r * r));
    }

    /** Axial stress under plane strain: nu times the sum of the in-plane stresses. */
    public static double axialStressPlaneStrain(double a, double b, double p, double nu) {
        double a2 = a * a, b2 = b * b;
        return 2.0 * nu * p * a2 / (b2 - a2);
    }

    /** Radial displacement with the ends fully restrained, eps_z = 0. */
    public static double radialDisplacementPlaneStrain(double r, double a, double b, double p,
                                                       double E, double nu) {
        double a2 = a * a, b2 = b * b;
        return (1.0 + nu) * p * a2 / (E * (b2 - a2)) * ((1.0 - 2.0 * nu) * r + b2 / r);
    }

    /**
     * Radial displacement with the ends open and unloaded, sigma_z = 0.
     *
     * <p>Worth keeping both. Note that {@link #radialStress} and {@link #hoopStress} contain
     * neither E nor nu: the Lame stresses are <em>identical</em> under plane stress and
     * plane strain. Only the displacements differ, by around 7% for steel. A validation
     * harness that checks stress alone therefore cannot detect a wrong end condition, and
     * will report a clean pass while solving a different problem. The displacement check is
     * the one that has teeth.
     */
    public static double radialDisplacementPlaneStress(double r, double a, double b, double p,
                                                       double E, double nu) {
        double a2 = a * a, b2 = b * b;
        return p * a2 / (E * (b2 - a2)) * ((1.0 - nu) * r + (1.0 + nu) * b2 / r);
    }

    /**
     * Barlow's thin-wall hoop stress, for the M1 comparison. Valid above roughly
     * d/t = 20; below that it under-predicts and Lame should be used instead.
     */
    public static double barlowHoopStress(double meanDiameter, double wallThickness, double p) {
        return p * meanDiameter / (2.0 * wallThickness);
    }
}
