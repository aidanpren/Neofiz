package org.neofiz.core;

/**
 * Material constants. SI throughout: Pa, dimensionless, kg/m^3, K.
 *
 * <p>Carries the elastic constants, a linear yield surface with combined isotropic and
 * kinematic hardening, and optionally a {@link JohnsonCook} flow stress that supersedes the
 * linear one. Damage (Johnson-Cook D1-D5) is not here yet.
 *
 * <p><b>The two hardening descriptions are exclusive, not additive.</b> When a Johnson-Cook
 * law is attached it <em>is</em> the flow stress, and {@link #isotropicHardening} is ignored
 * -- so {@link #withJohnsonCook} zeroes it rather than leaving a stale value that reads as if
 * it were still contributing. {@link #yieldStress} is kept in step with the law's quasi-static
 * yield {@code A}, because it is what the elastic trial checks against and a disagreement
 * between the two would put the yield check and the return map on different surfaces.
 */
public record Material(String name, double youngsModulus, double poissonRatio, double density,
                       double yieldStress, double isotropicHardening, double kinematicHardening,
                       JohnsonCook johnsonCook) {

    public Material {
        if (youngsModulus <= 0) throw new IllegalArgumentException("E must be positive");
        if (poissonRatio <= -1.0 || poissonRatio >= 0.5) {
            throw new IllegalArgumentException("nu must lie in (-1, 0.5); got " + poissonRatio);
        }
        if (density <= 0) throw new IllegalArgumentException("density must be positive");
        if (yieldStress <= 0) throw new IllegalArgumentException("yield stress must be positive");
        if (isotropicHardening < 0 || kinematicHardening < 0) {
            throw new IllegalArgumentException("hardening moduli must be non-negative");
        }
        if (johnsonCook != null && yieldStress != johnsonCook.a()) {
            throw new IllegalArgumentException(
                    "yield stress " + yieldStress + " disagrees with the Johnson-Cook A of "
                            + johnsonCook.a() + "; the elastic check and the return map would "
                            + "be working on different surfaces");
        }
    }

    /** Rate-independent linear hardening. */
    public Material(String name, double youngsModulus, double poissonRatio, double density,
                    double yieldStress, double isotropicHardening, double kinematicHardening) {
        this(name, youngsModulus, poissonRatio, density,
                yieldStress, isotropicHardening, kinematicHardening, null);
    }

    /**
     * A purely elastic material, expressed as a plastic one whose yield stress is infinite.
     *
     * <p>This is not a trick to save a field. It means the return-mapping kernel has exactly
     * one path: an infinite yield surface is never crossed, so the elastic case falls out of
     * the same arithmetic rather than being special-cased around it. One path is one thing
     * to validate, and the elastic gates keep testing the plastic code.
     */
    public Material(String name, double youngsModulus, double poissonRatio, double density) {
        this(name, youngsModulus, poissonRatio, density, Double.POSITIVE_INFINITY, 0.0, 0.0);
    }

    public boolean isElastic() {
        return Double.isInfinite(yieldStress);
    }

    /**
     * Whether the flow stress depends on how fast the material is deforming. This is what
     * decides whether the return map has a closed form or needs a scalar solve.
     */
    public boolean isRateDependent() {
        return johnsonCook != null;
    }

    /** The same material with a yield surface and linear isotropic hardening attached. */
    public Material yielding(double yieldStress, double isotropicHardening) {
        return new Material(name, youngsModulus, poissonRatio, density,
                yieldStress, isotropicHardening, 0.0);
    }

    /** The same material with linear hardening split between isotropic and kinematic. */
    public Material yielding(double yieldStress, double isotropic, double kinematic) {
        return new Material(name, youngsModulus, poissonRatio, density,
                yieldStress, isotropic, kinematic);
    }

    /**
     * The same material with a Johnson-Cook flow stress, which replaces the linear hardening
     * rather than adding to it. Kinematic hardening is carried through: it translates the
     * surface, which is orthogonal to how Johnson-Cook sizes it.
     */
    public Material withJohnsonCook(JohnsonCook law) {
        return new Material(name, youngsModulus, poissonRatio, density,
                law.a(), 0.0, kinematicHardening, law);
    }

    /** First Lame parameter. */
    public double lambda() {
        return youngsModulus * poissonRatio / ((1.0 + poissonRatio) * (1.0 - 2.0 * poissonRatio));
    }

    /** Second Lame parameter, the shear modulus G. */
    public double mu() {
        return youngsModulus / (2.0 * (1.0 + poissonRatio));
    }

    /**
     * Dilatational (p-wave) speed, sqrt((lambda + 2*mu) / rho). This is the speed that
     * sets the CFL limit: the timestep must not let a pressure wave cross an element.
     */
    public double dilatationalWaveSpeed() {
        return Math.sqrt((lambda() + 2.0 * mu()) / density);
    }

    /** Shear wave speed, sqrt(mu / rho). */
    public double shearWaveSpeed() {
        return Math.sqrt(mu() / density);
    }

    /**
     * Bar (extensional) wave speed, sqrt(E / rho). This is the speed that governs
     * <em>structural</em> response -- a ring breathing, a rod ringing -- as distinct from
     * the dilatational speed, which governs the stability limit.
     *
     * <p>The two are not interchangeable and they diverge exactly where it hurts. As nu
     * approaches 0.5, lambda runs away and the dilatational speed with it, while this one
     * barely moves: for steel at nu = 0.29 they are 5850 and 5111 m/s, but at nu = 0.499
     * they are 66000 and 5111. Anything structural that is timed off the dilatational speed
     * is therefore correct for ordinary metals and silently wrong for anything
     * near-incompressible, which is the worst way for a constant to be wrong.
     */
    public double barWaveSpeed() {
        return Math.sqrt(youngsModulus / density);
    }

    /** AISI 4340, the Taylor-impact validation alloy at M1. */
    public static final Material STEEL_4340 = new Material("AISI 4340", 205.0e9, 0.29, 7850.0);

    /** Annealed OFHC copper, the cylinder-expansion validation material at M5. */
    public static final Material COPPER_OFHC = new Material("OFHC Cu", 117.0e9, 0.34, 8960.0);
}
