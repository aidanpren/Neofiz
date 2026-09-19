package org.neofiz.core;

/**
 * Johnson-Cook flow stress: strain hardening, strain-rate hardening and thermal softening,
 * multiplied together.
 *
 * <pre>  sigma_y = (A + B eps^n) (1 + C ln(epsdot*)) (1 - T*^m)</pre>
 *
 * <p>with {@code epsdot* = epsdot / epsdot_0} the normalised equivalent plastic strain rate and
 * {@code T* = (T - T_room) / (T_melt - T_room)} the homologous temperature.
 *
 * <h2>Why all three, and why not one at a time</h2>
 *
 * The three terms are not independent refinements that can be added in any order and each
 * called an improvement. Rate hardening and thermal softening <b>act in opposite
 * directions</b>, and on a Taylor impact they are comparable in size: at 10^4 to 10^5 per
 * second the rate term is worth roughly 16 % on the flow stress, while the temperature rise
 * from the plastic work doing the deforming takes several per cent back. Shipping the rate
 * term alone would over-correct the mushroom and look like a regression; shipping the thermal
 * term alone would make it worse. Their partial cancellation is the physics, not a
 * convenience.
 *
 * <p>The multiplicative form is what makes this tractable: each term depends on one variable,
 * so the derivative needed by the return map is a sum of three products rather than a
 * coupled mess.
 *
 * <h2>Where the temperature comes from</h2>
 *
 * There is no heat equation here and there should not be. A Taylor impact lasts tens of
 * microseconds and the thermal diffusion length in steel over that time is under a
 * micrometre -- far below one element -- so the deformation is <b>adiabatic</b> and each
 * material point heats only itself:
 *
 * <pre>  T = T_room + beta W_p / (rho c_p)</pre>
 *
 * <p>where {@code W_p} is the plastic work per unit volume already accumulated at that point
 * and {@code beta} is the Taylor-Quinney fraction, the share of plastic work that becomes heat
 * rather than stored defect energy -- typically 0.9. This is why plastic work is stored per
 * Gauss point: the temperature field needs no storage of its own, because it is a function of
 * a quantity the return map was already integrating.
 *
 * <h2>Two guards that are not cosmetic</h2>
 *
 * <b>The rate term is clamped at the reference rate.</b> Below {@code epsdot_0} the logarithm
 * is negative, so the raw expression <em>softens</em> the material as it slows and passes
 * through zero strength at {@code epsdot* = exp(-1/C)}, which for these constants is about
 * 1e-31 per second. That is unreachable in a real run and immediately reachable in the first
 * Newton iteration, where the trial rate starts at zero. Holding the factor at 1 for
 * {@code epsdot < epsdot_0} is the standard treatment and it also makes the rate term's
 * derivative vanish there rather than blow up.
 *
 * <p><b>Thermal softening is clamped at zero.</b> Past the melt temperature {@code T*^m}
 * exceeds one and the raw expression returns a <em>negative</em> flow stress -- a material
 * that resists deformation backwards. A molten point should carry no deviatoric strength, so
 * the factor floors at zero.
 *
 * @param a                   quasi-static yield stress, Pa
 * @param b                   strain hardening coefficient, Pa
 * @param n                   strain hardening exponent
 * @param c                   strain-rate sensitivity
 * @param m                   thermal softening exponent
 * @param referenceStrainRate {@code epsdot_0}, 1/s
 * @param roomTemperature     K
 * @param meltTemperature     K
 * @param specificHeat        J/(kg K)
 * @param taylorQuinney       fraction of plastic work that becomes heat
 */
public record JohnsonCook(double a, double b, double n, double c, double m,
                          double referenceStrainRate, double roomTemperature,
                          double meltTemperature, double specificHeat,
                          double taylorQuinney) {

    public JohnsonCook {
        if (a <= 0) throw new IllegalArgumentException("A must be positive");
        if (b < 0) throw new IllegalArgumentException("B must be non-negative");
        if (n <= 0 || n > 1) throw new IllegalArgumentException("n must lie in (0, 1]");
        if (c < 0) throw new IllegalArgumentException("C must be non-negative");
        if (m <= 0) throw new IllegalArgumentException("m must be positive");
        if (referenceStrainRate <= 0) {
            throw new IllegalArgumentException("reference strain rate must be positive");
        }
        if (meltTemperature <= roomTemperature) {
            throw new IllegalArgumentException("melt temperature must exceed room temperature");
        }
        if (specificHeat <= 0) throw new IllegalArgumentException("specific heat must be positive");
        if (taylorQuinney < 0 || taylorQuinney > 1) {
            throw new IllegalArgumentException("Taylor-Quinney fraction must lie in [0, 1]");
        }
    }

    /**
     * AISI 4340 steel, the original Johnson-Cook 1983 constants and the build plan's
     * Taylor-impact validation alloy.
     */
    public static final JohnsonCook STEEL_4340 = new JohnsonCook(
            792.0e6, 510.0e6, 0.26, 0.014, 1.03,
            1.0, 293.0, 1793.0, 477.0, 0.9);

    /** OFHC copper, for the cylinder expansion test at M5. */
    public static final JohnsonCook COPPER_OFHC = new JohnsonCook(
            90.0e6, 292.0e6, 0.31, 0.025, 1.09,
            1.0, 293.0, 1356.0, 383.0, 0.9);

    /**
     * The same law stripped to its strain term alone: {@code C = 0} and {@code beta = 0}, so
     * the flow stress is exactly {@code A + B eps^n} at every rate and every accumulated
     * plastic work.
     *
     * <p>This is not a simplification for convenience, it is the right material model for a
     * slow test. A hydrostatic burst test runs over minutes: the strain rate is 1e-4 per
     * second or less, and the wall has ample time to conduct away the few kelvin the
     * deformation generates. Neither of the two terms this drops is <em>small</em> under
     * those conditions -- they are inapplicable, being an extrapolation of a logarithm far
     * below its reference rate and an adiabatic assumption over a minute.
     *
     * <p>There is a second reason, and it is the one that matters for validation. The burst
     * instability condition is a statement about a flow <em>curve</em>: stress as a function
     * of strain alone. Leave the rate term in and the instability strain acquires a dependence
     * on how fast the pump runs, so there is no longer a single number for the closed form to
     * predict. Turning the terms off is what makes {@link org.neofiz.validate.Burst} an exact
     * oracle rather than an approximate one, and it is enforced rather than hoped for:
     * {@code C = 0} makes {@link #rateTerm} identically one at any rate, and {@code beta = 0}
     * holds every point at room temperature, so no run can quietly violate it.
     */
    public JohnsonCook quasiStatic() {
        return new JohnsonCook(a, b, n, 0.0, m, referenceStrainRate,
                roomTemperature, meltTemperature, specificHeat, 0.0);
    }

    /** Whether this law reduces to a flow curve: no rate term, no self-heating. */
    public boolean isQuasiStatic() {
        return c == 0.0 && taylorQuinney == 0.0;
    }

    // ------------------------------------------------------------------ the three terms

    /** Strain hardening, {@code A + B eps^n}. Pa. */
    public double strainTerm(double plasticStrain) {
        return plasticStrain <= 0.0 ? a : a + b * Math.pow(plasticStrain, n);
    }

    /**
     * Derivative of the strain term, {@code B n eps^(n-1)}. Pa.
     *
     * <p>Singular at zero strain, because {@code n < 1} makes the very first increment of
     * plastic strain infinitely stiff. That is a property of the law rather than a defect,
     * and it is why the return map's Newton solve is bracketed: an unguarded Newton step
     * computed from an infinite slope goes nowhere.
     */
    public double strainSlope(double plasticStrain) {
        return plasticStrain <= 0.0 ? Double.POSITIVE_INFINITY
                : b * n * Math.pow(plasticStrain, n - 1.0);
    }

    /**
     * The same derivative, recovered from a {@link #strainTerm} already computed at the same
     * strain, with no second power. Since {@code B eps^n} is the strain term less {@code A},
     *
     * <pre>  B n eps^(n-1) = n (B eps^n) / eps = n (S - A) / eps</pre>
     *
     * <p>which is algebra rather than an approximation. It exists because the return map's
     * Newton solve needs both at every iteration, and a {@code pow} is tens of nanoseconds
     * against a handful for a divide -- the Johnson-Cook path costs several times the closed
     * form and two thirds of that is transcendentals.
     *
     * <p>The subtraction cancels when {@code B eps^n} is small against {@code A}, which is the
     * regime just past first yield. It cannot do any harm there: the true slope is enormous, a
     * cancelled one is merely differently enormous or zero, and either produces a Newton step
     * the bracket rejects in favour of a bisection. The guarantee that the solve converges
     * does not depend on the slope being accurate, only on the bracket being real.
     *
     * @param strainTerm    the value {@link #strainTerm} returned at this strain
     * @param plasticStrain strictly positive; the caller's bracket guarantees it
     */
    public double strainSlopeFrom(double strainTerm, double plasticStrain) {
        return n * (strainTerm - a) / plasticStrain;
    }

    /** Rate hardening, {@code 1 + C ln(epsdot/epsdot_0)}, held at 1 below the reference rate. */
    public double rateTerm(double plasticStrainRate) {
        if (plasticStrainRate <= referenceStrainRate) return 1.0;
        return 1.0 + c * Math.log(plasticStrainRate / referenceStrainRate);
    }

    /** Whether the rate term is on its clamp, and so contributing no derivative. */
    public boolean rateIsClamped(double plasticStrainRate) {
        return plasticStrainRate <= referenceStrainRate;
    }

    /** Thermal softening, {@code 1 - T*^m}, floored at zero for a molten point. */
    public double thermalTerm(double temperature) {
        final double homologous =
                (temperature - roomTemperature) / (meltTemperature - roomTemperature);
        if (homologous <= 0.0) return 1.0;
        if (homologous >= 1.0) return 0.0;
        return 1.0 - Math.pow(homologous, m);
    }

    // ------------------------------------------------------------------ assembled

    /** Flow stress, Pa. */
    public double flowStress(double plasticStrain, double plasticStrainRate, double temperature) {
        return strainTerm(plasticStrain) * rateTerm(plasticStrainRate) * thermalTerm(temperature);
    }

    /** Flow stress at the reference rate, which is the surface the yield check uses. */
    public double quasiStaticFlowStress(double plasticStrain, double temperature) {
        return strainTerm(plasticStrain) * thermalTerm(temperature);
    }

    /**
     * Adiabatic temperature from plastic work density, K. See the class comment: no heat
     * equation, because nothing conducts anywhere on this timescale.
     *
     * @param plasticWorkDensity J/m^3
     * @param density            kg/m^3
     */
    public double temperature(double plasticWorkDensity, double density) {
        return roomTemperature
                + taylorQuinney * plasticWorkDensity / (density * specificHeat);
    }

    /** The coefficient turning plastic work density into a temperature rise, K m^3/J. */
    public double thermalCoefficient(double density) {
        return taylorQuinney / (density * specificHeat);
    }
}
