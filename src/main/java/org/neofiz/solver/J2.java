package org.neofiz.solver;

import org.neofiz.core.JohnsonCook;
import org.neofiz.core.Material;

/**
 * J2 (von Mises) plasticity by radial return, with combined isotropic and kinematic
 * hardening. The flow stress is either linear in plastic strain and rate-independent, or
 * Johnson-Cook.
 *
 * <p>This is the first piece of the code that has <b>memory</b>. Everything in M0 was a
 * function of the current displacement: wipe the state, re-impose the displacement, get the
 * same stress. A plastic material does not work that way -- where it is now depends on how
 * it got here -- so stress, back stress and accumulated plastic strain become stored state
 * that the solver integrates along the strain path. That is the real architectural change
 * plasticity brings, and it is why the solver has to become incremental rather than total.
 *
 * <h2>The algorithm</h2>
 *
 * Operator split. Take a fully elastic trial step, then, if that trial lands outside the
 * yield surface, project it back. The projection is <em>radial</em>: it moves along the
 * normal to the surface, which for von Mises means straight toward the surface in deviatoric
 * space, leaving the direction of the deviator unchanged.
 *
 * <pre>
 *   s      = dev(sigma_trial)
 *   xi     = s - alpha                      relative stress
 *   f      = ||xi|| - sqrt(2/3) * K(epsP)   yield function
 *   dgamma = f / (2mu + (2/3)(H_iso + H_kin))
 *   sigma  = sigma_trial - 2 mu dgamma n    where n = xi / ||xi||
 * </pre>
 *
 * <p>Two properties fall out of this exactly rather than approximately, and both are
 * asserted to machine precision in the tests rather than to a tolerance:
 *
 * <ul>
 *   <li><b>For linear hardening the return is exact, not iterative.</b> The consistency
 *       condition is linear in dgamma, so the closed form above lands precisely on the
 *       updated yield surface. No Newton loop, no convergence tolerance, no failure mode
 *       where a badly conditioned element quietly stops converging.
 *
 *       <p>Johnson-Cook breaks that, because {@code A + B eps^n} is nonlinear in plastic
 *       strain and the rate term depends on dgamma itself -- the multiplier appears on both
 *       sides. Only the solve for dgamma changes; everything after it is the same arithmetic.
 *       The guarantee is restored rather than abandoned: see {@link #solve}, which brackets
 *       the root and falls back to bisection, so the iteration cannot fail to converge on a
 *       badly conditioned point. It can only get there in more steps.</li>
 *   <li><b>Plastic flow is exactly volume-preserving.</b> n is deviatoric by construction,
 *       so the return alters no part of the pressure. Metals do not change volume when they
 *       yield, and a scheme that leaks even slightly is a scheme that will manufacture or
 *       destroy energy over 200,000 steps.</li>
 * </ul>
 *
 * <h2>Why both hardening kinds</h2>
 *
 * Isotropic hardening grows the yield surface; kinematic hardening translates it. They are
 * indistinguishable under monotonic loading and completely different under reversal, which
 * is the Bauschinger effect: a metal pulled into tension yields <em>early</em> when pushed
 * back into compression. Anything that cycles -- a barrel breathing shot after shot, a case
 * wall flexing -- is a reversal problem, and a purely isotropic model will over-predict the
 * reverse yield point every time. The two moduli are carried separately so the split is a
 * material property rather than a hidden assumption.
 *
 * <h2>Component ordering</h2>
 *
 * Four components, {rr, zz, theta-theta, rz}, which covers axisymmetric and plane strain.
 * The shear is the tensor component, so the norm of a deviator carries the factor of two:
 * {@code ||s||^2 = s_rr^2 + s_zz^2 + s_tt^2 + 2 s_rz^2}. Strain increments arrive with the
 * engineering shear, twice the tensor shear, which is why {@code mu * dGamRZ} and not
 * {@code 2 mu * dGamRZ} appears below.
 *
 * <p>Plane <em>stress</em> is deliberately absent. Enforcing sigma_zz = 0 through a return
 * map is a genuinely different algorithm, not a special case of this one, and nothing in the
 * validation program needs it: the open-ended cylinder is solved as a true axisymmetric
 * problem with a free axial direction, not as plane stress.
 */
public final class J2 {

    /** Components per stress tensor. */
    public static final int COMPONENTS = 4;

    private static final double SQRT_2_3 = Math.sqrt(2.0 / 3.0);
    private static final double TWO_THIRDS = 2.0 / 3.0;

    /** Residual tolerance for the Johnson-Cook solve, relative to the quasi-static yield. */
    private static final double SOLVE_TOLERANCE = 1.0e-10;
    /** Insurance only. Newton reaches the tolerance in a handful of steps; see {@link #solve}. */
    private static final int SOLVE_ITERATIONS = 25;

    private J2() {
    }

    /**
     * The plasticity parameters, hoisted out of the element loop.
     *
     * <p>Bundled rather than passed as eight more scalars, and built once per solver rather
     * than per element: these are <em>uniforms</em>, constant across the whole mesh, which is
     * exactly what becomes a constant buffer when this kernel is ported. The per-point state
     * stays in the flat arrays where it belongs.
     *
     * @param thermalCoefficient {@code beta / (rho c_p)}, precomputed because it is the only
     *                           thing the temperature needs from the density, and zero when
     *                           there is no Johnson-Cook law to soften
     */
    public record Flow(double sy0, double hIso, double hKin,
                       JohnsonCook jc, double thermalCoefficient) {

        public static Flow of(Material material) {
            JohnsonCook law = material.johnsonCook();
            return new Flow(material.yieldStress(), material.isotropicHardening(),
                    material.kinematicHardening(), law,
                    law == null ? 0.0 : law.thermalCoefficient(material.density()));
        }

        /** Linear hardening, no rate or temperature dependence. */
        public static Flow linear(double sy0, double hIso, double hKin) {
            return new Flow(sy0, hIso, hKin, null, 0.0);
        }
    }

    /**
     * Advances one element's stress and history through a strain increment, in place.
     *
     * @param sig     stress, stride 4, {rr, zz, tt, rz}
     * @param back    back stress (deviatoric), stride 4
     * @param epsP    accumulated equivalent plastic strain, stride 1
     * @param work     accumulated plastic work per unit volume, J/m^3, stride 1
     * @param e        element index
     * @param dGamRZ   <em>engineering</em> shear strain increment, twice the tensor component
     * @param flow     the flow stress law and its constants
     * @param strainDt the interval the strain increment was taken over, which is what turns
     *                 a plastic strain increment into a plastic strain <em>rate</em>. Zero
     *                 assembles without advancing the material.
     * @return the plastic multiplier, zero if the step was elastic
     */
    public static double update(double[] sig, double[] back, double[] epsP, double[] work, int e,
                                double dRR, double dZZ, double dTT, double dGamRZ,
                                double lambda, double mu, Flow flow, double strainDt) {
        return update(sig, back, epsP, work, e, dRR, dZZ, dTT, dGamRZ,
                lambda, mu, flow, strainDt, null);
    }

    /**
     * The same update with a damage state attached.
     *
     * <p>{@link Damage} enters as a <b>negative isotropic hardening modulus</b> and nothing
     * else. Once a point has softened, its frozen flow stress and its softening slope take the
     * places that {@code sy0 + hIso * epsP} and {@code hIso} occupy for an undamaged one, and
     * every line after that is unchanged -- including the consistency condition, which stays
     * linear in dgamma and therefore stays exact. A softening Johnson-Cook point does not even
     * reach {@link #solve}: freezing the flow stress at onset removes the strain, rate and
     * temperature dependence that made the solve necessary, so the branch that needed an
     * iteration is the one branch that no longer has one.
     *
     * <p>The one thing that has to be got right is the sign of the denominator. It is
     * {@code 2 mu + (2/3)(H_soft + H_kin)}, and if softening were steep enough to make it
     * negative the return would push the stress the wrong way. That is why {@link Damage}
     * clamps the slope at {@code -E}: the worst case is then {@code 2 mu - (2/3) E}, which is
     * positive for every Poisson ratio below 0.5.
     *
     * @param damage the damage state, or null for a material that cannot fail
     */
    public static double update(double[] sig, double[] back, double[] epsP, double[] work, int e,
                                double dRR, double dZZ, double dTT, double dGamRZ,
                                double lambda, double mu, Flow flow, double strainDt,
                                Damage damage) {
        final int s = COMPONENTS * e;

        // Elastic trial. This is the whole step if the surface is not reached.
        final double lt = lambda * (dRR + dZZ + dTT);
        final double tRR = sig[s] + lt + 2.0 * mu * dRR;
        final double tZZ = sig[s + 1] + lt + 2.0 * mu * dZZ;
        final double tTT = sig[s + 2] + lt + 2.0 * mu * dTT;
        final double tRZ = sig[s + 3] + mu * dGamRZ;

        // Relative stress: how far the deviator sits from the centre of the yield surface.
        final double p = (tRR + tZZ + tTT) / 3.0;
        final double xRR = tRR - p - back[s];
        final double xZZ = tZZ - p - back[s + 1];
        final double xTT = tTT - p - back[s + 2];
        final double xRZ = tRZ - back[s + 3];

        final double norm = Math.sqrt(xRR * xRR + xZZ * xZZ + xTT * xTT + 2.0 * xRZ * xRZ);

        // Temperature is a function of the plastic work already stored at this point, not a
        // field of its own: a Taylor impact is over in tens of microseconds and the thermal
        // diffusion length in steel over that time is well under one element, so nothing
        // conducts anywhere and each point heats only itself. Taken at the start of the step
        // rather than solved for alongside dgamma -- over a step of tens of nanoseconds the
        // temperature moves by millikelvin, and making it implicit would couple the scalar
        // solve to the energy update for no accuracy worth having.
        final JohnsonCook jc = flow.jc();
        final double temperature = jc == null ? 0.0
                : jc.roomTemperature() + flow.thermalCoefficient() * work[e];

        // The surface the elastic trial is tested against. For Johnson-Cook this is the
        // quasi-static surface, at the reference rate, which is deliberate: it is the same
        // surface the solve's residual starts from, so "yielding" and "the solve has a
        // positive root" mean exactly the same thing and cannot disagree. Rate hardening
        // then raises the surface as the material flows, which is what limits how much it
        // flows.
        //
        // A softening point overrides all of that with its frozen flow stress. The floor at
        // zero is what stops a fully damaged point from being handed a negative yield radius:
        // past that strain it has no deviatoric strength left and carries pressure only.
        final boolean softening = damage != null && damage.hasBegun(e);
        double hardening = 0.0;
        double flowStress;
        if (softening) {
            flowStress = damage.strength[e]
                    + damage.softening[e] * (epsP[e] - damage.failureStrain[e]);
            if (flowStress > 0.0) {
                hardening = damage.softening[e];
            } else {
                flowStress = 0.0;
            }
        } else if (jc == null) {
            flowStress = flow.sy0() + flow.hIso() * epsP[e];
            hardening = flow.hIso();
        } else {
            flowStress = jc.quasiStaticFlowStress(epsP[e], temperature);
        }
        final double radius = SQRT_2_3 * flowStress;

        if (!(norm > radius)) {
            // Elastic, including the elastic material where radius is infinite. Written as a
            // negated > so that a NaN norm fails loudly downstream instead of being swallowed.
            sig[s] = tRR;
            sig[s + 1] = tZZ;
            sig[s + 2] = tTT;
            sig[s + 3] = tRZ;
            return 0.0;
        }

        final double stiffness = 2.0 * mu + TWO_THIRDS * flow.hKin();
        double dGamma;
        if (jc == null || softening) {
            dGamma = (norm - radius) / (stiffness + TWO_THIRDS * hardening);
            // The last increment of a softening branch would otherwise carry the yield radius
            // straight through zero: the consistency condition is linear and has no idea the
            // surface has a floor, so it solves for a negative radius and the return
            // overshoots the origin, reflecting the deviator instead of exhausting it. The
            // stress that comes back is small and has the wrong sign, which is exactly the
            // kind of wrong that no energy audit notices. Land on zero instead, which is the
            // state the material is actually in.
            if (softening && flowStress + hardening * SQRT_2_3 * dGamma < 0.0) {
                dGamma = norm / stiffness;
            }
        } else {
            dGamma = solve(norm, stiffness, epsP[e], temperature, strainDt, jc);
        }
        final double inv = 1.0 / norm;
        final double nRR = xRR * inv, nZZ = xZZ * inv, nTT = xTT * inv, nRZ = xRZ * inv;

        final double pull = 2.0 * mu * dGamma;
        sig[s] = tRR - pull * nRR;
        sig[s + 1] = tZZ - pull * nZZ;
        sig[s + 2] = tTT - pull * nTT;
        sig[s + 3] = tRZ - pull * nRZ;

        // Plastic work per unit volume, sigma : deps^p, taken with the back stress as it was
        // before this increment shifted it. Expanding sigma : n over the return gives
        // ||xi_trial|| - 2 mu dgamma + alpha : n, so no second pass over the deviator is
        // needed. The doubled shear term is the tensor inner product in this storage.
        //
        // Two reasons this is worth carrying rather than inferring later. It is the term that
        // closes a dynamic energy balance -- in an impact almost all the input energy ends up
        // here, so without it the balance is an identity between small residuals. And
        // Johnson-Cook thermal softening needs it per point, because the temperature rise that
        // does the softening is adiabatic: beta * W_p / (rho c_p), with beta the
        // Taylor-Quinney fraction.
        final double alphaDotN = back[s] * nRR + back[s + 1] * nZZ + back[s + 2] * nTT
                + 2.0 * back[s + 3] * nRZ;
        work[e] += dGamma * (norm - pull + alphaDotN);

        final double shift = TWO_THIRDS * flow.hKin() * dGamma;
        back[s] += shift * nRR;
        back[s + 1] += shift * nZZ;
        back[s + 2] += shift * nTT;
        back[s + 3] += shift * nRZ;

        epsP[e] += SQRT_2_3 * dGamma;

        // Onset. Checked after the strain advances rather than before, so the flow stress
        // frozen here is the one the point actually reached rather than the one it had a step
        // earlier. The increment that crosses was taken on the undamaged surface, which makes
        // the very first step of softening late by one step and everything after it exact --
        // the same order of error the whole central-difference scheme runs at.
        if (damage != null && !softening && epsP[e] >= damage.failureStrain[e]) {
            damage.begin(e, jc == null
                            ? flow.sy0() + flow.hIso() * epsP[e]
                            : jc.quasiStaticFlowStress(epsP[e], temperature),
                    nRR, nZZ, nRZ);
        }
        return dGamma;
    }

    /**
     * Solves the consistency condition for the plastic multiplier when the flow stress is
     * Johnson-Cook. Safeguarded Newton: Newton where it behaves, bisection where it does not.
     *
     * <pre>  g(dg) = ||xi|| - (2mu + (2/3)H_kin) dg - sqrt(2/3) sigma_y(eps, epsdot, T)</pre>
     *
     * <p>with {@code eps = eps_p + sqrt(2/3) dg} and {@code epsdot = sqrt(2/3) dg / dt}. The
     * multiplier appears in three places -- the elastic pull-back, the accumulated strain the
     * hardening term reads, and the rate the rate term reads -- which is what makes this a
     * solve rather than a formula.
     *
     * <h2>Why a root always exists, and why it is unique</h2>
     *
     * {@code g} is strictly decreasing in {@code dg}: the linear term falls, and both the
     * strain and rate terms raise {@code sigma_y}. At {@code dg = 0} the rate is zero, so the
     * rate factor sits on its clamp at 1 and {@code g(0) = ||xi|| - sqrt(2/3) sigma_y} at the
     * quasi-static surface, which the caller has already established is positive. At
     * {@code dg = ||xi|| / (2mu + (2/3)H_kin)} -- the multiplier a material with no strength
     * at all would need -- the linear term has consumed the whole of {@code ||xi||} and
     * {@code g = -sqrt(2/3) sigma_y < 0}. So {@code [0, that]} is a genuine bracket, not a
     * guess at one, and the root inside it is unique.
     *
     * <p>That bracket is the whole point. It is what lets this keep the property the exact
     * return had: <b>there is no non-convergence failure mode.</b> Every iteration either
     * takes a Newton step that stays inside the bracket or halves it, so the answer is
     * pinned regardless of conditioning. An unguarded Newton here would not merely be slower
     * -- {@code A + B eps^n} has an <em>infinite</em> hardening slope at zero plastic strain,
     * because {@code n < 1}, so the first yielding increment at any point produces a Newton
     * step of zero and the iteration stalls exactly where every plastic zone begins.
     *
     * <p>Started from the top of the bracket, which is the perfectly-plastic estimate and
     * therefore above the root. Approaching from above keeps every iterate at a plastic
     * strain where the hardening slope is finite.
     */
    private static double solve(double norm, double stiffness, double plasticStrain,
                                double temperature, double dt, JohnsonCook jc) {
        final double thermal = jc.thermalTerm(temperature);
        final double tolerance = SOLVE_TOLERANCE * jc.a();

        double lo = 0.0;
        double hi = norm / stiffness;
        double dGamma = hi;

        for (int iteration = 0; iteration < SOLVE_ITERATIONS; iteration++) {
            final double dEps = SQRT_2_3 * dGamma;
            final double eps = plasticStrain + dEps;
            final double rate = dt > 0.0 ? dEps / dt : 0.0;

            final double strain = jc.strainTerm(eps);
            final double rateFactor = jc.rateTerm(rate);
            final double g = norm - stiffness * dGamma
                    - SQRT_2_3 * strain * rateFactor * thermal;

            if (g > 0.0) lo = dGamma; else hi = dGamma;
            if (Math.abs(g) <= tolerance) return dGamma;

            // dg' of the residual. The rate term contributes S C Theta / dg, because the
            // logarithm's derivative against the rate cancels the dt and the sqrt(2/3)
            // introduced by the rate itself -- so no timestep appears here.
            // The hardening slope comes out of the strain term already computed rather than
            // from a second pow: B n eps^(n-1) is n (S - A) / eps exactly. eps is strictly
            // positive here because dGamma is, being confined to the open bracket below.
            double slope = -stiffness
                    - TWO_THIRDS * jc.strainSlopeFrom(strain, eps) * rateFactor * thermal;
            if (!jc.rateIsClamped(rate)) {
                slope -= SQRT_2_3 * strain * jc.c() * thermal / dGamma;
            }

            double next = dGamma - g / slope;
            // Also catches the infinite slope at zero plastic strain, where the step is zero
            // and Newton would sit still, and any NaN.
            if (!(next > lo && next < hi)) next = 0.5 * (lo + hi);
            dGamma = next;
        }
        // Bisection has halved the bracket on every iteration Newton did not resolve, so this
        // is accurate to the bracket width whether or not the tolerance was met.
        return dGamma;
    }

    /** Von Mises equivalent stress, sqrt(3/2) ||dev(sigma)||. */
    public static double vonMises(double[] sig, int e) {
        final int s = COMPONENTS * e;
        final double p = (sig[s] + sig[s + 1] + sig[s + 2]) / 3.0;
        final double a = sig[s] - p, b = sig[s + 1] - p, c = sig[s + 2] - p, d = sig[s + 3];
        return Math.sqrt(1.5 * (a * a + b * b + c * c + 2.0 * d * d));
    }

    /** Hydrostatic pressure, positive in compression: -tr(sigma)/3. */
    public static double pressure(double[] sig, int e) {
        final int s = COMPONENTS * e;
        return -(sig[s] + sig[s + 1] + sig[s + 2]) / 3.0;
    }

    /**
     * Distance of the relative stress from the yield surface. Zero while flowing, negative
     * while elastic. Never positive after an update: that would mean a stress state outside
     * the surface, which is the one thing the return map exists to prevent.
     */
    public static double yieldFunction(double[] sig, double[] back, double[] epsP, int e,
                                       double sy0, double hIso) {
        final int s = COMPONENTS * e;
        final double p = (sig[s] + sig[s + 1] + sig[s + 2]) / 3.0;
        final double xRR = sig[s] - p - back[s];
        final double xZZ = sig[s + 1] - p - back[s + 1];
        final double xTT = sig[s + 2] - p - back[s + 2];
        final double xRZ = sig[s + 3] - back[s + 3];
        final double norm = Math.sqrt(xRR * xRR + xZZ * xZZ + xTT * xTT + 2.0 * xRZ * xRZ);
        return norm - SQRT_2_3 * (sy0 + hIso * epsP[e]);
    }
}
