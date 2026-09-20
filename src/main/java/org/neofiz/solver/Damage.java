package org.neofiz.solver;

import org.neofiz.core.DefectField;

/**
 * Progressive softening, regularised by fracture energy per unit crack area.
 *
 * <p>This is the half of M2 that lets the solver read {@link DefectField}. The field gives
 * every element a failure strain drawn from the weakest-link distribution; this decides what
 * happens once accumulated plastic strain reaches it. The charter forbids binary failure, so
 * nothing is deleted and nothing switches off: the flow stress falls linearly to zero and the
 * element carries less and less of the load until its neighbours are carrying all of it.
 *
 * <h2>Damage is negative hardening</h2>
 *
 * At onset the flow stress is <b>frozen</b> at whatever value the material had reached, and
 * from then on the yield radius is
 *
 * <pre>  sigma_y = sigma_f + H_soft * (eps_p - eps_f)      H_soft &lt; 0, floored at zero</pre>
 *
 * which is the linear hardening the return map already knows how to solve, with the sign
 * turned around. That is not a presentational trick, it is the whole implementation: there is
 * no damage variable in the state and no second constitutive branch. {@link #at} reconstructs
 * D for reporting, and that is the only place it exists.
 *
 * <p>Freezing buys three things at once. The softening branch becomes <b>exactly</b>
 * triangular, so the energy it dissipates has a closed form rather than an accumulator; the
 * return map stays closed-form even for Johnson-Cook, because a frozen flow stress removes
 * the nonlinearity that {@link J2#solve} exists to handle; and the fracture energy stops
 * depending on strain rate and temperature, which is the point -- G_f is meant to be a
 * material constant, and parameterising the softening branch by anything else would put the
 * rate sensitivity back into the one quantity the regularisation is trying to preserve. The
 * softening branch is a cohesive law, not a flow law. It has one job.
 *
 * <h2>Why the slope has to depend on the element size</h2>
 *
 * Softening makes the tangent negative, and a negative tangent localises: the first element to
 * soften unloads its neighbours, so they stop deforming and it takes everything. In an
 * explicit solve the band is therefore exactly one element wide, always, whatever the mesh.
 * If the softening slope were a material constant, the energy dissipated per unit
 * <em>volume</em> would be fixed, the band volume would fall with the mesh, and the energy it
 * takes to break the part would go to zero as the mesh refined. Refining the mesh would make
 * the tube weaker without bound -- the same shape of failure as the one
 * {@link DefectField#effectiveVolume} exists to stop, pointing the other way.
 *
 * <p>The fix is Hillerborg's: fix the energy per unit <b>crack area</b>, not per unit volume,
 * and let the slope follow the element.
 *
 * <pre>  g_f = sigma_f * (eps_u - eps_f) / 2        energy per unit volume, the triangle
 *  G_f = g_f * h                             energy per unit crack area
 *  =&gt;  H_soft = -sigma_f^2 * h / (2 G_f)</pre>
 *
 * <h2>h is not the square root of the area</h2>
 *
 * h is the width of the band, which is the element's extent <em>normal to the crack</em>. The
 * common shortcut of taking {@code sqrt(A)} is that width only for a square element, and it
 * drifts as the square root of the aspect ratio otherwise -- which is not a corner case here,
 * because a tube wall is meshed thin in r and long in z and 10:1 is ordinary. Getting it wrong
 * by a factor of three in the flattering direction would be indistinguishable from the
 * regularisation working.
 *
 * <p>So the band orientation is taken from the physics. At onset the plastic flow direction
 * {@code n} is known, and its in-plane part has a largest principal direction: that is where
 * the material is stretching hardest and therefore where it will open. Writing t for the unit
 * vector along the band and u, v for the element's two mean edge vectors, the crack line
 * through the element has length {@code L = |u.t| + |v.t|} -- exact for any parallelogram --
 * and
 *
 * <pre>  h = A / L</pre>
 *
 * For a band across a rectangle this returns the other side exactly, for a diagonal band
 * across a square it returns {@code a/sqrt(2)}, and it never returns the support width, which
 * is the mistake that looks right until the band is at 45 degrees.
 *
 * <h2>When the in-plane direction is genuinely undetermined</h2>
 *
 * A thin open-ended tube under internal pressure has {@code sigma_z = sigma_r = 0} and all the
 * load in the hoop. The flow direction is then {@code (-1/2, -1/2, 1)/||.||}: its in-plane
 * part is isotropic, and there is no preferred direction in the meridian plane at all. That is
 * not a numerical degeneracy, it is the model telling the truth -- the crack that wants to
 * form is a longitudinal split, and an axisymmetric formulation has no such thing. The ring
 * just thins. When the in-plane anisotropy falls below 1e-9 of the flow direction's norm,
 * {@code sqrt(A)} is used, because with no orientation to prefer the isotropic answer is the
 * only defensible one.
 *
 * <h2>The 2 pi r cancels, exactly</h2>
 *
 * An element here is a ring, so it is worth checking that any of this is still a length. A
 * ring element of meridian area A at radius r has volume {@code 2 pi r A}. A crack through it
 * is a full circumferential ring -- axisymmetry admits no other kind -- whose area is
 * {@code 2 pi r L}. So
 *
 * <pre>  G_f = g_f * (2 pi r A) / (2 pi r L) = g_f * A / L = g_f * h</pre>
 *
 * and the radius drops out on both sides. A ring at r = 100 mm and a ring at r = 1 mm with the
 * same cross-section take the same energy per unit area to break, which is what a material
 * constant has to mean.
 *
 * <h2>The band is one element wide, not one Gauss point wide</h2>
 *
 * Under full integration each point owns a quarter of the element's volume but is handed the
 * <b>whole element's</b> area and edge vectors, so all four compute the same h. Giving each
 * point its own quarter-sized band would be the obvious thing and it is wrong: a bilinear quad
 * cannot represent a strain discontinuity in its interior, so the narrowest band the mesh can
 * actually resolve is one element across however many points sample it. Sizing the band off
 * the point instead would let one element host two parallel cracks in each direction and
 * dissipate twice the energy, which is a quadrature rule changing the fracture toughness.
 *
 * <p>Because the point volumes still sum to the element volume, the total works out to
 * {@code G_f * V_element / h_element} either way, and the full and reduced rules agree on the
 * energy exactly rather than approximately.
 *
 * <h2>The size limit, and what is done about it</h2>
 *
 * There is a largest element this can work on. The element's total response -- elastic
 * unloading plus softening -- has to stay monotone, or the softening branch snaps back and the
 * element cannot be driven through it quasi-statically at all:
 *
 * <pre>  eps_u - eps_f &gt;= sigma_f / E    =&gt;    h &lt;= 2 E G_f / sigma_f^2</pre>
 *
 * Above that size the regularisation is not merely inaccurate, it is asking for a response the
 * material cannot have. It cannot be checked when the solver is built, because sigma_f is not
 * known until the element gets there and it includes however much the material hardened on the
 * way. So it is checked at onset, the slope is clamped at {@code H_soft = -E}, and the clamp is
 * <em>counted</em>: see {@link #clampedPoints}. A clamped point dissipates more than G_f per
 * unit area and its mesh independence is gone. Reporting the count is the difference between
 * knowing that and not.
 *
 * <p>The clamp also happens to be what keeps the return map non-singular. The plastic
 * multiplier divides by {@code 2 mu + (2/3)(H_soft + H_kin)}, which at the clamp is
 * {@code 2 mu - (2/3) E}: positive for every nu &lt; 0.5, and vanishing exactly at nu = 0.5,
 * which {@link org.neofiz.core.Material} already refuses. It is positive but it is not large
 * -- 0.108 E for steel at nu = 0.29, 0.071 E for copper at nu = 0.34 -- so a clamped point is
 * numerically delicate as well as wrong. Both problems have the same cure, which is a mesh
 * fine enough not to clamp.
 *
 * <h2>What this model does not have</h2>
 *
 * <p><b>Triaxiality.</b> The failure strain here is a property of the material and the defect
 * field, and nothing else. Real ductile failure strain falls roughly exponentially with stress
 * triaxiality, and a pressurised wall spans a wide range of it between bore and outside
 * surface, so the bore should fail earlier than this model makes it. That is the
 * Johnson-Cook D1-D5 term and it is the next thing to add; until then the <em>location</em> of
 * failure is set by the defect field alone, where it ought to be set by the defect field and
 * the stress state together.
 *
 * <p><b>Loss of hydrostatic strength.</b> Only the flow stress is degraded, so a fully damaged
 * element has no shear strength but still resists volume change: it behaves as a fluid, not as
 * a crack. For a bursting wall that is enough, because the wall fails by losing its hoop
 * capacity and thinning, and hoop tension is deviatoric. For anything that needs a
 * traction-free surface to open, it is not.
 */
public final class Damage {

    /** Stride of {@link #geometry}: meridian area, then the two mean edge vectors. */
    public static final int GEOMETRY_STRIDE = 5;

    /** Below this in-plane anisotropy the band has no preferred direction. See the class note. */
    private static final double ISOTROPIC = 1.0e-9;

    /**
     * Fracture energy per unit crack area, J/m^2, one per integration point.
     *
     * <p>Per point rather than per body because a mesh may now hold more than one substance,
     * and G_f is the property that differs most between them -- four orders of magnitude from
     * a glass to a tough steel. Broadcasting a single value across the array is the uniform
     * case and costs the same arithmetic.
     */
    private final double[] fractureEnergy;
    /** Young's modulus, for the snap-back clamp. One per point, for the same reason. */
    private final double[] youngsModulus;

    /** Plastic strain at which softening begins, one per integration point. */
    final double[] failureStrain;
    /**
     * Per point, stride {@link #GEOMETRY_STRIDE}: the meridian area of the point's
     * <em>element</em> and that element's two mean edge vectors {@code (uR, uZ)} and
     * {@code (vR, vZ)}, in the reference configuration. Enough to turn a band direction into a
     * band width and nothing more. Held per point rather than per element because that is the
     * indexing everything else in the return map uses.
     */
    private final double[] geometry;

    /**
     * Flow stress frozen at onset, Pa. Zero until the point gets there, which is what
     * {@link #hasBegun} reads -- so "damaged" is a state of the material, not a flag beside it.
     */
    final double[] strength;
    /** Softening modulus, Pa per unit plastic strain, negative. Zero until onset. */
    final double[] softening;
    /** Band width chosen at onset, metres. Zero until then. */
    private final double[] width;

    /**
     * @param fractureEnergy energy per unit crack area, J/m^2. For a ductile steel this is
     *                       tens to hundreds of kJ/m^2; for a brittle ceramic, tens of J/m^2.
     * @param youngsModulus  the material's E, which sets the largest element this can be
     *                       regularised on
     * @param failureStrain  plastic strain at onset, one per integration point
     * @param geometry       stride {@link #GEOMETRY_STRIDE} per point: area, uR, uZ, vR, vZ
     */
    public Damage(double fractureEnergy, double youngsModulus,
                  double[] failureStrain, double[] geometry) {
        this(fill(fractureEnergy, failureStrain.length),
                fill(youngsModulus, failureStrain.length), failureStrain, geometry);
    }

    private static double[] fill(double value, int points) {
        final double[] out = new double[points];
        java.util.Arrays.fill(out, value);
        return out;
    }

    /**
     * The per-point form, for a mesh holding more than one substance.
     *
     * @param fractureEnergy energy per unit crack area, J/m^2, one per integration point
     * @param youngsModulus  the E of the material at each point
     * @param failureStrain  plastic strain at onset, one per integration point
     * @param geometry       stride {@link #GEOMETRY_STRIDE} per point: area, uR, uZ, vR, vZ
     */
    public Damage(double[] fractureEnergy, double[] youngsModulus,
                  double[] failureStrain, double[] geometry) {
        if (fractureEnergy.length != failureStrain.length
                || youngsModulus.length != failureStrain.length) {
            throw new IllegalArgumentException(
                    "fracture energy, Young's modulus and failure strain must agree in length");
        }
        for (int p = 0; p < failureStrain.length; p++) {
            if (!(fractureEnergy[p] > 0.0)) {
                throw new IllegalArgumentException(
                        "fracture energy must be positive, was " + fractureEnergy[p]);
            }
            if (!(youngsModulus[p] > 0.0)) {
                throw new IllegalArgumentException(
                        "Young's modulus must be positive, was " + youngsModulus[p]);
            }
        }
        if (geometry.length != GEOMETRY_STRIDE * failureStrain.length) {
            throw new IllegalArgumentException(
                    "failure strain has " + failureStrain.length + " points but the geometry "
                            + "has room for " + (geometry.length / (double) GEOMETRY_STRIDE));
        }
        for (int p = 0; p < failureStrain.length; p++) {
            if (!(failureStrain[p] > 0.0)) {
                throw new IllegalArgumentException(
                        "failure strain must be positive; point " + p + " has "
                                + failureStrain[p]);
            }
            if (!(geometry[GEOMETRY_STRIDE * p] > 0.0)) {
                throw new IllegalArgumentException(
                        "area must be positive; point " + p + " has "
                                + geometry[GEOMETRY_STRIDE * p]);
            }
        }
        this.fractureEnergy = fractureEnergy.clone();
        this.youngsModulus = youngsModulus.clone();
        this.failureStrain = failureStrain;
        this.geometry = geometry;
        this.strength = new double[failureStrain.length];
        this.softening = new double[failureStrain.length];
        this.width = new double[failureStrain.length];
    }

    /** Fracture energy at one point, J/m^2. */
    public double fractureEnergy(int point) {
        return fractureEnergy[point];
    }

    public int points() {
        return failureStrain.length;
    }

    /** Whether softening has started at a point. */
    public boolean hasBegun(int point) {
        return strength[point] > 0.0;
    }

    /** The band width chosen at onset, metres, or zero if the point has not softened. */
    public double bandWidth(int point) {
        return width[point];
    }

    /** The flow stress frozen at onset, Pa, or zero if the point has not softened. */
    public double strengthAtOnset(int point) {
        return strength[point];
    }

    /** The softening modulus chosen at onset, Pa, negative, or zero before then. */
    public double softeningModulus(int point) {
        return softening[point];
    }

    /**
     * Freezes the flow stress, orients the band, and fixes the softening slope. Called once per
     * point, by {@link J2#update}, on the increment that carries the plastic strain past the
     * failure strain.
     *
     * <p>The increment that crosses is taken on the undamaged surface, so the first step of
     * softening is late by one step. Everything after it is exact.
     *
     * @param nRR in-plane components of the plastic flow direction, which is where the band
     *            orientation comes from. The hoop component is not passed because a
     *            hoop-normal crack is a longitudinal split and this formulation has none.
     */
    void begin(int point, double flowStress, double nRR, double nZZ, double nRZ) {
        final int g = GEOMETRY_STRIDE * point;
        final double area = geometry[g];

        // In-plane principal direction of the flow: the eigenvector of [[nRR, nRZ], [nRZ, nZZ]]
        // with the larger eigenvalue, which is the direction the material is stretching fastest
        // and therefore the normal to the band. t is the band's own direction, at right angles.
        final double diff = nRR - nZZ;
        final double cross = 2.0 * nRZ;
        final double anisotropy = Math.hypot(diff, cross);

        final double h;
        if (anisotropy < ISOTROPIC) {
            h = Math.sqrt(area);
        } else {
            final double phi = 0.5 * Math.atan2(cross, diff);
            final double tR = -Math.sin(phi), tZ = Math.cos(phi);
            final double line = Math.abs(geometry[g + 1] * tR + geometry[g + 2] * tZ)
                    + Math.abs(geometry[g + 3] * tR + geometry[g + 4] * tZ);
            h = area / line;
        }
        width[point] = h;

        strength[point] = flowStress;
        final double slope = -flowStress * flowStress * h / (2.0 * fractureEnergy[point]);
        // Assigned rather than computed through the range so that a clamped point holds
        // exactly -E, which is what makes clampedPoints() an equality test rather than a
        // tolerance.
        final double floor = -youngsModulus[point];
        softening[point] = slope < floor ? floor : slope;
    }

    /**
     * The damage variable, 0 before onset and 1 once the flow stress has reached zero. Derived
     * rather than stored: it is a way of reading the state, not part of it.
     */
    public double at(int point, double plasticStrain) {
        final double sf = strength[point];
        if (sf <= 0.0) return 0.0;
        final double d = -softening[point] * (plasticStrain - failureStrain[point]) / sf;
        return d <= 0.0 ? 0.0 : (d >= 1.0 ? 1.0 : d);
    }

    /**
     * Energy dissipated by softening at a point, per unit volume, J/m^3. Closed form: the area
     * of the triangle traversed so far.
     *
     * <p>At full damage this is exactly {@code G_f / h}, so multiplied by the point's volume it
     * is exactly {@code G_f} times the crack area -- which is the property the whole
     * construction exists to have, and is asserted as an identity rather than a tolerance.
     *
     * <p>This is a <em>closed form</em>, not a running total, so it is the analytic answer
     * rather than what the solver actually took out. The solver's own figure, which lands in
     * the plastic work density alongside the pre-damage plastic work, is a right Riemann sum of
     * a falling stress and so undershoots this by one strain increment out of the softening
     * range. Comparing the two is how the discretisation error is measured instead of assumed.
     */
    public double dissipation(int point, double plasticStrain) {
        final double sf = strength[point];
        if (sf <= 0.0) return 0.0;
        final double slope = softening[point];
        final double full = -sf / slope;                 // excess at which strength reaches zero
        double excess = plasticStrain - failureStrain[point];
        if (excess <= 0.0) return 0.0;
        if (excess > full) excess = full;
        return excess * (sf + 0.5 * slope * excess);
    }

    /** Largest damage anywhere, given the current plastic strains. */
    public double maximum(double[] plasticStrain) {
        double max = 0.0;
        for (int p = 0; p < strength.length; p++) {
            final double d = at(p, plasticStrain[p]);
            if (d > max) max = d;
        }
        return max;
    }

    /** Fraction of points that have started to soften. */
    public double begunFraction() {
        int n = 0;
        for (double s : strength) if (s > 0.0) n++;
        return (double) n / strength.length;
    }

    /** Fraction of points whose flow stress has reached zero. */
    public double failedFraction(double[] plasticStrain) {
        int n = 0;
        for (int p = 0; p < strength.length; p++) {
            if (at(p, plasticStrain[p]) >= 1.0) n++;
        }
        return (double) n / strength.length;
    }

    /**
     * Points whose softening slope hit the snap-back clamp, and whose fracture energy is
     * therefore larger than G_f and mesh-dependent again. Should be zero. If it is not, the
     * mesh is too coarse for this material's fracture energy and
     * {@link #largestRegularisableElement} says by how much.
     */
    public int clampedPoints() {
        int n = 0;
        for (int p = 0; p < softening.length; p++) {
            if (softening[p] == -youngsModulus[p]) n++;
        }
        return n;
    }

    /**
     * The largest band width this fracture energy can regularise at a given flow stress,
     * {@code 2 E G_f / sigma_f^2}. An element wider than this snaps back.
     *
     * <p>Static and taking every argument because it is the number to consult when choosing a
     * mesh, which happens before any of this exists.
     */
    public static double largestRegularisableElement(double fractureEnergy,
                                                     double youngsModulus, double flowStress) {
        return 2.0 * youngsModulus * fractureEnergy / (flowStress * flowStress);
    }

    /**
     * Plastic strain range over which the flow stress falls from {@code flowStress} to zero,
     * for a band of this width. The inverse of the slope formula, for tests and for reporting.
     */
    public static double softeningRange(double fractureEnergy, double flowStress, double width) {
        return 2.0 * fractureEnergy / (flowStress * width);
    }
}
