package org.neofiz.solver;

/**
 * Incremental co-rotation of stored tensors, the Hughes-Winget construction.
 *
 * <h2>The problem</h2>
 *
 * Stress components are stored in the fixed global frame, but they describe a state that
 * belongs to the material. When the material rotates, those numbers become wrong even though
 * nothing about the physical stress changed -- a bar under tension, spun through ninety
 * degrees, is still under the same tension, but the component that was sigma_zz is now
 * sigma_rr.
 *
 * <p>A constitutive model driven by the rate of deformation cannot fix this on its own,
 * because rigid rotation produces <em>zero</em> rate of deformation. There is no stretching,
 * so the model is handed a zero strain increment and correctly returns an unchanged stress --
 * unchanged in the global frame, which is precisely the error. Left uncorrected, the stress
 * a spinning element reports drifts with accumulated rotation rather than with load. It grows
 * without any corresponding force, and it does so in folding and spalling regions, which are
 * exactly the places nobody has a closed form to check against.
 *
 * <h2>The fix, and why 2D makes it cheap</h2>
 *
 * Rotate the stored tensors by the incremental material rotation before applying the strain
 * increment. In three dimensions that means a polar decomposition or a Rodrigues formula. In
 * the axisymmetric half-plane the spin has exactly one independent component, so the rotation
 * is a single angle:
 *
 * <pre>  phi = dt * (dv_z/dr - dv_r/dz) / 2</pre>
 *
 * <p>and the hoop direction does not participate at all -- it is normal to the r-z plane, so
 * an in-plane rotation leaves sigma_theta untouched. What would be a nine-component problem
 * in general reduces to one 2x2 similarity transform.
 *
 * <h2>Which objective rate this is</h2>
 *
 * Applying the exact finite rotation over the increment rather than a first-order rate term
 * makes the update <em>algorithmically</em> objective: a rigid rotation of any size, taken in
 * any number of steps, returns the stress tensor exactly rotated and its invariants exactly
 * unchanged. That is asserted to machine precision in the tests.
 *
 * <p>This sidesteps the well-known pathology of the naive Jaumann rate, which in simple shear
 * makes the shear stress <em>oscillate sinusoidally</em> with increasing strain -- a material
 * that periodically unloads itself as you keep shearing it, which no metal does. The
 * oscillation only appears past shear strains of order one, which is late enough to survive
 * every small test and arrive during the run that matters.
 */
public final class Corotational {

    private Corotational() {
    }

    /**
     * Rotates one stored tensor in place through angle phi in the r-z plane.
     *
     * @param t     tensor storage, stride 4: {rr, zz, tt, rz}
     * @param index tensor index, not the array offset
     */
    public static void rotate(double[] t, int index, double cos, double sin) {
        final int s = J2.COMPONENTS * index;
        final double rr = t[s], zz = t[s + 1], rz = t[s + 3];

        final double cc = cos * cos, ss = sin * sin, cs = cos * sin;

        t[s] = cc * rr - 2.0 * cs * rz + ss * zz;
        t[s + 1] = ss * rr + 2.0 * cs * rz + cc * zz;
        t[s + 3] = cs * (rr - zz) + (cc - ss) * rz;
        // t[s + 2], the hoop component, is normal to the plane of rotation and invariant.
    }

    /**
     * Cosine of the incremental rotation, given the spin times the interval.
     *
     * <p>This is the Cayley transform, {@code (I - W dt/2)^-1 (I + W dt/2)}, not the
     * exponential. In the plane it reduces to the tangent half-angle identity, so with
     * {@code h = spinDt/2} the rotation is exactly {@code 2 arctan(h)}: no trigonometry, two
     * multiplies and a divide.
     *
     * <p>The choice of Cayley over exponential is not cosmetic and it is not an
     * approximation to it. Paired with a velocity gradient evaluated at the <em>mid-step</em>
     * configuration, this is the rotation that makes rigid body motion produce exactly zero
     * rate of deformation -- the Cayley transform of a rotation is precisely skew, so the
     * symmetric part vanishes identically rather than to some order in the step. Use the
     * exponential instead and rigid rotation leaks a spurious strain of order phi squared per
     * step, which is small, one-signed, and accumulates over a long spin into stress that was
     * never applied.
     */
    public static double cos(double spinDt) {
        final double h = 0.5 * spinDt;
        return (1.0 - h * h) / (1.0 + h * h);
    }

    /** Sine of the incremental rotation. See {@link #cos(double)}. */
    public static double sin(double spinDt) {
        final double h = 0.5 * spinDt;
        return 2.0 * h / (1.0 + h * h);
    }
}
