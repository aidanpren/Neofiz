package org.neofiz.solver;

/**
 * Quadrature rule for the element internal-force kernel.
 *
 * <p>This is a solver-wide flag rather than a rewrite for the same reason
 * {@link org.neofiz.core.Formulation} is: the two rules must be runnable against the same
 * mesh, the same material and the same validation harness, because the only convincing
 * evidence that the cheap one is trustworthy is that it reproduces the expensive one.
 */
public enum Integration {

    /**
     * Full 2x2 Gauss. Four stress evaluations per element. No spurious modes, but it is
     * four times the work and it shear-locks in bending.
     */
    FULL,

    /**
     * One Gauss point at the element centroid, with Flanagan-Belytschko hourglass control.
     *
     * <p>A quarter of the stress work, and it does not shear-lock. The price is that the
     * single point cannot see the hourglass deformation mode -- the one that alternates
     * sign around the element -- so that mode produces exactly zero internal force and is
     * free to grow without bound. Hourglass control supplies the missing stiffness.
     */
    REDUCED;

    /** Stress evaluations per element per step. */
    public int gaussPoints() {
        return this == FULL ? 4 : 1;
    }
}
