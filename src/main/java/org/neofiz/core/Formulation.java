package org.neofiz.core;

/**
 * Which two-dimensional idealisation an element uses.
 *
 * <p>All three share the mesher, the time integrator, the constitutive models and the
 * batching. Axisymmetric adds the hoop strain component eps_theta = u_r / r, giving the
 * strain-displacement matrix one extra row, and integrates element volume over 2*pi*r
 * rather than a constant out-of-plane thickness. That is the entire difference.
 *
 * <p>The flag exists from M0 because choosing one formulation and designing around it
 * costs a rewrite. Axisymmetric covers the cannon, pressure vessel, rocket, firearm and
 * hydraulic cylinder; planar covers the truss.
 */
public enum Formulation {
    /** Body of revolution about the z axis. In-plane coordinates are (r, z). */
    AXISYMMETRIC,
    /** Long prismatic body; eps_z = 0. Out-of-plane thickness is nominal. */
    PLANE_STRAIN,
    /** Thin plate loaded in its own plane; sigma_z = 0. */
    PLANE_STRESS;

    public boolean isAxisymmetric() {
        return this == AXISYMMETRIC;
    }

    /** Number of independent strain components carried by the B matrix. */
    public int strainComponents() {
        return this == AXISYMMETRIC ? 4 : 3;
    }
}
