package org.neofiz.solver;

/**
 * Whether the solver measures deformation against the original mesh or the current one.
 *
 * <p>Production codes are always finite strain; the small-strain path is kept for the same
 * reason Abaqus keeps NLGEOM off, and for one more: the small-strain gates have exact closed
 * forms behind them, so being able to run both proves the finite-strain machinery agrees
 * where it must before it is trusted where no closed form exists.
 */
public enum Kinematics {

    /**
     * Geometry frozen at the reference configuration. Strain is the linearised symmetric
     * gradient of displacement, and stress is never rotated.
     *
     * <p>Correct to first order in strain, which is enough for the pressurised cylinder
     * (0.13% strain) and worthless for anything that mushrooms, buckles or spins.
     */
    SMALL_STRAIN,

    /**
     * Updated Lagrangian. The gradient operator is rebuilt on the deformed configuration
     * every step, the pressure load follows the surface it acts on, the stable timestep
     * tracks the shrinking elements, and the stress is co-rotated so the update is objective.
     *
     * <p>The objectivity is the part that is easy to omit and catastrophic to omit. Rigid
     * rotation involves no stretching at all, so a constitutive model driven by the rate of
     * deformation sees nothing happening -- while the stored stress components, which are
     * written in the fixed global frame, go on describing a stress state that has rotated
     * away underneath them. The result is stress manufactured out of pure rotation, growing
     * with total rotation rather than with load, in exactly the spinning, folding regions
     * where nobody is checking.
     */
    FINITE_STRAIN;

    public boolean isFinite() {
        return this == FINITE_STRAIN;
    }
}
