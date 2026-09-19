package org.neofiz.solver;

/**
 * A rigid, frictionless, immovable plane normal to the axis, which the body may not pass
 * through. The anvil of the Taylor impact test.
 *
 * <h2>Kinematic, not penalty</h2>
 *
 * The build plan's eventual contact model is penalty-based with Coulomb friction, because
 * that is what deformable-on-deformable pairs need -- a projectile in a bore, a sealing
 * surface, a bolted joint. For a wall that is rigid and cannot move, the kinematic form is
 * not a simplification of that, it is strictly better, and it is worth being explicit about
 * why rather than reaching for the penalty method out of habit:
 *
 * <ul>
 *   <li><b>Non-penetration is exact.</b> A penalty wall is a stiff spring, so it is always
 *       penetrated -- that is the mechanism by which it generates force. The penetration is
 *       small if the stiffness is large, which is the trade the next point is about.</li>
 *   <li><b>There is no timestep to lose.</b> A penalty spring adds a frequency to the system,
 *       so the stable step becomes {@code min(CFL, 2/omega_penalty)}. Stiff enough to keep
 *       penetration negligible is stiff enough to halve the timestep, and with
 *       <em>no mass scaling</em> available to buy it back, that is a direct doubling of run
 *       cost. The kinematic constraint adds no frequency at all.</li>
 *   <li><b>There is no stiffness to pick.</b> Penalty stiffness is a free parameter that
 *       changes the answer, and the honest way to use one is a sensitivity sweep. This has
 *       none.</li>
 * </ul>
 *
 * <h2>Two states, and why it is not one</h2>
 *
 * A node interacting with the plane is in one of exactly two situations, and they need
 * different treatment. Collapsing them into one is the obvious implementation and it carries
 * an error that does not converge.
 *
 * <p><b>Arrival.</b> A node above the plane whose new velocity would carry it across. The
 * axial velocity is replaced by the one that lands it exactly on the plane,
 * {@code v_z = (plane - z) / dt}, so it never penetrates -- not by a little, not for one
 * step. This is the impact, and it is the one event that removes energy: the node's
 * approach velocity is destroyed in a single step, where the continuum would have destroyed
 * the velocity of an infinitesimal sliver. Each node arrives once, so the cost is
 * {@code m v^2 / 2} per node of the impact face, which is first order in element size and
 * converges away.
 *
 * <p><b>Resting.</b> A node already on the plane with the material still pressing it there.
 * Here the wall does not impose a velocity at all -- it supplies a <em>reaction force</em>
 * that exactly cancels the internal axial force. The node does not accelerate, does not
 * move, and because it does not move the reaction does no work. Nothing is lost.
 *
 * <p>The distinction matters far more than it looks. Treating a resting node as an arrival
 * every step -- letting the internal force accelerate it for a full {@code dt} and then
 * zeroing the result -- discards {@code dt^2 f^2 / 2m} per node per step. That is
 * {@code O(h^2)}, which sounds harmless, until it is summed over the {@code O(1/h)} steps of
 * a fixed physical time and the {@code O(1/h)} nodes of a refined contact face. The total is
 * {@code O(1)}: <b>refining the mesh does not reduce it.</b> On the reference Taylor case it
 * is 13 % of the impact energy at every resolution, silently removed from the plastic work
 * that is supposed to be forming the mushroom.
 *
 * <p><b>Release needs no criterion beyond the sign of the force.</b> The usual trap in
 * kinematic contact is a node that sticks -- held to the surface after the material has
 * begun to pull it away, so the wall develops a tensile grip it does not have. A resting
 * node is held only while the internal force presses it into the plane; the moment that
 * force reverses, the node is integrated normally and leaves. There is no contact set aged
 * across steps and no penetration tolerance to tune.
 *
 * <p><b>Frictionless means the radial velocity is never touched.</b> The impact face slides
 * freely outward, which is the standard idealisation for the Taylor test and the reason the
 * mushroom is allowed to form at all. A fully stuck face would be a different problem with a
 * different answer, and the two bracket the real one.
 *
 * <h2>The two audits</h2>
 *
 * <b>Momentum is the strong one, and it is exact.</b> The axial internal forces of every
 * element sum to zero -- the shape function derivatives do, since they are the gradient of a
 * partition of unity -- so nothing inside the mesh can change the body's total axial
 * momentum. With no pressure load, no gravity and no axial restraint, this wall is the only
 * thing that can, and therefore
 *
 * <pre>  sum(m_i v_i)(t) - sum(m_i v_i)(0)  ==  impulse()</pre>
 *
 * to floating point, at every step. That is a machine-precision identity rather than a
 * tolerance, and it is sensitive to the mistake that matters: an impulse applied to a node
 * whose mass is not the mass used to accelerate it.
 *
 * <p><b>Energy is the honest one, and it is not exact.</b> Zeroing the axial velocity of a
 * node removes {@code m v^2 / 2} from the system, and that energy goes nowhere -- it is not
 * stored, not dissipated by the material, not transferred to the wall, which is immovable and
 * so does no work. In the continuum there is nothing to lose, because the material arriving
 * at the wall in any instant is an infinitesimal sliver; in the discrete body the first row
 * of nodes carries real mass, and it stops in one step. The loss is therefore first order in
 * element size and converges away, which is exactly the kind of numerical prop that gets
 * shipped silently. {@link #energyLoss()} against initial kinetic energy is what stops that,
 * on the same principle as the hourglass energy audit: the question is never whether a prop
 * exists, it is how much of the answer it is carrying.
 */
public final class RigidWall {

    /** Axial position of the plane, metres. Material is admissible at or above it. */
    private final double plane;

    /** Coulomb friction coefficient. Zero slides freely; infinite is welded. */
    private final double friction;

    private double impulse;
    private double energyLoss;
    private double frictionDissipation;
    private int contactCount;
    private int stuckCount;

    private RigidWall(double plane, double friction) {
        if (friction < 0.0 || Double.isNaN(friction)) {
            throw new IllegalArgumentException("friction coefficient must be >= 0");
        }
        this.plane = plane;
        this.friction = friction;
    }

    /**
     * A frictionless wall at {@code z}, with the body on the +z side of it. A node is
     * admissible while {@code z_node >= z}.
     */
    public static RigidWall atZ(double z) {
        return new RigidWall(z, 0.0);
    }

    /**
     * A wall at {@code z} with Coulomb friction.
     *
     * <p>{@code Double.POSITIVE_INFINITY} is the welded limit and is exact rather than
     * approached: an infinite cone is never left, so the tangential freedom is held at zero
     * through the same arithmetic every other coefficient uses. That is the same idiom as an
     * elastic material being a plastic one with infinite yield stress -- the limit case falls
     * out of the general path instead of being special-cased beside it.
     */
    public static RigidWall atZ(double z, double friction) {
        return new RigidWall(z, friction);
    }

    /** Axial position of the plane, metres. */
    public double plane() {
        return plane;
    }

    /** Coulomb friction coefficient. */
    public double friction() {
        return friction;
    }

    boolean hasFriction() {
        return friction > 0.0;
    }

    /** Clears the per-step counts. Accumulated impulse and energy are not reset. */
    void beginStep() {
        contactCount = 0;
        stuckCount = 0;
    }

    /**
     * Whether a node resting on the plane stays held for this step: true while the internal
     * axial force presses it into the wall, false the moment that force reverses, which is
     * the release.
     */
    boolean holds(double forceZ) {
        return forceZ <= 0.0;
    }

    /**
     * Accounts for a node held at rest on the plane. The wall's reaction exactly cancels
     * {@code forceZ}, so the node neither accelerates nor moves and no work is done -- only
     * the impulse is recorded.
     *
     * @param dtEffective the interval the velocity update would have applied the force over,
     *                    which is the centred mean of the current and previous steps
     */
    void hold(double forceZ, double dtEffective) {
        impulse += -forceZ * dtEffective;
        contactCount++;
    }

    /**
     * Whether {@code velocityZ} over a step of {@code dt} would carry a node at
     * {@code zCurrent} across the plane.
     */
    boolean crosses(double zCurrent, double velocityZ, double dt) {
        return zCurrent + velocityZ * dt < plane;
    }

    /**
     * Accounts for a node arriving at the plane. The caller places it exactly on the plane
     * and sets its velocity to zero; this records the impulse that took the whole approach
     * momentum away and the kinetic energy that went with it.
     *
     * <p>Destroying the entire approach velocity is not a simplification: for a wall that is
     * rigid and cannot move there is no restitution in the normal direction, so a node that
     * reaches it stops. What is approximate is the timing -- the continuum stops an
     * infinitesimal sliver of material, this stops a node's worth of lumped mass in one step.
     *
     * @param trialVelocity velocity the internal force produced, which would have penetrated
     * @param nodeMass      lumped mass of the node, for the impulse and energy audits
     */
    void arrive(double trialVelocity, double nodeMass) {
        impulse += -nodeMass * trialVelocity;
        energyLoss += 0.5 * nodeMass * trialVelocity * trialVelocity;
        contactCount++;
    }

    /**
     * Coulomb friction on a node resting against the plane: returns its new tangential
     * velocity, accounting for the energy the friction dissipated.
     *
     * <h2>This is a return map, and saying so is not an analogy</h2>
     *
     * It has exactly the structure of the plasticity kernel next door, for the same reason --
     * a constraint that holds until a threshold and then flows:
     *
     * <ul>
     *   <li><b>Trial.</b> Assume the node sticks, so its new tangential velocity is zero, and
     *       solve for the reaction that would achieve it.</li>
     *   <li><b>Yield surface.</b> The friction cone, {@code |T| <= mu N}. Inside it the trial
     *       stands: static friction is whatever it needs to be, exactly like an elastic step
     *       is whatever the strain says.</li>
     *   <li><b>Return.</b> Outside it, the reaction is scaled back onto the cone --
     *       {@code T = mu N} in the direction the trial pointed, which is the direction
     *       opposing the motion the node would otherwise have -- and the node slides under
     *       the remainder.</li>
     * </ul>
     *
     * <p>The trial reaction is recovered from the velocity update the caller would otherwise
     * have performed, {@code v_new = decay v_old + gain (f + T) / m} with {@code v_new = 0},
     * so the damping terms are carried rather than assumed away and the cone is tested
     * against the force the integrator would actually have applied.
     *
     * <h2>Where the energy goes</h2>
     *
     * <b>Sliding dissipates and that is physical</b> -- {@code T} opposes the motion, so its
     * work is negative and the magnitude is heat. It is accumulated separately from
     * {@link #energyLoss()} because it is not an artefact: a frictional interface is supposed
     * to absorb energy, and an energy balance that omitted it would not close.
     *
     * <p><b>Sticking destroys the tangential kinetic energy</b> of a node that was sliding,
     * in the same way and for the same reason arrival destroys its normal kinetic energy: the
     * transition happens in one step rather than over the instant the continuum would take.
     * It is charged to {@link #energyLoss()}, it happens once per node per sticking event
     * rather than every step, and it converges away with element size.
     *
     * @param forceR   internal tangential force on the node
     * @param forceZ   internal normal force; the wall's reaction is its negative
     * @param velocityR tangential velocity at the start of the step
     * @param decay,gain the integrator's damping coefficients, so the cone sees the real force
     * @param dt       the interval the resulting velocity will move the node over
     */
    double slide(double forceR, double forceZ, double velocityR,
                 double decay, double gain, double nodeMass, double dt) {
        // The wall pushes, so the normal reaction is the negative of the internal force that
        // is pressing into it. holds() has already established this is not positive.
        final double normal = -forceZ;

        // The tangential force that would bring this node exactly to rest this step.
        final double trial = -nodeMass * decay * velocityR / gain - forceR;

        // Infinity is handled before the multiply: mu * N with mu infinite and N zero is NaN,
        // and a welded node momentarily carrying no normal force would silently become free.
        final double limit = Double.isInfinite(friction)
                ? Double.POSITIVE_INFINITY : friction * normal;

        if (Math.abs(trial) <= limit) {
            // Inside the cone: static friction holds. Energy is destroyed only on the step
            // the node actually stops, because after that its velocity is already zero.
            energyLoss += 0.5 * nodeMass * velocityR * velocityR;
            stuckCount++;
            return 0.0;
        }

        // On the cone: the wall supplies all it can, in the direction the trial pointed.
        final double tangential = Math.copySign(limit, trial);
        final double slid = decay * velocityR + gain * (forceR + tangential) / nodeMass;
        frictionDissipation += -tangential * slid * dt;
        return slid;
    }

    /**
     * Net axial impulse this wall has delivered to the body, N s. Positive, since the wall
     * can only push. The reaction on the wall is the negative of it.
     */
    public double impulse() {
        return impulse;
    }

    /**
     * Kinetic energy destroyed by arrival, joules. Resting contact contributes nothing. A
     * discretisation artefact, first order in element size, and the honest measure of what
     * the constraint is costing the answer. See the class comment.
     */
    public double energyLoss() {
        return energyLoss;
    }

    /**
     * Energy dissipated by sliding friction, joules. Physical rather than numerical, so it
     * belongs in the energy balance as a term of its own and not alongside
     * {@link #energyLoss()}.
     */
    public double frictionDissipation() {
        return frictionDissipation;
    }

    /** Nodes in contact during the last step. On the Taylor test, the contact patch. */
    public int contactCount() {
        return contactCount;
    }

    /** Nodes held by static friction during the last step, a subset of the contact patch. */
    public int stuckCount() {
        return stuckCount;
    }
}
