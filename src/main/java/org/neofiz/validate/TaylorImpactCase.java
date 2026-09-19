package org.neofiz.validate;

import org.neofiz.core.Formulation;
import org.neofiz.core.JohnsonCook;
import org.neofiz.core.Material;
import org.neofiz.mesh.QuadMesh;
import org.neofiz.solver.ExplicitSolver;
import org.neofiz.solver.Integration;
import org.neofiz.solver.Kinematics;
import org.neofiz.solver.RigidWall;

/**
 * The Taylor cylinder impact test: a solid cylinder fired flat-on into a rigid anvil,
 * measured afterwards. The M1 gate.
 *
 * <h2>Why this case and not another</h2>
 *
 * Every gate up to here had a closed form behind it, and each was passed by getting one thing
 * right. Lame needs a correct stiffness assembly. The elastic-plastic cylinder needs the
 * yield criterion, the flow rule, the return map and equilibrium to be right simultaneously,
 * but it is still a quasi-static problem in which nothing rotates and nothing moves far.
 *
 * <p>This one has no closed form at all, and that is the point: it is the first case whose
 * answer depends on the <em>material model</em> rather than on the assembly. It is also the
 * first that exercises inertia as physics rather than as something to damp away -- there is no
 * relaxation here, no pressure ramp, no settling. A plastic wave runs up the specimen from the
 * impact face and the mushroom forms behind it.
 *
 * <h2>The reference case</h2>
 *
 * 4340 steel, 37.97 mm long and 7.595 mm in diameter, at 181 m/s. The measured mushroom
 * diameter is 9.5 mm; published simulations land at 9.80 to 9.83 mm. The gate is 5 %.
 *
 * <p>The published simulations sitting 3 % above the measurement, consistently, is worth
 * noticing rather than averaging away: it is the width of the band within which a
 * well-posed code and a real specimen agree, and it is a floor on what any tolerance here
 * can mean.
 *
 * <h2>The case reports two outputs and they are not the same kind of number</h2>
 *
 * <b>Final length is a constitutive measure.</b> It is converged -- one part in ten thousand
 * between the two finest meshes, with nothing left for an extrapolation to work on -- agrees
 * between quadrature rules, and barely moves when the interface condition is changed from free
 * sliding to welded, 0.9079 against 0.9107. It is shortening driven by the plastic wave, and it is a
 * clean read on the flow stress.
 *
 * <p><b>Mushroom diameter is an interface measure</b>, and it behaves differently in every
 * respect. It is set by how freely the impact face may slide, and the limits are far apart:
 * 10.13 mm sliding freely against 9.04 mm at {@code mu = 0.20}, mesh-converged, with the same
 * material. The measured 9.5 mm is crossed at {@code mu} of about 0.09.
 *
 * <p>It also converges at <em>first</em> order rather than second, because it is read at the
 * contact edge where the free surface meets the anvil, and that corner is a geometric
 * singularity. So it needs {@link MeshConvergence} rather than a finer mesh -- and that
 * matters for more than tidiness. Fitting the coefficient on an unconverged mesh asks for
 * {@code mu = 0.047}, which is implausibly low for steel on steel; fitting it on the
 * extrapolated value asks for 0.09, which is ordinary. The discretisation error was going
 * straight into the material parameter.
 *
 * <p>So a mushroom diameter quoted without both its friction assumption and its mesh
 * convergence is not a result, and tuning a constitutive model to hit one is fitting the flow
 * stress to a boundary condition. The coefficient is a parameter of the case, with
 * {@link #FRICTIONLESS} and {@link #WELDED} as its two limits, and the sweep across it is the
 * result rather than any single value drawn from it.
 */
public final class TaylorImpactCase {

    /** Reference specimen length, metres. */
    public static final double LENGTH = 37.97e-3;
    /** Reference specimen diameter, metres. */
    public static final double DIAMETER = 7.595e-3;
    /** Reference impact speed, m/s. */
    public static final double IMPACT_SPEED = 181.0;
    /** Measured final mushroom diameter, metres. */
    public static final double MEASURED_MUSHROOM_DIAMETER = 9.5e-3;
    /** Centre of the published simulation band, metres. */
    public static final double PUBLISHED_MUSHROOM_DIAMETER = 9.815e-3;

    /** How much time the deformation needs. Everything is frozen well before this. */
    public static final double DURATION = 80e-6;

    /** Free sliding: the widest mushroom, and the standard idealisation for this test. */
    public static final double FRICTIONLESS = 0.0;
    /** Welded: the narrowest mushroom, and the other end of the bracket. */
    public static final double WELDED = Double.POSITIVE_INFINITY;

    /**
     * Plastic strain the rate-independent surrogate is fitted at, for
     * {@link #surrogate4340()}. A Taylor mushroom at this speed reaches equivalent plastic
     * strains of roughly this order, so a secant through {@code A + B eps^n} here is the
     * least-bad single slope available.
     */
    public static final double SURROGATE_FIT_STRAIN = 0.6;

    public record Result(
            Integration integration,
            double friction,
            boolean rateDependent,
            int elementsAcrossRadius,
            int elementCount,
            double elementSize,
            double impactSpeed,
            double finalLength,
            double lengthRatio,
            double mushroomDiameter,
            double mushroomErrorVsMeasuredPercent,
            double mushroomErrorVsPublishedPercent,
            double maxPlasticStrain,
            double maxTemperature,
            double yieldedFraction,
            int contactNodes,
            double momentumBalanceError,
            double energyBalanceError,
            double plasticFraction,
            double contactEnergyLossFraction,
            double frictionDissipationFraction,
            double hourglassFraction,
            double kineticEnergyFraction,
            double duration,
            int steps,
            double wallClockSeconds) {
    }

    private TaylorImpactCase() {
    }

    /** 4340 steel with the full Johnson-Cook flow stress. What the gate is run against. */
    public static Material johnsonCook4340() {
        return Material.STEEL_4340.withJohnsonCook(JohnsonCook.STEEL_4340);
    }

    /**
     * 4340 steel as a rate-independent linear-hardening solid: the Johnson-Cook strain term
     * {@code A + B eps^n} replaced by a secant through it at {@link #SURROGATE_FIT_STRAIN},
     * with the rate and thermal terms simply absent.
     *
     * <p>Kept after Johnson-Cook landed, not superseded by it. Running both is what turns
     * "Johnson-Cook changed the answer" into a measurement of <em>how much</em>, and on this
     * case the honest answer is: remarkably little, for a reason worth knowing. See the
     * report from {@code ./gradlew run}.
     *
     * <p>The secant rather than the initial slope because {@code eps^0.26} does most of its
     * rising in the first few percent of strain: the tangent at the origin is infinite and
     * the tangent at large strain is nearly flat, so neither describes the range this test
     * spends its time in. Hardening is put entirely in the isotropic term because the Taylor
     * test is monotonic -- nothing reverses, so the kinematic split is unobservable here and
     * choosing it would be asserting something this case cannot check.
     */
    public static Material surrogate4340() {
        JohnsonCook jc = JohnsonCook.STEEL_4340;
        double secant = jc.b() * Math.pow(SURROGATE_FIT_STRAIN, jc.n()) / SURROGATE_FIT_STRAIN;
        return Material.STEEL_4340.yielding(jc.a(), secant);
    }

    /** The reference case, Johnson-Cook, frictionless, at the given radial refinement. */
    public static Result reference(int elementsAcrossRadius, Integration integration) {
        return referenceAt(elementsAcrossRadius, integration, FRICTIONLESS);
    }

    /** The reference case at the other end of the interface bracket. */
    public static Result referenceWelded(int elementsAcrossRadius, Integration integration) {
        return referenceAt(elementsAcrossRadius, integration, WELDED);
    }

    /** The reference case at a given Coulomb friction coefficient. */
    public static Result referenceAt(int elementsAcrossRadius, Integration integration,
                                     double friction) {
        return run(johnsonCook4340(), LENGTH, DIAMETER, IMPACT_SPEED,
                elementsAcrossRadius, integration, friction, DURATION);
    }

    /**
     * Mesh convergence of the mushroom diameter, over three meshes each half the element size
     * of the last.
     *
     * <p>The mushroom needs this and the final length does not, which is the whole point of
     * measuring rather than assuming. Length converges at the second order a bilinear quad
     * promises and is done by the second refinement; the mushroom converges at <em>first</em>
     * order, because it is read at the contact edge where the free surface meets the anvil and
     * that corner is a geometric singularity. No amount of refinement restores the missing
     * order, so the right response is to extrapolate and quote a band rather than to keep
     * buying mesh.
     *
     * @param coarsest elements across the radius on the coarsest of the three meshes
     */
    public static MeshConvergence mushroomConvergence(Integration integration, double friction,
                                                      int coarsest) {
        return MeshConvergence.halving(
                referenceAt(coarsest, integration, friction).mushroomDiameter(),
                referenceAt(2 * coarsest, integration, friction).mushroomDiameter(),
                referenceAt(4 * coarsest, integration, friction).mushroomDiameter());
    }

    /** Mesh convergence of the final length ratio, for contrast with the mushroom. */
    public static MeshConvergence lengthConvergence(Integration integration, double friction,
                                                    int coarsest) {
        return MeshConvergence.halving(
                referenceAt(coarsest, integration, friction).lengthRatio(),
                referenceAt(2 * coarsest, integration, friction).lengthRatio(),
                referenceAt(4 * coarsest, integration, friction).lengthRatio());
    }

    /** The same case with the rate-independent surrogate, for comparison. */
    public static Result referenceSurrogate(int elementsAcrossRadius, Integration integration) {
        return run(surrogate4340(), LENGTH, DIAMETER, IMPACT_SPEED,
                elementsAcrossRadius, integration, FRICTIONLESS, DURATION);
    }

    /**
     * @param elementsAcrossRadius elements from the axis to the free surface; the axial count
     *                             is chosen to keep elements square
     * @param duration             physical time to run, seconds
     */
    public static Result run(Material material, double length, double diameter, double speed,
                             int elementsAcrossRadius, Integration integration, double friction,
                             double duration) {
        if (material.isElastic()) {
            // An elastic cylinder bounces off the anvil with no permanent set at all, so
            // "final length" and "mushroom diameter" are not measurements of anything. The
            // run would succeed and the numbers would be meaningless.
            throw new IllegalArgumentException(
                    "the Taylor test measures permanent deformation; the material must yield");
        }

        final double radius = 0.5 * diameter;
        final int nr = elementsAcrossRadius;
        final double h = radius / nr;
        // Square elements. The mushroom is a shear-dominated region and an element stretched
        // along the axis resolves it badly in exactly the direction the material is flowing.
        final int nz = Math.max(1, (int) Math.round(length / h));

        QuadMesh mesh = QuadMesh.solidCylinder(radius, length, nr, nz);
        ExplicitSolver solver = new ExplicitSolver(mesh, material, Formulation.AXISYMMETRIC,
                integration, Kinematics.FINITE_STRAIN, 1.0, 0.5);

        // u_r = 0 on the axis. A symmetry condition, not a choice: material on the centreline
        // cannot leave it without opening the body along its length.
        for (int n : QuadMesh.axisNodes(nr, nz)) solver.fixRadial(n);

        // The anvil, with its Coulomb cone. The welded limit is an infinite coefficient rather
        // than a pinned degree of freedom, which matters for more than tidiness: a pinned face
        // stays pinned after the specimen has rebounded and left the anvil, so it would go on
        // constraining material that is no longer touching anything.
        RigidWall plate = RigidWall.atZ(0.0, friction);
        solver.setRigidWall(plate);

        // A projectile in flight: every node carries the same velocity. No damping -- this is
        // the first case where inertia is the physics rather than something to settle out.
        solver.setUniformVelocity(0.0, -speed);

        final double momentum0 = solver.axialMomentum();
        final double kinetic0 = solver.kineticEnergy();

        long t0 = System.nanoTime();
        int steps = 0;
        while (solver.time() < duration) {
            solver.step();
            steps++;
        }
        double wallClock = (System.nanoTime() - t0) * 1e-9;

        // Geometry is measured on the deformed body, so the current node positions are what
        // count -- not the reference mesh.
        final double[] ur = solver.radialDisplacement();
        final double[] uz = solver.axialDisplacement();

        // Length is the extent of the body, max minus min -- not the height above the anvil.
        // The specimen separates from the anvil once deformation stops and drifts away from
        // it, so measuring from the plane reports the flight, not the length: it goes on
        // climbing for as long as the run does, and the number quietly becomes a function of
        // the duration rather than of the physics.
        double zMax = Double.NEGATIVE_INFINITY, zMin = Double.POSITIVE_INFINITY;
        for (int n = 0; n < mesh.nodeCount; n++) {
            final double zc = mesh.z[n] + uz[n];
            zMax = Math.max(zMax, zc);
            zMin = Math.min(zMin, zc);
        }
        final double finalLength = zMax - zMin;

        // The widest point of the body, wherever it is. On a freely sliding anvil that is the
        // impact face itself; with the face held it moves just above it, so measuring only the
        // face would report the undeformed diameter and call the stuck limit a specimen that
        // never spread.
        double mushroomRadius = 0.0;
        for (int n = 0; n < mesh.nodeCount; n++) {
            mushroomRadius = Math.max(mushroomRadius, mesh.r[n] + ur[n]);
        }
        final double mushroom = 2.0 * mushroomRadius;

        // Exact to floating point: nothing inside the mesh can change axial momentum, so the
        // whole of its history is the anvil's impulse.
        final double momentumError =
                Math.abs((solver.axialMomentum() - momentum0) - plate.impulse())
                        / Math.abs(momentum0);

        // Where the impact energy went. Everything that leaves the kinetic account has to
        // arrive somewhere nameable, and on this case roughly 98 % of it arrives as plastic
        // work -- which is why this balance cannot be written without it. Note the
        // denominators: every fraction here is against the energy that went in. Hourglass
        // energy over recoverable strain energy, the ratio the quasi-static cylinder gates
        // use, is meaningless on an impact: by the time the specimen has unloaded, strain
        // energy is a small residual rather than the scale of the problem, and that ratio
        // reads about 1.0 on a run whose hourglass control is in fact carrying two parts in
        // a thousand.
        // Sliding friction is a term of its own here, not folded into the arrival loss: one is
        // physics the interface is supposed to dissipate and the other is a discretisation
        // artefact, and an audit that adds them cannot tell a frictional anvil from a leaky
        // constraint.
        final double accounted = solver.plasticDissipation() + solver.kineticEnergy()
                + solver.strainEnergy() + plate.energyLoss() + solver.hourglassEnergy()
                + plate.frictionDissipation();

        return new Result(integration, friction, material.isRateDependent(),
                nr, mesh.elementCount, h, speed,
                finalLength, finalLength / length,
                mushroom,
                100.0 * (mushroom - MEASURED_MUSHROOM_DIAMETER) / MEASURED_MUSHROOM_DIAMETER,
                100.0 * (mushroom - PUBLISHED_MUSHROOM_DIAMETER) / PUBLISHED_MUSHROOM_DIAMETER,
                solver.maxPlasticStrain(),
                material.isRateDependent() ? solver.maxTemperature() : Double.NaN,
                solver.yieldedFraction(),
                plate.contactCount(),
                momentumError,
                accounted / kinetic0 - 1.0,
                solver.plasticDissipation() / kinetic0,
                plate.energyLoss() / kinetic0,
                plate.frictionDissipation() / kinetic0,
                solver.hourglassEnergy() / kinetic0,
                solver.kineticEnergy() / kinetic0,
                solver.time(), steps, wallClock);
    }
}
