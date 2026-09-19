package org.neofiz.solver;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.EnumSource;
import org.neofiz.core.Formulation;
import org.neofiz.core.Material;
import org.neofiz.mesh.QuadMesh;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Tests for the rigid frictionless anvil.
 *
 * <p>Contact has no closed form to check against, so what can be asserted are the
 * <em>conservation properties</em> and the discretisation behaviour. Those turn out to be
 * sharper than a tolerance on a displacement: the momentum identity below is exact to
 * floating point and it is sensitive to precisely the bookkeeping mistakes that are otherwise
 * invisible, and the convergence rate of the energy loss distinguishes an artefact that goes
 * away under refinement from one that does not.
 */
class ContactTest {

    /** 4340 steel with a yield surface. Linear hardening; see TaylorImpactCase for the fit. */
    private static final Material STEEL = Material.STEEL_4340.yielding(792.0e6, 850.0e6);

    private static final double DIAMETER = 7.595e-3;
    private static final double SPEED = 181.0;

    /**
     * A stubby projectile aimed at an anvil at z = 0. Deliberately short so that the tests
     * run in milliseconds while still doing everything the reference case does: impact,
     * spread, and separation.
     */
    private record Shot(ExplicitSolver solver, RigidWall anvil, QuadMesh mesh, int nr, int nz,
                        double momentum0, double kinetic0) {

        void runTo(double seconds) {
            while (solver.time() < seconds) solver.step();
        }

        /** Lowest current axial coordinate anywhere in the body. */
        double lowestZ() {
            double lo = Double.POSITIVE_INFINITY;
            double[] uz = solver.axialDisplacement();
            for (int n = 0; n < mesh.nodeCount; n++) lo = Math.min(lo, mesh.z[n] + uz[n]);
            return lo;
        }

        /** Largest current radius on the impact face. */
        double faceRadius() {
            double max = 0.0;
            double[] ur = solver.radialDisplacement();
            for (int n : QuadMesh.nearFaceNodes(nr, nz)) max = Math.max(max, mesh.r[n] + ur[n]);
            return max;
        }

        /** |change in axial momentum - anvil impulse| / |initial momentum|. */
        double momentumError() {
            return Math.abs((solver.axialMomentum() - momentum0) - anvil.impulse())
                    / Math.abs(momentum0);
        }
    }

    private static Shot shot(int nr, double length, Integration integration, double speed,
                             double wallZ, double friction) {
        final double radius = 0.5 * DIAMETER;
        final double h = radius / nr;
        final int nz = Math.max(1, (int) Math.round(length / h));

        QuadMesh mesh = QuadMesh.solidCylinder(radius, length, nr, nz);
        ExplicitSolver solver = new ExplicitSolver(mesh, STEEL, Formulation.AXISYMMETRIC,
                integration, Kinematics.FINITE_STRAIN, 1.0, 0.5);
        for (int n : QuadMesh.axisNodes(nr, nz)) solver.fixRadial(n);

        RigidWall anvil = RigidWall.atZ(wallZ, friction);
        solver.setRigidWall(anvil);
        solver.setUniformVelocity(0.0, -speed);

        return new Shot(solver, anvil, mesh, nr, nz,
                solver.axialMomentum(), solver.kineticEnergy());
    }

    private static Shot shot(int nr, double length, Integration integration, double speed,
                             double wallZ) {
        return shot(nr, length, integration, speed, wallZ, 0.0);
    }

    private static Shot shot(int nr, Integration integration) {
        return shot(nr, 10.0e-3, integration, SPEED, 0.0);
    }

    private static Shot shot(int nr, Integration integration, double friction) {
        return shot(nr, 10.0e-3, integration, SPEED, 0.0, friction);
    }

    @ParameterizedTest
    @EnumSource(Integration.class)
    @DisplayName("no node ever crosses the anvil, at any step")
    void nothingPenetrates(Integration integration) {
        Shot s = shot(3, integration);

        // Checked every step rather than at the end. A constraint that lets a node dip below
        // the plane and pushes it back next step would pass a final-state check while having
        // solved a different problem, and the recovery would be invisible in the answer.
        while (s.solver().time() < 40e-6) {
            s.solver().step();
            assertTrue(s.lowestZ() >= 0.0,
                    "a node reached z = " + s.lowestZ() + " at t = " + s.solver().time()
                            + "; non-penetration is meant to be exact, not a tolerance");
        }
        assertTrue(s.anvil().impulse() > 0.0, "the anvil must actually have been hit");
    }

    @ParameterizedTest
    @EnumSource(Integration.class)
    @DisplayName("axial momentum change equals the anvil impulse, to floating point")
    void momentumBalanceIsExact(Integration integration) {
        Shot s = shot(3, integration);
        s.runTo(40e-6);

        // The strong audit, and the reason it is exact: every element's axial internal forces
        // sum to zero, because the shape function derivatives are the gradient of a partition
        // of unity. So nothing inside the mesh can move the body's axial momentum and the
        // anvil's impulse is its entire history.
        //
        // This is the assertion that catches contact bookkeeping errors nothing else does.
        // Splitting arrival into "land with the velocity that reaches the plane" and "stop on
        // the following step" leaves the momentum destroyed by the second half belonging to no
        // impulse; it is silent at first impact, where nodes already touching the plane land
        // with zero velocity, and only appears once the contact patch spreads to nodes
        // arriving from above. It reads as 1e-3 here.
        assertTrue(s.momentumError() < 1e-11,
                "momentum balance broken by " + s.momentumError()
                        + " relative; the anvil is applying an impulse it is not accounting "
                        + "for, or accounting for one it is not applying");
    }

    @Test
    @DisplayName("an anvil that is not touched changes nothing at all")
    void aWallBehindTheBodyIsANoOp() {
        // Fired away from the anvil, which sits below the body and is never reached.
        Shot away = shot(3, 10.0e-3, Integration.REDUCED, -SPEED, -1.0);
        away.runTo(20e-6);

        QuadMesh mesh = QuadMesh.solidCylinder(0.5 * DIAMETER, 10.0e-3, 3, away.nz());
        ExplicitSolver free = new ExplicitSolver(mesh, STEEL, Formulation.AXISYMMETRIC,
                Integration.REDUCED, Kinematics.FINITE_STRAIN, 1.0, 0.5);
        for (int n : QuadMesh.axisNodes(3, away.nz())) free.fixRadial(n);
        free.setUniformVelocity(0.0, SPEED);
        while (free.time() < 20e-6) free.step();

        assertEquals(0.0, away.anvil().impulse(), 0.0,
                "an anvil nothing touched must deliver exactly zero impulse");
        assertEquals(0.0, away.anvil().energyLoss(), 0.0,
                "an anvil nothing touched must cost exactly zero energy");
        assertEquals(0, away.anvil().contactCount(),
                "an anvil nothing touched must report no contact");

        // Not merely "small": the presence of an untouched constraint must not perturb the
        // trajectory by a single bit, or the constraint is doing something when idle.
        for (int n = 0; n < mesh.nodeCount; n++) {
            assertEquals(free.axialDisplacement()[n], away.solver().axialDisplacement()[n], 0.0,
                    "an untouched anvil altered the motion at node " + n);
        }
    }

    @ParameterizedTest
    @EnumSource(Integration.class)
    @DisplayName("contact releases on its own once the material pulls away")
    void contactReleases(Integration integration) {
        Shot s = shot(3, integration);

        // Tracked across the run rather than probed at a chosen instant. This specimen is
        // short, so its wave transit is about 2 us and it has finished deforming and let go
        // well inside 15 us -- any fixed probe time is really an assumption about the
        // specimen's length.
        int peakContact = 0;
        double lastTouch = 0.0;
        while (s.solver().time() < 120e-6) {
            s.solver().step();
            if (s.anvil().contactCount() > 0) {
                peakContact = Math.max(peakContact, s.anvil().contactCount());
                lastTouch = s.solver().time();
            }
        }

        assertTrue(peakContact > 1, "the anvil was barely touched; nothing to release from");
        assertTrue(lastTouch < 100e-6,
                "contact was still live at " + lastTouch * 1e6 + " us, right up to the end of "
                        + "the run; this test cannot tell release from a run that stopped too "
                        + "early");

        // Release is not a criterion here, it is a consequence: a node is held only while the
        // internal force presses it into the plane. A stuck contact set would show up as a
        // body still gripped long after the deformation has finished, which is the anvil
        // developing a tensile grip it does not have.
        assertEquals(0, s.anvil().contactCount(),
                "the specimen is still gripped by the anvil after deformation stopped");
        assertTrue(s.lowestZ() > 0.0,
                "having released, the body must have left the plane; lowest z is "
                        + s.lowestZ());
        assertTrue(s.solver().axialMomentum() > 0.0,
                "the specimen must be travelling away from the anvil after it lets go");
    }

    @Test
    @DisplayName("the energy arrival destroys is first order in element size")
    void arrivalEnergyLossConverges() {
        // The one place this contact costs energy: a node's worth of lumped mass stopped in a
        // single step, where the continuum stops an infinitesimal sliver. It is therefore
        // proportional to the mass of the impact face, which is proportional to element size.
        //
        // The rate is what matters, not the magnitude. Letting the internal force accelerate a
        // resting node for a full dt and then zeroing it also "loses a little energy per
        // step", but summed over O(1/h) steps and O(1/h) contact nodes that version is O(1) --
        // 13 % of the impact energy at every resolution, refinement or not. This test is the
        // one that tells the two apart.
        double coarse = lossFraction(3);
        double fine = lossFraction(6);
        double ratio = fine / coarse;

        assertTrue(ratio > 0.4 && ratio < 0.62,
                "halving the element size should roughly halve the arrival loss; got "
                        + coarse + " then " + fine + ", a ratio of " + ratio
                        + ". A ratio near 1 means resting contact is bleeding energy every "
                        + "step and no mesh will fix it");
    }

    private static double lossFraction(int nr) {
        Shot s = shot(nr, Integration.REDUCED);
        s.runTo(40e-6);
        return s.anvil().energyLoss() / s.kinetic0();
    }

    @ParameterizedTest
    @EnumSource(Integration.class)
    @DisplayName("frictionless: the impact face slides outward freely")
    void theFaceIsFreeToSpread(Integration integration) {
        Shot s = shot(3, integration);
        final double before = s.faceRadius();
        s.runTo(40e-6);

        // The anvil constrains the axial degree of freedom and nothing else. A face that
        // could not slide would not mushroom at all, and the Taylor test would be measuring
        // the wrong boundary condition -- fully stuck and frictionless bracket the real one,
        // so which has been implemented has to be observable.
        assertTrue(s.faceRadius() > before * 1.1,
                "the impact face barely spread: " + before * 1e3 + " mm to "
                        + s.faceRadius() * 1e3 + " mm. A frictionless anvil must leave the "
                        + "radial degree of freedom completely alone");
    }

    // ------------------------------------------------------------------ Coulomb friction

    @Test
    @DisplayName("a zero coefficient is the frictionless path, bit for bit")
    void zeroFrictionChangesNothing() {
        // Friction was threaded through the same branch that previously had none, so the
        // guarantee worth having is not "close" but "identical". Anything else means the
        // frictionless results this project has already reported were quietly restated.
        Shot without = shot(3, Integration.REDUCED);
        Shot zero = shot(3, Integration.REDUCED, 0.0);
        without.runTo(40e-6);
        zero.runTo(40e-6);

        for (int n = 0; n < without.mesh().nodeCount; n++) {
            assertEquals(without.solver().radialDisplacement()[n],
                    zero.solver().radialDisplacement()[n], 0.0, "radial drift at node " + n);
            assertEquals(without.solver().axialDisplacement()[n],
                    zero.solver().axialDisplacement()[n], 0.0, "axial drift at node " + n);
        }
        assertEquals(0.0, zero.anvil().frictionDissipation(), 0.0,
                "a frictionless anvil must dissipate exactly nothing");
    }

    @ParameterizedTest
    @EnumSource(Integration.class)
    @DisplayName("friction restrains the spread and the axial audit survives it")
    void frictionRestrainsSpreading(Integration integration) {
        Shot free = shot(3, integration, 0.0);
        Shot gripped = shot(3, integration, 0.2);
        free.runTo(40e-6);
        gripped.runTo(40e-6);

        assertTrue(gripped.faceRadius() < free.faceRadius(),
                "a frictional face spread to " + gripped.faceRadius() * 1e3 + " mm against the "
                        + "frictionless " + free.faceRadius() * 1e3 + " mm; friction must "
                        + "oppose sliding, not assist it");
        assertTrue(gripped.anvil().frictionDissipation() > 0.0,
                "sliding against friction must dissipate energy");

        // The tangential constraint must not disturb the normal one. These are separate
        // directions and the momentum identity is about the normal one alone, so it has to
        // stay exact with friction switched on.
        assertTrue(gripped.momentumError() < 1e-11,
                "friction broke the axial momentum balance by " + gripped.momentumError());
    }

    @Test
    @DisplayName("an infinite coefficient is the welded limit, and a large one converges to it")
    void weldedIsTheLimitOfLargeFriction() {
        // Infinity goes through the same cone test as every other coefficient rather than
        // down a pinned-degree-of-freedom path beside it, so this is the check that the
        // general arithmetic really does reach its own limit.
        Shot welded = shot(3, Integration.REDUCED, Double.POSITIVE_INFINITY);
        Shot large = shot(3, Integration.REDUCED, 5.0);

        // Tracked across the run, not read at the end: the stuck count is per-step, and by
        // 40 us this specimen has finished deforming and left the anvil, so the final value
        // is zero however well the friction worked.
        int peakStuck = 0;
        while (welded.solver().time() < 40e-6) {
            welded.solver().step();
            peakStuck = Math.max(peakStuck, welded.anvil().stuckCount());
        }
        large.runTo(40e-6);

        assertEquals(welded.faceRadius(), large.faceRadius(), welded.faceRadius() * 1e-3,
                "a large coefficient gave " + large.faceRadius() * 1e3 + " mm against the "
                        + "welded " + welded.faceRadius() * 1e3 + " mm");
        assertEquals(0.0, welded.anvil().frictionDissipation(), 0.0,
                "a welded face never slides, so it cannot dissipate by sliding");
        assertTrue(peakStuck > 1, "the welded face never reported itself stuck");
    }

    @Test
    @DisplayName("frictional dissipation vanishes at both limits and peaks between them")
    void dissipationPeaksInTheMiddle() {
        // The structural check on the whole model, and the one that would catch a sign error
        // or a missing normal force. Dissipation is force times sliding distance, so it must
        // go to zero at mu = 0 where there is no force and again as mu grows large where
        // there is no sliding. A monotonic curve would mean one of the two factors is not
        // being applied.
        double none = dissipation(0.0);
        double middle = dissipation(0.2);
        double welded = dissipation(Double.POSITIVE_INFINITY);

        assertEquals(0.0, none, 0.0);
        assertEquals(0.0, welded, 0.0);
        assertTrue(middle > 0.0, "a partly gripped face must dissipate something");
        assertTrue(middle > dissipation(0.02) && middle > dissipation(2.0),
                "dissipation should peak at moderate friction; got " + dissipation(0.02)
                        + ", " + middle + ", " + dissipation(2.0));
    }

    private static double dissipation(double friction) {
        Shot s = shot(3, Integration.REDUCED, friction);
        s.runTo(40e-6);
        return s.anvil().frictionDissipation();
    }

    @Test
    @DisplayName("the anvil's impulse is the initial momentum once the body is arrested")
    void impulseAccountsForTheWholeApproach() {
        Shot s = shot(4, Integration.REDUCED);
        s.runTo(120e-6);

        // By the end the specimen has stopped and rebounded slightly, so the impulse delivered
        // is a little more than the momentum it arrived with. Bracketing rather than
        // asserting equality: the excess is the rebound, which is physical.
        double arrival = Math.abs(s.momentum0());
        assertTrue(s.anvil().impulse() > arrival,
                "the anvil must have removed at least the momentum the specimen arrived with; "
                        + "impulse " + s.anvil().impulse() + " against " + arrival);
        assertTrue(s.anvil().impulse() < 1.2 * arrival,
                "the anvil delivered " + s.anvil().impulse() + " against an arrival momentum "
                        + "of " + arrival + "; a rigid anvil has no restitution, so anything "
                        + "near twice this would mean it is throwing the specimen back");
    }
}
