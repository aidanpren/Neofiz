package org.neofiz.validate;

import org.neofiz.core.Formulation;
import org.neofiz.core.Material;
import org.neofiz.mesh.QuadMesh;
import org.neofiz.solver.ExplicitSolver;
import org.neofiz.solver.Integration;

/**
 * The M1 bridge case: the same thick-walled cylinder as {@link ThickWallCylinderCase}, now
 * pressurised past yield and compared against {@link ElasticPlastic}.
 *
 * <p>Same geometry, same mesh, same harness, one new ingredient. That is deliberate. When a
 * validation case changes several things at once, a failure tells you nothing about which
 * one broke; here the elastic gates still run alongside, so a regression can be attributed.
 *
 * <p><b>Path dependence changes what a careless ramp costs.</b> In the elastic case an
 * overshoot during loading was harmless -- the material forgot it on the way back down. A
 * plastic material does not forget: any transient that pushes an element past yield leaves
 * permanent strain behind, and the final answer carries it. So the ramp has to be slow
 * enough that inertia never overshoots, and the only way to know it was is the energy audit.
 */
public final class PlasticCylinderCase {

    public record Result(
            Integration integration,
            int elementsThroughWall,
            double elementSize,
            double pressure,
            double pressureOverElasticLimit,
            double frontRadiusFem,
            double frontRadiusExact,
            double maxYieldConditionErrorPercent,
            double maxRadialStressErrorPercent,
            double maxHoopStressErrorPercent,
            double boreDisplacement,
            double maxPlasticStrain,
            double yieldedFraction,
            double kineticOverStrainEnergy,
            int steps,
            double wallClockSeconds) {
    }

    private PlasticCylinderCase() {
    }

    public static Result run(Material material, double innerRadius, double outerRadius,
                             double pressure, int elementsThroughWall, Integration integration) {
        // Above the limit pressure the plastic zone reaches the outer surface and nothing
        // is left to contain it: the cylinder accelerates outward without bound and there
        // is no static answer to converge to. Refusing is the honest response. Left alone
        // the run returns a plausible-looking stress field with a quietly growing kinetic
        // energy, which is the failure mode this whole audit exists to catch.
        double limit = ElasticPlastic.limitPressure(innerRadius, outerRadius,
                material.yieldStress());
        if (pressure >= limit) {
            throw new IllegalArgumentException(String.format(
                    "pressure %.1f MPa is at or above the limit pressure %.1f MPa; the "
                            + "cylinder is fully plastic and has no equilibrium",
                    pressure / 1e6, limit / 1e6));
        }

        final int nr = elementsThroughWall;
        final int nz = 4;
        final double wall = outerRadius - innerRadius;
        final double h = wall / nr;

        QuadMesh mesh = QuadMesh.cylinderWall(innerRadius, outerRadius, h * nz, nr, nz);
        ExplicitSolver solver = new ExplicitSolver(mesh, material, Formulation.AXISYMMETRIC,
                integration, 1.0, 0.5);
        solver.fixAllAxial();   // eps_z = 0, which is what the closed form assumes

        final double rMid = 0.5 * (innerRadius + outerRadius);
        final double omega = material.barWaveSpeed() / rMid;
        final double period = 2.0 * Math.PI / omega;

        // Longer than the elastic case. A yielded cylinder is softer than an elastic one, so
        // it rings more slowly and settles more slowly, and any plasticity caused by
        // overshoot on the way up is permanent.
        solver.setPressureRamp(pressure, 20.0 * period);
        solver.setRelaxationDamping(0.7, omega);
        final int steps = (int) Math.ceil(40.0 * period / solver.timestep());

        long t0 = System.nanoTime();
        solver.run(steps);
        double wallClock = (System.nanoTime() - t0) * 1e-9;

        final double sy = material.yieldStress();
        final double exactFront = ElasticPlastic.frontRadius(pressure, innerRadius, outerRadius, sy);
        final double yieldGap = ElasticPlastic.MISES_PLANE_STRAIN * sy;

        double maxYieldErr = 0.0, maxRadialErr = 0.0, maxHoopErr = 0.0;
        for (int e = 0; e < mesh.elementCount; e++) {
            double rc = solver.centroidRadius(e);
            double[] s = solver.centroidStress(e);

            // Inside the plastic zone the yield condition is an identity, not an
            // approximation: sigma_theta - sigma_r must equal 2*sy/sqrt(3) everywhere.
            //
            // The set to check is the elements that actually yielded, not those the closed
            // form says should have. Scoring by the exact front instead makes the metric
            // depend on the mesh in a way that has nothing to do with accuracy: on a coarse
            // mesh no element sits far enough inside the front to qualify and the check goes
            // vacuous, while on a fine mesh it picks up exactly the elements straddling the
            // boundary, so the reported error grows under refinement while the solution
            // improves.
            if (solver.elementPlasticStrain(e) > 0.0) {
                maxYieldErr = Math.max(maxYieldErr, Math.abs((s[2] - s[0]) - yieldGap));
            }
            maxRadialErr = Math.max(maxRadialErr, Math.abs(s[0]
                    - ElasticPlastic.radialStress(rc, innerRadius, outerRadius, pressure, sy)));
            maxHoopErr = Math.max(maxHoopErr, Math.abs(s[2]
                    - ElasticPlastic.hoopStress(rc, innerRadius, outerRadius, pressure, sy)));
        }

        return new Result(integration, nr, h, pressure,
                pressure / ElasticPlastic.elasticLimitPressure(innerRadius, outerRadius, sy),
                solver.plasticFrontRadius(), exactFront,
                100.0 * maxYieldErr / yieldGap,
                100.0 * maxRadialErr / pressure,
                100.0 * maxHoopErr / yieldGap,
                solver.radialDisplacement()[0],
                solver.maxPlasticStrain(),
                solver.yieldedFraction(),
                solver.kineticEnergy() / solver.strainEnergy(),
                steps, wallClock);
    }
}
