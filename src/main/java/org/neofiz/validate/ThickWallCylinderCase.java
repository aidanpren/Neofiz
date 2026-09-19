package org.neofiz.validate;

import org.neofiz.core.Formulation;
import org.neofiz.core.Material;
import org.neofiz.mesh.QuadMesh;
import org.neofiz.solver.ExplicitSolver;
import org.neofiz.solver.Integration;

/**
 * The M0 validation case: a thick-walled cylinder under internal pressure, relaxed to a
 * static answer and compared against Lame.
 *
 * <p>Lame is a static solution and this is an explicit dynamic code, so the two are only
 * comparable once the dynamics have settled. Getting there is dynamic relaxation: ramp the
 * pressure over several breathing periods with a smoothstep so inertia barely participates,
 * then hold with mass-proportional damping until the ringing dies. The ratio of kinetic to
 * strain energy is reported so that "it settled" is an audited claim rather than an
 * assumption -- if that number is not small, the comparison below means nothing.
 *
 * <p>Both end conditions are run, because they are two different exact answers from the
 * same solver and the same mesh, separated only by which degrees of freedom are pinned.
 * Passing both is a far stronger statement than passing either: the Lame <em>stresses</em>
 * are identical in the two cases, so only the displacements distinguish them.
 */
public final class ThickWallCylinderCase {

    /** Which axial restraint the slice carries, and therefore which closed form applies. */
    public enum EndCondition {
        /** Ends fully restrained; eps_z = 0 everywhere. */
        RESTRAINED,
        /** Ends open and unloaded; sigma_z = 0 everywhere. */
        OPEN
    }

    public record Result(
            EndCondition endCondition,
            Integration integration,
            int elementsThroughWall,
            double elementSize,
            int elementCount,
            double timestep,
            int steps,
            double boreDisplacementFem,
            double boreDisplacementExact,
            double maxDisplacementErrorPercent,
            double maxHoopStressErrorPercent,
            double kineticOverStrainEnergy,
            double hourglassOverStrainEnergy,
            double wallClockSeconds) {
    }

    private ThickWallCylinderCase() {
    }

    public static Result run(Material material, double innerRadius, double outerRadius,
                             double pressure, int elementsThroughWall, EndCondition end) {
        return run(material, innerRadius, outerRadius, pressure, elementsThroughWall, end,
                Integration.FULL);
    }

    public static Result run(Material material, double innerRadius, double outerRadius,
                             double pressure, int elementsThroughWall, EndCondition end,
                             Integration integration) {
        final int nr = elementsThroughWall;
        final int nz = 4;
        final double wall = outerRadius - innerRadius;
        final double h = wall / nr;
        final double height = h * nz;                 // square elements

        QuadMesh mesh = QuadMesh.cylinderWall(innerRadius, outerRadius, height, nr, nz);
        ExplicitSolver solver = new ExplicitSolver(mesh, material, Formulation.AXISYMMETRIC,
                integration, 1.0, 0.5);

        if (end == EndCondition.RESTRAINED) {
            // u_z identically zero reduces the axisymmetric formulation to plane strain.
            solver.fixAllAxial();
        } else {
            // Pin only the z = 0 row. That removes the axial rigid-body mode without
            // introducing an end effect: the exact open-ended solution has a uniform
            // eps_z, so u_z = eps_z * z, which is already zero on that row.
            for (int i = 0; i <= nr; i++) solver.fixAxial(i);
        }

        final double rMid = 0.5 * (innerRadius + outerRadius);

        // The breathing mode of a ring is extensional, not dilatational: a hoop fibre
        // stretches, it does not compress. So the frequency is the bar wave speed over the
        // mean radius. Timing this off the dilatational speed instead is right to within a
        // factor for ordinary metals and catastrophically wrong as nu approaches 0.5, where
        // that speed runs away and the ramp collapses to a fraction of a real period. The
        // KE/SE audit is what exposes it; without that number the run just returns a wrong
        // answer with every other check still green.
        final double omega = material.barWaveSpeed() / rMid;
        final double breathingPeriod = 2.0 * Math.PI / omega;

        solver.setPressureRamp(pressure, 8.0 * breathingPeriod);
        solver.setRelaxationDamping(0.5, omega);

        final int steps = (int) Math.ceil(13.0 * breathingPeriod / solver.timestep());

        long t0 = System.nanoTime();
        solver.run(steps);
        double wallClock = (System.nanoTime() - t0) * 1e-9;

        // Displacement error over every node. The answer is a function of r alone, so
        // sweeping all nodes also checks that the solution is z-independent.
        double[] ur = solver.radialDisplacement();
        double maxDispErr = 0.0;
        double peak = 0.0;
        for (int i = 0; i < mesh.nodeCount; i++) {
            double exact = exactDisplacement(mesh.r[i], innerRadius, outerRadius, pressure,
                    material, end);
            peak = Math.max(peak, Math.abs(exact));
            maxDispErr = Math.max(maxDispErr, Math.abs(ur[i] - exact));
        }
        double maxDispErrPct = 100.0 * maxDispErr / peak;

        // Hoop stress at element centroids, the optimal stress point for a bilinear quad.
        double maxHoopErr = 0.0;
        double peakHoop = Lame.hoopStress(innerRadius, innerRadius, outerRadius, pressure);
        for (int e = 0; e < mesh.elementCount; e++) {
            double rc = solver.centroidRadius(e);
            double exact = Lame.hoopStress(rc, innerRadius, outerRadius, pressure);
            maxHoopErr = Math.max(maxHoopErr, Math.abs(solver.centroidStress(e)[2] - exact));
        }
        double maxHoopErrPct = 100.0 * maxHoopErr / peakHoop;

        double boreExact = exactDisplacement(innerRadius, innerRadius, outerRadius, pressure,
                material, end);

        double se = solver.strainEnergy();
        double keRatio = se > 0 ? solver.kineticEnergy() / se : Double.NaN;
        double hgRatio = se > 0 ? solver.hourglassEnergy() / se : Double.NaN;

        return new Result(end, integration, nr, h, mesh.elementCount, solver.timestep(), steps,
                ur[0], boreExact, maxDispErrPct, maxHoopErrPct, keRatio, hgRatio, wallClock);
    }

    private static double exactDisplacement(double r, double a, double b, double p,
                                            Material m, EndCondition end) {
        return end == EndCondition.RESTRAINED
                ? Lame.radialDisplacementPlaneStrain(r, a, b, p, m.youngsModulus(), m.poissonRatio())
                : Lame.radialDisplacementPlaneStress(r, a, b, p, m.youngsModulus(), m.poissonRatio());
    }
}
