package org.neofiz.validate;

import org.neofiz.core.DefectField;
import org.neofiz.core.Formulation;
import org.neofiz.core.JohnsonCook;
import org.neofiz.core.Material;
import org.neofiz.core.Weibull;
import org.neofiz.mesh.QuadMesh;
import org.neofiz.solver.ExplicitSolver;
import org.neofiz.solver.Integration;
import org.neofiz.solver.Kinematics;

/**
 * The same tube, pressurised twice, bursting at two different pressures in two different
 * places.
 *
 * <p>This is what M2 was for. {@link BurstCase} bursts a tube whose every element is identical,
 * so it has one answer and it is the closed form's; the only interesting question there was how
 * to <em>measure</em> a burst once equilibrium has stopped existing. Here every element has its
 * own failure strain, drawn from {@link DefectField}, and the tube fails where the field is
 * weakest rather than everywhere at once.
 *
 * <h2>What makes a tube long enough to have a "where"</h2>
 *
 * {@link BurstCase} runs on two elements of axial extent, because the reference solution is
 * uniform along the axis and a longer tube would have been a waste of a hundred thousand
 * timesteps. That is exactly the wrong shape for this: a two-element slice has nowhere to
 * localise. So this runs a genuine length of tube -- forty millimetres, thirteen correlation
 * lengths of it -- and asks which axial station gives way.
 *
 * <p>The axial restraint is kept. Every node has {@code eps_z = 0}, which is not a convenience
 * but the stress state of a thin pressurised tube: the 1 : 1/2 : 0 state has no deviatoric
 * axial component, so von Mises flow produces no axial strain whether the ends are held or
 * free. What it does cost is the axial component of necking -- the wall can still bulge
 * unevenly along z, because neighbouring rings may expand by different amounts, but it cannot
 * draw material along the tube into the neck. The localisation this shows is a bulge, not a
 * full neck.
 *
 * <h2>Two mechanisms, and which one wins</h2>
 *
 * A ductile tube has two ways to stop holding pressure. It can reach the <b>plastic
 * instability</b>, where the wall thins faster than the material hardens and the load maximum
 * arrives with the material still perfectly healthy -- that is the whole of M1, and it has no
 * scatter at all, because it is a property of the flow curve. Or it can <b>fail</b>, at a
 * strain drawn from the defect field, which does have scatter.
 *
 * <p>Which arrives first is not a modelling choice, it is a number: compare the mean failure
 * strain to {@link Burst#burstHoopStrain}. Annealed copper's real ductility is several times
 * its instability strain, so the instability wins outright, and the measured answer is that a
 * defect field knocks the burst pressure down by <b>0.6 %</b> with a scatter of <b>0.08 %</b>
 * -- a true answer, a dull one, and the one most engineering metals actually give. The
 * interesting regime needs a failure strain at or below the instability, which is what
 * {@link #SCALE} and the sweeps in the report set up. That is a statement about what is being
 * demonstrated, not about copper.
 *
 * <h2>A pressurised tube is not a chain</h2>
 *
 * The result worth having here was not expected. Weibull weakest-link theory -- all of
 * {@link org.neofiz.core.Weibull}, all of {@link DefectField} -- describes a body that fails
 * when its <em>weakest</em> piece does, and it predicts two things about size: a longer tube is
 * weaker, and its scatter is size-independent, because the minimum of n Weibulls is Weibull
 * with the same modulus.
 *
 * <p>Neither happens. Over an eightfold change in modelled length the mean burst pressure is
 * <b>flat</b> -- 0.9354, 0.9280, 0.9293, 0.9298 of the defect-free value -- while the scatter
 * <b>falls</b>, 0.933 %, 0.742 %, 0.526 %, 0.408 %. That is {@code L^-0.40}, against
 * {@code L^-0.5} for an average over independent patches and {@code L^0} for a weakest link.
 *
 * <p>The reason is {@link #shearLagLength}. A thin shell shares load along its axis over
 * {@code sqrt(R t)}, seven millimetres here, and the defect field's correlation length is
 * three. Every weak patch is shorter than the distance over which its neighbours can carry it,
 * so the tube <em>averages</em> its defects instead of failing at them. Pushing the correlation
 * length up through {@code sqrt(R t)} turns the behaviour over exactly as that argument
 * predicts: on an 80 mm tube, at 2, 4, 8, 16 and 32 mm the tube holds 0.9373, 0.9262, 0.9140,
 * 0.9019 and 0.8905 of the defect-free pressure, monotonically weaker, with the knee where the
 * shell theory says to put it. The scatter climbs with it, 0.46 %, 0.60 %, 0.88 %, 1.59 % and
 * then 1.33 % -- the last figure is five draws from a tube containing two and a half patches
 * and is not worth reading as a trend. The mean is.
 *
 * <p>So the weakest-link apparatus is right about the material and describes only part of the
 * structure it is put into, and which part depends on a length scale that belongs to the
 * <em>geometry</em> rather than to the statistics. Along its axis this tube is closer to a
 * bundle of parallel rings than to a chain. Through its wall it is a chain, which is why
 * refining the mesh through the wall changes nothing at all.
 */
public final class DefectBurst {

    /** Mean radius, metres. Same tube as {@link BurstCase}. */
    public static final double MEAN_RADIUS = BurstCase.MEAN_RADIUS;
    /** Mean diameter over wall thickness. */
    public static final double SLENDERNESS = BurstCase.REFERENCE_SLENDERNESS;
    /** Elements through the wall. */
    public static final int THROUGH_WALL = 4;
    /** Modelled length of tube, metres. */
    public static final double TUBE_LENGTH = 40.0e-3;
    /** Elements along the axis. */
    public static final int ALONG_AXIS = 40;

    /** Correlation length of the defect field, metres. Thirteen of them fit in the tube. */
    public static final double CORRELATION = 3.0e-3;

    /**
     * Weibull modulus. Twelve is a wrought metal with real but not dramatic scatter: a
     * coefficient of variation of 10 % in the failure strain.
     */
    public static final double MODULUS = 12.0;

    /**
     * Weibull scale: the failure strain of a tube of {@link #REFERENCE_VOLUME}, which is the
     * nominal tube below. Quoted that way so it is directly comparable to
     * {@link Burst#burstHoopStrain} rather than being a per-unit-volume figure that has to be
     * scaled before it means anything.
     *
     * <p>Set well below the instability, so that failure controls and the run has something to
     * say. See the class note on why that is a choice about what to demonstrate rather than a
     * property of copper.
     */
    public static final double SCALE = 0.04;

    /**
     * Fracture energy per unit crack area, J/m^2.
     *
     * <p>Twelve kilojoules is an order below annealed copper's real toughness, and is set with
     * {@link #SCALE} rather than independently of it: a material that fails at 4 % strain is
     * not one that tears at 120 kJ/m^2, and pairing a low failure strain with a high tearing
     * energy would produce a tube that starts failing early and then refuses to finish.
     */
    public static final double FRACTURE_ENERGY = 12.0e3;

    /**
     * The volume the Weibull scale refers to: the ring volume of the nominal tube above.
     *
     * <p>A constant rather than whatever mesh happens to be running, so that changing the
     * modelled length changes the number of independent draws and <em>nothing else</em>. Tying
     * it to the mesh would have made a longer tube weaker per element as well as more
     * numerous, and the two effects are exactly the pair this is meant to separate.
     */
    public static final double REFERENCE_VOLUME =
            2.0 * Math.PI * MEAN_RADIUS * (2.0 * MEAN_RADIUS / SLENDERNESS) * TUBE_LENGTH;

    /**
     * Ring-breathing periods to keep running after the peak, so that the tube has time to
     * choose a station. Stops early the moment anything reaches full damage.
     */
    public static final double POST_PEAK_PERIODS = 60.0;

    private static final double CFL_SAFETY = 0.5;

    /** One firing of one tube. */
    public record Shot(long seed, boolean traversed,
                       double peakCapacity, double peakStrain, double peakApplied,
                       double weakestStrain, double meanStrain,
                       double failureZ, double maxDamage, double damagedFraction,
                       int clampedPoints, double fractureEnergy, double plasticEnergy,
                       double bulgeRatio, int steps, double wallClock) {

        /** Peak capacity as a fraction of the closed-form, defect-free burst pressure. */
        public double knockdown(double reference) {
            return peakCapacity / reference;
        }
    }

    /**
     * A population of nominally identical tubes, which is the only thing a defect model can
     * really be asked about. One tube's burst pressure is a number; a hundred tubes' is a
     * distribution, and the distribution is the claim.
     */
    public record Population(Setup setup, int shots, double meanCapacity, double sdCapacity,
                             double meanRatio, double locationSpread, double meanBulge,
                             double meanDamaged, int notTraversed) {

        /** Coefficient of variation of the burst pressure, per cent. */
        public double scatterPercent() {
            return 100.0 * sdCapacity / meanCapacity;
        }

        /**
         * Spread of the failure location as a fraction of the tube length. A uniform
         * distribution over the whole tube gives {@code 1/sqrt(12) = 0.2887}; anything much
         * below that means the tube is picking the same place every time, which would mean the
         * geometry is choosing rather than the field.
         */
        public double locationSpreadFraction() {
            return locationSpread / setup.length();
        }
    }

    private DefectBurst() {
    }

    /** The closed-form burst pressure of the same tube with no defects anywhere, Pa. */
    public static double referencePressure() {
        return referencePressure(SLENDERNESS);
    }

    /** The closed-form burst pressure of a defect-free tube of the given slenderness, Pa. */
    public static double referencePressure(double slenderness) {
        return Burst.burstPressure(BurstCase.COPPER.johnsonCook(),
                2.0 * MEAN_RADIUS / slenderness, MEAN_RADIUS);
    }

    /**
     * The shear-lag length of this shell, {@code sqrt(R t)}, metres.
     *
     * <p>The distance over which a thin shell shares load along its axis. A weak patch shorter
     * than this is carried by its neighbours and barely shows up in the burst pressure; a patch
     * longer than this stands alone. It is therefore the length the defect field's correlation
     * length has to be compared against, and the comparison decides whether the tube behaves as
     * a weakest-link chain or as a bundle of parallel rings. See the class note.
     */
    public static double shearLagLength() {
        return shearLagLength(SLENDERNESS);
    }

    /** The shear-lag length of a tube of the given slenderness, metres. */
    public static double shearLagLength(double slenderness) {
        return Math.sqrt(MEAN_RADIUS * 2.0 * MEAN_RADIUS / slenderness);
    }

    /** Fires {@code shots} nominally identical tubes, seeded 1..shots, and pools the result. */
    public static Population volley(int shots, Setup setup) {
        final double[] capacity = new double[shots];
        final double[] where = new double[shots];
        double bulge = 0.0, damaged = 0.0;
        int notTraversed = 0;

        for (int s = 0; s < shots; s++) {
            final Shot shot = fire(s + 1, setup);
            capacity[s] = shot.peakCapacity();
            where[s] = shot.failureZ();
            bulge += shot.bulgeRatio();
            damaged += shot.damagedFraction();
            if (!shot.traversed()) notTraversed++;
        }

        final double mean = mean(capacity);
        return new Population(setup, shots, mean, sd(capacity, mean),
                mean / referencePressure(setup.slenderness()), sd(where, mean(where)),
                bulge / shots, damaged / shots, notTraversed);
    }

    private static double mean(double[] a) {
        double m = 0.0;
        for (double x : a) m += x;
        return m / a.length;
    }

    private static double sd(double[] a, double mean) {
        if (a.length < 2) return 0.0;
        double v = 0.0;
        for (double x : a) v += (x - mean) * (x - mean);
        return Math.sqrt(v / (a.length - 1));
    }

    /**
     * Everything about the tube and its defects except which realisation of them. Bundled so
     * that a sweep varies one named thing and carries the rest forward unchanged, rather than
     * threading six positional doubles through every call and hoping.
     *
     * @param failureScale    Weibull scale at {@link #REFERENCE_VOLUME}, or zero for a
     *                        defect-free tube, which is how the same code produces the control
     * @param fractureEnergy  J/m^2
     * @param correlation     correlation length of the defect field, metres
     * @param length          modelled length of tube, metres
     * @param slenderness     mean diameter over wall thickness; the wall is the only
     *                        dimension a design map varies, and the mean radius is held at
     *                        {@link #MEAN_RADIUS} because two free lengths would make the
     *                        size effect and the thickness effect inseparable
     * @param throughWall     elements through the wall
     * @param alongAxis       elements along the axis
     */
    public record Setup(double failureScale, double fractureEnergy, double correlation,
                        double length, double slenderness, int throughWall, int alongAxis) {

        public static Setup nominal() {
            return new Setup(SCALE, FRACTURE_ENERGY, CORRELATION, TUBE_LENGTH,
                    SLENDERNESS, THROUGH_WALL, ALONG_AXIS);
        }

        /** The same tube with no defects in it. */
        public Setup perfect() {
            return new Setup(0.0, fractureEnergy, correlation, length, slenderness,
                    throughWall, alongAxis);
        }

        public Setup withFailureScale(double s) {
            return new Setup(s, fractureEnergy, correlation, length, slenderness,
                    throughWall, alongAxis);
        }

        public Setup withFractureEnergy(double g) {
            return new Setup(failureScale, g, correlation, length, slenderness,
                    throughWall, alongAxis);
        }

        public Setup withCorrelation(double c) {
            return new Setup(failureScale, fractureEnergy, c, length, slenderness,
                    throughWall, alongAxis);
        }

        /** A different length of the same tube, meshed at the same element size. */
        public Setup withLength(double l) {
            return new Setup(failureScale, fractureEnergy, correlation, l, slenderness,
                    throughWall, (int) Math.round(alongAxis * l / length));
        }

        /**
         * A different wall thickness on the same mean radius, meshed with the same number of
         * elements through the wall.
         *
         * <p>Holding the element <em>count</em> rather than the element <em>size</em> is the
         * choice that keeps a thickness sweep comparable: the through-wall discretisation
         * error is then the same at every point, so what moves between points is the physics
         * and not the mesh. It does mean a thin wall costs more, because the radial element
         * shrinks with the wall and the CFL timestep follows it.
         */
        public Setup withSlenderness(double s) {
            return new Setup(failureScale, fractureEnergy, correlation, length, s,
                    throughWall, alongAxis);
        }

        /** The same tube on a mesh refined by an integer factor in both directions. */
        public Setup refined(int factor) {
            return new Setup(failureScale, fractureEnergy, correlation, length, slenderness,
                    throughWall * factor, alongAxis * factor);
        }

        /** The same tube refined along the axis only, which is where the band width lives. */
        public Setup refinedAxially(int factor) {
            return new Setup(failureScale, fractureEnergy, correlation, length, slenderness,
                    throughWall, alongAxis * factor);
        }

        /** Axial element size, metres, which is the band width for a hoop-driven burst. */
        public double bandWidth() {
            return length / alongAxis;
        }

        /** Wall thickness, metres. */
        public double thickness() {
            return 2.0 * MEAN_RADIUS / slenderness;
        }
    }

    /**
     * Called at every sample of a run, for anything that wants to watch rather than wait.
     *
     * <p>A shot reports one pressure and one place, and a run that produced them took a few
     * hundred thousand steps to do it. This is the hook that lets something -- a film, a
     * trace, a live plot -- see the states in between without the physics having to know what
     * is looking at it.
     */
    public interface Observer {

        /** @param capacity the wall's current load capacity, Pa, as the run measures it */
        void sample(ExplicitSolver solver, QuadMesh mesh, double capacity);
    }

    public static Shot fire(long seed) {
        return fire(seed, Setup.nominal());
    }

    public static Shot fire(long seed, Setup setup) {
        return fire(seed, setup, null);
    }

    /**
     * @param seed     picks the realisation of the defect field, and nothing else
     * @param observer watches every sample, or null
     */
    public static Shot fire(long seed, Setup setup, Observer observer) {
        final double scale = setup.failureScale();
        final double fractureEnergy = setup.fractureEnergy();
        final int nr = setup.throughWall();
        final int nz = setup.alongAxis();
        final double length = setup.length();

        final Material copper = BurstCase.COPPER;
        final JohnsonCook law = copper.johnsonCook();

        final double thickness = setup.thickness();
        final double inner = MEAN_RADIUS - 0.5 * thickness;
        final double outer = MEAN_RADIUS + 0.5 * thickness;

        final QuadMesh mesh = QuadMesh.cylinderWall(inner, outer, length, nr, nz);
        final ExplicitSolver solver = new ExplicitSolver(mesh, copper, Formulation.AXISYMMETRIC,
                Integration.REDUCED, Kinematics.FINITE_STRAIN, 1.0, CFL_SAFETY);
        solver.fixAllAxial();

        double weakest = 0.0;
        double mean = 0.0;
        if (scale > 0.0) {
            final DefectField field = new DefectField(new Weibull(scale, MODULUS),
                    REFERENCE_VOLUME, setup.correlation(), seed);
            final double[] failure = field.sampleOnto(mesh);
            weakest = Double.MAX_VALUE;
            for (double f : failure) {
                if (f < weakest) weakest = f;
                mean += f;
            }
            mean /= failure.length;
            solver.setDamage(fractureEnergy, failure);
        }

        final double exact = Burst.burstPressure(law, thickness, MEAN_RADIUS);
        final double ceiling = BurstCase.OVERSHOOT * exact;

        final double omega = copper.barWaveSpeed() / MEAN_RADIUS;
        final double period = 2.0 * Math.PI / omega;
        solver.setPressureRamp(ceiling, BurstCase.RAMP_PERIODS * period);
        solver.setRelaxationDamping(BurstCase.DAMPING, omega);

        boolean traversed = false;
        double peakCapacity = 0.0, peakStrain = 0.0, peakApplied = 0.0;
        double failureZ = Double.NaN, maxDamage = 0.0, damaged = 0.0, bulge = 1.0;
        double fracture = 0.0, plastic = 0.0;
        int clamped = 0;

        final double deadline = BurstCase.MAX_PERIODS * period;
        final long t0 = System.nanoTime();
        while (solver.time() < deadline) {
            solver.run(BurstCase.SAMPLE_STEPS);
            final double capacity = BurstCase.loadCapacity(solver, mesh);
            if (observer != null) observer.sample(solver, mesh, capacity);

            if (capacity > peakCapacity) {
                peakCapacity = capacity;
                peakStrain = BurstCase.hoopStrain(solver, mesh, MEAN_RADIUS);
                peakApplied = solver.borePressure();
            } else if (solver.borePressure() >= ceiling
                    && capacity < peakCapacity * (1.0 - BurstCase.PEAK_MARGIN)) {
                traversed = true;
                break;
            }
        }

        // Past the peak, and kept running. Everything above answers "at what pressure"; the
        // peak is where the tube started to give way, and a tube that has only just started is
        // still very nearly uniform. "Where" is a question about what happens next, and it
        // needs the run to continue into a regime with no equilibrium in it at all -- which is
        // legitimate, because past the peak the tube really is accelerating, and the only
        // reason BurstCase stops at the margin is that it is measuring a pressure rather than
        // watching a failure.
        final double runOn = solver.time() + POST_PEAK_PERIODS * period;
        while (traversed && solver.time() < runOn && solver.maxDamage() < 1.0) {
            solver.run(BurstCase.SAMPLE_STEPS);
            // Watched as well, and this is the half worth watching: everything before the
            // peak is a tube swelling evenly, and everything that makes it a failure rather
            // than a bulge happens here.
            if (observer != null) {
                observer.sample(solver, mesh, BurstCase.loadCapacity(solver, mesh));
            }
        }
        final double wallClock = (System.nanoTime() - t0) * 1e-9;

        // Read after the traverse rather than at the peak: the question "where did it fail" is
        // only answerable once it has, and the peak is where it started to, not where it ended.
        if (solver.damage() != null) {
            maxDamage = solver.maxDamage();
            damaged = solver.damagedFraction();
            clamped = solver.damage().clampedPoints();
            final int worst = solver.mostDamagedElement();
            if (worst >= 0) failureZ = mesh.centroidZ(worst);
            fracture = solver.fractureDissipation();
        }
        plastic = solver.plasticDissipation();
        bulge = bulgeRatio(solver, nr, nz);

        return new Shot(seed, traversed, peakCapacity, peakStrain, peakApplied,
                weakest, mean, failureZ, maxDamage, damaged, clamped, fracture, plastic,
                bulge, (int) solver.steps(), wallClock);
    }

    /**
     * Largest bore expansion along the tube over the mean of it. One for a tube that bulged
     * uniformly; larger once the deformation has picked a station and gone there.
     *
     * <p>This is the measurement that separates "failed" from "failed <em>somewhere</em>". A
     * burst pressure on its own cannot tell the two apart -- a uniformly weakened tube and a
     * tube with one bad ring both give way early -- and localisation is the thing the defect
     * field was built to produce.
     */
    private static double bulgeRatio(ExplicitSolver solver, int nr, int nz) {
        final double[] ur = solver.radialDisplacement();
        final int nodesR = nr + 1;
        double max = 0.0, sum = 0.0;
        for (int j = 0; j <= nz; j++) {
            final double u = ur[j * nodesR];        // the bore node at this station
            if (u > max) max = u;
            sum += u;
        }
        final double average = sum / (nz + 1);
        return average > 0.0 ? max / average : 1.0;
    }
}
