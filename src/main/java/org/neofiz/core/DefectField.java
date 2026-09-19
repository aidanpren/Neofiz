package org.neofiz.core;

import org.neofiz.mesh.QuadMesh;

/**
 * A spatially correlated random field of local failure strain, sampled onto a mesh.
 *
 * <p>One extra number per element. It is the cheapest object in the codebase and it is the one
 * that buys the most, because several things that would otherwise each need a rule become
 * consequences of it: failure happens somewhere plausible rather than exactly where the
 * continuum stress peaks, bigger parts are weaker, nominally identical parts burst at
 * different pressures, and a weld fails first because its elements draw from a worse
 * distribution rather than because a rule said welds fail.
 *
 * <h2>The trap this class exists to avoid</h2>
 *
 * A per-element independent draw from a fixed distribution is pathologically mesh-dependent.
 * Refine the mesh, take more samples, find a worse extreme, and the part gets weaker without
 * bound -- a simulation that changes its answer when the player changes a setting they did not
 * think was physical. Two fixes, and both are load-bearing:
 *
 * <ol>
 *   <li><b>A correlation length set by a physical scale</b> -- grain size, inclusion spacing,
 *       weld bead width -- not by the element size. The field is generated on a grid of that
 *       scale, anchored to the origin of the (r, z) frame, and sampled onto elements wherever
 *       they happen to be. Refining the mesh then <em>resolves the same field better</em>
 *       instead of drawing a new one. Two meshes of the same body see the same field at the
 *       same physical point, bit for bit.</li>
 *   <li><b>Volume normalisation.</b> An element holding volume V is given the distribution of
 *       the weakest of {@code V/V0} reference volumes, which is again Weibull -- see
 *       {@link Weibull#minimumOf}. Halving the element size halves its hazard and doubles the
 *       element count, and the body's failure distribution does not move. That is exact
 *       algebra, not a calibration.</li>
 * </ol>
 *
 * <h2>Why those two fixes are not enough on their own</h2>
 *
 * Applied naively together they trade one mesh dependence for another, in the opposite
 * direction, and the second one is the more flattering of the two so it is the more dangerous.
 *
 * <p>Volume normalisation says an element of volume V had {@code V/V0} independent chances to
 * be bad. That is true only while those chances really are independent. Refine past the
 * correlation length and neighbouring elements read almost the same value of the field -- they
 * are one chance shared out, not many -- while each is separately credited with the full
 * strengthening of its own small volume. The body gets stronger as the mesh is refined, by
 * about {@code k^(1/m)} per {@code k}-fold refinement, and it does not stop.
 *
 * <p>Measured before the fix, on a cylinder at m = 20: the weakest element went 0.183, 0.209,
 * 0.213, 0.223, 0.233, 0.251 over a thousandfold refinement. Every one of those numbers looks
 * reasonable. The sequence is a 37 % drift with nothing physical behind it.
 *
 * <p>So the count of independent chances is <b>floored at what the field can actually
 * resolve</b>. See {@link #effectiveVolume}. Above the correlation length nothing changes and
 * the effective volume is exactly the ring volume; below it, the element stops being credited
 * with in-plane chances it does not have.
 *
 * <h2>The volume is a ring volume, and that is not a detail</h2>
 *
 * An axisymmetric element is not a small chunk of material with a flaw in it. It is a full
 * 360-degree torus, {@code 2 pi r A}, and its stored failure strain has to stand for the
 * weakest of every flaw around that circumference. An element at r = 60 mm with a 1 mm square
 * section holds 377 cubic millimetres of steel; at a 50 micrometre correlation length that is
 * some 7500 independent chances to be bad, not one.
 *
 * <p>Normalising by {@link QuadMesh#ringVolume} <em>is</em> that correction. There is no
 * separate circumferential term to add, because the weakest of n samples and the distribution
 * of a volume n times larger are the same Weibull. The visible consequence is that a wide tube
 * is weaker than a narrow one of the same wall and length -- strength falls as
 * {@code r^(-1/m)} -- which is real, and which nothing in the solver was told to do.
 *
 * <p>The assumption being made is that the circumferential correlation length equals the
 * in-plane one. That is the only self-consistent choice for an isotropic defect population and
 * it is stated here rather than buried, because a drawn or rolled part has a texture and this
 * field does not know about it.
 *
 * <h2>Construction</h2>
 *
 * Standard Gaussian copula. White noise on the physical grid, smoothed by a separable Gaussian
 * kernel into a field with correlation {@code exp(-d^2 / 2 l^2)}, normalised to <em>exactly</em>
 * unit variance at every point by dividing by the root sum of squared weights, then pushed
 * through {@link Normal#cdf} to a uniform and through {@link Weibull#quantile} to the marginal.
 * Because the normalisation is exact, the marginal is exactly the Weibull asked for, which is
 * what lets the size effect be tested against closed form rather than against itself.
 *
 * <p>The white noise is drawn from a <b>stateless hash</b> of the seed and the integer grid
 * coordinates, not from a sequential generator. A sequential generator would make the value at
 * a grid node depend on how many nodes had been visited first, so the field would change when
 * the mesh bounding box changed, when the traversal order changed, or when the work was
 * threaded -- reintroducing the mesh dependence this class exists to remove, in a form that
 * would be far harder to see.
 */
public final class DefectField {

    /** Grid nodes per correlation length. Sets how well the smoothing kernel is resolved. */
    private static final int CELLS_PER_LENGTH = 3;
    /** Kernel truncation, in kernel standard deviations. Beyond 4 the weight is 3e-4. */
    private static final double KERNEL_RADII = 4.0;
    /**
     * Where the Gaussian field is clamped, in standard deviations. A guard for totality, not
     * physics: past this the copula's uniform rounds to 0 or 1 and the Weibull quantile
     * returns zero or infinity -- an element that fails instantly or never fails at all. At
     * nine sigma the probability excluded is 1e-19, so nothing measurable is changed.
     */
    private static final double CLAMP_SIGMA = 9.0;
    /** Largest white-noise grid worth memoising. Beyond it, hash on demand. */
    private static final int MAX_CACHE_NODES = 1 << 24;

    private final Weibull base;
    private final double referenceVolume;
    private final double correlationLength;
    private final long seed;

    private final double cell;
    private final double sigmaCells;
    private final int radius;

    /**
     * @param base              the failure-strain distribution of one reference volume
     * @param referenceVolume   the volume {@code base} was calibrated on, cubic metres
     * @param correlationLength the physical scale of the defect population, metres
     * @param seed              selects one realisation; the same seed is the same part
     */
    public DefectField(Weibull base, double referenceVolume, double correlationLength, long seed) {
        if (!(referenceVolume > 0.0)) {
            throw new IllegalArgumentException(
                    "reference volume must be positive, was " + referenceVolume);
        }
        if (!(correlationLength > 0.0)) {
            throw new IllegalArgumentException(
                    "correlation length must be positive, was " + correlationLength);
        }
        this.base = base;
        this.referenceVolume = referenceVolume;
        this.correlationLength = correlationLength;
        this.seed = seed;

        this.cell = correlationLength / CELLS_PER_LENGTH;
        // A Gaussian kernel of standard deviation s has an autocorrelation of standard
        // deviation s*sqrt(2), so s = l/sqrt(2) gives the field the correlation asked for.
        this.sigmaCells = (correlationLength / Math.sqrt(2.0)) / cell;
        this.radius = (int) Math.ceil(KERNEL_RADII * sigmaCells);
    }

    public Weibull base() {
        return base;
    }

    public double referenceVolume() {
        return referenceVolume;
    }

    public double correlationLength() {
        return correlationLength;
    }

    public long seed() {
        return seed;
    }

    /** The correlation this field is built to have between two points {@code d} apart. */
    public double correlation(double d) {
        double q = d / correlationLength;
        return Math.exp(-0.5 * q * q);
    }

    /**
     * The underlying unit-variance Gaussian field at a point. Exposed because it is the thing
     * whose statistics can be checked directly -- the marginal and the correlation both live
     * here, before the copula, where they are a standard normal and a known function rather
     * than something that has to be inferred from Weibull samples.
     */
    public double gaussianAt(double r, double z) {
        return smooth(r, z, this::white);
    }

    /**
     * The volume of material an element's single stored value actually stands for: the ring
     * it sweeps, but never resolved finer than the correlation length in any direction.
     *
     * <pre>  V_eff = l^3 * max(1, A / l^2) * max(1, 2 pi r / l)</pre>
     *
     * <p>Read the two factors separately, because they are floored for different reasons.
     *
     * <p><b>In plane</b>, the field is a real function of position and an element bigger than a
     * correlation cell covers {@code A/l^2} of them while only being read once -- so it is
     * credited with that many chances, which is the standard local-averaging correction. An
     * element <em>smaller</em> than a cell covers one, and the reading at its centroid is that
     * one. Crediting it with a fraction of a chance would make it spuriously strong, and that
     * is the drift documented above.
     *
     * <p><b>Around the circumference</b>, the field is not a function of position at all --
     * there is no third coordinate -- so the {@code 2 pi r / l} cells there are never resolved
     * at any mesh density and the element is always credited with all of them. This is where
     * the ring effect comes from, and why it does not wash out under refinement.
     *
     * <p>The floors are a {@code max} rather than a smooth blend. A blend would be defensible
     * -- Vanmarcke's variance function is the usual one -- but it would put a fitted shape in
     * the middle of a chain that is otherwise exact, to smooth a corner each element crosses
     * once.
     *
     * @param centroidRadius   radius of the element's area centroid, metres
     * @param crossSectionArea its area in the (r, z) half-plane, square metres
     */
    public double effectiveVolume(double centroidRadius, double crossSectionArea) {
        double inPlane = Math.max(1.0, crossSectionArea / (correlationLength * correlationLength));
        double hoop = Math.max(1.0, 2.0 * Math.PI * centroidRadius / correlationLength);
        return correlationLength * correlationLength * correlationLength * inPlane * hoop;
    }

    /**
     * Local failure strain for an axisymmetric element, sampled at its centroid.
     *
     * <p>Two arguments that look independent and are not: the position picks which realisation
     * of the field is being read, and the size decides how many independent chances that much
     * material had to be worse than the reading suggests.
     */
    public double failureStrainAt(double r, double z, double crossSectionArea) {
        return transform(gaussianAt(r, z), effectiveVolume(r, crossSectionArea));
    }

    /**
     * Local failure strain for a caller that has worked out the effective volume itself --
     * the planar formulations, where the out-of-plane direction is a nominal thickness rather
     * than a circumference, and the tests that need to vary the volume directly.
     */
    public double failureStrainForVolume(double r, double z, double volume) {
        return transform(gaussianAt(r, z), volume);
    }

    /**
     * Samples the field onto every element of a mesh, returning failure strain per element.
     *
     * <p>Each element is evaluated at its own area centroid and normalised by its own
     * effective volume, so nothing accumulates and the result does not depend on the order
     * elements are visited. Memoises the white noise over the mesh's bounding box, which is
     * arithmetic only: the memoised value is the same hash the point-wise path computes, so
     * the two agree bit for bit.
     */
    public double[] sampleOnto(QuadMesh mesh) {
        Noise noise = cacheFor(mesh);
        double[] out = new double[mesh.elementCount];
        for (int e = 0; e < mesh.elementCount; e++) {
            double rc = mesh.centroidRadius(e);
            double g = smooth(rc, mesh.centroidZ(e), noise);
            out[e] = transform(g, effectiveVolume(rc, mesh.signedArea(e)));
        }
        return out;
    }

    /**
     * The failure-strain distribution a body of this volume draws from -- the weakest link
     * over the whole thing. The closed form the sampled field has to reproduce.
     */
    public Weibull forBody(double volume) {
        return base.forVolume(volume, referenceVolume);
    }

    // ------------------------------------------------------------- the transform

    /** Gaussian to Weibull, through the copula, with the volume's worth of extra chances. */
    private double transform(double g, double volume) {
        Weibull local = base.forVolume(volume, referenceVolume);
        double clamped = Math.max(-CLAMP_SIGMA, Math.min(CLAMP_SIGMA, g));
        // Each branch keeps the digits in the tail it is responsible for. Going through the
        // cdf on the strong side would round the uniform to exactly one and hand back an
        // infinite failure strain; going through the survival on the weak side would lose the
        // weak tail, which is the half that decides where the part breaks.
        return clamped <= 0.0
                ? local.quantile(Normal.cdf(clamped))
                : local.quantileFromSurvival(Normal.survival(clamped));
    }

    // ------------------------------------------------------------- the field

    /** A source of standard normal white noise on the integer grid. */
    private interface Noise {
        double at(int i, int j);
    }

    /**
     * Weighted sum of the white noise around a point, normalised to exactly unit variance.
     *
     * <p>The kernel is separable, so the weights cost {@code 2(2R+2)} exponentials rather than
     * {@code (2R+2)^2}, and the normaliser factorises exactly. Dividing by the root sum of
     * squared weights -- rather than by a constant -- is what makes the variance exactly one
     * at every point instead of only away from the grid nodes, and an exact variance is what
     * makes the marginal exactly Weibull.
     */
    private double smooth(double r, double z, Noise noise) {
        final double gr = r / cell;
        final double gz = z / cell;
        final int i0 = (int) Math.floor(gr) - radius;
        final int j0 = (int) Math.floor(gz) - radius;
        final int taps = 2 * radius + 2;
        final double twoSigmaSq = 2.0 * sigmaCells * sigmaCells;

        double[] wr = new double[taps];
        double[] wz = new double[taps];
        double sumRSq = 0.0;
        double sumZSq = 0.0;
        for (int k = 0; k < taps; k++) {
            double dr = gr - (i0 + k);
            double dz = gz - (j0 + k);
            wr[k] = Math.exp(-dr * dr / twoSigmaSq);
            wz[k] = Math.exp(-dz * dz / twoSigmaSq);
            sumRSq += wr[k] * wr[k];
            sumZSq += wz[k] * wz[k];
        }

        double sum = 0.0;
        for (int b = 0; b < taps; b++) {
            double row = 0.0;
            for (int a = 0; a < taps; a++) {
                row += wr[a] * noise.at(i0 + a, j0 + b);
            }
            sum += wz[b] * row;
        }
        return sum / Math.sqrt(sumRSq * sumZSq);
    }

    /**
     * Standard normal white noise at a grid node, as a pure function of the seed and the node.
     * Box-Muller over two SplitMix64 hashes.
     */
    private double white(int i, int j) {
        long h = mix(seed + 0x9E3779B97F4A7C15L);
        h = mix(h ^ (i * 0xBF58476D1CE4E5B9L));
        h = mix(h ^ (j * 0x94D049BB133111EBL));
        double u1 = unitOpen(h);
        double u2 = unitOpen(mix(h ^ 0xD6E8FEB86659FD93L));
        return Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
    }

    private static long mix(long z) {
        z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
        z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
        return z ^ (z >>> 31);
    }

    /**
     * 53 bits to the open interval (0, 1). Open at zero because Box-Muller takes a logarithm
     * of it, and a closed interval produces an infinity once every few billion nodes -- which
     * is rarely enough to survive every test and to appear in someone's long run.
     */
    private static double unitOpen(long bits) {
        return ((bits >>> 11) + 0.5) * 0x1.0p-53;
    }

    // ------------------------------------------------------------- memoisation

    private Noise cacheFor(QuadMesh mesh) {
        double rMin = Double.MAX_VALUE, rMax = -Double.MAX_VALUE;
        double zMin = Double.MAX_VALUE, zMax = -Double.MAX_VALUE;
        for (int n = 0; n < mesh.nodeCount; n++) {
            rMin = Math.min(rMin, mesh.r[n]);
            rMax = Math.max(rMax, mesh.r[n]);
            zMin = Math.min(zMin, mesh.z[n]);
            zMax = Math.max(zMax, mesh.z[n]);
        }
        int i0 = (int) Math.floor(rMin / cell) - radius;
        int j0 = (int) Math.floor(zMin / cell) - radius;
        int ni = (int) Math.floor(rMax / cell) + radius + 2 - i0;
        int nj = (int) Math.floor(zMax / cell) + radius + 2 - j0;

        long nodes = (long) ni * nj;
        if (nodes <= 0 || nodes > MAX_CACHE_NODES) return this::white;

        double[] v = new double[(int) nodes];
        for (int b = 0; b < nj; b++) {
            for (int a = 0; a < ni; a++) {
                v[b * ni + a] = white(i0 + a, j0 + b);
            }
        }
        final int fi0 = i0, fj0 = j0, fni = ni, fnj = nj;
        return (i, j) -> {
            int a = i - fi0;
            int b = j - fj0;
            // Outside the box only when a mesh node sits exactly on the far boundary and
            // rounding pushes a window one cell past it; falling back keeps the value
            // identical to the uncached path rather than clamping to an edge.
            if (a < 0 || a >= fni || b < 0 || b >= fnj) return white(i, j);
            return v[b * fni + a];
        };
    }
}
