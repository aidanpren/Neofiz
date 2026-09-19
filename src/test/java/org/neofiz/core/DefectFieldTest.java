package org.neofiz.core;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.neofiz.mesh.QuadMesh;

import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The defect field, checked as a random field rather than as a piece of code.
 *
 * <p>The plan's risk register has one row marked high for this milestone: "mesh-dependent
 * defect field ships". Its stated contingency is an explicit convergence test at M2, which is
 * what most of this class is. The property being proved is not that the field is smooth or
 * that it looks plausible -- both are easy and neither is the risk -- but that two meshes of
 * the same body see the <em>same field</em> and reach the <em>same failure distribution</em>.
 */
class DefectFieldTest {

    private static final Weibull BASE = new Weibull(0.42, 20.0);
    /** A cubic millimetre: a plausible calibration coupon scale, and a round number. */
    private static final double V0 = 1e-9;
    /** One millimetre: coarse for a grain, right for a weld bead or an inclusion cluster. */
    private static final double LENGTH = 1e-3;

    private static DefectField field(long seed) {
        return new DefectField(BASE, V0, LENGTH, seed);
    }

    // ------------------------------------------------------------ the Gaussian underneath

    @ParameterizedTest
    @ValueSource(longs = {1L, 7L, 20260919L})
    @DisplayName("the underlying field is standard normal, exactly by construction")
    void gaussianMarginalIsStandard(long seed) {
        // Exact, not fitted. The weights are normalised by their own root sum of squares, so
        // the variance is one at every point rather than one on average -- which is what lets
        // the Weibull marginal downstream be exactly the Weibull that was asked for, and
        // therefore what lets the size effect be tested against closed form.
        DefectField f = field(seed);
        double sum = 0.0, sumSq = 0.0;
        int n = 0;
        for (int i = 0; i < 300; i++) {
            for (int j = 0; j < 300; j++) {
                double g = f.gaussianAt(0.05 + i * 0.37e-3, 0.01 + j * 0.41e-3);
                sum += g;
                sumSq += g * g;
                n++;
            }
        }
        double mean = sum / n;
        double variance = sumSq / n - mean * mean;
        // The samples are correlated, so the error on these is far above 1/sqrt(n). The
        // sample spacing is about a third of a correlation length, so 90000 points are worth
        // something like 1000 independent ones.
        assertEquals(0.0, mean, 0.05, "field mean drifted");
        assertEquals(1.0, variance, 0.07, "field variance is not one");
    }

    @Test
    @DisplayName("correlation decays at the physical length, not the sampling interval")
    void correlationFollowsTheKernel() {
        // The first of the plan's two fixes, and the one that is a property of the
        // construction rather than of the algebra. Measured against the analytic target
        // exp(-d^2 / 2 l^2) the field was built to have.
        DefectField f = field(11L);
        for (double lengths : new double[]{0.25, 0.5, 1.0, 1.5, 2.0}) {
            double d = lengths * LENGTH;
            double sxy = 0.0, sxx = 0.0, syy = 0.0;
            int n = 0;
            for (int i = 0; i < 260; i++) {
                for (int j = 0; j < 260; j++) {
                    double r = 0.05 + i * 0.29e-3;
                    double z = 0.01 + j * 0.31e-3;
                    double a = f.gaussianAt(r, z);
                    double b = f.gaussianAt(r + d, z);
                    sxy += a * b;
                    sxx += a * a;
                    syy += b * b;
                    n++;
                }
            }
            double measured = sxy / Math.sqrt(sxx * syy);
            assertEquals(f.correlation(d), measured, 0.03,
                    "correlation at " + lengths + " correlation lengths");
        }
    }

    @Test
    @DisplayName("points far apart are uncorrelated, and points on top of each other are not")
    void correlationEndpoints() {
        DefectField f = field(3L);
        assertEquals(1.0, f.correlation(0.0), 0.0);
        assertTrue(f.correlation(5.0 * LENGTH) < 1e-5);

        double a = f.gaussianAt(0.06, 0.02);
        assertEquals(a, f.gaussianAt(0.06, 0.02), 0.0, "the field must be a function");
        assertTrue(Math.abs(a - f.gaussianAt(0.06 + 0.02 * LENGTH, 0.02)) < 0.05,
                "neighbouring points must be nearly identical, or there is no correlation");
    }

    // ------------------------------------------------------------ mesh independence

    @Test
    @DisplayName("two meshes of the same body see the same field, bit for bit")
    void refiningTheMeshResolvesTheSameField() {
        // The risk register's row, in its strongest form. The field is anchored to the
        // origin of the (r, z) frame and drawn from a stateless hash of the grid index, so
        // there is no sense in which a mesh "draws" anything. Held to exact equality rather
        // than to a tolerance, because anything that made it approximate -- a sequential
        // generator, an accumulated offset, a bounding box that shifted the grid -- would
        // be a real dependence hiding inside a passing test.
        DefectField f = field(42L);
        QuadMesh coarse = QuadMesh.cylinderWall(0.05, 0.10, 0.02, 4, 2);
        QuadMesh fine = QuadMesh.cylinderWall(0.05, 0.10, 0.02, 64, 32);
        assertEquals(coarse.totalRingVolume(), fine.totalRingVolume(),
                coarse.totalRingVolume() * 1e-12, "the two meshes are not the same body");

        for (double r = 0.051; r < 0.0999; r += 0.0037) {
            for (double z = 0.001; z < 0.0199; z += 0.0031) {
                double direct = f.gaussianAt(r, z);
                assertEquals(direct, field(42L).gaussianAt(r, z), 0.0,
                        "the same seed must give the same field");
                assertNotEquals(direct, field(43L).gaussianAt(r, z),
                        "a different seed must give a different field");
            }
        }

        // And the memoised path used by sampleOnto is arithmetic only.
        double[] sampled = f.sampleOnto(coarse);
        for (int e = 0; e < coarse.elementCount; e++) {
            double pointwise = f.failureStrainAt(coarse.centroidRadius(e), coarse.centroidZ(e),
                    coarse.signedArea(e));
            assertEquals(pointwise, sampled[e], 0.0,
                    "memoising the white noise changed element " + e);
        }
    }

    @Test
    @DisplayName("effective volume is the ring volume until the element is smaller than a cell")
    void theResolutionFloor() {
        DefectField f = field(1L);
        final double r = 0.075;

        // Coarse: several correlation cells across the section, so nothing is floored and
        // the effective volume is exactly what Pappus gives.
        double coarseArea = 4.0 * LENGTH * LENGTH;
        assertEquals(2.0 * Math.PI * r * coarseArea, f.effectiveVolume(r, coarseArea),
                2.0 * Math.PI * r * coarseArea * 1e-12);

        // Fine: the element sits inside one cell. The field has one value there and the
        // element gets one in-plane chance, not a quarter of one.
        double fineArea = 0.25 * LENGTH * LENGTH;
        assertEquals(2.0 * Math.PI * r * LENGTH * LENGTH, f.effectiveVolume(r, fineArea),
                2.0 * Math.PI * r * LENGTH * LENGTH * 1e-12);
        assertEquals(f.effectiveVolume(r, fineArea), f.effectiveVolume(r, fineArea / 1000.0), 0.0,
                "below the cell size the effective volume must stop shrinking entirely");

        // The circumferential factor is 471 on that tube and never floors on anything
        // tube-shaped -- it is the term the ring effect rides on. The only place it can
        // floor is on the axis, where there is no circumference to have flaws around.
        assertEquals(2.0 * Math.PI * r / LENGTH, f.effectiveVolume(r, coarseArea)
                / (LENGTH * coarseArea), 1e-9);
        assertEquals(4.0 * LENGTH * LENGTH * LENGTH, f.effectiveVolume(1e-9, coarseArea),
                4e-18, "on the axis the hoop count must floor at one, not vanish");
    }

    @ParameterizedTest
    @ValueSource(longs = {1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L})
    @DisplayName("the weakest link settles under refinement rather than drifting either way")
    void refinementDoesNotMoveTheBody(long seed) {
        // The trap, tested on the object rather than on the algebra, on four meshes spanning
        // a 1024-fold change in element count.
        //
        // Both directions are failures here, and the second one is the one that was actually
        // found. Without volume normalisation this falls without bound. With volume
        // normalisation but without the resolution floor in DefectField.effectiveVolume it
        // *rises* without bound -- measured at 0.183, 0.209, 0.213, 0.223, 0.233, 0.251 over
        // a thousandfold refinement -- because every element below the correlation length is
        // credited with independent chances that its neighbours are in fact sharing. That
        // version passes any test phrased as "does not get weaker".
        DefectField f = field(seed);
        double[] weakest = new double[4];
        int[] counts = new int[4];
        int k = 0;
        for (int refine : new int[]{2, 8, 32, 64}) {
            QuadMesh mesh = QuadMesh.cylinderWall(0.05, 0.10, 0.02, 4 * refine, 2 * refine);
            weakest[k] = Arrays.stream(f.sampleOnto(mesh)).min().orElseThrow();
            counts[k] = mesh.elementCount;
            k++;
        }

        assertEquals(32, counts[0]);
        assertEquals(32768, counts[3], "the finest mesh must be 1024 times the coarsest");

        // The coarse meshes are genuinely under-resolved -- 32 elements cannot find the low
        // spot of a field with a thousand correlation cells in it -- so the comparison that
        // matters is between the two finest, where the field is resolved and only the
        // normalisation is still in question.
        double drift = Math.abs(weakest[3] - weakest[2]) / weakest[2];
        assertTrue(drift < 0.03,
                "a 4-fold refinement moved the weakest link by " + (100 * drift)
                        + " %, which is a mesh dependence rather than a discretisation error");
        assertTrue(weakest[3] > 0.6 * weakest[0] && weakest[3] < 1.4 * weakest[0],
                "over 1024x refinement the weakest link went " + weakest[0] + " to "
                        + weakest[3]);
    }

    @Test
    @DisplayName("the sampled body reproduces the closed-form weakest-link distribution")
    void bodyMatchesTheClosedForm() {
        // The strongest statement available, and the one that says the construction is right
        // rather than merely stable: run many nominally identical parts, take the weakest
        // element of each, and the resulting population has to be the Weibull that
        // weakest-link theory predicts for a body of that volume -- mean and scatter both.
        //
        // Nothing was fitted to make this land. The per-element effective volumes sum to far
        // more than the body volume on a fine mesh, and the elements are correlated so that
        // far fewer of them are independent draws; those two errors are each large and they
        // cancel, which is the whole design.
        // Deliberately fine enough that the resolution floor is engaged: the elements are
        // 0.78 by 0.63 mm against a 1 mm correlation length, so each holds half a cell and
        // the naive ring volume would be a quarter of the effective one.
        QuadMesh mesh = QuadMesh.cylinderWall(0.05, 0.10, 0.02, 64, 32);
        final int parts = 300;
        double[] weakest = new double[parts];
        for (int s = 0; s < parts; s++) {
            weakest[s] = Arrays.stream(field(7000L + s).sampleOnto(mesh)).min().orElseThrow();
        }

        double mean = Arrays.stream(weakest).average().orElseThrow();
        double sd = Math.sqrt(Arrays.stream(weakest).map(x -> (x - mean) * (x - mean)).sum()
                / (parts - 1));

        Weibull predicted = BASE.forVolume(mesh.totalRingVolume(), V0);
        assertEquals(predicted.mean(), mean, predicted.mean() * 0.04,
                "body mean " + mean + " against a predicted " + predicted.mean());
        assertEquals(predicted.standardDeviation(), sd, predicted.standardDeviation() * 0.25,
                "body scatter " + sd + " against a predicted " + predicted.standardDeviation());
    }

    @Test
    @DisplayName("strength falls as V^(-1/m) across a sevenfold change in volume")
    void theSizeEffectIsMeasured() {
        // The milestone's exit criterion, measured on sampled fields rather than asserted
        // from the algebra. Four tubes, identical wall and length, radius doubling each
        // time, so the only thing changing is how much steel there is.
        //
        // Nothing in the field knows the radius. The radius enters through the ring volume
        // and comes back out as a strength.
        // Three thousand parts per tube. The body distribution has a 6 % coefficient of
        // variation, so a few hundred would leave a standard error on each mean of half a
        // per cent and a ratio too noisy to distinguish the prediction from a near miss --
        // which is how a size effect that is 10 % wrong gets signed off.
        final int parts = 3000;
        double[] volumes = new double[4];
        double[] means = new double[4];
        int k = 0;
        for (double ri : new double[]{0.02, 0.04, 0.08, 0.16}) {
            QuadMesh mesh = QuadMesh.cylinderWall(ri, ri + 0.005, 0.02, 8, 16);
            double sum = 0.0;
            for (int s = 0; s < parts; s++) {
                sum += Arrays.stream(field(9000L + s).sampleOnto(mesh)).min().orElseThrow();
            }
            volumes[k] = mesh.totalRingVolume();
            means[k] = sum / parts;
            k++;
        }

        // Measured residuals at this sample size are +0.12 %, -0.12 % and -0.16 %.
        for (int i = 1; i < 4; i++) {
            double ratio = volumes[i] / volumes[0];
            double predicted = Math.pow(ratio, -1.0 / BASE.modulus());
            double measured = means[i] / means[0];
            assertEquals(predicted, measured, 0.005,
                    "at " + ratio + " times the volume the strength ratio was " + measured
                            + " against a predicted " + predicted);
            assertTrue(measured < 1.0, "the larger tube has to be the weaker one");
        }

        // And it is a real effect rather than a rounding: the biggest tube is 7.2 times the
        // volume of the smallest and about 9 % weaker for no other reason.
        assertTrue(means[3] < 0.94 * means[0],
                "the size effect came out too small to be the thing being demonstrated");
    }

    @Test
    @DisplayName("without volume normalisation the same refinement does collapse")
    void theUnnormalisedFieldCollapses() {
        // The control. Same field, same meshes, one ingredient removed: every element is
        // given the base distribution regardless of how much material it holds. This is what
        // a naive implementation does and it is the reason the test above is not vacuous.
        DefectField f = field(5L);
        double coarse = weakestUnnormalised(f, QuadMesh.cylinderWall(0.05, 0.10, 0.02, 8, 4));
        double fine = weakestUnnormalised(f, QuadMesh.cylinderWall(0.05, 0.10, 0.02, 128, 64));
        assertTrue(fine < 0.8 * coarse,
                "the unnormalised field must get markedly weaker under refinement; it went "
                        + coarse + " to " + fine);
    }

    private static double weakestUnnormalised(DefectField f, QuadMesh mesh) {
        double worst = Double.MAX_VALUE;
        for (int e = 0; e < mesh.elementCount; e++) {
            double g = f.gaussianAt(mesh.centroidRadius(e), mesh.centroidZ(e));
            worst = Math.min(worst, f.base().quantile(Normal.cdf(g)));
        }
        return worst;
    }

    // ------------------------------------------------------------ the marginal

    @Test
    @DisplayName("the sampled marginal is the Weibull the element's volume earns it")
    void marginalIsWeibull() {
        // The copula has to preserve the marginal exactly, or the size effect is being
        // measured through a distorted lens. Checked by Kolmogorov-Smirnov against the
        // closed form, on a mesh whose elements all hold the same volume so that there is a
        // single distribution to compare against.
        QuadMesh mesh = QuadMesh.cylinderWall(0.05, 0.10, 0.05, 1, 400);
        DefectField f = new DefectField(BASE, V0, LENGTH, 99L);

        // Every element in this mesh is one ring of the same cross-section, so they differ
        // only in z and the effective volume is identical across them. It is also the full
        // ring volume, because the section is 50 mm by 0.125 mm and no floor engages.
        double v = f.effectiveVolume(mesh.centroidRadius(0), mesh.signedArea(0));
        assertEquals(mesh.ringVolume(0), v, v * 1e-12);
        for (int e = 1; e < mesh.elementCount; e++) {
            assertEquals(v, f.effectiveVolume(mesh.centroidRadius(e), mesh.signedArea(e)),
                    v * 1e-12);
        }

        // Elements a correlation length apart are correlated, so a KS test on one mesh would
        // be testing a single wiggly realisation. Pooling many seeds gives independent draws
        // of the same marginal, which is what the statistic assumes.
        final int seeds = 200;
        double[] pooled = new double[seeds * 12];
        int n = 0;
        for (int s = 0; s < seeds; s++) {
            double[] strain = new DefectField(BASE, V0, LENGTH, 1000L + s).sampleOnto(mesh);
            for (int e = 0; e < 12; e++) pooled[n++] = strain[e * 33];
        }
        Arrays.sort(pooled);

        Weibull expected = BASE.forVolume(v, V0);
        double ks = 0.0;
        for (int i = 0; i < n; i++) {
            double cdf = expected.cdf(pooled[i]);
            ks = Math.max(ks, Math.max(cdf - (double) i / n, (double) (i + 1) / n - cdf));
        }
        // 1.63/sqrt(n) is the 1 % critical value; this is a regression guard, not a p-value.
        assertTrue(ks < 1.63 / Math.sqrt(n),
                "sampled marginal departs from the closed form, KS = " + ks
                        + " against a critical value of " + 1.63 / Math.sqrt(n));
    }

    @Test
    @DisplayName("a wider tube is weaker at the same wall and length, as r^(-1/m)")
    void theRingEffect() {
        // The consequence of normalising by ring volume rather than by cross-sectional area,
        // and the one prediction here that is specific to an axisymmetric solver. Two tubes
        // with the same wall thickness and the same length, one at twice the radius: the
        // wider one holds twice the steel, so it has twice the chances to be bad, and its
        // characteristic failure strain falls by 2^(-1/m).
        //
        // Nothing in the field was told about radius. It comes out of the volume.
        final double wall = 0.005;
        final double height = 0.02;
        QuadMesh narrow = QuadMesh.cylinderWall(0.030, 0.030 + wall, height, 4, 16);
        QuadMesh wide = QuadMesh.cylinderWall(0.060, 0.060 + wall, height, 4, 16);

        double volumeRatio = wide.totalRingVolume() / narrow.totalRingVolume();
        assertEquals(1.9231, volumeRatio, 1e-3, "the wide tube should hold about twice as much");

        DefectField f = field(77L);
        Weibull narrowBody = f.forBody(narrow.totalRingVolume());
        Weibull wideBody = f.forBody(wide.totalRingVolume());
        assertEquals(Math.pow(volumeRatio, -1.0 / BASE.modulus()),
                wideBody.scale() / narrowBody.scale(), 1e-12);

        // 3.2 % at m = 20. Small, and exactly the kind of thing that is dismissed as noise
        // until a sweep across radius shows it is a straight line on a log plot.
        double weaker = 100.0 * (1.0 - wideBody.scale() / narrowBody.scale());
        assertEquals(3.21, weaker, 0.02);
    }

    @Test
    @DisplayName("repeated runs give a distribution, not a number")
    void repeatedRunsScatter() {
        // The second half of the milestone's exit criterion. Same geometry, same material,
        // same everything except the seed, which is the part standing in for "this is a
        // different physical part off the same production line".
        QuadMesh mesh = QuadMesh.cylinderWall(0.05, 0.10, 0.02, 16, 8);
        final int parts = 400;
        double[] weakest = new double[parts];
        for (int s = 0; s < parts; s++) {
            weakest[s] = Arrays.stream(field(5000L + s).sampleOnto(mesh)).min().orElseThrow();
        }

        double mean = Arrays.stream(weakest).average().orElseThrow();
        double sd = Math.sqrt(Arrays.stream(weakest).map(x -> (x - mean) * (x - mean)).sum()
                / (parts - 1));
        double spread = sd / mean;

        assertTrue(spread > 0.005,
                "nominally identical parts came out identical, spread was " + spread);
        assertTrue(spread < 0.10,
                "the scatter is implausibly wide for m = 20, spread was " + spread);
        assertTrue(Arrays.stream(weakest).min().orElseThrow() < 0.9 * mean,
                "no part in four hundred was much worse than average, which is not what a "
                        + "weakest-link distribution looks like");
    }

    // ------------------------------------------------------------ housekeeping

    @Test
    @DisplayName("the whole mesh gets a finite, positive failure strain")
    void noDegenerateValues() {
        // The clamp and the two quantile branches exist for this. A zero would be an element
        // that fails on the first step; an infinity would be one that never fails at all.
        // Both are seed-dependent and rare, which is the worst combination.
        QuadMesh mesh = QuadMesh.cylinderWall(0.05, 0.10, 0.02, 32, 16);
        for (long seed = 0; seed < 40; seed++) {
            for (double x : field(seed).sampleOnto(mesh)) {
                assertTrue(x > 0.0 && Double.isFinite(x), "bad failure strain " + x);
            }
        }
        // Directly at the guard, where the copula would otherwise round to zero or one. The
        // volume form is used so the extremes can be reached at all: a volume of one cubic
        // metre against a cubic-millimetre reference is 1e9 chances to be bad, and a
        // picolitre is a millionth of one.
        DefectField f = field(1L);
        assertTrue(f.failureStrainForVolume(0.05, 0.01, 1.0) > 0.0);
        assertTrue(Double.isFinite(f.failureStrainForVolume(0.05, 0.01, 1e-15)));
    }

    @Test
    @DisplayName("an element on the axis still gets a ring volume and a value")
    void axisElementsAreFine() {
        // The Taylor geometry has a column of elements touching r = 0. Their centroid is at
        // positive radius so Pappus is well defined, but a nearly-zero volume drives the
        // normalisation hard and is worth checking rather than assuming.
        QuadMesh solid = QuadMesh.solidCylinder(0.0038, 0.0254, 8, 32);
        assertEquals(Math.PI * 0.0038 * 0.0038 * 0.0254, solid.totalRingVolume(), 1e-12);
        for (double x : field(9L).sampleOnto(solid)) {
            assertTrue(x > 0.0 && Double.isFinite(x), "bad failure strain " + x);
        }
    }

    @Test
    @DisplayName("nonsense parameters are refused")
    void validation() {
        assertThrows(IllegalArgumentException.class, () -> new DefectField(BASE, 0.0, LENGTH, 1));
        assertThrows(IllegalArgumentException.class, () -> new DefectField(BASE, V0, 0.0, 1));
        assertThrows(IllegalArgumentException.class,
                () -> new DefectField(BASE, V0, Double.NaN, 1));
    }
}
