package org.neofiz.core;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

import java.util.SplittableRandom;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The weakest-link algebra, tested away from any field or mesh.
 *
 * <p>Almost everything here is an identity rather than an approximation, so almost everything
 * is held to round-off. That matters: the size effect is the one part of M2 whose answer is
 * known in advance, and if it is checked to three digits then a construction that is quietly
 * wrong by two per cent passes and the mesh-dependence test downstream inherits the error.
 */
class WeibullTest {

    private static final Weibull STEEL = new Weibull(0.42, 30.0);
    private static final Weibull IRON = new Weibull(0.012, 8.0);

    @Test
    @DisplayName("the scale is the value at which survival is 1/e")
    void scaleIsTheCharacteristicValue() {
        // The definition, and the reason "scale" is not a fitting constant with no meaning.
        assertEquals(Math.exp(-1.0), STEEL.survival(STEEL.scale()), 1e-15);
        assertEquals(Math.exp(-1.0), IRON.survival(IRON.scale()), 1e-15);
        assertEquals(1.0, STEEL.cumulativeHazard(STEEL.scale()), 1e-15);
    }

    @Test
    @DisplayName("survival and cdf are complementary, including in the weak tail")
    void survivalAndCdfAgree() {
        for (double x = 0.01; x < 1.0; x *= 1.3) {
            assertEquals(1.0, STEEL.survival(x) + STEEL.cdf(x), 1e-15);
        }
        // Where expm1 earns its place: at a hazard of 1e-14 the naive 1 - exp(-h) has lost
        // two digits, and this is exactly the regime a volume-normalised element sits in.
        Weibull tiny = new Weibull(1.0, 1.0);
        assertEquals(1e-14, tiny.cdf(1e-14), 1e-28);
    }

    @Test
    @DisplayName("quantile inverts the cdf over the whole range, both tails included")
    void quantileInvertsCdf() {
        for (double p = 1e-12; p < 1.0; p *= 3.0) {
            assertEquals(p, STEEL.cdf(STEEL.quantile(p)), Math.abs(p) * 1e-12);
        }
        // Approaching one from below is where quantile(p) would die if it used log(1 - p):
        // at p = 1 - 1e-15 the subtraction leaves a single significant bit.
        for (double s = 1e-15; s < 0.5; s *= 3.0) {
            assertEquals(s, STEEL.survival(STEEL.quantileFromSurvival(s)), Math.abs(s) * 1e-12);
        }

        // The two forms agree only where the caller could have written either, which is
        // where 1 - s still has digits left. Below about 1e-9 it does not, and the pair
        // parts company by more than the tolerance -- which is the reason both exist and is
        // therefore asserted rather than avoided.
        for (double s = 1e-9; s < 0.5; s *= 3.0) {
            assertEquals(STEEL.quantileFromSurvival(s), STEEL.quantile(1.0 - s),
                    STEEL.scale() * 1e-9,
                    "the two quantile forms must agree where both are still accurate");
        }
        assertTrue(Math.abs(STEEL.quantileFromSurvival(1e-15) - STEEL.quantile(1.0 - 1e-15))
                        > STEEL.scale() * 1e-8,
                "if these now agree in the far tail, 1 - s has stopped losing precision and "
                        + "quantileFromSurvival has lost its reason to exist");
    }

    @Test
    @DisplayName("the weakest of n samples is Weibull with the scale cut by n^(-1/m)")
    void minimumOfNIsWeibull() {
        // The identity the whole milestone rests on, checked in the two forms it gets used
        // in. The hazard form is the exact one -- n independent chances is n times the
        // hazard, with no exponential in the way -- so it is held to round-off.
        for (Weibull w : new Weibull[]{STEEL, IRON}) {
            for (double n : new double[]{1.0, 2.0, 7.5, 1000.0, 1e6}) {
                Weibull weakest = w.minimumOf(n);
                assertEquals(w.modulus(), weakest.modulus(), 0.0,
                        "the modulus must survive; only the scale moves");
                for (double x : new double[]{0.2 * w.scale(), w.scale(), 1.5 * w.scale()}) {
                    double expected = n * w.cumulativeHazard(x);
                    assertEquals(expected, weakest.cumulativeHazard(x), expected * 1e-13,
                            "minimum of " + n + " at x = " + x);
                }
            }
        }

        // And the probability form, which is the same statement read out loud: the chance
        // all n survive is the chance one survives, to the n. Looser on purpose -- raising a
        // probability to the millionth power multiplies its relative error by a million, and
        // that is arithmetic in the test rather than error in the distribution.
        for (double n : new double[]{1.0, 2.0, 7.5, 1000.0, 1e6}) {
            double x = 0.9 * STEEL.scale();
            double expected = Math.pow(STEEL.survival(x), n);
            assertEquals(expected, STEEL.minimumOf(n).survival(x), expected * 1e-8,
                    "survival of the weakest of " + n);
        }
    }

    @Test
    @DisplayName("one sample is the distribution itself")
    void minimumOfOneIsIdentity() {
        assertEquals(STEEL.scale(), STEEL.minimumOf(1.0).scale(), 0.0);
        assertEquals(STEEL.modulus(), STEEL.minimumOf(1.0).modulus(), 0.0);
    }

    @Test
    @DisplayName("taking the weakest twice is taking the weakest of the product")
    void minimumComposes() {
        // Needed because the field applies volume normalisation once per element while the
        // body-level closed form applies it once over the whole body. If these did not
        // compose, the two would be answering different questions.
        Weibull twice = STEEL.minimumOf(12.0).minimumOf(5.0);
        Weibull once = STEEL.minimumOf(60.0);
        assertEquals(once.scale(), twice.scale(), once.scale() * 1e-14);
        assertEquals(once.modulus(), twice.modulus(), 0.0);
    }

    @Test
    @DisplayName("volume normalisation is the weakest-of-n calculation, not a second one")
    void forVolumeIsMinimumOf() {
        // Exact equality, because forVolume is implemented as a call to minimumOf. If someone
        // later writes the formula out a second time for speed, this fails on the last bit
        // and says so, which is the point.
        assertEquals(STEEL.minimumOf(8.0).scale(), STEEL.forVolume(8e-6, 1e-6).scale(), 0.0);
    }

    @Test
    @DisplayName("eight times the volume at m = 10 is 18.8 % weaker, which the plan rounds to 19")
    void theHeadlineNumber() {
        // The one number from M2 a reader is most likely to check by hand, so it is pinned
        // rather than left implied. The plan says "about 19 %"; it is 18.77, and the figure
        // written down here is the computed one.
        Weibull w = new Weibull(1.0, 10.0);
        assertEquals(Math.pow(8.0, -0.1), w.sizeFactor(8.0, 1.0), 1e-15);
        assertEquals(18.775, 100.0 * (1.0 - w.sizeFactor(8.0, 1.0)), 5e-3);

        // And the reason the modulus is the parameter that decides a material's character:
        // the same eight-fold size change is worth under 7 % in tight ductile steel and
        // nearly 23 % in a brittle material, from one number and nothing else.
        assertEquals(6.697, 100.0 * (1.0 - STEEL.sizeFactor(8.0, 1.0)), 5e-3);
        assertEquals(22.889, 100.0 * (1.0 - IRON.sizeFactor(8.0, 1.0)), 5e-3);
    }

    @Test
    @DisplayName("the body's failure distribution does not depend on how it is partitioned")
    void hazardIsAdditiveOverAnyPartition() {
        // Mesh independence, in its purest form and before any mesh exists. Cut a body into
        // uneven pieces, give each the distribution its own volume earns it, and the hazards
        // must sum to the hazard of the whole. This is what makes refining a mesh not change
        // the answer: it is a different partition of the same volume, and nothing else.
        final double v0 = 1e-9;
        final double body = 3.7e-4;
        SplittableRandom rng = new SplittableRandom(20260919L);

        for (Weibull w : new Weibull[]{STEEL, IRON}) {
            for (int parts : new int[]{2, 17, 500, 20000}) {
                double[] cut = new double[parts];
                double total = 0.0;
                for (int i = 0; i < parts; i++) {
                    cut[i] = rng.nextDouble(0.1, 1.0);
                    total += cut[i];
                }

                final double x = 0.8 * w.scale();
                double summed = 0.0;
                for (int i = 0; i < parts; i++) {
                    summed += w.forVolume(body * cut[i] / total, v0).cumulativeHazard(x);
                }
                double whole = w.forVolume(body, v0).cumulativeHazard(x);

                assertEquals(whole, summed, whole * 1e-11,
                        "hazard changed when the body was cut into " + parts + " pieces");
            }
        }
    }

    @Test
    @DisplayName("refining a partition without volume normalisation destroys the body")
    void theTrapItself() {
        // The failure mode this milestone exists to avoid, written down so that the fix has
        // something to be a fix for. Draw one value per element from the unnormalised base
        // law and take the weakest: the body's characteristic strength falls as the mesh is
        // refined, without bound and without anything in the setup having changed.
        final double v0 = 1e-9;
        final double body = 3.7e-4;

        double naiveCoarse = STEEL.minimumOf(100).scale();
        double naiveFine = STEEL.minimumOf(100000).scale();
        assertTrue(naiveFine < 0.8 * naiveCoarse,
                "an unnormalised per-element draw has to get weaker under refinement, "
                        + "otherwise this test is not describing the trap");

        // With normalisation, the same two meshes agree exactly, because both are the same
        // body volume however it was cut.
        double fixedCoarse = STEEL.forVolume(body / 100, v0).minimumOf(100).scale();
        double fixedFine = STEEL.forVolume(body / 100000, v0).minimumOf(100000).scale();
        assertEquals(fixedCoarse, fixedFine, fixedCoarse * 1e-14);
        assertEquals(STEEL.forVolume(body, v0).scale(), fixedCoarse, fixedCoarse * 1e-14);
    }

    @ParameterizedTest
    @ValueSource(doubles = {1.0, 2.0, 5.0, 8.0, 30.0, 50.0})
    @DisplayName("the mean matches a numerical integral of the survival function")
    void meanMatchesIntegratedSurvival(double m) {
        // An independent check on the gamma function, which otherwise has nothing to be
        // wrong against. E[X] = integral of the survival function over the positive axis.
        Weibull w = new Weibull(0.42, m);
        final int n = 2_000_000;
        // Far enough out that the survival function has fallen to 1e-12, which for a small
        // modulus is a long way: an m = 1 Weibull is an exponential and still has 3e-4 of
        // its probability left at eight times the scale.
        final double top = w.scale() * Math.max(8.0, Math.pow(27.6, 1.0 / m));
        final double h = top / n;
        double area = 0.5 * (w.survival(0.0) + w.survival(top));
        for (int i = 1; i < n; i++) area += w.survival(i * h);
        assertEquals(w.mean(), area * h, w.mean() * 2e-6, "mean disagrees at m = " + m);
    }

    @Test
    @DisplayName("the mean and scatter match the two closed forms that are known exactly")
    void knownSpecialCases() {
        // m = 1 is the exponential distribution and m = 2 is the Rayleigh. Both have gamma
        // values that are exact numbers rather than series evaluations, so they pin the
        // Lanczos approximation against something that cannot itself be slightly wrong.
        Weibull exponential = new Weibull(3.0, 1.0);
        assertEquals(3.0, exponential.mean(), 1e-14);
        assertEquals(3.0, exponential.standardDeviation(), 1e-13);
        assertEquals(1.0, exponential.coefficientOfVariation(), 1e-13);

        Weibull rayleigh = new Weibull(3.0, 2.0);
        assertEquals(3.0 * Math.sqrt(Math.PI) / 2.0, rayleigh.mean(), 1e-14);
        assertEquals(3.0 * Math.sqrt(1.0 - Math.PI / 4.0), rayleigh.standardDeviation(), 1e-14);
    }

    @Test
    @DisplayName("the modulus alone sets the scatter, and sets it to the quoted figures")
    void scatterDependsOnModulusOnly() {
        assertEquals(new Weibull(1.0, 8.0).coefficientOfVariation(),
                new Weibull(1e6, 8.0).coefficientOfVariation(), 1e-15,
                "scatter must be scale-free, or the modulus is not the character parameter");

        // What a modulus means in plain terms, which is the translation the design map needs:
        // m = 8 is 15 % scatter and m = 30 is 4 %.
        assertEquals(0.148369, IRON.coefficientOfVariation(), 1e-5);
        assertEquals(0.041769, STEEL.coefficientOfVariation(), 1e-5);
    }

    @Test
    @DisplayName("gamma is right where its value is known in closed form")
    void gammaSpotChecks() {
        assertEquals(1.0, Weibull.gamma(1.0), 1e-14);
        assertEquals(1.0, Weibull.gamma(2.0), 1e-14);
        assertEquals(6.0, Weibull.gamma(4.0), 1e-13);
        assertEquals(362880.0, Weibull.gamma(10.0), 1e-8);
        assertEquals(Math.sqrt(Math.PI), Weibull.gamma(0.5), 1e-14);
        assertEquals(Math.sqrt(Math.PI) / 2.0, Weibull.gamma(1.5), 1e-14);
        assertEquals(-2.0 * Math.sqrt(Math.PI), Weibull.gamma(-0.5), 1e-13,
                "gamma is negative between -1 and 0; the reflection has to carry the sign");
    }

    @Test
    @DisplayName("nonsense parameters are refused, NaN included")
    void validation() {
        assertThrows(IllegalArgumentException.class, () -> new Weibull(0.0, 10.0));
        assertThrows(IllegalArgumentException.class, () -> new Weibull(-1.0, 10.0));
        assertThrows(IllegalArgumentException.class, () -> new Weibull(1.0, 0.0));
        // A NaN modulus passes the obvious "m <= 0" guard and then turns every sampled
        // failure strain into a NaN, which reads downstream as an element that never fails.
        assertThrows(IllegalArgumentException.class, () -> new Weibull(1.0, Double.NaN));
        assertThrows(IllegalArgumentException.class, () -> new Weibull(Double.NaN, 10.0));
        assertThrows(IllegalArgumentException.class, () -> STEEL.minimumOf(0.0));
        assertThrows(IllegalArgumentException.class, () -> STEEL.forVolume(1.0, 0.0));
        assertThrows(IllegalArgumentException.class, () -> STEEL.quantile(1.5));
    }
}
