package org.neofiz.validate;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.neofiz.validate.DefectBurst.Population;
import org.neofiz.validate.DefectBurst.Setup;
import org.neofiz.validate.DefectBurst.Shot;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The end of M2: a tube that does not have one answer.
 *
 * <p>Everything up to here has been asserted as an identity, because everything up to here had
 * one. A population does not. So these assert the <em>shape</em> of the distribution -- that it
 * has spread, that the spread moves the way a stated mechanism says it should, and that the
 * failure location is not quietly pinned to a corner of the mesh -- and they are written so
 * that the thing being claimed is the mechanism rather than the number.
 *
 * <p>Short tubes and coarse meshes throughout, because these are checking that a behaviour
 * exists and has a sign, not measuring it. The measured figures are in the report.
 */
class DefectBurstTest {

    /** Small enough to run a whole population in a test suite. */
    private static Setup small() {
        return Setup.nominal().withLength(20.0e-3);
    }

    @Test
    @DisplayName("the same tube, fired twice, bursts at two pressures in two places")
    void theSameTubeBurstsDifferently() {
        // The plan's actual claim, and the reason any of M2 exists.
        Setup setup = small();
        double[] pressure = new double[5];
        double[] where = new double[5];
        for (int s = 1; s <= 5; s++) {
            Shot shot = DefectBurst.fire(s, setup);
            assertTrue(shot.traversed(), "seed " + s + " never reached a load maximum");
            pressure[s - 1] = shot.peakCapacity();
            where[s - 1] = shot.failureZ();
            assertEquals(0, shot.clampedPoints(),
                    "seed " + s + " hit the snap-back clamp; the mesh is too coarse for G_f");
            assertTrue(shot.maxDamage() > 0.5,
                    "seed " + s + " never really failed, max damage " + shot.maxDamage());
        }

        for (int i = 0; i < pressure.length; i++) {
            for (int j = i + 1; j < pressure.length; j++) {
                assertNotEquals(pressure[i], pressure[j],
                        "seeds " + (i + 1) + " and " + (j + 1) + " burst identically");
            }
        }

        long places = java.util.Arrays.stream(where).distinct().count();
        assertTrue(places >= 3, "five tubes must fail in at least three places, got " + places);
    }

    @Test
    @DisplayName("the same seed is the same tube, bit for bit")
    void theSameSeedIsTheSameTube() {
        // Without this the scatter above is not a defect field, it is a race condition.
        Shot a = DefectBurst.fire(7, small());
        Shot b = DefectBurst.fire(7, small());
        assertEquals(a.peakCapacity(), b.peakCapacity(), 0.0);
        assertEquals(a.failureZ(), b.failureZ(), 0.0);
        assertEquals(a.steps(), b.steps());
        assertEquals(a.weakestStrain(), b.weakestStrain(), 0.0);
    }

    @Test
    @DisplayName("a tube with no defects in it still bursts where M1 said it would")
    void aPerfectTubeReproducesTheClosedForm() {
        // The control. If this drifts, the scatter below is being measured against a moving
        // reference and none of it means anything.
        Shot perfect = DefectBurst.fire(1, small().perfect());
        assertTrue(perfect.traversed(), "the control must reach a load maximum");
        assertEquals(0.0, perfect.maxDamage(), 0.0, "and must not damage anything");
        assertEquals(1.0, perfect.bulgeRatio(), 1e-6,
                "a defect-free tube must bulge perfectly uniformly");

        double ratio = perfect.peakCapacity() / DefectBurst.referencePressure();
        assertEquals(1.0, ratio, 0.01,
                "the defect-free tube must land on the closed form, got " + ratio);
    }

    @Test
    @DisplayName("a weaker material bursts at a lower pressure, monotonically")
    void theWeakerTheMaterialTheLowerTheBurst() {
        double previous = Double.MAX_VALUE;
        for (double scale : new double[] {0.08, 0.05, 0.03}) {
            double p = DefectBurst.fire(1, small().withFailureScale(scale)).peakCapacity();
            assertTrue(p < previous,
                    "a failure strain of " + scale + " must not burst higher than the one above");
            previous = p;
        }
        assertTrue(previous < 0.95 * DefectBurst.referencePressure(),
                "the weakest material must knock the burst pressure down by more than 5 %");
    }

    @Test
    @DisplayName("the tube averages its defects instead of failing at its weakest link")
    void theTubeAveragesInsteadOfFailingAtItsWeakestLink() {
        // Weakest-link theory predicts a longer body is weaker and scatters by the same
        // fraction, because the minimum of n Weibulls is Weibull with the same modulus. A tube
        // does neither: the mean holds still and the scatter falls like an average, because a
        // shell shares load along its axis over sqrt(R t) and every weak patch here is shorter
        // than that. This is the finding, so it is asserted rather than described.
        Population shortTube = DefectBurst.volley(6, small().withLength(10.0e-3));
        Population longTube = DefectBurst.volley(6, small().withLength(40.0e-3));

        assertEquals(0, shortTube.notTraversed() + longTube.notTraversed(),
                "every tube in both populations must reach a load maximum");

        assertEquals(shortTube.meanRatio(), longTube.meanRatio(), 0.02,
                "a fourfold longer tube must not be meaningfully weaker: "
                        + shortTube.meanRatio() + " then " + longTube.meanRatio());

        assertTrue(longTube.scatterPercent() < 0.8 * shortTube.scatterPercent(),
                "and its scatter must fall like an average, not hold like a weakest link: "
                        + shortTube.scatterPercent() + " % then " + longTube.scatterPercent()
                        + " %");
    }

    @Test
    @DisplayName("the correlation length decides which statistics apply")
    void theCorrelationLengthDecidesWhichStatisticsApply() {
        // The mechanism behind the test above, exercised directly. Below sqrt(R t) the
        // neighbours carry a weak patch; above it they cannot, and weakest-link behaviour
        // reappears -- weaker and more variable, both.
        double lag = DefectBurst.shearLagLength();
        assertEquals(7.07e-3, lag, 0.01e-3, "sqrt(R t) for this tube");

        Setup base = small().withLength(40.0e-3);
        Population carried = DefectBurst.volley(6, base.withCorrelation(0.3 * lag));
        Population alone = DefectBurst.volley(6, base.withCorrelation(3.0 * lag));

        assertTrue(alone.meanRatio() < carried.meanRatio() - 0.01,
                "patches longer than the shear-lag length must weaken the tube: "
                        + carried.meanRatio() + " then " + alone.meanRatio());
        assertTrue(alone.scatterPercent() > 1.5 * carried.scatterPercent(),
                "and must scatter it: " + carried.scatterPercent() + " % then "
                        + alone.scatterPercent() + " %");
    }

    @Test
    @DisplayName("the failure location is spread along the tube, not pinned by the mesh")
    void theFailureLocationIsSpreadAlongTheTube() {
        // A location that always lands at an end, or always in the middle, would mean the
        // boundary conditions were choosing and the field was decorative.
        Population p = DefectBurst.volley(8, small().withLength(40.0e-3));
        assertTrue(p.locationSpreadFraction() > 0.15,
                "failure locations are bunched: spread is "
                        + p.locationSpreadFraction() + " of the tube, uniform would be 0.289");
    }

    @Test
    @DisplayName("the burst pressure does NOT converge under axial refinement")
    void theBurstPressureDoesNotConverge() {
        // The honest one. Crack-band regularisation fixes the energy per unit crack area --
        // DamageTest proves that to 1e-10 -- but it fixes it by making the softening modulus
        // depend on the element size, and that is only the right thing to do once the band has
        // localised into one element. Here the load maximum arrives while damage is still
        // diffuse, so a finer mesh softens more slowly at the same strain and the tube comes
        // out stronger, monotonically, with no plateau. Asserted so that a future fix has to
        // change this test rather than quietly pass it.
        Setup base = small().withLength(40.0e-3);
        double[] peak = new double[3];
        for (int i = 0; i < 3; i++) {
            peak[i] = DefectBurst.fire(1, base.refinedAxially(1 << i)).peakCapacity();
        }

        assertTrue(peak[0] < peak[1] && peak[1] < peak[2],
                "refinement must move the answer one way: " + java.util.Arrays.toString(peak));

        double drift = (peak[2] - peak[0]) / peak[0];
        assertTrue(drift > 0.001 && drift < 0.02,
                "and it drifts by about half a per cent over a fourfold refinement, not by "
                        + (100.0 * drift) + " %");

        // Heading for the defect-free answer, which is where it ends up if the mesh is refined
        // without bound: the knockdown is what evaporates.
        double perfect = DefectBurst.fire(1, base.perfect()).peakCapacity();
        assertTrue(peak[2] < 0.98 * perfect,
                "but it must still be a long way from having evaporated");
    }

    // ------------------------------------------------------------------ wall thickness

    @Test
    @DisplayName("a thicker wall holds more, in proportion, with defects in it")
    void slendernessScales() {
        // The design map's x axis. Burst pressure goes as t/r, so halving the slenderness at
        // a fixed mean radius should double the capacity -- and it should keep doing that
        // with a defect field switched on, because the knockdown is a fraction rather than a
        // pressure. If it did not, the map's columns would not be comparable.
        Setup thin = small().withSlenderness(40.0);
        Setup thick = small().withSlenderness(20.0);

        double thinPeak = DefectBurst.fire(1, thin).peakCapacity();
        double thickPeak = DefectBurst.fire(1, thick).peakCapacity();

        assertEquals(2.0, thickPeak / thinPeak, 0.05,
                "halving D/t should about double the burst pressure; got "
                        + (thickPeak / thinPeak));

        // Both knocked down from their own defect-free reference by a similar fraction, which
        // is what makes a derating rule a rule rather than a table.
        double thinRatio = thinPeak / DefectBurst.referencePressure(40.0);
        double thickRatio = thickPeak / DefectBurst.referencePressure(20.0);
        assertTrue(thinRatio < 0.98 && thickRatio < 0.98,
                "the defect field should cost both tubes something: "
                        + thinRatio + " and " + thickRatio);
        assertEquals(thinRatio, thickRatio, 0.02,
                "the knockdown should be a fraction, not a pressure: "
                        + thinRatio + " against " + thickRatio);
    }

    @Test
    @DisplayName("the geometry helpers agree with the setup they describe")
    void geometryHelpersAgree() {
        Setup setup = small().withSlenderness(20.0);
        assertEquals(2.0 * DefectBurst.MEAN_RADIUS / 20.0, setup.thickness(), 0.0);
        // To a rounding error rather than to the bit: the two expressions multiply the same
        // three numbers in a different order.
        assertEquals(Math.sqrt(DefectBurst.MEAN_RADIUS * setup.thickness()),
                DefectBurst.shearLagLength(20.0), 1e-15);

        // The no-argument forms are the nominal tube and must not drift away from it.
        assertEquals(DefectBurst.referencePressure(DefectBurst.SLENDERNESS),
                DefectBurst.referencePressure(), 0.0);
        assertEquals(DefectBurst.shearLagLength(DefectBurst.SLENDERNESS),
                DefectBurst.shearLagLength(), 0.0);
        assertEquals(DefectBurst.SLENDERNESS, Setup.nominal().slenderness(), 0.0);
    }
}
