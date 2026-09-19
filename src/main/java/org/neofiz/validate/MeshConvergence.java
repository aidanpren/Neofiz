package org.neofiz.validate;

/**
 * Richardson extrapolation and the Grid Convergence Index over three systematically refined
 * meshes.
 *
 * <h2>What this is for</h2>
 *
 * Every quantity this project reports from a mesh is wrong by a discretisation error, and the
 * honest question is never "is the mesh fine enough" but "how big is that error and which way
 * does it point". Three meshes at a fixed refinement ratio answer both, provided the
 * differences are shrinking:
 *
 * <pre>
 *   p     = ln|d1 / d2| / ln r          observed order of convergence
 *   f_ext = f_fine + d2 / (r^p - 1)     the zero-mesh-size value
 *   GCI   = Fs |d2 / f_fine| / (r^p - 1)   an uncertainty band on f_fine
 * </pre>
 *
 * <p>with {@code d1 = medium - coarse} and {@code d2 = fine - medium}. This is the standard
 * verification procedure (Roache; ASME V&amp;V 20), and the safety factor {@code Fs = 1.25} is
 * the usual one for a three-grid study where the order is observed rather than assumed.
 *
 * <p><b>The observed order is the interesting output, not the extrapolated value.</b> A
 * quantity that converges at the rate the element promises is merely under-resolved and will
 * come good with refinement. One that converges more slowly is telling you something about
 * the problem: a bilinear quad should give second order on a smooth field, so a measured
 * order near one says the quantity is being read next to a singularity, and no amount of
 * refinement will restore the missing order. That is a reason to extrapolate and quote a
 * band, not a reason to buy a bigger machine.
 *
 * <h2>What it refuses to do</h2>
 *
 * The arithmetic above will produce a number from any three values at all, including
 * sequences that are oscillating or diverging, and that number looks exactly as authoritative
 * as a real one. So the guards are not conveniences:
 *
 * <ul>
 *   <li>If the two differences have opposite signs the sequence is <b>oscillatory</b>, there
 *       is no order to observe, and extrapolating is meaningless.</li>
 *   <li>If the second difference is not smaller than the first the sequence is <b>not
 *       converging</b>; {@code p} comes out zero or negative and the extrapolation runs the
 *       wrong way, confidently.</li>
 * </ul>
 *
 * <p>Both throw rather than returning a plausible-looking value.
 *
 * @param coarse          value on the coarsest mesh
 * @param medium          value on the intermediate mesh
 * @param fine            value on the finest mesh
 * @param refinementRatio mesh size ratio between successive levels, 2 for halving
 */
public record MeshConvergence(double coarse, double medium, double fine,
                              double refinementRatio) {

    /** Roache's safety factor for a three-grid study with an observed order. */
    public static final double SAFETY = 1.25;

    public MeshConvergence {
        if (refinementRatio <= 1.0) {
            throw new IllegalArgumentException("refinement ratio must exceed 1");
        }
    }

    /** Three meshes each half the element size of the last. */
    public static MeshConvergence halving(double coarse, double medium, double fine) {
        return new MeshConvergence(coarse, medium, fine, 2.0);
    }

    private double d1() {
        return medium - coarse;
    }

    private double d2() {
        return fine - medium;
    }

    /** Whether both steps move the same way, so there is a trend rather than a wobble. */
    public boolean isMonotone() {
        return d1() * d2() > 0.0;
    }

    /** Whether the steps are shrinking, which is what makes an order meaningful. */
    public boolean isConverging() {
        return isMonotone() && Math.abs(d1()) > Math.abs(d2());
    }

    /**
     * Whether the sequence is in a usable asymptotic range: the correction Richardson wants to
     * apply must be no larger than the movement refinement has actually produced.
     *
     * <p>This is the guard that shrinking steps alone do not give. As the two differences
     * approach each other the denominator {@code r^p - 1} goes to zero and the extrapolation
     * runs away, so a quantity that has <em>already converged</em> -- whose remaining
     * differences are run-to-run residual rather than a mesh trend -- produces a large, stable
     * looking, entirely fictional correction. The Taylor final length does exactly this: it
     * moves one part in ten thousand per refinement, {@code |d1/d2|} comes out at 1.04, and
     * the extrapolation confidently relocates it twelve times further than the whole refinement
     * study moved it.
     *
     * <p>A converged quantity is not a failure and should not be dressed as one. The right
     * report for it is the finest value and the size of the last change, which is what
     * {@link #summary} falls back to.
     */
    public boolean isAsymptotic() {
        if (!isConverging()) return false;
        return Math.abs(d2() / (Math.abs(d1() / d2()) - 1.0)) <= Math.abs(fine - coarse);
    }

    /**
     * The ratio {@code r^p}, which is where every formula here actually gets its work done.
     *
     * <p>Computed directly as {@code |d1/d2|} rather than by raising the refinement ratio to
     * the observed order: the order was derived from this ratio in the first place, so going
     * back through a logarithm and an exponential only adds round-trip error to a quantity
     * that is already exact.
     */
    private double ratio() {
        if (!isMonotone()) {
            throw new IllegalStateException(String.format(
                    "the sequence %.6g, %.6g, %.6g is oscillatory, not converging: the steps "
                            + "are %+.3g then %+.3g. There is no order to observe and "
                            + "extrapolating would invent one",
                    coarse, medium, fine, d1(), d2()));
        }
        if (!isConverging()) {
            throw new IllegalStateException(String.format(
                    "the sequence %.6g, %.6g, %.6g is not converging: the step grew from "
                            + "%.3g to %.3g under refinement. Extrapolation would run the "
                            + "wrong way",
                    coarse, medium, fine, Math.abs(d1()), Math.abs(d2())));
        }
        if (!isAsymptotic()) {
            throw new IllegalStateException(String.format(
                    "the sequence %.6g, %.6g, %.6g is not in an asymptotic range: its steps "
                            + "are %.3g then %.3g, barely shrinking, so extrapolation would "
                            + "move the answer by %.3g against the %.3g the whole refinement "
                            + "produced. This quantity has most likely already converged and "
                            + "the remainder is run-to-run residual; report the finest value "
                            + "and the last change instead",
                    coarse, medium, fine, d1(), d2(),
                    d2() / (Math.abs(d1() / d2()) - 1.0), fine - coarse));
        }
        return Math.abs(d1() / d2());
    }

    /**
     * Observed order of convergence. Compare against what the element promises -- second
     * order for a bilinear quad on a smooth field -- rather than against any fixed threshold.
     */
    public double observedOrder() {
        return Math.log(ratio()) / Math.log(refinementRatio);
    }

    /** The Richardson-extrapolated value: what an infinitely fine mesh would report. */
    public double extrapolated() {
        return fine + d2() / (ratio() - 1.0);
    }

    /**
     * Grid Convergence Index on the finest mesh, as a fraction. An uncertainty band, not an
     * error bar: the claim is that the exact value lies within this of {@link #fine}, not
     * that the error is exactly this.
     */
    public double gci() {
        return SAFETY * Math.abs(d2() / fine) / (ratio() - 1.0);
    }

    /** The uncertainty band in the same units as the values themselves. */
    public double uncertainty() {
        return gci() * Math.abs(fine);
    }

    /** Relative size of the last refinement step, for a quantity that has already converged. */
    public double lastChange() {
        return Math.abs(d2() / fine);
    }

    /**
     * One line, for a report. Extrapolated value and band where that is meaningful; for a
     * quantity that has already converged, the finest value and the size of the last step,
     * which is the honest statement rather than a manufactured extrapolation.
     */
    public String summary(String units) {
        final String suffix = units.isEmpty() ? "" : " " + units;
        if (!isAsymptotic()) {
            return String.format("%.4g%s, converged (last step %.3f%%)",
                    fine, suffix, lastChange() * 100.0);
        }
        return String.format("%.4g%s +/- %.2f%% (observed order %.2f)",
                extrapolated(), suffix, gci() * 100.0, observedOrder());
    }
}
