package org.neofiz.render;

import org.neofiz.mesh.QuadMesh;
import org.neofiz.report.Runs;
import org.neofiz.solver.ExplicitSolver;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * Deformed geometry captured over a run, so that a solve can be watched instead of summarised.
 *
 * <h2>Why this is a file and not a window</h2>
 *
 * The viewer will change -- it is a browser page today and will be the game's own renderer
 * eventually -- and the capture will not. What a frame needs to contain is settled by the
 * physics: where every node is, what every element is carrying, and when. So that is what
 * this holds, in the layout the solver already stores it in, and the drawing happens
 * somewhere else entirely. Writing the viewer first and hanging the capture off it would get
 * this backwards and would have to be undone.
 *
 * <p>It is also the cheapest thing in the pipeline by a long way. A run's frames are a few
 * hundred kilobytes; the run that produced them is hundreds of thousands of timesteps.
 *
 * <h2>Coordinates</h2>
 *
 * Positions are the <b>deformed</b> ones -- {@code r + ur}, {@code z + uz} -- because that is
 * the shape the body is actually in, and it is the entire point. Metres, unscaled. The viewer
 * fits them to its canvas; nothing here knows about pixels.
 *
 * <p>The axis convention is the solver's: {@code r} across, {@code z} along. For an
 * axisymmetric body {@code r} is a radius rather than a coordinate, so a viewer that wants to
 * show a whole tube rather than half of one mirrors the drawing about {@code r = 0}. That is
 * a decision about what the picture means and it belongs to the viewer, which is why the
 * capture stores the half the solver solved and says so rather than doubling it here.
 *
 * <h2>Precision</h2>
 *
 * Six significant figures, which on a 40 mm tube is a resolution of forty nanometres -- far
 * below anything a screen can show, and roughly a quarter the size of the file that full
 * {@code double} precision would produce. Traces that need to be recomputed from are
 * {@link org.neofiz.report.Csv}'s job and are written at seventeen.
 */
public final class Film {

    /**
     * Passed as a field's bounds to have the ramp fitted to what the run actually produced.
     *
     * <p>Worth preferring over a guessed ceiling nearly always. A hand-picked range that is
     * too wide compresses the whole field into a few adjacent colours and the picture reads
     * as uniform when it is not -- which is a failure mode that looks like a physical result
     * rather than like a bug. The exception is a quantity whose absolute value means
     * something to the reader: damage runs 0 to 1 because 0 is intact and 1 is gone, and
     * rescaling that to the range one run happened to reach would make an undamaged element
     * look damaged.
     *
     * <p>Fitted over <b>every</b> frame, never per frame. A ramp that rescaled as the film
     * played would make two frames incomparable, which is most of what a film is for.
     */
    public static final double AUTO = Double.NaN;

    /** What each element is coloured by. The name travels with the data to the viewer. */
    public record Field(String name, String units, double low, double high) {

        boolean isAuto() {
            return Double.isNaN(low) || Double.isNaN(high);
        }
    }

    private record Frame(double time, double[] r, double[] z, double[][] values,
                         double[] scalars, String live) {
    }

    private final QuadMesh mesh;
    /** False for a cross-section; see {@link #planar()}. */
    private boolean axisymmetric = true;

    private final List<Field> fields = new ArrayList<>();
    private final List<String> readouts = new ArrayList<>();
    private final List<Frame> frames = new ArrayList<>();
    private final List<ElementValue> sources = new ArrayList<>();
    private final List<Readout> meters = new ArrayList<>();
    private ElementAlive alive;

    /** A per-element quantity to colour by. */
    public interface ElementValue {
        double at(ExplicitSolver solver, int element);
    }

    /** A single number about the whole body, shown beside the picture. */
    public interface Readout {
        double of(ExplicitSolver solver);
    }

    /** Whether an element is still part of the solve. See {@link #aliveIf}. */
    public interface ElementAlive {
        boolean at(ExplicitSolver solver, int element);
    }

    public Film(QuadMesh mesh) {
        this.mesh = mesh;
    }

    /**
     * Adds a field the viewer can colour elements by.
     *
     * @param low  the value drawn as the cold end of the ramp
     * @param high the value drawn as the hot end; values outside are clamped by the viewer
     */
    public Film colourBy(String name, String units, double low, double high,
                         ElementValue value) {
        fields.add(new Field(name, units, low, high));
        sources.add(value);
        return this;
    }

    /**
     * Records which elements are still there, for a scene that can erode.
     *
     * <p>Without this a deleted element keeps being drawn in the last shape it held, which on
     * a run that erodes is a picture of debris that is not there. One character per element
     * per frame -- a few tens of kilobytes on a scene of this size, against several megabytes
     * of coordinates -- and nothing at all when it is not set.
     */
    public Film aliveIf(ElementAlive test) {
        this.alive = test;
        return this;
    }

    /** Adds a number to show beside the picture as it plays. */
    public Film readout(String name, Readout value) {
        readouts.add(name);
        meters.add(value);
        return this;
    }

    /**
     * Captures one frame. Cheap enough to call every sample of an ordinary run: it is a copy
     * of two node arrays and one pass over the elements.
     */
    public void capture(ExplicitSolver solver) {
        final double[] ur = solver.radialDisplacement();
        final double[] uz = solver.axialDisplacement();
        final double[] r = new double[mesh.nodeCount];
        final double[] z = new double[mesh.nodeCount];
        for (int i = 0; i < mesh.nodeCount; i++) {
            r[i] = mesh.r[i] + ur[i];
            z[i] = mesh.z[i] + uz[i];
        }

        final double[][] values = new double[sources.size()][mesh.elementCount];
        for (int f = 0; f < sources.size(); f++) {
            for (int e = 0; e < mesh.elementCount; e++) {
                values[f][e] = sources.get(f).at(solver, e);
            }
        }

        final double[] scalars = new double[meters.size()];
        for (int m = 0; m < meters.size(); m++) scalars[m] = meters.get(m).of(solver);

        String live = null;
        if (alive != null) {
            final char[] flags = new char[mesh.elementCount];
            for (int e = 0; e < mesh.elementCount; e++) {
                flags[e] = alive.at(solver, e) ? '1' : '0';
            }
            live = new String(flags);
        }

        frames.add(new Frame(solver.time(), r, z, values, scalars, live));
    }

    public int size() {
        return frames.size();
    }

    /**
     * Writes the film as a JavaScript file that assigns one global.
     *
     * <p>A {@code .js} assignment rather than a {@code .json} fetch because the viewer is
     * served under a content policy that allows a script tag and refuses an XHR, and because
     * a file that loads by being included cannot fail halfway through and leave the page
     * showing an empty canvas with no explanation.
     */
    public Path write(String name, String title) {
        final Path file = Runs.file(name);
        try (BufferedWriter out = Files.newBufferedWriter(file, StandardCharsets.UTF_8)) {
            // Two globals, not one. The named global is what a viewer built around a single
            // film reads; the array is what a viewer showing several of them reads, and a page
            // that includes two film files ends up with two entries in it without either file
            // knowing the other exists.
            out.write("window.NEOFIZ_FILM = ");
            out.write(json(title));
            out.write(";\n");
            out.write("(window.NEOFIZ_FILMS = window.NEOFIZ_FILMS || [])"
                    + ".push(window.NEOFIZ_FILM);\n");
        } catch (IOException e) {
            throw new UncheckedIOException("cannot write " + file, e);
        }
        return file;
    }

    /**
     * Marks the scene as a cross-section rather than a solid of revolution.
     *
     * <p>A viewer cannot tell from the mesh alone: the same rectangle of nodes is a tube wall
     * under one formulation and a flat bar under the other, and the difference is whether the
     * picture should be mirrored about r = 0. Getting it wrong draws a plane-strain bracket as
     * an impossible pair of brackets facing each other.
     */
    public Film planar() {
        this.axisymmetric = false;
        return this;
    }

    /** The film as JSON. Public so that a different viewer can be handed the same bytes. */
    public String json(String title) {
        final StringBuilder out = new StringBuilder(1 << 16);
        out.append("{\"title\":\"").append(escape(title)).append('"');
        out.append(",\"axisymmetric\":").append(axisymmetric);
        out.append(",\"nodes\":").append(mesh.nodeCount);
        out.append(",\"elements\":").append(mesh.elementCount);

        out.append(",\"conn\":[");
        for (int k = 0; k < mesh.conn.length; k++) {
            if (k > 0) out.append(',');
            out.append(mesh.conn[k]);
        }
        out.append(']');

        final int[] boundary = boundaryEdges();
        out.append(",\"boundary\":[");
        for (int k = 0; k < boundary.length; k++) {
            if (k > 0) out.append(',');
            out.append(boundary[k]);
        }
        out.append(']');

        out.append(",\"fields\":[");
        for (int f = 0; f < fields.size(); f++) {
            final Field field = fields.get(f);
            final double[] range = field.isAuto()
                    ? observedRange(f)
                    : new double[] {field.low(), field.high(), field.high()};
            if (f > 0) out.append(',');
            out.append("{\"name\":\"").append(escape(field.name()))
                    .append("\",\"units\":\"").append(escape(field.units()))
                    .append("\",\"low\":").append(number(range[0]))
                    .append(",\"high\":").append(number(range[1]))
                    .append(",\"peak\":").append(number(range[2]))
                    .append(",\"auto\":").append(field.isAuto()).append('}');
        }
        out.append(']');

        out.append(",\"readouts\":[");
        for (int m = 0; m < readouts.size(); m++) {
            if (m > 0) out.append(',');
            out.append('"').append(escape(readouts.get(m))).append('"');
        }
        out.append(']');

        out.append(",\"frames\":[");
        for (int i = 0; i < frames.size(); i++) {
            if (i > 0) out.append(',');
            final Frame frame = frames.get(i);
            out.append("{\"t\":").append(number(frame.time()));
            append(out, ",\"r\":", frame.r());
            append(out, ",\"z\":", frame.z());
            out.append(",\"v\":[");
            for (int f = 0; f < frame.values().length; f++) {
                if (f > 0) out.append(',');
                append(out, "", frame.values()[f]);
            }
            out.append(']');
            append(out, ",\"s\":", frame.scalars());
            if (frame.live() != null) {
                out.append(",\"live\":\"").append(frame.live()).append('"');
            }
            out.append('}');
        }
        out.append("]}");
        return out.toString();
    }

    /**
     * Fraction of samples allowed to sit above the top of an automatic ramp.
     *
     * <p>Not a rendering nicety. A single element can hold a value far outside what the rest of
     * the body ever reaches -- the one cell where a corner met a rigid plane carried 7.07 of
     * plastic strain in a run whose other 695 elements all stayed under 1.0 -- and a ramp
     * stretched to that maximum paints the entire body at the bottom of its scale. The picture
     * then says nothing happened, which is the opposite of what the film recorded.
     *
     * <p>Half a per cent is chosen against the thing that must survive: a localised band. A
     * band ten elements wide, present for the whole film, is about one and a half per cent of
     * the samples, so it is not clipped away -- and when a genuine band <em>is</em> wider than
     * this, clipping puts those elements at the top of the ramp, where they still read as the
     * hottest thing in the frame. What clipping costs is the ability to read an exact value off
     * the legend, which is why the true maximum is reported beside the range rather than
     * discarded.
     */
    private static final double RAMP_CLIP = 0.005;

    /**
     * The range an automatic ramp should span, and the largest value the field actually
     * reaches: {@code {low, high, peak}}.
     *
     * <p>{@code high} is the {@value #RAMP_CLIP} upper quantile rather than the maximum; see
     * {@link #RAMP_CLIP}. {@code peak} is the true maximum, so a viewer can say when the ramp
     * is clipped and by how much.
     *
     * <p>A field that never moves is given a little width, so the viewer divides by a range
     * rather than by zero and draws a uniform body instead of nothing at all.
     */
    private double[] observedRange(int field) {
        int n = 0;
        for (Frame frame : frames) {
            for (double v : frame.values()[field]) if (Double.isFinite(v)) n++;
        }
        if (n == 0) return new double[] {0.0, 1.0, 1.0};

        final double[] all = new double[n];
        int k = 0;
        for (Frame frame : frames) {
            for (double v : frame.values()[field]) if (Double.isFinite(v)) all[k++] = v;
        }
        java.util.Arrays.sort(all);

        final double low = all[0];
        final double peak = all[n - 1];
        // Counted as samples allowed above the top rather than as a quantile position, so that
        // a film with a handful of samples in it clips nothing at all: half a per cent of ten
        // values is zero values, and a rule that rounded that up to one would throw away the
        // maximum of every small film to guard against an outlier that cannot be there.
        final int above = (int) Math.floor(n * RAMP_CLIP);
        double high = all[Math.max(0, n - 1 - above)];
        if (high <= low) high = peak;

        if (high - low < 1e-12 * Math.max(1.0, Math.abs(high))) {
            return new double[] {low, low + Math.max(1e-9, Math.abs(low) * 1e-6), peak};
        }
        return new double[] {low, high, peak};
    }

    /**
     * Node pairs for every edge that belongs to exactly one element: the free surface.
     *
     * <p>A viewer needs this to draw an outline, and drawing an outline is not decoration --
     * it is what separates a body from a grid of coloured boxes. Computed from the
     * connectivity rather than from how a particular mesh generator happens to order its
     * nodes, so it keeps working for shapes {@link QuadMesh} cannot yet make.
     *
     * <p>An interior edge is shared by two elements and appears twice; a surface edge appears
     * once. That is the whole algorithm, and it is also how a body with a hole in it, or two
     * bodies in one mesh, come out right without being special-cased.
     */
    private int[] boundaryEdges() {
        return mesh.surface().edges();
    }

    private static void append(StringBuilder out, String key, double[] values) {
        out.append(key).append('[');
        for (int i = 0; i < values.length; i++) {
            if (i > 0) out.append(',');
            out.append(number(values[i]));
        }
        out.append(']');
    }

    /**
     * Six significant figures with the trailing zeros and the exponent's padding removed,
     * which is worth the trouble only because it is done a few million times per film.
     */
    private static String number(double value) {
        if (!Double.isFinite(value)) return "0";
        if (value == 0.0) return "0";
        String text = String.format(Locale.ROOT, "%.6g", value);
        if (text.contains("e")) {
            final String[] parts = text.split("e");
            return trim(parts[0]) + "e" + Integer.parseInt(parts[1]);
        }
        return text.contains(".") ? trim(text) : text;
    }

    private static String trim(String text) {
        int end = text.length();
        while (end > 0 && text.charAt(end - 1) == '0') end--;
        if (end > 0 && text.charAt(end - 1) == '.') end--;
        return text.substring(0, end);
    }

    private static String escape(String text) {
        return text.replace("\\", "\\\\").replace("\"", "\\\"");
    }
}
