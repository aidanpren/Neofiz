package org.neofiz.render;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.neofiz.core.Formulation;
import org.neofiz.mesh.QuadMesh;
import org.neofiz.solver.ExplicitSolver;
import org.neofiz.solver.Integration;
import org.neofiz.solver.Kinematics;
import org.neofiz.validate.BurstCase;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The frame capture.
 *
 * <p>A film is the one artefact in this project that nothing downstream can check for itself:
 * a plot with a wrong axis looks wrong, and a picture of a body with the wrong nodes in it
 * just looks like a body. So the properties that make it a faithful record -- deformed
 * positions rather than reference ones, a boundary that really is the free surface, and a
 * colour ramp that spans what the run produced -- are asserted here rather than eyeballed
 * there.
 */
class FilmTest {

    private static ExplicitSolver solver(QuadMesh mesh) {
        return new ExplicitSolver(mesh, BurstCase.COPPER, Formulation.AXISYMMETRIC,
                Integration.REDUCED, Kinematics.FINITE_STRAIN, 1.0, 0.5);
    }

    /** A four-by-three annulus: 12 elements, and a boundary anyone can count by hand. */
    private static QuadMesh wall() {
        return QuadMesh.cylinderWall(24.0e-3, 26.0e-3, 3.0e-3, 4, 3);
    }

    @Test
    @DisplayName("the free surface is every edge that belongs to one element")
    void boundaryIsTheFreeSurface() {
        // A 4-by-3 annulus has 3 elements' worth of edge on the bore, 3 on the outside, and
        // 4 on each end: 14 edges. Counting it by hand is the point -- the alternative is
        // trusting a traversal that happens to agree with itself.
        Film film = new Film(wall());
        ExplicitSolver s = solver(wall());
        film.capture(s);

        String json = film.json("t");
        int[] boundary = ints(json, "\"boundary\":[");
        assertEquals(2 * 14, boundary.length,
                "expected 14 surface edges, got " + boundary.length / 2);

        // And no edge appears twice, which is what would happen if interior edges leaked in.
        java.util.Set<Long> seen = new java.util.HashSet<>();
        for (int i = 0; i < boundary.length; i += 2) {
            long key = (long) Math.min(boundary[i], boundary[i + 1]) << 32
                    | Math.max(boundary[i], boundary[i + 1]);
            assertTrue(seen.add(key), "an edge was reported twice");
        }
    }

    @Test
    @DisplayName("positions are where the body is, not where it started")
    void capturesDeformedGeometry() {
        // The whole reason a film exists. A capture that stored reference coordinates would
        // produce a perfectly plausible movie of a tube doing nothing.
        QuadMesh mesh = wall();
        ExplicitSolver s = solver(mesh);
        Film film = new Film(mesh);
        film.capture(s);

        s.setPressureRamp(40.0e6, 1.0e-6);
        s.run(400);
        film.capture(s);

        String json = film.json("t");
        assertEquals(2, count(json, "{\"t\":"), "expected two frames");

        double[] first = doubles(json, "\"r\":[");
        double[] second = doublesAfter(json, "\"r\":[", 1);
        boolean moved = false;
        for (int i = 0; i < first.length; i++) {
            if (Math.abs(second[i] - first[i]) > 1e-9) moved = true;
        }
        assertTrue(moved, "the tube was pressurised and the capture did not notice");
    }

    @Test
    @DisplayName("an auto ramp spans what the run actually produced")
    void autoRangeFitsTheData() {
        QuadMesh mesh = wall();
        ExplicitSolver s = solver(mesh);
        Film film = new Film(mesh)
                .colourBy("fixed", "", 0.0, 10.0, (solver, e) -> e)
                .colourBy("fitted", "", Film.AUTO, Film.AUTO, (solver, e) -> e);
        film.capture(s);

        String json = film.json("t");
        assertTrue(json.contains("\"name\":\"fixed\",\"units\":\"\","
                        + "\"low\":0,\"high\":10,\"peak\":10,\"auto\":false"),
                "a fixed range must be passed through untouched:\n" + json);

        // Element index as the value, so the fitted range has to come back 0 to 11 -- the
        // largest index, not the declared ceiling of the field beside it. Twelve samples are
        // too few for the outlier clip to remove any of them, so the top of the ramp is the
        // maximum and the reported peak agrees with it.
        assertTrue(json.contains("\"name\":\"fitted\",\"units\":\"\","
                        + "\"low\":0,\"high\":11,\"peak\":11,\"auto\":true"),
                "the fitted range should reach the largest element index:\n" + json);
    }

    @Test
    @DisplayName("one wild element cannot flatten the ramp for every other one")
    void autoRangeClipsOutliers() {
        // The failure this exists to stop: a single cell where a corner met a rigid plane held
        // seven of plastic strain while the other 695 stayed under one, and the picture came
        // back uniformly black. The peak still has to be reported, or the viewer cannot say
        // that the hot end means "at least".
        // 400 elements over 10 frames is 4000 samples, of which the one wild element is 10 --
        // a quarter of a per cent, which is what makes it an outlier rather than a feature.
        // A band wide enough to exceed the allowance is meant to survive, and does: clipping
        // leaves it at the top of the ramp rather than removing it from the picture.
        QuadMesh mesh = QuadMesh.cylinderWall(24.0e-3, 26.0e-3, 3.0e-3, 40, 10);
        ExplicitSolver s = solver(mesh);
        Film film = new Film(mesh)
                .colourBy("spiky", "", Film.AUTO, Film.AUTO, (v, e) -> e == 0 ? 500.0 : e);
        for (int i = 0; i < 10; i++) film.capture(s);

        String json = film.json("t");
        double high = number(json, "\"name\":\"spiky\"", "\"high\":");
        double peak = number(json, "\"name\":\"spiky\"", "\"peak\":");
        assertEquals(500.0, peak, 0.0, "the true maximum must survive:\n" + json);
        assertTrue(high < 450.0, "the ramp should not stretch to the outlier, got " + high);
        assertTrue(high >= 380.0, "the ramp must still cover the ordinary values, got " + high);
    }

    private static double number(String json, String after, String key) {
        int at = json.indexOf(key, json.indexOf(after));
        int end = at + key.length();
        while (end < json.length() && ",}]".indexOf(json.charAt(end)) < 0) end++;
        return Double.parseDouble(json.substring(at + key.length(), end));
    }

    @Test
    @DisplayName("a field that never moves gets a range rather than a division by zero")
    void flatFieldIsSurvivable() {
        // An unyielded body has zero plastic strain everywhere, which is the first frame of
        // every run this will ever be pointed at.
        QuadMesh mesh = wall();
        Film film = new Film(mesh).colourBy("flat", "", Film.AUTO, Film.AUTO, (s, e) -> 0.0);
        film.capture(solver(mesh));

        String json = film.json("t");
        assertFalse(json.contains("NaN"), "a flat field produced a NaN range:\n" + json);
        assertFalse(json.contains("Infinity"), "a flat field produced an infinite range");
    }

    @Test
    @DisplayName("readouts and fields keep the order they were declared in")
    void orderIsDeclarationOrder() {
        // The viewer addresses both by index, so a reordering here would silently relabel
        // every number on the screen.
        QuadMesh mesh = wall();
        Film film = new Film(mesh)
                .colourBy("alpha", "", 0.0, 1.0, (s, e) -> 0.0)
                .colourBy("beta", "MPa", 0.0, 1.0, (s, e) -> 0.0)
                .readout("first", s -> 1.0)
                .readout("second", s -> 2.0);
        film.capture(solver(mesh));

        String json = film.json("t");
        assertTrue(json.indexOf("\"alpha\"") < json.indexOf("\"beta\""), json);
        assertTrue(json.indexOf("\"first\"") < json.indexOf("\"second\""), json);
        assertTrue(json.contains("\"s\":[1,2]"), "readouts in order:\n" + json);
    }

    @Test
    @DisplayName("the title survives a quote in it")
    void escapesTitles() {
        Film film = new Film(wall());
        film.capture(solver(wall()));
        assertTrue(film.json("a \"quoted\" title").startsWith("{\"title\":\"a \\\"quoted\\\" title\""));
    }

    // ------------------------------------------------------------------ helpers

    private static int count(String text, String needle) {
        int n = 0, at = 0;
        while ((at = text.indexOf(needle, at)) >= 0) { n++; at += needle.length(); }
        return n;
    }

    private static String slice(String json, String key, int occurrence) {
        int at = -1;
        for (int i = 0; i <= occurrence; i++) at = json.indexOf(key, at + 1);
        int start = at + key.length();
        return json.substring(start, json.indexOf(']', start));
    }

    private static int[] ints(String json, String key) {
        String body = slice(json, key, 0);
        if (body.isEmpty()) return new int[0];
        String[] parts = body.split(",");
        int[] out = new int[parts.length];
        for (int i = 0; i < parts.length; i++) out[i] = Integer.parseInt(parts[i].trim());
        return out;
    }

    private static double[] doubles(String json, String key) {
        return doublesAfter(json, key, 0);
    }

    private static double[] doublesAfter(String json, String key, int occurrence) {
        String[] parts = slice(json, key, occurrence).split(",");
        double[] out = new double[parts.length];
        for (int i = 0; i < parts.length; i++) out[i] = Double.parseDouble(parts[i].trim());
        return out;
    }
}
