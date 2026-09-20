package org.neofiz;

import org.neofiz.report.RunLog;

import java.time.ZoneId;
import java.time.format.DateTimeFormatter;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * Every run that has been recorded, most recent last.
 *
 * <p>The companion to {@link org.neofiz.report.RunLog}: something has to read the file back or
 * writing it is a ritual. Filtered by program name if one is given.
 *
 * <pre>  ./gradlew history --args="smash"</pre>
 */
public final class History {

    private static final DateTimeFormatter WHEN =
            DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm").withZone(ZoneId.systemDefault());

    private History() {
    }

    public static void main(String[] args) {
        final String program = args.length > 0 ? args[0] : null;
        final List<RunLog.Entry> entries =
                program == null ? RunLog.all() : RunLog.historyOf(program);

        System.out.println("NEOFIZ  run log");
        System.out.println("=".repeat(76));
        if (entries.isEmpty()) {
            System.out.println(program == null
                    ? "Nothing recorded yet. Run a scene and it will log itself."
                    : "No runs of " + program + " recorded.");
            return;
        }

        for (RunLog.Entry e : entries) {
            System.out.printf(Locale.ROOT, "%s  %s%n", WHEN.format(e.when()), e.program());
            System.out.printf(Locale.ROOT, "  in   %s%n", pairs(e.parameters()));
            System.out.printf(Locale.ROOT, "  out  %s%n", pairs(e.results()));
        }
        System.out.printf(Locale.ROOT, "%n%d run%s recorded%n",
                entries.size(), entries.size() == 1 ? "" : "s");
    }

    private static String pairs(Map<String, String> map) {
        if (map.isEmpty()) return "(none)";
        final StringBuilder out = new StringBuilder();
        for (Map.Entry<String, String> e : map.entrySet()) {
            if (out.length() > 0) out.append("   ");
            out.append(e.getKey()).append(' ').append(e.getValue());
        }
        return out.toString();
    }
}
