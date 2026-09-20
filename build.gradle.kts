plugins {
    id("java")
    id("application")
}

group = "org.neofiz"
version = "0.1.0-M0"

repositories {
    mavenCentral()
}

// Bytecode targets 21; compiled by whatever JDK runs the build. No toolchain block, so
// Gradle will not go hunting for a specific JDK install.
tasks.withType<JavaCompile>().configureEach {
    options.release = 21
}

dependencies {
    testImplementation(platform("org.junit:junit-bom:5.11.4"))
    testImplementation("org.junit.jupiter:junit-jupiter")
    testRuntimeOnly("org.junit.platform:junit-platform-launcher")
}

application {
    mainClass = "org.neofiz.M0Cylinder"
}

tasks.test {
    useJUnitPlatform()
    // Forward -Dneofiz.* through to the test JVM, so a gate can be re-run at a different
    // setting from the command line rather than by editing it and putting it back.
    System.getProperties().forEach { k, v ->
        if (k.toString().startsWith("neofiz.")) systemProperty(k.toString(), v.toString())
    }
    testLogging {
        events("passed", "failed", "skipped")
        showStandardStreams = true
    }
}

// The pressure-expansion curve of the reference tube, plotted in the terminal and written
// to runs/burst-curve.csv. Separate from `run` because it is an output, not a gate.
tasks.register<JavaExec>("burstCurve") {
    group = "application"
    description = "Traces the reference tube's pressure-expansion curve through burst."
    mainClass = "org.neofiz.BurstCurve"
    classpath = sourceSets["main"].runtimeClasspath
}

// Two sweeps and the exponents they produce: the M3 exit criterion in miniature.
tasks.register<JavaExec>("barlowSweep") {
    group = "application"
    description = "Sweeps wall thickness and bore diameter, and fits the scaling laws."
    mainClass = "org.neofiz.BarlowSweep"
    classpath = sourceSets["main"].runtimeClasspath
}

// The design map: wall thickness against applied pressure, cells coloured by how often the
// tube gives way. --args="<thicknesses> <samples per thickness>" to resize it.
tasks.register<JavaExec>("designMap") {
    group = "application"
    description = "Maps wall thickness against applied pressure, with the defect scatter."
    mainClass = "org.neofiz.DesignMap"
    classpath = sourceSets["main"].runtimeClasspath
}

// Captures a defect tube's burst frame by frame into runs/film.js, for a viewer to play.
tasks.register<JavaExec>("burstFilm") {
    group = "application"
    description = "Films a pressurised tube giving way. --args=\"<seed>\" picks the tube."
    mainClass = "org.neofiz.BurstFilm"
    classpath = sourceSets["main"].runtimeClasspath
}

// The first scene that is not a solid of revolution: plane strain, a polygon outline and
// two materials in one mesh. --args="<speed in m/s>".
tasks.register<JavaExec>("impact") {
    group = "application"
    description = "Throws a two-material L bracket at an anvil, in plane strain."
    mainClass = "org.neofiz.Impact"
    classpath = sourceSets["main"].runtimeClasspath
}

// Thirteen bodies in one solve: a steel slug through a stacked copper wall. The first
// scene here that is a scene rather than a specimen. --args="<speed in m/s>".
tasks.register<JavaExec>("smash") {
    group = "application"
    description = "Fires a slug through a stacked wall, in plane strain with contact."
    mainClass = "org.neofiz.Smash"
    classpath = sourceSets["main"].runtimeClasspath
}

// Where the time goes: element kernel, damage, wall, contact, threads, and what a second of
// simulated time would cost at each mesh size.
tasks.register<JavaExec>("bench") {
    group = "application"
    description = "Times the solver piece by piece, in element-steps per second."
    mainClass = "org.neofiz.Bench"
    classpath = sourceSets["main"].runtimeClasspath
}

// A leaning stack left alone under gravity, which mass scaling is what makes affordable.
tasks.register<JavaExec>("topple") {
    group = "application"
    description = "Lets a leaning stack of blocks fall over, in plane strain with contact."
    mainClass = "org.neofiz.Topple"
    classpath = sourceSets["main"].runtimeClasspath
}

// Reads the run log back. --args="<program>" narrows it.
tasks.register<JavaExec>("history") {
    group = "application"
    description = "Prints every run recorded in runs/log.csv."
    mainClass = "org.neofiz.History"
    classpath = sourceSets["main"].runtimeClasspath
}
