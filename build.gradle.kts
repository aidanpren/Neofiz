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
    testLogging {
        events("passed", "failed", "skipped")
        showStandardStreams = true
    }
}
