# Neofiz

A physics sandbox with no scripted objects. Materials, geometry, and physical parameters —
the content is the physics.

Built against the *Everything Is Parameters* charter and the
[build plan](https://claude.ai/code/artifact/fb241530-93c1-45ee-a03b-47dd540ecae2).

## Where this is

**M0 complete. Both halves of M1 are met** — the Taylor mushroom mesh-converged at a
physically plausible friction coefficient, and hydrostatic burst to 0.28 % against a 10 %
criterion. Two-dimensional Lagrangian explicit FEM, 4-node quads, central-difference
integration, full or reduced quadrature, J2 plasticity with combined isotropic and kinematic
hardening, Johnson–Cook flow stress with adiabatic thermal softening, updated Lagrangian
kinematics with an objective stress update, a rigid anvil with a Coulomb cone, a threaded
kernel that is bit-for-bit identical at any thread count, a mesh-independent Weibull defect
field, and progressive softening regularised by fracture energy per unit crack area. 313 tests
green.

**You can now watch a run instead of reading its summary.** `./gradlew burstFilm` keeps the
deformed geometry at every sample and writes it out; `viewer/burst-viewer.html` plays it back
with a scrub track that is the wall-capacity curve, the burst marked on it. Two things the
first render made obvious that no number had. A true-scale section of this tube is unreadable
— the wall is 4 % of the diameter, so it draws as two hairlines around an enormous hole, and
the bore has to be cut away the way a sectional drawing would. And **there is no crack**: the
bulge ratio is 1.009 and damage at the last frame runs 0.77 to 1.00 across the *whole* wall,
so the tube softens almost uniformly rather than necking. That is the parallel-rings result
and the non-converging burst pressure, seen directly for the first time.

**The first scene that is not a solid of revolution.** `./gradlew impact` throws an L-shaped
bracket — steel foot, copper leg, welded at the interface — onto a rigid anvil in plane
strain, and `viewer/scene-viewer.html` plays it. Three restrictions had to go to make that
sentence expressible: one material per mesh, two hardcoded mesh generators, and a solver
whose `PLANE_STRAIN` flag had existed since M0 with nothing built on it. None of them was the
physics. The constitutive code is untouched, and a single-material mesh now reaches the
kernel through `Materials` **bit-for-bit** identically to the way it did through a constant,
which is what lets every gate above keep testing what it was written to test.

**Bodies can now touch each other.** `./gradlew smash` fires a steel slug through a stacked
wall of twelve copper blocks — thirteen separate bodies in one solve, resting on one another
and on the floor, with a Coulomb cone between every pair. That took penalty contact
(`Contact`), a mesher that can produce bodies which *touch* rather than weld
(`Outline.assemble`), and gravity, which the solver had never had. Momentum through a
collision is exact to round-off; energy closes to a few per cent, and the whole shape of that
number is written up below.

**And things now come apart.** An element softened to nothing is deleted and the free surface
closes around the hole, so the slug perforates the wall instead of pushing a permanent dent
ahead of it — 32 of 652 elements at 700 m/s, none of them inverted. Deletion is bookkeeping
rather than physics, and the division matters: the *softening* that gets an element to zero
strength is regularised by fracture energy and is mesh-independent; deletion only disposes of
something that can no longer be integrated. Mass and nodes stay, so momentum does not move,
and the elastic energy that leaves with the debris is reported rather than dropped.

**And things can now simply fall over.** `./gradlew topple` leans five 40 mm blocks into a
staircase and lets go: nothing is thrown, gravity does all of it, and 320 ms of physics costs
4 s of wall clock. That was unaffordable a day ago — an explicit solve runs in microseconds
and gravity works in milliseconds, so every scene before this one had to *begin* already
stacked. `setMassScaling` buys the step back, and the lattice mesher makes the cheap form of
it the right form: on a mesh where every element is the same size, uniform scaling **is**
selective scaling.

Every conservation audit on the Taylor case is exact or converging. Three findings shaped it,
and none of them were the one expected going in:

- **Mushroom diameter is an interface measure, not a constitutive one.** It runs from 10.13 mm
  frictionless to 9.04 mm at μ = 0.20 with the same material, and crosses the measured 9.5 mm
  at **μ ≈ 0.09** — an ordinary coefficient for steel on steel.
- **It converges at first order, not the element's second**, because it is read at the contact
  edge, which is a geometric singularity. That is why it needed Richardson extrapolation
  rather than a finer mesh — and getting it converged roughly *doubled* the fitted μ, out of
  the implausible 0.047 an un-converged mesh had suggested.
- **Johnson–Cook did not do what it was added to do**: rate hardening (+15 %) and adiabatic
  thermal softening (−14 %) very nearly cancel here. See
  [gate M1c](#gate-m1c--taylor-impact).

Final length, by contrast, is converged — one part in ten thousand between the two finest
meshes — as well as rule-independent and interface-insensitive, at L/L₀ = 0.9079. It is the
clean read on the flow stress; the mushroom never was.

**Burst was the opposite experience: it converged immediately and the difficulty was
measuring it at all.** Past the load maximum there is no equilibrium, so the applied pressure
cannot be the measurement; the wall's own capacity is, read off an exact equilibrium identity
that is defined on both sides of the peak. What is left over — 0.28 % — turns out to be two
elastic terms of opposite sign that a rigid-plastic closed form omits by construction, and
they predict it to 0.004 % from E and ν alone. See [gate M1d](#gate-m1d--hydrostatic-burst).

**Gate 2 is now measured rather than projected, and the projection it replaces was wrong by a
factor of two in the flattering direction.** The worked example runs in 15.8 s on eight cores,
not the 4.1 s this file used to claim — and the 2–10 s it was being compared against is the
plan's wall clock for *one GPU*, so a CPU was never supposed to reach it. The cost model
itself survives intact: seconds, not minutes. See [gate 2](#gate-2--the-cost-model).

**M2 is complete, and the two most interesting things in it were not in the plan.** The defect
field is mesh-independent and its size effect matches weakest-link theory to 0.2 %; the solver
reads it; the same tube now bursts at six pressures in six places. Along the way:

- **Volume normalisation on its own makes the part *stronger* without bound as the mesh
  refines** — 35 % over a thousandfold — because below the correlation length neighbouring
  elements are not independent chances. The plan names one trap here; this is a second one,
  facing the flattering way, which survives any test phrased as "must not get weaker".
- **A pressurised shell averages its defects rather than failing at its weakest link.** Over an
  eightfold change in tube length the mean burst pressure is flat while the scatter falls as
  L^(−0.40) — neither of the two things Weibull theory predicts. The governing length is the
  shell's shear-lag length `√(Rt)`, and pushing the field's correlation length through it turns
  the behaviour over exactly as that argument says. Along its axis a tube is a bundle of
  parallel rings, not a chain.
- **The burst pressure does not converge under mesh refinement**: +0.22 % per doubling,
  monotone, heading for the defect-free answer. Crack-band regularisation makes the fracture
  *energy* exact to 1e-10 and it does that by making the softening *modulus* mesh-dependent,
  which is only right once the band has localised. Here the load maximum arrives first. The
  number is reported rather than smoothed.

See [M2](#m2--the-defect-field) and
[M2's second half](#m2s-second-half--what-the-solver-does-with-it).

**M3 has started, and the first sweep already answers its exit criterion.** The burst
harness no longer throws its curve away, and twelve tubes swept across a decade of wall
thickness return `P ∝ t^0.9989` with a worst residual of 0.02 % — Barlow's formula, measured
off the solver, which was never told what Barlow's formula is. The bore-diameter sweep is the
more interesting of the two, because it comes back *wrong*: `d^-0.9666`, with 1.3 % of
structure the exponent is flattening out. Plot the identical runs against the **mean**
diameter and the bend disappears and the exponent snaps to −0.9989. Equilibrium is written on
the mean radius; a wall does not know what its bore diameter is. Nobody was told that either —
it arrived as curvature in a residual plot. See [the first sweep](#the-first-sweep--barlow-measured).

**And the design map exists.** 108 tubes, 29 seconds, and the boundary between *holds* and
*bursts* comes out as a band rather than a line — which is the image the build plan calls the
single most valuable one in the product. Two things about it were not expected. The band is
**narrower than a pixel** on axes wide enough to show the whole design range, because the
scatter is 0.8 % and the range is a factor of three; it only resolves once the pressure axis
is divided through by each column's own defect-free burst pressure. And once it resolves it
**slopes**: 6.90 % knockdown at a 1.25 mm wall against 7.75 % at 4.00 mm. That is a size
effect — more material, more chances at a bad patch — and it is *ten times weaker* than
weakest-link theory's `V^(1/m)` asks for over the same range. Which is the M2 parallel-rings
result, arrived at from a different experiment and not looked for. See
[the design map](#the-design-map).

```
./gradlew run           # the M0 report
./gradlew test          # the validation gates
./gradlew burstCurve    # the pressure-expansion curve, through the maximum
./gradlew barlowSweep   # two sweeps, and the exponents they produce
./gradlew designMap     # wall thickness against applied pressure, with scatter
./gradlew burstFilm     # captures a burst frame by frame into runs/film.js
./gradlew impact        # a two-material bracket on an anvil, into runs/bracket.js
./gradlew smash         # a slug through a stacked wall, into runs/smash.js
./gradlew topple        # a leaning stack falls over, into runs/topple.js
./gradlew bench         # where the time goes, in element-steps per second
./gradlew history       # every run recorded, from runs/log.csv
node viewer/serve.js    # plays either film at http://localhost:8731
```

### Gate 1 — correctness

Elastic axisymmetric FEM against Lamé's thick-walled cylinder. Both end conditions against
both quadrature rules, at 0.5 mm elements:

| End condition | Rule | Displacement error | Hoop stress error |
| --- | --- | --- | --- |
| Restrained (ε_z = 0, plane strain) | full 2×2 | 0.0074 % | 0.0013 % |
| Restrained | reduced + hourglass | 0.0000 % | 0.0057 % |
| Open (σ_z = 0, plane stress) | full 2×2 | 0.0124 % | 0.0374 % |
| Open | reduced + hourglass | 0.0000 % | 0.0057 % |

Gate was 1 %. Stress error falls ~4× per mesh halving, which is the second-order
convergence a bilinear quad should show.

The exact zeros are real, not a formatting artefact: reduced integration is **nodally
exact** here, to seven digits at every node, at every ν, and out to a 10:1 wall ratio where
full integration is 30 % off. That is the classical Galerkin result for a one-dimensional
two-point boundary value problem, which is what this cylinder is. Stresses are still only
second-order accurate — exact nodes do not mean exact gradients.

**Two end conditions, deliberately.** The Lamé stresses contain neither E nor ν, so σ_r and
σ_θ are identical under plane stress and plane strain. A stress-only harness cannot detect a
wrong end restraint — it reports a clean pass while solving a different problem. This
happened during M0: the first run matched hoop stress to 0.04 % while displacement was 7 %
off, because the plane-strain constraint was documented but never applied. Only the
displacement check caught it.

### Volumetric locking — why reduced integration is not optional

Displacement error as ν → 0.5, restrained ends, fixed coarse 3-element wall:

| ν | full 2×2 | reduced + hourglass |
| --- | --- | --- |
| 0.290 | 0.082 % | 0.000 % |
| 0.450 | 0.321 % | 0.000 % |
| 0.490 | 1.550 % | 0.000 % |
| 0.499 | 13.537 % | 0.000 % |

Full integration locks — two orders of error growth on a fixed mesh. Reduced integration is
flat.

This matters beyond elasticity. **Plastic flow is incompressible**, so the tangent response
of a fully yielded element approaches exactly this limit. Reduced integration entered the
plan as a 4× speed lever; it is really a correctness requirement for M1. Full integration
would have locked inside the plastic zone of the Taylor cylinder, where no closed form is
standing by to catch it.

Finding it required fixing a different bug first. The relaxation schedule timed its pressure
ramp and its damping off the **dilatational** wave speed, which runs away as ν → 0.5 while
the ring's actual breathing mode — an extensional, bar-wave mode — barely moves. At ν = 0.499
the ramp had collapsed to a fraction of a real period and the damping was ~13× too high, so
the run simply had not settled, and it reported 64 % error that looked exactly like locking.
KE/SE was 2.6e-2 where a settled run reaches 1e-20. The audit is the only reason the two
were distinguishable. That threshold has since been tightened from 1e-3 to 1e-9.

### Gate M1a — J2 plasticity

Elastic-plastic thick-walled cylinder against the closed form. Thick wall, b/a = 2, σ_y =
800 MPa, elastic-perfectly-plastic, 20 elements through the wall:

| p / p_elastic | Front (FEM / exact) | Yield condition | σ_r | σ_θ |
| --- | --- | --- | --- | --- |
| 1.10 | 26.25 / 26.31 mm | 0.134 % | 0.072 % | 0.107 % |
| 1.35 | 30.00 / 30.19 mm | 0.245 % | 0.061 % | 0.227 % |
| 1.60 | 36.25 / 35.60 mm | 0.727 % | 0.051 % | 0.611 % |

*(reduced integration; full integration is comparable — see `./gradlew run`)*

A far sharper instrument than Lamé. The elastic gate can be passed by any code that
assembles a stiffness matrix correctly. This one needs the yield criterion, the flow rule,
the return map and equilibrium to be right *simultaneously*, and it is sensitive to the
2/√3 that separates von Mises from Tresca.

The **yield-condition** column is the one that matters: inside the plastic zone
σ_θ − σ_r must equal 2σ_y/√3 as an identity, not an approximation. It is the return map's own
consistency condition, read back out of the assembled solution rather than checked at a
material point.

Two things this exposed that are worth knowing rather than smoothing over:

- **Yield onset carries an O(h) bias, always non-conservative.** A centroid-sampled element
  cannot see the stress peak at the bore, only the value half an element inside it, so the
  discrete cylinder yields *late* — by exactly ((a + h/2)/a)², which is 10 % high at 10
  elements and 2.5 % at 40. It converges away; it is measured rather than tolerated.
- **Above the limit pressure the case refuses to run.** A fully plastic cylinder has no
  static equilibrium to relax onto. Returning a stress field anyway is the dangerous
  outcome: it looks entirely reasonable and the kinetic energy grows too slowly to be
  obvious. That was observed before the guard went in.

Thin walls are nearly useless for this. At the charter's b/a = 1.2 the *limit* pressure — where
the plastic front reaches the outer surface — is only **1.19×** the elastic limit, so there is
almost no contained-plastic regime between "entirely elastic" and "no equilibrium exists" in
which to validate anything. For an elastic-perfectly-plastic material that limit pressure is
also the burst pressure, and [gate M1d](#gate-m1d--hydrostatic-burst) reaches it from an
entirely separate derivation; thin walls are where *that* case wants to live, for the same
reason they are useless here.

### Gate M1b — finite strain

A stressed element rotated rigidly. No stretching anywhere, so the equivalent stress must not
move at all. Relative change:

| Rotation | Finite strain | Small strain |
| --- | --- | --- |
| 0.25 turn | 2.7e-13 | **4.99e+02** |
| 0.50 turn | 6.3e-14 | **9.99e+02** |
| 1.00 turn | 5.5e-13 | 2.1e-12 ← |
| 3.30 turn | 1.6e-13 | **6.54e+02** |

Machine precision at every angle, independent of step count. Objectivity has no closed form
to miss and no convergence study to run — it either holds exactly or the update is
manufacturing stress out of rotation.

**Why it's exact rather than merely accurate.** The rotation increment is the **Cayley
transform**, not the exponential, and the velocity gradient is taken at the **mid-step**
configuration. In that pairing, rigid motion produces an exactly skew velocity gradient, so
the symmetric part — the strain — vanishes identically rather than to some order in the step.
Either piece alone leaks a spurious strain of order φ² per step: small, one-signed, and
accumulating over a long spin into stress nobody applied. In the axisymmetric half-plane the
spin has one component, so this costs two multiplies and a divide, with no trigonometry.

**The 1.00-turn row is the warning.** The non-objective error is *periodic in rotation, not
monotonic*: it's 500× at a quarter turn and 1e-12 at a whole one, because the linearised
increments of a rigid rotation cancel over a complete revolution. An objectivity test that
happened to use 360° would certify a broken update as perfect. That row is pinned by a test
so the angle in the others doesn't get "tidied up" later.

**Large-deformation check with an exact answer.** Isochoric plane-strain compression to 50 %
height — squash in z by s, stretch in r by 1/s, so det F = 1 at every instant:

| | Measured | Exact |
| --- | --- | --- |
| Pressure | 0.0000 MPa | 0, since ln(det F) = 0 |
| σ_zz | −461.88 MPa | −461.88 = −σ_y/√3 |

Two independent things have to be right for both to land: the volumetric strain has to be
ln(det F) rather than a linearised trace, and plastic flow has to leak no volume over 50 %
deformation. Small strain fails the first.

Also in: **adaptive CFL** (elements that compress carry shorter wave transit, and holding dt
fixed while a mesh crushes is how an explicit run goes unstable late), and a **follower
pressure load** that acts on the surface where it now is, at the area it now has.

### Gate M1c — Taylor impact

A solid 4340 cylinder, 37.97 mm long and 7.595 mm in diameter, fired flat-on into a rigid
anvil at 181 m/s. Measured mushroom diameter 9.5 mm; published simulations 9.80–9.83 mm. No
damping, no relaxation, no pressure ramp — the first case where inertia is the physics rather
than something to settle out, and the first with no closed form behind it at all.

| Law | Rule | n | Mushroom | vs measured | L/L₀ | max ε_p | T max |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Johnson–Cook | reduced | 4 | 9.472 mm | −0.30 % | 0.9079 | 0.50 | 443 K |
| rate-indep. | reduced | 4 | 9.590 mm | +0.94 % | 0.8827 | 0.54 | — |
| Johnson–Cook | reduced | 8 | 9.776 mm | +2.90 % | 0.9081 | 0.67 | 495 K |
| rate-indep. | reduced | 8 | 9.855 mm | +3.74 % | 0.8843 | 0.72 | — |
| Johnson–Cook | reduced | 16 | 9.965 mm | +4.90 % | 0.9080 | 0.75 | 517 K |
| rate-indep. | reduced | 16 | 10.005 mm | +5.31 % | 0.8845 | 0.83 | — |
| Johnson–Cook | full 2×2 | 8 | 9.448 mm | −0.55 % | 0.9076 | 0.56 | 459 K |

#### Johnson–Cook did not do what it was added to do

The rate-independent surrogate over-predicted the mushroom and the obvious diagnosis was the
missing rate term: at the 10⁴–10⁵ s⁻¹ this test imposes, 4340 steel is much stronger than the
quasi-static curve, and Johnson–Cook's rate term is worth **+15 %** on flow stress. It is.
But the same deformation heats the specimen by **220 K** adiabatically — it is its own heat
source, 99 % of the impact energy arrives as plastic work — and thermal softening at that
temperature is worth **−14 %**.

```
rate hardening   1.1478
thermal softening 0.8551
product          0.9814     <- two per cent, in the other direction
```

So the mushroom moved by under 1 %. **Either term shipped alone would have looked like a
large effect and a large error**, which is why they were built together rather than in the
order that would have been easier to validate.

What did move is the length: L/L₀ goes 0.884 → 0.908, far outside the mesh convergence of
either model. That is the honest measure of what the better flow stress bought, and the
rate-independent surrogate is kept runnable so the comparison stays available rather than
becoming a claim about a version that no longer exists.

#### The mushroom is an interface measure

Coulomb friction sweep, n = 8, reduced integration:

| μ | Mushroom | vs measured | friction dissipation | L/L₀ |
| --- | --- | --- | --- | --- |
| 0 (frictionless) | 9.776 mm | +2.90 % | 0 | 0.9081 |
| 0.02 | 9.651 mm | +1.59 % | 0.46 % | 0.9082 |
| **0.05** | **9.476 mm** | **−0.25 %** | 1.02 % | 0.9085 |
| 0.10 | 9.259 mm | −2.54 % | 1.67 % | 0.9088 |
| 0.20 | 8.987 mm | −5.40 % | 2.21 % | 0.9093 |
| 0.50 | 8.921 mm | −6.09 % | 1.48 % | 0.9102 |
| ∞ (welded) | 8.928 mm | −6.02 % | 0 | 0.9106 |

The two limits are **a millimetre apart** — an order of magnitude more than Johnson–Cook was
worth — and the measurement is crossed at μ ≈ 0.05. Across the whole sweep the length ratio
moves 0.3 %. Length is set by the flow stress; the mushroom is set by the boundary condition
at the anvil. They are not the same kind of number, and a mushroom diameter quoted without its
friction assumption is not a result.

Two structural checks are visible in that table, and each would catch a plausible bug:

- **Dissipation vanishes at both ends and peaks between them.** At μ = 0 there is no
  tangential force; at the welded limit there is no sliding. Dissipation is the product of the
  two, so a monotonic curve would mean one factor isn't being applied.
- **μ = ∞ reaches the welded answer through the same cone test** every other coefficient uses,
  rather than down a pinned-DOF path beside it — and μ = 0.5 converges onto it from below.
  Same idiom as elasticity being plasticity with an infinite yield stress: the limit case
  falls out of the general arithmetic instead of being special-cased next to it.

#### Why the mushroom would not converge

The fitted μ drifted with mesh — 0.047 at n = 8, 0.073 at n = 16 — because the frictionless
baseline itself kept climbing. Three meshes at n = 8 / 16 / 32 say why:

| μ | n=8 | n=16 | n=32 | observed order | extrapolated | GCI |
| --- | --- | --- | --- | --- | --- | --- |
| 0 | 9.7758 | 9.9653 | 10.0545 | **1.09** | 10.134 mm | ±0.99 % |
| 0.10 | 9.2587 | 9.3551 | 9.4029 | **1.01** | 9.450 mm | ±0.63 % |
| 0.20 | 8.9872 | 9.0263 | 9.0373 | **1.83** | 9.042 mm | ±0.06 % |

**A bilinear quad gives second order on a smooth field. The mushroom gives first.** It is read
at the contact edge, where the free surface meets the anvil, and that corner is a geometric
singularity — no amount of refinement restores the missing order. The response is Richardson
extrapolation and an uncertainty band, not a bigger mesh.

The μ = 0.20 row is the confirmation. Grip the face and the order climbs back towards two
while the band collapses by a factor of sixteen, because material is no longer sliding past
the corner that was generating the singular field. The slow convergence is a property of
**sliding contact at an edge**, not of the mesh or of the mushroom as a feature.

It is also not the measurement's fault: the outer-surface profile near the face is smooth, and
the widest node sits only 0.01 mm above its neighbour at n = 32. There is no rogue folded
node — the whole near-face geometry is converging together, just slowly.

**Gate status: extrapolated frictionless mushroom 10.134 mm ± 1.0 %, which is +6.7 % against
the measured 9.5 mm. The extrapolated sweep crosses the measurement at μ ≈ 0.09** — confirmed
by a run at that value, which extrapolates to 9.511 mm, +0.11 % against measurement and well
inside its own ±0.69 % band.

That is roughly **double** the 0.047 the un-converged n = 8 mesh asked for, and it matters
which one you get: 0.09 is an ordinary coefficient for a steel-on-steel interface, and 0.047
is not. Converging the mesh did not merely tighten the fit — it moved the answer out of the
range where it would have had to be explained away.

#### The audits, which do hold

| | Value at finest mesh | |
| --- | --- | --- |
| Momentum balance | 9.6e-15 | exact |
| Energy balance | +0.245 % | closes |
| Plastic work | 98.9 % of impact energy | where it all goes |
| Hourglass | 1.2e-3 of impact energy | gated at 1e-2 |
| Arrival loss | 0.36 % of impact energy | O(h), halves per refinement |

**Momentum balance is exact, and it earned its keep.** Every element's axial internal forces
sum to zero — the shape function derivatives are the gradient of a partition of unity — so
nothing inside the mesh can move the body's axial momentum, and the anvil's impulse is its
entire history. That identity caught the one real bug in this work: arrival was implemented
as *land on the plane with the velocity that reaches it*, then *stop on the following step*,
and the momentum destroyed by the second half belonged to no impulse. It reads as 1e-3, and
it is **invisible at first impact** — a node already touching the anvil lands with zero
velocity anyway — appearing only once the contact patch spreads to nodes arriving from above.
Arrival is now one event: land exactly on the plane and stop.

**The energy balance needed a new term and a different denominator.** Plastic work is now
accumulated per Gauss point, because on an impact ~99 % of the input energy ends up there and
without it the balance is an identity between small residuals that closes whether or not the
run is sound. It is also what Johnson–Cook's thermal softening will read: the temperature
rise that does the softening is adiabatic, β·W_p/(ρc_p).

The hourglass audit needed the denominator changed outright. Hourglass over *recoverable
strain energy* is the right ratio for the quasi-static cylinder, where it runs at 1e-29. On
this case it reads **0.98** — not because hourglass control is carrying the answer, but
because by the time an impacted specimen has unloaded, strain energy is a small residual
rather than the scale of the problem. Against the energy that actually went in, hourglass
control holds 0.2 %. A wrong denominator here would have condemned a sound run.

### Contact

A rigid frictionless immovable plane. The plan's eventual contact model is penalty-based with
Coulomb friction, which is what deformable-on-deformable pairs need; for a wall that cannot
move, the kinematic form is not a simplification of that but strictly better. Non-penetration
is exact rather than a consequence of stiffness. There is no penalty frequency, so the stable
timestep is untouched — and with **no mass scaling** available to buy it back, a stiffness
that halves dt is a direct doubling of run cost. There is no stiffness parameter to sweep.

**A node in contact is in one of two states, and collapsing them costs an error that does not
converge.** A node *arriving* has its approach velocity destroyed in one step; each node
arrives once, so the cost is first order in element size and converges away. A node *resting*
is held by a reaction force that exactly cancels the internal force — it does not accelerate,
does not move, and because it does not move the reaction does no work, so nothing is lost.

Treating a resting node as an arrival every step is the obvious implementation. It discards
`dt²f²/2m` per node per step, which is O(h²) and sounds harmless, until it is summed over the
O(1/h) steps of a fixed physical time and the O(1/h) nodes of a refined contact face. The
total is **O(1): refining the mesh does not reduce it.** That version ran at 13 % of the
impact energy at every resolution, quietly removed from the plastic work forming the mushroom.
A test asserts the convergence *rate* rather than the magnitude, because the rate is what
tells the two apart.

Release needs no criterion beyond the sign of the force, so there is no contact set aged
across steps and no penetration tolerance.

#### Coulomb friction is a return map

The tangential half has exactly the structure of the plasticity kernel next door, for the
same reason — a constraint that holds until a threshold and then flows:

| | plasticity | friction |
| --- | --- | --- |
| trial | elastic step | assume the node sticks |
| surface | von Mises yield | the friction cone, `\|T\| ≤ μN` |
| inside | the stress stands | static friction is whatever it needs to be |
| return | scale onto the surface | cap at μN, opposing the slip |

The trial reaction is recovered from the velocity update the integrator *would* have done, so
the damping coefficients are carried rather than assumed away and the cone is tested against
the force that would actually have been applied.

The energy splits in two, and keeping them apart is the point. **Sliding dissipates, and that
is physical** — a frictional interface is supposed to absorb energy, and a balance omitting it
would not close. **Sticking destroys the tangential kinetic energy** of a node that was
sliding, in the same way and for the same reason arrival destroys its normal kinetic energy,
so it is charged to the artefact account instead. An audit that added the two could not tell a
frictional anvil from a leaky constraint.

The contact state is settled once per node, before either component is touched — the axial
branch rewrites the contact set, so asking again afterwards would apply friction using this
step's answer for the normal direction and the previous step's for the tangential one.

### Gate M1d — hydrostatic burst

An OFHC copper tube, mean radius 25 mm, wall 2 mm, pressurised until it stops being able to
hold it. σ = 90 + 292 ε^0.31 MPa, rate-independent and isothermal, because a burst test takes
minutes.

**Nothing in the model breaks.** There is no failure stress, no critical strain, no criterion
of any kind — the charter forbids binary failure and this case does not need one. A pressurised
tube is *geometrically self-weakening*: as it expands the wall thins and the radius grows, and
both raise the hoop stress that the same pressure produces. Strain hardening pushes back.
Burst is the crossing point, where the material stops being able to keep up with its own
geometry:

```
  P(ε̄) = (2/√3)(t₀/r₀) σ̄(ε̄) exp(−√3 ε̄)       maximised where   dσ̄/dε̄ = √3 σ̄
```

For a pure power law that collapses to **ε̄\* = n/√3** exactly — the burst strain *is* the
hardening exponent, scaled by the biaxiality. The same material in a tension test necks at
ε = n (Considère), so a tube goes unstable at 58 % of the strain it would reach in uniaxial
tension, because biaxial tension weakens the geometry twice as fast. `Burst` reproduces
n/√3 to nine digits, which is what certifies the bisection behind the general case.

| | Burst pressure | vs closed form |
| --- | --- | --- |
| closed form, rigid-plastic thin wall | 18.1186 MPa | — |
| measured, ceiling at 1.03× burst | 18.0491 MPa | −0.384 % |
| measured, ceiling at 1.12× burst | 17.9942 MPa | −0.688 % |
| **measured, ceiling extrapolated out** | **18.0674 MPa** | **−0.283 %** |

#### You cannot walk up to a load maximum under load control

Past burst there is no equilibrium at all: every state the tube can reach carries less than
what is applied, so it accelerates and keeps accelerating. Reporting the pressure at which a
run visibly diverges would work, but it is biased high and the bias depends on the ramp rate
and on the divergence detector — three knobs, one answer.

So the pressure is not what gets measured. **The wall's capacity is**, read straight off the
stress field through the exact axisymmetric equilibrium identity

```
  d(r σ_r)/dr = σ_θ      ⟹      P·a = ∫σ_θ dr
```

which is a function of state, and so is defined on *both* sides of the peak. No approximation
anywhere — not thin-wall, not elastic, not small-strain. It holds in the current configuration
with Cauchy stress, which is exactly what the solver stores.

Two things make the traverse gentle enough to resolve. The ceiling sits a few per cent above
burst, so the imbalance driving the runaway starts near zero; and **the relaxation damping that
exists to reach static answers turns the runaway into a creep**, doing the job of a
displacement-control device with none of the machinery. The drag does not even have to be
modelled out: it enters the radial balance as a body force, so `P·a = ∫σ_θ dr + ∫ρcv·r dr` and
the hoop resultant is the applied pressure *minus* the drag — which is to say, the pressure the
wall would hold if it were standing still.

The ceiling is the one knob, and it is **extrapolated away rather than chosen**. The
perturbation is linear in it with no constant term, so a line through two ceilings gives the
zero-ceiling intercept; the third measurement is withheld from the fit as a witness and lands
on the line to **0.00002 %**, against a perturbation being removed 10,000× larger. Quadrupling
the ramp length leaves the answer unchanged in four decimal places.

#### The 0.28 % is the elasticity the oracle leaves out, and it takes two terms

Write the rigid-plastic curve as `P_rp(x) = C σ(κx) exp(−2x)` in total hoop strain `x`. Two
things change when the material is allowed to be elastic, and **they pull in opposite
directions**:

| | |
| --- | --- |
| elastic hoop strain thins the wall without hardening it | −0.385 % |
| elastic dilatation resists that thinning | +0.105 % |
| **net, from E and ν alone, before the run** | **−0.279 %** |
| measured | −0.283 % |

The first is the obvious one: flow stress is set by the plastic strain `x − e` while the
geometry thins by the total, worth `−2e`. The second is easy to miss. The `exp(−2x)` factor
*comes from* incompressibility, and a real material under a mean stress of +p swells by
`D = (1−2ν)3p/E`, which with ε_z = 0 goes entirely into ε_r + ε_θ — so the wall is thicker
than incompressible kinematics predicts. Together `P(x) = exp(D − 2e)·P_rp(x − e)`: the
rigid-plastic curve **translated by e and scaled by exp(D − 2e)**. Keeping only the first term
would over-predict the deficit by a third and leave a residual with nowhere to live.

Honest caveat on that 0.004 %: the thickness remainder crosses zero near D/t = 25, which is
where the reference tube happens to sit. At D/t = 80 the same two terms leave −0.034 %, and
that is the fair figure for how well they explain the deficit.

#### Freeze the geometry and there is no burst at all

The same tube in **small strain**, held at 1.5× its burst pressure, finds no load maximum
whatever. It settles — KE/SE = 6×10⁻⁶ — carrying half again what it can actually hold, having
expanded by 24 % of its radius while insisting nothing changed shape.

This is the first case in the project with no small-strain answer to fall back on, which makes
it the sharpest available test that the finite-strain path is doing work rather than
reproducing the small-strain one to three digits. Burst is a property of the *kinematics*.

#### A coarse mesh that gets the right answer for the wrong reason

Full 2×2 integration at 4 elements through the wall reports **−0.032 %** — better than the
converged answer's −0.283 %. It is not better. Volumetric locking under near-incompressible
plastic flow stiffens the wall and delays the instability, and the two errors very nearly
cancel *in the pressure*.

The peak strain does not cancel. It implies a 0.612 % elastic shift where 0.385 % is predicted;
the sound mesh reads 0.364 %. One number hides the locking completely and the other gives it
away, which is the whole argument for measuring a peak twice.

#### Other things measured here

- **Mesh: done at four elements**, second order, 0.004 % from coarsest to finest with a 4 ppm
  GCI band. The contrast with the Taylor mushroom is the point — that one is read at a singular
  corner and crawls at first order, while this is an integral of a smooth field over the whole
  wall and has no corner to be spoiled by.
- **The equilibrium identity holds to 2×10⁻⁶ %** on a settled sub-burst run, with KE/SE at
  5×10⁻¹⁶. That is what licenses trusting the capacity integral past the peak, where there is
  nothing left to compare it against.
- **Thin-wall limit**: the remainder is +0.21 % at D/t = 10 and −0.03 % at D/t = 80. Two tenths
  of a per cent at a wall far thicker than any thin-wall rule of thumb would allow, because
  carrying the radial stress through the equilibrium integral puts the **mean** radius in the
  formula, and that absorbs the entire first-order thickness term. **So the claim previously
  made here — that Lamé has to take over below D/t ≈ 20 — was wrong about what fails.** What
  fails there is the thin-wall stress *distribution*, not the burst pressure.
- **Barlow's formula** on the engineering UTS (218.9 MPa, from Considère on the same curve):
  17.509 MPa, −3.36 %, conservative. The reason is worth knowing rather than trusting:
  biaxiality raises the pressure a given flow stress can hold by 15 %, and also drops the
  instability strain to 1/√3 of the tensile value, so most of the first error is cancelled by
  the second. Two large mistakes nearly annihilating is why the rule of thumb survives — and
  why this gate is not scored against it.
- **A non-hardening tube** is unstable from first yield, and `Burst` then returns
  (2/√3)·A·t/r, which agrees with `ElasticPlastic.limitPressure` to (t/r)²/12 — third order,
  from two derivations sharing nothing. The gap quarters for every halving of the wall, which
  is asserted rather than admired.

### Gate 2 — the cost model

The charter's worked example: 50 mm bore, 5 mm wall, 500 mm long, 0.5 mm mesh,
10 ms of physics. 10 × 1000 = 10,000 elements, 42.7 ns timestep (CFL, safety 0.5,
**no mass scaling**), 233,998 steps = 2.34 × 10⁹ element-steps.

| Threads | element-steps/s | Speedup | Full run |
| --- | --- | --- | --- |
| 1, full 2×2 | 1.5 × 10⁷ | — | 2.7 min |
| 1, reduced | 4.0 × 10⁷ | 1.00× | 58.2 s |
| 2 | 7.1 × 10⁷ | 1.76× | 33.1 s |
| 4 | 1.15 × 10⁸ | 2.85× | 20.4 s |
| **8** | **1.48 × 10⁸** | **3.68×** | **15.8 s** |
| 8 independent runs on 8 threads | 1.81 × 10⁸ | 4.50× | 12.9 s each |

*(One run of the report; these move a few per cent between runs, and the ten-thread rows move
more than that because ten spinning threads on ten cores leave nothing for anything else.)*

**Seconds, not minutes.** That is the claim M0 exists to test, and it holds. Everything else
in this section is a correction to how it used to be reported.

#### The old projection was wrong twice

This gate used to print *serial rate × core count* and conclude 4.1 s. Both halves of that
were wrong.

**The measured ceiling is 4.5–5.0×, not 10×.** The last row above is the hardware speaking:
N independent runs on N threads, which synchronise never and balance perfectly. No threading
of one run can beat it, and it is half the core count. Two reasons, neither of them the code.
The ten cores are **four performance and six efficiency**, so an equal split runs at the speed
of the slowest; and the kernel is **memory-bound**, so single-thread throughput falls as the
mesh outgrows cache (4.06 × 10⁷ at 10,000 elements, 3.50 × 10⁷ at 300,000) while the threaded
rate plateaus at the same figure whatever the mesh size.

**And the target it was being compared against is a GPU number.** The plan's table is headed
*"Wall clock, one GPU"*, and the same row that budgets 2–10 s for this case budgets 11–58 min
for full 3D at the same element size. The old line had a CPU beating a GPU budget on its own,
which should have been the tell.

So gate 2 closes on CPU at 15.5–15.9 s — eight and ten threads land within a few per cent of
each other and which one wins varies run to run — and what the GPU has to deliver is now a
number rather than a hope: **1.6× to 7.8× over the CPU**, on a kernel whose limit is bandwidth.
That has not been attempted, and the structure-of-arrays layout that exists to make it cheap
has never been tested by doing it.

#### Threading without changing the answer by a bit

Threading recovers 82 % of that ceiling, and the standard it is held to is **bit-for-bit
identity at any thread count** — asserted over 400–600 steps of a plastic cylinder and of a
Taylor impact, at 1, 2, 3, 5 and 8 threads, and across a thread count changed mid-run.

A tolerance would not do, and the reason is specific. The two obvious ways to thread an
explicit kernel — an atomic add into the nodal force, or a private force vector per thread
reduced afterwards — both change the *order* a node's contributions are summed in.
Floating-point addition is not associative, so both give an answer that depends on the thread
count; the difference feeds a nonlinear return map and grows over 200,000 steps. A test with a
tolerance would pass at every thread count while the answers drifted apart, and the tolerance
would have to be widened as runs got longer, which is indistinguishable from a physics bug.

So the kernel **does not scatter at all.** Every source — element or pressure edge — writes
its own slot, and a second pass has each node sum the sources incident on it in ascending
source order, which is the order the serial scatter used. Identical arithmetic, not equivalent
arithmetic. It costs about 10 % of the serial rate and buys three things:

- **The schedule becomes free.** Since no chunk produces a sum, the result is independent of
  the partition entirely — equal chunks, unequal chunks, one thread or twelve, same bits. That
  is what licenses handing chunks out from a shared cursor rather than dividing them in
  advance, and on a machine with two kinds of core that is not a tuning detail:

  | Threads | Equal split | Shared cursor |
  | --- | --- | --- |
  | 2 | **0.87×** | 2.03× |
  | 8 | 1.79× | 3.74× |
  | 10 | 2.20× | 3.94× |

  The first row is the one to keep. **Two threads with an equal split are slower than one
  thread**, because the second half of the work lands on a core that cannot finish half of it
  in the time the first core takes to do all of it, and the barrier waits. None of that is
  visible in a profile of the kernel — it is entirely a property of the schedule.
- **The pressure load rides along.** Sharing the buffer puts it in the same parallel region
  as the elements. Kept separate, a thousand edges are too little work to be worth a barrier
  and too much to leave on one thread — a third of the threaded step, measured.
- **It is the layout a GPU kernel wants anyway**, which is the port this gate is still open on.

What is *not* threaded is the nodal update when a rigid wall is present. The wall carries
running totals — impulse, destroyed kinetic energy, frictional dissipation — and splitting
those across a dynamic schedule would make the audits vary between runs. Those audits caught
the one real bug in the contact work; they are not worth trading for speed on a case whose
mesh is two orders of magnitude smaller than this one.

#### What is left on the table

The remaining 18 % is the barrier: three per step, and 10,000 elements is too little work to
hide them. A wider mesh threads visibly better — 3.6× at 10,000 elements against 4.3× at
100,000, same eight threads. That is an argument **for** the plan's batched ensemble rather
than against threading: a sweep synchronises once per run instead of three times per step, it
is already the fastest thing measured here, and it is what the product is for. The same
argument gets stronger on a GPU, where a per-step dispatch at this mesh size costs more than a
per-step barrier does here.

Reduced integration itself delivers **2.6–2.8×, not the 4× it was budgeted at**. The stress
kernel drops to a quarter of the work, but the per-element gather and scatter do not, and
hourglass control adds some back — Amdahl applies inside the element too.

### M2 — the defect field

Real things do not fail where the continuum stress is highest. They fail at defects. The plan's
construction for this is one spatially correlated random field on local failure strain, one
float per element, and it is meant to make several separate-looking behaviours into
consequences of the same number: failure lands somewhere plausible, bigger parts are weaker,
nominally identical parts burst at different pressures, and a weld fails first because its
elements draw from a worse distribution rather than because a rule said welds fail.

The exit criterion is that the size effect reproduces `(V₂/V₁)^(1/m)` and that repeated runs
give a distribution. Both hold. Measured on four tubes of identical wall and length with the
radius doubling each time, 3000 parts each:

| V/V₁ | measured strength ratio | (V/V₁)^(−1/m) | error |
|---|---|---|---|
| 1.889 | 0.96983 | 0.96870 | +0.12 % |
| 3.667 | 0.93598 | 0.93710 | −0.12 % |
| 7.222 | 0.90444 | 0.90587 | −0.16 % |

Nothing in the field knows the radius.

#### The one identity everything here rests on

The minimum of *n* independent Weibull samples is **again Weibull**, same modulus, scale
multiplied by `n^(−1/m)`. That is exact, and it collapses three things that look like three
features into one calculation:

- **Bigger parts are weaker.** A body of volume V contains V/V₀ reference volumes, so it is the
  minimum of that many. Eight times the volume at m = 10 is 18.8 % weaker — the plan says
  "about 19 %", and the figure in the test is the computed one.
- **Volume normalisation.** An element is given the distribution of the minimum over *its*
  volume, not a fresh draw from the base law.
- **An axisymmetric element is a ring.** A real flaw is a point; an element here is a full 360°
  torus, so its one stored number has to stand for the weakest of every flaw around the
  circumference. Normalising by the true ring volume 2πrA *is* that statement. There is no
  second correction, and the r^(−1/m) prediction this README flagged as a hypothesis last
  milestone turns out to be a restatement of volume normalisation rather than an addition to
  it.

Because those are one calculation, `Weibull.forVolume` is implemented as a call to
`Weibull.minimumOf` rather than as the same formula written twice. A test asserting they agree
would otherwise compare two roundings of the same algebra, and pass long after the intent had
drifted.

#### The mesh dependence the plan names, and the one it does not

The plan's trap: a naive per-element draw gets weaker without bound under refinement — more
samples, worse extreme. Its two fixes are a correlation length set by a physical scale and
volume normalisation per Weibull theory. Both are implemented and both are necessary.

**They are not sufficient, and what is left over drifts the other way.** Volume normalisation
says an element of volume V had V/V₀ independent chances to be bad. That is true only while the
chances really are independent. Refine past the correlation length and neighbouring elements
read almost the same value of the field — one chance shared out, not many — while each is
separately credited with the strengthening its own small volume earns. Measured, before the
fix, on the same part and the same seed:

```
elements   ring volume only   with the resolution floor
      32        0.27466              0.27466
     512        0.30095              0.27736
    8192        0.34525              0.27699
   32768        0.36991              0.27690
```

35 % of drift over a thousandfold refinement, and every individual number in that column looks
reasonable. It drifts **upward**, which is why it survives any test phrased as "the part must
not get weaker" — and upward is the flattering direction, so it is the one more likely to ship.

The fix is a floor on the count of independent chances at what the field can actually resolve:

```
V_eff = ℓ³ · max(1, A/ℓ²) · max(1, 2πr/ℓ)
```

The two factors are floored for different reasons, and the difference is the whole of the
axisymmetric story. **In plane**, the field is a real function of position: an element bigger
than a correlation cell covers A/ℓ² of them while being read once, so it is credited with that
many — the standard local-averaging correction — and an element smaller than a cell covers one.
**Around the circumference**, the field is not a function of position at all, because there is
no third coordinate; those 2πr/ℓ cells are never resolved at any mesh density, so the element is
always credited with all of them. That is where the ring effect comes from and why it does not
wash out under refinement. Above the correlation length the expression is exactly the ring
volume and nothing changes.

The floors are a `max` rather than a smooth blend. Vanmarcke's variance function would be the
usual choice and would be defensible; it would also put a fitted shape in the middle of a chain
that is otherwise exact, to smooth a corner each element crosses once.

#### Two errors that are each large and cancel

The strongest available statement is not that the field is stable but that it is *right*: run
many nominally identical parts, take the weakest element of each, and the population is the
Weibull that weakest-link theory predicts for a body of that volume — mean and scatter both.
Across a 144-fold change in element count the mean runs 1.0009, 1.0026, 0.9917, 0.9720, 0.9677
of the predicted value: a 3.2 % residual at the finest mesh, still settling, and reported as
that rather than as agreement.

That is not trivially true. On a fine mesh the per-element effective volumes sum to **far more**
than the body volume, because each element is credited with a full correlation cell it only
partly occupies. In the other direction the elements are strongly correlated, so far fewer of
them are independent draws than there are elements. Both errors are large. They cancel, and the
cancellation is the design rather than a coincidence — the over-credited volume per element is
exactly the volume per *independent* element.

#### Construction

Gaussian copula. White noise on a grid at the physical scale, smoothed by a separable Gaussian
kernel into a field with correlation `exp(−d²/2ℓ²)`, normalised to **exactly** unit variance at
every point by dividing by the root sum of squared weights, then pushed through the normal CDF
to a uniform and through the Weibull quantile to the marginal. The exact normalisation is what
makes the marginal exactly the Weibull that was asked for, which is what lets the size effect
be tested against closed form rather than against itself.

Two details that are not decoration:

**The white noise comes from a stateless hash of the seed and the integer grid coordinates**,
not from a sequential generator. A sequential generator would make a node's value depend on how
many nodes had been drawn first, so the field would move when the bounding box changed, when
the traversal order changed, or when the work was threaded — reintroducing exactly the mesh
dependence this object exists to remove, in a form far harder to see. As it is, two meshes of
the same body see the same field at the same physical point **bit for bit**, and that is
asserted with exact equality rather than a tolerance.

**The normal CDF is Cody's, not Abramowitz–Stegun 7.1.26.** The usual approximation is good to
1.5e-7, which looks perfect in any plot and is a million times worse than double precision. It
sits in the middle of the copula, so its error lands directly on the Weibull marginal and comes
out as a size effect that misses by a per cent — indistinguishable from the construction being
wrong, and absorbed by relaxing the tolerance on the test that is supposed to be proving it.
Measured against the platform's own `erfc` at 7681 points across ±30, this one is within three
units in the last place.

One consequence worth knowing before it is mistaken for a bug: the far tail is steep enough
that forming `x/√2` by multiplying by the rounded reciprocal, rather than dividing, moves the
answer by 1e-13 relative at 30σ. The logarithmic derivative of erfc is about −2x. That is
argument sensitivity, not approximation error, and there is a test that says so.

### M2's second half — what the solver does with it

The field is one number per element; this is the part that lets it matter. Damage enters as
**progressive softening** — the charter forbids binary failure, so nothing is deleted and
nothing switches off. Once accumulated plastic strain reaches an element's failure strain, the
flow stress is frozen at whatever it had reached and then falls linearly to zero.

**Damage is negative hardening, and that is the whole implementation.** There is no damage
variable in the state and no second constitutive branch: the frozen flow stress and the
softening slope take the places that `sy0 + hIso·εp` and `hIso` occupy for an undamaged point,
and every line of the return map after that is unchanged. Freezing buys three things at once —
the softening branch becomes exactly triangular, so its energy has a closed form rather than an
accumulator; the return stays closed-form **even for Johnson–Cook**, because a frozen flow
stress removes the nonlinearity that the bracketed Newton solve exists for; and the fracture
energy stops depending on strain rate and temperature, which is the point, since G_f is meant
to be a material constant.

#### The band width is not the square root of the area

Softening localises — that is what a negative tangent does — and in an explicit solve the band
is exactly one element wide whatever the mesh. So a fixed softening slope would make the energy
to break the part go to zero as the mesh refined. Hillerborg's fix is to fix the energy per unit
**crack area** and let the slope follow the element: `H_soft = −σ_f²·h / 2G_f`.

That makes `h` the whole question, and the usual shortcut of `sqrt(A)` is the element's band
width only if it is square. A tube wall is meshed thin in r and long in z; 10:1 is ordinary, and
`sqrt(A)` is out by the square root of that — a factor of three, in the flattering direction,
indistinguishable from the regularisation working. So the orientation comes from the physics
instead: at onset the plastic flow direction's largest in-plane principal direction is where
the material is stretching hardest and therefore where it will open, and with `t` along the band
and `u, v` the element's mean edge vectors,

```
h = A / (|u·t| + |v·t|)
```

exact for any parallelogram. Across a rectangle it returns the other side; diagonally across a
square it returns `a/√2`, not the support width `a√2`, which is the mistake that looks right
until the band is at 45°.

An element here is a ring, so it is worth checking this is still a length. Ring volume is
`2πrA` and a circumferential crack through it has area `2πrL`, so `G_f = g_f·A/L` and **the
radius cancels on both sides**. A ring at 100 mm and a ring at 1 mm with the same cross-section
cost the same per unit area to break, which is what a material constant has to mean.

One case is genuinely undetermined and it is worth knowing why. A thin open-ended tube has all
its load in the hoop, so the flow direction's in-plane part is isotropic and there is no
preferred meridian direction at all. That is the model telling the truth: the crack that wants
to form is a longitudinal split, and an axisymmetric formulation has none. Below 1e-9 of
anisotropy it falls back to `sqrt(A)`.

Under full integration all four Gauss points are handed the **element's** geometry, not a
quarter each. Giving each point its own quarter-sized band is the obvious thing and it is wrong:
a bilinear quad cannot represent a strain discontinuity in its interior, so the narrowest band
the mesh can resolve is one element across however many points sample it. Sizing off the point
would let one element host two parallel cracks in each direction and dissipate twice — a
quadrature rule changing the fracture toughness. Full and reduced integration agree on the
energy to 1e-10, and that is the test.

#### The size limit, clamped and counted

The element's total response — elastic unloading plus softening — has to stay monotone, or the
branch snaps back and the element cannot be driven through it at all: `h ≤ 2E·G_f/σ_f²`. It
cannot be checked at construction, because σ_f includes however much the material hardened on
the way, so it is checked at onset, clamped at `H_soft = −E`, and **counted**. A clamped point
dissipates more than G_f and its mesh independence is gone; reporting the count is the
difference between knowing that and not. The clamp is also what keeps the return map
non-singular — the worst denominator is `2μ − (2/3)E`, positive for every ν < 0.5 and vanishing
exactly at 0.5, which `Material` already refuses.

#### What holds exactly

A through-wall ring crack costs `G_f × π(r_o² − r_i²)` — to **1e-10**, at 2, 4 and 8 elements
through the wall, and across an eightfold axial refinement. The 2πr cancels, the annulus area
does not depend on how many elements span it, and neither does the energy. The solver's own
running total undershoots the closed form by exactly one step in n, a right Riemann sum of a
straight line, and that is asserted as `1/n` rather than absorbed into a tolerance.

#### The result, and the two things that were not in the plan

Six nominally identical tubes, seed the only difference: **16.648 to 16.920 MPa, failing at six
different axial stations.** That is the plan's claim and it holds.

**A pressurised shell averages its defects instead of failing at its weakest link.** Weibull
weakest-link theory — all of `Weibull`, all of `DefectField` — predicts that a longer body is
weaker and that its scatter is size-independent, because the minimum of n Weibulls is Weibull
with the same modulus. Neither happens:

| Tube length | L/√(Rt) | of defect-free | scatter | failure location, sd/L |
| --- | --- | --- | --- | --- |
| 10 mm | 1.41 | 0.9354 | 0.933 % | 0.226 |
| 20 mm | 2.83 | 0.9280 | 0.742 % | 0.377 |
| 40 mm | 5.66 | 0.9293 | 0.526 % | 0.284 |
| 80 mm | 11.31 | 0.9298 | 0.408 % | 0.270 |

The mean holds still and the scatter falls as **L^(−0.40)** — against L^(−1/2) for an average
over independent patches and L^0 for a weakest link, which is much nearer the first. (The last
column is the spread of the failure location as a fraction of the tube; 0.289 is what a uniform
distribution over the whole length gives, so the field is choosing and the boundaries are not.)
The reason is the **shear-lag length**
`√(Rt)` = 7.07 mm — the distance over which a thin shell shares load along its axis. Every weak
patch here is 3 mm long, so the neighbours carry it. Pushing the correlation length up through
`√(Rt)` turns the behaviour over exactly as that argument predicts: on an 80 mm tube, at 2, 4,
8, 16 and 32 mm it holds 0.9373, 0.9262, 0.9140, 0.9019 and 0.8905 of the defect-free pressure,
monotonically weaker, with the scatter climbing 0.46 %, 0.60 %, 0.88 %, 1.59 % and then 1.33 %
— that last figure is five draws from a tube containing two and a half patches and is noise,
not a turn. **Along its axis this tube is a bundle of parallel
rings, not a chain.** Through its wall it is a chain, which is why refining through the wall
changes nothing. The weakest-link apparatus is right about the material and describes only part
of the structure it is put into, and which part depends on a length scale that belongs to the
geometry.

#### The burst pressure does not converge

| Elements | Δz | Burst | of defect-free |
| --- | --- | --- | --- |
| 160 | 1.000 mm | 16.880 MPa | 0.9316 |
| 320 | 0.500 mm | 16.923 MPa | 0.9340 |
| 640 | 0.250 mm | 16.955 MPa | 0.9358 |
| *no field* | — | 18.031 MPa | 0.9952 |

**+0.22 % per mesh doubling, monotone, no plateau, heading for the last line.** What evaporates
under refinement is the knockdown itself.

Crack-band regularisation fixes the energy per unit crack area — exactly, to 1e-10, and that
identity holds — but it fixes it by making the softening *modulus* depend on the element size,
and that is only the right thing to do once the band has localised into one element. Here the
load maximum arrives while damage is still diffuse, so a finer mesh softens more slowly at the
same strain and the tube comes out stronger. One doubling moves the answer by about 40 % of the
0.53 % scatter it is meant to be measuring: the ranking between two tubes survives refinement,
the absolute knockdown does not.

This is the second mesh dependence in M2 pointing the flattering way, and the mirror of the
first. Volume normalisation made the part stronger as the mesh refined because it credited
correlated elements with independent chances; this makes the part stronger as the mesh refines
because it credits a narrower band with a longer softening range. The fix is a softening modulus
that is a material constant until localisation is detected and regularised only after it, or a
nonlocal damage model. Neither is in this increment, and the number is printed rather than
smoothed.

#### What this damage model does not have

**Triaxiality.** The failure strain is a property of the material and the defect field and
nothing else. Real ductile failure strain falls roughly exponentially with stress triaxiality,
and a pressurised wall spans a wide range of it between bore and outside, so the bore should
fail earlier than this makes it. Until that is in, *where* a part fails is set by the defect
field alone where it ought to be set by the defect field and the stress state together.

**Loss of hydrostatic strength.** Only the flow stress is degraded, so a fully damaged element
has no shear strength but still resists volume change: it is a fluid, not a crack. For a
bursting wall that is enough — the wall fails by losing hoop capacity and thinning, and hoop
tension is deviatoric. It also keeps plastic flow exactly isochoric, which is what the plastic
dissipation audit rests on. For anything that needs a traction-free surface to open, it is not
enough.

**And for a real ductile metal this whole mechanism is a footnote.** Annealed copper's ductility
is several times its instability strain, so the instability wins outright and the measured
effect of a defect field is a 0.6 % knockdown with 0.06 % scatter. The runs above use a failure
strain set below the instability on purpose, to have something to show. That is a statement
about what is being demonstrated, not about copper.

#### The axisymmetric price

One thing this cannot represent, and it is the thing a real tube actually does. A burst tube
splits **longitudinally**. Axisymmetry admits only circumferential cracks, so what this models
is the wall losing its hoop capacity and bulging, not the split that follows. The localisation
is a bulge and not a neck for a second reason as well: with ε_z = 0 the wall can expand unevenly
along z, because neighbouring rings may move by different amounts, but it cannot draw material
along the tube into the neck.

### The first sweep — Barlow, measured

The M3 exit criterion in the build plan is that *a player rediscovers Barlow from a sweep
without being told*. `./gradlew barlowSweep` is the machine half of that, and it costs 24 runs
and 11 seconds on eight cores. Nothing in it consults `Burst.barlow`; the exponents come off
the finite element solver, and the rule of thumb is compared against them afterwards.

| Sweep | Fit | Worst residual |
| --- | --- | --- |
| wall thickness, D/t from 100 to 10, mean radius fixed | `P ∝ t^0.9989` | 0.020 % |
| bore diameter, factor of 8, wall fixed at 2 mm | `P ∝ d^-0.9666` | 1.279 % |
| **the same runs against the mean diameter** | `P ∝ D^-0.9989` | 0.019 % |

Three things came out of it, and the second was not expected:

- **The exponents are the shape of the law; the coefficient is the material in it.** The
  thickness sweep ran at a fixed mean diameter, so solving its fitted coefficient for σ in
  `P = 2σt/D` gives **225.54 MPa** against the flow curve's own effective stress of 226.48 MPa
  — 0.42 % low, which is exactly the single-run creep bias every point carries. The solver
  recovered how strong the copper is from twelve tubes of different thicknesses bursting.
- **The bore-diameter sweep is wrong in a way that teaches something.** With the wall held
  fixed, burst pressure goes as `1/(d + t)`, which is not a power law in `d` at all. The fit
  returns −0.9666 and hides the rest in a residual, and the residual plot is a clean arch
  rather than noise. That arch *is* the distinction between Barlow written on the bore, the
  mean and the outside diameter, and it is most of why the rule of thumb is conservative.
  A fit reported without its residual would have called this a success.
- **`r²` is useless here and is reported anyway because people expect it.** It is 0.999907 on
  the bent fit — indistinguishable, at a glance, from the 1.000000 on the clean one. Over a
  decade in `x`, a power law fit has an `r²` of three nines almost regardless of how badly it
  is wrong in the middle. `PowerLaw.maxResidualPercent` is the number to read first, and the
  tests assert both halves of that: that the residual notices the structure and that `r²`
  does not.

Both thickness and diameter exponents land at 0.9989 in magnitude, and the residual arch on
the *thickness* sweep is the same shape at ±0.02 %. That is the second-order thick-wall term
that `Burst`'s mean-radius formulation was chosen to leave behind, showing up as the only
thing left once the first-order error is gone.

### The design map

`./gradlew designMap` fires nine wall thicknesses, twelve nominally identical tubes each, and
colours every cell of a wall-thickness-against-applied-pressure grid by how many of them gave
way. 108 runs, 29 s on eight cores, 2.1 s of solver per tube.

**The pressure axis is free, and that is worth saying out loud.** A map of *n* thicknesses by
*m* pressures looks like `n·m` runs and is not. A tube's burst pressure is a property of the
tube — the maximum of the wall's own load capacity, read off an equilibrium identity, measured
without reference to what ceiling the harness used to find it. So a tube bursts under applied
pressure `P` exactly when `P` exceeds its capacity, one population per thickness answers every
pressure at once, and the pressure resolution is limited by nothing but how many rows fit on a
terminal. The obvious implementation is two orders of magnitude more expensive and produces
the same picture.

| wall, mm | D/t | mean burst | defect-free | knockdown | scatter |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1.250 | 40.0 | 10.542 MPa | 11.324 MPa | −6.90 % | 0.81 % |
| 2.236 | 22.4 | 18.801 MPa | 20.257 MPa | −7.19 % | 0.77 % |
| 4.000 | 12.5 | 33.430 MPa | 36.237 MPa | −7.75 % | 0.72 % |

- **The band is sub-pixel on the axes anyone would choose first.** The widest column spans
  0.800 MPa from its weakest tube to its strongest, against a row height of 1.213 MPa on a
  map covering the whole design range. On those axes the boundary is a clean staircase and it
  would be easy — and wrong — to read that as the band not existing. Divide the pressure axis
  through by each column's own defect-free burst pressure and it opens up across seventeen
  rows. The report draws both and says which is which.
- **The band slopes, and the slope is the interesting number.** Knockdown runs from 6.90 % at
  1.25 mm of wall to 7.75 % at 4.00 mm — a size effect over a factor of 3.2 in material.
  Weakest-link theory says `V^(−1/m)` over that range, which at `m = 12` is a **9.24 %** fall.
  The measured fall is **0.91 %**, smaller by a factor of ten. The size effect is real and it
  is nothing like as strong as a chain of links would make it, which is exactly the M2 finding
  — a pressurised shell averages its defects over the shear-lag length — reproduced from an
  experiment that was not designed to test it.
- **Knockdown and scatter are both fractions, not pressures.** That is what makes the
  normalised band nearly horizontal, and it is what makes a derating rule a rule rather than a
  table. It is also a safety factor in the literal sense, derived from tubes rather than looked
  up.

`DefectBurst.Setup` gained a `slenderness` so the wall could be swept at all; it previously
hardcoded the reference tube's. Holding the element *count* through the wall rather than the
element *size* is deliberate — the through-wall discretisation error is then identical at every
column, so what moves between columns is physics. It does mean a thin wall costs more, because
the radial element shrinks with the wall and the CFL timestep follows it.

### What the curve was hiding

`./gradlew burstCurve` traces the reference tube's `P(ε_θ)` through the maximum instead of
reporting one point off the top of it. The peak was never the interesting part:

- **Nothing breaks and no criterion is consulted.** The curve rises while strain hardening
  wins, flattens as the wall thins and the radius grows, and turns over where geometry starts
  winning. The turnover is a 0.5 % feature on a 19 MPa curve, which is why the report draws it
  twice — once on axes that show the elastic rise and the yield knee, and once with the
  pressure axis opened up, where it is the only thing on the page.
- **The gap between the two curves past the peak is the instability, drawn.** Applied bore
  pressure sits on its ceiling; the wall's capacity falls away beneath it; the difference is
  the unbalanced force accelerating the wall outward. Up the stable branch they track each
  other to within the relaxation drag, and that agreement is what licenses reading the
  capacity past the peak at all.

## Layout

```
core/     Formulation (the axisym / plane-strain / plane-stress flag), Material
          Materials — which material each element is; the uniform case is one
              reference, so a single-substance mesh is unchanged bit for bit
          JohnsonCook — flow stress: strain, rate and thermal terms
          Weibull — the weakest-link algebra; minimumOf(n) is the whole of it,
              and forVolume is a call to it rather than the same formula twice
          Normal — Cody erfc, because the copula's error lands on the marginal
          DefectField — correlated failure strain from a stateless hash;
              effectiveVolume is where the ring and the resolution floor live
mesh/     QuadMesh — flat primitive arrays, the layout that ports to a GPU kernel
              cylinderWall for the annulus, solidCylinder for the Taylor specimen
              ringVolume — Pappus, 2*pi*rc*A, the volume a defect draw normalises by
          Outline — any polygon, holes included, rasterised onto a square lattice
              because a fitted mesher makes slivers and slivers set the timestep
              rotated/translated/restingOn, and one lattice for many shapes
              mesh welds shapes that touch into one body; assemble keeps them
              separate, which is the difference between a clad plate and a stack
solver/   ExplicitSolver — the hot loop
          Parallel — spin-barrier worker pool; chunks from a shared cursor,
              which is safe only because nothing accumulates per chunk
          Integration — full 2x2 or reduced + hourglass
          Kinematics — small strain or updated Lagrangian
          J2 — radial-return plasticity; closed form for linear hardening,
              bracketed Newton for Johnson-Cook. J2.Flow bundles the uniforms
          Damage — progressive softening as a NEGATIVE hardening modulus, so the
              return map needs no second branch; the band width comes from the
              flow direction, and the snap-back clamp is counted rather than hidden
          Corotational — objective incremental rotation
          RigidWall — kinematic anvil with a Coulomb cone; stick/slip is a
              return map, same shape as J2
          Contact — penalty between deformable surfaces. Stiffness from the pair
              mass and the reference step, so it is material- and formulation-
              agnostic; three tests decide a pair, and the one that is easy to
              miss is that the two surfaces have to FACE each other
          ExplicitSolver also carries gravity, erosion and mass scaling, all
              off by default so that every validation gate runs the arithmetic
              it always did
validate/ Lamé and ElasticPlastic closed forms, the cylinder cases,
              TaylorImpactCase — the M1 gate, which has no closed form,
              Burst — the instability condition, the P(ε) curve, Considère,
                  Barlow, and the elastic correction to a rigid-plastic oracle
              BurstCase — traverses the load maximum and measures the wall's
                  capacity rather than the pressure applied to it
              DefectBurst — the same tube many times over; Setup varies one named
                  thing, volley pools a population, shearLagLength is the number
                  the correlation length has to be compared against
              MeshConvergence — Richardson extrapolation and the GCI band
sweep/    Sweep — one parameter varied, one run per thread, results in parameter
              order; linear and logarithmic spacings with exact endpoints
          PowerLaw — the log-log fit. maxResidualPercent is the statistic that
              matters, because r² stays above 0.999 on data the exponent is
              visibly wrong about
report/   Plot — a line chart in characters, for keeping a curve visible before
              any of the product's plotting exists
          Heatmap — the two-parameter map. Blank for zero, because on a map
              whose point is where the outcome changes, the quiet region
              should carry no ink and the boundary should be the thing with it
          Csv — the same curve at seventeen digits, for plotting it properly
          RunLog — one appended line per run: what went in, what came out. Text,
              because a log nobody can read is not a log; parameters are the
              identity, so lastMatching answers "have I run this already"
          Runs — where traces land. Not out/, which the IDE already claims
render/   Film — deformed geometry over a run, plus per-element fields and
              readouts. boundaryEdges finds the free surface from connectivity
              alone, so it keeps working for shapes QuadMesh cannot yet make.
              Film.AUTO fits a colour ramp to the 99.5th percentile and reports
              the peak beside it, so one wild element cannot flatten the picture
viewer/   burst-viewer.html — the tube. The scrub track is the wall capacity
              curve, so the burst is a point on the control itself
          scene-viewer.html — any planar film: deformed mesh, four fields, the
              anvil, and the bodies at rest as an outline behind them. On a scene
              that erodes it rebuilds the outline per frame from the live set,
              because the free surface moves as elements go. Every film
              file appends itself to one global array, so a page that includes
              several of them gets a scene selector and none of them knows it
          serve.js — pages from viewer/, film data from runs/
```

### Knowing how much to believe a number

`MeshConvergence` takes three systematically refined meshes and returns the **observed order**,
the Richardson-extrapolated value, and a Grid Convergence Index band on the finest result. It
is what turns "the mushroom hasn't converged" into "10.134 mm ± 1.0 %", and the observed order
is the more useful of the two outputs: a quantity converging at the rate the element promises
is merely under-resolved, while one converging more slowly is telling you something about the
*problem*.

It refuses more than it computes, because the arithmetic will produce an authoritative-looking
number from any three values at all. Oscillatory sequences have no order to observe.
Sequences whose steps are growing extrapolate the wrong way.

And one subtler trap, which the final length walked straight into: **a quantity that has
already converged.** L/L₀ moves one part in ten thousand per refinement — residual dynamics,
not a mesh trend — but the steps *are* very slightly shrinking, so a shrinking-steps test lets
it through, the denominator `r^p − 1` goes near zero, and Richardson confidently relocates it
twelve times further than the entire refinement study moved it. The guard is that the proposed
correction must not exceed the movement refinement actually produced. A converged quantity is
not a failure and shouldn't be dressed as one: it gets reported as its finest value and the
size of its last step.

The solid-cylinder mesh puts a column of nodes **on the axis**, at r = 0, which looks fatal
for a formulation that divides by radius to form the hoop strain. It is not: that quotient is
only ever evaluated at a quadrature point, and every quadrature point of an element with
non-zero radial extent lies strictly inside it. What the axis does require is that its radial
degree of freedom be pinned — u_r(0, z) = 0 is a symmetry condition, not a boundary condition
the caller may choose, and forgetting it does not crash. It quietly lets the axis open into a
hole down the middle of the specimen while every other number stays plausible.

The small-strain path is kept for the reason Abaqus keeps NLGEOM off, and for one more: the
small-strain gates have exact closed forms behind them, so running both proves the
finite-strain machinery agrees where it must before it is trusted where no closed form exists.
All 38 pre-existing tests passed unchanged after the finite-strain work landed; all 45 did
again after contact and plastic-work accumulation; all 69 after Johnson–Cook; all 93 after
Coulomb friction, where a test asserts that a zero coefficient reproduces the previous
frictionless path bit for bit; and all 134 after the assembly was turned inside out for
threading. Each time the only edits to existing tests were mechanical — threading a new array
or a parameter bundle through a helper signature — which is the evidence that the older paths
are really untouched rather than merely still passing. The threading change made no edits to
them at all, which is the strongest version of that evidence available: the hot loop was
restructured and 134 assertions about physics did not notice.

**The solver is incremental, not total.** M0 computed stress from total displacement — state
in, stress out, no memory. That stops working the moment the material can yield, because a
plastic stress depends on the path taken and not merely where the material ended up. Stress,
back stress and equivalent plastic strain are stored per Gauss point and integrated along the
strain path. This is the real architectural change plasticity brings, and it is another point
for reduced integration: one history point per element instead of four.

An elastic material is expressed as a plastic one with **infinite yield stress**. Not a trick
to save a field — it means the return map has exactly one code path, so the elastic gates
keep testing the plastic code.

Radial return is **exact, not iterative** for linear hardening: the consistency condition is
linear in the plastic multiplier, so the closed form lands precisely on the updated yield
surface. No Newton loop, no convergence tolerance, no failure mode where a badly conditioned
element quietly stops converging.

Johnson–Cook breaks that. `A + Bεⁿ` is nonlinear in plastic strain and the rate term depends
on the multiplier itself, so the multiplier appears on both sides and has to be solved for.
Only that step changes; everything after it is the same arithmetic, and the linear path keeps
its closed form rather than being routed through the solve.

The guarantee is **restored rather than abandoned**, which matters because the naive version
fails exactly where plasticity starts. Since n < 1, the hardening slope `Bnε^(n-1)` is
*infinite* at zero plastic strain, so an unguarded Newton step at a point's first yielding
increment is zero and the iteration sits still — at every point of every plastic zone, on the
first step it matters. The solve is therefore bracketed:

- `g(0) > 0` is exactly the condition the elastic check already tested, because the yield
  check uses the same quasi-static surface the residual starts from. The two cannot disagree
  about whether a point is flowing.
- `g(‖ξ‖/(2μ + ⅔H_kin)) < 0`, that being the multiplier a material with no strength at all
  would need.
- `g` is strictly decreasing, so the root between them is unique.

Every iteration either takes a Newton step that stays inside the bracket or halves it. There
is no non-convergence mode; a badly conditioned point costs more iterations, not a wrong
answer. A test asserts that with `B = C = β = 0` — Johnson–Cook degenerated to perfect
plasticity — the solve reproduces the closed form to machine precision, which is the only way
to know the iteration is landing where the algebra says it should.

It costs **5.2× the closed form** at the same mesh — 3.5 s against 0.67 s at n = 16, reduced —
because each yielding point runs a handful of iterations carrying transcendentals. It used to
cost 6–8×, and the difference is that the hardening slope is now taken from the flow stress
already computed rather than from a second `pow`: `Bnε^(n-1) = n(S − A)/ε` exactly. The
cancellation in that subtraction is real and harmless, because it only bites where `Bεⁿ` is
small against `A` — just past first yield, where the true slope is enormous, and where any
Newton step the bracket does not like is replaced by a bisection. The convergence guarantee
never rested on the slope being accurate, only on the bracket being real.

**Temperature is explicit and needs no storage.** It is a function of the plastic work already
accumulated per Gauss point, `T = T_room + βW_p/(ρc_p)` — no heat equation, because over the
tens of microseconds an impact lasts the thermal diffusion length in steel is well under one
element, so no point exchanges heat with any other. Taking it at the start of the step rather
than solving for it alongside the multiplier lets the stress end a step fractionally outside
the *end-of-step* surface, by the amount the point heated while flowing. That is measured
rather than asserted away: at the increment the reference case actually takes it is ~1e-4
relative, against a gate of 5e-2.

Both hardening kinds are carried because they are **indistinguishable under monotonic loading
and completely different under reversal** — the Bauschinger effect. Anything that cycles (a
barrel breathing shot after shot, a case wall flexing) is a reversal problem, and a purely
isotropic model over-predicts the reverse yield point every time.

Hourglass control is Flanagan–Belytschko, **stiffness form rather than viscous**. The viscous
form resists hourglass *velocity*, which is useless under dynamic relaxation: velocity is
driven to zero, at which point a viscous term stops resisting and whatever hourglass
displacement has accumulated is frozen in. The stiffness form resists displacement, holds at
equilibrium, and being conservative has an energy that can be audited. That audit
(`hourglassEnergy()` over `strainEnergy()`) runs at 1e-29 to 1e-26 on the cylinder and is
gated at 1e-4 — hourglass stiffness is a numerical prop, not physics, and the only honest
question is how much of the answer it is carrying.

The shape vector is the projected one, γ_i = h_i − (h·r)b_i^r − (h·z)b_i^z, not the raw
alternating pattern. On a distorted element the raw pattern is contaminated by genuine linear
deformation, so resisting it would resist real strain. The projection makes γ orthogonal to
every linear displacement field, so the correction vanishes identically on anything the
element is supposed to represent. That property is asserted **to machine precision** rather
than to a tolerance, which is what makes it a correction rather than a fudge.

Two more choices that look like details and are not:

- **Structure-of-arrays over flat primitives**, not an object graph. This is what ports to
  CUDA or HLSL without restructuring, and what makes batched ensembles cheap later: N lanes
  are N strides into the same buffers.
- **B is recomputed every step**, though M0's small-strain elastic physics would let it be
  cached for a large speedup. The solver this becomes recomputes B from the deformed
  configuration every step, and a timing number from a cached B would flatter the cost model
  that M0 exists to test.

And one that is a stated principle: **no mass scaling, ever.** It is the standard trick for
making explicit dynamics affordable and it corrupts exactly the inertial behaviour this
project is about. Offline solve is what buys the right to refuse it.

## Direction

**The target is a general 2D destruction sandbox**, physically accurate rather than exactly
real, in which a thing you cannot do should be a scope decision — *the game does not grow
trees* — and never a solver limitation. Everything above is the validation programme that
says the material model is right. It is the test rig, not the game, and from here it is
frozen as a regression suite rather than extended.

What that target refuses today, and none of it is the physics:

| Needed | Today |
| --- | --- |
| Many materials in one scene | **Done** — `Materials`, one `Material` per element |
| Arbitrary shapes | **Done** — `Outline`, any polygon with holes, rasterised to quads |
| Cross-sections, not lathe shapes | **Done** — `PLANE_STRAIN`, and `Impact` is built on it |
| Objects that touch each other | **Done** — `Contact`, penalty, node against segment |
| Things falling | **Done** — `setGravity`, a uniform body force, off by default |
| Things breaking into pieces | **Done** -- elements softened to nothing are deleted |
| Mixed-scale scenes | **Done** -- `setMassScaling`, exact on a lattice; see below |

What is left is speed, and the estimate this section used to carry was wrong by an order of
magnitude. `./gradlew bench` measures it instead, on two blocks meeting at 120 m/s:

| configuration | element-steps/s |
| --- | --- |
| bare elements, finite strain | 9.03 × 10⁶ |
| + damage and erosion | 7.83 × 10⁶ |
| + rigid floor | 6.19 × 10⁶ |
| + deformable contact | 3.60 × 10⁶ |
| bare, 4 threads | 2.42 × 10⁷ |
| + contact, 4 threads | 9.21 × 10⁶ |

The **kernel scales**: quadrupling the mesh leaves the per-element rate at 0.97× of the small
one, so the element loop is doing the same work per element however big the scene gets.
Contact does not — at the 4× mesh it drops another third.

The contact **search** is threaded and the accumulation is not, which is the only split that
keeps the answer reproducible: a search writes nothing but its own node's master, while a sum
into shared force arrays has an order that is part of the result. That, plus tightening the
broad-phase cell from twice the longest segment to 1.5 times it — the smallest radius that
still reaches every segment a node could be touching, since a touching node is at most
`√(0.5² + 1²) = 1.12` segment lengths from that segment's ends — took the threaded sandbox
case from 8.22 to 9.21 × 10⁶, with every number in the slug scene bit-identical.

The gap against the M0 cost model above (4.0 × 10⁷ single-threaded) is not a regression: that
was a small-strain elastic run with nothing touching anything. Finite strain rebuilds the
geometry every step, and contact is 2.5× the whole element kernel on its own. Both are what a
sandbox scene actually needs.

What that costs, for a 200 mm square of copper:

| cell | elements | step | hours per simulated second |
| --- | --- | --- | --- |
| 8 mm | 625 | 892 ns | 0.1 |
| 4 mm | 2 500 | 446 ns | 0.4 |
| 2 mm | 10 000 | 223 ns | 3.5 |
| 1 mm | 40 000 | 112 ns | 27.6 |
| 0.5 mm | 160 000 | 56 ns | 221 |

Cost goes as the **cube** of linear resolution — the mesh grows as the square and the timestep
shrinks with the cell — which is why quoting a frame rate for one mesh says nothing about
another. Six minutes of compute for one second of physics at a coarse 8 mm, and a day at
1 mm, **unscaled**. Mass scaling divides those by the factor it buys — which is how a 320 ms
topple fits into four seconds of wall clock — and threading and a GPU are what is left after
that. Commit-and-watch is not a design preference here; it is what the arithmetic allows.

### What the mesher refuses to do, and why

`Outline` rasterises a polygon onto a square lattice and keeps the cells whose centres fall
inside, so boundaries come out stair-stepped. A fitted mesher that followed the outline is the
obvious thing to want, and it is the wrong trade here. An explicit solve is paced by the
**shortest edge anywhere on the mesh**; trimming cells against an outline makes slivers at
every corner, unpredictably, as a function of where the shape happened to be drawn, and one
sliver in ten thousand elements slows every other element by the same factor. On a lattice
every edge is the cell size, so the stable step is known before the mesh exists. Stair-steps
are a resolution error that halves when the cell halves. A sliver is a performance cliff that
does not.

The lattice is anchored at the origin rather than at each shape's bounding box, so shapes
meshed together land on the same cells and share nodes where they touch. Shared nodes are a
*welded* joint, which is right for a clad plate and is not a substitute for contact.

### Contact, and two things that went wrong on the way

Penalty contact has one genuinely free parameter — the stiffness — and the usual way to build
it, from a bulk modulus and an element size, makes it a material property and therefore
ambiguous the moment steel meets copper. The mass-based form avoids that:

```
k = SCALE · m* / dt₀²        ⇒        ω = √SCALE / dt₀
```

The frequency it adds depends on nothing but `SCALE`. Central difference is stable to
`ω·dt < 2`, and every surface here is both slave and master so a pair is found twice, which
puts the bound at `SCALE < 2`. It is also formulation-agnostic for free: an axisymmetric
node's mass is already a ring mass, so `k` comes out as a ring stiffness with no special case.

Being stable is not being usable. Measured on the head-on collision the gate runs:

| SCALE | energy closes to | deepest penetration |
| --- | --- | --- |
| 0.01 | −3.2 % | 25.5 % of a cell |
| 0.05 | +2.8 % | 8.7 % |
| **0.10** | **+2.2 %** | **4.6 %** |
| 0.20 | +6.0 % | 3.5 % |
| 0.80 | +4.1 % | 1.6 % |
| 1.50 | +553 % | 0.4 % |
| 3.00 | diverges | |

Soft, and surfaces sink far enough in that nodes slide out of the projection windows they were
resolving — the spring vanishes without doing the work of releasing, and energy is *lost*.
Through the middle two decades it is a few per cent with no clean trend: it is the size of the
discrete release impulse, and it does not fall under mesh refinement. A **resting** contact is
much better than that suggests — a block dropped on another and slid across it closes to
0.15 % over six thousand steps.

**The audit was lying, and not about contact.** The collision appeared to gain 4 % and the
obvious suspect was the model that had just been added. It was not. Central difference carries
velocity at half steps and displacement at whole ones, so `KE + SE` adds two numbers from two
different instants, and the error is `(dt/2)·d(KE)/dt` — *first* order in the step. It costs
nothing on a body in rigid translation and it is real for anything vibrating, which is exactly
what a collision produces. A bar with nothing to hit, given a standing wave, already moves the
budget by **1.36 %**. The fix is exact and free: the velocity half a step further on is
`v + dt·f/m` with the force already assembled, so the centred velocity is `v + (dt/2)·f/m`.
Measured on that bar, the audit goes 1.36 % → 0.16 %, and the convergence order goes from 1.05
to 2.05. `centredKineticEnergy()` is now what every energy gate here reads.

**One spurious pair was worth more than the impact.** A block sliding across another gained
245 % of its energy, which read exactly like friction pumping. The cause was a contact between
two *perpendicular* surfaces — a top-face node of the lower block, normal straight up, against
the vertical side face of the upper one, normal straight out. Their dot product is exactly
zero, so an epsilon test admitted it the instant a hair of elastic deformation tipped the node
normal negative, and it came in with a depth equal to the whole sideways overlap rather than to
any penetration. One pair, 2.4 J, against a 2.5 J impact. The region behind a segment is a
half-plane and not a box, and that is the part that is easy to get wrong.

The fix is a real angle rather than an epsilon: surfaces have to face each other within about
104°. It costs nothing, because a contact seen at a glancing angle from one side is seen
head-on from the other. With it in place the sliding case closes to 0.15 %, and the contact
damper that had been added to suppress the symptom turned out to make every other measurement
slightly worse — it is still there, defaulted off, with a note saying why.

### Mass scaling, and the one motion it does not touch

Mass scaling trades honesty for wall clock: inflate every mass by `f`, every wave speed falls
by `√f`, and the stable step grows by `√f`. Four hundred times the mass is twenty times the
step. The usual objection is that it wrecks inertia, and it does — a body at a given speed
carries `f` times the momentum, so every impact is wrong by that factor.

What is easy to miss is that **the motion it is wanted for does not notice at all**. A block
rotating about its corner under its own weight turns at

```
α = torque / inertia = m·g·d / (m·k²) = g·d / k²
```

and the mass cancels exactly. Every purely gravity-driven motion has that property, because
the force and the inertia are the same mass. A toppling stack under 400× mass falls in exactly
the time the unscaled one would, which is what makes `topple` a legitimate scene rather than a
cartoon of one.

What it does change, besides momentum, is the acoustic impedance `ρ·c`, which grows as `√f`.
A block *landing* hits `√f` times harder than it should. That, not stability, is what sets the
factor: 20× on a 200 mm steel stack keeps the landing stresses under yield, and 100× would
not.

The textbook version of this is *selective* — inflate only the elements whose own CFL step is
below target, so the mesh's inertia is disturbed as little as possible. That matters on a mesh
with a few small elements among many large ones, which is exactly what a fitted mesher makes
and exactly what `Outline` refuses to make. **On a lattice, selective and uniform scaling are
the same thing**, and the blocky mesher pays for itself a second time.

### What the bracket showed

The scene was built on the theory that a stiff foot under a soft, off-centre leg would tip on
its own. Measured off the film, it does not: landing flat, the top of the leg turned **1.08°**
over the whole event and the section's centroid moved sideways by 16 µm. It squats, spreads by
17 %, and stays upright. The whole foot reaches the anvil in the same instant, so the reaction
is distributed under the entire base and there is no moment arm — mass being off-centre is not
a torque when every part of the contact is supported. Landing the same bracket on a corner at
18° turns it **27.4°** and kicks it off the anvil.

One element paid for that corner: 7.07 of plastic strain in the single cell that met the
plane, against 695 others all under 1.0. Nothing inverted and the mesh stayed valid, but a
contact carried by one cell is exactly the artefact that element deletion exists to resolve,
and it moved separation up the list. It also forced a change in `Film`: an automatic colour
ramp now spans the 99.5th percentile rather than the maximum and reports the true peak beside
it, because a ramp stretched to one wild element painted the entire body at the cold end and
said nothing had happened.

The constitutive core does not care about any of this. Plasticity, hardening, rate and
thermal terms, finite-strain kinematics, the defect field and the fracture-energy
regularisation are all formulation-agnostic and survive intact. This is a rewrite of the
harness around the engine, not of the engine.

**"Close, not exact" is doing real work here.** It licenses three things the validation code
deliberately refuses: deleting elements at full damage, which is how debris happens and which
every commercial hydrocode does; mass scaling, which buys back timestep at the cost of
slightly wrong inertia; and penalty contact, which the rigid anvil rejects on purity grounds
and which is the right choice for deformable pairs anyway. The validation suite is what then
says how much each of those cost.

## Next

1. ~~**Plane strain, a polygon mesher, and per-element material.**~~ Done. `Outline` meshes
   any polygon with holes, `Materials` paints one material per element, and `Impact` is the
   first scene built on `PLANE_STRAIN`. What is still missing from the *mesher* rather than
   the solver: pressure edges on a rasterised outline, and a way to name a boundary so a
   constraint can be attached to it without counting node indices.
2. ~~**Contact between deformable bodies.**~~ Done. `Contact` is penalty with a Coulomb cone,
   momentum-exact, and characterised above. `Outline.assemble` makes bodies that touch rather
   than weld, and `setGravity` makes them fall. What contact still lacks: it runs serial, and
   it is node-to-segment rather than segment-to-segment, which is where the residual few per
   cent of energy lives.
3. ~~**Separation.**~~ Done. `setErosion` deletes elements softened to nothing; nodes and mass
   stay, so momentum is untouched, and `erodedEnergy()` reports what left with the debris.
   `QuadMesh.surface(skip)` rebuilds the free surface around the hole so contact sees the new
   faces. What it cannot do is **shatter**: fragments need a brittle material, and the
   crack-band limit `h ≤ 2·E·G_f/σ_f²` puts the largest regularisable element for an alumina
   at a few microns against the millimetre a scene can afford. That is a real constraint, not
   a missing feature, and it is why the demo wall is ductile.
4. **Speed, and then interaction.** Mass scaling is done, the `topple` scene is what it
   bought, and the contact search is threaded. What is left, in order of what it buys:
   **thread the nodal update**, which any scene with a floor currently serialises because the
   wall carries running totals a dynamic schedule would reorder — per-chunk accumulators
   summed in chunk order would fix that; **segment-to-segment contact**, which is where the
   residual few per cent of collision energy lives and which would scale better than the
   node-to-segment search does; then the **GPU port**, which the structure-of-arrays layout
   was chosen for. Placing, dragging and launching come after, because commit-and-watch is
   what the numbers allow.
5. ~~**Runs have to accumulate.**~~ Done. `RunLog` appends one line per run — what went in,
   what came out — and `./gradlew history` reads it back. Parameters are the identity rather
   than a hash, so `lastMatching` answers "have I already run this" and the file stays
   readable. `impact`, `smash` and `topple` log themselves. What it is not yet: a way to
   *diff* two runs, which is the thing you actually want when a number moves.

Frozen, not abandoned — all three are below the accuracy bar the sandbox sets, and all three
are written up above: the burst pressure not converging under refinement, triaxiality missing
from the failure strain, and the unexplained 0.04 % in the burst deficit.

Known gaps, deliberate: plane-stress plasticity is not implemented (enforcing σ_zz = 0 through
a return map is a different algorithm, and nothing in the validation program needs it — the
solver refuses the combination rather than silently solving plane strain). Contact is
node-to-segment rather than segment-to-segment, and serial, so the residual few per cent of
energy across a collision stays where it is. Johnson–Cook **damage** (D1–D5) is not in — the
failure strain comes from the defect field and carries no triaxiality dependence, which is item
2 above. Plastic dissipation is integrated against reference element volumes, which is exact
only because this return map is exactly isochoric — and `Damage` does not break that, because
softening scales the *size* of the yield surface and leaves the flow direction alone. A model
with volumetric plasticity, or one that degraded the hydrostatic response so a crack could open
traction-free, would break it.

## Build notes

- Bytecode targets Java 21; any JDK ≥ 21 can build it. No dependencies outside JUnit.
- On the Windows machine this currently builds on there is no standalone JDK installed. The
  JetBrains Runtime shipped with IntelliJ is a JDK 25 and works:
  `export JAVA_HOME="/c/Program Files/JetBrains/IntelliJ IDEA 2026.2.3/jbr"`
- Run output goes to `runs/`, which is gitignored. Deliberately **not** `out/`: the IDE
  configuration claims that as a compiler output root, and a directory something else empties
  is not a place to keep traces. `-Dneofiz.runs=<dir>` moves it.
- The viewer needs a server: `node viewer/serve.js` serves pages from `viewer/` and film data
  from `runs/`. Opening the HTML as a `file://` URL will not work, because the film loads as a
  script beside it.
- Sweeps run one solver per thread and each solver single-threaded, which is the default.
  Do not call `ExplicitSolver.setThreads` inside a swept case: the solver's workers spin
  rather than sleeping between parallel regions, so eight solvers each asking for eight
  workers puts sixty-four runnable threads on eight cores and most of the machine goes into
  waiting for a barrier held by a thread that is not scheduled.
