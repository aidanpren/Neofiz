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
kernel that is bit-for-bit identical at any thread count, and a mesh-independent Weibull defect
field. 197 tests green.

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

**M2's first half is in: the defect field exists, it is mesh-independent, and its size effect
matches weakest-link theory to 0.2 %.** One float per element. The plan names one trap here and
it is real, but there is a **second one facing the other way** that only shows up if you
actually run the refinement study: volume normalisation on its own makes the part *stronger*
without bound as the mesh refines, 35 % over a thousandfold, because below the correlation
length neighbouring elements are not independent chances. Nothing in the solver reads the field
yet. See [M2](#m2--the-defect-field).

```
./gradlew run     # the M0 report
./gradlew test    # the validation gates
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

#### What is not done

Nothing in the solver reads this number. Failure and the scattered burst pressure it should
produce are M2's second half, and the softening that carries them is where the mesh dependence
comes back — see [Next](#next).

## Layout

```
core/     Formulation (the axisym / plane-strain / plane-stress flag), Material
          JohnsonCook — flow stress: strain, rate and thermal terms
          Weibull — the weakest-link algebra; minimumOf(n) is the whole of it,
              and forVolume is a call to it rather than the same formula twice
          Normal — Cody erfc, because the copula's error lands on the marginal
          DefectField — correlated failure strain from a stateless hash;
              effectiveVolume is where the ring and the resolution floor live
mesh/     QuadMesh — flat primitive arrays, the layout that ports to a GPU kernel
              cylinderWall for the annulus, solidCylinder for the Taylor specimen
              ringVolume — Pappus, 2*pi*rc*A, the volume a defect draw normalises by
solver/   ExplicitSolver — the hot loop
          Parallel — spin-barrier worker pool; chunks from a shared cursor,
              which is safe only because nothing accumulates per chunk
          Integration — full 2x2 or reduced + hourglass
          Kinematics — small strain or updated Lagrangian
          J2 — radial-return plasticity; closed form for linear hardening,
              bracketed Newton for Johnson-Cook. J2.Flow bundles the uniforms
          Corotational — objective incremental rotation
          RigidWall — kinematic anvil with a Coulomb cone; stick/slip is a
              return map, same shape as J2
validate/ Lamé and ElasticPlastic closed forms, the cylinder cases,
              TaylorImpactCase — the M1 gate, which has no closed form,
              Burst — the instability condition, the P(ε) curve, Considère,
                  Barlow, and the elastic correction to a rigid-plastic oracle
              BurstCase — traverses the load maximum and measures the wall's
                  capacity rather than the pressure applied to it
              MeshConvergence — Richardson extrapolation and the GCI band
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

## Next

1. **M2's second half — let the solver read the field.** The field is built, tested and
   mesh-independent; nothing in `ExplicitSolver` looks at it. What it needs is a damage variable
   per element driven by the ratio of accumulated plastic strain to that element's local failure
   strain, and the charter's ban on binary failure means it has to enter as **progressive
   softening**, not deletion — an element that has reached its failure strain stops carrying
   deviatoric stress over some finite strain rather than vanishing between two timesteps.
   Softening is the part that needs care: it makes the tangent negative, which is exactly the
   regime where an explicit solve localises into one element and the answer becomes the mesh.
   The standard fix is to regularise by fracture energy per unit area rather than per unit
   volume, so the softening slope depends on the element size — and that lands in the same
   territory as the resolution floor already in `DefectField.effectiveVolume`. The payoff is
   the plan's actual claim: the same tube, run twice, bursts at two different pressures and in
   two different places.
2. **The 0.04 % nobody has accounted for in burst.** It is independent of wall thickness and
   of loading rate — both checked — so it is neither the thin-wall assumption nor the harness.
   The standing suspect is that the elastic response is *hypoelastic*, a rate form integrated
   along the path, whose departure from an exact hyperelastic law is of order elastic strain ×
   total strain: 0.002 × 0.1, which is the size observed. That is a hypothesis with the right
   magnitude and the right invariances, not a measurement. Settling it means a hyperelastic
   volumetric split, which is wanted anyway before anything is trusted past ~50 % strain.
3. **The pressure–volume curve as an output, not just its peak.** The burst harness already
   traces `P(ε_θ)` through the maximum and throws all of it away but one point. That curve *is*
   the bulge-as-positive-feedback story the charter promises, and it is the first thing in this
   project a person could be shown rather than told.
4. **The GPU port**, which is the only thing still holding gate 2 open, and which now has a
   target rather than an aspiration: 1.6–7.8× over ten CPU threads on a bandwidth-bound kernel.

Known gaps, deliberate: plane-stress plasticity is not implemented (enforcing σ_zz = 0 through
a return map is a different algorithm, and nothing in the validation program needs it — the
solver refuses the combination rather than silently solving plane strain). Contact is a rigid
plane -- Coulomb friction is in, but the plan's penalty model for deformable-on-deformable
pairs is not, and nothing before M5 needs it. Johnson–Cook **damage** (D1–D5) is not in — no case before
M5 reaches it, and the charter's ban on binary failure means it has to arrive as progressive
softening rather than deletion. Plastic dissipation is integrated against reference element
volumes, which is exact only because this return map is exactly isochoric — that argument
stops holding for volumetric plasticity or damage.

## Build notes

- Bytecode targets Java 21; any JDK ≥ 21 can build it.
- On this machine the JDKs under `~/Library/Java` cannot be executed by non-interactive
  processes (macOS privacy protection on `~/Library`). The system JDK at
  `/Library/Java/JavaVirtualMachines/jdk-25.jdk` works:
  `export JAVA_HOME=/Library/Java/JavaVirtualMachines/jdk-25.jdk/Contents/Home`
- `git` needs the Xcode licence accepted first: `sudo xcodebuild -license`.
