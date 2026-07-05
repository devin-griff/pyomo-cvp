# pyomo-cvp roadmap

Control vector parameterization for `pyomo.dae`. Replaces the
constraint-based `reduce_collocation_points` mechanism with variable
elimination: after any DAE discretization, each declared control keeps one
free value per finite element (component replaced under its own name), and
every other copy is substituted out of the model. No extra variables, no
linking constraints, correct model before the solve.

Registry: `TransformationFactory('cvp.parameterize')`. Declaration API:
`declare_profile(var, wrt=contset, profile=...)`.

## Phase 0 — scaffold
Repo, packaging (hatchling, src layout, Apache-2.0), CI, registration smoke
test. Acceptance: package installs; `cvp.parameterize` resolves.

## Phase 1 — core transformation, explicit API
`apply_to(m, var=..., contset=..., profile='piecewise_constant')` for
Lagrange-Radau collocation. Representative point per finite element,
component replacement under the same name, `replace_expressions` over all
active constraints and objectives (including expanded Integrals). Guards:
not-yet-discretized, DerivativeVar attached, double application.
Acceptance: race car model matches the `reduce_collocation_points` optimum
(IPOPT), with control copies 46 -> 15 and tie constraints 30 -> 0.

## Phase 2 — scheme coverage
Lagrange-Legendre (element-boundary copies substituted too, eliminating the
dangling-variable wart) and `dae.finite_difference` (uniform API; t0
ownership). Controls with extra (non-time) indices. Acceptance: same-optimum
tests per scheme; a test that Legendre leaves no unreferenced control copies.

## Phase 3 — declare_profile
Build-time annotation (inert metadata) + no-argument `apply_to(m)` discovery;
multiple controls, per-control profiles. Acceptance: annotated race car
solves identically to the explicit-args path.

## Phase 4 — more profiles
`'piecewise_linear'` (representatives at element boundaries, interior copies
replaced by interpolation expressions) and `('reduced_collocation', k)` for
full parity with `reduce_collocation_points(ncp=k)`. Acceptance: pwl matches
an independently hand-built pwl model; reduced-collocation matches rcp's
optimum with fewer variables.

## Phase 5 — polish and release
`u_at(t)` reporting helper, README with the race car before/after size and
solution table, example notebook, CHANGELOG, PyPI trusted publishing, v0.1.0.

## Milestone validation
The ndcbe race car three ways — reduce_collocation_points baseline, explicit
cvp, declared cvp — identical optimal final time, with the model-size table.
