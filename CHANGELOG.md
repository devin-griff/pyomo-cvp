# Changelog

All notable changes to this project are documented here. The format is based
on [Keep a Changelog](https://keepachangelog.com/), and the project adheres
to [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.6.3] - 2026-07-17

### Changed

- Piecewise-constant control references resolve by what contains them.
  Model equations (constraints carrying a DerivativeVar or indexed over
  the ContinuousSet, and Expressions likewise) take the value before the
  jump at an element boundary and the last value at the final time.
  Objectives and cost constraints charge the decision made at that
  instant, and a cost reference at the final time errors, since no move
  starts there.

### Removed

- The `final_node` option on `declare_profile` and `cvp.parameterize`:
  with references classified by their containing expression, nothing is
  left for a flag to select.

## [0.6.2] - 2026-07-17

### Changed

- `final_node` moves to the declaration: `declare_profile(..., final_node=...)`
  records it per control, so one `parameterize` call can mix conventions.
  Re-declaring a variable replaces its pending declaration (last wins), the
  sanctioned way to update a convention before it is applied. The call-level
  option remains for the explicit form only (`var=`, `contset=`); passing it
  in declaration mode errors. The final-node tests run under Radau,
  Legendre, and finite differences: the convention comes from the profile,
  not the scheme.

## [0.6.1] - 2026-07-17

### Added

- `final_node` option on `cvp.parameterize` (piecewise-constant only):
  `'remove'` (default, the terminal-horizon convention) errors on a control
  reference at the final grid point, where no move exists; `'keep'` (the
  horizon-continues convention) defines the control there as the held last
  move, so a DAE's algebraic equations at the final collocation point
  resolve like the collocation equation beside them.

## [0.6.0] - 2026-07-16

### Changed

- Renamed the `('reduced_collocation', k)` profile to `('collocation', k)`:
  with k = ncp nothing is reduced, and the profile's meaning is the control
  as the element's collocation polynomial. The plain string `'collocation'`
  is now accepted and resolves k to the discretization's ncp at transform
  time. Breaking rename with no deprecation shim: the old name raises
  `ValueError`.

## [0.5.0] - 2026-07-14

### Changed

- Relicensed from Apache-2.0 to BSD-3-Clause, matching the Pyomo and
  scientific-Python ecosystem this builds on.
- Minimum Pyomo is now 6.8.1 and minimum Python is now 3.10. The transform
  stores its per-block state through `Block.private_data` (added in Pyomo
  6.8.1) instead of ad-hoc block attributes.
- The `cvp.parameterize` transformation validates its options through a
  `ConfigDict`, so an unknown keyword option now raises `ValueError` instead
  of being silently ignored.

### Internal

- Aligned with Pyomo's contribution conventions: NumPy-style docstrings on the
  public and private API, Black formatting with Pyomo's settings, per-file BSD
  license headers, a shared test-helper module (tests no longer import one
  another), and coverage, a minimum-dependency check, and a spell-check step in
  CI.
- CI cuts a GitHub Release from the changelog on each version tag.

## [0.4.0] - 2026-07-13

### Fixed

- Control references in algebraic expressions (objectives, path constraints,
  non-differential named expressions) no longer shift for piecewise-constant
  controls. A reference to element k's control `u[fe[k]]` stays element k's,
  instead of being rewritten to element k-1's `u[fe[k-1]]`, which silently
  penalized the wrong element. The rewrite now splits by expression kind:
  differential equations (those carrying a `DerivativeVar`: the user's ODEs
  and pyomo.dae's `disc_eq`) are unchanged; algebraic expressions resolve a
  control reference to the element it sits in. piecewise-linear and
  reduced-collocation are unaffected, their maps already coincide.

### Changed

- **Breaking**: an algebraic reference to the final-time piecewise-constant
  control (`u[fe[N]]`, which does not exist: N elements carry N controls
  indexed by the element starts 0..N-1) now raises `ValueError` instead of
  being silently rewritten to the last element and double-counted. Reference
  controls over the elements, not the final boundary, and penalize the
  terminal state on its own.

## [0.3.0] - 2026-07-12

### Changed

- `declare_profile` accepts one or more control Vars per call
  (`declare_profile(m.v1, m.v2, wrt=m.t)`), with the profile settings
  applying to every variable listed; `wrt` is now keyword-only, and the
  old positional form raises a TypeError naming the fix.

## [0.2.0] - 2026-07-10

### Changed

- **Breaking**: piecewise-constant controls now index by each element's
  START time (the standard CVP convention, u(t) = u_k on [t_k, t_k+1)),
  not its end. `u[t0]` now exists and the final time carries no control;
  code written for 0.1.0 that indexed controls by element end times must
  shift its labels. Optima, profiles, and `control_value` evaluations are
  numerically unchanged: this is a relabeling of the free variables only.
- `control_value` now raises `ValueError` for `t` outside the control's
  domain instead of silently holding the end values.

## [0.1.0] - 2026-07-05

First release (alpha).

### Added
- `cvp.parameterize` transformation: reduce a control variable over a
  discretized ContinuousSet to its profile's free values by substitution,
  replacing the component under its own name --- no linking constraints, no
  leftover copies. Works after any `pyomo.dae` discretization
  (Lagrange-Radau, Lagrange-Legendre, finite difference); under Legendre the
  otherwise-dangling element-boundary copies are eliminated as well.
  Controls may carry additional non-time indices.
- Profiles: `'piecewise_constant'` (one value per finite element),
  `'piecewise_linear'` (continuous, one value per element boundary), and
  `('reduced_collocation', k)` (elimination form of
  `reduce_collocation_points(ncp=k)`).
- `declare_profile(var, wrt=..., profile=...)`: declare profiles at model
  build time; `apply_to(model)` with no arguments applies all declarations.
- `control_value(var, t)`: evaluate a parameterized control at any time.
- Guards: not-yet-discretized set, DerivativeVar-attached control, double
  application, reduced collocation without collocation or with k > ncp.

[Unreleased]: https://github.com/devin-griff/pyomo-cvp/compare/v0.6.3...HEAD
[0.6.3]: https://github.com/devin-griff/pyomo-cvp/compare/v0.6.2...v0.6.3
[0.6.2]: https://github.com/devin-griff/pyomo-cvp/compare/v0.6.1...v0.6.2
[0.6.1]: https://github.com/devin-griff/pyomo-cvp/compare/v0.6.0...v0.6.1
[0.6.0]: https://github.com/devin-griff/pyomo-cvp/compare/v0.5.0...v0.6.0
[0.5.0]: https://github.com/devin-griff/pyomo-cvp/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/devin-griff/pyomo-cvp/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/devin-griff/pyomo-cvp/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/devin-griff/pyomo-cvp/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/devin-griff/pyomo-cvp/releases/tag/v0.1.0
