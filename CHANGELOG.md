# Changelog

All notable changes to this project are documented here. The format is based
on [Keep a Changelog](https://keepachangelog.com/), and the project adheres
to [Semantic Versioning](https://semver.org/).

## [Unreleased]

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

[Unreleased]: https://github.com/devin-griff/pyomo-cvp/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/devin-griff/pyomo-cvp/releases/tag/v0.1.0
