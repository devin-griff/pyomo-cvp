# pyomo-cvp

[![PyPI](https://img.shields.io/pypi/v/pyomo-cvp.svg)](https://pypi.org/project/pyomo-cvp/)
[![Python versions](https://img.shields.io/pypi/pyversions/pyomo-cvp.svg)](https://pypi.org/project/pyomo-cvp/)
[![CI](https://github.com/devin-griff/pyomo-cvp/actions/workflows/ci.yml/badge.svg)](https://github.com/devin-griff/pyomo-cvp/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-BSD%203--Clause-blue.svg)](LICENSE)

Control vector parameterization for [pyomo.dae](https://pyomo.readthedocs.io/en/6.8.0/modeling_extensions/dae.html#collocation-transformation).

`pyomo.dae` makes control profiles piecewise-constant by adding linking
equality constraints (`reduce_collocation_points`), which keeps every
collocation-point copy of the control in the model plus one equation per
tied copy. `pyomo-cvp` does it by **elimination**: after any DAE
discretization, each declared control keeps only its profile's free values,
every other copy is substituted out of the model, and the component is
replaced under its own name. The model you solve is the model you meant:
no extra variables, no linking constraints.

On the classic race car problem (nfe=15, ncp=3, Lagrange-Radau):

|                                | control vars | linking constraints |
|--------------------------------|-------------:|--------------------:|
| `reduce_collocation_points`    |           46 |                  30 |
| `cvp.parameterize`             |           15 |                   0 |

## Install

```bash
pip install pyomo-cvp
```

## Usage

```python
import pyomo.environ as pyo
from pyomo_cvp import declare_profile, control_value

# ... build a pyomo.dae model with control m.u over ContinuousSet m.tau ...
declare_profile(m.u, wrt=m.tau, profile="piecewise_constant")

pyo.TransformationFactory("dae.collocation").apply_to(
    m, nfe=15, ncp=3, scheme="LAGRANGE-RADAU")
pyo.TransformationFactory("cvp.parameterize").apply_to(m)

# m.u now has exactly nfe members, one per finite element
pyo.SolverFactory("ipopt").solve(m)
control_value(m.u, 0.5)   # evaluate the profile at any time
```

The explicit form (no declaration) is equivalent:

```python
pyo.TransformationFactory("cvp.parameterize").apply_to(
    m, var=m.u, contset=m.tau, profile="piecewise_constant")
```

Works with any `pyomo.dae` discretization: Lagrange-Radau,
Lagrange-Legendre (where it also eliminates the dangling element-boundary
copies the constraint-based approach leaves unconstrained), or finite
difference. Controls may carry additional (non-time) indices.

## Profiles

- `'piecewise_constant'` --- one free value per finite element, indexed by
  the element's start time (`u[t0]` exists; the final time carries no
  control).
- `'piecewise_linear'` --- one free value per element boundary, continuous,
  interior points interpolated.
- `'collocation'` or `('collocation', k)`: the control is the element's
  collocation polynomial, with k free values per element (the last k
  collocation points, and k = ncp for the plain form) and Lagrange
  interpolation elsewhere. The elimination form of
  `reduce_collocation_points(ncp=k)`.

## Examples

Worked notebooks under [examples/](examples/):

- [racecar_cvp.ipynb](examples/racecar_cvp.ipynb): minimum-time race car,
  both invocation forms and all three profiles.
- [hicks_cvp.ipynb](examples/hicks_cvp.ipynb): the Hicks-Ray CSTR.
- [Quad_tank_cvp.ipynb](examples/Quad_tank_cvp.ipynb): the quadruple-tank
  process.

Install their dependencies with `pip install pyomo-cvp[examples]`.

## Citing

If you use this package, please also cite the pyomo.dae framework it builds on:

> Nicholson, B., Siirola, J.D., Watson, J.-P., Zavala, V.M., Biegler, L.T.
> (2018). pyomo.dae: a modeling and automatic discretization framework for
> optimization with differential and algebraic equations. *Mathematical
> Programming Computation* 10(2), 187-223.
> [doi:10.1007/s12532-017-0127-0](https://doi.org/10.1007/s12532-017-0127-0)

```bibtex
@article{nicholson2018pyomodae,
  author  = {Nicholson, Bethany and Siirola, John D. and Watson, Jean-Paul
             and Zavala, Victor M. and Biegler, Lorenz T.},
  title   = {pyomo.dae: a modeling and automatic discretization framework
             for optimization with differential and algebraic equations},
  journal = {Mathematical Programming Computation},
  volume  = {10},
  number  = {2},
  pages   = {187--223},
  year    = {2018},
  doi     = {10.1007/s12532-017-0127-0}
}
```

## Part of the DRTO stack

pyomo-cvp stands alone, but it is also the control-parameterization layer of
[DRTO](https://github.com/devin-griff/drto), a unified framework for dynamic
real-time optimization (NMPC, moving horizon estimation, and steady-state
RTO) built on Pyomo. In DRTO, `declare_control(m.u, profile=...)` delegates
to this package, so a declared model gets its control profiles without
calling pyomo-cvp directly. If you are parameterizing controls for a
receding-horizon controller, DRTO may be the layer you actually want.

## Maintainer

Maintained by [@devin-griff](https://github.com/devin-griff). Issues and pull
requests welcome.

## License

BSD 3-Clause License. See [LICENSE](LICENSE).
