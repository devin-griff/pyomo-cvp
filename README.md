# pyomo-cvp

[![PyPI](https://img.shields.io/pypi/v/pyomo-cvp.svg)](https://pypi.org/project/pyomo-cvp/)
[![Python versions](https://img.shields.io/pypi/pyversions/pyomo-cvp.svg)](https://pypi.org/project/pyomo-cvp/)
[![CI](https://github.com/devin-griff/pyomo-cvp/actions/workflows/ci.yml/badge.svg)](https://github.com/devin-griff/pyomo-cvp/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

Control vector parameterization for [pyomo.dae](https://pyomo.readthedocs.io/en/stable/explanation/modeling_extensions/dae.html).

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

This matters for NLP solvers such as IPOPT, which have no presolve to strip
redundant variables and equalities.

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

- `'piecewise_constant'` --- one free value per finite element.
- `'piecewise_linear'` --- one free value per element boundary, continuous,
  interior points interpolated.
- `('reduced_collocation', k)` --- k free values per element (the last k
  collocation points), Lagrange interpolation elsewhere; the elimination
  form of `reduce_collocation_points(ncp=k)`.

See [examples/racecar_cvp.ipynb](examples/racecar_cvp.ipynb) for a
complete worked example showing both forms and all three profiles.

## License

Apache License 2.0. See [LICENSE](LICENSE).
