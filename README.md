# pyomo-cvp

Control vector parameterization for [pyomo.dae](https://pyomo.readthedocs.io/en/stable/explanation/modeling_extensions/dae.html).

`pyomo.dae` makes control profiles piecewise-constant by adding linking
equality constraints (`reduce_collocation_points`), which leaves every
collocation-point copy of the control in the model plus one equation per tied
copy. `pyomo-cvp` does it by **elimination**: after any DAE discretization,
each declared control keeps one free value per finite element, every other
copy is substituted out, and the component is replaced under its own name.
The model you get is the model you meant: no extra variables, no linking
constraints.

> **Status: pre-alpha scaffold.** The `cvp.parameterize` transformation is
> registered but not yet implemented. See [ROADMAP.md](ROADMAP.md).

## Planned usage

```python
import pyomo.environ as pyo
from pyomo_cvp import declare_profile

# ... build a pyomo.dae model with control m.u over ContinuousSet m.tau ...
declare_profile(m.u, wrt=m.tau, profile="piecewise_constant")

pyo.TransformationFactory("dae.collocation").apply_to(
    m, nfe=15, ncp=3, scheme="LAGRANGE-RADAU")
pyo.TransformationFactory("cvp.parameterize").apply_to(m)
# m.u now has exactly nfe members, one per finite element
```

Works with any `pyomo.dae` discretization: Lagrange-Radau, Lagrange-Legendre,
or finite difference.

## License

Apache License 2.0. See [LICENSE](LICENSE).
