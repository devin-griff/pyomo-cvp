"""The explicit-arguments form of cvp.parameterize (no declare_profile),
shown with a finite-difference discretization: the transformation is the
same call regardless of scheme."""
import pyomo.environ as pyo
from pyomo.dae import ContinuousSet, DerivativeVar

from pyomo_cvp import control_value

L, R, NFE = 100.0, 0.001, 30

m = pyo.ConcreteModel()
m.tau = ContinuousSet(bounds=(0, 1))
m.x = pyo.Var(m.tau)
m.v = pyo.Var(m.tau, bounds=(0, None))
m.u = pyo.Var(m.tau, bounds=(-3, 1), initialize=0)
m.tf = pyo.Var(bounds=(1, None), initialize=10)
m.dx = DerivativeVar(m.x, wrt=m.tau)
m.dv = DerivativeVar(m.v, wrt=m.tau)


@m.Constraint(m.tau)
def ode_x(m, t):
    if t == m.tau.first():
        return pyo.Constraint.Skip
    return m.dx[t] == m.tf * m.v[t]


@m.Constraint(m.tau)
def ode_v(m, t):
    if t == m.tau.first():
        return pyo.Constraint.Skip
    return m.dv[t] == m.tf * (m.u[t] - R * m.v[t] ** 2)


m.ic_x = pyo.Constraint(expr=m.x[0] == 0)
m.ic_v = pyo.Constraint(expr=m.v[0] == 0)
m.fc_x = pyo.Constraint(expr=m.x[1] == L)
m.fc_v = pyo.Constraint(expr=m.v[1] == 0)
m.obj = pyo.Objective(expr=m.tf)

pyo.TransformationFactory("dae.finite_difference").apply_to(
    m, nfe=NFE, scheme="BACKWARD")

# explicit form: name the variable, the set, and the profile directly
pyo.TransformationFactory("cvp.parameterize").apply_to(
    m, var=m.u, contset=m.tau, profile="piecewise_constant")

print(f"u members after parameterization: {len(m.u)} (one per element)")
res = pyo.SolverFactory("ipopt").solve(m)
print(f"status: {res.solver.termination_condition}, "
      f"tf = {pyo.value(m.tf):.4f} s")
print(f"u(0.5) = {control_value(m.u, 0.5):.4f}")
