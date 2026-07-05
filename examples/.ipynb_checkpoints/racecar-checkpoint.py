"""Minimum-time race car (pyomo.dae) with a piecewise-constant control via
pyomo-cvp. Compare model sizes with the reduce_collocation_points approach."""
import pyomo.environ as pyo
from pyomo.dae import ContinuousSet, DerivativeVar

from pyomo_cvp import control_value, declare_profile

L, R, NFE, NCP = 200.0, 0.001, 15, 3

m = pyo.ConcreteModel()
m.tau = ContinuousSet(bounds=(0, 1))
m.x = pyo.Var(m.tau)
m.v = pyo.Var(m.tau, bounds=(0, None))
m.u = pyo.Var(m.tau, bounds=(-3, 1), initialize=0)
m.tf = pyo.Var(bounds=(1, None), initialize=10)
m.dx = DerivativeVar(m.x, wrt=m.tau)
m.dv = DerivativeVar(m.v, wrt=m.tau)

declare_profile(m.u, wrt=m.tau, profile="piecewise_constant")


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

pyo.TransformationFactory("dae.collocation").apply_to(
    m, nfe=NFE, ncp=NCP, scheme="LAGRANGE-RADAU")
pyo.TransformationFactory("cvp.parameterize").apply_to(m)

print(f"control members after parameterization: {len(m.u)} (one per element)")
res = pyo.SolverFactory("ipopt").solve(m)
print(f"status: {res.solver.termination_condition}, tf = {pyo.value(m.tf):.3f} s")
print("u(0.5) =", control_value(m.u, 0.5))
