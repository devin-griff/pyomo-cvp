"""Phase 3: declare_profile discovery path."""
import pytest
import pyomo.environ as pyo
from pyomo.dae import ContinuousSet, DerivativeVar

import pyomo_cvp
from pyomo_cvp import declare_profile

from test_parameterize import racecar, NFE, NCP, R

ipopt_available = pyo.SolverFactory("ipopt").available(False)
needs_ipopt = pytest.mark.skipif(not ipopt_available, reason="ipopt not available")


def discretize(m, **kw):
    kw.setdefault("nfe", NFE)
    kw.setdefault("ncp", NCP)
    kw.setdefault("scheme", "LAGRANGE-RADAU")
    pyo.TransformationFactory("dae.collocation").apply_to(m, **kw)
    return m


@needs_ipopt
def test_declared_matches_explicit():
    m1 = discretize(racecar())
    pyo.TransformationFactory("cvp.parameterize").apply_to(
        m1, var=m1.u, contset=m1.tau
    )
    r1 = pyo.SolverFactory("ipopt").solve(m1)
    assert r1.solver.termination_condition == pyo.TerminationCondition.optimal

    m2 = racecar()
    declare_profile(m2.u, wrt=m2.tau, profile="piecewise_constant")
    discretize(m2)
    pyo.TransformationFactory("cvp.parameterize").apply_to(m2)
    assert len(m2.u) == NFE
    r2 = pyo.SolverFactory("ipopt").solve(m2)
    assert r2.solver.termination_condition == pyo.TerminationCondition.optimal
    assert pyo.value(m2.tf) == pytest.approx(pyo.value(m1.tf), rel=1e-8)


def two_controls():
    m = pyo.ConcreteModel()
    m.tau = ContinuousSet(bounds=(0, 1))
    m.x = pyo.Var(m.tau)
    m.v = pyo.Var(m.tau, bounds=(0, None))
    m.throttle = pyo.Var(m.tau, bounds=(0, 1), initialize=0)
    m.brake = pyo.Var(m.tau, bounds=(0, 3), initialize=0)
    m.tf = pyo.Var(bounds=(1, None), initialize=10)
    m.dx = DerivativeVar(m.x, wrt=m.tau)
    m.dv = DerivativeVar(m.v, wrt=m.tau)

    def ode_x(m, t):
        if t == m.tau.first():
            return pyo.Constraint.Skip
        return m.dx[t] == m.tf * m.v[t]

    def ode_v(m, t):
        if t == m.tau.first():
            return pyo.Constraint.Skip
        return m.dv[t] == m.tf * (
            m.throttle[t] - m.brake[t] - R * m.v[t] ** 2
        )

    m.ode_x = pyo.Constraint(m.tau, rule=ode_x)
    m.ode_v = pyo.Constraint(m.tau, rule=ode_v)
    m.ic_x = pyo.Constraint(expr=m.x[0] == 0)
    m.ic_v = pyo.Constraint(expr=m.v[0] == 0)
    m.fc_x = pyo.Constraint(expr=m.x[1] == 100.0)
    m.fc_v = pyo.Constraint(expr=m.v[1] == 0)
    m.obj = pyo.Objective(expr=m.tf)
    return m


@needs_ipopt
def test_two_declared_controls():
    m = two_controls()
    declare_profile(m.throttle, wrt=m.tau)
    declare_profile(m.brake, wrt=m.tau)
    discretize(m)
    pyo.TransformationFactory("cvp.parameterize").apply_to(m)
    assert len(m.throttle) == NFE and len(m.brake) == NFE
    r = pyo.SolverFactory("ipopt").solve(m)
    assert r.solver.termination_condition == pyo.TerminationCondition.optimal


def test_no_declarations_error():
    m = discretize(racecar())
    with pytest.raises(RuntimeError, match="no control profile declarations"):
        pyo.TransformationFactory("cvp.parameterize").apply_to(m)


def test_declarations_consumed():
    m = racecar()
    declare_profile(m.u, wrt=m.tau)
    discretize(m)
    xf = pyo.TransformationFactory("cvp.parameterize")
    xf.apply_to(m)
    with pytest.raises(RuntimeError, match="already applied"):
        xf.apply_to(m)
