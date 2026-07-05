"""Phase 1: cvp.parameterize vs reduce_collocation_points on the race car."""
import pytest
import pyomo.environ as pyo
from pyomo.dae import ContinuousSet, DerivativeVar

import pyomo_cvp  # noqa: F401

ipopt_available = pyo.SolverFactory("ipopt").available(False)
needs_ipopt = pytest.mark.skipif(not ipopt_available, reason="ipopt not available")

NFE, NCP = 15, 3
L, R = 100.0, 0.001


def racecar():
    m = pyo.ConcreteModel()
    m.tau = ContinuousSet(bounds=(0, 1))
    m.x = pyo.Var(m.tau)
    m.v = pyo.Var(m.tau, bounds=(0, None))
    m.u = pyo.Var(m.tau, bounds=(-3, 1), initialize=0)
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
        return m.dv[t] == m.tf * (m.u[t] - R * m.v[t] ** 2)

    m.ode_x = pyo.Constraint(m.tau, rule=ode_x)
    m.ode_v = pyo.Constraint(m.tau, rule=ode_v)
    m.ic_x = pyo.Constraint(expr=m.x[0] == 0)
    m.ic_v = pyo.Constraint(expr=m.v[0] == 0)
    m.fc_x = pyo.Constraint(expr=m.x[1] == L)
    m.fc_v = pyo.Constraint(expr=m.v[1] == 0)
    m.obj = pyo.Objective(expr=m.tf)
    return m


def discretize(m):
    pyo.TransformationFactory("dae.collocation").apply_to(
        m, nfe=NFE, ncp=NCP, scheme="LAGRANGE-RADAU"
    )
    return m


def counts(m):
    nv = sum(
        1 for v in m.component_data_objects(pyo.Var, active=True) if not v.fixed
    )
    nc = sum(1 for _ in m.component_data_objects(pyo.Constraint, active=True))
    return nv, nc


def test_model_size_reduction():
    # rcp baseline must use the same transformation instance that discretized
    m_rcp2 = racecar()
    d2 = pyo.TransformationFactory("dae.collocation")
    d2.apply_to(m_rcp2, nfe=NFE, ncp=NCP, scheme="LAGRANGE-RADAU")
    d2.reduce_collocation_points(m_rcp2, var=m_rcp2.u, ncp=1, contset=m_rcp2.tau)

    m_cvp = discretize(racecar())
    pyo.TransformationFactory("cvp.parameterize").apply_to(
        m_cvp, var=m_cvp.u, contset=m_cvp.tau
    )

    assert len(m_cvp.u) == NFE
    v_rcp, c_rcp = counts(m_rcp2)
    v_cvp, c_cvp = counts(m_cvp)
    assert v_cvp < v_rcp
    assert c_cvp < c_rcp
    # exactly one u value per element; no interpolation constraint list
    assert m_cvp.find_component("u_interpolation_constraints") is None


@needs_ipopt
def test_same_optimum_as_rcp():
    m1 = racecar()
    d = pyo.TransformationFactory("dae.collocation")
    d.apply_to(m1, nfe=NFE, ncp=NCP, scheme="LAGRANGE-RADAU")
    d.reduce_collocation_points(m1, var=m1.u, ncp=1, contset=m1.tau)
    r1 = pyo.SolverFactory("ipopt").solve(m1)
    assert r1.solver.termination_condition == pyo.TerminationCondition.optimal

    m2 = discretize(racecar())
    pyo.TransformationFactory("cvp.parameterize").apply_to(
        m2, var=m2.u, contset=m2.tau
    )
    r2 = pyo.SolverFactory("ipopt").solve(m2)
    assert r2.solver.termination_condition == pyo.TerminationCondition.optimal

    assert pyo.value(m2.tf) == pytest.approx(pyo.value(m1.tf), rel=1e-6)
    # element controls agree with the rcp element values (the kept points)
    for t in m2.u:
        assert pyo.value(m2.u[t]) == pytest.approx(pyo.value(m1.u[t]), abs=1e-6)


def test_guard_not_discretized():
    m = racecar()
    with pytest.raises(RuntimeError, match="has not been discretized"):
        pyo.TransformationFactory("cvp.parameterize").apply_to(
            m, var=m.u, contset=m.tau
        )


def test_guard_derivative_var():
    m = discretize(racecar())
    with pytest.raises(ValueError, match="DerivativeVar"):
        pyo.TransformationFactory("cvp.parameterize").apply_to(
            m, var=m.v, contset=m.tau
        )


def test_guard_double_apply():
    m = discretize(racecar())
    xf = pyo.TransformationFactory("cvp.parameterize")
    xf.apply_to(m, var=m.u, contset=m.tau)
    with pytest.raises(ValueError, match="not indexed by"):
        xf.apply_to(m, var=m.u, contset=m.tau)
