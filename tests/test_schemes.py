"""Phase 2: Legendre, finite difference, and multi-indexed controls."""
import pytest
import pyomo.environ as pyo
from pyomo.dae import ContinuousSet, DerivativeVar

import pyomo_cvp  # noqa: F401

from test_parameterize import racecar, NFE, NCP, R

ipopt_available = pyo.SolverFactory("ipopt").available(False)
needs_ipopt = pytest.mark.skipif(not ipopt_available, reason="ipopt not available")


@needs_ipopt
def test_legendre_matches_rcp():
    m1 = racecar()
    d = pyo.TransformationFactory("dae.collocation")
    d.apply_to(m1, nfe=NFE, ncp=NCP, scheme="LAGRANGE-LEGENDRE")
    d.reduce_collocation_points(m1, var=m1.u, ncp=1, contset=m1.tau)
    r1 = pyo.SolverFactory("ipopt").solve(m1)
    assert r1.solver.termination_condition == pyo.TerminationCondition.optimal

    m2 = racecar()
    pyo.TransformationFactory("dae.collocation").apply_to(
        m2, nfe=NFE, ncp=NCP, scheme="LAGRANGE-LEGENDRE"
    )
    pyo.TransformationFactory("cvp.parameterize").apply_to(
        m2, var=m2.u, contset=m2.tau
    )
    assert len(m2.u) == NFE
    r2 = pyo.SolverFactory("ipopt").solve(m2)
    assert r2.solver.termination_condition == pyo.TerminationCondition.optimal
    assert pyo.value(m2.tf) == pytest.approx(pyo.value(m1.tf), rel=1e-6)


@needs_ipopt
def test_finite_difference():
    # baseline: plain backward finite difference, no reduction available
    m1 = racecar()
    pyo.TransformationFactory("dae.finite_difference").apply_to(
        m1, nfe=NFE, scheme="BACKWARD"
    )
    r1 = pyo.SolverFactory("ipopt").solve(m1)
    assert r1.solver.termination_condition == pyo.TerminationCondition.optimal

    m2 = racecar()
    pyo.TransformationFactory("dae.finite_difference").apply_to(
        m2, nfe=NFE, scheme="BACKWARD"
    )
    pyo.TransformationFactory("cvp.parameterize").apply_to(
        m2, var=m2.u, contset=m2.tau
    )
    assert len(m2.u) == NFE          # nfe+1 copies -> nfe (t0 copy eliminated)
    r2 = pyo.SolverFactory("ipopt").solve(m2)
    assert r2.solver.termination_condition == pyo.TerminationCondition.optimal
    assert pyo.value(m2.tf) == pytest.approx(pyo.value(m1.tf), rel=1e-6)


def twin_engine():
    """Race car with a two-component control u[t, engine]."""
    m = pyo.ConcreteModel()
    m.tau = ContinuousSet(bounds=(0, 1))
    m.E = pyo.Set(initialize=[1, 2])
    m.x = pyo.Var(m.tau)
    m.v = pyo.Var(m.tau, bounds=(0, None))
    m.u = pyo.Var(m.tau, m.E, bounds=(-3, 1), initialize=0)
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
            0.6 * m.u[t, 1] + 0.4 * m.u[t, 2] - R * m.v[t] ** 2
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
def test_multi_indexed_control():
    m1 = twin_engine()
    d = pyo.TransformationFactory("dae.collocation")
    d.apply_to(m1, nfe=NFE, ncp=NCP, scheme="LAGRANGE-RADAU")
    d.reduce_collocation_points(m1, var=m1.u, ncp=1, contset=m1.tau)
    r1 = pyo.SolverFactory("ipopt").solve(m1)
    assert r1.solver.termination_condition == pyo.TerminationCondition.optimal

    m2 = twin_engine()
    pyo.TransformationFactory("dae.collocation").apply_to(
        m2, nfe=NFE, ncp=NCP, scheme="LAGRANGE-RADAU"
    )
    pyo.TransformationFactory("cvp.parameterize").apply_to(
        m2, var=m2.u, contset=m2.tau
    )
    assert len(m2.u) == NFE * 2
    r2 = pyo.SolverFactory("ipopt").solve(m2)
    assert r2.solver.termination_condition == pyo.TerminationCondition.optimal
    assert pyo.value(m2.tf) == pytest.approx(pyo.value(m1.tf), rel=1e-6)
