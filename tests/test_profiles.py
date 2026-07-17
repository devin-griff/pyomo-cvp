# Copyright (c) 2026 Devin Griffith
# SPDX-License-Identifier: BSD-3-Clause
"""Phase 4: piecewise_linear and collocation profiles."""
import pytest
import pyomo.environ as pyo

import pyomo_cvp  # noqa: F401
from helpers import racecar, NFE, NCP, needs_ipopt


def discretize(m, scheme="LAGRANGE-RADAU"):
    pyo.TransformationFactory("dae.collocation").apply_to(
        m, nfe=NFE, ncp=NCP, scheme=scheme
    )
    return m


@needs_ipopt
def test_piecewise_linear_matches_handbuilt():
    # hand-built baseline: full u plus explicit interpolation constraints
    m1 = discretize(racecar())
    fe = list(m1.tau.get_finite_elements())
    interior = [t for t in m1.tau if t not in fe]

    def interp(m, t):
        import bisect

        i = bisect.bisect_left(fe, t) - 1
        a, b = fe[i], fe[i + 1]
        w = (t - a) / (b - a)
        return m.u[t] == (1 - w) * m.u[a] + w * m.u[b]

    m1.pwl_tie = pyo.Constraint(interior, rule=interp)
    r1 = pyo.SolverFactory("ipopt").solve(m1)
    assert r1.solver.termination_condition == pyo.TerminationCondition.optimal

    m2 = discretize(racecar())
    pyo.TransformationFactory("cvp.parameterize").apply_to(
        m2, var=m2.u, contset=m2.tau, profile="piecewise_linear"
    )
    assert len(m2.u) == NFE + 1
    r2 = pyo.SolverFactory("ipopt").solve(m2)
    assert r2.solver.termination_condition == pyo.TerminationCondition.optimal
    assert pyo.value(m2.tf) == pytest.approx(pyo.value(m1.tf), rel=1e-6)


@needs_ipopt
@pytest.mark.parametrize("scheme", ["LAGRANGE-RADAU", "LAGRANGE-LEGENDRE"])
def test_collocation_matches_rcp(scheme):
    k = 2
    m1 = racecar()
    d = pyo.TransformationFactory("dae.collocation")
    d.apply_to(m1, nfe=NFE, ncp=NCP, scheme=scheme)
    d.reduce_collocation_points(m1, var=m1.u, ncp=k, contset=m1.tau)
    r1 = pyo.SolverFactory("ipopt").solve(m1)
    assert r1.solver.termination_condition == pyo.TerminationCondition.optimal

    m2 = discretize(racecar(), scheme=scheme)
    pyo.TransformationFactory("cvp.parameterize").apply_to(
        m2, var=m2.u, contset=m2.tau, profile=("collocation", k)
    )
    assert len(m2.u) == NFE * k
    # variable bounds on eliminated copies survive as profile-bound rows
    assert m2.find_component("u_profile_bounds") is not None
    r2 = pyo.SolverFactory("ipopt").solve(m2)
    assert r2.solver.termination_condition == pyo.TerminationCondition.optimal
    if scheme == "LAGRANGE-RADAU":
        assert pyo.value(m2.tf) == pytest.approx(pyo.value(m1.tf), rel=1e-6)
    else:
        # Under Legendre, rcp leaves the element-boundary control copies as
        # free independent variables (its dangling-variable wart); we tie them
        # to the element polynomial, so the models differ slightly by design.
        assert pyo.value(m2.tf) == pytest.approx(pyo.value(m1.tf), rel=2e-3)


@needs_ipopt
def test_collocation_respects_bounds():
    from pyomo_cvp import control_value

    m = discretize(racecar())
    pyo.TransformationFactory("cvp.parameterize").apply_to(
        m, var=m.u, contset=m.tau, profile=("collocation", 2)
    )
    r = pyo.SolverFactory("ipopt").solve(m)
    assert r.solver.termination_condition == pyo.TerminationCondition.optimal
    for t in sorted(m.tau):
        u = control_value(m.u, t)
        assert -3 - 1e-6 <= u <= 1 + 1e-6


def test_collocation_guards():
    m = racecar()
    pyo.TransformationFactory("dae.finite_difference").apply_to(
        m, nfe=NFE, scheme="BACKWARD"
    )
    with pytest.raises(RuntimeError, match="requires a collocation"):
        pyo.TransformationFactory("cvp.parameterize").apply_to(
            m, var=m.u, contset=m.tau, profile=("collocation", 2)
        )

    m2 = discretize(racecar())
    with pytest.raises(ValueError, match="exceeds"):
        pyo.TransformationFactory("cvp.parameterize").apply_to(
            m2, var=m2.u, contset=m2.tau, profile=("collocation", NCP + 1)
        )


@needs_ipopt
def test_plain_collocation_resolves_to_ncp():
    m1 = discretize(racecar())
    pyo.TransformationFactory("cvp.parameterize").apply_to(
        m1, var=m1.u, contset=m1.tau, profile="collocation"
    )
    assert len(m1.u) == NFE * NCP
    r1 = pyo.SolverFactory("ipopt").solve(m1)
    assert r1.solver.termination_condition == pyo.TerminationCondition.optimal

    m2 = discretize(racecar())
    pyo.TransformationFactory("cvp.parameterize").apply_to(
        m2, var=m2.u, contset=m2.tau, profile=("collocation", NCP)
    )
    r2 = pyo.SolverFactory("ipopt").solve(m2)
    assert r2.solver.termination_condition == pyo.TerminationCondition.optimal
    assert pyo.value(m1.tf) == pytest.approx(pyo.value(m2.tf), rel=1e-8)


def test_unknown_tuple_profile():
    m = discretize(racecar())
    with pytest.raises(ValueError, match="unknown profile"):
        pyo.TransformationFactory("cvp.parameterize").apply_to(
            m, var=m.u, contset=m.tau, profile=("spline", 3)
        )


def test_control_value_helper():
    from pyomo_cvp import control_value

    m = discretize(racecar())
    pyo.TransformationFactory("cvp.parameterize").apply_to(m, var=m.u, contset=m.tau)
    for t in m.u:
        m.u[t] = 0.25
    assert control_value(m.u, 0.5) == pytest.approx(0.25)

    m2 = discretize(racecar())
    pyo.TransformationFactory("cvp.parameterize").apply_to(
        m2, var=m2.u, contset=m2.tau, profile="piecewise_linear"
    )
    knots = sorted(m2.u)
    m2.u[knots[0]] = 0.0
    m2.u[knots[1]] = 1.0
    mid = (knots[0] + knots[1]) / 2
    assert control_value(m2.u, mid) == pytest.approx(0.5)
