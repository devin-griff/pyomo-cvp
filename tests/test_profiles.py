"""Phase 4: piecewise_linear and reduced_collocation profiles."""
import pytest
import pyomo.environ as pyo

import pyomo_cvp  # noqa: F401

from test_parameterize import racecar, NFE, NCP

ipopt_available = pyo.SolverFactory("ipopt").available(False)
needs_ipopt = pytest.mark.skipif(not ipopt_available, reason="ipopt not available")


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
def test_reduced_collocation_dominates_rcp(scheme):
    # The problem is a nonconvex NLP, so structurally different (but
    # mathematically equivalent) models can land in different local optima.
    # The correct equivalence check is warm-start dominance: started from the
    # rcp optimum, the elimination model must do at least as well.
    k = 2
    m1 = racecar()
    d = pyo.TransformationFactory("dae.collocation")
    d.apply_to(m1, nfe=NFE, ncp=NCP, scheme=scheme)
    d.reduce_collocation_points(m1, var=m1.u, ncp=k, contset=m1.tau)
    r1 = pyo.SolverFactory("ipopt").solve(m1)
    assert r1.solver.termination_condition == pyo.TerminationCondition.optimal

    m2 = discretize(racecar(), scheme=scheme)
    pyo.TransformationFactory("cvp.parameterize").apply_to(
        m2, var=m2.u, contset=m2.tau, profile=("reduced_collocation", k)
    )
    assert len(m2.u) == NFE * k

    # warm start every m2 variable from the rcp solution (skip unset values)
    for t in m2.tau:
        for src, dst in ((m1.x, m2.x), (m1.v, m2.v), (m1.dx, m2.dx),
                         (m1.dv, m2.dv)):
            if src[t].value is not None:
                dst[t] = src[t].value
    for t in m2.u:
        if m1.u[t].value is not None:
            m2.u[t] = m1.u[t].value
    m2.tf = m1.tf.value

    r2 = pyo.SolverFactory("ipopt").solve(m2)
    assert r2.solver.termination_condition == pyo.TerminationCondition.optimal
    assert pyo.value(m2.tf) <= pyo.value(m1.tf) * (1 + 1e-6)


def test_reduced_collocation_guards():
    m = racecar()
    pyo.TransformationFactory("dae.finite_difference").apply_to(
        m, nfe=NFE, scheme="BACKWARD"
    )
    with pytest.raises(RuntimeError, match="requires a collocation"):
        pyo.TransformationFactory("cvp.parameterize").apply_to(
            m, var=m.u, contset=m.tau, profile=("reduced_collocation", 2)
        )

    m2 = discretize(racecar())
    with pytest.raises(ValueError, match="exceeds"):
        pyo.TransformationFactory("cvp.parameterize").apply_to(
            m2, var=m2.u, contset=m2.tau, profile=("reduced_collocation", NCP + 1)
        )


def test_unknown_tuple_profile():
    m = discretize(racecar())
    with pytest.raises(ValueError, match="unknown profile"):
        pyo.TransformationFactory("cvp.parameterize").apply_to(
            m, var=m.u, contset=m.tau, profile=("spline", 3)
        )


def test_control_value_helper():
    from pyomo_cvp import control_value

    m = discretize(racecar())
    pyo.TransformationFactory("cvp.parameterize").apply_to(
        m, var=m.u, contset=m.tau
    )
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
