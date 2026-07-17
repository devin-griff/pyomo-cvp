# Copyright (c) 2026 Devin Griffith
# SPDX-License-Identifier: BSD-3-Clause
"""Algebraic vs differential control references (piecewise_constant shift).

The bug: parameterize rewrote control references uniformly, so an algebraic
reference to element k's control u[fe[k]] became element k-1's u[fe[k-1]],
silently penalizing the wrong element in objectives and path constraints.

The fix splits the rewrite by expression kind. Differential equations (those
carrying a DerivativeVar: the user's ODEs and pyomo.dae's disc_eq) keep the
old map, which is correct for the collocation node there. Algebraic
expressions resolve a control reference to the element the reference sits in.
The final node starts no element, so an algebraic reference to it errors
instead of double-counting the last element. Only the discontinuous
piecewise_constant profile is affected; piecewise_linear and
collocation are unchanged (their maps already coincide).
"""
import re

import pytest
import pyomo.environ as pyo
from pyomo.dae import ContinuousSet, DerivativeVar

import pyomo_cvp  # noqa: F401
from pyomo_cvp import declare_profile

N = 4


def build(profile, obj="elements", path=False):
    m = pyo.ConcreteModel()
    m.i = ContinuousSet(initialize=pyo.RangeSet(0, N, 1))
    m.z = pyo.Var(m.i, initialize=0.6)
    m.zdot = DerivativeVar(m.z, wrt=m.i)
    m.u = pyo.Var(m.i, initialize=0.5)
    declare_profile(m.u, wrt=m.i, profile=profile)

    @m.Constraint(m.i)
    def ode(mm, t):
        return mm.zdot[t] == mm.u[t] - mm.z[t]

    pts = sorted(m.i)
    idxs = pts[:-1] if obj == "elements" else pts  # element starts vs all
    m.obj = pyo.Objective(expr=sum(m.u[i] for i in idxs))
    if path:

        @m.Constraint(pts[:-1])
        def cap(mm, t):
            return mm.u[t] <= 1.0

    return m


def discretize(m, scheme="LAGRANGE-RADAU"):
    if scheme == "FD":
        pyo.TransformationFactory("dae.finite_difference").apply_to(
            m, nfe=N, scheme="BACKWARD"
        )
    else:
        pyo.TransformationFactory("dae.collocation").apply_to(
            m, nfe=N, ncp=3, scheme=scheme
        )
    return m


def urefs(expr):
    return re.findall(r"u\[[0-9.]+\]", str(expr))


@pytest.mark.parametrize("scheme", ["LAGRANGE-RADAU", "LAGRANGE-LEGENDRE", "FD"])
def test_pwc_algebraic_no_shift(scheme):
    # the correctness fix: an objective over the element starts references each
    # element's control once, at its own index -- no shift, no double-count.
    m = discretize(build("piecewise_constant", obj="elements"), scheme)
    pyo.TransformationFactory("cvp.parameterize").apply_to(m)
    assert urefs(m.obj.expr) == ["u[0]", "u[1]", "u[2]", "u[3]"]


@pytest.mark.parametrize("scheme", ["LAGRANGE-RADAU", "LAGRANGE-LEGENDRE", "FD"])
def test_pwc_ode_unchanged(scheme):
    # differential regression: the ODE at boundary fe[k] uses element k-1's
    # control, and the final node uses the last element's -- unchanged.
    m = discretize(build("piecewise_constant", obj="elements"), scheme)
    pyo.TransformationFactory("cvp.parameterize").apply_to(m)
    assert urefs(m.ode[1].body) == ["u[0]"]
    assert urefs(m.ode[2].body) == ["u[1]"]
    assert urefs(m.ode[N].body) == ["u[3]"]


def test_pwc_path_constraint_no_shift():
    # the fix applies to any algebraic constraint, not only the objective.
    m = discretize(build("piecewise_constant", obj="elements", path=True))
    pyo.TransformationFactory("cvp.parameterize").apply_to(m)
    for k in range(N):
        assert urefs(m.cap[k].body) == [f"u[{k}]"]


@pytest.mark.parametrize("scheme", ["LAGRANGE-RADAU", "LAGRANGE-LEGENDRE", "FD"])
def test_pwc_final_node_reference_errors(scheme):
    # the terminal-horizon default: no move exists at the final node, so an
    # algebraic reference to the control there is a modeling error. The
    # convention comes from the profile, not the scheme.
    m = discretize(build("piecewise_constant", obj="all"), scheme)
    with pytest.raises(ValueError, match="final time"):
        pyo.TransformationFactory("cvp.parameterize").apply_to(m)


@pytest.mark.parametrize("scheme", ["LAGRANGE-RADAU", "LAGRANGE-LEGENDRE", "FD"])
def test_pwc_final_node_keep_resolves_to_the_held_move(scheme):
    # the horizon-continues convention, declared per profile: the control at
    # the final node is the held last move, so a DAE's algebraic equations
    # there (the last element's equations) resolve like the collocation
    # equation beside them.
    m = build("piecewise_constant", obj="all")
    declare_profile(m.u, wrt=m.i, profile="piecewise_constant", final_node="keep")
    m = discretize(m, scheme)
    pyo.TransformationFactory("cvp.parameterize").apply_to(m)
    assert urefs(m.obj.expr)[-1] == f"u[{N - 1}]"  # the last element start


def test_final_node_keep_via_the_explicit_form():
    m = discretize(build("piecewise_constant", obj="all"))
    pyo.TransformationFactory("cvp.parameterize").apply_to(
        m, var=m.u, contset=m.i, profile="piecewise_constant", final_node="keep"
    )
    assert urefs(m.obj.expr)[-1] == f"u[{N - 1}]"


def test_redeclaration_replaces_the_pending_declaration():
    # last declaration wins: the sanctioned way to update a declaration's
    # convention before it is applied.
    m = build("piecewise_constant", obj="all")
    declare_profile(m.u, wrt=m.i, profile="piecewise_constant", final_node="keep")
    m = discretize(m)
    pyo.TransformationFactory("cvp.parameterize").apply_to(m)  # one pass only
    assert urefs(m.obj.expr)[-1] == f"u[{N - 1}]"


def test_final_node_call_option_errors_in_declaration_mode():
    m = discretize(build("piecewise_constant", obj="elements"))
    with pytest.raises(ValueError, match="declaration mode"):
        pyo.TransformationFactory("cvp.parameterize").apply_to(m, final_node="keep")


def test_final_node_is_validated():
    m = build("piecewise_constant")
    with pytest.raises(ValueError, match="or 'keep'"):
        declare_profile(m.u, wrt=m.i, profile="piecewise_constant", final_node="drop")


@pytest.mark.parametrize("scheme", ["LAGRANGE-RADAU", "LAGRANGE-LEGENDRE", "FD"])
def test_pwl_algebraic_unchanged(scheme):
    # piecewise_linear is continuous: no shift, and u[fe[N]] is a real knot the
    # last segment interpolates toward, so referencing it must NOT error.
    m = discretize(build("piecewise_linear", obj="all"), scheme)
    pyo.TransformationFactory("cvp.parameterize").apply_to(m)
    assert urefs(m.obj.expr) == ["u[0]", "u[1]", "u[2]", "u[3]", "u[4]"]


@pytest.mark.parametrize("scheme", ["LAGRANGE-RADAU", "LAGRANGE-LEGENDRE"])
def test_rcp_algebraic_unchanged(scheme):
    # collocation: the control is a per-element polynomial, so an
    # algebraic reference resolves to the polynomial value and does not error.
    m = discretize(build(("collocation", 2), obj="elements"), scheme)
    pyo.TransformationFactory("cvp.parameterize").apply_to(m)  # must not raise
    assert urefs(m.obj.expr)
