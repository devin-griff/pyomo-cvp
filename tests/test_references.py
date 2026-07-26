# Copyright (c) 2026 Devin Griffith
# SPDX-License-Identifier: BSD-3-Clause
"""References and Ports through the control rebuild (gh #3).

The rebuild replaces the control container. Two Reference shapes interact
with that. A Reference TO a replaced flat control (the IDAES ``heat_duty``
idiom) pointed at deleted members and read pre-solve values; it is now
rebuilt under its own name as a view of the profile, and Ports are
re-pointed at the rebuilt view. A control that IS a Reference (the IDAES
inlet idiom, referents inside ``Block(time)`` members) left its referents
orphaned; each referent is now tied to its profile value by an equality
row, so it stays live for Ports, Arcs, and reporting, and the tie rows
against the otherwise-unconstrained referents keep the transform
dof-neutral.
"""
import pyomo.environ as pyo
import pytest
from pyomo.core.expr.visitor import identify_variables
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.network import Port

import pyomo_cvp  # noqa: F401
from pyomo_cvp import declare_profile

N = 4


def flat_with_alias(profile="piecewise_constant"):
    """A flat control plus a Reference alias and a Port holding it."""
    m = pyo.ConcreteModel()
    m.t = ContinuousSet(initialize=pyo.RangeSet(0, N, 1))
    m.z = pyo.Var(m.t, initialize=0.6)
    m.zdot = DerivativeVar(m.z, wrt=m.t)
    m.u = pyo.Var(m.t, initialize=0.5)
    m.u_alias = pyo.Reference(m.u[:])
    m.port = Port()
    m.port.add(m.u_alias, "u")
    declare_profile(m.u, wrt=m.t, profile=profile)

    @m.Constraint(m.t)
    def ode(mm, t):
        return mm.zdot[t] == mm.u[t] - mm.z[t]

    pyo.TransformationFactory("dae.collocation").apply_to(
        m, wrt=m.t, nfe=N, ncp=2, scheme="LAGRANGE-RADAU"
    )
    return m


def member_control():
    """The IDAES inlet idiom: the control is a Reference into Block members."""
    m = pyo.ConcreteModel()
    m.t = ContinuousSet(initialize=pyo.RangeSet(0, N, 1))
    m.props = pyo.Block(m.t, rule=lambda b, t: setattr(b, "f", pyo.Var(initialize=1.0)))
    m.fin = pyo.Reference(m.props[:].f)
    m.port = Port()
    m.port.add(m.fin, "f")
    m.z = pyo.Var(m.t, initialize=0.6)
    m.zdot = DerivativeVar(m.z, wrt=m.t)
    declare_profile(m.fin, wrt=m.t)

    @m.Constraint(m.t)
    def ode(mm, t):
        return mm.zdot[t] == mm.props[t].f - mm.z[t]

    pyo.TransformationFactory("dae.collocation").apply_to(
        m, wrt=m.t, nfe=N, ncp=2, scheme="LAGRANGE-RADAU"
    )
    return m


def _dof(m):
    from pyomo.util.model_size import build_model_size_report

    r = build_model_size_report(m)
    return r.activated.variables - r.activated.constraints


def test_reference_to_a_control_is_rebuilt_as_a_view():
    m = flat_with_alias()
    pyo.TransformationFactory("cvp.parameterize").apply_to(m)
    assert m.u_alias.is_reference()
    # piecewise_constant: every entry follows its interval's free value, so
    # the view keeps the full index and reads the applied profile
    m.u[0].set_value(0.7)
    t_in = next(t for t in sorted(m.t) if 0 < t < 1)
    assert m.u_alias[0] is m.u[0]
    assert m.u_alias[t_in] is m.u[0]  # interior point: the element's value
    assert pyo.value(m.u_alias[t_in]) == 0.7


def test_port_is_repointed_at_the_rebuilt_view():
    m = flat_with_alias()
    pyo.TransformationFactory("cvp.parameterize").apply_to(m)
    assert m.port.vars["u"] is m.u_alias
    assert m.port.vars["u"].parent_block() is m


def test_interpolated_entries_are_dropped_from_the_view():
    m = flat_with_alias(profile="piecewise_linear")
    pyo.TransformationFactory("cvp.parameterize").apply_to(m)
    # boundaries stay free and keep their view; interior collocation points
    # substitute to interpolations, which no Reference member can hold
    kept = sorted(m.u_alias.keys())
    assert kept == sorted(m.u.keys())
    assert m.u_alias[1] is m.u[1]


def test_reference_control_referents_are_tied_and_live():
    m = member_control()
    before = _dof(m)
    pyo.TransformationFactory("cvp.parameterize").apply_to(m)
    # the ties absorb the referents one-for-one, so the reference route
    # loses exactly what the flat route loses: the eliminated copies
    assert _dof(m) == before - (len(sorted(m.t)) - N)
    ties = m.component("fin_profile_ties")
    assert ties is not None and len(ties) == len(sorted(m.t))
    # every referent appears in a tie row, bound to the move variables
    tied = set()
    for row in ties.values():
        tied.update(id(v) for v in identify_variables(row.body))
    assert all(id(m.props[t].f) in tied for t in m.t)


def test_reference_control_port_reads_live_referents():
    m = member_control()
    pyo.TransformationFactory("cvp.parameterize").apply_to(m)
    # the port's entry views the referents, which the ties keep live: the
    # data it yields are exactly the member variables
    assert m.port.vars["f"][0] is m.props[0].f
    live = set(
        id(v)
        for c in m.component_data_objects(pyo.Constraint, active=True)
        for v in identify_variables(c.body)
    )
    assert all(id(m.props[t].f) in live for t in m.t)


def test_flat_model_without_references_is_unchanged():
    m = flat_with_alias()
    m2 = flat_with_alias()
    m2.del_component(m2.port)
    m2.del_component(m2.u_alias)
    pyo.TransformationFactory("cvp.parameterize").apply_to(m)
    pyo.TransformationFactory("cvp.parameterize").apply_to(m2)
    c1 = sorted(
        str(c.expr) for c in m.component_data_objects(pyo.Constraint, active=True)
    )
    c2 = sorted(
        str(c.expr) for c in m2.component_data_objects(pyo.Constraint, active=True)
    )
    assert c1 == c2
