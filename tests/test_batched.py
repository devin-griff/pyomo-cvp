# Copyright (c) 2026 Devin Griffith
# SPDX-License-Identifier: BSD-3-Clause
"""Batched parameterization: one pass equals sequential passes.

Controls sharing a time set parameterize in a single substitution pass
(the list form of the explicit call, and automatic grouping in declaration
mode). The transformed model must be identical to parameterizing the same
controls one call at a time.
"""
import pyomo.environ as pyo
from pyomo.dae import ContinuousSet, DerivativeVar

import pyomo_cvp  # noqa: F401
from pyomo_cvp import declare_profile

N = 4


def build(declare=False):
    # two controls with different profiles, and a cost over the elements
    m = pyo.ConcreteModel()
    m.i = ContinuousSet(initialize=pyo.RangeSet(0, N, 1))
    m.z = pyo.Var(m.i, initialize=0.6)
    m.zdot = DerivativeVar(m.z, wrt=m.i)
    m.u = pyo.Var(m.i, initialize=0.5, bounds=(0, 1))
    m.w = pyo.Var(m.i, initialize=0.25)

    @m.Constraint(m.i)
    def ode(mm, t):
        return mm.zdot[t] == -mm.z[t] + mm.u[t] + mm.w[t]

    m.obj = pyo.Objective(
        expr=sum(
            (m.z[t] - 0.5) ** 2 + m.u[t] ** 2 + m.w[t] ** 2 for t in sorted(m.i)[:-1]
        )
    )
    if declare:
        declare_profile(m.u, wrt=m.i, profile="piecewise_constant")
        declare_profile(m.w, wrt=m.i, profile="piecewise_linear")
    pyo.TransformationFactory("dae.collocation").apply_to(
        m, wrt=m.i, nfe=N, ncp=3, scheme="LAGRANGE-RADAU"
    )
    return m


def snapshot(m):
    lines = [
        f"{c.name}: {c.expr}"
        for c in m.component_data_objects(pyo.Constraint, active=True)
    ]
    lines += [
        f"{o.name}: {o.expr}"
        for o in m.component_data_objects(pyo.Objective, active=True)
    ]
    return lines


def test_the_list_form_equals_sequential_calls():
    a = build()
    pyo.TransformationFactory("cvp.parameterize").apply_to(
        a,
        var=[a.u, a.w],
        contset=a.i,
        profile=["piecewise_constant", "piecewise_linear"],
    )
    b = build()
    pyo.TransformationFactory("cvp.parameterize").apply_to(
        b, var=b.u, contset=b.i, profile="piecewise_constant"
    )
    pyo.TransformationFactory("cvp.parameterize").apply_to(
        b, var=b.w, contset=b.i, profile="piecewise_linear"
    )
    assert snapshot(a) == snapshot(b)


def test_declaration_mode_makes_one_pass(monkeypatch):
    # both declared controls arrive in a single _parameterize call
    calls = []
    cls = type(pyo.TransformationFactory("cvp.parameterize"))
    orig = cls._parameterize

    def counting(self, model, var, contset, profile):
        calls.append(var)
        return orig(self, model, var, contset, profile)

    monkeypatch.setattr(cls, "_parameterize", counting)
    m = build(declare=True)
    pyo.TransformationFactory("cvp.parameterize").apply_to(m)
    assert len(calls) == 1
    assert [v.local_name for v in calls[0]] == ["u", "w"]
