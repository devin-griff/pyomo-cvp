# Copyright (c) 2026 Devin Griffith
# SPDX-License-Identifier: BSD-3-Clause
"""Shared model builders and markers for the pyomo-cvp test suite.

Imported by the individual test modules so that no test module has to import
another (which would couple collection order and re-run the imported module's
top-level solver probe).
"""
import pytest
import pyomo.environ as pyo
from pyomo.dae import ContinuousSet, DerivativeVar

ipopt_available = pyo.SolverFactory("ipopt").available(False)
needs_ipopt = pytest.mark.skipif(not ipopt_available, reason="ipopt not available")

NFE, NCP = 15, 3
L, R = 100.0, 0.001


def racecar():
    """Build the minimum-time race-car model used across the suite."""
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
