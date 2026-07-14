# Copyright (c) 2026 Devin Griffith
# SPDX-License-Identifier: BSD-3-Clause
"""Phase 0: registration smoke tests."""
import pytest
import pyomo.environ as pyo

import pyomo_cvp  # noqa: F401  (registers the plugin)


def test_transformation_registers():
    xf = pyo.TransformationFactory("cvp.parameterize")
    assert xf is not None


def test_no_args_without_declarations():
    m = pyo.ConcreteModel()
    with pytest.raises(RuntimeError, match="no control profile declarations"):
        pyo.TransformationFactory("cvp.parameterize").apply_to(m)


def test_partial_args_rejected():
    from pyomo.dae import ContinuousSet

    m = pyo.ConcreteModel()
    m.tau = ContinuousSet(bounds=(0, 1))
    m.u = pyo.Var(m.tau)
    with pytest.raises(TypeError, match="both var= and contset="):
        pyo.TransformationFactory("cvp.parameterize").apply_to(m, var=m.u)


def test_declare_profile_records():
    from pyomo.dae import ContinuousSet
    from pyomo_cvp.parameterize import _cvp_data

    m = pyo.ConcreteModel()
    m.tau = ContinuousSet(bounds=(0, 1))
    m.u = pyo.Var(m.tau)
    pyomo_cvp.declare_profile(m.u, wrt=m.tau)
    decls = _cvp_data(m)["profiles"]
    assert len(decls) == 1 and decls[0]["var"] is m.u


def test_unknown_profile_rejected():
    from pyomo.dae import ContinuousSet

    m = pyo.ConcreteModel()
    m.tau = ContinuousSet(bounds=(0, 1))
    m.u = pyo.Var(m.tau)
    with pytest.raises(ValueError, match="unknown profile"):
        pyomo_cvp.declare_profile(m.u, wrt=m.tau, profile="cubic_spline")
