"""The ``cvp.parameterize`` transformation and the ``declare_profile`` API.

Reduces a control variable to one free value per finite element of a
discretized ContinuousSet, by substitution: the component is replaced under
its own name by a Var indexed by the element representative points, and every
other copy is substituted out of all constraints, objectives, and named
expressions. No linking constraints, no leftover fixed copies.

Ownership rule: finite element i owns the half-open interval
(fe[i], fe[i+1]]; the representative point is fe[i+1] (the element's last
point). The initial point fe[0] belongs to the first element. For
Lagrange-Radau collocation the representative coincides with the free point
kept by ``reduce_collocation_points``.
"""
from bisect import bisect_left

from pyomo.core import (
    Constraint,
    Expression,
    Objective,
    Transformation,
    TransformationFactory,
    Var,
)
from pyomo.core.expr import replace_expressions
from pyomo.dae import ContinuousSet, DerivativeVar

PROFILES = ("piecewise_constant",)

_DECLARATION_ATTR = "_pyomo_cvp_profiles"


def declare_profile(var, wrt, profile="piecewise_constant"):
    """Declare a control profile for ``var`` over ContinuousSet ``wrt``.

    Inert metadata recorded on the variable's parent block; applied later by
    ``TransformationFactory('cvp.parameterize').apply_to(model)``.
    """
    if profile not in PROFILES:
        raise ValueError(
            f"pyomo-cvp: unknown profile {profile!r}; supported: {PROFILES}"
        )
    block = var.parent_block()
    decls = getattr(block, _DECLARATION_ATTR, None)
    if decls is None:
        decls = []
        setattr(block, _DECLARATION_ATTR, decls)
    decls.append({"var": var, "wrt": wrt, "profile": profile})


@TransformationFactory.register(
    "cvp.parameterize",
    doc="Parameterize controls over finite elements by variable elimination "
    "(pyomo-cvp).",
)
class ParameterizeTransformation(Transformation):
    """Reduce a control to one free value per finite element, by substitution."""

    def _apply_to(self, model, var=None, contset=None,
                  profile="piecewise_constant", **kwds):
        if profile not in PROFILES:
            raise ValueError(
                f"pyomo-cvp: unknown profile {profile!r}; supported: {PROFILES}"
            )
        if var is None or contset is None:
            raise TypeError(
                "pyomo-cvp: pass var= and contset= explicitly (discovery of "
                "declare_profile declarations arrives in Phase 3)."
            )
        if contset.ctype is not ContinuousSet:
            raise TypeError(
                f"pyomo-cvp: contset must be a ContinuousSet, got "
                f"{contset.ctype.__name__}."
            )
        if not contset.get_discretization_info():
            raise RuntimeError(
                f"pyomo-cvp: ContinuousSet '{contset.name}' has not been "
                f"discretized; apply a dae.* transformation first."
            )
        if var.ctype is not Var:
            raise TypeError("pyomo-cvp: var must be a Var.")
        if var.index_set() is not contset:
            raise ValueError(
                f"pyomo-cvp: '{var.name}' is not indexed by (exactly) "
                f"'{contset.name}'. Either it was already parameterized, or "
                f"it has additional index sets (Phase 2)."
            )
        # DerivativeVar registers with ctype Var, so filter by isinstance
        for dv in model.component_objects(Var, active=True):
            if isinstance(dv, DerivativeVar) and dv.get_state_var() is var:
                raise ValueError(
                    f"pyomo-cvp: '{var.name}' has DerivativeVar '{dv.name}'; "
                    f"a piecewise-constant control cannot be differentiated."
                )

        fe = list(contset.get_finite_elements())
        points = sorted(contset)
        reps = fe[1:]                       # representative of element i: fe[i+1]

        def rep_of(t):
            if t <= fe[0]:
                return fe[1]
            return fe[bisect_left(fe, t)]

        # capture the old copies and their attributes
        old = {t: var[t] for t in points}
        rep_attrs = {
            r: (old[r].domain, old[r].lb, old[r].ub, old[r].value)
            for r in reps
        }

        name = var.local_name
        parent = var.parent_block()
        parent.del_component(var)
        new_var = Var(
            reps,
            domain=rep_attrs[reps[0]][0],
            bounds=lambda m, t: (rep_attrs[t][1], rep_attrs[t][2]),
            initialize=lambda m, t: rep_attrs[t][3],
        )
        parent.add_component(name, new_var)

        submap = {id(old[t]): new_var[rep_of(t)] for t in points}

        for c in model.component_data_objects(
            Constraint, active=True, descend_into=True
        ):
            c.set_value(replace_expressions(c.expr, submap))
        for o in model.component_data_objects(
            Objective, active=True, descend_into=True
        ):
            o.set_value(replace_expressions(o.expr, submap))
        for e in model.component_data_objects(
            Expression, active=True, descend_into=True
        ):
            e.set_value(replace_expressions(e.expr, submap))

        return model
