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
        if var is None and contset is None:
            # discovery mode: apply every declare_profile() declaration
            found = 0
            for block in model.block_data_objects(active=True,
                                                  descend_into=True):
                decls = getattr(block, _DECLARATION_ATTR, None)
                if not decls:
                    continue
                while decls:
                    d = decls.pop(0)
                    self._parameterize(model, d["var"], d["wrt"], d["profile"])
                    found += 1
            if found == 0:
                raise RuntimeError(
                    "pyomo-cvp: no control profile declarations found on the "
                    "model (either none were made with declare_profile(), or "
                    "they were already applied)."
                )
            return model
        return self._parameterize(model, var, contset, profile)

    def _parameterize(self, model, var, contset, profile):
        if profile not in PROFILES:
            raise ValueError(
                f"pyomo-cvp: unknown profile {profile!r}; supported: {PROFILES}"
            )
        if var is None or contset is None:
            raise TypeError(
                "pyomo-cvp: pass both var= and contset=, or neither (to apply "
                "declare_profile declarations)."
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
        idxset = var.index_set()
        if idxset is contset:
            subsets, pos = [contset], 0
        else:
            subsets = list(idxset.subsets())
            matches = [i for i, s in enumerate(subsets) if s is contset]
            if len(matches) != 1:
                raise ValueError(
                    f"pyomo-cvp: '{var.name}' is not indexed by "
                    f"'{contset.name}' (exactly once). Either it was already "
                    f"parameterized, or it is not a control over this set."
                )
            pos = matches[0]
        # DerivativeVar registers with ctype Var, so filter by isinstance
        for dv in model.component_objects(Var, active=True):
            if isinstance(dv, DerivativeVar) and dv.get_state_var() is var:
                raise ValueError(
                    f"pyomo-cvp: '{var.name}' has DerivativeVar '{dv.name}'; "
                    f"a piecewise-constant control cannot be differentiated."
                )

        fe = list(contset.get_finite_elements())
        reps = fe[1:]                       # representative of element i: fe[i+1]

        def rep_of(t):
            if t <= fe[0]:
                return fe[1]
            return fe[bisect_left(fe, t)]

        def as_tuple(i):
            return i if isinstance(i, tuple) else (i,)

        def to_rep(full):
            full = as_tuple(full)
            return full[:pos] + (rep_of(full[pos]),) + full[pos + 1:]

        # capture the old copies and, for representative indices, attributes
        old = {as_tuple(i): var[i] for i in var}
        rep_attrs = {}
        for full, vd in old.items():
            if full[pos] in reps:
                rep_attrs[full] = (vd.domain, vd.lb, vd.ub, vd.value)

        name = var.local_name
        parent = var.parent_block()
        parent.del_component(var)
        newsets = subsets[:pos] + [reps] + subsets[pos + 1:]
        any_dom = next(iter(rep_attrs.values()))[0]
        new_var = Var(
            *newsets,
            domain=any_dom,
            bounds=lambda m, *i: (rep_attrs[i][1], rep_attrs[i][2]),
            initialize=lambda m, *i: rep_attrs[i][3],
        )
        parent.add_component(name, new_var)

        submap = {
            id(vd): new_var[to_rep(full)] for full, vd in old.items()
        }

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
