"""The ``cvp.parameterize`` transformation and the ``declare_profile`` API.

Reduces a control variable to a small set of free values over a discretized
ContinuousSet, by substitution: the component is replaced under its own name
by a Var indexed by the profile's free time points, and every other copy is
substituted out of all constraints, objectives, and named expressions as a
linear combination of the free values. No linking constraints, no leftover
copies.

Profiles
--------
``'piecewise_constant'``
    One free value per finite element, at the element's last point (equals
    the point ``reduce_collocation_points`` keeps free under Lagrange-Radau).
    Element i owns the half-open interval (fe[i], fe[i+1]]; the initial point
    belongs to the first element.
``'piecewise_linear'``
    One free value per element boundary (nfe + 1 values), continuous, with
    interior points substituted by linear interpolation between the two
    surrounding boundary values.
``('reduced_collocation', k)``
    k free values per element, at the element's last k collocation points,
    with the remaining points substituted by Lagrange interpolation --- the
    elimination form of ``reduce_collocation_points(ncp=k)``. Requires a
    collocation discretization with k <= ncp.
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

PROFILES = ("piecewise_constant", "piecewise_linear")

_DECLARATION_ATTR = "_pyomo_cvp_profiles"


def declare_profile(var, wrt, profile="piecewise_constant"):
    """Declare a control profile for ``var`` over ContinuousSet ``wrt``.

    Inert metadata recorded on the variable's parent block; applied later by
    ``TransformationFactory('cvp.parameterize').apply_to(model)``.
    """
    _validate_profile(profile)
    block = var.parent_block()
    decls = getattr(block, _DECLARATION_ATTR, None)
    if decls is None:
        decls = []
        setattr(block, _DECLARATION_ATTR, decls)
    decls.append({"var": var, "wrt": wrt, "profile": profile})


def _validate_profile(profile):
    if isinstance(profile, (tuple, list)):
        if (
            len(profile) == 2
            and profile[0] == "reduced_collocation"
            and int(profile[1]) >= 1
        ):
            return "reduced_collocation", int(profile[1])
        raise ValueError(
            f"pyomo-cvp: unknown profile {profile!r}; tuple profiles must be "
            f"('reduced_collocation', k) with k >= 1."
        )
    if profile not in PROFILES:
        raise ValueError(
            f"pyomo-cvp: unknown profile {profile!r}; supported: "
            f"{PROFILES + (('reduced_collocation', 'k'),)}"
        )
    return profile, None


def control_value(var, t, index=()):
    """Evaluate a parameterized control at any time ``t``.

    ``var`` is the (replaced) control component; ``index`` supplies any
    non-time index components, in their original order. Works for every
    profile: piecewise-constant lookup, piecewise-linear interpolation, or
    the element's Lagrange polynomial for reduced collocation.
    """
    from pyomo.core import value as _value

    info_map = getattr(var.parent_block(), "_pyomo_cvp_info", {})
    if var.local_name not in info_map:
        raise ValueError(
            f"pyomo-cvp: '{var.name}' has not been parameterized."
        )
    info = info_map[var.local_name]
    fe, mode, pairs = info["fe"], info["mode"], info["pairs"]
    if t <= fe[0]:
        elem = 0
    else:
        elem = bisect_left(fe, min(t, fe[-1])) - 1

    if mode == "piecewise_constant":
        plist = [(fe[elem + 1], 1.0)]
    elif mode == "piecewise_linear":
        a, b = fe[elem], fe[elem + 1]
        w = (t - a) / (b - a)
        plist = [(a, 1.0 - w), (b, w)]
    else:  # reduced_collocation: the element's free knots
        knots = sorted(
            {ft for tt, pl in pairs.items()
             for ft, _ in pl
             if fe[elem] < tt <= fe[elem + 1] or (elem == 0 and tt <= fe[0])}
        )
        knots = [ft for ft in knots if fe[elem] < ft <= fe[elem + 1]]
        plist = list(zip(knots, _lagrange_coeffs(t, knots)))

    pos = info["pos"]
    index = tuple(index)

    def member(ft):
        full = index[:pos] + (ft,) + index[pos:]
        return var[full if len(full) > 1 else full[0]]

    return sum(c * _value(member(ft)) for ft, c in plist)


def _lagrange_coeffs(t, knots):
    coeffs = []
    for i in knots:
        c = 1.0
        for j in knots:
            if i != j:
                c *= (t - j) / (i - j)
        coeffs.append(c)
    return coeffs


@TransformationFactory.register(
    "cvp.parameterize",
    doc="Parameterize controls over finite elements by variable elimination "
    "(pyomo-cvp).",
)
class ParameterizeTransformation(Transformation):
    """Reduce a control to a profile's free values, by substitution."""

    def _apply_to(self, model, var=None, contset=None,
                  profile="piecewise_constant", **kwds):
        if var is None and contset is None:
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
        mode, k = _validate_profile(profile)
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
        disc_info = contset.get_discretization_info()
        if not disc_info:
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
                    f"a parameterized control cannot be differentiated."
                )

        fe = list(contset.get_finite_elements())
        points = sorted(contset)
        pairs = self._profile_pairs(mode, k, fe, points, disc_info, contset)
        free_times = sorted({ft for plist in pairs.values() for ft, _ in plist})

        def as_tuple(i):
            return i if isinstance(i, tuple) else (i,)

        def with_t(full, t):
            return full[:pos] + (t,) + full[pos + 1:]

        old = {as_tuple(i): var[i] for i in var}
        free_attrs = {}
        for full, vd in old.items():
            if full[pos] in free_times:
                free_attrs[full] = (vd.domain, vd.lb, vd.ub, vd.value)

        name = var.local_name
        parent = var.parent_block()
        parent.del_component(var)
        newsets = subsets[:pos] + [free_times] + subsets[pos + 1:]
        any_dom = next(iter(free_attrs.values()))[0]
        new_var = Var(
            *newsets,
            domain=any_dom,
            bounds=lambda m, *i: (free_attrs[i][1], free_attrs[i][2]),
            initialize=lambda m, *i: free_attrs[i][3],
        )
        parent.add_component(name, new_var)

        submap = {}
        for full, vd in old.items():
            plist = pairs[full[pos]]
            if len(plist) == 1 and plist[0][1] == 1.0:
                node = new_var[with_t(full, plist[0][0])]
            else:
                node = sum(
                    c * new_var[with_t(full, ft)] for ft, c in plist
                )
            submap[id(vd)] = node

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

        info = getattr(parent, "_pyomo_cvp_info", None)
        if info is None:
            info = {}
            setattr(parent, "_pyomo_cvp_info", info)
        info[name] = {"mode": mode, "k": k, "fe": fe, "pos": pos,
                      "pairs": pairs}

        return model

    @staticmethod
    def _profile_pairs(mode, k, fe, points, disc_info, contset):
        """Map each time point to [(free_time, coeff), ...]."""
        pairs = {}

        def element_of(t):
            if t <= fe[0]:
                return 0
            return bisect_left(fe, t) - 1

        if mode == "piecewise_constant":
            for t in points:
                pairs[t] = [(fe[element_of(t) + 1], 1.0)]
            return pairs

        if mode == "piecewise_linear":
            for t in points:
                if t in fe:
                    pairs[t] = [(t, 1.0)]
                else:
                    i = element_of(t)
                    a, b = fe[i], fe[i + 1]
                    w = (t - a) / (b - a)
                    pairs[t] = [(a, 1.0 - w), (b, w)]
            return pairs

        # reduced_collocation
        ncp = disc_info.get("ncp")
        if ncp is None:
            raise RuntimeError(
                "pyomo-cvp: ('reduced_collocation', k) requires a collocation "
                "discretization (dae.collocation), not finite differences."
            )
        if k > ncp:
            raise ValueError(
                f"pyomo-cvp: reduced_collocation k={k} exceeds the "
                f"discretization's ncp={ncp}."
            )
        scheme = disc_info.get("scheme", "")
        legendre = "LEGENDRE" in str(scheme).upper()
        # collocation points of element i: the points strictly inside
        # (fe[i], fe[i+1]) plus, for Radau, the right boundary
        elem_pts = []
        for i in range(len(fe) - 1):
            inside = [t for t in points if fe[i] < t < fe[i + 1]]
            if not legendre:
                inside = inside + [fe[i + 1]]
            elem_pts.append(inside)
        free = [pts[-k:] for pts in elem_pts]
        for t in points:
            i = element_of(t)
            knots = free[i]
            if t in knots:
                pairs[t] = [(t, 1.0)]
            else:
                pairs[t] = list(zip(knots, _lagrange_coeffs(t, knots)))
        return pairs
