# Copyright (c) 2026 Devin Griffith
# SPDX-License-Identifier: BSD-3-Clause
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
    One free value per finite element, indexed by the element's START time
    (u[t0] exists; the final time carries no control). Element i owns the
    half-open interval (fe[i], fe[i+1]], so the profile is left-continuous
    at interior boundaries; the initial point belongs to the first element.
``'piecewise_linear'``
    One free value per element boundary (nfe + 1 values), continuous, with
    interior points substituted by linear interpolation between the two
    surrounding boundary values.
``'collocation'``
    The control is the element's collocation polynomial: one free value per
    collocation point (k resolves to the discretization's ncp), with the
    element-boundary copies substituted by the polynomial's value there.
``('collocation', k)``
    k free values per element, at the element's last k collocation points,
    with the remaining points substituted by Lagrange interpolation: the
    elimination form of ``reduce_collocation_points(ncp=k)``. Requires a
    collocation discretization with k <= ncp.
"""
from bisect import bisect_left, bisect_right

from pyomo.common.config import ConfigDict, ConfigValue
from pyomo.core import (
    Constraint,
    ConstraintList,
    Expression,
    Objective,
    Transformation,
    TransformationFactory,
    Var,
    value as _value,
)
from pyomo.core.expr import replace_expressions
from pyomo.core.expr.visitor import identify_variables
from pyomo.dae import ContinuousSet, DerivativeVar

PROFILES = ("piecewise_constant", "piecewise_linear", "collocation")

#: Scope key under which pyomo-cvp stashes its per-block state through
#: :meth:`Block.private_data`: the profile declarations before the transform
#: runs, and the profile metadata that ``control_value`` reads afterward.
_CVP_SCOPE = "pyomo_cvp"


def _cvp_data(block, create=False):
    """Return ``block``'s pyomo-cvp private-data namespace.

    Parameters
    ----------
    block : BlockData
        Block whose namespace is requested.
    create : bool, optional
        When True, allocate the namespace via :meth:`Block.private_data`.
        When False (default), peek without side effect: return ``None`` if
        the block holds no pyomo-cvp data, so scanning many blocks does not
        stamp an empty namespace onto each.

    Returns
    -------
    dict or None
        The namespace, or ``None`` when ``create`` is False and none exists.
    """
    if create:
        return block.private_data(_CVP_SCOPE)
    store = getattr(block, "_private_data", None)
    return store.get(_CVP_SCOPE) if store else None


def declare_profile(*variables, wrt=None, profile="piecewise_constant"):
    """Declare a control-vector-parameterization profile for one or more Vars.

    The declaration is inert metadata recorded on each variable's parent
    block; it is applied later by
    ``TransformationFactory('cvp.parameterize').apply_to(model)``. The same
    profile settings apply to every variable in the call.

    Parameters
    ----------
    *variables : Var
        One or more control Vars, each indexed (once) by the ContinuousSet
        ``wrt``. Example: ``declare_profile(m.v1, m.v2, wrt=m.t)``.
    wrt : ContinuousSet
        The time set the controls are parameterized over. Keyword-only.
    profile : str or tuple, optional
        ``'piecewise_constant'`` (default), ``'piecewise_linear'``,
        ``'collocation'``, or ``('collocation', k)``. See the module
        docstring.

    Raises
    ------
    TypeError
        If no control Var is given, if ``wrt`` is missing or passed
        positionally, or if a ContinuousSet is passed as a control.
    ValueError
        If ``profile`` is not a recognized profile.
    """
    if not variables:
        raise TypeError(
            "pyomo-cvp: declare_profile needs at least one control Var: "
            "declare_profile(m.u, wrt=m.tau)"
        )
    if wrt is None:
        if len(variables) >= 2 and isinstance(variables[1], ContinuousSet):
            raise TypeError(
                "pyomo-cvp: wrt is keyword-only; write "
                "declare_profile(m.u, wrt=m.tau) instead of "
                "declare_profile(m.u, m.tau)"
            )
        raise TypeError("pyomo-cvp: wrt is required: declare_profile(m.u, wrt=m.tau)")
    _validate_profile(profile)
    for var in variables:
        if isinstance(var, ContinuousSet):
            raise TypeError(
                f"pyomo-cvp: {var.name} is a ContinuousSet, not a control "
                "Var; wrt is keyword-only: declare_profile(m.u, wrt=m.tau)"
            )
        block = var.parent_block()
        store = _cvp_data(block, create=True)
        store.setdefault("profiles", []).append(
            {"var": var, "wrt": wrt, "profile": profile}
        )


def _validate_profile(profile):
    """Normalize and validate a profile specifier.

    Parameters
    ----------
    profile : str or tuple
        A profile name or a ``('collocation', k)`` tuple.

    Returns
    -------
    tuple
        ``(mode, k)``: the profile name and the collocation order. ``k`` is
        ``None`` for the piecewise profiles, and for plain ``'collocation'``,
        where it resolves to the discretization's ncp at transform time.

    Raises
    ------
    ValueError
        If ``profile`` is not a recognized profile.
    """
    if isinstance(profile, (tuple, list)):
        if len(profile) == 2 and profile[0] == "collocation" and int(profile[1]) >= 1:
            return "collocation", int(profile[1])
        raise ValueError(
            f"pyomo-cvp: unknown profile {profile!r}; tuple profiles must be "
            f"('collocation', k) with k >= 1."
        )
    if profile not in PROFILES:
        raise ValueError(
            f"pyomo-cvp: unknown profile {profile!r}; supported: "
            f"{PROFILES + (('collocation', 'k'),)}"
        )
    return profile, None


def control_value(var, t, index=()):
    """Evaluate a parameterized control at any time ``t``.

    Works for every profile: piecewise-constant lookup, piecewise-linear
    interpolation, or the element's Lagrange polynomial for reduced
    collocation.

    Parameters
    ----------
    var : Var
        The (already parameterized) control component.
    t : float
        A time in the control's domain ``[fe[0], fe[-1]]``.
    index : tuple, optional
        Any non-time index members, in their original order. Empty for a
        control indexed only by time.

    Returns
    -------
    float
        The control value interpolated at ``t``.

    Raises
    ------
    ValueError
        If ``var`` has not been parameterized, or ``t`` is out of domain.
    """
    store = _cvp_data(var.parent_block())
    info_map = (store or {}).get("info", {})
    if var.local_name not in info_map:
        raise ValueError(f"pyomo-cvp: '{var.name}' has not been parameterized.")
    info = info_map[var.local_name]
    fe, mode, pairs = info["fe"], info["mode"], info["pairs"]
    if not fe[0] <= t <= fe[-1]:
        raise ValueError(
            f"pyomo-cvp: t={t} is outside the control's domain " f"[{fe[0]}, {fe[-1]}]."
        )
    if t <= fe[0]:
        elem = 0
    else:
        elem = bisect_left(fe, t) - 1

    if mode == "piecewise_constant":
        plist = [(fe[elem], 1.0)]
    elif mode == "piecewise_linear":
        a, b = fe[elem], fe[elem + 1]
        w = (t - a) / (b - a)
        plist = [(a, 1.0 - w), (b, w)]
    else:  # collocation: the element's free knots
        knots = sorted(
            {
                ft
                for tt, pl in pairs.items()
                for ft, _ in pl
                if fe[elem] < tt <= fe[elem + 1] or (elem == 0 and tt <= fe[0])
            }
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
    """Lagrange interpolation weights at ``t`` for the given ``knots``.

    Parameters
    ----------
    t : float
        The evaluation point.
    knots : sequence of float
        The interpolation nodes.

    Returns
    -------
    list of float
        One weight per knot, in knot order; ``sum(w_i * f(knot_i))`` is the
        interpolated value.
    """
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
    """Reduce a control to a profile's free values, by substitution.

    Invoke it two ways. With no options it applies every profile declared
    with :func:`declare_profile`. With ``var`` and ``contset`` it
    parameterizes a single control directly::

        TransformationFactory('cvp.parameterize').apply_to(model)
        TransformationFactory('cvp.parameterize').apply_to(
            model, var=m.u, contset=m.t, profile='piecewise_linear')
    """

    CONFIG = ConfigDict("cvp.parameterize")
    CONFIG.declare(
        "var",
        ConfigValue(
            default=None,
            description="Control Var to parameterize (paired with contset). "
            "Omit both var and contset to apply declare_profile declarations.",
        ),
    )
    CONFIG.declare(
        "contset",
        ConfigValue(
            default=None, description="ContinuousSet the control is indexed over."
        ),
    )
    CONFIG.declare(
        "profile",
        ConfigValue(
            default="piecewise_constant",
            description="'piecewise_constant', 'piecewise_linear', "
            "'collocation', or ('collocation', k). Ignored in declaration "
            "mode.",
        ),
    )

    def _apply_to(self, model, **kwds):
        """Apply the parameterization in place; see the class docstring.

        Parameters
        ----------
        model : Block
            The already-discretized model to transform.
        **kwds
            ``var``, ``contset``, and ``profile`` (see ``CONFIG``). Unknown
            options raise ``ValueError``.
        """
        config = self.CONFIG(kwds)
        var, contset, profile = config.var, config.contset, config.profile
        if var is None and contset is None:
            found = 0
            for block in model.block_data_objects(active=True, descend_into=True):
                store = _cvp_data(block)
                decls = store.get("profiles") if store else None
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
        """Parameterize one control ``var`` over ``contset`` in place.

        Replaces ``var`` with a Var indexed by the profile's free time
        points, substitutes every eliminated copy out of all constraints,
        objectives, and named expressions (the differential map on
        expressions carrying a DerivativeVar, the algebraic map elsewhere),
        and adds explicit bound constraints for eliminated copies whose
        substitution is a non-convex combination of the free values. Records
        profile metadata for :func:`control_value`.

        Parameters
        ----------
        model : Block
            The model being transformed.
        var : Var
            The control to parameterize.
        contset : ContinuousSet
            The discretized time set ``var`` is indexed over.
        profile : str or tuple
            The profile specifier (see :func:`declare_profile`).

        Returns
        -------
        Block
            ``model``, transformed in place.

        Raises
        ------
        TypeError
            If ``var`` or ``contset`` is missing or of the wrong type.
        ValueError
            If ``var`` is not indexed by ``contset`` exactly once, carries a
            DerivativeVar, or a piecewise-constant control is referenced at
            the final time.
        RuntimeError
            If ``contset`` is not discretized, or a generated bound-
            constraint name collides with an existing component.
        """
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
            return full[:pos] + (t,) + full[pos + 1 :]

        old = {as_tuple(i): var[i] for i in var}
        free_attrs = {}
        for full, vd in old.items():
            if full[pos] in free_times:
                free_attrs[full] = (vd.domain, vd.lb, vd.ub, vd.value)

        name = var.local_name
        parent = var.parent_block()
        parent.del_component(var)
        newsets = subsets[:pos] + [free_times] + subsets[pos + 1 :]
        any_dom = next(iter(free_attrs.values()))[0]
        new_var = Var(
            *newsets,
            domain=any_dom,
            bounds=lambda m, *i: (free_attrs[i][1], free_attrs[i][2]),
            initialize=lambda m, *i: free_attrs[i][3],
        )
        parent.add_component(name, new_var)

        submap = {}
        # (expr, lb, ub) for eliminated copies whose substitution is a
        # NON-convex combination of the free values: there the original
        # variable bounds are not implied by the knot bounds and must be kept
        # as explicit constraints (matching reduce_collocation_points, which
        # enforces bounds at every collocation point).
        bound_rows = []
        for full, vd in old.items():
            plist = pairs[full[pos]]
            if len(plist) == 1 and plist[0][1] == 1.0:
                node = new_var[with_t(full, plist[0][0])]
            else:
                node = sum(c * new_var[with_t(full, ft)] for ft, c in plist)
                convex = all(-1e-12 <= c <= 1 + 1e-12 for _, c in plist)
                if not convex and (vd.lb is not None or vd.ub is not None):
                    bound_rows.append((node, vd.lb, vd.ub))
            submap[id(vd)] = node

        # submap is the DIFFERENTIAL map: at a node it gives the control the
        # ODE there integrates with. It is applied to constraints that carry a
        # DerivativeVar (the user's ODEs and pyomo.dae's disc_eq).
        #
        # submap_alg is the ALGEBRAIC map, for objectives / path constraints /
        # non-differential expressions: a control reference names the control
        # of the element it sits in (the decision the user wrote). The two maps
        # differ only for the discontinuous piecewise_constant profile, and
        # only at a finite-element boundary: the differential map takes the
        # element ENDING at the node (left-continuous, correct for the
        # collocation equation there), the algebraic map the element STARTING
        # at it. For the continuous profiles the maps coincide. The final node
        # starts no element, so an algebraic reference to it (which cannot mean
        # a real control) is collected in `forbidden` and raised, rather than
        # silently mapped to the last element and double-counted.
        submap_alg = submap
        forbidden = {}
        if mode == "piecewise_constant":
            submap_alg = {}
            for full, vd in old.items():
                j = bisect_right(fe, full[pos]) - 1
                if j >= len(fe) - 1:
                    forbidden[id(vd)] = vd.name
                else:
                    submap_alg[id(vd)] = new_var[with_t(full, fe[j])]

        deriv_ids = {
            id(d)
            for comp in model.component_objects(Var, active=True)
            if isinstance(comp, DerivativeVar)
            for d in comp.values()
        }

        def is_differential(expr):
            return any(id(v) in deriv_ids for v in identify_variables(expr))

        def check_final(expr, where):
            for v in identify_variables(expr):
                nm = forbidden.get(id(v))
                if nm is not None:
                    raise ValueError(
                        f"pyomo-cvp: {where} references '{nm}', the control at "
                        f"the final time. A piecewise-constant control has one "
                        f"value per element (indexed by the element starts, "
                        f"0..N-1); the final node carries no control. Reference "
                        f"controls over the elements, not the final boundary."
                    )

        for c in model.component_data_objects(
            Constraint, active=True, descend_into=True
        ):
            if is_differential(c.expr):
                c.set_value(replace_expressions(c.expr, submap))
            else:
                if forbidden:
                    check_final(c.expr, f"constraint '{c.name}'")
                c.set_value(replace_expressions(c.expr, submap_alg))
        for o in model.component_data_objects(
            Objective, active=True, descend_into=True
        ):
            if forbidden:
                check_final(o.expr, f"objective '{o.name}'")
            o.set_value(replace_expressions(o.expr, submap_alg))
        for e in model.component_data_objects(
            Expression, active=True, descend_into=True
        ):
            if is_differential(e.expr):
                e.set_value(replace_expressions(e.expr, submap))
            else:
                if forbidden:
                    check_final(e.expr, f"expression '{e.name}'")
                e.set_value(replace_expressions(e.expr, submap_alg))

        if bound_rows:
            clname = name + "_profile_bounds"
            if parent.find_component(clname) is not None:
                raise RuntimeError(f"pyomo-cvp: component '{clname}' already exists.")
            cl = ConstraintList()
            parent.add_component(clname, cl)
            for node, lb, ub in bound_rows:
                if lb is not None and ub is not None:
                    cl.add((lb, node, ub))
                elif lb is not None:
                    cl.add(node >= lb)
                else:
                    cl.add(node <= ub)

        info = _cvp_data(parent, create=True).setdefault("info", {})
        info[name] = {"mode": mode, "k": k, "fe": fe, "pos": pos, "pairs": pairs}

        return model

    @staticmethod
    def _profile_pairs(mode, k, fe, points, disc_info, contset):
        """Map each time point to its substitution ``[(free_time, coeff), ...]``.

        For each discretized time point, gives the linear combination of free
        control values the eliminated copy there is replaced by: a single
        ``(free_time, 1.0)`` pair when the value stays free, or several pairs
        (the interpolation weights) when it is eliminated.

        Parameters
        ----------
        mode : str
            The normalized profile name.
        k : int or None
            Collocation order. ``None`` for the piecewise profiles, and for
            plain ``'collocation'``, where it resolves to ncp.
        fe : list of float
            The finite-element boundaries.
        points : list of float
            All discretized time points, sorted.
        disc_info : dict
            ``contset.get_discretization_info()``.
        contset : ContinuousSet
            The discretized time set.

        Returns
        -------
        dict
            Maps each time point to a list of ``(free_time, coeff)`` pairs.

        Raises
        ------
        RuntimeError
            If a collocation profile is used without a collocation
            discretization.
        ValueError
            If the collocation order exceeds the discretization ncp.
        """
        pairs = {}

        def element_of(t):
            if t <= fe[0]:
                return 0
            return bisect_left(fe, t) - 1

        if mode == "piecewise_constant":
            # u(t) = u_k on (t_k, t_{k+1}]: the free value is anchored at the
            # element's START, so controls index by fe[:-1] and u[t0] exists.
            # Boundary points map to the element they terminate (as the Radau
            # and backward-difference equations require), so the evaluated
            # profile is left-continuous at interior boundaries.
            for t in points:
                pairs[t] = [(fe[element_of(t)], 1.0)]
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

        # collocation
        ncp = disc_info.get("ncp")
        if ncp is None:
            raise RuntimeError(
                "pyomo-cvp: the 'collocation' profile requires a collocation "
                "discretization (dae.collocation), not finite differences."
            )
        if k is None:
            k = ncp
        if k > ncp:
            raise ValueError(
                f"pyomo-cvp: collocation k={k} exceeds the "
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
