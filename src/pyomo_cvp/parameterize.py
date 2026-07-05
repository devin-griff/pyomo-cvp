"""The ``cvp.parameterize`` transformation and the ``declare_profile`` API.

Phase 0 scaffold: registration and argument validation only; the substitution
mechanics land in Phase 1.
"""
from pyomo.core import Transformation, TransformationFactory

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
        raise NotImplementedError(
            "pyomo-cvp: the substitution mechanics arrive in Phase 1."
        )
