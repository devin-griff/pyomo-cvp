"""pyomo-cvp: control vector parameterization for pyomo.dae.

Importing this package registers the ``cvp.parameterize`` transformation.
"""
from pyomo_cvp import parameterize  # noqa: F401  (registers the plugin)
from pyomo_cvp.parameterize import declare_profile  # noqa: F401

try:
    from importlib.metadata import version

    __version__ = version("pyomo-cvp")
except Exception:  # noqa: BLE001
    __version__ = "unknown"
