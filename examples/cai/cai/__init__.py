"""
CompDam_DGD Compression-After-Impact (CAI) model generator.

The package keeps pure-Python logic, which you can test without Abaqus, apart
from the builders that call the Abaqus Scripting Interface:

Pure Python (unit tested with ``uv run pytest``)
    :mod:`cai.settings`, :mod:`cai.compdam`, :mod:`cai.layup`,
    :mod:`cai.damage_map`, :mod:`cai.inp_patcher`, :mod:`cai.validation`

Abaqus kernel only
    :mod:`cai.geometry`, :mod:`cai.material_section`, :mod:`cai.damage`,
    :mod:`cai.interactions`, :mod:`cai.loading`, :mod:`cai.orchestrator`,
    :mod:`cai.postprocess`
"""

__version__ = "0.1.0"


class CAIError(Exception):
    """Base class for all errors raised by the CAI model generator."""


class ConfigurationError(CAIError):
    """Raised when the configuration is incomplete or inconsistent."""


class ModelBuildError(CAIError):
    """Raised when an Abaqus model-building operation fails."""
