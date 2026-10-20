"""
Standalone CompDam_DGD Compression-After-Impact (CAI) model generator.

Plain Python (2.7 or 3.x, no third-party packages, no Abaqus/CAE): it builds
fibre-aligned ply meshes, sizes the delamination zone with Turon's
cohesive-zone equations and writes a complete Abaqus/Explicit input deck,
plus the files to run it with a pre-compiled CompDam library on a cluster.

Modules
-------
:mod:`cai.settings`     configuration loading, validation, overrides
:mod:`cai.compdam`      CompDam property vector, flags, state variables
:mod:`cai.layup`        plies, ply blocks, interfaces, CLT
:mod:`cai.damage_map`   seeded delaminations and matrix cracks
:mod:`cai.turon`        cohesive-zone length, element size, strength, penalty
:mod:`cai.mesh2d`       fibre-aligned structured quad mesh of one ply
:mod:`cai.model`        3-D extrusion, numbering, sets, surfaces
:mod:`cai.inp_writer`   Abaqus/Explicit input deck
:mod:`cai.jobfiles`     abaqus_v6.env, SLURM/PBS scripts, postprocess.py
:mod:`cai.builder`      orchestration (CAIJob)
:mod:`cai.validation`   NCAMP comparison helpers
"""

__version__ = "0.2.0"


class CAIError(Exception):
    """Base class for all errors raised by the CAI model generator."""


class ConfigurationError(CAIError):
    """Raised when the configuration is incomplete or inconsistent."""


class ModelBuildError(CAIError):
    """Raised when an Abaqus model-building operation fails."""
