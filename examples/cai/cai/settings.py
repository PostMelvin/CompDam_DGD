"""
Loading and validation of the configuration module.

The configuration is a plain Python module (``config.py`` by default). This
module turns it into a dictionary, applies command-line overrides and checks
that the values are consistent before any Abaqus object is created, so
mistakes show up in seconds rather than after meshing.
"""

from __future__ import division, print_function

import copy
import os
import runpy

from cai import ConfigurationError
from cai import compdam

#: Top-level sections every configuration must define.
REQUIRED_SECTIONS = (
    "MODEL", "GEOMETRY", "LAYUP", "MATERIALS", "NCAMP_BASELINE", "COMPDAM",
    "MESH", "IMPACT", "DAMAGE", "BOUNDARY", "STEP", "INTERFACE", "JOB",
)


def load_config(path=None):
    """Load a configuration module into a dictionary.

    Parameters
    ----------
    path : str, optional
        Path to a configuration file. Defaults to ``config.py`` next to the
        package.

    Returns
    -------
    dict
        Deep copy of the upper-case sections, plus ``"_config_dir"`` (the
        directory that holds the file, used to resolve relative paths).

    Raises
    ------
    ConfigurationError
        If the file does not exist or a required section is missing.
    """
    if path is None:
        path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                            "config.py")
    path = os.path.abspath(path)
    if not os.path.isfile(path):
        raise ConfigurationError("configuration file not found: %s" % path)
    namespace = runpy.run_path(path)
    cfg = dict((k, copy.deepcopy(v)) for k, v in namespace.items() if k.isupper())
    missing = [s for s in REQUIRED_SECTIONS if s not in cfg]
    if missing:
        raise ConfigurationError("configuration %s lacks sections: %s" % (path, ", ".join(missing)))
    cfg["_config_dir"] = os.path.dirname(path)
    return cfg


def resolve_path(cfg, path):
    """Resolve a path from the configuration relative to the config file.

    Parameters
    ----------
    cfg : dict
        Loaded configuration.
    path : str
        Absolute path, or path relative to the configuration file.

    Returns
    -------
    str
        Absolute, normalised path.
    """
    if os.path.isabs(path):
        return os.path.normpath(path)
    return os.path.normpath(os.path.join(cfg.get("_config_dir", os.getcwd()), path))


def active_material(cfg):
    """Return ``(name, preset)`` of the material selected in ``LAYUP``.

    Raises
    ------
    ConfigurationError
        If the material name is not in ``MATERIALS``.
    """
    name = cfg["LAYUP"]["material"]
    if name not in cfg["MATERIALS"]:
        raise ConfigurationError("material %r not in MATERIALS (%s)" % (
            name, ", ".join(sorted(cfg["MATERIALS"]))))
    return name, cfg["MATERIALS"][name]


def impact_center(cfg):
    """Impact location ``(x, y)``; defaults to the plate centre."""
    g = cfg["GEOMETRY"]
    c = g.get("impact_center")
    if c is None:
        return (0.5 * g["length"], 0.5 * g["width"])
    return (float(c[0]), float(c[1]))


def validate(cfg):
    """Check the configuration for consistency.

    Parameters
    ----------
    cfg : dict
        Loaded configuration.

    Returns
    -------
    list of str
        Non-fatal warnings (e.g. a coarse mesh that may snap back).

    Raises
    ------
    ConfigurationError
        On any fatal inconsistency.
    """
    warnings = []
    g, m, b = cfg["GEOMETRY"], cfg["MESH"], cfg["BOUNDARY"]
    L, W = float(g["length"]), float(g["width"])
    if L <= 0 or W <= 0:
        raise ConfigurationError("plate length and width must be positive")

    xc, yc = impact_center(cfg)
    lz, wz = g["refine_zone_length"], g["refine_zone_width"]
    es, inset = b["end_support_length"], b["side_support_inset"]
    if not (es < xc - lz / 2.0 and xc + lz / 2.0 < L - es):
        raise ConfigurationError(
            "the refined zone (x = %.1f..%.1f) must lie between the end supports "
            "(%.1f and %.1f)" % (xc - lz / 2.0, xc + lz / 2.0, es, L - es))
    if not (inset < yc - wz / 2.0 and yc + wz / 2.0 < W - inset):
        raise ConfigurationError(
            "the refined zone (y = %.1f..%.1f) must lie between the side supports"
            % (yc - wz / 2.0, yc + wz / 2.0))
    if es <= 0.0:
        raise ConfigurationError("end_support_length must be positive")
    if inset < 0.0:
        raise ConfigurationError("side_support_inset must be >= 0")

    if not 0 < m["fine_size"] <= m["coarse_size"]:
        raise ConfigurationError("need 0 < fine_size <= coarse_size")
    if int(m["elements_per_ply"]) < 1:
        raise ConfigurationError("elements_per_ply must be >= 1")
    if m["element_type"] != "C3D8R":
        raise ConfigurationError("only C3D8R is implemented (got %r)" % m["element_type"])

    name, mat = active_material(cfg)
    flags = cfg["COMPDAM"]["feature_flags"]
    compdam.build_property_vector(mat, flags)  # raises on missing properties
    f = compdam.parse_feature_flags(flags)
    if f["fiber_compression"] in (3, 4, 5):
        warnings.append(
            "FKT fibre kinking is enabled: CompDam expects a fibre-aligned, uniform mesh; "
            "the orthogonal mesh built here is fibre-aligned only for 0/90 plies")

    limits = compdam.max_element_size(mat)
    for mode, lmax in sorted(limits.items()):
        if m["fine_size"] > lmax:
            warnings.append(
                "fine_size %.2f mm exceeds the %s snap-back limit %.2f mm for %s"
                % (m["fine_size"], mode, lmax, name))
    if limits and m["coarse_size"] > min(limits.values()):
        warnings.append(
            "coarse_size %.2f mm exceeds the snap-back limit %.2f mm: damage must stay "
            "inside the refined zone" % (m["coarse_size"], min(limits.values())))

    if cfg["INTERFACE"]["type"] not in ("cohesive", "tie"):
        raise ConfigurationError("INTERFACE type must be 'cohesive' or 'tie'")
    if cfg["STEP"]["time_period"] <= 0:
        raise ConfigurationError("step time_period must be positive")
    if b["applied_displacement"] <= 0:
        raise ConfigurationError("applied_displacement is an end shortening and must be > 0")
    if cfg["IMPACT"]["reference_energy"] <= cfg["IMPACT"].get("threshold_energy", 0.0):
        raise ConfigurationError("IMPACT reference_energy must exceed threshold_energy")
    return warnings


def apply_overrides(cfg, material=None, energy=None, submit=None, interface=None,
                    job_name=None, cpus=None):
    """Apply command-line overrides to a loaded configuration in place.

    Parameters
    ----------
    cfg : dict
        Loaded configuration.
    material : str, optional
        Material preset name.
    energy : float, optional
        Impact energy [J] (sets ``IMPACT["default_energy"]``).
    submit : bool, optional
        Submit the job after writing the input deck.
    interface : str, optional
        ``"cohesive"`` or ``"tie"``.
    job_name : str, optional
        Base name of the model and job.
    cpus : int, optional
        Number of CPUs and domains.

    Returns
    -------
    dict
        The same ``cfg`` object.
    """
    if material is not None:
        cfg["LAYUP"]["material"] = material
    if energy is not None:
        cfg["IMPACT"]["default_energy"] = float(energy)
    if submit is not None:
        cfg["JOB"]["submit"] = bool(submit)
    if interface is not None:
        cfg["INTERFACE"]["type"] = interface
    if job_name is not None:
        cfg["MODEL"]["job_name"] = job_name
        cfg["MODEL"]["model_name"] = job_name
    if cpus is not None:
        cfg["JOB"]["num_cpus"] = int(cpus)
        cfg["JOB"]["num_domains"] = int(cpus)
    return cfg


def run_name(cfg, energy):
    """Model/job name for one impact energy, e.g. ``CAI_QI_AS4_8552_30J``.

    Parameters
    ----------
    cfg : dict
        Loaded configuration.
    energy : float
        Impact energy [J].

    Returns
    -------
    str
        Abaqus-safe name (letters, digits, underscores; starts with a letter).
    """
    e = ("%g" % energy).replace(".", "p")
    name = "%s_%s_%sJ" % (cfg["MODEL"]["job_name"], cfg["LAYUP"]["material"], e)
    return "".join(c if c.isalnum() or c == "_" else "_" for c in name)[:38]
