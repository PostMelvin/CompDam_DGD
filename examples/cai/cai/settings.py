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
    "MESH", "COHESIVE_ZONE", "IMPACT", "DAMAGE", "BOUNDARY", "STEP", "INTERFACE",
    "JOB", "CLUSTER",
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
        directory that holds the file).

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
        Non-fatal warnings.

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
    es, inset = float(b["end_support_length"]), float(b["side_support_inset"])
    if not 0.0 < es < 0.5 * L:
        raise ConfigurationError("end_support_length must lie in (0, length/2)")
    if not 0.0 <= inset < 0.5 * W:
        raise ConfigurationError("side_support_inset must lie in [0, width/2)")
    if not (es < xc < L - es and inset < yc < W - inset):
        raise ConfigurationError("the impact centre must lie inside the supports")

    if m["element_type"] != "C3D8R":
        raise ConfigurationError("only C3D8R is implemented (got %r)" % m["element_type"])
    if float(m["coarse_size"]) <= 0 or float(m.get("growth_rate", 1.2)) <= 1.0:
        raise ConfigurationError("need coarse_size > 0 and growth_rate > 1")
    if int(m.get("elements_per_ply", 1)) < 1:
        raise ConfigurationError("elements_per_ply must be >= 1")
    if float(m.get("zone_margin", 0.0)) < 0.0:
        raise ConfigurationError("zone_margin must be >= 0")
    hg = str(m.get("hourglass_control", "ENHANCED")).upper()
    if hg not in ("ENHANCED", "RELAX_STIFFNESS", "STIFFNESS", "DEFAULT"):
        raise ConfigurationError("unknown hourglass_control %r" % hg)

    name, mat = active_material(cfg)
    flags = cfg["COMPDAM"]["feature_flags"]
    compdam.property_values(mat, flags)  # raises on missing properties
    f = compdam.parse_feature_flags(flags)
    if f["fiber_compression"] in (3, 4, 5):
        warnings.append(
            "FKT fibre kinking is enabled: the random misalignment options (phi0 = 1/2/3) "
            "need a uniform fibre-aligned mesh, which holds only inside the damage zones")
        if not m.get("fiber_aligned", True):
            raise ConfigurationError("FKT needs MESH fiber_aligned = True")

    cz = cfg["COHESIVE_ZONE"]
    if float(cz.get("elements_in_zone", 3)) < 1:
        raise ConfigurationError("COHESIVE_ZONE elements_in_zone must be >= 1")
    if cz.get("min_element_size") is not None and float(cz["min_element_size"]) <= 0:
        raise ConfigurationError("COHESIVE_ZONE min_element_size must be positive")

    c = cfg["COMPDAM"]
    lib, src = c.get("library_dir"), c.get("user_subroutine")
    if not lib and not src:
        raise ConfigurationError("set COMPDAM library_dir (pre-compiled) or user_subroutine")
    if lib and "/path/to/" in str(lib):
        warnings.append("COMPDAM library_dir is still the placeholder %r; set it to the "
                        "directory of the compiled CompDam library on the cluster" % lib)

    if cfg["INTERFACE"]["type"] not in ("cohesive", "tie"):
        raise ConfigurationError("INTERFACE type must be 'cohesive' or 'tie'")
    if cfg["STEP"]["time_period"] <= 0:
        raise ConfigurationError("step time_period must be positive")
    if b["applied_displacement"] <= 0:
        raise ConfigurationError("applied_displacement is an end shortening and must be > 0")
    if cfg["IMPACT"]["reference_energy"] <= cfg["IMPACT"].get("threshold_energy", 0.0):
        raise ConfigurationError("IMPACT reference_energy must exceed threshold_energy")
    sched = str(cfg["CLUSTER"].get("scheduler") or "none").lower()
    if sched not in ("slurm", "pbs", "none"):
        raise ConfigurationError("CLUSTER scheduler must be 'slurm', 'pbs' or 'none'")
    return warnings


def apply_overrides(cfg, material=None, energy=None, interface=None, job_name=None,
                    cpus=None, output_dir=None, library_dir=None, min_element_size=None):
    """Apply command-line overrides to a loaded configuration in place.

    Parameters
    ----------
    cfg : dict
        Loaded configuration.
    material : str, optional
        Material preset name.
    energy : float, optional
        Impact energy [J] (sets ``IMPACT["default_energy"]``).
    interface : str, optional
        ``"cohesive"`` or ``"tie"``.
    job_name : str, optional
        Base job name.
    cpus : int, optional
        CPUs (and Explicit domains) in the cluster scripts.
    output_dir : str, optional
        Parent folder of the job folders.
    library_dir : str, optional
        Directory of the compiled CompDam library on the cluster.
    min_element_size : float, optional
        Floor on the damage-zone element size [mm].

    Returns
    -------
    dict
        The same ``cfg`` object.
    """
    if material is not None:
        cfg["LAYUP"]["material"] = material
    if energy is not None:
        cfg["IMPACT"]["default_energy"] = float(energy)
    if interface is not None:
        cfg["INTERFACE"]["type"] = interface
    if job_name is not None:
        cfg["MODEL"]["job_name"] = job_name
    if cpus is not None:
        cfg["CLUSTER"]["cpus"] = int(cpus)
    if output_dir is not None:
        cfg["JOB"]["output_dir"] = output_dir
    if library_dir is not None:
        cfg["COMPDAM"]["library_dir"] = library_dir
    if min_element_size is not None:
        cfg["COHESIVE_ZONE"]["min_element_size"] = float(min_element_size)
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
