"""
Read a finished CAI output database and compare the result with NCAMP (Abaqus Python).

Run with ``abaqus python cai_postprocess.py -- <job>.odb``. Extracts the
load-displacement history of ``RP_LOAD``, the peak load, the CAI strength
(ASTM D7137) and the kinetic-to-internal energy ratio (quasi-static check),
writes a CSV and prints the comparison with ``config.NCAMP_BASELINE``.
"""

from __future__ import division, print_function

import csv
import os

from cai import validation


def _find_region(step, required):
    for key, region in step.historyRegions.items():
        if all(v in region.historyOutputs.keys() for v in required):
            return key, region
    return None, None


def extract(odb_path):
    """Extract load-displacement and energy histories from an ODB.

    Parameters
    ----------
    odb_path : str
        Path to the ``.odb`` file.

    Returns
    -------
    dict
        ``time``, ``displacement`` (end shortening, positive), ``load``
        (compressive, positive) and ``ke_ie_ratio_max`` (max ALLKE/ALLIE over
        the second half of the step, or None if unavailable).

    Raises
    ------
    RuntimeError
        If the RP history output is missing.
    """
    from odbAccess import openOdb

    odb = openOdb(odb_path, readOnly=True)
    try:
        step = odb.steps[list(odb.steps.keys())[-1]]
        _, rp = _find_region(step, ("U1", "RF1"))
        if rp is None:
            raise RuntimeError("no history region with U1 and RF1 (RP_LOAD) in %s" % odb_path)
        u = rp.historyOutputs["U1"].data
        f = rp.historyOutputs["RF1"].data
        out = {
            "time": [d[0] for d in u],
            "displacement": [-d[1] for d in u],
            "load": [-d[1] for d in f],
            "ke_ie_ratio_max": None,
        }
        _, en = _find_region(step, ("ALLKE", "ALLIE"))
        if en is not None:
            ke = en.historyOutputs["ALLKE"].data
            ie = en.historyOutputs["ALLIE"].data
            half = len(ie) // 2
            ratios = [k[1] / i[1] for k, i in zip(ke[half:], ie[half:]) if abs(i[1]) > 1e-12]
            out["ke_ie_ratio_max"] = max(ratios) if ratios else None
        return out
    finally:
        odb.close()


def write_csv(path, data):
    """Write time, displacement and load columns to a CSV file.

    Parameters
    ----------
    path : str
        Output CSV path.
    data : dict
        Result of :func:`extract`.
    """
    mode = "w" if str is not bytes else "wb"  # Python 3 / Python 2.7 (Abaqus 2023)
    kwargs = {"newline": ""} if str is not bytes else {}
    with open(path, mode, **kwargs) as fh:
        w = csv.writer(fh)
        w.writerow(["time_s", "end_shortening_mm", "load_N"])
        for row in zip(data["time"], data["displacement"], data["load"]):
            w.writerow(["%.6e" % v for v in row])


def summarize(odb_path, cfg, energy=None, thickness=None):
    """Extract results, write ``<job>_load_displacement.csv`` and return a report.

    Parameters
    ----------
    odb_path : str
        Path to the ``.odb`` file.
    cfg : dict
        Loaded configuration (for width, material and NCAMP baseline).
    energy : float, optional
        Impact energy used in the model [J].
    thickness : float, optional
        Laminate thickness [mm].

    Returns
    -------
    str
        Human-readable report.
    """
    data = extract(odb_path)
    csv_path = os.path.splitext(odb_path)[0] + "_load_displacement.csv"
    write_csv(csv_path, data)
    peak = max(data["load"]) if data["load"] else 0.0
    width = cfg["GEOMETRY"]["width"]
    lines = ["Results for %s" % odb_path, "Peak compressive load: %.0f N" % peak]
    if thickness:
        strength = validation.cai_strength(peak, width, thickness)
        baseline = cfg["NCAMP_BASELINE"].get(cfg["LAYUP"]["material"])
        lines.append(validation.cai_report(strength, baseline, energy, thickness))
    r = data["ke_ie_ratio_max"]
    if r is not None:
        lines.append("Max ALLKE/ALLIE (2nd half of step): %.3f %s" % (
            r, "(OK, quasi-static)" if r < 0.05 else "(> 5 %: slow the loading down)"))
    lines.append("Load-displacement CSV: %s" % csv_path)
    return "\n".join(lines)
