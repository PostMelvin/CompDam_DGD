"""
Comparison of model inputs and results with NCAMP baseline data.

Pure Python. :func:`material_report` lists the CompDam inputs against the
NCAMP lamina values in ``config.NCAMP_BASELINE``; :func:`cai_report`
compares a predicted CAI strength with the NCAMP laminate value.
"""

from __future__ import division, print_function

from cai import compdam

#: Mapping of NCAMP baseline keys to CompDam property names.
NCAMP_TO_COMPDAM = {
    "E1t": "E1", "E1c": "E1", "E2t": "E2", "G12": "G12", "nu12": "nu12",
    "XT": "XT", "XC": "XC", "YT": "YT", "YC": "YC", "SL": "SL",
}


def material_report(material_name, material, baseline):
    """Tabulate CompDam inputs against NCAMP lamina baseline values.

    Parameters
    ----------
    material_name : str
        Preset name.
    material : dict
        Material preset from ``config.MATERIALS``.
    baseline : dict or None
        Entry of ``config.NCAMP_BASELINE`` for the material.

    Returns
    -------
    str
        Fixed-width text table with the relative difference of each property.
    """
    if not baseline or not baseline.get("lamina"):
        return "No NCAMP lamina baseline defined for %s." % material_name
    flat = compdam.flatten_material(material)
    flat["ply_thickness"] = material.get("ply_thickness")
    lines = [
        "NCAMP comparison for %s: %s, %s %s%s" % (
            material_name, baseline.get("report", "?"), baseline.get("condition", ""),
            baseline.get("basis", ""),
            "" if baseline.get("verified") else "  [baseline NOT verified]"),
        "%-14s %12s %12s %9s" % ("property", "NCAMP", "model", "delta"),
    ]
    for key in sorted(baseline["lamina"]):
        ref = baseline["lamina"][key]
        model_key = NCAMP_TO_COMPDAM.get(key, key)
        val = flat.get(model_key)
        if ref is None or val is None:
            lines.append("%-14s %12s %12s %9s" % (key, ref, val, "-"))
            continue
        delta = 100.0 * (val - ref) / ref if ref else float("nan")
        lines.append("%-14s %12.4g %12.4g %8.1f%%" % (key, ref, val, delta))
    return "\n".join(lines)


def cai_strength(peak_load, width, thickness):
    """CAI strength per ASTM D7137: ``F_CAI = P_max / (w * h)``.

    Parameters
    ----------
    peak_load : float
        Maximum compressive load magnitude [N].
    width, thickness : float
        Specimen width and thickness [mm].

    Returns
    -------
    float
        CAI strength [MPa].
    """
    return abs(peak_load) / (width * thickness)


def cai_report(predicted_strength, baseline, energy=None, thickness=None):
    """Compare a predicted CAI strength with the NCAMP baseline.

    Parameters
    ----------
    predicted_strength : float
        Model CAI strength [MPa].
    baseline : dict or None
        Entry of ``config.NCAMP_BASELINE``.
    energy : float, optional
        Impact energy the model used [J].
    thickness : float, optional
        Laminate thickness [mm]; with ``energy``, used to report J/mm.

    Returns
    -------
    str
        Short summary.
    """
    lines = ["Predicted CAI strength: %.1f MPa" % predicted_strength]
    if energy is not None and thickness:
        lines.append("Impact energy: %.1f J (%.2f J/mm)" % (energy, energy / thickness))
    cai = (baseline or {}).get("cai", {})
    ref = cai.get("cai_strength")
    if ref:
        lines.append("NCAMP CAI strength: %.1f MPa at %.2f J/mm -> delta %.1f%%%s" % (
            ref, cai.get("impact_energy_per_thickness", float("nan")),
            100.0 * (predicted_strength - ref) / ref,
            "" if baseline.get("verified") else "  [baseline NOT verified]"))
    else:
        lines.append("No NCAMP CAI reference value configured.")
    return "\n".join(lines)
