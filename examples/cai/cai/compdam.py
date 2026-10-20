"""
CompDam_DGD interface definitions.

Maps material presets onto the 40-entry ``*User Material`` property vector,
parses feature flags, sizes the state-variable array and formats the
``*Initial Conditions, type=SOLUTION`` data CompDam needs.

All definitions follow the CompDam_DGD readme (sections "Material properties",
"Controlling which features are enabled", "State variables" and "Initial
conditions").
"""

from __future__ import division, print_function

from cai import ConfigurationError

#: Order of the 40 CompDam material properties (1-based position = index + 1).
#: Positions 2 and 4 are reserved; 5-8 are fatigue properties.
PROPERTY_ORDER = (
    "feature_flags", "reserved_2", "thickness", "reserved_4",
    "fatigue_gamma", "fatigue_epsilon", "fatigue_eta", "fatigue_p_mod",
    "E1", "E2", "G12", "nu12", "nu23", "YT", "SL", "GYT",
    "GSL", "eta_BK", "YC", "alpha0", "E3", "G13", "G23", "nu13",
    "alpha11", "alpha22", "alpha_PL", "n_PL", "XT", "fXT", "GXT", "fGXT",
    "XC", "fXC", "GXC", "fGXC", "cl", "w_kb", "T_sf", "mu",
)

#: CompDam state variable names (continuum elements), 1-based in Abaqus.
STATE_VARIABLE_NAMES = (
    "CDM_d2", "CDM_Fb1", "CDM_Fb2", "CDM_Fb3", "CDM_B", "CDM_Lc1", "CDM_Lc2",
    "CDM_Lc3", "CDM_FIm", "CDM_alpha", "CDM_STATUS", "CDM_Plas12",
    "CDM_Inel12", "CDM_FIfT", "CDM_slide1", "CDM_slide2", "CDM_FIfC",
    "CDM_d1T", "CDM_d1C", "CDM_Plas13", "CDM_Inel13", "CDM_phi0_12",
    "CDM_gamma_12", "CDM_phi0_13", "CDM_gamma_13", "CDM_reserved", "CDM_Ep1",
    "CDM_Ep2", "CDM_Ep3", "CDM_Ep4", "CDM_Ep5", "CDM_Ep6", "CDM_fp1",
)

#: 1-based state-variable positions used when seeding damage.
SDV_D2 = 1
SDV_ALPHA = 10
SDV_STATUS = 11
SDV_D1T = 18
SDV_D1C = 19
SDV_PHI0_12 = 22
SDV_PHI0_13 = 24

#: Properties that must be non-blank for each enabled feature.
FEATURE_REQUIREMENTS = {
    "elastic": ("E1", "E2", "G12", "nu12", "nu23"),
    "matrix_damage": ("YT", "SL", "GYT", "GSL", "eta_BK", "YC", "alpha0"),
    "shear_ro": ("alpha_PL", "n_PL"),
    "fiber_tension": ("XT", "fXT", "GXT", "fGXT"),
    "fiber_compression_bl": ("XC", "fXC", "GXC", "fGXC"),
    "fiber_compression_fkt": ("XC", "YC", "w_kb", "alpha0", "alpha_PL", "n_PL"),
    "friction": ("mu",),
}

_FLAG_NAMES = (
    "matrix_damage", "shear_nonlinearity", "fiber_tension",
    "fiber_compression", "energy_output", "friction",
)


def parse_feature_flags(flags):
    """Parse a CompDam feature-flag integer/string into its six positions.

    Parameters
    ----------
    flags : int or str
        Feature flags, e.g. ``"111101"`` or ``111101``. Missing leading zeros
        are restored the same way CompDam does it (``11000`` -> ``"011000"``).

    Returns
    -------
    dict
        Mapping of feature name to integer option, with keys
        ``matrix_damage``, ``shear_nonlinearity``, ``fiber_tension``,
        ``fiber_compression``, ``energy_output`` and ``friction``.

    Raises
    ------
    ConfigurationError
        If the flags are not a non-negative integer of at most six digits or
        an option is outside the range CompDam supports.
    """
    text = str(flags).strip()
    if not text.isdigit() or len(text.lstrip("0")) > 6:
        raise ConfigurationError(
            "feature_flags must be an integer of up to 6 digits, got %r" % (flags,))
    text = text[-6:].zfill(6)
    parsed = dict(zip(_FLAG_NAMES, [int(c) for c in text]))

    allowed = {
        "matrix_damage": (0, 1),  # 2 = cohesive elements, not used here
        "shear_nonlinearity": (0, 1, 2, 3, 4, 5),
        "fiber_tension": (0, 1),
        "fiber_compression": (0, 1, 3, 4, 5),
        "energy_output": (0, 1, 2),
        "friction": (0, 1),
    }
    for name, value in parsed.items():
        if value not in allowed[name]:
            raise ConfigurationError(
                "feature flag %s=%d is not supported for solid elements "
                "(allowed: %s)" % (name, value, allowed[name]))
    return parsed


def feature_flags_as_float(flags):
    """Return the feature flags as the float stored in the property vector.

    Parameters
    ----------
    flags : int or str
        Feature flags.

    Returns
    -------
    float
        Numeric value written as material property 1.
    """
    parse_feature_flags(flags)
    return float(int(str(flags).strip()))


def number_of_state_variables(flags):
    """Number of state variables (``*Depvar``) needed for the given features.

    Parameters
    ----------
    flags : int or str
        Feature flags.

    Returns
    -------
    int
        18 minimum; 19 for fibre compression model 1; 21 for Ramberg-Osgood
        3-D/1-3 shear nonlinearity; 25 for fibre kinking (FKT); 33 for Schaefer.
    """
    f = parse_feature_flags(flags)
    n = 18
    if f["fiber_compression"] == 1:
        n = max(n, 19)
    if f["shear_nonlinearity"] in (3, 4):
        n = max(n, 21)
    if f["fiber_compression"] in (3, 4, 5):
        n = max(n, 25)
    if f["shear_nonlinearity"] == 5:
        n = max(n, 33)
    return n


def flatten_material(material):
    """Collapse the grouped material preset into a single property dictionary.

    Parameters
    ----------
    material : dict
        Preset from ``config.MATERIALS`` with ``elastic``, ``strength``,
        ``fracture`` and ``other`` groups.

    Returns
    -------
    dict
        ``{property_name: value}`` for every CompDam property found.

    Raises
    ------
    ConfigurationError
        If a property name appears in more than one group.
    """
    flat = {}
    for group in ("elastic", "strength", "fracture", "other"):
        for key, value in material.get(group, {}).items():
            if key in flat:
                raise ConfigurationError("property %r defined twice in material" % key)
            flat[key] = value
    return flat


def required_properties(flags):
    """List the properties that must be defined for the enabled features.

    Parameters
    ----------
    flags : int or str
        Feature flags.

    Returns
    -------
    list of str
        Property names, without duplicates, in first-required order.
    """
    f = parse_feature_flags(flags)
    groups = ["elastic"]
    if f["matrix_damage"]:
        groups.append("matrix_damage")
    if f["shear_nonlinearity"] in (1, 3, 4):
        groups.append("shear_ro")
    if f["fiber_tension"]:
        groups.append("fiber_tension")
    if f["fiber_compression"] == 1:
        groups.append("fiber_compression_bl")
    elif f["fiber_compression"] in (3, 4, 5):
        groups.append("fiber_compression_fkt")
    if f["friction"]:
        groups.append("friction")
    names = []
    for g in groups:
        for p in FEATURE_REQUIREMENTS[g]:
            if p not in names:
                names.append(p)
    return names


def property_values(material, flags, thickness=None):
    """The 40 CompDam properties in input-deck order, ``None`` for blanks.

    Parameters
    ----------
    material : dict
        Material preset from ``config.MATERIALS``.
    flags : int or str
        Feature flags (material property 1).
    thickness : float, optional
        Material property 3. Leave as ``None`` for solid elements; CompDam
        gets the characteristic lengths from VUCHARLENGTH.

    Returns
    -------
    list
        40 entries; ``None`` marks a blank field (feature disabled).

    Raises
    ------
    ConfigurationError
        If a property needed by an enabled feature is missing or zero, or
        a property is unknown to CompDam.
    """
    flat = flatten_material(material)
    unknown = [k for k in flat if k not in PROPERTY_ORDER]
    if unknown:
        raise ConfigurationError("unknown CompDam properties: %s" % ", ".join(sorted(unknown)))

    missing = [p for p in required_properties(flags) if not flat.get(p)]
    if missing:
        raise ConfigurationError(
            "feature flags %s need these properties, which are missing or zero: %s"
            % (flags, ", ".join(missing)))

    values = dict(flat)
    values["feature_flags"] = int(feature_flags_as_float(flags))
    values["thickness"] = thickness
    return [values.get(name) for name in PROPERTY_ORDER]


def build_property_vector(material, flags, thickness=None):
    """Assemble the 40-entry property vector as floats (blank -> 0.0).

    Parameters
    ----------
    material : dict
        Material preset from ``config.MATERIALS``.
    flags : int or str
        Feature flags (material property 1).
    thickness : float, optional
        Material property 3.

    Returns
    -------
    list of float
        Exactly 40 values; 0.0 is what Abaqus passes to CompDam for a blank.
    """
    return [0.0 if v is None else float(v)
            for v in property_values(material, flags, thickness)]


def format_user_material(material, flags):
    """Data lines of ``*User Material, constants=40`` (8 values per line).

    Parameters
    ----------
    material : dict
        Material preset.
    flags : int or str
        Feature flags.

    Returns
    -------
    list of str
        Five data lines; blank fields are left empty, exactly as in the
        CompDam example decks.
    """
    vals = property_values(material, flags)
    out = []
    for i in range(0, 40, 8):
        chunk = []
        for v in vals[i:i + 8]:
            if v is None:
                chunk.append("")
            elif isinstance(v, int):
                chunk.append("%d" % v)
            else:
                chunk.append("%.10g" % v)
        out.append(", ".join(chunk) + ("," if i < 32 else ""))
    return out


def max_element_size(material):
    """Largest element length that avoids constitutive snap-back per mode.

    For a bilinear/linear softening law the regularised element length must
    satisfy ``Lc < 2 E G / sigma^2``. Uses E2 for the matrix modes and E1 for
    the fibre modes.

    Parameters
    ----------
    material : dict
        Material preset.

    Returns
    -------
    dict
        ``{mode: max_length_mm}`` for each mode with complete data.
    """
    p = flatten_material(material)
    checks = {
        "matrix_tension": ("E2", "GYT", "YT"),
        "matrix_shear": ("G12", "GSL", "SL"),
        "fiber_tension": ("E1", "GXT", "XT"),
        "fiber_compression": ("E1", "GXC", "XC"),
    }
    out = {}
    for mode, (e, g, s) in checks.items():
        if p.get(e) and p.get(g) and p.get(s):
            out[mode] = 2.0 * p[e] * p[g] / p[s] ** 2
    return out


def initial_state(nstatev, d2=0.0, alpha=0, d1T=0.0, d1C=0.0, phi0=0.0):
    """Build the initial state-variable vector for one element set.

    Parameters
    ----------
    nstatev : int
        Number of state variables.
    d2 : float, optional
        Initial matrix damage ``CDM_d2`` (0 pristine, 1 fully cracked).
    alpha : int, optional
        Initial matrix-crack normal angle ``CDM_alpha`` [deg].
    d1T, d1C : float, optional
        Initial fibre tension/compression damage (needs ``nstatev >= 19``
        for ``d1C``).
    phi0 : float, optional
        Initial fibre misalignment for FKT (``nstatev >= 25``).

    Returns
    -------
    list
        Values in state-variable order. ``CDM_STATUS`` is always 1 and
        ``CDM_alpha`` is an integer, as CompDam requires.

    Raises
    ------
    ConfigurationError
        If a damage value lies outside [0, 1] or the state-variable array is
        too short for the requested value.
    """
    if nstatev < 18:
        raise ConfigurationError("CompDam needs at least 18 state variables")
    for name, v in (("d2", d2), ("d1T", d1T), ("d1C", d1C)):
        if not 0.0 <= v <= 1.0:
            raise ConfigurationError("%s must lie in [0, 1], got %r" % (name, v))
    if d1C and nstatev < SDV_D1C:
        raise ConfigurationError("d1C seeding needs fibre compression model 1 (19 SDVs)")

    state = [0.0] * nstatev
    state[SDV_D2 - 1] = float(d2)
    state[SDV_ALPHA - 1] = int(alpha)
    state[SDV_STATUS - 1] = 1
    state[SDV_D1T - 1] = float(d1T)
    if nstatev >= SDV_D1C:
        state[SDV_D1C - 1] = float(d1C)
    if nstatev >= SDV_PHI0_13:
        state[SDV_PHI0_12 - 1] = float(phi0)
        state[SDV_PHI0_13 - 1] = float(phi0)
    return state


def _fmt(value):
    if isinstance(value, int):
        return "%d" % value
    return "%.8g" % value


def format_initial_conditions(sets_and_states):
    """Format an ``*Initial Conditions, type=SOLUTION`` keyword block.

    Parameters
    ----------
    sets_and_states : list of tuple
        ``(element_set_name, state_vector)`` pairs. Set names refer to
        assembly-level element sets.

    Returns
    -------
    str
        Keyword block, terminated with a newline. The first data line of each
        set holds the set name plus 7 values; continuation lines hold 8.
    """
    lines = ["*Initial Conditions, type=SOLUTION"]
    for elset, state in sets_and_states:
        values = [_fmt(v) for v in state]
        lines.append(", ".join([elset] + values[:7]) + ("," if len(values) > 7 else ""))
        rest = values[7:]
        for i in range(0, len(rest), 8):
            chunk = rest[i:i + 8]
            trailing = "," if i + 8 < len(rest) else ""
            lines.append(", ".join(chunk) + trailing)
    return "\n".join(lines) + "\n"


def format_depvar_names(nstatev):
    """Return the ``*Depvar`` name data lines for ``nstatev`` variables.

    Parameters
    ----------
    nstatev : int
        Number of state variables.

    Returns
    -------
    list of str
        Lines of the form ``"  1, CDM_d2"``.
    """
    names = list(STATE_VARIABLE_NAMES[:nstatev])
    names += ["SDV%d" % i for i in range(len(names) + 1, nstatev + 1)]
    return ["%3d, %s" % (i + 1, n) for i, n in enumerate(names)]


def format_parameter_tables(float_params, int_params):
    """Format the COMPDAM_PARAM table collection for non-default parameters.

    Parameters
    ----------
    float_params : dict
        Floating-point model parameters (e.g. ``{"gamma_max": 4.0}``).
    int_params : dict
        Integer model parameters (e.g. ``{"logLevel": 2}``).

    Returns
    -------
    str
        Keyword block, or an empty string when both dictionaries are empty.
    """
    if not float_params and not int_params:
        return ""
    lines = [
        "** CompDam model parameters (non-default values only)",
        "*Parameter Table Type, name=KEY_FLOAT, parameters=2",
        "STRING",
        "FLOAT",
        "*Parameter Table Type, name=KEY_INT, parameters=2",
        "STRING",
        "INTEGER",
        "*Table Collection, name=COMPDAM_PARAM",
    ]
    if float_params:
        lines.append("*Parameter Table, type=KEY_FLOAT, label=CDP_FLOAT")
        for k in sorted(float_params):
            lines.append("%s, %.10g" % (k, float(float_params[k])))
    if int_params:
        lines.append("*Parameter Table, type=KEY_INT, label=CDP_INT")
        for k in sorted(int_params):
            lines.append("%s, %d" % (k, int(int_params[k])))
    return "\n".join(lines) + "\n"
