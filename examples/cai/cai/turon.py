"""
Cohesive-zone based mesh sizing and interface properties (Turon et al. 2007).

References
----------
.. [Turon2007] A. Turon, C.G. Davila, P.P. Camanho, J. Costa. An engineering
   solution for mesh size effects in the simulation of delamination using
   cohesive zone models. Eng Fract Mech 74 (2007) 1665-1682.
.. [Harper2008] P.W. Harper, S.R. Hallett. Cohesive zone length in numerical
   simulations of composite delamination. Eng Fract Mech 75 (2008) 4774-4792.

Equations
---------
From [Turon2007]_; the mode II cohesive-zone length may use G12 as in
[Harper2008]_:

* cohesive zone length: ``l_cz = M E Gc / tau0^2`` with ``M = 9 pi / 32``
* element size needed for ``N_e`` elements in the zone: ``l_e = l_cz / N_e``
* strength that gives ``N_e`` elements for a chosen ``l_e``:
  ``tau0' = sqrt(9 pi E Gc / (32 N_e l_e))``
* interface penalty stiffness: ``K = alpha E3 / t`` with ``alpha >= 50``
"""

from __future__ import division, print_function

import math

from cai import ConfigurationError
from cai import compdam

M_RICE = 9.0 * math.pi / 32.0


def cohesive_zone_length(modulus, toughness, strength, M=M_RICE):
    """Length of the cohesive (process) zone.

    Parameters
    ----------
    modulus : float
        Representative modulus E [MPa].
    toughness : float
        Fracture toughness Gc [N/mm].
    strength : float
        Interface strength tau0 [MPa].
    M : float, optional
        Model parameter (9 pi / 32 by default, as in Turon et al.).

    Returns
    -------
    float
        ``M E Gc / tau0^2`` [mm].
    """
    return M * modulus * toughness / strength ** 2


def reduced_strength(modulus, toughness, n_elements, element_size, M=M_RICE):
    """Interface strength for which ``n_elements`` span the cohesive zone.

    Parameters
    ----------
    modulus, toughness : float
        E [MPa] and Gc [N/mm].
    n_elements : float
        Required number of elements in the cohesive zone, N_e.
    element_size : float
        Element length l_e [mm].
    M : float, optional
        Model parameter.

    Returns
    -------
    float
        ``sqrt(M E Gc / (N_e l_e))`` [MPa]; equals Turon's
        ``sqrt(9 pi E Gc / (32 N_e l_e))`` for the default M.
    """
    return math.sqrt(M * modulus * toughness / (n_elements * element_size))


def penalty_stiffness(e3, thickness, alpha=50.0):
    """Turon's interface penalty stiffness ``K = alpha E3 / t``.

    Parameters
    ----------
    e3 : float
        Through-thickness modulus [MPa].
    thickness : float
        Thickness of the adjacent sublaminate [mm] (one ply here).
    alpha : float, optional
        Must be >= 50 for a < 2 % loss of laminate compliance.

    Returns
    -------
    float
        Penalty stiffness [N/mm^3].
    """
    return alpha * e3 / thickness


class CohesiveSizing(object):
    """Damage-zone element size and effective interface properties.

    Parameters
    ----------
    material : dict
        Material preset (``elastic`` and ``interface`` groups are used).
    zone_cfg : dict
        ``config.COHESIVE_ZONE``.
    ply_thickness : float
        Ply thickness [mm] for the penalty stiffness.

    Attributes
    ----------
    l_cz : dict
        Cohesive zone length per mode, ``{"I": ..., "II": ...}`` [mm].
    turon_size : float
        ``min(l_cz) / N_e``: the element size Turon's criterion asks for [mm].
    element_size : float
        Element size actually used in the damage zone [mm].
    strength_normal, strength_shear : float
        Interface strengths written to the deck (possibly reduced) [MPa].
    strength_reduced : bool
        True when Turon's strength reduction was applied.
    penalty : float
        Interface penalty stiffness [N/mm^3].
    notes : list of str
        Human-readable explanation of the decisions taken.
    """

    def __init__(self, material, zone_cfg, ply_thickness):
        flat = compdam.flatten_material(material)
        itf = material.get("interface")
        if not itf:
            raise ConfigurationError("material preset lacks 'interface' properties")
        e3 = flat.get("E3") or flat["E2"]
        g13 = flat.get("G13") or flat["G12"]
        mode2 = zone_cfg.get("mode_II_modulus", "G12").upper()
        if mode2 not in ("G12", "E3"):
            raise ConfigurationError("mode_II_modulus must be 'G12' or 'E3'")
        e_mode2 = g13 if mode2 == "G12" else e3
        M = float(zone_cfg.get("M", M_RICE))
        ne = float(zone_cfg.get("elements_in_zone", 3))
        if ne < 1:
            raise ConfigurationError("elements_in_zone must be >= 1")

        self.M = M
        self.moduli = {"I": e3, "II": e_mode2}
        self.toughness = {"I": itf["GIc"], "II": itf["GIIc"]}
        self.nominal = {"I": itf["strength_normal"], "II": itf["strength_shear"]}
        self.l_cz = dict((m, cohesive_zone_length(self.moduli[m], self.toughness[m],
                                                  self.nominal[m], M)) for m in ("I", "II"))
        self.turon_size = min(self.l_cz.values()) / ne
        self.notes = []

        size = self.turon_size
        lo = zone_cfg.get("min_element_size")
        if lo and size < lo:
            self.notes.append("Turon size %.3f mm is below min_element_size; using %.3f mm"
                              % (size, lo))
            size = float(lo)
        hi = zone_cfg.get("max_element_size")
        if hi is None:
            limits = compdam.max_element_size(material)
            matrix = [limits[k] for k in ("matrix_tension", "matrix_shear") if k in limits]
            hi = min(matrix) if matrix else None
        if hi and size > hi:
            self.notes.append("capped at %.3f mm (max_element_size / CompDam snap-back limit)"
                              % hi)
            size = float(hi)
        self.element_size = size

        strengths = dict(self.nominal)
        self.strength_reduced = False
        if zone_cfg.get("reduce_strength", True):
            for m in ("I", "II"):
                tau = reduced_strength(self.moduli[m], self.toughness[m], ne, size, M)
                if tau < strengths[m]:
                    strengths[m] = tau
                    self.strength_reduced = True
            if self.strength_reduced:
                self.notes.append(
                    "interface strengths reduced (Turon) to %.1f / %.1f MPa so that %d "
                    "elements span the cohesive zone" % (strengths["I"], strengths["II"], ne))
        elif size > self.turon_size:
            self.notes.append("WARNING: only %.1f elements in the cohesive zone (target %d)"
                              % (min(self.l_cz.values()) / size, ne))
        self.strength_normal = strengths["I"]
        self.strength_shear = strengths["II"]

        pen = zone_cfg.get("penalty", "turon")
        if isinstance(pen, str):
            if pen.lower() != "turon":
                raise ConfigurationError("penalty must be 'turon' or a number")
            self.penalty = penalty_stiffness(e3, ply_thickness,
                                             float(zone_cfg.get("penalty_alpha", 50.0)))
        else:
            self.penalty = float(pen)

    def elements_in_zone(self):
        """Number of elements spanning the (possibly enlarged) cohesive zone.

        Returns
        -------
        float
            ``min(l_cz') / l_e`` using the strengths actually written.
        """
        M = self.M
        lI = cohesive_zone_length(self.moduli["I"], self.toughness["I"], self.strength_normal, M)
        lII = cohesive_zone_length(self.moduli["II"], self.toughness["II"],
                                   self.strength_shear, M)
        return min(lI, lII) / self.element_size

    def summary(self):
        """Readable multi-line summary.

        Returns
        -------
        str
        """
        lines = [
            "Cohesive zone (Turon et al. 2007): l_cz,I = %.3f mm, l_cz,II = %.3f mm"
            % (self.l_cz["I"], self.l_cz["II"]),
            "  Turon element size %.3f mm -> damage-zone element size %.3f mm"
            % (self.turon_size, self.element_size),
            "  interface strengths %.1f / %.1f MPa (nominal %.1f / %.1f), "
            "%.1f elements in zone, K = %.3g N/mm^3"
            % (self.strength_normal, self.strength_shear, self.nominal["I"],
               self.nominal["II"], self.elements_in_zone(), self.penalty),
        ]
        lines += ["  " + n for n in self.notes]
        return "\n".join(lines)
