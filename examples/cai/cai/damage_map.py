"""
Geometric description of the pre-defined impact damage state.

The CAI analysis starts from an idealised damage state that represents the
barely-visible impact damage (BVID) left by an ASTM D7136 drop-weight
impact:

* **Delaminations**: one ellipse per ply-block interface, with its major axis
  along the fibre direction of the ply below the interface. The equivalent
  radius grows linearly from the impacted face to the back face (conical
  damage distribution).
* **Matrix cracks**: circular regions per ply, also growing toward the back
  face, seeded with a fully formed CompDam matrix crack.
* **Fibre damage**: optional small region in the plies under the impactor.

Footprints scale with impact energy (see :func:`energy_scale`). Membership
tests take global coordinates, so the same functions serve the Abaqus element
selection and the pure-Python tests.
"""

from __future__ import division, print_function

import math

from cai import ConfigurationError


def energy_scale(energy, reference_energy, threshold_energy=0.0):
    """Radius scale factor for an impact energy.

    Projected damage area is taken to grow linearly with impact energy above
    a threshold, so radii scale with the square root of the energy ratio.

    Parameters
    ----------
    energy : float
        Impact energy [J].
    reference_energy : float
        Energy at which the configured reference radii apply [J].
    threshold_energy : float, optional
        Damage threshold energy [J]; no damage below it.

    Returns
    -------
    float
        Factor applied to the reference radii (0 at or below the threshold).

    Raises
    ------
    ConfigurationError
        If ``reference_energy <= threshold_energy``.
    """
    if reference_energy <= threshold_energy:
        raise ConfigurationError("reference_energy must exceed threshold_energy")
    if energy <= threshold_energy:
        return 0.0
    return math.sqrt((energy - threshold_energy) / (reference_energy - threshold_energy))


def interpolate_radius(z, thickness, radius_front, radius_back):
    """Linearly interpolate a radius through the thickness.

    Parameters
    ----------
    z : float
        Global z-coordinate (0 = back face, ``thickness`` = impacted face).
    thickness : float
        Laminate thickness.
    radius_front, radius_back : float
        Radii at the impacted (z = t) and back (z = 0) faces.

    Returns
    -------
    float
    """
    s = min(max(z / thickness, 0.0), 1.0)
    return radius_back + (radius_front - radius_back) * s


def inside_ellipse(x, y, center, semi_major, semi_minor, angle_deg):
    """Test whether a point lies inside a rotated ellipse.

    Parameters
    ----------
    x, y : float
        Point coordinates.
    center : tuple of float
        Ellipse centre ``(xc, yc)``.
    semi_major, semi_minor : float
        Semi-axes; zero or negative means an empty ellipse.
    angle_deg : float
        Orientation of the major axis from global X [deg].

    Returns
    -------
    bool
    """
    if semi_major <= 0.0 or semi_minor <= 0.0:
        return False
    t = math.radians(angle_deg)
    dx, dy = x - center[0], y - center[1]
    u = dx * math.cos(t) + dy * math.sin(t)
    v = -dx * math.sin(t) + dy * math.cos(t)
    return (u / semi_major) ** 2 + (v / semi_minor) ** 2 <= 1.0


class DamageMap(object):
    """Pre-defined impact damage for a given laminate and impact energy.

    Parameters
    ----------
    laminate : cai.layup.Laminate
        Laminate definition.
    damage : dict
        ``config.DAMAGE``.
    impact : dict
        ``config.IMPACT`` (threshold and reference energies).
    center : tuple of float
        Impact location ``(x, y)`` [mm].
    energy : float
        Impact energy [J].

    Attributes
    ----------
    scale : float
        Radius scale factor from :func:`energy_scale`.
    """

    def __init__(self, laminate, damage, impact, center, energy):
        self.laminate = laminate
        self.cfg = damage
        self.center = (float(center[0]), float(center[1]))
        self.energy = float(energy)
        self.enabled = bool(damage.get("enabled", False))
        self.scale = energy_scale(energy, impact["reference_energy"],
                                  impact.get("threshold_energy", 0.0))
        ar = damage["delamination"].get("aspect_ratio", 1.0)
        if ar < 1.0:
            raise ConfigurationError("delamination aspect_ratio must be >= 1")

    # ------------------------------------------------------------------ delam
    def delamination_axes(self, interface):
        """Semi-axes and orientation of the delamination at an interface.

        Parameters
        ----------
        interface : cai.layup.Interface

        Returns
        -------
        tuple
            ``(semi_major, semi_minor, angle_deg)``; zero axes mean no
            delamination at this interface.
        """
        d = self.cfg["delamination"]
        if not (self.enabled and d.get("enabled", False)) or self.scale <= 0.0:
            return 0.0, 0.0, 0.0
        r = self.scale * interpolate_radius(
            interface.z, self.laminate.thickness, d["radius_front"], d["radius_back"])
        root_ar = math.sqrt(d.get("aspect_ratio", 1.0))
        return r * root_ar, r / root_ar, interface.lower_ply.angle

    def is_delaminated(self, interface, x, y):
        """Return True if point ``(x, y)`` on ``interface`` is pre-delaminated."""
        a, b, ang = self.delamination_axes(interface)
        return inside_ellipse(x, y, self.center, a, b, ang)

    def delamination_area(self, interface):
        """Area of the seeded delamination at an interface [mm^2]."""
        a, b, _ = self.delamination_axes(interface)
        return math.pi * a * b

    def projected_damage_radius(self):
        """Largest seeded delamination semi-major axis over all interfaces [mm]."""
        radii = [self.delamination_axes(i)[0] for i in self.laminate.interfaces]
        return max(radii) if radii else 0.0

    # ----------------------------------------------------------- intralaminar
    def matrix_crack_radius(self, ply):
        """Radius of the seeded matrix-crack region in a ply [mm]."""
        m = self.cfg["matrix_cracking"]
        if not (self.enabled and m.get("enabled", False)) or self.scale <= 0.0:
            return 0.0
        return self.scale * interpolate_radius(
            ply.z_mid, self.laminate.thickness, m["radius_front"], m["radius_back"])

    def has_fiber_damage(self, ply):
        """Return True if fibre damage is seeded in ``ply``."""
        f = self.cfg["fiber_damage"]
        if not (self.enabled and f.get("enabled", False)) or self.scale <= 0.0:
            return False
        n = len(self.laminate.plies)
        return ply.index > n - int(f.get("plies_from_front", 0))

    def classify_element(self, ply, x, y):
        """Damage category of an element from its centroid.

        Parameters
        ----------
        ply : cai.layup.Ply
            Ply that contains the element.
        x, y : float
            Element centroid (global in-plane coordinates).

        Returns
        -------
        str
            ``"fiber"``, ``"matrix"`` or ``"pristine"`` (fibre takes
            precedence; fibre-damaged elements also get the matrix crack).
        """
        r2 = (x - self.center[0]) ** 2 + (y - self.center[1]) ** 2
        if self.has_fiber_damage(ply) and r2 <= self.cfg["fiber_damage"]["radius"] ** 2:
            return "fiber"
        rm = self.matrix_crack_radius(ply)
        if rm > 0.0 and r2 <= rm * rm:
            return "matrix"
        return "pristine"

    def summary(self):
        """Readable description of the seeded damage.

        Returns
        -------
        str
        """
        lines = ["Impact damage: E = %.1f J, radius scale = %.3f, enabled = %s" % (
            self.energy, self.scale, self.enabled)]
        for i in self.laminate.interfaces:
            a, b, ang = self.delamination_axes(i)
            if a > 0.0:
                lines.append("  %s z=%.3f: ellipse %.1f x %.1f mm at %+g deg (%.0f mm^2)" % (
                    i.name, i.z, 2 * a, 2 * b, ang, math.pi * a * b))
        return "\n".join(lines)
