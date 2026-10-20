"""
Laminate stacking, ply-block grouping and classical lamination theory checks.

Pure Python (no Abaqus import) so you can unit test it and reuse it in pre- and
post-processing.
"""

from __future__ import division, print_function

import math

from cai import ConfigurationError


class Ply(object):
    """A single ply in the laminate.

    Parameters
    ----------
    index : int
        1-based ply number counted from the back face (z = 0).
    angle : float
        Fibre angle in degrees from global X, positive about +Z.
    z_bottom, z_top : float
        Global through-thickness coordinates of the ply faces [mm].
    block : int
        1-based index of the ply block the ply belongs to.

    Attributes
    ----------
    thickness : float
        ``z_top - z_bottom``.
    z_mid : float
        Mid-plane coordinate of the ply.
    name : str
        Set name used in Abaqus, e.g. ``"PLY_07"``.
    """

    def __init__(self, index, angle, z_bottom, z_top, block):
        self.index = int(index)
        self.angle = float(angle)
        self.z_bottom = float(z_bottom)
        self.z_top = float(z_top)
        self.block = int(block)

    @property
    def thickness(self):
        return self.z_top - self.z_bottom

    @property
    def z_mid(self):
        return 0.5 * (self.z_bottom + self.z_top)

    @property
    def name(self):
        return "PLY_%02d" % self.index

    def fiber_direction(self):
        """Unit vector of the fibre direction in the global X-Y plane.

        Returns
        -------
        tuple of float
            ``(cos(theta), sin(theta), 0.0)``.
        """
        t = math.radians(self.angle)
        return (math.cos(t), math.sin(t), 0.0)

    def __repr__(self):
        return "Ply(%d, %+.1f deg, z=[%.4f, %.4f], block %d)" % (
            self.index, self.angle, self.z_bottom, self.z_top, self.block)


class PlyBlock(object):
    """Group of consecutive plies meshed as one Abaqus part.

    Parameters
    ----------
    index : int
        1-based block number from the back face.
    plies : list of Ply
        Plies in the block, ordered bottom to top.

    Attributes
    ----------
    z_bottom, z_top, thickness : float
        Through-thickness extent of the block.
    name : str
        Part name, e.g. ``"BLOCK_03"``.
    """

    def __init__(self, index, plies):
        if not plies:
            raise ConfigurationError("a ply block needs at least one ply")
        self.index = int(index)
        self.plies = list(plies)

    @property
    def z_bottom(self):
        return self.plies[0].z_bottom

    @property
    def z_top(self):
        return self.plies[-1].z_top

    @property
    def thickness(self):
        return self.z_top - self.z_bottom

    @property
    def name(self):
        return "BLOCK_%02d" % self.index

    @property
    def instance_name(self):
        return self.name + "-1"

    def __repr__(self):
        return "PlyBlock(%d, plies %d-%d)" % (
            self.index, self.plies[0].index, self.plies[-1].index)


class Interface(object):
    """Interface between two adjacent ply blocks.

    Parameters
    ----------
    index : int
        1-based interface number from the back face.
    lower, upper : PlyBlock
        Blocks below and above the interface.

    Attributes
    ----------
    z : float
        Global z-coordinate of the interface.
    lower_ply, upper_ply : Ply
        Plies touching the interface.
    name : str
        e.g. ``"INT_05"``.
    """

    def __init__(self, index, lower, upper):
        self.index = int(index)
        self.lower = lower
        self.upper = upper

    @property
    def z(self):
        return self.lower.z_top

    @property
    def lower_ply(self):
        return self.lower.plies[-1]

    @property
    def upper_ply(self):
        return self.upper.plies[0]

    @property
    def name(self):
        return "INT_%02d" % self.index

    @property
    def angle_mismatch(self):
        """Absolute fibre-angle difference across the interface [deg], in [0, 90]."""
        d = abs(self.upper_ply.angle - self.lower_ply.angle) % 180.0
        return min(d, 180.0 - d)


def expand_sequence(sequence, symmetric=False):
    """Expand a (half-)stacking sequence.

    Parameters
    ----------
    sequence : list of float
        Ply angles from the back face upward.
    symmetric : bool, optional
        If True, append the mirrored sequence.

    Returns
    -------
    list of float
        Full stacking sequence.

    Raises
    ------
    ConfigurationError
        If the sequence is empty.
    """
    angles = [float(a) for a in sequence]
    if not angles:
        raise ConfigurationError("the stacking sequence is empty")
    if symmetric:
        angles = angles + angles[::-1]
    return angles


class Laminate(object):
    """Stacking sequence with ply coordinates and ply-block grouping.

    Parameters
    ----------
    angles : list of float
        Full stacking sequence from the back face (z = 0) upward [deg].
    ply_thickness : float
        Cured ply thickness [mm].
    plies_per_block : int, optional
        Number of consecutive plies per block (ignored if ``blocks`` given).
    blocks : list of int, optional
        Explicit ply counts per block, bottom to top. Must sum to the number
        of plies.

    Attributes
    ----------
    plies : list of Ply
    blocks : list of PlyBlock
    interfaces : list of Interface
    thickness : float
    """

    def __init__(self, angles, ply_thickness, plies_per_block=1, blocks=None):
        if ply_thickness is None or ply_thickness <= 0.0:
            raise ConfigurationError("ply thickness must be positive")
        n = len(angles)
        if blocks is None:
            k = int(plies_per_block)
            if k < 1:
                raise ConfigurationError("plies_per_block must be >= 1")
            blocks = [k] * (n // k) + ([n % k] if n % k else [])
        blocks = [int(b) for b in blocks]
        if any(b < 1 for b in blocks) or sum(blocks) != n:
            raise ConfigurationError(
                "block sizes %s must be positive and sum to %d plies" % (blocks, n))

        self.ply_thickness = float(ply_thickness)
        self.plies = []
        self.blocks = []
        ply_index = 0
        for b_idx, size in enumerate(blocks, start=1):
            block_plies = []
            for _ in range(size):
                z0 = ply_index * self.ply_thickness
                ply = Ply(ply_index + 1, angles[ply_index], z0, z0 + self.ply_thickness, b_idx)
                self.plies.append(ply)
                block_plies.append(ply)
                ply_index += 1
            self.blocks.append(PlyBlock(b_idx, block_plies))
        self.interfaces = [
            Interface(i + 1, self.blocks[i], self.blocks[i + 1])
            for i in range(len(self.blocks) - 1)
        ]

    @classmethod
    def from_config(cls, layup, material):
        """Build a laminate from ``config.LAYUP`` and a material preset.

        Parameters
        ----------
        layup : dict
            ``config.LAYUP``.
        material : dict
            Material preset (provides the nominal ply thickness).

        Returns
        -------
        Laminate
        """
        angles = expand_sequence(layup["sequence"], layup.get("symmetric", False))
        t = layup.get("ply_thickness") or material.get("ply_thickness")
        return cls(angles, t, layup.get("plies_per_block", 1), layup.get("blocks"))

    @property
    def thickness(self):
        return self.plies[-1].z_top

    @property
    def angles(self):
        return [p.angle for p in self.plies]

    def is_symmetric(self, tol=1e-9):
        """Return True if the stacking sequence is mirror-symmetric."""
        a = self.angles
        return all(abs(x - y) < tol for x, y in zip(a, a[::-1]))

    def is_balanced(self, tol=1e-9):
        """Return True if every off-axis ply has a matching ply of opposite sign."""
        off = sorted(((a + 90.0) % 180.0) - 90.0 for a in self.angles
                     if abs(math.sin(math.radians(2 * a))) > tol)
        return all(abs(x + y) < tol for x, y in zip(off, off[::-1]))

    def summary(self):
        """Readable multi-line description of the stack.

        Returns
        -------
        str
        """
        lines = ["Laminate: %d plies, t = %.3f mm, %d blocks, %d interfaces" % (
            len(self.plies), self.thickness, len(self.blocks), len(self.interfaces))]
        for b in self.blocks:
            lines.append("  %s: %s" % (b.name, ", ".join(
                "%s(%+g)" % (p.name, p.angle) for p in b.plies)))
        return "\n".join(lines)


# -----------------------------------------------------------------------------
# Classical lamination theory (in-plane membrane stiffness only)
# -----------------------------------------------------------------------------
def _reduced_stiffness(E1, E2, G12, nu12):
    nu21 = nu12 * E2 / E1
    d = 1.0 - nu12 * nu21
    return (E1 / d, nu12 * E2 / d, E2 / d, G12)


def _transformed_q(q, angle_deg):
    Q11, Q12, Q22, Q66 = q
    t = math.radians(angle_deg)
    c, s = math.cos(t), math.sin(t)
    c2, s2 = c * c, s * s
    Qb11 = Q11 * c2 * c2 + 2 * (Q12 + 2 * Q66) * s2 * c2 + Q22 * s2 * s2
    Qb22 = Q11 * s2 * s2 + 2 * (Q12 + 2 * Q66) * s2 * c2 + Q22 * c2 * c2
    Qb12 = (Q11 + Q22 - 4 * Q66) * s2 * c2 + Q12 * (s2 * s2 + c2 * c2)
    Qb66 = (Q11 + Q22 - 2 * Q12 - 2 * Q66) * s2 * c2 + Q66 * (s2 * s2 + c2 * c2)
    Qb16 = (Q11 - Q12 - 2 * Q66) * s * c2 * c + (Q12 - Q22 + 2 * Q66) * s2 * s * c
    Qb26 = (Q11 - Q12 - 2 * Q66) * s2 * s * c + (Q12 - Q22 + 2 * Q66) * s * c2 * c
    return ((Qb11, Qb12, Qb16), (Qb12, Qb22, Qb26), (Qb16, Qb26, Qb66))


def _inv3(m):
    a, b, c = m[0]
    d, e, f = m[1]
    g, h, i = m[2]
    det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
    if abs(det) < 1e-30:
        raise ConfigurationError("singular laminate stiffness matrix")
    return (
        ((e * i - f * h) / det, (c * h - b * i) / det, (b * f - c * e) / det),
        ((f * g - d * i) / det, (a * i - c * g) / det, (c * d - a * f) / det),
        ((d * h - e * g) / det, (b * g - a * h) / det, (a * e - b * d) / det),
    )


def membrane_properties(laminate, elastic):
    """Effective in-plane laminate engineering constants from CLT.

    Parameters
    ----------
    laminate : Laminate
        Laminate definition.
    elastic : dict
        ``E1``, ``E2``, ``G12``, ``nu12`` of the ply [MPa].

    Returns
    -------
    dict
        ``Ex``, ``Ey``, ``Gxy``, ``nuxy`` [MPa] and the membrane stiffness
        matrix ``A`` [N/mm].

    Notes
    -----
    Use ``Ex`` to check the initial slope of the CAI load-displacement curve:
    ``K = Ex * width * thickness / free_length``.
    """
    q = _reduced_stiffness(elastic["E1"], elastic["E2"], elastic["G12"], elastic["nu12"])
    A = [[0.0] * 3 for _ in range(3)]
    for ply in laminate.plies:
        qb = _transformed_q(q, ply.angle)
        for i in range(3):
            for j in range(3):
                A[i][j] += qb[i][j] * ply.thickness
    a = _inv3(A)
    h = laminate.thickness
    return {
        "Ex": 1.0 / (h * a[0][0]),
        "Ey": 1.0 / (h * a[1][1]),
        "Gxy": 1.0 / (h * a[2][2]),
        "nuxy": -a[0][1] / a[0][0],
        "A": tuple(tuple(r) for r in A),
    }
