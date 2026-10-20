"""
In-plane partitioning and edge-seeding plan for the CAI plate.

The plate is cut by planes at

* the end-support (clamp) lines ``x = s`` and ``x = L - s``,
* the side knife-edge lines ``y = inset`` and ``y = W - inset``,
* the edges of the refined zone around the impact site,
* the specimen centre line ``y = y_c`` (gives a node line for the rigid-body
  restraint),

which leaves a grid of hexahedral blocks that can all be meshed with the
STRUCTURED technique. Each edge gets one of three seeds:

* ``fine`` inside the refined zone,
* a single bias from ``fine`` to ``coarse`` in the transition band between
  the refined zone and the support lines,
* ``coarse`` everywhere else (under the clamps and outside the knife edges).

Opposite edges of every block get identical seeds, which keeps the
structured mesh valid. Every ply block gets the same plan, so the in-plane
meshes of all blocks match and their interface nodes coincide.

Pure Python: the Abaqus builder only queries this plan.
"""

from __future__ import division, print_function

from cai import ConfigurationError


def _unique_sorted(values, lo, hi, tol):
    out = []
    for v in sorted(values):
        if lo + tol < v < hi - tol and (not out or v - out[-1] > tol):
            out.append(v)
    return out


class MeshPlan(object):
    """Partition coordinates and seeding rules for the plate.

    Parameters
    ----------
    geometry : dict
        ``config.GEOMETRY``.
    mesh : dict
        ``config.MESH``.
    boundary : dict
        ``config.BOUNDARY``.
    center : tuple of float
        Impact location ``(x, y)``.

    Attributes
    ----------
    x_cuts, y_cuts : list of float
        Interior partition-plane coordinates.
    zone_x, zone_y : tuple of float
        Bounds of the refined zone.
    band_x, band_y : tuple of float
        Outer bounds of the transition band (support lines).
    """

    def __init__(self, geometry, mesh, boundary, center):
        self.L = float(geometry["length"])
        self.W = float(geometry["width"])
        self.fine = float(mesh["fine_size"])
        self.coarse = float(mesh["coarse_size"])
        self.n_through = int(mesh.get("elements_per_ply", 1))
        self.tol = float(mesh.get("tolerance", 1e-4))
        xc, yc = center
        lz = float(geometry["refine_zone_length"])
        wz = float(geometry["refine_zone_width"])
        es = float(boundary["end_support_length"])
        inset = float(boundary["side_support_inset"])

        self.center = (float(xc), float(yc))
        self.zone_x = (xc - lz / 2.0, xc + lz / 2.0)
        self.zone_y = (yc - wz / 2.0, yc + wz / 2.0)
        self.band_x = (es, self.L - es)
        self.band_y = (inset, self.W - inset)
        if not (self.band_x[0] < self.zone_x[0] and self.zone_x[1] < self.band_x[1]):
            raise ConfigurationError("refined zone must lie inside the end supports")
        if not (self.band_y[0] < self.zone_y[0] and self.zone_y[1] < self.band_y[1]):
            raise ConfigurationError("refined zone must lie inside the side supports")

        self.x_cuts = _unique_sorted([es, self.zone_x[0], self.zone_x[1], self.L - es],
                                     0.0, self.L, self.tol)
        self.y_cuts = _unique_sorted([inset, self.zone_y[0], yc, self.zone_y[1], self.W - inset],
                                     0.0, self.W, self.tol)

    def _axis_rule(self, lo, hi, zone, band):
        t = self.tol
        if lo >= zone[0] - t and hi <= zone[1] + t:
            return ("size", self.fine)
        if abs(hi - zone[0]) < t and lo >= band[0] - t:
            return ("bias", self.fine, self.coarse, hi)
        if abs(lo - zone[1]) < t and hi <= band[1] + t:
            return ("bias", self.fine, self.coarse, lo)
        return ("size", self.coarse)

    def edge_seed(self, p0, p1):
        """Seeding rule for a straight, axis-aligned edge.

        Parameters
        ----------
        p0, p1 : sequence of float
            End-point coordinates ``(x, y, z)`` of the edge.

        Returns
        -------
        tuple
            One of:

            * ``("number", n)``: through-thickness edges (``n`` elements),
            * ``("size", h)``: uniform element size ``h``,
            * ``("bias", h_min, h_max, c_fine)``: single bias, smallest
              elements at the end whose axial coordinate is ``c_fine``.

        Raises
        ------
        ConfigurationError
            If the edge is not parallel to a global axis.
        """
        d = [abs(p1[i] - p0[i]) for i in range(3)]
        axis = d.index(max(d))
        others = [d[i] for i in range(3) if i != axis]
        if max(others) > self.tol:
            raise ConfigurationError("edge %s-%s is not axis aligned" % (p0, p1))
        if axis == 2:
            return ("number", self.n_through)
        lo, hi = sorted((p0[axis], p1[axis]))
        if axis == 0:
            return self._axis_rule(lo, hi, self.zone_x, self.band_x)
        return self._axis_rule(lo, hi, self.zone_y, self.band_y)

    def estimated_elements_per_layer(self):
        """Rough element count per through-thickness element layer.

        Returns
        -------
        int
        """
        def count(cuts, length, zone):
            pts = [0.0] + cuts + [length]
            n = 0.0
            for a, b in zip(pts[:-1], pts[1:]):
                rule = self._axis_rule(a, b, zone, (cuts[0], cuts[-1]))
                h = self.fine if rule[0] == "size" and rule[1] == self.fine else (
                    0.5 * (self.fine + self.coarse) if rule[0] == "bias" else self.coarse)
                n += max(1.0, round((b - a) / h))
            return n
        return int(count(self.x_cuts, self.L, self.zone_x) * count(self.y_cuts, self.W, self.zone_y))
