"""
Fibre-aligned structured quadrilateral meshing of one ply (block) of the plate.

Concept
-------
For a ply with fibre angle ``theta`` the mesh is generated in a working frame
in which the fibres are "steep" (``|theta_g| >= 45 deg``). For
``|theta| < 45 deg`` (e.g. 0 deg plies) the plate is transposed (x <-> y)
first, so a single algorithm covers every angle.

In the working frame the skewed coordinate ``s = x - c y`` with
``c = cot(theta_g)`` is constant along a fibre. The **damage zone** is the
parallelogram ``s0 <= s <= s1, y0 <= y <= y1``, with two edges along the
fibres and two edges horizontal. It is meshed with a uniform grid of
spacing ``h`` along s and ``h |sin theta_g|`` along y, so every element edge
in the zone is either parallel to the fibres (length h) or horizontal
(length h): the mesh is **fibre-aligned**. For 0/90 deg plies the
parallelogram degenerates to a rectangle.

The rest of the plate is filled without breaking the structured
``Nx x Ny`` grid:

* rows in the zone's y-range continue left and right of the zone, graded
  geometrically from ``h`` up to the coarse size, with a node forced on the
  support line (end clamp or knife edge),
* below and above the zone, rows are graded in y (with a row forced on the
  other support line) and the columns are vertical continuations of the
  zone's bottom/top row.

All quads are counter-clockwise in the original (x, y) frame.

Pure Python: no Abaqus, no numpy.
"""

from __future__ import division, print_function

import math

from cai import ConfigurationError


# -----------------------------------------------------------------------------
# 1-D spacing helpers
# -----------------------------------------------------------------------------
def normalize_angle(angle):
    """Map an angle in degrees to the interval (-90, 90].

    Parameters
    ----------
    angle : float
        Angle [deg].

    Returns
    -------
    float
        Equivalent fibre direction angle in (-90, 90].
    """
    a = math.fmod(float(angle), 180.0)
    if a <= -90.0:
        a += 180.0
    elif a > 90.0:
        a -= 180.0
    return a


def count_cells(length, first, max_size, growth):
    """Number of cells to cover ``length`` growing from ``first`` to ``max_size``.

    Parameters
    ----------
    length : float
        Segment length.
    first : float
        Size of the first cell.
    max_size : float
        Size at which growth stops.
    growth : float
        Ratio between consecutive cells (> 1).

    Returns
    -------
    int
        Number of cells (0 for a zero-length segment).
    """
    if length <= 0.0:
        return 0
    if first >= length:
        return 1
    n, total, size = 0, 0.0, first
    while total < length * (1.0 - 1e-9):
        total += size
        n += 1
        size = min(size * growth, max_size)
    return n


def graded_sizes(length, n, first, max_size=None):
    """Split ``length`` into ``n`` cells growing geometrically from ``first``.

    Parameters
    ----------
    length : float
        Segment length.
    n : int
        Number of cells.
    first : float
        Desired size of the first cell.
    max_size : float, optional
        Cells stop growing at this size (the growth ratio is solved so that
        the capped sizes still sum to ``length``). Ignored if it cannot be
        honoured with ``n`` cells.

    Returns
    -------
    list of float
        ``n`` sizes summing to ``length``. Uniform when ``n * first >= length``.
    """
    if n <= 0:
        return []
    if length <= n * first * (1.0 + 1e-12):
        return [length / n] * n
    cap = max_size if (max_size and first + (n - 1) * max_size >= length) else float("inf")

    def sizes_for(r):
        out, size = [], first
        for _ in range(n):
            out.append(min(size, cap))
            size *= r
        return out

    lo, hi = 1.0, 2.0
    while sum(sizes_for(hi)) < length and hi < 1e3:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if sum(sizes_for(mid)) < length:
            lo = mid
        else:
            hi = mid
    sizes = sizes_for(0.5 * (lo + hi))
    scale = length / sum(sizes)
    return [x * scale for x in sizes]


def segment(start, end, n, first=None, fine_at="start", max_size=None):
    """Node positions on a segment (both ends included).

    Parameters
    ----------
    start, end : float
        Segment end points (``end > start``).
    n : int
        Number of cells.
    first : float, optional
        Size of the finest cell; uniform spacing when None.
    fine_at : {"start", "end"}
        Which end carries the finest cell.
    max_size : float, optional
        Largest cell size of the grading.

    Returns
    -------
    list of float
        ``n + 1`` increasing positions.
    """
    length = end - start
    if n <= 0:
        return [start] if length <= 0 else [start, end]
    sizes = [length / n] * n if first is None else graded_sizes(length, n, first, max_size)
    if fine_at == "end":
        sizes = sizes[::-1]
    pos = [start]
    for s in sizes:
        pos.append(pos[-1] + s)
    pos[-1] = end
    return pos


def join(*parts):
    """Concatenate position lists, dropping the duplicated junction points."""
    out = []
    for p in parts:
        if not p:
            continue
        if out and abs(out[-1] - p[0]) < 1e-9:
            out.extend(p[1:])
        else:
            out.extend(p)
    return out


# -----------------------------------------------------------------------------
# Working frame and damage zone
# -----------------------------------------------------------------------------
class Shape(object):
    """Ellipse (or circle) the damage zone of a ply must contain.

    Parameters
    ----------
    center : tuple of float
        Centre (x, y) in the original frame.
    semi_major, semi_minor : float
        Semi-axes [mm].
    angle : float
        Major-axis angle from global X [deg].
    """

    def __init__(self, center, semi_major, semi_minor=None, angle=0.0):
        self.center = (float(center[0]), float(center[1]))
        self.a = float(semi_major)
        self.b = float(semi_major if semi_minor is None else semi_minor)
        self.angle = float(angle)


class PlyFrame(object):
    """Working frame for meshing a ply of fibre angle ``angle``.

    Parameters
    ----------
    length, width : float
        Plate dimensions in the original frame.
    angle : float
        Fibre (mesh) angle [deg].
    x_lines : tuple of float
        Original-frame x-coordinates of the support lines at the ends. The
        damage zone stays between them.
    y_lines : tuple of float
        Original-frame y-coordinates of the side support lines.
    x_required, y_required : bool
        Whether the x / y support lines must be exact node lines. Forcing a
        line that crosses the slanted rows of an off-axis zone creates sliver
        elements, so only the lines that need exact nodes (knife edges)
        should be required.

    Attributes
    ----------
    transposed : bool
        True when the working frame swaps x and y.
    Lg, Wg : float
        Plate dimensions in the working frame.
    angle_g : float
        Fibre angle in the working frame, ``|angle_g| >= 45``.
    c : float
        ``cot(angle_g)``.
    sin : float
        ``|sin(angle_g)|``.
    xs, ys : tuple of float
        Support lines in the working frame (columns, rows).
    force_xs, force_ys : bool
        Whether those lines are exact node lines in the working frame.
    """

    def __init__(self, length, width, angle, x_lines, y_lines, x_required=False,
                 y_required=True):
        a = normalize_angle(angle)
        self.angle = a
        self.transposed = abs(a) < 45.0
        if self.transposed:
            self.Lg, self.Wg = float(width), float(length)
            self.angle_g = normalize_angle(90.0 - a)
            self.xs, self.ys = tuple(y_lines), tuple(x_lines)
            self.force_xs, self.force_ys = y_required, x_required
        else:
            self.Lg, self.Wg = float(length), float(width)
            self.angle_g = a
            self.xs, self.ys = tuple(x_lines), tuple(y_lines)
            self.force_xs, self.force_ys = x_required, y_required
        t = math.radians(self.angle_g)
        self.sin = abs(math.sin(t))
        self.c = math.cos(t) / math.sin(t)
        if abs(self.c) < 1e-12:
            self.c = 0.0

    def to_gen(self, x, y):
        """Original -> working frame."""
        return (y, x) if self.transposed else (x, y)

    def to_orig(self, x, y):
        """Working -> original frame."""
        return (y, x) if self.transposed else (x, y)

    def angle_to_gen(self, angle):
        """Direction angle original -> working frame [deg]."""
        return 90.0 - angle if self.transposed else angle


class Zone(object):
    """Fibre-aligned damage zone in the working frame.

    Attributes
    ----------
    s0, y0 : float
        Lower-left corner in skewed coordinates.
    ns, ny : int
        Number of cells along s and y.
    h, dy : float
        Cell sizes along s (horizontal) and y.
    clipped : list of str
        Sides that had to be clipped to stay inside the supports.
    """

    def __init__(self, s0, y0, ns, ny, h, dy, clipped=None):
        self.s0, self.y0, self.ns, self.ny, self.h, self.dy = s0, y0, ns, ny, h, dy
        self.clipped = clipped or []

    @property
    def s1(self):
        return self.s0 + self.ns * self.h

    @property
    def y1(self):
        return self.y0 + self.ny * self.dy


def _extent(cx, cy, a, b, phi, gx, gy):
    """Centre value and half-range of the linear map g.(x, y) over an ellipse."""
    u = (math.cos(phi), math.sin(phi))
    v = (-math.sin(phi), math.cos(phi))
    gu = gx * u[0] + gy * u[1]
    gv = gx * v[0] + gy * v[1]
    return gx * cx + gy * cy, math.sqrt((a * gu) ** 2 + (b * gv) ** 2)


def _fit(lo, hi, cell, lim_lo, lim_hi, name, clipped):
    """Snap [lo, hi] to whole cells, centred, and keep it inside the limits."""
    if lo < lim_lo:
        lo = lim_lo
        clipped.append(name + "_low")
    if hi > lim_hi:
        hi = lim_hi
        clipped.append(name + "_high")
    n = max(1, int(math.ceil((hi - lo) / cell - 1e-9)))
    start = max(0.5 * (lo + hi) - 0.5 * n * cell, lim_lo)
    if start + n * cell > lim_hi + 1e-9:
        n = int(math.floor((lim_hi - start) / cell + 1e-9))
    return start, n


def compute_zone(frame, shapes, h, margin):
    """Smallest fibre-aligned zone containing ``shapes`` plus ``margin``.

    Parameters
    ----------
    frame : PlyFrame
        Working frame of the ply.
    shapes : list of Shape
        Regions to contain (original frame).
    h : float
        Element size in the zone.
    margin : float
        Extra distance around the shapes [mm], measured perpendicular to the
        zone edges.

    Returns
    -------
    Zone
        Zone snapped to the cell grid and clipped to the support lines.

    Raises
    ------
    ConfigurationError
        If no shapes are given or the clipped zone is empty.
    """
    if not shapes:
        raise ConfigurationError("a damage zone needs at least one shape")
    c, sn = frame.c, frame.sin
    dy = h * sn
    s_lo = y_lo = float("inf")
    s_hi = y_hi = -float("inf")
    for sh in shapes:
        cx, cy = frame.to_gen(*sh.center)
        phi = math.radians(frame.angle_to_gen(sh.angle))
        cs, hs = _extent(cx, cy, sh.a, sh.b, phi, 1.0, -c)
        cyy, hy = _extent(cx, cy, sh.a, sh.b, phi, 0.0, 1.0)
        s_lo = min(s_lo, cs - hs - margin / sn)
        s_hi = max(s_hi, cs + hs + margin / sn)
        y_lo = min(y_lo, cyy - hy - margin)
        y_hi = max(y_hi, cyy + hy + margin)

    clipped = []
    p, q = frame.ys
    # Rows: keep at least one cell between the zone and the row lines/edges.
    y0, ny = _fit(y_lo, y_hi, dy, p + dy, q - dy, "y", clipped)
    y1 = y0 + ny * dy
    # Columns: the zone edges are slanted, so check both ends of the y-range.
    a, b = frame.xs
    s_min = a + h - min(c * y0, c * y1)
    s_max = b - h - max(c * y0, c * y1)
    s0, ns = _fit(s_lo, s_hi, h, s_min, s_max, "x", clipped)
    if ns < 1 or ny < 1:
        raise ConfigurationError("damage zone of the %+g deg ply does not fit between "
                                 "the supports" % frame.angle)
    return Zone(s0, y0, ns, ny, h, dy, clipped)


# -----------------------------------------------------------------------------
# Quad mesh
# -----------------------------------------------------------------------------
class QuadMesh(object):
    """Structured quadrilateral mesh of the plate for one ply angle.

    Attributes
    ----------
    nodes : list of tuple
        Node coordinates ``(x, y)`` in the original frame.
    quads : list of tuple
        Node indices (0-based), counter-clockwise in the original frame.
    nx, ny : int
        Grid dimensions (nodes per row, rows).
    angle : float
        Mesh (fibre) angle [deg].
    zone : Zone
        Fibre-aligned damage zone (working frame).
    frame : PlyFrame
    zone_quads : set of int
        Indices of quads inside the fibre-aligned zone.
    """

    def __init__(self, nodes, quads, nx, ny, frame, zone, zone_quads):
        self.nodes = nodes
        self.quads = quads
        self.nx = nx
        self.ny = ny
        self.frame = frame
        self.angle = frame.angle
        self.zone = zone
        self.zone_quads = zone_quads

    def centroid(self, q):
        """Centroid (x, y) of quad ``q``."""
        pts = [self.nodes[i] for i in self.quads[q]]
        return (sum(p[0] for p in pts) / 4.0, sum(p[1] for p in pts) / 4.0)

    def quality(self):
        """Mesh quality indicators.

        Returns
        -------
        dict
            ``min_jacobian`` (smallest corner cross product, must be > 0),
            ``min_angle`` [deg], ``max_aspect`` (longest/shortest edge),
            ``area`` (total, should equal L*W).
        """
        min_j, min_ang, max_ar, area = float("inf"), 180.0, 1.0, 0.0
        for quad in self.quads:
            p = [self.nodes[i] for i in quad]
            edges = []
            for k in range(4):
                x0, y0 = p[k]
                x1, y1 = p[(k + 1) % 4]
                edges.append((x1 - x0, y1 - y0))
            lengths = [math.hypot(*e) for e in edges]
            max_ar = max(max_ar, max(lengths) / min(lengths))
            for k in range(4):
                e_in = edges[k - 1]
                e_out = edges[k]
                cross = e_in[0] * e_out[1] - e_in[1] * e_out[0]
                min_j = min(min_j, cross)
                dot = -(e_in[0] * e_out[0] + e_in[1] * e_out[1])
                ang = math.degrees(math.atan2(abs(cross), dot))
                min_ang = min(min_ang, ang)
            area += 0.5 * ((p[0][0] - p[2][0]) * (p[1][1] - p[3][1])
                           - (p[1][0] - p[3][0]) * (p[0][1] - p[2][1]))
        return {"min_jacobian": min_j, "min_angle": min_ang, "max_aspect": max_ar,
                "area": area}


def mesh_ply(length, width, angle, shapes, h, coarse, growth, x_lines, y_lines, margin,
             x_required=False, y_required=True):
    """Generate the fibre-aligned structured mesh of one ply block.

    Parameters
    ----------
    length, width : float
        Plate dimensions [mm].
    angle : float
        Mesh (fibre) angle [deg]; 0 gives an orthogonal mesh.
    shapes : list of Shape
        Regions the fine zone must contain (seeded damage).
    h : float
        Element size inside the damage zone [mm].
    coarse : float
        Target element size far from the zone [mm].
    growth : float
        Geometric growth ratio of the grading.
    x_lines : tuple of float
        ``(x_a, x_b)``: end-support lines; the damage zone stays between them.
    y_lines : tuple of float
        ``(y_p, y_q)``: side-support (knife-edge) lines.
    margin : float
        Margin around the shapes [mm].
    x_required, y_required : bool, optional
        Whether ``x_lines`` / ``y_lines`` must be exact node lines
        (default: knife edges yes, end supports no).

    Returns
    -------
    QuadMesh
    """
    frame = PlyFrame(length, width, angle, x_lines, y_lines, x_required, y_required)
    zone = compute_zone(frame, shapes, h, margin)
    c, Lg, Wg = frame.c, frame.Lg, frame.Wg
    a, b = frame.xs
    p, q = frame.ys
    dy = zone.dy

    def uniform(start, end):
        n = int(math.ceil((end - start) / coarse - 1e-9)) if end - start > 1e-9 else 0
        return segment(start, end, n) if n else []

    # Without a required line the graded part runs all the way to the edge.
    if not frame.force_ys:
        p, q = 0.0, Wg
    if not frame.force_xs:
        a, b = 0.0, Lg

    # Rows (y positions in the working frame)
    nb = count_cells(zone.y0 - p, dy, coarse, growth)
    nt = count_cells(q - zone.y1, dy, coarse, growth)
    y_mid = [zone.y0 + j * dy for j in range(zone.ny + 1)]
    below = join(uniform(0.0, p), segment(p, zone.y0, nb, dy, "end", coarse))
    Y = join(below, y_mid, segment(zone.y1, q, nt, dy, "start", coarse), uniform(q, Wg))
    k0 = len(below) - 1  # row index of the zone's bottom edge
    k1 = k0 + zone.ny

    # Columns: common cell counts for the graded parts of all middle rows
    left_len = [zone.s0 + c * y - a for y in y_mid]
    right_len = [b - (zone.s1 + c * y) for y in y_mid]
    if min(left_len) <= 0.0 or min(right_len) <= 0.0:
        raise ConfigurationError("damage zone of the %+g deg ply crosses a support line"
                                 % frame.angle)
    nl = count_cells(max(left_len), h, coarse, growth)
    nr = count_cells(max(right_len), h, coarse, growth)

    rows = []
    i_zone0 = None
    for y in y_mid:
        xl = zone.s0 + c * y
        xr = zone.s1 + c * y
        left = join(uniform(0.0, a), segment(a, xl, nl, h, "end", coarse))
        i_zone0 = len(left) - 1  # column index of the zone's left edge
        rows.append(join(left,
                         [xl + i * h for i in range(zone.ns + 1)],
                         segment(xr, b, nr, h, "start", coarse),
                         uniform(b, Lg)))
    nx = len(rows[0])
    if any(len(r) != nx for r in rows):
        raise ConfigurationError("inconsistent row lengths while meshing the %+g deg ply"
                                 % frame.angle)
    X = [rows[0]] * k0 + rows + [rows[-1]] * (len(Y) - k1 - 1)

    nodes = []
    for k, y in enumerate(Y):
        for x in X[k]:
            nodes.append(frame.to_orig(x, y))

    quads, zone_quads = [], set()
    for k in range(len(Y) - 1):
        for i in range(nx - 1):
            n00 = k * nx + i
            n10 = n00 + 1
            n01 = n00 + nx
            n11 = n01 + 1
            if frame.transposed:  # the x<->y swap mirrors the orientation
                quads.append((n00, n01, n11, n10))
            else:
                quads.append((n00, n10, n11, n01))
            if k0 <= k < k1 and i_zone0 <= i < i_zone0 + zone.ns:
                zone_quads.add(len(quads) - 1)
    return QuadMesh(nodes, quads, nx, len(Y), frame, zone, zone_quads)
