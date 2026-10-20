"""
3-D CAI finite element model: extruded ply-block meshes, numbering, sets and surfaces.

Each ply block gets its own fibre-aligned quad mesh (:mod:`cai.mesh2d`),
extruded through the block thickness with ``elements_per_ply`` C3D8R layers
per ply, so element faces coincide with every ply interface. The blocks do
not share nodes; general contact (cohesive or tie) joins them.

Numbering is block by block and layer by layer: block ``b`` owns the
contiguous node range ``node_offset + 1 ...`` and element range
``elem_offset + 1 ...``, with element label
``elem_offset + layer * n_quads + quad + 1``. Each ply therefore owns one
contiguous element range, written with ``generate``.

Pure Python: no Abaqus, no numpy.
"""

from __future__ import division, print_function

from cai import ConfigurationError
from cai.mesh2d import Shape, mesh_ply

#: C3D8 face id of each quad edge (n1-n2, n2-n3, n3-n4, n4-n1) and of the
#: bottom (z-) and top (z+) faces, for the node ordering used here.
EDGE_FACE = ("S3", "S4", "S5", "S6")
BOTTOM_FACE = "S1"
TOP_FACE = "S2"

#: Element-set names of the three damage categories.
CATEGORY_SETS = (
    ("pristine", "CDM_PRISTINE"),
    ("matrix", "CDM_MATRIX_CRACK"),
    ("fiber", "CDM_FIBER_DAMAGE"),
)


class BlockMesh(object):
    """Extruded hexahedral mesh of one ply block.

    Parameters
    ----------
    block : cai.layup.PlyBlock
        Ply block.
    mesh2d : cai.mesh2d.QuadMesh
        In-plane mesh of the block.
    layers_per_ply : int
        Elements through each ply.
    node_offset, elem_offset : int
        Label offsets (labels start at offset + 1).

    Attributes
    ----------
    z : list of float
        z-coordinate of every node level, bottom to top.
    n2d, nq, nlayers : int
        In-plane node and quad counts, and element layers.
    """

    def __init__(self, block, mesh2d, layers_per_ply, node_offset, elem_offset):
        self.block = block
        self.mesh2d = mesh2d
        self.layers_per_ply = int(layers_per_ply)
        self.node_offset = node_offset
        self.elem_offset = elem_offset
        self.z = [block.z_bottom]
        for ply in block.plies:
            dz = ply.thickness / self.layers_per_ply
            for k in range(1, self.layers_per_ply + 1):
                self.z.append(ply.z_bottom + k * dz)
        self.z[-1] = block.z_top
        self.n2d = len(mesh2d.nodes)
        self.nq = len(mesh2d.quads)
        self.nlayers = len(self.z) - 1

    @property
    def n_nodes(self):
        return self.n2d * len(self.z)

    @property
    def n_elements(self):
        return self.nq * self.nlayers

    def node_label(self, i2d, level):
        """Global label of in-plane node ``i2d`` at node level ``level``."""
        return self.node_offset + level * self.n2d + i2d + 1

    def elem_label(self, q, layer):
        """Global label of quad ``q`` extruded in element layer ``layer``."""
        return self.elem_offset + layer * self.nq + q + 1

    def iter_nodes(self):
        """Yield ``(label, x, y, z)`` for every node of the block."""
        for level, z in enumerate(self.z):
            base = self.node_offset + level * self.n2d
            for i, (x, y) in enumerate(self.mesh2d.nodes):
                yield base + i + 1, x, y, z

    def iter_elements(self):
        """Yield ``(label, n1..n8)`` for every C3D8R element of the block."""
        n2d = self.n2d
        for layer in range(self.nlayers):
            lo = self.node_offset + layer * n2d + 1
            hi = lo + n2d
            base = self.elem_offset + layer * self.nq + 1
            for q, (a, b, c, d) in enumerate(self.mesh2d.quads):
                yield (base + q, lo + a, lo + b, lo + c, lo + d, hi + a, hi + b, hi + c, hi + d)

    def ply_layers(self, ply):
        """Element layer range of ``ply`` inside the block."""
        k = self.block.plies.index(ply)
        return range(k * self.layers_per_ply, (k + 1) * self.layers_per_ply)

    def ply_range(self, ply):
        """``(first, last)`` element labels of ``ply`` (contiguous)."""
        layers = self.ply_layers(ply)
        return self.elem_label(0, layers[0]), self.elem_label(self.nq - 1, layers[-1])


class CAIMesh(object):
    """Complete 3-D mesh of the CAI plate with all sets and surfaces.

    Parameters
    ----------
    cfg : dict
        Loaded configuration.
    laminate : cai.layup.Laminate
        Laminate definition.
    damage_map : cai.damage_map.DamageMap
        Seeded impact damage.
    element_size : float
        Damage-zone element size from :class:`cai.turon.CohesiveSizing`.
    center : tuple of float
        Impact location.

    Attributes
    ----------
    blocks : list of BlockMesh
    warnings : list of str
    """

    def __init__(self, cfg, laminate, damage_map, element_size, center):
        self.cfg = cfg
        self.laminate = laminate
        self.dmap = damage_map
        self.h = float(element_size)
        self.center = center
        g, m, bc = cfg["GEOMETRY"], cfg["MESH"], cfg["BOUNDARY"]
        self.L, self.W = float(g["length"]), float(g["width"])
        self.tol = float(m.get("tolerance", 1e-6))
        self.es = float(bc["end_support_length"])
        self.inset = float(bc["side_support_inset"])
        self.warnings = []
        self.blocks = []
        node_off = elem_off = 0
        for block in laminate.blocks:
            angle = self._mesh_angle(block)
            q = mesh_ply(self.L, self.W, angle, self.zone_shapes(block), self.h,
                         float(m["coarse_size"]), float(m.get("growth_rate", 1.2)),
                         (self.es, self.L - self.es), (self.inset, self.W - self.inset),
                         float(m.get("zone_margin", 0.0)))
            if q.zone.clipped:
                self.warnings.append(
                    "%s: damage zone clipped at %s to stay inside the supports"
                    % (block.name, ", ".join(q.zone.clipped)))
            bm = BlockMesh(block, q, m.get("elements_per_ply", 1), node_off, elem_off)
            self.blocks.append(bm)
            node_off += bm.n_nodes
            elem_off += bm.n_elements
        self.n_nodes = node_off
        self.n_elements = elem_off
        self.rp_label = node_off + 1

    # ------------------------------------------------------------ meshing
    def _mesh_angle(self, block):
        if not self.cfg["MESH"].get("fiber_aligned", True):
            return 0.0
        angle = block.uniform_angle()
        if angle is None:
            raise ConfigurationError(
                "%s mixes ply angles %s; fibre-aligned meshing needs one angle per "
                "block (use LAYUP blocks='auto' or MESH fiber_aligned=False)"
                % (block.name, block.angles))
        return angle

    def zone_shapes(self, block):
        """Regions the fine damage zone of ``block`` must contain.

        Includes the seeded delaminations at the block's interfaces, the
        seeded matrix cracks of its plies and a minimum square around the
        impact site.

        Parameters
        ----------
        block : cai.layup.PlyBlock

        Returns
        -------
        list of cai.mesh2d.Shape
        """
        shapes = [Shape(self.center, float(self.cfg["MESH"].get("min_zone_half_size", 10.0)))]
        for itf in self.laminate.interfaces:
            if itf.lower is block or itf.upper is block:
                a, b, ang = self.dmap.delamination_axes(itf)
                if a > 0.0:
                    shapes.append(Shape(self.center, a, b, ang))
        for ply in block.plies:
            r = self.dmap.matrix_crack_radius(ply)
            if r > 0.0:
                shapes.append(Shape(self.center, r))
        return shapes

    def quality(self):
        """Per-block in-plane mesh quality.

        Returns
        -------
        list of tuple
            ``(block_name, angle, n_elements, n_zone_quads, quality_dict)``.
        """
        return [(bm.block.name, bm.mesh2d.angle, bm.n_elements, len(bm.mesh2d.zone_quads),
                 bm.mesh2d.quality()) for bm in self.blocks]

    # -------------------------------------------------------- element sets
    def ply_ranges(self):
        """``[(ply, first_label, last_label)]`` for every ply."""
        out = []
        for bm in self.blocks:
            for ply in bm.block.plies:
                first, last = bm.ply_range(ply)
                out.append((ply, first, last))
        return out

    def damage_categories(self):
        """Element labels of each damage category.

        Returns
        -------
        dict
            ``{"pristine": [...], "matrix": [...], "fiber": [...]}``.
        """
        cats = dict((c, []) for c, _ in CATEGORY_SETS)
        for bm in self.blocks:
            centroids = [bm.mesh2d.centroid(q) for q in range(bm.nq)]
            for ply in bm.block.plies:
                kinds = [self.dmap.classify_element(ply, x, y) for x, y in centroids]
                for layer in bm.ply_layers(ply):
                    base = bm.elem_label(0, layer)
                    for q, kind in enumerate(kinds):
                        cats[kind].append(base + q)
        return cats

    # ------------------------------------------------------------ surfaces
    def interface_surfaces(self):
        """Element faces on both sides of every interface, split by seeded damage.

        Returns
        -------
        list of dict
            One entry per interface with keys ``interface`` and
            ``lower``/``upper``, each ``{"bonded": (face, labels),
            "delam": (face, labels)}``.
        """
        out = []
        for itf in self.laminate.interfaces:
            lower = self.blocks[itf.lower.index - 1]
            upper = self.blocks[itf.upper.index - 1]
            entry = {"interface": itf}
            for side, bm, layer, face in (("lower", lower, lower.nlayers - 1, TOP_FACE),
                                          ("upper", upper, 0, BOTTOM_FACE)):
                groups = {"bonded": [], "delam": []}
                for q in range(bm.nq):
                    x, y = bm.mesh2d.centroid(q)
                    key = "delam" if self.dmap.is_delaminated(itf, x, y) else "bonded"
                    groups[key].append(bm.elem_label(q, layer))
                entry[side] = dict((k, (face, v)) for k, v in groups.items())
            out.append(entry)
        return out

    def block_faces(self, bm, top):
        """``(face, labels)`` of the top or bottom face of a block."""
        layer = bm.nlayers - 1 if top else 0
        return (TOP_FACE if top else BOTTOM_FACE,
                [bm.elem_label(q, layer) for q in range(bm.nq)])

    def end_faces(self, x_end):
        """Element faces on the plate end ``x = x_end``.

        Returns
        -------
        dict
            ``{face_id: [labels]}``.
        """
        faces = {}
        for bm in self.blocks:
            nodes = bm.mesh2d.nodes
            for q, quad in enumerate(bm.mesh2d.quads):
                for e in range(4):
                    a, b = quad[e], quad[(e + 1) % 4]
                    if abs(nodes[a][0] - x_end) < self.tol and abs(nodes[b][0] - x_end) < self.tol:
                        labels = faces.setdefault(EDGE_FACE[e], [])
                        labels.extend(bm.elem_label(q, layer) for layer in range(bm.nlayers))
        return faces

    # ------------------------------------------------------------ node sets
    def _select(self, bm, levels, predicate):
        idx = [i for i, (x, y) in enumerate(bm.mesh2d.nodes) if predicate(x, y)]
        return [bm.node_label(i, lv) for lv in levels for i in idx]

    def node_sets(self):
        """Boundary-condition node sets.

        Returns
        -------
        dict
            ``END_X0``, ``END_SUPPORTS``, ``KNIFE_EDGES`` and ``ANCHOR_Y``
            label lists.

        Raises
        ------
        ConfigurationError
            If a support line has no nodes (should not happen: the mesher puts
            node lines on every support).
        """
        L, W, tol, es, inset = self.L, self.W, self.tol, self.es, self.inset
        bottom, top = self.blocks[0], self.blocks[-1]
        faces = ((bottom, [0]), (top, [len(top.z) - 1]))
        sets = {"END_X0": [], "END_SUPPORTS": [], "KNIFE_EDGES": [], "ANCHOR_Y": []}
        for bm in self.blocks:
            all_levels = range(len(bm.z))
            sets["END_X0"] += self._select(bm, all_levels, lambda x, y: abs(x) < tol)
            # Rigid-body Y restraint: x = 0 node line closest to the impact y
            ys = [y for x, y in bm.mesh2d.nodes if abs(x) < tol]
            y_anchor = min(ys, key=lambda y: abs(y - self.center[1]))
            sets["ANCHOR_Y"] += self._select(
                bm, all_levels, lambda x, y: abs(x) < tol and abs(y - y_anchor) < tol)
        for bm, levels in faces:
            sets["END_SUPPORTS"] += self._select(
                bm, levels, lambda x, y: x <= es + tol or x >= L - es - tol)
            sets["KNIFE_EDGES"] += self._select(
                bm, levels, lambda x, y: es - tol <= x <= L - es + tol and
                (abs(y - inset) < tol or abs(y - (W - inset)) < tol))
        for name, labels in sets.items():
            if not labels:
                raise ConfigurationError("no nodes found for node set %s" % name)
        return sets

    def summary(self):
        """Readable mesh summary.

        Returns
        -------
        str
        """
        lines = ["Mesh: %d nodes, %d C3D8R elements, damage-zone size %.3f mm, %s" % (
            self.n_nodes, self.n_elements, self.h,
            "fibre-aligned" if self.cfg["MESH"].get("fiber_aligned", True) else "orthogonal")]
        for name, angle, ne, nz, qd in self.quality():
            lines.append("  %s %+6.1f deg: %7d elements (%6d in fibre-aligned zone), "
                         "min angle %.0f deg, max aspect %.1f" % (
                             name, angle, ne, nz, qd["min_angle"], qd["max_aspect"]))
            if qd["min_jacobian"] <= 0.0:
                raise ConfigurationError("%s has an inverted element" % name)
            if abs(qd["area"] - self.L * self.W) > 1e-6 * self.L * self.W:
                raise ConfigurationError("%s does not cover the plate" % name)
        return "\n".join(lines)
