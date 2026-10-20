"""
Geometry, partitioning, ply-aligned solid meshing and assembly (Abaqus).

Each ply block of the laminate becomes its own 3-D deformable part:

* the plate is extruded to the block thickness,
* cells are cut at every ply boundary inside the block, so element faces
  line up exactly with the ply interfaces (``elements_per_ply`` elements
  through each ply),
* in-plane cuts from :class:`cai.mesh_plan.MeshPlan` give a structured
  hexahedral mesh, refined around the impact site and coarse near the
  supports,
* elements are C3D8R (Explicit library) with enhanced hourglass control and
  distortion control, as CompDam requires.

Blocks are then instanced and stacked in z. The blocks are joined by
cohesive contact or tie constraints (:mod:`cai.interactions`).

Must run inside the Abaqus/CAE kernel.
"""

from __future__ import division, print_function

from abaqusConstants import (THREE_D, DEFORMABLE_BODY, YZPLANE, XZPLANE, XYPLANE,
                             SINGLE, FINER, FIXED, HEX, STRUCTURED, ON, OFF,
                             EXPLICIT, C3D8R, C3D6, C3D4, ENHANCED, RELAX_STIFFNESS,
                             STIFFNESS, DEFAULT, AVERAGE_STRAIN, ORTHOGONAL, CENTROID)
import mesh

from cai import ModelBuildError

_HOURGLASS = {
    "ENHANCED": ENHANCED,
    "RELAX_STIFFNESS": RELAX_STIFFNESS,
    "STIFFNESS": STIFFNESS,
    "DEFAULT": DEFAULT,
}

_KINEMATIC_SPLIT = {
    "AVERAGE_STRAIN": AVERAGE_STRAIN,
    "ORTHOGONAL": ORTHOGONAL,
    "CENTROID": CENTROID,
}


def box(L, W, z0, z1, tol, x=None, y=None):
    """Keyword arguments for ``getByBoundingBox`` around a region of the plate.

    Parameters
    ----------
    L, W : float
        Plate length and width.
    z0, z1 : float
        Through-thickness range.
    tol : float
        Tolerance added on every side.
    x, y : tuple of float, optional
        Restrict the in-plane range; defaults to the whole plate.

    Returns
    -------
    dict
        ``xMin``, ``yMin``, ``zMin``, ``xMax``, ``yMax``, ``zMax``.
    """
    x = x or (0.0, L)
    y = y or (0.0, W)
    return dict(xMin=x[0] - tol, yMin=y[0] - tol, zMin=z0 - tol,
                xMax=x[1] + tol, yMax=y[1] + tol, zMax=z1 + tol)


class PlyBlockPartBuilder(object):
    """Create, partition and mesh one Abaqus part per ply block.

    Parameters
    ----------
    model : Model
        Abaqus model object.
    cfg : dict
        Loaded configuration.
    laminate : cai.layup.Laminate
        Laminate definition.
    plan : cai.mesh_plan.MeshPlan
        Partition and seeding plan.

    Notes
    -----
    Part-level sets created on each part (local coordinates, block bottom at
    z = 0):

    * ``ALL``: all cells,
    * ``PLY_nn``: cells of ply ``nn`` (global ply number),

    and surfaces ``TOP`` (z = block thickness) and ``BOTTOM`` (z = 0).
    """

    def __init__(self, model, cfg, laminate, plan):
        self.model = model
        self.cfg = cfg
        self.laminate = laminate
        self.plan = plan
        self.L = plan.L
        self.W = plan.W
        self.tol = plan.tol

    def build(self):
        """Create and mesh all ply-block parts.

        Returns
        -------
        dict
            ``{block_name: Part}``.

        Raises
        ------
        ModelBuildError
            If partitioning or meshing fails for a block.
        """
        parts = {}
        for block in self.laminate.blocks:
            try:
                p = self._create_part(block)
                self._partition(p, block)
                self._create_sets(p, block)
                self._seed(p)
                self._mesh(p, block)
            except ModelBuildError:
                raise
            except Exception as exc:  # Abaqus raises plain exceptions
                raise ModelBuildError("failed to build part %s: %s" % (block.name, exc))
            parts[block.name] = p
            print("  %s: %d elements, %d nodes" % (block.name, len(p.elements), len(p.nodes)))
        return parts

    # ------------------------------------------------------------------ steps
    def _create_part(self, block):
        sketch = self.model.ConstrainedSketch(name="__profile__",
                                              sheetSize=2.0 * max(self.L, self.W))
        sketch.rectangle(point1=(0.0, 0.0), point2=(self.L, self.W))
        p = self.model.Part(name=block.name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
        p.BaseSolidExtrude(sketch=sketch, depth=block.thickness)
        del self.model.sketches["__profile__"]
        return p

    def _cut(self, p, principal_plane, offset):
        datum = p.DatumPlaneByPrincipalPlane(principalPlane=principal_plane, offset=offset)
        p.PartitionCellByDatumPlane(datumPlane=p.datums[datum.id], cells=p.cells)

    def _partition(self, p, block):
        for x in self.plan.x_cuts:
            self._cut(p, YZPLANE, x)
        for y in self.plan.y_cuts:
            self._cut(p, XZPLANE, y)
        # Ply boundaries inside the block (local z)
        for ply in block.plies[:-1]:
            self._cut(p, XYPLANE, ply.z_top - block.z_bottom)

    def _create_sets(self, p, block):
        L, W, tol, h = self.L, self.W, self.tol, block.thickness
        p.Set(cells=p.cells, name="ALL")
        for ply in block.plies:
            z0 = ply.z_bottom - block.z_bottom
            cells = p.cells.getByBoundingBox(**box(L, W, z0, z0 + ply.thickness, tol))
            if len(cells) == 0:
                raise ModelBuildError("no cells found for %s in %s" % (ply.name, block.name))
            p.Set(cells=cells, name=ply.name)
        p.Surface(side1Faces=p.faces.getByBoundingBox(**box(L, W, h, h, tol)), name="TOP")
        p.Surface(side1Faces=p.faces.getByBoundingBox(**box(L, W, 0.0, 0.0, tol)), name="BOTTOM")

    def _seed(self, p):
        for edge in p.edges:
            v0, v1 = edge.getVertices()
            p0 = p.vertices[v0].pointOn[0]
            p1 = p.vertices[v1].pointOn[0]
            rule = self.plan.edge_seed(p0, p1)
            edges = p.edges[edge.index:edge.index + 1]
            if rule[0] == "number":
                p.seedEdgeByNumber(edges=edges, number=rule[1], constraint=FIXED)
            elif rule[0] == "size":
                p.seedEdgeBySize(edges=edges, size=rule[1], deviationFactor=0.1,
                                 constraint=FINER)
            else:
                _, h_min, h_max, c_fine = rule
                axis = 0 if abs(p1[0] - p0[0]) > abs(p1[1] - p0[1]) else 1
                # end1Edges puts the smallest elements at curve parameter 0,
                # i.e. at the edge's first vertex.
                fine_at_start = abs(p0[axis] - c_fine) < abs(p1[axis] - c_fine)
                key = "end1Edges" if fine_at_start else "end2Edges"
                p.seedEdgeByBias(biasMethod=SINGLE, minSize=h_min, maxSize=h_max,
                                 constraint=FINER, **{key: edges})

    def _element_types(self):
        m = self.cfg["MESH"]
        hg = _HOURGLASS[m.get("hourglass_control", "ENHANCED").upper()]
        ks = _KINEMATIC_SPLIT[m.get("kinematic_split", "AVERAGE_STRAIN").upper()]
        dc = ON if m.get("distortion_control", True) else OFF
        hex_type = mesh.ElemType(elemCode=C3D8R, elemLibrary=EXPLICIT,
                                 kinematicSplit=ks, secondOrderAccuracy=OFF,
                                 hourglassControl=hg, distortionControl=dc)
        wedge_type = mesh.ElemType(elemCode=C3D6, elemLibrary=EXPLICIT, distortionControl=dc)
        tet_type = mesh.ElemType(elemCode=C3D4, elemLibrary=EXPLICIT, distortionControl=dc)
        return (hex_type, wedge_type, tet_type)

    def _mesh(self, p, block):
        p.setMeshControls(regions=p.cells, elemShape=HEX, technique=STRUCTURED)
        p.setElementType(regions=(p.cells,), elemTypes=self._element_types())
        p.generateMesh()
        if len(p.elements) == 0:
            raise ModelBuildError("meshing produced no elements for %s" % block.name)
        for ply in block.plies:
            if len(p.sets[ply.name].elements) == 0:
                raise ModelBuildError("%s of %s has no elements" % (ply.name, block.name))


class AssemblyBuilder(object):
    """Instance the ply-block parts and stack them through the thickness.

    Parameters
    ----------
    model : Model
        Abaqus model object.
    laminate : cai.layup.Laminate
        Laminate definition.
    """

    def __init__(self, model, laminate):
        self.model = model
        self.laminate = laminate

    def build(self, parts):
        """Create dependent instances translated to their global z-position.

        Parameters
        ----------
        parts : dict
            ``{block_name: Part}`` from :class:`PlyBlockPartBuilder`.

        Returns
        -------
        Assembly
            The root assembly, with an ``ALL_PLIES`` element set.
        """
        a = self.model.rootAssembly
        for block in self.laminate.blocks:
            a.Instance(name=block.instance_name, part=parts[block.name], dependent=ON)
            if block.z_bottom > 0.0:
                a.translate(instanceList=(block.instance_name,),
                            vector=(0.0, 0.0, block.z_bottom))
        elements = None
        for block in self.laminate.blocks:
            e = a.instances[block.instance_name].elements
            elements = e if elements is None else elements + e
        a.Set(elements=elements, name="ALL_PLIES")
        return a
