"""
CompDam material, solid section and ply orientations (Abaqus).

* The material is an Abaqus/Explicit user material (VUMAT) with the
  40-constant CompDam property vector, a density and ``*Depvar``. The state
  variable names and ``*Characteristic Length, definition=USER`` have no
  CAE support, so :mod:`cai.inp_patcher` adds them to the input deck.
* Every ply gets a rectangular local coordinate system rotated about the
  stacking direction (global Z) by the ply angle: material direction 1 is
  the fibre, 2 the in-plane transverse direction, 3 the ply normal.

Must run inside the Abaqus/CAE kernel.
"""

from __future__ import division, print_function

import math

from abaqusConstants import (MECHANICAL, CARTESIAN, SYSTEM, AXIS_3, ROTATION_NONE,
                             STACK_3, MIDDLE_SURFACE, FROM_SECTION)
import regionToolset

from cai import compdam
from cai import ModelBuildError


class MaterialLayupBuilder(object):
    """Define the CompDam material and assign sections and ply orientations.

    Parameters
    ----------
    model : Model
        Abaqus model object.
    cfg : dict
        Loaded configuration.
    laminate : cai.layup.Laminate
        Laminate definition.
    material_name : str
        Name of the material preset (also used as the Abaqus material name).
    material : dict
        Material preset from ``config.MATERIALS``.

    Attributes
    ----------
    nstatev : int
        Number of state variables allocated for CompDam.
    props : list of float
        The 40 user-material constants.
    section_name : str
        Name of the homogeneous solid section.
    """

    def __init__(self, model, cfg, laminate, material_name, material):
        self.model = model
        self.cfg = cfg
        self.laminate = laminate
        self.material_name = material_name
        self.material = material
        flags = cfg["COMPDAM"]["feature_flags"]
        self.nstatev = compdam.number_of_state_variables(flags)
        self.props = compdam.build_property_vector(material, flags)
        self.section_name = "SEC_" + material_name

    def build(self, parts):
        """Create the material and section, then assign them to all parts.

        Parameters
        ----------
        parts : dict
            ``{block_name: Part}``.

        Raises
        ------
        ModelBuildError
            If a ply set is missing on a part.
        """
        self._create_material()
        self.model.HomogeneousSolidSection(name=self.section_name,
                                           material=self.material_name, thickness=None)
        for block in self.laminate.blocks:
            p = parts[block.name]
            p.SectionAssignment(region=p.sets["ALL"], sectionName=self.section_name,
                                offset=0.0, offsetType=MIDDLE_SURFACE, offsetField="",
                                thicknessAssignment=FROM_SECTION)
            for ply in block.plies:
                self._orient_ply(p, ply)

    def _create_material(self):
        mat = self.model.Material(name=self.material_name,
                                  description=self.material.get("description", ""))
        mat.Density(table=((self.material["density"],),))
        delete = self.cfg["COMPDAM"].get("element_deletion", False)
        if delete:
            mat.Depvar(n=self.nstatev, deleteVar=compdam.SDV_STATUS)
        else:
            mat.Depvar(n=self.nstatev)
        mat.UserMaterial(type=MECHANICAL, mechanicalConstants=tuple(self.props))
        print("  material %s: %d constants, %d state variables, flags %s" % (
            self.material_name, len(self.props), self.nstatev,
            self.cfg["COMPDAM"]["feature_flags"]))

    def _orient_ply(self, p, ply):
        if ply.name not in p.sets.keys():
            raise ModelBuildError("set %s missing on part %s" % (ply.name, p.name))
        t = math.radians(ply.angle)
        csys = p.DatumCsysByThreePoints(
            name="CSYS_" + ply.name, coordSysType=CARTESIAN, origin=(0.0, 0.0, 0.0),
            point1=(math.cos(t), math.sin(t), 0.0),  # local 1 = fibre direction
            point2=(-math.sin(t), math.cos(t), 0.0))  # in the local 1-2 plane
        region = regionToolset.Region(cells=p.sets[ply.name].cells)
        p.MaterialOrientation(region=region, orientationType=SYSTEM, axis=AXIS_3,
                              localCsys=p.datums[csys.id], fieldName="",
                              additionalRotationType=ROTATION_NONE, angle=0.0,
                              additionalRotationField="", stackDirection=STACK_3)
