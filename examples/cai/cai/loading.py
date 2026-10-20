"""
Explicit compression step, ASTM D7137 fixture boundary conditions and outputs (Abaqus).

Boundary conditions (global axes, X = loading direction):

=====================  ========================================================
Region                 Constraint
=====================  ========================================================
Loaded end (x = L)     Kinematic coupling (U1) to reference point ``RP_LOAD``;
                       the RP moves by ``-applied_displacement`` along X with a
                       smooth-step amplitude (displacement control).
Reacting end (x = 0)   U1 = 0 (rigid base plate).
End supports           U3 = 0 on both faces over ``end_support_length`` at each
                       end (anti-rotation clamping plates of the fixture).
Side knife edges       U3 = 0 along y = inset and y = W - inset on both faces,
                       between the end supports (simple support).
Rigid-body Y           U2 = 0 on the node line x = 0, y = y_c.
=====================  ========================================================

Must run inside the Abaqus/CAE kernel.
"""

from __future__ import division, print_function

from abaqusConstants import (ON, OFF, STEP, SEMI_AUTOMATIC, MODEL, AT_BEGINNING,
                             BELOW_MIN, KINEMATIC, WHOLE_SURFACE)

from cai import ModelBuildError
from cai.geometry import box

AMPLITUDE = "AMP_SMOOTH_RAMP"


def _union(arrays):
    """Concatenate Abaqus entity arrays (possibly from several instances)."""
    result = None
    for arr in arrays:
        if arr is None or len(arr) == 0:
            continue
        result = arr if result is None else result + arr
    return result


class LoadingBuilder(object):
    """Create the explicit step, fixture boundary conditions and output requests.

    Parameters
    ----------
    model : Model
        Abaqus model object.
    cfg : dict
        Loaded configuration.
    laminate : cai.layup.Laminate
        Laminate definition.
    center : tuple of float
        Impact location ``(x, y)``; its y-coordinate sets the rigid-body
        restraint line.
    """

    def __init__(self, model, cfg, laminate, center):
        self.model = model
        self.cfg = cfg
        self.laminate = laminate
        self.center = center
        self.L = float(cfg["GEOMETRY"]["length"])
        self.W = float(cfg["GEOMETRY"]["width"])
        self.t = laminate.thickness
        self.tol = float(cfg["MESH"].get("tolerance", 1e-4))
        self.step = cfg["STEP"]["name"]

    def build(self, assembly):
        """Create step, amplitude, sets, constraints, BCs and outputs.

        Parameters
        ----------
        assembly : Assembly
            Root assembly with the stacked instances.
        """
        self._step()
        self._sets(assembly)
        self._load_introduction(assembly)
        self._boundary_conditions(assembly)
        self._outputs(assembly)

    # ------------------------------------------------------------------- step
    def _step(self):
        s = self.cfg["STEP"]
        kwargs = dict(name=self.step, previous="Initial", timePeriod=s["time_period"],
                      nlgeom=ON,  # CompDam requires geometric nonlinearity
                      linearBulkViscosity=s["linear_bulk_viscosity"],
                      quadBulkViscosity=s["quadratic_bulk_viscosity"],
                      description="CAI compression, displacement control")
        dt = s.get("mass_scaling_target_dt")
        if dt:
            kwargs["massScaling"] = ((SEMI_AUTOMATIC, MODEL, AT_BEGINNING, 0.0, dt,
                                      BELOW_MIN, 0, 0, 0.0, 0.0, 0, None),)
        self.model.ExplicitDynamicsStep(**kwargs)
        self.model.SmoothStepAmplitude(name=AMPLITUDE, timeSpan=STEP,
                                       data=((0.0, 0.0), (s["time_period"], 1.0)))

    # ------------------------------------------------------------------- sets
    def _instances(self, a):
        return [a.instances[b.instance_name] for b in self.laminate.blocks]

    def _nodes(self, a, instances, z, x=None, y=None):
        L, W, tol = self.L, self.W, self.tol
        return _union([inst.nodes.getByBoundingBox(**box(L, W, z, z, tol, x=x, y=y))
                       for inst in instances])

    def _require(self, arr, what):
        if arr is None or len(arr) == 0:
            raise ModelBuildError("no mesh entities found for %s" % what)
        return arr

    def _sets(self, a):
        L, W, t, tol = self.L, self.W, self.t, self.tol
        b = self.cfg["BOUNDARY"]
        es, inset = b["end_support_length"], b["side_support_inset"]
        insts = self._instances(a)
        bottom, top = [insts[0]], [insts[-1]]

        x_min = _union([i.faces.getByBoundingBox(**box(L, W, 0.0, t, tol, x=(0.0, 0.0)))
                        for i in insts])
        x_max = _union([i.faces.getByBoundingBox(**box(L, W, 0.0, t, tol, x=(L, L)))
                        for i in insts])
        a.Set(faces=self._require(x_min, "x = 0 end face"), name="END_X0")
        a.Surface(side1Faces=self._require(x_max, "x = L end face"), name="END_XL")

        clamp = _union([
            self._nodes(a, bottom, 0.0, x=(0.0, es)), self._nodes(a, bottom, 0.0, x=(L - es, L)),
            self._nodes(a, top, t, x=(0.0, es)), self._nodes(a, top, t, x=(L - es, L)),
        ])
        a.Set(nodes=self._require(clamp, "end-support strips"), name="END_SUPPORTS")

        knife = []
        for y in (inset, W - inset):
            for group, z in ((bottom, 0.0), (top, t)):
                knife.append(self._nodes(a, group, z, x=(es, L - es), y=(y, y)))
        a.Set(nodes=self._require(_union(knife), "side knife edges"), name="KNIFE_EDGES")

        yc = self.center[1]
        anchor = _union([i.nodes.getByBoundingBox(**box(L, W, 0.0, t, tol,
                                                        x=(0.0, 0.0), y=(yc, yc)))
                         for i in insts])
        a.Set(nodes=self._require(anchor, "rigid-body anchor line"), name="ANCHOR_Y")

    # ----------------------------------------------------------------- loads
    def _load_introduction(self, a):
        rp = a.ReferencePoint(point=(self.L + 5.0, 0.5 * self.W, 0.5 * self.t))
        a.Set(referencePoints=(a.referencePoints[rp.id],), name="RP_LOAD")
        self.model.Coupling(name="CPL_LOADED_END", controlPoint=a.sets["RP_LOAD"],
                            surface=a.surfaces["END_XL"], influenceRadius=WHOLE_SURFACE,
                            couplingType=KINEMATIC, localCsys=None,
                            u1=ON, u2=OFF, u3=OFF, ur1=OFF, ur2=OFF, ur3=OFF)

    def _boundary_conditions(self, a):
        m = self.model
        d = float(self.cfg["BOUNDARY"]["applied_displacement"])
        m.DisplacementBC(name="BC_END_X0", createStepName="Initial",
                         region=a.sets["END_X0"], u1=0.0)
        m.DisplacementBC(name="BC_END_SUPPORTS", createStepName="Initial",
                         region=a.sets["END_SUPPORTS"], u3=0.0)
        m.DisplacementBC(name="BC_KNIFE_EDGES", createStepName="Initial",
                         region=a.sets["KNIFE_EDGES"], u3=0.0)
        m.DisplacementBC(name="BC_ANCHOR_Y", createStepName="Initial",
                         region=a.sets["ANCHOR_Y"], u2=0.0)
        # RP: all DOFs prescribed (no free massless DOFs in Explicit)
        m.DisplacementBC(name="BC_RP", createStepName="Initial", region=a.sets["RP_LOAD"],
                         u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0)
        m.boundaryConditions["BC_RP"].setValuesInStep(stepName=self.step, u1=-d,
                                                      amplitude=AMPLITUDE)

    # --------------------------------------------------------------- outputs
    def _outputs(self, a):
        s = self.cfg["STEP"]
        m = self.model
        field = ["S", "LE", "U", "RF", "SDV"]
        if self.cfg["COMPDAM"].get("element_deletion"):
            field.append("STATUS")
        if self.cfg["INTERFACE"]["type"] == "cohesive" and self.laminate.interfaces:
            field += ["CSTRESS", "CSDMG"]
        try:
            m.fieldOutputRequests["F-Output-1"].setValues(
                variables=tuple(field), numIntervals=s["field_output_intervals"])
        except Exception as exc:  # unsupported variable in this Abaqus version
            print("  WARNING: field output %s rejected (%s); using basic set" % (field, exc))
            m.fieldOutputRequests["F-Output-1"].setValues(
                variables=("S", "LE", "U", "RF", "SDV"),
                numIntervals=s["field_output_intervals"])

        m.historyOutputRequests["H-Output-1"].setValues(
            variables=("ALLAE", "ALLDMD", "ALLFD", "ALLIE", "ALLKE", "ALLPD", "ALLSE",
                       "ALLVD", "ALLWK", "ETOTAL"),
            numIntervals=s["history_output_intervals"])
        m.HistoryOutputRequest(name="H-RP_LOAD", createStepName=self.step,
                               variables=("U1", "RF1"), region=a.sets["RP_LOAD"],
                               numIntervals=s["history_output_intervals"])
