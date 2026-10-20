"""
Ply-block interfaces: cohesive contact with seeded delaminations, or ties (Abaqus).

``INTERFACE["type"] == "cohesive"``
    Each interface gets two Abaqus/Explicit surface-to-surface contact pairs
    (penalty enforcement, small sliding):

    * the *bonded* region uses a surface-based cohesive behaviour with
      quadratic-traction initiation and Benzeggagh-Kenane energy-based
      evolution (delamination growth), plus hard normal contact and Coulomb
      friction once debonded,
    * the *delaminated* region (pre-defined impact damage) uses hard,
      frictional contact only, so the sublaminates can separate and
      buckle locally during compression.

``INTERFACE["type"] == "tie"``
    Perfectly bonded blocks; use it for debugging and pristine baselines.

Must run inside the Abaqus/CAE kernel.
"""

from __future__ import division, print_function

from abaqusConstants import (ON, OFF, PENALTY, FRACTION, HARD, DEFAULT, QUAD_TRACTION,
                             ENERGY, BK, SMALL, COMPUTED, INITIAL_NODES, UNCOUPLED,
                             ISOTROPIC)

from cai import ModelBuildError

COHESIVE_PROPERTY = "IP_COHESIVE"
FRICTION_PROPERTY = "IP_FRICTION"


class InterfaceBuilder(object):
    """Create interaction properties and interface contact pairs or ties.

    Parameters
    ----------
    model : Model
        Abaqus model object.
    cfg : dict
        Loaded configuration.
    laminate : cai.layup.Laminate
        Laminate definition.
    material : dict
        Material preset (its ``interface`` group is used).
    seeder : cai.damage.DamageSeeder
        Holds the bonded/delaminated surface names per interface.
    """

    def __init__(self, model, cfg, laminate, material, seeder):
        self.model = model
        self.cfg = cfg
        self.laminate = laminate
        self.props = material.get("interface")
        self.seeder = seeder
        self.kind = cfg["INTERFACE"]["type"]

    def build(self, assembly):
        """Create all interface definitions.

        Parameters
        ----------
        assembly : Assembly
            Root assembly with the stacked instances.

        Raises
        ------
        ModelBuildError
            If cohesive interfaces are requested but the material preset has
            no ``interface`` properties.
        """
        if not self.laminate.interfaces:
            return
        if self.kind == "tie":
            self._ties(assembly)
            return
        if not self.props:
            raise ModelBuildError("material preset lacks 'interface' properties")
        self._create_properties()
        n_pairs = 0
        for itf in self.laminate.interfaces:
            n_pairs += self._contact_pairs(assembly, itf)
        print("  %d contact pairs on %d interfaces" % (n_pairs, len(self.laminate.interfaces)))

    # ------------------------------------------------------------- properties
    def _add_friction(self, prop):
        prop.TangentialBehavior(formulation=PENALTY, directionality=ISOTROPIC,
                                slipRateDependency=OFF, pressureDependency=OFF,
                                temperatureDependency=OFF, dependencies=0,
                                table=((self.props["friction"],),),
                                shearStressLimit=None, maximumElasticSlip=FRACTION,
                                fraction=0.005, elasticSlipStiffness=None)
        prop.NormalBehavior(pressureOverclosure=HARD, allowSeparation=ON,
                            constraintEnforcementMethod=DEFAULT)

    def _create_properties(self):
        ip = self.props
        coh = self.model.ContactProperty(COHESIVE_PROPERTY)
        self._add_friction(coh)
        coh.CohesiveBehavior(eligibility=INITIAL_NODES, defaultPenalties=OFF,
                             coupling=UNCOUPLED,
                             table=((ip["penalty_normal"], ip["penalty_shear"],
                                     ip["penalty_shear"]),))
        coh.Damage(criterion=QUAD_TRACTION,
                   initTable=((ip["strength_normal"], ip["strength_shear"],
                               ip["strength_shear"]),),
                   useEvolution=ON, evolutionType=ENERGY, useMixedMode=ON,
                   mixedModeType=BK, exponent=ip["eta_BK"],
                   evolTable=((ip["GIc"], ip["GIIc"], ip["GIIc"]),))
        fric = self.model.ContactProperty(FRICTION_PROPERTY)
        self._add_friction(fric)

    # ------------------------------------------------------------------ pairs
    def _contact_pairs(self, a, itf):
        names = self.seeder.interface_surfaces.get(itf.name)
        if names is None:
            raise ModelBuildError("interface surfaces for %s were not created" % itf.name)
        low = a.instances[itf.lower.instance_name]
        up = a.instances[itf.upper.instance_name]
        n = 0
        for key, prop in (("bonded", COHESIVE_PROPERTY), ("delam", FRICTION_PROPERTY)):
            s_low, s_up = names["lower"][key], names["upper"][key]
            if s_low is None and s_up is None:
                continue
            if s_low is None or s_up is None:
                raise ModelBuildError(
                    "%s %s region exists on one side only; the meshes of the blocks "
                    "do not match" % (itf.name, key))
            # The upper (impact side) block is main; both meshes match.
            self.model.SurfaceToSurfaceContactExp(
                name="CP_%s_%s" % (itf.name, key.upper()), createStepName="Initial",
                main=up.surfaces[s_up], secondary=low.surfaces[s_low],
                mechanicalConstraint=PENALTY, sliding=SMALL,
                interactionProperty=prop)
            n += 1
        return n

    def _ties(self, a):
        for itf in self.laminate.interfaces:
            low = a.instances[itf.lower.instance_name]
            up = a.instances[itf.upper.instance_name]
            self.model.Tie(name="TIE_%s" % itf.name, main=up.surfaces["BOTTOM"],
                           secondary=low.surfaces["TOP"],
                           positionToleranceMethod=COMPUTED, adjust=ON,
                           tieRotations=ON, thickness=ON)
        print("  %d tie constraints" % len(self.laminate.interfaces))
