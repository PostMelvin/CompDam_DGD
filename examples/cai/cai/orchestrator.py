"""
Top-level orchestration of the CAI model build, input-deck patching and job submission.

:class:`CAIModel` runs the builders in order:

1. validate the configuration (:mod:`cai.settings`),
2. build the laminate, mesh plan and damage map (pure Python),
3. create and mesh the ply-block parts (:mod:`cai.geometry`),
4. define the CompDam material, sections and orientations
   (:mod:`cai.material_section`),
5. seed the impact damage (:mod:`cai.damage`),
6. assemble, then create the interfaces (:mod:`cai.interactions`),
7. create the explicit step, BCs and outputs (:mod:`cai.loading`),
8. write the input deck and add the CompDam keywords
   (:mod:`cai.inp_patcher`),
9. submit the job if requested.

Must run inside the Abaqus/CAE kernel.
"""

from __future__ import division, print_function

import os
import shutil

from abaqus import mdb
from abaqusConstants import ANALYSIS, DOUBLE_PLUS_PACK, FULL, DOMAIN, DEFAULT, OFF

from cai import ConfigurationError, ModelBuildError
from cai import compdam, settings, validation
from cai.layup import Laminate, membrane_properties
from cai.mesh_plan import MeshPlan
from cai.damage_map import DamageMap
from cai.geometry import PlyBlockPartBuilder, AssemblyBuilder
from cai.material_section import MaterialLayupBuilder
from cai.damage import DamageSeeder
from cai.interactions import InterfaceBuilder
from cai.loading import LoadingBuilder


class CAIModel(object):
    """Build, write and (optionally) submit one CAI analysis.

    Parameters
    ----------
    cfg : dict
        Configuration from :func:`cai.settings.load_config`.
    energy : float, optional
        Impact energy [J]; defaults to ``IMPACT["default_energy"]``.

    Attributes
    ----------
    name : str
        Model and job name.
    laminate : cai.layup.Laminate
    damage_map : cai.damage_map.DamageMap
    model : Model
        The Abaqus model (after :meth:`build`).
    inp_path : str
        Path of the patched input deck (after :meth:`write_input`).

    Examples
    --------
    Inside ``abaqus cae noGUI=cai_model.py``::

        cfg = settings.load_config()
        cai = CAIModel(cfg, energy=30.0)
        cai.build()
        cai.write_input()
        cai.submit()
    """

    def __init__(self, cfg, energy=None):
        self.cfg = cfg
        self.energy = float(cfg["IMPACT"]["default_energy"] if energy is None else energy)
        self.warnings = settings.validate(cfg)
        self.material_name, self.material = settings.active_material(cfg)
        self.name = settings.run_name(cfg, self.energy)
        self.center = settings.impact_center(cfg)
        self.laminate = Laminate.from_config(cfg["LAYUP"], self.material)
        self.plan = MeshPlan(cfg["GEOMETRY"], cfg["MESH"], cfg["BOUNDARY"], self.center)
        self.damage_map = DamageMap(self.laminate, cfg["DAMAGE"], cfg["IMPACT"],
                                    self.center, self.energy)
        self.nstatev = compdam.number_of_state_variables(cfg["COMPDAM"]["feature_flags"])
        self.model = None
        self.seeder = None
        self.ic_sets = []
        self.inp_path = None
        self.user_subroutine = settings.resolve_path(cfg, cfg["COMPDAM"]["user_subroutine"])

    # ----------------------------------------------------------------- report
    def report(self):
        """Print a pre-build summary: laminate, CLT stiffness, damage and checks."""
        print("=" * 78)
        print("CAI model %s" % self.name)
        print(self.laminate.summary())
        clt = membrane_properties(self.laminate, self.material["elastic"])
        print("CLT membrane: Ex = %.0f MPa, Ey = %.0f MPa, Gxy = %.0f MPa, nuxy = %.3f" % (
            clt["Ex"], clt["Ey"], clt["Gxy"], clt["nuxy"]))
        print(self.damage_map.summary())
        print(validation.material_report(self.material_name, self.material,
                                         self.cfg["NCAMP_BASELINE"].get(self.material_name)))
        for w in self.warnings:
            print("WARNING: %s" % w)
        print("=" * 78)

    # ------------------------------------------------------------------ build
    def build(self):
        """Create the complete Abaqus model.

        Returns
        -------
        Model
            The Abaqus model object.

        Raises
        ------
        ModelBuildError
            If any build stage fails.
        """
        self.report()
        if self.name in mdb.models.keys():
            del mdb.models[self.name]
        self.model = mdb.Model(name=self.name,
                               description=self.cfg["MODEL"].get("description", ""))

        print("[1/6] parts, partitions and mesh")
        parts = PlyBlockPartBuilder(self.model, self.cfg, self.laminate, self.plan).build()

        print("[2/6] CompDam material, section and ply orientations")
        MaterialLayupBuilder(self.model, self.cfg, self.laminate,
                             self.material_name, self.material).build(parts)

        print("[3/6] pre-defined impact damage")
        tol = float(self.cfg["MESH"].get("tolerance", 1e-4))
        self.seeder = DamageSeeder(self.laminate, self.damage_map, tol)
        self.seeder.classify_elements(parts)
        if self.cfg["INTERFACE"]["type"] == "cohesive":
            self.seeder.split_interface_surfaces(parts)

        print("[4/6] assembly")
        assembly = AssemblyBuilder(self.model, self.laminate).build(parts)
        self.ic_sets = self.seeder.create_assembly_sets(assembly)

        print("[5/6] interfaces (%s)" % self.cfg["INTERFACE"]["type"])
        InterfaceBuilder(self.model, self.cfg, self.laminate, self.material,
                         self.seeder).build(assembly)

        print("[6/6] step, boundary conditions and outputs")
        LoadingBuilder(self.model, self.cfg, self.laminate, self.center).build(assembly)
        if "Model-1" in mdb.models.keys() and len(mdb.models.keys()) > 1:
            del mdb.models["Model-1"]
        return self.model

    # -------------------------------------------------------------- job/deck
    def _job_kwargs(self):
        j = self.cfg["JOB"]
        if not os.path.isfile(self.user_subroutine):
            raise ConfigurationError("CompDam subroutine not found: %s" % self.user_subroutine)
        version = os.path.join(os.path.dirname(self.user_subroutine), "version.for")
        if not os.path.isfile(version):
            print("WARNING: %s is missing; run 'python setup.py' in the CompDam root or "
                  "rename for/version.for.nogit before submitting." % version)
        return dict(type=ANALYSIS, explicitPrecision=DOUBLE_PLUS_PACK,
                    nodalOutputPrecision=FULL, userSubroutine=self.user_subroutine,
                    numCpus=int(j["num_cpus"]), numDomains=int(j["num_domains"]),
                    parallelizationMethodExplicit=DOMAIN, multiprocessingMode=DEFAULT)

    def initial_conditions_block(self):
        """``*Initial Conditions, type=SOLUTION`` text for all damage categories.

        Returns
        -------
        str
        """
        d = self.cfg["DAMAGE"]
        phi0 = self.cfg["COMPDAM"].get("phi0", 0.0)
        states = {
            "pristine": compdam.initial_state(self.nstatev, phi0=phi0),
            "matrix": compdam.initial_state(
                self.nstatev, d2=d["matrix_cracking"]["d2"],
                alpha=d["matrix_cracking"]["alpha"], phi0=phi0),
            "fiber": compdam.initial_state(
                self.nstatev, d2=d["matrix_cracking"]["d2"],
                alpha=d["matrix_cracking"]["alpha"], d1T=d["fiber_damage"]["d1T"],
                d1C=d["fiber_damage"]["d1C"], phi0=phi0),
        }
        return compdam.format_initial_conditions(
            [(set_name, states[cat]) for set_name, cat in self.ic_sets])

    def write_input(self):
        """Write the input deck and add the CompDam-specific keywords.

        Returns
        -------
        str
            Absolute path of the patched ``.inp`` file.

        Raises
        ------
        ModelBuildError
            If the model was not built or the deck cannot be patched.
        """
        from cai.inp_patcher import InputDeckPatcher

        if self.model is None:
            raise ModelBuildError("call build() before write_input()")
        if self.name in mdb.jobs.keys():
            del mdb.jobs[self.name]
        job = mdb.Job(name=self.name, model=self.name, **self._job_kwargs())
        job.writeInput(consistencyChecking=OFF)
        path = os.path.abspath(self.name + ".inp")

        c = self.cfg["COMPDAM"]
        patcher = InputDeckPatcher.from_file(path)
        patcher.name_state_variables(self.material_name,
                                     compdam.format_depvar_names(self.nstatev))
        patcher.add_characteristic_length(self.material_name,
                                          c.get("characteristic_length_components", 3))
        patcher.insert_before_first_step(
            compdam.format_parameter_tables(c.get("parameters_float", {}),
                                            c.get("parameters_int", {})))
        patcher.insert_before_first_step(
            "**\n** COMPDAM INITIAL STATE (pristine + seeded impact damage)\n**\n"
            + self.initial_conditions_block())
        patcher.write(path)
        self.inp_path = path
        print("Patched input deck: %s" % path)
        print("Command-line equivalent:\n  abaqus job=%s input=%s user=%s double=both "
              "cpus=%d interactive" % (self.name, os.path.basename(path),
                                       self.user_subroutine, self.cfg["JOB"]["num_cpus"]))
        return path

    def submit(self, wait=True):
        """Submit the patched input deck as an Abaqus/Explicit job.

        Parameters
        ----------
        wait : bool, optional
            Block until the analysis finishes.

        Returns
        -------
        Job
            The Abaqus job object.
        """
        if self.inp_path is None:
            self.write_input()
        if self.name in mdb.jobs.keys():
            del mdb.jobs[self.name]
        job = mdb.JobFromInputFile(name=self.name, inputFileName=self.inp_path,
                                   **self._job_kwargs())
        job.submit(consistencyChecking=OFF)
        print("Submitted %s" % self.name)
        if wait:
            job.waitForCompletion()
            print("Job %s finished with status %s" % (self.name, job.status))
        return job


def prepare_work_dir(cfg):
    """Create and enter the working directory; copy the Abaqus environment file.

    Parameters
    ----------
    cfg : dict
        Loaded configuration.

    Returns
    -------
    str
        Absolute path of the working directory.
    """
    work = os.path.abspath(cfg["JOB"].get("work_dir") or os.getcwd())
    if not os.path.isdir(work):
        os.makedirs(work)
    env_src = cfg["COMPDAM"].get("env_file")
    if env_src:
        env_src = settings.resolve_path(cfg, env_src)
        env_dst = os.path.join(work, "abaqus_v6.env")
        if os.path.isfile(env_src) and not os.path.exists(env_dst):
            shutil.copy(env_src, env_dst)
            print("Copied %s -> %s" % (env_src, env_dst))
    os.chdir(work)
    return work


def run(cfg, energies):
    """Build (and optionally submit) one CAI model per impact energy.

    Parameters
    ----------
    cfg : dict
        Loaded configuration.
    energies : list of float
        Impact energies [J].

    Returns
    -------
    list of CAIModel
        The built models.
    """
    prepare_work_dir(cfg)
    built = []
    for energy in energies:
        cai = CAIModel(cfg, energy)
        cai.build()
        if cfg["JOB"].get("write_input", True) or cfg["JOB"].get("submit", False):
            cai.write_input()
        built.append(cai)
    if cfg["JOB"].get("save_cae", True):
        cae = os.path.abspath(cfg["MODEL"]["job_name"] + ".cae")
        mdb.saveAs(pathName=cae)
        print("Saved %s" % cae)
    if cfg["JOB"].get("submit", False):
        for cai in built:
            cai.submit(wait=True)
    return built
