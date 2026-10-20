"""
Orchestration: configuration -> laminate, damage, cohesive sizing, mesh -> job folder.

:class:`CAIJob` runs the stages for one impact energy:

1. validate the configuration (:mod:`cai.settings`),
2. build the laminate and ply blocks (:mod:`cai.layup`),
3. build the seeded impact damage (:mod:`cai.damage_map`),
4. size the damage-zone mesh and interface with Turon's equations
   (:mod:`cai.turon`),
5. mesh every ply block along its fibre direction (:mod:`cai.model`,
   :mod:`cai.mesh2d`),
6. write ``<job>.inp`` (:mod:`cai.inp_writer`) and the cluster files
   (:mod:`cai.jobfiles`) into a self-contained job folder.

Pure Python: runs anywhere (Python 2.7 or 3), including ``abaqus python``.
"""

from __future__ import division, print_function

import io
import os

from cai import compdam, jobfiles, settings, validation
from cai.damage_map import DamageMap
from cai.inp_writer import DeckWriter
from cai.layup import Laminate, membrane_properties
from cai.model import CAIMesh
from cai.turon import CohesiveSizing


class CAIJob(object):
    """One CAI analysis (one impact energy).

    Parameters
    ----------
    cfg : dict
        Configuration from :func:`cai.settings.load_config`.
    energy : float, optional
        Impact energy [J]; defaults to ``IMPACT["default_energy"]``.

    Attributes
    ----------
    name : str
        Job name.
    laminate : cai.layup.Laminate
    damage_map : cai.damage_map.DamageMap
    sizing : cai.turon.CohesiveSizing
    mesh : cai.model.CAIMesh
        Available after :meth:`build_mesh`.
    """

    def __init__(self, cfg, energy=None):
        self.cfg = cfg
        self.energy = float(cfg["IMPACT"]["default_energy"] if energy is None else energy)
        self.warnings = settings.validate(cfg)
        self.material_name, self.material = settings.active_material(cfg)
        self.name = settings.run_name(cfg, self.energy)
        self.center = settings.impact_center(cfg)
        self.laminate = Laminate.from_config(cfg["LAYUP"], self.material)
        self.damage_map = DamageMap(self.laminate, cfg["DAMAGE"], cfg["IMPACT"],
                                    self.center, self.energy)
        self.sizing = CohesiveSizing(self.material, cfg["COHESIVE_ZONE"],
                                     self.laminate.ply_thickness)
        self.nstatev = compdam.number_of_state_variables(cfg["COMPDAM"]["feature_flags"])
        self.mesh = None
        self.ic_sets = []
        self.contact_pairs = []
        self.category_counts = {}

    def build_mesh(self):
        """Mesh all ply blocks.

        Returns
        -------
        cai.model.CAIMesh
        """
        self.mesh = CAIMesh(self.cfg, self.laminate, self.damage_map,
                            self.sizing.element_size, self.center)
        self.warnings += self.mesh.warnings
        limit = self.cfg["MESH"].get("max_elements_warning")
        if limit and self.mesh.n_elements > limit:
            self.warnings.append("%d elements exceed max_elements_warning (%d); raise "
                                 "COHESIVE_ZONE min_element_size or reduce zone_margin"
                                 % (self.mesh.n_elements, limit))
        return self.mesh

    def report(self):
        """Build summary (laminate, CLT, damage, Turon sizing, mesh, NCAMP).

        Returns
        -------
        str
        """
        clt = membrane_properties(self.laminate, self.material["elastic"])
        parts = [
            "CAI job %s (impact %.1f J, material %s)" % (self.name, self.energy,
                                                        self.material_name),
            self.laminate.summary(),
            "CLT membrane: Ex = %.0f MPa, Ey = %.0f MPa, Gxy = %.0f MPa, nuxy = %.3f"
            % (clt["Ex"], clt["Ey"], clt["Gxy"], clt["nuxy"]),
            self.damage_map.summary(),
            self.sizing.summary(),
        ]
        if self.mesh is not None:
            parts.append(self.mesh.summary())
        if self.category_counts:
            parts.append("Seeded elements: %s" % ", ".join(
                "%s=%d" % (k, self.category_counts.get(k, 0))
                for k in ("pristine", "matrix", "fiber")))
        parts.append(validation.material_report(
            self.material_name, self.material, self.cfg["NCAMP_BASELINE"].get(self.material_name)))
        parts += ["WARNING: %s" % w for w in self.warnings]
        return "\n".join(parts)

    def write(self, output_dir):
        """Write the self-contained job folder.

        Parameters
        ----------
        output_dir : str
            Parent directory; the job goes to ``output_dir/<job name>``.

        Returns
        -------
        str
            Path of the job folder.
        """
        if self.mesh is None:
            self.build_mesh()
        folder = os.path.abspath(os.path.join(output_dir, self.name))
        if not os.path.isdir(folder):
            os.makedirs(folder)
        DeckWriter(self).write(os.path.join(folder, self.name + ".inp"))

        files = {
            "abaqus_v6.env": jobfiles.env_file(self.cfg, self.name),
            "run.sh": jobfiles.run_script(self.cfg, self.name),
            "postprocess.py": jobfiles.postprocess_script(
                self.cfg, self.name, self.laminate.thickness, self.energy),
            "summary.txt": self.report() + "\n",
        }
        sched = jobfiles.scheduler_script(self.cfg, self.name)
        if sched:
            files[sched[0]] = sched[1]
        for fname, text in files.items():
            with io.open(os.path.join(folder, fname), "w", encoding="utf-8") as fh:
                fh.write(u"" + text)
        for fname in ("run.sh", "submit.slurm", "submit.pbs"):
            path = os.path.join(folder, fname)
            if os.path.exists(path):
                os.chmod(path, 0o755)
        return folder


def run(cfg, energies, output_dir=None):
    """Generate one job folder per impact energy.

    Parameters
    ----------
    cfg : dict
        Loaded configuration.
    energies : list of float
        Impact energies [J].
    output_dir : str, optional
        Defaults to ``JOB["output_dir"]`` (relative to the current directory).

    Returns
    -------
    list of str
        Job folder paths.
    """
    out = output_dir or cfg["JOB"].get("output_dir", "jobs")
    folders = []
    for energy in energies:
        job = CAIJob(cfg, energy)
        job.build_mesh()
        folder = job.write(out)
        print(job.report())
        print("Job folder: %s" % folder)
        print("-" * 78)
        folders.append(folder)
    return folders
