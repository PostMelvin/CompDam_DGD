"""
Generate standalone CompDam_DGD Compression-After-Impact (CAI) jobs.

Plain Python (2.7 or 3.x, no packages, no Abaqus/CAE needed)::

    python cai_model.py                          # default energy -> jobs/<job>/
    python cai_model.py --energy 40
    python cai_model.py --sweep                  # every energy in IMPACT["energies"]
    python cai_model.py --library-dir /apps/compdam/lib --cpus 32
    python cai_model.py --check                  # validate + sizing summary, no mesh
    abaqus python cai_model.py --energy 30       # also works with Abaqus' Python

Each job folder holds ``<job>.inp``, ``abaqus_v6.env`` (usub_lib_dir -> your
compiled CompDam library), ``submit.slurm``/``submit.pbs``, ``run.sh``,
``postprocess.py`` and ``summary.txt``. Copy it to the cluster and submit.

All engineering inputs live in ``config.py``.
"""

from __future__ import division, print_function

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from cai import CAIError  # noqa: E402
from cai import settings  # noqa: E402


def parse_args(argv):
    """Parse command-line arguments (a leading ``--`` is ignored).

    Parameters
    ----------
    argv : list of str
        Arguments without the program name.

    Returns
    -------
    argparse.Namespace
    """
    argv = [a for a in argv if a != "--"]
    p = argparse.ArgumentParser(description="Generate CompDam CAI job folders.")
    p.add_argument("--config", default=os.path.join(HERE, "config.py"),
                   help="configuration file (default: config.py next to this script)")
    p.add_argument("--material", help="material preset, e.g. AS4_8552 or IM7_8552")
    p.add_argument("--energy", type=float, help="impact energy [J]")
    p.add_argument("--sweep", action="store_true",
                   help="one job per energy in IMPACT['energies']")
    p.add_argument("--interface", choices=("cohesive", "tie"), help="ply-block interface type")
    p.add_argument("--job-name", help="base job name")
    p.add_argument("--cpus", type=int, help="CPUs / Explicit domains in the cluster scripts")
    p.add_argument("--output-dir", help="parent folder of the job folders")
    p.add_argument("--library-dir", help="compiled CompDam library directory on the cluster")
    p.add_argument("--min-element-size", type=float,
                   help="floor on the damage-zone element size [mm] (COHESIVE_ZONE)")
    p.add_argument("--check", action="store_true",
                   help="validate the configuration and print the sizing summary only")
    return p.parse_args(argv)


def main(argv=None):
    """Entry point; returns the process exit code."""
    args = parse_args(sys.argv[1:] if argv is None else argv)
    try:
        cfg = settings.load_config(args.config)
        settings.apply_overrides(cfg, material=args.material, energy=args.energy,
                                 interface=args.interface, job_name=args.job_name,
                                 cpus=args.cpus, output_dir=args.output_dir,
                                 library_dir=args.library_dir,
                                 min_element_size=args.min_element_size)
        energies = list(cfg["IMPACT"]["energies"]) if args.sweep \
            else [cfg["IMPACT"]["default_energy"]]
        from cai.builder import CAIJob, run
        if args.check:
            for e in energies:
                print(CAIJob(cfg, e).report())
                print("-" * 78)
            return 0
        run(cfg, energies)
    except CAIError as exc:
        print("ERROR: %s" % exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
