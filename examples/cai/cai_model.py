"""
Build the CompDam_DGD Compression-After-Impact (CAI) model in Abaqus/CAE.

Usage (from ``examples/cai``)::

    abaqus cae noGUI=cai_model.py                        # default energy, write .inp only
    abaqus cae noGUI=cai_model.py -- --energy 40         # one impact energy
    abaqus cae noGUI=cai_model.py -- --sweep             # all IMPACT["energies"]
    abaqus cae noGUI=cai_model.py -- --material IM7_8552 --submit --cpus 8
    abaqus cae noGUI=cai_model.py -- --check             # validate config only

Inside an open CAE session: *File > Run Script...* (arguments cannot be
passed; edit ``config.py`` instead).

All engineering inputs live in ``config.py``; this script only parses
command-line overrides and calls :func:`cai.orchestrator.run`.
"""

from __future__ import division, print_function

import argparse
import inspect
import os
import sys


def _script_dir():
    """Directory of this script (``__file__`` is not always defined under CAE)."""
    try:
        return os.path.dirname(os.path.abspath(__file__))
    except NameError:
        return os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))


HERE = _script_dir()
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from cai import CAIError  # noqa: E402
from cai import settings  # noqa: E402


def parse_args(argv):
    """Parse the arguments that follow ``--`` on the Abaqus command line.

    Parameters
    ----------
    argv : list of str
        Full ``sys.argv``.

    Returns
    -------
    argparse.Namespace
    """
    user = argv[argv.index("--") + 1:] if "--" in argv else []
    p = argparse.ArgumentParser(prog="abaqus cae noGUI=cai_model.py --",
                                description="Build the CompDam CAI model.")
    p.add_argument("--config", default=os.path.join(HERE, "config.py"),
                   help="configuration file (default: config.py next to this script)")
    p.add_argument("--material", help="material preset name, e.g. AS4_8552 or IM7_8552")
    p.add_argument("--energy", type=float, help="impact energy [J]")
    p.add_argument("--sweep", action="store_true",
                   help="build one model per energy in IMPACT['energies']")
    p.add_argument("--interface", choices=("cohesive", "tie"), help="ply-block interface type")
    p.add_argument("--job-name", help="base model/job name")
    p.add_argument("--cpus", type=int, help="number of CPUs/domains")
    p.add_argument("--submit", action="store_true", help="submit the job(s) after writing")
    p.add_argument("--check", action="store_true",
                   help="validate the configuration and print the summary, no Abaqus build")
    return p.parse_args(user)


def main(argv=None):
    """Entry point.

    Parameters
    ----------
    argv : list of str, optional
        Defaults to ``sys.argv``.

    Returns
    -------
    int
        Exit code (0 on success).
    """
    args = parse_args(sys.argv if argv is None else argv)
    try:
        cfg = settings.load_config(args.config)
        settings.apply_overrides(cfg, material=args.material, energy=args.energy,
                                 submit=True if args.submit else None,
                                 interface=args.interface, job_name=args.job_name,
                                 cpus=args.cpus)
        energies = list(cfg["IMPACT"]["energies"]) if args.sweep \
            else [cfg["IMPACT"]["default_energy"]]
        if args.check:
            for w in settings.validate(cfg):
                print("WARNING: %s" % w)
            print("Configuration OK for energies %s" % energies)
            return 0
        from cai.orchestrator import run  # needs the Abaqus kernel
        run(cfg, energies)
    except CAIError as exc:
        print("ERROR: %s" % exc)
        return 1
    return 0


if __name__ == "__main__":
    code = main()
    if code:
        sys.exit(code)
