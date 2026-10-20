"""
Post-process a finished CAI job and compare it with the NCAMP baseline.

Usage (from the Abaqus working directory)::

    abaqus python ../cai_postprocess.py -- CAI_QI_AS4_8552_30J.odb --energy 30

Writes ``<job>_load_displacement.csv`` and prints the peak load, the CAI
strength (ASTM D7137: P_max / (w h)), the NCAMP comparison and the
quasi-static energy check.
"""

from __future__ import division, print_function

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from cai import settings  # noqa: E402
from cai.layup import Laminate  # noqa: E402


def main():
    """Entry point; returns the process exit code."""
    argv = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else sys.argv[1:]
    p = argparse.ArgumentParser(prog="abaqus python cai_postprocess.py --")
    p.add_argument("odb", help="path to the .odb file")
    p.add_argument("--config", default=os.path.join(HERE, "config.py"))
    p.add_argument("--material", help="material preset used for the run")
    p.add_argument("--energy", type=float, help="impact energy used for the run [J]")
    args = p.parse_args(argv)

    from cai.postprocess import summarize  # needs odbAccess (Abaqus Python)

    cfg = settings.load_config(args.config)
    settings.apply_overrides(cfg, material=args.material)
    _, material = settings.active_material(cfg)
    thickness = Laminate.from_config(cfg["LAYUP"], material).thickness
    print(summarize(args.odb, cfg, energy=args.energy, thickness=thickness))
    return 0


if __name__ == "__main__":
    sys.exit(main())
