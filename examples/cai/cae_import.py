"""
Optional: import a generated CAI deck into Abaqus/CAE for inspection.

Not needed to run the analysis. Usage::

    abaqus cae noGUI=cae_import.py -- jobs/CAI_QI_AS4_8552_30J/CAI_QI_AS4_8552_30J.inp

Creates ``<job>.cae`` next to the input file. Open it in CAE to check the
fibre-aligned mesh, sets, orientations and contact assignments. Abaqus/CAE
may warn about the CompDam-specific keywords it does not support; the
analysis always runs from the ``.inp`` file, not from the ``.cae``.
"""

from __future__ import print_function

import os
import sys

from abaqus import mdb


def main():
    """Import the input deck given on the command line and save a .cae file."""
    args = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else []
    if not args:
        print("usage: abaqus cae noGUI=cae_import.py -- <job>.inp")
        return
    inp = os.path.abspath(args[0])
    name = os.path.splitext(os.path.basename(inp))[0]
    mdb.ModelFromInputFile(name=name, inputFileName=inp)
    if "Model-1" in mdb.models.keys() and len(mdb.models.keys()) > 1:
        del mdb.models["Model-1"]
    cae = os.path.splitext(inp)[0] + ".cae"
    mdb.saveAs(pathName=cae)
    print("Saved %s" % cae)


main()
