"""
Post-processing of the Abaqus input deck written by CAE.

Abaqus/CAE has no GUI or scripting support for several keywords CompDam
needs, so they are added to the ``.inp`` text after ``Job.writeInput()``:

* named state variables under ``*Depvar``,
* ``*Characteristic Length, definition=USER`` (activates CompDam's
  VUCHARLENGTH),
* ``*Initial Conditions, type=SOLUTION`` (pristine and seeded damage states),
* the ``COMPDAM_PARAM`` parameter table collection.

Editing the text gives a transparent result you can diff, and it avoids
fragile edits by index to ``model.keywordBlock``.
"""

from __future__ import division, print_function

import io
import re

from cai import ModelBuildError


class InputDeckPatcher(object):
    """Insert CompDam-specific keywords into an Abaqus input deck.

    Parameters
    ----------
    text : str
        Contents of the ``.inp`` file.

    Examples
    --------
    >>> patcher = InputDeckPatcher.from_file("CAI.inp")   # doctest: +SKIP
    >>> patcher.name_state_variables("AS4_8552", ["  1, CDM_d2"])  # doctest: +SKIP
    >>> patcher.write("CAI.inp")   # doctest: +SKIP
    """

    _KEYWORD = re.compile(r"^\*(?!\*)")

    def __init__(self, text):
        self.lines = text.splitlines()

    @classmethod
    def from_file(cls, path):
        """Read an input deck from disk.

        Parameters
        ----------
        path : str
            Path to the ``.inp`` file.

        Returns
        -------
        InputDeckPatcher
        """
        with io.open(path, "r", encoding="utf-8", errors="replace") as f:
            return cls(f.read())

    @property
    def text(self):
        """Current deck contents as a single string."""
        return "\n".join(self.lines) + "\n"

    def write(self, path):
        """Write the patched deck to disk.

        Parameters
        ----------
        path : str
            Destination path.
        """
        with io.open(path, "w", encoding="utf-8") as f:
            f.write(u"" + self.text)

    # ------------------------------------------------------------- helpers
    def _is_keyword(self, line):
        return bool(self._KEYWORD.match(line))

    def _find_material(self, material_name):
        pattern = re.compile(r"^\*Material\s*,\s*name\s*=\s*%s\s*$" % re.escape(material_name),
                             re.IGNORECASE)
        for i, line in enumerate(self.lines):
            if pattern.match(line.strip()):
                return i
        raise ModelBuildError("material %r not found in the input deck" % material_name)

    def _material_end(self, start):
        """Index of the first line after the material block starting at ``start``."""
        for i in range(start + 1, len(self.lines)):
            line = self.lines[i].strip()
            if self._is_keyword(line):
                head = line.split(",")[0].strip().lower()
                if head in ("*material", "*step", "*boundary", "*initial conditions",
                            "*surface interaction", "*amplitude", "*contact pair",
                            "*end assembly", "*part", "*assembly"):
                    return i
            elif line.startswith("** ") and line[3:].strip().isupper() and len(line) > 4:
                # CAE section banners such as "** INTERACTION PROPERTIES"
                return i
        return len(self.lines)

    def _find_in_block(self, start, end, keyword):
        kw = keyword.lower()
        for i in range(start, end):
            line = self.lines[i].strip()
            if self._is_keyword(line) and line.split(",")[0].strip().lower() == kw:
                return i
        return None

    # ------------------------------------------------------------ patches
    def name_state_variables(self, material_name, name_lines):
        """Add state-variable names after the ``*Depvar`` count line.

        Parameters
        ----------
        material_name : str
            Name of the CompDam material.
        name_lines : list of str
            Lines from :func:`cai.compdam.format_depvar_names`.

        Raises
        ------
        ModelBuildError
            If the material or its ``*Depvar`` keyword is missing, or names
            were already present.
        """
        start = self._find_material(material_name)
        end = self._material_end(start)
        idx = self._find_in_block(start, end, "*depvar")
        if idx is None:
            raise ModelBuildError("*Depvar not found in material %r" % material_name)
        count_line = idx + 1
        nxt = self.lines[count_line + 1].strip() if count_line + 1 < len(self.lines) else ""
        if re.match(r"^\d+\s*,\s*[A-Za-z]", nxt):
            raise ModelBuildError("state variable names already present for %r" % material_name)
        self.lines[count_line + 1:count_line + 1] = list(name_lines)

    def add_characteristic_length(self, material_name, components=3):
        """Add ``*Characteristic Length, definition=USER`` to a material.

        Parameters
        ----------
        material_name : str
            Name of the CompDam material.
        components : int, optional
            3 for C3D8R; 6 with the random fibre misalignment options.
        """
        start = self._find_material(material_name)
        end = self._material_end(start)
        if self._find_in_block(start, end, "*characteristic length") is not None:
            return
        line = "*Characteristic Length, definition=USER, components=%d" % int(components)
        self.lines.insert(start + 1, line)

    def insert_before_first_step(self, block):
        """Insert model-data keywords just before the first ``*Step``.

        Parameters
        ----------
        block : str
            Keyword text (may span several lines). Empty strings are ignored.

        Raises
        ------
        ModelBuildError
            If no ``*Step`` keyword exists.
        """
        if not block:
            return
        for i, line in enumerate(self.lines):
            if self._is_keyword(line.strip()) and line.strip().lower().startswith("*step"):
                # Keep the CAE "** STEP: ..." banner directly above *Step.
                j = i
                while j > 0 and self.lines[j - 1].startswith("**"):
                    j -= 1
                self.lines[j:j] = block.rstrip("\n").splitlines()
                return
        raise ModelBuildError("no *Step keyword found in the input deck")

    def has_keyword(self, keyword):
        """Return True if ``keyword`` (e.g. ``"*Initial Conditions"``) occurs."""
        kw = keyword.lower()
        return any(l.strip().lower().startswith(kw) for l in self.lines)
