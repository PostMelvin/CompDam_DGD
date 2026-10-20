"""
Seed the pre-defined impact damage state on the meshed parts (Abaqus).

Uses :class:`cai.damage_map.DamageMap` to

* sort every element into a pristine / matrix-cracked / fibre-damaged
  category (element centroid test). The categories become assembly element
  sets that get different CompDam initial state variables,
* split each ply-block interface surface into a *bonded* region (cohesive
  contact) and a *delaminated* region (frictional contact only), using
  mesh-based element-face surfaces.

Must run inside the Abaqus/CAE kernel, after meshing.
"""

from __future__ import division, print_function

from cai import ModelBuildError

#: Assembly element-set names of the three damage categories.
CATEGORY_SETS = {
    "pristine": "CDM_PRISTINE",
    "matrix": "CDM_MATRIX_CRACK",
    "fiber": "CDM_FIBER_DAMAGE",
}

_FACE_KEYWORD = dict(("FACE%d" % i, "face%dElements" % i) for i in range(1, 7))


def _centroid(nodes):
    n = float(len(nodes))
    x = sum(nd.coordinates[0] for nd in nodes) / n
    y = sum(nd.coordinates[1] for nd in nodes) / n
    z = sum(nd.coordinates[2] for nd in nodes) / n
    return x, y, z


class DamageSeeder(object):
    """Create damage element sets and interface surfaces.

    Parameters
    ----------
    laminate : cai.layup.Laminate
        Laminate definition.
    damage_map : cai.damage_map.DamageMap
        Pre-defined damage description.
    tolerance : float
        Geometric tolerance [mm].

    Attributes
    ----------
    labels : dict
        ``{category: {instance_name: [element labels]}}``.
    interface_surfaces : dict
        ``{interface_name: {"lower": {"bonded": name, "delam": name},
        "upper": {...}}}``; a value of ``None`` means that region is empty.
    """

    def __init__(self, laminate, damage_map, tolerance):
        self.laminate = laminate
        self.dmap = damage_map
        self.tol = tolerance
        self.labels = dict((c, {}) for c in CATEGORY_SETS)
        self.interface_surfaces = {}

    # ------------------------------------------------------------ elements
    def classify_elements(self, parts):
        """Sort the elements of every ply into damage categories.

        Parameters
        ----------
        parts : dict
            ``{block_name: Part}``.
        """
        counts = dict((c, 0) for c in CATEGORY_SETS)
        for block in self.laminate.blocks:
            p = parts[block.name]
            per_cat = dict((c, []) for c in CATEGORY_SETS)
            for ply in block.plies:
                for elem in p.sets[ply.name].elements:
                    x, y, _ = _centroid(elem.getNodes())
                    per_cat[self.dmap.classify_element(ply, x, y)].append(elem.label)
            for cat, labels in per_cat.items():
                if labels:
                    self.labels[cat][block.instance_name] = labels
                    counts[cat] += len(labels)
        print("  damage seeding: %s" % ", ".join(
            "%s=%d" % (c, counts[c]) for c in ("pristine", "matrix", "fiber")))

    def create_assembly_sets(self, assembly):
        """Create the assembly-level element sets of the damage categories.

        Parameters
        ----------
        assembly : Assembly
            Root assembly (instances must exist).

        Returns
        -------
        list of tuple
            ``(set_name, category)`` for every non-empty category.

        Raises
        ------
        ModelBuildError
            If no element ends up in any category.
        """
        created = []
        for cat in ("pristine", "matrix", "fiber"):
            by_instance = self.labels[cat]
            if not by_instance:
                continue
            spec = tuple((inst, tuple(lbls)) for inst, lbls in sorted(by_instance.items()))
            assembly.SetFromElementLabels(name=CATEGORY_SETS[cat], elementLabels=spec)
            created.append((CATEGORY_SETS[cat], cat))
        if not created:
            raise ModelBuildError("no elements were classified; was the mesh generated?")
        return created

    # ------------------------------------------------------------ surfaces
    def split_interface_surfaces(self, parts):
        """Create bonded / delaminated mesh surfaces on both sides of each interface.

        Parameters
        ----------
        parts : dict
            ``{block_name: Part}``.
        """
        for itf in self.laminate.interfaces:
            lower = parts[itf.lower.name]
            upper = parts[itf.upper.name]
            self.interface_surfaces[itf.name] = {
                "lower": self._split_side(lower, itf, itf.lower_ply, itf.lower.thickness,
                                          "%s_LOWER" % itf.name),
                "upper": self._split_side(upper, itf, itf.upper_ply, 0.0,
                                          "%s_UPPER" % itf.name),
            }

    def _split_side(self, p, itf, ply, z_local, prefix):
        groups = {"bonded": {}, "delam": {}}
        for elem in p.sets[ply.name].elements:
            for face in elem.getElementFaces():
                nodes = face.getNodes()
                if any(abs(nd.coordinates[2] - z_local) > self.tol for nd in nodes):
                    continue
                x, y, _ = _centroid(nodes)
                key = "delam" if self.dmap.is_delaminated(itf, x, y) else "bonded"
                groups[key].setdefault(str(face.face), []).append(elem.label)
                break  # an element has one face on the interface
        names = {}
        for key, by_face in groups.items():
            if not by_face:
                names[key] = None
                continue
            kwargs = {}
            for face_id, labels in by_face.items():
                kwargs[_FACE_KEYWORD[face_id]] = p.elements.sequenceFromLabels(labels)
            name = "%s_%s" % (prefix, key.upper())
            p.Surface(name=name, **kwargs)
            names[key] = name
        if names["bonded"] is None and names["delam"] is None:
            raise ModelBuildError("no element faces found on %s of %s" % (prefix, p.name))
        return names
