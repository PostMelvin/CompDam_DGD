"""Unit tests for the standalone CAI generator (no Abaqus required)."""

import math
import os
import re

import pytest

from cai import ConfigurationError
from cai import compdam, jobfiles, settings, validation
from cai.builder import CAIJob
from cai.damage_map import DamageMap, energy_scale, inside_ellipse
from cai.layup import Laminate, auto_blocks, expand_sequence, membrane_properties
from cai.mesh2d import Shape, count_cells, graded_sizes, mesh_ply, normalize_angle
from cai.turon import CohesiveSizing, cohesive_zone_length, reduced_strength

HERE = os.path.dirname(os.path.abspath(__file__))
CONFIG = os.path.join(os.path.dirname(HERE), "config.py")


@pytest.fixture
def cfg():
    return settings.load_config(CONFIG)


@pytest.fixture
def small_cfg(cfg):
    """Fast variant: 8 plies, coarse damage-zone mesh."""
    cfg["LAYUP"]["sequence"] = [45.0, 0.0, -45.0, 90.0]
    cfg["COHESIVE_ZONE"]["min_element_size"] = 1.0
    cfg["COHESIVE_ZONE"]["max_element_size"] = 1.0
    cfg["MESH"]["coarse_size"] = 6.0
    cfg["COMPDAM"]["library_dir"] = "/apps/compdam/lib"
    return cfg


# ----------------------------------------------------------------- compdam
def test_feature_flags_parse_and_leading_zeros():
    f = compdam.parse_feature_flags("111101")
    assert f == {"matrix_damage": 1, "shear_nonlinearity": 1, "fiber_tension": 1,
                 "fiber_compression": 1, "energy_output": 0, "friction": 1}
    assert compdam.parse_feature_flags(11000)["shear_nonlinearity"] == 1


@pytest.mark.parametrize("bad", ["2abc", "1234567", "200000", "100200"])
def test_feature_flags_reject_invalid(bad):
    with pytest.raises(ConfigurationError):
        compdam.parse_feature_flags(bad)


@pytest.mark.parametrize("flags,n", [("100000", 18), ("101100", 19), ("130000", 21),
                                     ("101300", 25), ("150000", 33)])
def test_number_of_state_variables(flags, n):
    assert compdam.number_of_state_variables(flags) == n


def test_property_vector_matches_compdam_order(cfg):
    mat = cfg["MATERIALS"]["IM7_8552"]
    props = compdam.build_property_vector(mat, "111101")
    assert len(props) == 40 and props[0] == 111101.0
    expected = {9: 171420.0, 10: 9080.0, 14: 62.3, 19: 199.8, 20: 0.925,
                29: 2326.2, 33: 1200.1, 36: 0.432, 40: 0.3}
    for pos, val in expected.items():
        assert props[pos - 1] == pytest.approx(val)


def test_user_material_lines_keep_blanks(cfg):
    lines = compdam.format_user_material(cfg["MATERIALS"]["IM7_8552"], "111101")
    assert len(lines) == 5
    assert lines[0].startswith("111101, , ,")
    assert lines[2].split(",")[4].strip() == ""  # E3 blank
    assert sum(len(l.rstrip(",").split(",")) for l in lines) == 40


def test_property_vector_missing_required():
    mat = {"elastic": {"E1": 1.0, "E2": 1.0, "G12": 1.0, "nu12": 0.3, "nu23": 0.4}}
    with pytest.raises(ConfigurationError, match="YT"):
        compdam.build_property_vector(mat, "100000")


def test_initial_state_and_format():
    s = compdam.initial_state(19, d2=1.0, alpha=0)
    assert s[0] == 1.0 and s[9] == 0 and s[10] == 1
    lines = compdam.format_initial_conditions([("CDM_MATRIX_CRACK", s)]).strip().splitlines()
    assert lines[0] == "*Initial Conditions, type=SOLUTION"
    assert lines[1].startswith("CDM_MATRIX_CRACK, 1, 0, 0")
    assert lines[2].split(",")[3].strip() == "1"  # CDM_STATUS (var 11)
    with pytest.raises(ConfigurationError):
        compdam.initial_state(19, d2=1.5)


# ------------------------------------------------------------------- layup
def test_auto_blocks_merge_equal_angles():
    angles = expand_sequence([45, 0, -45, 90], True)
    assert auto_blocks(angles) == [1, 1, 1, 2, 1, 1, 1]
    assert auto_blocks([90, -90, 0]) == [2, 1]  # -90 and 90 are the same direction


def test_laminate_from_config(cfg):
    lam = Laminate.from_config(cfg["LAYUP"], cfg["MATERIALS"]["AS4_8552"])
    assert len(lam.plies) == 24 and lam.thickness == pytest.approx(24 * 0.188)
    assert len(lam.blocks) == 23 and len(lam.interfaces) == 22
    assert lam.blocks[11].uniform_angle() == 90.0
    assert lam.is_symmetric() and lam.is_balanced()


def test_quasi_isotropic_clt():
    lam = Laminate(expand_sequence([45, 0, -45, 90], True), 0.2)
    c = membrane_properties(lam, {"E1": 135000.0, "E2": 9600.0, "G12": 5100.0, "nu12": 0.32})
    assert c["Ex"] == pytest.approx(c["Ey"], rel=1e-9)
    assert c["Gxy"] == pytest.approx(c["Ex"] / (2 * (1 + c["nuxy"])), rel=1e-6)


# ------------------------------------------------------------------- turon
def test_turon_equations():
    l = cohesive_zone_length(9600.0, 0.28, 62.3)
    assert l == pytest.approx(9 * math.pi / 32 * 9600 * 0.28 / 62.3 ** 2)
    # Reduced strength gives exactly N_e elements of size l_e
    tau = reduced_strength(9600.0, 0.28, 3, 0.5)
    assert cohesive_zone_length(9600.0, 0.28, tau) == pytest.approx(3 * 0.5)


def test_cohesive_sizing_floor_and_reduction(cfg):
    mat = cfg["MATERIALS"]["AS4_8552"]
    sz = CohesiveSizing(mat, cfg["COHESIVE_ZONE"], 0.188)
    assert sz.turon_size == pytest.approx(min(sz.l_cz.values()) / 3)
    assert sz.element_size == pytest.approx(0.3)  # min_element_size floor
    assert sz.strength_reduced
    assert sz.elements_in_zone() == pytest.approx(3.0, rel=1e-6)
    assert sz.penalty == pytest.approx(50 * 9600 / 0.188)
    cz = dict(cfg["COHESIVE_ZONE"], min_element_size=None)
    full = CohesiveSizing(mat, cz, 0.188)
    assert full.element_size == pytest.approx(full.turon_size)
    assert not full.strength_reduced


# ------------------------------------------------------------------ damage
def test_energy_scale():
    assert energy_scale(30, 30, 5) == pytest.approx(1.0)
    assert energy_scale(4, 30, 5) == 0.0
    assert energy_scale(55, 30, 5) == pytest.approx(math.sqrt(2))


def test_inside_rotated_ellipse():
    assert inside_ellipse(10.0, 10.0, (0, 0), 15.0, 2.0, 45.0)
    assert not inside_ellipse(10.0, -10.0, (0, 0), 15.0, 2.0, 45.0)


def test_damage_map_cone(cfg):
    lam = Laminate.from_config(cfg["LAYUP"], cfg["MATERIALS"]["AS4_8552"])
    dm = DamageMap(lam, cfg["DAMAGE"], cfg["IMPACT"], (75.0, 50.0), 30.0)
    assert dm.delamination_area(lam.interfaces[0]) > dm.delamination_area(lam.interfaces[-1])
    assert dm.classify_element(lam.plies[0], 75.0, 50.0) == "matrix"
    assert dm.classify_element(lam.plies[0], 5.0, 5.0) == "pristine"


# ------------------------------------------------------------------ mesh2d
def test_spacing_helpers():
    assert normalize_angle(135) == -45 and normalize_angle(-90) == 90
    sizes = graded_sizes(40.0, 15, 0.4, 4.0)
    assert sum(sizes) == pytest.approx(40.0)
    assert sizes[0] == pytest.approx(0.4, rel=0.05) and max(sizes) <= 4.0 + 1e-9
    assert all(b >= a - 1e-12 for a, b in zip(sizes, sizes[1:]))
    assert count_cells(0.0, 1, 4, 1.2) == 0


@pytest.mark.parametrize("angle", [45.0, -45.0, 0.0, 90.0, 60.0, -30.0])
def test_ply_mesh_valid_and_fiber_aligned(angle):
    shapes = [Shape((75.0, 50.0), 20.0, 12.0, angle), Shape((75.0, 50.0), 10.0)]
    m = mesh_ply(150.0, 100.0, angle, shapes, 0.6, 4.0, 1.2, (8.0, 142.0), (5.0, 95.0), 5.0)
    q = m.quality()
    assert q["min_jacobian"] > 0.0
    assert q["area"] == pytest.approx(150.0 * 100.0, rel=1e-9)
    assert q["max_aspect"] < 25.0
    f = (math.cos(math.radians(angle)), math.sin(math.radians(angle)))
    for qi in m.zone_quads:
        pts = [m.nodes[i] for i in m.quads[qi]]
        # one pair of opposite edges is parallel to the fibres
        e = [(pts[(k + 1) % 4][0] - pts[k][0], pts[(k + 1) % 4][1] - pts[k][1]) for k in range(4)]
        cross = [abs(v[0] * f[1] - v[1] * f[0]) / math.hypot(*v) for v in e]
        assert min(cross[0], cross[1]) < 1e-9
    # knife-edge lines are exact node rows
    assert any(abs(y - 5.0) < 1e-9 for x, y in m.nodes)
    assert any(abs(y - 95.0) < 1e-9 for x, y in m.nodes)
    # the zone contains every shape centre and its ellipse extremities
    xs = [m.nodes[i][0] for qi in m.zone_quads for i in m.quads[qi]]
    assert min(xs) < 75.0 - 10.0 and max(xs) > 75.0 + 10.0


def test_zone_rejects_impossible_fit():
    with pytest.raises(ConfigurationError):
        mesh_ply(150.0, 100.0, 45.0, [Shape((75.0, 50.0), 5.0)], 0.5, 4.0, 1.2,
                 (74.0, 76.0), (5.0, 95.0), 1.0)


# ------------------------------------------------------------ full model
def test_small_job_mesh_and_sets(small_cfg):
    job = CAIJob(small_cfg, 30.0)
    mesh = job.build_mesh()
    assert len(mesh.blocks) == 7
    for _, angle, ne, nz, q in mesh.quality():
        assert q["min_jacobian"] > 0 and nz > 0
    ns = mesh.node_sets()
    assert ns["KNIFE_EDGES"] and ns["END_SUPPORTS"] and ns["ANCHOR_Y"]
    cats = mesh.damage_categories()
    assert sum(len(v) for v in cats.values()) == mesh.n_elements
    surfaces = mesh.interface_surfaces()
    assert len(surfaces) == 6
    assert surfaces[0]["lower"]["delam"][1] and surfaces[0]["upper"]["bonded"][1]
    # positive hexahedron volume (check first element of each block)
    coords = {}
    for bm in mesh.blocks:
        for n in bm.iter_nodes():
            coords[n[0]] = n[1:]
    for bm in mesh.blocks:
        e = next(iter(bm.iter_elements()))
        p = [coords[i] for i in e[1:]]
        a = [p[1][k] - p[0][k] for k in range(3)]
        b = [p[3][k] - p[0][k] for k in range(3)]
        c = [p[4][k] - p[0][k] for k in range(3)]
        det = (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
               + a[2] * (b[0] * c[1] - b[1] * c[0]))
        assert det > 0


def test_written_job_folder_is_consistent(small_cfg, tmp_path):
    job = CAIJob(small_cfg, 30.0)
    folder = job.write(str(tmp_path))
    files = set(os.listdir(folder))
    assert {job.name + ".inp", "abaqus_v6.env", "run.sh", "postprocess.py",
            "summary.txt", "submit.slurm"} <= files
    text = open(os.path.join(folder, job.name + ".inp")).read()
    for kw in ("*Characteristic Length, definition=USER, components=3",
               "*Initial Conditions, type=SOLUTION", "*Contact Property Assignment",
               "*Cohesive Behavior", "*Coupling", "*User Material, constants=40",
               "*Section Controls, name=EC_COMPDAM, distortion control=YES",
               "*Fixed Mass Scaling", "*Table Collection, name=COMPDAM_PARAM"):
        assert kw in text, kw
    # every referenced set/surface is defined
    defined = set(re.findall(r"^\*(?:Elset|Nset|Node), (?:elset|nset)=(\w+)", text, re.M))
    defined |= set(re.findall(r"^\*Surface, type=ELEMENT, name=(\w+)", text, re.M))
    refs = set(re.findall(r"elset=(\w+), orientation", text))
    refs |= set(re.findall(r"surface=(\w+)", text))
    block = text.split("*Contact Property Assignment")[1].split("**")[0]
    for line in block.strip().splitlines()[1:]:
        refs |= set(t.strip() for t in line.split(",")[:2])
    ic = text.split("*Initial Conditions, type=SOLUTION")[1].split("*")[0]
    refs |= set(re.findall(r"^(CDM_\w+),", ic, re.M))
    bc = text.split("*Boundary\n")[1].split("*")[0]
    refs |= set(l.split(",")[0] for l in bc.strip().splitlines())
    assert refs <= defined, refs - defined
    # element connectivity references existing nodes only
    n_nodes = job.mesh.n_nodes + 1
    elems = re.findall(r"^\d+(?:, \d+){8}$", text, re.M)
    assert len(elems) == job.mesh.n_elements
    assert max(max(int(v) for v in e.split(",")[1:]) for e in elems[-50:]) <= n_nodes
    env = open(os.path.join(folder, "abaqus_v6.env")).read()
    assert "usub_lib_dir = '/apps/compdam/lib'" in env
    slurm = open(os.path.join(folder, "submit.slurm")).read()
    assert "double=both" in slurm and "user=" not in slurm
    compile(open(os.path.join(folder, "postprocess.py")).read(), "postprocess.py", "exec")


def test_source_compile_mode(small_cfg):
    small_cfg["COMPDAM"]["library_dir"] = None
    small_cfg["COMPDAM"]["user_subroutine"] = "/home/me/CompDam_DGD/for/CompDam_DGD.for"
    assert "compile_fortran" in jobfiles.env_file(small_cfg, "J")
    assert "user=/home/me/CompDam_DGD/for/CompDam_DGD.for" in jobfiles.abaqus_command(small_cfg, "J")
    small_cfg["CLUSTER"]["scheduler"] = "pbs"
    name, text = jobfiles.scheduler_script(small_cfg, "J")
    assert name == "submit.pbs" and "#PBS -l select=1:ncpus=16:mem=64gb" in text


# ---------------------------------------------------------------- settings
def test_config_validates(cfg):
    warnings = settings.validate(cfg)
    assert any("placeholder" in w for w in warnings)


def test_config_requires_library_or_source(cfg):
    cfg["COMPDAM"]["library_dir"] = None
    with pytest.raises(ConfigurationError):
        settings.validate(cfg)


def test_mixed_angle_blocks_need_orthogonal_mesh(small_cfg):
    small_cfg["LAYUP"]["blocks"] = [2, 2, 2, 2]
    with pytest.raises(ConfigurationError, match="mixes ply angles"):
        CAIJob(small_cfg, 30.0).build_mesh()
    small_cfg["MESH"]["fiber_aligned"] = False
    assert CAIJob(small_cfg, 30.0).build_mesh().n_elements > 0


def test_run_name(cfg):
    assert settings.run_name(cfg, 22.5) == "CAI_QI_AS4_8552_22p5J"


# -------------------------------------------------------------- validation
def test_validation_reports(cfg):
    txt = validation.material_report("AS4_8552", cfg["MATERIALS"]["AS4_8552"],
                                     cfg["NCAMP_BASELINE"]["AS4_8552"])
    assert "XT" in txt and "NOT verified" in txt
    assert validation.cai_strength(-90000.0, 100.0, 4.5) == pytest.approx(200.0)
