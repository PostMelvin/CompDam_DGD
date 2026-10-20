"""Unit tests for the Abaqus-independent parts of the CAI generator."""

import math
import os

import pytest

from cai import ConfigurationError, ModelBuildError
from cai import compdam, settings, validation
from cai.damage_map import DamageMap, energy_scale, inside_ellipse
from cai.inp_patcher import InputDeckPatcher
from cai.layup import Laminate, expand_sequence, membrane_properties
from cai.mesh_plan import MeshPlan

HERE = os.path.dirname(os.path.abspath(__file__))
CONFIG = os.path.join(os.path.dirname(HERE), "config.py")


@pytest.fixture
def cfg():
    return settings.load_config(CONFIG)


# ----------------------------------------------------------------- compdam
def test_feature_flags_parse_and_leading_zeros():
    f = compdam.parse_feature_flags("111101")
    assert f == {"matrix_damage": 1, "shear_nonlinearity": 1, "fiber_tension": 1,
                 "fiber_compression": 1, "energy_output": 0, "friction": 1}
    assert compdam.parse_feature_flags(11000)["matrix_damage"] == 0
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
    assert len(props) == 40
    assert props[0] == 111101.0
    # Spot-check against tests/IM7-8552-C3D.inp positions (1-based)
    expected = {9: 171420.0, 10: 9080.0, 14: 62.3, 19: 199.8, 20: 0.925,
                29: 2326.2, 33: 1200.1, 36: 0.432, 40: 0.3}
    for pos, val in expected.items():
        assert props[pos - 1] == pytest.approx(val)
    assert props[20] == 0.0  # E3 blank -> transverse isotropy


def test_property_vector_missing_required():
    mat = {"elastic": {"E1": 1.0, "E2": 1.0, "G12": 1.0, "nu12": 0.3, "nu23": 0.4}}
    with pytest.raises(ConfigurationError, match="YT"):
        compdam.build_property_vector(mat, "100000")
    assert len(compdam.build_property_vector(mat, "000000")) == 40


def test_initial_state_and_format():
    s = compdam.initial_state(19, d2=1.0, alpha=0)
    assert s[0] == 1.0 and s[9] == 0 and s[10] == 1 and len(s) == 19
    block = compdam.format_initial_conditions([("CDM_MATRIX_CRACK", s)])
    lines = block.strip().splitlines()
    assert lines[0] == "*Initial Conditions, type=SOLUTION"
    assert lines[1].startswith("CDM_MATRIX_CRACK, 1, 0, 0")
    values = ",".join(lines[1:]).replace("CDM_MATRIX_CRACK,", "").split(",")
    assert len([v for v in values if v.strip()]) == 19
    assert lines[2].split(",")[3].strip() == "1"  # CDM_STATUS (var 11)


def test_initial_state_validation():
    with pytest.raises(ConfigurationError):
        compdam.initial_state(19, d2=1.5)
    with pytest.raises(ConfigurationError):
        compdam.initial_state(18, d1C=0.5)


def test_parameter_tables():
    assert compdam.format_parameter_tables({}, {}) == ""
    text = compdam.format_parameter_tables({"gamma_max": 4.0}, {"logLevel": 2})
    assert "*Table Collection, name=COMPDAM_PARAM" in text
    assert "gamma_max, 4" in text and "logLevel, 2" in text


# ------------------------------------------------------------------- layup
def test_laminate_from_config(cfg):
    lam = Laminate.from_config(cfg["LAYUP"], cfg["MATERIALS"]["AS4_8552"])
    assert len(lam.plies) == 24
    assert lam.thickness == pytest.approx(24 * 0.188)
    assert lam.is_symmetric() and lam.is_balanced()
    assert len(lam.blocks) == 24 and len(lam.interfaces) == 23
    assert lam.interfaces[0].z == pytest.approx(0.188)


def test_laminate_blocks():
    lam = Laminate(expand_sequence([45, 0, -45, 90], True), 0.2, blocks=[2, 4, 2])
    assert [len(b.plies) for b in lam.blocks] == [2, 4, 2]
    assert lam.interfaces[1].lower_ply.index == 6
    with pytest.raises(ConfigurationError):
        Laminate([0, 90], 0.2, blocks=[1, 2])


def test_quasi_isotropic_clt():
    lam = Laminate(expand_sequence([45, 0, -45, 90], True), 0.2)
    el = {"E1": 135000.0, "E2": 9600.0, "G12": 5100.0, "nu12": 0.32}
    c = membrane_properties(lam, el)
    assert c["Ex"] == pytest.approx(c["Ey"], rel=1e-9)
    assert c["Gxy"] == pytest.approx(c["Ex"] / (2 * (1 + c["nuxy"])), rel=1e-6)
    assert 45000 < c["Ex"] < 55000


# ------------------------------------------------------------------ damage
def test_energy_scale():
    assert energy_scale(30, 30, 5) == pytest.approx(1.0)
    assert energy_scale(4, 30, 5) == 0.0
    assert energy_scale(55, 30, 5) == pytest.approx(math.sqrt(2))


def test_inside_rotated_ellipse():
    assert inside_ellipse(10.0, 10.0, (0, 0), 15.0, 2.0, 45.0)
    assert not inside_ellipse(10.0, -10.0, (0, 0), 15.0, 2.0, 45.0)


def test_damage_map_cone_and_orientation(cfg):
    mat = cfg["MATERIALS"]["AS4_8552"]
    lam = Laminate.from_config(cfg["LAYUP"], mat)
    dm = DamageMap(lam, cfg["DAMAGE"], cfg["IMPACT"], (75.0, 50.0), 30.0)
    back, front = lam.interfaces[0], lam.interfaces[-1]
    assert dm.delamination_area(back) > dm.delamination_area(front)
    a, b, ang = dm.delamination_axes(back)
    assert ang == lam.interfaces[0].lower_ply.angle
    assert a / b == pytest.approx(cfg["DAMAGE"]["delamination"]["aspect_ratio"])
    assert dm.classify_element(lam.plies[0], 75.0, 50.0) == "matrix"
    assert dm.classify_element(lam.plies[0], 5.0, 5.0) == "pristine"
    off = DamageMap(lam, cfg["DAMAGE"], cfg["IMPACT"], (75.0, 50.0), 2.0)
    assert off.classify_element(lam.plies[0], 75.0, 50.0) == "pristine"


# --------------------------------------------------------------- mesh plan
def test_mesh_plan_seeds(cfg):
    plan = MeshPlan(cfg["GEOMETRY"], cfg["MESH"], cfg["BOUNDARY"], (75.0, 50.0))
    assert plan.x_cuts == [8.0, 45.0, 105.0, 142.0]
    assert plan.y_cuts == [5.0, 20.0, 50.0, 80.0, 95.0]
    fine = cfg["MESH"]["fine_size"]
    coarse = cfg["MESH"]["coarse_size"]
    assert plan.edge_seed((45, 0, 0), (105, 0, 0)) == ("size", fine)
    assert plan.edge_seed((8, 0, 0), (45, 0, 0)) == ("bias", fine, coarse, 45.0)
    assert plan.edge_seed((142, 0, 0), (105, 0, 0)) == ("bias", fine, coarse, 105.0)
    assert plan.edge_seed((0, 0, 0), (8, 0, 0)) == ("size", coarse)
    assert plan.edge_seed((0, 0, 0), (0, 0, 0.188)) == ("number", 1)
    with pytest.raises(ConfigurationError):
        plan.edge_seed((0, 0, 0), (1, 1, 0))


# ---------------------------------------------------------------- settings
def test_config_validates(cfg):
    warnings = settings.validate(cfg)
    assert isinstance(warnings, list)


def test_config_rejects_zone_outside_supports(cfg):
    cfg["GEOMETRY"]["refine_zone_length"] = 140.0
    with pytest.raises(ConfigurationError):
        settings.validate(cfg)


def test_run_name(cfg):
    assert settings.run_name(cfg, 22.5) == "CAI_QI_AS4_8552_22p5J"


# ------------------------------------------------------------- inp patcher
SAMPLE_INP = """*Heading
*Part, name=BLOCK_01
*End Part
**
** MATERIALS
**
*Material, name=AS4_8552
*Density
 1.59e-09,
*Depvar
     19,
*User Material, constants=40
111101.,0.,0.,0.,0.,0.,0.,0.
**
** INTERACTION PROPERTIES
**
*Surface Interaction, name=IP_FRICTION
** ----------------------------------------------------------------
**
** STEP: Compression
**
*Step, name=Compression, nlgeom=YES
*Dynamic, Explicit
, 0.02
*End Step
"""


def test_patcher_inserts_compdam_keywords():
    p = InputDeckPatcher(SAMPLE_INP)
    p.name_state_variables("AS4_8552", compdam.format_depvar_names(19))
    p.add_characteristic_length("AS4_8552", 3)
    p.insert_before_first_step(compdam.format_initial_conditions(
        [("CDM_PRISTINE", compdam.initial_state(19))]))
    lines = p.text.splitlines()
    i_dep = lines.index("*Depvar")
    assert lines[i_dep + 1].strip() == "19,"
    assert lines[i_dep + 2] == "  1, CDM_d2"
    assert lines[i_dep + 20] == " 19, CDM_d1C"
    i_mat = lines.index("*Material, name=AS4_8552")
    assert lines[i_mat + 1].startswith("*Characteristic Length, definition=USER, components=3")
    i_ic = lines.index("*Initial Conditions, type=SOLUTION")
    i_step = [i for i, l in enumerate(lines) if l.startswith("*Step")][0]
    assert i_ic < i_step
    assert lines[i_step - 1] == "**" and lines[i_step - 2] == "** STEP: Compression"
    with pytest.raises(ModelBuildError):
        p.name_state_variables("AS4_8552", compdam.format_depvar_names(19))


def test_patcher_missing_material():
    with pytest.raises(ModelBuildError):
        InputDeckPatcher(SAMPLE_INP).add_characteristic_length("NOPE")


# -------------------------------------------------------------- validation
def test_validation_reports(cfg):
    txt = validation.material_report("AS4_8552", cfg["MATERIALS"]["AS4_8552"],
                                     cfg["NCAMP_BASELINE"]["AS4_8552"])
    assert "XT" in txt and "NOT verified" in txt
    assert validation.cai_strength(-90000.0, 100.0, 4.5) == pytest.approx(200.0)
    rep = validation.cai_report(180.0, cfg["NCAMP_BASELINE"]["AS4_8552"], 30.0, 4.5)
    assert "delta" in rep
