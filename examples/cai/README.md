# Compression After Impact (CAI) with CompDam_DGD

A modular, object-oriented Abaqus/CAE Python generator for an **ASTM D7137 Compression-After-Impact** model of a carbon/epoxy laminate. Intralaminar damage uses the NASA **CompDam_DGD** VUMAT; delaminations use surface-based cohesive contact.

The model starts from an idealised impact damage state (seeded delaminations and matrix cracks) and simulates the quasi-static compression phase in Abaqus/Explicit.

```
examples/cai/
├── config.py            # ALL inputs: geometry, layup, NCAMP-based materials, mesh, damage, BCs, solver, job
├── cai_model.py         # entry point: abaqus cae noGUI=cai_model.py -- [options]
├── cai_postprocess.py   # entry point: abaqus python cai_postprocess.py -- <job>.odb
├── cai/
│   ├── settings.py          # load / validate / override the config         (pure Python)
│   ├── compdam.py           # 40-constant property vector, flags, SDVs, ICs (pure Python)
│   ├── layup.py             # plies, ply blocks, interfaces, CLT checks     (pure Python)
│   ├── mesh_plan.py         # partition planes + edge-seeding rules         (pure Python)
│   ├── damage_map.py        # delamination / matrix-crack footprints        (pure Python)
│   ├── inp_patcher.py       # adds CompDam keywords CAE cannot write        (pure Python)
│   ├── validation.py        # NCAMP comparison, CAI strength                (pure Python)
│   ├── geometry.py          # parts, partitions, ply-aligned C3D8R mesh, assembly   (Abaqus)
│   ├── material_section.py  # VUMAT material, section, ply orientations             (Abaqus)
│   ├── damage.py            # damage element sets, bonded/delaminated surfaces      (Abaqus)
│   ├── interactions.py      # cohesive contact pairs or ties                        (Abaqus)
│   ├── loading.py           # explicit step, D7137 fixture BCs, outputs             (Abaqus)
│   ├── orchestrator.py      # CAIModel: build -> write .inp -> patch -> submit      (Abaqus)
│   └── postprocess.py       # ODB -> load/displacement CSV, CAI strength, KE/IE     (Abaqus)
├── tests/test_pure.py   # pytest suite for the pure-Python modules
├── docs/                # Sphinx (numpydoc via napoleon; Abaqus modules mocked)
└── pyproject.toml       # uv project (tests and docs only)
```

## Prerequisites

1. **Abaqus 2023 or newer** with a working Fortran setup for user subroutines (CompDam ≥ 2.7 needs Abaqus 2023+). See [`usersubroutine-prerequisites.md`](../../usersubroutine-prerequisites.md).
2. **CompDam set up**: from the repository root run `python setup.py`, or rename `for/version.for.nogit` to `for/version.for`.
3. An `abaqus_v6.env` with the CompDam compiler flags (`-free -fpp -mkl=sequential`). The build copies [`tests/abaqus_v6.env`](../../tests/abaqus_v6.env) into the working directory if none exists there.
4. Optional: [`uv`](https://docs.astral.sh/uv/) for the unit tests and documentation.

## Running

All commands run from `examples/cai/`.

```bash
# 0. Check the configuration without Abaqus (any Python 3)
python cai_model.py -- --check

# 1. Build the model, save CAI_QI.cae and write the patched input deck (no solve)
abaqus cae noGUI=cai_model.py

# 2. Other options
abaqus cae noGUI=cai_model.py -- --energy 40               # one impact energy [J]
abaqus cae noGUI=cai_model.py -- --sweep                   # every energy in IMPACT["energies"]
abaqus cae noGUI=cai_model.py -- --material IM7_8552
abaqus cae noGUI=cai_model.py -- --interface tie           # bonded blocks (debug/baseline)
abaqus cae noGUI=cai_model.py -- --submit --cpus 8         # build, write, submit, wait
```

Results go to `JOB["work_dir"]` (default `abaqus_work/`). To solve the written deck yourself:

```bash
cd abaqus_work
abaqus job=CAI_QI_AS4_8552_30J user=../../../for/CompDam_DGD.for double=both cpus=4 interactive
```

`double=both` is required: CompDam only runs in double precision. Jobs submitted from the script use `DOUBLE_PLUS_PACK` for the same reason.

From an open CAE session, use *File → Run Script… → cai_model.py*. Arguments cannot be passed there, so edit `config.py` instead.

### Post-processing

```bash
cd abaqus_work
abaqus python ../cai_postprocess.py -- CAI_QI_AS4_8552_30J.odb --energy 30
```

This writes `<job>_load_displacement.csv` and prints:

* the peak load,
* the CAI strength σ_CAI = P_max / (w·h) per ASTM D7137,
* the difference from the NCAMP value,
* the quasi-static check: max ALLKE/ALLIE should stay below about 5 %.

## What the model contains

| Aspect | Implementation |
|---|---|
| Geometry | 150 × 100 mm plate (ASTM D7137). X = loading direction, Z = thickness; ply 1 at z = 0 (back face), impact face at z = t |
| Layup | `[45/0/-45/90]_3s` with 0.188 mm plies: 24 plies, 4.51 mm thick (D7137 asks for 4–6 mm) |
| Ply blocks | One part per block (default: one ply per block, so 24 parts and 23 interfaces). Change with `plies_per_block` or an explicit `blocks` list |
| Elements | C3D8R, Explicit library, **enhanced hourglass control**, **distortion control** (CompDam requirement), average-strain kinematic split |
| Ply-aligned through-thickness mesh | Cells are cut at every ply boundary, so element faces coincide with ply interfaces; `elements_per_ply` (default 1) elements through each ply |
| Ply orientation | One rotated rectangular CSYS per ply (1 = fibre, 3 = stacking direction) assigned as a material orientation |
| In-plane mesh | Structured hex mesh: uniform `fine_size` in the refined zone around the impact, single-bias grading to `coarse_size` toward the clamps; identical in every block, so interface nodes coincide |
| CompDam material | `*User Material, constants=40` built from `config.MATERIALS` in CompDam order; `*Depvar` (19 SDVs for flags `111101`) with CompDam names; `*Characteristic Length, definition=USER, components=3` |
| Impact damage | Conical, peanut-shaped delaminations (ellipse along the fibre direction of the ply below each interface), circular matrix-crack zones (`CDM_d2 = 1`, `CDM_alpha = 0`), optional fibre damage. Radii scale with √((E − E_th)/(E_ref − E_th)) |
| Interfaces | Bonded regions: cohesive contact (quadratic-traction initiation, BK energy evolution) + friction. Seeded delaminations: frictional hard contact only. Or ties (`--interface tie`) |
| Fixture BCs | Loaded end kinematically coupled (U1) to `RP_LOAD`; reacting end U1 = 0; anti-rotation end strips U3 = 0 on both faces; side knife edges U3 = 0 on both faces; U2 = 0 on the centre line at x = 0 |
| Step | Explicit dynamic, NLGEOM on (CompDam requirement), smooth-step displacement ramp, semi-automatic mass scaling to a target Δt, bulk viscosity |
| Outputs | S, LE, U, RF, SDV, CSTRESS, CSDMG fields; energy histories; RP U1/RF1 |

### Keywords added after `writeInput()`

Abaqus/CAE cannot write some keywords CompDam needs, so `cai/inp_patcher.py` adds them to the `.inp` text. The patched deck is plain text you can diff:

* the state-variable names under `*Depvar`,
* `*Characteristic Length, definition=USER, components=3` (calls CompDam's VUCHARLENGTH),
* `*Initial Conditions, type=SOLUTION` for the `CDM_PRISTINE`, `CDM_MATRIX_CRACK` and `CDM_FIBER_DAMAGE` sets (`CDM_STATUS = 1` everywhere),
* the `COMPDAM_PARAM` table collection for non-default model parameters (`COMPDAM["parameters_float"/"parameters_int"]`).

## NCAMP validation workflow

`config.py` holds two kinds of material data on purpose:

* **`MATERIALS`**: the values fed to CompDam (lamina elastic and strength data plus fracture properties, which NCAMP does not report).
* **`NCAMP_BASELINE`**: the NCAMP reference values. These are the lamina RTD properties and the laminate CAI strength at 1500 in-lbf/in (6.67 J/mm).

Every build prints a table of model inputs against the NCAMP lamina values with the % difference. The post-processor compares the predicted CAI strength with `NCAMP_BASELINE[...]["cai"]["cai_strength"]`.

**Before you use the numbers quantitatively:**

1. Open the NCAMP report revision you validate against (the AS4/8552 tape report is NCP-RP-2010-002 / CAM-RP-2010-002; see the NIAR/NCAMP website). Enter its values in `NCAMP_BASELINE`, choosing the condition (RTD) and basis (mean vs. B-basis, normalised vs. as-measured), then set `"verified": True`. The shipped values are representative and flagged `NOT verified`.
2. Make `MATERIALS["AS4_8552"]` consistent with that baseline. Calibrate the fibre fracture toughnesses (`GXT`, `GXC`); they are currently the IM7/8552 values from CompDam's `tests/IM7-8552-C3D.inp`.
3. Calibrate the damage seed with your own C-scans. Use the projected damage area and depth distribution at the impact energies in `IMPACT["energies"]` to set `radius_front`, `radius_back` and `aspect_ratio`.
4. Run `--sweep` and compare CAI strength vs. impact energy with the NCAMP/test data.

The CLT membrane modulus printed at build time gives a quick check of the initial slope of the load-displacement curve: K ≈ E_x·w·t / L_free.

## Practical notes and limitations

* **Model size.** The defaults (0.8 mm fine mesh, one ply per block) give roughly 10 k elements per ply and about 235 k C3D8R elements in total. For quick checks, set `plies_per_block = 2` (or `blocks = [...]`) and `fine_size = 1.0`.
* **Snap-back limit.** `validate()` warns when an element exceeds 2EG/σ² for a mode. The coarse outer mesh exceeds the matrix limit on purpose; damage must stay inside the refined zone, so enlarge it if damage reaches its edge.
* **In-plane fibre alignment.** The mesh is orthogonal in-plane, so it is fibre-aligned for 0° and 90° plies only. CompDam's matrix-crack DGD and max-strain fibre models work on non-aligned meshes, because VUCHARLENGTH computes the characteristic lengths in material directions. The **fibre-kinking (FKT) models with random misalignment (φ0 = 1/2/3) need a uniform, fibre-aligned mesh**, which this generator does not build. It warns if FKT flags are set.
* **Ties and coupling.** With `--interface tie`, tie secondary nodes on the loaded end are also in the kinematic coupling. Abaqus may report over-constraints there. Tie mode is meant for debugging and baselines.
* **Hourglass control.** `ENHANCED` is used as requested. If your Abaqus release rejects it with the VUMAT, set `MESH["hourglass_control"] = "RELAX_STIFFNESS"`.
* **Quasi-static loading.** Check ALLKE/ALLIE after each run. If it is too high, lengthen `STEP["time_period"]` or reduce mass scaling.
* **Element deletion** is off, following the CompDam guidance (`COMPDAM["element_deletion"]`). `terminate_on_no_convergence = 0` deletes the rare element where DGD fails to converge instead of aborting the run.

## Tests and documentation (uv)

```bash
uv run pytest                                                    # pure-Python unit tests
uv run --group docs sphinx-build -b html docs docs/_build/html   # API docs (Abaqus mocked)
```

The Abaqus builders follow Python 2.7/3 compatible syntax because Abaqus 2023 ships Python 2.7 and Abaqus 2024+ ships Python 3.
