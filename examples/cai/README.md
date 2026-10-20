# Compression After Impact (CAI) with CompDam_DGD

A **standalone** generator for an ASTM D7137 Compression-After-Impact model of a carbon/epoxy laminate. Intralaminar damage uses the NASA **CompDam_DGD** VUMAT; delaminations use surface-based cohesive contact.

* **No Abaqus/CAE needed to build the model.** Plain Python (2.7 or 3.x, standard library only) writes the complete Abaqus/Explicit input deck.
* **Fibre-aligned meshes.** Each ply block is meshed along its own fibre direction: parallelogram elements for ±45°, rectangles for 0/90°.
* **Cohesive-zone driven refinement.** The delamination-zone element size comes from Turon et al. (2007). Interface strengths are lowered automatically when you choose a coarser mesh.
* **Cluster-ready job folders.** Each folder holds `<job>.inp`, an `abaqus_v6.env` that points `usub_lib_dir` at your **pre-compiled CompDam library**, SLURM/PBS scripts and a standalone post-processor. Copy the folder to the cluster and submit; nothing refers back to this repository.

```
examples/cai/
├── config.py          # ALL inputs: geometry, layup, materials (NCAMP), mesh, Turon sizing, damage, BCs, step, cluster
├── cai_model.py       # CLI: python cai_model.py [--energy E | --sweep] [--library-dir DIR] ...
├── cae_import.py      # optional: import a generated deck into Abaqus/CAE for viewing
├── cai/
│   ├── settings.py    # load / validate / override the configuration
│   ├── compdam.py     # 40-constant property vector, feature flags, SDVs, initial states
│   ├── layup.py       # plies, ply blocks ("auto" = group equal angles), interfaces, CLT
│   ├── damage_map.py  # seeded delaminations (conical, peanut-shaped) and matrix cracks
│   ├── turon.py       # cohesive-zone length, element size, strength reduction, penalty
│   ├── mesh2d.py      # fibre-aligned structured quad mesh of one ply
│   ├── model.py       # 3-D extrusion, numbering, element/node sets, surfaces
│   ├── inp_writer.py  # Abaqus/Explicit input deck
│   ├── jobfiles.py    # abaqus_v6.env, submit.slurm / submit.pbs, run.sh, postprocess.py
│   ├── builder.py     # CAIJob: config -> mesh -> job folder
│   └── validation.py  # NCAMP comparison helpers
├── tests/test_pure.py # pytest suite (meshes, Turon, deck consistency, job files)
├── docs/              # Sphinx (numpydoc via napoleon)
└── pyproject.toml     # uv project (tests and docs only; the generator has no dependencies)
```

## Quick start

```bash
cd examples/cai

# 1. Point the jobs at your compiled CompDam library on the cluster (or set
#    COMPDAM["library_dir"] in config.py), then generate
python cai_model.py --library-dir /apps/compdam/lib            # -> jobs/CAI_QI_AS4_8552_30J/
python cai_model.py --sweep --library-dir /apps/compdam/lib    # one job per IMPACT["energies"]
python cai_model.py --check                                    # validation + Turon sizing only

# 2. Copy the job folder to the cluster and submit
scp -r jobs/CAI_QI_AS4_8552_30J cluster:/scratch/me/
ssh cluster 'cd /scratch/me/CAI_QI_AS4_8552_30J && sbatch submit.slurm'
```

The batch script runs:

```bash
abaqus job=CAI_QI_AS4_8552_30J input=CAI_QI_AS4_8552_30J.inp double=both cpus=16 \
       parallel=domain domains=16 mp_mode=mpi interactive
abaqus python postprocess.py CAI_QI_AS4_8552_30J.odb
```

There is no `user=`: Abaqus loads CompDam from `usub_lib_dir` in the job's `abaqus_v6.env`. `double=both` is mandatory for CompDam.

Other CLI options: `--energy`, `--material IM7_8552`, `--interface tie`, `--cpus`, `--output-dir`, `--min-element-size`, `--job-name`, `--config`. The generator also runs under `abaqus python cai_model.py ...`.

### Compiled library

Build the library once on the cluster, as described in CompDam's [readme](../../readme.md#building-a-shared-library):

```bash
cd CompDam_DGD/for && cp ../tests/abaqus_v6.env . && abaqus make library=CompDam_DGD
```

Point `COMPDAM["library_dir"]` (or `--library-dir`) at that directory.

To compile at run time instead, set `library_dir = None` and `user_subroutine = "/abs/path/CompDam_DGD.for"`. The env file then carries the Intel Fortran flags and the command gets `user=...`.

### Job folder contents

| File | Purpose |
|---|---|
| `<job>.inp` | complete flat Abaqus/Explicit deck (nodes, C3D8R elements, sets, orientations, CompDam material, contact, BCs, step, outputs) |
| `abaqus_v6.env` | `usub_lib_dir = '<library_dir>'` (plus `MPI_FLUSH_FCACHE`, from CompDam's sample env) |
| `submit.slurm` / `submit.pbs` | batch script built from `CLUSTER` (cpus, memory, walltime, partition, account, modules) |
| `run.sh` | the same commands for an interactive node |
| `postprocess.py` | standalone `abaqus python` script: load–displacement CSV, peak load, σ_CAI = P_max/(w·h), NCAMP delta, ALLKE/ALLIE check |
| `summary.txt` | laminate, CLT stiffness, damage seed, Turon sizing, per-block mesh statistics, NCAMP table, warnings |

## Fibre-aligned meshing

Each ply block gets its own structured quad mesh, extruded with `elements_per_ply` C3D8R layers per ply:

* **Damage zone (fibre-aligned).** Inside a parallelogram with two edges along the fibres and two along an axis, the grid is uniform. One edge family follows the fibre direction exactly and the other is the axis most nearly perpendicular to it. For 0° and 90° the zone is a rectangle; for ±45° the elements are rhombi with 45°/135° corners. Edges are `h` long in both directions.
* **Zone size per block.** The zone covers the seeded delaminations of the block's own interfaces, its matrix-crack disks and `zone_margin` (room for delamination growth). Plies near the impact face therefore get smaller zones than plies near the back face (the cone).
* **Outside the zone.** The rows and columns of the same structured grid are graded geometrically (`growth_rate`) up to `coarse_size`. The knife-edge lines are exact node lines. The end-clamp strips take every node with x ≤ `end_support_length`. Forcing that line through the slanted rows of a ±45° zone would create sliver elements, so it is not a node line.
* **Ply blocks.** `LAYUP["blocks"] = "auto"` groups consecutive plies with the same angle. The default `[45/0/-45/90]_3s` gives 23 blocks (the mid-plane 90/90 pair is one block) and 22 cohesive interfaces. Blocks share no nodes; general contact joins them, so their meshes need not match.
* **Quality checks.** Every block is checked for positive Jacobians and exact plate coverage. `summary.txt` lists element counts, the fibre-aligned element count, the minimum corner angle and the maximum aspect ratio. The highest aspect ratios (≈ coarse/fine) are in the coarse cross-bands outside the damage zone, which stay elastic.

## Delamination-zone element size (Turon et al. 2007)

`cai/turon.py` implements Turon et al.'s engineering solution for mesh-size effects in cohesive zone models (Eng. Fract. Mech. 74, 2007):

| Quantity | Equation |
|---|---|
| cohesive zone length | l_cz = M · E · G_c / τ₀², with M = 9π/32 |
| element size | l_e = l_cz / N_e, with N_e ≥ 3 |
| strength for a chosen l_e | τ₀' = √(9π · E · G_c / (32 · N_e · l_e)) |
| penalty stiffness | K = α · E₃ / t, with α = 50 |

l_cz is evaluated for mode I (E₃) and mode II (G₁₂ following Harper & Hallett 2008, or E₃ via `mode_II_modulus`). The smaller value governs. For the AS4/8552 defaults:

```
l_cz,I = 0.61 mm, l_cz,II = 0.42 mm  ->  Turon element size l_e = 0.14 mm (N_e = 3)
min_element_size = 0.30 mm           ->  strengths reduced 62.3/92.3 -> 51.4/62.9 MPa (3 elements in zone)
K = 50 * 9600 / 0.188 = 2.55e6 N/mm^3
```

The damage-zone element size is max(Turon size, `min_element_size`), capped at the CompDam snap-back limit. When the floor is active and `reduce_strength = True`, the strengths are lowered with Turon's formula: the fracture toughness is unchanged and N_e elements still span the cohesive zone. Set `min_element_size = None` to mesh at the full Turon size of 0.14 mm. That costs about 5× the elements of the 0.3 mm default. The same size is used for CompDam's intralaminar elements, well below their snap-back limit (0.82 mm).

**Model size.** The defaults (0.3 mm damage zone, 23 blocks) give about 1.26 M C3D8R elements and 2.5 M nodes. Generation takes ~30 s and the `.inp` is ~200 MB. For scoping runs, use `--min-element-size 0.5` (about 0.5 M elements). Note that a large K can lower the stable time increment; check the `.sta` file and reduce `penalty_alpha` towards 50 if needed.

## Rest of the model

| Aspect | Implementation |
|---|---|
| Geometry | 150 × 100 mm (ASTM D7137). X = loading, Z = thickness; ply 1 at z = 0 (back face), impact face at z = t |
| Layup | `[45/0/-45/90]_3s`, 0.188 mm plies, 4.51 mm |
| Elements | C3D8R, `*Section Controls`: distortion control (CompDam requirement), enhanced hourglass, average-strain kinematic split; one `*Orientation` per ply (1 = fibre, 3 = stacking direction) |
| CompDam material | `*User Material, constants=40` in CompDam order with blank fields for unused properties; named `*Depvar` (19 SDVs for flags `111101`); `*Characteristic Length, definition=USER, components=3`; non-default parameters in the `COMPDAM_PARAM` table collection |
| Impact damage | Conical, peanut-shaped delaminations (ellipse along the fibre direction of the ply below each interface) and matrix-crack disks (`CDM_d2 = 1`, `CDM_alpha = 0`) set via `*Initial Conditions, type=SOLUTION`; optional fibre damage. Radii scale with √((E − E_th)/(E_ref − E_th)) |
| Interfaces | General contact (`ALL EXTERIOR`). Bonded regions: cohesive behaviour (QUADS initiation, BK energy evolution) + friction. Seeded delaminations: friction only. Optional `*Tie` (`--interface tie`) |
| Fixture BCs | Loaded end kinematically coupled (U1) to `RP_LOAD`, which is driven by a smooth-step displacement; reacting end U1 = 0; end clamp strips U3 = 0 on both faces; knife edges U3 = 0 on both faces at y = inset and y = W − inset; U2 = 0 on one node line at x = 0 |
| Step | Explicit, NLGEOM, fixed mass scaling to a target Δt, bulk viscosity; outputs S, LE, SDV, U, RF, CSTRESS, CSDMG, energies, RP U1/RF1 |

## NCAMP validation workflow

`config.py` separates the values fed to CompDam (`MATERIALS`) from the NCAMP reference values (`NCAMP_BASELINE`: lamina RTD properties and the laminate CAI strength at 1500 in-lbf/in, i.e. 6.67 J/mm). `summary.txt` tabulates the model inputs against the NCAMP lamina values. The job's `postprocess.py` compares the predicted CAI strength with the NCAMP value.

Before you use the numbers quantitatively:

1. Enter the values from the NCAMP report revision you validate against (AS4/8552 tape: NCP-RP-2010-002 / CAM-RP-2010-002), then set `"verified": True`. The shipped values are representative and flagged `NOT verified`.
2. Calibrate `GXT` and `GXC` for AS4/8552. They are currently the IM7/8552 values from CompDam's `tests/IM7-8552-C3D.inp`.
3. Calibrate the damage seed (`radius_front`, `radius_back`, `aspect_ratio`) against your C-scans, and keep `zone_margin` large enough to hold the delamination growth before failure.
4. Run `--sweep` and compare CAI strength vs. impact energy with the NCAMP/test data.

## Limitations

* The fibre-aligned region is each block's damage zone. Outside it the grid is graded and no longer fibre-aligned, which is acceptable where the material stays elastic. CompDam's FKT models with random misalignment (φ0 = 1/2/3) assume a uniform fibre-aligned mesh, which therefore holds only inside the zones.
* Near-perpendicular delaminations in a ±45° block (e.g. a 90° ellipse) make the parallelogram large. If it would reach the supports, it is clipped and `summary.txt` lists a warning; check those blocks' margins.
* Not yet run in Abaqus from this environment (no Abaqus here). The deck is checked structurally by the tests: every referenced set and surface is defined, every element references existing nodes, every hex has a positive volume, and every keyword the model needs is present. Submit a short trial first (e.g. `STEP["time_period"]` 1e-4 and `--min-element-size 1.0`).

## Tests and documentation (uv)

```bash
uv run pytest                                                    # 38 tests, ~3 s
uv run --group docs sphinx-build -b html docs docs/_build/html   # API docs
```

Optional visual check in CAE: `abaqus cae noGUI=cae_import.py -- jobs/<job>/<job>.inp` writes `<job>.cae`.
