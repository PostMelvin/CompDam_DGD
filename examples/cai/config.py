"""
Configuration for the CompDam Compression-After-Impact (CAI) model.

Every modelling input (geometry, layup, material properties, mesh, cohesive
zone sizing, pre-defined impact damage, boundary conditions, solver, cluster
and job settings) lives in this module. The generator in :mod:`cai` reads
these settings and holds no hard-coded engineering values.

The generator is standalone: it needs only Python (2.7 or 3.x, no third-party
packages, no Abaqus/CAE) to write a complete Abaqus/Explicit input deck plus
the files to run it with a pre-compiled CompDam library on a cluster.

Units
-----
Consistent ``N - mm - s - tonne`` units, the system CompDam recommends:

* length: mm
* force: N
* stress / modulus: MPa
* fracture toughness: N/mm (kJ/m^2)
* density: tonne/mm^3
* energy: impact energies are given in J.

Coordinate system
-----------------
* Global X: loading (length) direction of the CAI coupon (0 deg fibre direction).
* Global Y: width direction.
* Global Z: thickness direction. Ply 1 is at z = 0 (back/non-impacted face),
  and the impacted face is at z = t.
* Ply angles are measured from global X, positive counter-clockwise about +Z.
"""

# =============================================================================
# Model / job identification
# =============================================================================
MODEL = {
    # Base name of the job; the material and impact energy are appended,
    # e.g. CAI_QI_AS4_8552_30J.
    "job_name": "CAI_QI",
    "description": "CompDam_DGD CAI model, ASTM D7137 configuration",
}

# =============================================================================
# Geometry (ASTM D7137/D7137M standard specimen)
# =============================================================================
GEOMETRY = {
    "length": 150.0,  # X, loading direction [mm]
    "width": 100.0,  # Y [mm]
    # Centre of the impact site (defaults to the plate centre when None).
    "impact_center": None,  # (x, y) [mm]
}

# =============================================================================
# Laminate definition
# =============================================================================
LAYUP = {
    # Key into MATERIALS below.
    "material": "AS4_8552",
    # Half-stack (if symmetric) listed from the BACK face (z = 0) upward.
    # [45/0/-45/90]_3s, 24 plies, ~4.5 mm thick (ASTM D7137 asks for 4-6 mm).
    "sequence": [45.0, 0.0, -45.0, 90.0] * 3,
    "symmetric": True,
    # Overrides the material's nominal cured ply thickness when not None [mm].
    "ply_thickness": None,
    # Grouping of plies into separately meshed "ply blocks". Delamination
    # interfaces (cohesive contact) exist only BETWEEN blocks.
    #   "auto" -> consecutive plies with the same angle form one block (e.g. the
    #             90/90 pair at the mid-plane); every angle change is an
    #             interface. Needed for fibre-aligned meshes, because each block
    #             is meshed along a single fibre direction.
    #   [n1, n2, ...] -> explicit ply counts per block (bottom to top). Blocks
    #             with mixed angles require MESH["fiber_aligned"] = False.
    "blocks": "auto",
}

# =============================================================================
# Material presets
# =============================================================================
# Each preset feeds the 40-entry CompDam property vector (see
# ``cai.compdam.PROPERTY_ORDER``). ``None`` is written as a blank field in the
# input deck, which disables the corresponding CompDam feature.
#
# IMPORTANT: The AS4/8552 set is assembled from NCAMP-style RTD lamina data
# plus literature values for the fracture properties. Check every value
# against the revision of the NCAMP report you are validating against
# (see NCAMP_BASELINE below) before drawing quantitative conclusions.
MATERIALS = {
    "AS4_8552": {
        "description": "Hexcel AS4/8552 unidirectional tape, 190 g/m^2 FAW",
        "density": 1.59e-9,  # tonne/mm^3
        "ply_thickness": 0.188,  # nominal cured ply thickness [mm]
        "elastic": {
            "E1": 135000.0,  # mean of 0 deg tension (~141 GPa) and compression (~128 GPa)
            "E2": 9600.0,
            "G12": 5100.0,
            "nu12": 0.32,
            "nu23": 0.487,  # literature (not an NCAMP-reported property)
            # Orthotropic extras: None -> transverse isotropy assumed by CompDam
            "E3": None,
            "G13": None,
            "G23": None,
            "nu13": None,
        },
        "strength": {
            "XT": 2207.0,  # 0 deg tension
            "XC": 1531.0,  # 0 deg compression
            "YT": 81.0,  # 90 deg tension
            "YC": 267.0,  # 90 deg compression
            "SL": 91.6,  # in-plane shear
        },
        "fracture": {
            "GYT": 0.28,  # mode I, matrix (8552 resin)
            "GSL": 0.79,  # mode II, matrix (8552 resin)
            "eta_BK": 1.634,  # Benzeggagh-Kenane exponent
            "GXT": 205.0,  # fibre tension; IM7/8552 value assumed -> calibrate
            "fXT": 0.5,  # 0.5/0.5 -> linear softening
            "fGXT": 0.5,
            "GXC": 61.0,  # fibre compression; IM7/8552 value assumed -> calibrate
            "fXC": 0.5,
            "fGXC": 0.5,
        },
        "other": {
            "alpha0": 0.925,  # fracture-plane angle in pure transverse compression [rad]
            "alpha_PL": 4.06e-9,  # Ramberg-Osgood shear nonlinearity (8552 resin, IM7/8552 fit)
            "n_PL": 5.4,
            "mu": 0.3,  # friction on matrix-crack surfaces
            "alpha11": None,
            "alpha22": None,
            "cl": None,
            "w_kb": None,  # kink-band width, only needed for FKT (flag pos. 4 = 3/4/5)
            "T_sf": None,
        },
        # Interlaminar properties for the surface-based cohesive contact.
        # Strengths may be lowered automatically (Turon et al. 2007) when the
        # mesh is coarser than the cohesive-zone requirement, see COHESIVE_ZONE.
        "interface": {
            "strength_normal": 62.3,  # MPa
            "strength_shear": 92.3,  # MPa
            "GIc": 0.28,  # N/mm
            "GIIc": 0.79,  # N/mm
            "eta_BK": 1.634,
            "friction": 0.3,
        },
    },
    # IM7/8552 set taken verbatim from CompDam's tests/IM7-8552-C3D.inp (refs:
    # Camanho et al. 2007, 2013; Catalanotti et al. 2014). 134 g/m^2 tape.
    "IM7_8552": {
        "description": "Hexcel IM7/8552 unidirectional tape, 134 g/m^2 FAW",
        "density": 1.57e-9,
        "ply_thickness": 0.131,
        "elastic": {
            "E1": 171420.0,
            "E2": 9080.0,
            "G12": 5290.0,
            "nu12": 0.32,
            "nu23": 0.52,
            "E3": None,
            "G13": None,
            "G23": None,
            "nu13": None,
        },
        "strength": {
            "XT": 2326.2,
            "XC": 1200.1,
            "YT": 62.3,
            "YC": 199.8,
            "SL": 92.30,
        },
        "fracture": {
            "GYT": 0.277,
            "GSL": 0.788,
            "eta_BK": 1.634,
            "GXT": 205.0,
            "fXT": 0.63,
            "fGXT": 0.469,
            "GXC": 61.0,
            "fXC": 0.614,
            "fGXC": 0.432,
        },
        "other": {
            "alpha0": 0.925,
            "alpha_PL": 4.06e-9,
            "n_PL": 5.4,
            "mu": 0.3,
            "alpha11": None,
            "alpha22": None,
            "cl": None,
            "w_kb": None,
            "T_sf": None,
        },
        "interface": {
            "strength_normal": 62.3,
            "strength_shear": 92.3,
            "GIc": 0.277,
            "GIIc": 0.788,
            "eta_BK": 1.634,
            "friction": 0.3,
        },
    },
}

# =============================================================================
# NCAMP baseline data for validation / comparison
# =============================================================================
# Reference values for comparison with the model inputs and outputs. The
# build summary lists the deltas between MATERIALS and these entries, and the
# generated postprocess.py compares the predicted CAI strength with "cai".
#
# "verified": False means the numbers are representative values that have
# NOT been checked line by line against the cited report. Replace them with
# the values (mean, B-basis, normalised or as-measured) from the report
# revision you validate against, then set "verified": True.
NCAMP_BASELINE = {
    "AS4_8552": {
        "report": "NCAMP NCP-RP-2010-002 / CAM-RP-2010-002 (Hexcel 8552 AS4 unitape)",
        "condition": "RTD",
        "basis": "mean",
        "verified": False,
        "lamina": {
            "E1t": 141000.0,  # MPa (0 deg tension modulus)
            "E1c": 128000.0,  # MPa (0 deg compression modulus)
            "E2t": 9600.0,
            "G12": 5100.0,
            "nu12": 0.32,
            "XT": 2207.0,  # F1tu
            "XC": 1531.0,  # F1cu
            "YT": 81.0,  # F2tu
            "YC": 267.0,  # F2cu
            "SL": 91.6,  # F12 (5 % strain)
            "ply_thickness": 0.188,
        },
        # ASTM D7136/D7137 CAI on a quasi-isotropic laminate,
        # 1500 in-lbf/in (6.67 J/mm) impact.
        "cai": {
            "impact_energy_per_thickness": 6.67,  # J/mm
            "cai_strength": 176.3,  # MPa (25.57 ksi); verify against the report
            "dent_depth": None,  # mm
            "damage_area": None,  # mm^2 (C-scan)
        },
    },
    "IM7_8552": {
        "report": "NCAMP NCP-RP-2009-028 (Hexcel 8552 IM7 unitape)",
        "condition": "RTD",
        "basis": "mean",
        "verified": False,
        "lamina": {},
        "cai": {
            "impact_energy_per_thickness": 6.67,
            "cai_strength": None,
            "dent_depth": None,
            "damage_area": None,
        },
    },
}

# =============================================================================
# CompDam_DGD VUMAT settings
# =============================================================================
COMPDAM = {
    # Directory on the CLUSTER that holds the pre-compiled CompDam shared
    # library (built with "abaqus make library=CompDam_DGD"). It is written to
    # usub_lib_dir in each job's abaqus_v6.env, so the job runs without
    # "user=". Use an absolute path; the generator does not check that it
    # exists, since it may run on another machine.
    "library_dir": "/path/to/compdam/lib",
    # Alternative to library_dir: compile at run time from source. Set to the
    # absolute path of CompDam_DGD.for on the cluster (and library_dir to
    # None); the generated env file then includes the Intel Fortran flags.
    "user_subroutine": None,
    # Feature flags (positions):
    #   1 matrix damage        (1 = intralaminar DGD cracking)
    #   2 shear nonlinearity   (0 off, 1 RO 1-2, 2 Schapery, 3 RO 3-D, 4 RO 1-3, 5 Schaefer)
    #   3 fibre tension        (1 = max strain, bilinear softening)
    #   4 fibre compression    (1 = max strain BL, 3/4/5 = FKT 12/13/3-D)
    #   5 energy output        (0 = all mechanisms)
    #   6 friction             (1 = on)
    "feature_flags": "111101",
    # Element deletion based on CDM_STATUS (SDV 11). CompDam recommends
    # leaving this off; see its readme, "Element deletion".
    "element_deletion": False,
    # Number of components for *Characteristic Length (3 for C3D8R;
    # 6 when a random fibre misalignment phi0 = 1/2/3 is used with FKT).
    "characteristic_length_components": 3,
    # Initial fibre misalignment phi0 for FKT (rad). Ignored unless FKT is on.
    "phi0": 0.0,
    # Non-default model parameters written to the COMPDAM_PARAM table
    # collection. Leave empty to use the CompDam defaults.
    "parameters_float": {
        # "gamma_max": 4.0,
    },
    "parameters_int": {
        "logLevel": 2,
        # Delete (rather than abort on) the rare element where DGD fails to
        # converge, so one bad element does not end a long CAI run.
        "terminate_on_no_convergence": 0,
    },
}

# =============================================================================
# Mesh
# =============================================================================
MESH = {
    "element_type": "C3D8R",
    # Fibre-aligned in-plane meshing: inside each ply's damage zone the element
    # edges follow the fibre direction (parallelogram elements for +/-45 deg,
    # rectangles for 0/90 deg). Every block gets its own mesh; blocks are
    # joined through general contact, so their meshes need not match.
    # False -> orthogonal mesh in every block (also allows mixed-angle blocks).
    "fiber_aligned": True,
    # The damage-zone element size comes from COHESIVE_ZONE (Turon et al.).
    # Outside the damage zone the mesh grades geometrically up to coarse_size.
    "coarse_size": 4.0,  # mm
    "growth_rate": 1.2,  # max size ratio of neighbouring elements in the grading
    # The fine zone of each block covers the seeded delaminations at the
    # block's interfaces and its matrix cracks, plus this margin for
    # delamination growth during compression [mm].
    "zone_margin": 8.0,
    # Minimum half-size of the fine zone around the impact site [mm]
    # (also used when the impact energy is below the damage threshold).
    "min_zone_half_size": 10.0,
    # Elements through the thickness of each ply (1 is the CompDam practice).
    "elements_per_ply": 1,
    # Element controls (*Section Controls)
    "hourglass_control": "ENHANCED",  # ENHANCED | RELAX_STIFFNESS | STIFFNESS | DEFAULT
    "distortion_control": True,  # required by CompDam ("distortion control=YES")
    "kinematic_split": "AVERAGE_STRAIN",
    # Warn when the total element count exceeds this value.
    "max_elements_warning": 1500000,
    # Geometric tolerance used when selecting nodes by coordinates [mm].
    "tolerance": 1.0e-6,
}

# =============================================================================
# Cohesive-zone driven element size in the delamination zone (Turon et al.)
# =============================================================================
# Turon A., Davila C.G., Camanho P.P., Costa J. "An engineering solution for
# mesh size effects in the simulation of delamination using cohesive zone
# models", Eng Fract Mech 74 (2007) 1665-1682.
#
#   cohesive zone length   l_cz = M * E * Gc / tau0^2,    M = 9*pi/32
#   element size           l_e  = l_cz / N_e              (N_e >= 3 elements)
#   strength reduction     tau0' = sqrt(9*pi*E*Gc / (32 * N_e * l_e))
#   penalty stiffness      K = alpha * E3 / t,            alpha >= 50
#
# l_cz is evaluated for mode I (E = E3) and mode II (E per "mode_II_modulus")
# and the smaller value governs. The resulting l_e is often ~0.1-0.2 mm, which
# is expensive for a 24-ply CAI plate. Set "min_element_size" to the smallest
# size you can afford: if it is larger than the Turon size, the interface
# strengths are reduced with Turon's formula (when "reduce_strength" is True)
# so that N_e elements still span the cohesive zone, and the fracture
# toughness is unchanged.
COHESIVE_ZONE = {
    "M": 9.0 * 3.141592653589793 / 32.0,  # Rice / Falk et al. parameter used by Turon
    "elements_in_zone": 3,  # N_e; Turon et al. recommend at least 3
    # Modulus for the mode II cohesive-zone length: "G12" (Harper & Hallett
    # 2008, more conservative) or "E3" (Turon et al. 2007).
    "mode_II_modulus": "G12",
    "min_element_size": 0.3,  # mm; floor on the damage-zone element size
    "max_element_size": None,  # mm; optional ceiling (defaults to the CompDam snap-back limit)
    "reduce_strength": True,
    # Interface penalty stiffness: "turon" -> alpha*E3/t_ply, or a number [N/mm^3].
    "penalty": "turon",
    "penalty_alpha": 50.0,
}

# =============================================================================
# Impact energies and pre-defined impact damage state
# =============================================================================
# The CAI model starts from an idealised damage state that represents the
# impact (no impact step is simulated). Damage footprints are scaled from a
# reference energy with an empirical rule: projected damage area grows
# linearly with impact energy above a threshold, so radii scale with
# sqrt((E - E_th) / (E_ref - E_th)). Calibrate against C-scan data when
# available.
IMPACT = {
    # Energies [J] to generate. cai_model.py builds one job per value with
    # --sweep, or uses "default_energy" otherwise.
    "energies": [15.0, 22.5, 30.0, 40.0],
    "default_energy": 30.0,  # ~6.67 J/mm * 4.5 mm (ASTM D7136 standard)
    "threshold_energy": 5.0,  # J, below which no damage is seeded
    "reference_energy": 30.0,  # J, energy the reference radii below belong to
}

DAMAGE = {
    "enabled": True,
    # --- Delamination (interlaminar) ------------------------------------
    # Elliptical footprint per interface. The equivalent radius varies
    # linearly from the impacted face (front) to the back face, giving the
    # classic conical ("pine tree") distribution. Each ellipse's major axis
    # follows the fibre direction of the ply below the interface (the
    # "peanut" shape seen in C-scans).
    "delamination": {
        "enabled": True,
        "radius_front": 5.0,  # mm, at reference energy
        "radius_back": 20.0,  # mm, at reference energy
        "aspect_ratio": 1.6,  # major / minor axis (area preserving)
    },
    # --- Matrix cracking (intralaminar) ---------------------------------
    # Elements whose centroid lies within the (circular) radius are seeded
    # with a fully formed matrix crack: CDM_d2 = 1, CDM_alpha = "alpha".
    "matrix_cracking": {
        "enabled": True,
        "radius_front": 3.0,
        "radius_back": 15.0,
        "d2": 1.0,
        "alpha": 0,  # crack normal angle [deg, integer]; 0 -> through-thickness crack
    },
    # --- Fibre damage under the impactor --------------------------------
    # Off by default. Seeding CDM_d1T / CDM_d1C changes the initial fibre
    # stiffness; check the CompDam documentation on how the damage
    # variables evolve before turning this on.
    "fiber_damage": {
        "enabled": False,
        "radius": 3.0,  # mm (not scaled with energy)
        "plies_from_front": 4,  # number of plies from the impacted face
        "d1T": 0.0,
        "d1C": 0.9,
    },
}

# =============================================================================
# Boundary conditions (ASTM D7137 support fixture)
# =============================================================================
BOUNDARY = {
    # Length of the top/bottom clamping plates along X. Out-of-plane motion
    # (Uz) of both faces is restrained over this strip at each loaded end.
    "end_support_length": 8.0,  # mm
    # Side knife edges: lines of Uz = 0 on both faces at y = inset and
    # y = width - inset, spanning between the end supports.
    "side_support_inset": 5.0,  # mm (0 -> on the free edges)
    # End shortening applied through a reference point coupled to the
    # loaded end (x = length) [mm].
    "applied_displacement": 1.5,
}

# =============================================================================
# Abaqus/Explicit step
# =============================================================================
STEP = {
    "name": "Compression",
    # Long enough to stay quasi-static: check ALLKE << ALLIE afterwards.
    "time_period": 0.02,  # s
    # Fixed mass scaling to a target stable increment [s], applied to
    # elements below it (None -> no mass scaling).
    "mass_scaling_target_dt": 2.0e-7,
    "linear_bulk_viscosity": 0.06,
    "quadratic_bulk_viscosity": 1.2,
    "field_output_intervals": 50,
    "history_output_intervals": 500,
}

# =============================================================================
# Interfaces between ply blocks
# =============================================================================
INTERFACE = {
    # "cohesive": general contact with surface-based cohesive behaviour on the
    #             bonded regions and frictional contact inside the seeded
    #             delaminations (delamination growth).
    # "tie":      perfectly bonded blocks (useful for debugging/baselines).
    "type": "cohesive",
}

# =============================================================================
# Output location
# =============================================================================
JOB = {
    # One sub-folder per job is created here, holding everything needed to
    # run it: <job>.inp, abaqus_v6.env, the scheduler script and
    # postprocess.py. Copy the folder to the cluster as-is.
    "output_dir": "jobs",
}

# =============================================================================
# Cluster execution
# =============================================================================
CLUSTER = {
    "scheduler": "slurm",  # "slurm" | "pbs" | "none" (plain shell script only)
    "abaqus_command": "abaqus",
    # Shell lines run before Abaqus (e.g. module loads); one string per line.
    "setup_commands": ["module load abaqus"],
    "cpus": 16,  # also the number of Explicit domains
    "nodes": 1,
    "memory": "64G",
    "walltime": "48:00:00",
    "partition": None,
    "account": None,
    "mp_mode": "mpi",
    # Extra scheduler directives, e.g. ["--qos=long"] (SLURM) or ["-q long"] (PBS).
    "extra_directives": [],
    # Run postprocess.py with "abaqus python" after the analysis finishes.
    "postprocess": True,
}
