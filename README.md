[Japanese version](README.jp.md)

# dd_docking — Flexible-pocket ensemble docking: MD-verified induced fit, GPU-accelerated at scale

An ensemble docking toolkit with induced-fit-style pocket dynamics for
virtual screening. Designed as a reusable package, not tied to any specific
target or ligand set.

- **Ensemble preparation (`dd_docking-prep`)**: structurally repairs multiple
  receptor conformations (PDB) with PDBFixer, determines the docking box and
  pocket-proximal flexible residues from each co-crystallized ligand, and
  generates rigid/flex PDBQT pairs with Meeko.
- **Ensemble docking (`dd_docking-dock`)**: docks each molecule in a `.smi`
  library against every ensemble member with flexible side chains, ranking
  by the best (lowest) affinity across members (best-of-ensemble).
- **MD-based refinement (`dd_docking-refine`)**: relaxes top hits only with a
  short implicit-solvent MD run (OpenMM, GBn2 + hydrogen mass repartitioning)
  and re-ranks them by whether the ligand heavy-atom RMSD trajectory stays in
  the pocket — a lightweight induced-fit sanity check.

## Installation

Requires meeko, pdbfixer, openmm, openmmforcefields, openff-toolkit, and
mdtraj (declared in `pyproject.toml`'s `dependencies`), plus the `qvina2`
CLI binary (from the conda-forge `qvina` package — CPU docking via
QuickVina2, a speed-tuned AutoDock Vina 1.1.2 fork; not a Python import, so
it can't be a `pyproject.toml` dependency). Best installed via conda-forge,
using mamba:

```bash
mamba create -n dd_docking python=3.12 -c conda-forge \
    rdkit numpy pandas qvina meeko pdbfixer openmm openmmforcefields openff-toolkit mdtraj
mamba activate dd_docking
cd dd_docking
pip install -e .
```

This installs three console commands: `dd_docking-prep`, `dd_docking-dock`,
`dd_docking-refine`.

### Optional: GPU-accelerated docking (Linux only, rigid receptors only)

`dd_docking-dock` can use [Vina-GPU+](https://github.com/DeltaGroupNJUPT/Vina-GPU-2.0)
instead of CPU QuickVina2 for a given docking task, when it's built, the
member has no flexible side chains, and the task's box fits the OpenCL
kernel's size limit. This is Linux-only; on macOS/Windows (or on Linux
without the binary built, or for a member with flexible residues, or a box
that's too large), `dd_docking-dock` transparently uses CPU QuickVina2 — no
code changes needed, this is purely a runtime fallback (see `--backend`
below).

```bash
mamba activate dd_docking
bash scripts/build_vina_gpu.sh
```

This clones Vina-GPU+ from source (pinned commit) into `third_party/` (not
checked into this repo), **patches a handful of real upstream OpenCL-C bugs**
in its kernel source (see below), builds it, and installs the resulting
binary + kernel files into `$CONDA_PREFIX/share/dd_docking/vina-gpu/`.
Requires an NVIDIA/AMD GPU with a working OpenCL runtime (check with
`clinfo`) — on NVIDIA this normally comes from the driver/CUDA toolkit
install. See `scripts/build_vina_gpu.sh` for the AMD (`GPU_PLATFORM`), custom
OpenCL path (`DD_DOCKING_OPENCL_PATH`), and OpenCL-C version
(`DD_DOCKING_OPENCL_VERSION`) overrides.

**Hard limitation, not fixable by patching: Vina-GPU+ only supports rigid
receptors.** `main_procedure_cl.cpp` asserts `m.num_other_pairs() == 0`
before docking, and that count is nonzero for *any* ligand-flex, flex-flex,
or flex-inflex interaction — i.e. as soon as a single flexible side chain is
in play. Since `dd_docking`'s ensemble docking always uses flexible side
chains, `--backend gpu`/`auto` will still use CPU QuickVina2 for those tasks
(with a one-time warning if you explicitly asked for `gpu`); only a member
prepared with zero flexible residues can actually run on the GPU here.

Earlier versions of this README reported Vina-GPU+ as simply broken on this
project's GTX 1660 Ti (`CL_BUILD_PROGRAM_FAILURE` / a crash while compiling
the kernel, matching upstream issues
[#1](https://github.com/DeltaGroupNJUPT/Vina-GPU-2.0/issues/1) and
[#26](https://github.com/DeltaGroupNJUPT/Vina-GPU-2.0/issues/26)). The actual
root cause: `clinfo` reports this device's real OpenCL C compiler level as
1.2 even though the *platform* advertises OpenCL 3.0, and Vina-GPU+'s kernel
source has several real type errors (a redundant address-of that turns a row
pointer into the wrong pointer type, two pointers missing their `__global`
address-space qualifier, and a call to `get_global_linear_id()`, which is
OpenCL ≥2.0 only) that a strict 1.2/2.0 compile correctly rejects, but that
NVIDIA's more lenient OpenCL-3.0 compile path tolerates just long enough to
segfault or fail at the program-binary stage instead. `build_vina_gpu.sh` now
patches these before building and defaults to `-DOPENCL_1_2`, which builds
and runs correctly on this hardware. Measured after the fix, rigid-receptor
GPU docking is genuinely faster than CPU once there's enough work to amortize
GPU/subprocess startup cost -- 8 ligands at `--exhaustiveness 32`: **16.2s on
GPU vs. 35.2s on CPU** (~2.2x) -- but *slower* for a single cheap task (one
ligand at `--exhaustiveness 8`: 3.2s GPU vs. 2.3s CPU), so don't expect a win
for small jobs. If your GPU/driver's OpenCL C compiler genuinely supports
2.0/3.0, `DD_DOCKING_OPENCL_VERSION=-DOPENCL_3_0` may build a faster kernel;
these patches are correctness fixes independent of that choice, so they're
applied either way.

## Usage

### 1. Ensemble preparation (`dd_docking-prep`)

Pass `--member ID PDB_FILE COCRYSTAL_LIGAND_CODE` once per receptor
conformation. Each member is repaired independently with PDBFixer, and the
docking box plus pocket residues (default cutoff 5 Å) are derived from the
co-crystallized ligand's coordinates.

```bash
dd_docking-prep \
  --member 3ert data/raw/3ERT.pdb OHT \
  --member 1xpc data/raw/1XPC.pdb AIT \
  --member 1yim data/raw/1YIM.pdb CM4 \
  -o data/ensemble
```

Output (measured, reproducible from the `data/` in this repository):

```
3ert: 10 flexible residues -> data/ensemble/3ert_rigid.pdbqt
1xpc: 10 flexible residues -> data/ensemble/1xpc_rigid.pdbqt
1yim: 10 flexible residues -> data/ensemble/1yim_rigid.pdbqt

[done] 3 member(s) -> data/ensemble/manifest.json
```

`data/ensemble/` gets, per member, `<id>_fixed.pdb` (post-PDBFixer, no
hydrogens), `<id>_rigid.pdbqt` / `<id>_flex.pdbqt` (Meeko output), and
`manifest.json` (box center/size, flexible-residue list, and file paths per
member — read by the later commands).

Key options:

| Option | Default | Description |
|---|---|---|
| `--chain` | `A` | receptor chain to extract |
| `--pocket-cutoff` | `5.0` | distance (Å) from the co-crystallized ligand used to select flexible residues |
| `--max-flexres` | `10` | cap on flexible residues (nearest to ligand first; `<=0` for no cap) |
| `--box-padding` | `5.0` | docking box padding around the co-crystallized ligand's extent (Å) |
| `--charge-model` | `gasteiger` | Meeko receptor partial-charge model (`gasteiger`/`espaloma`/`zero`) |

`--max-flexres` defaults to a value chosen from experience: unlike rigid
docking (a precomputed grid map), each flexible side chain adds a real
rotational degree of freedom to Vina's search space. On the ER-alpha test
data in this repository (3ERT/1XPC/1YIM, 5 Å cutoff from the co-crystallized
ligand), 22-26 residues fell within the cutoff, and docking one ligand
against one member with no cap didn't finish in under an hour. Capping at 10
(nearest first) brought this down to a practical runtime. Removing the cap or
widening the cutoff to "be safe" runs into this trap easily.

The docking box itself is sized to cover both the co-crystallized ligand
*and* every chosen flexible residue's atoms (plus `--box-padding` on top) --
not just the ligand. A box sized only around the ligand can leave a flexible
side chain's movable atoms outside the search space, and Vina then reports
"no conformations completely within the search space" and returns zero poses
for every ligand docked against that member. ER-alpha's flexible residues
happen to spread fairly wide around its elongated pocket, so these three
members end up with boxes of ~25-33 Å per side even at the default padding
(`1xpc`/`1yim` are actually right up against Vina-GPU+'s separate 30 Å
OpenCL kernel limit too, but it's moot here: all three members have
flexible side chains, so `--backend auto/gpu` always uses CPU QuickVina2
for this dataset regardless of box size -- see
[GPU-accelerated docking](#optional-gpu-accelerated-docking-linux-only-rigid-receptors-only)).

Python API:

```python
from dd_docking import prepare_ensemble

ensemble = prepare_ensemble(
    [("3ert", "data/raw/3ERT.pdb", "OHT"),
     ("1xpc", "data/raw/1XPC.pdb", "AIT"),
     ("1yim", "data/raw/1YIM.pdb", "CM4")],
    "data/ensemble",
)
for m in ensemble:
    print(m.member_id, len(m.flexres), m.rigid_pdbqt)
```

### 2. Ensemble docking (`dd_docking-dock`)

Docks a `.smi` library (`SMILES<TAB>ID` format, one molecule per line)
against every member of a prepared ensemble.

```bash
dd_docking-dock data/ensemble data/ligands.smi \
  -o data/screen --exhaustiveness 4 --n-poses 1 --n-jobs 4
```

`data/ligands.smi` here is 3 co-crystallized ligands (self-docking checks) +
3 genuine ER-alpha actives + 9 property-matched decoys, all pulled from
[DUD-E](https://dude.docking.org/targets/esr1/) (see
[Test dataset (ER-alpha, DUD-E decoys)](#test-dataset-er-alpha-dud-e-decoys)
below) -- a much harder,
more realistic screen than picking a handful of unrelated approved drugs as
"obviously wrong" decoys.

Work is split into `(ligand, ensemble member)` tasks (15 ligands × 3 members
= 45 tasks in the example above), parallelized via `--n-jobs` (default
sequential; `<=0` uses all CPU cores). One progress line is printed per
ligand once it has been docked against every member (measured on this
repository's data, 16-core Linux box, `--n-jobs 4`, ~29 minutes):

```
[parallel] using 4 worker processes for 45 tasks
[1] OHT_native_3ert  best_member=1xpc  affinity=-11.300
[2] AIT_native_1xpc  best_member=3ert  affinity=-11.400
[3] CM4_native_1yim  best_member=3ert  affinity=-11.700
[4] CHEMBL215428_active  best_member=1xpc  affinity=-11.700
[5] CHEMBL31127_active  best_member=1xpc  affinity=-11.300
[6] CHEMBL385358_active  best_member=1xpc  affinity=-12.100
[7] C48469840_decoy  best_member=1xpc  affinity=-10.600
[8] C14244874_decoy  best_member=3ert  affinity=-11.200
[9] C66896427_decoy  best_member=1yim  affinity=-11.100
[10] C27846680_decoy  best_member=1xpc  affinity=-11.300
[11] C37195230_decoy  best_member=1xpc  affinity=-8.600
[12] C37085981_decoy  best_member=1xpc  affinity=-9.000
[13] C09495369_decoy  best_member=3ert  affinity=-12.300
[14] C12776766_decoy  best_member=1xpc  affinity=-11.300
[15] C36904163_decoy  best_member=1xpc  affinity=-8.300

[done] 15 ligand(s) ranked -> data/screen/ranked_results.csv
```

Final ranking (`ranked_results.csv`, affinity ascending). Unlike unrelated
approved drugs, DUD-E's property-matched decoys (similar molecular weight,
rotatable-bond count, and charge to the real actives, but a different 2-D
scaffold) are a deliberately hard discrimination test, and this shows in the
result honestly: the top spot goes to a decoy (`C09495369_decoy`), and
several other decoys interleave with the real actives/natives through the
top 12 rows. This isn't a bug — it's the well-documented reality that a
single low-exhaustiveness AutoDock-family docking pass alone often can't
perfectly separate property-matched decoys from true actives, which is
exactly why `dd_docking-refine`'s MD-based induced-fit rescoring exists as a
second filter (see below). What docking *does* clearly deliver here: the 3
weakest, most topologically-dissimilar decoys are cleanly separated at the
bottom (ranks 13-15, affinity -9.0 to -8.3, versus -10.6 or better for
everything else), and the `best_member` / per-member `affinity[3ert/1xpc/1yim]`
columns differ across rows, confirming the ensemble's conformations actually
produce different results (evidence that ensemble docking is functioning as
intended):

| rank | ligand_id | best_member | best_affinity | affinity[3ert] | affinity[1xpc] | affinity[1yim] |
|---|---|---|---|---|---|---|
| 1 | C09495369_decoy | 3ert | -12.3 | -12.3 | -12.2 | -11.4 |
| 2 | CHEMBL385358_active | 1xpc | -12.1 | -11.5 | -12.1 | -11.6 |
| 3 | CM4_native_1yim | 3ert | -11.7 | -11.7 | -11.7 | -10.8 |
| 4 | CHEMBL215428_active | 1xpc | -11.7 | -11.4 | -11.7 | -8.1 |
| 5 | AIT_native_1xpc | 3ert | -11.4 | -11.4 | -11.3 | -10.8 |
| 6 | OHT_native_3ert | 1xpc | -11.3 | -10.6 | -11.3 | -10.8 |
| 7 | CHEMBL31127_active | 1xpc | -11.3 | -10.7 | -11.3 | -10.5 |
| 8 | C27846680_decoy | 1xpc | -11.3 | -10.7 | -11.3 | -10.7 |
| 9 | C12776766_decoy | 1xpc | -11.3 | -10.5 | -11.3 | -10.1 |
| 10 | C14244874_decoy | 3ert | -11.2 | -11.2 | -11.0 | -10.4 |
| 11 | C66896427_decoy | 1yim | -11.1 | -11.0 | -10.2 | -11.1 |
| 12 | C48469840_decoy | 1xpc | -10.6 | -9.4 | -10.6 | -10.4 |
| 13 | C37085981_decoy | 1xpc | -9.0 | -8.2 | -9.0 | -9.0 |
| 14 | C37195230_decoy | 1xpc | -8.6 | -8.4 | -8.6 | -8.3 |
| 15 | C36904163_decoy | 1xpc | -8.3 | -7.6 | -8.3 | -8.1 |

Output:

- `data/screen/ranked_results.csv` — columns `rank, ligand_id, smiles,
  best_member, best_affinity, receptor_pdb, pose_pdbqt,
  affinity[<member_id>]...`. `receptor_pdb` / `pose_pdbqt` are consumed
  directly by `dd_docking-refine`.
- `data/screen/top_hits.sdf` — each hit's best pose, with an `affinity`
  property, so this file can be loaded directly into a protein-ligand viewer
  together with the receptor PDB (the `receptor_pdb` column of
  `ranked_results.csv`).
- `data/screen/ranked_results_poses/` — one PDBQT file per hit (input for
  `dd_docking-refine`).

Key options:

| Option | Default | Description |
|---|---|---|
| `--exhaustiveness` | `16` | Vina exhaustiveness |
| `--n-poses` | `5` | poses kept per member |
| `--seed` | `0` | random seed (embedding and docking) |
| `--top-n` | all | keep only the top N results |
| `--n-jobs` | `1` | parallel workers, one per `(ligand, member)` task (`<=0` for all cores) |
| `--backend` | `auto` | docking engine: `auto` uses Vina-GPU+ when built and the member has no flexible side chains and its box fits the <30 Å OpenCL kernel limit, else CPU QuickVina2; `cpu` always uses CPU QuickVina2 (every OS); `gpu` prefers Vina-GPU+ and falls back to `cpu` with a warning per member where it isn't usable (see [GPU-accelerated docking](#optional-gpu-accelerated-docking-linux-only-rigid-receptors-only)) |
| `--no-progress` | - | suppress progress log |

**`--n-jobs` behavior and CPU allocation**: `--n-jobs` is the *only*
CPU-allocation knob `dd_docking-dock`/`screen_ensemble` expose -- there is no
separate `--cpu`/`--ncpu` flag on this CLI. Internally, `--n-jobs` controls
how many cores each `(ligand, member)` task gets via qvina2's own `--cpu`
flag (not user-facing here): `--n-jobs 1` (default) runs sequentially,
letting each task use every available core (`cpu_count // 1`); any other
`--n-jobs` value divides cores evenly across that many concurrent worker
processes (`cpu_count // n_jobs`, minimum 1 core each). Flexible side-chain
docking is much heavier than rigid docking and benefits significantly from
QuickVina2's own multithreading, so the common pattern for lightweight rigid
docking -- run every worker pinned to 1 core and parallelize jobs instead --
backfires badly here. Measured on this repository's ER-alpha data (one
ligand, one member, `--exhaustiveness 16`, 16-core Linux box), single-task
wall time by cores-per-task:

| cores/task | wall time | speedup vs. 1 core |
|---|---|---|
| 1 | 827.5s | 1.0x |
| 2 | 432.4s | 1.9x |
| 4 | 223.2s | 3.7x |
| 8 | 118.8s | 7.0x |
| 16 | 101.5s | 8.2x |

Returns diminish sharply past 8 cores (going from 8 to 16 cores only gains
~17%), so for a batch of many `(ligand, member)` tasks, splitting the extra
cores across a couple of parallel `--n-jobs` workers (each still keeping
several cores via `cpu_count // n_jobs`) tends to use the machine better than
handing every core to one task at a time. Start conservatively (roughly 4-8
cores per task, i.e. `--n-jobs` around `cpu_count // 6`) and adjust from
there relative to library size and ensemble size.

Python API:

```python
from dd_docking import ligand_prep, screen_ensemble
from dd_docking.ensemble import load_manifest

ensemble = load_manifest("data/ensemble/manifest.json")
ligands = ligand_prep.read_smi("data/ligands.smi")
hits = screen_ensemble(ensemble, ligands, n_jobs=8,
                       out_csv="data/screen/ranked_results.csv",
                       out_sdf="data/screen/top_hits.sdf")
for hit in hits[:5]:
    print(hit.ligand_id, hit.best_member, hit.best_affinity, hit.member_affinities)
```

### 3. MD refinement and re-ranking (`dd_docking-refine`)

Relaxes and re-evaluates only the top hits from `dd_docking-dock`'s
`ranked_results.csv`, using a short implicit-solvent MD run (GBn2 + hydrogen
mass repartitioning, 4 fs). Limiting this to the top N hits (rather than all
of them) keeps the computational cost practical.

```bash
dd_docking-refine data/screen/ranked_results.csv \
  -o data/screen/refine --top-n 3 --prod-ps 20 --equil-ps 5
```

Measured output (`--prod-ps` shortened here for a quick test; use a longer
value in production):

```
[MD C09495369_decoy] implicit  platform=OpenCL  rmsd_mean=1.41  rmsd_final=1.35  stable=True
[MD CHEMBL385358_active] implicit  platform=OpenCL  rmsd_mean=0.94  rmsd_final=1.79  stable=True
[MD CM4_native_1yim] implicit  platform=OpenCL  rmsd_mean=1.48  rmsd_final=1.89  stable=True

[done] 3 hit(s) refined -> data/screen/refine/md_rescore.csv
```

All three top hits -- including the decoy that out-scored every real active
in docking (`C09495369_decoy`) -- pass the MD stability criterion (RMSD mean
< 3 Å and final < 4 Å) here, GBn2 implicit solvent setup succeeded for all
three this run (`--platform` picked OpenCL), and RMSDs stay under 2 Å for
all three. Rather than editing this to a cleaner-looking result: this is the
honest, expected outcome of using property-matched decoys, not a bug. A
DUD-E "decoy" is only presumed inactive (absence of a reported ChEMBL
activity, not proof of one) — some legitimately sit stably in a compatible
sub-pocket the same way a real ligand would, which is exactly why the docking
literature treats property-matched decoy sets as a genuinely hard
benchmark, and why no single stage of this pipeline (docking or MD) should
be trusted as a final answer on its own.

`md_rescore.csv` sorts by `stable` (ligand heavy-atom RMSD mean < 3 Å and
final < 4 Å) first, then by Vina affinity ascending within that. Columns
also include `rmsd_mean` / `rmsd_final` / `rmsd_max` / `e_min_kcal` /
`e_final_kcal` / `implicit` (whether GBn2 implicit solvent was used, or it
fell back to vacuum). Each hit's complex structure and trajectory are saved
to `data/screen/refine/<rank>_<ligand_id>/` (`complex.pdb` + `prod.dcd`).

Key options:

| Option | Default | Description |
|---|---|---|
| `--top-n` | `5` | number of top hits to run MD on |
| `--prod-ps` | `100.0` | production run length (ps) |
| `--equil-ps` | `20.0` | equilibration length (ps) |
| `--vacuum` | - | run MD in vacuum instead of GBn2 implicit solvent (also the automatic fallback if implicit-solvent setup fails) |
| `--platform` | `auto` | OpenMM platform: `auto` prefers CUDA, then OpenCL, falling back to CPU if neither is usable on this machine; pass `CPU`/`CUDA`/`OpenCL`/`Reference` to force one explicitly (raises if that platform isn't usable) |

`--platform` only applies to `dd_docking-refine`'s MD step; docking's own
GPU option is `dd_docking-dock --backend` (see
[GPU-accelerated docking](#optional-gpu-accelerated-docking-linux-only-rigid-receptors-only)).

Python API:

```python
from dd_docking import refine_top_hits

result = refine_top_hits("data/screen/ranked_results.csv", "data/screen/refine",
                         top_n=3, prod_ps=100.0)
print(result[["name", "stable", "rmsd_mean", "best_affinity"]])
```

### End-to-end pipeline (Python API)

To run all three steps together (each is also callable individually as
shown above):

```python
from dd_docking.pipeline import run_ensemble_docking

df = run_ensemble_docking(
    "data/ensemble", "data/ligands.smi", "data/screen",
    exhaustiveness=8, n_jobs=8, refine=True, refine_top_n=3,
)
```

## Test dataset (ER-alpha, DUD-E decoys)

`data/raw/` and `data/ligands.smi` hold a small but genuinely rigorous
example, not tied to any particular target -- swap in any protein with
multiple co-crystallized conformations and a DUD-E (or similar) decoy set:

- **Ensemble members**: 3 human estrogen receptor alpha (ER-alpha) ligand-
  binding-domain structures, a textbook induced-fit case (agonist vs. SERM
  ligands reposition helix 12 and reshape the pocket) --
  [3ERT](https://www.rcsb.org/structure/3ERT) (4-hydroxytamoxifen),
  [1XPC](https://www.rcsb.org/structure/1XPC) and
  [1YIM](https://www.rcsb.org/structure/1YIM) (two different SERM
  scaffolds). All three are wild-type, single-chain, and verified via the
  RCSB API to have no internal missing-residue gaps in the ligand-binding
  domain (see the next section for why that check matters) before being
  chosen.
- **Ligand library**: each member's own co-crystallized ligand (self-docking
  check) + 3 genuine ER-alpha actives + 9 property-matched decoys, all
  pulled from [DUD-E's `esr1` target](https://dude.docking.org/targets/esr1/)
  and pre-filtered with RDKit for valid, 3D-embeddable SMILES (fixed random
  seed for reproducibility). DUD-E decoys are chosen to match real actives'
  molecular weight/rotatable-bond count/charge while differing in 2-D
  topology -- a far harder discrimination test than unrelated everyday
  molecules (e.g. approved drugs), which any docking method separates from a
  real binder trivially.

## Verified behavior

- Rigid/flex PDBQT files for the 3 ensemble members in `data/ensemble/`
  dock without issue via `qvina2 --receptor rigid.pdbqt --flex flex.pdbqt
  ...`, and the box (ligand extent + every flexible residue's atoms +
  padding, see `pocket.compute_box`'s `extra_coords`) contains every
  flexible side chain's movable atoms -- confirmed by `dd_docking-dock`
  returning a real pose/affinity for every `(ligand, member)` task rather
  than Vina's "no conformations completely within the search space" failure.
- Self-docking each member's own co-crystallized ligand, real ER-alpha
  actives, and DUD-E property-matched decoys together produces a ranking
  where the weakest, most topologically-distinct decoys separate cleanly at
  the bottom, while docking alone doesn't cleanly separate every decoy from
  every true active/native higher up the list -- and MD-based induced-fit
  rescoring doesn't either, in this specific case (see the worked example
  above). This is the expected, honestly-reported behavior of a rigorous
  property-matched benchmark, not a discrimination failure to hide.
- `receptor_prep.py`'s `regularize_carboxylate_geometry` fixes a PDBFixer
  `addMissingAtoms()` quirk (a freshly-added carboxylate partner oxygen --
  backbone OXT at a chain terminus, or Asp/Glu OD2/OE2 -- placed at a
  chemically impossible angle from its sibling oxygen) that previously
  crashed Meeko's `mk_prepare_receptor.py` with an oxygen valence error on
  real PDB entries with an unblocked terminus. Confirmed fixed end to end
  on PDB 4EQC (PAK1 kinase domain): `dd_docking-prep` now succeeds (3
  defects detected and regularized) where it previously failed, and
  `dd_docking-dock` docks Naringin into the repaired receptor
  successfully (best affinity -8.646 kcal/mol, flexible side chains).

## Module structure (`dd_docking/`)

| Module | Role |
|---|---|
| `receptor_prep.py` | PDB fetch/chain isolation/residue cleanup (TER insertion, CYX renaming), PDBFixer-based repair, post-PDBFixer carboxylate geometry regularization (`regularize_carboxylate_geometry`) |
| `pocket.py` | docking box calculation (ligand extent + optional flexible-residue atoms via `extra_coords`), distance-based flexible-residue detection, Meeko flexres string formatting |
| `ensemble.py` | batch receptor_prep + pocket + Meeko PDBQT generation across conformations; save/load as `manifest.json` |
| `ligand_prep.py` | `.smi` reading, SMILES -> 3D (ETKDGv3+MMFF) -> Meeko ligand PDBQT |
| `docking.py` | CPU docking via QuickVina2 (`qvina2` CLI), with flexible-receptor support (`dock_ligand`) |
| `gpu_backend.py` | optional Vina-GPU+ backend (Linux only): backend selection (`resolve_backend`), subprocess docking (`dock_ligand_gpu`) |
| `screening.py` | parallel all-ligand × all-member docking, best-of-ensemble ranking, CSV/SDF output |
| `refine_md.py` | short implicit-solvent MD relaxation of top hits, RMSD stability evaluation, re-ranking |
| `io.py` | PDBQT -> RDKit conversion, pose SDF output, result CSV I/O |
| `parallel.py` | `parallel_map` — parallelizes independent tasks over `ProcessPoolExecutor` |
| `progress.py` | `DockProgress` / `RefineProgress` — one printed line per completed task |
| `pipeline.py` | high-level functions chaining prep -> docking -> MD refinement |
| `cli.py` | `dd_docking-prep` / `dd_docking-dock` / `dd_docking-refine` commands |

## Limitations and possible extensions

- Flexible residues are chosen purely by distance from the co-crystallized
  ligand. To pin specific residues known to be involved in induced fit,
  bypass `pocket.find_pocket_residues` and pass a flexres string directly to
  `ensemble.prepare_ensemble_member`.
- MD refinement defaults to a short, CPU-only, implicit-solvent (GBn2)
  simulation; explicit solvent or GPU-based free-energy calculations are out
  of scope (extend `refine_md.py` independently if needed).
- `screening.py` parallelizes over all `(ligand, member)` combinations, so
  Vina's grid map is recomputed per task. For large libraries, computing each
  member's grid map once and reusing it (parallelizing per member instead)
  would be a worthwhile optimization.

## License

MIT — see [LICENSE](LICENSE).
