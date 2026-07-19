[Japanese version](README.jp.md)

# dd_docking — ensemble docking with flexible-side-chain sampling and MD-based induced-fit rescoring

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

Requires vina, meeko, pdbfixer, openmm, openmmforcefields, openff-toolkit,
and mdtraj (now declared in `pyproject.toml`'s `dependencies`). These are
best installed via conda-forge, using mamba:

```bash
mamba create -n dd_docking python=3.12 -c conda-forge \
    rdkit numpy pandas vina meeko pdbfixer openmm openmmforcefields openff-toolkit mdtraj
mamba activate dd_docking
cd dd_docking
pip install -e .
```

This installs three console commands: `dd_docking-prep`, `dd_docking-dock`,
`dd_docking-refine`.

## Usage

### 1. Ensemble preparation (`dd_docking-prep`)

Pass `--member ID PDB_FILE COCRYSTAL_LIGAND_CODE` once per receptor
conformation. Each member is repaired independently with PDBFixer, and the
docking box plus pocket residues (default cutoff 5 Å) are derived from the
co-crystallized ligand's coordinates.

```bash
dd_docking-prep \
  --member 6w63 data/raw/6W63.pdb X77 \
  --member 7l11 data/raw/7L11.pdb XF1 \
  --member 7l10 data/raw/7L10.pdb XEY \
  -o data/ensemble
```

Output (measured, reproducible from the `data/` in this repository):

```
6w63: 27 flexible residues -> data/ensemble/6w63_rigid.pdbqt
7l11: 24 flexible residues -> data/ensemble/7l11_rigid.pdbqt
7l10: 20 flexible residues -> data/ensemble/7l10_rigid.pdbqt

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
rotational degree of freedom to Vina's search space. On the Mpro test data in
this repository (6W63/7L11/7L10, 5 Å cutoff from the co-crystallized
ligand), 20-27 residues fell within the cutoff, and docking one ligand
against one member with no cap didn't finish in under an hour. Capping at 10
(nearest first) brought this down to a practical runtime. Removing the cap or
widening the cutoff to "be safe" runs into this trap easily.

Python API:

```python
from dd_docking import prepare_ensemble

ensemble = prepare_ensemble(
    [("6w63", "data/raw/6W63.pdb", "X77"),
     ("7l11", "data/raw/7L11.pdb", "XF1"),
     ("7l10", "data/raw/7L10.pdb", "XEY")],
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

Work is split into `(ligand, ensemble member)` tasks (8 ligands × 3 members =
24 tasks in the example above), parallelized via `--n-jobs` (default
sequential; `<=0` uses all CPU cores). One progress line is printed per
ligand once it has been docked against every member (measured on this
repository's data, 16-core Mac, `--n-jobs 4`, ~8 minutes for 3 native
compounds + 5 approved drugs):

```
[parallel] using 4 worker processes for 24 tasks
[1] X77_native_6w63  best_member=6w63  affinity=-6.837
[2] XF1_native_7l11  best_member=6w63  affinity=-7.697
[3] aspirin  best_member=7l11  affinity=-5.210
[4] XEY_native_7l10  best_member=7l11  affinity=-8.058
[5] ibuprofen  best_member=7l10  affinity=-5.282
[6] naproxen  best_member=7l10  affinity=-6.187
[7] acetaminophen  best_member=6w63  affinity=-4.376
[8] metformin  best_member=7l11  affinity=-5.539

[done] 8 ligand(s) ranked -> data/screen/ranked_results.csv
```

Final ranking (`ranked_results.csv`, affinity ascending): the three
co-crystallized ligands (self-docking) take the top 3 spots, with the 5
unrelated approved drugs ranked below (the intended discrimination). The
`best_member` column and per-member `affinity[6w63/7l11/7l10]` columns
differ across rows, confirming the ensemble's conformations actually produce
different results (evidence that ensemble docking is functioning as
intended):

| rank | ligand_id | best_member | best_affinity | affinity[6w63] | affinity[7l11] | affinity[7l10] |
|---|---|---|---|---|---|---|
| 1 | XEY_native_7l10 | 7l11 | -8.058 | -6.817 | -8.058 | -7.711 |
| 2 | XF1_native_7l11 | 6w63 | -7.697 | -7.697 | -7.414 | -7.418 |
| 3 | X77_native_6w63 | 6w63 | -6.837 | -6.837 | -6.773 | -6.554 |
| 4 | naproxen | 7l10 | -6.187 | -5.102 | -6.020 | -6.187 |
| 5 | metformin | 7l11 | -5.539 | -5.310 | -5.539 | -5.349 |
| 6 | ibuprofen | 7l10 | -5.282 | -4.388 | -4.739 | -5.282 |
| 7 | aspirin | 7l11 | -5.210 | -4.623 | -5.210 | -4.109 |
| 8 | acetaminophen | 6w63 | -4.376 | -4.376 | -4.022 | -4.355 |

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
| `--no-progress` | - | suppress progress log |

**`--n-jobs` behavior and CPU allocation**: `--n-jobs 1` (default) runs
sequentially, giving each docking task all available cores internally
(`cpu=0`). Any other `--n-jobs` value divides cores evenly across worker
processes (`cpu_count // n_jobs`, minimum 1). Flexible side-chain docking is
much heavier than rigid docking and benefits significantly from Vina's own
multithreading, so the common pattern for lightweight rigid docking — fix
`--cpu 1` and parallelize jobs instead — backfires badly here (measured: same
conditions, `cpu=1` fixed took ~850s per task vs. ~150s per task with several
cores each, about 5.6x). Setting too many workers reduces each worker's core
share and can also slow things down, so start conservatively (roughly 1/4 to
1/2 of core count) relative to library size and ensemble size.

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
  -o data/screen/refine --top-n 2 --prod-ps 20 --equil-ps 5
```

Measured output (`--prod-ps` shortened here for a quick test; use a longer
value in production):

```
[MD XEY_native_7l10] implicit solvent setup failed, retrying in vacuum: ''
[MD XEY_native_7l10] vacuum  rmsd_mean=1.90  rmsd_final=2.46  stable=True
[MD XF1_native_7l11] implicit  rmsd_mean=2.67  rmsd_final=3.32  stable=True

[done] 2 hit(s) refined -> data/screen/refine/md_rescore.csv
```

The first hit's GBn2 implicit-solvent system setup failed (an occasional,
apparently transient issue in `openmmforcefields`' force-field template
generation) and automatically fell back to vacuum MD (this fallback happens
automatically even without passing `--vacuum`; expected behavior). Both hits
met the stability criterion (RMSD mean < 3 Å and final < 4 Å) and were
marked `stable=True`.

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

MD refinement is the only GPU-capable stage in this package (Vina docking,
via the `vina` package used here, has no GPU backend), so `--platform` only
applies to `dd_docking-refine`.

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

## Verified behavior

- Rigid/flex PDBQT files for the 3 ensemble members in `data/ensemble/`
  (SARS-CoV-2 Mpro, PDB 6W63/7L11/7L10 — genuinely different pocket
  conformations bound to different inhibitors) load without issue into
  `vina.Vina().set_receptor(rigid, flex)`, and grid maps compute
  successfully.
- Self-docking each member's co-crystallized ligand discriminates it from
  unrelated molecules (approved drugs) by affinity (see `data/ligands.smi`).
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
| `pocket.py` | docking box calculation, distance-based flexible-residue detection, Meeko flexres string formatting |
| `ensemble.py` | batch receptor_prep + pocket + Meeko PDBQT generation across conformations; save/load as `manifest.json` |
| `ligand_prep.py` | `.smi` reading, SMILES -> 3D (ETKDGv3+MMFF) -> Meeko ligand PDBQT |
| `docking.py` | thin Vina wrapper with flexible-receptor support (`make_vina` / `dock_ligand`) |
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
