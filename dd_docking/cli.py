"""Command-line entry points:
  dd_docking-prep    --member ID PDB LIG_RESNAME ... -o ensemble_dir
  dd_docking-dock     ensemble_dir ligands.smi -o out_dir
  dd_docking-refine   out_dir/ranked_results.csv -o out_dir/refine
"""
import argparse
from pathlib import Path

from . import ensemble as ensemble_mod
from . import ligand_prep, refine_md, screening


def build_prep_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dd_docking-prep",
        description=(
            "Fix N receptor conformations with PDBFixer and prepare Vina-ready "
            "rigid/flexible-side-chain PDBQT files (an 'ensemble') from them."
        ),
    )
    parser.add_argument(
        "--member", nargs=3, action="append", required=True, metavar=("ID", "PDB", "LIG_RESNAME"),
        help="One ensemble member: unique id, raw co-crystal PDB path, and its "
             "bound ligand's 3-letter residue name (used to define the docking "
             "box and pocket flexible residues). Repeat for each conformation.",
    )
    parser.add_argument("-o", "--out-dir", required=True, help="Output directory for the prepared ensemble")
    parser.add_argument("--chain", default="A", help="Receptor chain to keep (default: A)")
    parser.add_argument("--pocket-cutoff", type=float, default=5.0, help="Flexible-residue distance cutoff from the reference ligand, in Angstrom")
    parser.add_argument("--max-flexres", type=int, default=10, help="Cap on flexible residues per member, closest-to-ligand first (0 or negative = no cap; keep small, each one adds real search dimensions to Vina)")
    parser.add_argument("--box-padding", type=float, default=5.0, help="Docking box padding around the reference ligand, in Angstrom")
    parser.add_argument("--charge-model", default="gasteiger", choices=["gasteiger", "espaloma", "zero"], help="Meeko receptor partial-charge model")
    return parser


def main_prep(argv=None):
    args = build_prep_parser().parse_args(argv)
    max_flexres = args.max_flexres if args.max_flexres and args.max_flexres > 0 else None
    ensemble = ensemble_mod.prepare_ensemble(
        [(mid, Path(pdb), lig) for mid, pdb, lig in args.member],
        args.out_dir,
        chain=args.chain, pocket_cutoff=args.pocket_cutoff, max_flexres=max_flexres,
        box_padding=args.box_padding, charge_model=args.charge_model,
    )
    for m in ensemble:
        flex_note = f"{len(m.flexres)} flexible residues" if m.flexres else "rigid (no flexible residues found)"
        print(f"{m.member_id}: {flex_note} -> {m.rigid_pdbqt}", flush=True)
    print(f"\n[done] {len(ensemble)} member(s) -> {args.out_dir}/manifest.json")


def build_dock_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dd_docking-dock",
        description="Dock a ligand library against every member of a prepared ensemble and rank by best-of-ensemble affinity.",
    )
    parser.add_argument("ensemble_dir", help="Directory containing manifest.json from dd_docking-prep")
    parser.add_argument("ligand_smi", help=".smi file: 'SMILES  ID' per line")
    parser.add_argument("-o", "--out-dir", required=True, help="Output directory for ranked_results.csv / top_hits.sdf")
    parser.add_argument("--exhaustiveness", type=int, default=16)
    parser.add_argument("--n-poses", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--top-n", type=int, default=None, help="Only keep the top N ranked ligands")
    parser.add_argument("--no-progress", action="store_true")
    parser.add_argument("--n-jobs", type=int, default=1, help="Parallel worker processes for the (ligand, member) task grid (<=0 = all CPU cores, default 1 = sequential). Each worker is auto-pinned to (CPU count / n-jobs) cores -- see README for why flexible docking needs several cores per task, unlike cheap rigid docking")
    parser.add_argument(
        "--backend", default="auto", choices=["auto", "cpu", "gpu"],
        help="Docking backend: auto (default) uses Vina-GPU+ when built and the box fits its "
             "OpenCL kernel (<30 A per side), else CPU vina; cpu always uses CPU vina (every "
             "platform); gpu requires Vina-GPU+ (Linux only, see scripts/build_vina_gpu.sh) and "
             "falls back to cpu with a warning if it isn't usable for a given member",
    )
    return parser


def main_dock(argv=None):
    args = build_dock_parser().parse_args(argv)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ensemble = ensemble_mod.load_manifest(Path(args.ensemble_dir) / "manifest.json")
    ligands = ligand_prep.read_smi(args.ligand_smi)

    hits = screening.screen_ensemble(
        ensemble, ligands,
        exhaustiveness=args.exhaustiveness, n_poses=args.n_poses, seed=args.seed,
        n_jobs=args.n_jobs, show_progress=not args.no_progress, top_n=args.top_n,
        backend=args.backend,
        out_csv=str(out_dir / "ranked_results.csv"), out_sdf=str(out_dir / "top_hits.sdf"),
    )
    print(f"\n[done] {len(hits)} ligand(s) ranked -> {out_dir}/ranked_results.csv")


def build_refine_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dd_docking-refine",
        description="Short implicit-solvent OpenMM MD rescoring of the top-ranked docking hits (induced-fit-style pocket relaxation).",
    )
    parser.add_argument("ranked_csv", help="ranked_results.csv from dd_docking-dock")
    parser.add_argument("-o", "--out-dir", required=True, help="Output directory for per-hit MD workdirs and md_rescore.csv")
    parser.add_argument("--top-n", type=int, default=5)
    parser.add_argument("--prod-ps", type=float, default=100.0, help="Production MD length per hit, in picoseconds")
    parser.add_argument("--equil-ps", type=float, default=20.0)
    parser.add_argument("--vacuum", action="store_true", help="Use vacuum instead of GBn2 implicit solvent")
    parser.add_argument("--no-progress", action="store_true")
    parser.add_argument(
        "--platform", default="auto", choices=["auto", "CPU", "CUDA", "OpenCL", "Reference"],
        help="OpenMM platform for MD (default: auto -- prefers CUDA, then OpenCL, falls back to "
             "CPU if no GPU platform is usable on this machine; pass CPU/CUDA/OpenCL/Reference "
             "to force one explicitly)",
    )
    return parser


def main_refine(argv=None):
    args = build_refine_parser().parse_args(argv)
    out = refine_md.refine_top_hits(
        args.ranked_csv, args.out_dir, top_n=args.top_n,
        implicit=not args.vacuum, prod_ps=args.prod_ps, equil_ps=args.equil_ps,
        show_progress=not args.no_progress, platform=args.platform,
    )
    print(f"\n[done] {len(out)} hit(s) refined -> {args.out_dir}/md_rescore.csv")
