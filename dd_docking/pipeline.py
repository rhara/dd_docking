"""High-level convenience functions chaining ensemble prep -> screening ->
(optional) MD refinement. Each stage remains independently usable via its
own module (`ensemble`, `screening`, `refine_md`) and CLI subcommand.
"""
from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import pandas as pd

from . import ensemble as ensemble_mod
from . import ligand_prep, refine_md, screening


def prepare_ensemble(members: Sequence[Tuple[str, str, str]], out_dir: str, **kwargs) -> List[ensemble_mod.EnsembleMember]:
    """Thin re-export of `ensemble.prepare_ensemble` for pipeline symmetry."""
    return ensemble_mod.prepare_ensemble(
        [(mid, Path(p), lig) for mid, p, lig in members], Path(out_dir), **kwargs
    )


def run_ensemble_docking(
    ensemble_dir: str,
    ligand_smi: str,
    out_dir: str,
    *,
    exhaustiveness: int = 16,
    n_poses: int = 5,
    seed: int = 0,
    n_jobs: int = 1,
    top_n: Optional[int] = None,
    show_progress: bool = True,
    refine: bool = False,
    refine_top_n: int = 5,
    refine_prod_ps: float = 100.0,
) -> pd.DataFrame:
    """Load a prepared ensemble (from `dd_docking-prep`'s `manifest.json`),
    dock a `.smi` ligand library against it, and optionally MD-refine the
    top hits. Returns the ranked-results DataFrame (post-refine re-ranking
    if `refine=True`, otherwise the raw Vina ranking).
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ensemble = ensemble_mod.load_manifest(Path(ensemble_dir) / "manifest.json")
    ligands = ligand_prep.read_smi(ligand_smi)

    ranked_csv = out_dir / "ranked_results.csv"
    hits = screening.screen_ensemble(
        ensemble, ligands,
        exhaustiveness=exhaustiveness, n_poses=n_poses, seed=seed, n_jobs=n_jobs,
        show_progress=show_progress, top_n=top_n,
        out_csv=str(ranked_csv), out_sdf=str(out_dir / "top_hits.sdf"),
    )
    if not refine:
        return pd.read_csv(ranked_csv)

    return refine_md.refine_top_hits(
        str(ranked_csv), str(out_dir / "refine"),
        top_n=refine_top_n, prod_ps=refine_prod_ps, show_progress=show_progress,
    )
