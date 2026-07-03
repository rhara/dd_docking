"""Progress reporting: one concise line per completed item, printed to
stdout with flush=True (mirrors overlaylib/progress.py and embodylib's
pipeline.py console conventions)."""


class DockProgress:
    """Prints one line per docked ligand:
    `<rank> <ligand_id> best=<member_id> affinity=<kcal/mol>`
    """

    def __init__(self, enabled: bool = True):
        self.enabled = enabled
        self._count = 0

    def update(self, ligand_id: str, best_member: str, best_affinity: float) -> None:
        self._count += 1
        if not self.enabled:
            return
        print(
            f"[{self._count}] {ligand_id}  best_member={best_member}  "
            f"affinity={best_affinity:.3f}",
            flush=True,
        )


class RefineProgress:
    """Prints one line per MD-refined hit:
    `<name> <implicit/vacuum>  rmsd_mean=.. rmsd_final=.. stable=..`
    """

    def __init__(self, enabled: bool = True):
        self.enabled = enabled

    def update(self, name: str, implicit: bool, rmsd_mean: float, rmsd_final: float, stable: bool) -> None:
        if not self.enabled:
            return
        print(
            f"[MD {name}] {'implicit' if implicit else 'vacuum'}  "
            f"rmsd_mean={rmsd_mean:.2f}  rmsd_final={rmsd_final:.2f}  stable={stable}",
            flush=True,
        )
