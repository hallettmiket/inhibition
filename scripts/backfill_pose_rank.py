"""
Purpose: stamp `pose_rank` onto v2 SDFs written before the writer set it.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06
Input: 00_outputs/blacksmith/nac_v2_poses/*.sdf
Output: the same files, with a pose_rank property per record

`bpmd_run.read_pose` selects a pose by its `pose_rank` PROPERTY and refuses to
take one by file position -- correctly, since a pose identified by where it sits
is a pose a re-sort silently redefines. The v2 screen's first writer did not
stamp it, so poses already on disk cannot be addressed by rank.

Rank is assigned by FILE ORDER here, which is safe only because the writer wrote
them in energy order and nothing has re-sorted them since. That assumption is
stated rather than hidden, and the file is skipped if it already carries ranks.
"""
from __future__ import annotations
import sys
from pathlib import Path
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

POSES = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith/nac_v2_poses")


def main() -> None:
    files = sorted(POSES.glob("*.sdf"))
    done = skipped = failed = 0
    for f in files:
        mols = [m for m in Chem.SDMolSupplier(str(f), removeHs=False) if m is not None]
        if not mols:
            failed += 1
            continue
        if mols[0].HasProp("pose_rank"):
            skipped += 1
            continue
        tmp = f.with_suffix(".sdf.tmp")
        w = Chem.SDWriter(str(tmp))
        for rank, m in enumerate(mols, 1):
            m.SetProp("pose_rank", str(rank))
            m.SetProp("energy_rank", str(rank))
            w.write(m)
        w.close()
        tmp.replace(f)
        done += 1
    print(f"  stamped {done}, already had ranks {skipped}, unreadable {failed}, "
          f"of {len(files)} SDFs")


if __name__ == "__main__":
    sys.exit(main())
