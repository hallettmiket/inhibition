"""
Purpose: dock one ligand pool as N parallel chunks across several GPUs.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-01
Input: an experiment whose ligands are already prepared under docking/ligands_<tag>/
Output: poses in the experiment's shared pose directory + one merged frame

Run:  python scripts/dock_chunked.py --experiment 02b_t2_liu_c3_crem \
          --approach t2 --frame-prefix D2 --gpus 1,3,4,6,7

WHY THIS EXISTS. Vina-GPU in virtual-screening mode is ALL-OR-NOTHING: it
writes every pose at the END of the run. On 2026-08-01 the liu_2024_c3 pool
(16,806 ligands, mean 10.65 rotatable bonds) ran for 24 h on one card, hit the
flat `timeout=86400` in `run_vina_gpu`, and was killed with **0 of 16,806 poses
written**. A full day of GPU time produced nothing, and no frame was written.

Two things follow, and this script is the second:

1. The timeout now scales with the pool (`vina_timeout_s`). That stops a run
   being killed by a constant that was sized for a different pool.
2. A pool too large for one card is split. Five chunks on five cards turn a
   ~38 h serial run into ~8 h, and — more importantly — a failure loses one
   chunk rather than the entire pool.

ONE POSE DIRECTORY, NOT N. Every chunk writes into the experiment's single
`poses_<tag>/`. Vina names each output for its ligand and `candidate_id` is
unique across the pool, so concurrent writers cannot collide, and the existing
`collect_modes` reads the directory afterwards with no idea it was filled by
five processes. Per-chunk pose directories would have needed a merge step whose
only job was to undo the split.

CHUNK DIRECTORIES ARE SYMLINKS, IN THE SCRATCH TREE. Vina takes a directory,
so the split has to be expressed as directories — but copying 400 MB of PDBQTs
five ways to express it would be absurd, and the append-only tree is the wrong
place for a working artefact that exists for one run. They are symlinks under
`/data/lab_vm/modifiable/`, which is what that tree is for.

THE MERGE IS NOT REIMPLEMENTED HERE. It calls
`noncovalent_dock_run.merge_poses_onto_frame`, the same function `run()` uses,
with the same derived drop list and the same `_x`/`_y` assertion. Two code
paths that both "merge the docking results" is how the covalent and GROMACS
frames acquired suffixed columns nobody noticed.
"""

from __future__ import annotations

import argparse
import logging
import shutil
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import noncovalent_dock_run as ncd        # noqa: E402

log = logging.getLogger("dock-chunked")

SCRATCH = Path("/data/lab_vm/modifiable/inhibition/dock_chunks")


def chunk_ligands(ligand_dir: Path, dest_root: Path, n: int) -> list[Path]:
    """Symlink the pool into `n` chunk directories, round-robin.

    ROUND-ROBIN, NOT CONTIGUOUS BLOCKS. Ligands arrive in candidate order,
    which correlates with the CReM edit that produced them and therefore with
    size and flexibility. Contiguous blocks would hand one card a chunk of
    systematically heavier molecules and leave the others idle waiting for it;
    interleaving evens the cost out without needing to measure it.
    """
    ligands = sorted(ligand_dir.glob("*.pdbqt"))
    if not ligands:
        raise SystemExit(f"no prepared ligands under {ligand_dir}")
    if dest_root.exists():
        shutil.rmtree(dest_root)          # scratch tree: safe to rebuild
    chunks = []
    for i in range(n):
        d = dest_root / f"chunk{i:02d}"
        d.mkdir(parents=True)
        chunks.append(d)
    for i, lig in enumerate(ligands):
        (chunks[i % n] / lig.name).symlink_to(lig)
    log.info("split %d ligands into %d chunks of ~%d",
             len(ligands), n, len(ligands) // n)
    return chunks


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--experiment", required=True)
    ap.add_argument("--approach", required=True)
    ap.add_argument("--frame-prefix", required=True)
    ap.add_argument("--gpus", required=True,
                    help="comma-separated device ids, e.g. 1,3,4,6,7")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s %(name)s: %(message)s")

    gpus = [int(g) for g in args.gpus.split(",") if g.strip()]
    work = ncd.DATA_ROOT / args.experiment / "docking"
    ligand_dir = work / f"ligands_{ncd.LIGAND_PREP_TAG}"
    out_dir = work / f"poses_{ncd.LIGAND_PREP_TAG}"

    already = len(list(out_dir.glob("*.pdbqt"))) if out_dir.is_dir() else 0
    if already:
        log.warning("%d pose file(s) already present in %s; Vina will rewrite "
                    "the ones it re-docks", already, out_dir)

    chunks = chunk_ligands(ligand_dir, SCRATCH / args.experiment, len(gpus))
    for gpu, c in zip(gpus, chunks):
        n = len(list(c.glob("*.pdbqt")))
        log.info("  GPU %d <- %s (%d ligands, deadline %.1f h)",
                 gpu, c.name, n, ncd.vina_timeout_s(n) / 3600)
    if args.dry_run:
        log.info("--dry-run: not launching")
        return

    t0 = time.time()
    # One thread per chunk; each blocks on its own Vina subprocess, so the
    # threads are just a way to wait on five processes at once.
    with ThreadPoolExecutor(max_workers=len(gpus)) as ex:
        futures = {ex.submit(ncd.run_vina_gpu, c, out_dir, g): (g, c)
                   for g, c in zip(gpus, chunks)}
        failed = []
        for f, (g, c) in futures.items():
            try:
                dt = f.result()
                log.info("GPU %d finished %s in %.1f s", g, c.name, dt)
            except Exception as exc:  # noqa: BLE001
                # A chunk failing must not discard the ones that succeeded --
                # their poses are already on disk and belong in the frame.
                log.error("GPU %d FAILED on %s: %s", g, c.name, exc)
                failed.append((g, c.name, str(exc)[:200]))
    elapsed = time.time() - t0
    log.info("all chunks done in %.1f h (%d failed)", elapsed / 3600, len(failed))

    merged, out, survivors, n_docked, _ = ncd.merge_poses_onto_frame(
        experiment=args.experiment, approach=args.approach,
        frame_prefix=args.frame_prefix, out_dir=out_dir, elapsed=elapsed,
        gpu=f"chunked:{args.gpus}")

    print(f"\nChunked docking -> {out}")
    print(f"  docked successfully {n_docked} / {len(survivors)}")
    print(f"  elapsed             {elapsed / 3600:.2f} h across {len(gpus)} GPUs")
    if failed:
        print(f"  CHUNKS FAILED       {len(failed)} — {failed}")
        print("  The frame carries what succeeded; re-run to fill the rest.")


if __name__ == "__main__":
    main()
