#!/usr/bin/env python3
"""
Purpose: cluster each molecule's pose cloud in 3N coordinate space with HDBSCAN, and count what comes out.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-19
Input: a screen topic's persisted pose clouds, <topic>_allposes/<ident>.sdf
Output: append_only/00_outputs/blacksmith/coord_modes/ -- per-mode and per-molecule tables

@tt8804: *"we need to use HDBSCAN on only the molecules dimensions in 3d space
(3 x atoms dimensions) so that we generate clusters that are essentially the same
poses being recreated."*

WHAT THIS ANSWERS, AND WHAT IT DELIBERATELY DOES NOT. Three numbers: how many
modes a molecule gets, how many poses are in one, and how many poses belong to no
mode at all (ORPHANS -- HDBSCAN's noise label). Nothing here scores a mode.
D0088's whole point is that the shipped splitter formed groups along the axis it
then graded them on, so a rule that is judged by the grade it produces cannot be
compared against it. This measures the CLUSTERING, alone.

HEAVY ATOMS, NO SUPERPOSITION. Every pose is already in the receptor's frame.
Fitting poses onto each other would ask whether they are the same SHAPE; the
question is whether they are in the same PLACE. Hydrogens are dropped because a
rotating methyl is not a different binding mode, and meeko's polar-hydrogen
placement is not a measurement.

WIDTHS ARE REPORTED IN RMSD, NOT IN THE CLUSTERING METRIC. The clustering runs on
raw Euclidean distance in 3N space, which is `sqrt(n_heavy)` times RMSD -- about
7x for a 50-atom ligand. Reporting the clustering's own units would put a 1 A
group on the page as 7 A and invite exactly the comparison-against-the-old-rule
that this experiment exists to make. Every width printed here is heavy-atom RMSD,
the same unit D0086 and D0088 measured mode span in.

WHAT THE INPUT IS, AND THE ONE CAVEAT ON IT. `<topic>_allposes/` holds the poses
the SHIPPED rule assigned to a mode -- `nac_screen_v2` writes
`labels[i] in mode_ids` -- so the ~16% of each cloud that the old DBSCAN called
noise never reached disk. The orphan rate below is therefore measured on a cloud
that was already filtered once, by a different rule, and is a LOWER BOUND on what
this rule would call orphan given all 500 poses. Fixing that needs the screen to
persist its noise too, which is a change to the screen, not to this experiment.
"""

from __future__ import annotations

import argparse
import glob
import json
import logging
import multiprocessing as mp
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from shared import compute as cp                           # noqa: E402
from shared import outputs as sout                         # noqa: E402
from shared import pose_cluster as pclust                  # noqa: E402
from shared import run_paths as rp                         # noqa: E402

log = logging.getLogger("coord-modes")


def heavy_coords(sdf: Path) -> np.ndarray:
    """(n_poses, n_heavy, 3) from one all-poses SDF.

    Read by RECORD rather than as a multi-conformer molecule: `write_sdf` emits
    the same molecule once per pose, so the supplier yields one single-conformer
    mol per pose. A pose whose record will not sanitise is dropped and counted,
    never silently replaced by its neighbour.

    NOT keyed on `pose_idx`. The clouds written before #76 do not carry it, and
    nothing here needs to join back to the per-pose table -- coordinates are the
    whole input. A loader that required the key would refuse the library this
    was asked about.
    """
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    out = []
    for m in Chem.SDMolSupplier(str(sdf), removeHs=False, sanitize=True):
        if m is None or m.GetNumConformers() != 1:
            continue
        pos = m.GetConformer().GetPositions()
        keep = [a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum() > 1]
        out.append(pos[keep])
    if not out:
        return np.empty((0, 0, 3))
    n = {len(c) for c in out}
    if len(n) != 1:
        raise ValueError(f"{sdf.name}: {sorted(n)} heavy atoms across records")
    return np.array(out)


def widths(coords: np.ndarray, labels: np.ndarray) -> dict[int, float]:
    """Widest heavy-atom RMSD between any two poses of a mode, per mode.

    The DIAMETER, not the mean spread: "essentially the same pose" is a claim
    about the worst pair in the group, and a mean hides one outlier in twenty.
    """
    out = {}
    for c in sorted({int(x) for x in labels if x >= 0}):
        m = coords[labels == c]
        if len(m) < 2:
            out[c] = 0.0
            continue
        d = np.sqrt(((m[:, None] - m[None, :]) ** 2).sum(axis=3).mean(axis=2))
        out[c] = float(d.max())
    return out


def one(args) -> dict | None:
    sdf, min_cluster_size, selection = args
    cp.pin_to_one_thread()
    ident = sdf.stem
    try:
        coords = heavy_coords(sdf)
    except Exception as exc:                                   # noqa: BLE001
        log.warning("%s: %s", ident, exc)
        return None
    if not len(coords):
        return None
    lab = pclust.cluster_coords(coords, min_cluster_size=min_cluster_size,
                                selection=selection)
    w = widths(coords, lab)
    sizes = [int((lab == c).sum()) for c in sorted(w)]
    return {"ident": ident, "n_poses": int(len(coords)),
            "n_heavy": int(coords.shape[1]),
            "n_modes": len(sizes),
            "n_orphan": int((lab < 0).sum()),
            "mode_sizes": sizes,
            "mode_widths": [w[c] for c in sorted(w)]}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--topic", default=rp.topic(),
                    help="screen topic whose <topic>_allposes/ clouds to read")
    ap.add_argument("--min-cluster-size", type=int,
                    default=pclust.MIN_CLUSTER_SIZE)
    ap.add_argument("--selection", default=pclust.SELECTION,
                    choices=["leaf", "eom"])
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--workers", type=int, default=0)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        datefmt="%H:%M:%S")

    cloud_dir = rp.BLACKSMITH / f"{args.topic}_allposes"
    sdfs = sorted(Path(p) for p in glob.glob(str(cloud_dir / "*.sdf")))
    if args.limit:
        sdfs = sdfs[:args.limit]
    if not sdfs:
        log.error("no clouds under %s", cloud_dir)
        return 1
    log.info("%d clouds under %s", len(sdfs), cloud_dir)

    # HOW MUCH OF THE DOCKED CLOUD REACHED DISK. Reported, not assumed: the
    # screen persists only mode-assigned poses, and the orphan rate below is
    # conditional on that filter. If the table is missing the run still works
    # and the coverage line says so rather than printing a fabricated 100%.
    docked = None
    tabs = glob.glob(str(rp.BLACKSMITH / args.topic / "poses_s*_*.csv"))
    if tabs:
        t = pd.concat([pd.read_csv(f, usecols=["ident", "mode"]) for f in tabs],
                      ignore_index=True)
        docked = t.groupby("ident").size()

    n_workers = cp.available_workers(args.workers or None)
    os.nice(cp.NICE)
    payload = [(p, args.min_cluster_size, args.selection) for p in sdfs]
    with mp.Pool(n_workers) as pool:
        rows = [r for r in pool.imap_unordered(one, payload, chunksize=4)
                if r is not None]

    per_mode = pd.DataFrame([
        {"ident": r["ident"], "mode": i, "size": s, "width_rmsd": w}
        for r in rows for i, (s, w) in
        enumerate(zip(r["mode_sizes"], r["mode_widths"]))])
    per_mol = pd.DataFrame([{k: v for k, v in r.items()
                             if k not in ("mode_sizes", "mode_widths")}
                            for r in rows])
    per_mol["orphan_frac"] = per_mol.n_orphan / per_mol.n_poses
    if docked is not None:
        per_mol["n_docked"] = per_mol.ident.map(docked)
        per_mol["cloud_coverage"] = per_mol.n_poses / per_mol.n_docked

    out = sout.Topic("blacksmith", "coord_modes")
    p_mode = out.write(f"modes_{args.topic}", ".csv")
    p_mol = out.write(f"molecules_{args.topic}", ".csv")
    per_mode.to_csv(p_mode, index=False)
    per_mol.to_csv(p_mol, index=False)

    tot_poses = int(per_mol.n_poses.sum())
    tot_orphan = int(per_mol.n_orphan.sum())
    tot_modes = int(per_mol.n_modes.sum())
    summary = {
        "topic": args.topic, "min_cluster_size": args.min_cluster_size,
        "selection": args.selection,
        "molecules": len(per_mol), "poses_clustered": tot_poses,
        "modes": tot_modes,
        "modes_per_molecule_mean": float(per_mol.n_modes.mean()),
        "modes_per_molecule_median": float(per_mol.n_modes.median()),
        "poses_per_mode_mean": (tot_poses - tot_orphan) / max(tot_modes, 1),
        "poses_per_mode_median": float(per_mode["size"].median()),
        "poses_per_mode_max": int(per_mode["size"].max()),
        "orphans": tot_orphan,
        "orphan_frac": tot_orphan / max(tot_poses, 1),
        "width_rmsd_median": float(per_mode.width_rmsd.median()),
        "width_rmsd_p90": float(per_mode.width_rmsd.quantile(0.90)),
        "width_rmsd_max": float(per_mode.width_rmsd.max()),
        "cloud_coverage": (float(per_mol.cloud_coverage.mean())
                           if "cloud_coverage" in per_mol else None),
    }
    (out.dir / p_mode.name.replace("modes_", "summary_")
     .replace(".csv", ".json")).write_text(json.dumps(summary, indent=2))

    print(f"\n  HDBSCAN in 3N coordinate space -- {args.topic}, "
          f"min_cluster_size {args.min_cluster_size}, {args.selection}\n")
    print(f"  molecules                 {summary['molecules']:>10,}")
    print(f"  poses clustered           {tot_poses:>10,}"
          + (f"   ({100*summary['cloud_coverage']:.0f}% of the docked cloud"
             " reached disk)" if summary["cloud_coverage"] else ""))
    print(f"  modes                     {tot_modes:>10,}")
    print(f"  modes per molecule        {summary['modes_per_molecule_mean']:>10.1f}"
          f"   median {summary['modes_per_molecule_median']:.0f}")
    print(f"  poses per mode            {summary['poses_per_mode_mean']:>10.1f}"
          f"   median {summary['poses_per_mode_median']:.0f}"
          f"   max {summary['poses_per_mode_max']}")
    print(f"  ORPHAN poses              {tot_orphan:>10,}"
          f"   {100*summary['orphan_frac']:.1f}% of the cloud")
    print(f"\n  mode width (heavy-atom RMSD, widest pair)")
    print(f"    median {summary['width_rmsd_median']:.2f} A"
          f"   p90 {summary['width_rmsd_p90']:.2f} A"
          f"   max {summary['width_rmsd_max']:.2f} A")
    print(f"\n  -> {p_mode}\n  -> {p_mol}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
