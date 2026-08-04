"""
Purpose: cluster each candidate's docked modes and pick a REAL representative pose.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-04
Input: a shortlist's pose files + the prepared receptor's pocket residues
Output: 00_outputs/blacksmith/pose_clusters/pose_clusters_<N>.csv

Issue #14, the clustering half. Groups a candidate's modes by pocket-contact
profile and returns the **medoid** of each cluster -- an index into the modes
that were actually generated, never a blended structure. That was #14's own
conclusion once the averaging hazard was raised: weight WHICH real pose is
chosen, never the geometry.

WHAT THE CLUSTER POPULATIONS ARE NOT. Within ONE Vina run they are not a
confidence measure and must never be read as one. Vina diversifies its reported
modes with a minimum-RMSD floor before writing them, and our own median
nearest-neighbour RMSD of 1.176 A sits right on that floor -- so counting how
many of the nine agree measures the OUTPUT FORMATTER, not agreement. That is
the trap flagged in #10 and it is the reason `bpmd_confidence` in the output is
NA rather than filled with something plausible.

The honest confidence sources, neither of which exists yet:

* **BPMD** (#14) -- needs PLUMED, which is not installed on this machine.
  GROMACS 2026.3 carries the native `-plumed` flag but there is no libplumed to
  load, and `openmmplumed` is absent.
* **Replicate-run convergence** (#10 item 2) -- cluster population across
  INDEPENDENT runs, which is meaningful because the seeds differ. Vina-GPU
  already draws a fresh seed per invocation, so this needs no code change, only
  compute.

Until one of those lands, the columns here are DESCRIPTIVE GEOMETRY. Per #13's
pre-registration rule, nothing here ranks, gates or filters anything.

ONE THING THE DESCRIPTIVE COLUMNS ALREADY SHOW. `top_mode_is_medoid` records
whether Vina's best-scoring mode is also the geometric consensus of its own
nine. Where it is False, the pose the pipeline has been carrying is not the one
the rest of the modes agree on -- which is worth seeing next to D0046's finding
that top-1 recovers the crystal pose 22.5% of the time while best-of-9 reaches
55%.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "integration" / "app"))

from shared import io as dio                       # noqa: E402
from shared import outputs as sout                 # noqa: E402
from shared import pose_vector as pv               # noqa: E402
import pose3d as p3d                               # noqa: E402

log = logging.getLogger("cluster-poses")

OUT = sout.Topic("blacksmith", "pose_clusters")
DATA = Path("/data/lab_vm/append_only/inhibition")

# Contact-profile distance below which two modes are the same binding mode.
# 2.0 A of summed per-residue contact difference. Chosen BEFORE looking at any
# ranking, and reported alongside the result at 1.0 and 3.0 so a reader can see
# how sensitive the cluster count is rather than taking one number on trust.
DEFAULT_THRESHOLD = 2.0
REPORT_THRESHOLDS = (1.0, 2.0, 3.0)


def receptor_residues() -> tuple[tuple[int, ...], dict[int, np.ndarray]]:
    """Pocket residues' heavy-atom coordinates, from the prepared receptor."""
    if not p3d.receptor_readable():
        raise SystemExit(
            f"receptor not readable: {p3d.RECEPTOR}. Presence is not "
            "readability on this data root — see D0054.")
    resi = tuple(p3d.pocket_resi())
    want, out = set(resi), {}
    for ln in p3d.RECEPTOR.read_text(errors="replace").splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        el = (ln[76:78].strip() or ln[12:16].strip()[:1]).upper()
        if el == "H":
            continue
        try:
            r = int(ln[22:26])
            xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
        except ValueError:
            continue
        if r in want:
            out.setdefault(r, []).append(xyz)
    got = {k: np.array(v) for k, v in out.items()}
    missing = [r for r in resi if r not in got]
    if missing:
        # Recorded, not silently tolerated: a shorter basis changes every
        # distance and nothing downstream could tell.
        log.warning("%d pocket residue(s) absent from the receptor: %s",
                    len(missing), missing)
    return resi, got


def vectors_for(pose_file: Path, resi, rec) -> list[pv.PoseVector]:
    out = []
    for pose in p3d.read_poses(pose_file):
        atoms = p3d.pose_atoms(pose)
        if not atoms:
            continue
        xyz = np.array([(x, y, z) for _, x, y, z in atoms], dtype=float)
        out.append(pv.contact_vector(xyz, rec, resi))
    return out


def cluster_one(pose_file: Path, resi, rec,
                threshold: float = DEFAULT_THRESHOLD) -> dict | None:
    vs = vectors_for(pose_file, resi, rec)
    if len(vs) < 2:
        return None
    labels = pv.cluster(vs, threshold)
    medoid = pv.representative(vs)                 # index into REAL modes
    sizes = pd.Series(labels).value_counts()
    largest = int(sizes.idxmax())
    members = [i for i, l in enumerate(labels) if l == largest]
    row = {
        "candidate_id": pose_file.stem.replace("_out", "").replace("_docked", ""),
        "n_modes": len(vs),
        "n_clusters": len(set(labels)),
        "cluster_threshold": threshold,
        # 1-BASED, to match how the GUI and Vina both number modes. Recorded as
        # a mode NUMBER, never as a blended structure.
        "medoid_mode": medoid + 1,
        "largest_cluster_size": int(sizes.max()),
        "largest_cluster_modes": ",".join(str(i + 1) for i in members),
        "medoid_in_largest_cluster": bool(labels[medoid] == largest),
        # Vina sorts by affinity, so mode 1 is its best. Is it the geometric
        # consensus of its own nine?
        "top_mode_is_medoid": bool(medoid == 0),
        "top_mode_cluster_size": int(sizes.get(labels[0], 0)),
        # No confidence yet, and NOT invented — see the module docstring.
        "bpmd_confidence": pd.NA,
        "confidence_source": pd.NA,
    }
    for t in REPORT_THRESHOLDS:
        row[f"n_clusters_at_{t:g}"] = len(set(pv.cluster(vs, t)))
    return row


def resolve_poses(approach: str, df, ids: set[str]) -> list[tuple[str, Path]]:
    """(candidate_id, pose file) using the frame's OWN answer.

    THE FILENAME CANNOT BE DERIVED FROM THE CANDIDATE ID. The covalent arms
    name pose files by a separate `dock_id` sharing no hash with
    `candidate_id`; `pose3d.find_pose` says so in as many words, and building
    the name by hand here found 0 of 25 T_3 shortlist candidates -- exactly the
    "overlap is exactly zero across 4080 files" its docstring records. Use the
    resolver that owns the convention, and pass the frame's `dock_id` and
    `pose_path` so the lookup is the frame's answer rather than a second one.
    """
    sub = df[df["candidate_id"].isin(ids)]
    out = []
    for _, r in sub.iterrows():
        p = p3d.find_pose(approach, r["candidate_id"],
                          dock_id=r.get("dock_id"),
                          pose_path=r.get("pose_path"))
        if p is not None:
            out.append((r["candidate_id"], p))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--approach", default=None, help="t1|t2|t3|t4 (default: all)")
    ap.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD)
    ap.add_argument("--limit", type=int, default=None,
                    help="cap candidates per approach (smoke runs)")
    ap.add_argument("--all-candidates", action="store_true",
                    help="every ranked candidate, not just the shortlist")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    resi, rec = receptor_residues()
    log.info("pocket basis: %d residues", len(resi))

    import data as appdata                          # noqa: PLC0415
    approaches = [args.approach] if args.approach else list(appdata.APPROACHES)

    rows = []
    for a in approaches:
        df, name = appdata.load_frame(a)
        if df is None:
            log.warning("%s: %s", a, name)
            continue
        if args.all_candidates:
            ids = set(df["candidate_id"])
        else:
            col = appdata.shortlist_column(df)
            ids = set(df.loc[df[col] == True, "candidate_id"])  # noqa: E712
        files = resolve_poses(a, df, ids)
        if args.limit:
            files = files[:args.limit]
        log.info("%s: %d candidates selected, %d pose files resolved",
                 a, len(ids), len(files))
        for cid, f in files:
            try:
                r = cluster_one(f, resi, rec, args.threshold)
                if r:
                    r["candidate_id"] = cid      # from the frame, not the path
            except Exception as exc:  # noqa: BLE001
                log.warning("%s: %s", f.name, str(exc)[:100])
                continue
            if r:
                r["approach"] = a
                rows.append(r)

    if not rows:
        raise SystemExit("no poses clustered — check the shortlist and pose dirs")

    out = pd.DataFrame(rows)
    dest = OUT.write("pose_clusters", ".csv")
    out.to_csv(dest, index=False)

    print(f"\npose clusters -> {dest}")
    print(f"  candidates            {len(out)}")
    print(f"  median clusters/cand  {out['n_clusters'].median():.1f} "
          f"of {out['n_modes'].median():.0f} modes")
    top_is_medoid = out["top_mode_is_medoid"].mean()
    print(f"  top-scored mode IS the geometric medoid: "
          f"{top_is_medoid:.1%} of candidates")
    print(f"  medoid inside the largest cluster:       "
          f"{out['medoid_in_largest_cluster'].mean():.1%}")
    print("\n  cluster count by threshold (sensitivity, not a result):")
    for t in REPORT_THRESHOLDS:
        print(f"    {t:g} A: median {out[f'n_clusters_at_{t:g}'].median():.1f}")
    print("\n  bpmd_confidence is NA: PLUMED is not installed on this machine, "
          "and intra-run\n  cluster population measures Vina's diversification "
          "filter, not agreement (#10).")


if __name__ == "__main__":
    main()
