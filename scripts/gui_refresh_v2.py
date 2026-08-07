"""
Purpose: put the 2.1.0 ranking and its poses into the frames the GUI reads.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06
Input: rank_v2/rank_v2_<TIER>_<N>.csv + the persisted v2 poses + the latest frames
Output: new versioned D3_/D4_ parquet frames carrying nac2_* columns

FOLLOWS `merge_nac_into_frames`'s RULE AND EXTENDS IT. That script wrote `nac_*`
columns and deliberately left `shortlist`/`shortlist_synth` alone, because
silently redefining what the GUI calls "the shortlist" makes two different
rankings indistinguishable on screen. Same reasoning one level up: the 2.1.0
scores go in as **`nac2_*`**, beside the 2.0.0 `nac_*` rather than over them.

That matters more than usual here. The two rankings genuinely disagree --
topn_viable_frac vs the old enrichment ranks at rho ~= +0.33 -- so overwriting
would silently move every molecule on screen with no way to see what changed or
to go back. Both are on the frame; the GUI can show either and the difference is
inspectable.

COLUMNS ADDED

  nac2_topn_viable_frac   fraction of the top-10 poses BY ENERGY in attack
                          geometry. The metric D0068 argued for. 41.6% of
                          molecules score 0 here, which is the population the
                          2.0.0 ranking could not see.
  nac2_enrichment_cond    P(angle | distance in window) / isotropic null --
                          the sampler-bias-cancelled form.
  nac2_enrichment_joint   the 2.0.0 quantity, recomputed on the v2 run so the
                          comparison is like-for-like.
  nac2_consensus_gnina    top-10 agreement under gnina's CNN ordering.
  nac2_consensus_autodock top-10 agreement under AutoDock's ordering.
  nac2_class_rank         rank WITHIN its warhead class among filter survivors.
  nac2_passes             cleared the within-class consensus quota.
  nac2_pose_path          the persisted top-20 poses for this molecule.
  nac2_frac_in_range      how often the sampler put the warhead in range -- the
                          quantity the joint score was silently multiplying in.

`nac2_class_rank` is per class BY DESIGN and there is no global rank column. A
single cross-class ordering re-imports the rigidity bias the per-class quota
exists to remove (D0073), so the frame does not carry one to sort on.
"""

from __future__ import annotations

import argparse
import glob
import logging
import re
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

log = logging.getLogger("gui-refresh")
DATA = Path("/data/lab_vm/append_only/inhibition")
RANK = DATA / "00_outputs/blacksmith/rank_v2"
POSES = DATA / "00_outputs/blacksmith/nac_v2_poses"

FRAMES = {"T3": ("03_t3_reinvent", "D3"), "T4": ("04_t4_combinatorial", "D4")}

COLS = {
    "topn_viable_frac": "nac2_topn_viable_frac",
    "enrichment_conditional": "nac2_enrichment_cond",
    "enrichment_joint": "nac2_enrichment_joint",
    "consensus_gnina": "nac2_consensus_gnina",
    "consensus_autodock": "nac2_consensus_autodock",
    "class_rank": "nac2_class_rank",
    "passes": "nac2_passes",
    "frac_in_range": "nac2_frac_in_range",
    "topn_best_dist": "nac2_topn_best_dist",
}


def latest(pattern: str, key) -> Path | None:
    fs = glob.glob(pattern)
    return Path(max(fs, key=key)) if fs else None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    wrote = 0
    for tier, (sub, stem) in FRAMES.items():
        r = latest(str(RANK / f"rank_v2_{tier}_*.csv"),
                   key=lambda p: int(p.rsplit("_", 1)[1].split(".")[0]))
        if r is None:
            log.warning("%s: no ranking yet, skipping", tier)
            continue
        rk = pd.read_csv(r)

        fdir = DATA / sub
        f = latest(str(fdir / f"{stem}_*.parquet"),
                   key=lambda p: int(re.search(r"_(\d+)\.parquet$", p).group(1)))
        if f is None:
            log.warning("%s: no frame under %s", tier, fdir)
            continue
        df = pd.read_parquet(f)
        n = int(re.search(r"_(\d+)\.parquet$", f.name).group(1))

        keep = {k: v for k, v in COLS.items() if k in rk.columns}
        m = rk[["ident"] + list(keep)].rename(columns=keep)
        m["nac2_pose_path"] = m.ident.map(
            lambda i: str(POSES / f"{i}.sdf") if (POSES / f"{i}.sdf").is_file() else None)

        # drop any stale nac2_* from a previous refresh BEFORE merging, or pandas
        # suffixes them to _x/_y and the GUI silently finds neither. That exact
        # bug shipped once already in export_nac_poses.
        stale = [c for c in df.columns if c.startswith("nac2_")]
        if stale:
            df = df.drop(columns=stale)
            log.info("%s: dropped %d stale nac2_ columns", tier, len(stale))

        before = len(df)
        out = df.merge(m, left_on="candidate_id", right_on="ident", how="left")
        out = out.drop(columns=["ident"])
        assert len(out) == before, f"{tier}: merge changed row count {before} -> {len(out)}"
        scored = out.nac2_topn_viable_frac.notna().sum() if "nac2_topn_viable_frac" in out else 0
        posed = out.nac2_pose_path.notna().sum()
        log.info("%s: %d rows, %d scored, %d with poses", tier, len(out), scored, posed)

        if args.dry_run:
            continue
        dest = fdir / f"{stem}_{n + 1}.parquet"
        if dest.exists():
            log.warning("%s exists; refusing to overwrite (append_only)", dest.name)
            continue
        out.to_parquet(dest, index=False)
        log.info("%s -> %s", tier, dest.name)
        wrote += 1

    print(f"\n  {wrote} frame(s) written. The GUI reads the highest-numbered "
          f"D3_/D4_ parquet.")
    print("  2.0.0's nac_* columns are untouched and sit beside the new nac2_*, "
          "because\n  the two rankings disagree (rho ~ +0.33) and overwriting "
          "would move every\n  molecule on screen with no way to see what changed.")


if __name__ == "__main__":
    main()
