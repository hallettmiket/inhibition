"""
Purpose: re-derive affinity_kcal as the BEST mode, not gnina's row 0 (D0047).
Author: Mike Hallett (with Claude Code)
Date: 2026-07-31
Input: existing gnina pose SDFs under append_only/<experiment>/docking/
Output: a new D3/D4 frame with corrected affinity_kcal / cnn_* / selected_mode

Run:  /data/lab_vm/envs/dwi_cheminf/bin/python3 \
        scripts/reextract_covalent_affinity.py [--approach t3 --approach t4] [--dry-run]

NO RE-DOCKING. Every pose already exists; only the choice of WHICH pose to read
was wrong. This re-reads the SDFs and rewrites the columns.

WHY. `covalent_protocol.py` took `rows[0]` from gnina's results table and read
`affinity` off it. That table is ordered by CNN POSE SCORE, so row 0 is the
CNN-best pose, not the affinity-best one. D0011 had already demoted
`cnn_affinity` to advisory because gnina reports CNN scoring is uncalibrated
for covalent docking — so the project rejected the CNN as a ranking signal and
then let it pick the pose whose affinity it ranked on.

ALL THREE SCORES MOVE TOGETHER. Reading affinity from one pose and the CNN
values from another would describe two geometries in one row. `selected_mode`
records which pose won, so the disagreement is visible rather than inferred.

THE OLD COLUMNS ARE KEPT. `affinity_kcal_rows0` preserves what the frame said
before, because a correction that erases the thing it corrected cannot be
audited. The rank columns are NOT rewritten here — re-ranking is a separate
step, so this script cannot silently reorder a shortlist as a side effect.
"""

from __future__ import annotations

import argparse
import gzip
import logging
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import io as dio                       # noqa: E402

log = logging.getLogger("reextract-affinity")

DATA = Path("/data/lab_vm/append_only/inhibition")
APPROACHES = {
    "t3": {"experiment": "03_t3_reinvent", "prefix": "D3"},
    "t4": {"experiment": "04_t4_combinatorial", "prefix": "D4"},
}

#: SDF tag -> frame column. All read from the SAME record.
TAGS = {
    "minimizedAffinity": "affinity_kcal",
    "CNNscore": "cnn_score",
    "CNNaffinity": "cnn_affinity",
}


def read_modes(path: Path) -> list[dict]:
    """Every mode in a gnina SDF, in file order, as {tag: float}."""
    op = gzip.open if path.suffix == ".gz" else open
    out: list[dict] = []
    cur: dict = {}
    pending: str | None = None
    try:
        with op(path, "rt", errors="ignore") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith("$$$$"):
                    if cur:
                        out.append(cur)
                    cur, pending = {}, None
                    continue
                if line.startswith("> "):
                    for tag in TAGS:
                        if f"<{tag}>" in line:
                            pending = tag
                            break
                    else:
                        pending = None
                    continue
                if pending:
                    try:
                        cur[pending] = float(line.strip())
                    except ValueError:
                        pass
                    pending = None
    except OSError as exc:
        log.warning("unreadable %s: %s", path.name, exc)
        return []
    if cur:
        out.append(cur)
    return out


def best_mode(modes: list[dict]) -> tuple[dict | None, int | None]:
    """The affinity-best mode and its 1-based index."""
    scored = [(i, m) for i, m in enumerate(modes)
              if m.get("minimizedAffinity") is not None]
    if not scored:
        return None, None
    i, m = min(scored, key=lambda t: t[1]["minimizedAffinity"])
    return m, i + 1


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--approach", action="append", choices=sorted(APPROACHES))
    ap.add_argument("--dry-run", action="store_true",
                    help="report the change without writing a frame")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s %(name)s: %(message)s")

    for a in (args.approach or sorted(APPROACHES)):
        cfg = APPROACHES[a]
        df, frame = dio.latest_frame(cfg["experiment"], a)
        dock_dir = DATA / cfg["experiment"] / "docking"

        key = "dock_id" if "dock_id" in df.columns else "candidate_id"
        by_id: dict[str, tuple[dict, int, int]] = {}
        for f in sorted(list(dock_dir.glob("*.sdf"))
                        + list(dock_dir.glob("*.sdf.gz"))):
            did = f.name.split("_docked")[0]
            modes = read_modes(f)
            m, idx = best_mode(modes)
            if m is not None:
                by_id[did] = (m, idx, len(modes))

        log.info("[%s] %s: %d rows, %d pose files parsed",
                 a, frame.name, len(df), len(by_id))

        out = df.copy()
        # Preserve what the frame said before. A correction that erases the
        # thing it corrected cannot be audited afterwards.
        if "affinity_kcal" in out.columns:
            out["affinity_kcal_rows0"] = out["affinity_kcal"]
        for col in ("cnn_score", "cnn_affinity"):
            if col in out.columns:
                out[f"{col}_rows0"] = out[col]

        hit = out[key].map(lambda k: by_id.get(str(k)))
        out["selected_mode"] = [t[1] if t else pd.NA for t in hit]
        out["n_modes"] = [t[2] if t else pd.NA for t in hit]
        for tag, col in TAGS.items():
            out[col] = [t[0].get(tag) if t else pd.NA for t in hit]
        out["affinity_selection"] = "min_affinity_over_modes"

        matched = int(hit.notna().sum())
        changed = int((out["selected_mode"].fillna(1) != 1).sum())
        delta = (out["affinity_kcal"] - out.get("affinity_kcal_rows0")).dropna()
        log.info("[%s] matched %d/%d; row0 was NOT affinity-best for %d (%.1f%%)",
                 a, matched, len(out), changed, 100 * changed / max(1, matched))
        if len(delta):
            log.info("[%s] affinity change: median %+.2f, best %+.2f kcal/mol",
                     a, delta.median(), delta.min())

        if args.dry_run:
            log.info("[%s] dry run — no frame written", a)
            continue

        p = dio.write_full_frame(
            out, approach=a, experiment=cfg["experiment"],
            stage=f"{a}_reextract_affinity",
            params={"correction": "D0047",
                    "affinity_selection": "min_affinity_over_modes",
                    "was": "rows[0] of gnina's CNN-score-ordered table",
                    "n_matched": matched,
                    "n_rows0_not_affinity_best": changed,
                    "median_affinity_change_kcal":
                        float(delta.median()) if len(delta) else None,
                    "note": "scores re-read from existing poses; NO re-docking. "
                            "Ranks are NOT rewritten here — re-rank separately."},
            inputs={"frame": frame})
        print(f"\n{a}: {frame.name} -> {p.name}")
        print(f"   row0 was not affinity-best for {changed}/{matched} "
              f"({100 * changed / max(1, matched):.1f}%)")
        if len(delta):
            print(f"   median improvement {delta.median():+.2f} kcal/mol")
        print("   ranks unchanged — run 04_rank.py to rebuild the shortlist")


if __name__ == "__main__":
    main()
