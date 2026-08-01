"""
Purpose: re-measure how much each arm's ranking is a molecular-size sort (D0043).
Author: Mike Hallett (with Claude Code)
Date: 2026-08-01
Input: the latest D1-D4 frames, plus each T_2 seed neighbourhood
Output: append_only/inhibition/00_outputs/blacksmith/size_correlation/
        size_correlation_<N>.csv + a printed table

WHY RE-MEASURE. D0043 established that our rankings are partly a size sort
(Spearman rho = -0.617 for T_1, -0.479 for T_3 against heavy-atom count) and
that ligand efficiency OVER-corrects (-0.938), making it a smallness sort
rather than a fix. Issue #9 item 4 asks us to normalise to LE anyway. That was
measured and rejected -- but three inputs have genuinely changed since D0043
was written on 2026-07-30, so the number deserves a re-test rather than a
citation:

  1. D0047 -- `affinity_kcal` was the CNN-selected pose's affinity, not the
     best affinity. 89% of covalent candidates were affected. T_3 and T_4's
     size correlations were computed on the pre-fix column.
  2. Ligands are now protonated at pH 7.4 rather than docked as drawn. Heavy
     atom count is unchanged; the scores are not.
  3. T_2 is five seed neighbourhoods, not one. D0043's T_2 number is ATRA only.

For T_3 and T_4 the pre-fix column survives as `affinity_kcal_rows0`, so (1)
is measured directly here rather than inferred: both are reported side by side.

THE THIRD COLUMN IS THE POINT. `residual` is the metric after regressing out
heavy-atom count -- the size-decorrelated score D0043 left open as the thing to
rank on instead. Its rho against HAC is ~0 BY CONSTRUCTION; it is reported as a
control confirming the decorrelation did what it claims, not as evidence.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import rank_shortlist as rs                      # noqa: E402
from shared import outputs as sout           # noqa: E402

DATA_ROOT = Path("/data/lab_vm/append_only/inhibition")
# Analysis outputs live under the GOVERNED root, not in the repo
# (rules/data-storage.md). See shared/outputs.py for why, and for the
# versioned-write / resolve-latest policy the append-only tree needs.
OUT = sout.Topic("blacksmith", "size_correlation")
OUT_DIR = OUT.dir
SIZE_COL = "HAC"

# (label, experiment, frame prefix, rank metric). The T_2 seed neighbourhoods
# are separate experiments; liu_2024_c3 and potter_astex are still docking and
# are skipped with a note rather than silently omitted.
POOLS = [
    ("T_1 de novo",        "01_t1_de_novo",         "D1", "vina_affinity"),
    ("T_2 atra",           "02_t2_atra_crem",       "D2", "vina_affinity"),
    ("T_2 du_xu",          "02d_t2_duxu_crem",      "D2", "vina_affinity"),
    ("T_2 guo_pfizer",     "02e_t2_guo_crem",       "D2", "vina_affinity"),
    ("T_2 liu_2024_c3",    "02b_t2_liu_c3_crem",    "D2", "vina_affinity"),
    ("T_2 potter_astex",   "02c_t2_potter_crem",    "D2", "vina_affinity"),
    ("T_2 atra degree-2",  "02_t2_atra_crem_degree2", "D2", "vina_affinity"),
    ("T_3 decoration",     "03_t3_reinvent",        "D3", "affinity_kcal"),
    ("T_4 combinatorial",  "04_t4_combinatorial",   "D4", "affinity_kcal"),
]

# The pre-D0047 column, kept on the covalent frames so the fix's effect on the
# size correlation can be measured rather than assumed.
PREFIX_COL = "affinity_kcal_rows0"


def _latest(experiment: str, prefix: str) -> Path | None:
    d = DATA_ROOT / experiment
    if not d.is_dir():
        return None
    best, best_n = None, -1
    for p in d.glob(f"{prefix}_*.parquet"):
        tail = p.stem[len(prefix) + 1:]
        if tail.isdigit() and int(tail) > best_n:
            best, best_n = p, int(tail)
    return best


def _spearman(x: pd.Series, y: pd.Series) -> tuple[float, int]:
    ok = x.notna() & y.notna()
    if ok.sum() < 3:
        return float("nan"), int(ok.sum())
    return float(stats.spearmanr(x[ok], y[ok]).statistic), int(ok.sum())


def size_decorrelated(metric: pd.Series, size: pd.Series) -> pd.Series:
    """The production decorrelation, so this measures what ranking actually does.

    Deliberately NOT a second implementation. An analysis script that
    re-derives the thing it is auditing measures its own copy and can report a
    clean result while the pipeline ships a different one -- the same shape as
    the hand-maintained drop list in `how_this_project_breaks.md`.
    """
    return rs.size_decorrelated_score(metric, size).astype(float)


def main() -> None:
    rows = []
    skipped = []
    for label, experiment, prefix, metric in POOLS:
        path = _latest(experiment, prefix)
        if path is None:
            skipped.append(f"{label}: no {prefix} frame under {experiment}")
            continue
        df = pd.read_parquet(path)
        if metric not in df.columns or SIZE_COL not in df.columns:
            skipped.append(f"{label}: {path.name} has no {metric}/{SIZE_COL} "
                           "(not docked yet)")
            continue
        # Rank sees survivors only; measuring on rejects would describe a
        # population no shortlist is drawn from.
        d = df[df["rejected_at"].isna()] if "rejected_at" in df.columns else df
        if d[metric].notna().sum() < 3:
            skipped.append(f"{label}: {path.name} carries no {metric} values yet")
            continue

        rho_raw, n = _spearman(d[SIZE_COL], d[metric])
        rho_le, _ = _spearman(d[SIZE_COL], d.get("ligand_efficiency", pd.Series(dtype=float)))
        rho_res, _ = _spearman(d[SIZE_COL], size_decorrelated(d[metric], d[SIZE_COL]))
        rho_pre = float("nan")
        if PREFIX_COL in d.columns:
            rho_pre, _ = _spearman(d[SIZE_COL], d[PREFIX_COL])

        rows.append({"pool": label, "frame": path.name, "n": n,
                     "metric": metric,
                     "rho_size_vs_metric": rho_raw,
                     "rho_size_vs_metric_preD0047": rho_pre,
                     "rho_size_vs_ligand_efficiency": rho_le,
                     "rho_size_vs_residual": rho_res,
                     "mean_HAC": float(d[SIZE_COL].mean())})

    out = pd.DataFrame(rows)
    dest = OUT.write("size_correlation", ".csv")
    out.to_csv(dest, index=False)

    print("\nSpearman rho against heavy-atom count (HAC)")
    print("D0043 baseline: T_1 raw -0.617, T_3 raw -0.479, LE -0.938\n")
    hdr = (f"{'pool':22s} {'n':>6s} {'mean':>6s} {'raw':>8s} {'preFix':>8s} "
           f"{'LE':>8s} {'resid':>8s}")
    print(hdr)
    print("-" * len(hdr))
    for _, r in out.iterrows():
        def f(v):
            return "     --" if pd.isna(v) else f"{v:8.3f}"
        print(f"{r['pool']:22s} {r['n']:6d} {r['mean_HAC']:6.1f} "
              f"{f(r['rho_size_vs_metric'])} {f(r['rho_size_vs_metric_preD0047'])} "
              f"{f(r['rho_size_vs_ligand_efficiency'])} {f(r['rho_size_vs_residual'])}")
    for s in skipped:
        print(f"  SKIPPED  {s}")
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
