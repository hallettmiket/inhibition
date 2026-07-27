"""
Purpose: The D^i column contract every approach must satisfy, and its validator.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: an approach's candidate frame
Output: validation errors, or a conforming top-10 hand-off frame

WHY A CONTRACT AT ALL. Four approaches are written independently and their
outputs are pooled by a GUI none of them know about. Without a shared column
contract the integration phase becomes a negotiation with four separate authors
about what `score` meant — and the whole design rests on the pooled shortlist
being comparable on the axes where comparison is legitimate.

WHAT IS AND IS NOT COMPARABLE (this is the load-bearing part):

  COMPARABLE across all four — the RDKit physicochemical axes plus SAscore and
  novelty. Same tool, same call, same reference set, so the numbers mean the
  same thing everywhere. These are what the GUI plots across the pooled 40.

  NOT COMPARABLE — every approach's rank metric. T_1 and T_2 rank on Vina
  kcal/mol (lower better), T_3 on gnina affinity (kcal/mol, lower better, but a
  different program and a covalent complex), T_4 on MM-GBSA. So `rank_metric` is
  carried as a NAME plus a VALUE plus a DIRECTION, never as a bare number, and
  the GUI is required to display the name and direction beside the value.

  INTERNAL ONLY — T_2's `rank_score` is a weight-vector artifact and must never
  cross the join.

STAMP, DO NOT DELETE. `Di.parquet` keeps every candidate ever enumerated,
including those stamped `rejected_at`. Only `Di_top10.csv` is filtered. A
validator that accepted a trimmed full frame would silently permit the one thing
the design forbids.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd
import yaml

log = logging.getLogger(__name__)

_REPO_ROOT = Path(__file__).resolve().parent.parent

# --- the top-10 hand-off contract ------------------------------------------
# Exactly these columns, in this order. The GUI pools four of these frames.
TOP10_COLUMNS: tuple[str, ...] = (
    "candidate_id",
    "canonical_smiles",          # THE join key
    "approach",                  # t1 | t2 | t3 | t4
    "mechanism",                 # covalent | non_covalent
    "rank",                      # 1..N within the approach
    "rank_metric_name",          # e.g. vina_affinity, gnina_affinity_kcal
    "rank_metric_value",
    "rank_metric_direction",     # lower_is_better | higher_is_better
    "uncertainty",
    "novelty_external",
    # cross-approach comparable axes — identical RDKit computation everywhere
    "MW", "cLogP", "TPSA", "QED", "HBD", "HBA", "rot_bonds", "frac_sp3", "SAscore",
    "alert_flags",
    "conditions_resolved",       # e.g. "i,iii,iv"
    "inhibition_proxy_strength", # weak | strong
    "pose_path",
    "manifest_ref",
)

REQUIRED_FULL_COLUMNS: tuple[str, ...] = (
    "candidate_id", "canonical_smiles", "approach", "rejected_at",
)

VALID_APPROACHES = {"t1", "t2", "t3", "t4"}
VALID_MECHANISMS = {"covalent", "non_covalent"}
VALID_DIRECTIONS = {"lower_is_better", "higher_is_better"}
VALID_PROXY = {"weak", "strong"}

# Never permitted to cross the join, whatever an approach computes internally.
INTERNAL_ONLY_COLUMNS = ("rank_score", "axes_covered", "panel_cycle", "w_vector")


class SchemaError(RuntimeError):
    """A frame does not satisfy the contract."""


@dataclass
class ValidationReport:
    """Everything wrong with a frame, gathered rather than raised one at a time."""

    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors

    def raise_if_bad(self, context: str = "") -> None:
        if self.errors:
            raise SchemaError(
                f"{context or 'frame'} violates the D^i contract:\n  - "
                + "\n  - ".join(self.errors))

    def __str__(self) -> str:
        parts = []
        if self.errors:
            parts.append("ERRORS:\n  - " + "\n  - ".join(self.errors))
        if self.warnings:
            parts.append("WARNINGS:\n  - " + "\n  - ".join(self.warnings))
        return "\n".join(parts) or "OK"


def _metric_direction(name: str) -> str | None:
    """Look a metric's declared direction up in choreography.yaml."""
    cfg = yaml.safe_load(
        (_REPO_ROOT / "config" / "choreography.yaml").read_text(encoding="utf-8"))
    m = (cfg.get("metrics") or {}).get(name)
    return (m or {}).get("direction")


def validate_full(df: pd.DataFrame, approach: str) -> ValidationReport:
    """Validate an approach's FULL frame (``Di.parquet``).

    The full frame must retain rejected candidates. This is checked rather than
    trusted: reweighting in T_2's panel loop can resurrect a molecule an earlier
    tier set aside, which is impossible if the row was deleted.
    """
    r = ValidationReport()
    for col in REQUIRED_FULL_COLUMNS:
        if col not in df.columns:
            r.errors.append(f"missing required column {col!r}")
    if r.errors:
        return r

    if approach not in VALID_APPROACHES:
        r.errors.append(f"approach {approach!r} not in {sorted(VALID_APPROACHES)}")
    if (df["approach"] != approach).any():
        r.errors.append(f"rows present whose 'approach' is not {approach!r}")

    dupes = df["canonical_smiles"].duplicated().sum()
    if dupes:
        r.errors.append(f"{dupes} duplicate canonical_smiles — the join key must be unique")

    if df["canonical_smiles"].isna().any():
        r.errors.append("null canonical_smiles present; the join key cannot be null")

    # Stamp-don't-delete: a frame with no stamped rejections at all is suspicious
    # for every approach that runs gates.
    if "rejected_at" in df.columns and df["rejected_at"].notna().sum() == 0:
        r.warnings.append(
            "no rows carry 'rejected_at'. Every approach here runs at least one "
            "gate, so an entirely unstamped frame usually means rows were dropped "
            "rather than stamped.")
    return r


def validate_top10(df: pd.DataFrame, approach: str) -> ValidationReport:
    """Validate the GUI hand-off frame (``Di_top10.csv``)."""
    r = ValidationReport()

    missing = [c for c in TOP10_COLUMNS if c not in df.columns]
    if missing:
        r.errors.append(f"missing column(s): {missing}")
    leaked = [c for c in INTERNAL_ONLY_COLUMNS if c in df.columns]
    if leaked:
        r.errors.append(
            f"internal-only column(s) present: {leaked}. These are approach "
            "artifacts and must never cross the join.")
    if r.errors:
        return r

    if len(df) > 10:
        r.warnings.append(f"{len(df)} rows; the contract is a top-10 hand-off")
    if df["canonical_smiles"].duplicated().any():
        r.errors.append("duplicate canonical_smiles in the hand-off")

    bad_ap = set(df["approach"].unique()) - VALID_APPROACHES
    if bad_ap:
        r.errors.append(f"invalid approach value(s): {sorted(bad_ap)}")
    bad_mech = set(df["mechanism"].unique()) - VALID_MECHANISMS
    if bad_mech:
        r.errors.append(f"invalid mechanism value(s): {sorted(bad_mech)}")
    bad_dir = set(df["rank_metric_direction"].unique()) - VALID_DIRECTIONS
    if bad_dir:
        r.errors.append(f"invalid rank_metric_direction value(s): {sorted(bad_dir)}")
    bad_proxy = set(df["inhibition_proxy_strength"].unique()) - VALID_PROXY
    if bad_proxy:
        r.errors.append(f"invalid inhibition_proxy_strength value(s): {sorted(bad_proxy)}")

    # The declared direction must match choreography.yaml. A frame claiming
    # higher-is-better for a kcal/mol metric would invert the GUI's ranking
    # silently — this is exactly what changed under D0015.
    for name in df["rank_metric_name"].dropna().unique():
        declared = _metric_direction(str(name))
        used = set(df.loc[df["rank_metric_name"] == name, "rank_metric_direction"])
        if declared and used - {declared}:
            r.errors.append(
                f"metric {name!r} is declared {declared!r} in choreography.yaml "
                f"but the frame uses {sorted(used)}")

    # Rank must be a dense 1..N ordering consistent with the stated direction.
    ranks = df["rank"].tolist()
    if sorted(ranks) != list(range(1, len(df) + 1)):
        r.errors.append(f"rank must be 1..{len(df)} with no gaps; got {sorted(ranks)}")
    else:
        d = df.sort_values("rank")
        vals = pd.to_numeric(d["rank_metric_value"], errors="coerce")
        if vals.notna().all() and len(vals) > 1:
            direction = d["rank_metric_direction"].iloc[0]
            ordered = (vals.is_monotonic_increasing if direction == "lower_is_better"
                       else vals.is_monotonic_decreasing)
            if not ordered:
                r.errors.append(
                    f"rank order contradicts rank_metric_direction={direction!r}; "
                    "the best-ranked candidate must have the best metric value")

    for axis in ("MW", "cLogP", "TPSA", "QED", "SAscore"):
        if df[axis].isna().any():
            r.warnings.append(
                f"null values in {axis!r}; the GUI pools this axis across all four "
                "approaches and gaps will show as missing points")
    return r


def build_top10(df: pd.DataFrame, *, approach: str, mechanism: str,
                rank_metric_name: str, conditions_resolved: str,
                inhibition_proxy_strength: str, n: int = 10) -> pd.DataFrame:
    """Assemble a conforming top-10 frame from an approach's scored candidates.

    Sorts by the metric in the direction declared in choreography.yaml rather
    than a direction passed in — so an approach cannot rank the wrong way round
    by supplying the wrong argument.
    """
    direction = _metric_direction(rank_metric_name)
    if direction not in VALID_DIRECTIONS:
        raise SchemaError(
            f"metric {rank_metric_name!r} has no valid direction in "
            "choreography.yaml; declare it before ranking on it")

    work = df[df[rank_metric_name].notna()].copy()
    work = work.sort_values(rank_metric_name,
                            ascending=(direction == "lower_is_better")).head(n)
    work = work.reset_index(drop=True)

    out = pd.DataFrame({
        "candidate_id": work.get("candidate_id"),
        "canonical_smiles": work["canonical_smiles"],
        "approach": approach,
        "mechanism": mechanism,
        "rank": range(1, len(work) + 1),
        "rank_metric_name": rank_metric_name,
        "rank_metric_value": work[rank_metric_name],
        "rank_metric_direction": direction,
        "uncertainty": work.get("uncertainty"),
        "novelty_external": work.get("novelty_external"),
    })
    for axis in ("MW", "cLogP", "TPSA", "QED", "HBD", "HBA",
                 "rot_bonds", "frac_sp3", "SAscore"):
        out[axis] = work.get(axis)
    out["alert_flags"] = work.get("whole_alert_names", work.get("alert_flags"))
    out["conditions_resolved"] = conditions_resolved
    out["inhibition_proxy_strength"] = inhibition_proxy_strength
    out["pose_path"] = work.get("pose_path")
    out["manifest_ref"] = work.get("manifest_ref")
    return out[list(TOP10_COLUMNS)]
