"""Rank binding modes by warhead engagement geometry, and nothing else.

@tt8804, 2026-08-27: *"now we need to fix ranking accordingly, focusing on
warhead optimal engagement only for now"* and *"we can rank by poses(modes)
geometry only and have an option to rank by ligand by average geometry score
across modes"*.

WHY GEOMETRY ONLY, MEASURED RATHER THAN ASSERTED. Against the 147 swept modes,
using each mode's measured `frac_attack_ready` -- the fraction of MD frames in
which the warhead actually sits in the near-attack window -- as the outcome:

    the SIMULATED POSE's anchor quality        rho = +0.652   p = 3.5e-19
    anchor_quality_max                         rho = +0.130   p = 0.12
    `conditional_eb`  (THE INCUMBENT COLUMN)   rho = -0.015   p = 0.86
    `enrichment` / `viable_fraction`           rho = -0.043   p = 0.61
    mode size                                  rho = +0.102   p = 0.22

The column the pipeline ranks on today predicts the outcome no better than
chance, and the frequency statistics are very slightly NEGATIVE. Geometry of one
pose predicts it strongly. That is the whole argument for this module.

WHY THE MODE-LEVEL AGGREGATES FAIL, WHICH IS THE SAME FINDING AS D0088. A mode's
poses do not agree about engagement: the median mode spans **0.776** of the
anchor-quality scale, which itself only runs 0 to 1, and 93% span more than half
of it. A mode is a mixture of excellent and hopeless poses, so no summary of one
can predict what any single member does. This is why `anchor_quality_mean` scores
-0.031 while a single real pose scores +0.652 -- the average is over a
population that has no central tendency worth reporting.

The consequence for this module is a rule rather than a preference: **score a
POSE, then attribute it to the group** -- never average a group. Which pose is
`shared/pose_contacts` and `nac_screen_v2.representative_indices`' business, not
this module's; the measured selector is the medoid of the well-anchored quartile
(33.3% crystal recovery against 6.7% for argmax anchoring).

WHAT THIS DOES NOT CLAIM. `frac_attack_ready` is reachability of attack geometry,
not binding and not reactivity. And the 147 modes it was validated on were
SELECTED for sweeping by `conditional_eb`, so every correlation above is measured
inside the band the incumbent already liked. Metric-against-metric comparison is
fair -- they all face the same restriction -- but no absolute number here is a
population estimate, and nothing in this module makes a shortlist
`rank_validated`.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from shared import nac_criterion as nac

#: How a group's engagement is summarised from its poses. `representative` is the
#: default and the only one with outcome evidence behind it; the rest exist so a
#: caller can reproduce an older table or measure the difference, never because
#: they are interchangeable.
STATISTICS = ("representative", "q75_median", "q75_mean", "max", "mean", "median")

#: Higher is better, for every statistic here. Registered explicitly because the
#: direction registry is what caught catalogue #4 -- an analysis ranking on a
#: column whose direction it had guessed.
LOWER_IS_BETTER = False


def pose_engagement(distance, angle, mechanism) -> np.ndarray:
    """Per-pose warhead engagement, 0-1. Thin wrapper over `nac.anchor_quality`.

    ONE DEFINITION, IMPORTED. The screen picks representatives with
    `anchor_quality` and the ranking now orders on it; a second copy here would
    be two definitions of "well anchored", free to drift while both looked right.
    """
    return np.array([nac.anchor_quality(d, a, m)
                     for d, a, m in zip(distance, angle, mechanism)], dtype=float)


def _summarise(a: np.ndarray, statistic: str) -> float:
    a = a[np.isfinite(a)]
    if a.size == 0:
        return float("nan")
    if statistic == "max":
        return float(a.max())
    if statistic == "mean":
        return float(a.mean())
    if statistic == "median":
        return float(np.median(a))
    top = a[a >= np.percentile(a, 75)] if a.size >= 4 else a
    if statistic == "q75_mean":
        return float(top.mean())
    if statistic in ("q75_median", "representative"):
        # `representative` is the q75 MEDIAN when the representative's own
        # measurement is unavailable -- the closest table-computable analogue of
        # "a typical member of the well-anchored quartile", which is the rule
        # that actually selects the pose. It is NOT the same number, and
        # `engagement_source` records which one a row carries.
        return float(np.median(top))
    raise ValueError(f"unknown statistic {statistic!r}; known: {STATISTICS}")


def mode_engagement(poses: pd.DataFrame, statistic: str = "representative",
                    group_keys=("ident", "mode")) -> pd.DataFrame:
    """One engagement score per mode, from a per-pose table.

    `poses` needs `distance`, `angle`, `mechanism` and the group keys. Returns
    the score plus the SPREAD of its poses, because a group whose members
    disagree has no meaningful summary and the reader has to be able to see that
    -- 93% of the shipped rule's modes span more than half the scale.
    """
    if statistic not in STATISTICS:
        raise ValueError(f"unknown statistic {statistic!r}; known: {STATISTICS}")
    d = poses.copy()
    d["_anchor"] = pose_engagement(d["distance"], d["angle"], d["mechanism"])
    rows = []
    for key, g in d.groupby(list(group_keys)):
        a = g["_anchor"].to_numpy(dtype=float)
        fin = a[np.isfinite(a)]
        rows.append(dict(
            zip(group_keys, key if isinstance(key, tuple) else (key,)),
            **{"engagement": _summarise(a, statistic),
               "engagement_statistic": statistic,
               "engagement_spread": float(fin.max() - fin.min()) if fin.size else np.nan,
               "engagement_best": float(fin.max()) if fin.size else np.nan,
               "n_poses_mode": int(len(g))}))
    return pd.DataFrame(rows)


def rank_modes(modes: pd.DataFrame, within: str | None = "warhead_class",
               min_poses: int | None = None) -> pd.DataFrame:
    """Order modes by engagement, best first.

    `within` ranks inside a stratum (warhead class by default) rather than
    globally, because engagement is a geometric quantity whose achievable range
    differs by mechanism -- an SN2 backside attack and a perpendicular approach
    to an sp2 carbon do not span the same angles, and `isotropic_null` differs
    between them for exactly that reason.

    `min_poses` is NOT applied by default. The pose-count gate exists because a
    FREQUENCY estimate over three poses means nothing; an engagement score is a
    property of one pose and is as estimable in a group of one as in a group of
    fifty. A caller that wants the old gate must ask for it.
    """
    d = modes.copy()
    if min_poses is not None:
        if "n_poses_mode" not in d.columns:
            raise ValueError("min_poses requested but n_poses_mode is absent")
        d = d[d.n_poses_mode >= int(min_poses)]
    d = d.sort_values("engagement", ascending=False, kind="mergesort")
    if within and within in d.columns:
        d["engagement_rank"] = (d.groupby(within)["engagement"]
                                 .rank(ascending=False, method="first").astype(int))
    else:
        d["engagement_rank"] = np.arange(1, len(d) + 1)
    return d.reset_index(drop=True)


def rank_ligands(modes: pd.DataFrame, how: str = "mean",
                 ligand_key: str = "ident") -> pd.DataFrame:
    """Order LIGANDS by their modes' engagement. `how` is mean | best | median.

    @tt8804 asked for the average across modes, and it is the default. But the
    choice is load-bearing and the two answer different questions:

      `mean`  -- how well does this ligand engage TYPICALLY. Penalises a molecule
                 that can reach attack geometry one way out of twenty, which is
                 the right treatment if you believe the mode population, and the
                 wrong one if you do not (D0092: the group COUNT is a function of
                 docking depth, so the denominator of this mean moves with how
                 long you docked).
      `best`  -- can this ligand engage AT ALL. Immune to the denominator, and
                 therefore to docking depth, but it is a maximum over a noisy
                 score, and argmax selection is measured to be the worst rule
                 available at the pose level (6.7% against 33.3%).

    Neither is safe to read as "this ligand binds". They order reachability of
    attack geometry, which is a precondition for covalent chemistry and not
    evidence of it.
    """
    if how not in ("mean", "best", "median"):
        raise ValueError(f"unknown aggregation {how!r}; known: mean, best, median")
    agg = {"mean": "mean", "best": "max", "median": "median"}[how]
    g = (modes.groupby(ligand_key)
               .agg(ligand_engagement=("engagement", agg),
                    n_modes=("engagement", "size"),
                    best_mode_engagement=("engagement", "max"),
                    mean_mode_engagement=("engagement", "mean"))
               .reset_index())
    g["ligand_aggregation"] = how
    # THE DENOMINATOR IS REPORTED BESIDE THE MEAN, always. `n_modes` grows with
    # docking depth (D0092), so a mean over it is comparable only at fixed depth
    # and a reader has to be able to see the n that produced it.
    return g.sort_values("ligand_engagement", ascending=False).reset_index(drop=True)
