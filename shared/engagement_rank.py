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
                 ligand_key: str = "ident", cutoff: float | None = None,
                 score_col: str = "engagement",
                 cut_col: str | None = None) -> pd.DataFrame:
    """Order LIGANDS by their modes' engagement.

    `how` is mean | best | median | fraction_above.

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
    if how not in ("mean", "best", "median", "fraction_above"):
        raise ValueError(f"unknown aggregation {how!r}; known: mean, best, "
                         "median, fraction_above")
    if how == "fraction_above":
        # @tt8804: "rank by mol based either on the proportion of a rankscore
        # above a cutoff, or a molecule mode mean rankscore."
        #
        # WHAT IT ASKS: of everything this molecule can do, how much of it is
        # reaction-competent? A molecule with 200 groups of which 40 clear the
        # cutoff is a better bet than one with 200 of which 2 do, even if both
        # have the same best group -- which `best` cannot see and `mean` blurs,
        # because the median group scores ~0 and drags every mean toward it.
        #
        # STILL DEPTH-DEPENDENT, and it has to be said: both numerator and
        # denominator grow with docking depth, and the groups added at depth are
        # mostly sparse singletons (D0092: 61% of groups hold <= 3 poses, and the
        # count never saturates). So a fraction falls as you dock longer. It is
        # comparable across molecules AT FIXED DEPTH and nowhere else; `best` is
        # the only depth-immune aggregation of the four.
        if cutoff is None:
            raise ValueError("fraction_above needs an explicit cutoff -- there is "
                             "no defensible default, and a silent one would set "
                             "the selection rule without anybody choosing it")
        # THE CUTOFF MAY BE APPLIED TO A DIFFERENT COLUMN THAN THE ORDERING
        # (`cut_col`). The threshold is a GEOMETRIC criterion -- anchor >= 0.5 is
        # exactly "inside the NAC window", because the score is a product of two
        # [0,1] factors and so >= 0.5 forces both >= 0.5, which are precisely
        # 2.8-4.2 A and <= 30 degrees. `rank_score` multiplies that anchor by a
        # support bonus of up to 1.15, so thresholding IT at 0.5 admits poses at
        # anchor 0.44 -- outside the window, promoted by a hydrogen bond. 1,448
        # modes on nac_v6 do exactly that. The bonus belongs in the ORDER, never
        # in the criterion.
        cut_on = cut_col or score_col
        if cut_on not in modes.columns:
            raise KeyError(
                f"fraction_above needs {cut_on!r} to apply its cutoff and the "
                f"frame does not carry it. Falling back to {score_col!r} would "
                f"silently change what the threshold MEANS, which is the whole "
                f"defect this argument exists to prevent.")
        g = (modes.assign(_hit=modes[cut_on] >= float(cutoff))
                  .groupby(ligand_key)
                  .agg(ligand_engagement=("_hit", "mean"),
                       n_modes=("_hit", "size"),
                       n_above=("_hit", "sum"),
                       best_mode_engagement=(score_col, "max"),
                       mean_mode_engagement=(score_col, "mean"))
                  .reset_index())
        g["ligand_aggregation"] = f"fraction_above({cut_on} >= {cutoff})"
        # MEAN BREAKS THE TIE (@tt8804, 2026-08-29). A fraction over a few
        # hundred modes is coarse -- two molecules landing on the same value is
        # a real event, not a rounding artefact -- and leaving the tie to the
        # frame's row order is selection by position, which is the shape this
        # project breaks on (`how_this_project_breaks.md`, disguise #2). A
        # secondary key makes the order TOTAL and reproducible, and `mean` is
        # the right one: it is the same quantity read at finer resolution, so
        # the tiebreak never contradicts the primary key's logic.
        #
        # MEASURED on nac_v6: at cutoff 0.4 the fraction takes 1,288 distinct
        # values over 1,684 molecules and the top 50 hold 49 of them, so the
        # tiebreak decides roughly one placement in fifty. It is cheap
        # insurance, not a second ranking.
        return (g.sort_values(["ligand_engagement", "mean_mode_engagement"],
                              ascending=False)
                 .reset_index(drop=True))
    agg = {"mean": "mean", "best": "max", "median": "median"}[how]
    g = (modes.groupby(ligand_key)
               .agg(ligand_engagement=(score_col, agg),
                    n_modes=(score_col, "size"),
                    best_mode_engagement=(score_col, "max"),
                    mean_mode_engagement=(score_col, "mean"))
               .reset_index())
    g["ligand_aggregation"] = how
    # THE DENOMINATOR IS REPORTED BESIDE THE MEAN, always. `n_modes` grows with
    # docking depth (D0092), so a mean over it is comparable only at fixed depth
    # and a reader has to be able to see the n that produced it.
    return g.sort_values("ligand_engagement", ascending=False).reset_index(drop=True)


def summarise_anchor(anchor, statistic: str | None = None) -> float:
    """One engagement number from an array of per-pose anchor qualities.

    The screen's entry point: it already holds `anchor` for every pose of a
    group, so it needs the summary and not the whole `mode_engagement` frame.
    Shares `_summarise` with that function rather than reimplementing it -- two
    definitions of "this group's engagement" is exactly the drift this project
    keeps paying for.
    """
    import numpy as _np
    a = _np.asarray(anchor, dtype=float)
    return _summarise(a, statistic or "representative")


def anchor_spread(anchor) -> float:
    """Range of per-pose engagement inside one group, on the same 0-1 scale.

    Reported beside every engagement score. A group spanning most of the scale
    is a mixture and its summary means nothing; 93% of the shipped splitter's
    modes were in that state.
    """
    import numpy as _np
    a = _np.asarray(anchor, dtype=float)
    a = a[_np.isfinite(a)]
    return float(a.max() - a.min()) if a.size else float("nan")


# --------------------------------------------------------------------------- #
#  Which metrics are frequencies, and therefore need a population gate
# --------------------------------------------------------------------------- #
#: A FREQUENCY metric estimates "how often does this group reach attack
#: geometry", so it is unestimable on a handful of poses and needs a population
#: floor. A POSE-PROPERTY metric asks "how good is the geometry", which one pose
#: answers as well as fifty.
#:
#: AN ALLOWLIST OF FREQUENCIES, not a denylist of the others. A metric nobody
#: registered is treated as a frequency and KEEPS the gate, because the failure
#: of guessing wrong in that direction is a smaller shortlist, while guessing
#: wrong the other way silently publishes ranks computed from three poses.
#: #14 is the record of a denylist admitting a value nobody anticipated.
FREQUENCY_METRICS = frozenset({
    "enrichment", "enrichment_joint", "enrichment_conditional",
    "conditional_eb", "conditional_lcb", "conditional_x_consensus",
    "viable_fraction", "consensus", "weighted_score",
})

#: Metrics that are a property of geometry rather than of a count.
POSE_PROPERTY_METRICS = frozenset({
    "engagement", "engagement_best", "anchor_quality_max",
    "anchor_quality_mean", "anchor_quality_p90",
})


def needs_population_gate(metric: str) -> bool:
    """Does ranking on `metric` require the pose-count gate?

    The gate is a property of the METRIC, not of the pipeline. Applying a
    frequency's estimability threshold to a geometry score is what would have
    made the nac_v6 re-screen useless: under contact grouping only 2.8% of
    groups hold 12 poses, against 66.7% under the rule the floor was calibrated
    on, so a 12-pose gate would have discarded 97% of the run's own output while
    reporting a full-looking shortlist.
    """
    if metric in POSE_PROPERTY_METRICS:
        return False
    return True


# --------------------------------------------------------------------------- #
#  Supporting contacts — conservative, bounded, and multiplicative
# --------------------------------------------------------------------------- #
#: The residues whose side chains may earn a pose credit, and the atom on each
#: that does it. DELIBERATELY SHORT, and the omissions are the design.
#:
#: @tt8804: "we do not want to influence the physical nature of molecules we
#: select, only adding basic derivative/supporting rules that would be essential
#: for cys113 engagement."
#:
#: Measured on the prepared 3IKD, six residues sit inside the 4.2 A near-attack
#: window: CYS113, HIS59, SER114, SER115, LEU122, LEU61. Only the two serines are
#: here, because only they can hold the warhead at the reaction centre without
#: expressing a preference about what KIND of molecule wins:
#:
#:   SER114 (3.3 A) / SER115 (3.4 A)  -- INCLUDED. Both donor and acceptor, both
#:       flanking the sulfur. A polar contact here holds the warhead in place at
#:       the reaction centre; it does not favour any scaffold class.
#:   HIS59 (3.2 A)  -- EXCLUDED. Rewarding contact with it rewards pi-stacking,
#:       which is a shape preference and would change what gets selected.
#:   LEU122 / LEU61 -- EXCLUDED. Rewarding hydrophobic burial pushes selection
#:       toward greasier molecules. That is precisely the influence being ruled
#:       out.
#:   The Arg loop (LYS63, ARG68, ARG69, 7-10 A) -- EXCLUDED. It is where Pin1's
#:       own phosphopeptide affinity lives, and it is the least conservative
#:       thing that could be added: it would select for anionic ligands.
SUPPORT_ATOMS = (("A", 114, "SER", "OG"), ("A", 115, "SER", "OG"))

#: An H-bond, generously: heavy-atom donor/acceptor separation. Wider than the
#: 2.6-3.2 A ideal because the pose is a docked prediction, not a crystal.
HBOND_MIN_A = 2.4
HBOND_MAX_A = 3.6

#: The most a supporting contact may add. A TIE-BREAKER, not a weight -- at 0.15
#: a fully supported pose outranks an unsupported one of equal geometry, and
#: cannot outrank a pose with meaningfully better attack geometry.
MAX_SUPPORT = 0.15


def support_factor(ligand_xyz, ligand_elements, support_xyz,
                   max_support: float = MAX_SUPPORT) -> float:
    """1.0 to 1+max_support: does the pose hold itself at the reaction centre?

    BONUS-ONLY AND MULTIPLICATIVE, and both properties are load-bearing:

      * bonus-only -- an unsupported pose keeps its geometry score untouched, so
        adding this rule can never DEMOTE a pose that reaches attack geometry.
      * multiplicative -- `rank_score = anchor_quality * support`, so a pose with
        anchor_quality 0 scores 0 however well supported it is. Support cannot
        promote a pose that is not reaction-competent. An additive term could,
        and that is the whole risk being avoided.

    Only N and O ligand atoms count: a polar contact is the claim, and letting
    carbon satisfy it would make this a proximity bonus, which is a shape
    preference by another name.
    """
    import numpy as _np
    xyz = _np.asarray(ligand_xyz, dtype=float)
    els = [str(e).upper() for e in ligand_elements]
    polar = _np.array([i for i, e in enumerate(els) if e in ("N", "O")], dtype=int)
    if polar.size == 0 or len(support_xyz) == 0:
        return 1.0
    P = xyz[polar]
    hits = 0
    for s in support_xyz:
        d = _np.sqrt(((P - _np.asarray(s, dtype=float)) ** 2).sum(-1))
        if ((d >= HBOND_MIN_A) & (d <= HBOND_MAX_A)).any():
            hits += 1
    return 1.0 + max_support * (hits / len(support_xyz))


def receptor_support_atoms(receptor_pdb=None) -> list:
    """Coordinates of `SUPPORT_ATOMS`, matched by (chain, resi, resname, atom).

    RAISES on a missing one rather than scoring without it. A support term that
    silently loses one of its two serines would return a systematically lower
    factor for every pose and look exactly like a chemistry result.
    """
    import numpy as _np
    from shared import run_paths as _rp
    path = receptor_pdb or _rp.receptor_prep()
    want = {(c, i, rn, an): None for c, i, rn, an in SUPPORT_ATOMS}
    for ln in open(path, errors="replace"):
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        k = ((ln[21:22].strip() or "A"), int(ln[22:26]), ln[17:20].strip(),
             ln[12:16].strip())
        if k in want:
            want[k] = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
    missing = [k for k, v in want.items() if v is None]
    if missing:
        raise ValueError(f"support atoms absent from {path}: {missing}")
    return [_np.array(v) for v in want.values()]


def rank_score(anchor: float, support: float = 1.0) -> float:
    """The per-mode rank score. A mode and a single pose are the same thing here
    -- a pose is a mode of one -- so this is the only score either needs.
    """
    a = float(anchor)
    if a != a:
        return float("nan")
    return a * float(support)
