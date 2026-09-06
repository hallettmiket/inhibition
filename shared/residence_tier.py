#!/usr/bin/env python3
"""
Purpose: the three-tier verdict for a 100 ns production run.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-17

@tt8804: "so there should be held in pocket, held but not optimal and below
max .35 is optimal" -- and "red green yellow".

WHY THREE TIERS AND NOT A BAR. Two thresholds were already in play and were
being read as rivals: 1.0 nm (`mdprio_report.BOUND_NM`, the pre-registered
residence criterion -- did the ligand ever leave and not come back) and 0.35 nm
(`tc.md_survivor_rmsd_nm()`, the triage-sweep survivor bar). Scoring a
production run against 0.35 calls a run that never left the pocket a failure;
scoring it against 1.0 alone cannot distinguish a ligand pinned at the warhead
from one rattling around the site. They are not rivals -- 0.35 is the tighter
tier inside 1.0, so a run gets one of three verdicts, not a pass/fail.

  optimal   never dissociated AND max ligand RMSD < 0.35 nm
  held      never dissociated, but wandered past 0.35 nm
  left      dissociated -- a frame after which it never came back

COLOUR IS THE SECOND SIGNAL, NEVER THE ONLY ONE. report_theme's own note says
these tables are read by people with red-green deficiency, so every tier
carries its label and colour only reinforces it.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shared import target_config as tc          # noqa: E402

# The residence criterion: the ligand is "in the pocket" at RMSD <= this.
BOUND_NM = 1.0

#: (key, label, theme colour token)
#:
#: FOUR TIERS (@tt8804: "it should be yellow held vs green held and green
#: optimal"). "Never dissociated" was doing too much work: a run that travels
#: 5.765 nm, comes back under 1.0 and stays counts as held under the
#: pre-registered rule, and t4_cc678a20a3d0_m2 did exactly that -- sitting in
#: the held band with 23% of its frames outside the pocket, next to runs that
#: never left at all. Both are "held"; only one is a candidate. The residence
#: fraction separates them, so held splits by how much of the run was actually
#: spent bound.
#:
#: The two green tiers are deliberately both green -- optimal is a tighter
#: reading of the same good outcome -- and every tier still carries a distinct
#: WORD, because report_theme's own note is that these tables are read by people
#: with red-green deficiency and colour can never be the only signal.
#: FIVE TIERS since 2026-09-06. `inert` was added because two 100 ns runs held
#: the pocket beautifully and never once put their warhead in reach of Cys113:
#: `t4_f431658edd9a_m37` scored the best residence in the project (mean ligand
#: RMSD 0.124 nm, 6.7 contacts, n_independent 462 -- genuinely still, not slowly
#: drifting) with the warhead engaged **0.0% of 100 ns**, parked at 5.65 A and
#: 66 deg off-normal from the first quarter to the last. It rendered as
#: **"Optimal"**.
#:
#: That is the project's signature defect in a label: a verdict about the LIGAND
#: read as a verdict about the CANDIDATE. In a covalent screen "optimal" has to
#: mean optimal for the reaction, not for sitting still. `inert` is yellow, not
#: green, because a molecule that cannot reach the cysteine is not a hit however
#: still it sits -- and it keeps its own WORD, since colour is never the only
#: signal here.
TIERS = (
    ("optimal", "optimal", "good"),                    # green
    ("held", "held", "good"),                          # green
    ("inert", "held, warhead never engaged", "warn"),  # yellow
    ("unstable", "held, unstable", "warn"),            # yellow
    ("left", "left", "bad"),                           # red
)
_LABEL = {k: (lab, col) for k, lab, col in TIERS}


def optimal_nm() -> float:
    """The tight tier, from config.

    NOT the sweep survivor bar. This was `md_survivor_rmsd_nm()` (0.35 nm), the
    number that decides which 8 ns poses earn a 100 ns run -- and reusing it
    here made the tier unreachable by construction rather than by chemistry: a
    ligand explores more in 100 ns than in 8, and 13 finished runs bottomed out
    at 0.410 nm, leaving the tier empty. @tt8804 set the production bar to
    0.45 nm; 0.35 remains the sweep bar and is unchanged.
    """
    return float(tc.md_production_optimal_rmsd_nm())


def residence_floor() -> float:
    """Fraction of frames a run must spend bound to count as a clean hold.

    Below it the run is `unstable`: it came back, so it did not dissociate, but
    it was out of the pocket for a real part of the trajectory.
    """
    return float(tc.md_held_residence_floor())


#: A warhead engaged in less than this fraction of the run never got into
#: position at all. Deliberately near zero: this separates "never" from "rarely",
#: not "often" from "sometimes", and the two runs that motivated it sat at
#: exactly 0.000 and 0.080.
ENGAGED_FLOOR = 0.01


def tier(rmsd_max_nm: float | None, dissociated: bool | None,
         residence_frac: float | None = None,
         engaged_frac: float | None = None) -> str:
    """Return the tier key for one run.

    Raises on unusable input rather than defaulting. A missing RMSD is "not
    measured", and this project has been bitten repeatedly by an unmeasured
    value entering a table as a passing one. `residence_frac` is required to
    separate held from unstable, and omitting it raises for the same reason --
    silently treating an unknown residence as a clean hold is precisely the
    "value taken by label not identity" failure.
    """
    if dissociated is None or rmsd_max_nm is None:
        raise ValueError("residence tier needs both rmsd_max_nm and dissociated; "
                         f"got {rmsd_max_nm!r} and {dissociated!r}")
    if dissociated:
        return "left"
    if residence_frac is None:
        raise ValueError("residence tier needs residence_frac to tell a clean "
                         "hold from an excursion that came back")
    if float(residence_frac) < residence_floor():
        return "unstable"
    # THE WARHEAD DECIDES WHETHER A CLEAN HOLD IS A CANDIDATE.
    #
    # Optional rather than required, because this function is also called on
    # runs where the geometry series is unavailable and raising there would lose
    # a real residence verdict. But WHEN IT IS KNOWN it outranks the RMSD: a
    # molecule pinned 5.7 A from Cys113 for 100 ns is not "optimal", whatever
    # its ligand RMSD says.
    if engaged_frac is not None and float(engaged_frac) < ENGAGED_FLOOR:
        return "inert"
    return "optimal" if float(rmsd_max_nm) < optimal_nm() else "held"


def label(key: str) -> str:
    return _LABEL[key][0]


def colour(key: str) -> str:
    """Theme token -- 'good' | 'warn' | 'bad', resolved by report_theme."""
    return _LABEL[key][1]


def badge(key: str) -> str:
    """Inline HTML badge: colour AND text, never colour alone."""
    lab, col = _LABEL[key]
    return (f'<span class="tier tier-{key}" '
            f'style="color:var(--{col});border-color:var(--{col})">{lab}</span>')


TIER_CSS = """
.tier{display:inline-block;padding:.05rem .45rem;border:1px solid;border-radius:.7rem;
      font-size:.78rem;font-weight:600;letter-spacing:.01em;white-space:nowrap}
"""
