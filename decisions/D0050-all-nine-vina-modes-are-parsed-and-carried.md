---
id: D0050
title: All nine Vina modes are parsed and carried as descriptive labels, not as a confidence score
date: 2026-08-01
status: accepted
approach: shared
decided_by: '@mhallet'
origin: user
supersedes: []
superseded_by: null
affects:
  - shared/noncovalent_dock_run.py
  - tests/test_vina_modes.py
evidence:
  - '400 completed du_xu ligands re-parsed: every one reports exactly 9 modes — the cap, always saturated'
  - 'median affinity gap between mode 2 and mode 1: 0.100 kcal/mol'
  - '98% of ligands have mode 2 within 0.5 kcal/mol of mode 1'
  - "Vina's own reported RMSE is ~2-3 kcal/mol — 20-30x the margin mode 1 is selected by"
  - 'median nearest-neighbour RMSD from another mode back to the best: 1.176 nm, piling up at the ~1 A diversification floor'
  - 'best-affinity values are byte-identical to the previous single-value parse across all 400'
  - 'D0046 headroom: top-1 pose accuracy 22.5% vs best-of-9 55.0%; median RMSD 4.04 A vs 1.80 A'
---

# Stop discarding eight ninths of what we compute

## Context

`collect_scores` regex-searched for the **first** `REMARK VINA RESULT` line and
returned one float per ligand. Vina writes every mode it reports, each with its
own affinity and its RMSD back to the best mode, so 8 of every 9 poses we paid
to compute were written to disk and never read — across all ~56,000 candidates.

Issue #10 (@tt8804): *"each pose is very meaningless on its own and we want to
look for multiple similar poses"*. Measured on our own output, that is right in
a stronger sense than asserted. The margin by which mode 1 beats mode 2 is a
**median 0.100 kcal/mol**, and 98% of ligands have mode 2 within 0.5 kcal/mol,
against Vina's own reported RMSE of ~2–3 kcal/mol. **We select one pose over
another on a margin twenty to thirty times smaller than the error bar of the
function producing it.** Picking mode 1 is, statistically, picking arbitrarily.

## Decision

Parse **all** modes. `collect_modes` returns one row per candidate carrying
`vina_n_modes`, `vina_mode2_gap`, `vina_mode_rmsd_nn` and
`vina_affinity_spread` alongside `vina_affinity`. `collect_scores` is retained
as the single-value accessor, implemented in terms of the full parse so the two
cannot drift.

**These are descriptive labels. Nothing ranks on them.**

## The trap this deliberately does not walk into

**Vina's 9 modes are not independent samples.** Vina diversifies its reported
modes with a minimum-RMSD floor before writing them out, and our own data shows
it: the median nearest-neighbour RMSD back to the best mode is **1.176**,
piling up right at the ~1 A floor.

So a cluster population computed over these nine would measure **the output
diversification filter, not pose convergence**. It would produce a number, it
would look like a consensus statistic, and it would be an artefact of the
formatter — precisely the shape catalogued in `how_this_project_breaks.md`: a
populated, plausible value computed from the wrong thing. The honest version is
replicate runs with independent seeds, which Vina-GPU supports with no code
change since it draws a fresh seed per invocation.

`vina_mode_rmsd_nn` therefore excludes mode 1 (whose RMSD to itself is 0.0 by
construction and would read as perfect convergence for every ligand ever
docked), and reports **NaN rather than 0.0** for a single-mode ligand. Both are
pinned by tests.

## Consequences

No frame changes value: best-affinity is byte-identical to the previous parse
on all 400 ligands checked, so nothing downstream moves. The new columns are
added by a merge whose drop list is **derived from the merge** rather than
hand-maintained — catalogue entry #5 — with a post-merge assertion that fails
loudly on any `_x`/`_y`.

Frames docked before this landed lack the columns. Nothing needs re-docking:
the poses are on disk, so `collect_modes` can backfill them.

**The headroom is large and measurable.** D0046 puts top-1 pose accuracy at
22.5% against **best-of-9 at 55.0%** (median RMSD 4.04 A → 1.80 A), so a
selection rule that chose well among the nine could more than double pose
accuracy. That is the prize; this decision only stops throwing away the data
needed to chase it.

**The rule that gates the next step:** no consensus metric enters the pipeline
until it has been scored on D0046's 80 redock cases against crystal ground
truth. If it does not beat 22.5% there, it does not get to label anything. And
consensus measures precision, not accuracy — a systematically wrong scoring
function will be reproducibly, consensually wrong.
