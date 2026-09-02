---
id: D0102
title: Running one outside molecule through the screen found three guards that each stopped one step short
date: 2026-08-31
status: accepted
approach: shared
decided_by: '@tt8804'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - scripts/dock_reference_modes.py
  - scripts/nac_screen_v2.py
  - scripts/shortlist_report.py
  - scripts/pose_modes_report.py
  - shared/reference_set.py
  - shared/ionisation.py
  - integration/pose_group_viewer.py
  - tests/test_topic_paths.py
  - tests/test_version_resolution_is_numeric.py
  - tests/test_ionisation_azole_anion.py
evidence:
  - 'issue #81, @tt8804: a SMILES for a known positive control, asked to be run "from docking up"'
  - 'DEFECT 1: dock_reference_modes set `nsv.OUT` alone; `one()` reads POSE_DIR and ALL_POSE_DIR off the module at write time, so a control run would have written its representative and its whole 500-pose cloud into the PRODUCTION topic nac_v6_poses/ and nac_v6_allposes/'
  - 'tests/test_topic_paths.py existed for exactly this defect and could not fail on it: every assertion read nac_screen_v2.py, and the caller is a different file'
  - 'DEFECT 2: `sorted(glob("warhead_classes_*.csv"))[-1]` sorts LEXICALLY, so `_9` beat `_10`. `_10` adds exactly one class, `cinnamamide` -- the class this molecule needs -- so it was unreachable'
  - 'two call sites: dock_reference_modes.py:55 (raises on a miss, so it failed loudly) and shortlist_report.py:307 (returns None on a miss, so it was silent)'
  - 'tests/test_reference_version.py could not catch it: it walks the AST for pinned version LITERALS, and there is no literal here -- the bug is in the comparison'
  - 'DEFECT 3: obabel deprotonates the tetrazole correctly (pKa ~4.9) then writes the anion as `[N-]1=C(NN=N1)R`, three bonds on an anionic N, which RDKit rejects. `ionisation.protonate` guaranteed IDENTITY and never checked VALIDITY'
  - 'obabel only mis-writes an anion it CREATES: handed the aromatic `c1nnn[n-]1` it round-trips it perfectly, which is why every neutral molecule screened to date passed'
  - 'RESULT after the fixes: 500 poses, 424 PoseBusters-valid (84.8%), 94 modes, best viable_fraction 0.200, best engagement 0.596'
  - 'ON THE CONFIGURED RANKING (ranking.ligand_agg = fraction_above, cut anchor_quality_max >= 0.5, score rank_score): control frac_above = 0.1596 (15 of 94 modes), mean score 0.1203, rank 213 of 1,685 -- top 12.6%, 7.9x better than random'
  - 'library frac_above median 0.0849, mean 0.0942, max 0.3197; the control is 1.88x the median'
  - 'THE AGGREGATIONS DISAGREE ABOUT THIS MOLECULE: fraction_above 213/1685 (top 12.6%) against best 1,645/1,685 (bottom 2.4%), from one pose cloud'
  - 'CONFOUND CHECKED: frac_above is negatively correlated with mode count in the library (rho = -0.197, p = 3.4e-16) and the control has 94 modes against a library median of 192, so the denominator flatters it. Against the 37 molecules within 75-125% of its mode count it still beats 86.5%, against 87.4% overall -- the placement is not an artefact of the denominator'
  - 'RECOVERY: a top-25, top-50 or top-100 shortlist MISSES it; top-250 finds it'
  - '297 of 500 poses put the BETA-CARBON in the 2.8-4.2 A window, but only 14 of 500 (2.8%) are `viable` -- the off-normal angle rejects them, median 58.1 deg against a 45 deg bar'
  - '@tt8804, 2026-08-31, settling the mechanism: "the alpha,beta-unsaturated ketone (enone) is the warhead and is supposed to target the cys113"'
runbook: null
---

## Context

Issue #81 supplied one SMILES for a known positive control and asked for it to
be run through the pipeline from docking up. The molecule is
`O=C(C=Cc1ccccc1CCc1nnn[nH]1)n1c(-c2ccccc2)cc2ccc(Br)cc21` — a cinnamoyl amide
on a 5-bromo-2-phenylindole, with a tetrazole on the other arm. It is in none of
the reference binder sets, so nothing in this repo had seen it.

`scripts/dock_reference_modes.py` already exists for exactly this — "would our
screen have found it" — so the expectation was that this was a one-command job.
It was not. Three separate guards let it through, and **all three are the same
shape**: a check that exists, runs, and is scoped one step short of the case in
front of it. That is disguise #4 in
[`how_this_project_breaks.md`](../docs/how_this_project_breaks.md), and it is
still the fastest-growing group.

## Decision

Fix all three as classes, and move each guard out to the scope that would have
caught it.

**1. Output paths move together, or not at all.** `nac_screen_v2.one()` reads
`POSE_DIR` and `ALL_POSE_DIR` off the module when it writes, so any caller
driving `one()` directly has to move them itself. `dock_reference_modes` set
`OUT` alone. Tables would have gone to `ian_ctrl/` while a molecule that **is
not in the library** put its representative and its 500-pose cloud into
`nac_v6_poses/` and `nac_v6_allposes/` — where anything globbing the production
cloud picks it up as a candidate.

The replacement is `nac_screen_v2.use_topic(topic)`, which rebinds all four
globals and is the only supported way to move them. The test now bans assigning
**any** of the four individually, anywhere in the repo — not "check all four are
assigned", because that is a condition a later edit deletes one line from and
nothing complains.

**2. Version resolution compares integers.** `latest_reference()` exists so
nobody hand-pins a version, and two scripts rolled their own with
`sorted(glob(...))[-1]`. That *looks* like dynamic resolution and is — to the
wrong file, for as long as the version count has two digits.
`shared/reference_set.warhead_library()` / `load_warhead_row()` are now the one
resolver, and `load_warhead_row` **raises** on an unregistered class rather than
returning a default, following `canonical_class()`.

**3. `protonate()` validates what it returns.** It guaranteed the right string
reached the right id — with recursive-split machinery above it built precisely
so a short return could not slide results onto the wrong candidate — and said
nothing about whether the string was a molecule. It now parses obabel's output,
repairs this one well-defined serialisation fault, and **drops** anything it
cannot repair, keeping the existing contract that a missing id means no species
could be built. The repair neutralises the over-valent nitrogen, lets RDKit
perceive the aromatic ring, then removes the ring NH again: obabel decides
*which* site ionises, RDKit decides what a valid structure looks like, and the
charge obabel chose is re-asserted rather than recomputed.

## Why the wrong thing looked right

**Defect 1** looked right because the fix for this exact bug was already in the
file, three functions up, with a comment explaining it: *"Tables and poses are
the same claim about the same run; a code path that can separate them will
eventually separate them."* `main()` does it correctly. The test file
`test_topic_paths.py` is entirely about it. Everything was in place except that
every assertion opened `nac_screen_v2.py` and the offending line is in a
different file. **A guard scoped to one file cannot see a caller.**

**Defect 2** looked right because it is not a pin. The project's whole defence
against stale versions — `test_reference_version.py` walking the AST over every
directory — hunts for version *literals*, and there is no literal here. The line
reads as the correct idiom. It also could not fail until the library passed nine
versions, so it was correct for its first nine and silently wrong afterwards,
with no edit in between to blame. And its two call sites failed differently:
`dock_reference_modes` raises on an unknown class, so it stopped dead with a
clear message; `shortlist_report` returns `None`, so there an aryl Michael
acceptor simply had no reactive atom marked and nothing said why.

**Defect 3** looked right because obabel is right about the chemistry. At pH 7.4
the tetrazole *is* deprotonated, the charge is correct, and the string is a
perfectly plausible SMILES. It is wrong only in valence, and only for an anion
obabel itself created — the same converter round-trips the aromatic form
faultlessly. Every molecule screened to this point was neutral at every azole,
so the first ionisable heteroaryl through the door found it. The failure then
surfaced two stages downstream in `prepare_ligand` as "unparseable pH 7.4
SMILES for ian_ctrl_issue81", which reads like a property of the molecule
rather than of the converter.

## Consequences

The control ran: 500 poses, 424 PoseBusters-valid, 94 modes, into topic
`ian_ctrl`, at the library's own settings (500 runs, `contact_linkage`, seed 42)
so its numbers are comparable to nac_v6's.

**CORRECTION, same day.** This section first reported the control as ranking
above only 3.0% of the library. That number was computed on the per-molecule
**maximum** engagement — a metric nobody ranks on. `config/target.yaml` declares
`ranking.ligand_agg: fraction_above`, and on the metric that actually orders the
library the control is **213rd of 1,685, top 12.6%**. The error is this
project's own house defect committed while documenting it: a value selected
because it was plausible and adjacent, rather than because the config named it.

The four aggregations do not merely differ in degree — they disagree about
whether this molecule is interesting at all, from one pose cloud:

| aggregation | control | rank | |
|---|---:|---:|---|
| **`fraction_above`** (configured) | 0.1596 | **213 / 1,685** | top 12.6% |
| `mean` | 0.1203 | 693 / 1,685 | top 41.1% |
| `median` | 0.0000 | 1,149 / 1,685 | top 68.2% |
| `best` | 0.5962 | 1,645 / 1,685 | **bottom 2.4%** |

`best` and `fraction_above` place the same molecule at opposite ends. That is
not a tie-break-sized disagreement, and it means **the aggregation choice, not
the geometry, decides this molecule's fate.** `rank_ligands`' own docstring says
mean and best "answer different questions"; this is the first measurement of how
far apart the answers can be on a molecule someone cares about.

**SETTLED, 2026-08-31.** @tt8804: *"the α,β-unsaturated ketone (enone) is the
warhead and is supposed to target the cys113."* The reactive dock was therefore
the right measurement and the ranking result stands. (The question was worth
asking: the tetrazole is a phosphate mimic for the Arg68/Arg69 basic triad and
is the signature of a *non-covalent* Pin1 inhibitor, and scoring a reversible
binder on a near-attack criterion measures nothing — that is
`elevate_reference.py`'s stated reason for plain-docking ATRA.)

**AND THE SCREEN FAILS IT ON DIRECTION, NOT DISTANCE.** Measured over the
persisted cloud, against Cys113 SG located by residue identity:

| | |
|---|---:|
| poses with the β-carbon in the 2.8–4.2 Å window | **297 / 500 (59.4%)** |
| poses `viable` (distance **and** ≤45° off-normal) | **14 / 500 (2.8%)** |
| best-energy pose (−10.99 kcal/mol), β-C···SG | **3.21 Å** |
| in-window among the best 25% by energy | 82.4% |
| in-window among the rest | 51.7% |

So the warhead gets to the sulfur at the right separation in most poses, and
energy and distance *agree* about it — the best-scoring quartile is in the
window 82.4% of the time against 51.7% for the rest. What fails is the approach
**direction**: median 58.1° off the sp2 normal where the criterion wants ≤45°.

Two readings, and they are not distinguishable from this run alone: either the
molecule genuinely cannot present its enone face-on in this pocket, or the
rigid-receptor docking cannot produce that geometry. The distance half is also
partly manufactured — the reactive potential is a *sampler* biased toward
warhead–sulfur contact (D0064), which is exactly why `conditional_eb` conditions
on distance. **The 2.8% viable rate is the number to carry forward, not the
59.4%.**

Also worth recording: the four warhead classes whose reactive SMARTS match this
molecule (`naphthoquinone_c2`, `naphthoquinone_benzo`, `acrylamide`,
`cinnamamide`) all share the identical pattern `[CX3]=[CX3][CX3]=O` and the
identical mechanism, so the class label changes nothing that is measured. It was
chosen as `cinnamamide` on chemistry, not by taking the first match.

`integration/pose_group_viewer.py` now lists **every topic's** clouds, not only
`run.topic`'s, so viewing a targeted run no longer requires bumping `run.topic`
— which is global state that detached supervisors poll. It also highlights a
named class's reactive atoms and draws the anchor they must reach, with the
anchor located by residue **identity**: Pin1 has two cysteines and Cys57 is the
one a positional shortcut finds.
