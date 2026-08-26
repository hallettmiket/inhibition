---
id: D0095
title: Contact groups reproduce across independent dockings, and D0093's filter does not invalidate exp/16
date: 2026-08-26
status: proposed
approach: shared
decided_by: '@tt8804'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - exp/16_contact_clustering/run_all.py
  - exp/19_contact_reproducibility/run_all.py
  - scripts/persist_raw_clouds.py
  - decisions/D0093-the-file-named-allposes-is-not-all-poses-it-is-dbscan-cleaned.md
evidence:
  - 'REPRODUCIBILITY, 5 independent 500-run dockings of t4_716800c125a7 at tolerance 0.73 A: groups of >=5 poses match pairwise at 98.5% (range 91-100%) and 31 of 34 (91%) are present in ALL FIVE'
  - 'including singletons the rates fall to 79% pairwise and 59% in all five -- the disagreement is the singleton tail, which is what a singleton is'
  - 'group counts per replicate 112, 112, 115, 118, 109 -- stable to +-4%'
  - 'CONTRAST: D0088/#78 measured HDBSCAN keeping only 1 of 3 modes across an independent draw, and losing the MD-validated pose in 3 of 30 replicates'
  - 'THE DBSCAN BASELINE WAS REFUSED, NOT REPORTED: it returned ONE mode in every replicate and therefore reproduced at 100%. A rule that always answers "one mode" is perfectly reproducible and discriminates nothing; it is also circular on these clouds, which were already DBSCAN-cleaned'
  - 'RAW vs FILTERED, matched on the same 20 molecules: 10,000 raw poses against 8,376 filtered (16% removed), groups per molecule 217 vs 174, median within-group Cartesian RMSD 1.12 A vs 1.08 A, 90th percentile 1.63 A vs 1.62 A'
  - 'the worst single group degrades on raw: 3.76 A against 2.71 A -- above the 3.5 A sweep bar, and the one place the filtered clouds flattered the method'
  - 'a matched pair confirms D0093 directly: t4_9e43c2c37bce docks 500 poses raw and its production cloud holds 433 (13.4% absent)'
  - 'scripts/persist_raw_clouds.py now EXITS NON-ZERO when nothing is persisted; its first run failed all 20 molecules on `No module named gemmi` (wrong environment), printed "0 ok, 20 failed" and returned success, and the chained step after it ran against nothing and also reported success'
runbook: python exp/19_contact_reproducibility/run_all.py; ~/.micromamba/envs/dwi_reactive/bin/python scripts/persist_raw_clouds.py --n-molecules 20 --gpu 2; python exp/16_contact_clustering/run_all.py --raw
---

# D0095 — the two tests exp/16 had not passed

D0092 closed with two open items: contact groups had never been checked against
an **independent docking**, and every measurement had been made on clouds the
shipped DBSCAN had already cleaned (D0093). Both are now run.

## 1. Reproducibility across independent dockings

Five independent 500-run dockings, grouped at the same tolerance, groups matched
by the same rule that defines membership:

| | pairwise | present in all five |
|---|---|---|
| groups of ≥ 5 poses | **98.5%** (91–100%) | **31 of 34 — 91%** |
| all groups incl. singletons | 79.0% | 66 of 112 — 59% |

Group counts per replicate: 112, 112, 115, 118, 109 — stable to ±4%.

The singleton tail is where the disagreement lives, which is what a singleton is:
one pose that happened to land somewhere, with no claim to being a mode. The
groups that carry population reproduce.

**The contrast is the point.** HDBSCAN kept **1 of 3** modes across an
independent draw and lost the MD-validated pose in 3 of 30 replicates (D0088,
#78). exp/17 showed contact groups survive *deeper sampling of one cloud*, which
is a weak claim — a fixed tolerance carves fixed regions, so a deeper draw cannot
move them. An independent docking can, and does not.

### The baseline was refused rather than reported

A DBSCAN baseline was run on the same five clouds and returned **one mode in every
replicate**, reproducing at 100%. That number is not a comparison: a rule that
always answers "one mode" is perfectly reproducible and useless. It is also
circular — these clouds have already been DBSCAN-cleaned, so it asks whether a
filter agrees with itself. The experiment now prints a degeneracy refusal in place
of the rate. A fair head-to-head needs raw clouds on both sides and is still open.

## 2. Does D0093's filter invalidate exp/16?

Twenty molecules re-docked with no clustering between the dock and the write, and
exp/16 run on both sets — the same molecules, the same draw, one difference:

| | raw | filtered |
|---|---:|---:|
| poses | 10,000 | 8,376 |
| groups per molecule (median) | 217 | 174 |
| largest group | 14 (max 36) | 13 (max 40) |
| singleton groups | 48% | 47% |
| median within-group RMSD | **1.12 Å** | **1.08 Å** |
| 90th percentile | 1.63 Å | 1.62 Å |
| **worst group anywhere** | **3.76 Å** | **2.71 Å** |

**exp/16's conclusion survives.** The filter removes 16% of poses and adds 25% to
the group count — the removed poses are scattered ones that become singletons —
but group *quality*, the thing exp/16 exists to measure, is unchanged: 1.12 Å
against 1.08 Å at the median, 1.63 against 1.62 at the 90th percentile. Against
the 9.3 Å bag D0088 condemned, both are the same answer.

**One number is genuinely worse on raw data,** and it should be quoted: the worst
single group widens from 2.71 Å to 3.76 Å, which is above the 3.5 Å sweep bar.
The filtered clouds flattered the method exactly there, in the tail, which is
where a filter that removes outliers would flatter it.

## 3. A guard, from the way this nearly went wrong

`persist_raw_clouds.py`'s first run failed all 20 molecules on `No module named
'gemmi'` — the wrong environment — printed **"0 ok, 20 failed"**, and **exited 0**.
The chained step after it then ran exp/16 against zero raw clouds and also
reported success, in six seconds. Nothing in either output said the comparison had
not happened.

It now raises when nothing is persisted, naming the environment it needs, and the
chain aborts on a count below the minimum rather than proceeding. Catalogue
disguise #4: a stage that cannot distinguish an empty result from a clean one.

## What is still open

* A head-to-head against the shipped rule on **raw clouds for both**, which is the
  only comparison that is neither degenerate nor circular.
* Reproducibility is one molecule (`t4_716800c125a7`), the only one with five
  independent dockings.
* The tolerance these groups were built at is itself under question — D0094.
