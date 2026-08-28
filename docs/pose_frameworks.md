# Pose handling — which module is live, and for what

*Written 2026-08-27, because five modules in this repo group poses and nothing
said which one runs.*

This repo has accumulated **five** modules that turn poses into groups. They were
each written for a real reason and several are still load-bearing, so this is a
map rather than a cull. **The live splitter is named in
`config/target.yaml: splitting.method` — not in code, and not here.**

## The five

| module | represents a pose as | groups by | status |
|---|---|---|---|
| [`pose_modes`](../shared/pose_modes.py) | reactive-atom position + warhead direction | DBSCAN, eps 3.0 | **LIVE** — `splitting.method: warhead_dbscan` |
| [`pose_subsplit`](../shared/pose_subsplit.py) | whole-molecule heavy-atom RMSD | complete linkage, 2 Å, max 5 | **LIVE** — stage 2 of the above |
| [`pose_contacts`](../shared/pose_contacts.py) | **atom × residue** contact tensor, RMSF-weighted | complete linkage, per-molecule tolerance | **BUILT, NOT DEFAULT** — `splitting.method: contact_linkage` |
| [`pose_vector`](../shared/pose_vector.py) | **one number per residue** (min distance) | single linkage | **LIVE, DIFFERENT JOB** — fit score + representative over ~9 modes |
| [`pose_cluster`](../shared/pose_cluster.py) | heavy-atom RMSD | HDBSCAN, leaf selection | **SUPERSEDED** — experiments only, never adopted |

## The two that are easy to confuse

`pose_vector` and `pose_contacts` both describe a pose by what it touches, and
they were written three weeks apart without either knowing about the other. They
are not duplicates, and the differences are the whole point:

| | `pose_vector` (2026-08-04) | `pose_contacts` (2026-08-26) |
|---|---|---|
| granularity | **one number per residue** — min over all ligand atoms | **one number per (atom, residue)** |
| orientation | **lost** — a flipped pose touching the same residues is identical | **kept** — a flip puts a different atom near a given residue |
| weighting | none | per-atom, inverse predicted RMSF |
| linkage | **single** | **complete** |
| sized for | ~9 Vina modes (an O(n³) Python loop) | 500–6,000 poses (`pdist`) |
| built for | a fit score against a reference profile | splitting a cloud into groups |

### The linkage rationales contradict each other, and both are right

`pose_vector` says complete linkage is wrong: *"two poses that differ by a small
rotation form a chain of near-neighbours and belong together, and
complete-linkage would split them on the widest pair."*

`pose_contacts` says single linkage is wrong: it chains, which is precisely
D0088's defect — the shipped rule produced a 137-pose "mode" spanning 9.3 Å by
chaining.

**Both hold, for different n.** Over 9 Vina modes already spread by a
minimum-RMSD floor, chaining cannot run away and single linkage recovers
rotational near-neighbours. Over 500 poses filling a continuous cloud, chaining
is exactly what happens, and only complete linkage bounds a group's width. The
rule to carry forward is **the linkage must match the density of the cloud**, and
neither module's rationale transfers to the other's input.

## Which to use

* **Splitting a docked cloud into groups** → `pose_contacts.split_poses`, once
  `splitting.method` is switched. Until then, `pose_modes` + `pose_subsplit` via
  the screen.
* **Scoring one pose's arrangement against observed binders** → `pose_vector`
  (`reference_profile`, `fit_score`). `pose_contacts` has no reference profile
  and makes no claim about affinity.
* **Anything new** → `pose_contacts`. Do not add a sixth representation.

## Why `contact_linkage` is not the default yet

Every frame, ranking, sweep result and 100 ns run on disk was produced under
`warhead_dbscan`. Switching the default re-groups every cloud and makes them
incomparable in a way none of the artefacts would announce. The switch belongs
with the re-screen (#79).

What is measured and holds (D0092, D0095): median within-group RMSD **1.12 Å**
against the 9.3 Å bags of both predecessors; **91%** of groups holding ≥5 poses
found in all five independent dockings, against HDBSCAN's 1 of 3; **100%**
persistence under 12× deeper sampling.

What is open: the count does not saturate and **must never be reported as a
number of binding modes** (b = +0.69, no plateau); the per-molecule tolerance is
close to a constant and does not beat writing one number down (**D0094**); and
the framework has never been tested against the MD-validated pose, nor against
**SIFt** (Deng 2004), which is the published prior art for describing a pose by
its interactions.

## `pose_cluster` is superseded

HDBSCAN over heavy-atom RMSD. It has **no production caller** — only
`exp/4,5,6,7,8,9,10` and their tests, which are the record of why it was
rejected: it discards 29% of a cloud as noise, lost the MD-validated pose in 3 of
30 replicates, kept only 1 of 3 modes across an independent draw (#78), and its
cluster count grows linearly with sampling because it has no length scale
(D0090). It stays on disk so those experiments still run; it is listed in
`data/ready_to_delete.md` and should go when they are archived.

---

## Before you bump `run.topic`

`config/target.yaml` is global state, and **detached processes from other Claude
Code sessions read it.** Long-running supervisors survive the session that
started them, poll every few minutes, and resolve every path from `run.topic`.
Bumping it therefore redirects work already in flight, and a resume check that
finds the new topic empty reads as *"nothing has been done yet"* — which is how a
supervisor comes to launch a full sweep against a screen that has produced no
modes.

Measured 2026-08-27: two `overnight.sh` supervisors, 14 and 9 days old, started
from a session in `~/repos/murmurent`, were keeping the nac_v5 campaign alive.
Bumping the topic to `nac_v6` made one of them build an empty report tree for it
within five minutes.

```bash
ps -eo pid,etime,cmd | grep -E '[o]vernight|[s]weep_worker|[p]ipeline.py serve|[p]romote_to'
```

**Stop the supervisor before its children** — killing a watched child only makes
it respawn. And do not treat an idle process as finished: a `promote_to_bpmd` had
been resident for 10 days 22 hours with zero GPU usage.
