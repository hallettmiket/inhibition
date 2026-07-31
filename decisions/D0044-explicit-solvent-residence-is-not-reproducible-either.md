---
id: D0044
title: Residence is not reproducible in explicit solvent either — 5 replicates settle it
date: 2026-07-31
status: accepted
approach: shared
decided_by: '@mhallet'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - decisions/D0038-implicit-solvent-residence-is-a-property-of-the-solvent-model.md
  - shared/gromacs_explicit.py
  - scripts/merge_gromacs_results.py
  - integration/app/pose3d.py
evidence:
  - '235 of 240 replicates completed: 48 candidates x 5 replicates x 10 ns explicit TIP3P'
  - 'all 5 failures are ONE candidate, t1_1224c0ee20c2, which exploded in every replicate'
  - 'T_1: replicate max/min ligand RMSD ratio median 3.53, 90th pct 5.71, max 13.95'
  - 'T_2: median 2.16, 90th pct 3.28, max 3.75'
  - 'T_1: 22 of 24 candidates disagree between replicates by more than 2x'
  - 'T_2: 15 of 24 disagree by more than 2x'
  - 'T_1 mean ligand RMSD across candidates 0.19-1.08 nm; T_2 0.17-0.49 nm'
  - '286 of 288 trajectories analysed (the 2 failures are the exploded candidate)'
---

# Five replicates, and the answer is that one run tells you nothing

## What this was run to settle

D0038 withdrew a claim. Two T_1 candidates had been reported as leaving the
pocket under implicit solvent (ligand RMSD 9.00 and 7.30 nm), and that was
attributed to the water model. It did not reproduce: re-running under the SAME
model gave 1.75 and 0.59 nm, a 5.1x and 12.5x swing. The corrected lesson was
that a per-candidate residence claim needs replicates in either solvent model,
and the obvious hypothesis was that explicit water would be better behaved --
it has structure, it must be displaced, it should hold a ligand in place.

240 runs were launched to test that: 48 candidates, 5 replicates, 10 ns of
TIP3P each, differing only in their velocity seed.

## It is not better behaved

| | candidates | median replicate max/min | 90th pct | worst | disagree > 2x |
|---|---|---|---|---|---|
| T_1 | 24 | **3.53x** | 5.71x | 13.95x | **22 / 24** |
| T_2 | 24 | **2.16x** | 3.28x | 3.75x | 15 / 24 |

Five 10 ns trajectories of the *same molecule*, in the *same solvent model*,
differing only in initial velocities, give mean ligand RMSDs that span a factor
of **3.5 on a typical T_1 candidate and 14 on the worst**. The spread between
replicates of one molecule is comparable to the spread *between molecules*
(T_1 means run 0.19-1.08 nm across candidates).

So D0038's finding was not about implicit solvent. It was about 10 ns of
molecular dynamics. Explicit water does not rescue reproducibility; it narrows
it somewhat (T_2's median 2.16x against T_1's 3.53x) and leaves it far too wide
to rank on.

## What this retires

**Any per-candidate residence claim from a single 10 ns trajectory, in either
solvent model.** Not "treat with caution" -- there is no claim to make. A
number whose replicate spread is 3.5x cannot distinguish two candidates whose
means differ by less than that, and almost none of ours differ by more.

This is now the second time residence has been reported and withdrawn. The
first withdrawal (D0038) blamed the solvent model; that was too generous to the
method, and the replicates say so.

## What survives, and it is worth something

**As a QC flag, on the aggregate.** A candidate whose replicates ALL show large
displacement is different from one where a single run wandered. The merge now
carries `explicit_rmsd_replicate_min/max/sd` and the max/min ratio beside the
mean, so a reader can see the spread rather than infer it.

**One candidate is genuinely unstable.** `t1_1224c0ee20c2` failed all five
replicates with LINCS constraint violations -- the system explodes. It
completed in the single-run campaign, so the replicates converted an apparent
success into a real finding: that pose or those parameters are not physical.
That is exactly the kind of thing residence measurement is good for, and it is
categorical rather than a ranking.

**The comparison between solvent models stands and is now better founded.**
Explicit ligand RMSD is lower than implicit by a mean of 0.29 nm (T_1) and 0.47
nm (T_2), computed against replicate-averaged values rather than single runs.

## Two silent failures caught on the way

Both would have "succeeded", and neither would have been visible in a log line.

**The analysis discovered no replicates.** `discover()` looked only for
`wd/prod.xtc`, the flat layout of the original campaign; the replicates are in
`wd/repN/prod.xtc`. It would have re-analysed the 48 old single runs and
reported them as the 240-replicate result. Found by checking the analysis path
BEFORE the campaign finished rather than after.

**The merge then discarded them anyway.** `drop_duplicates("candidate_id")`
keeps the first row per candidate, and `discover()` emits the flat run first --
so even after the discovery was fixed, the merge reported the old single-run
numbers. The log read "24 rows carry explicit-solvent results" in both cases.
It was caught only because the count of 24 was inspected against an expectation
of five values per candidate.

The pattern is the one this project keeps meeting: the correction was written,
and the consumer one layer up kept reading the old thing.

## The general lesson

Sampling was the answer to a question nobody had asked precisely enough. The
project ran replicates to find out whether specific candidates dissociate, and
the replicates answered a prior question instead: whether this measurement is
reproducible at all at this trajectory length. It is not. Anything built on top
of it -- rankings, per-molecule verdicts, the "does the pose survive water"
framing -- had to wait on that answer and nobody had asked for it.

A measurement should be shown to be reproducible before it is used to
discriminate. That is cheaper than any of the analysis that came after it here.
