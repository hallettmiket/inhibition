---
id: D0112
title: The pose movie measured an atom named C10 rather than the warhead, and the sweep rail now ranks on the 100 ns gate
date: 2026-09-02
status: accepted
approach: shared
decided_by: '@twu383'
origin: user
supersedes: []
superseded_by: null
affects:
  - scripts/elevation_report.py
  - scripts/mdprio_report.py
  - scripts/sweep_report.py
  - scripts/shortlist_report.py
  - scripts/sweep_combine.py
  - shared/results_shell.py
  - tests/test_movie_reactive_atom.py
evidence:
  - '`elevation_report.surface_payload` selected the reactive atom as `next(i for i, a in enumerate(f0) if a[0] == "MOL" and a[2] == "C10")` -- by PDB atom name'
  - 'measured over 98 finished nac_v8 sweeps: the atom named C10 was the reactive atom in 0 of 98; movie frame 0 disagreed with the plot beneath it by a median of 3.11 A and up to 10.08 A; only 2 of 98 agreed within 0.05 A'
  - 'after the fix the two agree in 98 of 98, max |diff| 0.0049 A, which is the viewer 2-dp rounding'
  - '`tests/test_crystal_pose_audit.py` already recorded the reason: "6VAJ calls it C10; the other five covalent Pin1 entries call the equivalent atom C19, C14, C24, C12 and C3"'
  - '`nac_series` picked the warhead class by the FIRST matching SMARTS in warhead_classes_10.csv: 62 acrylamides reported as naphthoquinone_c2 (row 5 beats row 7), 34 bdhi_c4 as bdhi_c5 (row 3 beats row 4)'
  - 'that mislabelling was HARMLESS to every number: the same reactive atom in 98 of 98 and the paired classes share a mechanism, so distance, angle and occupancy are bit-identical (max |diff| 0.0000)'
  - '`nac_series` took Cys113 SG by matching resname CYS and keeping the LAST match; 3IKD has Cys57 (resid 7) and Cys113 (resid 63), so it was correct only because 63 sorts after 7'
  - 'all 98 modes were selected under 3.0 A docked; 0 of 98 are still under 3.0 A at the first production frame; median drift over the 300 ps equilibration is +1.51 A'
runbook: null
---

# The movie was measuring the wrong atom

@twu383, 2026-09-02, looking at a mode selected under 3.0 A whose movie opened
at **5.35 A**: *"audit why the movie shows a starting warhead distance higher
than 3.0 and does not match the rmsd plots"*.

Two separate reasons, and only the first is a defect.

## 1. The defect: the warhead was selected by atom NAME

The plot and the movie sat on the same page, described the same quantity, and
asked two different questions.

| | how it found the reactive atom |
|---|---|
| the plot (`mdprio_report.nac_series`) | the molecule's warhead **SMARTS** |
| the movie (`elevation_report.surface_payload`) | the atom **named `C10`** |

```python
i_c10 = next(i for i, a in enumerate(f0) if a[0] == "MOL" and a[2] == "C10")
```

`C10` is sulfopin's name for its reactive carbon in 6VAJ -- the one molecule
this report was written for. On a combinatorial library the names come from
antechamber and mean nothing.

**Measured over 98 finished nac_v8 sweeps: `C10` was the reactive atom in 0 of
them.** The viewer's "warhead->SG" readout was a median of **3.11 A** from the
truth and up to **10.08 A**; only 2 of 98 happened to agree within 0.05 A. After
the fix all 98 agree, with a maximum difference of 0.0049 A -- the viewer's own
2-decimal rounding.

### Why it looked right

The label said "warhead->SG" and the number was a plausible distance in a
plausible range, sitting under a plot of the same named quantity. Nothing
raised, because `C10` exists in most of these ligands -- just not as the warhead.

**The project had already written down why this fails**, in a test, in as many
words: *"6VAJ calls it C10; the other five covalent Pin1 entries call the
equivalent atom C19, C14, C24, C12 and C3. Anything that selected it by name
would pass on 6VAJ and quietly pick a different atom on the next structure."*

And the containing function is the one with the residue-numbering guard: forty
lines above, it **raises** if a protein residue does not match the name it is
about to label. The protein's identity was verified and the ligand's was
assumed.

### The fix is the class

`surface_payload` no longer has a default. It takes `reactive_idx` (from the
shared resolver) or, to ask for the old behaviour explicitly, `reactive_name`.
Passing neither raises. The `--reactive-idx` / `--reactive-atom-name` CLI flags
turn a pin that could not announce itself into a choice visible in the command
line -- catalogue #32/#35, the keyword argument whose non-default value no
caller ever passes.

`mdprio_report.reactive_atom(cand, rep)` is now the single resolver, and
`nac_series` was refactored onto it, so the plot and the movie cannot drift
apart again.

## 2. Not a defect: the movie starts after equilibration

Even measured correctly, that mode's frame 0 is **3.58 A**, not the sub-3 A it
was selected at. That is real, and it is D0110:

| | |
|---|---|
| docked, at selection | min 2.83, **median 2.89**, max 3.00 A |
| first production frame | min 3.15, **median 4.45**, max 6.46 A |

**All 98 modes were selected under 3.0 A. None is still under 3.0 A when the
movie opens.** The median warhead drifts **+1.51 A** during the 300 ps of
unrestrained NVT/NPT that precedes production. The movie has never shown the
docked pose and was never going to.

The sweep page now says so under the legend, because "selected under 3 A" and a
movie opening at 3.6 A otherwise read as a contradiction.

## 3. Two latent traps found in the same function, both closed

**The warhead CLASS was taken by file order.** `nac_series` walked
`warhead_classes_10.csv` and took the first class whose SMARTS matched: 62
acrylamides came back as `naphthoquinone_c2` (row 5 beats row 7) and 34
`bdhi_c4` as `bdhi_c5` (row 3 beats row 4).

**It changed no number, and that was checked before it was claimed.** The paired
classes match the same reactive atom (98 of 98) and share a mechanism, so every
distance, angle and occupancy is bit-identical -- maximum difference 0.0000.
What was wrong was the label, and what was fragile was that its harmlessness
depended entirely on which rows happen to sit above which in a CSV. The resolver
now asks the candidate frame what the molecule *is*.

**Cys113's sulfur was taken as "the last CYS SG in the frame."** 3IKD has Cys57
(resid 7) as well as Cys113 (resid 63), so this was correct only because 63
sorts after 7 -- file order again, and exactly what catalogue #38 says not to
lean on. Now selected by residue number through `md_movie.PIN1_OFFSET`, agreeing
with `surface_payload`, which had always done it properly.

## 4. The rail now ranks on the gate

@twu383: *"can you update the gui to show this priority"*.

The sweep rail ranked on **max ligand RMSD**, which is the pose-stability half
of D0111's gate -- so a molecule that sat perfectly still facing the wrong way
topped the list. It now ranks on the gate itself:

| tier | |
|---|---|
| 0 | clears the gate -- engaged >= 60% **and** pose held (green rule) |
| 1 | pose held, under the engagement bar -- the near misses |
| 2 | pose left the site |

Engagement descending within each tier, then max RMSD. The headline is
engagement because engagement leads the gate; stability sits beside it; and the
progress bar is now distance-to-the-gate rather than "lowest RMSD so far". The
legend states the thresholds by **reading them from the gate**, so the caption
cannot describe a rule the order does not follow.

Current standing, 105 finished: **0 clear it**, 76 held but under the
engagement bar, 29 left.

## 5. Guards

`tests/test_movie_reactive_atom.py`, 6 tests, each verified against a mutation.
The fixture places `C10` 20 A from the sulfur and the real warhead at 3.0 A, so
a regression reports an absurd number rather than a plausible one:

* an index selects the right atom; an index from another molecule is **refused**
  rather than clamped
* there is **no implicit default** -- calling with neither argument raises
* selection by name still works when explicitly asked for (sulfopin needs it)
* `nac_series` calls the shared resolver and has no SMARTS scan of its own
* no file in `scripts/` or `shared/` compares a ligand atom against a `C<n>`
  name -- **walked over the AST, not the source text**, because the first
  version of that test flagged the docstring in `surface_payload` that quotes
  the removed line as an example. That is the same false positive the version-pin
  test produced on four docstrings; prose describing a defect is not the defect.
