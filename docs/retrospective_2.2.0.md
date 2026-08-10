# 2.2.0 “Chalcopyrite” retrospective

*Closed 2026-08-10. Companion to
[`outline_2.2.0.md`](outline_2.2.0.md) (what was promised),
[`framework_2.2.0.md`](framework_2.2.0.md) (what was built, canonically) and
[`retrospective_2.1.0.md`](retrospective_2.1.0.md) (the version it follows).*

**The one-line version: 2.1.0 ended by discovering that its score ranked poses on
a quantity carrying no information. 2.2.0 replaced the unit of work — a mode, not
a molecule — found that the crystallographic pose was in our output all along and
the selection step was throwing it away, and then, on its last night, ran the
positive control and discovered that the chemistry had been sound the whole time
and it was the READOUT rejecting it, in three independent ways.**

Forty commits. The version's own founding premise did not survive it.

---

## 1. What was built

| | |
|---|---|
| **Pose splitting** | Poses cluster into binding **modes** on the reactive atom's position and the direction the warhead faces. Deliberately *not* whole-molecule RMSD (D0062), *not* docking energy (#23/#30 — would re-import the defect this version exists to remove), and *not* distance-to-SG or the NAC angle, because those **are** the score |
| **A mode is a candidate row** | @tt8804's design ruling. Ranking, selection and the GUI consume modes unchanged, so the change removed downstream work rather than adding it |
| **Screen** | One row per binding mode; default raised to **500 runs** after measuring that 95% pose coverage needs ~300 |
| **Reactive docking** | Protonates at pH 7.4 (D0074) |
| **Score selection** | Chosen on evidence — Wilson lower confidence bound, reference molecules, and a convergence check |
| **Attack-geometry gate** | Free pose check, then a 10 ns sweep (#32); the 10 ns sweep predicts 100 ns attack geometry |
| **MD priority** | Pre-registered before the six 100 ns runs, then reported against the pre-registration |
| **Throughput** | Weekend worklist, concurrent sweep/MD worker pools, asymmetric per-class allocation (priority 50, the rest 5) |
| **GUI** | The catalogue viewer, adopted as the foundation every future results interface follows |
| **Docs** | A canonical framework description with every claim measured, and a per-stage ownership map |

## 2. What 2.2.0 got right

**Pre-registration stopped being a gesture.** Four pre-registrations were written
before the work they governed — MD priority, co-folding, the attack sweep, score
selection. The payoff is visible in the one that came back negative: **BPMD does
not predict residence** (ρ = +0.213, n = 6). A verdict written down in advance is
allowed to disappoint, and this one did, on the version's own headline tool.

**Auditing before spending, twice, and both times it paid.** Six defects were
found by auditing *before* the library re-dock rather than after. Separately, the
2.2.0 gate was found to pass **zero molecules** — before the run, not after a
weekend of empty output. Neither audit was prompted by a failure; both were
prompted by the habit.

**The version killed its own premise rather than quietly re-planning.** The
outline's central selling point was that steps 1–3 needed *no new simulation* —
the dividend from 2.1.0 persisting its poses. That premise turned out to be
wrong, and the build plan says so in its title. Recording a dead premise in the
open is what makes the outline trustworthy next time.

**A design ruling that removed work.** "A mode IS a candidate row" could have been
a new parallel structure with its own ranking, its own selection and its own
panel. Making it a row meant every existing consumer worked unchanged. The best
architectural decisions in this project keep having this shape.

**Measure on this target, then act on it.** The 93.3% finding (§4) was not
inherited from a benchmark paper. It was measured on 15 deposited Pin1 complexes
at the production protocol, and the pipeline changed the same week.

## 3. What went wrong

### 3.1 The recurring defect, still recurring

The class named in `recap_2.0.0.md` §3 and again in `retrospective_2.1.0.md` §3.1
— **a value taken by position, name, or inheritance rather than by identity**.
2.2.0's instances:

| defect | the value taken wrongly | cost |
|---|---|---|
| **`md_movie.py` Cys113** | residue **number** 113, but the MD system renumbers from 1 — Cys113 is residue 63 | a **glutamate** rendered in sticks and labelled as the target cysteine, in a published picture |
| **anchoring selection** | the **ARGMAX** of anchoring | selected strained poses |
| **`frac_in_range`** | divided by the wrong denominator | a rate that was not the rate it was named |
| **the 2.2.0 gate** | thresholds inherited from the prior definition | passed **zero** molecules |
| **co-folding harness** | hydrogens lost in conversion | blocked **339 poses** |
| **`noncovalent_dock_run.py`** | `RECEPTOR_PDBQT` still hardcoded to `6VAJ_prepared.pdbqt` | see §3.3 — **still live** |
| **persistence docstring** | described a window that had been removed | documentation asserting a behaviour the code no longer had |

The `md_movie` instance deserves attention because of *how* it got through:
`elevation_report.py` already had the offset right and **refuses to label a
structure whose residue types do not match**. The guard existed, was correct, and
the new styling code simply did not go through it. A guard only protects the path
that calls it.

It affected no computed value — the warhead→SG series locates the SG atom
independently and was right throughout. Only the picture was wrong. That is a
patch by the versioning rule, and it is also the most quietly dangerous kind of
error this project produces, because a picture is exactly what a reader trusts
without checking.

### 3.2 Silent failures, again

`retrospective_2.1.0.md` §3.2 concluded: *a stage that produces no output should
fail the run, not log a non-zero exit and continue.* 2.2.0 produced two more of
the same shape:

- **The gate that passed zero molecules.** It did not error. It returned an empty
  set, which is a valid result, and would have been read as "nothing qualified."
- **`mdprio` plots that never rendered**, warhead geometry never plotted, PBC
  silently breaking the calculation.

Both were caught by looking, not by the system complaining. The 2.1.0 lesson was
recorded and not yet implemented — that is the honest reading.

### 3.3 The receptor is not settled, and the code disagrees with itself

**D0059 makes 3IKD the receptor, replacing 6VAJ**, on the chemist's instruction
relayed by @tt8804, with measured support: cross-docking into 6VAJ ranked the
crystal pose #1 in **0 of 82** cases against 5/82 self-docked, top-3 3.7% vs
22.0%. 6VAJ's pocket is induced-fit around sulfopin, so it was always the wrong
frame to dock strangers into.

But:

- **D0059's status is still `proposed`**, not accepted.
- **`config/receptor.yaml` still pins `pdb_id: 6VAJ`**, including the prepared
  `.pdbqt` the docking box points at.
- **`shared/noncovalent_dock_run.py:61` still hardcodes `6VAJ_prepared.pdbqt`**
  and defines a `SIX_VAJ` receptor.
- Meanwhile the benchmark and reference-screen paths call `resolve_3ikd_ian()`,
  which *refuses to run against the wrong 3IKD*.

So one half of the codebase guards fiercely against using the wrong receptor
while the other half defaults to the receptor that was replaced. Both are
populated, both are plausible, and which one you get depends on which entry point
you came through. This is the catalogue defect at the level of the project's most
load-bearing artifact, and `config/receptor.yaml` warns in its own header that
switching receptors mid-choreography invalidates every join.

**This is the first thing the next version should close.**

### 3.4 The new GUI foundation dropped the controls

`retrospective_2.1.0.md` §1 records the 2.1.0 GUI as *"a Ranking 2.1.0 panel with
references beside candidates"*, and §2 names **references as a yardstick, not a
gate** among the things that version got right. The catalogue viewer — adopted in
2.2.0 as the pattern every future results interface follows — has none of it.

`mdprio_combine.py` builds the rail from three globs (`attack_sweep/`,
`md_residence/`, `rank_v2/`) and contains no reference to controls, references,
crystal poses, Sulfopin or ATRA. `crystal_controls.py` writes to
`00_outputs/blacksmith/crystal_controls/`, and its only consumer is
`crystal_reactant.py`. The output is otherwise terminal.

This is worse than a missing panel because of **#47**: the warhead classes with
crystal structures and measured kinetics score **last** on the near-attack
criterion, while a class with no measured Pin1 activity scores **first**.
`crystal_controls.py` exists precisely to separate the two readings of that — bad
docking versus a bad criterion. Until it is on screen beside the candidates, the
catalogue presents a ranking whose validity is under active falsification, with
the falsifying evidence built, runnable, and invisible.

A foundation is the right thing to declare; declaring one that quietly drops a
control the previous version was praised for keeping is how a good pattern
inherits a bad default.

## 4. The finding that ends the version

From #30, measured at the production protocol — AutoDock-GPU reactive, `--nrun
200`, 3IKD, flexible Cys113 — against **15 deposited Pin1 complexes**, with the
ligand superposed into 3IKD's frame and scored by element-aware symmetry-corrected
RMSD:

| the crystallographic pose is… | |
|---|---:|
| present anywhere in the 200 poses | **93.3%** (median best **1.03 Å**) |
| still present after `KEEP_TOP = 20` by energy | 46.7% |
| present in the top 10 the score actually reads | **33.3%** |

The single miss is 9INP at 2.14 Å — a near-miss on a 2.0 Å bar, not a failure to
find the site.

**The search is not the problem. It never was.** We find the right answer nine
times in ten and then discard it during selection, because selection was ordered
by energy — the quantity 2.1.0 already proved carries no information about
reaction competence (ρ = +0.009 across 115,300 poses).

And the replacement works. Mode-0 selection, calibrated against the same 15
crystal complexes at 500 runs with each molecule docked twice, recovers the
crystal pose **93.3%** of the time against energy's **60.0%**.

That is the version in one exchange: 2.1.0 found the score was noise; 2.2.0 found
that the cost of that noise was not a worse answer but a *discarded correct
answer*, and that changing the unit of selection recovers it.

**What this says about the version.** The three headline numbers — 93.3% present,
46.7% surviving the cut, 33.3% surviving the score — were all computable from the
production protocol at any point in 2.1.0. Nothing new had to be built to learn
that we were throwing away the right pose. It took someone asking #30.

## 5. What carries forward

**Keep:** modes as candidate rows; pre-registration before measurement; auditing
before spending; the 500-run default and the coverage measurement behind it; the
catalogue viewer as the interface pattern; measure-on-this-target.

**Fix, in order:**

0. **Run the docked Sulfopin at 100 ns.** D0077 shows the `rx_*` crystal-reactant
   controls model an adduct as a Michaelis complex, so they cannot answer the
   positive-control question. The docked pose enters at a valid near-attack
   geometry and can. This is the first experiment of the next version.
1. **The dwell filter (D0076).** Re-derive `n_visits` at the save interval and
   report the raw count beside it. Costs no GPU — the raw number is already in
   every row. Pre-register the new reading first: it is being changed after
   seeing which molecules it rejected.
2. **Never rank across mechanisms.** The 58-survivor list is a cross-mechanism
   ranking and 56 of the 58 come from the two laxer criteria. Within SN2 the
   controls top their own chemistry.
3. **The receptor split (§3.3).** Accept or reject D0059, then make
   `config/receptor.yaml` and `noncovalent_dock_run.py` agree with the decision.
   Nothing else should be run until this is settled.
4. **A stage that yields nothing must fail the run.** Twice recorded, not yet
   implemented.
5. **Route every renderer through the guard that already knows the residue
   offset**, rather than trusting each new styling path.

**Test:** whether mode assignment is stable across re-runs of the same molecule —
the outline named unstable mode counts as a failure condition, and 500 runs was
chosen for coverage, not yet demonstrated for reproducibility.

**Unresolved, carried in and still carried:** the SN2 150° threshold — now with
evidence on both sides, since the controls clear it repeatedly and never hold it; the
secondary pocket Reddi 2023 reports; N-activated acrylamides at 97% of T₃
(D0066), still without a chemist's ruling; the within-class rigidity confound.

## 5b. The positive control, run on the last night — and what it found

Added 2026-08-10, after the retrospective's first draft. @tt8804 asked the
question the version had not yet answered: **would our own screen have caught
Sulfopin?**

Answering it produced three records in a chain, each qualifying the one before,
and the chain is the most useful thing 2.2.0 produced.

**D0075 — it says no.** Sulfopin through the production protocol gives 1 mode,
465 poses, 47 reaction-competent, then fails the 10 ns sweep with **zero sustained
visits**, ranked 104 of 234. Liu-2022-ZL-Pin13 and Juglone fail it too. Meanwhile
58 of 233 candidates pass. On the face of it, the criterion admits a quarter of
our generated matter and none of the chemistry known to react with Cys113.

**D0076 — why, and it is a defect.** `rx_7F0M` entered attack geometry **13
separate times** and scored zero visits, because `MIN_DWELL_PS = 100` requires
each excursion to *last* 100 ps and the sweep saves every 19.96 ps — five
consecutive frames. The pre-registration chose visits over occupancy in as many
words: *"a covalent reaction needs ONE good approach, not sustained occupancy"*.
The implementation then filtered on persistence. Re-derived on raw visits and
compared within mechanism, as #47 requires: **Liu-2022 beats 20 of 20 SN2
candidates**, Sulfopin's crystal form 18 of 20, its docked form 17 of 20, against
an SN2 candidate median of **0**.

**D0077 — and the control itself was the wrong shape.** 6VAJ's ligand is a
covalent adduct. Cleaving the bond leaves the reactive carbon **1.98 Å** from SG,
*below* the 2.8 Å window floor, and the 180° angle is *constructed* by placing the
halogen along the S→C vector rather than measured. Equilibration relaxes it to
3.57 Å / 100.9° before the first production frame. Our own docked pose enters at
**3.36 Å / 156.8°** — inside the window, over the SN2 bar. The 5.01 Å RMSD between
them at 1.45 Å centroid separation is a pivot, not a mislocation.

**What the chain amounts to.** The screen's chemistry was sound. Three separate
readout defects — a persistence filter on an observable chosen for the opposite
reason, cross-mechanism ranking, and a control built in the wrong state — made it
look otherwise. Each was found by asking what a number meant rather than by a test
failing, which is the same route §3 records for every other defect this project
has caught.

The 100 ns results stand and are worth quoting on their own: **Liu-2022 99.95%
engaged and held**, Sulfopin 78.9%, Juglone 47.5%.

## 6. Judged against its own failure criteria

`outline_2.2.0.md` §5 stated in advance what would make 2.2.0 a failure. Scoring
honestly:

| stated failure condition | verdict |
|---|---|
| A score that ranks well and still gives Sulfopin a zero | **Avoided.** Score selection was made on evidence including reference molecules |
| Mode counts that change between re-runs of the same molecule | **Not yet tested.** Coverage was measured; reproducibility was not |
| Elevating a mode that BPMD then shows was the wrong one, repeatedly | **Open, and complicated** — BPMD turned out not to predict residence at all (n = 6) |
| **Another silent stage** | **Failed.** Two: the gate that passed zero, and plots that never rendered |

One of four failed outright and one is untested. Writing the criteria down in
advance is what makes that sentence possible.

---

## 7. Numbers worth remembering

**93.3%** — how often the crystallographic pose is somewhere in our 200.
**33.3%** — how often it is somewhere the score would let anyone see it.

2.1.0's numbers were about how often docking *finds* the pose. 2.2.0's are about
how often we *keep* it. The gap moved from the instrument to our own handling of
its output, which is the more embarrassing place for it to be and the easier one
to fix.

And two more, from the last night, which say the same thing again:

**13** — independent approaches Liu-2022-ZL-Pin13 made into attack geometry.
**0** — the number our own readout scored it.

## 8. Where 2.2.0 leaves the project

Every defect this version found was in the reading, not the chemistry. Pose
splitting works; the screen puts a known covalent inhibitor into a valid
near-attack geometry; the crystal pose is in the output. What kept failing was
each successive attempt to *summarise* that into a number — an energy window
(#23/#30), a dwell filter (D0076), a control built in the wrong state (D0077),
and a ranking pooled across incomparable mechanisms.

That is a better position than 2.1.0 ended in, and it is worth saying plainly
because the negative results make it easy to read the opposite. 2.1.0 closed not
knowing whether anything in the pipeline discriminated. 2.2.0 closes knowing the
pipeline finds the right chemistry and that the instruments pointed at it have
been measuring the wrong quantity — which is a fixable problem, and four of the
fixes cost no GPU at all.

**Not stamped with a successor.** Whether the next version is 2.3.0 or 3.0.0 is
[#46](https://github.com/hallettmiket/inhibition/issues/46) and turns on whether
"the ranking predicts nothing" invalidates measurements or their interpretation.
Nothing in this retrospective settles it.
