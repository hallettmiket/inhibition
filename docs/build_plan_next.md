# Build plan — the next release

*A LIVING DOCUMENT. Written with @tt8804 as each pipeline step is discussed, one
at a time. Sections appear as we settle them; a section marked **DRAFT** has not
been agreed, and a section marked **OPEN** is a question we have not answered.*

*Started 2026-08-17. Companion to
[`pipeline_in_plain_language.md`](pipeline_in_plain_language.md), which
describes the pipeline as it is TODAY. This one describes what we are changing
it to, and why.*

---

## 0. Status, and the version question — **OPEN**

| step | change under discussion | section | state |
|---|---|---|---|
| 1. Generate | — | — | not discussed |
| **2. Dock** | **add PoseBusters as a validity gate; quota on valid poses** | **§1** | **measured — D0089** |
| **3. Group into modes** | **replace the two-stage splitter — D0086, D0088** | **§2** | **in discussion** |
| 4. Score (near-attack) | — | — | not discussed |
| 5. Rank | — | — | not discussed |
| 6. Shortlist | — | — | not discussed |
| 7. Triage sweep | — | — | not discussed |
| 8. Production MD | — | — | not discussed |
| 9. BPMD | — | — | not discussed |

**This is probably a MAJOR release, and we should say so before we start.** By
the CHANGELOG's own rule, MAJOR means *previously reported numbers are invalid
and must be re-measured*. Changing which poses exist changes the mode split, the
per-mode counts, `viable_fraction`, `enrichment`, the ranking and the shortlist —
so nothing from 3.1.0 stays quotable beside it. 3.1.0 is itself still in
progress, so the numbering (3.2.0 vs 4.0.0) is a decision for @tt8804, not a
consequence of this document.

**Consequence worth stating early:** every change we agree here costs a full
re-screen, and they should therefore be batched into one release rather than
landed one at a time. That is an argument for discussing all nine steps before
running anything.

---

## 1. Docking — PoseBusters as a validity gate — **DRAFT** · recorded as **D0089**

### 1.1 What we do now

Each in-scope molecule is docked with AutoDock-GPU at a **fixed `--nrun 500`**.
Whatever poses come back are what the molecule gets. The pose count is therefore
*variable and unremarked*: measured on `nac_v5`, molecules return anywhere from
~418 to ~486 poses from the same 500 runs, and nothing downstream records that
one molecule was scored on 418 observations and another on 486.

No physical validity check is applied to a pose at any point. A pose that
clashes into the protein is scored on its near-attack geometry exactly like one
that does not.

### 1.2 What @tt8804 proposes

> *"Add PoseBusters to this step since it is quite lightweight. Change 500 rounds
> of docking → 500 poses, to running docking checked by PoseBusters to get up to
> 500 PoseBusters-stable poses as a quota."*

So: **the quota moves from runs to valid poses.** Dock until 500 poses pass
PoseBusters, rather than docking 500 times and keeping whatever appears.

Two things this buys:

1. **Physical validity.** Poses that clash into the protein, or that carry
   absurd internal energy, stop entering the scoring population.
2. **An equal denominator.** Every molecule is scored on 500 poses rather than
   on however many its docking happened to return. `viable_fraction` is a
   proportion, and it is currently measured on populations that differ by ~15%
   between molecules with nothing recording it.

The second is arguably the bigger win and was not the stated motivation.

### 1.3 Measured: what PoseBusters actually rejects

`posebusters 0.6.5` is already installed in `dwi_cheminf`. No new environment.

Sample: **1,500 poses**, 25 molecules drawn at random from each of the three
warhead families, run in `dock` mode (ligand + protein, no crystal reference)
against the prepared 3IKD.

| | pass rate |
|---|---:|
| **overall** | **90.60%** |
| acrylamide | 91.6% |
| bdhi_c4 | 88.2% |
| bdhi_c5 | 92.0% |

Per molecule: **min 70%, median 90%, max 100%**; 65 of 75 molecules have at
least one failing pose.

**Only 2 of the 22 checks ever fail:**

| check | failure rate |
|---|---:|
| `minimum_distance_to_protein` | 7.53% |
| `internal_energy` | 2.00% |
| the other 20 | 0.00% |

### 1.4 Why 20 of the 22 checks can never fail here

This is structural and worth writing down, because it decides how much of
PoseBusters is actually doing work.

AutoDock-GPU varies **rigid-body placement and torsion angles only**. Bond
lengths, bond angles, ring conformations and stereochemistry are inherited
unchanged from the input conformer, which meeko/RDKit generated and
energy-minimised before docking began.

So `bond_lengths`, `bond_angles`, `aromatic_ring_flatness`,
`double_bond_flatness`, `all_atoms_connected`, `sanitization` and the rest are
**tests of our conformer generator, not of our docking** — and they are the same
answer for all 500 poses of a molecule. They will pass at 100% until we change
how ligands are prepared.

The two that vary with the pose are the two that fail:
`minimum_distance_to_protein` (where the pose sits relative to the protein) and
`internal_energy` (which torsions can change).

**Implication for the design:** we are not adding "22 physical checks". We are
adding **two**, plus twenty that certify the conformer generator once per
molecule. That is still worth having — but the honest description of this change
is *a protein-ligand clash gate and an internal-energy gate*, and the build
should say so rather than claim broader coverage than it has.

### 1.5 Compute cost

**The gate itself is cheap.** 232 ms/pose, measured sustained over 1,500 poses
(not a warm-up artefact — the rate held from the first batch to the last).

A full screen at the new quota is ~561 × 552 ≈ **310,000 poses ≈ 20 CPU-hours**,
or about **25 minutes wall** on the project's 50-worker cap
(`shared/compute.py`). Negligible against everything else in the pipeline.

**The docking increase is the real cost.**

| | |
|---|---:|
| runs needed for 500 valid poses at 90.6% | **~552** |
| GPU increase, on the mean | **+10.4%** |
| runs needed for the worst molecule sampled (70%) | **714** |
| GPU increase, worst case | **+43%** |

The spread matters more than the mean. **A fixed over-request sized on the
average will undershoot the quota for roughly half the library** — which would
reintroduce the unequal-denominator problem this change is partly meant to fix,
while charging us the extra GPU time anyway.

### 1.6 The risk we have not resolved — **OPEN**

PoseBusters flags a clash when any ligand heavy atom is closer to a protein
heavy atom than **0.75 × the sum of their van der Waals radii**
(`clash_cutoff: 0.75`, `radius_type: vdw`, confirmed in the installed
`dock.yml`).

For the reactive carbon against Cys113's sulfur:

```
0.75 × (1.70 + 1.80)  =  2.625 Å
```

**Our near-attack window starts at 2.8 Å.**

So a pose sitting at the ideal SN2 approach distance is **0.175 Å above being
declared a clash** — and `minimum_distance_to_protein` is precisely the check
doing three quarters of the rejecting. Compounding it: we dock into a
*reactive* receptor whose van der Waals parameters were deliberately softened
(`R_EQ_12 = 3.2`) to pull the warhead toward the sulfur. We built a receptor
that encourages close approach and are now proposing a filter that penalises it.

**If PoseBusters rejects attack-ready poses at a higher rate than ordinary
ones, this gate spends GPU time buying reaction-incompetent poses.** That is
the same shape as D0082, where a filter rejects the one molecule known to work.

**Measurement in flight.** All poses of the 18 molecules with the most viable
poses, joined to the screen's own `viable` flag by the established
(mode, pose_idx) correspondence, cross-tabulated against the PoseBusters
verdict. Result goes in §1.7.

### 1.6a A data defect this measurement surfaced — already fixed in code

The join above refused two molecules for a pose-count mismatch (461 SDF vs 457
rows; 474 vs 456), so I measured how widespread it is:

> **105 of 561 molecules (18.7%) have a persisted pose cloud that does not match
> their own scored rows.** Differences run from −42 to +39 poses; the SDF holds
> more in 43 cases and fewer in 62.

For those molecules the cloud cannot be joined to its own measurements at all,
which is #44's rule — *the cloud must come from the same run that produced the
scores* — broken in practice.

**The cause is found and already fixed on this branch** (`a7cca45`,
`nac_screen_v2.py`): the all-poses SDF was written under `if not
adest.exists()`, so a re-screened molecule kept the **previous** run's cloud
beside the current run's table. It is now always rewritten.

**But the `nac_v5` data on disk predates the fix**, so the 18.7% stands for
every number currently derived from those clouds. It resolves only on a
re-screen — which this build requires anyway. Worth stating explicitly because
it means *the current pose clouds are not a safe baseline to compare the
PoseBusters gate against*: §1.10's criterion 4 must be measured on a re-screen
with the fix in, not against what is on disk today.

### 1.7 Result of the attack-ready cross-tabulation — **MEASURED**

**The risk in §1.6 does not exist. PoseBusters protects the poses we want.**

5,834 poses across 13 molecules (5 refused by the join guard, §1.6a), each pose
carrying the screen's own `viable` flag:

| | PoseBusters pass |
|---|---:|
| attack-ready (viable) | **98.57%** (n = 1,819) |
| not attack-ready | **92.80%** (n = 4,015) |
| odds ratio | **5.35** in favour of viable poses |
| Fisher exact | **p = 1.9 × 10⁻²³** |

The clash check — the one I expected to do the damage — fails **0.93%** of
attack-ready poses against **6.48%** of the rest. It rejects seven times more of
what we do not want.

**Why the §1.6 arithmetic was right and the conclusion was wrong.** The 2.625 Å
clash threshold sits *below* our 2.8 Å window, so it bites the poses that are
**too close** — the ones our own criterion already rejects for being inside a
formed bond rather than approaching one:

| reactive C to SG | poses | PB pass | clash fails |
|---|---:|---:|---:|
| **< 2.8 Å** (below our window) | 42 | **45.2%** | **54.8%** |
| 2.8–3.2 Å (window floor) | 580 | 89.5% | 10.0% |
| 3.2–3.6 Å | 1,872 | 96.2% | 3.5% |
| 3.6–4.2 Å | 1,267 | 96.5% | 2.6% |
| 4.2–6 Å | 835 | 92.6% | 6.4% |
| > 6 Å | 1,238 | 95.7% | 3.6% |

PoseBusters and the near-attack criterion **agree**, independently, about which
poses are physically real. That is a small mutual validation of both and worth
recording as such.

Residual: the 2.8–3.2 Å band loses 10%, which is the one place the two rules
genuinely overlap. 58 poses of 580 in this sample.

### 1.7a But it does not change the answer — **criterion 4 fires**

Per-mode `viable_fraction` recomputed with the failing poses removed, over the
75 modes holding ≥ 12 poses:

| | |
|---|---:|
| median change in `viable_fraction` | **+0.00 pp** |
| mean change | +0.67 pp |
| largest change | 5.03 pp |
| rank correlation, before vs after | **ρ = 0.9989** |
| modes whose rank moves at all | 19 of 75 (median move **0** places, max 5) |
| top-10 overlap | **10 / 10** |
| top-25 overlap | **25 / 25** |

**§1.10's criterion 4 is met: the gate does not change what the pipeline
recommends.** ρ = 0.9989 is above the ≳ 0.99 threshold written down before the
measurement.

**Sample limitation, stated plainly.** These 13 molecules were chosen as *those
with the most viable poses* — the best-behaved end of the library — because that
is where the §1.6 risk would show most sharply. It is the right sample for the
question it was built to answer and the **wrong** sample for estimating a
library-wide effect. The §1.3 pass rate (90.6%, random across families) is the
unbiased number; 94.6% here reflects the selection. A library-wide rerank
comparison belongs in the re-screen, not here.

Two of 75 modes fell below the 12-pose estimability gate after filtering. Under
a true quota that cannot happen, since the run tops up to 500 *valid* poses.

Decision rule, **set before the number is known**:

* **If the odds ratio is ≈ 1** (PoseBusters rejects attack-ready and ordinary
  poses at indistinguishable rates) — adopt the gate as proposed, all checks on.
* **If attack-ready poses fail materially more often** — do **not** adopt
  `minimum_distance_to_protein` as-is. Two options, in preference order:
  1. exempt the anchor atom pair (reactive atom ↔ Cys113 SG) from the clash
     check and keep the check for every other atom, since the close approach we
     want is *specifically* at that pair;
  2. drop `minimum_distance_to_protein` from the gate entirely and keep the
     other 21 — which, per §1.4, means keeping `internal_energy` and twenty
     conformer certifications.
* **If attack-ready poses fail materially LESS often** — that is evidence the
  criterion and the clash check agree, and worth recording as a small
  independent validation of both.

### 1.8 Implementation options — **DRAFT**

The quota loop interacts with two known AutoDock-GPU properties.

**(a) The `--nrun` ceiling.** Measured on this build: `--nrun 5000` corrupts its
own stack *and still writes a `.dlg`*; `--nrun 10000` exits 0 with no output.
`dock()` already refuses all three failure shapes. Any quota loop must therefore
cap the per-call `--nrun` well below 2,000 and top up with additional calls
rather than one large request.

**(b) Seeding.** Docking is seeded for reproducibility (#77) — the same seed
returns the same poses, so a top-up call must use a *different* seed or it
returns the identical cloud. The seed sequence has to be deterministic and
recorded, or we lose the reproducibility that #77 bought.

Three shapes:

| option | how | cost | risk |
|---|---|---|---|
| **A. Flat over-request** | one call at `--nrun ceil(500 / 0.906) = 552`, keep the first 500 valid | simplest, one call | undershoots for ~half the library — the molecules with the *worst* pose quality get the *smallest* populations, which is backwards |
| **B. Iterate to quota** | call at 500, gate, top up in batches of ~100 with successive seeds until 500 valid or a hard cap | exact quota | more calls, each with AutoDock start-up cost; needs the seed ledger |
| **C. Adaptive** | call at 500, measure this molecule's pass rate, size one top-up call from it | one extra call, near-exact | slight overshoot/undershoot; still needs the seed ledger |

**Recommendation: C**, with B's hard cap as a backstop. It costs one extra
AutoDock call per molecule, hits the quota for nearly everything, and its
failure mode (a molecule that still falls short after the top-up) is *recorded*
rather than silent — see §1.9.

### 1.9 What must be recorded — **DRAFT**

The current screen records `nrun: 500` as a constant. Under a quota it becomes a
per-molecule outcome, and a constant in the config would then be a lie.

Per molecule, the run must persist:

- `n_runs_requested` — total AutoDock runs actually spent
- `n_poses_returned` — before the gate
- `n_poses_valid` — after the gate
- `quota_met` — bool; **false is a result, not an error**
- `pb_fail_reasons` — counts per failing check
- `seeds` — the seed sequence used, for reproducibility

**A molecule that cannot reach 500 valid poses must be stamped, not dropped and
not silently short.** That is the same rule as `docked_species_ok` (#58): the
substitution we are guarding against is a molecule quietly entering the ranking
on a smaller population than everything it is ranked against.

### 1.10 Pre-registered failure criteria

This change is **rejected** if any of the following is true after
implementation:

1. **It selects against reaction competence.** Attack-ready poses are rejected
   materially more often than ordinary poses and the fix in §1.7 does not
   resolve it.
2. **The quota is not actually met.** More than 10% of molecules end the screen
   below the 500-valid-pose quota, leaving the unequal denominators the change
   was meant to remove.
3. **It costs more than +25% GPU on the library mean.** At that point the same
   GPU time is better spent on the sweep, where §"What we can and cannot claim"
   in the plain-language doc says the real uncertainty lies.
4. **It changes nothing measurable.** If the ranking produced with the gate is
   indistinguishable from the ranking without it (rank correlation ≳ 0.99 on the
   same seeds), then we have spent GPU time and added a dependency for no
   change in what the pipeline recommends — and the honest thing is to record
   that and not adopt it.

Criterion 4 deserves emphasis. **It is the most likely outcome**: 90.6% of poses
already pass, the rejects are concentrated in a check that may be firing on the
poses we want, and mode assignment is a clustering over hundreds of poses that
is unlikely to move much when 9% are removed. We should measure it explicitly
rather than assume the gate helped because it sounds like it should.

---

## 1.11 Recommendation — **DRAFT, for @tt8804**

**Adopt it — but for the second reason, not the first, and record that it does
not change the ranking.**

The case *against* the stated rationale: it was proposed to improve pose
quality, and measured, it does not improve the answer. ρ = 0.9989, top-25
unchanged. Spending **+10.4% GPU** to buy an identical shortlist is not a good
trade on its own.

The case *for* adopting anyway, three reasons that survive §1.7a:

1. **The equal denominator is a real methodological fix.** Molecules are
   currently ranked against each other on 418–486 poses with nothing recording
   it, and `viable_fraction` is a proportion. The quota fixes that. *Note this
   benefit comes from the QUOTA, not from PoseBusters* — quota-ing on raw pose
   count would deliver it for zero CPU. If GPU is tight, that is the cheap
   version of this change.
2. **It is a guard that can fail.** Twenty of the 22 checks certify the
   conformer generator and currently pass 100%. They would fire if ligand prep
   regressed — which is exactly the class of silent defect this project keeps
   finding. Cheap standing insurance at ~20 CPU-h per screen.
3. **It is the field's standard, and #66 needs cheap credibility.** "All poses
   PoseBusters-validated" pre-empts a reviewer question. And the measurement
   itself is reportable: **9.4% of raw AutoDock poses on this target are
   physically invalid**, and they are disproportionately the non-reactive ones.

**What must go in the release notes, or this becomes a claim we cannot support:**
that the gate was adopted for validity and reproducibility, **not** because it
improved the ranking — because it measurably did not.

---

## 2. Pose splitting — replace the two-stage splitter — **DRAFT**

Recorded position: **D0086** (two proposed fixes measured and rejected) and
**D0088** (modes come from pose similarity alone, HDBSCAN) — the latter still
`proposed`. This section adds what @tt8804 asked for on 2026-08-17:
reproducibility, and whether the noise matters.

### 2.1 There is no answer key for this state, and that decides the method

@tt8804: *"sulf is not relevant here. it is already bound in an induced fit. we
are trying to model a transient pose here that is pre-covalent and a naive
protein."*

Sulfopin's crystal is the **covalent adduct** in an **induced-fit** pocket. We
model a **transient pre-covalent** encounter against a **naive** receptor —
the wrong state on both axes. And a transient complex is transient, which is why
it is not in the PDB. **No experimental ground truth exists for what we are
clustering.**

Two consequences:

1. D0088's "the pose we know is right" **overstates what those three references
   are**. `exp/4_election`'s `REFERENCES` table holds three T_4 candidates —
   never synthesised, never assayed — whose poses our own docking produced, our
   own ranking elected, and our own 100 ns MD did not dislodge. The code is
   honest (*"a mode a run elected AND a trajectory confirmed"*); the prose is
   not. They are the best available reference **because nothing better exists**,
   not because they are known-correct.
2. **Accuracy is unavailable, so the case must rest on internal consistency** —
   width, purity, reproducibility, scaling. Those need no answer key, and they
   are most of what D0088 already measured.

### 2.2 The argument D0088 should lead with

The purity test labels each pose by whether it reaches attack geometry, then
asks whether a mode's poses agree — a mode at 0.4 is two populations under one
label.

**The shipped rule clusters on reactive-atom position and warhead direction,
which ARE the distance and angle terms of that label.** It groups along the
exact axes the purity test measures, so it should score artificially well. It
has every advantage on this metric.

It still loses: widest mode 9.30 Å against 3.91, largest 137 poses against 14,
p90 width 7.86 Å against 2.52. **A rule that cheated on the metric and still
came last** is a far harder result to argue with than a single MD-confirmed
pose landing in a big group.

### 2.3 Reproducibility — **MEASURED** (`exp/8_hdbscan_reproducibility`)

@tt8804: *"is it reproducable. lets check the reproducability first."*

A mode that does not survive an independent draw of the pose cloud is not a
binding mode; it is a partition of one sample. Five independent 500-run dockings
of `t4_716800c125a7`, distinct seeds, each clustered on its own; modes matched
between replicates by medoid heavy-atom RMSD ≤ 2.0 Å — the same in-place metric
the clustering itself uses.

| | modes per replicate | noise | largest mode | pairwise recovery | in ALL 5 |
|---|---|---:|---:|---:|---:|
| **HDBSCAN** | 54–63 | 26–33% | **11–20** | **88.6%** (79.7–94.7) | **41 of 59 (69%)** |
| shipped | 2–4 | 0% | **266–349** | 65.8% (25.0–100) | **1 of 3 (33%)** |

**HDBSCAN is markedly the more reproducible rule**, on both measures.

Two things about the shipped row deserve their own sentence:

* **Its largest mode holds 266–349 of ~400 poses** — 65–86% of the entire cloud
  in one "binding mode". That is not a mode; it is the cloud with a label.
* **Its mode count swings 2 → 4 across replicates**, a 2× change in how many
  ways the molecule is said to bind, from re-running the same docking. Only one
  of its three modes survives all five draws.

@tt8804's proposed proxy — *"as long as HDBSCAN scales logarithmically with
poses we could consider it reproducible, since we generate about the same number
of poses"* — is sound reasoning, and this measures the thing directly rather
than through the proxy. Stable count is **necessary** for reproducibility, not
sufficient: the same *number* of modes could be different modes each time, which
is close to what the shipped rule does. **Scaling still matters for a different
reason** and is §2.5.

### 2.4 The noise is a fringe, and it is mildly protective — **MEASURED**

The open worry in D0088 was that 29% of poses become noise, and that a rare
transient pose might be exactly what lives there. Measured across the same five
replicates, against the screen's own `viable` flag:

| | |
|---|---:|
| overall noise rate | 28.2% |
| **P(noise \| attack-ready)** | **23.5%** |
| **P(noise \| not attack-ready)** | **29.3%** |
| risk ratio | **0.80** |

**An attack-ready pose is *less* likely to be discarded as noise than an
ordinary one.** The noise is a sparse fringe, not where the good poses hide.

Caveat worth carrying: the per-replicate spread is wide — 13.7% of viable poses
were noise in r5 against 38.0% in r2. The mean is reassuring; the variance is
not yet explained, and on a single molecule it should not be over-read.

### 2.5 Still open

- **Scaling.** No HDBSCAN saturation run exists — `exp/5` covered `shipped` and
  a `fine` recipe only. It matters because `consensus = mode_size / n_poses`: if
  mode count climbs with depth, the denominator moves under every per-mode
  score. Needs a deep cloud, which is not on disk (deepest persisted is 439
  poses) and would need docking.
- **One molecule.** Everything in §2.3 and §2.4 is `t4_716800c125a7`. D0088's
  quality figures cover three. Neither is a library.
- **What replaces `consensus`.** At 54–63 modes over ~400 poses, a mode holds
  ~2% of the cloud. `consensus` and the 12-pose estimability gate were both
  sized against 2–5 modes per molecule; under HDBSCAN nearly every mode fails
  both. **This is the real cost of adopting D0088, and it is not the clustering
  — it is everything downstream of it.**

---

## 3-9. Remaining steps — **not yet discussed**

Placeholders, so the order of work is visible. Each gets its own section, with
the same structure: what we do now, what changes, what it costs, what would make
us reject it.

---

## Appendix — measurements behind this document

| what | where |
|---|---|
| PoseBusters pass rate, 1,500 poses, 3 families | §1.3; raw CSV in the session scratchpad |
| per-check failure breakdown | §1.3 |
| 232 ms/pose timing | §1.5 |
| clash threshold 2.625 Å vs the 2.8 Å window | §1.6, from the installed `dock.yml` |
| attack-ready cross-tabulation | §1.7, pending |

**Reproduce:** the sampling scripts are in the session scratchpad and should be
promoted to `exp/` before this document is used to justify a decision. A
measurement that cannot be re-run is a claim, not a result.
