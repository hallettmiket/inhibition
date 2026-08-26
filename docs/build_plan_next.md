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

### 2.4a Frequency does not predict outcome, at any clustering tightness — **MEASURED**

@tt8804: *"I feel like a pose that shows up many times is more likely to be the
real pose but I guess not."*

It is the standard intuition and it is what `consensus` encodes. The project had
measured against it twice — **D0071** (neither ranking metric predicts pose
stability) and **D0073** (consensus *depletes* validated mechanisms) — but both
under the **shipped** clustering, where a mode routinely holds 65–86% of the
cloud. A group that large has a consensus near 1 whatever the molecule does, so
those tests may have been measuring the clustering rather than the intuition.

`exp/9_consensus_vs_outcome` re-asks it under tight clustering. For each of the
147 swept modes, the pose actually simulated is located inside an HDBSCAN
re-clustering of its own cloud, and that group's size is correlated against what
the trajectory then did.

| predictor | vs 5 ns attack-ready fraction (n = 147) | vs 100 ns max ligand RMSD (n = 15) |
|---|---:|---:|
| consensus, shipped clustering | ρ = +0.102 (p = 0.22) | ρ = +0.320 (p = 0.25) |
| consensus, HDBSCAN (tight) | **ρ = +0.012** (p = 0.89) | ρ = +0.121 (p = 0.67) |
| group size, HDBSCAN (tight) | **ρ = −0.004** (p = 0.96) | ρ = +0.130 (p = 0.65) |

**Null everywhere, and tighter clustering makes it *more* null, not less.** The
median group holding a simulated pose was 74 poses under the shipped rule and
**6** under HDBSCAN — a 12× change in what "shows up many times" means, with no
change in the (absent) relationship to outcome.

At n = 147 this has ~80% power to detect ρ = 0.2, so what is ruled out is
anything beyond a weak effect. Note also the sign against 100 ns RMSD is
*positive* — higher consensus, more drift — though far from significant.

**So D0071 and D0073 were not artefacts of loose clustering.** How often docking
returns to a pose carries no information about whether that pose survives
dynamics, and no amount of tightening the clusters creates any.

### 2.4b Consequence: clustering is de-duplication, not ontology — **@tt8804's reframe**

> *"we should just treat the 'noise' as singleton poses and reframe clustering as
> a cost saving mechanism to collapse the number of poses."*

This follows directly from §2.4a, and it resolves the noise question by
dissolving it:

* A pose HDBSCAN calls noise is a **group of one**, not a pose that failed to
  exist. `exp/9` already scores it that way — `tight_size = 1` for noise — and
  under that reading "how many poses agree with this one" stays defined for
  every pose in the cloud.
* If group size predicts nothing (§2.4a), then a mode is **not** a claim about
  how the molecule binds. It is a statement that these N poses are close enough
  that simulating all N would be buying the same answer N times.
* Clustering's job is therefore to **choose what to simulate under a budget**,
  and the property that matters is *coverage* — does the set of representatives
  span the cloud? — not *population*.

**What this changes downstream, and it is most of the cost of D0088:**

| | under "modes are real" | under "clustering is de-duplication" |
|---|---|---|
| `consensus` = mode_size / n_poses | a score term | **has no meaning as a score** — §2.4a |
| the 12-pose rank gate (D0084) | estimability of a proportion | **wrong question**; a singleton is a legitimate candidate |
| `viable_fraction` per mode | a property of a binding mode | a property of a *cluster*, i.e. of the clustering |
| what to simulate | the highest-scoring modes | a **covering set** of the cloud, budget-sized |

That is a bigger change than swapping the clusterer, and it is the honest
consequence of the measurement rather than a preference. **Open: what selects
the covering set, if not a score over populations.**

### 2.4c The gate is not what stops singletons — the SCORE is — **MEASURED**

`exp/10_gate_vs_all`, bdhi only, 312 molecules (62 skipped for the cloud/table
mismatch of §1.6a). Same HDBSCAN grouping and the same scores on both sides;
**only the gate differs.**

| | modes | scorable | score ≥ 4.0 |
|---|---:|---:|---:|
| gated (≥ 12 poses) | **639** | 444 | **137** |
| all groups, singletons kept | **55,816** | 31,055 | **2,767** |

Dropping the gate multiplies the candidate count **87×** — and adds **zero**
singletons above the floor. All 2,767 that clear it are multi-pose groups.

**Why, and this is the finding.** A singleton has 0 or 1 poses in the distance
window. With 0 it is unscorable. With 1, the conditional proportion is 0/1 or
1/1, and after empirical-Bayes shrinkage that is **exactly two numbers**:

| `conditional_eb` | singletons |
|---|---:|
| 2.8505 (in window, not viable) | 14,844 |
| 3.9642 (in window, viable) | 4,929 |
| NaN (never entered the window) | 18,350 |

**Two distinct values across 19,773 scorable singletons.** One of them — 2.8505
— is shared by 50.1% of *all* scorable groups. That is the flat step across the
right-hand panel: not a property of the molecules, the score's floor of
resolution.

So the answer to §2.4b's open question is that **`conditional_eb` cannot select
a covering set, by construction**. It is a shrunken proportion, and a proportion
estimated from one observation carries no information to order on. The 12-pose
gate is not the obstacle; removing it changes nothing because the score was
never able to rank these groups in the first place.

**What that leaves.** If clustering is de-duplication and the population score
is dead (§2.4a), the selection rule has to come from something other than
counting poses. Three candidates, none measured:

1. **Geometry directly** — rank a representative on its own approach distance
   and angle, with no population term at all. Every pose has these; a singleton
   is as measurable as a 200-pose group.
2. **Coverage** — choose representatives to span the cloud (a facility-location
   or max-min-distance pick) and spend the budget on breadth rather than depth.
   Selection stops being a score at all.
3. **Energy** — the one per-pose quantity we have and currently do not rank on,
   deliberately (#23/#30 found it uninformative for mode choice). Worth
   re-testing under tight grouping, since that finding has the same
   loose-clustering caveat §2.4a just cleared for consensus.
4. **Cheap intrinsic chemistry, and rigidity first** (@tt8804, raised in earlier
   issues and again here). Rotatable-bond count, ring fraction, heavy-atom
   count — free to compute, defined for every molecule, and **rigidity is a
   direct mechanistic proxy for the thing the sweep measures**: a pose with few
   rotatable bonds has fewer ways to leave the geometry it was placed in, so
   "will it stay as shown" is partly answerable before any simulation runs.
   Unlike the population terms it is a property of the MOLECULE rather than of
   the clustering, so §2.4a's null does not touch it.

   **Tested, and null — but the test is weak and should not be read as a
   refutation.** Over the 147 swept modes (119 molecules), against the 5 ns
   attack-ready fraction:

   | descriptor | ρ | p |
   |---|---:|---:|
   | rotatable bonds | +0.004 | 0.97 |
   | rotatable bonds / heavy atom | −0.002 | 0.98 |
   | ring count | −0.006 | 0.94 |
   | fraction sp3 | −0.068 | 0.41 |
   | heavy atoms | −0.016 | 0.85 |
   | TPSA | −0.072 | 0.39 |
   | clogP | +0.029 | 0.73 |

   **The predictor barely varies: rotatable bonds run 2–8, median 4.** T_4 is one
   core with R-group substitution, so the library is chemically homogeneous and
   there is not enough spread in rigidity to detect an effect that exists. This
   is range restriction, the same limitation #71 records for enrichment against
   outcome. The hypothesis needs a library that varies — it has not been tested
   here, only failed to be visible.

### 2.4d Granularity: is the pose space finite, and does HDBSCAN bound it? — **PREDICTION, recorded before the measurement**

@tt8804: *"given that our search area is bound by posebusters, is there still
infinite granularity within that area and does hdbscan compensate for infinite
granularity"* — and, on the objective: *"at the end of the day we just need to
catch the real pose with confidence."*

**Granularity is infinite, and that is a property of the space.** A pose is a
point in a continuous space — six rigid-body degrees of freedom plus torsions.
PoseBusters bounds the *volume*; within that volume there is always a pose
0.01 A from any pose already held. "How many distinct poses exist" therefore has
no physical answer, only an answer relative to a stated **resolution**.

**HDBSCAN does not supply one.** Its only size control is `min_cluster_size = 3`,
an absolute COUNT, not a distance. As a bounded volume densifies, local density
variations resolve at ever finer scales and it keeps subdividing; in the limit it
can return n/3 clusters. Nothing in it says *stop at 1.5 A*.

| method | length scale | count as density → ∞ |
|---|---|---|
| DBSCAN at fixed `eps` | yes | bounded |
| complete linkage at fixed cut | yes | bounded |
| **HDBSCAN** | **none** | **unbounded** |
| covering number at radius r | yes | bounded |

**This cuts against D0088, and the same property is responsible for both.**
HDBSCAN was adopted because it never produces a bag (widest mode 3.91 A against
the shipped rule's 9.30) — and it never produces one *because* it subdivides
until density says stop rather than until distance says stop. A virtue at 500
poses; potentially pathological at 50,000.

**Recorded before the numbers** (`exp/5` at 6,000 PoseBusters-filtered poses,
ladder to 5,000, both metrics):

* HDBSCAN cluster count **grows roughly linearly in n**, no plateau.
* Covering number at fixed r **plateaus** — bounded volume at fixed resolution
  has a finite answer, and deeper sampling only re-covers it.

If both plateau, the prediction is wrong and HDBSCAN is self-limiting on real
clouds, which is worth knowing. If HDBSCAN climbs while the cover flattens, then
**the covering number is the instrument for counting and HDBSCAN is not** —
though it may remain fine for *choosing* representatives.

### 2.4e What the count is actually for — **@tt8804's reframe of the objective**

> *"we just need to catch the real pose with confidence"* … *"if the counts
> taper off and are related logarithmically then we can find a 95% confidence
> spec and that's good enough for us"*

This makes the count instrumental rather than interesting. The quantity we care
about is **P(the cloud holds a pose within r A of the real one)**, and the
project already reasons this way one stage earlier: `docking.n_runs: 500` is
justified in config as *"500 runs give >= 95% probability of sampling at least
one pose within 2 A of the true one, provided the per-run hit rate exceeds
0.597%"*.

**Log growth is sufficient — a true plateau is not required.** If the covering
number is a + b·ln(N), it never flattens, but the marginal new coverage per added
pose falls as 1/N, so there is a finite depth at which 95% of the covering set
attainable at infinite depth is already held. That depth is a measured answer to
"how deep should we dock", replacing a justification that currently asserts a
sampling property it does not measure.

The spec would read: **cover the PoseBusters-allowed volume at resolution r,
with 95% confidence, at depth N** — and N becomes `docking.n_runs`.

**THE LIMIT OF THE CLAIM, and it must travel with it.** Coverage saturation says
we have covered *what docking can reach*, not what is real. If the true pose sits
where docking never samples, depth does not help — and two reasons that is live
here are already recorded: the receptor is rigid where the real complex would
relax (§2.1), and the scoring function steers the search. So the defensible
statement is *"we have covered the accessible space at r A with 95%
confidence"* — necessary, not sufficient, and notably stronger than anything the
pipeline claims today.

### 2.4f Saturation: measured, and the prediction was half right — **D0090**

`exp/5`, 6,000 poses docked (3 × 2,000, distinct seeds), **5,390 PoseBusters-valid**
(89.8% — a third independent confirmation of the ~90% rate in §1.3), ladder to
5,390.

| metric | exponent *b* in *a·n^b* | R² |
|---|---:|---:|
| HDBSCAN modes | **0.977** | 0.999 |
| modes + singletons | **1.025** | 1.000 |
| covering number @ 1.0 Å | 0.879 | 0.998 |
| covering number @ 1.5 Å | 0.753 | 0.995 |
| covering number @ 2.0 Å | **0.628** | 0.992 |

| poses | modes | modes+singletons | cover @ 2.0 Å |
|---:|---:|---:|---:|
| 500 | 58 | 237 | 275 |
| 2,000 | 237 | 985 | 634 |
| 5,390 | 657 | 2,794 | 1,035 |

**§2.4d's prediction was half right.** HDBSCAN grows linearly — b = 0.977, and a
linear fit at R² 1.000 against 0.816 logarithmic — confirming it has no length
scale and does not compensate for granularity. **But the covering number does not
plateau either.** It is sublinear, yet a power law with a positive exponent has no
asymptote: at 2 Å, doubling the poses still returns ~1.5× the distinct places.

**So §2.4e's 95% depth does not exist.** Log growth would have yielded one;
n^0.63 does not. The honest replacement names the resolution as a choice: *"at
2 Å and 500 runs we hold 275 distinct placements; doubling the depth adds ~50%
more."* Diminishing returns, not completeness.

**And de-duplication is weaker than the reframe assumed.** At production depth,
500 poses need **275** representatives at 2 Å — 45% collapse — and at 1 Å,
**97% of poses are their own representative**. Clustering is not mostly removing
repeats.

### 2.4g Why it does not saturate: the cloud fits inside the score's error bar — **D0090**

@tt8804: *"I don't understand how so many different poses can be such low
energy??"* **They are not low energy. They are indistinguishable.**

| | |
|---|---:|
| energy span within one molecule's cloud | **3.96 kcal/mol** (median) |
| poses within 0.5 kcal/mol of that molecule's best | 3% |
| within 1.0 | 15% |
| within 2.0 | **63%** |
| within 3.0 | **95%** |

An empirical docking score carries ~2–3 kcal/mol of error, and this project
already treats that as disqualifying — `state_of_the_project` §4 rejects
alchemical FEP partly for a blinded **2.44 kcal/mol RMSE**, and AutoDock's
function is not better than FEP. **63% of a molecule's poses sit within the
tool's own uncertainty of the best one.**

Two corollaries, measured:

* **Energy does not predict attack geometry**: ρ = **−0.093** over 236,313
  poses; viable poses average −5.76 against −5.51, a **0.25 kcal/mol** gap.
* **Distinct places are not energetically distinct**: ten HDBSCAN modes of
  `t4_716800c125a7` span **1.7 kcal/mol** in median energy.

**This explains the non-saturation.** A search saturates when the landscape has
basins to fall into. A flat landscape has none, so every run ends somewhere
slightly new and the count grows with the looking. It also kills candidate 3 in
§2.4c before it was tested: energy cannot select a covering set either.

**Checked, and the reactive receptor is NOT the cause** (`exp/12`, 7 molecules,
2,000 runs into each of the reactive and the plain 3IKD). The objection was that
our softened van der Waals parameters (`R_EQ_12 = 3.2`, `EPS_12 = 1.0`) flatten
the landscape by construction. **The plain receptor is flatter, in 7 of 7:**

| | reactive − plain (mean) | reading |
|---|---:|---|
| energy span across the cloud | **+1.89 kcal/mol** | reactive spans MORE |
| fraction within 2 kcal/mol of best | **−0.307** | reactive discriminates MORE |
| saturation exponent b @ 3.5 Å | +0.072 | essentially unchanged |

On the plain receptor **25–75% of poses (mean 53%)** lie within 2 kcal/mol of
the best, against 3.5–42% on the reactive one. Stock AutoDock on this pocket is
*less* discriminating than our modified setup. **§2.4g stands and is
strengthened**, and saturation is untouched — b is 0.28–0.45 in both arms.

(The experiment compares the reactive SETUP against the plain one; the reactive
arm also carries reactive typing and a flexible sidechain. It supports "our setup
is not what flattens the landscape", not a claim about the vdW term alone.)

### 2.4h The resolution is already chosen by the pipeline — **@tt8804**

> *"we already tolerate 0.35 RMSD, we just need to decide how the molecule should
> fit in different regions of the cloud."*

`md.sweep_survivor_rmsd_nm = 0.35` nm — **3.5 Å of ligand RMSD still counts as
"held"** at the next stage. Two docked poses closer than that are within the
tolerance of the very thing that will judge them, so resolving finer is resolving
distinctions the pipeline then discards. HDBSCAN's median mode is ~1.5 Å wide:
**it works four times finer than the tolerance downstream applies.**

At the pipeline's own tolerance the numbers become tractable, and the exponent
falls sharply with resolution:

| resolution | exponent b | centres for 6,000 poses | centres at 500 |
|---|---:|---:|---:|
| 1.0 Å | 0.885 | 3,632 (61%) | 454 |
| 2.0 Å | 0.649 | 1,104 (18%) | 275 |
| **3.5 Å** (MD tolerance) | **0.417** | **254 (4%)** | **96** |
| 5.0 Å | 0.257 | 56 (1%) | 32 |

Still no asymptote, but at 1 Å the count doubles with every doubling of poses,
while at 3.5 Å it takes ~10× the poses to double. **And a covering set IS the
"partition the volume" proposal**: the greedy centres are the representatives,
it carries an explicit length scale where HDBSCAN carries none, and singletons
stop being a special case — a pose alone in its ball is a ball with one member.

### 2.4i Is the covering set reproducible? — **MEASURED**

If a covering set is to replace mode clustering, it must survive the test the
incumbent already passed (§2.3). Same five independent dockings,
`exp/8/cover_reproducibility.py`:

| method | centres / replicate | pairwise recovery | in **all 5** |
|---|---:|---:|---:|
| HDBSCAN modes | 54–63 | **88.6%** | 41/59 (69%) |
| cover @ 2.0 Å | 88–97 | 69.5% | 37/92 (40%) |
| **cover @ 3.5 Å** (MD tolerance) | **34–43** | **83.4%** | **25/37 (68%)** |
| cover @ 5.0 Å | 12–16 | **93.8%** | 12/15 (80%) |

**HDBSCAN is not magically more reproducible.** Its 88.6% sits between the
3.5 Å and 5.0 Å covers — it is operating at *some* effective resolution and
delivering the reproducibility that resolution buys. It just does not let you
choose which.

At the pipeline's own tolerance the cover matches HDBSCAN's core-set fraction
(68% vs 69%) with **34–43 centres instead of 54–63**, and at 5 Å it beats it
outright on both measures with a quarter of the representatives.

**That is the argument for the covering approach, and it is not "it is more
reproducible".** It is that **reproducibility becomes a dial**. Resolution,
count and reproducibility move together and explicitly; under HDBSCAN they move
together and silently.

**The arbitrary-start caveat is resolved, and it resolved by not mattering.**
Greedy farthest-point is deterministic only after its first centre, and the
implementation took index 0 — whatever AutoDock wrote first, which is a
positional choice and this project's defining defect shape. Restarting from the
**medoid** (a property of the cloud, so two independent dockings begin at the
same place in the pocket) changes almost nothing:

| radius | start = index 0 | start = medoid |
|---|---:|---:|
| 2.0 Å | 69.5% | 68.2% |
| 3.5 Å | 83.4% | 82.6% |
| 5.0 Å | 93.8% | 92.9% |

Within noise, and marginally *lower*. So the recovery rates above are **not** a
lower bound depressed by write order — **the limit is intrinsic to the cloud and
the resolution.** The medoid start is kept as the default regardless, because
depending on file order is not a property worth having even when it costs
nothing; `--start first` reproduces the old behaviour.

### 2.4j Does chemistry predict saturation? No — and there is no cutoff to place — **MEASURED**

`exp/11`, 10 molecules stratified on rigidity and size, 2,000 runs each,
PoseBusters-filtered, covering number fitted at the pipeline's 3.5 Å tolerance.

| ident | rotb | heavy | fsp3 | PB % | **b** | centres@500 |
|---|---:|---:|---:|---:|---:|---:|
| t4_98951476e4f8 | 2 | 22 | 0.38 | 88.7 | 0.368 | 73 |
| t4_e6cf2d8e26d2 | 2 | 24 | 0.50 | 91.7 | 0.385 | 66 |
| t4_b2b9dab376eb | 3 | 19 | 0.60 | 87.0 | 0.438 | 58 |
| t4_61d3ee480cc3 | 3 | 24 | 0.59 | 88.3 | 0.384 | 46 |
| t4_7625f3e0d7f3 | 4 | 21 | 0.54 | 94.4 | 0.396 | 74 |
| t4_df8e56d56792 | 4 | 27 | 0.45 | 87.2 | 0.479 | 86 |
| t4_04fb46a6929a | 5 | 22 | 0.57 | 92.0 | 0.452 | 64 |
| t4_45aa6c98bf0e | 5 | 28 | 0.50 | 87.3 | 0.429 | 110 |
| t4_f0b05f412e7f | **6** | 23 | 0.56 | 89.2 | **0.431** | 74 |
| t4_d6cd64168a1c | **8** | 29 | 0.67 | 88.2 | **0.431** | 99 |

**Nothing predicts it.** Across 14 tests, **zero** reach p < 0.05 — against 0.7
expected by chance alone. Rotatable bonds come closest at ρ = +0.55, p = 0.102.

**And the interim signals did not survive the extremes**, which is the reason to
have run all ten:

| | at n = 8 | at n = 10 |
|---|---|---|
| b vs rotatable bonds | +0.683, p = 0.062 | **+0.546, p = 0.102** |
| centres@500 vs fsp3 | −0.719, **p = 0.045** | **−0.293, p = 0.412** |

The two most flexible molecules — **6 and 8 rotatable bonds** — both returned
**b = 0.431**, dead on the median of 0.430. If flexibility drove saturation they
should have been the worst; they were average. The apparent trend at n = 8 was
the middle of the range doing the work.

**The practical answer, which stands regardless of significance: b spans
0.368–0.479 across a 4× range in flexibility.** Every molecule sits in a narrow
band around 0.43. There is no subpopulation that saturates and another that does
not, so **there is no cutoff to place** — a rigidity filter would exclude
molecules at b = 0.45 to keep ones at b = 0.39, discarding real chemistry to buy
essentially nothing.

**The bias trap this experiment was written to avoid did not need to be sprung.**
§2.4c's worry was that b would track rigidity and tempt a rotatable-bond filter
that bought a saturating search by discarding flexible chemistry. It does not
track it, so the question does not arise — and the honest reading is that
**non-saturation is a property of the pocket and the scoring function, not of the
ligands.** That is consistent with §2.4g (the landscape is flat to within the
score's error) and with `exp/12` (the plain receptor is flatter still).

**Fifth independent confirmation of the PoseBusters rate**: 87.0–94.4%, mean
**89.4%** across these ten.

### 2.4k The volume partition is refuted — **D0091**

@tt8804 proposed partitioning the pose-occupied volume at 3 Å instead of
clustering by similarity. I measured it, reported that it held, and an
adversarial audit demolished it. All three load-bearing findings verified
independently:

**The cloud fills the box.** Heavy-atom extent **25.48 × 25.42 × 25.16 Å** in a
26 Å cube — 97–98% of every axis. AutoDock hard-clips to the grid, so the volume
question had one possible answer before any pose was docked. Catalogue entry #12
in a new coordinate system: the **search box** standing in for the **receptor
cavity**. PoseBusters, credited with the bounding, removes ~10%; the box removes
the rest.

**The partition rebuilds the bag.** 49 cells, largest holds **1,758 of 6,000
poses**, median **6.08 Å** and max **9.13 Å** RMSD inside a cell **3 Å** wide.
D0088 condemned the shipped rule at 137 poses spanning 9.3 Å — this is 13× the
membership at the same width. Unfixable by tuning: two poses can share a centroid
and be 180° flips, so a partition of ℝ³ cannot bound a distance in configuration
space.

**The exponent contrast dissolves.** One OLS slope was fitted across a curved
log-log plot, and every rung was a subsample of one pooled cloud rather than an
independent docking. That bias inflates the volume ladder's low rungs ~19% while
leaving the covering number untouched (177 vs 177) — it **manufactured** the
contrast. Corrected: volume **b ≈ 0.32**, cover **b ≈ 0.325**.

**Two numbers in §2.4 were wrong and are corrected here.** "14% of the box" used
point occupancy (verified 14.4%), the column `exp/13`'s own docstring disowns;
sphere occupancy is ~34%. "1.71× poses → 1.029× cells" was one draw; over 20 the
mean is 1.102×.

**What survives:** at 6,000 poses the cloud has swept ~a third of its own box and
is still adding territory at b ≈ 0.32 by every metric. Diminishing returns, not a
bound.

**The one informative experiment:** re-dock at 1.5× and 2× box. If the envelope
expands proportionally the bound was always the box; if not, the intuition was
right and I tested it in a way that could not show it.

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

---

## 2.5 Does the group count taper? (exp/17)

*2026-08-26. @tt8804: "lets see the groups as a function of poses generated like
before. if it tapers off we just accept the huge number and go to ranking" and
"also can we check if this residue contact space grows with poses".*

Full record in [D0092](../decisions/D0092-contact-space-is-fixed-the-group-count-climbs-because-6000-poses-undersample-it.md).

### 2.5a The answer to the first question is no

Ladder on the raw 6,000-pose cloud, tolerance held fixed at the molecule's own
0.73 Å:

```
    poses   groups   per 1k new   largest   singletons
      100     79.2                    4.2        80%
      500    285.0        429.3       9.4        57%
    1,000    464.2        327.2      13.8        47%
    3,000    948.7        200.0      27.7        33%
    6,000  1,438.0        140.5      37.0        26%
```

`b = +0.693`, and the species-accumulation fit implies no finite plateau. The
marginal rate is falling — 429 new groups per 1,000 poses at n=500, 140 at
n=6,000 — but it is falling toward a ceiling we never reach. Replicated across 12
production molecules: median `b = +0.668`, range +0.589 to +0.780, none above 0.9.

**And loosening the tolerance is not a way out.** It buys a lower count by
rebuilding the bag:

| tolerance | groups at n=6,000 | largest group |
|---:|---:|---:|
| 0.73 Å (RMSF) | 1,438 | 37 |
| 1.50 Å | 209 | 194 |
| 2.00 Å | 64 | 433 |
| 3.00 Å | 9 | 2,004 |
| 3.50 Å (the sweep bar) | 3 | 3,147 |

There is no tolerance at which the count saturates *and* the groups stay tight.

### 2.5b The answer to the second question is what rescues it

The **space** does not grow. Diameter exponent `+0.019`, mean pairwise separation
exponent `+0.001`; 60× more poses widened the diameter by 1.10× — and the maximum
of a larger sample is larger by construction, so that is an upper bound on the
effect. The 99th percentile of pairwise distance is 3.11 Å at every depth from 500
to 6,000.

So the climb is **undersampling of a fixed region**, not expansion of the region.
That distinction is the whole decision: expansion would mean the pose set is an
artefact of runtime, undersampling means the region is a property of the molecule.

Two supporting numbers:

* **Effective dimension 3.54** out of 420 coordinates offered (28 atoms × 15
  residues), from `N(ε) ~ ε^−d`. That is the right order for rigid-body placement,
  and is the evidence the metric tracks pose rather than noise. Still rising with
  n, so it is a floor.
* **Groups never move.** 100% of n=500 groups have an n=6,000 counterpart inside
  the tolerance, median displacement 0.254 Å against 0.73 Å, over 5 draws and
  1,460 shallow groups — and 100% for non-singletons and for groups of ≥5 alone.
  Growth is entirely tail: 26% of deep groups are singletons, 61% hold ≤3 poses,
  while 65% of the cloud lives in the 431 groups of ≥5.

This is exactly the property D0088's rule lacked, and for a nameable reason: an
absolute, molecule-owned tolerance makes group identity a property of the region,
while HDBSCAN's density criterion makes it a property of the draw.

### 2.5c What it decides

Proceed to ranking, under two conditions:

1. **Artefacts say *groups*, never *modes*,** and never report the count as a
   property of the molecule. It is a monotone function of docking depth with no
   plateau; reporting it as a mode count invites precisely the reading D0088 found
   in the shipped pipeline.
2. **Rankings are compared at fixed docking depth,** because the group population
   is depth-dependent even though each group is not.

### 2.5d The defect this turned up — read before citing exp/14–16

[D0093](../decisions/D0093-the-file-named-allposes-is-not-all-poses-it-is-dbscan-cleaned.md),
catalogue entry #26. `<topic>_allposes/` is **not** all poses:
[`nac_screen_v2.py:501`](../scripts/nac_screen_v2.py#L501) writes only poses whose
DBSCAN label is in `mode_ids`, so 21% of every production cloud — the scattered
poses — is missing. exp/14, exp/15 and exp/16 all read that path, which means
**every candidate replacement for DBSCAN so far has been measured on clouds
DBSCAN had already cleaned.** At n≈400 the raw cloud gives 241–254 groups where
the filtered clouds give 109–118.

exp/17's conclusions are unaffected — it reads the raw 6,000-pose cloud
throughout, which is why the discrepancy was visible at all. exp/16's headline
numbers need re-running on raw clouds before they are cited.

### 2.5e Still open

* Re-run exp/14–16 on raw clouds (D0093).
* Reproducibility of contact groups across the five independent dockings —
  persistence within one cloud is established, agreement between clouds is not.
* Extent, dimension and persistence are **one molecule**; only the exponent is
  replicated (12 molecules, ≤450 poses each).
