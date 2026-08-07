# 2.2.0 “Chalcopyrite” — build plan

*What gets implemented, why, what could go wrong with each, and how each one is
accepted or rejected. Written 2026-08-07, after the measurements in #23/#29/#30
changed what is possible. Companion to [`outline_2.2.0.md`](outline_2.2.0.md),
which stated the intent; this states the build.*

---

## 0. Read this first — the outline's logistics premise is dead

`outline_2.2.0.md` §4 says:

> **Steps 1–3 need no new simulation.** That is the dividend from 2.1.0
> persisting its working.

**That is not true, and the whole order of work depends on it.** `nac_screen_v2`
persists per-pose rows for the **top 20 by energy only** — verified: every
molecule has exactly 20 rows, `energy_rank` 1–20. The other 180 poses are
discarded at write time and no longer exist.

What survives per molecule:

| | persisted? | consequence |
|---|---|---|
| aggregate counts over all 200 (`n_in_range`, `n_viable_given_in_range`) | **yes** | `enrichment_conditional` **is** computable from existing data |
| per-pose geometry (distance, angle, viable) for all 200 | **no** — 20 only | `anchor_quality` over the full population is **not** computable |
| pose coordinates for all 200 | **no** — 20 only | **pose splitting is not computable.** Clustering 20 poses selected *by energy* into modes clusters the survivors of the exact filter we are removing |

So the central feature of 2.2.0 requires re-docking the library. That is not a
disaster — measured from the v2 run's own frames, ~7.8 s/molecule/shard, so 5,765
molecules on 4 shards is **≈3 GPU-hours** — but it reorders everything and it
must be done **once**, carrying every change that needs a re-dock.

---

## 1. The build

### 1.1 Persist the population, not the survivors *(new — forced by the above)*

**Change.** `nac_screen_v2` stops writing "the 20 lowest-energy poses" and starts
writing a pose set that can be split. Per-pose geometry for **all 200**;
coordinates for a mode-diverse subset.

**Considerations.**

- **Storage is the reason `KEEP_TOP` exists, and it is not actually binding.**
  Per-pose geometry for 200 poses × 5,765 molecules is ~1.1 M rows — trivial.
  Coordinates are the expensive part (~30 heavy atoms × 200 × 5,765 × 4 bytes ≈
  1.4 GB as float32), which is affordable but should be a `.npy` per molecule
  rather than a column in a CSV.
- **Write coordinates for a mode-diverse subset, geometry for everything.** The
  geometry table is what scoring and clustering need; full coordinates are needed
  only for the poses that might get elevated.
- **This lands in `append_only/`**, so the old frames stay. That is the rule, and
  it means the change is non-destructive by construction.

**Accepted when:** a re-run of the 15 ground-truth molecules reproduces #30's
93.3% containment *from the persisted files alone*, with no docking.

---

### 1.2 Pose splitting — the core feature

**@tt8804, 2026-08-07 — the shape of this is settled, and it is simpler than what
I first proposed:**

> *"within the 200 to split/group poses by consensus and we treat all those split
> poses as separate candidates and we filter by our current methods"*

**A mode is a candidate row.** Not a nested structure hanging off a molecule — a
row in the same table, scored by the same functions, filtered by the same
pipeline.

```
   BEFORE                          AFTER
   t4_abc123   score 0.31          t4_abc123_m1   n=64   score 0.71
                                   t4_abc123_m2   n=51   score 0.12
                                   t4_abc123_m3   n=18   score 0.44
```

**This is the right call and it is worth saying why.** My version made a molecule
own a ranked set of hypotheses, which meant every downstream consumer — ranking,
class stratification, selection, the GUI, the elevation queue — needed to learn a
new nested shape. Flattening a mode into a row means **none of them change**. The
candidate table gets longer and everything that reads it keeps working. The
5,765-molecule table becomes perhaps 12,000–20,000 rows, which is nothing.

It also makes modes compete *directly*: a strong mode of a mediocre molecule can
outrank a weak mode of a good one, which is exactly the comparison the current
design cannot express.

**The stage sequence, confirmed by @tt8804:**

```
1. chemical-space generation
2. pose generation          dock 200, keep all 200
3. POSE SPLITTING           ← new stage: group by consensus, emit one row per mode
4. ranking                  unchanged methods, longer table
5. selection
6. elevation
```

Splitting is a **stage**, not a scoring option. It sits between generation and
ranking, and its output is the candidate table ranking already consumes. That is
why nothing downstream needs rewriting.

**"By consensus" now means something different, and better.** Consensus today is
*"do the top-10 by energy agree?"* — which #23/#30 killed, because the top-10 is
an energy-selected sample and energy is uninformative. Under splitting, consensus
becomes **how populated a mode is over all 200 poses**: 64 of 200 poses agreeing
on a geometry *is* the consensus for that mode. Same word, no energy anywhere in
it, and it stops being a whole-molecule scalar that punishes genuine multi-mode
binders.

**Design constraints, each from something measured:**

| constraint | why |
|---|---|
| cluster on **reactive-atom position + approach vector**, not whole-molecule RMSD | D0062 — two poses that place the warhead identically and differ in a distal ring are one mode for a covalent question |
| clustering **must not use docking energy** | #23/#30 — energy places the correct pose at a rank indistinguishable from uniform (KS p = 0.666). Using it to define or order modes re-imports the defect this version exists to remove |
| **mode count is measured, not a parameter** | fixing *k* decides the answer in advance. Density/agreement criterion, count reported with its stability across re-runs |
| a mode below a **minimum population is labelled noise**, not dropped | same rule the consensus floor follows |
| modes get **identities that survive a re-dock** — geometric, not rank-based | a mode named by rank is a mode a re-run silently redefines |

**Considerations that are not in the outline:**

- **Docking is stochastic, and mode identity has to survive that.** Re-running
  the 15 cases reproduced containment exactly (93.3% both times) but moved
  `KEEP_TOP`-survival from 46.7% → 40.0% and top-10 from 33.3% → 26.7%. If mode
  *membership* is that unstable, mode-level scores inherit the instability. The
  stability check is therefore not optional garnish — it is the acceptance test.
- **The 150° SN2 threshold interacts with clustering — settled in §3.1.** Sulfopin's
  crystal mode peaks at 146.4°, so a hard cut calls the crystallographic mode of a
  nanomolar inhibitor dead by 3.6°. If viability were a per-mode gate, one
  unjustified number would silently kill whole modes. **The angular criterion
  therefore becomes a graded term inside the mode score and gates nothing.**
- **What happens when a molecule has one mode?** Most will. The output shape must
  not privilege multi-mode molecules, or ranking acquires a bias toward
  promiscuous binders — the exact failure mode consensus had, inverted.

- **`topn_viable_frac` cannot survive into a mode.** It is defined as "the
  fraction of the top-10 *by energy* that are reaction-competent". Inside a mode
  that is still an energy-ordered sample of the mode's members. It is replaced by
  the mode's own `viable_fraction` over **all** its poses — no energy anywhere in
  the score. This retires the last consumer of the energy ordering.

- **⚠ THE TRAP: the deliverable is molecules, not modes.** Selection takes the top
  *n* per warhead class. Once modes are rows, one molecule contributing four modes
  can occupy four slots — and the top-5 becomes five modes of two compounds. The
  chemists synthesise **one compound a week** and the deliverable is a **top-5 of
  distinct molecules**, so a top-5 that is really a top-2 is a broken deliverable
  even though every row is legitimate.

  **The rule:** modes compete as rows through ranking, and **collapse to distinct
  parent molecules at selection.** A molecule enters the top-*n* once, on its best
  mode. Its other modes are retained and reported — because "this molecule has a
  second credible mode" is exactly the thing worth knowing at elevation, and
  testing two modes of one molecule is a stated goal — but they do not consume
  deliverable slots.

- **Class stratification survives untouched.** Each mode inherits its parent's
  warhead class, so `class_rank` and per-class selection (D0073 /
  `docs/class_stratification.md`) keep working on the longer table with no change.

**Accepted when:** (a) mode assignments are stable across two independent
re-docks of the same molecule; (b) on the 15 ground-truth cases, the crystal pose
falls inside a named mode ≥90% of the time; (c) elevating the top mode beats
33.3% (the current top-10 number), with 93.3% as the ceiling; (d) a top-5 request
returns **5 distinct molecules**.

---

### 1.3 The score — conditional enrichment and anchor quality, per mode

**Change.** Primary score becomes **conditional enrichment** — P(angle viable |
distance in window) — computed **per mode**. `anchor_quality` recomputed over a
mode's full population rather than averaged over the top-20 by energy.

**Considerations.**

- **`enrichment_conditional` is the one thing computable today**, from persisted
  aggregate counts. It can and should be tested before the re-dock.
- **The current `WEIGHTS = {anchor_quality: 0.5, topn_viable_frac: 0.5}` is a
  placeholder, not a finding.** #30 showed `topn_viable_frac` reads a window that
  holds the crystal pose 26.7–33.3% of the time. **Half the current weight rests
  on a mostly-wrong window.** Re-derive after the re-dock; do not carry the 0.5/0.5
  forward as if it meant something.
- **Test against Sulfopin first, not last.** Any candidate score that gives the
  parent compound a zero is wrong. Its crystal pose is produced at 1.00 Å at
  energy rank 160/200 — so a correct score *can* find it, and there is no excuse
  for a repeat.
- **Convergence check before it ranks anything** — 200 vs 2,000 runs, per D0068.
  No score in this project has passed one yet.

**Accepted when:** Sulfopin scores non-zero; the score converges between 200 and
2,000 runs; and it is not rank-correlated with docking energy (if it is, it has
re-imported the defect).

---

### 1.4 PoseBusters as a validity gate

**Change.** Installed and unused today. Gate poses **before** clustering.

**Consideration — this ordering is the whole point.** A reproducibly-invalid pose
is *reproducible*, so it will form a confident, well-populated cluster and look
exactly like a real mode. Gating after clustering would let physical nonsense
become a named hypothesis with a good score.

**Accepted when:** the fraction gated is reported per warhead class, and no class
loses so much that its ranking becomes unpopulated.

---

### 1.5 D0074 — the reactive path protonates at pH 7.4

**Change.** `nac_screen.prepare_ligand` adopts `obabel -p 7.4`, matching the
non-covalent path. **Accepted 2026-08-07; not yet implemented.**

**Considerations.**

- **594 of 1,782 T₄ (33.3%) and 331 of 5,370 T₃ are currently docked as the wrong
  species.** Every score for those describes a molecule that does not exist at
  pH 7.4. This is a correctness blocker, not a refinement.
- **It rides the same re-dock as pose splitting.** Separately, the library gets
  docked twice.
- **`charge_ph74` becomes the single charge annotation** once both paths agree —
  the reason for the divergence disappears with it.
- **The T₃ bidirectional disagreement is a separate defect** (142 rows where
  `charge_ph74` is *lower* than the SMILES charge, which protonation cannot
  cause). Not fixed by this and must not be assumed fixed.
- **Two 100 ns runs launched with `--net-charge +1`** (`t4_da2e98512d02`,
  `t4_9a973be6b946`) simulate the unprotonated species. Internally consistent,
  but the species this decision rules against. They must be redone.

---

### 1.6 Boltz-2 — adopt as a cross-check, hold on triage

**Measured (pre-registered in `prereg_cofolding.md`, cutoff verified 2023-06-01):**

| | n | within 2 Å | median |
|---|---:|---:|---:|
| **held out** | 10 | **60.0%** | 1.93 Å |
| in training (control) | 5 | 80.0% | 1.36 Å |
| **our docking, same 10** | 10 | **10.0%** | 2.82 Å |

**Change.** Adopt as an **independent cross-check on every elevated pose** — it
never sees our poses, so it is the only orthogonal signal we have.

**Considerations.**

- **Non-covalent is the right register, not a limitation — @tt8804, 2026-08-07.**
  I repeatedly flagged "Boltz-2 cannot model the covalent bond" as though it were
  a defect specific to co-folding. It is not, and singling it out was wrong.
  **The entire pipeline through MD is non-covalent**: the docking is
  reactive-*biased* (`rec.reactive_config`) but produces non-covalent poses, the
  NAC criterion is a geometric proxy measured on those poses, and
  `md_residence_3ikd` parameterises through `mmgbsa_noncovalent`. The covalent
  topology has existed since 2.0.0 and has never been run. **Covalent chemistry
  enters at FEP, on the few candidates that earn it.** Boltz-2 predicting the
  pre-reaction complex is therefore predicting exactly the same thing every other
  stage predicts — it is in register with the pipeline, not short of it.
- **The inferential trap is real but general, and belongs to the whole funnel.**
  "It binds, therefore it reacts" is unsupported — by co-folding, by our docking,
  by the NAC geometry and by 100 ns residence alike. That caveat attaches to the
  *stage*, not to any one tool, and it is discharged at FEP.
- **Its affinity head is largely pose-independent** (established in 2.1.0). Use
  structure and confidence; never affinity. *This* is a Boltz-2-specific limit.
- **The MSA cache is the premise, not an optimisation.** One MSA per construct
  makes a 300-molecule triage a forward pass each. Without it the cost argument
  collapses.
- **The triage use (T4) is unproven** — it needs the bornite runs. Do not build
  selection on it yet.

---

### 1.7 Carried-forward fixes

| fix | why it matters now |
|---|---|
| **`mmgbsa.RECEPTOR_PDB` required, not defaulted** | still defaults to **6VAJ**, which is 48.6 Å from the 3IKD pocket, and every covalent path takes it |
| **A chain stage that produces no output stops the run** | the overnight chain logged `exit 1` twice and continued, turning a good ranking into zero elevation |
| **Flexible Cys113 sidechain** | the anchor distance is measured to one arbitrary rotamer of the residue the criterion is about |
| **Covalent MD** | topology built and verified since 2.0.0, never run |
| **Manual elevate button** (#22) | selection is automatic-only; no way to push a molecule onto the queue by hand |
| **Mode-aware GUI** | once a molecule has modes, the viewer and tables must show them, or the feature is invisible to the person using it |

---

## 2. Revised order of work

The outline's order assumed no re-dock. This one does not.

| # | step | GPU | blocks |
|---|---|---|---|
| 1 | **Design + validate the keep rule** on the 15 ground-truth cases (200 poses each, already dumped) | none | everything |
| 2 | **`enrichment_conditional` + Sulfopin test + convergence**, on existing aggregates | none | — |
| 3 | **PoseBusters gate** wired in | none | 4 |
| 4 | **ONE re-dock**: D0074 protonation + new persistence + new keep rule | ~3 h | 5–7 |
| 5 | **Pose splitting** on the new pose sets; stability check | light | 6 |
| 6 | **Re-rank per mode**; re-derive weights | none | 7 |
| 7 | **Re-select and re-elevate** | heavy | — |
| — | *in parallel:* MD-priority verdict, chain robustness, receptor default, GUI | — | — |

**Steps 1–3 genuinely need no GPU**, which is the real dividend — and step 1 is
already running.

---

## 3. Decisions — @tt8804, 2026-08-07

**The re-dock does not increase sampling.** Worth stating plainly because it is
easy to read the plan as proposing more work: production *already* runs
`--nrun 200`. The change is that we stop deleting 180 of the 200 results. Same
docking, different write.

```
today:     generate 200 → sort by energy → keep 20 → delete 180 → score the 20
proposed:  generate 200 → keep all 200's geometry → split into modes → score each
```

| # | question | ruling |
|---|---|---|
| 1 | re-dock now or after the MD verdict | **after** — the bornite runs hold 5 cards and a 3 GPU-hour job is not worth delaying the verdict for |
| 2 | the 150° SN2 threshold | **not escalating — too in the weeds.** Resolved in design instead: see below |
| 3 | N-activated acrylamides, 97% of T₃ | **irrelevant — T₃ stays in the re-dock.** Question dropped |
| 4 | FEP licensing | **not now.** Nothing found yet that justifies FEP. Stays in the plan as a terminal, optional step; revisit when there is a molecule worth it |

### 3.1 The 150° threshold, resolved in design rather than escalated

Since this is not a question to put to a chemist mid-build, it is settled by
**not making it a decision at all**: the angular criterion becomes a **graded term
inside the mode score**, never a binary gate on a mode's members.

The reason is that a hard cut has to be *exactly* right or it silently destroys
information, and we already know it is not exactly right — it calls Sulfopin's
crystallographically-confirmed mode dead by 3.6°. A graded term degrades smoothly
instead: a 146° pose scores slightly below a 155° pose rather than scoring zero.
Nothing is thrown away on the strength of a number nobody has justified.

The binary `viable` flag stays in the output for continuity with every earlier
measurement, but it **stops gating anything**.

### 3.2 How many poses for 95% coverage? — **answer: ~300. Recommend 500.**

Measured, not assumed. All 15 ground-truth molecules were docked at `--nrun 2000`
and every pose scored against its crystal structure, which gives each molecule's
**hit rate** — the fraction of docking runs that land within 2 Å.

| | hit rate | 1 correct pose every |
|---|---:|---:|
| 9INN (16 heavy atoms) | 22.30% | 4 runs |
| 6VAJ / sulfopin (17) | 15.55% | 6 runs |
| *median of the 15* | *5.70%* | *18 runs* |
| 7F0M (31) | 1.70% | 59 runs |
| 9INP (23) | 1.05% | 95 runs |

**Difficulty is driven by molecular size, strongly and exponentially.**

- ρ(heavy atoms, hit rate) = **−0.683**, p = 0.005
- ρ(rotatable bonds, hit rate) = −0.587, p = 0.021
- fit: **each additional heavy atom multiplies the hit rate by 0.883** — a 12%
  penalty per atom, compounding
- molecule-to-molecule spread at fixed size: ~1.9×

**This is why the crystal benchmark cannot be read off directly.** Those 15 have a
median of 22 heavy atoms. Our library is bigger: T₃ median 25, **T₄ median 28**,
90th percentile 34. We sit at the hard end of the fitted range and partly beyond
it.

Runs needed so a molecule is covered 95% of the time (`n = ln0.05 / ln(1−p)`):

| heavy atoms | predicted hit rate | n for 95% | n if unlucky (−1 sd) |
|---:|---:|---:|---:|
| 22 (benchmark median) | 5.79% | 50 | 95 |
| 25 (T₃ median) | 3.99% | 73 | 139 |
| **28 (T₄ median)** | **2.75%** | **107** | 202 |
| 32 | 1.68% | 177 | 333 |
| 35 | 1.16% | 258 | 484 |

Covering **95% of the pool** means covering the 5th-percentile-hardest molecule,
not the median: **T₃ needs 227 runs, T₄ needs 258.**

**This agrees with what we observed directly.** Two independent 200-run dockings
of the 15 both gave 93.3% — just under 95%, exactly where ~250–300 being the right
number would put them.

| | |
|---|---|
| **~300 runs** | the point estimate for 95% coverage of our library |
| **500 runs — recommended** | margin for an extrapolation fitted on 15 molecules and pushed past its size range. ≈7.5 GPU-h for the library |
| 2000 runs | reaches 100%, but 10× the cost for the last few percent |
| 200 runs (today) | **marginally short**, and never justified — inherited |

**Recommend 500.** The library is docked once; the difference between 300 and 500
is ~3 GPU-hours, and discovering afterwards that 300 was short means docking it
all again.

**What more sampling does NOT fix.** At 2000 runs the crystal pose was present for
every molecule — and the energy cut still discarded 60% of them. Sampling raises
the ceiling; only the keep rule (§1.2) reaches it.

### 3.3 Original convergence question — resolved by the above

Your question about pose counts exposed one nobody had asked: **200 is inherited,
not justified.** D0068 requires every number to carry the parameter that defines
it, and this one does not.

#30 makes it directly testable, so it is running: the same 15 ground-truth
molecules at `--nrun 2000`. The readout is whether containment rises above the
93.3% measured at 200.

| observation | conclusion |
|---|---|
| containment materially above 93.3% | 200 under-samples; the re-dock uses more runs, and the cost estimate rises with it |
| containment ≈ 93.3% | 200 is sufficient and now **justified rather than inherited**; re-dock proceeds at 200, ~3 GPU-h |
| mode structure changes with sampling depth | more serious than either — "mode" would be partly a sampling artefact, and §1.2's stability test has to cover sampling depth as well as re-runs |

Answering this **before** the re-dock costs under an hour on one card and stops us
committing 5,765 molecules to an unjustified parameter.

---

## 4. What would make 2.2.0 a failure

Stated in advance:

- **A score that ranks well and still gives Sulfopin a zero.**
- **Mode counts that change between re-runs** — then a "mode" is an artefact of
  the clustering, not a property of the ligand.
- **Elevating a mode BPMD then shows was the wrong one, repeatedly** — would mean
  modes are ranked on something as uninformative as energy was.
- **Another silent stage.** If a 2.2.0 run can produce zero output and report
  success, nothing else here matters.
- **Shipping the re-dock without D0074** — a third of T₄ would stay the wrong
  molecule, and every downstream number would inherit it.
