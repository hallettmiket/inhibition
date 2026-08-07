# The 2.2.0 framework — how a molecule gets to synthesis

*Canonical description of the pipeline as it now stands, with every claim
measured. Written from [#31](https://github.com/hallettmiket/inhibition/issues/31)
(@tt8804's framing) plus the measurements of 2026-08-07. Successor to
[`framework_2.1.0.md`](framework_2.1.0.md); the build itself is in
[`build_plan_2.2.0.md`](build_plan_2.2.0.md).*

---

## 0. The problem in one paragraph

We have ~5,769 generated molecules and want to know which few are worth
synthesising, at one compound a week. Every stage before FEP is a **cheap,
non-covalent proxy** for a covalent event, arranged so that the expensive
measurements are only spent on molecules that survived the cheap ones. The
deliverable is a **top-5 of distinct molecules**, ranked, with the whole list
kept queryable — nothing is deleted.

## 1. The pipeline

```
chemical-space generation
        ↓
pose generation            500 independent runs        ~4 h wall-clock, whole library
        ↓
POSE SPLITTING             consensus → one candidate ROW per binding mode
        ↓
ranking                    anchoring × consensus, per mode, stratified by warhead class
        ↓                  collapse to distinct parent molecules
Boltz-2 pose confirmation  ~50 s/molecule              SHORTLIST ONLY
        ↓                  gate: does the prediction fall inside this mode?
[BPMD]                     ~1 GPU-h                    kept only if it earns it
        ↓
100 ns MD                  residence AND attack geometry   ~4 GPU-h each
        ↓
FEP                        covalent, on the few that earn it
        ↓
synthesis
```

**Everything above FEP is non-covalent.** The docking is reactive-*biased*
(`rec.reactive_config`) but yields non-covalent poses; the NAC criterion is a
geometric proxy computed on them; `md_residence_3ikd` parameterises through
`mmgbsa_noncovalent`. The covalent topology has existed since 2.0.0 and has never
been run. Covalent chemistry enters once, late, on the shortlist.

The consequence is a **stage-level** caveat that belongs to the whole funnel and
to no single tool: *"it binds, therefore it reacts"* is unsupported by docking,
by the NAC geometry, by co-folding and by 100 ns residence alike. It is
discharged at FEP.

---

## 2. Pose generation — the tool is fine, the selector was not

This is the correction that reorganised the version.

| runs | correct pose present in the set |
|---:|---|
| 200 | 93.3% of molecules |
| **500** | **100%** |

The generator finds the right answer for every molecule. It is *stochastic* — a
single run has a ~5.7% median hit rate — but 500 independent searches are not
"the 500 best" of a larger pool; each is one search's own answer, which is why
hit rates compose the way they do.

**What failed was the step immediately after.** Sorting those poses by docking
energy and keeping 20:

| how 20 poses are chosen | crystal pose within 2 Å |
|---|---:|
| all 500 (ceiling) | 100.0% |
| **mode 0 — the new rule** | **93.3%** |
| naive densest 20, no clustering | 73.3% |
| **random 20** | **80.0%** |
| lowest-energy 20 — *what 2.1.0 did* | **60.0%** |

Energy is **worse than random** at this. It places the correct pose at a rank
indistinguishable from uniform — KS *p* = 0.666 against uniform(1,200) — and
Sulfopin's crystal pose is produced at 1.00 Å at **energy rank 160 of 200**
(#23, #30).

**Energy still generates the poses; it no longer selects them.** Each run is an
independent Lamarckian GA optimising AutoDock's scoring function, so energy
decides where a run *lands* — that cannot be removed without removing the
docking. What is removed is energy as a *selection* criterion. `mean_energy` is
reported per mode and consumed by nothing.

### Why 500

200 was inherited and never defended, which D0068 forbids. Per-run hit rate falls
exponentially with molecule size: ρ(heavy atoms, hit rate) = **−0.683**
(*p* = 0.005), each additional heavy atom multiplying it by **0.883**.

Covering 95% of the *pool* means covering the 5th-percentile-hardest molecule,
not the median: **227 runs for T₃, 258 for T₄**. Two independent 200-run dockings
both returned 93.3% — just under 95%, exactly where that estimate puts them. 500
is that plus margin for a fit built on 15 molecules and extrapolated past its
size range; the library is docked once.

**It costs ~25% more than 200, not 2.5×.** AutoDock-GPU runs its LGA instances
concurrently, so ~3.6 s per molecule is fixed cost (ligand prep, meeko rebuild,
gnina rescoring, I/O) and only ~0.0032 s is per-run. We were paying the fixed
cost 5,769 times and under-using a GPU that was already spun up.

### Persist the population, not the survivors

2.1.0 wrote per-pose rows for the top 20 by energy and discarded the other 180 at
write time. That made pose splitting impossible retrospectively — clustering the
kept 20 would cluster the survivors of the exact filter being removed. 2.2.0
persists every pose's geometry, and coordinates for the mode representatives.

---

## 3. Pose splitting — a mode is a candidate

@tt8804: *"within the 200 to split/group poses by consensus and we treat all those
split poses as separate candidates and we filter by our current methods."*

```
   BEFORE                          AFTER
   t4_abc123   score 0.31          t4_abc123_m0   n=303  consensus 0.606
                                   t4_abc123_m1   n=75   consensus 0.150
                                   t4_abc123_m2   n=32   consensus 0.064
```

A mode is a **row**, not a nested structure. Ranking, class stratification,
selection and the GUI all read rows, so the table simply gets longer and nothing
downstream was rewritten. It also lets modes compete directly — a strong mode of
a mediocre molecule can outrank a weak mode of a good one, which the 2.1.0 shape
could not express.

**What is clustered on:** the reactive atom's position and the direction its
warhead faces.

- **Not whole-molecule RMSD** — D0062: two poses that place the warhead
  identically and differ in a distal ring are one mode for a covalent question.
- **Not docking energy, at all** — using it to define or order modes re-imports
  the defect this version exists to remove.
- **Not the NAC geometry itself** — that is the score. A mode defined partly by
  its own score is guaranteed to look internally consistent.

**"Consensus" is redefined and improved.** It was *"do the top-10 by energy
agree?"*, which read an energy-selected sample of a uniformly uninformative
ordering. It is now **how populated a mode is over all 500 poses**. Same word, no
energy in it, and it stops penalising molecules with a genuine second mode.

`topn_viable_frac` does not survive: inside a mode it is still an energy-ordered
sample. Replaced by the mode's own `viable_fraction` over all its poses — which
retires the last consumer of the energy ordering.

**Both parameters are calibrated, not chosen (D0068):**

| `eps` | modes/mol | crystal in a named mode | crystal in the top mode |
|---:|---:|---:|---:|
| 1.0 | 5.5 | 60% | 33% |
| **3.0** | **2.9** | **100%** | **93%** |
| 4.0 | 1.5 | 100% | 100% |

4.0–5.0 looks perfect only because it degenerates to ~1 mode, making "the crystal
is in the top mode" a tautology. `min_population_frac = 0.05` likewise: at 0.02
the mode *count* reproduces only 47% across independent dockings.

**Measured on the library:** ~2.1 modes per molecule, range 1–8.

**Honest limit.** The dominant mode reproduces across independent dockings 86.7%
of the time, the mode *count* only 73%. **Ranking the dominant mode is supported;
treating every minor mode as an independent candidate is not yet.**

---

## 4. Ranking — consensus works, anchoring-within-a-mode does not (yet)

**Consensus picks the right mode 93.3% of the time, which *is* the ceiling.**
Nothing improves on it at this n, including Boltz-2.

**Picking one pose out of that mode is unsolved:**

| rule | within 2 Å |
|---|---:|
| ceiling — best pose in the mode | 93.3% |
| medoid | 26.7% |
| lowest energy in the mode | 20.0% |
| best `anchor_quality` | 6.7% |

**That last number is measured against the wrong target and must not be read as
"anchoring is broken."** The deposited ligand is a *post-reaction adduct* whose
leaving group is gone, so RMSD-to-crystal matches 16 of 17 atoms and leaves the
leaving-group direction **unconstrained** — while the SN2 angle depends entirely
on it. Measured: the crystal-matching pose sits in the NAC distance window 80% of
the time but clears the 150° criterion for **1 of 9** chloroacetamides.

An adduct records where a ligand *ended up*, not the trajectory it took. So the
crystal set arbitrates **where a molecule binds** and cannot arbitrate **which
pose reacts**. MD is the only arbiter we have for the latter.

### The score

Primary is **conditional enrichment** — P(angle viable | distance in window) —
which uses no ordering and so cannot inherit a bad one. On the 2.1.0 data it
already passes the test the old default fails catastrophically:

| | Sulfopin | percentile of 5,765 |
|---|---:|---:|
| `topn_viable_frac` (2.1.0 default) | **0.000** | **0.0%** |
| `enrichment_conditional` | 2.590 | **78.1%** |

The score we had been ranking on put the crystallographically-confirmed nanomolar
parent **dead last in the library**. Juglone — the known promiscuous
naphthoquinone — scores 0.000 on conditional enrichment, which is the right
direction.

**The `0.5/0.5` weights are a placeholder, not a finding**, and half that weight
is `topn_viable_frac`. Re-derive after the re-dock; do not carry them forward.

### The angular criterion gates nothing

The 150° SN2 threshold becomes a **graded term inside the mode score**. A hard cut
has to be exactly right or it destroys information silently, and it is not: it
calls Sulfopin's crystallographically-confirmed mode dead by **3.6°**. Graded, a
146° pose scores just below a 155° one instead of scoring zero. The binary
`viable` flag stays for continuity with earlier measurements and stops gating.

### Class stratification

Ranking and selection are both **within warhead class** — `class_rank` in
`rank_v2`, `nsmallest(per_class)` in `select_elevate`. A molecule competes only
against others carrying the same warhead.

This is what contains D0073: consensus *depletes* validated chemistry (library
90.3% validated, pool 77.8%, Fisher OR 0.34, *p* = 2.3×10⁻¹⁴), and BDHI outscores
the validated chemistry on anchoring too (Cliff's *d* = +0.632). Under stratified
selection the mechanism mix of the elevated set is fixed by `--per-class`, not by
the score, so a composition effect of any size cannot change what gets elevated.
Full argument in [`class_stratification.md`](class_stratification.md).

### The deliverable collapses to molecules

Once modes are rows, one molecule with four modes can occupy four slots in a
top-5. Every row would be legitimate, and the deliverable would be broken —
**you synthesise a molecule, not a mode.**

**Rule:** modes compete as rows through ranking, then **collapse to distinct
parent molecules at selection**. A molecule enters the top-*n* once, on its best
mode. Its other modes are retained and reported — "this one has a second credible
mode" is worth knowing at elevation — but they do not consume deliverable slots.

---

## 5. Boltz-2 — a pose source, not an arbiter

Pre-registered in [`prereg_cofolding.md`](prereg_cofolding.md) before installation.
Training cutoff verified at **2023-06-01**, then deposition dates read out of the
mmCIF files: all 10 held-out entries 2024-07→2025-05, all 5 controls 2019→2021.
No entry changed sides.

| | n | within 2 Å | median |
|---|---:|---:|---:|
| **held out** | 10 | **60.0%** | 1.93 Å |
| in training (control) | 5 | 80.0% | 1.36 Å |
| **our docking, same 10 molecules** | 10 | **10.0%** | 2.82 Å |

The control behaves as required, so T1 is interpretable. The comparator had to be
*measured*: four held-out entries appear in the 82-case benchmark, but every one
of those cases docked `A1ERA` — a second ligand present in all four — not the
covalent inhibitor. Overlap with the ten under test was **zero**.

**As a mode arbiter it fails** — 86.7% against consensus's 93.3%, and 77.8%
against 88.9% on the multi-mode subset where it would matter. Consensus is at
ceiling. Dropped.

**As a pose source it succeeds:** its prediction lands **inside** the consensus
mode for **100%** of molecules (median 0.89 Å), and its pose is within 2 Å of the
crystal **67%** of the time against our medoid's **27%**.

**So consensus picks the mode, Boltz-2 picks the pose, MD decides.**

**Placement:** ~50 s/molecule with the MSA cached is **~80 GPU-hours across
5,769** — more than the docking it validates — and **~1.4 GPU-h across ~100**. It
is cheap *relative to MD*, not cheap absolutely, so it sits on the shortlist. Its
*affinity* head is largely pose-independent and is never used.

**The gate** — @tt8804's *"if it fails the transition we move down the list and
don't waste MD"* — is Boltz's prediction falling inside the mode, with
`dir_coherence` and `spread_a` as secondary confidence. **It never fired on the 15
benchmark molecules**, so its false-positive rate is unmeasured. Adopt as a
tripwire, log every trigger, measure the rate on real candidates.

### "Probabilistic pose" cannot mean averaging

| | median nearest-atom distance |
|---|---:|
| a real docked pose | 1.52 Å |
| the mode's coordinate average | **0.43 Å** |

A C–C bond is ~1.5 Å. Averaging many orientations of a flexible molecule
collapses it toward its centroid; the result is not a structure and cannot be
parameterised or simulated. A consensus pose must be an actual ensemble member or
an independently predicted structure. What the distribution legitimately gives is
**uncertainty** — `consensus`, `spread_a`, `dir_coherence` — not geometry.

---

## 6. MD — residence and attack geometry are different readings

Five of the six bornite runs:

| molecule | BPMD | residence | outcome | frames angle-competent |
|---|---:|---:|---|---:|
| `t4_da2e98512d02` | 0.365 | 0.791 | left 81 ns | — |
| **`t4_7e86b677bb2d`** | 0.189 | **1.000** | **held** | **0.8%** |
| `t4_9a973be6b946` | 0.161 | 0.774 | left 79 ns | — |
| **`t4_4e608398fd6a`** | 0.125 | **1.000** | **held** | **0.0%** |
| `t4_9265b4bff789` | 0.108 | 0.383 | left 39 ns | — |

**Both molecules that held perfectly are angle-competent in under 1% of frames.**
They sit in the pocket beautifully and essentially never present the warhead
correctly.

If that survives the sixth molecule, **residence measures something real but not
reaction competence**, and the funnel needs a second gate rather than a better
residence threshold. Note that the warhead→SG *distance and angle* had never been
plotted for a 100 ns trajectory before 2026-08-07 — the runs produce
`rmsd`/`mindist`/`numcont`, and `mindist` is the closest approach by *any* atom
pair, which is neither the warhead nor the sulfur.

**BPMD's fate is pre-registered.** It asks a physics question Boltz-2 cannot, so
they are not substitutes — but it is a ~1 GPU-hour *proxy* for a 4 GPU-hour
measurement, and a proxy that does not predict its target has no role. If
occupancy does not rank-correlate with residence, drop it and go straight from
the Boltz-2-confirmed pose to MD.

---

## 7. Species, cost, and the shared box

**D0074 — the reactive path protonates at pH 7.4.** It previously docked the
SMILES as drawn while the non-covalent path ran `obabel -p 7.4`, so the two arms
docked **different species of the same molecule** for **594 of 1,782 T₄ (33.3%)**
and 331 of 5,370 T₃. Now both use the same call. It refuses rather than falling
through to the unprotonated species, and refuses if the reactive SMARTS stops
matching after protonation. **It drops ~0.4% of the library** (~23 molecules)
where obabel emits SMILES RDKit cannot parse — deprotonated tetrazoles and
enolates. @tt8804: that is the filter doing its job.

**Cost, measured:**

| stage | cost |
|---|---|
| docking, whole library at 500 runs | ~4 h wall-clock on 3–4 shards |
| pose splitting | **~3 min**, one core — 28 ms/molecule |
| ranking | minutes |
| Boltz-2, shortlist of ~100 | ~1.4 GPU-h |
| 100 ns MD | ~4 GPU-h **each** |

The box is 8× A100-SXM4-80GB, 500 W cap, drawing ~69 W idle and ~190 W under MD.
A 100 ns run is ~0.6 kWh — about twelve cents. **The binding constraint is not
electricity, it is the ~35 people sharing the machine**: hence `nice -n 19`
throughout, GPUs 0 and 4 refused by name, and launchers that wait for a card to
be genuinely idle (both <10% utilisation *and* <1.5 GB resident — utilisation
alone reads 0% between kernels of a job that still owns the card).

---

## 8. Tooling, settled

| tool | status |
|---|---|
| **GNINA** | **adopted** — CNN rescoring lifted pose recovery 18.3% → 26.8% on this receptor |
| **Boltz-2** | **adopted as a shortlist pose source**, §5 |
| **Desmond** | **declined** — not on engine quality. Desmond is well regarded and BPMD's canonical implementation is Desmond's (Clark 2016). FEP+ for Academic Research *exists*, but excludes *"drug discovery or other IP generating activities"*, which is what this project is. A licensing decision, not a technical one |
| **FEP** | **deferred** — the right terminal gate; nothing found yet that justifies it |
| **PoseBusters** | installed, **still unused**. Must gate **before** clustering: a reproducibly-invalid pose is *reproducible*, so it forms a confident cluster and looks exactly like a real mode |

---

## 9. Open, and honest about it

- **Mode count reproduces only 73%.** Minor modes as independent candidates are
  not yet supported.
- **The `0.5/0.5` weights are a placeholder**, and half rests on a window
  measured as ~30% correct.
- **No score has passed a convergence check** (200 vs 2,000 runs, D0068).
- **`consensus_gnina` runs higher in T₄ acrylamides than T₃ acrylamides at
  matched flexibility** (*p* = 1.8×10⁻¹¹, #29). Unexplained, and a live ranking
  input.
- **Why unvalidated chemistry outscores validated chemistry** is unresolved —
  either BDHI genuinely forms better near-attack geometry, or `anchor_quality` is
  partly measuring conformational rigidity under another name. No measured BDHI
  activity exists to distinguish them, which is why the class is unvalidated, and
  that circularity is why the answer is stratification rather than re-weighting.
- **PoseBusters, the `mmgbsa` 6VAJ default, and covalent MD** are all still
  outstanding from the carried-forward list.

## 10. What would make 2.2.0 a failure

- A score that ranks well and still gives Sulfopin a zero.
- Mode counts that change between re-runs — then a "mode" is an artefact of the
  clustering, not a property of the ligand.
- Elevating a mode that MD repeatedly shows was the wrong one.
- **Another silent stage.** If a run can produce zero output and report success,
  nothing else here matters.
