# Which module does what

*The four stages — generation → ranking → selection → elevation — and the module
that owns each. 86 scripts live in `scripts/`; this says which are the pipeline,
which measure it, and which are finished business. 2026-08-07.*

---

## The stages at a glance

| # | stage | question it answers | owns |
|---|---|---|---|
| 1 | **generation** | what molecules exist? | enumeration from the T₂ seeds |
| 2 | **ranking** | in what order should we look at them? | docking, pose splitting, scoring |
| 3 | **selection** | which few do we spend money on? | quotas, stratification, the queue |
| 4 | **elevation** | do these actually work? | the validation ladder, cheap → expensive |

A stage boundary is where the **unit changes**. Generation emits molecules;
ranking emits *scored binding modes*; selection emits a *queue*; elevation emits
*evidence about individual molecules*.

---

## 1. Generation — molecules

| module | role |
|---|---|
| `sample_t2_degree2.py` | uniform random sample of T₂'s degree-2 neighbourhood — **a sample, not a truncation** |
| `audit_degree2.py` | checks that it *is* a sample: uniform, complete, reproducible |
| *(T₃ / T₄ frames)* | REINVENT and combinatorial arms, consumed as `D3_*.parquet` / `D4_*.parquet` |

**Emits:** ~5,769 molecules with SMILES and a warhead class.

**Not the constraint.** Enumeration delivered 4,082 T₃ + 1,683 T₄; the cost is
entirely downstream (#15). More molecules is not the bottleneck.

---

## 2. Ranking — scored binding modes

This is where 2.2.0 changed most. Three steps inside one stage:

| step | module | what it does |
|---|---|---|
| 2a. pose generation | `nac_screen_v2.py` | 500 independent reactive dockings per molecule, ligand protonated at pH 7.4 (D0074) |
| 2b. **pose splitting** | `shared/pose_modes.py` | clusters poses into binding modes on the **reactive atom's position and warhead direction** |
| 2c. scoring | `rank_v2.py --topic nac_v3` | conditional enrichment, anchoring, consensus — **per mode** |

**Supporting:**

| module | role |
|---|---|
| `nac_rank.py` | builds the candidate pool from the frames (`load_candidates`) |
| `shared/nac_criterion.py` | the near-attack criterion **and** `anchor_quality` — one definition, two consumers |
| `shared/warhead_library.py` | mechanism and reactive SMARTS per class; never a per-frame column |
| `shared/ionisation.py` | the single pH 7.4 protonation both arms now use |

**Emits:** one row per **binding mode**, ident `<parent>_m<k>`, into `nac_v3`.
~2.1 modes per molecule, so ~12k rows from ~5.7k molecules.

**Key invariants:**
- **Energy generates the poses and never selects them.** Each run optimises
  AutoDock's function, so energy decides where a run lands; nothing downstream
  reads it. `mean_energy` is reported and consumed by nothing.
- **Consensus = mode population** over all 500 poses. Not "do the top-10 by
  energy agree".
- **Class stratification lives here** (`class_rank`), which is what contains
  D0073.

**Superseded:** `rank_2_1.py`, `nac_screen.py`'s own `main` (its helpers are
still used), `nac_stage4.py`, `topn_vs_fraction.py` — all 2.0.0/2.1.0 scoring
experiments, kept for provenance.

---

## 3. Selection — a queue

| module | role |
|---|---|
| `select_elevate.py` | walks the ranked list, **re-measures the pose it is about to elevate**, writes one queue |
| `cluster_poses.py` | picks a real representative pose (superseded by `pose_modes` inside the screen) |

**Emits:** an elevation queue — molecule, mode, pose, and why it was chosen.

**Two rules that live here and nowhere else:**

1. **Per-class quotas.** Top *n* within each warhead class, never a global top-*n*.
   This is what makes D0073's depletion unable to reach the elevated set.
2. **Collapse to distinct parent molecules.** Modes compete as rows through
   ranking; at selection a molecule enters **once**, on its best mode. You
   synthesise a molecule, not a mode, so a top-5 must be five compounds.
   ⚠️ **NOT YET IMPLEMENTED** — `select_elevate` has no `parent_ident` handling.

**Selection never re-scores.** It applies quotas to an existing order and
re-measures geometry as a guard. If selection is choosing, ranking has failed.

---

## 4. Elevation — evidence, cheap before expensive

The ladder, with measured cost per molecule:

| rung | module | cost | asks |
|---|---|---|---|
| 4a. pose confirmation | `cofold_bench.py` (Boltz-2) | ~50 s | does an independent method put it here? |
| 4b. physical validity | **PoseBusters** | seconds | is this pose even physically possible? ⚠️ **NOT BUILT** |
| 4c. attack sweep | `attack_sweep.py` | ~0.4 GPU-h | does the warhead reach and hold attack geometry? |
| 4d. pose metadynamics | `bpmd_run.py` | ~1 GPU-h | does the pose survive a biased perturbation? |
| 4e. residence | `md_residence_3ikd.py` | ~4 GPU-h | does it stay for 100 ns, and stay competent? |
| 4f. covalent workup | `covalent_workup_one.py` | — | adduct, covalent docking, MM-GBSA |
| 4g. developability | `medchem_workup.py`, `medchem_admet.py`, `medchem_retro_report.py` | — | is it a drug, and can it be made? |
| 4h. free energy | **FEP** | — | deferred; the only covalent-aware rung |

**Drivers and reporting:**

| module | role |
|---|---|
| `elevate_queue.py` | launches the suite over a selection queue |
| `elevate_reference.py` | puts Sulfopin/ATRA through the identical path — the yardstick |
| `mdprio_report.py` | one molecule → one report, the moment it lands |
| `elevation_report.py` | the 2.1.0 cohort report (tier-1/tier-2 design) |

**The ordering inside elevation is not arbitrary** — realism before geometry:

```
consensus (where it sits) → Boltz-2 (independent support) → PoseBusters (valid?)
                                  → attack geometry (RANKS what survived)
```

Attack geometry never *selects* a pose; it orders a set already established as
realistic. Reversing that selects strained poses that happen to point the warhead
— measured, and the reason the representative is now a top-quartile medoid rather
than an argmax.

---

## Not the pipeline — the instruments that measure it

These never run in production. They exist to tell us whether a stage works, and
they are why most of this document's numbers exist.

| module | measures |
|---|---|
| `pose_split_validation.py` | is the crystal pose **in** our pose set? (93.3% at 200, 100% at 500) |
| `mode_arbitration.py` | can Boltz-2 pick the right mode? (no — consensus is at ceiling) |
| `attack_sweep_check.py` | does 10 ns predict 100 ns? (ρ = +0.83) |
| `redock_3ikd_benchmark.py` | pose recovery on 82 crystal cases |
| `cofold_docking_comparator.py` | our docking on the co-folding benchmark's own ligands |
| `pose_selection_bench.py` | does any cheap pose rule beat random? |
| `nac_decoy_site.py` | does the criterion measure **Cys113** or just "a warhead can point somewhere"? |
| `nac_robustness.py` | is the AUC real or an artefact of which negatives were drawn? |
| `consensus_convergence.py` | does consensus converge where viable-NAC frequency did not? |
| `screen_references.py` | known binders through the **identical** criterion |

---

## Where the current work sits

| | stage | status |
|---|---|---|
| library re-dock, 500 runs | **2a** | running, ETA ~22:30 |
| automated ranking chain | **2c** | armed, fires on completion |
| collapse-to-parent | **3** | not built |
| PoseBusters | **4b** | not built — and it must gate *before* 2b clustering, so it needs a re-dock |
| `attack_sweep.py` | **4c** | built, never run end to end |
| 6th bornite molecule | **4e** | ~19:15 |

## Two boundary rules worth keeping

**Ranking scores; selection allocates.** If a quota decision appears in
`rank_v2`, or a scoring decision in `select_elevate`, the boundary has been
crossed and the two will drift.

**Elevation never re-orders.** It produces evidence about individual molecules.
If elevation results start feeding back into the order, that is a new ranking
input and belongs in stage 2 with a pre-registration — not an ad-hoc re-sort.
