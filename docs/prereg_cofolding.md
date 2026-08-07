# Pre-registration — is a co-folding model worth putting in this pipeline?

*Written and committed **before** Boltz-2 is installed or run. @tt8804,
2026-08-07. Companion to `docs/prereg_md_priority.md`.*

---

## Why bother at all

Every signal this project has tested — docking energy, enrichment, consensus,
top-N viability, anchor quality — is computed from **the same AutoDock pose set**.
They therefore share failure modes, and #23 showed the shared failure mode is
serious: the ordering they all sit on carries no information about reaction
geometry.

A co-folding model builds the complex from sequence and ligand alone. **It never
sees our poses.** Whatever it says is orthogonal, and orthogonal is the one thing
a 300-molecule triage does not currently have.

Two candidate uses:

1. **Pose recheck** — does an independent method put the ligand where we did?
2. **MD-priority triage** — does its confidence predict which molecules hold their
   pose over 100 ns? At **seconds per complex** against BPMD's **~1 GPU-hour**,
   this is the one that would change what we can afford.

## The contamination problem, and why this target escapes it

Co-folding models are trained on the PDB. Asking one to reproduce a crystal pose
it was trained on measures memorisation, not prediction — so most obvious tests
here are worthless.

**Our crystallographic positives split by deposition era:**

| era | PDB codes | n | status |
|---|---|---:|---|
| 2024+ | 9INN 9INO 9INP 9INQ 9JF6 9JFH 9V6G 9V6I 9V6P 9V6W | **10** | **almost certainly held out** — AF3's cutoff is 2021-09-30, Boltz's comparable |
| 2020–23 | 7EFJ 7EFX 7EKV 7F0M | 4 | likely in training |
| pre-2020 | 6VAJ (sulfopin) | 1 | in training |

So there is a **genuine held-out benchmark of 10 Pin1 complexes**, on our exact
target, spanning our exact chemistry (chloroacetamide, naphthoquinone, SNAr).
That is the rare thing that makes this testable rather than arguable.

**Deposition date is a proxy for exclusion, not proof.** Before any result is
quoted, the model's stated training cutoff will be checked and any 2024+ entry
that turns out to be included will be moved to the contaminated set.

## Tests, fixed now

### T1 — held-out pose accuracy *(primary)*

On the **10 held-out** complexes: fraction where the predicted ligand is within
**2 Å symmetry-corrected RMSD** of the crystal ligand.

**Our docking is measured on the same 10 molecules**, not compared against the
82-case benchmark. Different case sets are not comparable and quoting 18.3%
against a 10-case number would be exactly the sort of thing this project keeps
catching.

### T2 — contamination control

Same measurement on the **5 in-training** complexes. Accuracy there should be
**higher** than T1. If it is not, either the split is wrong or the pipeline is
broken, and T1 cannot be interpreted until that is resolved.

### T3 — mode arbitration on sulfopin

Sulfopin's docked poses form two modes (issue #26): sulfolane in the **proline**
pocket (crystal-correct, ranks 8+) or against the **basic cluster** (ranks 1–7).
Docking energy prefers the wrong one by 0.16 kcal/mol.

**Does the co-folded pose land in the proline mode?**

6VAJ is in training, so a *correct* answer proves little. A **wrong** answer is
informative and would rule out mode arbitration immediately. Asymmetric test,
recorded as such.

### T4 — MD-priority prediction *(the brief)*

On our **generated** molecules — uncontaminated by construction, since they do
not exist outside this project — does interface confidence predict 100 ns
residence?

Set: Sulfopin and ATRA (held), `t4_72f5671e89cb` (left at 54 ns), and the six
MD-priority molecules as their runs land. Note Sulfopin and ATRA are *not*
generated and carry the contamination caveat for the pose, though not for the
residence outcome.

**Primary metric:** ligand-interface confidence (`ligand_iptm`, or the closest
the model reports). **Secondary:** complex pLDDT and interface PAE.

## Readings, fixed in advance

| observation | conclusion |
|---|---|
| **T1 ≥ our docking + 15 points**, T2 higher still | Co-folding predicts poses better than our docking on held-out Pin1 chemistry. Adopt as a pose source for the non-covalent leg and as a cross-check on every elevated pose |
| **T1 within ±15 points of our docking** | No better at posing, but still *orthogonal*. Value rests entirely on T4; do not adopt as a pose source |
| **T1 ≪ our docking**, T2 high | Memorisation without generalisation. Useless here. Report and stop |
| **T4 separates held from left** | Adopt as the MD-priority triage. Three orders of magnitude cheaper than BPMD |
| **T4 does not separate** | Not a triage signal at this n. Fall back to BPMD occupancy; do not spend more on co-folding |
| **T3 lands in the basic mode** | Cannot arbitrate modes. Pose splitting must be settled by physics, not by co-folding |

## What this cannot settle, whatever it returns

- **Co-folding does not model the covalent bond.** Every prediction is of the
  *non-covalent* pre-reaction complex. It can never speak to whether a molecule
  reacts — only to whether it sits in the pocket. That is a hard limit, not a
  caveat to be worked around.
- **n = 10 held out** supports only large effects. A null means "not
  demonstrated at n = 10".
- **T4 will have n ≤ 9** and mixes generated molecules with two references. No
  p-value will be quoted from it.
- **Boltz-2's affinity head is largely pose-independent** (established during
  2.1.0). Its *structure* and *confidence* are usable for pose questions; its
  affinity is not, and will not be used as one.
- **A confidence score is not a stability measurement.** If T4 succeeds it means
  confidence *correlates with* residence on nine molecules, not that it measures
  it.

## Cost, and the one trick that makes it affordable

The protein is **always Pin1**. Its MSA is computed **once** and reused for every
ligand, so per-candidate cost collapses to a single forward pass. That is what
makes a 300+ molecule triage plausible; without it, MSA generation would dominate
and the whole argument for co-folding as a cheap filter fails.

Budget for this experiment: one MSA, ~25 predictions. Everything else is analysis
against ground truth we already hold.
