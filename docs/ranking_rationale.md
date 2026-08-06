# Ranking rationale

*Agreed with @tt8804, 2026-08-05. This is the reasoning behind how candidates are
ordered, and the measurements that force it. Read it before proposing a scoring
function.*

---

## The one sentence

**Rank on whether a molecule can orient to form the bond, not on how good the
bond would be.**

## Why

Covalent inhibition happens in two steps:

1. **Recognition** — the molecule binds non-covalently in an orientation that
   presents its warhead to Cys113.
2. **Chemistry** — the bond forms.

Step 2's rate is a property of the **warhead class**, not of the individual
molecule. All chloroacetamides react at broadly similar intrinsic rates; that is
what makes a warhead a class at all. So essentially all molecule-to-molecule
variation lives in step 1.

Rank on step 1.

## What that choice discards, and why each deserves it

| discarded | reason |
|---|---|
| **the docking score** | it estimates affinity, a proxy for step 1 that has been measured as carrying no signal here — five independent levels have failed (D0041, D0046, D0036, D0038/D0044, D0057) |
| **covalent docking** | it assumes step 2 already happened. It can say how an adduct sits once formed; it cannot say whether the molecule would ever react |
| **MM-GBSA, MD residence** | they measure the stability of *whatever* pose the molecule adopts, not whether the warhead is presented |

What survives is **geometry** — the one thing crystal structures directly
measure, so ground truth exists.

## The measurements that force it

3IKD, 82 crystal cases, 2026-08-05:

| | |
|---|---|
| a correct pose is **findable** (best-of-9) | **41.5%** |
| the score **picks** it (top-1) | **18.3%** |
| **choosing at random** among the nine | **19.8%** |

The search works. The score is indistinguishable from a coin flip. **Selection is
the bottleneck**, and the bar any method must clear is **random, not the score** —
those are effectively the same number.

---

## The pipeline

### The decomposition that keeps it falsifiable

Two questions, deliberately separate. **A must be settled before B**, because a
near-attack criterion is computed *on a pose* — pick the wrong pose and its
geometry means nothing.

| | question | ground truth |
|---|---|---|
| **A** | does the rule pick the right *pose*? | the 85-case redock benchmark |
| **B** | does near-attack geometry identify molecules that *react*? | 17 verified Cys113 adducts, against the 6-of-16 warhead classes that collapse to an inert amide |

Note the benchmark's 85 cases are **almost entirely non-covalent ligands** — only
4 covalent PDB entries appear, and for those the benchmark selected the
non-covalent ligand. So it validates A and cannot validate B; B needs its own set
built from the covalent complexes.

### Stage 1 — pose generation

3IKD, **free form** (never the adduct), non-covalent, box centred on the reactive
sub-pocket. Two cheap knobs: more modes than nine (does best-of-20 beat 41.5%?),
and **k replicate seeds** — Vina-GPU already draws a fresh seed per invocation.

### Stage 2 — pose selection *(must beat 19.8%)*

Validated on the 85 cases. Cheapest first:

1. **Replicate consensus** — dock k times, cluster all 9k poses, take the mode
   recurring across the most *independent* runs. The primary candidate: it is new
   information rather than a re-reading of the score, and it is what #10
   identifies as the honest version of consensus.
2. **Contact-profile medoid** — `shared/pose_vector.representative`.
3. **Largest-cluster representative.**
4. **BPMD stability** — expensive; only if the free rules fail.

**This runs first and decides the rest.** If a free rule beats 19.8%, that is an
immediate improvement and a calibrated bar for BPMD. If nothing free beats
random, that is the strongest argument that BPMD is the only remaining option —
worth knowing before committing a week of GPU.

### Stage 3 — near-attack screen *(T_3/T_4 only)*

**Mechanism-specific, never a distance cutoff.** `reactive_atom_smarts` and
`mechanism` are already in the warhead library:

| mechanism | approach required |
|---|---|
| `sn2_displacement`, `sn2_ring_opening` | backside, anti to the leaving group (S···C–LG ≈ 180°) |
| `michael_addition` | perpendicular to the alkene plane |
| `snar_displacement` | perpendicular to the aromatic ring |

A chloroacetamide sitting 3.5 Å from the sulfur with its chlorine pointing at it
is chemically dead and passes any distance-only filter.

### Stage 4 — rank

**Binary gate, then continuous rank.** A molecule either reaches near-attack
geometry reproducibly or it does not. Among those that do, rank by the
free-energy cost of pulling the warhead out of the near-attack window — BPMD
biased along **d(reactive atom, SG)**, not whole-ligand RMSD. A molecule whose
scaffold drifts while the warhead stays locked is a *good* answer, and
whole-ligand RMSD would penalise it.

**Scope:** T_3 + T_4 = 7,178 molecules. T_1 and T_2 have no warhead — report
*not applicable*, never ranked last. Pooling different kinds of object on one
axis is what D0043 is a cautionary tale about.

---

## How this fails

Stated up front so none of them is a surprise:

- **Nothing beats random.** Then poses within one Vina ensemble are genuinely
  indistinguishable by any cheap signal and only BPMD remains, at 200–1,000
  poses per week on 6 GPUs.
- **The near-attack gate does not discriminate.** If nearly everything reaches
  it, there is no funnel and stage 4 has nothing to rank.
- **The validation set for B is 17 molecules.** Enough to support "this clearly
  works" or "this clearly does not"; not enough for a marginal claim.

## Pre-registration

Fixed before any candidate is scored, per D0045: the near-attack window per
mechanism, the acceptance bar on the 17 positives, and the rule combining stage-2
reproducibility with stage-4 stability. Each is otherwise a knob that can be
turned until the ranking looks reasonable.
