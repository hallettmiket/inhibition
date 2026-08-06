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

---

# Measured, 2026-08-05

The pipeline above was built and run. What follows replaces the plan for stages
1–3; stage 4 is still unbuilt.

## Stage 1 changed: bias the search, don't filter afterwards

Reactive docking (D0063) replaced dock-then-filter. It biases *sampling* toward
reaction-competent geometry, which addresses the sampling problem rather than
working around it, and it costs ~2 seconds for 40 independent runs.

**It solves the distance and not the chemistry** (D0064). On a test
chloroacetamide, 20/20 poses put the electrophilic carbon 1.55 Å from Cys113's
sulfur — and every one was chemically dead, S–C–Cl median **97.6°** where SN2
needs ~180°. The modified term is a potential in *distance*, isotropic by
construction, so no parameter value encodes an angle. Two changes follow:

- **`r_eq_12` = 3.2 Å, not the published 1.8 Å.** A NAC is a van der Waals
  *contact* — the reactant state. 1.8 Å is a bond distance, past the transition
  state, where the free molecule cannot be while it still holds its leaving group.
- **Neighbour radii unscaled**, so the leaving group keeps its real bulk.

## Stage 3 is now the load-bearing one

`shared/nac_criterion.py`. **Mechanism matters, and not as a refinement**: SNAr
attacks along the ring normal with the leaving group in-plane, so its S–C–LG sits
near **90°** — the exact value that is dead for SN2. One rule across mechanisms
inverts the verdict for one of them.

**Raw viable fractions are not comparable across mechanisms.** The SN2 window
(≥150°) covers 6.7% of approach directions; the perpendicular window (≤45°, two
faces) covers 29.3% — **4.4× wider by solid angle alone**. `isotropic_null()`
computes that baseline exactly and `enrichment()` divides it out. This is D0020's
"rank within warhead class, not globally" arriving from a new direction.

## Stage 2's bar is met, for one class

3IKD, 200 runs per molecule. Positives are ligands **crystallographically bonded
to Cys113** (verified against `_struct_conn`; issue #12 §A). Negatives are
warhead-matched molecules **measured** inactive in AID 504891, shuffled rather
than taken in file order.

| class | mechanism | positives | negatives | AUC | p |
|---|---|---|---|---|---|
| **chloroacetamide** | SN2 | 9 | 30 | **0.872** | **0.0004** |
| snar_chloroazine | SNAr | 2 | 30 | 0.317 | 0.81 |

Chloroacetamide positives enrich **2.39×** over chance against **0.82×** for
negatives; 8 of 9 beat ≥73% of same-class negatives. **This is the first
measurement in the project that separates actives from inactives** — against five
levels of theory that did not (D0041, D0046, D0036, D0038/D0044, D0057, D0061).

It is also one class of two, and chloroacetamide is where the literature already
converged (#12 §A). The SNAr result is not a failure of the idea so much as a
failure of the *window*: negatives score 57% median against a 29.3% baseline, and
a criterion two thirds of random molecules pass is not a gate.

**AID 504891's 34 actives were rejected as positives.** Read as chemistry rather
than as labels, the 11 warhead-bearing ones are frequent hitters — two
rhodanines, an azlactone, an arylidene barbiturate, a furfurylidene indandione,
an embelin-like dihydroxyquinone, two naphthoquinone sulfonylimines — at 3–75 µM
in a 387,000-compound qHTS, plus a cephalosporin whose warhead match is spurious.
Validating a geometric criterion against compounds that hit everything would
confirm nothing.

## What is still open

- **Tighten the SNAr/Michael window**, then re-measure. It is pre-registered, so
  it is changed once, on stated stereoelectronic grounds, and never tuned against
  the positives.
- **Symmetric warheads are skipped** — a molecule whose reactive SMARTS matches
  twice (fumarate/maleate esters read the alkene from both carbonyls) is dropped
  rather than docked once per reactive centre. This currently removes the whole
  Michael class from validation.
- **Stage 4 is unbuilt.** Nothing yet ranks *among* the molecules that pass.
- **The pocket basis is still 6VAJ's** (D0062), unused by this criterion but
  wrong for anything contact-profile based.
