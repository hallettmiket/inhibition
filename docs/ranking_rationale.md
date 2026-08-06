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

| class | mechanism | positives | negatives | pos enrich | neg enrich | AUC | p |
|---|---|---|---|---|---|---|---|
| **chloroacetamide** | SN2 | 9 | 30 | **2.39×** | 0.82× | **0.822** | **0.0020** |
| **naphthoquinone_c2** | Michael | 4 | 30 | **2.39×** | 1.01× | **0.800** | **0.0288** |
| snar_chloroazine | SNAr | 2 | 30 | 2.66× | 2.39× | 0.558 | 0.41 |
| **pooled on enrichment** | — | 15 | 90 | **2.39×** | 1.14× | **0.722** | **0.0031** |

**Two independent mechanisms separate.** This is the first measurement in the
project that distinguishes actives from inactives at all — against five levels of
theory that did not (D0041, D0046, D0036, D0038/D0044, D0057, D0061). Pooling is
legitimate here *only* because `enrichment()` has divided each mechanism's own
baseline out; the raw fractions are not comparable and pooling them gives
AUC ≈ 0.5, which is an artefact of the windows rather than a result.

**It replicates.** AutoDock-GPU reseeds every run, so each screen is an
independent measurement. Four runs of the chloroacetamide arm gave AUC **0.872,
0.881, 0.852, 0.822** (mean 0.857). Quote the spread, not a single run.

**SNAr does not separate, and cannot be settled here.** Only two SNAr ligands
have ever been crystallised at Cys113, so n = 2 is underpowered by construction —
which is exactly issue #12 §A's finding of 3 verified chemotypes against a
statistical floor of 6. Its negatives also enrich 2.39× over chance on their own,
suggesting the pocket steers most chloroazines into perpendicular approach
regardless of whether they bind.

**AID 504891's 34 actives were rejected as positives.** Read as chemistry rather
than as labels, the 11 warhead-bearing ones are frequent hitters — two
rhodanines, an azlactone, an arylidene barbiturate, a furfurylidene indandione,
an embelin-like dihydroxyquinone, two naphthoquinone sulfonylimines — at 3–75 µM
in a 387,000-compound qHTS, plus a cephalosporin whose warhead match is spurious.
Validating a geometric criterion against compounds that hit everything would
confirm nothing.

### The obvious confound: is this just measuring molecular size?

A small, floppy molecule might reach a near-attack geometry trivially, which
would make the whole result an artefact of the positives being smaller. Tested,
and it is not:

| | |
|---|---|
| enrichment vs heavy-atom count | ρ = **−0.09** (p = 0.36) |
| enrichment vs molecular weight | ρ = −0.07 (p = 0.48) |
| enrichment vs rotatable bonds | ρ = −0.14 (p = 0.17) |
| positives vs negatives: HAC / MW / RotB | p = 0.97 / 0.33 / 0.64 — indistinguishable |
| chloroacetamide, negatives restricted to the positives' size range (n = 21) | **AUC 0.812, p = 0.0040** |

All three correlations are null and *negative* in sign — if anything larger
molecules score slightly lower, the opposite of the failure mode. Size-matching
the negatives leaves the result essentially unchanged (0.812 vs 0.822).

## What is still open

- **Stage 4 is unbuilt.** Nothing yet ranks *among* the molecules that pass the
  gate. This is now the main gap.
- **SNAr is underpowered, not disproven.** n = 2 cannot settle it; it needs
  either more solved SNAr structures or an orthogonal source of positives.
- **The negatives are HTS inactives**, which is weak evidence per molecule. The
  result rests on 15 crystallographic positives against 90 of them; it would be
  stronger against measured *non-covalent* binders, which would test recognition
  rather than warhead presence.
- **The pocket basis is still 6VAJ's** (D0062), unused by this criterion but
  wrong for anything contact-profile based.

*Closed since the first draft of this section: the perpendicular window was
completed with the Burgi-Dunitz constraint (once, on stereoelectronic grounds);
symmetric warheads are now docked once per reactive centre, restoring the whole
Michael class to validation.*
