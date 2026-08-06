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

### Is enrichment just the docking score in new clothes?

No. Measured on 3,602 T_3 candidates carrying both, Spearman ρ against everything
the previous pipeline computed:

| against | ρ |
|---|---|
| Vina `affinity_kcal` | **+0.155** |
| GNINA `cnn_affinity` | **−0.239** |
| GNINA `cnn_score` | −0.020 |
| `size_decorrelated_score` | −0.075 |
| SAscore | +0.055 |

**Max |ρ| = 0.255 — at most ~6% shared variance.** And the *signs* are the
interesting part: on both affinity predictors, higher enrichment goes with
slightly **worse** predicted binding. A molecule that presents its warhead well
is not the same molecule as one that binds tightly, which is precisely what this
framework asserts and what the docking score's five failures (D0041, D0046,
D0036, D0038/D0044, D0061) would predict.

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

### Was one draw of 30 negatives lucky? No — it was conservative

The headline rested on a single draw of 30 warhead-matched inactives from a pool
of 642. `scripts/nac_robustness.py` re-drew **300 per class with a different
seed**, split them into **ten disjoint subsets of 30**, and scored each — ten
independent answers to "what would we have concluded from a different 30", rather
than a bootstrap resampling the one pool we happened to have.

| class | full-pool AUC (300 neg) | ten disjoint draws of 30 | verdict |
|---|---|---|---|
| **chloroacetamide** | **0.908** [0.857, 0.954] | 0.846 – 0.974, median 0.914 | **robust** |
| **naphthoquinone / Michael** | **0.734** [0.593, 0.839] | 0.679 – 0.875, median 0.742 | **robust** |
| snar_chloroazine | 0.451 [0.151, 0.758] | 0.375 – 0.525 | **fragile — no signal** |

**Every disjoint draw of both validated classes clears chance.** The
chloroacetamide worst case (0.846) sits *above* the original point estimate of
0.822, so the first draw was mildly unlucky rather than fortunate.

SNAr straddles 0.5 in every single draw. That settles it as **no signal**, not
merely underpowered — a stronger and more useful statement than n = 2 allowed.

*Caveat: 11 of 915 molecules (1.2%) failed to prepare, all of them negatives —
nine where meeko's reactive-SMARTS match did not reproduce after fragment
stripping. The reasons are unrelated to what they would have scored, but dropping
negatives can only help an AUC, so the asymmetry is recorded rather than
buried.*

### Stage 4, and an uncomfortable result the pre-registration did not anticipate

The five stage-4 rules were fixed and committed **before** the labelled set
existed. Run against it:

| class | incumbent (enrichment) | best pre-registered rule | verdict |
|---|---|---|---|
| chloroacetamide | 0.908 | **C2 enrichment × energy, 0.953** | delta +0.045 [+0.018, +0.076] — **beats it** |
| Michael | 0.734 | C1 0.776 | delta +0.043 [−0.157, +0.182] — indistinguishable |
| SNAr | 0.451 | E1 0.839 | unstable at n = 2; not read as a result |

So **stage 4 exists, for chloroacetamide**: combining how often a molecule
reaches a near-attack conformation with how good that pose is beats frequency
alone, and the improvement clears a pre-registered bar.

**But E1 was mis-specified, and fixing it produces the real finding.** E1 was
labelled "energy alone, expected to fail" and used `best_viable_dg` — which is
conditioned on the geometric gate, so it was never the clean control. The
unconditioned score is `best_dg`:

| class | **`best_dg` alone** | enrichment |
|---|---|---|
| chloroacetamide | **0.915** | 0.908 |
| naphthoquinone / Michael | 0.636 | **0.734** |
| snar_chloroazine | **0.807** | 0.451 |

**The plain docking energy from these runs discriminates about as well as the
geometry** — better for two classes of three. And it is not the geometry in
disguise: `best_dg` correlates with the viable fraction at only ρ = +0.13, −0.01
and −0.12.

**This does not contradict the five prior failures.** None of them measured this:
D0041 was Vina enrichment on 6VAJ, D0046 was pose *recovery*, D0036 was MM-GBSA,
D0038/D0044 was MD residence, D0061 was pose *selection*. This is AutoDock4's
function, on the corrected receptor, asked to separate molecules — a question
none of them put. It is a new measurement, not a reversal.

**What it does do is remove the framework's licence to assume geometry is
necessary.** The honest position is that geometry and energy carry comparable
signal here, that combining them is better than either for chloroacetamide, and
that the outstanding control is **unbiased docking** — these energies come from
runs whose *sampling* was biased toward the warhead–sulfur contact, so a clean
AutoDock run without the reactive potential is required before `best_dg` can be
called an affinity signal at all.

### Confidence ledger — what is established, and what is not

*Asked directly by @tt8804, 2026-08-06: "are we confident in the current rank
build schema?" Kept here rather than in a thread so it stays current. Three
claims sit at different confidence levels and should not be quoted as one.*

| # | claim | status |
|---|---|---|
| **1** | **The measurement.** Enrichment separates crystallographic Cys113 binders from warhead-matched measured inactives. | **Supported, and now robust to the negative draw.** On a 10× larger, independently-seeded negative pool: chloroacetamide **AUC 0.908** [0.857, 0.954], Michael **0.734** [0.593, 0.839]. **Ten disjoint draws of 30 all clear chance** (0.846–0.974 and 0.679–0.875). Replicates across 4 docking seeds. Not molecular size (ρ = −0.09; size-matched AUC 0.812). Not the docking score (max \|ρ\| = 0.255). |
| **2** | **The interpretation.** It measures *Cys113 recognition*, rather than generic warhead exposure. | **UNTESTED.** Everything was measured at one site, so both readings fit. The Cys57 decoy control decides it. If positives enrich equally at both sites, this interpretation is withdrawn — the AUCs would stand but would not mean what §"The one sentence" says. |
| **3** | **The application.** The ranking can order 5,769 candidates. | **Insufficient as built.** 95% CI ≈ 1.12× at 200 runs: the leader is separated, but 1,239 of 1,806 molecules have an interval reaching the top-25 band. Supports "these ~300 deserve a closer look", not "these are the top 25". `--refine-top` addresses it. |

**Structural limits that no amount of compute removes:**

- **15 positives.** Enough for "clearly works" or "clearly does not"; not for a
  marginal claim. Stated before the measurement, and unchanged by it.
- **The negatives are HTS inactives** — weak evidence per compound, since
  single-concentration inactivity has many causes besides failing to bind.
- **One draw of 30 negatives** underpins the headline; the robustness run tests
  exactly that and is pending.
- **SNAr cannot be validated at all** (n = 2). Issue #12 §A's chemotype shortage
  arriving as a hard statistical wall rather than a projection.

**The honest summary.** This is the best-supported ranking the project has
produced and the first thing to separate actives from inactives after five levels
of theory failed — and "best so far" is not "trustworthy enough to commit
synthesis to". The gap between them is claims 2 and 3, plus #12 §D3: measured
actives and inactives are the binding constraint, and no computational control
substitutes for them.

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
