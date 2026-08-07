# Recap of 2.0.0 “Azurite” — what was built, what it proved, and what carries forward

*Written at the open of the ranking rework. Branch `receptor/3ikd-chemist-prepared`,
58 commits, tagged `v2.0.0`. @tt8804, 2026-08-06.*

The one-line version: **2.0.0 built a geometric ranking for covalent candidates,
validated it against crystal structures, then disproved it against physics — and
the disproof is the release's most valuable output.**

---

## 1. Where 2.0.0 started

@mhallet handed over a pipeline that docked generated candidates into **6VAJ**
and ranked them on binding energy and MM-GBSA. Two things were wrong with it for
this problem, and 2.0.0 opened by establishing both.

**The receptor was wrong.** The chemists had prepared **3IKD**, and it is not a
variant of 6VAJ — the two place the pocket **48.6 Å apart**, with different box
centres and sizes. Every 6VAJ number in the project measures the wrong site
(**D0059**). This is what makes 2.0.0 a major version.

**The endpoint was wrong.** Whole-molecule RMSD to a reference pose answers "did
docking reproduce a crystal structure", which is not the covalent question
(**D0062**). And nothing cheap beat chance at picking the right pose — the bar to
clear is *random*, not the docking score, and no simple rule cleared it
(**D0061**).

---

## 2. The idea 2.0.0 was built on

> Docking is bad at *how tightly* a molecule binds. But a covalent inhibitor only
> needs to **present its warhead to Cys113 at the right geometry**. That is a
> question about shape, not energy — and shape is what docking is least bad at.

Three commitments followed, and they remain the right ones:

**Reactive docking, not dock-then-filter.** Bias the *search* toward near-attack
geometry with derivative atom types and a modified vdW pair potential, rather
than docking normally and filtering afterwards (**D0063**). Critically: the
reactive potential is a **sampler, not a criterion** — it changes what poses you
see, so it cannot also be what scores them (**D0064**). The approach angle is
scored separately, which is why the criterion is honest.

**Mechanism-specific geometry.** A thiolate attacking an sp³ carbon needs a
backside approach ≥150° anti to the leaving group. Attacking an sp² carbon needs
approach ≤30° off the plane normal, in the Bürgi–Dunitz window. These are
different criteria and cannot share a threshold.

**An exact isotropic null.** Each mechanism's viable solid angle is computed in
closed form (SN2 6.70%, perpendicular 8.16%), so `enrichment = viable_fraction ÷
isotropic_null` is comparable across mechanisms. Without it, "40% of poses were
viable" means different things for different warheads.

**It validated.** Against crystallographic Cys113 binders vs warhead-matched
measured inactives: chloroacetamide **AUC 0.908**, Michael **0.734** (**D0065**).
Robust across ten disjoint negative draws; size and docking-score confounds ruled
out.

---

## 3. What went wrong, and what each error taught

Every significant defect in 2.0.0 is the same mistake wearing a different
costume: **a value taken by position, name, label, or inheritance rather than by
identity.** Worth stating as a class, because it will recur.

| defect | the value taken wrongly |
|---|---|
| **D0067** — BDHI scored with sp³ geometry at an sp² carbon | the *mechanism name* (`sn2_ring_opening`) over the actual hybridisation. **374 candidates read 0.00×**; after the fix BDHI became the top two classes of nine |
| pose atoms keyed on the PDBQT name field | the atom's *position in a list* over its chemical identity — every carbon is named `C`, so the wrong carbon got measured (2.2 Å vs the true 1.54 Å) |
| reactive typing guarded by a literal type string | a *type name* over the fact of being retyped. Silently deleted the entire SNAr class — 30 negatives and 2 positives — from the validation |
| GUI drew every pose into 6VAJ | an *inherited default* receptor over the receptor the pose was actually docked into. Poses rendered floating in solvent, 48.6 Å from the pocket |
| adduct re-embedded from SMILES | *connectivity* over stereochemistry — drew the opposite configuration at the warhead centre, worth **1.35–2.02 kcal/mol** |
| PLUMED "installed" | the *absence of an error* over a positive check. GROMACS accepts `-plumed` with no PLUMED present; the kernel is `dlopen`ed at runtime |
| residue labels by offset | a *numeric offset* over the residue identity at that position |

The general fix, applied throughout: **verify by identity, and make the check
fail loudly.** Reactive typing is now guarded by diffing against a plain
preparation rather than looking for a named type. Receptor identity is read back
out of the structure. Residue labels are checked against the residue actually
there.

---

## 4. The two metrics, and how they fell apart

**Enrichment** — the fraction of docking runs reaching a near-attack
conformation, over the isotropic null.

**It does not converge (D0068).** The same molecules fall from 2.91× to 0.96× at
10× search effort; rank correlation across efforts **ρ = −0.117**. The cause is
the *window*, not the search: dividing by *every* run puts every mediocre pose in
the denominator, and more searching only adds more of them. In the limit the
fraction approaches the pocket's background rate for anything.

@tt8804 diagnosed this correctly where the first analysis had not — *"the tool
makes the best poses in order; you are over-diluting it"*. Scoring the **top-N by
energy** instead is rank-stable (**ρ = +0.568**, and +0.688 on the direct test),
because a metric defined on a molecule's own best poses cannot be diluted by
adding runs.

**Consensus** — do a molecule's top-10 poses by energy agree with each other
(**D0070**). Rank-stable where frequency is not. Chosen over largest-cluster
because the distribution is bimodal.

**And then a harder result: D0069.** Plain docking energy on 3IKD separates
covalent binders from measured inactives *better than the geometric criterion
does*. The receptor was **not** the explanation for the earlier weak showing. The
geometric criterion is not obviously buying anything over the score it was
supposed to improve on.

---

## 5. The elevation experiment — the centrepiece

Two metrics ranked the same molecules and disagreed almost completely. @tt8804's
sanity check exposed how completely: **397 molecules have a single binding mode
at ≥0.90 pose agreement, and only 4 of them clear enrichment > 5.70.** The
enrichment cut captures **1%** of the well-aligned molecules and misses five
chloroacetamides — the one class with a validated criterion.

So: which is selecting for something physical?

**Design.** Four groups crossing enrichment × consensus within BDHI, matched on
rotatable-bond count (consensus correlates with rigidity at ρ = −0.259), plus the
validated chloroacetamide class, plus **crystallographic positives as an anchor**.

**Pre-registered.** `docs/elevation_prereg.md`, committed to git before any
simulation ran, fixing the groups, the readouts, and a table mapping each
possible outcome to the conclusion it forces.

**Result (D0071).** 37 molecules, 111 tier-1 runs, 0 failures.

| group | median \|Δd\| (nm) |
|---|---:|
| A · high enrichment, high consensus | 0.277 |
| B · low enrichment, high consensus | 0.198 |
| D · low enrichment, low consensus | 0.204 |
| V · chloroacetamide, high consensus | 0.203 |
| **REF · crystallographic** | **0.102** |

All three pre-registered contrasts: **null** (Holm p 0.39–0.88). Every anchor
contrast: **significant** (δ −0.59 to −0.78, p 0.007–0.050). Tier 2 independently
reproduced the null on the completed BDHI groups.

**Neither metric predicts pose stability. The assay that shows this is
demonstrably working.**

Two things make this a result rather than a disappointment:

1. **The anchor.** Without REF, "A ≈ B ≈ D" is indistinguishable from an assay
   that can't separate anything. With it, the null is a statement about the
   metrics.
2. **Tier 1 was free.** `gromacs_explicit` applies no position restraints during
   NVT/NPT, so 300 ps of unrestrained dynamics already ran before every
   production job. The signal was being computed and discarded.

---

## 6. The worked molecule

`t4_72f5671e89cb` — top of T_4, enrichment 6.86, consensus 1.000, 10/10
near-attack poses at 3.04 Å, QED 0.796. Full literature, med chem, MD, BPMD and
covalent workup.

- Holds its pose **54.45 ns** in 100 ns of explicit water, then dissociates.
- Warhead in the near-attack window for **7.5%** of the bound phase — visited,
  not inhabited.
- **Rank 37 of 37** on tier-1 warhead stability, reproduced across two
  independently built systems.
- **D0072 — NO GO.**

The decision record is explicit that dissociation at 54 ns is *not* the reason:
for a covalent binder the warhead need only reach the window. The reasons are
that its rank came from metrics that predict nothing, it is the cohort's worst on
the readout that works, and its warhead class has **zero crystallographic
positives**.

---

## 7. What carries forward into the rework

### Keep — these earned their place

- **3IKD, and receptor labelling on every number.** Non-negotiable.
- **The mechanism-specific criterion and the isotropic null.** The chemistry is
  right; D0067 was a defect in applying it, not in the idea.
- **Reactive docking as a sampler, scored separately.**
- **Consensus over whole-population frequency.** ρ = +0.568 vs −0.117.
- **Top-N windows over whole-population fractions**, everywhere.
- **The elevation suite** (`docs/elevation_example.md`) — cohort, anchor,
  pre-registration, tiers 1–4.
- **Pre-registration as standard practice.** It cost us the answer we wanted and
  that is precisely when it paid.

### Fix — known-unfixed, inherited

- `mmgbsa.RECEPTOR_PDB` **defaults to 6VAJ**, and every covalent path in the repo
  takes that default. Should be required, not defaulted.
- `nac_rank.refine()` counts `failed:` rows as done when resuming.
- **Covalent MD never ran.** Topology is built and verified; only the trajectory
  is missing.
- **No tier-3 baseline on a crystallographic positive** — so no 100 ns residence
  number in this project is interpretable, including the lead's 54 ns.
- The **accessible-decoy-cysteine control** is still untested: does the criterion
  measure Cys113 *recognition* or merely warhead exposure? Cys57 turned out to be
  buried, so the first attempt was inconclusive.
- Chemist ruling still outstanding on **N-activated acrylamides** — 97% of T_3
  (D0066).

### Open questions the rework has to answer

1. **What replaces enrichment?** D0068 says the top-N fraction is the right
   shape. If it is renamed or redefined, that is a **major** by our own rule
   (`docs/versioning.md`).
2. **Does the geometric criterion beat plain docking energy at all?** D0069 says
   not on 3IKD. If it does not, the honest move is to say so.
3. **What does the shortlist get built on?** The pool is **397 single-mode
   molecules**. Validation attaches to the *mechanism*, not the class name, so
   **309 of 397 (78%) sit in a mechanism that cleared validation** — Michael
   addition 293 (acrylamide 245, naphthoquinone 48) and SN2 displacement 16.
   **But see D0073**: consensus *depletes* validated chemistry relative to the
   library (90.3% → 77.8%, OR 0.34), and the best-validated class has the worst
   pass rate. "Consensus plus the validated class" are pulling against each
   other, and the conflict is larger than it looked.
4. **What would validate BDHI?** It is the largest well-aligned non-acrylamide
   group and it has zero positives. Nothing in 2.x can rank it honestly.

---

## 8. The standing lesson

2.0.0's most reusable output is not a ranking. It is a way of working that caught
its own errors:

- **Anchor every measurement** to something known to be true. The tier-1 null is
  interpretable only because crystal structures went through the same pipeline.
- **Pre-register** when the outcome has more than one plausible reading.
- **Verify by identity, never by position, name, or inherited default.** Every
  significant defect in this release was that one mistake.
- **A null is a result** when the assay is shown to work.
