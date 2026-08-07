# Ranking 2.1.0 — anchoring, nulls, and where the freed compute should go

*Design note, not a decision record. @tt8804's framework (#18) with the pieces
worked out. Nothing here is committed until reviewed. 2026-08-06.*

---

## 0. The proposal being designed against

> Consensus as a filter → rank the surviving poses by an anchoring score
> (covalent case additionally uses the warhead's angle/distance) → selection
> (automatic + manual) → elevation suite. Since we are no longer running
> library-scale GROMACS/BPMD, spend the freed compute on more robust docking.

Two additions asked for: **a null for anchoring pairs that are only non-covalent
interactors**, and **a review of pose-generation options**.

---

## 1. What just landed — and a claim withdrawn

**Tier 2 completed while this was being written — 111/111 replicas, 0 failures.**

### 1.1 The core null holds, on both tiers independently

| contrast | tier 1 δ / p(Holm) | tier 2 δ / p(Holm) |
|---|---|---|
| A vs B | −0.469 / 0.391 | −0.156 / 1.000 |
| B vs D | −0.062 / 0.884 | −0.094 / 1.000 |
| A vs D | −0.250 / 0.884 | −0.250 / 1.000 |
| **A vs REF** | −0.781 / **0.0070** | −0.781 / **0.0070** |
| **B vs REF** | −0.750 / **0.0104** | −0.719 / **0.0148** |
| **D vs REF** | −0.594 / **0.0499** | −0.688 / **0.0207** |

Two independent physics readouts, same answer: **neither ranking metric predicts
stability, and the anchor separates from every candidate group on both.** The
tiers agree at Spearman ρ = 0.475, p = 0.003 (n = 37).

### 1.2 V ≈ REF was claimed prematurely, and is withdrawn

An earlier draft of this document, and a summary sent to @tt8804, said the
prereg's fifth reading (**V ≈ REF**) had fired on the strength of tier-2
occupancy — V 0.172 vs REF 0.163, not distinguished. **That was an overreach and
does not survive checking the other readout.**

| readout | V | REF | δ (+ = V more stable) | p |
|---|---:|---:|---:|---:|
| **tier 1 \|Δd\|** — the *pre-registered primary* | 0.203 | 0.102 | **−0.300** | 0.435 |
| tier 2 frac-in-window | 0.172 | 0.163 | **+0.200** | 0.622 |
| tier 1 still-in-window at production start | 0.53 | 0.54 | ~0 | — |

**The two tiers disagree in sign on the same comparison.** On the readout the
pre-registration actually named, V is *worse* than REF, not equivalent to it.

Three separate reasons the claim cannot stand:

1. **Failing to distinguish is not equivalence.** At n = 5 vs n = 8 the test has
   almost no power; a non-significant p is the expected outcome whether or not
   the groups differ. Demonstrating equivalence requires an equivalence test
   against a pre-specified margin, which was never specified.
2. **The prereg forbids it explicitly** — *"no significance claim will be made
   from n = 5"*. That constraint binds the favourable reading exactly as much as
   the unfavourable ones.
3. **The supporting readouts were selected after seeing them.** Tier-2 occupancy
   and still-in-window both favour V; the pre-registered primary does not.
   Choosing the two that agree is the failure the pre-registration exists to
   prevent.

**What can honestly be said:** V is not *distinguished* from REF on two of three
readouts and is worse on the third, at a sample size that cannot resolve any of
them. That is weaker than "behaves like known binders", and it is not a basis for
a synthesis shortlist.

**Consequence for D0073.** The reframing built on this — *"the depletion is a
problem of quantity, not quality"* — is withdrawn with it. D0073 stands as
written: consensus depletes validated chemistry, and whether the survivors are
nonetheless good is **not established**. Testing that properly needs more
chloroacetamides through the suite, not a re-reading of five.

### 1.3 Two protocol facts worth carrying forward

- **`bias_at_exit` separated nothing** (all p ≥ 0.08) and tracks occupancy at
  ρ = 0.974, so the escape-cost term is nearly inert at 3 ns. Occupancy is the
  discriminating tier-2 readout. Do not build a score on the escape cost.
- **Tier 2 is short and unconverged by design** — 108/111 replicas escaped,
  median bias at exit 0.11 kJ/mol. The comparison rests on every molecule
  receiving the identical protocol, not on any single molecule's score being its
  stability.

## 2. Anchoring, generalised

### 2.1 The abstraction

An **anchor** is a triple:

```
anchor := (ligand atom selector, receptor atom selector, interaction type)
```

- **ligand atom selector** — a SMARTS match, as today. Returns atom indices on
  *this* molecule, so indices address the same atoms the conformers hold.
- **receptor atom selector** — residue(s) + atom name(s), resolved **by identity**
  (`shared/receptors.describe`), never by an offset.
- **interaction type** — determines both the geometric predicate and, critically,
  **which null applies**.

Today's covalent criterion is the single instance
`(warhead reactive atom, Cys113 SG, sn2_displacement | perpendicular)`.

### 2.2 The type table

| interaction type | predicate | angular content? |
|---|---|---|
| `sn2_displacement` | d ∈ [2.8, 4.2] Å ∧ S–C–LG ≥ 150° | yes |
| `perpendicular` (Michael / SNAr / BDHI) | d ∈ [2.8, 4.2] ∧ ≤30° off normal ∧ Bürgi–Dunitz 85–125° | yes |
| `hbond` | d(D,A) ∈ [2.6, 3.5] ∧ D–H···A ≥ 120° | yes |
| `salt_bridge` | d(charged, charged) ≤ 4.0 Å | **no** |
| `pi_stack` | centroid d ∈ [3.4, 4.5] ∧ plane angle ≤ 30° or ≥ 150° | yes |
| `hydrophobic` | d(C, C) ≤ 4.5 Å | **no** |
| `proximity` | d ≤ *d*₀ | **no** |

The rows with **no angular content are exactly the ones the current null cannot
serve**, and they are the general non-covalent case.

---

## 3. The null — the part that needs the most care

### 3.1 A defect in the current score, found while designing this

`nac_criterion.viable_fraction` returns the **joint** rate:

```python
viable = in_range and angle >= SN2_ANGLE_MIN        # distance AND angle
```

but `isotropic_null` is **purely orientational** — the solid-angle fraction of
approach directions clearing the window (SN2 6.70%, perpendicular 8.16%).

So today:

```
enrichment = P(distance ∧ angle) / P(angle | isotropic)
```

The numerator carries the distance hit rate; the denominator does not. **The two
have different dimensionality, and the distance term is the one the reactive
potential deliberately biases** (D0064: the reactive potential is a *sampler*,
not a criterion). Enrichment therefore partly measures how well the sampler
worked on a given molecule rather than the molecule's own ability to orient.

This has never been separated, and it is a candidate explanation for D0071 that
does not require the geometric idea itself to be wrong.

**Two consistent repairs:**

**(a) Conditional — recommended for covalent.**
```
score = P(angle viable | distance in window) / P(angle viable | isotropic)
```
Both terms orientational; the sampler's distance bias cancels because it applies
equally to numerator and conditioning set. This is the version D0064's own
reasoning implies and it was never actually implemented.

**(b) Joint with a joint null.**
```
score = P(dist ∧ angle) / [ P(dist | null) × P(angle | null) ]
```
Needs a *positional* null for the distance term — which is exactly what the
non-covalent case needs anyway.

### 3.2 The general principle

> **The null must be produced by the same pipeline as the measurement.**

An analytic null does not know the sampler exists. Where the sampler biases the
quantity being scored, an analytic null cannot cancel that bias, and the score
inherits it. This is the same defect class as everything in `recap_2.0.0.md` §3 —
a value taken from something adjacent to the thing you meant to measure.

### 3.3 Three tiers of null, and every score declares which it used

| tier | null | applies to | cost |
|---|---|---|---|
| **N1 · analytic** | closed-form solid angle (angular) or volume fraction (positional) | separable, purely geometric predicates | free |
| **N2 · empirical decoy panel** | hit rate of the predicate over a fixed decoy set docked through the **identical** pipeline | **any** predicate — the general fallback, and the answer for non-covalent anchors | one-time panel docking per anchor |
| **N3 · within-molecule site control** | same predicate evaluated at matched decoy sites in the same pocket | controls molecule size and flexibility | cheap, per molecule |

**N2 is the mechanism for non-covalent anchors.** It absorbs pocket shape, box
size, sampler bias and ligand-size effects in one measured quantity, and it needs
no analytic derivation for interaction types that have none.

Design constraints on N2, each of which is a way to get it wrong:

1. **The panel must go through the identical sampler**, including any reactive
   bias. A null docked plainly against a measurement docked reactively does not
   cancel anything.
2. **The panel is fixed and versioned**, not resampled per anchor, or the null
   becomes a source of variance in the score.
3. **The panel must not contain molecules that engage the anchor.** For Cys113
   this is knowable; for an arbitrary residue it is an assumption that has to be
   stated.

### 3.4 Why a volumetric N1 is not enough on its own

The naive positional null — sphere volume ÷ box volume — makes the score depend
on the **docking box size**, which is a configuration parameter, not a property
of the molecule or the pocket. A 20 Å box and a 26 Å box would give different
"enrichment" for identical poses. That is precisely the failure mode this project
keeps hitting. If an analytic positional null is used at all it must be
normalised to **pocket-accessible volume**, and even then N2 is the safer answer.

### 3.5 Anchors need a declared dynamic range

A degenerate case that will otherwise waste a campaign: if the anchor residue
lines the pocket and the predicate is loose (`proximity`, 5 Å), then nearly every
pose satisfies it, the null approaches 1, and enrichment approaches 1 for
everything. The anchor discriminates nothing.

So at **anchor-definition time**, compute and store the null rate, and refuse or
flag anchors whose null sits outside a usable band (proposal: reject π₀ > 0.5 or
π₀ < 0.01). An anchor's usefulness is a property that can be checked before any
molecule is ranked, and it should be.

### 3.6 Combining multiple anchors

Once anchors are plural, a molecule has a vector of scores. Do **not** average
them into a single number — that repeats the composite-ranking mistake of letting
a strong component hide a missing one. `shared/composite_rank.py` already has the
right shape: an unmeasured component contributes the interval [0, 1] rather than
an imputed value. Extend that, and keep per-anchor scores visible in the GUI.

---

## 4. Pose generation — what the evidence supports

The reason this matters more in 2.1.0 than 2.0.0: **both the consensus filter and
the anchoring score are computed on the top-N poses.** Pose *ranking* quality is
now the input to everything downstream, not just a display detail.

### 4.1 The number that should shape the design

Best-in-class top-1 pose accuracy in **cross-docking** — the regime we are
actually in, since our molecules are generated and novel — is about **37%**
([GNINA benchmark](https://www.mdpi.com/1420-3049/30/16/3361): Vina 27% → GNINA
37% within 2 Å; redocking 58% → 73%).

**The single top pose is more likely wrong than right.** That is a strong
independent justification for @tt8804's design: filtering on agreement across the
top poses, and scoring over a pose *ensemble*, is the correct response to a
regime where no individual pose can be trusted.

### 4.2 Recommendations

**Adopt — evidence-backed, both already available:**

1. **AutoDock-GPU for sampling + GNINA CNN as pose re-ranker.** The
   [LIT-PCBA benchmark](https://arxiv.org/abs/2605.01681) (15 targets, 578,295
   pairs) found AutoDock-GPU + GNINA rescoring the best single combination
   (median EF1% 2.14). This attacks our actual bottleneck — which pose sits at
   the top — rather than adding sampling we may not need.

   **Carry the D0011 distinction explicitly.** D0011 demoted GNINA's `cnn_*` for
   **covalent affinity**, on gnina's own warning that CNN scoring is uncalibrated
   for covalent docking. That demotion does **not** extend to non-covalent pose
   ranking, which is the CNN's validated strength. If we reuse the D0011 verdict
   by name rather than by scope we will have made a value-taken-by-label error
   about a decision record *about* a value-taken-by-label error.

2. **PoseBusters as a hard validity gate, before consensus.** Already installed in
   `dwi_cheminf`. Checks bond geometry, internal clashes, stereochemistry, ligand
   strain (energy ratio < 100) and protein overlap (< 7.5% vdW volume). The
   [PoseBusters paper](https://pubs.rsc.org/sc/article/15/9/3130/827511/PoseBusters-AI-based-docking-methods-fail-to)
   is explicit that a method can score well on RMSD while producing physically
   impossible poses.

   **It must gate before consensus, not after.** An invalid pose can be
   *reproducibly* invalid — the same clash every run — which inflates consensus.
   Filtering after consensus would let physical nonsense set the agreement score.

3. **Flexible sidechain on the anchor residue.** The anchor distance is measured
   to a specific atom of a specific rotamer. With a rigid receptor we are
   measuring one arbitrary rotamer of Cys113 and calling it geometry. This is
   anchoring-specific and is not a general docking nicety.

**Do not adopt:**

4. **Cross-program consensus rescoring.** The LIT-PCBA study found consensus
   ranking improved *consistency* but **did not beat the single best scorer**. It
   doubles cost for no accuracy. Note this is different from *pose* consensus
   within one program, which is what we use and which remains justified.

5. **DL docking (DiffDock and relatives) as the pose source.** PoseBusters:
   no DL method outperforms classical docking once physical plausibility is
   scored alongside RMSD, and DiffDock underperformed on LIT-PCBA.

**Evaluate separately, not as a pose source:**

6. **Co-folding (Boltz-2 / Chai-1 / AF3).** Prospectively these reproduce >50% of
   novel-ligand poses under 2 Å, and Boltz-2 is reported to discriminate true from
   false positives among docking hits better than any scoring function tested.
   **But its affinity is derived from features largely independent of the final
   ligand pose**, so it cannot feed an anchoring score — anchoring is a statement
   about the pose. Candidate as an orthogonal *molecule-level* filter at
   selection, never as the geometry source. It also does not model the covalent
   bond.

### 4.3 Where the freed compute goes, in priority order

1. GNINA CNN re-ranking of the top-N — directly improves the input to consensus.
2. PoseBusters validity — cheap, and removes a class of false agreement.
3. More runs per molecule — improves top-N stability. **Test, do not assume:**
   D0068 found more search produced lower-energy, *less* reaction-competent poses.
   A properly conditioned score (§3.1a) should be immune to that, which is itself
   a test of the repair.
4. Flexible anchor sidechain.
5. MM-GBSA rescoring of survivors — now affordable, but keep it out of the
   ranking until it has passed §5.

---

## 5. How this gets validated before it ranks anything

**2.0.0's most valuable asset is not a ranking — it is a working stability assay
and a cohort already run through it.** Tier 1 separates crystallographic binders
from generated candidates at p = 0.007, and tier 2 reproduces it on occupancy.

So the new score does not need a new experiment. It needs the *existing* one:

1. Recompute the reworked anchoring score on the **elevation cohort's existing
   poses** — 37 molecules, already docked, already measured.
2. Test it against tier-1 |Δd| and tier-2 occupancy, the same contrasts D0071
   pre-registered.
3. **A score that predicts stability where enrichment and consensus did not is a
   real advance. A score that does not is another null**, and cheaper than the
   first one because the cohort exists.

Pre-register it, for the same reason it was pre-registered last time: the
outcome has more than one plausible reading and the least convenient one has to
stay readable.

**One caution.** The elevation cohort was *selected on* enrichment and consensus.
A new score correlated with either inherits that selection. The cleanest test
adds a fresh stratum selected on the new score — otherwise the cohort can only
falsify, not confirm.

---

## 6. Open questions

- **Is the conditional repair (§3.1a) enough to rescue enrichment?** Cheap to
  test: per-pose distances and angles are not currently persisted, so this needs
  one screen re-run with per-pose output. That output should be persisted anyway.
- **Which non-covalent anchors are worth defining for Pin1?** The basic cluster
  (Lys63, Arg68, Arg69) binds the phosphate of the pSer/pThr-Pro substrate and is
  the obvious candidate. Requires a decoy panel per anchor (§3.3).
- **Does the consensus bar want to be per-mechanism at a fixed pass rate**
  rather than a fixed 0.90? D0073's rigidity skew argues yes; tier 2's V ≈ REF
  argues the bar is doing something right where the chemistry is validated.
- **The chloroacetamide supply problem.** Five molecules clear the bar and they
  behave like crystal structures. That argues for generating a more rigid
  chloroacetamide series rather than lowering the bar — a stage-1 request, not a
  ranking change.

---

## Sources

- [PoseBusters: AI-based docking methods fail to generate physically valid poses or generalise to novel sequences](https://pubs.rsc.org/sc/article/15/9/3130/827511/PoseBusters-AI-based-docking-methods-fail-to) — Chemical Science
- [Benchmarking Single-Pose Docking, Consensus Rescoring, and Supervised ML on the LIT-PCBA Library](https://arxiv.org/abs/2605.01681)
- [High-Throughput, High-Quality: Benchmarking GNINA and AutoDock Vina for Precision Virtual Screening](https://www.mdpi.com/1420-3049/30/16/3361)
- [Boltz-2: Towards Accurate and Efficient Binding Affinity Prediction](https://pubmed.ncbi.nlm.nih.gov/40667369/)
- [Large scale prospective evaluation of co-folding across 557 Mac1-ligand complexes](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12776374/)
