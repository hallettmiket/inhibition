# What kind of molecule does each T_i actually produce?

*Chemotype characterisation of the four generation approaches, Dance with
Inhibition (Pin1). Written for Mike; describes chemical matter, not merit —
see the caveats in §0 before reading any number as "this approach is
better."*

Data: latest frame per approach — `D1_19` (T_1, n=4,803), `D2_19` (T_2,
n=1,882), `D3_16` (T_3, n=5,396), `D4_29` (T_4, n=1,782). Descriptors read
from existing columns (not recomputed) except Murcko scaffolds, ring-system
macrocycle flags, and Tanimoto similarity, which were computed here with
RDKit 2025.09.5 (ECFP4 = Morgan radius 2, 2048-bit, ignoring stereochemistry)
via `/data/lab_vm/envs/dwi_cheminf/bin/python3.11`.

## 0. Caveats that bear on every number below

- **The shortlists are size-selected, not chemotype-selected.** D0043:
  spearman(heavy atoms, own rank metric) is −0.617 (T_1), −0.479 (T_3),
  −0.230 (T_2), +0.181 (T_4). T_3's shortlist median is 39 heavy atoms
  against 25 in its pool. Every table below reports **POOL** numbers as the
  description of what the *method* produces, and calls out shortlist figures
  only where they add something the pool doesn't show.
- **T_1 is structurally the least sound of the four.** 26.7% of its pool
  (1,284/4,803, `synth_fail=True`) fails a structural rule in
  `shared/synthesizability.py` — dominated by `strained_small_ring_fusion`
  (643 of the 1,284 failures) and `geminal_diol_or_hemiketal` (368), with a
  smaller `peroxide_or_higher` tail (255) — against 0.8% T_2 (15/1,882),
  0.8% T_3 (44/5,396), 0.0% T_4 (0/1,782). This is discussed as a chemotype
  fact in §1, not smoothed over.
- **Covalent molecules have two forms.** `canonical_smiles` is the parent
  (pre-reaction); `adduct_smiles` is the species actually docked, after the
  leaving group departs. Section 2 and the warhead discussion are explicit
  about which is being described. All Tanimoto/scaffold work below uses the
  **parent** form for a fair four-way comparison (T_1/T_2 have no adduct
  form at all).
- **No ranking here is validated** (non-covalent gate AUC 0.599, EF1% 0.0;
  covalent gate underpowered — D0016). Nothing below should be read as "T_i
  found a better/worse molecule." This is a description of chemical space
  occupied, full stop.

---

## 1. Character sketch — one paragraph each

**T_1 (DiffSBDD de novo)** reads as an unconstrained pocket-filling exercise
rather than a chemical series. Median MW 221 (IQR 141–285), median 16 heavy
atoms — smaller than any other approach, and bimodally so: 1,442 pool
members (30%) sit in a genuine fragment band (10–17 heavy atoms) and 1,903
(40%) in a lead-like band (18–45), with 1,458 more generated below/above the
window entirely and retained-but-labelled `too_large`/`degenerate_too_small`
rather than discarded. It is the most saturated set (median fraction sp3
0.75, lowest median cLogP 0.59, virtually no aromatic rings at the median)
and carries the most stereocentres per heavy atom (median 2, up to 32 in one
outlier). It is also the most scattered internally — mean pairwise Tanimoto
0.083, median nearest-neighbour similarity 0.31 — meaning most T_1 molecules
don't structurally resemble even their closest sibling in the same pool. And
it is the least synthetically sound: 26.7% (1,284/4,803) trips a structural
rule, dominated by `strained_small_ring_fusion` (643 hits),
`geminal_diol_or_hemiketal` (368), and `peroxide_or_higher` (255) — fused
strained rings and carbonyl-hydrate/O–O motifs a chemist would never draw as
a target, not the fused-bicyclic congestion rule (`adjacent_quaternary_
ring_carbons`, only 20 hits — that rule was deliberately narrowed after an
earlier version over-fired on ordinary fused bicyclics, per
`shared/synthesizability.py`). 822 molecules (17%)
carry no ring system at all, and Brenk-type reactive/unstable-group flags
hit 68% of the pool (isolated alkenes, non-ring acetals/ketals, aldehydes,
thiols, even a phosphorus-containing tail) — the signature of a generator
with no scaffold prior and no chemist's instinct for what's isolable.

**T_2 (CReM/ATRA neighbourhood)** is, structurally, almost a single
compound family wearing different substituents — which is exactly what
"derivative neighbourhood of a named seed" should look like. Median MW 399
(IQR 373–412, the tightest MW spread of any set), median cLogP 6.28 (by far
the greasiest — driven by ATRA's own polyene/terpenoid character), and a
single Murcko scaffold (`C1=CCCCC1`, ATRA's cyclohexenyl ring) accounts for
576/1,882 molecules (31%); the top 10 scaffolds cover 45%. It is the
tightest internally diverse set by a wide margin — mean pairwise Tanimoto
0.477, median nearest-neighbour similarity 0.84 — and it is almost entirely
achiral (87% have zero stereocentres, vs. 0–19% for every other approach),
because CReM edits ATRA's periphery and rarely touches the one existing
stereocentre. 86.5% of the pool trips a Lipinski rule (essentially all on
cLogP > 5 alone — none trip two rules, so this is a single-axis liability,
not general non-drug-likeness), and 90.5% carry a `Michael_acceptor|polyene`
Brenk flag inherited directly from ATRA's own conjugated tail, not from
anything CReM added. It is non-covalent, achiral, greasy, and structurally
monotonous by construction.

**T_3 (REINVENT4 LibInvent on sulfopin+acrylamide)** is a fixed
covalent core wearing a genuinely diverse wardrobe of attachments — the
opposite failure mode from T_2's near-single-scaffold problem, despite also
starting from one fixed seed. Median MW 380 (IQR 347–421), median cLogP 1.39
(the most polar of the four, TPSA median 98), median QED 0.66 (highest of
any approach at the pool level). 2,354 distinct Murcko scaffolds across
5,396 molecules — LibInvent varies not just the terminal R-group but the
*linker chemistry* connecting it to the fixed sulfolane-acrylamide core
(benzamide, phenylsulfonamide, phenylurea, phenoxyacetamide all appear among
the top five scaffolds), so "decoration" here means new bonds and new
functional classes at the attachment point, not just aryl swaps. It is a
near-obligate carrier of the D0026 `acyclic_imide` flag (58%, 3,128/5,396)
because the scaffold nitrogen already bears the acrylamide carbonyl and
LibInvent frequently adds a second acyl there — a structural consequence of
*where* it decorates, carried as an excused flag rather than filtered. 24.3%
fails the alert gate outright (`rejected_at='alerts'`) on flags beyond the
excused imide (aliphatic long chains, nitro groups, beta-keto/anhydrides).
Every molecule, by construction, is a Michael-acceptor acrylamide.

**T_4 (warhead × R-group combinatorial)** is not one series but nine
parallel series sharing a core and an R-group library, one per warhead
class, and the differences between classes are large enough that "T_4's
character" is really nine characters. Median MW ranges from 348
(acrylamide) to 508 (sulfamate_acetamide) across classes; median cLogP from
0.5 (sulfonate_acetamide) to 3.2 (naphthoquinone_benzo); TPSA from 70
(acrylamide/chloroacetamide) to 125 (sulfamate_acetamide, snar_chloroazine).
It is the most ring-dense set (median 4 total rings, vs. 1–2 for the other
three) because the fixed sulfopin core plus an aryl/heteroaryl R-group
plus, for four of the nine classes, a warhead that is itself ring-based
(naphthoquinone, BDHI-isoxazoline) stacks up rings quickly. 924 distinct
Murcko scaffolds across 1,782 molecules, but many recur in clusters of
exactly 9 — the same R-group decorated with all nine warheads collapses to
the same Murcko scaffold once the (largely acyclic) warhead is stripped,
which is a clean structural signature of "same R-group library, warhead
swapped." It has the highest formal-charge incidence of any approach (9.1%
non-zero, from the nitro-bearing snar_chloroazine and sulfonate/sulfamate
classes) and, at 24.9%, by far the highest PAINS incidence (quinone
substructures in the two naphthoquinone classes are classic PAINS hits).
0.0% synthesizability-rule failures — unsurprising, since every molecule is
assembled from a curated, literature-anchored warhead and R-group library
rather than generated de novo.

---

## 2. Comparison table (POOL, i.e. all generated candidates, not shortlists)

| Axis | T_1 (de novo) | T_2 (ATRA/CReM) | T_3 (LibInvent) | T_4 (combinatorial) |
|---|---|---|---|---|
| n (pool) | 4,803 | 1,882 | 5,396 | 1,782 |
| Mechanism | non-covalent | non-covalent | covalent (acrylamide, fixed) | covalent (9 warhead classes) |
| MW, median (IQR) | 221 (141–285) | 399 (373–412) | 380 (347–421) | 439 (390–466), range by class 348–508 |
| Heavy atoms, median (IQR) | 16 (9–19) | 29 (27–30) | 25 (23–28) | 28 (25–32) |
| cLogP, median (IQR) | 0.59 (−0.52–1.78) | 6.28 (5.48–6.93) | 1.39 (0.74–2.09) | 2.08 (1.22–2.88), range by class 0.5–3.2 |
| TPSA, median | 70 | 50 | 98 | 87, range by class 70–125 |
| frac sp3, median | 0.75 | 0.44 | 0.36 | 0.44 |
| Aromatic rings, median | 0 | 1 | 1 | 2 |
| Total rings, median | 1 | 2 | 2 | 4 |
| Stereocentres, median (frac. 0) | 2 (19%) | 0 (87%) | 1 (0%) | 1 (0%) |
| QED, median | 0.48 | 0.39 | 0.66 | 0.67 |
| Distinct Murcko scaffolds / n | 2,643 / 4,803 (55%) | 641 / 1,882 (34%) | 2,354 / 5,396 (44%) | 924 / 1,782 (52%) |
| Top-scaffold share | 822/4,803 (17%, *acyclic — no ring*) | 576/1,882 (31%) | 572/5,396 (11%) | 9/1,782 (0.5%, max repeat) |
| Macrocycle (≥8-ring) present | 2.6% | 0.2% | 0.2% | 0.0% |
| Lipinski ≥1 violation | 8.7% | 86.5% (cLogP alone) | 6.0% | 7.6% |
| Veber ≥1 violation | 12.1% | 1.4% | 5.3% | 4.4% |
| PAINS > 0 | 1.9% | 0.2% | 1.2% | 24.9% |
| Brenk > 0 | 67.9% | 99.5% | 100% | 46.7% |
| NIH > 0 | 45.0% | 98.6% | 100% | 78.8% |
| Structural synth-rule failure | 26.7% | 0.8% | 0.8% | 0.0% |
| Alert-gate rejection (own pipeline) | n/a (no gate) | n/a | 24.3% | 5.6% |
| Mean pairwise Tanimoto (within) | 0.083 | 0.477 | 0.421 | 0.289 |
| Median NN Tanimoto (within) | 0.31 | 0.84 | 0.76 | 0.85 |
| SAscore, median (not used as criterion) | 4.14 | 3.29 | 3.40 | 3.69 |

Axes that turned out **not** to discriminate: formal charge is 0 for
essentially everyone except T_3 (0.3%) and T_4 (9.1%, driven entirely by the
nitro/sulfonate-bearing warhead classes) — this is a warhead-choice
artefact, not a series-level property worth dwelling on. HBD/HBA medians
cluster in a narrow 0–8 band across all four and don't separate the sets
usefully on their own (TPSA, which folds them together with polarity,
does).

---

## 3. What is genuinely shared across all four

Less than might be expected, and most of what's shared is a design decision
rather than an emergent chemical property:

- **All four target the same pocket geometry**, so all four cluster in a
  MW/heavy-atom range compatible with a mid-size druggable-cleft occupant
  (roughly 200–500 Da at the population level) — no approach produces
  peptide-sized or PROTAC-sized matter.
- **Formal charge is neutral for the overwhelming majority** of every set
  (T_2 and T_1: 100%; T_3: 99.7%; T_4: 90.9%), i.e. none of the approaches
  is systematically generating ionisable/zwitterionic chemotypes.
- **PAINS incidence is low-to-moderate everywhere except T_4's quinone
  classes** — as a population, none of the four is dominated by classic
  frequent-hitter chemotypes.
- **Every approach's SAscore sits in the same narrow synthetic-accessibility
  band** (medians 3.3–4.1) — but per `shared/synthesizability.py`'s own
  rationale ("WHY NOT SA SCORE"), SAscore is a resemblance-to-known-chemistry
  statistic, not a route-existence claim, and it is the one axis Mike asked
  us to treat as a comparator only. It is
  worth noting it does **not** discriminate T_1 (worst by the structural
  rules) from the other three — SAscore's median (4.14) is close to T_4's
  (3.69) despite a 26.7% vs 0.0% gap on the rules that actually catch
  impossible structures. This is direct evidence for why the project uses
  the SMARTS rules instead of SAscore as the synthesizability gate.
- Beyond pocket-scale MW and neutrality, there is **no shared chemotype**:
  no functional group, ring system, or scaffold family appears as a dominant
  motif across even two of the four approaches, other than the covalent pair
  sharing their fixed core by construction (§4).

---

## 4. Where any two approaches are closest, and where most disjoint

Cross-approach max-Tanimoto per molecule (parent/`canonical_smiles`, ECFP4,
computed both directions; "combined median" = median of the per-molecule
max-similarity-to-other-set values, pooled across both directions):

| Pair | Combined median max-Tanimoto | 90th pctile (larger direction) | Max observed |
|---|---|---|---|
| T_3 – T_4 | **0.456** | 0.56 | up to fingerprint-identical (see below) |
| T_2 – T_3 | 0.212 | 0.30 | 0.449 |
| T_1 – T_3 | 0.224 | 0.29 | 0.447 |
| T_1 – T_2 | 0.176 | 0.26 | 0.400 |
| T_1 – T_4 | 0.157 | 0.26 | 0.333 |
| T_2 – T_4 | 0.148 | 0.18 | 0.277 |

**Closest pair, by a wide margin: T_3 and T_4.** This is not a coincidence —
both are built on the identical sulfopin core (`N(...)C1CCS(=O)(=O)C1`) and
T_4's warhead library includes `acrylamide` as one of its nine classes,
which is also T_3's single fixed warhead. Stripped of stereochemistry, **6
molecules are constitutionally identical** between T_3's pool and T_4's
198-member acrylamide subset (e.g. `C=CC(=O)N(Cc1ccc2ccccc2c1)C1CCS(=O)(=O)C1`
appears in both — T_3 generates it with unspecified stereochemistry, T_4's
enumeration carries the config-pinned `[C@@H]` centre, so the two are
fingerprint-identical under a 2D/non-chiral ECFP4 but are technically
distinct `canonical_smiles`, hence **zero exact-string overlap** even here —
consistent with what Mike already knew). That 6/198 (3.0% of T_4's
acrylamide subset) is a real, small, and mechanistically explicable
intersection: an independent generative decorator (LibInvent) and a
hand-curated R-group library occasionally reach for the same simple,
common substituent (here, a 2-naphthylmethyl group) even without sharing an
algorithm.

**Most disjoint pair: T_2 and T_4** (combined median 0.148, 90th percentile
0.18, max observed only 0.277 across 1,882 × 1,782 = 3.35M comparisons).
This tracks the character sketches directly — T_2 is a greasy, achiral,
non-covalent polyene neighbourhood built from ATRA, and T_4 is a
ring-dense, more polar, covalent, formally-charged-in-part combinatorial
set built from sulfopin. They share almost nothing: different mechanism,
different seed, different ring chemistry, different charge profile.
**T_1 and T_2 are nearly as disjoint** (0.176 combined median) despite both
being non-covalent — T_1's unconstrained pocket-filling and T_2's
ATRA-anchored neighbourhood simply don't converge on similar chemical
matter.

The two "seeded" approaches (T_2 from ATRA, T_3/T_4 from sulfopin) do not
resemble each other (T_2–T_3 combined median 0.212, T_2–T_4 0.148) — the
two seeds themselves are structurally unrelated (a retinoid vs. a
sulfolane-based covalent binder), so inheriting a neighbourhood from one
gives no proximity to a neighbourhood inherited from the other.

---

## 5. What surprised us, and what looks like a configuration artefact rather than a chemical fact

1. **T_3 and T_4 nearly-collide by construction, not by convergent
   discovery.** The 6-molecule constitutional overlap (§4) is genuine, but
   it exists *because* both approaches were handed the same fixed core and
   T_4 happens to include acrylamide as one of nine warhead choices — not
   because two independent methods "discovered" the same chemistry from
   different starting points. Read this as evidence the two approaches
   probe the *same local neighbourhood* around one scaffold, one from a
   generative-decoration angle and one from a library-enumeration angle,
   not as evidence of chemical convergence more broadly.
2. **T_2's Brenk/NIH alert rate (99.5%/98.6%) is almost entirely one
   inherited flag, not 1,882 independently-flawed molecules.** 1,703 of
   1,882 (90.5%) carry exactly `Michael_acceptor_1|polyene|polyene` —
   ATRA's own conjugated tail reads as a Michael acceptor to Brenk's SMARTS
   even though it isn't a designed electrophile and T_2 is a non-covalent
   approach. Treating this as "T_2 produces reactive molecules" would be a
   misreading; it is a property of the seed, present before CReM edits
   anything.
3. **SAscore does not track the structural-rule failure rate at all**
   (§3) — T_1's median SAscore (4.14) is unremarkable relative to the other
   three despite failing structural checks at 26.7% vs. ≤0.8% everywhere
   else. This is the clearest direct evidence in the data for the project's
   own stated reason (`shared/synthesizability.py` docstring) for not using
   SAscore as a synthesizability gate.
4. **T_4's per-warhead-class spread is wide enough that "T_4" as a single
   row in a table is misleading** — MW spans 348–508 and cLogP spans 0.5–3.2
   across its nine classes, a range comparable to the spread *between* the
   other three approaches entirely. Any future comparison should probably
   treat T_4 as nine series, consistent with how T_4 itself ranks
   internally (`per_class_quota`, D0020) rather than as one series.
5. **T_1's fragment/lead-like split is a filter artefact worth remembering
   when reading its physchem medians.** 1,458 of 4,803 generated candidates
   (30%) sit outside the reported size window entirely
   (`degenerate_too_small`/`too_large`) and are retained-but-labelled rather
   than removed (see `config/approaches/t1_de_novo.yaml`, "THE FLOOR KILLS
   DEGENERATE OUTPUT ONLY"). All pool statistics in this document include
   them (n=4,803 throughout) — excluding them would shift T_1's median MW
   and heavy-atom count upward and should be done deliberately, not by
   accident, in any downstream analysis.
6. **The synthesizability gap (26.7% vs. 0.0–0.8%) is the single most
   discriminating number in this whole document**, more so than any
   physchem axis — it says T_1's lack of a scaffold prior costs it in
   *makeability*, not just in the diffuse "greater diversity" sense one
   might assume from "de novo, no seed."

---

## Methods note (for reproducibility)

- Descriptors (MW, HAC, cLogP, TPSA, HBD, HBA, rot_bonds, ring counts,
  frac_sp3, formal_charge, n_stereocenters, QED, SAscore, alert flags) are
  read directly from the frame — not recomputed — except where noted.
- Murcko scaffolds: `rdkit.Chem.Scaffolds.MurckoScaffold.GetScaffoldForMol`
  on `canonical_smiles`, canonicalised with `Chem.MolToSmiles`.
- Macrocycle flag: any ring of ≥8 atoms in RDKit's `RingInfo.AtomRings()`.
- Tanimoto: ECFP4 = `AllChem.GetMorganFingerprintAsBitVect(mol, 2,
  nBits=2048)` on `canonical_smiles` (parent form, no stereochemistry
  encoded in the fingerprint bit-hashing by default), `DataStructs.
  BulkTanimotoSimilarity`. Within-approach mean pairwise excludes
  self-pairs; nearest-neighbour = max similarity to any other molecule in
  the same pool. Cross-approach = max similarity of each molecule in set A
  to every molecule in set B, both directions reported.
- Lipinski violation = any of MW>500, cLogP>5, HBD>5, HBA>10. Veber
  violation = rot_bonds>10 or TPSA>140.
- All computation via `/data/lab_vm/envs/dwi_cheminf/bin/python3.11`
  (RDKit 2025.09.5, pandas 3.0.3) against the four latest append_only
  frames listed at the top of this document.
