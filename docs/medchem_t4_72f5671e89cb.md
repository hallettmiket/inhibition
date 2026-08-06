# Computational med-chem and ADMET workup — candidate `t4_72f5671e89cb`

*Written 2026-08-06 by Blacksmith (murmurent) for the T_4 team, on request.
Scope: everything computable about this one molecule — ADMET, physicochemistry,
developability, retrosynthesis, stereochemistry, and a like-for-like comparison
with sulfopin. The literature side is
[`docs/lit_t4_72f5671e89cb.md`](lit_t4_72f5671e89cb.md) (Bookworm, same day);
this report does not repeat it and cross-references it where the two meet. A
third report on the same molecule, `docs/simulations_t4_72f5671e89cb.md`,
covers the docking and simulation side. All three were written independently;
where they overlap they agree, and any disagreement should be treated as a
finding rather than reconciled by hand.*

---

## Read this before anything else

**This molecule's rank is not a reason to make it, and one measurement made
today for this report says so directly.**

Three separate things qualify every number below. The first two were given to me
with the task; the third I measured.

**1. The rank comes from a criterion corrected today.** Decision
[D0067](../decisions/D0067-bdhi-was-scored-with-sp3-geometry-at-an-sp2-carbon.md):
`shared/nac_criterion.py` had mapped the mechanism label `sn2_ring_opening` to
sp3 backside-attack geometry, on the strength of the name. BDHI's attacked
carbon is the **sp2** carbon of a C=N, and a thiolate adds perpendicular to that
plane rather than from behind. Until this morning all 374 BDHI candidates read
**0.00× enrichment** — completely unreactive. The fix is chemically correct
independent of any data (and the literature agrees: Byun 2023's DFT puts the
thiolate addition at the sp2 C3, see the Bookworm report §2). But it means this
molecule went from "last in the ranking" to "first in the ranking" in one code
change, a few hours ago, with no intervening experiment.

**2. BDHI is the one warhead class in the ranking with zero crystallographic
Cys113 positives.** The reference set's only BDHI row,
`Byun-2023-BDHI-fragment` in
[`pin1_covalent_cys113_anchors_2.csv`](../data/reference/pin1_covalent_cys113_anchors_2.csv),
carries `smiles_status = UNVERIFIED` and `potency_kinetics = UNVERIFIED`, and
`shared/reference_set.covalent_anchors()` **refuses to return it by default**.
So the class cannot be validated at all in this project's own framework — unlike
chloroacetamide (AUC 0.908 [0.857, 0.954]) and Michael (0.734 [0.593, 0.839]).
Its warhead library row says `structure_status: DESIGNED_UNTESTED` and its
adduct is `INFERRED`, not observed.

**3. NEW, measured for this report: the enrichment does not survive the
convergence control.**
[D0068](../decisions/D0068-enrichment-depends-on-search-effort-and-energy-beats-it-at-convergence.md)
established that near-attack
enrichment is a function of search effort — the top 300 fell from 2.91× to 0.96×
at 2,000 runs, and the crystallographic positives fell identically, so it is not
winner's curse. **No BDHI molecule had ever been run at 2,000 runs**; the whole
2,000-run set is 15 positives and 60 negatives, all chloroacetamide,
naphthoquinone or SNAr. I ran it on this molecule and its two regiochemical /
warhead siblings today:

| molecule | class | 200 runs | **2,000 runs** | median S···C | median off-normal |
|---|---|---|---|---|---|
| **`t4_72f5671e89cb`** (this one) | bdhi_c5 | 6.74× | **0.60×** | 3.15 → 5.17 Å | 24.5° → 61.9° |
| `t4_9c44a3f8892a` (same R-group, C4) | bdhi_c4 | 4.72× | **0.94×** | 3.89 → 5.52 Å | 30.9° → 52.9° |
| `t4_45901f30d2a1` (same R-group, sulfopin's warhead) | chloroacetamide | 1.12× | **0.94×** | 3.45 → 3.58 Å | 121.7° → 105.3° |

My independent 200-run re-measurement (6.74×) reproduces the production value
(7.23×) inside the seed-to-seed spread, so the headline number is real and
reproducible — **it just does not survive more searching.** At 2,000 runs this
molecule is *below chance* (0.60×), and the entire 7.23-vs-1.19 gap that put it
at rank 1 of 5,765 collapses to nothing. More search finds lower-energy poses
that are less reaction-competent, exactly as D0068 described, and BDHI is not an
exception to it.

**The docking energy, by contrast, does converge — and it does not favour this
molecule much.** In the same runs, `best_dg` barely moved with 10× the search:
candidate **−8.05 → −7.94** kcal/mol, `bdhi_c4` sibling −7.68 → −7.65,
chloroacetamide sibling −7.22 → −7.21. So the quantity the ranking is built on
is the unstable one and the quantity the ranking discards is the stable one,
which is
[D0069](../decisions/D0069-plain-docking-on-3ikd-outperforms-the-geometric-criterion.md)'s
finding arriving on this molecule. On the converged energy the candidate leads
its chloroacetamide sibling by 0.7 kcal/mol — a real but ordinary margin, and
one that 20 heavy atoms against 17 partly explains.

`rank_validated = False` on this candidate's frame row, and it should be read
literally. **The honest description of this molecule is "an ordering the
pipeline produced", not "a predicted binder".**

One further piece of context a chemist should have: the corrected criterion did
not just promote this molecule, it promoted its whole class. BDHI is 374 of
5,765 scored candidates (6.5%) and now holds **73 of the top 100** (57 `bdhi_c5`
+ 16 `bdhi_c4`); the validated chloroacetamide class holds **one**. A fix that
moves one unvalidated class from last to first across the board deserves the
convergence control before it is acted on, and the control is above.

---

## 0. What was run, and with what

Every number in this report carries its tool. Where a tool was unavailable or
out of its applicability domain, that is stated rather than papered over.

| purpose | tool | version | env |
|---|---|---|---|
| ADMET prediction | **ADMET-AI** | **2.0.1** (chemprop 2.3.0, torch 2.13.0, RDKit 2026.03.4) | `/data/lab_vm/envs/dwi_admet` |
| descriptors, QED, SAscore, alerts | RDKit + `shared/descriptors.py`, `shared/alerts.py` | RDKit **2025.09.5** | `/data/lab_vm/envs/dwi_cheminf` |
| structural synthesizability rules | `shared/synthesizability.py` | repo | `dwi_cheminf` |
| retrosynthesis | **AiZynthFinder** | **4.4.1** (USPTO expansion + ringbreaker + filter ONNX models; ZINC stock, 17,422,831 compounds) | `~/.local/pin1_tools/aizynth_venv` |
| near-attack re-scoring | AutoDock-GPU, reactive 3IKD receptor (D0063 potential) | repo `scripts/nac_rank.py` | `~/.micromamba/envs/dwi_reactive`, GPU 7 |

**Two environment problems, recorded because they change what you can
reproduce.**

*ADMET-AI's console script `admet_predict` is broken in `dwi_admet`.* chemprop
2.3.0 does an unconditional `import cuik_molmaker` at
`chemprop/featurizers/molgraph/molecule.py:4`, and the installed
`cuik_molmaker_pin` wheel links against the auditwheel-mangled RDKit libraries
of the *pip* rdkit wheel (`libRDKitAbbreviations-a55b7b38.so.1`), which this
env does not have — its RDKit is the conda build. The import fails before any
model loads. [`scripts/medchem_admet.py`](../scripts/medchem_admet.py) works
around it by pre-inserting an empty stub module, which is sound only because
`cuik_molmaker` is referenced exclusively inside the alternative
`CuikmolmakerMolGraphFeaturizer` class and ADMET-AI's checkpoints carry
`SimpleMoleculeMolGraphFeaturizer`. The script does **not** assert that in a
comment: `--assert-featurizer` (on by default) walks every loaded model and
aborts if any featurizer is the stubbed one. It passed. **The env itself was
not modified** (D0010, and it belongs to another user).

*`nac_rank.py --refine-top` cannot re-measure a molecule that once failed.*
`refine()` builds its resume set from every ident present in a chunk file
**regardless of status**, so a row written as `failed: <reason>` counts as
done. My first convergence attempt ran from an env without `gemmi` and wrote
five `failed:` rows; the immediate re-run in the correct env reported
`5 assigned, 0 to do` and scored nothing. Append-only means those rows cannot
be removed, so those five molecules are permanently un-refinable through that
path. I worked around it with a standalone script
([`scripts/medchem_nac_converge_one.py`](../scripts/medchem_nac_converge_one.py))
rather than editing shared code mid-run. **This is a live defect and it is the
project's signature shape — a value (here, "done") taken by position in a file
rather than by identity of what happened.** It should be fixed in
`nac_rank._ids_in` by filtering on `status == "ok"`, with a test.

---

## 1. ADMET

ADMET-AI 2.0.1 ran successfully: 41 TDC endpoints plus 8 physicochemical, each
with a percentile against 2,845 DrugBank-approved drugs. Full output at
`00_outputs/blacksmith/medchem_t4_72f5671e89cb/admet_ai_1.csv`.

### How to read these numbers

**Read the percentiles, not the absolute regression values.** ADMET-AI's
regression heads are unconstrained, and on its *own* DrugBank reference set they
produce physically impossible values at a high rate — measured from
`admet_ai/resources/data/drugbank_approved.csv`:

| endpoint | negative in … of 2,845 approved drugs |
|---|---|
| Volume of distribution (L/kg) | **35.0%** |
| Half-life (hr) | **22.2%** |
| Clearance, microsome (µL/min/mg) | 24.4% |
| Clearance, hepatocyte (µL/min/10⁶ cells) | 15.1% |
| Plasma protein binding (%) | 2.6% |

A predicted VDss of −1.9 L/kg for this candidate is therefore not a number about
this candidate; it is the model saying "low", in a scale it does not respect.
The percentile is the defensible readout. Three of the five molecules I scored
got a **negative half-life**, including sulfopin's own chloroacetamide analogue.

**Two applicability-domain caveats that are specific to this molecule and are
not small.** (i) ADMET-AI is trained on TDC datasets that are overwhelmingly
non-covalent, drug-like molecules; a latent electrophile bearing a C(sp2)–Br is
outside that distribution, and the models have no mechanism for "reacts with the
target". (ii) The molecule that reaches tissue is not the molecule that was
scored: after reaction the bromide is gone. I therefore ran the **adduct form**
as well, and where the two disagree the disagreement is the interesting part.

### Absorption and permeability

| endpoint | candidate | pct | sulfopin | pct | reading |
|---|---|---|---|---|---|
| HIA (human intestinal absorption) | 0.999 | 57 | 0.942 | 30 | absorbed |
| Oral bioavailability (Ma) | 0.838 | 60 | 0.731 | 38 | favourable |
| PAMPA permeability | 0.845 | 57 | 0.751 | 48 | permeable |
| Caco-2 (log 10⁻⁶ cm/s) | −5.06 | 42 | −4.33 | 91 | **candidate ~5× lower than sulfopin** |
| P-gp inhibition | 0.012 | 35 | 0.056 | 47 | not a P-gp inhibitor |
| Aqueous solubility (log mol/L) | −1.84 | 73 | −1.87 | 73 | good, ~1–2 mM |
| BBB penetration | 0.941 | 76 | 0.776 | 55 | **CNS-penetrant predicted** |

Passive absorption looks fine and is consistent with the physicochemistry (MW
364, cLogP 1.12, TPSA 85, HBD 0, 4 rotatable bonds). The one real gap against
sulfopin is Caco-2: −5.06 vs −4.33, i.e. about a five-fold lower predicted
effective permeability, driven by the extra 31 Å² of TPSA and four extra HBAs.

The predicted BBB penetration deserves a note rather than a tick. The
class's clinically-tested relative, acivicin, failed on **dose-limiting CNS
toxicity** (Bookworm §2). A CNS-penetrant covalent electrophile from that
lineage is a combination worth raising before synthesis, not after.

### hERG, hepatotoxicity, mutagenicity, and the rest of the tox panel

| endpoint | candidate | pct | adduct | sulfopin | pct |
|---|---|---|---|---|---|
| **hERG blocking** | **0.134** | 36 | 0.143 | 0.174 | 40 |
| **DILI** | **0.761** | 70 | 0.687 | **0.319** | 45 |
| **AMES mutagenicity** | **0.945** | **98** | 0.791 | 0.729 | 92 |
| Carcinogenicity (Lagunin) | 0.502 | 94 | 0.389 | 0.404 | 88 |
| Skin reaction | 0.840 | 89 | 0.666 | 0.837 | 88 |
| ClinTox | 0.094 | 49 | 0.123 | 0.109 | 53 |
| Acute tox LD50 (log 1/(mol/kg)) | 2.70 | 66 | 2.28 | 2.44 | 50 |

**hERG is clean** — 0.134, below sulfopin, 36th percentile. That is the one
unambiguously good tox number here and it is consistent with the
physicochemistry (no basic-amine-plus-lipophile motif; cLogP 1.12).

**AMES is the standout liability: 0.945, the 98th percentile of approved
drugs.** It falls to 0.791 for the adduct, so some but not all of the signal is
carried by the C–Br. This is what an Ames model *should* do with an
alkyl/vinyl halide, and it does not distinguish "designed covalent warhead" from
"genotoxic alkylator". Take it as a flag to run a real Ames rather than as a
prediction: every covalent-warhead compound scores high here, sulfopin included
at 0.729/92nd percentile. But note the candidate is *higher than sulfopin*, and
the difference survives into the adduct (0.791 vs 0.729), so it is not only the
halide.

**DILI 0.761 vs sulfopin's 0.319 is the largest gap in the whole panel** — a
2.4-fold higher predicted hepatotoxicity risk than the reference compound, and
it stays at 0.687 for the adduct, so it is not the warhead. This is a property
of the scaffold as changed: the isoxazolylmethyl group and the dihydroisoxazole
ring replacing sulfopin's neopentyl and chloroacetamide.

The nuclear-receptor and stress-response panels are quiet: the highest of the
twelve is SR-ARE at 0.255, every other below 0.11.

### CYP inhibition and metabolism

| endpoint | candidate | pct | sulfopin |
|---|---|---|---|
| CYP1A2 inhibition | 0.040 | 44 | 0.006 |
| CYP2C9 inhibition | 0.030 | 48 | 0.052 |
| CYP2C19 inhibition | 0.141 | 58 | 0.173 |
| CYP2D6 inhibition | 0.014 | 29 | 0.043 |
| CYP3A4 inhibition | 0.048 | 49 | 0.041 |
| CYP3A4 substrate | 0.645 | 73 | 0.574 |
| CYP2C9 substrate | 0.216 | 73 | 0.167 |
| CYP2D6 substrate | 0.091 | 45 | 0.089 |

**No CYP inhibition liability**: all five inhibition probabilities are below
0.15, all in the middle of the approved-drug distribution. Predicted to be a
CYP3A4 substrate (0.645), which is unremarkable.

### Clearance and exposure

| endpoint | candidate | pct | adduct | sulfopin | pct |
|---|---|---|---|---|---|
| Clearance, hepatocyte (µL/min/10⁶ cells) | 15.9 | 32 | 19.8 | 60.8 | 70 |
| Clearance, microsome (µL/min/mg) | 7.1 | 31 | **−3.7** | 37.0 | 60 |
| Half-life (hr) | 1.63 | 25 | **−1.55** | 1.91 | 25 |
| Plasma protein binding (%) | 44.6 | 24 | 29.8 | 44.3 | 24 |
| VDss (L/kg) | **−1.90** | 18 | **−2.98** | −2.15 | 16 |

Predicted **low clearance** (32nd/31st percentile, roughly a third of sulfopin's)
and **low plasma protein binding** (~45%, 24th percentile — high free fraction).
Half-life and VDss both land in the bottom quartile, but two of the four adduct
values are negative, so read the direction and not the magnitude. For a covalent
inhibitor, systemic half-life matters less than target-residence anyway; what
matters is exposure long enough to react, and low clearance plus high free
fraction is the right direction for that.

### Structural alerts (RDKit FilterCatalog via `shared/alerts.py`)

Computed twice, per the module's two-tier rule — whole molecule is advisory
(it *will* flag a warhead, that is what a warhead is), decoration is the gate.

| | PAINS | BRENK | NIH | names |
|---|---|---|---|---|
| **candidate**, whole molecule | 0 | 0 | 1 | `halo_imino` |
| **candidate**, R-group only | 0 | 0 | 0 | — |
| **adduct** (post-reaction) | 0 | 0 | 0 | — |
| **sulfopin**, whole molecule | 0 | **1** | **2** | `alkyl_halide`, `alpha_halo_carbonyl`, `primary_halide_sulfate` |

The candidate is *cleaner* than sulfopin on published alerts — but only because
`halo_imino` is the single motif the catalogs carry for BDHI, where
chloroacetamide is covered three times over. **Do not read this as the candidate
being the safer electrophile.** The decoration passes the hard gate with zero
alerts in all three catalogs, which is the part that means something.

**The alert catalogs miss the liability I would actually raise.** See §2.

---

## 2. Physicochemistry and developability

Descriptors from `shared/descriptors.py` — the project's single source, so these
are the same numbers as the T_4 frame carries, by construction.

| | candidate | adduct | sulfopin |
|---|---|---|---|
| formula | C11H14BrN3O4S | C11H15N3O4S | C11H20ClNO3S |
| MW | 364.22 | 285.33 | 281.81 |
| heavy atoms | 20 | 19 | 17 |
| cLogP | 1.12 | 0.40 | 1.29 |
| TPSA (Å²) | 85.0 | 85.0 | 54.45 |
| HBD / HBA | 0 / 7 | 0 / 7 | 0 / 3 |
| rotatable bonds | 4 | 4 | 3 |
| rings (arom / aliph) | 1 / 2 | 1 / 2 | 0 / 1 |
| fraction sp3 | 0.64 | 0.64 | 0.91 |
| **QED** | **0.796** | 0.796 | 0.733 |
| **SAscore** | **4.50** | — | **3.19** |
| stereocentres (assigned/total) | **1 / 2** | 1 / 2 | 0 / 1 (seed) · 1 / 1 (PDB QT7) |
| novelty vs external Pin1 set | 0.737 | — | 0 (is the set) |

### Rule sets — all three pass, cleanly

| rule set | verdict | components |
|---|---|---|
| **Lipinski** | **pass, 0 violations** | MW 364 ≤ 500 ✓ · cLogP 1.12 ≤ 5 ✓ · HBD 0 ≤ 5 ✓ · HBA 7 ≤ 10 ✓ |
| **Veber** | **pass** | rot bonds 4 ≤ 10 ✓ · TPSA 85 ≤ 140 ✓ |
| **Ghose** | **pass, 0 violations** | MW 364 ∈ [160, 480] ✓ · cLogP 1.12 ∈ [−0.4, 5.6] ✓ · MR 74.96 ∈ [40, 130] ✓ · 34 atoms ∈ [20, 70] ✓ |

QED 0.796 sits at the 88th percentile of approved drugs (ADMET-AI's own
DrugBank comparison), above sulfopin's 0.733/80th. On the standard
drug-likeness axes this is a well-behaved small molecule and better-behaved than
the reference compound.

### PAINS and Brenk

Zero PAINS and zero Brenk on the whole molecule; zero on the decoration. See the
table in §1. `shared/synthesizability.py`'s structural-impossibility rules also
return **no violations** — which, per that module's own docstring, means only
that none of the named impossibilities is present, not that the molecule is easy
to make. §3 is where that question actually gets answered.

### The liability nothing in the pipeline flags

**The bond joining the sulfopin core to the warhead is an N,O-acetal.**

In `...N(Cc2cnoc2)C2CC(Br)=NO2`, the ring carbon bonded to the core nitrogen is
**C5 of the dihydroisoxazole, which also carries the ring oxygen O1**. Nitrogen
and an ether oxygen on the same sp3 carbon — and that oxygen is an oxime-ether
oxygen (O–N=C), so the departing group on ring-opening is an oximate, a
better leaving group than a plain alkoxide. This is a classic hydrolytic
liability: acid-catalysed opening to an iminium/oxocarbenium, giving the
secondary amine and a ring-opened oxime. **PAINS, BRENK and NIH do not carry
this motif; SAscore does not see it; the synthesizability rules do not test for
it.** I added it as an explicit SMARTS
(`[NX3;!$(N[CX3]=[OX1,SX1,NX2])][CX4;H1,H2]([OX2])`) and it fires once on the
candidate and once on the adduct — the adduct too, so **reaction with Cys113
does not remove it**.

It is not a property of this molecule. It is a property of the **class**,
created by the attachment choice itself. Scanned across all 1,782 T_4 molecules:

| warhead class | molecules with an N,O-acetal |
|---|---|
| **bdhi_c5** | **198 / 198** (mean 1.03 per molecule) |
| every other class | 6 / 198 (3.0%), all from R-groups that carry one themselves |

`bdhi_c4` — the same warhead attached at C4 instead of C5 — has **none**,
because C4 does not bear the ring oxygen. The library's own note on `bdhi_c5`
says the C5 attachment "follows the usual literature elaboration of BDHI
fragments"; that convention is about elaborating C5 with **carbon**. Putting a
*nitrogen* there is a different compound class, and the label carried over while
the identity did not. That is the shape of defect
[`how_this_project_breaks.md`](how_this_project_breaks.md) catalogues.

And the project's own decision of record already prefers the other arm:
**[D0024](../decisions/D0024-bdhi-c4-naphthoquinone-benzo-on-clean-poses.md)
says "carry `bdhi_c4`"**, reversing D0021, on adduct-form poses (median −3.79
vs −2.87 kcal/mol), verdict UNDERPOWERED. The molecule now sitting at rank 1 is
in the arm that decision deselected.

Two further observations, smaller:

- `charge_ph74 = 1`, `charge_class = cation` in the frame. That comes from
  `obabel -p 7.4`, which protonates aliphatic tertiary amines by rule. This
  nitrogen is an **α-alkoxy amine** flanked by two electron-poor rings; its pKa
  is very likely well below 7.4 and the molecule most likely neutral at
  physiological pH. `shared/ionisation.py` is explicit that it uses obabel to
  match what docking did, not because it is a pKa model — so this is a known
  limitation surfacing, not a new bug. It matters for permeability and for the
  acetal's stability.
- The reactivity window says `in_window: True` (LUMO −7.63 eV, window
  [−7.82, −6.21]). That window is anchored on **six chloroacetamide and
  sulfamate/sulfonate-acetamide actives and nothing else**; its own recorded
  caveat is "chloroacetamide-centric … do not over-trust a single global
  cutoff." A bromo-oxime-ether electrophile landing inside a window built from
  α-halo carbonyls is not corroboration.

### Ligand efficiency

Two energies exist for this molecule and they are not the same measurement.
Both are reported; neither is an affinity.

| metric | AutoDock4 reactive, free ligand | GNINA/Vina, adduct form |
|---|---|---|
| ΔG (kcal/mol) | **−8.06** | −5.03 |
| **LE** (kcal/mol per heavy atom) | **0.403** | 0.251 |
| implied pKd *if this were an affinity* | 5.91 | 3.68 |
| **LLE surrogate** (pKd − cLogP) | **4.79** | 2.57 |

- The −8.06 comes from AutoDock-GPU with the **D0063 reactive potential**
  (`r_eq_12` = 3.2 Å) biasing sampling toward the warhead–sulfur contact. It is
  not a plain dock and the ranking doc is explicit that an unbiased AutoDock run
  is the outstanding control before `best_dg` can be called an affinity signal.
- The −5.03 is Vina's score on the adduct (D0022), taken as the minimum over
  nine modes, in a stratum whose gate verdict is **UNDERPOWERED**
  (AUC 0.542 [0.350, 0.756]).
- There is **no measured activity**, so the LLE figures are surrogates built by
  converting a docking score to a pKd. They are included because a chemist will
  ask; they should not be quoted outside this table.

For scale: sulfopin's ensemble MM-GBSA is −7.58 ± 0.28 kcal/mol over 17 heavy
atoms, LE 0.446 — and in that same measurement 50 of 80 decoys scored better
than it (D0036), which is the clearest available statement of what these
energies are worth on this target.

### Where it sits in known Pin1 chemical space

ECFP4 Tanimoto against the 21 parseable binders in the frozen external
reference set (`master_sha256` `4c8f7d1b…`):

| nearest | Tanimoto |
|---|---|
| Tian chloropyrimidine covalent 6a | 0.263 |
| Reddi-2023-4g | 0.256 |
| **Sulfopin** | **0.250** |
| Reddi-2023-4d | 0.244 |
| Potter/Astex indole-furancarboxamide | 0.102 |

Nothing above 0.27. The project's own `novelty_external` is 0.737. This is
**outside** known Pin1 chemical space — which is the point of T_4, and also
means there is no analogue whose behaviour can be borrowed. The Bookworm report
independently found no PubChem/ChEMBL match for the molecule and no record of
this scaffold ever carrying a BDHI warhead.

---

## 3. Synthetic accessibility, properly

SAscore 4.50 is "moderate" and is a resemblance metric, not a claim that a route
exists. AiZynthFinder 4.4.1 was run against the USPTO expansion + ringbreaker
policies with the USPTO filter and the ZINC stock (17,422,831 compounds), at two
budgets: default (100 iterations, 120 s) and deep (3,000 iterations, 3,600 s,
max 9 transforms — the same config the project used previously).

### Result: no route found, at either budget

| target | solved | routes explored | **solved routes** | best depth | top score | time |
|---|---|---|---|---|---|---|
| **`t4_72f5671e89cb`** (bdhi_c5) | **NO** | 1,939 | **0** | 6 | 0.798 | 327 s |
| adduct form | **NO** | 1,596 | **0** | 4 | 0.738 | 258 s |
| `t4_9c44a3f8892a` (bdhi_c4, same R-group) | **NO** | 2,248 | **0** | 8 | 0.815 | 287 s |
| `t4_45901f30d2a1` (chloroacetamide, same R-group) | **YES** | 1,760 | **460** | 2 | 0.994 | 136 s |
| **sulfopin** | **YES** | 4,741 | **316** | 3 | 0.987 | 365 s |

The controls are the value here. **Change only the warhead, keep the core and
the R-group identical, and the target goes from solved-in-two-steps to
unsolvable in 1,939 routes.**

"Not solved" is a statement about the stock, not about chemistry — a route
ending at a real-but-uncatalogued intermediate is reported unsolved. So the
useful question is *where* it stalls.

### The BDHI ring and the C–Br are formed by recognised chemistry. The C–N bond is not.

The best unsolved route for the candidate, read outward from the target:

```
target  O=S1(=O)CC[C@@H](N(Cc2cnoc2)C2CC(Br)=NO2)C1
 └─ sulfide → sulfone oxidation                    m-CPBA                  [IN STOCK]
    └─ amide → amine reduction
       └─ N-acylation                              isoxazole-4-COCl
          ├─ acid chloride                         isoxazole-4-COOH + SOCl2 [IN STOCK]
          └─ ** [3+2] nitrile-oxide cycloaddition **
             ├─ dibromoformaldoxime  ON=C(Br)Br                            [IN STOCK]
             └─ N-vinyl amine  C=CN[C@@H]1CCSC1                        [NOT IN STOCK]
                └─ carbamate cleavage ← C=CN(C(=O)OCC)[C@@H]1CCSC1
                   └─ ← diethyl sulfate + C=CN(C(=O)O)[C@@H]1CCSC1    [DEAD END]
```

**The warhead is not the problem, and this is worth saying plainly.** The policy
proposed exactly the recognised route: bromonitrile oxide, generated in situ
from **dibromoformaldoxime** (in stock), undergoing [3+2] cycloaddition across
an alkene to give the 3-bromo-4,5-dihydroisoxazole. **The BDHI ring and the
C(sp2)–Br bond are formed together, in one step, from a catalogue reagent** —
this is the standard construction for the class, and AiZynthFinder found it
without prompting. Full marks for the warhead.

**The problem is the C5–N bond — the N,O-acetal from §2.** To install the amine
at C5 by cycloaddition you need an **enamine** dipolarophile
(`C=CN[C@@H]1CCSC1`), the N-vinyl derivative of the sulfolane amine. The search
could not buy it, and its own proposal for making it terminates at
`C=CN(C(=O)O)[C@@H]1CCSC1` — a free N-vinyl **carbamic acid**, which is not an
isolable compound. Depth 6 and a dead end.

The `bdhi_c4` sibling fails the same way, deeper (depth 8), stalling on
`O=C(O)N(C1CON=C1Br)[C@@H]1CCSC1` — again the C–N bond to the dihydroisoxazole.
So **the hard bond is amine-to-dihydroisoxazole regardless of regiochemistry**;
C5 additionally makes the resulting bond an N,O-acetal.

Policy confidence was low throughout: every step on the candidate's best route
except the acid-chloride formation scored below 0.03.

**Practical reading for a chemist.** This does not say the molecule cannot be
made. It says the retrosynthetic policy, given 3,000 iterations and 17.4 M
purchasable compounds, could not reach one, and that the obstacle is a specific,
nameable bond that a chemist can either solve (reductive amination onto a C5
carbonyl or C5-halide; addition of the amine to an isoxazoline-derived
oxocarbenium; a completely different cycloaddition partner) or design around by
moving to `bdhi_c4`, which D0024 already prefers. If the route to C5-amino BDHIs
is not routine, that is a cost that should be weighed against a rank that §0
shows does not survive its own convergence control.

---

## 4. The stereocentre — and the second one nobody specified

### Is `[C@@H]` the sulfopin-matching configuration? **Yes.**

The sulfolane carbon is **(R)**, and that matches sulfopin **as crystallised**.

Note which sulfopin: `config/seeds.yaml` and
`pin1_covalent_cys113_anchors_2.csv` both store sulfopin with **no
stereochemistry at all** (`CC(C)(C)CN(C1CCS(=O)(=O)C1)C(=O)CCl`). The PDB
chemical component **QT7 in 6VAJ** — the molecule actually solved bonded to
Cys113, recorded in this project's own
`00_outputs/blacksmith/pdb_covalent/covalent_links_3.csv` — carries it:
`CC(C)(C)CN([C@@H]1CCS(=O)(=O)C1)C(=O)CCl`, which is **(R)**.

I did not settle this by comparing R/S letters, because a CIP label is a
function of substituent priorities and the substituents differ. I measured the
handedness directly — the signed volume of the (N, CH₂–S, CH₂–CH₂) triple about
the sulfolane C3 in an embedded 3D structure, which is the configuration itself
with no convention in between:

| molecule | signed volume | sign | CIP |
|---|---|---|---|
| **candidate** | −2.007 | **−** | R |
| sulfopin **QT7 (6VAJ)** | −2.101 | **−** | R |
| T_4 protected core (Me-capped) | −2.027 | **−** | R |
| candidate, C3 epimer | +2.044 | + | S |
| sulfopin, enantiomer | +2.181 | + | S |

Same sign, therefore same configuration. The T_4 protected core
(`N([1*])([2*])[C@@H]1CCS(=O)(=O)C1`) carries it too, so **every** T_4 molecule
inherits the crystallographic configuration. The project got this right and did
it deliberately.

**Cost of the other enantiomer: a starting-material substitution, not a new
route.** The configuration is set by the sulfolane amine, which is a purchased
building block — AiZynthFinder's solved route for the chloroacetamide sibling
begins at `N[C@@H]1CCS(=O)(=O)C1`. I checked the ZINC stock (17.4 M purchasable
compounds) for both enantiomers by InChIKey:

| building block | InChIKey | in stock |
|---|---|---|
| (R)-3-aminosulfolane | `OVKIDXBGVUQFFC-SCSAIBSYSA-N` | **yes** |
| (S)-3-aminosulfolane | `OVKIDXBGVUQFFC-BYPYZUCNSA-N` | **yes** |
| (R)- / (S)-3-aminothiolane | `GBNRIMMKLMTDLW-SCSAIBSYSA-N` / `-BYPYZUCNSA-N` | **yes / yes** |
| racemic 3-aminosulfolane | `OVKIDXBGVUQFFC-UHFFFAOYSA-N` | no (only the resolved forms) |

So the (S) compound costs the price difference on one starting material plus a
chiral-purity check. Given §0, a matched wrong-enantiomer control is arguably
the most informative *second* compound to make: if the (S) compound is equally
active, the recognition story is wrong regardless of what the ranking said.

### The stereocentre that was never specified, and was decided by a random seed

**The candidate has two stereocentres, and only one is specified.**
`Chem.FindMolChiralCenters(..., includeUnassigned=True)` returns
`[(5, 'R'), (13, '?')]`. Atom 13 is **C5 of the dihydroisoxazole** — the
N,O-acetal carbon, the one bonded to both the core nitrogen and the ring oxygen.
The SMILES writes it `C2CC(Br)=NO2`, with no `@`.

**It was not left free during docking. It was fixed, arbitrarily, by the
embedder.** `scripts/nac_screen.prepare_ligand` calls
`AllChem.EmbedMolecule(mol, randomSeed=0xC0FFEE)` once on the molecule as
parsed, writes one PDBQT, and every AutoDock-GPU run then samples torsions —
chirality is not a torsion. Reconstructing that embedding with the same seed:

```
embedded:  O=S1(=O)CC[C@@H](N(Cc2cnoc2)[C@H]2CC(Br)=NO2)C1     centres: (5, R), (13, R)
```

So **all 200 production runs, and both of my convergence runs, scored one
diastereomer — the (R) one at the dihydroisoxazole C5 — chosen by an ETKDG
random seed rather than by design. The other diastereomer has never been
scored.** The 7.23× belongs to a molecule the SMILES does not uniquely name.

This is not a cosmetic point. The Bookworm report §3 item 5 records that Watts
et al. 2006 found a **~50-fold** potency difference between 5-(S) and 5-(R) BDHI
stereoisomers on TG2, and that Byun 2023 tested BDHI compounds **as racemates**
because enantiopure BDHI synthesis is nontrivial. I would not map their CIP
labels onto ours — the C5 substituent sets differ, so the letters are not
guaranteed to describe the same arrangement — but the *magnitude* of the SAR
transfers as a warning: at this warhead's ring carbon, configuration is
potency-determining in the one programme that measured it.

**Three consequences, in order of how much they cost to fix:**

1. **Cheap.** Score both diastereomers. It is two runs of
   `medchem_nac_converge_one.py` and it would say whether the ranking is even a
   property of a definite molecule.
2. **Cheap and structural.** Enumerate the unspecified centre at library-build
   time, or refuse to rank a molecule with an unassigned stereocentre, rather
   than letting `EmbedMolecule` decide. Every `bdhi_c5` and `bdhi_c4` molecule
   in T_4 has this same unspecified C5 — this is 374 molecules, not one.
3. **Expensive.** If the compound is made, it will most likely be made as a
   diastereomeric mixture at that centre (the cycloaddition in §3 sets it with
   little control), and separating it is a real cost that no number in the
   pipeline currently accounts for.

---

## 5. Comparison to sulfopin

### The shared scaffold

MCS (ring-matches-ring, complete rings): **9 atoms, 9 bonds** —

```
[#6&!R]-[#7&!R]-[#6]1-[#6]-[#6]-[#16](-[#6]-1)(=[#8])=[#8]
```

i.e. **the sulfolane ring, its exocyclic tertiary nitrogen, and one carbon on
that nitrogen** — the "protected core" that T_4 fixes by construction, and
nothing more. ECFP4 Tanimoto **0.250**.

### What differs

| | sulfopin | candidate |
|---|---|---|
| **warhead** | chloroacetamide, `N–C(=O)CH₂Cl` | BDHI, `N–C5H(ring)`, C3–Br on a C=N |
| **mechanism** | SN2 at sp3 C, backside | addition–elimination at sp2 C, perpendicular |
| **core–warhead bond** | **amide** (planar, robust) | **N,O-acetal** (hydrolytically labile) |
| **other N substituent** | neopentyl (`CH₂C(CH₃)₃`), aliphatic | isoxazol-4-ylmethyl, aromatic heterocycle |
| aromatic rings | 0 | 1 |
| fraction sp3 | 0.91 | 0.64 |
| stereocentres | 1 (defined in the crystal) | 2 (**one undefined**) |
| validation status | 9 crystallographic Cys113 positives, class AUC 0.908, measured k = 0.028 M⁻¹s⁻¹ | **0 crystallographic positives, class cannot be validated** |

### Head to head, on everything computed

| property | candidate | sulfopin | better |
|---|---|---|---|
| MW | 364.22 | 281.81 | sulfopin |
| heavy atoms | 20 | 17 | sulfopin |
| cLogP | 1.12 | 1.29 | ≈ |
| TPSA | 85.0 | 54.45 | sulfopin |
| HBA | 7 | 3 | sulfopin |
| **QED** | **0.796** | 0.733 | **candidate** |
| **SAscore** | 4.50 | **3.19** | **sulfopin** |
| Lipinski / Veber / Ghose | pass / pass / pass | pass / pass / pass | tie |
| PAINS / BRENK / NIH (whole) | 0 / 0 / 1 | 0 / 1 / 2 | candidate (see §1 caveat) |
| N,O-acetal | **1** | 0 | **sulfopin** |
| **retrosynthesis solved** | **NO** (0 of 1,939 routes) | **YES** (316 solved, 3 steps) | **sulfopin** |
| **hERG** | 0.134 | 0.174 | candidate |
| **AMES** | **0.945 (98th pct)** | 0.729 (92nd) | **sulfopin** |
| **DILI** | **0.761** | **0.319** | **sulfopin** |
| carcinogenicity | 0.502 | 0.404 | sulfopin |
| Caco-2 | −5.06 | −4.33 | sulfopin |
| PAMPA / HIA / bioavailability | 0.845 / 0.999 / 0.838 | 0.751 / 0.942 / 0.731 | candidate |
| BBB | 0.941 | 0.776 | depends (see §1) |
| CYP inhibition (max of 5) | 0.141 | 0.173 | ≈ (both clean) |
| clearance, hepatocyte | 15.9 | 60.8 | candidate |
| plasma protein binding | 44.6% | 44.3% | tie |
| **NAC enrichment, 200 runs** | **7.23×** (6.74× re-measured) | 2.54× | candidate |
| **NAC enrichment, 2,000 runs** | **0.60×** | not measured at 2,000 | **neither — candidate is below chance** |
| ligand efficiency | 0.403 (AD4 reactive) / 0.251 (Vina adduct) | 0.446 (MM-GBSA, 50/80 decoys better) | not comparable |
| **crystallographic Cys113 evidence** | **none** | **6VAJ, 1.42 Å** | **sulfopin** |

The candidate wins on drug-likeness (QED, clearance, passive permeability
predictors) and on the single unvalidated geometric metric at a fixed run count.
Sulfopin wins on everything that has been *measured* — synthesizability, hepatic
safety prediction, mutagenicity prediction, and the existence of a crystal
structure.

---

## 6. What I would do before committing this to synthesis

In cost order. Items 1–3 are cheap and would change the answer.

1. **Score both diastereomers at the unspecified C5** (§4). Two runs. If they
   differ materially, the ranking never described a definite molecule.
2. **Run the 2,000-run convergence control across the whole BDHI class**, not
   just the three molecules in §0. The measured collapse (6.74× → 0.60×) is one
   molecule; if it is general, the D0067 fix corrected a defect and did **not**
   produce a shortlist, and the top of the list needs re-reading.
3. **Fix `nac_rank._ids_in` to filter on `status == "ok"`** (§0) and add a test.
   Five molecules are currently un-refinable through the refine path because a
   failed row counts as done.
4. **Take the question of the C5-amino BDHI to a chemist before anything else.**
   Two independent computations point at the same bond: the N,O-acetal SMARTS
   flags it as a hydrolytic liability present in 198/198 `bdhi_c5` molecules and
   surviving into the adduct, and AiZynthFinder dead-ends on exactly that
   disconnection while forming the warhead itself from a catalogue reagent in
   one step. If a chemist says the C5-amino linkage is fine and routine, both
   flags relax. If not, `bdhi_c4` — which D0024 already prefers, avoids the
   acetal entirely, and has an identical formula — is the obvious substitute.
5. **Measure hydrolytic stability in buffer before any assay.** Byun 2023
   reports at least one BDHI analogue excluded for being unstable in phosphate
   buffer (Bookworm §3 item 2); the N,O-acetal gives this molecule a second,
   independent route to the same failure.
6. **A real Ames**, given 0.945 / 98th percentile — and a matched-pair Ames on
   the `bdhi_c4` sibling and the chloroacetamide sibling, which are the natural
   controls.

**Handoffs.** The Adversary should audit §0's convergence measurement (n = 3,
one GPU, one seed each — it wants replicates before it is quoted as settled) and
the `nac_rank` resume defect. The Artist has the tables for a
warhead-comparison figure. The Bookworm's report already covers the literature
and should be read alongside this one; the two agree on the stereocentre and on
the absence of a Pin1–BDHI structure, and were written independently.

---

## Provenance and reproduction

Code (this repo):

- [`scripts/medchem_admet.py`](../scripts/medchem_admet.py) — ADMET-AI driver + the featurizer guard
- [`scripts/medchem_workup.py`](../scripts/medchem_workup.py) — descriptors, rule sets, alerts, LE/LLE, stereochemistry, handedness, reference-set comparison
- [`scripts/medchem_retro_report.py`](../scripts/medchem_retro_report.py) — AiZynthFinder output summariser
- [`scripts/medchem_nac_converge_one.py`](../scripts/medchem_nac_converge_one.py) — per-molecule convergence check

Outputs (`/data/lab_vm/append_only/inhibition/00_outputs/blacksmith/medchem_t4_72f5671e89cb/`):

| file | contents |
|---|---|
| `admet_ai_1.csv` | ADMET-AI, candidate + adduct + sulfopin |
| `admet_ai_2.csv` | ADMET-AI, `bdhi_c4` + chloroacetamide siblings |
| `workup_2.json` | every physchem / alert / stereo / LE number in §2, §4, §5 |
| `aizynth_default_1.json.gz`, `aizynth_deep3000_1.json.gz` | retrosynthesis, candidate + adduct + sulfopin |
| `aizynth_siblings_deep3000_1.json.gz` | retrosynthesis, both siblings |
| `retro_summary_1.json`, `retro_summary_2.json` | parsed routes, leaves, stock status |
| `nac_converge_one_1.csv` | **the 200 vs 2,000-run measurement in §0** |
| `aizynth_config_default_1.yml`, `aizynth_config_deep_1.yml` | search configs |
| `smiles_medchem_1.smi`, `smiles_siblings_1.smi` | inputs |

Reference-set hashes: master `4c8f7d1be94bdf55dda4b641e9512cbeb596d7806a3a80fe1492fc339a25ff89`,
anchors `d832eda3b4f70d224c66c83c86af0872862f9f4e8ee69a4769e4a47bc77bcef2`.
Candidate frame: `04_t4_combinatorial/D4_43.parquet`, row `t4_72f5671e89cb`.
NAC ranking read with version-order deduplication, `keep="last"` (D0067).

To reproduce:

```bash
# ADMET (CPU is fine; pin a free GPU if one is used)
CUDA_VISIBLE_DEVICES=7 nice -n 19 /data/lab_vm/envs/dwi_admet/bin/python \
  scripts/medchem_admet.py --smiles <file>.smi --outdir <outdir> --version N

# physchem / alerts / stereo / LE
nice -n 19 /data/lab_vm/envs/dwi_cheminf/bin/python scripts/medchem_workup.py <outdir> N

# retrosynthesis  (the dwi_retro env's python is not readable by other users;
# ~/.local/pin1_tools/aizynth_venv is an equivalent AiZynthFinder 4.4.1 and
# reads the same model + stock files)
V=~/.local/pin1_tools/aizynth_venv
PATH=$V/bin:$PATH nice -n 19 $V/bin/aizynthcli --smiles <file>.smi \
  --config aizynth_config_deep_1.yml --output <out>.json.gz --nproc 3

# convergence check  (needs meeko + gemmi: dwi_reactive, not dwi_cheminf)
nice -n 19 ~/.micromamba/envs/dwi_reactive/bin/python \
  scripts/medchem_nac_converge_one.py --idents t4_72f5671e89cb \
  --nrun 200 2000 --gpu 7 --outdir <outdir> --version N
```
