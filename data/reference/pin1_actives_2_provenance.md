---
title: Provenance and chemotype assessment for expanded Pin1 actives set (v2)
date: 2026-07-27
project: dance_with_inhibition
sensitivity: standard
tags: [pin1, docking, actives-set, chembl, provenance, enrichment-control]
sources: ['@bookworm']
url: https://www.ebi.ac.uk/chembl/api/data/target/CHEMBL2288
---

# Pin1 reference binders v2 — provenance and skeptical assessment

Companion file to `pin1_reference_binders_2.csv`. Assembled to expand the
ACTIVES set for the "Dance with Inhibition" docking-enrichment control
(repo `inhibition`). v1 had n=11 rows (6 claimed covalent, 5 claimed
non-covalent); this version corrects a target-identity error found in v1
and adds 8 new, independently-sourced actives spanning several new
chemotypes.

## 0. Critical correction to v1 — READ FIRST

**`CHEMBL3391` is not human Pin1.** I queried it directly:

```
target_chembl_id: CHEMBL3391
pref_name:        Threonine--tRNA ligase 1, cytoplasmic (TARS1)
accession:        P26639
```

The correct ChEMBL target for human Pin1 (peptidyl-prolyl cis-trans
isomerase NIMA-interacting 1, UniProt `Q13526`) is **`CHEMBL2288`**.

v1 rows 11–12 ("reversible-Thr-sulfonamide-quinazoline",
"reversible-Thr-sulfamate-adenosine") cited `CHEMBL2311920` and
`CHEMBL1163068` "target CHEMBL3391". I traced both molecules' full
activity records:

- `CHEMBL2311920`: Ki 1.2 nM vs `CHEMBL3391` (TARS1) and vs
  `CHEMBL2311235`("Threonine--tRNA ligase"), plus antibacterial MIC
  data (*E. coli*, *H. influenzae*, *B. thailandensis*) — this is a
  **bacterial/human threonyl-tRNA synthetase inhibitor**, structurally
  and functionally unrelated to Pin1.
- `CHEMBL1163068`: IC50 15 nM vs `CHEMBL3391` (TARS1), plus a PBMC
  cytotoxicity assay and a second ThrRS record (`CHEMBL5465377`) —
  same story.

Neither compound has **any** activity record against `CHEMBL2288`
(Pin1) in ChEMBL. **Both rows are removed in v2.** This means v1's
claimed "5 non-covalent actives" was really only **3 valid non-covalent
Pin1 actives** (ATRA, EGCG, PiB) — and two of those three are already
flagged promiscuous. This makes the expansion in this file more
consequential than it might have looked: the non-covalent side of the
control was nearly empty of clean data.

## 1. Sources queried

- **ChEMBL REST API** (`www.ebi.ac.uk/chembl/api/data`), target
  `CHEMBL2288`: pulled all 661 quantitative activities
  (`pchembl_value` not null, Homo sapiens, IC50/Ki/Kd/EC50), filtered
  to ≤10 µM potency (479 records / 394 unique molecules), then
  triaged by document to find independent chemotypes. Raw JSON
  retained at `/tmp/.../scratchpad/chembl2288_quant.json` (session
  scratch, not committed).
- **ChEMBL document endpoint**: pulled title/authors/journal/DOI/PMID
  for every document cluster considered, to attribute compounds to
  primary literature rather than leaving them as bare ChEMBL IDs.
- **PubChem BioAssay** (`pug/bioassay/target/genesymbol/PIN1/aids`):
  target search returned 130+ AIDs. Spot-checked several; the sample
  (e.g. AID 1904, a GNF genome-wide **siRNA** circadian-rhythm screen
  that happens to include *PIN1* as a hit gene, not a chemical-binding
  assay) was predominantly off-target/phenotypic/gene-level data, not
  independent chemical Pin1-binder evidence beyond what ChEMBL already
  had. I did **not** exhaustively triage all 130+ AIDs — that is a
  scope limitation, noted honestly rather than silently skipped.
- **PubMed / WebSearch / WebFetch**: read abstracts for the four most
  recent (2023–2025) primary papers behind new chemotypes, to confirm
  mechanism claims rather than trust ChEMBL's bare `standard_type`
  label.
- **ChEMBL document-by-DOI** and **PubChem xref-by-DOI**: attempted
  to resolve the two pre-existing `UNVERIFIED` rows (Reddi 2023 JACS
  `10.1021/jacs.2c08853`; Byun 2023 JACS `10.1021/jacs.3c00598`).
  Neither DOI has a matching ChEMBL document, and PubChem's PUG-REST
  does not support DOI cross-reference lookup the way I tried it.
  **Left both rows `UNVERIFIED`, unchanged** — no invented SMILES.

## 2. New rows added (v2), with provenance

All SMILES below were parsed and canonicalized with RDKit
(`/data/lab_vm/envs/dwi_cheminf/bin/python`, RDKit `Chem.MolFromSmiles`
+ `MolToSmiles`) before being written to the CSV — 15/15 real
structures parsed cleanly (the 2 `UNVERIFIED` rows carried forward
from v1 are not structures and were not run through RDKit).

| name | ChEMBL mol ID | mechanism | potency | primary source |
|---|---|---|---|---|
| Wildemann-macrocyclic-peptide | CHEMBL219531 | non_covalent | Ki 1.2 nM | Wildemann 2006, *J Med Chem*, PMID 16570909, DOI 10.1021/jm060036n |
| Guo-Pfizer-benzothiophene-phosphonate | CHEMBL585917 | non_covalent | Ki 6 nM | Guo 2009, *Bioorg Med Chem Lett*, PMID 19729306, DOI 10.1016/j.bmcl.2009.08.034 |
| Potter-Astex-indole-furancarboxamide | CHEMBL595291 | non_covalent | IC50 25 nM | Potter 2010, *Bioorg Med Chem Lett*, PMID 19969456, DOI 10.1016/j.bmcl.2009.11.090 |
| Liu-Pei-cyclic-peptide | CHEMBL1089622 | non_covalent | IC50 31 nM / Kd 47 nM | Liu 2010, *J Med Chem*, PMID 20180533, DOI 10.1021/jm901778v |
| Jiang-Pei-bicyclic-CPP-peptide | CHEMBL3590308 | non_covalent | Kd 72 nM | Jiang & Pei 2015, *J Med Chem*, PMID 26196061, DOI 10.1021/acs.jmedchem.5b00411 |
| Du-Xu-naphthalenecarboxamide | CHEMBL1090715 | non_covalent | Ki 6 nM | Du 2021, *Bioorg Med Chem*, PMID 33246256, DOI 10.1016/j.bmc.2020.115878 |
| Tian-chloropyrimidine-covalent-6a | CHEMBL6152351 | covalent_cys113 | IC50 3.15 µM (X-ray confirmed) | Tian 2025, *ACS Med Chem Lett*, PMID 39811131, DOI 10.1021/acsmedchemlett.4c00477 |
| Pu-benzylguanine-API-series | CHEMBL6120754 | unknown | IC50 64.4 nM | Pu 2025, *J Med Chem*, PMID 39868498, DOI 10.1021/acs.jmedchem.4c02144 |

Every row's `target_organism` in the underlying ChEMBL activity record
is `Homo sapiens` against `CHEMBL2288` — verified per-row, not assumed
from the target-level search alone.

### Notes on selection reasoning (best record chosen per molecule)

- **Guo-Pfizer-benzothiophene-phosphonate (CHEMBL585917)** is the
  single best-corroborated compound in the whole set: it was
  independently re-assayed and re-reported potent (6–50 nM) in **five**
  separate papers spanning 2009–2021 (Guo 2009 original; Moore 2013
  review; two intervening Bioorg Med Chem papers in 2016/2018; Du
  2021). That kind of independent cross-lab reproduction is exactly
  the "known/validated" bar Bookworm should be applying, and it fills
  a real gap: after removing the two mistargeted rows, v1 had **zero**
  clean non-peptidic, non-covalent small-molecule Pin1 binders.
- **Wildemann-macrocyclic-peptide (CHEMBL219531)**: I found a second,
  near-identical compound (`CHEMBL5395865`, same Ki 1.2 nM, 1–2 atom
  differences) attributed **only** to a 2023 review document
  (`CHEMBL5325477`, He et al., *J Med Chem*, a survey of Pin1 biology
  and inhibitors) with no independently traceable original citation.
  ChEMBL's review-document activities turned out to be systematically
  useful for **discovery** (see §3) but not reliable as **primary
  provenance** — several are simply re-transcriptions of compounds
  from other papers under the review's own document ID. I used the
  review to find leads, then always re-attributed to the original
  paper where traceable, and excluded `CHEMBL5395865` specifically
  because I could not find its original source.
- **Liu-Pei-cyclic-peptide** and **Jiang-Pei-bicyclic-CPP-peptide** are
  from the same lab (Dixon "Pei" lab) and share the same phosphonate
  Ser/Thr-mimetic warhead logic as the Wildemann macrocycle — see the
  chemotype-clustering discussion below for how I handled this
  relatedness honestly rather than silently counting all three as
  fully independent.
- **Jiang-Pei-bicyclic-CPP-peptide** is MW ≈ 2372 Da (RDKit-confirmed
  exact mass). I flagged this directly in the CSV `citation` field:
  it is real, wet-lab-validated affinity data, but a molecule this
  large and this flexible may not be a fair test case for a
  standard small-molecule docking-enrichment pipeline. That's a call
  for the blacksmith/adversary, not me — I'm surfacing it, not
  deciding it.
- **Du-Xu-naphthalenecarboxamide (CHEMBL1090715)**: the source paper's
  title is about *thiazole-based* Pin1 inhibitors, but the specific
  SMILES ChEMBL attributes to that document has no thiazole ring in
  it. I could not fully resolve whether this is the paper's flagship
  compound or a reference/intermediate compound reused from earlier
  SAR — the Ki 6 nM record is real and traceable to a genuine
  biochemical Pin1 assay in that document, so I kept it, but flagged
  the discrepancy explicitly in the CSV rather than silently
  presenting it as "the thiazole lead."
- **Tian-chloropyrimidine-covalent-6a (CHEMBL6152351)**: I confirmed
  via the paper's abstract that this is "compound 6a," IC50 3.15 µM,
  with an **X-ray cocrystal structure showing a covalent Cys113
  adduct** via 2-chloro-5-nitropyrimidine SNAr chemistry — genuinely
  rigorous validation. But its recognition fragment
  (`CC(C)(C)CN(...)C1CCS(=O)(=O)C1`, tert-butyl + sulfolane) is
  **identical** to Sulfopin's. This is a warhead-swap analog of
  Sulfopin, not an independent chemotype — I say so explicitly in the
  CSV so it cannot silently inflate the chemotype count. (The same
  document also independently re-confirms the parent Sulfopin
  structure/potency, `CHEMBL5397880`, Ki 17 nM — a nice, if redundant,
  independent replication of the existing v1 row.)
- **Pu-benzylguanine-API-series (CHEMBL6120754)**: mechanism is
  genuinely uncertain. I could not get full-text access to confirm
  whether the "6-O-benzylguanine" naming reflects an intended
  AGT-like covalent transfer mechanism onto Cys113, or is coincidental
  scaffold nomenclature for a non-covalent PPIase-domain binder found
  by virtual screening (as a secondary source claims for the
  predecessor compound API-1, IC50 72 nM, itself not independently
  verified here). I marked `mechanism: unknown` rather than guess —
  this is exactly the kind of ambiguity the task asked me to surface,
  not silently resolve.

## 3. Exclusions — and why

| Item | Reason excluded |
|---|---|
| `CHEMBL2311920`, `CHEMBL1163068` (v1 rows 11–12) | **Wrong target.** Both bind Threonine-tRNA ligase (TARS1, CHEMBL3391), not Pin1 (CHEMBL2288). Zero Pin1 activity records exist for either in ChEMBL. Removed entirely — do not use these SMILES as Pin1 actives under any circumstances. |
| `CHEMBL5395865` (peptide near-identical to Wildemann's) | Only traceable to a 2023 review document, not an independently verifiable original paper. Provenance too thin to trust as distinct primary data. |
| `CHEMBL4762559`, `CHEMBL4526940` (naphthoquinone-sulfonylacetate E/Z isomers) | Same chemotype as the already-included KPT-6566 (self-immolative aryl-sulfonyl-acetate → naphthoquinone); adding them would just re-inflate the quinone-warhead bias the task explicitly wants reduced. |
| Organoselenium diesters, e.g. `CHEMBL5411977` (`CCOC(=O)C([Se]c1ccccc1)([Se]c1ccccc1)C(=O)OCC`) and analogs from the He 2023 review table | **Excluded as likely assay artifacts.** Organoselenium compounds are well-documented nonspecific thiol-reactive / redox-active species. Against an enzyme whose mechanism-of-action literature is centered on a single reactive active-site cysteine (Cys113), a diselenide/selenide ester binding with modest potency (440–780 nM in this dataset) is a textbook false-positive risk, not a validated chemotype. No independent replication or selectivity data found. |
| Indanone/benzophenone-skeleton nitrophenol-catechol series, e.g. `CHEMBL2017129` (Liu 2012, *Bioorg Med Chem*, PMID 22459212) | Considered but not promoted to primary picks. Nitrophenol/catechol-bearing fragments carry classic redox-cycling/PAINS liability; potency in this series was modest (Ki 200 nM–50 µM, mostly weaker end) with no independent replication found across other documents. Kept out of the actives set; noted here in case the user wants a separate "low-confidence / caution tier" later. |
| PiB, ATRA, EGCG, juglone, KPT-6566, Sulfopin, BJP-06-005-3 re-appearing under other ChEMBL document IDs (`CHEMBL43612`, `CHEMBL38`, `CHEMBL2409076`, `CHEMBL5397880`, `CHEMBL6148397`) | Not duplicated — these are the *same* molecules already in the file (SMILES-confirmed identical), re-cited by later review/comparator papers. Their reappearance is a useful **independent cross-validation** of the v1 entries (particularly reassuring for Sulfopin and PiB), not new data. |
| Reddi 2023 (JACS), Byun 2023 (JACS) — pre-existing `UNVERIFIED` rows | Attempted resolution via ChEMBL document-by-DOI and PubChem xref-by-DOI; both came back empty. Left as `UNVERIFIED`, unchanged, per the no-invented-SMILES rule. |
| ~120 of the ~130 PubChem AIDs tagged to gene symbol PIN1 | Not exhaustively mined (scope/time). Spot-checks suggested predominantly off-target or phenotypic/genomic screens (e.g. a GNF siRNA circadian-rhythm screen), not independent chemical-binder evidence. Flagged as a genuine limitation, not claimed as "checked and clean." |

## 4. Independent chemotype assessment (the honest count, not just headcount)

I clustered every structurally-resolved compound in v2 (15 of 17 rows;
the 2 `UNVERIFIED` rows can't be clustered without a structure) by
scaffold/pharmacophore family:

| Cluster | Members | Independent? |
|---|---|---|
| A1 — tert-butyl-sulfolane amine + electrophile | Sulfopin (chloroacetamide), Tian-6a (chloro-nitropyrimidine) | **1 chemotype**, 2 warhead variants — do not double-count |
| A2 — Trp-Pro-Phe-Arg peptidomimetic + N-Me-chloroacetamide | BJP-06-005-3 | independent |
| B — self-immolative aryl-sulfonyl-acetate → naphthoquinone | KPT-6566 | independent |
| C — simple 1,4-naphthoquinone Michael acceptor | Juglone | independent (distinct warhead kinetics from B, both "quinone" only in a loose sense) |
| D — retinoid | ATRA | independent (promiscuity-flagged) |
| E — flavan-3-ol gallate | EGCG | independent (promiscuity-flagged) |
| F — naphthalimide/pyromellitic diimide diester | PiB | independent (promiscuity-flagged) |
| G — disulfide-bridged bicyclic peptide + phosphonate | Wildemann-macrocyclic-peptide | independent library origin |
| J — head-to-tail cyclic hexapeptide + phosphonate | Liu-Pei-cyclic-peptide | related to G/K by pharmacophore logic (peptide macrocycle + phospho-mimetic), distinct macrocyclization chemistry/lab of origin |
| K — bicyclic CPP-fused macrocycle + phosphonate | Jiang-Pei-bicyclic-CPP-peptide | related to G/J (same lab lineage, follow-up design), but structurally much larger/different (added CPP tail + aromatic diamide bridge) |
| H — benzothiophene-carboxamide phosphonate | Guo-Pfizer-benzothiophene-phosphonate | independent, non-peptidic |
| I — indole-alanine + furancarboxamide, phosphate-free | Potter-Astex-indole-furancarboxamide | independent |
| L — naphthalenecarboxamide-cinnamylalanine, COOH warhead-mimic | Du-Xu-naphthalenecarboxamide | independent |
| M — 6-O-benzylguanine purine ether | Pu-benzylguanine-API-series | independent |

**Strict count** (merging G/J/K into one "phosphonate-peptide-macrocycle"
superfamily, and A1 counted once): **11 independent chemotypes** across
15 structurally-known actives.

**Loose count** (treating G, J, K as 3 distinct entries because their
macrocyclization chemistry and library of origin genuinely differ, even
though the core pharmacophore logic — phosphonate mimicking pSer/Thr-Pro
— is shared): **13 independent chemotypes**.

I'd report this honestly as **"11–13 independent chemotypes depending
on how strictly you split the phosphonate-peptide-macrocycle family,"**
not a single flattering number. Compare to v1, which — once the two
mistargeted rows are removed — really only had **7 independent
chemotypes** (A1{Sulfopin}, A2{BJP}, B{KPT-6566}, C{Juglone}, D{ATRA},
E{EGCG}, F{PiB}) behind 9 structurally-known compounds, 2 of which
were the *same* two chemotypes the user specifically complained about
(2 Sulfopin-family, 2 quinones). v2 roughly **doubles** both the
compound count and the independent-chemotype count, and does not add
any new quinone or Sulfopin-fragment representatives to the primary
picks (Tian-6a is flagged, not silently counted, as a Sulfopin analog).

## 5. Compounds whose reported Pin1 activity I consider doubtful

Answering the "assessment I specifically want" directly, beyond what's
already in the existing file's promiscuity flags (ATRA, EGCG, juglone,
KPT-6566, PiB are already flagged `y` in the schema and correctly so —
I did not change those):

1. **The organoselenium diesters** (§3) — doubtful as genuine Pin1
   engagement; more likely nonspecific thiol reactivity. Not included.
2. **`CHEMBL5395865`** (peptide near-identical to Wildemann's compound)
   — doubtful *provenance*, not doubtful *chemistry*; I can't verify
   who actually measured this Ki. Not included.
3. **Pu-benzylguanine-API-series mechanism claim** — I am not
   confident the compound is covalent despite the suggestive name; I
   also cannot confirm this exact structure is literally "API-32" as
   opposed to a related SAR compound ChEMBL attributes to the same
   document. Included, but both caveats are stated in the CSV
   citation field itself so they travel with the data.
4. **Du-Xu naphthalenecarboxamide's "thiazole" attribution** — the
   paper title doesn't match the deposited structure's substructure.
   Included with the discrepancy flagged.

None of these four rises to "exclude entirely" the way the
mistargeted TARS1 compounds or the organoselenides did, but all four
carry a caveat that should survive into any downstream docking
enrichment write-up — an enrichment control is only as trustworthy as
its weakest-provenance active.
