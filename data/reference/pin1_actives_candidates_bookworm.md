---
title: Candidate Pin1 actives for chemotype-diversity expansion (bookworm literature sweep)
date: 2026-07-30
project: dance_with_inhibition
sensitivity: standard
tags: [pin1, docking, actives-set, chembl, chemotype-diversity, enrichment-gate, provenance]
sources: ['@bookworm']
url: https://www.ebi.ac.uk/chembl/api/data/target/CHEMBL2288
related: ['[[pin1_reference_binders_2.csv]]', '[[pin1_covalent_cys113_anchors_2.csv]]', '[[warhead_classes_7.csv]]', '[[D0040-the-residual-auc-was-never-significant]]']
---

# Candidate Pin1 actives — chemotype-diversity expansion

Request: find published Pin1 actives that clear a genuinely **new** warhead/
scaffold chemotype (not already in `pin1_reference_binders_2.csv`'s 4
covalent chemotypes: chloroacetamide, naphthoquinone, sulfamate acetamide,
BDHI), are small molecules (no peptidic macrocycles), have a resolvable
SMILES, and a named-assay potency. Target: human Pin1 (Q13526), PPIase
domain, Cys113, PDB 6VAJ. ChEMBL target confirmed as **CHEMBL2288** (not
CHEMBL3391 — that is threonine-tRNA ligase, per this project's prior
correction in `pin1_actives_2_provenance.md`; re-confirmed directly against
the ChEMBL target endpoint before running any query below).

**Bottom line up front:** I found **one** genuinely new *covalent*
chemotype and **one** genuinely new *non-covalent* chemotype, both with
resolvable structures and named-assay potency. Two chloroacetamide analogs
worth keeping as anchors turned up too, but they do **not** count toward
chemotype diversity — same warhead as Sulfopin/BJP-06-005-3, just a
different recognition scaffold. Whether the non-covalent hit moves your
gate's chemotype counter at all depends on whether that counter is
warhead-keyed (as `warhead_classes_7.csv` suggests) or active-keyed — see
the note at the end.

## Method

1. Re-pulled all 661 quantitative Pin1 activities from ChEMBL
   (`target_chembl_id=CHEMBL2288`, `pchembl_value` not null), raw JSON at
   `/tmp/.../scratchpad/chembl2288_activities.json` (session scratch, not
   committed) — 396 unique molecules at pIC50/pKi/pKd ≥ 5.0 (≤10 µM).
2. Resolved all 35 source documents behind those activities to
   title/journal/year/DOI, and triaged every one **I had not already seen**
   in `pin1_actives_2_provenance.md` — 9 documents new to this project
   (2019–2025), by potency and by SMILES warhead pattern.
3. Cross-checked warhead assignments against `warhead_classes_7.csv`'s
   existing SMARTS definitions to confirm "new" vs "same chemotype, new
   scaffold."
4. Ran targeted WebSearch/WebFetch passes for: covalent-fragment /
   chemoproteomics papers with PIN1_C113 as a hit (Cravatt/Backus-style
   cysteinome datasets), boronic acid / sulfonyl fluoride / vinyl sulfone /
   epoxide / nitrile / aldehyde electrophiles at Cys113, and recent
   PDB depositions.
5. Verified every candidate SMILES parses in RDKit and computed MW/formula
   as a sanity check (peptidic-macrocycle screen + gross-error catch).

## Candidates table

| name | canonical_smiles | chemotype/warhead_class | covalent? | Cys113 engagement evidence | potency + assay | PDB | DOI or PMID | promiscuity/PAINS concern |
|---|---|---|---|---|---|---|---|---|
| **Ieda-(S)-2** (ChEMBL `CHEMBL4467081`) — **NEW chemotype** | `O=C(/C=C/c1ccc2ccccc2c1)N[C@@H](CCC(=O)N1CCOCC1)C(=O)O` | cinnamamide / aryl Michael acceptor ((E)-3-(naphthalen-2-yl)acrylamide conjugated to a glutamine-morpholinamide amino acid) | yes | ESI-MS-confirmed Cys113 adduct (intact-protein MS, Ieda et al.) | IC50 3.2 µM (α-chymotrypsin protease-coupled PPIase assay); Ki 1.37 µM; kinact 3.42×10⁻⁷ s⁻¹; kinact/Ki 0.249 M⁻¹s⁻¹ — comparable-to-better than Sulfopin's 0.028 M⁻¹s⁻¹ | none | Ieda et al. 2019, *Bioorg Med Chem Lett* 29:353-356, DOI 10.1016/j.bmcl.2018.12.044, PMID 30585173; ChEMBL doc `CHEMBL4371086` | moderate — cinnamamide/styryl-amide Michael acceptors are a known reactive pharmacophore family (structurally adjacent to PAINS-flagged cinnamaldehydes); this single paper did **not** run a proteome-wide chemoproteomic selectivity panel (unlike Sulfopin/Reddi), so off-target reactivity is unmeasured, not just "clean." Flag accordingly. |
| Ieda-6 (AM-ester prodrug of (S)-2) (`CHEMBL4518987`) — same chemotype as above, not independent | `CC(=O)OCOC(=O)[C@H](CCC(=O)N1CCOCC1)NC(=O)/C=C/c1ccc2ccccc2c1` | cinnamamide (acetoxymethyl-ester cell-permeable prodrug of Ieda-(S)-2) | yes (after intracellular esterase cleavage to (S)-2) | inferred — same warhead as (S)-2, not independently MS-confirmed as the prodrug | cytotoxic in PC-3 (20 µM MTT), HCT116 (28 µM), suppresses cyclin D1; no direct enzymatic IC50 for the ester itself | none | same as above | same caveat as (S)-2, plus prodrug-specific off-target esterase dependence |
| **Liu-2024-C3** (ChEMBL `CHEMBL5589996`) — **NEW chemotype** | `N#Cc1ccc(-c2ccc(-n3cc(CC(NC(=O)C4(c5ccccc5)CC4)C(=O)NC4CC(O)C4)nn3)cc2)cc1` | 1-arylcyclopropane-1-carboxamide + triazolylalanine + biphenyl-4-carbonitrile — a "nonacidic"/neutral (no phosphate/phosphonate/carboxylate) DEL-derived scaffold | **no** (non-covalent) | not applicable (non-covalent); functional evidence for **active-site** (not WW-domain) engagement: IC50 measured as competitive displacement of a canonical PPIase-domain phosphopeptide substrate probe (TAMRA-Bth-D-pThr-Pip-L-2Nal — the same substrate-mimetic pharmacophore family as the Guo/Wildemann/Liu-Pei phosphonate binders already in your reference set), by fluorescence polarization | Kd 130 nM (SPR, GST-Pin1); IC50 660 nM (FP substrate-competition assay) | **9INR** — co-crystal structure of full-length Pin1 with this exact compound ("C3" / ligand ID `A1D9K`, formula C32H30N6O3 — matches this SMILES) | Liu C, et al. 2024, *J Med Chem* 67(17):15780-15795, DOI 10.1021/acs.jmedchem.4c01412, PMID 39229909 | low-moderate — aryl nitrile + triazole are not classical PAINS motifs, but this is a DEL-derived hit (higher false-positive base rate for the modality generally); the paper's own headline finding is that PIN1 inhibition (including PROTAC-mediated degradation) showed **no meaningful antiproliferative effect** and siRNA knockdown gave unfavorable evidence for PIN1 as an oncology target — doesn't undermine the *binding* data, but worth knowing before citing this paper for biological relevance |
| Liu-2024-C10 (`CHEMBL5591180`) — same chemotype as C3, most potent analog, not independent | `Cc1ccccc1C1(C(=O)NC(Cc2cn(-c3ccc(-c4ccc(C#N)cc4)cc3)nn2)C(=O)NC2CCN(C)CC2)CC1` | same as C3 (2-methylphenyl-cyclopropanecarboxamide variant, N-methylpiperidinyl amide tail) | no | inferred from C3's assay family (same FP substrate-competition assay used across the series) | Kd 25 nM (SPR); IC50 150 nM (FP substrate-competition) — the most potent compound in the whole 396-molecule ChEMBL2288 potency pull I ran, after the sub-nM Wildemann/Guo/Liu-Pei phosphonate peptidomimetics | none (C3 is the crystallized analog; C10 has no deposited co-structure) | same document as C3 | same caveats as C3 |

### Not new chemotypes — noted so you don't re-search this ground

| name | canonical_smiles | chemotype/warhead_class | covalent? | Cys113 engagement evidence | potency + assay | PDB | DOI or PMID | promiscuity/PAINS concern |
|---|---|---|---|---|---|---|---|---|
| ZL-Pin13 (`CHEMBL5075017`) | `O=C(CCl)N1CCC2(CC1)SCC(=O)N2Cc1ccc(-c2cccc3ccccc23)o1` | **chloroacetamide** on a spiro-thiazolidinone-piperidine + furyl-naphthyl scaffold — different recognition fragment from Sulfopin's tert-butyl-sulfolane, but the identical reactive electrophile | yes | X-ray-confirmed Cys113 adduct; co-crystal shows a Gln129 conformational change induced by the inhibitor | IC50 0.067±0.03 µM (cell-active; culminated from hit ZL-Pin01 IC50 1.33±0.07 µM across a structure-guided SAR series) | **7F0M** (ZL-Pin13); series also deposited 7EFJ (ZL-Pin01), 7EFX (ZL-Pin03), 7EKV (ZL-Pin05) | Liu L, et al. 2022, *J Med Chem* 65(3):2174-2190, DOI 10.1021/acs.jmedchem.1c01686; ChEMBL doc `CHEMBL5038660` | low reported here, but same chloroacetamide reactivity class as Sulfopin — worth using as a **second, independently-discovered chloroacetamide anchor** (different recognition scaffold) if your project ever wants within-chemotype scaffold diversity for chloroacetamide docking/decoy work, but it does **not** move the chemotype counter |
| Triazine-DEL-covalent (`CHEMBL5207244`) | `CNc1nc(Nc2ccc3cn[nH]c3c2)nc(N2CC3CC2CN3C(=O)CCl)n1` | **chloroacetamide** appended to a triazine DNA-encoded-library scaffold | yes (by DEL-selection design; not independently MS-confirmed for this specific off-DNA resynthesized hit in the accessible abstract) | not confirmed beyond covalent-DEL selection pressure | IC50 1.69 µM | none | Li L, et al. 2022, *ACS Med Chem Lett* 13(10):1574-1581, DOI 10.1021/acsmedchemlett.2c00127, PMID 36262386; ChEMBL doc `CHEMBL5126642` | same chloroacetamide-class caveat; DEL provenance adds its own false-positive risk (no independent orthogonal confirmation surfaced) |

## Rejected / considered-and-dropped

- **3,5-diaminobenzoic acid Pin1 inhibitor patent** (US patent, "…and therapeutic agent for inflammatory diseases using same") — turned up in search, described as non-covalent. The document I retrieved is an **image-based scanned PDF**; no SMILES or text-searchable structure was extractable. Fails your criterion 3 outright (name/Markush-only is explicitly excluded). Not included. If you have access to a text-searchable version or the compound's CAS number, worth a second pass.
- **Cyclohexyl-ketone Pin1 inhibitors** (Namanja/Etzkorn-lineage work, *PLoS ONE* 2012, DOI 10.1371/journal.pone.0044226) — designed as an electrophilic Cys113 mimetic, but the paper's own conclusion is that **inhibition is non-covalent / substrate-analogue**, not nucleophilic addition. Weak (best IC50 61±8 µM), phosphoserine-peptide scaffold closely related to compounds already in your reference set (Wildemann/Liu-Pei family), and no SMILES was explicitly given (would need derivation from a synthetic scheme, which I did not do given the weak potency). Not included.
- **Chemoproteomics / cysteinome-profiling datasets** (Cravatt/Backus-style isoTOP-ABPP and plate-based reactive-cysteine platforms, e.g. Litwin/Crowley/Suciu/Boger/Cravatt *Tetrahedron Lett* 2021, and the ~38,450-cysteine / 192-electrophile HEK293T dataset) — I specifically went looking here because chemoproteomics papers can surface electrophile classes absent from the med-chem literature, per your prompt. Two relevant facts, both negative for new-chemotype purposes: (1) the Litwin et al. panel found **PIN1_C113 was liganded only by the α-chloroacetamide scout fragment (KB02)**, not by any of their 7 other candidate electrophiles — i.e., in that specific screen Cys113 looked chloroacetamide-selective, not evidence of a new warhead engaging it; (2) the large-scale plate-based/isoTOP-ABPP datasets report reactive-cysteine coverage and a second reactive cysteine (**C57**, lower potency than C113) but do not publish resolvable per-hit compound structures for Pin1 in the material I could access — a real scope limitation (large chemoproteomics resource papers report hit counts and residue IDs, not always structures for every hit protein), noted honestly rather than fabricated.
- **Boronic acid, sulfonyl fluoride/SuFEx, vinyl sulfone, cyanoacrylamide, fumarate, epoxide, aldehyde warheads** — targeted searches for each of these against Pin1/Cys113 came back empty of any primary-literature hit with a resolvable structure. I did not find evidence any of these electrophile classes has been published against Pin1 at all, covalent-fragment screen or otherwise.
- **Tian-2025 pyrimidine-series analogs beyond compound 6a** (already in your table) and **Pu-2025 benzylguanine-series analogs beyond the one row already present** — both documents (`CHEMBL6096711`, `CHEMBL6096731`) contain dozens of additional SAR analogs, all sharing their parent's chemotype. Not independent; not enumerated further here.

## Chemotype-count assessment

**Two genuinely new chemotypes surfaced, not the "2 or more" that would comfortably clear the gate, and one of the two may not even count depending on how your pipeline defines "chemotype."**

1. **Cinnamamide/aryl-acrylamide Michael acceptor (Ieda 2019)** is the clean win: covalent, MS-confirmed Cys113 adduct, kinact/Ki measured, structurally and mechanistically distinct from all 4 existing covalent chemotypes. This also happens to close a gap you already knew about — `warhead_classes_7.csv`'s `acrylamide` row is flagged `"NO validated covalent-Cys113 Pin1 precedent"` (it was T_3's expert-chosen, precedent-free warhead). Ieda's compound isn't plain acrylamide, it's an aryl-conjugated cinnamamide variant, but it is the closest thing to literature precedent for that warhead family that exists. Worth a decision record cross-referencing this back to whatever used that "no precedent" flag.
2. **Cyclopropanecarboxamide/triazole/biphenylnitrile non-covalent series (Liu et al. 2024)** is real, well-validated (SPR + competitive FP + a co-crystal structure, PDB 9INR), and structurally novel relative to everything else non-covalent in your set — but it is **non-covalent**, and your gate's chemotype bookkeeping (`warhead_classes_*.csv`, the R-group decoy-enumeration files) is built entirely around covalent warhead SMARTS. If chemotype-counting for the gate is warhead-keyed, this candidate doesn't move that counter at all, regardless of how good the binding data is. Flagging this explicitly rather than quietly assuming it counts — that determination is the blacksmith's/adversary's call, not mine, since it depends on how the enrichment pipeline itself is wired, not on the literature.

**Honest read on your actual question** — is the Pin1 covalent literature only 4 chemotypes deep: **essentially yes, with one addition.** After a systematic pull of all 396 sub-10-µM ChEMBL2288 actives across 35 source documents (9 of them new to this project) plus targeted chemoproteomics/patent/preprint searches, the covalent side of the literature is chloroacetamide (repeatedly, across at least 4 independent scaffolds: Sulfopin, BJP-06-005-3, ZL-Pin13, the triazine DEL hit), naphthoquinone (Michael acceptor), sulfamate acetamide, BDHI, an SNAr chloroazine that's explicitly a Sulfopin warhead-swap (not independent), and now this one cinnamamide/aryl-acrylamide compound. That's **5** independent covalent chemotypes, not 6. I would not spend further search budget expecting a 6th to appear in the small-molecule med-chem literature — the field has converged hard on chloroacetamide as the default Cys113 electrophile (it's what Dubiella's own 993-fragment screen selected for), and the exhaustive-feeling searches for boronic acid/SuFEx/vinyl sulfone/aldehyde/epoxide/nitrile came back completely empty. If 6 chemotypes is a hard floor, the realistic paths are: (a) accept the non-covalent Liu-2024 series if your gate can be adapted to count it, (b) treat this as the honest finding it is and pursue a different validation strategy per your own framing, or (c) commission covalent-fragment screening data outside the published literature (Dubiella's own 993-fragment dataset reportedly had 111 fragments labeling Pin1 >50%, only a handful of which became named, published chemotypes — the rest may sit in supplementary data I have not obtained).

## Zotero / reading list

`$ZOTERO_USER_ID` / `$ZOTERO_API_KEY` were not set in this environment, so I did not attempt Zotero POSTs. The four primary sources above (Ieda 2019, Liu 2024, the ZL-Pin13 paper, and the triazine-DEL paper) are worth adding to the project reading list and to Zotero (tag `inhibition`) when credentials are available.
