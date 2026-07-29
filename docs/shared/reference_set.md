# The frozen reference set

!!! info "Generated page"
    Rendered at build time from the repo's source of truth. Edit the underlying file, not this page.

The single source for two things and nothing else: the **novelty axis** for every approach (`1 - max Tanimoto ECFP4`, computed against this set and **never against the seed**), and **T_4's reactivity window**.

## Master set — novelty axis

`data/reference/pin1_reference_binders_1.csv`

Rows marked `UNVERIFIED` are excluded from the novelty computation.

| name | canonical_smiles | mechanism | warhead_class | potency | promiscuity_flag | pdb |
|---|---|---|---|---|---|---|
| `Sulfopin` | `CC(C)(C)CN(C1CCS(=O)(=O)C1)C(=O)CCl` | `covalent_cys113` | `chloroacetamide` | `nanomolar (covalent)` | `n` | `6VAJ` |
| `KPT-6566` | `CC(C)(C)C1=CC=C(C=C1)S(=O)(=O)N=C2C=C(C(=O)C…` | `covalent_cys113` | `aryl-sulfonyl-acetate (self-immolative -> 1;…` | `IC50 ~0.3-1.4 uM` | `y` |  |
| `Juglone` | `C1=CC2=C(C(=O)C=CC2=O)C(=C1)O` | `covalent_cys113` | `1;4-naphthoquinone (Michael acceptor)` | `kinact 5.3e-4-4.5e-3 s^-1` | `y` |  |
| `BJP-06-005-3` | `CCOC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H](Cc1c[nH…` | `covalent_cys113` | `chloroacetamide (N-methyl peptidomimetic)` | `potent/selective` | `n` |  |
| `Reddi-sulfamate-acetamide` | `UNVERIFIED` | `covalent_cys113` | `sulfamate acetamide` | `potent/selective` | `n` |  |
| `Byun-BDHI-fragment` | `UNVERIFIED` | `covalent_cys113` | `3-bromo-4;5-dihydroisoxazole (BDHI)` | `fragment-level` | `n` |  |
| `ATRA` | `CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(=C/C(=O)O)…` | `non_covalent` |  | `low-uM (binds + degrades Pin1)` | `y` |  |
| `EGCG` | `C1[C@H]([C@H](OC2=CC(=CC(=C21)O)O)C3=CC(=C(C…` | `non_covalent` |  | `direct binding; PPIase inhibition` | `y` |  |
| `PiB` | `CCOC(=O)CN1C(=O)c2ccc3c4c(ccc(c24)C1=O)C(=O)…` | `non_covalent` |  | `low-uM PPIase inhibition` | `y` |  |
| `reversible-Thr-sulfonamide-quinazoline` | `C[C@@H](O)[C@H](N)C(=O)NS(=O)(=O)c1cccc(-c2c…` | `non_covalent` |  | `Ki 1.2 nM (non-drug-like)` | `n` |  |
| `reversible-Thr-sulfamate-adenosine` | `C[C@@H](O)[C@H](N)C(=O)NS(=O)(=O)OC[C@H]1O[C…` | `non_covalent` |  | `IC50 15 nM (non-drug-like)` | `n` |  |

## Covalent Cys113 anchors — reactivity window

`data/reference/pin1_covalent_cys113_anchors_2.csv`

`reference_set.py` refuses `UNVERIFIED` rows into the window.

| anchor_id | name | warhead_class | canonical_smiles | smiles_status | potency_kinetics | promiscuity_flag |
|---|---|---|---|---|---|---|
| `1` | `Sulfopin` | `chloroacetamide` | `CC(C)(C)CN(C1CCS(=O)(=O)C1)C(=O)CCl` | `verified` | `nanomolar (covalent); k=0.028 M-1 s-1` | `n` |
| `2` | `BJP-06-005-3` | `chloroacetamide (N-methyl; peptidomimetic)` | `CCOC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H](Cc1c[nH…` | `verified` | `potent/selective` | `n` |
| `3` | `KPT-6566` | `aryl-sulfonyl-acetate -> 1;4-naphthoquinone` | `CC(C)(C)C1=CC=C(C=C1)S(=O)(=O)N=C2C=C(C(=O)C…` | `verified` | `IC50 ~0.3-1.4 uM` | `y` |
| `4` | `Juglone` | `1;4-naphthoquinone Michael acceptor` | `C1=CC2=C(C(=O)C=CC2=O)C(=C1)O` | `verified` | `kinact 5e-4-4.5e-3 s^-1` | `y` |
| `5` | `Reddi-2023-4d` | `sulfamate acetamide` | `CC(C)(C)CN(C1CCS(=O)(=O)C1)C(=O)COS(=O)(=O)N…` | `verified` | `k=0.023 M-1 s-1; 95% labeling; chemoproteomi…` | `n` |
| `6` | `Reddi-2023-4g` | `sulfamate acetamide` | `C1CCCCC1CN(C2CCS(=O)(=O)C2)C(=O)COS(=O)(=O)N…` | `verified` | `k=0.068 M-1 s-1; 93% labeling; most potent i…` | `n` |
| `7` | `Byun-2023-BDHI-fragment` | `3-bromo-4;5-dihydroisoxazole (BDHI)` | `UNVERIFIED` | `UNVERIFIED` | `fragment-level` | `n` |

## Warhead classes — T_4 enumeration

`data/reference/warhead_classes_2.csv`

`warhead_library.enumerable()` defaults to `VERIFIED` only.

| class_id | display_name | warhead_fragment_smiles | mechanism | reactive_atom_smarts | leaving_group | precedent | structure_status |
|---|---|---|---|---|---|---|---|
| `chloroacetamide` | `Chloroacetamide` | `[*]C(=O)CCl` | `sn2_displacement` | `[CH2][Cl]` | `chloride` | `anchored_verified` | `VERIFIED` |
| `sulfamate_acetamide` | `Sulfamate acetamide` | `[*]C(=O)COS(=O)(=O)N[*]` | `sn2_displacement` | `[CH2][OX2][SX4](=O)(=O)` | `sulfamate` | `anchored_verified` | `VERIFIED` |
| `sulfonate_acetamide` | `Sulfonate acetamide (mesylate)` | `[*]C(=O)COS(=O)(=O)C` | `sn2_displacement` | `[CH2][OX2][SX4](=O)(=O)` | `sulfonate` | `anchored_verified` | `VERIFIED` |
| `bdhi` | `3-bromo-4,5-dihydroisoxazole (BDHI)` | `[*]C1=NOCC1` | `sn2_ring_opening` | `[CX3]([Br])=[NX2]` | `bromide` | `anchored_verified` | `VERIFIED_CLASS_ONLY` |
| `naphthoquinone` | `1,4-naphthoquinone Michael acceptor` | `UNRESOLVED` | `michael_addition` | `[CX3]=[CX3][CX3]=O` | `none_addition` | `anchored_verified` | `NEEDS_DESIGN` |

## Measured reactivity kinetics

`data/reference/pin1_reactivity_kinetics_1.csv`

Digitized from a figure to ~1 significant figure. Bound a window with these; do not treat them as precise.

| compound | warhead_class | r_substituent | canonical_smiles | k_M-1_s-1 | pct_labeling | value_status |
|---|---|---|---|---|---|---|
| `Sulfopin` | `chloroacetamide` |  | `CC(C)(C)CN(C1CCS(=O)(=O)C1)C(=O)CCl` | `0.028` | `55` | `DIGITIZED_FROM_FIGURE` |
| `4a` | `sulfonate_acetamide` | `Me (sulfonate)` | `CC(C)(C)CN(C1CCS(=O)(=O)C1)C(=O)COS(=O)(=O)C` | `0.030` | `17` | `DIGITIZED_FROM_FIGURE` |
| `4b` | `sulfamate_acetamide` | `Me` | `CC(C)(C)CN(C1CCS(=O)(=O)C1)C(=O)COS(=O)(=O)NC` | `0.005` | `8` | `DIGITIZED_FROM_FIGURE` |
| `4c` | `sulfamate_acetamide` | `Bn` | `CC(C)(C)CN(C1CCS(=O)(=O)C1)C(=O)COS(=O)(=O)N…` | `0.022` | `57` | `DIGITIZED_FROM_FIGURE` |
| `4d` | `sulfamate_acetamide` | `Ph` | `CC(C)(C)CN(C1CCS(=O)(=O)C1)C(=O)COS(=O)(=O)N…` | `0.023` | `95` | `DIGITIZED_FROM_FIGURE` |
| `4e` | `sulfamate_acetamide` | `4-Br-Ph` | `CC(C)(C)CN(C1CCS(=O)(=O)C1)C(=O)COS(=O)(=O)N…` | `0.007` | `97` | `DIGITIZED_FROM_FIGURE` |
| `4f` | `sulfamate_acetamide` | `4-Me-Ph` | `CC(C)(C)CN(C1CCS(=O)(=O)C1)C(=O)COS(=O)(=O)N…` | `0.006` | `28` | `DIGITIZED_FROM_FIGURE` |
| `4g` | `sulfamate_acetamide` | `Ph (N-cyclohexylmethyl core)` | `C1CCCCC1CN(C2CCS(=O)(=O)C2)C(=O)COS(=O)(=O)N…` | `0.068` | `93` | `DIGITIZED_FROM_FIGURE` |

---

## Provenance


Source: issue #108, Rev 3 spec (comment 5083543621), Appendix A. Assembled by the
**bookworm**. Structures from PubChem/ChEMBL; literature confirmations from PubMed.
**No SMILES or DOIs were invented.** Structures that could not be verified against a
public record carry `UNVERIFIED` rather than a guess.

This set is the single source for two things, and nothing else:

1. the **novelty axis** for every approach — `1 - max Tanimoto (ECFP4)` against
   `pin1_reference_binders_1.csv`, **never against the seed** (closes adversary
   finding B4, novelty-axis circularity);
2. **T_4's covalent reactivity window** — anchored on
   `pin1_covalent_cys113_anchors_1.csv` (closes B5, n≈1 anchor).

The project's own computational leads are **excluded as anchors** by design.

## Files

| File | Rows | Feeds |
|---|---|---|
| `pin1_reference_binders_1.csv` | 11 validated binders (9 with verified SMILES) | novelty axis, all approaches |
| `pin1_covalent_cys113_anchors_1.csv` | 6 covalent-at-Cys113 (4 verified) | T_4 reactivity window |

Integer-versioned per the lab's data rule: `_1` is the first freeze. A revised set
becomes `_2`; the old file is retired via `src/ready_to_delete.md`, never edited in
place, so a pipeline run always pins an exact version.

## CSV escaping note

Two source values contain commas inside chemical names (`1,4-naphthoquinone`,
`3-bromo-4,5-dihydroisoxazole`). These are written with semicolons
(`1;4-naphthoquinone`) in the `warhead_class` column so the CSV parses unambiguously
regardless of reader. The chemistry is unchanged; only the delimiter-safe rendering
of the class label differs from Appendix A's prose.

## Open item — blocks T_4's reactivity window only

Two covalent anchors have **SI-only structures** and are marked `UNVERIFIED`:

- **Reddi 2023** sulfamate-acetamide — JACS 10.1021/jacs.2c08853
- **Byun 2023** BDHI fragment — JACS 10.1021/jacs.3c00598

`shared/reference_set.py` **refuses to feed an `UNVERIFIED` row into the reactivity
window**. Until the bookworm pulls their SMILES from the paper SI, T_4's LUMO window
is computed from the 4 verified anchors alone. Nothing else is blocked: the novelty
axis uses the 9 verified master-set structures and is unaffected.

## Honest caveat carried forward (Rev 3, B5 disposition)

By headcount there are 6 covalent anchors; by distinct electrophile chemistry roughly
4. But the *clean, selective* anchors (Sulfopin, BJP-06-005-3) are **both
chloroacetamides**, and the extra kinetic anchors (juglone, KPT-6566) are
**promiscuous quinones** (`promiscuity_flag = y`). The resulting LUMO window is
therefore **chloroacetamide-centric with a reactive-quinone upper tail — not a
chemotype-balanced distribution.** Do not over-trust a single global cutoff.

That scarcity of clean, drug-like, selective covalent-Cys113 chemistry is itself a
finding the choreography carries forward, not a defect in the set.

## Excluded, and why (no invented structures)

| Excluded | Reason |
|---|---|
| API-1 | name collides with an Akt inhibitor; no confirmed Pin1 binding |
| buparvaquone, DTM | no primary Pin1 citation found on this pass |
| borrelidin | not a bona fide Pin1 ligand |
| phosphopeptide / WW-domain peptides | genuine binders, but non-drug-like with no useful SMILES |

## Seeds — why only ATRA and sulfopin

The reference set confirms the seed choice. The other validated binders are either
promiscuous aggregators (EGCG, PiB, juglone, KPT-6566) or non-drug-like
peptidomimetics (BJP-06-005-3, the reversible active-site peptides) — poor starting
points. So candidates are **not** fanned across all binders as seeds. The seed is an
explicit parameter in `config/seeds.yaml`: running an approach from an additional
curated, drug-like, non-promiscuous seed is a config change, not a redesign. Strategy
diversity comes from the four approaches, not from multiplying seeds.

---

## 2026-07-27 update — Reddi 2023 sulfamate acetamide RESOLVED

Figure 5 of Reddi 2023 (JACS 10.1021/jacs.2c08853) was supplied directly and the
`UNVERIFIED` sulfamate-acetamide anchor is now resolved. Panel A gives the
structures of the sulfopin analogs **4a–4g**; panels B–E give adduct mass,
kinetics, cellular potency, and chemoproteomic selectivity.

**How the transcription was checked.** These structures were read from a figure
image, not pulled from a database, so they were cross-checked two ways before
being accepted:

1. **Core check.** The transcribed Sulfopin core returns C11H20ClNO3S / 281.80,
   the known Sulfopin formula — so the neopentyl-N-(sulfolan-3-yl) scaffold was
   read correctly.
2. **Adduct-mass check.** Panel B reports a **+272 Da** shift for Pin1 + **4g**.
   For SN2 displacement (Cys-SH + R-CH2-OSO2NHPh → Cys-S-CH2-R + HOSO2NHPh) the
   predicted shift from the transcribed 4g is **271.4 Da**. Match.

Both checks are recorded because a mis-read structure here would silently
corrupt the reactivity window, which is a control.

**What changed:**

| File | Change |
|---|---|
| `pin1_covalent_cys113_anchors_2.csv` | Reddi row (`UNVERIFIED`) → two verified anchors, **4d** and **4g**. Anchors 6 → 7; verified 4 → 6. |
| `warhead_classes_2.csv` | `sulfamate_acetamide` `UNVERIFIED` → `VERIFIED`; `sulfonate_acetamide` (4a) added. Enumerable classes 1 → 3. |
| `pin1_reactivity_kinetics_1.csv` | new — measured k and % labeling for Sulfopin and 4a–4g. |

**The kinetics values are DIGITIZED FROM A SCATTER PLOT (Figure 5C), by eye.**
They are good to roughly one significant figure and are flagged
`DIGITIZED_FROM_FIGURE` in the CSV. Use them to bound a window, not as precise
measurements. If exact values matter, pull them from the paper's SI tables.

**A finding worth carrying forward.** Across Sulfopin and 4a–4g, intrinsic
reactivity and Pin1 labeling efficiency are only weakly related — Pearson
**r = 0.396** over 8 compounds, spanning a 13.6× range of k:

- **4e** has nearly the lowest k (0.007) and the highest labeling (97%).
- **4a** has among the highest k (0.030) and the lowest labeling (17%).

So within the precedented range, intrinsic electrophilicity does **not** predict
target engagement — recognition does. That means the reactivity window is a
**safety filter for condition (ii)**, not a potency predictor, and T_4 must not
rank candidates by LUMO. It also means the window can be anchored on measured
kinetics (0.005–0.068 M⁻¹ s⁻¹ across real actives) rather than on computed LUMO
alone, which is a strictly better anchor.

**Still unresolved:** the Byun 2023 BDHI *fragment* (anchor 7). The BDHI warhead
*class* is verified from PubChem CID 21983498, which is what the window needs,
but the specific fragment structure is still SI-only.
