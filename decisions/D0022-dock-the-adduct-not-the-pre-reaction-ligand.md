---
id: D0022
title: Covalent docking must use the adduct form, not the pre-reaction ligand
date: 2026-07-27
status: accepted
approach: t4
decided_by: '@mhallet'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - shared/covalent_protocol.py
  - shared/warhead_library.py
  - data/reference/warhead_classes_3.csv
  - approaches/t4_combinatorial/03_covalent_dock.py
  - approaches/t4_combinatorial/05_regiochemistry_comparison.py
  - decisions/D0021-regiochemistry-resolved-by-paired-docking.md
evidence:
  - 'every docked reactive carbon sits 1.81 A from Cys113 SG while retaining its leaving group — 5 bonds on carbon'
  - 'chloroacetamide pose t4_5e235921c8c0: reactive C has 2 heavy neighbours + 2 H + the new S bond'
  - 'clash rates (leaving-group atom within 2.5 A of SG): bdhi_c5 67.9%, snar_chloroazine 55.6%, bdhi_c4 47.1%, chloroacetamide 0%'
  - 'shortest observed contact 0.89 A (Cl to SG, snar_chloroazine top-ranked pose)'
  - 'chloroacetamide, sulfamate_acetamide and sulfonate_acetamide yield ONE identical adduct — their docking spread is leaving-group artifact'
  - 'sulfamate carried a 12-heavy-atom leaving group through docking and ranked best by median (-5.47)'
  - 'the Michael acceptors have no leaving group and are unaffected (0-3.2% clash)'
runbook: null
---

## Context

gnina's `--covalent_lig_atom_pattern` does exactly what its help text says: it
picks a ligand atom and bonds it to the receptor atom. It does not perform
reaction chemistry, and it does not remove a leaving group. The gnina
documentation says nothing either way about how the ligand should be prepared —
so the convention had to be established from the output, and the output is
unambiguous.

We supplied the intact, pre-reaction ligand. The reactive carbon therefore
already had a full valence, and gnina added the S–C bond on top of it. **Every
one of the 1,683 docked complexes contains a pentavalent carbon.**

For a leaving group on a flexible sp3 carbon this is survivable. Chloroacetamide
puts its chlorine a median 2.93 Å from SG and not one of its 187 poses clashes —
the CH2 rotates and the chlorine finds somewhere to go.

For a leaving group on a rigid sp2 or aromatic carbon it is not survivable. The
ring fixes the C–halogen direction, so bonding S to that carbon drives the
halogen into the sulfur:

| class | median halogen···SG | poses clashing (<2.5 Å) |
|---|---|---|
| `bdhi_c5` | 1.63 Å | 67.9% |
| `snar_chloroazine` | 2.38 Å | 55.6% |
| `bdhi_c4` | 2.55 Å | 47.1% |
| `chloroacetamide` | 2.93 Å | 0% |

The shortest contact observed is **0.89 Å** — two atoms occupying the same
space — and it belongs to the pose that scored best in the entire run
(`snar_chloroazine`, −9.16 kcal/mol).

### The three acetamides are the same molecule

Worse than the clashes, and easy to miss: `chloroacetamide`,
`sulfamate_acetamide` and `sulfonate_acetamide` are all SN2 displacements at the
same CH2. They differ only in what leaves. **Their adducts are one identical
molecule.** Verified directly: applying each displacement to the same R-group
gives a single distinct product SMILES.

So any difference between those three classes in the docking results is
*entirely* leaving-group artifact — and the differences were large, spanning
−5.47 to −2.87 kcal/mol in median. Sulfamate topped the table while carrying a
12-heavy-atom phantom group that is not present in anything that binds Pin1.
What actually distinguishes those three warheads is kinetics, which is what the
reactivity window measures, and which docking cannot see.

## Decision

**Dock the adduct form: the post-reaction ligand with the leaving group removed
and the attachment atom left with an open valence for gnina to fill.**

This requires a second SMARTS per class. `reactive_atom_smarts` identifies the
reactive atom *in the pre-reaction molecule* and names the leaving group in
doing so (`[CH2][Cl]`, `[c]([Cl])[n]`); those patterns cannot match the adduct,
because the atom they key on is gone. Confirmed empirically: re-docking a
leaving-group-stripped ligand under the existing SMARTS produces no reactive
atom and no result. The library therefore needs a per-class adduct transform and
an `adduct_attachment_smarts` that matches the product.

Both forms stay in the library. The pre-reaction SMARTS is still what the 5b
validity gate and the warhead tests need — the question "is this a genuine
electrophile of its class" is a question about the *unreacted* warhead.

## Consequences

**T_4's docking must be re-run.** The current per-class table is not a valid
cross-class comparison, and D0015's choice of rank metric is unaffected but the
values it ranks are not.

**D0021 is partially withdrawn.** The BDHI regiochemistry call rested on pose
success, and both BDHI arms clash heavily — with the winning arm clashing *more*
(67.9% vs 47.1%). That call cannot be separated from the artifact. The
naphthoquinone call stands: both arms are Michael acceptors carrying no leaving
group, so their comparison was between two ligands of identical composition.

**The enrichment gate (D0015, D0016) needs re-checking.** It was measured on
6 actives and 294 warhead-bearing decoys through this same protocol. Actives and
decoys were treated identically, so the comparison is internally consistent and
ROC-AUC is unlikely to move much — but "unlikely" is a prediction, and it should
be re-measured rather than assumed.

**After the fix, the three acetamide classes should converge.** Their adducts
are identical, so their docking scores must agree to within the search's own
noise. That gives a free internal control: if they do not converge after
re-docking, something in the protocol is still wrong. Their reactivity
differences remain real and remain the reactivity window's business.

Michael acceptors gain an H on the α-carbon in the true adduct, which docking
has been ignoring. That is a one-atom difference on a flexible position and is
minor next to a 12-atom sulfamate, but the adduct transform should handle it for
consistency rather than special-casing.
