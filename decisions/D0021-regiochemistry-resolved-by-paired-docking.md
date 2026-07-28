---
id: D0021
title: BDHI and naphthoquinone attachment points resolved by paired docking
date: 2026-07-27
status: superseded
approach: t4
decided_by: '@mhallet'
origin: implementation
supersedes: []
superseded_by: D0024
affects:
  - approaches/t4_combinatorial/05_regiochemistry_comparison.py
  - config/approaches/t4_combinatorial.yaml
  - data/reference/warhead_classes_3.csv
  - integration/app/DECISIONS_TAB_SPEC.md
evidence:
  - 'naphthoquinone: benzo poses in 69 pairs where c2 does not, vs 5 the other way, McNemar p = 1.8e-15'
  - 'naphthoquinone_c2 finds NO favourable pose for 96.8% of its 187 R-groups; benzo 62.6%'
  - 'naphthoquinone paired median difference +4.68 kcal/mol favouring benzo, Wilcoxon p = 2.5e-4'
  - 'bdhi: c5 poses in 59 pairs where c4 does not, vs 34 the other way, McNemar p = 0.0124'
  - 'bdhi no-pose 50.3% (c4) vs 36.9% (c5); paired median difference +1.09 kcal/mol favouring c5, Wilcoxon p = 6.6e-3'
  - 'both endpoints — pose success and affinity — favour the same arm in both pairs'
  - 'all four classes are paired on the SAME 187 R-groups, so the comparison is matched'
runbook: null
---

## RESOLUTION AFTER RE-DOCK (2026-07-27, same day)

The re-dock on adduct-form ligands is done, and **the BDHI call reverses**.

| | old (pre-reaction) | new (adduct) |
|---|---|---|
| `bdhi_c4` median | +0.01 | **−3.79** |
| `bdhi_c5` median | −2.22 | −2.87 |
| `bdhi_c4` no-pose | 50% | **14%** |
| `bdhi_c5` no-pose | 37% | 20% |

**Carry `bdhi_c4`, drop `bdhi_c5`** — the opposite of what this record originally
concluded. Verdict **UNDERPOWERED**: the affinity difference is real
(median −0.58 kcal/mol, p = 0.0035) but modest, and the pose-success endpoint no
longer separates the arms at all (18 vs 29 discordant pairs, p = 0.14). With the
clash artifact gone, most poses now succeed in both arms, so the binary endpoint
has little left to discriminate. Treat as a weak preference, not a finding.

The original call was an artifact, and the artifact was large enough to invert
the answer. Note it did so through the *loser*: `bdhi_c4`'s retained bromine
blocked half its poses, which read as a geometric failure of the C4 attachment
when it was really a failure to remove a leaving group.

**Naphthoquinone is unchanged and remains STRONG** — benzo over C2, 69 discordant
pairs against 6, p = 1.2e-14. This is what the withdrawal notice predicted:
Michael acceptors carry no leaving group, so nothing about that comparison could
move. The prediction holding is mild evidence the diagnosis was right.

## WITHDRAWAL NOTICE (2026-07-27, same day)

**The BDHI half of this decision is withdrawn. The naphthoquinone half stands.**

Preparing the MM-GBSA input revealed that gnina was docking the **pre-reaction**
ligand: the leaving group was never removed, so gnina formed the S–C bond onto a
carbon whose valence was already full. Every docked reactive carbon has five
bonds. See D0022.

For a leaving group on a flexible sp3 carbon this is survivable — chloroacetamide
puts its chlorine a median 2.93 Å from SG and **0%** of its poses clash. For a
leaving group on a rigid sp2/aromatic carbon it is not: the ring geometry fixes
the C–halogen direction, so bonding S to that carbon drives the halogen into the
sulfur.

- `bdhi_c5`: halogen a median **1.63 Å** from SG, **68%** of poses clashing
- `bdhi_c4`: **47%** clashing
- `snar_chloroazine`: **56%** clashing, with contacts as short as 0.89 Å

**Both BDHI arms are contaminated, and the arm that won has the worse clash
rate.** The pose-success endpoint that decided it cannot be distinguished from a
clash artifact. The BDHI call is withdrawn pending a re-dock on adduct-form
ligands.

**The naphthoquinone call stands.** Both arms are Michael acceptors with no
leaving group at all, so neither is affected: 0% clash for `naphthoquinone_benzo`
and 3.2% for `naphthoquinone_c2`. The C2-versus-benzo comparison was between two
ligands of identical composition, and the 97%-no-pose result for C2 reflects the
geometry it was measuring.

## Context

Two chemotypes entered T_4 with a verified warhead class but an *untested*
attachment point, and the PI's instruction was explicit: "For BDHI and
Naphthoquinone, can we just try all the four cases? I just don't know."

So all four were enumerated as separate classes, and
`config/approaches/t4_combinatorial.yaml` declared in advance what would settle
it: `discriminated_by: [warhead_validity_gate_5b, covalent_docking_geometry,
lumo_window]`. All four passed 5b — coupling does not destroy the electrophile
in any of them — and the LUMO window did not separate the members of either
pair. That left the docking geometry, which is the physically apt test anyway:
the question is whether Cys113 can reach the reactive atom with the core in the
way.

Each regiochemistry was built against the same 187 R-groups, so the arms are
matched pair-for-pair on the only other varying factor. An R-group that docks
well does so in both arms and cancels.

## Decision

**Carry `naphthoquinone_benzo` and `bdhi_c5`. Stamp `naphthoquinone_c2` and
`bdhi_c4` as superseded regiochemistries.** Neither is deleted — both keep their
rows and their docking evidence, per stamp-don't-delete.

Two endpoints were measured on each pair:

- **Pose success** (paired binary, McNemar). A Vina-style score at or above zero
  is not a weak binding energy; it means the search found no favourable pose.
- **Affinity** (Wilcoxon signed-rank, with matched-pairs rank-biserial effect
  size computed from the signed rank sums, not from win counts).

For **naphthoquinone** the result is not close. Attachment at C2 fails to pose
for 97% of R-groups; in 69 pairs the benzo isomer poses where C2 does not,
against 5 the other way. This matches the chemistry: C2 sits adjacent to the
Michael-acceptor positions, so the core is placed directly in the path of the
approaching cysteine. Verdict **STRONG**.

For **BDHI** the same direction holds with a modest margin: 59 discordant pairs
favouring C5 against 34 favouring C4, p = 0.012. Verdict **WEAK**. C4 places the
core adjacent to the reactive carbon, which the config anticipated might hinder
Cys approach; it does, but not decisively.

### Where the primary endpoint came from, and a knife edge

When most of the affinity column is censored, a rank test on it compares noise.
The verdict function therefore treats pose success as primary above 50% no-pose
in either arm. **That branch was written after seeing that C2 fails 97% of the
time** — the exact circumstance in which a threshold change can become a way of
choosing the answer, so it is recorded rather than buried. It does not change
any conclusion: both endpoints favour the same arm in both pairs, so only the
confidence label moves, never the winner.

BDHI lands on the knife edge — 50.3% no-pose, against a 50% threshold. Under the
affinity branch it would grade **UNDERPOWERED** instead of **WEAK**. The winner
is C5 either way; only the strength of the claim is sensitive to a threshold
that a handful of R-groups could flip. Treat the BDHI call as provisional.

## Consequences

Step 9 (MM-GBSA on the true covalent adduct) is the most expensive stage in the
plan, and this halves the quinone and BDHI work it has to cover.

The reporting rule in `t4_combinatorial.yaml` — `keep_regiochemistries_separate:
true` — did its job and stays. Merging the arms would have averaged a working
geometry with a non-working one and reported a mediocre chemotype where there
are really one good and one bad attachment point.

The losing arms remain available: if MM-GBSA later contradicts the docking
geometry for the surviving arm, the alternative has not been thrown away.

None of this is a measured potency difference. It is a decision about which
attachment to carry forward, on docking geometry, which is the evidence the
config committed to in advance.
