---
id: D0108
title: NO GO on t4_80fbed3bdf1e — its warhead is 8th of 9 on its own scaffold, and two crystallographically validated warheads on the identical molecule reach the crystal band
date: 2026-09-02
status: accepted
approach: t4
decided_by: '@twu383'
origin: user
supersedes: []
superseded_by: null
affects:
  - scripts/warhead_panel_tier1.py
  - decisions/D0072-no-go-on-t4-72f5671e89cb.md
  - decisions/D0107-bpmd-and-non-covalent-residence-both-fail-their-own-positive-control.md
evidence:
  - 'the molecule: O=S1(=O)CC[C@@H](N(c2ccc3c(c2)CCCO3)C2CON=C2Br)C1 -- bdhi_c4 on sulfopin''s sulfolane with a chromane R-group (rg0126). MW 415.3, QED 0.757, SAscore 4.14, LUMO -7.673 eV (in window)'
  - 'TIER-1 |delta d| over 300 ps unrestrained equilibration, SAME R-GROUP, warhead the only variable (n=3 replicates each): snar_chloroazine 0.057, naphthoquinone_c2 0.057, naphthoquinone_benzo 0.082, sulfonate_acetamide 0.121, bdhi_c5 0.173, sulfamate_acetamide 0.214, acrylamide 0.217, bdhi_c4 (THIS MOLECULE) 0.281, chloroacetamide 0.437'
  - 'REF crystallographic median 0.102 nm; three of the nine BEAT it, so the chromane scaffold is not the limitation'
  - 'this molecule sits exactly on the A_hiEnr_hiCons_bdhi group median (0.281) -- the group D0072''s rejected molecule came from'
  - 'BDHI has ZERO crystallographic Cys113 positives: 15 curated structures give chloroacetamide 9, naphthoquinone_c2 4, snar_chloroazine 2. Unchanged since D0072 flagged it, and @twu383 confirmed 2026-09-02 that none exist'
  - 'warhead_status DESIGNED_UNTESTED; the library''s own note: "Places the core adjacent to the reactive carbon, which may hinder Cys approach"'
  - 'engagement rank on the same panel: 8th of 9 (0.725 vs 0.966 best); SAscore the WORST of the nine (4.14 vs 3.14-3.55)'
  - 'in nac_v6 it ranks 125/187 within bdhi_c4 and 1224/1684 across T_4 on best-mode engagement -- below its class median, on no shortlist'
  - 'Welch p vs the two leaders: 0.085 (snar) and 0.098 (naphthoquinone), n=3 each -- suggestive, not significant'
  - 'BPMD and 100 ns residence both fail their positive control (D0107) and contribute nothing to this decision in either direction'
runbook: null
---

# NO GO on `t4_80fbed3bdf1e`

```
O=S1(=O)CC[C@@H](N(c2ccc3c(c2)CCCO3)C2CON=C2Br)C1
C17H19BrN2O4S · MW 415.3 · QED 0.757 · SAscore 4.14 · BDHI-C4 on sulfopin's sulfolane
```

**Do not send this molecule for synthesis.** Two reasons, and the first is new
evidence rather than a restatement of D0072.

## 1. On its own scaffold, its warhead is 8th of 9

T_4 is combinatorial, so nine molecules already existed sharing this molecule's
exact R-group (`rg0126`, chromane) and differing **only** in warhead. All nine
were docked in one screen — one seed, one receptor, one splitter — and measured
on tier 1, the only readout on this project that has passed its own validation
(D0071: separates crystallographic binders from candidates at p = 0.007).

| warhead | xtal Cys113 positives | \|Δd\| nm |
|---|---|---|
| snar_chloroazine | 2 | **0.057** |
| naphthoquinone_c2 | 4 | **0.057** |
| naphthoquinone_benzo | 0 | 0.082 |
| — REF crystallographic median — | | *0.102* |
| sulfonate_acetamide | 0 | 0.121 |
| bdhi_c5 | 0 | 0.173 |
| sulfamate_acetamide | 0 | 0.214 |
| acrylamide | 0 | 0.217 |
| **bdhi_c4 — this molecule** | **0** | **0.281** |
| chloroacetamide | 9 | 0.437 |

**Three of the nine beat the crystallographic median, so the scaffold is not the
problem.** The chromane R-group can hold a warhead as stably as a deposited
binder does. This molecule's warhead does not, and it lands on the same 0.281
group median that D0072's rejected molecule came from.

## 2. Its warhead class still has no validation of any kind

D0072 flagged that BDHI has zero crystallographic Cys113 positives. Four weeks
on, the curated count is unchanged — chloroacetamide 9, naphthoquinone_c2 4,
snar_chloroazine 2, BDHI 0 — and @twu383 confirmed directly that no BDHI
positive exists to be found. `warhead_status` is `DESIGNED_UNTESTED` and the
library's own note anticipates the failure mode: *"places the core adjacent to
the reactive carbon, which may hinder Cys approach."*

Committing synthesis to an unvalidated warhead is defensible when the molecule
is otherwise the strongest thing available. It is 8th of 9 on its own scaffold
and has the worst SAscore of the nine.

## What is explicitly NOT the reason

- **Med chem does not object.** MW 415, QED 0.757, cLogP 2.11, `synth_fail`
  False, and LUMO −7.673 eV sits *inside* the reactivity safety window — which
  the two recommended alternatives do not.
- **Boltz-2 supports the site.** Independent co-folding places it 3.62 Å from
  Cys113 SG, at the catalytic cysteine and not Cys57, confidence 0.973. That is
  real orthogonal evidence for *where* it binds.
- **Its 100 ns trajectory is fine** — 80% engaged, mean RMSD 0.314 nm, better
  than sulfopin's. That readout discriminates nothing (D0107) and is not being
  used against it either.
- **BPMD says nothing.** 3/3 escaped, and so does the crystal pose (D0107).
- **The ranking is not the reason.** `rank_validated = False` throughout.

## The replacement

**`t4_b720c5a33d32`** — same chromane, same sulfolane, naphthoquinone_c2 warhead:

```
O=C1C=C(N(c2ccc3c(c2)CCCO3)[C@@H]2CCS(=O)(=O)C2)C(=O)c2ccccc21
MW 423.5 · QED 0.755 · SAscore 3.39 · |delta d| 0.057 nm · 4 crystallographic Cys113 positives
```

Tier 1 0.057 against 0.281, QED unchanged (0.755 vs 0.757), **easier to make**
(SAscore 3.39 vs 4.14), and its warhead class has four deposited structures
reacting at this exact cysteine. `t4_03974eb61946` (SNAr, 2 positives) ties it on
tier 1 but QED falls to 0.418.

**The trade, stated rather than buried:** both alternatives have LUMO far below
the reactivity window — naphthoquinone −9.18 eV, SNAr −9.56 eV, against a floor
of −7.82. They are predicted substantially more electrophilic than the
wet-lab-validated selective anchors that bound that window, which is an
off-target-reactivity risk this molecule does not carry. That is a chemist's
judgement, not a computational one, and it is the reason this record recommends
rather than decides.

## Limits of this decision

- **n = 3 replicates.** Welch p = 0.085 and 0.098 against the two leaders.
  Suggestive, not significant. Six more replicates would settle it in ~30 min.
- **One molecule per warhead class**, so "this warhead class" and "this
  particular molecule" cannot be separated. The result is about nine specific
  compounds sharing an R-group.
- **Chloroacetamide is worst (0.437) despite nine crystal positives.** Class
  validation does not straightforwardly predict tier-1 on this scaffold, which
  is the strongest single warning against over-reading the table — including
  against reading the two leaders as safe.
- **Tier 1 measures whether the docked pose survives equilibration.** It does
  not measure reactivity, selectivity, or whether anything binds.

## What would reverse it

Any BDHI positive against Cys113 — one structure, or one measured IC50 in the
series. That was D0072's condition, it is still the condition, and it is now
known not to exist.
