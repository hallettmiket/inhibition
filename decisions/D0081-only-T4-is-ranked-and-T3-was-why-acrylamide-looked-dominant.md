---
id: D0081
title: Only T_4 is ranked — and T_3 is why acrylamide looked like the library
date: 2026-08-12
status: accepted
approach: shared
decided_by: '@tt8804'
origin: user
supersedes: []
superseded_by: null
affects:
  - config/target.yaml
  - scripts/rank_v2.py
  - shared/mode_ranking.py
  - scripts/sweep_gap_worklist.py
evidence:
  - '@tt8804, 2026-08-12: "t_3???????????? we are only doing t_4" and "no wonder there are so many acrylamides"'
  - 'T_3 is 96% acrylamide: 4,062 of the 4,249 acrylamide molecules, against 187 from T_4'
  - 'every other warhead class is 100% T_4, at 162-188 molecules each'
  - 'T_3 contributes 19,888 of 34,059 ranked modes -- 58% of the ranking, one warhead class, one generator'
  - 'the acrylamide sweep depth of 85 modes was 84 T_3 rows and 1 T_4 row'
  - '83 of the 259 modes on the worklist were T_3, all acrylamide; 4 were swept (~84 min GPU) before this was caught'
  - 'the tier filter is applied BEFORE ranking, so class_rank counts only molecules in contention'
  - 'shared/mode_ranking reads run.tiers as well: dropping a tier stops a NEW file being written but does not remove the old one, and the reader takes the newest match per tier'
---

# D0081 — only T_4 is ranked, and T_3 is why acrylamide looked dominant

## The decision

`run.tiers: ["T4"]` in `config/target.yaml`. T_3 rows stay screened and stay on
disk; they are no longer **ranked**, which is what the GUI shows and what the
sweep can select from.

## What T_3 was doing to every class-composition number

T_3 is the REINVENT de novo output and it is essentially a single-warhead
library:

| class | T_3 | T_4 | % T_3 |
|---|---|---|---|
| acrylamide | 4,062 | 187 | **96%** |
| chloroacetamide | 0 | 188 | 0% |
| naphthoquinone_benzo | 0 | 187 | 0% |
| snar_chloroazine | 0 | 187 | 0% |
| naphthoquinone_c2 | 0 | 187 | 0% |
| sulfamate_acetamide | 0 | 187 | 0% |
| sulfonate_acetamide | 0 | 187 | 0% |
| bdhi_c4 | 0 | 162 | 0% |
| bdhi_c5 | 0 | 162 | 0% |

So "acrylamide is 76% of the in-scope modes" was never a statement about
chemistry. It was a statement about which generator ran, and it silently
propagated into every downstream argument:

- the **per-family quota** was justified as stopping acrylamide from crowding
  out the smaller families — the crowding was T_3;
- the **depth** for acrylamide came out at 85 modes clearing enrichment 4.0, of
  which **84 were T_3 and one was T_4**;
- the **Chalcopyrite-vs-Galena class shifts** were computed over a pool 58% of
  which is one generator's acrylamides.

The R-group profile is the exception and survives: it was computed on T_4 alone
(the top 100 was already 96% T_4), so the N-heteroaromatic enrichment of +11.5
points still stands.

## Why the filter goes before the ranking, not after

`class_rank` is a position within a warhead class. Ranking both tiers together
and dropping T_3 afterwards would leave acrylamide's T_4 rows carrying ranks
earned against 4,062 REINVENT molecules nobody intends to synthesise — a rank of
"90th" that is really "second, behind 88 things we are not making".

## Why the GUI reader needs the same key

Removing a tier from `rank_v2` stops it writing a **new** file for that tier. It
does not remove the old one, and `mode_ranking.gather` resolves the newest match
for each tier independently — so without the same filter the view would keep
showing a pre-decision T_3 table beside T_4 rows ranked without it. Same shape as
D0080: one decision, two readers, and the one that was not told carries on.

## Cost

83 of the 259 modes on the sweep worklist were T_3 acrylamide. Four were swept
before this was caught — about 84 GPU-minutes. The sweep was stopped, the
worklist rebuilt from the T_4-only ranking, and the T_3 results remain on disk
and are valid measurements of the molecules they were run on; they simply belong
to a tier that is no longer in contention.
