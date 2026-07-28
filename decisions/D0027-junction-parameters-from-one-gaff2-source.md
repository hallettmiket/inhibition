---
id: D0027
title: Derive every Cys-ligand junction parameter from one GAFF2 source
date: 2026-07-28
status: accepted
approach: t4
decided_by: '@mhallet'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - data/params/cys_gaff2_junction_2.frcmod
  - shared/mmgbsa.py
  - approaches/t4_combinatorial/06_mmgbsa.py
evidence:
  - '18 of 27 MM-GBSA builds failed; the 9 that succeeded were exactly the three SN2 acetamide classes'
  - 'missing terms: S-c2 (9 candidates), S-cc (6), S-ca (3), plus 9 distinct angle types'
  - 'junction v1 covered only sp3 carbon (c3); 6 of 9 warhead classes attach through sp2 or aromatic carbon'
  - 'v1 also lacked the S-ca BOND while carrying S-ca angles — an incomplete set'
  - 'GAFF2 provides every needed term under its thioether sulfur type `ss`: c2-ss 213.76/1.7842, ca-ss 213.48/1.7847, cc-ss 231.66/1.7538, plus 136 terminal-S angles'
runbook: null
---

## Context

The covalent bond joins ff19SB's protein sulfur `S` to a GAFF2 carbon, and no
force field covers that pair, so the junction terms are supplied by a hand-built
frcmod. Version 1 was written while debugging a chloroacetamide, whose
attachment carbon is sp3 (`c3`), and it covered that case only.

Six of the nine warhead classes do not attach through sp3 carbon. The Michael
acceptors and BDHI attach through sp2 (`c2`), the SNAr azine through aromatic
carbon (`ca`/`cc`). All 18 of their MM-GBSA builds failed on missing `S-c2`,
`S-cc` and `S-ca` terms. The nine that succeeded were exactly the three SN2
acetamide classes — the sp3 case v1 was built for.

Version 1 also carried `S-ca` *angles* without the `S-ca` *bond*, which is
incoherent on its own terms and would have failed the aromatic case even if the
class had been considered.

## Decision

**Every junction term comes from GAFF2's own parameter for the same geometry,
with GAFF2's thioether sulfur type `ss` substituted for the protein's `S`.**

That substitution is an argument, not a convenience: once Cys113's SG has bonded
the ligand it sits between two carbons, `CB-S-C_lig`, which is precisely the
chemistry `ss` describes. GAFF2 supplies all of it — five bonds across the
carbon types that occur, and 136 terminal-S angles — so nothing is hand-chosen
and each line cites the `ss` form it came from.

**All 27 candidates are rebuilt under v2, including the 9 that already
succeeded.** Their v1 results are set aside as
`result_SUPERSEDED_junction_v1.json` rather than deleted. Keeping them would
have meant the acetamide classes rested on parm19-derived parameters while every
other class rested on GAFF2-derived ones — and while D0020 and D0023 already
forbid comparing dG across warhead classes, a parameter set that differs by
class makes even the within-class numbers rest on different footings depending
on which class you are in. One source, uniformly applied, costs a re-run and
removes the question.

## Consequences

The junction remains the largest modelling assumption in the MM-GBSA module: it
is an approximation at the exact bond the calculation is about. What changed is
that it is now a *single, cited, uniform* approximation rather than an ad-hoc set
that silently covered one third of the classes.

The failure was loud, which is the only reason it was cheap. tleap refused to
build and named the missing atom types; had it substituted a default and
proceeded, 18 candidates would have carried quietly wrong energies into a
ranking.

If the junction is ever suspected of driving a result, the check is to re-run one
class with a deliberately perturbed junction parameter and confirm the
within-class ordering is unchanged. That has not been done.
