---
id: D0034
title: The gate token erased sibling metrics — the same defect as before, one level down
date: 2026-07-28
status: accepted
approach: shared
decided_by: '@mhallet'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - shared/enrichment_gate.py
  - tests/test_gate_token_merge.py
  - decisions/D0031-class-matched-decoys-remove-the-apparent-covalent-enrichment.md
  - decisions/D0032-mmgbsa-gate-and-the-power-floor-on-negative-verdicts.md
evidence:
  - "write_token popped the whole stratum before re-adding: `for s in {r.stratum for r in results}: by_stratum.pop(s, None)`"
  - 'run_mmgbsa_gate.py writes a single mmgbsa_dG result for the covalent stratum'
  - "so it deleted the covalent affinity_kcal verdict (ROC-AUC 0.537) that D0031 established"
  - "the live token was found carrying covalent metrics = ['mmgbsa_dG'] only"
  - "with docking gone, recommended_rank_metric became mmgbsa_dG (0.140) — a metric NO approach ranks on"
  - 'demonstrated directly: old logic leaves [mmgbsa_dG], new logic leaves [affinity_kcal, mmgbsa_dG]'
  - 'token rebuilt: covalent affinity_kcal 0.537 + mmgbsa_dG 0.260, recommended back to affinity_kcal'
---

# The gate token erased sibling metrics

## What happened

`write_token` merged at **stratum** granularity and replaced at **metric**
granularity:

```python
for s in {r.stratum for r in results}:
    by_stratum.pop(s, None)          # this run supersedes its own strata
for r in results:
    s = by_stratum.setdefault(r.stratum, {"metrics": {}})
    s["metrics"][r.metric] = r.to_dict()
```

`run_mmgbsa_gate.py` evaluates exactly one metric — `mmgbsa_dG` — and writes it
for the `covalent` stratum. The pop therefore discarded `affinity_kcal`, the
covalent docking verdict D0031 measured at ROC-AUC 0.537, and re-added only
MM-GBSA.

The live token was found in exactly this state: `covalent` held one metric.

## Why it matters more than a missing entry

`recommended_rank_metric` is chosen as the best-scoring surviving metric in the
stratum. With docking deleted, the only candidate left was `mmgbsa_dG` at
ROC-AUC 0.140 — so the token recommended ranking on a metric **no approach in
this project ranks on**. T_3 and T_4 both rank on `affinity_kcal`; T_1 and T_2
on `vina_affinity`. The token was advertising a rank metric that contradicted
what the pipeline actually did, and the verdict travelling on every ranked row
described the wrong measurement.

## This is the second occurrence of one defect

The docstring immediately above the faulty lines documents the first:

> **MERGE, DO NOT REPLACE.** This used to rebuild the payload from only the
> current run's results, so `run_enrichment_gate.py covalent` deleted the
> non_covalent verdict [...] Nothing errored; their ranking stage just reported
> UNGATED.

That fix raised the merge granularity from *token* to *stratum* and left the
destructive replace sitting one level down. The same class of run — one that
evaluates a subset — produced the same class of silent loss, against a comment
explaining the hazard.

Raising the granularity of a destructive operation is not the same as removing
it. The fix now supersedes exactly what a run computed: one metric.

## What was restored

The covalent stratum was rebuilt from data that already existed — D0031's
class-matched docking results and the D0033-corrected dG values:

| stratum | metric | ROC-AUC | verdict |
|---|---|---|---|
| non_covalent | vina_affinity | 0.535 | UNDERPOWERED |
| covalent | **affinity_kcal** | **0.537** | UNDERPOWERED |
| covalent | mmgbsa_dG | 0.260 | UNDERPOWERED |

`recommended_rank_metric` for `covalent` is `affinity_kcal` again. Note that
this restoration changes no scientific conclusion: docking at 0.537 was already
indistinguishable from chance (D0031), and MM-GBSA at 0.260 is still worse. The
damage was to the record, not to the finding.

## Guard

`tests/test_gate_token_merge.py` pins four properties: writing one metric keeps
its siblings, writing one stratum keeps the others, rewriting the same metric
still supersedes it, and the recommended metric is the best survivor. The
sibling test was demonstrated to fail against the pre-fix logic before being
committed, so it pins the actual defect rather than the current behaviour.

## The pattern, now stated generally

D0033 was `dict.get(key, 0.0)` making a parse failure indistinguishable from a
zero term. This is `dict.pop(key)` making a partial run indistinguishable from
a full one. Both are **silent-loss defaults**: an operation that quietly
substitutes an empty value where the correct response is either to preserve
what exists or to fail.

The project has now been bitten six times by silent substitution (D0025, D0028,
D0029, D0030, D0033, D0034). The recurring rule: **when a value is absent,
absence must be representable and visible — never coerced to a plausible
default.**
