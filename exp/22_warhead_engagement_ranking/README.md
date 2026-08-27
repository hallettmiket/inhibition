# exp/22 — which static metric predicts warhead engagement in MD?

**Verdict: the column we rank on today predicts nothing; the geometry of one pose
predicts strongly.** Against 147 swept modes, using measured `frac_attack_ready`
as the outcome:

| | ρ | p |
|---|---:|---:|
| the simulated pose's anchor quality | **+0.652** | 3.5e-19 |
| `anchor_quality_max` | +0.130 | 0.12 |
| **`conditional_eb` (incumbent)** | **−0.015** | **0.86** |
| `enrichment` / `viable_fraction` | −0.043 | 0.61 |

Record: [D0098](../../decisions/D0098-rank-on-warhead-engagement-geometry-the-incumbent-column-predicts-nothing.md).

## Run

```bash
E=/data/lab_vm/envs/dwi_cheminf/bin/python
$E run_all.py            # every candidate metric against the MD outcome
$E splitter_effect.py    # does contact grouping shrink within-mode spread?
```

`splitter_effect.py` needs two screens of the same molecules, one per method:

```bash
RX=$HOME/.micromamba/envs/dwi_reactive/bin/python
for m in contact_linkage warhead_dbscan; do
  $RX ../../scripts/nac_screen_v2.py --topic cmp_$m --only ids.txt \
     --nrun 300 --gpu 2 --no-gnina --split-method $m
done
```

## Reading the numbers

* **`lift` is the honest column.** A ρ of +0.13 sounds like something; a top-20
  hit rate of 60% against a 44% base rate is what it buys.
* **Always read `rho vs size` beside a metric.** `anchor_quality_max` looks
  second-best until you see it correlates with mode size at +0.638.
* **Mode-level aggregates fail because modes are mixtures** — median within-mode
  engagement spread is 0.776 on a 0–1 scale. Contact grouping cuts that to 0.288.

## Caveats

* **Range restriction.** These 147 were *selected* for sweeping by
  `conditional_eb`. Metric-vs-metric comparison is fair; absolute values are not
  population estimates.
* `frac_attack_ready` is **reachability of attack geometry**, not binding.
* The strongest metric needs the representative pose's geometry, which for
  `nac_v5` **cannot be recovered** — those files predate #76 and carry no
  `pose_idx`.
* No MD has been run under `contact_linkage`.
