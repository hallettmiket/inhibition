# exp/21 — pose generation audit

**Verdict: docking is sound; the analysis was blind to the scores.** Poses that
look like they sit outside the pocket are real, rare (2.6%), and sit at the
**88th energy percentile** — the scorer ranks them correctly (ρ = +0.446 with
exposure). What was wrong is that the persisted clouds carried **no energies at
all**, so four clustering experiments and the viewer weighted the best pose and
the 500th equally. Record:
[D0096](../../decisions/D0096-pose-generation-is-sound-the-clustering-was-measured-without-the-scores.md).

## Run

Needs the **docking** environment (gemmi/meeko), not `dwi_cheminf`:

```bash
RX=$HOME/.micromamba/envs/dwi_reactive/bin/python
$RX run_all.py --n-molecules 6 --gpu 2        # geometry vs energy, PoseBusters
/data/lab_vm/envs/dwi_admet/bin/python energy_filtered_grouping.py
```

## Reading the numbers

* **The energy↔pose pairing is solved, not assumed.** AutoDock reports a cluster
  ranking beside the run order, so "the order" is ambiguous. `energies_aligned`
  matches records to conformers by an order-invariant signature under a Hungarian
  assignment and requires the identity permutation at ~0 Å. Two earlier versions
  were wrong — one counted the flexible Cys113 side chain as ligand atoms, one
  compared coordinates element-wise across differing atom orders — and both
  refused to run rather than guessing. That was the right failure mode.
* **Every filtered statistic needs the size-matched control.** Keeping the best
  25% leaves a quarter of the poses, and a smaller sample concentrates its top
  group for purely arithmetic reasons. The control is a random subset of
  identical size; only the gap between them is a finding.
* **`exposed_frac` is the useful burial measure**, not `contacts` — contacts
  scale with molecule size, exposure does not.

## Caveats

* 6 molecules for the geometry panel, 21 for the energy-filter panel.
* The energy is AutoDock's, on the reactive receptor — the same function D0041
  and D0046 show does not rank actives on this target. That poses agree with each
  other says nothing yet about whether they are *right*.
* This fixed the audit path. **`nac_screen_v2` still writes production clouds
  without energies.**
