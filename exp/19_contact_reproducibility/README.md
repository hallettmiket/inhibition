# exp/19 — do contact groups reproduce across independent dockings?

**Verdict: yes, for the groups that carry population.** 98.5% pairwise, and 91%
of ≥5-pose groups present in all five independent dockings — against HDBSCAN's
1 of 3 modes surviving the same test (D0088, #78). Record:
[D0095](../../decisions/D0095-contact-groups-reproduce-across-independent-dockings-and-survive-the-raw-cloud.md).

## Run

```bash
/data/lab_vm/envs/dwi_admet/bin/python run_all.py     # needs sklearn
```

## Reading the numbers

* **Singletons drag the all-groups rate to 59%,** and that is what a singleton is
  — one pose that landed somewhere, with no claim to being a mode. Read the
  `>= 5 poses` row.
* **The DBSCAN baseline is refused, not reported.** It returns one mode per
  replicate and "reproduces" at 100%; a rule that always answers "one mode" is
  perfectly reproducible and discriminates nothing. It is also circular on these
  clouds, which were already DBSCAN-cleaned (D0093).
* This is **stronger than exp/17's persistence**: a deeper draw from one cloud
  cannot move a fixed-tolerance region, but an independent docking can.

## Caveats

* One molecule (`t4_716800c125a7`) — the only one with five independent dockings.
* The clouds are `_allposes` files, ~21% DBSCAN-filtered (D0093). Both methods see
  the same filter so the comparison is fair; the absolute rates are the easy case.
* The tolerance used is itself under question (D0094).
