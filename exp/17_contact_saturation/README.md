# exp/17 — does the contact-space group count taper with docking depth?

**Verdict: no, and it does not need to.** The count climbs as `n^0.69` with no
plateau at any tolerance that keeps groups tight — but the region it is counting
inside is fixed (diameter exponent `+0.019`), and the groups themselves never
move (100% of n=500 groups persist at n=6,000, displaced 0.254 Å against a 0.73 Å
tolerance). The climb is undersampling, not expansion. Record:
[D0092](../../decisions/D0092-contact-space-is-fixed-the-group-count-climbs-because-6000-poses-undersample-it.md).

## Run

Environment `/data/lab_vm/envs/dwi_admet`. Each script is a few minutes; the
distance matrix is a `pdist`, so depth is cheap.

```bash
python run_all.py           # the ladder, plus 12 shallow ladders and the independence check
python tolerance_sweep.py   # the same ladder at 7 tolerances, one distance matrix per rung
python persistence.py       # do n=500 groups survive at n=6,000?
python space_growth.py      # extent, effective dimension, coverage cost
python plots.py             # the six-panel figure (run the four above first)
```

Outputs: `00_outputs/blacksmith/contact_saturation/`, figure at `figures_*.png`.

To look at a grouping rather than read it off a table:
`bash integration/run_pose_group_viewer.sh` (port 8932).

## Reading the numbers

* **`b`** is the exponent in `groups ~ n^b`. `1.0` = every pose is a new group,
  `0.0` = flat.
* **Never quote the group count as a number of binding modes.** It is a monotone
  function of docking depth. Compare rankings only at fixed depth.
* **The bounds are entailed, not measured.** Contact coordinates are capped at
  `pose_contacts.CAP_A`, so the region is bounded before anything is docked, and
  complete linkage at a fixed tolerance cannot exceed the covering number. D0091
  is the record of reporting that kind of bound as a finding. The rates are the
  content.

## Caveats

* Extent, dimension and persistence are **one molecule** (`t4_716800c125a7`) —
  the only 6,000-pose cloud in existence. Only the exponent is replicated (12
  molecules, ladders to ~450 poses).
* Coverage estimates use `N·ln N` and assume uniform cells; 65% of the cloud sits
  in 30% of the groups. They are an order of magnitude, not a spec.
* The shallow ladders read `<topic>_allposes`, which is DBSCAN-cleaned — see
  [D0093](../../decisions/D0093-the-file-named-allposes-is-not-all-poses-it-is-dbscan-cleaned.md).
  The deep ladder does not, and every headline number above is from the raw cloud.
