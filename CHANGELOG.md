# Changelog

Versions describe **the discovery method**, not a public API. The rule and its
justification are in [`docs/versioning.md`](docs/versioning.md); in short:

- **MAJOR** — previously reported numbers are **invalid** and must be re-measured.
- **MINOR** — new capability or metric; existing numbers stay valid.
- **PATCH** — a defect fix that corrects numbers within an unchanged definition.

The load-bearing distinction is between invalidating a *measurement* and
invalidating an *interpretation*. Replacing the receptor (D0059) killed the
measurements — major. Discovering that a metric predicts nothing (D0071) killed
an interpretation while the measured values stayed correct — not a major.

Every entry below states whether prior numbers survive it.

---

## 3.1.0 — in progress

### 2026-09-02 — the sulfur moved and nobody followed it

**MAJOR.** Cys113 docks as a **flexible** sidechain, so every pose has its own
SG. `nac_screen.sg_position` returned the FIRST docked model's sulfur and
`measure_poses` broadcast that one value across all 640 conformers while taking
ligand coordinates per conformer — so poses 2..N were measured against where the
sulfur used to be (D0109). **Every distance, angle, `in_range`, `viable`,
`viable_fraction`, `enrichment`, `anchor_quality` and `engagement` in `nac_v5`,
`nac_v6` and `nac_v7` must be re-measured.**

The error reads SHORT, which is why it suppressed scores instead of inflating
them. Same cloud, same 208-group split, only the sulfur differing: `viable`
49 → 96, in-window 199 → 238, poses under 3.0 Å 96 → 7. Non-uniform, because it
scales with how far each run's sulfur wandered — a property of the docking, not
of the molecule. That makes it a ranking defect rather than an offset.

It surfaced only because a distance ranking with no lower bound put impossible
1.22 Å distances on top, shorter than a C–S bond. **PoseBusters had been
disagreeing all along and was right.**

Survives: the modes (contact linkage never touches distance-to-SG — verified
identical before and after), PoseBusters validity, and D0108's NO GO, which
measures in the MD frame through a separate code path.

- `sg_positions(dlg)` returns one sulfur per model; `sg_position` deprecated,
  not deleted, so stale callers fail by name.
- **`measure_poses` refuses a single SG for a multi-conformer molecule.** That
  is the fix — a bare `(3,)` raises unless `allow_static_sg=True` is passed
  explicitly, and a wrong-length array raises too. `tests/test_per_pose_sulfur.py`.

### 2026-09-02 — the screen samples to a fixed number of VALID poses

**MAJOR for `consensus` and anything built on it.** `consensus` = mode_size /
n_poses with n_poses = the number DOCKED — equal at 500, and verified equal
across all 34,059 `nac_v6` rows. But only PoseBusters-valid poses can join a
mode, and that fraction runs **0.812–0.982**: a molecule at 81% valid had a
consensus *ceiling* of 0.81 where one at 98% could reach 0.98. The denominator
was equal and the numerator's headroom was not (D0106).

`docking.target_pb_valid: 500`, `n_runs: 640` — sized from the measured pass-rate
distribution of the same 562 molecules, where the worst case observed is 0.812.
The kept set is the **first 500 valid in docking order**: the runs are
independent GA replicates, so that is an unbiased sample, where an energy cut
would make it "the best-scoring 500" and inflate `engagement` (exp/21: the best
25% of a cloud concentrates attack-ready poses 2.60x over a random 25%).

Nothing is deleted — all 640 poses keep their row and their place in the cloud.
`pb_kept` distinguishes "in the analysed 500" from "valid", which are now
different things.

### 2026-09-02 — two physics readouts fail their own positive control

**Not a major** — the measured values are correct; their interpretation was not.
BPMD's `escaped` is True for 7 of 7 runs ever completed at 10 ns, *including
sulfopin's crystal pose*; and the first 100 ns non-covalent control on 3IKD puts
the positive **worst of three** (mean ligand RMSD 0.803 nm against a candidate's
0.314). Neither may be used to rank or reject a covalent candidate (D0107).
Tier-1 warhead drift over 300 ps remains the one validated readout.

**First synthesis verdict on it: NO GO on `t4_80fbed3bdf1e`** (D0108). Nine
molecules already existed sharing its exact R-group and differing only in
warhead; on tier 1 it came 8th of 9 (0.281 nm) while two crystallographically
validated warheads on the identical molecule reached 0.057 nm — better than the
REF median of 0.102. BDHI still has zero crystallographic Cys113 positives.

- **BPMD ran `pose_rank=1` for every run ever made** (D0105). `read_pose`,
  `run_pose` and `prepare_pose` all supported the argument; neither call site in
  `main()` passed one and no CLI flag existed. `already_done` was keyed without
  it too, so adding the flag alone would have made a finished rank-1 run mark
  rank 11 as done.


### 2026-08-19 — the run completed, and the modes it ranked are mixtures

**The screen finished**: 147/147 triaged, 15/15 at 100 ns — 5 held, 1 held
unstably, 9 left. Then the modes underneath it did not survive inspection.

**MAJOR, pending.** Every per-mode number in `nac_v5` — `viable_fraction`,
`enrichment`, `conditional_eb` — is measured over a group that is not one pose.
Median mode spans 3.51 A, 42% have a viable fraction between 0.1 and 0.9, and
the largest holds 137 poses across 9.3 A. The cause is circular: the pipeline
clusters on the reactive atom's position and direction, which ARE the score's
distance and angle terms (D0088). The 100 ns trajectories stand; which pose
earned each one was chosen by this machinery.

Fixed outright:

- **The screen was not reproducible.** AutoDock-GPU was invoked with no `--seed`,
  so every run drew a different cloud — v4 and v5 ranked the same 504 molecules
  at rho = +0.43, agreeing on the winning sub-mode 22.6% of the time.
  `docking.seed: 42` (#77).
- **The persisted pose cloud could not be joined to its own table.** The SDF
  numbered poses by position and was never rewritten on a re-screen, so it
  described a different run than the measurements beside it. `pose_idx` is now
  written (#76).
- **AutoDock-GPU fails silently above ~2,000 runs**: at 5,000 it exits -6 with
  "stack smashing detected" *and still writes a .dlg*; at 10,000 it reports
  failure and exits 0. `dock()` now verifies exit code, output file and log.
- **The empirical-Bayes prior did nothing.** Method of moments fitted it at 2.17
  poses on this heterogeneous library — a prior worth two poses shrinks nothing.
  A floor of 10 moves rho(score, mode_size) from +0.143 to -0.016
  (`ranking.eb_prior_min_strength`).
- **The GUI reported the wrong run.** The rail took its headline number from the
  md CSV and its verdict from the trajectory sidecar; with replicates those were
  different runs, so one row read "6.324 nm max / HELD". Both now come from one
  place, ranked on mean RMSD with max beside it.

Spec changes:

- **Triage sweep 8 ns → 5 ns** (D0087). Inside D0085's own bootstrap CI
  (4.3–9.5 ns), and truncation is one-sided, so it can only admit extras.
- **The 100 ns "optimal" bar is 0.45 nm** and is now separate from the 0.35 nm
  sweep bar. They were the same number, which made "optimal" unreachable at
  100 ns by construction.
- **Four residence tiers**: optimal / held / held-unstable / left, the last split
  by a residence floor so a run that travels 5.8 nm and returns is not reported
  beside one that never moved.

Built and **not adopted**: `shared/pose_cluster.py` — one clustering step on pose
similarity alone (HDBSCAN), attack geometry used only to rank. The only rule
measured that never produces a bag (largest mode 14 poses, widest 3.91 A), and it
places the validated pose in an 8-pose group 1.5 A wide against the shipped
rule's 108 poses across 8 A. Held at `proposed`: 29% of poses become noise and
the validated pose was lost in 3 of 30 replicates (#78); adopting it requires a
re-screen (#79).


**MAJOR: no 3.0.0 number survives this release, and 3.0.0 never closed.** The
re-run is on a new topic (`nac_v5`) with a different molecule set, so nothing
from Galena is quotable beside it:

- **The protonation fix changes which molecules exist.** `protonate()` could
  only ever build a ±1 species — it protonated one site and tested the total —
  so every dication in the library was stamped un-dockable and dropped. All 60
  of D4's failures were BDHI (30 bdhi_c4, 30 bdhi_c5, zero acrylamide), so both
  BDHI arms had entered 3.0.0 **15% short against a full acrylamide arm**. The
  families now enter at 187/187/187 after the bis-electrophile gate. Any 3.0.0
  cross-family comparison ran on unequal denominators.
- **The screen is scoped at docking, not only at ranking.** 3.0.0 docked 1,783
  molecules to rank 594; this docks the 561 in scope.
- **Every stage writes to run-scoped directories** (`attack_sweep_<topic>`,
  `md_residence_<topic>`, `sweep_gaps_<topic>`, `mdprio_reports_<topic>`).
  Eleven instances were found where a page or a query read a previous screen's
  directory — five of them surfaced by the user rather than by a test. Guarded
  by `tests/test_runs_are_standalone.py`.
- **The spec is read from `config/target.yaml` by the code that runs.**
  `md.production_ps` had been read only by the pipeline *diagram* while the
  runner defaulted to 100 ps, and `SWEEP_PS` stayed a 10 ns literal after D0085
  moved the decision to 8 ns.

**What is unchanged from 3.0.0:** the receptor (3IKD), the near-attack
criterion, the two-stage pose splitting, and the cascade's gates
(8 ns / 0.35 nm / 100 ns, D0085).

**Version coined 2026-08-17** while writing the handover comparison; the work it
labels began with the topic bump to `nac_v5` on the same day.

---

## 3.0.0 “Galena” — superseded before closing

**MAJOR: no 2.2.0 number survives this release.** Four changes each alter what is
measured, and no combination of them leaves a prior value comparable:

- **Second-stage pose splitting** (#61) changes `consensus` = mode_size/n_poses,
  and every score computed from it. Sized by measurement: on the 82-case
  benchmark, carrying 4–5 representatives instead of 1 lifts crystal-pose
  recovery from 22.0% to 39–40% (14 cases gained, 0 lost, McNemar p = 1.2×10⁻⁴),
  and past 5 there is nothing left (k=5 vs k=9, p = 1.00).
- **Per-mode selection** (#53). 2.2.0 ranked per mode and then sent **mode 0 for
  242 of 242 molecules**; 5 modes ranking *first* in their warhead class were
  never simulated. The 2.2.0 shortlist is a shortlist of mode 0s.
- **Physiological ionic strength** (#57). Every 2.2.0 system was built with
  `addions … 0` — neutralise and stop — so **no simulation had salt**, at ~0 M
  against a cell's 0.15 M, with 726 of 1,782 T₄ candidates cationic at pH 7.4.
- **The pH 7.4 species is what gets docked** (#58, @tt8804). The pose set and the
  warhead library disagreed about protonation, which blocked the covalent workup
  on 41% of the library.

Prior numbers are not merely uncertain — they answer a different question.

### Also in this release
- **The rank gate counts poses, not a fraction of the cloud** (#65, D0084).
  `consensus >= 0.05` was exactly `n_poses_mode >= 25` on a 500-pose cloud — a
  number nothing measured, and one that would silently become 50 if
  `docking.n_runs` doubled. Replaced by the estimability rule already measured
  for the sweep, `>= 12` poses. T_4: 5,132 → 6,338 modes ranked, 0 lost; among
  modes clearing enrichment 4.0, 289 → 434. **No re-dock and no score changed** —
  only which modes may hold a position. The old rule stays reachable as
  `--gate consensus_fraction` and every row records which gate ran.
- **Modes are named by their own index in the ranking view** (#65, D0083). The
  lettered `0a`/`1b` labels implied that sub-splits of one first-stage mode are
  variants of each other; measured, 22% of split first-stage clusters hold modes
  whose median reactive-atom distance spans more than the criterion's entire
  2.8–4.2 Å window. Display-only — the identities were already flat.
- `config/target.yaml`: the target and every screen decision in one file, read
  only by `shared/target_config.py`, which **refuses** to supply a sweep floor
  that has not been measured for the target at hand (#59).
- `inchikey` and `docked_smiles` on every candidate frame — 1,783 T₄ and 5,396 T₃
  rows, 0 duplicate keys. A molecule whose pH 7.4 form cannot be built is stamped
  and excluded rather than quietly docked as the neutral.
- Pose clouds persist by rule, not by flag (#44, `CLAUDE.md`).
- Runs record which pose they simulated (#35, #36).

---

## 2.2.0 “Chalcopyrite” — closed 2026-08-11

Pose splitting and tooling upgrades. Outline: `docs/outline_2.2.0.md`.
Framework as built: `docs/framework_2.2.0.md`.
Retrospective on the release it follows: `docs/retrospective_2.1.0.md`.
Retrospective on this one: `docs/retrospective_2.2.0.md`, written 2026-08-09 and
judged against the outline's four pre-registered failure criteria — one failed
outright (another silent stage), one is untested (mode-count reproducibility).

### The positive control, and what it cost to read it (2026-08-10)

Sulfopin and Liu-2022-ZL-Pin13 were put through the pipeline as positive controls.
**No measured value changes for any candidate** — this is about what the readout
was doing to them.

Three records, in the order they were found, each qualifying the one before:

* **D0075** — the 10 ns sweep rejects *every* known active. Sulfopin is not
  rejected early (500 runs, 1 mode, 465 poses, 47 reaction-competent); it fails at
  the sweep, with zero sustained visits, ranked 104 of 234.
* **D0076** — *why*. `rx_7F0M` entered attack geometry **13 separate times** and
  scored zero, because `MIN_DWELL_PS = 100` requires each excursion to *last*
  100 ps at a 19.96 ps save interval. The pre-registration chose visits precisely
  because "a covalent reaction needs ONE good approach, not sustained occupancy",
  and the implementation filters on persistence. Re-derived on raw visits, within
  mechanism: **Liu-2022 beats 20 of 20 SN2 candidates**, Sulfopin 18 of 20.
* **D0077** — the `rx_*` controls were the wrong shape. 6VAJ's ligand is a
  covalent **adduct**; cleaving the bond leaves the reactive carbon 1.98 Å from
  SG, *below* the 2.8 Å window floor, and its 180° is constructed rather than
  measured. Equilibration relaxes it to 3.57 Å / 100.9° before production starts.
  Our own docked pose, by contrast, enters at 3.36 Å / 156.8° — inside the window
  and over the SN2 bar.

Net: the screen's chemistry was sound, and three separate readout defects made it
look otherwise. The 100 ns results stand — Liu-2022 **99.95%** engaged and held,
Sulfopin 78.9%, Juglone 47.5% — and the docked Sulfopin is the control that should
carry the claim.

**Deliberately not stamped released.** Two things are open and neither is the
retrospective's to decide: the receptor split (D0059 is still `proposed`, while
`config/receptor.yaml` and `noncovalent_dock_run.py` still default to 6VAJ), and
the next version number, which is the open question in
[#46](https://github.com/hallettmiket/inhibition/issues/46) — `2.3.0` or `3.0.0`
turns on whether "the ranking predicts nothing" invalidates measurements or only
their interpretation, the distinction this file opens with.

### The catalogue viewer — the new GUI foundation (2026-08-09)

@tt8804: *"we will use this as the new gui foundation."* Adopted as the pattern
every future results interface follows. **No measured value changes** — this is
presentation only, so every prior number survives it.

The shape: **a selector rail on the left, one viewer on the right.** Each row
carries a 2D structure, identifier, warhead class, sustained visits, max RMSD, an
engagement bar and a held/left tag. Selecting a row loads that molecule's full
report beside it — 3D pose, **MD movie**, **RMSD plots** — with everything in
collapsible panels that start closed, so the page paints immediately instead of
after ~9 MB of movie frames.

Two composable toggles, both of which encode a real distinction rather than a
preference:

- **all classes / by warhead class** — cross-class ranking is biased, because the
  SN2 angular criterion is far stricter than the perpendicular one (#47). The
  toggle makes that visible instead of something the reader must remember.
- **combined / split held-left** — engagement and residence are near-independent
  (rho = −0.007, #46), so a molecule can rank high and still leave the pocket.

Light/dark toggle stamps both the shell and the framed report, and persists.

Structural fixes worth carrying forward: the report writes its blocks straight
into `<body>` with no wrapper, so the measure has to be applied to those blocks;
`box-sizing: border-box` globally, or bordered cards and plain text resolve to
different left edges; and a 98%-opaque surface hides anything inside it, so the
ligand and Cys113 are **cut out of** the surface rather than merely drawn under
it.

Implementation: `scripts/mdprio_combine.py`, `scripts/mdprio_report.py`,
`shared/report_theme.py`, `shared/md_movie.py`. Full record in issue #49.

### PATCH — the movie drew the wrong residue as Cys113

`shared/md_movie.py` styled `resi: 113`. The MD system renumbers from 1, so the
crystal's Cys113 is residue **63** (`PIN1_OFFSET = 50`); residue 113 is a
**glutamate**, and it was rendered in sticks and labelled as the target cysteine.
`elevation_report.py` had the offset right and refuses to label a structure whose
residue types do not match; the new styling code bypassed that guard.

Same metric, same inputs, a previously wrong picture — a patch by the rule. It
affected **no computed value**: the warhead→SG distance series is built from the
SG atom located independently, and it was correct throughout. Only the rendering
was wrong. `CYS113_RESI = 113 - PIN1_OFFSET` now drives the surface cut-out, the
sticks and the sulfur sphere.

---

## 2.1.0 “Bornite” — 2026-08-07

The ranking rework. Screen re-run persisting per-pose geometry, gnina scores and
the poses themselves; consensus as a per-warhead-class quota (D0073); a weighted
anchoring score; selection that re-measures the pose it elevates; BPMD pose
ranking as a separate stage; references through the identical criterion, with
Sulfopin and ATRA through 100 ns.

**Numbered 2.1.0, not 3.0.0.** The open question when this work began was whether
the rework *redefined* `enrichment` — which would make old and new values share a
name while not being comparable, and would take a major. It did not: the 2.0.0
quantity is carried forward unchanged as `enrichment_joint` and the new
quantities sit beside it under new names. Existing numbers stay valid; what
changed is their standing in the ranking.

**Closed by its own finding (issue #23):** the pose window it scores is ordered by
docking energy, and energy carries no information about reaction geometry —
Sulfopin, which has a crystallographic Cys113 adduct, scores 0.000. The structure
carries forward; the score does not. Full account in
`docs/retrospective_2.1.0.md`.

Names come from the copper-mineral alphabet in issue #27 — see
[`docs/versioning.md`](docs/versioning.md).

---

## 2.0.0 “Azurite” — 2026-08-06

**MAJOR: every 1.0.0 number is invalid.** D0059 replaced the receptor, and 6VAJ
and 3IKD place the pocket 48.6 Å apart — prior values measure the wrong site.

The 3IKD receptor and the geometric (near-attack) ranking, end to end: screen,
rank, elevation suite, and the first full molecular workup.

*This release bundles what convention would have shipped as roughly one major,
four minors and several patches; the whole line was developed unreleased on one
branch. The audit is in [`docs/versioning.md`](docs/versioning.md) rather than
back-dated into tags that were never cut.*

**The headline result is a negative one, and it is the point of the release:**
neither of the two ranking metrics predicts whether a docked pose survives
physics (D0071), established on a pre-registered cohort with a crystallographic
anchor that separates at p = 0.007.

### Receptor
- **3IKD replaces 6VAJ**, used exactly as the chemist prepared it (D0059). Every
  6VAJ measurement is invalidated as receptor-current; the two are never pooled.
- Vina-GPU segfault when it cannot write its kernel cache (D0060).

### The ranking
- Whole-molecule RMSD is the wrong endpoint for a covalent question (D0062).
- No cheap pose-selection rule beats random (D0061).
- **Reactive docking replaces dock-then-filter** (D0063); the reactive potential
  is a *sampler*, not a criterion (D0064).
- **Mechanism-specific near-attack criterion** — SN2 backside vs perpendicular
  approach at sp² centres, each against an exact isotropic solid-angle null.
  Validated on two independent mechanisms: chloroacetamide **AUC 0.908**, Michael
  0.734 (D0065).
- Production ranking of all **5,769** candidates through the validated gate.
- Pose **consensus** — agreement among the top-N poses by energy — added as a
  second component (D0070), plus a composite that lets an unmeasured component
  contribute a `[0, 1]` interval rather than an imputed value.

### Defects found and fixed
- **D0067 — BDHI was scored with sp³ backside geometry at an sp² carbon.** The
  mechanism *name* was trusted over the hybridisation. 374 candidates read
  0.00×; after the fix BDHI became the top two classes of nine. Independently
  corroborated by Byun 2023's DFT.
- **D0066 — 97% of T_3 is an acryl*imide*, not an acryl*amide*** (corrected from
  an initial 77%), and the ranking does **not** prefer them.
- **D0068 — enrichment does not converge.** The same molecules fall 2.91× → 0.96×
  at 10× search effort; rank correlation across efforts ρ = −0.117. The cause is
  the *window*: dividing by every run puts every mediocre pose in the
  denominator. Scoring the top-N by energy instead is rank-stable (ρ = +0.568).
- **D0069 — plain docking on 3IKD still beats the geometric criterion.** The
  receptor was *not* the explanation for the earlier weak result.
- Pose atoms were keyed on the PDBQT name field (every carbon named `C`);
  reactive typing was guarded by a literal type name, silently deleting the whole
  SNAr class; the GUI drew every pose into 6VAJ regardless of where it was
  docked; the pose export was not idempotent.

### The elevation suite
- `docs/elevation_prereg.md`, committed **before** any simulation ran.
- Tiers 1–4 (`scripts/elevation_run.py`, `elevation_launch.sh`,
  `elevation_analysis.py`), fair use enforced in the launcher.
- **D0071 — neither metric predicts stability**, and the anchor is what makes
  that readable rather than merely disappointing.
- `docs/elevation_example.md` — the suite documented as a repeatable protocol.

### The worked molecule
- `t4_72f5671e89cb`: literature, ADMET/developability, 100 ns MD, BPMD, covalent
  docking and adduct topology, every docking score labelled by receptor.
- **D0072 — NO GO.** Rank 37 of 37 on tier-1 warhead stability, and its warhead
  class has zero crystallographic positives.

### Tooling
- BPMD (`shared/bpmd.py`) with the PLUMED kernel actually installed — GROMACS
  accepts `-plumed` with no PLUMED present, which had produced a false positive.
- GUI: near-attack ranking panel, 2D structure beside the 3D pose, per-pose
  receptor resolution, enrichment formula surfaced.
- `scripts/elevation_report.py` — self-contained HTML report; every figure and
  table recomputed from the shard CSVs at build time.

### Known-unfixed, carried forward
- `mmgbsa.RECEPTOR_PDB` still defaults to **6VAJ**, and every covalent path takes
  that default.
- `nac_rank.refine()` counts `failed:` rows as done when resuming.
- Covalent MD never ran.
- No tier-3 baseline on a crystallographic positive.

---

## 1.0.0 — up to the handoff from @mhallet

**Assigned retroactively; no 1.0.0 was ever cut.** It names everything before
@tt8804 took over — the 6VAJ receptor, the T_3/T_4 generative campaigns, the
original docking and MM-GBSA pipeline, decisions D0001–D0058 — so that
pre-handoff numbers have a pipeline to be attributed to.

**Numbers from 1.0.0 are 6VAJ measurements** and are not comparable to anything
in 2.x without re-measurement (D0059).
