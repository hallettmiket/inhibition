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

## 3.0.0 “Galena” — in progress

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
