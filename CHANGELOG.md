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

## 2.2.0 “Chalcopyrite” — in progress

Pose splitting and tooling upgrades. Outline: `docs/outline_2.2.0.md`.
Retrospective on the release it follows: `docs/retrospective_2.1.0.md`.

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
