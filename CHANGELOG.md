# Changelog

Versions describe **the method**, not the code — a version is the state of the
discovery pipeline you would have to name when quoting a number from it. The
boundary between 1.x and 2.x is the handoff from @mhallet to @tt8804.

Numbers produced under different major versions are **not comparable**: 2.0.0
changed the receptor, and 6VAJ and 3IKD place the pocket 48.6 Å apart (D0059).

---

## 2.1.0 — in progress

Ranking rework, opened after 2.0.0 established that the 2.0.0 ranking has no
physical support. See `docs/recap_2.0.0.md` for what carries forward.

---

## 2.0.0 — 2026-08-06

The 3IKD receptor and the geometric (near-attack) ranking, end to end: screen,
rank, elevation suite, and the first full molecular workup.

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

### Known-unfixed, carried into 2.1.0
- `mmgbsa.RECEPTOR_PDB` still defaults to **6VAJ**, and every covalent path takes
  that default.
- `nac_rank.refine()` counts `failed:` rows as done when resuming.
- Covalent MD never ran.
- No tier-3 baseline on a crystallographic positive.

---

## 1.0.0 — up to the handoff from @mhallet

Everything before @tt8804 took over the project: the 6VAJ receptor, the T_3/T_4
generative campaigns, the original docking and MM-GBSA pipeline, and decisions
D0001–D0058.

Retained for provenance. **Numbers from 1.0.0 are 6VAJ measurements** and are not
comparable to anything in 2.x without re-measurement (D0059).
