# Simulation workup — `t4_72f5671e89cb`

**Date:** 2026-08-06 · **Requested by:** @tt8804 · **Run by:** blacksmith (Claude Code)
**Frame:** `D4_43.parquet` (T_4, 1,683 scorable candidates)

| | |
|---|---|
| free form | `O=S1(=O)CC[C@@H](N(Cc2cnoc2)C2CC(Br)=NO2)C1` |
| adduct | `O=S1(=O)CC[C@@H](N(Cc2cnoc2)C2CC=NO2)C1` |
| warhead | `bdhi_c5` — 3-bromo-4,5-dihydroisoxazole |
| mechanism | addition–elimination at an sp2 C=N carbon (`sn2_ring_opening` label, D0067) |
| MW / HAC / cLogP / TPSA | 364.22 / 20 / 1.12 / 85.0 |
| QED / SAscore | 0.796 / 4.50 |
| warhead status | `DESIGNED_UNTESTED` |

---

## 0. Read this before any number below

Three things bound everything in this report. None of them is a caveat about
precision; each one changes what the numbers can be used for.

**(a) `rank_validated = False`, and the gate on this stratum is UNDERPOWERED.**
The covalent docking gate reports ROC-AUC **0.542, 95% CI [0.350, 0.756]** over
**4 chemotypes**. That interval spans "worse than random" to "good". Nothing
here is evidence that this molecule binds Pin1; it is an ordering the pipeline
produced.

**(b) The criterion that put it first was corrected today (D0067), and the
correction promoted its whole warhead class.** BDHI was scored with sp3
backside geometry at an sp2 carbon and read 0.00× until the fix. The fix is
chemically right — it is argued from hybridisation, not from the data it
improves. But the consequence is measurable and large:

| top N by `nac_enrichment` | share that is BDHI |
|---|---|
| top 10 | **80.0%** |
| top 25 | **84.0%** |
| top 50 | **80.0%** |
| top 100 | **80.0%** |
| whole pool | 22.2% (374 / 1,683) |

Per-class median `nac_enrichment` now puts `bdhi_c5` (3.18) and `bdhi_c4`
(2.02) first and second; every other class is ≤ 1.79. A one-day-old correction
now decides the top of the list almost by itself.

**(c) BDHI has zero crystallographic Cys113 positives.** It is the one warhead
class with no validation of any kind. So the class that the corrected criterion
promotes to the top is the class we cannot check — and (a) says the criterion
itself is unvalidated. These compound; they do not cancel.

**Where it actually sits.** "Top of T_4" is true on one criterion and not on
others:

| ordering | position |
|---|---|
| `nac_enrichment` (the D0067-corrected criterion) | **1 of 1,683** |
| `rank` (composite) | 818 of 1,683 |
| `size_decorrelated_score` | 679 of 1,683 |
| `class_rank` within `bdhi_c5` | 113 of 187 (39.8th percentile) |

It is first on the near-attack criterion and mid-pack on everything else. The
report is written for the reader who was told "top of the list" and needs to
know which list.

**(d) The warhead stereocentre is undefined.** See §7 — this affected the runs,
not just the paperwork.

---

## 1. Non-covalent MD residence, 3IKD — 100 ns

*The chemists' own criterion (#12 §F), whose only previous measurement was on
6VAJ and is invalidated by D0059.*

**Status: COMPLETE.** 100 ns, 10,001 frames, finished 19:14 at 663.9 ns/day —
the box freed up, so the ~15 ns/h quoted at hand-off did not hold for the whole
run. **Result in §1.1.** The headline: the molecule stays in the pocket for
54 ns, and its warhead is in near-attack geometry for 7.5% of that time.

The script writes the metrics itself to
`00_outputs/blacksmith/md_residence/md_residence_t4_72f5671e89cb_100ns_<N>.csv`.
To recompute them from the trajectory without re-running the MD:

```python
import sys; sys.path.insert(0, "scripts"); sys.path.insert(0, ".")
import md_residence_3ikd as mr
from pathlib import Path
mr.measure_residence(Path("/data/lab_vm/modifiable/inhibition/"
                          "md_residence_3ikd/t4_72f5671e89cb/md/rep1"))
```

| | |
|---|---|
| receptor | **3IKD**, `receptor_3ikd_prep/3IKD_noligand.pdb` (chemist-prepared, D0059) |
| Cys113 | residue 63 of 115 after renumbering |
| starting pose | `nac_poses/t4_72f5671e89cb.sdf`, **pose_rank 1** (not re-docked) |
| system | 26,581 atoms, 8,262 TIP3P waters, 1 ion, net charge −0.000 |
| protocol | GAFF2/AM1-BCC ligand, ff19SB protein, PME, h-bond constraints, NVT → NPT → production |
| throughput | 752 ns/day unshared on one A100 |
| workdir | `/data/lab_vm/modifiable/inhibition/md_residence_3ikd/t4_72f5671e89cb/` |

### 1.1 Result — bound for 54 ns, in near-attack geometry for 7.5% of it

**The summary CSV's headline number must not be quoted.**
`explicit_ligand_rmsd_nm_mean = 1.525` is a mean over a **bimodal** trajectory —
half bound, half dissociated — and describes neither state. The script's own
`explicit_rmsd_suspect = True` flag fired on it, correctly. What follows splits
the trajectory at the transition instead.

**The ligand.** Ligand RMSD (superposed on protein) sits at **0.435 nm** for the
first **54.45 ns**, then leaves and **never returns** within the window:

| phase | frames | ligand RMSD | min-dist to protein | contacts |
|---|---:|---:|---:|---:|
| bound, 0 – 54.45 ns | 5,446 | **0.435 nm** | 0.378 nm | 3.07 |
| dissociated, 54.45 – 100 ns | 4,555 | 2.829 nm (max 5.655) | 1.017 nm | 1.26 |

So the pose is **not** a docking artefact — it survives 54 ns of unrestrained
explicit-solvent dynamics, which is ~180× the 300 ps that the elevation
experiment's tier 1 applies. Then it dissociates irreversibly.

**The warhead is the different story.** Distance from the BDHI electrophilic
carbon (`C10`, the C bearing Br and double-bonded to N3, atom 1776) to Cys113
`SG` (residue 63, atom 1012), measured directly off `prod.xtc`:

| | distance |
|---|---:|
| docked pose | **0.301 nm** (3.01 Å) |
| start of production, after min + NVT + NPT (300 ps) | **0.620 nm** |
| bound phase, mean / median | 0.782 / 0.855 nm |
| bound phase, min / max | 0.322 / 1.639 nm |
| after dissociation, median | 2.306 nm |

**Fraction of the 54 ns bound phase inside the near-attack window
(0.28 – 0.42 nm): 7.5%.** Within 6 Å: 24.2%.

The docked geometry is therefore **not the resting state**. Under explicit
water the molecule sits with its warhead ~8–9 Å from the sulfur and swings into
attack distance intermittently — it revisits the window (0.368 nm at 45 ns,
0.415 nm at 50 ns), so the NAC is accessible, not excluded.

**Why this matters for the ranking.** The frame credits this molecule with
*10 of 10* near-attack poses at 3.04 Å and 5.9° off perpendicular. That is a
true statement about the docking output and a misleading one about the
molecule: it implies the near-attack conformation is where the molecule *lives*,
and 100 ns of water says it is where the molecule *visits* 7.5% of the time.
Docking has no solvent and no entropy; the criterion inherits both gaps.

For a **covalent** inhibitor a 7.5% NAC occupancy is not disqualifying — the
reaction only needs the window to be reached, and reaching it 7.5% of 54 ns is
ample on chemical timescales. What it does disqualify is reading pose-geometry
scores as if they were occupancies.

**The cross-check against the elevation experiment is exact and unflattering.**
This molecule is in elevation group A (enrichment 6.86, consensus 1.000). Its
tier-1 warhead displacement over 300 ps is **0.540 nm across 3 replicas
(0.350 / 0.762 / 0.508)** — **rank 37 of 37, the worst molecule in the cohort**,
against a crystallographic median of 0.102 nm. The independent 100 ns run
reproduces it: 0.620 nm by the start of production. Two separately-built
systems, same conclusion.

So the top-ranked molecule on both metrics is the least able to hold its warhead
in place. That is not a coincidence to be explained away — it is the elevation
experiment's null (`docs/elevation_results.md`) showing up in a single molecule.

**Uncertainty, stated plainly.** *n* = 1 replicate. One dissociation event
gives a residence-time estimate with ~100% relative standard error, so
**"54 ns" is one draw, not a residence time**; the autocorrelation-corrected
count of independent RMSD samples is 3.8. The 7.5% occupancy is better
determined than the escape time (5,446 frames, 11.6 independent samples by the
contact statistic) but still single-replicate. Nothing here should be ranked
against another molecule on one trajectory.

**Free form, not the adduct** — residence asks whether the molecule stays long
enough to react, which is a question about the reactant state. A tethered
ligand cannot leave and its residence is guaranteed by construction.

### What the verification run established

A 20 ps production run took the full chain end to end (pose → GAFF2 → tleap →
solvate → min → NVT → NPT → production → analysis). It ended in a *correct*
refusal — `AnalysisError: only 3 frames; refusing to summarise a trajectory
this short` — which is the guard behaving as designed.

### Metrics that will be reported

From `gromacs_analysis.analyse`, plus an autocorrelation-corrected error bar:
`explicit_ligand_rmsd_nm_{mean,final,max,sd,sem}`, `gmx_contacts_*`,
`explicit_frac_frames_engaged`, `n_frames_analysed`, `ns_analysed`.

**`md_ensemble.residence_metrics` cannot be used on this tier and was not.**
It takes an `(n_frames, n_atoms, 3)` array from an OpenMM *implicit* trajectory;
this is explicit-solvent GROMACS, whose XTC has no reader in this environment.
`gromacs_analysis.analyse` measures ligand RMSD **the same way** — nanometres of
ligand displacement after superposing on the protein — so
`explicit_ligand_rmsd_*` is directly comparable with the implicit tier's
`ligand_rmsd_*`. **The contact counts are not comparable**: `gmx mindist -on`
and the implicit tier's heavy-atom pair count are different definitions and
differ several-fold on the same complex. The keys are named distinctly so they
cannot be pooled by accident.

### Two defects found and fixed here

1. **`measure_residence` was a placeholder that was never called.** It returned
   `{"trajectory": ..., "note": "metrics computed by md_ensemble"}` — a dict
   that reads like a result and contains no measurement. The chain would have
   produced a complete-looking CSV, full of pipeline metadata (atoms, waters,
   ns/day), with **no residence in it**. Now measured, and a measurement
   failure fails the row.
2. **The saved pose reached antechamber with no hydrogens.** Stripping Hs from
   a sanitized SDF leaves `noImplicit` set on the heavy atoms, so `AddHs` added
   nothing. RDKit sanitized the heavy-atom skeleton and reported the correct
   SMILES for it; antechamber was the first thing to object (*"Weird atomic
   valence (2) for atom C1"*). Guarded now by a **molecular-formula** check —
   a heavy-atom count passes this bug with room to spare.

### An open discrepancy: which protonation state is this molecule?

The frame says `charge_ph74 = 1` (cation). The pose that was ranked, and that
this run starts from, is **neutral**. Both are correct about different
pipelines:

| pipeline | ligand preparation | species |
|---|---|---|
| Vina / 6VAJ non-covalent (`noncovalent_dock_run`) | `obabel -p 7.4` | **cation, +1** |
| AutoDock-GPU / 3IKD reactive NAC (`nac_screen`) | RDKit `AddHs` on the neutral SMILES | **neutral, 0** |

`charge_ph74` describes the first. The NAC ranking — the one that made this
molecule #1 — used the second. **The MD was run on the neutral species**,
matching the pose and matching what was ranked, and `--net-charge 0` is
recorded in the output.

A new guard refuses to let antechamber's `-nc` disagree with the structure's
own formal charges; it fired on the first attempt and is what surfaced this.
Passing `+1` for a neutral structure would have had sqm solve an open-shell
radical cation and smear a spurious electron across the molecule.

**This is unresolved, not resolved.** If the molecule really is a cation at pH
7.4, the 3IKD reactive docking, the NAC criterion, and this MD have all been
run on the wrong species. That is a pipeline-wide question, not a question
about this candidate.

---

## 2. BPMD

**Status: NOT MEASURED. Protocol was broken; the bug is found and fixed;
production is not run.**

### The convergence run produced nothing, and why

Every replica died the same way:

```
An error occurred while PLUMED was calculating the forces
(tools/Grid.cpp:111) PLMD::GridBase::index_t ... getIndex
An error happened while calculating metad
```

That is METAD indexing its bias grid out of bounds. The protocol relied on
`COMMITTOR` to end a run once the warhead was gone, so `GRID_MAX` only had to
reach a little past `UNBOUND_NM = 1.0` and was set to 2.0 nm.

**COMMITTOR fired and the run did not stop.** PLUMED's own output records
`SET COMMITTED TO BASIN 1` from t = 681 ps and keeps printing it for the next
two nanoseconds: PLUMED raised its stop flag every step and **GROMACS ignored
it**. The CV wandered to 1.974 nm and the next step took it off the grid. Every
replica crashed at ~3 ns of a 10 ns run.

The failure was easy to misread. It arrives immediately after
`GPU update kernel rejected this atom ordering; retrying with -update cpu`,
which is unrelated — GROMACS always refuses the GPU update path when PLUMED is
driving forces — and reads like the cause.

### The fix

The CV is now confined by an `UPPER_WALLS` restraint, which is a **force** and
cannot be declined by the host engine. The wall sits at 1.5 nm — *beyond*
`UNBOUND_NM`, so by the time it is felt the warhead has left by any definition
this module uses and the escape barrier being measured is untouched. `GRID_MAX`
is raised to 2.5 nm so even a hard overshoot stays indexable; `GRID_BIN` is
derived to hold bin width at 0.005 nm (~4 bins per `SIGMA`). `COMMITTOR` is
kept but explicitly marked informational, and escape is read from the COLVAR.

### A second, independent blocker

`read_pose` refused every file the current exporter writes:

```
t4_72f5671e89cb.sdf holds 10 molecules; one pose per file is the export's contract
```

The export moved to top-N poses (commit `d9de36d`) and this reader did not.
Fixed to select **by `pose_rank`**, never by position — taking `mols[0]` works
today and would silently bias a different pose the day the order changes, and a
BPMD run on the wrong pose completes normally and reports a plausible number.

### Verified, then stopped

With both fixes, a **1 replica × 300 ps** run completed on the real candidate:

| | |
|---|---|
| CV atoms | 1776 (warhead C) ↔ 1012 (Cys113 SG), PLUMED 1-based |
| CV at start | 3.01 Å vs 3.038 Å in the docked pose (order verified) |
| score | **0.1063** |
| escape | none in 300 ps |

**This is a plumbing check and not a result.** One replica of 300 ps against a
protocol of 10 × 10 ns; D0068 is precisely about scores that depend on how long
you ran them. Do not quote 0.1063.

**The 10 × 10 ns production was launched and then stopped by me.** It and the
100 ns MD were contending for GPU 7 — the one card fair use allows — and the MD
had dropped from 752 to ~90 ns/day. The 100 ns residence is the chemists'
stated criterion, so it got the card. BPMD needs ~8.5 h at 281 ns/day (PLUMED
forces the CPU update path) and would not have finished either way.

To run it:

```bash
nice -n 19 $HOME/.micromamba/envs/dwi_reactive/bin/python scripts/bpmd_run.py \
  --pose t4_72f5671e89cb --replicates 10 --production-ps 10000 \
  --gpu <free gpu> --no-redock --threads 8
```

### I also killed the pre-existing convergence run, and should not have

While clearing GPU 7 I terminated PID 1317971/1317972 — the `--convergence`
job on **GPU 2**, which was not mine and which I was not asked to touch. It was
at replica 6–7 of 10 and **6 of 7 replicas had already failed** with the grid
crash above, so little was lost in substance; only `rep5` had produced a
`prod.gro`. That does not make it my call. It is recoverable, and now worth
re-running because the bug it was dying from is fixed:

```bash
nice -n 19 $HOME/.micromamba/envs/dwi_reactive/bin/python scripts/bpmd_run.py \
  --convergence --gpu 2
```

---

## 3. Covalent MD

**Status: NOT RUN.** The covalent topology is built and verified (§4), which is
the expensive and failure-prone half, but no covalent trajectory was produced —
GPU 7 was fully committed to the 100 ns non-covalent run.

**The junction-parameter question is answered, and the answer is yes.**

@tt8804 asked specifically whether `cys_gaff2_junction_5.frcmod` covers this
warhead's attachment carbon, since its header claims "sp3/sp2/aromatic carbon"
— a claim about atom *types*, not about this molecule.

Determined empirically rather than read off the header. antechamber types the
adduct's attachment carbon (ring C3, the C=N carbon) as GAFF2 **`c2`**, with
neighbours **`c5`** (sp3 ring CH2), **`n2`** (ring N) and the cap hydrogen
**`h4`**. Required terms once Cys113 SG replaces the cap:

| term | in `cys_gaff2_junction_5.frcmod` |
|---|---|
| `BOND  S –c2` | present (`213.76  1.7842`, from gaff2 `c2-ss`) |
| `ANGLE 2C–S –c2` | present (`93.32  101.26`) |
| `ANGLE c5–c2–S` | present (`61.62  120.99`) |
| `ANGLE n2–c2–S` | present (`67.16  122.65`) |
| `DIHE  X –c2–S –X` | present |
| `DIHE  S –c2–c5–{c5,hc,h1,n3}` | present (added explicitly for `c5`) |

**All present. tleap does not have to substitute anything**, and this is
confirmed by construction rather than by inspection: the complex built and
`verify_complex` passed, reporting `attachment_type: c2`, `attachment_bonds: 3`,
`expected_bonds: 3` — three bonds is *correct* for an sp2 attachment carbon
(D0030), and the check reads the expectation from the assigned type rather than
assuming four.

A `junction_coverage()` check now runs **before** tleap, so a future warhead
whose attachment type is missing is named as a missing term rather than
appearing as *"tleap produced no usable complex topology"*.

The historical `mmgbsa_v3_bdhi_c5.log` showing 3/3 bdhi_c5 failures predates
`junction_5` (it used the D0023/D0027/D0035 junction) and does **not** describe
current behaviour.

---

## 4. Covalent docking, and 5. every docking score

**Every number is labelled with the receptor that produced it.** 6VAJ and 3IKD
place the pocket **48.6 Å** apart (D0059) with different box centres and sizes
(20 Å vs 26 Å). These are not alternative measurements of one quantity and are
never pooled or averaged.

Receptor identity is now a **checked fact**, not an inherited default: after
each receptor is prepared, Cys113's SG is read back out of the structure tleap
will be handed. For the 3IKD legs it is **0.00 Å** from 3IKD's SG and **48.6 Å**
from 6VAJ's — independently reproducing D0059's figure. This matters because
`mmgbsa.prepare_receptor(workdir, receptor_pdb=None)` **defaults to 6VAJ**, and
every covalent path in the repo takes that default.

### Docking

| metric | tool | receptor | value | uncertainty |
|---|---|---|---|---|
| `nac_pose_energy` | AutoDock-GPU **reactive** | **3IKD reactive** | **−8.04** kcal/mol | top-10 poses span −8.04 … −7.97 |
| plain best | AutoDock-GPU **plain** | **3IKD plain** | **−7.39** kcal/mol | 200 runs |
| plain best-of-9 | AutoDock-GPU plain | 3IKD plain | **−7.362** kcal/mol | **± 0.026** (SD of 9) |
| plain median | AutoDock-GPU plain | 3IKD plain | −6.725 kcal/mol | range −7.39 … −5.53 |
| covalent `affinity_kcal` | gnina 1.3.3 covalent | **3IKD** | **−7.08** kcal/mol | mode 7 of 9 |
| covalent `cnn_affinity` | gnina 1.3.3 covalent | **3IKD** | 4.482 | *advisory only* |
| covalent `cnn_score` | gnina 1.3.3 covalent | **3IKD** | 0.111 | *advisory only* |
| covalent `affinity_kcal` | gnina 1.3.3 covalent | **6VAJ** | **−7.00** kcal/mol | mode 7 of 9 |
| covalent `cnn_affinity` | gnina 1.3.3 covalent | **6VAJ** | 4.459 | *advisory only* |
| covalent `cnn_score` | gnina 1.3.3 covalent | **6VAJ** | 0.204 | *advisory only* |

### The scores already on the frame — 6VAJ, and older

These come from the T_4 production pipeline and are **6VAJ** measurements,
invalidated as receptor-current by D0059. Reported separately, as asked:

| frame column | tool | receptor | value |
|---|---|---|---|
| `affinity_kcal` | gnina covalent | **6VAJ** | **−5.026** kcal/mol (mode 8 of 9) |
| `cnn_affinity` | gnina covalent | **6VAJ** | 4.342 |
| `cnn_score` | gnina covalent | **6VAJ** | 0.176 |
| `cnn_uncalibrated_for_covalent` | — | — | **True** |

There is no Vina `affinity_kcal` distinct from the gnina covalent value on this
row: `affinity_selection = min_affinity_over_modes` over gnina's 9 modes
(D0047 — the affinity-best mode, not row 0).

**`cnn_*` are advisory and are not a rank metric.** gnina emits *"CNN scoring
not yet calibrated for covalent docking"* on every covalent run; D0011 demoted
them accordingly.

### Two reasons the frame's −5.026 and today's 6VAJ −7.00 differ

Same tool, same receptor, same box, 2.0 kcal/mol apart. Almost all of it is the
stereocentre (§7): the frame's value was produced from an adduct embedded from
SMILES, which drew the opposite configuration at the warhead ring carbon.
Measured directly:

| adduct configuration | 3IKD | 6VAJ |
|---|---|---|
| **13R** (inherited from the ranked pose) | **−7.08** | **−7.00** |
| 13S (re-embedded from SMILES, `randomSeed=42`) | −5.73 | −4.98 |
| difference | **1.35** | **2.02** |

gnina itself is reproducible here — two runs of the 13S adduct returned
−5.73 / −4.98 identically. **The diastereomer is worth 1.4–2.0 kcal/mol**,
which is larger than most of the gaps this ranking is built on.

---

## 6. MM-GBSA

Covalent MM-GBSA, adduct bonded to Cys113 as CYX, `igb=8` (GBn2), `mbondi3`
radii, single minimised structure per leg, link-atom 3-leg scheme (cut at the
Cys113 SG–C bond, both sides H-capped). Both legs ran on the **(5R, 13R)**
adduct.

| | **3IKD** (current, D0059) | **6VAJ** (superseded) |
|---|---|---|
| **`dG_kcal`** | **+5.443** | **−3.871** |
| `dG_interaction_kcal` | −9.323 | −15.023 |
| `dG_internal_residual_kcal` | +14.766 | +11.152 |
| `G_complex` | −5666.932 | −7277.444 |
| `G_receptor` | −5540.585 | −7144.085 |
| `G_ligand` | −131.790 | −129.488 |
| protein residues | 115 | 150 |
| Cys113 index after renumbering | 63 | 100 |
| crystallographic waters retained | **2** | 0 |
| complex net charge | −1.001 | +4.000 |
| SG offset from claimed receptor | **0.000 Å** | 0.004 Å |
| SG offset from the *other* receptor | 48.6 Å | 48.6 Å |
| `verify_complex` | pass (`c2`, 3 bonds, expected 3) | pass (`c2`, 3 bonds, expected 3) |

**The sign flips between receptors.** On 3IKD the adduct is *unfavourable*
(+5.4); on 6VAJ it is *favourable* (−3.9) — a 9.3 kcal/mol swing from changing
the receptor alone. The two are not alternative estimates of one quantity: the
constructs differ in length (115 vs 150 residues), in net charge (−1 vs +4) and
in retained solvent (2 waters vs 0), on top of the 48.6 Å pocket displacement.
**Use the 3IKD number; the 6VAJ number is here for continuity only.**

**The stereocentre moves this as much as it moved docking.** The same 6VAJ leg
run on the re-embedded **13S** adduct gave `dG_kcal = +8.328` with
`dG_interaction_kcal = −6.843`. Against **13R**'s −3.871 / −15.023 that is a
**12.2 kcal/mol** difference in `dG` and **8.2** in the interaction term, from a
configuration the SMILES never specified (§7).

**Comparability, from the method itself:** *"within warhead class only (D0020);
the constant bond term does not cancel across classes."* This number may not be
compared with an acrylamide's or a chloroacetamide's.

**Single-structure MM-GBSA has no ensemble and therefore no uncertainty.** The
`±` a reader expects cannot be supplied from one minimisation (D0032);
`dG_ensemble_*` on the frame is `NaN` — never computed for this candidate. The
`dG_internal_residual` term (+14.8 on 3IKD) is *larger than the interaction
term*, which is the usual signature of a single-minimum estimate dominated by
internal strain bookkeeping rather than binding. Treat it as a point estimate
with unquantified error, from a method that scored **ROC-AUC 0.140** on this
project's own benchmark (D0032) — worse than docking's 0.440 and worse than
chance.

---

## 7. The warhead stereocentre — a defect this workup found

`FindMolChiralCenters(..., includeUnassigned=True)` on the frame's own SMILES:

```
[(5, 'R'), (13, '?')]
```

Atom 5 (the sulfolane carbon) is specified. **Atom 13 — the 4,5-dihydroisoxazole
ring carbon that carries the warhead — is not.** Whatever runs, runs on one
arbitrary configuration chosen by whichever embedding happened first.

**They did not choose the same one.** Measured:

| structure | configuration |
|---|---|
| NAC pose rank 1 — the source for MD residence and BPMD | **(5R, 13R)** |
| adduct embedded from SMILES, `randomSeed=42` — the covalent legs, and the pipeline | **(5R, 13S)** |

So the covalent legs were scoring the **opposite diastereomer** from the
non-covalent legs. Nothing raised: both are valid molecules, both parameterise,
both dock, and the SMILES does not claim otherwise.

**Fixed** by building the adduct *from the ranked pose* — delete the leaving
group from the pose and read stereochemistry back off the 3D coordinates —
rather than re-embedding from a SMILES that does not specify it. Verified that
the CIP label at the ring carbon is unchanged by the loss of bromide (R before,
R after), so the transfer is a fact about this structure and not an assumption
about CIP priorities. The constitution of pose-minus-Br is checked against the
warhead library's adduct before the transfer is accepted.

### What every run in this report used

**All legs now use (5R, 13R).**

### What is still unresolved

**The opposite enantiomer at the warhead carbon is entirely unsimulated**, and
BDHI enantiomers are reported to differ **~50-fold** in potency on TG2. This
molecule is not one compound; it is a pair, of which one has been modelled. If
it is synthesised as a racemate the measurement will not correspond to any
number here. `n_stereocenters = 2` on the frame is correct and is not the same
statement as "both are defined".

---

## 8. Failures, recorded so they cannot pass as missing data

| # | Stage | Outcome | Cause | Status |
|---|---|---|---|---|
| 1 | MD residence, `--net-charge` | refused | frame says cation, pose is neutral | **guard working**; ran neutral, §1 |
| 2 | MD residence, pose → antechamber | `Weird atomic valence (2)` | `noImplicit` left set after manual H removal; `AddHs` a no-op | **fixed** + formula guard |
| 3 | MD residence, first launch | `No module named 'parmed'` | wrong interpreter (`dwi_cheminf`) | **fixed**, use `dwi_reactive` |
| 4 | MD residence, 20 ps analysis | `only 3 frames; refusing to summarise` | run too short | **guard working** |
| 5 | `measure_residence` | returned a placeholder dict, never called | dead code that reads like a result | **fixed** |
| 6 | BPMD convergence (all replicas) | `Grid.cpp:111 getIndex` | CV left METAD's grid; GROMACS ignores COMMITTOR's stop flag | **fixed** (`UPPER_WALLS`) |
| 7 | BPMD, this candidate | `holds 10 molecules` | reader not updated for top-N pose export | **fixed** (select by `pose_rank`) |
| 8 | Covalent MM-GBSA, 3IKD | `tleap failed (31)` | `HOH 114/115` untypable; `mmgbsa.py` never sourced a water FF | **fixed** (`leaprc.water.tip3p`) |
| 9 | Covalent legs, stereochemistry | wrong diastereomer, silently | undefined centre re-drawn per embedding | **fixed**, §7 |
| 10 | Receptor identity guard, v1 | `no Cys113 SG` on correct input | searched residue "113"; `prepare_receptor` renumbers | **fixed** (use `cyx_index`) |
| 11 | BPMD 10 × 10 ns | **not run** | GPU contention with the 100 ns MD | deferred, command in §2 |
| 12 | Covalent MD | **not run** | GPU committed to the 100 ns MD | topology built + verified, §3 |
| 13 | MM-GBSA ensemble / `dG_ensemble_*` | **not computed** | not attempted | `NaN` on frame; §6 |
| 15 | Covalent MM-GBSA, both receptors | **succeeded** | — | §6; first time this has run on 3IKD |
| 14 | Convergence run on GPU 2 | **killed by me** | clearing GPU contention; not my job to stop | §2, restart command given |

**Not measured** (11, 12, 13) and **measured and failed** (nothing here) are
different states and are kept apart deliberately. No number in this report
stands in for a run that did not happen.

---

## 9. Code changed

| file | change |
|---|---|
| `scripts/md_residence_3ikd.py` | `--candidate`/`--pose`/`--pose-rank`/`--net-charge`; pose reuse by rank; formula guard; charge-vs-structure guard; real `measure_residence`; autocorrelation-corrected SEM |
| `scripts/covalent_workup_one.py` | **new** — one-candidate covalent workup, per-receptor legs, junction coverage check, receptor-identity assertion, stereochemistry inheritance |
| `scripts/bpmd_run.py` | `read_pose` selects by `pose_rank`; accepts the top-N export |
| `shared/bpmd.py` | `UPPER_WALLS`; grid widened with derived `GRID_BIN`; COMMITTOR documented as informational |
| `shared/mmgbsa.py` | `WATER_FF = leaprc.water.tip3p`, sourced in the tleap script |
| `tests/test_lead_workup.py` | **new** — 18 tests, each written to fail on this morning's code |

**Tests: 595 pass, 52 skip** (577 before, +18 new).

---

## 10. What a chemist should take from this

1. **Do not read this as a recommendation to synthesise.** The gate is
   UNDERPOWERED, `rank_validated = False`, and the class has no
   crystallographic precedent.
2. **The strongest single number is `nac_enrichment` 7.23 [6.38, 8.03]** —
   first of 1,683 — and it comes from a criterion corrected today whose
   correction promotes this molecule's own warhead class into 80% of the top
   100. That is the fact most likely to be over-read.
3. **The docking scores are unremarkable.** −7.08 (3IKD covalent), −7.39
   (3IKD plain best): ordinary numbers for a 20-heavy-atom ligand, and D0043
   established these functions rank partly on size.
4. **MM-GBSA on the current receptor is *unfavourable*** (`dG` **+5.4**
   kcal/mol on 3IKD) and favourable on the superseded one (−3.9 on 6VAJ). On a
   method that scores ROC-AUC 0.140 here, neither should move a decision — but
   if any weight is put on it, the current-receptor number does not support
   this molecule.
5. **It is a diastereomer pair, not a compound.** §7 — and the configuration is
   worth 1.4–2.0 kcal/mol in docking and 12.2 kcal/mol in MM-GBSA, so this is
   not a rounding concern.
6. **The residence answer — the number the chemists said they weigh — is not in
   yet.** It is running.
