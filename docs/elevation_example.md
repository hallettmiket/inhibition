# The elevation suite — what it is, what we ran, and the worked example

*Pipeline **2.0.0**. @tt8804, 2026-08-06. Written after the first full elevation,
so that the next one is a repeat rather than an invention. The worked example is
`t4_72f5671e89cb`, and **its verdict is NO GO** (§8, D0072).*

---

## 0. What elevation is for

The screen ranks ~5,769 candidates on geometry alone: can the molecule present
its warhead to Cys113. **Elevation is the step between a rank and a synthesis
request.** It asks a different question — *does this pose survive contact with
physics* — using progressively more expensive simulation, and it exists because
the ranking is a docking artefact until something outside docking agrees with it.

The lab constraint that shapes everything here: **synthesis is one compound per
week.** Elevation is not a filter applied to hundreds of molecules. It is a
workup applied to a handful, and its output is a go/no-go with reasons.

**Elevation does not measure reactivity.** Every readout below is about whether
a *pose* is physically real. Whether the molecule then reacts with Cys113 is a
separate question that none of this touches.

---

## 1. The non-negotiable: a cohort, never a single molecule

The single most important thing this experiment established is procedural.

> **Never elevate one molecule alone. Elevate it inside a cohort that contains
> crystallographic positives.**

A stability number on one molecule cannot be interpreted. 0.54 nm of warhead
drift is neither good nor bad in isolation. It became meaningful only because
eight crystallographic Cys113 binders went through the identical pipeline in the
same batch and landed at 0.102 nm.

This is not a statistical nicety — it is the difference between a result and a
number. The tier-1 run produced a **null** on the question it was designed to
answer (neither ranking metric predicts stability, D0071), and that null is
*informative* purely because the anchor separated at p = 0.007. Without the
anchor, the same data would have been indistinguishable from a broken assay.

`crystal_positives()` in `scripts/nac_screen.py` supplies them. The membership
rule must be **fixed before the run** — ours was *sort by ident, take the first
8* — because "which 8 of the 15" is otherwise a knob that can be turned after
seeing the answer.

---

## 2. Pre-registration

`docs/elevation_prereg.md` was written and **committed to git before any
simulation started** (`efcfd8d`, before `73aed19`).

It fixes, in advance: the groups, the readouts, and — most importantly — a table
mapping each possible observation to the conclusion it forces. With four groups
and several plausible readouts, a result chosen after the fact can be made to say
almost anything, and D0045 exists in this project because that has happened.

When the outcome landed on the **least convenient row** of that table, the row
was reported as written. That is the whole point of the exercise; a
pre-registration only earns its keep on the run where it costs you something.

---

## 3. The suite, in the order it runs

| tier | what it asks | cost | verdict power |
|---|---|---|---|
| **0** | does the molecule exist, is it makeable, is the warhead sane | hours, no GPU | can veto outright |
| **1** | does the docked pose survive 300 ps of plain water | ~free — it is the equilibration you already run | strong, and cheap |
| **2** | how hard must you push the warhead out of the window | ~3 GPU-h per molecule | moderate |
| **3** | does the molecule stay in the pocket for 100 ns | ~4 GPU-h per molecule | descriptive at n=1 |
| **4** | covalent adduct: does the bonded complex behave | ~1 GPU-h + topology work | confirmatory |

**Tier 1 is the discovery of this experiment.** `gromacs_explicit` applies no
position restraints during NVT/NPT, so 300 ps of unrestrained dynamics *already
runs* before any production. The warhead's displacement across that window was
being computed and thrown away. It is the cheapest physical signal in the stack
and it is the one that separated the crystallographic anchor from every
generated group.

---

## 4. What we actually ran

All commands from `~/repos/inhibition`. Python is
`~/.micromamba/envs/dwi_reactive/bin/python` unless noted.

### 4.0 Prerequisites

```bash
# receptor — chemist-prepared 3IKD, NOT 6VAJ (D0059; they are 48.6 A apart)
python scripts/prepare_3ikd_receptor.py

# the screen the cohort is drawn from
bash scripts/run_consensus_all.sh          # consensus for all 5,769
python scripts/export_nac_poses.py         # top-10 poses per molecule -> SDF
```

### 4.1 Cohort selection

```bash
python scripts/elevation_cohort.py --per-group 8 --seed 0xE1E7A7
```

Writes `00_outputs/blacksmith/elevation_cohort/elevation_cohort_<N>.csv`.
Groups are **matched on rotatable-bond count**, not merely balanced, because
consensus correlates with rigidity at ρ = −0.259 (p = 4×10⁻⁸⁹) and an unmatched
draw would rediscover that instead of the thing being tested.

### 4.2 Pre-register, then commit

```bash
$EDITOR docs/elevation_prereg.md
git add docs/elevation_prereg.md && git commit    # BEFORE the next step
```

### 4.3 Tier 1 and tier 2

```bash
bash scripts/elevation_launch.sh 1              # equilibration survival
bash scripts/elevation_launch.sh 2 3000         # BPMD, 3 replicas x 3 ns
```

Four shards, one per GPU, in tmux session `elevate`. The launcher **refuses GPUs
0 and 7 by name** and puts `nice -n 19` on every worker — fair use is enforced in
the script rather than remembered, because this box carries ~35 users.

Logs tee to `/data/lab_vm/modifiable/inhibition/elevation_logs/`.
Tier 2 reuses tier 1's equilibration when the `.mdp` is identical, which is most
of why 111 BPMD replicas fit in an evening.

### 4.4 Analysis

```bash
python scripts/elevation_analysis.py --write
```

Mann-Whitney per contrast with Cliff's δ alongside, Holm across the three
pre-registered contrasts, every value carrying its replica spread.

### 4.5 The lead's deeper workup

```bash
# 100 ns non-covalent residence on 3IKD
bash scripts/run_md_residence.sh                     # --production-ps 100000

# covalent adduct + covalent docking + MM-GBSA, each receptor its own leg
python scripts/covalent_workup_one.py --candidate t4_72f5671e89cb \
       --receptors 3IKD 6VAJ --gpu 7

# ADMET, developability, retrosynthesis
python scripts/medchem_workup.py --candidate t4_72f5671e89cb
```

### 4.6 Report

```bash
python scripts/elevation_report.py \
       --movie <trjconv multi-model pdb> --out elevation_report.html
```

Recomputes every figure and table from the shard CSVs at build time rather than
transcribing them, so prose and data cannot silently diverge.

---

## 5. Reading each tier honestly

**Tier 1 — |Δd|, warhead-to-SG displacement over 300 ps.** Smaller is more
stable. Report the median with IQR and the per-replica spread. Do **not** read
the direction of a non-significant effect; n = 8 supports only large effects, and
"≈" means *not distinguished*, not *shown to be equal*.

**Tier 2 — BPMD.** Every replica escapes within 3 ns. **That is expected and is
not the readout** — well-tempered metadynamics is designed to force escape. The
readout is the *cost*: `bias_at_exit_kj` and `frac_in_window`. A replica
reporting `bias_at_exit = 0.000` did not resist at all; the warhead left before
the bias accumulated anything, which means the pose was not in a well.

Tier 2 starts from the **post-equilibration** pose, so tier 1 qualifies tier 2:
a molecule that drifted 5 Å during equilibration is not having its docked pose
tested at all.

**Tier 3 — long MD.** One dissociation event gives a residence estimate with
~100% relative standard error. Quote it as *one draw*, never as a residence time,
and never rank two molecules on single trajectories.

**Never quote a mean across a bimodal trajectory.** The lead's summary CSV
reported `explicit_ligand_rmsd_nm_mean = 1.525`, which averages a bound state and
a dissociated one and describes neither. Split at the transition. The
`explicit_rmsd_suspect` flag exists to catch exactly this.

---

## 6. Failure modes this run hit, so the next one does not

1. **Unit drift between gmx tools.** `gmx rms`, `mindist` and `numcont` write the
   xvg time column in **ns**; `gmx distance` writes **ps** — same trajectory,
   nothing in the file says which. Dividing everything by 1000 squashed three
   series into the first 0.1 ns and produced a plot panel that looked *empty*
   rather than *wrong*. `elevation_report.to_ns()` now anchors on the known
   production length instead of guessing.

2. **Residue-numbering offsets.** 3IKD is renumbered −50 from UniProt, so Cys113
   is residue 63. Verify every label against the **residue identity** at that
   position; an offset off by one produces labels that look entirely plausible.
   Note that AMBER writes histidine as HID/HIE/HIP by protonation state, so a
   naive `== "HIS"` check reports false mismatches.

3. **Receptor inheritance.** `mmgbsa.RECEPTOR_PDB` still defaults to **6VAJ**,
   and every covalent path in the repo takes that default. 6VAJ and 3IKD place
   the pocket 48.6 Å apart. Pass the receptor explicitly and label every score
   with the receptor that produced it; never pool or average the two legs.

4. **Stereocentre inheritance.** Re-embedding an adduct from SMILES can draw the
   *opposite* configuration at the warhead centre. On this molecule the
   diastereomer is worth **1.35 kcal/mol on 3IKD and 2.02 on 6VAJ** — larger than
   most of the gaps the ranking is built on. Carry the configuration from the
   ranked pose; do not re-embed.

Every one of these is the project's signature defect in a new costume: **a value
taken by position, name, label, or inheritance rather than by identity.**

---

## 7. What the suite produced this time

Full result in `docs/elevation_results.md`; decision in **D0071**.

37 molecules (A 8, B 8, D 8, V 5, REF 8), 3 replicas per tier, **111 tier-1 runs,
0 failures**.

| group | n | median \|Δd\| (nm) |
|---|---:|---:|
| A · high enrichment, high consensus | 8 | 0.277 |
| B · low enrichment, high consensus | 8 | 0.198 |
| D · low enrichment, low consensus | 8 | 0.204 |
| V · chloroacetamide, high consensus | 5 | 0.203 |
| **REF · crystallographic** | 8 | **0.102** |

The three pre-registered contrasts: **nothing** (Holm p 0.39–0.88). The anchor:
**everything** (δ −0.59 to −0.78, p 0.007–0.050). Tier 2 independently
reproduced the same null on the completed BDHI groups (p 0.13–0.88 across both
occupancy and bias cost).

**Conclusion: neither enrichment nor consensus predicts pose stability, and the
assay that shows this is demonstrably working.**

---

## 8. Worked example — `t4_72f5671e89cb`, and why it is a NO GO

```
O=S1(=O)CC[C@@H](N(Cc2cnoc2)C2CC(Br)=NO2)C1
C11H14BrN3O4S · MW 364.2 · QED 0.796 · cLogP 1.12 · TPSA 85.0 · SAscore 4.50
warhead: 3-bromo-4,5-dihydroisoxazole (BDHI) on sulfopin's sulfolane core
```

It entered elevation as the **top of the T_4 list**: enrichment 6.86, consensus
1.000 (a single binding mode), 10 of 10 near-attack poses at 3.04 Å and 5.9° off
perpendicular, and genuinely drug-like at 364 Da.

### What each tier returned

| tier | result |
|---|---|
| **0 · literature** | no close analogue; BDHI has thin clinical precedent; the one Pin1 BDHI precedent (Byun 2023) has no structure |
| **0 · med chem** | QED 0.796, SAscore 4.50 — **passes**; the stereocentre is a liability, worth 1.35–2.02 kcal/mol |
| **1 · equilibration** | \|Δd\| 0.350 / 0.762 / 0.508 nm, mean **0.540** — **rank 37 of 37, the worst molecule in the cohort** (REF median 0.102) |
| **2 · BPMD** | frac-in-window 0.081 / 0.041 / 0.150; bias-at-exit 0.010 / **0.000** / 0.464 kJ. Middling *within* group A (rank 3–4 of 8) — but group A is itself the least stable group |
| **3 · 100 ns MD** | stays in the pocket **54.45 ns**, then dissociates and does not return. Warhead is inside the near-attack window for **7.5%** of the bound phase, sitting 8–9 Å out the rest of the time |
| **4 · covalent** | gnina covalent on 3IKD **−7.08 kcal/mol**; adduct topology builds and verifies (GAFF2 `c2`, all junction terms present). **Covalent MD not run** |

### The verdict, and the honest version of the reasoning

**NO GO.** Not because any single number is disqualifying, but because the case
for it rested on a claim the physics does not support.

- **Its rank came from a metric that predicts nothing.** D0071 shows enrichment
  and consensus do not predict pose stability. This molecule was ranked first by
  both.
- **It is the cohort's worst on the one readout that does separate real binders.**
  Rank 37 of 37 on tier 1, against an anchor that separates at p = 0.007. Two
  independently built systems agree — the 100 ns run reproduces the drift
  (0.620 nm by the start of production).
- **The near-attack conformation is visited, not inhabited.** "10 of 10 poses at
  3.04 Å" is a true statement about the docking output and a misleading one about
  the molecule: 7.5% occupancy, concentrated in the 14 ns immediately before it
  left.
- **Its warhead class has no validation whatsoever.** The 17 crystallographic
  Cys113 structures are chloroacetamide (10), naphthoquinone (4) and SNAr (2) —
  **zero BDHI**. Unlike chloroacetamide (AUC 0.908), the BDHI criterion has never
  been tested against a single positive.

### What is *not* the reason

Three things that look like grounds for rejection and are not, recorded so the
no-go is not later remembered as stronger than it was:

- **Dissociating at 54 ns is not a failure.** A covalent inhibitor does not need
  its warhead parked on the sulfur; it needs to reach the window, and the bond is
  permanent once formed. 54 ns is a *long* non-covalent residence for a 364 Da
  fragment.
- **7.5% NAC occupancy is not disqualifying either.** On chemical timescales,
  reaching the window 7.5% of the time across countless rebinding events is
  ample. What it disqualifies is reading pose geometry as occupancy.
- **Tier 2 does not single it out.** It sits 3rd–4th of 8 within its own group.
  The tier-1 result and the missing class validation carry this decision; tier 2
  does not.

### The gap that would change the answer

**There is no baseline for the 54 ns.** We have never run 100 ns on a
crystallographic positive on this receptor, so we do not know whether 54 ns is
excellent or unremarkable. If a REF molecule sits for 500 ns, this result is
weak; if a decoy leaves in 5, it is strong. **The next elevation should put the
anchor through tier 3, not just tiers 1–2** — that is the single cheapest thing
that would make tier-3 numbers interpretable, and it is the same lesson §1
already learned at tier 1.

---

## 9. Standing rules for the next elevation

1. Cohort, never a single molecule. Crystallographic positives in every batch.
2. Pre-register the readings before the first job launches. Report the row that
   fires, including when it is the unwelcome one.
3. Run the anchor through **every** tier you intend to quote, tier 3 included.
4. Fair use lives in the launcher: named GPUs, `nice -n 19`, four cards.
5. Label every number with the receptor, the replica count, and the trajectory
   length that produced it.
6. Split bimodal trajectories before summarising them.
7. Prefer a verified identity over an inherited default — receptor, residue
   number, stereocentre, atom index. Every defect in §6 was that one mistake.
