# Branch `receptor/3ikd-chemist-prepared`

*Opened 2026-08-05 by @tt8804. Chemist's recommendation: use only the prepared
3IKD structure, not 6VAJ.*

This branch exists so the receptor change can be tried without disturbing what
is already measured on 6VAJ. **Nothing here is merged until the benchmarks are
re-run** — see "what switching invalidates" below.

---

## Why

6VAJ is Pin1 co-crystallised with sulfopin, so the pocket is in an **induced-fit
state shaped around that ligand**. Docking into it biases toward sulfopin-like
chemistry, which is the criticism raised in #14.

There is now independent evidence for the concern from our own data. Building
the contact-profile fit score (D0057), the same test run two ways gave:

| | crystal pose ranked #1 | top 3 |
|---|---|---|
| cross-docked into 6VAJ | 0 / 82 (0.0%) | 3.7% |
| self-docked into each ligand's own receptor | 5 / 82 (6.1%) | 22.0% |

A ligand's crystal pose transplanted into 6VAJ does not make 6VAJ-like contacts.
That is a receptor-transfer penalty measured on 82 structures, and it is a
separate line of evidence from the chemist's recommendation.

**The honest caveat, which the branch must not lose:** 3IKD is *also* a
ligand-bound structure ("Structure-Based Design of Novel PIN1 Inhibitors (I)",
2.0 Å, cognate ligand **J9Z**). Swapping 6VAJ for 3IKD trades sulfopin-induced
fit for J9Z-induced fit *unless the chemist's preparation specifically addresses
that*. "Not induced-fit" and "induced-fit by a different ligand" are different
claims and the prepared file should state which it is.

---

## Status

- [ ] **File not yet on the server.** It is at `/home/tt/Downloads/3ikd_well_prepared.pdb`
      on @tt8804's own machine. Transfer:
      `scp /home/tt/Downloads/3ikd_well_prepared.pdb twu383@129.100.24.200:/data/lab_vm/immutable/inhibition/receptor/`
- [ ] **Provenance unrecorded.** A chemist-modified PDB is not a deposited
      structure. Needed before it becomes load-bearing: what changed relative to
      deposited 3IKD (ligand stripped? waters? protonation and at what pH?
      rotamers or loops rebuilt?), and whether **J9Z is still present** — the
      box has to be derived from the entry's own reference ligand, because
      6VAJ's box is a set of coordinates in 6VAJ's frame and means nothing here.
- [ ] Verify it parses, Cys113 SG present and a **reactive thiol** (or T_3/T_4
      have nothing to attack).
- [ ] Run through `shared/receptor_prep.py`'s exact path (strip → `reduce
      -BUILD` → `obabel -xr`) so the receptor is the only thing that differs.
- [ ] Derive the box; record SHA-256 + provenance as a decision record.

## What switching invalidates

Every gate verdict and benchmark on record is a **6VAJ number**. Re-run before
any of them is quoted about 3IKD:

| measurement | record | why it does not transfer |
|---|---|---|
| non-covalent enrichment ROC-AUC 0.599 | D0016 / D0041 | measured on 6VAJ with `box_expanded.json` |
| **pose recovery 5% production / 55% best-of-9** | **D0046** | the number this whole redesign rests on |
| covalent enrichment at chance | D0031 | class-matched decoys docked into 6VAJ |
| size decorrelation ρ ≤ 0.034 | D0049 | residuals of a 6VAJ metric |

`attach_gate` looks a verdict up by (stratum, metric). The metric NAME does not
change when the receptor does, so an old verdict would attach silently to the
new receptor. D0051 makes an *unknown* metric fail closed — this is a *known*
metric on a different structure, which is worse. **Re-run the gate before
anything ranks.**

The receptor ensemble machinery already landed on `main` (`Receptor` registry,
per-receptor boxes, receptor-tagged pose directories), so 3IKD can be a
first-class receptor rather than a hard-coded swap.

---

## Parked work, to return to

### In flight when this branch opened

- **guo_pfizer degree-2 enumeration** — running in tmux `degree2`, 8,670
  parents, reservoir target 30,000, local `/dev/shm` fragment DB. atra is done
  (D2_7: 30,000 kept of a measured population of 4,063,427, which reproduces the
  earlier run exactly). potter_astex, du_xu and liu_2024_c3 were **dropped by
  @tt8804 to save time**.
- **ATRA degree-2 docking** finished (exit 0) with some failures to account for
  — "the frame carries what succeeded; re-run to fill the rest."

### The ranking method worked out in #14, not yet built

Design agreed with @tt8804 2026-08-05. Non-covalent only, free form, no covalent
docking at this stage: the question is **whether the molecule can orient to form
the bond, not how good the bond is.**

1. **The criterion is mechanism-specific, not a distance cutoff.**
   `warhead_classes_10.csv` carries `reactive_atom_smarts` and `mechanism`, and
   the approaches differ: `sn2_displacement` needs backside attack anti to the
   leaving group (S···C–LG ≈ 180°); `michael_addition` needs approach
   perpendicular to the alkene plane; `snar_displacement` perpendicular to the
   ring. A chloroacetamide 3.5 Å from the sulfur with its chlorine pointing at
   it is chemically dead and passes any distance filter.
2. **Binary gate then continuous rank.** Can it reach a near-attack conformation
   (NAC), reproducibly? Then among those that can, how stably is the warhead
   held?
3. **Funnel:** free NAC screen on existing poses → 5-seed replicate docking
   (~37 GPU-h for T_3+T_4) → BPMD only on survivors.
4. **Bias along d(reactive atom, SG), not whole-ligand RMSD.** A molecule whose
   scaffold drifts while the warhead stays locked is a *good* answer; standard
   BPMD would penalise it. A 1-D CV also converges faster, so the protocol can
   be shorter than Clark's 10 × 10 ns.
5. **Validation before candidates.** 17 ligands verified covalent at Cys113 are
   positives; the 6 of 16 warhead classes that collapse to an inert amide are
   built-in negatives. If the method cannot separate them, nothing downstream is
   worth running.
6. **Scope:** applies natively to T_3 + T_4 (7,178 molecules). T_1 and T_2 have
   no warhead — report as *not applicable*, not ranked last.

**Measured costs** (A100, a real system from this project):

| | rate |
|---|---|
| plain MD | 695.9 ns/day |
| under PLUMED bias | **485.8 ns/day** (30% overhead) |
| a week on 6 free GPUs | ~20,400 ns ≈ **200–1,000 poses** depending on protocol |

**PLUMED 2.10 is installed and working** at
`~/.micromamba/envs/dwi_plumed`; GROMACS 2026.3 reports `Plumed support:
enabled` and loads the kernel via `PLUMED_KERNEL`, verified end to end on a real
system. Proper home is `/data/lab_vm/envs/dwi_plumed`, which needs someone with
write access there to create.

### Open, unchanged by this branch

**#13** open technical problems · **#12** chemistry judgement out to the Lu lab ·
**#14** this framework · **#15** enumeration efficiency.
