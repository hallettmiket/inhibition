# 2.2.0 “Chalcopyrite” — pose splitting, and a score that survives its own test

*Opened 2026-08-07 at @tt8804's direction: "tooling upgrades and an upgraded
pose-splitting feature". Companion to `docs/retrospective_2.1.0.md`.*

**The thesis: 2.1.0 asked "how good is this molecule's pose?" of a molecule that
does not have *a* pose. It has a distribution of them, often several distinct
binding modes, and every score so far has flattened that into one number —
either by averaging over an arbitrary window or by measuring how much the poses
agree, which penalises exactly the molecules that have a real second mode.**

---

## 1. Pose splitting

### 1.1 What is wrong now

The current pipeline treats a molecule's 200 poses as one population and reduces
it to scalars: a viable fraction, a consensus, an average anchor quality. Three
consequences, all measured in 2.1.0:

1. **A molecule with two genuine binding modes is punished.** Consensus asks
   whether the top poses agree. Two well-defined modes score as disagreement,
   which is indistinguishable from a molecule whose poses are simply scattered.
   `consensus_modes` already records that the two cases differ; nothing acts on it.
2. **Averaging mixes modes that should never be averaged.** `anchor_quality` is a
   mean over the window. A molecule with one excellent reaction-competent mode
   and one inert mode gets a mediocre mean describing neither.
3. **42.7% of molecules have viable poses that the scoring window never sees**
   (#23). If those viable poses form a coherent mode, that mode is a real
   hypothesis being discarded because its members score poorly on energy.

### 1.2 The change

**Split each molecule's poses into distinct binding modes, then score every mode
on its own and carry them forward separately.**

```
molecule ─┬─ mode A   n=64   anchor 0.71   viable 0.55   ← elevate this
          ├─ mode B   n=51   anchor 0.12   viable 0.02
          └─ mode C   n=18   anchor 0.44   viable 0.30   ← and maybe this
```

A molecule stops having *a* score and starts having a **ranked set of
hypotheses**, which is what a docking run actually produces.

### 1.3 Design constraints, each from something already measured

- **Cluster on the reactive-atom position and approach vector, not whole-molecule
  RMSD.** D0062 established that whole-molecule RMSD is the wrong endpoint for a
  covalent question. Two poses that place the warhead identically and differ in a
  distal ring are the same mode for our purposes.
- **Clustering must not use docking energy.** #23: energy carries no signal about
  reaction geometry, so using it to define or order modes re-imports the defect
  2.2.0 exists to remove. Energy may be *reported* per mode; it must not shape
  them.
- **Mode count is a measured property, not a parameter.** Fixing *k* would decide
  the answer in advance. Use a density- or agreement-based criterion and report
  the count with its stability across re-runs, the way D0068 forces any number to
  carry its defining parameter.
- **A mode below a minimum population is noise, and is labelled as such** rather
  than dropped — the same rule the consensus floor already follows.
- **Modes get identities that survive a re-dock.** A mode named by its rank is a
  mode a re-run silently redefines; name it by its geometry (reactive-atom
  centroid + approach direction) so the same mode is recognisable across runs.

### 1.4 What it unlocks downstream

- **Selection elevates a *mode*, not a molecule.** The pose handed to MD is the
  representative of a named mode, and the record says which.
- **Elevation can test two modes of one molecule** and find out which is real —
  currently impossible to even express.
- **Consensus is replaced by something honest**: not "do the poses agree" but
  "how many modes, how populated, how reaction-competent is the best one".
- **#23's population becomes reachable.** Viable poses outside the energy window
  form modes that can be scored on their merits.

---

## 2. The score (issue #23)

Pose splitting does not by itself fix the score. Two changes, in order:

1. **Primary score moves to conditional enrichment** — P(angle viable | distance
   in window) — computed **per mode**. It uses no ordering, so it cannot inherit
   a bad one, and it is a ratio within the in-range subset so it should not
   dilute the way the per-run rate did in D0068.
2. **`anchor_quality` recomputed over the full population, per mode**, not
   averaged over the top-20 by energy.

**Both must pass the convergence check before they rank anything.** 200 vs 2,000
runs, the check D0068 forces and that no score in this project has yet passed.
Pre-register it, as with the elevation experiment.

**And test against Sulfopin first.** Any candidate score that gives the parent
compound a zero is wrong, and that check costs an afternoon. Not doing it in
2.1.0 is the single most avoidable thing in the retrospective.

---

## 3. Tooling — the four in issue #26, assessed

Availability checked first, because a tool we cannot run is not a decision.
**Neither Desmond nor Boltz-2/AF3 is installed here.** AlphaFold3 *is* reachable
in this lab: @emucaki runs it containerised on A100s for `sh2_domains`
(`/app/alphafold`, `/alphafold3_venv` in the job logs), so weights and a working
image exist under DeepMind's terms of use.

### 3.1 Desmond — recommend against, and the reason is your own FEP plan

| | |
|---|---|
| **Licence** | free for academics, **but the academic licence excludes free-energy calculations** |
| **Force field** | fixed OPLS in the academic build |
| **Our stack** | GAFF2/AM1-BCC + ff19SB, a hand-built covalent junction (`cys_gaff2_junction_5.frcmod`), PLUMED-driven BPMD |

Three problems, in descending order:

1. **It blocks the thing you want it for.** FEP is a free-energy calculation and
   the free academic licence does not permit them. Desmond + FEP means
   Schrödinger FEP+ commercially. So "upgrade MD to Desmond" and "do FEP before
   synthesis" do not compose on an academic licence — that has to be a purchase
   decision, not a tooling one.
2. **The covalent junction would have to be rebuilt in OPLS.** D0030 and the
   `junction_5` work exist because parameterising a Cys–warhead bond is
   failure-prone; redoing it in a force field we do not control, to gain what
   GROMACS already does, is cost without a measured benefit.
3. **BPMD is PLUMED-on-GROMACS.** Pose ranking — the stage that just proved it
   changes which pose gets elevated — would need reimplementing.

**What Desmond is actually good at is analysis and presentation**, which is the
one place we already solved the problem (`shared/report_theme.py`). If the appeal
is the reporting rather than the engine, that is worth saying out loud, because it
is a much cheaper thing to want.

**Recommendation:** stay on GROMACS. Revisit only if a licence covering FEP+ is
bought, in which case Desmond comes with it and the question answers itself.

### 3.2 Boltz-2 / AlphaFold3 — recommend testing, for a reason neither of us stated

The interesting property is not accuracy. It is that **Boltz-2 never sees our
docking poses.**

Every signal we have tested — docking energy, enrichment, consensus,
top-N viability, anchor quality — is computed from the same AutoDock pose set, so
they inherit the same failure modes and correlate with each other for reasons
that have nothing to do with the molecule (#23). A co-folding model builds the
complex from sequence and ligand alone. Whatever it says is **orthogonal**, and
orthogonal is precisely what a 300-molecule triage is missing.

Two distinct uses, worth separating:

**(a) Pose recheck / MD preparation** — what #26 asks for. Both models emit
per-prediction confidence (pLDDT, interface PAE, ipTM) usable as a filter, and
prospectively they reproduce >50% of novel-ligand poses under 2 Å. But note the
caveat from the 2.1.0 work: **Boltz-2's affinity head is largely pose-independent**,
so it is a weaker instrument for pose-level questions than its headline numbers
suggest. Use the *structure* and its confidence, not the affinity, for anything
about a pose.

**(b) MD-priority triage** — not in #26, and possibly the bigger prize. Interface
PAE is a **forward pass, seconds per complex**, against ~1 GPU-hour for BPMD. If
it predicts 100 ns residence even moderately, it triages 300+ molecules for
roughly the cost of one BPMD run.

**Boltz-2 also bears on the FEP question directly:** it is reported to approach
FEP accuracy while running ~1000× faster. That does not replace FEP before
synthesis, but it plausibly *triages* which molecules deserve one.

### 3.3 The test, which is cheap and uses molecules whose answers we know

Same design as `docs/prereg_md_priority.md`. Run Boltz-2 (or the lab's AF3
container) on the nine molecules whose 100 ns outcome we either have or will have
within a day:

- **Sulfopin** and **ATRA** — both held 100 ns
- **`t4_72f5671e89cb`** — left at 54 ns
- **the six MD-priority molecules** — outcomes landing now

Then ask whether interface confidence separates held from left. Nine points, four
of them with known answers today. **If it does, it is the cheapest filter
available by three orders of magnitude. If it does not, we have spent an
afternoon.**

Two things to get right, both learned the hard way here:

- **Predict before looking.** Pre-register as with the BPMD test.
- **The covalent question is not asked.** Co-folding models do not model the
  Cys113 bond, so this can only ever be about whether the molecule sits in the
  pocket — never about whether it reacts.

### 3.4 gnina and FEP — settled and deferred

**gnina** is adopted and running: CNN rescoring lifted top-1 pose accuracy from
18.3% to 26.8% on this receptor, and it scores every retained pose in the screen.
Standing filter, no further decision needed.

**FEP** stays where #26 puts it — an optional terminal step before synthesis, on
a handful of molecules. The licence question in §3.1 decides *how*, and the
Boltz-2 result may decide *which*.

### 3.5 Carried-forward fixes, unchanged

| | why |
|---|---|
| **PoseBusters as a validity gate** | installed and unused. A reproducibly-invalid pose is *reproducible*, so it survives mode-splitting as a confident cluster. Must gate **before** clustering |
| **`mmgbsa.RECEPTOR_PDB` required, not defaulted** | still defaults to **6VAJ** and every covalent path takes it |
| **Chain fails on a failed stage** | the overnight run logged `exit 1` twice and carried on, producing zero elevation from a good ranking |
| **`charge_ph74` vs `canonical_smiles`** | disagree in 33% of T₄ (#28). Feeds antechamber `-nc`, so a wrong value parameterises the wrong species silently |
| **Flexible Cys113 sidechain** | the anchor distance is measured to one arbitrary rotamer of the residue the criterion is about |
| **Covalent MD** | topology built and verified since 2.0.0; never run |

## 4. Order of work

1. **Score first, on existing data.** Conditional enrichment and full-population
   anchor quality are computable from the persisted v2 poses — no docking. Test
   against Sulfopin and the Reddi compounds, then convergence.
2. **Pose splitting**, on the same persisted poses. Also no re-dock.
3. **PoseBusters gate**, before splitting is trusted.
4. **Re-rank and re-elevate** once 1–3 hold.
5. Chain robustness and the receptor default alongside, since they are small.

**Steps 1–3 need no new simulation.** That is the dividend from 2.1.0 persisting
its working, and it is why the version was worth running even though its score
did not survive.

---

## 5. What would make 2.2.0 a failure

Stated in advance, in the spirit of the elevation pre-registration:

- **A score that ranks well and still gives Sulfopin a zero.** The first test, not
  the last.
- **Mode counts that change between re-runs of the same molecule.** Then a "mode"
  is an artefact of the clustering, not a property of the ligand.
- **Elevating a mode that BPMD then shows was the wrong one, repeatedly.** Would
  mean modes are being ranked on something as uninformative as energy was.
- **Another silent stage.** If a 2.2.0 run can produce zero output and report
  success, nothing else in this document matters.

---

## 6. Open questions carried in

- Is the **SN2 150° threshold** too strict? Sulfopin clears it 34 times in 200
  poses but never among its best-scoring. A chemist's question.
- Does the **secondary pocket** Reddi 2023 reports next to the active site
  ([10.1021/jacs.2c08853](https://doi.org/10.1021/jacs.2c08853)) give a second
  anchor worth defining? It is exactly the generalisation issue #17 asks for, on
  a real series with measured potency.
- **N-activated acrylamides** — 97% of T₃ (D0066), still no chemist ruling.
- The **within-class rigidity confound** (ρ = −0.312 inside T₃ alone) is measured
  and unaddressed.
