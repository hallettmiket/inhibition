# 2.2.0 “Chalcopyrite” — build plan

*What gets implemented, why, what could go wrong with each, and how each one is
accepted or rejected. Written 2026-08-07, after the measurements in #23/#29/#30
changed what is possible. Companion to [`outline_2.2.0.md`](outline_2.2.0.md),
which stated the intent; this states the build.*

---

## 0. Read this first — the outline's logistics premise is dead

`outline_2.2.0.md` §4 says:

> **Steps 1–3 need no new simulation.** That is the dividend from 2.1.0
> persisting its working.

**That is not true, and the whole order of work depends on it.** `nac_screen_v2`
persists per-pose rows for the **top 20 by energy only** — verified: every
molecule has exactly 20 rows, `energy_rank` 1–20. The other 180 poses are
discarded at write time and no longer exist.

What survives per molecule:

| | persisted? | consequence |
|---|---|---|
| aggregate counts over all 200 (`n_in_range`, `n_viable_given_in_range`) | **yes** | `enrichment_conditional` **is** computable from existing data |
| per-pose geometry (distance, angle, viable) for all 200 | **no** — 20 only | `anchor_quality` over the full population is **not** computable |
| pose coordinates for all 200 | **no** — 20 only | **pose splitting is not computable.** Clustering 20 poses selected *by energy* into modes clusters the survivors of the exact filter we are removing |

So the central feature of 2.2.0 requires re-docking the library. That is not a
disaster — measured from the v2 run's own frames, ~7.8 s/molecule/shard, so 5,765
molecules on 4 shards is **≈3 GPU-hours** — but it reorders everything and it
must be done **once**, carrying every change that needs a re-dock.

---

## 1. The build

### 1.1 Persist the population, not the survivors *(new — forced by the above)*

**Change.** `nac_screen_v2` stops writing "the 20 lowest-energy poses" and starts
writing a pose set that can be split. Per-pose geometry for **all 200**;
coordinates for a mode-diverse subset.

**Considerations.**

- **Storage is the reason `KEEP_TOP` exists, and it is not actually binding.**
  Per-pose geometry for 200 poses × 5,765 molecules is ~1.1 M rows — trivial.
  Coordinates are the expensive part (~30 heavy atoms × 200 × 5,765 × 4 bytes ≈
  1.4 GB as float32), which is affordable but should be a `.npy` per molecule
  rather than a column in a CSV.
- **Write coordinates for a mode-diverse subset, geometry for everything.** The
  geometry table is what scoring and clustering need; full coordinates are needed
  only for the poses that might get elevated.
- **This lands in `append_only/`**, so the old frames stay. That is the rule, and
  it means the change is non-destructive by construction.

**Accepted when:** a re-run of the 15 ground-truth molecules reproduces #30's
93.3% containment *from the persisted files alone*, with no docking.

---

### 1.2 Pose splitting — the core feature

**Change.** Each molecule stops having *a* score and gets a **ranked set of
binding-mode hypotheses**, each scored on its own.

```
molecule ─┬─ mode A   n=64   anchor 0.71   viable 0.55   ← elevate
          ├─ mode B   n=51   anchor 0.12   viable 0.02
          └─ mode C   n=18   anchor 0.44   viable 0.30   ← and maybe
```

**Design constraints, each from something measured:**

| constraint | why |
|---|---|
| cluster on **reactive-atom position + approach vector**, not whole-molecule RMSD | D0062 — two poses that place the warhead identically and differ in a distal ring are one mode for a covalent question |
| clustering **must not use docking energy** | #23/#30 — energy places the correct pose at a rank indistinguishable from uniform (KS p = 0.666). Using it to define or order modes re-imports the defect this version exists to remove |
| **mode count is measured, not a parameter** | fixing *k* decides the answer in advance. Density/agreement criterion, count reported with its stability across re-runs |
| a mode below a **minimum population is labelled noise**, not dropped | same rule the consensus floor follows |
| modes get **identities that survive a re-dock** — geometric, not rank-based | a mode named by rank is a mode a re-run silently redefines |

**Considerations that are not in the outline:**

- **Docking is stochastic, and mode identity has to survive that.** Re-running
  the 15 cases reproduced containment exactly (93.3% both times) but moved
  `KEEP_TOP`-survival from 46.7% → 40.0% and top-10 from 33.3% → 26.7%. If mode
  *membership* is that unstable, mode-level scores inherit the instability. The
  stability check is therefore not optional garnish — it is the acceptance test.
- **The 150° SN2 threshold interacts with clustering.** Sulfopin's crystal mode
  peaks at 146.4°, so a hard cut calls the crystallographic mode of a nanomolar
  inhibitor dead by 3.6°. If viability is a *per-mode* property, one threshold
  decision now silently kills whole modes. Prefer a graded angular term inside
  the mode score over a binary gate on its members.
- **What happens when a molecule has one mode?** Most will. The output shape must
  not privilege multi-mode molecules, or ranking acquires a bias toward
  promiscuous binders — the exact failure mode consensus had, inverted.

**Accepted when:** (a) mode assignments are stable across two independent
re-docks of the same molecule; (b) on the 15 ground-truth cases, the crystal pose
falls inside a named mode ≥90% of the time; (c) elevating the top mode beats
33.3% (the current top-10 number), with 93.3% as the ceiling.

---

### 1.3 The score — conditional enrichment and anchor quality, per mode

**Change.** Primary score becomes **conditional enrichment** — P(angle viable |
distance in window) — computed **per mode**. `anchor_quality` recomputed over a
mode's full population rather than averaged over the top-20 by energy.

**Considerations.**

- **`enrichment_conditional` is the one thing computable today**, from persisted
  aggregate counts. It can and should be tested before the re-dock.
- **The current `WEIGHTS = {anchor_quality: 0.5, topn_viable_frac: 0.5}` is a
  placeholder, not a finding.** #30 showed `topn_viable_frac` reads a window that
  holds the crystal pose 26.7–33.3% of the time. **Half the current weight rests
  on a mostly-wrong window.** Re-derive after the re-dock; do not carry the 0.5/0.5
  forward as if it meant something.
- **Test against Sulfopin first, not last.** Any candidate score that gives the
  parent compound a zero is wrong. Its crystal pose is produced at 1.00 Å at
  energy rank 160/200 — so a correct score *can* find it, and there is no excuse
  for a repeat.
- **Convergence check before it ranks anything** — 200 vs 2,000 runs, per D0068.
  No score in this project has passed one yet.

**Accepted when:** Sulfopin scores non-zero; the score converges between 200 and
2,000 runs; and it is not rank-correlated with docking energy (if it is, it has
re-imported the defect).

---

### 1.4 PoseBusters as a validity gate

**Change.** Installed and unused today. Gate poses **before** clustering.

**Consideration — this ordering is the whole point.** A reproducibly-invalid pose
is *reproducible*, so it will form a confident, well-populated cluster and look
exactly like a real mode. Gating after clustering would let physical nonsense
become a named hypothesis with a good score.

**Accepted when:** the fraction gated is reported per warhead class, and no class
loses so much that its ranking becomes unpopulated.

---

### 1.5 D0074 — the reactive path protonates at pH 7.4

**Change.** `nac_screen.prepare_ligand` adopts `obabel -p 7.4`, matching the
non-covalent path. **Accepted 2026-08-07; not yet implemented.**

**Considerations.**

- **594 of 1,782 T₄ (33.3%) and 331 of 5,370 T₃ are currently docked as the wrong
  species.** Every score for those describes a molecule that does not exist at
  pH 7.4. This is a correctness blocker, not a refinement.
- **It rides the same re-dock as pose splitting.** Separately, the library gets
  docked twice.
- **`charge_ph74` becomes the single charge annotation** once both paths agree —
  the reason for the divergence disappears with it.
- **The T₃ bidirectional disagreement is a separate defect** (142 rows where
  `charge_ph74` is *lower* than the SMILES charge, which protonation cannot
  cause). Not fixed by this and must not be assumed fixed.
- **Two 100 ns runs launched with `--net-charge +1`** (`t4_da2e98512d02`,
  `t4_9a973be6b946`) simulate the unprotonated species. Internally consistent,
  but the species this decision rules against. They must be redone.

---

### 1.6 Boltz-2 — adopt as a cross-check, hold on triage

**Measured (pre-registered in `prereg_cofolding.md`, cutoff verified 2023-06-01):**

| | n | within 2 Å | median |
|---|---:|---:|---:|
| **held out** | 10 | **60.0%** | 1.93 Å |
| in training (control) | 5 | 80.0% | 1.36 Å |
| **our docking, same 10** | 10 | **10.0%** | 2.82 Å |

**Change.** Adopt as an **independent cross-check on every elevated pose** — it
never sees our poses, so it is the only orthogonal signal we have.

**Considerations.**

- **It cannot model the covalent bond.** Every prediction is the *non-covalent*
  pre-reaction complex. It can never speak to whether a molecule reacts. Hard
  limit.
- **Its affinity head is largely pose-independent** (established in 2.1.0). Use
  structure and confidence; never affinity.
- **The MSA cache is the premise, not an optimisation.** One MSA per construct
  makes a 300-molecule triage a forward pass each. Without it the cost argument
  collapses.
- **The triage use (T4) is unproven** — it needs the bornite runs. Do not build
  selection on it yet.

---

### 1.7 Carried-forward fixes

| fix | why it matters now |
|---|---|
| **`mmgbsa.RECEPTOR_PDB` required, not defaulted** | still defaults to **6VAJ**, which is 48.6 Å from the 3IKD pocket, and every covalent path takes it |
| **A chain stage that produces no output stops the run** | the overnight chain logged `exit 1` twice and continued, turning a good ranking into zero elevation |
| **Flexible Cys113 sidechain** | the anchor distance is measured to one arbitrary rotamer of the residue the criterion is about |
| **Covalent MD** | topology built and verified since 2.0.0, never run |
| **Manual elevate button** (#22) | selection is automatic-only; no way to push a molecule onto the queue by hand |
| **Mode-aware GUI** | once a molecule has modes, the viewer and tables must show them, or the feature is invisible to the person using it |

---

## 2. Revised order of work

The outline's order assumed no re-dock. This one does not.

| # | step | GPU | blocks |
|---|---|---|---|
| 1 | **Design + validate the keep rule** on the 15 ground-truth cases (200 poses each, already dumped) | none | everything |
| 2 | **`enrichment_conditional` + Sulfopin test + convergence**, on existing aggregates | none | — |
| 3 | **PoseBusters gate** wired in | none | 4 |
| 4 | **ONE re-dock**: D0074 protonation + new persistence + new keep rule | ~3 h | 5–7 |
| 5 | **Pose splitting** on the new pose sets; stability check | light | 6 |
| 6 | **Re-rank per mode**; re-derive weights | none | 7 |
| 7 | **Re-select and re-elevate** | heavy | — |
| — | *in parallel:* MD-priority verdict, chain robustness, receptor default, GUI | — | — |

**Steps 1–3 genuinely need no GPU**, which is the real dividend — and step 1 is
already running.

---

## 3. What I need from you

1. **Does the re-dock go now or after the MD verdict?** It is ~3 GPU-hours and
   the bornite runs are using 5 cards. My inclination is after, so we are not
   competing for GPUs and so the verdict is not delayed.
2. **The 150° SN2 threshold.** A chemist's call, and it now has teeth: it decides
   whether whole modes are viable, and it currently calls Sulfopin's
   crystallographic mode dead by 3.6°.
3. **N-activated acrylamides — 97% of T₃** (D0066), still no ruling. If they are
   not credible, T₃ is 97% noise and the re-dock should skip it.
4. **FEP licensing** — a governance decision, not a technical one (see
   `outline_2.2.0.md` §3.1, corrected).

---

## 4. What would make 2.2.0 a failure

Stated in advance:

- **A score that ranks well and still gives Sulfopin a zero.**
- **Mode counts that change between re-runs** — then a "mode" is an artefact of
  the clustering, not a property of the ligand.
- **Elevating a mode BPMD then shows was the wrong one, repeatedly** — would mean
  modes are ranked on something as uninformative as energy was.
- **Another silent stage.** If a 2.2.0 run can produce zero output and report
  success, nothing else here matters.
- **Shipping the re-dock without D0074** — a third of T₄ would stay the wrong
  molecule, and every downstream number would inherit it.
