# The 2.1.0 “Bornite” framework

*Settled 2026-08-06 with @tt8804. This is the architecture; `CHANGELOG.md` will
carry the release when the version number is fixed (see §8 — whether this is
2.1.0 or 3.0.0 depends on a decision not yet made).*

The one-line version: **2.0.0 ranked molecules on a score nobody had tested and
elevated whichever pose docking happened to score best. 2.1.0 separates those
into two ranking problems, scores geometry continuously instead of as a
threshold, and makes every stage's inputs verified rather than inherited.**

---

## 1. The four stages

```
1  chemical space          T_3 (4,082, one warhead) · T_4 (1,683, 9 classes × 187)
        │
2  RANK MOLECULES          consensus filter WITHIN class → weighted anchoring score
        │                  → one ranked list per warhead class
3  SELECT                  automatic (in order, geometry re-measured) | manual (GUI)
        │
4  ELEVATE                 4a  RANK POSES within the molecule, by BPMD
                           4b  100 ns non-covalent MD on the winning pose
                           4c  gate: ≥90% residence AND majority in attack geometry
                           4d  covalent MD 50 ns + covalent docking score
                           → candidates list
```

**The separation at stages 2 and 4a is the structural change.** They are different
questions and were being answered by the same machinery:

| | question | evidence used |
|---|---|---|
| **stage 2** | which molecules deserve the compute | consensus + anchoring geometry, over a pose window |
| **stage 4a** | which of *this molecule's* poses is real enough to spend 4 GPU-hours on | BPMD occupancy per pose |

Stage 2 asks whether a molecule *can* present its warhead. Stage 4a asks whether
a specific pose *is physically there*. Docking energy answers neither well — on
this receptor AutoDock's own ordering puts the right pose first **18.3%** of the
time — so neither stage takes it on trust.

---

## 2. Ranking molecules (stage 2)

### The filter: consensus, as a quota within warhead class

Pose agreement among the top 10. Applied as a **fixed fraction within each
warhead class**, not as a single library-wide bar.

**Why** (D0073): a single 0.90 cut passed BDHI at 16.6% and chloroacetamide at
2.9%, because pose agreement is partly a statement about how many ways a molecule
*can* sit. Pass rate is monotone in rotatable bonds (19.5% at 0–2 down to 2.9% at
7+), and warhead classes differ systematically in flexibility. A library-wide bar
is therefore a rigidity ranking wearing a geometry label.

**What it does not fix.** The *within*-class part of the confound survives: T_3 is
100% acrylamide and consensus still correlates with rotatable bonds at
**ρ = −0.312** there. Rigidity is reported per class rather than claimed handled.

A **consensus floor** sits under the quota. A class whose 20% cut falls below
~0.5 has no well-determined poses, and the quota would otherwise pass molecules
whose poses simply disagree — `sulfonate_acetamide` cut at 0.24 and is held to 5
survivors instead of 37.

### The score: weighted, anchoring-based

```
weighted_score = w₁ · anchor_quality  +  w₂ · topn_viable_frac
```

**`anchor_quality`** is enrichment's question asked *continuously*. Enrichment
counts poses that clear a window; that throws away most of what was measured — a
pose at 3.5 Å and 2° off ideal and one at 4.19 Å and 29.9° both score 1, and 4.21 Å
scores 0. The window is a decision boundary, not a measurement.

So distance and angle each get a factor that is 1.0 at the ideal and decays, and
they **multiply**. A pose at perfect distance and hopeless angle is not
half-good; a sum would let one term hide the other.

**`topn_viable_frac`** is the fraction of the top 10 poses *by energy* that are
reaction-competent — the metric D0068 argued for and 2.0.0 never implemented. It
cannot be diluted by more searching, because it is defined on the molecule's own
best poses rather than as a per-run rate.

**Weights are equal and unvalidated on purpose.** Nothing has yet shown which
component predicts anything, so a tuned weight would be fitted to nothing.
`WEIGHTS` is one edit and the sensitivity sweep is the next experiment.

### Output shape

**One ranked list per warhead class. No global top-N.** A merged ordering
re-imports the bias the per-class quota exists to remove, so the frame does not
carry a column to sort on.

### What was removed

**The validation gate.** Earlier drafts restricted the shortlist to classes with
crystallographic positives. @tt8804: the 15 depositions are too few and too poor
to decide which chemistry to pursue, and the chloroacetamide series they mostly
come from is not of interest. Classes are ranked on their own terms.

The crystallographic set keeps a *different* job — calibrating the elevation
stability assay, where 8 molecules make one comparison readable. Choosing
chemistry and calibrating an instrument are not the same use.

---

## 3. Selection (stage 3)

**Automatic** (default): walk the ranked list in order, taking the top *n* per
class.

**The re-check measures the pose, not the score.** The ranking score is an
aggregate over a pose window; what gets elevated is one specific pose, and they
disagree often. Measured: two of the top three T_3 molecules had a rank-1 pose
outside attack geometry, one **8.03 Å** from the sulfur.

A molecule whose poses all fail is **queued with `geometry_ok = False` and a
reason**, not dropped — "the ranking liked it and none of its poses is
reaction-competent" is a finding about the ranking.

**Manual**: an elevate button in the GUI, recording who queued it and why, so a
queued molecule traces to a decision rather than looking like the pipeline chose
it.

---

## 4. Elevation (stage 4)

### 4a. Rank the molecule's own poses, by BPMD

Every reaction-competent pose gets BPMD; the most stable wins.

**Readout is occupancy (`frac_in_window`), not escape cost.** The completed
elevation run measured `bias_at_exit` separating nothing (all p ≥ 0.08) while
tracking occupancy at ρ = 0.974 — the escape-cost term is nearly inert at 3 ns.
Occupancy is what separated crystallographic positives from candidates
(p = 0.005–0.021).

**Stability decides among reaction-competent poses, it does not override them.**
Collapsing the two criteria would let a stable but badly-angled pose win.

Cost ~1 GPU-hour, protecting a 4 GPU-hour run. Worth it: `t3_0ac2aa9133ef` has
three viable poses at 12.0°, 4.4° and 22.5° off ideal, and nothing before this
stage chose between them on evidence.

### 4b–4d. The trajectory, the gate, the covalent leg

| | |
|---|---|
| **4b** | 100 ns unbiased explicit-solvent MD from the winning pose |
| **4c** | gate: **≥90% residence** AND **majority of time in attack geometry** |
| **4d** | covalent MD 50 ns + covalent docking score → candidates list |

**The residence half of the gate is a floor, not a stretch target, and needs no
calibration.** Residence time is 1/k_off, and k_off = K_D × k_on; at a
conservative k_on of 10⁷ M⁻¹s⁻¹ even a *millimolar* fragment should sit for
~100 µs. Dissociating inside 100 ns implies K_D ≈ 1 M — which is another way of
saying the molecule does not bind. (This is why `t4_72f5671e89cb` leaving at
54.45 ns is a strong statement, not a near miss.)

**The attack-geometry half is deliberately ambitious.** Reaction rate scales with
time spent in attack geometry, so demanding a majority demands a preorganised
molecule rather than one that merely visits.

Caveats kept on the record: force fields over-dissociate, so a pass is more
trustworthy than a fail; one trajectory is not a rate, so a single-replicate
failure is a flag rather than a verdict; periodic-boundary effects contribute.
Together those argue for 3 replicates on any molecule that matters, and for
treating 1 replicate as a screen.

**The covalent leg is not chained automatically.** Firing 50 ns off a
single-replicate measurement carrying ~100% relative standard error spends
compute on a coin flip.

---

## 5. What every stage stands on

### Pose generation

**AutoDock-GPU for sampling, gnina CNN for re-ranking.** Measured on 3IKD_ian
over 82 Pin1 crystal ligands, sampling held fixed:

| ordering | top-1 ≤2 Å |
|---|---:|
| sampling finds a good pose (ceiling) | 41.5% |
| AutoDock's ordering | 18.3% |
| **gnina CNNscore** | **26.8%** |

23.2 points of pure ranking failure; gnina closes 37% of it at zero sampling
cost. And **when a good pose exists it is in the top 10 in 100% of cases** —
which is why every score here is computed over a pose window rather than on
rank 1.

Rejected on evidence: cross-program consensus rescoring (improves consistency,
does not beat the best single scorer), DL docking as the pose source (fails
physical validity), co-folding as a geometry source (Boltz-2's affinity is
largely pose-independent).

### The receptor

**3IKD_ian**, resolved by pinned SHA-256 and refusing on mismatch. There are two
structures on this box answering to "3IKD" — the chemist's (1,807 atoms, 7
waters) and a deposited RCSB copy (2,741 atoms, 484 waters) — and the
non-covalent docking path was reading the deposited one.

### Data that is now kept

The 2.0.0 screen computed per-pose geometry, reduced it to one number, and
discarded the poses. 2.1.0 persists per-pose distance/angle/viability, gnina
scores on every retained pose, and **the poses themselves**, so a score can be
repaired without docking again.

### Reference molecules (#19)

The known binders go through the identical criterion, with warhead class assigned
**by SMARTS match, not by the prose in the file** (only 2 of 22 rows matched a
canonical class as a string). Sulfopin — the parent — scores **2.54**; ZL-Pin13
(IC50 67 nM) scores **3.43**.

They are on the screen so a chemist can see what the incumbent scores on the same
measurement. **They are not a validation set and nothing is calibrated against
them.**

---

## 6. The recurring defect this version is built against

Every significant error in 2.0.0, and several caught while building 2.1.0, is one
mistake: **a value taken by position, name, label, or inheritance rather than by
identity.**

| caught | taken wrongly |
|---|---|
| D0067 | mechanism *name* over hybridisation — 374 candidates read as dead |
| pose atoms | list *position* over chemical identity |
| SNAr deletion | a type *string* over the fact of retyping |
| GUI poses | an *inherited* default receptor |
| adduct | connectivity over stereochemistry — 1.35–2.02 kcal/mol |
| PLUMED | absence of an error over a positive check |
| Cys113 SG | a residue *number* that changes between pipeline stages |
| `nac.measure` | argument *position* — called as (sg, coords, mechanism) |
| reference classes | a prose *label* over a SMARTS match |
| BPMD resume | an incomplete *key* — ident+replicate, ignoring trajectory length |

The standing rule: **verify by identity, and make the check fail loudly.**

---

## 7. Known-unfixed, carried forward

- `mmgbsa.RECEPTOR_PDB` still defaults to **6VAJ**, and every covalent path takes
  that default.
- Covalent MD has never run; the topology is built and verified.
- No crystallographic positive has been through the 100 ns tier — not needed for
  the gate (§4c), but it means no tier-3 number has a peer.
- The accessible-decoy-cysteine control is still inconclusive (Cys57 is buried).
- The **within-class** rigidity confound is measured and unaddressed.
- Chemist ruling outstanding on N-activated acrylamides — 97% of T_3 (D0066).
- Weights in the composite are equal and untested.

---

## 8. What is not settled

- **The version number.** If the rework redefines `enrichment` rather than adding
  beside it, old and new values share a name and are not comparable, which takes
  a major bump. See `docs/versioning.md`.
- **Whether any of these scores predicts stability.** The elevation cohort is the
  instrument; the test is to recompute the new scores on its existing poses and
  run the pre-registered contrasts. Cheap, and it can falsify — though the cohort
  was *selected on* the old metrics, so it cannot confirm.
- **The quota and floor values** (20%, 0.5). Chosen to be sensible, not measured.
- **Whether the residual mechanism effect survives rigidity.** SNAr and Michael
  share a median rotatable-bond count and pass at 13.9% vs 6.3%.

---

## 9. Where the code is

| stage | script |
|---|---|
| screen, persisting everything | `scripts/nac_screen_v2.py` |
| rank molecules | `scripts/rank_v2.py` |
| GUI frames | `scripts/gui_refresh_v2.py` |
| select | `scripts/select_elevate.py` |
| rank poses (BPMD) | `scripts/rank_poses_bpmd.py` |
| elevate | `scripts/elevate_queue.py` |
| references | `scripts/screen_references.py` |
| re-scorer benchmark | `scripts/rescore_benchmark.py` |
| the whole chain | `scripts/overnight.sh` |
| receptor identity | `shared/receptors.py` |

Design rationale: `docs/ranking_2.1.0_design.md`. Elevation protocol:
`docs/elevation_example.md`. What 2.0.0 established and disproved:
`docs/recap_2.0.0.md`.
