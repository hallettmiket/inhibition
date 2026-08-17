# Handover delta

*Written 2026-08-17 for @mhallet's return. Compares the project against
[`state_of_the_project.md`](state_of_the_project.md) as it stood on main at the
handover (5b7b682, 2026-08-02). 253 commits since the fork.*

The short version: **the central finding survived the pivot.** §2 of the
handover document was *"we have ~72,000 candidates and no validated way to rank
any of them"*, over four failed levels of theory. The receptor changed, the
ranking was rebuilt around geometry instead of affinity, and the new ranking
measures as unpredictive too — so the finding now stands against **six** methods
rather than four, on two independent endpoints.

Most of the numbered roadmap was **superseded rather than executed**. Phases 2.1
and 2.2 exist to repair the affinity scorer; the work instead removed affinity
from the ranking entirely. That was deliberate and documented, but it means the
plan is not a fair checklist against what happened, and this document separates
the two.

---

## 1. The roadmap, item by item

From §5 of the handover document. Verified against code and outputs on disk.

| # | Step | Status | Evidence |
|---|---|---|---|
| 1 | Phase 0.3a — curate covalent PDB ligands | **done** | Completed before handover. 3 chemotypes vs a floor of 6; unchanged. |
| 2 | Phase 1 — ingest measured-inactive datasets | **done** | AID 504891: **34 actives, 361,354 measured inactives**, 22,578 inconclusive. Exactly the set specified. |
| 3 | Phase 2.1 — run docking through those gates | **done, underpowered** | Run on **3IKD, not 6VAJ**. 22 actives vs 257 inactives: **AUC 0.624**, *p* = 0.054, EF1% 4.23. |
| 4 | Phase 2.2 — pharmacophore as orthogonal scorer | **not started** | No `psearch`/`pmapper`/Pharmit code in the repo. |
| 5 | Phase 0.3c — ensemble docking | **built, not decisive** | `prepare_ensemble_receptors`, `dock_ensemble_shortlist`, `merge_ensemble_dg` exist. Overtaken by the single-receptor change, which addressed the same problem. |
| 6 | Purchasable-analogue mapping | **not started** | No SmallWorld/Arthor integration. |

**On §9 — "if you only do one thing".** Phase 2.1 was run, against measured
inactives, on the new receptor: **AUC 0.624 at *p* = 0.054 over 22 actives**.
That neither rescues docking nor cleanly confirms D0046's prediction. It is the
underpowered middle, and the active count is the binding constraint, not the
method.

---

## 2. Where the work left the roadmap

Four changes. The first two happened within 72 hours of the handover.

**2026-08-05 · D0059 — receptor 6VAJ → chemist-prepared 3IKD.** 6VAJ is
co-crystallised with sulfopin, so its pocket is induced-fit around that ligand.
Re-running D0046's benchmark on 3IKD: best-of-9 pose recovery **15.9% → 41.5%**,
top-1 6.1% → 18.3%. The right pose is in the ensemble 2.6× more often. **Every
6VAJ measurement was invalidated** — D0016, D0041, D0046, D0031, D0049.

**2026-08-05 · #14 — ranking rebuilt around geometry, not affinity.** The
question became *can the molecule orient to form the bond*, not how good the
bond is. A mechanism-specific near-attack criterion — SN2 needs backside
approach, Michael addition perpendicular — replaced the affinity score, which
was removed from the ranking entirely. This is what supersedes Phases 2.1/2.2.

**2026-08-12 · D0081 — scope narrowed to T_4, then to three warhead families.**
T_3 is REINVENT output and 96% acrylamide, so "acrylamide dominates the library"
was a statement about which generator ran. Ranking is T_4 only; scope is
acrylamide, bdhi_c4 and bdhi_c5 — the chemistry the Lu lab will synthesise.

**2026-08-16 · D0085 — a four-stage cascade with measured gates.** Docking +
near-attack rank → 8 ns triage MD → 100 ns production → BPMD. The 8 ns length
and the 0.35 nm survivor bar came from a documented experiment over a 0.1 ns
grid ([`sweep_length.md`](sweep_length.md)), not from a guess.

---

## 3. Then and now

| | 2026-08-02 handover | 2026-08-17 |
|---|---|---|
| receptor | 6VAJ | **3IKD**, chemist-prepared |
| arms in play | four (T_1–T_4) | **T_4 only** |
| molecules | ~72,000 docked and ranked | **561** in scope, all screened |
| ranking basis | affinity, size-decorrelated | **near-attack geometry**, empirical-Bayes shrunk |
| unit of selection | the molecule | **the binding mode** |
| pose handling | one representative per molecule | two-stage splitting, ≤5 sub-modes |
| MD | ad hoc, not reproducible | 3-stage cascade with measured gates |
| decision records | 57 | 85+ |
| releases | 1.0.0 | **3.1.0** — four since |
| running compute | nothing | pipeline live on 3 GPUs |

---

## 4. The finding that matters most

| method | result | record |
|---|---|---|
| Docking enrichment | AUC 0.599, EF1% 0.0 | D0041 |
| Docking pose recovery | 5% in production | D0046 |
| Ensemble MM-GBSA | below chance | D0036 |
| MD residence | not reproducible | D0038, D0044 |
| Contact-profile fit score | worse than chance; built and killed | D0057 |
| **Near-attack geometry ranking** | **ρ = +0.119, *p* = 0.33** vs sweep outcome | this run |

Measured on the 68 modes swept to date. `enrichment` gives ρ = +0.033.
`class_rank` gives **ρ = −0.256 at *p* = 0.035** — the only one that clears
significance, and it points the *wrong way*: better-ranked modes were less
stable. Median enrichment among survivors is **6.03**; among those that left,
**6.12**. The single highest-enrichment mode in the library (12.25) left the
pocket.

**The caveat that must travel with this:** every one of the 68 cleared the
enrichment floor, so this measures discrimination *above* the floor, not whether
the floor itself excludes anything. The stratified pilot that would test the
floor (#71) has never been run.

This replicates 3.0.0's ρ = +0.016 with four times the power and against a
*different endpoint* — that was attack-readiness, this is stability. Two
independent outcomes, neither predicted.

**The cascade's first full result.** The best mode either run has produced —
`t4_071099f4034c_m1`, 0.323 nm and 84.8% attack-ready over 8 ns — came back from
100 ns at **1.140 nm**. It moved 3.5× further and missed the bar. One data
point, but it is the best candidate in hand.

---

## 5. What needs a decision

| Question | Why it is open | Issue |
|---|---|---|
| **Is the enrichment floor real?** | It discards 3,700 of 4,432 modes on a parameter with no demonstrated relationship to outcome. The pilot is ~28 GPU-h, ~9 h on 3 cards. | [#71](https://github.com/hallettmiket/inhibition/issues/71) |
| **Should the sweep run best-ranked-first?** | If rank carries no information, best-first buys nothing and costs an unbiased sample if the campaign stops early. | — |
| **Does BPMD get a card?** | Stage four has no GPU. It consumed **1.3 TB** in 3.0.0 — 87% of the project's footprint — for no usable result, and nothing prunes it. | [#72](https://github.com/hallettmiket/inhibition/issues/72) |
| **Does naphthoquinone come back?** | It was in the original three-family ask and dropped when scope narrowed. Two complete 198-molecule arms sit unused; ~13 min of docking to add. | — |

---

## 6. Honest gaps in this comparison

**Phase 2.1's result is not a clean test of the D0046 prediction.** It ran on
3IKD against a different active set, with 22 actives. Reading it as confirming
or refuting D0046 would be over-reading.

**The 100 ns result is n = 1.** Five more survivors are queued; the conclusion
above is a calibration, not a rate.

**The old-vs-new ranking comparison is confounded.** `t4_716800c125a7`, 3.0.0's
headline molecule, ranks 89th of 774 in its class now — but the clustering was
rebuilt, so its old mode 4 and its current mode 2 are not the same object. It
has not yet been swept in this run.

**The §8 warning about discipline-not-enforcement was right, and it cost us.**
Cross-run contamination — pages and queries reading a previous screen's
directories — occurred **eleven times** in the last week, and five were caught
by @tt8804 rather than by a test. There is now a guard,
`tests/test_runs_are_standalone.py`, that fails on any pipeline module naming a
shared run directory, a hardcoded topic, or a dataset root. 103 files outside
the pipeline path remain unguarded (#74).
