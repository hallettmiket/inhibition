# State of the project

*Written 2026-07-31 at handover. Start here.*

This is the orientation document, not the README. The README says how to run
things. This says **what we are trying to find out, what we have established,
what is still open, and what I would do next.**

Read alongside:

* [`how_this_project_breaks.md`](how_this_project_breaks.md) — the one pattern
  behind every bug found here. **Read it second, before writing any code.**
* `decisions/` — 45 records. They document what was decided *and what was wrong
  and why it looked right*. They are the most valuable thing in the repo.
* GitHub **#4** (the plan) and **#6** (open decisions). Only two issues are open,
  deliberately.

---

## 1. What this project is

`inhibition` is a murmurent **choreography**: one problem attacked by four
independent approaches, with an integration layer that presents their
shortlists for a human to adjudicate.

The problem is finding an inhibitor of human **Pin1**, catalytic **Cys113**,
against **PDB 6VAJ**.

**The deliverable is the method, not the molecule.** Pin1 is the testbed. This
matters when prioritising: a result about whether the choreography works beats
a result about any individual compound. (One consequence: whether Pin1 is a
good oncology target is *not* this project's concern — that was briefly listed
as a blocking question and withdrawn.)

| | approach | seed | search | mechanism |
|---|---|---|---|---|
| **T_1** | de novo, structure-based | none — the pocket | DiffSBDD | non-covalent |
| **T_2** | derivative neighbourhood | 5 seeds (was ATRA only) | CReM | non-covalent |
| **T_3** | R-group decoration | sulfopin' | REINVENT LibInvent | covalent Cys113 |
| **T_4** | warhead × R-group | sulfopin core | combinatorial | covalent Cys113 |

**13,863 candidates** across the four arms, plus **42,588** more from the T_2
reseeding now docking, plus a ~30,000-molecule degree-2 sample generating.

---

## 2. The one thing to understand

> **We have ~56,000 candidates and no validated way to rank any of them.**

That is not pessimism; it is the measured position, and it is the project's
central finding so far. Four levels of theory have been tested and none
discriminates:

| level | result | record |
|---|---|---|
| Docking enrichment | AUC 0.599, CI [0.311, 0.874], **EF1% 0.0** | D0041 |
| **Docking pose recovery** | **5% in production**, 16% self-dock, vs 60–80% norm | **D0046** |
| Ensemble MM-GBSA | below chance | D0036 |
| Implicit + explicit MD residence | not reproducible | D0038, D0044 |

**D0046 is the newest and the most decisive**, because it explains the others.
Docking recovers a known Pin1 pose 5% of the time on the receptor T_1 and T_2
actually use. And the failure is **scoring, not sampling** — best-of-9 is 55%
while top-1 is 22.5%, so the search *finds* the crystallographic pose and the
scoring function then ranks it below a wrong one three times in four.

So the enrichment null is not a decoy-construction artefact. The same function
fails a strictly easier question than enrichment.

**This yields a live prediction.** Phase 2.1 of #4 re-runs Vina against gates
built from *experimentally measured* inactives. If decoys were the problem,
that rescues it. On D0046's evidence it should not. That experiment is now a
test of an explanation rather than an open question — run it.

---

## 3. What is established

Things I would defend in review.

* **Docking does not work on this pocket** — enrichment (D0041) *and* pose
  recovery (D0046), by independent measurements.
* **Ranking is partly a size sort.** ρ = −0.617 (T_1), −0.479 (T_3) between
  rank metric and heavy-atom count. Published elsewhere too, so not a local
  bug. Ligand efficiency over-corrects (ρ = −0.938) and is not the fix. D0043.
* **The covalent stratum is UNDERPOWERED and stays that way.** Counting
  chemotypes by warhead class gives **4** against a floor of 6. Structural
  clustering would have given exactly 6 — the decision was fixed *before* the
  counts were taken, precisely so the answer could not be chosen. D0045.
* **Explicit solvent stabilises poses relative to implicit** — T_1 −0.290 nm,
  T_2 −0.469 nm mean ligand-RMSD change, every T_2 candidate improving. The
  implicit-solvent dissociations were substantially a solvent-model artefact.
* **T_3's "protac-like" molecules come from the scorer, not the generator.**
  LibInvent runs once. Generated median 25 heavy atoms; *shortlist* median 39.
  D0043.
* **Synthesizability needs structural rules, not SAscore.** Seven rules, each a
  named impossibility. On a 5-molecule retrosynthesis test SAscore ordered
  correctly but its fragment-additive form structurally cannot represent "the
  combination is unprecedented" — which is how a molecule with an acyl
  phosphate scored *easiest* in T_1's top 10.
* **The four arms differ in how synthesizable-by-construction they are.** T_4's
  route is known (it is enumerated from a coupling); T_1 carries no synthetic
  information at all. Confirmed empirically: AiZynthFinder solved Sulfopin (4
  steps) and T_4's top (2 steps) and failed T_1, T_2 and T_3's tops.

---

## 4. What is ruled out — do not redo these

From the literature sweeps and the measurements. Full reasoning in #4.

* **Swapping DiffSBDD for a better SBDD model** (TargetDiff, DecompDiff,
  MolCRAFT, PocketFlow…). All are trained on CrossDocked and ranked by Vina —
  the function D0046 shows recovers poses 5% of the time here. And D0043 makes
  it worse than neutral: a model generating larger molecules scores better on
  our own ranking without binding better.
* **Boltz-2 as a gate scorer.** Its affinity head trained on ChEMBL/BindingDB/
  PDBbind, and every one of our actives came from ChEMBL2288. That is leakage,
  not signal.
* **Alchemical FEP.** Blinded benchmark 2.44 kcal/mol RMSE; our pocket is
  shallow/polar/water-mediated (its worst case), our ligands are charged, and
  they share no common core — forcing ABFE at ~52 GPU-h/compound.
* **ML rescoring** (RTMScore, KarmaDock, PIGNet2…). PDBbind-trained,
  DUD-E-validated; the critique literature says the wins are memorisation.
* **Genetic algorithms / REINVENT RL against the current objective.** They
  optimise the oracle harder, and the oracle is broken. They would amplify
  D0043 fast and it would look like progress.
* **More literature searching for a 6th covalent chemotype.** The field
  converged on chloroacetamide. The missing chemotypes are in *screening data*
  and *structures*, not papers.

---

## 5. What to do next

In order. Detail in **#4**.

1. **Phase 0.3a — curate the covalent PDB ligands.** The single highest-value
   task. Verify the covalent partner via `struct_conn` (RCSB says "covalent"
   but not to *what* — could be Cys57), map warheads via
   `warhead_library.canonical_class()`. The survey holds ≥4 chemistries the
   reference set lacks — **aryl aldehyde (11 structures), maleate ester (4),
   SuFEx (2), Mannich (1)** — taking the chemotype count 4 → 8 and unblocking
   the covalent gate legitimately. 1–2 days, no compute.
2. **Phase 1 — ingest the two measured-inactive datasets.** Dubiella's
   993-fragment screen (111 hits / 882 measured non-hits) and PubChem AID
   504891 (34 actives / **361,392 measured inactives**, verified). These give
   gates whose negatives were *assayed* rather than assumed.
3. **Phase 2.1 — run Vina through those gates.** Now a test of D0046's
   prediction.
4. **Phase 2.2 — pharmacophore as an orthogonal, gate-testable scorer.**
   `psearch`/`pmapper` (same lab as CReM) or Pharmit. Measure it on **pose
   recovery first** — D0046 gives a cheaper, stricter harness than enrichment.
5. **Phase 0.3c — ensemble docking.** If cross-docking into an ensemble of the
   163 X-ray receptors beats 5%, rigid-receptor error dominates and the
   ensemble is the fix.
6. **Purchasable-analogue mapping** (SmallWorld/Arthor, free) — converts an
   unrankable shortlist into orderable matter regardless of whether ranking is
   ever fixed.

---

## 6. Open decisions — #6

Eight items, each with question, evidence and a recommendation. Two are
**defects with a decision attached**, not preferences:

* **`affinity_kcal` is the CNN-best pose's affinity**, not the best affinity.
  89% of covalent candidates affected; **>50% of each covalent shortlist would
  change**. No re-docking needed.
* **T_1's structural alerts are computed and never acted on** —
  `alert_gate_pass` is `True` for all 4,803 rows.

The rest: adopt two new synthesizability rules (`acyl_phosphate`,
`stereogenic_phosphorus` — both kill 0 known binders); the N-hydroxylamine
rule; T_2 phosphate protect-vs-label; receptor ensemble-vs-fixed;
charge-stratified ranking; and whether there is synthesis/assay capacity for a
real covalent fragment screen.

---

## 7. Running compute, and what it needs

| | state |
|---|---|
| T_2 five-seed docking, GPUs 4 & 7 | liu + du_xu running; potter, guo, atra queued. **Chain ends ~06:00** |
| Degree-2 ATRA sample | ~900/1,882 parents, ~30,000 molecules; then needs docking |
| Explicit MD | **done** — 243/245 replicates, merged into D1_21/D2_21 |
| Redocking benchmark | **done** — D0046 |

All detached (`PPID 1`), survive SSH disconnects. **GPUs are shared** — three
other users are on this box; cut back if it is busy.

**One known defect in the merged MD frames:** `D1_21`/`D2_21` carry
`explicit_rmsd_replicate_sd_x`/`_y` and no canonical column, because
`merge_gromacs_results.py` drops a hand-maintained column list that omits the
aggregates it builds. The fix is written out in
`docs/session_state_2026-07-31.md`; the script is **untouched**, so the tree
is consistent.

---

## 8. Things that will bite you

* **Environments live outside the repo** — `/data/lab_vm/envs/dwi_{cheminf,gui,
  amber_md,gromacs_cuda,vinagpu,gnina,reinvent4,diffsbdd,admet,retro}`. Clone
  and run will not work without them.
* **`immutable/` is read-only, `append_only/` is append-only**, enforced by
  hooks, not convention. Frames are integer-versioned; retire superseded ones
  in `data/ready_to_delete.md`.
* **Streamlit does not re-import helper modules.** Editing `curate.py` and
  clicking Rerun gives you the old module. There is a guard that stops the page
  and says so; restart the process.
* **Reference files resolve by glob** (`reference_set.latest_reference`). Do
  not re-pin a version literal — a test walks the AST and will fail you.
* **Read `how_this_project_breaks.md` before trusting any number.**

---

## 9. If you only do one thing

Run **Phase 2.1**. It is cheap, it tests D0046's prediction, and either outcome
is publishable: a scorer that works against measured inactives, or a second
independent demonstration that docking does not work on this pocket.

The honest framing of this project's contribution, as it stands today, is not a
Pin1 inhibitor. It is a **worked demonstration that a multi-agent choreography
can measure its own methods honestly enough to report that they do not work** —
D0016 → D0031 → D0041 → D0043 → D0045 → D0046 is a coherent negative-result
methods paper, and very few groups measure their own null at their own power
floor and publish the number.
