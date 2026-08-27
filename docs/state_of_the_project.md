# State of the project

*Written 2026-07-31 at handover. Last updated 2026-08-19. Start here.*

> **The measured numbers in §7 are generated now** — `scripts/refresh_orientation.py`,
> and `tests/test_orientation_current.py` fails the suite if they drift (D0055,
> closing #11). It caught 800 molecules of drift on its first run.
> **The prose is still yours to keep true.** If you change what this describes,
> change this: it drifted badly within 24 h of being written and a new
> maintainer read it as fact.

This is the orientation document, not the README. The README says how to run
things. This says **what we are trying to find out, what we have established,
what is still open, and what I would do next.**

Read alongside:

* [`how_this_project_breaks.md`](how_this_project_breaks.md) — the one pattern
  behind every bug found here. **Read it second, before writing any code.**
* `decisions/` — they document what was decided *and what was wrong and why it
  looked right*. They are the most valuable thing in the repo. (Count in §7,
  generated.)
* GitHub, consolidated 2026-08-04 — **#12** (the chemistry judgement we cannot
  supply computationally, out to the Lu lab) and **#13** (every open technical
  problem, audited against the code rather than against the threads). **#4**
  remains the plan and reasoning of record. #2, #6 and #8 are closed into those
  two; read them for history, not for status.

---

## 0. Where 3.1.0 landed, and the thing it uncovered

*2026-08-19. Read this before §3 — it changes how much of §3 you should trust.*

**The run finished.** 561 molecules screened, 4,432 modes ranked, 147 triaged at
8 ns, 15 given 100 ns runs. Five held the pocket, one held unstably, nine left.
The GUI reports it at `http://localhost:8931/`, ranked on mean ligand RMSD with
max shown beside it, four tiers (optimal / held / held-unstable / left).

**And the modes it ranked are not modes.** @tt8804, looking at one: *"why would
only 20 out of 82 poses be good?? that suggests that they arent the same
poses."* They are not. Measured across nac_v5: the median mode spans **3.51 A**
in warhead-to-anchor distance, 87% span more than 2 A, and **42% have a viable
fraction between 0.1 and 0.9** — two populations under one label. The largest is
137 poses spanning 9.3 A.

**The cause is an ordering mistake, and it is circular** (D0088). The pipeline
clusters on the reactive atom's POSITION and the direction its warhead faces,
then scores each group by how often it reaches attack geometry. But attack
geometry *is* position and direction — SG is fixed within a run. It forms groups
along the axis it then grades them on. The code comment claiming otherwise
(*"never on the NAC geometry itself, which is the score"*) is wrong.

**So `viable_fraction`, `enrichment` and `conditional_eb` are all measured over
mixtures.** Every per-mode number in nac_v5 carries that caveat. The 100 ns
results do not — a trajectory is a trajectory — but *which* pose earned each one
was chosen by this machinery.

**Three things were fixed outright** and are on `main`-track code:

| | |
|---|---|
| the screen was not reproducible | `docking.seed` now set; it was seeded from the clock, so v4 and v5 ranked the same 504 molecules at rho = +0.43 (#77) |
| the pose cloud could not be joined to its own measurements | `pose_idx` is now written, and the cloud is rewritten with its run rather than cached (#76) |
| AutoDock-GPU failed silently above ~2,000 runs | at 5,000 it corrupts its stack **and still writes a .dlg**; at 10,000 it fails with exit 0. `dock()` now checks all three |

**One replacement is built and NOT adopted**: `shared/pose_cluster.py` — a single
clustering step on pose similarity alone (HDBSCAN over heavy-atom RMSD), with
attack geometry used only to rank afterwards. It is the only rule measured that
never produces a bag (largest mode 14 poses, widest 3.91 A), and it puts the
pose a 100 ns run validated into an **8-pose group 1.5 A wide** where the shipped
rule puts it in **108 poses spanning 8 A**. It stays `proposed` because it
discards 29% of the cloud as noise and lost the validated pose in 3 of 30
replicates (#78), and adopting it means a full re-screen (#79).

**What NOT to do next:** do not run the enrichment-floor pilot (#71) or BPMD
(#72) against nac_v5's modes. Both would calibrate against the artefact. The
sequence is #78, then #79, then re-screen, then those.

**Also settled since §3 was written:** the triage sweep is 5 ns (D0087, inside
D0085's own CI, and truncation is one-sided so it cannot drop a survivor); the
100 ns "optimal" bar is 0.45 nm and separate from the 0.35 nm sweep bar (they
were the same number, which made "optimal" unreachable); and mode count does not
saturate with sampling depth because the density threshold is 5% *of the sample*
(D0088).

---

## 1. What this project is

`inhibition` is a murmurent **choreography**: one problem attacked by four
independent approaches, with an integration layer that presents their
shortlists for a human to adjudicate.

The problem is finding an inhibitor of human **Pin1**, catalytic **Cys113**.

**The receptor is contested, and the code currently disagrees with itself.**
D0059 (2026-08-05, @tt8804, relaying the chemist) replaces **6VAJ** with the
prepared **3IKD**: 6VAJ is co-crystallised with sulfopin, so its pocket is
induced-fit around that ligand, and cross-docking into it ranked the crystal
pose #1 in **0 of 82** cases against 5/82 self-docked. But D0059's status is
still `proposed`, `config/receptor.yaml` still pins `pdb_id: 6VAJ`, and
`shared/noncovalent_dock_run.py:61` still hardcodes `6VAJ_prepared.pdbqt` —
while the benchmark and reference-screen paths call `resolve_3ikd_ian()`, which
refuses to run against the wrong 3IKD. **Which receptor you get depends on which
entry point you came through.** Settle this before running anything; see
[`retrospective_2.2.0.md`](retrospective_2.2.0.md) §3.3.

Measurements below that predate D0059 were made on 6VAJ and are labelled where
it matters — they are not silently reinterpreted as 3IKD results.

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

**13,863 candidates** across the four arms, plus **42,588** from the five-seed
T_2 reseeding and **15,653** from a degree-2 ATRA sample — **~72,000
molecules, all docked and ranked** as of 2026-08-02.

---

## 2. The one thing to understand

> **We have ~72,000 candidates and no validated way to rank any of them.**

That is not pessimism; it is the measured position, and it is the project's
central finding so far. Four levels of theory have been tested and none
discriminates:

| level | result | record |
|---|---|---|
| Docking enrichment | AUC 0.599, CI [0.311, 0.874], **EF1% 0.0** | D0041 |
| **Docking pose recovery** | **5% in production** (6VAJ, invalidated); on 3IKD **18.3% top-1 / 41.5% best-of-9** | **D0046**, #66 |
| Ensemble MM-GBSA | below chance | D0036 |
| Implicit + explicit MD residence | not reproducible | D0038, D0044 |

**D0046's framing has been corrected against the literature (#66), and the
correction matters.** This table used to compare our recovery against a
"60–80% norm". **That is the SELF-docking norm** — redocking a ligand into the
structure it was crystallised in. Our 82-case benchmark docks *non-cognate*
Pin1 ligands into one prepared 3IKD, which is **cross-docking**, where the
published baseline is roughly **41–50% top-1 for a single receptor** and 67–69%
only with best-structure selection. So:

* our **top-1 of 18.3%** is genuinely below the cross-docking norm — by about
  2×, not the ~3× the old framing implied;
* our **best-of-9 of 41.5% is *at* the single-structure cross-docking norm**,
  and the old framing presented it as a failure.

And the failure being **scoring rather than sampling** is not our finding: it is
documented and quantified for the exact program we run. Across docking programs,
sampling produces a near-native pose in 85–99% of cases while scoring ranks it
first in 35–73%, and **AutoDock Vina shows the largest gap of any program
tested — 93.4% sampling against 35–40% ranking**. Our numbers are that
phenomenon at lower absolute rates because the setting is harder. Treat it as a
**positive control** — evidence the measurement apparatus works — not as a
result. See [`publication_audit.md`](publication_audit.md) §2–3.

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
* **Ranking is partly a size sort — and the DIRECTION depends on the pool.**
  ρ = −0.617 (T_1), −0.695 (T_3) between rank metric and heavy-atom count, but
  **positive** in the heavier T_2 pools: +0.205 (liu, mean 45 heavy atoms),
  +0.119 (du_xu), +0.092 (potter). Consistent with Vina's size preference
  saturating and reversing once molecules outgrow the pocket. Ligand efficiency
  is NOT the fix — re-measured 2026-08-01 it is worse than the raw score in
  five of six pools (T_1 −0.938 vs −0.617). D0043, D0049.
* **Ranking is now size-decorrelated** (D0049): the metric's residual against
  heavy-atom count, taken as a local median within equal-population strata.
  Every arm lands at |ρ| ≤ 0.034. This makes the ordering *less wrong in one
  identified way*; it does not make it valid.
* **The covalent stratum is UNDERPOWERED and stays that way — now measured
  twice, independently.** Counting chemotypes by warhead class over the
  reference set's lead-tier actives gives **4** against a floor of 6; structural
  clustering would have given exactly 6, and the definition was fixed *before*
  the counts were taken precisely so the answer could not be chosen (D0045).
  The Pin1 PDB, curated against `struct_conn` on 2026-08-01, independently
  gives **3** — chloroacetamide, naphthoquinone, SNAr chloroazine, over 17
  ligands verified covalent at Cys113. Two different sources, both short of the
  floor, neither rescued by the other.
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
* **More literature searching for a 6th covalent chemotype** — and, as of
  2026-08-01, **the structural route too.** The field converged on
  chloroacetamide, and the PDB has now been curated against `struct_conn`: 17
  ligands verified covalent at Cys113, yielding 3 chemotypes. The four
  chemistries the survey appeared to offer do not survive — the SuFEx entry is
  bonded to TYR23, and aryl aldehyde and maleate ester have no covalently
  linked instance at all. **The missing chemotypes are in *screening data*
  only** (Dubiella's 993-fragment set, #4 Phase 1.1), or must be commissioned
  (#8 C3).

---

## 5. What to do next

In order. Detail in **#4**.

1. ~~**Phase 0.3a — curate the covalent PDB ligands.**~~ **DONE 2026-08-01,
   and it does NOT unblock the covalent gate.** Verified against `struct_conn`:
   31 ligand–cysteine covalent links across 18 entries and 17 distinct ligands,
   **all 31 at Cys113** (the Cys57 worry does not bite). Counted by warhead
   class: chloroacetamide 10, naphthoquinone 4, SNAr chloroazine 2, plus one
   unclassified Mannich aminoketone. **3 chemotypes against a floor of 6 —
   lower than D0045's 4, not higher.** The hoped-for chemistries do not
   survive: the SuFEx entry (8VZ3) is bonded to **TYR23**, not a cysteine, and
   aryl aldehyde and maleate ester have no covalently-linked instance at all
   (9KE5 has no `_struct_conn` section — it is non-covalent). Those counts came
   from entry titles, not modelled linkages. `scripts/curate_covalent_pdb.py`;
   full result in #4.
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

## 6. Open decisions — now #13

*Updated 2026-08-04. #6 is closed; every item below is tracked in **#13**,
audited against the code. Chemistry judgement moved to **#12**.*

**Closed since this was written:**

* `affinity_kcal` was the CNN-best pose's affinity — **fixed**, D0047. Worse
  than first reported: 25% of T_4 candidates were ranked on *positive*
  (clashing) affinities, up to +159.7.
* T_1's structural alerts computed and never acted on — **fixed**.
  `alert_gate_pass` is now NA when no gate ran, with `alert_gate_applied`
  beside it. Report-not-gate, per PI decision.
* Two synthesizability rules — **adopted**, D0048. Both reject zero known
  binders; 236/4,803 fire in T_1, zero elsewhere.
* Size-decorrelated ranking — **decided and implemented**, D0049. Ligand
  efficiency was measured and rejected: worse than the raw score in five of
  six pools.

Also closed since: `rank_validated` validated on a denylist (**D0051** — only
`STRONG` validates, unknown verdicts fail closed as `UNGATED`); the stale-pin
guard generalised to every versioned stem across `config/` too, with the three
stale warhead pins fixed.

**Still open:** the N-hydroxylamine rule and synthesis/assay capacity (both in
**#12**, out to the Lu lab, unanswered since 2026-07-31); charge-stratified
ranking; T_2 phosphate protect-vs-label (decided *label*, not yet implemented);
receptor ensemble — **partially built as of 2026-08-04**: the receptor registry,
per-receptor boxes, receptor-tagged pose directories and the median combination
rule have landed with tests (D0052, D0053); 3IKG/3IKD/9INR still need preparing
and the gate still has to be run on the ensemble metric before it ranks
anything. **Scope decided: shortlists first** (~125 molecules/seed, minutes)
rather than the ~265 GPU-hour full re-dock.

---

## 7. Running compute, and what it needs

**Nothing is running as of 2026-08-02 07:30.** The T_2 campaign is complete.

| | state |
|---|---|
| Explicit MD | **done** — 243/245 replicates, merged and clean from `D1_24`/`D2_24` |
| Redocking benchmark | **done** — D0046 |
| Covalent PDB curation | **done** — #4 Phase 0.3a; see §5 |

### Measured state of every frame

*Generated — do not hand-edit. `python3 scripts/refresh_orientation.py`.
`tests/test_orientation_current.py` fails the suite if this drifts, which is
issue #11 and the reason the hand-maintained version of this table was wrong
within 24 h of being written.*

<!-- AUTO:arms:BEGIN -->
| arm | frame | rows | docked | ranked | shortlist |
|---|---|---:|---:|---:|---:|
| T_1 de novo (DiffSBDD) | `D1_32.parquet` | 4,803 | 3,233 | 3,233 | 25 (`shortlist_synth`) |
| T_3 R-group (LibInvent) | `D3_38.parquet` | 5,396 | 4,080 | 4,080 | 25 (`shortlist_synth`) |
| T_4 warhead x R-group | `D4_52.parquet` | 1,783 | 1,683 | 1,684 | 27 (`shortlist_synth`) |
<!-- AUTO:arms:END -->

<!-- AUTO:t2:BEGIN -->
| T_2 seed | frame | docked | ranked | shortlist |
|---|---|---:|---:|---:|
| ATRA | `D2_33.parquet` | 1,882 | 1,882 | 25 |
| Liu-2024-C3 | `D2_5.parquet` | 16,806 | 16,806 | 25 |
| Potter-Astex | `D2_5.parquet` | 7,376 | 7,376 | 25 |
| Du-Xu | `D2_10.parquet` | 9,736 | 9,736 | 25 |
| Guo-Pfizer | `D2_10.parquet` | 8,670 | 8,670 | 25 |
| ATRA degree-2 | `D2_8.parquet` | 127 | 0 | 0 |
| **all six** | | **44,597** | | |
<!-- AUTO:t2:END -->

<!-- AUTO:decisions:BEGIN -->
**97** decision records.
<!-- AUTO:decisions:END -->

All six T_2 variants are ranked (size-decorrelated, D0049) and carry rebuilt
synthesizable shortlists.

**GPUs are shared** — `ysun2443` and `wzhan564` are also on this box. Check
`nvidia-smi` for other people's processes before taking a card.

**Two operational lessons, both paid for in GPU time:**

*Vina-GPU is all-or-nothing.* It writes every pose at the END of a
virtual-screening run. The liu pool ran 24 h on one card, hit a flat
`timeout=86400`, and was killed with **0 of 16,806 poses written**. The
timeout now scales with the pool (`vina_timeout_s`) and is logged before the
run starts.

*Split large pools.* `scripts/dock_chunked.py` runs one pool as N chunks
across N GPUs into one pose directory. The same liu pool then finished in
**7.4 h across five cards, 0 chunks failed** — and a failure now costs one
chunk rather than everything.

---

## 8. Things that will bite you

* **Environments live outside the repo** — `/data/lab_vm/envs/dwi_{cheminf,gui,
  amber_md,gromacs_cuda,vinagpu,gnina,reinvent4,diffsbdd,admet,retro}`. Clone
  and run will not work without them.
* **`immutable/` and `append_only/` are a DISCIPLINE, not enforcement.** This
  said "enforced by hooks, not convention", which is misleading in a way worth
  being exact about. Verified 2026-08-02: both trees are **writable at the
  filesystem level** — `test -w` succeeds, and @tt8804 created and removed a
  directory inside `immutable/`. The guarantee comes from
  `~/.claude/hooks/block-rm.sh`, a **per-user Claude Code hook**: it does
  nothing in a plain shell, nothing for a script run outside a CC session, and
  nothing for a different user until they install their own. Treat read-only as
  a property you maintain, not one the system maintains for you. Frames are
  integer-versioned; retire superseded ones in `data/ready_to_delete.md`.
* **`modifiable/` is the scratch tree** — logs, launch drivers, chunk
  directories. Nothing there may be cited by a manifest, decision record or
  frame, and it is the only place under `/data/lab_vm` where deletion is
  allowed. See its own README.
* **Analysis artefacts go to `append_only/inhibition/00_outputs/<agent>/`**,
  resolved by `shared/outputs.py`, which versions every write and resolves
  every read to the newest. Nothing derived belongs in the repo.
* **Permissions are governed by an Isilon ACL the client cannot see.** The
  POSIX mode `ls` shows is an approximation; `/data` is NFSv3, so
  `nfs4_getfacl` returns "not supported" and `chmod` may not change what is
  actually enforced. If a directory is unreadable despite a mode that plainly
  permits it, that is the ACL and it needs a storage admin.
* **Two receptors are live at once.** `config/receptor.yaml` and
  `noncovalent_dock_run.py` default to 6VAJ; the benchmark and reference-screen
  paths guard hard for 3IKD_ian. Both are populated and plausible, so nothing
  errors — you simply get whichever your entry point chose. See §1 and
  [`retrospective_2.2.0.md`](retrospective_2.2.0.md) §3.3.
* **A guard only protects the path that calls it.** `elevation_report.py` knows
  the MD system renumbers from 1 (Cys113 is residue 63, `PIN1_OFFSET = 50`) and
  refuses to mislabel a structure. New styling code in `md_movie.py` bypassed it
  and drew a **glutamate** labelled as the target cysteine. Route renderers
  through the guard rather than trusting each new path.
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
