# Decisions

!!! info "Generated page"
    Rendered at build time from the repo's source of truth. Edit the underlying file, not this page.

37 record(s). Each answers *why* a choice was made; the [runbooks](../runbooks/index.md) answer *how* to make that kind of choice again, and a run's `manifest.json` records *what* it actually consumed.

Records marked `origin: adversary` are decisions the adversarial review forced — the audit trail that review changed the design.

## Shared substrate

### D0001 — Use PDB 6VAJ as the single shared receptor

**:material-check: accepted** · `origin: spec` · 2026-07-27

All four approaches must dock into the identical prepared receptor or their scores are not comparable. A covalently-bound reference ligand also defines the box for free and hands over the exact attachment atom.

**Decision.** 6VAJ is the shared receptor for the whole choreography. Verified first-hand rather than on the spec's assertion (this closed adversary finding M6). It is hash-pinned; changing it mid-choreography is forbidden.

**Consequences.** Every cross-approach comparison depends on this file. Any future retarget swaps this one entry and re-runs receptor_prep, but all previously computed comparisons are invalidated and must be recomputed.

??? note "Evidence"
    - LINK record: SG CYS A 113 - C10 QT7 A 201 at 1.78 A (a real covalent bond length)
    - resolution 1.42 A
    - TITLE: CRYSTAL STRUCTURE ANALYSIS OF HUMAN PIN1
    - QT7 ligand present, 16 atoms
    - sha256 820fd5969131bef8... pinned in sources.lock.json

Affects: `config/receptor.yaml`, `shared/receptor_prep.py`, `config/sources.yaml`

Runbook: `docs/runbooks/receptor_selection.md`

---

### D0002 — Emit two docking boxes rather than one

**:material-check: accepted** · `origin: adversary` · 2026-07-27

Adversary finding M5: a box drawn around a covalent ligand is centred on the warhead sub-pocket. That is right for the covalent approaches, which attack that atom, and wrong for the non-covalent ones, which would be biased toward a sub-pocket they have no reason to prefer.

**Decision.** receptor_prep.py emits box.json (20 A, covalent, T_3/T_4) and box_expanded.json (26 A, full PPIase pocket, T_1/T_2). Each records which approaches use it.

**Consequences.** Non-covalent and covalent dock scores are computed in different volumes and are therefore NOT directly comparable even before the tool difference is considered. This reinforces the no-authoritative-cross-approach-join design.

??? note "Evidence"
    - QT7 is covalent at Cys113, so a tight box is centred on the warhead sub-pocket
    - covalent box 20 A (t3,t4); expanded box 26 A (t1,t2)
    - Cys113 SG is 4.26 A from the box centre, inside both

Affects: `config/receptor.yaml`, `shared/receptor_prep.py`

Runbook: `docs/runbooks/receptor_selection.md`

---

### D0003 — Retain and report unknown heteroatoms rather than stripping by default

**:material-check: accepted** · `origin: implementation` · 2026-07-27

Stripping every HETATM would silently remove structural cofactors and metals; keeping everything that is not water would silently leave cryoprotectants occupying pocket volume. Both failures are invisible - docking succeeds either way.

**Decision.** receptor_prep.py strips only an explicit list of solvent/buffer/cryoprotectant codes. Anything unrecognized is KEPT and counted in prep_log.json as other_het_atoms_retained. A non-zero count is a signal to go look.

**Consequences.** Preparing a new target requires a heteroatom inventory pass rather than a blind run; the runbook makes that a step. The mechanism proved itself immediately by catching PG4.

??? note "Evidence"
    - 6VAJ heteroatoms: 132 HOH, 31 PG4, 16 QT7, 10 SO4
    - PG4 (cryoprotectant) was retained by the first run and reported, then added to the strip list
    - PG4 nearest atom was 22.65 A from the box centre, so no result was affected
    - final prep: 1215 protein atoms kept, 0 unrecognized retained

Affects: `shared/receptor_prep.py`

Runbook: `docs/runbooks/receptor_selection.md`

---

### D0007 — Keep pins in a lockfile, never write back into hand-authored config

**:material-check: accepted** · `origin: implementation` · 2026-07-27

Acquisition needs to record observed hashes so they can be enforced on later runs. The obvious implementation - write them back into the source config - destroys anything a human wrote there.

**Decision.** config/sources.yaml is hand-authored and never written by code. Observed hashes and resolved commits go to config/sources.lock.json. Config is written by humans; pins are written by code; the two never share a file.

**Consequences.** Standard dependency-lockfile pattern. Adding a source means editing YAML; pinning happens automatically on first acquisition and is enforced thereafter.

??? note "Evidence"
    - the first sources.py round-tripped sources.yaml through yaml.safe_dump
    - this erased all 22 explanatory comments in the file
    - the file had not been committed, so git could not restore it
    - pin enforcement verified by tampering: mismatch raises SourceError

Affects: `config/sources.yaml`, `config/sources.lock.json`, `shared/sources.py`

---

### D0009 — Protonate with reduce, not obabel, and assert post-conditions

**:material-check: accepted** · `origin: implementation` · 2026-07-27

The first `receptor_prep.py` protonated with `obabel -p 7.4`. That is not a
receptor-preparation tool: on 6VAJ it renumbered residues from 1, **renamed**
them, invented four extra chain IDs, and silently dropped 28 of 150 residues —
including **Cys113**, the catalytic residue that T_3 and T_4 both target.

The output still looked like a protein. It would have docked without error and
produced plausible, meaningless scores. It was hash-pinned into a manifest and
described in published documentation before anyone noticed.

It was caught only incidentally: the DiffSBDD smoke test tried to look up
pocket residue 59 and got a `KeyError`.

The deeper failure was in verification, not tool choice. The module checked that
Cys113 existed in the **input** and never re-checked the **output**, and the
atom counts it logged were taken before protonation — so the log looked healthy.

**Decision.** Protonate with `reduce -BUILD` (AmberTools; Word et al. 1999), which preserves
the heavy-atom record. Convert to PDBQT with `obabel -xr` and **no** `-p` flag,
which is safe.

Add post-conditions checked on the **outputs**, not the inputs: residue count
preserved, chain set unchanged, no residue renamed, and the catalytic residue
present with the right name in both the PDB and the PDBQT. Failure raises and
refuses to produce a receptor.

**Consequences.** M1 had to be redone; the earlier prepared receptor and its recorded hashes were
invalid. `reduce` lives in the `amber_md` env rather than `cheminf`, so the
binary is searched for rather than assumed on PATH.

The general lesson is recorded in the receptor runbook: **verify the artifact
you produced, not the one you consumed.** Any preparation step that can silently
drop part of a structure needs a post-condition, not a log line.

??? note "Evidence"
    - obabel -p 7.4 on 6VAJ dropped 28 of 150 residues, including Cys113
    - it also renamed residues (LYS A 6 -> TRP) and invented chains A I J M N
    - reduce -BUILD preserves all 150 residues, chain A only, all names identical
    - obabel -xr without -p is fine for pdbqt: Cys113 SG at exact raw coordinates
    - guard verified against the corrupt output: catches all four failure modes

Affects: `shared/receptor_prep.py`, `config/receptor.yaml`

Runbook: `docs/runbooks/receptor_selection.md`

---

### D0010 — Isolate pip installs behind the target env's bin

**:material-check: accepted** · `origin: implementation` · 2026-07-27

The five-env design exists so three incompatible torch builds cannot collide.
That protects against *conda* cross-contamination, but not against a third-party
installer shelling out to a bare `pip`.

REINVENT's `install.py` builds `cmd = ["pip", "install", ...]`. Invoked as
`$ENV/bin/python install.py`, the Python is the env's but the `pip` is whatever
PATH resolves first — here base conda's. REINVENT and roughly thirty chemistry
packages went into the base environment, and pandas and pydantic were upgraded
out from under base tooling that pins them.

Nothing failed loudly. The install reported success; the env stayed empty; the
damage was in a different environment entirely.

**Decision.** Any third-party installer that may shell out to `pip` is run with the target
env's `bin` FIRST on PATH, and the result is verified by importing from that
env rather than trusting the installer's exit code.

Base conda is treated as off-limits. It belongs to the user's other projects and
is not part of this choreography's environment contract.

**Consequences.** `setup_envs.sh` now sets PATH explicitly for the REINVENT build and asserts the
CLI runs afterward. Base was repaired to pandas 2.3.3 / pydantic 1.10.26; the
remaining `pip check` failures there (requests, urllib3, tqdm, numpy) predate
this and are not ours.

The wider lesson matches D0009: verify the artifact you produced. An installer
exit code of 0 says the installer ran, not that the package landed where you
intended.

??? note "Evidence"
    - REINVENT's install.py shells out to a bare `pip`, not sys.executable -m pip
    - it resolved to base conda pip and installed reinvent + ~30 packages into base
    - pandas was upgraded 2.x -> 3.0.5, breaking streamlit 1.55.0 and anndata 0.12.8
    - pydantic was upgraded <2 -> 2.13.4, breaking anaconda-cloud-auth 0.1.3
    - repaired by pinning pandas<3 and pydantic<2 back into base
    - plan also stale: REINVENT needs python >=3.11 (built 3.10), and install.py first positional is processor_type not the dependency set

Affects: `scripts/setup_envs.sh`

---

### D0011 — gnina CNN scoring is uncalibrated for covalent docking

**:material-check: accepted** · `origin: implementation` · 2026-07-27

M2 pinned the shared gnina covalent protocol and ran it end-to-end on Sulfopin,
the anchor whose true covalent pose 6VAJ resolves. It works. But gnina prints,
on every covalent run:

> CNN scoring not yet calibrated for covalent docking. Recommend running with
> `--cnn_scoring none`

The Rev 3 spec makes `CNNaffinity` **T_3's rank metric** and T_4's secondary
metric, and §7 offers an optional within-covalent-stratum re-score built on it.
The tool's own authors are saying that number is not calibrated for the mode we
are using it in.

This is the same shape as the LUMO finding in D0005: a metric the spec treats as
a ranking signal, which the evidence says is not one.

**Decision.** **ACCEPTED 2026-07-27: option 2 below — let M3 decide empirically.**

The question is settled by measurement rather than judgement: the enrichment
gate runs the covalent stratum with BOTH candidate metrics and the one that
enriches known actives is the one that ranks. If neither enriches, docking is
demoted to a displayed label in T_3 and T_4, which is the gate's existing
FAIL branch.

Options as weighed:

1. **Rank covalent candidates by gnina's Vina-style `affinity` (kcal/mol,
   lower better) and carry `CNNaffinity` as an advisory annotation.** Keeps the
   shared-protocol parity (S3) intact, uses a metric that is at least
   calibrated for what it measures, and costs nothing to adopt.
2. **Keep CNNaffinity but run the enrichment gate (M3) on the covalent stratum
   specifically**, and let the measured result decide. This is what the gate
   exists for, and it turns the question into an empirical one.
3. Keep CNNaffinity as the rank metric and note the caveat. Weakest option: it
   ranks on a number the tool says is uncalibrated.

Option 2 subsumes option 1 and is the honest route, at the cost of needing M3
before T_3 can rank anything. **This is the accepted path.** Consequence: the
protocol's `cnn_scoring` pin — and therefore the protocol fingerprint — must be
settled by M3 BEFORE T_3 or T_4 produce any results, not after.

Interim, already implemented: `dock()` detects the warning, logs it, and
returns `cnn_uncalibrated_for_covalent` so it reaches the manifest and the GUI
rather than scrolling past in a log.

**Consequences.** If CNNaffinity is demoted, T_3's output contract changes (rank metric, and its
stated direction flips from higher-better to lower-better), and §7's
within-covalent re-score changes with it. The docs state the direction
explicitly in several places and would need updating together.

Nothing is blocked meanwhile: the protocol is pinned, parity holds, and both
numbers are recorded per dock regardless of which one ends up ranking.

??? note "Evidence"
    - gnina v1.3.3 prints: "CNN scoring not yet calibrated for covalent docking. Recommend running with --cnn_scoring none"
    - the Rev 3 spec makes gnina CNNaffinity T_3 rank metric and T_4 secondary metric
    - Sulfopin covalent dock: CNNaffinity 4.64, CNNscore 0.6286, Vina-style affinity -2.35 kcal/mol
    - the warning is emitted on every covalent run, not just edge cases

Affects: `shared/covalent_protocol.py`, `config/choreography.yaml`, `docs/approaches/t3.md`, `docs/approaches/t4.md`

---

### D0012 — Gates report evidence strength; they do not adjudicate

**:material-check: accepted** · `origin: user` · 2026-07-27

The enrichment gate was specified as a binary PASS/FAIL on ROC-AUC, EF1% and
BEDROC thresholds. Assessing the available actives showed those thresholds
cannot carry that weight: six actives, three or four independent chemotypes, and
an EF1% that is quantised to a handful of values.

The tempting responses were both wrong. Emitting a confident PASS from six
actives manufactures precision that is not there. Emitting FAIL, or refusing to
run until the statistics are strong, discards real signal — and the statistics
will *never* be strong here, because validated Pin1 chemistry is genuinely
scarce. That scarcity is a finding about the target, not a defect to engineer
around.

The PI's framing settles it: **this is not purely a statistical exercise.** The
choreography exists so PIs and researchers can bring their own priors and
expertise to bear on the ranking. False discovery is expected and accepted; the
job is to compute where the most promising leads are, not to prove significance.

**Decision.** Gates **report evidence strength; they do not adjudicate.** Concretely:

- The enrichment gate emits a graded verdict — `STRONG`, `WEAK`,
  `UNDERPOWERED`, `FAIL` — not a binary, and always alongside confidence
  intervals, the actives count, and the **independent chemotype count**.
- Evaluation is **per chemotype**, leave-one-chemotype-out, so analog bias
  cannot inflate a result. Six actives that are three chemotypes get reported as
  three.
- `UNDERPOWERED` does **not** mean "discard". It means the ranking is carried
  forward with its uncertainty displayed, for a human to weigh. It is a label on
  the evidence, not a veto.
- Only `FAIL` — docking demonstrably anti-correlated with known actives —
  demotes `dock_score` to a displayed label.

This makes the gate consistent with the choreography's existing stance in Rev 3
section 7: present the evidence and its limits, let the human adjudicate. A gate
that silently vetoed on thin statistics would be the one component doing the
opposite of what everything else does.

**Consequences.** The GUI must display gate verdicts with their power characteristics, not just a
green tick — a `WEAK` verdict shown as a pass would be worse than no gate. The
Open Questions panel already surfaces this class of limitation.

D0011 still stands: M3 decides the covalent rank metric empirically. But if the
comparison between CNNaffinity and Vina-style affinity comes back
`UNDERPOWERED`, that is a real answer — it means the choice should be made on
mechanistic grounds and the tool's own calibration warning, not on six data
points.

??? note "Evidence"
    - covalent stratum has 6 verified actives spanning ~3-4 independent chemotypes
    - ROC-AUC standard error at n=6 is roughly +/-0.2 — a 0.70 PASS sits within noise of FAIL
    - EF1% over 306 molecules resolves to the top 3 compounds, so it is quantised, not continuous
    - the choreography already refuses an authoritative cross-approach numeric join (Rev 3 section 7)

Affects: `config/gates.yaml`, `shared/enrichment_gate.py`, `integration/app/DECISIONS_TAB_SPEC.md`

---

### D0014 — Covalent decoys must carry a warhead

**:material-check: accepted** · `origin: implementation` · 2026-07-27

The first covalent decoy set was property-matched on size, greasiness,
polarity and charge — the standard DUD-E criteria — but never checked whether a
decoy could *react*. Covalent docking requires a reactive atom to bond to;
gnina matches a SMARTS pattern against the ligand and a molecule with no
electrophile cannot be scored at all.

Only 10.6% of that set carried an electrophile. Nine decoys in ten were
unrunnable, and the comparison that survived would have been "electrophiles
versus inert molecules" — which docking wins for reasons that have nothing to do
with Pin1.

**Decision.** Covalent decoys must carry a warhead motif, enforced at generation. Where the
pool allows, they are drawn from the SAME warhead class as the active they are
matched to.

That second part is what makes the control interesting rather than trivial. The
question becomes *does docking prefer our chloroacetamide over other
chloroacetamides* — discrimination within a chemistry — instead of *does docking
prefer electrophiles over inert molecules*, which is not in doubt and not
informative.

**Consequences.** 294 decoys, none of them unrunnable. But class-matching is only partial: 119 of
294. The ChEMBL pool holds 13 chloroacetamides and no sulfamate acetamides in
the relevant property range, so most decoys are Michael acceptors topped up from
outside the active's class. That fallback is recorded per decoy
(`class_matched`) rather than hidden, and it caps how strong a covalent verdict
can honestly be.

A targeted ChEMBL substructure search per warhead class would deepen the pool
and raise the class-matched fraction. Worth doing if the covalent gate turns out
to be the binding constraint on T_3/T_4.

The peptidomimetic BJP-06-005-3 still matches nothing and remains excluded by
the minimum-decoy filter, as before.

??? note "Evidence"
    - decoys_covalent_1: only 32 of 302 (10.6%) carried any electrophile
    - gnina covalent docking requires --covalent_lig_atom_pattern to MATCH the ligand
    - decoys_covalent_2: 294 decoys, 0 without an electrophile
    - 119 of 294 share their active warhead class; the ChEMBL pool holds 13 chloroacetamides and 0 sulfamate acetamides
    - peptidomimetic BJP-06-005-3 still matches 0 warhead-bearing decoys and is excluded by the >=10 filter

Affects: `shared/decoys.py`, `data/reference/warhead_classes_3.csv`

---

### D0015 — Covalent ranking uses gnina affinity, not CNNaffinity

**partially_withdrawn** · `origin: implementation` · 2026-07-27

D0011 recorded that the Rev 3 spec makes gnina `CNNaffinity` T_3's rank metric
and T_4's secondary, while gnina itself prints a warning on every covalent run
that CNN scoring is not calibrated for covalent docking. Rather than decide on
the warning alone, D0011 accepted the empirical route: run both metrics through
the enrichment gate and let the measurement decide.

M3 ran 6 actives against 294 warhead-bearing decoys.

**Decision.** **T_3 and T_4 rank on gnina's Vina-style `affinity` (kcal/mol, lower better).**
`CNNaffinity` is carried as an advisory annotation, never as a rank metric.

The two metrics separate cleanly on the evidence that matters. `affinity_kcal`
enriches: its confidence interval excludes 0.5 and it puts actives in the top
1% (EF1% 16.7). `CNNaffinity` does neither — its interval includes 0.5, which is
consistent with no enrichment at all, and an EF1% of 0.0 means **not one known
active reached the top 1%** of its ranking.

That is independent confirmation of the tool's own warning, arrived at without
relying on it.

**Consequences.** T_3's output contract changes: the rank metric is now kcal/mol and
**lower-is-better**, where the spec had it dimensionless and higher-is-better.
The direction is stated explicitly in several places and all of them move
together. Section 7's within-covalent re-score follows the same metric.

The protocol keeps `cnn_scoring: rescore` so the advisory number is still
produced; only its ROLE changes. The fingerprint is therefore unaffected by this
decision.

Both verdicts remain UNDERPOWERED at 4 independent chemotypes, so this is a
comparison between two metrics on the same data rather than a claim that
covalent docking is validated. Under D0012 the ranking carries forward with its
uncertainty displayed.

??? note "Evidence"
    - covalent affinity_kcal: ROC-AUC 0.815, CI [0.667, 0.931] EXCLUDES 0.5, EF1% 16.7
    - covalent CNNaffinity: ROC-AUC 0.707, CI [0.408, 0.921] INCLUDES 0.5, EF1% 0.0
    - EF1% of 0.0 means not one known active landed in the top 1% of the CNNaffinity ranking
    - gnina itself warns CNN scoring is not calibrated for covalent docking
    - both verdicts capped at UNDERPOWERED by the 4-chemotype floor

Affects: `config/choreography.yaml`, `shared/covalent_protocol.py`, `docs/approaches/t3.md`, `docs/approaches/t4.md`

---

### D0016 — Non-covalent docking barely enriches on Pin1

**:material-check: accepted** · `origin: implementation` · 2026-07-27

The enrichment gate ran the non-covalent stratum through the same Vina protocol
T_1 and T_2 will use, against property-matched decoys on the prepared 6VAJ
receptor.

**Decision.** Record that **non-covalent docking shows essentially no enrichment on this
target**, and treat T_1 and T_2 dock-based rankings as weakly supported until
more actives are available.

ROC-AUC 0.535 is a coin flip. The interval [0.215, 0.855] is so wide it cannot
distinguish "docking works" from "docking is actively misleading". EF1% 0.0 and
BEDROC 0.083 agree: no known binder reached the top of the ranking.

This does NOT trip FAIL, which requires the interval's upper bound to fall below
0.5, so under D0012 the ranking carries forward with its uncertainty displayed
rather than being vetoed.

**Consequences.** **T_1 and T_2 both rank on Vina.** On this evidence their dock-based shortlists
should not be presented as evidence-backed, and the GUI must show the gate
verdict beside them rather than a bare score. The covalent stratum is in
markedly better shape (AUC 0.815, CI excluding 0.5), so the two families are not
equally supported and should not be displayed as if they were.

This is a plausible property of the target rather than a defect in the setup.
Rev 3 section 3 already notes Pin1's PPIase pocket is shallow and
solvent-exposed — the regime where structure-based methods are weakest, and the
reason T_1's sanitise/filter stages were called load-bearing.

Before concluding docking cannot rank here, the honest next step is more
actives: 5 non-covalent actives over 5 chemotypes cannot separate a real null
from a small effect. If the result survives an expanded set, T_1 and T_2 should
lean on their other evidence and the human, exactly as the FAIL branch intends.

??? note "Evidence"
    - non-covalent Vina: ROC-AUC 0.535, CI [0.215, 0.855], EF1% 0.0, BEDROC 0.083
    - EF1% 0.0 — no known Pin1 binder reached the top 1% of the ranking
    - 5 actives, 243 decoys, 5 independent chemotypes
    - the CI spans from clearly-worse-than-chance to strongly-enriching
    - Rev 3 section 3 already flags the PPIase pocket as shallow and solvent-exposed

Affects: `docs/approaches/t1.md`, `docs/approaches/t2.md`, `integration/app/DECISIONS_TAB_SPEC.md`

---

### D0017 — Vina-GPU adopted at search_depth >= 20

**:material-check: accepted** · `origin: implementation` · 2026-07-27

AutoDock Vina is CPU-only, making it the throughput wall for T_1 (5-10k
candidates) and T_2. Vina-GPU 2.1 was built and validated against the 248
ligands M3 had already scored with CPU Vina, on the identical receptor and box.

Adoption required four checks: rank agreement, score agreement, that the
enrichment ROC-AUC reproduces within 0.10, and that the graded verdict matches.

**Decision.** **RESOLVED 2026-07-27 by sweeping `search_depth`: ADOPT at >= 20.**

The original run used `search_depth 10` and failed the AUC check by 0.001.
Sweeping the parameter settled it — the discrepancy was search **convergence**,
not implementation:

| search_depth | Spearman | mean abs diff | AUC drift | verdict |
|---|---|---|---|---|
| 10 | 0.869 | 0.229 | 0.101 | DO_NOT_ADOPT |
| 20 | 0.909 | 0.158 | **0.005** | **ADOPT** |
| 40 | 0.946 | 0.105 | 0.006 | **ADOPT** |

Every metric improves monotonically with depth, which is what a convergence
explanation predicts and an implementation difference would not. At depth 20 the
AUC drift is 0.005 — twentyfold inside the threshold — and all four checks pass.

The literature agrees independently: Tang et al. (Molecules 2022,
10.3390/molecules27093041) report CPU -8.9 vs GPU -8.7 kcal/mol with Pearson
0.965 and set 0.5 kcal/mol as their own agreement tolerance. Our +0.188 bias at
depth 10 reproduces their +0.2 almost exactly, and every mean absolute
difference here sits inside their tolerance.

The original reasoning below stands as the record of why the depth-10 failure
was not itself informative.

---

**(superseded reasoning, retained)** Not adopted on the depth-10 run — but the
failing check was not informative, and the score agreement was good.

Three of four checks pass. Rank correlation is 0.869 and scores differ by 0.229
kcal/mol on average with a small positive bias; both engines return the same
graded verdict (UNDERPOWERED) and the same EF1% of 0.0.

The failure is ROC-AUC drift: 0.101 against a 0.10 threshold, missing by 0.001.
That check is close to meaningless on this data. **The baseline AUC is 0.535 —
chance.** When there is no signal to preserve, ranking is decided by noise, so
small score changes move the AUC freely. The two confidence intervals overlap
across 0.238-0.624 and *both contain 0.5*: the engines are not disagreeing about
which molecules bind, they are agreeing that neither can tell.

Adopting on a knife-edge threshold would be as wrong as rejecting on one. The
honest position is that **score agreement is demonstrated and enrichment
agreement is untested**, because this data set cannot test it.

**Consequences.** Vina-GPU stays built, wrapped and documented but out of the pipeline. CPU Vina
remains the T_1/T_2 engine, so nothing downstream changes and D0016's baseline
stands.

To settle it properly, one of:

- **Re-validate on a set where docking actually enriches.** An AUC-reproduction
  test is only meaningful against a non-null baseline. If the expanded actives
  set lifts non-covalent enrichment above chance, re-run this script.
- **Raise `--search_depth`.** DONE — this is what resolved it. Adopt at >= 20.
- **Adopt on score agreement alone**, accepting that enrichment equivalence is
  unproven, if T_1's throughput becomes the binding constraint. A 3.5x speedup
  on 248 ligands should widen on 10,000, where kernel setup amortises.

The speedup is real but was measured on a small batch; it is not yet a
projection for T_1 scale.

??? note "Evidence"
    - Spearman rho 0.869 and Pearson 0.890 against CPU Vina over 248 ligands
    - mean absolute difference 0.229 kcal/mol, mean bias +0.188
    - ROC-AUC 0.535 (CPU) -> 0.433 (GPU); drift 0.101 against a 0.10 threshold
    - both CIs contain 0.5 and overlap over 0.238-0.624; EF1% is 0.0 for both
    - 346.9 s on GPU vs ~20 min on 20 CPU cores — roughly 3.5x on 248 ligands
    - SWEEP: search_depth 20 -> Spearman 0.909, mean |diff| 0.158, AUC drift 0.005, ADOPT
    - SWEEP: search_depth 40 -> Spearman 0.946, mean |diff| 0.105, AUC drift 0.006, ADOPT
    - so the discrepancy was search CONVERGENCE, not implementation

Affects: `scripts/build_vina_gpu.sh`, `scripts/validate_vina_gpu.py`

---

### D0025 — isolate_rgroup ignores its cap argument, producing false aldehyde alerts

**:material-check: accepted** · `origin: implementation` · 2026-07-27

T_3 decorates a fixed scaffold at the sulfopin nitrogen, and LibInvent
frequently attaches its decoration through a carbonyl — ureas, carbamates,
amides. That is ordinary medicinal chemistry.

`shared.alerts.isolate_rgroup` cuts the R-group away from the core so alerts are
scored on the decoration rather than on the warhead every candidate is required
to carry. It takes a `cap` argument, documented as the group used to satisfy the
open valence left by the cut, defaulting to `"C"`.

**The argument is ignored.** `cap="C"`, `cap="[H]"` and `cap="CC"` all return
byte-identical output. The cut valence is filled with hydrogen regardless.

For a decoration attached through carbon this is harmless. For one attached
through a carbonyl it is not: an amide `>N–C(=O)–R` becomes, once severed from
the nitrogen and H-capped, a formamide `H–C(=O)–R`. BRENK then flags an aldehyde
that exists only in the fragment, never in the molecule.

**Consequences.** **T_3's rejection rate barely moved — 97.7% to 77.2% — and that is the
interesting part.** The aldehyde and thiol artefacts are gone entirely, but the
dominant alert is now `acyclic_imide` (3,787 candidates), and it is REAL. The
T_3 scaffold nitrogen already carries the acrylamide carbonyl; LibInvent
frequently decorates it with a second acyl group, and an N with two acyl groups
is an acyclic imide — genuinely more electrophilic and hydrolytically labile
than the acrylamide alone.

So the finding is not "the gate was broken", it is **"LibInvent's preferred
decoration chemistry for this scaffold creates imides"**. That is a fact about
T_3 worth carrying to the panel, and it would have stayed hidden underneath the
aldehyde artefact. 1,233 of 5,396 now pass.

`max_rgroup_alerts = 0` remains the default and remains worth revisiting
separately: rejecting on a single alert is strict for an approach whose value is
proposing chemistry a person would not have picked.

T_4 is affected in principle but far less in practice: its R-groups come from a
frequency-derived ChEMBL library and attach through carbon, not carbonyl. The
count of T_4 candidates rejected specifically for `aldehyde` should be checked
before assuming its gate was unaffected.

The general shape is the same one D0022 had: a control that ran, reported a
plausible number, and was measuring an artefact of how its input was prepared.
Both were found by looking at a specific molecule rather than at a summary
statistic.

??? note "Evidence"
    - isolate_rgroup(cap="C"), cap="[H]" and cap="CC" all return the identical R-group SMILES — the argument has no effect
    - T_3 candidate C=CC(=O)N(C(=O)NCc1ccccc1OC)C1CCS(=O)(=O)C1 has a UREA decoration; its isolated R-group reads COc1ccccc1CNC=O, a formamide
    - BRENK correctly flags H-C(=O)-N as an aldehyde; the intact molecule has no aldehyde
    - 2657 of 5396 T_3 candidates were stamped with "aldehyde, aldehyde", the single largest rejection reason
    - two_tier defaults to max_rgroup_alerts=0, so a single artefactual alert rejects
    - 5270/5396 (97.7%) of T_3 stamped rejected at the alert gate

Affects: `shared/alerts.py`, `approaches/t3_reinvent/02_annotate.py`, `approaches/t4_combinatorial/01_enumerate.py`

---

### D0026 — Excuse named alerts rather than raising the tolerated count

**:material-check: accepted** · `origin: user` · 2026-07-28

The decoration gate rejected on a single alert, failing 77% of T_3. The natural
response is to raise the tolerated count, and the PI asked whether we could.

The distribution says a count is the wrong instrument. Of the 3,054 T_3
candidates carrying exactly one alert, **2,853 carry `acyclic_imide` and nothing
else**. So "tolerate one alert" would have been a decision to accept imides,
taken by accident, while simultaneously admitting every other one-off alert —
including thioesters and phenol esters, which hydrolyse, and thiocarbonyls,
which are reactive. Those are worse liabilities than the thing the change was
actually meant to permit.

**Decision.** **Keep the count at 0. Excuse alerts by NAME, and carry the excusal.**

For T_3, `acyclic_imide` is excused. The scaffold nitrogen already bears the
acrylamide carbonyl, and LibInvent frequently adds a second acyl group, which
makes a genuine imide. This is not a false positive: an N-acyl acrylamide is
more electrophilic and more hydrolytically labile than the acrylamide alone.

It is excused rather than rejected because it is a structural consequence of
*where T_3 decorates*, so gating on it discards 53% of the approach before
anything has been measured — the filter-before-you-measure trap that D0012 and
D0019 both refuse. The alert travels with the candidate as
`excused_alert_names`, and the GUI must display it: an imide shown without that
caveat implies a cleanliness it did not earn.

**Consequences.** T_3 goes from 1,233 to 4,086 passing. Every one of the 2,853 imide-bearing
candidates carries a visible flag rather than a silent pass.

This changes nothing for T_1 and T_2, and the record should be explicit because
it was briefly misremembered: the gate lives on the two-tier path, which runs
only for approaches that pass a `core_smarts`. T_1 and T_2 are non-covalent,
have no fixed core, and take the whole-molecule path where nothing is
disqualifying by default. **Neither has ever had an alert-gate rejection.** T_1's
own lever is `--alert-limit`, still unset, so its alerts are annotated only.

T_4 is barely sensitive to the count (1,683 of 1,782 at 0, 1,728 at 1) and gets
no excusals; its decoration alerts are genuine and few.

If docking and MM-GBSA end up favouring the imides, the honest resolution is a
hydrolytic-stability measurement, not a recomputed alert.

??? note "Evidence"
    - T_3: 3054 candidates carry exactly one attributable alert, and 2853 of those (93%) are acyclic_imide alone
    - raising the count to 1 would also admit thioester, phenol_ester, Thiocarbonyl_group and isolated_alkene one-offs
    - T_3 passing: 1233 (count 0) -> 4086 (count 0, imide excused); rejected 4163 -> 1310
    - T_4 is barely sensitive: 1683 pass at count 0, 1728 at count 1 of 1782
    - the gate only applies to approaches passing core_smarts — T_1 and T_2 have zero alert-gate rejections

Affects: `shared/alerts.py`, `shared/annotate.py`, `config/approaches/t3_reinvent.yaml`, `approaches/t3_reinvent/02_annotate.py`

---

### D0028 — Enrichment gate re-measured on adduct forms — D0015's decisive evidence does not survive

**:material-check: accepted** · `origin: implementation` · 2026-07-28

D0022 changed what gets docked: the adduct form rather than the pre-reaction
ligand. That record noted the enrichment gate had been measured through the old
protocol and said it "should be re-measured rather than assumed", since actives
and decoys were treated identically and ROC-AUC was unlikely to move much.

That prediction was wrong in the way that matters.

**Decision.** **D0015's conclusion is retained; its stated justification is not.**
`affinity_kcal` remains T_3 and T_4's rank metric, because it beats the
alternative on every statistic and because CNNaffinity is uncalibrated for
covalent docking by the tool's own admission. But the record must not keep
citing an interval that excludes 0.5, because it no longer does.

**The covalent gate's verdict stays UNDERPOWERED**, which is what it has always
said. The graded-verdict floor (6 independent chemotypes, gates.yaml) refused to
claim more than UNDERPOWERED even when the point estimates looked good, and that
refusal is now vindicated: the flattering point estimate did not survive a
change in the ligand form.

`decoys_covalent_3.csv` records the verified class assignment for future runs.

**Consequences.** **T_4's ranking is weakly supported and should be presented that way.** It was
already displayed with its gate verdict; the verdict has not changed, only the
confidence one should place in the numbers behind it.

**The decoy set needs regenerating before the gate can say anything stronger.**
Property matching alone is not enough — decoys must be matched on warhead class
too, or the gate cannot separate binding discrimination from chemotype
discrimination. That is D0014's territory and is not attempted here.

**Assigning a chemotype with a reactive-atom SMARTS is a category error** and it
appeared twice in this project: here, and in the alert attribution, where
excusing only the reactive atoms let `alpha_halo_carbonyl` straddle the boundary
(D0025). The reactive-atom pattern says where a bond forms. The whole-warhead
fragment says what the chemistry is. They are not interchangeable and the
library carries both for that reason.

??? note "Evidence"
    - affinity_kcal on adduct forms: ROC-AUC 0.718, CI [0.483, 0.944] — the CI now INCLUDES 0.5
    - D0015 measured the same metric at ROC-AUC 0.815, CI [0.667, 0.931], which EXCLUDED 0.5
    - cnn_affinity on adduct forms: ROC-AUC 0.392, CI [0.181, 0.645], EF1% 0.0 (was 0.707)
    - affinity EF1% rose 16.7 -> 19.0; BEDROC 0.333 vs cnn 0.146
    - both verdicts remain UNDERPOWERED: 4 independent chemotypes against a floor of 6
    - decoy warhead classes were assigned by the NARROW reactive-atom SMARTS; [CH2][Cl] matched nitrogen mustards (cyclophosphamide) and nitrosoureas (lomustine) as chloroacetamides
    - only 112 of 294 decoys carry a whole warhead group AND survive the adduct transform
    - verified decoys are 104 acrylamide, 5 naphthoquinone_c2, 3 chloroacetamide
    - two of six actives (sulfamate_acetamide) and one (snar_chloroazine) have NO same-class decoy

Affects: `scripts/run_enrichment_gate.py`, `decisions/D0015-covalent-ranking-uses-affinity-not-cnnaffinity.md`, `data/reference/decoys_covalent_3.csv`, `config/gates.yaml`

---

### D0030 — Acrylamide's adduct is saturated; the quinones' is not — one mechanism, two chemistries

**:material-check: accepted** · `origin: implementation` · 2026-07-28

**Decision.** 1. Acrylamide's adduct form saturates the acceptor C=C. Its attachment
   SMARTS becomes `[CH3][CH2][CX3](=O)[NX3]` — the terminal carbon of a
   propanamide, exactly parallel to how the acetamide classes bond the
   CH3 of `[*]C(=O)C`.
2. The quinones are unchanged; their note now states the
   re-aromatization assumption instead of claiming a missing hydrogen.
3. Whether a Michael acceptor saturates is declared in the library
   (`adduct_saturates_alkene`, `warhead_classes_5.csv`), not inferred.
   Both chemistries carry `mechanism: michael_addition` and need
   opposite treatment, so the mechanism label cannot decide it and
   neither can a SMARTS. It is chemistry, and it is written down.
4. The docking run itself moves to `shared/covalent_dock_run.py`. T_3
   and T_4 now execute one function rather than two scripts that began
   identical — the drift risk T_3's own config warns about, applied to
   the loop instead of just the transform.

**Consequences.** The protocol fingerprint changes, so **T_4's existing docks are no
longer at parity with T_3's** and the integration GUI will correctly
refuse the within-covalent comparison until T_4 re-docks. T_4's
non-acrylamide classes will reproduce exactly — the docking is seeded
and deterministic and their SMARTS are untouched — so the re-run's real
content is acrylamide's 187 rows.

`05_regiochemistry_comparison`'s naphthoquinone result is **not**
affected. Both of its arms are quinones, neither is transformed, and its
STRONG verdict (benzo over c2 on pose success, p = 1.2e-14) stands. The
likely cause of c2's 96.3% pose failure is that c2 places the R-group
and the Cys113 sulfur on adjacent carbons of one rigid ring while benzo
places the R-group on the distal ring — a real steric difference, not a
modelling artefact.

??? note "Evidence"
    - D0022 transformed only classes with a leaving group; all three michael_addition classes were passed through untouched
    - 561 of 1,782 T_4 rows carried the resulting approximation note (acrylamide, naphthoquinone_c2, naphthoquinone_benzo)
    - acrylamide is T_3's ONLY warhead, so 100% of T_3 was affected
    - docking the alkene gives Cys-S-CH=CH-C(=O)NR2, a planar vinyl thioether; the true adduct is Cys-S-CH2-CH2-C(=O)NR2 with two rotatable bonds
    - thiol addition to 1,4-naphthoquinone gives a hydroquinone that re-oxidizes to the 2-thio-quinone, so the quinone sulfur genuinely sits on an sp2 carbon
    - protocol fingerprint 67366274f425a371 -> a2854a6e6f7edc43
    - T_3 smoke test: 12/12 docked on the saturated adduct

Affects: `shared/covalent_adduct.py`, `shared/covalent_dock_run.py`, `data/reference/warhead_classes_5.csv`, `approaches/t3_reinvent/03_covalent_dock.py`, `approaches/t4_combinatorial/03_covalent_dock.py`, `decisions/D0022-dock-the-adduct-not-the-pre-reaction-ligand.md`

---

### D0031 — Class-matched decoys remove the apparent covalent enrichment

**:material-check: accepted** · `origin: implementation` · 2026-07-28

**Consequences.** - `decoys_covalent_2` must not be used for a covalent gate again.
- D0015's metric choice stands on mechanism, not measurement.
- The integration GUI must show the covalent ranking with this verdict
  attached, not the D0015 figure.
- Strengthening this gate requires more actives, which means the
  literature, not more compute.

??? note "Evidence"
    - class-matched gate: affinity_kcal ROC-AUC 0.537, CI [0.346, 0.728], EF1% 0.0, BEDROC 0.001
    - class-matched gate: cnn_affinity ROC-AUC 0.552, CI [0.358, 0.741], EF1% 0.0, BEDROC 0.001
    - the same metric measured 0.815 on decoys_covalent_2 (D0015) and 0.718 on adduct forms (D0028)
    - decoys_covalent_2 held 104 acrylamide decoys against ZERO acrylamide actives
    - new set: 90 decoys, every one carrying its chemotype whole reactive group and producing a valid adduct
    - Sulfopin 50/50 same-class decoys; Juglone 31; BJP-06-005-3 8; Tian-6a 0; Reddi-4d 0; Reddi-4g 1
    - ChEMBL holds 4,430 chloroacetamides, 3,963 naphthoquinones, 41 sulfonate acetamides, 6 sulfamate acetamides, and 3 non-Pin1 nitro-chloropyrimidines
    - max ECFP4 Tanimoto of any decoy to any active: 0.344 against a 0.35 cap
    - both verdicts remain UNDERPOWERED: 2 chemotypes against a floor of 6

Affects: `shared/decoys_classmatched.py`, `scripts/build_covalent_decoys.py`, `scripts/run_enrichment_gate.py`, `data/reference/decoy_chemotypes_2.csv`, `decisions/D0014-covalent-decoys-must-carry-a-warhead.md`, `decisions/D0015-covalent-ranking-uses-affinity-not-cnnaffinity.md`, `decisions/D0028-enrichment-gate-remeasured-on-adduct-forms.md`

---

### D0032 — MM-GBSA does not rescue the ranking — and a negative verdict needs as much power as a positive one

**:material-check: accepted** · `origin: implementation` · 2026-07-28

??? note "Evidence"
    - MM-GBSA on the class-matched gate set: ROC-AUC 0.140, CI [0.060, 0.240], EF1% 0.0, BEDROC 0.000
    - docking on the SAME 51 ligands: ROC-AUC 0.440, EF1% 0.0 — MM-GBSA is 0.300 WORSE
    - Sulfopin ranks 44 of 51 by MM-GBSA dG and 29 of 51 by docking affinity
    - Sulfopin dG -15.76; decoy dG median -21.90, best -46.45
    - dG vs heavy-atom count Spearman -0.291; affinity vs heavy-atom count -0.409, so the size artefact does NOT explain it
    - coverage: 51 of 83 ligands scored — all 32 naphthoquinone_c2 failed on the sp2 junction gap
    - only 1 active survives, so the result is UNDERPOWERED in both directions
    - the gate initially graded this FAIL, because its power floor did not govern the FAIL branch

Affects: `scripts/run_mmgbsa_gate.py`, `shared/enrichment_gate.py`, `config/gates.yaml`, `decisions/D0031-class-matched-decoys-remove-the-apparent-covalent-enrichment.md`

---

### D0033 — Every dG in the project was summed from a partial energy, and a plausible number hid it

**:material-check: accepted** · `origin: implementation` · 2026-07-28

??? note "Evidence"
    - ENERGY_TERMS asked for '1-4VDW'/'1-4EEL'; sander prints '1-4 VDW'/'1-4 EEL' with a space
    - the token regex [A-Z0-9\-]+ stopped at the space and stored 1-4 VDW's value under the key 'VDW', which nothing read
    - 1-4 EEL collided with the already-set 'EEL' key under setdefault and was discarded
    - CMAP was never in ENERGY_TERMS at all
    - net effect: three terms contributed exactly 0.0 to every leg total in the project
    - recomputing with the old logic reproduces every stored dG exactly, confirming this produced all of them
    - gate set: shift +17.00 kcal/mol mean, range +7.47 to +38.00, sd 6.63 -- large and non-constant
    - shift differs BY APPROACH: t3 -8.77 mean vs gate +17.00, a ~26 kcal/mol systematic gap
    - T_4 ranking near-inverts: Spearman(original, corrected) = -0.735, 0 of 5 top candidates retained
    - D0032 re-run: MM-GBSA ROC-AUC 0.140 -> 0.260, Sulfopin 44/51 -> 38/51, still below docking 0.440, verdict still UNDERPOWERED
    - 134 candidates recomputed by parsing alone; no minimisation was rerun

Affects: `shared/mmgbsa.py`, `scripts/recompute_mmgbsa_totals.py`, `tests/test_mmgbsa_energy_terms.py`, `decisions/D0032-mmgbsa-gate-and-the-power-floor-on-negative-verdicts.md`

---

### D0034 — The gate token erased sibling metrics — the same defect as before, one level down

**:material-check: accepted** · `origin: implementation` · 2026-07-28

??? note "Evidence"
    - write_token popped the whole stratum before re-adding: `for s in {r.stratum for r in results}: by_stratum.pop(s, None)`
    - run_mmgbsa_gate.py writes a single mmgbsa_dG result for the covalent stratum
    - so it deleted the covalent affinity_kcal verdict (ROC-AUC 0.537) that D0031 established
    - the live token was found carrying covalent metrics = ['mmgbsa_dG'] only
    - with docking gone, recommended_rank_metric became mmgbsa_dG (0.140) — a metric NO approach ranks on
    - demonstrated directly: old logic leaves [mmgbsa_dG], new logic leaves [affinity_kcal, mmgbsa_dG]
    - token rebuilt: covalent affinity_kcal 0.537 + mmgbsa_dG 0.260, recommended back to affinity_kcal

Affects: `shared/enrichment_gate.py`, `tests/test_gate_token_merge.py`, `decisions/D0031-class-matched-decoys-remove-the-apparent-covalent-enrichment.md`, `decisions/D0032-mmgbsa-gate-and-the-power-floor-on-negative-verdicts.md`

---

### D0035 — The sp2 junction gap that cost 32 ligands and an active was three missing angle lines

**:material-check: accepted** · `origin: implementation` · 2026-07-28

??? note "Evidence"
    - tleap error was always the same single line: 'Could not find angle parameter for atom types: 2C - S - cc'
    - junction v2 declared BONDS S-c2, S-c3, S-ca, S-cc, S-cd but only ANGLES 2C-S-c2 and 2C-S-c3
    - so every ligand attaching through aromatic (ca) or conjugated (cc/cd) carbon failed at the last build step
    - cost: 32 of 83 gate ligands, 5 of 7 T_4 adduct classes, and Juglone -- the second active
    - fix is 3 lines, taken from gaff2.dat by the analogue rule the file already documented
    - the 2 pre-existing entries match gaff2.dat exactly (c2-ss-c3, c3-ss-c3), confirming the rule before applying it
    - Juglone now builds: tleap Errors = 0, all three prmtops written

Affects: `data/params/cys_gaff2_junction_3.frcmod`, `shared/mmgbsa.py`, `decisions/D0032-mmgbsa-gate-and-the-power-floor-on-negative-verdicts.md`

---

### D0036 — Better sampling does not rescue MM-GBSA — the ensemble is precise and still below chance

**:material-check: accepted** · `origin: implementation` · 2026-07-29

??? note "Evidence"
    - 167 candidates rescored per-frame over 2 ns GB implicit-solvent MD, 0 failures
    - gate set complete: 82 of 83 ligands, 2 actives, 2 chemotypes
    - docking ROC-AUC 0.537; MM-GBSA single-structure 0.350; MM-GBSA ensemble 0.394
    - propagating each candidate MEASURED SEM into the metric: AUC 95% [0.356, 0.463], P(AUC>0.5) = 0.002
    - Sulfopin dG -7.58 +/- 0.28; 50 of 80 decoys score better
    - Juglone dG -7.76 +/- 0.51; 47 of 80 decoys score better
    - decoy dG median -8.91, sd 6.35; beating 80% of decoys needs dG < -11.38
    - one-frame vs 90-frame Spearman on the gate set is 0.283, so a SINGLE structure there is noise-dominated
    - the 90-frame mean is not: SEM 0.28-0.51 kcal/mol against a decoy spread of 6.35
    - verdict remains UNDERPOWERED: 2 actives < 3 floor, 2 chemotypes < 6 floor

Affects: `shared/mmgbsa_ensemble.py`, `scripts/run_mmgbsa_ensemble.py`, `decisions/D0031-class-matched-decoys-remove-the-apparent-covalent-enrichment.md`, `decisions/D0032-mmgbsa-gate-and-the-power-floor-on-negative-verdicts.md`

---

### D0037 — The junction dihedral was the sp3 analogue, and the reported dG was never an interaction energy

**:material-check: accepted** · `origin: adversary` · 2026-07-29

??? note "Evidence"
    - the frcmod header states every term is GAFF2's ss analogue; the whole DIHE block said 'from parm19' and used 1.00/3-fold/0deg for all five carbon types
    - GAFF2 gives X-c2-ss-X 2.200/2-fold/180, X-ca-ss-X 0.800/2-fold/180; only X-c3-ss-X matched what was in use
    - Juglone's built topology carried 2C-S-cc at PK 0.333, per 3, phase 0.0 -- the sp3 form on an sp2 attachment
    - affected 31 of 82 gate ligands: every cc/cd/ca/c2 attachment, i.e. exactly the ligands D0035 restored
    - GAFF2 has NO generic X-cc-ss-X or X-cd-ss-X, only specific 4-atom terms; cc/cd use cd-cc-ss-ca (2.430/2-fold/0) and are the least certain of the five
    - link-atom cap prevents bonded-term cancellation: residual 9.77 +/- 11.38 kcal/mol against a decoy spread of 6.35
    - gate ROC-AUC: full potential 0.425, standard interaction energy 0.181, the residual ALONE 0.831
    - Juglone interaction energy +9.08 kcal/mol, i.e. unfavourable, for a compound with published covalent activity
    - the measurement-error propagation quoted in D0036 existed in no code in the repository

Affects: `data/params/cys_gaff2_junction_5.frcmod`, `shared/mmgbsa.py`, `shared/mmgbsa_ensemble.py`, `shared/enrichment_gate.py`, `decisions/D0035-the-sp2-junction-gap-was-three-missing-angles.md`, `decisions/D0036-ensemble-mmgbsa-is-precise-and-still-below-chance.md`

---

### D0038 — The two solvent models disagree; but the dissociation I blamed on water was a single-trajectory artefact

**partially_withdrawn** · `origin: implementation` · 2026-07-30

??? note "Evidence"
    - 48 of 48 non-covalent candidates (T_1, T_2) run for 10 ns in explicit TIP3P, 0 failures
    - WITHDRAWN: t1_8a3f4861ac34 9.00 nm was NOT reproducible; implicit re-run gave 1.75 nm (engaged 0.51)
    - WITHDRAWN: t1_bd563e94c862 7.30 nm was NOT reproducible; implicit re-run gave 0.59 nm (engaged 0.91)
    - run-to-run divergence under the SAME model: 5.1x and 12.5x on mean ligand RMSD for those two
    - STANDS: spearman(implicit, explicit) = -0.102 run 1 and -0.144 run 2, i.e. two independent implicit runs both uncorrelated with explicit
    - Spearman(implicit RMSD, explicit RMSD) = -0.102 across 47 paired candidates
    - both models flag a similar COUNT as leaving (4 vs 5 of 47) but they are different candidates
    - candidates stable under GB drift furthest in water: t2_bc8a4b62eb0e 1.78 -> 4.86 nm, t1_c1ec9e35dba7 0.63 -> 4.66 nm
    - GROMACS 2026.3 CUDA sees all 8 A100s at ~740 ns/day; the shared OpenCL build refuses NVIDIA devices entirely

Affects: `shared/gromacs_explicit.py`, `shared/gromacs_analysis.py`, `scripts/run_gromacs_explicit.py`, `scripts/merge_gromacs_results.py`, `decisions/D0036-ensemble-mmgbsa-is-precise-and-still-below-chance.md`

---

## T_2 ATRA neighborhood

### D0018 — CReM fragment DB — ChEMBL33 SA<=2 primary, Enamine secondary, radius 2

**:material-check: accepted** · `origin: user` · 2026-07-27

T_2 was blocked on "choosing the fragment-DB radius". That framing was wrong and
worth correcting: the published CReM databases differ by **source**, **synthetic
-accessibility filter** and **minimum fragment frequency**. The context
**radius is a runtime argument** to `grow_mol`/`mutate_mol`, not baked into the
file. Two separate decisions were hiding inside one.

**Decision.** **Database — both staged, ChEMBL33 first.**

- **Primary: `chembl33_sa2_f5`** (281 MB). `f5` matches the spec's stated
  `--min-freq 5` default exactly. `sa2` restricts to synthetically accessible
  fragments, which serves condition (iii) *at enumeration time* rather than
  discovering unmakeable derivatives after labelling 10^4-10^5 of them.
- **Secondary: `enamine2025_sa2_f5`** (777 MB), staged now and unused for the
  first run. Enamine stock biases the neighbourhood toward fragments that can
  actually be bought — a different scientific stance from ChEMBL's "what has
  been published". Having it hash-pinned means that comparison is later a config
  change, not a fresh acquisition.

**Radius: 2** — revised from 3 on evidence, before any run.

Radius 3 was chosen first, on the reasoning that a larger radius demands more
surrounding context and so yields fewer but more chemically sensible
replacements. A smoke test against the actual seed refuted it:

| radius | ATRA mutate | ATRA grow |
|---|---|---|
| 1 | 45 | 50 |
| **2** | **43** | **38** |
| 3 | **0** | **0** |
| 4 | 0 | 0 |
| 5 | 0 | 0 |

**Radius 3 produces nothing at all for ATRA.** A control rules out a broken
setup: benzoic acid at radius 3 gives 47 mutations from the same database. The
cause is the seed itself — ATRA is a conjugated polyene with methyl branches, and
requiring three bonds of matching context finds no precedent for it anywhere in
ChEMBL33.

So radius 3 is not "conservative" here, it is inoperable: T_2 would have
enumerated an empty frontier and reported success. Radius 2 retains more context
than radius 1 while still being productive, and is the operative choice.

**Consequences.** The radius is the parameter that actually defines the neighbourhood, so changing
it changes what "degree-1 derivative of ATRA" means and invalidates comparison
with any earlier run.

**A general lesson for retargeting:** the usable radius is a property of the
SEED, not of the method. Any new seed needs this smoke test before its radius is
pinned — an unusual scaffold can silently produce an empty neighbourhood at a
radius that works fine for ordinary drug-like molecules. It is recorded in the T_2 config and pinned into each
run's manifest rather than passed ad hoc.

Running the Enamine database is a documented extension: swap the source in the
T_2 config, re-run, and compare. Both are hash-pinned so the comparison is
between two known inputs.

**Verification gap, stated plainly:** the file list and the SA/frequency
semantics come from the CReM download page. The radius guidance is a reading of
what the parameter does, not a quoted recommendation — the docs page gave none.
Worth checking against crem.readthedocs.io before any published result rests
on it.

??? note "Evidence"
    - published CReM DBs differ by SOURCE, SA filter and min fragment frequency — NOT by radius
    - radius is a runtime argument to grow_mol/mutate_mol, not a property of the file
    - chembl33_sa2_f5 is 281 MB; enamine2025_sa2_f5 is 777 MB
    - f5 matches the spec's stated --min-freq 5 default
    - the DB contains radius1..radius5 tables — all radii ship in one file
    - ATRA at radius 3, 4 and 5 yields ZERO mutations and ZERO grows
    - ATRA at radius 2 yields 43 mutations / 38 grows; radius 1 yields 45 / 50
    - control: benzoic acid at radius 3 yields 47 mutations — so radius 3 works in general, just not for this seed
    - CReM docstring: 'radius: radius of context which will be considered for replacement. Default: 3.'
    - benzoic acid replacements by radius: 1449 / 593 / 584 / 530 / 526 for r=1..5 — MONOTONICALLY DECREASING
    - benzoic acid by max_size at fixed radius 2: 40 / 204 / 593 / 593 for max_size=2/4/8/12 — INCREASING
    - max_size defaults to 10 and was never set; ATRA has 22 heavy atoms

Affects: `config/sources.yaml`, `config/approaches/t2_atra_crem.yaml`

Runbook: `docs/runbooks/adding_a_source.md`

---

## T_4 combinatorial

### D0004 — Build T_4's library fresh rather than reusing the prior run

**:material-check: accepted** · `origin: user` · 2026-07-27

A 7,104-member combinatorial library from a prior run existed and was of good quality (unique SMILES, core-verified, built on the graph engine). Reusing it would have saved real work.

**Decision.** Do not reuse it. Build T_4's warhead and R-group libraries fresh, grounded in the frozen reference set. The inherited library_size: 7104 was removed from gates.yaml; library size is now a derived, pinned output.

**Consequences.** Loses the prior work but gains a library that can actually satisfy the spec's grounding requirement, and removes a dependency on another group's directory tree (owner hemam, group ssmd-u-biodatsci-otherslab) that could move or vanish.

??? note "Evidence"
    - prior library: /data/lab_vm/refined/pin1_acr_screen, owner hemam, built 2026-07-21
    - reference set assembled 2026-07-26, five days later
    - so the prior library could not have been grounded in the reference set
    - prior library was 7104 rows = 16 warhead classes x 444 R-groups

Affects: `config/gates.yaml`, `data/reference/warhead_classes_2.csv`

---

### D0005 — Anchor the reactivity window on measured kinetics, and never rank by LUMO

**:material-check: accepted** · `origin: implementation` · 2026-07-27

The spec bounds T_4's reactivity window using computed LUMO. Figure 5C of Reddi 2023 supplies MEASURED second-order rate constants for Sulfopin and 4a-4g, which prompted checking whether reactivity predicts engagement.

**Decision.** The window is anchored on measured kinetics where available, with computed LUMO used to place NEW warheads on that calibrated scale rather than as the source of truth. T_4 must not rank candidates by LUMO: within the precedented range, electrophilicity does not predict engagement - recognition does.

**Consequences.** The reactivity window is a SAFETY filter for condition (ii), not a potency signal. Ranking by LUMO would have been an easy and invisible mistake. The kinetics values are figure-digitized and flagged as such; exact values need the SI tables.

??? note "Evidence"
    - Pearson r = 0.396 between intrinsic k and Pin1 labeling across 8 compounds
    - k spans 13.6x: 0.005 to 0.068 M-1 s-1
    - 4e: near-lowest k (0.007) but highest labeling (97%)
    - 4a: among highest k (0.030) but lowest labeling (17%)
    - values digitized from Reddi 2023 Figure 5C, ~1 significant figure

Affects: `data/reference/pin1_reactivity_kinetics_1.csv`, `config/gates.yaml`

Runbook: `docs/runbooks/resolving_unverified_structures.md`

---

### D0006 — Tier warhead structures by evidence and default to VERIFIED only

**:material-check: accepted** · `origin: implementation` · 2026-07-27

The warhead set must be data so the choreography can go wide later by adding rows. But breadth without provenance is how a library ends up half-dead, and a status field nothing enforces is just a comment.

**Decision.** warhead_classes_*.csv carries a structure_status per class; warhead_library.enumerable() defaults to VERIFIED only and logs a warning naming every class a caller widens to. window_anchor_classes() separately admits VERIFIED_CLASS_ONLY, because a window needs the chemotype, not the compound.

**Consequences.** T_4 currently has 3 enumerable classes. BDHI and naphthoquinone are blocked on attachment-regiochemistry design, which is a chemist's call and not a coding task. Going wide later is a CSV edit plus an explicit widening argument.

??? note "Evidence"
    - 4 tiers: VERIFIED, VERIFIED_CLASS_ONLY, NEEDS_DESIGN, UNVERIFIED
    - enumerable now: chloroacetamide, sulfamate_acetamide, sulfonate_acetamide
    - BDHI is VERIFIED_CLASS_ONLY (PubChem CID 21983498) - usable as a window anchor, not enumerable
    - naphthoquinone is NEEDS_DESIGN - juglone and KPT-6566 are intact quinones yielding no attachable fragment
    - prior run: 6 of 16 warhead classes collapsed to inert amides once attached

Affects: `data/reference/warhead_classes_2.csv`, `shared/warhead_library.py`

Runbook: `docs/runbooks/resolving_unverified_structures.md`

---

### D0013 — Enumerate competing attachment regiochemistries rather than choosing one

**:material-check: accepted** · `origin: user` · 2026-07-27

Two warhead classes had verified chemotypes but no established attachment
regiochemistry, and the honest answer from the PI was "I do not know". That is a
real state, not an oversight — the literature does not settle it for a sulfolane
core, and guessing would put an arbitrary choice underneath every T_4 result.

The constraint is tighter than it first looks: the reactive atom must stay free.
That reduces an open design question to four concrete candidates.

**Decision.** Enumerate **all four** as separate warhead classes and let the gates decide.
T_4 opts into `DESIGNED_UNTESTED` alongside `VERIFIED`.

The discriminating evidence costs nothing extra, because it comes from steps
already in the pipeline:

- **step 5b** — is the attached warhead still a genuine electrophile of its
  class, or did coupling kill it? A regiochemistry that fails this is refuted,
  not merely disfavoured.
- **step 6** — can Cys113 SG reach the reactive atom with the core in the way?
  A blocked approach shows up directly as poor covalent docking geometry.
- **step 7** — is its LUMO inside the window bounded by real actives?

Competing regiochemistries are reported **separately** through to ranking. Which
attachment works is a finding worth carrying, not an implementation detail to
collapse.

**Consequences.** Eight enumerable classes instead of four, so the enumerated library doubles.
That is affordable: the alert gate and 5b run before covalent docking, which is
the throughput wall, and MM-GBSA is still capped at per-class representatives.

DESIGNED_UNTESTED classes remain barred from anchoring the reactivity window
(control B5) however well they dock — surviving a gate is not the same as being
a validated Pin1 active.

If both members of a pair survive every gate, that is also an answer: the
attachment is not the discriminating variable, and both should carry forward.

??? note "Evidence"
    - BDHI: C3 bears the Br and is the Cys attack site, so attachment is limited to C4 or C5
    - 1,4-naphthoquinone: C2/C3 are the Michael acceptor positions, so attachment is C2 or the benzo ring
    - all four candidates parse and match their mechanism SMARTS
    - the discriminating evidence (5b validity, docking geometry, LUMO) comes from steps already in the pipeline

Affects: `config/approaches/t4_combinatorial.yaml`, `data/reference/warhead_classes_3.csv`, `shared/warhead_library.py`

---

### D0019 — Reactivity window flags classes, it does not reject them

**:material-check: accepted** · `origin: user` · 2026-07-27

The reactivity window as first built excluded nothing. Anchored on all six
verified covalent actives, it spanned ~4 eV — bounded at one end by mild
chloroacetamides and at the other by highly reactive quinones — and every one of
the nine warhead classes fell inside. A filter that does not filter, reported as
though everything had passed a safety check.

Two corrections were approved: drop the promiscuous quinone anchors (juglone and
KPT-6566 both carry `promiscuity_flag = y`, so bounding a SAFETY window with them
admits reactivity nobody would accept in a lead), and calibrate against the
measured rate constants rather than computed LUMO alone.

Both were applied. The window narrows from 3.96 eV to 1.61 eV and starts
discriminating — and the measured kinetics barely move the bounds set by the
clean anchors, which is reassuring: computed and measured chemistry agree on
where the range sits.

But it now excludes four of nine classes, **including acrylamide — the warhead
chosen by the PI for T_3**.

**Decision.** **The window FLAGS classes; it does not reject them.** Candidates outside the
window carry `reactivity_flag = OUTSIDE_WINDOW` and proceed to docking with
their evidence collected. `rejected_at` is not set.

Three reasons the exclusion is not strong enough to act on as a veto:

1. **The window is chemotype-narrow by construction.** Every clean anchor and
   every kinetics compound is chloroacetamide, sulfamate or sulfonate — one
   chemical family. A window built from one family will exclude other families
   almost by definition, whether or not they are genuinely unsafe. This is the
   "chloroacetamide-centric" caveat from the reference provenance, now biting in
   the opposite direction from before.
2. **LUMO is a weak proxy and we measured how weak.** D0005 found computed
   reactivity and Pin1 labelling correlate at r = 0.396. The window answers
   condition (ii) — is this electrophile in a precedented safety range — not
   whether the warhead will work.
3. **The margin is 0.106 eV**, against a 0.5 eV tolerance I chose myself. That
   is not a robust exclusion, and widening the tolerance after seeing which
   class it would readmit would be choosing the answer.

This is D0012 applied to a second gate: report evidence strength, let the human
adjudicate. The docking data for the flagged classes will exist regardless, so
the flag costs nothing and preserves the option.

**Consequences.** T_4 keeps all nine classes through docking. The GUI must display
`reactivity_flag` beside those candidates — a flagged candidate presented
without it would imply a safety assessment it did not pass.

Note this does not bear on T_3 directly: T_3's condition (ii) rests on the
up-front expert warhead choice, not on this window. It is informative about that
choice, not a contradiction of it.

If a flagged class subsequently ranks well on docking and MM-GBSA, that is
exactly the case worth a chemist's attention — a warhead outside the
precedented range but performing — and the honest resolution is measured
kinetics on the specific compound, not a recomputed LUMO.

??? note "Evidence"
    - window with ALL anchors spans 3.96 eV and excludes ZERO of 9 classes — a filter that does not filter
    - dropping the 2 promiscuous quinone anchors narrows it to 1.58 eV and excludes 4 of 9
    - adding the 8 measured-kinetics compounds barely moves it (1.61 eV, same 4 excluded)
    - excluded: acrylamide, naphthoquinone_c2, naphthoquinone_benzo, snar_chloroazine
    - acrylamide misses by 0.106 eV against a hand-chosen 0.5 eV tolerance
    - every clean anchor and every kinetics compound is chloroacetamide/sulfamate/sulfonate — ONE chemical family
    - D0005: computed reactivity vs measured Pin1 labelling correlate at only r = 0.396

Affects: `approaches/t4_combinatorial/02_reactivity_triage.py`, `config/approaches/t4_combinatorial.yaml`, `integration/app/DECISIONS_TAB_SPEC.md`

---

### D0020 — T_4 ranks within warhead class, with a per-class quota, not globally

**:material-check: accepted** · `origin: implementation` · 2026-07-27

T_4 fixes the sulfopin core and varies the warhead and the R-group. The question
it exists to answer is therefore *which warhead chemistry works on this core* —
a comparison across chemotypes.

The obvious implementation is to sort all 1,683 docked survivors on
`affinity_kcal` and take the best. That would answer a different question badly.

gnina's affinity is comparable only among molecules docked the same way, and
these are not. Each class is docked against its own
`covalent_lig_atom_pattern`, so the search is constrained differently per class
and the resulting pose ensembles are not samples from a common distribution. On
top of that, a Vina-style score grows with heavy-atom count, and the classes
differ substantially in size — the naphthoquinones carry roughly twice the heavy
atoms of the acrylamide. A global sort would largely report which warhead is
biggest and greasiest, and would hand most of the shortlist to one chemotype.

**Decision.** **Rank within warhead class; each class contributes a fixed quota (3) to the
shortlist.**

Three supporting choices:

1. **Ligand efficiency is computed but advisory.** `LE = -affinity / HAC` is the
   standard size correction and reviewers will ask for it, so it is reported.
   It is *not* the rank metric. D0015 fixed the rank metric on `affinity_kcal`
   because the enrichment gate measured it on this target (ROC-AUC 0.815, EF1%
   16.7); nothing comparable has been measured for LE here. Substituting an
   unmeasured metric because it is conventional would discard the only
   calibration the choreography has.

2. **Classes with few successful docks are flagged, not dropped.** Below 20
   docks, "best in class" is most of the class and the rank is not selective.
   Such rows carry `rank_is_selective = False` and say so in
   `shortlist_reason`.

3. **`OUTSIDE_WINDOW` classes still contribute** (D0019). Excluding them at
   ranking time would silently re-impose the kinetics filter that D0019
   deliberately relaxed — a veto reintroduced through the back door of a
   selection step.

Nothing is rejected here. This stage stamps and orders.

**Consequences.** The shortlist is a designed cross-chemotype comparison, not a league table. It
is deliberately *not* the set of nine best-docking molecules, and the GUI must
not present it as such: each row carries `class_rank`, `class_n_docked`,
`class_percentile` and `shortlist_reason` so a reader can see it was chosen as
best-of-its-class rather than best-overall.

Cross-class affinity comparison remains unlicensed downstream. If the
integration phase wants to rank warhead chemistries against each other, the
evidence for that is MM-GBSA on the true covalent adduct (step 9), which models
the bonded complex explicitly, not the docking score.

The per-class quota is in `config/approaches/t4_combinatorial.yaml`, so widening
the shortlist is a config change with a manifest record, not an edit to the
ranking code.

??? note "Evidence"
    - each warhead class docks against a DIFFERENT covalent_lig_atom_pattern, so pose ensembles are not drawn from one distribution
    - warhead heavy-atom count varies across classes (20-44 in the shortlist fixture) and a Vina-style score tracks size
    - on a synthetic frame where one class has the best raw affinities, a global top-9 excludes an entire class; the quota shortlist represents all three
    - D0015 fixed affinity_kcal as the rank metric after measuring it on 6 actives + 294 decoys; no equivalent measurement exists for ligand efficiency

Affects: `approaches/t4_combinatorial/04_rank_within_class.py`, `config/approaches/t4_combinatorial.yaml`, `tests/test_rank_within_class.py`, `integration/app/DECISIONS_TAB_SPEC.md`

---

### D0021 — BDHI and naphthoquinone attachment points resolved by paired docking

**:material-history: superseded** · `origin: implementation` · 2026-07-27 · superseded by **D0024**

Two chemotypes entered T_4 with a verified warhead class but an *untested*
attachment point, and the PI's instruction was explicit: "For BDHI and
Naphthoquinone, can we just try all the four cases? I just don't know."

So all four were enumerated as separate classes, and
`config/approaches/t4_combinatorial.yaml` declared in advance what would settle
it: `discriminated_by: [warhead_validity_gate_5b, covalent_docking_geometry,
lumo_window]`. All four passed 5b — coupling does not destroy the electrophile
in any of them — and the LUMO window did not separate the members of either
pair. That left the docking geometry, which is the physically apt test anyway:
the question is whether Cys113 can reach the reactive atom with the core in the
way.

Each regiochemistry was built against the same 187 R-groups, so the arms are
matched pair-for-pair on the only other varying factor. An R-group that docks
well does so in both arms and cancels.

**Decision.** **Carry `naphthoquinone_benzo` and `bdhi_c5`. Stamp `naphthoquinone_c2` and
`bdhi_c4` as superseded regiochemistries.** Neither is deleted — both keep their
rows and their docking evidence, per stamp-don't-delete.

Two endpoints were measured on each pair:

- **Pose success** (paired binary, McNemar). A Vina-style score at or above zero
  is not a weak binding energy; it means the search found no favourable pose.
- **Affinity** (Wilcoxon signed-rank, with matched-pairs rank-biserial effect
  size computed from the signed rank sums, not from win counts).

For **naphthoquinone** the result is not close. Attachment at C2 fails to pose
for 97% of R-groups; in 69 pairs the benzo isomer poses where C2 does not,
against 5 the other way. This matches the chemistry: C2 sits adjacent to the
Michael-acceptor positions, so the core is placed directly in the path of the
approaching cysteine. Verdict **STRONG**.

For **BDHI** the same direction holds with a modest margin: 59 discordant pairs
favouring C5 against 34 favouring C4, p = 0.012. Verdict **WEAK**. C4 places the
core adjacent to the reactive carbon, which the config anticipated might hinder
Cys approach; it does, but not decisively.

### Where the primary endpoint came from, and a knife edge

When most of the affinity column is censored, a rank test on it compares noise.
The verdict function therefore treats pose success as primary above 50% no-pose
in either arm. **That branch was written after seeing that C2 fails 97% of the
time** — the exact circumstance in which a threshold change can become a way of
choosing the answer, so it is recorded rather than buried. It does not change
any conclusion: both endpoints favour the same arm in both pairs, so only the
confidence label moves, never the winner.

BDHI lands on the knife edge — 50.3% no-pose, against a 50% threshold. Under the
affinity branch it would grade **UNDERPOWERED** instead of **WEAK**. The winner
is C5 either way; only the strength of the claim is sensitive to a threshold
that a handful of R-groups could flip. Treat the BDHI call as provisional.

**Consequences.** Step 9 (MM-GBSA on the true covalent adduct) is the most expensive stage in the
plan, and this halves the quinone and BDHI work it has to cover.

The reporting rule in `t4_combinatorial.yaml` — `keep_regiochemistries_separate:
true` — did its job and stays. Merging the arms would have averaged a working
geometry with a non-working one and reported a mediocre chemotype where there
are really one good and one bad attachment point.

The losing arms remain available: if MM-GBSA later contradicts the docking
geometry for the surviving arm, the alternative has not been thrown away.

None of this is a measured potency difference. It is a decision about which
attachment to carry forward, on docking geometry, which is the evidence the
config committed to in advance.

??? note "Evidence"
    - naphthoquinone: benzo poses in 69 pairs where c2 does not, vs 5 the other way, McNemar p = 1.8e-15
    - naphthoquinone_c2 finds NO favourable pose for 96.8% of its 187 R-groups; benzo 62.6%
    - naphthoquinone paired median difference +4.68 kcal/mol favouring benzo, Wilcoxon p = 2.5e-4
    - bdhi: c5 poses in 59 pairs where c4 does not, vs 34 the other way, McNemar p = 0.0124
    - bdhi no-pose 50.3% (c4) vs 36.9% (c5); paired median difference +1.09 kcal/mol favouring c5, Wilcoxon p = 6.6e-3
    - both endpoints — pose success and affinity — favour the same arm in both pairs
    - all four classes are paired on the SAME 187 R-groups, so the comparison is matched

Affects: `approaches/t4_combinatorial/05_regiochemistry_comparison.py`, `config/approaches/t4_combinatorial.yaml`, `data/reference/warhead_classes_3.csv`, `integration/app/DECISIONS_TAB_SPEC.md`

---

### D0022 — Covalent docking must use the adduct form, not the pre-reaction ligand

**:material-check: accepted** · `origin: implementation` · 2026-07-27

gnina's `--covalent_lig_atom_pattern` does exactly what its help text says: it
picks a ligand atom and bonds it to the receptor atom. It does not perform
reaction chemistry, and it does not remove a leaving group. The gnina
documentation says nothing either way about how the ligand should be prepared —
so the convention had to be established from the output, and the output is
unambiguous.

We supplied the intact, pre-reaction ligand. **Every one of the 1,683 docked
complexes therefore contains a molecule that does not exist**: the reactive
carbon sits 1.81 Å from Cys113's sulfur while still carrying the leaving group
it is supposed to have lost. For chloroacetamide that is an α-chloro thioether;
the real adduct has no chlorine on that carbon at all.

*Correction, same day.* An earlier version of this record asserted the docked
carbon was pentavalent. That is not established and should not be relied on.
gnina writes only polar hydrogens into its output — 1 of 24 survives for
chloroacetamide `t4_5e235921c8c0`, the imidazoline N–H — and it does not record
the S–C bond in the SDF at all, so the file cannot say whether a hydrogen was
displaced. Both readings are wrong molecules and both justify the same remedy,
but the specific pentavalence claim was not supported by the evidence cited for
it. What is certain is the retained leaving group and the geometry below.

For a leaving group on a flexible sp3 carbon this is survivable. Chloroacetamide
puts its chlorine a median 2.93 Å from SG and not one of its 187 poses clashes —
the CH2 rotates and the chlorine finds somewhere to go.

For a leaving group on a rigid sp2 or aromatic carbon it is not survivable. The
ring fixes the C–halogen direction, so bonding S to that carbon drives the
halogen into the sulfur:

| class | median halogen···SG | poses clashing (<2.5 Å) |
|---|---|---|
| `bdhi_c5` | 1.63 Å | 67.9% |
| `snar_chloroazine` | 2.38 Å | 55.6% |
| `bdhi_c4` | 2.55 Å | 47.1% |
| `chloroacetamide` | 2.93 Å | 0% |

The shortest contact observed is **0.89 Å** — two atoms occupying the same
space — and it belongs to the pose that scored best in the entire run
(`snar_chloroazine`, −9.16 kcal/mol).

### The three acetamides are the same molecule

Worse than the clashes, and easy to miss: `chloroacetamide`,
`sulfamate_acetamide` and `sulfonate_acetamide` are all SN2 displacements at the
same CH2. They differ only in what leaves. **Their adducts are one identical
molecule.** Verified directly: applying each displacement to the same R-group
gives a single distinct product SMILES.

So any difference between those three classes in the docking results is
*entirely* leaving-group artifact — and the differences were large, spanning
−5.47 to −2.87 kcal/mol in median. Sulfamate topped the table while carrying a
12-heavy-atom phantom group that is not present in anything that binds Pin1.
What actually distinguishes those three warheads is kinetics, which is what the
reactivity window measures, and which docking cannot see.

**Decision.** **Dock the adduct form: the post-reaction ligand with the leaving group removed
and the attachment atom left with an open valence for gnina to fill.**

This requires a second SMARTS per class. `reactive_atom_smarts` identifies the
reactive atom *in the pre-reaction molecule* and names the leaving group in
doing so (`[CH2][Cl]`, `[c]([Cl])[n]`); those patterns cannot match the adduct,
because the atom they key on is gone. Confirmed empirically: re-docking a
leaving-group-stripped ligand under the existing SMARTS produces no reactive
atom and no result. The library therefore needs a per-class adduct transform and
an `adduct_attachment_smarts` that matches the product.

Both forms stay in the library. The pre-reaction SMARTS is still what the 5b
validity gate and the warhead tests need — the question "is this a genuine
electrophile of its class" is a question about the *unreacted* warhead.

**Consequences.** **T_4's docking must be re-run.** The current per-class table is not a valid
cross-class comparison, and D0015's choice of rank metric is unaffected but the
values it ranks are not.

**D0021 is partially withdrawn.** The BDHI regiochemistry call rested on pose
success, and both BDHI arms clash heavily — with the winning arm clashing *more*
(67.9% vs 47.1%). That call cannot be separated from the artifact. The
naphthoquinone call stands: both arms are Michael acceptors carrying no leaving
group, so their comparison was between two ligands of identical composition.

**The enrichment gate (D0015, D0016) needs re-checking.** It was measured on
6 actives and 294 warhead-bearing decoys through this same protocol. Actives and
decoys were treated identically, so the comparison is internally consistent and
ROC-AUC is unlikely to move much — but "unlikely" is a prediction, and it should
be re-measured rather than assumed.

**After the fix, the three acetamide classes should converge.** Their adducts
are identical, so their docking scores must agree to within the search's own
noise. That gives a free internal control: if they do not converge after
re-docking, something in the protocol is still wrong. Their reactivity
differences remain real and remain the reactivity window's business.

Michael acceptors gain an H on the α-carbon in the true adduct, which docking
has been ignoring. That is a one-atom difference on a flexible position and is
minor next to a 12-atom sulfamate, but the adduct transform should handle it for
consistency rather than special-casing.

??? note "Evidence"
    - every docked reactive carbon sits 1.81 A from Cys113 SG while STILL CARRYING its leaving group — the docked molecule is not the adduct
    - clash rates (leaving-group atom within 2.5 A of SG): bdhi_c5 67.9%, snar_chloroazine 55.6%, bdhi_c4 47.1%, chloroacetamide 0%
    - shortest observed contact 0.89 A (Cl to SG, snar_chloroazine top-ranked pose)
    - chloroacetamide, sulfamate_acetamide and sulfonate_acetamide yield ONE identical adduct — their docking spread is leaving-group artifact
    - sulfamate carried a 12-heavy-atom leaving group through docking and ranked best by median (-5.47)
    - the Michael acceptors have no leaving group and are unaffected (0-3.2% clash)

Affects: `shared/covalent_protocol.py`, `shared/warhead_library.py`, `data/reference/warhead_classes_3.csv`, `approaches/t4_combinatorial/03_covalent_dock.py`, `approaches/t4_combinatorial/05_regiochemistry_comparison.py`, `decisions/D0021-regiochemistry-resolved-by-paired-docking.md`

---

### D0024 — Regiochemistry re-decided on adduct-form poses — bdhi_c4 and naphthoquinone_benzo

**:material-check: accepted** · `origin: implementation` · 2026-07-27

D0021 decided both regiochemistries on poses that turned out to have been docked
in the pre-reaction form (D0022). Its BDHI call was withdrawn; its naphthoquinone
call was argued to be unaffected. The re-dock on adduct-form ligands settles
both.

**Decision.** **Carry `bdhi_c4` and `naphthoquinone_benzo`.**

**BDHI reverses.** `bdhi_c4` now leads on median (−3.79 vs −2.87) where it
previously trailed (+0.01 vs −2.22). The reversal came through the loser: C4's
retained bromine had blocked half its poses, and that read as a geometric
failure of the C4 attachment when it was really a failure to remove a leaving
group. Verdict **UNDERPOWERED** — the affinity difference is real but modest
(p = 0.0035, rank-biserial −0.247) and the pose-success endpoint no longer
separates the arms (p = 0.14), because with the artifact gone most poses succeed
in both. This is a weak preference, and worth revisiting if MM-GBSA disagrees.

**Naphthoquinone is unchanged, STRONG.** Benzo over C2, 69 discordant pairs
against 6, p = 1.2e-14; C2 still finds no pose for 96% of R-groups. D0021's
withdrawal notice predicted exactly this on the grounds that Michael acceptors
carry no leaving group. The prediction holding is mild independent evidence the
diagnosis in D0022 was correct.

**Consequences.** The BDHI attachment is now a *weak* call rather than a decided one. It should be
carried as provisional and revisited at MM-GBSA, which models the bonded complex
explicitly rather than inferring geometry from a docking constraint.

The general lesson is recorded rather than left implicit: **a downstream
comparison inherits every defect of the poses it reads.** Both D0021 calls were
computed correctly, tested appropriately, and reported with calibrated
confidence — and one of them was still backwards, because the input was wrong in
a way no statistic on that input could reveal. The controls that caught it were
chemical (valence, interatomic distance), not statistical.

??? note "Evidence"
    - bdhi_c4 median -3.79 vs bdhi_c5 -2.87; paired median difference -0.58 kcal/mol, Wilcoxon p = 0.0035
    - bdhi pose success no longer separates the arms: 18 vs 29 discordant pairs, McNemar p = 0.14
    - bdhi no-pose collapsed once the bromine was removed: c4 50% -> 14%, c5 37% -> 20%
    - naphthoquinone_benzo over c2: 69 discordant pairs vs 6, McNemar p = 1.16e-14, unchanged from the pre-redock run
    - naphthoquinone_c2 still fails to pose for 96% of R-groups
    - convergence control passed: the three SN2 acetamides give identical best (-8.16) and median (-5.02)

Affects: `approaches/t4_combinatorial/05_regiochemistry_comparison.py`, `config/approaches/t4_combinatorial.yaml`, `data/reference/warhead_classes_4.csv`

---

### D0027 — Derive every Cys-ligand junction parameter from one GAFF2 source

**:material-check: accepted** · `origin: implementation` · 2026-07-28

The covalent bond joins ff19SB's protein sulfur `S` to a GAFF2 carbon, and no
force field covers that pair, so the junction terms are supplied by a hand-built
frcmod. Version 1 was written while debugging a chloroacetamide, whose
attachment carbon is sp3 (`c3`), and it covered that case only.

Six of the nine warhead classes do not attach through sp3 carbon. The Michael
acceptors and BDHI attach through sp2 (`c2`), the SNAr azine through aromatic
carbon (`ca`/`cc`). All 18 of their MM-GBSA builds failed on missing `S-c2`,
`S-cc` and `S-ca` terms. The nine that succeeded were exactly the three SN2
acetamide classes — the sp3 case v1 was built for.

Version 1 also carried `S-ca` *angles* without the `S-ca` *bond*, which is
incoherent on its own terms and would have failed the aromatic case even if the
class had been considered.

**Decision.** **Every junction term comes from GAFF2's own parameter for the same geometry,
with GAFF2's thioether sulfur type `ss` substituted for the protein's `S`.**

That substitution is an argument, not a convenience: once Cys113's SG has bonded
the ligand it sits between two carbons, `CB-S-C_lig`, which is precisely the
chemistry `ss` describes. GAFF2 supplies all of it — five bonds across the
carbon types that occur, and 136 terminal-S angles — so nothing is hand-chosen
and each line cites the `ss` form it came from.

**All 27 candidates are rebuilt under v2, including the 9 that already
succeeded.** Their v1 results are set aside as
`result_SUPERSEDED_junction_v1.json` rather than deleted. Keeping them would
have meant the acetamide classes rested on parm19-derived parameters while every
other class rested on GAFF2-derived ones — and while D0020 and D0023 already
forbid comparing dG across warhead classes, a parameter set that differs by
class makes even the within-class numbers rest on different footings depending
on which class you are in. One source, uniformly applied, costs a re-run and
removes the question.

**Consequences.** The junction remains the largest modelling assumption in the MM-GBSA module: it
is an approximation at the exact bond the calculation is about. What changed is
that it is now a *single, cited, uniform* approximation rather than an ad-hoc set
that silently covered one third of the classes.

The failure was loud, which is the only reason it was cheap. tleap refused to
build and named the missing atom types; had it substituted a default and
proceeded, 18 candidates would have carried quietly wrong energies into a
ranking.

If the junction is ever suspected of driving a result, the check is to re-run one
class with a deliberately perturbed junction parameter and confirm the
within-class ordering is unchanged. That has not been done.

??? note "Evidence"
    - 18 of 27 MM-GBSA builds failed; the 9 that succeeded were exactly the three SN2 acetamide classes
    - missing terms: S-c2 (9 candidates), S-cc (6), S-ca (3), plus 9 distinct angle types
    - junction v1 covered only sp3 carbon (c3); 6 of 9 warhead classes attach through sp2 or aromatic carbon
    - v1 also lacked the S-ca BOND while carrying S-ca angles — an incomplete set
    - GAFF2 provides every needed term under its thioether sulfur type `ss`: c2-ss 213.76/1.7842, ca-ss 213.48/1.7847, cc-ss 231.66/1.7538, plus 136 terminal-S angles

Affects: `data/params/cys_gaff2_junction_2.frcmod`, `shared/mmgbsa.py`, `approaches/t4_combinatorial/06_mmgbsa.py`

---

### D0029 — Three of T_4's nine warhead classes are one class after the reaction

**:material-check: accepted** · `origin: implementation` · 2026-07-28

**Decision.** 1. Post-reaction ranking, quotas and diversity counts group by
   **adduct class** (7), not warhead class (9).
2. Warhead class is retained on every row and carried into the GUI as
   the *route*, because the reactivity triage needs it and because
   three routes to one adduct is a result worth showing.
3. The shortlist is re-derived on adduct class once the current MM-GBSA
   run completes; the six freed slots go to molecules not yet scored.
   The in-flight results stay valid — they are correct energies for the
   molecules they name, merely redundant.
4. Docking should key its cache on `dock_id` so the 22% redundant work
   is not repeated on a re-run.

**Consequences.** **The class quota is triple-counting one chemotype.** Every class gets
three shortlist slots. The acetamide family therefore received nine
slots for three molecules, and the shortlist that MM-GBSA is scoring is
27 rows over 21 distinct molecules. Six of the 27 MM-GBSA runs
recompute an identical system.

**Any "diversity across warhead classes" claim overstates by two.**
T_4 covers **seven** post-reaction chemotypes, not nine. The
integration GUI must say seven.

**D0020 is not violated, it is confirmed.** D0020 forbids comparing
affinity across warhead classes. Here three classes are the same class,
so their scores are comparable — and are exactly equal. The one case
where cross-class comparison is legitimate is the case where there is
no cross-class comparison to make.

??? note "Evidence"
    - 1,683 docked rows collapse to 1,309 unique adducts (dock_id)
    - 100% of chloroacetamide, sulfamate_acetamide and sulfonate_acetamide share their adduct with the other two: 187 adducts, each reached by all three routes
    - the other six classes have 0% cross-class adduct sharing
    - affinity_kcal is identical to the last decimal across all 187 triples: max |chloro - sulfamate| = 0.0, max |chloro - sulfonate| = 0.0
    - the shortlist quota gave these three classes 9 slots for 3 unique molecules; the 27-candidate shortlist is 21 molecules
    - 374 of 1,683 covalent docks (22%) recomputed a pose already computed
    - 05_regiochemistry_comparison compared only bdhi_c4/c5 and naphthoquinone_c2/benzo, so no reported comparison is affected

Affects: `approaches/t4_combinatorial/04_rank_within_class.py`, `approaches/t4_combinatorial/06_mmgbsa.py`, `approaches/t4_combinatorial/05_regiochemistry_comparison.py`, `config/approaches/t4_combinatorial.yaml`, `integration/app/DECISIONS_TAB_SPEC.md`

---

## Integration

### D0008 — Files are the source of truth; the GUI is a view

**:material-check: accepted** · `origin: user` · 2026-07-27

Searching seven formats to answer 'why is the box 26 A' is untenable, and a single interface holding everything is genuinely more discoverable. But making the GUI the STORE would make reproduction depend on standing up Streamlit, and would lose git diff, code review, and bisect.

**Decision.** One store, many views, with files as the store. Decisions live in decisions/ as git-versioned records with machine-readable frontmatter. The GUI aggregates decisions + manifests + runbooks into a single Decisions pane, per approach and for the choreography.

**Consequences.** The GUI gains a real job beyond candidate display. If lab members ever need to AUTHOR decisions rather than read them, the GUI should write THROUGH to these files rather than into a database, so the repo stays authoritative.

??? note "Evidence"
    - provenance was scattered across 7 places: commit messages, .provenance.md, prep_log.json, manifest.json, runbooks, config comments, ready_to_delete.md
    - a published result must be reproducible from the repo alone, without running a web app
    - the GUI is built last (M5); M0-M4 need provenance now

Affects: `decisions/`, `shared/decisions.py`, `integration/app/`

---

