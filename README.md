# Dance with Inhibition

Four independent computational approaches to finding an **inhibitor** of human
Pin1 (peptidyl-prolyl cis-trans isomerase NIMA-interacting 1), plus an
integration layer that presents their shortlists for a human to adjudicate.

This is a murmurent **choreography**: a problem statement many people attack,
each with their own **approach**. The pipeline code lives here; the generic
machinery for *describing* choreographies belongs in the murmurent repo.

<!-- release-block:start -->
## Release 3.1.0 — in progress

> **Branch** `release/3.0.0` · **topic** `nac_v5` · **status** running
>
> | | |
> |---|---|
> | **Changelog** | [`CHANGELOG.md`](CHANGELOG.md) — no 3.0.0 number survives this release |
> | **Since the handoff** | [`docs/since_handoff.md`](docs/since_handoff.md) — what changed after @mhallet left, and why |
> | **Orientation** | [`docs/state_of_the_project.md`](docs/state_of_the_project.md) |
> | **How the project breaks** | [`docs/how_this_project_breaks.md`](docs/how_this_project_breaks.md) — read before writing code |
> | **Receptor change** | [`docs/branch_3ikd_receptor.md`](docs/branch_3ikd_receptor.md) — 6VAJ → 3IKD, and what it invalidated |
> | **Sweep length** | [`docs/sweep_length.md`](docs/sweep_length.md) — how 8 ns and 0.35 nm were chosen |
> | **GUI** | `python scripts/serve_reports.py` then `localhost:8931` |
> | **Pipeline** | `python scripts/pipeline.py status \| start <stage> \| stop <stage>` |
<!-- release-block:end -->

> ### New here? Read these two first
>
> 1. **[`docs/state_of_the_project.md`](docs/state_of_the_project.md)** — what
>    we are trying to find out, what is established, what is ruled out, and what
>    to do next. This README tells you how to *run* things; that tells you
>    *why*.
> 2. **[`docs/how_this_project_breaks.md`](docs/how_this_project_breaks.md)** —
>    every substantive bug found here has been the same bug. Read it before
>    trusting a number or writing code.
>
> Then `decisions/` (51 records — they record what was wrong and why it looked
> right), and the open issues: **#4** (the master plan), **#6** (open decisions
> and known defects), **#8** (questions put to the Lu lab — the chemistry
> judgement we cannot supply computationally), **#9** and **#10** (where the
> direction is currently being set, from a med-chemist review).

- **Spec:** [issue #108](https://github.com/hallettmiket/murmurent/issues/108),
  Rev 3 (comment `5083543621`) + implementation plan (comment `5083611947`).
- **Authors:** Mike Hallett (hallett.mike.t@gmail.com), with Claude Code.
- **Target:** Pin1, catalytic Cys113. Shared receptor **3IKD, curated by the
  project's medicinal chemist** — replaced 6VAJ on 2026-08-05 (D0059).

### Why the receptor changed, and what it bought

6VAJ is Pin1 co-crystallised with sulfopin, so its pocket is **induced-fit around
that ligand** and biases docking toward sulfopin-like chemistry. The curated 3IKD
is used *exactly as delivered* — no re-protonation, its 6 waters kept, nothing
rebuilt. Only the cognate ligand J9Z is removed (a receptor holding its own
ligand has an occupied pocket); its coordinates define the box. Cys113 arrives as
a **reactive thiol**, which the covalent arms require.

Re-running the pose-recovery benchmark on 82 crystal cases:

| is *any* pose in the top-k within 2 Å | 6VAJ | **3IKD** |
|---|---:|---:|
| top-1 — *what the pipeline carries* | 6.1% | **18.3%** |
| **best-of-9 — *the ceiling*** | 15.9% | **41.5%** |
| **random pick of the nine — *the floor*** | 5.3% | **19.8%** |

**The right pose is in the ensemble 2.6× more often — and the score is
indistinguishable from a coin flip at picking it.** Sampling works; selection is
the bottleneck. That is the finding the ranking design is built on.

*Confounded: 6VAJ was water-stripped and `reduce`-protonated, 3IKD keeps its
waters and the chemist's protonation, so part of the gain may be solvation rather
than conformer. The separating control has not been run.*

### How candidates are ranked

**On whether a molecule can orient to form the bond, not on how good the bond
would be.** Covalent inhibition is recognition then chemistry; the chemistry rate
belongs to the warhead *class*, so the molecule-to-molecule variation is almost
entirely in the recognition step. The docking score, covalent docking and
MM-GBSA are all excluded — each either measures the wrong thing or assumes the
bond already exists. Full reasoning, the pipeline, and how it fails:
**[`docs/ranking_rationale.md`](docs/ranking_rationale.md)**.

## The four approaches

| | Approach | Seed | Search | Covalent? | Inhibition proxy |
|---|---|---|---|---|---|
| **T_1** | de novo, structure-based generation | none — the pocket | DiffSBDD (pocket-conditioned diffusion) | non-covalent | weak (reversible occupancy) |
| **T_2** | derivative neighborhood of a known binder | 5 seeds (ATRA, Liu-2024-C3, Potter-Astex, Du-Xu, Guo-Pfizer) | CReM, degree-bounded | non-covalent | weak (reversible occupancy) |
| **T_3** | single-warhead R-group decoration | sulfopin (core + warhead fixed) | REINVENT 4 `libinvent` | covalent | strong (covalent Cys113) |
| **T_4** | warhead × R-group combinatorial | sulfopin (core fixed) | combinatorial, 9 × 198 = 1,782 | covalent | strong (covalent Cys113) |

Each approach emits a candidate frame `D^i` (rows = candidates, columns =
attributes, keyed on canonical SMILES) and hands its **top 10** to the
integration phase.

## What this deliberately does *not* do

**There is no authoritative cross-approach numeric ranking.** Cross-method
quantitative comparison is known a priori to be hard, and forcing it would be
the single easiest way to produce a confident wrong answer. Vina affinity
(kcal/mol, lower better) and gnina `CNNaffinity` (dimensionless, *higher*
better) are not the same axis; a non-covalent complex and a covalent adduct with
the leaving group removed are not the same physical quantity.

So the integration phase **presents rather than merges**: four shortlists side
by side, with structural convergence and the shared RDKit physicochemical axes
as the score-free cross-approach signals, an optional *within-stratum* re-score
clearly labelled as an aid, and the human making the call.

**Inhibition versus activation is not resolved computationally by any approach
here.** Catalytic-site occupancy is the proxy, and the proxy is not equally
strong across approaches — that asymmetry is displayed, not buried.

## Layout

```
config/       choreography.yaml, receptor.yaml, seeds.yaml, gates.yaml
shared/       BUILD FIRST — every approach imports this
data/reference/   the frozen reference set (the ONLY data in git)
approaches/t{1..4}_*/   per-approach stages
integration/app/  the artist's Streamlit GUI
scripts/      env setup, run wrappers (tmux-friendly)
tests/
```

Data lives outside git, under the governed roots:

- `/data/lab_vm/immutable/inhibition/` — read-only sources: 6VAJ, model weights,
  CReM fragment DB, decoys.
- `/data/lab_vm/append_only/inhibition/<exp>/` — derived, large,
  integer-versioned: frontiers, poses, trajectories, `Di.parquet`,
  `Di_top10.csv`.
- `/data/lab_vm/append_only/inhibition/00_outputs/<agent>/<topic>/` — analysis
  artefacts that belong to no single experiment: benchmark tables, retrosynthesis
  search trees, rendered poses, reading lists. Resolved by
  [`shared/outputs.py`](shared/outputs.py), which versions every write and
  resolves every read to the newest — the append-only tree needs both.

**Nothing derived goes in the repo.** There was an `outputs/` directory holding
1.1 MB of gzipped search trees and rendered pose HTML; it moved to `00_outputs/`
and is now `.gitignore`d. `data/` is for very small in-repo files only.

## Reading the results

Clone your own copy and run the GUI against the shared data — nothing is
written, so several people can read at once. Pick a port nobody else is using:

```bash
git clone https://github.com/hallettmiket/inhibition.git ~/inhibition
cd ~/inhibition
nice -n 19 /data/lab_vm/envs/dwi_gui/bin/python3.11 -m streamlit run \
  integration/app/app.py --server.port 8901 --server.address 127.0.0.1 \
  --server.headless true
```

Then from your own machine, forward that port and open
`http://localhost:8901`:

```bash
ssh -L 8901:127.0.0.1:8901 <you>@<host>
```

The interpreter and the data roots are shared, so there is no environment to
build. What you see is whatever the pipeline has written so far; panels for
stages that have not run are absent rather than empty.

## Controls (these are not optional)

Four controls came out of an adversary audit of the spec and are wired in as
gates, not conventions. Each blocks a failure that is invisible as a crash:

1. **Docking-enrichment gate (B2).** Docking is not trusted to *rank* anything
   until it enriches known actives over property-matched decoys **on 6VAJ**. On
   failure `dock_score` is demoted to a displayed label — the choreography does
   not stall. **This gate has fired.** Non-covalent enrichment is at chance
   (D0041: ROC-AUC 0.599, CI [0.311, 0.874], EF1% 0.0) and covalent enrichment
   goes to chance under class-matched decoys (D0031), so every shortlist is
   stamped `rank_validated = False`. Read a shortlist as an ordering the
   pipeline produced, **not** as evidence the molecules at the top bind.
2. **Frozen external reference set (B4/B5).** Novelty is `1 − max Tanimoto
   (ECFP4)` against the published binder set, **never against the seed** —
   seed-relative novelty is circular. The same set anchors T_4's reactivity
   window; the project's own leads are excluded as anchors.
3. **Warhead-validity gate (S2).** In the prior real run, 6 of 16 warhead
   classes collapsed to an inert amide once attached to the core. This gate
   catches that *before* covalent docking spends on it.
4. **Covalent-protocol parity (S3).** T_3 and T_4 import one pinned gnina setup
   unchanged; if their protocol hashes differ the within-covalent re-score is
   disabled rather than silently comparing incomparable numbers.

## Setup

```bash
bash scripts/setup_envs.sh all      # five isolated envs under /data/lab_vm/envs
```

Five envs, not one: DiffSBDD, REINVENT 4 and ADMET-AI each pull an incompatible
torch+CUDA build, and co-installing them is the most common way this stack
breaks. gnina is a binary/Docker call and needs no env.

```bash
conda activate /data/lab_vm/envs/dwi_cheminf   # the shared CPU workhorse
```

Long compute (CReM, DiffSBDD, REINVENT RL, covalent docking, xtb, MM-GBSA/MD)
runs under **tmux** and checkpoints, so a killed job resumes rather than
restarts.

## Status

*Current as of 2026-08-17. Numbers here are from the live run; the fuller
comparison against the handover document is in
[`docs/handover_delta.md`](docs/handover_delta.md).*

**The central finding has survived two changes of method.** @mhallet's handover
put it as: *we have candidates and no validated way to rank any of them*, over
four failed levels of theory. Since then the receptor changed, the ranking was
rebuilt around geometry instead of affinity, and **the new ranking measures as
unpredictive too** — so the finding now stands against six methods, on two
independent endpoints.

| method | result | record |
|---|---|---|
| Docking enrichment | AUC 0.599, EF1% 0.0 | D0041 |
| Docking pose recovery | 5% in production | D0046 |
| Ensemble MM-GBSA | below chance | D0036 |
| MD residence | not reproducible | D0038, D0044 |
| Contact-profile fit score | worse than chance; built and killed | D0057 |
| **Near-attack geometry ranking** | **ρ = +0.119, *p* = 0.33** vs sweep outcome | this run |

That last row is measured on 68 swept modes of the current run. `enrichment`
gives ρ = +0.033; `class_rank` gives ρ = −0.256 at *p* = 0.035 — the only one
that clears significance, and it points the **wrong way**. Median enrichment
among survivors is 6.03, among those that left 6.12.

**The caveat that must travel with it:** every one of those 68 cleared the
enrichment floor, so this measures discrimination *above* the floor, not whether
the floor excludes anything. The stratified pilot that would test the floor
([#71](https://github.com/hallettmiket/inhibition/issues/71)) has never run.

### Then and now

| | 2026-08-02 handover | 2026-08-17 |
|---|---|---|
| receptor | 6VAJ | **3IKD**, chemist-prepared (D0059) |
| arms in play | four (T_1–T_4) | **T_4 only** (D0081) |
| molecules | ~72,000 docked and ranked | **561** in scope, all screened |
| ranking basis | affinity, size-decorrelated | **near-attack geometry**, EB-shrunk |
| unit of selection | the molecule | **the binding mode** |
| pose handling | one representative | two-stage splitting, ≤5 sub-modes |
| MD | ad hoc, not reproducible | 3-stage cascade with measured gates |
| running compute | nothing | pipeline live on 3 GPUs |

### The current run

561 molecules screened → 4,432 binding modes → 2,019 ranked → **147 selected for
MD**. The 8 ns triage is in progress; 6 modes have held under 0.35 nm so far.
The first 100 ns run finished at **1.140 nm** — the best mode either run has
produced (0.323 nm and 84.8% attack-ready over 8 ns) moved 3.5× further over the
longer trajectory and missed the bar.

### What is open

- **Is the enrichment floor real?** It discards 3,700 of 4,432 modes on a
  parameter with no demonstrated relationship to outcome
  ([#71](https://github.com/hallettmiket/inhibition/issues/71)).
- **BPMD has no GPU** ([#72](https://github.com/hallettmiket/inhibition/issues/72)),
  and consumed 1.3 TB in 3.0.0 — 87% of the project's footprint — for no usable
  result.
- **103 files outside the pipeline path still hardcode the dataset root**
  ([#74](https://github.com/hallettmiket/inhibition/issues/74)).

