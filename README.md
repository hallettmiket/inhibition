# Dance with Inhibition

Four independent computational approaches to finding an **inhibitor** of human
Pin1 (peptidyl-prolyl cis-trans isomerase NIMA-interacting 1), plus an
integration layer that presents their shortlists for a human to adjudicate.

This is a murmurent **choreography**: a problem statement many people attack,
each with their own **approach**. The pipeline code lives here; the generic
machinery for *describing* choreographies belongs in the murmurent repo.

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
> Then `decisions/` (47 records — they record what was wrong and why it looked
> right), and the open issues: **#4** (the master plan), **#6** (open decisions
> and known defects), **#8** (questions put to the Lu lab — the chemistry
> judgement we cannot supply computationally), **#9** and **#10** (where the
> direction is currently being set, from a med-chemist review).

- **Spec:** [issue #108](https://github.com/hallettmiket/murmurent/issues/108),
  Rev 3 (comment `5083543621`) + implementation plan (comment `5083611947`).
- **Authors:** Mike Hallett (hallett.mike.t@gmail.com), with Claude Code.
- **Target:** Pin1, catalytic Cys113. Shared receptor **PDB 6VAJ**
  (Pin1 + sulfopin/QT7, covalent at Cys113, 1.42 Å).

## The four approaches

| | Approach | Seed | Search | Covalent? | Inhibition proxy |
|---|---|---|---|---|---|
| **T_1** | de novo, structure-based generation | none — the pocket | DiffSBDD (pocket-conditioned diffusion) | non-covalent | weak (reversible occupancy) |
| **T_2** | derivative neighborhood of a known binder | 5 seeds (ATRA, Liu-2024-C3, Potter-Astex, Du-Xu, Guo-Pfizer) | CReM, degree-bounded | non-covalent | weak (reversible occupancy) |
| **T_3** | single-warhead R-group decoration | sulfopin (core + warhead fixed) | REINVENT 4 `libinvent` | covalent | strong (covalent Cys113) |
| **T_4** | warhead × R-group combinatorial | sulfopin (core fixed) | combinatorial, 16 × 444 = 7,104 | covalent | strong (covalent Cys113) |

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

**All four arms have generated and been scored.** **13,863** candidates across
T_1–T_4 (T_1 4,803 · T_2 1,882 · T_3 5,396 · T_4 1,782), plus **42,588** from
the four new T_2 seeds and **15,653** from a degree-2 ATRA sample — about
**72,000 molecules generated**. Docking has completed on ~32,000 of them; the
Liu-2024-C3 and Potter-Astex pools and the degree-2 sample are still running or
queued.

The binding constraint is not generation. It is that **no scorer we have tested
discriminates on this target** — docking enrichment (D0041), docking *pose
recovery* (D0046: 5% in production against a 60–80% norm), ensemble MM-GBSA
(D0036) and MD residence (D0038, D0044) have each been measured and each
failed. That is the project's central finding so far, and it is why the current
plan (#4) spends its next phase on **measured inactives** rather than on more
generators.

For anything more specific than this paragraph, read
[`docs/state_of_the_project.md`](docs/state_of_the_project.md) — it is kept
current and this section deliberately is not.
