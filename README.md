# Dance with Inhibition

Four independent computational approaches to finding an **inhibitor** of human
Pin1 (peptidyl-prolyl cis-trans isomerase NIMA-interacting 1), plus an
integration layer that presents their shortlists for a human to adjudicate.

This is a murmurent **choreography**: a problem statement many people attack,
each with their own **approach**. The pipeline code lives here; the generic
machinery for *describing* choreographies belongs in the murmurent repo.

- **Spec:** [issue #108](https://github.com/hallettmiket/murmurent/issues/108),
  Rev 3 (comment `5083543621`) + implementation plan (comment `5083611947`).
- **Authors:** Mike Hallett (hallett.mike.t@gmail.com), with Claude Code.
- **Target:** Pin1, catalytic Cys113. Shared receptor **PDB 6VAJ**
  (Pin1 + sulfopin/QT7, covalent at Cys113, 1.42 Å).

## The four approaches

| | Approach | Seed | Search | Covalent? | Inhibition proxy |
|---|---|---|---|---|---|
| **T_1** | de novo, structure-based generation | none — the pocket | DiffSBDD (pocket-conditioned diffusion) | non-covalent | weak (reversible occupancy) |
| **T_2** | derivative neighborhood of ATRA | ATRA | CReM, degree-bounded | non-covalent | weak (reversible occupancy) |
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

## Controls (these are not optional)

Four controls came out of an adversary audit of the spec and are wired in as
gates, not conventions. Each blocks a failure that is invisible as a crash:

1. **Docking-enrichment gate (B2).** Docking is not trusted to *rank* anything
   until it enriches known actives over property-matched decoys **on 6VAJ**. On
   failure `dock_score` is demoted to a displayed label — the choreography does
   not stall.
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

See the milestone table in the implementation plan. Current: **M0 in progress.**
