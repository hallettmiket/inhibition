---
id: D0063
title: Reactive docking biases the search toward near-attack geometry, rather than filtering for it afterwards
date: 2026-08-05
status: proposed
approach: shared
decided_by: '@tt8804'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - docs/ranking_rationale.md
evidence:
  - 'Forli lab, Reactive Docking, JCIM 2023 (doi 10.1021/acs.jcim.3c00832): reactive parameters do not affect the length or complexity of the calculation'
  - 'AutoDock-GPU exposes --derivtype (derivative atom types) and --modpair (modified vdW pair params)'
  - "meeko.get_reactive_config(['C1'],['S4'],eps12=1.5,r12=2.0,...) -> derivtypes {'C':['C1'],'SA':['S4']}, modpairs {('C1','S4'): eps 1.5, r_eq 2.0, n 13, m 7}"
  - 'meeko MoleculePreparation takes reactive_smarts + reactive_smarts_idx; warhead_classes_10.csv already carries reactive_atom_smarts per class'
  - 'built locally: meeko 0.7.1, autogrid4 4.2.9, AutoDock-GPU CUDA build'
---

# Bias the search, not the score

## The change

`docs/ranking_rationale.md` proposed docking non-covalently and then **filtering**
poses for near-attack geometry. The established method does it the other way
round: **bias the sampling** toward reaction-competent geometry during docking,
with a restraint on the reactive centre that guides the warhead toward the
nucleophile while the rest of the ligand stays fully flexible.

That is strictly better here, and for a measured reason. Filtering can only find
a near-attack pose if the search produced one, and our search produces a correct
whole-molecule pose 41.5% of the time. Biasing addresses the sampling problem
rather than working around it.

**Biasing the SEARCH is legitimate; biasing the SCORE would not be.** We are
asking "can this molecule reach a reaction-competent geometry" — searching where
the answer must be is the question, not a thumb on the scale. Nothing about how
poses are then compared changes.

## How it works

Not a distance restraint bolted on. The reactive ligand atom is retyped to a
**derivative atom type** (`C1`), the receptor's nucleophile to another (`S4`),
and the **vdW pair potential between those two types alone** is replaced with one
whose well sits at the near-attack distance:

```
derivtypes  {'C': ['C1'], 'SA': ['S4']}
modpairs    {('C1','S4'): {'eps': 1.5, 'r_eq': 2.0, 'n': 13, 'm': 7}}
```

Every other interaction is untouched. Cost is unchanged — the published protocol
notes reactive parameters "do not affect the length nor the complexity of the
calculation".

## Why it fits what we already have

`warhead_classes_10.csv` carries **`reactive_atom_smarts` per class**, and Meeko's
`MoleculePreparation` takes `reactive_smarts` + `reactive_smarts_idx`. The
mapping is direct: each warhead class's SMARTS names the atom to retype. The
mechanism column then tells us which approach geometry to expect, which is the
part the pair potential does not encode.

## What it does not solve

The pair potential is **isotropic** — it rewards the reactive atoms being close,
not the *angle* of approach. An SN2 warhead pulled to 2 Å with its leaving group
pointing at the sulfur satisfies the potential and is chemically dead. So the
mechanism-specific angular criterion in `docs/ranking_rationale.md` stage 3 is
still required, as a filter on top of biased sampling rather than a replacement
for it.

## Installed

- `meeko` 0.7.1 — `~/.micromamba/envs/dwi_reactive`
- `autogrid4` 4.2.9 — same env
- **AutoDock-GPU**, CUDA build from `ccsb-scripps/AutoDock-GPU@89fd1c5` —
  `/data/lab_vm/modifiable/inhibition/autodock_gpu/bin/autodock_gpu_64wi`

AutoDock-GPU is a **different program from Vina-GPU** — AutoDock4's scoring
function and Lamarckian GA, not Vina's. So poses from it are not comparable to
anything measured so far, and it needs its own baseline on the benchmark before
any claim is made about it.
