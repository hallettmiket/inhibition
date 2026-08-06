---
id: D0062
title: Whole-molecule RMSD is the wrong endpoint for a covalent question, and the controls are not authoritative
date: 2026-08-05
status: proposed
approach: shared
decided_by: '@tt8804'
origin: user
supersedes: []
superseded_by: null
affects:
  - scripts/pose_selection_bench.py
  - docs/ranking_rationale.md
evidence:
  - '738 poses on 3IKD: spearman(whole-molecule RMSD, error in reactive-region placement) = +0.433'
  - 'reactive region within 1.0 A of the crystal in 59.8% of poses; whole-molecule RMSD <=2 A in 19.8%'
  - '52.7% of all poses place the reactive region correctly AND fail the RMSD test'
  - 'pocket_residues.json records derived_from 6VAJ.pdb, reference_ligand QT7, for_use_with 6VAJ_prepared.pdb — and was applied to 3IKD'
  - 'the 85-case benchmark contains 4 covalent PDB entries, and for those the benchmark selected the non-covalent ligand'
---

# The controls were inherited, not established

@tt8804, 2026-08-05: *"assume we are not sure of the validity of the known
dockers... heavily scrutinise each control you have been using. They are not
authoritative unless we determine them to be so."*

Auditing them produced one live error and one result that changes the endpoint.

## The endpoint was measuring the wrong thing

Every pose result so far — D0046's 5%, the 41.5% best-of-9, D0061's selection
rules — is scored on **whole-molecule RMSD to a crystal pose**. Our question is
whether a molecule can place its **warhead** to react. Those are not the same
question, and measured over 738 poses they are not even strongly related:

| | |
|---|---|
| spearman(whole-molecule RMSD, reactive-region placement error) | **+0.433** |
| reactive region within 1.0 Å of the crystal's | **59.8%** |
| whole-molecule RMSD ≤ 2 Å | 19.8% |
| **reactive region right, RMSD test failed** | **52.7% of all poses** |

**Docking places the reactive region correctly about three times more often than
the RMSD benchmark reports.** A pose at 3 Å whole-molecule RMSD with the reactive
carbon in the right place is *correct for our purpose* and is currently counted
as a failure.

So part of "docking does not work on this pocket" is an artefact of the endpoint.
Not all of it — the enrichment nulls (D0041, D0031) are about ranking, not pose —
but the pose half needs re-measuring against a criterion that matches the
question.

**Limits, stated so this is not over-read.** "Reactive region" here is the
ligand's closest heavy atom to the Cys113 sulfur: a one-dimensional proxy, and
two differently-oriented poses can share a minimum distance. Most benchmark
ligands are non-covalent and have no warhead. This establishes that the endpoint
matters; it does not establish the replacement.

## The pocket basis was 6VAJ's, applied to 3IKD

`pocket_residues.json` states its own provenance: `reference_ligand: QT7 A:201`,
`derived_from: 6VAJ.pdb`, `for_use_with: 6VAJ_prepared.pdb`. It is the 8 Å shell
around **sulfopin** in the **6VAJ conformer**, and every contact-profile rule in
D0061 used it against 3IKD. Residue numbering is shared so it is not nonsense,
but it is the wrong shell and must be re-derived from J9Z in 3IKD.

D0061's conclusion survives — the basis was consistent across rules, so "nothing
beat random" holds — but the contact rules were handicapped and deserve a re-run.

## The other controls, and what each is worth

| control | status |
|---|---|
| **crystal poses as truth** | each is that ligand's pose in **its own crystal's conformer**, superposed into 3IKD. For an induced-fit pocket some "failures" may be the benchmark being wrong, not the docking |
| **2.0 Å threshold** | conventional, not derived; weights all atoms equally so peripheral groups dominate |
| **Vina as sampler** | not authoritative — best-of-9 is 41.5%, i.e. it fails to find the answer more often than it finds it |
| **the 85-case benchmark** | almost entirely **non-covalent** ligands, used to validate a covalent method |
| **the 26 Å box** | inherited from `box_expanded`, chosen for 6VAJ's T_1/T_2 use, not derived for this purpose |

## What IS authoritative

**Experimentally measured reactivity.** The 17 verified covalent Cys113
complexes are molecules that demonstrably react — a fact about chemistry rather
than a docking convention. The reference set's measured potencies are the other.
Everything else is instrumentation adopted for convenience, and should be
labelled as such.

## What follows

1. Re-derive the pocket basis from J9Z in 3IKD; re-run D0061's rules on it.
2. Replace the endpoint with **reactive-atom placement plus approach angle**,
   defined per mechanism from `reactive_atom_smarts`, rather than whole-molecule
   RMSD.
3. Build the covalent validation set that the 85 cases cannot provide — the 17
   verified adducts, docked as free forms.
4. Re-state D0046 and the 41.5% as *whole-molecule pose reproduction*, which is
   what they measure, rather than as "docking works / does not work".
