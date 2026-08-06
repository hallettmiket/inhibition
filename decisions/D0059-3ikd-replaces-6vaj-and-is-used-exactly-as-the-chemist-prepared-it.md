---
id: D0059
title: 3IKD replaces 6VAJ as the receptor, used exactly as the chemist prepared it
date: 2026-08-05
status: proposed
approach: shared
decided_by: '@tt8804'
origin: user
supersedes: []
superseded_by: null
affects:
  - scripts/prepare_3ikd_receptor.py
  - config/receptor.yaml
  - shared/noncovalent_dock_run.py
evidence:
  - "chemist's recommendation, relayed by @tt8804 2026-08-05: use only the prepared 3IKD structure"
  - '6VAJ is Pin1 co-crystallised with sulfopin, so its pocket is induced-fit around that ligand (#14)'
  - 'measured on 82 crystal cases while building D0057: crystal pose ranked #1 in 0/82 cross-docked into 6VAJ vs 5/82 self-docked into each ligand own receptor; top-3 3.7% vs 22.0%'
  - 'uploaded file sha256 37cccd49b585c4bce1d1f8e6e8a68aa0178aae2f79f480271cac24e3ccacb68d, 1807 atoms, 888 hydrogens, 6 waters, J9Z present'
  - 'Cys113 SG at (13.385, 3.989, -2.040) carrying HG — a reactive thiol, which config/receptor.yaml requires'
  - 'J9Z centroid 3.57 A from the Cys113 sulfur'
---

# One receptor, taken as given

## The decision

**3IKD, exactly as `@tt8804`'s chemist prepared it, replaces 6VAJ.** The file is
authoritative: nothing is re-protonated, no waters are stripped, no residue or
rotamer is rebuilt.

## Why

6VAJ is Pin1 co-crystallised with **sulfopin**, so its pocket is in an
induced-fit state shaped around that ligand. Docking into it biases toward
sulfopin-like chemistry, which is the criticism raised in #14 and the reason the
chemist recommended the change.

Our own data gives an independent measurement of the transfer penalty. Building
the contact-profile fit score (D0057), the same test run two ways over 82 crystal
complexes:

| | crystal pose ranked #1 | top 3 |
|---|---|---|
| cross-docked into 6VAJ | 0 / 82 (0.0%) | 3.7% |
| self-docked into each ligand's own receptor | 5 / 82 (6.1%) | 22.0% |

A ligand's crystal pose transplanted into 6VAJ does not make 6VAJ-like contacts.

**The caveat this record must not lose: 3IKD is also a ligand-bound structure**
(cognate ligand J9Z, 2.0 A). Swapping 6VAJ for 3IKD trades sulfopin-induced fit
for J9Z-induced fit unless the chemist's preparation specifically addressed that,
and the file does not say. "Not induced-fit" and "induced-fit by a different
ligand" are different claims.

## What the file is

Verified on arrival rather than assumed:

| | |
|---|---|
| chain A, residues 51-163 | the PPIase domain; the WW domain is absent |
| **J9Z present**, 26 heavy atoms | the box comes from the entry's own reference ligand |
| **Cys113 SG at (13.385, 3.989, -2.040), with HG** | a **reactive thiol** — required, or T_3/T_4 have nothing to attack |
| Cys57 also present | the second ligandable cysteine |
| 888 hydrogens | already protonated |
| 6 waters | retained deliberately |

## The one change made, and why it is unavoidable

**J9Z is removed; its coordinates are kept.** A receptor still containing its own
cognate ligand has an occupied pocket, and every docked molecule would be scored
against a full site. The coordinates define the box. This is the identical rule
`config/receptor.yaml` states for 6VAJ's QT7: *"Stripped before docking, but its
coordinates DEFINE the box."*

CONECT records were dropped with it — they referenced the removed atoms.
Conversion to PDBQT is `obabel -xr` (rigid receptor), with **no `-p`**.

## What was deliberately NOT done

* **No `reduce -BUILD`.** The file arrives with 888 hydrogens assigned. 6VAJ went
  through `reduce` at pH 7.4; re-running it here would double-add or silently
  re-assign the chemist's choices.
* **Waters kept — all 6.** `6VAJ_prepared.pdb` has none. Pin1's site is
  water-mediated, which is cited in the FEP rule-out as a reason this pocket is
  hard, and six were retained on purpose. They are now part of the rigid
  receptor.
* No residue rebuilding, no rotamer changes.

Verified after conversion: 6 waters survived, SG is at identical coordinates, and
**HG survived as a polar `HD`** — so the thiol is intact in the same form 6VAJ
carries it.

## Two things left unresolved, recorded rather than assumed away

**1. 6VAJ and 3IKD are now prepared DIFFERENTLY.** 6VAJ: waters stripped,
protonated by `reduce` at pH 7.4. 3IKD: waters kept, chemist's protonation. So
when the redock benchmark is re-run and compared against D0046's 5% pose
recovery, **any difference conflates receptor with preparation.** Attributing it
to the receptor needs a control -- deposited 3IKD through 6VAJ's original path --
which is cheap and has not been run.

**2. The preparation pH is unknown.** `REMARK 200 PH: 8.00` in the file is the
**crystallisation** pH of the 2009 deposition, not the protonation target; the
file does not record what that was. Our ligands are prepared at pH 7.4
(`LIGAND_PREP_TAG`). His and Cys are exactly where a 7.4-vs-8.0 difference is not
automatically negligible, and Cys113 is the residue that matters. A question for
the chemist.

## What this invalidates

Every gate verdict and benchmark on record is a **6VAJ** number and does not
transfer:

| measurement | record |
|---|---|
| non-covalent enrichment ROC-AUC 0.599 | D0016 / D0041 |
| **pose recovery 5% production, 55% best-of-9** | **D0046** |
| covalent enrichment at chance | D0031 |
| size decorrelation \|rho\| <= 0.034 | D0049 |

`rank_shortlist.attach_gate` looks a verdict up by **(stratum, metric)**, and the
metric NAME does not change when the receptor does — so a 6VAJ verdict would
attach silently to 3IKD. D0051 makes an *unknown* metric fail closed; this is a
*known* metric on a different receptor, which is worse. **Re-run the gate before
anything ranks.**

Generation is unaffected: all 86,451 molecules are receptor-independent.

## First experiment

Re-run D0046's redock benchmark on 3IKD. The 83 crystal cases already exist and
cross-docking them is minutes of GPU. It answers the question that justifies the
switch — **does pose recovery beat 6VAJ's 5%?** — and either answer is worth
having before a week of BPMD is spent.
