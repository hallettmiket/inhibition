# Runbook — selecting and verifying the shared receptor

**Use this when:** starting the choreography on a new target, or when someone
proposes changing the receptor mid-project (usually: don't).

**Why it matters:** every approach docks into the artifact derived from this one
file. If it is wrong, or if two approaches use different preparations of it,
every cross-approach comparison downstream is meaningless — and it will still
*look* fine, because docking scores come out either way.

---

## Procedure

### 1. Decide what the structure must contain

Write this down before searching, so you are not talked into a near-miss:

- the target protein, **human** sequence unless there is a stated reason;
- the **binding site of interest** actually resolved (not a construct where it
  is disordered);
- for a **covalent** campaign, the catalytic residue present and, ideally, a
  ligand bonded to it — the bond geometry defines the attack vector;
- a resolution good enough that side-chain positions are real, not modelled.
  Under ~2.0 Å is comfortable; over ~2.5 Å, be suspicious of rotamers.

### 2. Prefer a holo structure with a relevant ligand

A structure with a ligand in the pocket gives you the **docking box for free**:
the reference ligand's coordinates define the centre. An apo structure forces
you to define the box by eye or by pocket detection, which is a judgment call
you then have to defend.

If the ligand is covalent, the `LINK` record gives you the exact receptor atom
and the ligand attachment atom — which is precisely what a covalent docking
protocol needs.

### 3. Verify the entry yourself; do not trust the citation

Download it and check, in this order. Each of these has caught a real error at
some point:

```bash
grep -E "^HEADER|^TITLE|^REMARK   2 RESOLUTION" <PDB>.pdb
grep -c "<LIGAND_HET_CODE>" <PDB>.pdb
grep "^LINK" <PDB>.pdb
grep -E "^ATOM.*<RESNAME> <CHAIN> *<RESID>" <PDB>.pdb | head -3
```

- **TITLE** — is it the protein you think, and the right organism?
- **RESOLUTION** — matches what the paper claimed?
- **LIGAND** — present, and the HET code you expect?
- **LINK** — for covalent: does a bond genuinely exist between the catalytic
  residue and the ligand, and at what distance? A real covalent bond is
  ~1.7–1.9 Å. A "covalent" complex with a 3 Å contact is not bonded.
- **CATALYTIC RESIDUE** — present in the coordinates at the expected number.

### 4. Inventory every heteroatom before stripping anything

This is the step most likely to be done carelessly.

```bash
grep "^HETATM" <PDB>.pdb | awk '{print substr($0,18,3)}' | sort | uniq -c | sort -rn
```

Classify each code into:

| Class | Action | Examples |
|---|---|---|
| The reference ligand | strip, **but keep its coordinates** for the box | target-specific |
| Solvent / cryoprotectant / buffer | strip | HOH, GOL, EDO, PEG/PG4/1PE, SO4, DMS, MPD, ACT, TRS |
| Structural cofactor or metal | **KEEP** | ZN, MG, FAD, HEM, NAD |
| Unrecognized | **stop and look it up** | — |

Never strip by "everything that is HETATM" and never keep by "everything that
is not water." A stripped structural zinc silently changes the site; a retained
cryoprotectant silently occupies pocket volume.

`shared/receptor_prep.py` implements this as retain-and-report: anything not on
the strip list is **kept and counted** in `prep_log.json` under
`other_het_atoms_retained`. **If that number is not zero, go look at what it
is** before trusting the prepared receptor.

### 5. Sanity-check the box against the catalytic residue

After preparation, confirm the catalytic atom actually falls inside the box:

```
distance(box_centre, catalytic_atom) << box_size / 2
```

If the catalytic residue is near the box edge, the box is centred on the wrong
thing.

### 6. Decide whether one box is enough

If the reference ligand is covalent and small, a tight box around it is centred
on the **warhead sub-pocket**, which is right for covalent approaches and wrong
for non-covalent ones. Emit two boxes rather than compromising on one. Record
which approaches use which.

### 7. Pin it and never change it

Add to `config/sources.yaml`, stage with `shared.sources`, let the hash land in
`sources.lock.json`. Changing the receptor mid-choreography invalidates every
cross-approach comparison already computed.

---

## Worked example — Pin1, 2026-07-27

**Requirements.** Human Pin1; PPIase catalytic site; Cys113 present; covalent
campaign, so a covalently-bound ligand preferred; good resolution.

**Chosen: PDB 6VAJ.** Verification performed:

```
HEADER    ISOMERASE/ISOMERASE INHIBITOR           17-DEC-19   6VAJ
TITLE     CRYSTAL STRUCTURE ANALYSIS OF HUMAN PIN1
REMARK   2 RESOLUTION.    1.42 ANGSTROMS.
LINK         SG  CYS A 113                 C10 QT7 A 201     1555   1555  1.78
```

Human Pin1 ✓, 1.42 Å ✓, QT7 (sulfopin) present with 16 atoms ✓, and the LINK
record shows a genuine covalent bond `Cys113 SG — QT7 C10` at **1.78 Å** ✓ —
a real bond length, and it also hands over the exact ligand attachment atom.

**Heteroatom inventory:**

```
    132 HOH     -> strip (water)
     31 PG4     -> strip (tetraethylene glycol, cryoprotectant)
     16 QT7     -> strip, retain coords for box (reference ligand)
     10 SO4     -> strip (buffer)
```

PG4 was **not** on the initial strip list, so the first run retained 31 atoms
and reported them. Checking placed its nearest atom 22.65 Å from the box
centre — outside both boxes, so no result was affected — but it has no business
in a docking receptor and the strip list was corrected. **This is the
retain-and-report design working; the number was non-zero and got looked at.**

**Result:** 1,215 protein atoms kept; 16 ligand + 173 solvent removed; 0
unrecognized retained. Box centre = QT7 centroid `(-12.610, -34.943, 12.344)`;
Cys113 SG is 4.26 Å away, comfortably inside.

**Two boxes emitted** because QT7 is covalent at Cys113: a 20 Å covalent box for
T_3/T_4, and a 26 Å expanded box for T_1/T_2 so non-covalent search is not
biased toward the warhead sub-pocket.

---

## Failure modes seen or guarded against

| Failure | How it shows up | Guard |
|---|---|---|
| Cryoprotectant left in the receptor | nothing — docking succeeds, pocket is subtly occluded | retain-and-report, then inspect |
| Structural metal stripped | poses look plausible, chemistry is wrong | never strip unrecognized codes |
| Box centred on the wrong ligand copy | catalytic residue near box edge | step 5 distance check |
| Non-covalent search biased by a covalent box | candidates all crowd the warhead sub-pocket | emit two boxes |
| Receptor swapped mid-project | cross-approach numbers quietly incomparable | hash pin in `sources.lock.json` |
