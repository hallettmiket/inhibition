# Runbook — resolving a reference compound whose structure is not in a database

**Use this when:** a compound you need as a reference or anchor has no PubChem
or ChEMBL entry — its structure exists only in a paper figure, an SI PDF, or
behind a paywall.

**Why it matters:** these compounds anchor controls. In this choreography the
covalent anchors bound T_4's reactivity window and the master set defines the
novelty axis. A wrong structure does not crash anything — it silently shifts a
control, and every downstream number inherits the error while looking normal.

**The rule that overrides convenience: never write a SMILES you cannot tie to a
record or verify by an independent check.** `UNVERIFIED` is a perfectly good
answer. A plausible guess is not.

---

## Procedure

### 1. Try the cheap sources first, in order

1. **PubChem by name** — `pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/<name>/property/SMILES/JSON`
2. **PubChem / ChEMBL by the paper's compound code** — often deposited.
3. **PubChem for the warhead or scaffold alone.** Frequently enough, and see
   step 2 — often it is all you actually need.
4. **Open-access mirrors** — bioRxiv preprint, PMC, institutional repository.
   A paywalled JACS paper often has a free preprint with the same figures.
5. **The paper itself**, if someone has access.

Record the identifier you found (CID, ChEMBL ID) — not just the SMILES.

### 2. Ask what you actually need: the compound, or the chemotype?

**This is the step that most often unblocks things, so do it before giving up.**

Look at how the structure is consumed downstream. If a quantity is a property
of a *substructure* rather than of the whole molecule, then the substructure's
structure is sufficient and may be far easier to verify.

Intrinsic electrophilicity is primarily a **warhead** property, not an R-group
property. So a reactivity window computed "per warhead class on a fixed
reference R-group" needs the **warhead chemotype**, which may sit in PubChem as
a simple heterocycle even when the paper's full fragment does not.

If you use a chemotype in place of a compound, say so explicitly in the status
field — it is a different claim, and a reader must be able to tell.

### 3. If you have a figure, transcribe carefully — then verify independently

Reading structures out of an image is legitimate and often the only route. It
is also error-prone in ways that do not announce themselves. So transcribe, and
then find **at least one independent number in the paper that your transcription
must reproduce.** Candidates, best first:

| Check | Uses | Strength |
|---|---|---|
| **Covalent adduct mass shift** | mass-spec panel | very strong — tests the whole structure and the reaction mechanism at once |
| **Molecular formula / exact mass** | any reported MW or formula | strong |
| **Elemental composition from HRMS** | SI characterization | strong |
| Degrees of unsaturation | formula | weak, but cheap |

**The adduct-mass check, generalized.** For a covalent inhibitor reacting with a
protein nucleophile:

```
predicted shift = MW(compound) - MW(leaving group as neutral)
```

Work out the leaving group from the mechanism (SN2 displacement loses the
leaving group; Michael addition loses nothing). Compare against the reported
shift. Agreement within ~1 Da over a 200–400 Da shift is a strong signal that
both the structure *and* your mechanistic assumption are right — it is two
checks for the price of one.

Do this in RDKit, not by hand:

```python
from rdkit import Chem
from rdkit.Chem import Descriptors
cmpd = Chem.MolFromSmiles("<transcribed>")
lg   = Chem.MolFromSmiles("<neutral leaving group>")
print(Descriptors.MolWt(cmpd) - Descriptors.MolWt(lg))   # vs reported shift
```

Also transcribe a **known** compound from the same figure (a control, a parent)
and confirm its formula matches the literature. That tests your reading of the
shared scaffold independently of the novel analogs.

### 4. Record status honestly, with a vocabulary that has teeth

Do not collapse everything to "verified"/"unverified". Use a tier that says
what is actually established:

| Status | Means |
|---|---|
| `VERIFIED` | traced to a public record, or derived from a verified structure, or confirmed by an independent numeric check |
| `VERIFIED_CLASS_ONLY` | the chemotype is verified; the specific literature compound is not |
| `NEEDS_DESIGN` | a real anchor exists but yields no usable fragment; someone must design the attachment |
| `UNVERIFIED` | unresolved. Refused by controls. |

Make the *code* enforce the tiers — `shared/reference_set.py` and
`shared/warhead_library.py` both refuse `UNVERIFIED` into the reactivity window.
A status field that nothing checks is a comment.

### 5. Write down what you did

In `data/reference/.provenance.md`: the source, the checks performed **and their
numbers**, and what remains unresolved. A future reader must be able to tell a
database lookup from a figure transcription.

---

## Worked example — Reddi 2023 sulfamate acetamide, 2026-07-27

**Problem.** `pin1_covalent_cys113_anchors_1.csv` carried the Reddi anchor as
`UNVERIFIED`. It was one of only two non-chloroacetamide chemotypes available,
so its absence left the reactivity window anchored on essentially one chemistry.

**Steps 1–2.** PubChem by name returned 404. The abstract
([PubMed 36738297](https://pubmed.ncbi.nlm.nih.gov/36738297/)) confirmed the
work but named no compound. It did reveal the Pin1 compounds are **4a–4g,
Figure 5** — a precise pointer, which is itself worth recording even when you
cannot yet act on it.

**Step 3.** Figure 5 was supplied directly. Panel A gives the structures.
Transcribed the series, then ran two independent checks:

1. **Control check.** Transcribed the shared Sulfopin scaffold → `C11H20ClNO3S`,
   MW 281.80 — the known Sulfopin formula. The neopentyl-N-(sulfolan-3-yl) core
   was read correctly.
2. **Adduct-mass check.** Panel B reports **+272 Da** for Pin1 + **4g**.
   Mechanism is SN2: `Cys-SH + R-CH2-OSO2NHPh -> Cys-S-CH2-R + HOSO2NHPh`, so
   the leaving group is N-phenylsulfamic acid (173.19). Predicted shift from the
   transcribed 4g = **271.4 Da**. Reported 272. **Match.**

Two checks, and the second validated the reaction mechanism as well as the
structure.

**Outcome.** Anchor promoted to `VERIFIED` as two compounds (4d, 4g). Verified
anchors 4 → 6. Enumerable warhead classes 1 → 3. A distinction the transcription
also surfaced: compound **4a** is a *sulfonate* (mesylate), not a sulfamate, and
was recorded as its own class rather than lumped in.

**A bonus finding.** Panel C gave measured rate constants, which turned out to
be more useful than the structures: across all 8 compounds, intrinsic reactivity
and Pin1 labeling correlate only weakly (**Pearson r = 0.396** over a 13.6×
range of k). That changed T_4's design — the reactivity window is a *safety*
filter, not a potency predictor, so T_4 must not rank by LUMO. **When you go to
a figure for one thing, read the other panels.**

**Still unresolved.** The Byun 2023 BDHI *fragment*. But step 2 applied: the
BDHI *warhead class* is [PubChem CID 21983498](https://pubchem.ncbi.nlm.nih.gov/compound/3-Bromo-4_5-dihydroisoxazole)
(`C1CON=C1Br`), and the window needs the class — so it is usable as
`VERIFIED_CLASS_ONLY` while the fragment stays `UNVERIFIED`.

---

## Failure modes

| Failure | Consequence | Guard |
|---|---|---|
| Inventing a plausible SMILES | control silently shifts; nothing errors | never write one you cannot tie to a record or a check |
| Transcribing a figure without verification | one wrong bond, no symptom | independent numeric check (step 3) |
| Treating a chemotype as the compound | overstates what is known | `VERIFIED_CLASS_ONLY` tier |
| Reading only the panel you came for | miss the data that changes the design | read all panels |
| Status field nothing enforces | `UNVERIFIED` rows reach the window anyway | loaders refuse them |
