# Versioning the method

This project versions **the discovery method**, not a public API. SemVer is
written for libraries, so the three fields are given an explicit meaning here
rather than borrowed by analogy — otherwise "is this a major?" becomes a matter
of taste, and the number stops carrying information.

## The rule

The question a version answers is: **can I quote a number from an earlier
version alongside one from this version?**

| bump | meaning | test |
|---|---|---|
| **MAJOR** | previously reported numbers are **invalid** and must be re-measured | would re-running the old pipeline give a different answer to the *same* question? |
| **MINOR** | new capability, metric, or workup. Existing numbers stay valid | does anything already reported change? If no → minor |
| **PATCH** | defect fix that corrects numbers **within an unchanged definition** | same metric, same inputs, previously wrong value |

The distinction that does the most work is between **invalidating a measurement**
and **invalidating an interpretation**. They are not the same bump:

- D0059 replaced the receptor. Every 6VAJ number is now a measurement of the
  wrong thing — 6VAJ and 3IKD place the pocket 48.6 Å apart. **The measurements
  died. That is a MAJOR.**
- D0071 showed the ranking metrics do not predict pose stability. The enrichment
  and consensus values are still correct measurements of what they measure; what
  died is the claim that they predict something physical. **An interpretation
  died. That is not a MAJOR.**

A metric **redefinition** is a MAJOR even without a receptor change, because the
old and new values share a name and are not comparable. D0068 identified exactly
this hazard: redefining `enrichment` from a whole-population fraction to a top-N
fraction would make every published enrichment value non-comparable while still
being called "enrichment". If 2.1.0 does that, it is **3.0.0**.

## Release names

Releases carry a name as well as a number, from a **copper-mineral alphabet**
(issue #27): every name is an IMA-recognised mineral with Cu in its formula,
taken in order.

| version | letter | mineral | formula | colour |
|---|---|---|---|---|
| 2.0.0 | **A** | **Azurite** | Cu₃(CO₃)₂(OH)₂ | deep azure blue |
| 2.1.0 | **B** | **Bornite** | Cu₅FeS₄ | iridescent peacock |
| 2.2.0 | **C** | **Chalcopyrite** | CuFeS₂ | brass yellow |

The letter advances **per release**, not per major — so the next release after
Chalcopyrite is Dioptase whether it is 2.3.0 or 3.0.0. A name is easier to say
in a meeting than "two point one", and unlike a number it cannot be confused
with a different quantity that also has dots in it.

**1.0.0 has no name.** It was assigned retroactively and never cut (below), so
giving it a letter would imply a release that never happened and would push every
real release one letter along.

## Applying it to what we have

**1.0.0** is assigned **retroactively** and is not a real release. It names
everything up to the handoff from @mhallet — a 6VAJ pipeline with no version
history of its own. It exists so that "which pipeline produced this?" has an
answer for pre-handoff numbers, not because a 1.0.0 was ever cut.

**2.0.0** is correct as a MAJOR: D0059 invalidated every number the project had.

But 2.0.0 **bundles what convention would have released separately**, because the
whole line was developed unreleased on one branch. Audited honestly, the branch
contains at least:

| would have been | change |
|---|---|
| MAJOR | D0059 — 3IKD replaces 6VAJ; every prior number invalid |
| MINOR | D0063/D0064 — reactive docking replaces dock-then-filter |
| MINOR | D0065 — the mechanism-specific near-attack criterion |
| MINOR | D0070 — pose consensus as a second component |
| MINOR | the elevation suite (tiers 1–4) |
| **PATCH** | **D0067 — BDHI scored with sp³ geometry at an sp² carbon.** Same metric, same inputs, 374 candidates previously wrong |
| PATCH | pose atoms keyed on the PDBQT name field; SNAr class silently deleted; GUI drew every pose into 6VAJ; non-idempotent pose export |
| — | D0061/D0062/D0066/D0068/D0069/D0071 — **findings, not code changes.** They change what we believe, not what the pipeline computes, and do not bump anything on their own |

That last row is the one most easily got wrong. **A decision record is not a
release.** D0071 is the most consequential thing in 2.0.0 and it bumps nothing,
because it did not change a single computed value.

Recording the bundle rather than back-dating a dozen tags: the tags would be
fiction, and the audit above carries the same information honestly.

## Going forward

Tag at the point work becomes quotable — when a number leaves the repo, into the
manuscript, a chemist's shortlist, or a report. Not on every merge.

`VERSION` holds the current version. `CHANGELOG.md` carries the history and
states, per entry, whether prior numbers survive it.
