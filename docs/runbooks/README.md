# Runbooks — the judgment steps, written down

Not everything in this choreography can or should be code. Choosing the right
PDB entry, deciding which heteroatoms are cryoprotectant and which are
structural, reading a compound out of a figure — these need judgment, and
having Claude Code (or a person) exercise that judgment is the point, not a
workaround.

**But the judgment has to leave a trace.** A one-shot decision made well and
then forgotten is worth very little the second time the problem appears with a
different structure, a different paper, or a different target. So each runbook
here records:

- **the decision procedure**, generalized — what to check, in what order, and
  what would make you reject the candidate;
- **the worked example** actually performed for Pin1, so the procedure is
  concrete rather than abstract;
- **the failure modes**, especially the silent ones.

Each is written to be followed by either a human or a future CC session. When
you hit one of these problems for a new target, read the runbook, follow it,
and **append your worked example** to it — the runbook improves by accretion.

| Runbook | Use it when |
|---|---|
| [`receptor_selection.md`](receptor_selection.md) | Choosing and verifying the shared receptor for a new target |
| [`resolving_unverified_structures.md`](resolving_unverified_structures.md) | A reference compound's structure is SI-only, paywalled, or in a figure |
| [`adding_a_source.md`](adding_a_source.md) | Bringing any new external file into the choreography |

## The division of labour

| Kind of step | Where it lives |
|---|---|
| Deterministic, parameterized, run repeatedly | code in `shared/`, config in `config/` |
| Acquisition of a declared file | `config/sources.yaml` + `shared/sources.py` |
| **Judgment, performed once per target** | **a runbook here** |
| What a specific run actually consumed | `manifest.json`, written by `shared/manifest.py` |

If you find yourself about to do something clever in a terminal that is not
covered by any of the above — that is the signal to write a fourth runbook.
