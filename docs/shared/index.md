# Shared substrate

Everything in `shared/` is imported by every approach. It exists so that four
independent pipelines cannot quietly diverge on things that must match.

| Module | Guarantees |
|---|---|
| `smiles.py` | one canonicalization, so the join key is consistent |
| `descriptors.py` | one RDKit call per axis, so physchem values are comparable |
| `novelty.py` | novelty against the external set only, never the seed |
| `reference_set.py` | `UNVERIFIED` rows cannot reach the reactivity window |
| `warhead_library.py` | warhead classes tiered by evidence; `VERIFIED` only by default |
| `receptor_prep.py` | one prepared receptor, two boxes |
| `sources.py` | declarative acquisition, hash-pinned |
| `manifest.py` | every stage records what it consumed |
| `decisions.py` | the decision log, validated |

!!! warning "Build order"
    No approach may dock or rank until the receptor, the reference set, the
    covalent protocol, and the enrichment-gate token exist. The gate is what
    licenses docking to *rank* anything at all.

## The comparability chain

Four approaches produce numbers that mostly cannot be compared. The few that
can are comparable *because* of specific choices:

1. **Same receptor** — one prepared 6VAJ, hash-pinned ([D0001](../decisions/index.md)).
2. **Same descriptors** — one `descriptors.py`, so MW/cLogP/TPSA/QED/SAscore
   mean the same thing everywhere.
3. **Same novelty reference** — one frozen external set.
4. **Same covalent protocol** — T_3 and T_4 import one pinned gnina setup; if
   their protocol hashes differ, the within-covalent re-score is disabled
   rather than silently comparing incomparable numbers.

Break any of these and the integration phase loses the little commensurability
it has.
