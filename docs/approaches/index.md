# The four approaches

Two complementary families. **T_1–T_2 generate broadly** (de novo, and the
derivative neighborhood of a weak non-covalent seed); **T_3–T_4 pursue covalent
inhibitors** off sulfopin under expert constraints.

| | Seed | Search | Rank metric | Direction |
|---|---|---|---|---|
| [T_1](t1.md) | none — the pocket | DiffSBDD | Vina `dock_score` | lower better |
| [T_2](t2.md) | ATRA | CReM degree-1 | `rank_score` (internal) | higher better |
| [T_3](t3.md) | sulfopin (core + warhead) | REINVENT `libinvent` | gnina `CNNaffinity` | **higher** better |
| [T_4](t4.md) | sulfopin (core) | combinatorial | MM-GBSA ΔG (within-approach) | lower better |

!!! danger "The directions genuinely differ"
    There is no global "more-negative-is-better" rule — that is false for the
    covalent tools. gnina `CNNaffinity` is dimensionless and **higher is
    better**; Vina affinity is kcal/mol and lower is better. Each output
    contract states its own direction, and the GUI reads it from
    `config/choreography.yaml`.

## What "protection" means

Encoding a substructure as a hard constraint the expansion cannot alter.
Nothing is physically installed or removed, unlike the synthetic-chemistry sense
of the word.

- **T_3** fixes core **and** warhead — one-point search, R-group only.
- **T_4** fixes only the core — two-point search, warhead **and** R-group. Which
  is exactly why T_4 needs a warhead-reactivity triage that T_3 does not.

## Shared obligations

Every approach must: key on canonical SMILES; use `shared/descriptors.py` for
the comparable physchem axes; compute novelty against the external set only;
stamp rather than delete rejected candidates; verify after expansion; and wait
for the enrichment-gate token before any dock-based ranking.

Each writes `Di.parquet` (full frame, stamped rejects retained), `Di_top10.csv`
(the GUI hand-off), and a `manifest.json`.
