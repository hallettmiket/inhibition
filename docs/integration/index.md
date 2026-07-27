# Integration — present, don't auto-rank

**Status: specified, not implemented.**

Cross-method quantitative comparison is hard and is **not attempted as an
authoritative ranking**. Each approach delivers its top 10, ranked by its own
internally-valid metric. Integration is a **presentation and human-decision
layer** — the artist's Streamlit GUI — showing the four shortlists (40
candidates) side by side.

## What it shows

### 1. Per-candidate dossier
2D depiction; an interactive 3D viewer for the predicted binding mode against
Pin1; the approach's own metrics **with their stated directions**; mechanism and
inhibition-proxy strength; and which of conditions (i)–(v) it actually resolved,
with coverage gaps shown rather than hidden.

### 2. Score-free cross-approach signals
- **Structural convergence** — cluster all 40 by scaffold and ECFP4 Tanimoto. A
  molecule surfaced independently by more than one approach is a soft
  cross-validation relying on **no** score commensurability. This is the most
  defensible cross-approach signal available.
- **Shared physicochemical axes** — MW, cLogP, TPSA, QED, SAscore, computed by
  the identical RDKit call, plotted across the pool.

### 3. Optional within-stratum re-score
Non-covalent (T_1+T_2) through one Vina protocol; covalent (T_3+T_4) through one
gnina protocol. **Two leaderboards, never one** — a comparison aid, never
authoritative, and cross-stratum ordering is not implied.

### 4. The human decides
The choreography surfaces and organizes the evidence; it does not output a
single winner. Honest limits are displayed: no authoritative cross-method
ranking, inhibition-versus-activation unresolved, no wet-lab ground truth.

## The Decisions pane

Beyond candidates, the GUI is the single place to answer "why is this the way it
is" — aggregating decision records, run manifests and runbooks. Spec:
`integration/app/DECISIONS_TAB_SPEC.md`.

The GUI **reads; it does not own**. Everything it shows is a rendering of files
in the repo ([D0008](../decisions/index.md)).
