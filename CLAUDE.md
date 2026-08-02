# inhibition — Claude Code context

`inhibition` ("Dance with Inhibition") is a murmurent **choreography**: four
independent computational approaches to finding an inhibitor of human **Pin1**
at catalytic **Cys113**, against **PDB 6VAJ**, plus an integration layer that
presents their shortlists for a human to adjudicate.

**The deliverable is the method, not the molecule.** Pin1 is the testbed. A
result about whether the choreography works beats a result about any individual
compound.

---

## Read these two before doing anything else

They are imported below so they are in context every session. They are not
background reading — the first tells you what is true right now, and the second
tells you how the numbers in this repo go wrong.

@docs/state_of_the_project.md

@docs/how_this_project_breaks.md

---

## What those two documents mean for how you work here

**Do not trust a number because it is populated and plausible.** Every
substantive bug found in this project has been the same bug — a value taken by
position, name, or default rather than by identity, failing silently because
both the right and the wrong candidate were populated and plausible. Before
using any column, ask what selected it and what that selection omits.

**Ranking is not validated on this target.** The docking-enrichment gate has
fired (D0041), and pose recovery is 5% in production (D0046). Every shortlist
carries `rank_validated = False`. Never describe a shortlist as evidence that
the molecules at the top bind — it is an ordering the pipeline produced.

**When you find a defect, fix the class and not the case**, write a decision
record in [`decisions/`](decisions/) including *why the wrong thing looked
right*, and add a guard that can actually fail.

## Where the rest of the context lives

- [`decisions/`](decisions/) — 51 records, the most valuable thing in the repo.
  Format and rules in [`decisions/README.md`](decisions/README.md).
- Open issues: **#4** (master plan), **#6** (open decisions + known defects),
  **#8** (questions out to the Lu lab), **#9** / **#10** (current direction,
  from a med-chemist review).
- [`README.md`](README.md) — how to *run* things, the four controls, setup.
- Spec: murmurent issue #108, Rev 3.

## Operational rules that bite

- **Data lives outside git.** `/data/lab_vm/immutable/inhibition/` is read-only;
  `/data/lab_vm/append_only/inhibition/<exp>/` is append-only and
  integer-versioned. Both are enforced by hooks, not convention. Retire
  superseded frames in `data/ready_to_delete.md` rather than deleting them.
- **Environments live outside the repo**, under `/data/lab_vm/envs/dwi_*`.
  Clone-and-run will not work without them. The shared CPU workhorse is
  `/data/lab_vm/envs/dwi_cheminf`.
- **Respect the compute budget** — `shared/compute.py` caps CPU workers (50).
  Run long jobs under `nice -n 19`. GPUs are shared with other users: check
  `nvidia-smi` for other people's processes before taking a card.
- **Reference files resolve by glob** (`shared/reference_set.latest_reference`).
  Never re-pin a version literal — a test walks the AST and will fail you.
- **Streamlit does not re-import helper modules.** Editing `curate.py` and
  clicking Rerun gives you the old module; restart the process.
