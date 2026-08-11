# inhibition — Claude Code context

`inhibition` ("Dance with Inhibition") is a murmurent **choreography**: four
independent computational approaches to finding an inhibitor of human **Pin1**
at catalytic **Cys113**, plus an integration layer that
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

**The receptor is 3IKD, not 6VAJ** (D0059, 2026-08-05) — the chemist's
prepared structure, on branch `receptor/3ikd-chemist-prepared`. **Every 6VAJ
measurement is invalidated until re-run**, including D0046's 5% pose recovery.
Re-measured on 3IKD: top-1 2.4%, best-of-9 **41.5%** (was 15.9% on 6VAJ by the
same metric) — the search improved 2.6x and the score still cannot exploit it.

**Ranking is not validated on this target.** The docking-enrichment gate has
fired (D0041) and five independent levels of theory have now been measured and
failed. Every shortlist carries `rank_validated = False`. Never describe a
shortlist as evidence that the molecules at the top bind — it is an ordering
the pipeline produced.

**When you find a defect, fix the class and not the case**, write a decision
record in [`decisions/`](decisions/) including *why the wrong thing looked
right*, and add a guard that can actually fail.

## Where the rest of the context lives

- **[`docs/since_handoff.md`](docs/since_handoff.md) — what has changed since
  @mhallet left, newest first.** Read it before trusting a number from before
  2026-08-02: the receptor changed, five defects were producing plausible wrong
  output, and every 6VAJ measurement is invalidated.

- [`decisions/`](decisions/) — the most valuable thing in the repo. Format and
  rules in [`decisions/README.md`](decisions/README.md). `origin` is an
  **allowlist** (`adversary`/`implementation`/`spec`/`user`) and will refuse a
  value you invent — that is deliberate, per D0051.
- Open issues, consolidated 2026-08-04: **#12** (chemistry judgement, out to the
  Lu lab) and **#13** (every open technical problem, audited against the code).
  **#4** is the plan and reasoning of record. #2, #6, #8 and #11 are closed into
  those; read them for history, not for status.
- [`README.md`](README.md) — how to *run* things, the four controls, setup.
- Spec: murmurent issue #108, Rev 3.

## Operational rules that bite

- **Data lives outside git.** `/data/lab_vm/immutable/inhibition/` is read-only;
  `/data/lab_vm/append_only/inhibition/<exp>/` is append-only and
  integer-versioned. **Both are a DISCIPLINE, not enforcement** — this line used
  to claim "enforced by hooks, not convention" and that was wrong. Verified
  2026-08-02: both trees are writable at the filesystem level, and the only
  guarantee is `~/.claude/hooks/block-rm.sh`, a **per-user Claude Code hook**
  that does nothing in a plain shell, nothing outside a CC session, and nothing
  for another user until they install their own. See
  [`docs/state_of_the_project.md`](docs/state_of_the_project.md) §8. Retire
  superseded frames in `data/ready_to_delete.md` rather than deleting them.
- **Read access is governed by an Isilon ACL the client cannot see**, and it
  does not match the POSIX mode. Frames under `append_only/inhibition/` show
  mode `rwxrwx---` with group `ssmd-ud-vmlab`; a member of that group can still
  get `Permission denied` on every one of them (measured for `@tt8804`,
  2026-08-04: 0 of 166 sampled files readable, and several experiment
  directories not even listable). If the GUI comes up with every panel absent,
  check `test -r` on a frame before debugging the app. `load_frame` carries the
  real reason (`no frame: [Errno 13] Permission denied: .../D1_30.parquet`), but
  **`seed_status` does not** — [`data.py:186`](integration/app/data.py#L186)
  catches every exception and reports the fixed string `"no frame written yet"`,
  so an unreadable seed is displayed as an unfinished one. That is a live
  instance of the catalogue's disguise #4 in
  [`how_this_project_breaks.md`](docs/how_this_project_breaks.md), and the
  docstring three lines above it warns against exactly this confusion.
- **Every docked pose is persisted, always** (@tt8804, 2026-08-11, #44). A
  screen run MUST write the whole pose cloud grouped by mode, not just each
  mode's representative. `nac_screen_v2` docks into a `tempfile.mkdtemp` and
  `rmtree`s it in a `finally`, so by default the 500 poses behind every mode are
  destroyed the moment the run ends — verified by a filesystem-wide search for
  `t4_716800c125a7`, which found no `.dlg` or `.pdbqt` anywhere. What survives
  is one medoid per mode and a per-pose table of measurements with no
  coordinates. That is why the ranking view can never show the cloud behind a
  mode, why mode membership cannot be re-derived, and why the pose that was
  simulated cannot be matched back to its energy rank. **Run with
  `--all-poses`.** It must come from the SAME run that produced the scores:
  docking is stochastic with no fixed seed, so a re-dock gives a different cloud
  and showing it beside the existing numbers puts structures on screen that the
  scores were not computed from.
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
