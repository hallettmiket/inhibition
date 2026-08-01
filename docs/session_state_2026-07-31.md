# Session state — 2026-07-31 12:55 EDT

Written before a planned client shutdown at 13:30. Everything here is on
**biodatsci**; nothing depends on the laptop that was connected to it.

## Resuming the Claude Code session

The transcript is persisted server-side and written continuously:

```
/home/UWO/mhallet/.claude/projects/-home-UWO-mhallet-repos-murmurent/6662da8c-133f-4ac0-b004-07f9c7feab6a.jsonl   (147 MB)
```

Two ways back in, in order of preference:

1. **Reconnect VSCode Remote-SSH to biodatsci.** The extension host has been up
   4 days; if `vscode-server` survives the disconnect the session is simply
   still there.
2. **If it does not**, from `/home/UWO/mhallet/repos/murmurent`:
   ```bash
   claude --resume 6662da8c-133f-4ac0-b004-07f9c7feab6a
   ```
   The project key is derived from the working directory, so start from
   `~/repos/murmurent` — not from `~/repos/inhibition` — or the session will
   not be found.

## Running compute — survives the shutdown

Verified detached: both docking chains have `PPID 1` and are their own session
leaders with no controlling terminal. The MD campaign's parent looks like a
Claude Code child, but `1598021` is *also* its own session leader
(`SESS=PGID=1598021`, TTY `?`), so it is reparented to init rather than
signalled when the extension host exits.

**UPDATED 13:12 — the neutral chains were stopped and everything relaunched at
pH 7.4.** @hallettmiket approved the recommendation below at 12:58.

| what | where | state at 13:12 | ETA |
|---|---|---|---|
| T_2 docking, liu_2024_c3 (16,806 mols) | GPU 4 | **relaunched 13:01, pH 7.4** | ~00:50 |
| T_2 docking, potter_astex (7,376) | GPU 4, queued | not started | ~06:00 |
| T_2 docking, du_xu (9,736) | GPU 7 | **relaunched 13:01, pH 7.4** | ~19:50 |
| T_2 docking, guo_pfizer (8,670) | GPU 7, queued | not started | ~02:00 |
| T_2 docking, **atra re-dock** (1,882) | GPU 7, queued | not started | ~03:20 |
| GROMACS explicit MD (T_1/T_2 synth shortlists) | GPUs 0,1,2,3,5,6 | 7 `gmx` left of 245 | soon |

ATRA is now in the chain too — its original D2_4 docking was neutral, and all
five seeds must be prepared identically for the cross-seed comparison to mean
anything.

Chain driver + logs:
```
/tmp/claude-162990128/-home-UWO-mhallet-repos-murmurent/6662da8c-133f-4ac0-b004-07f9c7feab6a/scratchpad/
    dock_reseeds.sh
    logs/chain_gpu4.log  logs/chain_gpu7.log  logs/dock_<seed>.log
```
**These live under `/tmp` and are not guaranteed across a host reboot.** They
are logs only — no result depends on them. Frames are written to
`/data/lab_vm/append_only/inhibition/02*_t2_*/` and only at the END of each
dock, so an interrupted run loses that seed's docking but nothing already
committed.

## The open decision that gates all of the above

`shared/noncovalent_dock_run.py` prepares every ligand with `Chem.AddHs()` on
neutral SMILES and `obabel <sdf> -O <pdbqt>` with **no `-p 7.4`**. Ligands are
docked exactly as drawn — neutral. Measured composition of the five T_2 pools:

| seed | ionizable group | fraction | net charge at pH 7.4 |
|---|---|---|---|
| guo_pfizer | phosphate | 83.4% | ≈ −2 |
| atra | carboxylic acid | 86.6% | ≈ −1 |
| du_xu | carboxylic acid | 85.5% | ≈ −1 |
| potter_astex | basic amine | 86.7% | ≈ +1 |
| liu_2024_c3 | — | — | ≈ 0 |

The cross-seed comparison is therefore across four charge states, all docked
neutral, into a pocket recognised by Lys63/Arg68/Arg69. Raised by Tim in issue
#5 for one seed; it generalises to the whole five-seed design. Plan posted as a
comment on #5.

### A cache trap that must be fixed with it

`_prepare_one` short-circuits on `if out.is_file()` and the ligand path is
keyed on `candidate_id` **alone** — nothing about protonation. Cached neutral
PDBQTs already on disk:

```
399M  02b_t2_liu_c3_crem/docking/ligands
231M  02d_t2_duxu_crem/docking/ligands
 45M  02_t2_atra_crem/docking/ligands
  0   02e_t2_guo_crem/docking/ligands   (not started)
```

Adding `-p 7.4` without invalidating these would silently re-dock the **same
neutral files** and report success. Same defect class as the class-pool cache
keyed on `class_id` alone (D0042 era): a cache keyed on less than its inputs.
The fix must write to a new ligand directory, not refresh the old one.

## DONE — chains stopped, fix landed, all five relaunched

Approved and executed 12:58–13:01.

1. Both neutral chains stopped (`kill -TERM` on the process groups); GPUs 4 and
   7 released.
2. `shared/noncovalent_dock_run.py`: `LIGAND_PH = 7.4`, `obabel ... -p 7.4`,
   and `LIGAND_PREP_TAG` in the ligand/pose directory names so the prep cache
   cannot serve neutral files under a pH-7.4 run. `ligand_ph` now appears in
   every manifest, so a frame states its own protonation instead of leaving a
   reader to infer it from the run date.
3. All five seeds relaunched, ATRA included.

**Verified live, not assumed.** `obabel -p 7.4` gives the right ionization at
the SMILES level (21b → −2, carboxylic acid → −1, primary amine → +1), and in
the running pipeline the polar-hydrogen (`HD`, H-bond donor) counts differ
between the old and new preparations of the same candidates: 3→2 and 2→1 where
acids deprotonate, 0→1 where an amine protonates.

### One correction to the memo's reasoning, worth carrying

The memo warns that a ranking over a phosphate-seeded library "will rank
substantially on formal charge" because Lys63/Arg68/Arg69 is a basic cluster.
The concern is right but **the mechanism is not electrostatic**: AutoDock Vina
has no electrostatic term, and obabel writes all-zero partial charges into the
PDBQT either way — verified here, with and without `-p`. Vina never reads the
charge.

What protonation actually changes is **atom typing**: deprotonating a phosphate
removes two polar hydrogens, and Vina's H-bond term counts donors from those
`HD` atoms. Docking the neutral acid presents H-bond donors that do not exist at
pH 7.4. So the fix matters, and charge stratification in ranking is still worth
doing — but nobody should expect it to remove an electrostatic bias, because
there was never an electrostatic term to produce one.

## Uncommitted work in the tree

`git status` at 12:55 — GUI work from issue #3 plus the stale-module guard:

```
M integration/app/app.py      (stale guard moved to top; PANEL_SCOPE via getattr)
M integration/app/curate.py   (PANEL_SCOPE, UNFILTERED_FACTS, sidebar scope)
M integration/app/pose3d.py   (surface + sub-pockets, all poses, PyMOL export)
? tests/test_subpockets.py  test_pose_modes.py  test_curation_scope.py
? tests/test_app_renders.py  test_stale_guard.py
? outputs/
```

292 tests pass under `dwi_cheminf`. Nothing committed. Last pushed commit is
`9c6a1d4` (branded footer).

## Open threads

- **#5** — T_2 extension plan posted; three decisions pending (stop-or-finish,
  phosphate protect-vs-label, receptor ensemble vs fixed).
- **#4** — measurement-problem plan. Asked whether to fold in the Pin1 PDB
  survey (190 entries, 120 drug-like ligands, 39 covalent spanning ≥5 warhead
  chemistries) and revise Phase 0.3. **Not yet answered.**
- **#3** — GUI work done and documented; not committed.
- **#2** — chemotype-counting definition, now on the critical path for #4.
- `affinity_kcal` is the CNN-best pose's affinity, not the best affinity —
  89% of covalent candidates affected, >50% shortlist churn. Asked whether to
  open an issue and dispatch the adversary. **Not yet answered.**
- Degree-2 T_2 sampler built and smoke-tested; queued for tonight.

---

# Addendum — 13:40, after the MD campaign completed

## Explicit-solvent MD: DONE, and it answers the question it was built for

`243 of 245` replicates succeeded; 285–764 ns/day; 49 distinct candidates, 48
with all five replicates. Merged into `D1_21.parquet` and `D2_21.parquet`.

**Explicit water stabilises the poses relative to implicit solvent.** Mean
change in ligand RMSD, explicit minus implicit:

- T_1: **−0.290 nm** (min −2.292, max +0.764) over 23 candidates
- T_2: **−0.469 nm** (min −1.340, max −0.006) over 24 candidates

That is the question `run_gromacs_explicit.py` was written to answer. Under
implicit solvent two T_1 shortlist candidates dissociated outright (RMSD 9.0 and
7.3 nm) and the docstring asked whether that was real or a GB artefact. On this
evidence it was substantially the solvent model: T_2 improved on **every**
candidate (max change is −0.006, i.e. nothing got worse).

Seven trajectories still exceed 3.3 nm ligand RMSD after PBC correction and are
flagged suspect — genuine dissociation or an analysis problem, not yet
distinguished.

Both failures are the same molecule, `t1_1224c0ee20c2` (rep2, rep3, `gmx`
SIGSEGV). Its rep1 is the run that produced ~1400 LINCS constraint-violation
dumps already listed in `ready_to_delete.md`. **One candidate is simply not
simulable under this protocol** — three independent failures by two different
mechanisms. Worth stating as a result rather than filing as noise.

## A defect in the merge, found in its own output — NOT YET FIXED

`D1_21` and `D2_21` carry `explicit_rmsd_replicate_sd_x` / `_y`,
`explicit_rmsd_replicate_min_x` / `_y`, and so on — **and no column under the
canonical name.**

Cause: `scripts/merge_gromacs_results.py` line ~193 drops stale columns from a
hardcoded list, `(*COLS, *IMPLICIT_COLS)`, which omits every aggregate column
the same function builds a few lines earlier — `explicit_rmsd_replicate_sd`,
`_min`, `_max`, `_ratio`, `n_replicates`, `explicit_rmsd_suspect_any`. On a
re-merge those already exist on the frame, survive the drop, and pandas suffixes
both copies. Anything downstream reading the plain name now gets nothing.

This is the same defect already documented for `affinity_kcal_x/_y` in
`covalent_dock_run.py`, reproduced because the drop list is maintained by hand
and the aggregate columns were added later. A list that must stay in sync with
code five lines away will not.

**The fix (not applied):** derive the drop list from the merge itself rather
than from a constant —

```python
incoming = [c for c in agg.columns if c != "candidate_id"]
df = df.drop(columns=[c for c in incoming if c in df.columns])
```

same for the implicit frame, plus a post-merge assertion that no column ends in
`_x`/`_y`, so a future miss fails loudly instead of shipping.

**Why it was not applied:** the client connection dropped mid-edit and the
Read/Edit tools began failing on their PreToolUse hook. Patching a file that
produces governed data without the normal read-and-verify path is not worth a
few hours' head start. `scripts/merge_gromacs_results.py` is **untouched** —
`git diff` on it is empty — so the tree is consistent. Re-running the merge
after the fix writes `D1_22`/`D2_22`; `D1_21`/`D2_21` should then be retired in
`ready_to_delete.md`.

The numbers reported above are unaffected: they come from
`explicit_ligand_rmsd_nm_mean`, which is in `COLS`, was dropped correctly, and
carries no suffix.

## Other state at 13:40

- Docking chains healthy on GPUs 4 and 7, ~40 min in, running `ligands_ph7.4`.
- **GPU 6 belongs to another user** (`wzhan564`, ~50 GB). Not ours; leave it.
  GPUs 0,1,2,3,5 are now idle since the MD finished.
- `D1_21`/`D2_21` manifests carry the DIRTY-tree warning — they were written
  against commit `9c6a1d4` plus the uncommitted GUI and protonation changes.
