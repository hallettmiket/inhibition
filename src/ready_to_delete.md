# ready_to_delete

Append-only outputs and obsolete code considered safe to garbage-collect.
Required by `rules/data-storage.md`: nothing under `append_only/` is ever
deleted or overwritten in place, so superseded versions are **listed here** and
swept by the reconciliation routine rather than removed ad hoc.

Each entry says what the file is, what supersedes it, and why it is safe — a
bare list of paths is not reviewable.

---

## 2026-08-07 — `blacksmith/pose_split/pose_split_validation_{1..14}.csv`

**Supersedes:** `pose_split_validation_15.csv` (the complete 15-row result).

**Why they exist:** `scripts/pose_split_validation.py` checkpointed after every
molecule through `OUT.write()`, which mints a *new numbered file* on each call.
Fourteen partial prefixes of one table were the result — rows 1..1, 1..2, …
1..14 — and "the result" became ambiguous, since a reader has to know that only
the highest number is complete.

**Fixed at source:** the script now checkpoints to a fixed
`_partial_pose_split_validation.csv` and writes the versioned artefact once, at
the end, removing the partial on success. This class of clutter cannot recur.

**Safe because:** every row in files 1–14 is byte-identical to the
correspondingly-numbered prefix of file 15. Nothing is lost. Confirmed by
`pd.read_csv(f_n).equals(pd.read_csv(f_15).head(n))`.
