"""
Purpose: Resolve agent analysis outputs under the governed append-only root.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-01
Input: an agent name, a topic, a file stem and suffix
Output: a versioned path under /data/lab_vm/append_only/inhibition/00_outputs/

WHY THIS EXISTS. Analysis artefacts were being written to `outputs/` INSIDE THE
REPO -- 1.1 MB of gzipped search trees, rendered pose HTML and benchmark CSVs
sitting in git. `rules/data-storage.md` and this repo's own README both say
derived data lives under the governed root and only very small files stay in
the repo, and the repo was already inconsistent about it:
`build_noncovalent_decoys_2.py` wrote to `append_only/` while the whole redock
benchmark wrote to `outputs/`. There is no version of that split that is
correct, so it is resolved toward the rule.

WRITES ARE VERSIONED BECAUSE THE ROOT IS APPEND-ONLY. Every filename under
`outputs/` was pinned `_1` and rewritten in place on each run -- so re-running
`redock_01_survey_pdb.py` silently replaced the survey that
`redock_02_build_cases.py` had already consumed, with nothing recording that
the input to a finished benchmark had changed underneath it. Under the governed
root the protected-paths hook refuses the overwrite outright, and
`write_path()` returns the next integer instead, so a re-run produces `_2`
beside `_1` rather than destroying it.

READS RESOLVE TO THE NEWEST. The counterpart: `latest_path()` globs rather than
naming a version, for exactly the reason `reference_set.latest_reference()`
does. A pinned `_1` read cannot announce that `_2` exists -- that is catalogue
entries #6 and #7 in `docs/how_this_project_breaks.md`, and five of the six
scripts here read a version literal written by another script.

LAYOUT. `00_outputs/<agent>/<topic>/<stem>_<N>.<suffix>`, mirroring the
`00_shared_substrate/` convention already used for the gate token and decoys:
a `00_` prefix sorts it above the numbered experiments and marks it as
belonging to no single one.
"""

from __future__ import annotations

from pathlib import Path

from . import io as dio

ROOT = Path("/data/lab_vm/append_only/inhibition/00_outputs")


def dir_for(agent: str, topic: str | None = None) -> Path:
    """The directory for one agent's outputs on one topic, created if absent."""
    d = ROOT / agent if topic is None else ROOT / agent / topic
    d.mkdir(parents=True, exist_ok=True)
    return d


def write_path(agent: str, topic: str | None, stem: str, suffix: str) -> Path:
    """The next unused version to write. Never returns an existing path."""
    return dio.next_version(dir_for(agent, topic), stem, suffix)


def latest_path(agent: str, topic: str | None, stem: str,
                suffix: str) -> Path:
    """The newest version of a stem.

    Raises rather than returning None: every caller here is reading an input it
    genuinely requires, and a None flowing into `pd.read_csv` produces a
    confusing failure some distance from the actual cause.
    """
    p = dio.latest(dir_for(agent, topic), stem, suffix)
    if p is None:
        raise FileNotFoundError(
            f"no {stem}_<N>{suffix} under {dir_for(agent, topic)} — "
            "has the stage that writes it been run?")
    return p


class Topic:
    """One agent's outputs on one topic, with the version policy attached.

    RESOLVE A WRITE PATH ONCE PER RUN, THEN REUSE IT. `redock_04_rmsd.py`
    writes its summary JSON twice -- a checkpoint, then the same payload with a
    cross-check folded in. Calling `write()` twice would allocate `_1` and `_2`
    and leave a half-populated file looking like a separate result. Bind the
    path to a name and write to that name:

        out = Topic("blacksmith", "redock_pin1")
        dest = out.write("redock_summary", ".json")
        dest.write_text(...)      # checkpoint
        dest.write_text(...)      # same run, same artefact

    Re-running the SCRIPT still produces a new version, which is the property
    the append-only root is for: one run, one artefact, nothing overwritten.
    """

    def __init__(self, agent: str, topic: str | None = None) -> None:
        self.agent, self.topic = agent, topic

    @property
    def dir(self) -> Path:
        return dir_for(self.agent, self.topic)

    def write(self, stem: str, suffix: str = ".csv") -> Path:
        return write_path(self.agent, self.topic, stem, suffix)

    def latest(self, stem: str, suffix: str = ".csv") -> Path:
        return latest_path(self.agent, self.topic, stem, suffix)
