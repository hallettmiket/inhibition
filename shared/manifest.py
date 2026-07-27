"""
Purpose: Run manifests — the provenance spine every stage writes through.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: stage name, config used, input/output paths
Output: manifest.json alongside the stage's outputs

WHY EVERY STAGE WRITES ONE. The choreography is not a one-shot script. It is
meant to be run repeatedly and parameterized differently — another target,
another seed, another checkpoint, another warhead library. That only works if
each run records what it consumed, so that two results can be compared with
confidence about *what actually differed between them*.

A manifest answers, without re-running anything:
  - which git commit of this repo produced it
  - which config files, at which content hash
  - which inputs, at which content hash
  - which tool versions
  - when, on which host, under which run id

That last point matters more than it looks. "The numbers changed" is only
diagnosable if you can tell whether the code changed, the config changed, or
the input changed. Without a manifest all three look identical.

Manifests are written into the append-only tree next to their outputs, never
into git — they reference large data and are themselves run artifacts.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import platform
import subprocess
import uuid
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

log = logging.getLogger(__name__)

_REPO_ROOT = Path(__file__).resolve().parent.parent


def sha256_file(path: Path | str) -> str | None:
    """SHA-256 of a file, streamed. None if the path is not a readable file.

    Streamed rather than read-whole because MD trajectories and generated
    libraries reach the GB range and a manifest must never be the thing that
    exhausts memory.
    """
    p = Path(path)
    if not p.is_file():
        return None
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_tree(path: Path | str, *, max_files: int = 5000) -> str | None:
    """Stable hash over a directory's contents (names + file hashes).

    Used for directory-shaped inputs such as a decoy set or a fragment DB.
    Returns None if the directory is missing or exceeds ``max_files``, rather
    than spending unbounded time hashing.
    """
    p = Path(path)
    if not p.is_dir():
        return None
    files = sorted(f for f in p.rglob("*") if f.is_file())
    if len(files) > max_files:
        log.warning("sha256_tree: %s has %d files (> %d); not hashed",
                    p, len(files), max_files)
        return None
    h = hashlib.sha256()
    for f in files:
        h.update(str(f.relative_to(p)).encode("utf-8"))
        fh = sha256_file(f)
        if fh:
            h.update(fh.encode("ascii"))
    return h.hexdigest()


def git_commit(repo: Path | None = None) -> dict[str, str | bool]:
    """Current git commit of the code, and whether the tree is dirty.

    A dirty tree is recorded rather than rejected — blocking a run because a
    file is unsaved would be worse than noting it — but ``dirty: true`` means
    the commit hash does NOT fully describe the code that ran, and any result
    carrying it should be treated as provisional.
    """
    root = repo or _REPO_ROOT
    out: dict[str, str | bool] = {"commit": "unknown", "dirty": False}
    try:
        out["commit"] = subprocess.check_output(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            stderr=subprocess.DEVNULL, text=True).strip()
        status = subprocess.check_output(
            ["git", "-C", str(root), "status", "--porcelain"],
            stderr=subprocess.DEVNULL, text=True)
        out["dirty"] = bool(status.strip())
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        log.warning("git provenance unavailable: %s", exc)
    if out["dirty"]:
        log.warning(
            "working tree is DIRTY — commit %s does not fully describe the "
            "code producing this run; treat its outputs as provisional.",
            str(out["commit"])[:12])
    return out


def tool_versions(names: list[str] | None = None) -> dict[str, str]:
    """Best-effort versions of the python tools in the active environment."""
    names = names or ["rdkit", "pandas", "numpy", "crem", "torch", "openmm"]
    versions: dict[str, str] = {"python": platform.python_version()}
    for n in names:
        try:
            mod = __import__(n)
            versions[n] = getattr(mod, "__version__", "unknown")
        except Exception:  # noqa: BLE001 - absent tool is data, not an error
            versions[n] = "not installed"
    return versions


@dataclass
class Manifest:
    """A single stage's provenance record.

    Attributes
    ----------
    stage : str
        Stage name, e.g. ``"t4_enumeration"``.
    approach : str
        Approach id (``t1``..``t4``, or ``shared``).
    params : dict
        The parameters this run was given — the thing that differs between
        two runs of the same code, and therefore the first thing to record.
    """

    stage: str
    approach: str = "shared"
    run_id: str = field(default_factory=lambda: uuid.uuid4().hex[:12])
    params: dict = field(default_factory=dict)
    inputs: dict = field(default_factory=dict)
    outputs: dict = field(default_factory=dict)
    configs: dict = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)

    def add_input(self, name: str, path: Path | str) -> "Manifest":
        """Record an input file or directory with its content hash."""
        p = Path(path)
        self.inputs[name] = {
            "path": str(p),
            "sha256": sha256_tree(p) if p.is_dir() else sha256_file(p),
            "exists": p.exists(),
        }
        return self

    def add_output(self, name: str, path: Path | str) -> "Manifest":
        """Record an output file or directory with its content hash."""
        p = Path(path)
        self.outputs[name] = {
            "path": str(p),
            "sha256": sha256_tree(p) if p.is_dir() else sha256_file(p),
            "bytes": p.stat().st_size if p.is_file() else None,
        }
        return self

    def add_config(self, path: Path | str) -> "Manifest":
        """Record a config file by content hash.

        Hashing the config rather than copying it means a manifest stays small
        while still detecting that a run used different thresholds.
        """
        p = Path(path)
        self.configs[p.name] = {"path": str(p), "sha256": sha256_file(p)}
        return self

    def note(self, text: str) -> "Manifest":
        """Attach a human-readable note (a caveat, a skipped step, a warning)."""
        self.notes.append(text)
        return self

    def to_dict(self) -> dict:
        return {
            "run_id": self.run_id,
            "stage": self.stage,
            "approach": self.approach,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "host": platform.node(),
            "user": os.environ.get("USER", "unknown"),
            "git": git_commit(),
            "tools": tool_versions(),
            "params": self.params,
            "configs": self.configs,
            "inputs": self.inputs,
            "outputs": self.outputs,
            "notes": self.notes,
        }

    def write(self, out_dir: Path | str, *, filename: str = "manifest.json") -> Path:
        """Write the manifest beside a stage's outputs.

        Existing manifests are never overwritten; a second run in the same
        directory writes ``manifest.<run_id>.json``. Overwriting would destroy
        the record of the run whose outputs may still be sitting there.
        """
        d = Path(out_dir)
        d.mkdir(parents=True, exist_ok=True)
        target = d / filename
        if target.exists():
            target = d / f"manifest.{self.run_id}.json"
        target.write_text(json.dumps(self.to_dict(), indent=2) + "\n", encoding="utf-8")
        log.info("manifest [%s] %s -> %s", self.run_id, self.stage, target)
        return target
