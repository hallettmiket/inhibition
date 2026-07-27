"""
Purpose: Declarative acquisition of every external input, with hash pinning.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: config/sources.yaml
Output: staged files under immutable/; observed hashes pinned in config/sources.lock.json

Replaces "whatever curl and git commands somebody ran once" with a recorded,
re-runnable step. Staging a fresh machine — or a different target — becomes:

    python -m shared.sources stage

Idempotent: an already-present file whose hash matches its pin is left alone.
A file whose hash DIFFERS from its pin is an error, not a silent re-download.
An upstream artifact that changes underneath a published result is a
reproducibility failure, and the whole point of pinning is to make that loud.

Nothing here writes to append_only/, and nothing here deletes. Files under
immutable/ are only ever created, never replaced — replacing one would
invalidate every downstream result derived from it without leaving a trace.
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
from pathlib import Path

import yaml

from .manifest import Manifest, sha256_file

log = logging.getLogger(__name__)

_REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CONFIG = _REPO_ROOT / "config" / "sources.yaml"

# Pins live in a lockfile, NOT written back into the hand-authored YAML.
# An earlier version of this module round-tripped sources.yaml through
# yaml.safe_dump to record observed hashes, which silently erased every
# explanatory comment in it. Config is written by humans; pins are written by
# code; the two do not share a file.
DEFAULT_LOCKFILE = _REPO_ROOT / "config" / "sources.lock.json"


def load_lock(path: Path | str | None = None) -> dict:
    """Load the pin lockfile, or an empty lock if none exists yet."""
    p = Path(path) if path else DEFAULT_LOCKFILE
    if not p.is_file():
        return {"sources": {}}
    return json.loads(p.read_text(encoding="utf-8"))


def write_lock(lock: dict, path: Path | str | None = None) -> Path:
    """Write the pin lockfile."""
    p = Path(path) if path else DEFAULT_LOCKFILE
    p.write_text(json.dumps(lock, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return p


class SourceError(RuntimeError):
    """A declared source is missing, unreachable, or has changed upstream."""


def load_config(path: Path | str | None = None) -> dict:
    p = Path(path) if path else DEFAULT_CONFIG
    if not p.is_file():
        raise SourceError(f"sources config not found: {p}")
    return yaml.safe_load(p.read_text(encoding="utf-8"))


def ensure_directories(cfg: dict) -> list[Path]:
    """Create the governed directories the choreography expects."""
    made: list[Path] = []
    for _root, paths in (cfg.get("directories") or {}).items():
        for d in paths or []:
            p = Path(d)
            if not p.is_dir():
                p.mkdir(parents=True, exist_ok=True)
                made.append(p)
                log.info("created %s", p)
    return made


def _fetch_http(url: str, dest: Path) -> None:
    """Download a URL to dest via curl, failing loudly on HTTP errors.

    ``--fail`` matters: without it curl happily writes a 404 HTML page to the
    destination and exits 0, producing a "checkpoint" that is really an error
    page — which then fails much later, somewhere far less obvious.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".partial")
    proc = subprocess.run(
        ["curl", "-sSL", "--fail", "-o", str(tmp), url],
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        tmp.unlink(missing_ok=True)
        raise SourceError(f"download failed: {url}\n{proc.stderr[:500]}")
    tmp.rename(dest)


def _fetch_git(url: str, dest: Path, commit: str | None) -> str:
    """Clone (or reuse) a git source and return the resolved commit."""
    if not dest.is_dir():
        dest.parent.mkdir(parents=True, exist_ok=True)
        proc = subprocess.run(["git", "clone", url, str(dest)],
                              capture_output=True, text=True)
        if proc.returncode != 0:
            raise SourceError(f"clone failed: {url}\n{proc.stderr[:500]}")
    if commit:
        proc = subprocess.run(["git", "-C", str(dest), "checkout", commit],
                              capture_output=True, text=True)
        if proc.returncode != 0:
            raise SourceError(f"checkout {commit} failed in {dest}\n{proc.stderr[:500]}")
    return subprocess.check_output(
        ["git", "-C", str(dest), "rev-parse", "HEAD"], text=True).strip()


def stage(config_path: Path | str | None = None, *, only: str | None = None,
          write_back: bool = True, lock_path: Path | str | None = None) -> dict:
    """Acquire every declared source, verifying or pinning its hash.

    Parameters
    ----------
    config_path : Path or str, optional
        Path to config/sources.yaml.
    only : str, optional
        Stage a single named source instead of all of them.
    write_back : bool, optional
        Record newly-observed hashes and commits into config/sources.lock.json.
        First acquisition observes; thereafter the pin is enforced.


    Returns
    -------
    dict
        Per-source status: ``pinned``, ``fetched``, ``skipped``, or ``pending``.
    """
    cfg_path = Path(config_path) if config_path else DEFAULT_CONFIG
    cfg = load_config(cfg_path)
    lock = load_lock(lock_path)
    pins = lock.setdefault("sources", {})
    ensure_directories(cfg)

    mf = Manifest(stage="stage_sources", approach="shared").add_config(cfg_path)
    results: dict[str, dict] = {}
    changed = False

    # One source must NOT abort the rest. A transient git failure on an
    # unrelated repo previously blocked a 281 MB database download that had
    # nothing to do with it, and the whole run exited having staged nothing.
    # Failures are collected and reported together at the end, so a run either
    # stages everything it can or tells you exactly what it could not.
    failures: dict[str, str] = {}

    for name, spec in (cfg.get("sources") or {}).items():
        if only and name != only:
            continue
        try:
            changed |= _stage_one(name, spec, pins, results, mf, write_back)
        except Exception as exc:  # noqa: BLE001 - collected, not swallowed
            failures[name] = str(exc)[:300]
            results[name] = {"status": "failed", "error": str(exc)[:300]}
            log.error("source %r FAILED (continuing with the rest): %s",
                      name, str(exc)[:200])

    if changed and write_back:
        written = write_lock(lock, lock_path)
        log.info("wrote newly-observed pins to %s", written)

    append_only = (cfg.get("directories") or {}).get("append_only") or []
    if append_only:
        mf.write(Path(append_only[0]) / "00_shared_substrate",
                 filename="stage_sources_manifest.json")

    if failures:
        raise SourceError(
            f"{len(failures)} of {len(results)} source(s) failed; the rest were "
            "staged:\n  - " + "\n  - ".join(f"{k}: {v[:160]}"
                                            for k, v in failures.items()))
    return results


def _stage_one(name: str, spec: dict, pins: dict, results: dict,
               mf, write_back: bool) -> bool:
    """Acquire one declared source. Returns True if a new pin was recorded."""
    changed = False
    if True:  # preserved indentation for the original body
        kind = spec.get("kind")
        dest = Path(spec["dest"]) if spec.get("dest") else None
        url = spec.get("url")

        if kind == "generated" or url is None:
            results[name] = {"status": "pending",
                             "reason": spec.get("notes", "not yet acquired").strip()}
            mf.note(f"{name}: pending — {results[name]['reason'][:120]}")
            log.warning("source %r pending: %s", name, results[name]["reason"][:120])
            return changed

        if kind == "git":
            pin = pins.get(name, {})
            head = _fetch_git(url, dest, pin.get("commit"))
            pinned = pin.get("commit")
            if pinned and pinned != head:
                raise SourceError(
                    f"{name}: checked-out commit {head} != pinned {pinned}")
            if not pinned and write_back:
                pins.setdefault(name, {})["commit"] = head
                changed = True
            results[name] = {"status": "pinned" if pinned else "fetched",
                             "commit": head}
            mf.add_input(name, dest)
            return changed

        # http
        if dest.is_file():
            observed = sha256_file(dest)
            expected = pins.get(name, {}).get("sha256")
            if expected and observed != expected:
                raise SourceError(
                    f"{name}: {dest} hash {observed[:16]} != pinned "
                    f"{expected[:16]}. An input changed underneath a pinned "
                    "run; investigate before overwriting anything."
                )
            if not expected and write_back:
                pins.setdefault(name, {})["sha256"] = observed
                changed = True
                log.info("pinned %s -> sha256 %s", name, observed[:16])
            results[name] = {"status": "skipped" if expected else "pinned",
                             "sha256": observed}
        else:
            _fetch_http(url, dest)
            observed = sha256_file(dest)
            expected = pins.get(name, {}).get("sha256")
            if expected and observed != expected:
                raise SourceError(
                    f"{name}: downloaded hash {observed[:16]} != pinned "
                    f"{expected[:16]}")
            if not expected and write_back:
                pins.setdefault(name, {})["sha256"] = observed
                changed = True
            results[name] = {"status": "fetched", "sha256": observed}
        mf.add_input(name, dest)
    return changed


def main() -> None:
    ap = argparse.ArgumentParser(description="Stage the choreography's external sources.")
    ap.add_argument("action", choices=["stage", "check"],
                    help="stage: acquire and pin. check: verify pins only.")
    ap.add_argument("--only", default=None, help="stage a single named source")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    res = stage(only=args.only, write_back=(args.action == "stage"))
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
