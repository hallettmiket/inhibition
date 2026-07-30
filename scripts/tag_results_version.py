"""
Purpose: Freeze the current complement of results as a named, checkable version.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-30
Input: the append-only data root, the git working tree, the gate token
Output: results_versions/<label>.json + an annotated git tag

Run:  python scripts/tag_results_version.py --label 2026.07.30.01 [--dry-run]
      python scripts/tag_results_version.py --list

WHY THIS IS CHEAP HERE. The data root is append-only: no derived file is ever
overwritten, and every approach's frame is integer-versioned. So a "version of
the results" does not need to copy anything. It is a MANIFEST that names the
exact file version each approach was at, with a SHA-256 for each, plus the git
commit of the code that produced them. Re-reading a version means reading the
files it names; nothing has to be preserved specially because nothing was ever
at risk of being overwritten.

WHAT A VERSION IS FOR. Two things, and they are different.

  1. "What did we believe on 30 July?" -- answerable later, exactly, including
     the decision records that were live and the gate verdicts that stood.
  2. "What changed between then and now?" -- a diff over two manifests, which
     is what makes a claim like "fixing X moved the verdict" checkable rather
     than remembered.

REFUSES A DIRTY TREE BY DEFAULT. A version whose git commit does not describe
the code that produced the files is not a version, it is a timestamp. The
manifest system already warns about this per-run (`working tree is DIRTY`);
here it is fatal, because a release is exactly the moment the provenance has
to be true.

INCOMPLETE COMPUTATIONS ARE RECORDED, NOT HIDDEN. A version tagged while
GROMACS is still running is legitimate -- but the manifest says so, with the
completed count against the expected one, so nobody later mistakes a partial
campaign for a finished one.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

log = logging.getLogger("tag-results")

DATA = Path("/data/lab_vm/append_only/inhibition")
VERSIONS = DATA / "00_shared_substrate" / "results_versions"
GATE_TOKEN = DATA / "00_shared_substrate" / "enrichment_gate.token"

EXPERIMENTS = {
    "t1": ("01_t1_de_novo", "D1"),
    "t2": ("02_t2_atra_crem", "D2"),
    "t3": ("03_t3_reinvent", "D3"),
    "t4": ("04_t4_combinatorial", "D4"),
}


def _sha256(p: Path, cap: int = 512 * 1024 * 1024) -> str:
    h = hashlib.sha256()
    n = 0
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
            n += len(chunk)
            if n > cap:
                return f"sha256-partial-{cap}:{h.hexdigest()}"
    return h.hexdigest()


def _git(*args: str) -> str:
    return subprocess.run(["git", "-C", str(REPO), *args],
                          capture_output=True, text=True).stdout.strip()


def latest_frame(experiment: str, prefix: str) -> Path | None:
    d = DATA / experiment
    if not d.is_dir():
        return None
    frames = sorted(d.glob(f"{prefix}_*.parquet"),
                    key=lambda p: int(p.stem.split("_")[-1]))
    return frames[-1] if frames else None


def gromacs_state() -> dict:
    """Completed 10 ns replicates against the campaign's intent."""
    done = 0
    for p in DATA.glob("0[12]_*/gromacs/*/rep*/gromacs_result.json"):
        try:
            if float(json.loads(p.read_text()).get("production_ps", 0)) >= 10000:
                done += 1
        except Exception:  # noqa: BLE001
            continue
    expected = 240
    return {"completed_10ns_replicates": done, "expected": expected,
            "complete": done >= expected,
            "note": ("48 candidates x 5 replicates x 10 ns explicit TIP3P. "
                     "A partial count is recorded rather than rounded up.")}


def decision_records() -> list[dict]:
    out = []
    for p in sorted((REPO / "decisions").glob("D*.md")):
        status = None
        for line in p.read_text(errors="ignore").splitlines()[:30]:
            if line.startswith("status:"):
                status = line.split(":", 1)[1].strip()
                break
        out.append({"id": p.stem.split("-")[0], "file": p.name, "status": status})
    return out


def gate_verdicts() -> dict:
    if not GATE_TOKEN.is_file():
        return {"error": "no gate token"}
    try:
        tok = json.loads(GATE_TOKEN.read_text())
    except Exception as exc:  # noqa: BLE001
        return {"error": str(exc)[:200]}
    out = {}
    for stratum, v in tok.items():
        if not isinstance(v, dict):
            continue
        for metric, m in v.items():
            if isinstance(m, dict) and "verdict" in m:
                out[f"{stratum}/{metric}"] = {
                    k: m.get(k) for k in
                    ("verdict", "roc_auc", "n_actives", "n_decoys",
                     "n_chemotypes", "ef_1pct")}
    return out


def build(label: str) -> dict:
    frames = {}
    for a, (exp, prefix) in EXPERIMENTS.items():
        f = latest_frame(exp, prefix)
        frames[a] = ({"frame": str(f), "sha256": _sha256(f),
                      "bytes": f.stat().st_size} if f else
                     {"frame": None, "note": "no frame found"})
    return {
        "label": label,
        "git_commit": _git("rev-parse", "HEAD"),
        "git_describe": _git("describe", "--always", "--dirty"),
        "clean_tree": _git("status", "--porcelain") == "",
        "frames": frames,
        "gate": gate_verdicts(),
        "gromacs": gromacs_state(),
        "decisions": decision_records(),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--label", help="e.g. 2026.07.30.01")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--allow-dirty", action="store_true",
                    help="tag anyway; the manifest records that it is untrue")
    ap.add_argument("--allow-incomplete", action="store_true",
                    help="tag with computations still running (recorded)")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    if args.list:
        VERSIONS.mkdir(parents=True, exist_ok=True)
        for p in sorted(VERSIONS.glob("*.json")):
            d = json.loads(p.read_text())
            g = d.get("gromacs", {})
            print(f"{d['label']:<18} {d.get('git_commit','')[:8]}  "
                  f"gromacs {g.get('completed_10ns_replicates','?')}/"
                  f"{g.get('expected','?')}  "
                  f"{len(d.get('decisions', []))} decisions")
        return

    if not args.label:
        raise SystemExit("--label is required (e.g. 2026.07.30.01)")

    man = build(args.label)

    if not man["clean_tree"] and not args.allow_dirty:
        raise SystemExit(
            "working tree is dirty. A version whose commit does not describe "
            "the code that produced the files is a timestamp, not a version. "
            "Commit first, or pass --allow-dirty to record it as untrue.")
    if not man["gromacs"]["complete"] and not args.allow_incomplete:
        g = man["gromacs"]
        raise SystemExit(
            f"GROMACS is at {g['completed_10ns_replicates']}/{g['expected']} "
            "replicates. Wait, or pass --allow-incomplete — the manifest will "
            "record the partial count either way.")

    dest = VERSIONS / f"{args.label}.json"
    if dest.exists():
        raise SystemExit(f"{dest} already exists; versions are never rewritten")

    print(json.dumps({k: v for k, v in man.items() if k != "decisions"},
                     indent=2)[:2000])
    print(f"\n... plus {len(man['decisions'])} decision records")

    if args.dry_run:
        print("\n--dry-run: nothing written, no tag created")
        return

    VERSIONS.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(man, indent=2), encoding="utf-8")
    print(f"\nwrote {dest}")

    tag = f"results/{args.label}"
    r = subprocess.run(
        ["git", "-C", str(REPO), "tag", "-a", tag, "-m",
         f"results {args.label}\n\n"
         f"gromacs {man['gromacs']['completed_10ns_replicates']}/"
         f"{man['gromacs']['expected']} replicates\n"
         f"manifest: {dest}"],
        capture_output=True, text=True)
    print(f"git tag {tag}: {'ok' if r.returncode == 0 else r.stderr.strip()}")
    if r.returncode == 0:
        print(f"push it with:  git push origin {tag}")


if __name__ == "__main__":
    main()
