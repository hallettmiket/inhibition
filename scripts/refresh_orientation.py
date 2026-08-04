"""
Purpose: regenerate the measured numbers in the two orientation docs, and prove they are current.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-04
Input: the latest frame per experiment, decisions/, the gate token
Output: rewritten AUTO blocks in docs/state_of_the_project.md + a drift check

ISSUE #11. `state_of_the_project.md` says of itself: *"it drifted badly within
24 h of being written and a new maintainer read it as fact."* That is the whole
problem. These two documents are the project's tier-II memory -- the thing a
fresh Claude Code session reads to know what is true -- so a stale number in
them is not a documentation nit, it is context poisoning.

WHY A REFRESH SCRIPT ALONE WOULD NOT HAVE FIXED IT. A script only helps if
somebody runs it, and nobody runs it precisely when things are moving fastest,
which is when the drift happens. So the load-bearing half of this is not here:
it is `tests/test_orientation_current.py`, which FAILS THE SUITE when a number
in the doc no longer matches the data. The script makes the fix easy; the test
makes skipping it impossible.

WHAT IS AUTOMATED AND WHAT IS DELIBERATELY NOT. Only the counts that drift are
generated -- frame names, row counts, docked/ranked/shortlisted, per-seed pool
sizes, the decision count. **The prose is not touched.** "What is established",
"what is ruled out" and "what to do next" are judgements, and a generator that
rewrote them would produce a confident document nobody decided. Those sections
stay hand-written and stay the maintainer's responsibility.

The generated region is fenced by AUTO markers. A missing marker is a hard
error rather than a silent skip -- a refresh that quietly updates nothing is
the same failure as never running it.
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import io as dio                       # noqa: E402

log = logging.getLogger("refresh-orientation")

STATE_DOC = REPO / "docs" / "state_of_the_project.md"
DECISIONS = REPO / "decisions"

BEGIN = "<!-- AUTO:{key}:BEGIN -->"
END = "<!-- AUTO:{key}:END -->"

# (approach, frame prefix, experiment dir, label). The T_2 seeds each have
# their own experiment directory; they are read from config rather than listed
# here so a new seed cannot be silently omitted -- see `_t2_experiments`.
ARMS = [("t1", "D1", "01_t1_de_novo", "T_1 de novo (DiffSBDD)"),
        ("t2", "D2", None, "T_2 neighbourhood (CReM)"),
        ("t3", "D3", "03_t3_reinvent", "T_3 R-group (LibInvent)"),
        ("t4", "D4", "04_t4_combinatorial", "T_4 warhead x R-group")]


def _t2_experiments() -> list[tuple[str, str]]:
    """(label, experiment) for every T_2 seed, from the app's own resolver.

    Read from `integration/app/data.py` rather than hardcoded: `04_rank.py`
    naming ONE T_2 experiment when five existed is catalogue entry #16, and a
    hardcoded list here would reproduce it in the document that exists to stop
    people being misled.
    """
    sys.path.insert(0, str(REPO / "integration" / "app"))
    import data as appdata                          # noqa: PLC0415
    return [(v["label"], v["experiment"])
            for v in appdata.variants("t2").values()]


def _summarise(prefix: str, experiment: str) -> dict | None:
    root = dio.DATA_ROOT if hasattr(dio, "DATA_ROOT") else Path(
        "/data/lab_vm/append_only/inhibition")
    path = dio.latest(root / experiment, prefix, ".parquet")
    if path is None:
        return None
    df = dio.read_frame(path)
    metric = "vina_affinity" if prefix in ("D1", "D2") else "affinity_kcal"
    short = ("shortlist_synth" if "shortlist_synth" in df.columns
             else "shortlist" if "shortlist" in df.columns else None)
    return {
        "frame": path.name,
        "rows": len(df),
        "docked": int(df[metric].notna().sum()) if metric in df.columns else 0,
        "ranked": int(df["rank"].notna().sum()) if "rank" in df.columns else 0,
        "shortlisted": int(df[short].sum()) if short else 0,
        "shortlist_col": short or "-",
    }


def measure() -> dict:
    """Every number the document quotes. One place, one read."""
    arms, t2 = {}, {}
    for approach, prefix, experiment, label in ARMS:
        if approach == "t2":
            for seed_label, exp in _t2_experiments():
                s = _summarise(prefix, exp)
                if s:
                    t2[seed_label] = s
            continue
        s = _summarise(prefix, experiment)
        if s:
            arms[label] = s
    return {
        "arms": arms,
        "t2": t2,
        "n_decisions": len(sorted(DECISIONS.glob("D0*.md"))),
        "total_t2": sum(v["docked"] for v in t2.values()),
    }


def render(m: dict) -> dict[str, str]:
    """AUTO block key -> markdown. Pure formatting; no judgement."""
    lines = ["| arm | frame | rows | docked | ranked | shortlist |",
             "|---|---|---:|---:|---:|---:|"]
    for label, s in m["arms"].items():
        lines.append(f"| {label} | `{s['frame']}` | {s['rows']:,} | "
                     f"{s['docked']:,} | {s['ranked']:,} | "
                     f"{s['shortlisted']} (`{s['shortlist_col']}`) |")
    arms_tbl = "\n".join(lines)

    lines = ["| T_2 seed | frame | docked | ranked | shortlist |",
             "|---|---|---:|---:|---:|"]
    for label, s in m["t2"].items():
        lines.append(f"| {label} | `{s['frame']}` | {s['docked']:,} | "
                     f"{s['ranked']:,} | {s['shortlisted']} |")
    lines.append(f"| **all six** | | **{m['total_t2']:,}** | | |")
    t2_tbl = "\n".join(lines)

    return {"arms": arms_tbl, "t2": t2_tbl,
            "decisions": f"**{m['n_decisions']}** decision records."}


def _shown(p: Path) -> str:
    """Repo-relative when it is under the repo, absolute otherwise.

    `relative_to` RAISES on a path outside the repo, so a bare call here made
    the script crash while formatting its own error message -- surfacing a
    `ValueError` from pathlib instead of the "this doc is stale" it was trying
    to report. Found by the test that corrupts a copy in a tmp directory.
    """
    try:
        return str(p.relative_to(REPO))
    except ValueError:
        return str(p)


def _replace(text: str, key: str, body: str) -> str:
    b, e = BEGIN.format(key=key), END.format(key=key)
    pat = re.compile(re.escape(b) + r".*?" + re.escape(e), re.DOTALL)
    if not pat.search(text):
        raise SystemExit(
            f"marker {b} ... {e} not found in {STATE_DOC.name}. A refresh that "
            "silently updates nothing is the failure this script exists to "
            "prevent; add the markers or fix the key.")
    return pat.sub(f"{b}\n{body}\n{e}", text)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--check", action="store_true",
                    help="exit non-zero if the doc is stale; write nothing")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    blocks = render(measure())
    text = original = STATE_DOC.read_text(encoding="utf-8")
    for key, body in blocks.items():
        text = _replace(text, key, body)

    if args.check:
        if text != original:
            raise SystemExit(f"{_shown(STATE_DOC)} is STALE. Run: "
                             "python3 scripts/refresh_orientation.py")
        print(f"{STATE_DOC.name} is current.")
        return

    if text == original:
        print(f"{STATE_DOC.name} already current; nothing written.")
        return
    STATE_DOC.write_text(text, encoding="utf-8")
    print(f"refreshed {_shown(STATE_DOC)}")


if __name__ == "__main__":
    main()
