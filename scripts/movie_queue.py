#!/usr/bin/env python3
"""
Purpose: which long-run poses have a movie OF THAT RUN, and which still need one.
Author: Timothy Wu (with Claude Code)
Date: 2026-09-09
Input: every completed run of at least --min-ns, across every results topic
Output: a status table; --build renders the missing movies

@twu383, 2026-09-09: *"can we start a list of the poses we want to have 100 ns
md movies ready for"*.

WHY THIS IS A SCRIPT AND NOT A LIST. A hand-kept list of poses goes stale the
moment a run finishes, and this project has a catalogue entry for exactly that
(#15, #19: a value that was right when written and cannot announce that it is
not any more). The set of poses wanting a movie is derivable -- it is every
completed long run -- so it is derived.

THE TRAP THIS EXISTS TO CATCH. Every one of the first 16 long runs HAD a viewer
page and a movie that played, and not one of those movies was of the 100 ns
trajectory: they were built from the 1.2 ns or 10 ns sweep and never rebuilt
when the longer run finished. A page that exists is not a page that is current,
and the movie carries no visible clue -- the slider reads frame counts. So the
check here is not "is there a page" but "does the movie's own TOTAL_PS match the
run's length", read out of the rendered page.
"""

from __future__ import annotations

import argparse
import glob
import logging
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import run_paths as rp                        # noqa: E402

log = logging.getLogger("movie-queue")
PY = Path.home() / ".micromamba/envs/dwi_reactive/bin/python"

#: Topics whose report trees may hold a viewer for a pose. Searched in order.
REPORT_TOPICS = ("holders_100ns", "thiadiazole_100ns", "nac_v8")
#: A movie is "of this run" if its length matches within this many ps.
LENGTH_TOL_PS = 500.0


def completed_runs(min_ns: float) -> pd.DataFrame:
    """Every finished run of at least `min_ns`, from every results topic.

    Two shapes reach here and both count: an `attack_sweep` row (the adaptive
    runs, whose `sweep_ps` is what they actually reached) and an
    `md_residence` row (the direct elevations, `production_ps`). Reading only
    one would silently omit half the long runs.
    """
    B = rp.BLACKSMITH
    rows = []
    for d in sorted(B.glob("attack_sweep_*")):
        for f in d.glob("*.csv"):
            try:
                t = pd.read_csv(f)
            except Exception:                             # noqa: BLE001
                continue
            if "ident" not in t.columns or "sweep_ps" not in t.columns:
                continue
            t = t[t.get("status", "").astype(str) == "ok"]
            for r in t.itertuples():
                ps = pd.to_numeric(getattr(r, "sweep_ps", None), errors="coerce")
                if ps == ps and ps / 1000 >= min_ns:
                    rows.append(dict(ident=str(r.ident), ps=float(ps),
                                     source=d.name))
    for f in B.glob("md_residence_*/*.csv"):
        try:
            t = pd.read_csv(f)
        except Exception:                                 # noqa: BLE001
            continue
        if "ident" not in t.columns or "production_ps" not in t.columns:
            continue
        if "status" in t.columns:
            t = t[t.status.astype(str) == "ok"]
        for r in t.itertuples():
            ps = pd.to_numeric(getattr(r, "production_ps", None), errors="coerce")
            if ps == ps and ps / 1000 >= min_ns:
                rows.append(dict(ident=str(r.ident), ps=float(ps),
                                 source=f.parent.name))
    if not rows:
        return pd.DataFrame(columns=["ident", "ps", "source"])
    d = pd.DataFrame(rows)
    # the LONGEST run of a pose is the one a movie should show
    return d.sort_values("ps").drop_duplicates("ident", keep="last")


def movie_length_ps(ident: str) -> tuple[str | None, float | None]:
    """(where the viewer lives, what length its movie covers).

    Read out of the rendered page rather than inferred from its existence --
    `const TOTAL_PS` is what the viewer actually plays.
    """
    for topic in REPORT_TOPICS:
        R = rp.reports_dir(topic)
        for p in (R / "sweep_pages" / f"{ident}.html", R / f"{ident}.html"):
            if not p.is_file():
                continue
            head = p.read_text(errors="replace")
            m = re.search(r"const TOTAL_PS = ([0-9.]+)", head)
            return f"{topic}/{p.parent.name}", (float(m.group(1)) if m else None)
    return None, None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--min-ns", type=float, default=50.0,
                    help="a run this long or longer wants a movie")
    ap.add_argument("--build", action="store_true",
                    help="render the missing ones (minutes each)")
    ap.add_argument("--topic", default="holders_100ns",
                    help="report tree to build into")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    d = completed_runs(args.min_ns)
    if d.empty:
        print(f"  no completed runs of {args.min_ns:.0f} ns or more")
        return
    where, have = [], []
    for i in d.ident:
        w, ps = movie_length_ps(i)
        where.append(w or "—")
        have.append(ps)
    d = d.assign(viewer=where, movie_ps=have)
    d["current"] = [
        (mp is not None and abs(mp - ps) <= LENGTH_TOL_PS)
        for mp, ps in zip(d.movie_ps, d.ps)]
    d = d.sort_values(["current", "ps"], ascending=[True, False])

    print(f"\n  POSES WITH A RUN >= {args.min_ns:.0f} ns: {len(d)}\n")
    print(f"  {'pose':<24}{'run':>8}{'movie':>10}  {'status':<26}viewer")
    for r in d.itertuples():
        mv = f"{r.movie_ps/1000:.1f} ns" if r.movie_ps else "none"
        st = ("current" if r.current
              else "NO MOVIE" if r.movie_ps is None
              else f"STALE — shows {r.movie_ps/1000:.1f} ns")
        print(f"  {r.ident:<24}{r.ps/1000:7.1f}{mv:>10}  {st:<26}{r.viewer}")
    need = d[~d.current]
    print(f"\n  current: {int(d.current.sum())}   NEED A MOVIE: {len(need)}")

    if not args.build:
        if len(need):
            print(f"\n  build them with:  {Path(__file__).name} --build "
                  f"--topic {args.topic}")
        return

    for r in need.itertuples():
        log.info("building %s (%.0f ns)", r.ident, r.ps / 1000)
        for step in (["scripts/sweep_assets.py", "--topic", args.topic,
                      "--idents", r.ident],
                     ["scripts/sweep_report.py", "--topic", args.topic,
                      "--ident", r.ident]):
            out = subprocess.run([str(PY)] + step, cwd=str(REPO),
                                 capture_output=True, text=True)
            if out.returncode != 0:
                log.warning("%s failed: %s", step[0],
                            (out.stderr or "").strip()[-200:])
                break


if __name__ == "__main__":
    main()
