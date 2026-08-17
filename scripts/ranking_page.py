#!/usr/bin/env python3
"""
Purpose: build the Ranking page (modes.html) from the ranking alone.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-17
Input: the current run.topic's rank_v2 table (via shared.mode_ranking.gather)
Output: <reports>/modes.html + mode_poses/ + mode_thumbs/

WHY THIS IS ITS OWN SCRIPT. `modes.html` was written by `mdprio_combine`, which
exists to combine 100 ns REPORTS -- and which exits "no reports found for any
requested candidate" before it gets to the ranking. So the Ranking page, whose
only input is stage 2, could not be built until stage 5 had produced something.

On a fresh run that is backwards in the most visible way: 2,019 ranked modes sat
on disk while the page showing them stayed an awaiting-stage placeholder, and
the first stage whose results a reader wants was the last one to appear.
@tt8804, with the ranking already finished: "so I should be seeing a ranked list
soon?"

Nothing here is new code -- `mode_ranking.build` and `mode_assets.write_assets`
were always independent of the reports. Only the call site was coupled.

`mdprio_combine` still writes the page too, so the two must not drift; both call
the same builder with the same arguments, which is why the 3Dmol payload and the
`no_pose` list are resolved here exactly as they are there.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import logging
import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import mode_assets as massets                # noqa: E402
from shared import mode_ranking as moderank              # noqa: E402
from shared import run_paths as rp                       # noqa: E402
from shared import target_config as tc                   # noqa: E402

log = logging.getLogger("ranking-page")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--title", default="DWI covalent screen")
    ap.add_argument("--force-assets", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    r = moderank.gather()
    if r.empty:
        # Not an error: before stage 2 there is genuinely nothing to rank, and
        # build_gui's placeholder already says so. Exiting non-zero here would
        # make a watch loop log a failure every minute of a normal screen.
        log.info("no ranked modes for topic %s yet — nothing to build", tc.topic())
        return

    out = rp.reports_dir()
    # A molecule whose representative predates this run must not have its pose
    # drawn from the older clustering: mode 2 of one run is not mode 2 of
    # another. `expected` lets write_assets detect that and report it as stale.
    need = {}
    if "parent_ident" in r.columns and "mode" in r.columns:
        for p, g in r.groupby("parent_ident"):
            need[str(p)] = int(g["mode"].max()) + 1
    a = massets.write_assets(
        out, moderank.idents(r),
        force=args.force_assets or os.environ.get("MODE_ASSETS_FORCE") == "1",
        expected=need)
    log.info("mode assets: +%d poses, +%d thumbs, %d molecules with no pose "
             "from this run", a["poses"], a["thumbs"], len(a.get("stale", [])))

    # The viewer library is inlined from the same cache mdprio_combine uses, so
    # the page works with no network and the two builders cannot ship different
    # versions of it.
    cache = REPO / "scripts" / ".cache_3dmol-min.js"
    three = cache.read_text() if cache.is_file() else ""
    if not three:
        log.warning("no cached 3Dmol payload at %s — the pose viewer will be "
                    "blank rather than absent", cache.name)

    (out / "modes.html").write_text(
        moderank.build(args.title, _dt.date.today().isoformat(), three,
                       no_pose=a.get("stale", [])))
    n_mol = r.parent_ident.nunique() if "parent_ident" in r.columns else 0
    print(f"\n  {len(r):,} modes over {n_mol:,} molecules ({tc.topic()}) "
          f"-> {out / 'modes.html'}")


if __name__ == "__main__":
    main()
