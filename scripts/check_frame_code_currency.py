"""
Purpose: does the code that WOULD produce a frame differ from the code that DID?
Author: Mike Hallett (with Claude Code)
Date: 2026-08-04
Input: each generation stage's latest frame manifest (git.commit) + git history
Output: a per-stage report; `--check` exits non-zero on unacknowledged staleness

THE DEFECT CLASS THIS EXISTS FOR. A fix lands in a generation stage, the stage
is never re-run, and the fix is inert on the data everybody is looking at.
Nothing announces it: the frame is present, populated and plausible, and the
code reads as though the fix is in force.

Found twice on 2026-08-04, in the same audit:

* **T_3's pocket size ceiling.** Added 2026-07-30 (`b776843`, "pocket-derived
  size ceiling for T_2 and T_3"). T_3's production frame `D3_2` was written
  2026-07-27 23:20 -- three days EARLIER -- and T_3 generation has never been
  re-run. `verify()` computes `heavy_atoms` and stamps `exceeds_pocket_ceiling`;
  no T_3 frame has ever carried that column. Measured impact: 2 of 5,396
  molecules exceed the 55-atom ceiling, both already rejected by `alerts` and on
  neither shortlist. Real, and nearly harmless.
* **T_2 ATRA, same commit.** ATRA was generated 2026-07-27; the four newer seeds
  on 2026-07-31, after the ceiling. Measured impact: ZERO -- no T_2 pool exceeds
  49 heavy atoms, so the ceiling would never have fired on any of them.

Both impacts are small. THE CLASS IS NOT. Neither was detectable by reading the
code or the frame, and the next instance may not be harmless.

WHY IT REPORTS RATHER THAN SIMPLY FAILING. Staleness is often the correct state:
re-running DiffSBDD or a 16,806-molecule CReM expansion to pick up a comment
change would be absurd. So known-stale stages are ACKNOWLEDGED here with a
reason and a MEASURED impact, and the test fails only on staleness nobody has
looked at. That is the allowlist discipline from D0051 -- name what passes, so
an unanticipated case is refused rather than admitted.

WHAT IT CANNOT SEE. A frame produced from a DIRTY working tree is not described
by its commit at all; `shared/manifest.py` already warns at write time, and this
reports it separately rather than pretending the comparison is meaningful.
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import io as dio                       # noqa: E402

log = logging.getLogger("frame-code-currency")

DATA_ROOT = Path("/data/lab_vm/append_only/inhibition")

# stage name -> (experiment dir, frame prefix, source files that define it).
# Scoped to GENERATION/ENUMERATION, the stages this audit covered. Downstream
# stages re-run cheaply and are re-run often; the generation stages are the ones
# expensive enough that a fix silently fails to reach the data.
STAGES = {
    "t1_generate": ("01_t1_de_novo", "D1",
                    ["approaches/t1_de_novo/01_generate.py"]),
    "t2_generate": ("02_t2_atra_crem", "D2",
                    ["approaches/t2_atra_crem/01_generate.py"]),
    "t3_generate": ("03_t3_reinvent", "D3",
                    ["approaches/t3_reinvent/01_generate.py"]),
    "t4_enumerate": ("04_t4_combinatorial", "D4",
                     ["approaches/t4_combinatorial/01_enumerate.py"]),
}

# Staleness somebody has looked at. Each entry MUST carry the measured impact,
# not merely an assertion that it is fine -- "we checked and it does not matter"
# is only worth something with the number attached.
#
# EVERY GENERATION FRAME WAS WRITTEN FROM A DIRTY TREE. Measured 2026-08-04:
# all four record `git.dirty = true`, so the commit they name does NOT describe
# the code that ran, and this comparison systematically OVER-reports. That is
# the right direction to be wrong in -- it asks a question rather than granting
# a pass -- but it means a STALE verdict is a prompt to measure, never a finding
# on its own. Two of the four flags below dissolved the moment they were
# measured; two were real. Read the impact line, not the tag.
ACKNOWLEDGED = {
    "t1_generate":
        "1c40ba4 (D0025 alerts attribution, T_1 floor 10 + size_class) is later "
        "than D1_3's recorded commit, but the tree was DIRTY and the change is "
        "demonstrably IN the data. MEASURED 2026-08-04: `size_class` is present "
        "on the frame and the floor is applied -- 1,376 rows stamped "
        "`degenerate_too_small`, 82 `too_large`. Zero impact.",
    "t4_enumerate":
        "72cf331 (dangling attachment points, colliding candidate ids, suffixed "
        "merge) is later than D4_6's recorded commit, tree DIRTY. MEASURED "
        "2026-08-04: 0 of 1,782 SMILES contain an unfilled attachment point, 0 "
        "fail to parse, and 0 candidate_id collisions across all nine pools "
        "(72,104 molecules, injective id->SMILES). All three fixes are "
        "reflected in the data. Zero impact.",
    "t3_generate":
        "Pocket ceiling (b776843, 2026-07-30) postdates D3_2 (2026-07-27). "
        "MEASURED 2026-08-04: 2 of 5,396 molecules exceed 55 heavy atoms "
        "(max 59); both already rejected by `alerts`, neither on `shortlist` "
        "or `shortlist_synth`. Re-running LibInvent to stamp two rows that are "
        "already excluded is not worth a GPU day.",
    "t2_generate":
        "Same commit; ATRA generated 2026-07-27, the four newer seeds "
        "2026-07-31 (after). MEASURED 2026-08-04: no T_2 pool exceeds 49 heavy "
        "atoms (atra max 30), so the 55-atom ceiling would not have fired on "
        "any molecule in any pool. Zero impact.",
}


def _git(*args: str) -> str:
    return subprocess.check_output(["git", "-C", str(REPO), *args],
                                   text=True).strip()


def frame_commit(experiment: str, prefix: str, stage: str) -> tuple[str, dict]:
    """The commit recorded by the newest frame written by `stage`."""
    newest, meta = None, {}
    d = DATA_ROOT / experiment
    if not d.is_dir():
        return "", {}
    for mf in d.glob(f"{prefix}_*_manifest.json"):
        try:
            m = json.loads(mf.read_text())
        except Exception:  # noqa: BLE001
            continue
        if m.get("stage") != stage:
            continue
        n = int(mf.name.split("_")[1])
        if newest is None or n > newest:
            newest, meta = n, m
    git = (meta.get("git") or {})
    return str(git.get("commit", "")), {"n": newest, "dirty": git.get("dirty")}


def check() -> list[dict]:
    out = []
    for stage, (experiment, prefix, sources) in STAGES.items():
        commit, meta = frame_commit(experiment, prefix, stage)
        if not commit or commit == "unknown":
            out.append({"stage": stage, "state": "no-provenance",
                        "detail": "no frame, or the manifest records no commit"})
            continue
        try:
            later = _git("log", "--oneline", f"{commit}..HEAD", "--", *sources)
        except subprocess.CalledProcessError:
            out.append({"stage": stage, "state": "unknown-commit",
                        "detail": f"{commit[:12]} is not in this repo's history"})
            continue
        rec = {"stage": stage, "frame": f"{prefix}_{meta['n']}",
               "commit": commit[:12], "dirty": bool(meta.get("dirty")),
               "later": [ln for ln in later.splitlines() if ln]}
        rec["state"] = "stale" if rec["later"] else "current"
        out.append(rec)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--check", action="store_true",
                    help="exit non-zero on UNACKNOWLEDGED staleness")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    rows = check()
    unacknowledged = []
    for r in rows:
        tag = r["state"].upper()
        print(f"\n{r['stage']}  [{tag}]")
        if r["state"] in ("no-provenance", "unknown-commit"):
            print(f"  {r['detail']}")
            continue
        print(f"  frame {r['frame']} produced at {r['commit']}"
              + ("  (DIRTY TREE — the commit does not describe the code)"
                 if r["dirty"] else ""))
        for ln in r["later"]:
            print(f"    since: {ln}")
        if r["state"] == "stale":
            note = ACKNOWLEDGED.get(r["stage"])
            if note:
                print(f"  ACKNOWLEDGED: {note}")
            else:
                unacknowledged.append(r["stage"])
                print("  NOT ACKNOWLEDGED — re-run the stage, or record the "
                      "measured impact in ACKNOWLEDGED.")

    if args.check and unacknowledged:
        raise SystemExit(
            "frames produced by code that has since changed, with no recorded "
            f"impact: {unacknowledged}. Re-run the stage, or measure the "
            "impact and add it to ACKNOWLEDGED in "
            "scripts/check_frame_code_currency.py.")
    if args.check:
        print("\nno unacknowledged staleness.")


if __name__ == "__main__":
    main()
