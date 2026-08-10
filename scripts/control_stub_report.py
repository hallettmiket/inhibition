"""
Purpose: a report page for a control whose 100 ns numbers survived but whose trajectory did not.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-10
Input: --candidate <ident> (its row in 00_outputs/blacksmith/md_residence/)
Output: 00_outputs/blacksmith/mdprio_reports/<ident>.html

WHY THIS EXISTS, AND WHAT IT IS NOT. The three 100 ns control runs of 2026-08-09
were launched without `--keep`, so md_residence cleaned up each GROMACS workdir on
completion. The measured row was written first and is intact; `prod.xtc` is gone.
`mdprio_report.py` needs the trajectory -- RMSD series, attack geometry, movie --
so it cannot build these, and without a report page `mdprio_combine` drops the
molecule from the GUI entirely.

That would be the worst outcome: a control that RAN, and produced the result the
whole positive-control exercise was for, missing from the ranking because a
directory was deleted. This writes the page from the surviving row so the control
takes its place in the ranking on its real number.

IT IS LABELLED AS A STUB, LOUDLY AND ON THE PAGE. It carries no movie, no plots
and no pose, because those inputs no longer exist -- and a page that looked like
a full report while silently lacking them is exactly the populated-and-plausible
failure this project keeps finding. Rebuild properly with `mdprio_report.py` once
the trajectory is re-run.
"""

from __future__ import annotations

import argparse
import glob
import logging
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import report_theme as rt                   # noqa: E402

log = logging.getLogger("control-stub")
B = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith")
REPORTS = B / "mdprio_reports"
BOUND_NM = 1.2


def md_row(ident: str) -> pd.Series | None:
    """The molecule's 100 ns row, preferring one that succeeded."""
    parts = []
    for f in glob.glob(str(B / "md_residence/*.csv")):
        try:
            parts.append(pd.read_csv(f))
        except Exception:                                # noqa: BLE001
            continue
    if not parts:
        return None
    d = pd.concat(parts, ignore_index=True)
    d = d[(d.ident.astype(str) == ident) & (d.get("production_ps", 0) >= 50000)]
    if "status" in d.columns:
        d = d[d.status.astype(str).str.startswith("ok")]
    e = "explicit_frac_frames_engaged"
    if e in d.columns:
        d = d[d[e].notna()]
    return None if d.empty else d.iloc[-1]


def sweep_row(ident: str) -> pd.Series | None:
    parts = []
    for f in sorted(glob.glob(str(B / "attack_sweep/attack_sweep_*.csv"))):
        try:
            parts.append(pd.read_csv(f))
        except Exception:                                # noqa: BLE001
            continue
    if not parts:
        return None
    d = pd.concat(parts, ignore_index=True)
    d = d[(d.parent_ident.astype(str) == ident) & (d.status == "ok")]
    return None if d.empty else d.sort_values("frac_attack_ready").iloc[-1]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--candidate", required=True)
    ap.add_argument("--note", default="", help="one line on this control's provenance")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    m = md_row(args.candidate)
    if m is None:
        raise SystemExit(f"{args.candidate}: no successful 100 ns row to build from")
    s = sweep_row(args.candidate)

    eng = float(m["explicit_frac_frames_engaged"])
    rmax = float(m["explicit_ligand_rmsd_nm_max"])
    held = rmax < BOUND_NM
    rows = [("100 ns target engagement", f"{eng * 100:.2f}%"),
            ("max ligand RMSD", f"{rmax:.3f} nm"),
            ("verdict", "held" if held else f"left (above the {BOUND_NM} nm bar)"),
            ("production", f"{float(m['production_ps']) / 1000:.0f} ns")]
    for k, lbl in (("explicit_ligand_rmsd_nm_mean", "mean ligand RMSD"),
                   ("explicit_ligand_rmsd_nm_final", "final ligand RMSD")):
        if k in m and pd.notna(m[k]):
            rows.append((lbl, f"{float(m[k]):.3f} nm"))
    if s is not None:
        rows.append(("10 ns sweep, attack-ready", f"{float(s.frac_attack_ready):.4f}"))
        rows.append(("10 ns sweep, sustained visits", f"{float(s.n_visits):.0f}"))

    facts = "".join(f"<tr><th>{k}</th><td>{v}</td></tr>" for k, v in rows)
    sweep_line = ""
    if s is not None and float(s.n_visits) == 0:
        sweep_line = (
            "<p class='bad'><strong>This control was rejected by the 10 ns sweep</strong> "
            f"— {float(s.frac_attack_ready):.4f} attack-ready, zero sustained visits — "
            "and would never have been elevated. It then produced the 100 ns result "
            "above. See <code>D0075</code>.</p>")

    body = f"""
<h1>{args.candidate}</h1>
<p class="lead">{args.note or 'Control, 100 ns.'}</p>

<div class="stub"><strong>Numbers only — this is not a full report.</strong>
The 100 ns run completed and its measured row is intact, but it was launched
without <code>--keep</code>, so the GROMACS working directory was cleaned up on
completion. <strong>There is no trajectory, so there is no movie, no RMSD plot and
no pose here.</strong> Everything below is read from the surviving row. Rebuild
with <code>mdprio_report.py</code> once the run is repeated with
<code>--keep</code>.</p></div>

{sweep_line}

<table class="facts">{facts}</table>

<p class="note">Held/left uses the same {BOUND_NM} nm maximum-ligand-RMSD bar as
every candidate, so this row is directly comparable in the ranking.</p>
"""
    html = ("<!doctype html><html><head><meta charset='utf-8'>"
            f"<title>{args.candidate}</title><style>{rt.CSS}"
            ".stub{border-left:4px solid #8a5a00;background:#fdf0dc;padding:10px 14px;"
            "border-radius:0 4px 4px 0;margin:1rem 0}"
            ".bad{border-left:4px solid #b3261e;background:#fbeae8;padding:10px 14px;"
            "border-radius:0 4px 4px 0}"
            "table.facts{border-collapse:collapse;margin-top:1rem}"
            "table.facts th{text-align:left;padding:5px 14px 5px 0;color:#5b6b80;"
            "font-weight:500}table.facts td{padding:5px 0;font-family:ui-monospace,monospace}"
            ".lead{color:#5b6b80}.note{color:#5b6b80;font-size:.9rem}"
            "</style></head><body>" + body + "</body></html>")

    dest = REPORTS / f"{args.candidate}.html"
    dest.write_text(html)
    print(f"  {args.candidate}: {eng*100:.2f}% engaged, "
          f"{'held' if held else 'left'} -> {dest}")


if __name__ == "__main__":
    main()
