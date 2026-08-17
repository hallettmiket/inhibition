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
from shared import run_paths as rp            # noqa: E402

log = logging.getLogger("control-stub")
B = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith")
REPORTS = rp.reports_dir()               # this run only (#74)
BOUND_NM = 1.2


def md_row(ident: str) -> pd.Series | None:
    """The molecule's 100 ns row, preferring one that succeeded."""
    parts = []
    for f in glob.glob(str(rp.residence_dir() / "*.csv")):  # this run only (#74)
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
    for f in sorted(glob.glob(str(rp.sweep_dir() / "attack_sweep_*.csv"))):  # this run only (#74)
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

    # THE SAME LAYOUT AS EVERY OTHER REPORT (@tt8804). A control that renders
    # differently from a candidate is harder to compare against, which is the only
    # reason it is on the rail. Same masthead, same verdict line, same panels --
    # what this page lacks is a trajectory, and that belongs in a panel of its
    # own rather than in a different page design.
    mast_facts = [(args.candidate, "molecule"),
                  ("control", "role"),
                  (f"{eng * 100:.2f}%", "100 ns engagement"),
                  (f"{rmax:.3f} nm", "max ligand RMSD")]
    if s is not None:
        mast_facts.append((f"{float(s.frac_attack_ready)*100:.1f}%  ·  "
                           f"{float(s.n_visits):.0f} visits", "attack-ready (10 ns)"))
    mast_facts.append((f"{float(m['production_ps']) / 1000:.0f} ns", "trajectory"))

    stand = (f"{'Held' if held else 'Left'}. Engaged the target in "
             f"{eng * 100:.2f}% of the 100 ns run.")

    detail = "".join(f"<tr><th>{k}</th><td>{v}</td></tr>" for k, v in rows)
    sweep_panel = ""
    if s is not None and float(s.n_visits) == 0:
        sweep_panel = (
            '<details class="panel"><summary>Rejected by the 10 ns sweep'
            '<span class="hint">and then produced the result above</span></summary>'
            '<div class="pbody"><p>'
            f"{float(s.frac_attack_ready):.4f} attack-ready, "
            f"<strong>zero sustained visits</strong>, so it would never have been "
            "elevated. D0076 shows the dwell filter discards exactly the brief "
            "approaches the observable was chosen to count; D0077 shows the "
            "crystal-reactant controls model an adduct as a Michaelis complex."
            "</p></div></details>")

    body = (
        rt.masthead(f"{args.candidate} — 100 ns residence", stand,
                    rt.eyebrow("CONTROL"), mast_facts)
        + f'<p>{rt.pill("Held" if held else "Left")} '
          f'Max ligand RMSD {rmax:.3f} nm against the {BOUND_NM} nm bar '
          f'&middot; {args.note}</p>'
        + '<details class="panel" open><summary>Measured values'
          '<span class="hint">the same readings every candidate carries</span>'
          f'</summary><div class="pbody"><table class="kv">{detail}</table></div></details>'
        + sweep_panel
        + '<details class="panel"><summary>No trajectory on this run'
          '<span class="hint">why there is no movie or RMSD plot</span></summary>'
          '<div class="pbody"><p>The 100 ns run completed and its measured row is '
          'intact, but it was launched without <code>--keep</code>, so the GROMACS '
          'working directory was cleaned up on completion. Every number on this '
          'page is read from the surviving row; none of it is inferred. Rebuild '
          'with <code>mdprio_report.py</code> once the run is repeated with '
          '<code>--keep</code> and this page becomes a full report.</p></div></details>')

    html = ("<!doctype html><html><head><meta charset='utf-8'>"
            f"<title>{args.candidate}</title><style>{rt.CSS}"
            "table.kv{border-collapse:collapse}"
            "table.kv th{text-align:left;padding:4px 16px 4px 0;color:var(--muted);"
            "font-weight:500}table.kv td{padding:4px 0;font-family:var(--mono)}"
            "</style></head><body>" + body + "</body></html>")

    dest = REPORTS / f"{args.candidate}.html"
    dest.write_text(html)
    print(f"  {args.candidate}: {eng*100:.2f}% engaged, "
          f"{'held' if held else 'left'} -> {dest}")


if __name__ == "__main__":
    main()
