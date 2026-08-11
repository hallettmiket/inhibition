#!/usr/bin/env python3
"""
Purpose: one molecule, its independent 100 ns replicates side by side.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-10
Input: --candidate <ident> --run <dir>:<label> ... (each a finished rep<N> dir)
Output: 00_outputs/blacksmith/shortlist/replicates_<N>.html

The shortlist format shows one trajectory. One trajectory cannot tell you whether
what you are looking at is the molecule or the seed. This shows the replicates
together: the comparison table first, then the three RMSD panels in a row, then
each run in full in the ordinary shortlist format.

WHAT A SPREAD ACROSS THESE COLUMNS MEANS, AND WHAT IT DOES NOT. These are
independent velocity seeds on ONE pose of ONE complex. Agreement says the
behaviour is not an artefact of a single seed. It does not make the number a
prediction of activity, it does not feed a gate, and it reorders no shortlist
(D0036, D0044). The uncertainty printed on a single run is within-trajectory and
autocorrelation-corrected; the spread across columns is the replicate spread, and
the two are different quantities that must not be mixed.
"""

from __future__ import annotations

import argparse
import glob
import logging
import sys
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

import shortlist_report as sr                              # noqa: E402
from shared import md_movie as mov                         # noqa: E402
from shared import outputs as sout                        # noqa: E402
from shared import report_theme as rt                      # noqa: E402

log = logging.getLogger("replicates")


def _mp():
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "mdprio_report", REPO / "scripts" / "mdprio_report.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules["mdprio_report"] = m
    spec.loader.exec_module(m)
    return m


def engagement_for(rep: Path) -> float | None:
    """Target engagement for THIS run, found by the directory it was measured from.

    Not by filename, and not by position in a glob. Every md_residence row
    carries the directory the numbers came from -- `equilibration_dir` for a run
    measured in place, `remeasured_from` for one rescued afterwards -- so the
    number can be tied to the trajectory rather than to a name that happens to
    contain "rep2". This project's recurring defect is a value taken by name or
    position instead of by identity; three near-identical replicate rows in one
    directory is exactly where that bites.
    """
    import pandas as pd
    want = str(rep).rstrip("/")
    best = None
    for f in sorted(glob.glob(str(sr.B / "md_residence" / "*.csv"))):
        try:
            d = pd.read_csv(f)
        except Exception:                                  # noqa: BLE001
            continue
        cols = [c for c in ("equilibration_dir", "remeasured_from", "trajectory")
                if c in d.columns]
        if not cols or "explicit_frac_frames_engaged" not in d.columns:
            continue
        for _, row in d.iterrows():
            if any(str(row[c]).rstrip("/").startswith(want) for c in cols):
                v = row["explicit_frac_frames_engaged"]
                if pd.notna(v):
                    best = float(v)
    return best


def summarise(rep: Path, er, mp) -> dict | None:
    """The measured values for one run, or None if it has no usable trajectory."""
    if not rep.is_dir():
        log.warning("no directory: %s", rep)
        return None
    total_ns = mp.prod_ns(rep)
    s = mp.series(rep, er, total_ns)
    res = mp.residence(s)
    if res.get("status") != "ok":
        log.warning("%s: %s", rep, res.get("status"))
        return None
    res["length_ns"] = res.get("length_ns", total_ns)
    return res


def comparison(runs: list[tuple[str, Path, dict]]) -> str:
    """The table that is the point of the document.

    Rows are metrics, columns are runs, so the eye travels ACROSS a row to see
    whether the runs agree -- which is the only comparison worth making here.
    Laid out the other way round the reader has to hold three numbers in their
    head to answer one question.
    """
    def cells(fn):
        return "".join(f"<td>{fn(r)}</td>" for _, _, r in runs)

    def held(r):
        left = r.get("left_at_ns")
        return ('<span class="bad">left</span>' if left is not None
                else '<span class="good">held</span>')

    def eng(r):
        v = r.get("engaged")
        return f"{float(v) * 100:.2f}%" if v is not None else "&mdash;"

    body = [
        ("target engagement, 100 ns", cells(eng)),
        ("verdict", cells(held)),
        ("mean ligand RMSD", cells(lambda r: f"{r['rmsd_mean_nm']:.3f} nm")),
        ("max ligand RMSD", cells(lambda r: f"{r['rmsd_max_nm']:.3f} nm")),
        ("final ligand RMSD", cells(lambda r: f"{r['rmsd_final_nm']:.3f} nm")),
        ("residence fraction", cells(lambda r: f"{r['residence_frac']:.3f}")),
        ("trajectory", cells(lambda r: f"{r['length_ns']:.1f} ns, "
                                       f"{r['n_frames']:,} frames")),
    ]
    head = "".join(f"<th>{lab}</th>" for lab, _, _ in runs)
    rows = "".join(f"<tr><th>{k}</th>{v}</tr>" for k, v in body)
    return (f'<table class="cmp"><thead><tr><th></th>{head}</tr></thead>'
            f"<tbody>{rows}</tbody></table>")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--candidate", required=True)
    ap.add_argument("--run", nargs="+", required=True,
                    help="<dir>:<label>, one per replicate, in the order to show")
    ap.add_argument("--name", default=None)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    import importlib.util as _u
    _sp = _u.spec_from_file_location("mdprio_combine",
                                     REPO / "scripts" / "mdprio_combine.py")
    _mc = _u.module_from_spec(_sp); _sp.loader.exec_module(_mc)
    ver, code = _mc._version()

    er, mp = sr._er(), _mp()
    three = (REPO / "scripts/.cache_3dmol-min.js").read_text()
    cls = sr.classes()

    runs, figs = [], []
    for spec in args.run:
        d, _, lab = spec.rpartition(":")
        rep = Path(d)
        res = summarise(rep, er, mp)
        if res is None:
            log.warning("%s: skipped", lab)
            continue
        res["engaged"] = engagement_for(rep)
        runs.append((lab, rep, res))

    if len(runs) < 2:
        raise SystemExit("need at least two usable runs to compare")

    # The three RMSD panels in a row -- the literal side-by-side. Same y-scales
    # come from mp.figure's own conventions; nothing here rescales them, because
    # three axes that differ silently is worse than three that are simply wide.
    for lab, rep, res in runs:
        total_ns = res["length_ns"]
        s = mp.series(rep, er, total_ns)
        mpdb = rep / "movie.pdb"
        nacs = mp.nac_series(args.candidate, rep, mpdb, total_ns) \
            if mpdb.is_file() else None
        img = mp.figure(args.candidate, s, res, er, nacs)
        figs.append(f'<figure><img src="data:image/png;base64,{img}" alt="">'
                    f"<figcaption>{lab}</figcaption></figure>")

    # Each run in full, in the ordinary shortlist format. Element ids are
    # namespaced per run: getElementById returns the FIRST match, so three
    # sections sharing ids means every control drives the first viewer.
    blocks = []
    for i, (lab, rep, _res) in enumerate(runs, 1):
        b = sr.block(args.candidate, er, three, cls,
                     rep_dir=rep, suffix=f"_r{i}",
                     heading=f"{args.candidate} &mdash; {lab}")
        if b:
            blocks.append(b)

    name = args.name or f"{args.candidate} — {len(runs)} × 100 ns replicates"
    title = f"{date.today().isoformat()} {name}"
    byline_ver = " ".join(x for x in (f"version {ver}" if ver else "",
                                      f"“{code}”" if code else "") if x)

    page = f"""<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title}</title><style>{rt.CSS}{mov.VIEWER_CSS}
body{{max-width:1180px;margin:0 auto;padding:0 30px 70px}}
body>*{{max-width:none;padding-left:0;padding-right:0}}
table.cmp{{border-collapse:collapse;width:100%;font-size:13.5px;margin:1rem 0 1.6rem}}
table.cmp th,table.cmp td{{padding:.45rem .9rem;text-align:left;
  border-bottom:1px solid var(--rule)}}
table.cmp thead th{{font:600 12px var(--sans);color:var(--muted);
  border-bottom:2px solid var(--rule);text-transform:uppercase;
  letter-spacing:.05em}}
table.cmp tbody th{{font-weight:500;color:var(--muted);white-space:nowrap}}
table.cmp td{{font-family:var(--mono)}}
.good{{color:#0f7a54;font-weight:600}}
.bad{{color:#b3261e;font-weight:600}}
.figrow{{display:grid;grid-template-columns:repeat({len(figs)},1fr);gap:14px;
  margin:1rem 0 1.8rem}}
.figrow figure{{margin:0}}
.figrow img{{width:100%;height:auto;border:1px solid var(--rule);
  border-radius:5px;background:#fff}}
.figrow figcaption{{font:600 12px var(--sans);color:var(--muted);
  margin-top:.4rem;text-align:center}}
</style></head><body>
<header class="mast"><h1>{title}</h1>
<p class="standfirst">Timothy Wu &middot; {byline_ver}</p></header>
<p>Independent velocity seeds on one pose of one complex. Agreement across the
columns says the behaviour is not an artefact of a single seed. It is not a
prediction of activity, it feeds no gate, and it reorders no shortlist.</p>
{comparison(runs)}
<div class="figrow">{''.join(figs)}</div>
<script>{three}</script>
{''.join(blocks)}
</body></html>"""

    dest = sout.Topic("blacksmith", "shortlist").write("replicates", ".html")
    dest.write_text(page)
    stable = dest.parent / "replicates.html"
    stable.write_text(page)
    print(f"\n  {len(runs)} runs -> {stable}  ({len(page)/1048576:.1f} MB)")
    print(f"  versioned: {dest}")


if __name__ == "__main__":
    main()
