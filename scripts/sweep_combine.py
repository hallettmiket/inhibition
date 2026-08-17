#!/usr/bin/env python3
"""
Purpose: the Sweep results page — the MD results layout, driven by sweep data.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-16
Input: --worklist <sweep_gaps_N.csv>
Output: mdprio_reports/sweep.html

@tt8804: "just copy over the design for MD results and add/modify the tables to
use the sweep result data", and "the sweep results should be ranked in the order
of 100 ns MD priority".

THE LAYOUT IS NOT REDESIGNED, IT IS THE SAME OBJECT. `shared/results_shell.CSS`
is the stylesheet `mdprio_combine` uses, moved there so both interpolate one
constant: a 376px rail of rows on the left, one report at a time in an iframe on
the right. Three bespoke sweep layouts were three too many.

RANKED ON WHETHER IT STAYED PUT. @tt8804: "I don't care about visits I care if
they stayed in a min RMSD range", and the case that settles it: a mode 75.8%
attack-ready that flies out of the pocket. Attack geometry measured on a molecule
that has left is geometry against a site it is no longer in.

So the rail ranks on MAX ligand RMSD over the sweep, lowest first -- the same
headline the 100 ns results page ranks on -- and marks held/left at
`md.sweep_survivor_rmsd_nm`, the bar that actually gates the next stage.
Attack-ready sits beside it as the triage reading it is.

This supersedes ranking on `n_visits` (prereg T4). T4's argument -- one good
approach beats sustained occupancy -- is about WHICH attack observable to prefer,
and it is sound; it presumes the molecule is still bound, which is the thing this
order tests first.

The rail therefore answers one question: what earns a 100 ns run next, best first.
"""

from __future__ import annotations

import argparse
import glob
import html
import logging
import os
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import gui_shell as gs                        # noqa: E402
from shared import results_shell as rs                    # noqa: E402
from shared import sweep_state as ss                      # noqa: E402
from shared import run_paths as rp                  # noqa: E402

log = logging.getLogger("sweep-combine")
B = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith")
REPORTS = rp.reports_dir()
PAGES = REPORTS / "sweep_pages"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--worklist", required=True)
    ap.add_argument("--title", default="DWI covalent screen")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    wl = Path(args.worklist)
    d = ss.state(wl)
    summ = ss.summary(d)
    pred = ss.predicts(d)
    ok = d[d.sweep_state == "ok"].copy()
    # THIS CAMPAIGN, not every sweep ever run. `state()` also carries results
    # from earlier worklists, other warhead classes and pre-tier-decision T_3
    # rows; listing them here would rank modes selected under rules that are no
    # longer in force, and the "qualify" count would describe a set nobody chose.
    if "_queued" in ok.columns and bool(ok["_queued"].any()):
        ok = ok[ok["_queued"]]

    # ---- MD PRIORITY ORDER: DID IT STAY PUT ---------------------------------
    # @tt8804: "I don't care about visits I care if they stayed in a min RMSD
    # range". So the rail ranks on MAX ligand RMSD over the 10 ns, lowest first --
    # how far the molecule ever got from where it started -- which is exactly the
    # headline the 100 ns results page ranks on. A molecule that wandered out of
    # the pocket has answered the question whatever its attack geometry did.
    #
    # THE BAR IS THE SWEEP'S OWN GATE, not the 100 ns one. This read
    # `md.bound_rmsd_nm` (1.2 nm) -- the old "did not dissociate" reading -- so
    # the page marked 134 of 168 modes "held" while the rule that actually
    # decides what earns a 100 ns run passes 12. A page whose badge disagrees
    # with the gate is worse than a page with no badge.
    import sweep_assets as sa
    from shared import target_config as tc
    bound = float(tc.get("md.sweep_survivor_rmsd_nm", default=0.35))
    # POSE RANK IS PART OF THE KEY. Without it every mode of a molecule resolves
    # to the same trajectory: three modes of t4_59100a17abd6 all reported max
    # RMSD 0.347 nm when they are 0.416, 0.663 and 6.62 -- the last of which
    # leaves the pocket entirely and was sitting in the held list.
    wlf = pd.read_csv(wl)
    prank = (dict(zip(wlf.ident.astype(str), wlf.pose_rank.astype(int)))
             if "pose_rank" in wlf.columns else {})
    mx, mean = {}, {}
    for ident in ok.ident.astype(str):
        rep = sa.rep_dir(ident.rsplit("_m", 1)[0], prank.get(ident))
        if rep is None:
            continue
        _t, y = sa._xvg(rep / "rmsd.xvg")
        if y is not None and len(y):
            mx[ident], mean[ident] = float(y.max()), float(y.mean())
    ok["rmsd_max"] = ok.ident.astype(str).map(mx)
    ok["rmsd_mean"] = ok.ident.astype(str).map(mean)
    for c in ("n_visits", "frac_attack_ready"):
        if c not in ok.columns:
            ok[c] = float("nan")
    # Lowest max RMSD first. A mode with no RMSD trace sorts last rather than
    # first, which a NaN would otherwise do silently under ascending order.
    ok = ok.sort_values("rmsd_max", ascending=True, na_position="last")
    log.info("ranked on max ligand RMSD; %d of %d have a trace, held bar %.2f nm",
             int(ok.rmsd_max.notna().sum()), len(ok), bound)
    pending = d[d.sweep_state.isin(("pending", "failed"))]

    thumbs = {}
    try:
        import importlib.util as _u
        sp = _u.spec_from_file_location("mc", REPO / "scripts/mdprio_combine.py")
        mc = _u.module_from_spec(sp); sp.loader.exec_module(mc)
        thumbs = mc._thumbs([str(i).rsplit("_m", 1)[0] for i in ok.ident])
    except Exception as exc:                               # noqa: BLE001
        log.warning("depictions unavailable: %s", exc)

    rows, tabs = [], []
    for i, r in enumerate(ok.itertuples()):
        ident = str(r.ident)
        parent = ident.rsplit("_m", 1)[0]
        if not (PAGES / f"{ident}.html").is_file():
            continue
        tabs.append(ident)
        vis = 0 if pd.isna(r.n_visits) else int(r.n_visits)
        ar = 0.0 if pd.isna(r.frac_attack_ready) else float(r.frac_attack_ready)
        rmx = None if pd.isna(getattr(r, "rmsd_max", float("nan"))) else float(r.rmsd_max)
        # THE HEADLINE IS THE NUMBER THE RAIL IS SORTED ON, exactly as on the MD
        # results page. Putting anything else here invites reading the order as
        # something it is not.
        headline = f"{rmx:.3f} nm max" if rmx is not None else "—"
        held = rmx is not None and rmx < bound
        meta = ((("held" if held else "left") + f" · {ar*100:.0f}% attack-ready")
                if rmx is not None else f"{ar*100:.0f}% attack-ready")
        cls = html.escape(str(getattr(r, "warhead_class", "") or "—"))
        th = thumbs.get(parent)
        img = (f"<img class='thumb' loading='lazy' alt='' src='{th}'>" if th
               else "<span class='thumb'></span>")
        rows.append(
            f"<button class='row{'' if held else ' left'}' data-cls=\"{cls}\" "
            f"id='b_{html.escape(ident)}' "
            f"onclick=\"show('{html.escape(ident)}')\">"
            f"<span class='rk'>{i+1}</span>{img}"
            f"<span class='body'><span class='l1'>"
            f"<span class='mid-id'>{html.escape(ident)}</span>"
            f"<span class='eng'>{headline}</span></span>"
            f"<span class='l2'><span class='wc'>{cls}</span>"
            f"<span class='meta'>{meta}</span></span>"
            f"<span class='bar'><i style='width:"
            f"{(100*max(0.0,1-(rmx/bound)) if rmx is not None else 0):.0f}%'></i></span>"
            f"</span></button>")

    first = tabs[0] if tabs else ""
    # AN EMPTY RAIL MUST NOT PRODUCE AN IFRAME. With no finished sweep the src
    # interpolated to `sweep_pages/.html` -- a request for the empty ident --
    # and the viewer pane rendered the server's 404 page inside the layout,
    # which reads as a broken GUI rather than as a run with no results yet.
    # @tt8804: "showing a 404".
    _viewer = (f'<iframe id="frame" src="sweep_pages/{html.escape(first)}.html"'
               f' title="sweep report"></iframe>' if first else
               '<div class="legend" style="margin:22px">No sweep has finished '
               'yet, so there is no report to show. Modes appear in the rail as '
               'they come back, and this panel loads the one you select.</div>')
    n_pri = int((ok.rmsd_max < bound).sum())
    rho = pred.get("rho")
    note = (f"Spearman(enrichment, attack-ready) = {rho:+.3f}, p = {pred['p']:.3f} "
            f"over {pred['n']} finished modes — the docked ranking does not "
            f"predict these results, so this order is the empirical one."
            if rho is not None else "")

    page = f"""<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(args.title)} — Sweep results</title><style>
{rs.CSS}
{gs.CSS}
</style></head><body>
<div id="topbar">
 <h1 title="Pick a mode on the left; its plots, movie and readings load beside it.">
 {html.escape(args.title)} — sweep results</h1>
 <span class="msep"></span>
 <span class="mhint" id="mhint">{len(tabs)} modes &middot; ranked by 100&nbsp;ns
 priority &middot; {n_pri} qualify</span>
 <span class="msep"></span>
 <a class="mbtn lnk" href="pipeline.html">how this works &#8599;</a>
 <button id="theme" class="mbtn tbtn" onclick="toggleTheme()">dark</button>
</div>
{gs.nav("sweep.html", {"sweep.html": f"{summ['ok']} ok · {summ['pending']} pending"})}
<main>
 <div id="rail">
  {rs.SEARCH_HTML}
  <div class="legend">ranked by <b>max ligand RMSD</b> over the 10&nbsp;ns run,
   lowest first &mdash; how far the molecule ever got from where it started, the
   same headline the 100&nbsp;ns results page ranks on. <b>held</b> means it never
   exceeded {bound:.2f}&nbsp;nm. Attack-ready % is shown beside it as the triage
   reading, not as the order. {note}</div>
  {''.join(rows)}
  {'<div class="legend">' + str(len(pending)) + ' modes pending or failed — not shown, they have no readings yet.</div>' if len(pending) else ''}
 </div>
 <div id="viewer">
  {_viewer}
 </div>
</main>
<script>
let CUR = "{html.escape(first)}";
function show(id){{
  if(!id) return;
  CUR = id;
  const f = document.getElementById('frame');
  if(!f) {{ location.reload(); return; }}
  f.src = 'sweep_pages/' + id + '.html';
  document.querySelectorAll('#rail .row').forEach(function(b){{
    b.classList.toggle('on', b.id === 'b_' + id); }});
}}
function toggleTheme(){{
  const d = document.documentElement.getAttribute('data-theme') === 'dark';
  document.documentElement.setAttribute('data-theme', d ? 'light' : 'dark'); }}
show(CUR);
{rs.SEARCH_JS}
railFilter();
</script>
</body></html>"""
    (REPORTS / "sweep.html").write_text(page)
    print(f"\n  {len(tabs)} sweep reports, ranked by max ligand RMSD "
          f"({n_pri} held under {bound:.2f} nm) -> {REPORTS / 'sweep.html'}")


if __name__ == "__main__":
    main()
