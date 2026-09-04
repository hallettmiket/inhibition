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

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import gui_shell as gs                        # noqa: E402
from shared import results_shell as rs                    # noqa: E402
from shared import sweep_state as ss                      # noqa: E402
from shared import run_paths as rp                  # noqa: E402

log = logging.getLogger("sweep-combine")
B = rp.BLACKSMITH
REPORTS = rp.reports_dir()
PAGES = REPORTS / "sweep_pages"
#: Sweep length, derived -- this legend said 10 ns while the sweep has run
#: at 8 ns since D0085 (@tt8804).
from shared import target_config as _tc                   # noqa: E402
_SWEEP_NS = int(round(_tc.md_sweep_ps() / 1000))


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
    # KEYED ON `task_id`, WHICH IS THE MODE. The block above says pose rank is
    # part of the key and then built the map from `ident`, which is the
    # MOLECULE (`t4_215b12bd9b34`) while every lookup below passes a mode
    # (`t4_215b12bd9b34_m184`). Measured 2026-09-02: 0 of 98 finished modes
    # resolved a pose rank, so `rep_dir` fell back to whichever sibling sorted
    # first for every single row -- the exact defect the comment was written to
    # prevent, reintroduced one column name away. `sweep_assets` and
    # `recompute_attack_ready` already key on `task_id`; this was the last
    # reader that did not.
    _k = "task_id" if "task_id" in wlf.columns else "ident"
    prank = (dict(zip(wlf[_k].astype(str), wlf.pose_rank.astype(int)))
             if "pose_rank" in wlf.columns else {})
    if not prank:
        log.warning("no pose_rank map from %s; RMSD falls back to the first "
                    "sibling directory per molecule", wl.name)

    # PREFER THE RESULTS TABLE. `attack_sweep` now writes `rmsd_max_a` /
    # `rmsd_mean_a` (ANGSTROM) on every completed row, computed by the run that
    # owns the trajectory. Recomputing them here was a second implementation of
    # the same quantity reading the same files by a different route -- which is
    # how a page and a gate come to disagree while both look right. The xvg
    # walk below is now only the fallback for rows written before that column
    # existed.
    mx, mean = {}, {}
    for c_a, dst in (("rmsd_max_a", mx), ("rmsd_mean_a", mean)):
        if c_a in ok.columns:
            for i, v in zip(ok.ident.astype(str), ok[c_a]):
                if pd.notna(v):
                    dst[i] = float(v) / 10.0              # Angstrom -> nm
    missing = [i for i in ok.ident.astype(str) if i not in mx]
    for ident in missing:
        rep = sa.rep_dir(ident.rsplit("_m", 1)[0], prank.get(ident))
        if rep is None:
            continue
        _t, y = sa._xvg(rep / "rmsd.xvg")
        if y is not None and len(y):
            mx[ident], mean[ident] = float(y.max()), float(y.mean())
    if missing:
        log.info("%d of %d rows had no stored RMSD; read from rmsd.xvg",
                 len(missing), len(ok))
    ok["rmsd_max"] = ok.ident.astype(str).map(mx)
    ok["rmsd_mean"] = ok.ident.astype(str).map(mean)
    for c in ("n_visits", "frac_attack_ready"):
        if c not in ok.columns:
            ok[c] = float("nan")
    # ---- PRIORITY = THE 100 ns GATE, not one of its two halves --------------
    # @twu383, 2026-09-02: *"when the 1.2 ns finishes this will decide elevating
    # to 100 ns"*. The rail used to rank on max ligand RMSD alone, which is the
    # POSE-STABILITY half; a molecule that sits perfectly still facing the wrong
    # way topped the list. The gate is engagement AND stability, so the order is
    # too:
    #
    #   tier 0  clears the gate           (engaged >= bar AND pose held)
    #   tier 1  pose held, not engaged    -- the near misses, best engagement first
    #   tier 2  pose left                 -- ranked by engagement anyway, so a
    #                                        mode that engaged before departing
    #                                        is visible rather than buried
    #
    # Within every tier: engagement descending, then max RMSD ascending. A mode
    # with no RMSD trace sorts last rather than first, which a NaN does silently
    # under ascending order.
    import attack_sweep as _asw
    _th = _asw.elevation_thresholds()
    occ_bar = _th["occupancy_min"]
    occ_window = float(getattr(_asw, "COMMON_WINDOW_PS", 1200.0))

    def _held(r):
        mxa, mna = r.get("rmsd_max_a"), r.get("rmsd_mean_a")
        if pd.isna(mxa) or pd.isna(mna):
            # fall back to the nm value the rail already has, max-bar only
            v = r.get("rmsd_max")
            return bool(pd.notna(v) and float(v) * 10.0 < _th["rmsd_max_a"])
        return (float(mxa) < _th["rmsd_max_a"]) or (float(mna) < _th["rmsd_mean_a"])

    # NOT `_held`. `DataFrame.itertuples` renames any column whose name
    # starts with an underscore to a POSITIONAL name (`_5`), so
    # `getattr(r, "_held")` in the row loop below silently returned the
    # default for every row and every mode was labelled "left" -- including
    # ones at 0.189 nm. The DataFrame-side logic was right; only the
    # rendering read a name that no longer existed.
    ok["gate_held"] = ok.apply(_held, axis=1)
    # THE COMPARABLE FIGURE RANKS, ALWAYS. `frac_attack_ready` is a fraction of
    # whatever the run actually ran, and since 2026-09-03 that varies per mode
    # (adaptive length: stop when the molecule leaves, cap 10 ns). Ranking on it
    # would put a 1.2 ns run's 33.9% above a 10 ns run's 26.7% as though they
    # were the same measurement -- two populations under one column, which is
    # the defect `frac_attack_ready_common` exists to prevent. That column is
    # the first 1.2 ns on EVERY row, and equals the full-run figure on the 702
    # fixed-length rows, so nothing already ranked moves for the wrong reason.
    if "frac_attack_ready_common" in ok.columns:
        eng = ok.frac_attack_ready_common.astype(float)
        # A row with no common-window reading falls back to the full-run figure
        # ONLY when the run is exactly one window long; otherwise it is left out
        # of the ranking rather than given a number from a different window.
        gap = eng.isna() & (ok.get("sweep_ps", pd.Series(index=ok.index)) == occ_window)
        eng = eng.where(~gap, ok.frac_attack_ready.astype(float))
        ok["gate_eng"] = eng.fillna(0.0)
        ok["eng_basis"] = "first 1.2 ns"
    else:
        ok["gate_eng"] = ok.frac_attack_ready.fillna(0.0).astype(float)
        ok["eng_basis"] = "full run"
    ok["gate_pass"] = ok["gate_held"] & (ok["gate_eng"] >= occ_bar)
    ok["gate_tier"] = np.where(ok["gate_pass"], 0, np.where(ok["gate_held"], 1, 2))
    ok = ok.sort_values(["gate_tier", "gate_eng", "rmsd_max"],
                        ascending=[True, False, True], na_position="last")
    log.info("ranked on the 100 ns gate: %d clear it, %d held but under the "
             "%.0f%% engagement bar, %d left the site",
             int((ok.gate_tier == 0).sum()), int((ok.gate_tier == 1).sum()),
             occ_bar * 100, int((ok.gate_tier == 2).sum()))
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
        # THE HEADLINE IS THE NUMBER THE RAIL IS SORTED ON. The rail is sorted
        # on the gate, whose leading term is engagement, so engagement is the
        # headline and stability sits beside it. Showing max RMSD here while
        # sorting on something else is how a reader infers the wrong order.
        # The headline is the COMPARABLE figure, because that is what the rail
        # is sorted on. Showing the full-run percentage here while ordering on
        # the common window is how a reader infers the wrong ranking.
        ar = float(getattr(r, "gate_eng", ar))
        headline = f"{ar*100:.0f}% engaged"
        held = bool(getattr(r, "gate_held", False))
        passes = bool(getattr(r, "gate_pass", False))
        rtxt = f"{rmx:.3f} nm max" if rmx is not None else "no trace"
        # HOW LONG IT RAN, AND WHETHER IT LEFT. Under adaptive length the run
        # length IS a result -- a mode that held for 10 ns and one that left at
        # 1.2 are not the same finding, and without this the page shows two
        # engagement figures with no way to tell which window produced them.
        sps = getattr(r, "sweep_ps", None)
        ltxt = ""
        if sps is not None and not pd.isna(sps):
            left_site = getattr(r, "left_site", None)
            if left_site is not None and not pd.isna(left_site):
                ltxt = (f" · left at {float(sps)/1000:.1f} ns" if bool(left_site)
                        else f" · held {float(sps)/1000:.0f} ns")
            elif float(sps) != occ_window:
                ltxt = f" · {float(sps)/1000:.1f} ns"
        meta = f"{'held' if held else 'left'} · {rtxt}{ltxt}"
        if passes:
            meta = "ELEVATE · " + meta
        cls = html.escape(str(getattr(r, "warhead_class", "") or "—"))
        th = thumbs.get(parent)
        img = (f"<img class='thumb' loading='lazy' alt='' src='{th}'>" if th
               else "<span class='thumb'></span>")
        rows.append(
            f"<button class='row{'' if held else ' left'}"
            f"{' pass' if passes else ''}' data-cls=\"{cls}\" "
            f"id='b_{html.escape(ident)}' "
            f"onclick=\"show('{html.escape(ident)}')\">"
            f"<span class='rk'>{i+1}</span>{img}"
            f"<span class='body'><span class='l1'>"
            f"<span class='mid-id'>{html.escape(ident)}</span>"
            f"<span class='eng'>{headline}</span></span>"
            f"<span class='l2'><span class='wc'>{cls}</span>"
            f"<span class='meta'>{meta}</span></span>"
            # THE BAR IS PROGRESS TOWARDS THE GATE, so a full bar means
            # elevated rather than "lowest RMSD seen so far".
            f"<span class='bar'><i style='width:"
            f"{min(100.0, 100.0*ar/occ_bar):.0f}%'></i></span>"
            f"</span></button>")

    # Values the legend states, taken from the gate itself rather than typed
    # into the prose -- a page that describes a threshold it does not read is
    # how the caption and the order come to disagree.
    hi_a = _asw.attack_ready_max_a()
    occ_pct = occ_bar * 100
    rmx_a, rmn_a = _th["rmsd_max_a"], _th["rmsd_mean_a"]
    n_pass = int((ok.gate_tier == 0).sum())
    n_near = int((ok.gate_tier == 1).sum())
    n_left = int((ok.gate_tier == 2).sum())
    n_ok = len(ok)
    # HOW MANY ROWS ARE OF EACH KIND, for the legend. A page that mixes
    # fixed-length and adaptive rows has to say so, and say how many.
    _sps = ok.get("sweep_ps")
    n_fixed = int((_sps == occ_window).sum()) if _sps is not None else n_ok
    n_adaptive = int(ok.get("adaptive", pd.Series(False, index=ok.index))
                     .fillna(False).astype(bool).sum())

    # A PAGE THAT CANNOT BE READ IS NOT A SWEEP THAT DID NOT RUN.
    #
    # `tabs` only collects idents whose `sweep_pages/<ident>.html` exists, so if
    # those are absent, being rewritten, or unreadable (the Isilon ACL does not
    # honour the POSIX mode it advertises -- see CLAUDE.md), this rendered "No
    # sweep has finished yet" over a run with hundreds of results. That is the
    # catalogue's disguise #4: a message naming the wrong cause, and it is the
    # same shape as `seed_status` reporting "no frame written yet" for a
    # permission error.
    #
    # Refusing to write is the right answer: the LAST GOOD page stays up, and
    # the exception says what is actually wrong. Writing a claim that the
    # campaign has produced nothing is the worst thing this page can get wrong.
    if not tabs and not ok.empty:
        missing = [str(r.ident) for r in ok.itertuples()
                   if not (PAGES / f"{str(r.ident)}.html").is_file()][:5]
        raise SystemExit(
            f"{len(ok)} completed sweeps but 0 readable per-mode pages in "
            f"{PAGES} — refusing to write 'no sweep has finished yet' over the "
            f"existing page. Check that the pages exist and are readable "
            f"(test -r), then re-run. First few expected but absent: {missing}")

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
    # QUALIFY MEANS THE GATE, not one half of it. This counted modes under the
    # 0.35 nm RMSD bar -- the pose-stability term -- and printed it in the header
    # as "N qualify" beside a rail ranked on the full gate: 62 qualify against 0
    # that actually clear it. A header that contradicts the order beneath it is
    # worse than no header.
    n_pri = int((ok.gate_tier == 0).sum())
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
{gs.nav("sweep.html", gs.step_counts())}
<main>
 <div id="rail">
  {rs.SEARCH_HTML}
  <div class="legend">ranked by the <b>100&nbsp;ns gate</b> &mdash; a mode is elevated
   when the warhead is within {hi_a:.1f}&nbsp;&#8491; of Cys113 SG for
   <b>&ge;{occ_pct:.0f}%</b> of the {int(_SWEEP_NS)}&nbsp;ns run <b>and</b> the pose held
   (max&nbsp;RMSD&nbsp;&lt;&nbsp;{rmx_a:.1f}&nbsp;&#8491; <i>or</i> mean&nbsp;&lt;&nbsp;{rmn_a:.1f}&nbsp;&#8491;,
   the allowance for brief spikes). Engagement is the headline because it leads
   the gate; stability sits beside it. <b>{n_pass}</b> of {n_ok} clear it,
   {n_near} held but under the engagement bar, {n_left} left the site.
   <br><b>Engagement is measured over the first {occ_window/1000:.1f}&nbsp;ns of every
   run</b>, not over the whole of it. Runs are no longer the same length &mdash;
   a sweep now continues while the molecule is still in the site and stops when it
   leaves, capped at 10&nbsp;ns &mdash; and a fraction of a 1.2&nbsp;ns run is not the
   same quantity as a fraction of a 10&nbsp;ns one. The common window is identical
   to the full-run figure for the {n_fixed} fixed-length rows, so nothing collected
   earlier changed meaning. <b>{n_adaptive}</b> rows are adaptive; each says beside
   it how long it ran and whether it left.
   <br><span class="muted">The movie's first frame is the start of
   <b>production</b>, after 300&nbsp;ps of unrestrained equilibration &mdash; not the
   docked pose. Every mode here was selected under 3.0&nbsp;&#8491; docked and the
   median has drifted +1.5&nbsp;&#8491; by then, so a movie opening above 3&nbsp;&#8491;
   is the equilibration, not a mismatch.</span> {note}</div>
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
    print(f"\n  {len(tabs)} sweep reports, ranked by the 100 ns gate "
          f"({int((ok.gate_tier == 0).sum())} clear it) -> {REPORTS / 'sweep.html'}")


if __name__ == "__main__":
    main()
