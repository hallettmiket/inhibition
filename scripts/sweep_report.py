#!/usr/bin/env python3
"""
Purpose: one swept MODE -> one report, the same shape as a 100 ns MD report.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-16
Input: --ident <parent>_m<k> (needs sweep_assets built for it)
Output: mdprio_reports/sweep_pages/<ident>.html

@tt8804: "just copy over the design for MD results and add/modify the tables to
use the sweep result data".

SO THIS IS `mdprio_report`'S PAGE, NOT A NEW ONE. Same `report_theme` CSS, same
masthead, same `<details class="panel">` blocks in the same order, same movie
viewer from `shared/md_movie`. What changes is the DATA: 10 ns instead of 100,
the sweep's own readings in the tables, and the pose that was simulated.

Three attempts at a bespoke sweep layout were three attempts too many. The design
question was already answered by the MD page; the only real work is which numbers
go in the tables.
"""

from __future__ import annotations

import argparse
import base64
import glob
import html
import logging
import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import md_movie as mov                        # noqa: E402
from shared import target_config as _tc                   # noqa: E402

#: Sweep length, derived -- these strings said 10 ns while the sweep has run
#: at 8 ns since D0085 (@tt8804: "update the gui to say 8 ns sweep not 10").
_SWEEP_NS = int(round(_tc.md_sweep_ps() / 1000))
from shared import report_theme as rt                     # noqa: E402
from shared import run_paths as rp                  # noqa: E402

log = logging.getLogger("sweep-report")
B = rp.BLACKSMITH
REPORTS = rp.reports_dir()
ASSETS = REPORTS / "sweep_assets"
PAGES = REPORTS / "sweep_pages"


def sweep_row(ident: str):
    """This mode's sweep result, newest attempt wins."""
    import pandas as pd
    best, best_t = None, -1.0
    for f in rp.sweep_result_files():                 # ordered, not raw glob
        try:
            d = pd.read_csv(f)
        except Exception:                                  # noqa: BLE001
            continue
        if "ident" not in d.columns:
            continue
        hit = d[d.ident.astype(str) == ident]
        if len(hit) and os.path.getmtime(f) > best_t:
            best, best_t = hit.iloc[0], os.path.getmtime(f)
    return best


def rank_row(ident: str):
    """What the docking ranking said about this mode, for the comparison table."""
    import pandas as pd
    # THE `*` WHERE THE TOPIC GOES MATCHED EVERY SCREEN. `rank_v2_T4_*_conditional_eb_*`
    # is satisfied by nac_v3, nac_v4 and nac_v5 alike, so this read whichever
    # sorted last -- a previous run's ranking on this run's report page.
    fs = glob.glob(str(rp.BLACKSMITH / "rank_v2" /
                       f"rank_v2_T4_{rp.topic()}_conditional_eb_*.csv"))
    if not fs:
        return None
    f = max(fs, key=lambda p: int(p.rsplit("_", 1)[1].split(".")[0]))
    d = pd.read_csv(f)
    hit = d[d.ident.astype(str) == ident]
    return hit.iloc[0] if len(hit) else None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--ident", required=True)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    import pandas as pd

    ident = args.ident
    s = sweep_row(ident)
    if s is None:
        raise SystemExit(f"{ident}: no sweep result")
    r = rank_row(ident)
    png, pdb = ASSETS / f"{ident}.png", ASSETS / f"{ident}.pdb"

    def num(v, fmt="{:.3f}"):
        try:
            if v is None or pd.isna(v):
                return "—"
            return fmt.format(float(v))
        except Exception:                                  # noqa: BLE001
            return str(v)

    ar = float(s.frac_attack_ready) if not pd.isna(s.get("frac_attack_ready")) else None
    vis = int(s.n_visits) if not pd.isna(s.get("n_visits")) else 0
    verdict = ("Reaches attack geometry" if (ar or 0) > 0.01 else "Never in position")

    # ---- the movie, from the asset the sweep build wrote --------------------
    movie_block = ""
    if pdb.is_file():
        try:
            # The SAME payload builder the 100 ns report uses, so the movie
            # is coloured and fitted identically -- it lives in
            # scripts/elevation_report, which is where mdprio_report gets it.
            import elevation_report as er
            import mdprio_report as mp
            import sweep_assets as sa
            # THE REACTIVE ATOM COMES FROM THE MOLECULE, not from an atom name.
            # `pose_rank` is on the sweep row, so the trajectory this movie was
            # built from is resolvable without guessing a sibling directory.
            _parent = ident.rsplit("_m", 1)[0]
            _pr = s.get("pose_rank")
            _rep = sa.rep_dir(_parent, None if pd.isna(_pr) else int(_pr))
            _ra = mp.reactive_atom(_parent, _rep) if _rep is not None else None
            if _ra is None:
                raise ValueError(
                    f"{ident}: cannot resolve the reactive atom; refusing to "
                    f"label an arbitrary atom as the warhead")
            pdb_txt, dsg, labels, lpos = er.surface_payload(
                pdb, reactive_idx=_ra["heavy_idx"])
            three = (REPO / "scripts/.cache_3dmol-min.js").read_text()
            movie_block = mov.viewer_html(pdb_txt, dsg, labels, lpos, three)
        except Exception as exc:                           # noqa: BLE001
            movie_block = rt.callout(
                "Movie unavailable",
                f"The {_SWEEP_NS} ns trajectory rendered no viewer: <code>{html.escape(str(exc))}</code>. "
                "The readings below are unaffected — they come from the trajectory "
                "directly.", "warn")

    img = base64.b64encode(png.read_bytes()).decode() if png.is_file() else ""

    # ---- the tables: SWEEP data in the MD report's shape --------------------
    sweep_tbl = (
        '<table class="kv"><tbody>'
        f'<tr><th>attack-ready</th><td>{num(ar)} '
        f'({(ar or 0)*100:.1f}% of 10 ns)</td></tr>'
        f'<tr><th>sustained visits</th><td>{vis}</td></tr>'
        f'<tr><th>fraction in the distance window</th><td>{num(s.get("frac_in_window"))}</td></tr>'
        f'<tr><th>median warhead–Sγ distance</th><td>{num(s.get("median_dist_a"), "{:.2f}")} Å</td></tr>'
        f'<tr><th>closest approach</th><td>{num(s.get("min_dist_a"), "{:.2f}")} Å</td></tr>'
        f'<tr><th>median attack angle</th><td>{num(s.get("median_angle_deg"), "{:.1f}")}°</td></tr>'
        f'<tr><th>attack-ready at frame 0</th><td>{s.get("start_attack_ready", "—")}</td></tr>'
        f'<tr><th>mechanism</th><td>{html.escape(str(s.get("mechanism", "—")))}</td></tr>'
        '</tbody></table>')

    sel_tbl = (
        '<table class="kv"><tbody>'
        f'<tr><th>enrichment (docked)</th><td>{num(r.get("enrichment"), "{:.2f}") if r is not None else "—"}</td></tr>'
        f'<tr><th>viable fraction</th><td>{num(r.get("viable_fraction")) if r is not None else "—"}</td></tr>'
        f'<tr><th>poses in this mode</th><td>{num(r.get("n_poses_mode"), "{:.0f}") if r is not None else "—"}</td></tr>'
        f'<tr><th>rank within warhead class</th><td>{num(r.get("class_rank"), "{:.0f}") if r is not None else "—"}</td></tr>'
        '</tbody></table>'
        '<p>These are the DOCKED numbers that selected this mode for simulation. '
        'Over the completed campaign they do not predict the readings above '
        '(Spearman +0.016, p = 0.84, n = 168), so they are shown as provenance — '
        'what chose this molecule — and not as evidence about it.</p>')

    # ---- the molecule, and the pose that was simulated ---------------------
    # @tt8804: "show the damn structure in the viewer". The report had numbers,
    # plots and a movie and never once said what the molecule IS. The depiction
    # and the pose are the assets the ranking view already writes, embedded here
    # so the page stays self-contained like every other report.
    parent = ident.rsplit("_m", 1)[0]
    mode_n = int(ident.rsplit("_m", 1)[1]) if "_m" in ident else -1
    # THE DEPICTION, NOT A SECOND 3D PANEL. @tt8804: "get rid of the pose I can
    # just see the pose in the movie, I want to see the structure". The movie
    # already shows the molecule in the pocket, from the trajectory this page is
    # about; a static pose beside it answered a question nothing had asked while
    # the one thing the page never showed was what the molecule IS.
    #
    # Drawn here rather than reusing the 96x64 rail thumbnail: that one exists to
    # be recognised at a glance in a list, and a chemist reading a report needs to
    # read the atoms.
    struct = ""
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import Draw, rdCoordGen
        RDLogger.DisableLog("rdApp.*")
        smi = None
        for sub, stem in (("04_t4_combinatorial", "D4"), ("03_t3_reinvent", "D3")):
            fs = [str(x) for x in rp.frames(stem)]
            if not fs:
                continue
            fr = pd.read_parquet(fs[-1])
            hit = fr[fr.candidate_id.astype(str) == parent]
            if len(hit):
                smi = str(hit.iloc[0].canonical_smiles); break
        m = Chem.MolFromSmiles(smi) if smi else None
        if m is not None:
            # rdCoordGen, not Compute2DCoords: the default layout puts visibly
            # wrong angles on substituted centres, and this is the figure a
            # chemist decides from.
            rdCoordGen.AddCoords(m)
            d2 = Draw.rdMolDraw2D.MolDraw2DSVG(560, 300)
            d2.drawOptions().bondLineWidth = 2
            d2.drawOptions().addStereoAnnotation = True
            Draw.rdMolDraw2D.PrepareAndDrawMolecule(d2, m)
            d2.FinishDrawing()
            b = base64.b64encode(d2.GetDrawingText().encode()).decode()
            struct = (f'<img src="data:image/svg+xml;base64,{b}" alt="structure" '
                      f'style="max-width:620px;width:100%;background:#fff;'
                      f'border:1px solid var(--rule);border-radius:4px;padding:10px">'
                      f'<p class="mono" style="font-size:11px;color:var(--muted);'
                      f'word-break:break-all">{html.escape(smi or "")}</p>')
    except Exception as exc:                               # noqa: BLE001
        struct = rt.callout("Structure unavailable",
                            f"<code>{html.escape(str(exc))}</code>", "warn")
    blocks = [
        rt.masthead(ident, f"{verdict} &middot; 10 ns attack-geometry sweep",
                    "sweep result",
                    [("attack-ready", f"{(ar or 0)*100:.1f}%"),
                     ("visits", str(vis)),
                     ("class", str(s.get("warhead_class", "—")))]),
        f'<p>{rt.pill("Held" if (ar or 0) > 0.01 else "Left")} '
        f'Attack-ready {(ar or 0)*100:.1f}% of the run over {vis} independent '
        f'approaches.</p>',
        # STRUCTURE FIRST. Everything below is about a molecule, and the page
        # used to never show which one.
        ('<details class="panel" open><summary>Structure'
         '<span class="hint">what the molecule is — the pose is in the movie'
         '</span></summary><div class="pbody">'
         f'{struct}</div></details>') if struct else "",
        (f'<details class="panel" open><summary>Trajectory plots'
         f'<span class="hint">ligand RMSD with its maximum, and warhead–Cys113 '
         f'distance against the attack window</span></summary><div class="pbody">'
         f'<img src="data:image/png;base64,{img}" alt="{_SWEEP_NS} ns trajectory"></div>'
         f'</details>') if img else "",
        # OPEN BY DEFAULT (@tt8804: "update the sweep results page to show the md
        # movies"). The movie was built for all 147 modes and embedded in every
        # page -- but inside a collapsed <details>, while Structure, plots and
        # readings all open. A panel nobody expands is a panel nobody knows is
        # there, which is indistinguishable from one that was never built.
        (f'<details class="panel" open><summary>{_SWEEP_NS} ns movie'
         f'<span class="hint">surface by charge, ligand in yellow, CA-fitted</span>'
         f'</summary><div class="pbody">{movie_block}</div></details>')
        if movie_block else "",
        '<details class="panel" open><summary>Sweep readings'
        f'<span class="hint">what the {_SWEEP_NS} ns run measured</span></summary>'
        f'<div class="pbody">{sweep_tbl}</div></details>',
        '<details class="panel"><summary>How it was selected'
        '<span class="hint">the docked numbers — provenance, not evidence</span>'
        f'</summary><div class="pbody">{sel_tbl}</div></details>',
    ]

    page = (f'<!doctype html><html><head><meta charset="utf-8">'
            f'<meta name="viewport" content="width=device-width,initial-scale=1">'
            f'<title>{html.escape(ident)} — 10 ns sweep</title>'
            f"<style>{rt.CSS}{mov.VIEWER_CSS}</style></head><body>\n"
            + "".join(blocks) + "</body></html>")
    PAGES.mkdir(parents=True, exist_ok=True)
    dest = PAGES / f"{ident}.html"
    dest.write_text(page)
    print(f"  {ident}: {verdict} -> {dest}")


if __name__ == "__main__":
    main()
