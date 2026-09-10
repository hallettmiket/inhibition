#!/usr/bin/env python3
"""
Purpose: one standalone viewer page per pose, for the presentation GUI.
Author: Timothy Wu (with Claude Code)
Date: 2026-09-10
Input: presentation_data/<ident>.json + .movie.pdb + .rmsd.png, written by
       scripts/presentation_build.py
Output: presentation_gui/<ident>.html

@twu383, 2026-09-10: "do not carry over any of the design and numbers from
the last guis." This page's CSS, layout and copy are written fresh for this
spec -- structure, a stats table (ligand RMSD mean/max, warhead distance mean/
max from Cys113), the RMSD plot, and the movie. The only REUSED pieces are
low-level chemistry/rendering utilities (`elevation_report.surface_payload`,
`shared.md_movie.viewer_html`) -- infrastructure this session fixed a real
clock defect in today, not page design, and re-deriving a 3D trajectory
viewer from scratch under presentation deadline would risk reintroducing
exactly that bug.
"""

from __future__ import annotations

import base64
import html
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import md_movie as mov                        # noqa: E402
from shared import outputs as sout                         # noqa: E402
import elevation_report as er                               # noqa: E402
import mdprio_report as mp                                  # noqa: E402

DATA = sout.Topic("blacksmith", "presentation_data").dir
OUT = sout.Topic("artist", "presentation_gui").dir

CSS = """
:root{
  --ink:#1a1f2b; --sub:#5b6472; --line:#e4e7ec; --bg:#ffffff; --panel:#f7f8fa;
  --accent:#1d6fd6; --bad:#c0392b; --good:#1c8a5a;
}
@media (prefers-color-scheme: dark){
  :root:not([data-force-light]){
    --ink:#e8ebf1; --sub:#9aa4b2; --line:#2a2f3a; --bg:#12151c; --panel:#181c25;
    --accent:#5b9df0; --bad:#e0685c; --good:#3fc088;
  }
}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
     font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}
.wrap{max-width:980px;margin:0 auto;padding:20px 24px 60px}
h1{font-size:19px;margin:0 0 2px}
.fam{color:var(--sub);font-size:12.5px;margin:0 0 18px}
.top{display:grid;grid-template-columns:1fr 1fr;gap:18px;margin-bottom:18px}
.card{background:var(--panel);border:1px solid var(--line);border-radius:10px;
      padding:14px 16px}
.card h2{font-size:11px;text-transform:uppercase;letter-spacing:.05em;
         color:var(--sub);margin:0 0 10px;font-weight:600}
.struct{display:flex;align-items:center;justify-content:center;min-height:190px}
.struct svg{max-width:100%;height:auto}
table.stats{width:100%;border-collapse:collapse}
table.stats td{padding:5px 0;border-bottom:1px dashed var(--line);font-size:13px}
table.stats td:first-child{color:var(--sub)}
table.stats td:last-child{text-align:right;font-variant-numeric:tabular-nums;
                          font-weight:600}
.plot{background:var(--panel);border:1px solid var(--line);border-radius:10px;
      padding:10px;margin-bottom:18px;text-align:center}
.plot img{max-width:100%;border-radius:6px}
.plots{display:grid;grid-template-columns:1fr 1fr;gap:18px;margin-bottom:18px}
.plots .plot{margin-bottom:0}
.plot h3{text-align:left;font-size:11px;text-transform:uppercase;letter-spacing:.05em;
         color:var(--sub);margin:0 0 8px;font-weight:600}
.movie{background:var(--panel);border:1px solid var(--line);border-radius:10px;
       padding:14px;overflow:hidden}
.movie h2{font-size:11px;text-transform:uppercase;letter-spacing:.05em;
          color:var(--sub);margin:0 0 10px;font-weight:600}
"""


def build_one(ident: str) -> Path | None:
    jf = DATA / f"{ident}.json"
    if not jf.is_file():
        print(f"  SKIP {ident}: no data ({jf})")
        return None
    d = json.loads(jf.read_text())
    movie_pdb = DATA / f"{ident}.movie.pdb"
    plot_png = DATA / f"{ident}.rmsd.png"
    warhead_png = DATA / f"{ident}.warhead.png"
    if not movie_pdb.is_file():
        print(f"  SKIP {ident}: no movie file")
        return None

    rep = Path(d["rep_dir"])
    ra = mp.reactive_atom(d["parent"], rep)
    reactive_idx = ra["heavy_idx"] if ra else None
    pdb_txt, dist, labels, lpos = er.surface_payload(
        movie_pdb, reactive_idx=reactive_idx)
    three = (REPO / "scripts/.cache_3dmol-min.js").read_text()
    elem = "gl_" + "".join(c if c.isalnum() else "_" for c in ident)
    movie_block = mov.viewer_html(
        pdb_txt, dist, labels, lpos, three, total_ps=d.get("total_ps"),
        elem_id=elem)

    def _img(png: Path, alt: str) -> str:
        if not png.is_file():
            return "<p style='color:var(--sub)'>not available</p>"
        b64 = base64.b64encode(png.read_bytes()).decode()
        return f'<img src="data:image/png;base64,{b64}" alt="{alt}">'

    rmsd_plot_img = _img(plot_png, "ligand RMSD vs time")
    warhead_plot_img = _img(warhead_png, "warhead-Cys113 distance vs time")

    # `depiction()` returns a data: URI (base64 SVG), not raw markup -- an
    # <img> is the correct container; dropping the URI string straight into
    # the DOM rendered as inert text with no visible structure.
    struct_uri = d.get("structure_svg")
    svg = (f'<img src="{struct_uri}" alt="{html.escape(d.get("smiles") or "")}">'
          if struct_uri else "<p style='color:var(--sub)'>no structure</p>")

    def fmt(k, unit=" &Aring;", nd=2):
        v = d.get(k)
        return f"{float(v):.{nd}f}{unit}" if v is not None else "&mdash;"

    stats = f"""
      <tr><td>ligand RMSD, mean</td><td>{fmt('rmsd_mean_a')}</td></tr>
      <tr><td>ligand RMSD, max</td><td>{fmt('rmsd_max_a')}</td></tr>
      <tr><td>warhead&ndash;Cys113, mean</td><td>{fmt('warhead_dist_mean_a')}</td></tr>
      <tr><td>warhead&ndash;Cys113, max</td><td>{fmt('warhead_dist_max_a')}</td></tr>
      <tr><td>warhead&ndash;Cys113, min</td><td>{fmt('warhead_dist_min_a')}</td></tr>
      <tr><td>mechanism</td><td>{html.escape(str(d.get('mechanism') or '&mdash;'))}</td></tr>
      <tr><td>trajectory</td><td>{d.get('total_ps',0)/1000:.0f} ns, {d.get('n_frames','?')} frames</td></tr>
    """

    # @twu383, 2026-09-10: "the only subtitle on the viewer should be the mol
    # names ... The only annotation should be the thiadiazole-benzene family."
    # The molecule's own id, plain -- and the ONE exception, appended only for
    # that one family, never a second family name invented for any other pose.
    sub = html.escape(d.get("parent", ""))
    if d.get("family") == "thiadiazole_benzene":
        sub += " &middot; thiadiazole-benzene family"

    page = f"""<!doctype html><html><head><meta charset="utf-8">
<title>{html.escape(ident)}</title><style>{CSS}{mov.VIEWER_CSS}</style></head><body>
<div class="wrap">
  <h1>{html.escape(ident)}</h1>
  <p class="fam">{sub}</p>
  <div class="top">
    <div class="card"><h2>Structure</h2><div class="struct">{svg}</div></div>
    <div class="card"><h2>100 ns summary</h2><table class="stats">{stats}</table></div>
  </div>
  <div class="plots">
    <div class="plot"><h3>Ligand RMSD vs time</h3>{rmsd_plot_img}</div>
    <div class="plot"><h3>Warhead&ndash;Cys113 distance vs time</h3>{warhead_plot_img}</div>
  </div>
  <div class="movie"><h2>100 ns trajectory</h2>{movie_block}</div>
</div>
</body></html>"""

    OUT.mkdir(parents=True, exist_ok=True)
    dest = OUT / f"{ident}.html"
    dest.write_text(page)
    return dest


def main() -> None:
    import argparse
    import presentation_build as pb
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--only", nargs="*", default=None)
    args = ap.parse_args()

    idents = [e["ident"] for e in pb.MANIFEST
             if args.only is None or e["ident"] in args.only]
    ok = 0
    for ident in idents:
        try:
            dest = build_one(ident)
        except Exception as exc:                          # noqa: BLE001
            print(f"  FAILED {ident}: {exc}")
            dest = None
        if dest:
            ok += 1
            print(f"  OK {ident} -> {dest}")
    print(f"\n  {ok}/{len(idents)} pages built -> {OUT}")


if __name__ == "__main__":
    main()
