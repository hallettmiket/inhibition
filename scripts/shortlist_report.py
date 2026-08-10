"""
Purpose: one self-contained HTML for a chemist — structure, SMILES, MD movie, plots.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-10
Input: --candidates <ident...> (each needs a finished 100 ns trajectory)
Output: 00_outputs/blacksmith/shortlist/shortlist_<N>.html

For sending outside the project. It carries no ranking, no gate verdict and no
pipeline commentary — a reader who does not work on this repo cannot check those
claims and does not need them to look at a trajectory. What it carries is the
molecule, its SMILES as selectable text, the 100 ns movie, the RMSD plots, and the
measured values behind them.

SELF-CONTAINED ON PURPOSE. The per-molecule reports are ~9.5 MB each because the
movie frames are base64 in the page; four of them inlined is large, and that is
the price of a file that opens by double-clicking with no server, no directory and
no network. 3Dmol.js is vendored ONCE for the whole document rather than per
viewer.
"""

from __future__ import annotations

import argparse
import base64
import glob
import logging
import re
import sys
from datetime import date
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import md_movie as mov                      # noqa: E402
from shared import outputs as sout                      # noqa: E402
from shared import report_theme as rt                   # noqa: E402

log = logging.getLogger("shortlist")
B = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith")
MD = Path("/data/lab_vm/modifiable/inhibition/md_residence_3ikd")
BOUND_NM = 1.2


def _er():
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "elevation_report", REPO / "scripts" / "elevation_report.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules["elevation_report"] = m
    spec.loader.exec_module(m)
    return m


def smiles_of(ident: str) -> str | None:
    for sub, stem in (("04_t4_combinatorial", "D4"), ("03_t3_reinvent", "D3")):
        fs = sorted(glob.glob(f"/data/lab_vm/append_only/inhibition/{sub}/{stem}_*.parquet"),
                    key=lambda q: int(q.rsplit("_", 1)[1].split(".")[0]))
        if not fs:
            continue
        fr = pd.read_parquet(fs[-1]).drop_duplicates("candidate_id").set_index("candidate_id")
        if ident in fr.index:
            return str(fr.loc[ident, "canonical_smiles"])
    p = B / f"pose_sidecars/{ident}.json"
    if p.is_file():
        import json
        return json.loads(p.read_text()).get("canonical_smiles")
    return None


def depiction(smi: str, w: int = 340, h: int = 210) -> str:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Draw, rdCoordGen
    RDLogger.DisableLog("rdApp.*")
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return ""
    rdCoordGen.AddCoords(m)
    d = Draw.rdMolDraw2D.MolDraw2DSVG(w, h)
    d.drawOptions().bondLineWidth = 1
    Draw.rdMolDraw2D.PrepareAndDrawMolecule(d, m)
    d.FinishDrawing()
    svg = re.sub(r"<\?xml.*?\?>", "", d.GetDrawingText(), flags=re.S)
    svg = re.sub(r"<!--.*?-->", "", svg, flags=re.S)
    return "data:image/svg+xml;base64," + base64.b64encode(svg.encode()).decode()


def md_row(ident: str):
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


def sweep_row(ident: str):
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


def classes() -> dict:
    out = {}
    for tier, score in (("T4", "conditional_eb"), ("T3", "enrichment_conditional")):
        fs = sorted(glob.glob(str(B / f"rank_v2/rank_v2_{tier}_{score}_*.csv")))
        if fs:
            d = pd.read_csv(fs[-1]).drop_duplicates("parent_ident")
            out.update(dict(zip(d.parent_ident, d.warhead_class)))
    return out


def block(ident: str, er, three: str, cls: dict) -> str:
    rep = MD / ident / "md" / "rep1"
    if not rep.is_dir():
        log.warning("%s: no trajectory at %s", ident, rep)
        return ""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "mdprio_report", REPO / "scripts" / "mdprio_report.py")
    mp = importlib.util.module_from_spec(spec)
    sys.modules["mdprio_report"] = mp
    spec.loader.exec_module(mp)

    total_ns = mp.prod_ns(rep)
    s = mp.series(rep, er, total_ns)
    res = mp.residence(s)
    if res["status"] != "ok":
        log.warning("%s: %s", ident, res["status"])
        return ""

    mpdb = rep / "movie.pdb"
    if not mpdb.is_file():
        mov.build_movie_pdb(rep, mpdb, total_ps=total_ns * 1000.0)
    movie = ""
    nacs = None
    if mpdb.is_file():
        pdb_txt, dsg, labels, lpos = er.surface_payload(mpdb)
        movie = mov.viewer_html(pdb_txt, dsg, labels, lpos, "", elem_id=f"gl_{ident}")
        nacs = mp.nac_series(ident, rep, mpdb, total_ns)
    img = mp.figure(ident, s, res, er, nacs)

    smi = smiles_of(ident) or ""
    svg = depiction(smi) if smi else ""
    m, sw = md_row(ident), sweep_row(ident)

    rows = [("trajectory", f"{res['length_ns']:.1f} ns, {res['n_frames']:,} frames"),
            ("warhead class", cls.get(ident, "unclassified")),
            ("mean ligand RMSD", f"{res['rmsd_mean_nm']:.3f} nm"),
            ("max ligand RMSD", f"{res['rmsd_max_nm']:.3f} nm"),
            ("final ligand RMSD", f"{res['rmsd_final_nm']:.3f} nm"),
            ("residence fraction", f"{res['residence_frac']:.3f}")]
    if m is not None and pd.notna(m.get("explicit_frac_frames_engaged")):
        rows.insert(2, ("target engagement, 100 ns",
                        f"{float(m['explicit_frac_frames_engaged'])*100:.2f}%"))
    if res.get("left_at_ns") is not None:
        rows.append(("left the pocket at", f"{res['left_at_ns']:.1f} ns"))
    if sw is not None:
        rows.append(("attack-ready, 10 ns sweep", f"{float(sw.frac_attack_ready):.4f}"))
        rows.append(("median C&ndash;S&gamma; distance",
                     f"{float(sw.median_dist_a):.2f} &Aring;"))
        rows.append(("median attack angle", f"{float(sw.median_angle_deg):.1f}&deg;"))
    facts = "".join(f"<tr><th>{k}</th><td>{v}</td></tr>" for k, v in rows)

    return f"""
<section class="mol">
  <h2>{ident}</h2>
  <div class="top">
    <div class="struct">{f'<img alt="" src="{svg}">' if svg else ''}</div>
    <div class="side">
      <label for="s_{ident}">SMILES</label>
      <textarea id="s_{ident}" readonly rows="3" onclick="this.select()">{smi}</textarea>
      <table class="kv">{facts}</table>
    </div>
  </div>
  <h3>100 ns MD</h3>
  {movie}
  <h3>Plots</h3>
  <img class="plots" src="data:image/png;base64,{img}" alt="RMSD, distance and angle traces">
</section>"""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--candidates", nargs="+", required=True)
    ap.add_argument("--name", default="Shortlist")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    er = _er()
    three = (REPO / "scripts/.cache_3dmol-min.js").read_text()
    cls = classes()
    blocks = [block(c, er, three, cls) for c in args.candidates]
    blocks = [b for b in blocks if b]
    if not blocks:
        raise SystemExit("nothing to report")

    title = f"{date.today().isoformat()} {args.name}"
    page = f"""<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title}</title><style>{rt.CSS}{mov.VIEWER_CSS}
body>*{{max-width:1180px}}
section.mol{{border-top:2px solid var(--rule);padding-top:1.4rem;margin-top:1.8rem}}
section.mol:first-of-type{{border-top:0}}
h2{{font-family:var(--mono);font-size:1.05rem;color:var(--navy);margin:0 0 .8rem}}
h3{{font-size:.66rem;font-family:var(--mono);letter-spacing:.14em;
  text-transform:uppercase;color:var(--blue);margin:1.4rem 0 .5rem}}
.top{{display:grid;grid-template-columns:360px 1fr;gap:22px;align-items:start}}
@media(max-width:820px){{.top{{grid-template-columns:1fr}}}}
.struct img{{width:100%;height:auto;background:#fff;border:1px solid var(--rule);
  border-radius:5px}}
.side label{{font:600 .6rem var(--mono);letter-spacing:.12em;text-transform:uppercase;
  color:var(--muted);display:block;margin-bottom:.3rem}}
textarea{{width:100%;font-family:var(--mono);font-size:12.5px;padding:.5rem;
  border:1px solid var(--rule);border-radius:4px;background:var(--card);
  color:var(--ink);resize:vertical;margin-bottom:.9rem}}
table.kv{{border-collapse:collapse;width:100%;font-size:13px}}
table.kv th{{text-align:left;font-weight:500;color:var(--muted);padding:3px 14px 3px 0;
  white-space:nowrap}}
table.kv td{{font-family:var(--mono);padding:3px 0}}
img.plots{{width:100%;height:auto;border:1px solid var(--rule);border-radius:5px;
  background:#fff}}
</style></head><body>
<header class="mast"><h1>{title}</h1>
<p class="standfirst">Timothy Wu</p></header>
<script>{three}</script>
{''.join(blocks)}
</body></html>"""

    dest = sout.Topic("blacksmith", "shortlist").write("shortlist", ".html")
    dest.write_text(page)
    side = B / "shortlist" / "shortlist.html"
    side.write_text(page)
    print(f"\n  {len(blocks)} molecules -> {side}  ({len(page)/1048576:.1f} MB)")
    print(f"  versioned: {dest}")


if __name__ == "__main__":
    main()
