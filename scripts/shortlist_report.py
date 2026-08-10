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


#: Heavy-atom contact cutoff, Angstrom. A residue counts as contacting the ligand
#: in a frame if any of its heavy atoms is within this of any ligand heavy atom.
CONTACT_A = 4.5
#: Polar contact: N/O to N/O within this. Called "polar", not "hydrogen bond" --
#: the movie PDB carries no hydrogens, so donor geometry cannot be checked and
#: calling it an H-bond would claim more than was measured.
POLAR_A = 3.5


def contacts(movie_pdb: Path, max_frames: int = 40) -> tuple[list, int]:
    """Per-residue contact frequency across the trajectory.

    Over FRAMES, not one snapshot: a single frame says which residues happened to
    be near the ligand at one instant, which for a molecule that moves is close to
    arbitrary. Frequency over the run is the thing a chemist can act on.

    Returns (rows, n_frames) where each row is
    (resname, resid, fraction_of_frames, polar_fraction).
    """
    import numpy as np
    text = movie_pdb.read_text()
    models = [m for m in text.split("MODEL")[1:]] or [text]
    step = max(1, len(models) // max_frames)
    models = models[::step]

    seen: dict[tuple, list] = {}
    for mdl in models:
        lig, prot = [], []
        for ln in mdl.splitlines():
            if not ln.startswith(("ATOM", "HETATM")):
                continue
            el = ln[76:78].strip() or ln[12:16].strip()[:1]
            if el == "H":
                continue
            try:
                xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
            except ValueError:
                continue
            rn, ri = ln[17:20].strip(), ln[22:26].strip()
            (lig if rn == "MOL" else prot).append((xyz, rn, ri, el))
        if not lig or not prot:
            continue
        L = np.array([a[0] for a in lig])
        P = np.array([a[0] for a in prot])
        d = np.linalg.norm(P[:, None, :] - L[None, :, :], axis=2)
        near = d.min(axis=1) <= CONTACT_A
        lig_pol = np.array([a[3] in ("N", "O") for a in lig])
        hit_res, pol_res = set(), set()
        for i, ok in enumerate(near):
            if not ok:
                continue
            _, rn, ri, el = prot[i]
            hit_res.add((rn, ri))
            if el in ("N", "O") and lig_pol.any():
                if d[i][lig_pol].min() <= POLAR_A:
                    pol_res.add((rn, ri))
        for k in hit_res:
            seen.setdefault(k, [0, 0])
            seen[k][0] += 1
            if k in pol_res:
                seen[k][1] += 1
    n = len(models)
    rows = [(rn, ri, c / n, p / n) for (rn, ri), (c, p) in seen.items()]
    rows.sort(key=lambda r: -r[2])
    return rows, n


def interaction_map(smi: str, rows: list, n_frames: int, cys_resi: int = 63,
                    offset: int = 50) -> str:
    """Ligand in the middle, contacting residues around it by frequency.

    A contact SUMMARY, not a LigPlot: the residues are placed for legibility, not
    at their real positions, and no line claims to join a particular ligand atom
    to a particular residue atom. Saying so on the figure is the point — a diagram
    that looks like a LigPlot will be read as one.

    NUMBERED AS THE CRYSTAL, NOT AS THE MD SYSTEM. GROMACS renumbers from 1, so
    the catalytic cysteine is residue 63 in the trajectory and Cys113 in every
    paper and every PDB entry. This report leaves the project, and a chemist
    reading "Cys63" would either not recognise it or would look up the wrong
    residue — the same offset that once had us draw a glutamate and label it
    Cys113.
    """
    import math
    keep = [r for r in rows if r[2] >= 0.20][:12]
    if not keep:
        return ""
    W, H, CX, CY = 720, 470, 360, 235
    core = depiction(smi, 250, 175) if smi else ""
    parts = []
    for i, (rn, ri, frac, pol) in enumerate(keep):
        a = -math.pi / 2 + 2 * math.pi * i / len(keep)
        x, y = CX + 250 * math.cos(a), CY + 158 * math.sin(a)
        is_cys = (rn == "CYS" and str(ri) == str(cys_resi))
        try:
            shown = f"{rn}{int(ri) + offset}"
        except ValueError:
            shown = f"{rn}{ri}"
        col = "#b3261e" if is_cys else ("#0f7a54" if pol > 0.2 else "#4a6885")
        x0, y0 = CX + 120 * math.cos(a), CY + 84 * math.sin(a)
        dash = " stroke-dasharray='4 3'" if pol > 0.2 else ""
        parts.append(
            f"<line x1='{x0:.0f}' y1='{y0:.0f}' x2='{x:.0f}' y2='{y:.0f}' "
            f"stroke='{col}' stroke-width='{0.8 + 2.6 * frac:.1f}' opacity='.5'"
            f"{dash}/>"
            f"<circle cx='{x:.0f}' cy='{y:.0f}' r='25' fill='{col}' fill-opacity='.12' "
            f"stroke='{col}' stroke-width='1.2'/>"
            f"<text x='{x:.0f}' y='{y - 2:.0f}' class='rl' fill='{col}'>"
            f"{shown}</text>"
            f"<text x='{x:.0f}' y='{y + 11:.0f}' class='rf' fill='{col}'>"
            f"{frac*100:.0f}%</text>")
    return f"""<svg viewBox="0 0 {W} {H}" class="imap" role="img"
 aria-label="residues contacting the ligand, by fraction of frames">
<style>.rl{{font:600 11px ui-monospace,monospace;text-anchor:middle}}
.rf{{font:9.5px ui-monospace,monospace;text-anchor:middle;opacity:.85}}
.ik{{font:10px Helvetica,Arial,sans-serif;fill:#5b6b80}}</style>
{''.join(parts)}
<image href="{core}" x="{CX-125}" y="{CY-88}" width="250" height="175"/>
<text x="10" y="{H-24}" class="ik">line width = fraction of frames in contact
&#183; dashed + green = polar contact (N/O within {POLAR_A} &#8491;)
&#183; red = catalytic Cys113 &#183; crystal numbering</text>
<text x="10" y="{H-10}" class="ik">Contact = any heavy atom within {CONTACT_A} &#8491;,
over {n_frames} frames of the run. Residues are placed for legibility, not at
their real positions.</text>
</svg>"""


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

    # THE STRUCTURE THE CHEMIST CAN OPEN. First model of the fitted movie: the
    # protein and ligand as simulated, PBC-repaired and CA-fitted, one frame.
    # Offered as a download rather than a path, because the recipient has no
    # access to this filesystem.
    imap, pdb_href, pdb_bytes = "", "", 0
    if mpdb.is_file():
        raw = mpdb.read_text()
        first = raw.split("ENDMDL")[0]
        if not first.lstrip().startswith(("MODEL", "ATOM", "HETATM", "TITLE", "REMARK")):
            first = raw
        frame1 = first.replace("MODEL", "REMARK MODEL", 1).rstrip() + "\nEND\n"
        pdb_bytes = len(frame1.encode())
        pdb_href = ("data:chemical/x-pdb;base64,"
                    + base64.b64encode(frame1.encode()).decode())
        try:
            rows_c, nfr = contacts(mpdb)
            imap = interaction_map(smiles_of(ident) or "", rows_c, nfr)
        except Exception as exc:                          # noqa: BLE001
            log.warning("%s: interaction map unavailable: %s", ident, exc)

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
    dl = (f'<a class="dl" download="{ident}_md.pdb" href="{pdb_href}">'
          f'Download the MD structure (PDB, {pdb_bytes/1024:.0f} KB)</a>'
          if pdb_href else "")

    return f"""
<section class="mol">
  <h2>{ident}</h2>
  <div class="top">
    <div class="struct">{f'<img alt="" src="{svg}">' if svg else ''}</div>
    <div class="side">
      <label for="s_{ident}">SMILES</label>
      <textarea id="s_{ident}" readonly rows="3" onclick="this.select()">{smi}</textarea>
      <table class="kv">{facts}</table>
      {dl}
    </div>
  </div>
  <details class="panel"><summary>2D interaction map
    <span class="hint">residues contacting the ligand, by fraction of frames</span></summary>
    <div class="pbody">{imap}</div></details>
  <details class="panel"><summary>MD movie
    <span class="hint">100 ns, surface by charge, ligand in yellow</span></summary>
    <div class="pbody">{movie}</div></details>
  <details class="panel"><summary>RMSD plots
    <span class="hint">ligand RMSD, warhead&ndash;Cys113 distance, attack angle</span></summary>
    <div class="pbody">
      <img class="plots" src="data:image/png;base64,{img}" alt="RMSD, distance and angle traces">
    </div></details>
</section>"""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--candidates", nargs="+", required=True)
    ap.add_argument("--name", default="Shortlist")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    # Version from the CHANGELOG, via the GUI's own parser rather than a second
    # copy of the logic -- one source, so the report and the GUI cannot disagree
    # about which release produced the numbers.
    import importlib.util as _u
    _sp = _u.spec_from_file_location("mdprio_combine", REPO / "scripts" / "mdprio_combine.py")
    _mc = _u.module_from_spec(_sp); _sp.loader.exec_module(_mc)
    ver, code = _mc._version()

    er = _er()
    three = (REPO / "scripts/.cache_3dmol-min.js").read_text()
    cls = classes()
    blocks = [block(c, er, three, cls) for c in args.candidates]
    blocks = [b for b in blocks if b]
    if not blocks:
        raise SystemExit("nothing to report")

    title = f"{date.today().isoformat()} {args.name}"
    byline_ver = " ".join(x for x in (f"version {ver}" if ver else "",
                                      f"\u201c{code}\u201d" if code else "") if x)
    page = f"""<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title}</title><style>{rt.CSS}{mov.VIEWER_CSS}
/* CENTRED ON THE BODY, not on each child. rt.CSS centres direct children
   individually, which leaves a page of differently-sized blocks looking
   left-anchored; giving the body itself the measure centres the whole column and
   every child then fills it. */
body{{max-width:1180px;margin:0 auto;padding:0 30px 70px}}
body>*{{max-width:none;padding-left:0;padding-right:0}}
section.mol{{border-top:4px solid var(--rule);padding-top:1.6rem;margin-top:2rem}}
section.mol:first-of-type{{border-top:0;margin-top:1rem}}
header.mast{{border-bottom:4px solid var(--rule);padding-bottom:1rem}}
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
svg.imap{{width:100%;height:auto;background:var(--card);border:1px solid var(--rule);
  border-radius:5px}}
a.dl{{display:inline-block;margin-top:.7rem;font:600 12px var(--sans);
  color:var(--blue);text-decoration:none;border:1px solid var(--blue);
  border-radius:4px;padding:.35rem .7rem}}
a.dl:hover{{background:var(--blue-pale)}}
</style></head><body>
<header class="mast"><h1>{title}</h1>
<p class="standfirst">Timothy Wu &middot; {byline_ver}</p></header>
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
