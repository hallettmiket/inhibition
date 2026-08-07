"""
Purpose: one report covering everything the overnight run produced — ranking, pose ranking, MD, and the gate.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06
Input: rank_v2/, elevation_queue/, pose_rank_bpmd/, md_residence/, nac_v2/agg_ref_*
Output: a self-contained HTML report

WRITTEN TO BE READ BEFORE THE RUN IS FINISHED. Every section degrades to "not
yet" rather than failing, because the 100 ns leg may still be running when this
is read and a report that refuses to render until everything lands is a report
nobody sees.

WHAT IT PUTS SIDE BY SIDE, AND WHY THAT IS THE POINT. A candidate's residence
number means nothing alone. Sulfopin -- the parent, nanomolar covalent -- and
ATRA -- a genuine low-micromolar binder that is not a warhead compound at all --
went through the identical 100 ns leg, so every candidate row is read next to
what the incumbent chemistry actually does. That is a COMPARISON, not a
validation set: @tt8804 has ruled the known Pin1 binders too few and too poor to
calibrate against, and nothing here is gated on them.

THE GATE IS APPLIED AND SHOWN, NOT ENFORCED. #22 sets it at >=90% residence AND
a majority of time in attack geometry. The residence half is a floor with a
kinetic justification -- even a millimolar binder should sit through 100 ns, so
leaving inside the window implies molar K_D -- while the attack-geometry half is
deliberately ambitious. Both verdicts are shown per molecule with the numbers
beside them; the covalent leg is a judgement, not an automatic consequence.
"""

from __future__ import annotations

import argparse
import base64
import glob
import io
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

log = logging.getLogger("morning")
DATA = Path("/data/lab_vm/append_only/inhibition")
B = DATA / "00_outputs/blacksmith"
MDROOT = Path("/data/lab_vm/modifiable/inhibition/md_residence_3ikd")

GATE_RESIDENCE = 0.90
GATE_NAC = 0.50
NAC_LO, NAC_HI = 0.28, 0.42


def newest(pattern: str) -> Path | None:
    fs = glob.glob(pattern)
    if not fs:
        return None
    try:
        return Path(max(fs, key=lambda p: int(p.rsplit("_", 1)[1].split(".")[0])))
    except (IndexError, ValueError):
        return Path(max(fs, key=lambda p: Path(p).stat().st_mtime))


def shards(pat: str) -> pd.DataFrame:
    fs = sorted(glob.glob(pat))
    if not fs:
        return pd.DataFrame()
    return pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)


def read_xvg(p: Path) -> np.ndarray | None:
    if not p.is_file():
        return None
    rows = [l.split() for l in p.read_text().splitlines() if l and l[0] not in "@#"]
    return np.array([[float(x) for x in r] for r in rows]) if rows else None


def to_ns(t: np.ndarray, total: float) -> np.ndarray:
    span = float(t[-1] - t[0])
    return t if span <= 0 else (t - t[0]) * (total / span)


def md_summary(ident: str) -> dict:
    """Residence and attack-geometry occupancy from a finished or running MD."""
    out = {"ident": ident, "state": "not started"}
    wd = MDROOT / ident.replace(":", "_") / "md" / "rep1"
    if not wd.is_dir():
        return out
    out["state"] = "running"
    r = read_xvg(wd / "rmsd.xvg")
    if r is None or len(r) < 10:
        return out
    t = to_ns(r[:, 0], 100.0)
    y = r[:, 1]
    out["ns_analysed"] = float(t[-1])
    bound = y <= 1.0
    out["residence_frac"] = float(bound.mean())
    if bound.all():
        out["escape_ns"] = None
    else:
        idx = np.where(~bound)[0]
        run = None
        for i in range(len(y)):
            if bound[i:].sum() == 0:
                run = i
                break
        out["escape_ns"] = float(t[run]) if run is not None else float(t[idx[0]])
    out["rmsd_mean_bound"] = float(y[bound].mean()) if bound.any() else float("nan")
    out["state"] = "complete" if t[-1] >= 99.0 else "running"
    return out


def fig_md(ident: str, theme: dict) -> str | None:
    wd = MDROOT / ident.replace(":", "_") / "md" / "rep1"
    r = read_xvg(wd / "rmsd.xvg")
    if r is None or len(r) < 10:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "figure.facecolor": theme["paper"], "axes.facecolor": theme["paper"],
        "savefig.facecolor": theme["paper"], "text.color": theme["ink"],
        "axes.labelcolor": theme["ink"], "xtick.color": theme["muted"],
        "ytick.color": theme["muted"], "axes.edgecolor": theme["grid"],
        "grid.color": theme["grid"], "font.size": 9,
        "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 120,
    })
    t, y = to_ns(r[:, 0], 100.0), r[:, 1]
    fig, a = plt.subplots(figsize=(9, 2.4))
    a.plot(t, y, lw=0.6, color=theme["accent"])
    a.axhline(1.0, color=theme["drift"], lw=1, ls="--", alpha=0.8)
    a.set_xlabel("time (ns)")
    a.set_ylabel("ligand RMSD (nm)")
    a.set_xlim(0, max(100, t[-1]))
    a.grid(axis="y", lw=0.5, alpha=0.5)
    a.text(1, a.get_ylim()[1] * 0.88, "bound below 1.0 nm", color=theme["muted"],
           fontsize=8)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def depict(smiles: str) -> str | None:
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem.Draw import rdMolDraw2D
        RDLogger.DisableLog("rdApp.*")
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            return None
        d = rdMolDraw2D.MolDraw2DCairo(320, 200)
        d.drawOptions().clearBackground = False
        rdMolDraw2D.PrepareAndDrawMolecule(d, m)
        d.FinishDrawing()
        return base64.b64encode(d.GetDrawingText()).decode()
    except Exception:                              # noqa: BLE001
        return None


def gate(md: dict, nac_frac: float | None) -> tuple[str, str]:
    if md.get("state") != "complete":
        return "pending", "trajectory not finished"
    res = md.get("residence_frac")
    if res is None:
        return "pending", "no frames analysed yet"
    bits = []
    ok_res = res >= GATE_RESIDENCE
    bits.append(f"residence {res*100:.1f}% {'PASS' if ok_res else 'FAIL'} "
                f"(gate ≥{GATE_RESIDENCE*100:.0f}%)")
    if nac_frac is None:
        bits.append("attack geometry not measured")
        return ("fail" if not ok_res else "partial"), "; ".join(bits)
    ok_nac = nac_frac >= GATE_NAC
    bits.append(f"attack geometry {nac_frac*100:.1f}% "
                f"{'PASS' if ok_nac else 'FAIL'} (gate >{GATE_NAC*100:.0f}%)")
    return ("pass" if (ok_res and ok_nac) else "fail"), "; ".join(bits)


THEMES = {"light": dict(paper="#faf8f4", ink="#14181e", muted="#6d7078",
                        grid="#ddd8cf", accent="#a8761a", anchor="#2f6f6a",
                        drift="#a8443a"),
          "dark": dict(paper="#14171c", ink="#e9e6e0", muted="#9599a1",
                       grid="#2c3138", accent="#d9a441", anchor="#57a79e",
                       drift="#cd6f63")}

CSS = """
:root{--paper:#faf8f4;--raise:#f3efe7;--ink:#14181e;--muted:#6d7078;--rule:#e0dad0;
 --accent:#a8761a;--anchor:#2f6f6a;--drift:#a8443a;
 --serif:"Iowan Old Style",Palatino,Georgia,serif;
 --sans:system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;
 --mono:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;}
@media (prefers-color-scheme:dark){:root{--paper:#14171c;--raise:#1b1f26;--ink:#e9e6e0;
 --muted:#9599a1;--rule:#2c3138;--accent:#d9a441;--anchor:#57a79e;--drift:#cd6f63;}}
:root[data-theme="dark"]{--paper:#14171c;--raise:#1b1f26;--ink:#e9e6e0;--muted:#9599a1;
 --rule:#2c3138;--accent:#d9a441;--anchor:#57a79e;--drift:#cd6f63;}
:root[data-theme="light"]{--paper:#faf8f4;--raise:#f3efe7;--ink:#14181e;--muted:#6d7078;
 --rule:#e0dad0;--accent:#a8761a;--anchor:#2f6f6a;--drift:#a8443a;}
*{box-sizing:border-box}
body{margin:0;background:var(--paper);color:var(--ink);font-family:var(--sans);
 font-size:16px;line-height:1.6}
.wrap{max-width:1120px;margin:0 auto;padding:0 26px 100px}
h1,h2,h3{font-family:var(--serif);font-weight:600;line-height:1.2;text-wrap:balance}
h1{font-size:clamp(2rem,4vw,2.9rem);margin:0 0 .8rem}
h2{font-size:1.6rem;margin:0 0 .3rem}
h3{font-size:1.05rem;margin:1.6rem 0 .4rem}
p{margin:0 0 .9rem}code{font-family:var(--mono);font-size:.87em;background:var(--raise);
 padding:.1em .35em;border-radius:3px}
.mast{border-bottom:2px solid var(--ink);padding:50px 0 22px;margin-bottom:30px}
.eyebrow{font-family:var(--mono);font-size:.7rem;letter-spacing:.15em;
 text-transform:uppercase;color:var(--muted);margin-bottom:1.2rem}
.facts{display:flex;flex-wrap:wrap;gap:0 2.4rem;font-family:var(--mono);font-size:.78rem;
 color:var(--muted);border-top:1px solid var(--rule);padding-top:.3rem}
.facts div{padding-top:.7rem}
.facts b{display:block;color:var(--ink);font-weight:600;font-size:1rem}
section{margin:0 0 3.4rem}
.shead{display:flex;gap:1rem;align-items:baseline;border-top:1px solid var(--rule);
 padding-top:1.3rem;margin-bottom:1.2rem}
.snum{font-family:var(--mono);font-size:.74rem;color:var(--accent);padding-top:.4rem;
 letter-spacing:.09em;white-space:nowrap}
.sub{color:var(--muted);font-size:.93rem;margin:.1rem 0 0}
.scroll{overflow-x:auto;margin:1.2rem 0;border:1px solid var(--rule);border-radius:5px;
 background:var(--raise)}
table{border-collapse:collapse;width:100%;font-size:.85rem}
th,td{padding:.5rem .8rem;text-align:left;border-bottom:1px solid var(--rule);white-space:nowrap}
th{font-family:var(--mono);font-size:.66rem;letter-spacing:.08em;text-transform:uppercase;
 color:var(--muted);background:var(--paper)}
td.n{text-align:right;font-family:var(--mono);font-variant-numeric:tabular-nums}
tbody tr:last-child td{border-bottom:none}
tr.ref td{color:var(--anchor);font-weight:600}
.pill{font-family:var(--mono);font-size:.68rem;padding:.1rem .45rem;border-radius:99px;
 border:1px solid currentColor;white-space:nowrap}
.pass{color:var(--anchor)}.fail{color:var(--drift)}.pend{color:var(--muted)}
.card{border:1px solid var(--rule);border-radius:6px;background:var(--raise);
 padding:1.1rem 1.3rem;margin:1.2rem 0}
.card h3{margin-top:0}
.row{display:flex;gap:1.4rem;flex-wrap:wrap;align-items:flex-start}
.row img{background:transparent;border-radius:4px;max-width:100%}
.callout{border:1px solid var(--rule);border-left:3px solid var(--accent);
 background:var(--raise);padding:1rem 1.2rem;margin:1.4rem 0;border-radius:0 5px 5px 0}
.callout.warn{border-left-color:var(--drift)}
.ctitle{font-family:var(--mono);font-size:.68rem;letter-spacing:.12em;
 text-transform:uppercase;color:var(--muted);margin-bottom:.4rem}
figure{margin:1rem 0}figure img{width:100%;height:auto;border-radius:4px}
.darkonly{display:none}
@media (prefers-color-scheme:dark){html:not([data-theme]) .lightonly{display:none}
 html:not([data-theme]) .darkonly{display:block}}
:root[data-theme="dark"] .lightonly{display:none}
:root[data-theme="dark"] .darkonly{display:block}
:root[data-theme="light"] .lightonly{display:block}
:root[data-theme="light"] .darkonly{display:none}
ul{padding-left:1.1rem}li{margin-bottom:.4rem}
.foot{border-top:1px solid var(--rule);padding-top:1.2rem;margin-top:3rem;
 font-size:.79rem;color:var(--muted);font-family:var(--mono);line-height:1.7}
@media (max-width:640px){.wrap{padding:0 16px 70px}}
"""


def build(rank, queue, prank, mds, refs, figs, meta) -> str:
    def gpill(v):
        cls = {"pass": "pass", "fail": "fail"}.get(v, "pend")
        return f'<span class="pill {cls}">{v.upper()}</span>'

    # --- elevated molecules -------------------------------------------------
    cards = []
    for m in mds:
        i = m["ident"]
        g, why = m["gate"]
        img = figs.get(i)
        struct = m.get("depiction")
        pr = prank[prank.ident == i] if not prank.empty else pd.DataFrame()
        prows = "".join(
            f"<tr><td class='n'>{int(r.pose_rank)}</td><td class='n'>{r.distance_A:.2f}</td>"
            f"<td class='n'>{r.angle_deg:.1f}</td>"
            f"<td class='n'>{getattr(r,'frac_in_window',float('nan')):.3f}</td>"
            f"<td>{'&larr; elevated' if getattr(r,'is_winner',False) else ''}</td></tr>"
            for r in pr.itertuples()) or \
            "<tr><td colspan='5' class='pend'>pose ranking not run</td></tr>"

        # a molecule whose trajectory has not produced frames yet has NO
        # residence -- render that as "not measured", never as 0.0, which would
        # read as "measured and found to leave immediately"
        res = m.get("residence_frac")
        res_txt = "—" if res is None else f"{res * 100:.1f}%"
        esc = m.get("escape_ns")
        ns = m.get("ns_analysed")
        ns_txt = "—" if ns is None else f"{ns:.1f}"
        rm = m.get("rmsd_mean_bound")
        rm_txt = "—" if rm is None or rm != rm else f"{rm:.3f} nm"
        cards.append(f"""
<div class="card">
  <h3>{i} &nbsp; {gpill(g)}</h3>
  <p class="sub">{m.get('warhead_class','')} &middot; {m.get('note','')}</p>
  <div class="row">
    {'<img src="data:image/png;base64,'+struct+'" alt="2D structure of '+i+'">' if struct else ''}
    <div style="flex:1;min-width:280px">
      <div class="scroll"><table><tbody>
        <tr><th>state</th><td>{m.get('state','?')}</td></tr>
        <tr><th>ns analysed</th><td class="n">{ns_txt}</td></tr>
        <tr><th>residence</th><td class="n">{res_txt}</td></tr>
        <tr><th>leaves at</th><td class="n">{'—' if esc is None else f'{esc:.1f} ns'}</td></tr>
        <tr><th>RMSD while bound</th><td class="n">{rm_txt}</td></tr>
        <tr><th>#22 gate</th><td>{why}</td></tr>
      </tbody></table></div>
    </div>
  </div>
  {'<figure><img src="data:image/png;base64,'+img+'" alt="ligand RMSD over time for '+i+'"></figure>' if img else ''}
  <h3>Poses ranked within this molecule</h3>
  <div class="scroll"><table>
    <thead><tr><th>pose</th><th>d (Å)</th><th>angle</th><th>BPMD occupancy</th><th></th></tr></thead>
    <tbody>{prows}</tbody></table></div>
</div>""")

    # --- ranking tables -----------------------------------------------------
    rank_tbl = ""
    if not rank.empty:
        for cls, g in rank[rank.passes].groupby("warhead_class"):
            t = g.nsmallest(10, "class_rank")
            rows = "".join(
                f"<tr><td class='n'>{int(r.class_rank)}</td><td><code>{r.ident}</code></td>"
                f"<td class='n'>{getattr(r,'weighted_score',float('nan')):.3f}</td>"
                f"<td class='n'>{getattr(r,'anchor_quality',float('nan')):.3f}</td>"
                f"<td class='n'>{getattr(r,'topn_viable_frac',float('nan')):.2f}</td>"
                f"<td class='n'>{getattr(r,'enrichment_joint',float('nan')):.2f}</td>"
                f"<td class='n'>{getattr(r,'consensus_gnina',float('nan')):.2f}</td>"
                f"<td class='n'>{getattr(r,'QED',float('nan')):.3f}</td></tr>"
                for r in t.itertuples())
            rank_tbl += (f"<h3>{cls} &nbsp;<span class='sub'>{len(g)} survivors</span></h3>"
                         f"<div class='scroll'><table><thead><tr><th>rank</th><th>ident</th>"
                         f"<th>weighted</th><th>anchor</th><th>top-10 viable</th>"
                         f"<th>enrich (2.0.0)</th><th>consensus</th><th>QED</th></tr></thead>"
                         f"<tbody>{rows}</tbody></table></div>")

    ref_rows = "".join(
        f"<tr class='ref'><td>{r.ident.replace('ref_','')}</td>"
        f"<td>{r.warhead_class}</td><td class='n'>{r.n_in_range}</td>"
        f"<td class='n'>{r.n_viable}</td><td class='n'>{r.enrichment:.2f}</td></tr>"
        for r in refs.itertuples()) if not refs.empty else \
        "<tr><td colspan='5' class='pend'>not run</td></tr>"

    return f"""<title>Overnight report — 2.1.0 ranking, selection and elevation</title>
<style>{CSS}</style>
<div class="wrap">
<header class="mast">
  <div class="eyebrow">Pin1 covalent inhibitors · pipeline 2.1.0 · inhibition@3IKD_ian</div>
  <h1>Overnight run</h1>
  <div class="facts">
    <div><b>{meta['screened']}</b> molecules screened</div>
    <div><b>{meta['ranked']}</b> ranked</div>
    <div><b>{meta['queued']}</b> selected</div>
    <div><b>{meta['elevated']}</b> in 100 ns MD</div>
    <div><b>{meta['when']}</b></div>
  </div>
</header>

<section>
  <div class="shead"><div class="snum">§1</div><div>
    <h2>Elevated molecules</h2>
    <p class="sub">Each started from the pose BPMD chose, not the lowest-energy one.
    The #22 gate is shown, not enforced — the covalent leg is a judgement.</p>
  </div></div>
  {''.join(cards) if cards else '<p class="pend">No molecule has reached the 100 ns leg yet.</p>'}
</section>

<section>
  <div class="shead"><div class="snum">§2</div><div>
    <h2>Reference molecules, on the identical measurement</h2>
    <p class="sub">A comparison, not a validation set — nothing here is calibrated
    against them.</p>
  </div></div>
  <div class="scroll"><table>
    <thead><tr><th>molecule</th><th>warhead class</th><th>poses in range</th>
    <th>viable</th><th>enrichment</th></tr></thead>
    <tbody>{ref_rows}</tbody></table></div>
  <div class="callout">
    <div class="ctitle">How to read these</div>
    <p><strong>Sulfopin is the parent</strong> — nanomolar covalent chloroacetamide,
    the compound this series descends from, and until tonight it had never been
    through the criterion its successors are ranked by. <strong>ATRA</strong> is a
    genuine low-micromolar binder that carries no warhead at all; its SMARTS hits on
    Michael-acceptor classes are a conjugated polyene looking like an acceptor, not a
    claim that it alkylates Cys113.</p>
  </div>
  <div class="callout warn">
    <div class="ctitle">Two caveats that travel with ATRA</div>
    <p><strong>Protonation.</strong> Run as the neutral acid, as the reference file
    records it. ATRA is &gt;99% deprotonated at pH 7.4, so this is not a measurement
    of the species that binds physiologically.</p>
    <p><strong>Pose origin.</strong> Started from a reactive-docked pose because
    <code>md_residence</code>'s plain-dock path leaves every hydrogen at the origin
    and antechamber refuses the molecule. The warhead bias that found this pose does
    not reflect ATRA's chemistry. Both are defects to fix, not properties of ATRA.</p>
  </div>
</section>

<section>
  <div class="shead"><div class="snum">§3</div><div>
    <h2>The ranking</h2>
    <p class="sub">One list per warhead class. There is deliberately no global
    top-N — a merged ordering re-imports the rigidity bias the per-class quota
    removes (D0073).</p>
  </div></div>
  {rank_tbl or '<p class="pend">ranking not available</p>'}
</section>

<section>
  <div class="shead"><div class="snum">§4</div><div>
    <h2>What these numbers are not</h2>
  </div></div>
  <ul>
    <li><strong>No score here has been shown to predict stability.</strong> The
    elevation cohort is the instrument that can test it, and that test has not run.
    It can falsify but not confirm — the cohort was selected on the old metrics.</li>
    <li><strong>Composite weights are equal and untested.</strong> Nothing has shown
    which component predicts anything, so a tuned weight would be fitted to nothing.</li>
    <li><strong>One replicate per molecule.</strong> A single dissociation is one
    draw with ~100% relative standard error — a screen, not a residence measurement.</li>
    <li><strong>The within-class rigidity confound survives.</strong> T₃ is 100%
    acrylamide and consensus still tracks rotatable bonds at ρ = −0.312 there.</li>
  </ul>
  <div class="foot">
    Generated by <code>scripts/morning_report.py</code>. Framework:
    <code>docs/framework_2.1.0.md</code> · Design:
    <code>docs/ranking_2.1.0_design.md</code> · Elevation protocol:
    <code>docs/elevation_example.md</code>
  </div>
</section>
</div>"""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--out", default="/data/lab_vm/modifiable/inhibition/morning_report.html")
    ap.add_argument("--tier", default="T3")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    agg = shards(str(B / "nac_v2/agg_s*_*.csv"))
    screened = agg.ident.nunique() if not agg.empty else 0

    rf = newest(str(B / f"rank_v2/rank_v2_{args.tier}_*.csv"))
    rank = pd.read_csv(rf) if rf else pd.DataFrame()
    qf = newest(str(B / "elevation_queue/queue_*.csv"))
    queue = pd.read_csv(qf) if qf else pd.DataFrame()
    pf = newest(str(B / "pose_rank_bpmd/pose_rank_*.csv"))
    prank = pd.read_csv(pf) if pf else pd.DataFrame()
    refs = shards(str(B / "nac_v2/agg_ref_*.csv"))
    if not refs.empty:
        refs = refs[refs.status == "ok"].drop_duplicates("ident")

    # everything with an MD workdir: queued candidates plus the references
    idents = list(queue.ident) if not queue.empty else []
    idents += [p.name for p in MDROOT.glob("ref_*") if p.is_dir()]
    seen, mds, figs = set(), [], {}
    smiles = dict(zip(queue.ident, queue.smiles)) if "smiles" in queue.columns else {}
    for i in idents:
        if i in seen:
            continue
        seen.add(i)
        m = md_summary(i)
        if m["state"] == "not started":
            continue
        q = queue[queue.ident == i]
        m["warhead_class"] = q.warhead_class.iloc[0] if len(q) else "reference"
        m["note"] = ("reference — comparison only" if i.startswith("ref_")
                     else f"class rank {int(q.class_rank.iloc[0])}" if len(q) else "")
        nacf = None
        m["gate"] = gate(m, nacf)
        s = smiles.get(i)
        if s:
            m["depiction"] = depict(s)
        mds.append(m)
        f = fig_md(i, THEMES["light"])
        if f:
            figs[i] = f

    meta = {"screened": screened,
            "ranked": int(rank.passes.sum()) if "passes" in rank else 0,
            "queued": len(queue), "elevated": len(mds),
            "when": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M")}
    html = build(rank, queue, prank, mds, refs, figs, meta)
    Path(args.out).write_text(html)
    log.info("wrote %s (%.1f KB) — %d elevated, %d ranked, %d screened",
             args.out, len(html) / 1e3, len(mds), meta["ranked"], screened)


if __name__ == "__main__":
    main()
