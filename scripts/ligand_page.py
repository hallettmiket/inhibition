#!/usr/bin/env python3
"""
Purpose: the per-LIGAND ranking page — one row per molecule, not per mode.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-28
Input: the run's rank_v2 engagement table (topic from config)
Output: <reports>/ligands.html

@tt8804: "can i see rank by mol in gui". The Ranking page is per MODE -- 327,167
rows for nac_v6 -- because a mode is the thing that gets simulated. But the thing
you order, buy and test is a MOLECULE, and nothing showed that.

BOTH AGGREGATIONS ARE SHOWN, SIDE BY SIDE, because they answer different
questions and picking one silently would answer the wrong one (D0098):

  best  -- can this ligand reach attack geometry AT ALL. Immune to how many
           groups it happens to have, and therefore to docking depth. But it is
           a maximum over a noisy score, and argmax selection is separately
           measured as the worst rule available at the pose level (6.7% crystal
           recovery against 33.3% for the medoid of the well-anchored quartile).
  mean  -- how well does it engage TYPICALLY. Penalises a molecule that can only
           reach attack geometry one way in twenty -- right if you believe the
           mode population, wrong if you do not, because `n_modes` grows with
           docking depth and never saturates (D0092, b = +0.69).

`n_modes` is therefore printed beside the mean on every row: a mean over a
denominator that moves with runtime is comparable only at fixed depth.

NOT A CLAIM THAT ANY OF THESE BIND. The ordering is reachability of attack
geometry, which is a precondition for covalent chemistry and not evidence of it.
`rank_validated` is False for this run as for every other.
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

from shared import engagement_rank as er           # noqa: E402
from shared import run_paths as rp                 # noqa: E402
from shared import target_config as tc             # noqa: E402

log = logging.getLogger("ligand-page")

_CSS = """
:root{--ink:#16202b;--muted:#5b6b7b;--rule:#dbe2e9;--bg:#f6f8fa;--card:#fff;--accent:#1f4e79}
@media(prefers-color-scheme:dark){:root{--ink:#eef3f8;--muted:#93a4b4;--rule:#2a3742;--bg:#101820;--card:#18222c;--accent:#7db2e0}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.5 system-ui,-apple-system,Segoe UI,sans-serif}
.wrap{max-width:1500px;margin:0 auto;padding:26px 22px 70px}
h1{font-size:22px;margin:0 0 6px}
.sub{color:var(--muted);font-size:13.5px;margin-bottom:18px;max-width:100ch}
.bar{display:flex;gap:10px;flex-wrap:wrap;align-items:center;margin-bottom:14px}
input,select{padding:7px 10px;border:1px solid var(--rule);border-radius:5px;background:var(--card);color:var(--ink);font-size:14px}
.tw{overflow-x:auto;border:1px solid var(--rule);border-radius:6px;background:var(--card)}
table{border-collapse:collapse;width:100%;font-size:14px}
th,td{padding:8px 12px;border-bottom:1px solid var(--rule);text-align:right;white-space:nowrap;font-variant-numeric:tabular-nums}
th:nth-child(2),td:nth-child(2),th:last-child,td:last-child{text-align:left}
thead th{position:sticky;top:0;background:var(--card);cursor:pointer;font-size:12px;
  text-transform:uppercase;letter-spacing:.06em;color:var(--muted);border-bottom:2px solid var(--rule)}
thead th:hover{color:var(--accent)}
tbody tr:hover{background:rgba(127,127,127,.08)}
code{font:12.5px ui-monospace,Menlo,monospace}
.smi{max-width:340px;overflow:hidden;text-overflow:ellipsis;display:inline-block;vertical-align:bottom;color:var(--muted)}
.note{margin-top:16px;color:var(--muted);font-size:13px;max-width:100ch}
"""

_JS = """
const tb=document.querySelector('tbody');
let dir={};
document.querySelectorAll('thead th').forEach((th,i)=>th.onclick=()=>{
  dir[i]=!dir[i];
  const rows=[...tb.rows];
  rows.sort((a,b)=>{
    const x=a.cells[i].dataset.v??a.cells[i].textContent, y=b.cells[i].dataset.v??b.cells[i].textContent;
    const nx=parseFloat(x), ny=parseFloat(y);
    const c=(!isNaN(nx)&&!isNaN(ny))?nx-ny:String(x).localeCompare(String(y));
    return dir[i]?c:-c;});
  rows.forEach(r=>tb.appendChild(r));});
const q=document.getElementById('q'), f=document.getElementById('f');
function filt(){const s=q.value.toLowerCase(),c=f.value;
  [...tb.rows].forEach(r=>{const hit=r.textContent.toLowerCase().includes(s)&&(!c||r.dataset.c===c);
    r.style.display=hit?'':'none';});}
q.oninput=filt; f.onchange=filt;
"""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--topic", default=None)
    ap.add_argument("--sort", default="best", choices=("best", "mean", "median"))
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    topic = a.topic or rp.topic()
    fs = sorted(glob.glob(str(rp.BLACKSMITH /
                              f"rank_v2/rank_v2_T4_{topic}_engagement_*.csv")),
                key=os.path.getmtime)
    if not fs:
        raise SystemExit(f"no engagement ranking for topic {topic!r}")
    d = pd.read_csv(fs[-1])
    d["family"] = d.warhead_class.map(tc.family_of())

    # one row per MOLECULE, both aggregations
    # KEYED ON parent_ident DIRECTLY. Renaming it to `ident` collided with the
    # mode-level `ident` already in the frame, and pandas groups on the ambiguous
    # name rather than raising anything a reader would understand.
    lig = er.rank_ligands(d, how=a.sort, ligand_key="parent_ident")
    meta = (d.groupby("parent_ident")
              .agg(warhead_class=("warhead_class", "first"),
                   family=("family", "first"),
                   smiles=("smiles", "first"),
                   QED=("QED", "first"),
                   best_mode=("engagement", "idxmax"))
              .reset_index())
    meta["best_mode"] = d.loc[meta.best_mode, "mode"].values
    t = lig.merge(meta, on="parent_ident", how="left").rename(
        columns={"parent_ident": "ident"})
    t = t.sort_values("ligand_engagement", ascending=False).reset_index(drop=True)
    t.insert(0, "rank", range(1, len(t) + 1))
    log.info("%s: %d molecules from %d modes", topic, len(t), len(d))

    cols = [("rank", "rank"), ("ident", "molecule"),
            ("best_mode_engagement", "best"), ("mean_mode_engagement", "mean"),
            ("n_modes", "modes"), ("best_mode", "best mode"),
            ("warhead_class", "warhead class"), ("QED", "QED"), ("smiles", "SMILES")]
    body = []
    for r in t.itertuples():
        fam = r.family if isinstance(r.family, str) else ""
        cells = []
        for k, _ in cols:
            v = getattr(r, k, "")
            if k == "smiles":
                cells.append(f'<td><span class="smi" title="{html.escape(str(v))}">'
                             f'{html.escape(str(v))}</span></td>')
            elif isinstance(v, float):
                cells.append(f'<td data-v="{v}">{v:.3f}</td>')
            elif k == "ident":
                cells.append(f'<td><code>{html.escape(str(v))}</code></td>')
            else:
                cells.append(f'<td data-v="{v}">{html.escape(str(v))}</td>')
        body.append(f'<tr data-c="{html.escape(fam)}">' + "".join(cells) + "</tr>")

    fams = sorted({x for x in t.family.dropna().unique()})
    opts = "".join(f'<option value="{html.escape(f)}">{html.escape(f)}</option>' for f in fams)
    head = "".join(f"<th>{html.escape(lbl)}</th>" for _, lbl in cols)
    page = f"""<!doctype html><html><head><meta charset="utf-8">
<title>{topic} — ligands by engagement</title><style>{_CSS}</style></head><body><div class="wrap">
<h1>{topic} — {len(t):,} ligands by warhead engagement</h1>
<p class="sub">One row per MOLECULE. The Ranking page is per mode ({len(d):,} of them);
this is the thing you would actually order. <b>best</b> is the ligand's strongest
group — can it reach attack geometry at all. <b>mean</b> is how well it engages
typically, and <b>modes</b> is the denominator behind that mean, which grows with
docking depth and never saturates, so a mean is comparable only at fixed depth.
Neither is evidence that anything binds: this orders reachability of attack
geometry, and <code>rank_validated</code> is False.</p>
<div class="bar"><input id="q" placeholder="filter — id, SMILES, warhead…" size="34">
<select id="f"><option value="">all families</option>{opts}</select>
<span class="sub" style="margin:0">click a header to sort</span></div>
<div class="tw"><table><thead><tr>{head}</tr></thead><tbody>
{chr(10).join(body)}
</tbody></table></div>
<p class="note">Sorted by <b>{a.sort}</b>. Engagement is
<code>anchor_quality</code> — distance and angle to Cys113 SG, each 0–1 and
multiplied, so a pose at perfect distance and hopeless angle scores near zero
rather than half. D0098 measured it against the MD outcome at rho = +0.652,
against −0.015 for the column this replaced.</p>
</div><script>{_JS}</script></body></html>"""

    dest = rp.reports_dir() / "ligands.html"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(page, encoding="utf-8")
    print(dest)
    print(f"  {len(t):,} ligands · top: {t.iloc[0].ident} "
          f"(best {t.iloc[0].best_mode_engagement:.3f}, "
          f"{int(t.iloc[0].n_modes)} modes)")


if __name__ == "__main__":
    main()
