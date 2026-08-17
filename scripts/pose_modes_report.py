#!/usr/bin/env python3
"""
Purpose: one molecule's docked poses, split into modes, medoids shown in 3D.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-11
Input: --candidate <ident> (needs rank_v2 rows and an exported pose SDF)
Output: 00_outputs/blacksmith/pose_modes/pose_modes_<N>.html

A mode IS a candidate row. This shows what that means for one molecule: how the
docked poses split, how many of each mode reach near-attack geometry, and where
each mode's medoid actually sits in the pocket.

WHAT IT CANNOT SHOW, AND SAYS SO ON THE PAGE. The full pose cloud is not on disk.
`export_nac_poses` writes one representative per mode, not every pose, so this
draws 4 medoids and not 500 dots. Persisting the whole cloud is #44, and it has
to come from the SAME run that produced the scores -- docking is stochastic with
no fixed seed, so a re-dock gives a different cloud and putting it beside these
numbers would show structures the scores were not computed from.

THE POSES AND THE RECEPTOR ARE IN ONE FRAME. The SDF holds poses docked into the
chemist-prepared 3IKD (D0059), and the receptor drawn here is that same file. No
superposition happens in this script, because a second implementation of "fit the
two frames" is how two answers to one question appear.
"""

from __future__ import annotations

import argparse
import base64
import glob
import json
import logging
import sys
from datetime import date
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import gui_shell as gs                       # noqa: E402
from shared import mode_assets as massets                # noqa: E402
from shared import outputs as sout                       # noqa: E402
from shared import report_theme as rt                    # noqa: E402

log = logging.getLogger("pose-modes")
B = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith")
RECEPTOR = Path("/data/lab_vm/modifiable/inhibition/receptor_3ikd_prep/3IKD_noligand.pdb")
#: The production run's poses, from `run.topic` -- one definition, in mode_assets.
POSES = massets.poses_dir()

#: Cys113 in the crystal numbering the receptor file uses.
CYS_RESI = 113
#: One colour per mode, extended from the pipeline schematic's palette so the
#: two figures agree about which mode is which.
MODE_COLS = ["#0072ce", "#7b5ea7", "#c2703d", "#0f7a54", "#b3261e", "#8a6d1f"]


def mode_rows(ident: str) -> pd.DataFrame:
    """Every mode of this molecule, from the newest rank table that holds it."""
    out = []
    for tier, score in (("T4", "conditional_eb"), ("T3", "enrichment_conditional")):
        fs = sorted(glob.glob(str(B / f"rank_v2/rank_v2_{tier}_{score}_*.csv")))
        if not fs:
            continue
        d = pd.read_csv(fs[-1])
        s = d[d.parent_ident == ident]
        if len(s):
            out.append(s.sort_values("mode"))
    if not out:
        return pd.DataFrame()
    return out[0]


def sweep_rows(ident: str) -> pd.DataFrame:
    """The triage sweep, per mode. Scored PER MODE, so it is joined per mode."""
    fs = sorted(glob.glob(str(B / "attack_sweep/attack_sweep_*.csv")))
    if not fs:
        return pd.DataFrame()
    d = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    key = "parent_ident" if "parent_ident" in d.columns else "ident"
    return d[d[key] == ident].drop_duplicates("ident", keep="last")


def medoids(ident: str):
    """(mol, mode, pose_rank, energy_rank) per exported pose, in file order."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    p = POSES / f"{ident}.sdf"
    if not p.is_file():
        raise SystemExit(f"no exported pose set at {p}")
    out = []
    for m in Chem.SDMolSupplier(str(p), removeHs=False):
        if m is None:
            continue
        g = lambda k: (m.GetProp(k) if m.HasProp(k) else None)   # noqa: E731
        out.append((m, g("mode"), g("pose_rank"), g("energy_rank")))
    return out


def pose_pdb(mol) -> str:
    """One pose as a PDB block, residue MOL, so 3Dmol can select it by resn."""
    from rdkit import Chem
    m = Chem.Mol(mol)
    for a in m.GetAtoms():
        ri = a.GetPDBResidueInfo()
        if ri is None:
            ri = Chem.AtomPDBResidueInfo(f" {a.GetSymbol():<3s}")
            a.SetMonomerInfo(ri)
        ri.SetResidueName("MOL")
        ri.SetResidueNumber(1)
        ri.SetIsHeteroAtom(True)
    return "\n".join(l for l in Chem.MolToPDBBlock(m).splitlines()
                     if l.startswith(("ATOM", "HETATM")))


def viewer(ident: str, poses, receptor: str, ran_mode: str | None) -> str:
    """All medoids in the pocket, one colour per mode, each toggleable.

    Lazy boot and an explicit render(): 3Dmol draws nothing until render() is
    called, and a viewer built into a closed <details> has a zero-height
    container. Both cost a blank panel and neither announces itself.
    """
    blocks, ctl = [], []
    for i, (mol, mode, prank, erank) in enumerate(poses):
        col = MODE_COLS[i % len(MODE_COLS)]
        blocks.append({"m": mode, "col": col})
        tag = " · ran 100 ns" if mode == ran_mode else ""
        # NOT `energy_rank`. backfill_pose_rank.py stamps energy_rank = file
        # position, so it reads 1,2,3,4 for any molecule and means nothing.
        # Showing it invited the reading "mode 0 holds the lowest-energy pose",
        # which is the opposite of how the representative is chosen.
        ctl.append(
            f'<label class="mk"><input type="checkbox" id="pm_{ident}_{i}" checked>'
            f'<i style="background:{col}"></i>mode {mode}'
            f'<em>medoid{tag}</em></label>')
    pdbs = "".join(
        f'<script type="text/plain" id="pm_{ident}_pdb{i}">{pose_pdb(m)}</script>'
        for i, (m, *_r) in enumerate(poses))
    return f"""
<div class="glwrap"><div class="glbox">
<div id="pm_{ident}_gl" style="position:absolute;inset:0"></div></div>
<div class="vctl"><label class="mk"><input type="checkbox"
 id="pm_{ident}_surf" checked><i style="background:#b9c7db"></i>pocket surface</label>
{''.join(ctl)}</div></div>
<script type="text/plain" id="pm_{ident}_rec">{receptor}</script>
{pdbs}
<script>
(function(){{
  let M = null;
  const N = {len(poses)}, COLS = {json.dumps([b["col"] for b in blocks])};
  let built = false, v = null, surf = null;
  function styleAll() {{
    v.setStyle({{}}, {{cartoon:{{color:'#c3ccd8', opacity:0.5}}}});
    const CC = Object.assign({{}}, (M.elementColors||{{}}).defaultColors||{{}},
                             {{C: 0xb3261e}});
    v.setStyle({{resi:[{CYS_RESI}]}},
               {{stick:{{radius:0.26, colorscheme:{{prop:'elem', map: CC}}}},
                cartoon:{{color:'#c3ccd8', opacity:0.5}}}});
    for (let i = 0; i < N; i++) {{
      const on = document.getElementById('pm_{ident}_' + i).checked;
      v.setStyle({{model: i + 1}},
                 on ? {{stick:{{radius:0.19, color: COLS[i]}}}} : {{}});
    }}
    if (surf) {{ try {{ v.removeSurface(surf.surfid); }} catch(e) {{}} surf = null; }}
    if (document.getElementById('pm_{ident}_surf').checked) {{
      surf = v.addSurface(M.SurfaceType.VDW, {{opacity:0.55, color:'#b9c7db'}},
                          {{model: 0, not: {{resi:[{CYS_RESI}]}}}});
    }}
    v.render();
  }}
  function boot() {{
    if (built) return; built = true;
    M = window.$3Dmol || window['3Dmol'];
    if (!M) {{ built = false; return; }}   // library not parsed yet; retry on load
    requestAnimationFrame(function(){{ requestAnimationFrame(function(){{
      v = M.createViewer(document.getElementById('pm_{ident}_gl'),
                         {{backgroundColor:'#eef1f6'}});
      v.addModel(document.getElementById('pm_{ident}_rec').textContent, 'pdb');
      for (let i = 0; i < N; i++)
        v.addModel(document.getElementById('pm_{ident}_pdb'+i).textContent, 'pdb');
      styleAll();
      v.zoomTo({{resn:'MOL'}}); v.zoom(0.5); v.resize(); v.render();
      for (let i = 0; i < N; i++)
        document.getElementById('pm_{ident}_'+i)
                .addEventListener('change', styleAll);
      document.getElementById('pm_{ident}_surf')
              .addEventListener('change', styleAll);
    }}); }});
  }}
  const host = document.getElementById('pm_{ident}_gl');
  const det = host && host.closest ? host.closest('details') : null;
  if (det) {{ if (det.open) window.addEventListener('load', boot);
              det.addEventListener('toggle', function(){{ if (det.open) boot(); }}); }}
  else {{ window.addEventListener('load', boot); }}
}})();
</script>"""


def per_pose(ident: str) -> pd.DataFrame:
    """Every docked pose's energy and geometry, from the run that scored them.

    The COORDINATES are gone -- nac_screen_v2 docks in a tempfile.mkdtemp and
    rmtree's it in a finally, so only each mode's representative survives as
    structure (#44). The per-pose measurements were written to nac_v3/poses_*
    and did survive, which is why energies can be shown when poses cannot.
    """
    for f in sorted(glob.glob(str(B / "nac_v3/poses_*.csv"))):
        try:
            d = pd.read_csv(f)
        except Exception:                                  # noqa: BLE001
            continue
        if "ident" in d.columns and (d.ident == ident).any():
            return d[d.ident == ident]
    return pd.DataFrame()


def energy_table(pp: pd.DataFrame, n: int = 10) -> str:
    """The n lowest-energy poses of the whole cloud, with their mode."""
    if pp.empty:
        return "<p class='na'>no per-pose record survives for this molecule.</p>"
    t = pp.nsmallest(n, "energy")
    head = ("<tr><th>energy rank</th><th>pose</th><th>mode</th><th>energy</th>"
            "<th>C&ndash;S&gamma;</th><th>angle</th><th>in range</th>"
            "<th>NAC viable</th></tr>")
    body = []
    for _, r in t.iterrows():
        i = int(r["mode"])
        col = MODE_COLS[i % len(MODE_COLS)] if i >= 0 else "#8a94a6"
        yes = "<strong>yes</strong>" if bool(r.viable) else "<span class='na'>no</span>"
        body.append(
            f"<tr><td>{int(r.energy_rank)}</td><td>#{int(r.pose_idx)}</td>"
            f"<td><i class='sw' style='background:{col}'></i>m{i}</td>"
            f"<td>{r.energy:.2f}</td><td>{r.distance:.2f} &Aring;</td>"
            f"<td>{r.angle:.1f}&deg;</td>"
            f"<td>{'yes' if bool(r.in_range) else '<span class=na>no</span>'}</td>"
            f"<td>{yes}</td></tr>")
    return f'<table class="modes"><thead>{head}</thead><tbody>{"".join(body)}</tbody></table>'


def table(mr: pd.DataFrame, sw: pd.DataFrame, poses) -> str:
    sw_by_mode = {}
    if len(sw) and "mode" in sw.columns:
        for _, r in sw.iterrows():
            sw_by_mode[str(int(r["mode"]))] = r
    head = ("<tr><th>mode</th><th>poses</th><th>consensus</th><th>in NAC range</th>"
            "<th>viable</th><th>viable fraction</th><th>enrichment</th>"
            "<th>spread</th><th>coherence</th><th>mean energy</th>"
            f"<th>gnina affinity</th><th>{gs.sweep_label()} sweep</th></tr>")
    body = []
    for i, (_m, mode, _pr, _er) in enumerate(poses):
        row = mr[mr["mode"].astype("Int64").astype(str) == str(mode)]
        if row.empty:
            continue
        r = row.iloc[0]
        s = sw_by_mode.get(str(mode))
        swept = (f"{float(s['frac_attack_ready']):.3f} attack-ready"
                 f" · {int(s['n_visits'])} visits" if s is not None
                 else "<span class='na'>not swept</span>")
        body.append(
            f"<tr><td><i class='sw' style='background:{MODE_COLS[i % len(MODE_COLS)]}'></i>"
            f"m{mode}</td>"
            f"<td>{int(r.n_poses_mode)}</td><td>{r.consensus:.3f}</td>"
            f"<td>{int(r.n_in_range)}</td><td>{int(r.n_viable)}</td>"
            f"<td>{r.viable_fraction*100:.1f}%</td><td>{r.enrichment:.2f}</td>"
            f"<td>{r.spread_a:.2f} &Aring;</td><td>{r.dir_coherence:.3f}</td>"
            f"<td>{r.mean_energy:.3f}</td><td>{r.Affinity:.2f}</td>"
            f"<td>{swept}</td></tr>")
    return f'<table class="modes"><thead>{head}</thead><tbody>{"".join(body)}</tbody></table>'


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--candidate", required=True)
    ap.add_argument("--name", default=None)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    ident = args.candidate

    mr, sw, poses = mode_rows(ident), sweep_rows(ident), medoids(ident)
    pp = per_pose(ident)
    if mr.empty:
        raise SystemExit(f"no rank_v2 rows for {ident}")
    rec = "\n".join(l for l in RECEPTOR.read_text().splitlines()
                    if l.startswith(("ATOM", "HETATM")))

    # WHICH MODE WAS ACTUALLY SIMULATED. Taken from the sweep row's own mode, not
    # assumed to be 0 -- the bare ident in the sweep table IS mode 0, and at
    # least one molecule in this project reached 100 ns on a minority mode.
    ran = None
    if len(sw) and "mode" in sw.columns:
        ran = str(int(sw.iloc[0]["mode"]))

    n_poses = int(mr.n_poses.max())
    in_modes = int(mr.n_poses_mode.sum())
    unassigned = n_poses - in_modes

    import importlib.util as _u
    _sp = _u.spec_from_file_location("mdprio_combine", REPO / "scripts" / "mdprio_combine.py")
    _mc = _u.module_from_spec(_sp); _sp.loader.exec_module(_mc)
    ver, code = _mc._version()
    three = (REPO / "scripts/.cache_3dmol-min.js").read_text()

    name = args.name or f"{ident} — poses and modes"
    title = f"{date.today().isoformat()} {name}"
    byline = " ".join(x for x in (f"version {ver}" if ver else "",
                                  f"“{code}”" if code else "") if x)

    page = f"""<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title}</title><style>{rt.CSS}
body{{max-width:1180px;margin:0 auto;padding:0 30px 70px}}
body>*{{max-width:none;padding-left:0;padding-right:0}}
.glwrap{{margin:1.1rem 0}}
.glbox{{position:relative;width:100%;height:520px;overflow:hidden;
  background:#eef1f6;border:1px solid var(--rule);border-radius:5px}}
.glbox canvas{{position:absolute;top:0;left:0}}
.vctl{{display:flex;flex-wrap:wrap;gap:.5rem 1.5rem;padding:.7rem .2rem 0}}
label.mk{{display:flex;align-items:center;gap:.45rem;font:600 12px var(--sans);
  cursor:pointer;user-select:none}}
label.mk i{{width:22px;height:4px;border-radius:2px;display:inline-block}}
label.mk em{{font-style:normal;font-weight:400;color:var(--muted)}}
table.modes{{border-collapse:collapse;width:100%;font-size:13px;margin:1rem 0 1.4rem}}
table.modes th,table.modes td{{padding:.42rem .7rem;text-align:right;
  border-bottom:1px solid var(--rule);white-space:nowrap}}
table.modes th:first-child,table.modes td:first-child{{text-align:left}}
table.modes th:last-child,table.modes td:last-child{{text-align:left}}
table.modes thead th{{font:600 11px var(--sans);color:var(--muted);
  text-transform:uppercase;letter-spacing:.05em;border-bottom:2px solid var(--rule)}}
table.modes td{{font-family:var(--mono)}}
i.sw{{width:11px;height:11px;border-radius:2px;display:inline-block;
  margin-right:.5rem;vertical-align:-1px}}
.na{{color:var(--muted)}}
.caveat{{border-left:3px solid var(--rule);padding:.1rem 0 .1rem 1rem;
  color:var(--muted);font-size:13.5px;margin:1.2rem 0}}
</style>
<!-- 3Dmol is vendored in the HEAD, before any viewer script runs. Loaded at the
     end of the body it is not defined yet when a viewer's module-level
     `window.$3Dmol` lookup executes, and the viewer silently draws nothing. -->
<script>{three}</script>
</head><body>
<header class="mast"><h1>{title}</h1>
<p class="standfirst">Timothy Wu &middot; {byline}</p></header>

<p><strong>{n_poses} docked poses split into {len(mr)} modes.</strong>
{in_modes} land in a mode; <strong>{unassigned}
({unassigned/n_poses*100:.1f}%)</strong> fall below the clustering threshold and
appear in no row at all &mdash; they are in the denominator of
<code>consensus</code> and never in the numerator.</p>

{table(mr, sw, poses)}

<p>Each mode's <strong>medoid</strong> in the pocket, in the frame it was docked
into. Toggle a mode to show or hide it. Cys113 is in red sticks; the surface
covers the receptor and is left off Cys113.</p>

{viewer(ident, poses, rec, ran)}

<h2>The ten lowest-energy poses of the {n_poses}</h2>
<p>Docking energy over the whole cloud, from the run that produced the scores
above. <strong>Their coordinates no longer exist</strong> &mdash; the screen docks
in a temporary directory and deletes it, so only each mode's representative
survives as structure (#44). The measurements survived; the structures did not.</p>

{energy_table(pp)}

<div class="caveat"><strong>The pose that was simulated is NOT the lowest-energy
one, by design.</strong> <code>nac_screen_v2</code>: <em>"the pose most central to
its own mode, not its lowest-energy member"</em> &mdash; the medoid of the
best-anchored quartile. Argmax of a noisy score recovered the crystal pose 6.7% of
the time against 26.7% for a typical member, so taking the best-scoring pose is
the thing this pipeline deliberately stopped doing. The energy table above is
therefore context, not a competing ranking.</div>

<div class="caveat"><strong>These are medoids, not the pose cloud.</strong>
<code>export_nac_poses</code> writes one representative per mode, so this shows
{len(poses)} structures and not {n_poses}. Persisting every pose per mode is
issue #44, and it has to come from the <em>same</em> run that produced these
scores &mdash; docking is stochastic with no fixed seed, so a re-dock gives a
different cloud and putting it beside these numbers would show structures the
scores were not computed from.</div>

<div class="caveat"><strong>Only mode {ran if ran else '?'} was ever
simulated.</strong> The {gs.sweep_label()} sweep ran on it, and the {gs.production_label()} MD and its
replicates all started from that same pose. The other modes have no dynamics of
any kind behind them &mdash; their columns above are docking-derived only.</div>

</body></html>"""

    dest = sout.Topic("blacksmith", "pose_modes").write(f"pose_modes_{ident}", ".html")
    dest.write_text(page)
    stable = dest.parent / f"pose_modes_{ident}.html"
    print(f"\n  {len(poses)} modes of {ident} -> {dest}  "
          f"({len(page)/1048576:.1f} MB)")


if __name__ == "__main__":
    main()
