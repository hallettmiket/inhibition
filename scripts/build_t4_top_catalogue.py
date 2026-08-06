"""
Purpose: a single, self-contained, downloadable HTML catalogue of the top T_4
         covalent candidates by near-attack enrichment -- structure and ALL TEN
         docking poses for every molecule, visible simultaneously -- for @tt8804.
Author: Artist (murmurent), with Claude Code
Date: 2026-08-06
Input: the latest 04_t4_combinatorial/D4_*.parquet frame, the exported NAC pose
       SDFs (00_outputs/blacksmith/nac_poses/<candidate_id>.sdf), the 3IKD
       receptor, shared/nac_criterion.py's own constants.
Output: one .html file. No CDN dependency -- 3Dmol.js is vendored inline so it
        renders offline.

WHY THE THRESHOLD IS RESOLVED HERE AND NOT HARDCODED AS A LIST OF IDS. The
brief changed twice while this was being built -- >=6.00 (5 candidates), then
>5.70 (6, adding an SNAr) -- because @tt8804 changed the threshold, not because
the frame changed. Filtering the LATEST frame by NAC_ENRICHMENT_MIN on every
run means a future re-run with a new frame or a further threshold change is one
constant edit, not a re-transcription of ids that can drift from what is
actually on disk.

WHY THIS PAGE CARRIES THREE WARNINGS A NORMAL RESULTS PAGE WOULD NOT.
Ordinary practice would be to show the top candidates and let the ranking speak
for itself. It cannot, here, for reasons specific to this exact list:

  1. Five of six are one warhead class (BDHI) with ZERO crystallographic
     Cys113 positives -- the class cannot be validated in this project's own
     framework (D0067's report).
  2. The sixth is a DIFFERENT warhead class (SNAr) that WAS tested against
     measured actives and inactives and FAILED (AUC 0.558/0.451, D0065/D0070)
     -- worse than "unvalidated", this one is measured to carry no signal.
  3. Every enrichment on this page is a 200-run value, and the one molecule
     independently re-measured at 2,000 runs collapsed from 6.74x to 0.60x --
     BELOW chance (docs/medchem_t4_72f5671e89cb.md Sec 0). No BDHI molecule had
     ever been run at convergence before that report. The entire gap that put
     these candidates at the top of 5,765 scored molecules may not survive more
     search, and it has not been checked for five of the six.

A page that showed six confident-looking molecule cards without these would
read as a shortlist. It is not one, and saying so is part of the deliverable,
not a caveat bolted on afterwards -- see rules/headline_first.md's spirit
applied to a report rather than a chat reply: the reader who scrolls no
further than the top must still get the true picture.

WHY 3Dmol.js IS FETCHED ONCE AND VENDORED INLINE RATHER THAN CDN-LINKED.
py3Dmol's own `_make_html()` always emits a
`loadScriptAsync('https://cdn.jsdelivr.net/...')` call -- there is no vendored
copy anywhere in this repo or its five envs (checked). A page that phones home
for its own rendering engine is not "self-contained" and will render blank for
anyone reading it offline or behind a firewall. The library (BSD-3, ~525 KB
minified) is downloaded once by this script and inlined as a <script> block;
each per-viewer snippet then has its CDN-loader boilerplate stripped so it
calls the already-loaded global $3Dmol directly instead of racing a network
request six times.

WHY THE RECEPTOR IS RESOLVED VIA `pose3d.receptor_for("nac_pose_path")` AND
NEVER HARDCODED. These poses were docked into 3IKD, not 6VAJ -- the two sit
48.6 A apart (pose3d.py's own docstring) -- and that exact mistake ("poses were
given without the pocket for some reason??") already happened once on this
project today. Asking the module that owns the mapping is the only way this
page cannot make it a second time.
"""

from __future__ import annotations

import re
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from integration.app import depict, pose3d          # noqa: E402
from shared import nac_criterion as nac              # noqa: E402
from shared import outputs                            # noqa: E402

DATA = Path("/data/lab_vm/append_only/inhibition")
T4_DIR = DATA / "04_t4_combinatorial"
POSES_DIR = DATA / "00_outputs" / "blacksmith" / "nac_poses"

# Set by @tt8804, twice (>=6.00 then >5.70) -- read the LATEST frame and filter
# live rather than pinning a list of ids that would silently go stale.
NAC_ENRICHMENT_MIN = 5.70
TOP_N = 12  # generous cap; the filter is expected to return ~6

THREEDMOL_CDN = "https://cdn.jsdelivr.net/npm/3dmol@2.5.5/build/3Dmol-min.js"
THREEDMOL_CACHE = REPO / "scripts" / ".cache_3dmol-min.js"


# ---------------------------------------------------------------------------
# frame resolution
# ---------------------------------------------------------------------------

def latest_frame() -> Path:
    """The newest D4_<N>.parquet, by integer suffix -- not by mtime or glob order."""
    fs = list(T4_DIR.glob("D4_*.parquet"))
    if not fs:
        raise FileNotFoundError(f"no D4_<N>.parquet under {T4_DIR}")
    return max(fs, key=lambda p: int(re.search(r"_(\d+)\.parquet$", p.name).group(1)))


def top_candidates(frame: Path) -> pd.DataFrame:
    df = pd.read_parquet(frame)
    sel = df[(df["nac_enrichment"] > NAC_ENRICHMENT_MIN)
             & df["nac_pose_path"].notna()].copy()
    sel = sel.sort_values("nac_enrichment", ascending=False).head(TOP_N)
    if sel.empty:
        raise ValueError(
            f"no candidates with nac_enrichment > {NAC_ENRICHMENT_MIN} and an "
            f"exported pose in {frame.name}")
    return sel


# ---------------------------------------------------------------------------
# per-pose properties, parsed straight from the SDF text
#
# NOT via pose3d.read_poses(): that parser only keeps FLOAT-parseable tags
# (pose3d.py's own docstring explains why -- gnina's blank header lines), which
# silently drops nac_viable ("True"/"False") and nac_angle_kind (a string). The
# per-pose table on this page needs exactly those, so it is parsed separately
# rather than asking pose3d to do something its contract does not cover.
# ---------------------------------------------------------------------------

_PROP_RE = re.compile(r">\s*<(\w+)>[^\n]*\n([^\n]*)\n")


def pose_properties(sdf_path: Path) -> list[dict]:
    text = sdf_path.read_text(errors="ignore")
    rows = []
    for rec in text.split("$$$$"):
        if not rec.strip():
            continue
        props = {m.group(1): m.group(2).strip() for m in _PROP_RE.finditer(rec)}
        if "pose_rank" in props:
            rows.append(props)
    rows.sort(key=lambda r: int(r["pose_rank"]))
    return rows


# ---------------------------------------------------------------------------
# 3Dmol.js vendoring
# ---------------------------------------------------------------------------

def vendored_3dmol_js() -> str:
    """The 3Dmol.js library text, fetched once and cached beside this script.

    Cached in `scripts/` (tiny relative to the repo, and not governed data) so
    a second run of this builder -- or a run with no network reachable -- does
    not re-fetch or fail. Delete the cache file to force a refresh.
    """
    if THREEDMOL_CACHE.is_file() and THREEDMOL_CACHE.stat().st_size > 100_000:
        return THREEDMOL_CACHE.read_text(errors="ignore")
    with urllib.request.urlopen(THREEDMOL_CDN, timeout=30) as resp:
        text = resp.read().decode("utf-8", errors="ignore")
    if len(text) < 100_000:
        raise RuntimeError(f"3Dmol.js fetch looks truncated: {len(text)} bytes")
    THREEDMOL_CACHE.write_text(text)
    return text


# Removes ONLY the `loadScriptAsync` function definition and the
# `if(typeof $3Dmolpromise === 'undefined') {...}` block that kicks off the CDN
# fetch. Deliberately narrow: an earlier version of this regex also swallowed
# the few lines just after it -- `var viewer_ID = null;` and the
# `warn.parentNode.removeChild(warn)` call that deletes the "3Dmol.js failed to
# load" placeholder <p> -- because it matched greedily up to the FIRST
# `$3Dmolpromise.then(function() {`, which comes after those lines, not before.
# The bug was invisible in a rendered screenshot (the canvas still drew, on
# top of the now-orphaned placeholder) and only showed up by checking the DOM:
# the failure warning was still sitting there, unremoved, describing a failure
# that had not happened.
_LOADER_BLOCK = re.compile(
    r"var loadScriptAsync = function\(uri\)\{.*?\nif\(typeof \$3Dmolpromise"
    r" === 'undefined'\) \{\n\$3Dmolpromise = null;\n  \$3Dmolpromise = "
    r"loadScriptAsync\('[^']*'\);\n\}\n\n",
    re.S)
# The `$3Dmolpromise.then(function() { ... });` wrapper around the per-viewer
# setup calls. Replaced, not merely stripped: the six calls inside it are
# rewritten into a named `function initViewer_UID(){...}` and deferred behind
# an IntersectionObserver (see `offline_viewer_html`) rather than left to run
# immediately in document order.
#
# WHY DEFERRED, MEASURED RATHER THAN ASSUMED NECESSARY. Six viewers each build
# a receptor VDW surface (three coloured sub-pockets + the rest of the pocket
# shell) plus ten overlaid ligand poses -- eleven models and several surface
# meshes per card, SIX times, all at page load. Verified headless
# (Chromium + SwiftShare software GL, this box has no free GPU): running all
# six simultaneously produced an intermittent "Error compiling shader" and a
# null `createViewer()` return on 1-2 of the six viewers per load, in 4 of 5
# repeated loads -- never the SAME candidate's DATA at fault (every viewer that
# DID initialise carried the correct 11 models, receptor + 10 poses), which
# points at simultaneous shader-compile contention rather than a bug in any
# one candidate's pose file. Deferring construction until each card is near
# the viewport (rather than firing all six at t=0) removes that contention on
# a constrained renderer and costs nothing on a capable one -- the intersection
# margin is generous enough that a reader scrolling at a normal pace never
# notices the deferral.
_PROMISE_OPEN = re.compile(r"\$3Dmolpromise\.then\(function\(\) \{\n")
_PROMISE_CLOSE = re.compile(r"\n\}\);\n</script>\s*$")


def offline_viewer_html(pose_html: str, uid: str) -> str:
    """One py3Dmol viewer snippet, made to run against an already-loaded
    global $3Dmol instead of racing its own CDN fetch, and deferred behind an
    IntersectionObserver so six heavy GL contexts do not all initialise in the
    same frame.

    `uid` replaces py3Dmol's own `time.time()`-derived id. NOT because that id
    collides in practice (microsecond resolution), but because "in practice"
    is not a guarantee worth trusting across six viewers built in one process,
    and a collision here means two candidates silently share one <div> and one
    renders nothing -- with no error, in the one place this page cannot afford
    that shape of failure again.
    """
    m = re.search(r"3dmolviewer_(\d+)", pose_html)
    if not m:
        raise RuntimeError("py3Dmol output did not carry the expected viewer id")
    html = pose_html.replace(m.group(1), uid)
    html, n1 = _LOADER_BLOCK.subn("", html)
    html, n2 = _PROMISE_OPEN.subn(f"function initViewer_{uid}() {{\n", html)
    html, n3 = _PROMISE_CLOSE.subn(
        "\n}\n"
        "if (window.IntersectionObserver) {\n"
        f"  var _io_{uid} = new IntersectionObserver(function(entries) {{\n"
        "    entries.forEach(function(e) {\n"
        f"      if (e.isIntersecting) {{ initViewer_{uid}(); _io_{uid}.disconnect(); }}\n"
        "    });\n"
        "  }, {rootMargin: '800px 0px'});\n"
        f"  _io_{uid}.observe(document.getElementById('3dmolviewer_{uid}'));\n"
        "} else {\n"
        f"  initViewer_{uid}();\n"
        "}\n"
        "</script>",
        html)
    if (n1, n2, n3) != (1, 1, 1):
        raise RuntimeError(
            f"py3Dmol's HTML shape changed (loader={n1}, promise-open={n2}, "
            f"promise-close={n3}); the offline-stripping regex needs updating")
    if "3dmolwarning_" not in html or "removeChild(warn)" not in html:
        raise RuntimeError(
            "the warning-placeholder removal code was stripped along with the "
            "loader -- the 'failed to load' message would stay visible")
    return html


# ---------------------------------------------------------------------------
# HTML assembly
# ---------------------------------------------------------------------------

CSS = """
:root {
  --bg: #0e1116; --panel: #161b22; --panel2: #1c2330; --text: #e6edf3;
  --muted: #9aa7b5; --border: #2a3240; --accent: #4fa3ff;
  --warn-bg: #3a1f10; --warn-border: #d97b2b; --warn-text: #ffd9b3;
  --fail-bg: #3a1414; --fail-border: #d94f4f; --fail-text: #ffc2c2;
  --ok-accent: #2bd97c;
}
@media (prefers-color-scheme: light) {
  :root { --bg:#f5f7fa; --panel:#ffffff; --panel2:#f0f3f7; --text:#1a2129;
          --muted:#5b6672; --border:#d8dee5; --accent:#0b63c5;
          --warn-bg:#fff2e2; --warn-border:#e08a2e; --warn-text:#6b3d05;
          --fail-bg:#fdecec; --fail-border:#d94f4f; --fail-text:#7a1717;
          --ok-accent:#137a3f; }
}
* { box-sizing: border-box; }
body { background: var(--bg); color: var(--text); font-family: -apple-system,
  Segoe UI, Helvetica, Arial, sans-serif; margin: 0; padding: 0 0 4em 0; }
.wrap { max-width: 1180px; margin: 0 auto; padding: 1.5em; }
h1 { font-size: 1.7em; margin-bottom: 0.1em; }
h2 { font-size: 1.25em; margin-top: 2em; border-bottom: 1px solid var(--border);
  padding-bottom: 0.3em; }
.subtitle { color: var(--muted); margin-top: 0; }
.box { border: 1px solid var(--border); border-radius: 10px; padding: 1.1em 1.3em;
  margin: 1.2em 0; background: var(--panel); }
.box.warn { background: var(--warn-bg); border-color: var(--warn-border);
  color: var(--warn-text); }
.box.fail { background: var(--fail-bg); border-color: var(--fail-border);
  color: var(--fail-text); }
.box h3 { margin-top: 0; }
table { border-collapse: collapse; width: 100%; font-size: 0.92em; margin: 0.6em 0; }
th, td { border: 1px solid var(--border); padding: 0.35em 0.6em; text-align: left; }
th { background: var(--panel2); }
code, .mono { font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; }
.smiles { word-break: break-all; font-size: 0.85em; }
.pill { display: inline-block; padding: 0.15em 0.6em; border-radius: 999px;
  font-size: 0.78em; font-weight: 600; margin-right: 0.4em; }
.pill.bdhi { background: #6a3d0f; color: #ffd9a0; }
.pill.snar { background: #6a1414; color: #ffb3b3; }
.pill.cation { background: #123a5e; color: #a9d4ff; }
.pill.viable { background: #123a1e; color: #9be8b0; }
.card { border: 1px solid var(--border); border-radius: 12px; background: var(--panel);
  margin: 1.6em 0; padding: 1.2em; }
.card.snar-card { border-color: var(--fail-border); box-shadow: 0 0 0 1px var(--fail-border) inset; }
.card-head { display: flex; flex-wrap: wrap; align-items: baseline; gap: 0.6em;
  justify-content: space-between; }
.rank-badge { font-size: 1.4em; font-weight: 800; color: var(--accent); }
.grid2 { display: grid; grid-template-columns: minmax(280px, 380px) 1fr; gap: 1.2em;
  margin-top: 0.8em; }
@media (max-width: 900px) { .grid2 { grid-template-columns: 1fr; } }
.struct-panel { background: var(--panel2); border-radius: 8px; padding: 0.6em;
  text-align: center; }
.struct-panel svg { max-width: 100%; height: auto; background: white; border-radius: 4px; }
.viewer-panel { background: #05070a; border-radius: 8px; overflow: hidden; }
.metrics { display: grid; grid-template-columns: repeat(auto-fit, minmax(150px,1fr));
  gap: 0.6em; margin: 0.8em 0; }
.metric { background: var(--panel2); border-radius: 8px; padding: 0.5em 0.7em; }
.metric .label { color: var(--muted); font-size: 0.75em; text-transform: uppercase;
  letter-spacing: 0.03em; }
.metric .value { font-size: 1.15em; font-weight: 700; }
.small { font-size: 0.85em; color: var(--muted); }
.controls-note { font-size: 0.8em; color: var(--muted); margin-top: 0.4em; }
.formula { font-family: ui-monospace, monospace; background: var(--panel2);
  padding: 0.8em 1em; border-radius: 8px; font-size: 0.95em; overflow-x: auto; }
footer { color: var(--muted); font-size: 0.85em; margin-top: 3em;
  border-top: 1px solid var(--border); padding-top: 1em; }
a { color: var(--accent); }
"""


def fmt_enrich(row) -> str:
    return (f"{row['nac_enrichment']:.2f}&times; "
            f"[{row['nac_enrichment_lo']:.2f}, {row['nac_enrichment_hi']:.2f}], "
            f"N={int(row['nac_run_count'])} runs")


def class_pill(warhead_class: str) -> str:
    cls = "bdhi" if warhead_class.startswith("bdhi") else "snar"
    return f'<span class="pill {cls}">{warhead_class}</span>'


def candidate_card(rank: int, row, uid: str) -> str:
    cid = row["candidate_id"]
    smi = row["canonical_smiles"]
    wclass = row["warhead_class"]
    is_snar = wclass == "snar_chloroazine"
    # STRUCTURAL charge, not the pH-rule flag. `charge_class == "cation"` fires
    # for four of the six candidates here, driven by obabel protonating any
    # aliphatic tertiary amine by rule at pH 7.4 -- and the med-chem workup on
    # the rank-1 candidate found that call is LIKELY WRONG for this scaffold
    # (docs/medchem_t4_72f5671e89cb.md Sec 2: the amine is alpha to two
    # electron-poor rings, its real pKa is probably well below 7.4, the
    # molecule is most likely neutral). Only `t4_9a973be6b946` carries a charge
    # that is not a rule call at all -- an aromatic pyridinium, `[nH+]`,
    # written directly into the canonical SMILES -- so that is the one flagged
    # as "charged" here. Using the frame's `charge_class` column instead would
    # have put an inaccurate "protonated pyridinium" label on three molecules
    # that have no pyridinium at all.
    has_explicit_pyridinium = "[nH+]" in smi

    smarts = depict.warhead_smarts(wclass)
    svg = depict.svg(smi, highlight_smarts=smarts, width=340, height=280,
                      legend=cid)

    pose_path = Path(row["nac_pose_path"])
    receptor = pose3d.receptor_for("nac_pose_path")
    raw_viewer = pose3d.pose_html(
        pose_path, show=tuple(range(1, 11)), width=None, height=460,
        surface=True, label_subpockets=True, show_covalent=True,
        cartoon=False, zoom_on="ligand", receptor=receptor)
    viewer = offline_viewer_html(raw_viewer, uid)

    poses = pose_properties(pose_path)
    n_viable = sum(1 for p in poses if p.get("nac_viable") == "True")
    pose_rows = "\n".join(
        f'<tr class="{"viable-row" if p["nac_viable"]=="True" else ""}">'
        f'<td>{p["pose_rank"]}</td>'
        f'<td>{"viable" if p["nac_viable"]=="True" else "&mdash;"}</td>'
        f'<td>{float(p["nac_energy_kcal"]):.2f}</td>'
        f'<td>{float(p["nac_distance_A"]):.2f}</td>'
        f'<td>{float(p["nac_angle_deg"]):.1f}&deg; ({p.get("nac_angle_kind","")})</td>'
        "</tr>"
        for p in poses)

    flags = []
    if has_explicit_pyridinium:
        flags.append(
            '<span class="pill cation">charged &mdash; aromatic pyridinium '
            '(<code>[nH+]</code>) written directly into the SMILES, not a pH-rule '
            "call</span>")
    elif row.get("charge_class") == "cation":
        flags.append(
            '<span class="pill" style="background:var(--panel2);color:var(--text)">'
            "obabel's pH 7.4 rule flags this molecule's aliphatic amine as "
            "protonated &mdash; likely an over-call for this scaffold (the amine sits "
            "alpha to two electron-poor rings), see "
            "<code class=\"mono\">docs/medchem_t4_72f5671e89cb.md</code> &sect;2, "
            "not independently confirmed as the physiological charge state</span>")
    if is_snar:
        flags.append('<span class="pill snar">nitro group + chloropyrimidine leaving group</span>')

    snar_caveat = ""
    if is_snar:
        snar_caveat = """
        <div class="box fail">
          <h3>This card's warhead class was TESTED and FAILED validation &mdash; worse than "unvalidated"</h3>
          <p>SNAr (<code>snar_chloroazine</code>) was measured against 2 crystallographic Cys113
          positives and warhead-matched inactives: <b>AUC 0.558, p = 0.41</b> (D0065). The
          robustness control redrew the negative pool nine independent times and every single
          draw straddled chance: <b>AUC 0.451</b>, range 0.375&ndash;0.525 (D0070/ranking_rationale.md).
          That is <b>no signal</b>, not merely an underpowered measurement.</p>
          <p>Worse: this class's <b>negatives enrich 2.39&times; over chance on their own</b> &mdash;
          the pocket appears to steer most chloroazines into a perpendicular approach whether or
          not they bind, so a high SNAr enrichment number may be measuring geometry of the pocket
          rather than anything about this molecule. Only two SNAr ligands have ever been
          crystallised at Cys113 (n=2), so this cannot be settled either way with more of the
          same kind of evidence.</p>
          <p><b>This molecule is on the page because it cleared the numeric threshold
          (&gt;5.70&times;), not because there is evidence the number means anything for this
          warhead class.</b> Weight its rank accordingly &mdash; it should not be read as
          equivalent to the BDHI entries above it, whose caveat (unvalidated, not failed) is a
          different and lesser problem.</p>
        </div>"""

    return f"""
    <div class="card {'snar-card' if is_snar else ''}" id="card-{cid}">
      <div class="card-head">
        <div><span class="rank-badge">#{rank}</span>
          &nbsp; <span class="mono" style="font-size:1.15em">{cid}</span></div>
        <div>{class_pill(wclass)}
          <span class="pill" style="background:var(--panel2);color:var(--text)">
            mechanism: {row['warhead_mechanism']}</span></div>
      </div>
      <div class="small">{' '.join(flags) if flags else ''}</div>

      <div class="metrics">
        <div class="metric"><div class="label">NAC enrichment (200 runs)</div>
          <div class="value">{row['nac_enrichment']:.2f}&times;</div>
          <div class="small">95% CI [{row['nac_enrichment_lo']:.2f}, {row['nac_enrichment_hi']:.2f}]</div></div>
        <div class="metric"><div class="label">median S&middot;&middot;&middot;C distance</div>
          <div class="value">{row['nac_median_dist']:.2f} &#8491;</div></div>
        <div class="metric"><div class="label">median approach angle</div>
          <div class="value">{row['nac_median_angle']:.1f}&deg;</div></div>
        <div class="metric"><div class="label">poses shown / NAC-viable</div>
          <div class="value">{len(poses)} / {n_viable}</div></div>
      </div>

      <div class="smiles mono"><b>SMILES</b> &nbsp;{smi}</div>

      {snar_caveat}

      <div class="grid2">
        <div class="struct-panel">
          <div class="small" style="margin-bottom:0.4em">2D structure &mdash; warhead
            highlighted (SMARTS <code class="mono">{smarts}</code>)</div>
          {svg}
        </div>
        <div class="viewer-panel">
          {viewer}
        </div>
      </div>
      <div class="controls-note">3D viewer: left-drag rotate &middot; scroll / right-drag
        zoom &middot; Ctrl+left-drag pan. Receptor: <b>3IKD</b> (chemist-prepared, reactive
        near-attack docking, D0059/D0063/D0064) &mdash; <b>not</b> 6VAJ. All ten saved poses
        overlaid; pose 1 (thickest sticks) is the best NAC-viable pose, the rest follow by
        energy. The dashed yellow line marks any pose within {pose3d.COVALENT_BOND_MAX_A} &#8491;
        of Cys113 SG.</div>

      <details style="margin-top:0.8em">
        <summary class="small">Per-pose geometry (all {len(poses)} poses)</summary>
        <table>
          <tr><th>pose</th><th>NAC</th><th>energy (kcal/mol)</th>
              <th>dist. to Cys113 SG (&#8491;)</th><th>approach angle</th></tr>
          {pose_rows}
        </table>
      </details>
    </div>
    """


def build(frame: Path, df: pd.DataFrame) -> str:
    js = vendored_3dmol_js()
    n_bdhi = int((df["warhead_class"].str.startswith("bdhi")).sum())
    n_snar = int((df["warhead_class"] == "snar_chloroazine").sum())
    n_c5 = int((df["warhead_class"] == "bdhi_c5").sum())
    n_c4 = int((df["warhead_class"] == "bdhi_c4").sum())

    cards = "\n".join(
        candidate_card(i + 1, row, f"cand{i}")
        for i, (_, row) in enumerate(df.iterrows()))

    sn2_null = nac.isotropic_null("sn2_displacement")
    perp_null = nac.isotropic_null("michael_addition")

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>T&#8324; top candidates &mdash; near-attack enrichment catalogue</title>
<style>{CSS}</style>
</head>
<body>
<script>
{js}
</script>
<div class="wrap">

<h1>T&#8324; top candidates by near-attack enrichment</h1>
<p class="subtitle">{len(df)} molecules with enrichment &gt; {NAC_ENRICHMENT_MIN:.2f}&times;,
  latest frame <code class="mono">{frame.name}</code> ({frame.parent.name}), read
  {pd.Timestamp.now(tz='UTC').strftime('%Y-%m-%d %H:%M UTC')}. All poses docked into the
  <b>3IKD</b> receptor. Built for @tt8804.</p>

<div class="box fail">
  <h3>Read this before the molecules below &mdash; the head of this ranking does not
    survive more searching</h3>
  <p>Every enrichment number on this page is measured at <b>200</b> independent docking
  runs. That is not a stylistic choice of run count: <a href="#formula">the formula
  below</a> shows the score is literally a fraction of those 200 runs, so it is a
  property of the search as well as of the molecule (D0068).</p>
  <p>The top candidate on this page and two matched siblings (same core, same R-group,
  only the warhead differs) were independently re-docked at <b>2,000</b> runs
  (10&times; the search) for the med-chem workup. The result:</p>
  <table>
    <tr><th>molecule</th><th>class</th><th>200 runs</th><th>2,000 runs</th>
        <th>median dist. to Cys113 SG</th></tr>
    <tr><td class="mono">t4_72f5671e89cb</td><td>bdhi_c5 (rank&nbsp;1, this page)</td>
        <td>6.74&times;</td><td><b>0.60&times;</b> &mdash; below chance</td>
        <td>3.15 &rarr; 5.17 &#8491;</td></tr>
    <tr><td class="mono">t4_9c44a3f8892a</td><td>bdhi_c4 (same R-group, not on this page)</td>
        <td>4.72&times;</td><td><b>0.94&times;</b></td>
        <td>3.89 &rarr; 5.52 &#8491;</td></tr>
    <tr><td class="mono">t4_45901f30d2a1</td><td>chloroacetamide (same R-group, control)</td>
        <td>1.12&times;</td><td>0.94&times;</td>
        <td>3.45 &rarr; 3.58 &#8491; (barely moved)</td></tr>
  </table>
  <p>The BDHI ligands drift roughly <b>2 &#8491; away</b> from the sulfur once the search is
  allowed to converge; the chloroacetamide control barely moves (0.12 &#8491;). At
  2,000 runs all three sit at or below chance, and the entire 200-run gap that put
  <code>t4_72f5671e89cb</code> at rank 1 of 5,765 scored T&#8324; molecules collapses to
  nothing.</p>
  <p><b>Why this had never been caught:</b> D0068 established that enrichment does not
  converge, but its convergence set (15 crystallographic positives + the top-300 refine)
  was entirely chloroacetamide, naphthoquinone and SNAr &mdash; <b>no BDHI molecule had
  ever been run at 2,000 runs before this check.</b> BDHI is also the class that now
  dominates the ranking: the T&#8324; top 100 is <b>64 bdhi_c5 + 16 bdhi_c4 + 17 SNAr, and
  only 2 chloroacetamide</b>, from a pool with exactly 187 of each class &mdash; i.e. the
  head of the list is ~97% from a class whose convergence was never checked (BDHI) or a
  class that failed validation outright (SNAr, see below).</p>
  <p><b>What follows for this page:</b> the six molecules below are real, their poses are
  real, and they are worth looking at &mdash; that is what was asked for. <b>They are not,
  on the strength of anything measured so far, a shortlist to synthesise.</b> They are the
  top of a ranking whose top is now known to be substantially an artefact of
  under-sampling for five of six molecules, and untested-for-signal for the sixth.</p>
</div>

<h2 id="formula">What "enrichment" means, exactly</h2>
<div class="box">
  <div class="formula">enrichment&nbsp; = &nbsp;(fraction of independent docking runs reaching a
viable near-attack conformation)
              &divide; (fraction a randomly-oriented approach would reach by chance)</div>
  <p><b>Numerator, measured.</b> The molecule is docked <b>N = 200</b> times independently
  into 3IKD (every number on this page is at N=200 &mdash; see the warning above for why that
  matters). A pose counts as <i>viable</i> only if <b>both</b> hold, from
  <code class="mono">shared/nac_criterion.py</code>:</p>
  <ul>
    <li><b>distance</b> &mdash; warhead atom to Cys113 SG, <b>{nac.NAC_DIST_MIN}&ndash;{nac.NAC_DIST_MAX} &#8491;</b>,
      a van der Waals <i>contact</i> (the reactant state), not a formed bond</li>
    <li><b>angle</b> &mdash; mechanism-specific:
      <table>
        <tr><th>mechanism</th><th>angular criterion</th><th>chance baseline (exact)</th></tr>
        <tr><td>SN2 (<code>sn2_displacement</code>)</td>
            <td>S&middot;&middot;&middot;C&ndash;LG &ge; {nac.SN2_ANGLE_MIN:.0f}&deg;, backside anti to
              the leaving group</td>
            <td><b>{sn2_null*100:.2f}%</b></td></tr>
        <tr><td>Michael / SNAr / BDHI ring-opening (all <code>perpendicular_to_plane</code>)</td>
            <td>&le; {nac.PERPENDICULAR_MAX_OFF_NORMAL:.0f}&deg; off the sp2 plane normal <b>and</b>
              Burgi&ndash;Dunitz approach {nac.APPROACH_WINDOW[0]:.0f}&ndash;{nac.APPROACH_WINDOW[1]:.0f}&deg;</td>
            <td><b>{perp_null*100:.2f}%</b></td></tr>
      </table>
    </li>
  </ul>
  <p><b>Denominator, computed exactly</b> (not sampled) &mdash; the fraction of all possible
  approach directions that would satisfy the angular window above, from its solid angle
  (<code class="mono">nac_criterion.isotropic_null()</code>).</p>
  <p><b>Reading:</b> how many times more often than chance this molecule reaches a
  chemically competent geometry. <b>1.0&times; = no better than a random approach.</b></p>
  <p><b>Why divide at all:</b> five of the six molecules below are BDHI and one is SNAr
  &mdash; two different angular windows, {perp_null*100:.1f}% vs {sn2_null*100:.2f}% permissive by
  solid angle alone (a {perp_null/sn2_null:.1f}&times; difference). Raw viable fractions are
  therefore <b>not comparable across mechanisms</b>; dividing each by its own chance
  baseline is what makes a BDHI 6.43&times; and an SNAr 5.76&times; describe the same kind of
  thing. Pooling the raw (undivided) fractions instead gives AUC &asymp; 0.5 &mdash; an artefact
  of the windows, not a result (<code>docs/ranking_rationale.md</code>).</p>
  <p class="small"><b>Reference points, same gate</b> (D0065): crystallographic Cys113
  binders score <b>1.6&ndash;4.3&times;</b>; random warhead-matched <i>measured</i> inactives
  score <b>&asymp;0.8&times;</b>. Placed here so the six numbers below can be read against
  something other than each other.</p>
</div>

<h2>Three things this list is, plainly</h2>
<div class="box">
  <ol>
    <li><b>The top of this ranking is dominated by one small part of chemical space, not
      six independent findings.</b> {n_bdhi} of {len(df)} candidates below are BDHI
      ({n_c5} <code>bdhi_c5</code> + {n_c4} <code>bdhi_c4</code>); the sixth is a single
      SNAr chloroazine. Across the full T&#8324; top 100 the same skew holds at far larger
      scale (64 bdhi_c5 + 16 bdhi_c4 + 17 SNAr + 2 chloroacetamide, from a pool with 187 of
      each class) &mdash; this page's composition is not a coincidence of six molecules, it
      is what the ranking currently does at scale.</li>
    <li><b>BDHI (5 of 6 candidates) has ZERO crystallographic Cys113 positives.</b> The 17
      solved covalent Cys113 structures split chloroacetamide (10), naphthoquinone (4),
      SNAr (2) &mdash; and BDHI, none. Unlike chloroacetamide (validation AUC 0.908) and
      Michael/naphthoquinone (0.734), this warhead class is <b>unvalidated</b> in this
      project's hands: not shown wrong, simply never checked against a known binder.</li>
    <li><b>These BDHI scores exist because of a defect fixed TODAY (D0067).</b> Until this
      morning, BDHI was scored with backside (sp3) attack geometry at what is actually an
      sp2 carbon (the C of its C=N), and every BDHI candidate read 0.00&times; &mdash;
      completely unreactive. The fix (perpendicular approach, matching addition&ndash;
      elimination chemistry at an sp2 centre) is independently corroborated by the
      literature (Byun 2023, DFT, addition&ndash;elimination at the sp2 C3 carbon) &mdash; so
      the geometry the fix now uses is very likely the chemically correct one. <b>What is
      NOT independently confirmed is these specific ranks</b>: no experiment, crystal
      structure, or (per the warning above) even a converged docking run backs the ordering
      the fix produced. It corrected a defect; it did not validate a shortlist.</li>
  </ol>
</div>

<h2>The six candidates</h2>
{cards}

<footer>
  <p><b>Source frame:</b> <code class="mono">{frame}</code> &mdash; resolved as the
  highest-numbered <code>D4_&lt;N&gt;.parquet</code> under <code>04_t4_combinatorial/</code>
  at build time, not hardcoded.</p>
  <p><b>Poses:</b> <code class="mono">00_outputs/blacksmith/nac_poses/&lt;candidate_id&gt;.sdf</code>,
  10 saved modes per candidate, best-NAC-viable-first (<code>scripts/export_nac_poses.py</code>).</p>
  <p><b>Receptor:</b> 3IKD, chemist-prepared, reactive near-attack docking protocol
  (D0059, D0063, D0064). This is <b>not</b> 6VAJ.</p>
  <p><b>Decisions referenced</b> (repo-relative, not linked &mdash; this page is meant
  to be read offline): <code class="mono">decisions/D0067</code> (BDHI sp2 geometry
  fix), <code class="mono">D0068</code> (enrichment does not converge, energy beats
  geometry at convergence), <code class="mono">D0069</code> (plain docking on 3IKD
  outperforms the geometric criterion), <code class="mono">D0070</code> (consensus
  vs. frequency under convergence) &mdash; see <code class="mono">decisions/</code>
  and <code class="mono">docs/ranking_rationale.md</code> in the
  <code class="mono">inhibition</code> repo. Med-chem convergence workup:
  <code class="mono">docs/medchem_t4_72f5671e89cb.md</code>.</p>
  <p>Built by the Artist agent (murmurent) with
  <code class="mono">scripts/build_t4_top_catalogue.py</code>, reusing
  <code class="mono">integration/app/pose3d.py</code> and
  <code class="mono">integration/app/depict.py</code> rather than re-implementing pose
  rendering or 2D depiction. 3Dmol.js 2.5.5 vendored inline (BSD-3) &mdash; no CDN
  dependency, renders offline.</p>
</footer>
</div>

</body>
</html>
"""


def main() -> None:
    frame = latest_frame()
    df = top_candidates(frame)
    print(f"latest frame: {frame}")
    print(f"{len(df)} candidates above {NAC_ENRICHMENT_MIN}x:")
    for _, r in df.iterrows():
        print(f"  {r['candidate_id']:<20} {r['nac_enrichment']:.2f}x  {r['warhead_class']}")

    html = build(frame, df)

    out = outputs.Topic("artist").write("t4_top_candidates_catalogue", ".html")
    out.write_text(html, encoding="utf-8")
    print(f"wrote {out}  ({out.stat().st_size/1e6:.2f} MB)")

    # Convenience copy for @tt8804 -- same shared-group storage as the rest of
    # this project's deliverables (/data/lab_vm/modifiable/inhibition/), but at
    # a fixed, unversioned name so it is easy to point at without knowing the
    # append-only version number. NOT a second source of truth: it is
    # overwritten every run and the append-only copy above is canonical.
    convenience_dir = Path("/data/lab_vm/modifiable/inhibition/deliverables")
    convenience_dir.mkdir(parents=True, exist_ok=True)
    convenience = convenience_dir / "t4_top_candidates_catalogue.html"
    convenience.write_text(html, encoding="utf-8")
    print(f"convenience copy: {convenience}")


if __name__ == "__main__":
    main()
