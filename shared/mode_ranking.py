"""The ranking view: every molecule and every mode, before anything is simulated.

TWO VIEWS, AND THIS IS THE FIRST ONE. `combined.html` shows the sweep and the
100 ns results -- the subset that was simulated. This shows what the screen
*scored*: every molecule, every mode, ranked, with its pose. In the pipeline's
real order this comes FIRST -- you read the ranked list, look at the poses, and
then choose what goes to the sweep. It was built retrospectively, after #53 found
that the sweep took mode 0 for 242 of 242 molecules while the ranking is per
mode, and that gap was invisible precisely because no view like this existed.

SCOPE IS ONE CONTROL, NOT TWO. A dropdown picks a warhead class, every class
ranked within itself, or one global order. Two separate toggles let a reader
combine "global" with a class filter, a combination with no meaning.

Within-class is the default. The SN2 angular criterion is far stricter than the
perpendicular one (#47), so a global order compares scores computed under
different bars. It is offered because "where does this sit overall" is a real
question, and refusing to answer it does not remove the bias -- the option names
it instead.

THE JOIN IS ON (parent_ident, mode). Never on `ident`: mode 0 is the bare ident
in the sweep table and `_m0` in the rank table, so a merge on the label silently
drops exactly the rows that were simulated (`shared/mode_key.py`).

Depictions and poses are FILES fetched on demand, not inlined. 8,096 rows of
base64 is tens of megabytes; the results GUI can inline its 59 and this cannot.
"""

from __future__ import annotations

import glob
import logging
import os
import html
import json
from pathlib import Path

from . import run_paths as rp

import pandas as pd

from shared import mode_key as mk

B = rp.BLACKSMITH
RECEPTOR = rp.receptor_prep()
CYS_RESI = 113
#: Molecules to draw the eye to. They are ordinary candidates -- screened,
#: ranked and swept through the identical path -- and the flag is presentational
#: only. `t4_sulfopin` is the literature parent of the chloroacetamide series
#: (Dubiella 2021), added to the candidate frame so the screen has to place it
#: without knowing what it is.
HIGHLIGHT = {"t4_sulfopin"}
#: Surface shell: residues with a heavy atom within this of Cys113's SG. A
#: whole-protein VDW mesh is the expensive call in a 3Dmol viewer and the far
#: side of the protein is not the subject.
SURF_SHELL_A = 12.0


def pocket_residues() -> list[int]:
    """Residue numbers forming the pocket wall, from the receptor itself."""
    import numpy as np
    if not RECEPTOR.is_file():
        return []
    sg, res = None, {}
    for ln in RECEPTOR.read_text().splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        ri, nm = ln[22:26].strip(), ln[12:16].strip()
        try:
            xyz = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
        except ValueError:
            continue
        res.setdefault(ri, []).append(xyz)
        if ri == str(CYS_RESI) and nm == "SG":
            sg = xyz
    if sg is None:
        return []
    return sorted(int(ri) for ri, xs in res.items()
                  if ri.lstrip("-").isdigit()
                  and min(float(np.linalg.norm(x - sg)) for x in xs) <= SURF_SHELL_A)


def _latest(pattern: str) -> Path | None:
    fs = sorted(glob.glob(str(B / pattern)))
    return Path(fs[-1]) if fs else None


def gather() -> pd.DataFrame:
    """One row per mode: rank, docking-derived scores, and what was simulated."""
    frames = []
    # WHICH RUN THIS IS A VIEW OF, STATED RATHER THAN INFERRED. The pattern used
    # to be `rank_v2_{tier}_{score}_*.csv` -- the un-suffixed name -- so the view
    # showed whichever screen last claimed that slot and could not report which
    # one it was. It now comes from `run.topic` in config/target.yaml, the same
    # key the ranking stamps into its filenames, so the GUI and the ranking
    # cannot end up describing different screens.
    from shared import target_config as tc
    topic = tc.get("run.topic")
    # TIER SCOPE IS READ HERE TOO, NOT ONLY IN THE RANKING. Dropping a tier from
    # `rank_v2` stops it writing a NEW file for that tier; it does not remove the
    # old one, and this reader takes the newest match for each tier
    # independently. Without the same filter the view would keep showing a T_3
    # table from before the decision, beside T_4 rows ranked without it.
    want = {str(t).upper() for t in tc.get("run.tiers", default=["T3", "T4"])}
    # THE SCORE IS NAMED IN CONFIG, NOT HARDCODED HERE. It was the literal pair
    # (T4, conditional_eb) / (T3, enrichment_conditional), so a run ranked on any
    # other score was invisible to the whole GUI: nac_v6's 327,167 groups ranked
    # on `engagement` produced "no ranked modes for topic nac_v6 yet -- nothing
    # to build", which reads as "the ranking has not run" rather than "this
    # reader is looking for a different filename". Same shape as the topic
    # literals D0080 records.
    _scores = tc.get("ranking.score_by_tier",
                     default={"T4": "conditional_eb",
                              "T3": "enrichment_conditional"})
    for tier, score in ((t, _scores[t]) for t in ("T4", "T3") if t in _scores):
        if tier not in want:
            continue
        f = _latest(f"rank_v2/rank_v2_{tier}_{topic}_{score}_*.csv")
        if f is None:
            continue
        d = pd.read_csv(f)
        if "mode" in d.columns:
            d = d[d["mode"].notna()]
        d["tier"] = d.get("tier", tier)
        frames.append(d)
    if not frames:
        return pd.DataFrame()
    r = pd.concat(frames, ignore_index=True)

    # THIS RUN'S SWEEP, NOT EVERY SWEEP EVER. These read the UNSCOPED
    # `attack_sweep/` and `md_residence/` directories, so every status badge on
    # the ranking page described the PREVIOUS screen: a mode showed "100 ns"
    # because it ran in 3.0.0, and a mode swept this morning showed "not run"
    # because nac_v5's directory was not read at all. @tt8804: "why are these
    # top ranked mols showing not run and some showing 100 ns??"
    #
    # Sorted by MTIME, not lexicographically: `keep="last"` over a lexicographic
    # sort makes `_10` older than `_9`, which decides which of two measurements
    # of the same mode is believed by string order.
    sf = sorted(glob.glob(str(_rp().sweep_dir() / "attack_sweep_*.csv")),
                key=os.path.getmtime)
    sweep = (pd.concat([pd.read_csv(x) for x in sf], ignore_index=True)
             .drop_duplicates("ident", keep="last") if sf else pd.DataFrame())
    if not sweep.empty:
        keep = [c for c in ("ident", "parent_ident", "mode", "frac_attack_ready",
                            "n_visits", "status") if c in sweep.columns]
        # ATTEMPTED is not SUCCEEDED. A row exists for every mode sent;
        # frac_attack_ready is null when the run failed. Counting only the
        # successful ones as "swept" would report a mode that was tried and
        # crashed as one nobody ever chose.
        sweep = sweep.assign(_sent=True)
        # A BARE IDENT HERE MEANS "THE MODE IS UNKNOWN", NOT "MODE 0".
        #
        # `bare_is_mode_zero` exists for sweep tables written before #53, when
        # mode 0 really was recorded bare. This topic is post-#53: an `ok` row
        # carries `<parent>_m<mode>`, and a row is bare only when the run DIED
        # BEFORE IT KNEW ITS MODE -- which is what all 24 of this run's failures
        # are, from a launcher that asked for pose ranks in the hundreds.
        #
        # Read as mode 0 they collided with the real mode-0 row of the same
        # molecule: the join duplicated those modes (4,432 rows in, 4,434 out)
        # and, worse, badged a real mode FAILED for a failure that belongs to no
        # mode at all. A status is an accusation; it has to be attributable.
        sw = mk.add_key(sweep[keep + ["_sent"]].rename(
            columns={"status": "sweep_status"}), bare_is_mode_zero=False)
        unattributable = int(sw.mode_key.isna().sum())
        if unattributable:
            logging.getLogger(__name__).info(
                "sweep: %d row(s) carry no mode and cannot label one "
                "(the run failed before it knew) — excluded", unattributable)
            sw = sw[sw.mode_key.notna()]
        r = mk.join(r, sw.drop(columns=[c for c in ("parent_ident", "mode")
                                        if c in sw.columns]),
                    suffixes=("", "_sw"))

    md_ids: set[str] = set()
    for f in glob.glob(str(_rp().residence_dir() / "*.csv")):
        try:
            d = pd.read_csv(f)
        except Exception:                                  # noqa: BLE001
            continue
        if "ident" not in d.columns or "production_ps" not in d.columns:
            continue
        d = d[(d.production_ps >= 50000)
              & d.status.astype(str).str.startswith("ok")]
        md_ids |= set(d.ident.astype(str))

    r["sent"] = (r["_sent"].fillna(False).astype(bool)
                 if "_sent" in r.columns else False)
    r["swept"] = (r["frac_attack_ready"].notna()
                  if "frac_attack_ready" in r.columns else False)
    # A mode gets the 100 ns badge only if its molecule ran AND this mode is the
    # one that was sent. The MD rows do not record their mode (#36), so anything
    # looser would badge a mode that never moved.
    r["ran_md"] = r.parent_ident.isin(md_ids) & r["sent"]

    # GLOBAL RANK IS COMPUTED HERE AND LABELLED BIASED WHEREVER IT IS SHOWN.
    # conditional_eb is not comparable across warhead classes (#47); this exists
    # so the page can answer "where does this sit overall" while saying so.
    # References carry no class_rank of their own -- they are excluded from the
    # candidates' per-class competition by design -- so give them one computed
    # the same way, against the candidates of their class, and mark it.
    # HIGHLIGHTED ROWS. Not a separate tier and not a separate score -- these ARE
    # candidates, screened, ranked and swept through the identical path, and the
    # flag only draws the eye. rank_v2's reference table is deliberately NOT read
    # here: references are scored on a different column (enrichment_conditional,
    # because conditional_eb is not computed for them) and carry no rank slot, so
    # putting them in this list would mix two scores in one ordering. ATRA is the
    # concrete reason to keep them out: it is `mechanism_declared: non_covalent`
    # and matched three covalent warhead classes on the loose Michael SMARTS
    # `[CX3]=[CX3][CX3]=O`, which its alpha,beta-unsaturated carboxylic ACID
    # satisfies. A conjugated carboxylate is not a warhead.
    r["is_control"] = r.parent_ident.isin(HIGHLIGHT)

    if "conditional_eb" in r.columns:
        r["global_rank"] = r["conditional_eb"].rank(
            ascending=False, method="min", na_option="bottom")
    return r


def idents(r: pd.DataFrame) -> set[str]:
    """The molecules the view will ask for assets for."""
    return set(r.parent_ident.astype(str)) if not r.empty else set()


def _rows_json(r: pd.DataFrame) -> str:
    out = []
    for _, x in r.iterrows():
        # UNRANKED MODES ARE EMITTED, NOT DROPPED. A mode with
        # viable_fraction == 0 has no enrichment_conditional and therefore no
        # class_rank, and this used to `continue` past it. That silently removed
        # 1,869 modes -- 23% of the library -- from the ranking view, and it
        # removed them NON-RANDOMLY: always the zero-viable ones, which are
        # exactly the stratum the sweep-depth question is about (#59).
        #
        # It also made the per-molecule table lie. t3_5b92831a5d23 has SIX modes
        # and the page said four, because m2 and m4 are zero-viable (@tt8804
        # spotted it). "Stamp, do not drop" -- gui_spec.md §6.2.
        pass
        st = ("md" if x.ran_md else "swept" if x.swept
              else "failed" if x.sent else "none")

        def num(k, nd=None):
            v = x.get(k)
            if pd.isna(v):
                return None
            return round(float(v), nd) if nd is not None else int(v)

        # `i` IS NOT SENT. It is exactly `p + "_m" + m` for all 34,076 rows
        # (checked, 0 exceptions) and is rebuilt once on load. At ~30 bytes a row
        # it was a megabyte of the payload restating two fields already present.
        # `ctl` is omitted when false for the same reason -- it is true for the
        # control alone, and `"ctl":false` on every other row cost 400 KB.
        row = {
            "p": str(x.parent_ident), "m": int(x["mode"]),
            "c": str(x.warhead_class),
            "cr": (int(x.class_rank) if pd.notna(x.get("class_rank")) else None),
            "gr": num("global_rank"), "n": num("n_poses_mode"),
            "np": num("n_poses"), "vf": num("viable_fraction", 4),
            "eb": num("conditional_eb", 3), "en": num("enrichment", 2),
            "sp": num("spread_a", 2), "dc": num("dir_coherence", 3),
            "fa": num("frac_attack_ready", 4), "s": st,
            # `mode_label` is 1a / 1b when a first-stage mode was subdivided
            # (#61) and a plain number otherwise. Absent on every frame screened
            # before sub-splitting existed, so it falls back to the number.
            "ml": (str(x["mode_label"]) if pd.notna(x.get("mode_label"))
                   else str(int(x["mode"]))),
            "pm": (int(x["parent_mode"]) if pd.notna(x.get("parent_mode"))
                   else int(x["mode"])),
        }
        if bool(x.get("is_control", False)):
            row["ctl"] = True
        out.append(row)
    return json.dumps(out, separators=(",", ":"))


_TPL = r"""<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>__TITLE__ — ranking</title>
<script>__THREE__</script>
<style>
/* The same palette and the same shell as the results GUI, so the two views read
   as one instrument rather than two pages that happen to link to each other. */
:root{--ink:#10233f;--navy:#003087;--blue:#0072ce;--blue-pale:#e8f1fb;
 --rule:#d6dee8;--muted:#5b6b80;--paper:#fff;--raise:#f5f8fc;--rail:#fafcfe;
 --good:#0f7a54;--bad:#b3261e;--warn:#8a6d1f;
 --sans:"Helvetica Neue",Helvetica,Arial,system-ui,sans-serif;
 --mono:ui-monospace,"SF Mono",Menlo,Consolas,monospace}
:root[data-theme="dark"]{--ink:#dfe7f0;--navy:#8ab4e8;--blue:#6aa9e0;
 --blue-pale:#16283a;--rule:#25333f;--muted:#93a3b4;--paper:#0e151c;
 --raise:#16202a;--rail:#121b24;--good:#4fc4a0;--bad:#e08a70;--warn:#d0ae5a}
*{box-sizing:border-box}
html,body{height:100%}
body{margin:0;background:var(--paper);color:var(--ink);font-family:var(--sans);
 font-size:14px;line-height:1.5;font-variant-numeric:tabular-nums;
 display:flex;flex-direction:column}
#topbar{display:flex;align-items:center;gap:8px;padding:6px 14px;min-height:38px;
 border-bottom:1px solid var(--rule);background:var(--raise);
 overflow-x:auto;white-space:nowrap;scrollbar-width:thin}
h1{margin:0;font-size:.86rem;font-weight:600;letter-spacing:-.01em;color:var(--navy)}
.mbtn{font:600 11px var(--sans);padding:3px 10px;border-radius:99px;cursor:pointer;
 border:1px solid var(--rule);background:var(--paper);color:var(--ink)}
.mbtn.on{background:var(--navy);color:#fff;border-color:var(--navy)}
.mbtn.lnk{text-decoration:none;color:var(--blue)}
select#scope{font:600 11px var(--sans);padding:3px 26px 3px 10px;border-radius:99px;
 border:1px solid var(--rule);background:var(--paper);color:var(--ink);cursor:pointer;
 appearance:none;background-image:linear-gradient(45deg,transparent 50%,var(--muted) 50%),
 linear-gradient(135deg,var(--muted) 50%,transparent 50%);
 background-position:calc(100% - 14px) 52%,calc(100% - 9px) 52%;
 background-size:5px 5px,5px 5px;background-repeat:no-repeat}
select#scope:focus{outline:2px solid var(--blue);outline-offset:1px}
.msep{width:1px;height:16px;background:var(--rule);flex:none}
.mhint{font-size:11px;color:var(--muted);margin-left:4px}
main{flex:1;display:grid;grid-template-columns:376px 1fr;min-height:0}
@media(max-width:880px){main{grid-template-columns:1fr;grid-template-rows:250px 1fr}}
/* The rail's COLUMN is the grid child; the banner and the scroller stack inside
   it. min-height:0 on both, or the scroller refuses to shrink and the whole
   column grows to the height of the list. */
#railcol{display:flex;flex-direction:column;min-height:0;min-width:0;
 border-right:1px solid var(--rule);background:var(--rail)}
/* `contain:layout paint`, NEVER `strict`. `strict` adds SIZE containment, which
   makes the element size itself without regard to its contents -- inside a grid
   row that is not explicitly sized, the rail then contributes zero height and
   collapses to a strip. Layout and paint containment are what is wanted here:
   they stop the virtualised window from invalidating the rest of the page. */
#rail{flex:1;overflow-y:auto;background:var(--rail);
 position:relative;contain:layout paint}
#railPad{position:relative;width:100%}
#railWin{position:absolute;top:0;left:0;right:0;will-change:transform}
/* Fixed geometry is what makes the window computable without measuring 34,076
   elements. Both are enforced, not assumed -- a row that overflowed its box
   would drift the whole list out of register with the scrollbar. */
.row{height:64px;box-sizing:border-box;overflow:hidden}
/* AND THIS PAGE'S OWN VOCABULARY, which the shared row does not carry and must
   not: the rail pages badge a mode held/left against an RMSD bar, this one
   badges HOW FAR DOWN THE CASCADE it got. Extracting ROW_CSS deleted these
   along with the duplicated rules, and the rows rendered their mode badge and
   status tag unstyled -- present, legible, and wrong, which is why it read as
   "the selectors look different" rather than as anything broken. */
.mtag{display:inline-block;padding:0 .3rem;margin-left:.3rem;border-radius:3px;
  color:#fff;font:600 10px var(--mono);vertical-align:1px}
.t-md{background:#e6f4ee;color:var(--good)}
.t-swept{background:#e8f1fb;color:var(--navy)}
.t-failed{background:#faf3e0;color:var(--warn)}
.t-none{background:var(--raise);color:var(--muted)}
/* The class banner that replaces sticky headers under virtualisation. */
.railbn{font-family:var(--mono);font-size:.6rem;letter-spacing:.14em;
 text-transform:uppercase;color:var(--blue);font-weight:600;padding:8px 14px 6px;
 background:var(--raise);border-bottom:1px solid var(--rule)}
.chd{font-family:var(--mono);font-size:.6rem;letter-spacing:.14em;
 text-transform:uppercase;color:var(--blue);font-weight:600;padding:10px 14px 6px;
 background:var(--raise);border-bottom:1px solid var(--rule);position:sticky;top:0;z-index:1}
__ROWCSS__
__RAILQCSS__
/* Only what VIRTUALISATION requires, on top of the shared row. Every item is
   exactly ROW_H tall, because the window offset is computed as i * ROW_H
   rather than measured -- a row free to size itself would put every row below
   it at the wrong offset. */
.row{height:64px;box-sizing:border-box;overflow:hidden}
#viewer{min-width:0;min-height:0;display:flex;flex-direction:column;background:var(--paper)}
#vhead{padding:9px 18px;border-bottom:1px solid var(--rule);background:var(--raise);
 display:flex;justify-content:space-between;align-items:center;gap:12px}
#vname{font-family:var(--mono);font-size:12.5px;font-weight:600;color:var(--navy)}
#vbody{flex:1;min-height:0;overflow-y:auto;padding:16px 18px 30px}
.facts{display:grid;grid-template-columns:repeat(auto-fit,minmax(128px,1fr));
 gap:10px 18px;margin-bottom:14px}
.fact b{display:block;font-family:var(--mono);font-size:1.05rem}
.fact span{font-size:11px;color:var(--muted)}
/* THE DEPICTION IS CONTAINED, NOT TRUSTED (#62). @tt8804 reported the structure
   overlapping the facts row above it. The thumbs themselves are clean -- checked
   over all 561, none draws outside its own viewBox -- so the collision is a
   layout property, and the fix is to make the box unable to be overrun rather
   than to chase one molecule's geometry:
     - `overflow:hidden` so nothing paints outside the frame, whatever the SVG;
     - `display:block` on the img, because an inline image sits on a text
       baseline and contributes descender space that shifts what follows;
     - `max-height` so a depiction with an extreme aspect ratio cannot grow the
       panel without bound;
     - `min-height` so a missing or slow SVG does not collapse the box to zero
       and let the 3D viewer jump up into the facts grid, which is the one
       mechanism that produces exactly the reported overlap. */
#vstructwrap{border:1px solid var(--rule);border-radius:4px;background:#fff;
 padding:6px;margin-bottom:10px;display:flex;justify-content:center;
 overflow:hidden;min-height:96px;align-items:center}
#vstruct{display:block;width:100%;max-width:420px;height:auto;max-height:300px;
 object-fit:contain}
.facts{align-items:start}
#glbox{position:relative;width:100%;height:440px;background:#eef1f6;
 border:1px solid var(--rule);border-radius:4px;overflow:hidden}
#glbox>div{position:absolute;inset:0}
#glbox canvas{position:absolute;top:0;left:0}
.vctl{display:flex;flex-wrap:wrap;gap:.4rem 1.2rem;padding:.6rem .1rem 0;font-size:12px}
.vctl label{display:flex;align-items:center;gap:.35rem;cursor:pointer;font-weight:600}
.note{font-size:12px;color:var(--muted);margin:.9rem 0 0;max-width:78ch}
#sibs h3{font:600 11px var(--sans);letter-spacing:.05em;text-transform:uppercase;
 color:var(--muted);margin:1.3rem 0 .4rem}
table.sib{border-collapse:collapse;width:100%;font-size:12.5px}
table.sib th,table.sib td{padding:.34rem .6rem;border-bottom:1px solid var(--rule);
 text-align:right;white-space:nowrap}
table.sib th{font:600 10px var(--sans);color:var(--muted);text-transform:uppercase;
 letter-spacing:.05em;border-bottom:2px solid var(--rule)}
table.sib th:first-child,table.sib td:first-child,
table.sib th:last-child,table.sib td:last-child{text-align:left}
table.sib td{font-family:var(--mono)}
tr.sibrow{cursor:pointer}
tr.sibrow:hover{background:var(--blue-pale)}
tr.sibrow.cur{background:var(--blue-pale);font-weight:700}
input.mchk{margin:0 .45rem 0 0;vertical-align:-1px;cursor:pointer}
/* The sub-split GROUP HEADER and its indent are gone (#65). They boxed a
   molecule's modes into one block per first-stage cluster, each with its own
   all/none, which read as a distinction between kinds of mode -- and it is not
   one: both clusters are first-stage clusters, and which one a mode came out of
   is provenance, not chemistry. The provenance survives as the `cluster` column
   and the row's hue. */
button.mini{font:600 10px var(--sans);margin-left:.5rem;padding:.1rem .42rem;
  border:1px solid var(--rule);border-radius:3px;background:var(--bg);
  color:var(--fg);cursor:pointer;text-transform:uppercase}
button.mini:hover{background:var(--blue-pale)}
span.gp{font:400 10px var(--mono);color:var(--muted);
  text-transform:none;letter-spacing:0}
i.sw{width:11px;height:11px;border-radius:2px;display:inline-block;margin-right:.45rem;
 vertical-align:-1px}
.win{color:var(--good);font-weight:700}
.na{color:var(--muted);font-style:italic}
.warnbox{border-left:3px solid var(--warn);background:#fdf8ea;padding:.6rem .9rem;
 font-size:12px;margin:0 0 12px;border-radius:0 3px 3px 0}
:root[data-theme="dark"] .warnbox{background:#241f12}
:root[data-theme="dark"] .thumb{background:#fff}
a{color:var(--blue)}
__STEPCSS__
</style></head><body>
<div id="topbar">
 <h1 title="Pick a mode on the left; its pose and scores load on the right.">__TITLE__ — ranking</h1>
 <span class="msep"></span>
 <select id="scope" onchange="setScope(this.value)"
   title="rank within one warhead class, within every class, or across all of them"></select>
 <span class="mhint" id="mhint"></span>
 <span class="msep"></span>
 <a class="mbtn lnk" href="pipeline.html" title="how a molecule becomes a row">how this works &#8599;</a>
 <button class="mbtn" onclick="toggleTheme()">dark</button>
</div>
<!-- THE SAME STEPPER AS EVERY OTHER PAGE (#63). The ad-hoc "results" link that
     sat in the topbar is gone: two routes to one page, only one of which knows
     where you are, is what made these read as separate tools. -->
__STEPNAV__
<main>
 <!-- ONE GRID CHILD. `main` is a two-column grid, so the banner must live INSIDE
      the rail's column rather than beside it -- as a third child it took column 1
      and pushed the viewer onto row 2. railPad carries the FULL scroll height so
      the scrollbar stays honest about how much library there is; railWin holds
      only the ~40 items actually on screen. -->
 <div id="railcol">
  __RAILSEARCH__
  <div class="railbn" id="railBanner" style="display:none"></div>
  <div id="rail"><div id="railPad"><div id="railWin"></div></div></div>
 </div>
 <div id="viewer">
  <div id="vhead"><span id="vname">select a mode</span>
   <span class="mhint">medoid pose, in the receptor it was docked into</span></div>
  <div id="vbody">
   <div id="vempty" class="note">This is the <strong>ranking</strong> view: every
   molecule and every mode the screen scored, simulated or not. In the pipeline's
   real order it comes first — you read the ranked list and look at the poses, then
   choose what goes to the triage sweep. It was built after the fact, which is
   how <a href="https://github.com/hallettmiket/inhibition/issues/53">#53</a> went
   unnoticed: nothing in the project had ever shown the per-mode ranking it
   computes.</div>
   <div id="vfull" style="display:none">
    <div id="gwarn" class="warnbox" style="display:none"></div>
    <div id="vmiss" class="warnbox" style="display:none"></div>
    <div class="facts" id="vfacts"></div>
    <div id="vstructwrap"><img id="vstruct" alt="2D structure"></div>
    <div id="glbox"><div id="gl"></div></div>
    <div class="vctl">
     <label><input type="checkbox" id="c-surf" checked> pocket surface</label>
    </div>
    <div id="sibs"></div>
    <p class="note" id="vnote"></p>
   </div>
  </div>
 </div>
</main>
<script type="text/plain" id="recpdb">__RECEPTOR__</script>
<script>
const ROWS = __ROWS__;
// `i` is rebuilt here rather than sent 34,076 times. Everything downstream keeps
// using `x.i`, so this is a wire-format saving and not a change of shape.
ROWS.forEach(function(x){ x.i = x.p + '_m' + x.m; });
// Molecules with NO representative pose from the production run. Their stored
// asset is from an earlier screen with a different clustering, so its models do
// not mean what this table's modes mean and NONE of them may be drawn.
const NOPOSE = new Set(__NOPOSE__);
const MODE_COLS = [0x0072ce, 0x7b5ea7, 0xc2703d, 0x0f7a54, 0xb3261e, 0x8a6d1f];
const MODE_CSS  = ['#0072ce','#7b5ea7','#c2703d','#0f7a54','#b3261e','#8a6d1f'];
// HUE BY FIRST-STAGE CLUSTER, LIGHTNESS BY POSITION WITHIN IT (#61).
//
// The colour records PROVENANCE -- which first-stage cluster a mode was cut out
// of -- so a reader can draw one cluster's rows together and see how they sit
// relative to each other. That is a genuinely useful comparison and it is why
// the grouping survives.
//
// IT IS NOT A STATEMENT THAT THEY ARE ALIKE. This comment used to say the
// second stage produces "alternative scaffold placements of one warhead
// geometry", i.e. rows differing only away from the reactive end. Measured on
// the 3.0.0 run that is false for a large minority: the first stage is a
// chaining cluster (DBSCAN, eps 3 A with 2 A per radian, so poses in one place
// join at up to 86 degrees apart and chains extend indefinitely), and 22% of
// split clusters hold modes whose median reactive-atom distance spans more than
// the criterion's entire 2.8-4.2 A window. Same hue, opposite sides of the
// window. The rows are therefore NAMED by their own mode index and the shared
// hue is labelled as provenance in the group header.
const MODE_HUES = [205, 268, 25, 158, 4, 45];
function subIx(m){
  // `ml` is '3' when a mode was never subdivided and '3b' when it was.
  const s = String((m && m.ml) || '');
  const c = s.charCodeAt(s.length - 1);
  return (c >= 97 && c <= 122) ? c - 97 : 0;
}
function modeHsl(m){
  const h = MODE_HUES[Math.abs((m && m.pm) || 0) % MODE_HUES.length];
  return [h, 54, 32 + Math.min(subIx(m), 4) * 9];
}
function modeCss(m){ const c = modeHsl(m); return 'hsl(' + c[0] + ',' + c[1] + '%,' + c[2] + '%)'; }
function modeHex(m){
  // 3Dmol wants an integer colour, so the same HSL is converted rather than a
  // second palette being kept in parallel -- two palettes is how the rail and
  // the viewer come to disagree about which mode is which colour.
  const c = modeHsl(m), s = c[1] / 100, l = c[2] / 100;
  const k = n => (n + c[0] / 30) % 12;
  const f = n => l - s * Math.min(l, 1 - l) *
                 Math.max(-1, Math.min(Math.min(k(n) - 3, 9 - k(n)), 1));
  return (Math.round(f(0) * 255) << 16) | (Math.round(f(8) * 255) << 8)
         | Math.round(f(4) * 255);
}
const POCKET = __POCKET__;
// CARBONS CARRY THE MODE COLOUR; EVERY OTHER ELEMENT KEEPS ITS CONVENTIONAL ONE.
// Colouring a whole molecule by mode hides its chemistry -- the sulfur, the
// halogen and the oxygens are what a chemist reads a pose by, and they must look
// the same in every mode so the only thing that changes is the carbon skeleton.
function carbonScheme(col){
  const M = lib();
  return {prop:'elem',
          map: Object.assign({}, (M.elementColors||{}).defaultColors||{}, {C: col})};
}
// SCOPE is one control: a warhead class name, '*' for every class ranked within
// itself, or '__global__' for one order across all of them. Two orthogonal
// toggles let a reader combine "global" with a class filter, which is a
// combination with no meaning.
let SCOPE = '*', SEL = null, V = null, SURF = null;
// SHOWN holds the modes currently drawn for the selected molecule, so several
// alternatives can be compared and any of them switched off again. SEL stays the
// PRIMARY -- the one the facts panel describes -- because a panel of numbers has
// to be about one mode, and "several are visible" is a different question from
// "which one am I reading".
let SHOWN = new Set(), PDBCACHE = {};
// mode number -> that mode's row, for the molecule on screen. The 3D layer only
// ever knows a mode NUMBER (it comes off the PDB models), but colour and label
// now depend on the row's `pm`/`ml`, so the lookup has to exist somewhere. One
// map, rebuilt in pick(), rather than a scan per model per redraw.
let MBY = {};

function lib(){ return window.$3Dmol || window['3Dmol']; }
function fmt(x, d){ return (x === null || x === undefined) ? '—' : (+x).toFixed(d); }

function isGlobal(){ return SCOPE === '__global__'; }
document.addEventListener('keydown', function(e){
  const el = document.getElementById('railq');
  if (!el) return;
  if (e.key === '/' && document.activeElement !== el){ e.preventDefault(); el.focus(); }
  else if (e.key === 'Escape' && document.activeElement === el){
    el.value = ''; setQuery(''); }
});

let QUERY = '';
// The shared box's clear button calls this by name.
function railClear(){
  const el = document.getElementById('railq');
  if (el){ el.value = ''; setQuery(''); el.focus(); }
}
function setQuery(q){
  QUERY = (q || '').trim().toLowerCase();
  // Rebuild rather than hide: THE RAIL IS VIRTUALISED, so only the rows
  // currently on screen exist in the DOM. Hiding those would filter the window
  // and leave everything below it unfiltered as soon as you scrolled.
  railHTML();
  const c = document.getElementById('mhint');
  if (c && QUERY) c.textContent = visible().length.toLocaleString()
    + ' of ' + ROWS.length.toLocaleString() + ' match';
  else if (c) c.textContent = ROWS.length.toLocaleString() + ' modes';
}
function matches(x){
  if (!QUERY) return true;
  // The ident is not stored per row (it is rebuilt from p + "_m" + m to save a
  // megabyte of payload), so it is rebuilt here too rather than assumed absent.
  const id = x.p + '_m' + x.m;
  return (id + ' ' + x.c + ' ' + (x.ml || '')).toLowerCase().indexOf(QUERY) !== -1;
}
function visible(){
  let r = ROWS.slice();
  if (SCOPE !== '*' && !isGlobal()) r = r.filter(x => x.c === SCOPE);
  if (QUERY) r = r.filter(matches);
  const K = x => (x === null || x === undefined) ? 1e9 : x;   // unranked sorts last
  if (isGlobal()) r.sort((a,b) => K(a.gr) - K(b.gr));
  else r.sort((a,b) => a.c.localeCompare(b.c) || K(a.cr) - K(b.cr));
  return r;
}

function buildScope(){
  const n = {};
  ROWS.forEach(x => { n[x.c] = (n[x.c] || 0) + 1; });
  const opts = ['<optgroup label="ranked within its own class">',
    '<option value="*">all classes</option>'];
  Object.keys(n).sort().forEach(c =>
    opts.push('<option value="' + c + '">' + c + ' (' + n[c].toLocaleString() + ')</option>'));
  opts.push('</optgroup><optgroup label="across classes">',
    '<option value="__global__">global — biased (#47)</option></optgroup>');
  const el = document.getElementById('scope');
  el.innerHTML = opts.join('');
  el.value = SCOPE;
}

// ---------------------------------------------------------------------------
// THE RAIL IS VIRTUALISED. Only what fits on screen exists in the DOM.
//
// It used to render every visible mode as a <button> holding an <img> and eight
// <span>s, in one innerHTML assignment. At 8,097 modes that was slow and
// tolerable. Sub-splitting (#61) took the library to 34,076, which is roughly
// 350,000 nodes and 34,000 lazy images -- and `pick()` called it, so EVERY row
// click and EVERY mode checkbox rebuilt the whole thing. The page became
// unusable, and it became unusable because the view scaled with the library
// instead of with the viewport.
//
// Now: one flat ITEMS list (headers and rows), a cumulative offset table, and a
// window of ~40 items rendered on scroll. Cost is constant in the library size.
// Selection no longer rebuilds anything -- it moves one CSS class.
// ONE CLASS LABEL, NOT TWO. The pre-virtualisation rail used a `position:sticky`
// class header: one element that both marked the boundary and stayed on screen as
// you scrolled past it. Sticky cannot survive a transformed window, so it was
// replaced by a banner above the rail -- but the inline header was left in place
// too, and the two said the same word one under the other (@tt8804: "says
// acrylamide twice"). The banner is strictly the better of the two: it is
// readable at every scroll position, not just at the boundary. So the inline
// headers are gone and every item in the list is a row.
let ITEMS = [], OFF = [], TOTAL = 0, VIS = [];
const ROW_H = 64, OVERSCAN = 8;

function rowHTML(x){
  const rank = isGlobal() ? (x.gr === null ? '—' : x.gr)
                         : (x.cr === null ? '—' : x.cr);
  const pct = x.vf === null ? 0 : Math.round(x.vf * 100);
  const badge = x.ctl ? 'known inhibitor'
              : x.s === 'md' ? '100 ns' : x.s === 'swept' ? 'swept'
              : x.s === 'failed' ? 'failed' : 'not run';
  return '<button class="row' + (SEL === x.i ? ' on' : '') + (x.ctl ? ' ctl' : '') +
    '" data-i="' + x.i + '" onclick="pick(\'' + x.i + '\')">' +
    '<span class="rk">' + rank + '</span>' +
    '<img class="thumb" loading="lazy" alt="" src="mode_thumbs/' + x.p + '.svg">' +
    '<span class="body"><span class="l1">' +
    '<span class="mid-id">' + x.p +
    // THE MODE IS NAMED BY ITS OWN INDEX, which is its identity everywhere else
    // -- `t4_x_m3` in the rank table, the sweep table and the pose file. This
    // printed the `0d` letter label for a while, on the reasoning that the
    // letter says which first-stage mode a row came from. It does, and that is
    // the problem: a lettered name reads as a variant of `0a`, and sub-modes of
    // one first-stage mode are NOT variants of each other. Measured on this
    // run, 22% of split first-stage modes hold sub-modes whose median
    // reactive-atom distance spans more than the criterion's entire 2.8-4.2 A
    // window, and 18% have some sub-modes inside the window and some outside.
    // The first-stage origin is still shown -- as provenance, in the group
    // header and the detail caption, not as the row's name.
    ' <span class="mtag" style="background:' + modeCss(x) + '">m' + x.m +
    '</span></span>' +
    '<span class="eng">' + fmt(x.eb, 2) + '</span></span>' +
    '<span class="l2"><span class="wc">' + x.c + '</span>' +
    '<span class="meta">' + (x.n === null ? '—' : x.n) + ' poses · ' + pct + '% viable</span>' +
    '<span class="tag ' + (x.ctl ? 't-ctl' : 't-' + x.s) + '">' + badge +
    '</span></span>' +
    '<span class="bar"><i style="width:' + pct + '%"></i></span></span></button>';
}

function buildItems(){
  // The sort is here and ONLY here. `visible()` slices and sorts 34,076 rows;
  // doing that inside the render made every scroll frame pay for it.
  VIS = visible();
  ITEMS = VIS;
  // Every item is one row of ROW_H, so the offset of item i IS i * ROW_H and the
  // window is arithmetic rather than a search. The cumulative table is kept
  // because it is what makes the geometry explicit and checkable.
  OFF = new Array(ITEMS.length + 1);
  for (let i = 0; i <= ITEMS.length; i++) OFF[i] = i * ROW_H;
  TOTAL = ITEMS.length * ROW_H;
  document.getElementById('railPad').style.height = TOTAL + 'px';
  document.getElementById('mhint').textContent =
    VIS.length.toLocaleString() + ' modes' +
    (isGlobal() ? ' · one order across classes scored under different bars (#47)'
                : ' · rank is within the warhead class');
}

function renderRail(){
  const el = document.getElementById('rail');
  const top = el.scrollTop, h = el.clientHeight || 600;
  // Binary search rather than a scan: at 34,076 items a linear search per scroll
  // frame is the very cost this is removing.
  let lo = 0, hi = ITEMS.length - 1, s = 0;
  while (lo <= hi){ const m = (lo + hi) >> 1;
    if (OFF[m + 1] <= top) lo = m + 1; else { s = m; hi = m - 1; } }
  let e = s;
  while (e < ITEMS.length && OFF[e] < top + h) e++;
  s = Math.max(0, s - OVERSCAN); e = Math.min(ITEMS.length, e + OVERSCAN);
  const out = [];
  for (let i = s; i < e; i++) out.push(rowHTML(ITEMS[i]));
  const win = document.getElementById('railWin');
  win.style.transform = 'translateY(' + (ITEMS.length ? OFF[s] : 0) + 'px)';
  win.innerHTML = out.join('');
  // THE ONLY CLASS LABEL IN THE RAIL. It names the class of the row at the top of
  // the viewport, so it both marks the boundary (the text changes as you cross
  // one) and answers "where am I" deep inside a class -- which is what the old
  // sticky header did and what an inline header alone cannot.
  const first = ITEMS[Math.min(Math.max(s, 0) + (s > 0 ? OVERSCAN : 0),
                               ITEMS.length - 1)];
  const bn = document.getElementById('railBanner');
  if (first){
    bn.textContent = isGlobal()
      ? 'global order — biased across classes (#47)'
      : first.c + ' — ranked within the class';
    bn.style.display = '';
  } else { bn.style.display = 'none'; }
}

function railHTML(){ buildItems(); renderRail(); }

// Selection moves a class; it does not rebuild the list. This is the other half
// of the fix -- virtualising the render is wasted if every click re-renders.
function markSel(){
  const win = document.getElementById('railWin');
  win.querySelectorAll('button.row.on').forEach(b => b.classList.remove('on'));
  const b = win.querySelector('button.row[data-i="' + SEL + '"]');
  if (b) b.classList.add('on');
}

function setScope(v){ SCOPE = v;
  const el = document.getElementById('rail'); if (el) el.scrollTop = 0;
  railHTML(); }

function showModes(on){
  // Draw or clear EVERY mode of the selected molecule at once.
  //
  // This used to be one control per first-stage cluster, so a molecule with two
  // clusters got two "all" buttons that each drew part of the cloud. That split
  // the one interaction anybody wants -- see this molecule's modes together --
  // along a boundary that is provenance rather than chemistry, and @tt8804 read
  // the two groups as a claim that the clusters are different kinds of thing:
  // "get rid of the separated selection modes, both are first stage". They are.
  const cur = SEL ? ROWS.find(r => r.i === SEL) : null;
  if (!cur) return;
  ROWS.filter(r => r.p === cur.p).forEach(function(m){
    if (on) SHOWN.add(m.m);
    else if (m.m !== cur.m) SHOWN.delete(m.m);   // the primary always stays drawn
  });
  pick(SEL);
}

function toggleMode(m){
  // The primary cannot be hidden -- the facts panel is describing it, and a
  // panel of numbers with nothing drawn beside it reads as a rendering failure.
  const cur = SEL ? ROWS.find(r => r.i === SEL) : null;
  if (cur && m === cur.m) return;
  if (SHOWN.has(m)) SHOWN.delete(m); else SHOWN.add(m);
  if (SEL) pick(SEL);
}
function toggleTheme(){
  const d = document.documentElement.getAttribute('data-theme') === 'dark';
  document.documentElement.setAttribute('data-theme', d ? 'light' : 'dark'); }

async function pick(id){
  const x = ROWS.find(r => r.i === id); if (!x) return;
  const prev = SEL ? ROWS.find(r => r.i === SEL) : null;
  if (!prev || prev.p !== x.p) SHOWN = new Set();   // new molecule, new selection
  SHOWN.add(x.m);                                    // the primary is always drawn
  SEL = id; markSel();
  document.getElementById('vempty').style.display = 'none';
  document.getElementById('vfull').style.display = '';
  // The ident names the mode; the first-stage cluster is stated after it as
  // PROVENANCE. "sub-mode 0c of first-stage mode 0" was the old wording and it
  // asserted a similarity the geometry does not support -- see the rail.
  document.getElementById('vname').textContent =
    x.i + (x.ml !== String(x.m)
           ? '   (from first-stage cluster ' + x.pm + ', ranked on its own)' : '');
  // The same depiction the rail uses, at panel size. It is an SVG, so one file
  // serves both; drawing a second at a larger size would be a second answer to
  // "what does this molecule look like".
  document.getElementById('vstruct').src = 'mode_thumbs/' + x.p + '.svg';
  document.getElementById('vfacts').innerHTML = [
    ['class rank', x.cr === null ? 'unranked' : x.cr],
    ['global rank', x.gr === null ? '—' : x.gr],
    ['poses in mode', x.n === null ? '—' : x.n + ' of ' + (x.np === null ? '?' : x.np)],
    ['viable fraction', x.vf === null ? '—' : (x.vf*100).toFixed(1) + '%'],
    ['enrichment', fmt(x.en, 2)],
    ['conditional_eb', fmt(x.eb, 3)],
    ['spread', fmt(x.sp, 2) + ' Å'],
    ['direction coherence', fmt(x.dc, 3)],
    ['__SWEEPLABEL__ sweep', x.fa === null ? (x.s === 'none' ? 'never sent' : 'no score')
                                  : (x.fa*100).toFixed(1) + '% ready'],
  ].map(kv => '<div class="fact"><b>' + kv[1] + '</b><span>' + kv[0] + '</span></div>').join('');

  const g = document.getElementById('gwarn');
  if (x.ctl){
    g.style.display = '';
    g.innerHTML = '<strong>Known inhibitor, run as a candidate.</strong> Sulfopin '
      + '(Dubiella 2021, 6VAJ) sits in the candidate frame and went through '
      + 'docking, mode splitting, ranking and the sweep by the identical path as '
      + 'every other row, scored on the same column. The highlight is '
      + 'presentational: nothing here was computed differently because we know '
      + 'the answer.';
  } else if (x.s === 'none'){
    g.style.display = '';
    g.innerHTML = '<strong>Never simulated.</strong> This mode was scored and ranked, '
      + 'and no sweep or MD was ever run from it. Every number above is docking-derived.';
  } else { g.style.display = 'none'; }

  // EVERY MODE OF THIS MOLECULE, and how they were ranked against each other.
  // The rail orders modes across the whole library; this is the comparison the
  // pipeline claims to make -- a molecule's modes competing as separate rows --
  // and it is the one place a reader can see whether the mode that was simulated
  // is the one that scored best.
  const sibs = ROWS.filter(r => r.p === x.p).sort((a,b) => a.m - b.m);
  MBY = {}; sibs.forEach(function(m){ MBY[m.m] = m; });
  const RK = m => (m.cr === null || m.cr === undefined) ? 1e9 : m.cr;
  const best = sibs.reduce((a,b) => (RK(b) < RK(a) ? b : a), sibs[0]);
  // ONE FLAT LIST. Every row is a mode; nothing is nested under anything.
  //
  // The rows used to be boxed under a "first-stage cluster N" header with its
  // own all/none pair, which put a molecule's modes into two selection groups
  // and made the grouping look like a distinction between kinds of mode. It is
  // not one: both clusters are first-stage clusters, and the split between them
  // is provenance -- which connected component of a CHAINING clusterer a mode
  // was cut out of -- not chemistry. @tt8804: "get rid of the separated
  // selection modes, both are first stage."
  //
  // The provenance is kept, demoted to a column and the row's hue, so a reader
  // who wants it can still see which modes came out of one cluster.
  const nClust = new Set(sibs.map(m => m.pm)).size;
  document.getElementById('sibs').innerHTML =
    '<h3>modes of ' + x.p + ' — ' + sibs.length + ' in total, '
    + sibs.filter(m => m.cr !== null).length + ' ranked'
    + ' · every row is its own mode, ranked and swept independently (#61)</h3>' +
    '<p class="note" style="margin:.2rem 0 .5rem">Tick to draw a mode; click the '
    + 'row to read it. Several can be shown at once. '
    + '<button class="mini" onclick="showModes(1)">all</button>'
    + '<button class="mini" onclick="showModes(0)">none</button>'
    + (nClust < sibs.length
       ? ' <b>cluster</b> is the first-stage group a mode was cut from, and is '
         + 'provenance only — the first stage chains, so one cluster can hold '
         + 'modes on opposite sides of the 2.8–4.2 Å window (#65).'
       : '') + '</p>'
    + '<table class="sib"><thead><tr><th>show</th><th>class rank</th><th>poses</th>' +
    '<th>viable</th><th>enrichment</th><th>conditional_eb</th><th>spread</th>' +
    '<th>coherence</th><th>cluster</th><th>simulated</th></tr></thead><tbody>' +
    sibs.map(function(m){
      const col = modeCss(m);
      const badge = m.s === 'md' ? '100 ns' : m.s === 'swept' ? 'swept'
                  : m.s === 'failed' ? 'sweep failed' : 'never';
      return '<tr class="sibrow' + (m.i === x.i ? ' cur' : '') + '"'
        + ' onclick="pick(\'' + m.i + '\')">'
        + '<td onclick="event.stopPropagation();toggleMode(' + m.m + ')">'
        + '<input type="checkbox" class="mchk"' + (SHOWN.has(m.m) ? ' checked' : '')
        + ' onclick="event.stopPropagation();toggleMode(' + m.m + ')">'
        + '<i class="sw" style="background:' + col + '"></i>m' + m.m + '</td>'
        + '<td' + (m.i === best.i ? ' class="win"' : '') + '>'
        + (m.cr === null ? '<span class="na">unranked</span>' : m.cr) + '</td>'
        + '<td>' + (m.n === null ? '—' : m.n) + '</td>'
        + '<td>' + (m.vf === null ? '—' : (m.vf*100).toFixed(1) + '%') + '</td>'
        + '<td>' + fmt(m.en, 2) + '</td><td>' + fmt(m.eb, 3) + '</td>'
        + '<td>' + fmt(m.sp, 2) + '</td><td>' + fmt(m.dc, 3) + '</td>'
        + '<td><span class="gp">' + m.pm + '</span></td>'
        + '<td><span class="tag t-' + m.s + '">' + badge + '</span></td></tr>';
    }).join('') + '</tbody></table>' +
    (sibs.length > 1 && best.s === 'none'
      ? '<p class="note"><strong>The best-ranked mode of this molecule was never '
        + 'simulated.</strong> m' + best.m + ' ranks ' + best.cr + ' in '
        + best.c + '; the sweep took mode 0 (#53).</p>'
      : '');

  document.getElementById('vnote').innerHTML =
    'The pose is this mode\'s <strong>medoid</strong> — the pose most central to '
    + 'the mode among its best-anchored quartile, not its lowest-energy member. The '
    + 'individual poses were not persisted '
    + '(<a href="https://github.com/hallettmiket/inhibition/issues/44">#44</a>).';

  // REFUSED, NOT WARNED ABOUT. The stored asset for these molecules belongs to
  // an earlier run whose mode numbering is unrelated to this table's, so drawing
  // model 2 for mode 2 would show a pose from a different clustering of a
  // different docking run and look entirely normal.
  if (NOPOSE.has(x.p)){
    const miss = document.getElementById('vmiss');
    miss.style.display = '';
    miss.innerHTML = '<strong>No pose from this run — nothing is drawn.</strong> '
      + x.p + ' was docked and scored by the current screen (the numbers above '
      + 'are its own), but its pose files were never written: the screen aimed '
      + 'them at the previous run\'s directory, where a file already existed, and '
      + 'the append-only guard skipped the write. The asset still on disk is the '
      + 'PREVIOUS run\'s, under a different mode numbering, so it is refused '
      + 'rather than shown. This molecule needs re-docking.';
    document.getElementById('gl').innerHTML =
      '<p class="note" style="padding:14px">no pose from this run</p>';
    return;
  }
  try {
    if (!PDBCACHE[x.p]){
      const res = await fetch('mode_poses/' + x.p + '.pdb');
      if (!res.ok) throw new Error(res.status);
      PDBCACHE[x.p] = await res.text();
    }
    draw(PDBCACHE[x.p], x);
  } catch (e) {
    // THE REASON, NOT A GUESS AT IT. This said "no pose file for X" for every
    // failure, including exceptions thrown inside draw() -- so a rendering bug
    // was reported as a missing file and looked like a data problem. The message
    // now carries what actually went wrong.
    document.getElementById('gl').innerHTML =
      '<p class="note" style="padding:14px"><strong>No pose drawn.</strong><br>'
      + String(e && e.message ? e.message : e) + '<br><span class="na">'
      + 'asset: mode_poses/' + x.p + '.pdb</span></p>';
  }
}

function draw(pdbTxt, x){
  const M = lib(); if (!M) return;
  if (!V) V = M.createViewer(document.getElementById('gl'), {backgroundColor:'#eef1f6'});
  V.clear(); SURF = null;
  V.addModel(document.getElementById('recpdb').textContent, 'pdb');
  // Each block carries its own mode in the MODEL record, so the model->mode
  // mapping is READ rather than counted. Counting positions is #53.
  const blocks = pdbTxt.split('ENDMDL').filter(b => b.indexOf('MODEL') >= 0);
  const modes = [];
  blocks.forEach(b => {
    const m = /MODEL\s+(-?\d+)/.exec(b);
    modes.push(m ? parseInt(m[1], 10) : -1);
    V.addModel(b.replace(/MODEL[^\n]*\n/, ''), 'pdb');
  });
  V.setStyle({}, {cartoon:{color:'#c3ccd8', opacity:0.5}});
  // CYS113 IN CONVENTIONAL ELEMENT COLOURS, so its sulfur reads as sulfur --
  // it is the atom the whole screen is aimed at. SG additionally gets a sphere,
  // because at stick radius a single S is easy to lose against the cartoon.
  V.setStyle({resi:[__CYS__]},
             {stick:{radius:0.28, colorscheme:'default'},
              cartoon:{color:'#c3ccd8', opacity:0.5}});
  V.addStyle({resi:[__CYS__], atom:'SG'}, {sphere:{radius:0.62}});
  // Draw every mode the reader has ticked. The PRIMARY is thicker and fully
  // opaque so it stays identifiable in a stack of alternatives; the rest are
  // thinner and translucent, which is the difference between comparing poses and
  // producing one unreadable object out of several.
  modes.forEach(function(mo, i){
    if (!SHOWN.has(mo)){ V.setStyle({model: i+1}, {}); return; }
    const primary = (mo === x.m);
    V.setStyle({model: i+1}, {stick:{
      radius: primary ? 0.22 : 0.14,
      opacity: primary ? 1 : 0.6,
      colorscheme: carbonScheme(modeHex(MBY[mo]))}});
  });
  if (document.getElementById('c-surf').checked){
    // The pocket shell only, and never over Cys113 or a ligand: a mesh drawn on
    // top of those hides the two things the panel exists to show. `and` rather
    // than a bare `not`, so the selection cannot leak onto the pose models.
    SURF = V.addSurface(M.SurfaceType.VDW, {opacity:0.62, color:'#b9c7db'},
      {and:[{model:0}, {resi:POCKET}, {not:{resi:[__CYS__]}}]});
  }
  const sel = modes.indexOf(x.m);
  V.zoomTo(sel >= 0 ? {model: sel+1} : {resn:'MOL'});
  V.zoom(0.5); V.resize();
  V.render();                       // 3Dmol draws NOTHING without this.
  // THE ASSET MAY NOT HOLD THIS MODE, AND SILENCE LOOKS IDENTICAL TO A BUG.
  // 195 molecules carry an asset from an earlier run with a single pose in it
  // while the ranking lists five modes; asking for model 3 then selects nothing
  // and the panel renders an empty box. Say which it is.
  const miss = document.getElementById('vmiss');
  if (sel < 0){
    miss.style.display = '';
    miss.innerHTML = '<strong>This mode has no pose in the stored asset.</strong> '
      + 'The asset for ' + x.p + ' holds ' + modes.length + ' pose(s) ('
      + modes.map(q => 'm' + q).join(', ') + ') and this row is m' + x.m
      + '. It predates the current screen, so the molecule needs re-docking '
      + 'before its poses can be shown. The numbers above are unaffected — they '
      + 'come from the tables, not from this file.';
  } else { miss.style.display = 'none'; }
}

document.getElementById('c-surf').addEventListener('change', function(){ if (SEL) pick(SEL); });
// One render per animation frame at most. A scroll event can fire far more often
// than the screen refreshes, and rendering per event is how a virtualised list
// ends up no faster than the thing it replaced.
(function(){
  const el = document.getElementById('rail');
  let queued = false;
  el.addEventListener('scroll', function(){
    if (queued) return;
    queued = true;
    requestAnimationFrame(function(){ queued = false; renderRail(); });
  }, {passive: true});
  window.addEventListener('resize', function(){ renderRail(); });
})();
buildScope();
railHTML();
// The rail's clientHeight is 0 until layout settles, so the first window would
// be sized from the fallback and could come up short. Re-render once the browser
// has actually laid the page out -- the same double-rAF the 3D panel needs.
requestAnimationFrame(function(){ requestAnimationFrame(renderRail); });
</script>
</body></html>"""


def build(title: str, date_str: str, three: str = "",
          no_pose: list | None = None) -> str:
    """`no_pose` names molecules with NO representative from the production run.

    THEIR ASSET MUST NOT BE DRAWN AT ALL, and a warning is not enough. 196
    molecules kept an asset from 2.2.0 holding that run's five modes while the
    3.0.0 table lists twelve. Asking for mode 7 drew nothing and was reported;
    asking for mode 2 drew 2.2.0's mode 2 -- a pose from a DIFFERENT clustering
    of a different docking run -- silently, beside 3.0.0's numbers. That is the
    exact defect this whole release exists to remove, reappearing one layer out.

    So the page is given the list and refuses those molecules outright. Refusing
    is not a limitation to work around: there is no correct pose to show until
    the molecule is re-docked, and showing a plausible wrong one is worse than
    showing none.
    """
    r = gather()
    if r.empty:
        return "<!doctype html><p>no rank tables found</p>"
    rec = ("\n".join(l for l in RECEPTOR.read_text().splitlines()
                     if l.startswith(("ATOM", "HETATM")))
           if RECEPTOR.is_file() else "")
    return (_TPL
            .replace("__ROWS__", _rows_json(r))
            .replace("__NOPOSE__", json.dumps(sorted(no_pose or [])))
            .replace("__RECEPTOR__", rec)
            .replace("__THREE__", three)
            .replace("__CYS__", str(CYS_RESI))
            .replace("__POCKET__", json.dumps(pocket_residues()))
            .replace("__ROWCSS__", _rs().ROW_CSS)
            .replace("__RAILQCSS__", _rs().SEARCH_CSS)
            .replace("__RAILSEARCH__",
                     _rs().search_html("setQuery(this.value)",
                                       "filter — id, class, mode"))
            .replace("__SWEEPLABEL__", _gs().sweep_label())
            .replace("__STEPCSS__", _gs().CSS)
            .replace("__STEPNAV__", _gs().nav("modes.html", _step_counts(r)))
            .replace("__TITLE__", html.escape(title)))


def _rp():
    from shared import run_paths
    return run_paths


def _gs():
    from shared import gui_shell
    return gui_shell


def _rs():
    """The shared rail-row stylesheet. Imported lazily for the same reason
    `_gs` is: this module is imported by the builders that own it."""
    from shared import results_shell
    return results_shell


def _step_counts(r) -> dict:
    """Delegates to `gui_shell.step_counts()` -- the one source.

    THIS IS WHERE THE 447 CAME FROM. It computed its own counts and read
    `mdprio_reports/sweep_state.json`: the UNSCOPED path, so a page built for
    nac_v5 showed 3.0.0's 447 swept modes in its nav while the Sweep page next
    to it showed this run's 34. Both files existed, so nothing failed.

    Kept as a shim rather than deleted because two call sites use it, one of them
    inside a template replacement.
    """
    return _gs().step_counts()
