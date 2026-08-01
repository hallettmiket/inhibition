"""
Purpose: The integration GUI — present the four shortlists, do not auto-rank them.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-28
Input: the append-only frames, run manifests and decision records
Output: a Streamlit app

Run:  /data/lab_vm/envs/dwi_gui/bin/streamlit run integration/app/app.py

PRESENT, DON'T AUTO-RANK. Each approach delivers a shortlist ranked by its own
internally-valid metric. `vina_affinity` and `affinity_kcal` are different
quantities produced under different protocols, so there is no defensible way to
sort all four together, and this app never does. It shows them side by side and
lets a human decide.

EVERY RANK CARRIES ITS GATE. The non-covalent gate has issued a verdict —
WEAK, ROC-AUC 0.599, CI [0.311, 0.874], EF1% 0.0, over 6 actives and 6
independent chemotypes (D0041). It is the first non-UNDERPOWERED verdict in the
project and it says docking does not demonstrably enrich for known binders. The
covalent gate is still UNDERPOWERED. So no shortlist here is evidence that the
molecules at the top bind, and a rank displayed without that verdict would imply
a confidence nothing supports.

THE RANKINGS ARE PARTLY SIZE RANKINGS (D0043). spearman(heavy_atoms, rank
metric) is -0.617 for T_1, -0.479 for T_3 and -0.230 for T_2, all
lower-is-better, so larger molecules score better in three of four approaches.
T_4 is the exception at +0.181. That is a mechanism for the WEAK verdict, since
the decoys are matched on molecular weight. Ligand efficiency is not the fix —
it over-corrects to -0.938 and becomes a smallness ranking.

THE SCORE-FREE SIGNALS ARE THE DEFENSIBLE ONES. Structural convergence and the
shared physicochemical axes need no commensurability between metrics, which is
why they get their own panels while the cross-approach "leaderboard" does not
exist.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

APP_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(APP_DIR))
sys.path.insert(0, str(APP_DIR.parent.parent))

import data as D                                  # noqa: E402
import depict                                     # noqa: E402
import curate                                     # noqa: E402
import pose3d as p3d                              # noqa: E402
from shared import synthesizability as syn        # noqa: E402

st.set_page_config(page_title="Dance with Inhibition — integration",
                   layout="wide")


# --------------------------------------------------------------------------
# STALE-MODULE GUARD — MUST RUN BEFORE ANY HELPER ATTRIBUTE IS TOUCHED
# --------------------------------------------------------------------------
# Streamlit re-runs THIS file on every interaction but does NOT re-import local
# helper modules: they stay in sys.modules from process start. So after an edit
# to curate.py, the NEW app.py executes against the OLD curate module and dies
# on the first attribute the old one lacks:
#
#     AttributeError: module 'curate' has no attribute 'PANEL_SCOPE'
#
# This check used to sit at the BOTTOM of the file next to the version badge,
# which meant the very error it exists to explain was raised ~30 lines before it
# ever ran. A diagnostic that only fires after the crash it diagnoses is not a
# diagnostic. It runs first now, and it STOPS the page with an instruction
# rather than letting an AttributeError reach the user as a traceback.
#
# Each helper freezes its own mtime at import (LOADED_MTIME); comparing that to
# the file's CURRENT mtime is the only thing that answers "is what I imported
# still what is on disk". Comparing current mtimes to each other -- the first
# version of this check -- flags staleness whenever any file is newer than
# another, which is always, and it never knew when anything was imported.
_HELPERS = (D, depict, curate, p3d)


def stale_modules() -> list[str]:
    """Helper modules whose file changed since this process imported them."""
    out = []
    for m in _HELPERS:
        loaded = getattr(m, "LOADED_MTIME", None)
        if loaded is None:
            continue
        try:
            if Path(m.__file__).stat().st_mtime > loaded + 1:
                out.append(m.__name__)
        except OSError:
            continue
    return out


_STALE = stale_modules()
if _STALE:
    st.error(
        "### Restart required — this process is running stale code\n\n"
        f"`{'`, `'.join(_STALE)}` changed on disk since this process imported "
        "them. Streamlit re-runs `app.py` on every interaction but **does not "
        "re-import helper modules**, so the new `app.py` is calling into old "
        "ones. Rerunning will not fix it, and neither will a browser "
        "refresh.\n\n"
        "Stop the process and start it again:\n\n"
        "```bash\n"
        "/data/lab_vm/envs/dwi_gui/bin/python3 -m streamlit run "
        "integration/app/app.py\n"
        "```\n\n"
        "The page is halted here deliberately. Continuing would raise an "
        "`AttributeError` from whichever old module lacks the newest "
        "attribute, which reads like a code bug rather than a restart.")
    st.stop()


# --------------------------------------------------------------------------
# shared furniture
# --------------------------------------------------------------------------

def gate_badge(verdict: str) -> str:
    return {"STRONG": "🟢", "MODERATE": "🟡", "UNDERPOWERED": "🟠",
            "FAIL": "🔴", "UNGATED": "⚪"}.get(str(verdict).upper(), "⚪")


# --------------------------------------------------------------------------
# the curation filter, applied everywhere (issue #3.2)
#
# IT LIVES IN THE SIDEBAR BECAUSE THE SIDEBAR IS THE ONLY THING ON EVERY PAGE.
# The previous version was a text box inside the Shortlists panel; it wrote the
# spec to session state and nothing else ever read it, so a chemist who excluded
# chlorines met chlorinated molecules again the moment they clicked "Candidate
# dossier". A filter that is invisible from the panel it is not affecting cannot
# be noticed to be off, which is what made this misleading rather than merely
# incomplete.
#
# EVERY PANEL DECLARES ITS SCOPE IN `curate.PANEL_SCOPE`, and an undeclared
# panel raises. Defaulting an unknown panel to "unfiltered" is precisely the
# behaviour being fixed, and defaulting it to "filtered" would silently curate a
# provenance table.
# --------------------------------------------------------------------------

CURATE_KEY = "_curate_spec"


def curation_spec() -> str:
    """The constraint text the chemist has entered, or ""."""
    return str(st.session_state.get(CURATE_KEY, "") or "")


def curation_rules() -> tuple[list, str | None]:
    """(parsed rules, error). An unparseable spec yields ([], the message).

    Parsed here rather than inside each panel so that a typo produces ONE
    message in the sidebar instead of five identical ones down the page — and
    so that every panel agrees about whether a filter is currently active.
    """
    spec = curation_spec()
    if not spec.strip():
        return [], None
    try:
        return curate.parse(spec), None
    except curate.ConstraintError as exc:
        return [], str(exc)


def curation_sidebar() -> None:
    """The global constraint box plus a live statement of what it is doing."""
    st.sidebar.divider()
    st.sidebar.subheader("🧪 Curate")
    st.sidebar.caption(
        "Chemistry you have ruled out. Applies to **every** candidate view.")
    st.sidebar.text_area(
        "constraints", key=CURATE_KEY, height=100,
        placeholder="no chlorine\nmw < 450",
        help="One per line. `no chlorine` · `no [Cl]` (SMARTS works too) · "
             "`require sulfonamide` · `mw < 450` · `sp2_bonds <= 8` for less "
             "conjugation · `rotatable_bonds <= 6`. Known group names: "
             + ", ".join(sorted(curate.NAMED_GROUPS)))

    rules, err = curation_rules()
    if err:
        # REFUSED, NOT IGNORED. A mis-parsed constraint that filters nothing is
        # indistinguishable from a working one that finds nothing, and a
        # chemist would act on it.
        st.sidebar.error(f"**Filter OFF — constraint not understood.** {err}")
        return
    if not rules:
        st.sidebar.caption("No filter active — every panel shows all candidates.")
        return

    # The count is over the pooled shortlists, so the sidebar states the size of
    # the effect before the reader has scrolled to any particular panel.
    pool = D.all_shortlists()
    kept = pool
    if len(pool) and "canonical_smiles" in pool.columns:
        try:
            kept, rules = curate.apply(pool, curation_spec())
        except curate.ConstraintError:
            kept = pool
    st.sidebar.success(
        f"**Filter ON — {len(kept)} of {len(pool)} shortlisted candidates.**\n\n"
        + "\n".join(f"- `{r.text}` −{r.removed}" for r in rules))
    st.sidebar.caption(
        "Not applied to: "
        + ", ".join(s.panel for s in curate.PANEL_SCOPE if not s.filtered)
        + " — see *what the filter never touches* on any curated panel.")


def curation_header(panel: str) -> list:
    """State this panel's relationship to the filter, and return its rules.

    Rendered at the TOP of every panel, filtered or not. An unfiltered panel
    says so explicitly: "no banner" is exactly how the original bug read to the
    user, and silence is not a claim they can check.
    """
    scope = curate.scope_for(panel)
    rules, err = curation_rules()
    if err:
        st.error(f"**Curation filter is OFF — constraint not understood.** {err}"
                 "  \nEverything below is UNFILTERED. Fix the constraint in the "
                 "sidebar.")
        return []
    if not rules:
        return []
    if not scope.filtered:
        st.info(f"**Curation filter is active but not applied here.** {scope.why}")
        return []
    return rules


def curated(df: pd.DataFrame, panel: str, *, note: bool = True,
            label: str = "") -> tuple[pd.DataFrame, list]:
    """Filter one frame for `panel`, rendering the persistent indicator.

    Returns the frame unchanged when the panel's declared scope says the filter
    does not belong there, so a caller cannot accidentally curate a provenance
    view by routing it through this helper.
    """
    scope = curate.scope_for(panel)
    rules, err = curation_rules()
    if err or not rules or not scope.filtered or df.empty:
        return df, []
    if "canonical_smiles" not in df.columns:
        st.warning(f"Curation not applied{' to ' + label if label else ''} — "
                   "this frame carries no `canonical_smiles` column to match on.")
        return df, []
    n_before = len(df)
    try:
        kept, applied = curate.apply(df, curation_spec())
    except curate.ConstraintError as exc:
        st.error(f"**Curation not applied — {exc}**")
        return df, []
    if note:
        st.info(curate.banner(applied, n_before, len(kept), label or None))
    return kept, applied


def unfiltered_facts_note() -> None:
    """What the filter deliberately never touches, and why.

    Sits on the curated panels rather than in documentation because the reader
    who needs it is the one currently looking at a curated table beside an
    uncurated count.
    """
    with st.expander("what the curation filter never touches — and why"):
        st.caption(
            "Curation says which molecules **you** are willing to consider. It "
            "says nothing about what the pipeline generated, docked or "
            "measured, so anything reporting on that population is left alone. "
            "Filtering these would not produce a curated fact; it would "
            "produce a false one.")
        for name, why in curate.UNFILTERED_FACTS:
            st.markdown(f"- **{name}** — {why}")


# THE MURMURENT BRANDED FOOTER, ported from the site's own definition:
# `overrides/partials/copyright.html` + the `.mm-brand-footer` rules in
# `docs/stylesheets/extra.css` in the murmurent repo. Markup, class names,
# colours and logo order all follow that source rather than an approximation of
# the rendered page, so the two stay recognisably the same object.
#
# THE COLOURS ARE THE BRAND, not decoration: #201436 (--purple-deep) with a 3px
# #F0A757 (--tiger) top accent. Dropping them makes this a generic grey footer.
#
# THE LAND ACKNOWLEDGEMENT IS PART OF THE FOOTER, reproduced verbatim including
# the diacritics in "Lūnaapéewak" -- a transliteration stripped of its accents
# is a different word.
#
# LOGOS ARE EMBEDDED AS DATA URIs. Streamlit serves no static directory for
# arbitrary app assets, and an <img src> pointing at a repo path renders as a
# broken icon. They are downscaled copies (3x their CSS height, for retina) of
# the originals in murmurent's docs/assets/logos, kept in this repo so the GUI
# does not depend on the murmurent checkout being present next to it.
ASSETS = APP_DIR / "assets" / "logos"

_LOGO_LINKS = [
    ("lab-logo.png", "https://mikehallett.science", "Hallett Lab",
     "mm-logo-lab"),
    ("schulich.png", "https://www.schulich.uwo.ca/",
     "Schulich School of Medicine &amp; Dentistry", "mm-logo-schulich"),
    ("western.png", "https://www.uwo.ca/", "Western University",
     "mm-logo-western"),
]


def _data_uri(path: Path) -> str | None:
    """base64 data URI for an image, or None if it is missing."""
    try:
        import base64
        return ("data:image/png;base64,"
                + base64.b64encode(path.read_bytes()).decode("ascii"))
    except OSError:
        return None


@st.cache_data(show_spinner=False)
def _footer_html() -> str:
    """Build the footer once; the base64 payload does not change per rerun."""
    imgs = []
    for fname, href, alt, cls in _LOGO_LINKS:
        uri = _data_uri(ASSETS / fname)
        if uri is None:
            # A missing logo must not take the footer -- and must not silently
            # look like a design choice either. The text link stands in.
            imgs.append(f'<a href="{href}" target="_blank" rel="noopener" '
                        f'class="mm-logo-missing">{alt}</a>')
            continue
        imgs.append(
            f'<a href="{href}" target="_blank" rel="noopener">'
            f'<img class="{cls}" src="{uri}" alt="{alt}"></a>')
    return f"""
<style>
.mm-brand-footer {{
  width: 100%;
  background: #201436;                 /* --purple-deep */
  border-top: 3px solid #F0A757;       /* --tiger */
  color: rgba(255, 255, 255, 0.85);
  padding: 16px 24px;
  margin-top: 2.5rem;
  font-size: 0.72rem;
  line-height: 1.5;
  border-radius: 2px;
}}
.mm-brand-bar {{
  max-width: 1220px;
  margin: 0 auto;
  display: flex;
  align-items: center;
  gap: 20px;
  flex-wrap: wrap;
}}
.mm-brand-bar a {{ border-bottom: 0; line-height: 0; }}
.mm-brand-footer img {{ display: block; width: auto !important; max-width: none; }}
.mm-brand-footer img.mm-logo-lab      {{ height: 30px !important; border-radius: 2px; }}
.mm-brand-footer img.mm-logo-schulich {{ height: 22px !important; filter: brightness(0) invert(1); opacity: 0.85; }}
.mm-brand-footer img.mm-logo-western  {{ height: 18px !important; }}
.mm-brand-footer .mm-dept {{ color: rgba(255, 255, 255, 0.72); flex-grow: 1; }}
.mm-brand-footer .mm-logo-missing {{ color: rgba(255,255,255,0.6); line-height: 1.4; }}
.mm-brand-footer .mm-ack {{
  max-width: 1220px;
  margin: 10px auto 0;
  font-style: italic;
  font-size: 0.68rem;
  color: rgba(255, 255, 255, 0.7);
}}
</style>
<div class="mm-brand-footer">
  <div class="mm-brand-bar">
    {"".join(imgs)}
    <span class="mm-dept">Built by the Hallett Lab &middot; Department of
    Biochemistry &middot; London, ON, Canada</span>
  </div>
  <div class="mm-ack">
    We acknowledge that Western University is located on the traditional lands
    of the Anishinaabek, Haudenosaunee, L&#363;naap&#233;ewak and Attawandaron
    peoples.
  </div>
</div>
"""


def site_footer() -> None:
    """Render the murmurent branded footer at the bottom of the page."""
    st.markdown(_footer_html(), unsafe_allow_html=True)


def honest_limits() -> None:
    st.warning(
        "**Honest limits.** No authoritative cross-method ranking exists here — "
        "the four metrics are different quantities.\n\n"
        "**The non-covalent gate has now issued a verdict: WEAK** "
        "(ROC-AUC 0.599, 95% CI [0.311, 0.874], **EF1% 0.0**, over 6 actives / "
        "6 independent chemotypes; D0041). The point estimate is above chance "
        "and the interval contains chance, so docking does not demonstrably "
        "enrich for known Pin1 binders. The covalent gate remains UNDERPOWERED.\n\n"
        "**The rankings are partly size rankings** (D0043). Spearman between "
        "heavy-atom count and each approach's OWN rank metric (all "
        "lower-is-better): −0.617 T₁, −0.479 T₃, −0.230 T₂ — larger "
        "molecules score better in all three. T₄ runs the other way at "
        "+0.181. T₃'s shortlist has a median 39 heavy atoms against 25 "
        "generated. Read every rank with that in mind.\n\n"
        "Inhibition versus activation is unresolved, and there is no wet-lab "
        "ground truth for any candidate.")


def show_gate(stratum: str, metric: str) -> None:
    tok = D.gate_verdicts()
    try:
        g = tok["strata"][stratum]["metrics"][metric]
    except Exception:  # noqa: BLE001
        st.info(f"No gate verdict recorded for {stratum}/{metric} — treat this "
                "ranking as UNGATED.")
        return
    ci = g.get("roc_auc_ci") or [None, None]
    cols = st.columns(4)
    cols[0].metric("gate verdict", f"{gate_badge(g.get('verdict'))} {g.get('verdict')}")
    cols[1].metric("ROC-AUC", f"{g.get('roc_auc', float('nan')):.3f}")
    cols[2].metric("95% CI", f"[{ci[0]:.3f}, {ci[1]:.3f}]" if ci[0] is not None else "—")
    cols[3].metric("EF 1%", f"{g.get('ef_1pct', 0):.1f}")
    for r in g.get("reasons", []):
        st.caption(f"· {r}")


# --------------------------------------------------------------------------
# panels
# --------------------------------------------------------------------------

def row_value(frame: pd.DataFrame, row: pd.Series, name: str) -> str | None:
    """A row's value for `name` as a string, or None when it is absent or NaN.

    Written as a plain function rather than a closure inside `locate_pose`
    because `tests/test_app_names.py` walks each function's own scope and does
    not model enclosing ones — a nested helper reading `frame` from outside
    trips the guard that exists to catch exactly that mistake in panels.
    """
    if name not in frame.columns:
        return None
    val = row.get(name)
    return str(val) if pd.notna(val) else None


def locate_pose(approach: str, frame: pd.DataFrame, row: pd.Series) -> Path | None:
    """The pose file for one row, using the frame's own answer when it has one.

    `pose_path` is recorded by T_3/T_4 and is authoritative; `dock_id` is the
    fallback, because covalent pose files are named by a hash that shares
    nothing with `candidate_id` — deriving the name from the candidate id found
    nothing and rendered as missing data rather than as a lookup bug.
    """
    return p3d.find_pose(approach, str(row["candidate_id"]),
                         dock_id=row_value(frame, row, "dock_id"),
                         pose_path=row_value(frame, row, "pose_path"))


def render_pose_viewer(approach: str, approach_name: str,
                       frame: pd.DataFrame, row: pd.Series, *,
                       key: str, height: int = 620) -> None:
    """The docked-pose viewer: labelled surface, every mode, and an export.

    3Dmol.js has no on-screen controls at all -- no buttons, no hint that
    right-drag zooms or that Ctrl+drag pans. Someone who tries left-drag alone
    finds the view rotates and never centres, and reasonably concludes it is
    broken. The control table is not decoration; it is the only affordance.

    ALL MODES ARE AVAILABLE, AND ONE IS SELECTED BY DEFAULT. Overlaying nine
    sticks by default reads as a single impossible molecule, which is why the
    old viewer showed one — but the fix for that is a control, not a truncation.
    """
    cid = str(row["candidate_id"])
    pose_file = locate_pose(approach, frame, row)
    if pose_file is None:
        st.info(f"no docked pose file found for {cid}")
        return

    poses = p3d.read_poses(pose_file)
    if not poses:
        st.warning(f"`{pose_file.name}` holds no readable binding mode.")
        return

    # WHICH SCORE THE FILE IS SORTED BY, SAID BEFORE THE PICKER. gnina returns
    # modes in CNNscore order and the frame's `affinity_kcal` is read off the
    # FIRST record, so "pose 1" is the CNN's choice and not the affinity's.
    order_score = p3d.FILE_ORDER_SCORE.get(pose_file.suffix.lower().lstrip("."), "")
    aff_score = ("minimizedAffinity" if pose_file.suffix.lower() == ".sdf"
                 else "vina_affinity")
    better = p3d.hidden_better_pose(poses, shown=1, score=aff_score)

    c1, c2, c3 = st.columns([3, 2, 2])
    mode = c1.radio("which modes", ["best (pose 1)", "pick", "overlay all"],
                    horizontal=True, key=f"posemode_{key}",
                    help="'overlay all' answers the question a single pose "
                         "cannot: do the modes agree about where the ligand "
                         "sits, or scatter across the site?")
    if mode == "pick":
        chosen = c2.multiselect("pose(s)", [p.index for p in poses], default=[1],
                                key=f"posepick_{key}")
        show = tuple(chosen) or (1,)
    elif mode == "overlay all":
        show = tuple(p.index for p in poses)
    else:
        show = (1,)
    framing = c3.selectbox(
        "framing", ["ligand", "pocket", "all"], key=f"posezoom_{key}",
        help="Changing this re-renders and RESETS the view — it is also how "
             "you recover from having rotated or panned somewhere unhelpful.")

    if better is not None:
        idx, gap = better
        st.warning(
            f"**Pose 1 is not the best-scoring mode in this file.** Pose {idx} "
            f"has a `{aff_score}` better by **{gap:.2f} kcal/mol**. The file is "
            f"ordered by `{order_score}` and the frame's rank metric is read "
            "off the first record, so the two gnina scores disagree about this "
            "candidate. Switch to *pick* or *overlay all* to see it.")

    components.html(
        p3d.pose_html(pose_file, show=show, zoom_on=framing, height=height),
        height=height + 20)

    st.caption(
        f"**{approach_name}** · `{pose_file.name}` · {len(poses)} binding mode(s), "
        f"showing {len(show)}. The receptor is a **surface** with the three "
        "sub-pockets coloured and labelled; **Cys113 is in yellow sticks**; the "
        "first shown pose is cyan and thicker than the rest. A dashed yellow "
        "bond appears only when a pose is actually within bonding distance of "
        "the Cys113 SG — T₁/T₂ are non-covalent and correctly get none.")

    with st.expander("the three sub-pockets — what each one is", expanded=False):
        st.markdown(p3d.subpocket_legend())
        st.caption(
            "Residue numbers are checked against "
            "`6VAJ_prepared.pdb` at test time (`verify_subpockets`), not "
            "transcribed — a mislabelled pocket is worse than an unlabelled "
            "one. Regions are kept disjoint: the sulfopin paper's wider "
            "proline-pocket description also includes Thr152/His157, which are "
            "coloured with the catalytic tetrad here.")

    with st.expander(f"all {len(poses)} modes and their scores"):
        st.caption(
            "★ marks the best pose on each score. Note the direction differs: "
            "`minimizedAffinity` and `vina_affinity` are kcal/mol and **lower "
            "is better**; `CNNscore`, `CNNaffinity` and `CNN_VS` are "
            "dimensionless and **higher is better**.")
        st.dataframe(pd.DataFrame(p3d.pose_score_table(poses)),
                     use_container_width=True, hide_index=True)

    # THE HAND-OFF (issue #3.1). The embedded viewer cannot measure a distance,
    # mutate a residue or render a figure, and anything a chemist wants to DO
    # with a pose happens in PyMOL or ChimeraX. The alternative to this button
    # is that they hunt for the file by hand — which for T_3/T_4 means knowing
    # that pose files are named by dock_id, the exact lookup that has already
    # caught this project out once.
    is_cov = D.APPROACHES[approach]["stratum"] == "covalent"
    try:
        bundle = p3d.pose_bundle(pose_file, covalent=is_cov)
    except OSError as exc:  # noqa: BLE001 - an unreadable receptor is not fatal
        bundle = None
        st.caption(f"export unavailable: {exc}")
    if bundle is not None:
        st.download_button(
            "⤓ open in PyMOL / ChimeraX — pose file + receptor + session script",
            data=bundle, file_name=f"{cid}_poses.zip", mime="application/zip",
            key=f"posedl_{key}", use_container_width=True,
            help="A zip with all binding modes, the shared prepared receptor, "
                 "and a .pml that reproduces this exact view — sub-pocket "
                 "selections generated from the same definitions the GUI uses, "
                 "so the session and this page cannot disagree. Run "
                 "`pymol view_pin1.pml`.")

    with st.expander("how to move the view", expanded=False):
        st.markdown(p3d.CONTROLS)


@st.dialog("Docked poses", width="large")
def show_pose_dialog(approach: str, approach_name: str, top: pd.DataFrame) -> None:
    """Full-width pose viewer, opened from a shortlist column."""
    ids = list(top["candidate_id"])
    pick = st.selectbox(
        "candidate", ids, key=f"dlg_pose_{approach}",
        format_func=lambda c: (
            f"#{int(top.loc[top.candidate_id == c, 'rank'].iloc[0])} · {c}"))
    prow = top.loc[top.candidate_id == pick].iloc[0]
    render_pose_viewer(approach, approach_name, top, prow,
                       key=f"dlg_{approach}", height=560)


def panel_candidates() -> None:
    st.header("Shortlists — four approaches, side by side")
    st.caption("Each column is ranked by ITS OWN metric. The columns are not "
               "comparable with each other and are deliberately not merged.")
    honest_limits()

    # THE CURATOR FILTER (issue #1 note 1; issue #3.2). Removes rows a chemist
    # has ruled out and preserves the existing order. It deliberately does NOT
    # re-score: the gate has measured that the underlying ranking does not
    # demonstrably enrich (D0041) and is partly a size ranking (D0043), so
    # inventing a second ordering on top of an unvalidated one would compound
    # the problem. Dropping molecules you have ruled out is safe either way.
    #
    # The control itself now lives in the SIDEBAR so it reaches every panel;
    # this panel only consumes it.
    curation_header("Shortlists")

    list_choice = st.radio(
        "which shortlist",
        ["synthesizable (rule failures replaced)", "raw metric top-25"],
        index=0, horizontal=True,
        help="DEFAULT: rule failures are REMOVED from the quota and the "
             "next-best passing candidates take their slots — a different set "
             "of molecules, not a reordering. Switch to the raw list to see "
             "what the metric alone selected.")
    prefer_synth = list_choice.startswith("synthesizable")

    # SAY IT ONCE, AT THE TOP. Which list a reader is looking at changes which
    # molecules exist on the page, so it cannot be something they infer from a
    # per-approach caption further down.
    # NOT CURATED, ON PURPOSE. These counts describe what the synthesizable
    # rebuild did to each approach's whole quota, over the full generated
    # population. Curation happens downstream of the rebuild and cannot change
    # what it did, so filtering them would report a number that was never true.
    _deltas = {k: D.shortlist_delta(k) for k in D.APPROACHES}
    _dropped = sum(d["dropped"] for d in _deltas.values())
    _promoted = sum(d["promoted"] for d in _deltas.values())

    # AN APPROACH WHOSE FRAME HAS NO REBUILD IS NOT A FILTERED APPROACH.
    # `shortlist_delta` reports `available: False` when the frame carries no
    # `shortlist_synth` column, and that flag used to be computed and then
    # ignored — so a T₂ seed that had never been through
    # `reshortlist_synthesizable.py` contributed `dropped: 0` and sat under a
    # green "non-synthesizable compounds have been filtered out" banner while
    # showing the RAW list. Computed-and-unused is the same defect as T_1's
    # `alert_gate_pass`, and it reads as a filter to anyone looking.
    _unfiltered = [k for k, d in _deltas.items() if not d.get("available")]
    if prefer_synth and _unfiltered:
        st.error(
            "### ⚠️ The synthesizability rebuild has NOT been run for "
            + ", ".join(D.display_name(k).split("·")[0].strip()
                        + (f" (seed {D.variant_label(k)})" if D.variant_label(k)
                           else "")
                        for k in _unfiltered)
            + "\nThose columns show the **raw metric top-25 and may contain "
              "compounds that fail a structural synthesizability rule.** Every "
              "other column is filtered, so the lists below are not directly "
              "comparable. Run `scripts/reshortlist_synthesizable.py` to fix.")

    if prefer_synth and _dropped:
        st.success(
            f"### ✅ Non-synthesizable compounds have been filtered out\n"
            f"**{_dropped} candidates that fail a structural synthesizability "
            f"rule were removed** from these shortlists, and **{_promoted} "
            f"next-best passing candidates were promoted** into the freed "
            f"slots. Per approach — "
            + " · ".join(f"{D.APPROACHES[k]['name'].split('·')[0].strip()} "
                         f"−{d['dropped']}"
                         for k, d in _deltas.items() if d["dropped"])
            + ". These are different molecules from the raw metric top-25, not "
              "a re-ordering of it. Switch to *raw metric top-25* above to see "
              "what the score alone selected.")
    elif prefer_synth:
        st.info("Synthesizability filter active — no candidate in any "
                "approach's quota failed a rule, so these lists are identical "
                "to the raw metric top-25.")
    else:
        st.warning(
            f"### ⚠️ Showing the RAW metric list — filter OFF\n"
            f"These lists still contain {_dropped} compound(s) that fail a "
            "structural synthesizability rule (highlighted below).")

    cols = st.columns(len(D.APPROACHES))
    for col, (key, cfg) in zip(cols, D.APPROACHES.items()):
        with col:
            st.subheader(D.display_name(key))
            s = D.shortlist(key, prefer_synth=prefer_synth)
            if s.empty:
                st.info("no shortlist yet")
                continue
            s, _applied = curated(s, "Shortlists",
                                  label=D.display_name(key).split("·")[0].strip())
            if s.empty:
                st.warning("every candidate in this approach was filtered out "
                           "by the curation constraints")
                continue
            verdict = str(s["gate_verdict"].iloc[0]) if "gate_verdict" in s else "UNGATED"
            st.markdown(
                f"{gate_badge(verdict)} **{verdict}** · metric "
                f"`{cfg['metric']}` (lower better) · {cfg['mechanism']} · "
                # The ACTIVE seed, not the config's static one. `cfg["seed"]`
                # is "ATRA" for every T₂ variant, so the caption read "seed:
                # ATRA" under a Guo-Pfizer header.
                f"seed: {D.variant_label(key) or cfg['seed']}")
            # SYNTHESIZABILITY IS SHOWN, NEVER RANKED ON (issue #1). The Lu lab
            # would not make a compound that fails these, so a failure has to be
            # VISIBLE in the list a chemist reads — but the rules are structural
            # red flags, not a score, and folding them into the ordering would
            # add a second unvalidated ranking to one the gate has already
            # measured as not demonstrably enriching (D0041).
            #
            # SAscore is displayed beside them deliberately: the two disagree
            # often, and seeing where is how the rules earn or lose trust.
            s = s.copy()
            s["synth"] = [
                ("⚠ " + ", ".join(r.name for r in syn.violations(str(x)))
                 if syn.violations(str(x)) else "✓")
                for x in s["canonical_smiles"]]
            n_fail = int((s["synth"] != "✓").sum())
            delta = D.shortlist_delta(key)
            using_synth = (s["shortlist_column"].iloc[0] == D.SYNTH_COLUMN
                           if "shortlist_column" in s.columns else False)
            if using_synth:
                if delta["dropped"]:
                    st.success(
                        f"**synthesizable list** — {delta['dropped']} rule "
                        f"failures removed, {delta['promoted']} replacements "
                        "promoted into their slots. These are different "
                        "molecules from the raw top-25, not a reordering.")
                else:
                    st.success(
                        "**synthesizable list** — no candidate in this "
                        "approach's quota failed a rule, so it is identical to "
                        "the raw top-25.")
            elif delta["available"]:
                st.warning(
                    f"**raw metric list** — includes {n_fail} rule "
                    f"failure(s), highlighted. The synthesizable list drops "
                    f"{delta['dropped']} and promotes {delta['promoted']} "
                    "replacements.")
            if n_fail:
                st.error(
                    f"**{n_fail} of {len(s)} fail a structural synthesizability "
                    "rule** and are highlighted below — the `synth` column "
                    "names the rule.")
                st.caption(
                    "`SAscore` sits beside it for comparison (lower = the "
                    "heuristic thinks it is easier to make). The two agree only "
                    "partly — SAscore separates rule-failures from rule-passes "
                    "at AUC 0.74 (T₁), 0.83 (T₂), 0.75 (T₃), where 0.5 would "
                    "be no information. So neither replaces the other: SAscore "
                    "misses specific impossibilities, and the rules say nothing "
                    "about how hard a route is.")

            show = [c for c in ("display_rank", "rank", "candidate_id",
                                "synth", "SAscore",
                                cfg["metric"],
                                "ligand_efficiency", "dG_kcal",
                                "dG_ensemble_interaction_kcal",
                                "dG_ensemble_interaction_sem_kcal",
                                "dG_ensemble_internal_residual_kcal",
                                "dG_ensemble_kcal", "dG_ensemble_sem_kcal",
                                "warhead_class")
                    if c in s.columns]
            # FILTERED, NOT REORDERED (PI decision, issue #1). The synthesizable
            # list is a genuinely different SET: failures are out of the quota
            # and the next-best passers are in. Sorting is therefore just by
            # that list's own rank.
            #
            # Failures are still demoted and highlighted here as a SAFETY NET.
            # On the synthesizable list nothing should trip it -- if a row does
            # light up, the rebuild is stale relative to the rules and that must
            # be visible rather than assumed away.
            s["_synth_fail"] = s["synth"].str.startswith("⚠")
            sort_by = ["_synth_fail"] + (["display_rank"]
                                         if "display_rank" in s.columns
                                         else ["rank"])
            s = s.sort_values(sort_by)
            show_cols = show
            _t = s[show_cols].head(25)
            st.dataframe(
                _t.style.apply(
                    lambda row: ["background-color: rgba(255,75,75,0.16)"
                                 if str(row.get("synth", "✓")).startswith("⚠")
                                 else "" for _ in row], axis=1),
                use_container_width=True, hide_index=True)
            if using_synth and n_fail:
                st.error(
                    f"{n_fail} row(s) on the SYNTHESIZABLE list still fail a "
                    "rule. That should be impossible — rerun "
                    "`scripts/reshortlist_synthesizable.py`; the rules have "
                    "changed since the list was rebuilt.")
            if using_synth and delta["dropped"]:
                st.caption(
                    "`display_rank` is this list's own 1..N order; `rank` is "
                    "the original metric rank, untouched — gaps in it are the "
                    "removed failures.")
            if "dG_kcal" in s.columns and s["dG_kcal"].notna().any():
                st.caption(f"MM-GBSA dG on {int(s['dG_kcal'].notna().sum())} of "
                           f"{len(s)} — an INDEPENDENT estimate, not confirmation "
                           "of the docking rank.")
            if ("dG_ensemble_kcal" in s.columns
                    and s["dG_ensemble_kcal"].notna().any()):
                # The two dG columns are DIFFERENT ESTIMATORS, not two
                # precisions of one (D0036). Said plainly here, because a
                # reader who sees a value beside an uncertainty will otherwise
                # read the pair as one number and its error bar.
                st.caption(
                    f"Ensemble dG on "
                    f"{int(s['dG_ensemble_kcal'].notna().sum())} of {len(s)}, "
                    "from 2 ns of implicit-solvent MD; the uncertainty is the "
                    "SEM widened by the trajectory's statistical inefficiency. "
                    "It is a DIFFERENT estimator from `dG_kcal` (one "
                    "trajectory vs three independent minimisations), so it "
                    "replaces rather than refines it.")

            # WHICH COLUMN IS THE BINDING ENERGY (D0037). Single-trajectory
            # MM-GBSA is defined by its bonded terms cancelling between the
            # legs. For the covalent approaches the decomposition cuts the
            # Cys113 SG-C bond and caps both sides with hydrogen, and those caps
            # do not cancel -- so `dG_ensemble_kcal`, the full potential
            # difference, carries a contamination that has nothing to do with
            # binding. The size of it is measured per approach and shown here
            # rather than described, because "small" and "three times larger
            # than the signal" are both possible and the reader cannot tell
            # which without the number.
            ic, rc = ("dG_ensemble_interaction_kcal",
                      "dG_ensemble_internal_residual_kcal")
            if ic in s.columns and s[ic].notna().any():
                m = s[s[ic].notna()]
                ratio = (m[rc].abs()
                         / m[ic].abs().replace(0, float("nan"))).median()
                if ratio >= 0.10:
                    st.warning(
                        f"**Read `{ic}`, not `dG_ensemble_kcal`.** The "
                        f"link-atom caps leave a residual of "
                        f"{m[rc].median():+.2f} kcal/mol against an interaction "
                        f"energy of {m[ic].median():+.2f} — a median "
                        f"|residual|/|interaction| of **{ratio:.2f}**. The full "
                        "potential column is not an interaction energy here.")
                else:
                    st.caption(
                        f"Bonded terms cancel: median "
                        f"|residual|/|interaction| = {ratio:.2f}, so "
                        f"`{ic}` and `dG_ensemble_kcal` agree. This approach is "
                        "non-covalent, so there is no link atom to contaminate "
                        "the decomposition.")
            with st.expander("structures"):
                top = s.sort_values("rank").head(9)
                smi_col = ("adduct_smiles"
                           if "adduct_smiles" in top.columns
                           and top["adduct_smiles"].notna().any()
                           else "canonical_smiles")
                if smi_col == "adduct_smiles":
                    st.caption("Showing the ADDUCT form — the post-reaction "
                               "species that was docked (D0022/D0030), not the "
                               "molecule as synthesised.")
                st.markdown(
                    depict.grid(list(top[smi_col]),
                                [f"#{int(r)} · {v:.2f}"
                                 for r, v in zip(top["rank"], top[cfg["metric"]])],
                                width=150, height=120),
                    unsafe_allow_html=True)

            # THE POSE OPENS IN A DIALOG, NOT IN THIS COLUMN. A ligand rendered
            # alone is a conformer, not a pose -- but a pose rendered into a
            # quarter-width column is barely better. The 3D canvas was 700 px
            # inside a ~300 px column, so the scene sat off-frame and could not
            # be recovered with the mouse. st.dialog gives it the full width.
            top = s.sort_values("rank").head(9)
            if st.button("view docked poses ⤢", key=f"posebtn_{key}",
                         use_container_width=True):
                show_pose_dialog(key, D.display_name(key), top)

    unfiltered_facts_note()


def panel_dossier() -> None:
    st.header("Per-candidate dossier")
    st.caption("Everything recorded about one candidate, including the SMILES "
               "in a form you can copy into any other tool.")

    # THE DOSSIER IS THE PANEL THE BUG WAS REPORTED AGAINST (issue #3.2). The
    # curation filter was typed into the shortlist panel and stopped there, so
    # a chemist who had excluded a whole chemotype still found it offered here
    # as a candidate in good standing.
    rules = curation_header("Candidate dossier")

    c1, c2 = st.columns([1, 2])
    approach = c1.selectbox("approach", list(D.APPROACHES),
                            format_func=D.display_name)
    s_all = D.shortlist(approach)
    if s_all.empty:
        st.info("no shortlist for this approach yet")
        return
    s_all = s_all.sort_values("rank")

    # AN EXPLICIT ESCAPE HATCH, DEFAULT OFF. Hiding the excluded molecules is
    # the point, but "why exactly did my constraint drop this one" is a fair
    # question and answering it must not require clearing the filter.
    show_excluded = False
    if rules:
        show_excluded = st.checkbox(
            "also list candidates the curation filter excluded", value=False,
            key=f"dossier_show_excluded_{approach}",
            help="Off by default. Excluded candidates are labelled ✗ in the "
                 "picker so an inspected molecule is never mistaken for one "
                 "still in contention.")
    s, _applied = curated(s_all, "Candidate dossier",
                          label=D.display_name(approach).split("·")[0].strip())
    excluded_ids = set(s_all["candidate_id"]) - set(s["candidate_id"])
    if s.empty and not show_excluded:
        st.warning("Every candidate in this approach was excluded by the "
                   "curation constraints. Tick the box above to inspect one "
                   "anyway, or relax the filter in the sidebar.")
        return

    pick_from = s_all if show_excluded else s
    labels = {r["candidate_id"]:
              ("✗ " if r["candidate_id"] in excluded_ids else "")
              + f"#{int(r['rank'])} · {r['candidate_id']}"
              for _, r in pick_from.iterrows()}
    cid = c2.selectbox("candidate", list(labels), format_func=lambda k: labels[k])
    row = pick_from[pick_from["candidate_id"] == cid].iloc[0]
    cfg = D.APPROACHES[approach]
    if cid in excluded_ids:
        st.error(
            "**This candidate is EXCLUDED by your current curation "
            "constraints.** It is shown because you asked to see excluded "
            "candidates; it is not part of the curated shortlist.")

    # Everything below reads columns off the frame the row actually came from.
    # Binding it once here keeps a curated-to-empty frame from silently losing
    # its schema and hiding sections that the candidate does have data for.
    s = pick_from

    # --- structure(s) -----------------------------------------------------
    is_covalent = "adduct_smiles" in s.columns and pd.notna(row.get("adduct_smiles"))
    hl = depict.warhead_smarts(str(row.get("warhead_class", ""))) \
        if "warhead_class" in s.columns else None

    if is_covalent:
        a, b = st.columns(2)
        with a:
            st.markdown("**As synthesised** (pre-reaction)")
            _img = depict.png(row["canonical_smiles"], highlight_smarts=hl,
                              width=420, height=320)
            if _img:
                st.image(_img)
            st.caption("Warhead highlighted. This is the molecule a chemist makes.")
        with b:
            st.markdown("**As docked** (adduct form)")
            _img = depict.png(row["adduct_smiles"], width=420, height=320)
            if _img:
                st.image(_img)
            lg = row.get("leaving_group_smiles")
            st.caption(
                f"Post-reaction species (D0022/D0030). Leaving group `{lg}` is "
                "gone." if pd.notna(lg) and lg else
                "Post-reaction species (D0022/D0030). Nothing leaves in a "
                "Michael addition.")
    else:
        _img = depict.png(row["canonical_smiles"], width=520, height=400)
        if _img:
            st.image(_img)

    # --- SMILES -----------------------------------------------------------
    st.markdown("**SMILES**")
    st.code(row["canonical_smiles"], language="text")
    if is_covalent:
        st.markdown("**Adduct SMILES** (what was actually docked)")
        st.code(row["adduct_smiles"], language="text")
    if "protonated_smiles" in s.columns and pd.notna(row.get("protonated_smiles")):
        if str(row["protonated_smiles"]) != str(row["canonical_smiles"]):
            st.markdown("**At pH 7.4** (what MM-GBSA parameterised)")
            st.code(row["protonated_smiles"], language="text")
            st.caption(f"Formal charge {int(row.get('protonated_charge', 0)):+d} — "
                       "the generator emitted a neutral form.")

    # --- numbers ----------------------------------------------------------
    st.divider()
    verdict = str(row.get("gate_verdict", "UNGATED"))
    cols = st.columns(5)
    cols[0].metric(cfg["metric"], f"{row.get(cfg['metric'], float('nan')):.2f}")
    # NEVER RE-DERIVED UNDER CURATION. "3rd of 1,204 docked" is a fact about
    # the docked population; recomputing the denominator over the curated
    # shortlist would give a number that was never true of anything.
    cols[1].metric("rank", f"{int(row['rank'])} of {int(row.get('group_n_docked', 0))}",
                   help="Out of everything this approach docked — not out of "
                        "the curated shortlist. Curation does not change what "
                        "was ranked.")
    le = row.get("ligand_efficiency")
    cols[2].metric("ligand efficiency", f"{le:.3f}" if pd.notna(le) else "—")
    dg = row.get("dG_kcal")
    cols[3].metric("MM-GBSA dG", f"{dg:.2f}" if pd.notna(dg) else "not scored",
                   help="Single-structure, three independent minimisations.")
    # The ensemble value is shown WITH its uncertainty or not at all: a mean
    # from 100 correlated frames quoted bare invites exactly the false
    # precision the ensemble tier exists to remove.
    dge, sem = row.get("dG_ensemble_kcal"), row.get("dG_ensemble_sem_kcal")
    cols[4].metric(
        "ensemble dG",
        f"{dge:.2f} ± {sem:.2f}" if pd.notna(dge) and pd.notna(sem)
        else "not scored",
        help="2 ns implicit-solvent MD, single-trajectory three-leg rescoring. "
             "± is the SEM widened by the statistical inefficiency. A "
             "DIFFERENT estimator from MM-GBSA dG, not a refinement of it "
             "(D0036).")

    st.markdown(f"{gate_badge(verdict)} **Gate: {verdict}** — this rank is an "
                "ordering the pipeline produced, not evidence of binding (D0031).")

    # Does the docked pose survive the solvent? Every label here names the
    # solvent model, the run length and the tool, because the two rows are
    # produced by different engines with different definitions and a bare
    # side-by-side table invites the reader to subtract them.
    imp_r = row.get("ligand_rmsd_nm_mean")
    exp_r = row.get("explicit_ligand_rmsd_nm_mean")
    if pd.notna(exp_r) or pd.notna(imp_r):
        st.markdown("#### Does the docked pose hold up in solvent?")
        st.caption(
            "Ligand RMSD = how far the ligand moved from its docked pose, in "
            "nanometres, after superposing on the protein backbone so protein "
            "tumbling is removed. Smaller = stayed put. This column IS the same "
            "quantity in both rows and may be compared directly."
        )
        rows_ = []
        if pd.notna(imp_r):
            rows_.append({
                "solvent model": "implicit GB (no water molecules)",
                "engine / length": "OpenMM · 2 ns · 90 frames",
                "ligand RMSD (nm)": round(float(imp_r), 3),
                "frames engaged": (round(float(row["frac_frames_engaged"]), 3)
                                   if pd.notna(row.get("frac_frames_engaged"))
                                   else None),
                "replicates": 1,
            })
        if pd.notna(exp_r):
            rows_.append({
                "solvent model": "explicit TIP3P (~9,200 waters)",
                "engine / length": (
                    f"GROMACS · {row.get('ns_analysed', 10):g} ns · "
                    f"{int(row['n_frames_analysed']) if pd.notna(row.get('n_frames_analysed')) else '?'} frames"),
                "ligand RMSD (nm)": round(float(exp_r), 3),
                "frames engaged": (round(float(row["explicit_frac_frames_engaged"]), 3)
                                   if pd.notna(row.get("explicit_frac_frames_engaged"))
                                   else None),
                "replicates": 1,
            })
        st.dataframe(pd.DataFrame(rows_), use_container_width=True,
                     hide_index=True)

        # PER-REPLICATE RMSD (issue #1: "quickly see MD movies, RMSD plots").
        # Plotted as separate traces, never averaged. D0038 found that two runs
        # of ONE molecule under ONE model diverged 5x and 12x, so the spread
        # between replicates is the result; a mean line would hide precisely
        # what the replicates were run to expose.
        series = p3d.rmsd_series(approach, str(cid))
        if series:
            st.markdown("##### Ligand RMSD, per replicate")
            wide = pd.DataFrame({f"rep{k}": pd.Series([v for _, v in rows],
                                                      index=[t for t, _ in rows])
                                 for k, rows in sorted(series.items())})
            wide.index.name = "time (ns)"
            st.line_chart(wide, height=260)
            finals = {k: rows[-1][1] for k, rows in series.items()}
            means = {k: sum(v for _, v in rows) / len(rows)
                     for k, rows in series.items()}
            spread = (max(means.values()) / min(means.values())
                      if min(means.values()) > 0 else float("inf"))
            st.caption(
                f"{len(series)} replicates · mean RMSD "
                f"{min(means.values()):.2f}–{max(means.values()):.2f} nm "
                f"(**{spread:.1f}× spread between replicates**) · final "
                f"{min(finals.values()):.2f}–{max(finals.values()):.2f} nm. "
                "A large spread means this candidate's residence is not "
                "reproducible and no per-molecule verdict should be read from "
                "any single run.")
        else:
            st.warning(
                "**Single run each — do not read a per-molecule verdict from "
                "this.** Velocities are drawn afresh each run, and re-running "
                "the implicit simulation changed one molecule's RMSD by 5x and "
                "another's by 12x (D0038). Five replicates per candidate are "
                "needed before \"this one is unstable\" is a claim about the "
                "molecule rather than about one trajectory.")

        # THE MOVIE. Built on demand: converting a trajectory takes seconds and
        # writes a file, so it is not done for every candidate a reader clicks.
        traj = p3d.find_trajectory(approach, str(cid))
        if traj is not None:
            with st.expander("MD movie (explicit solvent, PBC-corrected)"):
                st.caption(
                    "Protein and ligand only — water and ions are ~90% of the "
                    "30,000-atom system and would show nothing you need. This "
                    "is the **PBC-corrected** trajectory: the raw file shows "
                    "the ligand jumping across the box whenever it crosses a "
                    "periodic boundary, which reads as dissociation and is not "
                    "(D0038).")
                if st.button("build / show movie", key=f"movie_{cid}"):
                    try:
                        html = p3d.trajectory_html(*traj)
                    except p3d.MovieTooLarge as exc:
                        html = None
                        st.error(str(exc))
                    if html:
                        components.html(html, height=640)
                        st.markdown(p3d.CONTROLS)
                    else:
                        st.info("could not build the movie for this trajectory")

        if pd.notna(exp_r) and pd.notna(imp_r):
            st.caption(
                "**The two models disagree.** Across the 47 candidates run "
                "under both, ligand RMSD correlates at Spearman **−0.120** — "
                "effectively not at all. Implicit solvent also reports "
                "systematically more movement (mean 0.78 nm against 0.35 nm), "
                "which is expected: it has no water to hold a ligand in place. "
                "In explicit water nothing dissociates — the largest mean "
                "displacement across all 47 is 1.16 nm.")

        if bool(row.get("explicit_rmsd_suspect", False)):
            st.error(
                "**Suspect measurement.** Displacement approaches half the "
                "periodic box, the signature of a ligand crossing a boundary "
                "and being measured against its own image. Inspect the "
                "trajectory before reporting this value.")

        cm = row.get("gmx_contacts_mean")
        if pd.notna(cm):
            st.caption(
                f"Contacts (GROMACS `mindist`, 0.45 nm): **{cm:g}**. NOT "
                "comparable with the implicit tier's contact count, which "
                "counts heavy-atom pairs and returns several-fold larger "
                "numbers on the same complex. Different definition, not a "
                "solvent effect. Explicit-solvent trajectories are "
                "periodic-boundary corrected before measurement.")

    axes = [a for a in D.SHARED_AXES if a in s.columns]
    if axes:
        st.markdown("**Shared physicochemical axes** (identical RDKit call "
                    "across all four approaches)")
        st.dataframe(pd.DataFrame([{a: row[a] for a in axes}]),
                     use_container_width=True, hide_index=True)

    # --- the docked pose, in the pocket ------------------------------------
    # THE DOSSIER IS WHERE "EVERYTHING RECORDED ABOUT ONE CANDIDATE" LIVES, and
    # until now the pose was reachable only from a dialog on another panel. It
    # is full width here, which is the width the viewer actually needs.
    st.divider()
    st.markdown("#### Docked pose in the Pin1 pocket")
    render_pose_viewer(approach, D.display_name(approach), s, row,
                       key=f"dossier_{approach}_{cid}", height=600)

    # Structural convergence: did any OTHER approach reach this molecule?
    # Reported as a rank when the identical molecule was ranked elsewhere, and
    # as the nearest shortlisted analogue otherwise -- an empty panel would let
    # "nobody else found it" read as "the lookup is broken".
    #
    # DELIBERATELY NOT CURATED. "Was this molecule also ranked by another
    # approach" is a fact about the pipeline's output; a chemist's constraints
    # cannot un-rank it, and filtering the lookup would turn a real convergence
    # into an apparent absence of one.
    smi = row.get("canonical_smiles")
    if pd.notna(smi):
        conv = D.cross_approach_ranks(str(smi), approach)
        if conv:
            st.markdown("**Also found by other approaches?**")
            exact = [c for c in conv if c["exact"]]
            st.dataframe(pd.DataFrame([{
                "approach": c["name"],
                "same molecule?": "YES" if c["exact"] else "no",
                "rank there": (f"{int(c['rank'])} of {c['n_ranked']}"
                               if c["exact"] and pd.notna(c["rank"]) else "—"),
                "nearest in their shortlist": c["candidate_id"] or "—",
                "Tanimoto": c["similarity"],
            } for c in conv]), use_container_width=True, hide_index=True)
            if exact:
                st.success(
                    f"Surfaced independently by {len(exact)} other approach(es)."
                    " Convergence is a soft cross-validation that needs no "
                    "shared metric — the approaches report incommensurable "
                    "numbers, but agreeing on a molecule requires no units.")
            else:
                st.caption(
                    "No other approach reached this molecule. That is the norm "
                    "here, not a gap: exact overlap between every pair of "
                    "approaches is zero, and the closest cross-approach "
                    "shortlist pair in the whole build is T₃~T₄ at Tanimoto "
                    "0.455. The four searches are effectively disjoint, so "
                    "convergence currently provides no cross-validation.")

    flags = {k: row[k] for k in ("shortlist_reason", "reactivity_flag",
                                 "adduct_approximation", "excused_alert_names",
                                 "rgroup_alert_names", "size_class")
             if k in s.columns and pd.notna(row.get(k)) and str(row.get(k)).strip()}
    if flags:
        st.markdown("**Flags carried with this candidate**")
        for k, v in flags.items():
            st.markdown(f"- `{k}`: {v}")

    unfiltered_facts_note()


def panel_convergence() -> None:
    st.header("Structural convergence")
    st.markdown(
        "A molecule surfaced independently by more than one approach is a soft "
        "cross-validation that relies on **no** score commensurability — the "
        "most defensible cross-approach signal available here.")
    st.caption(
        "Read it with care: all four approaches dock into the same receptor "
        "under related protocols, so their errors are correlated and agreement "
        "may report a shared bias rather than a real signal. T₁ is seed-free, "
        "so its agreement with a seeded approach at least is not ancestral.")
    curation_header("Convergence")

    pool = D.all_shortlists()
    if pool.empty or "canonical_smiles" not in pool.columns:
        st.info("no shortlists to compare yet")
        return
    # Pairs are drawn FROM the shortlists, so a pair whose members the chemist
    # has ruled out is not evidence about anything still in contention.
    pool, _cur = curated(pool, "Convergence", label="pooled shortlists")
    if pool.empty:
        st.warning("every shortlisted candidate was excluded by the curation "
                   "constraints — no pairs left to compare")
        return

    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import rdFingerprintGenerator as fpg
    RDLogger.DisableLog("rdApp.*")
    gen = fpg.GetMorganGenerator(radius=2, fpSize=2048)

    pool = pool.reset_index(drop=True)
    fps, keep = [], []
    for i, smi in enumerate(pool["canonical_smiles"]):
        m = Chem.MolFromSmiles(str(smi))
        if m is not None:
            fps.append(gen.GetFingerprint(m))
            keep.append(i)
    sub = pool.loc[keep].reset_index(drop=True)

    thresh = st.slider("ECFP4 Tanimoto threshold", 0.4, 1.0, 0.7, 0.05)
    pairs = []
    for i in range(len(fps)):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[i + 1:])
        for j, sim in enumerate(sims, start=i + 1):
            if sim >= thresh and sub.at[i, "approach"] != sub.at[j, "approach"]:
                pairs.append({
                    "approach A": sub.at[i, "approach"],
                    "candidate A": sub.at[i, "candidate_id"],
                    "approach B": sub.at[j, "approach"],
                    "candidate B": sub.at[j, "candidate_id"],
                    "Tanimoto": round(sim, 3)})
    if pairs:
        st.dataframe(pd.DataFrame(pairs).sort_values("Tanimoto", ascending=False),
                     use_container_width=True, hide_index=True)
    else:
        st.info(f"No cross-approach pair reaches Tanimoto {thresh:.2f}. "
                "That is itself informative — the approaches are exploring "
                "genuinely different regions.")

    unfiltered_facts_note()


def panel_axes() -> None:
    st.header("Shared physicochemical axes")
    st.caption("Computed by the identical RDKit call for all four approaches "
               "(shared/descriptors.py), which is what makes pooling these — "
               "and only these — legitimate.")
    curation_header("Shared axes")

    pool = D.all_shortlists()
    if pool.empty:
        st.info("no shortlists yet")
        return
    pool, applied = curated(pool, "Shared axes", label="pooled shortlists")
    if pool.empty:
        st.warning("every shortlisted candidate was excluded by the curation "
                   "constraints — no distribution left to summarise")
        return

    axes = [a for a in D.SHARED_AXES if a in pool.columns]
    if not axes:
        st.info("no descriptor columns on the frames")
        return
    axis = st.selectbox("axis", axes)

    # A NUMERIC CONSTRAINT TRUNCATES THE AXIS IT NAMES, so the median below is
    # partly the filter's doing rather than the approach's. `curate.PROPERTIES`
    # and `shared/descriptors.py` make the identical RDKit call for each of
    # these pairs, so this is an exact truncation, not an approximate one — and
    # a reader comparing "T_1 median MW" across two filter settings would
    # otherwise read their own constraint as a property of T_1.
    bounded = curate.bounded_axes(applied)
    if axis in bounded:
        st.warning(
            f"**`{axis}` is bounded by your own constraint `{bounded[axis]}`.** "
            "The median below describes the molecules the filter left, not what "
            "the approach produced. Clear that constraint before reading this "
            "axis as a property of the approach.")
    other = {k: v for k, v in bounded.items() if k != axis}
    if other:
        st.caption("Also bounded by your constraints, in the table below: "
                   + ", ".join(f"`{k}` ({v})" for k, v in other.items()))

    st.bar_chart(pool.pivot_table(index="approach", values=axis, aggfunc="median"))
    st.dataframe(
        pool.groupby("approach")[axes].median().round(2),
        use_container_width=True)

    unfiltered_facts_note()


def panel_within_stratum() -> None:
    st.header("Within-stratum re-score — two leaderboards, never one")
    st.caption("Cross-stratum ordering is not implied and is not offered.")
    curation_header("Within-stratum")

    fps = D.protocol_fingerprints()
    t3, t4 = fps.get("t3", set()), fps.get("t4", set())
    if t3 and t4 and t3 != t4:
        st.error(
            f"**Within-covalent comparison DISABLED.** T₃ and T₄ recorded "
            f"different protocol fingerprints ({sorted(t3)} vs {sorted(t4)}). "
            "They did not dock under identical rules, so their affinities are "
            "not comparable. Re-dock the lagging approach before comparing — "
            "showing the numbers anyway is exactly what the fingerprint exists "
            "to prevent.")
    else:
        st.success(f"Covalent protocol fingerprints agree: {sorted(t3 or t4)}")

    for label, keys, stratum, metric in (
            ("Non-covalent (T₁ + T₂)", ("t1", "t2"), "non_covalent", "vina_affinity"),
            ("Covalent (T₃ + T₄)", ("t3", "t4"), "covalent", "affinity_kcal")):
        st.subheader(label)
        # The gate is measured on known actives against property-matched
        # decoys. No candidate of ours is in that calculation, so curation has
        # nothing to say about it and it is shown unfiltered.
        show_gate(stratum, metric)
        if stratum == "covalent" and t3 and t4 and t3 != t4:
            st.info("leaderboard withheld — see the fingerprint mismatch above")
            continue
        frames = [D.shortlist(k) for k in keys]
        frames = [f for f in frames if len(f) and metric in f.columns]
        if not frames:
            st.info("no candidates yet")
            continue
        pooled = pd.concat(frames, ignore_index=True)
        pooled, _cur = curated(pooled, "Within-stratum", label=label)
        if pooled.empty:
            st.warning("every candidate in this stratum was excluded by the "
                       "curation constraints")
            continue
        cols = [c for c in ("approach", "candidate_id", metric,
                            "ligand_efficiency", "dG_kcal") if c in pooled.columns]
        st.dataframe(pooled.sort_values(metric)[cols].head(20),
                     use_container_width=True, hide_index=True)

    unfiltered_facts_note()


def panel_decisions() -> None:
    st.header("Choreography decision log")
    curation_header("Decisions")
    ds = D.decisions_all()
    if not ds:
        st.info("no decision records found")
        return
    c1, c2 = st.columns(2)
    appr = c1.multiselect("approach", sorted({d.get("approach", "?") for d in ds}))
    orig = c2.multiselect("origin", sorted({d.get("origin", "?") for d in ds}))
    st.caption("`origin: adversary` records are the audit trail that adversarial "
               "review actually changed the design.")

    rows = [d for d in ds
            if (not appr or d.get("approach") in appr)
            and (not orig or d.get("origin") in orig)]
    rows.sort(key=lambda d: str(d.get("id")), reverse=True)

    for d in rows:
        superseded = bool(d.get("superseded_by"))
        title = f"~~{d['id']} · {d['title']}~~" if superseded else f"{d['id']} · {d['title']}"
        badge = "⛔" if superseded else {"accepted": "✅", "proposed": "❓"}.get(
            str(d.get("status")), "•")
        with st.expander(f"{badge} {title}  ·  {d.get('approach')} / {d.get('origin')}"):
            if superseded:
                st.warning(f"Superseded by {d['superseded_by']}. Shown, not "
                           "hidden — why the answer changed is usually more "
                           "informative than the current answer.")
            for sec in ("context", "decision", "consequences"):
                if d.get(sec):
                    st.markdown(f"**{sec.title()}**")
                    st.markdown(d[sec])
            if d.get("evidence"):
                st.markdown("**Evidence** (numeric by design — the numbers are "
                            "the argument)")
                for e in d["evidence"]:
                    st.markdown(f"- `{e}`")
            if d.get("affects"):
                st.caption("affects: " + ", ".join(d["affects"]))


def panel_why_this_file() -> None:
    st.header("Why is this file like this?")
    st.caption("The panel that replaces grepping four formats.")
    curation_header("Why this file?")
    frag = st.text_input("path fragment", placeholder="warhead_classes / receptor.yaml")
    if not frag:
        return
    hits = D.decisions_affecting(frag)
    if not hits:
        st.info(f"No decision record names anything matching {frag!r}.")
        return
    for d in hits:
        st.markdown(f"**{d['id']} · {d['title']}**  ({d.get('status')})")
        if d.get("decision"):
            st.markdown(d["decision"])
        st.divider()


def panel_provenance() -> None:
    st.header("Run provenance")
    curation_header("Provenance")
    approach = st.selectbox("approach", list(D.APPROACHES),
                            format_func=D.display_name)
    ms = D.manifests(approach)
    if not ms:
        st.info("no manifests recorded for this approach")
        return
    dirty = [m for m in ms if (m.get("git") or {}).get("dirty")]
    if dirty:
        st.error(
            f"**{len(dirty)} of {len(ms)} runs were made from a DIRTY working "
            "tree.** The recorded commit does not fully describe the code that "
            "ran, so those outputs are provisional.")
    for m in ms:
        g = m.get("git") or {}
        flag = "⚠️ DIRTY" if g.get("dirty") else "clean"
        with st.expander(f"{m.get('stage', '?')} · {m.get('run_id', m['_file'])} · {flag}"):
            st.json(m, expanded=False)


def panel_open_questions() -> None:
    st.header("Open questions")
    curation_header("Open questions")
    q = D.open_questions()
    if not q:
        st.success("Nothing is marked proposed, pending or unverified.")
        return
    st.dataframe(pd.DataFrame(q)[["id", "title", "status", "approach"]],
                 use_container_width=True, hide_index=True)


# --------------------------------------------------------------------------

def panel_seed_comparison() -> None:
    """Each T₂ seed's own top-25, side by side — presented, never pooled.

    THIS IS THE QUESTION THE RESEEDING WAS RUN TO ANSWER: does the starting
    point determine what comes out? Five CReM neighbourhoods, everything else
    held fixed.

    WHY THE SHORTLISTS ARE SHOWN SIDE BY SIDE RATHER THAN MERGED. Each seed's
    top-25 is ranked WITHIN its own pool, and the pools run 1,882-16,806
    molecules. A merged top-25 over unmatched N would mostly report which seed
    generated most, and a merged sort on `vina_affinity` would compare four
    charge states and four chemotypes on a scoring function measured at chance
    for this pocket (D0041) and at 5% pose recovery (D0046). That is exactly
    the cross-approach merge the choreography refuses at the top level, one
    scale down — so the same answer applies: present, do not merge.

    WHAT IS COMPARABLE, and is shown below the shortlists: the score-free axes
    (size, lipophilicity, novelty against the external reference set) and
    whether different seeds ever reach the SAME molecule. Those need no shared
    metric and no matched N.
    """
    st.header("T₂ seed comparison")
    curation_header("T₂ seed comparison")
    status = D.variant_status("t2")
    ready = [v for v in status if v["ready"]]
    pending = [v for v in status if not v["ready"]]

    if not ready:
        st.info("No T₂ seed has a ranked frame yet.")
        return
    st.caption(
        "Each column is that seed's own top-25, ranked within its own pool. "
        "**The scores are not comparable across columns** — different pool "
        "sizes, different chemotypes, and a rank metric the enrichment gate "
        "has not validated (D0041). Read the columns, not a winner.")
    if pending:
        st.warning("Still docking, so absent below: "
                   + ", ".join(f"**{v['label']}**" for v in pending))

    lists = {v["key"]: D.shortlist_variant("t2", v["key"]) for v in ready}

    # WRAPPED AT THREE. With all five seeds plus the degree-2 sample there are
    # six columns, and st.columns(6) makes each too narrow to read a SMILES or
    # a candidate id in. Three per row keeps them legible.
    PER_ROW = 3
    for start in range(0, len(ready), PER_ROW):
        row = ready[start:start + PER_ROW]
        cols = st.columns(PER_ROW)          # fixed width, so a short last row
        for col, v in zip(cols, row):       # does not stretch its columns
            s = lists[v["key"]]
            with col:
                st.subheader(v["label"])
                # THE DEGREE-2 SAMPLE IS NOT A FIFTH SEED. It is a second CReM
                # edit from ATRA, so it answers "does another edit help?" and
                # not "does the starting point matter?". Read next to ATRA it is
                # a comparison; read next to Du-Xu it is a category error.
                if v["key"] == "atra_degree2":
                    st.caption("⚠️ derived from ATRA — compare with the ATRA "
                               "column, not with the other seeds")
                st.caption(f"{v['n']:,} docked · top {len(s)}")
                show = [c for c in ("display_rank", "candidate_id",
                                    "metric_value", "HAC", "SAscore")
                        if c in s.columns]
                st.dataframe(s.sort_values("display_rank")[show],
                             hide_index=True, use_container_width=True)

    st.divider()
    st.subheader("Score-free axes — these ARE comparable")
    st.caption(
        "Medians over each seed's top-25. No shared metric is involved, so "
        "unmatched pool sizes do not distort them. This is where a real "
        "seed effect would show: a seed whose neighbourhood is systematically "
        "larger, greasier or less novel than another's.")
    rows = []
    for v in ready:
        s = lists[v["key"]]
        row = {"seed": v["label"], "n docked": v["n"], "shortlist": len(s)}
        for ax in D.SHARED_AXES:
            if ax in s.columns:
                row[ax] = round(float(pd.to_numeric(s[ax], errors="coerce")
                                      .median()), 3)
        rows.append(row)
    st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)

    st.divider()
    st.subheader("Do two seeds ever reach the same molecule?")
    st.caption(
        "Exact identity on InChIKey across the shortlists. CReM edits a seed, "
        "so overlap would mean two different starting points converged on one "
        "structure — a genuine cross-validation that needs no scoring "
        "function. Absence is a result too, and is reported as one.")
    keys = {}
    for v in ready:
        s = lists[v["key"]]
        if "canonical_smiles" not in s.columns:
            continue
        for _, r in s.iterrows():
            k = _inchikey(str(r["canonical_smiles"]))
            if k:
                keys.setdefault(k, []).append((v["label"], r.get("candidate_id")))
    shared = {k: hits for k, hits in keys.items()
              if len({h[0] for h in hits}) > 1}
    if shared:
        st.dataframe(pd.DataFrame(
            [{"InChIKey": k, "seeds": ", ".join(sorted({h[0] for h in hits})),
              "ids": ", ".join(str(h[1]) for h in hits)}
             for k, hits in shared.items()]),
            hide_index=True, use_container_width=True)
    else:
        st.info(
            f"**No molecule appears in more than one seed's top-25** "
            f"({len(keys)} distinct structures across {len(ready)} seeds). "
            "The neighbourhoods are disjoint at the shortlist level, which is "
            "what degree-bounded editing of chemically distinct seeds should "
            "produce — it means the seed choice, not the search, determines "
            "what you get.")


@st.cache_data(show_spinner=False)
def _inchikey(smiles: str) -> str | None:
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    m = Chem.MolFromSmiles(smiles)
    return Chem.MolToInchiKey(m) if m is not None else None


PANELS = {
    "Shortlists": panel_candidates,
    "T₂ seed comparison": panel_seed_comparison,
    "Candidate dossier": panel_dossier,
    "Convergence": panel_convergence,
    "Shared axes": panel_axes,
    "Within-stratum": panel_within_stratum,
    "Decisions": panel_decisions,
    "Why this file?": panel_why_this_file,
    "Provenance": panel_provenance,
    "Open questions": panel_open_questions,
}

st.sidebar.title("Dance with Inhibition")
st.sidebar.caption("Integration is a presentation and human-decision layer. "
                   "It surfaces and organises evidence; it does not output a "
                   "winner.")
choice = st.sidebar.radio("panel", list(PANELS))


def seed_sidebar() -> None:
    """Choose which seed neighbourhood a multi-variant approach shows.

    T₂ is five CReM neighbourhoods plus a degree-2 sample, each in its own
    experiment directory. Until this existed the app read a hardcoded ATRA
    path, so seeds that had finished and been ranked were simply unreachable.

    SEEDS THAT ARE NOT READY ARE NAMED, NOT OMITTED. A selector that silently
    lists three of six reads as "T₂ has three seeds", which is a different and
    false statement. They are listed below the control with the reason.

    The control sets ONE seed for the whole app rather than pooling them.
    Pooling is not a display choice here: the pools run 1,882-16,806 molecules,
    and a top-N over unmatched N would be dominated by whichever seed generated
    most. That comparison needs a decision (#6), not a widget.
    """
    for approach in D.APPROACHES:
        status = D.variant_status(approach)
        if not status:
            continue
        ready = [v for v in status if v["ready"]]
        pending = [v for v in status if not v["ready"]]
        name = D.APPROACHES[approach]["name"].split("·")[0].strip()
        st.sidebar.divider()
        st.sidebar.subheader(f"🌱 {name} seed")
        if not ready:
            st.sidebar.warning(f"No {name} seed has a ranked frame yet.")
            continue
        keys = [v["key"] for v in ready]
        labels = {v["key"]: f"{v['label']} — {v['n']:,} docked" for v in ready}
        current = D.active_variant(approach)
        chosen = st.sidebar.selectbox(
            "seed neighbourhood", keys,
            index=keys.index(current) if current in keys else 0,
            format_func=lambda k: labels[k], key=f"seed_{approach}")
        D.set_variant(approach, chosen)
        if pending:
            st.sidebar.caption(
                "Not yet available — " + "; ".join(
                    f"**{v['label']}** ({v['why']})" for v in pending))
        st.sidebar.caption(
            "Scores are NOT comparable across seeds: the pools differ ~9× in "
            "size and each inherits its seed's chemotype. Compare property, "
            "novelty and chemotype distributions instead.")


seed_sidebar()

# EVERY PANEL NAME MUST HAVE A DECLARED CURATION SCOPE. Checked here at start-up
# rather than when someone clicks the tab: a panel added without a scope should
# fail immediately and visibly, not filter (or fail to filter) by accident three
# clicks into a session. Same class of defect as the shortlist column that was
# selected by the name already there rather than by the one that answered the
# question.
_declared = {s.panel for s in getattr(curate, "PANEL_SCOPE", ())}
_undeclared = [p for p in PANELS if p not in _declared]
if _undeclared:
    st.sidebar.error(
        f"**Panels with no curation scope declared:** {', '.join(_undeclared)}. "
        "Add them to `curate.PANEL_SCOPE`.")

# The curation control lives HERE, not inside a panel, because the sidebar is
# the only thing present on every page (issue #3.2).
curation_sidebar()
st.sidebar.divider()

# WHICH CODE IS ACTUALLY LOADED. Streamlit re-runs this script on every
# interaction but does NOT re-import local helper modules -- they stay in
# sys.modules from process start. So editing pose3d.py and clicking around
# gives a TypeError about an argument the file on disk plainly has, and the
# obvious conclusion ("the fix didn't work") is wrong. Showing the loaded
# commit next to the working-tree commit makes a stale process visible.
try:
    import subprocess as _sp

    _repo = str(APP_DIR.parent.parent)
    _head = _sp.run(["git", "-C", _repo, "rev-parse", "--short", "HEAD"],
                    capture_output=True, text=True).stdout.strip()
    st.sidebar.caption(f"code: `{_head}`")

    # The stale-module check itself now runs at the TOP of this file, before any
    # helper attribute is dereferenced -- see STALE-MODULE GUARD. Reaching this
    # line means it already passed, so there is nothing to repeat here.
except Exception:  # noqa: BLE001 - a version badge must never break the page
    pass

st.sidebar.caption("The GUI reads; it does not own (D0008). Everything shown is "
                   "a rendering of files in the repo or the append-only tree.")
PANELS[choice]()

# Last, so it sits below whichever panel rendered.
site_footer()
