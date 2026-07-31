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
# shared furniture
# --------------------------------------------------------------------------

def gate_badge(verdict: str) -> str:
    return {"STRONG": "🟢", "MODERATE": "🟡", "UNDERPOWERED": "🟠",
            "FAIL": "🔴", "UNGATED": "⚪"}.get(str(verdict).upper(), "⚪")


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
        "lower-is-better): −0.617 T_1, −0.479 T_3, −0.230 T_2 — larger "
        "molecules score better in all three. T_4 runs the other way at "
        "+0.181. T_3's shortlist has a median 39 heavy atoms against 25 "
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

@st.dialog("Docked pose", width="large")
def show_pose_dialog(approach: str, approach_name: str, top: pd.DataFrame) -> None:
    """Full-width pose viewer with its controls spelled out.

    3Dmol.js has no on-screen controls at all -- no buttons, no hint that
    right-drag zooms or that Ctrl+drag pans. Someone who tries left-drag alone
    finds the view rotates and never centres, and reasonably concludes it is
    broken. The control table is not decoration; it is the only affordance.
    """
    ids = list(top["candidate_id"])
    c1, c2 = st.columns([3, 2])
    pick = c1.selectbox(
        "candidate", ids, key=f"dlg_pose_{approach}",
        format_func=lambda c: (
            f"#{int(top.loc[top.candidate_id == c, 'rank'].iloc[0])} · {c}"))
    framing = c2.selectbox(
        "framing", ["ligand", "pocket", "all"], key=f"dlg_zoom_{approach}",
        help="Changing this re-renders and RESETS the view — it is also how "
             "you recover from having rotated or panned somewhere unhelpful.")

    prow = top.loc[top.candidate_id == pick].iloc[0]
    # dock_id, not candidate_id: covalent pose files are named by a separate
    # hash that shares nothing with the candidate id.
    pose = p3d.find_pose(
        approach, str(pick),
        dock_id=(str(prow["dock_id"]) if "dock_id" in top.columns
                 and pd.notna(prow.get("dock_id")) else None))
    if pose is None:
        st.info(f"no docked pose file found for {pick}")
        return

    components.html(p3d.pose_html(pose, zoom_on=framing, height=620), height=640)
    st.caption(
        f"**{approach_name}** · `{pose.name}` · receptor 6VAJ as a cartoon, "
        "**Cys113 in yellow sticks**, ligand in cyan. The grey surface is the "
        "pocket wall within 12 residues of Cys113 — a whole-protein surface "
        "would hide the ligand."
        + (" Vina writes all nine binding modes into one file; only the "
           "top-scoring one is drawn." if pose.suffix == ".pdbqt" else ""))
    with st.expander("how to move the view", expanded=True):
        st.markdown(p3d.CONTROLS)


def panel_candidates() -> None:
    st.header("Shortlists — four approaches, side by side")
    st.caption("Each column is ranked by ITS OWN metric. The columns are not "
               "comparable with each other and are deliberately not merged.")
    honest_limits()

    # THE CURATOR FILTER (issue #1, general note 1). Removes rows a chemist has
    # ruled out and preserves the existing order. It deliberately does NOT
    # re-score: the gate has measured that the underlying ranking does not
    # demonstrably enrich (D0041) and is partly a size ranking (D0043), so
    # inventing a second ordering on top of an unvalidated one would compound
    # the problem. Dropping molecules you have ruled out is safe either way.
    with st.expander("curate — filter these lists by chemistry you have ruled out"):
        st.caption(
            "One constraint per line. `no chlorine` · `no [Cl]` (SMARTS works "
            "too) · `require sulfonamide` · `mw < 450` · `sp2_bonds <= 8` for "
            "less conjugation · `rotatable_bonds <= 6`. "
            f"Known group names: {', '.join(sorted(curate.NAMED_GROUPS))}.")
        spec = st.text_area("constraints", value="", height=90,
                            placeholder="no chlorine\nmw < 450")
    st.session_state["_curate_spec"] = spec

    demote = st.checkbox(
        "sort candidates failing a synthesizability rule to the bottom",
        value=True,
        help="On by default: a compound the lab would not make should not hold "
             "a top-25 slot. Nothing is deleted and the original metric rank is "
             "kept in the `rank` column. Uncheck to see the raw metric order.")

    cols = st.columns(len(D.APPROACHES))
    for col, (key, cfg) in zip(cols, D.APPROACHES.items()):
        with col:
            st.subheader(cfg["name"])
            s = D.shortlist(key)
            if s.empty:
                st.info("no shortlist yet")
                continue
            if spec.strip():
                try:
                    n_before = len(s)
                    s, applied = curate.apply(s, spec)
                    st.info(
                        f"**curated: {len(s)} of {n_before}** · "
                        + " · ".join(f"`{r.text}` −{r.removed}" for r in applied))
                except curate.ConstraintError as exc:
                    # Refused, not silently ignored: a mis-parsed constraint
                    # that filters nothing is indistinguishable from a working
                    # one that finds nothing, and a chemist would act on it.
                    st.error(f"constraint not applied — {exc}")
                if s.empty:
                    st.warning("every candidate was filtered out")
                    continue
            verdict = str(s["gate_verdict"].iloc[0]) if "gate_verdict" in s else "UNGATED"
            st.markdown(
                f"{gate_badge(verdict)} **{verdict}** · metric "
                f"`{cfg['metric']}` (lower better) · {cfg['mechanism']} · "
                f"seed: {cfg['seed']}")
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
            if n_fail:
                st.error(
                    f"**{n_fail} of {len(s)} fail a structural synthesizability "
                    "rule** and are highlighted below. They are NOT removed and "
                    "NOT reranked — the `synth` column names the rule.")
                st.caption(
                    "`SAscore` sits beside it for comparison (lower = the "
                    "heuristic thinks it is easier to make). The two agree only "
                    "partly — SAscore separates rule-failures from rule-passes "
                    "at AUC 0.74 (T_1), 0.83 (T_2), 0.75 (T_3), where 0.5 would "
                    "be no information. So neither replaces the other: SAscore "
                    "misses specific impossibilities, and the rules say nothing "
                    "about how hard a route is.")

            show = [c for c in ("rank", "candidate_id", "synth", "SAscore",
                                cfg["metric"],
                                "ligand_efficiency", "dG_kcal",
                                "dG_ensemble_interaction_kcal",
                                "dG_ensemble_interaction_sem_kcal",
                                "dG_ensemble_internal_residual_kcal",
                                "dG_ensemble_kcal", "dG_ensemble_sem_kcal",
                                "warhead_class")
                    if c in s.columns]
            # DEMOTE, DO NOT DELETE (PI decision, issue #1). A compound the Lu
            # lab would not synthesise should not occupy a top-25 slot, so
            # failures sort BELOW every passing candidate while keeping their
            # metric order within each group. They stay visible and keep their
            # original `rank`, because the frame's rank is provenance and the
            # demotion is a reading order -- deleting them would make the list
            # look clean while hiding why it changed.
            s["_synth_fail"] = s["synth"].str.startswith("⚠")
            if demote:
                s = s.sort_values(["_synth_fail", "rank"])
                s["shown_as"] = range(1, len(s) + 1)
                show_cols = ["shown_as"] + [c for c in show if c != "shown_as"]
            else:
                s = s.sort_values("rank")
                show_cols = show
            _t = s[show_cols].head(25)
            st.dataframe(
                _t.style.apply(
                    lambda row: ["background-color: rgba(255,75,75,0.16)"
                                 if str(row.get("synth", "✓")).startswith("⚠")
                                 else "" for _ in row], axis=1),
                use_container_width=True, hide_index=True)
            if demote and n_fail:
                promoted = int((~s["_synth_fail"].head(25)).sum())
                st.caption(
                    f"`shown_as` is the reading order after demotion; `rank` is "
                    f"the original metric rank, untouched. {promoted} candidates "
                    "now occupy the top slots that failures held.")
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
                show_pose_dialog(key, cfg["name"], top)


def panel_dossier() -> None:
    st.header("Per-candidate dossier")
    st.caption("Everything recorded about one candidate, including the SMILES "
               "in a form you can copy into any other tool.")

    c1, c2 = st.columns([1, 2])
    approach = c1.selectbox("approach", list(D.APPROACHES),
                            format_func=lambda k: D.APPROACHES[k]["name"])
    s = D.shortlist(approach)
    if s.empty:
        st.info("no shortlist for this approach yet")
        return
    s = s.sort_values("rank")
    labels = {r["candidate_id"]: f"#{int(r['rank'])} · {r['candidate_id']}"
              for _, r in s.iterrows()}
    cid = c2.selectbox("candidate", list(labels), format_func=lambda k: labels[k])
    row = s[s["candidate_id"] == cid].iloc[0]
    cfg = D.APPROACHES[approach]

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
    cols[1].metric("rank", f"{int(row['rank'])} of {int(row.get('group_n_docked', 0))}")
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

    # Structural convergence: did any OTHER approach reach this molecule?
    # Reported as a rank when the identical molecule was ranked elsewhere, and
    # as the nearest shortlisted analogue otherwise -- an empty panel would let
    # "nobody else found it" read as "the lookup is broken".
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
                    "shortlist pair in the whole build is T_3~T_4 at Tanimoto "
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


def panel_convergence() -> None:
    st.header("Structural convergence")
    st.markdown(
        "A molecule surfaced independently by more than one approach is a soft "
        "cross-validation that relies on **no** score commensurability — the "
        "most defensible cross-approach signal available here.")
    st.caption(
        "Read it with care: all four approaches dock into the same receptor "
        "under related protocols, so their errors are correlated and agreement "
        "may report a shared bias rather than a real signal. T_1 is seed-free, "
        "so its agreement with a seeded approach at least is not ancestral.")

    pool = D.all_shortlists()
    if pool.empty or "canonical_smiles" not in pool.columns:
        st.info("no shortlists to compare yet")
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


def panel_axes() -> None:
    st.header("Shared physicochemical axes")
    st.caption("Computed by the identical RDKit call for all four approaches "
               "(shared/descriptors.py), which is what makes pooling these — "
               "and only these — legitimate.")
    pool = D.all_shortlists()
    if pool.empty:
        st.info("no shortlists yet")
        return
    axes = [a for a in D.SHARED_AXES if a in pool.columns]
    if not axes:
        st.info("no descriptor columns on the frames")
        return
    axis = st.selectbox("axis", axes)
    st.bar_chart(pool.pivot_table(index="approach", values=axis, aggfunc="median"))
    st.dataframe(
        pool.groupby("approach")[axes].median().round(2),
        use_container_width=True)


def panel_within_stratum() -> None:
    st.header("Within-stratum re-score — two leaderboards, never one")
    st.caption("Cross-stratum ordering is not implied and is not offered.")

    fps = D.protocol_fingerprints()
    t3, t4 = fps.get("t3", set()), fps.get("t4", set())
    if t3 and t4 and t3 != t4:
        st.error(
            f"**Within-covalent comparison DISABLED.** T_3 and T_4 recorded "
            f"different protocol fingerprints ({sorted(t3)} vs {sorted(t4)}). "
            "They did not dock under identical rules, so their affinities are "
            "not comparable. Re-dock the lagging approach before comparing — "
            "showing the numbers anyway is exactly what the fingerprint exists "
            "to prevent.")
    else:
        st.success(f"Covalent protocol fingerprints agree: {sorted(t3 or t4)}")

    for label, keys, stratum, metric in (
            ("Non-covalent (T_1 + T_2)", ("t1", "t2"), "non_covalent", "vina_affinity"),
            ("Covalent (T_3 + T_4)", ("t3", "t4"), "covalent", "affinity_kcal")):
        st.subheader(label)
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
        cols = [c for c in ("approach", "candidate_id", metric,
                            "ligand_efficiency", "dG_kcal") if c in pooled.columns]
        st.dataframe(pooled.sort_values(metric)[cols].head(20),
                     use_container_width=True, hide_index=True)


def panel_decisions() -> None:
    st.header("Choreography decision log")
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
    approach = st.selectbox("approach", list(D.APPROACHES),
                            format_func=lambda k: D.APPROACHES[k]["name"])
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
    q = D.open_questions()
    if not q:
        st.success("Nothing is marked proposed, pending or unverified.")
        return
    st.dataframe(pd.DataFrame(q)[["id", "title", "status", "approach"]],
                 use_container_width=True, hide_index=True)


# --------------------------------------------------------------------------

PANELS = {
    "Shortlists": panel_candidates,
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

    # COMPARE EACH MODULE'S IMPORT-TIME MTIME WITH ITS CURRENT ONE. The first
    # version of this compared the helper modules' CURRENT mtimes against the
    # newest file in the directory -- so editing app.py alone made it cry stale,
    # every time, including immediately after a restart. It never knew when
    # anything was imported. Each helper now freezes its own mtime at import
    # (LOADED_MTIME), which is the only value that answers the question.
    _stale = [m.__name__ for m in (D, depict, curate, p3d)
              if getattr(m, "LOADED_MTIME", None) is not None
              and Path(m.__file__).stat().st_mtime > m.LOADED_MTIME + 1]
    if _stale:
        st.sidebar.error(
            f"**Running stale code:** `{'`, `'.join(_stale)}` changed on disk "
            "since this process imported them. Streamlit does not re-import "
            "helper modules — restart the app.")
except Exception:  # noqa: BLE001 - a version badge must never break the page
    pass

st.sidebar.caption("The GUI reads; it does not own (D0008). Everything shown is "
                   "a rendering of files in the repo or the append-only tree.")
PANELS[choice]()
