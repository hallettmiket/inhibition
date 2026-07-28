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

EVERY RANK CARRIES ITS GATE (D0031). On class-matched decoys the covalent gate
reads ROC-AUC 0.537 and the non-covalent 0.535 — both indistinguishable from
chance. So no shortlist here is evidence that the molecules at the top bind, and
a rank displayed without that verdict would imply a confidence nothing supports.
The verdict is shown beside every ranking, not buried in a methods note.

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

APP_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(APP_DIR))
sys.path.insert(0, str(APP_DIR.parent.parent))

import data as D                                  # noqa: E402

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
        "the four metrics are different quantities. On class-matched decoys "
        "docking is indistinguishable from chance on this receptor (covalent "
        "ROC-AUC 0.537, non-covalent 0.535, EF1% 0.0 for both; D0031), so no "
        "shortlist is evidence of binding. Inhibition versus activation is "
        "unresolved, and there is no wet-lab ground truth for any candidate.")


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

def panel_candidates() -> None:
    st.header("Shortlists — four approaches, side by side")
    st.caption("Each column is ranked by ITS OWN metric. The columns are not "
               "comparable with each other and are deliberately not merged.")
    honest_limits()

    cols = st.columns(len(D.APPROACHES))
    for col, (key, cfg) in zip(cols, D.APPROACHES.items()):
        with col:
            st.subheader(cfg["name"])
            s = D.shortlist(key)
            if s.empty:
                st.info("no shortlist yet")
                continue
            verdict = str(s["gate_verdict"].iloc[0]) if "gate_verdict" in s else "UNGATED"
            st.markdown(
                f"{gate_badge(verdict)} **{verdict}** · metric "
                f"`{cfg['metric']}` (lower better) · {cfg['mechanism']} · "
                f"seed: {cfg['seed']}")
            show = [c for c in ("rank", "candidate_id", cfg["metric"],
                                "ligand_efficiency", "dG_kcal", "warhead_class")
                    if c in s.columns]
            st.dataframe(s.sort_values("rank")[show].head(25),
                         use_container_width=True, hide_index=True)
            if "dG_kcal" in s.columns and s["dG_kcal"].notna().any():
                st.caption(f"MM-GBSA dG on {int(s['dG_kcal'].notna().sum())} of "
                           f"{len(s)} — an INDEPENDENT estimate, not confirmation "
                           "of the docking rank.")


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
st.sidebar.caption("The GUI reads; it does not own (D0008). Everything shown is "
                   "a rendering of files in the repo or the append-only tree.")
PANELS[choice]()
