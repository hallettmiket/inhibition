"""
Purpose: see every pose of one ligand, coloured by contact-space group, and open a group
Author: Timothy Wu (with Claude Code)
Date: 2026-08-26
Input: a persisted pose cloud (raw deep cloud or a production cloud) + the prepared receptor
Output: a Streamlit page; nothing is written

WHY THIS EXISTS. Until now the ranking view could show a mode's representative and
nothing else -- #44's rule (persist every pose) was added precisely because the
cloud behind a mode was being destroyed. Now that the cloud survives and groups
are built from it, the only way to check a grouping is to LOOK at it: a group
whose members are the same pose looks like one molecule drawn twice, and a bag
looks like a bag. D0088's 137-pose group spanning 9.3 A would have been obvious
on screen years before it was measured.

THE RAW CLOUD IS THE DEFAULT SOURCE. `<topic>_allposes` is not all poses -- it
holds only poses whose DBSCAN label survived, so ~21% of it is missing (D0093).
Opening a viewer named "all poses" onto a filtered file is how that defect stayed
invisible, so the source is named on screen and the filtered one is labelled.

THE RECEPTOR IS ALWAYS DRAWN, AS A SURFACE. A ligand rendered alone is a
conformer, not a pose (pose3d.py's founding note, issue #1) -- and a CARTOON is
barely better. pose3d.py says so in as many words: "a spectrum cartoon tells a
reader where the chain runs; it does not tell them where the ligand IS". The
first version of this page defaulted to cartoon with the surface OFF, and
@tt8804 reported poses "literally outside of the pocket". They are not; Pin1's
catalytic site is a shallow surface groove and a ribbon cannot show that. The
surface is now on by default.

ENERGY IS SHOWN, BECAUSE THE ABSENCE OF IT WAS THE REAL DEFECT. The clouds this
page was written against carried no energies at all, so the best-scoring pose and
the 500th were drawn identically. Audited (exp/21): poses with >30% of their atoms
uncontacted are real but rare (2.6%) and sit at the 88th ENERGY PERCENTILE, with
zero of them in the best decile -- the scoring function ranks them correctly and
the viewer simply never said so. A pose set displayed without its scores invites
exactly the conclusion it produced.
"""

from __future__ import annotations

import glob
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st

REPO = Path(__file__).resolve().parents[1]
for p in (REPO, REPO / "exp" / "15_rmsf_predictor", REPO / "exp" / "16_contact_clustering"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from shared import pose_contacts as pc              # noqa: E402
from shared import run_paths as rp                  # noqa: E402

#: Distinct enough to tell apart on a dark surface, and stable across reruns so a
#: group keeps its colour while you click through the table.
PALETTE = [
    "#e6194b", "#3cb44b", "#4363d8", "#f58231", "#911eb4", "#46f0f0",
    "#f032e6", "#bcf60c", "#fabebe", "#008080", "#e6beff", "#9a6324",
    "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000075",
    "#a9a9a9", "#ff4500", "#00ced1", "#7cfc00", "#dc143c", "#1e90ff",
]
GREY = "#5a5a5a"


def _import_by_path(name: str, path: Path):
    """Import by FILE. exp/15 and exp/16 both expose `run_all`."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


_RMSF = _import_by_path("pgv_rmsf", REPO / "exp" / "15_rmsf_predictor" / "run_all.py")
_CONTACT = _import_by_path("pgv_contact", REPO / "exp" / "16_contact_clustering" / "run_all.py")


# --------------------------------------------------------------------------- #
#  sources
# --------------------------------------------------------------------------- #
def clouds() -> list[tuple[str, str, Path]]:
    """(label, ident, path) for every readable pose cloud, raw ones first.

    Read access is governed by an ACL the client cannot see, so a cloud is listed
    only if `os.access` says it can actually be opened -- an unreadable file
    presented as available is the seed_status confusion in a new place.
    """
    out = []
    # RAW CLOUDS FIRST: these are the only ones carrying energies, and a page
    # that opens by default on an energy-less cloud is the defect this page had.
    for d in sorted(rp.BLACKSMITH.glob("raw_cloud_*")):
        fs = sorted(d.glob("cloud_*.sdf"), key=os.path.getmtime)
        if fs and os.access(fs[-1], os.R_OK):
            ident = d.name[len("raw_cloud_"):]
            out.append((f"RAW + energies — {ident}", ident, fs[-1]))
    for d in sorted(rp.BLACKSMITH.glob("deep_cloud_*")):
        for f in sorted(d.glob("cloud_*.sdf")):
            if os.access(f, os.R_OK):
                out.append((f"RAW deep cloud, no energies — {d.name[11:]}",
                            d.name[11:], f))
    ap = rp.allposes_dir()
    if ap.is_dir():
        for f in sorted(ap.glob("*.sdf")):
            if os.access(f, os.R_OK):
                out.append((f"mode-assigned only (D0093) — {f.stem}", f.stem, f))
    return out


@st.cache_data(show_spinner=False)
def load_cloud(path_str: str, mtime: float):
    """Heavy-atom coordinates and the template molecule.

    KEYED ON THE FILE AND ITS MTIME, not on the ident. A cache keyed on less than
    its inputs is catalogue entry #18; a rewritten cloud must invalidate this.
    """
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    ms = [m for m in Chem.SDMolSupplier(path_str, removeHs=False, sanitize=True)
          if m is not None]
    if not ms:
        return None, None, None
    hv = [a.GetIdx() for a in ms[0].GetAtoms() if a.GetAtomicNum() > 1]
    xyz = np.array([m.GetConformer().GetPositions()[hv] for m in ms])
    sym = [ms[0].GetAtomWithIdx(i).GetSymbol() for i in hv]
    # NaN, not a default. A cloud written before energies were persisted has no
    # score, and inventing one (0.0, or the median) would put a number on screen
    # that no docking produced.
    en = np.array([float(m.GetProp("free_energy_kcal"))
                   if m.HasProp("free_energy_kcal") else np.nan for m in ms])
    return xyz, sym, ms[0], en


@st.cache_data(show_spinner=False)
def residue_landmarks(n_res: int):
    names = _CONTACT.key_residues(n_res)
    return names, _CONTACT.receptor_coords(names)


@st.cache_data(show_spinner=False)
def rmsf_for(path_str: str, mtime: float, n_conf: int, seed: int):
    xyz, sym, tmpl, _en = load_cloud(path_str, mtime)
    hv = [a.GetIdx() for a in tmpl.GetAtoms() if a.GetAtomicNum() > 1]
    return _RMSF.predict_rmsf(tmpl, hv, n_conf, seed)


@st.cache_data(show_spinner=False)
def grouping(path_str: str, mtime: float, n_res: int, tol: float,
             n_conf: int, seed: int, max_poses: int, keep_frac: float = 1.0):
    """Labels, medoid index per group, and a per-group summary table."""
    xyz, sym, _t, en = load_cloud(path_str, mtime)
    # ENERGY FILTER FIRST, so the grouping describes the poses on screen.
    if keep_frac < 1.0 and np.isfinite(en).any():
        k = max(2, int(round(len(xyz) * keep_frac)))
        best = np.argsort(np.where(np.isfinite(en), en, np.inf))[:k]
        best.sort()
        xyz, en = xyz[best], en[best]
    else:
        best = np.arange(len(xyz))
    if len(xyz) > max_poses:
        idx = np.random.default_rng(seed).choice(len(xyz), max_poses, replace=False)
        idx.sort()
        xyz, en = xyz[idx], en[idx]
    else:
        idx = np.arange(len(xyz))
    _, res = residue_landmarks(n_res)
    rmsf = rmsf_for(path_str, mtime, n_conf, seed)
    w = pc.atom_weights(rmsf)
    D = pc.pose_distances(pc.contact_tensor(xyz, res), w)
    lab = pc.group(D, tol)
    rows, med = [], []
    for k in range(lab.max() + 1):
        mem = np.flatnonzero(lab == k)
        sub = D[np.ix_(mem, mem)]
        m = int(mem[sub.sum(1).argmin()])
        med.append(m)
        c = xyz[mem]
        if len(mem) > 1:
            s = c if len(mem) <= 120 else c[:120]
            dd = np.array([np.sqrt(((s - s[i]) ** 2).sum(-1).mean(-1))
                           for i in range(len(s))])
            iu = np.triu_indices(len(s), 1)
            cm, cx = float(np.median(dd[iu])), float(dd[iu].max())
        else:
            cm = cx = 0.0
        rows.append(dict(group=k, poses=len(mem),
                         best_energy=float(np.nanmin(en[mem])) if np.isfinite(en[mem]).any() else np.nan,
                         contact_width=float(sub.max()),
                         rmsd_median=cm, rmsd_max=cx))
    t = pd.DataFrame(rows).sort_values("poses", ascending=False).reset_index(drop=True)
    return xyz, sym, lab, np.array(med), t, idx, en


# --------------------------------------------------------------------------- #
#  rendering
# --------------------------------------------------------------------------- #
def receptor_block() -> str:
    return rp.receptor_prep().read_text(encoding="utf-8", errors="replace")


def xyz_block(coords: np.ndarray, sym: list) -> str:
    """One pose as an XYZ record — no bonds, so RDKit never has to re-perceive."""
    lines = [str(len(coords)), "pose"]
    lines += [f"{s} {c[0]:.3f} {c[1]:.3f} {c[2]:.3f}" for s, c in zip(sym, coords)]
    return "\n".join(lines)


def render(xyz, sym, draw: list[tuple[int, str]], style: str, opacity: float,
           surface: bool, height: int) -> str:
    """`draw` is (pose index, colour). Receptor first, ligands into it."""
    import py3Dmol
    v = py3Dmol.view(width="100%", height=height)
    v.addModel(receptor_block(), "pdb")
    v.setStyle({"model": 0}, {"cartoon": {"color": "#d8d8d8", "opacity": 0.55}})
    if surface:
        v.addSurface("VDW", {"opacity": opacity, "color": "#b9c6d4"}, {"model": 0})
    for i, colour in draw:
        v.addModel(xyz_block(xyz[i], sym), "xyz")
        n = v.getModel().__dict__.get("id", None)
        spec = ({"stick": {"colorscheme": {"prop": "elem", "map": {}},
                           "color": colour, "radius": 0.13}}
                if style == "stick" else
                {"line": {"color": colour, "linewidth": 2.0}})
        v.setStyle({"model": -1}, spec)
    v.zoomTo({"model": 0})
    return v._make_html()


# --------------------------------------------------------------------------- #
def main() -> None:
    st.set_page_config(page_title="Pose groups", layout="wide")
    st.title("Pose groups — every pose of one ligand, coloured by contact-space group")

    cs = clouds()
    if not cs:
        st.error("No readable pose cloud found under "
                 f"`{rp.BLACKSMITH}`. Read access is governed by an Isilon ACL the "
                 "client cannot see — check `test -r` on a cloud SDF before "
                 "debugging this page.")
        return

    with st.sidebar:
        st.header("Source")
        labels = [c[0] for c in cs]
        pick = st.selectbox("pose cloud", range(len(cs)), format_func=lambda i: labels[i])
        label, ident, path = cs[pick]
        if "no energies" in label:
            st.warning("This cloud has no per-pose energies. Prefer a "
                       "**RAW + energies** entry so the score filter works.")
        if label.startswith("mode-assigned"):
            st.warning("This file holds only poses whose DBSCAN label survived — "
                       "~21% of the cloud is absent (D0093). Prefer a RAW cloud.")
        mtime = os.path.getmtime(path)

        st.header("Grouping")
        n_res = st.slider("landmark residues", 5, 25, 15, help="exp/14's greedy pick")
        max_poses = st.select_slider("poses to load", [500, 1000, 2000, 4000, 6000],
                                     value=2000)
        rmsf = rmsf_for(str(path), mtime, 50, 7)
        auto = float(np.median(rmsf) / pc.RMSF_CALIBRATION)
        use_auto = st.checkbox(f"tolerance from predicted RMSF ({auto:.2f} Å)", True)
        tol = auto if use_auto else st.slider("tolerance (Å)", 0.3, 3.5, auto, 0.05)
        st.caption("The tolerance is NOT meaningfully per-molecule — the ensemble "
                   "ranks atoms within a molecule, not molecules against each "
                   "other (D0094).")

        st.header("Energy")
        _x, _s, _t2, _en_probe = load_cloud(str(path), mtime)
        has_energy = bool(np.isfinite(_en_probe).any())
        if has_energy:
            keep_pct = st.select_slider(
                "keep the best N% by docking energy",
                [10, 25, 50, 75, 100], value=100,
                help="Poses far down the energy list are real docking output but "
                     "are NOT what the score favours. exp/21: poses with >30% of "
                     "atoms uncontacted sit at the 88th energy percentile.")
            keep_frac = keep_pct / 100.0
        else:
            keep_frac = 1.0
            st.warning("This cloud carries **no energies** — it was written before "
                       "`persist_raw_clouds.py` recorded them. Every pose is drawn "
                       "identically regardless of score. Re-persist it with "
                       "`--force` to filter here.")

        st.header("Display")
        style = st.radio("ligand style", ["line", "stick"], horizontal=True)
        # ON BY DEFAULT. A ribbon cannot show a shallow surface groove, and the
        # first version of this page defaulted it off -- see the module docstring.
        surface = st.checkbox("pocket surface (shows where the pocket IS)", True)
        opacity = st.slider("surface opacity", 0.1, 1.0, 0.65, 0.05) if surface else 0.65
        height = st.slider("viewer height (px)", 400, 1100, 700, 50)

    with st.spinner("grouping poses…"):
        xyz, sym, lab, med, table, _, energy = grouping(
            str(path), mtime, n_res, tol, 50, 7, max_poses, keep_frac)

    n_g = len(table)
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("poses", f"{len(xyz):,}")
    c2.metric("groups", f"{n_g:,}")
    c3.metric("largest group", f"{int(table.poses.max()):,}")
    c4.metric("tolerance", f"{tol:.2f} Å")
    if np.isfinite(energy).any():
        st.caption(
            f"**Docking energy** best {np.nanmin(energy):.2f}, median "
            f"{np.nanmedian(energy):.2f}, worst {np.nanmax(energy):.2f} kcal/mol"
            + (f" · showing the best {keep_pct}%" if keep_frac < 1 else "")
            + f" · {np.isfinite(energy).mean() * 100:.0f}% of poses scored.")
    st.caption(
        f"`{ident}` · {label} · {n_res} landmark residues · "
        f"{int((table.poses == 1).sum())} singleton groups "
        f"({(table.poses == 1).mean() * 100:.0f}%). "
        "**The group count is a function of docking depth and is not a number of "
        "binding modes** (D0092) — change *poses to load* and watch it move.")

    left, right = st.columns([2, 1], gap="medium")

    with right:
        st.subheader("Groups")
        st.caption("Ranked by size. Pick one or more to see the poses inside them; "
                   "leave empty to see one representative per group.")
        top_n = st.slider("groups to colour in the overview", 4, 24, 12)
        opts = table.group.tolist()
        chosen = st.multiselect(
            "open group(s)", opts, default=[],
            format_func=lambda g: (f"#{g} — {int(table.loc[table.group == g, 'poses'].iloc[0])} poses, "
                                   f"{table.loc[table.group == g, 'rmsd_max'].iloc[0]:.2f} Å wide"))
        show = table.copy()
        show["colour"] = [PALETTE[i % len(PALETTE)] if i < top_n else GREY
                          for i in range(len(show))]
        st.dataframe(
            show.head(200).style.apply(
                lambda r: [f"background-color: {r.colour}33"] * len(r), axis=1),
            column_config={
                "group": "group", "poses": "poses",
                "best_energy": st.column_config.NumberColumn("best kcal", format="%.2f"),
                "contact_width": st.column_config.NumberColumn("contact Å", format="%.2f"),
                "rmsd_median": st.column_config.NumberColumn("RMSD med", format="%.2f"),
                "rmsd_max": st.column_config.NumberColumn("RMSD max", format="%.2f"),
                "colour": None},
            hide_index=True, use_container_width=True, height=340)
        st.caption("`contact Å` is bounded by the tolerance **structurally** — "
                   "complete linkage defines group distance by the farthest pair. "
                   "`RMSD max` is the Cartesian width, which is *not* guaranteed "
                   "and is the number worth watching.")

    with left:
        colour_of = {g: (PALETTE[i % len(PALETTE)] if i < top_n else GREY)
                     for i, g in enumerate(table.group)}
        if chosen:
            draw = []
            for g in chosen:
                c = colour_of.get(g, GREY)
                for i in np.flatnonzero(lab == g):
                    draw.append((int(i), c))
            st.subheader(f"{len(chosen)} group(s) — {len(draw)} poses")
        else:
            draw = [(int(med[g]), colour_of[g]) for g in table.group.head(top_n)]
            draw += [(int(med[g]), GREY) for g in table.group.iloc[top_n:]]
            st.subheader(f"Overview — one representative per group ({len(draw)})")

        cap = 600
        if len(draw) > cap:
            st.warning(f"{len(draw):,} poses selected; drawing the first {cap:,}. "
                       "Narrow the selection or lower *poses to load*.")
            draw = draw[:cap]
        st.components.v1.html(
            render(xyz, sym, draw, style, opacity, surface, height),
            height=height + 24)

        if chosen:
            sub = table[table.group.isin(chosen)]
            st.markdown("**Selected groups**")
            for _, r in sub.iterrows():
                st.markdown(
                    f"<span style='color:{colour_of.get(r.group, GREY)};font-size:20px'>&#9632;</span> "
                    f"**#{int(r.group)}** — {int(r.poses)} poses · best "
                    f"{r.best_energy:.2f} kcal/mol · contact width "
                    f"{r.contact_width:.2f} Å (≤ {tol:.2f}) · Cartesian RMSD "
                    f"median {r.rmsd_median:.2f} Å, max {r.rmsd_max:.2f} Å",
                    unsafe_allow_html=True)


if __name__ == "__main__":
    main()
