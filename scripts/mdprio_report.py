"""
Purpose: one 100 ns MD-priority molecule -> one self-contained report, the moment it lands.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-07
Input: --candidate <ident> (its trajectory directory under md_residence_3ikd/)
Output: 00_outputs/blacksmith/mdprio_reports/<ident>.html

WHY THIS EXISTS. `md_residence_3ikd.py` already runs its gmx analysis inline and
writes that molecule's row as soon as it finishes, so the NUMBERS were never
batched. The REPORT was: `elevation_report.py` is hard-wired to the 2.1.0
elevation experiment (tier-1 cohort, tier-2 replicas, one lead), so there was no
way to look at a single MD-priority molecule without waiting for all six and
hand-running an analysis. This closes that gap -- six molecules that finish hours
apart produce six reports hours apart.

THE PRE-REGISTERED VERDICT IS DELIBERATELY NOT COMPUTED HERE. `docs/prereg_md_priority.md`
fixes the readout as a RANK CORRELATION between BPMD occupancy and 100 ns
residence across all six. A rank correlation over a partial set is not a
preliminary version of the final number -- it is a different number, and quoting
it early is how a null becomes a trend. Each report states its own molecule's
residence and where that sits against the BPMD prediction, and stops there.

UNITS ARE READ, NOT ASSUMED. `gmx rms`/`mindist`/`numcont` write nanoseconds and
`gmx distance` writes picoseconds. Dividing everything by 1000 once squashed
three series into the first 0.1 ns of a 100 ns plot, and the panel looked EMPTY
rather than wrong. `to_ns` is imported from `elevation_report` rather than
rewritten, so there is one rule for this.
"""

from __future__ import annotations

import argparse
import importlib.util
import logging
import sys
from pathlib import Path

import numpy as np
import re
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import outputs as sout                  # noqa: E402
from shared import report_theme as rt               # noqa: E402
from shared import md_movie as mov                  # noqa: E402

log = logging.getLogger("mdprio-report")
MD = Path("/data/lab_vm/modifiable/inhibition/md_residence_3ikd")
OUT = sout.Topic("blacksmith", "mdprio_reports")
BOUND_NM = 1.0            # the residence criterion: ligand RMSD <= 1.0 nm

#: BPMD occupancy each molecule was selected on (docs/prereg_md_priority.md).
#: Carried here so a report can state what was PREDICTED for this molecule
#: without re-deriving it, and without pooling across molecules.
PREDICTED = {
    "t4_da2e98512d02": (0.365, "bdhi_c5", 1),
    "t4_7e86b677bb2d": (0.189, "acrylamide", 6),
    "t4_9a973be6b946": (0.161, "bdhi_c4", 2),
    "t4_28f5ea16adeb": (0.152, "acrylamide", 1),
    "t4_4e608398fd6a": (0.125, "bdhi_c4", 1),
    "t4_9265b4bff789": (0.108, "acrylamide", 8),
}
REF_MEDIAN = 0.163        # crystallographic BPMD median — a yardstick, not a control


def nac_lo() -> float:
    import shared.nac_criterion as nac
    return nac.NAC_DIST_MIN


def nac_hi() -> float:
    import shared.nac_criterion as nac
    return nac.NAC_DIST_MAX


def _er():
    spec = importlib.util.spec_from_file_location(
        "elevation_report", REPO / "scripts" / "elevation_report.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def series(rep: Path, er) -> dict:
    """Every gmx series this run produced, each on its own true time axis."""
    out = {}
    for key, fname in (("rmsd", "rmsd.xvg"), ("mindist", "mindist.xvg"),
                       ("contacts", "numcont.xvg"), ("dist", "dist.xvg")):
        p = rep / fname
        if not p.is_file():
            continue
        a = er.read_xvg(p)
        if a.size == 0:
            continue
        out[key] = (er.to_ns(a[:, 0]), a[:, 1])
    return out


def residence(s: dict) -> dict:
    """Residence fraction and whether the ligand ever left.

    The readout the pre-registration names: the fraction of frames with ligand
    RMSD <= 1.0 nm, plus whether dissociation happened at all. A single
    dissociation event is ONE draw -- reported as a screen, never as a rate.
    """
    if "rmsd" not in s:
        return {"status": "no rmsd.xvg — run incomplete"}
    t, r = s["rmsd"]
    bound = r <= BOUND_NM
    left = None
    if (~bound).any():
        # First frame after which it never comes back within the criterion.
        idx = np.where(~bound)[0]
        for i in idx:
            if not bound[i:].any():
                left = float(t[i])
                break
    return {"status": "ok", "n_frames": int(len(r)),
            "length_ns": float(t[-1]), "residence_frac": float(bound.mean()),
            "rmsd_mean_nm": float(r.mean()), "rmsd_max_nm": float(r.max()),
            "rmsd_final_nm": float(r[-1]),
            "left_at_ns": left, "dissociated": left is not None}


def nac_series(cand: str, rep: Path, movie: Path) -> dict | None:
    """Warhead->Cys113 SG distance AND the near-attack angle, per frame.

    NEITHER IS PRODUCED BY THE RUN. `md_residence` writes rmsd/mindist/numcont
    only: `mindist` is the closest ligand-protein approach by ANY atom pair,
    which is not the warhead and not the sulfur. So the two quantities the whole
    project is about were never plotted for a 100 ns trajectory.

    Computed here from the fitted movie frames, with the reactive atom located by
    the SAME SMARTS the screen uses and the geometry measured by the SAME
    `nac_criterion.measure`. A second definition of "the attack angle" is how the
    plot and the ranking would come to disagree while both looked right.

    Atom correspondence: the MD ligand was parameterised FROM `ligand_pose.sdf`,
    so the MOL residue's atoms are that file's atoms in order. Hydrogens are
    stripped from both sides before pairing, because the movie carries heavy
    atoms only.
    """
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    import shared.nac_criterion as nac

    sdf = rep.parent.parent / "ligand_pose.sdf"
    if not (sdf.is_file() and movie.is_file()):
        return None
    mol = Chem.SDMolSupplier(str(sdf), removeHs=False)[0]
    if mol is None:
        return None
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    pos = {old: new for new, old in enumerate(heavy)}     # sdf idx -> heavy idx

    wh = pd.read_csv(REPO / "data/reference/warhead_classes_10.csv")
    match = mech = None
    for r in wh.itertuples():
        patt = Chem.MolFromSmarts(r.reactive_atom_smarts)
        if patt is None:
            continue
        ms = mol.GetSubstructMatches(patt)
        if ms:
            match, mech = ms[0], r.mechanism
            break
    if match is None or any(i not in pos for i in match):
        log.warning("%s: no reactive SMARTS match on the MD ligand", cand)
        return None

    # walk the movie: MOL heavy atoms in order, and Cys113's SG
    frames_lig, frames_sg, cur_lig, sg = [], [], [], None
    for ln in movie.read_text().splitlines():
        if ln.startswith("MODEL"):
            cur_lig, sg = [], None
        elif ln.startswith("ENDMDL"):
            if cur_lig and sg is not None:
                frames_lig.append(np.array(cur_lig)); frames_sg.append(sg)
        elif ln.startswith(("ATOM", "HETATM")):
            resn, name = ln[17:20].strip(), ln[12:16].strip()
            xyz = [float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]
            if resn == "MOL":
                cur_lig.append(xyz)
            elif resn == "CYS" and name == "SG":
                sg = np.array(xyz)
    if not frames_lig:
        return None

    idx = [pos[i] for i in match]
    dist, ang = [], []
    for lig, s_ in zip(frames_lig, frames_sg):
        if max(idx) >= len(lig):
            return None
        r = nac.measure(mech, lig[idx], s_)
        dist.append(r.distance); ang.append(r.angle)
    t = np.linspace(0.0, 100.0, len(dist))
    return {"t": t, "dist": np.array(dist), "angle": np.array(ang),
            "kind": nac.MECHANISMS.get(mech, ""), "mechanism": mech}


def figure(ident: str, s: dict, res: dict, er, nacs: dict | None = None) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update(rt.MPL)

    panels = [k for k in ("rmsd", "mindist", "contacts") if k in s]
    extra = ["nac_dist", "nac_angle"] if nacs else []
    panels = panels[:1] + extra + panels[1:]
    fig, axes = plt.subplots(len(panels), 1, figsize=(9, 2.2 * len(panels)),
                             sharex=True)
    axes = np.atleast_1d(axes)
    titles = {"rmsd": "Ligand RMSD after superposing on the protein (nm)",
              "nac_dist": "Warhead → Cys113 SG distance (Å)",
              "nac_angle": "Near-attack angle (°)",
              "mindist": "Minimum ligand–protein distance (nm)",
              "contacts": "Ligand–protein contacts"}
    for ax, k in zip(axes, panels):
        if k == "nac_dist":
            t, y = nacs["t"], nacs["dist"]
        elif k == "nac_angle":
            t, y = nacs["t"], nacs["angle"]
        else:
            t, y = s[k]
        ax.plot(t, y, lw=0.7, color=rt.SERIES["accent"])
        ax.set_title(titles[k], loc="left", fontsize=9)
        ax.margins(x=0)
        if k == "rmsd":
            ax.axhline(BOUND_NM, ls="--", lw=1, color=rt.SERIES["alert"])
            ax.text(0.995, BOUND_NM, f" bound ≤ {BOUND_NM} nm", va="bottom",
                    ha="right", transform=ax.get_yaxis_transform(),
                    fontsize=8, color=rt.SERIES["alert"])
            if res.get("left_at_ns"):
                ax.axvline(res["left_at_ns"], lw=1, color=rt.SERIES["alert"])
        if k == "nac_dist":
            # The attack window, so "close enough to react" is on the plot
            # rather than in the reader's head.
            ax.axhspan(nac_lo(), nac_hi(), color=rt.SERIES["ref"], alpha=0.13)
            ax.text(0.995, nac_hi(), " attack window", va="bottom", ha="right",
                    transform=ax.get_yaxis_transform(), fontsize=8,
                    color=rt.SERIES["ref"])
        if k == "nac_angle":
            # The competent band depends on the MECHANISM, and getting it wrong
            # would draw a target the molecule was never aiming at: SN2 wants a
            # backside approach at >=150 deg, the perpendicular mechanisms want
            # <=30 deg off the sp2 plane normal. Read from the criterion rather
            # than hardcoded per panel.
            import shared.nac_criterion as nacm
            if "anti" in (nacs.get("kind") or ""):
                ax.axhspan(nacm.SN2_ANGLE_MIN, 180, color=rt.SERIES["ref"], alpha=0.13)
                lbl, at = f" ≥{nacm.SN2_ANGLE_MIN:.0f}° backside", nacm.SN2_ANGLE_MIN
            else:
                ax.axhspan(0, nacm.PERPENDICULAR_MAX_OFF_NORMAL,
                           color=rt.SERIES["ref"], alpha=0.13)
                lbl = f" ≤{nacm.PERPENDICULAR_MAX_OFF_NORMAL:.0f}° off plane normal"
                at = nacm.PERPENDICULAR_MAX_OFF_NORMAL
            ax.text(0.995, at, lbl, va="bottom", ha="right",
                    transform=ax.get_yaxis_transform(), fontsize=8,
                    color=rt.SERIES["ref"])
    axes[-1].set_xlabel("time (ns)")
    fig.suptitle(ident, x=0.005, ha="left", fontsize=11, weight="bold")
    fig.tight_layout()
    return er._png(fig, plt)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--candidate", required=True)
    ap.add_argument("--rep", default="rep1")
    ap.add_argument("--out", default=None)
    ap.add_argument("--no-movie", action="store_true",
                    help="skip the 3D trajectory movie (it costs a trjconv pass)")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    er = _er()
    rep = MD / args.candidate / "md" / args.rep
    if not rep.is_dir():
        raise SystemExit(f"no trajectory directory {rep}")

    s = series(rep, er)
    res = residence(s)
    if res["status"] != "ok":
        raise SystemExit(f"{args.candidate}: {res['status']}")
    occ, cls, pose = PREDICTED.get(args.candidate, (None, None, None))

    # PREDICTED only covers the six BPMD-priority molecules. Everything swept
    # since then fell through to "—" for class and pose, so the masthead of every
    # newer report advertised missing data it could have looked up. Class comes
    # from the ranking, the sweep readings from the sweep, and the SMILES gives
    # the depiction -- all keyed on the molecule, all already on disk.
    import glob as _g
    smiles = None
    for _sub, _stem in (("04_t4_combinatorial", "D4"), ("03_t3_reinvent", "D3")):
        _fs = sorted(_g.glob(f"/data/lab_vm/append_only/inhibition/{_sub}/{_stem}_*.parquet"),
                     key=lambda q: int(q.rsplit("_", 1)[1].split(".")[0]))
        if not _fs:
            continue
        _fr = pd.read_parquet(_fs[-1]).drop_duplicates("candidate_id").set_index("candidate_id")
        if args.candidate in _fr.index:
            smiles = _fr.loc[args.candidate].get("canonical_smiles")
            if not cls:
                cls = _fr.loc[args.candidate].get("warhead_class")
            break
    if not cls:
        for _t, _sc in (("T4", "conditional_eb"), ("T3", "enrichment_conditional")):
            _fs = sorted(_g.glob("/data/lab_vm/append_only/inhibition/00_outputs/"
                                 f"blacksmith/rank_v2/rank_v2_{_t}_{_sc}_*.csv"))
            if not _fs:
                continue
            _rk = pd.read_csv(_fs[-1]).drop_duplicates("parent_ident").set_index("parent_ident")
            if args.candidate in _rk.index:
                cls = _rk.loc[args.candidate, "warhead_class"]; break
    sweep_ar = sweep_v = None
    _fs = sorted(_g.glob("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith/"
                         "attack_sweep/attack_sweep_*.csv"),
                 key=lambda q: int(q.rsplit("_", 1)[1].split(".")[0]))
    if _fs:
        _sw = pd.concat([pd.read_csv(f) for f in _fs], ignore_index=True)
        _sw = _sw[(_sw.get("sweep_ps", 0) > 1000) & (_sw.status == "ok")
                  & (_sw.parent_ident == args.candidate)]
        if len(_sw):
            _b = _sw.sort_values("frac_attack_ready").iloc[-1]
            sweep_ar, sweep_v = float(_b.frac_attack_ready), int(_b.n_visits)

    struct_svg = ""
    if isinstance(smiles, str):
        try:
            from rdkit import Chem as _C, RDLogger as _R
            from rdkit.Chem import Draw as _D, AllChem as _A
            _R.DisableLog("rdApp.*")
            _m = _C.MolFromSmiles(smiles)
            if _m is not None:
                _A.Compute2DCoords(_m)
                _d = _D.rdMolDraw2D.MolDraw2DSVG(300, 190)
                _D.rdMolDraw2D.PrepareAndDrawMolecule(_d, _m); _d.FinishDrawing()
                _svg = re.sub(r"<\?xml.*?\?>", "", _d.GetDrawingText(), flags=re.S)
                struct_svg = re.sub(r"<!--.*?-->", "", _svg, flags=re.S)
        except Exception:                                  # noqa: BLE001
            struct_svg = ""

    log.info("%s: %.1f ns, residence %.3f, dissociated=%s", args.candidate,
             res["length_ns"], res["residence_frac"], res["dissociated"])

    # ---- the movie: surface, charge colouring, labelled key residues -------
    # Built from the SAME trajectory the numbers come from, so the panel and the
    # structure cannot describe different runs.
    movie_block = ""
    mpdb = rep / "movie.pdb"
    if not args.no_movie:
        if not mpdb.is_file():
            mov.build_movie_pdb(rep, mpdb, total_ps=res["length_ns"] * 1000.0)
        if mpdb.is_file():
            try:
                pdb_txt, dsg, labels, lpos = er.surface_payload(mpdb)
                three = (REPO / "scripts/.cache_3dmol-min.js").read_text()
                movie_block = mov.viewer_html(pdb_txt, dsg, labels, lpos, three)
                log.info("movie embedded: %d frames", len(dsg))
            except Exception as exc:                   # noqa: BLE001
                # A failed movie must not cost the numbers. Recorded in the page
                # rather than dropped, so a missing viewer is visible as a
                # failure and not mistaken for "this run had no movie".
                log.warning("movie failed: %s", exc)
                movie_block = rt.callout(
                    "Movie unavailable",
                    f"The trajectory rendered no viewer: <code>{exc}</code>. "
                    "The figures and residence numbers above are unaffected — "
                    "they are computed from the trajectory directly.", "warn")
    verdict = ("Held" if not res["dissociated"] else
               f"Left at {res['left_at_ns']:.0f} ns")

    # The NAC series needs the fitted movie frames, so it comes after the movie.
    nacs = nac_series(args.candidate, rep, mpdb) if not args.no_movie else None
    if nacs is None:
        log.warning("no warhead->SG distance/angle series for %s", args.candidate)
    img = figure(args.candidate, s, res, er, nacs)

    facts = [("molecule", args.candidate), ("warhead class", cls or "unclassified")]
    if sweep_ar is not None:
        facts.append(("attack-ready (10 ns)", f"{sweep_ar*100:.1f}%  ·  {sweep_v} visits"))
    if occ is not None:
        facts.append(("BPMD occupancy", rt.num(occ, "{:.3f}")))
    if pose:
        facts.append(("pose elevated", f"rank {pose}"))
    facts.append(("trajectory", f"{res['length_ns']:.1f} ns, {res['n_frames']:,} frames"))
    body = [
        rt.masthead(
            f"{args.candidate} — 100 ns residence",
            f"{verdict}. Residence fraction {res['residence_frac']:.3f} "
            f"(frames with ligand RMSD ≤ {BOUND_NM} nm).",
            "MD-PRIORITY · 2.2.0 BORNITE", facts),
        (f'<div class="structrow"><div class="structbox">{struct_svg}</div>'
         f'<div class="structnote"><b>{cls or "unclassified"}</b>'
         + (f' &middot; attack-ready {sweep_ar*100:.1f}% of the 10 ns sweep'
            f' over {sweep_v} sustained visits' if sweep_ar is not None else '')
         + '</div></div>') if struct_svg else "",
        f'<p>{rt.pill("Held" if not res["dissociated"] else "Left")} '
        f'Mean RMSD {res["rmsd_mean_nm"]:.3f} nm &middot; max '
        f'{res["rmsd_max_nm"]:.3f} nm &middot; final {res["rmsd_final_nm"]:.3f} nm.</p>',
        f'<details class="panel"><summary>RMSD plots'
        f'<span class="hint">RMSD, warhead&ndash;Cys113 distance, attack angle</span>'
        f'</summary><div class="pbody">'
        f'<img src="data:image/png;base64,{img}" alt="trajectory plots"></div></details>',
        (f'<details class="panel"><summary>MD movie'
         f'<span class="hint">surface by charge, ligand in yellow, CA-fitted</span>'
         f'</summary><div class="pbody">{movie_block}</div></details>')
        if movie_block else "",
        '<details class="panel"><summary>How it was selected, and what that is worth'
        '<span class="hint">pre-registration context</span></summary>'
        '<div class="pbody">',
        f"<p>It was selected on a BPMD occupancy of "
        f"<strong>{rt.num(occ, '{:.3f}')}</strong>, against a crystallographic "
        f"median of {REF_MEDIAN:.3f}. That number is what the pre-registration "
        f"is testing — it is <em>not</em> evidence about this molecule.</p>",
        rt.callout(
            "One molecule cannot answer the question",
            "<code>docs/prereg_md_priority.md</code> fixes the readout as a rank "
            "correlation between BPMD occupancy and 100 ns residence across all "
            "six molecules. A correlation computed on the subset that has "
            "finished is not an early version of that number — it is a different "
            "number, and reading it as a trend is exactly how a null becomes a "
            "result. This report deliberately stops at one molecule. "
            "<strong>n = 6 supports only large effects</strong>, and a single "
            "dissociation event is one draw with ~100% relative standard error, "
            "so residence is a screen and not a rate.",
            "warn"),
        "</div></details>",
    ]
    dest = Path(args.out) if args.out else OUT.dir / f"{args.candidate}.html"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(
        "<!doctype html><html lang='en'><head><meta charset='utf-8'>"
        "<meta name='viewport' content='width=device-width,initial-scale=1'>"
        f"<title>{args.candidate} — 100 ns residence</title>"
        f"<style>{rt.CSS}{mov.VIEWER_CSS}</style></head><body>\n"
        + "\n".join(body) + "\n</body></html>")
    print(f"{args.candidate}: {verdict}, residence {res['residence_frac']:.3f} "
          f"-> {dest}")

    row = {"ident": args.candidate, "bpmd_occupancy": occ, "warhead_class": cls,
           "pose_rank": pose, **{k: v for k, v in res.items() if k != "status"}}
    pd.DataFrame([row]).to_csv(OUT.write(f"mdprio_{args.candidate}", ".csv"),
                               index=False)


if __name__ == "__main__":
    main()
