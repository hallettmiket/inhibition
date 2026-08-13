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
import json
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


def prod_ns(rep: Path, default: float = 100.0) -> float:
    """The production length of THIS run, from the mdp that produced it.

    `elevation_report.to_ns` rescales an xvg time column onto a known total,
    because different gmx tools write the column in different units and nothing
    in the file says which. That is sound -- but its total defaults to 100.0, and
    a default is a pin. Pointed at the 10 ns sweep trajectories the controls were
    run under, it stretched a 10 ns axis onto a 100 ns one and every time
    coordinate in the report came out 10x too large, while the residence fraction
    and the RMSD statistics -- which are per-frame and carry no time -- stayed
    correct. Populated, plausible, and wrong: the shape in
    `how_this_project_breaks.md`.

    `nsteps * dt` is the run's own statement of its length, so ask the run.
    """
    for name in ("prod.mdp", "mdout.mdp"):
        p = rep / name
        if not p.is_file():
            continue
        vals = {}
        for line in p.read_text().splitlines():
            if "=" not in line or line.lstrip().startswith(";"):
                continue
            k, _, v = line.partition("=")
            k = k.strip().lower().replace("-", "_")
            if k in ("nsteps", "dt"):
                try:
                    vals[k] = float(v.split(";")[0].strip())
                except ValueError:
                    pass
        if "nsteps" in vals and "dt" in vals and vals["nsteps"] > 0:
            return vals["nsteps"] * vals["dt"] / 1000.0
    log.warning("no prod.mdp under %s — assuming %.0f ns", rep, default)
    return default


def series(rep: Path, er, total_ns: float = 100.0) -> dict:
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
        out[key] = (er.to_ns(a[:, 0], total_ns), a[:, 1])
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


def nac_series(cand: str, rep: Path, movie: Path,
               total_ns: float = 100.0) -> dict | None:
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
    # The movie frames span the production run, whatever its length -- see
    # prod_ns(). Hardcoding 100.0 put the controls' attack-geometry trace on a
    # 100 ns axis for a 10 ns trajectory.
    t = np.linspace(0.0, total_ns, len(dist))
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


def _key3_css() -> str:
    """The interaction key's styles, from the module that owns the key.

    Imported rather than copied, and imported LAZILY so a report can still be
    built if `shortlist_report`'s heavier dependencies are unavailable -- the
    page then loses the key's styling, not the run's numbers.
    """
    try:
        import shortlist_report as sr
        return sr.KEY3_CSS
    except Exception:                                      # noqa: BLE001
        return ""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--candidate", required=True)
    ap.add_argument("--rep", default="rep1")
    ap.add_argument("--out", default=None)
    ap.add_argument("--no-movie", action="store_true",
                    help="skip the 3D trajectory movie (it costs a trjconv pass)")
    # The 100 ns runs and the 10 ns sweep live under DIFFERENT roots (the split
    # md_residence_3ikd.py --work-root introduced). The controls were swept, not
    # elevated, so their trajectories sit under the sweep root -- without this
    # they could have no report at all, and the viewer had nothing to show.
    ap.add_argument("--work-root", default=None,
                    help=f"trajectory root holding <ident>/md/<rep> (default {MD})")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    er = _er()
    root = Path(args.work_root) if args.work_root else MD
    rep = root / args.candidate / "md" / args.rep
    if not rep.is_dir():
        raise SystemExit(f"no trajectory directory {rep}")

    total_ns = prod_ns(rep)
    log.info("%s: production length %.1f ns (from the mdp)", args.candidate, total_ns)
    s = series(rep, er, total_ns)
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
    if not isinstance(smiles, str):
        # Controls are not rows in D3/D4 -- they come from crystal structures -- so
        # the generated-candidate frames cannot describe them. The pose sidecar is
        # written next to the pose itself by whatever produced it, which makes it
        # the right place to ask: it is keyed on the pose, not on a generator.
        _sc = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith/"
                   f"pose_sidecars/{args.candidate}.json")
        if _sc.is_file():
            try:
                smiles = json.loads(_sc.read_text()).get("canonical_smiles")
            except Exception as exc:                       # noqa: BLE001
                log.warning("sidecar unreadable for %s: %s", args.candidate, exc)
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
    # ---- WHAT THE MOLECULE ACTUALLY TOUCHED (@tt8804, #63) -----------------
    # The interaction map, from the SAME movie the viewer above shows, so the
    # picture and the contact list cannot describe different frames.
    #
    # IMPORTED, NOT REIMPLEMENTED. `shortlist_report` already defines contact
    # occupancy, the representative frame, the 3D interaction view and the
    # colour key -- and the key is the part that must not fork: a colour meaning
    # "polar" on one page and something else on another is worse than no colour
    # at all. Importing across scripts/ is how nac_rank and nac_screen are
    # already used here.
    inter_block = ""
    if not args.no_movie and mpdb.is_file():
        try:
            import shortlist_report as sr
            rows_c, nfr = sr.contacts(mpdb)
            rx = None
            try:
                rx = sr.reactive_atom_index(
                    sr.ligand_mol(mpdb, str(s.get("smiles", ""))),
                    s.get("warhead_class"))
            except Exception:                          # noqa: BLE001
                pass                    # the view degrades to no anchor, not to nothing
            i3d = sr.interaction_3d(mpdb, rows_c, f"i3_{args.candidate}",
                                    rx_atom=rx)
            if i3d:
                occ = "".join(
                    f"<tr><td>{r[0]}{r[1]}</td><td>{r[2]*100:.0f}%</td>"
                    f"<td>{r[3]*100:.0f}%</td></tr>"
                    for r in rows_c[:12] if r[2] >= 0.20)
                inter_block = (
                    '<h2 class="ih">Interactions</h2>'
                    '<p class="note">Contacts over ' + str(nfr) +
                    ' frames of the same trajectory. Occupancy is the fraction '
                    'of frames in which the residue is in contact; polar is the '
                    'fraction with an N/O–N/O pair inside the cutoff. Every line '
                    'in the 3D view joins the actual closest atom pair in the '
                    'representative frame — nothing here is projected.</p>'
                    + i3d
                    + ('<table class="occ"><thead><tr><th>residue</th>'
                       '<th>occupancy</th><th>polar</th></tr></thead><tbody>'
                       + occ + "</tbody></table>" if occ else ""))
                log.info("interaction map: %d residues over %d frames",
                         len([r for r in rows_c if r[2] >= 0.20]), nfr)
        except Exception as exc:                       # noqa: BLE001
            # Same rule as the movie: a failed figure is RECORDED, never dropped,
            # so an absent map reads as a failure and not as "nothing touched".
            log.warning("interaction map failed: %s", exc)
            inter_block = rt.callout(
                "Interaction map unavailable",
                f"Contacts could not be computed: <code>{exc}</code>. The "
                "residence numbers above are unaffected.", "warn")

    verdict = ("Held" if not res["dissociated"] else
               f"Left at {res['left_at_ns']:.0f} ns")

    # The NAC series needs the fitted movie frames, so it comes after the movie.
    nacs = (nac_series(args.candidate, rep, mpdb, total_ns)
            if not args.no_movie else None)
    if nacs is None:
        log.warning("no warhead->SG distance/angle series for %s", args.candidate)
    img = figure(args.candidate, s, res, er, nacs)

    # THE SWEEP IS NOT A RESULT AND NO LONGER SITS BESIDE ONE (#55). It is
    # triage: it decided whether this molecule earned a 100 ns run. Presented in
    # the headline facts next to the engagement number it read as a second,
    # competing score, so it has left the selector rail entirely and appears here
    # in its own table, under a heading that says what it is for.
    facts = [("molecule", args.candidate), ("warhead class", cls or "unclassified")]
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
        # Directly after the movie, because it is the same trajectory read a
        # different way: the movie shows what moved, this shows what it touched.
        (f'<details class="panel"><summary>Interactions'
         f'<span class="hint">contact occupancy over the run, drawn on the '
         f'representative frame</span>'
         f'</summary><div class="pbody">{inter_block}</div></details>')
        if inter_block else "",
        '<details class="panel"><summary>How it was selected, and what that is worth'
        '<span class="hint">the 10 ns triage sweep — not a result</span></summary>'
        '<div class="pbody">',
        (('<table class="kv"><tbody>'
          f'<tr><th>attack-ready, 10 ns sweep</th><td>{sweep_ar*100:.1f}%</td></tr>'
          f'<tr><th>sustained visits</th><td>{sweep_v}</td></tr>'
          '</tbody></table>'
          '<p>These are the numbers that decided this molecule earned a 100 ns '
          'run. They are <strong>triage, not a measurement of it</strong> — the '
          'result is the 100&nbsp;ns engagement and RMSD above. D0075 and D0076 '
          'record what this sweep does and does not order correctly.</p>')
         if sweep_ar is not None else
         '<p>No 10 ns sweep reading for this molecule.</p>'),
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
        # The interaction key's styles come from the module that OWNS the key,
        # so the swatches cannot say one thing here and another on the shortlist.
        f"<style>{rt.CSS}{mov.VIEWER_CSS}{_key3_css()}</style></head><body>\n"
        + "\n".join(body) + "\n</body></html>")
    print(f"{args.candidate}: {verdict}, residence {res['residence_frac']:.3f} "
          f"-> {dest}")

    row = {"ident": args.candidate, "bpmd_occupancy": occ, "warhead_class": cls,
           "pose_rank": pose, **{k: v for k, v in res.items() if k != "status"}}
    pd.DataFrame([row]).to_csv(OUT.write(f"mdprio_{args.candidate}", ".csv"),
                               index=False)


if __name__ == "__main__":
    main()
