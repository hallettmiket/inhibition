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
import os
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
from shared import mode_key                         # noqa: E402
from shared import run_paths as rp        # noqa: E402
from shared import residence_tier as rtier            # noqa: E402
from shared import target_config as _tc               # noqa: E402

#: The triage sweep's length, DERIVED. Every one of these strings said "10 ns"
#: while the sweep has run at 8 ns since D0085 (@tt8804: "update the gui to say
#: 8 ns sweep not 10"). A length written into prose is a number that stops
#: tracking the run the first time the spec changes.
_SWEEP_NS = int(round(_tc.md_sweep_ps() / 1000))

log = logging.getLogger("mdprio-report")
MD = rp.residence_work()
OUT = sout.Topic("blacksmith", rp.reports_topic())
# The residence criterion: ligand RMSD <= 1.0 nm. Imported rather than redefined
# so the tier rule and this page agree by construction.
BOUND_NM = rtier.BOUND_NM

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
    # `prot_rmsd` IS COMPUTED ON DEMAND, not read passively: it did not exist
    # when these runs were analysed, and it is one `gmx rms` over the corrected
    # trajectory that is already on disk. Failure is silent by design -- an old
    # run without `whole.xtc` keeps the single trace it has.
    try:
        from shared import gromacs_analysis as ga
        ga.protein_rmsd(rep)
    except Exception as exc:                                   # noqa: BLE001
        log.debug("%s: no protein RMSD (%s)", rep.name, exc)
    for key, fname in (("rmsd", "rmsd.xvg"), ("prot_rmsd", "rmsd_protein.xvg"),
                       ("mindist", "mindist.xvg"),
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


def reactive_atom(cand: str, rep: Path) -> dict | None:
    """The candidate's reactive atom, BY THE CLASS THE MOLECULE IS RECORDED AS.

    Returns {"heavy_idx", "match", "mechanism", "class_id", "n_heavy"} where
    `heavy_idx` indexes the ligand's HEAVY atoms in sdf order -- which is the
    order the MOL residue appears in a movie frame.

    TWO SELECTION BUGS LIVE HERE AND BOTH ARE FIXED BY ASKING THE MOLECULE.

    1. This used to walk `warhead_classes_10.csv` and take the FIRST class whose
       SMARTS matched. Measured 2026-09-02 over 98 finished nac_v8 sweeps: 62
       acrylamides were reported as `naphthoquinone_c2` (row 5 beats row 7) and
       34 `bdhi_c4` as `bdhi_c5` (row 3 beats row 4). It happened to be harmless
       -- the same reactive ATOM in 98 of 98, and the paired classes share a
       mechanism, so every distance and angle was identical -- but it is
       selection by file order, and the next class added to that CSV decides
       silently whether it stays harmless.
    2. `elevation_report.surface_payload` selected the reactive atom as the one
       NAMED `C10`, which is sulfopin's name for it in 6VAJ and nothing else's.
       `tests/test_crystal_pose_audit.py` already says so in as many words:
       "6VAJ calls it C10; the other five covalent Pin1 entries call the
       equivalent atom C19, C14, C24, C12 and C3."

    Falls back to a first-match scan ONLY when the candidate's own class cannot
    be resolved, and says so in the log rather than doing it silently.
    """
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")

    sdf = rep.parent.parent / "ligand_pose.sdf"
    if not sdf.is_file():
        return None
    mol = Chem.SDMolSupplier(str(sdf), removeHs=False)[0]
    if mol is None:
        return None
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    pos = {old_i: new_i for new_i, old_i in enumerate(heavy)}

    wh = pd.read_csv(REPO / "data/reference/warhead_classes_10.csv")

    # THE MOLECULE'S OWN CLASS FIRST, from the candidate frame that defines it.
    recorded = None
    try:
        frames = sorted((rp.DATA / "04_t4_combinatorial").glob("D4_*.parquet"),
                        key=lambda f: int(f.stem.split("_")[1]))
        if frames:
            df = pd.read_parquet(frames[-1], columns=["candidate_id",
                                                      "warhead_class"])
            hit = df[df.candidate_id == cand]
            if len(hit):
                recorded = str(hit.iloc[0].warhead_class)
    except Exception as exc:                              # noqa: BLE001
        log.debug("%s: could not read the recorded warhead class (%s)", cand, exc)

    order = []
    if recorded is not None:
        order = [r for r in wh.itertuples() if r.class_id == recorded]
        if not order:
            log.warning("%s: recorded warhead class %r is not in the library",
                        cand, recorded)
    if not order:
        log.warning("%s: falling back to a first-SMARTS-match scan; the class "
                    "this reports is the first one that matched, not the "
                    "molecule's own", cand)
        order = list(wh.itertuples())

    for r in order:
        patt = Chem.MolFromSmarts(r.reactive_atom_smarts)
        if patt is None:
            continue
        ms = mol.GetSubstructMatches(patt)
        if not ms:
            continue
        match = ms[0]
        if any(i not in pos for i in match):
            continue
        return {"heavy_idx": pos[match[0]],
                "match": tuple(pos[i] for i in match),
                "mechanism": r.mechanism, "class_id": r.class_id,
                "n_heavy": len(heavy)}
    if recorded is not None:
        log.warning("%s: recorded as %r but its SMARTS does not match the MD "
                    "ligand", cand, recorded)
    return None


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
    import shared.nac_criterion as nac

    if not movie.is_file():
        return None
    # ONE RESOLVER, shared with the movie payload. Two implementations of "the
    # reactive atom" is how a plot and the structure beneath it came to disagree
    # by a median of 3.11 A while both looked right.
    ra = reactive_atom(cand, rep)
    if ra is None:
        log.warning("%s: no reactive SMARTS match on the MD ligand", cand)
        return None
    match, mech = ra["match"], ra["mechanism"]

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
            # BY RESIDUE NUMBER. 3IKD has Cys57 as well as Cys113, so matching
            # on the residue NAME alone took whichever CYS SG appeared last in
            # the frame. It is Cys113 today only because 63 sorts after 7 --
            # correct by file order, which is exactly the guarantee catalogue
            # #38 says not to rely on. `surface_payload` selects by number and
            # this now agrees with it.
            elif (resn == "CYS" and name == "SG"
                  and int(ln[22:26]) == 113 - mov.PIN1_OFFSET):
                sg = np.array(xyz)
    if not frames_lig:
        return None

    idx = list(match)                                    # already heavy indices
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


def _promote_bar() -> float | None:
    """The max-ligand-RMSD gate, from config -- one definition, not a constant
    retyped here that would drift the moment the bar is retuned."""
    try:
        from shared import target_config as tc
        return float(tc.get("md.sweep_survivor_rmsd_nm", default=0.35))
    except Exception:                                          # noqa: BLE001
        return None


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
    titles = {"rmsd": "RMSD after superposing on protein CA (nm) — "
                      "ligand vs. the protein itself",
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
            # THE PROTEIN'S OWN DISPLACEMENT, on the same axis. This panel is
            # measured AFTER superposing on protein CA, so it rises both when
            # the ligand slides out of a rigid pocket and when the protein
            # relaxes around a ligand that never let go. @tt8804: "if they
            # change tgt then they are fine". Drawn behind the ligand and in a
            # neutral grey: it is the reference, not the reading.
            if "prot_rmsd" in s:
                tp, yp = s["prot_rmsd"]
                # The ligand trace is the one already drawn above; relabelled
                # rather than redrawn, so there is exactly one line per series.
                ax.lines[-1].set_label("ligand")
                ax.lines[-1].set_zorder(2)
                ax.plot(tp, yp, lw=0.7, color="#8a94a6", alpha=0.85, zorder=1,
                        label="protein CA")
                ax.legend(loc="lower left", fontsize=7.5, frameon=False,
                          ncol=2, handlelength=1.4, borderaxespad=0.2)
            # THE LINE THAT ACTUALLY DECIDES, drawn alongside the prereg one.
            # `BOUND_NM` (1.0) is the residence criterion -- did it dissociate
            # -- and essentially everything passes it. Promotion to 100 ns, and
            # onward to BPMD, is gated at `md.sweep_survivor_rmsd_nm` (0.35), so
            # a panel showing only the 1.0 line shows only the test nothing
            # fails. @tt8804: "md results should be updated to our new specs".
            _bar = _promote_bar()
            if _bar:
                ax.axhline(_bar, ls="--", lw=1, color="#8a6d1f")
                ax.text(0.995, _bar, f" promotes ≤ {_bar:.2f} nm ", va="bottom",
                        ha="right", transform=ax.get_yaxis_transform(),
                        fontsize=8, color="#8a6d1f")
            # AND THE AXIS IS SCALED TO THE DECISION, not to the prereg line.
            # Forcing the top to 1.0 nm squashed a held complex into the bottom
            # quarter of the panel, where a ligand at 0.15 and one at 0.30 --
            # promote and do-not-promote -- are the same two pixels. A run that
            # really does leave auto-scales past 1.0 and shows that line anyway.
            _top = max(float(np.nanmax(y)),
                       *( [float(np.nanmax(s["prot_rmsd"][1]))]
                          if "prot_rmsd" in s else [] ),
                       (_bar or 0.0) * 1.15)
            _fits = _top >= BOUND_NM
            if not _fits:
                ax.set_ylim(0, _top * 1.15)
            # The dissociation line is drawn only when it is ON the axis. Drawn
            # at y=1.0 under a 0.4 nm axis, `get_yaxis_transform` places its
            # LABEL outside the axes -- text floating in the figure margin
            # pointing at an invisible line. When it is off-scale the fact that
            # matters is that it was never approached, so the title says so.
            if _fits:
                ax.axhline(BOUND_NM, ls="--", lw=1, color=rt.SERIES["alert"])
                ax.text(0.995, BOUND_NM, f" bound ≤ {BOUND_NM} nm", va="bottom",
                        ha="right", transform=ax.get_yaxis_transform(),
                        fontsize=8, color=rt.SERIES["alert"])
            else:
                ax.set_title(f"{titles[k]} — never approached the "
                             f"{BOUND_NM} nm dissociation line",
                             loc="left", fontsize=9)
            if res.get("left_at_ns"):
                ax.axvline(res["left_at_ns"], lw=1, color=rt.SERIES["alert"])
        if k == "nac_dist":
            # The attack window, so "close enough to react" is on the plot
            # rather than in the reader's head.
            # The GATE's band, not the screen's window -- see
            # nac_criterion.attack_ready_window.
            ax.axhspan(*nacm.attack_ready_window(), color=rt.SERIES["ref"],
                       alpha=0.15)
            if nacm.NAC_DIST_MAX > nacm.attack_ready_window()[1]:
                ax.axhline(nacm.NAC_DIST_MAX, color=rt.SERIES["ref"], ls=":",
                           lw=0.8, alpha=0.55)
            ax.text(0.995, nacm.attack_ready_window()[1],
                    f" attack ready ≤{nacm.attack_ready_window()[1]:.1f} Å",
                    va="bottom", ha="right",
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
    # THE TRAJECTORY IS PER RUN; EVERYTHING ELSE IS PER MOLECULE. Once a run is
    # identified `<parent>_m<mode>`, the workdir and this report's filename follow
    # the RUN, while the SMILES, the warhead class, the pose sidecar and the sweep
    # row are all facts about the MOLECULE. Looking those up under the run ident
    # would miss every one of them and render a report about a molecule the
    # project appears to know nothing about. Legacy idents have no suffix, so
    # `split_ident` returns them unchanged.
    parent, run_mode = mode_key.split_ident(args.candidate)
    rep = root / args.candidate / "md" / args.rep
    if not rep.is_dir():
        raise SystemExit(f"no trajectory directory {rep}")

    total_ns = prod_ns(rep)
    log.info("%s: production length %.1f ns (from the mdp)", args.candidate, total_ns)
    s = series(rep, er, total_ns)
    res = residence(s)
    if res["status"] != "ok":
        raise SystemExit(f"{args.candidate}: {res['status']}")
    occ, cls, pose = PREDICTED.get(parent, (None, None, None))

    # PREDICTED only covers the six BPMD-priority molecules. Everything swept
    # since then fell through to "—" for class and pose, so the masthead of every
    # newer report advertised missing data it could have looked up. Class comes
    # from the ranking, the sweep readings from the sweep, and the SMILES gives
    # the depiction -- all keyed on the molecule, all already on disk.
    import glob as _g
    smiles = None
    for _sub, _stem in (("04_t4_combinatorial", "D4"), ("03_t3_reinvent", "D3")):
        _fs = [str(x) for x in rp.frames(_stem)]
        if not _fs:
            continue
        _fr = pd.read_parquet(_fs[-1]).drop_duplicates("candidate_id").set_index("candidate_id")
        if parent in _fr.index:
            smiles = _fr.loc[parent].get("canonical_smiles")
            if not cls:
                cls = _fr.loc[parent].get("warhead_class")
            break
    if not isinstance(smiles, str):
        # Controls are not rows in D3/D4 -- they come from crystal structures -- so
        # the generated-candidate frames cannot describe them. The pose sidecar is
        # written next to the pose itself by whatever produced it, which makes it
        # the right place to ask: it is keyed on the pose, not on a generator.
        _sc = rp.sidecars() / f"{parent}.json"
        if _sc.is_file():
            try:
                smiles = json.loads(_sc.read_text()).get("canonical_smiles")
            except Exception as exc:                       # noqa: BLE001
                log.warning("sidecar unreadable for %s: %s", args.candidate, exc)
    if not cls:
        for _t, _sc in (("T4", "conditional_eb"), ("T3", "enrichment_conditional")):
            # THIS RUN'S RANKING. This glob carried no topic at all, so it
            # matched every screen's tables and `[-1]` picked by string order --
            # a warhead class read off whichever run sorted last.
            _fs = sorted(_g.glob(str(rp.BLACKSMITH / "rank_v2" /
                                     f"rank_v2_{_t}_{rp.topic()}_{_sc}_*.csv")))
            if not _fs:
                continue
            _rk = pd.read_csv(_fs[-1]).drop_duplicates("parent_ident").set_index("parent_ident")
            if parent in _rk.index:
                cls = _rk.loc[parent, "warhead_class"]; break
    sweep_ar = sweep_v = None
    # THIS RUN'S SWEEP. The unscoped `attack_sweep/` here put a previous
    # screen's triage readings on this run's report page.
    _fs = sorted(_g.glob(str(rp.sweep_dir() / "attack_sweep_*.csv")),
                 key=os.path.getmtime)
    if _fs:
        _sw = pd.concat([pd.read_csv(f) for f in _fs], ignore_index=True)
        _sw = _sw[(_sw.get("sweep_ps", 0) > 1000) & (_sw.status == "ok")
                  & (_sw.parent_ident == parent)]
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
                _ra = reactive_atom(cand, rep)
                if _ra is None:
                    raise ValueError(
                        f"{cand}: cannot resolve the reactive atom for the "
                        f"movie; refusing to label an arbitrary atom")
                pdb_txt, dsg, labels, lpos = er.surface_payload(
                    mpdb, reactive_idx=_ra["heavy_idx"])
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

    # THREE TIERS, NOT A PASS/FAIL (@tt8804: "held in pocket, held but not
    # optimal and below max .35 is optimal"). Against 1.0 nm alone this page
    # could not tell a ligand pinned at the warhead from one rattling around the
    # site; against 0.35 alone it called a run that never left a failure.
    # shared/residence_tier owns the rule, so this page and the combined rail
    # cannot drift on what "optimal" means.
    tier_key = rtier.tier(res.get("rmsd_max_nm"), res.get("dissociated"),
                          res.get("residence_frac"))
    verdict = (f"Left at {res['left_at_ns']:.0f} ns" if tier_key == "left"
               else rtier.label(tier_key).capitalize())

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
    facts = [("molecule", parent),
             ("binding mode",
              f"m{run_mode}" if run_mode is not None else "not recorded"), ("warhead class", cls or "unclassified")]
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
            # THE STAGE AND THE RELEASE, BOTH DERIVED. This was the literal
            # "MD-PRIORITY · 2.2.0 BORNITE" -- wrong twice over, since 2.2.0 is
            # Chalcopyrite and Bornite is 2.1.0, so every 100 ns report named a
            # release that never existed, two releases back. "MD-PRIORITY" was
            # also the 2.2.0 name for a prereg experiment, not the stage this
            # page now reports.
            rt.eyebrow(f"PRODUCTION MD · {int(total_ns)} NS"), facts),
        (f'<div class="structrow"><div class="structbox">{struct_svg}</div>'
         f'<div class="structnote"><b>{cls or "unclassified"}</b>'
         + (f' &middot; attack-ready {sweep_ar*100:.1f}% of the {_SWEEP_NS} ns sweep'
            f' over {sweep_v} sustained visits' if sweep_ar is not None else '')
         + '</div></div>') if struct_svg else "",
        f'<style>{rtier.TIER_CSS}</style>'
        f'<p>{rtier.badge(tier_key)} '
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
        f'<span class="hint">the {_SWEEP_NS} ns triage sweep — not a result</span></summary>'
        '<div class="pbody">',
        (('<table class="kv"><tbody>'
          f'<tr><th>attack-ready, {_SWEEP_NS} ns sweep</th><td>{sweep_ar*100:.1f}%</td></tr>'
          f'<tr><th>sustained visits</th><td>{sweep_v}</td></tr>'
          '</tbody></table>'
          '<p>These are the numbers that decided this molecule earned a 100 ns '
          'run. They are <strong>triage, not a measurement of it</strong> — the '
          'result is the 100&nbsp;ns engagement and RMSD above. D0075 and D0076 '
          'record what this sweep does and does not order correctly.</p>')
         if sweep_ar is not None else
         f'<p>No {_SWEEP_NS} ns sweep reading for this molecule.</p>'),
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
    # THE TIER TRAVELS WITH THE PAGE. The md100 row carries rmsd max but not
    # whether the ligand ever came back, and max alone cannot separate "crossed
    # 1.0 nm and returned" from "left" -- t4_071099f4034c_m1 is 4 frames out of
    # 10,001 above the line and is held. Rather than have the rail re-parse
    # rmsd.xvg and re-derive the rule, this run writes what it already computed
    # and the rail reads it. A run with no sidecar is shown as not scored, not
    # as left.
    (dest.with_suffix(".tier.json")).write_text(json.dumps({
        "ident": args.candidate, "tier": tier_key,
        "label": rtier.label(tier_key), "colour": rtier.colour(tier_key),
        "rmsd_mean_nm": res.get("rmsd_mean_nm"),
        "rmsd_max_nm": res.get("rmsd_max_nm"),
        "residence_frac": res.get("residence_frac"),
        "dissociated": res.get("dissociated"),
        "left_at_ns": res.get("left_at_ns"),
        "bound_nm": rtier.BOUND_NM, "optimal_nm": rtier.optimal_nm(),
    }, indent=2))
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
