#!/usr/bin/env python3
"""
Purpose: build every artefact the presentation GUI needs, fresh, from the
         trajectories themselves -- no numbers copied from an older report.
Author: Timothy Wu (with Claude Code)
Date: 2026-09-10
Input: the 14-run manifest below, each entry naming a rep directory this
       script resolves and verifies rather than assumes
Output: presentation_data/<ident>.json (+ .rmsd.png, movie via md_movie)
        under append_only/inhibition/00_outputs/blacksmith/

@twu383, 2026-09-10: "I want to make a new gui from scratch for presentation
specifcally ... do not carry over any of the design and numbers from the last
guis." This script is the "no carried-over numbers" half of that: every
figure the page shows -- RMSD mean/max, warhead distance mean/max, the movie
itself -- is computed here, now, from the trajectory on disk. Nothing is read
from an existing .tier.json sidecar or an older report's cached PNG.

WHY MEAN/MAX WARHEAD DISTANCE IS COMPUTED HERE AND NOT READ FROM A TABLE.
`attack_sweep.geometry_stats` records median_dist_a and min_dist_a -- never
mean or max -- because the campaign's question was occupancy of a window, not
central tendency. The user asked specifically for mean and max here, which is
a different pair of statistics nothing on disk carries, so it is computed
fresh from the same per-frame distance array `nac_series` already produces
(the SAME resolver sweep_report.py uses for its own plots -- one reactive-atom
definition, not a second one invented for this page).
"""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import md_movie as mov                        # noqa: E402
from shared import outputs as sout                         # noqa: E402
import attack_sweep as asw                                 # noqa: E402
import mdprio_report as mp                                 # noqa: E402
import sweep_assets as sa                                  # noqa: E402
import shortlist_report as sr                               # noqa: E402

log = logging.getLogger("presentation-build")

MOD = Path("/data/lab_vm/modifiable/inhibition")

# ---------------------------------------------------------------------------
# THE MANIFEST. One entry per 100 ns run this presentation shows. `rep` is
# resolved once, here, and is either `sa.rep_dir(parent, pose_rank)` -- the
# one resolver the sweep path uses -- or a literal path for the two runs that
# went through direct elevation (`elevate_nac_v8`, not `attack_sweep_nac_v8`)
# and therefore have no pose_rank-keyed sibling directory to disambiguate.
# Resolved and printed by `python3 presentation_build.py --list` before any
# heavy work runs, so a bad path is caught before 10 minutes of GROMACS.
# ---------------------------------------------------------------------------
MANIFEST = [
    # -- family: thiadiazole-benzene. @twu383, 2026-09-10: "the lead result
    # is part of the tda benzne family" -- m50 groups here, not on its own. --
    dict(ident="t4_7b02d1dc4fd2_m50", label="t4_7b02d1dc4fd2 · m50",
         family="thiadiazole_benzene", parent="t4_7b02d1dc4fd2",
         rep=lambda: sa.rep_dir("t4_7b02d1dc4fd2", 51)),

    # t4_7b02d1dc4fd2_m50_rep2 is DELIBERATELY NOT HERE. @twu383, 2026-09-10:
    # "get rid of m50 rep 2 its not even the same starting pose just keep it
    # as its own mol." The pre-equilibration ligand_pose.sdf fed to both runs
    # was md5-identical (D0118); what differed is the pose each independent
    # equilibration handed to PRODUCTION (3.86 A vs 5.39 A at frame 1) -- one
    # molecule, one shown pose here, by the same editorial call applied to
    # m184's departure: this GUI shows the pose, not the replicate spread.

    dict(ident="t4_33b6dd2ff31f_m200", label="t4_33b6dd2ff31f · m200",
         family="thiadiazole_benzene", parent="t4_33b6dd2ff31f",
         rep=lambda: sa.rep_dir("t4_33b6dd2ff31f", 201)),
    dict(ident="t4_a7d09b6da355_m21", label="t4_a7d09b6da355 · m21",
         family="thiadiazole_benzene", parent="t4_a7d09b6da355",
         rep=lambda: sa.rep_dir("t4_a7d09b6da355", 22)),
    dict(ident="t4_75642786dc33_m195", label="t4_75642786dc33 · m195",
         family="thiadiazole_benzene", parent="t4_75642786dc33",
         rep=lambda: sa.rep_dir("t4_75642786dc33", 196)),
    dict(ident="t4_6ff34366c747_m255", label="t4_6ff34366c747 · m255",
         family="thiadiazole_benzene", parent="t4_6ff34366c747",
         rep=lambda: sa.rep_dir("t4_6ff34366c747", 256)),
    dict(ident="t4_e0b03662d460_m130", label="t4_e0b03662d460 · m130",
         family="thiadiazole_benzene", parent="t4_e0b03662d460",
         rep=lambda: sa.rep_dir("t4_e0b03662d460", 131)),

    # -- every other group is ITS OWN MOLECULE. `family` == `parent` for all
    # of these, which is what makes the header-derivation rule in
    # presentation_index.py correct with no per-molecule entry needed: the
    # only group that is not one molecule is thiadiazole_benzene above, and
    # that is the one case FAMILY_TITLE overrides. ----------------------------
    dict(ident="t4_7606519fce77_m218", label="m218",
         family="t4_7606519fce77", parent="t4_7606519fce77",
         rep=lambda: sa.rep_dir("t4_7606519fce77", 219)),
    dict(ident="t4_7606519fce77_m121", label="m121",
         family="t4_7606519fce77", parent="t4_7606519fce77",
         rep=lambda: sa.rep_dir("t4_7606519fce77", 122)),
    dict(ident="t4_7606519fce77_m186", label="m186",
         family="t4_7606519fce77", parent="t4_7606519fce77",
         rep=lambda: sa.rep_dir("t4_7606519fce77", 187)),
    dict(ident="t4_7606519fce77_m216", label="m216",
         family="t4_7606519fce77", parent="t4_7606519fce77",
         rep=lambda: sa.rep_dir("t4_7606519fce77", 217)),

    dict(ident="t4_b2092e63baa3_m67", label="m67",
         family="t4_b2092e63baa3", parent="t4_b2092e63baa3",
         rep=lambda: sa.rep_dir("t4_b2092e63baa3", 68)),

    dict(ident="t4_88ae42890e2d_m239", label="m239",
         family="t4_88ae42890e2d", parent="t4_88ae42890e2d",
         rep=lambda: MOD / "elevate_nac_v8/t4_88ae42890e2d_m239_100ns"
                          "/t4_88ae42890e2d_m239/md/rep1"),

    # t4_215b12bd9b34_m184 is DELIBERATELY NOT HERE. @twu383, 2026-09-10:
    # "if it left the site I dont want to see it." Mean ligand RMSD 18.5 A,
    # confirmed left (docs/state_of_the_project.md SS0c) -- not presentation
    # material. If it belongs anywhere it is a decision record, not this GUI.
]

FAMILY_ORDER = ["thiadiazole_benzene", "t4_7606519fce77", "t4_b2092e63baa3",
                "t4_88ae42890e2d"]
# THE ONLY OVERRIDE. @twu383, 2026-09-10: "the only subtitle on the viewer
# should be the mol names ... The only annotation should be the
# thiadiazole-benzene family." Every group not listed here shows its own
# `family` value as its header verbatim -- which, since `family == parent`
# for every group except this one, IS the molecule's name. No second title
# to keep in sync with the manifest by hand.
FAMILY_TITLE = {
    "thiadiazole_benzene": "Thiadiazole-benzene family",
}

OUT_TOPIC = sout.Topic("blacksmith", "presentation_data")


def build_one(entry: dict, *, force: bool) -> dict | None:
    ident = entry["ident"]
    rep = entry["rep"]()
    if rep is None or not (rep / "prod.tpr").is_file():
        log.error("%s: no resolvable rep dir (got %s)", ident, rep)
        return None

    movie = rep / "movie_pres.pdb"
    stamp = mov.movie_total_ps(movie) if movie.is_file() else None
    if force or stamp is None or abs(stamp - 100_000.0) > 500.0:
        log.info("%s: building movie from %s", ident, rep)
        got = mov.build_movie_pdb(rep, movie, n_frames=180)
        if got is None:
            log.error("%s: movie build failed", ident)
            return None
        stamp = mov.movie_total_ps(movie)
    if stamp is None or abs(stamp - 100_000.0) > 500.0:
        log.error("%s: movie does not cover 100 ns (got %s ps) -- refusing "
                  "to present it as one", ident, stamp)
        return None

    # Ligand RMSD, Angstrom -- attack_sweep's own converter, one place nm
    # becomes Angstrom for a trajectory (see its "UNITS ARE THE TRAP" note).
    rs = asw.rmsd_stats(rep)
    if not rs:
        log.error("%s: no rmsd.xvg readable at %s", ident, rep)
        return None

    # Warhead distance, mean/max, Angstrom -- fresh, over THIS movie's frames,
    # by the molecule's own recorded reactive atom (mp.reactive_atom / the
    # SAME resolver sweep_report.py uses, so this page and the sweep pages can
    # never disagree about which atom is "the warhead").
    sdf_ident = entry.get("sdf_ident", entry["parent"])
    # `total_ns` is NOT left at its keyword default. `nac_series` defaults to
    # 100.0, which happens to be right for every run here -- but "happens to
    # be right" is exactly the shape catalogue #32/#35 warns about: a default
    # nobody overrides, correct by coincidence rather than by construction.
    # `stamp` is the verified length (checked against the movie two lines
    # above), so it is passed explicitly.
    nac = mp.nac_series(sdf_ident, rep, movie, total_ns=stamp / 1000.0)
    if nac is None:
        log.error("%s: no reactive-atom series (nac_series returned None)",
                  ident)
        return None
    dist = np.asarray(nac["dist"], dtype=float)
    t_ns = np.asarray(nac["t"], dtype=float)

    # SMILES IS A MOLECULE PROPERTY, NOT A POSE PROPERTY -- always the PARENT
    # ident, even when `sdf_ident` (used above for the geometry lookup) names a
    # specific mode. `smiles_of()` resolves against `candidate_id`, which is
    # always the parent form; passing it a mode-suffixed ident (as this did for
    # m50's replicate 2 on the first pass) returns nothing, silently, and the
    # structure panel renders empty with no error anywhere in the chain.
    smi = sr.smiles_of(entry["parent"]) or ""
    svg = sr.depiction(smi) if smi else ""

    data = {
        "ident": ident, "label": entry["label"], "family": entry["family"],
        "parent": entry["parent"],
        "rmsd_mean_a": rs.get("rmsd_mean_a"), "rmsd_max_a": rs.get("rmsd_max_a"),
        "warhead_dist_mean_a": round(float(dist.mean()), 3),
        "warhead_dist_max_a": round(float(dist.max()), 3),
        "warhead_dist_min_a": round(float(dist.min()), 3),
        "n_frames": int(len(dist)), "total_ps": stamp,
        "mechanism": nac.get("mechanism"),
        "smiles": smi, "structure_svg": svg,
        "rep_dir": str(rep),
    }

    out_dir = OUT_TOPIC.dir
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / f"{ident}.json").write_text(json.dumps(data, indent=2))

    plot = plot_rmsd(rep, rs, out_dir / f"{ident}.rmsd.png")
    data["rmsd_plot"] = plot.name if plot else None
    wplot = plot_warhead(t_ns, dist, out_dir / f"{ident}.warhead.png")
    data["warhead_plot"] = wplot.name if wplot else None
    (out_dir / f"{ident}.json").write_text(json.dumps(data, indent=2))

    movie_dst = out_dir / f"{ident}.movie.pdb"
    movie_dst.write_bytes(movie.read_bytes())

    log.info("%s: OK  rmsd %.2f/%.2f A  warhead %.2f/%.2f A  (%d frames)",
             ident, data["rmsd_mean_a"], data["rmsd_max_a"],
             data["warhead_dist_mean_a"], data["warhead_dist_max_a"],
             data["n_frames"])
    return data


def plot_rmsd(rep: Path, rs: dict, dest: Path) -> Path | None:
    """A fresh RMSD-vs-time figure -- not one carried over from an older
    report's cache. Reads the same rmsd.xvg `rmsd_stats` already parsed."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    f = rep / "rmsd.xvg"
    if not f.is_file():
        return None
    xs, ys = [], []
    for ln in f.read_text(errors="replace").splitlines():
        ln = ln.strip()
        if not ln or ln[0] in "#@":
            continue
        parts = ln.split()
        if len(parts) >= 2:
            xs.append(float(parts[0]))
            ys.append(float(parts[1]) * 10.0)          # nm -> Angstrom
    if not xs:
        return None
    xs = np.asarray(xs)
    # `.xvg` time axis is in ns already (gmx rms -tu ns, per gromacs_analysis).
    fig, ax = plt.subplots(figsize=(5.2, 2.2), dpi=140)
    ax.plot(xs, ys, lw=0.7, color="#2563eb")
    ax.axhline(10.0, color="#94a3b8", lw=0.8, ls="--")   # 1.0 nm departure bar
    ax.set_xlabel("time (ns)", fontsize=9)
    ax.set_ylabel("ligand RMSD (Å)", fontsize=9)
    ax.tick_params(labelsize=8)
    fig.tight_layout()
    fig.savefig(dest)
    plt.close(fig)
    return dest


def plot_warhead(t_ns: np.ndarray, dist_a: np.ndarray, dest: Path) -> Path | None:
    """Warhead-to-Cys113 distance vs time -- the geometry readout, plotted
    beside the ligand RMSD rather than left as a table of summary statistics.
    Shares its time axis and reactive-atom resolver with `nac_series`, so this
    plot and the 3D viewer's own distance readout cannot disagree.
    """
    if not len(dist_a):
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from shared import nac_criterion as nac

    lo, hi = nac.attack_ready_window()          # the 2.8-3.5 A gate, on config
    fig, ax = plt.subplots(figsize=(5.2, 2.2), dpi=140)
    ax.plot(t_ns, dist_a, lw=0.7, color="#c2410c")
    ax.axhspan(lo, hi, color="#16a34a", alpha=0.12, lw=0)
    ax.axhline(lo, color="#16a34a", lw=0.7, ls=":")
    ax.axhline(hi, color="#16a34a", lw=0.7, ls=":")
    ax.set_xlabel("time (ns)", fontsize=9)
    ax.set_ylabel("warhead–Cys113 (Å)", fontsize=9)
    ax.tick_params(labelsize=8)
    fig.tight_layout()
    fig.savefig(dest)
    plt.close(fig)
    return dest


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--list", action="store_true",
                    help="resolve and print every rep dir, build nothing")
    ap.add_argument("--force", action="store_true",
                    help="rebuild movies even if already stamped at 100 ns")
    ap.add_argument("--only", nargs="*", default=None,
                    help="restrict to these idents")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    entries = MANIFEST if not args.only else [
        e for e in MANIFEST if e["ident"] in args.only]

    if args.list:
        bad = 0
        for e in entries:
            rep = e["rep"]()
            ok = rep is not None and (rep / "prod.tpr").is_file()
            if not ok:
                bad += 1
            print(f"  {'OK ' if ok else 'BAD'}  {e['ident']:<32s} {rep}")
        print(f"\n  {len(entries)} entries, {bad} unresolved")
        return

    results = []
    for e in entries:
        try:
            d = build_one(e, force=args.force)
        except Exception as exc:                          # noqa: BLE001
            log.error("%s: FAILED: %s", e["ident"], exc)
            d = None
        results.append((e["ident"], d is not None))

    ok = sum(1 for _, s in results if s)
    print(f"\n  {ok}/{len(results)} built")
    for ident, s in results:
        if not s:
            print(f"    FAILED: {ident}")


if __name__ == "__main__":
    main()
