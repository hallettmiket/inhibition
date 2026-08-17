#!/usr/bin/env python3
"""
Purpose: per-mode 10 ns assets — movie and trajectory plots — for the Sweep page.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-16
Input: --worklist <sweep_gaps_N.csv> [--limit N]
Output: mdprio_reports/sweep_assets/<ident>.{pdb,png}

@tt8804: "sweep results should show rmsd and 10 ns md movie, it should look just
like MD results but with the sweep info also in a table on the viewer side".

NOTHING IS RE-SIMULATED, AND ALMOST NOTHING IS RE-COMPUTED. Every finished sweep
already left `rmsd.xvg`, `mindist.xvg` and `numcont.xvg` beside its trajectory --
1,001 rows each, written by the run itself. The plots read those. The movie is
`md_movie.build_movie_pdb`, the SAME function the 100 ns page uses, pointed at
the sweep's rep directory: 3 seconds per mode, so the whole campaign is minutes.

THE DISTANCE PANEL CARRIES THE CRITERION'S OWN WINDOW. A distance trace without
the 2.8-4.2 A band is a line that a reader has to hold a threshold against in
their head; with it, "did this molecule ever get into position" is answered by
looking. The band is read from `nac_criterion`, not typed in here, so the picture
and the score cannot disagree about what attack range means.

ASSETS ARE PER MODE, NOT PER MOLECULE. A molecule can have several modes swept
and they are different trajectories -- naming a file by the parent would let one
mode's movie stand for another's, which is the mistake `mode_key` exists to stop.
"""

from __future__ import annotations

import argparse
import base64
import glob
import io
import logging
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

log = logging.getLogger("sweep-assets")
OUT = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith/"
           "mdprio_reports/sweep_assets")
SWEEP_ROOT = Path("/data/lab_vm/modifiable/inhibition/attack_sweep_10ns")


def rep_dir(parent: str, pose_rank: int | None = None) -> Path | None:
    """The rep directory for THIS MODE's finished 10 ns run.

    KEYED ON POSE RANK, NOT JUST THE MOLECULE. The runner writes
    `rank<pose_rank>_<ps>ps/<parent>/md/rep1`, so a molecule with four swept
    modes has four sibling directories. Matching on the parent alone returned
    whichever sorted first and handed EVERY mode of that molecule the same
    trajectory -- one mode's RMSD, one mode's plots, one mode's movie, shown
    under four different idents.

    It was visible in the ranking as identical `rmsd_max` to six decimal places
    for modes that are different binding poses, which is not something two
    trajectories do. Caught before 12 GPU-days were committed on that order.

    `pose_rank` is required for a correct answer; without it this still resolves
    to the first finished run, which is right only for single-mode molecules.
    A finished run is one whose prod.log says so -- the same completeness test
    the sweep's resume guard uses.
    """
    pats = ([f"rank{int(pose_rank)}_*ps/{parent}/md/rep1"] if pose_rank is not None
            else [f"*/{parent}/md/rep1"])
    for pat in pats:
        for p in sorted(SWEEP_ROOT.glob(pat)):
            log_f = p / "prod.log"
            if log_f.is_file() and "Finished mdrun" in log_f.read_text(errors="replace"):
                return p
    return None


def _xvg(path: Path):
    """(x, y) from a GROMACS .xvg, with x converted to NANOSECONDS.

    THE UNIT IS READ FROM THE FILE, NOT ASSUMED. These traces declare
    `@ xaxis label "Time (ns)"` and run 0 -> 10.0; the first version of this
    divided by 1000 on the assumption they were picoseconds, and every plot came
    out with a 10 ns run spanning 0.010 ns. The axis was wrong by three orders of
    magnitude and the shape of the trace looked entirely normal, which is why it
    survived a review -- @tt8804 caught it by reading the tick labels.
    """
    if not path.is_file():
        return None, None
    scale = 1.0
    xs, ys = [], []
    for ln in path.read_text(errors="replace").splitlines():
        if ln.startswith("@") and "xaxis" in ln and "label" in ln:
            lab = ln.lower()
            if "(ps)" in lab:
                scale = 1e-3
            elif "(fs)" in lab:
                scale = 1e-6
        if not ln or ln[0] in "#@":
            continue
        f = ln.split()
        if len(f) >= 2:
            try:
                xs.append(float(f[0]) * scale); ys.append(float(f[1]))
            except ValueError:
                continue
    return np.array(xs), np.array(ys)


def plots(rep: Path, ident: str) -> str:
    """RMSD and warhead-to-sulfur distance over the 10 ns, as base64 PNG."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from shared import nac_criterion as nac

    t_r, rmsd = _xvg(rep / "rmsd.xvg")
    t_d, dist = _xvg(rep / "mindist.xvg")
    if t_r is None and t_d is None:
        return ""
    fig, ax = plt.subplots(2, 1, figsize=(7.2, 3.9), dpi=150, sharex=True)
    if t_r is not None:
        ax[0].plot(t_r, rmsd, lw=1.0, color="#0072ce")
        ax[0].set_ylabel("ligand RMSD (nm)")
        # MAX RMSD MARKED, because it is the number the 100 ns page RANKS on --
        # how far the molecule ever got from where it started. A trace without it
        # makes the reader eyeball the peak, and the eye reads a noisy maximum
        # low. Drawn as a line plus its value so the plot states it outright.
        mx = float(np.nanmax(rmsd))
        ax[0].axhline(mx, color="#8a6d1f", ls="--", lw=1.0)
        ax[0].text(0.995, mx, f" max {mx:.3f} nm ", ha="right", va="bottom",
                   fontsize=7.5, color="#8a6d1f",
                   transform=ax[0].get_yaxis_transform())
        ax[0].set_ylim(0, max(mx * 1.18, 0.05))
    if t_d is not None:
        # THE WINDOW THE SCORE USES, drawn from the criterion itself.
        ax[1].axhspan(nac.NAC_DIST_MIN, nac.NAC_DIST_MAX, color="#0f7a54",
                      alpha=0.13, lw=0)
        ax[1].plot(t_d, dist * 10.0, lw=1.0, color="#b3261e")
        ax[1].set_ylabel("warhead–S$\\gamma$ (Å)")
        ax[1].text(0.995, nac.NAC_DIST_MAX, " attack range", ha="right",
                   va="bottom", fontsize=7, color="#0f7a54",
                   transform=ax[1].get_yaxis_transform())
    ax[1].set_xlabel("time (ns)")
    for a in ax:
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
    fig.suptitle(ident, fontsize=8, color="#5b6b80", y=0.995)
    fig.tight_layout()
    buf = io.BytesIO(); fig.savefig(buf, format="png"); plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--worklist", required=True)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--plots-only", action="store_true",
                    help="rebuild the figures and leave the movies alone -- "
                         "plots are ~1 s, movies ~3 s and rarely change")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    import pandas as pd
    from shared import md_movie as mov

    wl = pd.read_csv(args.worklist)
    OUT.mkdir(parents=True, exist_ok=True)
    idents = list(wl.ident.astype(str))[:args.limit]
    # ident -> pose_rank, so each mode's OWN trajectory is found.
    prank = dict(zip(wl.ident.astype(str), wl.pose_rank.astype(int))) \
        if "pose_rank" in wl.columns else {}
    n_mov = n_png = 0
    miss: dict = {}
    for i, ident in enumerate(idents, 1):
        parent = ident.rsplit("_m", 1)[0]
        pdb, png = OUT / f"{ident}.pdb", OUT / f"{ident}.png"
        want_mov = not args.plots_only
        # `--plots-only` REGENERATES the figures. Skipping when the PNG exists
        # made the flag a no-op in exactly the case it is for -- rebuilding the
        # plots after changing how they are drawn -- and the run printed
        # "0 plots" while looking like it had succeeded.
        if not (args.force or args.plots_only) and png.is_file() and pdb.is_file():
            continue
        rep = rep_dir(parent, prank.get(ident))
        if rep is None:
            miss[parent] = "no finished 10 ns run"
            continue
        try:
            if want_mov and (args.force or not pdb.is_file()):
                mov.build_movie_pdb(rep, pdb, total_ps=10_000.0)
                n_mov += 1
            if args.force or args.plots_only or not png.is_file():
                b64 = plots(rep, ident)
                if b64:
                    png.write_bytes(base64.b64decode(b64)); n_png += 1
        except Exception as exc:                           # noqa: BLE001
            # Recorded, never dropped: a molecule with no movie must read as a
            # failure to build one, not as a run that did not happen.
            miss[ident] = str(exc)[:90]
        if i % 25 == 0:
            log.info("  %d/%d  movies +%d  plots +%d", i, len(idents), n_mov, n_png)
    print(f"\n  {n_mov} movies, {n_png} plots -> {OUT}")
    if miss:
        print(f"  {len(miss)} without assets:")
        for k, v in list(miss.items())[:6]:
            print(f"    {k}: {v}")


if __name__ == "__main__":
    main()
