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

from shared import run_paths as rp                  # noqa: E402

log = logging.getLogger("sweep-assets")
OUT = rp.reports_dir() / "sweep_assets"
SWEEP_ROOT = rp.sweep_work()


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
    # THE PROTEIN'S OWN DISPLACEMENT, on the same axis as the ligand's. Ligand
    # RMSD is measured after superposing on protein CA, so it rises both when
    # the ligand leaves a rigid pocket and when the protein relaxes around a
    # still-bound ligand. @tt8804: "if they change tgt then they are fine".
    # Computed on demand from the persisted corrected trajectory; a mode without
    # one simply plots the single trace it has.
    t_p = p_rmsd = None
    try:
        from shared import gromacs_analysis as ga
        xvg = ga.protein_rmsd(rep)
        if xvg is not None:
            t_p, p_rmsd = _xvg(xvg)
    except Exception as exc:                                   # noqa: BLE001
        log.debug("%s: no protein RMSD (%s)", ident, exc)
    if t_r is None and t_d is None:
        return ""
    fig, ax = plt.subplots(2, 1, figsize=(7.2, 3.9), dpi=150, sharex=True)
    if t_p is not None and len(t_p):
        ax[0].plot(t_p, p_rmsd, lw=0.9, color="#8a94a6", alpha=0.85,
                   label="protein CA", zorder=1)
    if t_r is not None:
        ax[0].plot(t_r, rmsd, lw=1.0, color="#0072ce", label="ligand", zorder=2)
        ax[0].set_ylabel("RMSD (nm)")
        # MAX RMSD MARKED, because it is the number the 100 ns page RANKS on --
        # how far the molecule ever got from where it started. A trace without it
        # makes the reader eyeball the peak, and the eye reads a noisy maximum
        # low. Drawn as a line plus its value so the plot states it outright.
        mx = float(np.nanmax(rmsd))
        ax[0].axhline(mx, color="#8a6d1f", ls="--", lw=1.0)
        ax[0].text(0.995, mx, f" ligand max {mx:.3f} nm ", ha="right",
                   va="bottom", fontsize=7.5, color="#8a6d1f",
                   transform=ax[0].get_yaxis_transform())
        # The top of the axis must cover BOTH traces, or a protein that moved
        # further than the ligand is silently clipped -- which reads as a
        # protein that stayed still, the opposite of what happened.
        top = mx if t_p is None or not len(t_p) else max(mx, float(np.nanmax(p_rmsd)))
        ax[0].set_ylim(0, max(top * 1.18, 0.05))
        if t_p is not None and len(t_p):
            ax[0].legend(loc="upper left", fontsize=7, frameon=False,
                         ncol=2, handlelength=1.4, borderaxespad=0.2)
    if t_d is not None:
        # THE BAND THE GATE USES, drawn from the gate itself. This shaded
        # 2.8-4.2 A, the screen's NAC window, while the sweep is judged at
        # 2.8-3.5 (D0111) -- so a trace could sit inside the green zone for most
        # of the run and still score 0% engaged on the same page.
        _lo, _hi = nac.attack_ready_window()
        ax[1].axhspan(_lo, _hi, color="#0f7a54", alpha=0.15, lw=0)
        # The wider NAC window is kept as a hairline rather than deleted: it is
        # still the criterion the SCREEN scored these poses with, and dropping
        # it would make the two pages describe different physics.
        if nac.NAC_DIST_MAX > _hi:
            ax[1].axhline(nac.NAC_DIST_MAX, color="#0f7a54", ls=":", lw=0.8,
                          alpha=0.55)
            ax[1].text(0.005, nac.NAC_DIST_MAX, " NAC window ", ha="left",
                       va="bottom", fontsize=6.5, color="#0f7a54", alpha=0.8,
                       transform=ax[1].get_yaxis_transform(),
                       bbox=dict(fc="white", ec="none", alpha=0.7, pad=1.0))
        ax[1].plot(t_d, dist * 10.0, lw=1.0, color="#b3261e")
        ax[1].set_ylabel("warhead–S$\\gamma$ (Å)")
        # BOXED, because the label sits at the band edge where the trace
        # spends most of its time -- unboxed it was overprinted by the red line
        # and read as noise.
        ax[1].text(0.995, _hi, f" attack ready ≤{_hi:.1f} Å", ha="right",
                   va="bottom", fontsize=7, color="#0f7a54",
                   transform=ax[1].get_yaxis_transform(),
                   bbox=dict(fc="white", ec="none", alpha=0.78, pad=1.0))
    ax[1].set_xlabel("time (ns)")
    for a in ax:
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
    # THE LENGTH AND THE FATE, ON THE FIGURE. Two traces that both end are not
    # the same result if one stopped because the molecule left at 5.2 ns and the
    # other because it was still there at the 10 ns cap.
    fate = _sweep_fate(ident)
    ttl = ident
    if fate.get("sweep_ps"):
        ns = fate["sweep_ps"] / 1000.0
        if fate.get("left") is True:
            ttl = f"{ident}  —  left the site at {ns:.1f} ns"
        elif fate.get("left") is False:
            ttl = f"{ident}  —  held to the {ns:.0f} ns cap"
        else:
            ttl = f"{ident}  —  {ns:.1f} ns"
    fig.suptitle(ttl, fontsize=8, color="#5b6b80", y=0.995)
    # And mark WHERE it stopped on the distance panel, so the end of the trace
    # reads as an event rather than as the edge of the plot.
    if fate.get("sweep_ps") and t_d is not None and len(t_d):
        end = float(t_d[-1])
        left = fate.get("left")
        if left is not None:
            col = "#b3261e" if left else "#0f7a54"
            ax[1].axvline(end, color=col, ls="--", lw=0.9, alpha=0.75)
            ax[1].text(end, 0.02, (" left " if left else " held "),
                       ha="right" if left else "right", va="bottom",
                       fontsize=6.5, color=col, rotation=90,
                       transform=ax[1].get_xaxis_transform(),
                       bbox=dict(fc="white", ec="none", alpha=0.75, pad=0.8))
    fig.tight_layout()
    buf = io.BytesIO(); fig.savefig(buf, format="png"); plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()



def _sweep_fate(ident: str) -> dict:
    """How long this mode ran and whether it left, from its own row.

    UNDER ADAPTIVE LENGTH THE RUN LENGTH IS A RESULT. A sweep now continues
    while the molecule is in the site and stops when it leaves, capped at 10 ns
    -- so "5.2 ns, left" and "10 ns, held" are different findings about
    different molecules, and a plot that shows only a trace makes them look like
    the same experiment run for different amounts of time by accident.

    Returns {} when the row says nothing, and the caller then annotates nothing
    rather than asserting a fate it does not have.
    """
    try:
        import pandas as _pd
        fs = [str(f) for f in rp.sweep_result_files()]
        if not fs:
            return {}
        d = _pd.concat([_pd.read_csv(f) for f in fs], ignore_index=True)
        hit = d[d.ident.astype(str) == ident]
        if not len(hit):
            return {}
        r = hit.iloc[-1]
        out = {}
        if _pd.notna(r.get("sweep_ps")):
            out["sweep_ps"] = float(r["sweep_ps"])
        if "left_site" in r.index and _pd.notna(r.get("left_site")):
            out["left"] = bool(r["left_site"])
        if "left_at_ps" in r.index and _pd.notna(r.get("left_at_ps")):
            out["left_at_ps"] = float(r["left_at_ps"])
        if "adaptive" in r.index and _pd.notna(r.get("adaptive")):
            out["adaptive"] = bool(r["adaptive"])
        return out
    except Exception:                                      # noqa: BLE001
        return {}


def _sweep_ps(ident: str, default: float = 1200.0) -> float:
    """This mode's sweep length, from its own row.

    Falls back to the CONFIG value rather than a literal, and only then to
    `default` -- so a run whose rows are missing still gets the length the
    campaign was configured with rather than one from two campaigns ago.
    """
    try:
        import glob as _g
        import pandas as _pd
        fs = [str(f) for f in rp.sweep_result_files()]   # ordered, one resolver
        if fs:
            d = _pd.concat([_pd.read_csv(f) for f in fs], ignore_index=True)
            hit = d[d.ident.astype(str) == ident]
            if len(hit) and _pd.notna(hit.iloc[-1].get("sweep_ps")):
                return float(hit.iloc[-1]["sweep_ps"])
    except Exception:                                      # noqa: BLE001
        pass
    try:
        from shared import target_config as _tc
        v = _tc.get("md.sweep_ps", default=None)
        if v:
            return float(v)
    except Exception:                                      # noqa: BLE001
        pass
    return default



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
    # THE MODE IDENT, WHATEVER THE WORKLIST CALLS IT.
    #
    # This read `wl.ident`, which in the older worklists WAS `<parent>_m<mode>`.
    # `sweep_supervisor`'s worklists put the MOLECULE in `ident` and the mode
    # ident in `task_id`, so every lookup here resolved a molecule instead of a
    # mode -- and `prank`, keyed on that, collapsed a molecule's several modes
    # onto whichever pose_rank happened to be last. The symptom was "no
    # finished run" for all 474 molecules while 16 trajectories sat on disk.
    #
    # Prefer `task_id`; fall back to `ident` so the old worklists still work.
    key = "task_id" if "task_id" in wl.columns else "ident"
    idents = list(wl[key].astype(str))[:args.limit]
    # mode ident -> pose_rank, so each mode's OWN trajectory is found. Keyed on
    # the MODE: keying on the molecule is what #23 was, and it handed every mode
    # of a molecule the same trajectory.
    prank = dict(zip(wl[key].astype(str), wl.pose_rank.astype(int))) \
        if "pose_rank" in wl.columns else {}
    n_mov = n_png = 0
    miss: dict = {}
    for i, ident in enumerate(idents, 1):
        parent = ident.rsplit("_m", 1)[0]
        pdb, png = OUT / f"{ident}.pdb", OUT / f"{ident}.png"
        src = OUT / f"{ident}.src"
        want_mov = not args.plots_only
        # RESOLVED BEFORE THE CACHE IS CONSULTED, because the cache is keyed on
        # WHICH trajectory drew the asset, not merely on the asset existing.
        rep = rep_dir(parent, prank.get(ident))
        if rep is None:
            miss[parent] = f"no finished {_sweep_ps(ident)/1000:.1f} ns run"
            continue
        # THE STALENESS BUG. Assets built before `rep_dir` took `pose_rank`
        # came from whichever run sorted first, so a multi-mode molecule's
        # plots and movie belonged to a sibling pose. Fixing the resolver did
        # not fix the files: t4_2f88a2f534fd_m1 was ranked on its own 0.255 nm
        # trace while its plot showed rank13's 0.857 nm, and nothing on the page
        # could have revealed the disagreement. @tt8804: "why is the selector
        # showing 0.255 nm max but the rmsd plots show max 0.857 nm".
        #
        # Recording the source path makes the cache self-invalidating: change
        # how a mode resolves and its assets rebuild without anyone having to
        # remember `--force`. A missing sidecar means the asset predates this,
        # i.e. exactly the suspect population, so it rebuilds.
        stale = (not src.is_file()) or src.read_text().strip() != str(rep)
        if not (args.force or args.plots_only or stale) \
                and png.is_file() and pdb.is_file():
            continue
        try:
            if want_mov and (args.force or stale or not pdb.is_file()):
                # THE SWEEP'S OWN LENGTH, NOT A LITERAL. This was
                # `total_ps=10_000.0`, written when the triage sweep was 10 ns.
                # The triage is 1.2 ns now (config `md.sweep_ps`), so every
                # movie built for this campaign would have carried a timeline
                # 8x too long -- frame 100 labelled 6.7 ns when it is 0.8 ns.
                # A number that was right when written and cannot announce that
                # it is not: `how_this_project_breaks` disguise #3, the same
                # shape as the timeout sized for a smaller pool (#19).
                mov.build_movie_pdb(rep, pdb, total_ps=_sweep_ps(ident))
                n_mov += 1
            if args.force or args.plots_only or stale or not png.is_file():
                b64 = plots(rep, ident)
                if b64:
                    png.write_bytes(base64.b64decode(b64)); n_png += 1
            # Written LAST, so an asset that failed half-way keeps its old
            # (or absent) sidecar and is retried rather than recorded as current.
            src.write_text(str(rep))
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
