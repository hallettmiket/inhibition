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
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import outputs as sout                  # noqa: E402
from shared import report_theme as rt               # noqa: E402

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


def figure(ident: str, s: dict, res: dict, er) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update(rt.MPL)

    panels = [k for k in ("rmsd", "dist", "mindist", "contacts") if k in s]
    fig, axes = plt.subplots(len(panels), 1, figsize=(9, 2.2 * len(panels)),
                             sharex=True)
    axes = np.atleast_1d(axes)
    titles = {"rmsd": "Ligand RMSD after superposing on the protein (nm)",
              "dist": "Warhead → Cys113 SG distance (nm)",
              "mindist": "Minimum ligand–protein distance (nm)",
              "contacts": "Ligand–protein contacts"}
    for ax, k in zip(axes, panels):
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
    axes[-1].set_xlabel("time (ns)")
    fig.suptitle(ident, x=0.005, ha="left", fontsize=11, weight="bold")
    fig.tight_layout()
    return er._png(fig, plt)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--candidate", required=True)
    ap.add_argument("--rep", default="rep1")
    ap.add_argument("--out", default=None)
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
    occ, cls, pose = PREDICTED.get(args.candidate, (None, "—", None))
    log.info("%s: %.1f ns, residence %.3f, dissociated=%s", args.candidate,
             res["length_ns"], res["residence_frac"], res["dissociated"])

    img = figure(args.candidate, s, res, er)
    verdict = ("Held" if not res["dissociated"] else
               f"Left at {res['left_at_ns']:.0f} ns")

    facts = [("molecule", args.candidate), ("class", cls),
             ("pose elevated", f"rank {pose}" if pose else "—"),
             ("BPMD occupancy", rt.num(occ, "{:.3f}")),
             ("trajectory", f"{res['length_ns']:.1f} ns, "
                            f"{res['n_frames']:,} frames")]
    body = [
        rt.masthead(
            f"{args.candidate} — 100 ns residence",
            f"{verdict}. Residence fraction {res['residence_frac']:.3f} "
            f"(frames with ligand RMSD ≤ {BOUND_NM} nm).",
            "MD-PRIORITY · 2.2.0 BORNITE", facts),
        rt.section("1", "What this molecule did"),
        f'<p>{rt.pill("Held" if not res["dissociated"] else "Left")} '
        f'Mean ligand RMSD {res["rmsd_mean_nm"]:.3f} nm, max '
        f'{res["rmsd_max_nm"]:.3f} nm, final {res["rmsd_final_nm"]:.3f} nm.</p>',
        f'<img src="{img}" style="max-width:100%">',
        rt.section("2", "Against what was predicted for it"),
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
    ]
    dest = Path(args.out) if args.out else OUT.dir / f"{args.candidate}.html"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(
        "<!doctype html><html lang='en'><head><meta charset='utf-8'>"
        "<meta name='viewport' content='width=device-width,initial-scale=1'>"
        f"<title>{args.candidate} — 100 ns residence</title>"
        f"<style>{rt.CSS}</style></head><body>\n"
        + "\n".join(body) + "\n</body></html>")
    print(f"{args.candidate}: {verdict}, residence {res['residence_frac']:.3f} "
          f"-> {dest}")

    row = {"ident": args.candidate, "bpmd_occupancy": occ, "warhead_class": cls,
           "pose_rank": pose, **{k: v for k, v in res.items() if k != "status"}}
    pd.DataFrame([row]).to_csv(OUT.write(f"mdprio_{args.candidate}", ".csv"),
                               index=False)


if __name__ == "__main__":
    main()
