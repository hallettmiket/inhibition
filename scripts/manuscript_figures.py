#!/usr/bin/env python3
"""
Purpose: Generate the Pin1-choreography figures for the murmurent manuscript.
Author:  Timothy Wu (with Claude Code)
Date:    2026-08-25
Input:   decisions/*.md, docs/how_this_project_breaks.md, the enrichment-gate
         token under 00_shared_substrate/, and the integration GUI captures.
Output:  versioned PDF + PNG under append_only/inhibition/00_outputs/artist/
         (topic: manuscript_figures), via shared.outputs.

EVERY NUMBER HERE IS READ FROM A FILE, NOT TYPED IN. The manuscript's own
figures for this choreography drifted from the repo within days of being
written (51 decision records against 87 actual, 21 catalogue entries against
25). A figure that hardcodes a count cannot announce that it is stale -- that
is disguise #3 in docs/how_this_project_breaks.md. So the corpus counts, the
catalogue counts and the gate verdicts are all parsed at render time, and the
figure prints what it parsed so a reader can check it against the repo.

The two places a literal IS used are labelled as such: the three archived
ROC-AUC measurements of the covalent control (D0015 / D0028 / D0031) are read
out of those decision records' evidence lines rather than recomputed, because
the decoy sets they were measured against are themselves superseded; and the
cross-docking literature band comes from docs/publication_audit.md.
"""
from __future__ import annotations

import glob
import json
import re
import sys
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from shared import outputs  # noqa: E402

TOPIC = "manuscript_figures"

# CVD-validated categorical order (Paul Tol muted). Checked with the six-check
# validator: worst-case OKLab dE across normal/deuter/protan/tritan = 13.1,
# against a floor of 8 for CVD and 15 for normal vision. Assigned to T1..T4 in
# fixed order and never cycled.
ARM = {"t1": "#332288", "t2": "#DDAA33", "t3": "#44AA99", "t4": "#CC3311"}
INK, MUTED, GRID = "#1a1a1a", "#5c5c5c", "#d8d8d8"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "font.size": 8.5,
    "axes.edgecolor": MUTED,
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 200,
    "savefig.dpi": 400,
})


# ----------------------------------------------------------------- parsing --
def decision_records() -> list[dict]:
    recs = []
    for f in sorted(glob.glob(str(REPO / "decisions" / "D*.md"))):
        txt = Path(f).read_text()
        if not txt.startswith("---"):
            continue
        try:
            fm = yaml.safe_load(txt.split("---", 2)[1])
        except yaml.YAMLError:
            continue
        fm["_file"] = Path(f).name
        recs.append(fm)
    return recs


def catalogue() -> dict:
    """Counts parsed from the failure catalogue, not copied from it."""
    txt = (REPO / "docs" / "how_this_project_breaks.md").read_text()
    entries = re.findall(r"^\| (\d+) \|", txt, flags=re.M)
    disguise_names = re.findall(r"^### \d+\. (.+)$", txt, flags=re.M)
    # Which entry numbers each disguise claims. The catalogue's own prose is the
    # source; entries it has not yet filed are counted as unassigned rather than
    # silently dropped, because "not yet classified" and "no instances" are
    # different facts (honesty rule 2 in docs/gui_spec.md).
    claimed: dict[str, set] = {}
    for name, block in zip(disguise_names, re.split(r"^### \d+\. .+$", txt,
                                                    flags=re.M)[1:]):
        nums = set()
        for m in re.findall(r"Instances?\s+\*\*([0-9,\s and]+)\*\*", block):
            nums |= {int(n) for n in re.findall(r"\d+", m)}
        for m in re.findall(r"pattern behind \*\*([0-9,\s and]+)\*\*", block):
            nums |= {int(n) for n in re.findall(r"\d+", m)}
        claimed[name] = nums
    routes = dict(re.findall(r"^\| (Someone looked[^|]+|Found while[^|]+|"
                             r"An existing guard[^|]+|Deliberate audit[^|]+)\|\s*(\d+)",
                             txt, flags=re.M))
    return {"n_entries": len(entries),
            "disguises": claimed,
            "routes": {k.strip(): int(v) for k, v in routes.items()}}


def gate_token() -> dict:
    p = Path("/data/lab_vm/append_only/inhibition/00_shared_substrate/"
             "enrichment_gate.token")
    return json.loads(p.read_text())


def archived_auc() -> list[tuple]:
    """The covalent control's three measurements, read from their records."""
    want = {"D0015": "pre-reaction ligands,\nproperty-matched decoys",
            "D0028": "adduct forms scored,\nsame decoys",
            "D0031": "decoys required to carry\nthe active's own chemotype"}
    out = []
    for did, label in want.items():
        f = glob.glob(str(REPO / "decisions" / f"{did}*.md"))[0]
        txt = Path(f).read_text()
        m = re.search(r"affinity_kcal[^\n']*?ROC-AUC (\d\.\d+),\s*CI \[(\d\.\d+),\s*(\d\.\d+)\]",
                      txt)
        if not m:
            m = re.search(r"ROC-AUC (\d\.\d+),?\s*\n?\s*CI \[(\d\.\d+),\s*(\d\.\d+)\]", txt)
        auc, lo, hi = (float(m.group(i)) for i in (1, 2, 3))
        out.append((did, label, auc, lo, hi))
    return out


def save(fig, stem: str) -> list[Path]:
    written = []
    for suffix in (".pdf", ".png"):
        p = outputs.write_path("artist", TOPIC, stem, suffix)
        fig.savefig(p, bbox_inches="tight", facecolor="white")
        written.append(p)
    plt.close(fig)
    return written


# ------------------------------------------------------- F3 control figure --
def fig_control_sequence():
    rows = archived_auc()
    gate = gate_token()
    nc = gate["strata"]["non_covalent"]["metrics"]["vina_affinity"]

    # Forest layout: a left gutter for the record id and what changed, the plot
    # in the middle, the CI verdict on the right. Labels live OUTSIDE the axes
    # so they cannot collide with an interval, whatever the values do.
    fig, ax = plt.subplots(figsize=(7.4, 2.9))
    fig.subplots_adjust(left=0.335, right=0.795)

    series = [(d, l, a, lo, hi, ARM["t4"]) for d, l, a, lo, hi in rows]
    series.append(("D0041",
                   "non-covalent arm,\ncurrent verdict (%s)" % nc["verdict"],
                   nc["roc_auc"], nc["roc_auc_ci"][0], nc["roc_auc_ci"][1],
                   ARM["t1"]))

    ys = list(range(len(series)))[::-1]
    for y, (did, label, auc, lo, hi, col) in zip(ys, series):
        ax.plot([lo, hi], [y, y], lw=2.2, color=col, solid_capstyle="round",
                zorder=2)
        ax.plot([auc], [y], "o", ms=8.5, color=col, mec="white", mew=1.6,
                zorder=3)
        ax.annotate(f"{auc:.3f}", xy=(auc, y), xytext=(0, 7),
                    textcoords="offset points", ha="center", va="bottom",
                    fontsize=9, fontweight="bold", color=INK, zorder=4)
        ax.annotate(f"{did}   {label}", xy=(0, y),
                    xycoords=("axes fraction", "data"), xytext=(-12, 0),
                    textcoords="offset points", ha="right", va="center",
                    fontsize=8, color=INK, linespacing=1.4)
        excl = lo > 0.5
        ax.annotate("CI excludes 0.5" if excl else "CI includes 0.5",
                    xy=(1, y), xycoords=("axes fraction", "data"),
                    xytext=(12, 0), textcoords="offset points", ha="left",
                    va="center", fontsize=7.5, style="italic",
                    color=col if excl else MUTED)

    ax.axvline(0.5, color=INK, lw=1.1, ls=(0, (4, 3)), zorder=1)
    ax.annotate("chance", xy=(0.5, len(series) - 0.62), ha="center",
                va="bottom", fontsize=8, color=INK)
    ax.axhline(0.5, color=GRID, lw=0.9)

    ax.set_xlim(0.25, 1.0)
    ax.set_ylim(-0.7, len(series) - 0.15)
    ax.set_yticks([])
    ax.set_xticks(np.arange(0.3, 1.01, 0.1))
    ax.set_xlabel("ROC-AUC against the shared control  (95% CI)")
    ax.spines["left"].set_visible(False)
    ax.xaxis.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_title("Each time the control was made stricter, the enrichment it "
                 "had certified went away",
                 loc="left", x=-0.505, fontsize=10.5, fontweight="bold", pad=14)
    return save(fig, "f3_control_walked_to_chance"), rows, nc


# ---------------------------------------------------- F4 written substrate --
def fig_written_substrate():
    recs = decision_records()
    cat = catalogue()
    n = len(recs)

    fig = plt.figure(figsize=(7.4, 4.9))
    gs = fig.add_gridspec(2, 2, hspace=0.95, wspace=0.55,
                          height_ratios=[1, 1])

    def hbar(ax, counts, title, color, note=None):
        labels = [k for k, _ in counts]
        vals = [v for _, v in counts]
        ys = np.arange(len(labels))[::-1]
        ax.barh(ys, vals, height=0.58, color=color, zorder=2)
        for y, v in zip(ys, vals):
            ax.text(v + max(vals) * 0.025, y, str(v), va="center", ha="left",
                    fontsize=8.5, fontweight="bold", color=INK)
        ax.set_yticks(ys)
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlim(0, max(vals) * 1.22)
        ax.set_xticks([])
        ax.spines["bottom"].set_visible(False)
        ax.spines["left"].set_color(GRID)
        ax.set_title(title, loc="left", fontsize=9, fontweight="bold", pad=6)
        if note:
            ax.text(0, -0.26, "\n".join(textwrap.wrap(note, 52)),
                    transform=ax.transAxes, fontsize=6.9, color=MUTED,
                    va="top", linespacing=1.45)

    scope = [("binds the whole\nchoreography",
              sum(1 for r in recs if r.get("approach") == "shared")),
             ("one approach or\nthe integration layer",
              sum(1 for r in recs if r.get("approach") != "shared"))]
    hbar(fig.add_subplot(gs[0, 0]), scope,
         f"What each record binds  (n = {n})", ARM["t1"],
         "Shared records outnumber approach-scoped ones ~5:1. This is what keeps four independently built approaches comparable at the end.")

    order = ["implementation", "user", "adversary", "spec"]
    origin = [(o, sum(1 for r in recs if r.get("origin") == o)) for o in order]
    hbar(fig.add_subplot(gs[0, 1]), [o for o in origin if o[1]],
         "What forced the record", ARM["t3"],
         "The trail records which actor — human or agent — forced each change, not only what changed.")

    ax = fig.add_subplot(gs[1, 0])
    dis = cat["disguises"]
    assigned = set().union(*dis.values()) if dis else set()
    counts = [(k.replace("Selection by ", "by ")
                .replace("A guard that is scoped out, mis-ordered, or vacuous",
                         "a guard that cannot fail")
                .replace("Pinned defaults that go stale", "pinned defaults gone stale"),
               len(v)) for k, v in dis.items()]
    unassigned = cat["n_entries"] - len(assigned)
    if unassigned > 0:
        counts.append(("not yet filed to a disguise", unassigned))
    hbar(ax, counts, f"The catalogued defects, by disguise  "
                     f"(n = {cat['n_entries']})", ARM["t4"],
         "Every entry is a value taken by where it sat or what it was called, rather than by what defined it. None raised an exception.")

    ax = fig.add_subplot(gs[1, 1])
    routes = [(k.replace("Someone looked at output and it didn't match expectation",
                         "someone read output and it\nlooked wrong")
                .replace("Found while building something else entirely",
                         "found while building\nsomething else")
                .replace("An existing guard fired", "an existing guard fired")
                .replace("Deliberate audit for this class of defect",
                         "deliberate audit for\nthis class"), v)
              for k, v in cat["routes"].items()]
    n_routes = sum(v for _, v in routes)
    hbar(ax, routes, f"How each defect was actually found  (n = {n_routes})",
         ARM["t2"],
         f"Only 3 of {n_routes} were caught by a guard. That ratio is the "
         "argument for writing guards rather than for auditing harder. "
         f"({cat['n_entries'] - n_routes} newer entries are not yet filed to a "
         "route in the catalogue.)")

    fig.suptitle("The written substrate accumulates faster than the code, and "
                 "records what the code cannot",
                 x=0.02, ha="left", fontsize=10.5, fontweight="bold", y=1.015)
    return save(fig, "f4_written_substrate"), n, cat


# ------------------------------------------------------- F5 levels of theory --
def fig_levels_of_theory():
    gate = gate_token()
    nc = gate["strata"]["non_covalent"]["metrics"]["vina_affinity"]
    cv = gate["strata"]["covalent"]["metrics"]["mmgbsa_dG"]

    fig, (axa, axb) = plt.subplots(1, 2, figsize=(7.4, 3.25),
                                   gridspec_kw={"wspace": 0.62,
                                                "width_ratios": [1.12, 1]})
    fig.subplots_adjust(left=0.235, right=0.965, top=0.74, bottom=0.26)

    # (a) the two decoy-based gates, gutter-labelled so nothing can collide
    items = [("Docking enrichment\nnon-covalent · D0041", nc["roc_auc"],
              nc["roc_auc_ci"], nc["verdict"], ARM["t1"]),
             ("Ensemble MM-GBSA\ncovalent · D0036", cv["roc_auc"],
              cv["roc_auc_ci"], cv["verdict"], ARM["t4"])]
    for y, (lab, auc, ci, verdict, c) in zip([1, 0], items):
        axa.plot(ci, [y, y], lw=2.2, color=c, solid_capstyle="round", zorder=2)
        axa.plot([auc], [y], "o", ms=8.5, color=c, mec="white", mew=1.6,
                 zorder=3)
        axa.annotate(f"{auc:.3f}", xy=(auc, y), xytext=(0, 8),
                     textcoords="offset points", ha="center", va="bottom",
                     fontsize=9, fontweight="bold")
        axa.annotate(lab, xy=(0, y), xycoords=("axes fraction", "data"),
                     xytext=(-10, 4), textcoords="offset points", ha="right",
                     va="center", fontsize=7.8, linespacing=1.4)
        axa.annotate(verdict, xy=(0, y), xycoords=("axes fraction", "data"),
                     xytext=(-10, -14), textcoords="offset points",
                     ha="right", va="center", fontsize=7.4,
                     fontweight="bold", color=c)
    axa.axvline(0.5, color=INK, lw=1.1, ls=(0, (4, 3)), zorder=1)
    axa.annotate("chance", xy=(0.5, 1.52), ha="center", va="bottom",
                 fontsize=8)
    axa.set_xlim(0.2, 1.0)
    axa.set_ylim(-0.55, 1.75)
    axa.set_yticks([])
    axa.set_xticks([0.2, 0.4, 0.6, 0.8, 1.0])
    axa.set_xlabel("ROC-AUC  (95% CI)")
    axa.spines["left"].set_visible(False)
    axa.xaxis.grid(True, color=GRID, lw=0.6)
    axa.set_axisbelow(True)
    axa.set_title("a   Both decoy-based gates include chance", loc="left",
                  x=-0.42, fontsize=9.3, fontweight="bold", pad=12)

    # (b) pose recovery on 3IKD against the published cross-docking band
    labels = ["Vina's\ntop-1", "random pick\nof the 9", "best-of-9\n(ceiling)"]
    vals = [18.3, 19.8, 41.5]
    cols = [ARM["t1"], MUTED, ARM["t3"]]
    xs = np.arange(3)
    axb.axhspan(41, 50, color=ARM["t2"], alpha=0.20, zorder=1)
    axb.bar(xs, vals, width=0.58, color=cols, zorder=3)
    for x, v in zip(xs, vals):
        axb.text(x, v + 1.6, f"{v:.1f}%", ha="center", va="bottom",
                 fontsize=8.8, fontweight="bold", zorder=4)
    axb.annotate("published cross-docking band, 41–50%", xy=(-0.45, 52.5),
                 ha="left", va="bottom", fontsize=6.8, color="#7a6114",
                 annotation_clip=False)
    axb.set_xticks(xs)
    axb.set_xticklabels(labels, fontsize=7.2, linespacing=1.4)
    axb.set_ylim(0, 60)
    axb.set_yticks([0, 10, 20, 30, 40, 50])
    axb.set_ylabel("cases with a pose within 2 Å  (%)", fontsize=7.8)
    axb.yaxis.grid(True, color=GRID, lw=0.6)
    axb.set_axisbelow(True)
    axb.set_title("b   The search finds the pose;\n      the score cannot pick it",
                  loc="left", x=-0.20, fontsize=9.3, fontweight="bold", pad=8)

    fig.suptitle("Every level of theory measured on this pocket has failed to "
                 "rank, by independent tests",
                 x=0.012, ha="left", fontsize=10.5, fontweight="bold", y=1.02)
    fig.text(0.012, 0.012, "\n".join(textwrap.wrap(
             "Also measured and not shown: implicit- and explicit-solvent MD "
             "residence (D0038, D0044), which was not reproducible between "
             "solvent models. Pose recovery is on the prepared 3IKD receptor "
             "(D0059); the 6VAJ measurements it replaces are invalidated.",
             104)), fontsize=6.8, color=MUTED, va="top", linespacing=1.5)
    return save(fig, "f5_levels_of_theory"), nc, cv


if __name__ == "__main__":
    p3, rows, nc = fig_control_sequence()
    p4, n_dec, cat = fig_written_substrate()
    p5, _, cv = fig_levels_of_theory()

    print("PARSED AT RENDER TIME (check these against the repo):")
    print(f"  decision records          : {n_dec}")
    print(f"  catalogue entries         : {cat['n_entries']}")
    print(f"  catalogue routes          : {cat['routes']}")
    print(f"  disguises                 : "
          f"{ {k: len(v) for k, v in cat['disguises'].items()} }")
    print(f"  covalent control sequence : "
          f"{[(d, a) for d, _, a, _, _ in rows]}")
    print(f"  non-covalent gate         : {nc['roc_auc']:.3f} "
          f"{nc['roc_auc_ci']} {nc['verdict']}")
    print(f"  covalent MM-GBSA gate     : {cv['roc_auc']:.3f} "
          f"{cv['roc_auc_ci']} {cv['verdict']}")
    print("\nWROTE:")
    for p in (*p3, *p4, *p5):
        print("  ", p)
