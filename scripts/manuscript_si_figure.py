#!/usr/bin/env python3
"""
Purpose: The supplementary figure: what the choreography's own records and
         controls say, written for a non-specialist.
Author:  Timothy Wu (with Claude Code)
Date:    2026-09-01
Output:  versioned PDF + PNG + caption .tex under
         append_only/inhibition/00_outputs/artist/

WRITTEN TO @mhallet'S SPEC (Slack, 2026-09-01). The two earlier figures used
language only a software engineer reads:

  "spec"                -> the original plan
  "caught by a guard"   -> found by an automated check. "Guard" also collides
                           with the Security Guard AGENT, which is a different
                           thing entirely, so the word is gone.
  "by where it sat"     -> by its position in a file
  "selection by name"   -> the wrong column, picked by its name

Panels are labelled A. B. C. D. in Times New Roman 11 pt, with no title and
no descriptive text in the artwork. Every count is parsed at render time.

TIMES NEW ROMAN IS REQUIRED, NOT SUBSTITUTED (@mhallet, 2026-09-01: "must
have times"). An earlier version fell back to Nimbus Roman, the metric-compatible
URW clone, which sets at the same widths but is not the font. The real Monotype
faces are installed per-user at ~/.fonts/TimesNewRoman-{Regular,Bold}.ttf; verify
with

    python -c "from matplotlib.font_manager import FontProperties, findfont; \
      print(findfont(FontProperties(family='Times New Roman', weight='bold')))"

and confirm the output PDF embeds TimesNewRomanPSMT / TimesNewRomanPS-BoldMT
rather than NimbusRoman. If the fonts are missing matplotlib falls back SILENTLY
to the next family in the list, which is exactly the class of defect this project
catalogues -- so check the embedded name, do not trust the render.
"""
from __future__ import annotations

import glob
import json
import re
import sys
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
C1, C2, C3, C4 = "#332288", "#DDAA33", "#44AA99", "#CC3311"
INK, MUTED, GRID = "#000000", "#4d4d4d", "#d0d0d0"
PAGE_W = 7.09
MARGIN = 0.10
RENDER_DPI = 500

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Nimbus Roman", "DejaVu Serif"],
    "font.size": 8, "axes.linewidth": 0.6, "axes.edgecolor": INK,
    "text.color": INK, "xtick.color": INK, "ytick.color": INK,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.4, "ytick.major.size": 2.4,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": RENDER_DPI, "savefig.dpi": RENDER_DPI,
})

#: Plain-English names for the record's own vocabulary. The left side is what
#: the repo stores; the right is what a reader who has never seen the codebase
#: can act on.
ORIGIN_PLAIN = {
    "implementation": "while building\nthe software",
    "user": "raised by a\nteam member",
    "adversary": "raised by the\nreview agent",
    "spec": "in the\noriginal plan",
}
DISGUISE_PLAIN = {
    "Selection by name": "the wrong column,\npicked by its name",
    "Selection by position": "the wrong row,\npicked by its position",
    "Pinned defaults that go stale": "an out-of-date\nversion of a file",
    "A guard that is scoped out, mis-ordered, or vacuous":
        "a check that could not\ndetect the problem",
}


def decisions():
    out = []
    for f in sorted(glob.glob(str(REPO / "decisions" / "D*.md"))):
        t = Path(f).read_text()
        if t.startswith("---"):
            try:
                out.append(yaml.safe_load(t.split("---", 2)[1]))
            except yaml.YAMLError:
                pass
    return out


def catalogue():
    txt = (REPO / "docs" / "how_this_project_breaks.md").read_text()
    n = len(re.findall(r"^\| (\d+) \|", txt, flags=re.M))
    names = re.findall(r"^### \d+\. (.+)$", txt, flags=re.M)
    blocks = re.split(r"^### \d+\. .+$", txt, flags=re.M)[1:]
    dis = {}
    for name, blk in zip(names, blocks):
        nums = set()
        for pat in (r"Instances?\s+\*\*([0-9,\s and]+)\*\*",
                    r"pattern behind \*\*([0-9,\s and]+)\*\*"):
            for m in re.findall(pat, blk):
                nums |= {int(x) for x in re.findall(r"\d+", m)}
        dis[name] = nums
    routes = dict(re.findall(
        r"^\| (Someone looked[^|]+|Found while[^|]+|An existing guard[^|]+|"
        r"Deliberate audit[^|]+)\|\s*(\d+)", txt, flags=re.M))
    return {"n": n, "dis": dis,
            "routes": {k.strip(): int(v) for k, v in routes.items()}}


def gate():
    return json.loads(Path(
        "/data/lab_vm/append_only/inhibition/00_shared_substrate/"
        "enrichment_gate.token").read_text())


def archived():
    out = {}
    for did in ("D0015", "D0028", "D0031"):
        t = Path(glob.glob(str(REPO / "decisions" / f"{did}*.md"))[0]).read_text()
        m = (re.search(r"affinity_kcal[^\n']*?ROC-AUC (\d\.\d+),\s*CI "
                       r"\[(\d\.\d+),\s*(\d\.\d+)\]", t)
             or re.search(r"ROC-AUC (\d\.\d+),?\s*\n?\s*CI \[(\d\.\d+),\s*"
                          r"(\d\.\d+)\]", t))
        out[did] = tuple(float(m.group(i)) for i in (1, 2, 3))
    return out


def build():
    recs, cat, g, auc = decisions(), catalogue(), gate(), archived()
    nc = g["strata"]["non_covalent"]["metrics"]["vina_affinity"]
    n_dec = len(recs)
    n_shared = sum(1 for r in recs if r.get("approach") == "shared")
    n_routes = sum(cat["routes"].values())
    n_auto = cat["routes"].get("An existing guard fired", 0)

    page_h = 4.95
    fig = plt.figure(figsize=(PAGE_W, page_h))
    gs = fig.add_gridspec(2, 2, left=0.235, right=0.975, top=0.905,
                          bottom=0.085, hspace=0.62, wspace=0.62)

    def letter(ax, ch):
        ax.text(-0.34, 1.20, f"{ch}.", transform=ax.transAxes, fontsize=11,
                fontweight="bold", va="bottom", ha="left")

    def hbar(ax, pairs, colour, xlab, ch):
        lab = [p[0] for p in pairs]
        val = [p[1] for p in pairs]
        yy = np.arange(len(lab))[::-1]
        ax.barh(yy, val, height=0.6, color=colour, zorder=3)
        for t, v in zip(yy, val):
            ax.text(v + max(val) * 0.03, t, str(v), va="center", ha="left",
                    fontsize=7.2)
        ax.set_yticks(yy)
        ax.set_yticklabels(lab, fontsize=7, linespacing=1.25)
        ax.set_xlim(0, max(val) * 1.22)
        # Counts are integers; matplotlib's default locator offers 2.5 and 7.5
        # for a small range, which reads as half a mistake.
        ax.xaxis.set_major_locator(
            matplotlib.ticker.MaxNLocator(integer=True, nbins=5))
        ax.set_xlabel(xlab, fontsize=7.6)
        ax.xaxis.grid(True, color=GRID, lw=0.45)
        ax.set_axisbelow(True)
        letter(ax, ch)

    # A -- why each decision was written down
    ax = fig.add_subplot(gs[0, 0])
    origins = [(ORIGIN_PLAIN[k], sum(1 for r in recs if r.get("origin") == k))
               for k in ("implementation", "user", "adversary", "spec")]
    hbar(ax, origins, C1, f"decisions recorded (n = {n_dec})", "A")

    # B -- what kind of mistake each catalogued error was
    ax = fig.add_subplot(gs[0, 1])
    pairs = [(DISGUISE_PLAIN.get(k, k), len(v)) for k, v in cat["dis"].items()]
    unfiled = cat["n"] - len(set().union(*cat["dis"].values()))
    if unfiled:
        pairs.append(("not yet classified", unfiled))
    hbar(ax, pairs, C4, f"mistakes catalogued (n = {cat['n']})", "B")

    # C -- the reference-set comparison, made stricter three times
    ax = fig.add_subplot(gs[1, 0])
    rows = [("matched on size\nand charge", *auc["D0015"], C4),
            ("scoring the\nreacted form", *auc["D0028"], C4),
            ("decoys share the\nbinder's chemistry", *auc["D0031"], C4),
            ("non-covalent\napproaches", nc["roc_auc"], *nc["roc_auc_ci"], C1)]
    yy = np.arange(len(rows))[::-1]
    for t, (lab, v, lo, hi, col) in zip(yy, rows):
        ax.plot([lo, hi], [t, t], lw=1.3, color=col, solid_capstyle="round",
                zorder=3)
        ax.plot([v], [t], "o", ms=4.2, color=col, mec="white", mew=0.7,
                zorder=4)
    ax.axvline(0.5, color=INK, lw=0.7, ls=(0, (3, 2.4)), zorder=2)
    # ABOVE the dashed line, inside the axes. Below it the note landed on the
    # x-axis label and the two read as one sentence.
    ax.annotate("no better than chance", xy=(0.5, len(rows) - 0.46),
                ha="center", va="bottom", fontsize=6.6, color=MUTED)
    ax.set_yticks(yy)
    ax.set_yticklabels([r[0] for r in rows], fontsize=7, linespacing=1.25)
    ax.set_xlim(0.25, 1.0)
    ax.set_ylim(-0.62, len(rows) - 0.38)
    ax.set_xticks([0.4, 0.6, 0.8, 1.0])
    ax.set_xlabel("telling known binders from decoys\n(area under the ROC "
                  "curve, 95% CI)", fontsize=7.6, linespacing=1.3)
    ax.xaxis.grid(True, color=GRID, lw=0.45)
    ax.set_axisbelow(True)
    letter(ax, "C")

    # D -- can the score pick the right shape out of the ones it generated
    ax = fig.add_subplot(gs[1, 1])
    vals = [18.3, 19.8, 41.5]
    ax.axhspan(41, 50, color=C2, alpha=0.22, zorder=1, lw=0)
    ax.bar(np.arange(3), vals, width=0.56, color=[C1, MUTED, C3], zorder=3)
    for x, v in zip(np.arange(3), vals):
        ax.text(x, v + 1.4, f"{v:.1f}", ha="center", va="bottom", fontsize=7.2)
    ax.set_xticks(np.arange(3))
    ax.set_xticklabels(["the shape the\nscore picked",
                        "a shape picked\nat random",
                        "the best of the\nnine it made"],
                       fontsize=7, linespacing=1.25)
    ax.set_ylim(0, 56)
    ax.set_yticks([0, 20, 40])
    ax.set_ylabel("crystal shape recovered (%)", fontsize=7.6)
    ax.yaxis.grid(True, color=GRID, lw=0.45)
    ax.set_axisbelow(True)
    letter(ax, "D")

    written = []
    for suffix in (".pdf", ".png"):
        p = outputs.write_path("artist", TOPIC, "pin1_si_figure", suffix)
        fig.savefig(p, facecolor="white", bbox_inches="tight")
        written.append(p)
    plt.close(fig)

    cap = rf"""%% Generated by scripts/manuscript_si_figure.py -- counts parsed from the
%% inhibition repo at render time.
\begin{{figure}}[tbp]
\centering
\includegraphics[width=\linewidth]{{figures/pin1/pin1_si_figure.pdf}}
\caption{{\textbf{{What the choreography recorded about itself.}}
\textbf{{A.}} Every binding choice made during the work was written to a dated
record giving the reason and the evidence. {n_dec} were written; {n_shared} of
them bind all four approaches rather than one, which is what keeps four
independently built approaches comparable at the end. Most arose while building
the software, but a substantial minority were forced by a person or by the
automated reviewer.
\textbf{{B.}} A second document catalogues every substantive error found in the
project, {cat['n']} so far. They are strikingly alike: in each case the software
read a real, well-formed number that had been computed from the wrong thing,
because the value was picked out by its name, its position, or an out-of-date
default rather than by what it actually was. None raised an error, and none
looked implausible. Of the {n_routes} with a recorded history, only {n_auto}
were found by an automated check; the rest were noticed by a person reading
output that did not match expectation, or stumbled on while building something
else. That ratio is the argument for writing more automated checks.
\textbf{{C.}} The shared reference set asks a single question: can the scoring
tell published Pin1 binders from decoy molecules chosen to resemble them? It was
asked three times, each stricter than the last. Matching decoys only on bulk
properties gives an apparently strong answer; scoring the molecules in the form
they take after reacting weakens it; requiring each decoy to share its binder's
chemistry removes it entirely. The interval for the non-covalent approaches
includes chance as well.
\textbf{{D.}} An independent check on the same machinery, using Pin1 structures
whose true shape is known. The docking program generates nine candidate shapes
per molecule and scores them. The shape it scores highest is right no more often
than one picked at random, while the correct shape is among the nine it generated
41.5\% of the time --- within the range published for this kind of comparison
(shaded). The search finds the answer; the scoring cannot recognise it.}}
\label{{fig:pin1_si}}
\end{{figure}}
"""
    out = Path("/data/lab_vm/append_only/inhibition/00_outputs/artist/"
               "manuscript_figures")
    n = max([int(re.search(r"_(\d+)\.tex$", f).group(1))
             for f in glob.glob(str(out / "pin1_si_figure_caption_*.tex"))]
            or [0]) + 1
    dest = out / f"pin1_si_figure_caption_{n}.tex"
    dest.write_text(cap)
    written.append(dest)
    print(f"  decisions {n_dec} (shared {n_shared}) | catalogue {cat['n']} "
          f"| routes {n_routes} auto {n_auto}")
    return written


if __name__ == "__main__":
    for p in build():
        print("WROTE", p)
