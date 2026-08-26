#!/usr/bin/env python3
"""
Purpose: Compose the manuscript's "output of the Pin1 search" figure from real
         captures of the integration GUI.
Author:  Timothy Wu (with Claude Code)
Date:    2026-08-25
Input:   full-page PNG captures of the running integration app (see below)
Output:  versioned PDF + PNG under append_only/inhibition/00_outputs/artist/

THESE ARE CAPTURES, NOT A MOCK-UP. Every pixel below came from
`integration/app/app.py` served by `scripts/serve_gui.sh` and photographed with
headless chromium at device_scale_factor 2.5. Nothing is redrawn, retouched or
rearranged inside a panel; the only edits are the crop, the panel letters and
the margin annotations, all of which sit OUTSIDE the captured area. A figure of
an interface that has been redrawn is a figure of the redrawing, and this paper
is making a claim about what the interface actually shows.

To regenerate the captures:
    bash scripts/serve_gui.sh 8907                  # in one shell
    python scripts/manuscript_gui_capture.py        # in another
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
from PIL import Image

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from shared import outputs  # noqa: E402

TOPIC = "manuscript_figures"
INK, MUTED = "#1a1a1a", "#5c5c5c"
ACCENT = "#CC3311"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "text.color": INK,
    "figure.dpi": 200,
    "savefig.dpi": 400,
})


def crop(path: Path, box: tuple[int, int, int, int]) -> Image.Image:
    im = Image.open(path).convert("RGB")
    return im.crop(box)


def build(shots: Path):
    """Two captures stacked, with a right-hand margin reserved for callouts.

    The panels occupy PANEL_W of the figure width and the annotations live in
    the remaining margin with leader lines into the capture. Nothing is drawn
    over the interface itself -- an annotation that covers the thing it
    describes makes the figure unreadable and, worse, unfalsifiable.
    """
    a = crop(shots / "f2_shortlists_wide.png", (150, 250, 4060, 3380))
    b = crop(shots / "f2_dossier_wide.png", (150, 250, 4060, 2520))

    PANEL_W = 0.775                      # fraction of figure width
    wa, ha = a.size
    wb, hb = b.size
    ar_a = ha / wa                       # panel heights in figure-width units
    ar_b = (hb / wb)

    TITLE_H, CAP_H, GAP_H = 0.085, 0.062, 0.050   # in figure-width units
    total = (TITLE_H + CAP_H + ar_a * PANEL_W + GAP_H + CAP_H
             + ar_b * PANEL_W)

    fig_w = 7.4
    fig = plt.figure(figsize=(fig_w, fig_w * total))

    def f(y):                            # width-units from top -> figure frac
        return 1 - y / total

    fig.text(0.0, f(0.012), "What the Pin1 choreography actually outputs: "
             "four shortlists a human\nadjudicates, and never a winner",
             fontsize=10.6, fontweight="bold", va="top", ha="left",
             linespacing=1.35)

    y = TITLE_H

    def place(img, ar, letter, title, note):
        nonlocal y
        fig.text(0.0, f(y + 0.004), letter, fontsize=10, fontweight="bold",
                 va="top", ha="left")
        fig.text(0.026, f(y + 0.004), title, fontsize=9.2,
                 fontweight="bold", va="top", ha="left")
        fig.text(0.026, f(y + 0.026), note, fontsize=7.0, color=MUTED,
                 va="top", ha="left", linespacing=1.5)
        h = ar * PANEL_W
        ax = fig.add_axes([0.0, f(y + CAP_H + h), PANEL_W, h / total])
        ax.imshow(img, aspect="auto", interpolation="lanczos")
        # FREEZE THE LIMITS. annotate() with a data-coordinate xytext outside
        # the image autoscales the axes, so a second callout would measure a
        # different height from the first and land somewhere else entirely.
        # This bit me: anchors were resolving at roughly half their intended
        # depth, and the arrows still looked plausible, which is exactly the
        # failure mode docs/how_this_project_breaks.md is about.
        ax.set_xlim(0, img.size[0])
        ax.set_ylim(img.size[1], 0)
        ax.set_autoscale_on(False)
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_edgecolor("#c4c4c4"); sp.set_linewidth(0.7)
        y += CAP_H + h + GAP_H
        return ax

    axa = place(a, ar_a, "a",
                "The four approaches, side by side — and never merged",
                "Each column is ranked by its own metric. The banner above "
                "them carries the shared\ncontrol's verdict, so an ordering "
                "cannot be read without the warrant for it.")
    axb = place(b, ar_b, "b",
                "Every candidate opens into the evidence behind it",
                "One shared descriptor module for all four approaches, the "
                "docking rank, two\nindependent free-energy estimates, and "
                "the gate verdict attached to the rank.")

    # Callouts: anchor in the capture, text in the reserved right margin.
    def callout(ax, img, xf, yf, text, ytext_f):
        W, H = img.size
        ax.annotate(text, xy=(xf * W, yf * H),
                    xytext=(1.035 * W, ytext_f * H),
                    fontsize=6.8, color=ACCENT, ha="left", va="center",
                    linespacing=1.5, annotation_clip=False, zorder=6,
                    arrowprops=dict(arrowstyle="-|>", color=ACCENT, lw=0.85,
                                    shrinkA=1, shrinkB=2,
                                    connectionstyle="arc3,rad=0.16"))

    # Anchors are fractions of the CROPPED capture, calibrated by rendering a
    # ruled overlay of the crop rather than by eye.
    callout(axa, a, 0.60, 0.155, "the control's verdict\ntravels with the\n"
            "ranking it qualifies", 0.045)
    callout(axa, a, 0.40, 0.385, "the rebuild reports what\nit removed and\n"
            "what it promoted", 0.300)
    callout(axa, a, 0.845, 0.805, "T\u2084 ranks WITHIN warhead\nclass, so "
            "several\nrows share rank 1", 0.735)
    callout(axb, b, 0.090, 0.982, "every rank is stamped an\nordering the "
            "pipeline\nproduced, not binding", 0.910)

    written = []
    for suffix in (".pdf", ".png"):
        p = outputs.write_path("artist", TOPIC, "f2_gui_output", suffix)
        fig.savefig(p, bbox_inches="tight", facecolor="white")
        written.append(p)
    plt.close(fig)
    return written


if __name__ == "__main__":
    shots = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    if shots is None or not shots.is_dir():
        raise SystemExit("usage: manuscript_gui_figure.py <dir-with-captures>")
    for p in build(shots):
        print("WROTE", p)
