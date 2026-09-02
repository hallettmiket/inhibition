#!/usr/bin/env python3
"""
Purpose: The main-text full-page figure: the choreography's pipeline and the
         GUI that presents its output, panels A-D.
Author:  Timothy Wu (with Claude Code)
Date:    2026-09-01
Input:   Hagar Emam's flowchart PDF and captures of the running GUI.
Output:  versioned PDF + PNG under append_only/inhibition/00_outputs/artist/
         plus the caption as a .tex fragment beside it.

WRITTEN TO @mhallet'S SPEC (Slack, 2026-09-01):
  - one full-page figure carrying both the pipeline diagram and the GUI summary
  - panels labelled A. B. C. D. in the manuscript's own style: Times New Roman
    11 pt, bold, with the full stop
  - NO figure title, NO annotation arrows, NO per-panel descriptive text in the
    image. All of it moves to the caption, which is emitted separately as a
    .tex fragment rather than baked into the artwork -- a caption rendered into
    a PDF cannot be edited in Overleaf, and this manuscript is Overleaf-synced.

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

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from shared import outputs  # noqa: E402

TOPIC = "manuscript_figures"
PAGE_W, PAGE_H = 8.27, 11.69          # A4 portrait
MARGIN = 0.55
CW = PAGE_W - 2 * MARGIN
RENDER_DPI = 500

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Nimbus Roman", "DejaVu Serif"],
    "figure.dpi": RENDER_DPI,
    "savefig.dpi": RENDER_DPI,
})


def place(fig, img: Image.Image, x, y, w, letter, letter_x=None):
    """One panel: the capture, and its letter set clear of the artwork.

    `letter_x` overrides where the letter sits. A panel narrower than the text
    block is centred, but its letter still belongs on the left margin with
    every other letter -- a column of A. B. C. that steps in and out is harder
    to scan than the panels it labels.
    """
    h = w * img.size[1] / img.size[0]
    ax = fig.add_axes([x / PAGE_W, 1 - (y + h) / PAGE_H, w / PAGE_W, h / PAGE_H])
    ax.imshow(img, aspect="auto", interpolation="lanczos")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    # The letter sits ABOVE the panel, not into it. Anchored va="top" at the
    # panel's own y it hangs down over the artwork, and "B." landed on the
    # first column heading.
    lx = MARGIN if letter_x is None else letter_x
    fig.text(lx / PAGE_W, 1 - (y - 0.055) / PAGE_H, f"{letter}.",
             fontsize=11, fontweight="bold", va="bottom", ha="left")
    return h


def build(shots: Path):
    flow = Image.open(shots / "flowchart.png").convert("RGB")
    short = Image.open(shots / "panel_shortlists.png").convert("RGB")
    cand = Image.open(shots / "panel_candidate.png").convert("RGB")
    pose = Image.open(shots / "pose_interaction.png").convert("RGB")

    fig = plt.figure(figsize=(PAGE_W, PAGE_H))
    y = MARGIN + 0.20

    # A -- the pipeline. Narrower than the text block so the boxes do not set
    # smaller than the caption they are read alongside.
    fw = 6.05
    y += place(fig, flow, MARGIN + (CW - fw) / 2, y, fw, "A",
               letter_x=MARGIN) + 0.30

    # B -- the four shortlists, side by side
    y += place(fig, short, MARGIN, y, CW, "B", letter_x=MARGIN) + 0.30

    # C and D share a row: the candidate's evidence, and where it sits
    cw_ = 4.42
    dw = CW - cw_ - 0.30
    hc = place(fig, cand, MARGIN, y + 0.22, cw_, "C", letter_x=MARGIN)
    hd = place(fig, pose, MARGIN + cw_ + 0.30, y, dw, "D",
                letter_x=MARGIN + cw_ + 0.30)
    y += max(hc + 0.22, hd)

    written = []
    for suffix in (".pdf", ".png"):
        p = outputs.write_path("artist", TOPIC, "pin1_main_figure", suffix)
        fig.savefig(p, facecolor="white")
        written.append(p)
    plt.close(fig)
    print(f"  content ends at {y:.2f} in of {PAGE_H} (margin {MARGIN})")
    return written


if __name__ == "__main__":
    for p in build(Path(sys.argv[1])):
        print("WROTE", p)
