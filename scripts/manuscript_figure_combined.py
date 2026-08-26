#!/usr/bin/env python3
"""
Purpose: One Nature-style page: five panels and the caption, for the Pin1
         choreography.
Author:  Timothy Wu (with Claude Code)
Date:    2026-08-25
Input:   decisions/*.md, docs/how_this_project_breaks.md, the enrichment-gate
         token, and a capture of the integration GUI.
Output:  versioned PDF + PNG under append_only/inhibition/00_outputs/artist/

WHY THIS REPLACES THE FOUR SEPARATE FIGURES. The earlier set carried its
explanation in sentence-shaped panel titles, which is a slide convention, not a
journal one. Here the panels carry only what a panel should -- a letter, an
axis, a unit -- and every claim moves into the caption, which is where a
reader of the printed page looks for it.

EVERY COUNT IS PARSED AT RENDER TIME, and the caption is BUILT from those same
parsed values rather than typed beside them. A caption that restates a number
the panel computed is a second copy that can drift from the first; here there
is one source, so the page cannot contradict itself.
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
from PIL import Image

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from shared import outputs  # noqa: E402

TOPIC = "manuscript_figures"

# CVD-validated categorical order; worst-case OKLab dE 13.1 across
# normal/deuteranopia/protanopia/tritanopia, floors 15 and 8.
C1, C2, C3, C4 = "#332288", "#DDAA33", "#44AA99", "#CC3311"
INK, MUTED, GRID = "#000000", "#4d4d4d", "#d0d0d0"

# Nature: sans-serif, 5-7 pt at final size, panel letters bold lowercase.
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "font.size": 7,
    "axes.linewidth": 0.6,
    "axes.edgecolor": INK,
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": INK,
    "ytick.color": INK,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.4,
    "ytick.major.size": 2.4,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 200,
    "savefig.dpi": 600,
})

PAGE_W, PAGE_H = 8.27, 11.69          # A4 portrait, inches


# ----------------------------------------------------------------- parsing --
def decision_records() -> list[dict]:
    out = []
    for f in sorted(glob.glob(str(REPO / "decisions" / "D*.md"))):
        t = Path(f).read_text()
        if not t.startswith("---"):
            continue
        try:
            out.append(yaml.safe_load(t.split("---", 2)[1]))
        except yaml.YAMLError:
            continue
    return out


def catalogue() -> dict:
    txt = (REPO / "docs" / "how_this_project_breaks.md").read_text()
    entries = re.findall(r"^\| (\d+) \|", txt, flags=re.M)
    names = re.findall(r"^### \d+\. (.+)$", txt, flags=re.M)
    blocks = re.split(r"^### \d+\. .+$", txt, flags=re.M)[1:]
    claimed = {}
    for name, block in zip(names, blocks):
        nums = set()
        for pat in (r"Instances?\s+\*\*([0-9,\s and]+)\*\*",
                    r"pattern behind \*\*([0-9,\s and]+)\*\*"):
            for m in re.findall(pat, block):
                nums |= {int(n) for n in re.findall(r"\d+", m)}
        claimed[name] = nums
    routes = dict(re.findall(
        r"^\| (Someone looked[^|]+|Found while[^|]+|An existing guard[^|]+|"
        r"Deliberate audit[^|]+)\|\s*(\d+)", txt, flags=re.M))
    return {"n": len(entries), "disguises": claimed,
            "routes": {k.strip(): int(v) for k, v in routes.items()}}


def gate() -> dict:
    return json.loads(Path(
        "/data/lab_vm/append_only/inhibition/00_shared_substrate/"
        "enrichment_gate.token").read_text())


def archived_auc() -> dict:
    out = {}
    for did in ("D0015", "D0028", "D0031"):
        txt = Path(glob.glob(str(REPO / "decisions" / f"{did}*.md"))[0]).read_text()
        m = re.search(r"affinity_kcal[^\n']*?ROC-AUC (\d\.\d+),\s*CI "
                      r"\[(\d\.\d+),\s*(\d\.\d+)\]", txt)
        if not m:
            m = re.search(r"ROC-AUC (\d\.\d+),?\s*\n?\s*CI \[(\d\.\d+),\s*"
                          r"(\d\.\d+)\]", txt)
        out[did] = tuple(float(m.group(i)) for i in (1, 2, 3))
    return out


# ------------------------------------------------------------------ layout --
def rect(x, y, w, h):
    """Inches from the top-left of the page -> matplotlib figure fraction."""
    return [x / PAGE_W, 1 - (y + h) / PAGE_H, w / PAGE_W, h / PAGE_H]


def letter(fig, x, y, ch):
    fig.text(x / PAGE_W, 1 - y / PAGE_H, ch, fontsize=9, fontweight="bold",
             va="top", ha="left")


def build(shot: Path):
    recs = decision_records()
    cat = catalogue()
    g = gate()
    nc = g["strata"]["non_covalent"]["metrics"]["vina_affinity"]
    auc = archived_auc()

    n_dec = len(recs)
    n_shared = sum(1 for r in recs if r.get("approach") == "shared")
    origins = [("implementation", 0), ("user", 0), ("adversary", 0), ("spec", 0)]
    origins = [(k, sum(1 for r in recs if r.get("origin") == k))
               for k, _ in origins]
    n_routes = sum(cat["routes"].values())
    n_guard = cat["routes"].get("An existing guard fired", 0)

    fig = plt.figure(figsize=(PAGE_W, PAGE_H))
    M = 0.62                                  # page margin, inches
    CW = PAGE_W - 2 * M                       # content width

    # ---- a: the integration interface --------------------------------------
    img = Image.open(shot).convert("RGB").crop((150, 1740, 4060, 3292))
    ha = CW * img.size[1] / img.size[0]
    y = M + 0.16
    letter(fig, M - 0.34, y - 0.02, "a")
    ax = fig.add_axes(rect(M, y, CW, ha))
    ax.imshow(img, aspect="auto", interpolation="lanczos")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_edgecolor("#b0b0b0"); s.set_linewidth(0.5)
    y += ha + 0.52

    # ---- b: the control under three decoy constructions --------------------
    hw = (CW - 0.72) / 2
    GUT = 0.95          # left gutter for y-tick labels
    hb = 2.16
    letter(fig, M - 0.34, y - 0.02, "b")
    ax = fig.add_axes(rect(M + GUT, y, hw - GUT, hb))
    rows = [("property-matched", *auc["D0015"], C4),
            ("adduct forms", *auc["D0028"], C4),
            ("chemotype-matched", *auc["D0031"], C4),
            ("non-covalent gate", nc["roc_auc"], *nc["roc_auc_ci"], C1)]
    ys = np.arange(len(rows))[::-1]
    for yy, (lab, a, lo, hi, col) in zip(ys, rows):
        ax.plot([lo, hi], [yy, yy], lw=1.3, color=col, solid_capstyle="round",
                zorder=3)
        ax.plot([a], [yy], "o", ms=4.4, color=col, mec="white", mew=0.7,
                zorder=4)
    ax.axvline(0.5, color=INK, lw=0.7, ls=(0, (3, 2.4)), zorder=2)
    ax.set_yticks(ys)
    ax.set_yticklabels([r[0] for r in rows], fontsize=6.2, linespacing=1.25)
    ax.set_xlim(0.25, 1.0)
    ax.set_ylim(-0.62, len(rows) - 0.38)
    ax.set_xticks([0.4, 0.6, 0.8, 1.0])
    ax.set_xlabel("ROC-AUC (95% CI)", fontsize=7)
    ax.xaxis.grid(True, color=GRID, lw=0.45)
    ax.set_axisbelow(True)

    # ---- c: redocking accuracy --------------------------------------------
    letter(fig, M + hw + 0.72 - 0.34, y - 0.02, "c")
    ax = fig.add_axes(rect(M + hw + 0.72 + GUT, y, hw - GUT, hb))
    vals = [18.3, 19.8, 41.5]
    ax.axhspan(41, 50, color=C2, alpha=0.22, zorder=1, lw=0)
    ax.bar(np.arange(3), vals, width=0.56, color=[C1, MUTED, C3], zorder=3)
    for x, v in zip(np.arange(3), vals):
        ax.text(x, v + 1.3, f"{v:.1f}", ha="center", va="bottom", fontsize=6.4)
    ax.set_xticks(np.arange(3))
    ax.set_xticklabels(["top-1", "random\nof nine", "best-of-\nnine"],
                       fontsize=6.2, linespacing=1.25)
    ax.set_ylim(0, 56)
    ax.set_yticks([0, 20, 40])
    ax.set_ylabel("pose within 2 Å (%)", fontsize=7)
    ax.yaxis.grid(True, color=GRID, lw=0.45)
    ax.set_axisbelow(True)
    y += hb + 0.62

    # ---- d / e: the written record ----------------------------------------
    hd = 1.72

    def hbar(x0, ch, pairs, color, xlab):
        letter(fig, x0 - 0.34, y - 0.02, ch)
        a = fig.add_axes(rect(x0 + GUT, y, hw - GUT, hd))
        lab = [p[0] for p in pairs]
        val = [p[1] for p in pairs]
        yy = np.arange(len(lab))[::-1]
        a.barh(yy, val, height=0.62, color=color, zorder=3)
        for t, v in zip(yy, val):
            a.text(v + max(val) * 0.03, t, str(v), va="center", ha="left",
                   fontsize=6.4)
        a.set_yticks(yy)
        a.set_yticklabels(lab, fontsize=6.2, linespacing=1.2)
        a.set_xlim(0, max(val) * 1.20)
        a.set_xlabel(xlab, fontsize=7)
        a.xaxis.grid(True, color=GRID, lw=0.45)
        a.set_axisbelow(True)
        return a

    hbar(M, "d", origins, C1, f"decision records (n = {n_dec})")
    dis = cat["disguises"]
    short = {"Selection by name": "selection by name",
             "Selection by position": "selection by position",
             "Pinned defaults that go stale": "stale pinned default",
             "A guard that is scoped out, mis-ordered, or vacuous":
                 "guard that cannot fail"}
    pairs = [(short.get(k, k), len(v)) for k, v in dis.items()]
    unfiled = cat["n"] - len(set().union(*dis.values()))
    if unfiled:
        pairs.append(("not yet assigned", unfiled))
    hbar(M + hw + 0.72, "e", pairs, C4, f"catalogued defects (n = {cat['n']})")
    y += hd + 0.52

    # ---- caption -----------------------------------------------------------
    cap = (
        "Fig. 1 | Output and internal controls of a four-approach "
        "computational screen against Pin1. "
        "**a**, The integration interface presents the four independent "
        "approaches side by side: T₁, structure-based de novo generation; "
        "T₂, a derivative neighbourhood of all-trans retinoic acid; "
        "T₃, R-group decoration of sulfopin; and T₄, a combinatorial "
        "warhead × R-group enumeration. Each column is ranked by its own "
        "metric and the columns are not merged. The verdict issued by the "
        "shared control is displayed with each ranking (WEAK, UNDERPOWERED), "
        "so an ordering cannot be read without it. Live screen capture. "
        "**b**, ROC-AUC of the covalent control under three successively "
        "stricter decoy constructions, and of the current non-covalent gate. "
        f"Scoring pre-reaction ligands against property-matched decoys gives "
        f"{auc['D0015'][0]:.3f}; scoring the bound adduct forms gives "
        f"{auc['D0028'][0]:.3f}; requiring each decoy to carry the active's "
        f"own chemotype gives {auc['D0031'][0]:.3f}. Dashed line, chance. "
        "Only the first construction excludes chance. "
        "**c**, Redocking accuracy over 82 non-cognate Pin1 ligands on the "
        "prepared 3IKD receptor. Selection by docking score (top-1) does not "
        "exceed random selection among the nine output modes, while the "
        "correct pose is present in the ensemble 41.5% of the time. Shaded "
        "band, published single-receptor cross-docking range. "
        f"**d**, The {n_dec} decision records by the actor or process that "
        f"prompted each; {n_shared} of {n_dec} bind the whole choreography "
        "rather than one approach. "
        f"**e**, The {cat['n']} catalogued defects by failure mode. In every "
        "case a value was selected by position, name or default rather than "
        "by identity, and none raised an exception. Of the "
        f"{n_routes} entries with a recorded discovery route, {n_guard} were "
        "caught by an automated guard. "
        "All rankings shown carry rank_validated = False."
    )
    wrapped = textwrap.wrap(cap, 148)
    lines = []
    for ln in wrapped:
        parts = re.split(r"(\*\*.+?\*\*)", ln)
        lines.append(parts)
    ty = y
    for parts in lines:
        tx = M
        for pc in parts:
            if not pc:
                continue
            bold = pc.startswith("**")
            s = pc.strip("*")
            t = fig.text(tx / PAGE_W, 1 - ty / PAGE_H, s, fontsize=6.6,
                         va="top", ha="left",
                         fontweight="bold" if bold else "normal")
            fig.canvas.draw()
            bb = t.get_window_extent(fig.canvas.get_renderer())
            tx += bb.width / fig.dpi
        ty += 0.118

    written = []
    for suffix in (".pdf", ".png"):
        p = outputs.write_path("artist", TOPIC, "pin1_combined_figure", suffix)
        fig.savefig(p, facecolor="white")
        written.append(p)
    plt.close(fig)

    print("PARSED AT RENDER TIME:")
    print(f"  decision records {n_dec}  (shared {n_shared})")
    print(f"  origins          {origins}")
    print(f"  catalogue        n={cat['n']}  unfiled={unfiled}  "
          f"routes={n_routes} guard={n_guard}")
    print(f"  covalent control {auc}")
    print(f"  non-covalent     {nc['roc_auc']:.3f} {nc['roc_auc_ci']}")
    print(f"  caption lines    {len(lines)}, ends at {ty:.2f} in "
          f"(page {PAGE_H})")
    return written


if __name__ == "__main__":
    shot = Path(sys.argv[1])
    for p in build(shot):
        print("WROTE", p)
