#!/usr/bin/env python3
"""
Purpose: The two submission figures for the Pin1 choreography.
Author:  Timothy Wu (with Claude Code)
Date:    2026-08-25
Input:   GUI captures, a stored 3IKD pose render, decisions/*.md,
         docs/how_this_project_breaks.md, the enrichment-gate token.
Output:  versioned PDF + PNG under append_only/inhibition/00_outputs/artist/

TWO FIGURES, SPLIT BY THEME. Figure 1 is what the choreography presents to a
human; Figure 2 is what its own controls say about the ordering that
presentation carries. Forcing both onto one page made every panel too small to
read, and the two answer different questions.

PANEL TITLES ARE DELIBERATELY ABSENT. A panel carries a letter, an axis and a
unit; the caption carries the claim. The captions are BUILT from the same
parsed values the panels are drawn from, so the page cannot contradict itself.

THE POSE IN FIG 1c IS DRAWN AGAINST 3IKD, NOT 6VAJ. `pose3d.receptor_for()`
still defaults to `6VAJ_prepared.pdb` and the frame carries a `pose_path` from
that production run, so the easy render is the superseded one. D0059 replaced
the receptor on 2026-08-05 and every 6VAJ measurement is invalidated, so this
uses `nac3_pose_path` and the receptor registered against it. Both files exist
and both render a plausible picture, which is the whole problem.
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
from matplotlib.patches import Rectangle
from PIL import Image

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "integration" / "app"))
from shared import outputs  # noqa: E402

TOPIC = "manuscript_figures"
C1, C2, C3, C4 = "#332288", "#DDAA33", "#44AA99", "#CC3311"
INK, MUTED, GRID = "#000000", "#4d4d4d", "#d0d0d0"

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
    "font.size": 7, "axes.linewidth": 0.6, "axes.edgecolor": INK,
    "axes.labelcolor": INK, "text.color": INK, "xtick.color": INK,
    "ytick.color": INK, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.4, "ytick.major.size": 2.4,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 200, "savefig.dpi": 600,
})

PAGE_W = 8.27
MARGIN = 0.62
CW = PAGE_W - 2 * MARGIN
GUT = 0.95

#: The compound in Fig 1b and 1c -- ONE molecule across both panels, so the
#: reader sees the same candidate as a chemist draws it, as it was docked, and
#: as it sits in the site. T4's top-ranked candidate.
POSE_ID = "t4_c24c106bd005"


# ----------------------------------------------------------------- parsing --
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


def pose_facts():
    """What the compound in Fig 1b/1c is, and what its pose does -- measured.

    Delegates to the interaction renderer rather than recomputing contacts, so
    the numbers in the caption and the geometry in the panel cannot disagree.
    """
    sys.path.insert(0, str(REPO / "scripts"))
    import data as D
    import pose3d as p3d
    from pose_interaction_render import analyse, CONTACT_A, POLAR_A

    a = analyse(POSE_ID)
    row = a["row"]
    f, _ = D.load_frame("t4")
    same = f[f["warhead_class"] == row["warhead_class"]].dropna(
        subset=["affinity_kcal"])
    aff = float(row["affinity_kcal"])
    # TWO RANKS, AND THEY DISAGREE. `rank_synth` is what the shortlist shows
    # and it is 1; by RAW affinity_kcal this molecule is 10th of 187 in its
    # class, because the shortlist orders on the size-decorrelated residual of
    # that score (D0049), not on the score. Nine chloroacetamides score better
    # and are absent -- all of them pass the alert and docked-species gates, so
    # the ordering, not a filter, is what puts this one first. Reporting the
    # raw-score rank as "the rank" would be this project's own defect: two
    # populated, plausible rank columns and the wrong one named.
    rank_raw = int((same["affinity_kcal"] < aff).sum()) + 1
    rank_shortlist = row.get("rank_synth", row.get("rank"))
    polar = [c for c in a["contacts"] if c["polar"]]
    return {
        "aff": aff, "rank_raw": rank_raw,
        "rank_shortlist": int(rank_shortlist) if rank_shortlist == rank_shortlist
        else None, "n_class": len(same),
        "receptor": a["receptor"].name, "d_sg": a["d_sg"],
        "n_contacts": len(a["contacts"]), "n_polar": len(polar),
        "polar": [(f"{c['resn'].title()}{c['num']}", c["polar"][0])
                  for c in polar],
        "subpockets": [(sp.label, sp.colour) for sp in p3d.SUBPOCKETS],
        "warhead": row["warhead_class"], "adduct": row["adduct_class"],
        "contact_A": CONTACT_A, "polar_A": POLAR_A,
        "covalent_max": p3d.COVALENT_BOND_MAX_A,
    }


# ------------------------------------------------------------------ layout --
#: Text in the captions is positioned by MEASURING each run and advancing by
#: its width. get_window_extent reports at the figure's dpi, so measuring at
#: figure.dpi 200 and rasterising at savefig.dpi 600 made every run render
#: wider than it had been measured -- each bold panel letter was drawn on top
#: of the preceding full stop. Measure and rasterise at ONE dpi.
RENDER_DPI = 600


def mk(page_h):
    return plt.figure(figsize=(PAGE_W, page_h), dpi=RENDER_DPI)


def rect(page_h, x, y, w, h):
    return [x / PAGE_W, 1 - (y + h) / page_h, w / PAGE_W, h / page_h]


def letter(fig, page_h, x, y, ch):
    fig.text(x / PAGE_W, 1 - y / page_h, ch, fontsize=9, fontweight="bold",
             va="top", ha="left")


def _space_w(fig, size):
    """Width of one space, measured.

    matplotlib's text extent ignores leading and trailing whitespace, so a
    caption assembled from bold and plain runs loses exactly the space before
    every bold panel letter -- "presents**a**, ROC-AUC". Measured as the
    difference between two strings that differ only by a space.
    """
    r = fig.canvas.get_renderer()
    w = []
    for probe in ("i i", "ii"):
        t = fig.text(0, 0, probe, fontsize=size)
        fig.canvas.draw()
        w.append(t.get_window_extent(r).width)
        t.remove()
    return (w[0] - w[1]) / fig.dpi


def caption(fig, page_h, y, text, width=148, size=6.6):
    sw = _space_w(fig, size)
    r = fig.canvas.get_renderer()
    for ln in textwrap.wrap(text, width):
        tx = MARGIN
        for pc in re.split(r"(\*\*.+?\*\*)", ln):
            if not pc:
                continue
            bold = pc.startswith("**")
            body = pc.strip("*")
            lead = len(body) - len(body.lstrip())
            trail = len(body) - len(body.rstrip())
            body = body.strip()
            tx += lead * sw
            if body:
                t = fig.text(tx / PAGE_W, 1 - y / page_h, body, fontsize=size,
                             va="top", ha="left",
                             fontweight="bold" if bold else "normal")
                fig.canvas.draw()
                tx += t.get_window_extent(r).width / fig.dpi
            tx += trail * sw
        y += 0.118
    return y


def save(fig, stem):
    out = []
    for suffix in (".pdf", ".png"):
        p = outputs.write_path("artist", TOPIC, stem, suffix)
        fig.savefig(p, facecolor="white", dpi=RENDER_DPI)
        out.append(p)
    plt.close(fig)
    return out


# ------------------------------------------------------------------ fig 1 --
def figure1(shots: Path):
    pf = pose_facts()

    a = Image.open(shots / "f2_shortlists_wide.png").convert("RGB") \
             .crop((150, 1740, 4060, 3292))
    b = Image.open(shots / "f_dossier_t4.png").convert("RGB") \
             .crop((770, 752, 4060, 1768))
    c = Image.open(shots / "pose_interaction.png").convert("RGB")

    ha = CW * a.size[1] / a.size[0]
    hb = CW * b.size[1] / b.size[0]
    cw = 4.05
    hc = cw * c.size[1] / c.size[0]

    page_h = MARGIN + 0.16 + ha + 0.46 + hb + 0.50 + hc + 0.46 + 1.55 + MARGIN
    fig = mk(page_h)
    y = MARGIN + 0.16

    for img, h, ch, w in ((a, ha, "a", CW), (b, hb, "b", CW),
                          (c, hc, "c", cw)):
        letter(fig, page_h, MARGIN - 0.34, y - 0.02, ch)
        ax = fig.add_axes(rect(page_h, MARGIN, y, w, h))
        ax.imshow(img, aspect="auto", interpolation="lanczos")
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_edgecolor("#b0b0b0"); s.set_linewidth(0.5)
        if ch == "c":
            gx = MARGIN + cw + 0.28
            gy = y + 0.08
            fig.text(gx / PAGE_W, 1 - gy / page_h, "sub-pockets",
                     fontsize=6.8, fontweight="bold", va="top")
            gy += 0.19
            for label, colour in pf["subpockets"]:
                sw = fig.add_axes(rect(page_h, gx, gy - 0.052, 0.105, 0.070))
                sw.add_patch(Rectangle((0, 0), 1, 1, facecolor=colour,
                                       edgecolor="#7a7a7a", lw=0.5))
                sw.set_xlim(0, 1); sw.set_ylim(0, 1)
                sw.set_xticks([]); sw.set_yticks([])
                for sp_ in sw.spines.values():
                    sp_.set_visible(False)
                # Advance by the number of lines the label ACTUALLY wrapped
                # to, not by a length threshold guessed against a different
                # wrap width -- "proline-binding pocket" is 22 characters, wrapped
                # to two lines at width 21, and a `len < 24` test advanced one
                # line, so its second line landed on the next entry.
                wrapped = textwrap.wrap(label, 21)
                fig.text((gx + 0.155) / PAGE_W, 1 - gy / page_h,
                         "\n".join(wrapped), fontsize=6.3,
                         va="top", linespacing=1.3)
                gy += 0.115 + 0.105 * len(wrapped)
            sw = fig.add_axes(rect(page_h, gx, gy - 0.052, 0.105, 0.070))
            sw.add_patch(Rectangle((0, 0), 1, 1, facecolor="white",
                                   edgecolor="#7a7a7a", lw=0.5))
            sw.set_xlim(0, 1); sw.set_ylim(0, 1)
            sw.set_xticks([]); sw.set_yticks([])
            for sp_ in sw.spines.values():
                sp_.set_visible(False)
            fig.text((gx + 0.155) / PAGE_W, 1 - gy / page_h,
                     "claimed by none", fontsize=6.3, va="top")
            gy += 0.30
            fig.text(gx / PAGE_W, 1 - gy / page_h,
                     "— — polar contact", fontsize=6.3, va="top",
                     color="#c000c0")
            gy += 0.30
            polar_txt = ", ".join(f"{n} {d:.2f}" for n, d in pf["polar"])
            facts = (f"{POSE_ID}\n"
                     f"{pf['warhead'].replace('_', ' ')} \u00b7 rank "
                     f"{pf['rank_shortlist']} in class, T\u2084 shortlist\n"
                     f"docking score {pf['aff']:.2f} kcal mol\u207b\u00b9 "
                     f"({pf['rank_raw']}th of {pf['n_class']} by raw score \u2014\n"
                     f"  the shortlist orders on the size-decorrelated\n"
                     f"  residual of it, D0049)\n"
                     f"{pf['n_contacts']} residues within "
                     f"{pf['contact_A']} \u00c5\n"
                     f"{pf['n_polar']} polar contacts \u2264 "
                     f"{pf['polar_A']} \u00c5:\n  {polar_txt}\n"
                     f"warhead C to Cys113 SG, {pf['d_sg']:.2f} \u00c5\n"
                     f"receptor {pf['receptor']}")
            fig.text(gx / PAGE_W, 1 - gy / page_h, facts, fontsize=6.1,
                     va="top", color=MUTED, linespacing=1.5)
        y += h + (0.46 if ch == "a" else 0.50 if ch == "b" else 0.46)

    cap = (
        "Fig. 1 | The integration interface and what it presents. "
        "**a**, The four approaches side by side: T\u2081, structure-based de "
        "novo generation; T\u2082, a derivative neighbourhood of all-trans "
        "retinoic acid; T\u2083, R-group decoration of sulfopin; T\u2084, a "
        "combinatorial warhead \u00d7 R-group enumeration. Each column is "
        "ranked by its own metric and the columns are not merged. The verdict "
        "issued by the shared control is shown with each ranking (WEAK, "
        "UNDERPOWERED). "
        "**b**, A candidate is presented in both the form a chemist "
        "synthesises and the form that was docked. For a covalent candidate "
        "these differ: the leaving group is lost on reaction, so the docked "
        "species is the adduct. Three of T\u2084's nine warhead classes \u2014 "
        "chloroacetamide, sulfamate acetamide and sulfonate acetamide \u2014 "
        "converge on the same adduct and therefore dock as one species, so the "
        "enumeration explores nine classes but seven distinct docked "
        "molecules. "
        f"**c**, The same candidate's docked pose in the Pin1 site. It is "
        f"rank {pf['rank_shortlist']} of its warhead class on the T\u2084 "
        f"shortlist; by raw docking score it is {pf['rank_raw']}th of "
        f"{pf['n_class']} chloroacetamides, because the shortlist orders on "
        f"the size-decorrelated residual of that score (D0049) rather than the "
        f"score itself. The protein is drawn as a solvent-excluded surface "
        f"coloured by sub-pocket and the ligand in cyan sticks, so the panel "
        f"shows the shape the ligand occupies rather than burying it in side "
        f"chains. Every residue with a heavy atom within {pf['contact_A']} "
        f"\u00c5 is labelled, and the {pf['n_polar']} polar contacts \u2264 "
        f"{pf['polar_A']} \u00c5 are dashed. The pose contacts "
        f"{pf['n_contacts']} residues and spans all three sub-pockets. Its "
        f"nearest atom sits {pf['d_sg']:.2f} \u00c5 from "
        f"the Cys113 SG \u2014 outside the {pf['covalent_max']} \u00c5 bond "
        f"window, so this is a near-attack pose and not a formed adduct. "
        "Residue labels are shown for residues a sub-pocket claims or that "
        "make a polar contact. Drawn against the prepared 3IKD receptor "
        "adopted in D0059; poses from the superseded 6VAJ run are also present "
        "in the frames and render equally plausibly. Screen captures of the "
        "live interface, except c, which is the interface's own renderer run "
        "headless."
    )
    caption(fig, page_h, y, cap)
    return save(fig, "pin1_fig1_interface"), pf


# ------------------------------------------------------------------ fig 2 --
def figure2():
    recs, cat, g, auc = decisions(), catalogue(), gate(), archived()
    nc = g["strata"]["non_covalent"]["metrics"]["vina_affinity"]
    n_dec = len(recs)
    n_shared = sum(1 for r in recs if r.get("approach") == "shared")
    origins = [(k, sum(1 for r in recs if r.get("origin") == k))
               for k in ("implementation", "user", "adversary", "spec")]
    n_routes = sum(cat["routes"].values())
    n_guard = cat["routes"].get("An existing guard fired", 0)

    hw = (CW - 0.72) / 2
    hb, hd = 2.16, 1.72
    page_h = MARGIN + 0.16 + hb + 0.62 + hd + 0.52 + 1.35 + MARGIN
    fig = mk(page_h)
    y = MARGIN + 0.16

    # a: the control under three decoy constructions
    letter(fig, page_h, MARGIN - 0.34, y - 0.02, "a")
    ax = fig.add_axes(rect(page_h, MARGIN + GUT, y, hw - GUT, hb))
    rows = [("property-matched", *auc["D0015"], C4),
            ("adduct forms", *auc["D0028"], C4),
            ("chemotype-matched", *auc["D0031"], C4),
            ("non-covalent gate", nc["roc_auc"], *nc["roc_auc_ci"], C1)]
    ys = np.arange(len(rows))[::-1]
    for yy, (lab, v, lo, hi, col) in zip(ys, rows):
        ax.plot([lo, hi], [yy, yy], lw=1.3, color=col,
                solid_capstyle="round", zorder=3)
        ax.plot([v], [yy], "o", ms=4.4, color=col, mec="white", mew=0.7,
                zorder=4)
    ax.axvline(0.5, color=INK, lw=0.7, ls=(0, (3, 2.4)), zorder=2)
    ax.set_yticks(ys); ax.set_yticklabels([r[0] for r in rows], fontsize=6.2)
    ax.set_xlim(0.25, 1.0); ax.set_ylim(-0.62, len(rows) - 0.38)
    ax.set_xticks([0.4, 0.6, 0.8, 1.0])
    ax.set_xlabel("ROC-AUC (95% CI)")
    ax.xaxis.grid(True, color=GRID, lw=0.45); ax.set_axisbelow(True)

    # b: redocking accuracy
    letter(fig, page_h, MARGIN + hw + 0.72 - 0.34, y - 0.02, "b")
    ax = fig.add_axes(rect(page_h, MARGIN + hw + 0.72 + GUT, y, hw - GUT, hb))
    vals = [18.3, 19.8, 41.5]
    ax.axhspan(41, 50, color=C2, alpha=0.22, zorder=1, lw=0)
    ax.bar(np.arange(3), vals, width=0.56, color=[C1, MUTED, C3], zorder=3)
    for x, v in zip(np.arange(3), vals):
        ax.text(x, v + 1.3, f"{v:.1f}", ha="center", va="bottom", fontsize=6.4)
    ax.set_xticks(np.arange(3))
    ax.set_xticklabels(["top-1", "random\nof nine", "best-of-\nnine"],
                       fontsize=6.2, linespacing=1.25)
    ax.set_ylim(0, 56); ax.set_yticks([0, 20, 40])
    ax.set_ylabel("pose within 2 Å (%)")
    ax.yaxis.grid(True, color=GRID, lw=0.45); ax.set_axisbelow(True)
    y += hb + 0.62

    def hbar(x0, ch, pairs, colour, xlab):
        letter(fig, page_h, x0 - 0.34, y - 0.02, ch)
        a = fig.add_axes(rect(page_h, x0 + GUT, y, hw - GUT, hd))
        lab = [p[0] for p in pairs]; val = [p[1] for p in pairs]
        yy = np.arange(len(lab))[::-1]
        a.barh(yy, val, height=0.62, color=colour, zorder=3)
        for t, v in zip(yy, val):
            a.text(v + max(val) * 0.03, t, str(v), va="center", ha="left",
                   fontsize=6.4)
        a.set_yticks(yy); a.set_yticklabels(lab, fontsize=6.2)
        a.set_xlim(0, max(val) * 1.20); a.set_xlabel(xlab)
        a.xaxis.grid(True, color=GRID, lw=0.45); a.set_axisbelow(True)

    hbar(MARGIN, "c", origins, C1, f"decision records (n = {n_dec})")
    short = {"Selection by name": "selection by name",
             "Selection by position": "selection by position",
             "Pinned defaults that go stale": "stale pinned default",
             "A guard that is scoped out, mis-ordered, or vacuous":
                 "guard that cannot fail"}
    pairs = [(short.get(k, k), len(v)) for k, v in cat["dis"].items()]
    unfiled = cat["n"] - len(set().union(*cat["dis"].values()))
    if unfiled:
        pairs.append(("not yet assigned", unfiled))
    hbar(MARGIN + hw + 0.72, "d", pairs, C4,
         f"catalogued defects (n = {cat['n']})")
    y += hd + 0.52

    cap = (
        "Fig. 2 | Internal controls on the ordering the interface presents. "
        "**a**, ROC-AUC of the covalent control under three successively "
        "stricter decoy constructions, and of the current non-covalent gate. "
        f"Scoring pre-reaction ligands against property-matched decoys gives "
        f"{auc['D0015'][0]:.3f}; scoring the bound adduct forms gives "
        f"{auc['D0028'][0]:.3f}; requiring each decoy to carry the active's "
        f"own chemotype gives {auc['D0031'][0]:.3f}. Dashed line, chance. "
        "Only the first construction excludes chance, and each tightening was "
        "made in response to an audit rather than after the fact. "
        "**b**, Redocking accuracy over 82 non-cognate Pin1 ligands on the "
        "prepared 3IKD receptor. Selection by docking score does not exceed "
        "random selection among the nine output modes, while a pose within "
        "2 Å is present in the ensemble 41.5% of the time: the search finds "
        "the pose and the score cannot identify it. Shaded band, published "
        "single-receptor cross-docking range. "
        f"**c**, The {n_dec} decision records by the actor or process that "
        f"prompted each; {n_shared} of {n_dec} bind the whole choreography "
        "rather than one approach. "
        f"**d**, The {cat['n']} catalogued defects by failure mode. In every "
        "case a value was selected by position, name or default rather than "
        "by identity, and none raised an exception. Of the "
        f"{n_routes} entries with a recorded discovery route, {n_guard} were "
        "caught by an automated guard. Every ranking in Fig. 1 carries "
        "rank_validated = False."
    )
    caption(fig, page_h, y, cap)
    return save(fig, "pin1_fig2_controls"), (n_dec, n_shared, cat, auc, nc)


if __name__ == "__main__":
    shots = Path(sys.argv[1])
    p1, pf = figure1(shots)
    p2, meta = figure2()
    n_dec, n_shared, cat, auc, nc = meta
    print("PARSED AT RENDER TIME:")
    print(f"  pose compound  {POSE_ID} {pf['warhead']} aff {pf['aff']:.2f} | "
          f"shortlist rank {pf['rank_shortlist']} | raw-score rank "
          f"{pf['rank_raw']}/{pf['n_class']}")
    print(f"  pose receptor  {pf['receptor']}  warhead-SG {pf['d_sg']:.2f} A "
          f"(bond window <= {pf['covalent_max']})")
    print(f"  contacts       {pf['n_contacts']} within {pf['contact_A']} A, "
          f"{pf['n_polar']} polar: {pf['polar']}")
    print(f"  decisions      {n_dec} (shared {n_shared})")
    print(f"  catalogue      n={cat['n']} routes={sum(cat['routes'].values())}")
    print(f"  control        {auc}")
    print(f"  non-covalent   {nc['roc_auc']:.3f}")
    for p in (*p1, *p2):
        print("WROTE", p)
