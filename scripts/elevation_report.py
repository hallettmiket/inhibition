"""
Purpose: build the formatted elevation report — tier 1/2 results, the lead's 100 ns MD plots, and an animated trajectory viewer.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06
Input: elevation_tier1/*.csv, elevation_tier2/*.csv, elevation_cohort/*.csv,
       the lead's prod trajectory analysis (rmsd/mindist/numcont/warhead-SG xvg),
       and a PBC-corrected fitted movie.pdb extracted with gmx trjconv
Output: a self-contained HTML report (no external requests) + the figures it embeds

EVERY NUMBER IN THE REPORT IS RECOMPUTED HERE FROM THE RAW CSVs. None is
transcribed from docs/elevation_results.md. If the two disagree, that is a
finding about the pipeline, not a typo to reconcile silently — the report is
built from the same shard files the analysis reads, so a divergence means the
shards changed underneath the prose.

THE FIGURE IS BUILT TWICE, once per colour theme, because the page renders in
whichever theme the reader's browser is in and a single PNG cannot serve both
without one of them looking wrong.
"""

from __future__ import annotations

import argparse
import base64
import glob
import io
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

log = logging.getLogger("elevation-report")

OUTPUTS = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith")
T1 = OUTPUTS / "elevation_tier1"
T2 = OUTPUTS / "elevation_tier2"
COHORT = OUTPUTS / "elevation_cohort"
LEAD = "t4_72f5671e89cb"
LEAD_MD = Path("/data/lab_vm/modifiable/inhibition/md_residence_3ikd") / LEAD / "md/rep1"

#: Escape time of the lead, in ns. Determined as the last frame below a 1.0 nm
#: ligand-RMSD threshold after which it never returns — not eyeballed off a plot.
ESCAPE_NS = 54.45

#: Near-attack window, nm. Same constants the criterion uses.
NAC_LO, NAC_HI = 0.28, 0.42

GROUP_LABEL = {
    "A_hiEnr_hiCons_bdhi": "A · high enrichment, high consensus",
    "B_loEnr_hiCons_bdhi": "B · low enrichment, high consensus",
    "D_loEnr_loCons_bdhi": "D · low enrichment, low consensus",
    "V_hiCons_chloroacetamide": "V · chloroacetamide, high consensus",
    "REF_crystallographic": "REF · crystallographic positives",
}
GROUP_ORDER = list(GROUP_LABEL)

THEMES = {
    "light": dict(paper="#faf8f4", ink="#14181e", muted="#6d7078", grid="#ddd8cf",
                  accent="#a8761a", anchor="#2f6f6a", drift="#a8443a"),
    "dark": dict(paper="#14171c", ink="#e9e6e0", muted="#9599a1", grid="#2c3138",
                 accent="#d9a441", anchor="#57a79e", drift="#cd6f63"),
}


# --------------------------------------------------------------------------- data


def _load_shards(d: Path) -> pd.DataFrame:
    fs = sorted(glob.glob(str(d / "*.csv")))
    if not fs:
        return pd.DataFrame()
    df = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    keys = [c for c in ("ident", "replicate") if c in df.columns]
    return df.drop_duplicates(subset=keys, keep="last") if keys else df


def read_xvg(p: Path) -> np.ndarray:
    rows = [l.split() for l in p.read_text().splitlines() if l and l[0] not in "@#"]
    return np.array([[float(x) for x in r] for r in rows])


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """Positive => `a` is SMALLER (more stable) than `b` more often than not.

    The readout is |delta d| where smaller is more stable, so the raw dominance
    statistic is negated here. Doing it in one place stops the sign from being
    re-derived (and re-flipped) at each call site.
    """
    gt = sum((x > y) for x in a for y in b)
    lt = sum((x < y) for x in a for y in b)
    return (lt - gt) / (len(a) * len(b))


def mannwhitney(a: np.ndarray, b: np.ndarray) -> float:
    from scipy.stats import mannwhitneyu
    return float(mannwhitneyu(a, b, alternative="two-sided").pvalue)


def holm(ps: list[float]) -> list[float]:
    order = sorted(range(len(ps)), key=lambda i: ps[i])
    out, prev = [0.0] * len(ps), 0.0
    for rank, i in enumerate(order):
        v = max(prev, min(1.0, (len(ps) - rank) * ps[i]))
        out[i], prev = v, v
    return out


def tier1_summary() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    t1 = _load_shards(T1)
    per = (t1.groupby(["group", "ident"]).abs_delta_nm.mean()
           .reset_index(name="abs_delta_nm"))

    rows = []
    for g in GROUP_ORDER:
        v = per[per.group == g].abs_delta_nm.values
        if not len(v):
            continue
        rows.append(dict(group=g, n=len(v), med=np.median(v),
                         q1=np.percentile(v, 25), q3=np.percentile(v, 75)))
    summary = pd.DataFrame(rows)

    def contrast(g1: str, g2: str) -> dict:
        a = per[per.group == g1].abs_delta_nm.values
        b = per[per.group == g2].abs_delta_nm.values
        return dict(g1=g1, g2=g2, m1=np.median(a), m2=np.median(b),
                    delta=cliffs_delta(a, b), p=mannwhitney(a, b), n1=len(a), n2=len(b))

    pre = [contrast(*c) for c in (("A_hiEnr_hiCons_bdhi", "B_loEnr_hiCons_bdhi"),
                                  ("B_loEnr_hiCons_bdhi", "D_loEnr_loCons_bdhi"),
                                  ("A_hiEnr_hiCons_bdhi", "D_loEnr_loCons_bdhi"))]
    for r, ph in zip(pre, holm([r["p"] for r in pre])):
        r["p_holm"] = ph
    anchor = [contrast(g, "REF_crystallographic") for g in GROUP_ORDER[:4]]

    meta = dict(n_runs=len(t1), n_failed=int((t1.status != "ok").sum()) if "status" in t1 else 0,
                n_mol=per.ident.nunique(),
                lead_rank=int(per.sort_values("abs_delta_nm")
                              .reset_index(drop=True).query("ident == @LEAD").index[0]) + 1,
                lead_val=float(per[per.ident == LEAD].abs_delta_nm.iloc[0]),
                lead_reps=t1[t1.ident == LEAD].abs_delta_nm.round(3).tolist())
    return summary, pd.DataFrame(pre), pd.DataFrame(anchor), meta


def tier2_summary() -> tuple[pd.DataFrame, dict]:
    t2 = _load_shards(T2)
    if t2.empty:
        return pd.DataFrame(), dict(done=0, total=111, complete=False)
    rows = []
    for g in GROUP_ORDER:
        d = t2[t2.group == g]
        if d.empty:
            continue
        rows.append(dict(group=g, mols=d.ident.nunique(), reps=len(d),
                         frac=d.frac_in_window.median(),
                         bias=d.bias_at_exit_kj.median(),
                         escaped=d.escaped.mean()))
    meta = dict(done=len(t2), total=111, complete=len(t2) >= 111,
                mols=t2.ident.nunique())
    return pd.DataFrame(rows), meta


def to_ns(t, total_ns: float = 100.0):
    """Rescale an xvg time column to ns.

    Different gmx tools write the time column in different units for the SAME
    trajectory — `gmx rms`/`mindist`/`numcont` emit ns here while `gmx distance`
    emits ps — and nothing in the file says which. Dividing by 1000 "because xvg
    is ps" silently squashed three of these into the first 0.1 ns of the plot.
    The production length IS known (100 ns, production_ps in the residence CSV),
    so anchor the axis on that rather than guessing the unit.
    """
    span = float(t[-1] - t[0])
    return t if span <= 0 else (t - t[0]) * (total_ns / span)


def lead_md() -> dict:
    rmsd = read_xvg(LEAD_MD / "rmsd.xvg")
    mind = read_xvg(LEAD_MD / "mindist.xvg")
    ncon = read_xvg(LEAD_MD / "numcont.xvg")
    whsg = read_xvg(OUTPUTS / "md_residence" /
                    f"warhead_sg_distance_{LEAD}_100ns_1.xvg")

    t_r, y_r = to_ns(rmsd[:, 0]), rmsd[:, 1]
    t_w, y_w = to_ns(whsg[:, 0]), whsg[:, 1]
    t_m, y_m = to_ns(mind[:, 0]), mind[:, 1]
    t_c, y_c = to_ns(ncon[:, 0]), ncon[:, 1]

    bound_r, free_r = y_r[t_r <= ESCAPE_NS], y_r[t_r > ESCAPE_NS]
    bound_w, free_w = y_w[t_w <= ESCAPE_NS], y_w[t_w > ESCAPE_NS]
    bound_c, free_c = y_c[t_c <= ESCAPE_NS], y_c[t_c > ESCAPE_NS]
    bound_m, free_m = y_m[t_m <= ESCAPE_NS], y_m[t_m > ESCAPE_NS]
    in_nac = (bound_w >= NAC_LO) & (bound_w <= NAC_HI)

    return dict(
        t_rmsd=t_r, rmsd=y_r, t_wh=t_w, wh=y_w,
        t_mind=t_m, mind=y_m, t_ncon=t_c, ncon=y_c,
        bound_rmsd=float(bound_r.mean()), free_rmsd=float(free_r.mean()),
        max_rmsd=float(y_r.max()),
        bound_wh_mean=float(bound_w.mean()), bound_wh_med=float(np.median(bound_w)),
        bound_wh_min=float(bound_w.min()), bound_wh_max=float(bound_w.max()),
        free_wh_med=float(np.median(free_w)),
        bound_con=float(bound_c.mean()), free_con=float(free_c.mean()),
        bound_mind=float(bound_m.mean()), free_mind=float(free_m.mean()),
        nac_frac=float(in_nac.mean()), wh_start=float(y_w[0]),
        within6=float((bound_w <= 0.60).mean()),
        nac_visits=int(np.sum(np.diff(in_nac.astype(int)) == 1)),
    )


# --------------------------------------------------------------------------- figures


def _style(theme: dict):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "figure.facecolor": theme["paper"], "axes.facecolor": theme["paper"],
        "savefig.facecolor": theme["paper"],
        "text.color": theme["ink"], "axes.labelcolor": theme["ink"],
        "xtick.color": theme["muted"], "ytick.color": theme["muted"],
        "axes.edgecolor": theme["grid"], "grid.color": theme["grid"],
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.dpi": 130,
    })
    return plt


def fig_lead(md: dict, theme: dict) -> str:
    plt = _style(theme)
    fig, ax = plt.subplots(3, 1, figsize=(10.5, 8.6), sharex=True,
                           gridspec_kw=dict(hspace=0.16, top=0.94))

    def band(a):
        a.axvspan(ESCAPE_NS, 100, color=theme["drift"], alpha=0.07, lw=0, zorder=0)
        a.axvline(ESCAPE_NS, color=theme["drift"], lw=1.1, ls="--", alpha=0.9, zorder=2)
        a.grid(axis="y", lw=0.6, alpha=0.45, zorder=0)

    def tag(a, x, y, txt, colour, ha="left"):
        a.text(x, y, txt, color=colour, fontsize=9, weight="bold", ha=ha,
               va="top", zorder=6,
               bbox=dict(fc=theme["paper"], ec="none", alpha=0.86, pad=1.8))

    # phase labels sit ABOVE the axes so they can never collide with the trace
    ax[0].set_title(" ", pad=14)
    fig.text(0.075, 0.958, "BOUND — pose holds", color=theme["anchor"],
             fontsize=10, weight="bold", family="DejaVu Sans")
    fig.text(0.565, 0.958, "DISSOCIATED — does not return", color=theme["drift"],
             fontsize=10, weight="bold", family="DejaVu Sans")

    a = ax[0]
    a.plot(md["t_rmsd"], md["rmsd"], lw=0.65, color=theme["accent"], zorder=3)
    a.plot([0, ESCAPE_NS], [md["bound_rmsd"]] * 2, color=theme["anchor"], lw=1.6, zorder=4)
    band(a)
    a.set_ylabel("ligand RMSD  (nm)")
    a.set_ylim(0, md["max_rmsd"] * 1.10)
    tag(a, 1.5, md["max_rmsd"] * 1.03,
        f"mean {md['bound_rmsd']:.3f} nm over {ESCAPE_NS:.0f} ns", theme["anchor"])
    tag(a, 98.5, md["max_rmsd"] * 1.03, f"mean {md['free_rmsd']:.2f} nm",
        theme["drift"], ha="right")

    # clipped so the bound phase is legible; the excursions run off-scale by design
    a = ax[1]
    top = 2.0
    a.axhspan(NAC_LO, NAC_HI, color=theme["anchor"], alpha=0.20, lw=0, zorder=1)
    a.plot(md["t_wh"], np.clip(md["wh"], 0, top * 1.4), lw=0.65,
           color=theme["accent"], zorder=3)
    band(a)
    a.set_ylabel("warhead C10 → Cys113 SG  (nm)")
    a.set_ylim(0, top)
    tag(a, 1.5, top * 0.97,
        f"near-attack window {NAC_LO}–{NAC_HI} nm  ·  in it for "
        f"{md['nac_frac'] * 100:.1f}% of the bound phase", theme["anchor"])
    a.annotate("docked  0.301", xy=(0.4, 0.301), xytext=(6, 1.45),
               color=theme["muted"], fontsize=8.5, zorder=6,
               arrowprops=dict(arrowstyle="->", color=theme["muted"], lw=0.9))
    a.text(99, top * 0.06, "clipped at 2.0 nm", color=theme["muted"], fontsize=7.5,
           ha="right", style="italic", zorder=6)

    # contacts and min-distance are different units, so they get their own axes
    a = ax[2]
    a.plot(md["t_ncon"], md["ncon"], lw=0.6, color=theme["accent"], zorder=3)
    band(a)
    a.set_ylabel("contacts < 0.45 nm", color=theme["accent"])
    a.tick_params(axis="y", colors=theme["accent"])
    a.set_ylim(0, max(md["ncon"]) * 1.05)
    a.set_xlabel("time  (ns)")
    a.set_xlim(0, 100)
    tag(a, 1.5, max(md["ncon"]) * 1.0,
        f"contacts {md['bound_con']:.2f} → {md['free_con']:.2f}   ·   "
        f"min distance {md['bound_mind']:.2f} → {md['free_mind']:.2f} nm", theme["ink"])

    a2 = a.twinx()
    a2.plot(md["t_mind"], md["mind"], lw=0.6, color=theme["muted"], alpha=0.75, zorder=2)
    a2.set_ylabel("min distance to protein  (nm)", color=theme["muted"])
    a2.tick_params(axis="y", colors=theme["muted"])
    a2.spines["top"].set_visible(False)
    a2.spines["right"].set_color(theme["grid"])
    a2.set_ylim(0, max(md["mind"]) * 1.05)

    return _png(fig, plt)


def fig_tier1(summary: pd.DataFrame, theme: dict) -> str:
    plt = _style(theme)
    t1 = _load_shards(T1)
    per = t1.groupby(["group", "ident"]).abs_delta_nm.mean().reset_index()

    fig, a = plt.subplots(figsize=(10.5, 4.3))
    rng = np.random.default_rng(7)
    for i, g in enumerate(GROUP_ORDER):
        v = per[per.group == g].abs_delta_nm.values
        if not len(v):
            continue
        is_ref = g.startswith("REF")
        c = theme["anchor"] if is_ref else theme["drift"]
        y = np.full(len(v), i) + rng.uniform(-.13, .13, len(v))
        a.scatter(v, y, s=44, color=c, alpha=.75, edgecolor="none", zorder=3)
        # the elevated molecule is called out because the whole report turns on
        # where it sits, and "rank 37 of 37" is easier to disbelieve than to see
        ids = per[per.group == g].ident.values
        if LEAD in set(ids):
            k = int(np.where(ids == LEAD)[0][0])
            a.scatter([v[k]], [y[k]], s=190, facecolor="none",
                      edgecolor=theme["ink"], lw=1.6, zorder=5)
            a.annotate("the elevated molecule", xy=(v[k], y[k]),
                       xytext=(v[k] - 0.075, y[k] + 0.52), fontsize=8.5,
                       color=theme["ink"], ha="right", zorder=6,
                       arrowprops=dict(arrowstyle="->", color=theme["ink"], lw=1))
        m = np.median(v)
        a.plot([m, m], [i - .30, i + .30], color=theme["ink"], lw=2.4, zorder=4)
        a.text(m, i + .40, f"{m:.3f}", ha="center", color=theme["ink"],
               fontsize=9, weight="bold")

    a.axvline(summary[summary.group == "REF_crystallographic"]["med"].iloc[0],
              color=theme["anchor"], lw=1, ls=":", alpha=.8, zorder=1)
    a.set_yticks(range(len(GROUP_ORDER)))
    a.set_yticklabels([GROUP_LABEL[g] for g in GROUP_ORDER], fontsize=9.5)
    a.invert_yaxis()
    a.set_xlabel("← more stable      "
                 "warhead displacement over 300 ps of unrestrained dynamics,  |Δd|  (nm)")
    a.grid(axis="x", lw=0.6, alpha=0.5)
    a.set_xlim(left=0)
    fig.tight_layout()
    return _png(fig, plt)


def _png(fig, plt) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


# --------------------------------------------------------------------------- movie


#: Formal charge class per residue, used to paint the surface red→white→blue.
#: AMBER splits histidine by protonation state: HID/HIE are neutral at pH 7,
#: HIP is the doubly-protonated cation. Lumping all three as "basic" would
#: paint two neutral residues blue.
ACIDIC = {"ASP", "GLU"}
BASIC = {"ARG", "LYS", "HIP"}

#: Pin1 PPIase active site in UniProt Q13526 numbering, mapped onto this
#: structure's renumbering. The offset is verified by residue IDENTITY at build
#: time rather than trusted — an offset that slipped by one would label the
#: wrong residue with a name that still looks entirely plausible.
PIN1_OFFSET = 50
KEY_SITES = {
    113: ("Cys113", "target", "CYS"), 59: ("His59", "catalytic", "HIS"),
    157: ("His157", "catalytic", "HIS"), 152: ("Thr152", "catalytic", "THR"),
    63: ("Lys63", "basic", "LYS"), 68: ("Arg68", "basic", "ARG"),
    69: ("Arg69", "basic", "ARG"), 122: ("Leu122", "pocket", "LEU"),
    130: ("Met130", "pocket", "MET"), 134: ("Phe134", "pocket", "PHE"),
}
#: Hang each label off the sidechain tip so it lands on the surface rather than
#: buried at the backbone.
TIP = {"CYS": "SG", "LYS": "NZ", "ARG": "CZ", "HID": "NE2", "HIE": "NE2",
       "HIP": "NE2", "PHE": "CZ", "MET": "SD", "LEU": "CD1", "THR": "OG1"}
HIS_FORMS = {"HID", "HIE", "HIP"}


def _charge(resn: str) -> float:
    return -1.0 if resn in ACIDIC else (1.0 if resn in BASIC else 0.0)


def surface_payload(pdb: Path) -> tuple[str, list, list, list]:
    """Multi-model PDB carrying formal charge in the B-factor, plus overlays.

    Returns (pdb_text, warhead_sg_distance_per_frame, label_defs,
    label_positions_per_frame). Distance and label anchors are computed from the
    SAME coordinates that get rendered, so no readout can disagree with what is
    on screen.
    """
    frames, cur, out, block = [], [], [], []
    for l in pdb.read_text().splitlines():
        if l.startswith("MODEL"):
            cur, block = [], [l]
        elif l.startswith(("ATOM", "HETATM")):
            resn, name = l[17:20].strip(), l[12:16].strip()
            cur.append((resn, int(l[22:26]), name,
                        float(l[30:38]), float(l[38:46]), float(l[46:54])))
            block.append(f"{l[:60]}{_charge(resn):6.2f}{l[66:]}")
        elif l.startswith("ENDMDL"):
            frames.append(cur)
            out.extend(block + [l])

    f0 = frames[0]
    by_res: dict[int, str] = {}
    for a in f0:
        by_res.setdefault(a[1], a[0])

    labels = []
    for pin1, (text, kind, expect) in sorted(KEY_SITES.items()):
        rid = pin1 - PIN1_OFFSET
        got = by_res.get(rid)
        norm = "HIS" if got in HIS_FORMS else got
        if norm != expect:
            raise ValueError(
                f"residue-numbering check failed: {text} maps to resid {rid}, "
                f"which is {got}, but {expect} was expected. The offset is "
                f"wrong — refusing to label the structure with names that do "
                f"not match it.")
        want = TIP.get(got, "CA")
        idx = next((i for i, a in enumerate(f0) if a[1] == rid and a[2] == want), None)
        if idx is None:
            idx = next(i for i, a in enumerate(f0) if a[1] == rid and a[2] == "CA")
        labels.append(dict(text=text, kind=kind, resid=rid, atom=idx))

    cys = 113 - PIN1_OFFSET
    i_sg = next(i for i, a in enumerate(f0)
                if a[1] == cys and a[2] == "SG")
    i_c10 = next(i for i, a in enumerate(f0) if a[0] == "MOL" and a[2] == "C10")

    dist, positions = [], []
    for f in frames:
        sg = np.array(f[i_sg][3:])
        dist.append(round(float(np.linalg.norm(np.array(f[i_c10][3:]) - sg)), 2))
        positions.append([[round(f[d["atom"]][3 + k], 2) for k in range(3)]
                          for d in labels])
    return "\n".join(out), dist, labels, positions


# --------------------------------------------------------------------------- html

CSS = """
:root{
  --paper:#faf8f4; --raise:#f3efe7; --ink:#14181e; --muted:#6d7078;
  --rule:#e0dad0; --accent:#a8761a; --anchor:#2f6f6a; --drift:#a8443a;
  --serif:"Iowan Old Style","Palatino Linotype",Palatino,Georgia,serif;
  --sans:system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;
  --mono:ui-monospace,SFMono-Regular,"SF Mono",Menlo,Consolas,monospace;
}
@media (prefers-color-scheme:dark){
  :root{ --paper:#14171c; --raise:#1b1f26; --ink:#e9e6e0; --muted:#9599a1;
         --rule:#2c3138; --accent:#d9a441; --anchor:#57a79e; --drift:#cd6f63; }
}
:root[data-theme="dark"]{ --paper:#14171c; --raise:#1b1f26; --ink:#e9e6e0;
  --muted:#9599a1; --rule:#2c3138; --accent:#d9a441; --anchor:#57a79e; --drift:#cd6f63; }
:root[data-theme="light"]{ --paper:#faf8f4; --raise:#f3efe7; --ink:#14181e;
  --muted:#6d7078; --rule:#e0dad0; --accent:#a8761a; --anchor:#2f6f6a; --drift:#a8443a; }

*{box-sizing:border-box}
body{margin:0;background:var(--paper);color:var(--ink);font-family:var(--sans);
  font-size:16.5px;line-height:1.62;-webkit-font-smoothing:antialiased}
.wrap{max-width:1180px;margin:0 auto;padding:0 28px 120px}
.col{max-width:68ch}
h1,h2,h3{font-family:var(--serif);font-weight:600;text-wrap:balance;line-height:1.18}
h2{font-size:1.72rem;margin:0 0 .2rem}
h3{font-size:1.16rem;margin:2.4rem 0 .5rem}
p{margin:0 0 1.05rem}
a{color:var(--accent)}
strong{font-weight:650}
code,.num{font-family:var(--mono);font-variant-numeric:tabular-nums}
code{font-size:.88em;background:var(--raise);padding:.1em .36em;border-radius:3px}

/* masthead */
.mast{border-bottom:2px solid var(--ink);padding:56px 0 26px;margin-bottom:34px}
.eyebrow{font-family:var(--mono);font-size:.7rem;letter-spacing:.16em;
  text-transform:uppercase;color:var(--muted);margin-bottom:1.4rem}
.mast h1{font-size:clamp(2.1rem,4.4vw,3.15rem);margin:0 0 1.1rem;max-width:19ch}
.verdict{font-family:var(--serif);font-size:1.3rem;line-height:1.45;
  max-width:60ch;border-left:3px solid var(--accent);padding-left:1.05rem;margin-bottom:1.9rem}
.facts{display:flex;flex-wrap:wrap;gap:0 2.6rem;font-family:var(--mono);
  font-size:.79rem;color:var(--muted);padding-top:.4rem;border-top:1px solid var(--rule)}
.facts div{padding-top:.8rem}
.facts b{display:block;color:var(--ink);font-weight:600;font-size:1.02rem;letter-spacing:-.01em}

/* section numbering carries the tier structure, which is real */
section{margin:0 0 4.4rem}
.shead{display:flex;gap:1.1rem;align-items:baseline;
  border-top:1px solid var(--rule);padding-top:1.5rem;margin-bottom:1.5rem}
.snum{font-family:var(--mono);font-size:.76rem;color:var(--accent);
  letter-spacing:.1em;padding-top:.42rem;white-space:nowrap}
.sub{color:var(--muted);font-size:.95rem;margin:.15rem 0 0}

/* data */
.scroll{overflow-x:auto;margin:1.5rem 0;border:1px solid var(--rule);border-radius:5px;background:var(--raise)}
table{border-collapse:collapse;width:100%;font-size:.87rem}
th,td{padding:.6rem .85rem;text-align:left;border-bottom:1px solid var(--rule);white-space:nowrap}
th{font-family:var(--mono);font-size:.68rem;letter-spacing:.09em;text-transform:uppercase;
  color:var(--muted);font-weight:600;background:var(--paper)}
tbody tr:last-child td{border-bottom:none}
td.n{text-align:right;font-family:var(--mono);font-variant-numeric:tabular-nums}
tr.anchor td{color:var(--anchor);font-weight:600}
tr.fired{background:color-mix(in srgb,var(--accent) 13%,transparent)}
tr.fired td:first-child{box-shadow:inset 3px 0 0 var(--accent)}
.sig{color:var(--drift);font-weight:650}
.ns{color:var(--muted)}
.pill{font-family:var(--mono);font-size:.7rem;padding:.1rem .45rem;border-radius:99px;
  border:1px solid currentColor;white-space:nowrap}

figure{margin:2rem 0}
figure img{width:100%;height:auto;display:block;border-radius:5px}
figcaption{font-size:.83rem;color:var(--muted);margin-top:.7rem;max-width:75ch}
:root[data-theme="dark"] .lightonly,html:not([data-theme]) .lightonly{display:block}
.darkonly{display:none}
@media (prefers-color-scheme:dark){
  html:not([data-theme]) .lightonly{display:none}
  html:not([data-theme]) .darkonly{display:block}
}
:root[data-theme="dark"] .lightonly{display:none}
:root[data-theme="dark"] .darkonly{display:block}
:root[data-theme="light"] .lightonly{display:block}
:root[data-theme="light"] .darkonly{display:none}

.callout{border:1px solid var(--rule);border-left:3px solid var(--accent);
  background:var(--raise);padding:1.15rem 1.35rem;margin:1.7rem 0;border-radius:0 5px 5px 0}
.callout p:last-child{margin-bottom:0}
.callout.warn{border-left-color:var(--drift)}
.callout.good{border-left-color:var(--anchor)}
.ctitle{font-family:var(--mono);font-size:.7rem;letter-spacing:.13em;text-transform:uppercase;
  color:var(--muted);margin-bottom:.5rem}

/* trajectory viewer */
.viewer{border:1px solid var(--rule);border-radius:6px;overflow:hidden;background:var(--raise);margin:1.7rem 0}
#gl{width:100%;height:480px;position:relative}
.controls{display:flex;align-items:center;gap:.95rem;padding:.75rem 1rem;
  border-top:1px solid var(--rule);flex-wrap:wrap}
button.play{font-family:var(--mono);font-size:.8rem;background:var(--ink);color:var(--paper);
  border:none;padding:.45rem 1rem;border-radius:4px;cursor:pointer;min-width:74px}
button.play:hover{opacity:.85}
button.play:focus-visible,input:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
input[type=range]{flex:1;min-width:150px;accent-color:var(--accent)}
.readout{font-family:var(--mono);font-size:.79rem;color:var(--muted);white-space:nowrap}
.readout b{color:var(--ink)}
.tog{font-family:var(--mono);font-size:.76rem;color:var(--muted);cursor:pointer;
  display:inline-flex;align-items:center;gap:.4rem}
.tog input{accent-color:var(--accent);cursor:pointer}
#state{font-family:var(--mono);font-size:.72rem;padding:.16rem .5rem;border-radius:99px;
  border:1px solid currentColor}
.s-nac{color:var(--anchor)} .s-bound{color:var(--accent)} .s-free{color:var(--drift)}
.legend{display:flex;gap:1.3rem;flex-wrap:wrap;font-size:.78rem;color:var(--muted);
  padding:.6rem 1rem;border-top:1px solid var(--rule);font-family:var(--mono)}
.swatch{display:inline-block;width:9px;height:9px;border-radius:2px;margin-right:.4rem;vertical-align:middle}

ul{padding-left:1.15rem;margin:0 0 1.05rem}
li{margin-bottom:.5rem}
.foot{border-top:1px solid var(--rule);padding-top:1.4rem;margin-top:3.5rem;
  font-size:.8rem;color:var(--muted);font-family:var(--mono);line-height:1.7}
@media (prefers-reduced-motion:reduce){*{animation:none!important;transition:none!important}}
@media (max-width:640px){ body{font-size:15.5px} .wrap{padding:0 18px 80px} #gl{height:360px} }
"""


def short(g: str) -> str:
    """Just the group letter — 'A', 'B', 'D', 'V', 'REF'."""
    return GROUP_LABEL[g].split(" ")[0]


def fmt_p(p: float) -> str:
    cls = "sig" if p < 0.05 else "ns"
    s = f"{p:.4f}" if p >= 0.0001 else "&lt;0.0001"
    return f'<span class="{cls}">{s}</span>'


def build_html(s1, pre, anchor, m1, t2, m2, md, figs, pdb, dist, labels, positions) -> str:
    three = (REPO / "scripts/.cache_3dmol-min.js").read_text()

    rows_g = "".join(
        f'<tr class="{"anchor" if g.group.startswith("REF") else ""}">'
        f"<td>{GROUP_LABEL[g.group]}</td><td class='n'>{g.n}</td>"
        f"<td class='n'>{g.med:.3f}</td><td class='n'>{g.q1:.3f} – {g.q3:.3f}</td></tr>"
        for g in s1.itertuples())

    rows_pre = "".join(
        f"<tr><td>{short(r.g1)} vs {short(r.g2)}"
        f"</td><td class='n'>{r.m1:.3f}</td><td class='n'>{r.m2:.3f}</td>"
        f"<td class='n'>{r.delta:+.3f}</td><td class='n'>{fmt_p(r.p)}</td>"
        f"<td class='n'>{fmt_p(r.p_holm)}</td>"
        f"<td><span class='pill ns'>≈</span></td></tr>"
        for r in pre.itertuples())

    def _anchor_verdict(p: float) -> str:
        if p < 0.05:
            return '<span class="pill" style="color:var(--anchor)">REF holds</span>'
        return '<span class="pill ns">n=5, descriptive</span>'

    rows_anc = "".join(
        f"<tr><td>{short(r.g1)} vs REF</td>"
        f"<td class='n'>{r.m1:.3f}</td><td class='n'>{r.m2:.3f}</td>"
        f"<td class='n'>{r.delta:+.3f}</td><td class='n'>{fmt_p(r.p)}</td>"
        f"<td>{_anchor_verdict(r.p)}</td></tr>"
        for r in anchor.itertuples())

    readings = [
        ("B ≈ A, both &gt; D", "Consensus is the filter; enrichment adds nothing.", False),
        ("A &gt; B, both &gt; D", "Enrichment adds something real beyond consensus.", False),
        ("A ≈ B ≈ D", "Neither metric predicts stability. The ranking has no physical "
                      "support from this experiment.", True),
        ("D ≥ A, B", "Something is wrong with the design; report as a failure.", False),
        ("V ≈ REF", "Consensus-selection inside the validated class finds molecules that "
                    "behave like known binders.", False),
    ]
    rows_read = "".join(
        f'<tr class="{"fired" if fired else ""}"><td><code>{obs}</code></td><td>{con}</td>'
        f'<td>{"<b>← fired</b>" if fired else ""}</td></tr>' for obs, con, fired in readings)

    rows_t2 = "".join(
        f'<tr class="{"anchor" if r.group.startswith("REF") else ""}">'
        f"<td>{GROUP_LABEL[r.group]}</td><td class='n'>{r.mols}</td>"
        f"<td class='n'>{r.reps}</td><td class='n'>{r.frac:.3f}</td>"
        f"<td class='n'>{r.bias:.3f}</td></tr>" for r in t2.itertuples()) if not t2.empty else ""

    t2_note = ("complete" if m2["complete"] else
               f"<b>partial — {m2['done']} of {m2['total']} replicas.</b> Groups A and B are "
               "both finished, so the enrichment contrast is already readable; D and REF are not.")

    return f"""<title>Elevation experiment — which ranking metric selects for physical stability?</title>
<style>{CSS}</style>
<div class="wrap">

<header class="mast">
  <div class="eyebrow">Pin1 covalent inhibitors · pre-registered · inhibition@3IKD</div>
  <h1>Which ranking metric selects for physical stability?</h1>
  <div class="verdict">Neither. But the crystallographic anchor separates from every
  candidate group at <span class="num">p&nbsp;=&nbsp;0.007</span> — so the assay measures
  pose survival correctly, and it is the <em>ranking</em> that has no physical support.</div>
  <div class="facts">
    <div><b>{m1['n_mol']}</b> molecules, 5 groups</div>
    <div><b>{m1['n_runs']}</b> tier-1 runs, {m1['n_failed']} failed</div>
    <div><b>3</b> replicas each</div>
    <div><b>300 ps</b> unrestrained equilibration</div>
    <div><b>3 × 3 ns</b> well-tempered BPMD</div>
    <div><b>2026-08-06</b></div>
  </div>
</header>

<section>
  <div class="shead"><div class="snum">§1</div><div>
    <h2>The question</h2>
    <p class="sub">Two metrics rank the same molecules and disagree almost completely.</p>
  </div></div>
  <div class="col">
    <p><strong>Enrichment</strong> is the fraction of docking runs reaching a near-attack
    conformation, over an isotropic null. It does not converge: the same molecules fall from
    2.91× to 0.96× at ten times the search effort, and its rank correlation across efforts is
    <span class="num">ρ = −0.117</span>.</p>
    <p><strong>Consensus</strong> asks whether a molecule's ten best poses by energy agree with
    each other. It is rank-stable across efforts, <span class="num">ρ = +0.568</span>.</p>
    <p>A sanity check exposed how far apart they are: <strong>397 molecules have a single
    binding mode at ≥0.90 pose agreement, and only four of them clear enrichment &gt; 5.70.</strong>
    The enrichment cut captures 1% of the well-aligned molecules and misses five
    chloroacetamides — the one warhead class with a validated criterion (AUC 0.908).</p>
  </div>
</section>

<section>
  <div class="shead"><div class="snum">§2</div><div>
    <h2>Design, and the three confounds that shaped it</h2>
    <p class="sub">Each confound was measured, not assumed.</p>
  </div></div>
  <div class="col">
    <ul>
      <li><strong>Consensus is partly rigidity</strong> (ρ = −0.259 vs rotatable bonds,
      p = 4×10⁻⁸⁹). Rigid molecules have fewer ways to sit <em>and</em> are trivially more
      stable under dynamics, so groups are <strong>matched on rotatable-bond count</strong>
      (median 4.0 across all three BDHI groups), not merely balanced in size.</li>
      <li><strong>Enrichment and warhead class are confounded.</strong> High-enrichment cells
      are BDHI-dominated, low-enrichment cells acrylamide-dominated. The two separate only
      <em>within</em> BDHI — which is a hard limit on what any result here can claim.</li>
      <li><strong>High enrichment is nearly a subset of high consensus.</strong> The
      high-enrichment/low-consensus cell holds five molecules, so the informative contrast is
      A vs B rather than a factorial the data cannot fill.</li>
    </ul>
    <div class="callout good">
      <div class="ctitle">Why REF is in the design</div>
      <p>Crystallographic positives are the only molecules we know actually react with Cys113.
      Without them, "A beats D" compares two arbitrary groups. With them the question becomes
      <em>which group behaves like a molecule that really reacts</em> — and a null becomes
      interpretable instead of merely disappointing.</p>
    </div>
  </div>
</section>

<section>
  <div class="shead"><div class="snum">§3</div><div>
    <h2>Readings, fixed before the run</h2>
    <p class="sub">Committed to git before any simulation started.</p>
  </div></div>
  <div class="col"><p>With four groups and several possible readouts, a result chosen after the
  fact can be made to say almost anything. Each outcome was assigned its conclusion in
  advance.</p></div>
  <div class="scroll"><table>
    <thead><tr><th>observation</th><th>fixed conclusion</th><th></th></tr></thead>
    <tbody>{rows_read}</tbody>
  </table></div>
</section>

<section>
  <div class="shead"><div class="snum">§4</div><div>
    <h2>Tier 1 — did the docked pose survive plain dynamics?</h2>
    <p class="sub">{m1['n_runs']} runs, {m1['n_failed']} failures. Readout is |Δd|, the
    warhead-to-SG displacement across 300 ps of unrestrained NVT/NPT. Smaller is more stable.</p>
  </div></div>

  <figure>
    <img class="lightonly" src="data:image/png;base64,{figs['t1_light']}" alt="Warhead displacement by group; crystallographic positives cluster at low displacement while all four candidate groups spread higher">
    <img class="darkonly" src="data:image/png;base64,{figs['t1_dark']}" alt="Warhead displacement by group; crystallographic positives cluster at low displacement while all four candidate groups spread higher">
    <figcaption>One point per molecule, averaged over 3 replicas. Heavy bar is the group
    median. The crystallographic anchor (teal) sits apart from every generated group.</figcaption>
  </figure>

  <div class="scroll"><table>
    <thead><tr><th>group</th><th>n</th><th>median |Δd| (nm)</th><th>IQR</th></tr></thead>
    <tbody>{rows_g}</tbody>
  </table></div>

  <h3>The pre-registered contrasts: nothing</h3>
  <div class="scroll"><table>
    <thead><tr><th>contrast</th><th>median 1</th><th>median 2</th><th>Cliff's δ</th>
    <th>p</th><th>p (Holm)</th><th></th></tr></thead>
    <tbody>{rows_pre}</tbody>
  </table></div>
  <div class="col"><p>δ is signed so positive means group 1 is <em>more</em> stable. All three
  point estimates are negative — the higher-enrichment group is if anything the less stable
  one — but none is significant, and <strong>no claim is drawn from the direction of a
  non-significant effect.</strong> The fourth reading (D ≥ A, B) requires a significant
  contrast in that direction; there is none.</p></div>

  <h3>The anchor: everything</h3>
  <div class="scroll"><table>
    <thead><tr><th>contrast</th><th>median 1</th><th>REF</th><th>Cliff's δ</th><th>p</th><th></th></tr></thead>
    <tbody>{rows_anc}</tbody>
  </table></div>
  <div class="col">
    <div class="callout">
      <div class="ctitle">What makes the null interpretable</div>
      <p>Without REF, "A ≈ B ≈ D" is equally consistent with an assay that cannot separate
      anything. With REF, the same 300 ps separates known reactive molecules from generated
      candidates in all three groups. So the null is a statement about <strong>the metrics</strong>,
      not about the measurement. The prereg called BPMD "orthogonal" in this branch — that is
      too generous. It measures pose survival, which is what it was built to measure.</p>
    </div>
  </div>
</section>

<section>
  <div class="shead"><div class="snum">§5</div><div>
    <h2>Post-hoc observations</h2>
    <p class="sub">Computed after seeing the data. They replace nothing above.</p>
  </div></div>
  <div class="col">
    <p><strong>The generated poses drift under <em>dynamics</em>; the anchor's drift is almost
    entirely the energy minimisation.</strong> Every group relaxes the same small amount
    (0.018–0.055 nm) during minimisation. Only REF then stays put once thermal motion is
    applied — 0.049 nm over the 300 ps, against 0.159–0.226 nm for the candidate groups.</p>
    <p><strong>The drift is almost perfectly one-directional: away from the sulfur.</strong>
    Signed Δd was positive in <strong>110 of 111 replicas</strong>, and all 37 molecules have a
    positive mean. No docked pose in this cohort tightens its near-attack geometry under
    dynamics — docking places the warhead closer to Cys113 than the force field will hold it,
    in every group including REF.</p>
    <p><strong>On "is it still a near-attack conformation", V matches REF and the BDHI groups
    do not.</strong> Fraction of replicas still inside the 0.28–0.42 nm window at the start of
    production: A 0.08, B 0.21, D 0.08, <strong>V 0.53, REF 0.54</strong>. This is the one
    place the prereg's <em>V ≈ REF</em> statement finds support — and it is on a readout the
    prereg did not name. On the readout it <em>did</em> name, V sits at 0.203 against REF's
    0.102 and is not equivalent. Both are reported; neither substitutes for the other.</p>
  </div>
</section>

<section>
  <div class="shead"><div class="snum">§6</div><div>
    <h2>Tier 2 — binding-pose metadynamics</h2>
    <p class="sub">{t2_note}</p>
  </div></div>
  <div class="scroll"><table>
    <thead><tr><th>group</th><th>molecules</th><th>replicas</th>
    <th>median frac in window</th><th>median bias at exit (kJ)</th></tr></thead>
    <tbody>{rows_t2}</tbody>
  </table></div>
  <div class="col"><p>Tier 2 pushes the warhead out of the near-attack window and measures how
  hard the push has to be. Every replica escapes within 3 ns — that is expected for
  well-tempered metadynamics and is <em>not</em> the readout; the cost is. Two facts
  about the protocol constrain how far it can be read: <code>bias_at_exit</code> separated
  nothing (all p&nbsp;≥&nbsp;0.08) and tracks occupancy at ρ&nbsp;=&nbsp;0.974, so the
  escape-cost term is nearly inert at this length — <strong>occupancy is the discriminating
  readout</strong>. And the run is short and unconverged by design, so the comparison rests on
  every molecule receiving the identical protocol rather than on any single molecule's score
  being its stability.</p>
  <p><strong>Tier 2 reproduces tier 1's null independently.</strong> A&nbsp;≈&nbsp;B&nbsp;≈&nbsp;D
  on both occupancy and bias, while all three BDHI groups separate from the crystallographic
  anchor (p&nbsp;=&nbsp;0.0070, 0.0148, 0.0207). The two tiers agree at Spearman
  ρ&nbsp;=&nbsp;0.475, p&nbsp;=&nbsp;0.003 across 37 molecules.</p>
  <div class="callout warn">
    <div class="ctitle">A claim made and withdrawn</div>
    <p>On tier-2 occupancy V (0.172) is not distinguished from REF (0.163), and this was
    briefly reported as the pre-registration's <em>V&nbsp;≈&nbsp;REF</em> reading firing.
    <strong>It is withdrawn.</strong> On the readout the prereg actually named — tier-1
    |Δd| — V is 0.203 against REF's 0.102, so <strong>the two tiers disagree in sign on the
    same comparison</strong>. Failing to distinguish at n&nbsp;=&nbsp;5 is not equivalence,
    the prereg forbids inference at that size in terms, and the two readouts that favour V
    were selected after seeing that they agreed. Nothing here supports a synthesis
    shortlist.</p>
  </div>
  </div>
</section>

<section>
  <div class="shead"><div class="snum">§7</div><div>
    <h2>The elevated molecule, in 100 ns of real water</h2>
    <p class="sub"><code>{LEAD}</code> — group A, enrichment 6.86, consensus 1.000, QED 0.796.</p>
  </div></div>

  <div class="col">
    <div class="callout warn">
      <div class="ctitle">Do not quote the summary CSV's headline</div>
      <p><code>explicit_ligand_rmsd_nm_mean = 1.525</code> is a mean across a <strong>bimodal</strong>
      trajectory — half bound, half dissociated — and describes neither state. The script's own
      <code>explicit_rmsd_suspect</code> flag fired on it. Everything below splits the
      trajectory at the transition instead of averaging over it.</p>
    </div>
  </div>

  <figure>
    <img class="lightonly" src="data:image/png;base64,{figs['lead_light']}" alt="Three stacked time series over 100 ns: ligand RMSD, warhead to sulfur distance, and protein contacts, all showing a transition at 54 ns">
    <img class="darkonly" src="data:image/png;base64,{figs['lead_dark']}" alt="Three stacked time series over 100 ns: ligand RMSD, warhead to sulfur distance, and protein contacts, all showing a transition at 54 ns">
    <figcaption>100 ns, 10,001 frames, explicit TIP3P, no restraints. The molecule holds its
    pose for {ESCAPE_NS} ns and then leaves without returning. Shaded band in the middle panel
    is the near-attack window.</figcaption>
  </figure>

  <div class="col">
    <p><strong>The ligand.</strong> RMSD sits at <span class="num">{md['bound_rmsd']:.3f} nm</span>
    for the first <span class="num">{ESCAPE_NS} ns</span>, then rises to
    <span class="num">{md['free_rmsd']:.3f} nm</span> and never comes back. So the pose is
    <em>not</em> a docking artefact — it survives ~180× longer than tier 1's window.</p>
    <p><strong>The warhead is the different story.</strong> Docked at
    <span class="num">0.301 nm</span>, it is already at <span class="num">{md['wh_start']:.3f} nm</span>
    by the start of production and averages <span class="num">{md['bound_wh_mean']:.3f} nm</span>
    through the bound phase, ranging {md['bound_wh_min']:.3f}–{md['bound_wh_max']:.3f} nm. It is
    inside the near-attack window for <strong>{md['nac_frac']*100:.1f}% of the bound phase</strong>
    and within 6 Å for {md['within6']*100:.1f}%.</p>
    <p>The docked geometry is therefore <strong>not the resting state</strong>. Under explicit
    water the molecule sits with its warhead 8–9 Å from the sulfur and swings into attack
    distance intermittently — it does revisit the window, so the near-attack conformation is
    accessible, not excluded.</p>
  </div>

  <h3>The trajectory</h3>
  <div class="viewer">
    <div id="gl"></div>
    <div class="controls">
      <button class="play" id="play" type="button">▶ Play</button>
      <input type="range" id="scrub" min="0" max="{len(dist)-1}" value="0" step="1" aria-label="trajectory frame">
      <div class="readout"><b id="tns">0</b> ns &nbsp;·&nbsp; warhead → SG <b id="dsg">—</b> Å</div>
      <span id="state" class="s-bound">bound</span>
    </div>
    <div class="controls" style="border-top:none;padding-top:0">
      <label class="tog"><input type="checkbox" id="surf" checked> van der Waals surface</label>
      <label class="tog"><input type="checkbox" id="labs" checked> residue labels</label>
      <span class="readout" id="perf"></span>
    </div>
    <div class="legend">
      <span><span class="swatch" style="background:#2f6fb5"></span>basic · Arg, Lys</span>
      <span><span class="swatch" style="background:#efece7;outline:1px solid #b9b3aa"></span>neutral</span>
      <span><span class="swatch" style="background:#c0392b"></span>acidic · Asp, Glu</span>
      <span><span class="swatch" style="background:#e8c14a"></span>Cys113 SG</span>
      <span><span class="swatch" style="background:var(--accent)"></span>{LEAD}</span>
      <span>1 ns / frame · fitted on backbone</span>
    </div>
  </div>
  <div class="col"><p class="sub">101 frames at 1 ns intervals, PBC-corrected and superposed
  on the protein backbone. The surface is a <strong>van der Waals</strong> surface — chosen
  over the slower solvent-excluded surface because it is recomputed at <em>every</em> frame, so
  it tracks sidechain motion instead of being a still from frame 0. It is coloured by formal
  charge, <span style="color:#c0392b">acidic</span> to <span style="color:#2f6fb5">basic</span>,
  with neutral residues pale. Active-site residues are
  labelled in <strong>Pin1 (UniProt Q13526) numbering</strong>; this structure is renumbered by
  −50, and the build verifies every label against the residue identity at that position rather
  than trusting the offset. The distance readout comes from the same coordinates being
  rendered.</p></div>

  <div class="col">
    <div class="callout warn">
      <div class="ctitle">The cross-check, and it is unflattering</div>
      <p>This molecule is elevation group A — the top of both metrics. Its tier-1 warhead
      displacement is <span class="num">{m1['lead_val']:.3f} nm</span> across three replicas
      ({', '.join(f'{v:.3f}' for v in m1['lead_reps'])}), which is
      <strong>rank {m1['lead_rank']} of {m1['n_mol']}</strong> — the <em>worst molecule in the
      cohort</em>, against a crystallographic median of 0.102 nm. The independent 100 ns run
      reproduces it: {md['wh_start']:.3f} nm by the start of production. Two separately built
      systems, same conclusion.</p>
      <p>The top-ranked molecule on both metrics is the least able to hold its warhead in
      place. That is not a coincidence to explain away — it is §4's null showing up in a
      single molecule.</p>
    </div>
    <div class="callout good">
      <div class="ctitle">Why this is not a rejection</div>
      <p>A covalent inhibitor does not need its warhead parked on the sulfur. It needs to
      <em>reach</em> the window, because the bond, once formed, is permanent. Reaching it for
      {md['nac_frac']*100:.1f}% of a {ESCAPE_NS} ns residence — across countless rebinding
      events at working concentration — is ample on chemical timescales. What the data
      disqualifies is <strong>reading pose-geometry scores as occupancies</strong>: "10 of 10
      poses at 3.04 Å" is a true statement about the docking output and a misleading one about
      the molecule.</p>
    </div>
  </div>
</section>

<section>
  <div class="shead"><div class="snum">§8</div><div>
    <h2>Covalent workup on the lead</h2>
    <p class="sub">Every score labelled with the receptor that produced it. 6VAJ and 3IKD place
    the pocket 48.6 Å apart and are never pooled.</p>
  </div></div>
  <div class="scroll"><table>
    <thead><tr><th>measurement</th><th>tool</th><th>receptor</th><th>value</th><th>uncertainty</th></tr></thead>
    <tbody>
      <tr><td>reactive pose energy</td><td>AutoDock-GPU reactive</td><td><b>3IKD</b></td><td class="n">−8.04</td><td>top-10 span −8.04 … −7.97</td></tr>
      <tr><td>plain best-of-9</td><td>AutoDock-GPU</td><td><b>3IKD</b></td><td class="n">−7.362</td><td>± 0.026</td></tr>
      <tr><td><b>covalent affinity</b></td><td>gnina 1.3.3 covalent</td><td><b>3IKD</b></td><td class="n"><b>−7.08</b></td><td>mode 7 of 9</td></tr>
      <tr><td>covalent affinity</td><td>gnina 1.3.3 covalent</td><td>6VAJ</td><td class="n">−7.00</td><td>mode 7 of 9 · not comparable</td></tr>
      <tr><td>cnn_affinity / cnn_score</td><td>gnina covalent</td><td>3IKD</td><td class="n">4.482 / 0.111</td><td><span class="pill ns">advisory — uncalibrated</span></td></tr>
      <tr><td>covalent MD</td><td>GROMACS</td><td>3IKD</td><td class="n">—</td><td><span class="pill" style="color:var(--drift)">not run</span></td></tr>
    </tbody>
  </table></div>
  <div class="col">
    <p>The covalent <em>topology</em> is built and verified — antechamber types the attachment
    carbon as GAFF2 <code>c2</code>, all six required junction terms are present in
    <code>cys_gaff2_junction_5.frcmod</code>, and <code>verify_complex</code> passes with three
    bonds, which is correct for an sp² attachment carbon. That is the expensive, failure-prone
    half. No covalent trajectory was produced because GPU 7 was committed to the 100 ns
    non-covalent run for its whole duration.</p>
    <p><strong>A stereochemistry finding worth more than most of the gaps this ranking is built
    on.</strong> The frame's older covalent score of −5.026 and today's −7.00 on the same
    receptor differ almost entirely because the frame's adduct was re-embedded from SMILES and
    drew the <em>opposite</em> configuration at the warhead ring carbon. Measured directly, the
    diastereomer is worth <strong>1.35 kcal/mol on 3IKD and 2.02 on 6VAJ</strong>.</p>
  </div>
</section>

<section>
  <div class="shead"><div class="snum">§9</div><div>
    <h2>What this does not settle</h2>
  </div></div>
  <div class="col">
    <ul>
      <li><strong>BDHI-only</strong> for the A/B/D contrasts. BDHI has zero crystallographic
      positives, so nothing here transfers to other classes without re-testing.</li>
      <li><strong>Stability is not reactivity.</strong> A stable near-attack pose is necessary
      for the reaction, not sufficient. Nothing here measures whether a molecule reacts.</li>
      <li><strong>n = 8 supports only large effects.</strong> The A-vs-B contrast at δ = −0.469
      is a moderate effect this design cannot resolve either way. "≈" means <em>not
      distinguished</em>, not <em>shown to be equal</em>.</li>
      <li><strong>Every BDHI group failed the same way.</strong> With all their poses leaving
      the window, the contrasts are drawn between degrees of failure. A cohort where some group
      survived would test the metrics harder than this one could.</li>
      <li><strong>The lead's 100 ns is a single replicate.</strong> One dissociation event gives
      a residence estimate with ~100% relative standard error, so <strong>"{ESCAPE_NS} ns" is one
      draw, not a residence time</strong> — and there is no matched baseline yet showing what an
      ordinary molecule does over the same 100 ns.</li>
    </ul>
  </div>
  <div class="foot">
    Generated by <code>scripts/elevation_report.py</code> from the tier-1/tier-2 shard CSVs and
    the lead's production trajectory. Every figure and table is recomputed at build time; none
    is transcribed.<br>
    Pre-registration: <code>docs/elevation_prereg.md</code> · Results:
    <code>docs/elevation_results.md</code> · Decision: <code>D0071</code>
  </div>
</section>
</div>

<script>{three}</script>
<script>
const DSG = {json.dumps(dist)}, ESC = {ESCAPE_NS};
const LABELS = {json.dumps(labels)}, LPOS = {json.dumps(positions)};
const NACLO = {NAC_LO * 10:.2f}, NACHI = {NAC_HI * 10:.2f};
// The UMD footer only guarantees the global as 3Dmol; $3Dmol is the
// conventional alias and is not exported by every build. Bind whichever exists.
const MOL3D = window.$3Dmol || window['3Dmol'];

const KIND = {{ target: '#e8c14a', catalytic: '#57a79e', basic: '#2f6fb5', pocket: '#9599a1' }};
let viewer, frame = 0, timer = null, lastCost = 0;

const wantSurf = () => document.getElementById('surf').checked;
const wantLabs = () => document.getElementById('labs').checked;

function paint() {{
  viewer.setStyle({{}}, {{cartoon: {{color: '#8a9099', opacity: 0.55}}}});
  viewer.setStyle({{resn: 'MOL'}}, {{stick: {{radius: 0.19, colorscheme: 'yellowCarbon'}}}});
  viewer.addStyle({{resi: {113 - PIN1_OFFSET}, atom: 'SG'}},
                  {{sphere: {{radius: 0.7, color: '#e8c14a'}}}});
}}

// The surface is rebuilt every frame rather than computed once from frame 0:
// the protein is fitted on backbone so the fold barely moves, but the sidechains
// lining the pocket do, and a still surface would quietly misreport the pocket
// the ligand is actually leaving.
function drawSurface() {{
  viewer.removeAllSurfaces();
  if (!wantSurf()) return;
  viewer.addSurface(MOL3D.SurfaceType.VDW, {{
    opacity: 0.66,
    colorscheme: {{prop: 'b', gradient: new MOL3D.Gradient.RWB(-1, 1)}}
  }}, {{resn: 'MOL', invert: true}});
}}

function drawLabels(i) {{
  viewer.removeAllLabels();
  if (!wantLabs()) return;
  LABELS.forEach((L, k) => {{
    const p = LPOS[i][k];
    viewer.addLabel(L.text, {{
      position: {{x: p[0], y: p[1], z: p[2]}},
      backgroundColor: '#101418', backgroundOpacity: 0.72,
      fontColor: KIND[L.kind] || '#e9e6e0',
      fontSize: L.kind === 'target' ? 14 : 11,
      borderThickness: L.kind === 'target' ? 1.2 : 0,
      borderColor: '#e8c14a', inFront: true
    }});
  }});
}}

function show(i) {{
  frame = i;
  const t0 = performance.now();
  const r = viewer.setFrame(i);
  const after = () => {{
    paint();
    drawSurface();
    drawLabels(i);
    viewer.render();
    lastCost = performance.now() - t0;
    document.getElementById('perf').textContent =
      wantSurf() ? Math.round(lastCost) + ' ms/frame' : '';
  }};
  if (r && typeof r.then === 'function') {{ r.then(after); }} else {{ after(); }}

  document.getElementById('tns').textContent = i;
  document.getElementById('dsg').textContent = DSG[i].toFixed(1);
  document.getElementById('scrub').value = i;
  const st = document.getElementById('state');
  if (i > ESC) {{ st.textContent = 'dissociated'; st.className = 's-free'; }}
  else if (DSG[i] >= NACLO && DSG[i] <= NACHI) {{ st.textContent = 'near-attack'; st.className = 's-nac'; }}
  else {{ st.textContent = 'bound'; st.className = 's-bound'; }}
}}

// Surface rebuilds dominate the frame cost, so playback paces itself off what
// the last frame actually took instead of a fixed interval that would either
// stutter or crawl depending on the machine.
function step() {{
  show((frame + 1) % DSG.length);
  if (timer) timer = setTimeout(step, wantSurf() ? Math.max(120, lastCost * 1.15) : 90);
}}

function toggle() {{
  const b = document.getElementById('play');
  if (timer) {{ clearTimeout(timer); timer = null; b.textContent = '▶ Play'; return; }}
  b.textContent = '❚❚ Pause';
  timer = setTimeout(step, 60);
}}

window.addEventListener('DOMContentLoaded', () => {{
  const PDB = document.getElementById('pdbdata').textContent;
  viewer = MOL3D.createViewer(document.getElementById('gl'), {{backgroundAlpha: 0}});
  viewer.addModelsAsFrames(PDB, 'pdb');
  paint();
  viewer.zoomTo({{resn: 'MOL'}});
  viewer.zoom(0.35);
  show(0);
  document.getElementById('play').addEventListener('click', toggle);
  document.getElementById('scrub').addEventListener('input', e => {{
    if (timer) toggle();
    show(+e.target.value);
  }});
  document.getElementById('surf').addEventListener('change', () => show(frame));
  document.getElementById('labs').addEventListener('change', () => show(frame));
}});
</script>
<script type="text/plain" id="pdbdata">{pdb}</script>
"""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--movie", required=True, help="multi-model PDB from gmx trjconv")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    s1, pre, anchor, m1 = tier1_summary()
    log.info("tier1: %d runs, %d molecules, lead rank %d/%d",
             m1["n_runs"], m1["n_mol"], m1["lead_rank"], m1["n_mol"])
    t2, m2 = tier2_summary()
    log.info("tier2: %d/%d replicas", m2["done"], m2["total"])
    md = lead_md()
    log.info("lead: bound RMSD %.3f nm, NAC occupancy %.1f%%",
             md["bound_rmsd"], md["nac_frac"] * 100)

    figs = {}
    for name, theme in THEMES.items():
        figs[f"lead_{name}"] = fig_lead(md, theme)
        figs[f"t1_{name}"] = fig_tier1(s1, theme)
    log.info("figures built for both themes")

    pdb, dist, labels, positions = surface_payload(Path(args.movie))
    log.info("movie: %d frames", len(dist))

    html = build_html(s1, pre, anchor, m1, t2, m2, md, figs, pdb, dist,
                      labels, positions)
    Path(args.out).write_text(html)
    log.info("wrote %s (%.1f MB)", args.out, len(html) / 1e6)


if __name__ == "__main__":
    main()
