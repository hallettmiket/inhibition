"""
Purpose: combine several per-molecule MD reports into one browsable page.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-08
Input: --candidates <ident...> (their reports under mdprio_reports/)
Output: 00_outputs/blacksmith/mdprio_reports/combined_<N>.html

@tt8804: *"combine them into one. I want to see the movies and plots."*

WHY NOT ONE CONCATENATED FILE. Each report is ~9.5 MB because the movie frames
and every plot are embedded as base64 -- self-contained by design, so a report can
be copied anywhere and still work. Four of them inlined into a single document is
~38 MB of HTML, which the browser must parse before it shows anything and which no
artefact host will accept.

So this builds a FRAME-BASED index instead: one page, a molecule picker across the
top with each molecule's headline numbers, and the selected report rendered whole
in an iframe beneath. Every movie and every plot is the real one from the original
report -- nothing is regenerated or downsampled -- and only the report being looked
at is loaded. The originals stay individually openable.

THE COMPARISON TABLE IS THE POINT. Flipping between four reports to remember which
molecule was 50% attack-ready is the work this is meant to remove, so the numbers
that decide the shortlist sit above the viewer where they can be read together:
the 10 ns sweep readings, the 100 ns engagement, and which binding MODE was
actually elevated -- that last one because a molecule promoted on its minority
mode is a different claim from one promoted on its dominant mode.
"""

from __future__ import annotations

import argparse
import glob
import html
import json
import logging
import os
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import outputs as sout                     # noqa: E402
from shared import pipeline_schematic as schematic     # noqa: E402
from shared import mode_ranking as moderank        # noqa: E402
from shared import mode_key                        # noqa: E402
from shared import run_paths as rp                  # noqa: E402

log = logging.getLogger("mdprio-combine")
#: Set from config in main(); see D0085.
_HELD_BAR = 0.35
B = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith")
REPORTS = rp.reports_dir()


def _sweep() -> pd.DataFrame:
    fs = sorted(glob.glob(str(rp.sweep_dir() / "attack_sweep_*.csv")),
                key=lambda p: int(p.rsplit("_", 1)[1].split(".")[0]))
    if not fs:
        return pd.DataFrame()
    d = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    d = d[(d.get("sweep_ps", 0) > 1000) & (d.status == "ok")]
    # the mode that was ELEVATED is the best-scoring surviving mode, which is how
    # the worker chose it -- not necessarily mode 0
    return d.sort_values("frac_attack_ready", ascending=False) \
            .drop_duplicates("parent_ident")


def _thumbs(idents) -> dict:
    """A small 2D depiction per molecule, for the selector.

    Base64 rather than inline SVG markup: RDKit emits an XML declaration and an
    HTML comment, and both break parsing when they land inside markup the browser
    is already mid-way through. Encoding removes the question entirely.
    """
    import base64
    import re as _re
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Draw, AllChem
    RDLogger.DisableLog("rdApp.*")
    smi = {}
    for d in ("04_t4_combinatorial/D4", "03_t3_reinvent/D3"):
        sub, stem = d.split("/")
        fs = sorted(glob.glob(f"/data/lab_vm/append_only/inhibition/{sub}/{stem}_*.parquet"),
                    key=lambda q: int(q.rsplit("_", 1)[1].split(".")[0]))
        if not fs:
            continue
        fr = pd.read_parquet(fs[-1]).drop_duplicates("candidate_id")
        smi.update(dict(zip(fr.candidate_id, fr.canonical_smiles)))
    # Controls are not rows in D3/D4 -- they come from crystal structures. Their
    # SMILES is in the pose sidecar written beside the pose. Without this they
    # rendered as a text badge while every candidate showed a structure, which
    # made the control look like a different sort of object rather than the same
    # sort of object with a different provenance.
    import json as _json
    for p in sorted(glob.glob(str(B / "pose_sidecars/*.json"))):
        try:
            v = _json.loads(Path(p).read_text()).get("canonical_smiles")
        except Exception:                                  # noqa: BLE001
            continue
        if isinstance(v, str):
            smi.setdefault(Path(p).stem, v)

    # And the reference set itself, for controls with no sidecar -- ref_ATRA had
    # none, so it sat in the rail as a blank tile while every candidate beside it
    # showed a structure. A control that looks like a different KIND of object is
    # harder to compare against, which is the whole reason it is on the rail.
    # Keyed on `ref_<Name>` and on `ref_<Name>__<warhead>`, since the screen
    # writes the mechanism into the ident.
    rs = sorted(glob.glob(str(REPO / "data/reference/pin1_reference_binders_*.csv")),
                key=lambda q: int(q.rsplit("_", 1)[1].split(".")[0]))
    if rs:
        try:
            rd = pd.read_csv(rs[-1])
            name_col = rd.columns[0]
            for _, r in rd.iterrows():
                got = [v for v in r.values
                       if isinstance(v, str) and len(v) > 8 and Chem.MolFromSmiles(v)]
                if not got:
                    continue
                key = f"ref_{str(r[name_col]).strip()}"
                smi.setdefault(key, got[0])
                for i in idents:
                    if str(i).startswith(key + "__"):
                        smi.setdefault(i, got[0])
        except Exception as exc:                           # noqa: BLE001
            log.warning("reference depictions unavailable: %s", exc)

    out = {}
    for i in idents:
        v = smi.get(i)
        m = Chem.MolFromSmiles(v) if isinstance(v, str) else None
        if m is None:
            continue
        AllChem.Compute2DCoords(m)
        d2 = Draw.rdMolDraw2D.MolDraw2DSVG(96, 64)
        d2.drawOptions().bondLineWidth = 1
        Draw.rdMolDraw2D.PrepareAndDrawMolecule(d2, m)
        d2.FinishDrawing()
        svg = _re.sub(r"<\?xml.*?\?>", "", d2.GetDrawingText(), flags=_re.S)
        svg = _re.sub(r"<!--.*?-->", "", svg, flags=_re.S)
        out[i] = "data:image/svg+xml;base64," + base64.b64encode(svg.encode()).decode()
    return out


def _qualifying_controls() -> set[str]:
    """Reference molecules that earn a place on the rail as controls.

    @tt8804: *"for controls we just need confirmed pin1 inhibitors with the
    warheads used"*. Two conditions, both read from the reference table rather
    than hand-listed here:

    * ``mechanism == covalent_cys113`` — a confirmed covalent inhibitor of the
      residue this project targets. That drops ATRA, EGCG, PiB and every
      non-covalent peptide/phosphonate: they are Pin1 binders, but they cannot be
      compared against a near-attack criterion because they have nothing to
      attack with.
    * a ``warhead_class`` naming chemistry the screen actually enumerates. A
      control carrying a warhead we never make cannot tell us whether our screen
      would have found our own chemistry.

    Matched on the class VOCABULARY, not on the prose string: the table's
    warhead_class is written for a human ("1;4-naphthoquinone (Michael
    acceptor)"), and only a couple of rows would match a class id verbatim.
    """
    ours = {"chloroacetamide", "sulfamate", "sulfonate", "bdhi", "naphthoquinone",
            "acrylamide", "cinnamamide", "snar", "chloropyrimidine", "isoxazole"}
    fs = sorted(glob.glob(str(REPO / "data/reference/pin1_reference_binders_*.csv")),
                key=lambda q: int(q.rsplit("_", 1)[1].split(".")[0]))
    if not fs:
        return set()
    keep = set()
    try:
        d = pd.read_csv(fs[-1])
        for _, r in d.iterrows():
            if str(r.get("mechanism", "")).strip() != "covalent_cys113":
                continue
            w = str(r.get("warhead_class", "")).lower()
            if w in ("", "nan", "unverified") or not any(k in w for k in ours):
                continue
            keep.add(f"ref_{str(r['name']).strip()}")
    except Exception as exc:                               # noqa: BLE001
        log.warning("could not read the reference table: %s", exc)
    return keep


def _classes() -> dict:
    """Warhead class per molecule, for the within-class ranking view."""
    out = {}
    for tier, score in (("T4", "conditional_eb"), ("T3", "enrichment_conditional")):
        fs = sorted(glob.glob(str(B / f"rank_v2/rank_v2_{tier}_{score}_*.csv")))
        if not fs:
            continue
        d = pd.read_csv(fs[-1]).drop_duplicates("parent_ident")
        out.update(dict(zip(d.parent_ident, d.warhead_class)))
    return out


def _md() -> pd.DataFrame:
    """The 100 ns rows, preferring a run that SUCCEEDED over one that did not.

    `drop_duplicates(keep="last")` over an unsorted glob decided this by file
    order. A molecule can hold several 100 ns rows -- a launch that died in setup
    writes one too -- so for rx_6VAJ the GUI showed either a real result or a
    blank depending on which filename `glob` happened to return last. The failed
    row is not a later measurement; it is not a measurement at all.

    Rows are ordered by whether they carry a result, and only then deduped.
    """
    fs = glob.glob(str(rp.residence_dir() / "*.csv"))
    if not fs:
        return pd.DataFrame()
    d = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    d = d[d.get("production_ps", 0) >= 50000].copy()
    if "status" in d.columns:
        d = d[d.status.astype(str).str.startswith("ok")]
    eng = "explicit_frac_frames_engaged"
    if eng in d.columns:
        d["_has"] = d[eng].notna().astype(int)
        d = d.sort_values("_has")
    # DEDUPED PER (MOLECULE, MODE), NOT PER MOLECULE. Two modes of one molecule
    # are two different complexes with two different trajectories; collapsing
    # them on `ident` keeps one arbitrarily and silently discards the other --
    # the same "value taken by label, not identity" failure that gave every mode
    # of a molecule the same sweep trajectory. No molecule has two 100 ns modes
    # yet, so this changes nothing today; it is here so that the first one does
    # not have to be noticed.
    d["_mode"] = (d["pose_mode"] if "pose_mode" in d.columns
                  else pd.Series(index=d.index, dtype="float64")).fillna(-1)
    return d.drop_duplicates(["ident", "_mode"], keep="last")


def _version() -> tuple[str, str]:
    """Current version and its codename, read from the CHANGELOG's first entry.

    Not a constant in this file. A version literal here would be a pin, and pins
    in this repo go stale silently -- the CHANGELOG is where the release is
    actually declared, so ask it.
    """
    import re
    p = REPO / "CHANGELOG.md"
    if p.is_file():
        m = re.search(r'^##\s+([0-9]+\.[0-9]+\.[0-9]+)\s+[""“”"]?([^""“”"\n—-]*)',
                      p.read_text(), re.M)
        if m:
            return m.group(1), m.group(2).strip().strip('"“”')
    return "", ""


def _sweep_all() -> pd.DataFrame:
    """Every sweep row, INCLUDING the failures.

    ``_sweep()`` drops ``status != "ok"``, which is right for candidates and wrong
    for controls: a control that could not be swept is a result, and the reason is
    the interesting part. Keeps the ok row per molecule when there is one, so a
    control that succeeded on a re-run is not represented by its earlier failure.
    """
    fs = sorted(glob.glob(str(rp.sweep_dir() / "attack_sweep_*.csv")),
                key=lambda p: int(p.rsplit("_", 1)[1].split(".")[0]))
    if not fs:
        return pd.DataFrame()
    d = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    d = d[d.get("sweep_ps", 0) > 1000].copy()
    d["_ok"] = (d.status == "ok").astype(int)
    return (d.sort_values(["_ok", "frac_attack_ready"], ascending=[False, False])
             .drop_duplicates("parent_ident"))


def _controls() -> list[dict]:
    """The experimentally-determined poses, put through our own criterion.

    Two kinds, and the distinction is the entire point of the experiment (#47,
    ``crystal_controls.py``):

    * ``xtal_*`` — the deposited geometry, still BONDED to Cys113 SG at ~2.0 Å.
      That is the reaction PRODUCT. The near-attack window is 2.8–4.2 Å, so a
      bonded pose cannot be attack-ready **by construction**. These carry no
      sweep value, and the reason is shown instead of a number.
    * ``rx_*`` — the same crystal pose with the bond cleaved and the leaving group
      rebuilt. That IS a pre-reaction geometry, and it sweeps.

    Only the ``rx_*`` forms are commensurate with a candidate, so only those are
    ranked beside them. Stamping the rest rather than dropping them is deliberate:
    "could not be swept" and "swept badly" are different facts, and a control that
    silently vanished from the page would be indistinguishable from one that was
    never run.
    """
    T = sout.Topic("blacksmith", "crystal_controls")
    out: list[dict] = []
    swa = _sweep_all()
    swi = swa.set_index("parent_ident") if not swa.empty else pd.DataFrame()

    def _row(ident):
        return swi.loc[ident] if ident in getattr(swi, "index", []) else None

    for stem, kind in (("crystal_controls", "bonded"), ("crystal_reactant", "reactant")):
        try:
            p = T.latest(stem, ".csv")
        except Exception:                                  # noqa: BLE001
            log.warning("no %s output yet", stem)
            continue
        d = pd.read_csv(p)
        for _, r in d.iterrows():
            ident = str(r.get("ident", "") or "")
            if not ident:
                continue
            s = _row(ident)
            ar = None
            if s is not None and "frac_attack_ready" in s and not pd.isna(s["frac_attack_ready"]):
                ar = float(s["frac_attack_ready"])
            out.append({
                "ident": ident,
                "kind": kind,
                "pdb": str(r.get("pdb", "")),
                "label": str(r.get("name", "") or r.get("comp_id", "") or ""),
                "prep_status": str(r.get("status", "")),
                "sweep_status": (str(s["status"]) if s is not None else "not swept"),
                "frac_attack_ready": ar,
                "n_visits": (float(s["n_visits"]) if s is not None
                             and "n_visits" in s and not pd.isna(s["n_visits"]) else None),
                "dist_a": (float(r["built_x_to_sg_a"]) if not pd.isna(r.get("built_x_to_sg_a", float("nan")))
                           else (float(r["linked_atom_to_sg_a"])
                                 if not pd.isna(r.get("linked_atom_to_sg_a", float("nan"))) else None)),
                "fit_rmsd_a": (float(r["fit_rmsd_a"]) if not pd.isna(r.get("fit_rmsd_a", float("nan"))) else None),
            })
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--candidates", nargs="+", required=True)
    ap.add_argument("--title", default="DWI Derivative Screen")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    # Read BEFORE anything uses them. These were set further down next to the
    # page template, which put the `global` statement after the first read of
    # _HELD_BAR in the rail loop -- a SyntaxError, so the page silently did not
    # rebuild while the sweep page beside it did.
    from shared import target_config as _tc
    global _HELD_BAR
    _HELD_BAR = float(_tc.get("md.sweep_survivor_rmsd_nm", default=0.35))
    _SWEEP_NS = float(_tc.get("md.sweep_ps", default=8000)) / 1000.0
    log.info("qualifying bar for a 100 ns hold: %.2f nm; triage %.0f ns (config)",
             _HELD_BAR, _SWEEP_NS)

    sw, md, cls_of = _sweep(), _md(), _classes()
    # Controls are resolved BEFORE the depictions so their structures are drawn
    # in the same pass as the candidates', from the same function.
    ctl = _controls()
    thumbs = _thumbs(list(args.candidates) + [c["ident"] for c in ctl])
    # The crystal set from _controls(), PLUS the reference molecules. ref_* are
    # known binders put through the identical criterion -- they are controls by
    # construction, and tagging only the crystal ones left Juglone ranked but
    # absent from the controls tab.
    ctl_idents = {c['ident'] for c in ctl} | {t for t in args.candidates if t.startswith('ref_')}
    # (args.candidates is already filtered to qualifying controls above)
    _ver, _code = _version()
    # The version AND its codename belong in the title (@tt8804): a page saved,
    # screenshotted or pasted into a thread carries its release with it rather
    # than losing it to a chip someone cropped out.
    _full_title = " ".join(x for x in (
        args.title,
        f"version {_ver}" if _ver else "",
        f"“{_code}”" if _code else "") if x)
    swi = sw.set_index("parent_ident") if not sw.empty else pd.DataFrame()
    mdi = md.set_index("ident") if not md.empty else pd.DataFrame()

    # WHAT EARNS A ROW. Three kinds qualify and nothing else:
    #   * a candidate (t3_/t4_) that has a measurement -- a sweep or a 100 ns run;
    #   * a crystal control (rx_/xtal_);
    #   * a reference that is a confirmed covalent Cys113 inhibitor carrying a
    #     warhead we enumerate (_qualifying_controls).
    #
    # This drops two things that were cluttering the rail. Reference molecules
    # that are not comparable -- ATRA and the non-covalent binders, which have no
    # warhead to aim -- and rows carrying no number at all, like a 190 ns variant
    # whose trajectory is long gone: it showed a blank tile, an em-dash and
    # "awaiting 100 ns" for a run that is never coming.
    qual = _qualifying_controls()
    swept_ids = set(swi.index.astype(str)) if len(swi.index) else set()
    md_ids = set(mdi.index.astype(str)) if len(mdi.index) else set()

    def keep(c: str) -> tuple[bool, str]:
        if c.startswith(("rx_", "xtal_")):
            return True, ""
        if c.startswith("ref_"):
            base = c.split("__")[0]
            return (base in qual), "not a covalent control with a warhead we use"
        if c in swept_ids or c in md_ids:
            return True, ""
        return False, "no sweep and no 100 ns run"

    dropped = []
    kept = []
    for c in args.candidates:
        ok, why = keep(c)
        (kept if ok else dropped).append(c if ok else (c, why))
    if dropped:
        for c, why in dropped:
            log.info("dropped %s — %s", c, why)
    args.candidates = kept

    # `rows` (a second, per-candidate <table>) was built here and never
    # rendered -- the page has used `rows_html` since the rail replaced it.
    # Dead code that formats data is worse than absent: it reads as the
    # thing on screen, so a fix applied to it changes nothing visible.
    tabs, missing = [], []
    for c in args.candidates:
        f = REPORTS / f"{c}.html"
        if not f.is_file():
            missing.append(c)
            log.warning("%s: no report at %s", c, f.name)
            continue
        # This loop selects; the rail below formats. It used to do both, and the
        # formatting half fed a <table> that stopped being rendered when the rail
        # replaced it -- so the lookups and the `g()` here were computed every
        # run for nothing.
        tabs.append(c)

    if not tabs:
        raise SystemExit("no reports found for any requested candidate")

    btns = "".join(
        f"<button onclick=\"show('{html.escape(t)}')\" id='b_{html.escape(t)}'>"
        f"{html.escape(t)}</button>" for t in tabs)
    miss = (f"<p class='warn'>No report yet for: {', '.join(map(html.escape, missing))}"
            " — still running, or the trajectory is incomplete.</p>") if missing else ""

    # LEFT RAIL SELECTOR, RIGHT VIEWER (@tt8804) -- the same shape as the GUI's
    # ranking panel, because that is the layout the reading actually happens in:
    # you scan the list, click, and the pose/movie/plots replace themselves beside
    # it. A row of buttons across the top pushed the viewer below the fold and
    # made comparison a scroll.
    rows_html = []
    for k, t in enumerate(tabs):
        # THE IDENT IS A DISPLAY LABEL; THE JOINS ARE ON THE PARENT. Once the
        # runner keys a run `<parent>_m<mode>`, everything looked up by molecule
        # -- the sweep row, the warhead class, the depiction -- must be looked up
        # under the parent or it silently misses and the row renders as a
        # molecule nobody has any information about. Legacy idents have no
        # suffix, so `split_ident` returns them unchanged and they behave exactly
        # as before.
        par, mode_of_ident = mode_key.split_ident(t)
        s_ = swi.loc[par] if par in getattr(swi, "index", []) else None
        m_ = mdi.loc[t] if t in getattr(mdi, "index", []) else None
        # `.loc` on a duplicated index returns a FRAME, and every `g()` below
        # would format a Series into a cell. `_md` now keeps one row per
        # (molecule, mode), so a molecule with two 100 ns modes lands here.
        if isinstance(m_, pd.DataFrame):
            _k = "explicit_ligand_rmsd_nm_max"
            m_ = (m_.sort_values(_k).iloc[0] if _k in m_.columns else m_.iloc[0])
        if isinstance(s_, pd.DataFrame):
            s_ = s_.iloc[0]
        # WHICH POSE THIS RUN STARTED FROM -- from the run's own row. `pose_mode`
        # is written by the runner out of --pose-rank and is the only record of
        # it. The rail showed the molecule alone, so two modes of one molecule
        # were indistinguishable on screen. @tt8804: "it doesnt show sub modes
        # only molecules".
        #
        # Present on 3 of 70 rows. Everything before today was launched per
        # MOLECULE, and the mode cannot be recovered after the fact either --
        # matching each run's ligand_pose.sdf against the current representatives
        # puts only 10 of 83 within 0.05 A, the rest having been built from poses
        # that no longer exist. So the others read "mode not recorded", which is
        # the true state; borrowing a label from a same-molecule sweep row would
        # make an unknown look like a measurement.
        # Preferring the IDENT over the column: once the mode is part of the
        # identity it is the run's own name, whereas `pose_mode` is a value
        # beside it that a later rewrite could leave stale. They agree by
        # construction; when only one exists, use whichever does.
        mode = mode_of_ident
        if mode is None and m_ is not None and "pose_mode" in m_ \
                and pd.notna(m_["pose_mode"]):
            mode = int(m_["pose_mode"])
        def g(src, key, fmt="{:.3f}"):
            if src is None or key not in src or pd.isna(src[key]):
                return "\u2014"
            try:
                return fmt.format(src[key])
            except Exception:                          # noqa: BLE001
                return str(src[key])
        # THE 10 ns SWEEP IS TRIAGE, NOT THE RESULT. It exists to choose which
        # molecules earn a 100 ns run. Ranking on it ranks the SELECTION FILTER
        # and not the endpoint -- the same shape as ranking on docking energy,
        # one stage further down. The endpoint is 100 ns target engagement, so
        # that is the sort key, and the sweep is carried beside it as the triage
        # reading it is.
        ar = None
        if s_ is not None and "frac_attack_ready" in s_ and not pd.isna(s_["frac_attack_ready"]):
            ar = float(s_["frac_attack_ready"])
        eng = None
        if (m_ is not None and "explicit_frac_frames_engaged" in m_
                and not pd.isna(m_["explicit_frac_frames_engaged"])):
            eng = float(m_["explicit_frac_frames_engaged"])
        rmax = None
        if m_ is not None and "explicit_ligand_rmsd_nm_max" in m_ and not pd.isna(m_["explicit_ligand_rmsd_nm_max"]):
            rmax = float(m_["explicit_ligand_rmsd_nm_max"])
        # THE QUALIFYING BAR, FROM CONFIG (D0085). 0.35 nm -- the same number
        # the 8 ns sweep uses and the same one BPMD promotion reads, because it
        # is the same question at a longer timescale. It was hardcoded 1.2 here,
        # which is the old "did not dissociate" reading and three times looser
        # than what now earns a molecule anything.
        held = rmax is not None and rmax < _HELD_BAR
        has_md = eng is not None
        wcls = str(cls_of.get(par, cls_of.get(t, "unclassified")))
        # A molecule with no 100 ns run cannot be placed on the ranked axis at
        # all. It goes in its own band rather than being given a 0, which would
        # read as "measured and engaged nothing".
        # THE HEADLINE IS WHAT THE RAIL IS SORTED ON (@tt8804, #55): max ligand
        # RMSD over the 100 ns run, lowest first. Engagement moves to the meta
        # line -- still shown, no longer the thing being ordered.
        headline = (f"{rmax:.3f} nm max" if has_md and rmax is not None
                    else ("—" if not has_md else f"{eng*100:.0f}% engaged"))
        # THE SELECTOR CARRIES 100 ns FACTS ONLY (@tt8804, #55): max ligand RMSD,
        # held/left, engaged %. The 10 ns sweep is triage for deciding what earns
        # a 100 ns run -- it is not a result, and sitting in the rail beside the
        # engagement number it read as a second, competing score. It moves to a
        # table in the viewer, where it is clearly labelled as what selected the
        # molecule rather than what was found.
        meta = (f"{eng*100:.0f}% engaged" if has_md and eng is not None
                else "awaiting 100 ns")
        # A CONTROL THAT NOW HAS A 100 ns RUN IS ONE ROW, NOT TWO. It was being
        # emitted here as a ranked candidate AND again in the controls block as an
        # unranked control -- two rows, the same DOM id twice, so clicking either
        # resolved to whichever the browser found first. It stays a single row
        # that is ranked on its own number and still answers the controls tab.
        is_ctl = t in ctl_idents
        rows_html.append(
            f"<button class='row{' ctl' if is_ctl else ''}' "
            f"data-cls=\"{'control' if is_ctl else html.escape(wcls)}\" "
            + ("data-ctl='1' " if is_ctl else "")
            + f"data-eng='{(eng if has_md else -1):.6f}' "
            # THE SORT KEY, ASCENDING -- lower is better, so an unranked row
            # cannot be given 0. It gets a sentinel that sorts last, and is in
            # the unranked band anyway.
            f"data-rmax='{(rmax if has_md and rmax is not None else 9999):.6f}' "
            f"data-sweep='{(ar if ar is not None else -1):.6f}' "
            f"data-md='{1 if has_md else 0}' "
            f"data-held='{1 if held else 0}' "
            f"id='b_{html.escape(t)}' onclick=\"show('{html.escape(t)}')\">"
            f"<span class='rk'>{k+1}</span>"
            + (f"<img class='thumb' alt='' src=\"{thumbs.get(t) or thumbs[par]}\">"
               if (t in thumbs or par in thumbs) else "<span class='thumb'></span>")
            + f"<span class='body'>"
            f"<span class='l1'><span class='mid-id'>{html.escape(t)}"
            + (f"<span class='mode'>m{mode}</span>" if mode is not None
               else "<span class='mode unk' title='this run was launched per "
                    "molecule; which pose it started from was not recorded and "
                    "cannot be recovered'>mode ?</span>")
            + f"</span>"
            f"<span class='eng{'' if has_md else ' pend'}' "
            f"title='{'largest ligand RMSD over the 100 ns run — the sort key, lowest first' if has_md else 'no 100 ns run yet — 10 ns triage sweep only'}'>"
            f"{headline}</span></span>"
            f"<span class='l2'><span class='wc'>{html.escape(wcls)}</span>"
            f"<span class='meta'>{meta}</span>"
            + (f"<span class='tag {'t-held' if held else 't-left'}'>"
               f"{'held' if held else 'left'}</span>"
               if has_md else "<span class='tag t-pend'>swept</span>")
            + "</span>"
            f"<span class='bar'><i style='width:{max(1.5,(eng if has_md else 0)*100):.1f}%'></i></span>"
            f"</span></button>")

    # CONTROLS (#47/#48). The catalogue previously showed candidates only, so a
    # ranking under active falsification had none of the falsifying evidence on
    # screen. Controls that swept are ranked WITH the candidates on the identical
    # axis; controls that could not be swept appear in the controls tab carrying
    # the reason. Neither is given a held/left tag: no control has 100 ns, so
    # "left" would be a fact we do not have.
    ranked_ctl = [c for c in ctl if c["frac_attack_ready"] is not None]
    ctl_rows_html = []
    for c in ctl:
        if c["ident"] in tabs:
            continue      # already a ranked row above
        ar = c["frac_attack_ready"]
        has = ar is not None
        cid = c["ident"]
        nm = c["label"] or c["pdb"]
        why = ("bonded product ~{:.2f} A from SG — outside the 2.8-4.2 A "
               "near-attack window by construction".format(c["dist_a"])
               if c["kind"] == "bonded" and c["dist_a"] else c["sweep_status"])
        # A control with its own report opens in the SAME viewer as a candidate --
        # pose, movie, RMSD plots, identical layout. Only the ones with no report
        # fall back to the controls page.
        has_rep = (REPORTS / f"{cid}.html").is_file()
        ctl_rows_html.append(
            f"<button class='row ctl' data-ctl='1' data-kind='{c['kind']}' "
            f"data-cls='control' data-eng='-1' "
            f"data-sweep='{(ar if has else -1):.6f}' data-md='0' "
            f"data-noval='{0 if has else 1}' data-held='0' "
            + ("" if has_rep else "data-src='controls.html' ")
            + f"id='b_{html.escape(cid)}' onclick=\"show('{html.escape(cid)}')\">"
            f"<span class='rk'>&middot;</span>"
            + (f"<img class='thumb' alt='' src=\"{thumbs[cid]}\">"
               if cid in thumbs
               else f"<span class='thumb tctl'>{'RX' if c['kind']=='reactant' else 'XT'}</span>")
            + f"<span class='body'>"
            f"<span class='l1'><span class='mid-id'>{html.escape(nm)}</span>"
            f"<span class='eng pend' title='control — 10 ns only, no 100 ns run'>"
            f"&mdash;</span></span>"
            f"<span class='l2'><span class='wc'>control &middot; {html.escape(c['pdb'])}</span>"
            f"<span class='meta'>"
            + ("awaiting 100 ns" if c["n_visits"] is not None else html.escape(why))
            + "</span>"
            f"<span class='tag t-ctl'>{'reactant' if c['kind']=='reactant' else 'bonded'}</span>"
            f"</span>"
            f"<span class='bar'><i style='width:{max(1.5, (ar or 0)*100):.1f}%'></i></span>"
            f"</span></button>")

    # The controls page the iframe loads. Self-contained, same palette, and it
    # states the interpretation rather than leaving the reader to infer it.
    def _crow(c):
        def num(v, fmt="{:.4f}"):
            return fmt.format(v) if v is not None else "&mdash;"
        return ("<tr>"
                f"<td class='id'>{html.escape(c['label'] or c['pdb'])}</td>"
                f"<td>{html.escape(c['pdb'])}</td>"
                f"<td>{html.escape(c['kind'])}</td>"
                f"<td class='n'>{num(c['dist_a'], '{:.2f}')}</td>"
                f"<td class='n'>{num(c['frac_attack_ready'])}</td>"
                f"<td class='n'>{num(c['n_visits'], '{:.0f}')}</td>"
                f"<td class='st'>{html.escape(c['sweep_status'])}</td>"
                "</tr>")

    # Read the candidate range off the DATA, not by parsing it back out of the
    # markup we just wrote -- a string splice on generated HTML is precisely the
    # take-it-by-position defect this project keeps finding.
    cand_ar = []
    for t in tabs:
        s_ = swi.loc[t] if t in getattr(swi, "index", []) else None
        if s_ is not None and "frac_attack_ready" in s_ and not pd.isna(s_["frac_attack_ready"]):
            cand_ar.append(float(s_["frac_attack_ready"]))
    best_cand = max(cand_ar, default=0.0)
    worst_cand = min(cand_ar, default=0.0)
    # Look the controls up BY IDENT. Indexing ranked_ctl[0]/[-1] happened to give
    # the right pair only because of the order the two source files are read in --
    # a positional read of an identity, which is the defect this repo is named for.
    by_ident = {c["ident"]: c for c in ctl}
    ctl_best = max((c["frac_attack_ready"] for c in ranked_ctl), default=None)
    above = sum(1 for a in cand_ar if ctl_best is not None and a > ctl_best)

    def _ctl_ar(ident):
        c = by_ident.get(ident)
        v = c["frac_attack_ready"] if c else None
        return f"{v:.4f}" if v is not None else "&mdash;"
    ctl_table = "".join(_crow(c) for c in ctl)
    ctl_page = f"""<!doctype html><html><head><meta charset="utf-8">
<title>controls</title><style>
:root{{--ink:#10233f;--navy:#003087;--blue:#0072ce;--rule:#d6dee8;--muted:#5b6b80;
 --paper:#fff;--raise:#f5f8fc;--bad:#b3261e;
 --sans:"Helvetica Neue",Helvetica,Arial,system-ui,sans-serif;
 --mono:ui-monospace,"SF Mono",Menlo,Consolas,monospace}}
:root[data-theme="dark"]{{--ink:#dfe7f0;--navy:#8ab4e8;--blue:#6aa9e0;--rule:#25333f;
 --muted:#93a3b4;--paper:#0e151c;--raise:#16202a;--bad:#e08a70}}
*{{box-sizing:border-box}}
body{{margin:0;padding:22px 26px;background:var(--paper);color:var(--ink);
 font-family:var(--sans);font-size:14px;line-height:1.55;
 font-variant-numeric:tabular-nums;max-width:62rem}}
h1{{font-size:1.15rem;color:var(--navy);margin:0 0 2px}}
h2{{font-size:.62rem;font-family:var(--mono);letter-spacing:.14em;text-transform:uppercase;
 color:var(--blue);margin:26px 0 8px}}
p{{margin:.5em 0}}
table{{border-collapse:collapse;width:100%;margin-top:6px;font-size:13px}}
th,td{{text-align:left;padding:6px 9px;border-bottom:1px solid var(--rule)}}
th{{font-family:var(--mono);font-size:.58rem;letter-spacing:.1em;text-transform:uppercase;
 color:var(--muted);font-weight:700}}
td.n{{text-align:right;font-family:var(--mono)}}
td.id{{font-family:var(--mono);font-weight:600}}
td.st{{color:var(--muted);font-size:12px}}
.lead{{color:var(--muted);margin-bottom:4px}}
.box{{border-left:3px solid var(--bad);background:var(--raise);padding:10px 14px;
 margin:14px 0;border-radius:0 4px 4px 0}}
code{{font-family:var(--mono);font-size:12.5px}}
</style></head><body>
<h1>Controls — the experimentally-determined poses, through our own criterion</h1>
<p class="lead">Deposited Pin1 complexes with a covalent bond to Cys113, scored by the
same near-attack criterion that ranks every candidate in this catalogue.</p>

<h2>Why this exists</h2>
<p><a href="https://github.com/hallettmiket/inhibition/issues/47">#47</a> measured that
the warhead classes with crystal structures and measured kinetics score <em>last</em>
on our near-attack criterion, while a class with no measured Pin1 activity scores
first. Two explanations fit that equally well: our <strong>docking</strong> puts
those molecules in the wrong place so the criterion never sees the real geometry,
or our <strong>criterion</strong> is wrong and would reject the real geometry too.
Only one experiment separates them — take the pose the crystallographer determined
and score it.</p>

<h2>Two forms, and the difference matters</h2>
<p><strong>bonded</strong> — the deposited geometry, still attached to Cys113 SG at
~2 Å. That is the reaction <em>product</em>. The near-attack window is 2.8–4.2 Å, so
a bonded pose cannot be attack-ready <em>by construction</em>; it produces no
attack-geometry series and no number. That is a property of the experiment, not a
failure of the molecule.</p>
<p><strong>reactant</strong> — the same crystal pose with the bond cleaved and the
leaving group rebuilt. This is a genuine pre-reaction geometry and it sweeps, so it
is directly comparable to a candidate and is ranked beside them.</p>

<table><tr><th>molecule</th><th>pdb</th><th>form</th><th>d(X&rarr;SG) &Aring;</th>
<th>attack-ready</th><th>visits</th><th>sweep status</th></tr>
{ctl_table}</table>

<div class="box">
<p><strong>What the reactant controls say.</strong> The two rebuildable controls are
Sulfopin (6VAJ) and Liu-2022-ZL-Pin13 (7F0M) — both crystallographically bound to
Cys113, both with measured potency. Through our criterion they score
<code>{_ctl_ar('rx_6VAJ')}</code> and <code>{_ctl_ar('rx_7F0M')}</code>
attack-ready, with <strong>zero sustained visits</strong>.</p>
<p>The candidates in this catalogue run from <code>{best_cand:.4f}</code> down to
<code>{worst_cand:.4f}</code>, and <strong>{above} of {len(cand_ar)}</strong> of them
score above the better of the two controls. So the criterion ranks most of our
generated matter ahead of chemistry that is known to react with Cys113.</p>
<p>Read against #47, that points away from "docking mislocates these molecules" and
toward the criterion itself rejecting geometry that is known to react. It does not
settle it — the reactant forms are <em>rebuilt</em>, not observed, and the rebuild
places the leaving group. But a criterion that scores the answer key near zero is
not yet evidence that the molecules above it are better.</p>
</div>

<h2>What would settle it</h2>
<p>More rebuildable controls. Only 2 of the 6 covalent complexes carry a halogen
leaving group the reactant builder can restore; the other four are recorded above
with their bonded distance and no sweep. Extending the builder to the remaining
chemistries is the cheapest way to turn two points into a distribution.</p>
</body></html>"""
    (REPORTS / "controls.html").write_text(ctl_page)
    # The schematic carries the SAME title as the GUI it explains, rather than a
    # heading of its own that drifts the moment the release name changes.
    # pipeline.html is written to a STABLE name, not through outputs.py's
    # versioned writer, so stamping the build date here does not spawn a new
    # versioned file per day the way it would for the frames.
    from datetime import date as _date
    (REPORTS / "pipeline.html").write_text(
        schematic.build(_full_title, _date.today().isoformat()))
    # EVERY MODE, ranked individually (#53). The rail is one row per MOLECULE
    # because it indexes the sweep by parent_ident, so the per-mode ranking the
    # pipeline computes was never visible anywhere. Same stable-name treatment
    # as the schematic.
    # Assets first: the ranking view fetches a depiction and a pose per molecule
    # rather than inlining 8,096 of each. Existing files are left alone.
    from shared import mode_assets as massets
    # Vendored ONCE into the ranking page's <head>, before any viewer script --
    # loaded at the end of the body it is not defined yet when the viewer looks
    # for it, and the viewer then draws nothing without erroring.
    _three_js = (REPO / "scripts" / ".cache_3dmol-min.js")
    three = _three_js.read_text() if _three_js.is_file() else ""
    _mr = moderank.gather()
    # HOW MANY MODELS EACH ASSET MUST HOLD. An asset left by an earlier run with
    # fewer poses than the ranking has modes is STALE, and skipping it on mere
    # existence is what left 195 molecules showing an empty viewer for every mode
    # above 0.
    _need = {}
    if not _mr.empty and "parent_ident" in _mr.columns:
        for _p, _g in _mr.groupby("parent_ident"):
            _need[str(_p)] = int(_g["mode"].max()) + 1
    _a = massets.write_assets(REPORTS, moderank.idents(_mr),
                              force=os.environ.get("MODE_ASSETS_FORCE") == "1",
                              expected=_need)
    log.info("mode assets: +%d poses, +%d thumbs, %d molecules with no pose "
             "from this run", _a["poses"], _a["thumbs"], len(_a.get("stale", [])))
    (REPORTS / "modes.html").write_text(
        moderank.build(_full_title, _date.today().isoformat(), three,
                       no_pose=_a.get("stale", [])))

    # The stepper (#63), from the one module that defines it. Counts come from
    # the same places the pages do, so the nav cannot disagree with the page it
    # sits on.
    from shared import gui_shell as _gs
    from shared import results_shell as _rs
    # The shell's layout, from the module that owns it. The sweep page uses the
    # same constant, so the two cannot drift into two different designs.
    _shellcss = _rs.CSS
    _stepcss = _gs.CSS
    _stepnav = _gs.nav("combined.html", moderank._step_counts(_mr))

    page = f"""<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(_full_title)}</title><style>
{_shellcss}
{_stepcss}
</style></head><body>
<div id="topbar">
 <h1 title="Pick a molecule on the left; its pose, movie and plots load on the right.">{html.escape(_full_title)}</h1>
 <span class="msep"></span>
 <button id="m-all" class="mbtn on" onclick="setMode('all')">all classes</button>
 <button id="m-cls" class="mbtn" onclick="setMode('cls')">by warhead class</button>
 <span class="msep"></span>
 <button id="o-mix" class="mbtn on" onclick="setSplit(0)">combined</button>
 <button id="o-spl" class="mbtn" onclick="setSplit(1)">split held / left</button>
 <span class="msep"></span>
 <button id="c-tab" class="mbtn" onclick="setTab()" title="the crystallographic controls, scored by the same criterion">controls</button>
 <span class="mhint" id="mhint"></span>
 <a class="mbtn lnk" href="pipeline.html" target="_blank" rel="noopener"
    title="how a molecule becomes a row: docking, modes, criteria, ranking, sweep, MD">how this works &#8599;</a>
 <button id="theme" class="mbtn tbtn" onclick="toggleTheme()" title="light / dark">dark</button>
</div>
<!-- THE SAME STEPPER AS EVERY OTHER PAGE (#63). The "every mode" link is gone
     from the topbar: it opened the Ranking page in a NEW TAB, which is the
     opposite of one instrument with four steps. -->
{_stepnav}
<main>
 <div id="rail">{_rs.SEARCH_HTML}<div class="legend">ranked by <b>max ligand RMSD</b> over the
  100&nbsp;ns run, lowest first &mdash; how far the molecule ever got from where it
  started. <b>held</b> means it never exceeded {_HELD_BAR:.2f}&nbsp;nm, the bar
  that earns a molecule BPMD (D0085) and the same one the {int(_SWEEP_NS)}&nbsp;ns
  triage applies. Engagement is shown beside it. The triage sweep chooses what
  earns a 100&nbsp;ns run and is not the result; it is on each molecule&rsquo;s
  own page.</div>
 {''.join(rows_html)}{''.join(ctl_rows_html)}</div>
 <div id="viewer">
  <div id="vhead"><span id="vname">&mdash;</span>
   <a id="vopen" href="#" target="_blank" rel="noopener">open full report &#8599;</a></div>
  <iframe id="v" title="molecule report" src="{html.escape(tabs[0])}.html"></iframe>
 </div>
</main>
<script>
var RAIL=document.getElementById('rail');
var ROWS=Array.prototype.slice.call(RAIL.querySelectorAll('.row'));
var MODE='all', SPLIT=0, TAB=false;
function renumber(l){{l.forEach(function(b,i){{b.querySelector('.rk').textContent=i+1}});}}
function hdr(cls,txt){{var h=document.createElement('div');h.className=cls;h.textContent=txt;
  RAIL.appendChild(h);}}
// ASCENDING on max ligand RMSD: lowest excursion first. Every other ordering in
// this file is descending-is-better, so the direction is stated rather than left
// to be inferred from a minus sign.
function byRank(a,b){{return parseFloat(a.dataset.rmax)-parseFloat(b.dataset.rmax)}}
function bySweep(a,b){{return parseFloat(b.dataset.sweep)-parseFloat(a.dataset.sweep)}}
function layoutGroup(rows){{
  if(MODE==='all'){{rows.forEach(function(b){{RAIL.appendChild(b)}});renumber(rows);return;}}
  var g={{}};
  rows.forEach(function(b){{(g[b.dataset.cls]=g[b.dataset.cls]||[]).push(b)}});
  Object.keys(g).sort(function(x,y){{return byRank(g[x][0],g[y][0])}}).forEach(function(n){{
    hdr('chd',n+'  ('+g[n].length+')');
    g[n].forEach(function(b){{RAIL.appendChild(b)}}); renumber(g[n]);
  }});
}}
function relayout(){{
  RAIL.querySelectorAll('.chd,.ohd').forEach(function(h){{h.remove()}});
  // A control with no sweep value cannot be RANKED -- there is no number to sort
  // it by -- but it must still be visible somewhere, so it lives in the controls
  // tab. Dropping it entirely would make "could not be swept" look identical to
  // "was never run".
  var pool=ROWS.filter(function(b){{
    return TAB ? b.dataset.ctl==='1' : b.dataset.noval!=='1';
  }});
  ROWS.forEach(function(b){{b.style.display='none'}});
  pool.forEach(function(b){{b.style.display=''}});
  var all=pool.slice().sort(byRank);
  // ONLY A 100 ns RUN CAN BE RANKED. The sweep is the triage that decides what
  // earns one, so a swept-but-not-yet-run molecule has no position on this axis
  // -- it gets its own band, ordered by its sweep reading, rather than a zero
  // that would read as "measured, engaged nothing".
  var done=all.filter(function(b){{return b.dataset.md==='1'}});
  var pend=all.filter(function(b){{return b.dataset.md!=='1'}})
              .sort(bySweep);
  if(TAB){{ all.sort(bySweep).forEach(function(b){{RAIL.appendChild(b)}}); renumber(all); }}
  else{{
    if(!SPLIT){{ layoutGroup(done); }}
    else{{
      var held=done.filter(function(b){{return b.dataset.held==='1'}});
      var gone=done.filter(function(b){{return b.dataset.held!=='1'}});
      if(held.length){{hdr('ohd o-held','held the pocket  ('+held.length+')'); layoutGroup(held);}}
      if(gone.length){{hdr('ohd o-left','dissociated  ('+gone.length+')'); layoutGroup(gone);}}
    }}
    // Controls sit here too: none has a 100 ns trajectory, so none can be ranked
    // and none carries a held/left verdict.
    if(pend.length){{
      var nctl=pend.filter(function(b){{return b.dataset.ctl==='1'}}).length;
      hdr('ohd o-pend','10 ns sweep only \u2014 not yet ranked  ('+pend.length
          +(nctl?', incl. '+nctl+' controls':'')+')');
      pend.forEach(function(b){{RAIL.appendChild(b)}}); renumber(pend);
    }}
  }}
  var bits;
  if(TAB){{
    bits=[pool.length+' controls','crystallographic poses through the same criterion'];
  }} else {{
    bits=[done.length+' ranked on 100 ns engagement'];
    if(pend.length) bits.push(pend.length+' swept only');
    bits.push(MODE==='all'?'one ranking across all classes'
      :'ranked within warhead class \u2014 cross-class comparison is biased (#47)');
    if(SPLIT) bits.push('held and dissociated shown separately');
  }}
  document.getElementById('mhint').textContent=bits.join(' \u00b7 ');
}}
function setTab(){{TAB=!TAB;
  document.getElementById('c-tab').classList.toggle('on',TAB);
  ['m-all','m-cls','o-mix','o-spl'].forEach(function(i){{
    document.getElementById(i).disabled=TAB;
  }});
  relayout();
  if(TAB){{
    var f=document.getElementById('v');
    f.onload=function(){{applyTheme(f.contentDocument)}};
    f.src='controls.html';
    document.getElementById('vname').textContent='controls';
    document.getElementById('vopen').href='controls.html';
  }}
}}
function setMode(m){{MODE=m;
  document.getElementById('m-all').classList.toggle('on',m==='all');
  document.getElementById('m-cls').classList.toggle('on',m==='cls');
  relayout();}}
function setSplit(v){{SPLIT=v;
  document.getElementById('o-mix').classList.toggle('on',!v);
  document.getElementById('o-spl').classList.toggle('on',!!v);
  relayout();}}
function applyTheme(doc){{
  try{{ doc.documentElement.setAttribute('data-theme',
        document.documentElement.getAttribute('data-theme')||'light'); }}catch(e){{}}
}}
function toggleTheme(){{
  var cur=document.documentElement.getAttribute('data-theme')==='dark'?'light':'dark';
  document.documentElement.setAttribute('data-theme',cur);
  document.getElementById('theme').textContent=cur==='dark'?'light':'dark';
  try{{localStorage.setItem('cat-theme',cur)}}catch(e){{}}
  var f=document.getElementById('v');
  if(f&&f.contentDocument) applyTheme(f.contentDocument);
}}
(function(){{
  var saved='light';
  try{{saved=localStorage.getItem('cat-theme')||'light'}}catch(e){{}}
  document.documentElement.setAttribute('data-theme',saved);
  document.addEventListener('DOMContentLoaded',function(){{
    document.getElementById('theme').textContent=saved==='dark'?'light':'dark';
  }});
}})();
function show(t){{
  var f=document.getElementById('v');
  var el0=document.getElementById('b_'+t);
  // Controls have no per-molecule report; route them to the controls page rather
  // than letting the iframe 404 into a blank pane.
  var src=(el0&&el0.dataset.src)?el0.dataset.src:t+'.html';
  f.onload=function(){{applyTheme(f.contentDocument)}};
  f.src=src;
  document.getElementById('vname').textContent=t;
  document.getElementById('vopen').href=src;
  document.querySelectorAll('.row').forEach(function(b){{b.classList.remove('on')}});
  var el=document.getElementById('b_'+t); if(el){{el.classList.add('on');}}
}}
relayout();
show({json.dumps(tabs[0])});
{_rs.SEARCH_JS}
railFilter();
</script></body></html>"""

    dest = sout.Topic("blacksmith", rp.reports_topic()).write("combined", ".html")
    dest.write_text(page)
    # Also drop a stable name beside the reports so the iframes resolve relatively.
    side = REPORTS / "combined.html"
    side.write_text(page)
    print(f"\n  {len(tabs)} reports combined -> {side}")
    print(f"  versioned copy: {dest}")


if __name__ == "__main__":
    main()
