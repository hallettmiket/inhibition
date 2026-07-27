"""
Purpose: T_4 step 7 — warhead reactivity triage against an externally-anchored window.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: D4 frame; the warhead library; the verified covalent anchors
Output: per-class LUMO, a window verdict per class, and D4 with reactivity columns

WHAT THIS IS AND IS NOT FOR (D0005). The reactivity window is a SAFETY filter
for condition (ii) — is this electrophile inside the range spanned by real,
wet-lab-validated covalent Cys113 actives? It is NOT a potency signal, and T_4
must never rank on it. Across Sulfopin and the Reddi 4a-4g series, measured
intrinsic reactivity and Pin1 labelling correlate at only r = 0.396 over a 13.6x
range of k: 4e has nearly the lowest reactivity and the highest labelling. Within
the precedented range, electrophilicity does not predict engagement.

REACTIVITY IS A WARHEAD PROPERTY, NOT AN R-GROUP PROPERTY. So LUMO is computed
ONCE PER CLASS on a fixed reference R-group, not per candidate — 9 calculations
rather than 1,782. That is also why this stage is cheap enough to run BEFORE
covalent docking even though Rev 3 section 6 lists it after: section 9's own
funnel puts triage ahead of docking, and spending GPU-hours docking a class that
the window will reject afterwards is pure waste.

THE WINDOW IS CHEMOTYPE-SKEWED, AND THAT IS CARRIED FORWARD. The clean selective
anchors are chloroacetamides and the extra kinetic anchors are promiscuous
quinones, so the window is chloroacetamide-centric with a reactive-quinone upper
tail. It is usable, not chemotype-balanced, and no single global cutoff should be
over-trusted.
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from shared import io as dio                     # noqa: E402
from shared import reference_set as rs           # noqa: E402
from shared import warhead_library as wl         # noqa: E402
from shared.manifest import Manifest             # noqa: E402

log = logging.getLogger("t4-reactivity")

EXPERIMENT = "04_t4_combinatorial"
XTB = Path("/data/lab_vm/envs/dwi_cheminf/bin/xtb")
OUT_DIR = Path("/data/lab_vm/append_only/inhibition") / EXPERIMENT / "reactivity"

# The fixed reference R-group every warhead is scored on. Methyl: the smallest
# group that leaves the electrophile chemically intact, so differences between
# classes are the warhead's and not the decoration's.
REFERENCE_RGROUP = "C"

# How far outside the anchors' observed LUMO span a class may sit. The window is
# anchored on real actives, so this is a tolerance around measured chemistry
# rather than an arbitrary cutoff.
WINDOW_TOLERANCE_EV = 0.5


class ReactivityError(RuntimeError):
    """LUMO could not be computed."""


def build_reference_molecule(warhead_fragment: str) -> str | None:
    """Attach the warhead to the reference R-group, giving a scoreable molecule."""
    smiles = warhead_fragment.replace("[*]", REFERENCE_RGROUP)
    m = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(m) if m is not None else None


def lumo_ev(smiles: str, *, charge: int = 0) -> dict | None:
    """Single-point GFN2-xTB HOMO/LUMO for one molecule.

    Returns
    -------
    dict or None
        ``{"homo_ev", "lumo_ev", "gap_ev"}``, or None if the calculation failed.
    """
    if not XTB.is_file():
        raise ReactivityError(f"xtb not found at {XTB} — activate the cheminf env")
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    m = Chem.AddHs(m)
    if AllChem.EmbedMolecule(m, randomSeed=42) != 0:
        return None
    try:
        AllChem.MMFFOptimizeMolecule(m, maxIters=500)
    except Exception:  # noqa: BLE001
        pass

    with tempfile.TemporaryDirectory() as td:
        xyz = Path(td) / "mol.xyz"
        Chem.MolToXYZFile(m, str(xyz))
        proc = subprocess.run(
            [str(XTB), xyz.name, "--gfn", "2", "--sp", "--chrg", str(charge)],
            cwd=td, capture_output=True, text=True, timeout=600)
        if proc.returncode != 0:
            log.warning("xtb failed for %s: %s", smiles, proc.stderr[-200:])
            return None
        homo = lumo = None
        for line in proc.stdout.splitlines():
            # xtb marks the frontier orbitals in its orbital listing
            if "(HOMO)" in line:
                parts = line.split()
                homo = float(parts[-2]) if len(parts) >= 2 else None
            elif "(LUMO)" in line:
                parts = line.split()
                lumo = float(parts[-2]) if len(parts) >= 2 else None
        if homo is None or lumo is None:
            return None
        return {"homo_ev": homo, "lumo_ev": lumo, "gap_ev": lumo - homo}


def anchor_window() -> dict:
    """LUMO span of the VERIFIED covalent anchors — the window's bounds.

    Control B5: the window is bounded by real wet-lab-validated actives, and the
    project's own leads are excluded as anchors. `covalent_anchors` already
    refuses UNVERIFIED rows.
    """
    anchors = rs.covalent_anchors()
    results = {}
    for _, a in anchors.iterrows():
        r = lumo_ev(str(a["canonical_smiles"]))
        if r:
            results[str(a["name"])] = r["lumo_ev"]
            log.info("anchor %-28s LUMO %+.3f eV", a["name"], r["lumo_ev"])
        else:
            log.warning("anchor %s: LUMO failed", a["name"])
    if len(results) < 2:
        raise ReactivityError(
            f"only {len(results)} anchor LUMO(s) computed; a window bounded by "
            "fewer than two real actives is the n~1 problem control B5 exists to "
            "prevent.")
    lo, hi = min(results.values()), max(results.values())
    return {"anchors": results, "lumo_min": lo, "lumo_max": hi,
            "window_lo": lo - WINDOW_TOLERANCE_EV,
            "window_hi": hi + WINDOW_TOLERANCE_EV,
            "caveat": rs.window_caveat()}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--frame", default=None, help="D4 parquet (default: latest)")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    d = Path("/data/lab_vm/append_only/inhibition") / EXPERIMENT
    frame_path = Path(args.frame) if args.frame else dio.latest(d, "D4", ".parquet")
    if frame_path is None:
        raise SystemExit("no D4 frame found — run 01_enumerate.py first")
    df = pd.read_parquet(frame_path)
    log.info("loaded %s (%d candidates)", frame_path.name, len(df))

    log.info("computing the anchor window from VERIFIED covalent actives")
    window = anchor_window()
    log.info("anchor LUMO span %.3f .. %.3f eV; window %.3f .. %.3f eV",
             window["lumo_min"], window["lumo_max"],
             window["window_lo"], window["window_hi"])

    log.info("computing per-class LUMO on the fixed reference R-group (%s)",
             REFERENCE_RGROUP)
    lib = wl.load()
    per_class = {}
    for _, w in lib.iterrows():
        cls = str(w["class_id"])
        if cls not in set(df["warhead_class"]):
            continue
        ref = build_reference_molecule(str(w["warhead_fragment_smiles"]))
        if ref is None:
            log.warning("class %s: reference molecule could not be built", cls)
            continue
        r = lumo_ev(ref)
        if r is None:
            log.warning("class %s: LUMO failed", cls)
            continue
        inside = window["window_lo"] <= r["lumo_ev"] <= window["window_hi"]
        per_class[cls] = {**r, "reference_smiles": ref, "in_window": inside}
        log.info("%-22s LUMO %+.3f eV  gap %.3f  in_window=%s",
                 cls, r["lumo_ev"], r["gap_ev"], inside)

    df["lumo_ev"] = df["warhead_class"].map(
        lambda c: per_class.get(c, {}).get("lumo_ev"))
    df["homo_lumo_gap_ev"] = df["warhead_class"].map(
        lambda c: per_class.get(c, {}).get("gap_ev"))
    df["reactivity_in_window"] = df["warhead_class"].map(
        lambda c: per_class.get(c, {}).get("in_window"))

    # Stamp only rows not already stamped — a candidate keeps the reason it was
    # FIRST set aside. And note this stamps, it does not delete.
    stamp = df["reactivity_in_window"].eq(False) & df["rejected_at"].isna()
    df.loc[stamp, "rejected_at"] = "reactivity_window"

    out_json = OUT_DIR / "reactivity_window.json"
    out_json.write_text(json.dumps({"window": window, "per_class": per_class},
                                   indent=2) + "\n", encoding="utf-8")

    out = dio.write_full_frame(
        df, approach="t4", experiment=EXPERIMENT, stage="t4_reactivity_triage",
        params={"reference_rgroup": REFERENCE_RGROUP,
                "window_tolerance_ev": WINDOW_TOLERANCE_EV,
                "window_lo": window["window_lo"], "window_hi": window["window_hi"]},
        inputs={"d4_frame": frame_path})

    (Manifest(stage="t4_reactivity_window", approach="t4")
     .add_output("window", out_json)
     .note(window["caveat"])
     .write(OUT_DIR, filename="reactivity_manifest.json"))

    print(f"\nT_4 reactivity triage -> {out}")
    print(f"  anchor window   {window['window_lo']:+.3f} .. {window['window_hi']:+.3f} eV")
    print(f"  anchors used    {len(window['anchors'])}: "
          f"{', '.join(window['anchors'])}")
    print("\n  per warhead class:")
    for cls, v in sorted(per_class.items(), key=lambda kv: kv[1]["lumo_ev"]):
        mark = "in " if v["in_window"] else "OUT"
        print(f"    {mark} {cls:22s} LUMO {v['lumo_ev']:+.3f} eV")
    outside = [c for c, v in per_class.items() if not v["in_window"]]
    print(f"\n  classes outside the window: {outside or 'none'}")
    print(f"  candidates surviving        : {int(df['rejected_at'].isna().sum())}")
    print(f"\n  {window['caveat'][:200]}")


if __name__ == "__main__":
    main()
