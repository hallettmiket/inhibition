"""
Purpose: dock the co-folding benchmark's 15 ligands ourselves, so T1 has a comparator.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-07
Input: the set built by cofold_prep_set.py (truth PDB + SMILES per entry)
Output: 00_outputs/blacksmith/cofold/cofold_dock_<N>.csv

`docs/prereg_cofolding.md` fixes this in advance: *"Our docking is measured on
the same 10 molecules, not compared against the 82-case benchmark. Different case
sets are not comparable."* This script produces that number.

IT TURNS OUT NONE OF THE TEN WERE EVER DOCKED BY US. Four of the held-out
entries (9V6G/9V6I/9V6P/9V6W) do appear in the 82-case benchmark, which looked
at first like a head start. They do not help: every one of those cases docked
**A1ERA**, a second ligand present in all four entries, not the covalent
inhibitor that makes the entry interesting. The overlap with the ten molecules
under test is zero, so all fifteen are docked here from scratch.

THE COMPARISON IS DELIBERATELY LIKE-FOR-LIKE, AND THAT COSTS US THE HOME FIELD.
Both methods are denied the holo receptor of the complex they are predicting.
Boltz builds a structure from sequence; we dock into 3IKD_ian, a DIFFERENT Pin1
crystal. Docking into each entry's own receptor would measure self-docking, a
much easier problem, and would flatter our side of a comparison we are running to
inform a real decision.

EVERY MEASUREMENT FUNCTION IS IMPORTED, NOT REIMPLEMENTED. The superposition,
the pose parsing and the assignment RMSD all come from
`redock_02_build_cases` / `redock_3ikd_benchmark`. A second implementation of
"recovered the crystal pose" is how two different answers to one question end up
in the same report.

RESOLUTION IS A CAVEAT ON THE GROUND TRUTH, NOT ON EITHER METHOD. The held-out
entries run 1.99-2.96 A, against a median of 1.81 A across the 82-case set --
they are newer and harder. That blurs the target both methods aim at, equally.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import logging
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import outputs as sout                  # noqa: E402

log = logging.getLogger("cofold-dock")
RB = Path("/data/lab_vm/append_only/inhibition/05_redock_benchmark")
WORK = Path("/data/lab_vm/modifiable/inhibition/cofold/docking")
OUT = sout.Topic("blacksmith", "cofold")
SUCCESS_A = 2.0


def _load(name: str):
    spec = importlib.util.spec_from_file_location(
        name, REPO / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def prepare_ligand(smiles: str, out: Path) -> bool:
    """SMILES -> 3D -> pH 7.4 -> pdbqt, the production reactive-path preparation.

    Protonation at pH 7.4 is D0074: the species docked is the one that exists at
    physiological pH, not the one the SMILES happens to be drawn as.
    """
    obabel = "/data/lab_vm/envs/dwi_cheminf/bin/obabel"
    out.parent.mkdir(parents=True, exist_ok=True)
    smi = out.with_suffix(".smi")
    smi.write_text(smiles + "\n")
    r = subprocess.run(
        [obabel, str(smi), "-O", str(out), "--gen3d", "-p", "7.4"],
        capture_output=True, text=True, timeout=600)
    return out.is_file() and out.stat().st_size > 0 and "0 molecules" not in r.stderr


def build_ref(rd02, truth: Path, receptor_pdbqt: Path, out: Path) -> dict:
    """Put this entry's crystal ligand into 3IKD's coordinate frame."""
    ref_atoms, _ = rd02.read_pdb(receptor_pdbqt)
    ref_ca = rd02.ca_map(ref_atoms)
    atoms, _ = rd02.load_structure(truth)
    chains = {a["chain"] for a in atoms if a["record"] == "ATOM"}
    best = None
    for ch in sorted(chains):
        mine = rd02.ca_map(atoms, ch)
        shared = sorted(set(mine) & set(ref_ca))
        if len(shared) < 50:
            continue
        rot, trans, fit = rd02.kabsch(
            np.array([mine[r][1] for r in shared]),
            np.array([ref_ca[r][1] for r in shared]))
        if best is None or fit < best[2]:
            best = (rot, trans, fit, len(shared), ch)
    if best is None:
        return {"status": "no chain fit"}
    rot, trans, fit, n_fit, ch = best
    lig = [a for a in atoms if a["record"] == "HETATM"
           and a["resname"] not in ("HOH", "WAT")]
    if not lig:
        return {"status": "no ligand in truth"}
    xyz = np.array([a["xyz"] for a in lig])
    rd02.write_ligand_pdb(lig, out, xyz @ rot + trans)
    return {"status": "ok", "ca_fit_rmsd": round(fit, 3), "n_ca": n_fit,
            "chain": ch, "n_lig_atoms": len(lig)}


def elems_pdbqt(p: Path) -> list[tuple[np.ndarray, np.ndarray]]:
    """Every docked mode as (heavy-atom coords, elements)."""
    out, cur_x, cur_e = [], [], []
    for ln in p.read_text(errors="replace").splitlines():
        if ln.startswith("MODEL"):
            cur_x, cur_e = [], []
        elif ln.startswith("ENDMDL"):
            if cur_x:
                out.append((np.array(cur_x), np.array(cur_e)))
            cur_x, cur_e = [], []
        elif ln.startswith(("ATOM", "HETATM")):
            # PDBQT columns 78-79 hold an AUTODOCK TYPE, not an element: 'A' is
            # aromatic carbon, 'NA'/'OA' are H-bond-accepting N/O. Mapping those
            # to elements is required before any element-aware comparison, and
            # getting it wrong would silently forbid correct atom pairings.
            t = (ln[77:79].strip() or ln[12:16].strip()[:1]).upper()
            el = {"A": "C", "NA": "N", "OA": "O", "SA": "S",
                  "HD": "H", "N": "N", "C": "C", "O": "O", "S": "S"}.get(t, t[:1])
            if el == "H":
                continue
            try:
                cur_x.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
                cur_e.append(el)
            except ValueError:
                pass
    if cur_x:
        out.append((np.array(cur_x), np.array(cur_e)))
    return out


def elems_pdb(p: Path) -> tuple[np.ndarray, np.ndarray]:
    xs, es = [], []
    for ln in p.read_text(errors="replace").splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        el = (ln[76:78].strip() or ln[12:16].strip()[:1]).upper()
        if el == "H":
            continue
        xs.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
        es.append(el)
    return np.array(xs), np.array(es)


def score_matched(cb, refs: Path, poses: Path) -> pd.DataFrame:
    """Re-score every docked pose with the EXACT metric used on Boltz's poses.

    `redock_3ikd_benchmark.rmsd` truncates the longer atom list by index when the
    counts differ. Across the 82-case set the counts match, so this never bit.
    Here they NEVER match: the crystal ligand is the adduct and has lost its
    leaving group, while the docked ligand is the intact compound. Truncating by
    file order would drop whichever atom happened to be written last.

    The two methods must be judged by one ruler or the comparison means nothing,
    so both go through `cofold_bench.sym_rmsd` — rectangular, element-aware
    optimal assignment.
    """
    rows = []
    for ref in sorted(refs.glob("*_ref.pdb")):
        case = ref.name.replace("_ref.pdb", "")
        pose = poses / f"{case}_out.pdbqt"
        if not pose.is_file():
            continue
        rx, re_ = elems_pdb(ref)
        modes = elems_pdbqt(pose)
        if not len(rx) or not modes:
            continue
        vals = [cb.sym_rmsd(mx, me, rx, re_) for mx, me in modes]
        rows.append({"ident": case, "n_modes_m": len(vals),
                     "top1_rmsd_m": vals[0], "best_rmsd_m": float(np.nanmin(vals)),
                     "best_rank_m": int(np.nanargmin(vals)) + 1,
                     "n_atoms_ref": len(rx), "n_atoms_docked": len(modes[0][0])})
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--set", required=True)
    ap.add_argument("--gpu", type=int, default=3)
    ap.add_argument("--skip-dock", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    rb = _load("redock_3ikd_benchmark")
    rd02 = rb._rd02()
    rec3 = sout.Topic("blacksmith", "receptor_3ikd")
    receptor = rb.latest(rec3, "3IKD_prepared", ".pdbqt")
    box = json.loads(rb.latest(rec3, "box_3IKD", ".json").read_text())
    log.info("receptor %s — the SAME receptor production docks into", receptor.name)

    df = pd.read_csv(args.set)
    refs, ligs = WORK / "refs", WORK / "ligands"
    refs.mkdir(parents=True, exist_ok=True)
    ligs.mkdir(parents=True, exist_ok=True)

    rows = []
    for r in df.itertuples():
        rec = {"ident": r.ident, "held_out": r.held_out, "comp_id": r.comp_id}
        rec.update(build_ref(rd02, Path(r.truth), receptor,
                             refs / f"{r.ident}_ref.pdb"))
        rec["prepared"] = prepare_ligand(r.smiles, ligs / f"{r.ident}.pdbqt")
        rows.append(rec)
        log.info("%s  ref=%s  ligand=%s", r.ident, rec.get("status"),
                 "ok" if rec["prepared"] else "FAILED")

    prep = pd.DataFrame(rows)
    usable = prep[(prep.status == "ok") & prep.prepared]
    log.info("%d/%d entries usable (reference transformed AND ligand prepared)",
             len(usable), len(prep))

    poses = WORK / "poses"
    if not args.skip_dock:
        rb.dock(ligs, poses, receptor, box, args.gpu, timeout_s=3600)

    cb = _load("cofold_bench")
    sc = rb.score(refs, poses, "3IKD").rename(columns={"case": "ident"})
    out = (prep.merge(sc, on="ident", how="left")
                .merge(score_matched(cb, refs, poses), on="ident", how="left"))
    dest = OUT.write("cofold_dock", ".csv")
    out.to_csv(dest, index=False)

    print("\nOur docking into 3IKD, on the co-folding benchmark's own ligands")
    print(f"  -> {dest}\n")
    print("  Scored with the SAME element-aware assignment RMSD as Boltz's poses")
    print("  (`_m` columns). The benchmark's own metric is kept alongside for")
    print("  continuity with the 82-case number, but is NOT the comparator.\n")
    for grp, g in out.dropna(subset=["best_rmsd_m"]).groupby("held_out"):
        lab = "HELD OUT (T1 comparator)" if grp else "in training"
        top1 = (g.top1_rmsd_m <= SUCCESS_A).mean() * 100
        best = (g.best_rmsd_m <= SUCCESS_A).mean() * 100
        print(f"  {lab:<26} n={len(g):<3} top-1 {top1:5.1f}%   "
              f"best-of-N {best:5.1f}%   median best {g.best_rmsd_m.median():.2f} A")
    print("\n  Same receptor, same box, same 2.0 A bar as the 82-case benchmark —")
    print("  but measured on THESE molecules, which is the only comparison the")
    print("  pre-registration allows.")


if __name__ == "__main__":
    main()
