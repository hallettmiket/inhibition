#!/usr/bin/env python3
"""
Purpose: audit pose generation -- are low-energy poses actually in the pocket?
Author: Timothy Wu (with Claude Code)
Date: 2026-08-26
Input: docks molecules fresh, keeping the .dlg so per-pose energies survive
Output: 00_outputs/blacksmith/pose_generation_audit/

@tt8804, looking at the pose viewer: "there are literally poses outside of the
pocket, how is that possibly the lowest energy".

THE PERSISTED CLOUDS CARRY NO ENERGIES. `rebuild_and_match` returns conformers
and nothing else, `persist_raw_clouds` writes them with no SD tags, and every
clustering experiment so far -- exp/16, 17, 19, 20 -- therefore treated the
best-scoring pose and the 500th identically. So did the viewer. That is the gap
this closes: the question "how is that the lowest energy" could not be asked of
any artefact the project had written.

WHAT THIS TESTS, IN ORDER OF WHAT WOULD BE WORST:
  1. Is the ENERGY-GEOMETRY relationship inverted or absent? If the least buried
     poses score best, the scoring function is broken on this target and every
     downstream stage inherits it.
  2. Is the BOX wrong? `config/receptor.yaml` gives T_4 a 20 A covalent box;
     `nac_screen.py:130` hardcodes 26. Both are populated and plausible -- the
     receptor-disagreement shape CLAUDE.md already warns about.
  3. Does PoseBusters reject the poses that look wrong, and do those poses carry
     competitive energies?

BURIAL IS MEASURED THREE WAYS because no single one is decisive: contacts (how
many receptor atoms are near), exposure (what fraction of the LIGAND touches
nothing), and enclosure (how much protein surrounds each ligand atom). A ligand
lying in a shallow surface groove scores differently on each, and "outside the
pocket" has to mean something specific before it can be checked.
"""

from __future__ import annotations

import argparse
import logging
import os
import re
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

import nac_screen as ns                             # noqa: E402
import nac_rank as nr                               # noqa: E402
from shared import outputs as sout                  # noqa: E402
from shared import run_paths as rp                  # noqa: E402

log = logging.getLogger("pose-audit")


def receptor_atoms():
    """(heavy-atom coords, Cys113 SG) from the receptor the SCREEN docks into."""
    rec, sg = [], None
    for ln in rp.receptor_prep().read_text().splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        if ((ln[76:78].strip() or ln[12:16].strip()[:1]).upper()) == "H":
            continue
        x = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
        rec.append(x)
        if ln[22:26].strip() == "113" and ln[12:16].strip() == "SG":
            sg = np.array(x)
    return np.array(rec), sg


def _pose_signature(P: np.ndarray, ref: np.ndarray) -> np.ndarray:
    """Sorted distances from each ligand heavy atom to a fixed reference point.

    An identity for a POSE that survives the two things that broke the first two
    attempts at this check:
      * atom ORDER differs between the PDBQT record and the rebuilt RDKit mol,
        so any element-wise comparison comes out at chance (measured: best match
        1.839 A against a runner-up of 1.869 A -- no discrimination at all);
      * sorting makes it order-invariant, and taking distances to an EXTERNAL
        point keeps it sensitive to rigid-body placement, which an internal
        distance signature would not be -- two poses that differ only by a
        rotation have identical internal distances.
    """
    return np.sort(np.sqrt(((P - ref) ** 2).sum(-1)))


def _record_ligand_heavy(rec_str: str) -> np.ndarray:
    """Ligand heavy atoms from one DLG pose record, EXCLUDING the flexible residue.

    The docking is run with `--flexres`, so every record holds the ligand AND the
    movable Cys113 side chain -- 19 ligand atoms plus 4 protein atoms on the
    molecule measured here. Counting all 23 was what made the first version of
    this guard fire on all six molecules: the centroids disagreed by ~1 A for a
    reason that had nothing to do with pose order. The residue is delimited by
    BEGIN_RES/END_RES, so it is excluded structurally rather than by count.
    """
    xs, in_res = [], False
    for ln in rec_str.splitlines():
        if ln.startswith("BEGIN_RES"):
            in_res = True
        elif ln.startswith("END_RES"):
            in_res = False
        elif ln.startswith(("ATOM", "HETATM")) and not in_res:
            t = (ln[77:79].strip() or ln[12:16].strip()[:1] or "X")
            if not t.upper().startswith("H"):
                xs.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    return np.array(xs)


def energies_aligned(dlg: Path, mol, ref: np.ndarray) -> np.ndarray:
    """Per-pose free energy, PROVEN to line up with `mol`'s conformers.

    Energies and geometry come from two different objects and AutoDock reports a
    cluster ranking beside the run order, so "the order" is genuinely ambiguous
    and pairing by position is the shape this project fails on. The pairing is
    therefore SOLVED, not assumed: every record is matched to every conformer by
    pose signature under a Hungarian assignment, and the result must come back as
    the identity permutation at essentially zero error. Measured on this build it
    does -- 0.00000 A median. Anything else raises.
    """
    from meeko import PDBQTMolecule
    from scipy.optimize import linear_sum_assignment
    pm = PDBQTMolecule.from_file(str(dlg), is_dlg=True, skip_typing=True)
    fe = np.array(pm._pose_data["free_energies"], dtype=float)
    strings = list(pm._pose_data["pdbqt_string"])
    n = mol.GetNumConformers()
    if not (len(fe) == len(strings) == n):
        raise ValueError(f"{len(fe)} energies, {len(strings)} records, {n} conformers")
    hv = [k for k in range(mol.GetNumAtoms())
          if mol.GetAtomWithIdx(k).GetAtomicNum() > 1]
    SR = np.array([_pose_signature(_record_ligand_heavy(s), ref) for s in strings])
    SC = np.array([_pose_signature(mol.GetConformer(c).GetPositions()[hv], ref)
                   for c in range(n)])
    if SR.shape[1] != SC.shape[1]:
        raise ValueError(f"record has {SR.shape[1]} ligand heavy atoms, "
                         f"conformer has {SC.shape[1]}")
    D = np.sqrt(((SR[:, None, :] - SC[None, :, :]) ** 2).sum(-1))
    r, c = linear_sum_assignment(D)
    if not (c == np.arange(n)).all():
        bad = int((c != np.arange(n)).sum())
        raise ValueError(f"{bad} of {n} records do not map to their own conformer; "
                         "the energies are NOT in conformer order")
    err = float(np.median(D[r, c]))
    if err > 0.01:
        raise ValueError(f"pose matching is not exact (median {err:.4f} A)")
    return fe


def burial(xyz: np.ndarray, rec: np.ndarray, sg: np.ndarray) -> pd.DataFrame:
    out = []
    for i, p in enumerate(xyz):
        dm = np.sqrt(((p[:, None, :] - rec[None, :, :]) ** 2).sum(-1))
        near = (dm < 4.5).sum(1)
        out.append(dict(pose=i,
                        contacts=int((dm.min(0) < 4.5).sum()),
                        exposed_frac=float((near == 0).mean()),
                        enclosure=float((dm < 8.0).sum(1).mean()),
                        sg_dist=float(np.sqrt(((p - sg) ** 2).sum(-1)).min())))
    return pd.DataFrame(out)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--molecules", default="", help="comma-separated idents")
    ap.add_argument("--n-molecules", type=int, default=6)
    ap.add_argument("--n-runs", type=int, default=500)
    ap.add_argument("--gpu", default="2")
    ap.add_argument("--posebusters", action="store_true", default=True)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        datefmt="%H:%M:%S")
    from scipy.stats import spearmanr

    rec, sg = receptor_atoms()
    log.info("receptor: %d heavy atoms; Cys113 SG at %s", len(rec), sg is not None)

    cands = {c.ident: c for c in nr.load_candidates()}
    if a.molecules:
        want = [x.strip() for x in a.molecules.split(",") if x.strip()]
    else:
        have = [d.name[len("raw_cloud_"):] for d in rp.BLACKSMITH.glob("raw_cloud_*")]
        want = [m for m in have if m in cands][:a.n_molecules]
    log.info("auditing %d molecules", len(want))

    rec_dir = ns.build_reactive_receptor(ns.RX_RECEPTOR)
    allrows = []
    for n, ident in enumerate(want, 1):
        work = Path(tempfile.mkdtemp(prefix=f"audit_{ident[:10]}_"))
        try:
            ligs = list(ns.prepare_ligand(cands[ident], work / "lig.pdbqt"))
            dlg = ns.dock(ligs[0], rec_dir, work / "d", a.n_runs, a.gpu, seed=a.seed)
            mol, match = ns.rebuild_and_match(dlg, cands[ident])
            e = energies_aligned(dlg, mol, sg)
            hv = [x.GetIdx() for x in mol.GetAtoms() if x.GetAtomicNum() > 1]
            xyz = np.array([mol.GetConformer(c).GetPositions()[hv]
                            for c in range(mol.GetNumConformers())])
            b = burial(xyz, rec, sg)
            b["energy"] = e
            b["ident"] = ident
            b["energy_rank"] = b.energy.rank(method="first").astype(int)
            if a.posebusters:
                try:
                    from posebusters import PoseBusters
                    from rdkit import Chem
                    tmp = work / "poses.sdf"
                    w = Chem.SDWriter(str(tmp))
                    for c in range(mol.GetNumConformers()):
                        w.write(mol, confId=c)
                    w.close()
                    pbres = PoseBusters(config="mol").bust([tmp], None, None)
                    ok = pbres.all(axis=1).values
                    if len(ok) == len(b):
                        b["pb_valid"] = ok
                except Exception as exc:                       # noqa: BLE001
                    log.warning("  PoseBusters skipped: %s", str(exc)[:70])
            allrows.append(b)
            log.info("  %d/%d %s: %d poses, energy %.2f to %.2f kcal/mol",
                     n, len(want), ident, len(b), b.energy.min(), b.energy.max())
        except Exception as exc:                                # noqa: BLE001
            log.error("  %s FAILED: %s", ident, str(exc)[:140])
        finally:
            import shutil
            shutil.rmtree(work, ignore_errors=True)

    if not allrows:
        raise SystemExit("nothing audited")
    d = pd.concat(allrows, ignore_index=True)
    t = sout.Topic("blacksmith", "pose_generation_audit")
    d.to_csv(t.write("per_pose", ".csv"), index=False)

    P = print
    P("\n" + "=" * 82)
    P("  POSE GENERATION AUDIT — do the low-energy poses sit in the pocket?")
    P("=" * 82)
    P(f"\n  {d.ident.nunique()} molecules, {len(d):,} poses, energies recovered from "
      f"the .dlg\n")

    P("  1. THE ENERGY-GEOMETRY RELATIONSHIP (per molecule, then pooled)\n")
    P(f"    {'quantity':<22}{'rho vs energy':>15}   interpretation")
    for col, good in (("contacts", "negative = more contacts score better"),
                      ("enclosure", "negative = more buried scores better"),
                      ("exposed_frac", "POSITIVE = more exposed scores WORSE"),
                      ("sg_dist", "positive = closer to Cys113 scores better")):
        rr = [spearmanr(g[col], g.energy)[0] for _, g in d.groupby("ident")]
        P(f"    {col:<22}{np.median(rr):+15.3f}   {good}")
    P("\n    (energy is negative-is-better, so a NEGATIVE rho with `contacts`")
    P("     means better-scoring poses are more buried -- the expected sign.)")

    P("\n  2. ARE THE BEST-SCORING POSES THE BURIED ONES?\n")
    P(f"    {'energy decile':>14}{'contacts':>10}{'exposed':>10}{'enclosure':>11}{'SG dist':>9}")
    d["dec"] = d.groupby("ident").energy.transform(
        lambda s: pd.qcut(s.rank(method="first"), 10, labels=False))
    for k, g in d.groupby("dec"):
        tag = "best 10%" if k == 0 else ("worst 10%" if k == 9 else f"{int(k) + 1}")
        P(f"    {tag:>14}{g.contacts.median():10.0f}{g.exposed_frac.median() * 100:9.0f}%"
          f"{g.enclosure.median():11.0f}{g.sg_dist.median():8.1f}A")

    P("\n  3. THE ENERGY SPREAD — how much does the score separate poses?\n")
    for ident, g in list(d.groupby("ident"))[:6]:
        P(f"    {ident:<20} best {g.energy.min():7.2f}  median {g.energy.median():7.2f}"
          f"  worst {g.energy.max():7.2f}   span {g.energy.max() - g.energy.min():5.2f} kcal/mol")
    P(f"\n    poses within 1 kcal/mol of the best: "
      f"{d.groupby('ident').apply(lambda g: (g.energy <= g.energy.min() + 1).mean()).median() * 100:.0f}% (median molecule)")
    P(f"    poses within 2 kcal/mol of the best: "
      f"{d.groupby('ident').apply(lambda g: (g.energy <= g.energy.min() + 2).mean()).median() * 100:.0f}%")

    P("\n  4. THE WORST-LOOKING POSES — are they scored badly?\n")
    bad = d[d.exposed_frac > 0.3]
    if len(bad):
        P(f"    poses with >30% of atoms uncontacted: {len(bad):,} ({len(bad) / len(d) * 100:.1f}%)")
        P(f"    their median energy percentile within their molecule: "
          f"{d.assign(pct=d.groupby('ident').energy.rank(pct=True))[d.exposed_frac > 0.3].pct.median() * 100:.0f}%")
        P(f"    how many are in their molecule's BEST decile: "
          f"{(bad.dec == 0).sum()} of {len(bad)}")
    else:
        P("    none")

    if "pb_valid" in d:
        P("\n  5. POSEBUSTERS\n")
        P(f"    valid: {d.pb_valid.mean() * 100:.1f}%")
        P(f"    invalid poses in the best energy decile: "
          f"{int(((~d.pb_valid) & (d.dec == 0)).sum())}")
        P(f"    median energy percentile of INVALID poses: "
          f"{d.assign(pct=d.groupby('ident').energy.rank(pct=True))[~d.pb_valid].pct.median() * 100:.0f}%")

    P("\n" + "=" * 82)
    P(f"  written to {t.dir}")
    P("=" * 82 + "\n")


if __name__ == "__main__":
    main()
