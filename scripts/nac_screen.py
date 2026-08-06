"""
Purpose: does the viable-NAC fraction discriminate between molecules? The go/no-go for the ranking framework.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-05
Input: a set of warhead-bearing molecules + the reactive 3IKD receptor
Output: 00_outputs/blacksmith/nac_screen/nac_screen_<N>.csv + a report

THE EXPERIMENT THE WHOLE FRAMEWORK RESTS ON.

`docs/ranking_rationale.md` proposes ranking on whether a molecule can orient to
present its warhead. `shared/nac_criterion.py` turns that into a number: the
fraction of independent docking runs that reach a mechanism-appropriate
near-attack conformation. This script asks the only question that matters before
anything is built on top of it:

    DOES THAT NUMBER VARY BETWEEN MOLECULES?

If every molecule scores alike, there is no ranking, and the honest outcome is to
say so and stop -- that is the rationale's second stated failure mode, "the
near-attack gate does not discriminate". This runs before the expensive stages
for the same reason `pose_selection_bench` did.

THE CONTROL, AND ITS LIMITS. Negatives are drawn from AID 504891 molecules
MEASURED inactive, filtered to those carrying the SAME warhead SMARTS as the
positives. That matching is the point: a molecule with no warhead fails the NAC
test trivially for lack of a reactive atom, so an unmatched control would measure
nothing but the presence of a warhead (D0014 makes the same argument for covalent
decoys).

    These labels are NOT authoritative and are not treated as such. Single-
    concentration HTS inactivity is weak evidence -- compounds read inactive for
    solubility, aggregation, or assay interference as readily as for failing to
    bind. The primary readout here is therefore SPREAD, which needs no labels at
    all; the positive/negative comparison is secondary and reported with that
    caveat attached.

WHY SPREAD IS THE PRIMARY READOUT. Spread is a property of the measurement, not
of anyone's annotation. A flat distribution kills the framework no matter how
good the labels are, and a wide one earns the right to go and find better ones.
"""

from __future__ import annotations

import argparse
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import nac_criterion as nac       # noqa: E402
from shared import outputs as sout            # noqa: E402

log = logging.getLogger("nac-screen")

ENV = Path("/data/lab_vm/envs" if Path("/data/lab_vm/envs").is_dir() else "")
RX_ENV = Path.home() / ".micromamba/envs/dwi_reactive"
AUTODOCK = Path("/data/lab_vm/modifiable/inhibition/autodock_gpu/bin/autodock_gpu_64wi")
RECEPTOR_PDB = Path("/data/lab_vm/modifiable/inhibition/receptor_3ikd_prep/3IKD_noligand.pdb")
RX_RECEPTOR = Path("/data/lab_vm/modifiable/inhibition/receptor_3ikd_reactive")

OUT = sout.Topic("blacksmith", "nac_screen")

# The NAC parameterisation, justified in shared/nac_criterion.py: a near-attack
# conformation is a van der Waals CONTACT, so the well sits at contact distance
# and the neighbouring atoms keep their real steric radii. The published
# covalent-distance defaults (1.8 A, radii scaled 0.5x) produced 0/40 chemically
# viable poses.
R_EQ_12 = 3.2
EPS_12 = 1.0
SCALING = 1.0


@dataclass(frozen=True)
class Candidate:
    ident: str
    smiles: str
    warhead_class: str
    mechanism: str
    reactive_smarts: str
    label: str            # "positive" | "negative"


# --------------------------------------------------------------------------

def largest_fragment(smiles: str):
    """Parse a SMILES and keep its largest fragment.

    HTS tables carry salts and co-formers ("...CCl.Cl", sodium counterions).
    Docking a counterion alongside the ligand is meaningless, and meeko refuses
    a multi-fragment molecule outright, so a salt would otherwise be dropped as
    UNMEASURABLE — which silently biases the negative set toward whichever
    compounds happen to be supplied as free bases.
    """
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    return max(frags, key=lambda f: f.GetNumHeavyAtoms()) if len(frags) > 1 else mol


def _atom_types(pdbqt: str) -> list[str]:
    """The AutoDock type column, in file order."""
    return [ln[77:79].strip() for ln in pdbqt.splitlines()
            if ln.startswith(("ATOM", "HETATM"))]


def build_reactive_receptor(dest: Path) -> Path:
    """Prepare 3IKD with Cys113 SG as the reactive atom. Cached — it is deterministic."""
    cfg = dest / "rec.reactive_config"
    if cfg.is_file() and (dest / "rec_rigid.maps.fld").is_file():
        log.info("reactive receptor already built at %s", dest)
        return dest
    dest.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [str(RX_ENV / "bin/mk_prepare_receptor.py"), "--read_pdb", str(RECEPTOR_PDB),
         "-o", "rec", "-p", "-g", "-r", "A:113", "--reactive_name", "CYS:SG",
         "--box_center_off_reactive_res", "--box_size", "26", "26", "26",
         "--r_eq_12", str(R_EQ_12), "--eps_12", str(EPS_12),
         "--r_eq_13_scaling", str(SCALING), "--r_eq_14_scaling", str(SCALING)],
        cwd=dest, check=True, capture_output=True, text=True)
    subprocess.run([str(RX_ENV / "bin/autogrid4"), "-p", "rec_rigid.gpf",
                    "-l", "rec_rigid.glg"], cwd=dest, check=True, capture_output=True)
    log.info("built reactive receptor -> %s", dest)
    return dest


def prepare_ligand(cand: Candidate, path: Path) -> list[Path]:
    """Write one reactive PDBQT PER REACTIVE CENTRE, with the warhead atom retyped.

    Raises on failure rather than writing an untyped ligand: a ligand silently
    prepared WITHOUT reactive types docks fine and lands nowhere near the
    sulfur, which would read as "this molecule cannot reach a NAC" when the
    truth is that it was never asked to.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    RDLogger.DisableLog("rdApp.*")

    mol = largest_fragment(cand.smiles)
    if mol is None:
        raise ValueError(f"unparseable SMILES for {cand.ident}")

    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=0xC0FFEE) != 0:
        raise ValueError(f"{cand.ident}: could not embed a 3D conformer")
    AllChem.MMFFOptimizeMolecule(mol)

    # ONE SETUP PER REACTIVE CENTRE. meeko returns a setup for every match of
    # the reactive SMARTS, and taking [0] silently picks one -- a value chosen
    # by default rather than by identity, the defect this project keeps
    # rediscovering. A fumarate has two genuinely distinct electrophilic
    # carbons and needs only one of them to work, so each is docked and the
    # molecule takes the best. Previously a second match was refused outright,
    # which dropped every symmetric Michael acceptor: 4 crystallographic
    # positives, the whole class, removed from validation by a guard.
    setups = MoleculePreparation(reactive_smarts=cand.reactive_smarts,
                                 reactive_smarts_idx=0)(mol)
    if not setups:
        raise ValueError(f"{cand.ident}: reactive SMARTS did not match")

    plain, ok_plain, _ = PDBQTWriterLegacy.write_string(MoleculePreparation()(mol)[0])
    out = []
    for i, setup in enumerate(setups):
        txt, ok, err = PDBQTWriterLegacy.write_string(setup)
        if not ok:
            log.debug("%s centre %d: PDBQT write failed: %s", cand.ident, i, err)
            continue
        # Confirm reactive typing actually happened, by DIFFING against a plain
        # preparation. Never by looking for a named type: meeko derives it from
        # the BASE type, so an aromatic carbon becomes "A1" and an aliphatic one
        # "C1". Hardcoding "C1" silently rejected every SNAr ligand -- 30
        # negatives and 2 positives, a whole warhead class deleted by my check
        # rather than by meeko.
        if ok_plain and _atom_types(txt) == _atom_types(plain):
            continue
        dest = path.with_name(f"{path.stem}_{i}{path.suffix}")
        dest.write_text(txt)
        out.append(dest)
    if not out:
        raise ValueError(f"{cand.ident}: reactive typing had no effect")
    return out


def dock(lig: Path, rec_dir: Path, work: Path, nrun: int, gpu: str) -> Path:
    work.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, CUDA_VISIBLE_DEVICES=gpu)
    subprocess.run(
        [str(AUTODOCK), "-C", "1", "--import_dpf", "rec.reactive_config",
         "--flexres", "rec_flex.pdbqt", "-L", str(lig.resolve()),
         "--nrun", str(nrun), "--resnam", str((work / "out").resolve())],
        cwd=rec_dir, check=True, capture_output=True, env=env)
    return work / "out.dlg"


def sg_position(dlg: Path) -> np.ndarray:
    """Cys113 SG, taken from the FLEXIBLE residue in the first docked model.

    The flexible sidechain moves during docking, so its position must come from
    the pose being measured -- reading it from the rigid receptor would measure
    the approach to where the sulfur STARTED.
    """
    for ln in dlg.read_text(errors="replace").splitlines():
        if "DOCKED: ATOM" not in ln and "DOCKED: HETATM" not in ln:
            continue
        rec = ln.split("DOCKED: ", 1)[1]
        if rec[17:20].strip() == "CYS" and rec[12:16].strip() == "SG":
            return np.array([float(rec[30:38]), float(rec[38:46]), float(rec[46:54])])
    raise ValueError(f"no flexible Cys SG in {dlg}")


def _reactive_xyz(dlg: Path) -> np.ndarray:
    """Coordinates of the retyped reactive atom, from the first docked pose.

    meeko gives the reactive atom an order-1 derivative type — the base type
    with `1` appended, so `C` becomes `C1` and aromatic `A` becomes `A1`. That
    suffix is read rather than any particular type being looked for, since
    hardcoding one deleted a whole warhead class once already.
    """
    for ln in dlg.read_text(errors="replace").splitlines():
        if "DOCKED: ATOM" not in ln and "DOCKED: HETATM" not in ln:
            continue
        rec = ln.split("DOCKED: ", 1)[1]
        if rec[17:20].strip() == "CYS":            # the flexible sidechain
            continue
        t = rec[77:79].strip()
        if len(t) == 2 and t[1] == "1":
            return np.array([float(rec[30:38]), float(rec[38:46]), float(rec[46:54])])
    raise ValueError(f"no reactive-typed ligand atom in {dlg.name}")


def measure_dlg(dlg: Path, cand: Candidate) -> list[nac.NACResult]:
    """Rebuild the ligand from the docking output and score every pose.

    Rebuilding via meeko preserves atom ordering and bonding, so the SMARTS
    match addresses the same atoms the conformers hold. Reading coordinates by
    PDBQT atom name would not: every carbon in a PDBQT is named "C".
    """
    from meeko import PDBQTMolecule, RDKitMolCreate
    from rdkit import Chem

    pm = PDBQTMolecule.from_file(str(dlg), is_dlg=True, skip_typing=True)
    mols = [m for m in RDKitMolCreate.from_pdbqt_mol(pm) if m is not None]
    if not mols:
        raise ValueError(f"{cand.ident}: nothing rebuilt from {dlg.name}")
    mol = mols[0]

    matches = mol.GetSubstructMatches(Chem.MolFromSmarts(cand.reactive_smarts))
    if not matches:
        raise ValueError(f"{cand.ident}: reactive SMARTS does not match the rebuilt mol")
    if len(matches) == 1:
        match = matches[0]
    else:
        # Several reactive centres exist, and exactly ONE of them was docked.
        # Ask the docking output which, rather than assuming the first: the
        # reactive atom is the one meeko retyped to an order-1 derivative type,
        # and it is located by COORDINATE against pose 0. Coordinates are the
        # same numbers in both files, so this is an identity match — not a
        # positional one, which is what would silently pick the wrong carbon.
        rx = _reactive_xyz(dlg)
        pos = mol.GetConformer(0).GetPositions()
        idx = int(np.argmin(np.linalg.norm(pos - rx, axis=1)))
        if np.linalg.norm(pos[idx] - rx) > 0.05:
            raise ValueError(f"{cand.ident}: cannot locate the docked reactive atom")
        hits = [m for m in matches if m[0] == idx]
        if not hits:
            raise ValueError(f"{cand.ident}: no SMARTS match centres on the "
                             f"docked reactive atom")
        # SEVERAL MATCHES MAY SHARE ONE REACTIVE ATOM, and that is not ambiguity.
        # A chloroazine's ipso carbon sits between two ring nitrogens, so
        # `[c]([Cl])[n]` matches twice with the same attacked atom and a
        # different nitrogen. Both triples are coplanar with the ring, so they
        # define the SAME plane and the same criterion; only the reactive atom
        # has to be unambiguous, and it is, because the docking output named it.
        # Refusing these cost both SNAr positives on the previous run.
        # Verified, not assumed: the planes must actually be parallel. If a
        # molecule ever presents two genuinely different planes at one reactive
        # atom, the criterion means two different things and picking either
        # would be arbitrary.
        match = hits[0]
        if len(hits) > 1:
            normals = []
            for h in hits:
                a, b = pos[h[1]] - pos[h[0]], pos[h[2]] - pos[h[0]]
                n = np.cross(a, b)
                if np.linalg.norm(n) > 1e-6:
                    normals.append(n / np.linalg.norm(n))
            for other in normals[1:]:
                if abs(float(normals[0] @ other)) < 0.98:      # ~11 degrees
                    raise ValueError(f"{cand.ident}: matches at one reactive atom "
                                     f"define non-parallel planes")
    return nac.measure_poses(mol, match, cand.mechanism, sg_position(dlg))


# --------------------------------------------------------------------------

def screen(cands: list[Candidate], rec_dir: Path, nrun: int, gpu: str
           ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-molecule summary, and EVERY POSE's raw geometry.

    The per-pose frame is the point. `shared/nac_criterion.py` promises that a
    window can be redrawn without re-docking, and that promise is only real if
    the angles survive the run — a summary keeping just the median cannot
    support a different window, and re-docking to try one invites tuning the
    window against the answer.
    """
    rows, poses = [], []
    for i, c in enumerate(cands, 1):
        work = Path(tempfile.mkdtemp(prefix="nac_"))
        try:
            # Each reactive centre is docked separately; the molecule scores as
            # its BEST centre, because it only needs one of them to react.
            per_centre = []
            for j, lig in enumerate(prepare_ligand(c, work / "lig.pdbqt")):
                dlg = dock(lig, rec_dir, work / f"c{j}", nrun, gpu)
                per_centre.append(measure_dlg(dlg, c))
            res = max(per_centre, key=nac.viable_fraction)
            poses.extend({"ident": c.ident, "label": c.label,
                          "warhead_class": c.warhead_class, "mechanism": c.mechanism,
                          "pose": k, "angle": r.angle, "distance": r.distance}
                         for k, r in enumerate(res))
            angles = np.array([r.angle for r in res])
            dists = np.array([r.distance for r in res])
            rows.append({
                "ident": c.ident, "label": c.label, "warhead_class": c.warhead_class,
                "mechanism": c.mechanism, "n_poses": len(res),
                "viable_fraction": nac.viable_fraction(res),
                "best_angle": angles.max() if c.mechanism.startswith("sn2") else angles.min(),
                "median_angle": float(np.median(angles)),
                "median_dist": float(np.median(dists)),
                "smiles": c.smiles,
            })
            log.info("[%d/%d] %-28s %s viable=%.1f%%", i, len(cands), c.ident[:28],
                     c.label[:3], 100 * rows[-1]["viable_fraction"])
        except Exception as exc:                       # noqa: BLE001
            # Recorded, never silently dropped: a molecule that could not be
            # measured is not a molecule that failed.
            log.warning("[%d/%d] %-28s SKIPPED: %s", i, len(cands), c.ident[:28],
                        str(exc)[:90])
            rows.append({"ident": c.ident, "label": c.label,
                         "warhead_class": c.warhead_class, "mechanism": c.mechanism,
                         "n_poses": 0, "viable_fraction": np.nan,
                         "error": str(exc)[:200], "smiles": c.smiles})
        finally:
            shutil.rmtree(work, ignore_errors=True)
    return pd.DataFrame(rows), pd.DataFrame(poses)


def report(df: pd.DataFrame) -> None:
    ok = df.dropna(subset=["viable_fraction"])
    print(f"\n=== NAC screen: {len(ok)} measured, {len(df) - len(ok)} unmeasurable ===")
    if ok.empty:
        print("  nothing measured"); return

    v = ok.viable_fraction.values
    print(f"\n  THE PRIMARY READOUT — does the score vary at all?")
    print(f"    range   {v.min():.1%} .. {v.max():.1%}")
    print(f"    median  {np.median(v):.1%}   IQR {np.percentile(v,25):.1%}-{np.percentile(v,75):.1%}")
    print(f"    exactly zero: {(v == 0).sum()}/{len(v)} molecules reach NO viable NAC")
    if v.max() - v.min() < 0.05:
        print("\n    FLAT — the criterion does not separate molecules. "
              "This kills the framework as specified.")
    else:
        print(f"\n    SPREAD of {v.max() - v.min():.1%} — there is something to rank.")

    print("\n  SECONDARY — measured actives vs measured inactives")
    print("  (HTS inactivity is weak evidence; read as a hint, not a validation)")
    for lab, g in ok.groupby("label"):
        print(f"    {lab:<9} n={len(g):>3}  viable {g.viable_fraction.mean():.1%} "
              f"(median {g.viable_fraction.median():.1%})")
    pos, neg = ok[ok.label == "positive"], ok[ok.label == "negative"]
    if len(pos) >= 2 and len(neg) >= 2:
        from scipy.stats import mannwhitneyu
        u, p = mannwhitneyu(pos.viable_fraction, neg.viable_fraction,
                            alternative="greater")
        auc = u / (len(pos) * len(neg))
        print(f"    AUC {auc:.3f}, Mann-Whitney p={p:.3f}  "
              f"({'separates' if p < 0.05 else 'not distinguishable'})")

    print("\n  top 10 by viable fraction")
    for r in ok.nlargest(10, "viable_fraction").itertuples():
        print(f"    {r.viable_fraction:>6.1%}  {r.label:<9} {r.ident[:44]}")


# --------------------------------------------------------------------------

COVALENT_LINKS = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith/"
                      "pdb_covalent/covalent_links_3.csv")


def crystal_positives(meta: dict, classes: list[str] | None) -> list[Candidate]:
    """Ligands CRYSTALLOGRAPHICALLY bonded to Cys113 — the strongest positives available.

    17 distinct ligands across 31 solved links. These are not annotations: each
    one has an observed covalent bond to the catalytic cysteine at 1.6-1.9 A, so
    "this molecule reacts with Cys113" is a structural fact rather than a
    label.

    THEY REPLACE THE HTS ACTIVES, WHICH DID NOT SURVIVE SCRUTINY. AID 504891
    yields 34 actives, 11 carrying a warhead, and read as chemistry rather than
    as labels they are a catalogue of frequent hitters: two rhodanines, an
    azlactone, an arylidene barbiturate, a furfurylidene indandione, an
    embelin-like dihydroxyquinone, two naphthoquinone sulfonylimines -- at 3-75
    uM in a 387,000-compound qHTS. One is a cephalosporin whose warhead match is
    spurious. Validating a geometric criterion against compounds that hit
    everything would confirm nothing.

    The stored SMILES are the FREE forms (`warhead-as-drawn`, leaving group
    present), which is what stage 1 docks. The mechanism is taken from the
    warhead library rather than the file's `warhead_class`, which mislabels the
    fumarate/maleate esters as `naphthoquinone_c2` because those classes share a
    reactive SMARTS (D0029).
    """
    from rdkit import Chem
    links = pd.read_csv(COVALENT_LINKS).drop_duplicates("comp_id")
    out = []
    for r in links.dropna(subset=["smiles"]).itertuples():
        m = largest_fragment(r.smiles)
        if m is None:
            continue
        for cid, (mech, smarts) in meta.items():
            if classes and cid not in classes:
                continue
            p = Chem.MolFromSmarts(smarts)
            if p is not None and m.HasSubstructMatch(p):
                out.append(Candidate(f"xtal:{r.pdb_id}:{r.comp_id}", r.smiles,
                                     cid, mech, smarts, "positive"))
                break
    return out


def load_candidates(n_neg: int, classes: list[str] | None) -> list[Candidate]:
    wh = pd.read_csv(REPO / "data/reference/warhead_classes_10.csv")
    meta = {r.class_id: (r.mechanism, r.reactive_atom_smarts) for r in wh.itertuples()}

    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")

    pos = crystal_positives(meta, classes)

    want = {c.warhead_class for c in pos}
    inact = pd.read_csv("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith/"
                        "measured_inactives/aid504891_inactives_1.csv")
    neg, per = [], max(1, n_neg // max(1, len(want)))
    for cid in sorted(want):
        mech, smarts = meta[cid]
        p = Chem.MolFromSmarts(smarts)
        # SHUFFLED, not taken in file order. A PubChem datatable is ordered by
        # submission, which tracks depositor and therefore chemical series -- so
        # the first N matches can be N members of one series, and the negative
        # set would describe that series rather than the library.
        pool = inact.dropna(subset=["canonical_smiles"]).sample(
            frac=1.0, random_state=0xC0FFEE)
        taken = 0
        for r in pool.itertuples():
            if taken >= per:
                break
            m = largest_fragment(r.canonical_smiles)
            if m is None or not m.HasSubstructMatch(p):
                continue
            neg.append(Candidate(f"inactive:{int(r.PUBCHEM_CID)}", r.canonical_smiles,
                                 cid, mech, smarts, "negative"))
            taken += 1
    log.info("candidates: %d positive, %d negative, classes %s",
             len(pos), len(neg), sorted(want))
    return pos + neg


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--nrun", type=int, default=200,
                    help="independent docking runs per molecule")
    ap.add_argument("--n-neg", type=int, default=20)
    ap.add_argument("--classes", nargs="*", default=None,
                    help="restrict to these warhead class_ids")
    ap.add_argument("--gpu", default="1")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    rec = build_reactive_receptor(RX_RECEPTOR)
    cands = load_candidates(args.n_neg, args.classes)
    if not cands:
        raise SystemExit("no candidates")
    df, poses = screen(cands, rec, args.nrun, args.gpu)
    report(df)
    dest = OUT.write("nac_screen", ".csv")
    df.to_csv(dest, index=False)
    pdest = OUT.write("nac_poses", ".csv")
    poses.to_csv(pdest, index=False)
    print(f"\nwritten -> {dest}\n         -> {pdest}  ({len(poses)} poses, "
          f"so a window can be redrawn without re-docking)")


if __name__ == "__main__":
    main()
