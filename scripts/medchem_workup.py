"""
Purpose: single-molecule med-chem workup -- physchem, rule sets, alerts, ligand
         efficiency, stereochemistry, and a like-for-like sulfopin comparison.
Author: @tt8804 (with Claude Code)
Date: 2026-08-06
Input: nothing on disk except the project's own reference data and the T_4 frame
Output: <outdir>/workup_<N>.json  -- every number quoted in
        docs/medchem_t4_72f5671e89cb.md, with the tool that produced it

DESCRIPTORS COME FROM `shared.descriptors` AND NOWHERE ELSE. That module is the
project's single source for the physicochemical axes precisely so that a number
in a report and the same number in a candidate frame cannot disagree. Recomputing
MW here with a different call would be the project's signature defect in its
mildest form.

LIGAND EFFICIENCY IS COMPUTED FROM A SCORE THAT IS NOT VALIDATED ON THIS TARGET.
LE and LLE are quoted because a chemist expects them, and they are labelled with
the energy they came from, because the two available energies mean different
things: `affinity_kcal` is Vina's score for the ADDUCT under GNINA (D0022), and
`best_dg` is AutoDock4's score for the FREE ligand under the reactive potential
of D0063 -- a biased-sampling run, not a plain dock. Neither is an affinity
measurement, and `rank_validated` is False for both.
"""

from __future__ import annotations

import glob
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Crippen, Descriptors, rdFMCS, rdMolDescriptors
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.EnumerateStereoisomers import (EnumerateStereoisomers,
                                               StereoEnumerationOptions)

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import alerts as al            # noqa: E402
from shared import descriptors as desc     # noqa: E402
from shared import reference_set as rs     # noqa: E402
from shared import synthesizability as syn # noqa: E402

RDLogger.DisableLog("rdApp.*")

CAND_ID = "t4_72f5671e89cb"
CAND = "O=S1(=O)CC[C@@H](N(Cc2cnoc2)C2CC(Br)=NO2)C1"
ADDUCT = "O=S1(=O)CC[C@@H](N(Cc2cnoc2)C2CC=NO2)C1"
SULFOPIN = "CC(C)(C)CN(C1CCS(=O)(=O)C1)C(=O)CCl"

#: Sulfopin AS CRYSTALLISED -- PDB component QT7 in 6VAJ, from this project's
#: own `pdb_covalent/covalent_links_3.csv`. It carries the stereocentre that
#: `config/seeds.yaml` omits, so it is the only admissible reference for the
#: configuration question.
QT7 = "CC(C)(C)CN([C@@H]1CCS(=O)(=O)C1)C(=O)CCl"
QT7_ENANT = "CC(C)(C)CN([C@H]1CCS(=O)(=O)C1)C(=O)CCl"
CAND_C3_EPIMER = "O=S1(=O)CC[C@H](N(Cc2cnoc2)C2CC(Br)=NO2)C1"
#: The T_4 protected core with its dummies replaced by methyls, so it can be
#: embedded and its handedness measured on the same footing as the others.
CORE_ME_CAPPED = "N(C)(C)[C@@H]1CCS(=O)(=O)C1"

#: The sulfolane core both molecules share, as the T_4 build defines it.
CORE_SMARTS = "C1CCS(=O)(=O)C1"

T4_FRAMES = "/data/lab_vm/append_only/inhibition/04_t4_combinatorial/D4_*.parquet"
NAC_RANK = "/data/lab_vm/append_only/inhibition/00_outputs/blacksmith/nac_rank"

#: Liability motifs that no public alert catalog carries but a chemist will ask
#: about on sight. Each is a SMARTS plus the reason it matters; a hit is a
#: question for a chemist, never an automatic rejection.
LIABILITY_SMARTS = {
    "N_O_acetal_aminal_ether": (
        "[NX3;!$(N[CX3]=[OX1,SX1,NX2])][CX4;H1,H2]([OX2])",
        "An amine and an ether oxygen on the SAME sp3 carbon. Hydrolyses via "
        "the oxocarbenium/iminium to an aldehyde plus the amine; the rate is "
        "acid-catalysed and can be minutes at gastric pH. Neither PAINS, "
        "BRENK, nor NIH carries this motif.",
    ),
    "vinyl_or_imidoyl_halide": (
        "[Br,Cl,I][CX3]=[NX2,CX3]",
        "A halogen on an sp2 carbon. This is the BDHI warhead itself, so the "
        "hit is expected -- recorded so the count is not read as a surprise.",
    ),
    "tertiary_aliphatic_amine": (
        "[NX3;H0;!$(N[#6]=[O,S,N]);!$(N-[!#6])](-[#6])(-[#6])-[#6]",
        "Basic centre; drives the pH-7.4 charge state and hence permeability.",
    ),
    "sulfone": ("[SX4](=[OX1])(=[OX1])", "Metabolically robust, strongly polar."),
    "isoxazole": (
        "c1cc[nX2]o1",
        "Aromatic isoxazole. Reductive N-O cleavage is a known metabolic route "
        "(CYP and gut flora); the ring is not inert.",
    ),
    "dihydroisoxazole_NO": (
        "[NX2]=[CX3][CX4][CX4][OX2]",
        "The 4,5-dihydroisoxazoline N-O bond. Reducible, and the ring is the "
        "warhead -- listed so the count is attributable.",
    ),
}

#: The sulfolane C3 and the atoms whose arrangement about it defines the
#: configuration. Written so exactly one match exists on both molecules.
SULFOLANE_PATT = "[NX3][CX4;H1;R]1[CH2][S](=O)(=O)[CH2][CH2]1"


def _fp_gen():
    return rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def _tanimoto(a: str, b: str) -> float:
    from rdkit import DataStructs
    g = _fp_gen()
    ma, mb = Chem.MolFromSmiles(a), Chem.MolFromSmiles(b)
    return float(DataStructs.TanimotoSimilarity(g.GetFingerprint(ma),
                                                g.GetFingerprint(mb)))


def rule_sets(smiles: str) -> dict:
    """Lipinski / Veber / Ghose, each as its own components plus a verdict.

    The component values are reported, not just the pass/fail, because a chemist
    reading "Ghose fail" needs to know WHICH bound and by how much.
    """
    m = Chem.MolFromSmiles(smiles)
    mw = Descriptors.MolWt(m)
    logp = Crippen.MolLogP(m)
    hbd = rdMolDescriptors.CalcNumHBD(m)
    hba = rdMolDescriptors.CalcNumHBA(m)
    tpsa = rdMolDescriptors.CalcTPSA(m)
    rotb = rdMolDescriptors.CalcNumRotatableBonds(m)
    mr = Crippen.MolMR(m)
    nat = m.GetNumAtoms()  # Ghose counts ALL atoms, hydrogens included
    nat_h = Chem.AddHs(m).GetNumAtoms()

    lip = {"MW<=500": mw <= 500, "cLogP<=5": logp <= 5,
           "HBD<=5": hbd <= 5, "HBA<=10": hba <= 10}
    veb = {"rot_bonds<=10": rotb <= 10, "TPSA<=140": tpsa <= 140}
    gho = {"160<=MW<=480": 160 <= mw <= 480,
           "-0.4<=cLogP<=5.6": -0.4 <= logp <= 5.6,
           "40<=MR<=130": 40 <= mr <= 130,
           "20<=n_atoms<=70": 20 <= nat_h <= 70}
    return {
        "values": {"MW": mw, "cLogP": logp, "HBD": hbd, "HBA": hba,
                   "TPSA": tpsa, "rot_bonds": rotb, "molar_refractivity": mr,
                   "n_atoms_heavy": nat, "n_atoms_with_H": nat_h},
        "lipinski": {"components": lip, "violations": sum(not v for v in lip.values()),
                     "pass": sum(not v for v in lip.values()) <= 1},
        "veber": {"components": veb, "pass": all(veb.values())},
        "ghose": {"components": gho, "violations": sum(not v for v in gho.values()),
                  "pass": all(gho.values())},
    }


def liabilities(smiles: str) -> dict:
    m = Chem.MolFromSmiles(smiles)
    out = {}
    for name, (sma, why) in LIABILITY_SMARTS.items():
        patt = Chem.MolFromSmarts(sma)
        if patt is None:
            out[name] = {"error": "unparseable SMARTS", "why": why}
            continue
        out[name] = {"n_matches": len(m.GetSubstructMatches(patt)), "why": why}
    return out


def stereo(smiles: str) -> dict:
    """Assigned and unassigned stereocentres, plus what the embedder would do.

    THE POINT OF THE `embedded_*` FIELDS. `scripts/nac_screen.prepare_ligand`
    calls `EmbedMolecule(mol, randomSeed=0xC0FFEE)` on the molecule exactly as
    parsed, once, and every AutoDock-GPU run then samples torsions only. So an
    UNSPECIFIED centre is not sampled -- it is fixed to whatever ETKDG happened
    to build, deterministically, and that configuration is what all 200 runs
    scored. Reconstructing it here is the only way to say which molecule the
    enrichment number describes.
    """
    m = Chem.MolFromSmiles(smiles)
    Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    centres = Chem.FindMolChiralCenters(m, includeUnassigned=True,
                                        useLegacyImplementation=False)
    opts = StereoEnumerationOptions(onlyUnassigned=True, unique=True)
    isomers = sorted(Chem.MolToSmiles(i) for i in EnumerateStereoisomers(
        Chem.MolFromSmiles(smiles), options=opts))

    mh = Chem.AddHs(Chem.MolFromSmiles(smiles))
    rc = AllChem.EmbedMolecule(mh, randomSeed=0xC0FFEE)
    embedded = None
    if rc == 0:
        Chem.AssignStereochemistryFrom3D(mh)
        embedded = {
            "smiles": Chem.MolToSmiles(Chem.RemoveHs(mh)),
            "centres": [(i, str(t)) for i, t in Chem.FindMolChiralCenters(
                Chem.RemoveHs(mh), includeUnassigned=True,
                useLegacyImplementation=False)],
        }
    return {
        "centres_as_written": [(i, str(t)) for i, t in centres],
        "n_assigned": sum(1 for _, t in centres if t in ("R", "S")),
        "n_unassigned": sum(1 for _, t in centres if t not in ("R", "S")),
        "n_stereoisomers_from_unassigned": len(isomers),
        "stereoisomers_from_unassigned": isomers,
        "embed_returncode": rc,
        "embedded": embedded,
    }


def sulfolane_handedness(smiles: str) -> dict:
    """Handedness at the sulfolane C3, measured in 3D, not read off a CIP label.

    WHY NOT JUST COMPARE THE R/S LETTERS. A CIP label is a function of
    substituent priorities, and priorities change when substituents change --
    so two molecules can share a letter and differ in space, or differ in
    letter and agree. Here the priority order provably cannot differ (at C3 the
    four neighbours are N, a CH2 bonded to S, a CH2 bonded to C, and H, and
    every comparison is settled inside the sulfolane ring), but "provably"
    is an argument and this is a measurement: the signed volume of the
    (N, CH2-S, CH2-CH2) triple about C3 is the handedness itself, with no
    convention in between. Same sign = same configuration.
    """
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if AllChem.EmbedMolecule(m, randomSeed=0xC0FFEE) != 0:
        return {"error": "embedding failed"}
    AllChem.MMFFOptimizeMolecule(m)
    matches = m.GetSubstructMatches(Chem.MolFromSmarts(SULFOLANE_PATT))
    if len(matches) != 1:
        return {"error": f"{len(matches)} sulfolane matches, expected 1"}
    n_i, c3_i, c_s_i, _s, _o1, _o2, c_c_i, _c5 = matches[0]
    conf = m.GetConformer()
    pos = lambda i: np.array(conf.GetAtomPosition(i))  # noqa: E731
    c3 = pos(c3_i)
    vol = float(np.dot(np.cross(pos(n_i) - c3, pos(c_s_i) - c3), pos(c_c_i) - c3))
    return {"signed_volume": vol, "sign": "+" if vol > 0 else "-"}


def mcs_with(a: str, b: str) -> dict:
    ma, mb = Chem.MolFromSmiles(a), Chem.MolFromSmiles(b)
    res = rdFMCS.FindMCS([ma, mb], ringMatchesRingOnly=True,
                         completeRingsOnly=True, timeout=30)
    return {"smarts": res.smartsString, "n_atoms": res.numAtoms,
            "n_bonds": res.numBonds, "canceled": res.canceled}


def nac_latest() -> pd.DataFrame:
    """The NAC ranking with the newest measurement of each candidate winning.

    Reimplements `scripts.nac_rank.load_scored` deliberately: importing that
    module pulls in meeko and the docking stack, which this workup does not
    need. The ONE thing that must match is the version-order sort -- lexical
    order puts `_10` before `_2` and would make "newest" wrong for exactly the
    BDHI rows D0067 re-scored.
    """
    def key(p: str) -> tuple[int, int]:
        m = re.search(r"nac_rank_s(\d+)_(\d+)\.csv$", p)
        return (int(m.group(1)), int(m.group(2)))

    fs = sorted(glob.glob(f"{NAC_RANK}/nac_rank_s*.csv"), key=key)
    df = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    return df.drop_duplicates("ident", keep="last")


def main() -> None:
    outdir = Path(sys.argv[1] if len(sys.argv) > 1
                  else "/data/lab_vm/append_only/inhibition/00_outputs/"
                       "blacksmith/medchem_t4_72f5671e89cb")
    version = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    out = outdir / f"workup_{version}.json"
    if out.exists():
        raise SystemExit(f"{out} exists -- bump the version rather than overwriting")

    mols = {"candidate": CAND, "adduct": ADDUCT, "sulfopin": SULFOPIN}
    rec: dict = {"rdkit_version": Chem.rdBase.rdkitVersion, "molecules": {}}

    for name, s in mols.items():
        entry = {
            "smiles": s,
            "formula": rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(s)),
            "descriptors_shared": desc.compute(s),
            "rules": rule_sets(s),
            "alerts_whole_molecule": al.screen(s).__dict__,
            "synthesizability_violations": [r.name for r in syn.violations(s)],
            "liability_smarts": liabilities(s),
            "stereo": stereo(s),
        }
        rec["molecules"][name] = entry

    # --- two-tier alerts: the warhead is the mechanism, the decoration is not
    for name in ("candidate", "sulfopin"):
        try:
            tt = al.two_tier(mols[name], CORE_SMARTS)
            rec["molecules"][name]["alerts_two_tier"] = {
                k: (v if not hasattr(v, "__dict__") else v.__dict__)
                for k, v in tt.__dict__.items()
            }
        except Exception as e:  # noqa: BLE001 - record, never crash the workup
            rec["molecules"][name]["alerts_two_tier"] = {"error": repr(e)}

    # --- the stereocentre question, measured against the CRYSTALLISED sulfopin.
    # config/seeds.yaml stores sulfopin WITHOUT stereochemistry; the PDB
    # component QT7 in 6VAJ carries it, and that is the molecule that was
    # solved bonded to Cys113. Comparing against the seed SMILES would compare
    # against a drawing that omits the very thing being asked about.
    rec["stereo_comparison"] = {
        "qt7_6vaj_smiles": QT7,
        "seed_yaml_sulfopin_has_stereo": False,
        "handedness": {
            "candidate": sulfolane_handedness(CAND),
            "candidate_C3_epimer": sulfolane_handedness(CAND_C3_EPIMER),
            "sulfopin_qt7_6vaj": sulfolane_handedness(QT7),
            "sulfopin_enantiomer": sulfolane_handedness(QT7_ENANT),
            "t4_protected_core_me_capped": sulfolane_handedness(CORE_ME_CAPPED),
        },
    }

    # --- shared scaffold
    rec["mcs_candidate_vs_sulfopin"] = mcs_with(CAND, SULFOPIN)
    rec["tanimoto_candidate_vs_sulfopin_ecfp4"] = _tanimoto(CAND, SULFOPIN)

    # --- where it sits in known Pin1 chemical space
    refset = rs.load()
    rec["reference_set_sha256"] = {"master": refset.master_sha256,
                                   "anchors": refset.anchors_sha256}
    sims = []
    for _, r in refset.master.iterrows():
        s = r.get("canonical_smiles")
        if not isinstance(s, str) or Chem.MolFromSmiles(s) is None:
            continue
        sims.append({"name": r.get("name"), "tanimoto_ecfp4": _tanimoto(CAND, s)})
    sims.sort(key=lambda d: -d["tanimoto_ecfp4"])
    rec["nearest_reference_binders"] = sims[:10]
    rec["n_reference_binders_compared"] = len(sims)

    # verified_only=False deliberately: the point here is to SHOW that the BDHI
    # anchor is the one row the gate refuses, not to use it.
    rec["covalent_anchors"] = rs.covalent_anchors(
        refset, verified_only=False).to_dict("records")

    # --- the T_4 frame row: everything the pipeline already measured
    frames = sorted(glob.glob(T4_FRAMES),
                    key=lambda p: int(re.search(r"D4_(\d+)", p).group(1)))
    d4 = pd.read_parquet(frames[-1])
    row = d4[d4.candidate_id == CAND_ID]
    rec["t4_frame"] = {"path": frames[-1],
                       "row": json.loads(row.to_json(orient="records"))[0]}

    # --- rank position, on the corrected criterion
    nac = nac_latest()
    ok = nac[nac.status == "ok"].copy()
    t4 = ok[ok.approach == "T_4"]
    cls = ok[ok.warhead_class == "bdhi_c5"]
    me = ok[ok.ident == CAND_ID].iloc[0]
    rec["nac_rank_position"] = {
        "enrichment": float(me.enrichment),
        "viable_fraction": float(me.viable_fraction),
        "n_poses": int(me.n_poses),
        "best_dg": float(me.best_dg),
        "median_viable_dg": float(me.median_viable_dg),
        "rank_within_T4": int((t4.enrichment > me.enrichment).sum() + 1),
        "n_T4_scored": int(len(t4)),
        "rank_within_bdhi_c5": int((cls.enrichment > me.enrichment).sum() + 1),
        "n_bdhi_c5_scored": int(len(cls)),
        "T4_enrichment_median": float(t4.enrichment.median()),
        "all_scored": int(len(ok)),
        "rank_overall": int((ok.enrichment > me.enrichment).sum() + 1),
        "per_class_median": ok.groupby("warhead_class").enrichment.median().to_dict(),
        "per_class_n": ok.groupby("warhead_class").size().to_dict(),
    }

    # --- ligand efficiency, from each energy separately and labelled as such
    hac = int(desc.compute(CAND)["HAC"])
    clogp = float(desc.compute(CAND)["cLogP"])
    le = {}
    for label, dg, provenance in (
        ("autodock4_reactive_best_dg", float(me.best_dg),
         "AutoDock-GPU, 3IKD, FREE ligand, reactive potential of D0063 "
         "(r_eq_12=3.2 A). Biased sampling: not a plain dock, not an affinity."),
        ("vina_gnina_adduct_affinity", float(row.affinity_kcal.iloc[0]),
         "GNINA/Vina, ADDUCT form (D0022), min affinity over 9 modes. "
         "gate_verdict UNDERPOWERED, rank_validated False."),
    ):
        pkd = -dg / 1.364  # dG = -RT ln K, 298 K, kcal/mol -> log10 units
        le[label] = {
            "dG_kcal_per_mol": dg,
            "LE_kcal_per_mol_per_heavy_atom": -dg / hac,
            "implied_pKd_if_this_were_an_affinity": pkd,
            "LLE_surrogate_pKd_minus_cLogP": pkd - clogp,
            "provenance": provenance,
        }
    rec["ligand_efficiency"] = {"HAC": hac, "cLogP": clogp, "by_energy": le}

    # sulfopin's LE against its own MM-GBSA number, for scale only
    rec["sulfopin_reference_energies"] = {
        "mmgbsa_ensemble_dG_kcal": -7.58,
        "mmgbsa_ensemble_sd": 0.28,
        "source": "D0036 (ensemble MM-GBSA); 50 of 80 decoys scored better",
        "HAC": int(desc.compute(SULFOPIN)["HAC"]),
        "LE_kcal_per_mol_per_heavy_atom": 7.58 / int(desc.compute(SULFOPIN)["HAC"]),
    }

    outdir.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rec, indent=2, default=str))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
