"""
Purpose: Turn the Pin1 PDB survey into redocking CASES -- per (entry, ligand)
         a prepared receptor, a crystal-frame reference pose, a self-docking
         box, and the same reference transformed into 6VAJ's frame.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-31
Input: 00_outputs/blacksmith/redock_pin1/pdb_ligand_survey_<latest>.csv
Output: append_only/inhibition/05_redock_benchmark/cases_1/...
        00_outputs/blacksmith/redock_pin1/redock_cases_<N>.csv (incl. every DROP)

RECEPTOR PREPARATION IS 6VAJ's, NOT A NEW ONE. Each case's receptor goes
through `shared.receptor_prep`'s exact path -- strip solvent, `reduce -BUILD`,
`obabel -xr` -- because a benchmark prepared differently from the production
receptor measures the preparation, not the docking. D0009 is the reason
`reduce` and not `obabel -p` places the hydrogens: obabel renumbered and
RENAMED residues on 6VAJ and silently dropped Cys113.

EVERY LIGAND IN THE ENTRY IS STRIPPED, NOT JUST THE TARGET. A fragment-screen
entry can hold two drug-like components; leaving the other one in the pocket
would dock the target against a receptor the production protocol would never
present. Solvent goes by `receptor_prep.STRIPPABLE_HET`, unrecognised
heteroatoms are RETAINED and counted (D0003).

COVALENCY IS DETECTED FROM GEOMETRY, NOT ONLY FROM LINK RECORDS. Depositors are
inconsistent about LINK. A ligand heavy atom within COVALENT_MAX_A of a protein
heavy atom is bonded whatever the header says, so both are checked and the
union flags the case. Redocking a covalent adduct non-covalently is a different
experiment (D0022) and these are reported as their own stratum, never pooled
into the headline number.

SUPERPOSITION IS ON PPIase CA ATOMS ONLY. Pin1 is two domains joined by a
flexible linker, and the WW domain's position relative to the PPIase domain
varies between entries. Fitting on all CAs would let WW-domain motion push the
active-site ligand off by an amount that has nothing to do with docking. The
fit is therefore restricted to PPIASE_RANGE, matched residue-by-residue with an
identity check so a point mutant contributes only its unmutated positions.
"""

from __future__ import annotations

import json
import logging
import os
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

# receptor_prep calls a bare `obabel`, which is only on PATH inside the cheminf
# env. noncovalent_dock_run pins the absolute path for the same binary; put the
# env on PATH so BOTH resolve to the one interpreter this benchmark runs under.
CHEMINF_BIN = "/data/lab_vm/envs/dwi_cheminf/bin"
os.environ["PATH"] = CHEMINF_BIN + os.pathsep + os.environ.get("PATH", "")

from shared import receptor_prep as rp              # noqa: E402
from shared import outputs as sout           # noqa: E402

log = logging.getLogger("redock-cases")

# Analysis outputs live under the GOVERNED root, not in the repo
# (rules/data-storage.md). See shared/outputs.py for why, and for the
# versioned-write / resolve-latest policy the append-only tree needs.
OUT = sout.Topic("blacksmith", "redock_pin1")
OUT_DIR = OUT.dir
WORK = Path("/data/lab_vm/append_only/inhibition/05_redock_benchmark/cases_1")
REF_6VAJ = Path("/data/lab_vm/immutable/inhibition/receptor/6VAJ_prepared.pdb")
BOX_EXPANDED = Path("/data/lab_vm/immutable/inhibition/receptor/box_expanded.json")
PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
CIF_URL = "https://files.rcsb.org/download/{pdb_id}.cif"

# Pin1's PPIase (catalytic) domain in author numbering. The WW domain is 1-39
# and the linker runs to ~50.
PPIASE_RANGE = (50, 163)
MIN_FIT_RESIDUES = 40
MAX_FIT_RMSD_A = 3.0

# A ligand-to-protein heavy-atom contact at or below this is a covalent bond,
# not a close contact. C-S is 1.81 A, C-C 1.54 A; the shortest genuine
# non-bonded contact in a crystal structure is ~2.6 A.
COVALENT_MAX_A = 1.95

ALLOWED_ALTLOC = {"", "A"}


# --------------------------------------------------------------------------
# PDB parsing. Fixed-column slicing, never split() -- coordinates run together
# when a value reaches -100.000 (receptor_prep makes the same point).
# --------------------------------------------------------------------------

def _atom(line: str) -> dict:
    return {
        "record": line[:6].strip(),
        "name": line[12:16].strip(),
        "altloc": line[16].strip(),
        "resname": line[17:20].strip(),
        "chain": line[21],
        "resseq": line[22:26].strip(),
        "icode": line[26].strip(),
        "xyz": (float(line[30:38]), float(line[38:46]), float(line[46:54])),
        "element": line[76:78].strip().upper() if len(line) >= 78 else "",
        "line": line,
    }


def is_hydrogen(a: dict) -> bool:
    """True for H/D. The element column is authoritative; the name is a fallback.

    HYDROGENS MUST NEVER REACH THE ATOM COUNT OR THE DISTANCE TEST. Entries
    refined at high resolution deposit explicit H -- 7OQ9's ligand 0AW is 46
    atoms on a 27-heavy-atom component, exactly 27 + 19 H. Counting them made
    complete ligands look like the wrong molecule. Worse, an H-bond hydrogen
    sits 1.8-2.0 A from its acceptor, which is inside COVALENT_MAX_A: leaving H
    in the distance test would have flagged ordinary hydrogen-bonded ligands as
    COVALENT and moved them out of the headline set.
    """
    e = a.get("element", "")
    if e:
        return e in ("H", "D")
    nm = a.get("name", "").lstrip("0123456789")
    return bool(nm) and nm[0] in ("H", "D")


def read_pdb(path: Path) -> tuple[list[dict], list[str]]:
    """Return (atom records, LINK records) for MODEL 1 only.

    ONLY THE FIRST MODEL. 27 of the 190 Pin1 entries are solution NMR and carry
    a 20-model ensemble. Without the ENDMDL guard every model's atoms merge
    into one residue -- 2M9F's NAG came back as 560 atoms on a 15-atom
    component -- and the CA arrays used for superposition are 20x duplicated.
    """
    atoms, links = [], []
    seen_model = False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("MODEL "):
            if seen_model:
                break
            seen_model = True
        elif line.startswith("ENDMDL"):
            break
        elif line.startswith(("ATOM  ", "HETATM")):
            try:
                atoms.append(_atom(line))
            except (ValueError, IndexError):
                continue
        elif line.startswith("LINK"):
            links.append(line)
    return atoms, links


# --------------------------------------------------------------------------
# mmCIF fallback.
#
# MOST OF THIS BENCHMARK IS ONLY IN mmCIF. 66 of the 127 ligand-bearing entries
# -- the entire 2025 Xiao fragment series (9KX*, 9V6*) -- return HTTP 404 for
# `.pdb`; RCSB no longer emits legacy PDB format for new depositions. Falling
# back to mmCIF is not a nicety here, it is half the dataset, and a benchmark
# that quietly skipped them would be reporting on the pre-2024 chemistry only.
# --------------------------------------------------------------------------

def _cif_tokens(line: str) -> list[str]:
    """Split an mmCIF data line, honouring single/double quotes."""
    out, cur, quote = [], "", ""
    for ch in line:
        if quote:
            if ch == quote:
                quote = ""
                out.append(cur)
                cur = ""
            else:
                cur += ch
        elif ch in "'\"":
            quote = ch
        elif ch.isspace():
            if cur:
                out.append(cur)
                cur = ""
        else:
            cur += ch
    if cur:
        out.append(cur)
    return out


def _pdb_line(rec: str, serial: int, name: str, altloc: str, resname: str,
              chain: str, resseq: str, icode: str, xyz: tuple, element: str
              ) -> str:
    """Render one atom as a fixed-column PDB line."""
    # Atom names of <=3 chars for single-letter elements start in column 14.
    nm = name if len(name) >= 4 else (
        f" {name:<3s}" if len(element) <= 1 else f"{name:<4s}")
    return (f"{rec:<6s}{serial:5d} {nm[:4]}{altloc[:1]:1s}{resname[:3]:>3s} "
            f"{chain[:1]:1s}{resseq:>4s}{icode[:1]:1s}   "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
            f"{1.0:6.2f}{0.0:6.2f}          {element[:2]:>2s}")


def parse_cif(path: Path) -> tuple[list[dict], list[str]]:
    """Minimal mmCIF reader: `_atom_site` (model 1) and covalent `_struct_conn`.

    Only the fields this benchmark needs are extracted. Author numbering is
    used throughout so the CIF and PDB paths produce interchangeable records.
    """
    text = path.read_text(encoding="utf-8", errors="replace").splitlines()
    atoms: list[dict] = []
    links: list[str] = []
    i, n = 0, len(text)
    while i < n:
        if text[i].strip() != "loop_":
            i += 1
            continue
        j = i + 1
        tags: list[str] = []
        while j < n and text[j].lstrip().startswith("_"):
            tags.append(text[j].strip())
            j += 1
        cat = tags[0].split(".")[0] if tags else ""
        if cat not in ("_atom_site", "_struct_conn"):
            i = j
            continue
        idx = {t.split(".", 1)[1]: k for k, t in enumerate(tags)}
        rows: list[list[str]] = []
        while j < n:
            s = text[j].strip()
            if not s or s.startswith("#") or s.startswith("loop_") or s.startswith("_"):
                break
            rows.append(_cif_tokens(text[j]))
            j += 1

        if cat == "_atom_site":
            def g(r: list[str], key: str, default: str = "") -> str:
                k = idx.get(key)
                v = r[k] if k is not None and k < len(r) else default
                return "" if v in (".", "?") else v
            for r in rows:
                if len(r) < len(tags):
                    continue
                if g(r, "pdbx_PDB_model_num", "1") not in ("", "1"):
                    continue     # NMR ensembles: model 1 only
                try:
                    xyz = (float(g(r, "Cartn_x")), float(g(r, "Cartn_y")),
                           float(g(r, "Cartn_z")))
                except ValueError:
                    continue
                rec = g(r, "group_PDB", "ATOM")
                name = g(r, "auth_atom_id") or g(r, "label_atom_id")
                resname = g(r, "auth_comp_id") or g(r, "label_comp_id")
                chain = g(r, "auth_asym_id") or g(r, "label_asym_id")
                resseq = g(r, "auth_seq_id") or g(r, "label_seq_id")
                alt = g(r, "label_alt_id")
                icode = g(r, "pdbx_PDB_ins_code")
                elem = g(r, "type_symbol").upper()
                atoms.append({
                    "record": rec, "name": name, "altloc": alt,
                    "resname": resname, "chain": chain, "resseq": resseq,
                    "icode": icode, "xyz": xyz, "element": elem,
                    "line": _pdb_line(rec, len(atoms) + 1, name, alt, resname,
                                      chain, resseq, icode, xyz, elem),
                })
        else:  # _struct_conn -- keep only covalent links, as pseudo-LINK lines
            for r in rows:
                if len(r) < len(tags):
                    continue
                k = idx.get("conn_type_id")
                if k is None or k >= len(r) or r[k] != "covale":
                    continue
                a = idx.get("ptnr1_auth_comp_id") or idx.get("ptnr1_label_comp_id")
                b = idx.get("ptnr2_auth_comp_id") or idx.get("ptnr2_label_comp_id")
                c1 = r[a] if a is not None and a < len(r) else ""
                c2 = r[b] if b is not None and b < len(r) else ""
                links.append(f"LINK{'':13s}{c1:>3s}{'':27s}{c2:>3s}")
        i = j
    return atoms, links


def fetch_structure(pdb_id: str, dest_dir: Path) -> Path | None:
    """Download an entry, preferring PDB format and falling back to mmCIF."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    for ext, url in (("pdb", PDB_URL), ("cif", CIF_URL)):
        p = dest_dir / f"{pdb_id}.{ext}"
        if p.is_file() and p.stat().st_size > 1000:
            return p
        try:
            with urllib.request.urlopen(url.format(pdb_id=pdb_id), timeout=180) as fh:
                data = fh.read()
        except Exception:  # noqa: BLE001 - try the next format
            continue
        if len(data) > 1000:
            p.write_bytes(data)
            return p
    log.warning("%s: neither PDB nor mmCIF could be downloaded", pdb_id)
    return None


def load_structure(path: Path) -> tuple[list[dict], list[str]]:
    """Read either format into the same (atoms, links) shape."""
    return parse_cif(path) if path.suffix == ".cif" else read_pdb(path)


# --------------------------------------------------------------------------
# Superposition (Kabsch). Implemented here rather than pulling in gemmi or
# biopython -- neither is in the shared envs and this is 15 lines of numpy.
# --------------------------------------------------------------------------

def kabsch(mobile: np.ndarray, target: np.ndarray
           ) -> tuple[np.ndarray, np.ndarray, float]:
    """Rigid transform taking `mobile` onto `target`. Returns (R, t, rmsd)."""
    mc, tc = mobile.mean(axis=0), target.mean(axis=0)
    p, q = mobile - mc, target - tc
    v, _, wt = np.linalg.svd(p.T @ q)
    d = np.sign(np.linalg.det(v @ wt))
    rot = v @ np.diag([1.0, 1.0, d]) @ wt
    fitted = p @ rot + tc
    rmsd = float(np.sqrt(((fitted - target) ** 2).sum(axis=1).mean()))
    return rot, tc - mc @ rot, rmsd


def ca_map(atoms: list[dict], chain: str | None = None) -> dict[int, tuple]:
    """{resseq: (resname, xyz)} for CA atoms in the PPIase range."""
    out: dict[int, tuple] = {}
    lo, hi = PPIASE_RANGE
    for a in atoms:
        if a["record"] != "ATOM" or a["name"] != "CA":
            continue
        if a["altloc"] not in ALLOWED_ALTLOC or a["icode"]:
            continue
        if chain is not None and a["chain"] != chain:
            continue
        try:
            n = int(a["resseq"])
        except ValueError:
            continue
        if lo <= n <= hi and n not in out:
            out[n] = (a["resname"], a["xyz"])
    return out


def superpose_to_6vaj(atoms: list[dict], chain: str, ref_ca: dict[int, tuple]
                      ) -> tuple[np.ndarray, np.ndarray, float, int]:
    """Fit one entry's PPIase CA trace onto 6VAJ's. Raises on a bad fit."""
    mine = ca_map(atoms, chain)
    common = [n for n in sorted(set(mine) & set(ref_ca))
              if mine[n][0] == ref_ca[n][0]]   # identity check: skip mutations
    if len(common) < MIN_FIT_RESIDUES:
        raise ValueError(f"only {len(common)} matched PPIase CA residues "
                         f"(need {MIN_FIT_RESIDUES})")
    mob = np.array([mine[n][1] for n in common])
    tgt = np.array([ref_ca[n][1] for n in common])
    rot, trans, rmsd = kabsch(mob, tgt)
    if rmsd > MAX_FIT_RMSD_A:
        raise ValueError(f"CA fit RMSD {rmsd:.2f} A over {len(common)} residues "
                         f"exceeds {MAX_FIT_RMSD_A} A")
    return rot, trans, rmsd, len(common)


# --------------------------------------------------------------------------
# Case construction
# --------------------------------------------------------------------------

def write_ligand_pdb(inst: list[dict], path: Path,
                     xyz: np.ndarray | None = None) -> None:
    """Write a ligand instance as PDB, optionally with substituted coords."""
    lines = []
    for i, a in enumerate(inst):
        c = a["xyz"] if xyz is None else tuple(xyz[i])
        # Blank the altloc so RDKit does not see two conformers of one atom.
        lines.append(f"{a['line'][:16]} {a['line'][17:30]}"
                     f"{c[0]:8.3f}{c[1]:8.3f}{c[2]:8.3f}{a['line'][54:]}")
    path.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")


def build_receptor(atoms: list[dict], strip_codes: set[str], dst_dir: Path
                   ) -> tuple[Path, dict]:
    """Strip -> reduce -> pdbqt, following receptor_prep's exact path."""
    dst_dir.mkdir(parents=True, exist_ok=True)
    stripped = dst_dir / "receptor_stripped.pdb"
    kept, n_solvent, n_lig, n_other = [], 0, 0, 0
    for a in atoms:
        if a["record"] == "HETATM":
            if a["resname"] in rp.STRIPPABLE_HET:
                n_solvent += 1
                continue
            if a["resname"] in strip_codes:
                n_lig += 1
                continue
            if is_hydrogen(a):
                continue
            n_other += 1
            kept.append(a["line"])
            continue
        if a["altloc"] not in ALLOWED_ALTLOC or is_hydrogen(a):
            continue
        kept.append(a["line"])
    stripped.write_text("\n".join(kept) + "\nEND\n", encoding="utf-8")

    protonated = dst_dir / "receptor_h.pdb"
    rp.protonate(stripped, protonated, 7.4)
    pdbqt = dst_dir / "receptor.pdbqt"
    rp.to_pdbqt(protonated, pdbqt)
    return pdbqt, {"receptor_atoms_kept": len(kept),
                   "solvent_atoms_removed": n_solvent,
                   "ligand_atoms_removed": n_lig,
                   "other_het_retained": n_other}


def main() -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    WORK.mkdir(parents=True, exist_ok=True)
    box_size = json.loads(BOX_EXPANDED.read_text())

    survey = pd.read_csv(OUT.latest("pdb_ligand_survey", ".csv"))
    ligands = survey[survey.classification == "ligand"].copy()
    # Every drug-like component in an entry is stripped from that entry's
    # receptor, not just the one being redocked.
    strip_by_entry = ligands.groupby("pdb_id")["comp_id"].apply(set).to_dict()

    ref_atoms, _ = read_pdb(REF_6VAJ)
    ref_ca = ca_map(ref_atoms, "A")
    log.info("6VAJ reference: %d PPIase CA atoms in %d-%d",
             len(ref_ca), *PPIASE_RANGE)

    rows: list[dict] = []
    for pdb_id, grp in ligands.groupby("pdb_id"):
        path = fetch_structure(pdb_id, WORK / "structures")
        if path is None:
            for _, r in grp.iterrows():
                rows.append({**r.to_dict(), "status": "drop",
                             "drop_reason": "structure_download_failed"})
            continue
        atoms, links = load_structure(path)
        prot = [a for a in atoms if a["record"] == "ATOM" and not is_hydrogen(a)]
        if not prot:
            for _, r in grp.iterrows():
                rows.append({**r.to_dict(), "status": "drop",
                             "drop_reason": "no_protein_atoms"})
            continue
        prot_xyz = np.array([a["xyz"] for a in prot])

        receptor_cache: dict[str, tuple[Path, dict]] = {}
        for _, r in grp.iterrows():
            comp = r["comp_id"]
            base = {**r.to_dict()}
            # ---- pick one instance of this ligand -------------------------
            inst_groups: dict[tuple, list[dict]] = {}
            n_altloc = 0
            for a in atoms:
                if a["record"] != "HETATM" or a["resname"] != comp:
                    continue
                if a["altloc"] not in ALLOWED_ALTLOC:
                    n_altloc += 1
                    continue
                if is_hydrogen(a):
                    continue
                inst_groups.setdefault((a["chain"], a["resseq"], a["icode"]), []).append(a)
            if not inst_groups:
                rows.append({**base, "status": "drop",
                             "drop_reason": "ligand_absent_from_coordinates"})
                continue
            key = sorted(inst_groups, key=lambda k: (k[0] != "A", k[0], int(k[1])))[0]
            inst = inst_groups[key]
            chain, resseq = key[0], key[1]
            case_id = f"{pdb_id}_{comp}_{chain}{resseq}"
            lig_xyz = np.array([a["xyz"] for a in inst])

            # ---- covalency: LINK header OR bond-length geometry -----------
            d = np.linalg.norm(lig_xyz[:, None, :] - prot_xyz[None, :, :], axis=2)
            min_d = float(d.min())
            has_link = any(comp in (l[17:20].strip(), l[47:50].strip())
                           for l in links)
            is_covalent = bool(min_d <= COVALENT_MAX_A or has_link)

            base.update({
                "case_id": case_id, "lig_chain": chain, "lig_resseq": resseq,
                "n_ligand_atoms_modelled": len(inst),
                "n_altloc_atoms_skipped": n_altloc,
                "min_dist_to_protein_a": round(min_d, 3),
                "link_record": has_link, "is_covalent": is_covalent,
            })

            # ---- covalent adducts are their own stratum -------------------
            #
            # A COVALENT LIGAND IS SHORT ITS LEAVING GROUP, AND THAT IS NOT
            # DISORDER. 6VAJ's QT7 is deposited with 16 heavy atoms against a
            # chem-comp SMILES of 17: the chloride left when the
            # chloroacetamide alkylated Cys113. The first version of this
            # check read that as incomplete density and dropped it -- and the
            # same "N-1 of N" signature appeared across the whole covalent set,
            # which is what gave the misreading away.
            #
            # These are excluded from the headline for the reason D0022 gives:
            # the species that sits in the crystal is the ADDUCT, while what we
            # would dock is the pre-reaction compound, so a symmetry-corrected
            # RMSD between them is not even well defined -- the two molecular
            # graphs differ. Excluded with a reason and counted, never silent.
            expected = int(r["heavy_atoms"])
            if is_covalent:
                base["status"] = "excluded_covalent"
                base["drop_reason"] = (
                    f"covalent_at_{min_d:.2f}A"
                    + (f"_leaving_group_absent_{len(inst)}_of_{expected}"
                       if len(inst) < expected else "")
                    + ("_link_record" if has_link else ""))
                rows.append(base)
                continue

            # An incompletely modelled ligand cannot support an RMSD against a
            # fully built docked pose; caught here rather than at the compare.
            if len(inst) != expected:
                base["status"], base["drop_reason"] = "drop", (
                    f"atom_count_mismatch_{len(inst)}_modelled_vs_"
                    f"{expected}_in_chem_comp")
                rows.append(base)
                continue

            # ---- receptor (one per entry, shared across its ligands) ------
            try:
                if pdb_id not in receptor_cache:
                    receptor_cache[pdb_id] = build_receptor(
                        atoms, strip_by_entry[pdb_id], WORK / "receptors" / pdb_id)
                pdbqt, counts = receptor_cache[pdb_id]
            except Exception as exc:  # noqa: BLE001
                base["status"], base["drop_reason"] = "drop", (
                    f"receptor_prep_failed: {str(exc)[:120]}")
                rows.append(base)
                continue
            base.update(counts)
            base["receptor_pdbqt"] = str(pdbqt)

            # ---- crystal-frame reference + self-docking box ---------------
            ref_dir = WORK / "refs"
            ref_dir.mkdir(parents=True, exist_ok=True)
            ref_pdb = ref_dir / f"{case_id}_ref.pdb"
            write_ligand_pdb(inst, ref_pdb)
            base["ref_pdb"] = str(ref_pdb)
            centroid = lig_xyz.mean(axis=0)
            box = {"center_x": float(centroid[0]), "center_y": float(centroid[1]),
                   "center_z": float(centroid[2]),
                   "size_x": box_size["size_x"], "size_y": box_size["size_y"],
                   "size_z": box_size["size_z"],
                   "derived_from": f"{comp} centroid in {pdb_id}"}
            box_dir = WORK / "boxes"
            box_dir.mkdir(parents=True, exist_ok=True)
            box_path = box_dir / f"{case_id}.json"
            box_path.write_text(json.dumps(box, indent=2), encoding="utf-8")
            base["box_json"] = str(box_path)

            # ---- 6VAJ-frame reference (for the cross-docking arm) ---------
            # The chain fitted is the one the ligand actually sits on, not
            # chain A by assumption -- in a two-copy asymmetric unit the
            # ligand may belong to either.
            near = {}
            for a, dist in zip(prot, d.min(axis=0)):
                near[a["chain"]] = min(near.get(a["chain"], 1e9), dist)
            host_chain = min(near, key=near.get)
            base["host_chain"] = host_chain
            try:
                rot, trans, fit_rmsd, n_fit = superpose_to_6vaj(
                    atoms, host_chain, ref_ca)
                moved = lig_xyz @ rot + trans
                ref6_dir = WORK / "refs_6vaj"
                ref6_dir.mkdir(parents=True, exist_ok=True)
                ref6 = ref6_dir / f"{case_id}_ref6vaj.pdb"
                write_ligand_pdb(inst, ref6, xyz=moved)
                base.update({"ref_6vaj_pdb": str(ref6),
                             "ca_fit_rmsd_a": round(fit_rmsd, 3),
                             "n_fit_residues": n_fit,
                             "cross_dock_ready": True})
            except Exception as exc:  # noqa: BLE001
                # A CHAIN THAT WILL NOT ALIGN TO PIN1 IS NOT PIN1, AND THE CASE
                # DIES HERE RATHER THAN BECOMING A SELF-DOCKING-ONLY ENTRY.
                # 7OQ9, 7OQA and 8C3C match the Q13526 search because they
                # contain a five-residue Pin1 PHOSPHOPEPTIDE -- but chain A is a
                # ~230-residue 14-3-3 protein, and the ligand (fusicoccin in
                # 8C3C) binds 14-3-3, not Pin1. Self-docking them would have
                # quietly benchmarked a different target and folded the result
                # into a Pin1 success rate.
                base["status"] = "drop"
                base["drop_reason"] = f"not_pin1_ppiase_domain: {str(exc)[:120]}"
                rows.append(base)
                continue

            # Can the PRODUCTION box even reach this site? A ligand whose
            # crystallographic site falls outside the 26 A box centred on QT7
            # cannot be recovered by cross-docking however good the engine is --
            # Vina is never offered that volume. Recorded rather than dropped:
            # "the production box does not cover this site" is a finding about
            # the protocol, and the Xiao 2025 Site 2 fragments are expected here.
            half = np.array([box_size["size_x"], box_size["size_y"],
                             box_size["size_z"]]) / 2.0
            prod_c = np.array([box_size["center_x"], box_size["center_y"],
                               box_size["center_z"]])
            inside = (np.abs(moved - prod_c) <= half).all(axis=1)
            base.update({
                "ref_frac_atoms_in_prod_box": round(float(inside.mean()), 3),
                "ref_centroid_dist_to_prod_centre_a": round(
                    float(np.linalg.norm(moved.mean(axis=0) - prod_c)), 2),
                "ref_in_prod_box": bool(inside.all()),
            })

            base["status"] = "case"
            rows.append(base)
        log.info("%s: %d component(s) processed", pdb_id, len(grp))

    df = pd.DataFrame(rows)
    df.to_csv(OUT.write("redock_cases", ".csv"), index=False)
    log.info("=" * 60)
    log.info("status: %s", df["status"].value_counts().to_dict())
    if "drop_reason" in df:
        log.info("drops: %s", df.drop_reason.dropna().value_counts().to_dict())
    cases = df[df.status == "case"]
    if cases.empty:
        log.error("NO CASES BUILT -- nothing downstream can run.")
        return
    log.info("cases: %d over %d entries; cross-dock ready: %d",
             len(cases), cases.pdb_id.nunique(),
             int(cases.get("cross_dock_ready", pd.Series(dtype=bool))
                 .fillna(False).sum()))
    log.info("tiers: %s", cases.tier.value_counts().to_dict())
    log.info("excluded covalent: %d", int((df.status == "excluded_covalent").sum()))
    log.info("reference inside the production box: %d of %d "
             "(median centroid offset %.1f A)",
             int(cases.ref_in_prod_box.sum()), len(cases),
             cases.ref_centroid_dist_to_prod_centre_a.median())
    if "ca_fit_rmsd_a" in cases:
        log.info("CA fit RMSD to 6VAJ: median %.2f A, max %.2f A",
                 cases.ca_fit_rmsd_a.median(), cases.ca_fit_rmsd_a.max())


if __name__ == "__main__":
    main()
