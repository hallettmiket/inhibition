"""
Purpose: Show docked poses IN THE POCKET against a labelled surface, and animate MD.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-31
Input: docked SDF/PDBQT per candidate, the prepared receptor, GROMACS trajectories
Output: py3Dmol HTML for embedding in the Streamlit dossier, plus an export bundle

WHY THE RECEPTOR IS NOT OPTIONAL (issue #1, T_4 note: "poses were given without
the pocket for some reason??"). A ligand rendered alone is a conformer, not a
pose. Every claim a docked pose makes -- that a warhead reaches Cys113, that a
substituent occupies a subpocket, that a molecule is too large -- is a claim
about the ligand RELATIVE TO the protein, and none of it is checkable without
the protein on screen. So the receptor is drawn first and the ligand into it;
there is no ligand-only path in this module.

THE RECEPTOR IS A LABELLED SURFACE, NOT A RIBBON (issue #3.1). A spectrum
cartoon tells a reader where the chain runs; it does not tell them where the
ligand IS. The question a docked pose exists to answer is "which sub-pocket does
this occupy", and that is a question about shape and about three named regions,
so the pocket is drawn as a surface with those three regions coloured and
labelled: the proline-binding pocket, the phosphate-binding Arg loop, and the
Cys113 pocket. Residue numbers are VERIFIED against the prepared receptor at
import-check time (`verify_subpockets`) rather than transcribed -- a mislabelled
pocket is worse than an unlabelled one because it is confidently wrong.

ALL POSES, NOT JUST THE FIRST (issue #3.1). Both engines write nine or more
modes per ligand and the viewer used to draw model 1 only. For Vina that is at
least the best-scoring mode. For gnina it is NOT: the SDF comes back in CNNscore
order, and on this build 11 of 25 T_3 shortlist entries and 6 of 27 T_4 entries
have a pose in the same file with a >0.5 kcal/mol BETTER minimizedAffinity than
the one that was drawn (max 2.86 kcal/mol, T_3). Showing one pose therefore hid
a real disagreement between the two gnina scores. `read_poses` returns every
model with its own scores and `pose_html` can draw any subset or overlay them.

THE COVALENT LINK IS DRAWN, NOT ASSUMED. T_3/T_4 dock the adduct with
`--covalent_rec_atom A:113:SG`, so the attachment atom sits ~1.81 A from the
Cys113 SG in every mode. That bond is the entire mechanistic claim of those two
approaches and it is invisible in a stick rendering of two separate models, so
it is drawn explicitly with its measured length.

THE MD ANIMATION IS THE PBC-CORRECTED TRAJECTORY, NOT prod.xtc. `whole.xtc` is
what `trjconv -pbc mol -center` produced. Animating the raw trajectory shows the
ligand teleporting across the box every time it crosses a periodic boundary --
which is what made the uncorrected RMSD 24x too large (D0038). A viewer reading
the raw file would see that artefact and reasonably conclude the ligand
dissociated.
"""

from __future__ import annotations

#: Mtime of THIS file at the moment it was imported. Frozen at import, so
#: comparing it with the file's current mtime is the only reliable way to tell
#: that a running process is executing stale code -- Streamlit re-runs the
#: script on every interaction but never re-imports helper modules.
LOADED_MTIME = __import__("os").stat(__file__).st_mtime

import io
import logging
import math
import os
import subprocess
import zipfile
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Sequence

RECEPTOR = Path("/data/lab_vm/immutable/inhibition/receptor/6VAJ_prepared.pdb")

# THE RECEPTOR MUST MATCH THE POSE'S COORDINATE FRAME, AND 6VAJ IS NOT UNIVERSAL.
#
# 6VAJ and the chemist-prepared 3IKD sit **48.6 A apart** in space: their Cys113
# sulfurs are at [-12.53, -35.87, 8.19] and [13.38, 3.99, -2.04]. A pose docked
# into one, drawn against the other, renders as a molecule floating in open
# solvent — which is exactly how the near-attack poses first appeared, and it
# reads as "the docking is broken" rather than "the viewer loaded the wrong
# protein".
#
# Keyed by the pose SOURCE, because the two coexist: the production poses were
# docked into 6VAJ and the near-attack poses into 3IKD (D0059), and both are
# valid against their own receptor.
RECEPTOR_3IKD = Path(
    "/data/lab_vm/modifiable/inhibition/receptor_3ikd_prep/3IKD_noligand.pdb")

RECEPTOR_FOR_POSE_COLUMN = {
    "pose_path": RECEPTOR,           # production docking, 6VAJ
    "nac_pose_path": RECEPTOR_3IKD,  # reactive near-attack docking, 3IKD
}


def receptor_for(pose_column: str = "pose_path") -> Path:
    """The receptor a pose from `pose_column` must be drawn against.

    Unknown columns fall back to the historical default rather than raising:
    a new pose source should render something recognisable while its receptor is
    wired up, and the frames themselves say which column was used.
    """
    return RECEPTOR_FOR_POSE_COLUMN.get(pose_column, RECEPTOR)
POCKET_RESIDUES = Path(
    "/data/lab_vm/immutable/inhibition/receptor/pocket_residues.json")
DATA = Path("/data/lab_vm/append_only/inhibition")
GMX = Path("/data/lab_vm/envs/dwi_gromacs_cuda/bin/gmx")

def receptor_readable() -> bool:
    """Can this process actually READ the prepared receptor?

    NOT `RECEPTOR.is_file()`. The data root is governed by an Isilon ACL the
    client cannot see, so a file can be present, `stat`-able and listed while
    every read of it raises `PermissionError` -- measured on this box for
    @tt8804 on 2026-08-04, where `is_file()` is True and `os.access(R_OK)` is
    False for every file under `immutable/inhibition/receptor/`.

    Guards written as "is it present" therefore do not fire, and the failure
    arrives later as an assertion about Arg-loop residues rather than as "you
    cannot read the receptor". Thirteen tests failed that way before this
    existed. Presence and readability are different questions; ask the one you
    actually depend on.
    """
    return os.access(RECEPTOR, os.R_OK)


# The protonation tag the docking run writes under. Imported from the module
# that OWNS it rather than restated, so the read path cannot drift from the
# write path -- restating it as a literal here is how they came apart before.
try:
    from shared.noncovalent_dock_run import LIGAND_PREP_TAG as _LIGAND_PREP_TAG
except Exception:  # noqa: BLE001 - the GUI must still start without the stack
    _LIGAND_PREP_TAG = "ph7.4"

DOCKING_DIRS = {
    "t1": DATA / "01_t1_de_novo" / "docking",
    "t2": DATA / "02_t2_atra_crem" / "docking",
    "t3": DATA / "03_t3_reinvent" / "docking",
    "t4": DATA / "04_t4_combinatorial" / "docking",
}
GROMACS_DIRS = {
    "t1": DATA / "01_t1_de_novo" / "gromacs",
    "t2": DATA / "02_t2_atra_crem" / "gromacs",
}

LIGAND_RESNAMES = ("MOL", "LIG", "UNL")
CATALYTIC_RESI = 113
CATALYTIC_CHAIN = "A"
#: The atom gnina is told to bond to (`config/receptor.yaml`:
#: `--covalent_rec_atom A:113:SG`). Named, not hardcoded as coordinates -- the
#: coordinates are read out of the receptor so they cannot drift from it.
CATALYTIC_ATOM = "SG"

#: A carbon-sulfur single bond is ~1.8 A. Anything inside this window is the
#: covalent link gnina was told to form; anything outside it is a contact and
#: must not be drawn as a bond, or a non-covalent pose acquires a bond it does
#: not have.
COVALENT_BOND_MAX_A = 2.4


# --------------------------------------------------------------------------
# the three sub-pockets
#
# WHY THESE ARE DECLARED WITH THEIR RESIDUE NAMES. A sub-pocket label is a claim
# about which residues line it, and the cheapest way to get it wrong is to
# transcribe a number from a paper that used a different construct's numbering.
# Each entry carries the three-letter name the residue MUST have in the prepared
# receptor, and `verify_subpockets` checks every one against the file. A test
# runs that check against the real receptor, so a wrong number fails the suite
# rather than mislabelling a figure a chemist then acts on.
#
# The numbering below was read out of 6VAJ chain A AND checked against the
# primary literature, because the structure alone cannot tell you what a region
# is conventionally called:
#
#   Ranganathan, Lu, Hunter & Noel, Cell 89:875 (1997), PMID 9200606
#     -- the founding structure. Names the basic triad Lys63/Arg68/Arg69 and the
#        hydrophobic proline pocket Leu122/Met130/Phe134.
#   Behrsin et al., J Mol Biol 365:1143 (2007), PMID 17113106
#     -- unigenic evolution, corroborates both sets functionally.
#   Dubiella et al., Nat Chem Biol 17:954 (2021), PMID 33972797
#     -- the sulfopin paper, i.e. the source of THIS receptor (6VAJ). Describes
#        the pocket in inhibitor-bound terms as Met130/Gln131/Phe134/Thr152/
#        His157.
#   Vohringer-Martinez, Verstraelen & Ayers, J Phys Chem B 118:9871 (2014),
#     PMID 25059768 -- Cys113's role in NATIVE catalysis is electrostatic
#     stabilisation, not nucleophilic attack. See the note on `cys113` below.
#
# Cross-checked against the reference ligand's own contacts in 6VAJ: QT7
# (sulfopin) touches Met130/Phe134/Leu122/Gln131 with its aryl end, Arg69/His59
# with its sulfonyl end, and its C10 sits 1.78 A from Cys113 SG.
#
# THE THREE REGIONS ARE KEPT DISJOINT. Dubiella's inhibitor-era description of
# the proline pocket also sweeps in Thr152 and His157, which the 1997/2007
# definition assigns to the catalytic site. Putting a residue in two coloured
# regions would make the surface unreadable, so each is claimed once, by the
# older and narrower definition, and the overlap is stated in `why` rather than
# resolved silently. Phe125 and Ser154 were considered and Phe125 dropped: no
# source names it as pocket-defining, and it sits 10.5 A from the SG.
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Subpocket:
    """One named region of the Pin1 active site, with its lining residues."""

    key: str
    label: str
    #: (residue number, expected three-letter name) in the prepared receptor.
    residues: tuple[tuple[int, str], ...]
    colour: str
    why: str

    @property
    def resi(self) -> list[int]:
        return [i for i, _ in self.residues]


SUBPOCKETS: tuple[Subpocket, ...] = (
    Subpocket(
        key="proline",
        label="proline-binding pocket",
        residues=((122, "LEU"), (130, "MET"), (131, "GLN"), (134, "PHE")),
        # Warm orange = hydrophobic, matching depict.WARHEAD_HIGHLIGHT's warmth
        # so the 2D and 3D panels do not use the same colour for two meanings.
        colour="orange",
        why=("The hydrophobic slot that accepts the prolyl ring of a "
             "pSer/pThr-Pro substrate. Leu122/Met130/Phe134 are the canonical "
             "lining (Ranganathan 1997; Behrsin 2007); Gln131 caps the rim and "
             "hydrogen-bonds sulfopin's sulfone at 2.9 A here. The sulfopin "
             "paper's wider description adds Thr152/His157, which this view "
             "colours with the catalytic site instead."),
    ),
    Subpocket(
        key="phosphate",
        label="phosphate-binding Arg loop",
        residues=((63, "LYS"), (68, "ARG"), (69, "ARG")),
        # Blue for basic, the near-universal convention for positive charge.
        colour="royalblue",
        why=("The basic triad on the catalytic loop (residues 63-80) that grips "
             "the substrate phosphate — the reason Pin1 is phospho-specific, "
             "and the reason an anionic group on a ligand has somewhere to go. "
             "Ser71 sits in the same loop but is a phosphosite, not part of the "
             "triad, so it is not coloured here."),
    ),
    Subpocket(
        key="cys113",
        label="Cys113 pocket (catalytic tetrad)",
        residues=((59, "HIS"), (113, "CYS"), (154, "SER"), (157, "HIS")),
        # Yellow = sulfur, and it is the colour Cys113's sticks were already
        # drawn in, so the region and the residue read as the same thing.
        colour="yellow",
        why=("The catalytic tetrad Cys113/His59/His157/Ser154. This is where "
             "T₃ and T₄ form their covalent bond and where sulfopin's "
             "chloroacetamide alkylates the SG (1.78 A in 6VAJ). NOTE: Cys113 "
             "is alkylated by covalent INHIBITORS, but in native catalysis the "
             "current reading is that it stabilises the cis conformer "
             "electrostatically rather than attacking as a nucleophile "
             "(Vohringer-Martinez 2014) — so 'catalytic cysteine' names the "
             "residue, not a covalent step in the enzyme's own mechanism."),
    ),
)

SUBPOCKETS_BY_KEY = {s.key: s for s in SUBPOCKETS}


class ReceptorMismatch(RuntimeError):
    """A declared sub-pocket residue is not what the receptor file says it is."""


@lru_cache(maxsize=2)
def _receptor_residues(pdb: str | None = None) -> dict[int, tuple[str, list]]:
    """{resi: (resname, [(atom_name, x, y, z), ...])} for the catalytic chain.

    Parsed by column offsets rather than by whitespace: `reduce` writes atom
    names that abut the residue name in some records, and a split() parse
    silently mis-columns exactly those.
    """
    path = Path(pdb) if pdb else RECEPTOR
    out: dict[int, tuple[str, list]] = {}
    if not path.is_file():
        return out
    for line in path.read_text(errors="ignore").splitlines():
        if not line.startswith(("ATOM", "HETATM")) or len(line) < 54:
            continue
        if line[21] != CATALYTIC_CHAIN:
            continue
        try:
            resi = int(line[22:26])
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        entry = out.setdefault(resi, (line[17:20].strip(), []))
        entry[1].append((line[12:16].strip(), *xyz))
    return out


def verify_subpockets(pdb: str | None = None) -> list[str]:
    """Every declared residue that the receptor does not agree with.

    Returns a list of human-readable mismatches; empty means every sub-pocket
    residue exists in the prepared receptor with the name declared here. This is
    the check that stops a sub-pocket label from being confidently wrong.
    """
    res = _receptor_residues(pdb)
    if not res:
        return [f"receptor not readable: {pdb or RECEPTOR}"]
    problems: list[str] = []
    for sp in SUBPOCKETS:
        for resi, expected in sp.residues:
            found = res.get(resi)
            if found is None:
                problems.append(
                    f"{sp.key}: residue {resi} absent from chain "
                    f"{CATALYTIC_CHAIN}")
            elif found[0] != expected:
                problems.append(
                    f"{sp.key}: residue {resi} is {found[0]}, declared "
                    f"{expected}")
    return problems


def subpocket_centroid(key: str, pdb: str | None = None) -> tuple[float, float, float] | None:
    """Centre of a sub-pocket's side-chain atoms, for placing its label.

    Side chains only (backbone N/CA/C/O dropped): the backbone of a loop runs
    away from the cavity it lines, so a whole-residue centroid puts the label
    behind the surface instead of on it.
    """
    sp = SUBPOCKETS_BY_KEY.get(key)
    res = _receptor_residues(pdb)
    if sp is None or not res:
        return None
    pts = [(x, y, z)
           for resi, _ in sp.residues
           for name, x, y, z in res.get(resi, ("", []))[1]
           if name not in ("N", "CA", "C", "O")]
    if not pts:
        return None
    n = len(pts)
    return (sum(p[0] for p in pts) / n,
            sum(p[1] for p in pts) / n,
            sum(p[2] for p in pts) / n)


def catalytic_sg(pdb: str | None = None) -> tuple[float, float, float] | None:
    """Coordinates of the Cys113 SG that gnina bonds to, read from the receptor."""
    entry = _receptor_residues(pdb).get(CATALYTIC_RESI)
    if entry is None:
        return None
    for name, x, y, z in entry[1]:
        if name == CATALYTIC_ATOM:
            return (x, y, z)
    return None


@lru_cache(maxsize=2)
def pocket_resi(pdb: str | None = None) -> tuple[int, ...]:
    """Every residue lining the site, for the grey surface the labels sit on.

    Read from `pocket_residues.json` -- the 8 A shell around the reference
    ligand that receptor_prep wrote -- and falls back to the union of the three
    sub-pockets when that file is unavailable. The earlier version surfaced
    `resi 101..125`, a SEQUENCE window around Cys113, which is not the same set
    as the residues that are spatially near it: it included residues pointing
    away from the site and excluded the whole Arg loop.
    """
    declared = {i for sp in SUBPOCKETS for i in sp.resi}
    try:
        import json
        raw = json.loads(POCKET_RESIDUES.read_text(encoding="utf-8"))
        listed = {int(str(r).split(":")[-1]) for r in raw.get("resi_list", [])}
        return tuple(sorted(listed | declared)) if listed else tuple(sorted(declared))
    except Exception as exc:  # noqa: BLE001 - a missing shell file must not blank the view
        # MISSING AND UNREADABLE ARE DIFFERENT FACTS. The fallback is right for
        # both -- a blank surface is worse than a small one -- but only the
        # first is benign. When the file is PRESENT and we cannot read it, the
        # view silently shrinks from the measured 8 A shell (>25 residues) to
        # the 11 declared sub-pocket residues, and looks like a normal picture.
        # That is the shape catalogued in `how_this_project_breaks.md`: a
        # populated, plausible value computed from the wrong thing. Say so.
        if POCKET_RESIDUES.exists():
            logging.getLogger(__name__).warning(
                "pocket shell %s exists but could not be read (%s); the grey "
                "surface falls back to the %d declared sub-pocket residues and "
                "is SMALLER than the measured 8 A shell",
                POCKET_RESIDUES.name, type(exc).__name__, len(declared))
        return tuple(sorted(declared))


# --------------------------------------------------------------------------
# finding files
# --------------------------------------------------------------------------

def find_pose(approach: str, candidate_id: str,
              dock_id: str | None = None,
              pose_path: str | Path | None = None) -> Path | None:
    """The docked pose file for one candidate, or None.

    `pose_path` IS THE FRAME'S OWN ANSWER AND IS TRIED FIRST. T_3 and T_4 record
    the exact file they docked into; deriving it again from ids is a second
    source of truth that can only ever agree or be wrong.

    THE FILENAME CANNOT BE DERIVED FROM THE CANDIDATE ID, AND GUESSING IT FAILED
    SILENTLY. The covalent approaches name pose files by a SEPARATE `dock_id`
    (`d_31fe9e96d8bb_docked.sdf`) that shares no hash with `candidate_id`
    (`t3_31d6edc305b0`) — the overlap between the two sets is exactly zero
    across 4080 files. An earlier version built the name from `candidate_id`,
    found nothing, and rendered "no docked pose file found", which reads as
    missing data rather than as a lookup bug. Pass `dock_id` from the frame.

    THE TWO STRATA ALSO USE DIFFERENT LAYOUTS AND FORMATS:
      T_3, T_4 (gnina, covalent)  docking/d_<dock_id>_docked.sdf
      T_1, T_2 (Vina, non-covalent)  docking/poses/<candidate_id>_out.pdbqt
    Both are returned; `read_poses` dispatches on the suffix.
    """
    if pose_path:
        p = Path(pose_path)
        if p.is_file() and p.stat().st_size > 0:
            return p

    d = DOCKING_DIRS.get(approach)
    if d is None or not d.is_dir():
        return None

    # Non-covalent: Vina writes PDBQT into a poses/ subdirectory.
    #
    # THE PREP-TAGGED DIRECTORY FIRST, AND IT IS NOT COSMETIC. Ligands are
    # docked protonated for pH 7.4 and land in `poses_ph7.4/`; the untagged
    # `poses/` holds the SUPERSEDED neutral-form run from before 2026-07-31.
    # Both are populated, both parse, and both contain a file with exactly the
    # right name -- so preferring `poses/` returned a pose that was NOT the one
    # the frame's `vina_affinity` was computed from, and the viewer drew it
    # beside that score as though they matched. Measured 2026-08-04: every T_2
    # candidate had both, and the stale one was being shown.
    #
    # This is the same defect `LIGAND_PREP_TAG` was introduced to prevent, one
    # layer up: the write path was correctly tagged, the READ path was not.
    for sub in (f"poses_{_LIGAND_PREP_TAG}", "poses"):
        poses = d / sub
        if not poses.is_dir():
            continue
        for name in (f"{candidate_id}_out.pdbqt", f"{candidate_id}.pdbqt"):
            p = poses / name
            if p.is_file() and p.stat().st_size > 0:
                if sub == "poses" and (d / f"poses_{_LIGAND_PREP_TAG}").is_dir():
                    logging.getLogger(__name__).warning(
                        "%s: falling back to the untagged pose set for %s while "
                        "poses_%s exists — the pose shown may not be the one "
                        "the score came from", approach, candidate_id,
                        _LIGAND_PREP_TAG)
                return p

    # Covalent: gnina writes SDF named by dock_id.
    if dock_id:
        for name in (f"{dock_id}_docked.sdf", f"{dock_id}.sdf"):
            p = d / name
            if p.is_file() and p.stat().st_size > 0:
                return p
    return None


def find_trajectory(approach: str, candidate_id: str,
                    replicate: int = 1) -> tuple[Path, Path] | None:
    """(PBC-corrected trajectory, its tpr) for one replicate, or None."""
    root = GROMACS_DIRS.get(approach)
    if root is None:
        return None
    wd = root / candidate_id
    for rep_dir in (wd / f"rep{replicate}", wd):
        whole, tpr = rep_dir / "whole.xtc", rep_dir / "prod.tpr"
        if whole.is_file() and tpr.is_file():
            return whole, tpr
    return None


# --------------------------------------------------------------------------
# reading every pose out of a multi-model file
# --------------------------------------------------------------------------

#: Which score in each format is the affinity, and which way is better. Named
#: here so the viewer never has to guess: `minimizedAffinity` and
#: `vina_affinity` are kcal/mol and lower is better, `CNNaffinity` and
#: `CNNscore` are dimensionless and HIGHER is better. Mixing the directions is
#: the single easiest way to present the worst pose as the best.
SCORE_DIRECTION = {
    "vina_affinity": "lower",
    "minimizedAffinity": "lower",
    "CNNaffinity": "higher",
    "CNNscore": "higher",
    "CNN_VS": "higher",
}

#: The score each format is ORDERED by in the file, which is not necessarily the
#: score the frame ranks on. Vina sorts by affinity, so model 1 is the best
#: affinity. gnina sorts by CNNscore, so record 1 can be — and on this build
#: often is — a worse minimizedAffinity than a later record.
FILE_ORDER_SCORE = {"pdbqt": "vina_affinity", "sdf": "CNNscore"}


@dataclass(frozen=True)
class Pose:
    """One binding mode out of a multi-model docking output."""

    index: int                       # 1-based, in file order
    fmt: str                         # "pdb" (from PDBQT) or "sdf"
    text: str                        # a complete single-model block
    scores: dict[str, float] = field(default_factory=dict)

    def score(self, name: str) -> float | None:
        return self.scores.get(name)


def _pdbqt_poses(text: str) -> list[Pose]:
    """Split a Vina PDBQT into its MODEL blocks with their REMARK scores."""
    poses: list[Pose] = []
    cur: list[str] | None = None
    scores: dict[str, float] = {}
    for line in text.splitlines():
        if line.startswith("MODEL"):
            cur, scores = [line], {}
            continue
        if line.startswith("REMARK VINA RESULT") and cur is not None:
            parts = line.split(":", 1)[1].split()
            for key, val in zip(("vina_affinity", "rmsd_lb", "rmsd_ub"), parts):
                try:
                    scores[key] = float(val)
                except ValueError:
                    pass
        if cur is not None:
            cur.append(line)
        if line.startswith("ENDMDL") and cur is not None:
            poses.append(Pose(len(poses) + 1, "pdb", "\n".join(cur) + "\n",
                              dict(scores)))
            cur = None
    if not poses and text.strip():
        # A single-model file with no MODEL/ENDMDL wrapper is still one pose.
        poses.append(Pose(1, "pdb", text, {}))
    return poses


def _sdf_poses(text: str) -> list[Pose]:
    """Split a multi-record SDF into records with their property tags.

    ACCUMULATED LINE BY LINE RATHER THAN `text.split("$$$$")`. An SDF record
    opens with three header lines that gnina leaves BLANK, so a split-then-strip
    parse eats the title line and shifts the counts line into its place; 3Dmol
    then reads the atom count out of the wrong row and draws nothing.
    """
    poses: list[Pose] = []
    cur: list[str] = []
    scores: dict[str, float] = {}
    tag: str | None = None
    for line in text.splitlines():
        if line.startswith("$$$$"):
            if cur:
                poses.append(Pose(len(poses) + 1, "sdf",
                                  "\n".join(cur) + "\n$$$$\n", dict(scores)))
            cur, scores, tag = [], {}, None
            continue
        cur.append(line)
        stripped = line.strip()
        if stripped.startswith(">") and "<" in stripped and ">" in stripped[1:]:
            tag = stripped[stripped.index("<") + 1:stripped.rindex(">")]
            continue
        if tag is not None and stripped:
            try:
                scores[tag] = float(stripped)
            except ValueError:
                pass
            tag = None
    if cur and any(ln.strip() for ln in cur):
        poses.append(Pose(len(poses) + 1, "sdf", "\n".join(cur) + "\n", dict(scores)))
    return poses


def read_poses(path: Path) -> list[Pose]:
    """Every binding mode in a docked-pose file, in file order.

    Text parsing rather than RDKit on purpose: these are covalent ADDUCTS with
    an open valence at the attachment carbon, so RDKit needs `sanitize=False` to
    read them at all, and a sanitisation failure would drop poses silently. The
    viewer only needs the block and its scores, and both are in the text.
    """
    if not path.is_file():
        return []
    text = path.read_text(errors="ignore")
    if path.suffix.lower() in (".pdbqt", ".pdb"):
        return _pdbqt_poses(text)
    return _sdf_poses(text)


def pose_atoms(pose: Pose) -> list[tuple[str, float, float, float]]:
    """(element, x, y, z) for one pose, for measuring the covalent link."""
    out: list[tuple[str, float, float, float]] = []
    if pose.fmt == "pdb":
        for line in pose.text.splitlines():
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
                try:
                    out.append((line[76:78].strip() or line[12:16].strip()[:1],
                                float(line[30:38]), float(line[38:46]),
                                float(line[46:54])))
                except ValueError:
                    continue
        return out
    lines = pose.text.splitlines()
    if len(lines) < 4:
        return out
    try:
        n_atoms = int(lines[3][0:3])
    except (ValueError, IndexError):
        return out
    for line in lines[4:4 + n_atoms]:
        try:
            out.append((line[31:34].strip(), float(line[0:10]),
                        float(line[10:20]), float(line[20:30])))
        except (ValueError, IndexError):
            continue
    return out


def covalent_link(pose: Pose, pdb: str | None = None
                  ) -> tuple[tuple[float, float, float], tuple[float, float, float], float] | None:
    """(SG coords, nearest ligand atom coords, distance) if this pose is bonded.

    Returns None when the nearest ligand atom is further than a bond away, which
    is the correct answer for T_1/T_2 — they are non-covalent and drawing a bond
    into them would assert a mechanism they do not have.
    """
    sg = catalytic_sg(pdb)
    if sg is None:
        return None
    atoms = pose_atoms(pose)
    if not atoms:
        return None
    best = min(atoms, key=lambda a: math.dist((a[1], a[2], a[3]), sg))
    d = math.dist((best[1], best[2], best[3]), sg)
    if d > COVALENT_BOND_MAX_A:
        return None
    return sg, (best[1], best[2], best[3]), d


def pose_score_table(poses: Sequence[Pose]) -> list[dict]:
    """One row per pose, ready for `st.dataframe`, with the best of each marked.

    The point of the table is the DISAGREEMENT between the scores. On this build
    gnina's file order (CNNscore) and its own minimizedAffinity pick a different
    best pose for 11 of 25 T_3 shortlist entries, so a reader who sees only pose
    1 is seeing the CNN's choice presented as the affinity's.
    """
    if not poses:
        return []
    names = [n for n in SCORE_DIRECTION if any(n in p.scores for p in poses)]
    best_by = {}
    for n in names:
        vals = [(p.scores[n], p.index) for p in poses if n in p.scores]
        if not vals:
            continue
        best_by[n] = (min(vals) if SCORE_DIRECTION[n] == "lower" else max(vals))[1]
    rows = []
    for p in poses:
        row: dict = {"pose": p.index}
        for n in names:
            v = p.scores.get(n)
            row[n] = (f"{v:.3f} ★" if best_by.get(n) == p.index else
                      (f"{v:.3f}" if v is not None else None))
        rows.append(row)
    return rows


def hidden_better_pose(poses: Sequence[Pose], shown: int = 1,
                       score: str = "minimizedAffinity") -> tuple[int, float] | None:
    """(pose index, improvement) if a pose beats the shown one, else None.

    Exists so the caption can state the size of what one-pose rendering was
    hiding rather than merely offering more poses and leaving the reader to
    notice.
    """
    vals = [(p.index, p.scores[score]) for p in poses if score in p.scores]
    if not vals or all(i != shown for i, _ in vals):
        return None
    shown_val = next(v for i, v in vals if i == shown)
    better = (min(vals, key=lambda t: t[1]) if SCORE_DIRECTION.get(score) == "lower"
              else max(vals, key=lambda t: t[1]))
    gap = abs(better[1] - shown_val)
    return (better[0], gap) if better[0] != shown and gap > 0 else None


#: How to drive a 3Dmol.js canvas. Written out because the viewer offers no
#: on-screen affordances whatsoever -- no buttons, and nothing indicates that
#: right-drag zooms or Ctrl+drag pans. A reader who tries only left-drag
#: concludes the view is stuck off-centre and unusable.
CONTROLS = """
| action | mouse | trackpad |
|---|---|---|
| **rotate** | left-click + drag | one-finger drag |
| **zoom** | scroll wheel, or right-click + drag | two-finger scroll / pinch |
| **pan** — shift the molecule sideways | middle-click + drag, or **Ctrl** + left-drag | **Ctrl** + one-finger drag |
| **slab** — cut away the front | **Ctrl** + right-drag | — |
| **reset** | change the *framing* control and back | same |
"""

#: Colours for overlaid poses, in draw order. The first is the cyan the single
#: pose has always been, so a reader who has seen this viewer before still
#: recognises "the pose"; the rest are chosen to stay distinguishable against
#: the orange/blue/yellow sub-pocket surfaces.
POSE_COLOURS = ("cyanCarbon", "magentaCarbon", "greenCarbon", "whiteCarbon",
                "purpleCarbon", "salmonCarbon", "blueCarbon", "greyCarbon",
                "brownCarbon")


def pose_html(pose_file: Path, *, show: Sequence[int] = (1,),
              width: int | None = None, height: int = 620,
              surface: bool = True, label_subpockets: bool = True,
              show_covalent: bool = True, cartoon: bool = False,
              zoom_on: str = "ligand",
              receptor: Path | None = None) -> str:
    """Docked pose(s) inside the receptor, on a labelled sub-pocket surface.

    `show` is a list of 1-based pose indices. Passing more than one overlays
    them, each in its own colour, which is how a reader sees whether the nine
    modes agree about where the ligand sits or scatter across the site — the
    question a single pose cannot answer.

    `width=None` makes the canvas fill its container. This matters more than it
    sounds: the viewer was fixed at 700 px inside a one-of-four column roughly
    300 px wide, so most of the scene sat off-canvas and no amount of mouse work
    brought it back. The complaint was that the image could not be centred; the
    cause was that the canvas was wider than the space it was drawn into.

    `zoom_on` sets the opening framing -- "ligand" fills the view with the
    binding site, "pocket" backs off to show the labelled sub-pockets, "all"
    shows the whole protein for orientation.
    """
    import py3Dmol

    poses = read_poses(Path(pose_file))
    wanted = [p for p in poses if p.index in set(show)] or poses[:1]

    v = py3Dmol.view(width=width or "100%", height=height)
    rec = receptor or RECEPTOR
    v.addModel(rec.read_text(), "pdb")

    # THE SURFACE IS THE RECEPTOR REPRESENTATION, the cartoon is optional
    # context. A ribbon answers "where does the chain run"; the question a
    # docked pose asks is "what shape is it sitting in", and only a surface
    # answers that.
    # An EMPTY style, not a transparent one: `{"line": {"opacity": 0}}` still
    # builds and uploads every bond in the protein for nothing.
    v.setStyle({"model": 0}, {"cartoon": {"color": "spectrum", "opacity": 0.25}}
               if cartoon else {})

    lining = list(pocket_resi())
    claimed = {i for sp in SUBPOCKETS for i in sp.resi}
    if surface:
        # The rest of the site in grey FIRST, so the three coloured patches are
        # drawn over it and read as regions of one surface rather than as three
        # floating blobs.
        rest = [i for i in lining if i not in claimed]
        if rest:
            v.addSurface("VDW", {"opacity": 0.72, "color": "lightgrey"},
                         {"model": 0, "resi": rest})
        for sp in SUBPOCKETS:
            v.addSurface("VDW", {"opacity": 0.78, "color": sp.colour},
                         {"model": 0, "resi": sp.resi})

    # Cys113 as sticks on top of its own surface patch: the SG is the atom the
    # covalent approaches attack and it has to be locatable by eye, not merely
    # inside a coloured region.
    v.addStyle({"model": 0, "resi": CATALYTIC_RESI},
               {"stick": {"colorscheme": "yellowCarbon", "radius": 0.22}})

    if label_subpockets:
        for sp in SUBPOCKETS:
            c = subpocket_centroid(sp.key)
            if c is None:
                continue
            v.addLabel(sp.label,
                       {"position": {"x": c[0], "y": c[1], "z": c[2]},
                        "backgroundColor": sp.colour,
                        "backgroundOpacity": 0.75,
                        "fontColor": "black", "fontSize": 12,
                        "borderThickness": 0.5, "inFront": True})

    for n, pose in enumerate(wanted):
        v.addModel(pose.text, pose.fmt)
        model_i = n + 1
        # The FIRST requested pose is drawn thickest. With nine overlaid modes
        # the reader needs to know which one the frame's score refers to, and
        # colour alone does not carry that.
        v.setStyle({"model": model_i},
                   {"stick": {"colorscheme": POSE_COLOURS[n % len(POSE_COLOURS)],
                              "radius": 0.16 if n == 0 else 0.10}})
        if show_covalent:
            link = covalent_link(pose)
            if link is not None:
                sg, lig, dist = link
                v.addCylinder({
                    "start": {"x": sg[0], "y": sg[1], "z": sg[2]},
                    "end": {"x": lig[0], "y": lig[1], "z": lig[2]},
                    "radius": 0.09, "color": "yellow", "dashed": True,
                    "fromCap": 1, "toCap": 1})
                if n == 0:
                    mid = [(a + b) / 2 for a, b in zip(sg, lig)]
                    v.addLabel(f"covalent Cys113 SG · {dist:.2f} Å",
                               {"position": {"x": mid[0], "y": mid[1], "z": mid[2]},
                                "backgroundColor": "yellow",
                                "backgroundOpacity": 0.8, "fontColor": "black",
                                "fontSize": 11, "inFront": True})

    if zoom_on == "all":
        v.zoomTo()
    elif zoom_on == "pocket":
        v.zoomTo({"model": 0, "resi": lining})
    else:
        v.zoomTo({"model": 1})
        # Back off slightly: zoomTo on a small ligand crops the pocket walls
        # out of frame, which is the context the pose exists to show.
        v.zoom(0.7)
    return v._make_html()


def subpocket_legend() -> str:
    """Markdown legend naming the three regions and what each one is for."""
    rows = ["| region | colour | residues | what it does |",
            "|---|---|---|---|"]
    for sp in SUBPOCKETS:
        resi = ", ".join(f"{name.title()}{num}" for num, name in sp.residues)
        rows.append(f"| **{sp.label}** | {sp.colour} | {resi} | {sp.why} |")
    return "\n".join(rows)


# --------------------------------------------------------------------------
# hand-off to external software (issue #3.1c)
#
# WHY OFFER AN EXPORT AT ALL. The embedded viewer is a reading aid, not a
# workbench: it cannot measure a distance, cannot mutate a residue, cannot
# render for a figure, and silently gives up above a few MB. Anything a chemist
# wants to DO with a pose happens in PyMOL or ChimeraX, and the alternative to
# an export button is that they go and find the file by hand -- which for T_3/T_4
# means knowing that pose files are named by `dock_id` and not by candidate id,
# the exact lookup that has already caught this project out once.
#
# THE SCRIPT IS GENERATED FROM `SUBPOCKETS`, NOT WRITTEN OUT. If the exported
# session named its own residues, the GUI and the session could disagree about
# where the proline pocket is, and the export exists precisely to be trusted
# more than the embedded view.
#
# NOTHING IS WRITTEN TO DISK. The bundle is built in memory and handed to the
# browser: the GUI reads, it does not own (D0008), and the pose files live under
# the append-only root.
# --------------------------------------------------------------------------

def pymol_script(pose_name: str, *, receptor_name: str = "6VAJ_prepared.pdb",
                 covalent: bool = False) -> str:
    """A .pml that reproduces this view in PyMOL, sub-pockets and all."""
    lines = [
        "# Pin1 docked poses — generated by the Dance with Inhibition GUI.",
        "# Sub-pocket definitions come from integration/app/pose3d.py:SUBPOCKETS,",
        "# so this session and the GUI cannot disagree about where they are.",
        "",
        f"load {receptor_name}, receptor",
        f"load {pose_name}, poses",
        "",
        "hide everything",
        "show surface, receptor",
        "set transparency, 0.25",
        "color grey80, receptor",
        "show sticks, poses",
        "color cyan, poses and elem C",
        "set all_states, on          # every docked mode at once",
        "",
    ]
    for sp in SUBPOCKETS:
        sel = "+".join(str(i) for i in sp.resi)
        lines += [
            f"# {sp.label}: {sp.why}",
            f"select {sp.key}, receptor and resi {sel}",
            f"color {sp.colour}, {sp.key}",
            f"label {sp.key} and name CA and resi {sp.resi[0]}, \"{sp.label}\"",
            "",
        ]
    lines += [
        f"select cys113_sg, receptor and resi {CATALYTIC_RESI} and name {CATALYTIC_ATOM}",
        "show spheres, cys113_sg",
        "color yellow, cys113_sg",
    ]
    if covalent:
        lines += [
            "",
            "# The covalent link gnina was told to form "
            f"(--covalent_rec_atom A:{CATALYTIC_RESI}:{CATALYTIC_ATOM}).",
            "# Measured, not drawn as a bond: the adduct and the receptor are "
            "separate objects here.",
            f"distance covalent_link, cys113_sg, poses within "
            f"{COVALENT_BOND_MAX_A} of cys113_sg",
            "color yellow, covalent_link",
        ]
    lines += ["", "orient poses", "zoom poses, 6", "deselect"]
    return "\n".join(lines) + "\n"


def pose_bundle(pose_file: Path, *, covalent: bool = False,
                include_receptor: bool = True) -> bytes:
    """A zip of {pose file, receptor, view.pml, README} for external software.

    Built entirely in memory — the GUI never writes to the governed tree.
    """
    pose_file = Path(pose_file)
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr(pose_file.name, pose_file.read_bytes())
        if include_receptor and RECEPTOR.is_file():
            z.writestr(RECEPTOR.name, RECEPTOR.read_bytes())
        z.writestr("view_pin1.pml",
                   pymol_script(pose_file.name, receptor_name=RECEPTOR.name,
                                covalent=covalent))
        z.writestr("README.txt", _bundle_readme(pose_file, covalent=covalent))
    return buf.getvalue()


def _bundle_readme(pose_file: Path, *, covalent: bool) -> str:
    n = len(read_poses(pose_file))
    order = FILE_ORDER_SCORE.get(pose_file.suffix.lower().lstrip("."), "?")
    return (
        f"Docked poses for one candidate — {pose_file.name}\n"
        f"{'=' * 60}\n\n"
        f"Contents:\n"
        f"  {pose_file.name:<28} {n} binding mode(s), in FILE order\n"
        f"  {RECEPTOR.name:<28} the prepared receptor all four approaches used\n"
        f"  view_pin1.pml{'':<16} PyMOL session: surface + the three sub-pockets\n\n"
        f"Open with:  pymol view_pin1.pml\n\n"
        f"FILE ORDER IS BY {order}, WHICH MAY NOT BE THE AFFINITY ORDER.\n"
        + ("gnina returns modes sorted by CNNscore. The candidate frame's\n"
           "`affinity_kcal` is taken from the FIRST record, so a later record\n"
           "in this file can have a better minimizedAffinity — on the current\n"
           "build that happens for 11 of 25 T₃ shortlist entries, by up to\n"
           "2.86 kcal/mol. Check every mode before quoting one.\n\n"
           if pose_file.suffix.lower() == ".sdf" else
           "AutoDock Vina returns modes sorted by affinity, so model 1 is the\n"
           "best-scoring mode in this file.\n\n")
        + ("This candidate is a COVALENT adduct: it was docked with\n"
           f"--covalent_rec_atom A:{CATALYTIC_RESI}:{CATALYTIC_ATOM}, and the\n"
           "attachment atom sits ~1.8 A from the Cys113 SG in every mode. The\n"
           "SMILES in this file is the POST-reaction species, not the molecule\n"
           "as synthesised (D0022/D0030).\n\n" if covalent else "")
        + "No score here is evidence of binding: the non-covalent enrichment\n"
          "gate returned WEAK (ROC-AUC 0.599, CI [0.311, 0.874], EF1% 0.0;\n"
          "D0041) and the covalent gate is UNDERPOWERED.\n")


#: Streamlit embeds the viewer HTML inline. Past roughly this size the browser
#: renders nothing at all and the click looks like a dead button -- which is
#: exactly what a 12 MB movie did. Guarded rather than hoped about.
MAX_EMBED_BYTES = 4_000_000


def trajectory_html(xtc: Path, tpr: Path, *, n_frames: int = 15,
                    width: int | None = None, height: int = 620) -> str | None:
    """An animated MD trajectory, or None if it cannot be built or is too big.

    WHY SO FEW FRAMES. `whole.xtc` already contains ONLY the solute -- the
    PBC-correction step wrote protein + ligand and dropped ~27,600 waters and
    ions. What remains is 2,391 atoms, and at 60 frames that is 12 MB of inline
    HTML. Streamlit embeds the viewer directly in the page and the browser
    silently renders nothing above a few MB, so the button appeared dead. Frame
    count is the only lever that does not require re-deriving an index.

    A BACKBONE SELECTION WAS TRIED AND IS WRONG HERE. `gmx select` indexes
    against `prod.tpr` (30,037 atoms) while `whole.xtc` holds 2,391, so every
    index past the solute overruns the trajectory and trjconv aborts with a
    mismatch. Any future selection must be built against the solute subset, not
    the full topology.
    """
    ndx = xtc.parent / "analysis.ndx"
    if not ndx.is_file():
        return None
    out = xtc.parent / f"movie_{n_frames}.pdb"
    if not out.is_file():
        try:
            subprocess.run(
                [str(GMX), "trjconv", "-s", str(tpr), "-f", str(xtc),
                 "-n", str(ndx), "-o", str(out), "-skip",
                 str(max(1, 1000 // n_frames))],
                input="0\n", capture_output=True, text=True, timeout=900,
                cwd=xtc.parent,
                env={"GMX_MAXBACKUP": "-1", "PATH": "/usr/bin:/bin"}, check=True)
        except Exception:  # noqa: BLE001 - a missing movie is not a page failure
            return None
    if not out.is_file() or out.stat().st_size == 0:
        return None
    if out.stat().st_size > MAX_EMBED_BYTES:
        log_size = out.stat().st_size / 1e6
        raise MovieTooLarge(
            f"{out.name} is {log_size:.1f} MB; above ~{MAX_EMBED_BYTES/1e6:.0f} MB "
            "the browser renders nothing and the control looks broken. Reduce "
            "n_frames.")

    import py3Dmol

    text = out.read_text()
    # SELECT THE LIGAND BY RESIDUE NAME, NOT BY hetflag. `gmx trjconv` writes
    # every atom as an ATOM record -- there is not one HETATM in the file -- so
    # {"hetflag": True} matches nothing. The ligand was therefore never styled
    # as sticks AND zoomTo had an empty selection, which is why the animation
    # opened somewhere arbitrary and would not centre on the molecule.
    resn = next((r for r in LIGAND_RESNAMES
                 if f" {r} " in text or f"\n{r}" in text
                 or any(line[17:20].strip() == r
                        for line in text.splitlines()[:4000]
                        if line.startswith(("ATOM", "HETATM")))), None)
    lig_sel = {"resn": resn} if resn else {"hetflag": True}

    v = py3Dmol.view(width=width or "100%", height=height)
    v.addModelsAsFrames(text, "pdb")
    v.setStyle({"cartoon": {"color": "spectrum", "opacity": 0.55}})
    v.addStyle(lig_sel, {"stick": {"colorscheme": "cyanCarbon", "radius": 0.22}})
    v.addStyle({"resi": CATALYTIC_RESI},
               {"stick": {"colorscheme": "yellowCarbon", "radius": 0.2}})
    # Centre on the ligand, then back off so the pocket it sits in is visible.
    # zoomTo on a ~40-atom ligand alone fills the canvas with the ligand and
    # crops away the protein it is supposed to be staying inside.
    v.zoomTo(lig_sel)
    v.zoom(0.6)
    v.animate({"loop": "forward", "interval": 80})
    return v._make_html()


class MovieTooLarge(RuntimeError):
    """The rendered movie would exceed what a browser will embed."""


def rmsd_series(approach: str, candidate_id: str) -> dict[int, list[tuple]]:
    """Per-replicate (time_ns, ligand RMSD nm) from each rep's rmsd.xvg.

    Returned per replicate rather than averaged. D0038's whole lesson was that
    one trajectory cannot separate "this ligand leaves" from "this trajectory
    wandered", so the spread between replicates IS the result and collapsing it
    to a mean would hide exactly what the replicates were run to show.
    """
    root = GROMACS_DIRS.get(approach)
    if root is None:
        return {}
    wd = root / candidate_id
    if not wd.is_dir():
        return {}
    out: dict[int, list[tuple]] = {}
    for rep_dir in sorted(wd.glob("rep*")):
        xvg = rep_dir / "rmsd.xvg"
        if not xvg.is_file():
            continue
        rows = []
        for line in xvg.read_text(errors="ignore").splitlines():
            line = line.strip()
            if not line or line[0] in "@#":
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    rows.append((float(parts[0]), float(parts[1])))
                except ValueError:
                    continue
        if rows:
            try:
                out[int(rep_dir.name[3:])] = rows
            except ValueError:
                continue
    return out
