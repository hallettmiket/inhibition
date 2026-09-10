"""
Purpose: a ligand-based pharmacophore scorer, orthogonal to docking, that can be gated.
Author: @tt8804 (with Claude Code)
Date: 2026-09-09
Input: the frozen Pin1 reference binder set (resolved by glob) + candidate SMILES
Output: `ph2d_max_similarity` per candidate, and the model set it was scored against

WHY THIS EXISTS. Every level of theory this project has measured needs either
the receptor (docking, MM-GBSA, MD) or a GPU, and #4 Phase 2.2 asks for an
orthogonal, GATE-TESTABLE scorer. A 2D pharmacophore fingerprint needs neither:
no receptor, no pose, no scoring function, no GPU. That makes it the one scorer
that can be measured against AID 504891's 361,354 ASSAYED inactives on CPU
alone, which is Phase 2.1's harness pointed at a scorer that can actually run
today.

WHAT IT IS, STATED PLAINLY. This is pharmacophore SIMILARITY to known actives,
not a match against a pharmacophore HYPOTHESIS. It asks "does this molecule
carry the same arrangement of donors, acceptors, aromatics and charges as a
known Pin1 binder", by Gobbi 2D feature-pair fingerprint. It does NOT ask
"does this molecule place a feature where the pocket needs one". The second
question needs the co-crystals and 3D conformers and is deliberately deferred;
do not let this module's output be described as the second thing.

THE MODEL SET IS THE WHOLE DESIGN, AND IT HAS TWO LEAKS TO CLOSE.

1.  THE SEED. T_2 pools are CReM derivatives of a seed, and both degree-2
    seeds -- ATRA and Guo-Pfizer -- are THEMSELVES in the reference binder
    set. Scoring Guo's derivatives against a model containing Guo asks "how
    much do you still look like your parent", which is circular in exactly the
    way `reference_set`'s control B4 exists to prevent for the novelty axis.
    Measured 2026-09-09 on 2,000 molecules per pool:

        pool             to own seed   best to any other binder   seed nearest
        ATRA degree-2        0.053              0.107                16.8%
        Guo  degree-2        0.291              0.269                58.8%

    It bites for Guo and not for ATRA -- a polyene carries almost no
    pharmacophore features, so its products drift away from it fast. The
    asymmetry is the point: a rule applied only where it looked necessary
    would have been applied to neither. `hold_out` is therefore MANDATORY at
    the call site, not a default, because the caller is the only one that
    knows which pool it is scoring.

2.  THE PEPTIDES. The reference set's `lead` tier includes four peptides and
    cell-penetrating peptides carrying 28 to 103 pharmacophore features. A
    ~30-heavy-atom CReM product cannot approach them, so they contribute a
    ceiling nothing can reach and no signal at all.

THE MODEL-SET RULE, WRITTEN DOWN BEFORE IT WAS APPLIED. D0045 fixed its
chemotype definition before taking counts precisely so the answer could not be
chosen; the same discipline applies here, because "drug-like" is exactly the
kind of word that gets retrofitted to whatever produced the nicer number. A
molecule enters the model set iff:

    (i)   its SMILES parses -- which excludes the `UNVERIFIED` placeholder
          still sitting in the Byun-BDHI-fragment row;
    (ii)  it carries FEWER THAN 4 backbone amide bonds; and
    (iii) it has AT MOST 45 heavy atoms.

(ii) is the classical peptide criterion and does the real work. (iii) is a
backstop for a large non-peptidic macrocycle. Liu-2024-C3, the most potent
non-covalent binder in the set, sits at 41 heavy atoms and 2 amides and is
retained -- the thresholds were chosen to sit clear of it on both axes rather
than to trim it.

HIGHER IS BETTER, AND THE REGISTRY NOW SAYS SO. `ph2d_max_similarity` is a
Tanimoto, not an energy. It is the first metric in this project whose direction
differs from the rest, and registering it is what turned
`rank_shortlist.LOWER_IS_BETTER` from a set of names into a map that carries a
direction -- see the note there.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem.Pharm2D import Generate, Gobbi_Pharm2D

from . import reference_set as rs

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger(__name__)

#: The rank metric this module produces. Registered in
#: `rank_shortlist.RANK_DIRECTION` as higher-is-better before anything sorts on
#: it; that registry is the guard that caught catalogue #4.
METRIC = "ph2d_max_similarity"
LOWER_IS_BETTER = False

#: The model-set rule (see the module docstring). Declared as constants so the
#: thresholds appear in the manifest rather than only in prose.
MAX_AMIDE_BONDS = 4          # strictly fewer than this
MAX_HEAVY_ATOMS = 45         # at most this

#: A backbone amide. Deliberately excludes the N-H-free tertiary case only when
#: it is not on a carbonyl; this is the standard peptide-bond pattern.
_AMIDE = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")


class PharmacophoreError(ValueError):
    """The scorer cannot be built or applied as configured."""


def fingerprint(smiles: str):
    """Gobbi 2D pharmacophore fingerprint, or None if the SMILES will not parse.

    Returns None rather than raising because a candidate pool of 30,000 CReM
    products legitimately contains a few unparseable rows, and the caller
    stamps those as unscored. A model-set molecule that fails to parse is a
    different matter and `build_model_set` refuses it loudly.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Generate.Gen2DFingerprint(mol, Gobbi_Pharm2D.factory)


def _peptide_like(mol) -> tuple[bool, int, int]:
    """(excluded?, n_amide, n_heavy) under the declared rule."""
    n_amide = len(mol.GetSubstructMatches(_AMIDE))
    n_heavy = mol.GetNumHeavyAtoms()
    return (n_amide >= MAX_AMIDE_BONDS or n_heavy > MAX_HEAVY_ATOMS,
            n_amide, n_heavy)


def build_model_set(hold_out: str | None) -> pd.DataFrame:
    """The reference binders that pass the model-set rule, minus ``hold_out``.

    ``hold_out`` is the canonical SMILES of the pool's own seed, and it is a
    REQUIRED argument with no default. A default of None would mean "score
    against everything" -- correct for an unseeded pool like T_1 and circular
    for every T_2 pool, with nothing at the call site to say which was meant.
    Pass None explicitly for an unseeded pool.

    The reference set is resolved by glob (`reference_set.latest_reference`),
    never pinned: `tests/test_reference_version.py` walks the AST for pins and
    five have gone stale in this project already.
    """
    ref = pd.read_csv(rs.latest_reference("pin1_reference_binders"))
    rows, dropped = [], []
    for rec in ref.to_dict("records"):
        name, smi = rec.get("name"), rec.get("canonical_smiles")
        mol = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
        if mol is None:
            dropped.append((name, "unparseable SMILES"))
            continue
        excluded, n_amide, n_heavy = _peptide_like(mol)
        if excluded:
            dropped.append((name, f"peptide-like ({n_amide} amides, "
                                  f"{n_heavy} heavy atoms)"))
            continue
        if hold_out is not None and smi == hold_out:
            dropped.append((name, "HELD OUT — this pool's own seed"))
            continue
        rows.append({"name": name, "canonical_smiles": smi,
                     "tier": rec.get("tier"), "n_amide": n_amide,
                     "n_heavy": n_heavy})

    for name, why in dropped:
        log.info("model set excludes %-36s %s", name, why)
    if not rows:
        raise PharmacophoreError(
            "the model set is empty after applying the rule; refusing to score "
            "against nothing")

    out = pd.DataFrame(rows)
    out["fp"] = [fingerprint(s) for s in out.canonical_smiles]
    # Every model molecule parsed above, so a None here would be a silent
    # fingerprinting failure rather than a bad input. Fail rather than score
    # against a shorter set than the log just reported.
    if out.fp.isna().any():
        bad = out.loc[out.fp.isna(), "name"].tolist()
        raise PharmacophoreError(f"fingerprinting failed for {bad}")
    log.info("model set: %d molecules (held out: %s)", len(out),
             "none" if hold_out is None else hold_out)
    return out


def score(smiles_iter, model: pd.DataFrame) -> pd.DataFrame:
    """Score candidates against a model set.

    Returns one row per input in input order, carrying the maximum and mean
    Tanimoto to the model set and the name of the nearest model molecule.

    THE NEAREST MODEL MOLECULE IS RETURNED ON PURPOSE. A max over a set is a
    number with no provenance -- it says how similar without saying to what.
    Carrying `ph2d_nearest_active` makes "this pool scores well because
    everything is nearest to the one binder that resembles its seed" visible in
    a value_counts rather than requiring someone to suspect it and recompute.
    """
    fps = list(model.fp)
    names = list(model.name)
    if not fps:
        raise PharmacophoreError("empty model set")

    max_sim, mean_sim, nearest = [], [], []
    for smi in smiles_iter:
        fp = fingerprint(smi) if isinstance(smi, str) else None
        if fp is None:
            max_sim.append(np.nan); mean_sim.append(np.nan); nearest.append(None)
            continue
        sims = np.array(DataStructs.BulkTanimotoSimilarity(fp, fps))
        j = int(sims.argmax())
        max_sim.append(float(sims[j]))
        mean_sim.append(float(sims.mean()))
        nearest.append(names[j])
    return pd.DataFrame({METRIC: max_sim,
                         "ph2d_mean_similarity": mean_sim,
                         "ph2d_nearest_active": nearest})
