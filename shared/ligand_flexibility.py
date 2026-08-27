"""Per-atom flexibility from a cheap conformer ensemble, and what it may be used for.

MOVED OUT OF `exp/15_rmsf_predictor` BECAUSE PRODUCTION NEEDS IT. It lived in an
experiment directory while `shared/pose_contacts` depended on it, which meant the
splitter could only run from code that had put an `exp/` path on `sys.path` -- and
two modules named `run_all` then had to be imported by file to avoid shadowing
each other.

WHAT IT IS FOR, AND WHAT IT IS NOT FOR. This is the distinction D0094 exists to
record, and it was got wrong once:

  VALIDATED -- ranking ATOMS WITHIN one molecule. rho = +0.657 against measured
  MD RMSF over 147 swept modes, 100% positive. That is what `atom_weights` needs:
  a wagging tail should count for less than a rigid core.

  NOT VALIDATED -- ranking MOLECULES against each other on absolute scale. Across
  molecules the same prediction correlates with measurement at rho = +0.112, CI
  [-0.06, +0.27], crossing zero. The prediction varies at CV 0.15 where the truth
  varies at 0.45, and `median(rmsf)/2.21` does not beat writing one number down
  for every molecule (Wilcoxon p = 0.515).

So `predict_rmsf` is safe to use for weights and unsafe to use for a per-molecule
LENGTH SCALE. `pose_contacts.tolerance_for` is the one caller that needs the
second thing, and it says so.

BOND ORDERS COME FROM THE TEMPLATE, NOT FROM COORDINATES. Perceiving bonds from a
pose's coordinates produced [S+2], [C+] and [N-] where the real molecule has none
-- and still returned a correlation of +0.51, which is the trap.

NO HEAVY-ONLY MOLECULE IS EVER BUILT. Deleting hydrogens to get a heavy-atom
molecule makes an aromatic N-H unkekulizable; only an INDEX MAPPING is carried,
and conformers are embedded on the intact molecule and read at the mapped
positions.
"""

from __future__ import annotations

import numpy as np

#: Conformers in the ensemble. Converged: 50 sits within 1.4-2.1% of 100 and 200,
#: and seed-to-seed movement is 3.5% (exp/18). More does not help.
DEFAULT_CONFORMERS = 50
DEFAULT_SEED = 7


def predict_rmsf(template, heavy_idx: list[int], n_conf: int = DEFAULT_CONFORMERS,
                 seed: int = DEFAULT_SEED) -> np.ndarray:
    """Per-atom spread over an aligned conformer ensemble, in template atom order.

    `template` supplies the chemistry (an RDKit mol read from a pose SDF);
    `heavy_idx` selects which atoms come back, in the order given.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolAlign
    mh = Chem.Mol(template)
    cids = AllChem.EmbedMultipleConfs(mh, numConfs=n_conf, randomSeed=seed)
    if len(cids) == 0:
        cids = AllChem.EmbedMultipleConfs(mh, numConfs=n_conf, randomSeed=seed,
                                          useRandomCoords=True)
    if len(cids) < 5:
        raise ValueError(f"only {len(cids)} conformers embedded")
    AllChem.MMFFOptimizeMoleculeConfs(mh, maxIters=300)
    amap = [(i, i) for i in heavy_idx]
    for c in list(cids)[1:]:
        rdMolAlign.AlignMol(mh, mh, prbCid=c, refCid=cids[0], atomMap=amap)
    X = np.array([mh.GetConformer(c).GetPositions()[heavy_idx] for c in cids])
    return np.sqrt(((X - X.mean(0)) ** 2).sum(-1).mean(0))
