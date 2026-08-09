"""
Purpose: prove the pose written to disk for a mode is the pose the selection chose.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-07
Input: none — synthetic pose clouds with known structure
Output: pass/fail

WHY THIS EXISTS. @tt8804, on the 2.2.0 viewer: *"medoids, which I am not
confident are being shown."* That doubt was well placed — nothing checked it.

The screen asserts that the SDF coordinates match the conformer that was
measured, which catches a corrupted write. It does **not** check that the
conformer chosen is the one the selection rule identified, and it does not check
that record *i* carries mode *i*. Those are the two ways the viewer could show a
confidently-labelled wrong pose, and neither would raise anything: the file would
be valid, the geometry self-consistent, and every number downstream unchanged.

The checks below are on synthetic clouds with a known answer, so a regression
fails here rather than being noticed months later in a figure.

  1. mode labels are assigned by POPULATION, largest first
  2. the representative of mode k is a MEMBER of mode k
  3. the representative is the medoid of the top anchoring quartile, which is
     NOT in general the best-anchored pose -- the distinction @tt8804 asked for
     when argmax was rejected as unrealistic
  4. `write_sdf` stamps record i with the mode it was given, so a reader
     selecting by the `mode` property gets the pose the selection chose
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import pose_modes as pmod                # noqa: E402


def _cloud(sizes, sep=25.0, spread=0.6, seed=0):
    """Well-separated spherical clusters of (position, direction) features."""
    rng = np.random.default_rng(seed)
    feats, truth = [], []
    for k, n in enumerate(sizes):
        centre = np.array([k * sep, 0.0, 0.0])
        pos = centre + rng.normal(0, spread, (n, 3))
        d = np.tile(np.array([0.0, 0.0, 1.0]), (n, 1))
        feats.append(np.hstack([pos, d]))
        truth += [k] * n
    return np.vstack(feats), np.array(truth)


def test_modes_are_numbered_by_population():
    # deliberately built smallest-first, so scan order and population disagree
    feat, truth = _cloud([40, 300, 120])
    lab = pmod.split(feat)
    sizes = {k: int((lab == k).sum()) for k in sorted(set(lab)) if k >= 0}
    assert sizes, "no modes found in a cleanly separated cloud"
    ordered = [sizes[k] for k in sorted(sizes)]
    assert ordered == sorted(ordered, reverse=True), (
        f"modes must be numbered by population, largest first; got {sizes}")
    # the 300-member cluster must be mode 0 regardless of where it sat
    assert sizes[0] == max(sizes.values())


def _select_rep(feat, labels, anchor, k):
    """The selection rule as `nac_screen_v2` applies it."""
    dmat = pmod.distances(feat)
    idx = np.flatnonzero(labels == k)
    a = anchor[idx]
    sub = dmat[np.ix_(idx, idx)]
    if np.all(np.isnan(a)) or len(idx) < 4:
        return int(idx[np.argmin(sub.mean(axis=1))])
    keep = np.flatnonzero(a >= np.nanpercentile(a, 75))
    if len(keep) < 2:
        return int(idx[np.argmin(sub.mean(axis=1))])
    s2 = sub[np.ix_(keep, keep)]
    return int(idx[keep[np.argmin(s2.mean(axis=1))]])


def test_representative_belongs_to_its_own_mode():
    feat, _ = _cloud([200, 90, 60], seed=3)
    lab = pmod.split(feat)
    rng = np.random.default_rng(1)
    anchor = rng.random(len(feat))
    for k in sorted({int(x) for x in lab if x >= 0}):
        rep = _select_rep(feat, lab, anchor, k)
        assert lab[rep] == k, (
            f"representative of mode {k} is labelled mode {lab[rep]} -- the "
            "viewer would draw a pose from a different cluster")


def test_representative_is_not_simply_the_best_anchored():
    """argmax was rejected as unrealistic; the rule must actually differ from it.

    A cloud with one well-anchored OUTLIER far from the cluster centre: argmax
    picks the outlier, the top-quartile medoid must not.
    """
    feat, _ = _cloud([120], seed=7)
    lab = np.zeros(len(feat), dtype=int)
    anchor = np.full(len(feat), 0.2)
    anchor[np.arange(30)] = 0.8              # a well-anchored quartile
    outlier = 0
    feat[outlier, :3] += 12.0                # ...one of which is far away
    anchor[outlier] = 1.0                    # and is the single best
    rep = _select_rep(feat, lab, anchor, 0)
    assert rep != outlier, (
        "the representative is the best-anchored pose -- this is argmax, the "
        "behaviour rejected for showing an unrepresentative best case")
    assert anchor[rep] >= 0.8 - 1e-9, (
        "the representative fell outside the top anchoring quartile")


def test_write_sdf_stamps_the_mode_it_was_given():
    """Record i must carry mode_ids[i], so a reader selects by identity."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "nac_screen_v2", REPO / "scripts" / "nac_screen_v2.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules["nac_screen_v2"] = m
    spec.loader.exec_module(m)

    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMultipleConfs(mol, numConfs=5, randomSeed=42)
    reps, mode_ids = [4, 2, 0], [0, 1, 2]        # deliberately not in order
    dest = Path("/tmp/_test_write_sdf.sdf")
    m.write_sdf(mol, reps, dest, modes=mode_ids)

    got = [int(x.GetProp("mode")) for x in Chem.SDMolSupplier(str(dest))
           if x is not None]
    assert got == mode_ids, f"records carry modes {got}, expected {mode_ids}"
    dest.unlink(missing_ok=True)


if __name__ == "__main__":
    fails = 0
    for name, fn in sorted(globals().items()):
        if not name.startswith("test_") or not callable(fn):
            continue
        try:
            fn()
            print(f"  PASS  {name}")
        except AssertionError as exc:
            fails += 1
            print(f"  FAIL  {name}\n        {exc}")
        except Exception as exc:                       # noqa: BLE001
            fails += 1
            print(f"  ERROR {name}\n        {type(exc).__name__}: {exc}")
    print(f"\n{'all checks passed' if not fails else f'{fails} FAILED'}")
    sys.exit(1 if fails else 0)
