"""
Purpose: Test whether the internal residual's enrichment AUC reflects binding or
         molecular size, and what 2 actives can support either way.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-30
Input: the gate's MM-GBSA workdirs (ensemble_dg.json + ligand.mol2 per candidate)
Output: AUC per score, plus the AUC of ligand heavy-atom count used alone

THE QUESTION. D0037 measured an AUC of 0.831 for the internal residual against
0.181 for the standard interaction energy -- the artefact separating actives
from decoys better than the physics. Two explanations were left open: the
residual encodes something real, or it encodes a molecular property that happens
to correlate with being an active. Heavy-atom count is the obvious candidate,
because internal energy is extensive: more atoms, more bonded terms, larger
residual, with no reference to the pocket at all.

THE CONTROL is heavy-atom count scored on its own. If size alone reaches a
similar AUC, the residual is a size proxy and the 0.831 says nothing about
binding.
"""

from __future__ import annotations

import json
import pathlib
import itertools

import numpy as np

GATE = pathlib.Path(
    "/data/lab_vm/append_only/inhibition/00_shared_substrate/mmgbsa_gate")


def heavy_atoms(mol2: pathlib.Path) -> int | None:
    """Heavy-atom count from a mol2 ATOM block, hydrogens excluded."""
    if not mol2.is_file():
        return None
    n, inside = 0, False
    for line in mol2.read_text(errors="ignore").splitlines():
        s = line.strip()
        if s.startswith("@<TRIPOS>ATOM"):
            inside = True
            continue
        if s.startswith("@<TRIPOS>") and inside:
            break
        if inside and s:
            parts = s.split()
            if len(parts) >= 6 and not parts[5].split(".")[0].upper() == "H":
                n += 1
    return n or None


def auc(scores: np.ndarray, labels: np.ndarray) -> float:
    """Mann-Whitney AUC. Lower score = more active, matching dG convention."""
    pos, neg = scores[labels == 1], scores[labels == 0]
    if not len(pos) or not len(neg):
        return float("nan")
    wins = sum((p < n) + 0.5 * (p == n) for p in pos for n in neg)
    return float(wins) / (len(pos) * len(neg))


def exact_ci(scores: np.ndarray, labels: np.ndarray) -> tuple[float, float]:
    """The full range of AUCs obtainable by relabelling which molecules are the
    actives. With 2 actives this is the honest interval: it says what the
    statistic could have been had a different pair of the same molecules been
    the known binders."""
    n_pos = int((labels == 1).sum())
    idx = range(len(scores))
    vals = []
    for combo in itertools.combinations(idx, n_pos):
        lab = np.zeros(len(scores), dtype=int)
        lab[list(combo)] = 1
        vals.append(auc(scores, lab))
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


rows = []
for wd in sorted(GATE.iterdir()):
    if not wd.is_dir():
        continue
    f = wd / "ensemble" / "ensemble_dg.json"
    if not f.is_file():
        continue
    d = json.loads(f.read_text())
    ha = heavy_atoms(wd / "ligand.mol2")
    if ha is None:
        continue
    rows.append({
        "id": wd.name,
        "label": 1 if wd.name.startswith("act") else 0,
        "interaction": d.get("dG_interaction_kcal"),
        "residual": d.get("dG_internal_residual_kcal"),
        "full": d.get("dG_kcal_ensemble"),
        "heavy_atoms": ha,
    })

lab = np.array([r["label"] for r in rows])
print(f"candidates: {len(rows)}   actives: {int(lab.sum())}   "
      f"decoys: {int((lab == 0).sum())}\n")

print(f"{'score':<28} {'AUC':>7}   {'95% relabelling range':>24}")
print("-" * 64)
for key, name, sign in [
        ("interaction", "interaction energy", 1.0),
        ("residual", "internal residual", 1.0),
        ("full", "full potential (as shown)", 1.0),
        ("heavy_atoms", "heavy-atom count ALONE", -1.0),
]:
    s = np.array([r[key] if r[key] is not None else np.nan for r in rows],
                 dtype=float)
    ok = ~np.isnan(s)
    a = auc(sign * s[ok], lab[ok])
    lo, hi = exact_ci(sign * s[ok], lab[ok])
    print(f"{name:<28} {a:>7.3f}   [{lo:>7.3f}, {hi:>7.3f}]")

# Does the residual track size directly?
s = np.array([r["residual"] for r in rows], dtype=float)
h = np.array([r["heavy_atoms"] for r in rows], dtype=float)
ok = ~np.isnan(s)
r_pearson = float(np.corrcoef(s[ok], h[ok])[0, 1])
order_s, order_h = np.argsort(np.argsort(s[ok])), np.argsort(np.argsort(h[ok]))
r_spearman = float(np.corrcoef(order_s, order_h)[0, 1])
print(f"\ncorr(residual, heavy atoms): Pearson {r_pearson:+.3f}   "
      f"Spearman {r_spearman:+.3f}")

act = [r for r in rows if r["label"] == 1]
print(f"\nthe two actives: " + ", ".join(
    f"{r['id']} ({r['heavy_atoms']} heavy atoms)" for r in act))
print(f"decoy heavy atoms: median {np.median(h[lab == 0]):.0f}, "
      f"range {h[lab == 0].min():.0f}-{h[lab == 0].max():.0f}")
