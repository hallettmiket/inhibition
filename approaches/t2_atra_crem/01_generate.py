"""
Purpose: T_2 step 1 — degree-bounded CReM neighbourhood of ATRA.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: config/seeds.yaml (ATRA), config/approaches/t2_atra_crem.yaml
Output: the D2 frame — every generated derivative with its provenance

WHAT T_2 IS. Not de novo design and not a virtual library: the chemical
neighbourhood reachable from ATRA by single fragment swaps that ChEMBL33 has
actually seen. That makes its candidates synthetically plausible by
construction, at the cost of never leaving ATRA's chemotype.

DEGREE-BOUNDED, NOT EXHAUSTIVE (adversary control B3). `max_degree` is how many
successive rounds of mutation are applied. Degree 1 gives 10^4-10^5 molecules;
degree 2 gives 10^8-10^10, which no downstream stage could dock and which would
mostly be minor decorations of decorations. The cap is a scientific choice about
what neighbourhood means, not a performance workaround, so it is declared in
config and recorded in the manifest.

RADIUS IS A PROPERTY OF THE SEED (D0018). `radius` is how much surrounding
context must MATCH before CReM will swap a fragment — larger radius means FEWER
and more conservative replacements, which is the opposite of most people's first
guess. ATRA yields ZERO at radius 3 and 43 at radius 2, because a conjugated
polyene's 3-bond context has no precedent in ChEMBL33. Any new seed needs the
smoke test in `--probe` before its radius is pinned; inheriting another seed's
radius is how an approach silently generates nothing and reports success.

DEDUP ON InChIKey. Path multiplicity is high — the same molecule is reachable by
many edit sequences — so an undeduped frontier double-counts and, at degree > 1,
explodes combinatorially on the next round.

THE FRONTIER CAP IS ANNOUNCED. If generation hits `frontier_cap`, the run says
so in the log and in the manifest. A truncated enumeration reported as a
complete one is the quiet failure this stage is most likely to have.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
import yaml

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from shared import io as dio                      # noqa: E402
from shared import smiles as smi                  # noqa: E402

log = logging.getLogger("t2-generate")

EXPERIMENT = "02_t2_atra_crem"
APPROACH = "t2"
CONFIG = REPO / "config" / "approaches" / "t2_atra_crem.yaml"
SEEDS = REPO / "config" / "seeds.yaml"


def load_config() -> dict:
    cfg = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    exp = cfg["expansion"]
    return {
        "db": exp["fragment_db"],
        "radius": int(exp["radius"]),
        "max_degree": int(exp.get("max_degree", 1)),
        "frontier_cap": int(exp.get("frontier_cap", 200_000)),
        "mutate": exp.get("mutate") or {},
        "grow": exp.get("grow") or {},
        "seed_name": cfg["approach"]["seed"],
        "mechanism": cfg["approach"]["mechanism"],
    }


def load_seed(name: str) -> str:
    seeds = yaml.safe_load(SEEDS.read_text(encoding="utf-8"))["seeds"]
    s = seeds[name]["canonical_smiles"]
    canon = smi.canonical(s, strict=True)
    log.info("seed %s: %s", name, canon)
    return canon


def probe_radius(seed_smiles: str, db: str, radii=(1, 2, 3, 4, 5)) -> dict:
    """How many derivatives each radius yields for THIS seed.

    D0018's lesson as a runnable check: the usable radius belongs to the seed.
    A seed whose counts collapse to zero at the configured radius would
    enumerate an empty frontier and report success.
    """
    from crem.crem import grow_mol, mutate_mol
    from rdkit import Chem

    mol = Chem.MolFromSmiles(seed_smiles)
    out = {}
    for r in radii:
        try:
            m = len(set(mutate_mol(mol, db, radius=r, max_size=8, min_size=1)))
            g = len(set(grow_mol(mol, db, radius=r, max_atoms=8, min_atoms=1)))
        except Exception as exc:  # noqa: BLE001
            out[r] = {"error": str(exc)[:120]}
            continue
        out[r] = {"mutate": m, "grow": g}
        log.info("radius %d: %d mutations, %d grows", r, m, g)
    return out


def expand_once(smiles_list: list[str], cfg: dict) -> list[dict]:
    """One degree of expansion over a set of parents."""
    from crem.crem import grow_mol, mutate_mol
    from rdkit import Chem

    mut, grw = cfg["mutate"], cfg["grow"]
    rows = []
    for parent in smiles_list:
        mol = Chem.MolFromSmiles(parent)
        if mol is None:
            continue
        try:
            for child in mutate_mol(mol, cfg["db"], radius=cfg["radius"],
                                    min_size=int(mut.get("min_size", 1)),
                                    max_size=int(mut.get("max_size", 8)),
                                    max_inc=int(mut.get("max_inc", 4))):
                rows.append({"smiles": child, "parent_smiles": parent,
                             "operation": "mutate"})
            for child in grow_mol(mol, cfg["db"], radius=cfg["radius"],
                                  min_atoms=int(grw.get("min_atoms", 1)),
                                  max_atoms=int(grw.get("max_atoms", 8))):
                rows.append({"smiles": child, "parent_smiles": parent,
                             "operation": "grow"})
        except Exception as exc:  # noqa: BLE001
            log.warning("expansion failed for %s: %s", parent[:40], str(exc)[:120])
    return rows


def generate(cfg: dict, seed_smiles: str) -> tuple[pd.DataFrame, dict]:
    """Degree-bounded expansion with InChIKey dedup and an announced cap."""
    seen: set[str] = set()
    seed_key = smi.inchikey(seed_smiles)
    if seed_key:
        seen.add(seed_key)

    records: list[dict] = []
    frontier = [seed_smiles]
    truncated = False

    for degree in range(1, cfg["max_degree"] + 1):
        produced = expand_once(frontier, cfg)
        log.info("degree %d: %d raw products from %d parent(s)",
                 degree, len(produced), len(frontier))
        next_frontier: list[str] = []
        for row in produced:
            canon = smi.canonical(row["smiles"])
            if canon is None:
                continue
            key = smi.inchikey(canon)
            if not key or key in seen:
                continue
            seen.add(key)
            records.append({"canonical_smiles": canon,
                            "candidate_id": smi.candidate_id(canon, prefix=APPROACH),
                            "approach": APPROACH,
                            "parent_smiles": row["parent_smiles"],
                            "operation": row["operation"],
                            "degree": degree})
            next_frontier.append(canon)
            if len(records) >= cfg["frontier_cap"]:
                truncated = True
                break
        log.info("degree %d: %d new unique molecules (%d cumulative)",
                 degree, len(next_frontier), len(records))
        if truncated:
            log.warning("FRONTIER CAP %d reached at degree %d — enumeration is "
                        "TRUNCATED, not complete", cfg["frontier_cap"], degree)
            break
        if not next_frontier:
            log.warning("degree %d produced nothing new; stopping early", degree)
            break
        frontier = next_frontier

    df = pd.DataFrame(records)
    stats = {"n_unique": len(df), "truncated": truncated,
             "degrees_run": int(df["degree"].max()) if len(df) else 0}
    return df, stats


def main() -> None:
    ap = argparse.ArgumentParser(description="T_2: CReM neighbourhood of ATRA.")
    ap.add_argument("--probe", action="store_true",
                    help="report derivative counts per radius and exit (D0018)")
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

    cfg = load_config()
    seed = load_seed(cfg["seed_name"])

    if not Path(cfg["db"]).is_file():
        raise SystemExit(f"fragment DB not staged: {cfg['db']}\n"
                         "  python -m shared.sources stage")

    if args.probe:
        counts = probe_radius(seed, cfg["db"])
        print(f"\nCReM yield per radius for seed {cfg['seed_name']!r}:")
        for r, v in counts.items():
            print(f"  radius {r}: {v}")
        print(f"\nconfigured radius is {cfg['radius']}")
        if counts.get(cfg["radius"], {}).get("mutate", 0) == 0:
            print("  WARNING: the configured radius yields NOTHING for this seed")
        return

    df, stats = generate(cfg, seed)
    if df.empty:
        raise SystemExit(
            f"generation produced nothing at radius {cfg['radius']}. Run with "
            "--probe: the usable radius is a property of the seed (D0018).")

    if args.limit:
        df = df.head(args.limit)

    # Contract columns the schema requires of every approach's full frame.
    df["rejected_at"] = pd.NA
    df["mechanism"] = cfg["mechanism"]
    df["seed"] = cfg["seed_name"]

    out = dio.write_full_frame(
        df, approach=APPROACH, experiment=EXPERIMENT, stage="t2_generate",
        params={"engine": "crem", "fragment_db": cfg["db"],
                "radius": cfg["radius"], "max_degree": cfg["max_degree"],
                "mutate": cfg["mutate"], "grow": cfg["grow"],
                "frontier_cap": cfg["frontier_cap"],
                "frontier_truncated": stats["truncated"],
                "seed_smiles": seed},
        inputs={"fragment_db": Path(cfg["db"])})

    print(f"\nT_2 generation -> {out}")
    print(f"  {stats['n_unique']} unique derivatives, degree <= {stats['degrees_run']}")
    print(f"  radius {cfg['radius']} against {Path(cfg['db']).name}")
    if stats["truncated"]:
        print(f"  TRUNCATED at the {cfg['frontier_cap']} frontier cap — this is "
              "not the complete neighbourhood")
    print("\n  by operation:")
    for op, n in df["operation"].value_counts().items():
        print(f"    {op:8s} {n}")


if __name__ == "__main__":
    main()
