"""
Purpose: T_2 step 1 — degree-bounded CReM neighbourhood of ATRA.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: config/seeds.yaml (ATRA), config/approaches/t2_atra_crem.yaml
Output: the D2 frame — every generated derivative with its provenance

WHAT T_2 IS. Not de novo design and not a virtual library: the chemical
neighbourhood reachable from a KNOWN BINDER by single fragment swaps that
ChEMBL33 has actually seen. That makes its candidates synthetically plausible by
construction, at the cost of never leaving the seed's chemotype.

ONE APPROACH, MANY SEEDS (2026-07-31). `--seed` selects any T_2-admissible seed
from config/seeds.yaml; each writes to its OWN experiment directory while
keeping approach id `t2`, so the schema, gates and decision vocabulary are
untouched and no two seeds' frames can interleave. Running the same operator
from five starting points is what turns "never leaves the seed's chemotype"
from a limitation into a controlled comparison.

BUT NEIGHBOURHOOD SIZE IS NOT COMPARABLE ACROSS SEEDS. Measured at radius 2:
ATRA 1,882, Guo-Pfizer 8,674, Potter-Astex 9,192, Du-Xu 9,744, Liu-2024-C3
16,817 — a 9x spread. More draws yields a better best-of-pool for extreme-value
reasons alone, so comparing top scores between seeds is an artefact unless N is
matched. Compare distributions, or best-of-N at matched N.

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
from shared import pocket_size as ps              # noqa: E402
from shared import seeds as sd                    # noqa: E402
from shared import smiles as smi                  # noqa: E402

log = logging.getLogger("t2-generate")

APPROACH = "t2"
CONFIG = REPO / "config" / "approaches" / "t2_atra_crem.yaml"
SEEDS = REPO / "config" / "seeds.yaml"


def load_config(seed_name: str | None = None, *,
                require_radius: bool = True) -> dict:
    cfg = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    exp = cfg["expansion"]
    seed_name = seed_name or cfg["approach"]["seed"]
    try:
        rec = sd.resolve(APPROACH, seed_name, require_radius=require_radius)
    except sd.SeedError as exc:
        raise SystemExit(str(exc)) from exc
    return {
        "db": exp["fragment_db"],
        # RADIUS COMES FROM THE SEED, not from the approach config. The
        # approach-level `radius` remains as ATRA's historical value and is
        # deliberately not read here.
        "radius": int(rec["radius"]) if rec.get("radius") is not None else None,
        "max_degree": int(exp.get("max_degree", 1)),
        "frontier_cap": int(exp.get("frontier_cap", 200_000)),
        "max_heavy_atoms": int(exp.get("max_heavy_atoms",
                                       ps.MAX_HEAVY_ATOMS)),
        "mutate": exp.get("mutate") or {},
        "grow": exp.get("grow") or {},
        "seed_name": seed_name,
        "experiment": rec["t2_experiment"],
        "mechanism": cfg["approach"]["mechanism"],
    }


def load_seed(name: str) -> str:
    seeds = yaml.safe_load(SEEDS.read_text(encoding="utf-8"))["seeds"]
    s = seeds[name]["canonical_smiles"]
    canon = smi.canonical(s, strict=True)
    log.info("seed %s: %s", name, canon)
    return canon


def probe_radius(seed_smiles: str, db: str, cfg: dict,
                 radii=(1, 2, 3, 4, 5)) -> dict:
    """How many derivatives each radius yields for THIS seed.

    D0018's lesson as a runnable check: the usable radius belongs to the seed.
    A seed whose counts collapse to zero at the configured radius would
    enumerate an empty frontier and report success.

    THE PROBE MUST RUN THE SAME ENUMERATION AS `expand_once`, OR IT IS
    MEASURING A DIFFERENT QUESTION. It previously hardcoded `max_size=8,
    min_size=1` / `max_atoms=8, min_atoms=1` and omitted `max_inc` entirely,
    while `expand_once` reads all four from config and passes `max_inc`.
    Measured on ATRA at radius 2 (2026-08-04): the probe reported 113 mutate +
    1,679 grow = **1,792**, and the run produces 203 + 1,679 = **1,882**. The
    check that decides whether a radius is usable was reporting a stricter
    enumeration than the one about to happen, and any config change to the size
    parameters would have widened the gap silently.
    """
    from crem.crem import grow_mol, mutate_mol
    from rdkit import Chem

    mut, grw = cfg["mutate"], cfg["grow"]
    mol = Chem.MolFromSmiles(seed_smiles)
    out = {}
    for r in radii:
        try:
            m = len(set(mutate_mol(mol, db, radius=r,
                                   min_size=int(mut.get("min_size", 1)),
                                   max_size=int(mut.get("max_size", 8)),
                                   max_inc=int(mut.get("max_inc", 4)))))
            g = len(set(grow_mol(mol, db, radius=r,
                                 min_atoms=int(grw.get("min_atoms", 1)),
                                 max_atoms=int(grw.get("max_atoms", 8)))))
        except Exception as exc:  # noqa: BLE001
            # AN ERROR IS NOT A MEASURED ZERO. Callers used to read this with
            # `.get("mutate", 0)`, so a DB failure at radius 3 became "radius 3
            # yields nothing" -- the same shape as the mmCIF parser that
            # reported zero covalent entries across all 190 and read as a
            # finding. `total` is left absent so it cannot be summed.
            out[r] = {"error": str(exc)[:120]}
            log.warning("radius %d: PROBE FAILED (%s) — this is not a zero",
                        r, str(exc)[:80])
            continue
        out[r] = {"mutate": m, "grow": g, "total": m + g}
        log.info("radius %d: %d mutations, %d grows, %d total", r, m, g, m + g)
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
        pruned_oversize = 0
        for row in produced:
            canon = smi.canonical(row["smiles"])
            if canon is None:
                continue
            key = smi.inchikey(canon)
            if not key or key in seen:
                continue
            seen.add(key)
            # THE GOVERNOR (issue #1, note 2). Pruned BETWEEN enumeration and
            # reseeding, not afterwards: an oversized molecule that becomes a
            # seed produces an entire oversized lineage, so keeping it costs
            # more at every subsequent degree. Its InChIKey stays in `seen`, so
            # a pruned molecule is not re-enumerated by another parent later.
            #
            # Deliberately weak. 55 heavy atoms is ~1.6x what the 1018 A^3
            # pocket admits at tight packing and clears every known
            # non-peptidic Pin1 binder (shared/pocket_size.py). It removes what
            # cannot fit under any pose; it does not shape the chemistry.
            if not ps.fits_pocket(canon, max_heavy=cfg["max_heavy_atoms"]):
                pruned_oversize += 1
                continue
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
        log.info("degree %d: %d new unique molecules (%d cumulative); "
                 "governor pruned %d oversize (> %d heavy atoms)",
                 degree, len(next_frontier), len(records), pruned_oversize,
                 cfg["max_heavy_atoms"])
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
    sd.add_seed_argument(ap, APPROACH)
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

    cfg = load_config(args.seed, require_radius=not args.probe)
    seed = load_seed(cfg["seed_name"])
    log.info("seed %s -> experiment %s (radius %s)",
             cfg["seed_name"], cfg["experiment"], cfg["radius"])

    if not Path(cfg["db"]).is_file():
        raise SystemExit(f"fragment DB not staged: {cfg['db']}\n"
                         "  python -m shared.sources stage")

    if args.probe:
        counts = probe_radius(seed, cfg["db"], cfg)
        print(f"\nCReM yield per radius for seed {cfg['seed_name']!r}:")
        for r, v in counts.items():
            print(f"  radius {r}: {v}")
        if cfg["radius"] is None:
            best = max((r for r, v in counts.items()
                        if v.get("mutate", 0) + v.get("grow", 0) > 0),
                       default=None)
            print(f"\nseed {cfg['seed_name']!r} has NO pinned radius yet.")
            print(f"  largest radius that still yields anything: {best}")
            print(f"  pin `radius:` under seeds.yaml -> {cfg['seed_name']} "
                  "before running generation.")
        else:
            print(f"\nconfigured radius is {cfg['radius']}")
            # BOTH OPERATORS, NOT JUST mutate. This checked `mutate` alone
            # while `grow` supplies 1,679 of ATRA's 1,882 products -- 89% of
            # the frontier. A seed whose `grow` collapsed to zero passed the
            # check silently, and the `best` calculation eight lines above
            # already used mutate+grow, so two adjacent code paths disagreed
            # about what "yields anything" means.
            at_radius = counts.get(cfg["radius"], {})
            if "error" in at_radius:
                # Distinguished from a zero on purpose -- see probe_radius.
                print(f"  WARNING: the probe FAILED at the configured radius "
                      f"({at_radius['error']}). This is NOT a measurement that "
                      "the radius yields nothing; it is a missing measurement.")
            elif at_radius.get("mutate", 0) + at_radius.get("grow", 0) == 0:
                print("  WARNING: the configured radius yields NOTHING "
                      "for this seed (neither mutate nor grow)")
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
        df, approach=APPROACH, experiment=cfg["experiment"],
        stage="t2_generate",
        params={"engine": "crem", "fragment_db": cfg["db"],
                "radius": cfg["radius"], "max_degree": cfg["max_degree"],
                "mutate": cfg["mutate"], "grow": cfg["grow"],
                "frontier_cap": cfg["frontier_cap"],
                "frontier_truncated": stats["truncated"],
                "seed_name": cfg["seed_name"],
                "seed_smiles": seed},
        inputs={"fragment_db": Path(cfg["db"])})

    print(f"\nT_2 generation (seed {cfg['seed_name']}) -> {out}")
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
