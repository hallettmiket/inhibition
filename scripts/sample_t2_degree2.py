"""
Purpose: a UNIFORM RANDOM SAMPLE of T_2's degree-2 neighbourhood, not a truncation.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-31
Input: a seed's existing degree-1 D2 frame
Output: a sampled degree-2 D2 frame in its own experiment directory

Run:  /data/lab_vm/envs/dwi_cheminf/bin/python3 scripts/sample_t2_degree2.py \
        --seed atra --target 30000 [--workers 24] [--rng-seed 20260731]

WHY SAMPLE RATHER THAN ENUMERATE.

Measured 2026-07-31 on 30 random degree-1 parents of ATRA: mean 4,145 new
unique children each, union 99,885 over 30 parents with only 1.24x overlap.
Extrapolated over all 1,882 parents that is ~6.3M unique molecules -- NOT the
10^8-10^10 the approach's docstring estimated, which was wrong by two to four
orders of magnitude.

Enumeration is therefore cheap: ~7 CPU-hours single-core. DOCKING is not. At the
measured 2.53 s/molecule on one card (T_2's own D2_4 manifest: 1,882 molecules
in 4,763 s) the full degree-2 set is ~184 GPU-days on one GPU, ~31 days on six.
That is what makes exhaustive degree 2 infeasible, and nothing about the
enumeration.

WHY NOT JUST RAISE frontier_cap. Because `frontier_cap` TRUNCATES IN PARENT
ORDER. At 200,000 it stops around parent 48 of 1,882 -- the complete degree-2
expansion of the first 2.5% of the frontier, which is not a degree-2
neighbourhood and is not a sample of one. It is announced (the log says so, and
the manifest carries `frontier_truncated`), so it is not silent; it is simply
not interpretable. A uniform sample is.

HOW THE SAMPLE IS DRAWN. Every parent is expanded in full; each NEW unique child
that also clears the pocket governor enters a RESERVOIR of exactly `target`
items (Algorithm R). Every eligible molecule ends up equally likely to be in the
reservoir, regardless of which parent produced it or where that parent sat in
the frame — and the sample is exactly the size asked for.

THIS REPLACED BERNOULLI-AGAINST-AN-ESTIMATE, WHICH WAS WRONG BY 2x. The old
draw kept each child with p = target / (n_parents * MEAN_CHILDREN_PER_PARENT),
extrapolating linearly from a 30-parent probe. The deduplicated union grows
SUBLINEARLY, because parents' children overlap more and more as the pool grows,
so the denominator was inflated and p was too small. Measured on the ATRA run:
estimated 7,800,890 against a realised 4,063,427 — 1.92x high, so the frame kept
15,653 where the target was 30,000. The sample was still UNBIASED; only its SIZE
was wrong, and nothing detected it. Worse, MEAN_CHILDREN_PER_PARENT was measured
on ATRA and would have been applied unchanged to seeds whose degree-1
neighbourhoods differ by up to 9x. The reservoir needs no population estimate at
all, so the whole class of error is gone.

DEGREE 2 IS TERMINAL HERE. Nothing sampled is reseeded -- this measures whether
a second edit from the seed finds anything better, and a third would need its
own justification.

THE OUTPUT IS LABELLED A SAMPLE. `sampling_fraction`, `rng_seed` and
`estimated_population` go in the manifest and `is_sample` on every row, because
a sampled frame that later gets read as a complete enumeration would make every
count derived from it wrong by a factor of ~200.
"""

from __future__ import annotations

import argparse
import logging
import os
import random
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import io as dio                      # noqa: E402
from shared import pocket_size as ps              # noqa: E402
from shared import seeds as sd                    # noqa: E402
from shared import smiles as smi                  # noqa: E402

log = logging.getLogger("t2-degree2-sample")

APPROACH = "t2"

# Measured, not assumed. Mean new-unique children per degree-1 parent, from the
# 30-parent probe described above. Used only to set the Bernoulli probability;
# the realised population is reported afterwards so the estimate can be checked.
MEAN_CHILDREN_PER_PARENT = 4145.0


def _expand_one(args: tuple) -> list[str]:
    """Expand ONE parent and return (canonical_smiles, inchikey, fits) per child.

    THE PER-MOLECULE RDKit WORK HAPPENS HERE, NOT IN THE PARENT. Global dedup
    still happens in the parent -- a per-worker set would dedup within one
    parent only, and the whole point is that the same molecule is reachable
    from many parents -- but the parent now only needs a SET LOOKUP on a
    precomputed key.

    MEASURED, WHICH IS WHY THIS CHANGED. Returning raw SMILES left the parent
    doing three RDKit calls per child (canonical, inchikey, fits_pocket) on one
    thread while 48 workers waited: the parent sat at 99.7% CPU and the workers
    at 59%. On potter_astex that was ~1,940 s per 100 parents, dead linear --
    compute-bound on those calls, not on the growing dedup set -- projecting
    ~40 h for one seed and ~330 h for the four remaining. The expensive part is
    embarrassingly parallel and was running serially.

    `fits` is computed here too: the pocket governor is another RDKit parse,
    and it is applied AFTER dedup in the parent, so the flag travels with the
    molecule rather than being recomputed.
    """
    parent, cfg = args
    os.nice(19)
    sys.path.insert(0, str(REPO))
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    from crem.crem import grow_mol, mutate_mol

    from shared import pocket_size as _ps
    from shared import smiles as _smi

    mol = Chem.MolFromSmiles(parent)
    if mol is None:
        return []
    mut, grw = cfg["mutate"], cfg["grow"]
    raw: list[str] = []
    try:
        for child in mutate_mol(mol, cfg["db"], radius=cfg["radius"],
                                min_size=int(mut.get("min_size", 1)),
                                max_size=int(mut.get("max_size", 8)),
                                max_inc=int(mut.get("max_inc", 4))):
            raw.append(child)
        for child in grow_mol(mol, cfg["db"], radius=cfg["radius"],
                              min_atoms=int(grw.get("min_atoms", 1)),
                              max_atoms=int(grw.get("max_atoms", 8))):
            raw.append(child)
    except Exception as exc:  # noqa: BLE001 - one bad parent must not end the run
        log.warning("expansion failed for %s: %s", parent[:40], str(exc)[:120])

    out: list[tuple[str, str, bool]] = []
    max_heavy = cfg["max_heavy_atoms"]
    for s in raw:
        canon = _smi.canonical(s)
        if canon is None:
            continue
        key = _smi.inchikey(canon)
        if not key:
            continue
        out.append((canon, key, _ps.fits_pocket(canon, max_heavy=max_heavy)))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seed", default="atra",
                    help="which T_2 seed's degree-1 frame to expand")
    ap.add_argument("--target", type=int, default=30000,
                    help="approximate number of degree-2 molecules to keep")
    ap.add_argument("--workers", type=int, default=24,
                    help="parent-expansion processes (CPU budget is 50 cores)")
    ap.add_argument("--rng-seed", type=int, default=20260731,
                    help="pinned so the sample is reproducible")
    ap.add_argument("--experiment", default=None,
                    help="output experiment dir (default: <seed dir>_degree2)")
    ap.add_argument("--limit-parents", type=int, default=None,
                    help="smoke testing only: expand just N parents")
    ap.add_argument("--fragment-db-local", default=None,
                    help="read the CReM fragment DB from this LOCAL copy "
                         "instead of the governed path (e.g. /dev/shm/...). "
                         "The manifest still records the canonical path, and "
                         "the copy is SHA-256 verified against it first.")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s %(name)s: %(message)s")

    # Reuse the approach's own generation config so the operator is identical to
    # degree 1 -- same radius, same fragment DB, same size pins. A degree-2
    # sample drawn with different pins would not be comparable to the degree-1
    # frame it is being compared against, which is the entire experiment.
    sys.path.insert(0, str(REPO / "approaches" / "t2_atra_crem"))
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "t2gen", REPO / "approaches" / "t2_atra_crem" / "01_generate.py")
    t2gen = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(t2gen)

    cfg = t2gen.load_config(args.seed)
    rec = sd.resolve(APPROACH, args.seed)
    src_experiment = rec["experiment"]
    out_experiment = args.experiment or f"{src_experiment}_degree2"

    # THE FRAGMENT DB IS A 2 GB SQLite FILE ON NFS, AND 160 WORKERS DOING SMALL
    # RANDOM READS AGAINST IT OVER NFS IS A PATHOLOGICAL PATTERN. Measured on
    # potter_astex at 160 workers: 79 of 210 processes sat in `D` state
    # (uninterruptible I/O wait) and only 137 of the 160 requested cores were
    # busy, so the campaign was I/O-bound, not CPU-bound. Reading the same bytes
    # from a local copy removes the wait entirely.
    #
    # PROVENANCE IS NOT ALLOWED TO MOVE WITH THE BYTES. The manifest keeps
    # recording the canonical governed path, because that is the input the run
    # consumed; a manifest naming `/dev/shm/...` would record a location that
    # does not exist by the next reboot and cannot be checked by anyone. The
    # copy is SHA-256 verified against the original BEFORE use, so "same bytes"
    # is established rather than assumed -- a stale or truncated scratch copy
    # would otherwise produce a perfectly plausible, quietly different
    # neighbourhood.
    canonical_db = cfg["db"]
    if args.fragment_db_local:
        local = Path(args.fragment_db_local)
        if not local.is_file():
            raise SystemExit(f"--fragment-db-local not found: {local}")
        import hashlib

        def _sha(p: Path) -> str:
            h = hashlib.sha256()
            with open(p, "rb") as fh:
                for chunk in iter(lambda: fh.read(1 << 22), b""):
                    h.update(chunk)
            return h.hexdigest()

        log.info("verifying local fragment DB against the governed copy…")
        a, b = _sha(Path(canonical_db)), _sha(local)
        if a != b:
            raise SystemExit(
                f"local fragment DB does not match the governed one\n"
                f"  {canonical_db}  {a[:16]}\n  {local}  {b[:16]}\n"
                "Refusing to enumerate a different neighbourhood than the "
                "manifest would claim.")
        log.info("local DB verified (sha256 %s); reading from %s", a[:16], local)
        cfg["db"] = str(local)

    df1, frame1 = dio.latest_frame(src_experiment, APPROACH)
    d1 = df1[df1["degree"] == 1].drop_duplicates("canonical_smiles")
    parents = list(d1["canonical_smiles"])
    if args.limit_parents:
        parents = parents[:args.limit_parents]

    # RESERVOIR, NOT BERNOULLI-AGAINST-AN-ESTIMATE. The population no longer has
    # to be guessed at all, which removes the defect that halved the ATRA run.
    #
    # `p = target / (n_parents * MEAN_CHILDREN_PER_PARENT)` extrapolated
    # LINEARLY from a 30-parent probe and ignored that the deduplicated union
    # grows SUBLINEARLY as parents' children overlap. Measured on the ATRA run
    # that used it: estimated 7,800,890, realised 4,063,427 -- 1.92x high, so p
    # was half what it should have been and the frame kept 15,653 against a
    # target of 30,000. The sample was still UNBIASED (every molecule had the
    # same p); only its SIZE was wrong, and nothing flagged it. The docstring
    # claimed the realised size "varies slightly around target".
    #
    # Worse for the other seeds: 4145.0 was measured on ATRA and would be
    # applied by name to seeds whose degree-1 neighbourhoods differ by up to
    # 9x. A constant measured in one context, used in another.
    #
    # Reservoir sampling returns EXACTLY `target` items, uniformly, in one pass,
    # in O(target) memory, with no estimate of any kind. If the population is
    # smaller than the target it returns all of it -- the honest answer rather
    # than a scaled-down one.
    log.info("seed %s: %d degree-1 parents from %s", args.seed, len(parents),
             frame1.name)
    log.info("RESERVOIR sampling to exactly %d (the old estimator would have "
             "guessed a population of %.0f; it is no longer used)",
             args.target, len(parents) * MEAN_CHILDREN_PER_PARENT)

    # Seed the dedup set with degree 1 AND the seed itself, so "new at degree 2"
    # means genuinely new rather than rediscovered one edit later.
    seen: set[str] = set()
    seed_smiles = t2gen.load_seed(cfg["seed_name"])
    if (sk := smi.inchikey(seed_smiles)):
        seen.add(sk)
    for s in d1["canonical_smiles"]:
        k = smi.inchikey(s)
        if k:
            seen.add(k)

    rng = random.Random(args.rng_seed)
    kept: list[dict] = []
    n_raw = n_unique = n_oversize = n_eligible = 0
    t0 = time.time()

    payload = [(p, cfg) for p in parents]
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(_expand_one, item): item[0] for item in payload}
        for i, fut in enumerate(as_completed(futs), 1):
            parent = futs[fut]
            children = fut.result()
            n_raw += len(children)
            # The worker already canonicalised, keyed and governor-checked each
            # child; the parent's only serial job is the GLOBAL dedup set, which
            # is a hash lookup. See `_expand_one` for the measurement that moved
            # the RDKit work out of this loop.
            for canon, key, fits in children:
                if key in seen:
                    continue
                seen.add(key)
                n_unique += 1
                # The governor still applies. It prunes nothing at degree 1 for
                # any declared seed, but a degree-2 product can carry two
                # growths and is exactly what the ceiling exists for.
                if not fits:
                    n_oversize += 1
                    continue
                # RESERVOIR, AFTER DEDUP AND AFTER THE GOVERNOR. Sampling
                # before dedup would over-represent molecules reachable by many
                # edit paths; sampling before the governor would spend
                # reservoir slots on molecules that cannot fit the pocket.
                #
                # `n_eligible` counts everything that reached this point, which
                # is the population the sample is drawn FROM -- so the realised
                # fraction is a measured quantity afterwards rather than an
                # assumption beforehand.
                n_eligible += 1
                row = {
                    "canonical_smiles": canon,
                    "candidate_id": smi.candidate_id(canon, prefix=APPROACH),
                    "approach": APPROACH,
                    "parent_smiles": parent,
                    "degree": 2,
                }
                if len(kept) < args.target:
                    kept.append(row)
                else:
                    # Algorithm R: the j-th eligible item replaces a uniformly
                    # chosen slot with probability target/j, which leaves every
                    # item seen so far equally likely to be in the reservoir.
                    j = rng.randrange(n_eligible)
                    if j < args.target:
                        kept[j] = row
            if i % 100 == 0:
                log.info("[%d/%d parents] %d unique, %d kept, %.0fs",
                         i, len(parents), n_unique, len(kept), time.time() - t0)

    if not kept:
        raise SystemExit("sample is empty — check --target and the parent frame")

    out_df = pd.DataFrame(kept)
    out_df["rejected_at"] = pd.NA
    out_df["mechanism"] = cfg["mechanism"]
    out_df["seed"] = cfg["seed_name"]
    # EVERY ROW CARRIES THE FACT THAT IT IS SAMPLED. A frame that loses this in
    # the manifest alone is one join away from being counted as a census.
    out_df["is_sample"] = True
    # MEASURED, not assumed: the fraction is now a RESULT.
    out_df["sampling_fraction"] = len(kept) / n_eligible if n_eligible else 0.0

    realised_population = n_eligible
    out = dio.write_full_frame(
        out_df, approach=APPROACH, experiment=out_experiment,
        stage="t2_generate_degree2_sample",
        # THE CANONICAL PATH, never the scratch copy it was read from. The
        # bytes are identical (SHA-256 verified above); the governed path is
        # the one a reader can check.
        params={"engine": "crem", "fragment_db": canonical_db,
                "fragment_db_read_from":
                    args.fragment_db_local or canonical_db,
                "radius": cfg["radius"], "degree": 2,
                "mutate": cfg["mutate"], "grow": cfg["grow"],
                "seed_name": cfg["seed_name"],
                "source_experiment": src_experiment,
                "n_parents": len(parents),
                "sampling": "reservoir (Algorithm R) after global inchikey "
                            "dedup and the pocket governor",
                "sampling_fraction": (len(kept) / n_eligible) if n_eligible else 0.0,
                "rng_seed": args.rng_seed,
                "target": args.target,
                # No estimate is used any more. Recorded only so a reader can
                # see how far the retired estimator would have been off.
                "legacy_estimator_would_have_said":
                    len(parents) * MEAN_CHILDREN_PER_PARENT,
                "realised_population": realised_population,
                "n_kept": len(out_df),
                "governor_pruned_oversize": n_oversize,
                "IS_A_SAMPLE": True},
        inputs={"degree1_frame": frame1})

    print(f"\nT_2 degree-2 SAMPLE (seed {cfg['seed_name']}) -> {out}")
    print(f"  parents expanded      {len(parents)}")
    print(f"  raw products          {n_raw}")
    print(f"  new unique (degree 2) {n_unique}")
    print(f"  governor pruned       {n_oversize}")
    print(f"  population (post-gov) {realised_population}")
    legacy = len(parents) * MEAN_CHILDREN_PER_PARENT
    frac = (len(out_df) / realised_population) if realised_population else 0.0
    print(f"  kept                  {len(out_df)} of {realised_population} "
          f"(fraction {frac:.5f}, MEASURED not assumed)")
    print(f"  retired estimator     would have said {legacy:.0f} "
          f"({legacy / realised_population:.2f}x the truth)"
          if realised_population else "")
    print(f"  elapsed               {time.time() - t0:.0f}s")
    print("\n  THIS IS A SAMPLE, NOT AN ENUMERATION. Any count derived from it "
          f"must be scaled by 1/{frac:.5f}.")


if __name__ == "__main__":
    main()
