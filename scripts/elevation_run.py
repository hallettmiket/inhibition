"""
Purpose: run the pre-registered elevation experiment — which ranking metric selects for physical stability?
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06
Input: 00_outputs/blacksmith/elevation_cohort/elevation_cohort_<N>.csv + the crystallographic positives
Output: 00_outputs/blacksmith/elevation_tier1/elevation_t1_s<shard>_<N>.csv
        00_outputs/blacksmith/elevation_tier2/  (tier 2 rows land in the bpmd topic)

THE DESIGN IS FIXED IN `docs/elevation_prereg.md` AND THIS SCRIPT DOES NOT GET A
VOTE ON IT. Groups, readouts and readings were committed before any simulation
ran. Anything this script computes that the prereg did not name is labelled
post-hoc where it is reported, and never substituted for what was registered.

WHAT THIS ADDS TO `bpmd_run.py`, WHICH ALREADY RUNS BPMD:

1. **The cohort and the anchor.** `bpmd_run` runs whatever poses are on disk.
   This runs the stratified cohort plus the crystallographic positives, and
   carries the GROUP LABEL through to every row, because the whole experiment is
   a between-group comparison and a result table without the grouping is not one.

2. **Tier 1, which was being thrown away.** `gromacs_explicit` applies NO
   position restraints during NVT/NPT, so 300 ps of unrestrained dynamics runs
   before any bias is applied. The warhead's displacement across that window is a
   stability measurement that every previous run paid for and discarded. It is
   read here from `npt.gro` -- the exact frame production starts from -- against
   the docked pose in the SAME frame.

3. **Tier 1 qualifies tier 2 rather than merely preceding it.** With
   `reuse_equilibration`, tier 2 starts from the very frames tier 1 measured, so
   the two tiers describe one trajectory. A molecule whose warhead moved 3 A
   during equilibration is not having its DOCKED pose tested by BPMD, and the
   tier-1 column is what makes that visible instead of leaving a stability score
   standing on its own.

THE RECEPTOR IS ASSERTED, NOT ASSUMED. `mmgbsa.RECEPTOR_PDB` still defaults to
6VAJ, which D0059 invalidated, and 6VAJ's Cys113 SG sits 48.6 A from 3IKD's. A
run against the wrong receptor completes normally and reports plausible
distances, so the SG's coordinates are checked against the 3IKD value before any
GPU time is spent -- see `assert_receptor`.

FAIR USE. One GPU per shard, `nice -n 19`, and never GPU 0 or 7 (other people's
jobs). Four shards is the whole allowance:

    for s in 0 1 2 3; do
      tmux new-window -t elevate -d "nice -n 19 \\
        ~/.micromamba/envs/dwi_reactive/bin/python scripts/elevation_run.py \\
        --tier 1 --shard $s --n-shards 4 --gpu ${GPUS[$s]} 2>&1 | tee t1_s$s.log"
    done
"""

from __future__ import annotations

import argparse
import glob
import logging
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import bpmd                            # noqa: E402
from shared import gromacs_explicit as gx          # noqa: E402
from shared import mmgbsa as mg                    # noqa: E402
from shared import outputs as sout                 # noqa: E402
import bpmd_run as br                              # noqa: E402
import nac_screen as ns                            # noqa: E402

log = logging.getLogger("elevation")

T1 = sout.Topic("blacksmith", "elevation_tier1")
COHORT = sout.Topic("blacksmith", "elevation_cohort")

# TIER 2 GETS ITS OWN TOPIC RATHER THAN SHARING `bpmd`, AND THAT IS NOT TIDINESS.
# The bpmd topic already holds `status == ok` replicates for two molecules in
# this experiment -- t4_72f5671e89cb at 300 ps and xtal:6VAJ:QT7 at 10,000 ps --
# left by earlier protocol work. `bpmd_run.already_done()` keys on (ident,
# replicate) and knows nothing about trajectory length, so running the cohort
# through it would have SKIPPED those two molecules and quietly seated a 300 ps
# replica and a 10 ns replica inside a between-group comparison of 3 ns ones.
# A between-group comparison needs protocol consistency above all else, so the
# namespace is separated and `tier2_done` additionally matches on length.
T2 = sout.Topic("blacksmith", "elevation_tier2")

# Cys113's SG in the chemist-prepared 3IKD receptor, in the frame every docked
# pose on this branch lives in. 6VAJ's is at (-12.53, -35.87, 8.19), 48.6 A away,
# so this single check separates the two receptors unambiguously.
SG_XYZ_3IKD = np.array([13.385, 3.989, -2.040])
SG_TOL_A = 0.01

# THE ANCHOR'S MEMBERSHIP RULE, FIXED HERE BEFORE THE RUN. `crystal_positives`
# returns 15 ligands and the prereg budgets <= 8, so a rule is needed and it must
# not be "the ones that worked". Sorted by ident and take the first 8: it depends
# on nothing but the PDB accession codes, and it is reproducible by inspection.
# It happens to yield 5 chloroacetamides, 2 naphthoquinones and 1 chloroazine.
N_REF = 8
REF_GROUP = "REF_crystallographic"


class ElevationError(RuntimeError):
    """The experiment could not be run as pre-registered."""


# --------------------------------------------------------------------------
# the receptor, checked before anything is spent
# --------------------------------------------------------------------------

def assert_receptor(tmp: Path) -> dict:
    """Refuse to run unless Cys113's SG is where 3IKD's is.

    This is the cheapest possible guard against the defect that has cost this
    branch the most: `mmgbsa.prepare_receptor` falls back to 6VAJ when no
    receptor is passed, and everything downstream -- parameterisation, solvation,
    the CV, the stability score -- completes normally against the wrong protein.
    The two receptors' catalytic sulfurs are 48.6 A apart, so no tolerance
    question arises; either the coordinates match to file precision or the run is
    not the experiment that was registered.
    """
    tmp.mkdir(parents=True, exist_ok=True)
    _, cys, idx, nres = mg.prepare_receptor(tmp, receptor_pdb=br.RECEPTOR_3IKD)
    xyz = br.sg_in_pose_frame(cys, idx)
    off = float(np.linalg.norm(xyz - SG_XYZ_3IKD))
    if off > SG_TOL_A:
        raise ElevationError(
            f"Cys113 SG is at {xyz.round(3).tolist()}, {off:.2f} A from 3IKD's "
            f"{SG_XYZ_3IKD.tolist()}. This is not the chemist-prepared 3IKD "
            f"receptor (6VAJ's SG is 48.6 A away). Receptor: {br.RECEPTOR_3IKD}")
    log.info("receptor OK: %s, Cys113 = residue %d of %d, SG at %s (%.4f A off)",
             br.RECEPTOR_3IKD.name, idx, nres, xyz.round(3).tolist(), off)
    return {"receptor": str(br.RECEPTOR_3IKD), "cys113_residue": idx,
            "n_residues": nres, "sg_offset_A": round(off, 4)}


# --------------------------------------------------------------------------
# the cohort
# --------------------------------------------------------------------------

def latest_cohort() -> Path:
    fs = glob.glob(str(COHORT.dir / "elevation_cohort_*.csv"))
    if not fs:
        raise ElevationError(f"no cohort under {COHORT.dir}")
    return Path(max(fs, key=lambda f: int(re.search(r"_(\d+)\.csv$", f).group(1))))


def reference_positives(by_id: dict[str, ns.Candidate]) -> list[ns.Candidate]:
    """The anchor: crystallographic Cys113 positives, by the rule fixed above."""
    xtal = sorted((c for i, c in by_id.items() if i.startswith("xtal:")),
                  key=lambda c: c.ident)
    if not xtal:
        raise ElevationError(
            "no crystallographic positives — the prereg says the between-group "
            "comparison is uninterpretable without the anchor, so this is a "
            "refusal and not a warning")
    return xtal[:N_REF]


def cohort() -> list[tuple[str, ns.Candidate]]:
    """(group, candidate) for every molecule in the experiment, anchor included.

    Ordered cohort-first, anchor-last, and NOT shuffled: shard k takes every
    k-th molecule, so an interleaved order spreads all five groups across all
    four GPUs. If a shard dies it takes a slice of every group with it rather
    than an entire group, which is the difference between a thinner comparison
    and no comparison.
    """
    src = latest_cohort()
    df = pd.read_csv(src)
    log.info("cohort %s: %d molecules in %d groups", src.name, len(df),
             df.group.nunique())
    by_id = br.candidate_index()

    out: list[tuple[str, ns.Candidate]] = []
    missing = []
    for r in df.itertuples():
        c = by_id.get(r.ident)
        if c is None:
            missing.append(r.ident)
            continue
        out.append((r.group, c))
    if missing:
        raise ElevationError(
            f"{len(missing)} cohort molecules are not in the candidate index "
            f"({missing[:5]}); the cohort and the candidate set disagree about "
            "what these molecules are")

    ref = reference_positives(by_id)
    log.info("anchor: %d crystallographic positives (%s)", len(ref),
             ", ".join(c.ident for c in ref))
    out += [(REF_GROUP, c) for c in ref]
    return out


# --------------------------------------------------------------------------
# tier 1: did the docked pose survive unrestrained equilibration?
# --------------------------------------------------------------------------

def read_gro(path: Path) -> tuple[list[tuple[str, str]], np.ndarray, np.ndarray]:
    """((residue name, atom name) per atom, coordinates in nm, box in nm).

    Fixed-width GROMACS .gro, parsed by COLUMN and not by `split()`: atom names
    and residue names run together once a residue index reaches five digits, and
    a whitespace split then silently shifts every coordinate on the line.
    """
    lines = path.read_text().splitlines()
    if len(lines) < 3:
        raise ElevationError(f"{path} is not a .gro file")
    n = int(lines[1].strip())
    names, xyz = [], np.empty((n, 3), dtype=float)
    for i in range(n):
        l = lines[2 + i]
        names.append((l[5:10].strip(), l[10:15].strip()))
        xyz[i] = (float(l[20:28]), float(l[28:36]), float(l[36:44]))
    box = np.array([float(x) for x in lines[2 + n].split()])
    return names, xyz, box


def distance_nm(gro: Path, i0: int, j0: int, want: tuple[str, str]) -> dict:
    """The i0-j0 distance in nm, with the atom identities checked and PBC handled.

    TWO THINGS ARE CHECKED RATHER THAN ASSUMED.

    *The serials still name the right atoms.* They were derived from
    `solv.prmtop` via parmed; this file was written by GROMACS several stages
    later. The ordering does survive -- but "does survive" is a claim, and the
    atom names are right there in the file, so it is verified per read. An
    integer carried across a file-format boundary is exactly the shape of defect
    this project keeps finding.

    *Periodic images.* GROMACS may wrap a molecule across the box between
    frames, and a wrapped ligand's raw distance to SG is the box length minus the
    real one -- a number that is large, finite and completely wrong. The minimum
    image is taken, and any disagreement with the raw distance is REPORTED, since
    over 300 ps of equilibration a genuine wrap means the ligand left the pocket
    and the reader should see that rather than a quietly corrected number.
    """
    names, xyz, box = read_gro(gro)
    for k, idx in ((0, i0), (1, j0)):
        if not 0 <= idx < len(names):
            raise ElevationError(f"serial {idx} outside {gro.name}'s {len(names)} atoms")
        if names[idx][1] != want[k]:
            raise ElevationError(
                f"{gro.name} atom {idx} is named {names[idx][1]!r}, expected "
                f"{want[k]!r} — the serials do not name the CV's atoms in this file")
    if box.size >= 6 and float(np.abs(box[3:]).max()) > 1e-6:
        raise ElevationError(
            f"{gro.name} has a triclinic box {box.tolist()}; the orthorhombic "
            "minimum image used here does not apply to it")
    v = xyz[i0] - xyz[j0]
    raw = float(np.linalg.norm(v))
    mic = float(np.linalg.norm(v - box[:3] * np.round(v / box[:3])))
    return {"distance_nm": mic, "raw_distance_nm": raw,
            "pbc_wrapped": bool(abs(raw - mic) > 1e-4)}


def tier1_pose(group: str, cand: ns.Candidate, *, replicates: int, gpu: int,
               threads: int, nrun: int, dock_gpu: str, allow_redock: bool,
               on_row) -> list[dict]:
    """Equilibrate `replicates` replicas of one pose and measure the warhead's drift.

    Production is never reached: `stop_after="npt"` ends each replica at the
    frame production would have started from. The .mdp files it leaves behind are
    the ones tier 2 reuses, so this is not a separate run that tier 2 repeats.
    """
    base = {"group": group, "ident": cand.ident,
            "warhead_class": cand.warhead_class, "mechanism": cand.mechanism,
            "approach": cand.label}
    try:
        prep = br.prepare_pose(cand, nrun=nrun, gpu=dock_gpu,
                               allow_redock=allow_redock)
    except Exception as exc:                                  # noqa: BLE001
        row = {**base, "replicate": 0,
               "status": f"setup failed: {type(exc).__name__}: {str(exc)[:200]}"}
        log.warning("  SETUP FAILED: %s", row["status"])
        on_row(row)
        return [row]

    meta = {k: v for k, v in prep.items()
            if isinstance(v, (int, float, str, bool, type(None)))}
    # The docked pose's distance IN THE MD FRAME. Not `nac_distance_A`, which was
    # measured against a docking-time FLEXIBLE Cys113 sidechain and is a distance
    # to a different sulfur position; that is carried alongside, for reference.
    d_dock = float(prep["start_distance_nm"])
    want = (prep["warhead_atom_name"], prep["sg_atom_name"])
    rows = []
    for k in range(1, replicates + 1):
        row = {**base, **meta, "replicate": k, "status": "ok",
               "docked_distance_nm": round(d_dock, 4)}
        try:
            t0 = time.time()
            res = gx.run_pipeline(prep["wd"], prep["md_wd"], gpu_id=gpu,
                                  threads=threads, replicate=k,
                                  candidate_id=cand.ident, stop_after="npt",
                                  reuse_equilibration=True)
            rep = Path(res["equilibration_dir"])
            d = distance_nm(rep / "npt.gro", prep["warhead_serial0"],
                            prep["sg_serial0"], want)
            delta = d["distance_nm"] - d_dock
            row.update({
                "npt_distance_nm": round(d["distance_nm"], 4),
                "delta_nm": round(delta, 4),
                "abs_delta_nm": round(abs(delta), 4),
                "pbc_wrapped": d["pbc_wrapped"],
                "start_in_window": bool(bpmd.NAC_MIN_NM <= d_dock <= bpmd.NAC_MAX_NM),
                "npt_in_window": bool(bpmd.NAC_MIN_NM <= d["distance_nm"] <= bpmd.NAC_MAX_NM),
                "equilibration_ps": gx.NVT_PS + gx.NPT_PS,
                "velocity_seed": res["velocity_seed"],
                "reused": all(s.get("reused") for n, s in res["stages"].items()
                              if n in ("nvt", "npt")),
                "seconds": round(time.time() - t0, 1),
            })
            log.info("  rep%d: %.3f -> %.3f nm (delta %+.3f), %s",
                     k, d_dock, d["distance_nm"], delta,
                     "still in the NAC window" if row["npt_in_window"] else "outside")
        except Exception as exc:                              # noqa: BLE001
            row["status"] = f"failed: {type(exc).__name__}: {str(exc)[:200]}"
            log.warning("  rep%d FAILED: %s", k, row["status"])
        rows.append(row)
        on_row(row)
    return rows


# --------------------------------------------------------------------------
# resumability
# --------------------------------------------------------------------------

def _chunks(d: Path, pattern: str) -> list[str]:
    def key(f: str) -> tuple[int, int]:
        m = re.search(r"_s(\d+)_(\d+)\.csv$", f)
        return (int(m.group(1)), int(m.group(2))) if m else (0, 0)
    return sorted(glob.glob(str(d / pattern)), key=key)


def tier1_done() -> set[tuple[str, int]]:
    done = set()
    for f in _chunks(T1.dir, "elevation_t1_s*.csv"):
        try:
            df = pd.read_csv(f)
        except Exception as exc:                              # noqa: BLE001
            log.warning("unreadable chunk %s: %s", Path(f).name, exc)
            continue
        if not {"ident", "replicate", "status"} <= set(df.columns):
            continue
        ok = df[df.status == "ok"]
        done.update(zip(ok.ident.astype(str), ok.replicate.astype(int)))
    return done


def load_tier1() -> pd.DataFrame:
    fs = _chunks(T1.dir, "elevation_t1_s*.csv")
    if not fs:
        return pd.DataFrame()
    df = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    return df.drop_duplicates(["ident", "replicate"], keep="last")


def tier2_done(production_ps: float) -> set[tuple[str, int]]:
    """Replicates already run AT THIS TRAJECTORY LENGTH, and no others.

    Length is part of the identity of a replicate, not a detail of it. A replica
    run for 300 ps and one run for 3 ns are not two samples of the same quantity
    -- the shorter one was never given the time to escape -- and pooling them
    would make the comparison partly a record of which replicates happened to be
    lying around.
    """
    done = set()
    for f in _chunks(T2.dir, "elevation_t2_s*.csv"):
        try:
            df = pd.read_csv(f)
        except Exception as exc:                              # noqa: BLE001
            log.warning("unreadable chunk %s: %s", Path(f).name, exc)
            continue
        if not {"ident", "replicate", "status", "production_ps"} <= set(df.columns):
            continue
        ok = df[(df.status == "ok")
                & (df.production_ps.astype(float) - production_ps).abs().lt(1e-6)]
        done.update(zip(ok.ident.astype(str), ok.replicate.astype(int)))
    return done


# --------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--tier", type=int, choices=(1, 2), default=1)
    ap.add_argument("--replicates", type=int, default=3)
    ap.add_argument("--production-ps", type=float, default=3000.0,
                    help="tier 2 only; SHORT by design — see the report")
    ap.add_argument("--gpu", type=int, default=1,
                    help="ONE GPU, and never 0 or 7 (other people's jobs)")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--dock-gpu", default=None)
    ap.add_argument("--nrun", type=int, default=200,
                    help="docking runs, for a molecule with no exported pose")
    ap.add_argument("--no-redock", action="store_true")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--limit", type=int, default=None,
                    help="how many molecules; 1 verifies the chain end to end")
    ap.add_argument("--only", default=None, metavar="IDENT")
    ap.add_argument("--chunk", type=int, default=6)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format=f"%(levelname)s [t{args.tier} s{args.shard}] %(message)s")

    if args.gpu in (0, 7):
        raise SystemExit(f"GPU {args.gpu} belongs to someone else's job")

    assert_receptor(Path("/tmp") / f"elev_receptor_check_{args.shard}")
    if args.tier == 2:
        log.info("PLUMED kernel: %s", gx.plumed_kernel())

    mols = cohort()
    if args.only:
        mols = [(g, c) for g, c in mols if c.ident == args.only]
        if not mols:
            raise SystemExit(f"{args.only} is not in the cohort")
    mols = [m for k, m in enumerate(mols) if k % args.n_shards == args.shard]
    if args.limit:
        mols = mols[:args.limit]
    log.info("%d molecules on this shard, %d replicas each", len(mols),
             args.replicates)

    dock_gpu = args.dock_gpu or str(args.gpu)

    if args.tier == 1:
        done = tier1_done()
        chunk = br._ChunkWriter(T1, f"elevation_t1_s{args.shard}", args.chunk)
        for i, (group, cand) in enumerate(mols, 1):
            todo = [k for k in range(1, args.replicates + 1)
                    if (cand.ident, k) not in done]
            if not todo:
                log.info("[%d/%d] %s: all replicas already measured", i,
                         len(mols), cand.ident)
                continue
            log.info("[%d/%d] %s (%s, %s)", i, len(mols), cand.ident, group,
                     cand.mechanism)
            tier1_pose(group, cand, replicates=args.replicates, gpu=args.gpu,
                       threads=args.threads, nrun=args.nrun, dock_gpu=dock_gpu,
                       allow_redock=not args.no_redock, on_row=chunk.add)
        chunk.flush()
        return

    # Tier 2 delegates to the BPMD driver rather than reimplementing it, and
    # reuses the equilibration tier 1 already ran, so the biased trajectory
    # continues the one tier 1 measured instead of starting from a fresh one.
    chunk = br._ChunkWriter(T2, f"elevation_t2_s{args.shard}", args.chunk)
    done = tier2_done(args.production_ps)
    for i, (group, cand) in enumerate(mols, 1):
        todo = [k for k in range(1, args.replicates + 1)
                if (cand.ident, k) not in done]
        if not todo:
            log.info("[%d/%d] %s: all replicas already done", i, len(mols),
                     cand.ident)
            continue
        log.info("[%d/%d] %s (%s, %s)", i, len(mols), cand.ident, group,
                 cand.mechanism)
        br.run_pose(cand, replicates=args.replicates,
                    production_ps=args.production_ps, gpu=args.gpu,
                    threads=args.threads, nrun=args.nrun, dock_gpu=dock_gpu,
                    allow_redock=not args.no_redock,
                    on_row=lambda r, g=group: chunk.add({**r, "group": g}),
                    reuse_equilibration=True)
    chunk.flush()


if __name__ == "__main__":
    main()
