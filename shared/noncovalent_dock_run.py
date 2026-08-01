"""
Purpose: The non-covalent docking run, shared by T_1 and T_2.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-28
Input: an approach's latest frame (survivors carry `canonical_smiles`)
Output: the frame with `vina_affinity` merged back; ligand + pose PDBQTs

THE NON-COVALENT COUNTERPART TO `covalent_dock_run`. T_1 and T_2 are the two
reversible approaches and the integration phase pools them on one plot, so they
must dock through one engine, one box and one search depth for the same reason
T_3 and T_4 must. The loop lives here; each approach supplies only its identity.

ENGINE: VINA-GPU AT search_depth 20 (D0017). Not CPU Vina, and not gnina. D0017
adopted Vina-GPU only after a search_depth sweep showed the depth-10
disagreement with CPU Vina was search CONVERGENCE rather than implementation —
every metric improved monotonically with depth, and at 20 the enrichment ROC-AUC
drift is 0.005 against a 0.10 threshold. Below 20 the adoption evidence does not
hold, so the depth is pinned here rather than exposed as a tuning knob.

THE GATE ONLY TRANSFERS IF THE ENGINE MATCHES. The non-covalent enrichment gate
measured ROC-AUC 0.535 (D0016) using this receptor and this box. That verdict
describes T_1 and T_2 only insofar as they dock the same way, which is why the
expanded box is read from the shared receptor directory rather than restated.

WHY THE EXPANDED BOX. `box_expanded.json` (26 A) is declared `used_by: [t1, t2]`
against the covalent box's 20 A. T_1 and T_2 place whole molecules with no
anchor at Cys113, so they need room the covalent approaches do not.

GPUS ARE ALLOCATED EXPLICITLY, NOT INFERRED. Vina-GPU runs in virtual-screening
mode: one process, one GPU, a whole ligand directory. The caller passes the
device id. This is deliberate — `covalent_dock_run.select_gpus` infers idleness
from a 1024 MiB memory threshold, and gnina occupies only ~500 MiB, so a job
that auto-selected would land on GPUs already running a covalent dock while
idle ones sat unused.
"""

from __future__ import annotations

import json
import logging
import os
import re
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

from . import io as dio

RDLogger.DisableLog("rdApp.*")

log = logging.getLogger(__name__)

DATA_ROOT = Path("/data/lab_vm/append_only/inhibition")
RECEPTOR_ROOT = Path("/data/lab_vm/immutable/inhibition/receptor")
RECEPTOR_PDBQT = RECEPTOR_ROOT / "6VAJ_prepared.pdbqt"
BOX_EXPANDED = RECEPTOR_ROOT / "box_expanded.json"
VINA_GPU = Path("/data/lab_vm/envs/dwi_vinagpu/bin/vina-gpu")
OBABEL = "/data/lab_vm/envs/dwi_cheminf/bin/obabel"

SEARCH_DEPTH = 20      # D0017 — the adoption evidence does not hold below this
THREADS = 8000
NICE = 19

# LIGANDS ARE PROTONATED FOR pH 7.4, NOT DOCKED AS DRAWN (issue #5, @tt8804).
#
# Until 2026-07-31 every ligand was built with `Chem.AddHs()` on the neutral
# SMILES and converted with a bare `obabel` call -- docked exactly as written.
# That is wrong for most of what this project generates. Measured over the five
# T_2 seed neighbourhoods:
#
#   guo_pfizer    83.4% carry a phosphate    -> dianion at pH 7.4
#   atra          86.6% a carboxylic acid    -> anion
#   du_xu         85.5% a carboxylic acid    -> anion
#   potter_astex  86.7% a basic amine        -> cation
#   liu_2024_c3   largely neutral
#
# Each neighbourhood inherits its seed's ionizable group at ~85%, so comparing
# neutral-form scores across seeds compares four charge states all pretending to
# be neutral.
#
# WHAT THIS CHANGES AND WHAT IT DOES NOT. Vina carries no electrostatic term,
# and obabel writes all-zero partial charges into the PDBQT either way -- both
# verified here. So this does NOT make Vina "see" the charge, and the concern
# that a basic Lys63/Arg68/Arg69 subsite electrostatically rewards anions is not
# reachable through Vina's scoring function at all. What it changes is ATOM
# TYPING: deprotonating a phosphate drops two polar hydrogens (30 -> 28 atoms on
# 21b), and Vina's H-bond term counts donors from those `HD` atoms. Docking the
# neutral acid presents H-bond donors that do not exist at physiological pH.
LIGAND_PH = 7.4

# THE PREP CACHE MUST BE KEYED ON THE PROTONATION, NOT JUST THE CANDIDATE.
#
# `_prepare_one` short-circuits on `if out.is_file()`, and the path was
# `ligands/<candidate_id>.pdbqt` -- carrying nothing about how it was prepared.
# Adding `-p 7.4` without changing the path would silently re-dock the SAME
# neutral files and report success; 399 MB (liu), 231 MB (du_xu) and 45 MB
# (atra) of stale neutral PDBQTs were already on disk when this landed.
#
# Same defect as the class-pool cache keyed on `class_id` alone, which served a
# stale 3-molecule pool after its query was relaxed. A cache keyed on less than
# its inputs is a cache that lies. The tag goes in the directory name so the old
# and new sets coexist -- the tree is append-only and nothing is deleted.
LIGAND_PREP_TAG = f"ph{LIGAND_PH:g}"

# Ligand prep is embarrassingly parallel and would happily take all 224 cores.
# The project budget lives in shared/compute.py (raised to 50 by @mhallet on
# 2026-07-28).
from . import compute                             # noqa: E402
MAX_PREP_WORKERS = compute.MAX_CPU_WORKERS

_RESULT_RE = re.compile(r"REMARK VINA RESULT:\s*(-?\d+\.\d+)")

# ALL NINE MODES, NOT JUST THE FIRST (issue #10, @tt8804).
#
# Vina writes every mode it reports into the output PDBQT, each with its own
# affinity and its RMSD back to the best mode. `collect_scores` regex-searched
# for the FIRST `REMARK VINA RESULT` line and returned one float per ligand, so
# 8 of every 9 poses we computed were written to disk and never read -- across
# all ~56,000 candidates.
#
# That matters because mode 1 is not meaningfully better than mode 2. Measured
# over 400 completed du_xu ligands: the median affinity gap between mode 2 and
# mode 1 is 0.10 kcal/mol, and 99% of ligands have mode 2 within 0.5 kcal/mol.
# Vina's own reported RMSE is ~2-3 kcal/mol, so selecting mode 1 over mode 2
# happens on a margin twenty to thirty times smaller than the error bar of the
# function producing it.
#
# WHAT THESE COLUMNS ARE NOT. `mode_rmsd_nn` is the smallest RMSD from any
# other mode back to the best one, and it is NOT a convergence statistic. Vina
# diversifies its reported modes with a minimum-RMSD floor before writing them,
# so a cluster population computed over these nine would measure the output
# filter rather than agreement. The honest version is replicate runs with
# independent seeds. These columns are DESCRIPTIVE LABELS; nothing downstream
# should rank on them until a rule is scored against D0046's 80 redock cases.
_MODE_RE = re.compile(
    r"REMARK VINA RESULT:\s*(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)")

# Columns this module owns on the frame. Listed so the merge can drop exactly
# what it is about to supply -- see the drop-list defect in
# `docs/how_this_project_breaks.md` (catalogue #5), where a hand-maintained
# list that had to stay in sync with code five lines away did not.
MODE_COLS = ("vina_affinity", "vina_n_modes", "vina_mode2_gap",
             "vina_mode_rmsd_nn", "vina_affinity_spread")


def _prepare_one(job: dict) -> dict:
    """Embed one ligand and convert it to PDBQT. Runs in a worker process."""
    out = Path(job["ligand_dir"]) / f"{job['candidate_id']}.pdbqt"
    if out.is_file():
        return {"candidate_id": job["candidate_id"], "ok": True}
    sdf = out.with_suffix(".sdf")
    try:
        m = Chem.MolFromSmiles(job["smiles"])
        if m is None:
            return {"candidate_id": job["candidate_id"], "ok": False,
                    "error": "unparseable"}
        m = Chem.AddHs(m)
        if AllChem.EmbedMolecule(m, randomSeed=42) != 0:
            return {"candidate_id": job["candidate_id"], "ok": False,
                    "error": "embed_failed"}
        try:
            AllChem.MMFFOptimizeMolecule(m, maxIters=500)
        except Exception:  # noqa: BLE001 - MMFF can fail; the pose is still usable
            pass
        w = Chem.SDWriter(str(sdf))
        w.write(m)
        w.close()
        # -p protonates for the given pH: it strips the acidic hydrogens and
        # adds basic ones, so the PDBQT carries the physiological ionization
        # rather than the drawn one. See LIGAND_PH.
        conv = subprocess.run([OBABEL, str(sdf), "-O", str(out),
                               "-p", str(LIGAND_PH)],
                              capture_output=True, text=True, timeout=120)
        sdf.unlink(missing_ok=True)
        if conv.returncode != 0 or not out.is_file():
            return {"candidate_id": job["candidate_id"], "ok": False,
                    "error": "pdbqt_conversion_failed"}
        return {"candidate_id": job["candidate_id"], "ok": True}
    except Exception as exc:  # noqa: BLE001 - one bad ligand must not end prep
        return {"candidate_id": job["candidate_id"], "ok": False,
                "error": str(exc)[:200]}


def prepare_ligands(survivors: pd.DataFrame, ligand_dir: Path,
                    smiles_col: str = "canonical_smiles") -> list[dict]:
    """Embed + convert every survivor, capped at the lab's CPU limit."""
    ligand_dir.mkdir(parents=True, exist_ok=True)
    jobs = [{"candidate_id": r["candidate_id"], "smiles": r[smiles_col],
             "ligand_dir": str(ligand_dir)}
            for _, r in survivors.iterrows()]
    results = []
    with ProcessPoolExecutor(max_workers=MAX_PREP_WORKERS) as ex:
        for i, res in enumerate(ex.map(_prepare_one, jobs, chunksize=16), 1):
            results.append(res)
            if i % 500 == 0:
                log.info("prepared %d/%d ligands", i, len(jobs))
    bad = [r for r in results if not r["ok"]]
    if bad:
        log.warning("%d ligand(s) failed preparation and will not be docked; "
                    "first: %s", len(bad), bad[0])
    return results


# THE TIMEOUT SCALES WITH THE POOL, BECAUSE A FLAT ONE DESTROYED A DAY OF WORK.
#
# `run_vina_gpu` passed `timeout=86400` — a flat 24 h, correct when a pool was
# ATRA's 1,882 molecules and silently wrong once it was liu_2024_c3's 16,806.
# On 2026-08-01 that seed ran for exactly 24 h, raised `TimeoutExpired`, and
# lost EVERYTHING: Vina-GPU writes its poses only at the end of a
# virtual-screening run, so the output directory held 0 of 16,806 files and no
# frame was written. The projection for that pool was ~38 h, so it could never
# have finished under the cap — and nothing said so before the day was spent.
#
# Measured rate on this box: du_xu docked 9,736 ligands in 9 h 56 m and
# guo_pfizer 8,670 in 9 h 04 m, i.e. ~3.7 s/ligand. The safety factor is large
# because the rate is NOT uniform: cost climbs steeply with rotatable bonds,
# and liu's pool averages 10.65 against du_xu's 4.81.
SECONDS_PER_LIGAND = 3.7        # measured, du_xu + guo_pfizer
TIMEOUT_SAFETY_FACTOR = 4       # covers the flexibility spread, not just count
MIN_TIMEOUT_S = 3600


def vina_timeout_s(n_ligands: int) -> int:
    """A deadline proportional to the work, floored for tiny smoke runs."""
    return max(MIN_TIMEOUT_S,
               int(n_ligands * SECONDS_PER_LIGAND * TIMEOUT_SAFETY_FACTOR))


def run_vina_gpu(ligand_dir: Path, out_dir: Path, gpu: int) -> float:
    """Vina-GPU in virtual-screening mode on ONE explicitly chosen GPU."""
    out_dir.mkdir(parents=True, exist_ok=True)
    b = json.loads(BOX_EXPANDED.read_text())
    cmd = [str(VINA_GPU),
           "--receptor", str(RECEPTOR_PDBQT),
           "--ligand_directory", str(ligand_dir),
           "--output_directory", str(out_dir),
           "--center_x", str(b["center_x"]), "--center_y", str(b["center_y"]),
           "--center_z", str(b["center_z"]),
           "--size_x", str(b["size_x"]), "--size_y", str(b["size_y"]),
           "--size_z", str(b["size_z"]),
           "--thread", str(THREADS), "--search_depth", str(SEARCH_DEPTH)]
    env = dict(os.environ)
    env["CUDA_VISIBLE_DEVICES"] = str(gpu)
    env["GPU_DEVICE_ORDINAL"] = str(gpu)      # OpenCL honours this, not CUDA_*
    n_ligands = len(list(ligand_dir.glob("*.pdbqt")))
    timeout_s = vina_timeout_s(n_ligands)
    # STATED UP FRONT, NOT DISCOVERED AT THE END. The run that was killed by
    # the old flat cap gave no indication it was on a deadline it could not
    # meet until 24 h had already been spent.
    log.info("Vina-GPU on GPU %d, depth %d, %d ligands; deadline %.1f h "
             "(all-or-nothing: poses are written only on completion)",
             gpu, SEARCH_DEPTH, n_ligands, timeout_s / 3600)
    t0 = time.time()
    try:
        p = subprocess.run(cmd, capture_output=True, text=True,
                           timeout=timeout_s, env=env)
    except subprocess.TimeoutExpired:
        log.error("Vina-GPU exceeded its %.1f h deadline on %d ligands and was "
                  "killed; NOTHING was written, because poses land only at the "
                  "end of a virtual-screening run. Split the pool or raise "
                  "TIMEOUT_SAFETY_FACTOR before retrying.",
                  timeout_s / 3600, n_ligands)
        raise
    dt = time.time() - t0
    (out_dir / "vina_gpu_stdout.log").write_text(
        p.stdout + "\n--- stderr ---\n" + p.stderr)
    if p.returncode != 0:
        raise RuntimeError(f"Vina-GPU failed ({p.returncode}); see "
                           f"{out_dir/'vina_gpu_stdout.log'}\n"
                           f"{(p.stderr or p.stdout)[-1500:]}")
    log.info("Vina-GPU finished in %.1f s", dt)
    return dt


def collect_scores(out_dir: Path) -> dict[str, float]:
    """Best score per ligand from Vina-GPU's output PDBQTs.

    Kept as the single-value accessor; `collect_modes` is the full parse.
    """
    return {r["candidate_id"]: r["vina_affinity"]
            for r in collect_modes(out_dir).to_dict("records")}


def parse_modes(text: str) -> list[tuple[float, float, float]]:
    """(affinity, rmsd_lb, rmsd_ub) for every mode in one output PDBQT.

    Vina reports the RMSDs relative to the BEST mode, so mode 1's are 0.0 by
    construction and carry no information.
    """
    return [(float(a), float(lb), float(ub))
            for a, lb, ub in _MODE_RE.findall(text)]


def collect_modes(out_dir: Path) -> pd.DataFrame:
    """Every mode of every ligand, summarised one row per candidate.

    See MODE_COLS for why the other eight poses are worth reading and what
    these numbers may NOT be used for.
    """
    rows = []
    for f in out_dir.glob("*.pdbqt"):
        modes = parse_modes(f.read_text(errors="replace"))
        if not modes:
            continue
        aff = [m[0] for m in modes]
        # Nearest neighbour back to the best mode: the smallest lower-bound
        # RMSD among the OTHER modes. NaN when Vina reported a single mode --
        # absent, not zero, because zero would read as "perfectly converged".
        nn = min((m[1] for m in modes[1:]), default=float("nan"))
        rows.append({
            "candidate_id": f.stem.replace("_out", ""),
            "vina_affinity": aff[0],
            "vina_n_modes": len(modes),
            "vina_mode2_gap": (aff[1] - aff[0]) if len(aff) > 1 else float("nan"),
            "vina_mode_rmsd_nn": nn,
            "vina_affinity_spread": max(aff) - min(aff),
        })
    df = pd.DataFrame(rows, columns=list(("candidate_id", *MODE_COLS)))
    if not df.empty:
        log.info("parsed %d ligands; %.2f modes each on average, median "
                 "mode-2 gap %.3f kcal/mol", len(df),
                 df["vina_n_modes"].mean(), df["vina_mode2_gap"].median())
    return df


def run(*, experiment: str, approach: str, frame_prefix: str, gpu: int,
        limit: int | None = None):
    """Dock one approach's survivors and merge the result onto its frame."""
    os.nice(NICE)
    work = DATA_ROOT / experiment / "docking"
    # Ligand dir carries the prep tag so a protonation change cannot be served
    # from a cache built under different settings -- see LIGAND_PREP_TAG.
    ligand_dir = work / f"ligands_{LIGAND_PREP_TAG}"
    out_dir = work / f"poses_{LIGAND_PREP_TAG}"

    frame_path = dio.latest(DATA_ROOT / experiment, frame_prefix, ".parquet")
    if frame_path is None:
        raise SystemExit(f"no {frame_prefix} frame for {experiment}")
    df = dio.read_frame(frame_path)

    survivors = df[df["rejected_at"].isna()].copy()
    if limit:
        survivors = survivors.head(limit)
    log.info("[%s] %d survivors (of %d in the frame)",
             approach, len(survivors), len(df))

    prep = prepare_ligands(survivors, ligand_dir)
    n_ready = sum(1 for r in prep if r["ok"])
    log.info("[%s] %d/%d ligands ready to dock", approach, n_ready, len(prep))
    if not n_ready:
        raise SystemExit("no ligand survived preparation")

    elapsed = run_vina_gpu(ligand_dir, out_dir, gpu)
    return merge_poses_onto_frame(
        experiment=experiment, approach=approach, frame_prefix=frame_prefix,
        out_dir=out_dir, elapsed=elapsed, gpu=gpu, limit=limit,
        df=df, survivors=survivors)


def merge_poses_onto_frame(*, experiment: str, approach: str,
                           frame_prefix: str, out_dir: Path, elapsed: float,
                           gpu, limit: int | None = None,
                           df: pd.DataFrame | None = None,
                           survivors: pd.DataFrame | None = None):
    """Parse a pose directory, merge onto the frame, write it.

    SPLIT OUT SO A CHUNKED RUN USES THE SAME CODE. A pool too large for one
    GPU is docked as N chunks writing into ONE pose directory, and the frame is
    written once afterwards. That collection step must be the identical merge
    -- with the same derived drop list and the same `_x`/`_y` assertion -- and
    not a second implementation in a driver script, because two code paths that
    both "merge the docking results" is precisely how the covalent and GROMACS
    frames ended up with suffixed columns nobody noticed.

    `df`/`survivors` are passed when the caller already loaded them; a chunked
    run has not, and reads the latest frame itself.
    """
    if df is None:
        frame_path = dio.latest(DATA_ROOT / experiment, frame_prefix, ".parquet")
        if frame_path is None:
            raise SystemExit(f"no {frame_prefix} frame for {experiment}")
        df = dio.read_frame(frame_path)
    if survivors is None:
        survivors = df[df["rejected_at"].isna()] if "rejected_at" in df else df

    scored = collect_modes(out_dir)

    # Drop stale columns BEFORE merging, so a re-run does not produce
    # vina_affinity_x / _y and silently lose the column downstream. The list is
    # DERIVED FROM THE MERGE rather than written out, because a hand-maintained
    # drop list is exactly how the same defect reached D1_21/D2_21 via
    # `merge_gromacs_results.py` -- it omitted columns added five lines above
    # it. Suffixed survivors go too, so an already-damaged frame heals.
    incoming = [c for c in scored.columns if c != "candidate_id"]
    stale = {*incoming, *(f"{c}{s}" for c in incoming for s in ("_x", "_y"))}
    drop = [c for c in stale if c in df.columns]
    if drop:
        log.info("dropping %d stale column(s) before merge: %s",
                 len(drop), ", ".join(sorted(drop)))
        df = df.drop(columns=drop)
    merged = df.merge(scored, on="candidate_id", how="left")
    suffixed = [c for c in merged.columns if c.endswith(("_x", "_y"))]
    if suffixed:
        raise RuntimeError(
            f"merge produced suffixed columns {suffixed}; the drop list did "
            "not cover what the merge supplied")
    if len(merged) != len(df):
        raise RuntimeError(
            f"merge changed row count {len(df)} -> {len(merged)}; duplicate "
            "candidate_id in the frame or the results")

    n_docked = int(merged["vina_affinity"].notna().sum())

    # See covalent_dock_run: a --limit run must not become the latest frame.
    if limit:
        log.warning("--limit %d: NOT writing a frame. A partial run must not "
                    "become the latest frame for the next stage to read.", limit)
        return merged, None, survivors, n_docked, elapsed

    out = dio.write_full_frame(
        merged, approach=approach, experiment=experiment,
        stage=f"{approach}_noncovalent_dock",
        params={"engine": "Vina-GPU 2.1 (D0017)",
                "search_depth": SEARCH_DEPTH,
                "threads": THREADS,
                "box": BOX_EXPANDED.name,
                "gpu": gpu,
                # In the manifest so a frame states its own protonation rather
                # than leaving a reader to infer it from the run date.
                "ligand_ph": LIGAND_PH,
                "ligand_prep_tag": LIGAND_PREP_TAG,
                "rank_metric": "vina_affinity (kcal/mol, lower better)",
                # Descriptive only. `vina_mode_rmsd_nn` reflects Vina's mode
                # diversification floor, not pose convergence -- see MODE_COLS.
                "mode_columns": list(MODE_COLS),
                "enrichment_caveat":
                    "D0016: non-covalent enrichment on this pocket is ROC-AUC "
                    "0.535 (CI 0.215-0.855, EF1% 0.0). These rankings are "
                    "weakly supported and must be shown with that verdict.",
                "elapsed_s": round(elapsed, 1),
                "n_docked": n_docked},
        inputs={"frame": frame_path, "receptor": RECEPTOR_PDBQT})
    return merged, out, survivors, n_docked, elapsed
