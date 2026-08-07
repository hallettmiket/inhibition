"""
Purpose: re-run the geometric screen PERSISTING per-pose geometry and gnina scores, so the score can be repaired without docking again.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06
Input: all warhead-bearing candidates (T_3 + T_4), the chemist's 3IKD
Output: 00_outputs/blacksmith/nac_v2/{poses,agg}_s<shard>_<N>.csv + top-20 poses on disk

WHY THIS RE-RUN EXISTS. The 2.0.0 screen computed per-pose distances and angles,
used them to make one summary number per molecule, and threw the poses away. 79
molecules out of 5,765 have poses on disk. So two known repairs to the ranking
score are both blocked on the same missing data, and neither can be tested
without docking again:

  1. THE SCORE IS PARTLY MEASURING OUR OWN THUMB ON THE SCALE.
     `viable_fraction` is the JOINT rate P(distance in window AND angle in
     window), while `isotropic_null` is purely orientational, P(angle). The
     numerator therefore carries a distance hit rate -- and distance is exactly
     what the reactive potential is designed to bias (D0064: the reactive
     potential is a SAMPLER, not a criterion). Part of the score is how well the
     sampler worked on a molecule rather than the molecule's own ability to
     orient. The conditional form, P(angle | distance in window), cancels it:
     the bias applies equally to the numerator and to the conditioning set.

  2. WE ARE SCORING THE WRONG TEN POSES. Measured on this receptor tonight over
     82 Pin1 crystal ligands, AutoDock's own ordering puts a sub-2 A pose first
     18.3% of the time and gnina's CNN gets 26.8%, with sampling held fixed. Both
     consensus and the geometric score are computed on AutoDock's top ten, so a
     better ordering changes which ten they see.

WHAT IS PERSISTED, chosen so no third re-run is needed:

  poses_*.csv   one row per pose for EVERY pose -- not the top 20 by energy.
                energy, distance, angle, approach, viable, and the MODE the pose
                was assigned to.

                The old top-20-by-energy window was the single largest loss in
                the pipeline. #30 measured that the crystallographic pose is
                present in the pose set 93.3% of the time at 200 runs (100% at
                500) and survives the energy cut under half the time, because
                energy places the correct pose at a rank indistinguishable from
                uniform (KS p = 0.666). Persisting only the survivors also meant
                pose splitting could not be done retrospectively at all.

                ENERGY STILL GENERATES THE POSES AND NO LONGER SELECTS THEM.
                Each of the N runs is an independent Lamarckian GA optimising
                AutoDock's scoring function, so energy decides WHERE a run lands
                -- that cannot be removed without removing the docking. What is
                removed is energy as a SELECTION criterion. `mean_energy` is
                reported per mode and consumed by nothing.
  agg_*.csv     counts over ALL poses -- n_poses, n_in_range, n_viable,
                n_viable_given_in_range -- so the conditional score can also be
                computed over the whole population, not only the retained window.
  poses/        the top-20 poses themselves, as SDF. The thing 2.0.0 discarded.

POSE ORDER IS THE JOIN KEY between geometry and gnina scores, so it is asserted
rather than assumed at every hop, exactly as `rescore_benchmark` does: counts
must match, and SDF coordinates are checked against the conformer they came from.
This project's defect catalogue is values taken by position instead of identity.
"""

from __future__ import annotations

import argparse
import glob
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import nac_criterion as nac           # noqa: E402
from shared import outputs as sout                # noqa: E402
from shared import covalent_protocol as cp        # noqa: E402
from shared import receptors as R                 # noqa: E402
from shared import pose_modes as pmod             # noqa: E402
import nac_screen as ns                           # noqa: E402
import nac_rank as nr                             # noqa: E402

log = logging.getLogger("nac-v2")
#: 2.2.0 writes to its OWN topic, not nac_v2.
#:
#: `rank_v2` globs `agg_s*_*.csv` out of whichever topic it is pointed at. The
#: 2.1.0 frames there are one row per MOLECULE at nrun=200; these are one row per
#: MODE at nrun=500, with idents `<parent>_m<k>`. Because the idents differ, a
#: dedup would NOT collapse them -- the two tables would simply concatenate into
#: one frame with two incompatible meanings and no error. Separate topics make
#: that impossible, and leave the 2.1.0 data intact for comparison.
OUT = sout.Topic("blacksmith", "nac_v3")
#: 2.2.0 poses go to their OWN directory, and this is not cosmetic.
#:
#: `write_sdf` is guarded by `if not sdf.exists()` -- correct under the
#: append-only rule, never overwrite. But every molecule already has a file in
#: `nac_v2_poses` from the 200-run screen, so pointing 2.2.0 at that directory
#: would have written NONE of the new poses and left every downstream stage
#: reading the old energy-selected, un-protonated ones. The re-dock would have
#: completed, reported success, and changed nothing.
POSE_DIR = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith/nac_v3_poses")
#: RETIRED. Every pose is persisted now; this survives only so the docstring
#: below and #23/#30 can refer to what was removed. Nothing reads it.
KEEP_TOP = 20        # noqa: F841  (retired -- see nac_screen_v2 docstring)
_COORD_TOL = 0.05


def write_sdf(mol, order: list[int], dest: Path,
              modes: list[int] | None = None) -> int:
    """Write the top-`order` conformers to one SDF, stamping `pose_rank`.

    The rank is written as a PROPERTY, not left to file position, because
    `bpmd_run.read_pose` selects a pose by its own `pose_rank` and refuses to
    take one by index -- a pose identified by where it sits in a file is a pose
    that a re-sort silently redefines.
    """
    from rdkit import Chem
    w = Chem.SDWriter(str(dest))
    n = 0
    for rank, i in enumerate(order, 1):
        mol.SetProp("pose_rank", str(rank))
        mol.SetProp("energy_rank", str(rank))
        # The mode this pose represents, so a downstream reader selects by
        # IDENTITY rather than by file position -- the same rule pose_rank
        # already follows for bpmd_run.
        if modes is not None:
            mol.SetProp("mode", str(modes[rank - 1]))
        w.write(mol, confId=i)
        n += 1
    w.close()
    return n


def gnina_scores(receptor: Path, sdf: Path, gpu: str) -> pd.DataFrame:
    env = cp.gnina_env()
    env["CUDA_VISIBLE_DEVICES"] = str(gpu)
    r = subprocess.run(
        [str(cp.GNINA_BIN), "--receptor", str(receptor), "--ligand", str(sdf),
         "--score_only", "--cnn_scoring", "rescore", "--seed", "42"],
        capture_output=True, text=True, env=env, timeout=1800)
    if r.returncode != 0:
        raise RuntimeError(f"gnina rc={r.returncode}: {r.stderr[-200:]}")
    rows, cur = [], {}
    for line in r.stdout.splitlines():
        m = re.match(r"^(Affinity|CNNscore|CNNaffinity|CNNvariance):\s+(-?[\d.]+)", line)
        if not m:
            continue
        if m.group(1) == "Affinity" and cur:
            rows.append(cur)
            cur = {}
        cur[m.group(1)] = float(m.group(2))
    if cur:
        rows.append(cur)
    return pd.DataFrame(rows)


def one(cand, rec_dir: Path, plain_rec: Path, nrun: int, gpu: str,
        do_gnina: bool) -> tuple[pd.DataFrame, list[dict]]:
    """Dock one candidate; return (per-pose rows, ONE AGGREGATE ROW PER MODE).

    2.2.0 (@tt8804): a binding mode is a candidate, not a property of one. The
    second element used to be a single dict describing the molecule; it is now a
    list, one entry per mode, each carrying an ident of `<parent>_m<k>`. Ranking,
    class stratification, selection and the GUI all read rows, so a longer table
    is all they see -- nothing downstream needed rewriting for this.
    """
    from rdkit import Chem
    work = Path(tempfile.mkdtemp(prefix="nacv2_"))
    agg = {"ident": cand.ident, "warhead_class": cand.warhead_class,
           "mechanism": cand.mechanism, "nrun": nrun}
    try:
        best = None
        for j, lig in enumerate(ns.prepare_ligand(cand, work / "lig.pdbqt")):
            dlg = ns.dock(lig, rec_dir, work / f"c{j}", nrun, gpu)
            # ONE rebuild, ONE reactive-atom match, shared by the criterion and
            # the clustering, so they can never describe different atoms.
            mol_j, match_j = ns.rebuild_and_match(dlg, cand)
            res = nac.measure_poses(mol_j, match_j, cand.mechanism,
                                    ns.sg_position(dlg))
            en = ns.pose_energies(dlg)
            if len(en) != len(res):
                raise ValueError(f"{len(en)} energies vs {len(res)} poses")
            frac = nac.viable_fraction(res)
            if best is None or frac > best[0]:
                best = (frac, res, en, dlg, mol_j, match_j)
        if best is None:
            raise ValueError("no usable poses")
        frac, res, en, dlg, mol, match = best
        # `pose_energies` returns a LIST. Every per-mode aggregate indexes it
        # with a boolean mask, which a list cannot take -- that failed 6 of 6
        # molecules in the smoke test with "only integer scalar arrays can be
        # converted to a scalar index".
        en = np.asarray(en, dtype=float)

        in_rng = np.array([nac.NAC_DIST_MIN <= r.distance <= nac.NAC_DIST_MAX
                           for r in res])
        viable = np.array([bool(r.viable) for r in res])

        # ---- split the pose cloud into binding modes ----------------------
        # Clusters on the reactive atom's position and the direction its warhead
        # faces -- never on energy (#23/#30: energy places the crystal pose at a
        # rank indistinguishable from uniform) and never on the NAC geometry
        # itself, which is the score.
        feat = pmod.features(mol, match)
        labels = pmod.split(feat)
        mode_ids = sorted(set(int(l) for l in labels) - {-1})

        # ---- EVERY pose is persisted, not the 20 best-scoring --------------
        # #30: the crystal pose is present in the pose set 93.3% of the time and
        # survives `KEEP_TOP = 20 by energy` under half the time. The 180 poses
        # 2.1.0 discarded are exactly the ones splitting needs.
        rows = pd.DataFrame([{
            "ident": cand.ident, "warhead_class": cand.warhead_class,
            "mechanism": cand.mechanism, "pose_idx": i,
            "mode": int(labels[i]),
            "energy": en[i], "energy_rank": int(np.argsort(np.argsort(en))[i]) + 1,
            "distance": res[i].distance, "angle": res[i].angle,
            "approach_angle": res[i].approach_angle,   # None for SN2 by design
            "angle_kind": res[i].angle_kind,
            "viable": viable[i], "in_range": bool(in_rng[i]),
        } for i in range(len(res))])

        # Anchoring quality per pose, computed ONCE. Used twice below: to score
        # each mode, and to choose the representative pose within it. Defined
        # here rather than beside the representative block, which is where it
        # was and which put its use before its definition.
        anchor = np.array([nac.anchor_quality(r.distance, r.angle, cand.mechanism)
                           for r in res])
        dmat = pmod.distances(feat)

        # ---- one aggregate row PER MODE -----------------------------------
        aggs = []
        for k in mode_ids:
            sel = labels == k
            ident_k = f"{cand.ident}_m{k}"
            ident_row = dict(agg)
            ident_row.update({
                "ident": ident_k, "parent_ident": cand.ident, "mode": k,
                "n_poses": int(len(res)),
                "n_poses_mode": int(sel.sum()),
                # CONSENSUS IS NOW MODE POPULATION. Not "do the top-10 by energy
                # agree" -- that read an energy-selected sample of a uniformly
                # uninformative ordering. How much of the pose cloud lands on
                # this geometry is the same idea with the energy removed, and it
                # picks the crystal pose 93.3% of the time against 60.0% for the
                # old energy window.
                "consensus": float(sel.mean()),
                "n_in_range": int((in_rng & sel).sum()),
                "n_viable": int((viable & sel).sum()),
                "n_viable_given_in_range": int((viable & in_rng & sel).sum()),
                "viable_fraction": float(viable[sel].mean()) if sel.any() else 0.0,
                "enrichment": (float(viable[sel].mean()) /
                               nac.isotropic_null(cand.mechanism)) if sel.any() else 0.0,
                "isotropic_null": nac.isotropic_null(cand.mechanism),
                "mean_energy": float(np.nanmean(en[sel])) if sel.any() else np.nan,
                "anchor_quality_max": float(np.nanmax(anchor[sel])) if sel.any() else np.nan,
                "anchor_quality_mean": float(np.nanmean(anchor[sel])) if sel.any() else np.nan,
                "status": "ok",
            })
            ident_row.update(pmod.identity(feat, labels, k))
            aggs.append(ident_row)

        if not aggs:
            # Every pose was noise. Recorded as its own status rather than as a
            # molecule with zero modes, because those are different failures: one
            # is a molecule that does nothing reproducible, the other is a bug.
            agg["status"] = "no mode above the population floor"
            agg["n_poses"] = int(len(res))
            return rows, [agg]

        # ---- persist a representative of EVERY mode ------------------------
        # The pose most central to its own mode, not its lowest-energy member.
        POSE_DIR.mkdir(parents=True, exist_ok=True)
        sdf = POSE_DIR / f"{cand.ident}.sdf"
        # THE REPRESENTATIVE IS A TYPICAL POSE FROM THE WELL-ANCHORED QUARTILE,
        # NOT THE BEST-ANCHORED ONE.
        #
        # @tt8804 warned against prioritising attack geometry over realistic
        # poses. Measured on 15 crystal complexes, picking one pose out of the
        # dominant mode:
        #
        #     ceiling: best pose in the mode              93.3%
        #     medoid of the top-25% by anchoring          33.3%   <- adopted
        #     medoid of the whole mode                    26.7%
        #     argmax anchoring                             6.7%   <- was here
        #
        # Anchoring is informative -- across the whole pose population it
        # correlates with being CLOSER to both the crystal and Boltz-2's
        # independent prediction (median rho -0.14 for each). Its ARGMAX is not:
        # the maximum of a noisy score is an outlier, typically a strained pose
        # that happens to point the warhead well. Narrowing on anchoring and then
        # taking a TYPICAL member of what survives beats either alone.
        #
        # n = 15, so 33.3% against 26.7% is one molecule and the quartile width
        # is not tuned. What is not one molecule is 6.7% against 26.7%: argmax is
        # the thing to stop doing.
        reps = []
        for k in mode_ids:
            idx = np.flatnonzero(labels == k)
            a = anchor[idx]
            sub = dmat[np.ix_(idx, idx)]
            if np.all(np.isnan(a)) or len(idx) < 4:
                reps.append(int(idx[np.argmin(sub.mean(axis=1))]))
                continue
            keep = np.flatnonzero(a >= np.nanpercentile(a, 75))
            if len(keep) < 2:
                reps.append(int(idx[np.argmin(sub.mean(axis=1))]))
                continue
            s2 = sub[np.ix_(keep, keep)]
            reps.append(int(idx[keep[np.argmin(s2.mean(axis=1))]]))
        if not sdf.exists():                       # append_only: never overwrite
            write_sdf(mol, reps, sdf, modes=mode_ids)

        if do_gnina:
            tmp = work / "modes.sdf"
            write_sdf(mol, reps, tmp, modes=mode_ids)
            g = gnina_scores(plain_rec, tmp, gpu)
            if len(g) != len(reps):
                raise ValueError(f"gnina returned {len(g)} scores for {len(reps)} modes")
            supp = Chem.SDMolSupplier(str(tmp), removeHs=False)
            for kk, (m, i) in enumerate(zip(supp, reps)):
                if m is None:
                    raise ValueError(f"mode rep {kk} unreadable in the SDF handed to gnina")
                a = mol.GetConformer(i).GetPositions()
                b = m.GetConformer().GetPositions()
                if a.shape != b.shape or float(np.abs(a - b).max()) > _COORD_TOL:
                    raise ValueError(f"mode rep {kk}: SDF is not the conformer measured")
            for c in ("Affinity", "CNNscore", "CNNaffinity", "CNNvariance"):
                for kk in range(len(aggs)):
                    aggs[kk][c] = float(g[c].iloc[kk]) if c in g else np.nan
        return rows, aggs
    except Exception as exc:                       # noqa: BLE001
        # The traceback goes to the log, not just the status string. A one-line
        # message ("only integer scalar arrays...") is not enough to locate a
        # failure that hits every molecule, and finding that out after a
        # four-hour library run is the expensive way to learn it.
        log.debug("%s failed", cand.ident, exc_info=True)
        agg["status"] = f"failed: {str(exc)[:160]}"
        return pd.DataFrame(), [agg]
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--gpu", default="1")
    # 500, NOT 200 (@tt8804, 2026-08-07). 200 was inherited and never justified;
    # D0068 requires a number to carry its defining parameter. Measured on 15
    # crystal complexes docked at 2000 runs, each molecule's per-run hit rate
    # falls exponentially with size (rho = -0.683, p = 0.005; each extra heavy
    # atom multiplies it by 0.883). Covering 95% of the POOL means covering the
    # 5th-percentile-hardest molecule, which needs 227 runs for T_3 and 258 for
    # T_4 -- consistent with two independent 200-run dockings both returning
    # 93.3%, just under 95%. 500 is that estimate plus margin, since the fit is
    # on 15 molecules extrapolated past its size range and the library is docked
    # once. Costs ~25% more wall-clock than 200, not 2.5x: AutoDock-GPU runs the
    # LGA instances concurrently, so ~3.6 s of the per-molecule cost is fixed and
    # only ~0.0032 s is per-run.
    ap.add_argument("--nrun", type=int, default=500)
    ap.add_argument("--chunk", type=int, default=100)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--no-gnina", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=(logging.DEBUG if os.environ.get("NACV2_DEBUG") else logging.INFO),
                        format=f"%(levelname)s [s{args.shard}] %(message)s")
    os.nice(19)

    R.resolve_3ikd_ian()
    rec_dir = ns.build_reactive_receptor(ns.RX_RECEPTOR)
    plain_rec = sout.latest_path("blacksmith", "receptor_3ikd", "3IKD_prepared", ".pdbqt")
    log.info("reactive receptor ready; plain %s (3IKD_ian verified)", Path(plain_rec).name)

    cands = nr.load_candidates()
    cands = [c for i, c in enumerate(cands) if i % args.n_shards == args.shard]
    if args.limit:
        cands = cands[:args.limit]

    # Resume. A shard that dies at 90% must not restart from zero overnight.
    # Two rules, both learned the hard way:
    #   - only `ok` rows count as done. nac_rank.refine() counts `failed:` rows
    #     as done and therefore never retries a transient failure.
    #   - the run count is part of the identity. Skipping on ident alone would
    #     seat a 200-run measurement inside a re-run at a different effort,
    #     which is the same defect the BPMD resume had (catalogue #22).
    done = set()
    for f in sorted(glob.glob(str(OUT.dir / f"agg_s{args.shard}_*.csv"))):
        try:
            d = pd.read_csv(f)
        except Exception:                          # noqa: BLE001
            continue
        if {"ident", "status", "nrun"} <= set(d.columns):
            ok = d[(d.status == "ok") & (d.nrun == args.nrun)]
            # RESUME KEYS ON THE MOLECULE, NOT THE MODE. `ident` is now
            # `<parent>_m<k>`, so matching on it against a candidate's ident
            # would never hit and the resume would silently re-dock the entire
            # library -- four hours of GPU, reported as normal progress. Older
            # frames have no `parent_ident`, so fall back to `ident` for those.
            done |= set(ok["parent_ident"] if "parent_ident" in ok.columns
                        else ok["ident"])
    if done:
        before = len(cands)
        cands = [c for c in cands if c.ident not in done]
        log.info("resume: %d already complete at nrun=%d, %d remain of %d",
                 len(done), args.nrun, len(cands), before)
    log.info("%d candidates on this shard", len(cands))

    pose_buf, agg_buf, done = [], [], 0
    for i, c in enumerate(cands, 1):
        rows, aggs = one(c, rec_dir, Path(plain_rec), args.nrun, args.gpu,
                         not args.no_gnina)
        # `one` returns ONE ROW PER MODE now, so a molecule contributes several
        # candidates. `done` still counts MOLECULES -- resume and progress are
        # per-molecule, and counting modes would make the chunk size depend on
        # how multi-modal the shard happened to be.
        agg_buf.extend(aggs)
        if not rows.empty:
            pose_buf.append(rows)
        done += 1
        bad = [a for a in aggs if a["status"] != "ok"]
        if bad:
            log.warning("%s: %s", c.ident, bad[0]["status"])
        if done % args.chunk == 0 or i == len(cands):
            if pose_buf:
                pd.concat(pose_buf, ignore_index=True).to_csv(
                    OUT.write(f"poses_s{args.shard}", ".csv"), index=False)
            pd.DataFrame(agg_buf).to_csv(
                OUT.write(f"agg_s{args.shard}", ".csv"), index=False)
            ok = sum(a["status"] == "ok" for a in agg_buf)
            mols = len({a.get("parent_ident", a["ident"]) for a in agg_buf})
            log.info("[%d/%d] flushed; %d modes from %d molecules",
                     i, len(cands), ok, mols)
            pose_buf, agg_buf = [], []

    log.info("shard %d complete", args.shard)


if __name__ == "__main__":
    main()
