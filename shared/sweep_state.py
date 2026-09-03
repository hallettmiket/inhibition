"""What the sweep has done, is doing, and has not started (#63).

@tt8804: "show results as pending what is sweeped and show sweep results as it
goes". A sweep of 170 modes takes two days, so a page that shows only finished
rows is blank for the first half hour and silent about whether anything is
happening at all.

FOUR STATES, AND THE DISTINCTIONS ARE LOAD-BEARING:

    ok       a triage trajectory finished and was measured
    failed   it was attempted and did not produce a result, WITH the reason
    pending  on the active worklist, no result row yet -- queued or in flight
    (absent) ranked, never selected for this campaign

`pending` is the one that needs care. It is defined as "on the worklist and not
in the results", NOT as "a process is running" -- a page cannot see the process
table, and a worker that died would otherwise leave a row claiming to be in
flight forever. A mode that is queued and one that is mid-trajectory are both
honestly "not finished yet", and the page says exactly that.

THE WORKLIST IS NAMED, NOT INFERRED. Two sessions wrote different worklists on
2026-08-12 -- one scoped to the campaign's three warhead classes, one spanning
all nine -- and they overlap in a single mode. "Newest file wins" would make this
page report on whichever list was written last, which may not be the list the
running workers are executing. The caller passes the path and the page prints it.
"""

from __future__ import annotations

import glob
import os
from pathlib import Path

import pandas as pd

from . import run_paths as rp

B = rp.BLACKSMITH


def results() -> pd.DataFrame:
    """Every sweep row THIS RUN has written, newest wins per mode.

    "Ever written" is what it used to mean, and that was the defect: the sweep
    tables were a flat directory shared by every screen, so a freshly bumped
    topic still listed 554 rows from three superseded runs.
    """
    fs = [str(f) for f in rp.sweep_result_files()]
    if not fs:
        return pd.DataFrame()
    out = []
    for f in fs:
        try:
            d = pd.read_csv(f)
        except Exception:                                  # noqa: BLE001
            continue
        d["_t"] = os.path.getmtime(f)
        out.append(d)
    if not out:
        return pd.DataFrame()
    d = pd.concat(out, ignore_index=True)
    # A LATER ATTEMPT SUPERSEDES AN EARLIER ONE. A mode that failed, was
    # unblocked and re-run must read as ok, not carry its old failure.
    d = d.sort_values("_t").drop_duplicates("ident", keep="last")

    # BACKFILL THE COMPARABLE COLUMN FOR PRE-ADAPTIVE ROWS.
    #
    # `frac_attack_ready_common` (engagement over the first 1.2 ns, whatever the
    # run length) exists so variable-length adaptive rows can be ranked against
    # the fixed-length ones. Rows written before it existed do not carry it --
    # but for a run that IS exactly one common window long the two figures are
    # the same number by construction, so it can be filled rather than left
    # null and treated as missing data.
    #
    # ONLY where the length matches. A row from some other fixed length is left
    # NaN: inventing a common-window value for a window that was never observed
    # is exactly the kind of plausible-but-unmeasured number this project keeps
    # being bitten by.
    if "frac_attack_ready" in d.columns:
        try:
            import attack_sweep as _asw
            win = float(_asw.COMMON_WINDOW_PS)
        except Exception:                                 # noqa: BLE001
            win = 1200.0
        if "frac_attack_ready_common" not in d.columns:
            d["frac_attack_ready_common"] = float("nan")
        exact = d.get("sweep_ps", pd.Series(index=d.index, dtype=float)) == win
        fill = d.frac_attack_ready_common.isna() & exact
        d.loc[fill, "frac_attack_ready_common"] = d.loc[fill, "frac_attack_ready"]
    return d


def state(worklist: Path | None = None) -> pd.DataFrame:
    """One row per mode with a `sweep_state` column, joined to the ranking.

    Rows are the union of (everything on the worklist) and (everything with a
    result), so a mode swept under an earlier campaign still appears -- it is a
    measurement of this library and hiding it would misreport what is known.
    """
    from shared import mode_ranking as mr
    rk = mr.gather()
    res = results()
    wl = pd.DataFrame()
    if worklist and Path(worklist).is_file():
        wl = pd.read_csv(worklist)
        # `ident` IN A WORKLIST IS THE MOLECULE; `task_id` IS THE MODE, and every
        # result is keyed on the mode. Merging on `ident` therefore matched
        # NOTHING: 4,295 worklist rows and 261 result rows produced 4,556 rows
        # with zero overlap, so the whole worklist read as `pending` however many
        # sweeps had finished, and not one finished mode carried `_queued` --
        # which is the flag `sweep_combine` uses to mean "this campaign".
        # Same defect as catalogue #23, one file further along; `sweep_assets`,
        # `sweep_combine` and `recompute_attack_ready` already key on `task_id`.
        if "task_id" in wl.columns:
            wl = wl.assign(ident=wl["task_id"].astype(str))
        wl["_queued"] = True

    # A ROW FOR A MODE NOBODY SELECTED IS NOT A RESULT OF THIS CAMPAIGN.
    #
    # 24 rows in this topic come from a launcher that read the worklist
    # positionally and asked for `global_rank` as `pose_rank` -- ranks in the
    # hundreds against molecules with a handful of poses. The runner refused
    # each one against the pose SDF, which is the guard working; but the refusal
    # rows are keyed on (molecule, a pose_rank no real run will ever request),
    # so nothing supersedes them and Home read `failed: 24` permanently.
    #
    # They cannot be deleted -- the outputs root is append-only. So they are not
    # COUNTED: a row whose (parent, pose_rank) is not on the worklist describes
    # a job this campaign never asked for. Requires a worklist to judge against,
    # so with none passed nothing is dropped.
    if not res.empty and not wl.empty and \
            {"parent_ident", "pose_rank"} <= set(res.columns) and \
            {"parent_ident", "pose_rank"} <= set(wl.columns):
        asked = set(zip(wl.parent_ident.astype(str), wl.pose_rank.astype(int)))
        keep_row = [(p, r) in asked for p, r in
                    zip(res.parent_ident.astype(str), res.pose_rank.astype(int))]
        dropped = len(res) - sum(keep_row)
        if dropped:
            import logging as _lg
            _lg.getLogger(__name__).info(
                "%d sweep row(s) name a (molecule, pose_rank) this campaign "
                "never selected — not counted", dropped)
        res = res[keep_row]

    idents = set()
    if not res.empty:
        idents |= set(res.ident.astype(str))
    if not wl.empty:
        idents |= set(wl.ident.astype(str))
    if not idents:
        return pd.DataFrame()

    keep = [c for c in ("ident", "parent_ident", "warhead_class", "mode_label",
                        "enrichment", "viable_fraction", "conditional_eb",
                        "class_rank", "n_poses_mode") if c in rk.columns]
    base = (rk[rk.ident.astype(str).isin(idents)][keep].copy()
            if not rk.empty else pd.DataFrame({"ident": sorted(idents)}))
    # A mode on the worklist but absent from the ranking still gets a row: that
    # combination means the two disagree, which is worth SEEING rather than
    # dropping silently (it is how the T_3 rows stayed on a list after the tier
    # decision).
    missing = idents - set(base.ident.astype(str))
    if missing:
        base = pd.concat([base, pd.DataFrame({"ident": sorted(missing)})],
                         ignore_index=True)

    # DERIVED FROM THE RESULTS TABLE, NOT LISTED BY HAND. This was an allowlist
    # of nine column names, so every measurement added to a sweep row after it
    # was written -- `rmsd_max_a`, `rmsd_mean_a`, `elevate`, `pose_held`,
    # `attack_ready_max_a` -- was silently dropped on the way to the page, and
    # the page then RECOMPUTED one of them by a second route (catalogue #5, and
    # the reason `sweep_combine` was reading RMSD off the wrong trajectory).
    # Carry everything the results table has; drop only what `base` already
    # provides, so the merge cannot produce `_x`/`_y` pairs.
    _own = set(base.columns) - {"ident"}
    rescols = ["ident"] + [c for c in res.columns
                           if c != "ident" and c not in _own]
    d = base.merge(res[rescols], on="ident", how="left") if not res.empty else base
    if not wl.empty:
        d = d.merge(wl[["ident", "_queued"]], on="ident", how="left")
    if "_queued" not in d.columns:
        d["_queued"] = False
    d["_queued"] = d["_queued"].fillna(False).astype(bool)

    def _st(r):
        # `pd.isna` EXPLICITLY, NOT `or ""`. A missing status arrives as NaN, and
        # NaN is TRUTHY -- `nan or ""` evaluates to nan, `str(nan)` is "nan", and
        # every mode that had simply not been swept yet was reported as FAILED.
        # 162 of them, on the first run of this page.
        v = r.get("status")
        s = "" if v is None or (isinstance(v, float) and pd.isna(v)) else str(v).strip()
        if s == "ok":
            return "ok"
        # A GIVE-UP IS A RESULT, NOT A FAILURE. The pose was measured after
        # equilibration and had left; that is the campaign's own triage working,
        # and counting it as a failure makes a deliberate saving look like a
        # crash. It is reported separately.
        if s.startswith("aborted"):
            return "left early"
        # A PLACEHOLDER IS NOT AN ATTEMPT. `stage0 only` rows carry no
        # measurement -- they are the free geometry probe -- and `invalidated`
        # rows were withdrawn on purpose. Both were counted as `failed`, which
        # is how Home came to read "failed: 16" on a run with zero failures.
        if s.startswith("stage0") or s.startswith("invalidated"):
            return "not sent"
        if s:                       # any recorded non-ok status is an attempt
            return "failed"
        return "pending" if r["_queued"] else "not sent"

    d["sweep_state"] = d.apply(_st, axis=1)
    return d


def predicts(d: pd.DataFrame) -> dict:
    """Does the docked ranking predict what the trajectory did?

    THE SWEEP'S OWN HEADLINE, so it belongs on the sweep's page rather than in a
    message. Spearman of `enrichment` against `frac_attack_ready` over the modes
    that came back, plus the split between productive and not.

    IT IS RANGE-RESTRICTED AND THE PAGE MUST SAY SO. Every swept mode cleared the
    enrichment floor, so this measures whether enrichment discriminates ABOVE the
    floor -- not whether the floor works. Those are different claims, and the
    second one needs modes sampled from BELOW the floor, which is exactly what
    `sweep_rule.pilot` is for and has not been run. Reporting r ~ 0 as "the
    criterion is uninformative" would be the range-restriction mistake in print.
    """
    out = {"n": 0}
    if d.empty or "frac_attack_ready" not in d.columns or "enrichment" not in d.columns:
        return out
    # THIS CAMPAIGN ONLY. `state()` also carries every mode ever swept, and those
    # come from other worklists, other warhead classes and (before the tier
    # decision) other generation tiers. Pooling them lowers the apparent range
    # restriction by mixing in modes that were selected under different rules --
    # which changes the correlation without making it mean more. Measured on
    # 2026-08-13: pooled gave rho = +0.17 over 128, the campaign alone +0.07 over
    # 101, and only the second is a statement about the rule now in force.
    m = d[(d.sweep_state == "ok") & d.enrichment.notna()
          & d.frac_attack_ready.notna()]
    if "_queued" in d.columns and bool(d["_queued"].any()):
        m = m[m["_queued"]]
    if len(m) < 8:
        return {"n": int(len(m))}
    try:
        from scipy import stats
        r, p = stats.spearmanr(m.enrichment, m.frac_attack_ready)
    except Exception:                                      # noqa: BLE001
        return {"n": int(len(m))}
    prod = m.frac_attack_ready > 0.01
    return {"n": int(len(m)), "rho": float(r), "p": float(p),
            "productive": int(prod.sum()),
            "enr_prod": float(m[prod].enrichment.median()) if prod.any() else None,
            "enr_not": float(m[~prod].enrichment.median()) if (~prod).any() else None,
            "floor": float(m.enrichment.min())}


def summary(d: pd.DataFrame) -> dict:
    """Counts by state, for the step nav and the page header."""
    if d.empty:
        return {"ok": 0, "failed": 0, "pending": 0, "not sent": 0,
                "left early": 0, "productive": 0}
    c = d.sweep_state.value_counts().to_dict()
    out = {k: int(c.get(k, 0))
           for k in ("ok", "failed", "pending", "not sent", "left early")}
    ok = d[d.sweep_state == "ok"]
    out["productive"] = (int((ok.frac_attack_ready > 0.01).sum())
                         if "frac_attack_ready" in ok.columns else 0)
    return out
