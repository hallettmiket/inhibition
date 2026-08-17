"""What the sweep has done, is doing, and has not started (#63).

@tt8804: "show results as pending what is sweeped and show sweep results as it
goes". A sweep of 170 modes takes two days, so a page that shows only finished
rows is blank for the first half hour and silent about whether anything is
happening at all.

FOUR STATES, AND THE DISTINCTIONS ARE LOAD-BEARING:

    ok       a 10 ns trajectory finished and was measured
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

B = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith")


def results() -> pd.DataFrame:
    """Every sweep row THIS RUN has written, newest wins per mode.

    "Ever written" is what it used to mean, and that was the defect: the sweep
    tables were a flat directory shared by every screen, so a freshly bumped
    topic still listed 554 rows from three superseded runs.
    """
    fs = sorted(glob.glob(str(rp.sweep_dir() / "attack_sweep_*.csv")),
                key=os.path.getmtime)
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
        wl["_queued"] = True

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

    rescols = [c for c in ("ident", "status", "frac_attack_ready", "n_visits",
                           "frac_in_window", "median_dist_a", "median_angle_deg",
                           "min_dist_a", "sweep_ps", "_t") if c in res.columns]
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
        return {"ok": 0, "failed": 0, "pending": 0, "not sent": 0, "productive": 0}
    c = d.sweep_state.value_counts().to_dict()
    out = {k: int(c.get(k, 0)) for k in ("ok", "failed", "pending", "not sent")}
    ok = d[d.sweep_state == "ok"]
    out["productive"] = (int((ok.frac_attack_ready > 0.01).sum())
                         if "frac_attack_ready" in ok.columns else 0)
    return out
