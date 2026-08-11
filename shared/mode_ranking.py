"""Every binding mode, ranked individually, with what was and was not simulated.

WHY (#53). The pipeline scores modes independently and ranks them in one pooled
list -- "a mode IS a candidate row". The GUI's rail then shows one row per
MOLECULE, because it indexes the sweep by `parent_ident`. So the per-mode ranking
is computed and then never seen, and the fact that the sweep took mode 0 for 242
of 242 molecules -- including four cases where a DIFFERENT mode ranked first in
its warhead class -- was invisible in every view the project had.

This page is the missing view. Every mode, ordered by its own rank within its
warhead class, each stamped with whether it was swept and whether it reached
100 ns. What it is for is reading the gap: high-ranked modes with nothing behind
them are the ones the screen may have missed.

RANK IS WITHIN A WARHEAD CLASS, NEVER ACROSS. The SN2 angular criterion is far
stricter than the perpendicular one (#47), so a cross-class ordering compares
scores computed under different bars. The page offers no global rank column.

THE JOIN IS ON (parent_ident, mode). Never on `ident`: mode 0 is the bare ident
in the sweep table and `_m0` in the rank table, so a merge on the label drops
exactly the rows that were simulated (`shared/mode_key.py`).
"""

from __future__ import annotations

import glob
import html
from pathlib import Path

import pandas as pd

from shared import mode_key as mk

B = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith")


def _latest(pattern: str) -> Path | None:
    fs = sorted(glob.glob(str(B / pattern)))
    return Path(fs[-1]) if fs else None


def gather() -> pd.DataFrame:
    """One row per mode: its rank, its docking-derived scores, its simulations."""
    frames = []
    for tier, score in (("T4", "conditional_eb"), ("T3", "enrichment_conditional")):
        f = _latest(f"rank_v2/rank_v2_{tier}_{score}_*.csv")
        if f is None:
            continue
        d = pd.read_csv(f)
        d = d[d.get("mode").notna()] if "mode" in d.columns else d
        d["tier"] = d.get("tier", tier)
        frames.append(d)
    if not frames:
        return pd.DataFrame()
    r = pd.concat(frames, ignore_index=True)

    sf = sorted(glob.glob(str(B / "attack_sweep/attack_sweep_*.csv")))
    sweep = (pd.concat([pd.read_csv(x) for x in sf], ignore_index=True)
             .drop_duplicates("ident", keep="last") if sf else pd.DataFrame())

    if not sweep.empty:
        keep = [c for c in ("ident", "parent_ident", "mode", "frac_attack_ready",
                            "n_visits", "status") if c in sweep.columns]
        # ATTEMPTED is not the same as SUCCEEDED. A sweep row exists for every
        # mode that was sent; `frac_attack_ready` is null when the run failed.
        # Counting only the successful ones as "swept" would report a mode that
        # was tried and crashed as one nobody ever looked at.
        sweep = sweep.assign(_sent=True)
        keep = keep + ["_sent"]
        # bare_is_mode_zero: these rows were written before #53, when the sweep
        # wrote the bare ident for mode 0. Stated explicitly, not assumed.
        r = mk.join(r, sweep[keep].rename(columns={"status": "sweep_status"}),
                    right_bare_is_mode_zero=True, suffixes=("", "_sw"))
    else:
        r["frac_attack_ready"] = None

    # 100 ns runs are keyed by MOLECULE -- the rows do not record which mode was
    # simulated (fixed going forward, not retrofittable). So a mode is marked
    # "ran 100 ns" only when its molecule ran AND that molecule's swept mode is
    # this one. Anything else would put an MD badge on a mode that never moved.
    md_ids: set[str] = set()
    for f in glob.glob(str(B / "md_residence/*.csv")):
        try:
            d = pd.read_csv(f)
        except Exception:                                  # noqa: BLE001
            continue
        if "ident" not in d.columns or "production_ps" not in d.columns:
            continue
        d = d[(d.production_ps >= 50000)
              & d.status.astype(str).str.startswith("ok")]
        md_ids |= set(d.ident.astype(str))

    r["sent"] = r["_sent"].fillna(False).astype(bool) if "_sent" in r else False
    r["swept"] = (r["frac_attack_ready"].notna()
                  if "frac_attack_ready" in r else False)
    r["ran_md"] = r.parent_ident.isin(md_ids) & r["_sent"].fillna(False).astype(bool)
    return r


def build(title: str, date_str: str) -> str:
    r = gather()
    if r.empty:
        return "<!doctype html><p>no rank tables found</p>"

    n_modes = len(r)
    n_sent = int(r.sent.sum())
    n_swept = int(r.swept.sum())
    n_md = int(r.ran_md.sum())
    n_mol = r.parent_ident.nunique()
    swept_mode0 = int(((r["mode"] == 0) & r.sent).sum())

    # The finding, computed rather than asserted: modes that lead their class and
    # were never simulated.
    # Never SENT, not merely never scored — a mode that was tried and failed is
    # a different fact from one that was never chosen.
    miss = r[(r.class_rank == 1) & (~r.sent)] if "class_rank" in r.columns \
        else r.iloc[0:0]

    rows = []
    for _, x in r.sort_values(
            ["warhead_class", "class_rank"], na_position="last").iterrows():
        if pd.isna(x.get("class_rank")):
            continue
        if x.ran_md:
            badge, cls = "100 ns", "ok"
        elif x.swept:
            badge, cls = "swept", "mid"
        elif x.sent:
            badge, cls = "sweep failed", "mid"
        else:
            badge, cls = "not simulated", "no"
        rows.append(
            f"<tr class='{cls}'><td>{int(x.class_rank)}</td>"
            f"<td class='id'>{html.escape(str(x.ident))}</td>"
            f"<td>{html.escape(str(x.warhead_class))}</td>"
            f"<td>m{int(x['mode'])}</td>"
            f"<td>{int(x.n_poses_mode) if pd.notna(x.get('n_poses_mode')) else '—'}</td>"
            f"<td>{x.viable_fraction*100:.1f}%</td>"
            f"<td>{x.get('conditional_eb', float('nan')):.3f}</td>"
            f"<td><span class='b {cls}'>{badge}</span></td></tr>")

    missrows = "".join(
        f"<li><code>{html.escape(str(x.ident))}</code> &mdash; ranked "
        f"<strong>1st</strong> in <code>{html.escape(str(x.warhead_class))}</code>, "
        f"never simulated</li>" for _, x in miss.iterrows())

    return f"""<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(title)} — every mode</title><style>
:root{{--ink:#12181f;--muted:#5b6b80;--rule:#dfe4ea;--paper:#fff;--raise:#f6f8fa;
 --ok:#0f7a54;--mid:#8a6d1f;--no:#b3261e;--blue:#0b5fa5}}
*{{box-sizing:border-box}}
body{{margin:0;background:var(--paper);color:var(--ink);
 font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif}}
.wrap{{max-width:1220px;margin:0 auto;padding:28px 30px 80px}}
h1{{font-size:1.5rem;margin:0 0 .2rem}}
p.sub{{color:var(--muted);margin:0 0 1.4rem}}
.stats{{display:flex;flex-wrap:wrap;gap:1.6rem;padding:1rem 0;
 border-top:2px solid var(--rule);border-bottom:2px solid var(--rule);margin-bottom:1.4rem}}
.stat b{{display:block;font-size:1.5rem;font-family:ui-monospace,monospace}}
.stat span{{color:var(--muted);font-size:12px}}
.callout{{border-left:3px solid var(--no);background:#fdf5f4;padding:.8rem 1rem;
 margin:1.2rem 0;border-radius:0 4px 4px 0}}
.callout ul{{margin:.5rem 0 0;padding-left:1.2rem}}
table{{border-collapse:collapse;width:100%;font-size:13px;margin-top:1rem}}
th,td{{padding:.36rem .7rem;border-bottom:1px solid var(--rule);text-align:right;
 white-space:nowrap}}
th{{position:sticky;top:0;background:var(--paper);font:600 11px/1.4 inherit;
 color:var(--muted);text-transform:uppercase;letter-spacing:.05em;
 border-bottom:2px solid var(--rule)}}
td.id,th:nth-child(2),th:nth-child(3),td:nth-child(3),th:last-child,td:last-child
 {{text-align:left}}
td.id{{font-family:ui-monospace,monospace}}
tr.no td{{background:#fffaf9}}
.b{{font:600 11px inherit;padding:.1rem .45rem;border-radius:3px}}
.b.ok{{background:#e6f4ee;color:var(--ok)}}
.b.mid{{background:#faf3e0;color:var(--mid)}}
.b.no{{background:#fdeceb;color:var(--no)}}
a{{color:var(--blue)}}
p.nav{{margin:0 0 .8rem;font-size:13px}}
p.nav span{{color:var(--muted);margin:0 .4rem}}
</style></head><body><div class="wrap">
<p class="nav"><a href="combined.html">&#8592; sweep &amp; MD results</a>
 <span>&middot;</span> <a href="pipeline.html">how this works</a></p>
<h1>{html.escape(title)} &mdash; every mode, ranked individually</h1>
<p class="sub">{html.escape(date_str)} &middot; <strong>the ranking view</strong>
&mdash; every molecule and every mode the screen scored, simulated or not. The
<a href="combined.html">other view</a> carries the sweep and MD results, for the
subset that was simulated. Rank is <strong>within a warhead
class</strong>, never across &mdash; the S<sub>N</sub>2 angular criterion is
stricter than the perpendicular one (<a
href="https://github.com/hallettmiket/inhibition/issues/47">#47</a>), so a global
order would compare scores computed under different bars.</p>

<div class="stats">
<div class="stat"><b>{n_modes:,}</b><span>modes ranked</span></div>
<div class="stat"><b>{n_mol:,}</b><span>molecules</span></div>
<div class="stat"><b>{n_sent}</b><span>sent to the 10 ns sweep</span></div>
<div class="stat"><b>{n_swept}</b><span>of those returned a score</span></div>
<div class="stat"><b>{n_md}</b><span>reached 100 ns</span></div>
<div class="stat"><b>{swept_mode0}</b><span>of the swept are mode 0</span></div>
</div>

<p>The pipeline scores modes independently and ranks them in one pooled list
&mdash; <em>a mode is a candidate row</em>. Selection for simulation did not
follow that: it took <strong>mode 0, once per molecule</strong>. This page is the
view that was missing, so the gap can be read directly. Rows in red were never
simulated at any stage.</p>

{f'''<div class="callout"><strong>Modes that lead their warhead class and were
never simulated.</strong><ul>{missrows}</ul></div>''' if missrows else ''}

<table><thead><tr><th>class rank</th><th>mode</th><th>warhead class</th>
<th>m</th><th>poses</th><th>viable</th><th>conditional_eb</th><th>status</th>
</tr></thead><tbody>{''.join(rows)}</tbody></table>
</div></body></html>"""
