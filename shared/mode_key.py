"""The one way to key a binding mode, and to join two tables of modes.

WHY THIS EXISTS (#53). A mode was written two different ways:

    rank_v2      t4_716800c125a7_m0      always `_m<mode>`
    attack_sweep t4_716800c125a7         BARE for mode 0, `_m<mode>` otherwise

So the obvious join -- `sweep.merge(rank, on="ident")` -- silently drops every
mode-0 row, which is every row that was actually simulated. Nothing errors. The
merge returns a smaller frame that still looks like a frame.

That collision cost two wrong answers in one day: it hid #53 (the sweep took mode
0 for 242 of 242 molecules while the ranking is per mode) and it produced a wrong
sweep-vs-MD correlation in #36 that had to be retracted.

THE KEY IS `(parent_ident, mode)`, NEVER `ident`. `ident` is a display label. Two
tables agreeing about a mode is a fact about the pair, and this module is the
only place that pair is constructed.
"""

from __future__ import annotations

import re

import logging
import pandas as pd

log = logging.getLogger(__name__)

#: `t4_x_m3` -> ("t4_x", 3). Anchored at the end so a molecule whose own name
#: contains `_m` followed by digits is not truncated.
_SUFFIX = re.compile(r"^(?P<parent>.+)_m(?P<mode>\d+)$")


def split_ident(ident: str) -> tuple[str, int | None]:
    """('t4_x_m3') -> ('t4_x', 3); ('t4_x') -> ('t4_x', None).

    None, not 0. A bare ident means *the mode was not stated*, which is a
    different claim from "mode 0" — and reading it as 0 is precisely the
    assumption that produced #53's invisible collision. Callers that know the
    table's convention can substitute 0 explicitly and be seen doing it.
    """
    m = _SUFFIX.match(str(ident))
    if not m:
        return str(ident), None
    return m.group("parent"), int(m.group("mode"))


def add_key(df: pd.DataFrame, bare_is_mode_zero: bool = False) -> pd.DataFrame:
    """Add `parent_ident`, `mode` and `mode_key` columns, without guessing.

    Uses the frame's own `parent_ident`/`mode` columns when it has them, and
    falls back to parsing `ident` only for what is missing.

    `bare_is_mode_zero` is the ONE place the historical convention is applied,
    and it must be passed explicitly. Pass it for `attack_sweep` tables written
    before #53, where a bare ident did mean mode 0. Do not pass it for anything
    else: in `rank_v2` a bare ident is a molecule-level row, not mode 0.
    """
    d = df.copy()
    parsed = [split_ident(i) for i in d.get("ident", pd.Series(dtype=str))]
    if "parent_ident" not in d.columns:
        d["parent_ident"] = [p for p, _ in parsed]
    else:
        d["parent_ident"] = d["parent_ident"].fillna(
            pd.Series([p for p, _ in parsed], index=d.index))
    if "mode" not in d.columns:
        d["mode"] = [m for _, m in parsed]
    else:
        d["mode"] = d["mode"].where(d["mode"].notna(),
                                    pd.Series([m for _, m in parsed], index=d.index))
    if bare_is_mode_zero:
        d["mode"] = d["mode"].fillna(0)
    # A ROW THAT CANNOT BE KEYED IS DROPPED FROM THE KEY, NOT RAISED ON.
    #
    # `int(r["mode"])` raised ValueError on a single corrupt row whose columns
    # had shifted (an `elevate_why` string sitting in `pose_rank`, booleans in
    # `mode`), and that one row took down `build_gui` AND `sweep_combine` --
    # so the whole GUI froze at its last good build while 1,544 perfectly good
    # results sat on disk. One unparseable record must not be able to do that.
    #
    # It is a WARNING and a null key, not a silent zero: a row with no mode
    # cannot be joined to anything, and giving it mode 0 would attach it to a
    # real mode of the same molecule -- which is the failure this module exists
    # to prevent.
    def _key(r):
        m = r["mode"]
        if pd.isna(m):
            return None
        try:
            return f"{r.parent_ident}|{int(m)}"
        except (TypeError, ValueError):
            return None

    d["mode_key"] = d.apply(_key, axis=1)
    unkeyed = int(d["mode_key"].isna().sum() - d["mode"].isna().sum())
    if unkeyed > 0:
        log.warning("%d row(s) have a mode that is not an integer and cannot be "
                    "keyed; they are excluded from mode-level joins", unkeyed)
    return d


def join(left: pd.DataFrame, right: pd.DataFrame, how: str = "left",
         left_bare_is_mode_zero: bool = False,
         right_bare_is_mode_zero: bool = False,
         suffixes: tuple[str, str] = ("", "_r")) -> pd.DataFrame:
    """Merge two mode-level tables on (parent_ident, mode).

    Never on `ident`. That is the whole point.
    """
    l = add_key(left, left_bare_is_mode_zero)
    r = add_key(right, right_bare_is_mode_zero)
    r = r.drop(columns=[c for c in ("parent_ident", "mode") if c in r.columns])
    return l.merge(r, on="mode_key", how=how, suffixes=suffixes)
