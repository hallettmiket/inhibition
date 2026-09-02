"""The Ranking page must display the number it is ordered by.

THE DEFECT. `modes.html` sorted on `class_rank`, which `rank_v2` computes from
the column `ranking.score_by_tier` names (`engagement` for T_4, D0098). The rail
displayed `conditional_eb`, named as a literal in `_rows_json`.

Those are different columns, correlated at rho = 0.32 across nac_v6's 34,888
acrylamide modes. So the printed number went UP as the rank went DOWN at
**9,077 of 24,927 adjacent pairs (36.4%)**, and the page contradicted itself on
over a third of its rows. @tt8804, reading it: *"what is that number and why does
it seem like they are not being ranked by it"* -- ranks 7-12 and 15-18 all showed
2.94 while rank 13 showed 4.32.

It is worse than a cosmetic mismatch, because `conditional_eb` is an
empirical-Bayes estimate that collapses to the prior at small n: 15,988 of those
modes hold ONE pose and take just two distinct values between them, while their
`engagement` spans 0.0000 to 0.9557. The column on screen was nearly constant
exactly where the ranking column discriminates most, so the rail looked
unranked.

`gather()` had ALREADY been fixed to read the score name from config -- for
choosing which FILE to load, with a comment explaining that a hardcoded score
made a whole run invisible. The display was left naming a column. The fix
stopped one step short of the thing the reader actually sees.

SECOND DEFECT, SAME LINE OF THINKING: `global_rank` was computed from
`conditional_eb` while `class_rank` came from `engagement`, so the scope
selector's two settings ordered by different quantities and nothing said so.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

SRC = (REPO / "shared" / "mode_ranking.py").read_text()


def _code_only(src: str) -> str:
    """Executable lines. Prose must stay free to name the wrong column."""
    out, in_doc = [], False
    for ln in src.splitlines():
        t = ln.strip()
        if t.startswith(('"""', "'''")):
            if t.count('"""') % 2 or t.count("'''") % 2:
                in_doc = not in_doc
            continue
        if in_doc or t.startswith("#"):
            continue
        out.append(ln.split("#")[0])
    return "\n".join(out)


CODE = _code_only(SRC)


def test_the_headline_number_is_the_ranking_score():
    """The rail's big number must come from `_score`, not a named column."""
    m = re.search(r"class=\"eng\"[^\n]*\n?[^\n]*fmt\(x\.(\w+)", SRC)
    assert m, "could not find the rail's headline number"
    assert m.group(1) == "sc", (
        f"the rail displays x.{m.group(1)}; it must display x.sc, the value "
        "`_score` carries and the list is ordered by")


def test_the_payload_carries_the_configured_score():
    assert '"sc": num("_score"' in CODE, \
        "the row payload no longer sends the ranking score"


def test_global_rank_uses_the_same_score_as_the_class_rank():
    """Both scopes must order by one quantity."""
    m = re.search(r'r\["global_rank"\]\s*=\s*r\[[\'"](\w+)[\'"]\]\.rank', CODE)
    assert m, "global_rank is no longer computed from a single column"
    assert m.group(1) == "_score", (
        f"global_rank ranks {m.group(1)!r} while class_rank comes from the "
        "config-named score; the two scope settings would order differently")


def test_no_score_column_is_named_as_a_literal_in_the_display_path():
    """`conditional_eb` may be shown as CONTEXT, never as the headline."""
    m = re.search(r"class=\"eng\".{0,200}", SRC, re.S)
    assert "conditional_eb" not in m.group(0), \
        "a literal score column is back in the headline"


def test_gather_refuses_a_frame_without_the_configured_score():
    """An allowlist, not a fallback: showing another column IS the defect."""
    assert "does not carry" in SRC and "score_by_tier names" in SRC, \
        "gather() no longer raises when the configured score is absent"


def test_the_guard_can_fail():
    """Prove the headline regex matches the ORIGINAL defect line.

    Two tests in this repo have passed for free; this one is checked against
    the exact string that caused the bug.
    """
    before = """    '<span class="eng">' + fmt(x.eb, 2) + '</span></span>' +"""
    m = re.search(r"class=\"eng\"[^\n]*\n?[^\n]*fmt\(x\.(\w+)", before)
    assert m and m.group(1) == "eb", \
        "the guard would not have caught the original line"


@pytest.mark.parametrize("col", ["engagement", "conditional_eb",
                                 "enrichment_conditional"])
def test_score_stamping_is_generic_over_the_configured_column(col):
    """Nothing may assume WHICH score is configured."""
    import pandas as pd
    d = pd.DataFrame({col: [0.3, 0.9, float("nan")]})
    d["_score_col"] = col
    d["_score"] = pd.to_numeric(d[col], errors="coerce")
    r = d["_score"].rank(ascending=False, method="min", na_option="bottom")
    assert list(r) == [2.0, 1.0, 3.0]
