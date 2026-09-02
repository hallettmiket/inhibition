"""Resolving `<stem>_N` by a LEXICAL sort is a stale pin wearing a glob's clothes.

THE DEFECT THIS PINS. `data/reference/` is integer-versioned per the lab data
rule, and `reference_set.latest_reference` exists so no one hand-pins a version.
Two scripts still rolled their own:

    fs = sorted(glob.glob(".../warhead_classes_*.csv"))
    d  = pd.read_csv(fs[-1])

which reads like dynamic resolution and is not. `sorted` on strings puts
`warhead_classes_9.csv` AFTER `warhead_classes_10.csv`, so the moment the
library reached two digits the newest file became unreachable. `_10` adds
exactly one class -- `cinnamamide` -- so `dock_reference_modes` could not run an
aryl Michael acceptor at all, and `shortlist_report`, whose lookup returns None
on a miss, marked no reactive atom for one and said nothing.

`tests/test_reference_version.py` could not catch it: that test hunts for pinned
version LITERALS, and there is no literal here. The bug is in the comparison.

WHY A SOURCE SCAN RATHER THAN A BEHAVIOUR TEST. The wrong ordering only bites
once a stem passes nine versions, so a behaviour test over today's files passes
for every stem still in single digits -- it would have gone green on
`pin1_reference_binders` (4 versions) while the warhead library was already
broken. Banning the construct fails now, at nine versions or at one.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

SCAN_DIRS = ("scripts", "shared", "integration", "approaches", "tests")

#: A glob whose pattern ends in the integer-version wildcard.
_VERSIONED_GLOB = re.compile(r"""glob\s*\(\s*[^)]*_\*\.(csv|parquet)""")

#: ...and the result being INDEXED to one element. That narrowing is the whole
#: point of the rule. `sorted(glob("agg_s*_*.csv"))` gathers every shard of a
#: run and its order does not encode a version, so flagging it would be noise --
#: and a guard that cries wolf on twenty correct lines gets deleted, which is
#: how this project loses guards. Selecting `[-1]` from a sorted versioned glob
#: is the ONLY construct that means "the newest one", and it is the one that is
#: wrong.
_PICKS_ONE = re.compile(r"\[\s*-?\s*(?:1|0)\s*\]")


def _code_lines(path: Path):
    """Executable lines only -- prose must stay free to describe the defect.

    The version-pin test in this repo once flagged four "offenders" that were
    all docstrings. Comments and docstrings are stripped here for that reason.
    """
    out, in_doc = [], False
    for n, ln in enumerate(path.read_text(errors="ignore").splitlines(), 1):
        s = ln.strip()
        if s.startswith(('"""', "'''")):
            if s.count('"""') % 2 or s.count("'''") % 2:
                in_doc = not in_doc
            continue
        if in_doc or s.startswith("#"):
            continue
        out.append((n, ln.split("#")[0]))
    return out


def _files():
    for d in SCAN_DIRS:
        p = REPO / d
        if p.is_dir():
            yield from sorted(p.rglob("*.py"))


def test_no_lexical_sort_over_an_integer_versioned_glob():
    """`sorted(glob("..._*.csv"))` must not decide which version is newest."""
    offenders = []
    for f in _files():
        if f.name == Path(__file__).name:
            continue
        lines = _code_lines(f)
        for i, (n, ln) in enumerate(lines):
            if not _VERSIONED_GLOB.search(ln) or "sorted" not in ln:
                continue
            # THE WHOLE STATEMENT, not the one line the glob sits on. These
            # calls wrap, and `key=` -- the thing that makes an ordering
            # correct -- is usually on the continuation line. Judging line by
            # line reported eleven correct call sites as defects, which is how
            # a guard gets switched off.
            stmt, depth = "", 0
            for m, nxt in lines[i:i + 8]:
                stmt += nxt + "\n"
                depth += nxt.count("(") - nxt.count(")")
                if depth <= 0 and stmt.strip():
                    break
            if "key=" in stmt:
                continue
            # ...and only if ONE element is taken out of it. `sorted(glob(
            # "agg_s*_*.csv"))` gathers every shard of a run; its order encodes
            # no version and it is concatenated, not indexed. Selecting `[-1]`
            # from a sorted versioned glob is the only construct that means
            # "the newest one", and it is the one that is wrong.
            var = None
            head = ln.split("=")[0]
            if "=" in ln and "glob" not in head and "(" not in head:
                var = head.strip()
            window = stmt
            for m, nxt in lines[i + 1:i + 9]:
                if var is None or re.search(rf"\b{re.escape(var)}\s*\[", nxt):
                    window += nxt + "\n"
            if _PICKS_ONE.search(window):
                offenders.append(f"{f.relative_to(REPO)}:{n}: {ln.strip()}")
    assert not offenders, (
        "these order integer-versioned files LEXICALLY, so `_9` beats `_10`. "
        "Use shared.reference_set.latest_reference() (or sort on the parsed "
        "int):\n  " + "\n  ".join(offenders))


def test_the_resolver_disagrees_with_a_lexical_sort_and_is_right():
    """The resolver must compare integers, not strings.

    BOTH ORDERINGS ARE COMPUTED HERE, and no version literal is written down.
    Naming the expected file would make this test the very thing it polices --
    `tests/test_reference_version.py` walks the AST for exactly that and
    flagged the first draft of this test for hardcoding the older version.
    """
    from shared import reference_set as rs
    names = [p.name for p in (REPO / "data/reference").glob("warhead_classes_*.csv")]
    assert names, "no warhead library on disk"
    lexical = sorted(names)[-1]
    numeric = rs.latest_reference("warhead_classes").name
    biggest = max(int(n.rsplit("_", 1)[1].split(".")[0]) for n in names)
    assert numeric.endswith(f"_{biggest}.csv"), \
        f"latest_reference picked {numeric}, but the highest version is {biggest}"
    if lexical != numeric:
        # The interesting case, and the one live today: the two orderings
        # disagree, and the resolver takes the numerically-newest file.
        assert int(numeric.rsplit("_", 1)[1].split(".")[0]) > \
               int(lexical.rsplit("_", 1)[1].split(".")[0])


def test_the_class_that_only_exists_in_the_newest_library_is_reachable():
    """The concrete consequence, not just the ordering.

    `cinnamamide` exists ONLY in `_10`. If resolution regresses, this fails
    with the same symptom the run did.
    """
    from shared import reference_set as rs
    row = rs.load_warhead_row("cinnamamide")
    assert row.mechanism == "michael_addition"
    assert row.reactive_atom_smarts == "[CX3]=[CX3][CX3]=O"


def test_the_guard_can_fail():
    """Prove the regex matches the real defect line, not just clean code."""
    bad = ('fs = sorted(glob.glob(str(REPO / "data/reference/warhead_classes_*.csv")))\n'
           '    d = pd.read_csv(fs[-1])')
    assert _VERSIONED_GLOB.search(bad) and "sorted" in bad and "key=" not in bad, \
        "the regex no longer matches the line that caused this"
    assert _PICKS_ONE.search(bad), "the [-1] selection is no longer detected"

    # And the two shapes that must NOT be flagged, so the guard stays trusted.
    keyed = ('rs = sorted(glob.glob(str(REPO / "data/reference/x_*.csv")), '
             'key=lambda q: int(q.rsplit("_", 1)[1].split(".")[0]))\n    d = rs[-1]')
    assert "key=" in keyed, "a keyed sort must be exempt"
    collection = ('fs = sorted(glob.glob(str(DATA / topic / "agg_s*_*.csv")))\n'
                  '    d = pd.concat([pd.read_csv(f) for f in fs])')
    assert not _PICKS_ONE.search(collection), \
        "a shard collection must not be flagged -- nothing selects one element"
