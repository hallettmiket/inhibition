"""A run's tables and its poses must land in the same topic.

THE DEFECT THIS PINS. `--topic` used to rebind the table topic ONLY, and only
when the topic differed from the default; the two pose directories were module
constants naming `nac_v3` outright. So `--topic nac_v4` sent tables to nac_v4 and
poses to nac_v3_poses -- where every molecule already had a file, so the
`if not sdf.exists()` append-only guard skipped all of them. The screen ran to
completion, reported success, and wrote no poses.

It is silent by construction: the guard is doing exactly what it should, the run
exits 0, and the only evidence is an mtime. Measured mid-3.0.0-run, 5,772 of
5,774 representative files were still the Aug-07 2.1.0 ones -- one pose each --
while the tables described modes from fresh 500-run clouds. That mismatch is the
same one 2.2.0 was written off for, so the cost of missing it once is a screen.

These tests are cheap and the failure is expensive, which is the whole trade.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))


def _load():
    """Import the screen without running it."""
    sys.path.insert(0, str(REPO / "scripts"))
    spec = importlib.util.spec_from_file_location(
        "nac_screen_v2", REPO / "scripts" / "nac_screen_v2.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def m():
    return _load()


@pytest.mark.parametrize("topic", ["nac_v3", "nac_v4", "nac_v9", "pilot_x"])
def test_all_three_outputs_follow_the_topic(m, topic):
    out, poses, allposes = m.topic_paths(topic)
    assert poses.name == f"{topic}_poses"
    assert allposes.name == f"{topic}_allposes"
    # The table topic must agree too, whatever attribute carries it.
    assert topic in repr(out.__dict__) or topic in str(out.__dict__.values())


def test_the_default_topic_is_not_special_cased(m):
    """The old bug only bit on a NON-default topic, so test the boundary."""
    a = m.topic_paths("nac_v3")
    b = m.topic_paths("nac_v4")
    assert a[1] != b[1] and a[2] != b[2], \
        "two topics resolved to the same pose directory"


def test_poses_and_tables_cannot_be_given_different_topics(m):
    """Derivation from one argument is the invariant -- not three call sites."""
    import inspect
    src = inspect.getsource(m.topic_paths)
    assert src.count("topic") >= 3, "a path stopped being derived from `topic`"


def test_no_hardcoded_topic_directory_survives_in_the_source():
    """The literal that caused it must not come back anywhere in the screen."""
    src = (REPO / "scripts" / "nac_screen_v2.py").read_text()
    # Prose is allowed to name the directories -- the docstrings above explain
    # the defect and must be able to. Code is not. Comment lines are stripped
    # and the remainder is searched for the BARE substring, not for a
    # quote-prefixed one: the constant that caused this was
    # `Path(".../nac_v3_poses")`, where the topic sits mid-string and a
    # `'"' + name` test sails straight past it.
    code, in_doc = [], False
    for ln in src.splitlines():
        s = ln.strip()
        if s.startswith(('"""', "'''")):
            # naive but sufficient: docstrings here open and close on their own
            # lines or are single-line
            in_doc = not in_doc if s.count('"""') % 2 else in_doc
            continue
        if in_doc or s.startswith(("#", "#:")):
            continue
        code.append(ln)
    code = "\n".join(code)
    for bad in ("nac_v3_poses", "nac_v3_allposes", "nac_v2_poses"):
        assert bad not in code, \
            f"{bad} is hardcoded in code again -- this is the 2.2.0 defect"


# ---------------------------------------------------------------------------
# The same defect, one caller out.
#
# Everything above reads `nac_screen_v2.py` and nothing else. But `one()` reads
# POSE_DIR and ALL_POSE_DIR off the MODULE when it writes, so any script that
# drives `one()` directly has to move them itself -- and `dock_reference_modes`
# set `nsv.OUT` alone, sending a reference compound's tables to its own topic
# and its 500-pose cloud into the production run's `<topic>_allposes`.
#
# The guard above could not fail on that, because it never looked outside one
# file. This is that scope hole, closed: a partial rebind is now a test failure
# wherever it is written.
# ---------------------------------------------------------------------------

_PATHS = ("OUT", "POSE_DIR", "ALL_POSE_DIR", "TOPIC")


def _callers():
    """Every repo file that imports the screen as a module."""
    out = []
    for d in ("scripts", "shared", "integration", "approaches"):
        for f in (REPO / d).rglob("*.py") if (REPO / d).is_dir() else []:
            if f.name == "nac_screen_v2.py":
                continue
            src = f.read_text(errors="ignore")
            if "nac_screen_v2" in src:
                out.append((f, src))
    return out


def test_no_caller_rebinds_the_screens_output_paths_piecemeal():
    """`nsv.OUT = ...` without the other three is the whole bug.

    Assigning ANY of the four individually is banned outright rather than
    checking that all four are assigned together: 'all four are present' is a
    condition a future edit deletes one line from and nothing complains. There
    is exactly one supported way to move them, `use_topic()`, so anything else
    is the defect regardless of how complete it looks.
    """
    import re
    offenders = []
    for f, src in _callers():
        alias = re.findall(r"import\s+nac_screen_v2\s+as\s+(\w+)", src)
        alias.append("nac_screen_v2")
        for ln_no, ln in enumerate(src.splitlines(), 1):
            code = ln.split("#")[0]
            for a in alias:
                for attr in _PATHS:
                    if re.search(rf"\b{a}\.{attr}\s*=(?!=)", code):
                        offenders.append(f"{f.relative_to(REPO)}:{ln_no}: {ln.strip()}")
    assert not offenders, (
        "these assign the screen's output paths directly; call "
        "nac_screen_v2.use_topic(topic) instead so all four move together:\n  "
        + "\n  ".join(offenders))


def test_use_topic_moves_all_four(m, tmp_path):
    """The replacement must actually be the thing that cannot be half-done."""
    m.use_topic("probe_topic_a")
    assert m.TOPIC == "probe_topic_a"
    assert m.POSE_DIR.name == "probe_topic_a_poses"
    assert m.ALL_POSE_DIR.name == "probe_topic_a_allposes"
    before = (m.OUT, m.POSE_DIR, m.ALL_POSE_DIR)
    m.use_topic("probe_topic_b")
    assert m.POSE_DIR.name == "probe_topic_b_poses"
    assert m.ALL_POSE_DIR.name == "probe_topic_b_allposes"
    assert (m.OUT, m.POSE_DIR, m.ALL_POSE_DIR) != before, \
        "use_topic() did not move the paths"


def test_the_guard_can_fail(tmp_path):
    """A vacuous guard is this project's most common bug -- prove this one bites.

    Written because two tests in this repo passed for free: one counted zero
    accesses, the other flagged only docstrings.
    """
    import re
    bad = "import nac_screen_v2 as nsv\nnsv.OUT = something\n"
    alias = re.findall(r"import\s+nac_screen_v2\s+as\s+(\w+)", bad)
    hit = any(re.search(rf"\b{a}\.OUT\s*=(?!=)", ln.split("#")[0])
              for a in alias for ln in bad.splitlines())
    assert hit, "the pattern the guard searches for does not match the defect"
