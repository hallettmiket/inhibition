"""
Purpose: reference files must resolve to the NEWEST version, never a hand-written pin.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-31
Input: shared.reference_set, data/reference/
Output: pass/fail

WHY THIS EXISTS. This project has now written a stale version pin three times,
and every instance was silent because the pinned file still existed and still
parsed:

* `warhead_library.DEFAULT_LIBRARY` stayed on `warhead_classes_7.csv` after
  `_8` was written, so the acrylamide precedent correction in `_8` was loaded
  by nothing at all.
* `gen_docs.py` named `pin1_reference_binders_1.csv` literally and carried a
  load error on the rendered status page for weeks.
* `reference_set.DEFAULT_MASTER` pinned `_3` while `_4` -- adding two potent
  covalent actives the gate had been calling UNDERPOWERED without -- sat
  unread next to it.

A pin cannot announce that it is out of date. Only a comparison against the
directory can.
"""

from __future__ import annotations

import pandas as pd
import pytest

from shared import reference_set as rs

STEMS = ["pin1_reference_binders", "pin1_covalent_cys113_anchors"]

# WHERE TO SCAN. The first version of the stale-pin guard looked only at
# `shared/*.py` and `scripts/*.py`, and matched only the stem
# `pin1_reference_binders`. Both limits were invisible until #9 found three
# stale `warhead_classes` pins by hand -- in `config/`, `approaches/` and
# `integration/`, none of which the guard could see, naming a stem it did not
# check. A guard that cannot reach the code it is guarding is the "scoped out"
# disguise in `docs/how_this_project_breaks.md`, and it is why this now
# discovers its stems from the directory and walks the whole tree.
SCAN_DIRS = ("shared", "scripts", "approaches", "integration", "tests", "config")

# A builder necessarily names its own output. Keyed on filename so that moving
# the file does not silently carry the exemption somewhere it is not warranted.
EXEMPT_FILES = {"build_reference_binders_4.py"}


def _newest_reference_versions() -> dict[str, int]:
    """Every versioned stem under data/reference/, mapped to its highest N.

    Discovered from the directory rather than listed here, for the same reason
    `latest_reference()` globs: a hand-maintained list of stems goes stale in
    exactly the way this module exists to catch.
    """
    newest: dict[str, int] = {}
    for p in rs._REF_DIR.glob("*_*.csv"):
        stem, _, tail = p.stem.rpartition("_")
        if stem and tail.isdigit():
            newest[stem] = max(newest.get(stem, 0), int(tail))
    return newest


def _executable_strings(path):
    """(lineno, text) for literals code can actually evaluate.

    Python is walked by AST so docstrings are exempt -- prose legitimately
    cites old versions, and this module's own docstring cites `_1` and `_7`.
    YAML has no docstrings, so comments are stripped instead.
    """
    import ast

    if path.suffix == ".py":
        tree = ast.parse(path.read_text())
        docstrings = set()
        for node in ast.walk(tree):
            if isinstance(node, (ast.Module, ast.FunctionDef,
                                 ast.AsyncFunctionDef, ast.ClassDef)):
                body = getattr(node, "body", None)
                if (body and isinstance(body[0], ast.Expr)
                        and isinstance(body[0].value, ast.Constant)
                        and isinstance(body[0].value.value, str)):
                    docstrings.add(id(body[0].value))
        for node in ast.walk(tree):
            if (isinstance(node, ast.Constant) and isinstance(node.value, str)
                    and id(node) not in docstrings):
                yield node.lineno, node.value
    else:
        for i, line in enumerate(path.read_text().splitlines(), 1):
            yield i, line.split("#", 1)[0]


@pytest.mark.parametrize("stem", STEMS)
def test_resolver_picks_the_highest_version(stem):
    got = rs.latest_reference(stem)
    versions = []
    for p in (rs._REF_DIR).glob(f"{stem}_*.csv"):
        tail = p.stem[len(stem) + 1:]
        if tail.isdigit():
            versions.append(int(tail))
    assert versions, f"no versioned {stem}_N.csv on disk"
    assert got.stem == f"{stem}_{max(versions)}", (
        f"{stem} resolved to {got.name} but {max(versions)} is the highest "
        "version present")


def test_defaults_are_the_resolved_latest():
    """The module constants must agree with the resolver, not drift from it."""
    assert rs.DEFAULT_MASTER == rs.latest_reference("pin1_reference_binders")
    assert rs.DEFAULT_ANCHORS == rs.latest_reference(
        "pin1_covalent_cys113_anchors")


def test_no_source_file_hardcodes_an_older_reference_version():
    """No evaluated literal may name a superseded `data/reference/` file.

    Covers EVERY versioned stem in `data/reference/` -- not just
    `pin1_reference_binders` -- across Python and YAML in all of SCAN_DIRS.
    Python docstrings and YAML comments are exempt: prose legitimately names
    old versions, and only a literal code evaluates can actually be opened.

    A literal naming the newest version is still latent (it goes stale the
    moment the next version lands), but it is not wrong today, so this fails
    only on genuinely outdated ones.
    """
    import re
    from pathlib import Path

    repo = Path(rs._REPO_ROOT)
    newest = _newest_reference_versions()
    assert newest, "no versioned reference files found; the guard would be vacuous"
    pat = re.compile(r"([a-z0-9_]+?)_(\d+)\.csv")

    offenders = []
    for d in SCAN_DIRS:
        for path in sorted((repo / d).rglob("*")):
            if (path.suffix not in {".py", ".yaml", ".yml"}
                    or path.name in EXEMPT_FILES):
                continue
            for lineno, text in _executable_strings(path):
                for stem, ver in pat.findall(text):
                    if stem in newest and int(ver) < newest[stem]:
                        offenders.append(
                            f"{path.relative_to(repo)}:{lineno} -> "
                            f"{stem}_{ver}.csv (newest is _{newest[stem]})")
    assert not offenders, (
        "hard-coded reference version(s) superseded on disk:\n  "
        + "\n  ".join(offenders)
        + "\nUse reference_set.latest_reference() / warhead_library's default.")


def test_the_master_set_carries_the_structured_potency_columns():
    df = pd.read_csv(rs.DEFAULT_MASTER)
    for col in ("potency", "potency_value_nM", "potency_unit", "potency_type",
                "tier"):
        assert col in df.columns, f"{rs.DEFAULT_MASTER.name} lacks {col}"
    # The free-text source is preserved: the parse is additive, never a
    # replacement, so a mis-parse cannot destroy what the paper said.
    assert df["potency"].notna().any()


def test_tiers_are_a_closed_vocabulary():
    df = pd.read_csv(rs.DEFAULT_MASTER)
    allowed = {"lead", "fragment", "historical_promiscuous"}
    bad = set(df["tier"].dropna()) - allowed
    assert not bad, f"unknown tier value(s) {bad}; allowed: {sorted(allowed)}"


def test_no_promiscuous_compound_is_tiered_as_a_lead():
    """The whole point of the tier column.

    juglone, EGCG, ATRA, KPT-6566 and PiB are promiscuity-flagged; carrying
    them as gate actives plausibly depresses the AUC rather than supporting it.
    """
    df = pd.read_csv(rs.DEFAULT_MASTER)
    flagged = df["promiscuity_flag"].astype(str).str.lower().isin(
        {"y", "yes", "true"})
    bad = df[flagged & (df["tier"] == "lead")]["name"].tolist()
    assert not bad, f"promiscuity-flagged but tiered lead: {bad}"


def test_parsed_potency_is_never_invented():
    """A parsed value must be traceable to the free-text string it came from."""
    df = pd.read_csv(rs.DEFAULT_MASTER)
    got = df[df["potency_value_nM"].notna()]
    assert len(got) > 0
    for _, r in got.iterrows():
        assert isinstance(r["potency"], str) and r["potency"].strip(), (
            f"{r['name']}: parsed a value with no source text")
        assert str(r["potency_type"]).lower() in {"ki", "kd", "ic50", "ec50"}
        assert r["potency_unit"] == "nM"


def test_the_two_recovered_covalent_actives_are_present():
    """Both were verified against RCSB before being added."""
    df = pd.read_csv(rs.DEFAULT_MASTER)
    for name, pdb in (("Liu-2022-ZL-Pin13", "7F0M"),
                      ("Alboreggia-2024-164A10", "8VJG")):
        row = df[df["name"] == name]
        assert len(row) == 1, f"{name} missing from {rs.DEFAULT_MASTER.name}"
        assert str(row["pdb"].iloc[0]).upper() == pdb


def test_every_smiles_still_parses():
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    df = pd.read_csv(rs.DEFAULT_MASTER)
    bad = [r["name"] for _, r in df.iterrows()
           if str(r["canonical_smiles"]) != rs.UNVERIFIED
           and Chem.MolFromSmiles(str(r["canonical_smiles"])) is None]
    assert not bad, f"unparseable SMILES: {bad}"
