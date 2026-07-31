"""
Purpose: resolve a seed name to its run identity, once, for every stage that needs it.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-31
Input: config/seeds.yaml
Output: a validated seed record (SMILES, radius, experiment directory)

WHY THIS IS SHARED. T_2 has four stages and every one of them needs to know
which experiment directory a seed writes to. Four copies of that lookup is four
chances for one stage to resolve a seed differently from the next -- which for
this pipeline means stage 3 docking the frames stage 2 did not annotate, with
nothing raising. The codebase already makes this argument for
`noncovalent_dock_run` and `covalent_dock_run`; the same reasoning applies here.

TWO THINGS ARE REFUSED RATHER THAN DEFAULTED:

* **An unpinned radius.** The usable radius is a property of the SEED, not the
  method (D0018). ATRA yields ZERO at radius 3 while every other declared seed
  yields thousands. A seed that inherits another's radius is how this approach
  enumerates an empty frontier and reports success.

* **A missing experiment directory.** Two seeds sharing one directory interleave
  their frames in a single D2 series, and `latest_frame` -- which keys on the
  directory -- would then hand a later stage whichever seed happened to run
  last. Silent, and indistinguishable from a correct run.
"""

from __future__ import annotations

from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parent.parent
SEEDS_YAML = REPO / "config" / "seeds.yaml"

EXPERIMENT_KEY = {"t2": "t2_experiment"}


class SeedError(RuntimeError):
    """A seed is unknown, or not admissible for the approach asking for it."""


def load_all() -> dict:
    return yaml.safe_load(SEEDS_YAML.read_text(encoding="utf-8"))["seeds"]


def declared_for(approach: str) -> list[str]:
    """Every seed name this approach may run from, in file order."""
    return [k for k, v in load_all().items()
            if approach in (v.get("used_by") or [])]


def resolve(approach: str, name: str, *, require_radius: bool = True) -> dict:
    """One seed's block, checked for admissibility under `approach`.

    Parameters
    ----------
    approach:
        The approach id asking, e.g. ``"t2"``.
    name:
        Key in ``config/seeds.yaml``.
    require_radius:
        ``False`` only for the ``--probe`` path, which is the step that
        DETERMINES the radius and so cannot demand one already be pinned.

    Returns
    -------
    dict
        The seed block, with ``experiment`` added.

    Raises
    ------
    SeedError
        If the seed is unknown, not declared for this approach, has no
        experiment directory, or (unless probing) no pinned radius.
    """
    seeds = load_all()
    if name not in seeds:
        raise SeedError(
            f"unknown seed {name!r}. Declared for {approach}: "
            f"{declared_for(approach)}")
    rec = dict(seeds[name])

    if approach not in (rec.get("used_by") or []):
        raise SeedError(
            f"seed {name!r} is not declared for {approach} "
            f"(used_by = {rec.get('used_by')}). Add it to used_by "
            "deliberately: T_3/T_4 seeds carry a protected core that T_2 would "
            "ignore, and T_2 seeds carry none that T_3/T_4 require.")

    key = EXPERIMENT_KEY.get(approach)
    if key is None:
        raise SeedError(
            f"{approach} has no per-seed experiment convention. Reseeding it "
            "means authoring a core decomposition, not editing config — see "
            "the extension_note in seeds.yaml.")
    if not rec.get(key):
        raise SeedError(
            f"seed {name!r} has no {key}. Give it its own directory; sharing "
            "one with another seed interleaves their frames and latest_frame() "
            "will silently return the wrong seed's work.")
    rec["experiment"] = rec[key]

    if require_radius and rec.get("radius") is None:
        raise SeedError(
            f"seed {name!r} has no pinned radius. Run "
            f"`01_generate.py --seed {name} --probe` and pin the result "
            "(D0018) — do not inherit another seed's radius.")
    return rec


def add_seed_argument(parser, approach: str = "t2") -> None:
    """The identical ``--seed`` flag on every stage of an approach."""
    parser.add_argument(
        "--seed", default=None,
        help=f"which seed in config/seeds.yaml to operate on. Declared for "
             f"{approach}: {', '.join(declared_for(approach))}. Each seed has "
             "its own experiment directory; default is the approach config's "
             "seed.")
