"""Read `config/target.yaml`, and refuse to invent anything it does not say.

The tool is target-agnostic, so every target-specific number has to come from one
file rather than from a default buried in whichever script needed it first. This
module is that file's only reader.

IT REFUSES RATHER THAN SUBSTITUTES. `sweep_rule.floor` is null until a pilot has
measured it for the target at hand (#59), and asking for it raises instead of
returning a number inherited from Pin1. A floor is the one parameter here that
cannot travel: it is a property of the chemistry and the criterion, and a
plausible default would be indistinguishable from a measured one in every
artefact downstream.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parent.parent
CONFIG = REPO / "config" / "target.yaml"


class ConfigError(RuntimeError):
    """The config cannot answer this. Named so it cannot pass as a value."""


def load(path: Path | None = None) -> dict:
    import yaml
    p = Path(path or CONFIG)
    if not p.is_file():
        raise ConfigError(f"no target config at {p}")
    with p.open() as fh:
        return yaml.safe_load(fh) or {}


def get(dotted: str, cfg: dict | None = None, default: Any = "__raise__") -> Any:
    """`get("md.salt_molar")`. Missing keys raise unless a default is given."""
    cur = cfg if cfg is not None else load()
    for part in dotted.split("."):
        if not isinstance(cur, dict) or part not in cur:
            if default != "__raise__":
                return default
            raise ConfigError(f"{dotted}: not in the target config")
        cur = cur[part]
    return cur


def sweep_floor(cfg: dict | None = None) -> float:
    """The enrichment floor, or a refusal explaining what has to happen first.

    Never falls back. A screen that silently used Pin1's floor on another target
    would produce a shortlist whose selection rule nobody chose.
    """
    v = get("sweep_rule.floor", cfg, default=None)
    if v is None:
        pilot = get("sweep_rule.pilot", cfg, default={})
        raise ConfigError(
            "sweep_rule.floor is not set for this target. It is measured, not "
            "inherited: run the stratified pilot "
            f"({pilot.get('n_per_stratum', '?')} per stratum over "
            f"{pilot.get('strata', '?')}), fit P(productive | enrichment), and "
            "set the floor at the capture_target. See docs/sweep_depth.md.")
    return float(v)


def summary(cfg: dict | None = None) -> str:
    """One block a run can print so its settings are in its own log."""
    c = cfg if cfg is not None else load()
    t, s, m = c.get("target", {}), c.get("splitting", {}), c.get("md", {})
    st2 = s.get("stage2", {})
    floor = c.get("sweep_rule", {}).get("floor")
    return "\n".join([
        f"target      {t.get('name')} / {t.get('domain')} / {t.get('pdb')} "
        f"/ {t.get('anchor')}",
        f"docking     {c.get('docking', {}).get('n_runs')} runs, "
        f"persist_all_poses={c.get('docking', {}).get('persist_all_poses')}",
        f"splitting   stage2={st2.get('enabled')} cut={st2.get('cut_diameter_a')} A "
        f"max_sub={st2.get('max_sub')}",
        f"sweep rule  {c.get('sweep_rule', {}).get('parameter')} >= "
        f"{floor if floor is not None else 'UNSET (pilot required)'}, "
        f"select_by={c.get('sweep_rule', {}).get('select_by')}",
        f"md          sweep {m.get('sweep_ps')} ps, prod {m.get('production_ps')} ps, "
        f"salt {m.get('salt_molar')} M",
        f"chemistry   docked species {c.get('chemistry', {}).get('docked_species')}",
    ])
