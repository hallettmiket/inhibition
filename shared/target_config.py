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


def sweep_families(cfg: dict | None = None) -> dict:
    """{family: [warhead_class, ...]} -- the chemistry that earns a simulation.

    Raises if the scope is absent rather than defaulting to "everything". A run
    that silently widened its own scope would spend the budget on classes nobody
    chose, and the shortlist would not say so.
    """
    fam = get("sweep_rule.scope.families", cfg, default=None)
    if not fam:
        raise ConfigError(
            "sweep_rule.scope.families is not set. The sweep scope is a campaign "
            "decision -- which warhead chemistry the lab will synthesise -- and "
            "it has no safe default. See config/target.yaml.")
    return {str(k): [str(c) for c in v] for k, v in fam.items()}


def sweep_classes(cfg: dict | None = None) -> set:
    """Every warhead class in scope, flattened."""
    return {c for cs in sweep_families(cfg).values() for c in cs}


def family_of(cfg: dict | None = None) -> dict:
    """{warhead_class: family}. Inverted once so callers never re-derive it."""
    out = {}
    for fam, classes in sweep_families(cfg).items():
        for c in classes:
            if c in out:
                raise ConfigError(f"{c} appears in two families: {out[c]}, {fam}")
            out[c] = fam
    return out


def sweep_budget_floor(cfg: dict | None = None) -> float:
    """The enrichment below which we decline to spend GPU time.

    DISTINCT FROM `sweep_floor`, and the distinction is the point. `sweep_floor`
    is the capture-validated threshold the pilot measures and it refuses until
    that pilot has run. This one is a spending rule chosen from the ranking's own
    enrichment distribution and the GPU budget. Sharing a name would let a
    budget decision be reported as a chemistry result.
    """
    v = get("sweep_rule.budget_floor", cfg, default=None)
    if v is None:
        raise ConfigError("sweep_rule.budget_floor is not set; see config/target.yaml")
    return float(v)


def sweep_max_depth(cfg: dict | None = None) -> int:
    """Modes per FAMILY that may be swept. A budget ceiling, not a threshold.

    Distinct from `sweep_floor`, and both apply: the floor is a statement about
    chemistry (this mode is not worth simulating), the cap is a statement about
    money (there is no more GPU time). A caller must be able to report which one
    stopped it, so they are separate numbers with separate readers.
    """
    v = get("sweep_rule.max_depth", cfg, default=None)
    if v is None:
        raise ConfigError("sweep_rule.max_depth is not set; see config/target.yaml")
    return int(v)


def _scope_line(c: dict) -> str:
    """One line naming the families in scope and the depth cap, or UNSET."""
    sr = c.get("sweep_rule", {}) or {}
    fam = (sr.get("scope", {}) or {}).get("families", {}) or {}
    shown = ", ".join(f"{k} ({len(v)})" for k, v in fam.items()) or "UNSET"
    return f"{shown}  max_depth={sr.get('max_depth', 'UNSET')}/family"


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
        f"sweep scope {_scope_line(c)}",
        f"md          sweep {m.get('sweep_ps')} ps, prod {m.get('production_ps')} ps, "
        f"salt {m.get('salt_molar')} M",
        f"chemistry   docked species {c.get('chemistry', {}).get('docked_species')}",
    ])
