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

import copy as _copy
import functools as _functools
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parent.parent
CONFIG = REPO / "config" / "target.yaml"


class ConfigError(RuntimeError):
    """The config cannot answer this. Named so it cannot pass as a value."""


def load(path: Path | None = None) -> dict:
    """The target config. CACHED ON (path, mtime, size).

    `_cfg` in the screen reads single keys, and each read re-parsed the whole
    YAML: 196 loads and 8 s per molecule, profiled 2026-08-27.

    KEYED ON THE FILE'S MTIME AND SIZE, not just its path. A cache on the path
    alone would serve a stale config to a long run whose file was edited
    mid-flight -- and two different configs used inside one run is exactly the
    D0080 defect. Returns a DEEP COPY so a caller mutating the result cannot
    poison every later reader.
    """
    p = Path(path or CONFIG)
    if not p.is_file():
        raise ConfigError(f"no target config at {p}")
    st = p.stat()
    return _copy.deepcopy(_load_cached(str(p.resolve()), st.st_mtime_ns, st.st_size))


@_functools.lru_cache(maxsize=8)
def _load_cached(path_str: str, mtime_ns: int, size: int) -> dict:
    import yaml
    with open(path_str) as fh:
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


def sweep_min_mode_poses(cfg: dict | None = None) -> int:
    """Poses a mode needs before its enrichment is worth acting on.

    Asked BEFORE the floor, and it is a different question: the floor asks
    whether the value is high, this asks whether the value means anything. A
    two-pose mode scores the arithmetic maximum whenever both poses are viable.
    """
    v = get("sweep_rule.min_mode_poses", cfg, default=None)
    if v is None:
        raise ConfigError("sweep_rule.min_mode_poses is not set; see config/target.yaml")
    return int(v)


def rank_min_mode_poses(cfg: dict | None = None) -> int:
    """Poses a mode needs before it may hold a RANK (#65).

    THE SAME QUESTION AS `sweep_min_mode_poses`, ASKED ONE STAGE EARLIER, and
    kept as its own key so the two can be answered differently without either
    silently following the other. Ranking asks whether a mode may be ORDERED
    against its class; sweeping asks whether it may be SIMULATED. That they
    currently share the value 12 is a fact about this target, asserted by a test
    rather than achieved by aliasing.

    It replaces `consensus >= 0.05`. `consensus` is mode_size / n_poses, so a
    fraction floor is a size floor divided by the cloud -- on a 500-pose cloud
    exactly "at least 25 poses", a number nothing measured, and one that moves
    to 50 the moment `docking.n_runs` doubles without any output saying so.
    """
    v = get("ranking.mode_gate.min_poses", cfg, default=None)
    if v is None:
        raise ConfigError(
            "ranking.mode_gate.min_poses is not set. The rank gate is a "
            "measured estimability threshold, not a default: see "
            "config/target.yaml and D0084.")
    return int(v)


def rank_gate_parameter(cfg: dict | None = None) -> str:
    """The column the rank gate is applied to. Named so it cannot drift silently."""
    v = get("ranking.mode_gate.parameter", cfg, default=None)
    if v is None:
        raise ConfigError("ranking.mode_gate.parameter is not set; see config/target.yaml")
    return str(v)


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


def topic(cfg: dict | None = None) -> str:
    """The topic every output of a run is keyed on -- ONE source of truth.

    THE DEFAULT WAS THE BUG. Five scripts each carried their own literal --
    `nac_screen_v2` and `score_selection` said "nac_v3", `rank_v2` said
    "nac_v2", `sweep_gap_worklist` said "nac_v4" -- while `config/target.yaml`
    said something else again, with a comment reading "bump with the screen,
    never on its own". So a re-run launched without an explicit --topic wrote
    the SCREEN to nac_v3 and the WORKLIST to nac_v4, and nothing announced the
    split; D0080 exists because that already happened once.

    A default that lives in five places is five defaults. This is the one, and
    bumping `run.topic` is now the whole ceremony for starting a fresh screen.
    """
    return str(get("run.topic", cfg, default="nac_v4"))


def md_production_ps(cfg: dict | None = None) -> float:
    """Production length of a full MD run, in ps.

    Read by `pipeline_schematic` -- the DIAGRAM -- and by nothing that runs. The
    runner's own default was 100.0 ps against a config saying 100_000, so the
    page stated 100 ns while the code would have produced 0.1 ns, and the real
    value came from a flag typed into a scratch shell script. A spec the runner
    does not read is documentation.
    """
    return float(get("md.production_ps", cfg, default=100_000.0))


def md_sweep_ps(cfg: dict | None = None) -> float:
    """Triage-sweep length, in ps (D0085: 8 ns)."""
    return float(get("md.sweep_ps", cfg, default=8_000.0))


def md_replicates(cfg: dict | None = None) -> int:
    """Replicates per production run. Had no reader at all."""
    return int(get("md.replicates", cfg, default=1))


def md_survivor_rmsd_nm(cfg: dict | None = None) -> float:
    """Max ligand RMSD in the 8 ns TRIAGE SWEEP that earns a 100 ns run
    (D0085: 0.35 nm).

    This used to gate the sweep, the 100 ns verdict and BPMD promotion alike,
    "because it is the same question at three timescales". It is not: a ligand
    explores more in 100 ns than in 8, so the same number that selects tight
    8 ns poses makes the 100 ns "optimal" tier unreachable. The production bar
    is `md_production_optimal_rmsd_nm`.
    """
    return float(get("md.sweep_survivor_rmsd_nm", cfg, default=0.35))


def md_production_optimal_rmsd_nm(cfg: dict | None = None) -> float:
    """Max ligand RMSD over a finished 100 ns run that counts as OPTIMAL
    (@tt8804, 2026-08-18: 0.45 nm).

    Held vs left is decided separately, by the 1.0 nm dissociation criterion in
    `shared.residence_tier`; this only splits the runs that never left.
    """
    return float(get("md.production_optimal_rmsd_nm", cfg, default=0.45))


def md_held_residence_floor(cfg: dict | None = None) -> float:
    """Fraction of frames bound that separates a clean hold from an excursion
    that came back (@tt8804, 2026-08-18: 0.95).

    Dissociation (`bound_rmsd_nm`) says whether the ligand ever left for good;
    this says whether it stayed put while it was there.
    """
    return float(get("md.held_residence_floor", cfg, default=0.95))


def md_tier_ps(tier: str, cfg: dict | None = None) -> float:
    """Length in ps of one MD tier: `triage` | `screen` | `production` (D0101).

    AN ALLOWLIST OF THREE, and it RAISES on anything else. The three lengths are
    a filter cascade at ONE bar, and a caller that invents a fourth name would
    otherwise get a default and run a length nobody chose -- the shape #14
    records. `md.sweep_ps` and `md.production_ps` remain the readers for tiers 1
    and 3 so existing callers are untouched; this is the one place that knows the
    cascade as a whole.
    """
    names = {"triage": "triage_ps", "screen": "screen_ps",
             "production": "production_ps"}
    if tier not in names:
        raise ConfigError(f"unknown MD tier {tier!r}; known: {sorted(names)}")
    v = get(f"md.tiers.{names[tier]}", cfg, default=None)
    if v is None:
        raise ConfigError(
            f"md.tiers.{names[tier]} is not set. The three tiers are a measured "
            "cascade (D0101), not defaults: see config/target.yaml.")
    return float(v)
