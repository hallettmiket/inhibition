"""
Purpose: guard the machinery the pre-registered elevation experiment added.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06

THREE THINGS ARE UNDER TEST, and each one has a failure mode that produces a
COMPLETE, PLAUSIBLE result rather than an error.

1. **Stage reuse.** Tier 1 measures unrestrained equilibration and tier 2 runs
   the bias from that same frame, which is only true if tier 2 REUSES tier 1's
   NVT/NPT. Reuse keyed on "a .gro exists" would pool a 100 ps verification run
   into a production protocol; reuse offered to a biased stage would reuse a
   trajectory biased along a DIFFERENT pair of atoms, since two PLUMED inputs
   can share an .mdp. Both are refused, and that refusal is what is tested.

2. **Reading a distance out of a .gro.** The serials come from parmed via
   `solv.prmtop`; the file was written by GROMACS several stages later. Column
   parsing, the atom-name identity check and the periodic minimum image are each
   tested on hand-built files, because an integer carried across a file-format
   boundary is this project's signature defect.

3. **The wall's headroom.** The convergence run died with METAD indexing off its
   own grid. The fix bounds the CV with UPPER_WALLS, and a harmonic wall does not
   forbid an overshoot -- it charges for one. The overshoot is `sqrt(2B/kappa)`
   for a standing bias B, so the constants are only safe if that stays inside
   GRID_MAX. That is arithmetic, and it is checked here rather than rediscovered
   at 3 ns.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

from shared import bpmd                                      # noqa: E402
from shared import gromacs_explicit as gx                    # noqa: E402


# --------------------------------------------------------------------------
# 1. stage reuse
# --------------------------------------------------------------------------

@pytest.fixture
def no_gromacs(monkeypatch):
    """Record the stages that would have been executed, and execute none of them."""
    ran: list[str] = []

    def fake_run(cmd, cwd, log_name, timeout=86400, env=None):
        ran.append(log_name)
        if log_name.startswith("mdrun_"):
            name = log_name[len("mdrun_"):].removesuffix(".log")
            (Path(cwd) / f"{name}.gro").write_text("run\n")
        (Path(cwd) / log_name).write_text("")

    monkeypatch.setattr(gx, "_run", fake_run)
    monkeypatch.setattr(gx, "_bin", lambda env, tool: tool)
    return ran


def test_stage_reuses_only_on_an_identical_mdp(tmp_path, no_gromacs):
    gx._stage(tmp_path, "nvt", "nsteps = 50000\n", "min.gro", None, None, 8,
              reuse=True)
    assert no_gromacs == ["grompp_nvt.log", "mdrun_nvt.log"]

    no_gromacs.clear()
    r = gx._stage(tmp_path, "nvt", "nsteps = 50000\n", "min.gro", None, None, 8,
                  reuse=True)
    assert r.reused and no_gromacs == [], "an identical protocol was re-run"

    no_gromacs.clear()
    r = gx._stage(tmp_path, "nvt", "nsteps = 999\n", "min.gro", None, None, 8,
                  reuse=True)
    assert not r.reused and no_gromacs, \
        "a DIFFERENT protocol was served from disk: a short verification run " \
        "would be promoted into a production one"


def test_stage_never_reuses_a_biased_run(tmp_path, no_gromacs):
    """Two PLUMED inputs can share an .mdp and bias different atoms."""
    mdp = "nsteps = 1500000\n"
    gx._stage(tmp_path, "prod", mdp, "npt.gro", None, None, 8,
              plumed="d: DISTANCE ATOMS=10,20\n", reuse=True)
    no_gromacs.clear()
    r = gx._stage(tmp_path, "prod", mdp, "npt.gro", None, None, 8,
                  plumed="d: DISTANCE ATOMS=99,20\n", reuse=True)
    assert not r.reused and no_gromacs, \
        "a run biased along a different atom pair was reused because the .mdp matched"


def test_stage_reuse_is_off_by_default(tmp_path, no_gromacs):
    gx._stage(tmp_path, "nvt", "nsteps = 50000\n", "min.gro", None, None, 8)
    no_gromacs.clear()
    r = gx._stage(tmp_path, "nvt", "nsteps = 50000\n", "min.gro", None, None, 8)
    assert not r.reused and no_gromacs


def test_stop_after_refuses_to_stop_at_production():
    """Production is where the bias lives; a run that skipped it is not a BPMD replicate."""
    with pytest.raises(gx.GromacsError, match="only the equilibration stages"):
        gx.run_pipeline(Path("/nonexistent"), Path("/nonexistent"),
                        stop_after="prod")
    with pytest.raises(gx.GromacsError):
        gx.run_pipeline(Path("/nonexistent"), Path("/nonexistent"),
                        stop_after="min")


def test_stopped_pipeline_never_claims_it_was_biased(tmp_path, monkeypatch,
                                                     no_gromacs):
    """`bpmd_run.run_replicate` refuses a replicate on `plumed`; it must not read True here."""
    wd = tmp_path / "md"
    (wd / "rep1").mkdir(parents=True)
    for f in ("sys.top", "sys.gro", "min.gro"):
        (wd / f).write_text("x\n")
    monkeypatch.setattr(gx, "solvate", lambda src, w: {})
    monkeypatch.setattr(
        gx, "run_pipeline", gx.run_pipeline)          # keep the real one

    class FakeParmed:
        atoms, residues = [], []
        box = [1.0, 1.0, 1.0]
    monkeypatch.setitem(sys.modules, "parmed",
                        type("m", (), {"load_file": staticmethod(
                            lambda *a, **k: FakeParmed())}))

    res = gx.run_pipeline(tmp_path, wd, production_ps=3000.0, replicate=1,
                          candidate_id="x", plumed="d: DISTANCE ATOMS=1,2\n",
                          stop_after="npt")
    assert res["plumed"] is False
    assert res["stopped_after"] == "npt"
    assert res["trajectory"] is None
    assert res["production_ps"] == 0.0
    assert "mdrun_prod.log" not in no_gromacs, "production ran despite stop_after"


# --------------------------------------------------------------------------
# 2. reading a distance out of a .gro
# --------------------------------------------------------------------------

import elevation_run as er                                   # noqa: E402


def write_gro(path: Path, atoms, box=(5.0, 5.0, 5.0)) -> Path:
    """A .gro in the real fixed-width format: %5d%-5s%5s%5d%8.3f%8.3f%8.3f."""
    lines = ["test system", f"{len(atoms):5d}"]
    for i, (resn, name, xyz) in enumerate(atoms, 1):
        lines.append(f"{i:5d}{resn:<5s}{name:>5s}{i:5d}"
                     f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}")
    lines.append("".join(f"{b:10.5f}" for b in box))
    path.write_text("\n".join(lines) + "\n")
    return path


def test_distance_is_read_by_column_and_checks_the_atom_names(tmp_path):
    g = write_gro(tmp_path / "npt.gro", [
        ("LIG", "C7", (1.000, 1.000, 1.000)),
        ("CYS", "SG", (1.300, 1.400, 1.000)),
    ])
    r = er.distance_nm(g, 0, 1, ("C7", "SG"))
    assert r["distance_nm"] == pytest.approx(0.5, abs=1e-6)
    assert not r["pbc_wrapped"]


def test_distance_refuses_when_the_serial_names_another_atom(tmp_path):
    """The exact failure the serials are exposed to: right integer, wrong file."""
    g = write_gro(tmp_path / "npt.gro", [
        ("LIG", "C7", (1.0, 1.0, 1.0)),
        ("CYS", "CB", (1.3, 1.4, 1.0)),
        ("CYS", "SG", (1.3, 1.4, 1.2)),
    ])
    with pytest.raises(er.ElevationError, match="do not name the CV's atoms"):
        er.distance_nm(g, 0, 1, ("C7", "SG"))


def test_distance_takes_the_minimum_image_and_says_so(tmp_path):
    """A wrapped ligand's raw distance is the box length minus the real one."""
    g = write_gro(tmp_path / "npt.gro", [
        ("LIG", "C7", (0.100, 1.0, 1.0)),
        ("CYS", "SG", (4.900, 1.0, 1.0)),
    ], box=(5.0, 5.0, 5.0))
    r = er.distance_nm(g, 0, 1, ("C7", "SG"))
    assert r["raw_distance_nm"] == pytest.approx(4.8)
    assert r["distance_nm"] == pytest.approx(0.2)
    assert r["pbc_wrapped"] is True


def test_distance_refuses_a_triclinic_box(tmp_path):
    g = write_gro(tmp_path / "npt.gro", [("LIG", "C7", (1.0, 1.0, 1.0)),
                                         ("CYS", "SG", (1.5, 1.0, 1.0))])
    lines = g.read_text().splitlines()
    lines[-1] = "   5.00000   5.00000   5.00000   0.00000   0.00000   1.20000"
    g.write_text("\n".join(lines) + "\n")
    with pytest.raises(er.ElevationError, match="triclinic"):
        er.distance_nm(g, 0, 1, ("C7", "SG"))


def test_gro_columns_survive_a_five_digit_residue_number(tmp_path):
    """`split()` would merge residue number and name here and shift every coordinate."""
    p = tmp_path / "big.gro"
    rows = "".join(
        f"{n:5d}{resn:<5s}{name:>5s}{n:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n"
        for n, resn, name, (x, y, z) in
        ((12345, "WAT", "OW", (1.0, 2.0, 3.0)),
         (12346, "LIG", "C7", (1.3, 2.4, 3.0))))
    # No whitespace survives between the residue-name and atom-name fields at
    # this width -- "12346LIG     C712346" -- which is exactly what breaks split().
    assert "C712346" in rows
    p.write_text("t\n    2\n" + rows + "   9.00000   9.00000   9.00000\n")
    names, xyz, box = er.read_gro(p)
    assert names == [("WAT", "OW"), ("LIG", "C7")]
    assert xyz[1] == pytest.approx([1.3, 2.4, 3.0])


# --------------------------------------------------------------------------
# 3. the wall's headroom
# --------------------------------------------------------------------------

def test_the_wall_sits_past_unbound_so_it_cannot_alter_the_barrier():
    assert bpmd.WALL_NM > bpmd.UNBOUND_NM, \
        "a wall inside the unbound threshold would be part of the escape barrier " \
        "being measured, not a guard on the grid"


def test_the_grid_can_absorb_any_overshoot_well_tempered_can_pay_for():
    """A harmonic wall charges for an overshoot rather than forbidding one.

    Measured on this branch: with 100x the production deposition rate the CV sat
    0.206 nm past a wall at 0.5 nm under a standing bias of 93.5 kJ/mol, against
    the harmonic bound sqrt(2B/kappa) = 0.306 nm. So the bound holds and is
    conservative, and the question is only whether the grid can absorb it.
    """
    headroom = bpmd.GRID_MAX_NM - bpmd.WALL_NM
    affordable = 0.5 * bpmd.WALL_KAPPA_KJ_PER_NM2 * headroom ** 2

    # Well-tempered plateaus at (gamma-1)/gamma * dF. Even at an implausible
    # 200 kJ/mol escape barrier the bias cannot reach `affordable`.
    plateau = (bpmd.BIASFACTOR - 1) / bpmd.BIASFACTOR * 200.0
    assert plateau < affordable / 4, (
        f"a {plateau:.0f} kJ/mol bias overshoots the wall by "
        f"{np.sqrt(2*plateau/bpmd.WALL_KAPPA_KJ_PER_NM2):.2f} nm against "
        f"{headroom:.2f} nm of grid — this is how the convergence run died")


def test_a_measured_overshoot_stays_within_the_harmonic_bound():
    b, kappa, observed = 93.5, bpmd.WALL_KAPPA_KJ_PER_NM2, 0.206
    assert observed <= np.sqrt(2 * b / kappa) + 1e-9


# --------------------------------------------------------------------------
# the anchor
# --------------------------------------------------------------------------

class _C:
    def __init__(self, ident):
        self.ident = ident


def test_the_anchor_is_deterministic_and_capped():
    by_id = {f"xtal:{p}:L": _C(f"xtal:{p}:L")
             for p in ("9V6W", "6VAJ", "7EFJ", "9INN", "7F0M", "9INO",
                       "7EFX", "7EKV", "9INP", "9INQ")}
    by_id["t4_abc"] = _C("t4_abc")
    got = [c.ident for c in er.reference_positives(by_id)]
    assert len(got) == er.N_REF
    assert got == sorted(got), "the anchor's membership must not depend on dict order"
    assert all(i.startswith("xtal:") for i in got)
    assert er.reference_positives(by_id)[0].ident == got[0]


def test_the_anchor_is_a_refusal_not_a_warning():
    """The prereg says the between-group comparison is uninterpretable without it."""
    with pytest.raises(er.ElevationError, match="uninterpretable without the anchor"):
        er.reference_positives({"t4_abc": _C("t4_abc")})


# --------------------------------------------------------------------------
# the analysis: direction, the n = 5 rule, and refusing to read an unfinished run
# --------------------------------------------------------------------------

import elevation_analysis as ea                              # noqa: E402


def _per(**groups) -> "pd.DataFrame":
    import pandas as pd
    rows = []
    for g, vals in groups.items():
        for i, v in enumerate(vals):
            rows.append({"group": g, "ident": f"{g}_{i}", "tier1": v,
                         "tier2": v, "n_replicas": 3})
    return pd.DataFrame(rows)


A, B, D = ea.GROUPS[0], ea.GROUPS[1], ea.GROUPS[2]
V, R = ea.GROUPS[3], ea.GROUPS[4]


def test_cliffs_delta_is_the_pairwise_statement_it_claims_to_be():
    assert ea.cliffs_delta(np.array([4., 5., 6.]), np.array([1., 2., 3.])) == 1.0
    assert ea.cliffs_delta(np.array([1., 2., 3.]), np.array([4., 5., 6.])) == -1.0
    assert ea.cliffs_delta(np.array([1., 2.]), np.array([1., 2.])) == 0.0


def test_tier1_direction_is_flipped_so_smaller_displacement_reads_as_more_stable():
    """Tier 1 is a displacement and tier 2 is a score; they point opposite ways."""
    per = _per(**{A: [0.1] * 8, B: [0.9] * 8})
    t1 = ea.compare(per, "tier1", A, B, "tier1")
    t2 = ea.compare(per, "tier2", A, B, "tier2")
    assert t1["cliffs_delta_more_stable"] == 1.0, \
        "the group that moved LESS must read as more stable on tier 1"
    assert t2["cliffs_delta_more_stable"] == -1.0, \
        "the group with the LOWER score must read as less stable on tier 2"


def test_no_p_value_is_produced_for_the_n_equals_5_group():
    """The prereg forbids it, and the way that survives is the number not existing."""
    per = _per(**{V: [0.1] * 5, R: [0.9] * 8})
    r = ea.compare(per, "tier1", V, R, "tier1")
    assert not np.isfinite(r["p"])
    assert "n = 5" in r["note"]
    # The effect size IS reported -- the prereg calls that arm descriptive, not absent.
    assert np.isfinite(r["cliffs_delta_more_stable"])


def test_a_normal_contrast_still_gets_its_p_value():
    per = _per(**{A: [0.1] * 8, B: [0.9] * 8})
    r = ea.compare(per, "tier1", A, B, "tier1")
    assert np.isfinite(r["p"]) and r["p"] < 0.05


def test_holm_is_monotone_and_leaves_nans_alone():
    assert ea.holm([0.01, 0.04, 0.5]) == pytest.approx([0.03, 0.08, 0.5])
    out = ea.holm([0.01, np.nan])
    assert out[0] == pytest.approx(0.01) and not np.isfinite(out[1])


def test_an_unfinished_run_cannot_emit_a_reading():
    """Empty groups make every contrast "~", which is a substantive conclusion."""
    per = _per(**{A: [0.2] * 8})
    con = ea.contrast_table(per, "tier1", "tier1")
    r = ea.readings(con, per, "tier1", "tier1")
    assert r["reading"].startswith("NOT READABLE")


def test_the_wrong_direction_is_reported_as_a_failure_not_reinterpreted():
    """D more stable than A and B is the prereg's fourth reading."""
    per = _per(**{A: [0.9] * 8, B: [0.9] * 8, D: [0.1] * 8})
    con = ea.contrast_table(per, "tier1", "tier1")
    r = ea.readings(con, per, "tier1", "tier1")
    assert r["reading"] == "D >= A, B"
    assert "do not reinterpret" in r["conclusion"]


def test_consensus_is_the_filter_reading():
    per = _per(**{A: [0.1] * 8, B: [0.1] * 8, D: [0.9] * 8})
    con = ea.contrast_table(per, "tier1", "tier1")
    r = ea.readings(con, per, "tier1", "tier1")
    assert r["reading"] == "B ~ A, both > D"
    assert "Consensus is the filter" in r["conclusion"]


def test_enrichment_adds_something_reading():
    per = _per(**{A: [0.1] * 8, B: [0.5] * 8, D: [0.9] * 8})
    con = ea.contrast_table(per, "tier1", "tier1")
    r = ea.readings(con, per, "tier1", "tier1")
    assert r["reading"] == "A > B, both > D"


def test_no_difference_reading():
    rng = np.random.default_rng(0)
    per = _per(**{A: rng.normal(0.5, .01, 8), B: rng.normal(0.5, .01, 8),
                  D: rng.normal(0.5, .01, 8)})
    con = ea.contrast_table(per, "tier1", "tier1")
    assert ea.readings(con, per, "tier1", "tier1")["reading"] == "A ~ B ~ D"


def test_a_large_effect_without_significance_is_not_promoted_to_a_difference():
    """Both halves of the rule are required; either alone is misreadable at n = 8."""
    per = _per(**{A: [0.1, 0.2, 0.3], B: [0.4, 0.5, 0.6],
                  D: [0.4, 0.5, 0.6]})
    con = ea.contrast_table(per, "tier1", "tier1")
    ab = con[(con.group_1 == "A") & (con.group_2 == "B")].iloc[0]
    assert ab.cliffs_delta_more_stable == 1.0, "the effect is as large as it gets"
    assert ab.verdict == "~", "n = 3 cannot reach significance and must not claim it"
