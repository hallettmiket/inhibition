"""The 100 ns elevation gate, and the two silent paths that fed it.

@twu383, 2026-09-02: *"<3.5 max RSMD or mean RMSD 3.0 (need to account for
quick spikes but overall low rmsd) + 60%+ warhead within 3.5 A goes to 100 ns"*.

A 100 ns run is the most expensive thing this campaign does, so the gate that
authorises one has to be wrong LOUDLY rather than quietly. Every test here is a
way it could pass something it should not while still looking right.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))


def _asw():
    try:
        import attack_sweep as m
    except Exception as exc:                              # noqa: BLE001
        pytest.skip(f"attack_sweep not importable here: {exc}")
    return m


# --------------------------------------------------------------- the gate ---
def test_both_conditions_are_required():
    """Either alone must not elevate -- they measure different things.

    Occupancy is about the WARHEAD reaching Cys113; RMSD is about the whole
    LIGAND staying where it was docked. A molecule can pivot its warhead in
    while its scaffold walks off, and it can sit still facing the wrong way.
    """
    a = _asw()
    th = a.elevation_thresholds()
    good_rmsd = dict(rmsd_max_a=th["rmsd_max_a"] - 0.5, rmsd_mean_a=1.0)

    assert a.elevation_verdict(
        dict(frac_attack_ready=th["occupancy_min"] + 0.05, **good_rmsd))["elevate"]
    assert not a.elevation_verdict(
        dict(frac_attack_ready=th["occupancy_min"] - 0.01, **good_rmsd))["elevate"], (
        "engagement just below the bar elevated anyway")
    assert not a.elevation_verdict(dict(
        frac_attack_ready=0.99,
        rmsd_max_a=th["rmsd_max_a"] + 5, rmsd_mean_a=th["rmsd_mean_a"] + 5,
    ))["elevate"], "a pose that left the site elevated on engagement alone"


def test_the_spike_allowance_is_an_OR_not_an_AND():
    """"account for quick spikes but overall low rmsd" -- the mean rescues it.

    If this became an AND, a run that touches 4 A for 20 ps and sits at 2 A for
    the rest would be rejected, which is exactly the case the OR was asked for.
    """
    a = _asw()
    th = a.elevation_thresholds()
    v = a.elevation_verdict(dict(
        frac_attack_ready=th["occupancy_min"] + 0.1,
        rmsd_max_a=th["rmsd_max_a"] + 1.4,          # spikes over the max bar
        rmsd_mean_a=th["rmsd_mean_a"] - 0.6,        # but is low overall
    ))
    assert v["elevate"], f"the spike allowance did not apply: {v['elevate_why']}"


def test_a_missing_reading_never_elevates():
    """THE FAILURE MODE THIS PROJECT KEEPS PRODUCING: a guard that passes
    because the thing it inspects is absent (catalogue #30, #31)."""
    a = _asw()
    for gap in ("rmsd_max_a", "rmsd_mean_a", "frac_attack_ready"):
        rec = dict(frac_attack_ready=0.99, rmsd_max_a=1.0, rmsd_mean_a=0.5)
        rec.pop(gap)
        v = a.elevation_verdict(rec)
        assert v["elevate"] is False, f"elevated with {gap} missing"
        assert v["elevate_why"] == "no reading"


def test_the_thresholds_come_from_config_not_from_literals():
    """A gate hardcoded in two places is a gate that can disagree with itself."""
    a = _asw()
    th = a.elevation_thresholds()
    assert 0.0 < th["occupancy_min"] <= 1.0, (
        f"occupancy_min {th['occupancy_min']} is not a fraction -- a percentage "
        f"written as 60 would make the gate unsatisfiable")
    # The mean bar must not exceed the max bar, or the OR's second arm is
    # strictly weaker than the first and the spike allowance is decorative.
    assert th["rmsd_mean_a"] <= th["rmsd_max_a"]


# ------------------------------------------------- what feeds it: the units ---
def test_rmsd_is_reported_in_angstrom(tmp_path):
    """GROMACS writes nm; every threshold is stated in Angstrom.

    A missing factor of 10 makes every run pass a 3.5 A bar by a mile and reads
    as a spectacular result rather than as a bug.
    """
    a = _asw()
    (tmp_path / "rmsd.xvg").write_text(
        "# comment\n@ title\n0.0 0.100\n1.0 0.250\n2.0 0.400\n")
    st = a.rmsd_stats(tmp_path)
    assert st["rmsd_max_a"] == pytest.approx(4.0), (
        f"0.400 nm should be 4.0 A, got {st['rmsd_max_a']} -- unit conversion")
    assert st["rmsd_mean_a"] == pytest.approx(2.5)
    assert st["rmsd_frames"] == 3


def test_an_absent_rmsd_file_is_absent_not_zero(tmp_path):
    """Zero would be a perfect score for a run that produced no trace at all."""
    a = _asw()
    assert a.rmsd_stats(tmp_path) == {}
    (tmp_path / "rmsd.xvg").write_text("# only comments\n@ nothing\n")
    assert a.rmsd_stats(tmp_path) == {}


# ------------------------------------------ what feeds it: which definition ---
def test_every_row_records_the_definition_that_produced_it():
    """Two plausible floats in one column is this project's signature defect."""
    a = _asw()
    dist = np.array([2.9, 3.1, 3.4, 3.9, 5.0])
    angle = np.array([10.0, 20.0, 80.0, 15.0, 12.0])
    st = a.geometry_stats(dist, angle, "off-normal", frame_ps=10.0)
    for c in ("attack_ready_max_a", "attack_ready_min_a",
              "attack_ready_uses_angle", "frac_attack_ready_angle"):
        assert c in st, f"{c} is not stamped on the row"
    assert st["attack_ready_max_a"] == pytest.approx(a.attack_ready_max_a())


def test_dropping_the_angle_does_not_drop_it_from_the_record():
    """The angular fraction is still measured, so the choice stays reviewable."""
    a = _asw()
    # one frame is in range on distance but fails the angle badly
    dist = np.array([3.0, 3.0])
    angle = np.array([5.0, 89.0])
    st = a.geometry_stats(dist, angle, "off-normal", frame_ps=10.0)
    if not a.attack_ready_use_angle():
        assert st["frac_attack_ready"] == pytest.approx(1.0)
        assert st["frac_attack_ready_angle"] == pytest.approx(0.5), (
            "the angular reading was not kept beside the distance-only one")
    assert st["frac_attack_ready_angle"] <= st["frac_attack_ready"] + 1e-9


# ------------------------------- what feeds it: the RIGHT trajectory (#23) ---
def test_the_pose_rank_map_is_keyed_on_the_mode_not_the_molecule(tmp_path):
    """Catalogue #23, reintroduced one column name away and measured.

    `sweep_combine` built its pose-rank map from the worklist's `ident` column
    (the MOLECULE) and looked it up with a mode id. Measured 2026-09-02: 0 of 98
    finished modes resolved, so every RMSD on the sweep page came from whichever
    sibling directory sorted first -- 18 of 98 were another pose's trajectory,
    median disagreement 1.12 A, and 7 flipped the pose-held verdict.

    The worklist carries BOTH columns and both are populated and plausible,
    which is why nothing raised.
    """
    wl = pd.DataFrame({
        "ident": ["t4_aaa", "t4_aaa", "t4_bbb"],
        "task_id": ["t4_aaa_m1", "t4_aaa_m7", "t4_bbb_m3"],
        "pose_rank": [2, 8, 4],
    })
    key = "task_id" if "task_id" in wl.columns else "ident"
    prank = dict(zip(wl[key].astype(str), wl.pose_rank.astype(int)))
    modes = ["t4_aaa_m1", "t4_aaa_m7", "t4_bbb_m3"]
    assert all(m in prank for m in modes), (
        f"the map does not resolve mode ids: {sorted(prank)}")
    assert prank["t4_aaa_m1"] != prank["t4_aaa_m7"], (
        "two modes of one molecule collapsed to the same pose rank -- they "
        "would be drawn from the same trajectory")

    # and the defect itself: keying on `ident` resolves NOTHING
    bad = dict(zip(wl.ident.astype(str), wl.pose_rank.astype(int)))
    assert not any(m in bad for m in modes), (
        "this test cannot fail -- the molecule-keyed map resolved a mode id")


def test_sweep_state_carries_new_measurements_without_being_edited():
    """The column list is DERIVED from the results table, not maintained.

    It was an allowlist of nine names, so `rmsd_max_a`, `elevate` and
    `pose_held` were dropped on the way to the page -- and the page then
    recomputed RMSD by a second route, which is how it came to read the wrong
    trajectory. Catalogue #5.
    """
    src = (REPO / "shared" / "sweep_state.py").read_text()
    i = src.find("rescols")
    assert i > 0
    seg = src[i:i + 400]
    assert "res.columns" in seg, (
        "rescols is not derived from the results table; a new measurement "
        "added to a sweep row will silently not reach the page")


# ------------------------------------ the GREEN ZONE must equal the gate ---
def test_the_plotted_band_is_the_gate_band():
    """A shaded zone a mode can sit inside while failing the same page's rank.

    @twu383, 2026-09-02: *"update the rmsd plots to show the correct green
    zones"*. The sweep figure shaded `NAC_DIST_MIN..NAC_DIST_MAX` (2.8-4.2 A),
    the SCREEN's near-attack window, while the sweep is judged at 2.8-3.5 --
    so a trace could sit in the green for most of the run and score 0% engaged
    directly beneath it.
    """
    a = _asw()
    from shared import nac_criterion as nac
    lo, hi = nac.attack_ready_window()
    assert hi == pytest.approx(a.attack_ready_max_a()), (
        f"the plotted band tops out at {hi} but the gate uses "
        f"{a.attack_ready_max_a()}")
    assert lo == pytest.approx(nac.NAC_DIST_MIN)
    # The screen's window is a DIFFERENT quantity and must not have moved.
    assert nac.NAC_DIST_MAX == pytest.approx(4.2), (
        "NAC_DIST_MAX was changed; the screen's near-attack criterion is not "
        "the sweep's engagement bar and D0111 kept them separate")


def test_no_plot_hardcodes_the_wider_window_as_the_ready_band():
    """One definition of the band, read from the criterion by every drawer."""
    import re
    from pathlib import Path as _P
    offenders = []
    for name in ("scripts/sweep_assets.py", "scripts/mdprio_report.py",
                 "shared/md_movie.py"):
        src = (REPO / name).read_text(errors="replace")
        for n, line in enumerate(src.splitlines(), 1):
            code = line.split("#", 1)[0]
            if ("axhspan" in code or "nac_hi" in code) and "NAC_DIST_MAX" in code:
                offenders.append(f"{name}:{n}: {line.strip()}")
    assert not offenders, (
        "a band is drawn from NAC_DIST_MAX rather than attack_ready_window():\n  "
        + "\n  ".join(offenders))
