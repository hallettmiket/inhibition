"""Make the GUI coverage gap LOUD instead of a number in a summary line.

`tests/test_app_renders.py` is the only suite that executes every GUI panel, and
`tests/test_pose_modes.py` needs py3Dmol. Both dependencies live in `dwi_gui`,
which has no pytest, so under the suite's own interpreter they `importorskip`
and the run reports `... N skipped`. Two real bugs shipped inside that number
(#45): a panel that claimed to filter and did not, and a selectbox that crashed
when a filter removed the selection.

A skipped test that reads as covered is worse than no test. This prints, at the
end of every run that skipped them, what was NOT executed and the one command
that executes it.

It does not fail the run. Installing pytest into `dwi_gui` is not this repo's to
do -- the env is shared and read-only to most of us -- so the runner script is
the supported route, and the reminder is what keeps it from being forgotten.
"""

from __future__ import annotations

from pathlib import Path

#: Suites whose dependencies live in the GUI env, with what each one covers.
GUI_ENV_SUITES = {
    "test_app_renders.py": "every GUI panel, rendered headlessly (needs streamlit)",
    "test_pose_modes.py": "pose-mode viewer output (needs py3Dmol)",
}

RUNNER = "tests/run_gui_env.sh"

#: The same gap, one env over: scikit-learn lives in `dwi_admet`, which has no
#: pytest, and the suite's own interpreter has pytest and no sklearn. Kept as a
#: separate registry because the message names a different runner -- folding it
#: into GUI_ENV_SUITES would print "GUI COVERAGE" over a clustering test.
SKLEARN_ENV_SUITES = {
    "test_pose_cluster_coords.py":
        "HDBSCAN in 3N coordinate space, and that it agrees with the RMSD-matrix "
        "path (needs scikit-learn)",
}

SKLEARN_RUNNER = ("/data/lab_vm/envs/dwi_admet/bin/python "
                  "tests/run_pose_cluster_env.py")


def pytest_terminal_summary(terminalreporter, exitstatus, config):
    skipped = terminalreporter.stats.get("skipped", [])
    hit = {}
    for rep in skipped:
        name = Path(str(getattr(rep, "nodeid", ""))).name.split("::")[0]
        for suite, what in GUI_ENV_SUITES.items():
            if name == suite:
                hit.setdefault(suite, [what, 0])
                hit[suite][1] += 1
    sk = {}
    for rep in skipped:
        name = Path(str(getattr(rep, "nodeid", ""))).name.split("::")[0]
        if name in SKLEARN_ENV_SUITES:
            sk.setdefault(name, [SKLEARN_ENV_SUITES[name], 0])
            sk[name][1] += 1
    w = terminalreporter.write_line
    if hit:
        w("")
        w("GUI COVERAGE NOT EXECUTED IN THIS RUN (#45)", bold=True)
        for suite, (what, n) in sorted(hit.items()):
            w(f"  {suite}: {n} skipped — {what}")
        w(f"  These do not run under this interpreter. Run them with:  bash {RUNNER}")
    if sk:
        w("")
        w("CLUSTERING COVERAGE NOT EXECUTED IN THIS RUN", bold=True)
        for suite, (what, n) in sorted(sk.items()):
            w(f"  {suite}: {n} skipped — {what}")
        w(f"  These do not run under this interpreter. Run them with:  {SKLEARN_RUNNER}")
