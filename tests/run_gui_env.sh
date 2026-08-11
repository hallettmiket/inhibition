#!/bin/bash
# Run the suites whose dependencies live in the GUI env (#45).
#
# streamlit and py3Dmol live in `dwi_gui`; `dwi_gui` has no pytest and is
# read-only to most of us, so these cannot run under the main suite's
# interpreter and report as skips there. This is the supported way to execute
# them. Two real bugs shipped inside those skips -- see #45.
#
#     bash tests/run_gui_env.sh
#
# Exits non-zero if any panel fails. Cases the AppTest harness cannot drive are
# reported as NOT COVERED by name and do not count as passes.
set -uo pipefail

GUI_PY=${GUI_PY:-/data/lab_vm/envs/dwi_gui/bin/python3}
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [ ! -x "$GUI_PY" ]; then
  echo "no GUI interpreter at $GUI_PY — set GUI_PY to one with streamlit" >&2
  exit 2
fi

rc=0

echo "=== GUI panel render tests ($GUI_PY) ==="
"$GUI_PY" "$REPO/tests/run_app_renders_gui_env.py" || rc=1

# test_pose_modes only needs py3Dmol, which dwi_gui has. It is a pytest file, so
# it runs here only if the GUI env ever gains pytest; until then say so rather
# than let it look covered.
echo
echo "=== pose-mode viewer tests ==="
if "$GUI_PY" -c "import pytest" 2>/dev/null; then
  "$GUI_PY" -m pytest -q "$REPO/tests/test_pose_modes.py" || rc=1
else
  echo "  NOT RUN: $GUI_PY has no pytest, so test_pose_modes.py cannot execute"
  echo "  here either. It needs py3Dmol (present) AND pytest (absent)."
fi

exit $rc
