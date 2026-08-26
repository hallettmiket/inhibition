#!/usr/bin/env bash
# Launch the pose-group viewer. Port is overridable: PORT=8932 bash <this>
set -euo pipefail
cd "$(dirname "$0")/.."
exec /data/lab_vm/envs/dwi_gui/bin/streamlit run integration/pose_group_viewer.py \
  --server.port "${PORT:-8932}" --server.address 127.0.0.1 \
  --server.headless true --browser.gatherUsageStats false
