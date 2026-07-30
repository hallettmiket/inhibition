#!/usr/bin/env bash
# Purpose: Serve the integration GUI on the loopback interface.
# Author:  Mike Hallett (with Claude Code)
# Date:    2026-07-30
#
# Project members read the RESULTS through this interface rather than the data
# tree. The application runs as the data owner and a viewer needs no permissions
# on the tree at all.
#
# Binds the loopback interface only, and refuses to start if the port is already
# held. Remote viewers connect over an SSH tunnel, the same model murmurent uses
# for its own dashboard. Do not override the address.
#
# The interface is read-only by construction: it never writes to the data tree.
# It also has no per-user authorisation, so it suits panels whose contents are
# meant to be shared with the whole project.
#
# Usage:  bash scripts/serve_gui.sh [PORT]        (default 8899)

set -euo pipefail

PORT="${1:-8899}"
ENV_PY="/data/lab_vm/envs/dwi_gui/bin/python3.11"
APP="$(cd "$(dirname "$0")/.." && pwd)/integration/app/app.py"

[ -x "$ENV_PY" ] || { echo "no interpreter at $ENV_PY" >&2; exit 1; }
[ -f "$APP" ]    || { echo "no app at $APP" >&2; exit 1; }

if ss -ltn 2>/dev/null | grep -q ":${PORT}\b"; then
  echo "port ${PORT} is already in use. Check with:  ss -ltnp | grep ${PORT}" >&2
  exit 1
fi

echo "Serving the integration GUI on the loopback interface, port ${PORT}."
echo "Connect over an SSH tunnel; see the project lead for the recipe."
echo

exec nice -n 19 "$ENV_PY" -m streamlit run "$APP" \
  --server.port "$PORT" \
  --server.address 127.0.0.1 \
  --server.headless true \
  --browser.gatherUsageStats false
