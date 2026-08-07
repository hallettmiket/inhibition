#!/usr/bin/env bash
# Purpose: launch the v2 geometric screen across four GPUs, under fair use.
# Author: Mike Hallett (with Claude Code)
# Date: 2026-08-06
#
# Persists per-pose geometry + gnina scores + the poses themselves, so the two
# known score repairs can be tested without docking a third time.
#
#   scripts/nac_v2_launch.sh
#
# GPUs 0, 4 and 7 are refused by name: 0 and 4 are carrying other people's jobs
# and 7 is left free. nice -n 19 on every worker; this box has ~35 users.
set -euo pipefail
GPUS=(1 2 3 5)
REPO=/home/UWO/twu383/repos/inhibition
PY=~/.micromamba/envs/dwi_reactive/bin/python
LOGS=/data/lab_vm/modifiable/inhibition/nac_v2_logs
for g in "${GPUS[@]}"; do
  if [ "$g" -eq 0 ] || [ "$g" -eq 4 ] || [ "$g" -eq 7 ]; then
    echo "GPU $g is spoken for; refusing" >&2; exit 1
  fi
done
mkdir -p "$LOGS"
tmux has-session -t nacv2 2>/dev/null || tmux new-session -d -s nacv2 -n idle
for s in 0 1 2 3; do
  win="s$s"; log="$LOGS/s$s.log"
  cmd="cd $REPO && nice -n 19 $PY scripts/nac_screen_v2.py --shard $s --n-shards 4 \
       --gpu ${GPUS[$s]} --chunk 50 2>&1 | tee $log"
  tmux kill-window -t "nacv2:$win" 2>/dev/null || true
  tmux new-window -t nacv2 -n "$win" "bash -lc '$cmd; exec bash'"
  echo "shard $s -> GPU ${GPUS[$s]}, window nacv2:$win"
done
echo; echo "watch: tmux attach -t nacv2"
