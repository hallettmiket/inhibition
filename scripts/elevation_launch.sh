#!/usr/bin/env bash
# Purpose: launch one tier of the elevation experiment across four GPUs, under fair use.
# Author: Mike Hallett (with Claude Code)
# Date: 2026-08-06
#
# FAIR USE IS ENFORCED HERE RATHER THAN REMEMBERED. GPUs 0 and 7 carry other
# people's jobs and are refused by name; the allowance is four GPUs and the
# script takes exactly the four it is given. `nice -n 19` on every worker.
#
#   scripts/elevation_launch.sh 1          # tier 1: equilibration survival
#   scripts/elevation_launch.sh 2 3000     # tier 2: BPMD, 3 replicas x 3 ns
#
# Each shard runs in its own window of the `elevate` tmux session and tees to
# /data/lab_vm/modifiable/inhibition/elevation_logs/.
set -euo pipefail

TIER="${1:?usage: elevation_launch.sh <tier 1|2> [production_ps]}"
PROD_PS="${2:-3000}"
GPUS=(1 2 3 5)
REPO=/home/UWO/twu383/repos/inhibition
PY=~/.micromamba/envs/dwi_reactive/bin/python
LOGS=/data/lab_vm/modifiable/inhibition/elevation_logs
REPLICATES=3

for g in "${GPUS[@]}"; do
  if [ "$g" -eq 0 ] || [ "$g" -eq 7 ]; then
    echo "GPU $g belongs to someone else's job; refusing" >&2; exit 1
  fi
done

mkdir -p "$LOGS"
tmux has-session -t elevate 2>/dev/null || tmux new-session -d -s elevate -n idle

for s in 0 1 2 3; do
  win="t${TIER}s${s}"
  log="$LOGS/t${TIER}_s${s}.log"
  cmd="cd $REPO && nice -n 19 $PY scripts/elevation_run.py --tier $TIER \
       --shard $s --n-shards 4 --gpu ${GPUS[$s]} --replicates $REPLICATES \
       --production-ps $PROD_PS 2>&1 | tee $log"
  tmux kill-window -t "elevate:$win" 2>/dev/null || true
  tmux new-window -t elevate -n "$win" "bash -lc '$cmd; exec bash'"
  echo "shard $s -> GPU ${GPUS[$s]}, window elevate:$win, log $log"
done

echo
echo "tier $TIER launched: 4 shards, $REPLICATES replicas each, ${PROD_PS} ps production"
echo "watch:  tmux attach -t elevate"
