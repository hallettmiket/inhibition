#!/usr/bin/env bash
# Purpose: chain the whole pipeline overnight — wait for the screen, rank, select, elevate, report.
# Author: Mike Hallett (with Claude Code)
# Date: 2026-08-06
#
#   tmux new-session -d -s overnight 'bash scripts/overnight.sh'
#
# WHY A CHAIN RATHER THAN FOUR MANUAL STEPS. @tt8804 is away and asked that the
# ranking be finished by morning at minimum. Each stage here is cheap except the
# screen and the MD, so the whole thing is gated on the screen finishing and then
# runs unattended. Every stage logs to $LOGS and writes its own timestamped
# output, so a stage that fails leaves the earlier ones intact.
#
# ORDER IS THE PRIORITY ORDER @tt8804 GAVE:
#   1. ranking finished           <- must land
#   2. GUI carrying the new list
#   3. selection + elevation reports on the top few
#
# Fair use: GPUs 0, 4 and 7 are never touched. nice -n 19 throughout.
set -uo pipefail

REPO=/home/UWO/twu383/repos/inhibition
PY=~/.micromamba/envs/dwi_reactive/bin/python
PYC=/data/lab_vm/envs/dwi_cheminf/bin/python
LOGS=/data/lab_vm/modifiable/inhibition/overnight_logs
GPUS=(1 2 3 5)
MDGPUS=(1 2 3 5)
mkdir -p "$LOGS"
cd "$REPO" || exit 1

say() { echo "$(date -Is) | $*" | tee -a "$LOGS/chain.log"; }

# ---------------------------------------------------------------- 1. the screen
say "waiting for the v2 screen (4 shards)"
while pgrep -u "$USER" -f "nac_screen_v2.py" >/dev/null 2>&1; do
  n=$($PYC - <<'EOF' 2>/dev/null || echo 0
import glob, pandas as pd
fs = glob.glob('/data/lab_vm/append_only/inhibition/00_outputs/blacksmith/nac_v2/agg_s*_*.csv')
d = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True).drop_duplicates('ident')
print(len(d))
EOF
)
  say "screen at ${n}/5769"
  sleep 600
done
say "screen finished"

# ------------------------------------------------------------------ 2. ranking
say "ranking (topn_viable_frac, filtered on consensus_gnina)"
nice -n 19 $PY scripts/rank_v2.py --score topn_viable_frac --top 15 \
  > "$LOGS/rank.log" 2>&1
say "ranking exit $?  -> $LOGS/rank.log"

# also emit the alternative orderings so the morning has the comparison
for sc in enrichment_conditional enrichment_joint; do
  nice -n 19 $PY scripts/rank_v2.py --score "$sc" --top 10 \
    > "$LOGS/rank_$sc.log" 2>&1
done
say "alternative orderings written"

# ------------------------------------------------------------------ 3. the GUI
say "refreshing the GUI's ranked list"
nice -n 19 $PY scripts/gui_refresh_v2.py > "$LOGS/gui.log" 2>&1
say "gui exit $?  -> $LOGS/gui.log"

# ---------------------------------------------------------------- 4. selection
say "automatic selection"
for tier in T4 T3; do
  nice -n 19 $PY scripts/select_elevate.py --tier "$tier" --per-class 2 \
    --require-geometry > "$LOGS/select_$tier.log" 2>&1
  say "select $tier exit $? -> $LOGS/select_$tier.log"
done

# ------------------------------------------------- 4b. pose ranking by BPMD
# Ranking molecules and ranking a molecule's own poses are DIFFERENT problems.
# Stage 2 chose the molecules; this chooses which of a molecule's poses is real
# enough to spend 4 GPU-hours on. ~1 GPU-hour to protect a 4 GPU-hour run.
say "ranking poses within each selected molecule (BPMD)"
# poses written before the writer stamped pose_rank cannot be addressed by rank
nice -n 19 $PY scripts/backfill_pose_rank.py > "$LOGS/backfill.log" 2>&1
nice -n 19 $PY scripts/rank_poses_bpmd.py --max-poses 3 --replicates 2 \
  --production-ps 3000 --gpu "${GPUS[0]}" > "$LOGS/poserank.log" 2>&1
say "pose ranking exit $? -> $LOGS/poserank.log"
WINNERS=$(ls -t /data/lab_vm/append_only/inhibition/00_outputs/blacksmith/pose_rank_bpmd/pose_rank_*.csv 2>/dev/null | head -1)
say "winners: ${WINNERS:-none}"

# ---------------------------------------------------------------- 5. elevation
# 100 ns is ~4 GPU-hours per molecule, so only as many as there are cards, and
# ONE replicate. That is a screen, not a residence measurement -- a single
# dissociation event carries ~100% relative standard error (#22).
say "launching elevation on the queue"
nice -n 19 $PY scripts/elevate_queue.py --n "${#MDGPUS[@]}" \
  --gpus "${MDGPUS[@]}" --production-ps 100000 \
  ${WINNERS:+--winners "$WINNERS"} \
  > "$LOGS/elevate.log" 2>&1
say "elevation launcher exit $? -> $LOGS/elevate.log"

say "chain complete; MD continues in tmux session 'elevate100'"
