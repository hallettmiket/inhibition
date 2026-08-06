#!/bin/bash
# 100 ns MD residence on 3IKD — the chemists' own criterion (#12 §F), whose only
# measurement (D0038/D0044, "not reproducible") was made on 6VAJ and invalidated
# by D0059.
#
# Runs under tmux session `mdres`, following the project's convention
# (degree2, dockd2, repdock, gui).
#
# FAIR USE: two GPUs only, nice -n 19, on a box with ~35 active users. GPUs 0 and
# 4 are deliberately avoided — other people's jobs live there.
set -u
P=$HOME/.micromamba/envs/dwi_reactive/bin/python
cd "$HOME/repos/inhibition" || exit 1

GPUS=(2 3)
N=${#GPUS[@]}
for i in "${!GPUS[@]}"; do
  g=${GPUS[$i]}
  echo "$(date -Is) shard $i -> GPU $g"
  nice -n 19 timeout 259200 "$P" scripts/md_residence_3ikd.py \
    --limit 0 --n-neg 5 --production-ps 100000 --nrun 200 \
    --gpu "$g" --shard "$i" --n-shards "$N" --keep \
    > "/tmp/mdres_s$i.log" 2>&1 &
done
wait
echo "$(date -Is) all shards complete"
