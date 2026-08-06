#!/bin/bash
# Consensus + enrichment for every candidate, from ONE docking run each.
# Fair use: 4 GPUs, nice -n 19. Avoids 0 (other users), 2 (my BPMD run) and 4
# (someone holding memory). Checked with nvidia-smi before scaling, not assumed.
set -u
P=$HOME/.micromamba/envs/dwi_reactive/bin/python
cd "$HOME/repos/inhibition" || exit 1
GPUS=(1 3 5 6)
for i in "${!GPUS[@]}"; do
  g=${GPUS[$i]}
  echo "$(date -Is) shard $i -> GPU $g"
  nice -n 19 timeout 86400 "$P" scripts/nac_consensus_all.py \
    --shard "$i" --n-shards "${#GPUS[@]}" --nrun 200 --top-n 10 \
    --gpu "$g" --chunk 100 > "/tmp/consall_s$i.log" 2>&1 &
done
wait
echo "$(date -Is) consensus-for-all complete"
