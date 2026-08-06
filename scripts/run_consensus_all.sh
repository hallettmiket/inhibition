#!/bin/bash
# Consensus + enrichment for every candidate, from ONE docking run each.
# Fair use: 2 GPUs, nice -n 19, avoiding 0/2/4 (other jobs + the BPMD run).
set -u
P=$HOME/.micromamba/envs/dwi_reactive/bin/python
cd "$HOME/repos/inhibition" || exit 1
GPUS=(5 6)
for i in "${!GPUS[@]}"; do
  g=${GPUS[$i]}
  echo "$(date -Is) shard $i -> GPU $g"
  nice -n 19 timeout 86400 "$P" scripts/nac_consensus_all.py \
    --shard "$i" --n-shards "${#GPUS[@]}" --nrun 200 --top-n 10 \
    --gpu "$g" --chunk 100 > "/tmp/consall_s$i.log" 2>&1 &
done
wait
echo "$(date -Is) consensus-for-all complete"
