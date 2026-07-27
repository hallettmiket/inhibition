#!/usr/bin/env bash
#
# Purpose: Build the five isolated conda envs for Dance with Inhibition.
# Author:  Mike Hallett (with Claude Code)
# Date:    2026-07-27
# Input:   none (channels + pins are declared inline)
# Output:  envs under $ENV_ROOT, one per tool family; a log per env
#
# WHY FIVE ENVS: three of these consume torch (DiffSBDD, REINVENT 4, ADMET-AI)
# with mutually incompatible CUDA/PyG pins. Co-installing them is the single
# most common way this stack breaks. gnina is a binary/Docker call, so it needs
# no env at all and can never conflict.
#
# Usage:
#   bash scripts/setup_envs.sh cheminf      # one env
#   bash scripts/setup_envs.sh all          # everything (long; run under tmux)

set -euo pipefail

ENV_ROOT="${DWI_ENV_ROOT:-/data/lab_vm/envs}"
LOG_DIR="${DWI_LOG_DIR:-/data/lab_vm/append_only/inhibition/_env_logs}"
mkdir -p "$ENV_ROOT" "$LOG_DIR"

# Solving a torch/CUDA env with the classic solver takes hours; libmamba is the
# difference between a coffee and an afternoon.
CONDA_SOLVER_ARGS=(--yes --solver=libmamba)

_log() { echo "[$(date -Is)] $*"; }

build_cheminf() {
  # The shared CPU workhorse. Every approach imports RDKit from here, so there
  # is exactly ONE RDKit in the choreography and descriptor values cannot drift
  # between approaches (which is what makes the cross-approach physchem axes
  # comparable at all).
  local p="$ENV_ROOT/dwi_cheminf"
  _log "building cheminf -> $p"
  conda create --prefix "$p" "${CONDA_SOLVER_ARGS[@]}" \
    -c conda-forge \
    python=3.11 rdkit pandas pyarrow numpy scipy pyyaml \
    openbabel xtb click tqdm
  # CReM's fragment DB is a separate ~GB download (see fetch_fragment_db.sh).
  "$p/bin/pip" install --no-input crem meeko vina posebusters
}

build_diffsbdd() {
  # T_1's generator. torch-scatter/torch-cluster against the wrong CUDA is
  # DiffSBDD's #1 install failure — the PyG wheel index must match the torch
  # build, so these are installed from the pinned index, not plain PyPI.
  local p="$ENV_ROOT/dwi_diffsbdd"
  _log "building diffsbdd -> $p"
  conda create --prefix "$p" "${CONDA_SOLVER_ARGS[@]}" \
    -c conda-forge python=3.10 openbabel rdkit numpy scipy pyyaml biopython
  "$p/bin/pip" install --no-input torch --index-url https://download.pytorch.org/whl/cu121
  local tv
  tv="$("$p/bin/python" -c 'import torch; print(torch.__version__.split("+")[0])')"
  "$p/bin/pip" install --no-input torch-scatter torch-cluster \
    -f "https://data.pyg.org/whl/torch-${tv}+cu121.html"
  "$p/bin/pip" install --no-input pytorch-lightning wandb
}

build_reinvent4() {
  # T_3's generator. Three things the implementation plan got wrong, all because
  # upstream moved: (1) requirements-linux-64.lock no longer exists; (2) REINVENT
  # now needs python >= 3.11, not 3.10; (3) install.py's FIRST positional is the
  # PyTorch processor type (cu124/cpu/mac), NOT the optional-dependency set —
  # passing "none" there yields a bogus --extra-index-url .../whl/none.
  #
  # CRITICAL: install.py shells out to a BARE `pip`. Without this env's bin
  # first on PATH it resolves to whatever pip comes first — which on this
  # machine was base conda's, installing REINVENT and ~30 chemistry packages
  # into the base environment and upgrading pandas and pydantic out from under
  # streamlit, anndata and anaconda-cloud-auth.
  local p="$ENV_ROOT/dwi_reinvent4"
  local src="$ENV_ROOT/_src/REINVENT4"
  _log "building reinvent4 -> $p"
  conda create --prefix "$p" "${CONDA_SOLVER_ARGS[@]}" -c conda-forge python=3.11
  mkdir -p "$ENV_ROOT/_src"
  [ -d "$src" ] || git clone https://github.com/MolecularAI/REINVENT4.git "$src"
  ( cd "$src" && PATH="$p/bin:$PATH" "$p/bin/python" install.py cu124 -d none )
  # Not declared by REINVENT's own metadata, but its CLI imports it.
  "$p/bin/pip" install --no-input scipy
  "$p/bin/reinvent" --help >/dev/null && _log "reinvent CLI OK"
}

build_amber_md() {
  # T_4 step 9 + T_2 step 11. gmx_MMPBSA drags an old AmberTools pin along with
  # GROMACS and mpi4py; this env exists to quarantine that.
  local p="$ENV_ROOT/dwi_amber_md"
  _log "building amber_md -> $p"
  conda create --prefix "$p" "${CONDA_SOLVER_ARGS[@]}" \
    -c conda-forge \
    python=3.10 ambertools=23 gromacs openmm mpi4py numpy pandas
  "$p/bin/pip" install --no-input gmx_MMPBSA || \
    _log "WARN: gmx_MMPBSA pip install failed — may need ambertools=21; see plan 2.x"
}

build_gui() {
  # The artist's Streamlit app. Read-only over the four Di_top10.csv, so it
  # stays light and never needs a docking or MD dependency.
  local p="$ENV_ROOT/dwi_gui"
  _log "building gui -> $p"
  conda create --prefix "$p" "${CONDA_SOLVER_ARGS[@]}" \
    -c conda-forge python=3.11 rdkit pandas pyarrow numpy matplotlib
  "$p/bin/pip" install --no-input streamlit py3Dmol stmol
}

build_admet() {
  # Kept separate: chemprop pulls its own torch, which must not land in cheminf.
  local p="$ENV_ROOT/dwi_admet"
  _log "building admet -> $p"
  conda create --prefix "$p" "${CONDA_SOLVER_ARGS[@]}" \
    -c conda-forge python=3.11 rdkit pandas numpy
  "$p/bin/pip" install --no-input admet-ai
}

main() {
  local what="${1:-}"
  case "$what" in
    cheminf)   build_cheminf   2>&1 | tee "$LOG_DIR/cheminf.log" ;;
    diffsbdd)  build_diffsbdd  2>&1 | tee "$LOG_DIR/diffsbdd.log" ;;
    reinvent4) build_reinvent4 2>&1 | tee "$LOG_DIR/reinvent4.log" ;;
    amber_md)  build_amber_md  2>&1 | tee "$LOG_DIR/amber_md.log" ;;
    gui)       build_gui       2>&1 | tee "$LOG_DIR/gui.log" ;;
    admet)     build_admet     2>&1 | tee "$LOG_DIR/admet.log" ;;
    all)
      for e in cheminf gui admet diffsbdd reinvent4 amber_md; do
        main "$e"
      done ;;
    *)
      echo "usage: $0 {cheminf|diffsbdd|reinvent4|amber_md|gui|admet|all}" >&2
      exit 2 ;;
  esac
  _log "done: $what"
}

main "$@"
