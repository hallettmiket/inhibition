#!/usr/bin/env bash
#
# Purpose: Build Vina-GPU 2.1 (OpenCL) against the conda Boost + NVIDIA OpenCL.
# Author:  Mike Hallett (with Claude Code)
# Date:    2026-07-27
# Input:   the Vina-GPU-2.1 checkout under $ENV_ROOT/_src
# Output:  AutoDock-Vina-GPU-2-1 binary, plus a wrapper that sets the runtime env
#
# WHY THIS IS NOT JUST `make`:
#
# 1. Upstream's Makefile COMPILES Boost thread sources directly
#    ($BOOST/libs/thread/src/pthread/thread.cpp), which needs a Boost SOURCE
#    tree. Conda ships headers plus prebuilt shared libraries, so those sources
#    are dropped and we link -lboost_thread instead.
# 2. Upstream ships `MACRO=-DAMD_PLATFORM` while also setting
#    GPU_PLATFORM=-DNVIDIA_PLATFORM. On an NVIDIA box the AMD macro must go.
# 3. The README recommends -DOPENCL_3_0 on Linux and warns off OPENCL_1_2.
# 5. -DSMALL_BOX is required, not optional: `option` in wrapcl.cpp is declared
#    ONLY inside the LARGE_BOX/SMALL_BOX ifdefs, so omitting both fails to
#    compile. SMALL_BOX is upstream's default and covers boxes up to 100 A;
#    ours are 20 and 26 A.
# 4. Conda's OpenCL loader does not read /etc/OpenCL/vendors, so the NVIDIA ICD
#    is invisible and clinfo reports ZERO platforms. OCL_ICD_VENDORS must point
#    at a vendors dir — without it the binary builds fine and then finds no GPU.
#
# Vina-GPU is a DIFFERENT ENGINE from CPU Vina. It must be validated against
# already-scored ligands before it replaces anything (see validate_vina_gpu.py):
# swapping the docking engine changes every number downstream.

set -euo pipefail

ENV_ROOT="${DWI_ENV_ROOT:-/data/lab_vm/envs}"
P="$ENV_ROOT/dwi_vinagpu"
SRC="$ENV_ROOT/_src/Vina-GPU-2.1/AutoDock-Vina-GPU-2.1"
ICD_DIR="$P/etc/OpenCL/vendors"

_log() { echo "[$(date -Is)] $*"; }

[ -d "$SRC" ] || { echo "Vina-GPU source not found at $SRC" >&2; exit 1; }
[ -f "$P/include/boost/version.hpp" ] || { echo "boost not in $P" >&2; exit 1; }

# --- OpenCL ICD so the loader can see the NVIDIA platform ------------------
mkdir -p "$ICD_DIR"
if [ ! -s "$ICD_DIR/nvidia.icd" ]; then
  if [ -f /etc/OpenCL/vendors/nvidia.icd ]; then
    cp /etc/OpenCL/vendors/nvidia.icd "$ICD_DIR/nvidia.icd"
  else
    echo "/lib/x86_64-linux-gnu/libnvidia-opencl.so.1" > "$ICD_DIR/nvidia.icd"
  fi
fi
export OCL_ICD_VENDORS="$ICD_DIR"
export LD_LIBRARY_PATH="$P/lib:${LD_LIBRARY_PATH:-}"

_log "OpenCL platforms visible:"
"$P/bin/clinfo" 2>/dev/null | grep -m1 "Number of platforms" || true

# --- build ----------------------------------------------------------------
cd "$SRC"
_log "building Vina-GPU 2.1"
nice -n 19 g++ -o AutoDock-Vina-GPU-2-1 \
  -I"$P/include" -I"$P/include/boost" \
  -I"$SRC/lib" -I"$SRC/OpenCL/inc" \
  ./main/main.cpp -O3 \
  ./lib/*.cpp ./OpenCL/src/wrapcl.cpp \
  -lboost_program_options -lboost_system -lboost_filesystem -lboost_thread \
  -lOpenCL -lstdc++ -lstdc++fs -lm -lpthread \
  -L"$P/lib" \
  -DOPENCL_3_0 -DNVIDIA_PLATFORM -DSMALL_BOX -DNDEBUG \
  -DBOOST_TIMER_ENABLE_DEPRECATED \
  -DBUILD_KERNEL_FROM_SOURCE \
  2>&1 | tail -25

[ -x ./AutoDock-Vina-GPU-2-1 ] || { echo "build produced no binary" >&2; exit 1; }
_log "binary built: $SRC/AutoDock-Vina-GPU-2-1"

# --- wrapper so callers need not know about the ICD ------------------------
cat > "$P/bin/vina-gpu" <<WRAP
#!/usr/bin/env bash
# Vina-GPU 2.1 with the OpenCL ICD and library path preset.
export OCL_ICD_VENDORS="$ICD_DIR"
export LD_LIBRARY_PATH="$P/lib:\${LD_LIBRARY_PATH:-}"
cd "$SRC"
exec nice -n 19 "$SRC/AutoDock-Vina-GPU-2-1" "\$@"
WRAP
chmod +x "$P/bin/vina-gpu"
_log "wrapper installed: $P/bin/vina-gpu"
