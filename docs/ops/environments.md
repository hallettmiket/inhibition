# Environments

Five isolated conda envs under `/data/lab_vm/envs`, **not** in `$HOME` — `/home`
is at 89% and three of these carry separate torch+CUDA builds.

```bash
bash scripts/setup_envs.sh all        # or one at a time
conda activate /data/lab_vm/envs/dwi_cheminf
```

| Env | Contents | Status |
|---|---|---|
| `dwi_cheminf` | RDKit, CReM, Meeko, Vina, PoseBusters, xtb, OpenBabel | ✓ verified |
| `dwi_diffsbdd` | torch 2.5.1+cu121, PyG, DiffSBDD deps | ✓ 8 GPUs, scatter/cluster OK |
| `dwi_reinvent4` | REINVENT 4 | rebuilding |
| `dwi_amber_md` | AmberTools, GROMACS, OpenMM, gmx_MMPBSA | ⚠ wheel-build failure |
| `dwi_gui` | Streamlit, py3Dmol, stmol, MkDocs | ✓ |
| `dwi_admet` | ADMET-AI | ✓ built |
| gnina | **binary/Docker, no env** | not yet pinned |

!!! danger "Never merge the torch consumers"
    DiffSBDD, REINVENT 4 and ADMET-AI pull mutually incompatible torch/CUDA
    pins. Co-installing them is the single most common way this stack breaks.
    `torch-scatter`/`torch-cluster` against the wrong CUDA is DiffSBDD's #1
    install failure — use the matching PyG wheel index.

## Verified on biodatsci

RDKit 2025.09.5 · CReM 0.3.1 · Vina 1.2.7 · Meeko 0.7.1 · PoseBusters 0.6.5 ·
torch 2.5.1+cu121 · 8× A100-80GB · CUDA 12.0 · 224 cores.

SAscore resolves from RDKit Contrib at
`$ENV/share/RDKit/Contrib/SA_Score` — all four approaches report that axis, so
if it ever goes missing, vendor `sascorer.py` rather than substituting a
different metric.

## Known issues

- **`dwi_amber_md`** built but with a wheel-build failure inside it, likely
  `gmx_MMPBSA`. Needs checking before T_4 step 9 or T_2 step 11.
- **REINVENT 4** upstream dropped `requirements-linux-64.lock`, which the
  implementation plan named. Now `pyproject` + `uv.lock` + `install.py`.
