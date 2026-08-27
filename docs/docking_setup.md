# Docking setup — AutoDock-GPU, gnina, Vina-GPU

*Written 2026-08-27. Every parameter below was read out of the code, the config
or the built receptor on disk, not from memory. Where the code and the config
disagree, that is stated rather than reconciled.*

## 0. There are three docking paths, not one

This is the first thing to understand, because "our docking settings" is not a
single answer and the three paths do not share a receptor, a box, or a scoring
function.

| | engine | used by | mechanism | receptor |
|---|---|---|---|---|
| **A** | **AutoDock-GPU, reactive** | T_3, T_4 — the production screen | covalent, Cys113 | 3IKD prepared, reactive |
| **B** | **gnina** | T_3, T_4 — rescoring + a separate covalent path | covalent | 3IKD plain (rescore) / 6VAJ (protocol default) |
| **C** | **Vina-GPU** | T_1, T_2 | non-covalent | 6VAJ prepared |

Path A generates the poses everything downstream reads. Path B rescores them.
Path C is a different arm of the choreography entirely.

---

## 1. Path A — AutoDock-GPU, reactive docking

The production covalent screen. `scripts/nac_screen.py` builds it,
`scripts/nac_screen_v2.py` drives it.

### 1.1 Binary and environment

| | |
|---|---|
| binary | `/data/lab_vm/modifiable/inhibition/autodock_gpu/bin/autodock_gpu_64wi` |
| version | `89fd1c5e6b4639c22e9a2bea4cc805c42347fffb` |
| prep env | `~/.micromamba/envs/dwi_reactive` (gemmi + meeko + AutoDockTools) |
| receptor source | `/data/lab_vm/modifiable/inhibition/receptor_3ikd_prep/3IKD_noligand.pdb` |
| built receptor | `/data/lab_vm/modifiable/inhibition/receptor_3ikd_reactive/` |

The prep environment matters: `persist_raw_clouds.py` run under `dwi_cheminf`
fails every molecule with `No module named 'gemmi'` (D0095).

### 1.2 Receptor preparation

```bash
mk_prepare_receptor.py --read_pdb 3IKD_noligand.pdb -o rec -p -g \
  -r A:113 --reactive_name CYS:SG \
  --box_center_off_reactive_res --box_size 26 26 26 \
  --r_eq_12 3.2 --eps_12 1.0 \
  --r_eq_13_scaling 1.0 --r_eq_14_scaling 1.0
autogrid4 -p rec_rigid.gpf -l rec_rigid.glg
```

Cached — it is deterministic, and `build_reactive_receptor` returns early if
`rec.reactive_config` and `rec_rigid.maps.fld` both exist.

**The flexible residue is Cys113's side chain, 4 atoms**: `CA`, `CB`, `SG`, `HG`
(`rec_flex.pdbqt`). It moves during docking, which is why `sg_position()` reads
the sulfur out of the *pose* rather than out of the rigid receptor — reading it
from the rigid structure would measure the approach to where the sulfur started.

### 1.3 The reactive parameterisation, and why the published defaults fail

```python
R_EQ_12 = 3.2      # equilibrium distance, reactive atom <-> SG
EPS_12  = 1.0      # well depth
SCALING = 1.0      # 1-3 and 1-4 radii left at full size
```

meeko retypes the ligand's reactive atom to an *order-1 derivative* of its base
type — an aliphatic carbon becomes `C1`, an aromatic one `A1` — and its
neighbours to order 2 and 3. `r_eq_12` is the well for the 1–2 pair. The
`ligand_types` line in `rec.reactive_config` carries the whole derivative
alphabet: `S1…S9`, `C1…C6`, `A1…A3`, `N1…N6`, `H4…H6`, and so on.

**The published covalent defaults (1.8 Å, radii scaled 0.5×) produced 0 of 40
chemically viable poses.** A near-attack conformation is a van der Waals
*contact*, not a bond, so the well sits at contact distance and the neighbouring
atoms keep their real steric radii. That is the whole reason these three numbers
are not the defaults.

### 1.4 The grid, as actually built

Read out of `rec_rigid.gpf`:

```
npts       68 68 68
spacing    0.375
gridcenter 14.036 7.186 -2.108
```

68 × 0.375 = **25.5 Å per side** — AutoGrid rounds the requested 26 Å to an even
point count. The centre is the reactive residue, not a ligand centroid.

> **Divergence worth knowing.** `config/receptor.yaml` gives T_3/T_4 a **20 Å**
> covalent box centred on the QT7 centroid in 6VAJ (−12.61, −34.94, 12.34).
> `nac_screen.py:130` hardcodes **26 Å** centred off the reactive residue in
> 3IKD. Both are populated and plausible; the production screen uses the
> hardcoded one. See §5.

### 1.5 Ligand preparation

`prepare_ligand()`, in order, and every step raises rather than falling through:

1. **Protonate at pH 7.4** via `shared.ionisation.protonate` — the *same* obabel
   call the non-covalent arm makes, reused rather than reimplemented. Before
   D0074 the reactive path docked the SMILES as drawn while the non-covalent path
   protonated, so the two arms docked different species of the same molecule
   (33.3% of T_4, 331 of 5,370 T_3).
2. **Largest fragment** — HTS tables carry salts and counterions.
3. **Re-check the reactive SMARTS still matches.** Warheads are not titratable,
   but an imidazole becoming imidazolium changes ring perception, and a SMARTS
   that quietly stopped matching would read as "no reactive centre" rather than
   as an error.
4. `AddHs` → `EmbedMolecule(randomSeed=0xC0FFEE)` → `MMFFOptimizeMolecule`.
5. **One PDBQT per reactive centre.** `MoleculePreparation(reactive_smarts=...,
   reactive_smarts_idx=0)` returns a setup for *every* match; taking `[0]` picks
   one by default rather than by identity. A fumarate has two genuinely distinct
   electrophilic carbons and needs only one to work, so each is docked and the
   molecule takes the best. (Refusing multiple matches outright previously
   deleted every symmetric Michael acceptor — 4 crystallographic positives.)
6. **Verify reactive typing happened by DIFFING against a plain preparation** —
   never by looking for a named type. meeko derives the type from the base type,
   so hardcoding `"C1"` silently rejected every SNAr ligand (30 negatives and 2
   positives, a whole warhead class).

### 1.6 The docking call

```bash
CUDA_VISIBLE_DEVICES=<gpu> autodock_gpu_64wi \
  -C 1 \
  --import_dpf rec.reactive_config \
  --flexres rec_flex.pdbqt \
  -L <ligand>.pdbqt \
  --nrun 500 \
  --seed 42 \
  --resnam <work>/out
```

| flag | value | note |
|---|---|---|
| `-C` / `--contact_analysis` | `1` | distance-based contact analysis in the DLG |
| `--import_dpf` | `rec.reactive_config` | carries the maps and the reactive types |
| `--flexres` | `rec_flex.pdbqt` | the Cys113 side chain |
| `--nrun` | **500** (`docking.n_runs`) | independent Lamarckian GA runs |
| `--seed` | **42** (`docking.seed`) | `null` ⇒ clock-seeded, for replicate draws |

**Why 500.** `sampling_floor_p: 0.00597` — 500 runs give ≥95% probability of
sampling at least one pose within 2 Å of the true one *provided* the per-run hit
rate exceeds 0.597% (≥99% at 0.917%). Sulfopin's measured rate is 20.8%
[17.4, 24.8], where 13 runs would do. 500 is far past saturation for a molecule
that behaves.

**Why the seed exists (#77).** AutoDock-GPU seeds its GA from the clock unless
told otherwise, so re-screens drew different clouds and produced different
rankings — v4 against v5 agreed at ρ = +0.43 over 504 shared molecules. A seed
makes a run *repeatable*; it does not make the answer *stable* (D0086).

### 1.7 The failure modes, and why `check=True` is not enough

Measured on this build, 2026-08-19:

| `--nrun` | behaviour |
|---|---|
| 5,000 | exit −6, `*** stack smashing detected ***` — **and it still writes a .dlg** |
| 10,000 | exit 0, "The job was not successful", **no .dlg at all** |

The second is the dangerous shape: a zero exit with no result, surfacing later as
a `FileNotFoundError` far from its cause. `dock()` therefore checks **all three**
of return code, DLG presence, and the string `not successful` in the output.
**Depth past ~2,000 runs must come from repeated calls with distinct seeds.**

---

## 2. Path B — gnina

| | |
|---|---|
| binary | `/data/lab_vm/immutable/inhibition/bin/gnina` |
| version | **v1.3.3** `master:6fe1ce2`, built Jun 30 2026 |
| minimum | **≥ 1.3** — where the covalent flags appeared; `assert_version_ok` raises below it |
| runtime env | `/data/lab_vm/envs/dwi_gnina` — exists *solely* to supply `libcudnn.so.9` |

The "static" build is not static. Invoked without `gnina_env()` it dies with
`error while loading shared libraries: libcudnn.so.9`.

### 2.1 Use one — CNN rescoring of mode representatives

`nac_screen_v2.gnina_scores()`, run on the **plain** 3IKD receptor
(`receptor_3ikd/3IKD_prepared.pdbqt`), not the reactive one:

```bash
gnina --receptor <plain>.pdbqt --ligand modes.sdf \
      --score_only --cnn_scoring rescore --seed 42
```

`--score_only` means **no searching** — it scores the poses AutoDock already
produced. Parsed out: `Affinity`, `CNNscore`, `CNNaffinity`, `CNNvariance`.
Timeout 1800 s.

> `cnn_affinity` is **not** the ranking column. Catalogue entry #4: an analysis
> once used it as T_3/T_4's rank metric, and the `LOWER_IS_BETTER` direction
> registry refused the column — the only reason it was caught. The ranking
> column is `affinity_kcal`.

### 2.2 Use two — covalent docking through the pinned protocol

`shared/covalent_protocol.dock()`. Unlike the rescoring path, this one *searches*.

```python
@dataclass(frozen=True)
class DockingParams:
    exhaustiveness: int = 16
    num_modes: int = 9
    cnn_scoring: str = "rescore"
    seed: int = 42
    covalent_optimize_lig: bool = True
    covalent_bond_order: int = 1
```

emitted as:

```bash
gnina -r <receptor>.pdbqt -l <ligand> -o <out>.sdf \
  --center_x/-y/-z ... --size_x/-y/-z ...        # from box.json
  --covalent_rec_atom A:113:SG \
  --covalent_lig_atom_pattern <per-warhead SMARTS> \
  --exhaustiveness 16 --num_modes 9 \
  --cnn_scoring rescore --seed 42 \
  --covalent_bond_order 1 --covalent_optimize_lig
```

The parameter set is *deliberately small* — "every knob here is one more way T_3
and T_4 can diverge without noticing" — and it is frozen, so changing any value
changes the protocol fingerprint recorded with the results.

### 2.3 The SMARTS hazard (D0022, catalogue #21)

`--covalent_lig_atom_pattern` takes **`adduct_attachment_smarts`**, not
`reactive_atom_smarts`. The reactive-atom pattern describes the *unreacted*
warhead and names its leaving group — `[CH2][Cl]`, `[c]([Cl])[n]` — so it cannot
match the post-reaction ligand. Docking through it scored a molecule that does
not exist: the reactive carbon bonded to Cys113 while still carrying its chloride.

The mirror-image error also happened: matching PDB components against the
*adduct* pattern reported **11 novel chemistries**, because the PDB deposits the
**free** ligand and the bond lives in `_struct_conn`. Sulfopin itself came back
"unclassified".

---

## 3. Path C — Vina-GPU, non-covalent (T_1 / T_2)

```python
SEARCH_DEPTH = 20     # D0017 — the adoption evidence does not hold below this
THREADS      = 8000
```

```bash
vina-gpu --receptor 6VAJ_prepared.pdbqt \
  --ligand_directory <in> --output_directory <out> \
  --center_x/-y/-z ... --size_x/-y/-z ...   # box_expanded.json, 26 A
  --thread 8000 --search_depth 20
```

**Why depth 20.** D0017 adopted Vina-GPU only after a sweep showed the depth-10
disagreement with CPU Vina was search *convergence*, not implementation — every
metric improved monotonically with depth, and at 20 the enrichment ROC-AUC drift
is 0.005 against a 0.10 threshold. Pinned, not exposed as a knob: the
non-covalent enrichment gate (ROC-AUC 0.535, D0016) only transfers if the engine
matches.

**The box comes from the receptor object, never from a module constant** — a box
is a set of coordinates in one structure's frame, and carrying it separately is
how a receptor and a box that do not belong together end up on the same command
line, docking into empty space beside the site and returning ordinary-looking
affinities.

**Timeout scales with the pool** (catalogue #19):

```python
SECONDS_PER_LIGAND    = 3.7     # measured, du_xu + guo_pfizer
TIMEOUT_SAFETY_FACTOR = 4
MIN_TIMEOUT_S         = 3600
timeout = max(3600, n_ligands * 3.7 * 4)
```

The old flat `timeout=86400` killed a 16,806-ligand run after 24 h with **0 poses
written** — Vina-GPU writes everything at the *end* of a virtual-screening run,
so a kill loses all of it. The deadline is now logged **before** the run starts.

**Two operational hazards:** the governed wrapper at
`/data/lab_vm/envs/dwi_vinagpu/bin/vina-gpu` **segfaults for anyone who is not
its owner** (it recompiles OpenCL kernels into a directory it cannot write), and
it fails *silently* — a driver reported exit 0 while all six chunks failed,
docking 127 of 30,000, with the previous run's 15,653 pose files still on disk
looking populated and plausible.

---

## 4. What consumes the poses

`shared/nac_criterion.py` — the geometric gate the ranking is actually built on:

```python
NAC_DIST_MIN = 2.8      # A, reactive atom to SG
NAC_DIST_MAX = 4.2
SN2_ANGLE_MIN = 150.0   # degrees, nucleophile-carbon-leaving group
PERPENDICULAR_MAX_OFF_NORMAL = 30.0
APPROACH_WINDOW = (85.0, 125.0)    # Burgi-Dunitz, for addition mechanisms
```

`enrichment = viable_fraction / isotropic_null(mechanism)`. AutoDock's energy is
**deliberately not used to rank molecules** — five measurements say it carries no
signal on this target (D0041, D0046, D0061). It is used only to compare poses
*of one molecule* that have already passed the geometric gate, which is a much
weaker thing to ask of it.

---

## 5. Divergences and hazards — read this part

1. **Two receptors are live.** `covalent_protocol.dock()` defaults to
   `6VAJ_prepared.pdbqt` and `box.json`; the production screen builds and uses
   3IKD reactive. `config/receptor.yaml` still pins `pdb_id: 6VAJ` while
   `config/target.yaml` says `pdb: 3IKD`. **Which receptor you get depends on
   which entry point you came through.**
2. **Two boxes.** Config says 20 Å (covalent, T_3/T_4) and 26 Å (expanded,
   T_1/T_2); `nac_screen.py` hardcodes 26 Å for the covalent screen. The grid on
   disk is 25.5 Å.
3. **The production screen does not persist energies** (D0096). `nac_screen_v2`
   writes coordinates and a mode label; `persist_raw_clouds.py` now writes
   `free_energy_kcal`, but the production path does not. Four clustering
   experiments were energy-blind because of this.
4. **`<topic>_allposes` is not all poses** (D0093) — it holds only poses whose
   DBSCAN label survived, about 79% of the cloud.
5. **`--nrun` above ~2,000 is unsafe** on this AutoDock-GPU build, and one of the
   two failure modes still writes a `.dlg`.
6. **Vina-GPU is all-or-nothing** and segfaults for non-owners.

## Provenance

Every value here was read on 2026-08-27 from: `scripts/nac_screen.py`,
`scripts/nac_screen_v2.py`, `shared/covalent_protocol.py`,
`shared/noncovalent_dock_run.py`, `shared/nac_criterion.py`,
`config/target.yaml`, `config/receptor.yaml`, and the built receptor at
`/data/lab_vm/modifiable/inhibition/receptor_3ikd_reactive/`. Versions were
queried from the binaries.
