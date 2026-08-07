# From ranked modes to one pose worth 4 GPU-hours

*@tt8804's proposal, worked through and measured. 2026-08-07.*

> *"we can have a probabilistic pose based on the poses used to generate each
> mode and try to validate that with maybe alphafold 3 or Boltz-2 and then do
> bpmd (or maybe dont even need bpmd) and then use that high confidence pose for
> non-cov 100 ns MD, if it fails the transition from ranked pose modes to single
> consensus pose we move down the list and dont waste md"*

The instinct is right and it points at the correct bottleneck. The measurements
change **which** tool does **which** job, and they kill one part of the proposal
outright.

---

## 1. What is already solved, and what actually isn't

| step | status |
|---|---|
| finding the right pose at all | **solved** — present in 500 poses for 100% of molecules (#30) |
| picking the right **mode** | **solved** — mode population picks it 93.3%, which *is* the ceiling |
| picking the right **pose within a mode** | **NOT solved** — best in-house rule is 26.7% |

So the bottleneck is exactly where you put it: the **mode → single pose**
transition. Everything upstream of it is working.

For completeness, the in-house rules for picking one pose out of the consensus
mode, scored against the crystal:

| rule | within 2 Å |
|---|---:|
| ceiling — best pose in the mode | 93.3% |
| medoid (most typical) | 26.7% |
| lowest energy in the mode | 20.0% |
| best `anchor_quality` | 6.7% |

---

## 2. "Probabilistic pose" cannot mean averaging — measured

Averaging the Cartesian coordinates of a mode's poses **destroys the molecule**:

| | median nearest-atom distance |
|---|---:|
| a real docked pose | 1.52 Å (6VAJ), 1.40 Å (7EKV) |
| the mode's coordinate average | **0.43 Å**, **0.68 Å** |

A C–C bond is ~1.5 Å. The average of many orientations of a flexible molecule
collapses toward its own centroid and is not a structure at all — it cannot be
parameterised, simulated, or drawn. Any "consensus pose" has to be **an actual
member of the ensemble**, or a real structure predicted independently.

What the distribution *can* legitimately give is **uncertainty**, not a
structure: `consensus` (population), `spread_a` (positional tightness),
`dir_coherence` (orientational agreement, 0–1). Those are the confidence numbers
your gate needs.

---

## 3. The measurement that reorganises the design

Boltz-2 was tested two ways on the 15 crystal complexes.

**As a mode ARBITER — it fails.** It does not pick the right mode better than
population does:

| arbiter | all 15 | multi-mode only (n=9) |
|---|---:|---:|
| most populated mode | **93.3%** | **88.9%** |
| closest to Boltz-2 | 86.7% | 77.8% |
| consensus × Boltz agreement | 93.3% | 88.9% |

Consensus is already at the ceiling, so there is nothing for an arbiter to add.
**Drop this leg.**

**As a POSE SOURCE — it succeeds, and it solves the open problem:**

| | |
|---|---:|
| Boltz-2's prediction lands **inside** the consensus mode | **100%** (15/15), median 0.89 Å |
| Boltz-2's pose within 2 Å of the crystal | **67%** |
| our medoid within 2 Å of the crystal | 27% |

The two methods agree completely about **where** — and Boltz-2 is **2.5× better
at the pose within that region**, which is precisely the step our own rules fail.

---

## 4. The design

```
  500 docked poses                                   ~4 GPU-h, whole library
        ↓  pose splitting (density on reactive atom + warhead direction)
  modes → candidate rows, ranked by consensus × anchoring
        ↓  collapse to distinct parent molecules → shortlist
  for each mode, in rank order:                      ~50 s/molecule, SHORTLIST ONLY
        Boltz-2 prediction
        ├─ falls inside this mode  →  ADOPT BOLTZ'S POSE, elevate it
        └─ falls outside           →  next mode down; no MD spent
        ↓
  [BPMD — kept only if tonight's verdict earns it]
        ↓
  100 ns non-covalent MD: residence + attack radius   ~4 GPU-h each
        ↓
  FEP (covalent) on the shortlist
        ↓
  synthesis
```

**Each tool does the thing it was measured to be good at.** Consensus picks the
mode; Boltz-2 picks the pose; MD decides.

### Why Boltz-2 sits on the shortlist, not the library

~50 s per molecule with the MSA cached. Across 5,765 molecules that is **~80
GPU-hours** — more than the docking it would be validating. On a shortlist of
~100 it is **~1.4 GPU-hours**. It is cheap *relative to MD*, not cheap absolutely,
and placing it before the ranking would invert the funnel.

### The gate, and its honest status

*"if it fails the transition… we move down the list and don't waste md"* — the
gate is **Boltz-2's prediction falling inside the mode**, with `dir_coherence`
and `spread_a` as secondary confidence.

**It never fired on these 15 — 100% landed inside.** So it costs nothing here,
and its false-positive rate is **unmeasured**. Our generated molecules are larger
and floppier than this benchmark (median 28 heavy atoms against 22), and pose
difficulty scales exponentially with size, so it should fire more often on real
candidates. Adopt it as a **cheap tripwire**, not as a validated filter, and
record every time it triggers so its rate can be measured on the molecules that
matter.

### BPMD

Tonight's pre-registered verdict decides it. BPMD asks whether a pose survives a
biased dynamics perturbation — a physics question that Boltz-2's pattern-matching
cannot answer, so they are not substitutes. But BPMD is a ~1 GPU-hour *proxy* for
a 4 GPU-hour measurement, and a proxy that does not predict its target has no
role. **If occupancy does not rank-correlate with 100 ns residence, drop it** and
go straight from the Boltz-2-confirmed pose to MD.

---

## 5. What this still does not solve

- **The crystal is the wrong target for reaction competence.** The deposited
  ligand is a post-reaction adduct with its leaving group gone, so RMSD-to-crystal
  leaves the leaving-group direction unconstrained — and the SN2 angle depends
  entirely on it. Every "within 2 Å" number here measures *where the molecule
  sits*, never *whether it can react*. MD is the only arbiter of the latter that
  we have.
- **Boltz-2 predicts the non-covalent complex** — which is in register with the
  whole pre-FEP pipeline, not a deficit, but it means the gate tests placement
  and not reactivity.
- **Mode count reproduces only 73%** across independent dockings (dominant mode
  86.7%). Ranking the dominant mode is supported; treating every minor mode as an
  independent candidate is not yet.
- **n = 15**, all crystallographic, median 22 heavy atoms. Everything here is
  extrapolated to a library that is bigger and harder.
