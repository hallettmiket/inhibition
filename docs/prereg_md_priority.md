# Pre-registration — what selects for MD priority?

*Written and committed **before** the six 100 ns runs launch. @tt8804, 2026-08-07.*

---

## The problem

The filter returns **300+ molecules**. A 100 ns MD is ~4 GPU-hours, so running
them all is ~50 GPU-days. We need one cheap measurement that predicts which
molecules will hold their pose, so the expensive step is spent on the ones that
will repay it. **Efficiency without cost to accuracy** — @tt8804.

## What is already ruled out

Nothing derived from docking predicts stability. Measured, not assumed:

| candidate signal | verdict | evidence |
|---|---|---|
| docking energy | **no signal** | ρ(energy_rank, viable) = +0.009 over 115,300 poses (#23) |
| enrichment | **no signal** | D0071, pre-registered cohort |
| consensus | **no signal** | ρ = +0.102 vs BPMD occupancy (n = 29) |
| `topn_viable_frac` | inherits the energy ordering | Sulfopin scores 0.000 (#23) |

## The two survivors, and why neither is proven

**Tier-1 equilibration drift** (300 ps, ~5 GPU-min) separated crystallographic
positives from every candidate group at **p = 0.007**. But that is a *group*
comparison. Whether it ranks *individual* molecules is untested, and the only
per-molecule check available — three finished 100 ns runs — gives drifts of
0.067 and 0.111 nm (held) against 0.120 nm (left). It orders correctly and the
margin is **0.009 nm at n = 3**, which is not evidence of anything.

**BPMD occupancy** (3 ns × replicas, ~1 GPU-h) separated REF from candidates at
p = 0.005–0.021, and is reference-calibrated: the crystallographic median is
**0.163**.

## The prediction, fixed now

Six molecules go to 100 ns from the pose BPMD chose. Their BPMD occupancies span
3.4×:

| molecule | class | pose | BPMD occupancy | vs REF median 0.163 |
|---|---|---:|---:|---|
| `t4_da2e98512d02` | bdhi_c5 | 1 | **0.365** | 2.2× above |
| `t4_7e86b677bb2d` | acrylamide | 6 | 0.189 | above |
| `t4_9a973be6b946` | bdhi_c4 | 2 | 0.161 | at |
| `t4_28f5ea16adeb` | acrylamide | 1 | 0.152 | at |
| `t4_4e608398fd6a` | bdhi_c4 | 1 | 0.125 | below |
| `t4_9265b4bff789` | acrylamide | 8 | 0.108 | below |

**Readout:** residence fraction over 100 ns — the fraction of frames with ligand
RMSD ≤ 1.0 nm — and whether the molecule dissociates at all.

**The prediction:** BPMD occupancy is **rank-correlated with 100 ns residence**.
Specifically, `t4_da2e98512d02` holds and `t4_9265b4bff789` does not.

## Readings, fixed in advance

| observation | conclusion |
|---|---|
| **ρ ≥ +0.7 and the two extremes behave as predicted** | BPMD occupancy is the MD-priority filter. Run it on the 300+, elevate the top by occupancy. ~300 GPU-h, against ~1,200 for elevating everything |
| **ρ positive but weak (+0.3 to +0.7)** | Useful for the extremes only. Elevate the top decile and the bottom decile to keep testing it; do not rank the middle on it |
| **ρ ≈ 0** | BPMD does not predict residence either. Fall back to tier-1 drift on all 300+ (~25 GPU-h) and test *that* per-molecule, since it is 12× cheaper again |
| **ρ negative** | Something is wrong with the protocol or the readout. Report as a failure; do not reinterpret |

**n = 6 supports only large effects.** A null here means "not demonstrated at
n = 6", never "shown to be absent". No p-value will be quoted from six points; the
rank correlation is reported with its scatter and that is all it can carry.

## What this cannot settle

- **All six are T₄ and five of six are BDHI or acrylamide.** A result does not
  transfer to chemistry not represented here.
- **One replicate per molecule at 100 ns.** A single dissociation event is one
  draw with ~100% relative standard error, so residence is a screen not a rate.
- **Poses were chosen by BPMD**, so the 100 ns starts from the pose BPMD liked
  best. That is the intended pipeline, but it means the test is of the pair
  (BPMD-chosen pose + its occupancy), not of occupancy alone.
- **The reference median 0.163 comes from a different cohort** run at the same
  protocol. It is a yardstick, not a control in this experiment.
