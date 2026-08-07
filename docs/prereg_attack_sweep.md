# Pre-registration — the attack-geometry gate, and what BPMD actually measures

*Written before the sweep is built or run. @tt8804, 2026-08-07. Raised in
[#32](https://github.com/hallettmiket/inhibition/issues/32).*

---

## Why

100 ns residence and attack geometry are close to independent:

| molecule | residence | frames attack-ready |
|---|---:|---:|
| `t4_4e608398fd6a` | **1.000** | **0.4%** |
| `t4_7e86b677bb2d` | **1.000** | **1.7%** |
| `t4_da2e98512d02` | 0.791 | **55.2%** |
| `t4_9a973be6b946` | 0.774 | 7.3% |
| `t4_72f5671e89cb` | 0.541 | 0.7% |
| `t4_9265b4bff789` | 0.383 | 0.0% |

ρ(residence, attack-ready) = **+0.319**, *p* = 0.538. **The two molecules with
perfect residence essentially never present the warhead; the molecule that
dissociated at 81 ns is attack-ready more than half the time.** Ranking on
residence selects against the thing we want.

@tt8804: *"100% residence can still result in a warhead that does not stay near
the attack distance or angle... we should do 10 ns sweeps to first filter for
that ability... before comitting to 100 ns MD."*

## Three observations already in hand, and what they are worth

These come from truncating six existing trajectories. They motivate the design;
they are **not** the test, because the design was chosen after seeing them.

1. **A 10 ns window orders molecules like the full 100 ns** for window occupancy
   (ρ = +0.829, *p* = 0.042) and more weakly for attack-readiness (ρ = +0.600,
   *p* = 0.208).
2. **BPMD occupancy predicts attack-readiness (ρ = +0.900, *p* = 0.037) far
   better than it predicts residence (ρ = +0.410, *p* = 0.493)** — and the
   original pre-registration graded it on residence, where it looks like a null.
   Mechanistically unsurprising in hindsight: BPMD's collective variable **is**
   the warhead→SG distance.
3. **Only 1 of 6 elevated poses starts attack-ready**, and that one is the only
   molecule that got anywhere. Starting geometry costs nothing to check.

## The tests, fixed now

Run on a **fresh cohort of 8–10 molecules** drawn from the 2.2.0 mode-based
ranking, each taken through all four stages so the cheap ones can be scored
against the expensive one.

### T1 — does the 10 ns sweep predict the 100 ns attack geometry? *(primary)*

**Readout:** Spearman ρ between the sweep's `frac_attack_ready` and the 100 ns
`frac_attack_ready`, on molecules that were **not** used to design this.

### T2 — does BPMD beat the sweep at a fraction of the cost?

Both are proxies for the same quantity. BPMD is ~1 GPU-h, the sweep ~0.4 GPU-h,
the truth ~4 GPU-h. **Readout:** ρ for each against the 100 ns attack-readiness,
on the same molecules.

### T3 — is starting geometry a free filter?

**Readout:** fraction of molecules whose 100 ns attack-readiness exceeds 5%,
split by whether the elevated pose was attack-ready at frame 0. Costs no GPU at
all — it is a property of the pose.

### T4 — occupancy or visits?

`frac_attack_ready` measures how *long* a molecule is competent; `n_visits`
counts independent excursions into attack geometry. A covalent reaction needs
**one** good approach, not sustained occupancy, so visits is the more
mechanistically honest observable. On the six existing trajectories the two rank
identically (ρ = +1.000), which means the choice is currently free — and that
may not survive a cohort with more range.

## Readings, fixed in advance

| observation | conclusion |
|---|---|
| **T1 ρ ≥ +0.7** | Adopt the sweep as the gate before 100 ns. ~10× cheaper per rejection |
| **T1 ρ +0.3 to +0.7** | Use it to reject the bottom only, never to rank the middle |
| **T1 ρ ≈ 0** | The sweep filters on noise. Drop it; the early window is dominated by the starting pose |
| **T2: BPMD ρ > sweep ρ** | BPMD moves EARLIER in the funnel and the sweep is dropped — one proxy, the better one |
| **T2: sweep ρ > BPMD ρ** | BPMD retires. It would then predict neither residence nor attack better than a cheaper unbiased run |
| **T2: both ≈ equal** | Keep the cheaper one (sweep) and retire BPMD on cost |
| **T3 separates cleanly** | Gate on starting geometry FIRST, for free, before any simulation |
| **T4 diverge, visits wins** | Switch the observable to visits and re-derive every attack number on it |

## What this cannot settle

- **n = 8–10.** Supports large effects only. A null means "not demonstrated at
  this n", never "absent".
- **Every molecule starts from a pose chosen for its anchoring geometry**, so the
  early window is biased toward looking competent and a sweep could pass
  everything. On the six existing runs it did not — sweep values spanned 0.000 to
  0.832 — but that is the failure mode to watch, and T1's ρ is the thing that
  detects it.
- **Attack geometry is a proxy for reactivity, not a measurement of it.** Nothing
  here says a molecule reacts; it says the warhead reaches a geometry from which
  reaction is possible. That is discharged at FEP, on the covalent leg.
- **The BPMD result (ρ = +0.900) is a lead on five points**, found while checking
  something else. It is being tested here precisely because it was not predicted.

## The original MD-priority pre-registration stands

`docs/prereg_md_priority.md` fixed BPMD occupancy against **100 ns residence**,
and it will be reported that way when the sixth molecule lands — as written, and
as a null if that is what it shows. **This document does not re-read it.** The
attack-geometry hypothesis is a separate claim, tested on separate molecules,
with its readings fixed above before any of them are run.
