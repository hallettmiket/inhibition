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
2. ~~**BPMD occupancy predicts attack-readiness (ρ = +0.900, *p* = 0.037) far
   better than it predicts residence (ρ = +0.410, *p* = 0.493)**~~ — and the
   original pre-registration graded it on residence, where it looks like a null.
   Mechanistically unsurprising in hindsight: BPMD's collective variable **is**
   the warhead→SG distance.

   > **WITHDRAWN 2026-08-11 (#35). Do not quote the +0.900.** The two sides of
   > that correlation were computed on **different starting poses**. BPMD ran on
   > 2–3 poses per molecule and kept the best; the 100 ns run always used pose 1.
   > For 4 of 7 molecules those are not the same pose, and for two of them pose 1
   > was never put through BPMD at all. The number compares *best-of-N poses*
   > against *one different pose* — not what the sentence claims.
   >
   > #35 asks for it to be re-derived on matched poses. **That is not possible
   > from what is on disk.** Of the five molecules named, three have no surviving
   > elevation rows at all, and the remaining two
   > (`t4_72f5671e89cb`, `t4_9a973be6b946`) carry `warhead_pose_idx` — the
   > reactive-atom index — not the pose rank, so which pose each BPMD replicate
   > ran on cannot be recovered. A matched re-derivation needs new compute, not
   > new analysis.
   >
   > A separate defect in the same number, from #35: best-of is taken over three
   > tries for most molecules and fewer for at least one, and fewer tries
   > systematically gives a lower best-of.
   >
   > **Nothing downstream may inherit this figure**, including the downstream-BPMD
   > pre-registration #51 asks for. If BPMD-vs-attack-readiness is to be claimed,
   > it is a fresh measurement on matched poses at equal replicate count.
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

## The countervailing risk — @tt8804, and it changed the build

> *"be careful that we are not prioritizing attack geometry over realistic poses.
> we are going through a lot of trouble to translate consensus poses to more
> supported ones using much larger tools than gnina."*

The concern is that gating on attack geometry selects poses that *point the
warhead well* rather than poses that are *where the molecule actually sits* —
undoing the work Boltz-2 and mode consensus do.

**Tested, and the answer is more specific than the worry.** Across the whole pose
population, higher `anchor_quality` correlates with being **closer** to both
references — median ρ = **−0.143** against Boltz-2's independent prediction and
**−0.135** against the crystal, over 15 molecules (negative = closer). Anchoring
is not systematically selecting unrealistic poses.

**But its ARGMAX is.** Picking one pose out of the dominant mode:

| rule | within 2 Å of crystal |
|---|---:|
| ceiling — best pose in the mode | 93.3% |
| **medoid of the top-25% by anchoring** | **33.3%** |
| medoid of the whole mode | 26.7% |
| **argmax anchoring** | **6.7%** |

The maximum of a noisy score is an outlier, typically a strained pose that
happens to present the warhead. **Narrow on anchoring, then take a typical member
of what survives** — that beats either alone, and it is now what the screen does.

`n = 15`, so 33.3% against 26.7% is one molecule and the quartile width is
untuned. What is *not* one molecule is 6.7% against 26.7%: **argmax is the thing
to stop doing**, and it was what I had built.

**The ordering the concern implies is now the architecture.** Realism first,
attack geometry second, and never the reverse:

```
mode consensus          where the molecule actually sits        (93.3% at ceiling)
   ↓
Boltz-2 confirmation    independent structural support          (67% vs 27%)
   ↓
PoseBusters            physical validity, still unused
   ↓
attack geometry        RANKS what survives; never selects it
```

Attack geometry is the last term applied, never the first, and it ranks within a
set already established as realistic. **T5 below tests that this ordering is the
right one rather than assuming it.**

### T5 — does the realism-first ordering beat attack-geometry-first?

On the fresh cohort, rank the same molecules two ways: (a) realism gates then
attack ranking, (b) attack ranking alone. **Readout:** which ordering better
predicts 100 ns attack-readiness, and whether (b) elevates poses that PoseBusters
or Boltz-2 disagree with.

| observation | conclusion |
|---|---|
| (a) ≥ (b) | ordering confirmed; keep realism as a gate, attack as a rank |
| (b) > (a) | the realism gates are costing us; re-examine what they reject |
| (b) elevates PoseBusters failures | attack-geometry-first is selecting artefacts, exactly as @tt8804 warned |

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
- ~~**The BPMD result (ρ = +0.900) is a lead on five points**, found while checking
  something else.~~ **Withdrawn 2026-08-11 (#35)** — the two sides were computed on
  different starting poses and it cannot be re-derived on matched ones from data
  on disk. See the note under observation 2. It is not a lead; it is a number
  that was never measuring what its sentence said.

## The original MD-priority pre-registration stands

`docs/prereg_md_priority.md` fixed BPMD occupancy against **100 ns residence**,
and it will be reported that way when the sixth molecule lands — as written, and
as a null if that is what it shows. **This document does not re-read it.** The
attack-geometry hypothesis is a separate claim, tested on separate molecules,
with its readings fixed above before any of them are run.
