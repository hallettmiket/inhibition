# The shape of our problem, and the fields that already solved it

*Written 2026-08-17 by @twu383 (with Claude Code), at @tt8804's request: "we
generate a huge space of derivatives but to actually pick and commit to one to
test in real life is very expensive, so we rank and screen — what are some well
known analogous problems in CS or other fields?"*

> **Sourcing caveat, stated up front.** Unlike
> [`publication_audit.md`](publication_audit.md), this note was **not** produced
> by literature search. It is background knowledge, written to give #66 and #71
> a shared vocabulary. The named results (CAST, the Prentice criterion, the NAS
> random-baseline papers, Hyperband) are real and well known, but **verify any
> of them against the primary source before citing it in a manuscript.** The
> value here is the framing, not the references.

---

## The shape, stated without chemistry

> A cheap generator produces an enormous candidate set. The **true objective**
> can only be evaluated a handful of times, at high cost. So a **cheap proxy**
> is used to rank and filter — and the proxy's relationship to the true
> objective is **assumed rather than measured**.

Every difficulty this project has had is a consequence of the last clause, not
of the first two. We are good at generation (~72,000 candidates) and the
expensive step is genuinely expensive (synthesis + assay). What we do not have
is a demonstrated relationship between what we rank on and what we care about.
That is not a chemistry problem, and several fields have a mature methodology
for it.

**This reframing matters for publication.** "Our scorer did not work" is a null
result about one tool. "**We ran an unvalidated surrogate endpoint and measured
its validity**" is a methods contribution with an established literature to sit
in. Same experiments, different claim. See [`publication_audit.md`](publication_audit.md) §B.

---

## 1. Surrogate endpoints in clinical trials — the closest analogue

The most rigorous treatment of exactly our problem, because people died of
getting it wrong.

**The canonical failure — CAST (Cardiac Arrhythmia Suppression Trial, 1989).**
Antiarrhythmic drugs suppressed ventricular arrhythmia — the surrogate —
convincingly, and **increased mortality**. The surrogate moved the right way
while the endpoint moved the wrong way. Others in the same family: bone mineral
density vs fracture, CD4 count vs AIDS progression, tumour response rate vs
survival.

**The methodology it produced — the Prentice criterion.** A surrogate is valid
only if the effect on the surrogate *captures* the effect on the true endpoint,
and **this is something you measure, not something you argue**. Validation means
running both endpoints on a designed subset and reporting the proportion of the
effect explained.

**What it says about us.** `rank_validated = False` is, in this vocabulary, an
**unvalidated surrogate endpoint** — and that is a far more precise statement
than "the ranking is not validated". #71 is a textbook instance:
`sweep_rule.floor` selects on `enrichment`, and the measured relationship
between `enrichment` and the outcome it selects for is
**Spearman +0.016, p = 0.84**. The field's answer is not "find a better
surrogate"; it is "measure the surrogate against the endpoint on a stratified
subset", which is exactly what the pilot in #71 is. We already designed the
right experiment; we have not run it.

## 2. Neural architecture search — the closest CS analogue, including the crisis

NAS has our problem *and* our reckoning.

Enormous architecture space; the true objective (train to convergence and
measure accuracy) is expensive; so the field built cheap proxies — low-fidelity
training, learning-curve extrapolation, zero-cost proxies.

**The reckoning.** A run of papers around 2019–2021 showed that **random search
was a strong baseline**, that several proxies barely beat it, and that published
comparisons were not reproducible because baselines were weak and search spaces
did the work. "NAS evaluation is frustratingly hard" is a real title.

**What the field did next is the instructive part.** It did not abandon NAS. It
adopted **mandatory random baselines** and shared benchmarks (NAS-Bench-style
tabular datasets where every architecture's true score is precomputed, so a
proxy can be scored against ground truth cheaply).

**What it says about us.** Two things:

1. **We owe a random baseline.** See §6 — this is the recommendation.
2. Our 82-case redocking benchmark *is* our NAS-Bench: a set where the true
   answer is known, so a proxy can be scored cheaply. We should use it that way
   more aggressively than we do — D0082 is exactly that move, applied to one
   molecule.

## 3. Cascade ranking — the design rule we have backwards

Our pipeline is literally a **cascade ranker**: cheap candidate generation
(docking), then progressively more expensive re-rankers (8 ns sweep, 100 ns
production, BPMD). Web search and object detection (Viola–Jones) have used this
architecture for decades, and they settled a design rule:

> **An early stage is tuned for RECALL AT THE CUT, not precision at the top.**

The asymmetry is structural: a true positive lost at stage 1 is **unrecoverable**
— no downstream stage can retrieve it. A false positive merely costs one
expensive evaluation. So early stages are run at near-100% recall and are
*allowed* to be imprecise.

**What it says about us.** `sweep_rule.budget_floor = 4.0` discards **3,700 of
4,432 modes**, and it was chosen to fit the GPU budget rather than to retain
productive modes. That is an early stage tuned for precision. Worse, D0082
measured that the one molecule we know works — sulfopin — **does not clear it**,
at any achievable pose quality. In cascade terms we have a stage-1 filter with a
demonstrated false-negative on the only ground truth we have.

Note the config already encodes the right instinct and does not act on it:
`sweep_rule.capture_target: 0.95` says the floor must retain 95% of productive
modes, with the comment that missing one costs a candidate while a wasted sweep
costs ~0.5 GPU-h. That *is* the recall-at-cut rule. It is written down and
unenforced, because the floor it governs has never been measured.

## 4. Best-arm identification — the framework for spending the budget

The formal treatment of "many options, expensive evaluations, one fixed budget"
is **pure-exploration multi-armed bandits** / best-arm identification, and the
practical algorithms are **successive halving** and **Hyperband**: evaluate
everything cheaply, keep the top fraction, spend more on the survivors, repeat.

The theory is about how to split a budget between **breadth** (how many
candidates get a cheap look) and **depth** (how much each survivor gets), which
is precisely the `max_depth` vs `budget_floor` tension in `config/target.yaml` —
currently resolved by hand and by GPU arithmetic.

**What it says about us.** Our cascade is already successive halving with
hand-set fractions. The literature would have us set those fractions from the
measured *variance* of the cheap stage rather than from the GPU budget — and
would predict that a stage whose signal is near zero (see §1) should be given
**less** depth and more breadth, not the reverse.

## 5. Goodhart's law, and model exploitation in RL

The generic name. "When a measure becomes a target, it ceases to be a good
measure." In model-based reinforcement learning the specific failure is **model
exploitation**: the policy discovers and exploits errors in the learned model,
so performance *against the model* climbs while true performance falls.

**What it says about us.** This is the principled justification for a decision we
already made on instinct. [`state_of_the_project.md`](state_of_the_project.md) §4
rules out genetic algorithms and REINVENT RL against the current objective
because "they optimise the oracle harder, and the oracle is broken … it would
look like progress". That is model exploitation, it has a name and a literature,
and citing it converts a preference into a principle. D0043 is our own measured
instance: a model generating larger molecules scores better on our ranking
without binding better.

---

## 6. What I would actually do: the random baseline

Of everything above, one action is cheap, decisive, and currently missing.

**We have 147 modes selected by `enrichment ≥ 4.0` and swept. We have no
comparison against 147 modes selected at random.**

Without it we cannot answer the first question any reviewer of #66 will ask —
*does your ranking beat picking at random?* — and every other claim in the
project is downstream of that answer. It is the discipline NAS adopted after its
own crisis, and it is a much smaller experiment than #71's full stratified pilot
(which it complements rather than replaces: the random arm establishes whether
there is signal at all; the stratified pilot locates where the floor should sit).

Note the asymmetry in outcomes, which is what makes it worth running:

* if the ranking beats random, we have the project's first validated claim;
* if it does not, that is a stronger and more publishable negative result than
  anything currently in the repo — and it is the honest version of what §2's
  papers did to NAS.

---

## Vocabulary, for the manuscript

| what we say now | what the literature calls it |
|---|---|
| "the ranking is not validated" | an **unvalidated surrogate endpoint** (Prentice) |
| "the floor was chosen for GPU budget" | an early cascade stage tuned for **precision, not recall at the cut** |
| "optimising harder would look like progress" | **model exploitation** / Goodhart's law |
| "the criterion rejects the positive control" | a **demonstrated false negative against ground truth** |
| "we don't know if the ranking helps" | **no random baseline** |
| the 82-case redocking benchmark | a **tabular benchmark** for cheap proxy scoring (NAS-Bench-style) |

## Related

- [`publication_audit.md`](publication_audit.md) — what survives for publication;
  §B's reframing is this document's §1.
- **#71** — the enrichment floor has never been measured. §1 and §3.
- **#66** — publication audit. §2 and §6.
- **D0082** — the criterion rejects sulfopin at crystallographic poses. §3.
- **D0043** — larger molecules score better without binding better. §5.
- [`how_this_project_breaks.md`](how_this_project_breaks.md) — the failure
  catalogue. Orthogonal to this document: that one is about defects in what we
  compute, this one is about whether the thing we compute is the right thing.
