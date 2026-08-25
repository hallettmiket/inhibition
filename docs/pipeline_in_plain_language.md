# How the ranking pipeline works, in plain language

*Written 2026-08-17. For anyone who needs to understand what this project does
without reading the code. Numbers are from the current run (`nac_v5`, 3.1.0).*

> **The one-sentence version.** We build tens of thousands of candidate
> molecules, use cheap computer simulations to guess which ones might stick to
> our target protein, and narrow that down to a handful worth making in a lab —
> because making one in real life is slow and expensive.
>
> **The honest caveat, up front:** we have never proven that our ranking is
> better than picking at random. Everything below describes what the pipeline
> *does*, not evidence that its answer is right. See §"What we can and cannot
> claim".

---

## The problem in one picture

```
     ~72,000 candidate molecules
              │   (cheap to invent on a computer)
              ▼
       561 in this run's scope
              ▼
         4,432 binding modes
              ▼
         2,019 ranked
              ▼
       170 simulated briefly
              ▼
        a handful simulated properly
              ▼
       1-5 you would actually synthesise      (slow, expensive, real)
```

Each step throws molecules away. The whole design question is whether we are
throwing away the right ones.

---

## The target

**Pin1** is a human protein involved in cancer. We want a molecule that sticks
to it and shuts it off.

Pin1 has a specific spot we aim at: an amino acid called **Cys113**. Cysteine
has a sulfur atom that is chemically "sticky" — under the right conditions it
will form a permanent chemical bond with a molecule that reaches it correctly.
That permanent bond is the goal. A drug that bonds permanently is called a
**covalent inhibitor**.

So every candidate molecule has two parts that matter:

- a **warhead** — the reactive end that is supposed to bond with the sulfur;
- **the rest of the molecule** — which has to hold the warhead in the right
  place, like a handle holding a key to a lock.

---

## Step 1 — Make the candidates

We generate molecules computationally. In this run we use one family of
candidates (called **T_4**), built by combining a warhead with different
chemical decorations.

Three warhead types are in play, chosen because they are the chemistry the lab
can actually synthesise:

| warhead family | molecules |
|---|---|
| acrylamide | 187 |
| bdhi_c4 | 187 |
| bdhi_c5 | 187 |

**561 molecules** go into this run.

**Output:** a list of molecules, as chemical formulas.

---

## Step 2 — Docking: where does each molecule sit?

**Docking** is a simulation that asks: if you drop this molecule into the
protein's pocket, what shape does it curl into, and where does it sit?

The answer is not one shape. The software tries the molecule **500 separate
times**, and each attempt lands somewhere slightly different. Each attempt is
called a **pose** — a specific 3D position and orientation.

So after this step we have 561 molecules × 500 poses = **about 280,000 poses**.

Think of it as dropping a key into a lock 500 times and photographing where it
lands each time.

**Output:** 500 poses per molecule.

---

## Step 3 — Grouping poses into "binding modes"

500 poses is too many to look at, and most of them are near-duplicates. So we
group similar poses together. Each group is a **binding mode** — one distinct
way the molecule sits in the pocket.

This happens in two passes:

1. **First pass** — group poses by where the warhead is and which way it points.
2. **Second pass** — split those groups further if the *rest* of the molecule is
   arranged very differently, up to 5 sub-groups, and only for groups holding at
   least 12 poses.

A molecule usually ends up with a handful of modes. Across the run:
**4,432 binding modes** from 561 molecules — about 8 each.

**Why this matters:** a molecule is not simply "good" or "bad". It might have
one mode that would react beautifully and three that would do nothing. We score
the *modes*, not the molecules.

> ⚠️ **A known weakness.** The first pass uses a grouping method that can chain
> poses together — A is near B, B is near C, so A and C end up in one group even
> if they are far apart. Measured: 22% of first-pass groups contain modes that
> disagree about the key measurement by more than the entire acceptable range.
> Tracked as issue #65.

**Output:** 4,432 modes, each a cluster of poses.

---

## Step 4 — Scoring: could this pose actually react?

Now we ask of each pose: **is the warhead in a position where the chemistry
could happen?**

Two conditions must both hold:

1. **Distance** — the reactive atom must be between **2.8 and 4.2 Å** from the
   sulfur. (An ångström is a ten-billionth of a metre; this is a couple of atom
   widths.) Too far and it cannot reach; too close and the bond has already
   formed.
2. **Angle** — the molecule must be pointing the right way. For this chemistry
   the incoming sulfur has to approach from directly *behind* the part that
   leaves, at **150° or more** out of a possible 180°.

A pose that satisfies both is called **attack-ready**.

Then, for each mode, we count: **what fraction of this mode's poses are
attack-ready?** That fraction is the mode's core score.

Finally we divide by what you would get from a molecule pointing in a completely
random direction. So:

- **score above 1** = better than random
- **score below 1** = worse than random

This score is called **enrichment**.

**Output:** every mode has a score.

---

## Step 5 — Ranking

Modes are sorted best-first, **within their own warhead family**. We do not
compare across families, because the angle requirement is stricter for some
chemistries than others — comparing them directly would punish a chemistry for
being harder to satisfy rather than for being worse.

Two filters apply before a mode is allowed a rank:

- **Enough evidence.** A mode needs at least **12 poses**. A mode with 2 poses,
  both attack-ready, scores the maximum possible — from two observations. That
  is noise wearing a perfect score.
- **Actually measurable.** If no pose ever got close enough to measure the
  angle, the mode has no score. It is left unranked rather than scored zero —
  *not measured* is not the same as *measured and found bad*.

Result: **2,019 of 4,432 modes carry a rank.**

**Output:** a ranked list, viewable in the GUI's *Ranking* page.

---

## Step 6 — Choosing what to simulate

Ranking is cheap. The next step is not, so we pick a shortlist.

Two rules:

- **A score floor** — a mode must score **4.0 or higher** to qualify.
- **A budget cap** — at most **150 modes per warhead family**, so one family
  cannot consume the entire compute budget.

**379 modes** clear the floor; the current worklist holds **170**.

> ⚠️ **The floor has never been validated.** 4.0 was chosen to fit the available
> GPU time, not because anything showed that modes above 4.0 do better than
> modes below it. When we checked, the relationship between this score and the
> outcome it selects for was essentially zero. Tracked as issue #71.

**Output:** a worklist of modes to simulate.

---

## Step 7 — The short simulation ("triage sweep")

Docking is a still photograph. Real molecules move, and water is everywhere.
So we run a proper physics simulation — **molecular dynamics** — where every
atom moves under real forces, in water, at body-like salt concentration.

Each shortlisted mode gets a **5-nanosecond** simulation. That is a very short
time in human terms and a long time in atomic terms.

The question is simple: **does the molecule stay put, or does it drift away?**
We measure how far it wanders from its starting position. If it stays within
**0.35 nm**, it survives.

**Output:** survivors, viewable in the GUI's *Sweep* page.

---

## Step 8 — The long simulation

Survivors get a **100-nanosecond** run — twenty times longer. Same question,
much more demanding: does it *keep* holding on?

**Output:** viewable in the GUI's *MD results* page.

---

## Step 9 — BPMD (the stress test)

Finally, survivors of Step 8 can get a test that actively tries to *push the
molecule out* of the pocket, and measures how hard that is. A molecule that
resists is a better bet than one that merely sat still while nothing disturbed
it.

> This stage is currently not running — it needs its own GPU allocation, and it
> produced 1.3 TB of files last time for no usable result. Tracked as issue #72.

---

## The whole funnel, with current numbers

| step | what it asks | survivors |
|---|---|---:|
| 1. Generate | — | 561 molecules |
| 2. Dock | where does it sit? | ~280,000 poses |
| 3. Group | how many distinct ways? | 4,432 modes |
| 4. Score | could it react? | all scored |
| 5. Rank | which look best? | 2,019 ranked |
| 6. Shortlist | which earn simulation? | 170 |
| 7. Sweep (5 ns) | does it stay put? | in progress |
| 8. Production (100 ns) | does it keep holding? | a few |
| 9. BPMD | how hard to dislodge? | not running |

---

## What we can and cannot claim

This matters more than any number above.

**We can say:** the pipeline produces an ordering, and every step of that
ordering is documented and reproducible.

**We cannot say:** that the molecules at the top of the list will actually work.
Every shortlist this project produces carries the flag `rank_validated = False`,
and that is deliberate and honest.

Three specific reasons:

1. **The scoring step has never beaten a fair test.** When measured on molecules
   whose answers were already known, our scoring did not reliably separate the
   ones that work from the ones that do not.

2. **The one molecule we know works does not pass our own filter.** Sulfopin is
   a published Pin1 inhibitor. We put it through the pipeline disguised as an
   ordinary candidate. The docking *found* its correct shape — 99 of its 455
   poses match the real crystal structure — and the pipeline *kept* that shape.
   But our scoring rule still rates it below random, and it would never reach
   Step 7. No improvement to the docking could fix that; the rule itself is what
   rejects it. (Decision record D0082.)

3. **We have never compared against random selection.** We do not currently know
   whether ranking by our score beats picking 170 modes out of a hat. That is
   the cheapest and most important experiment still outstanding — see
   [`analogous_problems.md`](analogous_problems.md).

None of this means the work is wrong. It means the pipeline is a **candidate
generator that reports its own uncertainty honestly**, rather than a
prediction machine. Given how often computational screens quietly overstate
their confidence, that is a deliberate choice.

---

## Where to look

| you want | go to |
|---|---|
| the live pipeline and its four pages | the GUI on port 8931 |
| what the project currently believes | [`state_of_the_project.md`](state_of_the_project.md) |
| how the numbers go wrong | [`how_this_project_breaks.md`](how_this_project_breaks.md) |
| why each choice was made | [`../decisions/`](../decisions/) |
| whether this is publishable | [`publication_audit.md`](publication_audit.md) |
| every setting named above | [`../config/target.yaml`](../config/target.yaml) |
