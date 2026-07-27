# What are "decoys" and why do we use them?

*Plain-language explainer for the docking-enrichment gate, Dance with Inhibition (Pin1).*

## 1. The core idea

Decoys are the **negative controls** of docking. We already have a small set
of molecules we know bind Pin1 (the "actives" — our ~6 known binders). If we
only ever score those, we can't tell whether the docking program is smart or
just lucky — every method gives *some* score to *something*. So we hand the
docking program a big pile of molecules we don't think bind ("decoys") and
check whether it can tell them apart from the real binders.

Think of it like testing a metal detector: you don't just wave it over a
coin and declare victory. You bury the coin in a box of sand, gravel, and
other junk, and check that the detector beeps loudest over the coin and
stays quiet over everything else. The sand and gravel are the decoys.

## 2. Why decoys have to be "property-matched" but structurally different

This is the part that's easy to get wrong in two opposite directions.

**Decoys must look similar on the outside** (roughly the same size, the same
greasiness/oiliness, the same charge/polarity) **but be different on the
inside** (a different chemical scaffold, not just a tweak of the active
molecule).

- **If decoys are too different physically** — e.g. all our decoys are small
  and polar while the actives are big and greasy — then the docking score
  doesn't have to know anything about Pin1's binding pocket to "win." It can
  just learn the trivial rule "bigger, greasier molecules score better,"
  which is true of almost any docking program on almost any protein, real
  binder or not. The test would pass for the wrong reason and we'd have
  learned nothing about whether docking works *on this receptor*.

- **If decoys are too structurally similar to the actives** — e.g. we pick
  decoys that are just the actives with one atom changed — some of those
  decoys are probably real, undiscovered binders too. We would have
  accidentally "salted" our negative control with hidden positives, making
  the test artificially *harder* to pass, and any failure would be
  ambiguous (is docking bad, or did we just poison our own control set?).

The property-matching is what makes it a fair test: same "outside"
statistics as the actives, different "inside" chemistry. This is the same
logic behind decoy sets used in the field generally (e.g. the DUD, DUD-E,
and DEKOIS decoy-generation approaches), though our decoys here are
generated from ChEMBL for this specific project rather than pulled from one
of those published sets.

## 3. What the enrichment metrics mean, in plain terms

Once actives and decoys are both docked and scored, we ask: did the actives
end up ranked near the top, or scattered randomly among the decoys?

- **ROC-AUC** — the probability that, if you pick one active and one decoy
  at random, the docking score ranks the active better. 0.5 = coin flip
  (docking has no discriminating power on this target). 1.0 = perfect
  separation.
- **Enrichment factor (EF)** — look only at the top slice of the ranked list
  (say, the top 5%) and ask: how many actives showed up there compared to
  what you'd expect by chance? An EF of 10 means "10x better than random" in
  that top slice.
- **BEDROC** — like enrichment factor, but it specifically rewards actives
  that land at the very top of the list, and cares less about the rest.
  This matters because in practice we can only afford to synthesize and
  test a handful of top-ranked compounds — a method that gets 3 of our 6
  actives into the top 20 is far more useful to us than one that gets all 6
  actives scattered somewhere in the top half.

## 4. What it means if this gate fails

If docking can't separate our known Pin1 binders from the property-matched
decoys, that's evidence docking scores aren't trustworthy *for ranking
candidates on this receptor*. In that case, we would demote the docking
score from something we rank candidates by to something we merely *display*
alongside a candidate (for context, not for prioritization). The
downstream approaches would then lean more heavily on their other lines of
evidence — plus our own scientific judgement — rather than trusting the
docking number to pick winners.

## 5. The specific snag we just hit: covalent decoys need a "warhead"

Pin1 inhibitors can work in two modes: **non-covalent** (binds and lets go)
and **covalent** (binds and forms a permanent chemical bond to the protein,
usually to a specific cysteine). Docking each mode requires a different
kind of simulation.

For **covalent docking specifically**, a molecule has to physically be
*capable* of forming that bond — it needs a reactive chemical group called a
"warhead" (e.g. an acrylamide or similar electrophile) that can react with
the target cysteine. A molecule with no warhead simply cannot be scored by
covalent docking at all — the program has nothing to attach.

We just discovered that only about **11% of the decoys we generated carry a
warhead**. That means roughly 9 out of 10 of our covalent decoys can't even
be *run* through the covalent docking pipeline, let alone serve as a
negative control for it. This breaks the control the same way a metal
detector test breaks if 90% of your "junk" box is empty air instead of
sand and gravel — there's nothing there to fail to detect, so a clean-looking
result wouldn't actually tell us the detector works.

**Options going forward:**
- **Select (or generate more) decoys that do carry a warhead**, so the
  covalent gate has a real, size-matched pool of warhead-bearing
  non-binders to test against.
- **Treat the covalent gate differently** — e.g. evaluate it only on the
  warhead-bearing subset we do have (smaller, but honest), and be explicit
  that the covalent enrichment numbers are lower-confidence than the
  non-covalent ones until we've built out a proper covalent decoy set.

Either way, the honest read right now is: our non-covalent enrichment
numbers stand on reasonably solid ground, but our covalent enrichment
numbers are not yet backed by an adequate negative control and shouldn't be
over-interpreted until we fix the decoy set.
