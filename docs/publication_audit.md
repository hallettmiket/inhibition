# Publication audit — the pipeline's claims against the literature

*Written 2026-08-12 by @twu383 (with Claude Code), at @tt8804's request: audit
the main claims against existing literature to evaluate publication viability.*

> **Caveat on sourcing, stated up front.** The literature figures below come from
> search-result summaries and abstracts, not from my reading the full papers.
> They are good enough to decide strategy and **not** good enough to put in a
> manuscript. Every number marked ⚠ must be verified against the primary source
> before it is cited. I have flagged where a claim rests on one summary.

---

## The verdict in one paragraph

**The negative-result methods paper, as currently framed in
[`state_of_the_project.md`](state_of_the_project.md) §9, is not viable.** Its
three pillars are: docking enrichment fails here, pose recovery fails here, and
the failure is scoring rather than sampling. Measured against the literature,
the first is *unremarkable* (we are at the field median), the second is
*mis-benchmarked* (we compare cross-docking against a self-docking baseline),
and the third is *a textbook result already quantified for the exact tool we
used*. What survives is narrower, sharper, and still worth publishing — but it
is a different paper, and one of the project's declared dead ends has just been
reopened by someone else's data.

---

## Claim-by-claim

### 1. "Docking does not work on this pocket" — enrichment AUC 0.599, EF1% 0.0 (D0041)

**Literature.** On LIT-PCBA — the realistic, unbiased successor to DUD-E —
methods cluster around a **median ROC-AUC of 0.61–0.62** ⚠. All methods perform
worse on LIT-PCBA than on DUD-E, which is the point of LIT-PCBA.

**Verdict: NOT ANOMALOUS.** Our 0.599 is indistinguishable from the field median
on realistic data. Our own CI, [0.311, 0.874], includes both chance and
excellent — the measurement is underpowered, not decisive. "Docking does not
work on this pocket" is not supportable; "docking works here about as poorly as
it works everywhere" is, and that is not a finding.

**What to do.** Retire this as a headline. Keep it as a characterisation of the
target with the CI stated honestly.

### 2. "Pose recovery 22% top-1 / 41.5% best-of-9, vs a 60–80% norm" (D0046, since_handoff)

**This is the one that will get us rejected, and it is fixable.**

The 60–80% figure is the **self-docking** norm — redocking a ligand into the
structure it was crystallised in. Our 82-case benchmark docks *non-cognate* Pin1
ligands into a single prepared 3IKD. That is **cross-docking**, and the
literature baseline is entirely different ⚠:

| cross-docking setup | top-1 within 2 Å |
|---|---:|
| random receptor selection | 41.1% |
| DUD-E reference receptor | 44.1% |
| largest-pocket receptor | 49.9% |
| best-receptor selection / POSIT | 66.9–68.9% |
| general expectation, structure-based docking | 50–60% |

**Verdict: MIS-BENCHMARKED.** The stated baseline is roughly double the correct
one. Our 22% top-1 is genuinely below the cross-docking norm — that part
survives — but the gap is ~2× rather than ~3×, and our **best-of-9 of 41.5% is
squarely *at* the single-structure cross-docking norm**. A reviewer who knows
this field will catch it in the abstract.

**What to do.** Re-state every recovery number against a cross-docking baseline,
and say explicitly that it is cross-docking into one receptor. This is a writing
fix, not an experiment. **Do it before anything else.**

### 3. "The failure is scoring, not sampling" (D0046, called the most decisive finding)

**Literature.** This is quantified, published, and known specifically for our
tool ⚠: sampling generates a near-native pose in **85–99%** of cases while
scoring ranks it first in only **35–73%**. **AutoDock Vina shows the largest gap
of any program tested — 93.4% sampling against 35–40% ranking.**

**Verdict: KNOWN.** We rediscovered a textbook result about the exact program we
ran. Our own numbers (best-of-9 41.5% vs top-1 18–22%) are the same phenomenon
at lower absolute rates because the setting is harder.

**What to do.** Demote from "central finding" to **positive control**: the
pipeline reproduces a known, quantified effect, which is evidence the
measurement apparatus works. That is a legitimate and useful role for it — but
it is not the paper's contribution.

### 4. "Ranking is partly a size sort, and the direction depends on the pool" (D0043, D0049)

Vina's heavy-atom bias is well documented. The **direction reversal** in the
heavier T_2 pools (ρ = −0.617 in T_1 vs +0.205 in liu) is a nicer observation
than the base effect and I did not find it stated explicitly elsewhere — though
I did not search hard. The size-decorrelation fix (residual within
equal-population strata) is sound engineering, not a novel method.

**Verdict: KNOWN BASE, possibly novel detail.** Worth a figure, not a paper.

### 5. "The covalent stratum is underpowered — 3–4 chemotypes against a floor of 6" (D0045)

**This is the finding that has been overtaken, and it is the most actionable
item in this audit.**

According to PubMed, Xiao et al., *Eur J Med Chem* 2025 — [10.1016/j.ejmech.2025.118048](https://doi.org/10.1016/j.ejmech.2025.118048)
(PMID 40803165) — ran X-ray crystallographic fragment screening on Pin1 and
report:

* **~50 Pin1–fragment complex structures**, deposited (9JYO, 9JZ2, 9JZ4, 9JZ6,
  9JZU, 9KE7, 9KXG, 9KXN, 9KXO and others);
* **two druggable hotspots** — the catalytic centre (Site 1) and **a neighbouring
  site near Cys113 (Site 2)**;
* **a subset of fragments covalently reactive at Cys113**, confirmed by
  crystallography *and* mass spectrometry;
* **enzymatic inhibition assays** on several fragments.

Compare our position: [`state_of_the_project.md`](state_of_the_project.md) §4
declares "more literature searching for a 6th covalent chemotype" ruled out, and
§5 says the missing chemotypes "are in *screening data* only … or must be
commissioned".

**Verdict: SUPERSEDED.** The declared dead end is partly open, and the data is
public. This paper plausibly supplies (a) new covalent Cys113 chemotypes for
D0045's floor, (b) measured actives *and* inactives for the gate that #4 Phase 1
was waiting on, and (c) ~50 structures for the receptor ensemble of Phase 0.3c.

Note we may already have brushed past it: the catalogue records "9KE5 has no
`_struct_conn` section — it is non-covalent", so our curation saw at least one
9K\* entry and classified it non-covalent. **That curation should be re-run over
the whole deposition set**, and if our count is still 3, that is now a claim
about *our SMARTS*, not about the field.

### 6. The near-attack criterion, and that it rejects the positive control (D0065, D0075, D0082)

**Literature.** The gap is acknowledged: covalent docking tools "rank ligands in
their bound (adduct) state and neglect to model the kinetics of covalent
binding, ignoring both the reactivities of the reactants and the orientation of
the prereacted, or intermediate step" ⚠. So *that ground-state pre-reaction
geometry is unmodelled* is known and stated.

What I did **not** find precedent for is the specific measurement in D0082:
building a mechanism-specific pre-reaction geometric criterion, running a known
active blind through the whole pipeline, and showing that **the criterion scores
it below random orientation even at poses within 1.5 Å of its own crystal
structure** (enrichment 0.591 over the cloud → 2.764 at ≤1.5 Å, against a floor
of 4.0).

**Verdict: THE MOST DEFENSIBLE NOVEL ELEMENT.** With one weakness a reviewer
will name: it is a negative result about a criterion *we invented*. A negative
result about a method the field uses is much stronger than one about your own.

### 7. "A multi-agent choreography that measures its own methods honestly"

**Literature.** This space is crowded and moving fast ⚠. Contemporary systems
report *positive, experimentally-validated* discoveries: **Robin** (identified
ripasudil for dry AMD), **OriGene** (validated GPR160 for liver cancer, ARG2 for
colorectal). Plus AgentMol, PiFlow, Mimosa, Mozi, "Beyond SMILES: Evaluating
Agentic Systems for Drug Discovery", **HeurekaBench** (ICLR 2026), and surveys
specifically on collaboration and *failure attribution* in LLM multi-agent
systems.

**Verdict: THE FRAMING HAS BEEN OVERTAKEN.** "We built a multi-agent pipeline"
is no longer a contribution. "Our multi-agent pipeline honestly reported that it
did not work" competes against systems reporting validated hits, and will read
as a null result dressed as a methods claim unless the *self-auditing mechanism*
is the object of study rather than the pipeline.

There is a real contribution buried here — the **defect catalogue** in
[`how_this_project_breaks.md`](how_this_project_breaks.md): 22 instances of one
failure mode, classified into four disguises, with the honest tally of how each
was caught (9 by someone noticing output, 6 while building something else, only
**3 of 22 by a guard**). That is an empirical study of silent-failure modes in
agentic scientific pipelines, and I did not find its equal. It is a
research-practice / MLSys / *Nature Methods*-comment contribution, not a med-chem
one.

### 8. Boltz-2 ruled out for leakage

**Literature.** On a 2,172-complex, 9-warhead covalent benchmark, **Boltz-2 shows
the strongest pose reproduction** among AutoDock4, CovDock, GNINA and Boltz-2 ⚠.

**Verdict: DEFENSIBLE BUT UNDER-ARGUED.** The leakage objection (ChEMBL2288
overlap) is sound for the **affinity** head. It does not obviously transfer to
**pose** prediction, which is what our bottleneck actually is. Reviewers will
ask. Either run it pose-only or state the distinction explicitly.

### 9. Solvent model, synthesizability rules, AiZynthFinder, ensemble docking

Sound, standard, well executed. Supporting methods, not claims. No conflict
found.

---

## What this means for publication

### Option A — the negative-results methods paper, as framed
**Not viable.** Two pillars are known results, one is mis-benchmarked. Reviewers
in this field will recognise all three.

### Option B — reframe on what survives *(recommended, lowest cost)*
A short, honest methods paper on **covalent pose scoring**:

> A mechanism-specific pre-reaction geometric criterion, applied at scale to
> ~72,000 candidates against Pin1 Cys113, cannot rank a known covalent inhibitor
> above random orientation *even when scoring poses within 1.5 Å of that
> inhibitor's own crystal structure* — because the crystallographic reference is
> a post-reaction adduct whose reactive carbon sits below the near-attack
> window, and the reactive geometry is a freely-rotating torsion that a rigid
> docked snapshot samples at close to chance.

That is one sharp, quantitative, reproducible claim, with a blind positive
control and a public codebase. It needs: the recovery numbers re-baselined
against cross-docking (§2), and the pipeline's reproduction of the known
Vina sampling/scoring gap presented as a control (§3).

**Strengths to lead with:** the blind positive control carried through the entire
pipeline is genuinely good practice and rarely done; the pre-registration
documents (`prereg_*.md`); the decision records.

### Option C — convert to a positive-result paper using the new fragment data
Ingest Xiao et al.'s ~50 structures and covalent fragments. That potentially
unblocks the covalent floor (D0045), supplies assayed actives *and* inactives
for the enrichment gate, and gives a real receptor ensemble. Highest ceiling,
highest cost, and it depends on how much of that data is actually usable.

### Option D — split the methods contribution out
The defect catalogue as its own paper, on silent failure modes in agentic
scientific pipelines. Different venue, different audience, and it does not
depend on Pin1 working.

---

## What I would do next, in order

1. **Fix the cross-docking baseline everywhere** (§2). Writing only. Nothing else
   should be submitted until this is done — it is the single most likely cause of
   a desk rejection.
2. **Pull Xiao et al. and its depositions**, re-run `curate_covalent_pdb.py` over
   them, and re-count chemotypes (§5). Cheap, and it decides whether D0045 stands.
3. **Re-verify every ⚠ figure** against the primary source before drafting.
4. **Decide between B and C** on the evidence from step 2.
5. **Boltz-2 pose-only**, as a one-run answer to an objection we will otherwise
   receive (§8).

## Sources

- Xiao et al., *Eur J Med Chem* 299:118048 (2025), [10.1016/j.ejmech.2025.118048](https://doi.org/10.1016/j.ejmech.2025.118048) — Pin1 fragment screening (via PubMed, PMID 40803165)
- [Cross-docking benchmark for automated pose and ranking prediction](https://onlinelibrary.wiley.com/doi/full/10.1002/pro.3784) — Wierbowski et al., *Protein Sci* 2020
- [Benchmarking Cross-Docking Strategies in Kinase Drug Discovery](https://pubmed.ncbi.nlm.nih.gov/39558632/)
- [Efficient conformational sampling and weak scoring in docking programs](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-017-0227-x) — *J Cheminform* 2017
- [Benchmarking single-pose docking and rescoring on LIT-PCBA](https://arxiv.org/pdf/2605.01681)
- [Revealing the limits of covalent docking](https://pubs.rsc.org/en/content/articlelanding/2026/cp/d5cp04981d) — *PCCP* 2026
- [CovDocker: Benchmarking Covalent Drug Design](https://arxiv.org/abs/2506.21085) — KDD 2025
- [Comparative Evaluation of Covalent Docking Tools](https://pubs.acs.org/doi/abs/10.1021/acs.jcim.8b00228) — *JCIM* 2018
- [Re-Evaluating PIN1 as a Therapeutic Target Using Neutral Inhibitors and PROTACs](https://pubs.acs.org/doi/10.1021/acs.jmedchem.4c01412) — *J Med Chem* 2024
- [Discovery of Novel Pyrimidine Derivatives as Human Pin1 Covalent Inhibitors](https://pubs.acs.org/doi/10.1021/acsmedchemlett.4c00477)
- [Beyond SMILES: Evaluating Agentic Systems for Drug Discovery](https://arxiv.org/pdf/2602.10163)
- [Towards Scientific Intelligence: A Survey of LLM-based Scientific Agents](https://arxiv.org/pdf/2503.24047)
