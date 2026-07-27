---
title: Literature review — CPU AutoDock Vina vs GPU-accelerated reimplementations, score/ranking agreement
date: 2026-07-27
project: inhibition
sensitivity: standard
tags: [docking, vina-gpu, autodock-vina, gpu, benchmarking, literature-review]
sources: ['@bookworm']
---

# Does the literature document score/ranking discrepancies between CPU AutoDock Vina and GPU reimplementations?

**Bottom line up front:** yes — every Vina-GPU-family paper that reports a
head-to-head comparison against CPU AutoDock Vina finds a small, systematic,
GPU-less-negative score offset in the same direction and roughly the same
magnitude you measured (~0.2 kcal/mol, Pearson r in the 0.87–0.98 range
depending on benchmark). None of the primary Vina-GPU papers attribute this
to floating-point precision, despite using single-precision (FP32) GPU
kernels — they attribute it to the different (massively-parallel,
subspace-partitioned) search algorithm. **Direct CPU-vs-GPU enrichment/AUC
comparisons on an actives+decoys benchmark — the thing you actually did — are
essentially absent from the primary literature.** The enrichment numbers the
Vina-GPU 2.1 paper advertises are GPU-2.1-vs-GPU-2.0, not GPU-vs-CPU. That is
a real gap, and your finding (CPU AUC 0.535 → GPU AUC 0.433, both consistent
with chance, CIs overlapping) sits in a space nobody else seems to have
published data on.

---

## 1. Do the Vina-GPU papers report score agreement with CPU Vina?

### Vina-GPU (original), Tang et al., *Molecules* 2022
- Ding J, Tang S, et al. is **not** this one — the original is Tang S, Chen R,
  et al., "Accelerating AutoDock Vina with GPUs," *Molecules* 27(9):3041
  (2022). DOI: `10.3390/molecules27093041`. PMCID: `PMC9103882`.
- **140-complex redocking benchmark** (their curated set, described as
  "well pre-processed according to the original literature" — this is the
  same 140-complex set used by the AutoDock-GPU paper, not explicitly
  labeled PDBbind/CASF in the text I could retrieve):
  - Average docking score: CPU AutoDock Vina **−8.9** kcal/mol vs Vina-GPU
    **−8.7** kcal/mol → **mean bias ≈ +0.2 kcal/mol, GPU less negative**,
    matching your +0.188 kcal/mol bias almost exactly, both in sign and
    magnitude.
  - Pearson **r = 0.965**.
  - RMSD success (<2 Å from crystal pose): CPU **114/140**, GPU **107/140**;
    average RMSD 1.5 Å (CPU) vs 1.7 Å (GPU).
  - The paper explicitly states: *"most complexes lie around the diagonal
    line and fall into the lavender margin of a 0.5 kcal/mol difference,
    only with a few exceptions due to the randomness of the Vina-GPU
    algorithm."* — i.e., a **0.5 kcal/mol band is the authors' own stated
    tolerance for "agreement."** Your 0.229 kcal/mol mean |difference| is
    well inside that band.
  - **9,125-molecule DrugBank virtual screen**: r = **0.981**, average score
    −7.9 (CPU) vs −7.8 (GPU).
  - Implementation detail confirmed directly in the text: the GPU kernel
    runs in **single-precision floating point (FP32)**. The authors do
    **not** attribute the score gap to precision — they attribute it to
    algorithmic randomness/search-space partitioning (see Q2 below).

### Vina-GPU 2.0, Ding et al., *J. Chem. Inf. Model.* 2023
- Ding J, Tang S, Mei Z, Wang L, Huang Q, Hu H, Ling M, Wu J. "Vina-GPU 2.0:
  Further Accelerating AutoDock Vina and Its Derivatives with Graphics
  Processing Units." *J Chem Inf Model* 63(7):1982–1998 (2023).
  DOI: `10.1021/acs.jcim.2c01504`, PMID: `36941232`.
- I could **not** get past ACS's paywall/bot-check to pull primary-text
  quotes for this one (403 on the abstract page, CAPTCHA on a proxy fetch).
  A secondary aggregator summary reports Pearson correlation coefficients of
  **0.971, 0.967, 0.910, and 0.976** for the various Vina-GPU 2.0 method
  variants (Vina-GPU+, QuickVina2-GPU, QuickVina-W-GPU, etc.) against their
  respective CPU counterparts, and a **65.6-fold** average speedup on RIPK1/
  RIPK3 kinase virtual screens "while ensuring comparable docking accuracy."
  **I flag this as unverified against primary text** — I'm reporting what a
  secondary source states, not something I read in the paper myself. Treat
  the specific correlation numbers as provisional until you or I can get the
  ACS PDF directly (institutional access would fix this).
- The 140-complex redocking set is reused across the Vina-GPU 2.0 and 2.1
  papers (referred to in the 2.1 preprint as the "AutoDock-GPU 140 dataset").

### Vina-GPU 2.1, bioRxiv preprint (2023) → published *IEEE/ACM TCBB* 2024
- Preprint: bioRxiv `10.1101/2023.11.04.565429` (posted 2023-11-06,
  **preprint, not yet peer-reviewed** at that DOI).
- Published, peer-reviewed version: "Vina-GPU 2.1: Towards Further
  Optimizing Docking Speed and Precision of AutoDock Vina and Its
  Derivatives," *IEEE/ACM Trans Comput Biol Bioinform* 21(6):2382–2393
  (2024). DOI: `10.1109/TCBB.2024.3467127`.
- **On the 140-complex set**, the paper states: *"There is no obvious
  difference between the methods in Vina-GPU 2.1 and the corresponding
  method in Vina-GPU 2.0 in terms of the RMSD... and the docking score."*
  This is a **GPU-2.1-vs-GPU-2.0** comparison, not a fresh CPU comparison —
  the paper leans on the earlier papers' CPU-comparison numbers rather than
  re-reporting them.
- **Critical for your Q4/Q1 overlap**: the headline enrichment numbers in
  this paper's abstract — *"average 4.97-fold acceleration... and an average
  342% improvement in EF1%"* — are **Vina-GPU 2.1 vs Vina-GPU 2.0**, not vs
  CPU AutoDock Vina. Section III.C/III.D report, on DrugBank and Selleck
  compound libraries respectively: **321%/312%/209%** improvement in
  Hit(1%/5%/10%) and **319%/313%/210%** improvement in EF(1%/5%/10%) on
  DrugBank vs Vina-GPU 2.0; **328%/440%/231%** and **365%/445%/231%** on
  Selleck. I found **no ROC-AUC, BEDROC, or CPU-baseline enrichment numbers**
  anywhere in the retrievable text of this paper. This is GPU-generation
  improvement, structurally not the comparison you need.
- Search-depth: *"The search depth in AutoDock Vina-GPU 2.1 is set
  heuristically... 1.5 times larger than those in Vina-GPU and
  Vina-GPU+."* No absolute numeric default is given in the text I could
  retrieve.

---

## 2. Is a small systematic bias (GPU less negative) expected, and why?

**Yes, and your finding replicates it almost exactly.** The original
Vina-GPU paper's 140-complex benchmark shows CPU **−8.9** vs GPU **−8.7**
kcal/mol on average — a **+0.2 kcal/mol** GPU-less-negative bias, essentially
identical to your **+0.188 kcal/mol**. That's a strong, independently
converging replication of the direction and rough magnitude of the effect
you observed, on a completely different receptor/ligand set (their 140
diverse PDB complexes vs. your 248 Pin1 ligands).

**On attribution:** the GPU implementations do run in single-precision
(FP32) — confirmed directly in the original Vina-GPU paper's methods text —
whereas CPU AutoDock Vina uses double precision. That is a real, documented
implementation difference. **However, the authors themselves do not
attribute the score gap to precision.** Their own explanation is algorithmic:
Vina-GPU partitions the search space into thousands of concurrent docking
"lanes"/threads (`--thread`), each running a shallower local search
(`--search_depth`) than a single CPU run's full exhaustiveness-driven
optimization, and they explicitly call the residual scatter "randomness of
the Vina-GPU algorithm" rather than a precision artifact. I did not find any
paper in this family that runs a controlled FP32-vs-FP64 ablation to
separate precision from search-algorithm effects — so the "is it precision
or is it search convergence" question is asserted by the authors, not
empirically decomposed by them. Take the algorithmic explanation as the
authors' claim, not as something independently verified in the papers.

Separately, AutoDock Vina's own documentation confirms the CPU tool itself
is non-deterministic and platform-sensitive: the official FAQ states exact
reproducibility requires identical random seed *and* identical platform, and
that "even minor changes to the input can have an effect similar to a new
random seed." The peer-reviewed "1001 Ways to Run AutoDock Vina for Virtual
Screening" paper (Jaime Feinstein & Brylinski, PMC4801993, *J Comput Aided
Mol Des*) goes further and reports that **the same CPU Vina binary, same
seed, same inputs, on different operating systems/platforms, produced energy
differences up to 0.7 kcal/mol** — larger than the CPU/GPU gap you or the
original Vina-GPU paper observed. Their practical recommendation is to run
Vina multiple times with different seeds rather than push exhaustiveness
higher, because seed/platform variance dominates.

**Net read on your Q2:** your 0.229 kcal/mol mean |difference| and +0.188
kcal/mol bias are (a) in the same direction and same rough size as the only
directly comparable published number I found (+0.2 kcal/mol, Vina-GPU
original paper), and (b) smaller than the documented platform-to-platform
variance of CPU Vina itself (up to 0.7 kcal/mol). Your reading that this is
ordinary stochastic/convergence variation is supported by what's published —
though "supported by what's published" is not the same as "mechanistically
explained by what's published"; nobody has published a precision-vs-search
ablation that would settle the causal question definitively.

---

## 3. What `search_depth`/`thread` settings do the authors recommend?

- **`--thread`**: the original Vina-GPU 2.0 paper's own tuning selected
  **8000** as the empirically optimal value "across complexity levels" — the
  exact value you used. The Vina-GPU 2.1 GitHub README documents a default
  of 8000 for most method variants (5000 for QuickVina2-GPU 2.1
  specifically) and warns thread count is "preferably less than 10000"
  given GPU memory/scheduling limits. **Your `--thread 8000` matches the
  authors' own recommended/default setting.**
- **`--search_depth`**: this is the parameter I'd flag as worth
  double-checking. Left unspecified, Vina-GPU computes it **heuristically
  per-ligand** from the published formula (Vina-GPU 2.0 paper):
  `search_depth = max(1, floor(0.24·N_atom + 0.29·N_rot − 3.41))`,
  where N_atom and N_rot are the ligand's (heavy) atom count and rotatable
  bond count. Vina-GPU 2.1 sets its heuristic default **1.5× larger** than
  this. **By explicitly passing `--search_depth 10`, you fixed a single
  value for all 248 ligands rather than letting the tool compute a
  per-ligand value from this formula.** For a mid-sized drug-like ligand
  (e.g., ~25 heavy atoms, ~6 rotatable bonds) the 2.0-paper heuristic
  formula evaluates to roughly 4; for larger/more flexible ligands it would
  give more. Whether 10 is "low" therefore depends on your ligand size
  distribution relative to what the heuristic would have picked — I could
  not find published guidance on an absolute floor below which
  `search_depth` is known to under-converge; the papers only publish the
  heuristic formula itself, not a sensitivity curve of score-vs-search_depth.
  **I did not find any paper reporting that deeper search demonstrably
  converges scores toward CPU Vina values** (no search_depth ablation showing
  monotonic convergence) — this is a gap, not a documented result. If score
  parity with CPU Vina specifically (not just speed) is your goal, letting
  Vina-GPU auto-compute `search_depth` per ligand (i.e., omitting the flag)
  rather than fixing it at 10 would more closely track what the authors
  tuned for.

---

## 4. Has anyone systematically compared enrichment (ROC-AUC/EF/BEDROC) between CPU and GPU Vina specifically?

**Largely no — this looks like a genuine, publishable gap, consistent with
your suspicion.** What exists:

- **Vina-GPU 2.1** reports EF/Hit-rate improvements, but strictly
  **GPU-2.1-vs-GPU-2.0** (see §1) on DrugBank/Selleck compound libraries used
  as a generic "virtual screening speed" demo, not an actives-vs-decoys
  enrichment benchmark against a CPU baseline.
- **Uni-Dock** (Yu et al., *J Chem Theory Comput* 19(11):3336–3345, 2023,
  DOI `10.1021/acs.jctc.2c01145`) is the one paper whose abstract/README
  claims direct comparability to CPU Vina on "screening power" using DUD-E
  (102 targets) and docking power using CASF-2016 (285 complexes), and a
  secondary source states Uni-Dock "generally performs better than AutoDock
  Vina on EF20%." **I was blocked by ACS paywall/CAPTCHA and a
  garbled/undecodable PDF from every route I tried**, so **I could not
  independently verify a single exact EF/AUC/BEDROC number from this paper's
  primary text.** I'm flagging the claim's existence, not certifying its
  numbers — please treat "Uni-Dock ≈ or > CPU Vina on DUD-E EF20%" as
  reported-but-unverified-by-me until someone with ACS access pulls the
  actual table.
- I found **no paper** — Vina-GPU family, AutoDock-GPU, QuickVina2-GPU, or
  gnina — that runs the specific experiment you ran: same receptor, same
  ligand set including a modest actives+decoys panel, CPU Vina and a GPU
  reimplementation scored independently, ROC-AUC/EF1% compared head-to-head
  with confidence intervals. The closest analog (GNINA benchmarking papers)
  compares a *different scoring function* (CNN rescoring) to Vina, which is
  a different question (scoring-function accuracy, not GPU-vs-CPU
  implementation parity).
- **Plain statement, as requested**: nobody appears to have systematically
  checked whether GPU-vs-CPU Vina score noise translates into enrichment/AUC
  differences on a controlled actives+decoys benchmark. Your CPU 0.535 →
  GPU 0.433 AUC result, with both CIs straddling 0.5 and EF1% = 0 for both,
  is — as far as I could determine — not contradicted or confirmed by any
  published number, because no directly comparable published number exists.
  Given both engines are at chance-level enrichment on your 5-actives/243-
  decoys panel, the "AUC drop" reads as noise around an uninformative
  baseline rather than a GPU-specific degradation — that interpretation is
  consistent with everything above (documented CPU-Vina-alone run-to-run/
  platform variance of up to 0.7 kcal/mol, well capable of moving a 5-active
  AUC around) but I want to be explicit that this is inference from adjacent
  literature, not a directly cited confirmation.

---

## 5. Documented cases where GPU docking changed conclusions?

**I found none**, in the sense of a published virtual screen where switching
CPU Vina → GPU reimplementation changed which compounds were prioritized/
pursued experimentally, or reversed a scientific conclusion. Searches for
"GPU docking changed ranking/conclusion/hit list" surfaced only tangential
material: studies about *rescoring* changing top-1% hit rates (e.g., an
AK-Score2 rescoring study reporting top-1% success improving from 22.8% to
71.9% after CNN rescoring of AutoDock-GPU poses — a rescoring-function
effect, not a CPU-vs-GPU-Vina implementation effect), and large-scale
screening papers (Uni-Dock's 38.2M-compound KRAS G12D screen, HASTEN's 1.56B-
compound screen) that use GPU docking as a *stage* in a funnel but don't
report a parallel CPU-Vina run to compare against. **This is an absence of
evidence, not evidence of absence** — it plausibly reflects that CPU-scale
docking of the ultra-large libraries these tools are built for is simply
infeasible, so nobody runs the CPU control needed to make the comparison.

---

## Sources

| # | Citation | Type | Verified against primary text? |
|---|---|---|---|
| 1 | Tang S, Chen R, et al. "Accelerating AutoDock Vina with GPUs." *Molecules* 27(9):3041 (2022). DOI `10.3390/molecules27093041`, PMCID `PMC9103882` | Peer-reviewed | Yes — fetched full text |
| 2 | Ding J, Tang S, Mei Z, et al. "Vina-GPU 2.0..." *J Chem Inf Model* 63(7):1982–1998 (2023). DOI `10.1021/acs.jcim.2c01504`, PMID `36941232` | Peer-reviewed | Bibliographic details yes (PubMed); numeric results **no** (ACS paywalled/blocked — numbers reported are from a secondary aggregator, flagged as unverified) |
| 3 | "Vina-GPU 2.1: towards further optimizing docking speed and precision of AutoDock Vina and its derivatives." bioRxiv `10.1101/2023.11.04.565429` (2023) | **Preprint** | Yes — fetched full text via proxy |
| 4 | Same title, published version: *IEEE/ACM Trans Comput Biol Bioinform* 21(6):2382–2393 (2024). DOI `10.1109/TCBB.2024.3467127` | Peer-reviewed | Not independently fetched (IEEE paywalled); treated as equivalent in content to the bioRxiv preprint I did read — flagged as an assumption |
| 5 | Yu Y, Cai C, Wang J, Bo Z, Zhu Z, Zheng H. "Uni-Dock: GPU-Accelerated Docking Enables Ultralarge Virtual Screening." *J Chem Theory Comput* 19(11):3336–3345 (2023). DOI `10.1021/acs.jctc.2c01145` | Peer-reviewed | **No** — paywalled/CAPTCHA/garbled PDF on every route tried; claims reported here are from secondary summaries only, flagged throughout §1/§4 |
| 6 | Jaime Feinstein W, Brylinski M. "1001 Ways to Run AutoDock Vina for Virtual Screening." *J Comput Aided Mol Des* (2015). PMCID `PMC4801993` | Peer-reviewed | Yes — fetched full text |
| 7 | AutoDock Vina official FAQ, `autodock-vina.readthedocs.io/en/latest/faq.html` | Official documentation, not peer-reviewed literature | Yes — fetched directly |
| 8 | iwatobipen blog, "Comparison between Vina-GPU and Vina – Is life worth living?" (2024-05-18) | **Informal blog, not peer-reviewed** — single-target, author self-describes as insufficient | Yes — fetched directly; cited only as weak corroborating anecdote (scores "well correlated but not same"), not as evidence |
| 9 | Santos-Martins D, et al. "Accelerating AutoDock4 with GPUs and Gradient-Based Local Search." *J Chem Theory Comput* (2021). PMCID `PMC8063785` | Peer-reviewed | Located but not fetched for numeric content — listed for completeness on AutoDock-GPU, not used for any claim above |
| 10 | DeltaGroupNJUPT GitHub READMEs (`Vina-GPU-2.1`, `Vina-GPU-2.0`) | Software documentation, not literature | Yes — fetched directly; used only for parameter defaults (`--thread`, `--search_depth`) |

**Not added to Zotero** — `$ZOTERO_USER_ID`/`$ZOTERO_API_KEY` were not
checked/available in this session; flag to me if you want these filed and
I'll do it in a follow-up pass with credentials confirmed.

**Not filed to the Oracle** — this is a targeted literature answer for
`inhibition`, written directly to the docs path requested. Let me know if
you want a curated summary entry pushed into your personal or lab Oracle as
well (I'd file it as `sensitivity: standard`, `tags: [docking, vina-gpu,
gpu-cpu-agreement]`, `project: inhibition`).
