# Manuscript changes implied by the current repo

*@twu383 with Claude Code, 2026-08-25. Companion to the four Pin1 figures.*

Against `murmurent_manuscript` @ `main`, version string
**v2026.07.26-1254**. Every "actual" below is parsed from the repo at
2026-08-25, not read off a document, and the parse is reproducible with
`python scripts/manuscript_figures.py`, which prints what it parsed before it
draws anything.

**Why this file exists.** The manuscript quotes the Pin1 choreography's own
counts, and those counts have moved. They moved in the direction that flatters
nothing — the corpus is larger, the failure catalogue is longer, and the
receptor changed — but a reviewer who clones the repo the paper points at will
find different numbers than the paper states, and the paper's central claim is
that this project's records are trustworthy.

---

## A. Numbers that are now wrong

Each is stated in the manuscript as a bare figure. All are one-token edits.

| § | manuscript says | actual, 2026-08-25 | note |
|---|---|---|---|
| 1.5 | "generated **51** such records" | **87** | |
| 1.5 | "**37** of the **51** records bind the whole choreography" | **72** of **87** | ratio strengthens, ~2:1 → ~5:1 |
| 1.5 | "**34** … implementation, **10** from the PI, **6** from an Adversary audit" | implementation **54**, user **24**, adversary **8**, spec **1** | see note (i) |
| 1.5, 1.6, A | "the Pin1 exploration has **21** entries" / "**two of the twenty-one** defects" | **25** | four places |
| A, Table A1 | accepted **47** · partially withdrawn **2** · superseded **1** · proposed **1** | **61** · **2** · **1** · **23** | see note (ii) |
| A, Table A1 | "**356** separate pieces of evidence, a median of **seven**" | **645**, median **7** | median holds |
| A, Table A2 | 4 · 3 · 7 · 7 (sums to 21) | 4 · 3 · **8** · 7 = 22 assigned, **3 unfiled**, of **25** | see note (iii) |
| A, Table A3 | 9 · 6 · 3 · 3 (sums to 21) | unchanged, but now covers **21 of 25** | see note (iii) |

**(i) "from the PI" is not the same field as `origin: user`.** The record schema's
`origin` allowlist is `spec | adversary | implementation | user`; `decided_by` is
separate and currently reads `@mhallet` 50, `@tt8804` 35, `@twu383` 1, plus one
`@mhallett` — **a typo of the PI's own handle that splits his count**. Decide
which field the sentence means before changing the number, and fix the typo
either way; it is a one-character instance of the project's own selection-by-name
failure mode.

**(ii) `proposed` went from 1 to 23, and that is a finding, not bookkeeping.**
Several load-bearing decisions are currently `proposed`, not `accepted` —
including **D0059**, which replaces the receptor, and **D0088**, which says the
modes are mixtures. The manuscript's Table A1 currently implies a corpus that is
almost entirely settled. It is not, and the honest version of that table is more
interesting: it shows a project with live, unresolved calls in it.

**(iii) The catalogue's own sub-tables have not caught up with its entry list.**
`how_this_project_breaks.md` lists 25 entries, but the four "Instances **…**"
lines assign only 22 of them and the how-it-was-found table still sums to 21.
So Tables A2 and A3 cannot simply be re-typed — **the source document needs its
own counts refreshed first.** Figure 4 renders this honestly (a "not yet filed"
bar and an explicit *n* on the routes panel) rather than silently presenting a
partial total, but the underlying doc should be fixed.

---

## B. Statements that are no longer true

These are not number edits. Each needs a sentence rewritten.

### B1. The receptor — §1.2, §1.8, §A

> "docking makes use of one prepared Pin1 receptor (PDB 6VAJ)" (§1.2)
> "one hash-pinned receptor fixed at the outset and binding all four" (§1.8)

**D0059 (2026-08-05) replaces 6VAJ with a chemist-prepared 3IKD**, because 6VAJ
is co-crystallised with sulfopin and its pocket is induced-fit around that
ligand. Cross-docking into 6VAJ ranked the crystal pose #1 in **0 of 82** cases.
Re-measured on 3IKD, pose recovery went **6.1% → 18.3%** top-1 and
**15.9% → 41.5%** best-of-9.

Two consequences the manuscript should carry:

1. **Every 6VAJ measurement is invalidated until re-run.**
2. **The claim "one hash-pinned receptor binding all four" is, right now, false
   in a way that is worth reporting rather than hiding.** `config/receptor.yaml`
   still pins 6VAJ and `shared/noncovalent_dock_run.py:61` still hardcodes
   `6VAJ_prepared.pdbqt`, while the benchmark and reference-screen paths call
   `resolve_3ikd_ian()`. **Which receptor you get depends on which entry point
   you came through** — and D0059 is still `proposed`.

This is the paper's own thesis about shared substrate playing out live: one
shared artefact, changed by decision, not yet propagated to every caller. §A
already argues that shared infrastructure has a shared blast radius. This is
that argument with a current example, and I would use it rather than paper over
it.

### B2. "Roughly 72,000 candidates" — §1.3, §1.4, §A

The repo's **generated** orientation table (`scripts/refresh_orientation.py`,
guarded by `tests/test_orientation_current.py`) totals **56,579 rows** across
current frames: T₁ 4,803 · T₃ 5,396 · T₄ 1,783 · T₂ all six seeds 44,597.
Docked is lower still, ~**53,500**.

The ~72,000 traces to hand-maintained prose in `state_of_the_project.md` §1
("13,863 … plus 42,588 … plus 15,653"), which the generated table now
contradicts — the degree-2 ATRA sample it counts at 15,653 is **127 rows** in
the current frame. **I have not resolved which is right**, and it should not be
resolved by picking the nicer number. Re-run the generator and quote it, and fix
§1 of the state document to match.

### B3. The enforcement claim — §1.4 Tier 3, §1.6

> "two filesystem invariants which are enforced by hooks" (§1.4)
> "Neither guarantee depends on the agent's cooperation" (§1.6)

The choreography's own repo now says the opposite, and says it was verified:

> `immutable/` and `append_only/` are a **DISCIPLINE, not enforcement.** …
> Verified 2026-08-02: both trees are writable at the filesystem level …
> The guarantee comes from `~/.claude/hooks/block-rm.sh`, a **per-user Claude
> Code hook**: it does nothing in a plain shell, nothing for a script run
> outside a CC session, and nothing for a different user until they install
> their own. — `state_of_the_project.md` §8

§A already concedes most of this ("Both hooks are guardrails rather than a
security perimeter — a Bash-invoked script that shells out to a system call the
hook does not parse can bypass them"). **The main text should be brought into
line with the appendix and with the measurement**, because §1.6's phrasing is
the strongest claim in the paper and the easiest for a reviewer to test. The
Results section says the guarantee is unconditional; the supplement says it is
not; the repo says it was checked and it is not.

### B4. "Nine warhead classes" — §A, simultaneity subsection

> "T₄ enumerates nine warhead classes against a fixed 198-member R-group library"

Nine classes, but **seven distinct docked species**. Chloroacetamide, sulfamate
acetamide and sulfonate acetamide all collapse to `acetamide_adduct` and share
an identical `adduct_smiles`, so they dock as literally the same molecule and
score identically to five decimal places:

```
t4_617de1a16274  chloroacetamide       -5.43626
t4_97020e8242c4  sulfamate_acetamide   -5.43626
t4_0a52a08197fc  sulfonate_acetamide   -5.43626
```

The enumeration is nine; the search is seven. Both are worth saying, and the
collapse is a clean illustration of the adduct-form transform the same section
already discusses.

---

## C. Results the manuscript does not yet contain

### C1. The 3.1.0 campaign — absent entirely

561 molecules screened, 4,432 modes ranked, 147 triaged at 8 ns, **15 given
100 ns runs; five held the pocket, one held unstably, nine left.** This is the
closest thing the choreography has to an endpoint and the manuscript stops
before it.

**It must be reported with its caveat attached** — see C2 — but that pairing is
the paper's argument, not a weakness in it.

### C2. D0088 — the modes are mixtures

The pipeline **clustered on the reactive atom's position and warhead direction,
then scored each group by how often it reached attack geometry — which is
position and direction.** It formed groups along the axis it then graded them
on. The code comment asserting otherwise ("never on the NAC geometry itself,
which is the score") is wrong.

Measured over `nac_v5`: the median mode spans **3.51 Å**, 87% span more than
2 Å, **42% have a viable fraction between 0.1 and 0.9**, and the largest holds
137 poses across 9.3 Å. So `viable_fraction`, `enrichment` and `conditional_eb`
are all measured over mixtures. The 100 ns results survive as measurements —
a trajectory is a trajectory — but **which pose earned each one was chosen by
this machinery.**

This belongs in §1.9 beside the MM/GBSA audit, as a **third** instance of the
same pattern the paper already tells twice (D0015→D0028→D0031, and the
adduct-form correction): a number the project was ready to report, withdrawn by
its own audit. Reported that way it strengthens the paper's thesis. Omitted, it
is the largest thing the repo knows that the paper does not.

### C3. Pose recovery as a positive control — §1.9

The manuscript argues docking is weak on Pin1 from enrichment alone. The repo
has a second, stricter, cheaper measurement (D0046, corrected in #66) that says
**the failure is in the scoring, not the sampling** — and that this is the
documented behaviour of the exact program used. Across programs, sampling finds
a near-native pose in 85–99% of cases while scoring ranks it first in 35–73%;
**AutoDock Vina shows the largest gap of any program tested, 93.4% sampling
against 35–40% ranking.**

Two corrections travel with it, and #66 exists because the first framing was
wrong:

- The old comparison was against a **60–80% norm, which is the SELF-docking
  norm.** The 82-case benchmark is **cross-docking**, where the published
  baseline is ~**41–50% top-1** for a single receptor.
- So top-1 **18.3%** is below the cross-docking norm by about **2×, not 3×** —
  and **best-of-9 41.5% is *at* the single-structure cross-docking norm**, which
  the old framing presented as a failure.

Figure 5b shows this, including that Vina's top-1 (18.3%) is indistinguishable
from **randomly picking one of the nine modes (19.8%)**.

---

## D. The figures

Four, all regenerable, all writing versioned output under
`append_only/inhibition/00_outputs/artist/manuscript_figures/`.

| | file stem | script | what it carries |
|---|---|---|---|
| **2** | `f2_gui_output` | `scripts/manuscript_gui_figure.py` | **the missing "output of the Pin1 search" figure** — real captures of the integration GUI: four shortlists side by side, and one candidate's dossier |
| **3** | `f3_control_walked_to_chance` | `scripts/manuscript_figures.py` | D0015 → D0028 → D0031 with CIs against chance, plus the non-covalent gate |
| **4** | `f4_written_substrate` | ″ | decision corpus (scope, origin) and failure catalogue (disguise, discovery route) |
| **5** | `f5_levels_of_theory` | ″ | both decoy gates against chance; pose recovery against the cross-docking band |

**Figure 2 is captures, not a mock-up.** Every pixel inside a panel came from
`integration/app/app.py` served by `scripts/serve_gui.sh` and photographed with
headless chromium. The only edits are the crop, the panel letters and the margin
callouts, all outside the captured area.

**The figures parse their own numbers.** Nothing in Figure 3, 4 or 5 is typed
in; the counts come from `decisions/*.md`, `how_this_project_breaks.md` and the
enrichment-gate token at render time, and the script prints what it parsed. A
figure that hardcodes a count cannot announce that it is stale, which is
disguise #3 in the project's own catalogue — and is how the manuscript got here.

### Suggested placement

- **Figure 2** → §1.2, at "presented side by side along with supporting evidence
  … within a unified graphical user interface". This is the figure §1.2 has been
  describing in prose.
- **Figure 3** → §1.9, at the ROC-AUC 0.815 → 0.718 → 0.537 sequence.
- **Figure 4** → §1.5, replacing or accompanying Tables A1–A3. **It renders
  Tables A1 and A2 obsolete**; A3 could stay as a table if the source doc's
  routes are refreshed to cover all 25.
- **Figure 5** → §1.9, where docking's weakness on this pocket is asserted.

### Palette

T₁–T₄ are `#332288 · #DDAA33 · #44AA99 · #CC3311`, in fixed order, never cycled.
Validated for colour-vision deficiency: worst-case OKLab ΔE across
normal/deuteranopia/protanopia/tritanopia is **13.1** against a floor of 8 for
CVD and 15 for normal vision. If the figures are extended, keep the order and
re-run the check rather than adding a hue by eye.

---

## E. Smaller things

- **§1.9 "one prepared Pin1 receptor"** — also appears in the §1.2 approach
  table caption; same fix as B1.
- **§A "its 1,782 candidates"** — the T₄ frame is **1,783** rows, **1,683**
  docked. The 1,683 figure elsewhere in the same section is right.
- **Table 1 "fourteen reference agents"** — still correct (14 in `agents/`).
- **§1.5 "51 such records, with each carrying structured frontmatter"** — the
  schema block reproduced in §A is still accurate, including the `origin`
  allowlist. D0051 made `origin` refuse invented values; worth one clause, since
  the paper argues elsewhere for allowlists over denylists and this is the
  schema practising it.
- **Author list** — I am on it as Timothy Wu (`twu383@uwo.ca`); the repo work
  since 2026-08-02 is under `@tt8804` and `@twu383`. Worth checking the handle
  mapping in the acknowledgements if agent-attributed work is credited.

---

## F. What I did not touch

**No manuscript file was edited.** `~/repos/murmurent_manuscript` was cloned and
pulled per the pull-first rule, and left clean — Overleaf edits are
authoritative on conflict, and every change above is a judgement about the
argument, not a mechanical substitution. Say which of A–E you want and I will
make them as one coherent block and push, so Overleaf can fetch it.
