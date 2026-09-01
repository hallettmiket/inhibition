# How this project breaks

*Written 2026-07-31, at handover. Last updated 2026-09-01 (catalogue 31
entries; the discovery-route table below refreshed for entries 22-31, which it
had never been). Read this before the README.*

Every substantive bug found in this project has been the same bug.

> **A value taken by position, name, or default rather than by identity —
> failing silently, because both the right and wrong candidates are populated
> and plausible.**

Not one of them raised an exception. Not one produced an obviously wrong
number. Every single one produced a number that looked exactly like the number
it should have produced, and was computed from the wrong thing.

If you learn one thing from this document, learn to recognise that shape. The
next instance is already in the code; it just hasn't been looked at yet.

---

## The catalogue

Each entry: what was taken, what should have been taken, and how it surfaced.
The last column is the one to study — **almost nothing here was caught by a
test, because the code was doing exactly what it was written to do.**

| # | The value taken | What it should have been | How it surfaced |
|---|---|---|---|
| 1 | The shortlist column named `shortlist` | `shortlist_synth`, the rebuilt list with rule failures replaced | The PI noticed unsynthesizable compounds still on screen |
| 2 | Active ids numbered `act_{i:03d}` by **frame position** | An id derived from the molecule | One active counted twice, another never entered — found while re-reading gate output |
| 3 | `rows[0]` of gnina's results table | The **affinity-best** pose; gnina sorts by CNN score | Found while building a multi-pose viewer, not while auditing the metric |
| 4 | `cnn_affinity` as T_3/T_4's rank metric in an analysis | `affinity_kcal`, the actual ranking column | The `LOWER_IS_BETTER` direction registry refused the column |
| 5 | A **hand-maintained** list of columns to drop before a merge | The columns actually being merged in | `_x`/`_y` suffixes appeared in the output frame |
| 6 | `warhead_classes_7.csv`, pinned by hand | `_8`, which existed and contained the correction | Noticed while reading the library, by luck |
| 7 | `pin1_reference_binders_3.csv`, pinned by hand | `_4`, written minutes earlier | Caught only because I went looking for who reads the file |
| 8 | A cache keyed on `class_id` alone | Keyed on `class_id` **and the query** | Served a stale 3-molecule pool after the query was relaxed |
| 9 | A cache keyed on `candidate_id` alone | Keyed on candidate **and protonation** | Caught *before* it bit, by asking what the key omits |
| 10 | `nunique()` on free-text `warhead_class` | Canonical `class_id` — one warhead, two prose strings | Counting gave exactly the floor of 6; the correct count is 4 |
| 11 | The default H-DAB-style assumption that alerts were gated | `alert_gate_pass` is `True` for all 4,803 T_1 rows | A blacksmith checked a molecule's alerts and found them computed but unused |
| 12 | `resi 101..125`, a **sequence window** | The measured 8 Å pocket shell | Excluded the entire Arg loop from a figure nobody had checked |
| 13 | The stale-module detector placed at the **bottom** of `app.py` | Before any helper attribute is read | It ran 30 lines after the crash it exists to explain |
| 14 | `rank_validated = verdict not in (UNDERPOWERED, UNGATED, FAIL)` | An allowlist: only `STRONG` validates | The gate moved to `WEAK` and T_1/T_2 silently began claiming a **validated** ranking. Noticed in a log line during an unrelated re-rank |
| 15 | `shortlist_synth` left in place by a re-rank | Dropped — it is derived from the ranking | Its members reached rank 90 under the new order, and the GUI prefers that column. The PI saw unsynthesizable compounds in one seed |
| 16 | `reshortlist_synthesizable.py` naming ONE T_2 experiment | All five seeds | du_xu and guo were never rebuilt, under a banner saying the filter was on |
| 17 | `shortlist_delta`'s `available` flag, computed and never read | Used to warn | Same defect as #11, one year of code later |
| 18 | `@lru_cache(maxsize=1)` on a **zero-argument** function | Keyed on the frames it is built from | A cache keyed on *less than nothing*, in a Streamlit process that outlives its data |
| 19 | `timeout=86400`, sized when a pool was 1,882 molecules | Scaled to the pool | 16,806 ligands ran 24 h, were killed, and wrote **0 poses**. A day of GPU time |
| 20 | mmCIF parsed as `loop_` only | Also the single-record `_tag value` form | Reported **zero** covalent entries across all 190 — which reads as a finding |
| 21 | Covalent ligands matched against the **adduct** pattern | The warhead-as-drawn: PDB component SMILES is the FREE ligand | Reported **11 novel chemistries**. Sulfopin itself came back "unclassified" |
| 22 | `bpmd_run.already_done()` keyed on `(ident, replicate)` | Keyed on `(ident, replicate, **trajectory length**)` | `bpmd/` still held `status == ok` replicates for two elevation-cohort molecules at 300 ps and 10,000 ps. A 3 ns between-group run would have SKIPPED both molecules and seated a 300 ps replica beside 3 ns ones. Caught by reading the existing chunk files before launching — nothing in the code would have said |
| 23 | `sweep_assets.rep_dir()` matched the **molecule** | `(molecule, pose_rank)` — the runner writes one directory per swept mode | Every mode of a multi-mode molecule was drawn from whichever sibling sorted first. 114 of 168 modes had the wrong trajectory in their plot and movie, median disagreement 0.478 nm. The rail's number was right beside a plot from another pose; @tt8804 read the two together and asked why they differed |
| 24 | `protonate()` protonated **one site**, then tested `charge == want` | Protons accumulated across sites, in basicity order, until the recorded charge is reached | It could only ever build a ±1 species, so every dication was stamped `docked_species_ok = False` and dropped — all 60 of D4's failures, **all BDHI**, zero acrylamide. Both BDHI arms would have entered the screen 15% short against a full acrylamide arm, and the only symptom would have been BDHI underperforming |
| 25 | The GUI read **flat** run directories (`attack_sweep/`, `md_residence/`) | Directories scoped to `run.topic`, as docking and ranking already were | Bumping the topic emptied one page of three; the other two kept showing 554 sweep rows and 647 residence rows from four superseded screens *under the new run's title*. The report server had the same defect, serving the old directory from a literal path |
| 26 | The directory named `<topic>_allposes` | The RAW cloud; it holds only poses whose DBSCAN label is in `mode_ids`, so 21% is absent | Five independent 500-run dockings returned 109-118 contact-space groups where subsamples of the raw 6,000-pose cloud returned 241-254. The box was suspected first and ruled out (both paths share one cached receptor). Centroid extent 19.19 A raw against 7.1 A filtered gave it away: **every candidate replacement for DBSCAN had been measured on clouds DBSCAN already cleaned**. `exp/5`'s docstring had said so for weeks, and travelled with nothing |
| 27 | rho = 0.657, the WITHIN-molecule atom ranking, as licence for the tolerance | The ACROSS-molecule absolute scale, never measured -- rho = +0.112, CI [-0.06, +0.27], crossing zero | @tt8804 asked whether the experiment behind the one calibration constant was big enough. It was 147 modes, but of the wrong quantity: the predictor varies at CV 0.15 between molecules where the truth varies at 0.45, and `median(rmsf)/2.21` does not beat writing ONE number down for every molecule (Wilcoxon p = 0.515). Nothing was broken -- exp/15 is careful work reporting an honest number, for a different question |
| 28 | A pose cloud analysed and DISPLAYED with no scores attached | The poses paired with the energies that ranked them | @tt8804 saw poses "literally outside of the pocket" in the viewer and asked how they could be lowest-energy. They are not -- 2.6% of the cloud, 88th energy percentile, zero in the best decile, and the scorer ranks them correctly (rho = +0.446 with exposure). But `nac_screen_v2` and `persist_raw_clouds` wrote coordinates with NO energies, so exp/16, 17, 19 and 20 all weighted the best pose and the 500th equally. Filtering to the best 25% concentrates 2.60x more than a RANDOM 25% of the same cloud (21 of 21 molecules, p = 6e-05) -- signal that four experiments had been averaging away |
| 29 | `0.0` from `anchor_quality` for a mechanism name nobody registered | A raise -- the mechanism list is an allowlist of four | Writing a test fixture with `mechanism="sn2"` instead of `sn2_displacement`. Every pose scored 0.0, which ranks the molecule LAST in a metric where 0 is the worst LEGAL value, and nothing raised. The symptom would have been a warhead class that never appears near the top of a shortlist. Fixed by raising, the rule `canonical_class()` already follows |
| 30 | `cut = -inf` as a way to disable a filter | Not evaluating the filter at all | The engagement gate correctly logged "the mode_poses gate does not apply and is not being enforced" and then returned 0 of 93 groups. `consensus_gnina` is null whenever gnina did not run -- every nac_v6 shard -- and `NaN >= -inf` is False. A permissive threshold on an absent column is not permissive, it is unsatisfiable, and the log said the opposite |
| 31 | "0 of 500 poses valid" from PoseBusters | "the receptor could not be parsed" | The gate was handed `3IKD_prepared_1.pdbqt`, the receptor the docking itself consumes. PoseBusters WARNS rather than raising on a file type it cannot read, so every protein-ligand check failed for want of a receptor and the verdict was total rejection -- an infrastructure fault wearing a chemistry result's clothes. Now the suffix is checked, and a zero-valid verdict names its worst failing checks and says to suspect the receptor |

---

## The four disguises

### 1. Selection by name

Two columns exist. Both are populated. Both are plausible. The code names one.

Instances **1, 4, 10, 21, 26**. The tell: a column name that *describes* what you
want rather than being derived from the thing that defines it. `shortlist`
sounds like the shortlist. `cnn_affinity` sounds like an affinity.
`warhead_class` sounds like a class. `adduct_attachment_smarts` sounds like the
right pattern for a molecule crystallised as an adduct — and is wrong, because
the PDB deposits the FREE ligand and the bond lives in `_struct_conn`.

**Defence:** derive the name from the thing that owns it. `rank_shortlist`
refuses any metric not declared in `LOWER_IS_BETTER` *with its direction* —
that guard caught #4 and is the only reason it was caught.

### 2. Selection by position

The first row. The *i*th element. The order the file happens to be in.

Instances **2, 3, 5**. The tell: an index used where an identity was meant.
`rows[0]` is correct only if the row order encodes the property you want — and
in #3 it encoded a *different* score than the one being read off it.

**Defence:** if you index, write down what guarantees the order. If you can't,
you have this bug.

### 3. Pinned defaults that go stale

A constant naming a version. It was right when written. The file it names still
exists and still parses, so nothing ever complains.

Instances **6, 7, 15, 18, 19**, and the pattern behind **8, 9, 22** (a cache key
is a pin on its inputs). Reference files alone have done this **five** times:
generalising the version guard on 2026-08-01 found five stale pins at once, two
of them in `choreography.yaml` naming a binder set and a warhead library two
and eight versions behind what the code loaded — and read by nothing, so
nothing ever complained.

The family is wider than version literals. **#19** is a `timeout=86400` sized
when a pool was 1,882 molecules; at 16,806 it killed a 24-hour run that had
written nothing. **#18** is a cache keyed on *nothing at all*. **#15** is a
derived column that outlived the ranking it was derived from. All the same
shape: a value that was right when written, cannot announce that it is not any
more, and costs you before it raises anything.

**Defence:** `shared/reference_set.latest_reference()` resolves by glob;
`tests/test_reference_version.py` walks the AST across `shared/`, `scripts/`,
`approaches/`, `integration/`, `tests/` and `config/*.yaml`, over every stem in
`data/reference/`. A pin cannot announce that it is out of date — only a
comparison against the directory can. For non-version constants: scale them to
their input and **log the derived value before the work starts**, so a run that
cannot finish says so at launch rather than at the deadline.

### 4. A guard that is scoped out, mis-ordered, or vacuous

The check exists. It runs. It just doesn't cover the case, or it runs too late,
or it cannot fail.

Instances **11, 12, 13, 14, 16, 17, 20, 27, 28, 29, 30, 31** — the largest group, and growing
fastest. **#27** is the subtlest: the validation RAN, honestly, and reported a
real number -- for a different question than the one its result was used to
settle. A validation is scoped to the quantity it measured, and "the predictor
is validated" names a tool rather than a claim. **#14** is the sharpest: `rank_validated` was computed as `verdict not
in (UNDERPOWERED, UNGATED, FAIL)`, a DENYLIST, so when the gate produced a
verdict nobody had anticipated (`WEAK`) the ranking validated itself by
default. **#17** is #11 again in new clothing — a flag computed and never read.
**#20** is a parser that handled one of mmCIF's two serialisations and reported
zero covalent entries across all 190.

Plus two of my own tests during this session:

* `test_stale_guard`'s first version counted only dotted `curate.X` reads.
  Making the crashing line defensive with `getattr` left it **zero accesses to
  check** — it passed for free and would have kept passing through a
  reintroduced regression.
* The version-pin test flagged four "offenders" that were all **docstrings**.
  Fixing the prose was right; so was narrowing the test to executable code.

**#28** is the one where nothing was wrong with any number: the poses were real,
the geometry was fine, PoseBusters passed them, and the statistics were correct
over exactly the set they were given. The analysis answered "what does the whole
cloud look like" while every reader took it as "what does the docking think" --
and those differ by the 75% of poses the score rejects. A measurement is scoped
to the population it ran on, and a pose set displayed without its scores states
the wrong population silently.

**#29, #30 and #31 are all one week's worth of the same shape**: a guard that
reaches a verdict by a route nobody checked. #29 scores an unregistered mechanism
0.0 rather than raising, so a typo ranks a molecule last in a metric where 0 is
the worst LEGAL value. #30 disables a filter by setting its threshold to `-inf`,
which is unsatisfiable rather than permissive the moment the column is NaN -- and
the log announced the filter was off while it rejected everything. #31 reports
"0 of 500 poses valid" when the truth is "the receptor did not load", because the
library warns instead of raising. In all three the output is a legal value in the
right range, produced for a reason nobody intended.

**Defence:** for every guard, ask *what would make this pass when it should
fail?* If the answer is "the thing it inspects being absent," assert the thing
is present. **Name what passes, never what fails** — an allowlist refuses a
value nobody anticipated; a denylist admits it.

---

## Why this project produces this bug so reliably

Not carelessness. Three structural causes, all of which will still be there
after handover:

**The pipeline is wide and every stage has near-miss neighbours.** Four
approaches × several metrics × several selection columns × versioned reference
files. Almost every value has a plausible sibling one identifier away.

**Nearly every quantity is a plausible float.** A binding affinity of −5.85 and
one of −6.86 both look like binding affinities. There is no type error, no
range violation, nothing a schema catches. Only a *comparison against what the
value was supposed to mean* catches it.

**The outputs are read by eye at the end of a long chain.** By the time a
number reaches the GUI it has passed through generation, annotation, docking,
ranking, shortlisting and merging. A defect at stage two is indistinguishable
from a real result at stage six.

---

## How they actually got caught

Worth being honest, because it should change how you spend your attention.

| route | 1-21 | 22-31 | total | which of 22-31 |
|---|---:|---:|---:|---|
| Someone looked at output and it didn't match expectation | 9 | 7 | **16** | 23, 24, 25, 26, 28, 30, 31 |
| Found while building something else entirely | 6 | 1 | **7** | 29 |
| An existing guard fired | 3 | 0 | **3** | none |
| Deliberate audit for this class of defect | 3 | 2 | **5** | 22, 27 |

*Updated 2026-09-01 with entries 22-31. **The ratio got worse: still only
3 of 31 caught by a guard**, because not one of the ten new entries was. Seven
of the ten surfaced because someone read output that did not match
expectation, and two of those (#23, #28) came from @tt8804 seeing two things
presented side by side and asking why they disagreed — a route no guard
performs. The 22-31 column names which entry went where so the next person can
check the arithmetic instead of trusting it; the 1-21 column is left as
originally judged and is not re-derived here. This table went four weeks and
ten entries without being updated, which is itself a pinned default gone
stale (§3).*

*Updated 2026-08-02 with entries 14-21. Two of those (#20, #21) were caught by
checking against **6VAJ, whose covalent linkage D0001 already records at
1.78 Å** — ground truth the project already had, in a decision record, costing
one query to check. Both had produced output that read as a discovery: "zero
covalent entries in the PDB" and "11 novel chemistries". Neither would have
been caught by the number looking wrong.*

**If you take one operational habit from this document, take that one:** before
trusting a new pipeline, run it on something whose answer is already written
down in `decisions/`.

**The two that a guard caught are the two where a guard existed.** The direction
registry caught #4; the append-only manifest caught a provenance error. Every
other one needed a human noticing something looked off.

That ratio is the argument for writing guards — not for auditing harder.

---

## What to do when you find the next one

1. **Fix the class, not the case.** #7 wasn't fixed by editing a constant to
   `_4`; it was fixed by making version resolution impossible to pin.
2. **Write the guard so it fails loudly.** `canonical_class()` *raises* on an
   unmapped warhead rather than returning it unchanged, because returning it
   unchanged is how it becomes its own chemotype.
3. **Check the guard can fail.** See the two vacuous tests above.
4. **Write a decision record** — including what was wrong and *why it looked
   right*. The "why it looked right" is the part that helps the next person.
5. **Ask what else shares the shape.** #9 was caught *before* it bit, by asking
   what a cache key omits right after #8 taught us to ask.

---

## Where I would look next

*Updated 2026-08-02. **Four of the five leads below paid off within two days**,
which is the strongest argument this document makes: the shapes are
predictable, and looking for them deliberately works better than waiting.*

Resolved since first written:

* ~~Every other cache~~ → **#18**, `@lru_cache(maxsize=1)` on a zero-argument
  function. Now keyed on the frames it is built from.
* ~~The remaining hand-maintained lists~~ → **#5** fixed by deriving the drop
  list from the merge; the same pattern applied to the docking merge before it
  could bite there.
* ~~Anything else computed-and-unused~~ → **#17**, `shortlist_delta`'s
  `available` flag. Found by asking the question this lead poses.
* ~~Every `[0]` on a scored frame~~ → D0047: `affinity_kcal` was the CNN-best
  pose's affinity. 89% of covalent candidates affected.

Still untested:

* **`n_docked` and other denominators.** If any is computed after a filter a
  reader doesn't know about, every percentage derived from it is wrong. This is
  the one lead from the original list nobody has checked.
* **Every remaining denylist.** #14 was `verdict not in (...)`. Grep for
  `not in (` and `!=` in any predicate that decides whether something is
  valid, trusted, or safe, and ask what an unanticipated value does.
* **Every constant with a unit.** #19 was a timeout. Any literal seconds,
  bytes, counts or thresholds sized against a workload that has since grown —
  and which will fail by killing work rather than by raising.
* **Every parser that assumes one input shape.** #20 handled one of mmCIF's
  two serialisations. Anything reading an external format: what is the *other*
  valid way a producer could write this?
* **Every column pair where one name is the reacted/derived form of the
  other.** #21 was adduct-vs-warhead. `*_rows0` vs the fixed column, `rank` vs
  `rank_synth`, `shortlist` vs `shortlist_synth` — all live pairs today.

---

## The one-sentence version

If a number in this project is wrong, it is almost certainly because something
selected it by **where it sat** or **what it was called**, rather than by
**what it is** — and it will look completely normal.
