# Since the handoff

*What has changed since @mhallet left on 2026-08-02. Newest first. Kept by hand;
add to it when you change something a reader would be misled without.*

This is not a changelog — `git log` is the changelog. This records the things
that **change what someone should believe**: results that moved, defects that
were silently producing plausible numbers, and decisions that redirected the
work.

---

## The headline: the receptor was costing us 2.6× of the findable poses

**2026-08-05.** @tt8804's chemist recommended abandoning 6VAJ for a prepared
3IKD. Re-running D0046's pose-recovery benchmark on it:

| is *any* pose in the top-k within 2 Å | 6VAJ | **3IKD** |
|---|---:|---:|
| top-1 — *what the pipeline carries* | 6.1% | **18.3%** |
| top-3 | 9.8% | **29.3%** |
| top-5 | 12.2% | **34.1%** |
| **best-of-9 — *the ceiling*** | **15.9%** | **41.5%** |
| **random pick of the 9 — *the floor*** | 5.3% | **19.8%** |

**The right pose is now in the ensemble 2.6× more often.** 6VAJ is co-crystallised
with sulfopin, so its pocket is induced-fit around that ligand; the improvement is
the size of what that was costing.

**And the score is indistinguishable from chance at picking it.** Random selection
among the nine gives **19.8%**; Vina's score gives **18.3%**. The rank of the
closest-to-crystal mode is essentially uniform (5, 12, 9, 12, 11, 12, 9, 6, 6
across ranks 1–9).

So the headroom for pose selection is **18.3% → 41.5%** — more than doubling
recovery. But the bar any method must clear is **random (19.8%)**, not the score,
and those are effectively the same number.

**Correction, 2026-08-05:** an earlier version of this table reported top-1 as
2.4% and a 17× gap. That was a reporting bug of mine, not a docking result:
`top-k` was computed as *"the single closest pose is within the top k"* rather
than the standard *"a good pose exists within the top k"*, which excluded cases
where mode 1 was within 2 Å but some later mode was closer. Both arms were
understated. The direction and conclusion are unchanged; the magnitudes were
wrong. Catalogue shape, committed here.

**Caveats that travel with these numbers.** The comparison is confounded — 6VAJ
was water-stripped and `reduce`-protonated, 3IKD keeps 6 waters and the chemist's
protonation (D0059), so some of the gain may be solvation rather than conformer;
the control is one cheap docking run and has not been done. And the RMSD here is
an optimal-assignment proxy, not D0046's graph-matched symmetry-corrected metric,
so **these absolute rates must not be quoted against D0046's published 5%/55%** —
both arms used the same proxy, so only the comparison between them is sound.

---

## Direction changes

**The receptor is now 3IKD** (D0059), used exactly as the chemist prepared it —
no re-protonation, waters kept, nothing rebuilt. Only J9Z was removed, because a
receptor holding its own cognate ligand has an occupied pocket; its coordinates
define the box. Work is on branch `receptor/3ikd-chemist-prepared`.

**Every 6VAJ measurement is invalidated** and must be re-run before it is quoted
about 3IKD: D0016/D0041 enrichment, D0046 pose recovery, D0031, D0049.
`attach_gate` keys on *(stratum, metric)* and the metric NAME does not change
when the receptor does, so an old verdict would attach silently. Generation is
unaffected — all **86,451 molecules are receptor-independent**.

**Ranking is being rebuilt around geometry, not affinity** (#14). The framing is
@tt8804's: *can the molecule orient to form the bond, not how good the bond is.*
Non-covalent docking of the **free form** (never the adduct), then pose selection,
then a near-attack-geometry criterion that is **mechanism-specific** — SN2 needs
backside attack anti to the leaving group, Michael addition needs approach
perpendicular to the alkene plane. A single distance cutoff would be wrong for
every class.

**Docking score is out of the ranking.** Five independent levels have now been
measured and failed: docking enrichment (D0041), pose recovery (D0046), MM-GBSA
(D0036), MD residence (D0038/D0044), and the contact-profile fit score built and
killed here (D0057).

---

## Results that were measured and did not survive

**The contact-profile fit score is worse than chance** (D0057). Built to answer
#14's request for an agnostic fit score; scored on 82 crystal cases with a
leave-one-out reference profile before being allowed to label anything, per #13's
pre-registration rule.

| | crystal ranked #1 | top 3 |
|---|---|---|
| cross-docked into 6VAJ | 0/82 (0.0%) | 3.7% |
| self-docked into own receptor | 5/82 (6.1%) | 22.0% |
| chance | 10% | 30% |

Recorded as failed rather than tuned until it passed on the same 82 cases — that
would be the D0045 failure in a new costume. The self-vs-cross gap became
independent evidence for the receptor change months before the chemist's
recommendation was tested.

Two things survived it: `representative()` returns an **index into real poses**,
so a "weighted average pose" can never be a synthesised conformation; and the
contact vector itself, which is reusable.

---

## Defects found, and the one pattern behind them

Every one produced a populated, plausible, wrong result. None raised an
exception until someone went looking.

| what was wrong | how it presented | record |
|---|---|---|
| **Vina-GPU segfaulted for anyone but its owner** — it recompiles OpenCL kernels at startup and writes them into a directory only `@mhallet` can write | driver **exit 0** while all six chunks failed: **127 of 30,000** molecules docked, and the pose directory still held the previous run's 15,653 files | **D0060** |
| `merge_poses_onto_frame` raised `UnboundLocalError` whenever `df` was passed — the whole of T_1/T_2's ordinary dock route | failed at the manifest write, **after the entire GPU run was spent** (~1.3 h atra, ~10 h du_xu) | **D0053** |
| Degree-2 sampler estimated the population by **linear** extrapolation of a **sublinear** union | kept **15,653 against a target of 30,000** — unbiased, but half the intended size, and the docstring said the size "varies slightly" | — |
| **112 T_1 survivors silently failed ligand prep** | `rejected_at = NA`, `vina_affinity = NaN` — indistinguishable from "queued, not yet run" | — |
| `find_pose` preferred the untagged pose directory over `poses_ph7.4/` | the GUI drew the **superseded neutral-form pose** beside a pH-7.4 score | — |
| Guards testing `is_file()` where they needed readability | **13 tests failed with assertions about Arg-loop residues** when the cause was a filesystem ACL | **D0054** |
| T_1's `num_nodes_lig` declared in config, read by no code | setting it did nothing; it controls ligand **size**, directly upstream of D0043's size finding | — |
| T_2's radius probe omitted `max_inc` and hardcoded sizes | the safety check measured **1,792** where the run produces **1,882** | — |
| T_2's "yields NOTHING" warning checked `mutate` only | `grow` supplies **89%** of ATRA's frontier; a seed whose grow collapsed would pass silently | — |
| **mmCIF-only structures were never cached** | D0046 ran on **82 of 85** cases, and the 3 missing were systematically the *newest* entries — 5-character CCD codes the legacy PDB format cannot hold | — |
| A fix landing in a generation stage is **inert until the stage re-runs** | T_3's pocket ceiling postdates its frame by 3 days; **no T_3 frame has ever carried the column** | **D0056** |

**The pattern**, which is `how_this_project_breaks.md`'s thesis holding up: a
value taken by position, name, or default rather than by identity — failing
silently because the right and wrong candidates are both populated.

**The new one worth adding to the catalogue:** *a driver's exit code is not
evidence that the work happened.* `dock_chunked` printed `127 / 30000` and it was
not read.

---

## Guards added

Each fires on the *class*, not the instance.

- **`test_frame_code_currency`** (D0056) — fails when a generation stage's code
  changed after its frame was produced, unless the impact is **measured** and
  recorded. A second test asserts every acknowledgement contains the word
  MEASURED and a number.
- **`test_orientation_current`** (D0055, closes #11) — the orientation doc's
  numbers are generated, and drift **fails the suite**. Caught 800 molecules of
  drift on its first run, then fired again on its author's own change.
- **`test_dock_merge_provenance`** (D0053) — runtime test plus an **AST guard**
  that fails for any future caller passing `df` without `frame_path`.
- **`test_degree2_reservoir`** — exact sample size and uniformity over the whole
  stream, including a test that the *tail* is not starved.
- **`test_pose_vector`** — a deliberately bimodal fixture asserting the medoid
  lands on a real mode while the mean lands in the trough where nothing is.
- **`test_untested_modules{,_2}`** — coverage of decision-affected modules went
  **23/38 → 30/38**. Most load-bearing: `replicate_seed`, which **D0044's
  conclusion depends on** — five replicates are only independent if their seeds
  differ, and nothing checked it. They do.

---

## Infrastructure

- **Data-root access.** 0 of 166 files were readable by `@tt8804`; an Isilon ACL,
  invisible to the NFSv3 client, denying despite POSIX mode `2770` and correct
  group membership. Fixed by @jmucaki. `/data/lab_vm/envs/` was **not** covered —
  hence D0060.
- **PLUMED 2.10 installed** and verified end to end with GROMACS 2026.3
  (`Plumed support: enabled`). BPMD is unblocked. Proper home is
  `/data/lab_vm/envs/dwi_plumed`, which needs write access there.
- **Enumeration is ~10× faster.** The parent process was the bottleneck at 99.7%
  CPU with 48 workers idle at 59% — three RDKit calls per child on one thread.
  Moved into the workers. Then the fragment DB turned out to be a 2 GB SQLite
  file **on NFS** with 38% of workers in uninterruptible I/O wait; copied to
  `/dev/shm`, SHA-256 verified, manifest still recording the governed path (#15).
- **Measured BPMD cost:** 485.8 ns/day under bias on an A100 (30% PLUMED
  overhead), so ~20,400 ns/week on 6 GPUs — **200–1,000 poses**, depending on
  protocol.

---

## Issues

Consolidated 2026-08-04 from a tangle of seven into two, then two more opened:

- **#12** — chemistry judgement for the Lu lab. **Section B of the old #8 was
  wrong**: it promised four extra warhead chemistries, and verification against
  `_struct_conn` killed three (SuFEx bonded to TYR23; aryl aldehyde and maleate
  ester have no covalently-linked instance). The real position is **3 verified
  Cys113 chemotypes against a floor of 6** — fewer than we had.
- **#13** — every open technical problem, audited against the code.
- **#14** — the docking/MD framework redesign.
- **#15** — enumeration efficiency, with measured bottlenecks.

`#2`, `#4`, `#6`, `#8`, `#9`, `#10`, `#11` closed into those.

---

## The ranking framework now separates actives from inactives, for one class

Built and measured on 2026-08-05. Full reasoning in
[`ranking_rationale.md`](ranking_rationale.md); decisions D0063, D0064.

**Reactive docking was adopted, then cut down to size.** It biases docking
*sampling* toward reaction-competent geometry and it works — 20/20 poses put the
warhead carbon 1.55 Å from Cys113's sulfur, in ~2 seconds for 40 runs. Every one
of those poses was **chemically dead**: S–C–Cl median 97.6°, where SN2 needs the
nucleophile anti to the leaving group near 180°. The biasing term is a potential
in distance, isotropic by construction, so it cannot encode an angle at any
parameter value. It is a **sampler, not a criterion** (D0064).

**`shared/nac_criterion.py`** supplies the criterion it cannot: mechanism-specific
approach geometry, windows pre-registered from stereoelectronics before anything
was scored, raw angles always reported so a window can be redrawn without
re-docking. 26 tests, each of three mutations killing exactly the test that names
its claim.

**The result**, 3IKD, 200 runs per molecule, positives = ligands
**crystallographically** bonded to Cys113, negatives = warhead-matched molecules
**measured** inactive:

| class | positives | negatives | AUC | p |
|---|---|---|---|---|
| **chloroacetamide** | 9 | 30 | **0.872** | **0.0004** |
| snar_chloroazine | 2 | 30 | 0.317 | 0.81 |

Positives enrich 2.39× over chance, negatives 0.82×. **This is the first thing
measured on this project that separates actives from inactives** — after five
levels of theory that did not (D0041, D0046, D0036, D0038/D0044, D0057, D0061).

It is one class of two, and it is the class the literature had already converged
on. The SNAr arm failed for a diagnosable reason: its window admits 57% of
negatives against a 29.3% chance baseline, so it is not gating anything.

**Controls were scrutinised rather than trusted**, on @tt8804's instruction. The
HTS's own 34 *actives* were **rejected** as positives — the 11 warhead-bearing
ones are a catalogue of frequent hitters (two rhodanines, an azlactone, an
arylidene barbiturate, an embelin-like quinone…) at 3–75 µM in a 387,000-compound
screen, plus a cephalosporin matching the warhead SMARTS spuriously. Negatives
were re-drawn shuffled, because PubChem file order tracks depositor and therefore
chemical series. And raw viable fractions turned out **not comparable across
mechanisms** — the perpendicular window is 4.4× wider than the SN2 one by solid
angle alone, so an exact isotropic baseline is now divided out.

**Two bugs of mine, both the project's signature defect** — a value taken by
position or name rather than by identity. I keyed pose atoms on the PDBQT *name*
field, where every carbon is named `C`, so each overwrote the last and I measured
the wrong carbon (2.2 Å instead of 1.54 Å). And I guarded reactive typing by
looking for the literal atom type `C1`, where meeko derives it from the base type
— an aromatic carbon becomes `A1` — silently deleting **an entire warhead class**,
30 negatives and 2 positives, by my check rather than by the tool.

---

## Parked

- **guo_pfizer degree-2 enumeration** — running. atra is done (30,000 kept from a
  measured population of 4,063,427, reproducing the earlier run exactly). potter,
  du_xu and liu were dropped to save time.
- **ATRA degree-2 docking** — invalid, it was a D0060 casualty. Deliberately not
  being chased.
- **The NAC ranking build** — designed, costed, not started. See
  `docs/branch_3ikd_receptor.md`.
- **The preparation control** — deposited 3IKD through 6VAJ's path, to split
  receptor from waters in the 2.6× above. One docking run.
