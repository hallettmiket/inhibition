# The final GUI — an outline for edit

*@twu383 with Claude Code, 2026-08-25. Written to be marked up, not adopted.
Everything here is a proposal except §1, which is a measurement.*

You asked for an outline of how the final GUI should look, on the basis that
"we may need to update our current gui to show all 4 approaches". Before
proposing anything, §1 records what the two interfaces actually show today,
because the answer changes which of them needs the work.

---

## 1. What exists today — measured, 2026-08-25

**There are two GUIs, and only one of them is missing approaches.**

| | integration GUI | MD-priority GUI |
|---|---|---|
| built by | `integration/app/app.py` (Streamlit) | `scripts/mdprio_combine.py` (static HTML) |
| served by | `scripts/serve_gui.sh` (default :8899) | `scripts/serve_gui.sh`-style plain server, :8931 |
| what it ranks | four shortlists, each on its own metric | one rail, ranked on `explicit_ligand_rmsd_nm_max` |
| **approaches covered** | **all four** — T₁ 25, T₂ 25, T₃ 25, T₄ 27 | **T₄ only** |
| spec | `docs/gui_spec.md` §1–§9 | same file |

Verified by loading the data layer directly rather than by reading the page:

```
t1   T₁ · de novo (DiffSBDD)              rows=25  synth_avail=True
t2   T₂ · CReM neighbourhood · seed ATRA   rows=25  synth_avail=True
t3   T₃ · decoration (LibInvent)           rows=25  synth_avail=True
t4   T₄ · combinatorial                    rows=27  synth_avail=True
```

and, for the 3.1.0 report directory `mdprio_reports_nac_v5/`, by counting what
is in it: **2,054 `combined_*.html` fragments, 15 `t4_*` molecule reports, one
controls page — and no `t1_`, `t2_` or `t3_` report at all.** Every per-molecule
CSV in that directory is `mdprio_t4_*`.

**So the integration GUI is not the problem.** It has shown four approaches
side by side since 2.x, and it is the interface the manuscript describes in
§1.2 ("presented side by side … see structural convergence"). It is what
Figure 2 captures.

**The 3.1.0 MD-priority GUI is T₄-only**, and that is the gap. It is also the
one you have been looking at on :8931, which I think is where the question came
from.

### Why it is T₄-only, which matters for how much work this is

Not an oversight in the interface. The 100 ns campaign selected molecules by
near-attack-geometry triage, which is a **covalent** criterion — it measures
whether a warhead can present itself to Cys113 SG. T₁ and T₂ are non-covalent
arms with no warhead, so they were never eligible for that queue, and T₃ (also
covalent) had modes in earlier `nac_v3` runs but none in `nac_v5`.

**Adding T₁/T₂ to that rail is therefore not a display change. It is a decision
about what the rail ranks**, and it needs a criterion the non-covalent arms can
satisfy. That is §3 below and it is the part I most want you to rule on.

---

## 2. The shape I would propose

**One interface, three levels, and the level boundary is the honesty boundary.**

The current split is accidental — two GUIs because two campaigns — and it costs
a reader the thing the choreography is for: seeing four independent answers to
one question next to each other, with the evidence behind each.

```
┌────────────────────────────────────────────────────────────────┐
│ LEVEL 1 — THE FOUR SHORTLISTS          (integration GUI today) │
│ four columns, each ranked on its own metric, never merged      │
│ every column carries its gate verdict in the column header     │
└──────────────────────────────┬─────────────────────────────────┘
                               │ click a candidate
┌──────────────────────────────▼─────────────────────────────────┐
│ LEVEL 2 — THE CANDIDATE DOSSIER        (integration GUI today) │
│ structure · SMILES · rank · shared descriptor axes             │
│ two independent free-energy estimates · gate verdict           │
│ + NEW: "what happened to this molecule downstream"             │
└──────────────────────────────┬─────────────────────────────────┘
                               │ if it earned a trajectory
┌──────────────────────────────▼─────────────────────────────────┐
│ LEVEL 3 — THE TRAJECTORY REPORT        (MD-priority GUI today) │
│ pose · movie · RMSD per replicate · held/left verdict           │
│ today reachable only from the T₄ rail on :8931                 │
└────────────────────────────────────────────────────────────────┘
```

**The single highest-value change is the link from level 2 to level 3.** Right
now a molecule that went through a 100 ns run appears in *both* GUIs with no
route between them, and the two disagree about what "rank" means. A reader who
finds `t4_c24c106bd005` in the shortlist has no way to learn that it has a
trajectory, and a reader on the rail has no way back to the shortlist that
produced it.

### What each level owns

- **Level 1 owns comparison across approaches**, and owes the reader the
  reminder that the columns are incommensurable. It must never sort across
  columns.
- **Level 2 owns one molecule's whole evidence file**, and is the only level
  where quantities from different stages sit together. It therefore carries the
  heaviest labelling burden: every number needs its stage and its units.
- **Level 3 owns the endpoint**, and is the only level entitled to say "held"
  or "left".

---

## 3. The question I cannot answer for you

**What should the level-3 rail rank when a non-covalent molecule is on it?**

Three options, and they are not equivalent:

| | option | what it costs |
|---|---|---|
| **A** | **Keep the rail covalent-only, and say so in the masthead.** Level 3 stays T₃/T₄; T₁ and T₂ stop at level 2. | Honest and free. But the interface then never shows the two non-covalent arms at their endpoint, and the paper's "four approaches" claim is only ever demonstrated at levels 1–2. |
| **B** | **Rank on ligand RMSD alone**, which every arm can produce, and show near-attack geometry as a covalent-only column. | One rail for all four. But RMSD-only ranking drops the mechanistic argument that is precisely why T₃/T₄ are worth more than occupancy (manuscript §1.9), so the rail would rank a well-behaved non-covalent binder above a covalent one that engaged. |
| **C** | **Four bands, one per approach, each ranked on its own endpoint metric**, never cross-sorted — the level-1 rule applied at level 3. | Consistent with everything else in the interface, and the most defensible. Costs the most work, and means the rail can no longer answer "what is the single best molecule", which may be the question you actually want it to answer. |

**I would take C**, because it is the only one that does not either hide two
arms or invent a comparison the project has spent 87 decision records refusing
to invent. But B is much cheaper and A is free, and if the rail's job is triage
for your own next experiment rather than presentation, A may be right.

**This is the decision I would like marked up.** It should become a decision
record either way, because it binds what the interface can claim.

---

## 4. Specific changes, in the order I would make them

Each is small enough to do independently. **1–3 are worth doing whatever you
decide in §3.**

### 1. Join the two GUIs at the candidate

On the level-2 dossier, add a row: *"this candidate has a 100 ns trajectory —
open its report"*, linking to the level-3 page. Where there is none, say which:
*never triaged* / *triaged, did not earn a run* / *ran and failed*. Three
different facts (`gui_spec.md` §6 honesty rule 2), and today all three look the
same, which is absence.

### 2. Make the adduct collapse visible in T₄

T₄'s shortlist shows nine warhead classes, and **three of them —
chloroacetamide, sulfamate acetamide and sulfonate acetamide — collapse to the
same `acetamide_adduct` and therefore dock as literally the same species.**
Their scores are identical to five decimal places:

```
t4_617de1a16274  chloroacetamide       -5.43626
t4_97020e8242c4  sulfamate_acetamide   -5.43626
t4_0a52a08197fc  sulfonate_acetamide   -5.43626
   adduct_smiles for all three: CC(=O)N(c1ccccc1)[C@@H]1CCS(=O)(=O)C1
```

This is **correct** — the leaving group departs, so the docked adduct is shared
— but it reads as a bug, and I checked it as one before confirming otherwise.
Group those three under their adduct class, or add an `adduct_class` column, so
the reader is told rather than left to wonder. **Nine warhead classes are seven
distinct docked species**, and any count of "classes explored" that says nine is
overstating the search by two.

### 3. State the receptor in the masthead

`config/receptor.yaml` still pins 6VAJ while the benchmark paths guard for
3IKD, and which one you get depends on the entry point (`state_of_the_project.md`
§1, §8). The interface should say which receptor produced what is on screen. It
is one line and it closes a live instance of the project's own failure mode.

### 4. Carry D0088's caveat onto every per-mode number

The modes in `nac_v5` are mixtures — the splitter clustered on the score and
then scored the clusters. `viable_fraction`, `enrichment` and `conditional_eb`
are all measured over those mixtures. The 100 ns results are unaffected as
measurements, but **which** pose earned each one was chosen by that machinery.
Every per-mode number on level 3 should carry that, the way level 1 carries the
gate verdict. Until it does, the rail is presenting a number the project has
already recorded as compromised.

### 5. Then, whichever of §3 you chose

---

## 5. What not to change

From `gui_spec.md`, and worth restating because a rebuild is where these die:

- **The rail on the left, one viewer on the right.** Tabs across the top pushed
  the viewer below the fold.
- **Panels start closed.** A report is ~9 MB of embedded frames.
- **Never fabricate a value to make something sortable.** A molecule with no
  measurement on the ranking axis leaves the ranking; it does not get a zero.
- **Every ranking stays stamped `rank_validated = False`.** Five levels of
  theory have now failed on this pocket; nothing on screen may read as evidence
  of binding.
- **Controls appear twice** and `rx_*` and `xtal_*` stay distinguishable — the
  bonded crystal geometry cannot be attack-ready by construction, and must not
  be shown as though it scored badly.

---

## 6. What this outline does not cover

- **Whether the shortlists should be re-run.** D0088 makes the mode assignment
  circular, and #79 says adopting `pose_cluster.py` means a full re-screen. If
  that happens, level 1's contents change and this outline still holds.
- **Performance.** The level-3 directory is already 2,174 files; a four-approach
  version is roughly 4× that, and the current server is a plain `http.server`.
- **Anything about who may read it.** The integration GUI runs as the data owner
  precisely so a viewer needs no permissions on the tree, and that should
  survive any rebuild — the Isilon ACL means several project members cannot read
  a single frame directly.
