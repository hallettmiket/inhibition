# 2.1.0 “Bornite” retrospective

*Written at the close of the 2.1.0 ranking framework, 2026-08-07, before 2.2.0
opens. Companion to `docs/framework_2.1.0.md` (what was built) and
`docs/outline_2.2.0.md` (what comes next).*

**The one-line version: 2.1.0 rebuilt the ranking on a sound structure and then
discovered, on the last day, that the score inside that structure rests on an
ordering carrying no information. The structure is keepable. The score is not.**

---

## 1. What was built

| | |
|---|---|
| **v2 screen** | all 5,769 candidates re-docked, **5,765 ok**, 4 failed. Per-pose geometry, gnina scores and **the poses themselves** persisted — 115,300 pose rows and 5,765 SDFs |
| **Ranking** | consensus filter as a **quota within warhead class**, then a weighted anchoring score. One list per class, no global top-N |
| **Selection** | walks the ranked list, **re-measures the pose it is about to elevate** rather than re-reading the score |
| **Pose ranking** | new stage: BPMD every reaction-competent pose, elevate the most stable |
| **Elevation** | 100 ns non-covalent MD, then the #22 gate, then covalent |
| **References** | 25 known-binder jobs through the identical criterion; Sulfopin and ATRA through 100 ns |
| **Tooling** | gnina CNN re-ranking adopted on measurement; 3IKD_ian pinned by digest; one house report style |
| **GUI** | a Ranking 2.1.0 panel with references beside candidates, click-to-view, independent scrolling |

## 2. What 2.1.0 got right

**The screen persists its working.** 2.0.0 computed per-pose geometry, reduced it
to one number and threw the poses away — which is why two known repairs were both
blocked on the same missing data and needed a full re-dock. 2.1.0 keeps
everything, and that decision paid within a day: every analysis in §3 below was
possible only because the poses were on disk.

**Measure on this target, do not inherit.** gnina was adopted because it lifted
top-1 pose accuracy from 18.3% to 26.8% *on 3IKD over 82 Pin1 crystal ligands*,
not because a published benchmark liked it. The same discipline then produced the
finding that killed the score (§4) — the instrument that validated the re-scorer
also indicted the ordering.

**Separating the two ranking problems.** @tt8804 caught that ranking *molecules*
and ranking *a molecule's own poses* were being answered by the same machinery.
They are different questions and now have different stages. The dry run
vindicated it immediately: one molecule had three viable poses at 12.0°, 4.4° and
22.5° off ideal, and nothing before that stage chose between them on evidence.

**References as a yardstick, not a gate.** @tt8804's ruling — 15 crystallographic
depositions cannot decide which chemistry to pursue — removed a validation gate
that would have narrowed the shortlist to chemistry nobody wanted. Keeping the
references *on screen* while refusing to *calibrate* against them is the
distinction that made both possible.

**The kinetic argument for the residence gate.** Rather than calibrating ≥90%
residence against reference molecules, the threshold was derived: residence is
1/k_off, so even a millimolar binder should sit through 100 ns by four orders of
magnitude. Sulfopin and ATRA then both held for the full trajectory and the lead
left at 54 ns — the prediction and the measurement agreed, and the gate needed no
tuning.

## 3. What went wrong

### 3.1 The recurring defect, still recurring

Every significant error was **a value taken by position, name, or inheritance
rather than by identity** — the same class `recap_2.0.0.md` §3 named. Nine
instances in 2.1.0 alone:

| defect | the value taken wrongly | cost |
|---|---|---|
| **selection queue** | *newest file* rather than the queue with molecules in it | **17 candidates stranded; zero elevation overnight** |
| **ranking file** | *newest file* rather than the score by name | selection ran off `enrichment_joint`, the 2.0.0 quantity |
| `app.py` splice | region computed as `s[:a]+new+s[b:]` with **a > b** | every added function silently duplicated |
| `pose3d` | unregistered pose column **fell back to 6VAJ** | would have drawn every pose 48.6 Å from the pocket |
| `nac.measure` | argument **position** — called as `(sg, coords, mechanism)` | caught by the tool erroring |
| Cys113 SG | a residue **number** that differs between pipeline stages (113 vs 63) | caught by an assertion |
| reference classes | a prose **label** rather than a SMARTS match | would have mis-assigned 20 of 22 |
| BPMD resume | an incomplete **key** (ident+replicate, no trajectory length) | would have seated a 300 ps replica in a 3 ns comparison |
| reference table | no `smiles` — structures **absent** rather than wrong | the comparison the panel exists for could not render |

**Two of these reached the user.** The rest were caught by the repo's own guards,
by an assertion written the same hour, or by a test. That ratio is the argument
for the guards.

### 3.2 Silent failures beat loud ones

The three that cost the most all failed **quietly**:

- The overnight chain **exited 0 on every stage it skipped.** `pose ranking exit
  1` and `elevation launcher exit 1` were logged and nothing acted on them.
- `st.dataframe(on_select=...)` is **inert when handed a Styler** — no warning, no
  error, a table that renders perfectly and never selects.
- The `app.py` duplication left the file working, because Python takes the later
  definition.

The lesson is narrower than "add more tests": **a stage that produces no output
should fail the run, not log a non-zero exit and continue.** The chain had the
information and carried on.

### 3.3 Verification I did not do

I repeatedly verified that code *ran* and not that it *worked*. The Styler bug
shipped because I checked the panel rendered and the data path had no duplicate
columns — never that a click did anything, which I cannot do from here. The right
response is not to guess better; it is to say plainly which parts are unverified.
Where I did that (the movie fix), it was useful; where I did not, it wasted a
round trip.

## 4. The finding that ends the version

Prompted by @tt8804 asking *how were these poses selected*.

**They are the top 20 by AutoDock docking energy — and energy carries no
information about whether the warhead is aimed at Cys113.**

| top-10 window ordered by | mean viable fraction | zero viable |
|---|---:|---:|
| AutoDock energy *(what the score uses)* | 0.181 | 26.4% |
| gnina CNNscore | 0.183 | 24.2% |
| all 20, unordered | 0.185 | 0% |

ρ(energy_rank, viable) = **+0.009** across 115,300 poses — and *positive*, so
worse-scoring poses are marginally more reaction-competent. **42.7% of molecules
have viable poses among their 200 and none in their top 10.**

**Sulfopin — the parent, with a crystallographic Cys113 adduct — scores 0.000.**
34 of its 200 poses are reaction-competent; none of its 20 lowest-energy are.
Both Reddi-2023 sulfamates, more potent in cells than sulfopin, sit near the
floor too.

So `topn_viable_frac`, which I introduced and made the default, inherits an
ordering that is noise for the thing it counts. The top-N *idea* stands — D0068's
dilution argument is untouched — but the window was built on the wrong quantity.
Recorded as **issue #23**.

**What this says about the version.** The framework caught its own central error,
using an instrument the same version built, on a molecule chosen precisely
because the answer was already known. That is the system working. It is also,
plainly, a score that shipped into a GUI before it was tested against a known
positive — and testing against Sulfopin was a day's work available from the
start.

## 5. What carries into 2.2.0

**Keep:** the persisted screen; per-class quotas; the separation of molecule and
pose ranking; references as yardstick; the kinetically-derived gate; measure-on-
this-target discipline; the house report style.

**Fix:** the primary score (#23); `anchor_quality`, which is averaged over the
same energy window and is exposed to the same problem; `mmgbsa.RECEPTOR_PDB`
still defaulting to 6VAJ; a chain that continues past a failed stage.

**Test:** whether conditional enrichment converges at 200 vs 2,000 runs — the
check D0068 forces on any candidate score, and one no score has yet passed.

**Unresolved:** whether the SN2 150° threshold is too strict. Sulfopin clears it
34 times in 200 poses but never among its best-scoring, which is a question for a
chemist rather than a statistician.

---

## 6. Two numbers worth remembering

**41.5%** — how often docking finds a sub-2 Å pose for a Pin1 crystal ligand.
**18.3%** — how often it puts that pose first.

Almost everything 2.1.0 built sits in the gap between those, and almost every
mistake came from forgetting the gap was there.
