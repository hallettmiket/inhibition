# The GUI — specification

*The results interface, 2.2.0 "Chalcopyrite". @tt8804, 2026-08-09: "this GUI as it
stands right now is the perfect format." This records what that format is, so it
survives the next rebuild.*

Built by [`scripts/mdprio_combine.py`](../scripts/mdprio_combine.py). Everything
below is a rule with a reason; the reasons are the part worth keeping.

---

## 1. Shape

**A selector rail on the left, one viewer on the right.** The rail is the ranking.
Clicking a row loads that molecule's full report into the right pane, replacing
what was there.

Not a row of tabs across the top: that pushed the viewer below the fold and turned
comparison into scrolling.

Panels inside a report **start closed**. A report is ~9 MB of embedded movie
frames and plots; opening everything by default meant the page painted only after
all of it had parsed.

## 2. What ranks, and on what

**The sort key is `explicit_frac_frames_engaged` from the 100 ns MD run** — the
fraction of the run engaging the target.

The 10 ns sweep is **triage**, not a result. It decides which molecules earn a
100 ns run. Ranking on it ranks the selection filter rather than the endpoint,
which is the same defect as ranking on docking energy one stage earlier. It is
carried beside the headline number, in muted type, labelled as a sweep.

On the evidence, the sweep is a weak selector anyway: spearman(sweep, 100 ns
engagement) = **+0.240** over the molecules that have both.

**A molecule with no 100 ns run is not ranked.** It goes in its own band —
*10 ns sweep only — not yet ranked* — ordered by its sweep reading. It does not
get a zero: a zero would read as "measured, engaged nothing", which is a different
and false claim.

## 3. The row

| element | source | note |
|---|---|---|
| 2D structure | `canonical_smiles`, D3/D4 frame or pose sidecar | controls resolve via the sidecar |
| identifier | `parent_ident` | |
| headline number | `explicit_frac_frames_engaged` | shown as `NN% engaged` |
| headline (unranked) | `frac_attack_ready × 10` | shown as `sweep N.NN ns`, muted |
| warhead class | `rank_v2` | |
| max RMSD | `explicit_ligand_rmsd_nm_max` | |
| verdict tag | `rmsd_max < 1.2 nm` | `held` / `left`, or `swept` |
| bar | engagement | |

A **legend sits at the top of the rail** stating the ranking axis in words. A
number in a column has to be readable without knowing which of two similar
quantities it is; the same figure in ns previously meant the sweep, and nothing
said so.

## 4. Toggles, and why each exists

Each encodes a measured distinction, not a preference.

- **all classes / by warhead class** — cross-class ranking is biased because the
  S<sub>N</sub>2 angular criterion is far stricter than the perpendicular one
  (#47). The toggle makes that visible instead of something the reader must
  remember.
- **combined / split held-left** — engagement and residence are near-independent
  (ρ = −0.007, #46). A molecule can rank high and still leave the pocket.
- **controls** — the crystallographic controls, alone, with the controls page in
  the viewer.
- **light / dark** — stamps the shell *and* the framed report, and persists.

## 5. Controls

Controls appear **twice**: ranked among the candidates, and in their own tab.

Two kinds, treated differently because they are different experiments:

- **`rx_*` (reactant)** — the crystal pose with the covalent bond cleaved and the
  leaving group rebuilt; every other atom keeps its crystallographic coordinate.
  This is the non-covalent pre-reaction complex, it sweeps, and it is directly
  comparable to a candidate.
- **`xtal_*` (bonded)** — the deposited geometry, still bonded to SG at ~2 Å. That
  is the reaction *product*, and the near-attack window is 2.8–4.2 Å, so it cannot
  be attack-ready **by construction**. No number, and the reason is shown in its
  place.

**Never drop a control that produced no number.** "Could not be swept" and "was
never run" must not look the same.

**No control carries a held/left tag.** None has a 100 ns trajectory, so the
verdict does not exist for them.

A control with its own report opens in the **identical viewer** as a candidate —
pose, movie, RMSD plots. Only controls without a report fall back to the controls
page.

## 6. Honesty rules

These are the ones that will be quietly broken first.

1. **Never fabricate a value to make something sortable.** If a molecule has no
   measurement on the ranking axis, it leaves the ranking — it does not get a
   placeholder.
2. **Stamp, do not drop.** "Did not run", "ran and failed", and "ran and scored
   badly" are three different facts and must stay distinguishable.
3. **A rendering that describes a different frame must say so.** The MD surface is
   a mesh, rebuilt on release rather than per frame; while it lags it is labelled
   `surface: frame N`.
4. **Every ranking is stamped `rank_validated = False`** and is described as an
   ordering the pipeline produced, not as evidence of binding.
5. **Illustrative graphics are labelled in their own masthead**, not in a
   footnote. The pipeline schematic's geometry is drawn; its parameters are real;
   the page says both.

## 7. Top bar

One strip, ~38 px: title, **version + codename**, the toggles, the hint text, the
schematic link, the theme button. Everything below it is the work, and a two-line
masthead over a scrolling list is space the reader never gets back.

The version is **parsed from `CHANGELOG.md`**, never written here as a literal. A
version constant in a builder is a pin, and pins in this repo go stale silently.

## 8. The schematic

`how this works ↗` opens [`pipeline_schematic`](../shared/pipeline_schematic.py):
dock → split → criteria → rank → sweep → MD, in short sentences, each paired with
a diagram, plus **one interactive panel of real output** (491 real poses of one
screened molecule, coloured by the mode the screen assigned).

Seeded, so rebuilds are byte-identical — an unseeded page would produce a new
versioned file on every run differing only in dot positions.

## 9. Where it is served

A plain `http.server` over the report directory. Reads are concurrent-safe and
nothing is written, so several people can point at one instance:

```bash
ssh -L 8931:127.0.0.1:8931 <you>@<host>     # then http://localhost:8931/combined.html
```
