# The shortlist report — outward-facing format

*The format for anything leaving this project: a chemist, a collaborator, a
supervisor. Established 2026-08-10 for the 2.2.0 shortlist and adopted as the
template. Built by [`scripts/shortlist_report.py`](../scripts/shortlist_report.py).*

Build one:

```bash
/data/lab_vm/envs/dwi_cheminf/bin/python3 scripts/shortlist_report.py \
  --candidates t4_2f88a2f534fd t4_54c603efe816 t4_26002bfb953a t4_caf17775e15f \
  --name "T4 shortlist — 100 ns MD"
```

Output lands in `00_outputs/blacksmith/shortlist/` — `shortlist_<N>.html` is the
permanent versioned record, `shortlist.html` the stable name that the next build
overwrites. **Send the versioned one**, or a copy you have checked; the stable
name is not a stable artefact.

---

## 1. What goes in it

Per molecule, in this order:

1. **The identifier**, as a heading.
2. **2D structure** — RDKit via `rdCoordGen`, not `Compute2DCoords`; the default
   layout puts visibly wrong angles on substituted centres.
3. **SMILES in a selectable field.** A `<textarea readonly>` that selects on
   click. Not a code span: the recipient's first action is to paste it into their
   own software, and a span makes them drag-select across a line wrap.
4. **The measured values** — trajectory length, warhead class, 100 ns target
   engagement, mean/max/final ligand RMSD, residence fraction, where it left if it
   left, and the 10 ns sweep readings beside them.
5. **MD movie**, in a panel.
6. **RMSD plots**, in a panel.

## 2. What stays out

**No ranking, no gate verdict, no pipeline commentary.** A reader outside this
project cannot check those claims and does not need them to look at a trajectory.
Sending a rank invites it to be read as a prediction, and every ranking this
project produces is stamped `rank_validated = False` for reasons that take a
retrospective to explain.

Send the molecule and its measurements. Let the chemist do the chemistry.

## 3. Rules the format keeps

- **Self-contained.** One file, opens by double-clicking — no server, no
  directory, no network. That is worth the size.
- **3Dmol.js vendored once** per document, not once per viewer.
- **Every viewer control id namespaced by `elem_id`.** They were bare, which is
  fine for one movie per page and breaks the moment two share a document:
  `getElementById` returns the first match, so all the sliders drive the first
  movie and the rest are inert.
- **Panels start closed.** With four embedded movies this is the difference
  between painting immediately and painting after ~35 MB decodes.
- **Centred on the `body`**, not on each child — centring children individually
  leaves a page of differently-sized blocks reading as left-anchored.
- **4px rules between molecules**, so sections separate rather than scroll.

## 4. The masthead

```
<date> <report name>
<author> · version <X.Y.Z> “<codename>”
```

Nothing else. The version is read through the GUI's own CHANGELOG parser rather
than a second copy of the logic, so a report and the GUI cannot disagree about
which release produced the numbers.

## 5. Size, and getting it out

~9 MB per molecule, because the movie frames are base64 in the page. Four is
~35 MB, ~8.4 MB gzipped.

```bash
gzip -k shortlist.html          # then send the .gz
```

Past most mail limits raw; fine gzipped.

## 6. Before it goes

- **Read it yourself first.** It is leaving the project and no one downstream can
  check it against the tree.
- **Say what the numbers are not.** If a molecule's chemistry carries a caveat the
  recipient would want — a promiscuous warhead class, a control that failed a gate
  — put it in the covering message. The report is deliberately silent on
  interpretation, which means the interpretation has to travel with it.
- **Never send the stable filename's contents unverified.** Rebuilds overwrite it.
