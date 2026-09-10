# How the Pin1 manuscript figures are made

*@twu383 with Claude Code, 2026-09-01. Compact: what runs, in what order,
and the four things that will waste your afternoon if you don't know them.*

Two figures ship. Both are regenerable; neither has a number typed into it.

| | file | script |
|---|---|---|
| main text, full page | `pin1_main_figure.{pdf,png}` + `_caption.tex` | `scripts/manuscript_fullpage_figure.py` |
| supplementary | `pin1_si_figure.{pdf,png}` + `_caption.tex` | `scripts/manuscript_si_figure.py` |

Output lands under `append_only/inhibition/00_outputs/artist/manuscript_figures/`,
versioned by `shared/outputs.py`. Copy the highest integer of each stem into
`murmurent_manuscript/figures/pin1/` under the un-numbered name.

---

## The order

**1. Serve the GUI.**

```bash
bash scripts/serve_gui.sh 8907          # integration app, loopback only
```

**2. Capture it.** Headless chromium via playwright in a throwaway venv
(`pip install playwright && playwright install chromium`). Three shots:

- the Shortlists panel, sidebar collapsed, viewport 1680×1500 @ dsf 2.5
- the Candidate dossier with T₄ selected — pick the approach by **label**
  (`get_by_label("approach")`), not by index; with the sidebar open its
  selectbox is index 0 and you get the T₂ seed picker instead
- nothing else; panel D comes from step 3

**3. Render the pose.**

```bash
python scripts/pose_interaction_render.py t4_c24c106bd005 pose.html 0.55 22
```

Computes contacts (≤ 4.5 Å) and polar contacts (≤ 3.6 Å) from the pose and the
receptor, then writes a 3Dmol page: protein as a sub-pocket-coloured surface,
ligand in sticks, contacting residues labelled. Screenshot it, then crop to the
content bounding box — that is what centres it.

**4. Build.** Both scripts parse their own counts at render time from
`decisions/*.md`, `docs/how_this_project_breaks.md` and the enrichment-gate
token, print what they parsed, and emit the caption as a separate `.tex`
fragment. The caption is **not** baked into the artwork: the manuscript is
Overleaf-synced and a caption inside a PDF cannot be edited there.

---

## The four things that will cost you time

**The receptor defaults to the wrong one.** `pose3d.receptor_for()` returns
`6VAJ_prepared.pdb`, and every frame also carries a `pose_path` from that
superseded production run. D0059 replaced it with 3IKD. Both files exist and
both render a perfectly plausible picture. Use `nac3_pose_path` and the
receptor registered against it.

**Two rank columns disagree.** `rank_synth` is what the shortlist shows; raw
`affinity_kcal` gives a different order, because the shortlist ranks the
size-decorrelated residual (D0049). For `t4_c24c106bd005` that is rank 1 versus
10th of 187. Say which one you mean.

**The headless render goes blank, nondeterministically, and it is sticky per
page load.** Retrying inside one load fails 8/8; two near-identical scenes go
opposite ways. A *fourth* transparent surface fails every time, and full opacity
(1.0) fails while 0.85 is fine. Fix: discard the GL context between attempts —
new browser, new page — and keep taking shots until the frame has ink. This is
a software-GL limitation; on a real GPU interactively none of it bites.

**Fonts substitute silently.** Panel letters must be genuine Times New Roman
(`~/.fonts/TimesNewRoman-{Regular,Bold}.ttf`). If they are missing matplotlib
falls back to the next serif without a word and the figure looks right. Check
the *embedded* name, not the render:

```bash
python -c "import pymupdf; d=pymupdf.open('pin1_main_figure.pdf'); \
  print({x[3] for p in d for x in p.get_fonts(full=True)})"
# want TimesNewRomanPSMT / TimesNewRomanPS-BoldMT, not NimbusRoman
```

---

## Conventions the figures hold to

- Panel letters `A.` `B.` `C.` `D.`, Times New Roman 11 pt bold, set **above**
  the panel, aligned on the left margin.
- No figure title, no annotation arrows, no descriptive text inside the
  artwork. All of it lives in the caption.
- Captions are written for a bioscientist. No "spec", no "guard", no
  "selection by where it sat" — see the mapping in `manuscript_si_figure.py`.
- Categorical colours are a fixed CVD-validated order
  (`#332288 #DDAA33 #44AA99 #CC3311`, worst-case OKLab ΔE 13.1). Re-run the
  check before adding a hue.
- Counts are never typed. If a figure quotes a number, it parsed it.

---

## Every path

Absolute, so nothing depends on where you are standing.

### Scripts — `~/repos/inhibition/scripts/`

| | |
|---|---|
| `/home/UWO/twu383/repos/inhibition/scripts/manuscript_fullpage_figure.py` | builds the main figure |
| `/home/UWO/twu383/repos/inhibition/scripts/manuscript_si_figure.py` | builds the supplementary figure + its caption |
| `/home/UWO/twu383/repos/inhibition/scripts/pose_interaction_render.py` | computes contacts, writes the 3Dmol page for panel D |
| `/home/UWO/twu383/repos/inhibition/scripts/capture_3dmol_page.py` | screenshots a 3Dmol page, retrying on a fresh GL context |
| `/home/UWO/twu383/repos/inhibition/scripts/serve_gui.sh` | serves the integration GUI on loopback |

Retired — still run, nothing ships from them:
`manuscript_figures.py`, `manuscript_gui_figure.py`,
`manuscript_figure_combined.py`, `manuscript_figures_final.py`, same directory.

### Generated output — governed, versioned, append-only

`/data/lab_vm/append_only/inhibition/00_outputs/artist/manuscript_figures/`

Current highest version of each shipping stem:

```
pin1_main_figure_4.pdf          pin1_main_figure_4.png
pin1_main_figure_caption_1.tex
pin1_si_figure_3.pdf            pin1_si_figure_3.png
pin1_si_figure_caption_3.tex
```

### What the manuscript actually includes

`/home/UWO/twu383/repos/murmurent_manuscript/figures/pin1/`

```
pin1_main_figure.pdf   pin1_main_figure.png   pin1_main_figure_caption.tex
pin1_si_figure.pdf     pin1_si_figure.png     pin1_si_figure_caption.tex
pin1_fig2_gui_output.pdf/.png   <- superseded; kept only so main-article.tex
                                   line 720 still builds until it is repointed
CHANGES_IMPLIED.md
```

Panel A's source, Hagar's flowchart:
`/home/UWO/twu383/repos/murmurent_manuscript/figures/inhibition_flowchart_v3-2.pdf`

### Fonts — per-user, and required

```
/home/UWO/twu383/.fonts/TimesNewRoman-Regular.ttf
/home/UWO/twu383/.fonts/TimesNewRoman-Bold.ttf
```

After adding or replacing them, `fc-cache -f` and delete matplotlib's
`fontlist-*.json` under `python -c "import matplotlib;print(matplotlib.get_cachedir())"`.

### Intermediate captures — EPHEMERAL, regenerate them

The GUI captures and the rendered pose live in this session's scratchpad and
**will be deleted**:

```
<scratchpad>/shots/f2_shortlists_wide.png     full Shortlists panel capture
<scratchpad>/shots/f_dossier_t4.png           full Candidate-dossier capture
<scratchpad>/shots/flowchart.png              flowchart PDF at 420 dpi
<scratchpad>/shots/panel_shortlists.png       crop fed to panel B
<scratchpad>/shots/panel_candidate.png        composite fed to panel C
<scratchpad>/shots/pose_interaction.png       cropped render fed to panel D
```

Nothing is lost — steps 1–3 above rebuild all of them, and the build scripts
take the containing directory as their only argument. If you want them to
persist, write them under
`/data/lab_vm/append_only/inhibition/00_outputs/artist/` instead.

### Read alongside

```
/home/UWO/twu383/repos/inhibition/docs/manuscript_changes_2026-08-25.md
/home/UWO/twu383/repos/inhibition/docs/gui_final_outline.md
/home/UWO/twu383/repos/inhibition/data/ready_to_delete.md
```

Issue of record: https://github.com/hallettmiket/murmurent_manuscript/issues/1

## Retired

`manuscript_figures.py`, `manuscript_gui_figure.py`,
`manuscript_figure_combined.py` and `manuscript_figures_final.py` built the
earlier five-figure and two-figure sets. They still run; nothing ships from
them. Superseded renders are listed in `data/ready_to_delete.md`.
