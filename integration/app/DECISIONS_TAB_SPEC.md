# GUI spec — the Decisions pane

**For the artist.** This is the "one place to look" answer to a real problem:
provenance is currently spread across decision records, run manifests, runbooks,
and reference provenance, and nobody should grep four formats to learn why a
docking box is 26 Å.

**The GUI reads; it does not own.** Every panel below is a rendering of files in
the repo or the append-only tree (see `decisions/D0008`). Nothing here is the
source of truth, and the app must be safe to close and rebuild at any time.

## Data sources

| Panel input | Provided by |
|---|---|
| decision records | `shared.decisions.load()` → `Decision.to_dict()` |
| what a run consumed | `manifest.json` in `append_only/inhibition/<stage>/` |
| how to redo a judgment | `docs/runbooks/*.md` |
| reference-data origin | `data/reference/.provenance.md` |

`shared.decisions` already returns GUI-shaped dicts — frontmatter fields plus
the `context` / `decision` / `consequences` sections parsed out. No markdown
parsing belongs in the app.

## Panel 1 — Choreography decision log

All records, newest first. Columns: `id`, `title`, `approach`, `origin`,
`status`.

- Filter by `approach` (shared / t1–t4 / integration) and by `origin`
  (spec / adversary / implementation / user). **`origin: adversary` is worth
  surfacing prominently** — those records are the audit trail that the
  adversarial review actually changed the design, which is one of the project's
  stated goals.
- Superseded records are shown struck through, not hidden. *Why the answer
  changed* is usually more informative than the current answer.
- Expanding a row shows Context / Decision / Consequences, the evidence list,
  and links to the governing runbook and affected files.

## Panel 2 — Per-approach tab

On each approach's tab, `shared.decisions.by_approach("t4")` — the approach's
own records **plus the shared ones it inherits**. A reader on the T_4 tab needs
to know about the receptor and box decisions without going hunting.

Render the evidence list verbatim. It is deliberately numeric (`r = 0.396`,
`1.78 Å`, `13.6× span`) rather than adjectival, and those numbers are the
argument.

## Panel 3 — "Why is this file like this?"

A path box calling `shared.decisions.affecting(fragment)`. Typing
`warhead_classes` returns D0006; `receptor.yaml` returns D0001 + D0002. This is
the panel that replaces grepping.

## Panel 4 — Run provenance

For a selected approach, read its `manifest.json` files and show: run id, git
commit **and whether the tree was dirty**, config hashes, input/output hashes,
tool versions.

Surface `git.dirty == true` as a visible warning. It means the commit hash does
not fully describe the code that ran, and the outputs are provisional — that
must not be discoverable only by reading JSON.

## Panel 5 — Open questions

Records with `status: proposed`, plus known-blocked items pulled from source
status fields: sources marked `pending` in `sources.lock.json`, warhead classes
at `NEEDS_DESIGN` or `UNVERIFIED`, and reference rows marked `UNVERIFIED`.

The choreography's honest limits belong on screen next to its results, not in a
file nobody opens. As of writing that list includes: the CReM fragment DB radius
variant unchosen, the decoy set unbuilt, BDHI and naphthoquinone attachment
regiochemistry undesigned, and the Byun BDHI fragment structure unresolved.

## Non-goals

- **No editing in v1.** Records are authored in the repo and reviewed in a PR.
  If lab members later need to author decisions in the app, it must write
  *through* to `decisions/*.md` — never into a database that shadows them.
- **No scoring or ranking here.** This pane explains; the candidate panes rank.
