"""
Purpose: Generate the derived MkDocs pages from the repo's source-of-truth files.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: decisions/*.md, data/reference/*.csv, config/sources.lock.json
Output: virtual pages under the mkdocs docs tree (via mkdocs-gen-files)

Run automatically by mkdocs through the gen-files plugin; never invoked by hand.

WHY GENERATE. decisions/ and data/reference/ are the source of truth (see
decision D0008). Hand-copying them into the docs tree would create a second
copy that drifts, and a drifted decision log is worse than none. So the site
renders them at build time and the generated pages carry a banner saying so.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import mkdocs_gen_files

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

BANNER = (
    "!!! info \"Generated page\"\n"
    "    Rendered at build time from the repo's source of truth. "
    "Edit the underlying file, not this page.\n\n"
)

STATUS_BADGE = {
    "accepted": ":material-check: accepted",
    "superseded": ":material-history: superseded",
    "proposed": ":material-help: proposed",
    "rejected": ":material-close: rejected",
}


def _decisions():
    """Load decision records, tolerating an environment without RDKit."""
    from shared import decisions as dec  # noqa: PLC0415 - deferred for import cost
    return dec.load(REPO / "decisions")


def gen_decisions_index() -> None:
    """One page listing every decision, grouped by approach."""
    try:
        records = _decisions()
    except Exception as exc:  # noqa: BLE001 - a broken record must not kill the build
        with mkdocs_gen_files.open("decisions/index.md", "w") as f:
            f.write(f"# Decisions\n\n!!! danger \"Could not load\"\n    {exc}\n")
        return

    order = ["shared", "t1", "t2", "t3", "t4", "integration"]
    label = {"shared": "Shared substrate", "t1": "T_1 de novo",
             "t2": "T_2 ATRA neighborhood", "t3": "T_3 REINVENT R-groups",
             "t4": "T_4 combinatorial", "integration": "Integration"}

    with mkdocs_gen_files.open("decisions/index.md", "w") as f:
        f.write("# Decisions\n\n" + BANNER)
        f.write(
            f"{len(records)} record(s). Each answers *why* a choice was made; the "
            "[runbooks](../runbooks/index.md) answer *how* to make that kind of "
            "choice again, and a run's `manifest.json` records *what* it "
            "actually consumed.\n\n"
            "Records marked `origin: adversary` are decisions the adversarial "
            "review forced — the audit trail that review changed the design.\n\n"
        )
        for ap in order:
            rows = [r for r in records if r.approach == ap]
            if not rows:
                continue
            f.write(f"## {label[ap]}\n\n")
            for r in rows:
                f.write(f"### {r.id} — {r.title}\n\n")
                f.write(f"**{STATUS_BADGE.get(r.status, r.status)}** · "
                        f"`origin: {r.origin}` · {r.date}")
                if r.superseded_by:
                    f.write(f" · superseded by **{r.superseded_by}**")
                f.write("\n\n")
                if ctx := r.section("Context"):
                    f.write(f"{ctx}\n\n")
                if dcn := r.section("Decision"):
                    f.write(f"**Decision.** {dcn}\n\n")
                if csq := r.section("Consequences"):
                    f.write(f"**Consequences.** {csq}\n\n")
                if r.evidence:
                    f.write("??? note \"Evidence\"\n")
                    for e in r.evidence:
                        f.write(f"    - {e}\n")
                    f.write("\n")
                if r.affects:
                    f.write("Affects: " + ", ".join(f"`{a}`" for a in r.affects) + "\n\n")
                if r.runbook and r.runbook != "null":
                    f.write(f"Runbook: `{r.runbook}`\n\n")
                f.write("---\n\n")


def gen_reference_page() -> None:
    """Render the frozen reference set as browsable tables."""
    import csv

    ref = REPO / "data" / "reference"
    with mkdocs_gen_files.open("shared/reference_set.md", "w") as f:
        f.write("# The frozen reference set\n\n" + BANNER)
        f.write(
            "The single source for two things and nothing else: the **novelty "
            "axis** for every approach (`1 - max Tanimoto ECFP4`, computed "
            "against this set and **never against the seed**), and **T_4's "
            "reactivity window**.\n\n"
        )
        for name, title, note in [
            ("pin1_reference_binders_1.csv", "Master set — novelty axis",
             "Rows marked `UNVERIFIED` are excluded from the novelty computation."),
            ("pin1_covalent_cys113_anchors_2.csv", "Covalent Cys113 anchors — reactivity window",
             "`reference_set.py` refuses `UNVERIFIED` rows into the window."),
            ("warhead_classes_2.csv", "Warhead classes — T_4 enumeration",
             "`warhead_library.enumerable()` defaults to `VERIFIED` only."),
            ("pin1_reactivity_kinetics_1.csv", "Measured reactivity kinetics",
             "Digitized from a figure to ~1 significant figure. Bound a window with "
             "these; do not treat them as precise."),
        ]:
            p = ref / name
            if not p.is_file():
                continue
            f.write(f"## {title}\n\n`data/reference/{name}`\n\n{note}\n\n")
            rows = list(csv.reader(p.read_text(encoding="utf-8").splitlines()))
            if not rows:
                continue
            head, body = rows[0], rows[1:]
            keep = [i for i, h in enumerate(head)
                    if h not in ("citation", "structure_source", "notes", "source")]
            f.write("| " + " | ".join(head[i] for i in keep) + " |\n")
            f.write("|" + "---|" * len(keep) + "\n")
            for row in body:
                cells = [(row[i] if i < len(row) else "") for i in keep]
                cells = [(c[:44] + "…") if len(c) > 45 else c for c in cells]
                # Wrap every cell in backticks. SMILES contain [Br] and (=O),
                # which markdown parses as link syntax — an unwrapped structure
                # column silently renders as broken links instead of chemistry.
                f.write("| " + " | ".join(
                    f"`{c}`" if c else "" for c in
                    (c.replace("|", "\\|").replace("`", "'") for c in cells)
                ) + " |\n")
            f.write("\n")

        prov = ref / ".provenance.md"
        if prov.is_file():
            f.write("---\n\n## Provenance\n\n")
            f.write(prov.read_text(encoding="utf-8").split("\n", 1)[1])


def gen_status_page() -> None:
    """A live status page: what is pinned, what is pending, what is blocked."""
    with mkdocs_gen_files.open("overview/status.md", "w") as f:
        f.write("# Status\n\n" + BANNER)

        lock = REPO / "config" / "sources.lock.json"
        f.write("## Pinned external sources\n\n")
        if lock.is_file():
            pins = json.loads(lock.read_text())["sources"]
            f.write("| Source | Pin |\n|---|---|\n")
            for k, v in sorted(pins.items()):
                val = v.get("sha256") or v.get("commit") or "—"
                f.write(f"| `{k}` | `{val[:16]}…` |\n")
        else:
            f.write("No lockfile yet — run `python -m shared.sources stage`.\n")

        f.write("\n## Open questions\n\n")
        f.write(
            "The choreography's honest limits, kept next to its results rather "
            "than in a file nobody opens.\n\n"
        )
        try:
            from shared import warhead_library as wl
            lib = wl.load(REPO / "data" / "reference" / "warhead_classes_2.csv")
            blocked = lib[lib["structure_status"] != "VERIFIED"]
            f.write("### Warhead classes not enumerable\n\n")
            f.write("| Class | Status | Why |\n|---|---|---|\n")
            for _, r in blocked.iterrows():
                f.write(f"| `{r['class_id']}` | {r['structure_status']} | "
                        f"{str(r['notes'])[:110]} |\n")
            f.write("\n")
        except Exception as exc:  # noqa: BLE001
            f.write(f"*(warhead library unavailable in this environment: {exc})*\n\n")

        f.write(
            "### Not yet acquired\n\n"
            "- **CReM fragment DB** — the radius variant must be chosen and pinned "
            "before T_2 can enumerate; a different radius is a different "
            "neighbourhood definition, not a tuning knob.\n"
            "- **Decoy set** — built by the enrichment gate, which has not run.\n\n"
            "### Not yet verified\n\n"
            "- **Byun 2023 BDHI fragment** — SI-only. The BDHI *class* is verified "
            "(PubChem CID 21983498), which is what the reactivity window needs, so "
            "this blocks enumeration rather than the window.\n"
        )


gen_decisions_index()
gen_reference_page()
gen_status_page()
