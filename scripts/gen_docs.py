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

REPO_ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_))


def latest_reference(stem: str) -> str:
    """Highest integer-versioned `<stem>_N.csv` under data/reference.

    Resolved by GLOB, not by importing the loader that reads it. The docs are
    built in CI with only mkdocs + pandas installed -- importing
    `shared.reference_set` would pull in RDKit and break `mkdocs build
    --strict`. Globbing gets the same answer with no dependency, and it still
    cannot go stale: a new version is picked up the moment it is written.

    This page previously named `warhead_classes_2.csv` and
    `pin1_reference_binders_1.csv` as hard-coded strings, which is why the
    rendered status page carried a load error for weeks.
    """
    ref = REPO_ / "data" / "reference"
    hits = []
    for p in ref.glob(f"{stem}_*.csv"):
        tail = p.stem[len(stem) + 1:]
        if tail.isdigit():
            hits.append((int(tail), p.name))
    return max(hits)[1] if hits else f"{stem}_1.csv"

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
            (latest_reference("pin1_reference_binders"), "Master set — novelty axis",
             "Rows marked `UNVERIFIED` are excluded from the novelty computation."),
            (latest_reference("pin1_covalent_cys113_anchors"), "Covalent Cys113 anchors — reactivity window",
             "`reference_set.py` refuses `UNVERIFIED` rows into the window."),
            (latest_reference("warhead_classes"), "Warhead classes — T_4 enumeration",
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
            lib = wl.load()   # whatever warhead_library actually defaults to
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
            "### Where the gates stand\n\n"
            "The gate verdicts below are read from the live token, so this "
            "section cannot drift from what the pipeline actually recorded.\n\n"
        )
        try:
            import json as _json
            tok = _json.loads((Path("/data/lab_vm/append_only/inhibition/"
                                    "00_shared_substrate/enrichment_gate.token")
                               ).read_text())
            f.write("| stratum / metric | verdict | ROC-AUC | 95% CI | actives | "
                    "chemotypes | EF 1% |\n|---|---|---|---|---|---|---|\n")
            for stratum, sd in (tok.get("strata") or {}).items():
                for metric, m in (sd.get("metrics") or {}).items():
                    ci = m.get("roc_auc_ci") or [None, None]
                    ci_s = (f"[{ci[0]:.3f}, {ci[1]:.3f}]"
                            if ci and ci[0] is not None else "—")
                    f.write(f"| `{stratum}/{metric}` | **{m.get('verdict')}** | "
                            f"{m.get('roc_auc', float('nan')):.3f} | {ci_s} | "
                            f"{m.get('n_actives','?')} | "
                            f"{m.get('n_chemotypes','?')} | "
                            f"{m.get('ef_1pct', 0):.1f} |\n")
            f.write("\nThe floor for any verdict above UNDERPOWERED is 3 actives "
                    "AND 6 independent chemotypes.\n\n")
        except Exception as exc:  # noqa: BLE001
            f.write(f"*(gate token unreadable here: {exc})*\n\n")

        f.write(
            "### Known limits of the rankings\n\n"
            "- **The scores rank partly on molecular size** (D0043). Spearman "
            "between heavy-atom count and the ranking metric is +0.745 for T_3, "
            "−0.617 for T_1, +0.305 for T_4 and −0.230 for T_2 — every sign "
            "meaning larger molecules score better. Ligand efficiency "
            "over-corrects (−0.938 for T_1) and is not the fix.\n"
            "- **Inhibition versus activation is not resolved computationally** "
            "by any approach here. Catalytic-site occupancy is the proxy.\n"
            "- **No wet-lab ground truth** exists for any candidate.\n\n"
            "### Not yet verified\n\n"
            "- **Byun 2023 BDHI fragment** — SI-only. The BDHI *class* is verified "
            "(PubChem CID 21983498), which is what the reactivity window needs, so "
            "this blocks enumeration rather than the window.\n"
        )


gen_decisions_index()
gen_reference_page()
gen_status_page()
