#!/usr/bin/env python3
"""
Purpose: the presentation GUI's entry page -- left selector, right viewer.
Author: Timothy Wu (with Claude Code)
Date: 2026-09-10
Input: presentation_gui/<ident>.html, written by scripts/presentation_report.py
Output: presentation_gui/index.html

@twu383, 2026-09-10: "same left selector and viewer design [...] but do not
carry over any of the design and numbers from the last guis." Read as: keep
the two-pane PATTERN (a rail on the left switching what an iframe on the
right shows), rebuilt with fresh CSS/JS and grouped the way this page asks
for -- molecule labels, families adjacent -- not the prior pages' layout,
palette, or markup.

The rail-drives-an-iframe mechanism itself is the one piece carried over
deliberately: it is plumbing, not presentation, and it is what
`scripts/mdprio_combine.py` already proved handles many independently-built
viewer pages without the swapped-content script-execution problems a
fetch+innerHTML rail runs into.
"""

from __future__ import annotations

import html
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import outputs as sout                         # noqa: E402
import presentation_build as pb                             # noqa: E402

OUT = sout.Topic("artist", "presentation_gui").dir

CSS = """
:root{
  --ink:#1a1f2b; --sub:#5b6472; --line:#e4e7ec; --bg:#ffffff; --rail:#f7f8fa;
  --accent:#1d6fd6; --accent-bg:#e8f1fd;
}
@media (prefers-color-scheme: dark){
  :root:not([data-force-light]){
    --ink:#e8ebf1; --sub:#9aa4b2; --line:#2a2f3a; --bg:#12151c; --rail:#161a22;
    --accent:#5b9df0; --accent-bg:#1c2a3d;
  }
}
*{box-sizing:border-box}
html,body{height:100%;margin:0}
body{display:flex;background:var(--bg);color:var(--ink);
     font:14px/1.4 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}
#rail{width:280px;flex:0 0 280px;background:var(--rail);border-right:1px solid var(--line);
      overflow-y:auto;padding:14px 0}
#rail .brand{padding:0 16px 14px;font-size:13px;font-weight:700;
             border-bottom:1px solid var(--line);margin-bottom:8px}
#rail .brand small{display:block;font-weight:400;color:var(--sub);margin-top:2px}
.fam-head{padding:12px 16px 4px;font-size:10.5px;text-transform:uppercase;
          letter-spacing:.06em;color:var(--sub);font-weight:700}
.fam-head:first-of-type{padding-top:4px}
button.pose{display:block;width:100%;text-align:left;background:none;border:none;
            color:var(--ink);padding:7px 16px;font-size:13px;cursor:pointer;
            border-left:3px solid transparent}
button.pose:hover{background:var(--accent-bg)}
button.pose.active{background:var(--accent-bg);border-left-color:var(--accent);
                   font-weight:600}
#view{flex:1;min-width:0}
#view iframe{width:100%;height:100%;border:none;display:block}
"""

JS = """
function show(ident){
  document.getElementById('v').src = ident + '.html';
  document.querySelectorAll('button.pose').forEach(function(b){
    b.classList.toggle('active', b.dataset.ident === ident);
  });
}
"""


def main() -> None:
    rail = []
    first = None
    for fam in pb.FAMILY_ORDER:
        entries = [e for e in pb.MANIFEST if e["family"] == fam]
        if not entries:
            continue
        # Default to the family key itself -- which IS the molecule's own id
        # for every group except thiadiazole_benzene, the one override
        # FAMILY_TITLE carries. See its docstring in presentation_build.py.
        rail.append(f'<div class="fam-head">{html.escape(pb.FAMILY_TITLE.get(fam, fam))}</div>')
        for e in entries:
            ident = e["ident"]
            page = OUT / f"{ident}.html"
            if not page.is_file():
                continue                                   # not built yet
            first = first or ident
            rail.append(
                f'<button class="pose" data-ident="{html.escape(ident)}" '
                f'onclick="show(\'{html.escape(ident)}\')">'
                f'{html.escape(e["label"])}</button>')

    if first is None:
        raise SystemExit("no per-pose pages exist yet -- run "
                         "presentation_report.py first")

    page = f"""<!doctype html><html><head><meta charset="utf-8">
<title>Pin1 Cys113 — 100 ns evidence</title><style>{CSS}</style></head><body>
<div id="rail">
  <div class="brand">Pin1 Cys113 covalent candidates
    <small>100 ns MD, {sum(1 for e in pb.MANIFEST if (OUT/(e['ident']+'.html')).is_file())} poses</small>
  </div>
  {''.join(rail)}
</div>
<div id="view"><iframe id="v" src="{html.escape(first)}.html"></iframe></div>
<script>{JS}</script>
</body></html>"""

    OUT.mkdir(parents=True, exist_ok=True)
    dest = OUT / "index.html"
    dest.write_text(page)
    print(f"  index -> {dest}")


if __name__ == "__main__":
    main()
