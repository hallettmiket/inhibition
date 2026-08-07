"""
Purpose: turn a finished MD-priority report into a publishable page body.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-07
Input: --candidate <ident>  (reads that molecule's report HTML)
Output: <scratchpad>/reports/<ident>.html — <title> + <style> + body content

The Artifact host supplies its own `<!doctype>`, `<head>` and `<body>`, so a
complete document published as-is would nest one document inside another. This
strips the wrapper and keeps two things: the house `<style>` from
`shared/report_theme.py`, and the body content.

THE SINGLE-THEME LOOK IS DELIBERATE, NOT AN OMISSION. `report_theme` commits to a
white ground with navy/blue accents because that is the format that was asked for.
It is left alone rather than made theme-switching.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

SRC = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith/mdprio_reports")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--candidate", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    src = SRC / f"{args.candidate}.html"
    if not src.is_file():
        raise SystemExit(f"no report yet for {args.candidate} ({src})")
    h = src.read_text()
    style = re.search(r"<style>.*?</style>", h, re.S)
    body = re.search(r"<body>(.*)</body>", h, re.S)
    if not style or not body:
        raise SystemExit(f"{src} is not the expected report shape")

    out = Path(args.outdir) / f"{args.candidate}.html"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(f"<title>{args.candidate} — 100 ns residence</title>\n"
                   f"{style.group(0)}\n{body.group(1)}")
    print(f"{out}  {out.stat().st_size / 1024:.0f} KB")


if __name__ == "__main__":
    main()
