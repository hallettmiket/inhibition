#!/usr/bin/env python3
"""
Purpose: serve the current run's report pages, never from cache.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-17
Input: --port (default 8931)
Output: an HTTP server on loopback, rooted at this run's reports directory

WHY NOT `python -m http.server`. It sends `Last-Modified` and no
`Cache-Control`, so a browser is free to reuse a page it already has -- and
these pages are rebuilt every minute under a STABLE filename. The result is
that a rebuild lands on disk, the reader reloads, and sees the previous
version, with nothing anywhere saying so. That produced several rounds of "I
rebuilt it" / "it still shows the old thing", and the wrong conclusion each
time was that the build was broken.

`no-store` on every response, and the root is resolved from run_paths rather
than typed, so the server cannot end up serving a superseded topic's pages --
which it did once, for hours, from a literal path.
"""

from __future__ import annotations

import argparse
import sys
from functools import partial
from http.server import HTTPServer, SimpleHTTPRequestHandler
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import run_paths as rp                    # noqa: E402


class NoCache(SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header("Cache-Control", "no-store, must-revalidate")
        self.send_header("Pragma", "no-cache")
        self.send_header("Expires", "0")
        super().end_headers()

    def log_message(self, *a):                                 # noqa: A003
        pass


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--port", type=int, default=8931)
    args = ap.parse_args()
    root = rp.reports_dir()
    print(f"serving {root} on http://127.0.0.1:{args.port}  (no-store)")
    HTTPServer(("127.0.0.1", args.port),
               partial(NoCache, directory=str(root))).serve_forever()


if __name__ == "__main__":
    main()
