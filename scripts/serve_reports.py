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
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
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


def archive_root(name: str) -> Path:
    """A snapshot under `gui_archive/`, resolved by NAME rather than by path.

    THE ROOT IS STILL NOT TYPED, AND THAT IS THE POINT. This server exists
    because it once served a superseded topic's pages for hours from a literal
    path, so a general `--root` would hand that hazard straight back. An
    archived GUI is a real need -- a released run has to stay browsable after
    the topic moves on -- but an archive is a NAMED, FROZEN thing, so it is
    resolved inside `gui_archive/` and nowhere else.

    A live topic directory therefore cannot be reached through this flag, even
    by accident: `--archive nac_v4` does not resolve, because `nac_v4` is a run,
    not a snapshot.
    """
    base = rp.BLACKSMITH / "gui_archive"
    root = (base / name).resolve()
    if base.resolve() not in root.parents:
        raise SystemExit(
            f"{name!r} does not resolve inside {base}. This flag serves frozen "
            f"snapshots only; the live run is served with no argument.")
    if not root.is_dir():
        have = sorted(p.name for p in base.iterdir() if p.is_dir()) \
            if base.is_dir() else []
        raise SystemExit(
            f"no archived GUI called {name!r}. Available: {have or 'none'}")
    return root


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--port", type=int, default=8931)
    ap.add_argument("--archive", metavar="NAME", default=None,
                    help="serve a frozen snapshot from gui_archive/ instead of "
                         "the live run, e.g. --archive galena_3.0.0_20260816")
    args = ap.parse_args()
    if args.archive:
        root = archive_root(args.archive)
        # SAY WHICH RUN THIS IS, EVERY TIME. Two GUIs on two ports that look
        # identical is exactly how a superseded page gets read as current.
        print(f"serving ARCHIVED GUI {args.archive!r} — a frozen snapshot, NOT "
              f"the live run")
    else:
        root = rp.reports_dir()
    print(f"serving {root} on http://127.0.0.1:{args.port}  (no-store)")
    # THREADING, because `HTTPServer` serves ONE request at a time. A browser
    # holding a connection open, or a slow transfer of a 100 MB report, blocks
    # every other request -- the whole GUI stops answering and looks down.
    # `python -m http.server` has used ThreadingHTTPServer since 3.7; replacing
    # it with the plain class to add no-store headers silently downgraded that.
    ThreadingHTTPServer(("127.0.0.1", args.port),
                        partial(NoCache, directory=str(root))).serve_forever()


if __name__ == "__main__":
    main()
