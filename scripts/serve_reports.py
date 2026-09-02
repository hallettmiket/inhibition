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
import html as _html
import os
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


class LiveRun(NoCache):
    """The CURRENT run's reports, re-resolved on EVERY request.

    THE ROOT USED TO BE RESOLVED ONCE, AT STARTUP. The module docstring says
    the root comes from `run_paths` "so the server cannot end up serving a
    superseded topic's pages -- which it did once, for hours, from a literal
    path." That fixed the literal and left the lifetime: a server started under
    `nac_v5` keeps `nac_v5`'s directory for as long as the process lives, and
    these processes live for weeks.

    Measured 2026-08-31: the server on :8931 had been up 14 days and was serving
    `mdprio_reports_nac_v5/index.html` byte-for-byte, while `run.topic` had been
    `nac_v6` for days and `mdprio_reports_nac_v6/` was being rebuilt every few
    minutes beside it. Nothing was stale-looking: the pages render, the title
    says "DWI covalent screen", and the numbers are a real screen's numbers --
    just the previous one's. That is catalogue #25 exactly, in the component
    whose docstring claims immunity to it.

    `target_config.load` is cached on the file's mtime and size, so this costs a
    dict copy per request and picks up a topic change on the next one.
    """

    def translate_path(self, path):
        self.directory = str(rp.reports_dir())
        return super().translate_path(path)


def run_roots() -> dict:
    """{name: directory} for every browsable run, resolved fresh on each call.

    THREE KINDS OF THING END UP HERE AND THEY ARE NOT EQUIVALENT: the live run,
    a superseded run whose tree is still on disk, and a frozen archive of a
    released one. The chooser labels which is which, because the entire reason
    this project keeps a catalogue is that a superseded page reads as a current
    one (#25, and D0103 where a server served nac_v5 for 14 days).

    Resolved per call, never cached: `run.topic` moves, and a menu built once at
    startup is the D0103 defect with a nicer front end.
    """
    out = {}
    for d in sorted(rp.BLACKSMITH.glob("mdprio_reports_*")):
        if d.is_dir():
            out[d.name[len("mdprio_reports_"):]] = d
    arch = rp.BLACKSMITH / "gui_archive"
    if arch.is_dir():
        for d in sorted(arch.iterdir()):
            if d.is_dir() and d.name not in out:
                out[d.name] = d
    return out


def _chooser_html() -> bytes:
    """The landing page: pick a run. Generated per request, never stored."""
    cur = rp.topic()
    roots = run_roots()
    arch = rp.BLACKSMITH / "gui_archive"
    rows = []
    for name, d in roots.items():
        frozen = arch in d.parents
        live = (name == cur) and not frozen
        pages = sorted(x.name for x in d.glob("*.html"))
        # The landing page of a run tree is index.html where there is one; some
        # trees only ever got the individual pages.
        entry = "index.html" if "index.html" in pages else (pages[0] if pages else "")
        when = ""
        try:
            newest = max(x.stat().st_mtime for x in d.glob("*.html"))
            import datetime as _dt
            when = _dt.datetime.fromtimestamp(newest).strftime("%Y-%m-%d %H:%M")
        except ValueError:
            pass
        if live:
            tag, cls = "CURRENT RUN", "live"
        elif frozen:
            tag, cls = "ARCHIVED — frozen release", "arch"
        else:
            tag, cls = "SUPERSEDED — not the current run", "old"
        links = " · ".join(
            f'<a href="/{_html.escape(name)}/{p}">{p[:-5]}</a>'
            for p in pages if p in ("index.html", "modes.html", "ligands.html",
                                    "combined.html", "sweep.html", "controls.html"))
        rows.append(
            f'<tr class="{cls}"><td><a class="run" href="/{_html.escape(name)}/'
            f'{entry}">{_html.escape(name)}</a></td>'
            f'<td><span class="tag {cls}">{tag}</span></td>'
            f'<td class="n">{len(pages):,} pages</td>'
            f'<td class="n">{when}</td><td class="lk">{links}</td></tr>')
    body = f"""<!doctype html><meta charset="utf-8">
<title>DWI covalent screen — choose a run</title>
<style>
:root{{--fg:#16202b;--mut:#6b7a89;--rule:#dfe5ec;--bg:#f6f8fa;--card:#fff}}
*{{box-sizing:border-box}}
body{{margin:0;padding:2.5rem 1.5rem;background:var(--bg);color:var(--fg);
 font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif}}
.wrap{{max-width:1000px;margin:0 auto}}
h1{{font-size:20px;margin:0 0 .25rem}}
p.sub{{color:var(--mut);margin:0 0 1.75rem}}
table{{width:100%;border-collapse:collapse;background:var(--card);
 border:1px solid var(--rule);border-radius:6px;overflow:hidden}}
td{{padding:.7rem .85rem;border-top:1px solid var(--rule);vertical-align:middle}}
tr:first-child td{{border-top:0}}
a{{color:#1257a8;text-decoration:none}}a:hover{{text-decoration:underline}}
a.run{{font-weight:600;font-size:15px}}
.n{{color:var(--mut);white-space:nowrap;font-variant-numeric:tabular-nums}}
.lk{{font-size:12px;color:var(--mut)}}
.tag{{font-size:11px;font-weight:700;letter-spacing:.03em;padding:.15rem .5rem;
 border-radius:3px;white-space:nowrap}}
.tag.live{{background:#dff3e4;color:#17683a}}
.tag.old{{background:#fdeaea;color:#a12626}}
.tag.arch{{background:#eceff3;color:#556472}}
tr.old a.run,tr.arch a.run{{color:#44515e}}
.note{{margin-top:1.5rem;color:var(--mut);font-size:12.5px;max-width:70ch}}
@media(prefers-color-scheme:dark){{
:root{{--fg:#e6edf3;--mut:#8b98a5;--rule:#2a323c;--bg:#0d1117;--card:#161b22}}
a{{color:#6cb6ff}} .tag.live{{background:#12301f;color:#57d182}}
.tag.old{{background:#3a1b1b;color:#ff8a8a}} .tag.arch{{background:#232a33;color:#9aa7b4}}
tr.old a.run,tr.arch a.run{{color:#aab6c2}}}}
</style>
<div class="wrap">
<h1>DWI covalent screen</h1>
<p class="sub">Choose a run. The current run is <b>{_html.escape(cur)}</b>.</p>
<table>{''.join(rows)}</table>
<p class="note">Every run below is a complete, separately-screened library —
the numbers are <b>not</b> comparable across them without saying so. A
superseded run's pages render exactly like the current one's, which is how a
server left running on the default port served <b>nac_v5</b> for 14 days under
the current run's name (D0103). That is why each row is labelled.</p>
</div>"""
    return body.encode("utf-8")


class MultiRun(NoCache):
    """One server, every run, the first path segment selecting which.

    @tt8804: "can you combine them all and I just select which run I want to
    see". Previously each run needed its own process on its own port, which is
    how four servers came to be running with no way to tell from a page which
    was which.

    The run is resolved PER REQUEST from the directory listing, so a new run
    appears without a restart and the current-run label cannot go stale.
    """

    def _split(self):
        path = self.path.split("?", 1)[0].split("#", 1)[0]
        parts = [p for p in path.split("/") if p]
        return parts

    def do_GET(self):                                          # noqa: N802
        parts = self._split()
        if not parts or (len(parts) == 1 and parts[0] == "index.html"):
            body = _chooser_html()
            self.send_response(200)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
            return
        roots = run_roots()
        if parts[0] in roots and len(parts) == 1:
            # A bare run name must gain its trailing slash, or every RELATIVE
            # link inside the run's pages resolves one level too high.
            self.send_response(301)
            self.send_header("Location", f"/{parts[0]}/")
            self.end_headers()
            return
        super().do_GET()

    def translate_path(self, path):
        parts = [p for p in path.split("?", 1)[0].split("/") if p]
        roots = run_roots()
        if parts and parts[0] in roots:
            self.directory = str(roots[parts[0]])
            rest = "/" + "/".join(parts[1:])
            return super().translate_path(rest)
        # Unknown first segment: serve nothing rather than falling through to
        # some other run's tree, which would put one run's page under another
        # run's URL.
        self.directory = str(rp.BLACKSMITH / "__no_such_run__")
        return super().translate_path(path)


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
    ap.add_argument("--all-runs", action="store_true",
                    help="serve EVERY run behind a chooser, the first path "
                         "segment selecting one (e.g. /nac_v5/modes.html)")
    ap.add_argument("--archive", metavar="NAME", default=None,
                    help="serve a frozen snapshot from gui_archive/ instead of "
                         "the live run, e.g. --archive galena_3.0.0_20260816")
    args = ap.parse_args()
    if args.all_runs:
        root = rp.BLACKSMITH
        handler = partial(MultiRun, directory=str(root))
        names = run_roots()
        print(f"serving {len(names)} runs behind a chooser; current is "
              f"{rp.topic()!r}")
        for n, d in names.items():
            print(f"    /{n}/  -> {d}")
    elif args.archive:
        root = archive_root(args.archive)
        # SAY WHICH RUN THIS IS, EVERY TIME. Two GUIs on two ports that look
        # identical is exactly how a superseded page gets read as current.
        print(f"serving ARCHIVED GUI {args.archive!r} — a frozen snapshot, NOT "
              f"the live run")
        handler = partial(NoCache, directory=str(root))
    else:
        root = rp.reports_dir()
        # RE-RESOLVED PER REQUEST, not pinned here -- see `LiveRun`. The
        # topic this prints is the one at STARTUP; it is a starting point,
        # not a promise about what the server will still be serving next week.
        handler = partial(LiveRun, directory=str(root))
        print(f"live run: topic {rp.topic()!r} (re-resolved on every request)")
    print(f"serving {root} on http://127.0.0.1:{args.port}  (no-store)")
    # THREADING, because `HTTPServer` serves ONE request at a time. A browser
    # holding a connection open, or a slow transfer of a 100 MB report, blocks
    # every other request -- the whole GUI stops answering and looks down.
    # `python -m http.server` has used ThreadingHTTPServer since 3.7; replacing
    # it with the plain class to add no-store headers silently downgraded that.
    ThreadingHTTPServer(("127.0.0.1", args.port), handler).serve_forever()


if __name__ == "__main__":
    main()
