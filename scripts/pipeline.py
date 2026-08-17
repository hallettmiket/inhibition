#!/usr/bin/env python3
"""
Purpose: drive the screen — one command to see, start or stop any stage.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-17
Input: status | start <stage> | stop <stage> | auto | serve
Output: stage processes, and pipeline_state.json for the dashboard

@tt8804: "why are we still relying on you to do this cannot we code this as a
pipeline, we can use the gui as a dashboard to start each stage".

WHAT THIS REPLACES. Running 3.1.0 by hand needed roughly a dozen bespoke shell
scripts, and the stage code was never what broke -- the glue was, differently
each time. `shared/pipeline.py` records those failures as rules; this is the
interface to them.

  status   every stage, its progress counted from artefacts, and whether it can
           start. `unknown` is printed as `unknown`, never as zero.
  start    a stage, detached, refusing if its inputs are not done or if it is
           already running.
  stop     a stage and its children, PARENTS FIRST -- `md_residence_3ikd`
           relaunches `gmx`, so killing `gmx` first just makes a new one.
  auto     start every stage that is ready, then keep going as each completes.
  serve    a control endpoint the dashboard posts to.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import pipeline as pl                     # noqa: E402

log = logging.getLogger("pipeline")

_MARK = {"done": "[x]", "running": "[>]", "waiting": "[ ]",
         "stopped": "[!]", "unknown": "[?]"}


def _bar(done, total, w=22) -> str:
    if done is None or not total:
        return " " * w
    k = int(w * min(done, total) / total)
    return "#" * k + "." * (w - k)


def cmd_status(_args) -> int:
    st = pl.status()
    print(f"\n  3.1.0  topic {st['topic']}   "
          f"sweep {st['spec']['sweep_ns']:.0f} ns / production "
          f"{st['spec']['production_ns']:.0f} ns / bar "
          f"{st['spec']['survivor_rmsd_nm']} nm")
    print(f"  scope: {', '.join(st['spec']['families'])}  "
          f"tiers {','.join(st['spec']['tiers'])}\n")
    for s in st["stages"]:
        n = "?" if s["done"] is None else s["done"]
        t = "?" if s["total"] is None else s["total"]
        extra = ""
        if s["state"] == "running":
            extra = f"  {len(s['pids'])} proc"
        elif s["state"] == "waiting" and s["ready"]:
            extra = "  READY"
        if s["error"]:
            extra = f"  !! {s['error']}"
        print(f"  {_MARK[s['state']]} {s['title']:<24} {_bar(s['done'], s['total'])} "
              f"{n:>5}/{t:<5}{extra}")
    print()
    pl.write_status()
    return 0


def cmd_start(args) -> int:
    try:
        pid = pl.start(args.stage)
    except pl.StageError as exc:
        print(f"  refused: {exc}")
        return 1
    print(f"  started {args.stage} (pid {pid}); log "
          f"{pl.run_dir() / (args.stage + '.log')}")
    pl.write_status()
    return 0


def cmd_stop(args) -> int:
    n = pl.stop(args.stage)
    print(f"  stopped {args.stage}: {n} process(es) killed")
    pl.write_status()
    return 0


def cmd_auto(args) -> int:
    """Start whatever is ready, and keep starting as stages complete.

    Every stage is idempotent and resumable, so this is safe to run against a
    pipeline that is already partly done or already running -- it adopts what
    is there rather than restarting it.
    """
    while True:
        st = pl.status()
        for s in st["stages"]:
            if s["state"] in ("waiting", "stopped") and s["ready"] and s["total"]:
                if s["done"] and s["done"] >= s["total"]:
                    continue
                try:
                    pid = pl.start(s["name"])
                    log.info("started %s (pid %d)", s["name"], pid)
                except pl.StageError as exc:
                    log.info("not starting %s: %s", s["name"], exc)
        pl.write_status()
        if all(s["state"] == "done" for s in st["stages"]):
            log.info("every stage done")
            return 0
        if not args.watch:
            return 0
        time.sleep(args.poll)


class _Handler(BaseHTTPRequestHandler):
    """Tiny control endpoint. Loopback only, no auth -- same trust model as the
    report server it sits beside."""

    def _send(self, code, body, ctype="application/json"):
        b = body.encode()
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(b)))
        self.send_header("Access-Control-Allow-Origin", "*")
        self.end_headers()
        self.wfile.write(b)

    def do_GET(self):                                          # noqa: N802
        if self.path.rstrip("/") in ("", "/status"):
            self._send(200, json.dumps(pl.status(), default=str))
        else:
            self._send(404, json.dumps({"error": "not found"}))

    def do_POST(self):                                         # noqa: N802
        parts = [p for p in self.path.split("/") if p]
        if len(parts) != 2 or parts[0] not in ("start", "stop"):
            self._send(404, json.dumps({"error": "use /start/<stage> or /stop/<stage>"}))
            return
        action, stage = parts
        if stage not in pl.BY_NAME:
            self._send(404, json.dumps({"error": f"unknown stage {stage}"}))
            return
        try:
            if action == "start":
                pid = pl.start(stage)
                self._send(200, json.dumps({"ok": True, "pid": pid}))
            else:
                n = pl.stop(stage)
                self._send(200, json.dumps({"ok": True, "killed": n}))
        except pl.StageError as exc:
            # A refusal is a 409, not a 500: "inputs are not ready" is a normal
            # answer to a button press, and the dashboard should show the reason
            # rather than a stack trace.
            self._send(409, json.dumps({"ok": False, "error": str(exc)}))
        pl.write_status()

    def log_message(self, *a):                                 # noqa: A003
        pass


def cmd_serve(args) -> int:
    pl.write_status()
    # Threaded for the same reason the report server is: `status()` walks the
    # trajectory tree and can take seconds, and a serial server would make the
    # dashboard's poll block its own start/stop buttons.
    srv = ThreadingHTTPServer(("127.0.0.1", args.port), _Handler)
    log.info("pipeline control on http://127.0.0.1:%d  (GET /status, "
             "POST /start/<stage>, POST /stop/<stage>)", args.port)
    srv.serve_forever()
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("status")
    for name in ("start", "stop"):
        p = sub.add_parser(name)
        p.add_argument("stage", choices=[s.name for s in pl.STAGES])
    a = sub.add_parser("auto")
    a.add_argument("--watch", action="store_true")
    a.add_argument("--poll", type=int, default=300)
    s = sub.add_parser("serve")
    s.add_argument("--port", type=int, default=8932)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        datefmt="%H:%M:%S")
    return {"status": cmd_status, "start": cmd_start, "stop": cmd_stop,
            "auto": cmd_auto, "serve": cmd_serve}[args.cmd](args)


if __name__ == "__main__":
    raise SystemExit(main())
