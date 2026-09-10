#!/usr/bin/env python3
"""
Purpose: Screenshot a 3Dmol/py3Dmol HTML page reliably under headless chromium.
Author:  Timothy Wu (with Claude Code)
Date:    2026-09-01
Input:   an HTML file written by pose3d.pose_html or pose_interaction_render
Output:  a PNG of the rendered canvas

Usage:  python scripts/capture_3dmol_page.py <page.html> <out.png> [settle_s] [tries]

Needs playwright + chromium, which are NOT in the shared envs -- build a
throwaway venv rather than installing into /data/lab_vm/envs/*, which other
people use:

    python -m venv /tmp/shotenv
    /tmp/shotenv/bin/pip install playwright pillow
    /tmp/shotenv/bin/playwright install chromium
    /tmp/shotenv/bin/python scripts/capture_3dmol_page.py page.html out.png

WHY IT RETRIES ON A FRESH CONTEXT. Headless swiftshader returns an all-white
frame nondeterministically, and the outcome is STICKY PER PAGE LOAD -- retrying
inside one load failed 8/8, while two near-identical scenes went opposite ways
(12 labels without dashed contacts blank, the same 12 labels WITH them fine).
Adding objects sometimes fixed it, which rules out a complexity ceiling. So the
retry discards the GL context -- new browser, new page -- and keeps the frame
with the most ink.

Two scene properties fail DETERMINISTICALLY and no amount of retrying helps: a
fourth transparent surface, and a surface at opacity 1.0. Keep to three
surfaces at <= 0.85.
"""
from playwright.sync_api import sync_playwright
from PIL import Image
import io, sys, time, pathlib

src, out = sys.argv[1], sys.argv[2]
settle = float(sys.argv[3]) if len(sys.argv) > 3 else 12
tries = int(sys.argv[4]) if len(sys.argv) > 4 else 8
url = "file://" + str(pathlib.Path(src).resolve())


def ink(buf):
    im = Image.open(io.BytesIO(buf)).convert("RGB")
    return sum(1 for r, g, b in im.getdata() if r < 250 or g < 250 or b < 250)


best = None
with sync_playwright() as p:
    for i in range(tries):
        b = p.chromium.launch(args=["--use-gl=swiftshader",
                                    "--enable-unsafe-swiftshader"])
        pg = b.new_page(viewport={"width": 1200, "height": 1000},
                        device_scale_factor=2.5)
        pg.goto(url, wait_until="networkidle", timeout=120000)
        time.sleep(settle)
        pg.evaluate("""()=>{const k=Object.keys(window)
            .filter(n=>n.startsWith('viewer_'));
            if(k.length){try{window[k[0]].render();}catch(e){}}}""")
        time.sleep(2.5)
        buf = pg.screenshot()
        n = ink(buf)
        print(f"  attempt {i+1} (fresh context): {n} ink px")
        b.close()
        if best is None or n > best[0]:
            best = (n, buf)
        if n > 50000:
            break
pathlib.Path(out).write_bytes(best[1])
print("WROTE", out, "with", best[0], "ink px")
