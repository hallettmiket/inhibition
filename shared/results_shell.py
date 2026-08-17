"""The results-page shell: the layout the MD results page is built from.

@tt8804: "just copy over the design for MD results and add/modify the tables to
use the sweep result data".

SO IT IS MOVED HERE, NOT COPIED. `scripts/mdprio_combine` carried this CSS inline
and the sweep page needed the same layout -- a left rail of rows, a right viewer
holding one report at a time. Transcribing it would have produced two stylesheets
that start identical and drift, which is the failure this project keeps paying
for. Both pages now interpolate this constant.

The rail is 376px and the viewer takes the rest. That ratio is the design: the
thing being selected FROM is narrow, the thing being looked AT is wide. A sweep
page built with those reversed was the reason this module exists.
"""

from __future__ import annotations

#: The rail ROW, on its own, because three pages draw the same row and only two
#: of them could share a stylesheet. The ranking page keeps its own <style> (it
#: is virtualised and needs a fixed row height), so it used to hold a COPY of
#: these rules -- and the copy had already drifted: `.l2` lost `flex-wrap`, and
#: whitespace differences aside, nothing kept them equal. @tt8804: "the different
#: selectors look different on diff pages".
#:
#: Both now interpolate this. The ranking page adds only what virtualisation
#: requires, so a divergence has to be written deliberately rather than by
#: forgetting to copy an edit.
SEARCH_CSS = """\
/* Rail search. Sticky, because a filter you have to scroll back up to change is
   one you stop using on a 147-row list. */
.railq{position:sticky;top:0;z-index:3;display:flex;gap:6px;align-items:center;
 padding:6px 8px;background:var(--raise);border-bottom:1px solid var(--rule)}
.railq input{flex:1;min-width:0;font:12px var(--mono);padding:.3rem .45rem;
 border:1px solid var(--rule);border-radius:3px;background:var(--paper);
 color:var(--ink)}
.railq input:focus{outline:none;border-color:var(--blue)}
.railq .n{font:600 10px var(--sans);color:var(--muted);white-space:nowrap}
.railq button{font:600 11px var(--sans);border:1px solid var(--rule);
 border-radius:3px;background:var(--paper);color:var(--muted);cursor:pointer;
 padding:.25rem .4rem}
.railq button:hover{background:var(--blue-pale);color:var(--ink)}
"""

ROW_CSS = """\
.row{display:grid;grid-template-columns:22px 46px 1fr;gap:8px;align-items:start;
 width:100%;text-align:left;font:inherit;color:inherit;background:none;cursor:pointer;
 padding:9px 14px 8px;border:0;border-bottom:1px solid var(--rule)}
.row:hover{background:var(--blue-pale)}
.row.on{background:var(--blue-pale);box-shadow:inset 3px 0 0 var(--blue)}
.row:focus-visible{outline:2px solid var(--blue);outline-offset:-2px}
.rk{font-family:var(--mono);font-size:11px;font-weight:600;color:var(--muted);
 padding-top:3px}
.thumb{width:46px;height:32px;object-fit:contain;background:#fff;
 border:1px solid var(--rule);border-radius:3px;display:block}
.body{min-width:0;display:flex;flex-direction:column;gap:3px}
.l1{display:flex;align-items:baseline;justify-content:space-between;gap:8px}
.mid-id{font-family:var(--mono);font-size:12.5px;overflow:hidden;
 text-overflow:ellipsis;white-space:nowrap}
/* Which pose the run started from. `unk` is deliberately muted: an unrecorded
   mode is an absence, and styling it like a value would let the eye read it as
   one. */
.mode{font-family:var(--mono);font-size:10px;margin-left:5px;padding:0 4px;
 border-radius:3px;background:rgba(0,114,206,.12);color:var(--blue)}
.mode.unk{background:transparent;color:var(--muted,#8a94a6);font-style:italic}
.eng{font-family:var(--mono);font-size:11.5px;font-weight:600;color:var(--navy);
 flex:none}
.l2{display:flex;align-items:center;gap:6px;flex-wrap:nowrap;min-width:0}
.wc{font-size:10.5px;color:var(--blue);white-space:nowrap;overflow:hidden;
 text-overflow:ellipsis;max-width:11ch}
.meta{font-size:10.5px;color:var(--muted);white-space:nowrap;overflow:hidden;
 text-overflow:ellipsis;flex:1}
.tag{font-size:9px;letter-spacing:.06em;text-transform:uppercase;font-weight:700;
 padding:1px 6px;border-radius:99px;flex:none}
.t-held{background:#e6f4ee;color:var(--good)}
.t-left{background:#fbeae8;color:var(--bad)}
/* Controls read as a different KIND of row, not a better or worse one: an amber
   accent and a monospace badge instead of a structure thumbnail, because the
   thing to notice is where they land in the ranking, not what they look like. */
.t-ctl{background:#fdf0dc;color:#8a5a00}
.row.ctl{background:#fffaf2;box-shadow:inset 3px 0 0 #d99a2b}
.row.ctl:hover{background:#fdf3e4}
.thumb.tctl{display:flex;align-items:center;justify-content:center;
 font-family:var(--mono);font-size:10px;font-weight:700;color:#8a5a00;
 background:#fdf0dc;border-color:#e8cfa5}
:root[data-theme="dark"] .t-ctl{background:#3a2c14;color:#e0b070}
:root[data-theme="dark"] .row.ctl{background:#1b1710}
:root[data-theme="dark"] .row.ctl:hover{background:#241d13}
.bar{height:3px;background:var(--rule);border-radius:2px;overflow:hidden;margin-top:2px}
.bar i{display:block;height:100%;background:var(--blue)}
"""

CSS = """\
/* The per-molecule reports are light-only and use this exact palette
   (shared/report_theme.py). The shell inherits it so the frame and its
   contents read as one document rather than two. */
:root{--ink:#10233f;--navy:#003087;--blue:#0072ce;--blue-pale:#e8f1fb;
 --rule:#d6dee8;--muted:#5b6b80;--paper:#fff;--raise:#f5f8fc;
 --good:#0f7a54;--bad:#b3261e;
 --sans:"Helvetica Neue",Helvetica,Arial,system-ui,sans-serif;
 --mono:ui-monospace,"SF Mono",Menlo,Consolas,monospace}
*{box-sizing:border-box}
html,body{height:100%}
body{margin:0;background:var(--paper);color:var(--ink);font-family:var(--sans);
 font-size:14px;line-height:1.5;font-variant-numeric:tabular-nums;
 display:flex;flex-direction:column}
/* ONE STRIP. The title, both toggles and the hint share a single ~38px bar --
   roughly one selector row -- because everything below it is the actual work and
   a two-line masthead over a scrolling list is space the reader never gets back. */
#topbar{display:flex;align-items:center;gap:8px;padding:6px 14px;min-height:38px;
 border-bottom:1px solid var(--rule);background:var(--raise);
 overflow-x:auto;white-space:nowrap;scrollbar-width:thin}
h1{margin:0;font-size:.86rem;font-weight:600;letter-spacing:-.01em;color:var(--navy);
 flex:none}
.mbtn{font:11.5px var(--sans);padding:3px 10px;flex:none;border:1px solid var(--rule);
 background:var(--paper);color:var(--muted);border-radius:99px;cursor:pointer}
.mbtn.on{background:var(--navy);border-color:var(--navy);color:#fff;font-weight:600}
.mbtn:focus-visible{outline:2px solid var(--blue);outline-offset:2px}
.mhint{font-size:11px;color:var(--muted);flex:none;margin-left:4px}
.ver{font-family:var(--mono);font-size:10.5px;color:var(--muted);flex:none;
 border:1px solid var(--rule);border-radius:99px;padding:1px 8px;
 background:var(--paper)}
.msep{width:1px;height:16px;background:var(--rule);margin:0 2px;flex:none}
.ohd{font-family:var(--mono);font-size:.6rem;letter-spacing:.14em;text-transform:uppercase;
 font-weight:700;padding:11px 14px 7px;border-bottom:1px solid var(--rule);
 position:sticky;top:0;z-index:2}
.o-held{background:#e6f4ee;color:var(--good)}
.o-left{background:#fbeae8;color:var(--bad)}
.o-ctl{background:#fdf0dc;color:#8a5a00}
.o-pend{background:var(--raise);color:var(--muted)}
:root[data-theme="dark"] .o-ctl{background:#3a2c14;color:#e0b070}
.mbtn:disabled{opacity:.4;cursor:default}
a.mbtn.lnk{text-decoration:none;color:var(--blue);border-color:var(--blue);
 flex:none;line-height:1.5}
a.mbtn.lnk:hover{background:var(--blue-pale)}
.legend{font-size:11px;color:var(--muted);padding:8px 14px;background:var(--raise);
 border-bottom:1px solid var(--rule);line-height:1.45}
.legend b{color:var(--navy)}
/* A molecule that has only been triaged shows its sweep reading in muted type, so
   the ranked number and the not-yet-ranked number cannot be read as one column. */
.eng.pend{color:var(--muted);font-weight:500}
.t-pend{background:var(--raise);color:var(--muted)}
main{flex:1;display:grid;grid-template-columns:376px 1fr;min-height:0}
@media(max-width:880px){main{grid-template-columns:1fr;grid-template-rows:250px 1fr}}
#rail{overflow-y:auto;border-right:1px solid var(--rule);background:var(--rail)}
""" + SEARCH_CSS + """
.chd{font-family:var(--mono);font-size:.6rem;letter-spacing:.14em;text-transform:uppercase;
 color:var(--blue);font-weight:600;padding:10px 14px 6px;background:var(--raise);
 border-bottom:1px solid var(--rule);position:sticky;top:0;z-index:1}
""" + ROW_CSS + """
:root[data-theme="dark"] .thumb.tctl{background:#3a2c14;color:#e0b070;border-color:#4a3a1e}
#viewer{min-width:0;min-height:0;display:flex;flex-direction:column;background:var(--paper)}
#vhead{padding:9px 18px;border-bottom:1px solid var(--rule);background:var(--raise);
 display:flex;justify-content:space-between;align-items:center;gap:12px}
#vname{font-family:var(--mono);font-size:12.5px;font-weight:600;color:var(--navy)}
#vhead a{font-size:12px;color:var(--blue);text-decoration:none}
#vhead a:hover{text-decoration:underline}
iframe{flex:1;width:100%;border:0;background:var(--paper)}
.tbtn{margin-left:6px}
:root[data-theme="dark"]{--ink:#dfe7f0;--navy:#8ab4e8;--blue:#6aa9e0;
 --blue-pale:#16283a;--rule:#25333f;--muted:#93a3b4;--paper:#0e151c;
 --raise:#16202a;--rail:#121b24;--good:#4fc4a0;--bad:#e08a70}
:root[data-theme="dark"] .thumb{background:#fff}
{_stepcss}
"""


#: The rail's search box and its filter. Shared so the two rail pages cannot
#: drift into two behaviours -- which is the same reason `CSS` lives here.
#:
#: FILTERS THE DOM, NOT THE DATA. Both rails are static markup, so hiding rows
#: is enough and nothing has to be re-rendered or re-sorted. The rank numbers
#: keep their original values on purpose: renumbering 1..n over a filtered list
#: would make "rank 3" mean something different depending on what was typed.
#:
#: Matches the row's whole text -- ident, mode, warhead class, the readings --
#: so `bdhi_c5`, `held`, `0.3`, and a molecule id all work without the user
#: having to know which field they are searching.


def search_html(oninput: str = "railFilter()",
                placeholder: str = "filter — id, class, held/left, a number") -> str:
    """The rail's search box. Takes its handler because the ranking page cannot
    use the DOM filter -- its rail is virtualised, so it filters the DATA and
    rebuilds -- but it must still be the SAME box in the SAME place.

    IT LIVES IN THE RAIL, NOT THE TOPBAR. The first version put the ranking
    page's input in its topbar, which is a flex row with `overflow-x:auto`
    carrying a long title, the class select, a hint and two buttons -- so the
    box was served, well-formed, and scrolled out of the visible strip.
    @tt8804, three times: "still no search bar". A control the reader cannot
    find is a control that does not exist.
    """
    return (f'<div class="railq">\n'
            f' <input id="railq" type="search" placeholder="{placeholder}"\n'
            f'        oninput="{oninput}" autocomplete="off" spellcheck="false">\n'
            f' <span class="n" id="railn"></span>\n'
            f' <button type="button" onclick="railClear()" title="clear">&times;</button>\n'
            f'</div>')


#: Back-compat for the two pages that use the default handler.
SEARCH_HTML = search_html()

SEARCH_JS = """
function railFilter(){
  const el = document.getElementById('railq');
  const q = (el ? el.value : '').trim().toLowerCase();
  const rows = document.querySelectorAll('#rail .row');
  let shown = 0;
  rows.forEach(function(r){
    const hit = !q || r.textContent.toLowerCase().indexOf(q) !== -1;
    r.style.display = hit ? '' : 'none';
    if (hit) shown++;
  });
  // Section headers with nothing under them would otherwise float free.
  document.querySelectorAll('#rail .ohd').forEach(function(h){
    let n = 0;
    for (let s = h.nextElementSibling; s && !s.classList.contains('ohd');
         s = s.nextElementSibling) {
      if (s.classList.contains('row') && s.style.display !== 'none') n++;
    }
    h.style.display = n ? '' : 'none';
  });
  const c = document.getElementById('railn');
  if (c) c.textContent = q ? shown + ' / ' + rows.length : rows.length + '';
}
function railClear(){
  const el = document.getElementById('railq');
  if (el) { el.value = ''; railFilter(); el.focus(); }
}
// `/` focuses the box, Escape clears it -- the list is long enough that reaching
// for the mouse to filter is the thing that stops people filtering.
document.addEventListener('keydown', function(e){
  const el = document.getElementById('railq');
  if (!el) return;
  if (e.key === '/' && document.activeElement !== el) { e.preventDefault(); el.focus(); }
  else if (e.key === 'Escape' && document.activeElement === el) railClear();
});
"""
