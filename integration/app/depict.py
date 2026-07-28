"""
Purpose: 2D structure depiction for the integration GUI.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-28
Input: a SMILES string, optionally a SMARTS to highlight
Output: an inline SVG string

PNG THROUGH `st.image`, NOT INLINE SVG. The first version emitted SVG and passed
it to `st.markdown(..., unsafe_allow_html=True)`. Streamlit sanitises `<svg>` out
of markdown regardless of that flag, so every depiction rendered as nothing —
the page looked like it was showing SMILES only. `st.image` accepts PNG bytes
directly and is the reliable path. `svg()` is kept for anywhere HTML is genuinely
embedded (the grid uses a base64 data URI, which survives sanitisation).

HIGHLIGHTING IS THE INFORMATIVE PART. For the covalent approaches the reactive
warhead is highlighted, so a reader can see at a glance which atoms are the
mechanism and which are the decoration — the distinction that D0025's alert
attribution turns on, and the one a bare depiction hides.

BOTH FORMS ARE SHOWN FOR COVALENT CANDIDATES. `canonical_smiles` is the
pre-reaction molecule the chemist would synthesise; `adduct_smiles` is the
post-reaction species that was actually docked (D0022, D0030). Showing only one
invites exactly the confusion that cost a 1,683-ligand docking run: the chlorine
is present in the first and gone from the second, and that is not a rendering
artefact.
"""

from __future__ import annotations

import logging

import base64

from rdkit import Chem, RDLogger
from rdkit.Chem.Draw import rdMolDraw2D

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger(__name__)

WARHEAD_HIGHLIGHT = (0.95, 0.55, 0.25)     # warm orange — the reactive group
CORE_HIGHLIGHT = (0.55, 0.75, 0.95)        # cool blue — a matched scaffold


def _mol(smiles: str):
    if not smiles or not isinstance(smiles, str):
        return None
    return Chem.MolFromSmiles(smiles)


def svg(smiles: str, *, highlight_smarts: str | None = None,
        width: int = 340, height: int = 260,
        highlight_colour: tuple = WARHEAD_HIGHLIGHT,
        legend: str = "") -> str:
    """Render one molecule, optionally highlighting a substructure."""
    m = _mol(smiles)
    if m is None:
        return (f'<div style="padding:1em;border:1px dashed #999;'
                f'border-radius:6px;color:#999">unparseable SMILES</div>')

    atoms, bonds, colours, bond_colours = [], [], {}, {}
    if highlight_smarts:
        patt = Chem.MolFromSmarts(highlight_smarts)
        if patt is not None:
            match = m.GetSubstructMatch(patt)
            if match:
                atoms = list(match)
                colours = {i: highlight_colour for i in atoms}
                for b in m.GetBonds():
                    if b.GetBeginAtomIdx() in atoms and b.GetEndAtomIdx() in atoms:
                        bonds.append(b.GetIdx())
                        bond_colours[b.GetIdx()] = highlight_colour

    d = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = d.drawOptions()
    opts.addStereoAnnotation = True
    opts.clearBackground = False
    try:
        rdMolDraw2D.PrepareAndDrawMolecule(
            d, m, legend=legend,
            highlightAtoms=atoms or None, highlightBonds=bonds or None,
            highlightAtomColors=colours or None,
            highlightBondColors=bond_colours or None)
    except Exception as exc:  # noqa: BLE001 - a depiction must never break a page
        log.warning("depiction failed for %s: %s", smiles, exc)
        return '<div style="color:#999">could not depict</div>'
    d.FinishDrawing()
    return d.GetDrawingText()


def png(smiles: str, *, highlight_smarts: str | None = None,
        width: int = 340, height: int = 260,
        highlight_colour: tuple = WARHEAD_HIGHLIGHT,
        legend: str = "") -> bytes | None:
    """Render one molecule to PNG bytes, for `st.image`."""
    m = _mol(smiles)
    if m is None:
        return None

    atoms, bonds, colours, bond_colours = [], [], {}, {}
    if highlight_smarts:
        patt = Chem.MolFromSmarts(highlight_smarts)
        if patt is not None:
            match = m.GetSubstructMatch(patt)
            if match:
                atoms = list(match)
                colours = {i: highlight_colour for i in atoms}
                for b in m.GetBonds():
                    if b.GetBeginAtomIdx() in atoms and b.GetEndAtomIdx() in atoms:
                        bonds.append(b.GetIdx())
                        bond_colours[b.GetIdx()] = highlight_colour

    d = rdMolDraw2D.MolDraw2DCairo(width, height)
    d.drawOptions().addStereoAnnotation = True
    try:
        rdMolDraw2D.PrepareAndDrawMolecule(
            d, m, legend=legend,
            highlightAtoms=atoms or None, highlightBonds=bonds or None,
            highlightAtomColors=colours or None,
            highlightBondColors=bond_colours or None)
    except Exception as exc:  # noqa: BLE001 - a depiction must never break a page
        log.warning("PNG depiction failed for %s: %s", smiles, exc)
        return None
    d.FinishDrawing()
    return d.GetDrawingText()


def grid(smiles_list, legends=None, *, per_row: int = 5,
         width: int = 200, height: int = 160) -> str:
    """A responsive flex grid of depictions, as one HTML string."""
    legends = legends or [""] * len(smiles_list)
    cells = []
    for s, lg in zip(smiles_list, legends):
        # base64 data URI: an <img> survives Streamlit's sanitiser where a raw
        # inline <svg> element does not.
        raw = png(s, width=width, height=height)
        img = (f'<img src="data:image/png;base64,'
               f'{base64.b64encode(raw).decode()}" width="{width}"/>'
               if raw else '<div style="color:#999">could not depict</div>')
        cells.append(
            f'<div style="flex:0 0 {width}px;margin:4px;text-align:center;'
            f'font-size:0.75em">{img}'
            f'<div style="color:#888">{lg}</div></div>')
    return ('<div style="display:flex;flex-wrap:wrap;align-items:flex-start">'
            + "".join(cells) + "</div>")


def warhead_smarts(warhead_class: str) -> str | None:
    """The reactive group to highlight for a covalent candidate.

    Read from the warhead library rather than restated here — the same table
    that drives the adduct transform and the decoy class test. A depiction that
    highlighted a different set of atoms than the pipeline reacted would be
    worse than no highlight at all.
    """
    try:
        import sys
        from pathlib import Path
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
        from shared import warhead_library as wl

        lib = wl.load()
        row = lib[lib["class_id"] == warhead_class]
        if row.empty:
            return None
        frag = str(row.iloc[0]["warhead_fragment_smiles"])
        return frag or None
    except Exception:  # noqa: BLE001
        return None
