"""Per-molecule assets the ranking view fetches on demand.

The ranking view lists 8,096 modes. Inlining a depiction and a pose for each
would be tens of megabytes of base64 in one document, so they are written as
files beside the page and fetched when a row is shown or selected. The results
GUI can inline its 59 rows; this one cannot, and that is the only reason the two
differ.

One file per MOLECULE, not per mode: a molecule's modes share a depiction, and
its poses are models inside one PDB, so selecting a different mode of the same
molecule is a model switch rather than another fetch.
"""

from __future__ import annotations

import glob
import logging
from pathlib import Path

log = logging.getLogger("mode-assets")

B = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith")
POSES = B / "nac_v3_poses"


def write_assets(out_dir: Path, idents: set[str]) -> dict:
    """Write `<out>/mode_poses/<ident>.pdb` and `<out>/mode_thumbs/<ident>.svg`.

    Returns counts. Existing files are left alone -- these are derived, and a
    rebuild that rewrites 5,772 files every time is a rebuild nobody runs.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Draw, rdCoordGen
    RDLogger.DisableLog("rdApp.*")

    pd_dir, th_dir = out_dir / "mode_poses", out_dir / "mode_thumbs"
    pd_dir.mkdir(parents=True, exist_ok=True)
    th_dir.mkdir(parents=True, exist_ok=True)
    n_pose = n_thumb = 0

    for f in sorted(glob.glob(str(POSES / "*.sdf"))):
        ident = Path(f).stem
        if ident not in idents:
            continue
        pdb_out, svg_out = pd_dir / f"{ident}.pdb", th_dir / f"{ident}.svg"
        if pdb_out.exists() and svg_out.exists():
            continue
        mols = [m for m in Chem.SDMolSupplier(f, removeHs=False, sanitize=False)
                if m is not None]
        if not mols:
            continue

        if not pdb_out.exists():
            # MODEL n == the pose's own `mode`, so the viewer can select a model
            # by mode number instead of by position in the file. Position is what
            # #53 was about.
            parts = []
            for m in mols:
                mode = int(m.GetProp("mode")) if m.HasProp("mode") else -1
                body = "\n".join(l for l in Chem.MolToPDBBlock(m).splitlines()
                                 if l.startswith(("ATOM", "HETATM")))
                body = body.replace("UNL", "MOL")
                parts.append(f"MODEL     {mode}\n{body}\nENDMDL")
            pdb_out.write_text("\n".join(parts) + "\n")
            n_pose += 1

        if not svg_out.exists():
            try:
                flat = Chem.Mol(mols[0])
                flat.RemoveAllConformers()
                Chem.SanitizeMol(flat)
                rdCoordGen.AddCoords(flat)
                d = Draw.rdMolDraw2D.MolDraw2DSVG(92, 64)
                d.drawOptions().bondLineWidth = 1
                Draw.rdMolDraw2D.PrepareAndDrawMolecule(d, flat)
                d.FinishDrawing()
                svg_out.write_text(d.GetDrawingText())
                n_thumb += 1
            except Exception as exc:                       # noqa: BLE001
                log.debug("%s: no depiction (%s)", ident, exc)

    return {"poses": n_pose, "thumbs": n_thumb}
