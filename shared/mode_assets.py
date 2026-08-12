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


def poses_dir() -> Path:
    """The production run's representative poses, from `run.topic`.

    NEVER A LITERAL. Drawing poses from one run beside numbers from another is
    the defect that cost 2.2.0 and nearly cost 3.0.0, and a GUI is where it would
    be least visible: the page renders, the poses look like poses, and nothing
    says they came from a different screen than the table beside them.
    """
    from shared import target_config as tc
    return B / f"{tc.get('run.topic')}_poses"


POSES = poses_dir()


def write_assets(out_dir: Path, idents: set[str], force: bool = False,
                 expected: dict | None = None) -> dict:
    """Write `<out>/mode_poses/<ident>.pdb` and `<out>/mode_thumbs/<ident>.svg`.

    Returns counts. Existing files are left alone -- these are derived, and a
    rebuild that rewrites 5,772 files every time is a rebuild nobody runs.

    `expected` maps ident -> how many models the asset must hold, so a STALE
    asset is detected rather than kept. This is not hypothetical tidying: after
    the 3.0.0 topic fix, 195 molecules had no `nac_v4` representative yet and so
    kept the Aug-07 asset holding a SINGLE pose. The ranking listed five modes
    for each of them, the viewer asked for model 3, the model did not exist, and
    the panel rendered an empty box -- @tt8804: "no pose coming up". Skipping on
    mere existence cannot see that; skipping on existence AND sufficiency can.

    Assets it cannot fix are returned under `stale` rather than passed over in
    silence, because the honest statement is "this molecule has no pose from this
    run", and only the caller can say that on the page.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Draw, rdCoordGen
    RDLogger.DisableLog("rdApp.*")

    pd_dir, th_dir = out_dir / "mode_poses", out_dir / "mode_thumbs"
    pd_dir.mkdir(parents=True, exist_ok=True)
    th_dir.mkdir(parents=True, exist_ok=True)
    n_pose = n_thumb = 0
    exp = expected or {}
    have = set()

    def _models(p: Path) -> int:
        try:
            return p.read_text().count("MODEL ")
        except OSError:
            return 0

    for f in sorted(glob.glob(str(POSES / "*.sdf"))):
        ident = Path(f).stem
        if ident not in idents:
            continue
        have.add(ident)
        pdb_out, svg_out = pd_dir / f"{ident}.pdb", th_dir / f"{ident}.svg"
        # SUFFICIENT, not merely present. An asset with fewer models than the
        # ranking has modes is from an earlier run and must be rewritten.
        enough = (_models(pdb_out) >= exp.get(ident, 1)) if pdb_out.exists() else False
        if not force and enough and svg_out.exists():
            continue
        # HEAVY ATOMS ONLY, in both the depiction and the pose.
        #
        # Docked poses carry explicit hydrogens; drawing them gives a 2D
        # structure furred with H labels that no chemist wants to look at, and a
        # 3D pose whose sticks are mostly hydrogen. Everything else in this
        # project displays heavy atoms (the movie PDB has none at all), so these
        # were the odd ones out. `sanitize=False` on read means RemoveHs needs
        # its own sanitize first, or it silently returns the molecule unchanged.
        mols = []
        for m in Chem.SDMolSupplier(f, removeHs=False, sanitize=False):
            if m is None:
                continue
            try:
                Chem.SanitizeMol(m)
                m = Chem.RemoveHs(m)
            except Exception:                              # noqa: BLE001
                m = Chem.RemoveHs(m, sanitize=False)
            mols.append(m)
        if not mols:
            continue

        if force or not enough:
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

        if force or not svg_out.exists():
            try:
                flat = Chem.Mol(mols[0])
                flat.RemoveAllConformers()
                Chem.SanitizeMol(flat)
                flat = Chem.RemoveHs(flat)
                # rdCoordGen, not Compute2DCoords: the default layout puts
                # visibly wrong angles on substituted centres.
                rdCoordGen.AddCoords(flat)
                d = Draw.rdMolDraw2D.MolDraw2DSVG(92, 64)
                d.drawOptions().bondLineWidth = 1
                Draw.rdMolDraw2D.PrepareAndDrawMolecule(d, flat)
                d.FinishDrawing()
                svg_out.write_text(d.GetDrawingText())
                n_thumb += 1
            except Exception as exc:                       # noqa: BLE001
                log.debug("%s: no depiction (%s)", ident, exc)

    # NAMED, NOT SWALLOWED. A molecule the ranking lists but this run produced no
    # pose for cannot be fixed here -- there is nothing to write. Returning the
    # list lets the page say "no pose from this run" instead of showing an empty
    # viewer, which is indistinguishable from a broken one.
    stale = sorted(i for i in idents if i not in have)
    if stale:
        log.warning("%d molecules in the ranking have no pose under %s "
                    "(their viewer will report it, not blank)",
                    len(stale), POSES.name)
    return {"poses": n_pose, "thumbs": n_thumb, "stale": stale}
