"""
Purpose: The ONE gnina covalent docking protocol, shared unchanged by T_3 and T_4.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: a candidate SMILES or SDF, its warhead class, the prepared receptor
Output: a docked pose plus CNNaffinity / CNNscore, and a protocol fingerprint

WHY THIS MODULE EXISTS (control S3). T_3 and T_4 are the two covalent
approaches, and the integration phase offers an optional *within-covalent*
re-score comparing their candidates. That comparison is only defensible if both
approaches docked through the identical tool, the identical binary, and the
identical reactive-atom protocol.

"We both ran gnina" is not that claim. Versions differ, builds differ, and the
reactive-atom SMARTS is per-warhead-class and easy to diverge on. So the
protocol is pinned here — binary hash, version string, per-class SMARTS, and
docking parameters — and `protocol_fingerprint()` hashes the lot. T_3 and T_4
each record that fingerprint in their manifests; if the two differ, the GUI
disables the within-covalent re-score rather than silently comparing
incomparable numbers.

DIRECTION WARNING. gnina's `CNNaffinity` is **higher is better** and is NOT in
kcal/mol. It must never be mixed with Vina's affinity, which is kcal/mol and
lower-is-better. See config/choreography.yaml `metrics`.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from .manifest import sha256_file

log = logging.getLogger(__name__)

_REPO_ROOT = Path(__file__).resolve().parent.parent

GNINA_BIN = Path("/data/lab_vm/immutable/inhibition/bin/gnina")
# gnina's "static" build still needs libcudnn.so.9; dwi_gnina exists solely to
# provide it. Depending on another approach's env would break docking whenever
# that env was rebuilt.
GNINA_RUNTIME_ENV = Path("/data/lab_vm/envs/dwi_gnina")
MIN_VERSION = (1, 3)


class CovalentProtocolError(RuntimeError):
    """The covalent protocol is misconfigured, unpinned, or unusable."""


@dataclass(frozen=True)
class DockingParams:
    """Docking parameters, pinned. Changing any of these changes the fingerprint.

    Kept deliberately small: every knob here is one more way T_3 and T_4 can
    diverge without noticing.
    """

    exhaustiveness: int = 16
    num_modes: int = 9
    cnn_scoring: str = "rescore"
    seed: int = 42                    # deterministic; a docking run must replay
    covalent_optimize_lig: bool = True
    covalent_bond_order: int = 1

    def as_flags(self) -> list[str]:
        flags = [
            "--exhaustiveness", str(self.exhaustiveness),
            "--num_modes", str(self.num_modes),
            "--cnn_scoring", self.cnn_scoring,
            "--seed", str(self.seed),
            "--covalent_bond_order", str(self.covalent_bond_order),
        ]
        if self.covalent_optimize_lig:
            flags.append("--covalent_optimize_lig")
        return flags


def gnina_env() -> dict[str, str]:
    """Environment with the CUDA runtime libraries gnina needs on the path."""
    nvidia = GNINA_RUNTIME_ENV / "lib" / "python3.11" / "site-packages" / "nvidia"
    libdirs = [str(p) for p in nvidia.glob("*/lib") if p.is_dir()]
    if not libdirs:
        raise CovalentProtocolError(
            f"no NVIDIA runtime libs under {nvidia}. gnina needs libcudnn.so.9; "
            "build the dwi_gnina env (scripts/setup_envs.sh gnina)."
        )
    env = dict(os.environ)
    env["LD_LIBRARY_PATH"] = ":".join(libdirs + [env.get("LD_LIBRARY_PATH", "")])
    return env


def gnina_version() -> str:
    """Version string reported by the pinned binary."""
    if not GNINA_BIN.is_file():
        raise CovalentProtocolError(
            f"gnina not found at {GNINA_BIN}. Stage it: "
            "python -m shared.sources stage --only gnina"
        )
    proc = subprocess.run([str(GNINA_BIN), "--version"], capture_output=True,
                          text=True, env=gnina_env())
    out = (proc.stdout + proc.stderr).strip().splitlines()
    if not out:
        raise CovalentProtocolError("gnina --version produced no output")
    return out[0].strip()


def assert_version_ok(version_string: str | None = None) -> tuple[int, int, int]:
    """Require gnina >= 1.3, which is where the covalent flags appeared."""
    v = version_string or gnina_version()
    m = re.search(r"v(\d+)\.(\d+)\.(\d+)", v)
    if not m:
        raise CovalentProtocolError(f"cannot parse gnina version from {v!r}")
    parts = tuple(int(x) for x in m.groups())
    if parts[:2] < MIN_VERSION:
        raise CovalentProtocolError(
            f"gnina {parts[0]}.{parts[1]} < required {MIN_VERSION[0]}.{MIN_VERSION[1]}; "
            "the covalent flags do not exist in this build."
        )
    return parts  # type: ignore[return-value]


def assert_covalent_flags() -> list[str]:
    """Confirm the covalent flags actually exist in THIS binary.

    A version check alone is not enough — builds vary. This greps the binary's
    own help, so the guarantee is about the artifact rather than about a number.
    """
    proc = subprocess.run([str(GNINA_BIN), "--help"], capture_output=True,
                          text=True, env=gnina_env())
    text = proc.stdout + proc.stderr
    required = ["--covalent_rec_atom", "--covalent_lig_atom_pattern",
                "--covalent_optimize_lig"]
    missing = [f for f in required if f not in text]
    if missing:
        raise CovalentProtocolError(
            f"pinned gnina lacks covalent flag(s): {missing}. "
            "T_3 and T_4 cannot dock covalently with this build."
        )
    return required


def load_warhead_smarts() -> dict[str, str]:
    """Per-warhead-class reactive-atom SMARTS, from the warhead library.

    The reactive atom differs by MECHANISM — SN2 displacement, Michael addition,
    SN2 ring-opening each mark a different atom. One generic covalent constraint
    applied uniformly would be chemically wrong for most classes, so the SMARTS
    is carried per class in the reference data and read from there.
    """
    import pandas as pd

    lib = _REPO_ROOT / "data" / "reference" / "warhead_classes_2.csv"
    if not lib.is_file():
        raise CovalentProtocolError(f"warhead library not found: {lib}")
    df = pd.read_csv(lib)
    return {r["class_id"]: r["reactive_atom_smarts"] for _, r in df.iterrows()
            if str(r.get("reactive_atom_smarts", "")).strip()}


def receptor_atom() -> str:
    """The covalent attachment atom on the receptor, e.g. ``A:113:SG``."""
    cfg = yaml.safe_load(
        (_REPO_ROOT / "config" / "receptor.yaml").read_text(encoding="utf-8"))
    return cfg["covalent_attachment"]["receptor_atom"]


@dataclass
class CovalentProtocol:
    """The pinned protocol. Its fingerprint is what T_3 and T_4 must match."""

    version: str
    binary_sha256: str
    receptor_atom: str
    warhead_smarts: dict[str, str]
    params: DockingParams = field(default_factory=DockingParams)

    def fingerprint(self) -> str:
        """Stable hash over everything that could make T_3 and T_4 diverge.

        Deliberately includes the binary hash, not just the version string: two
        builds can report the same version and score differently.
        """
        payload = json.dumps({
            "version": self.version,
            "binary_sha256": self.binary_sha256,
            "receptor_atom": self.receptor_atom,
            "warhead_smarts": dict(sorted(self.warhead_smarts.items())),
            "params": self.params.__dict__,
        }, sort_keys=True)
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()

    def to_dict(self) -> dict:
        return {"version": self.version, "binary_sha256": self.binary_sha256,
                "receptor_atom": self.receptor_atom,
                "warhead_smarts": self.warhead_smarts,
                "params": self.params.__dict__,
                "fingerprint": self.fingerprint()}


def load(params: DockingParams | None = None) -> CovalentProtocol:
    """Build the pinned protocol, verifying the binary before returning it."""
    version = gnina_version()
    assert_version_ok(version)
    assert_covalent_flags()
    proto = CovalentProtocol(
        version=version,
        binary_sha256=sha256_file(GNINA_BIN) or "",
        receptor_atom=receptor_atom(),
        warhead_smarts=load_warhead_smarts(),
        params=params or DockingParams(),
    )
    log.info("covalent protocol pinned: %s | fingerprint %s",
             version, proto.fingerprint()[:16])
    return proto


def dock(ligand: Path, out_sdf: Path, warhead_class: str, *,
         receptor_pdbqt: Path | None = None, box: Path | None = None,
         protocol: CovalentProtocol | None = None,
         gpu: bool = True, timeout_s: int = 1800) -> dict:
    """Covalently dock one ligand through the pinned protocol.

    Parameters
    ----------
    ligand : Path
        Ligand file (SDF/MOL2/PDBQT) or a .smi file.
    out_sdf : Path
        Where the docked pose is written.
    warhead_class : str
        Key into the warhead library; selects the mechanism-specific SMARTS.
    gpu : bool
        gnina uses GPU when available; CPU is a fallback, not an equivalent.

    Returns
    -------
    dict
        ``CNNaffinity`` (higher better, dimensionless), ``CNNscore``, the
        protocol fingerprint, and the exact command run.
    """
    proto = protocol or load()
    smarts = proto.warhead_smarts.get(warhead_class)
    if not smarts:
        raise CovalentProtocolError(
            f"no reactive-atom SMARTS for warhead class {warhead_class!r}; "
            f"known: {sorted(proto.warhead_smarts)}"
        )

    rec = receptor_pdbqt or Path(
        "/data/lab_vm/immutable/inhibition/receptor/6VAJ_prepared.pdbqt")
    boxfile = box or Path("/data/lab_vm/immutable/inhibition/receptor/box.json")
    b = json.loads(Path(boxfile).read_text(encoding="utf-8"))

    out_sdf.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(GNINA_BIN),
        "-r", str(rec),
        "-l", str(ligand),
        "-o", str(out_sdf),
        "--center_x", str(b["center_x"]),
        "--center_y", str(b["center_y"]),
        "--center_z", str(b["center_z"]),
        "--size_x", str(b["size_x"]),
        "--size_y", str(b["size_y"]),
        "--size_z", str(b["size_z"]),
        "--covalent_rec_atom", proto.receptor_atom,
        "--covalent_lig_atom_pattern", smarts,
        *proto.params.as_flags(),
    ]
    if not gpu:
        cmd.append("--no_gpu")

    proc = subprocess.run(cmd, capture_output=True, text=True,
                          env=gnina_env(), timeout=timeout_s)
    if proc.returncode != 0:
        raise CovalentProtocolError(
            f"gnina failed ({proc.returncode}) on {ligand.name}:\n"
            f"{(proc.stderr or proc.stdout)[:1500]}")

    rows = parse_results_table(proc.stdout)
    if not rows:
        raise CovalentProtocolError(
            f"gnina exited 0 for {ligand.name} but produced no results table. "
            "A dock that returns no score is a failure, not a null result."
        )
    best = rows[0]

    # gnina emits this on every covalent run, and it matters: the plan makes
    # CNNaffinity T_3's RANK metric and T_4's secondary. Carried into the result
    # so it reaches the manifest and the GUI instead of scrolling past in a log.
    cnn_uncalibrated = "CNN scoring not yet calibrated for covalent docking" in (
        proc.stdout + proc.stderr)
    if cnn_uncalibrated:
        log.warning(
            "gnina reports CNN scoring is NOT calibrated for covalent docking; "
            "treat cnn_affinity as advisory, not as a rank metric (see D0011)")

    return {
        "cnn_affinity": best.get("cnn_affinity"),
        "cnn_score": best.get("cnn_score"),
        "affinity_kcal": best.get("affinity"),
        "cnn_uncalibrated_for_covalent": cnn_uncalibrated,
        "n_modes": len(rows),
        "protocol_fingerprint": proto.fingerprint(),
        "warhead_class": warhead_class,
        "reactive_atom_smarts": smarts,
        "command": " ".join(cmd),
        "pose_path": str(out_sdf),
    }


# gnina's results table has a TWO-LINE header, which a naive single-line parse
# silently misses — returning None for CNNaffinity while the dock itself
# succeeded:
#
#   mode |  affinity  |  intramol  |    CNN     |   CNN
#        | (kcal/mol) | (kcal/mol) | pose score | affinity
#   -----+------------+------------+------------+----------
#       1       -2.35        0.00       0.6286      4.640
#
# The column ORDER is stable across gnina 1.3.x, so rows are parsed positionally
# from the `-----+` rule rather than by matching header text.
_TABLE_RULE = re.compile(r"^-+\+-+")
_COLUMNS = ("mode", "affinity", "intramol", "cnn_score", "cnn_affinity")


def parse_results_table(stdout: str) -> list[dict[str, float]]:
    """Parse gnina's pose table into one dict per mode, best first.

    Returns
    -------
    list of dict
        Keys ``mode``, ``affinity``, ``intramol``, ``cnn_score``,
        ``cnn_affinity``. Empty when no table is present.
    """
    lines = stdout.splitlines()
    start = next((i for i, l in enumerate(lines) if _TABLE_RULE.match(l.strip())), None)
    if start is None:
        return []
    rows: list[dict[str, float]] = []
    for line in lines[start + 1:]:
        parts = line.split()
        if not parts or not parts[0].isdigit():
            if rows:
                break          # table ended
            continue
        if len(parts) < len(_COLUMNS):
            continue
        try:
            rows.append({c: float(v) for c, v in zip(_COLUMNS, parts)})
        except ValueError:
            continue
    return rows


def _parse_best(stdout: str, field_name: str) -> float | None:
    """Top-ranked pose's value for one column, or None if unavailable."""
    rows = parse_results_table(stdout)
    return rows[0].get(field_name) if rows else None
