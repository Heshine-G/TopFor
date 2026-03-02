# modules/run_antechamber.py
from __future__ import annotations

import os
import json
import subprocess
import time
from pathlib import Path

from modules.remove import process_mol2_file
from modules.utils import detect_formal_charge_from_mol2


def _run(cmd, cwd: Path, log_file: Path):
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(cwd))
    log_file.write_text(result.stdout + "\n" + result.stderr)
    return result.returncode == 0


def _read_net_charge_for_residue(residue_dir: Path, mol2_path: Path) -> tuple[int, str]:
    """
    Priority:
      1) residue_meta.json net_charge
      2) detect formal charge from MOL2 (RDKit) as a fallback
      3) 0 as last fallback
    """
    meta = residue_dir / "residue_meta.json"
    if meta.exists():
        try:
            data = json.loads(meta.read_text(encoding="utf-8"))
            if "net_charge" in data:
                return int(data["net_charge"]), "meta"
        except Exception:
            pass

    detected = detect_formal_charge_from_mol2(str(mol2_path))
    if detected is not None:
        return int(detected), "detected_formal"

    return 0, "fallback_0"


def run_antechamber_for_all(
    mol2_files: list[str],
    backbone: str = "ff19SB",
    sidechain: str = "gaff2",
    charge: str = "bcc",   # kept for compatibility, but AC step will use -c rc to preserve charges
    generate_gmx: bool = False,
):
    amberhome = os.environ.get("AMBERHOME")
    if not amberhome:
        print("AMBERHOME not set.")
        return

    for input_mol2 in mol2_files:
        start_time = time.perf_counter()

        mol2_path = Path(input_mol2).resolve()
        residue_dir = mol2_path.parent
        resname = mol2_path.stem.upper()

        log_file = residue_dir / f"{resname}.log"

        ac_file = residue_dir / f"{resname}.ac"
        mc_file = residue_dir / f"{resname}.mc"
        prepin_file = residue_dir / f"{resname}.prepin"
        lib_file = residue_dir / f"{resname}.lib"
        frcmod_gaff = residue_dir / f"{resname}_{sidechain}.frcmod"
        frcmod_backbone = residue_dir / f"{resname}_{backbone}.frcmod"

        net_charge, charge_source = _read_net_charge_for_residue(residue_dir, mol2_path)
        print(f"[{resname}] using net charge = {net_charge} (source: {charge_source})")

        # 1) AC generation
        
        cmd = [
            "antechamber",
            "-fi", "mol2",
            "-i", mol2_path.name,
            "-bk", resname,
            "-fo", "ac",
            "-o", ac_file.name,
                      # <-- preserve existing charges
            "-at", "amber",
            "-nc", str(net_charge),
        ]
        if not _run(cmd, residue_dir, log_file):
            print(f"{resname} failed at AC generation.")
            continue

        # 2) MC generation
        # IMPORTANT: CHARGE in the mc file must match net_charge (not hardcoded 0).
        process_mol2_file(
            str(mol2_path),
            str(mc_file),
            head_name="N",
            tail_name="C",
            main_chain=["CA"],
            charge=float(net_charge),   # <-- FIX
        )

        # 3) prepgen
        cmd = [
            "prepgen",
            "-i", ac_file.name,
            "-o", prepin_file.name,
            "-m", mc_file.name,
            "-rn", resname,
        ]
        if not _run(cmd, residue_dir, log_file):
            print(f"{resname} failed at prepgen.")
            continue

        # 4) tleap
        leap_script = residue_dir / "leap.in"
        leap_script.write_text(
            f"""
source leaprc.{sidechain}
source leaprc.protein.{backbone}
{resname} = loadmol2 {mol2_path.name}
saveoff {resname} {lib_file.name}
quit
""".strip()
            + "\n",
            encoding="utf-8",
        )

        if not _run(["tleap", "-f", leap_script.name], residue_dir, log_file):
            print(f"{resname} failed at tleap.")
            continue

        # 5) parmchk2
        _run(["parmchk2", "-i", ac_file.name, "-f", "ac", "-o", frcmod_gaff.name], residue_dir, log_file)
        _run(["parmchk2", "-i", ac_file.name, "-f", "ac", "-o", frcmod_backbone.name], residue_dir, log_file)

        total_time = time.perf_counter() - start_time

        print(f"\033[1m\n{resname} parametrization complete\033[0m")
        print(f"Time taken: {total_time:.2f} s\n")