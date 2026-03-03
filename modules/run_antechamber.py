from __future__ import annotations

import os
import json
import subprocess
import time
from pathlib import Path

from modules.remove import process_mol2_file
from modules.utils import detect_formal_charge_from_mol2


def _run(cmd, cwd: Path, log_file: Path) -> bool:
    """Run a command and append stdout/stderr to log file."""
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(cwd))

    with log_file.open("a", encoding="utf-8") as f:
        f.write("\n\n=== CMD ===\n")
        f.write(" ".join(cmd) if isinstance(cmd, list) else str(cmd))
        f.write("\n=== STDOUT ===\n")
        f.write(result.stdout)
        f.write("\n=== STDERR ===\n")
        f.write(result.stderr)
        f.write("\n=== RETURN CODE ===\n")
        f.write(str(result.returncode) + "\n")

    return result.returncode == 0


def _read_net_charge_for_residue(residue_dir: Path, mol2_path: Path) -> tuple[int, str]:
    """
    Priority:
        1) residue_meta.json
        2) RDKit formal charge detection
        3) fallback 0
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
    charge: str = "bcc",
    generate_gmx: bool = False,
):
    amberhome = os.environ.get("AMBERHOME")
    if not amberhome:
        print("AMBERHOME not set.")
        return

    # ------------------------------------------------------------
    # Proper mapping for parmchk2 base parameter files
    # ------------------------------------------------------------
    backbone_map = {
        "ff19SB": "parm19.dat",
        "ff14SB": "parm10.dat",
        "ff99SB": "parm99.dat",
    }

    sidechain_map = {
        "gaff": "gaff.dat",
        "gaff2": "gaff2.dat",
    }

    if backbone not in backbone_map:
        raise ValueError(f"Unsupported backbone: {backbone}")

    if sidechain not in sidechain_map:
        raise ValueError(f"Unsupported sidechain: {sidechain}")

    parm_file = backbone_map[backbone]
    gaff_file = sidechain_map[sidechain]

    backbone_parm_path = os.path.join(
        amberhome, "dat", "leap", "parm", parm_file
    )

    gaff_parm_path = os.path.join(
        amberhome, "dat", "leap", "parm", gaff_file
    )

    # ------------------------------------------------------------

    for input_mol2 in mol2_files:
        start_time = time.perf_counter()

        mol2_path = Path(input_mol2).resolve()
        residue_dir = mol2_path.parent
        resname = mol2_path.stem.upper()

        log_file = residue_dir / f"{resname}.log"
        log_file.write_text("", encoding="utf-8")

        ac_output = residue_dir / f"{resname}.ac"
        mc_file = residue_dir / f"{resname}.mc"
        prepin_file = residue_dir / f"{resname}.prepin"
        lib_file = residue_dir / f"{resname}.lib"

        backbone_frcmod_output = residue_dir / f"{resname}_{backbone}.frcmod"
        gaff_frcmod_output = residue_dir / f"{resname}_{sidechain}.frcmod"

        net_charge, charge_source = _read_net_charge_for_residue(residue_dir, mol2_path)
        print(f"[{resname}] using net charge = {net_charge} (source: {charge_source})")

        # ------------------------------------------------------------
        # AC generation (kept as requested)
        # ------------------------------------------------------------
        antechamber_cmd = [
            "antechamber",
            "-i", mol2_path.name,
            "-fi", "mol2",
            "-o", ac_output.name,
            "-fo", "ac",
            "-at", "amber",
            "-nc", str(net_charge),
        ]

        if not _run(antechamber_cmd, residue_dir, log_file):
            print(f"{resname} failed at AC generation.")
            continue

        # ------------------------------------------------------------
        # MC file generation
        # ------------------------------------------------------------
        process_mol2_file(
            str(mol2_path),
            str(mc_file),
            head_name="N",
            tail_name="C",
            main_chain=["CA"],
            charge=float(net_charge),
        )

        # ------------------------------------------------------------
        # PREPGEN
        # ------------------------------------------------------------
        prepgen_cmd = [
            "prepgen",
            "-i", ac_output.name,
            "-o", prepin_file.name,
            "-m", mc_file.name,
            "-rn", resname,
        ]

        if not _run(prepgen_cmd, residue_dir, log_file):
            print(f"{resname} failed at prepgen.")
            continue

        # ------------------------------------------------------------
        # FIXED parmchk2 section (ONLY corrected part)
        # ------------------------------------------------------------

        parmchk_backbone_cmd = [
            "parmchk2",
            "-i", ac_output.name,
            "-f", "ac",
            "-o", backbone_frcmod_output.name,
            "-a", "Y",
            "-p", backbone_parm_path,
        ]

        parmchk_sidechain_cmd = [
            "parmchk2",
            "-i", ac_output.name,
            "-f", "ac",
            "-o", gaff_frcmod_output.name,
            "-a", "Y",
            "-p", gaff_parm_path,
        ]

        if not _run(parmchk_backbone_cmd, residue_dir, log_file):
            print(f"{resname} failed at parmchk2 (backbone).")
            continue

        if not _run(parmchk_sidechain_cmd, residue_dir, log_file):
            print(f"{resname} failed at parmchk2 (sidechain).")
            continue

        # ------------------------------------------------------------
        # TLEAP
        # ------------------------------------------------------------
        leap_script = residue_dir / "leap.in"
        leap_script.write_text(
            (
                f"source leaprc.protein.{backbone}\n"
                f"source leaprc.{sidechain}\n"
                f"loadamberparams {backbone_frcmod_output.name}\n"
                f"loadamberparams {gaff_frcmod_output.name}\n"
                f"{resname} = loadmol2 {mol2_path.name}\n"
                f"saveoff {resname} {lib_file.name}\n"
                f"quit\n"
            ),
            encoding="utf-8",
        )

        if not _run(["tleap", "-f", leap_script.name], residue_dir, log_file):
            print(f"{resname} failed at tleap.")
            continue

        total_time = time.perf_counter() - start_time
        print(f"\033[1m\n{resname} parametrization complete\033[0m")
        print(f"Time taken: {total_time:.2f} s\n")