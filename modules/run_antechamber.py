from __future__ import annotations

import json
import os
import subprocess
import time
from pathlib import Path

from modules.remove import process_mol2_file
from modules.utils import (
    classify_residue_net_charge,
    normalize_resname,
    fix_backbone_atom_types_in_ac,
    fix_backbone_atom_types_in_prepin,
)


def _run(cmd, cwd: Path, log_file: Path) -> bool:
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


def _read_json_if_exists(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def _read_residue_meta(residue_dir: Path) -> dict:
    residue_meta = _read_json_if_exists(residue_dir / "residue_meta.json")
    capping_meta = _read_json_if_exists(residue_dir / "residue_capping_meta.json")

    if residue_meta:
        if capping_meta:
            residue_meta.setdefault("requested_head_name", capping_meta.get("requested_head_name"))
            residue_meta.setdefault("requested_tail_name", capping_meta.get("requested_tail_name"))
            residue_meta.setdefault("has_head", capping_meta.get("has_head"))
            residue_meta.setdefault("has_tail", capping_meta.get("has_tail"))
            residue_meta.setdefault("applied_caps", capping_meta.get("applied_caps", []))
        return residue_meta

    if capping_meta:
        return {
            "head_name": capping_meta.get("requested_head_name"),
            "tail_name": capping_meta.get("requested_tail_name"),
            "main_chain": None,
            "pre_head_type": "C",
            "post_tail_type": "N",
            "applied_caps": capping_meta.get("applied_caps", []),
        }

    return {}


def _read_net_charge_for_residue(residue_dir: Path, mol2_path: Path) -> tuple[int, str]:
    data = _read_residue_meta(residue_dir)
    if "net_charge" in data:
        return int(data["net_charge"]), "meta"

    classified, source = classify_residue_net_charge(str(mol2_path), resname=normalize_resname(mol2_path.stem))
    if classified is not None:
        return int(classified), source

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

    backbone_parm_path = os.path.join(amberhome, "dat", "leap", "parm", backbone_map[backbone])
    gaff_parm_path = os.path.join(amberhome, "dat", "leap", "parm", sidechain_map[sidechain])

    for input_mol2 in mol2_files:
        #start_time = time.perf_counter()

        mol2_path = Path(input_mol2).resolve()
        residue_dir = mol2_path.parent
        resname = normalize_resname(mol2_path.stem)
        meta = _read_residue_meta(residue_dir)

        log_file = residue_dir / f"{resname}.log"
        log_file.write_text("", encoding="utf-8")

        ac_output = residue_dir / f"{resname}.ac"
        mc_file = residue_dir / f"{resname}.mc"
        prepin_file = residue_dir / f"{resname}.prepin"
        lib_file = residue_dir / f"{resname}.lib"

        backbone_frcmod_output = residue_dir / f"{resname}_{backbone}.frcmod"
        gaff_frcmod_output = residue_dir / f"{resname}_{sidechain}.frcmod"

        net_charge, charge_source = _read_net_charge_for_residue(residue_dir, mol2_path)
        print(f"[{resname}] using net charge = {net_charge}")

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

        fix_backbone_atom_types_in_ac(str(ac_output))

        applied_caps = tuple(str(x).upper() for x in meta.get("applied_caps", []))
        process_mol2_file(
            str(mol2_path),
            str(mc_file),
            head_name=meta.get("head_name", "N"),
            tail_name=meta.get("tail_name", "C"),
            main_chain=meta.get("main_chain"),
            charge=float(net_charge),
            central_resname=resname,
            cap_resnames=applied_caps if applied_caps else ("ACE", "NME"),
            pre_head_type=str(meta.get("pre_head_type", "C")),
            post_tail_type=str(meta.get("post_tail_type", "N")),
            infer_mainchain_from_connectivity=True,
        )

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

        if fix_backbone_atom_types_in_prepin(str(prepin_file)):
            print(f"[{resname}] corrected backbone atom types in PREPIN file")

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

        #total_time = time.perf_counter() - start_time
        print(f"\033[1m\n{resname} parametrization complete\033[0m")
        #print(f"Time taken: {total_time:.2f} s\n")
