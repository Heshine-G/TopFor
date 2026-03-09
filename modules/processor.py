from __future__ import annotations

import json
import subprocess
import re
from pathlib import Path
from datetime import datetime

from modules.utils import (
    detect_formal_charge_from_mol2,
    renormalize_mol2_partial_charges_to_integer,
)

from modules.split_nonstandard_residues import extract_nonstandard_residues


class NonStandardAminoAcidProcessor:

    def __init__(
        self,
        input_file: str,
        charge_model: str = "bcc",
        residue_map: dict | None = None,
        default_net_charge: int = 0,
        net_charge_override: int | None = None,
        output_base: str = ".",
    ):

        self.input_file = input_file
        self.charge_model = charge_model
        self.residue_map = residue_map or {}
        self.default_net_charge = default_net_charge
        self.net_charge_override = net_charge_override
        self.output_base = Path(output_base).resolve()

    def _run(self, cmd: list[str], cwd: Path | None = None):

        return subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=str(cwd) if cwd else None,
        )

    # --------------------------------------------------
    # NEW ROBUST RESIDUE NAME DETECTION
    # --------------------------------------------------

    def _get_residue_name(self, file_path: str) -> str:
        """
        Determine residue name robustly.

        Priority:
        1. MOL2 SUBSTRUCTURE section
        2. atom subst_name
        3. filename fallback
        """

        p = Path(file_path)

        if p.suffix.lower() == ".mol2":

            with open(file_path, "r", encoding="utf-8", errors="ignore") as f:

                lines = f.readlines()

            in_sub = False

            for line in lines:

                if line.startswith("@<TRIPOS>SUBSTRUCTURE"):
                    in_sub = True
                    continue

                if line.startswith("@<TRIPOS>") and in_sub:
                    break

                if in_sub:

                    parts = line.split()

                    if len(parts) >= 2:
                        return parts[1].split(".")[0].upper()

            # fallback: read subst_name from atoms

            in_atoms = False

            for line in lines:

                if line.startswith("@<TRIPOS>ATOM"):
                    in_atoms = True
                    continue

                if line.startswith("@<TRIPOS>") and in_atoms:
                    break

                if in_atoms:

                    parts = line.split()

                    if len(parts) >= 8:
                        return parts[7].split(".")[0].upper()

        # filename fallback

        stem = Path(file_path).stem

        stem = re.sub(r"\d+$", "", stem)

        stem = stem.split("_")[0]

        return stem.upper()

    # --------------------------------------------------

    def _normalize_nullable_name(self, value: object, default: str | None) -> str | None:

        if value is None:
            return default

        s = str(value).strip()

        if not s:
            return default

        if s.upper() in {"NONE", "NULL", "0"}:
            return None

        return s

    # --------------------------------------------------

    def _get_residue_cfg(self, resname: str) -> dict:

        raw = self.residue_map.get(resname, {})

        head_name = self._normalize_nullable_name(raw.get("head", "N"), "N")
        tail_name = self._normalize_nullable_name(raw.get("tail", "C"), "C")

        main_chain = raw.get("mainchain")

        if not isinstance(main_chain, list):
            main_chain = None

        pre_head_type = str(raw.get("pre_head_type", "C"))
        post_tail_type = str(raw.get("post_tail_type", "N"))

        return {
            "head_name": head_name,
            "tail_name": tail_name,
            "main_chain": main_chain,
            "pre_head_type": pre_head_type,
            "post_tail_type": post_tail_type,
        }

    # --------------------------------------------------

    def _prepare_input_mol2(self, input_path: Path, residue_dir: Path) -> Path:

        residue_file = residue_dir / "residue.mol2"

        suffix = input_path.suffix.lower()

        if suffix == ".mol2":

            residue_file.write_text(
                input_path.read_text(encoding="utf-8", errors="replace"),
                encoding="utf-8",
            )

            return residue_file

        if suffix == ".pdb":

            convert_script = Path(__file__).parent / "pdb_to_mol2.py"

            result = self._run(
                ["python", str(convert_script), str(input_path), str(residue_file)],
                cwd=residue_dir,
            )

            if result.returncode != 0:
                raise RuntimeError("PDB conversion failed")

            return residue_file

        raise ValueError("Unsupported format")

    # --------------------------------------------------

    def _resolve_net_charge_after_capping(self, resname: str, capped_mol2_path: str):

        if resname in self.residue_map and "net_charge" in self.residue_map[resname]:
            return int(self.residue_map[resname]["net_charge"]), "map"

        if self.net_charge_override is not None:
            return int(self.net_charge_override), "cli"

        detected = detect_formal_charge_from_mol2(capped_mol2_path)

        if detected is not None:
            return int(detected), "detected_formal"

        return int(self.default_net_charge), "default"

    # --------------------------------------------------

    def _process_input_path(self, input_path: Path):

        resname = self._get_residue_name(str(input_path))

        residue_dir = self.output_base / resname

        residue_dir.mkdir(parents=True, exist_ok=True)

        cfg = self._get_residue_cfg(resname)

        residue_file = self._prepare_input_mol2(input_path, residue_dir)

        cap_script = Path(__file__).parent / "capping.py"

        head_arg = cfg["head_name"] if cfg["head_name"] else "NONE"
        tail_arg = cfg["tail_name"] if cfg["tail_name"] else "NONE"

        self._run(
            ["python", str(cap_script), str(residue_dir), str(head_arg), str(tail_arg)],
            cwd=residue_dir,
        )

        capped_file = residue_dir / "residue_capped.mol2"

        net_charge, charge_source = self._resolve_net_charge_after_capping(
            resname, str(capped_file)
        )

        print(f"[{resname}] net charge = {net_charge} (source: {charge_source})")

        charged_file = residue_dir / f"{resname}.mol2"

        cmd = [
            "antechamber",
            "-fi", "mol2",
            "-i", str(capped_file),
            "-bk", resname,
            "-fo", "mol2",
            "-o", str(charged_file),
            "-c", self.charge_model,
            "-at", "amber",
            "-nc", str(net_charge),
        ]

        self._run(cmd, cwd=residue_dir)

        renormalize_mol2_partial_charges_to_integer(str(charged_file), int(net_charge))

        return [str(charged_file)]

    # --------------------------------------------------

    def process_single_residue(self):

        return self._process_input_path(Path(self.input_file).resolve())

    # --------------------------------------------------

    def process_peptide(self):

        peptide_path = Path(self.input_file).resolve()

        split_dir = self.output_base / "_extracted_nonstandard_residues"

        split_dir.mkdir(parents=True, exist_ok=True)

        extracted_paths = extract_nonstandard_residues(
            str(peptide_path),
            str(split_dir),
        )

        all_charged = []

        for residue in extracted_paths:

            charged = self._process_input_path(Path(residue).resolve())

            all_charged.extend(charged)

        return all_charged