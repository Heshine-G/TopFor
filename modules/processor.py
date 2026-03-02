# modules/processor.py
from __future__ import annotations

import json
import subprocess
from pathlib import Path
from datetime import datetime

from modules.utils import (
    detect_formal_charge_from_mol2,
    renormalize_mol2_partial_charges_to_integer,
)


class NonStandardAminoAcidProcessor:
    """
    Step 1:
    Creates a clean residue folder:

      <OUT>/<RES>/

    and generates:
      residue.mol2
      residue_capped.mol2
      <RES>.mol2  (charged)
      residue_meta.json
    """

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
        return subprocess.run(cmd, capture_output=True, text=True, cwd=str(cwd) if cwd else None)

    def _get_residue_name(self, file_path: str) -> str:
        return Path(file_path).stem.split("_")[0].upper()

    def _resolve_net_charge_after_capping(self, resname: str, capped_mol2_path: str) -> tuple[int, str]:
        """
        Charge priority order:
          1) per-residue JSON map: residue_map[RES]["net_charge"]
          2) global CLI override: --net-charge
          3) auto-detect formal charge from capped MOL2 (RDKit)
          4) fallback default: --default-net-charge
        Returns: (net_charge, source_string)
        """
        overrides = self.residue_map.get(resname, {})
        if isinstance(overrides, dict) and "net_charge" in overrides:
            return int(overrides["net_charge"]), "map"

        if self.net_charge_override is not None:
            return int(self.net_charge_override), "cli"

        detected = detect_formal_charge_from_mol2(capped_mol2_path)
        if detected is not None:
            return int(detected), "detected_formal"

        return int(self.default_net_charge), "default"

    def process_single_residue(self) -> list[str]:
        input_path = Path(self.input_file).resolve()
        if not input_path.exists():
            print("Input not found.")
            return []

        resname = self._get_residue_name(str(input_path))
        residue_dir = self.output_base / resname
        residue_dir.mkdir(parents=True, exist_ok=True)

        # Copy original mol2
        residue_file = residue_dir / "residue.mol2"
        residue_file.write_text(input_path.read_text(encoding="utf-8", errors="replace"))

        # Run PyMOL capping
        cap_script = Path(__file__).parent / "capping.py"
        result = self._run(["python", str(cap_script), str(residue_dir)])
        if result.returncode != 0:
            print("Capping failed.")
            return []

        capped_file = residue_dir / "residue_capped.mol2"
        if not capped_file.exists():
            print("Capping failed: residue_capped.mol2 not found.")
            return []

        # Decide net charge (map > cli > detect_formal > default)
        net_charge, charge_source = self._resolve_net_charge_after_capping(resname, str(capped_file))
        print(f"[{resname}] net charge = {net_charge} (source: {charge_source})")

        # Charge generation (antechamber writes a new MOL2 with charges)
        charged_file = residue_dir / f"{resname}.mol2"

        cmd = [
            "antechamber",
            "-fi", "mol2",
            "-i", str(capped_file),
            "-bk", resname,
            "-fo", "mol2",
            "-o", str(charged_file),
            "-c", self.charge_model,   # bcc/gas
            "-at", "amber",
            "-nc", str(net_charge),
        ]
        result = self._run(cmd, cwd=residue_dir)
        if result.returncode != 0:
            print("Charge generation failed.")
            # show stderr for debugging
            if result.stderr:
                print(result.stderr)
            return []

        # Renormalize MOL2 partial charges so sum == integer net charge
        try:
            changed = renormalize_mol2_partial_charges_to_integer(str(charged_file), int(net_charge))
            if changed:
                print(f"[{resname}] partial charges renormalized to sum exactly {net_charge}.")
        except Exception as e:
            print(f"[{resname}] WARNING: charge renormalization failed: {e}")

        meta = {
            "resname": resname,
            "net_charge": int(net_charge),
            "net_charge_source": charge_source,
            "charge_model": self.charge_model,
            "timestamp": datetime.now().isoformat(),
        }
        (residue_dir / "residue_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

        return [str(charged_file)]