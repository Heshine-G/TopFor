from __future__ import annotations

import json
import subprocess
from pathlib import Path

from modules.split_nonstandard_residues import extract_nonstandard_residues
from modules.utils import (
    classify_residue_net_charge,
    normalize_resname,
    renormalize_mol2_partial_charges_to_integer,
    validate_molecule,
)


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

    def _run(self, cmd: list[str], cwd: Path | None = None) -> subprocess.CompletedProcess:
        return subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=str(cwd) if cwd else None,
        )

    def _get_residue_name(self, file_path: str) -> str:
        p = Path(file_path)

        if p.suffix.lower() == ".mol2":
            lines = p.read_text(encoding="utf-8", errors="ignore").splitlines()

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
                        return normalize_resname(parts[1])

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
                        return normalize_resname(parts[7])

        return normalize_resname(p.stem)

    def _normalize_nullable_name(self, value: object, default: str | None) -> str | None:
        if value is None:
            return default
        s = str(value).strip()
        if not s:
            return default
        if s.upper() in {"NONE", "NULL", "0"}:
            return None
        return s

    def _read_split_meta(self, input_path: Path) -> dict:
        candidates = [
            input_path.with_suffix(".split.json"),
            input_path.parent / f"{input_path.stem}.split.json",
        ]
        for meta_path in candidates:
            if not meta_path.exists():
                continue
            try:
                data = json.loads(meta_path.read_text(encoding="utf-8"))
                if isinstance(data, dict):
                    return data
            except Exception:
                continue
        return {}

    def _get_residue_cfg(self, resname: str, input_path: Path | None = None) -> dict:
        raw = dict(self.residue_map.get(resname, {}))
        split_meta = self._read_split_meta(input_path) if input_path else {}

        topology = str(split_meta.get("topology", "")).lower()
        default_head = split_meta.get("head_name", "N")
        default_tail = split_meta.get("tail_name", "C")

        if topology == "n_term_like":
            default_head = None
        elif topology == "c_term_like":
            default_tail = None

        head_raw = raw.get("head", default_head)
        tail_raw = raw.get("tail", default_tail)

        head_name = self._normalize_nullable_name(head_raw, default_head)
        tail_name = self._normalize_nullable_name(tail_raw, default_tail)

        main_chain = raw.get("mainchain")
        if not isinstance(main_chain, list):
            main_chain = split_meta.get("main_chain")
            if not isinstance(main_chain, list):
                main_chain = None

        pre_head_type = str(raw.get("pre_head_type", split_meta.get("pre_head_type", "C")))
        post_tail_type = str(raw.get("post_tail_type", split_meta.get("post_tail_type", "N")))

        return {
            "head_name": head_name,
            "tail_name": tail_name,
            "main_chain": main_chain,
            "pre_head_type": pre_head_type,
            "post_tail_type": post_tail_type,
            "split_meta": split_meta,
        }

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
                raise RuntimeError(
                    "PDB conversion failed\n"
                    f"STDOUT:\n{result.stdout}\n"
                    f"STDERR:\n{result.stderr}"
                )
            return residue_file

        raise ValueError(f"Unsupported format: {input_path.suffix}")

    def _resolve_input_net_charge(self, resname: str, input_mol2_path: str, split_meta: dict) -> tuple[int, str]:
        if resname in self.residue_map and "net_charge" in self.residue_map[resname]:
            return int(self.residue_map[resname]["net_charge"]), "map"

        if self.net_charge_override is not None:
            return int(self.net_charge_override), "cli"

        full_context = split_meta.get("full_context_net_charge")
        if full_context is not None:
            return int(full_context), str(split_meta.get("full_context_charge_source", "split_full_context"))

        classified, source = classify_residue_net_charge(input_mol2_path, resname=resname)
        if classified is not None:
            return int(classified), source

        return int(self.default_net_charge), "default"

    def _write_residue_meta(
        self,
        residue_dir: Path,
        *,
        resname: str,
        cfg: dict,
        capping_meta: dict,
        net_charge: int,
        charge_source: str,
        validation_warnings: list[str],
    ) -> None:
        split_meta = cfg.get("split_meta", {}) or {}
        meta = {
            "resname": resname,
            "head_name": cfg.get("head_name"),
            "tail_name": cfg.get("tail_name"),
            "main_chain": cfg.get("main_chain"),
            "pre_head_type": cfg.get("pre_head_type", "C"),
            "post_tail_type": cfg.get("post_tail_type", "N"),
            "requested_head_name": capping_meta.get("requested_head_name", cfg.get("head_name")),
            "requested_tail_name": capping_meta.get("requested_tail_name", cfg.get("tail_name")),
            "has_head": bool(capping_meta.get("has_head", False)),
            "has_tail": bool(capping_meta.get("has_tail", False)),
            "applied_caps": list(capping_meta.get("applied_caps", [])),
            "net_charge": int(net_charge),
            "net_charge_source": charge_source,
            "validation_warnings": validation_warnings,
            "source_input": split_meta.get("source_input"),
            "source_subst_id": split_meta.get("source_subst_id"),
            "source_subst_name": split_meta.get("source_subst_name"),
            "external_bonds": split_meta.get("external_bonds", []),
            "is_polymer_internal": bool(split_meta.get("is_polymer_internal", False)),
            "is_n_terminal_like": bool(split_meta.get("is_n_terminal_like", False)),
            "is_c_terminal_like": bool(split_meta.get("is_c_terminal_like", False)),
            "is_cyclic_like": bool(split_meta.get("is_cyclic_like", False)),
            "topology": split_meta.get("topology"),
            "full_context_net_charge": split_meta.get("full_context_net_charge"),
            "full_context_charge_source": split_meta.get("full_context_charge_source"),
        }
        (residue_dir / "residue_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    def _process_input_path(self, input_path: Path):
        resname = self._get_residue_name(str(input_path))
        residue_dir = self.output_base / resname
        residue_dir.mkdir(parents=True, exist_ok=True)

        cfg = self._get_residue_cfg(resname, input_path=input_path)
        input_mol2 = self._prepare_input_mol2(input_path, residue_dir)

        net_charge, charge_source = self._resolve_input_net_charge(
            resname,
            str(input_mol2),
            cfg.get("split_meta", {}) or {},
        )

        input_warnings = validate_molecule(str(input_mol2), expected_charge=net_charge)

        cap_script = Path(__file__).parent / "capping.py"
        head_arg = cfg["head_name"] if cfg["head_name"] else "NONE"
        tail_arg = cfg["tail_name"] if cfg["tail_name"] else "NONE"

        cap_result = self._run(
            ["python", str(cap_script), str(residue_dir), str(head_arg), str(tail_arg)],
            cwd=residue_dir,
        )
        if cap_result.returncode != 0:
            raise RuntimeError(
                f"Capping failed for {resname}\n"
                f"STDOUT:\n{cap_result.stdout}\n"
                f"STDERR:\n{cap_result.stderr}"
            )

        capped_file = residue_dir / "residue_capped.mol2"
        if not capped_file.exists() or capped_file.stat().st_size == 0:
            raise RuntimeError(f"Capped MOL2 was not created for {resname}: {capped_file}")

        capping_meta_file = residue_dir / "residue_capping_meta.json"
        try:
            capping_meta = json.loads(capping_meta_file.read_text(encoding="utf-8"))
        except Exception:
            capping_meta = {}

        validation_warnings = input_warnings + validate_molecule(str(capped_file), expected_charge=net_charge)
        print(f"[{resname}] net charge = {net_charge} (source: {charge_source})")
        for w in validation_warnings:
            print(f"[{resname}] warning: {w}")

        self._write_residue_meta(
            residue_dir,
            resname=resname,
            cfg=cfg,
            capping_meta=capping_meta,
            net_charge=net_charge,
            charge_source=charge_source,
            validation_warnings=validation_warnings,
        )

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

        charge_result = self._run(cmd, cwd=residue_dir)
        if charge_result.returncode != 0:
            raise RuntimeError(
                f"antechamber charge assignment failed for {resname}\n"
                f"STDOUT:\n{charge_result.stdout}\n"
                f"STDERR:\n{charge_result.stderr}"
            )

        renormalize_mol2_partial_charges_to_integer(str(charged_file), int(net_charge))

        from modules.utils import fix_backbone_atom_types
        fix_backbone_atom_types(str(charged_file))
        return [str(charged_file)]

    def process_single_residue(self):
        return self._process_input_path(Path(self.input_file).resolve())

    def process_peptide(self):
        peptide_path = Path(self.input_file).resolve()
        split_dir = self.output_base / "_extracted_nonstandard_residues"
        split_dir.mkdir(parents=True, exist_ok=True)

        extracted_paths = extract_nonstandard_residues(str(peptide_path), str(split_dir))

        all_charged = []
        for residue in extracted_paths:
            charged = self._process_input_path(Path(residue).resolve())
            all_charged.extend(charged)

        return all_charged
