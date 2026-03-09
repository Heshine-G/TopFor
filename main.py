from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

from modules.processor import NonStandardAminoAcidProcessor
from modules.run_antechamber import run_antechamber_for_all


def load_residue_map(path: str | None) -> dict:
    if not path:
        return {}

    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Mapping file not found: {p}")

    data = json.loads(p.read_text(encoding="utf-8")) or {}
    out: dict[str, dict] = {}
    for k, v in data.items():
        if isinstance(v, dict):
            out[str(k).upper()] = v
    return out


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="nsaa-paramgen: NSAA AMBER parameter generation (ff19SB/ff14SB + gaff/gaff2)."
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-i", "--input", help="Single non-standard residue (.mol2 or .pdb)")
    group.add_argument(
        "-p",
        "--peptide",
        help="Whole peptide input (.pdb or .mol2). The tool extracts each non-standard residue and parameterizes each separately.",
    )
    group.add_argument(
        "-b",
        "--batch",
        help=(
            "Batch mode: either a text file of paths, OR a directory (process all *.mol2/*.pdb), "
            "OR a glob/prefix pattern like 'de*' or '*.mol2'."
        ),
    )

    parser.add_argument("--backbone", "-bb", choices=["ff14SB", "ff19SB", "ff99SB"], default="ff19SB")
    parser.add_argument("--sidechain", "-sc", choices=["gaff", "gaff2"], default="gaff2")
    parser.add_argument("--charge", "-c", choices=["gas", "bcc"], default="bcc")
    parser.add_argument("--gmx", "-gmx", action="store_true", help="Also generate GROMACS files (best-effort).")
    parser.add_argument(
        "--map",
        default=None,
        help=(
            "JSON mapping residue -> overrides. Example: "
            "{'ACE': {'head': null, 'tail': 'C', 'mainchain': ['CH3', 'C']}, "
            " 'ABA': {'head': 'N', 'tail': 'C', 'mainchain': ['CB', 'CA']}}"
        ),
    )
    parser.add_argument(
        "--default-net-charge",
        type=int,
        default=0,
        help="Fallback integer net charge used if not specified in --map and --net-charge is not set. Default 0.",
    )
    parser.add_argument(
        "--net-charge",
        "-nc",
        type=int,
        default=None,
        help="Global integer net charge override (applied to all residues unless --map provides a per-residue value).",
    )
    parser.add_argument(
        "--out",
        "-o",
        default=".",
        help="Output base directory. Residue outputs go to <out>/<RES>/. Default: current directory.",
    )

    return parser


def expand_batch_argument(batch_arg: str) -> list[str]:
    p = Path(batch_arg)

    if p.exists() and p.is_file():
        items: list[str] = []
        with p.open("r", encoding="utf-8") as f:
            for line in f:
                s = line.strip()
                if s and not s.startswith("#"):
                    items.append(s)
        return items

    if p.exists() and p.is_dir():
        files = sorted(str(x) for x in p.glob("*.mol2"))
        files.extend(sorted(str(x) for x in p.glob("*.pdb")))
        return files

    return sorted(glob.glob(batch_arg))


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()

    residue_map = load_residue_map(args.map)
    out_base = Path(args.out).resolve()
    out_base.mkdir(parents=True, exist_ok=True)

    charged_files: list[str] = []

    def build_processor(path: str) -> NonStandardAminoAcidProcessor:
        return NonStandardAminoAcidProcessor(
            input_file=path,
            charge_model=args.charge,
            residue_map=residue_map,
            default_net_charge=args.default_net_charge,
            net_charge_override=args.net_charge,
            output_base=str(out_base),
        )

    def process_one(path: str) -> None:
        proc = build_processor(path)
        charged_files.extend(proc.process_single_residue())

    if args.input:
        process_one(args.input)

    elif args.peptide:
        proc = build_processor(args.peptide)
        charged_files.extend(proc.process_peptide())

    elif args.batch:
        batch_items = expand_batch_argument(args.batch)
        if not batch_items:
            print(f"Batch mode: nothing matched: {args.batch}")
            return

        for fp in batch_items:
            process_one(fp)

    if charged_files:
        run_antechamber_for_all(
            charged_files,
            backbone=args.backbone,
            sidechain=args.sidechain,
            charge=args.charge,
            generate_gmx=args.gmx,
        )
    else:
        print("No charged .mol2 files were generated. Skipping parameter generation.")


if __name__ == "__main__":
    main()