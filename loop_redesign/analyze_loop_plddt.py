#!/usr/bin/env python3
"""Score the fixed-length p35 loop (residues 147-160) from AlphaFold 3 models.

Unlike the helix campaign, this window does not move: ProteinMPNN redesigned
exactly B147-160 on the wild-type backbone.

Example:
  python analyze_loop_plddt.py \\
      --input-dir selected \\
      --output /tmp/loop_plddt.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path


def fix_long_path(path: Path) -> Path:
    text = str(path.resolve())
    if sys.platform == "win32" and len(text) > 260 and not text.startswith("\\\\?\\"):
        return Path("\\\\?\\" + text)
    return path


def extract_sequence(data_json: Path, chain_id: str = "B") -> str | None:
    data = json.loads(fix_long_path(data_json).read_text(encoding="utf-8"))
    for entry in data.get("sequences", []):
        protein = entry.get("protein", entry)
        if protein.get("id") == chain_id:
            return protein.get("sequence")
    return None


def extract_plddt(cif_file: Path, chain_id: str = "B", start_pos: int = 147, end_pos: int = 160):
    values = []
    with fix_long_path(cif_file).open(encoding="utf-8") as handle:
        in_atom = False
        for line in handle:
            if line.startswith("ATOM"):
                in_atom = True
            if not in_atom or not line.startswith("ATOM"):
                continue
            parts = line.split()
            if len(parts) < 15:
                continue
            try:
                if parts[6] != chain_id:
                    continue
                seq_id = int(parts[8])
                plddt = float(parts[14])
            except (ValueError, IndexError):
                continue
            if start_pos <= seq_id <= end_pos:
                values.append(plddt)
    if not values:
        return None, 0
    return sum(values) / len(values), len(values)


def extract_iptm(summary_json: Path) -> float | None:
    data = json.loads(fix_long_path(summary_json).read_text(encoding="utf-8"))
    return data.get("iptm")


def list_models(folder: Path) -> list[tuple[Path, Path]]:
    pairs = []
    sample_dirs = sorted(p for p in folder.iterdir() if p.is_dir() and p.name.startswith("seed-1_sample-"))
    search_roots = sample_dirs or [folder]
    for root in search_roots:
        cifs = list(root.glob("*_model.cif"))
        summaries = list(root.glob("*_summary_confidences.json"))
        if cifs and summaries:
            pairs.append((cifs[0], summaries[0]))
    return pairs


def analyze_folder(folder: Path, start_pos: int, end_pos: int) -> dict | None:
    data_files = list(folder.glob("*_data.json"))
    if not data_files:
        return None
    sequence = extract_sequence(data_files[0])
    if not sequence or len(sequence) < end_pos:
        return None

    plddts = []
    iptms = []
    atom_counts = []
    for cif, summary in list_models(folder):
        plddt, n_atoms = extract_plddt(cif, start_pos=start_pos, end_pos=end_pos)
        iptm = extract_iptm(summary)
        if plddt is not None:
            plddts.append(plddt)
            atom_counts.append(n_atoms)
        if iptm is not None:
            iptms.append(iptm)

    loop = sequence[start_pos - 1 : end_pos]
    return {
        "mutant_name": folder.name,
        "aa_count": end_pos - start_pos + 1,
        "start_pos": start_pos,
        "end_pos": end_pos,
        "avg_plddt": sum(plddts) / len(plddts) if plddts else None,
        "avg_iptm": sum(iptms) / len(iptms) if iptms else None,
        "avg_atom_count": sum(atom_counts) / len(atom_counts) if atom_counts else None,
        "target_subsequence": loop,
        "sequence": sequence,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True, help="Directory of AF3 design folders")
    parser.add_argument("--output", default="mIL12_redesignloop_147-160_score.csv")
    parser.add_argument("--start", type=int, default=147)
    parser.add_argument("--end", type=int, default=160)
    args = parser.parse_args()

    root = Path(args.input_dir)
    folders = sorted(p for p in root.iterdir() if p.is_dir())
    results = [analyze_folder(folder, args.start, args.end) for folder in folders]
    results = [row for row in results if row]

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "mutant_name", "aa_count", "start_pos", "end_pos",
        "avg_plddt", "avg_iptm", "avg_atom_count",
        "target_subsequence", "sequence",
    ]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print(f"Wrote {len(results)} rows to {output}")


if __name__ == "__main__":
    main()
