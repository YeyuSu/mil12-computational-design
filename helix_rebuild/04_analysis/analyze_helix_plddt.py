#!/usr/bin/env python3
"""Score designed p35 helices from AlphaFold 3 predictions.

For each design, find the DPYRVK motif after residue 147 on chain B and average
pLDDT over residues [147, residue before DPYRVK] across available models.

Supports:
  - full AF3 folders with seed-1_sample-* subdirectories
  - slim folders in this repository that only keep the top-ranked *_model.cif

Example:
  python analyze_helix_plddt.py \\
      --input-dir ../03_alphafold3/selected \\
      --output scores/selected_helix_plddt.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from multiprocessing import Pool, cpu_count
from pathlib import Path


def fix_long_path(path: Path) -> Path:
    text = str(path.resolve())
    if sys.platform == "win32" and len(text) > 260 and not text.startswith("\\\\?\\"):
        return Path("\\\\?\\" + text)
    return path


def find_helix_end(sequence: str, start_pos: int = 147, min_length: int = 14, max_length: int = 25):
    """Return 1-based last residue of the designed helix (residue before DPYRVK)."""
    motif = "DPYRVK"
    start_idx = start_pos - 1
    for end_idx in range(start_idx + min_length - 1, min(start_idx + max_length, len(sequence))):
        if sequence[end_idx + 1 : end_idx + 7] == motif:
            return end_idx + 1
    return None


def extract_sequence(data_json: Path, chain_id: str = "B") -> str | None:
    data = json.loads(fix_long_path(data_json).read_text(encoding="utf-8"))
    for entry in data.get("sequences", []):
        protein = entry.get("protein", entry)
        if protein.get("id") == chain_id:
            return protein.get("sequence")
    if data.get("sequences"):
        protein = data["sequences"][-1].get("protein", data["sequences"][-1])
        return protein.get("sequence")
    return None


def extract_plddt(cif_file: Path, chain_id: str = "B", start_pos: int = 147, end_pos: int | None = None):
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
            # AF3 mmCIF: residue name, chain, entity, seq id ... pLDDT near the end
            try:
                if parts[6] == chain_id:
                    seq_id = int(parts[8])
                else:
                    continue
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


def analyze_folder(folder: Path) -> dict | None:
    data_files = list(folder.glob("*_data.json"))
    if not data_files:
        return None
    sequence = extract_sequence(data_files[0])
    if not sequence:
        return None
    end_pos = find_helix_end(sequence)
    if end_pos is None:
        print(f"Skip {folder.name}: DPYRVK motif not found after residue 147")
        return None

    plddts = []
    iptms = []
    atom_counts = []
    for cif, summary in list_models(folder):
        plddt, n_atoms = extract_plddt(cif, end_pos=end_pos)
        iptm = extract_iptm(summary)
        if plddt is not None:
            plddts.append(plddt)
            atom_counts.append(n_atoms)
        if iptm is not None:
            iptms.append(iptm)

    helix = sequence[146:end_pos]
    return {
        "mutant_name": folder.name,
        "aa_count": end_pos - 147 + 1,
        "start_pos": 147,
        "end_pos": end_pos,
        "avg_plddt": sum(plddts) / len(plddts) if plddts else None,
        "avg_iptm": sum(iptms) / len(iptms) if iptms else None,
        "avg_atom_count": sum(atom_counts) / len(atom_counts) if atom_counts else None,
        "target_subsequence": helix,
        "sequence": sequence,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True, help="Directory of AF3 design folders")
    parser.add_argument("--output", default="mIL12_alphahelix_dynamic_score.csv")
    parser.add_argument("--ncores", type=int, default=1)
    parser.add_argument(
        "--prefix",
        default="",
        help="Only analyze subfolders whose names start with this prefix",
    )
    args = parser.parse_args()

    root = Path(args.input_dir)
    folders = sorted(
        p for p in root.iterdir()
        if p.is_dir() and (not args.prefix or p.name.startswith(args.prefix))
    )
    ncores = max(1, min(args.ncores, cpu_count()))
    print(f"Found {len(folders)} folders in {root} (ncores={ncores})")

    if ncores > 1:
        with Pool(processes=ncores) as pool:
            results = pool.map(analyze_folder, folders)
    else:
        results = [analyze_folder(folder) for folder in folders]
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
    plddts = [row["avg_plddt"] for row in results if row["avg_plddt"] is not None]
    if plddts:
        print(f"pLDDT range: {min(plddts):.2f} – {max(plddts):.2f}")


if __name__ == "__main__":
    main()
