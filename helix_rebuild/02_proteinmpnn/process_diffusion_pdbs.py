#!/usr/bin/env python3
"""Identify ProteinMPNN design positions from RFdiffusion B-factors.

RFdiffusion writes B-factor 0.00 on the inpainted helix and 1.00 on fixed
residues. This script checks that each backbone has 14-25 designable CA atoms
and writes ProteinMPNN-style JSONL dictionaries.

Matches the production preprocessing described in the lab notebook
(Mouse IL12 loop rebuild).
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path


def parse_ca_bfactors(pdb_file: Path) -> dict[str, dict[int, float]]:
    chains: dict[str, dict[int, float]] = defaultdict(dict)
    with pdb_file.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            if line[12:16].strip() != "CA":
                continue
            chain = line[21:22].strip() or "A"
            resid = int(line[22:26])
            bfactor = float(line[60:66])
            chains[chain][resid] = bfactor
    return chains


def design_info(chains: dict[str, dict[int, float]], min_length: int, max_length: int):
    info = {}
    for chain_id, residues in chains.items():
        sorted_resids = sorted(residues)
        designable = []
        fixed = []
        for idx, resid in enumerate(sorted_resids, start=1):
            if abs(residues[resid]) < 0.01:
                designable.append(idx)
            else:
                fixed.append(idx)
        n = len(designable)
        info[chain_id] = {
            "valid": min_length <= n <= max_length,
            "num_designable": n,
            "designable_positions": designable,
            "fixed_positions": fixed,
            "total_residues": len(sorted_resids),
        }
    return info


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb-dir", required=True, help="Directory of RFdiffusion PDBs")
    parser.add_argument("--outdir", default="parsed", help="Output directory for JSONL files")
    parser.add_argument("--min-length", type=int, default=14)
    parser.add_argument("--max-length", type=int, default=25)
    args = parser.parse_args()

    pdb_dir = Path(args.pdb_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    pdbs = sorted(pdb_dir.glob("*.pdb"))
    assigned = {}
    fixed = {}
    n_ok = 0
    for pdb in pdbs:
        info = design_info(parse_ca_bfactors(pdb), args.min_length, args.max_length)
        valid = [chain for chain, row in info.items() if row["valid"]]
        if not valid:
            print(f"SKIP {pdb.name}: no chain with {args.min_length}-{args.max_length} designable residues")
            continue
        n_ok += 1
        stem = pdb.stem
        assigned[stem] = [valid, [c for c in info if c not in valid]]
        fixed[stem] = {chain: row["fixed_positions"] for chain, row in info.items()}
        for chain, row in info.items():
            flag = "design" if row["valid"] else "fixed"
            print(f"{pdb.name} chain {chain}: {flag} {row['num_designable']}/{row['total_residues']}")

    (outdir / "assigned_chains.jsonl").write_text(json.dumps(assigned) + "\n", encoding="utf-8")
    (outdir / "fixed_positions.jsonl").write_text(json.dumps(fixed) + "\n", encoding="utf-8")
    print(f"Wrote dictionaries for {n_ok}/{len(pdbs)} PDBs to {outdir}")


if __name__ == "__main__":
    main()
