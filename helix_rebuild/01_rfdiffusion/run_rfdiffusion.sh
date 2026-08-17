#!/usr/bin/env bash
# Reconstruct murine IL-12 p35 loop (B147-160) as a 14-25 residue helix with RFdiffusion.
#
# Both production runs used the SAME contig. Folder names 14-25aa (20 designs) and
# 15-20aa (30 designs) are historical labels; they are not different contig ranges.
#
# Install RFdiffusion: https://github.com/RosettaCommons/RFdiffusion
# Run from an RFdiffusion environment that provides scripts/run_inference.py.

set -euo pipefail

INPUT_PDB="${INPUT_PDB:-../../inputs/AF3predicted_mIL12.pdb}"
OUTDIR="${OUTDIR:-./outputs}"
CONTIG='[A1-313/0 B1-146/14-25/B161-193]'

if [[ ! -f "${INPUT_PDB}" ]]; then
  echo "Input PDB not found: ${INPUT_PDB}" >&2
  exit 1
fi

if ! command -v python >/dev/null 2>&1; then
  echo "python not found" >&2
  exit 1
fi

if [[ ! -f scripts/run_inference.py ]]; then
  echo "Run this script from an RFdiffusion install that contains scripts/run_inference.py" >&2
  echo "Example:" >&2
  echo "  cd /path/to/RFdiffusion" >&2
  echo "  INPUT_PDB=${INPUT_PDB} bash /path/to/this/run_rfdiffusion.sh" >&2
  exit 1
fi

mkdir -p "${OUTDIR}"

# Campaign 1: 20 backbones (archived as 20251007_mIL12_diffusion_14-25aa_*.pdb)
python scripts/run_inference.py \
  "inference.output_prefix=${OUTDIR}/20251007_mIL12_diffusion_14-25aa/20251007_mIL12_diffusion_14-25aa" \
  inference.num_designs=20 \
  "inference.input_pdb=${INPUT_PDB}" \
  "contigmap.contigs=${CONTIG}"

# Campaign 2: 30 backbones (archived as 20251007_mIL12_diffusion_15-20aa_*.pdb)
python scripts/run_inference.py \
  "inference.output_prefix=${OUTDIR}/20251007_mIL12_diffusion_15-20aa/20251007_mIL12_diffusion_15-20aa" \
  inference.num_designs=30 \
  "inference.input_pdb=${INPUT_PDB}" \
  "contigmap.contigs=${CONTIG}"

echo "Wrote 50 backbone models under ${OUTDIR}"
echo "Diffusion trajectories (traj/) are not required for ProteinMPNN or AlphaFold 3."
