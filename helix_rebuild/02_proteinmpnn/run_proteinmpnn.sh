#!/usr/bin/env bash
# Sequence-design the RFdiffusion p35 helix with ProteinMPNN.
#
# Production settings from the lab notebook (Mouse IL12 loop rebuild):
#   model v_48_020, soluble_model_weights, T=0.2, seed=37
#   10 sequences / backbone, omit_AAs=C, batch_size=30
#   design positions = CA B-factor 0.00 (inpainted helix), 14-25 aa
#
# Install ProteinMPNN: https://github.com/dauparas/ProteinMPNN

set -euo pipefail

PDB_DIR="${PDB_DIR:-../01_rfdiffusion/outputs}"
OUT_DIR="${OUT_DIR:-.}"
PROTEINMPNN_DIR="${PROTEINMPNN_DIR:-}"
WEIGHTS_DIR="${WEIGHTS_DIR:-${PROTEINMPNN_DIR:+$PROTEINMPNN_DIR/soluble_model_weights}}"

if [[ -z "${PROTEINMPNN_DIR}" ]]; then
  echo "Set PROTEINMPNN_DIR to your ProteinMPNN clone." >&2
  exit 1
fi

mkdir -p "${OUT_DIR}/parsed" "${OUT_DIR}/seqs"

python "$(dirname "$0")/process_diffusion_pdbs.py" \
  --pdb-dir "${PDB_DIR}" \
  --outdir "${OUT_DIR}/parsed"

python "${PROTEINMPNN_DIR}/helper_scripts/parse_multiple_chains.py" \
  --input_path "${PDB_DIR}" \
  --output_path "${OUT_DIR}/parsed/parsed_pdbs.jsonl"

WEIGHTS_ARGS=()
if [[ -d "${WEIGHTS_DIR}" ]]; then
  WEIGHTS_ARGS+=(--path_to_model_weights "${WEIGHTS_DIR}")
fi

python "${PROTEINMPNN_DIR}/protein_mpnn_run.py" \
  --jsonl_path "${OUT_DIR}/parsed/parsed_pdbs.jsonl" \
  --chain_id_jsonl "${OUT_DIR}/parsed/assigned_chains.jsonl" \
  --fixed_positions_jsonl "${OUT_DIR}/parsed/fixed_positions.jsonl" \
  --out_folder "${OUT_DIR}" \
  --num_seq_per_target 10 \
  --sampling_temp "0.2" \
  --seed 37 \
  --batch_size 30 \
  --model_name v_48_020 \
  --omit_AAs C \
  --save_score 1 \
  --save_probs 1 \
  "${WEIGHTS_ARGS[@]}"

echo "Designed sequences are in ${OUT_DIR}/seqs"
echo "This public repository does not archive those FASTA files."
