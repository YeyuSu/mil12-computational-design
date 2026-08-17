#!/usr/bin/env bash
# Redesign the wild-type p35 loop (AF3 chain B, residues 147-160) with ProteinMPNN.
# No RFdiffusion: the backbone stays wild-type.
#
# Production settings from the lab notebook:
#   T=0.2 -> soluble_model_weights
#   T=0.3 -> vanilla_model_weights  (construct 6 is T=0.3 sample 358)
#   v_48_020, seed 37, 500 sequences, omit_AAs=C, batch_size=10
#   designed_chains=['B'], fixed_chains=['A']
#
# Install ProteinMPNN: https://github.com/dauparas/ProteinMPNN

set -euo pipefail

PDB="${PDB:-../inputs/AF3predicted_mIL12.pdb}"
OUT_DIR="${OUT_DIR:-.}"
PROTEINMPNN_DIR="${PROTEINMPNN_DIR:-}"

if [[ -z "${PROTEINMPNN_DIR}" ]]; then
  echo "Set PROTEINMPNN_DIR to your ProteinMPNN clone." >&2
  exit 1
fi

mkdir -p "${OUT_DIR}/parsed" "${OUT_DIR}/seqs" "${OUT_DIR}/input_pdb"
cp "${PDB}" "${OUT_DIR}/input_pdb/AF3predicted_mIL12.pdb"

python "${PROTEINMPNN_DIR}/helper_scripts/parse_multiple_chains.py" \
  --input_path "${OUT_DIR}/input_pdb" \
  --output_path "${OUT_DIR}/parsed/parsed_pdbs.jsonl"

python "${PROTEINMPNN_DIR}/helper_scripts/assign_fixed_chains.py" \
  --input_path "${OUT_DIR}/parsed/parsed_pdbs.jsonl" \
  --output_path "${OUT_DIR}/parsed/assigned_pdbs.jsonl" \
  --chain_list "B"

FIXED_B=$(python - <<'PY'
print(" ".join(str(i) for i in list(range(1, 147)) + list(range(161, 194))))
PY
)
python "${PROTEINMPNN_DIR}/helper_scripts/make_fixed_positions_dict.py" \
  --input_path "${OUT_DIR}/parsed/parsed_pdbs.jsonl" \
  --output_path "${OUT_DIR}/parsed/fixed_positions.jsonl" \
  --chain_list "B" \
  --position_list "${FIXED_B}"

run_one() {
  local T="$1"
  local WEIGHTS="$2"
  local WEIGHTS_ARGS=()
  if [[ -d "${WEIGHTS}" ]]; then
    WEIGHTS_ARGS+=(--path_to_model_weights "${WEIGHTS}")
  else
    echo "Warning: weights directory not found (${WEIGHTS}); using ProteinMPNN default." >&2
  fi
  python "${PROTEINMPNN_DIR}/protein_mpnn_run.py" \
    --jsonl_path "${OUT_DIR}/parsed/parsed_pdbs.jsonl" \
    --chain_id_jsonl "${OUT_DIR}/parsed/assigned_pdbs.jsonl" \
    --fixed_positions_jsonl "${OUT_DIR}/parsed/fixed_positions.jsonl" \
    --out_folder "${OUT_DIR}" \
    --num_seq_per_target 500 \
    --sampling_temp "${T}" \
    --seed 37 \
    --batch_size 10 \
    --model_name v_48_020 \
    --omit_AAs C \
    --save_score 1 \
    --save_probs 1 \
    "${WEIGHTS_ARGS[@]}"
}

run_one 0.2 "${PROTEINMPNN_DIR}/soluble_model_weights"
run_one 0.3 "${PROTEINMPNN_DIR}/vanilla_model_weights"

echo "Designed sequences are in ${OUT_DIR}/seqs"
echo "This public repository does not archive those FASTA files."
