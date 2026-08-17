# ProteinMPNN sequence design

Design amino acids onto the 50 RFdiffusion p35 helix backbones.

RFdiffusion marks the inpainted helix with CA B-factor `0.00` (fixed residues
are `1.00`). [`process_diffusion_pdbs.py`](process_diffusion_pdbs.py) keeps
backbones whose designable window is 14–25 residues and writes ProteinMPNN
`assigned_chains.jsonl` / `fixed_positions.jsonl`.

## Archived parameters

Taken from the original FASTA headers and the lab notebook (*Mouse IL12 loop rebuild*):

| Parameter | Value |
|---|---|
| Model | `v_48_020` |
| Weights | `soluble_model_weights` |
| Temperature | `0.2` |
| Sequences per backbone | 10 |
| Seed | 37 |
| Omit | Cys (`--omit_AAs C`) |
| Batch size | 30 |
| Designed positions | CA B-factor 0.00 on RFdiffusion chain **A** (p35) |
| Fixed chain | RFdiffusion chain **B** (p40) |

That yields **50 × 10 = 500** designed p35 sequences.

## Files

- `process_diffusion_pdbs.py` — B-factor preprocessor
- `run_proteinmpnn.sh` — production command

Designed FASTA (`seqs/`, `seqs_renamed/`) is **not** stored in this repository.
The first record in each ProteinMPNN FASTA is the poly-G backbone sequence;
designed samples are `T=0.2, sample=1` … `sample=10`.

## Reproduce

```bash
export PROTEINMPNN_DIR=/path/to/ProteinMPNN
bash run_proteinmpnn.sh
```

Install: [dauparas/ProteinMPNN](https://github.com/dauparas/ProteinMPNN) —
Dauparas et al., *Science* 2022.

## Next step

Convert designed p35 sequences to AlphaFold 3 JSON with
[`../03_alphafold3/generate_af3_json.py`](../03_alphafold3/generate_af3_json.py).
