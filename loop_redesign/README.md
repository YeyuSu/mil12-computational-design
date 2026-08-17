# Loop redesign (fixed-length p35 147–160)

This chapter is a **sequence-only** alternative to helix rebuilding. ProteinMPNN
redesigns the existing 14-residue p35 loop on the wild-type AlphaFold 3 backbone
instead of generating a new helix with RFdiffusion.

Wild-type p35 147–160: `NGETLRQKPPVGEA` (contains K). Residue 161 starts the
conserved `DPYRVK` motif.

**Experimental construct 6** is `mIL12_redesignloop_0_3_sample358`.
It is one of the six mutants that were cloned. Designed sequences are not
stored in this repository.

## Pipeline

1. ProteinMPNN on `../inputs/AF3predicted_mIL12.pdb`, redesign p35 147–160
   (AF3 chain B; chain A / p40 fixed). Example:

   ```bash
   export PROTEINMPNN_DIR=/path/to/ProteinMPNN
   bash run_proteinmpnn.sh
   ```

   Designed FASTA is not stored in this repository (500 designed + WT at each
   temperature).
2. Merge and rename:

   ```bash
   python process_fasta.py
   ```

   Output: a combined FASTA of 1000 sequences (local only).
3. Keep designs whose 147–160 segment has no R or K:

   ```bash
   python filter_sequences_no_rk.py
   ```

   Production split: 371 no-R/K and 629 with R/K (FASTA not stored here).
4. Build AF3 JSON locally:

   ```bash
   python generate_af3_json.py \
     --fasta seqs/mIL12_redesignloop_no_RK_147-160.fa \
     --template ../inputs/af3_sequence_template.json \
     --outdir json_inputs
   ```

   Original loop predictions used an MSA-containing wild-type AF3 job as the
   JSON template (`mil12_wt_data.json`). Those bulky MSA fields are not stored
   here. Scores in `scores/mIL12_redesignloop_147-160_score.csv` come from the
   original AF3 run (sequence columns removed).

## Selected structures

| Folder | Role |
|---|---|
| `selected/0_3_sample358/` | **Construct 6** (ProteinMPNN only; scores in the CSV) |
| `selected/wt/` | Wild-type mIL-12 AF3 model |
| `selected/0_2_sample1/` | Example T=0.2 AF3 job scores (not an experimental clone) |

Construct 6 files: [`selected/0_3_sample358/`](selected/0_3_sample358/).

## ProteinMPNN settings

From FASTA headers and the lab notebook (*Mouse IL12 loop rebuild*):

| Parameter | T = 0.2 | T = 0.3 (construct 6) |
|---|---|---|
| Model | `v_48_020` | `v_48_020` |
| Weights | `soluble_model_weights` | `vanilla_model_weights` |
| Designed chain | AF3 chain **B** (p35) | same |
| Fixed chain | AF3 chain **A** (p40) | same |
| Designed window | p35 147–160 | same |
| Sequences | 500 | 500 |
| Seed | 37 | 37 |
| Omit | Cys (`--omit_AAs C`) | Cys |
| Batch size | 10 | 10 |

Loop pLDDT over the fixed 147–160 window:

```bash
python analyze_loop_plddt.py \
  --input-dir selected \
  --output /tmp/loop_plddt.csv
```

