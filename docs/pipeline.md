# Computational pipeline

## Design goal

Murine IL-12 p35 residues 147–160 (`NGETLRQKPPVGEA`) form a solvent-exposed
loop. Mass spectrometry mapped cleavage near `NGETLRQK`. The C-terminal stretch
after this loop contacts both the IL-12 receptor and p40, so the redesign keeps
`DPYRVK` (starting at residue 161) and everything C-terminal of it. Two
computational strategies were used to remove the protease site while keeping the
p35/p40 heterodimer:

1. **Helix rebuild (main campaign).** RFdiffusion replaces B147–160 with a new
   14–25 residue backbone, ProteinMPNN designs the sequence, and AlphaFold 3
   scores helix confidence.
2. **Loop redesign (construct 6).** ProteinMPNN redesigns the existing
   14-residue loop on the wild-type backbone. No backbone diffusion.

```text
AF3predicted_mIL12.pdb
        |
        +--> RFdiffusion contig [A1-313/0 B1-146/14-25/B161-193]
        |         50 backbones
        |         |
        |         v
        |    ProteinMPNN v_48_020, T=0.2, 10 seq/backbone
        |         500 p35 sequences
        |         |
        |         v
        |    AlphaFold 3 (sequence-only JSON, seed=1)
        |         helix pLDDT + complex ipTM
        |         |
        |         v
        |    5 experimental helix candidates
        |
        +--> ProteinMPNN only on WT p35 B147-160 (T=0.2 and 0.3)
                  1000 sequences -> filter no R/K at 147-160
                  AlphaFold 3 -> experimental construct 6 (0_3_sample358)
```

## Helix rebuild details

### RFdiffusion

Input: [`inputs/AF3predicted_mIL12.pdb`](../inputs/AF3predicted_mIL12.pdb)
(chain A = p40, 313 aa; chain B = p35, 193 aa).

Contig used for **both** production jobs:

```text
contigmap.contigs=[A1-313/0 B1-146/14-25/B161-193]
```

| Historical output prefix | `inference.num_designs` |
|---|---|
| `20251007_mIL12_diffusion_14-25aa` | 20 |
| `20251007_mIL12_diffusion_15-20aa` | 30 |

The `15-20aa` folder name does **not** mean a 15–20 residue contig. Both jobs
sampled helix length 14–25. The lab notebook records one SLURM job with
`num_designs=50` and prefix `14-25aa`; the archived files are the 20 + 30 split
above. Commands: [`helix_rebuild/01_rfdiffusion/run_rfdiffusion.sh`](../helix_rebuild/01_rfdiffusion/run_rfdiffusion.sh).

RFdiffusion remaps chains in the output PDB (A = p35, B = p40). See
[chain_mapping.md](chain_mapping.md).

### ProteinMPNN

- Model `v_48_020`, **soluble** weights
- Temperature `0.2`, seed 37, 10 sequences per backbone, omit Cys, batch size 30
- Design positions = CA B-factor `0.00` on the RFdiffusion PDB (inpainted helix;
  14–25 aa). Preprocessor:
  [`helix_rebuild/02_proteinmpnn/process_diffusion_pdbs.py`](../helix_rebuild/02_proteinmpnn/process_diffusion_pdbs.py)
- Designed chain = RFdiffusion chain A (p35)

ProteinMPNN FASTA is **not** stored here. Reproduce locally with
[`helix_rebuild/02_proteinmpnn/run_proteinmpnn.sh`](../helix_rebuild/02_proteinmpnn/run_proteinmpnn.sh).

### AlphaFold 3

JSON dialect `alphafold3` version 3, chain A = p40, chain B = designed p35, no
MSA. The first AF3 JSON attempt reused the wild-type MSA template; that failed
because helix lengths (14–25) no longer match the WT MSA. Production JSON
therefore strips `unpairedMsa` / paired MSA fields. AF3 JSON for designed
sequences is **not** stored here; rebuild it locally with
[`helix_rebuild/03_alphafold3/generate_af3_json.py`](../helix_rebuild/03_alphafold3/generate_af3_json.py).

Scoring (`analyze_helix_plddt.py`):

- Locate `DPYRVK` after p35 residue 147
- Average atom pLDDT over residues `[147, residue before DPYRVK]`
- Average ipTM across the five AF3 samples
- Length of that window must fall in 14–25

All 500 scores (IDs and metrics only; no amino-acid strings):
[`helix_rebuild/04_analysis/scores/mIL12_alphahelix_dynamic_score.csv`](../helix_rebuild/04_analysis/scores/mIL12_alphahelix_dynamic_score.csv).

Two experimental designs (`20_2`, `5_8`) were chosen because the rebuilt helix
contains neither R nor K, not because they had the highest pLDDT.

## Loop redesign details (experimental construct 6)

ProteinMPNN redesigned only p35 B147–160 on the wild-type structure (500
sequences each; `v_48_020`; omit Cys; designed chain B, fixed chain A). T = 0.2
used **soluble** weights; T = 0.3 used **vanilla** weights (construct 6 is T =
0.3 sample 358). Commands:
[`loop_redesign/run_proteinmpnn.sh`](../loop_redesign/run_proteinmpnn.sh).
Fixed-window scorer: [`loop_redesign/analyze_loop_plddt.py`](../loop_redesign/analyze_loop_plddt.py).

Sequences were merged, then split by whether 147–160 contains R or K (371
without, 629 with). **Construct 6** is `mIL12_redesignloop_0_3_sample358`
(no R/K, loop pLDDT 77.16, ipTM 0.848). Designed sequences are not stored here.

Original AF3 jobs for this campaign used a wild-type JSON template that included
MSAs. This repository keeps the score table without sequence columns:
[`loop_redesign/scores/mIL12_redesignloop_147-160_score.csv`](../loop_redesign/scores/mIL12_redesignloop_147-160_score.csv).
The AF3 coordinate file for construct 6 was not kept in the working folder.

## What is not in this repository

- Designed amino-acid sequences (ProteinMPNN FASTA, AF3 JSON, designed CIF)
- RFdiffusion `traj/` frames (~1.3 GB)
- Atom-level AF3 `*_confidences.json`
- The other 495 helix AF3 models (IDs and scores are provided)
- Cloning maps, Prism files, and crosslinking proteomics
