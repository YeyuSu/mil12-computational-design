# Computational redesign of murine IL-12

Supporting code and design files for rebuilding the protease-sensitive p35 loop
of murine IL-12 as an alpha helix with **RFdiffusion**, **ProteinMPNN**, and
**AlphaFold 3**. A fixed-length loop redesign campaign is included as a second
chapter.

This repository is meant to accompany the manuscript. It is not a wrapper around
the three third-party models; install those from their official sources.

## Design idea

Wild-type mIL-12 p35 residues 147–160 (`NGETLRQKPPVGEA`) sit on an exposed loop
adjacent to the conserved `DPYRVK` motif. Mass spectrometry mapped cleavage near
`NGETLRQK`. Sequence C-terminal of the loop contacts the receptor and p40, so
both campaigns keep `DPYRVK` and the rest of p35. The main campaign **deletes
that loop and inpaints a 14–25 residue helix**, then designs and validates
sequences:

1. **RFdiffusion** — contig `[A1-313/0 B1-146/14-25/B161-193]` → 50 backbones
2. **ProteinMPNN** — soluble `v_48_020`, T = 0.2, omit Cys, 10 sequences / backbone → 500 designs
3. **AlphaFold 3** — sequence-only JSON (helix length ≠ WT MSA); rank by helix pLDDT and ipTM

Five helix designs (RFdiffusion + ProteinMPNN + AlphaFold 3) and **one loop
design (ProteinMPNN only, wild-type backbone)** were advanced experimentally.

Method details: [docs/pipeline.md](docs/pipeline.md).  
Chain IDs differ between RFdiffusion and AlphaFold 3: [docs/chain_mapping.md](docs/chain_mapping.md).

## Repository layout

```text
inputs/                         WT AF3 PDB + WT AF3 JSON template
helix_rebuild/
  01_rfdiffusion/               50 final backbones + run script
  02_proteinmpnn/               ProteinMPNN run script (sequences omitted)
  03_alphafold3/                AF3 JSON builder + selected-job scores
  04_analysis/                  helix pLDDT table (no sequences) and filters
loop_redesign/                  ProteinMPNN-only loop campaign + construct 6
designs/                        experimental candidate IDs and scores
docs/                           pipeline, chain map, AF3 output terms
```

## Experimental candidates (6)

Designed amino-acid sequences are **not** in this repository.

| Label | Campaign | ID | Length | pLDDT | ipTM | No R/K |
|---|---|---|---|---|---|---|
| 1 | helix (RFD+MPNN+AF3) | `26_7` | 21 | 90.1 | 0.85 | no |
| 5 | helix (RFD+MPNN+AF3) | `47_2` | 20 | 89.3 | 0.86 | no |
| 11 | helix (RFD+MPNN+AF3) | `38_10` | 19 | 88.4 | 0.85 | no |
| noRK1 | helix (RFD+MPNN+AF3) | `20_2` | 14 | 78.7 | 0.86 | yes |
| noRK2 | helix (RFD+MPNN+AF3) | `5_8` | 15 | 74.8 | 0.85 | yes |
| 6 | loop (ProteinMPNN only) | `0_3_sample358` | 14 | 77.2 | 0.85 | yes |

Construct 6 keeps the wild-type p35 backbone and only redesigns residues 147–160.

Score table (no sequences): [designs/selected_designs.csv](designs/selected_designs.csv).

## Reproduce the analysis scripts

```bash
conda env create -f environment.yml
conda activate mil12-computational-design

# After running ProteinMPNN locally, build AF3 JSON from FASTA
python helix_rebuild/03_alphafold3/generate_af3_json.py \
  --fasta /path/to/seqs_renamed \
  --template inputs/af3_sequence_template.json \
  --outdir /tmp/helix_json
```

RFdiffusion / ProteinMPNN / AlphaFold 3 themselves need GPU installations:

- [RFdiffusion](https://github.com/RosettaCommons/RFdiffusion)
- [ProteinMPNN](https://github.com/dauparas/ProteinMPNN)
- [AlphaFold 3](https://github.com/google-deepmind/alphafold3)

Example commands: `helix_rebuild/01_rfdiffusion/run_rfdiffusion.sh` and
`helix_rebuild/02_proteinmpnn/run_proteinmpnn.sh`.

## Data notes

- **Designed sequences are omitted** (ProteinMPNN FASTA, AF3 JSON inputs, CIF
  models of designed complexes, and sequence columns in score tables). Wild-type
  p35/p40 sequence is kept in `inputs/` because it is the public starting point.
- Final RFdiffusion PDBs are included (inpainted helix is poly-G); diffusion
  trajectories (`traj/`) are not.
- Score tables keep design IDs, lengths, pLDDT, and ipTM.
- AlphaFold 3 structures, if redistributed later, are subject to
  [docs/alphafold3_output_terms.md](docs/alphafold3_output_terms.md).

## Citation

Cite the manuscript (DOI to be added in `CITATION.cff`) and:

- Watson et al., *Nature* (2023) — RFdiffusion
- Dauparas et al., *Science* (2022) — ProteinMPNN
- Abramson et al., *Nature* (2024) — AlphaFold 3
