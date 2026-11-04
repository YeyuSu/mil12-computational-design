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

## System requirements and installation

The CPU preprocessing demo below uses only the Python standard library and was
tested with Python 3.9.6 on macOS 15.7.9. Other analysis scripts use modern
type annotations and need Python 3.10 or newer; use Python 3.12 for the local
analysis environment. `environment.yml` also installs pandas for local table
work. The complete design campaign has **not** been tested in this repository on a
standard desktop computer, and exact versions of the three historical model
installations were not archived.

Install the local analysis environment with Conda:

```bash
conda env create -f environment.yml
conda activate mil12-computational-design
```

`environment.yml` does **not** install RFdiffusion, ProteinMPNN, AlphaFold 3,
their model weights, or their CUDA dependencies. Install each separately using
its official instructions, in a compatible GPU environment:

- [RFdiffusion](https://github.com/RosettaCommons/RFdiffusion)
- [ProteinMPNN](https://github.com/dauparas/ProteinMPNN)
- [AlphaFold 3](https://github.com/google-deepmind/alphafold3)

The full backbone generation, sequence design, and local AlphaFold 3 prediction
workflow requires an NVIDIA GPU and software compatible with the respective
upstream model releases. GPU memory and runtime depend on the chosen model and
job size; this repository does not establish a tested minimum. AlphaFold Server
is an alternative submission route for prediction, subject to its own access
and terms. **Typical installation time and full-design runtime on a standard
desktop computer: Not applicable to a standard desktop computer; GPU required.**
No measured GPU installation or campaign time is available from the archived
run, so no numerical estimate is given.

## CPU demo with included data

From the repository root, the following step checks all 50 included
RFdiffusion PDB backbones and identifies the positions to be redesigned. It
does not require a GPU, model installation, or designed sequences:

```bash
python helix_rebuild/02_proteinmpnn/process_diffusion_pdbs.py \
  --pdb-dir helix_rebuild/01_rfdiffusion/outputs \
  --outdir /tmp/mil12-cpu-demo
```

Expected final line: `Wrote dictionaries for 50/50 PDBs to
/tmp/mil12-cpu-demo`. The output directory contains `assigned_chains.jsonl`
and `fixed_positions.jsonl`. On macOS 15.7.9 with Python 3.9.6, this command
took 0.08 seconds of wall time in one local run (`/usr/bin/time -p`); this is
an observed value, not a general desktop benchmark. The Conda installation
time was not measured.

## Full design workflow and local scripts

Run the model stages in the order described in [docs/pipeline.md](docs/pipeline.md).
Example production commands are in
`helix_rebuild/01_rfdiffusion/run_rfdiffusion.sh`,
`helix_rebuild/02_proteinmpnn/run_proteinmpnn.sh`, and
`loop_redesign/run_proteinmpnn.sh`. After ProteinMPNN creates a FASTA locally,
build the helix AlphaFold 3 JSON inputs with:

```bash
# After running ProteinMPNN locally, build AF3 JSON from FASTA
python helix_rebuild/03_alphafold3/generate_af3_json.py \
  --fasta /path/to/seqs_renamed \
  --template inputs/af3_sequence_template.json \
  --outdir /tmp/helix_json
```

The FASTA and designed AF3 inputs are intentionally absent from this public
repository. Thus, the JSON-generation command cannot run directly from the
checkout, and the historical AF3 scores cannot be fully regenerated from the
archived files alone. To run the workflow on another protein, supply your own
starting structure, adapt the hard-coded p35/p40 chain IDs and residue windows
in the scripts, and regenerate sequences and AF3 inputs; the published
parameters are specific to murine IL-12.

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
