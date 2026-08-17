# RFdiffusion backbone generation

Rebuild the protease-sensitive p35 loop (input chain B residues 147–160) as a
variable-length helix connector.

## Contig

```text
contigmap.contigs=[A1-313/0 B1-146/14-25/B161-193]
```

| Segment | Meaning |
|---|---|
| `A1-313` | Keep IL-12 p40 (input chain A) |
| `/0` | Chain break |
| `B1-146` | Keep p35 N-terminal to the loop |
| `14-25` | Diffuse a new 14–25 residue segment |
| `B161-193` | Keep p35 C-terminal motif starting at DPYRVK |

## Two runs, one contig

| Historical folder name | `num_designs` | Files in `outputs/` |
|---|---|---|
| `20251007_mIL12_diffusion_14-25aa` | 20 | `*_0.pdb` … `*_19.pdb` |
| `20251007_mIL12_diffusion_15-20aa` | 30 | `*_0.pdb` … `*_29.pdb` |

The `15-20aa` name is leftover from an earlier plan. Both jobs used `14-25`.
The lab notebook records a single job with `num_designs=50` and prefix
`14-25aa`; this repository archives the 20 + 30 files that were actually kept.

## Chain IDs in the output PDB

RFdiffusion remaps chains relative to the AlphaFold 3 input:

- output chain **A** = p35 (designed helix is on this chain)
- output chain **B** = p40

See [docs/chain_mapping.md](../../docs/chain_mapping.md).

## Reproduce

From an RFdiffusion installation:

```bash
export INPUT_PDB=/path/to/mil12-computational-design/inputs/AF3predicted_mIL12.pdb
bash /path/to/mil12-computational-design/helix_rebuild/01_rfdiffusion/run_rfdiffusion.sh
```

This repository already contains the 50 final backbones (`*.pdb`) and metadata
(`*.trb`). Intermediate diffusion frames (`traj/`) are omitted because they are
not needed to rerun ProteinMPNN or AlphaFold 3.

## Software

[RFdiffusion](https://github.com/RosettaCommons/RFdiffusion) — Watson et al., *Nature* 2023.
