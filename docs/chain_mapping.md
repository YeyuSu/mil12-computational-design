# Chain ID mapping

The AlphaFold 3 input/output convention is **not** the same as the RFdiffusion
output convention. ProteinMPNN headers refer to the RFdiffusion PDB.

The redesigned loop/helix is on **p35 (Il12a)**, the shorter subunit, not p40.

## AlphaFold 3 JSON and input PDB

File: [`inputs/AF3predicted_mIL12.pdb`](../inputs/AF3predicted_mIL12.pdb) and
[`inputs/af3_sequence_template.json`](../inputs/af3_sequence_template.json).

| AF3 chain | Subunit | Length (WT) | N-terminus |
|---|---|---|---|
| **A** | IL-12 p40 (Il12b) | 313 | `MWELEKDVYV...` |
| **B** | IL-12 p35 (Il12a) | 193 (WT); 193–204 after helix insertion | `RVIPVSGPAR...` |

The designed segment lives on **AF3 chain B (p35)**, residues 147 until the
residue before `DPYRVK`.

Wild-type p35 147–166:

```text
147        160 161
NGETLRQKPPVGEA DPYRVK
```

## RFdiffusion output PDB

RFdiffusion concatenates the contig and remaps chain IDs:

| RFdiffusion chain | Subunit | Typical residue count |
|---|---|---|
| **A** | p35, including the new helix | ~194–205 |
| **B** | p40 | 313 |

ProteinMPNN FASTA headers therefore say `designed_chains=['A']` and
`fixed_chains=['B']`. That **A** is p35, not AF3 chain A.

When sequences are converted to AF3 JSON, the designed p35 string is written
back onto **AF3 chain B**.

## Quick check

A designed helix string must appear in:

- ProteinMPNN FASTA (p35 sequence, RFdiffusion chain A)
- AF3 JSON `sequences[1].protein.sequence` (chain B)
- AF3 CIF chain B around residues 147–167

Those designed-sequence files are not redistributed in this repository.
