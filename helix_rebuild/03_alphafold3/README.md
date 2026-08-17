# AlphaFold 3 validation

Predict p35/p40 complexes for the 500 helix ProteinMPNN sequences. Designed
sequences and AF3 coordinate files are **not** stored here. The sixth
experimental mutant is the ProteinMPNN-only loop design
`mIL12_redesignloop_0_3_sample358`.

## JSON dialect

Production helix JSON files are sequence-only AlphaFold 3 inputs
(`dialect=alphafold3`, `version=3`):

- chain **A**: wild-type p40 (313 aa)
- chain **B**: designed p35
- `modelSeeds`: `[1]`
- no MSA fields

A first JSON attempt copied the wild-type AF3 MSA template. That failed because
the inpainted helix is 14–25 residues, not the wild-type 14, so the MSA length
no longer matches. Production inputs therefore drop `unpairedMsa` and paired MSA
fields and run as sequence-only jobs.

Designed JSON is not archived. Rebuild it locally from ProteinMPNN FASTA:

```bash
python generate_af3_json.py \
  --fasta /path/to/seqs_renamed \
  --template ../../../inputs/af3_sequence_template.json \
  --outdir ./json_inputs
```

Then submit each JSON to AlphaFold 3 or AlphaFold Server.

## Selected jobs

`selected/` keeps AlphaFold 3 **scores** for the five experimental helix designs
(summary confidences + ranking table). Coordinate files and `*_data.json` are
omitted because they encode designed sequences.

| Folder | Computational ID | No R/K |
|---|---|---|
| `selected/26_7/` | `mil12_alphahelix_14-25aa_26_7` | no |
| `selected/47_2/` | `mil12_alphahelix_14-25aa_47_2` | no |
| `selected/38_10/` | `mil12_alphahelix_14-25aa_38_10` | no |
| `selected/20_2_noRK/` | `mil12_alphahelix_14-25aa_20_2` | yes |
| `selected/5_8_noRK/` | `mil12_alphahelix_14-25aa_5_8` | yes |

Scores for all 500 designs (no amino-acid strings) are in
[`../04_analysis/scores/mIL12_alphahelix_dynamic_score.csv`](../04_analysis/scores/mIL12_alphahelix_dynamic_score.csv).

AlphaFold 3 outputs are subject to
[docs/alphafold3_output_terms.md](../../../docs/alphafold3_output_terms.md).

## Software

[AlphaFold 3](https://github.com/google-deepmind/alphafold3) — Abramson et al., *Nature* 2024.
