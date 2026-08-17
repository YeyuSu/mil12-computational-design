# Helix scoring

`mIL12_alphahelix_dynamic_score.csv` ranks all 500 ProteinMPNN sequences.

```bash
python analyze_helix_plddt.py \
  --input-dir ../03_alphafold3/selected \
  --output /tmp/selected_helix_plddt.csv

python filter_no_rk.py \
  --fasta /path/to/seqs_renamed \
  --pass-fasta /tmp/helix_no_RK.fa \
  --fail-fasta /tmp/helix_with_RK.fa
```

The helix window is located from `DPYRVK` (variable length 14–25). The
fixed-length loop campaign uses
[`../../loop_redesign/analyze_loop_plddt.py`](../../loop_redesign/analyze_loop_plddt.py)
on residues 147–160.

The archived CSV was computed from the full five-sample AF3 jobs. Sequence
columns (`sequence`, `target_subsequence`) have been removed. Re-running the
scorer needs local AF3 coordinate files.
