# Helix rebuild campaign

RFdiffusion → ProteinMPNN → AlphaFold 3 reconstruction of the p35 147–160 loop
as a 14–25 residue alpha helix.

| Step | Folder |
|---|---|
| 1. Backbone inpainting | [01_rfdiffusion](01_rfdiffusion) |
| 2. Sequence design | [02_proteinmpnn](02_proteinmpnn) |
| 3. Structure validation | [03_alphafold3](03_alphafold3) |
| 4. Helix pLDDT / no-R/K filters | [04_analysis](04_analysis) |

See [docs/pipeline.md](../docs/pipeline.md) for parameters and
[designs/selected_designs.csv](../designs/selected_designs.csv) for the five
experimental helix IDs and scores. Designed sequences are not stored here.
