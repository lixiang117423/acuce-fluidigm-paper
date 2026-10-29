# Acuce Fluidigm data and code

Repository for the analyses of:

**Rice intraspecific-neighbour interactions: disease, growth and leaf transcription**

**Status:** private development repository. The manuscript and response letter are withheld until publication. The repository may be made public after the paper is accepted or published.

## Contents

- `data/source_data/`: Source Data workbook (47 sheets) and its README.
- `data/processed_inputs/`: four processed analysis input tables.
- `data/raw_cq/`: three raw Fluidigm Cq exports (Supplementary Data 1-3).
- `code/analysis/`: data processing, inoculation and field disease, agronomic traits, expression and robustness analyses.
- `code/figure_generation/`: scripts for final Figures 1-5 and Supplementary Figures S1-S16.
- `figures/main/`: final main figures in PDF and PNG.
- `figures/supplementary/`: final Supplementary Figures S1-S16 in PDF and 600 dpi PNG.
- `tables/`: Supplementary Tables S1-S9.
- `environment/sessionInfo_v47.txt`: R session information.

## Reproducibility

- Run analysis and figure-generation commands from the repository root.
- R scripts use the repository root as the working directory.
- Processed input tables are provided in `data/processed_inputs/`.
- Figure-generation scripts read `data/source_data/Source_Data_v45.xlsx`.
- Script comments may retain Chinese text; all file and folder names are ASCII-only.

## Integrity

`MANIFEST_MD5_SHA256.tsv` records the byte count, MD5 and SHA-256 for the analysis files except the manifest itself.

## Contact

Repository owner: Xiang Li (`lixiang117423`).
