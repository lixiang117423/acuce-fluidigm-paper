# Acuce Fluidigm manuscript data and code

Public data and code repository for:

**Small and context-dependent effects of intra-specific neighbours on rice disease, growth and leaf transcription**

Current release: **v47** (29 October 2026). This release keeps exactly the same scientific content as v46 but uses English-only file and folder names for international collaboration and journal submission.

## Manuscript

- `manuscript/fluidigm-paper-v47_20261029.docx`
- `manuscript/fluidigm-paper-v47_20261029.pdf`
- `manuscript/response_second_round_v47_20261029.docx`
- `manuscript/response_second_round_v47_20261029.pdf`

## Contents

- `data/source_data/`: Source Data workbook (47 sheets) and its README.
- `data/processed_inputs/`: four processed analysis input tables.
- `data/raw_cq/`: three raw Fluidigm Cq exports (Supplementary Data 1-3).
- `code/analysis/`: data processing, inoculation and field disease, agronomic traits, expression and robustness analyses.
- `code/figure_generation/`: scripts for final Figures 1-5 and Supplementary Figures S1-S16.
- `code/reporting/`: archival manuscript, Source Data and response-generation scripts.
- `figures/main/`: final main figures in PDF and PNG.
- `figures/supplementary/`: final Supplementary Figures S1-S16 in PDF and 600 dpi PNG.
- `tables/`: Supplementary Tables S1-S9.
- `docs/review_audit_v47.md`: item-by-item audit against the second-round review.
- `environment/sessionInfo_v47.txt`: R session information.

## Reproducibility

- Run analysis and figure-generation commands from the repository root.
- R scripts use the repository root as the working directory.
- Processed input tables are provided in `data/processed_inputs/`.
- Figure-generation scripts read `data/source_data/Source_Data_v45.xlsx`.
- Script comments and historical change notes may retain Chinese text; all file and folder names are ASCII-only.

## Integrity

`MANIFEST_MD5_SHA256.tsv` records the byte count, MD5 and SHA-256 for every tracked file except the manifest itself.

## Contact

Repository owner: Xiang Li (`lixiang117423`).
