# Acuce Fluidigm manuscript data and code

Public data and code repository for:

**Small and context-dependent effects of intra-specific neighbours on rice disease, growth and leaf transcription**

Current release: **v46** (29 October 2026)

## Manuscript

- `manuscript/fluidigm-paper-v46_20261029.docx`
- `manuscript/fluidigm-paper-v46_20261029.pdf`
- `manuscript/response_second_round_v46_20261029.docx`
- `manuscript/response_second_round_v46_20261029.pdf`

## Contents

- `data/source_data/`: Source_Data_v45.xlsx (47 sheets) and its README.
- `data/processed_inputs/`: four processed analysis input tables.
- `data/raw_cq/`: three raw Fluidigm Cq exports (Supplementary Data 1-3).
- `code/analysis/`: data-processing, inoculation/field disease, agronomic trait, expression, and robustness scripts.
- `code/figure_generation/`: scripts for the final Figures 1-5 and Supplementary Figures S1-S16.
- `code/reporting/`: manuscript, Source Data and response generation scripts from v31 onwards.
- `figures/main/`: final main figures in PDF and PNG.
- `figures/supplementary/`: final Supplementary Figures S1-S16 in PDF and 600 dpi PNG.
- `tables/`: Supplementary Tables S1-S9.
- `docs/review_audit_v46.md`: item-by-item audit against the second-round review.
- `environment/sessionInfo_v46.txt`: R session information.

## Reproducibility notes

- Figure S2A-S2F are based on pot-level observations and pot-level Sidak letters.
- The three-stage Fluidigm panel covers tillering, flowering and grain filling.
- The root-exudate hormone figures from older drafts are not part of the current paper.
- A matched inoculation-plus-expression experiment was not performed; this is stated as a limitation in the manuscript.
- Analysis scripts retain the absolute project-root variables from the original working environment; update only those path variables when running elsewhere.

## Integrity

`MANIFEST_MD5_SHA256.tsv` records byte count, MD5 and SHA-256 for every file except the manifest itself.

## Contact

Repository owner: Xiang Li (`lixiang117423`).
