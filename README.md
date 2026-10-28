# Acuce Fluidigm manuscript code and Source Data

This repository contains the plotting code and Source Data for the Acuce
intra-specific interaction / Fluidigm manuscript.

## Contents

- `code/`: one plotting script per main or supplementary figure, plus the
  Source Data builder and validation scripts.
- `source_data/Source_Data.xlsx`: per-panel Source Data worksheets.
- `source_data/README.md`: sheet layout, figure-to-sheet mapping, and
  reproducibility notes.

## Reproduce figures

Run commands from the repository root:

```bash
Rscript code/Fig1.R
Rscript code/check_source_data.R
```

Figure scripts read `source_data/Source_Data.xlsx` and write rendered files to
`figures/`. Multi-panel assembly uses the VectorCraft CLI; set
`VECTORCRAFT_CLI` to its executable path when it is not installed at the
default location.

## Source Data verification

The workbook uses two header rows before each data table. Figure scripts read
worksheets with `readxl::read_excel(..., skip = 2)`. See
`source_data/README.md` for the complete panel-to-sheet mapping and validation
criteria.
