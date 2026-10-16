# Figures

Scripts that regenerate the figure panels of the paper from the Zenodo objects. Each script writes the panels as PDFs.

## Scripts

| Script | Panels |
|---|---|
| `fig1CD_figS2EF_pilot.R` | Fig 1C–D, Supp Fig 2E–F |
| `fig2_cell_types.R` | Fig 2A–C, Supp Fig 5, quasi-binomial tests of Fig 2B |
| `figS3_qc_histograms.R` | Supp Fig 3 (needs the GEO `cell_feature_matrix.h5` files) |
| `fig3_fig4B_plasma_cells.R` | Fig 3A–E, Fig 4B, Supp Figs 6, 8, 9, Supp Table 3 |
| `fig5_stromal_cells.R` | Fig 5A–E, Supp Fig 10 |
| `fig6_scider_neighbourhoods.R` | Fig 6A–B, Supp Figs 11–13 |

- `setup.R` loads the packages and defines the annotation directory, colour palettes and the paper's cluster labels;
- `config.R` defines paths, the sample sheet and helper functions.

Analysis scripts are in [`../Xenium_analysis`](../Xenium_analysis)

## Running the scripts

```r
setwd("path/to/this/repository/Figures")
source("fig1CD_figS2EF_pilot.R")
source("fig2_cell_types.R")
source("figS3_qc_histograms.R")
source("fig3_fig4B_plasma_cells.R")
source("fig5_stromal_cells.R")
source("fig6_scider_neighbourhoods.R")
```
